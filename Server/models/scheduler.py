# -*- coding: utf-8 -*-
import gzip, cStringIO, os, StringIO, zipfile, math
import boto
import traceback
from boto.s3.connection import S3Connection
from boto.s3.key import Key
import Eval
import sys
import datetime
import pytz

def log_error(userid, error_msg):
    user_row = db(db.user_data.user_id == userid).select().first()

    if db((db.error_reports.user_id == userid)&(db.error_reports.course_name == user_row.course_name)&(db.error_reports.course_term == user_row.course_term)).count() > 0:
        db((db.error_reports.user_id == userid)&(db.error_reports.course_name == user_row.course_name)&(db.error_reports.course_term == user_row.course_term)).update(error=str(error_msg)[:255])
    else:
        db.error_reports.insert(user_id=userid, course_name=user_row.course_name, course_term=user_row.course_term, error=str(error_msg)[:255])
    db.commit()

def clean_database():
    yesterday = datetime.datetime.now() - datetime.timedelta(days=1)
    db(db.auth_event.time_stamp < yesterday).delete()
    db(db.scheduler_task.start_time < yesterday).delete()
    db(db.scheduler_run.start_time < yesterday).delete()

def process_file(file_to_process):
    #only have 10,000 rows allowed for this heroku account, so this function is used to regularly delete unnecessary rows that accumulate in some of our tables

    conn = S3Connection(os.environ['S3_KEY'], os.environ['S3_SECRET'])

    #the genome id should start with e for easy, m for medium, h for hard, w for warmup, and a for assembly
    #diff_dict = {'e':'Easy', 'm':'Medium', 'h':'Hard', 'w':'Warmup', 'a':'Assembly'} #converts difficulty to easier to ready format

    '''
    #first process admin bucket
    bucket = conn.get_bucket(os.environ['ADMIN_UPLOAD_BUCKET'])
    rs = bucket.list() #lists the contents of our directory on S3

    #I'm assuming that only one course will be listed as active. The other courses should represent courses from previous terms
    course_row = db(db.course_data.current_course == True).select().first()
    course_term = course_row.course_quarter + ' ' + course_row.course_year

    #sort the keys by last modified time
    #process each key in order from oldest to newest
    rs = sorted(rs, key=lambda k: k.last_modified)
    for key in rs:
        f = StringIO.StringIO() #gives us a file-like object with which to receive the file from S3
        filename = key.name #name of the file

        # files that start with ans are used when scoring student submissions and should be stored in the genome table
        if filename.startswith('ans_'):
            genome_type = filename.split('_')[1].title() #file should be named ans_{type}_{difficulty}_{genome id #}_chr_{chromosome #} (e.g. ans_random_E_1_chr_1)
            genome_diff = filename.split('_')[2].lower()
            genome_id = filename.split('_')[3]
            chromosome_id = filename.split('_')[5]
            if '.' in chromosome_id:
                chromosome_id = chromosome_id.split('.')[0]

            diff = diff_dict[genome_diff]

            key.get_contents_to_file(f) #loads the file in our file-like object
            f.seek(0)
            #if there is already a record for this genome, then upload the new answer key
            if db((db.genome_table.genome_id == genome_id)&(db.genome_table.chromosome_id == chromosome_id)&(db.genome_table.course_name == course_row.course_name)&(db.genome_table.course_term == course_term)&
                    (db.genome_table.difficulty == diff)&(db.genome_table.genome_type == genome_type)).count() > 0:
                db((db.genome_table.genome_id == genome_id)&(db.genome_table.chromosome_id == chromosome_id)&(db.genome_table.course_name == course_row.course_name)&(db.genome_table.course_term == course_term)&
                   (db.genome_table.difficulty == diff)&
                   (db.genome_table.genome_type == genome_type)).update(answer_key_blob=f.read())
                db.commit()
            #else insert a new record
            else:
                db.genome_table.insert(genome_id=genome_id, chromosome_id=chromosome_id, difficulty=diff, genome_type=genome_type, course_name=course_row.course_name, course_term=course_term, answer_key_blob=f.read())
                db.commit()
            key.delete()
        elif filename.startswith('vote'):
            key.get_contents_to_file(f) #loads the file in our file-like object
            f.seek(0)

            zipped = zipfile.ZipFile(f, 'r')
            votes_file = zipped.open(zipped.infolist()[0].filename)
            zipped.close()
            del zipped
            #expected format of lines in votes.csv ==
            #{presenter username},{voter username},{YYYY-MM-DD},{Overall Score},{Clarity Score},{Difficulty Score}
            for line in votes_file:
                line = [item.strip() for item in line]
                presenter_username = line[0]
                voter_username = line[1]

                #retrieve the corresponding user_id for the student ids
                presenter_row = db(db.user_data.username == presenter_username).select().first()
                voter_row = db(db.user_data.username == voter_username).select().first()

                presenterid = presenter_row.user_id
                projectid = presenter_row.project_id
                voterid = voter_row.user_id
                present_date = line[2]
                overallscore = line[3]
                clarityscore = line[4]
                difficultyscore = line[5]

                #to make it easier to track missing votes, all students have their votes initialized as "Missing"
                #and this is changed to their real vote once it is received, initialization is done the first time
                #a vote is seen for a new presenter
                if db((db.student_votes.presenter_id == presenterid)&(db.student_votes.course_name == course_row.course_name)&(db.student_votes.course_term == course_term)&
                        (db.student_votes.vote_date == present_date)).count() == 0:
                    voter_rows = db((db.user_data.user_id != presenterid)&(db.user_data.course_name == course_row.course_name)&(db.user_data.course_term == course_term)).select()
                    for row in voter_rows:
                        db.student_votes.insert(presenter_id=presenterid,
                                            course_name=course_row.course_name,
                                            course_term=course_term,
                                            project_id=projectid,
                                            voter_id=row.user_id,
                                            vote_date=present_date,
                                            overall_score='Missing',
                                            clarity_score='Missing',
                                            difficulty_score='Missing')

                #update the rows for votes received
                db((db.student_votes.presenter_id == presenterid)&(db.student_votes.course_name == course_row.course_name)&(db.student_votes.course_term == course_term)&
                    (db.student_votes.voter_id == voterid)&
                    (db.student_votes.vote_date == present_date)).update(overall_score=overallscore,
                                                                         clarity_score=clarityscore,
                                                                         difficulty_score=difficultyscore)
                db.commit()

                #No Missing initialization for this section
                #if db((db.student_votes.presenter_id == presenterid)&
                #        (db.student_votes.voter_id == voterid)&
                #        (db.student_votes.vote_date == present_date)).count() > 0:
                #    db((db.student_votes.presenter_id == presenterid)&
                #        (db.student_votes.voter_id == voterid)&
                #        (db.student_votes.vote_date == present_date)).update(overall_score=overallscore,
                #                                                             clarity_score=clarityscore,
                #                                                             difficulty_score=difficultyscore)
                #    db.commit() #update is not done until commit
                #else insert a new record
                #else:
                #    db.student_votes.insert(presenter_id=presenterid,
                #                            project_id=projectid,
                #                            voter_id=voterid,
                #                            vote_date=present_date,
                #                            overall_score=overallscore,
                #                            clarity_score=clarityscore,
                #                            difficulty_score=difficultyscore)
                #    db.commit() #insert is not done until commit

            key.delete()
    '''
    #next process the student bucket
    bucket = conn.get_bucket(os.environ['STUDENT_ANS_BUCKET'])
    rs = bucket.list() #lists the contents of our directory on S3

    #sort the keys by last modified time
    #process each key in order from oldest to newest
    rs = sorted(rs, key=lambda k: k.last_modified)
    for key in rs:
        f = StringIO.StringIO() #gives us a file-like object with which to recieve the file from S3
        file_size = key.size #size of the file in bytes
        filename = key.name #name of the file

        try:
            queueid = int(filename.split('_')[0]) #queue id is the first part of the file name
            userid = int(filename.split('_')[1]) #user id is added to file name after an underscore
        except ValueError:
            key.delete() #deletes the file from s3, done when file name was incorrectly named
            break

        user_row = db(db.user_data.user_id == userid).select().first()
        auth_user_row = db(db.auth_user.id == userid).select().first()

        #log any oversize files that were found and delete them
        if (file_size / 1024 / 1024) > 15:
            #comment = 'Uploaded file of size ' + str(file_size)
            #if the user already has a record in the table, then update it
            #if db(db.user_uploads.user_id == userid).count() > 0:
            #    db(db.user_uploads.user_id == userid).update(comments=comment)
            #    db.commit()
            #else insert a new record
            #else:
            #    db.user_uploads.insert(user_id=userid, comments=comment)
            #    db.commit() #insert is not done until commit
            log_error(userid, 'Uploaded file of size: ' + str(file_size) + ' Exceeded max allowable size.')
            db(db.user_queue.id == queueid).delete() #remove from queue
            db.commit()
            key.delete()
        #for all files of appropriate size, conduct processing as needed
        else:
            key.get_contents_to_file(f) #loads the file in our file-like object
            f.seek(0) #without this reads from this file will not return anything
            ans_dict = None
            student_ans = None
            ans_key = None

            #catch the exception raised by non zip files and log it
            try:
                #files are submitted as zip files
                #1. read the file as a zipfile
                zipped = zipfile.ZipFile(f, 'r')
                #2. use infolist to get the names of the files contained in the zipfile and store them in
                #the db so they can be processed one at a time
                file_list = []
                for i in range(0, len(zipped.infolist())):
                    student_ans = zipped.open(zipped.infolist()[i].filename)
                    first_line = student_ans.readline()
                    if first_line.startswith('>'):
                        genomeid = first_line[1:].strip().lower()
                        if genomeid:
                            filename = genomeid
                            if db((db.upload_files.user_id == userid )&(db.upload_files.file_name==filename)).count() > 0:
                                db((db.upload_files.user_id == userid )&(db.upload_files.file_name==filename)).update(file_blob=student_ans.read())
                                db.commit()
                            else:
                                db.upload_files.insert(user_id=userid, file_name=filename, course_name=user_row.course_name, course_term=user_row.course_term, file_blob=student_ans.read())
                                db.commit()
                            file_list.append(filename)
                zipped.close()
                del f
                del zipped

                if len(file_list) == 0:
                    raise Exception("No answer files were found in the zip file.")

                for next_file in file_list:
                    student_row = db((db.upload_files.user_id == userid )&(db.upload_files.file_name==next_file)).select().first()
                    student_ans = StringIO.StringIO(student_row.file_blob)
                    student_ans.seek(0)

                    #format should be {genome type}_{difficulty}_{genome ID}_chr_{chr ID}
                    #genome_type = next_file.split('_')[0].title() #next_file should be {type}_{difficulty}_{genome id #}_chr_{chromosome #} (e.g. random_E_1_chr_1)
                    genome_type = next_file.split('_')[0]
                    genome_diff = next_file.split('_')[1].lower()
                    genome_id = next_file.split('_')[2]
                    chromosome_id = next_file.split('_')[4]

                    db((db.upload_files.user_id == userid )&(db.upload_files.file_name==next_file)).delete()
                    db.commit()

                    #answer key is in zip format stored in the db
                    #1. retrieve object
                    #2. set as a zip object
                    #3. retrieve file object from zip object
                    row = db((db.genome_table.genome_id == genome_id)&(db.genome_table.chromosome_id == chromosome_id)&(db.genome_table.genome_type == genome_type)&(db.genome_table.difficulty == genome_diff)&
                             (db.genome_table.course_term == user_row.course_term)&(db.genome_table.course_name == user_row.course_name)).select().first()
                    ans_key = StringIO.StringIO(row.answer_key_blob)
                    ans_key.seek(0)
                    zipped = zipfile.ZipFile(ans_key, 'r')
                    ans_key = zipped.open(zipped.infolist()[0].filename)
                    zipped.close()
                    del zipped
                    evaluator = Eval.Eval()

                    ans_dict = evaluator.eval(ans_key, student_ans)
                    #if ans_dict['A_COV'] != -1:
                    #    totalscore = (ans_dict['A_COV'] + ans_dict['A_ACC'] + ans_dict['A_CON']) / 3.0

                        #the following is for debugging purposes, assembly files are being stored so they can be used for later Eval development
                        ##UNCOMMENT TO SAVE FILES TO DEBUGGING BUCKET
                        ##conn2 = S3Connection(os.environ['S3_KEY'], os.environ['S3_SECRET'])
                        ##debugging_row = db(db.user_uploads.user_id == userid).select().first()
                        ##file_count = debugging_row.upload_count
                        ##debugging_bucket = conn2.get_bucket(os.environ['DEBUGGING_BUCKET'])
                        ##debugging_key = Key(debugging_bucket)
                        ##debugging_key.key = str(userid) + '_' + str(file_count)
                        ##student_ans = StringIO.StringIO(student_row.file_blob)
                        ##student_ans.seek(0)
                        ##debugging_key.set_contents_from_file(student_ans)

                    #else:
                    totalscore = 0
                    score_count = 0

                    for ans_dict_key in ans_dict:
                        if ans_dict[ans_dict_key] != -1:
                            score_count +=1
                            totalscore += ans_dict[ans_dict_key]
                        else:
                            ans_dict[ans_dict_key] = 0

                    if score_count > 0:
                        totalscore /= float(score_count)
                    else:
                        totalscore = 0

                    #round all the scores to two decimal points
                    for thiskey in ans_dict:
                        ans_dict[thiskey] = round(float(ans_dict[thiskey])*100, 2)
                    totalscore = round(totalscore*100, 2)

                    u_nickname = user_row.nickname
                    today_in_la = datetime.datetime.now(pytz.timezone('America/Los_Angeles')).date()
                    current_time_in_la = datetime.datetime.now(pytz.timezone('America/Los_Angeles')).time()

                    #keep a record of the last score saved for each day, lets us know what the student had at the end of each day for the assignment
                    if db((db.user_uploads.user_id == userid)&(db.user_uploads.course_term == user_row.course_term)&(db.user_uploads.course_name == user_row.course_name)&(db.user_uploads.genome_type == genome_type)&(db.user_uploads.difficulty == genome_diff)&(db.user_uploads.genome_id == genome_id)&(db.user_uploads.chromosome_id == chromosome_id)&(db.user_uploads.upload_date == today_in_la)).count() > 0:
                        db((db.user_uploads.user_id == userid)&(db.user_uploads.course_term == user_row.course_term)&(db.user_uploads.course_name == user_row.course_name)&(db.user_uploads.genome_type == genome_type)&(db.user_uploads.difficulty == genome_diff)&(db.user_uploads.genome_id == genome_id)&(db.user_uploads.chromosome_id == chromosome_id)&(db.user_uploads.upload_date == today_in_la)).update(user_nickname= u_nickname,
                                                                                                               upload_time=current_time_in_la,
                                                                                                               snp_score= ans_dict['SNP'],
                                                                                                               indel_score= ans_dict['INDEL'],
                                                                                                               copynumber_score= ans_dict['CNV'],
                                                                                                               inversion_score= ans_dict['INV'],
                                                                                                               str_score= ans_dict['STR'],
                                                                                                               alu_score= ans_dict['ALU'],
                                                                                                               assembly_coverage= ans_dict['A_COV'],
                                                                                                               assembly_accuracy= ans_dict['A_ACC'],
                                                                                                               assembly_contig_sizes= ans_dict['A_CON'],
                                                                                                               total_score= totalscore)
                        db.commit()
                    else:
                        db.user_uploads.insert(user_id=userid,
                                               user_nickname= u_nickname,
                                               course_name=user_row.course_name,
                                               course_term=user_row.course_term,
                                               upload_date=today_in_la,
                                               upload_time=current_time_in_la,
                                               genome_id=genome_id,
                                               chromosome_id=chromosome_id,
                                               genome_type=genome_type,
                                               difficulty=genome_diff,
                                               snp_score= ans_dict['SNP'],
                                               indel_score= ans_dict['INDEL'],
                                               copynumber_score= ans_dict['CNV'],
                                               inversion_score= ans_dict['INV'],
                                               str_score= ans_dict['STR'],
                                               alu_score= ans_dict['ALU'],
                                               assembly_coverage= ans_dict['A_COV'],
                                               assembly_contig_sizes= ans_dict['A_CON'],
                                               assembly_accuracy= ans_dict['A_ACC'],
                                               total_score= totalscore)
                        db.commit()

                    #if the student had a previous score, overwrite it
                    if db((db.scores_table.user_id == userid)&(db.scores_table.course_term == user_row.course_term)&(db.scores_table.course_name == user_row.course_name)&(db.scores_table.genome_type == genome_type)&(db.scores_table.difficulty == genome_diff)&(db.scores_table.genome_id == genome_id)&(db.scores_table.chromosome_id == chromosome_id)).count() > 0:
                        db((db.scores_table.user_id == userid)&(db.scores_table.genome_type == genome_type)&(db.scores_table.difficulty == genome_diff)&(db.scores_table.genome_id == genome_id)&(db.scores_table.chromosome_id == chromosome_id)).update(user_nickname= u_nickname,
                                                                                                               snp_score= ans_dict['SNP'],
                                                                                                               indel_score= ans_dict['INDEL'],
                                                                                                               copynumber_score= ans_dict['CNV'],
                                                                                                               inversion_score= ans_dict['INV'],
                                                                                                               str_score= ans_dict['STR'],
                                                                                                               alu_score= ans_dict['ALU'],
                                                                                                               assembly_coverage= ans_dict['A_COV'],
                                                                                                               assembly_accuracy= ans_dict['A_ACC'],
                                                                                                               assembly_contig_sizes= ans_dict['A_CON'],
                                                                                                               total_score= totalscore)
                        db((db.scores_history_table.user_id == userid)&(db.scores_history_table.genome_type == genome_type)&(db.scores_history_table.difficulty == genome_diff)&(db.scores_history_table.genome_id == genome_id)&(db.scores_history_table.chromosome_id == chromosome_id)).update(user_nickname= u_nickname,
                                                                                                               snp_score= ans_dict['SNP'],
                                                                                                               indel_score= ans_dict['INDEL'],
                                                                                                               copynumber_score= ans_dict['CNV'],
                                                                                                               inversion_score= ans_dict['INV'],
                                                                                                               str_score= ans_dict['STR'],
                                                                                                               alu_score= ans_dict['ALU'],
                                                                                                               assembly_coverage= ans_dict['A_COV'],
                                                                                                               assembly_accuracy= ans_dict['A_ACC'],
                                                                                                               assembly_contig_sizes= ans_dict['A_CON'],
                                                                                                               total_score= totalscore)
                        db.commit() #update is not done until commit
                    #else insert a new record
                    else:
                        db.scores_table.insert(user_id=userid,
                                               user_nickname= u_nickname,
                                               course_name=user_row.course_name,
                                               course_term=user_row.course_term,
                                               genome_id=genome_id,
                                               chromosome_id=chromosome_id,
                                               difficulty=genome_diff,
                                               genome_type=genome_type,
                                               snp_score= ans_dict['SNP'],
                                               indel_score= ans_dict['INDEL'],
                                               copynumber_score= ans_dict['CNV'],
                                               inversion_score= ans_dict['INV'],
                                               str_score= ans_dict['STR'],
                                               alu_score= ans_dict['ALU'],
                                               assembly_coverage= ans_dict['A_COV'],
                                               assembly_contig_sizes= ans_dict['A_CON'],
                                               assembly_accuracy= ans_dict['A_ACC'],
                                               total_score= totalscore)
                        db.scores_history_table.insert(user_id=userid,
                                               user_nickname= u_nickname,
                                               user_firstname= auth_user_row.first_name,
                                               user_lastname= auth_user_row.last_name,
                                               user_email= auth_user_row.email,
                                               course_name=user_row.course_name,
                                               course_term=user_row.course_term,
                                               genome_id=genome_id,
                                               chromosome_id=chromosome_id,
                                               difficulty=genome_diff,
                                               genome_type=genome_type,
                                               snp_score= ans_dict['SNP'],
                                               indel_score= ans_dict['INDEL'],
                                               copynumber_score= ans_dict['CNV'],
                                               inversion_score= ans_dict['INV'],
                                               str_score= ans_dict['STR'],
                                               alu_score= ans_dict['ALU'],
                                               assembly_coverage= ans_dict['A_COV'],
                                               assembly_contig_sizes= ans_dict['A_CON'],
                                               assembly_accuracy= ans_dict['A_ACC'],
                                               total_score= totalscore)
                        db.commit() #insert is not done until commit
                    #db.upload_files.insert(user_id='test',
                    #                            file_name='test',
                    #                            file_blob=f)
                    #db.upload_files.insert(user_id=userid,
                    #                        file_name=filename,
                    #                        file_blob=gzip.GzipFile(fileobj=f, mode='rb'))
                    #db.commit() #insert is not done until commit
                    db(db.user_queue.id <= queueid).delete() #clear the queue up to this point
                    db.commit()

                    if not student_ans.closed:
                        student_ans.close()
                    if not ans_key.closed:
                        ans_key.close()
                    key.delete() #deletes the file from S3

            except Exception as detail:
                db.rollback()
                exc_type, exc_value, exc_traceback = sys.exc_info()
                if "object has no attribute 'answer_key_blob" in str(detail):
                    detail = 'Answer key was not found, check that you used the correct genome name in your description line'
                #log exception and delete file/queue entry
                log_error(userid, 'An error occured while evaluating your file.' +
                              ' Please check your file to ensure that its format matches the format given in the documentation.' +
                              ' The exception was: ' + str(detail))
                traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stderr)
                db(db.user_queue.id == queueid).delete() #remove from queue
                db(db.scheduler_task.status == "COMPLETED").delete() #remove the old entries from the schedulers db so that it doesn't create too many rows
                db.commit()
                key.delete()
                if student_ans and not student_ans.closed:
                    student_ans.close()
                if ans_key and not ans_key.closed:
                    ans_key.close()

# scheduler is used to assign tasks to background processes
from gluon.scheduler import Scheduler
scheduler = Scheduler(db, dict(process_file=process_file), discard_results=True)