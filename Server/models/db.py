# -*- coding: utf-8 -*-

#########################################################################
## This scaffolding model makes your app work on Google App Engine too
## File is released under public domain and you can use without limitations
#########################################################################
## if SSL/HTTPS is properly configured and you want all HTTP requests to
## be redirected to HTTPS, uncomment the line below:
import os, shutil

# request.requires_https()
if not request.env.web2py_runtime_gae:
    ## if NOT running on Google App Engine use SQLite or other DB

    # HEROKU_POSTGRESQL_SOMECOLOR_URL should have our db color as SOMECOLOR,

    # The following code is a workaround for an error we encountered when using web2py on Heroku. Heroku doesn't let you write to the filesystem, which is the default way
    # in which web2py tracks the tables that exist in the database. So if migrate is enabled, then web2py will want to create any table that has not already been created.
    # But it only knows whether tables have been created based on what should be stored in the file system. It looks for the file that tells it what has been created,
    # finds nothing, and assumes that the tables have not been created yet. But when it issues the command to the DB to create the table, the DB finds that the table
    # already exists and throws an error, crashing the app. The following work around allows us to have tables created automatically when we first launch:
    # on initial launch: remove migrate_enabled=False from the following command, and git push, then immediately after the app has launched, replace migrate_enabled=False, and git push again
    # by doing so, the app launches and creates all the databases, then the second launch keeps it from crashing when it finds that tables have already been created. If you forget to replace
    # migrate_enabled=True, then the app will work until the first time it resets (git push or web thread resets), and then crash indefinitely until you replace the line or wipe the postgres database
    # After the initial launch, you need to manually add any additional tables and manually alter any tables when you make changes in this file (db.py).
    # All manual db work is done through the command line using heroku pg:psql (note: you will need to install postgres on your own system before this will work)
    db = DAL(os.environ['HEROKU_POSTGRESQL_ONYX_URL'], pool_size=10, check_reserved=['all'], migrate_enabled=False)
    # db = DAL('sqlite://storage.sqlite',pool_size=1,check_reserved=['all'])
else:
    ## connect to Google BigTable (optional 'google:datastore://namespace')

    db = DAL('google:datastore+ndb')
    ## store sessions and tickets there
    session.connect(request, response, db=db)
    ## or store session in Memcache, Redis, etc.
    ## from gluon.contrib.memdb import MEMDB
    ## from google.appengine.api.memcache import Client
    ## session.connect(request, response, db = MEMDB(Client()))

## by default give a view/generic.extension to all actions from localhost
## none otherwise. a pattern can be 'controller/function.extension'
response.generic_patterns = ['*'] if request.is_local else []
## (optional) optimize handling of static files
# response.optimize_css = 'concat,minify,inline'
# response.optimize_js = 'concat,minify,inline'
## (optional) static assets folder versioning
# response.static_version = '0.0.0'
#########################################################################
## Here is sample code if you need for
## - email capabilities
## - authentication (registration, login, logout, ... )
## - authorization (role based authorization)
## - services (xml, csv, json, xmlrpc, jsonrpc, amf, rss)
## - old style crud actions
## (more options discussed in gluon/tools.py)
#########################################################################

from gluon.tools import Auth, Crud, Service, PluginManager, prettydate

auth = Auth(db, secure=True)  # adding secure=True forces HTTPS
crud, service, plugins = Crud(db), Service(), PluginManager()

## create all tables needed by auth if not custom tables
auth.define_tables(username=True, signature=False)

## configure email
mail = auth.settings.mailer
mail.settings.server = 'smtp.sendgrid.net:587'
mail.settings.sender = 'cs122projects2016@gmail.com'
mail.settings.login = os.environ['SENDGRID_USERNAME'] + ':' + os.environ['SENDGRID_PASSWORD']

# This email will be sent to the admin's email address allowing them to verify new users
import sendgrid
import uuid


# this custom method is necessary because web2py only supports sending the email back to the user for them
# to verify it, it doesn't support sending an email to the admin notifying them that a user needs to be
# approved
def send_verify_message(form):
    # store a uuid with the record to match requests in emails to entries in the db
    key = uuid.uuid4()
    db(db.auth_user.username == form.vars.username).update(registration_id=key)
    db.commit()

    sg = sendgrid.SendGridClient(os.environ['SENDGRID_USERNAME'], os.environ['SENDGRID_PASSWORD'])
    verification_message = sendgrid.Mail()
    verification_message.add_to('cs122projects2016@gmail.com')
    verification_message.set_subject('Verify new user')
    verification_message.set_text('A new user has requested an account on the course Heroku app. ' +
                                  ' Their information is as follows:\n\n' +
                                  'First name: ' + form.vars.first_name + '\n' +
                                  'Last name: ' + form.vars.last_name + '\n' +
                                  'User name: ' + form.vars.username + '\n' +
                                  'Email: ' + form.vars.email + '\n\n' +
                                  'Approve this account with the following link: https://' + request.env.http_host +
                                  URL(r=request, c='default', f='approve_user') + '/' + str(key) + '\n\n\n\n' +
                                  'Reject this account with the following link: https://' + request.env.http_host +
                                  URL(r=request, c='default', f='reject_user') + '/' + str(key))
    verification_message.set_from('app24179024@heroku.com')
    sg.send(verification_message)


# auth.settings.register_onaccept.append(lambda form:
#    mail.send(to='ad...@emai.com',subject='New user registered for %s application' % (request.application),
#                message="new user email is %s" % (form.vars.email)))

## configure auth policy
auth.settings.registration_requires_verification = False
auth.settings.registration_requires_approval = True
auth.settings.reset_password_requires_verification = True
auth.settings.login_after_registration = False
auth.settings.register_onaccept = lambda form: send_verify_message(form)
# form = crud.create(db.comment, onaccept=lambda form,xyz=xyz:funcdone(form,xyz))

## turn off the user's ability to register and modify their account, auto created when logging on through UCLA email
auth.settings.actions_disabled = ['profile']  # 'register','change_password','request_reset_password','profile']

# Recaptcha doesn't work for us since we're using a shared SSL domain
# from gluon.tools import Recaptcha
# auth.settings.captcha = Recaptcha(request,
#                                  public_key='6LcBnvkSAAAAAAG8SbrAPaCRlKFdiWtN-aETFczY',
#                                  private_key='6LcBnvkSAAAAAMWvhCCjX5kUIEoOLIOIk1Q5p2CB')

## sets login authentication to use UCLA's SMTP (email) server
# from gluon.contrib.login_methods.email_auth import email_auth
# auth.settings.login_methods.append(email_auth("mail.ucla.edu:587", "@ucla.edu"))
# auth.settings.login_methods.append(email_auth("smtp.cs.ucla.edu:587", "@cs.ucla.edu"))
# auth.settings.login_methods.append(email_auth("smtp.gmail.com:587", "@gmail.com")) #this line allows anyone with a gmail account access

## if you need to use OpenID, Facebook, MySpace, Twitter, Linkedin, etc.
## register with janrain.com, write your domain:api_key in private/janrain.key
# from gluon.contrib.login_methods.janrain_account import use_janrain
# use_janrain(auth, filename='private/janrain.key')

#########################################################################
## Define your tables below (or better in another model file) for example
##
## >>> db.define_table('mytable',Field('myfield','string'))
##
## Fields can be 'string','text','password','integer','double','boolean'
##       'date','time','datetime','blob','upload', 'reference TABLENAME'
## There is an implicit 'id integer autoincrement' field
## Consult manual for more options, validators, etc.
##
## More API examples for controllers:
##
## >>> db.mytable.insert(myfield='value')
## >>> rows=db(db.mytable.myfield=='value').select(db.mytable.ALL)
## >>> for row in rows: print row.id, row.myfield
#########################################################################

# Used to track information about the different private genomes
db.define_table('genome_table',
                Field('genome_id', 'string'),
                Field('chromosome_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('answer_key_blob', 'blob'),
                Field('difficulty', 'string'),
                Field('genome_type', 'string'),
                #    Field('assignment_name', 'string'),
                #    Field('assignment_type', 'string'), #normal or assembly
                #    Field('is_practice', 'boolean'), #if True, then assignment submissions will not be backed up to the score_archive table or the user_uploads table
                #    Field('download_link', 'string'), #link to zip file containing data needed by the students
                primarykey=['genome_type', 'difficulty', 'genome_id', 'chromosome_id', 'course_name', 'course_term'])

db.define_table('assignment_table',
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('assignment_name', 'string'),
                Field('assignment_number', 'integer'),
                Field('deadline_date', 'date'),
                Field('deadline_time', 'time'),
                Field('graded_assignment', 'boolean'),
                Field('assembly', 'boolean'),
                Field('undergrad_instructions', 'text'),
                Field('undergrad_genome_id', 'string'),
                Field('undergrad_chromosome_id', 'string'),
                Field('grad_instructions', 'text'),
                Field('grad_genome_id', 'string'),
                Field('grad_chromosome_id', 'string'),
                primarykey=['course_name', 'course_term', 'assignment_name', 'assignment_number'])

# Used to track the student files that have been loaded onto S3 for each assignment
db.define_table('assignment_files_table',
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('assignment_name', 'string'),
                Field('assignment_number', 'integer'),
                Field('file_name', 'string'),
                Field('file_link', 'string'),
                primarykey=['course_name', 'course_term', 'assignment_name', 'assignment_number', 'file_name'])

# Used to track user scores for the different genomes
db.define_table('scores_table',
                Field('user_id', 'string'),
                Field('user_nickname', 'string'),
                Field('genome_id', 'string'),
                Field('chromosome_id', 'string'),
                Field('difficulty', 'string'),
                Field('genome_type', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('snp_score', 'double', default=0),
                Field('indel_score', 'double', default=0),  # based on their answers for both inserts and deletes
                Field('copynumber_score', 'double', default=0),
                Field('inversion_score', 'double', default=0),
                Field('str_score', 'double', default=0),
                Field('alu_score', 'double', default=0),
                Field('total_score', 'double', default=0),
                Field('assembly_coverage', 'double', default=0),
                Field('assembly_accuracy', 'double', default=0),
                Field('assembly_contig_sizes', 'double', default=0),
                primarykey=['user_id', 'course_name', 'course_term', 'genome_id', 'chromosome_id', 'difficulty',
                            'genome_type'])
# ALTER TABLE scores_table ADD COLUMN str_score double precision;
# ALTER TABLE scores_table ADD COLUMN alu_score double precision;
# UPDATE scores_table SET alu_score = 0;

# Used to keep a long term history across several classes
db.define_table('scores_history_table',
                Field('user_id', 'string'),
                Field('user_nickname', 'string'),
                Field('user_firstname', 'string'),
                Field('user_lastname', 'string'),
                Field('user_email', 'string'),
                Field('genome_id', 'string'),
                Field('chromosome_id', 'string'),
                Field('difficulty', 'string'),
                Field('genome_type', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('snp_score', 'double', default=0),
                Field('indel_score', 'double', default=0),  # based on their answers for both inserts and deletes
                Field('copynumber_score', 'double', default=0),
                Field('inversion_score', 'double', default=0),
                Field('str_score', 'double', default=0),
                Field('alu_score', 'double', default=0),
                Field('total_score', 'double', default=0),
                Field('assembly_coverage', 'double', default=0),
                Field('assembly_accuracy', 'double', default=0),
                Field('assembly_contig_sizes', 'double', default=0),
                primarykey=['user_id', 'course_name', 'course_term', 'genome_id', 'chromosome_id', 'difficulty',
                            'genome_type'])

# from the command line, you can manually create a new table using the following syntax
# heroku pg:psql
# DROP TABLE scores_table;
# CREATE SEQUENCE scores_table_seq;
# CREATE TABLE scores_table (
#    _id integer NOT NULL DEFAULT nextval('scores_table_seq'),
#    user_id varchar(512),
#    genome_id varchar(512),
#    snp_score double precision,
#    indel_score double precision,
#    copynumber_score double precision,
#    inversion_score double precision,
#    str_score double precision,
#    total_score double precision,
#    CONSTRAINT scores_key PRIMARY KEY(user_id, genome_id)
# );
# ALTER SEQUENCE scores_table_seq OWNED BY scores_table._id;


# Used to pass uploaded files between worker processes
# The writable=False, readable=False just keep these fields from showing up in SQLFORM forms
# The uploads field stores the files in a temp directory and tracks there name and location
# The blob field isn't necessary now, but put it there in case we later choose to store files in the db
db.define_table('upload_files',
                Field('user_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('file_name', 'string'),
                Field('file_blob', 'blob'))

# Used to track user uploads, allows us to check what score student had prior to the deadline
db.define_table('user_uploads',
                Field('user_id', 'string'),
                Field('user_nickname', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('upload_date', 'date'),
                Field('upload_time', 'time'),
                Field('genome_id', 'string', default='None'),
                Field('chromosome_id', 'string', default='None'),
                Field('difficulty', 'string', default='None'),
                Field('genome_type', 'string', default='None'),
                Field('snp_score', 'double', default=0),
                Field('indel_score', 'double', default=0),
                Field('copynumber_score', 'double', default=0),
                Field('inversion_score', 'double', default=0),
                Field('str_score', 'double', default=0),
                Field('alu_score', 'double', default=0),
                Field('total_score', 'double', default=0),
                Field('assembly_coverage', 'double', default=0),
                Field('assembly_accuracy', 'double', default=0),
                Field('assembly_contig_sizes', 'double', default=0))

# From the command line, you can manually alter a table using the following syntax
# heroku pg:psql
# ALTER TABLE user_uploads ADD COLUMN last_upload_size integer;
# ALTER TABLE user_uploads ADD COLUMN total_mb_uploaded integer;

# lists the overall data for each course, allows us to only display results from students in the current course and to easily delete data for old classes
db.define_table('course_data',
                Field('course_name', 'string', requires=IS_NOT_EMPTY(), default='Algorithms in Bioinformatics'),
                Field('course_quarter', requires=IS_IN_SET(['Fall', 'Winter', 'Spring', 'Summer'])),
                Field('course_year', requires=IS_IN_SET(
                    ['2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024', '2025', '2026',
                     '2027', '2028', '2029', '2030'])),
                Field('current_course', 'boolean', default=True),
                # set to true when you want a courses entries displayed in the scoreboard
                primarykey=['course_name', 'course_quarter', 'course_year'])

db.define_table('project_data',
                Field('project_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('project_name', 'string'),
                Field('comments', 'string'),
                primarykey=['project_id', 'course_name', 'course_term'])

db.define_table('user_data',
                Field('user_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('username', 'string', writable=False, readable=False),
                Field('nickname', 'string'),
                Field('presentation', 'string', default='Not Uploaded'),
                Field('presentation_count', 'integer', default=0),
                Field('project_id', 'string', default=''),
                Field('project_name', 'string', default='No Project Selected', writable=False),
                primarykey=['user_id', 'course_name', 'course_term'])

db.define_table('student_votes',
                Field('presenter_username', 'string'),
                Field('presenter_id', 'string'),
                Field('project_id', 'string'),
                Field('voter_username', 'string'),
                Field('voter_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('vote_date', 'string'),  # no date processing is needed so just storing the string representation
                Field('overall_score', 'string'),
                # using strings so that these can be recorded as "Not Submitted" as well
                Field('clarity_score', 'string'),
                Field('difficulty_score', 'string'))

# CREATE SEQUENCE student_votes_seq;
# CREATE TABLE student_votes (
#    _id integer NOT NULL DEFAULT nextval('student_votes_seq'),
#    presenter_id varchar(512),
#    project_id varchar(512),
#    voter_id varchar(512),
#    vote_datetime timestamp,
#    overall_score integer,
#    clarity_score integer,
#    difficulty_score integer,
# );
# ALTER SEQUENCE student_votes_seq OWNED BY student_votes._id;

# DROP TABLE user_nickname;
# CREATE SEQUENCE user_data_seq;
# CREATE TABLE user_data (
#    _id integer NOT NULL DEFAULT nextval('user_data_seq'),
#    user_id varchar(512),
#    nickname varchar(512),
#    project varchar(512),
#    CONSTRAINT user_data_key PRIMARY KEY(user_id)
# );
# ALTER SEQUENCE user_data_seq OWNED BY user_data._id;

# db.define_table('user_nickname',
#    Field('user_id', 'string'),
#    Field('nickname', 'string'),
#    primarykey=['user_id'])

# heroku pg:psql
# DROP TABLE user_nickname;
# CREATE SEQUENCE user_nickname_seq;
# CREATE TABLE user_nickname (
#    _id integer NOT NULL DEFAULT nextval('user_nickname_seq'),
#    user_id varchar(512),
#    nickname varchar(512),
#    CONSTRAINT user_nickname_key PRIMARY KEY(user_id)
# );
# ALTER SEQUENCE user_nickname_seq OWNED BY user_nickname._id;

db.define_table('user_queue',
                Field('user_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('upload_time', 'datetime'))

# heroku pg:psql
# DROP TABLE user_queue;
# CREATE SEQUENCE user_queue_seq;
# CREATE TABLE user_queue (
#    id integer NOT NULL DEFAULT nextval('user_queue_seq'),
#    user_id varchar(512),
#    upload_time timestamp
# );
# ALTER SEQUENCE user_queue_seq OWNED BY user_queue.id;

# db.define_table('sanity_checks',
#                Field('user_id', 'string'),
#                Field('last_score', 'double'),
#                Field('results_blob', 'blob'),
#                primarykey=['user_id'])

# heroku pg:psql
# DROP TABLE sanity_checks;
# CREATE SEQUENCE sanity_checks_seq;
# CREATE TABLE sanity_checks (
#    _id integer NOT NULL DEFAULT nextval('sanity_checks_seq'),
#    user_id varchar(512),
#    last_score double precision,
#    results_blob bytea,
#    CONSTRAINT sanity_checks_key PRIMARY KEY(user_id)
# );
# ALTER SEQUENCE sanity_checks_seq OWNED BY sanity_checks._id;

db.define_table('error_reports',
                Field('user_id', 'string'),
                Field('course_name', 'string'),
                Field('course_term', 'string'),
                Field('error', 'string'),
                primarykey=['user_id', 'course_name', 'course_term'])


# heroku pg:psql
# DROP TABLE error_reports;
# CREATE SEQUENCE error_reports_seq;
# CREATE TABLE error_reports (
#    _id integer NOT NULL DEFAULT nextval('error_reports_seq'),
#    user_id varchar(512),
#    error varchar(512),
#    CONSTRAINT error_reports_key PRIMARY KEY(user_id)
# );
# ALTER SEQUENCE error_reports_seq OWNED BY error_reports._id;

## after defining tables, uncomment below to enable auditing
# auth.enable_record_versioning(db)

def current_course_name():
    row = db(db.course_data.current_course == True).select().first()
    if row and row.course_name != "":
        return row.course_name
    else:
        return "Bioinformatics"