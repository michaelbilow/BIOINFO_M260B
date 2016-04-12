import pandas as pd
from os import listdir
from os.path import split, join, splitext
import datetime
import sys
import numpy as np

GRADED_GENOMES = {'hw1': 'PA1', 'hw2grad': 'PA2grad',
                  'hw3all': 'PA3', 'hw2undergrad': 'PA2ug'}

CUTOFF_DATES_W16 = {'PA1': pd.Timestamp('2016-01-15'), 'PA2': pd.Timestamp('2016-02-08'),
                    'PA3': pd.Timestamp('2016-02-26'), 'Final': pd.Timestamp('2016-03-16')}

SCORE_CUTOFFS = {'PA1': {'snp': {'grad': 65, 'ug': 55, 'floor': 40}},
                 'PA2': {'snp': {'grad': 65, 'ug': 55, 'floor': 40},
                         'indel': {'grad': 15, 'ug': 10, 'floor': 0}},
                 'PA3': {'coverage': {'grad': 90, 'ug': 85, 'floor': 75},
                         'accuracy': {'grad': 49, 'ug': 44, 'floor': 34},
                         'contig_sizes': {'grad': 4, 'ug': 2, 'floor': 0}},
                 'Final': {'snp': {'grad': 90, 'ug': 75, 'floor': 55},
                           'indel': {'grad': 50, 'ug': 30, 'floor': 15},
                           'str': {'grad': 70, 'ug': 65, 'floor': 40},  ## For future, raise the grad standard up to 80.
                           'inversion': {'grad': 20, 'ug': 20, 'floor': 0},
                           'copynumber': {'grad': 60, 'ug': 10, 'floor': 0}}
                 }

PROJECT_WEIGHTS = {'PA1': [],
                   'PA2': [{'snp': 70, 'indel': 30}],
                   'PA3': [],
                   'Final': [{'snp': 70, 'indel': 30}, {'str': 100}, {'inversion': 100}, {'copynumber': 100}]}


def score_df(project, score_df):
    cutoff_dict = SCORE_CUTOFFS[project]
    weights_dict = PROJECT_WEIGHTS[project]
    cutoff_date = CUTOFF_DATES_W16[project]

    ug_scores = []
    grad_scores = []
    late_days = []

    for ix in score_df.index:
        row = score_df.ix[ix]
        raw_score_dict = score_row(row, cutoff_dict)
        row_late_days = max((row['upload_date'] - cutoff_date).days, 0)
        max_row_score_dict = compute_max_score(raw_score_dict, weights_dict)
        ug_scores.append(max_row_score_dict['ug'])
        grad_scores.append(max_row_score_dict['grad'])
        late_days.append(row_late_days)

    score_df.ix[:, 'UG_RAW_SCORE'] = ug_scores
    score_df.ix[:, 'GRAD_RAW_SCORE'] = grad_scores
    score_df.ix[:, 'LATE_DAYS'] = late_days
    score_df.ix[:, 'UG_LATE_ADJUSTED'] = score_df.ix[:, 'UG_RAW_SCORE'] - 3 * score_df.ix[:, 'LATE_DAYS']
    score_df.ix[:, 'GRAD_LATE_ADJUSTED'] = score_df.ix[:, 'GRAD_RAW_SCORE'] - 3 * score_df.ix[:, 'LATE_DAYS']

    output_dict = {'{}_{}'.format(grade_level, score_type):
                       score_df['{}_{}'.format(grade_level, score_type)].max()
                   for grade_level in ('UG', 'GRAD') for score_type in ('RAW_SCORE', 'LATE_ADJUSTED')}
    output_dict = {k: v if v > 0 else 0 for k, v in output_dict.items()}

    output_dict['UG_LATE_DAYS'] = score_df.ix[score_df['UG_LATE_ADJUSTED'].idxmax(), 'LATE_DAYS']
    output_dict['GRAD_LATE_DAYS'] = score_df.ix[score_df['GRAD_LATE_ADJUSTED'].idxmax(), 'LATE_DAYS']
    output_dict['PROJECT'] = project
    assert len(set(score_df['project'])) == 1
    output_dict['DATASET'] = score_df['project'].min()
    # print output_dict
    # print score_df
    return output_dict


def score_row(row, cutoff_dict):
    # print row
    # print cutoff_dict
    row_score_dict = {}
    print project
    for k in cutoff_dict:
        component_cutoff_dict = cutoff_dict[k]
        floor = component_cutoff_dict['floor']
        ug_max = component_cutoff_dict['ug']
        grad_max = component_cutoff_dict['grad']
        raw_score = row[k]
        print k, raw_score,
        ug_score = min(100.0, max(100 * float((raw_score - floor)) / (ug_max - floor), 0.0))
        grad_score = min(100.0, max(100 * float((raw_score - floor)) / (grad_max - floor), 0.0))
        row_score_dict[k] = {'ug': ug_score, 'grad': grad_score}
    print row_score_dict
    return row_score_dict


def compute_max_score(scores_dict, weights_dict_list):
    score_lambda = lambda score_comps, grade_level, weights=None: \
        np.average([scores_dict[score_component][grade_level] for score_component in score_comps], weights=weights)

    if not weights_dict_list:
        ## Take the average of all the scored components
        output_scores = {grade_level:
                             score_lambda(scores_dict.keys(), grade_level)
                         for grade_level in ('ug', 'grad')}
        print output_scores
        return output_scores

    potential_scores = []
    for weights_dict in weights_dict_list:
        these_scores = {grade_level: score_lambda(weights_dict.keys(), grade_level, weights_dict.values())
                        for grade_level in ('ug', 'grad')}
        potential_scores.append(these_scores)

    output_scores = {grade_level: max([_[grade_level] for _ in potential_scores])
                     for grade_level in ('ug', 'grad')}
    print output_scores
    return output_scores


if __name__ == "__main__":
    input_folder = './submissions'
    input_fn = join(input_folder, listdir(input_folder)[0])
    df = pd.read_csv(input_fn)
    df.columns = [_.split('.')[1] for _ in df.columns]

    df['upload_date'] = pd.to_datetime(df['upload_date'])

    score_cols = [_ for _ in df.columns if any(_1 in _ for _1 in ('score', 'assembly'))]
    df = df.ix[df['genome_type'].isin(GRADED_GENOMES)]
    score_types = [GRADED_GENOMES[_] for _ in df['genome_type']]
    df['project'] = score_types
    df = df.ix[:, score_cols + ['user_id', 'project', 'upload_date']]
    df.columns = [_.split('_')[0] if _.split('_')[-1] == 'score' else '_'.join(_.split('_')[1:]) if _.split('_')[
                                                                                                        0] == 'assembly' else _
                  for _ in df.columns]

    grouped = df.groupby(['user_id', 'project'])

    output_scores = []
    for k, g_df in grouped:
        # print k
        # print g_df
        studentid, project = k
        if project.startswith('PA2'):
            thisproject = 'PA2'
            project_score = score_df('PA2', g_df)
            project_score['studentid'] = studentid
            print project_score
            output_scores.append(project_score)
            thisproject = 'Final'
            final_project_score = score_df('Final', g_df)
            final_project_score['studentid'] = studentid
            print final_project_score
            output_scores.append(final_project_score)
        else:
            project_score = score_df(project, g_df)
            project_score['studentid'] = studentid
            print project_score
            output_scores.append(project_score)

    output_df = pd.DataFrame.from_records(output_scores)
    output_df.to_excel('output_test.xlsx')

    # print df.columns
    # print df.head()
    # print set(df['genome_id'])
    # print set(df['difficulty'])
    # print set(df['chromosome_id'])
    # print set(df['genome_type'])
    # grouped = df.groupby()
