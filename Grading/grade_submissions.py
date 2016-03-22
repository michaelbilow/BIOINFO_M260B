import pandas as pd
from os import listdir
from os.path import split, join, splitext
import datetime
import sys

GRADED_GENOMES = {'hw1': 'PA1', 'hw2grad': 'PA2grad',
                  'hw3all': 'PA3', 'hw2undergrad': 'PA2ug'}


CUTOFF_DATES_W16 = {'PA1': datetime.date(2016, 1, 15), 'PA2': datetime.date(2016, 2, 8),
                    'PA3': datetime.date(2016, 2, 26), 'Final': datetime.date(2016, 3, 16)}

SCORE_CUTOFFS = {'PA1': {'snp': {'grad': 65, 'ug': 55, 'floor': 40}},
                 'PA2': {'snp': {'grad': 65, 'ug': 55, 'floor': 40},
                         'indel': {'grad': 15, 'ug': 10, 'floor': 0}},
                 'PA3': {'coverage': {'grad': 90, 'ug': 85, 'floor': 75},
                         'accuracy': {'grad': 49, 'ug': 44, 'floor': 34},
                         'contig_sizes': {'grad': 4, 'ug': 2, 'floor': 0}},
                 'Final': {'snp': {'grad': 90, 'ug': 75, 'floor': 55},
                           'indel': {'grad': 50, 'ug': 30, 'floor': 15},
                           'str': {'grad': 90, 'ug': 65, 'floor': 40},
                           'inversion': {'grad': 20, 'ug': 10, 'floor': 0},
                           'copynumber': {'grad': 20, 'ug': 10, 'floor': 0}}
                 }

PROJECT_WEIGHTS = {'PA1': [],
                   'PA2': [{'snp': 70, 'indel': 30}],
                   'PA3': [],
                   'Final': [{'snp': 70, 'indel': 30}, {'str': 100}, {'inversion': 100}, {'copynumber': 100}]}


def score_df(project, score_df):
    cutoff_dict = SCORE_CUTOFFS[project]
    weights_dict = PROJECT_WEIGHTS[project]
    cutoff_date = CUTOFF_DATES_W16[project]

    for ix in score_df.index:
        row = score_df.ix[ix]
        raw_score_dict = score_row(row, cutoff_dict, weights_dict, cutoff_date)
        late_days = (row['upload_date'] - cutoff_date).days
        max_row_score = compute_max_score(raw_score_dict, weights_dict)


def score_row(row, cutoff_dict, weights_dict, cutoff_date):
    print row
    print cutoff_dict
    row_score_dict = {}
    for k in cutoff_dict:
        component_cutoff_dict = cutoff_dict[k]
        floor = component_cutoff_dict['floor']
        ug_max = component_cutoff_dict['ug']
        grad_max = component_cutoff_dict['grad']
        raw_score = row[k]
        ug_score = min(100.0, max(100 * float((raw_score - ug_max))/(ug_max - floor), 0.0))
        grad_score = min(100.0, max(100 * float((raw_score - grad_max))/(grad_max - floor), 0.0))
        row_score_dict[k] = {'ug': ug_score, 'grad': grad_score}
        print row_score_dict
    return row_score_dict



if __name__ == "__main__":
    input_folder = './submissions'
    input_fn = join(input_folder, listdir(input_folder)[0])
    df = pd.read_excel(input_fn)
    df.columns = [_.split('.')[1] for _ in df.columns]

    df['upload_date'] = pd.to_datetime(df['upload_date'])

    score_cols = [_ for _ in df.columns if any(_1 in _ for _1 in ('score', 'assembly'))]
    df = df.ix[df['genome_type'].isin(GRADED_GENOMES)]
    score_types = [GRADED_GENOMES[_] for _ in df['genome_type']]
    df['project'] = score_types
    df = df.ix[:, score_cols + ['user_id', 'project', 'upload_date']]
    df.columns = [_.split('_')[0] if _.split('_')[-1] == 'score' else '_'.join(_.split('_')[1:]) if _.split('_')[0] == 'assembly' else _ for _ in df.columns]


    grouped = df.groupby(['user_id', 'project'])

    output_scores = []
    for k, g_df in grouped:
        # print k
        # print g_df
        studentid, project = k
        if project.startswith('PA2'):
            project_score = score_df('PA2', g_df)
            output_scores.append(project_score)
            final_project_score = score_df('Final', g_df)
            output_scores.append(final_project_score)

    # print df.columns
    # print df.head()
    # print set(df['genome_id'])
    # print set(df['difficulty'])
    # print set(df['chromosome_id'])
    # print set(df['genome_type'])
    # grouped = df.groupby()

