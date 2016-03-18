import pandas as pd
from os import listdir
from os.path import split, join, splitext
import datetime


CUTOFF_DATES_W16 = {'PA1': datetime.date(2016, 1, 15), 'PA2': datetime.date(2016, 2, 8),
                    'PA3': datetime.date(2016, 2, 26), 'Final': datetime.date(2016, 3, 16)}

SCORE_CUTOFFS = {'PA1': {'snp': {'grad': 65, 'ug': 55, 'floor': 40}},
                 'PA2': {'snp': {'grad': 65, 'ug': 55, 'floor': 40},
                         'indel': {'grad': 15, 'ug': 10, 'floor': 0}},
                 'PA3': {'coverage': {'grad': 90, 'ug': 85, 'floor': 75},
                         'accuracy': {'grad': 49, 'ug': 44, 'floor': 34},
                         'contig_sizes': {'grad': 4, 'ug': 2, 'floor': 0}},
                 'Final': {'snp': {'grad': },
                           'str': {'grad': 90, 'ug': }}
                 }

PROJECT_WEIGHTS = {'PA1': [],
                   'PA2': [{'snp': 70, 'indel': 30}],
                   'PA3': []
                   'Final': [{'snp': 70, 'indel': 30}, {'str': 100}, {'inversion':100}, {'copynumber':100}]}


def PA1_SCORE_FORMULA():

    return

def PA2_SCORE_FORMULA():
    return

def PA3_SCORE_FORMULA():
    return

def FINAL_PROJECT_SCORE_FORMULA():

    return

if __name__ == "__main__":
    input_folder = './submissions'
    input_fn = join(input_folder, listdir(input_folder)[0])
    df = pd.read_excel(input_fn)
    df.columns = [_.split('.')[1] for _ in df.columns]

    df['upload_date'] = pd.to_datetime(df['upload_date'])

    print df.columns
    print df.head()