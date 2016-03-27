import pandas as pd
from os.path import join
import numpy as np


if __name__ == "__main__":
    grades_df = pd.read_excel('output_test.xlsx')
    user_df = pd.read_excel(join('user_data', 'W16.xlsx'))
    user_df = user_df.ix[~np.isnan(user_df['UID'])]
    user_df['UID'] = ['{0:09d}'.format(int(_)) for _ in user_df['UID']]
    print grades_df.head()
    print user_df.head()
    merged_df = pd.merge(grades_df, user_df, left_on=['studentid'], right_on=['user_data.user_id'])
    print merged_df.head()
    merged_df.set_index(['Grad', 'Last Name', 'First Name', 'UID', 'DATASET', 'PROJECT'], inplace=True)
    merged_df.sort_index(inplace=True)
    merged_df.to_excel('homework_grades.xlsx')