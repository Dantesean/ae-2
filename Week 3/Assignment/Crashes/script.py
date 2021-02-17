import pandas as pd
import numpy as np
import glob
import tqdm

frame = pd.read_csv('data/large/frame.csv')

df1 = frame
df2 = pd.read_csv('data/large/intermmediatecsv.csv')
df1['point'] = [(x, y) for x,y in zip(df1['LAT'], df1['LON'])]

print('Loaded frame')

from worker import workaround
import multiprocessing
from multiprocessing import Pool, Manager
from multiprocessing import get_context
tester = df2.iloc[0:10]
print('Initiated tester')
if __name__ == '__main__':
    try:
        multiprocessing.set_start_method('spawn')
    except RuntimeError:
        pass
    print('In worker now')
    p = get_context("spawn").Pool()
    mgr = Manager()
    config = mgr.Namespace()
    config.df1 = df1
    config.df2 = tester
    result_list = list(tqdm.tqdm(p.map(workaround, list(zip(tester['closest'], tester['collision_date'])))))
    print(result_list)
    p.close()
    p.join()