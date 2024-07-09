from __future__ import annotations
from typing import Union
import numpy as np
import pandas as pd


def resample(data: Union[np.array, pd.DataFrame], nsamples: int):
    if isinstance(data, pd.DataFrame):
        data = data.to_numpy().flatten()
    np.random.seed(nsamples)
    NDRAWS = len(data) * 100
    resamples = [data.min(), data.max()]
    unq, ids = np.unique(data, return_inverse=True)
    all_ids = np.random.choice(ids, size=NDRAWS, replace=True)
    ar = np.bincount(all_ids) / NDRAWS
    resamples.extend(np.random.choice(a=unq, size=nsamples, p=ar))
    return resamples
