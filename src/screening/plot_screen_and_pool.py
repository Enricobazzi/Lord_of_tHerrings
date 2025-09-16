import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt

def get_color_dict(dataset: str):
    if dataset == 'baltic_v_atlantic_snps' or dataset == 'all_snps':
        col_dict = {
            'Baltic': 'orangered',
            'Atlantic': 'royalblue',
            'EnglishChannel': 'skyblue',
            'NorthSea': 'darkblue',
            'IrishSea': 'darkseagreen',
            'SeaOfJapan': 'gold',
            'WhiteSea': 'darkviolet',
            'PechoraSea': 'darkviolet',
            'KandalakshaBay': 'darkviolet',
            'Balsfjord': 'darkviolet',
            'BarentsSea': 'darkviolet',
            'Pacific': 'gold',
            'Germany': 'darkred',
            }
    elif dataset == 'spring_v_autumn_snps':
        col_dict = {
            'Autumn': 'royalblue', 'Spring': 'orangered',
            'Winter': 'skyblue',
            'Summer': 'gold',
            'Mixed': 'darkviolet',
            'AutumnWinter': 'darkseagreen',
            'NA': 'grey',
            }
    return col_dict

