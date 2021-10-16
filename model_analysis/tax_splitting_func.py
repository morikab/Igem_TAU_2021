#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import json

import random


# In[2]:




# load xl file to dataframe
data = pd.read_excel("rotem excel.xlsx", sheet_name='Sheet1')

# open organism dict
with open('data_for_analysis/org_name_to_dict.json') as f:
    org_dict = json.load(f)


# In[3]:


def splitting_org(org_dict, tax_con, x):
    s = data.sample(n = x )
    grouped = s.groupby(tax_con)         # group according to tax level
    l_grouped = list(grouped)
    
    # sort the list according to its length
    l_grouped.sort(key=lambda x: len(x[1]))

    # create organism list- each object is organisms with the same tax_name
    org_group_list = []
    for i in range(len(l_grouped)):
        org_group_list.append(l_grouped[i][1]['Organism'])

    # create tax_name list- matches org_list (not sure if needed)
    tax_list = []
    for i in range(len(l_grouped)):
        tax_list.append(l_grouped[i][0])
        
    # divide the organisms groups according to length
    org_group1 = org_group_list[0::2]     # take even indices
    org_group2 = org_group_list[1::2]     # take odd indices
    
    # extract the organism's names (in order to compare later to the json dict)
    org_list1 = []
    for i in range(len(org_group1)):
        for j in range(len(org_group1[i])):
            org_list1.append(org_group1[i].loc[org_group1[i].index[j]])
        
    org_list2 = []
    for i in range(len(org_group2)):
        for j in range(len(org_group2[i])):
            org_list2.append(org_group2[i].loc[org_group2[i].index[j]])
            
    # create two dictionaries - one for each subgroup, according to json dict
    org_group1_dict = {}
    for i in range(len(org_list1)):
        org_group1_dict[org_list1[i]] = org_dict[org_list1[i]]

    org_group2_dict = {}
    for i in range(len(org_list2)):
        org_group2_dict[org_list2[i]] = org_dict[org_list2[i]]
    
    return org_group1_dict,org_group2_dict



def choose_best_split(org_dict, tax_con, x, run_times=50):
    groups_tuple = []
    for i in range(run_times):
        groups_tuple.append(splitting_org(org_dict, tax_con, x))  # we get list of tuples, each run creates a tuple, tuple[i][0]=subgroup1(dict1), tuple[i][1]=subgroup2(dict2)

    running_dict = {}
    for i in range(run_times):
        try:
            optimized_dict = groups_tuple[i][0]
            deoptimized_dict = groups_tuple[i][1]
            groups_ratio = len(optimized_dict) / len(deoptimized_dict)
            running_dict[groups_ratio+ 1/groups_ratio] = \
                {'optimized_dict':optimized_dict, 'deoptimized_dict':deoptimized_dict}
        except:
            continue

    if len(running_dict) >0:
        best_run = min(running_dict.keys())
        optimized_dict = running_dict[best_run]['optimized_dict']
        deoptimized_dict = running_dict[best_run]['deoptimized_dict']
        return optimized_dict,deoptimized_dict
    else:
        return None

def format_output(seq, org_dict, tax_con, x, run_times=50):
    software_inp = {}
    software_inp['tuning_param'] = 0.5
    software_inp['sequence'] = seq

    software_inp['organisms'] = {}
    optimized_dict, deoptimized_dict = choose_best_split(org_dict, tax_con = tax_con, x=x , run_times=run_times)
    for org_name, org_dict in optimized_dict.items():
        org_dict['optimized'] = True
        org_dict['tai_profile'] = {}
        org_dict['tai_std'] = {}
        org_dict['tai_avg'] = {}
        software_inp['organisms'][org_name] = org_dict
    for org_name, org_dict in deoptimized_dict.items():
        org_dict['optimized'] = False
        org_dict['tai_profile'] = {}
        org_dict['tai_std'] = {}
        org_dict['tai_avg'] = {}
        software_inp['organisms'][org_name] = org_dict

    return software_inp



