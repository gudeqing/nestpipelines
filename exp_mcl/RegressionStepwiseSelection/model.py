# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 08:07:22 2019

@author: Sinan Talha Hascelik
modified by gudeqing
"""

# Import Libraries
import stepwiseSelection as ss
import pandas as pd


def var_selection(data,model_type='logistic', method='backward', y_col='y',
				  x_cols:list=None, drop_cols:list=None, significant=0.05,
				  varchar_process="dummy_dropfirst", factorize:tuple=None):
	# Read Data
	df = pd.read_csv(data, header=0, index_col=0, sep=None, engine='python')
	if factorize is not None:
		fac = {x:int(y) for x, y in zip(factorize[::2], factorize[1:][::2])}
		y_data = [fac[x] for x in df[y_col]]
	else:
		if df[y_col].dtype != int:
			y_data, rule = df[y_col].factorize()
			print(dict(zip(rule, range(len(rule)))))
		else:
			print('target col seems to be already factorized!')
			y_data = df[y_col]
	df[y_col] = y_data
	# Dependent and Independent Variables
	X = df.drop(columns= y_col)
	if drop_cols:
		for each in drop_cols:
			X = X.drop(columns=each)
	if x_cols:
		X = df[x_cols]
	y = df[y_col]
	corr = X.corr()
	corr.to_csv('corr.csv')
	tmp = X.copy()
	tmp[y_col] = y
	tmp.to_csv('input_variable.data.csv')

	# Magic Happens
	if method == 'backward':
		final_vars, iterations_logs = ss.backwardSelection(X,y, model_type=model_type, sl=significant, varchar_process=varchar_process)
	else:
		final_vars, iterations_logs = ss.forwardSelection(X,y, model_type=model_type, sl=significant, varchar_process=varchar_process)
	print(final_vars)

	# Write Logs To .txt
	iterations_file = open("Iterations.log","w")
	iterations_file.write(iterations_logs)
	iterations_file.close()

if __name__ == '__main__':
	from xcmds import xcmds
	xcmds.xcmds(locals(), include=['var_selection'])

