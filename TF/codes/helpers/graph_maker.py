#!/usr/bin/env python

import plotly.graph_objects as go
# Create random data with numpy
import numpy as np
import pandas as pd

df = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\plot\plot_test.xlsx', sheet_name='results',
                   index_col='index')
x_range = list(df.columns)
data = [(i[0], i[1]) for i in df.iterrows()]

xaxis = go.layout.XAxis(title='profundidade (um)', dtick=1, mirror=True, ticks='outside', showline=True, rangemode='tozero',
                        linecolor='black', titlefont=dict(family='Arial', size=10))
yaxis = go.layout.YAxis(title='concentração (at.%)', mirror=True, ticks='outside', showline=True, rangemode='tozero',
                        linecolor='black', titlefont=dict(family='Arial', size=10))
layout = go.Layout(plot_bgcolor='white', xaxis=xaxis, yaxis=yaxis, legend=dict(x=0.85, y=0.95))

fig = go.Figure(layout=layout)
# Add traces
for d in data:
    fig.add_trace(go.Scatter(x=x_range, y=d[1],
                             mode='lines',
                             name=str(d[0])))

fig.write_image(r'C:\Users\Julia\Documents\tcc\TF\codes\helpers\test7.png')
print("done")
