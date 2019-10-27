#!/usr/bin/env python
import argparse

import plotly.graph_objects as go
# Create random data with numpy
import numpy as np
import pandas as pd
from helpers.transform import from_vol_to_at


def plot(inputdir: str, outputdir: str):
    df = pd.read_excel(r'{}'.format(inputdir), sheet_name='results',
                       index_col='index')
    x_range = list(df.columns)
    data = [(i[0], from_vol_to_at(i[1])) for i in df.iterrows()]

    xaxis = go.layout.XAxis(title='profundidade (um)', dtick=2, mirror=True, ticks='outside', showline=True, rangemode='tozero',
                            linecolor='black', titlefont=dict(family='Arial', size=10))
    yaxis = go.layout.YAxis(title='concentração (at.%)', mirror=True, ticks='outside', showline=True, rangemode='tozero',
                            linecolor='black', titlefont=dict(family='Arial', size=10), tickformat=".2%")
    layout = go.Layout(plot_bgcolor='white', xaxis=xaxis, yaxis=yaxis, legend=dict(x=0.75, y=0.95))

    fig = go.Figure(layout=layout)
    # Add traces
    for d in data:
        fig.add_trace(go.Scatter(x=x_range, y=d[1],
                                 mode='lines',
                                 name=str(d[0])))

    fig.write_image(r'{}'.format(outputdir))
    print("done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Graph Maker')
    parser.add_argument('input', type=str,
                        help=r'Input Dir for data (e.g C:\Users\Julia\Documents\tcc\TF\codes\results\plot\plot_test2.xlsx)')
    parser.add_argument('outdir', type=str,
                        help=r'Output Dir for graph (e.g C:\Users\Julia\Documents\tcc\TF\codes\helpers\test8.png)')
    args = parser.parse_args()
    plot(args.input, args.outdir)
