#!/usr/bin/env python
import argparse
import os
import sys

import plotly.graph_objects as go
import pandas as pd
from helpers.transform import from_vol_to_at, from_yn_to_m3


def plot(inputdir: str, outputdir: str, exp: bool):
    df = pd.read_excel(r'{}'.format(inputdir), sheet_name='results',
                       index_col='index')
    x_range = list(df.columns)
    data = [(i[0], from_vol_to_at(i[1])) for i in df.iterrows()]

    xaxis = go.layout.XAxis(title='profundidade (um)', dtick=2, mirror=True, ticks='outside', showline=True,
                            rangemode='tozero',
                            linecolor='black', titlefont=dict(family='Arial', size=10))
    yaxis = go.layout.YAxis(title='concentração (at.%)', mirror=True, ticks='outside', showline=True,
                            rangemode='tozero',
                            linecolor='black', titlefont=dict(family='Arial', size=10), tickformat=".2%")
    layout = go.Layout(plot_bgcolor='white', xaxis=xaxis, yaxis=yaxis, width=700, height=350,

                       legend=dict(x=0.75, y=0.95), margin=dict(t=2, b=2, l=2, r=3))

    fig = go.Figure(layout=layout)
    dash_types = iter(["solid", "dot", "dash", "dashdot", "longdash", "longdashdot"])

    colors = iter(['#0048b3', '#e92123', '#13c23c', '#4c94ff', '#ffa500'])
    # Add traces
    for d in data:
        fig.add_trace(go.Scatter(x=x_range, y=d[1],
                                 mode='lines',
                                 name=str(d[0]),
                                 line=dict(dash=next(dash_types),
                                           color=next(colors)
                                           )
                                 )
                      )
    if exp:
        if exp == 'gas':
            df_exp = pd.read_excel(r'{}\ExpEPMAdata.xlsx'.format(os.path.dirname(os.path.abspath(__file__))),
                                   sheet_name='results',
                                   index_col='index')
            data = [(i[0], from_vol_to_at(from_yn_to_m3(i[1]))) for i in df_exp.iterrows()]
        elif exp == 'plasma':
                df_exp = pd.read_excel(r'{}\plasma_exp.xlsx'.format(os.path.dirname(os.path.abspath(__file__))),
                                       sheet_name='results',
                                       index_col='index')
                data = [(i[0], i[1]/100) for i in df_exp.iterrows()]

        x_range = list(df_exp.columns)
        fig.add_trace(go.Scatter(x=x_range, y=data[0][1],
                                 mode='markers',
                                 name=str(data[0][0]),
                                 marker=dict(symbol='circle',
                                             color='black'
                                             )
                                 )
                      )

    fig.write_image(r'{}'.format(outputdir))
    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Graph Maker')
    parser.add_argument('input', type=str,
                        help=r'Input Dir for data (e.g C:\Users\Julia\Documents\tcc\TF\codes\results\plot\20191029\plot_classicFickGasBeta.xlsx)')
    parser.add_argument('outdir', type=str,
                        help=r'Output Dir for graph (e.g C:\Users\Julia\Documents\tcc\TF\codes\helpers\test8.png)')
    parser.add_argument('--exp', type=str, help=r'Add experimental data (gas or plasma)')
    args = parser.parse_args()
    plot(args.input, args.outdir, args.exp)
