#!/usr/bin/env python

import os
import csv
import yaml
import pprint
from plotly import graph_objects as go
from plotly.subplots import make_subplots

def get_num_rows(config):
    n_ab = len(config["antibodies"])
    num_measurements = len(config["samples"]) * (n_ab + 1) + len(config["controls"])
    return num_measurements * config["num_replicates"]


def read_datafile(path, config, num_rows):
    with open(path, "r") as f:
        for _ in range(config["header_row"] - 1):
            f.readline()
        reader = csv.DictReader(f)
        data = [row for row in reader]
    return data[:num_rows]

def calc_ct_values(series, config):
    result = {}
    for sample in config["samples"]:
        ct_values = {}
        for ab in config["antibodies"]:
            avg = (float(next(series)["Ct"]) + float(next(series)["Ct"])) / 2
            ct_values[ab] = avg
        avg = (float(next(series)["Ct"]) + float(next(series)["Ct"])) / 2
        ct_values["input"] = avg
        result[sample] = ct_values
    return result

def take_IgG_diff(ct_values, config):
    result = {}
    for sample in config["samples"]:
        IgG_ct = ct_values[sample]["IgG"]
        result[sample] = {test:(IgG_ct - ct) for (test,ct) in ct_values[sample].items()}
    return result

def calc_final_ratio(IgG_diffs, config):
    result = {}
    for sample in config["samples"]:
        i = IgG_diffs[sample]["input"]
        result[sample] = {test:pow(2, ab) / ((100 / config["input_volume"]) * (pow(2, i))) for (test,ab) in IgG_diffs[sample].items()}
    return result

def create_table(final, config):
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'
    antibodies = config["antibodies"][:6]
    samples = config["samples"]
    registers = list(final.keys())
    cols = round(pow(len(antibodies), 0.5))
    rows = len(antibodies) // cols
    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=antibodies,
        vertical_spacing=0,
        horizontal_spacing=0.01,
        specs=[[{"type": "table"} for _ in range(cols)] for _ in range(rows)],)


    for i, ab in enumerate(antibodies):
        values = []
        for reg in registers:
            values.append([final[reg][sample][ab] for sample in samples])
        fig.add_trace(go.Table(
            header=dict(
                values=[""] + list(final.keys()),
                line_color='darkslategray',
                fill_color='grey',
                align=['left','center'],
                font=dict(color='white', size=12)
            ),
            cells=dict(
                values=[samples] + values,
                line_color='darkslategray',
                # 2-D list of colors for alternating rows
                fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor]*5],
                align = ['left', 'center'],
                font = dict(color = 'darkslategray', size = 11),
                format = ["", ".3f"],
                )),
            row=i//cols + 1, col=(i % cols) + 1)
    fig.show()

def create_graph(final, config):
    antibodies = config["antibodies"][:6]
    cols = round(pow(len(antibodies), 0.5))
    rows = len(antibodies) // cols
    colors = ["red", "green", "blue", "purple", "pink", "yellow"]
    fig = make_subplots(
        rows=rows, cols=cols,
        shared_xaxes=True,
        subplot_titles=antibodies,
        vertical_spacing=0.03,
        horizontal_spacing=0.03,)

    for i, ab in enumerate(antibodies):
        x = list(final.keys())
        for color, sample in enumerate(config["samples"]):
            y = [final[reg][sample][ab] for reg in x]
            fig.add_trace(go.Scatter(x=x, y=y, name=sample,
                                     line=go.scatter.Line(color=colors[color])),
                          row=i//cols + 1, col=(i % cols) + 1)
    fig.show()

def main():
    raw_data = {}
    config = yaml.load(open("config.yaml"))
    num_rows = get_num_rows(config)
    directory = os.fsencode("data")
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        series_name = os.path.basename(filename)[:-4]
        raw_data[series_name] = read_datafile(f"data/{filename}", config, num_rows)
    ct_values = {reg:calc_ct_values(iter(data), config) for (reg,data) in raw_data.items()}
    IgG_diffs = {reg:take_IgG_diff(data, config) for (reg,data) in ct_values.items()}
    final = {reg:calc_final_ratio(data, config) for (reg,data) in IgG_diffs.items()}
    # print(final)
    create_table(final, config)
    create_graph(final, config)



if __name__ == "__main__":
    main()
