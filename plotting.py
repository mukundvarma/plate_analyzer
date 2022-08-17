import plotly.express as px
import plotly.graph_objects as go
from analysis import fit_straight_line


def plot_rfu_panel(merged, ret=True):
    fig = px.scatter(
    merged, 
    x='Cycle', 
    y='RFU (smooth)', 
    color='Well',
    facet_col='group', 
    facet_col_wrap=4,
    hover_name='Well',
    hover_data={'RFU': ":.2f", 'RFU (smooth)': ":.2f", 'Linear Fit': ":.2f", "group": False},
    width=1000,
    height=1000,
    template='simple_white',
    color_discrete_sequence=px.colors.qualitative.T10[:6],


    )

    fig.update_traces(marker_size=2, showlegend=False)

    for i, group in enumerate(sorted(merged['group'].unique())):
        row = 4 -  (i // 4)
        col = 1 + (i % 4)
        each_df = merged.query('group == @group')
        fig.add_trace(
            go.Scatter(
                x=each_df['Cycle'],
                y=each_df['RFU'],
                mode="markers",
                ids=each_df.index,
                hoverinfo='skip',

                marker=dict(color='black', size=4, opacity=0.1)
            ),
            row=row, col=col
        )
        fig.add_trace(
            go.Scatter(
                x=each_df['Cycle'],
                y=each_df['Linear Fit'],
                mode="lines+markers",
                hoverinfo='skip',
                marker=dict(color='black', size=5, opacity=0.7),

            ),
            row=row, col=col
        )


    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(showlegend=False,)

    if ret:
        return fig


def plot_standard_curve(means):
    y_pred, x, m, c, r = fit_straight_line(means['Conc'], means['RFU (smooth)'])
    fig = px.scatter(means, 'Conc', 'RFU (smooth)',  title=f"Standard Curve")
    fig.add_trace(go.Scatter(x=x, y=y_pred, name='Linear Fit'))
    return m, r, fig


def plot_platemap(data_series, plate_map):
    fig = px.imshow(plate_map.applymap(lambda x: data_series.to_dict()[x] if x in data_series.index else 0, na_action=None),
                    aspect=1.2,
                    origin='top',
                    width=600,
                    height=400,
                    title=data_series.name
                   )

    fig.update_traces(text=plate_map, texttemplate="%{text}", hoverinfo='z', hovertemplate="Rate: %{z}")
    fig.update_yaxes(visible=False, showticklabels=False)
    fig.update_xaxes(visible=False, showticklabels=False)
    return fig