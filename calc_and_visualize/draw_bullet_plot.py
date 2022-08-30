import plotly.graph_objects as go

fig = go.Figure()

custom_snp_number = 154
overall_snp_number = 332

fig.add_trace(go.Indicator(
    value=17,
    gauge={
        'steps': [
            {'range': [0, 17], 'color': "blue"},
            {'range': [17, 100], 'color': "red"}],
        'bar': {'color': "blue"},
        'shape': "bullet",
        'bgcolor': "white",
        'axis': {'range': [None, 100], 'visible': False}},
    domain={'x': [0.1, 1], 'y': [0.5, 0.8]}))

fig.add_trace(go.Indicator(
    value=154,
    gauge={
        'steps': [
            {'range': [0, custom_snp_number], 'color': "limegreen"},
            {'range': [custom_snp_number, overall_snp_number], 'color': "red"}],
        'bar': {'color': "limegreen"},
        'shape': "bullet",
        'bgcolor': "white",
        'axis': {'range': [0, overall_snp_number]}},
    domain={'x': [0.1, 0.72], 'y': [0.3, 0.5]}))


fig.update_layout(
    grid={'rows': 2, 'columns': 1, 'pattern': "independent"},
    template={'data': {'indicator': [{
        'title': {'text': "<br><span style='color: black; font-size:0.8em'><b>Longevity<br>PRS<b></span>"},
        'mode': "number+delta+gauge",
        'delta': {'reference': 100}}, {
        'title': {'text': "<br><span style='color: black; font-size:0.8em'><b>SNPs<b></span>"},
        'mode': "number+delta+gauge",
        'delta': {'reference': overall_snp_number}}]
    }})

fig.update_layout(height=400)
fig.show()
fig.write_image("/home/alina_grf/progprojects/prs/calc_and_visualize/olga.svg")