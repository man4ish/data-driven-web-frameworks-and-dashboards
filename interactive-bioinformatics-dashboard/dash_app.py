import dash
from dash import dcc, html, Input, Output
import plotly.express as px
import pandas as pd

# Sample dataset
df = px.data.gapminder()

# Initialize Dash app
app = dash.Dash(__name__)

# Layout with a dropdown and a graph
app.layout = html.Div([
    html.H1("Dynamic Country GDP Plot"),
    
    dcc.Dropdown(
        id='country-dropdown',
        options=[{'label': c, 'value': c} for c in df['country'].unique()],
        value='India',
        clearable=False
    ),
    
    dcc.Graph(id='gdp-graph')
])

# Callback to update graph based on selected country
@app.callback(
    Output('gdp-graph', 'figure'),
    Input('country-dropdown', 'value')
)
def update_graph(selected_country):
    filtered_df = df[df['country'] == selected_country]
    fig = px.line(filtered_df, x='year', y='gdpPercap',
                  title=f"GDP per Capita Over Time for {selected_country}")
    return fig

# Run app
if __name__ == '__main__':
    app.run(debug=True)
