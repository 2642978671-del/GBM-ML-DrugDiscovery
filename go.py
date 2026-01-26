import pandas as pd
import numpy as np
import plotly.express as px

# Load the data
file_path = 'GO.txt'
try:
    data = pd.read_csv(file_path, sep="\t", encoding="utf-8")
except UnicodeDecodeError:
    data = pd.read_csv(file_path, sep="\t", encoding="gbk")

# Filter top 10 terms per ontology (BP, CC, MF) based on p-value
data['-log10(pvalue)'] = -data['pvalue'].apply(lambda x: -1 * np.log10(x))  # Calculate -log10(pvalue)
top_terms = data.groupby('ONTOLOGY').apply(lambda x: x.nsmallest(10, 'pvalue')).reset_index(drop=True)

# Custom color scale
custom_colors = ['#ccccff', '#b3d8b2', '#fea500', '#ffcccc']

# Plot using Plotly Express
fig = px.scatter(
    top_terms,
    x="GeneRatio",
    y="Description",
    color="-log10(pvalue)",  # Color by -log10(pvalue)
    size="Count",            # Bubble size by Count
    facet_col="ONTOLOGY",    # Create panels by ONTOLOGY
    color_continuous_scale=custom_colors,  # Use custom color scale
    title="GO Enrichment Analysis",
    labels={
        "GeneRatio": "Gene Ratio",
        "-log10(pvalue)": "-log10(pvalue)",
        "Description": "GO Term"
    },
    template="simple_white"
)

# Update layout for font styles and layout adjustments
fig.update_layout(
    title_x=0.5,
    title_font=dict(size=22, family='Times New Roman', color='black'),  # Title font
    legend_title=dict(font=dict(size=18, family='Times New Roman', color='black')),  # Legend font
    legend=dict(borderwidth=1),
    margin=dict(l=100, r=100, t=100, b=50),
    font=dict(family="Times New Roman", size=16, color="black")  # Global font settings
)

# Update axes for font styles
fig.update_xaxes(
    categoryorder="total ascending",  # Sort by total ascending for better readability
    title_font=dict(size=20, family='Times New Roman', color='black'),
    tickfont=dict(size=18, family='Times New Roman', color='black'),
)

fig.update_yaxes(
    title_font=dict(size=20, family='Times New Roman', color='black'),
    tickfont=dict(size=18, family='Times New Roman', color='black')
)

# Update color bar font styles
fig.update_coloraxes(
    colorbar_title_font=dict(size=18, family='Times New Roman', color='black'),
    colorbar_tickfont=dict(size=16, family='Times New Roman', color='black')
)

# Show the figure
fig.show()
