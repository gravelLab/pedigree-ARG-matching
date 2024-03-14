import plotly.graph_objects as go

# Sample data (replace with your own data)
vertices = [1, 2, 3, 4, 5]
edges = [(1, 2), (1, 3), (2, 4), (2, 5)]

# Create graph object
graph = go.Figure()

# Add nodes
graph.add_trace(go.Scatter(x=[1, 2, 3, 4, 5], y=[1, 2, 3, 4, 5],
                           mode='markers+text',
                           marker=dict(size=20),
                           text=vertices,
                           textposition='bottom center'))

# Add edges
for edge in edges:
    x0, y0 = vertices.index(edge[0]) + 1, vertices.index(edge[1]) + 1
    x1, y1 = vertices.index(edge[1]) + 1, vertices.index(edge[0]) + 1
    graph.add_trace(go.Scatter(x=[x0, x1], y=[y0, y1],
                               mode='lines'))

# Customize layout
graph.update_layout(
    showlegend=False,
    hovermode='closest',
    title='Interactive Graph Visualization',
    xaxis=dict(showgrid=False, zeroline=False),
    yaxis=dict(showgrid=False, zeroline=False)
)

# Show the interactive graph
graph.show()
