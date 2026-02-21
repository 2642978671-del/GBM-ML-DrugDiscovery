# import pandas as pd
# import plotly.graph_objects as go
# import random
#
# # 尝试读取文件，指定编码格式
# file_path = 'KEGG.txt'  # 文件路径
# try:
#     data = pd.read_csv(file_path, sep="\t", encoding="utf-8")  # 首先尝试utf-8
# except UnicodeDecodeError:
#     data = pd.read_csv(file_path, sep="\t", encoding="gbk")  # 如果失败，尝试gbk
#
# # 准备桑基图数据
# # 提取唯一的Description和geneID
# descriptions = data['Description'].unique()
# genes = set()
# for gene_list in data['geneID']:
#     genes.update(gene_list.split('/'))
# genes = list(genes)
#
# # 映射Description和gene到索引
# description_indices = {desc: i for i, desc in enumerate(descriptions)}
# gene_indices = {gene: i + len(descriptions) for i, gene in enumerate(genes)}
#
# # 定义颜色列表，初始颜色
# base_colors = ['#b3d8b2', '#ccccff', '#fea500', '#ffcccc']
# # 如果颜色不足，生成更多随机颜色
#
# while len(base_colors) < (len(descriptions) + len(genes)):
#     random_color = "#{:06x}".format(random.randint(0, 0xFFFFFF))  # 生成随机RGB颜色
#     base_colors.append(random_color)
#
# # 为Description和Gene分配颜色
# description_colors = {desc: base_colors[i] for i, desc in enumerate(descriptions)}
# gene_colors = {gene: base_colors[i + len(descriptions)] for i, gene in enumerate(genes)}
#
# # 创建桑基图的source和target索引
# sources = []
# targets = []
# values = []
# link_colors = []
#
# for _, row in data.iterrows():
#     desc_index = description_indices[row['Description']]
#     for gene in row['geneID'].split('/'):
#         gene_index = gene_indices[gene]
#         sources.append(desc_index)
#         targets.append(gene_index)
#         values.append(row['p.adjust'])  # 使用p.adjust作为权重
#         link_colors.append(description_colors[row['Description']])  # 设置线条颜色为来源节点颜色
#
# # 创建桑基图标签
# labels = list(descriptions) + genes
# node_colors = list(description_colors.values()) + list(gene_colors.values())  # 合并节点颜色
#
# # 绘制桑基图
# fig = go.Figure(data=[go.Sankey(
#     node=dict(
#         pad=50,  # 增加节点之间的间距
#         thickness=40,  # 节点厚度
#         line=dict(color="black", width=1),
#         label=labels,
#         color=node_colors,  # 节点颜色
#         hoverlabel=dict(font=dict(size=14))  # 调整悬停字体
#     ),
#     link=dict(
#         source=sources,
#         target=targets,
#         value=values,
#         color=link_colors  # 线条颜色
#     )
# )])
#
# # 设置布局
# fig.update_layout(
#     title_text="KEGG Pathway Enrichment Analysis",
#     font=dict(
#         family="Times New Roman",
#         size=18,
#         color='black'
#     ),
#     title_x=0.5,  # 标题居中
#     margin=dict(l=200, r=200, t=100, b=100),  # 增加上下左右边距
# )
#
# # 显示图表
# fig.show()

# import pandas as pd
# import plotly.graph_objects as go
#
# # 尝试读取文件，指定编码格式
# file_path = 'KEGG.txt'  # 文件路径
# try:
#     data = pd.read_csv(file_path, sep="\t", encoding="utf-8")  # 首先尝试utf-8
# except UnicodeDecodeError:
#     data = pd.read_csv(file_path, sep="\t", encoding="gbk")  # 如果失败，尝试gbk
#
# # 准备桑基图数据
# # 提取唯一的Description和geneID
# descriptions = data['Description'].unique()
# genes = set()
# for gene_list in data['geneID']:
#     genes.update(gene_list.split('/'))
# genes = sorted(list(genes))  # 对基因进行排序，确保右侧按顺序排列
#
# # 映射Description和gene到索引
# description_indices = {desc: i for i, desc in enumerate(descriptions)}
# gene_indices = {gene: i + len(descriptions) for i, gene in enumerate(genes)}
#
# # 使用循环颜色，基于目标节点分配颜色
# base_colors = ['#b3d8b2', '#ccccff', '#fea500', '#ffcccc']
#
# # 为每个目标节点（gene）分配颜色
# gene_colors = {gene: base_colors[i % len(base_colors)] for i, gene in enumerate(genes)}
#
# # 为每个Description分配颜色（循环使用颜色）
# description_colors = {desc: gene_colors[genes[0]] for desc in descriptions}  # 起始分配基因颜色一致
#
# # 创建桑基图的source和target索引
# sources = []
# targets = []
# values = []
# link_colors = []
#
# for _, row in data.iterrows():
#     desc_index = description_indices[row['Description']]
#     for gene in row['geneID'].split('/'):
#         gene_index = gene_indices[gene]
#         sources.append(desc_index)
#         targets.append(gene_index)
#         values.append(row['p.adjust'])  # 使用p.adjust作为权重
#         link_colors.append(gene_colors[gene])  # 线条颜色为目标节点的颜色
#
# # 创建桑基图标签
# labels = list(descriptions) + genes
# node_colors = list(description_colors.values()) + list(gene_colors.values())  # 合并节点颜色
#
# # 绘制桑基图
# fig = go.Figure(data=[go.Sankey(
#     node=dict(
#         pad=70,  # 调整节点之间的上下间距
#         thickness=40,  # 节点厚度
#         line=dict(color="black", width=1),
#         label=labels,
#         color=node_colors,  # 节点颜色
#         hoverlabel=dict(font=dict(size=14))  # 调整悬停字体
#     ),
#     link=dict(
#         source=sources,
#         target=targets,
#         value=values,
#         color=link_colors  # 线条颜色为目标节点颜色
#     )
# )])
#
# # 设置布局
# fig.update_layout(
#     title_text="KEGG Pathway Enrichment Analysis",
#     font=dict(
#         family="Times New Roman",
#         size=18,
#         color='black',
#         weight='bold'
#     ),
#     title_x=0.5,  # 标题居中
#     margin=dict(l=300, r=200, t=50, b=30),  # 增加左边距 l=300，调整为更大的间距
# )
#
# # 显示图表
# fig.show()




# import pandas as pd
# import plotly.graph_objects as go
#
# # 尝试读取文件，指定编码格式
# file_path = 'KEGG.txt'  # 文件路径
# try:
#     data = pd.read_csv(file_path, sep="\t", encoding="utf-8")  # 首先尝试utf-8
# except UnicodeDecodeError:
#     data = pd.read_csv(file_path, sep="\t", encoding="gbk")  # 如果失败，尝试gbk
#
# # 准备桑基图数据
# # 提取唯一的Description和geneID
# descriptions = data['Description'].unique()
# genes = set()
# for gene_list in data['geneID']:
#     genes.update(gene_list.split('/'))
# genes = sorted(list(genes))  # 对基因进行排序，确保右侧按顺序排列
#
# # 映射Description和gene到索引
# description_indices = {desc: i for i, desc in enumerate(descriptions)}
# gene_indices = {gene: i + len(descriptions) for i, gene in enumerate(genes)}
#
# # 使用循环颜色
# base_colors = ['#b3d8b2', '#ccccff', '#fea500', '#ffcccc']
#
# # 为Description和Gene分配颜色（循环使用颜色），并增加透明度
# description_colors = {desc: f"rgba({int(base_colors[i % len(base_colors)][1:3], 16)},"
#                             f"{int(base_colors[i % len(base_colors)][3:5], 16)},"
#                             f"{int(base_colors[i % len(base_colors)][5:7], 16)}, 0.6)"
#                       for i, desc in enumerate(descriptions)}
# gene_colors = {gene: f"rgba({int(base_colors[i % len(base_colors)][1:3], 16)},"
#                       f"{int(base_colors[i % len(base_colors)][3:5], 16)},"
#                       f"{int(base_colors[i % len(base_colors)][5:7], 16)}, 0.6)"
#                for i, gene in enumerate(genes)}
#
# # 创建桑基图的source和target索引
# sources = []
# targets = []
# values = []
# link_colors = []
#
# for _, row in data.iterrows():
#     desc_index = description_indices[row['Description']]
#     for gene in row['geneID'].split('/'):
#         gene_index = gene_indices[gene]
#         sources.append(desc_index)
#         targets.append(gene_index)
#         values.append(row['p.adjust'])  # 使用p.adjust作为权重
#         link_colors.append(description_colors[row['Description']])  # 设置线条颜色为来源节点颜色
#
# # 创建桑基图标签
# labels = list(descriptions) + genes
# node_colors = list(description_colors.values()) + list(gene_colors.values())  # 合并节点颜色
#
# # 绘制桑基图
# fig = go.Figure(data=[go.Sankey(
#     node=dict(
#         pad=70,  # 调整节点之间的上下间距
#         thickness=40,  # 节点厚度
#         line=dict(color="black", width=1),
#         label=labels,
#         color=node_colors,  # 节点颜色
#         hoverlabel=dict(font=dict(size=14))  # 调整悬停字体
#     ),
#     link=dict(
#         source=sources,
#         target=targets,
#         value=values,
#         color=[f"rgba({int(c[1:3], 16)}, {int(c[3:5], 16)}, {int(c[5:7], 16)}, 0.6)"
#                if c.startswith('#') else c for c in link_colors]  # 线条颜色增加透明度
#     )
# )])
#
# # 设置布局
# fig.update_layout(
#     title_text="KEGG Pathway Enrichment Analysis",
#     font=dict(
#         family="Times New Roman",
#         size=18,
#         color='black',
#         weight='bold'
#     ),
#     title_x=0.5,  # 标题居中
#     margin=dict(l=300, r=200, t=50, b=30),  # 增加左边距 l=300，调整为更大的间距
# )
#
# # 显示图表
# fig.show()

#
# import pandas as pd
# import plotly.graph_objects as go
#
# # 尝试读取文件，指定编码格式
# file_path = 'KEGG.txt'  # 文件路径
# try:
#     data = pd.read_csv(file_path, sep="\t", encoding="utf-8")  # 首先尝试utf-8
# except UnicodeDecodeError:
#     data = pd.read_csv(file_path, sep="\t", encoding="gbk")  # 如果失败，尝试gbk
#
# # 准备桑基图数据
# # 提取唯一的Description和geneID
# descriptions = data['Description'].unique()
# genes = set()
# for gene_list in data['geneID']:
#     genes.update(gene_list.split('/'))
# genes = sorted(list(genes))  # 对基因进行排序，确保右侧按顺序排列
#
# # 映射Description和gene到索引
# description_indices = {desc: i for i, desc in enumerate(descriptions)}
# gene_indices = {gene: i + len(descriptions) for i, gene in enumerate(genes)}
#
# # 使用循环颜色，并增加透明度
# base_colors = ['#b3d8b2', '#ccccff', '#fea500', '#ffcccc','#f26cab','#b3de68','#fcb89b']
#
# # 为右侧目标节点（gene）分配颜色
# gene_colors = {gene: f"rgba({int(base_colors[i % len(base_colors)][1:3], 16)},"
#                       f"{int(base_colors[i % len(base_colors)][3:5], 16)},"
#                       f"{int(base_colors[i % len(base_colors)][5:7], 16)}, 0.6)"
#                for i, gene in enumerate(genes)}
#
# # 初始化左侧节点颜色字典
# description_colors = {desc: None for desc in descriptions}
#
# # 创建桑基图的source和target索引
# sources = []
# targets = []
# values = []
# link_colors = []
#
# for _, row in data.iterrows():
#     desc_index = description_indices[row['Description']]
#     for gene in row['geneID'].split('/'):
#         gene_index = gene_indices[gene]
#         sources.append(desc_index)
#         targets.append(gene_index)
#         values.append(row['p.adjust'])  # 使用p.adjust作为权重
#         link_colors.append(gene_colors[gene])  # 线条颜色为目标节点颜色
#         # 动态分配左侧节点颜色
#         if description_colors[row['Description']] is None:
#             description_colors[row['Description']] = gene_colors[gene]
#
# # 如果某些左侧节点没有连线颜色，则设置默认颜色
# for desc in description_colors:
#     if description_colors[desc] is None:
#         description_colors[desc] = "rgba(0,0,0,0.6)"  # 设置默认颜色为黑色半透明
#
# # 创建桑基图标签和节点颜色
# labels = list(descriptions) + genes
# node_colors = [description_colors[desc] for desc in descriptions] + \
#               [gene_colors[gene] for gene in genes]  # 左右两侧节点颜色
#
# # 绘制桑基图
# fig = go.Figure(data=[go.Sankey(
#     node=dict(
#         pad=70,  # 调整节点之间的上下间距
#         thickness=40,  # 节点厚度
#         line=dict(color="black", width=1),
#         label=labels,
#         color=node_colors,  # 节点颜色
#         hoverlabel=dict(font=dict(size=14))  # 调整悬停字体
#     ),
#     link=dict(
#         source=sources,
#         target=targets,
#         value=values,
#         color=link_colors  # 线条颜色为目标节点颜色
#     )
# )])
#
# # 设置布局
# fig.update_layout(
#     title_text="KEGG Pathway Enrichment Analysis",
#     font=dict(
#         family="Times New Roman",
#         size=18,
#         color='black',
#         weight='bold'
#     ),
#     title_x=0.5,  # 标题居中
#     margin=dict(l=200, r=200, t=50, b=5),  # 增加左边距 l=200，调整为更大的间距
# )
#
# # 显示图表
# fig.show()



import pandas as pd
import plotly.graph_objects as go

# 尝试读取文件，指定编码格式
file_path = 'KEGG.txt'  # 文件路径
try:
    data = pd.read_csv(file_path, sep="\t", encoding="utf-8")  # 首先尝试utf-8
except UnicodeDecodeError:
    data = pd.read_csv(file_path, sep="\t", encoding="gbk")  # 如果失败，尝试gbk

# 根据p.adjust值排序，并取前20条记录
data = data.sort_values(by='p.adjust').head(20)

# 提取唯一的Description和geneID
descriptions = data['Description'].unique()
genes = set()
for gene_list in data['geneID']:
    genes.update(gene_list.split('/'))
genes = sorted(list(genes))  # 对基因进行排序，确保右侧按顺序排列

# 映射Description和gene到索引
description_indices = {desc: i for i, desc in enumerate(descriptions)}
gene_indices = {gene: i + len(descriptions) for i, gene in enumerate(genes)}

# 使用循环颜色，并增加透明度
base_colors = ['#b3d8b2', '#ccccff', '#fea500', '#ffcccc', '#f26cab', '#b3de68', '#fcb89b']
base_colors = ['#b3d8b2', '#ccccff', '#fea500', '#ffcccc', '#f26cab', '#b3de68', '#fcb89b']

# 为右侧目标节点（gene）分配颜色
gene_colors = {gene: f"rgba({int(base_colors[i % len(base_colors)][1:3], 16)},"
                      f"{int(base_colors[i % len(base_colors)][3:5], 16)},"
                      f"{int(base_colors[i % len(base_colors)][5:7], 16)}, 1)"
               for i, gene in enumerate(genes)}

line_colors = {gene: f"rgba({int(base_colors[i % len(base_colors)][1:3], 16)},"
                      f"{int(base_colors[i % len(base_colors)][3:5], 16)},"
                      f"{int(base_colors[i % len(base_colors)][5:7], 16)}, 0.3)"
               for i, gene in enumerate(genes)}

# 初始化左侧节点颜色字典
description_colors = {desc: None for desc in descriptions}

# 创建桑基图的source和target索引
sources = []
targets = []
values = []
link_colors = []

for _, row in data.iterrows():
    desc_index = description_indices[row['Description']]
    for gene in row['geneID'].split('/'):
        gene_index = gene_indices[gene]
        sources.append(desc_index)
        targets.append(gene_index)
        values.append(row['p.adjust'])  # 使用p.adjust作为权重
        # 动态分配左侧节点颜色
        link_colors.append(line_colors[gene])
        if description_colors[row['Description']] is None:
            description_colors[row['Description']] = gene_colors[gene]

# 如果某些左侧节点没有连线颜色，则设置默认颜色
for desc in description_colors:
    if description_colors[desc] is None:
        description_colors[desc] = "rgba(0,0,0,0.6)"  # 设置默认颜色为黑色半透明

# 创建桑基图标签和节点颜色
labels = list(descriptions) + genes
node_colors = [description_colors[desc] for desc in descriptions] + \
              [gene_colors[gene] for gene in genes]  # 左右两侧节点颜色

# 绘制桑基图
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=70,  # 调整节点之间的上下间距
        thickness=40,  # 节点厚度
        line=dict(color="black", width=1),
        label=labels,
        color=node_colors,  # 节点颜色
        hoverlabel=dict(font=dict(size=14))  # 调整悬停字体
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors  # 线条颜色为目标节点颜色
    )
)])

# 设置布局
fig.update_layout(
    title_text="KEGG Pathway Enrichment Analysis (Top 20 by p.adjust)",
    font=dict(
        family="Times New Roman",
        size=18,
        color='black',
        weight='bold'
    ),
    title_x=0.5,  # 标题居中
    margin=dict(l=250, r=250, t=50, b=50),  # 增加左边距 l=250，调整为更大的间距
)

# 显示图表
fig.show()

