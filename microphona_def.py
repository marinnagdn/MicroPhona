#!/usr/bin/env python3
#-*- coding : utf-8 -*-

import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import re 
import pickle
import math
import os



############## Functions ###############

"""load dictionary"""
def load_data(data): #Homo_sapiens for humans
    with open(data, 'rb') as f:
        mon_depickler = pickle.Unpickler(f)
        data = mon_depickler.load()
    return data

BIRMEC=load_data("database/irreversibleReactions_withExtrMet_byBacteria")
IRMEC=load_data("database/REACT_irreversibleReactions_withExtrMet")

#********************

def get_extracellular_products(bacteria):
	products = []
	irreversible_reactions=BIRMEC[bacteria]
	for ir in irreversible_reactions:
		for p in IRMEC[ir][1]:
			if re.match(r".*_e$",p):
				products.append(p)
	return products

#********************

def get_extracellular_reactants(bacteria):
	reactants = []
	irreversible_reactions=BIRMEC[bacteria]
	for ir in irreversible_reactions:
		for r in IRMEC[ir][0]:
			if re.match(r".*_e$",r):
				reactants.append(r)
	return reactants

#********************

def graph_threshold(G,threshold):
	list_interest_nodes = []
	for e in G.edges(data=True):
		bacteria1=e[0]
		bacteria2=e[1]
		data=e[2]
		if abs(data["weight"])>=threshold:
			list_interest_nodes.append(bacteria1)
			list_interest_nodes.append(bacteria2)
	subG = G.subgraph(list_interest_nodes)
	return subG

#********************

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]


"""graphical functions"""
def draw_labels_circular_layout(G, pos) :
	description = nx.draw_networkx_labels(G, pos, font_size=9, font_color='k')
	for node, t in description.items():
		position = (t.get_position()[0]*1.36, t.get_position()[1]*1.36)
		t.set_position(position)
		# t.set_rotation(math.atan2(y,x)/math.pi*180)
		angle_label = math.atan2(t.get_position()[1],t.get_position()[0])/math.pi*180
		if angle_label < -90 :
			t.set_rotation(angle_label+180)
		elif angle_label > 90 :
			t.set_rotation(angle_label-180)
		else :
			t.set_rotation(angle_label)

		t.set_clip_on(False)

#******************************

def format_margin_Graph(pos) :
	x_values, y_values = zip(*pos.values())
	x_max = max(x_values)
	x_min = min(x_values)
	x_margin = (x_max - x_min) * 0.4 #Extend the margin by 40%.

	y_max = max(y_values)
	y_min = min(y_values)
	y_margin = (y_max - y_min) * 0.4

	#Check for the longest distance from the center.
	if x_max+x_margin > y_max+y_margin :
	#we want a square, the longest distance becomes the distance for the X and Y axis.
		plt.xlim(x_min - x_margin, x_max + x_margin)
		plt.ylim(x_min - x_margin, x_max + x_margin)
	else :
		plt.xlim(y_min - y_margin, y_max + y_margin)
		plt.ylim(y_min - y_margin, y_max + y_margin)

	plt.axis('off') #helps remove the margins of the graph

#******************************

def save_graph(G,output,oriented):
	#remove grey nodes (nodes not representated in the database)
	nodelist=[u for u,v in G.nodes(data=True) if v["color"]!="#9C9C9C"]
	G=G.subgraph(nodelist)

	edge_color = [d['color'] for u,v,d in G.edges(data=True)]
	node_color = [d['color'] for u,d in G.node(data=True)]
	style = [d['style'] for u,v,d in G.edges(data=True)]

	###### weights as number of metabolites per edge
	nbr_metabolites_per_interaction = [len(d["metabolites"]) for u,v,d in G.edges(data=True)] #get nbr of metabolites per interaction
	nbr_metabolites_per_interaction = remove_values_from_list(nbr_metabolites_per_interaction,0) #remove all 0 that corresponds to empty list
	mini = min(nbr_metabolites_per_interaction)
	maxi =   max(nbr_metabolites_per_interaction)

	#to avoid 0 problems !!!
	if mini == maxi:
		maxi=0.000001

	#normalisation of weights
	weights = [((len(d["metabolites"])-mini)/(maxi-mini))*5+1 for u,v,d in G.edges(data=True)]
	weights = [1 if x<0 else x for x in weights] #to get 1 if no metabolites

	#position coordinates
	pos = nx.circular_layout(G, dim=2, scale= 5/np.sqrt(G.order()), center=None)
	pos["Homo_sapiens"]=np.array([0,0]) #if homo sapiens node
	fig = plt.figure(figsize=(10,10),dpi=300)

	if oriented:
		nx.draw(G, pos=pos, edge_color=edge_color, width=weights,node_color=node_color, arrowstyle=' -|>',arrowsize = 20,node_size=100,alpha=0.7)
	else:
		nx.draw(G, pos=pos, edge_color=edge_color, width=weights,node_color=node_color,node_size=100,alpha=0.7,style=style)
	
	draw_labels_circular_layout(G, pos)
	format_margin_Graph(pos)

	# UNCOMMENT BELOWTO GET CORRELATION SCORES ON EDGES
	# edge labels = correlation score
	# edge_labels = nx.get_edge_attributes(G, 'weight')
	# for k,v in edge_labels.items():
	# 	edge_labels[k]=round(v,3)
	# nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

	plt.draw()
	plt.savefig(output, bbox_inches='tight') #bbox_inches controls the margins
	print("Generated",output)
	plt.clf() #clean plot

