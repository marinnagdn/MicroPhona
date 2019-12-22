#!/usr/bin/env python3
#-*- coding : utf-8 -*-

from microphona_def import *


################# DATA ##################

BIRMEC=load_data("database/irreversibleReactions_withExtrMet_byBacteria")
IRMEC=load_data("database/REACT_irreversibleReactions_withExtrMet")
MKC=load_data("database/MET_idKegg_idCarveme")
BMEC=load_data("database/extrMet_idCarveme_byBacteria")
MCN=load_data("database/MET_idCarveme_wholeName")

############## Handle parameters ###############

if sys.argv!=1:
	print("Give parameters file please. See README.md for more informations")
	exit(1)

nameHandle = open(sys.argv[1], "r")
list_groups=[]
lines=nameHandle.readlines()
# remove spaces
lines = [line.replace(' ', '') for line in lines]


#default parameters
oriented=False
human=False
graph = 0
gml=False
list_index_bacteria = []
list_kegg_metabo = []
threshold = 0
suffix = ""


for line in lines :
	if line[0] == "#": #comment or blank line
		continue
	if re.search("groups",line):
		line = re.split("#",line)[0] #for comments
		list_groups=re.split(":",line[0:-1])[1]
		list_groups=re.split(",",list_groups)
	if re.search("human",line):
		line = re.split("#",line)[0] #for comments
		human=re.split(":",line[0:-1])[1]
		if human == "1":
			human = True
		else:
			human =False
	if re.search("oriented",line):
		line = re.split("#",line)[0] #for comments
		oriented=re.split(":",line[0:-1])[1]
		if oriented == "1":
			oriented = True
		else:
			oriented =False
	if re.search("gml",line):
		line = re.split("#",line)[0] #for comments
		gml=re.split(":",line[0:-1])[1]
		if gml == "1":
			gml = True
		else:
			gml=False
	if re.search("graph",line):
		line = re.split("#",line)[0] #for comments
		graph=re.split(":",line[0:-1])[1]
		if graph == "1":
			graph = True
		else:
			graph=False
	if re.search("suffix",line):
		line = re.split("#",line)[0] #for comments
		suffix=re.split(":",line[0:-1])[1]
	if re.search("threshold",line):
		line = re.split("#",line)[0] #for comments
		threshold=float(re.split(":",line)[1])
	if re.search("metabolites",line):
		line = re.split("#",line)[0] #for comments
		list_kegg_metabo=re.split(":",line)[1]
		list_kegg_metabo=re.split(",",list_kegg_metabo)
	if re.search("bacteria",line):
		line = re.split("#",line)[0] #for comments
		list_index_bacteria=re.split(":",line)[1]
		list_index_bacteria=re.split(",",list_index_bacteria)
		i = 0
		for i in range(len(list_index_bacteria)):
			list_index_bacteria[i]=int(list_index_bacteria[i])
			i+=1

nameHandle.close()
nameHandle = open("logo.txt","r")
lines = nameHandle.read().splitlines()

#print parameters
print("\nThanks for using MicroPhona !")
#logo
for line in lines:
	print(line)
nameHandle.close()

print("\nParameters:")
print("\n***Input gml files : ")
for i in range(len(list_groups)):
	print(list_groups[i],".gml",sep="")
print("\n***Output files : ")
if gml:
	for i in range(len(list_groups)):
		print(list_groups[i],"_",suffix,".gml",sep="")
if graph:
	for i in range(len(list_groups)):
		print(list_groups[i],"_",suffix,".png",sep="")
if oriented:
	print("\nOriented graph : YES (/!\\ only irreversible reactions data)")
else:
	print("Oriented graph : NO")
if human:
	print("Homo_sapiens node : YES")
else:
	print("Homo_sapiens node : NO")
print("Threshold :",threshold)

#bacteria
if list_index_bacteria != []:
	list_bacteria = []
	nameHandle = open("index_bacteria.txt","r")
	lines=nameHandle.read().splitlines()
	for line in lines:
		split = re.split(":",line)
		if int(split[0]) in list_index_bacteria:
			list_bacteria.append(split[1][1:])
	if list_bacteria != []:
		print("\n***Bacteria found :")
		for i in range(len(list_bacteria)):
			print(list_bacteria[i])

list_carveme_metabo = []
if list_kegg_metabo != []:
	for x in range(len(list_kegg_metabo)):
		if list_kegg_metabo[x] in MKC.keys():

			list_carveme_metabo.append(MKC[list_kegg_metabo[x]]) #get database id with kegg id

	if list_carveme_metabo != []:
		print("\n***Metabolites found :")
		for i in range(len(list_carveme_metabo)):
			print(MCN[list_carveme_metabo[i]])
	else :
		print("\n***No corresponding metabolite was found")


########### manipulation of data ############

#read gml files
list_graphs=[]
for i in range(len(list_groups)):
	list_graphs.append(nx.read_gml(list_groups[i]+".gml"))

#remove isolate nodes
for G in list_graphs : 
	 G.remove_nodes_from(list(nx.isolates(G)))

#rename the nodes
undesirable_character = '[]'
rename_dict={}
for G in list_graphs :
    for i in G.nodes():
        genre_espece_1 = re.split("g__",i)[1]
        genre_1 = re.split(";s__",genre_espece_1)[0]
        espece_1 = re.split(";s__",genre_espece_1)[1]
        bact1 = str(genre_1+'_'+espece_1)
        rename_dict[i] = bact1.translate({ord(i): None for i in undesirable_character})
for i in range(len(list_graphs)) :
    list_graphs[i]=nx.relabel_nodes(list_graphs[i],rename_dict)


########### create dict of edges ############
list_interactions = []
#preparation to get the weight
for i in range( len(list_groups) ) :
	list_interactions.append([])
	for  (a,b) in list_graphs[i].edges() :
		list_interactions[i].append( [ (round(d['weight'],2)) for (u,v,d) in list_graphs[i].edges(data=True) if (u == a and v == b) or (u == b and v ==a)][0])

dict_corr={}
nbr_graph = 0

for G in list_graphs :
	cpt = 0 

	for i in G.edges():
		bacterie1 = i[0]
		bacterie2 = i[1]
		temp_str = bacterie1+'-'+bacterie2
		temp_str_inverse = bacterie2+'-'+bacterie1
#ajout dans le dict
		if temp_str in dict_corr :
			dict_corr[temp_str][nbr_graph] = float(list_interactions[nbr_graph][cpt])
		elif temp_str_inverse in dict_corr :
			dict_corr[temp_str][nbr_graph] = float(list_interactions[nbr_graph][cpt])
		else :
			dict_corr[temp_str] = [0.0,0.0,0.0]
			dict_corr[temp_str][nbr_graph] = float(list_interactions[nbr_graph][cpt])

		cpt+=1
	nbr_graph+=1


########### design graphs ############
nbr_graph=0
present_bacteria = []

for G in list_graphs :
	
	#node color
	for n in G.nodes(data=True):
		bacterie = n[0]
		data = n[1]
		if bacterie in BIRMEC.keys():
			data["color"]="#345DFF" #default color of nodes
			present_bacteria.append(bacterie)
		else:
			data["color"]="#9C9C9C" #model not available in the database


	#edges color depending of correlation score
	for e in G.edges(data=True):
		bacterie1=e[0]
		bacterie2=e[1]
		data=e[2]
		corr_score = data["weight"]
		data["style"] = "solid" #by default
		if corr_score>0:
			data["color"]="#FF0101"
		elif corr_score<0:
			data["color"]="#21E8FF"
		# if corr_score>0 and sens==1:
		# 	 data["color"]="#FF0101"
		# elif corr_score>0 and sens==2:
		# 	data["color"]="#780303"
		# elif corr_score<0 and sens == 1:
		# 	data["color"]="#21E8FF"
		# elif corr_score<0 and sens == 2:
		# 	data["color"]="#325CBA"
		elif corr_score==0:
			data["color"]="#000000"
			data["weight"]=0 
		if bacterie1 not in present_bacteria or bacterie2 not in present_bacteria:
			data["color"]="#919191"
			data["weight"]=0 
			data["metabolites"]=""

		if oriented:
			new_edges = []
			data["sens"] = 1 #sens of the reaction

			if data["weight"]!=0: #remove 0 correlation
				if bacterie1 in BIRMEC.keys() and bacterie2 in BIRMEC.keys():
					#from bacteria one to bacteria two
					common = list(set(get_extracellular_products(bacterie1)).intersection(get_extracellular_reactants(bacterie2)))
					data["metabolites"]=common
					#from bacteria two to bacteria one 
					# ---> add new edge if not empty !
					common = list(set(get_extracellular_reactants(bacterie1)).intersection(get_extracellular_products(bacterie2)))
					if common != []:
						# create new edge to incorporate in list_graphs
						new_edges.append([bacterie2,bacterie1,data["weight"],common])
			else:
				data["metabolites"]=""

			# add new edges to list_graphs
			for ne in new_edges:
				G.add_edge(ne[0],ne[1],weight=ne[2],metabolites=ne[3],sens=2) #2 for the other sense between the same bacteria

		else : #not oriented
			if data["weight"] != 0:
				if bacterie1 in BMEC.keys() and bacterie2 in BMEC.keys():
					metabolites_in_common = list(set(BMEC[bacterie1]).intersection(BMEC[bacterie2]))
					data["metabolites"]=metabolites_in_common
				else:
					data["metabolites"]="" #bacteria not present in the dabatase
			else:
				data["metabolites"]="" #correlation 0

	nbr_graph+=1

####### add homo sapiens ########

if human:
	if oriented:
		for G in list_graphs :
			G.add_node("Homo_sapiens",color="#989898")
			for n in G.nodes(data=True):
				bacteria=n[0]
				if bacteria != "Homo_sapiens":
					if bacteria in BIRMEC.keys():
						#from human to bacteria
						common = list(set(get_extracellular_products("Homo_sapiens")).intersection(get_extracellular_reactants(bacteria)))
						if common != []:
							G.add_edge("Homo_sapiens",bacteria,metabolites=common,color="#989898",weight=0)
						#from bacteria to human
						common = list(set(get_extracellular_reactants("Homo_sapiens")).intersection(get_extracellular_products(bacteria)))
						if common != []:
							G.add_edge(bacteria,"Homo_sapiens",metabolites=common,color="#444444",weight=0,style="solid")
	else: # not oriented
		for G in list_graphs :
			G.add_node("Homo_sapiens",color="#989898")
			for i in G.nodes(data=True):
				if i[0] != "Homo_sapiens":
					if i[0] in BMEC.keys():
						metabolites_in_common = list(set(BMEC["Homo_sapiens"]).intersection(BMEC[i[0]]))
						G.add_edge("Homo_sapiens",i[0],metabolites=metabolites_in_common,color="#000000",weight=1,style="solid")

############# hightlight interest bacteria ###############

for G in list_graphs:
	for n in G.nodes(data=True):
		bacteria = n[0]
		data = n[1]
		if bacteria in list_bacteria:
			n[1]["color"]="#FFC400"

########### interest metabolites ############
print()
if list_carveme_metabo != [] :
	nbr_groupe=0
	for G in list_graphs:
		list_nodes = []
		for e in G.edges(data=True):
			bacteria1=e[0]
			bacteria2=e[1]
			data=e[2]
			#edge will be selected in any of the given metabolites appears
			if any(i in data["metabolites"] for i in list_carveme_metabo):
				list_nodes.append(bacteria1)
				list_nodes.append(bacteria2)
				# + we keep the neighbors connected nodes
				list_nodes.extend(list(G.neighbors(bacteria1)))
				list_nodes.extend(list(G.neighbors(bacteria2)))
				data["style"]="dashed" # dashed for metabolite annotated edges
		list_nodes = list(dict.fromkeys(list_nodes))
		subG=G.subgraph(list_nodes) #extract subgraph
		if len(subG.edges())>0: #if their is a subgraph
			subG=graph_threshold(subG,threshold)
			save_graph(subG,list_groups[nbr_groupe]+"_"+suffix+".png", oriented)
			if gml:
				nx.write_gml(subG,list_groups[nbr_groupe]+"_"+suffix+".gml")
				print("Generated",list_groups[nbr_groupe]+"_"+suffix+".gml")
		else:
			print("Any of the given metabolites were found on the graph ",list_groups[nbr_groupe],".gml",sep="")
		nbr_groupe+=1	

else: #no special metabolites
	nbr_groupe=0
	for G in list_graphs:
		subG=graph_threshold(G,threshold)
		save_graph(subG,list_groups[nbr_groupe]+"_"+suffix+".png",oriented)
		if gml:
				nx.write_gml(subG,list_groups[nbr_groupe]+"_"+suffix+".gml")
				print("Generated",list_groups[nbr_groupe]+"_"+suffix+".gml")
		nbr_groupe +=1
