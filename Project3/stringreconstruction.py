from collections import defaultdict
from sys import setrecursionlimit

def DeBruijn3(patterns):
	graph = []
	dict = defaultdict(list)
	for each in patterns:
		dict[each[:-1]].append(each[1:])
	for each in sorted(dict):
		string = each + ' -> '
		for item in dict[each]:
			string = string + item + ','
		graph.append(string[:-1])
	return graph

def EulerianPath(graph):
	dict = {}
	indegree_dict = {}
	bad_indegree = 0
	bad_outdegree = 0
	bad_in_degree_bool = False
	
	for line in graph:
		node_index = line.find(' ')
		connect_index = line.find('>')
		node = line[:node_index]
		connect_str = line[connect_index+2:]
		connect_list = connect_str.split(',')
		dict[node] = connect_list

	non_values = []

	for key in dict.keys():
		for value in dict[key]:
			if value not in dict.keys():
				non_values.append(value)

	for item in non_values:
		dict[item] = []
		bad_outdegree = item

	for key in dict:
		indegree_dict[key] = 0

	for key in dict:
		for value in dict[key]:
			indegree_dict[value] += 1

	for key in dict:
		if (len(dict[key]) != indegree_dict[key]):
			if key != bad_outdegree:
				bad_indegree = key
				bad_in_degree_bool = True

	if (len(non_values)!=0):
		dict[bad_outdegree].append(bad_indegree)

	path = FindCircuit(bad_indegree,dict)

	if (len(non_values) != 0):
		path = path[1:]

	path = path[::-1]
	
	if (not bad_in_degree_bool):
		brk = path.index(bad_outdegree)
		part1 = path[brk:]
		part2 = path[:brk]
		path = part1 + part2
	
	return path
			
def FindCircuit(node,dict):
	dict_copy = dict.copy()
	final_path = []
	if (len(dict_copy[node]) == 0):
		final_path.append(node)
		return final_path
	else:
		while (len(dict_copy[node]) != 0):
			next = dict_copy[node].pop()
			final_path += FindCircuit(next,dict_copy)
		final_path.append(node)
	return final_path

def combine(ls):
	combined = ls[0]
	for i in range(1,len(ls)):
		combined = combined + ls[i][-1]
	return combined	

#patts = ['CTTA','ACCA','TACC','GGCT','GCTT','TTAC']
#print combine(EulerianPath(DeBruijn3(patts)))

setrecursionlimit(5000)

patterns = []
input = open('dataset_203_6.txt','r')
for line in input:
	if (line.strip()[0] == '>'):
		continue
	patterns.append(line.rstrip())
input.close()
print len(patterns)
string = combine(EulerianPath(DeBruijn3(patterns)))

output = open('stringreconstruction.txt','w')
output.write(string)
output.close()