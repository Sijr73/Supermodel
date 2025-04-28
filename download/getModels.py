#!/usr/bin/python3


import urllib.request

import json

def getData(url):
	response = urllib.request.urlopen(url)
	data = response.read()
	text = data.decode('utf-8')
	return json.loads(text)



header = ["bigg_id", "organism", "gene_count", "metabolite_count", "reaction_count"]
print("\t".join(header))

if __name__ == "__main__":
	models = getData('http://bigg.ucsd.edu/api/v2/models')
	
	for m in models["results"]:
		values = []
		for i in header:
			values.append(str(m[i]))
		print("\t".join(values))
		
		response = urllib.request.urlopen("http://bigg.ucsd.edu/static/models/" + m["bigg_id"] + ".xml")
		data = response.read()
		f = open("../sourceData/models/"+ m["bigg_id"] + ".xml", "w")
		f.write(data.decode('utf-8'))
		f.close()
	
	
	
	
	
	
	
	
