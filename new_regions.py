with open("MonsterPlexRegionsFileSuperCont.txt",'r') as read:
	num_1 =0
	num_2 =0
	gate = True
	last = ""
	for line in read:
		if line.split("\t")[0] == last or gate:
			gate = False
			num_1 = 1 + num_2
			num_2 = int(line.split("\t")[2]) - int(line.split("\t")[1]) + num_2
			with open("MonsterPlexRegionsFileSuperCont_small_index.txt",'a') as write:
				write.write(line.split("\t")[0] + "\t" + str(num_1) + "\t" + str(num_2) + "\n")
			last = line.split("\t")[0]
		else:
			num_1 =0
			num_2 =0
			gate = True