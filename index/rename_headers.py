with open("70-15_small.fasta",'r') as read:
	with open("70-15_small2.fasta",'a') as write:
		used_sites =[]
		for line in read:
			if line[0] == ">":
				new_line = line.split(":")[0]
				if new_line not in used_sites:
					write.write("\n"+new_line+"\n")
					used_sites.append(new_line)
			else:
				write.write(line.rstrip())
			