with open("MonsterPlexRegionsFileFadix.txt",'r') as read_regions_names:
	with open("MonsterPlexSitesList.txt",'r') as read_sites:
		with open("MonsterPlexSitesList_small_index.txt",'a') as write_sites:
			site_names = []
			for line in read_regions_names:
				site_names.append(line.rstrip())
			print(site_names)
			for line in read_sites:
				for name in site_names:
					if line.split(' ')[0] == name.split(':')[0] and int(line.split(' ')[1].rstrip()) in range(int(name.split(':')[1].split('-')[0]),int(name.split(':')[1].split('-')[1])):
						new_pos = int(line.split(' ')[1].rstrip()) - int(name.split(':')[1].split('-')[0])
						write_sites.write(f"{name} {new_pos+1}\n")
