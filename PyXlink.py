from itertools import combinations

class xlink:

	def __init__(self, data, proteins, name):
		self.name = name
		self.data = []
		all_proteins = []
		#read data file into list of tuples
		with open(data, "rU") as data_handle:
			for line in data_handle:
				if line.startswith("PROTEIN") != True:
					self.data.append((line.split()[0], line.split()[1], line.split()[2], line.split()[3]))
					#make a list of all proteins (to be uniqified later)
					all_proteins.append(line.split()[0])
					all_proteins.append(line.split()[2])
		#make a dictionary of unique proteins (sequences added later)
		self.proteins = {}
		for item in all_proteins:
			if item in self.proteins:
				continue
			self.proteins[item] = ""
		#add sequences to self.proteins
		with open(proteins, "rU") as fasta_handle:
			fasta_seqids = []
			fasta = "".join(fasta_handle)
			fasta_split = fasta.split(">")
			for record in fasta_split:
				if len(record) > 0:
					record_split = record.split("\n")
					#check if record id is in list of ids derived from xlink data file
					if record_split[0] in self.proteins.keys():
						sequence = "".join(record_split[1:])
						sequence = sequence.strip("*")
						sequence = sequence.strip()
						self.proteins[record_split[0]] = sequence
						fasta_seqids.append(record_split[0])
					else:
						raise xlinkError("Sequence ID mismatch. Offending value: {} in {}".format(record_split[0], proteins))
			#check for ids in xlink data file that are not in fasta file
			if len(fasta_seqids) < len(self.proteins.keys()):
				mismatches = []
				for seqid in self.proteins.keys():
					if seqid not in fasta_seqids:
						mismatches.append(seqid)
				raise xlinkError("Sequence ID mismatch. Offending value(s): {} in {}".format(mismatches, data))

	def circos(self, colors, focus, linkedres = "K", show_intra = True):
		#create file name strings
		karyo_filename = self.name + "_karyo.txt"
		links_filename = self.name + "_links.txt"
		linkedres_filename   = self.name + "_linkedres.txt"
		conf_filename  = self.name + "_conf.txt"

		#create karyotype file
		with open(karyo_filename, "w") as karyo_handle:
			karyo_counter = 0
			karyo_string = ""
			for key in self.proteins.keys():
				end = str(len(self.proteins[key]) * 2)
				karyo_entry = " ".join(["chr -", key, key, "1", end, colors[karyo_counter]])
				karyo_entry += "\n"
				karyo_string += karyo_entry
				karyo_counter += 1
			karyo_handle.write(karyo_string)

		#create link file
		with open(links_filename, "w") as links_handle:
			links_string = ""
			for link in self.data:
				start1 = str(int(link[1]) * 2 - 1)
				end1   = str(int(link[1]) * 2)
				start2 = str(int(link[3]) * 2 - 1)
				end2   = str(int(link[3]) * 2)
				link_entry = " ".join([link[0], start1, end1, link[2], start2, end2])
				link_entry += "\n"
				links_string += link_entry
			links_handle.write(links_string)

		#create xlinked residues file
		with open(linkedres_filename, "w") as linkedres_handle:
			linkedres_string = ""
			for protid in self.proteins:
				linkedres_index = [index for index, res in enumerate(self.proteins[protid], start = 1) if res == linkedres]
				for index in linkedres_index:
					start = str(index * 2 - 1)
					end = str(index * 2)
					linkedres_entry = " ".join([protid, start, end])
					linkedres_entry += "\n"
					linkedres_string += linkedres_entry
			linkedres_handle.write(linkedres_string)

		#create conf file
		#create color dictionary
		color_dict = dict(zip(self.proteins, colors))
		with open(conf_filename, "w") as conf_handle:
			conf_string = ""
			conf_string += "karyotype = \"" + karyo_filename + "\"\n"
			conf_string += "<links>\n"
			conf_string +=     "<link>\n"
			conf_string +=         "file = \"" + links_filename + "\"\n"
			conf_string +=         "radius = 0.98r\nbezier_radius = 0r\nthickness = 5\n"
			conf_string +=         "<rules>\n"
			if show_intra == True:
				# conf_string +=         "<rule>\ncondition = var(intrachr)\nshow = yes\n</rule>\n"
				counter = 0
				for protid in self.proteins:
					conf_string +=     "<rule>\ncondition = between({},{})\ncolor = {}\n</rule>\n ".format(protid, protid, colors[counter])
					counter += 1
			else:
				conf_string +=         "<rule>\ncondition = var(intrachr)\nshow = no\n</rule>\n"
			protid_pairs = list(combinations(self.proteins, 2))
			inter_xlinks = [(xl[0], xl[2]) for xl in self.data if xl[0] != xl[2]]
			protid_pair_counter = []
			it = 0
			for protid_pair in protid_pairs:
				protid_pair_counter.append([protid_pair, 0])
				for inter_xlink in inter_xlinks:
					if inter_xlink == protid_pair or inter_xlink == (protid_pair[1], protid_pair[0]):
						protid_pair_counter[it][1] += 1
				it += 1
			protid_pair_counter = sorted(protid_pair_counter, key = lambda x: x[1])
			protid_pairs = [counter[0] for counter in protid_pair_counter]
			it = 0
			for protid_pair in protid_pairs:
				conf_string +=         "<rule>\ncondition = between({},{})\ncolor = {}\n</rule>\n".format(protid_pair[0], protid_pair[1], color_dict[protid_pair[0]])
			conf_string +=         "</rules>\n"
			conf_string +=     "</link>\n"
			conf_string += "</links>\n"
			conf_string += "<plots>\n"	
			conf_string +=     "type = tile\nlayers_overflow = hide\n"
			conf_string +=     "<plot>\n"
			conf_string +=         "file = \"" + linkedres_filename + "\"\n"
			conf_string +=         "r1 = 0.98r\nr0 = 0.97r\norientation = out\n"
			conf_string +=         "layers = 15\nmargin = 0.02u\nthickness = 15\npadding = 8\n"
			conf_string +=         "stroke_thickness = 5\nstroke_color = red\n"
			conf_string +=     "</plot>\n"
			conf_string += "</plots>\n"
			conf_string += "show_ticks = yes\nshow_tick_labels = yes\n"
			conf_string += "<ticks>\n"
			conf_string +=     "radius = 1r\ncolor = black\nthickness = 2p\n"
			conf_string +=     "multiplier = 0.5\nformat = %d\n"
			conf_string +=     "<tick>\n"
			conf_string +=         "spacing = 100\nsize = 10p\nshow_label = yes\n"
			conf_string +=         "label_size = 20p\nlabel_offset = 10p\nformat = %d\n"
			conf_string +=     "</tick>\n"
			conf_string += "</ticks>\n"
			conf_string += "<ideogram>\n"
			conf_string +=     "show_label = yes\nlabel_with_tag = yes\nlabel_font = light\n"
			conf_string +=     "label_radius = 1r +100p\nlabel_center = yes\nlabel_size = 48p\n"
			conf_string +=     "label_color = grey\nlabel_parallel = yes\nlabel_case = default\n"
			conf_string +=     "label_format = eval(sprintf(var(label)))\n"
			conf_string +=     "radius = 0.9r\nthickness = 40p\nfill = yes\n"
			conf_string +=     "<spacing>\n"
			conf_string +=         "default = 0.005r\n"
			conf_string +=     "</spacing>\n"
			conf_string += "</ideogram>\n"
			conf_string += "<image>\n<<include etc/image.conf>>\n</image>\n"
			conf_string += "<<include etc/colors_fonts_patterns.conf>>\n"
			conf_string += "<<include etc/housekeeping.conf>>\n"
			conf_handle.write(conf_string)

class xlinkError(Exception):

	def __init__(self, message):
		self.message = message
	def __str__(self):
		return repr(self.message)

			
