def generate_file(header='>hw1_W_2_chr_1',CNV=[],STR=[],INV=[],INS=[],DEL=[],SNP=[],ALU=[]):
	"""Creates output file"""
	file = open('answer.txt','w')
	file.write("{}\n".format(header))

	file.write('>CNV\n')
	for item in CNV:
		file.write("{}\n".format(item))

	file.write('>STR\n')
	for item in STR:
		file.write("{}\n".format(item))

	file.write('>INV\n')
	for item in INV:
		file.write("{}\n".format(item))

	file.write('>INS\n')
	for item in INS:
		file.write("{}\n".format(item))

	file.write('>DEL\n')
	for item in DEL:
		file.write("{}\n".format(item))

	file.write('>SNP\n')
	for item in SNP:
		file.write("{}\n".format(item))

	file.write('>ALU\n')
	for item in ALU:
		file.write("{}\n".format(item))

	file.close()