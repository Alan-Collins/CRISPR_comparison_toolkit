

def read_array_file(file):
	"""Read array spacers file into dict."""
	array_spacers_dict = {}
	with open(file, 'r') as fin:
		for line in fin.readlines():
			bits = line.split()
			array_spacers_dict[bits[0]] = bits[2:]
	return array_spacers_dict


