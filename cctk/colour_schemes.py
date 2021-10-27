Cols_tol = [
	"#332288",
	"#117733",
	"#44AA99",
	"#88CCEE",
	"#DDCC77",
	"#CC6677",
	"#AA4499",
	"#882255"
	]

Cols_hex_12 = [
	"#07001c",
	"#ff6f8d",
	"#4c62ff",
	"#92ffa9",
	"#810087",
	"#bcffe6",
	"#490046",
	"#00c8ee",
	"#b53900",
	"#ff8cf7",
	"#5b5800",
	"#14d625"
	]

Cols_hex_27 = [
	"#fd5925",
	"#dbc58e",
	"#008d40",
	"#304865",
	"#934270",
	"#f7b8a2",
	"#907500",
	"#45deb2",
	"#1f4195",
	"#d67381",
	"#8e7166",
	"#afb200",
	"#005746",
	"#a598ff",
	"#8f0f1b",
	"#b96000",
	"#667f42",
	"#00c7ce",
	"#9650f0",
	"#614017",
	"#59c300",
	"#1a8298",
	"#b5a6bd",
	"#ea9b00",
	"#bbcbb3",
	"#00b0ff",
	"#cd6ec6"
	]


def choose_col_scheme(ncolours):
	if ncolours > 8:
		if ncolours > 12: 
			if ncolours > 27:
				if ncolours > 40:
					if ncolours < 65:
						col_scheme = Cols_tol
					elif ncolours < 145:
						col_scheme = Cols_hex_12
					else:
						col_scheme = Cols_hex_27
					colours = []
					for i in range((ncolours+len(col_scheme)-1)//len(col_scheme)): # Repeat the same colour scheme.
						for j in col_scheme:
							colours += [(j, col_scheme[i])]

				else:
					colours = [(i, "#000000") for i in Cols_hex_40]
			else:
				colours = [(i, "#000000") for i in Cols_hex_27]
		else:
			colours = [(i, "#000000") for i in Cols_hex_12]
	else:
		colours = [(i, "#000000") for i in Cols_tol]
	return colours