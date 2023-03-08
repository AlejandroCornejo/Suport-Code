import matplotlib.pyplot as plt
import math
import numpy

###########################################
def ReadDataFromFile(name_of_the_file=""):
	x = []
	y = []
	file = open(name_of_the_file, "r")
	contents = file.readlines()
	counter = 1
	for line in contents:
		numbers = line.split("\t")
		counter += 1
		x.append(float(numbers[0]))
		y.append(float(numbers[1]))
	file.close()
	return x, y
###########################################

"""
Here we define the files to print
"""
## read first file data
x1, y1 = ReadDataFromFile("data/2w_1e-3.txt")
plt.plot(x1, y1, "g", label = "2 Way, dt=1e-3 s", linewidth=2)

x2, y2 = ReadDataFromFile("data/2w_1e-4.txt")
plt.plot(x2, y2, "blue", label = "2 Way, dt=1e-4 s", linewidth=2)

x3, y3 = ReadDataFromFile("data/2w_1e-5.txt")
plt.plot(x3, y3, "cyan", label = "2 Way, dt=1e-5 s", linewidth=2)

x3, y3 = ReadDataFromFile("data/2w_5e-6.txt")
plt.plot(x3, y3, "black", label = "2 Way, dt=5e-6 s", linewidth=2)

## read first file data
x1, y1 = ReadDataFromFile("data/1w_1e-3.txt")
plt.plot(x1, y1, "g+", label = "1 Way, dt=1e-3 s", linewidth=2)

x2, y2 = ReadDataFromFile("data/1w_1e-4.txt")
plt.plot(x2, y2, "b+", label = "1 Way, dt=1e-4 s", linewidth=2)

x3, y3 = ReadDataFromFile("data/1w_1e-5.txt")
plt.plot(x3, y3, "c+", label = "1 Way, dt=1e-5 s", linewidth=2)

x3, y3 = ReadDataFromFile("data/1w_5e-6.txt")
plt.plot(x3, y3, "k+", label = "1 Way, dt=5e-6 s", linewidth=2)

# markersize=10

# Constant values
# x = numpy.arange(10,20,0.2)
# y = 80*numpy.ones(len(x), dtype=float)
# plt.plot(x, y, "b-", label = "ref_value", linewidth=2)


# Error bar
#xe, ye = ReadDataFromFile("error.txt")
#plt.errorbar(xe, ye  , label="error_bar", yerr=0.5e4, ecolor="lightgrey", errorevery=1,
    #fmt='grey', capsize=5)


#####################################################################
#####################################################################
"""
Here we set the axes and limits
"""
x_max = 0.22
x_min = 0.1
y_max = -5e-2
y_min = -6.4e-2
x_interval = 0.02
y_interval = 0.2e-2

x_ticks = numpy.arange(x_min, x_max, x_interval)
plt.xticks(x_ticks,  fontsize = 12)

y_ticks = numpy.arange(y_min, y_max, y_interval)
plt.yticks(y_ticks, fontsize = 12)

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


#####################################################################
#####################################################################

"""
Here we put the labels and titles of the plot
"""
plt.xlabel('Time [s]', fontsize = 12)
plt.ylabel('Displacement [m]', fontsize = 12)
#plt.title('Three point bending', fontsize = 12)
plt.legend(fontsize = 13)
plt.grid()
plt.show()






#####################################################################
#####################################################################
# Settings
"""
'.'	point marker
','	pixel marker
'o'	circle marker
'v'	triangle_down marker
'^'	triangle_up marker
'<'	triangle_left marker
'>'	triangle_right marker
'1'	tri_down marker
'2'	tri_up marker
'3'	tri_left marker
'4'	tri_right marker
's'	square marker
'p'	pentagon marker
'*'	star marker
'h'	hexagon1 marker
'H'	hexagon2 marker
'+'	plus marker
'x'	x marker
'D'	diamond marker
'd'	thin_diamond marker
'|'	vline marker
'_'	hline marker
"""

"""
'-'	solid line style
'--'	dashed line style
'-.'	dash-dot line style
':'	dotted line style
"""

"""
'b'	blue
'g'	green
'r'	red
'c'	cyan
'm'	magenta
'y'	yellow
'k'	black
'w'	white
"""

#####################################################################
#####################################################################