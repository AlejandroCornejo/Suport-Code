import matplotlib.pyplot as plt
import math
import numpy
import scienceplots

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
#plt.rcParams.update({'figure.dpi': '100'})
plt.style.use(['science', 'high-vis']) # science  nature ieee // vibrant high-vis

## read first file data
x1, y1 = ReadDataFromFile("plot.txt")
plt.plot(x1, y1, label = "2 Way, dt=1e-3 s", linewidth=2)

x2, y2 = ReadDataFromFile("plot_2.txt")
plt.plot(x2, y2, label = "balcayde", linewidth=2)


#####################################################################
#####################################################################
"""
Here we set the axes and limits
"""
x_max = 8.001e-7
x_min = -8.001e-7
y_max = 0.06
y_min = -0.06
x_interval = 2.0e-7
y_interval = 0.02

x_ticks = numpy.arange(x_min, x_max, x_interval)
plt.xticks(x_ticks,  fontsize = 12)

y_ticks = numpy.arange(y_min, y_max, y_interval)
plt.yticks(y_ticks, fontsize = 12)

plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


#####################################################################
#####################################################################

"""
Here we put the labels and titles of the plot
"""
plt.xlabel('Displacement [m]', fontsize = 12)
plt.ylabel('Horizontal reaction [N]', fontsize = 12)
plt.title('Cyclic shear test', fontsize = 10)
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