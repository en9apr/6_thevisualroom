
# Compute the number of feet corresponding to a number of miles.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#miles = 13


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#miles = 57


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
miles = 82.67


###################################################
# Miles to feet conversion formula
# Student should enter formula on the next line.

feet = 5280 * miles

###################################################
# Test output
# Student should not change this code.

print str(miles) + " miles equals " + str(feet) + " feet."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#13 miles equals 68640 feet.

# Test 2 output:
#57 miles equals 300960 feet.

# Test 3 output:
#82.67 miles equals 436497.6 feet.



# Compute the number of seconds in a given number of hours, minutes, and seconds.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#hours = 7
#minutes = 21
#seconds = 37


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#hours = 10
#minutes = 1
#seconds = 7


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
hours = 1
minutes = 0
seconds = 1


###################################################
# Hours, minutes, and seconds to seconds conversion formula
# Student should enter formula on the next line.

total_seconds = hours * 60 * 60 + minutes * 60 + seconds

###################################################
# Test output
# Student should not change this code.

print str(hours) + " hours, " + str(minutes) + " minutes, and",
print str(seconds) + " seconds totals to " + str(total_seconds) + " seconds."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
# Test 1 output:
#7 hours, 21 minutes, and 37 seconds totals to 26497 seconds.

# Test 2 output:
#10 hours, 1 minutes, and 7 seconds totals to 36067 seconds.

# Test 3 output:
#1 hours, 0 minutes, and 1 seconds totals to 3601 seconds.



# Compute the length of a rectangle's perimeter, given its width and height.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#width = 4
#height = 7


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#width = 7
#height = 4


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
width = 10
height = 10


###################################################
# Rectangle perimeter formula
# Student should enter formula on the next line.

perimeter = 2 * width + 2 * height

###################################################
# Test output
# Student should not change this code.

print "A rectangle " + str(width) + " inches wide and " + str(height),
print "inches high has a perimeter of " + str(perimeter) + " inches."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#A rectangle 4 inches wide and 7 inches high has a perimeter of 22 inches.

# Test 2 output:
#A rectangle 7 inches wide and 4 inches high has a perimeter of 22 inches.

# Test 3 output:
#A rectangle 10 inches wide and 10 inches high has a perimeter of 40 inches.



# Compute the area of a rectangle, given its width and height.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#width = 4
#height = 7


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#width = 7
#height = 4


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
width = 10
height = 10


###################################################
# Rectangle area formula
# Student should enter formula on the next line.

area = width * height

###################################################
# Test output
# Student should not change this code.

print "A rectangle " + str(width) + " inches wide and " + str(height),
print "inches high has an area of " + str(area) + " square inches."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#A rectangle 4 inches wide and 7 inches high has an area of 28 square inches.

# Test 2 output:
#A rectangle 7 inches wide and 4 inches high has an area of 28 square inches.

# Test 3 output:
#A rectangle 10 inches wide and 10 inches high has an area of 100 square inches.



# Compute the circumference of a circle, given the length of its radius.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.
PI = 3.14

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#radius = 8


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#radius = 3


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
radius = 12.9


###################################################
# Circle circumference formula
# Student should enter formula on the next line.

circumference = 2 * PI * radius 

###################################################
# Test output
# Student should not change this code.

print "A circle with a radius of " + str(radius),
print "inches has a circumference of " + str(circumference) + " inches."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#A circle with a radius of 8 inches has a circumference of 50.24 inches.

# Test 2 output:
#A circle with a radius of 3 inches has a circumference of 18.84 inches.

# Test 3 output:
#A circle with a radius of 12.9 inches has a circumference of 81.012 inches.



# Compute the area of a circle, given the length of its radius.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.
PI = 3.14

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#radius = 8


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#radius = 3


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
radius = 12.9


###################################################
# Circle area formula
# Student should enter formula on the next line.

area = PI * radius ** 2

###################################################
# Test output
# Student should not change this code.

print "A circle with a radius of " + str(radius),
print "inches has an area of " + str(area) + " square inches."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#A circle with a radius of 8 inches has an area of 200.96 square inches.

# Test 2 output:
#A circle with a radius of 3 inches has an area of 28.26 square inches.

# Test 3 output:
#A circle with a radius of 12.9 inches has an area of 522.5274 square inches.



# Compute the future value of a given present value, annual rate, and number of years.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#present_value = 1000
#annual_rate = 7
#years = 10


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#present_value = 200
#annual_rate = 4
#years = 5


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
present_value = 1000
annual_rate = 3
years = 20


###################################################
# Future value formula
# Student should enter formula on the next line.

future_value = present_value * (1 + 0.01 * annual_rate) ** years

###################################################
# Test output
# Student should not change this code.

print "The future value of $" + str(present_value) + " in " + str(years),
print "years at an annual rate of " + str(annual_rate) + "% is $" + str(future_value) + "."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#The future value of $1000 in 10 years at an annual rate of 7% is $1967.15135729.

# Test 2 output:
#The future value of $200 in 5 years at an annual rate of 4% is $243.33058048.

# Test 3 output:
#The future value of $1000 in 20 years at an annual rate of 3% is $1806.11123467.



# Compute a name tag, given the first and last name.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#first_name = "Joe"
#last_name = "Warren"


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#first_name = "Scott"
#last_name = "Rixner"


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
first_name = "John"
last_name = "Greiner"


###################################################
# Name tag formula
# Student should enter formula on the next line.

name_tag = "My name is " + first_name + " " + last_name

###################################################
# Test output
# Student should not change this code.

print name_tag


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#My name is Joe Warren.

# Test 2 output:
#My name is Scott Rixner.

# Test 3 output:
#My name is John Greiner.



# Compute the statement about a person's name and age, given the person's name and age.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#name = "Joe Warren"
#age = 52


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#name = "Scott Rixner"
#age = 40


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
name = "John Greiner"
age = 46


###################################################
# Name and age formula
# Student should enter formula on the next line.

statement = name + " is " + str(age) + " years old"

###################################################
# Test output
# Student should not change this code.

print statement

###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#Joe Warren is 52 years old.

# Test 2 output:
#Scott Rixner is 40 years old.

# Test 3 output:
#John Greiner is 46 years old.



# Compute the distance between the points (x0, y0) and (x1, y1).

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#x0 = 2
#y0 = 2
#x1 = 5
#y1 = 6


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#x0 = 1
#y0 = 1
#x1 = 2
#y1 = 2


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
x0 = 0
y0 = 0
x1 = 3
y1 = 4


###################################################
# Distance formula
# Student should enter formula on the next line.

# Hint: You need to use the power operation ** .

distance = ((x0 - x1)**2+(y0 - y1)**2)**0.5


###################################################
# Test output
# Student should not change this code.

print "The distance from (" + str(x0) + ", " + str(y0) + ") to", 
print "(" + str(x1) + ", " + str(y1) + ") is " + str(distance) + "."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#The distance from (2, 2) to (5, 6) is 5.0.

# Test 2 output:
#The distance from (1, 1) to (2, 2) is 1.41421356237.

# Test 3 output:
#The distance from (0, 0) to (3, 4) is 5.0.


# Compute the area of a triangle (using Heron's formula),
# given its side lengths.

###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

# Test 1 - Select the following lines and use ctrl+shift+k to uncomment.
#x0, y0 = 0, 0
#x1, y1 = 3, 4
#x2, y2 = 1, 1


# Test 2 - Select the following lines and use ctrl+shift+k to uncomment.
#x0, y0 = -2, 4
#x1, y1 = 1, 6
#x2, y2 = 2, 1


# Test 3 - Select the following lines and use ctrl+shift+k to uncomment.
x0, y0 = 10, 0
x1, y1 = 0, 0
x2, y2 = 0, 10


###################################################
# Triangle area (Heron's) formula
# Student should enter formulas on the next lines.

a = ((x0 - x1)**2+(y0 - y1)**2)**0.5
b = ((x0 - x2)**2+(y0 - y2)**2)**0.5
c = ((x1 - x2)**2+(y1 - y2)**2)**0.5

s = 0.5 * (a + b + c)

area = (s * (s - a) * (s - b) * (s - c))**0.5

###################################################
# Test output
# Student should not change this code.

print "A triangle with vertices (" + str(x0) + "," + str(y0) + "),",
print "(" + str(x1) + "," + str(y1) + "), and",
print "(" + str(x2) + "," + str(y2) + ") has an area of " + str(area) + "."


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

# Test 1 output:
#A triangle with vertices (0,0), (3,4), and (1,1) has an area of 0.5.

# Test 2 output:
#A triangle with vertices (-2,4), (1,6), and (2,1) has an area of 8.5.

# Test 3 output:
#A triangle with vertices (10,0), (0,0), and (0,10) has an area of 50.




