
# Compute the number of feet corresponding to a number of miles.

###################################################
# Miles to feet conversion formula
# Student should enter function on the next lines.

def miles_to_feet(miles):
    return 5280 * miles


###################################################
# Tests
# Student should not change this code.

def test(miles):
	print str(miles) + " miles equals",
	print str(miles_to_feet(miles)) + " feet."

test(13)
test(57)
test(82.67)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#13 miles equals 68640 feet.
#57 miles equals 300960 feet.
#82.67 miles equals 436497.6 feet.


# Compute the number of seconds in a given number of hours, minutes, and seconds.

###################################################
# Hours, minutes, and seconds to seconds conversion formula
# Student should enter function on the next lines.

def total_seconds(hours, minutes, seconds):
    return hours * 60 * 60 + minutes * 60 + seconds

###################################################
# Tests
# Student should not change this code.

def test(hours, minutes, seconds):
	print str(hours) + " hours, " + str(minutes) + " minutes, and",
	print str(seconds) + " seconds totals to",
	print str(total_seconds(hours, minutes, seconds)) + " seconds."

test(7, 21, 37)
test(10, 1, 7)
test(1, 0, 1)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.


#7 hours, 21 minutes, and 37 seconds totals to 26497 seconds.
#10 hours, 1 minutes, and 7 seconds totals to 36067 seconds.
#1 hours, 0 minutes, and 1 seconds totals to 3601 seconds.


# Compute the length of a rectangle's perimeter, given its width and height.

###################################################
# Rectangle perimeter formula
# Student should enter function on the next lines.

def rectangle_perimeter(width, height):
    return 2 * width + 2 * height

###################################################
# Tests
# Student should not change this code.

def test(width, height):
	print "A rectangle " + str(width) + " inches wide and " + str(height),
	print "inches high has a perimeter of",
	print str(rectangle_perimeter(width, height)) + " inches."

test(4, 7)
test(7, 4)
test(10, 10)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#A rectangle 4 inches wide and 7 inches high has a perimeter of 22 inches.
#A rectangle 7 inches wide and 4 inches high has a perimeter of 22 inches.
#A rectangle 10 inches wide and 10 inches high has a perimeter of 40 inches.


# Compute the area of a rectangle, given its width and height.

###################################################
# Rectangle area formula
# Student should enter function on the next lines.

def rectangle_area(width, height):
    return width * height


###################################################
# Tests
# Student should not change this code.

def test(width, height):
	print "A rectangle " + str(width) + " inches wide and " + str(height),
	print "inches high has an area of",
	print str(rectangle_area(width, height)) + " square inches."

test(4, 7)
test(7, 4)
test(10, 10)

	
###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#A rectangle 4 inches wide and 7 inches high has an area of 28 square inches.
#A rectangle 7 inches wide and 4 inches high has an area of 28 square inches.
#A rectangle 10 inches wide and 10 inches high has an area of 100 square inches.


# Compute the circumference of a circle, given the length of its radius.

###################################################
# Circle circumference formula
# Student should enter function on the next lines.

import math

def circle_circumference(radius):
    return 2 * math.pi * radius

###################################################
# Tests
# Student should not change this code.

def test(radius):
	print "A circle with a radius of " + str(radius),
	print "inches has a circumference of",
	print str(circle_circumference(radius)) + " inches."

test(8)
test(3)
test(12.9)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#A circle with a radius of 8 inches has a circumference of 50.2654824574 inches.
#A circle with a radius of 3 inches has a circumference of 18.8495559215 inches.
#A circle with a radius of 12.9 inches has a circumference of 81.0530904626 inches.


# Compute the area of a circle, given the length of its radius.

###################################################
# Circle area formula
# Student should enter function on the next lines.

import math

def circle_area(radius):
    return math.pi * (radius ** 2)

###################################################
# Tests
# Student should not change this code.

def test(radius):
	print "A circle with a radius of " + str(radius),
	print "inches has an area of",
	print str(circle_area(radius)) + " square inches."

test(8)
test(3)
test(12.9)

###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#A circle with a radius of 8 inches has an area of 201.06192983 square inches.
#A circle with a radius of 3 inches has an area of 28.2743338823 square inches.
#A circle with a radius of 12.9 inches has an area of 522.792433484 square inches.


# Compute the future value of a given present value, annual rate, and number of years.

###################################################
# Future value formula
# Student should enter function on the next lines.

def future_value(present_value, annual_rate, years):
    return present_value * (1 + 0.01 * annual_rate) ** years

###################################################
# Tests
# Student should not change this code.

def test(present_value, annual_rate, years):
	"""Tests the future_value function."""
	
	print "The future value of $" + str(present_value) + " in " + str(years),
	print "years at an annual rate of " + str(annual_rate) + "% is",
	print "$" + str(future_value(present_value, annual_rate, years)) + "."


###################################################
# Tests
# Student should uncomment ONLY ONE of the following at a time.

test(1000, 7, 10)
test(200, 4, 5)
test(1000, 3, 20)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#The future value of $1000 in 10 years at an annual rate of 7% is $1967.15135729.
#The future value of $200 in 5 years at an annual rate of 4% is $243.33058048.
#The future value of $1000 in 20 years at an annual rate of 3% is $1806.11123467.


# Compute a name tag, given the first and last name.

###################################################
# Name tag formula
# Student should enter function on the next lines.

def name_tag(first_name, last_name):
    return "My name is " + first_name + " " + last_name + "."

###################################################
# Tests
# Student should not change this code.

def test(first_name, last_name):
	print name_tag(first_name, last_name)
	
test("Joe", "Warren")
test("Scott", "Rixner")
test("John", "Greiner")


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#My name is Joe Warren.
#My name is Scott Rixner.
#My name is John Greiner.


# Compute the statement about a person's name and age, given the person's name and age.

###################################################
# Name and age formula
# Student should enter function on the next lines.

def name_and_age(name, age):
    return name + " is " + str(age) + " years old."

###################################################
# Tests
# Student should not change this code.

def test(name, age):
	print name_and_age(name, age)
	
test("Joe Warren", 52)
test("Scott Rixner", 40)
test("John Greiner", 46)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#Joe Warren is 52 years old.
#Scott Rixner is 40 years old.
#John Greiner is 46 years old.


# Compute the distance between the points (x0, y0) and (x1, y1).

###################################################
# Distance formula
# Student should enter function on the next lines.

# Hint: You need to use the power operation ** .

def point_distance(x0, y0, x1, y1):
    return sqrt((x0 - x1)**2 + (y0 - y1)**2)

###################################################
# Tests
# Student should not change this code.

def test(x0, y0, x1, y1):
	print "The distance from (" + str(x0) + ", " + str(y0) + ") to",
	print "(" + str(x1) + ", " + str(y1) + ") is",
	print str(point_distance(x0, y0, x1, y1)) + "."

test(2, 2, 5, 6)
test(1, 1, 2, 2)
test(0, 0, 3, 4)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#The distance from (2, 2) to (5, 6) is 5.0.
#The distance from (1, 1) to (2, 2) is 1.41421356237.
#The distance from (0, 0) to (3, 4) is 5.0.


# Compute the area of a triangle (using Heron's formula),
# given its side lengths.

###################################################
# Triangle area (Heron's) formula
# Student should enter function on the next lines.
# Hint:  Also define point_distance as use it as a helper function.

def triangle_area(x0, y0, x1, y1, x2, y2):
    a = point_distance(x0, y0, x1, y1)
    b = point_distance(x1, y1, x2, y2)
    c = point_distance(x2, y2, x0, y0)
    s = (a + b + c) / 2
    return sqrt(s * (s - a) * (s - b) * (s - c))

###################################################
# Tests
# Student should not change this code.

def test(x0, y0, x1, y1, x2, y2):
	print "A triangle with vertices (" + str(x0) + "," + str(y0) + "),",
	print "(" + str(x1) + "," + str(y1) + "), and",
	print "(" + str(x2) + "," + str(y2) + ") has an area of",
	print str(triangle_area(x0, y0, x1, y1, x2, y2)) + "."

test(0, 0, 3, 4, 1, 1)
test(-2, 4, 1, 6, 2, 1)
test(10, 0, 0, 0, 0, 10)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#A triangle with vertices (0,0), (3,4), and (1,1) has an area of 0.5.
#A triangle with vertices (-2,4), (1,6), and (2,1) has an area of 8.5.
#A triangle with vertices (10,0), (0,0), and (0,10) has an area of 50.


# Compute and print tens and ones digit of an integer in [0,100).

###################################################
# Digits function
# Student should enter function on the next lines.

def print_digits(number):
    tens = number // 10
    ones = number % 10
    print "The tens digit is " + str(tens) + ",", 
    print "and the ones digit is " + str(ones) + "."

	
###################################################
# Tests
# Student should not change this code.
	
print_digits(42)
print_digits(99)
print_digits(5)


###################################################
# Expected output
# Student should look at the following comments and compare to printed output.

#The tens digit is 4, and the ones digit is 2.
#The tens digit is 9, and the ones digit is 9.
#The tens digit is 0, and the ones digit is 5.


# Compute and print powerball numbers.

###################################################
# Powerball function
# Student should enter function on the next lines.

import random

def different_numbers(no1, no2, no3, no4, no5, no6, no7):
    number1 = ((no1 != no2) and 
               (no1 != no3) and 
               (no1 != no4) and 
               (no1 != no5) and 
               (no1 != no6) and 
               (no1 != no7)) 
    number2 = ((no2 != no3) and 
               (no2 != no4) and 
               (no2 != no5) and 
               (no2 != no6) and 
               (no2 != no7))
    number3 = ((no3 != no4) and 
               (no3 != no5) and 
               (no3 != no6) and  
               (no3 != no7))
    number4 = ((no4 != no5) and 
               (no4 != no6) and 
               (no4 != no7))
    number5 = ((no5 != no6) and 
               (no5 != no7))
    number6 = (no6 != no7)
    
    return (number1 and number2 and number3 and
            number4 and number5 and number6)
    
def powerball():
    while True:
        ball_one = random.randrange(1, 60)
        ball_two = random.randrange(1, 60)
        ball_three = random.randrange(1, 60)
        ball_four = random.randrange(1, 60)
        ball_five = random.randrange(1, 60)
        ball_six = random.randrange(1, 60)
        ball_powerball = random.randrange(1, 36)
        
        if different_numbers(ball_one, ball_two, ball_three, 
                             ball_four, ball_five, ball_six,
                             ball_powerball):
            break

    print "Today's numbers are",
    print str(ball_one) + ",",
    print str(ball_two) + ",",
    print str(ball_three) + ",",
    print str(ball_four) + ",",
    print str(ball_five) + ",",
    print str(ball_six) + ".",
    print "The Powerball number is",
    print str(ball_powerball) + "."

	
###################################################
# Tests
# Student should not change this code.
	
powerball()
powerball()
powerball()




