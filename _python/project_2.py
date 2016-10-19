
# Testing template for name_to_number()

###################################################
# Copy and paste your definition of name_to_number() here

def name_to_number(name):
    if(name == "rock"):
        return 0
    elif(name == "Spock"):
        return 1
    elif(name == "paper"):
        return 2
    elif(name == "lizard"):
        return 3
    elif(name == "scissors"):
        return 4
    else:
        pass
  #      print "Sorry, your choice is not one of",
  #      print "rock, paper, scissors, lizard or Spock"

###################################################
# Test calls to name_to_number()
print name_to_number("rock")
print name_to_number("Spock")
print name_to_number("paper")
print name_to_number("lizard")
print name_to_number("scissors")
print name_to_number("scissfors")


###################################################
# Output to the console should have the form:
# 0
# 1
# 2
# 3
# 4

# Testing template for number_to_name()

###################################################
# Copy and paste your definition of number_to_name() here

def number_to_name(number):
    if(number == 0):
        return "rock"
    elif(number == 1):
        return "Spock"
    elif(number == 2):
        return "paper"
    elif(number == 3):
        return "lizard"
    elif(number == 4):
        return "scissors"
    else:
        pass

###################################################
# Test calls to number_to_name()
print number_to_name(0)
print number_to_name(1)
print number_to_name(2)
print number_to_name(3)
print number_to_name(4)
print number_to_name(5)


###################################################
# Output to the console should have the form:
# rock
# Spock
# paper
# lizard
# scissors



def verb(winning_choice, difference):
    if(winning_choice == "Spock"):
        if(abs(difference) == 3):
            return " smashes "
        elif(abs(difference) == 1):
            return " vaporises "
        
    elif(winning_choice == "lizard"):
        if(abs(difference) == 2):
            return " poisons "
        elif(abs(difference) == 1):
            return " eats "
        
    elif(winning_choice == "paper"):
        if(abs(difference) == 2):
            return " covers "
        elif(abs(difference) == 1):
            return " disproves "  
        
    elif(winning_choice == "rock"):
        if(abs(difference) == 4):
            return " blunts "
        elif(abs(difference) == 3):
            return " crushes "     
    elif(winning_choice == "scissors"):
        if(abs(difference) == 1):
            return " decapitate "
        elif(abs(difference) == 2):
            return " cuts "           


# Testing template for rpsls(player_choice)

###################################################

import random

def rpsls(player_choice):
    print
    print "You have chosen " + player_choice
    player_number = name_to_number(player_choice)
#    print "This is number " + str(player_number)
    
    if(player_number != None): 
        computer_number = random.randrange(0,5)
        computer_choice = number_to_name(computer_number)
        print "I have chosen " + computer_choice
#        print "This is number " + str(computer_number)
        difference = computer_number - player_number
        remainder = difference % 5
    
        if(remainder == 0):
            print "It's a draw!"
        elif(remainder == 1 or remainder == 2):
            print "I win -", 
            winning_verb = verb(computer_choice, difference)
            print computer_choice + winning_verb + player_choice + "!"
        elif (remainder == 3 or remainder == 4):
            print "You win -",
            winning_verb = verb(player_choice, difference)
            print player_choice + winning_verb + computer_choice + "!"
    else:
        print "Sorry, your choice was not one of",
        print "rock, paper, scissors, lizard or Spock."
        print "Please try again."
###################################################
# Test calls to number_to_name()
rpsls("rock")
rpsls("Spock")
rpsls("paper")
rpsls("lizard")
rpsls("scissors")
rpsls("scissfors")

###################################################
# Output to the console should have the form:

#You have chosen rock
#I have chosen ???
# You/I win!

#You have chosen Spock
#I have chosen ???
# You/I win!

#You have chosen paper
#I have chosen ???
# You/I win!

#You have chosen lizard
#I have chosen ???
# You/I win!

#You have chosen scissors
#I have chosen ???
# You/I win!


print 4 % 5 # player wins = 4

print 3 % 5 # player wins = 3
#Spock smashes scissors

print 2 % 5 # computer wins = 2
#Lizard poisons Spock

print 1 % 5 # computer wins = 1
#Spock vaporises rock
#Paper disproves Spock

print 0 % 5 #draw

print -1 % 5 # player wins = 4
#Spock vaporises rock
#Paper disproves Spock

print -2 % 5 # player wins = 3
#Lizard poisons Spock

print -3 % 5 # computer wins = 2 
#Spock smashes scissors

print -4 % 5 #computer wins = 1 

