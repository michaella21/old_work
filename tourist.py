"""
There is a lost tourist in Grandeville. The streets in Grandeville run east to west and go from
..., S. 2nd St., S. 1st St., Broadway St., N. 1st St., N. 2nd St., ...
The avenues run north to south and go from
..., E. 2nd Ave., E. 1st Ave., Broadway Ave., W. 1st Ave., W. 2nd Ave., ...
These streets form a square block grid. For each of the questions below, the tourist starts at the intersection of Broadway St. and Broadway Avenue and moves one block in each of the four cardinal directions with equal probability.

Q1 What is the probability that the tourist is at least 3 city blocks (as the crow flies) from Broadway and Broadway after 10 moves?
0.4539794922

 
Q2 What is the probability that the tourist is at least 10 city blocks (as the crow flies) from Broadway and Broadway after 60 moves?
 0.2065011181


Q3 What is the probability that the tourist is ever at least 5 city blocks (as the crow flies) from Broadway and Broadway within 10 moves?
 0.2156982422


Q4 What is the probability that the tourist is ever at least 10 city blocks (as the crow flies) from Broadway and Broadway within 60 moves?
 0.6252450324


Q5 What is the probability that the tourist is ever east of East 1st Avenue but ends up west of West 1st Avenue in 10 moves?
 0.8883504868


Q6 What is the probability that the tourist is ever east of East 1st Avenue but ends up west of West 1st Avenue in 30 moves?
 0.9996428368


Q7 What is the average number of moves until the first time the tourist is at least 10 city blocks (as the crow flies) from Broadway and Broadway.
 
Q8 What is the average number of moves until the first time the tourist is at least 60 city blocks (as the crow flies) from Broadway and Broadway.
"""

from sympy import *
import numpy as np
import math
from decimal import Decimal


"""
Let the intersection of Broadway and Broadway be an origin (0,0). 
Let the location of the trourist after 'n' moves be (X_n, Y_n). (thus (X_0, Y_0)=(0,0)). 
Every time at the 'n'th move, the tourist chooses 'a' horizontal moves and 'n-a' vertical moves. 
Among 'a' horizontal moves, 'k' moves to east(-1) and 'a-k' moves to west(+1).
Among 'n-a' vertical moves, 'r' movest to south(-1) and 'n-a-r' moves to north(+1).

"""

#compute the probability for a tourist to get (x,y) after n moves from (0,0)
def prob_to(x,y,n):
"""(n choose a)*(0.5)**n * (a choose k) * (0.5)** a * (n-a choose r) * (0.5)**(n-a)
=(n!/(k!(a-k)!r!(n-a-r)!)) * (0.25)**n
where a-2k=x, n-a-2r=y """
prob=0
for i in range(n+1):#from 0 moves to n moves
if ((i-x)/2.0).is_integer() and ((n-i-y)/2.0).is_integer():
k=(i-x)/2.0
r=(n-i-y)/2.0
if i>=k>=0 and n-i>=r>=0:
prob+= (0.25)**n*(factorial(n))/(math.factorial(k)*math.factorial(i-k)*math.factorial(r)*math.factorial(n-i-r))
else:
prob+=0
else: 
prob+=0
return prob 
#compute the probability for a tourist to visit (c,d) from (a,b) in exact n moves
def prob_from_to(a,b,c,d,n):
#moving from (a,b) to (c,d) is EQUIVALENT to moving from (0,0) to (c-a, d-b) 
return prob_to(c-a, d-b, n)

#find coordinates within m bloks (up to (m-1) blocks) from (0,0)=(Broadway, Broadway), excluding m, 

def within_blocks(m):
coordinates=[]
for i in range(int(-1*m), int(m+1)):#from -m blocks to m blocks
for j in range(int(-1*m), int(m+1)):
if np.sqrt(i**2+j**2)<m:
coordinates.append((i,j))
return coordinates
#whether it is reachable in one step from (a,b) to (x,y)
def one_step(a,b,x,y):
if (a-x)*(b-y)==0 and abs(a-x)+abs(b-y)==1:#can only choose either vertical or horizontal, and only one step
return True
else:
return False

#Q1. After 10 moves, i.e., n=10, farther than or equal to 3 blocks; 1-prob(closer than 3 (exlcusive) blocks afer 10 moves))
ans1=0
for (i,j) in within_blocks(3):
ans1+=prob_to(i,j,10)
print round(1-ans1, 10) 

#Q2. After 60 moves, i.e., n=60, farther than or equal to 10 blocks, 1-prob(closer than 10 (exlcusive) blocks afer 60 moves))
ans2=0
for (i,j) in within_blocks(10):
ans2+=prob_to(i,j,60)
print round(1-ans2, 10)

#Q3. Prob(one can ever be away for 5 or more blocks from (0,0) within 10 moves)
# = sum{Prob(one reaches at 5 blocks away from(0,0) exactly after k moves|one only stays within 4 blocks until (k-1)moves)\
# * Prob(only stays within 4 blocsk until (k-1) moves)} from k=1 to 10

#first list (x,y)'s that are exactly 5 blocks away from (0,0) 
target_5=list(set(within_blocks(6))^set(within_blocks(5)))
one_step_to_target5=[] # find all pairs of two consecutive locations that one can reach from 4 blocks away to 5 blocks away in one step

for (i,j) in target_5: 
for (a,b) in within_blocks(5):
if one_step(a,b,i,j):
one_step_to_target5.append((a,b,i,j))
ans3=0
for (a,b,c,d) in one_step_to_target5:
for i in range(11):#from 0 to 10 moves
ans3+=0.25*prob_to(a,b,i)

print round(ans3, 10) 

#Q4. similar to Q3, first ever reach to 10 blocks within 60 moves
target_10=list(set(within_blocks(11))^set(within_blocks(10)))
one_step_to_target10=[]

for (i,j) in target_10: 
for (a,b) in within_blocks(10):
if one_step(a,b,i,j):
one_step_to_target10.append((a,b,i,j))
ans4=0
for (a,b,c,d) in one_step_to_target10:
for i in range(61):#from 0 to 60 moves
ans4+=0.25*prob_to(a,b,i)

print round(ans4, 10)


#Q5.starting from (0,0), visit (1,y) and then visit (-1,y') in 10 moves 
# 1-Pr(never moves to east for 10 moves)-Pr(never moves to west for 10 moves)
# -Pr(moves only vertically for 10 moves)

never_east10=0
for i in range(1,11): #can go to only east for 1 step to 10 steps
for j in range(11-i): #rest of moves for vertical moves
never_east10+=\
(0.25)**10*(math.factorial(10)/(math.factorial(i)*math.factorial(j)*math.factorial(10-i-j)))

#pr(never moves to east)= pr(never moves to west), symmetric
never_west10=never_east10

only_vertical10=0
for i in range(1,11):
only_vertical10+=binomial(10,i)*(0.25)**10

print round(1-never_east10-never_west10-only_vertical10 ,10) 


#Q6, similar to Q5
never_east30=0
for i in range(1,31): #can go to only east for 1 step to 10 steps
for j in range(31-i): #rest of moves for vertical moves
never_east30+=\
(0.25)**30*(math.factorial(30)/(math.factorial(i)*math.factorial(j)*math.factorial(30-i-j)))

#pr(never moves to east)= pr(never moves to west), symmetric
never_west30=never_east30

only_vertical30=0
for i in range(1,31):
only_vertical30+=binomial(30,i)*(0.25)**30

print round(1-never_east30-never_west30-only_vertical30 ,10) 



