from random import randint
from random import random as r
from random import uniform
from random import sample
from random import gauss
import math
from sympy import symbols, Eq, solve
import time

number_of_generation=0
number_of_population=0
number_of_children=0
benchmark_function=0
n=0
C=0
sigma=0
f2_a=0
f2_b=0
f2_c=0
sigma_modifier=0  # oneFifth_Success_Rule (1)  OR  self_adaption (2)
xover_number=0
alpha=0
selection_method_number=0
tournament_k=0
scaling_number=0
pawer_Law_scaling_k=0
boltzmann_scaling_T=0
Sigma_Cutting_scaling_c=0
Linear_scaling_c=0
number_of_scaling=0


def read_from_file():
  global number_of_population, number_of_generation, number_of_children, benchmark_function, sigma_modifier, n, sigma, C, f2_a, f2_b, f2_c, xover_number, alpha, selection_method_number, tournament_k, scaling_number, pawer_Law_scaling_k, boltzmann_scaling_T, Sigma_Cutting_scaling_c, Linear_scaling_c, number_of_scaling
  reading_file = open("input.txt",'r')
  number_of_generation=int(reading_file.readline().split("=")[1])
  number_of_population=int(reading_file.readline().split("=")[1])
  number_of_children=int(reading_file.readline().split("=")[1])
  n=int(reading_file.readline().split("=")[1])
  sigma=float(reading_file.readline().split("=")[1])
  C=float(reading_file.readline().split("=")[1])
  for i in range(4):
    reading_file.readline()
  benchmark_function=int(reading_file.readline().split(":")[1])
  reading_file.readline()
  if benchmark_function==2:
    f2_a=float(reading_file.readline().split("=")[1])
    f2_b=float(reading_file.readline().split("=")[1])
    f2_c=float(reading_file.readline().split("=")[1])
    for i in range(3):
      reading_file.readline()
  else:
    for i in range(6):
      reading_file.readline()
  sigma_modifier=int(reading_file.readline().split(":")[1])
  for i in range(4):
    reading_file.readline()
  xover_number=int(reading_file.readline().split(":")[1])
  alpha=float(reading_file.readline().split("=")[1])
  for i in range(4):
    reading_file.readline()
  selection_method_number=int(reading_file.readline().split(":")[1])
  reading_file.readline()
  if selection_method_number==3:
    tournament_k=int(reading_file.readline().split("=")[1])
    for i in range(7):
      reading_file.readline()
  else:
    for i in range(8):
      reading_file.readline()
  scaling_number=int(reading_file.readline().split(":")[1])
  reading_file.readline()
  pawer_Law_scaling_k=float(reading_file.readline().split("=")[1])
  reading_file.readline()
  boltzmann_scaling_T=float(reading_file.readline().split("=")[1])
  reading_file.readline()
  Sigma_Cutting_scaling_c=float(reading_file.readline().split("=")[1])
  reading_file.readline()
  Linear_scaling_c=float(reading_file.readline().split("=")[1])
  reading_file.readline()
  number_of_scaling=int(reading_file.readline().split("=")[1])
  reading_file.close()

def f1(x):
  total=0
  for i in range(n):
    total = total + ( (x[i]**2) - 10*math.cos(2*math.pi*x[i]) )
  f = 10*n + total
  return f

def f2(x):
  total1=0
  total2=0
  for i in range(n):
    total1 = total1 + (x[i]**2)
  for i in range(n):
    total2 = total2 + math.cos(f2_c*x[i])
  f = ( (-1*f2_a) * math.exp((-1*f2_b) * math.sqrt((1/n) * total1)) ) - ( math.exp((1/n) * total2) ) + f2_a + math.exp(1)
  return f

def f3(x):
  total=0
  for i in range(n):
    if x[i]>5.12 or x[i]<-5.12:
      total = total + 10*(x[i]**2)
    if x[i]>=-5.12 and x[i]<=5.12:
      total = total + ( (x[i]**2) - 10*math.cos(2*math.pi*x[i]) )
  f = 10*n + total
  return f

def fitnessFunction(chromosome):
  f=func(chromosome)
  if f==0:
    fit=math.inf
  else:
    fit=1/f
  # fit=1/(1+f)
  return fit

class individual:
  def __init__(self, chromosome=[], fitness=-1):
    if fitness==-1:
      self.chromosome=chromosome
      self.fitness=fitnessFunction(self.chromosome)
    else:
      self.chromosome=chromosome
      self.fitness=fitness

def create_initial_population():
  population=[]
  for i in range(number_of_population):
    chromosome=[]
    for j in range(n):
      chromosome.append(uniform(-10,10))
    if sigma_modifier == 2:  #self_adaption
      for t in range(n):
        chromosome.append(sigma)
    population.append(individual(chromosome))
  return population

################################################## xover #############################################################

def single_xover(i1,i2):
  point=randint(1,n-1)
  x_temp=[]
  x_temp.append(i1.chromosome[point]*alpha + i2.chromosome[point]*(1-alpha))
  if sigma_modifier==2:  #self_adaption
    sigma_temp=[]
    sigma_temp.append(i1.chromosome[point+n]*alpha + i2.chromosome[point+n]*(1-alpha))
    chromosome1=i1.chromosome[:point]+x_temp[:]+i1.chromosome[point+1:point+n]+sigma_temp[:]+i1.chromosome[point+n+1:]
    chromosome2=i2.chromosome[:point]+x_temp[:]+i2.chromosome[point+1:point+n]+sigma_temp[:]+i2.chromosome[point+n+1:]
    return individual(chromosome1), individual(chromosome2)
  #oneFifth_Success_Rule
  chromosome1=i1.chromosome[:point]+x_temp[:]+i1.chromosome[point+1:]
  chromosome2=i2.chromosome[:point]+x_temp[:]+i2.chromosome[point+1:]
  return individual(chromosome1), individual(chromosome2)

def simple_xover(i1,i2):
  p=point=randint(1,n-1)
  x_temp=[]
  while(p<n):
    x_temp.append(i1.chromosome[p]*alpha + i2.chromosome[p]*(1-alpha))
    p=p+1
  if sigma_modifier==2:  #self_adaption
    p=point
    sigma_temp=[]
    while(p<n):
      sigma_temp.append(i1.chromosome[p+n]*alpha + i2.chromosome[p+n]*(1-alpha))
      p=p+1
    chromosome1=i1.chromosome[:point]+x_temp[:]+i1.chromosome[n:n+point]+sigma_temp[:]
    chromosome2=i2.chromosome[:point]+x_temp[:]+i2.chromosome[n:n+point]+sigma_temp[:]
    return individual(chromosome1), individual(chromosome2)
  #oneFifth_Success_Rule
  chromosome1=i1.chromosome[:point]+x_temp[:]
  chromosome2=i2.chromosome[:point]+x_temp[:]
  return individual(chromosome1), individual(chromosome2)

def whole_xover(i1,i2):
  x_temp=[]
  for p in range(n):
    x_temp.append(i1.chromosome[p]*alpha + i2.chromosome[p]*(1-alpha))
  if sigma_modifier==2:  #self_adaption
    p=0
    sigma_temp=[]
    while(p<n):
      sigma_temp.append(i1.chromosome[p+n]*alpha + i2.chromosome[p+n]*(1-alpha))
      p=p+1
    chromosome1=chromosome2=x_temp[:]+sigma_temp[:]
    return individual(chromosome1), individual(chromosome2)
  #oneFifth_Success_Rule
  chromosome1=chromosome2=x_temp[:]
  return individual(chromosome1), individual(chromosome2)

def LocalIntermediary_xover(i1,i2):
  chromosome=[]
  for p in range(n):
    chromosome.append( (i1.chromosome[p]+i2.chromosome[p])/2 )
  return individual(chromosome)

def GlobalIntermediary_xover(individuals):
  chromosome=[]
  for p in range(n):
    chromosome.append( (individuals[randint(0,len(individuals)-1)].chromosome[p]+individuals[randint(0,len(individuals)-1)].chromosome[p])/2 )
  return individual(chromosome)

def LocalDiscrete_xover(i1,i2):
  chromosome=[]
  for p in range(n):
    if r()<0.5:
      chromosome.append(i1.chromosome[p])
    else:
      chromosome.append(i2.chromosome[p])
  return individual(chromosome)

def GlobalDiscrete_xover(individuals):
  chromosome=[]
  for p in range(n):
    chromosome.append( (individuals[randint(0,len(individuals)-1)].chromosome[p]) )
  return individual(chromosome)

################################################## mutation #############################################################

def mutation(ind):
  if sigma_modifier==2:  #self_adaption
    xPrim=[]
    sigmaPrim=[]
    nn=n
    double_n=2*n
    while nn<double_n:
      sigmaPrim.append( ind.chromosome[nn]*math.exp((1/((2*n)**0.5))*gauss(0,1) + (1/((2*(n**0.5))**0.5))*gauss(0,1)) )
      nn=nn+1
    for i in range(n):
      total = ind.chromosome[i] + sigmaPrim[i]*gauss(0,1)
      if total>-10 and total<10:
        xPrim.append(total)
      else:
        xPrim.append(ind.chromosome[i])
    for s in sigmaPrim:
      xPrim.append(s)
    xPrim_fitness=fitnessFunction(xPrim)
    if xPrim_fitness > ind.fitness:
      return individual(xPrim,xPrim_fitness)
    return ind
  #oneFifth_Success_Rule
  xPrim=[]
  for x in ind.chromosome:
    total = x + gauss(0,sigma)
    if total>-10 and total<10:
      xPrim.append(total)
    else:
      xPrim.append(x)
    # xPrim.append(x)    
  xPrim_fitness=fitnessFunction(xPrim)
  if xPrim_fitness > ind.fitness:
    return individual(xPrim,xPrim_fitness) , 1
  return ind , 0

################################################## Selection methods #############################################################

def roulette_wheel(individuals, number):
  line=[(individuals[0], individuals[0].fitness)]
  for i in individuals[1:]:
    line.append((i, line[-1][1]+i.fitness))
  selected=[]
  for i in range(number):
    point=line[-1][1]*r()
    for j in line:
      if point<j[1]:
        selected.append(j[0])
        break
  return selected

# def sorted_by_fitness(e):
#   return e.fitness

# def sus(individuals, number):
#   individuals.sort(reverse=True, key=sorted_by_fitness)
#   line=[(individuals[0], individuals[0].fitness)]
#   for i in individuals[1:]:
#     line.append((i, line[-1][1]+i.fitness))
#   rand_Number=uniform(0,1/number)
#   selected=[]
#   for i in range(number):
#     point=line[-1][1]*rand_Number
#     for j in line:
#       if point<j[1]:
#         selected.append(j[0])
#         break
#     rand_Number=rand_Number + (1/number)
#   return selected

def sus(individuals, number=1):
  n=int(number)
  selected=[]
  sum_of_chances= individuals[0].fitness
  SS= [( individuals[0],individuals[0].fitness )]
  for indv in individuals[1:]:
      sum_of_chances += indv.fitness
      SS.append((indv,sum_of_chances) )
  step=sum_of_chances/n
  select=r()*step
  i=0
  while(n>0):
    if select < SS[i][1]:
      selected.append( SS[i][0] )
      select += step
      n=n-1
      i=i-1
    i= i+1
  return selected

def tournament(individuals, number):
  selected=[]
  for i in range(number):
    max_fitness=0
    updated_individuals=sample(individuals,tournament_k)
    for j in updated_individuals:
      if j.fitness > max_fitness:
        best=j
    selected.append(j)
  return selected

################################################## scaling #############################################################

def pawer_Law_scaling(individuals):
  new_individuals=[]
  for ind in individuals:
    new_individuals.append( individual(chromosome=ind.chromosome, fitness=(ind.fitness**pawer_Law_scaling_k) ) )
  return new_individuals

def panjere_bandi_scaling(individuals):
  min_fitness=math.inf
  new_individuals=[]
  for ind in individuals:
    if ind.fitness<min_fitness:
      min_fitness=ind.fitness
  for ind in individuals:
    new_individuals.append( individual(chromosome=ind.chromosome, fitness=ind.fitness-min_fitness ) )
  return new_individuals

def boltzmann_scaling(individuals):
  new_individuals=[]
  sum_e=0
  for ind in individuals:
    sum_e=sum_e+math.exp(ind.fitness/boltzmann_scaling_T)
  for ind in individuals:
    new_individuals.append( individual(chromosome=ind.chromosome, fitness= math.exp(ind.fitness/boltzmann_scaling_T)/sum_e ) )
  return new_individuals

def Sigma_Cutting_scaling(individuals):
  # 1<=c<=3
  sum_of_fitness=0
  new_individuals=[]
  for ind in individuals:
    sum_of_fitness=sum_of_fitness+ind.fitness
  avg_fitness=sum_of_fitness/len(individuals)
  for ind in individuals:
    new_individuals.append( individual(chromosome=ind.chromosome, fitness=ind.fitness-(avg_fitness - Sigma_Cutting_scaling_c*sigma) ) )
  return new_individuals

def Linear_scaling(individuals):
  sum_of_fitness=0
  max_fitness=0
  new_individuals=[]
  for ind in individuals:
    sum_of_fitness=sum_of_fitness+ind.fitness
    if ind.fitness>max_fitness:
      max_fitness=ind.fitness
  avg_fitness=sum_of_fitness/len(individuals)
  a,b=symbols('a b')
  eq1=Eq(sum_of_fitness*a + len(individuals)*b - sum_of_fitness , 0)
  eq2=Eq(a*max_fitness + b - Linear_scaling_c*avg_fitness , 0)
  sol_dict=solve((eq1,eq2), (a, b))
  # print(f'a = {sol_dict[a]}')
  # print(f'b = {sol_dict[b]}')
  for ind in individuals:
    fit=float(f'{sol_dict[a]}')*ind.fitness+float(f'{sol_dict[b]}')
    if fit>=0:
      new_individuals.append( individual(chromosome=ind.chromosome, fitness=fit) )
    else:
      new_individuals.append( individual(chromosome=ind.chromosome, fitness=0) )
  return new_individuals

################################################## ####### #############################################################

def create_generation(parents):
  children=[]
  number_of_mutation=0
  success=0
  for i in range(number_of_children):
    p1=randint(0,len(parents)-1)
    p2=randint(0,len(parents)-1)
    if r()<0.1:
      child1,child2=xover(parents[p1],parents[p2])
    else:
      child1,child2=parents[p1],parents[p2]
    if r()<0.9:
      if sigma_modifier==1:  #oneFifth_Success_Rule
        child1 , s = mutation(child1)
        number_of_mutation = number_of_mutation+1
        success = success+s
      elif sigma_modifier==2:  #self_adaption
        child1 = mutation(child1)
    if r()<0.9:
      if sigma_modifier==1:  #oneFifth_Success_Rule
        child2 , s = mutation(child2)
        number_of_mutation = number_of_mutation+1
        success = success+s
      elif sigma_modifier==2:  #self_adaption
        child2 = mutation(child2)
    children.append(child1)
    children.append(child2)
  if sigma_modifier==2:  #self_adaption
    return children
  #oneFifth_Success_Rule
  return children , number_of_mutation , success

################################################## ####### #############################################################

def oneFifth_Success_Rule(number_of_mutation, success):
  global sigma
  if success/number_of_mutation > 0.2:
    sigma = number_of_mutation/C
  elif success/number_of_mutation < 0.2:
    sigma = number_of_mutation*C


# main
average=[]
output_buffer=""
read_from_file()
if benchmark_function==1:
  func=f1
elif benchmark_function==2:
  func=f2
elif benchmark_function==3:
  func=f3
if xover_number==1:
  xover=single_xover
elif xover_number==2:
  xover=simple_xover
elif xover_number==3:
  xover=whole_xover
if selection_method_number==1:
  selection_method=roulette_wheel
elif selection_method_number==2:
  selection_method=sus
elif selection_method_number==3:
  selection_method=tournament
if scaling_number==1:
  scaling=pawer_Law_scaling
elif scaling_number==2:
  scaling=panjere_bandi_scaling
elif scaling_number==3:
  scaling=boltzmann_scaling
elif scaling_number==4:
  scaling=Sigma_Cutting_scaling
elif scaling_number==5:
  scaling=Linear_scaling

count=1
population=create_initial_population()
for i in range(number_of_generation):
  parents=selection_method(population, number_of_population)
  if sigma_modifier==1:  #oneFifth_Success_Rule
    children , number_of_mutation , success = create_generation(parents)
  elif sigma_modifier==2:  #self_adaption
    children = create_generation(parents)
  population=selection_method(parents+children, number_of_children)
  if scaling_number!=0: #scaling
    if count%(number_of_generation/number_of_scaling)==0:
      population=scaling(population)
  if sigma_modifier==1:  #oneFifth_Success_Rule
    oneFifth_Success_Rule(number_of_mutation, success)
  count=count+1
  avg=0
  for i in population:
    avg=avg+i.fitness
  # print(avg/len(population))
  average.append(avg/len(population))

best_in_population=0
for i in population:
  if best_in_population<i.fitness:
    best_in_population=i.fitness
    best=i


print(f2(best.chromosome))


output_buffer+="var chromosome = '"+str(best.chromosome)+"';\n"
output_buffer+="var best = "+str(best.fitness)+";\n"
output_buffer+="var valuesRy = ["
for i in range(len(average)):
  output_buffer+=str(average[i])+", "
# print(output_buffer)
output_buffer+="]\n"
output_buffer+="var propsRy = ["
for i in range(len(average)):
  output_buffer+='"'+str(i)+'", '
output_buffer+="]\n"
writing_file=open("files/variable.txt" ,'w')
writing_file.write(output_buffer)
writing_file.close()