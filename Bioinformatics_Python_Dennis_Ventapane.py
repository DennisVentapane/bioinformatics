#DETERMINAR PADROES

import sys
lines = sys.stdin.read().splitlines()
texto = lines[0]

"""Problema determinar padroes: retorna todas as posicoes de Genoma(texto) onde uma substring Padrao se (re)inicia"""

padrao = "ACACCA"
k = len(list(padrao))

posicoes = []

for i in range(len(texto)):
    if texto[i:i+k] == padrao:
   	 posicoes.append(i)

print posicoes



#CONTAR PADROES INEXATOS

def DistHam(p, q):
    d = 0
    for i in range(len(p)):
   	 if p[i] != q[i]:
   		 d += 1
    return d


texto = "TTTAGAGCCTTCAGAGG"
padrao = "GAGG"
d = 2

k = len(padrao)

ocorrencias = 0

for i in range((len(texto)-k)+1):
    pad = texto[i:i+k]
    if DistHam(pad, padrao) <= d:
   	 ocorrencias += 1


print ocorrencias




#DESEQUILIBRIO

import sys
lines = sys.stdin.read().splitlines()
genoma = lines[0]
 
desequilibrio = 0
 
for i in range(len(genoma)):
	 
	if genoma[i] == "C":
    	desequilibrio -= 1
    	 
	if genoma[i] == "G":
    	desequilibrio += 1
    	 
   #continue
 
 
   print (str(i), desequilibrio)



#DESCOBRIR PALAVRAS FREQUENTES

import sys
lines = sys.stdin.read().splitlines()

texto = lines[0]
k = 12
freq = [0]* (len(texto)-k)
PadraoFrequente = []



def ContarPadrao(texto,padrao):
	"""esta funcao retorna o numero de ocorrencias de padrao em texto"""
	
	l = len(padrao)
	ocorrencias = 0
	for c in range(len(texto)):
		if texto[c:c+l] == padrao:
			ocorrencias += 1
	return ocorrencias


for i in range(len(texto)-k):
    """ este for ira iterar sobre o modulo do texto - 12 e ira formar palavras de tamanho 12. Apos, ira tomar o numero de ocorrencias e colocar na lista freq"""

	palavra = texto[i:i+k]
	freq[i] = ContarPadrao(texto, palavra)
	

maxFreq = max(freq)


for j in range(len(texto)-k):
    """ este for ira iterar sobre o modulo do texto-12. Se o numero de ocorrencias for igual ao numero maximo de ocorrencias, a lista PadraoFreq sera acrescentada ao padrao mais frequente"""

	if freq[j] == maxFreq:
		PadraoFrequente.append(texto[j:j+k])

d = {}
for padrao in PadraoFrequente:
	""" como o dicionario tem chaves unicas podemos retornar d keys e ai teremos os padroes nao repetidos """
	if padrao not in d.keys():
		d[padrao] = maxFreq

print d.keys()
 


#REVERSO COMPLEMENTO

import sys
lines = sys.stdin.read().splitlines()

texto = lines[0]

complemento = ""

for i in range(len(texto)):
	if texto[i] == 'G':
		complemento += 'C'
	elif texto[i] == 'C':
		complemento += 'G'
	elif texto[i] == 'A':
		complemento += 'T'
	else:
		complemento += 'A'


reverso = ""

for i in range(len(complemento)):
	reverso = complemento[i] + reverso


print complemento
print reverso



#ALGORITMO FORÇA BRUTA PARA ENCONTRAR MOTIVOS (1.3 bitfun)

ltstr=["ATTTGGC","TGCCTTA","CGGTATC","GAAAATT"] 
k=3
d=1

def vizinhacad1(padrao): #gera vizinhaca do padrao com d=1
	nt=["A","C","T","G"]
	vizinhos=[]
	for i in range(len(padrao)):
		for j in nt:
			padrao2=list(padrao)
			if padrao2[i] != j:
				padrao2[i] = j
				vizinhos.append("".join(padrao2))
			elif padrao2[i] == j and "".join(padrao2) not in vizinhos:
				vizinhos.append("".join(padrao2))
	return vizinhos
	
def dham(padrao,padrao2): #calcula distancia de hamming entre 2 padroes
	disth=0
	for j in range(len(padrao)): 
		if padrao[j] != padrao2[j]:
			disth += 1
	return disth

listavizi=[] #aqui comeca a forca bruta de verdade
for i in range(len(ltstr)):
	for j in range(len(ltstr[i])-k+1):
		padrao=ltstr[i][j:j+k]
		vizi=vizinhacad1(padrao)
		for m in vizi:
			llmatch = []
			instr=0
			for l in range(len(ltstr)):
				llmatch.append(0)
				for n in range(len(ltstr[l])-k+1):
					padrao2 = ltstr[l][n:n+k]
					if dham(m,padrao2) <= d:
						llmatch[l] = 1
			for o in llmatch:
				if o == 1:
					instr+=1
			if instr == len(llmatch) and m not in listavizi:
				listavizi.append(m)

print listavizi


#ALGORITMO PARA PROCURAR TODOS OS PADRÕES DENTRO DE UMA SEQUENCIA (INPUTS SAO: SEQUENCIA, THRESHOLD DE MATCHES PARA EXIBIR O RESULTADO, TAMANHO DO K-MER (SE 0, PROCURA TODAS POSSIVEIS), DISTANCIA HAMMING)

import sys
lines = sys.stdin.read().splitlines()

def dham(padrao,padrao2):
	disth=0
	for j in range(len(padrao)): 
		if padrao[j] != padrao2[j]:
			disth += 1
	return disth

def fullsearch(seq,t,k,d):
	print "O tamanho da sequencia eh",len(seq),", o threshold eh",t,"e a distancia de hamming eh",d
	library=[]
	if k == 0:
		k = range(1,len(seq)+1)
		print "Procurando todas as combinacoes possiveis"
	else:
		k = [int(k)]
		print "Procurando",k,"-mers"
	for l in k:
		for c in range(len(seq)-l+1):
			occ = 0
			padrao = seq[c:c+l]
			if padrao not in library:
				for i in range(len(seq)-l+1):
					padrao2 = seq[i:i+l]
					if d > 0:
						if dham(padrao,padrao2) <= d:
							occ += 1	
					elif padrao == padrao2:
						occ += 1
				if occ >= t:
					library.append(padrao)
					print padrao, occ

fullsearch(lines[0],int(lines[1]),int(lines[2]),int(lines[3])









#GERAR COMBINACOES DE TAMANHO TAM COM ATE F MISMATCHES

"""
res = lista de posicoes
p = posicao
f = erros
tam = tamanho

"""

def gerar_combi(p, f, tam, res):

    #print ">> p=" + str(p) + ", f="+ str(f)+ ", tam="+ str(tam)+ ", res=" + str(res)
    #condicao de parada
    if f == 0:
   	 print res
   	 return
    
    for i in range(p, tam- f+ 1):
   	 res2 = res[0:len(res)]
   	 res2.append(i)
   	 #print str(len(res2)) + " " + str(i)
   	 gerar_combi(i+1, f-1, tam, res2)

gerar_combi(0, 3, 3, [])



#STRING MEDIANA 
    	
"""String mediana"""

import itertools
#gera combinacoes

def distHam(p,q):
    dH = 0
    for i in range(len(p)):
   	 if p[i] != q[i]:
   		 dH += 1
    return dH



def GerarPadroes(k):
    padroes_tuplas = list(itertools.product("AGCT", repeat = k))
    padroes_texto = []
    for padrao_tupla in padroes_tuplas:
   	 padroes_texto.append("".join(padrao_tupla))
    return padroes_texto


def StringMediana(Dna,k):
    possiveis = GerarPadroes(k)
    dist = k
    motivos= []
    distancias_hamming = []
    best = " "
    mediana = " "
    bestscore =k * len(Dna)

    for padrao in possiveis:
   	 for i in range(len(Dna)):
   		 for j in range(len(Dna[i])-k+1):
   			 kmer = Dna[i:i+k]
   			 dH = distHam(kmer, padrao)
   			 if (dH < dist):
   				 best = kmer
    
   		 distancias_hamming.append(dH)
   	 
   		 motivos.append(best)

   	 score = 0
   	 for dH in distancias_hamming:
   		 score += dH    

   	 if score < bestscore:
   		 bestscore = score
   		 
   		 

    for best in motivos:
   	 if dist > distHam(padrao, Dna):
   		 dist =distHam(padrao, Dna)
   		 mediana = padrao

    return mediana



Dna = "AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG",
"GCTGAGCACCGG" , "AGTTCGGGACAG"

k = 3
print(StringMediana(Dna, k))

	




#STRING MEDIANA 2 - ALGORITMO PARA ENCONTRAR O K-MER COM MENOR SCORE ENTRE VARIAS LINHAS DE DNA (POR HECTOR)

import itertools
import sys

lines = sys.stdin.read().splitlines()

k = int(lines[0])	#primeira linha do input eh o comprimento dos padroes que se quer encontrar
DNAS = []	#aqui cria um conjunto vazio onde entrarao as strings de dna
for i in range(1,len(lines)):	#este loop acrescenta cada linha de dna na lista DNAS
	DNAS.append(lines[i])


def padroes (k): #este algoritmo utiliza da ferramenta itertools para gerar uma lista com todos os padroes possiveis de tamanho k
	padroes_tuples = list(itertools.product("ACGT", repeat=k))
	padroes_strings = []
	for i in padroes_tuples:
		padroes_strings.append("".join(i))
	return padroes_strings
	
def dham(padrao,padrao2):	#bom e velho algoritmo que calcula a distancia hamming entre 2 padroes
	disth=0
	for j in range(len(padrao)): 
		if padrao[j] != padrao2[j]:
			disth += 1
	return disth
	
def best(DNA, padrao): #este algoritmo seleciona o trecho do dna com menor distancia hamming em relacao ao padrao e retorna ele
	distancia = len(padrao)+1
	best = ""
	for i in range (len(DNA)-len(padrao)+1):
		padrao2 = DNA[i:i+k]
		disth = dham(padrao, padrao2)
		if disth < distancia:
			distancia = disth
			best = padrao2
	return best

def best_matrix(DNAS, k): #este algoritmo monta uma matriz com todos os padroes e seus equivalente com menor distancia hamming em cada linha de dna (primeiro elemento eh o padrao, e os seguintes sao seus equivalente em cada fita)
	optimal_all = []
	for j in padroes (k):
		optimal_list = []
		padrao = j
		optimal_list.append(padrao)
		for i in range(len(DNAS)):
			optimal_list.append(best(DNAS[i], padrao))
		optimal_all.append(optimal_list)
	return optimal_all
	
def lowest_score(matrix): #esse calcula o score de cada padrao da matriz gerada e retorna somente o(s) padrao(oes) com menor score e seu(s) respectivo(s) score(s)
	lowest = ["", ]
	lowest_scr = len(matrix[0][0])*len(matrix[0])
	for i in range(len(matrix)):
		score = 0
		for j in range(1,len(matrix[i])):
			score += dham(matrix[i][0], matrix[i][j])
		if score < lowest_scr:
			lowest_scr = score
			lowest[0] = matrix[i][0]
	for k in range(len(matrix)): #esse aqui checa se ha algum outro padrao com mesmo score, para retornar
		score = 0
		for l in range(1,len(matrix[k])):
			score += dham(matrix[k][0], matrix[k][l])
		if score == lowest_scr and matrix[k][0] not in lowest:
			lowest.append(matrix[k][0])
	return lowest, lowest_scr
	
print lowest_score(best_matrix(DNAS, k))


#K-MER DE MAIOR PROBABILIDADE NUMA STRING, DADO UMA TABELA DE PERFIL

import sys

lines = sys.stdin.read().splitlines()

string = lines[0]	#aqui a primeira linha do input eh atribuida a uma variavel string
k = int(lines[1])	#a segunda linha eh o comprimento das sequencias a serem analisadas
prA = lines[2].split()	#aqui a terceira linha, que representa as probabilidades do nucleotideo A, em cada caracter do k-mer, eh dividido e adicionado a variavel prA, como uma lista onde cada indice eh uma posicao do k-mer
prC = lines[3].split()	#o mesmo aqui para o nucleotideo C
prG = lines[4].split() #e aqui para o G
prT = lines[5].split()	#e T


def prob_kmer(string, k, prA, prC, prG, prT): #aqui eh definida a funcao de probabilidade do kmer... os inputs sao a string, o k e o perfil 'destrinchado' em 4 listas, de probabilidades
	prob_save = 0	#aqui cria uma variavel vazia para futuramente salvar a maior probabilidade
	seq_save = ""	#aqui cria um string vazio para futuramente salvar a sequencia com maior probabilidade
	for i in range(len(string)-k):	#esse loop cria a janela deslizante que passa de 1 em 1 nucleotideo na string de input
		seq = string[i:i+k]	#aqui gera progressivamente as sequencias a partir da janela deslizante
		prstring = []	#aqui eh uma lista vazia onde o indice diz a probabilidade de cada posicao da sequencia analisada, baseado no perfil do input
		for j in range(len(seq)):	#aqui eh analisado de caracter a caracter da sequencia criada
			if seq[j] == "A":	#daqui pra baixo verifica o caracter que ta na posicao J e, baseado no caracter identificado, anexa a probabilidade daquele caracter, naquela posicao, a prstring
				prstring.append(prA[j])
			elif seq[j] == "C":
				prstring.append(prC[j])
			elif seq[j] == "G":
				prstring.append(prG[j])
			elif seq[j] == "T":
				prstring.append(prT[j])
		multprstring = float(prstring[0]) #aqui cria uma variavel cujo unico elemento eh o primeiro valor de prstring, ou seja, a probabilidade do primeiro caracter da sequencia gerada pela janela deslizante
		for l in range(1,len(prstring)):	#nesse loop sao multiplicados todos os valores contidos em prstring, a fin de dar a probabilidade final da sequencia gerada
			multprstring *= float(prstring[l])
		if multprstring > prob_save:	#entao, se a probabilidade for maior que a salva anteriormente, salva ela no lugar da anterior e salva a sequencia tambem
			prob_save = multprstring
			seq_save = seq
	return seq_save, prob_save	#entao, a sequencia com maior probabilidade eh retornada junto com a probabilidade
	

print prob_kmer(string, k, prA, prC, prG, prT)



#RASCUNHO DUDA (NAO ESTA COMPLETO)


"""GreedyBuscaMotivos(Dna,k)"""

def contagem(Dna, k):
    
    contagem = [[[0]*k]*4]
    t = len(Dna)
    #Dnai = Dna[i][j]
    contagem = [[[0]*k]*4]
    for i in range(t):
   	 for j in range(t):
   		 if Dna[i][j] == "A":
   			 contagem[0][i] += 1
   		 elif Dna[i][j] == "C":
   			 contagem[1][i] += 1
   		 elif Dna[i][j] == "G":
   			 contagem[2][i] += 1
   		 else:
   			 contagem[3][i] += 1

    return contagem


def
   			 


def GreedyBuscaMotivos(Dna,k):

    t = len(Dna)
    MelhoresMotivos = []
    
    contagem = [[0*k]*4]
    perfil = [[0*k]*4]


    for i in range(len(Dna)):
   	 for j in range(len(Dna[i])-k+1):
   		 

   			 
    
    motivos = []
    if score(motivos) < score(MelhoresMotivos):
   	 MelhoresMotivos = motivos

    return MelhoresMotivos


Dna = "GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"
 
t = len(Dna)
k = 3
print (Dna[2][1])
print ([[[0]*k]*4])
print (t)

#ALGORITMO GREEDY PARA A BUSCA DA STRING MEDIANA

#Esse codigo reaproveita muitos algoritmos de codigos anteriores que fiz

import sys

lines = sys.stdin.read().splitlines()

k = int(lines[0])
DNA = []
for i in range(1,len(lines)):
	DNA.append(lines[i])
	
	
def perfil_laplace(kmer):	#Esse algoritmo monta uma matriz de perfil para um unico dado kmer, corrigido por laplace. Como as probabilidades sao fixas para esse caso de um unico kmer, o calculo eh evitado ao se utilizar valores fixos pre-calculados.
	perfil = []
	prA = []
	prC = []
	prG = []
	prT = []
	for i in range(len(kmer)):	#Por exemplo: sao 4 nucleotideos diferentes, e em cada posicao de um kmer soh havera um unico nucleotideo, o que da 100% de chance pra esse e 0% pros demais. No entanto, com a correcao de laplace, as chances mudam para 40% do que esta na posicao k, e 20% para os que nao estao
		if kmer[i] == "A":
			prA.append(0.40)
			prC.append(0.20)
			prG.append(0.20)
			prT.append(0.20)
		elif kmer[i] == "C":
			prA.append(0.20)
			prC.append(0.40)
			prG.append(0.20)
			prT.append(0.20)
		elif kmer[i] == "G":
			prA.append(0.20)
			prC.append(0.20)
			prG.append(0.40)
			prT.append(0.20)
		elif kmer[i] == "T":
			prA.append(0.20)
			prC.append(0.20)
			prG.append(0.20)
			prT.append(0.40)
	perfil.append(prA)
	perfil.append(prC)
	perfil.append(prG)
	perfil.append(prT)
	return perfil
	
def prob_kmer(string, k, perfil):	#Este eh o mesmo algoritmo do ultimo codigo. A diferenca eh q ele n retorna a probabilidade, mas sim somente a string da sequencia
	prob_save = 0
	seq_save = ""
	for i in range(len(string)-k):
		seq = string[i:i+k]
		prstring = []
		for j in range(len(seq)):
			if seq[j] == "A":
				prstring.append(perfil[0][j])
			elif seq[j] == "C":
				prstring.append(perfil[1][j])
			elif seq[j] == "G":
				prstring.append(perfil[2][j])
			elif seq[j] == "T":
				prstring.append(perfil[3][j])
		multprstring = float(prstring[0])
		for l in range(1,len(prstring)):
			multprstring *= float(prstring[l])
		if multprstring > prob_save:
			prob_save = multprstring
			seq_save = seq
	return seq_save
	
def dham(padrao,padrao2):	#Novamente o bom e velho algoritmo que calcula a distancia hamming entre 2 padroes
	disth=0
	for j in range(len(padrao)): 
		if padrao[j] != padrao2[j]:
			disth += 1
	return disth
	
def merge_perfis(perfil1, perfil2):	#Esse algoritmo junta quaisquer duas tabelas de perfil de sequencias de mesmo tamanho k, ajustando as probabilidades de cada nucleotideo
	for i in range(len(perfil1)):
		for j in range(len(perfil1[i])):
			perfil1[i][j] = (perfil1[i][j]+perfil2[i][j])/2
	return perfil1

def lowest_score(matrix): #Usado no codigo da string mediana, que fiz anteriormente. Calcula o score de cada padrao de uma matriz e retorna somente o(s) padrao(oes) com menor score e seu(s) respectivo(s) score(s)
	lowest = ["", ]
	lowest_scr = len(matrix[0][0])*len(matrix[0])
	for i in range(len(matrix)):
		score = 0
		for j in range(1,len(matrix[i])):
			score += dham(matrix[i][0], matrix[i][j])
		if score < lowest_scr:
			lowest_scr = score
			lowest[0] = matrix[i][0]
	for k in range(len(matrix)):
		score = 0
		for l in range(1,len(matrix[k])):
			score += dham(matrix[k][0], matrix[k][l])
		if score == lowest_scr and matrix[k][0] not in lowest:
			lowest.append(matrix[k][0])
	return lowest, lowest_scr
	
def greedy_motivos_matrix(DNA, k):	#Esse aqui eh um algoritmo greedy que gera uma matriz onde o primeiro elemento de cada elemento eh uma sequencia gerada pela janela deslizante na primeira linha de dna, e os demais os kmers mais provaveis em relacao a tabela de perfil gerada do primeiro
	best_of_bests = []
	for j in range(len(DNA[0])-k):
		sequencia = DNA[0][j:j+k]
		perfil_sequencia = perfil_laplace(sequencia)
		best_kmers = [sequencia, ]
		for l in range(1,len(DNA)):
			best = prob_kmer(DNA[l], k, perfil_sequencia)
			best_kmers.append(best)
			perfil_sequencia = merge_perfis(perfil_sequencia, perfil_laplace(best))
		best_of_bests.append(best_kmers)
	return best_of_bests

print lowest_score(greedy_motivos_matrix(DNA, k))



#ROLAR OS DADOS

Motivos = []
for sequencia in Dna:
    posicao = rd.randint(0, len(sequencia)-k)
    #posicao eh o inicio da kmer

    Motivos.append(sequencia[posicao:posicao+k])

print Motivos


#segundo passo: matriz de contagem


def CreateContagem(lista):
    contagem = []

    for i in range(k):
   	 A = 0
   	 C = 0
   	 T = 0
   	 G = 0

   	 for j in range(t):
   		 x = lista[j][i]
   		 if x == "A":
   			 A += 1
   		 elif x == "C":
   			 C += 1
   		 elif x == "G":
   			 G += 1
   		 else:
   			 T += 1


   	 contagem.append([A, C, G, T])

    return contagem    


def ConvertProfile(matriz):
    
    perfil = []
    
    for i in range(len(matriz)):
   	 sublista = []
   	 for j in range(len(matriz[i])):
   		 soma = sum(matriz[i])
   		 freq = float(matriz[i][j])/float(soma)
   		 sublista.append(freq)
   	 
   	 perfil.append(sublista)

    return perfil



def pseudocount(matriz):

    pseudo_count = []

    for i in range(len(matriz)):
   	 sublista = []
   	 for j in range(len(matriz[i])):
   		 laplace = matriz[i][j] + 1
   		 sublista.append(laplace)

   	 pseudo_count.append(sublista)

    return pseudo_count


def give_motivos(perfil, Dna):
    #queremos uma lista de motivos como output

    lista_motivos = []    
    for seq in Dna:
   	 motivo = scan_motivo(perfil, seq)
   	 lista_motivos.append(motivo)
    
    return lista_motivos


def scan_motivo(perfil, seq):

    max_prob = 0
    for i in range(len(seq)-k):
   	 kmer = seq[i:i+k]

   	 prob = 1
   	 
   	 for j in range(len(kmer)):
   		 letra = kmer[j]
   		 if letra == "A":
   			 pos = 0
   		 elif letra == "C":
   			 pos = 1
   		 elif letra == "G":
   			 pos = 2
   		 else:
   			 pos = 3

   		 freq = perfil[j][pos]

   		 prob = freq * prob

   	 if prob > max_prob:
   		 max_prob = prob
   		 max_kmer = kmer    
    
    return max_kmer


def score_motivos(motifs, perfil):

    consensus = get_consensus(perfil)
    
    score = 0

    for i in range(len(motifs)):
   	 for j in range(len(motifs[i])):
   		 if motifs[i][j] != consensus[j]:
   			 score += 1
    return score
   	 


def get_consensus(matriz):
    consensus = ""
    
    for i in range(len(matriz)):
   	 
   	 max_prob = 0
   	 
   	 for j in range(len(matriz[i])):
   		 if matriz[i][j] > max_prob:
   			 max_prob = matriz[i][j]
   			 max_pos = j



   	 if max_pos == 0:
   		 consensus += "A"
   	 elif max_pos == 1:
   		 consensus += "C"
   	 elif max_pos == 2:
   		 consensus += "G"
   	 else:
   		 consensus += "T"

    return consensus

   			 
melhor_score = 999999999999999
lista_motivos = Motivos
   		 
i= 0
while 1 == 1:

    i= i+ 1
    #print "- number "+ str(i)+ " -"

    contagem = CreateContagem(lista_motivos)
    #print contagem
    pseudo = pseudocount(contagem)
    #print pseudo
    perfil = ConvertProfile(pseudo)
    #print perfil
    lista_motivos = give_motivos(perfil, Dna)
    #print lista_motivos
    score = score_motivos(lista_motivos, perfil)
    #print score

    if score < melhor_score:
   	 #print score, lista_motivos
   	 melhor_score = score
   	 melhor_motivo = lista_motivos
   	 melhor_perfil = perfil

    else:
   	 print "melhor score depois de "+ str(i)+ " iteracoes"
   	 print melhor_score, melhor_motivo
   	 print get_consensus(melhor_perfil)
   	 break



def give_motivos(perfil, Dna):
    #queremos uma lista de motivos como output

    lista_motivos = []    
    for seq in Dna:
   	 motivo = scan_motivo(perfil, seq)
   	 lista_motivos.append(motivo)
    
    return lista_motivos



#GIBBS SAMPLING

import random as rd

Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]  
k = 8
t = 5


#primeiro passo: selecionar motivos aleatorios

Motivos = []
for sequencia in Dna:
    posicao = rd.randint(0, len(sequencia)-k)
    #posicao eh o inicio da kmer

    Motivos.append(sequencia[posicao:posicao+k])

#print Motivos


#segundo passo: matriz de contagem


def CreateContagem(lista):
    contagem = []

    for i in range(k):
   	 A = 0
   	 C = 0
   	 T = 0
   	 G = 0

   	 for j in range(t):
   		 x = lista[j][i]
   		 if x == "A":
   			 A += 1
   		 elif x == "C":
   			 C += 1
   		 elif x == "G":
   			 G += 1
   		 else:
   			 T += 1


   	 contagem.append([A, C, G, T])

    return contagem    


def ConvertProfile(matriz):
    
    perfil = []
    
    for i in range(len(matriz)):
   	 sublista = []
   	 for j in range(len(matriz[i])):
   		 soma = sum(matriz[i])
   		 freq = float(matriz[i][j])/float(soma)
   		 sublista.append(freq)
   	 
   	 perfil.append(sublista)

    return perfil



def pseudocount(matriz):

    pseudo_count = []

    for i in range(len(matriz)):
   	 sublista = []
   	 for j in range(len(matriz[i])):
   		 laplace = matriz[i][j] + 1
   		 sublista.append(laplace)

   	 pseudo_count.append(sublista)

    return pseudo_count


def give_motivos(perfil, Dna):
    #queremos uma lista de motivos como output

    lista_motivos = []    
    for seq in Dna:
   	 motivo = scan_motivo(perfil, seq)
   	 lista_motivos.append(motivo)
    
    return lista_motivos

def subst_motivos(perfil, Dna, seq, motivos):
    #queremos uma lista de motivos como output

    lista_motivos = motivos[0:len(motivos)]    
    novo_motivo = scan_motivo(perfil, Dna[seq])
    lista_motivos[seq] = novo_motivo
   	 
    return lista_motivos

def scan_motivo(perfil, seq):

    max_prob = 0
    for i in range(len(seq)-k):
   	 kmer = seq[i:i+k]

   	 prob = 1
   	 
   	 for j in range(len(kmer)):
   		 letra = kmer[j]
   		 if letra == "A":
   			 pos = 0
   		 elif letra == "C":
   			 pos = 1
   		 elif letra == "G":
   			 pos = 2
   		 else:
   			 pos = 3

   		 freq = perfil[j][pos]

   		 prob = freq * prob

   	 if prob > max_prob:
   		 max_prob = prob
   		 max_kmer = kmer    
    
    return max_kmer


def score_motivos(motifs, perfil):

    consensus = get_consensus(perfil)
    
    score = 0

    for i in range(len(motifs)):
   	 for j in range(len(motifs[seq])):
   		 if motifs[i][j] != consensus[j]:
   			 score += 1
    return score
   	 


def get_consensus(matriz):
    consensus = ""
    
    for i in range(len(matriz)):
   	 
   	 max_prob = 0
   	 
   	 for j in range(len(matriz[i])):
   		 if matriz[i][j] > max_prob:
   			 max_prob = matriz[i][j]
   			 max_pos = j



   	 if max_pos == 0:
   		 consensus += "A"
   	 elif max_pos == 1:
   		 consensus += "C"
   	 elif max_pos == 2:
   		 consensus += "G"
   	 else:
   		 consensus += "T"

    return consensus

   			 
score_anterior = 999999999999999
lista_motivos = Motivos
   		 
i= 0
N = 10000
for i in range(N):
    i= i+ 1
    #print "- number "+ str(i)+ " -"

    contagem = CreateContagem(lista_motivos)
    #print contagem
    pseudo = pseudocount(contagem)
    #print pseudo
    perfil = ConvertProfile(pseudo)
    #print perfil
    seq = rd.randint(0, len(Dna)- 1)
    novo_motivos = subst_motivos(perfil, Dna, seq, lista_motivos)
    #print lista_motivos
    score = score_motivos(novo_motivos, perfil)
    #print score

    if score < score_anterior:
   	 lista_motivos= novo_motivos

    score_anterior= score

print lista_motivos
print score



#funcao para subamostragem enviesado pela lista de pesos fornecidos

import random

def randWeight(pesos):
    #converter pesos numa distribuicao de probabilidades
    soma = 0
    for i in range(len(pesos)):
   	 soma += pesos[i]
    pdist = []
    for i in range(len(pesos)):
   	 pdist.append(pesos[i]/soma)

    
    #criar uma distribuicao cumulativa da distribuicao cumulativa
    pcumu = []
    scumu = 0.0
    for i in range(len(pesos)):
   	 scumu += pdist[i]
   	 pcumu.append(scumu)

    #procurar uniform random na distribuicao cumulativa
    x = random.random()
    for i in range(len(pesos)):
   	 if x <= pcumu[i]:
   		 return i



pesos = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5]
faces = [1, 2, 3, 4, 5, 6]
for i in range(1000):
    print faces[randWeight(pesos)]



#CAMINHOS

"""caminhos"""

n = 16 + 1
m = 12 + 1

matriz = []

for i in range(n):
    matriz.append([0]*m)
    

#encher a primeira coluna com valores de 1:
for i in range(n):
    matriz[i][0] = 1

#encher a primeira linha com valores de 1:
for j in range(m):
    matriz[0][j] = 1

for i in range(1,n):
    for j in range(1,m):
   	 matriz[i][j] = matriz[i][j-1] + matriz[i-1][j]

print matriz
print "*---*"
print matriz[n-1][m-1]



#ALGORITIMO ALINHAMENTO EXAUSTIVO (RECURSIVO)
"""aula 16 de maio"""

from sys import maxint

maximum_int = maxint
minimum_int = -maxint - 2


def AlinhamentoExaustivo(v,w,i,j,peso):

    if i == len(v) and j == len(w):
   	 return peso

    x = y = z = minimum_int

    if i < len(v):
   	 x = AlinhamentoExaustivo(v,w,i+1, j, peso)
    if j < len(w):
   	 y = AlinhamentoExaustivo(v,w,i,j+1,peso)
    if i < len(v) and j < len(w):
   	 if v[i] == w[j]:
   		 peso = peso + 1
                    	#print peso
   	 z = AlinhamentoExaustivo(v,w,i+1,j+1,peso)
    #print x,y,z
    return max(x,y,z)



v = "ATGTTATA"
w = "ATCGTCC"



#i,j = 0,0 pois estamos na fonte
#peso = 0 no inicio
print AlinhamentoExaustivo(v,w,0,0,0)



#PROGRAMACAO DINAMICA
"""programacao dinamica 18 Maio"""


def ProgDin(v,w):

    
    n = len(v)
    m = len(w)

    cam = []
    direcao = []
    
    for i in range(n+1):
   	 sublista = [0] * (m+1)
   	 subdir = [""] * (m+1)
   	 cam.append(sublista)
   	 direcao.append(subdir)
   	 
   	 if i > 0:
   		 cam[i][0] = cam[i-1][0] + 0
   		 direcao[i][0] = "baixo"
   		 
    
    for j in range(1,m+1):
   	 
   	 cam[0][j] = cam[0][j-1] + 0
   	 direcao[0][j] = "lado"



    for i in range(1,n+1):
   	 
   	 for j in range(1,m+1):

   		 lado = cam[i][j-1] + 0

   		 baixo = cam[i-1][j] + 0

   		 x = v[i-1]
   		 y = w[j-1]

   		 if x == y:

   			 diag = cam[i-1][j-1] + 1
   		 else:
   			 diag = cam[i-1][j-1] + 0


   		 melhor_caminho = max(lado,baixo,diag)

   		 if melhor_caminho == lado:
   			 direcao[i][j] = "lado"

   		 elif melhor_caminho == baixo:
   			 direcao[i][j] = "baixo"

   		 else:
   			 direcao[i][j] = "diagonal"

   		 
   		 cam[i][j] = melhor_caminho

    print direcao

    
    i = n
    j = m

    v_alin = ""
    w_alin = ""

    while i != 0 or j != 0:


   	 if direcao[i][j] == "baixo":
   		 
   		 v_alin = v[i-1] + v_alin
   		 w_alin = "-" + w_alin
   		 i = i - 1

   	 elif direcao[i][j] == "lado":
   		 
   		 w_alin = w[j-1] + w_alin
   		 v_alin = "-" + v_alin
   		 j = j -1

   	 else:
   		 v_alin = v[i-1] + v_alin
   		 w_alin = w[j-1] + w_alin
   		 i = i - 1
   		 j = j - 1
   		 

    print v_alin
    print w_alin


    return cam[n-1][m-1]
   		 


v = "ATGTTATA"
w = "ATCGTCC"

print ProgDin(v,w)


#ALINHAMENTO GLOBAL 


mantes = []
with open("/home/aluno/Documentos/BLOSUM62.txt") as inputfile:
    for line in inputfile:
   	 mantes.append(line.split())


m = []
for i in  range(len(mantes)):
    if mantes[i][0] != "#":
   	 m.append(mantes[i])

v = m[0]
w = []
len_v = len(v)

dicionario = {}

j = 0
for i in range(len_v):
    dicionario[v[i]] = j
    j += 1



Q = dicionario["Q"]
R = dicionario["R"]

print Q
print R

print m[Q][R]
print m[5][1]




aa1 = "PLEASANTLY"
aa2 = "MEANLY"


"""
def AlinhamentoGlobal(seq1, seq2, dic, indel):

    dic = dicionario

    sequence1 = ""
    sequence1.append(seq1.split())

    sequence2 = ""
    sequence2.append(seq2.split())

    

    return



print AlinhamentoGlobal(aa1,aa2, dicionario, 5)

"""

#GAPS AFINOS

"""Gaps afinos

Objetivo: construir um alinhamento entre duas strings, com gap penalidades afinas

Input: 2 strings, uma matriz score, numeros sigma e epsilon
Output: o alinhamento com a pontuacao maxima entre essas strings, definido pela matriz score e pelas gap penalidades sigma e epsilon"""



def matriz_score(string1, string2):

    n = len(string1)
    m = len(string2)
    sigma = -1
    epsilon = -0.1

    diag = []
    for i in range(n+1):
   	 sublista = [0] * (m+1)
   	 diag.append(sublista)
   	 diag[i][0] = sigma + (epsilon * i)
    
    for j in range(m+1):
   	 diag[0][j] = sigma + (epsilon * j)
    

    horz = []
    for i in range(n+1):
   	 sublista = [0] * (m+1)
   	 horz.append(sublista)
   	 horz[i][0] = 2 * sigma + (epsilon * i)
   	 
    for j in range(m+1):   	 
   	 horz[0][j] = sigma + (epsilon * j)


    vert = []
    for i in range(n+1):
   	 sublista = [0] * (m+1)
   	 vert.append(sublista)
   	 vert[i][0] = sigma + (epsilon * i)

    for j in range(m+1):
   	 vert[0][j] = 2 * sigma + (epsilon * j)


    diag[0][0] = 0
    horz[0][0] = sigma
    vert[0][0] = sigma

    for i in range(1, n+1):
   	 for j in range(1, m+1):
   		 #passo horizontal   		 
   		 hh = horz[i][j-1] + epsilon
   		 
   		 #passo vertical
   		 vv = vert[i-1][j] + epsilon

   		 #passo diagonal
   		 var1 = string1[i-1]    		 
   		 var2 = string2[j-1]
   		 if var1 == var2:   	 
   			 dd = diag[i-1][j-1] + 1
   		 else:
   			 dd = diag[i-1][j-1] - 1
 
   		 diag[i][j] = max(hh, dd, vv)
   		 horz[i][j] = max(hh, diag[i][j] + sigma)
   		 vert[i][j] = max (vv, diag[i][j] + sigma)
   	 
    
    return max(diag[i][j], horz[i][j], vert[i][j])


string1 = "AGCTAAGC"
string2 = "GCTAAA"

print matriz_score(string1,string2)
    




