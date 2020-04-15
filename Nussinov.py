# Algoritmo de Nussinov para RNA folding

import numpy as np

# Lendo sequencia
arq = open('/home/aline/Documentos/Teste.txt', 'r')
seq = arq.read()
seq = seq.strip()
seq = list(seq)
print('Sequencia: {}'.format(seq))
arq.close()

# Criando matriz com zeros
def Cria_Matriz(tamanho):
    return np.zeros([tamanho, tamanho], dtype=int)

# Comparando para saber se as bases pareiam ou nao
def Pareia(y, x):

    bases = (y, x)
    if bases in (('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')):
        return 1
    return 0

# Calculando campos para preenchimento da matriz
def Calculo_matriz(M, i, j):

    # Condicoes de pareamento
    if (j < i) or (i == j) or (j == i + 1) or (j == i - 1):
        return 0
    else:
        k = i + 1
        list_k = []
        while (k < j):
            list_k.append(M[i, k] + M[k + 1, j])
            k = k + 1
        x0 = M[i + 1, j]
        x1 = M[i, j - 1]
        x2 = M[i + 1, j - 1] + (Pareia(seq[i], seq[j]))
        x3 = max(list_k)
    return max(x0, x1, x2, x3)

# Preenchendo a matriz
def Preenche_matriz(M, tam):
    for ck in range(2, tam):
        for ci in range(0, tam-ck):
           M[ci][ci+ck] = Calculo_matriz(M, ci, ci+ck)
    return M

#Tracebacking

#criando as pilhas
def empilha(M):
    pilha1 = []
    pilha2 = []
    i = 0
    j = len(seq) - 1

    #Inicializacao
    pilha1.append([i, j])

    while i < j:
        try:
            i, j = pilha1.pop()
        except IndexError:
            break

        if M[i+1, j] == M[i, j]:
            pilha1.append([i+1, j])

        elif M[i, j-1] == M[i, j]:
            pilha1.append([i, j-1])

        elif M[i+1, j-1] + Pareia(seq[i], seq[j]) == M[i, j]:
            pilha2.append([i, j])
            pilha1.append([i+1, j-1])

        else:
            k = i + 1
            while k < j:
                if M[i, k] + M[k + 1, j] == M[i, j]:
                    pilha1.append([k + 1, j])
                    pilha1.append([i, k])
                    break
                else:
                    k = k + 1

    return pilha2

# Chamando os metodos
Matriz = Cria_Matriz(len(seq))
Matriz = Preenche_matriz(Matriz, len(seq))
print(Matriz)
print("Resultado: {}".format(empilha(Matriz)))

########################Exemplo########################

'''
Entrada:

	Sequencia: ['G', 'G', 'G', 'A', 'A', 'A', 'U', 'C', 'C']

Saída:
	[[0 0 0 0 0 0 1 2 3]
	 [0 0 0 0 0 0 1 2 3]
	 [0 0 0 0 0 0 1 2 2]
	 [0 0 0 0 0 0 1 1 1]
	 [0 0 0 0 0 0 1 1 1]
	 [0 0 0 0 0 0 0 0 0]
	 [0 0 0 0 0 0 0 0 0]
	 [0 0 0 0 0 0 0 0 0]
	 [0 0 0 0 0 0 0 0 0]]

	Pares de bases: [[1, 8], [2, 7], [4, 6]]
'''
