#initialization
import matplotlib.pyplot as plt
import numpy as np

# importing Qiskit
from qiskit.tools.visualization import plot_bloch_multivector
from qiskit.quantum_info import Statevector
from qiskit import *
from qiskit import IBMQ, Aer, assemble, transpile
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

# import basic plot tools
from qiskit.visualization import plot_histogram
# função geradora dos binários:
def to_binary_with_sign(num):
    if int(num) >= 0:
        binary = bin(int(num))[2:]
        binary_with_sign = '0' + binary.zfill(len(binary))
    else:
        binary = bin(-int(num))[2:]
        binary_with_sign = '1' + binary.zfill(len(binary))
    return binary_with_sign
def QRAM(Dq,sizehvalue,I,has_negative): #Dq é o vetor de dado, sizehvalue é o número de elementos esperado do maior número em binário
    #Variáveis
    sizehvalue += 1
    size = len(Dq)
    qsize=int(np.ceil(len(Dq)/2))
    fsize = qsize+sizehvalue+1
    pvector = list(range(qsize))
    #Criando o circuito
    qr = QuantumRegister(fsize)
    cr = ClassicalRegister(qsize)
    circuit = QuantumCircuit(qr,cr)
    #Inicialização
    for h in range(qsize):
        circuit.h(h)
    circuit.h(fsize-1)
    circuit.z(fsize-1)
    circuit.barrier()
    #Memória
    for i in range(size):
        for k in range(qsize):
            if k < len(bin(i)):
                if bin(i)[len(bin(i))-1-k] != "1":
                    circuit.x(k)
            else:
                circuit.x(k)
        for j in range(sizehvalue):
            if len(to_binary_with_sign(int(Dq[i]))) == j+1 and I == 1:
                if to_binary_with_sign(int(Dq[i]))[len(to_binary_with_sign(int(Dq[i])))-1-j] == "1" and qsize!=1:
                    circuit.mct(pvector,qsize)
                    has_negative = True
                if to_binary_with_sign(int(Dq[i]))[len(to_binary_with_sign(int(Dq[i])))-1-j] == "1" and qsize == 1:
                    circuit.cx(pvector,qsize)
                    has_negative = True
            else:
                if len(to_binary_with_sign(int(Dq[i]))) == j and I == 1:
                    break
                if len(to_binary_with_sign(int(Dq[i]))) == j+1 and I != 1:
                    break
                if to_binary_with_sign(int(Dq[i]))[len(to_binary_with_sign(int(Dq[i])))-1-j] == "1" and qsize!=1:
                    circuit.mct(pvector,sizehvalue+qsize-1-j)
                if to_binary_with_sign(int(Dq[i]))[len(to_binary_with_sign(int(Dq[i])))-1-j] == "1" and qsize == 1:
                    circuit.cx(pvector,sizehvalue+qsize-1-j)
        for p in range(qsize):
            if p < len(bin(i)):
                if bin(i)[len(bin(i))-1-p] != "1":
                    circuit.x(p)
            else:
                circuit.x(p)
        circuit.barrier()
    #Operador P
    circuit.x(qsize)
    circuit.cx(qsize,fsize-1)
    circuit.x(qsize)
    circuit.barrier()
    #Memória "dagger"
    for i in range(size):
        for k in range(qsize):
            if k < len(bin(size-1-i)):
                if bin(size-1-i)[len(bin(size-1-i))-1-k] != "1":
                    circuit.x(k)
            else:
                circuit.x(k)
        for j in range(sizehvalue):
            if len(to_binary_with_sign(Dq[size-1-i])) == j+1 and I == 1:
                if to_binary_with_sign(Dq[size-1-i])[len(to_binary_with_sign(Dq[size-1-i]))-1-j] == "1" and qsize!=1:
                    circuit.mct(pvector,qsize)
                if to_binary_with_sign(Dq[size-1-i])[len(to_binary_with_sign(Dq[size-1-i]))-1-j] == "1" and qsize == 1:
                    circuit.cx(pvector,qsize)
            else:
                if len(to_binary_with_sign(Dq[size-1-i])) == j and I == 1:
                    break
                if len(to_binary_with_sign(Dq[size-1-i])) == j+1 and I != 1:
                    break
                if to_binary_with_sign(Dq[size-1-i])[len(to_binary_with_sign(Dq[size-1-i]))-1-j] == "1" and qsize!=1:
                    circuit.mct(pvector,sizehvalue+qsize-1-j)
                if to_binary_with_sign(Dq[size-1-i])[len(to_binary_with_sign(Dq[size-1-i]))-1-j] == "1" and qsize==1:
                    circuit.cx(pvector,sizehvalue+qsize-1-j)
        for p in range(qsize):
            if p < len(bin(size-1-i)):
                if bin(size-1-i)[len(bin(size-1-i))-1-p] != "1":
                    circuit.x(p)
            else:
                circuit.x(p)
        circuit.barrier()
    #Diffuser
    for n1 in range(qsize):
        circuit.h(n1)
    for n2 in range(qsize):
        circuit.x(n2)
    if qsize == 1:
        circuit.z(0)
    elif qsize == 2:
        circuit.cz(0,1)
    elif qsize != 1 and qsize != 2:
        circuit.h(qsize-1)
        circuit.mct(list(range(qsize-1)), qsize-1,[])
        circuit.h(qsize-1)
    for n3 in range(qsize):
        circuit.x(n3)
    if qsize != 1 or I==1:
        for n4 in range(qsize):
            circuit.h(n4)
    #Measure
    circuit.measure(pvector,cr)
    #Simulador
    backend = Aer.get_backend('aer_simulator')
    shots = 1024*(int(size))
    results = execute(circuit, backend=backend, shots=shots).result()
    #Post processing
    counts = list(results.get_counts().values())
    binaries = list(results.get_counts().keys())
    Dq_new = []
    if len(binaries)>2:
        for d in range(len(counts)):
            if counts[d]*10 < 1024*(int(size)) and int(binaries[d],2) < len(Dq) and I != 1 and has_negative == False:
                Dq_new.append(Dq[int(binaries[d],2)])
            elif counts[d]*10 > 1024*(int(size)) and int(binaries[d],2) < len(Dq) and I != 1 and has_negative == True:
                Dq_new.append(Dq[int(binaries[d],2)])
            elif counts[d]*10 > 1024*(int(size)) and int(binaries[d],2) < len(Dq) and I == 1:
                Dq_new.append(Dq[int(binaries[d],2)])
    elif len(binaries)==2:
        if I == 1:
            Dq_new.append(Dq[int(binaries[0],2)])
            Dq_new.append(Dq[int(binaries[1],2)])
        elif I != 1 and has_negative == False:
            del Dq[int(binaries[0],2)]
            if int(binaries[0],2) < int(binaries[1],2):
                del Dq[int(binaries[1],2)-1]
            else:
                del Dq[int(binaries[1],2)]
            Dq_new = Dq
        elif I != 1 and has_negative == True:
            Dq_new.append(Dq[int(binaries[0],2)])
            Dq_new.append(Dq[int(binaries[1],2)])
    elif len(binaries)==1:
        if I == 1:
            Dq_new.append(Dq[int(binaries[0],2)])
        elif I!=1 and has_negative == False:
            del Dq[int(binaries[0],2)]
            Dq_new = Dq
        elif I!=1 and has_negative == True:
            Dq_new.append(Dq[int(binaries[0],2)])
    if len(Dq_new) == 0:
        Dq_new = Dq
    print(Dq_new)
    return Dq_new, has_negative

#Previous-processing
V = list(input("Insira o vetor separado por espaço: ").split())
for i in range(len(V)):
    V[i] = int(V[i])
H = int(input("Número de elementos do binário do maior número em módulo: "))
I = 1
B = False
while len(V)>3 and H!=0 and (sum(V)/len(V)) != V[0]:
    V, B = QRAM(V,H,I,B)
    I += 1
    H -= 1