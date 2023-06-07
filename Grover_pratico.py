# initialization
import matplotlib.pyplot as plt
import numpy as np
import time

# importing Qiskit
from qiskit.tools.visualization import plot_bloch_multivector
from qiskit.quantum_info import Statevector
from qiskit import *
from qiskit import IBMQ, Aer, assemble, transpile
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

# import basic plot tools
from qiskit.visualization import plot_histogram
#Função manual de ccx
def ccx(circuit, control_qubit1, control_qubit2, target_qubit):
    circuit.h(target_qubit)
    circuit.cx(control_qubit2, target_qubit)
    circuit.tdg(target_qubit)
    circuit.cx(control_qubit1, target_qubit)
    circuit.t(target_qubit)
    circuit.cx(control_qubit2, target_qubit)
    circuit.tdg(target_qubit)
    circuit.cx(control_qubit1, target_qubit)
    circuit.t(control_qubit2)
    circuit.t(target_qubit)
    circuit.cx(control_qubit1, control_qubit2)
    circuit.h(target_qubit)
    circuit.t(control_qubit1)
    circuit.tdg(control_qubit2)
    circuit.cx(control_qubit1, control_qubit2)
    
#Função manual de mct
def mct(circuit, control_qubits, target_wire):
    num_ctrl = len(control_qubits)

    if num_ctrl == 1:
        circuit.cx(control_qubits[0], target_wire)
    elif num_ctrl == 2:
        ccx(circuit,control_qubits[0], control_qubits[1], target_wire)
    else:
        ancilla = QuantumRegister(num_ctrl - 2)
        circuit.add_register(ancilla)

        ccx(circuit,control_qubits[0], control_qubits[1], ancilla[0])

        for i in range(2, num_ctrl - 1):
            ccx(circuit,control_qubits[i], ancilla[i-2], ancilla[i-1])

        ccx(circuit,control_qubits[num_ctrl-1], ancilla[num_ctrl-3], target_wire)

        for i in range(num_ctrl - 2, 1, -1):
            ccx(circuit,control_qubits[i], ancilla[i-2], ancilla[i-1])

        ccx(circuit,control_qubits[0], control_qubits[1], ancilla[0])


# função geradora dos binários:
def to_binary_with_sign(num):
    if int(num) >= 0:
        binary = bin(int(num))[2:]
        binary_with_sign = '0' + binary.zfill(len(binary))
    else:
        binary = bin(-int(num))[2:]
        binary_with_sign = '1' + binary.zfill(len(binary))
    return binary_with_sign
def QRAM(Dq,sizehvalue,I,has_negative,A): #Dq é o vetor de dado, sizehvalue é o número de elementos esperado do maior número em binário
    #Variáveis
    sizehvalue += 1
    size = len(Dq)
    qsize=int(len(to_binary_with_sign(len(Dq)))-1)
    fsize = qsize+sizehvalue+1
    pvector = list(range(qsize))
    Dq_new = []
    numq = 0
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
            if len(to_binary_with_sign(int((Dq[i])/A))) == j+1 and I == 1:
                if to_binary_with_sign(int((Dq[i])/A))[len(to_binary_with_sign(int((Dq[i])/A)))-1-j] == "1" and qsize!=1:
                    mct(circuit,pvector,qsize)
                    has_negative = True
                if to_binary_with_sign(int((Dq[i])/A))[len(to_binary_with_sign(int((Dq[i])/A)))-1-j] == "1" and qsize == 1:
                    circuit.cx(pvector,qsize)
                    has_negative = True
            else:
                if len(to_binary_with_sign(int((Dq[i])/A))) == j and I == 1:
                    break
                if len(to_binary_with_sign(int((Dq[i])/A))) == j+1 and I != 1:
                    break
                if to_binary_with_sign(int((Dq[i])/A))[len(to_binary_with_sign(int((Dq[i])/A)))-1-j] == "1" and qsize!=1:
                    mct(circuit,pvector,sizehvalue+qsize-1-j)
                if to_binary_with_sign(int((Dq[i])/A))[len(to_binary_with_sign(int((Dq[i])/A)))-1-j] == "1" and qsize == 1:
                    circuit.cx(pvector,sizehvalue+qsize-1-j)
        for p in range(qsize):
            if p < len(bin(i)):
                if bin(i)[len(bin(i))-1-p] != "1":
                    circuit.x(p)
            else:
                circuit.x(p)
        circuit.barrier()
    if has_negative == False:
        return Dq_new, has_negative, numq
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
            if len(to_binary_with_sign((Dq[size-1-i])/A)) == j+1 and I == 1:
                if to_binary_with_sign((Dq[size-1-i])/A)[len(to_binary_with_sign((Dq[size-1-i])/A))-1-j] == "1" and qsize!=1:
                    mct(circuit,pvector,qsize)
                if to_binary_with_sign((Dq[size-1-i])/A)[len(to_binary_with_sign((Dq[size-1-i])/A))-1-j] == "1" and qsize == 1:
                    circuit.cx(pvector,qsize)
            else:
                if len(to_binary_with_sign((Dq[size-1-i])/A)) == j and I == 1:
                    break
                if len(to_binary_with_sign((Dq[size-1-i])/A)) == j+1 and I != 1:
                    break
                if to_binary_with_sign((Dq[size-1-i])/A)[len(to_binary_with_sign((Dq[size-1-i])/A))-1-j] == "1" and qsize!=1:
                    mct(circuit,pvector,sizehvalue+qsize-1-j)
                if to_binary_with_sign((Dq[size-1-i])/A)[len(to_binary_with_sign((Dq[size-1-i])/A))-1-j] == "1" and qsize==1:
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
        mct(circuit,list(range(qsize-1)),qsize-1)
        circuit.h(qsize-1)
    for n3 in range(qsize):
        circuit.x(n3)
    if qsize != 1 or I==1:
        for n4 in range(qsize):
            circuit.h(n4)
    #Measure
    circuit.measure(pvector,cr)
    numq = circuit.num_qubits
    #Simulador
    backend = Aer.get_backend('aer_simulator')
    shots = 50*(2**qsize)
    results = execute(circuit, backend=backend, shots=shots).result()
    #Post processing
    counts = list(results.get_counts().values())
    binaries = list(results.get_counts().keys())
    mean = shots/len(counts)
    if len(binaries)>2:
        for d in range(len(counts)):
            if counts[d] < mean + 20 and int(binaries[d],2) < len(Dq) and I != 1 and has_negative == False:
                Dq_new.append(Dq[int(binaries[d],2)])
            elif counts[d] > mean + 20 and int(binaries[d],2) < len(Dq) and I != 1 and has_negative == True:
                Dq_new.append(Dq[int(binaries[d],2)])
            elif counts[d] > mean + 20 and int(binaries[d],2) < len(Dq) and I == 1 and has_negative == True:
                Dq_new.append(Dq[int(binaries[d],2)])
            elif counts[d] < mean + 20 and int(binaries[d],2) < len(Dq) and I == 1 and has_negative == False:
                Dq_new.append(Dq[int(binaries[d],2)])
    elif len(binaries)==2:
        if I == 1 and int(binaries[0],2) < len(Dq) and int(binaries[1],2) < len(Dq):
            Dq_new.append(Dq[int(binaries[0],2)])
            Dq_new.append(Dq[int(binaries[1],2)])
        elif I != 1 and has_negative == False and int(binaries[0],2) < len(Dq) and int(binaries[1],2) < len(Dq):
            del Dq[int(binaries[0],2)]
            if int(binaries[0],2) < int(binaries[1],2):
                del Dq[int(binaries[1],2)-1]
            else:
                del Dq[int(binaries[1],2)]
            Dq_new = Dq
        elif I != 1 and has_negative == True and int(binaries[0],2) < len(Dq) and int(binaries[1],2) < len(Dq):
            Dq_new.append(Dq[int(binaries[0],2)])
            Dq_new.append(Dq[int(binaries[1],2)])
    elif len(binaries)==1:
        if I == 1 and int(binaries[0],2) < len(Dq):
            Dq_new.append(Dq[int(binaries[0],2)])
        elif I!=1 and has_negative == False and int(binaries[0],2) < len(Dq):
            del Dq[int(binaries[0],2)]
            Dq_new = Dq
        elif I!=1 and has_negative == True and int(binaries[0],2) < len(Dq):
            Dq_new.append(Dq[int(binaries[0],2)])
    if len(Dq_new) == 0:
        Dq_new = Dq
    return Dq_new, has_negative, numq

#Previous-processing
data = np.genfromtxt("results-ex4.txt",dtype=float,delimiter=";")
V_full = data[:,1]
H_inicial = int(input("Número de elementos do binário do maior número em módulo: "))
length = int(input("Quer dividir em blocos de que tamanho?: "))
acc = int((2**(H_inicial)))
print(f"Fator de precisão: {int((2**(H_inicial)))}")
Done = False
Q = []
Qaux = 1
start = time.time() #Início do tempo
while Done == False:
    H_inicial_linha = int(H_inicial - int(np.log2(acc)))
    V_final = []
    first = 0
    last = length
    i = 0
    while first != len(V_full):
        print(f"Bloco {i+1}")
        V = V_full[first:last]
        H = H_inicial_linha
        I = 1
        B = False
        A = acc
        while len(V)>3 and H!=-1 and int(sum(V)/len(V)) != int(V[0]):
            V,B,Qaux = QRAM(V,H,I,B,A)
            Q.append(Qaux)
            I += 1
            H -= 1
        if (len(V)>1 and int(sum(V)/len(V)) == int(V[0])) and B == True:
            V_final.append(V[0])
        elif len(V) > 1 and B == True:
            A = int(2**(np.log2(acc)-1))
            H = int(H_inicial - int(np.log2(A)))
            I = 1
            B = False
            while len(V)>3 and H!=-1 and int(sum(V)/len(V)) != int(V[0]):
                V,B,Qaux = QRAM(V,H,I,B,A)
                Q.append(Qaux)
                I += 1
                H -= 1
            if len(V)>=1:
                V_final.append(np.min(V))
        first = last
        if last + length > len(V_full):
            last = len(V_full)
        else:
            last += length
        i += 1
    #CASO USEMOS OS VALORES MODIFICADOS
    if acc != 1 and int(len(V_final)) != 0:
        V = np.array(V_final)
        V = np.min(V) #ajuste de precisão
    #CASO USEMOS OS VALORES IMACULADOS
    elif acc == 1 and int(len(V_final)) != 0:
        V = np.array(V_final)
        H = H_inicial
        I = 1
        B = False
        while len(V)>3 and H!=-1 and int(sum(V)/len(V)) != int(V[0]):
            V, B, Qaux = QRAM(V,H,I,B,A)
            Q.append(Qaux)
            I += 1
            H -= 1
        V = np.min(V) #ajuste de precisão
    if len(V_final)==0:
        print(f"Erro de precisão")
        print(f"Ajustarei o Fator de Precisão para {int(2**(np.log2(acc)-1))}")
        acc = int(2**(np.log2(acc)-1))
    else:
        end = time.time() #fim do tempo
        for j in range(len(V_full)):
            if V_full[j] == V:
                cen = j+1
        print(f"Pior retorno: {V}")
        print(f"Pior cenário: {cen}")
        print(f"Número máximo de Qubits usados: {np.max(Q)}")
        print(f"Tempo de processamento: {(end - start)} seg")
        Done = True