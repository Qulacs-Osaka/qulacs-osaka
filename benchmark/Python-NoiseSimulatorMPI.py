#
# Benchmark based on: https://github.com/keisukefujii/LXEB_qulacs/blob/master/LXEB_qulacs.ipynb
#
from mpi4py import MPI # MPI用ライブラリのmpi4pyが必要です。

import numpy as np
import time 
import random

from qulacs import NoiseSimulatorMPI

from qulacs import QuantumState
#GPU版をインストールしている場合のみ
#from qulacs import QuantumStateGpu
from qulacs import QuantumCircuit
from qulacs.gate import DenseMatrix
from qulacs.circuit import QuantumCircuitOptimizer
from qulacs.state import inner_product
from qulacs.gate import Identity, X,Y,Z #パウリ演算子
from qulacs.gate import H,S,Sdag, sqrtX,sqrtXdag,sqrtY,sqrtYdag #1量子ビット Clifford演算
from qulacs.gate import T,Tdag #1量子ビット 非Clifford演算
from qulacs.gate import RX,RY,RZ #パウリ演算子についての回転演算
from qulacs.gate import CNOT, CZ, SWAP #2量子ビット演算
from qulacs.gate import DephasingNoise,DepolarizingNoise,TwoQubitDepolarizingNoise


probability = 0.001

#xy座標から通し番号を出力
def GridtoId(x,y,length):
    return x+length*y

#通し番号からxy座標を出力
def IdtoGrid(id,length):
    return np.array([id%length,(id-id%length)/length],dtype = int)

#グリッド内にいるかどうか判定
def InGrid(x,y,length):
    if x<length and y<length and x>=0 and y>=0:
        return 0
    else:
        return 1
def grid_CZ_gate(qubit1,qubit2,length):
    if InGrid(qubit1[0],qubit1[1],length)==0 and InGrid(qubit2[0],qubit2[1],length)==0:
        circuit.add_gate(CZ(GridtoId(qubit1[0],qubit1[1],length),GridtoId(qubit2[0],qubit2[1],length)))

#グリッド上で2qubit gateを作用
def grid_2Q_gate(qubit1,qubit2,length,gate_array,circuit):
    if InGrid(qubit1[0],qubit1[1],length)==0 and InGrid(qubit2[0],qubit2[1],length)==0:
        circuit.add_noise_gate(DenseMatrix([GridtoId(qubit1[0],qubit1[1],length),GridtoId(qubit2[0],qubit2[1],length)],gate_array),"Depolarizing",probability)  

#グリッド上で1qubit gateを作用
def grid_1Q_gate(qubitId,circuit,gate_type):
    if gate_type == 0:
        circuit.add_gate(sqrtX(qubitId))
    elif gate_type == 1:
        circuit.add_gate(sqrtY(qubitId))
    elif gate_type == 2:
        circuit.add_gate(DenseMatrix(qubitId,[[1/np.sqrt(2.0),- np.sqrt(1.0j/2.0)],[np.sqrt(-1.0j/2.0),1/np.sqrt(2.0)]]))


#グリッド上でrandom 1Q gateを作用        
def random_1Q_gate(qubitId,circuit):
    hoge = 1
    prob = random.random()
    if prob < 1.0/3.0:
        circuit.add_noise_gate(sqrtX(qubitId),"Depolarizing",probability)
    elif prob < 2.0/3.0:
        circuit.add_noise_gate(sqrtY(qubitId),"Depolarizing",probability)
    else:
        circuit.add_noise_gate(DenseMatrix(qubitId,[[1/np.sqrt(2.0),- np.sqrt(1.0j/2.0)],[np.sqrt(-1.0j/2.0),1/np.sqrt(2.0)]]),"Depolarizing",probability)
        #circuit.add_gate(T(qubitId))

theta = np.pi/2.0
iSwapLikeGate = [[1,0,0,0],
                [0,np.cos(theta),-1.0j*np.sin(theta),0],
                [0,-1.0j*np.sin(theta),np.cos(theta),0],
                [0,0,0,1]]


#テスト用回路：アダマール演算を用いて全てのビット列の重ね合わせ状態を生成
def test_circuit(length,circuit):
    num_qubits = length**2
    for i in range(num_qubits):
        circuit.add_gate(H(i))

def Google_random_circuit(length,depth,circuit):
    num_qubits = length**2

    for k in range(depth):
        if k%4 == 0:
            #cycle A
            #ランダムに1量子ビット演算を追加
            for i in range(num_qubits):
                random_1Q_gate(i,circuit)

            #２量子ビットを追加
            for i in range(length):
                for j in range(length):
                    if (i+j)%2 == 0:
                        grid_2Q_gate([i,j],[i+1,j],length,iSwapLikeGate,circuit)
                        

        if k%4 == 1:
            #cycle B
            #ランダムに1量子ビット演算を追加
            for i in range(num_qubits):
                random_1Q_gate(i,circuit)

            #２量子ビットを追加
            for i in range(length):
                for j in range(length):
                    if (i+j)%2 == 0:
                        grid_2Q_gate([i,j],[i-1,j],length,iSwapLikeGate,circuit)

        if k%4 == 2:
            #cycle C
            #ランダムに1量子ビット演算を追加
            for i in range(num_qubits):
                random_1Q_gate(i,circuit)

            #２量子ビットを追加
            for i in range(length):
                for j in range(length):
                    if (i+j)%2 == 0:
                        grid_2Q_gate([i,j],[i,j+1],length,iSwapLikeGate,circuit)

        if k%4 == 3:
            #cycle D
            #ランダムに1量子ビット演算を追加
            for i in range(num_qubits):
                random_1Q_gate(i,circuit)

            #２量子ビットを追加
            for i in range(length):
                for j in range(length):
                    if (i+j)%2 == 0:
                        grid_2Q_gate([i,j],[i,j-1],length,iSwapLikeGate,circuit)


                        
num_samp = 100000
length = 4
nqubits = length**2
depth = 20

state = QuantumState(nqubits)
circuit = QuantumCircuit(nqubits)

Google_random_circuit(length,depth,circuit)

start = time.time()

# (コンストラクタ)NoiseSimulatorMPI(circuit,state) -> circuit: 量子回路 state: 初期状態
hoge = NoiseSimulatorMPI(circuit,state)

# NoiseSimulatorMPI.execute(num_samp) -> num_samp: サンプル数, 返り値: rank=0のノードにすべてのサンプル結果が返ってきます。それ以外は空配列が返ってきます。
ans = hoge.execute(num_samp)

elapsed_time = time.time() - start
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
if len(ans) != 0:
    for i in range(10):
        print(ans[i * int(len(ans) / 10)])
print(len(ans))
print(circuit)
