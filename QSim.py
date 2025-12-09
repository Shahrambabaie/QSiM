
# 3-Qubit Quantum Simulator

import numpy as np
import array as arr
import random
import matplotlib.pyplot as plt
import cmath
import math

class Circuit712():
  q0 = np.array([[1.+0.j],[0.+0.j]])
  q1 = np.array([[1.+0.j],[0.+0.j]])
  q2 = np.array([[1.+0.j],[0.+0.j]])
  MCounter = 0
  q3bit = np.array([[0.+0.j],[0.+0.j],[0.+0.j],[0.+0.j],[0.+0.j],[0.+0.j],[0.+0.j],[0.+0.j]])
  qubit0 = np.array([[1.+0.j],[0.+0.j]])
  qubit1 = np.array([[0.+0.j],[1.+0.j]])
  PMVector = np.array([0,0,0,0,0,0,0,0])
  q0GateList = []
  q1GateList = []
  q2GateList = []
  list_measurement = []
  counts = arr.array('i', [0, 0, 0, 0, 0, 0, 0, 0])

  def __init__(self, number_of_qubits: int):
    self._number_of_qubits = number_of_qubits
    if not isinstance(self._number_of_qubits, int):
      raise TypeError("The number_of_qubits parameter Qubit is not an :py:class:`int`")

    if self._number_of_qubits < 1:
      raise ValueError("The number_of_qubits parameter " + "must be greater than 0")

    Circuit712.MCounter = Circuit712.MCounter
    Circuit712.q0GateList = Circuit712.q0GateList
    Circuit712.q1GateList = Circuit712.q1GateList
    Circuit712.q2GateList = Circuit712.q2GateList
    Circuit712.q0GateList.append("q0:   |0>-----")
    Circuit712.q1GateList.append("q1:   |0>-----")
    Circuit712.q2GateList.append("q2:   |0>-----")
    Circuit712.counts = Circuit712.counts
    Circuit712.PMVector = Circuit712.PMVector
    Circuit712.list_measurement = Circuit712.list_measurement
    Circuit712.q3bit = Circuit712.q3bit
    M1 = np.kron(Circuit712.q0, Circuit712.q1)
    Circuit712.q3bit = np.kron(Circuit712.q2, M1)
    Circuit712.PMVector = Circuit712.q3bit

#----------------------------------------------------------
# 1-qubit Quantum Gates
#----------------------------------------------------------

  def id(self, qubit: int):
    idmatrix = np.array([[1.+0.j,0.+0.j],[0.+0.j,1.+0.j]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    id1 = np.kron(idmatrix,idmatrix)
    id2 = np.kron(id1, idmatrix)

    if qubit == 0:
      Circuit712.q3bit = np.dot(id2,Circuit712.q3bit)

    elif qubit == 1:
        Circuit712.q3bit = np.dot(id2,Circuit712.q3bit)

    elif qubit == 2:
        Circuit712.q3bit = np.dot(id2,Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit
    if qubit == 0:
      Circuit712.q0GateList.append("id")
      Circuit712.q1GateList.append("--")
      Circuit712.q2GateList.append("--")

    if qubit == 1:
      Circuit712.q1GateList.append("id")
      Circuit712.q0GateList.append("--")
      Circuit712.q2GateList.append("--")

    if qubit == 2:
      Circuit712.q2GateList.append("id")
      Circuit712.q1GateList.append("--")
      Circuit712.q0GateList.append("--")

  def x(self, qubit: int):
    Xmatrix = np.array([[0.+0.j,1.+0.j],[1.+0.j,0.+0.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      x0 = np.kron(Imatrix,Imatrix)
      x1 = np.kron(x0, Xmatrix)
      Circuit712.q3bit = np.dot(x1,Circuit712.q3bit)

    elif qubit == 1:
      x2 = np.kron(Imatrix,Xmatrix)
      x3 = np.kron(x2, Imatrix)
      Circuit712.q3bit = np.dot(x3, Circuit712.q3bit)

    elif qubit == 2:
      x4 = np.kron(Xmatrix, Imatrix)
      x5 = np.kron(x4, Imatrix)
      Circuit712.q3bit = np.dot(x5, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("X")
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 1:
      Circuit712.q1GateList.append("X")
      Circuit712.q0GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 2:
      Circuit712.q2GateList.append("X")
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")

  def z(self, qubit: int):
    Zmatrix = np.array([[1.+0.j,0.+0.j],[0.+0.j,-1.+0.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      z0 = np.kron(Imatrix,Imatrix)
      z1 = np.kron(z0, Zmatrix)
      Circuit712.q3bit = np.dot(z1,Circuit712.q3bit)

    elif qubit == 1:
      z2 = np.kron(Imatrix,Zmatrix)
      z3 = np.kron(z2, Imatrix)
      Circuit712.q3bit = np.dot(z3, Circuit712.q3bit)

    elif qubit == 2:
      z4 = np.kron(Zmatrix, Imatrix)
      z5 = np.kron(z4, Imatrix)
      Circuit712.q3bit = np.dot(z5, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("Z")
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 1:
      Circuit712.q1GateList.append("Z")
      Circuit712.q0GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 2:
      Circuit712.q2GateList.append("Z")
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")

  def y(self, qubit: int):
    Ymatrix = np.array([[0.+0.j,0.-1.j],[0.+1.j,0.+0.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      y0 = np.kron(Imatrix,Imatrix)
      y1 = np.kron(y0, Ymatrix)
      Circuit712.q3bit = np.dot(y1,Circuit712.q3bit)

    elif qubit == 1:
      y2 = np.kron(Imatrix,Ymatrix)
      y3 = np.kron(y2, Imatrix)
      Circuit712.q3bit = np.dot(y3, Circuit712.q3bit)

    elif qubit == 2:
      y4 = np.kron(Ymatrix, Imatrix)
      y5 = np.kron(y4, Imatrix)
      Circuit712.q3bit = np.dot(y5, Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("Y")
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 1:
      Circuit712.q1GateList.append("Y")
      Circuit712.q0GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 2:
      Circuit712.q2GateList.append("Y")
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")

  def s(self, qubit: int):
    Smatrix = np.array([[1.+0.j,0.+0.j],[0.+0.j,0.-1.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      s0 = np.kron(Imatrix,Imatrix)
      s1 = np.kron(s0, Smatrix)
      Circuit712.q3bit = np.dot(s1,Circuit712.q3bit)

    elif qubit == 1:
      s2 = np.kron(Imatrix,Smatrix)
      s3 = np.kron(s2, Imatrix)
      Circuit712.q3bit = np.dot(s3, Circuit712.q3bit)

    elif qubit == 2:
      s4 = np.kron(Smatrix, Imatrix)
      s5 = np.kron(s4, Imatrix)
      Circuit712.q3bit = np.dot(s5, Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("S")
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 1:
      Circuit712.q1GateList.append("S")
      Circuit712.q0GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 2:
      Circuit712.q2GateList.append("S")
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")


  def s_adjoint(self, qubit: int):
    SADmatrix = np.array([[1.+0.j,0.+0.j],[0.+0.j,0.-1.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      ST0 = np.kron(Imatrix,Imatrix)
      ST1 = np.kron(ST0, SADmatrix)
      Circuit712.q3bit = np.dot(ST1,Circuit712.q3bit)

    elif qubit == 1:
      ST2 = np.kron(Imatrix,SADmatrix)
      ST3 = np.kron(ST2, Imatrix)
      Circuit712.q3bit = np.dot(ST3, Circuit712.q3bit)

    elif qubit == 2:
      ST4 = np.kron(SADmatrix, Imatrix)
      ST5 = np.kron(ST4, Imatrix)
      Circuit712.q3bit = np.dot(ST5, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("S†")
      Circuit712.q1GateList.append("--")
      Circuit712.q2GateList.append("--")

    if qubit == 1:
      Circuit712.q1GateList.append("S†")
      Circuit712.q0GateList.append("--")
      Circuit712.q2GateList.append("--")

    if qubit == 2:
      Circuit712.q2GateList.append("S†")
      Circuit712.q1GateList.append("--")
      Circuit712.q0GateList.append("--")

  def sqrt_not(self, qubit: int):
    SQRTNOTmatrix = np.array([[1/2.+1/2.j,1/2.-1/2.j],[1/2.-1/2.j,1/2.+1/2.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      SQ0 = np.kron(Imatrix,Imatrix)
      SQ1 = np.kron(SQ0, SQRTNOTmatrix)
      Circuit712.q3bit = np.dot(SQ1,Circuit712.q3bit)

    elif qubit == 1:
      SQ2 = np.kron(Imatrix,SQRTNOTmatrix)
      SQ3 = np.kron(SQ2, Imatrix)
      Circuit712.q3bit = np.dot(SQ3, Circuit712.q3bit)

    elif qubit == 2:
      SQ4 = np.kron(SQRTNOTmatrix, Imatrix)
      SQ5 = np.kron(SQ4, Imatrix)
      Circuit712.q3bit = np.dot(SQ5, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("√NOT")
      Circuit712.q1GateList.append("----")
      Circuit712.q2GateList.append("----")

    if qubit == 1:
      Circuit712.q1GateList.append("√NOT")
      Circuit712.q0GateList.append("----")
      Circuit712.q2GateList.append("----")

    if qubit == 2:
      Circuit712.q2GateList.append("√NOT")
      Circuit712.q1GateList.append("----")
      Circuit712.q0GateList.append("----")

  def t(self, qubit: int):
    Tmatrix = np.array([[1.+0.j,0.+0.j],[0.+0.j,np.sqrt(2)/2+np.sqrt(2)/2.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")
    if (qubit < 0 or qubit >= self._number_of_qubits):

      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      t0 = np.kron(Imatrix,Imatrix)
      t1 = np.kron(t0, Tmatrix)
      Circuit712.q3bit = np.dot(t1,Circuit712.q3bit)

    elif qubit == 1:
      t2 = np.kron(Imatrix,Tmatrix)
      t3 = np.kron(t2, Imatrix)
      Circuit712.q3bit = np.dot(t3, Circuit712.q3bit)

    elif qubit == 2:
      t4 = np.kron(Tmatrix, Imatrix)
      t5 = np.kron(t4, Imatrix)
      Circuit712.q3bit = np.dot(t5, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("T")
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 1:
      Circuit712.q1GateList.append("T")
      Circuit712.q0GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 2:
      Circuit712.q2GateList.append("T")
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")

  def t_adjoint(self, qubit: int):
    TADmatrix = np.array([[1.+0.j,0.+0.j],[0.+0.j,np.sqrt(2)/2-np.sqrt(2)/2.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      TA0 = np.kron(Imatrix,Imatrix)
      TA1 = np.kron(TA0, TADmatrix)
      Circuit712.q3bit = np.dot(TA1,Circuit712.q3bit)

    elif qubit == 1:
      TA2 = np.kron(Imatrix,TADmatrix)
      TA3 = np.kron(TA2, Imatrix)
      Circuit712.q3bit = np.dot(TA3, Circuit712.q3bit)

    elif qubit == 2:
      TA4 = np.kron(TADmatrix, Imatrix)
      TA5 = np.kron(TA4, Imatrix)
      Circuit712.q3bit = np.dot(TA5, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("T†")
      Circuit712.q1GateList.append("--")
      Circuit712.q2GateList.append("--")

    if qubit == 1:
      Circuit712.q1GateList.append("T†")
      Circuit712.q0GateList.append("--")
      Circuit712.q2GateList.append("--")

    if qubit == 2:
      Circuit712.q2GateList.append("T†")
      Circuit712.q1GateList.append("--")
      Circuit712.q0GateList.append("--")


  def rx(self, qubit: int, phi: float):
    radian = (np.pi*phi/3.14)
    degree = int((np.pi/phi))
    RXmatrix = np.array([[cmath.cos(radian/2), complex(0,-cmath.sin(radian/2))],[complex(0,-cmath.sin(radian/2)),cmath.cos(radian/2)]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if not isinstance(phi, float):
      raise TypeError("PHI is not an :py:class:`float`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      RX0 = np.kron(Imatrix,Imatrix)
      RX1 = np.kron(RX0, RXmatrix)
      Circuit712.q3bit = np.dot(RX1,Circuit712.q3bit)

    elif qubit == 1:
      RX2 = np.kron(Imatrix,RXmatrix)
      RX3 = np.kron(RX2, Imatrix)
      Circuit712.q3bit = np.dot(RX3, Circuit712.q3bit)

    elif qubit == 2:
      RX4 = np.kron(RXmatrix, Imatrix)
      RX5 = np.kron(RX4, Imatrix)
      Circuit712.q3bit = np.dot(RX5, Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("RX(π/%i)"%degree)
      Circuit712.q1GateList.append("-------")
      Circuit712.q2GateList.append("-------")

    if qubit == 1:
      Circuit712.q1GateList.append("RX(π/%i)"%degree)
      Circuit712.q0GateList.append("-------")
      Circuit712.q2GateList.append("-------")

    if qubit == 2:
      Circuit712.q2GateList.append("RX(π/%i)"%degree)
      Circuit712.q1GateList.append("-------")
      Circuit712.q0GateList.append("-------")


  def ry(self, qubit: int, phi: float):
    radian = (np.pi*phi/3.14)
    degree = int((np.pi/phi))
    RYmatrix = np.array([[np.cos(radian/2), -np.sin(radian/2)],[np.sin(radian/2),np.cos(radian/2)]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if not isinstance(phi, float):
      raise TypeError("PHI is not an :py:class:`float`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      RY0 = np.kron(Imatrix,Imatrix)
      RY1 = np.kron(RY0, RYmatrix)
      Circuit712.q3bit = np.dot(RY1,Circuit712.q3bit)

    elif qubit == 1:
      RY2 = np.kron(Imatrix,RYmatrix)
      RY3 = np.kron(RY2, Imatrix)
      Circuit712.q3bit = np.dot(RY3, Circuit712.q3bit)

    elif qubit == 2:
      RY4 = np.kron(RYmatrix, Imatrix)
      RY5 = np.kron(RY4, Imatrix)
      Circuit712.q3bit = np.dot(RY5, Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("RY(π/%i)"%degree)
      Circuit712.q1GateList.append("-------")
      Circuit712.q2GateList.append("-------")
    if qubit == 1:
      Circuit712.q1GateList.append("RY(π/%i)"%degree)
      Circuit712.q0GateList.append("-------")
      Circuit712.q2GateList.append("-------")
    if qubit == 2:
      Circuit712.q2GateList.append("RY(π/%i)"%degree)
      Circuit712.q1GateList.append("-------")
      Circuit712.q0GateList.append("-------")


  def rz(self, qubit: int, phi: float):
    radian = (np.pi*phi/3.14)
    degree = int((np.pi/phi))
    RZmatrix = np.array([[cmath.exp(complex(0, -radian)), 0.+0.j],[0.+0.j,cmath.exp(complex(0, radian))]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if not isinstance(phi, float):
      raise TypeError("PHI is not an :py:class:`float`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      RZ0 = np.kron(Imatrix,Imatrix)
      RZ1 = np.kron(RZ0, RZmatrix)
      Circuit712.q3bit = np.dot(RZ1,Circuit712.q3bit)

    elif qubit == 1:
      RZ2 = np.kron(Imatrix,RZmatrix)
      RZ3 = np.kron(RZ2, Imatrix)
      Circuit712.q3bit = np.dot(RZ3, Circuit712.q3bit)

    elif qubit == 2:
      RZ4 = np.kron(RZmatrix, Imatrix)
      RZ5 = np.kron(RZ4, Imatrix)
      Circuit712.q3bit = np.dot(RZ5, Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("RZ(π/%i)"%degree)
      Circuit712.q1GateList.append("-------")
      Circuit712.q2GateList.append("-------")

    if qubit == 1:
      Circuit712.q1GateList.append("RZ(π/%i)"%degree)
      Circuit712.q0GateList.append("-------")
      Circuit712.q2GateList.append("-------")

    if qubit == 2:
      Circuit712.q2GateList.append("RZ(π/%i)"%degree)
      Circuit712.q1GateList.append("-------")
      Circuit712.q0GateList.append("-------")


  def h(self, qubit: int):
    Hmatrix = np.array([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      h0 = np.kron(Imatrix,Imatrix)
      h1 = np.kron(h0, Hmatrix)
      Circuit712.q3bit = np.dot(h1,Circuit712.q3bit)

    elif qubit == 1:
      h2 = np.kron(Imatrix,Hmatrix)
      h3 = np.kron(h2, Imatrix)
      Circuit712.q3bit = np.dot(h3, Circuit712.q3bit)

    elif qubit == 2:
      h4 = np.kron(Hmatrix, Imatrix)
      h5 = np.kron(h4, Imatrix)
      Circuit712.q3bit = np.dot(h5, Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if qubit == 0:
      Circuit712.q0GateList.append("H")
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 1:
      Circuit712.q1GateList.append("H")
      Circuit712.q0GateList.append("-")
      Circuit712.q2GateList.append("-")

    if qubit == 2:
      Circuit712.q2GateList.append("H")
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")

#----------------------------------------------------------
# 2-qubit Quantum Gates
#----------------------------------------------------------

  def ch(self, control_qubit: int, target_qubit: int):
    CHmatrix = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, np.sqrt(2)/2, np.sqrt(2)/2], [0.+0.j,0.+0.j, np.sqrt(2)/2, -np.sqrt(2)/2]])
    CHmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, np.sqrt(2)/2, 0.+0.j,np.sqrt(2)/2], [0.+0.j,0.+0.j, 1.+0j, 0.+0.j], [0.+0.j, np.sqrt(2)/2, 0.+0.j, -np.sqrt(2)/2]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(control_qubit, int):
      raise TypeError("Control_qubit is not an :py:class:`int`")

    if not isinstance(target_qubit, int):
      raise TypeError("Target_qubit is not an :py:class:`int`")

    if (control_qubit < 0 or control_qubit >= self._number_of_qubits):
      raise ValueError ("Control_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (target_qubit < 0 or target_qubit >= self._number_of_qubits):
      raise ValueError ("Target_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (control_qubit == target_qubit):
      raise ValueError ("Target_qubit and Cotrol_qubit can not be IDENTICAL ")

    if (control_qubit == 0 and target_qubit == 1):
        ch0 = np.kron(Imatrix, CHmatrix1)
        Circuit712.q3bit = np.dot(ch0,Circuit712.q3bit)

    if (control_qubit == 0 and target_qubit == 2):
        CHmatrix2 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, np.sqrt(2)/2, 0.+0.j, 0.+0.j, 0.+0.j, np.sqrt(2)/2, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, np.sqrt(2)/2, 0.+0.j, 0.+0.j, 0.+0.j,np.sqrt(2)/2],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, np.sqrt(2)/2, 0.+0.j, 0.+0.j, 0.+0.j,  -np.sqrt(2)/2, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, np.sqrt(2)/2, 0.+0.j, 0.+0.j, 0.+0.j,  -np.sqrt(2)/2 ]])
        Circuit712.q3bit = np.dot(CHmatrix2,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 2):
        ch0 = np.kron(CHmatrix1, Imatrix)
        Circuit712.q3bit = np.dot(ch0,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 0):
        ch0 = np.kron(Imatrix, CHmatrix)
        Circuit712.q3bit = np.dot(ch0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 0):
        CHmatrix3 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j ],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, np.sqrt(2)/2, np.sqrt(2)/2, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, np.sqrt(2)/2, -np.sqrt(2)/2, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, np.sqrt(2)/2, np.sqrt(2)/2,], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, np.sqrt(2)/2, -np.sqrt(2)/2, ]])
        Circuit712.q3bit = np.dot(CHmatrix3,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 1):
        ch0 = np.kron(CHmatrix, Imatrix)
        Circuit712.q3bit = np.dot(ch0,Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if (control_qubit == 0 and target_qubit == 1):
      Circuit712.q0GateList.append("O")
      Circuit712.q1GateList.append("H")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 0 and target_qubit == 2):
      Circuit712.q0GateList.append("O")
      Circuit712.q2GateList.append("H")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 1 and target_qubit == 2):
      Circuit712.q1GateList.append("O")
      Circuit712.q2GateList.append("H")
      Circuit712.q0GateList.append("-")

    elif (control_qubit == 1 and target_qubit == 0):
      Circuit712.q1GateList.append("O")
      Circuit712.q0GateList.append("H")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 2 and target_qubit == 0):
      Circuit712.q2GateList.append("O")
      Circuit712.q0GateList.append("H")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 2 and target_qubit == 1):
      Circuit712.q2GateList.append("O")
      Circuit712.q1GateList.append("H")
      Circuit712.q0GateList.append("-")

  def cnot(self, control_qubit: int, target_qubit: int):
    CNOTmatrix = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j ]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(control_qubit, int):
      raise TypeError("Control_qubit is not an :py:class:`int`")

    if not isinstance(target_qubit, int):
      raise TypeError("Target_qubit is not an :py:class:`int`")

    if (control_qubit < 0 or control_qubit >= self._number_of_qubits):
      raise ValueError ("Control_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (target_qubit < 0 or target_qubit >= self._number_of_qubits):
      raise ValueError ("Target_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (control_qubit == target_qubit):
      raise ValueError ("Target_qubit and Cotrol_qubit can not be IDENTICAL ")

    if (control_qubit == 0 and target_qubit == 1):
        CNOTmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j,0.+0.j, 1.+0.j], [0.+0.j, 0.+0.j, 1.+0.j,0.+0.j], [0.+0.j,1.+0.j,0.+0.j, 0.+0.j ]])
        CN0 = np.kron(Imatrix, CNOTmatrix1)
        Circuit712.q3bit = np.dot(CN0,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 2):
        CNOTmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j,0.+0.j, 1.+0.j], [0.+0.j, 0.+0.j, 1.+0.j,0.+0.j], [0.+0.j,1.+0.j,0.+0.j, 0.+0.j ]])
        CN0 = np.kron(CNOTmatrix1,Imatrix)
        Circuit712.q3bit = np.dot(CN0,Circuit712.q3bit)

    if (control_qubit == 0 and target_qubit == 2):
        CNOTmatrix3 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j,1.+0.j],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j ]])
        Circuit712.q3bit = np.dot(CNOTmatrix3,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 0):
        CN0 = np.kron(Imatrix, CNOTmatrix)
        Circuit712.q3bit = np.dot(CN0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 1):
        CN0 = np.kron(CNOTmatrix, Imatrix)
        Circuit712.q3bit = np.dot(CN0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 0):
        CNOTmatrix2 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j ],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,1.+0j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j ]])
        Circuit712.q3bit = np.dot(CNOTmatrix2,Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if (control_qubit == 0 and target_qubit == 1):
      Circuit712.q0GateList.append("O")
      Circuit712.q1GateList.append("X")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 0 and target_qubit == 2):
      Circuit712.q0GateList.append("O")
      Circuit712.q2GateList.append("X")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 1 and target_qubit == 2):
      Circuit712.q1GateList.append("O")
      Circuit712.q2GateList.append("X")
      Circuit712.q0GateList.append("-")

    elif (control_qubit == 1 and target_qubit == 0):
      Circuit712.q1GateList.append("O")
      Circuit712.q0GateList.append("X")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 2 and target_qubit == 0):
      Circuit712.q2GateList.append("O")
      Circuit712.q0GateList.append("X")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 2 and target_qubit == 1):
      Circuit712.q2GateList.append("O")
      Circuit712.q1GateList.append("X")
      Circuit712.q0GateList.append("-")

  def cs(self, control_qubit: int, target_qubit: int):
    CSmatrix = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j,], [0.+0.j,0.+0.j, 0.+0.j, 0.+1.j, ]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(control_qubit, int):
      raise TypeError("Control_qubit is not an :py:class:`int`")

    if not isinstance(target_qubit, int):
      raise TypeError("Target_qubit is not an :py:class:`int`")

    if (control_qubit < 0 or control_qubit >= self._number_of_qubits):
      raise ValueError ("Control_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (target_qubit < 0 or target_qubit >= self._number_of_qubits):
      raise ValueError ("Target_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (control_qubit == target_qubit):
      raise ValueError ("Target_qubit and Cotrol_qubit can not be IDENTICAL ")

    if (control_qubit == 0 and target_qubit == 1):
        CS0 = np.kron(Imatrix, CSmatrix)
        Circuit712.q3bit = np.dot(CS0,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 2):
        CS0 = np.kron(CSmatrix, Imatrix)
        Circuit712.q3bit = np.dot(CS0,Circuit712.q3bit)

    if (control_qubit == 0 and target_qubit == 2):
        CSmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 0.+1.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j,0.+1.j  ]])
        Circuit712.q3bit = np.dot(CSmatrix1,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 0):
        CS0 = np.kron(Imatrix, CSmatrix)
        Circuit712.q3bit = np.dot(CS0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 1):
        CS0 = np.kron(CSmatrix, Imatrix)
        Circuit712.q3bit = np.dot(CS0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 0):
        CSmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 0.+1.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j,0.+1.j  ]])
        Circuit712.q3bit = np.dot(CSmatrix1,Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if (control_qubit == 0 and target_qubit == 1):
      Circuit712.q0GateList.append("O")
      Circuit712.q1GateList.append("S")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 0 and target_qubit == 2):
      Circuit712.q0GateList.append("O")
      Circuit712.q2GateList.append("S")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 1 and target_qubit == 2):
      Circuit712.q1GateList.append("O")
      Circuit712.q2GateList.append("S")
      Circuit712.q0GateList.append("-")

    elif (control_qubit == 1 and target_qubit == 0):
      Circuit712.q1GateList.append("O")
      Circuit712.q0GateList.append("S")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 2 and target_qubit == 0):
      Circuit712.q2GateList.append("O")
      Circuit712.q0GateList.append("S")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 2 and target_qubit == 1):
      Circuit712.q2GateList.append("O")
      Circuit712.q1GateList.append("S")
      Circuit712.q0GateList.append("-")

  def cx(self, control_qubit: int, target_qubit: int):
    CXmatrix = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j ]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(control_qubit, int):
      raise TypeError("Control_qubit is not an :py:class:`int`")

    if not isinstance(target_qubit, int):
      raise TypeError("Target_qubit is not an :py:class:`int`")

    if (control_qubit < 0 or control_qubit >= self._number_of_qubits):
      raise ValueError ("Control_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (target_qubit < 0 or target_qubit >= self._number_of_qubits):
      raise ValueError ("Target_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (control_qubit == target_qubit):
      raise ValueError ("Target_qubit and Cotrol_qubit can not be IDENTICAL ")

    if (control_qubit == 0 and target_qubit == 1):
        CXmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j,0.+0.j, 1.+0.j], [0.+0.j, 0.+0.j, 1.+0.j,0.+0.j], [0.+0.j,1.+0.j,0.+0.j, 0.+0.j ]])
        CX0 = np.kron(Imatrix, CXmatrix1)
        Circuit712.q3bit = np.dot(CX0,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 2):
        CXmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j,0.+0.j, 1.+0.j], [0.+0.j, 0.+0.j, 1.+0.j,0.+0.j], [0.+0.j,1.+0.j,0.+0.j, 0.+0.j ]])
        CX0 = np.kron(CXmatrix1,Imatrix)
        Circuit712.q3bit = np.dot(CX0,Circuit712.q3bit)

    if (control_qubit == 0 and target_qubit == 2):
        CXmatrix3 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j,1.+0.j],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j ]])
        Circuit712.q3bit = np.dot(CXmatrix3,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 0):
        CX0 = np.kron(Imatrix, CXmatrix)
        Circuit712.q3bit = np.dot(CX0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 1):
        CX0 = np.kron(CXmatrix, Imatrix)
        Circuit712.q3bit = np.dot(CX0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 0):
        CXmatrix2 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j ],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,1.+0j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j ]])
        Circuit712.q3bit = np.dot(CXmatrix2,Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    if (control_qubit == 0 and target_qubit == 1):
      Circuit712.q0GateList.append("O")
      Circuit712.q1GateList.append("X")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 0 and target_qubit == 2):
      Circuit712.q0GateList.append("O")
      Circuit712.q2GateList.append("X")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 1 and target_qubit == 2):
      Circuit712.q1GateList.append("O")
      Circuit712.q2GateList.append("X")
      Circuit712.q0GateList.append("-")

    elif (control_qubit == 1 and target_qubit == 0):
      Circuit712.q1GateList.append("O")
      Circuit712.q0GateList.append("X")
      Circuit712.q2GateList.append("-")

    elif (control_qubit == 2 and target_qubit == 0):
      Circuit712.q2GateList.append("O")
      Circuit712.q0GateList.append("X")
      Circuit712.q1GateList.append("|")

    elif (control_qubit == 2 and target_qubit == 1):
      Circuit712.q2GateList.append("O")
      Circuit712.q1GateList.append("X")
      Circuit712.q0GateList.append("-")


  def swap(self, qubit1: int, qubit2: int):
    SWmatrix = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,1.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j,1.+0.j]])
    Imatrix = ([[1,0],[0,1]])

    if not isinstance(qubit1, int):
      raise TypeError("Origin qubit is not an :py:class:`int`")

    if not isinstance(qubit2, int):
      raise TypeError("Destination qubit is not an :py:class:`int`")

    if (qubit1 < 0 or qubit1 >= self._number_of_qubits):
      raise ValueError ("Origin qubit is less than zero or larger than the number of qubits in the circuit ")

    if (qubit2 < 0 or qubit2 >= self._number_of_qubits):
      raise ValueError ("Destination qubit is less than zero or larger than the number of qubits in the circuit ")

    if (qubit1 == qubit2):
      raise ValueError ("Origin qubit and Destination qubit can not be IDENTICAL ")

    if (qubit1 == 0 and qubit2 == 1):
        SW0 = np.kron(Imatrix, SWmatrix)
        Circuit712.q3bit = np.dot(SW0,Circuit712.q3bit)

    elif (qubit1 == 0 and qubit2 == 2):
        SWmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,1.+0.j, 0.+0.j],[ 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 1.+0.j ]])
        Circuit712.q3bit = np.dot(SWmatrix1,Circuit712.q3bit)

    elif (qubit1 == 1 and qubit2 == 2):
        SW0 = np.kron(SWmatrix, Imatrix)
        Circuit712.q3bit = np.dot(SW0,Circuit712.q3bit)

    elif (qubit1 == 1 and qubit2 == 0):
        SW0 = np.kron(Imatrix, SWmatrix)
        Circuit712.q3bit = np.dot(SW0,Circuit712.q3bit)

    elif (qubit1 == 2 and qubit2 == 0):
        SWmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,1.+0.j, 0.+0.j],[ 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 1.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, 1.+0.j ]])
        Circuit712.q3bit = np.dot(SWmatrix1,Circuit712.q3bit)

    elif (qubit1 == 2 and qubit2 == 1):
        SW0 = np.kron(SWmatrix, Imatrix)
        Circuit712.q3bit = np.dot(SW0,Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if (qubit1 == 0 and qubit2 == 1):
      Circuit712.q0GateList.append("X")
      Circuit712.q1GateList.append("X")
      Circuit712.q2GateList.append("-")

    elif (qubit1 == 0 and qubit2 == 2):
      Circuit712.q0GateList.append("X")
      Circuit712.q2GateList.append("X")
      Circuit712.q1GateList.append("|")

    elif (qubit1 == 1 and qubit2 == 2):
      Circuit712.q1GateList.append("X")
      Circuit712.q2GateList.append("X")
      Circuit712.q0GateList.append("-")

    elif (qubit1 == 1 and qubit2 == 0):
      Circuit712.q1GateList.append("X")
      Circuit712.q0GateList.append("X")
      Circuit712.q2GateList.append("-")

    elif (qubit1 == 2 and qubit2 == 0):
      Circuit712.q2GateList.append("X")
      Circuit712.q0GateList.append("X")
      Circuit712.q1GateList.append("|")

    elif (qubit1 == 2 and qubit2 == 1):
      Circuit712.q2GateList.append("X")
      Circuit712.q1GateList.append("X")
      Circuit712.q0GateList.append("-")


  def crz(self, control_qubit: int, target_qubit: int, phi: float):
    radian = (np.pi*phi/3.14)
    degree = int((np.pi/phi))
    CRZmatrix = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j,0.+0.j],[0.+0.j, 0.+0.j, cmath.exp(complex(0, -(radian))), 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, cmath.exp(complex(0, radian))]])
    Imatrix = ([[1,0],[0,1]])


    if not isinstance(control_qubit, int):
      raise TypeError("Control_qubit is not an :py:class:`int`")

    if not isinstance(target_qubit, int):
      raise TypeError("Target_qubit is not an :py:class:`int`")

    if (control_qubit < 0 or control_qubit >= self._number_of_qubits):
      raise ValueError ("Control_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (target_qubit < 0 or target_qubit >= self._number_of_qubits):
      raise ValueError ("Target_qubit is less than zero or larger than the number of qubits in the circuit ")

    if (control_qubit == target_qubit):
      raise ValueError ("Target_qubit and Cotrol_qubit can not be IDENTICAL ")

    if not isinstance(phi, float):
      raise TypeError("PHI is not an :py:class:`float`")

    if (control_qubit == 0 and target_qubit == 1):
        CRZmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, cmath.exp(complex(0, radian)), 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, cmath.exp(complex(0, radian)), 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, cmath.exp(complex(0, radian)), 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, cmath.exp(complex(0, radian)) ]])
        Circuit712.q3bit = np.dot(CRZmatrix1,Circuit712.q3bit)

    if (control_qubit == 0 and target_qubit == 2):
        CRZmatrix1 = np.array([[1.+0j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, cmath.exp(complex(0, radian)), 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j], [0.+0.j,0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, cmath.exp(complex(0, radian)), 0.+0.j, 0.+0.j, 0.+0.j,0.+0.j],[ 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, cmath.exp(complex(0, radian)), 0.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j], [0.+0.j,0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j,  0.+0.j, cmath.exp(complex(0, radian)) ]])
        Circuit712.q3bit = np.dot(CRZmatrix1,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 2):
        CRZ0 = np.kron(Imatrix, CRZmatrix)
        Circuit712.q3bit = np.dot(CRZ0,Circuit712.q3bit)

    if (control_qubit == 1 and target_qubit == 0):
        CRZ0 = np.kron(Imatrix, CRZmatrix)
        Circuit712.q3bit = np.dot(CRZ0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 1):
        CRZ0 = np.kron(CRZmatrix, Imatrix)
        Circuit712.q3bit = np.dot(CRZ0,Circuit712.q3bit)

    if (control_qubit == 2 and target_qubit == 0):
        CRZ0 = np.kron(CRZmatrix, Imatrix)
        Circuit712.q3bit = np.dot(CRZ0,Circuit712.q3bit)

    Circuit712.PMVector = Circuit712.q3bit

    if (control_qubit == 0 and target_qubit == 1):
      Circuit712.q0GateList.append("O")
      Circuit712.q0GateList.append("------")
      Circuit712.q1GateList.append("Z")
      Circuit712.q1GateList.append("R(π/%i)"%degree)
      Circuit712.q2GateList.append("--")
      Circuit712.q2GateList.append("-----")

    elif (control_qubit == 0 and target_qubit == 2):
      Circuit712.q0GateList.append("O")
      Circuit712.q0GateList.append("------")
      Circuit712.q2GateList.append("Z")
      Circuit712.q2GateList.append("R(π/%i)"%degree)
      Circuit712.q1GateList.append("|")
      Circuit712.q1GateList.append("------")

    elif (control_qubit == 1 and target_qubit == 2):
      Circuit712.q1GateList.append("O")
      Circuit712.q1GateList.append("------")
      Circuit712.q2GateList.append("Z")
      Circuit712.q2GateList.append("R(π/%i)"%degree)
      Circuit712.q0GateList.append("--")
      Circuit712.q0GateList.append("-----")

    elif (control_qubit == 1 and target_qubit == 0):
      Circuit712.q1GateList.append("O")
      Circuit712.q1GateList.append("------")
      Circuit712.q0GateList.append("Z")
      Circuit712.q0GateList.append("R(π/%i)"%degree)
      Circuit712.q2GateList.append("--")
      Circuit712.q2GateList.append("-----")

    elif (control_qubit == 2 and target_qubit == 0):
      Circuit712.q2GateList.append("O")
      Circuit712.q2GateList.append("------")
      Circuit712.q0GateList.append("Z")
      Circuit712.q0GateList.append("R(π/%i)"%degree)
      Circuit712.q1GateList.append("|")
      Circuit712.q1GateList.append("------")

    elif (control_qubit == 2 and target_qubit == 1):
      Circuit712.q2GateList.append("O")
      Circuit712.q2GateList.append("------")
      Circuit712.q1GateList.append("Z")
      Circuit712.q1GateList.append("R(π/%i)"%degree)
      Circuit712.q0GateList.append("--")
      Circuit712.q0GateList.append("-----")

#----------------------------------------------------------
# 3-qubit Quantum Gates
#----------------------------------------------------------
  def h_all (self):
    Hmatrix = np.array([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]])
    h_all = np.kron(Hmatrix, Hmatrix)
    h_all = np.kron(Hmatrix, h_all)
    Circuit712.q3bit = np.dot(h_all, Circuit712.q3bit)
    Circuit712.PMVector = Circuit712.q3bit

    Circuit712.q0GateList.append("-H")
    Circuit712.q1GateList.append("H")
    Circuit712.q2GateList.append("-H")

#----------------------------------------------------------
# Measurement and Simulation Quantum Gates
#----------------------------------------------------------

  def measure(self, qubit: int):
    if not isinstance(qubit, int):
      raise TypeError("Qubit is not an :py:class:`int`")

    if (qubit < 0 or qubit >= self._number_of_qubits):
      raise ValueError ("Qubit is less than zero or larger than the number of qubits in the circuit ")

    if qubit == 0:
      Circuit712.q0GateList.append("M(%i)"%Circuit712.MCounter)
      Circuit712.q1GateList.append("-")
      Circuit712.q2GateList.append("-")
      Circuit712.list_measurement.append(0)
      Circuit712.MCounter = Circuit712.MCounter+1

    if qubit == 1:
      Circuit712.q1GateList.append("M(%i)"%Circuit712.MCounter)
      Circuit712.q2GateList.append("-")
      Circuit712.q0GateList.append("-")
      Circuit712.list_measurement.append(1)
      Circuit712.MCounter = Circuit712.MCounter+1

    if qubit == 2:
      Circuit712.q2GateList.append("M(%i)"%Circuit712.MCounter)
      Circuit712.q1GateList.append("-")
      Circuit712.q0GateList.append("-")
      Circuit712.list_measurement.append(2)
      Circuit712.MCounter = Circuit712.MCounter+1


  def measure_all (self):
    Circuit712.q0GateList.append("M")
    Circuit712.q1GateList.append("M")
    Circuit712.q2GateList.append("M")

    # Collapsing the quantum state to a random state
    Circuit712.q3bit = [random.uniform(0, 1) for element in range(8)]

    Circuit712.list_measurement.append(0)
    Circuit712.list_measurement.append(1)
    Circuit712.list_measurement.append(2)


  def simulate (self, shots: int = 1):
    if not isinstance(shots, int):
      raise TypeError("The 'Shot Number' is not an :py:class:`int`")

    if (shots < 1):
      raise ValueError ("The shot Number' is smaller than 'one' ")

    for i in range(8):
        Circuit712.PMVector[i]= abs(Circuit712.PMVector[i] * Circuit712.PMVector[i])

    NVector = np.array([0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j])
    PMSVector =  np.array([0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j])

    Circuit712.counts = [0 for element in range(8)]
    for i in range(8):
        NVector[i] = Circuit712.PMVector[i][0]

    for i in range (8):
        PMSVector [i] = PMSVector [i-1] + NVector[i]

    for i in range (shots):
       RandomNumber = random.uniform(0, 1)
       if RandomNumber < PMSVector[0]:
            Circuit712.counts[0] = Circuit712.counts[0]+1
       if PMSVector[0] <= RandomNumber < PMSVector[1]:
            Circuit712.counts[1] = Circuit712.counts[1]+1
       if PMSVector[1] <= RandomNumber < PMSVector[2]:
            Circuit712.counts[2] = Circuit712.counts[2]+1
       if PMSVector[2] <= RandomNumber < PMSVector[3]:
            Circuit712.counts[3] = Circuit712.counts[3]+1
       if PMSVector[3] <= RandomNumber < PMSVector[4]:
            Circuit712.counts[4] = Circuit712.counts[4]+1
       if PMSVector[4] <= RandomNumber < PMSVector[5]:
            Circuit712.counts[5] = Circuit712.counts[5]+1
       if PMSVector[5] <= RandomNumber < PMSVector[6]:
            Circuit712.counts[6] = Circuit712.counts[6]+1
       if PMSVector[6] <= RandomNumber < PMSVector[7]:
            Circuit712.counts[7] = Circuit712.counts[7]+1

    if (Circuit712.list_measurement == [2]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[1] + Circuit712.counts[2] + Circuit712.counts[3]
        Circuit712.counts[1] = Circuit712.counts[4]+ Circuit712.counts[5] + Circuit712.counts[6] + Circuit712.counts[7]
        print("{'0': %d, '1': %d }" %(Circuit712.counts[0],Circuit712.counts[1]))

    if (Circuit712.list_measurement == [1]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[1] + Circuit712.counts[4] + Circuit712.counts[5]
        Circuit712.counts[1] = Circuit712.counts[2]+ Circuit712.counts[3] + Circuit712.counts[6] + Circuit712.counts[7]
        print("{'0': %d, '1': %d }" %(Circuit712.counts[0],Circuit712.counts[1]))

    if (Circuit712.list_measurement == [0]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[2] + Circuit712.counts[4] + Circuit712.counts[6]
        Circuit712.counts[1] = Circuit712.counts[1]+ Circuit712.counts[3] + Circuit712.counts[5] + Circuit712.counts[7]
        print("{'0': %d, '1': %d }" %(Circuit712.counts[0],Circuit712.counts[1]))

    if (Circuit712.list_measurement == [0,1]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[4]
        Circuit712.counts[1] = Circuit712.counts[1]+ Circuit712.counts[5]
        Circuit712.counts[2] = Circuit712.counts[2]+ Circuit712.counts[6]
        Circuit712.counts[3] = Circuit712.counts[3]+ Circuit712.counts[7]
        print("{'00': %d, '01': %d, '10': %d, '11': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[2], Circuit712.counts[3]))

    if (Circuit712.list_measurement == [1,0]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[4]
        Circuit712.counts[4] = Circuit712.counts[2]+ Circuit712.counts[6]
        Circuit712.counts[2] = Circuit712.counts[1]+ Circuit712.counts[5]
        Circuit712.counts[3] = Circuit712.counts[3]+ Circuit712.counts[7]
        print("{'00': %d, '01': %d, '10': %d, '11': %d }" %(Circuit712.counts[0],Circuit712.counts[4], Circuit712.counts[2], Circuit712.counts[3]))

    if (Circuit712.list_measurement == [0,2]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[2]
        Circuit712.counts[1] = Circuit712.counts[1]+ Circuit712.counts[3]
        Circuit712.counts[2] = Circuit712.counts[4]+ Circuit712.counts[6]
        Circuit712.counts[3] = Circuit712.counts[5]+ Circuit712.counts[7]
        print("{'00': %d, '01': %d, '10': %d, '11': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[2], Circuit712.counts[3]))

    if (Circuit712.list_measurement == [2,0]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[2]
        Circuit712.counts[2] = Circuit712.counts[1]+ Circuit712.counts[3]
        Circuit712.counts[1] = Circuit712.counts[4]+ Circuit712.counts[6]
        Circuit712.counts[3] = Circuit712.counts[5]+ Circuit712.counts[7]
        print("{'00': %d, '01': %d, '10': %d, '11': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[2], Circuit712.counts[3]))

    if (Circuit712.list_measurement == [1,2]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[1]
        Circuit712.counts[1] = Circuit712.counts[2]+ Circuit712.counts[3]
        Circuit712.counts[2] = Circuit712.counts[4]+ Circuit712.counts[5]
        Circuit712.counts[3] = Circuit712.counts[6]+ Circuit712.counts[7]
        print("{'00': %d, '01': %d, '10': %d, '11': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[2], Circuit712.counts[3]))

    if (Circuit712.list_measurement == [2,1]):
        Circuit712.counts[0] = Circuit712.counts[0]+ Circuit712.counts[1]
        Circuit712.counts[1] = Circuit712.counts[4]+ Circuit712.counts[5]
        Circuit712.counts[2] = Circuit712.counts[2]+ Circuit712.counts[3]
        Circuit712.counts[3] = Circuit712.counts[6]+ Circuit712.counts[7]
        print("{'00': %d, '01': %d, '10': %d, '11': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[2], Circuit712.counts[3]))

    if (Circuit712.list_measurement == [0,1,2]):
        print("{'000': %d, '001': %d, '010': %d, '011': %d, '100': %d, '101': %d, '110': %d, '111': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[2], Circuit712.counts[3], Circuit712.counts[4], Circuit712.counts[5], Circuit712.counts[6], Circuit712.counts[7]))

    if (Circuit712.list_measurement == [2,1,0]):
        print("{'000': %d, '001': %d, '010': %d, '011': %d, '100': %d, '101': %d, '110': %d, '111': %d }" %(Circuit712.counts[0],Circuit712.counts[4], Circuit712.counts[2], Circuit712.counts[6], Circuit712.counts[1], Circuit712.counts[5], Circuit712.counts[3], Circuit712.counts[7]))

    if (Circuit712.list_measurement == [1,2,0]):
        print("{'000': %d, '001': %d, '010': %d, '011': %d, '100': %d, '101': %d, '110': %d, '111': %d }" %(Circuit712.counts[0],Circuit712.counts[2], Circuit712.counts[4], Circuit712.counts[6], Circuit712.counts[1], Circuit712.counts[3], Circuit712.counts[5], Circuit712.counts[7]))

    if (Circuit712.list_measurement == [1,0,2]):
        print("{'000': %d, '001': %d, '010': %d, '011': %d, '100': %d, '101': %d, '110': %d, '111': %d }" %(Circuit712.counts[0],Circuit712.counts[2], Circuit712.counts[1], Circuit712.counts[3], Circuit712.counts[4], Circuit712.counts[6], Circuit712.counts[5], Circuit712.counts[7]))

    if (Circuit712.list_measurement == [2,0,1]):
        print("{'000': %d, '001': %d, '010': %d, '011': %d, '100': %d, '101': %d, '110': %d, '111': %d }" %(Circuit712.counts[0],Circuit712.counts[4], Circuit712.counts[1], Circuit712.counts[5], Circuit712.counts[2], Circuit712.counts[6], Circuit712.counts[3], Circuit712.counts[7]))

    if (Circuit712.list_measurement == [0,2,1]):
        print("{'000': %d, '001': %d, '010': %d, '011': %d, '100': %d, '101': %d, '110': %d, '111': %d }" %(Circuit712.counts[0],Circuit712.counts[1], Circuit712.counts[4], Circuit712.counts[5], Circuit712.counts[2], Circuit712.counts[3], Circuit712.counts[6], Circuit712.counts[7]))


  def simulation_pre_measurement_statevector(self):
      print("\n The Pre-Measurement Statevector:\n ",Circuit712.PMVector)


#----------------------------------------------------------
# Showing Functions
#----------------------------------------------------------


  def create_middle_line(self, top_list: list, bottom_list: list):
    source_target_indicator = ['O', 'X']
    source_target_indicator1 = ['O', 'S']
    source_target_indicator2 = ['O', 'Z']
    source_target_indicator3 = ['O', 'H']
    results = []

    for index in range(len(top_list)):
      top = top_list[index]
      bottom = bottom_list[index]

      if top == '|' or bottom == '|' or (top in source_target_indicator and bottom in source_target_indicator) or (top in source_target_indicator2 and bottom in source_target_indicator2) or (top in source_target_indicator1 and bottom in source_target_indicator1) or (top in source_target_indicator3 and bottom in source_target_indicator3):
        results.append('|')
      else:
        results.append(' ' * len(top))
    return results

  def draw_circuit(self):
    ar0 = np.array(Circuit712.q0GateList)
    ar1 = np.array(Circuit712.q1GateList)
    ar2 = np.array(Circuit712.q2GateList)
    array2 = np.stack((ar0, ar1), axis=0)
    array3 = np.r_[array2,[ar2]]
    item_decorator = '{}--'
    space_decorator = item_decorator.replace('-', ' ')
    input_count = len(array3)

    for index in range(len(array3)):
      input = array3[index]

      # Display input data
      formated_input = [item_decorator.format(item) for item in input]
      formated_input = ''.join(formated_input)
      print(formated_input)

      # Display the middle line if needed
      if input_count > index+1:
        middle_line = self.create_middle_line(input, array3[index+1])
        formated_space = [space_decorator.format(item) for item in middle_line]
        formated_space = ''.join(formated_space)
        print(formated_space)

  def draw_histogram (self, color: str = "blue"):
    if not isinstance(color, str):
      raise TypeError("The 'Color' is not an :py:class:`str")

    if (Circuit712.list_measurement == [1] or Circuit712.list_measurement == [2] or Circuit712.list_measurement == [0]):
        data = {'0':Circuit712.counts[0], '1':Circuit712.counts[1]}

    elif (Circuit712.list_measurement == [0,1]):
        data = {'00':Circuit712.counts[0], '01':Circuit712.counts[1], '10':Circuit712.counts[2], '11':Circuit712.counts[3]}

    elif (Circuit712.list_measurement == [1,0]):
        data = {'00':Circuit712.counts[0], '01':Circuit712.counts[4], '10':Circuit712.counts[2], '11':Circuit712.counts[3]}

    elif (Circuit712.list_measurement == [0,2]):
        data = {'00':Circuit712.counts[0], '01':Circuit712.counts[1], '10':Circuit712.counts[2], '11':Circuit712.counts[3]}

    elif (Circuit712.list_measurement == [2,0]):
        data = {'00':Circuit712.counts[0], '01':Circuit712.counts[1], '10':Circuit712.counts[2], '11':Circuit712.counts[3]}

    elif (Circuit712.list_measurement == [1,2]):
        data = {'00':Circuit712.counts[0], '01':Circuit712.counts[1], '10':Circuit712.counts[2], '11':Circuit712.counts[3]}

    elif (Circuit712.list_measurement == [2,1]):
        data = {'00':Circuit712.counts[0], '01':Circuit712.counts[1], '10':Circuit712.counts[2], '11':Circuit712.counts[3]}

    elif (Circuit712.list_measurement == [2,1,0]):
        data = {'000':Circuit712.counts[0], '001':Circuit712.counts[4], '010':Circuit712.counts[2], '011':Circuit712.counts[6], '100':Circuit712.counts[1], '101':Circuit712.counts[5], '110':Circuit712.counts[3], '111':Circuit712.counts[7]}

    elif (Circuit712.list_measurement == [1,2,0]):
        data = {'000':Circuit712.counts[0], '001':Circuit712.counts[2], '010':Circuit712.counts[4], '011':Circuit712.counts[6], '100':Circuit712.counts[1], '101':Circuit712.counts[3], '110':Circuit712.counts[5], '111':Circuit712.counts[7]}

    elif (Circuit712.list_measurement == [1,0,2]):
        data = {'000':Circuit712.counts[0], '001':Circuit712.counts[2], '010':Circuit712.counts[1], '011':Circuit712.counts[3], '100':Circuit712.counts[4], '101':Circuit712.counts[6], '110':Circuit712.counts[5], '111':Circuit712.counts[7]}

    else:
        data = {'000':Circuit712.counts[0], '001':Circuit712.counts[1], '010':Circuit712.counts[2], '011':Circuit712.counts[3], '100':Circuit712.counts[4], '101':Circuit712.counts[5], '110':Circuit712.counts[6], '111':Circuit712.counts[7]}

    values = list(data.keys())
    counts = list(data.values())
    res = [idx for idx, val in enumerate(counts) if val != 0]
    values1 = [0] * len(res)
    counts1 = [0] * len(res)
    for i in range(len(res)):
      values1[i] = values[res[i]]
      counts1[i] = counts[res[i]]

    fig = plt.figure(figsize = (7, 5))
    plt.bar(values1, counts1, color = color, width = 0.7)

    for i in range(len(res)):
        plt.text(i,counts1[i],counts1[i], ha = 'center')

    plt.ylabel("Counts")
    plt.title("Simulation Results")
    plt.show()