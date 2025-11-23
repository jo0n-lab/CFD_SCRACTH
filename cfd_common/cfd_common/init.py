# cfd_tdma/solver.py
import numpy as np

# 전역 변수들: 기존 코드 스타일 유지
IDIM = None
JDIM = None
IJDIM = None
NI = NJ = NIM = NJM = None
DX0 = DY0 = None
IR = 0

# 도메인 및 필드 배열들 (init()에서 생성)
XP = XU = None
DXP = DXU = None
YP = YV = None
DYP = DYV = None
RP = RU = None

T = TSOR = BT = None
U = V = USOR = VSOR = BU = BV = None
TKKP = ROCP = None
AW = AE = AN = AS = AP = APC = None
CORR = RES = None


def init():
    global IDIM, JDIM, IJDIM
    global NI, NIM, NJ, NJM
    global XP, XU, YP, YV
    global DXP, DXU, DYP, DYV
    global RP, RU
    global T, TSOR, BT
    global U, USOR, BU
    global V, VSOR, BV
    global TKKP, ROCP
    global AW, AE, AN, AS, AP, APC
    global CORR, RES

    # DOMAIN ARRAYS
    XP  = np.zeros(IDIM)
    XU  = np.zeros(IDIM)
    
    # CONTROL VOLUME SIZE ARRAYS
    DXP = np.zeros(IDIM) 
    DXU = np.zeros(IDIM)
    
    YP  = np.zeros(JDIM)
    YV  = np.zeros(JDIM)
    DYP = np.zeros(JDIM)
    DYV = np.zeros(JDIM)
    
    RP  = np.zeros(IDIM)
    RU  = np.zeros(IDIM)
    
    NIM = NI - 1
    NJM = NJ - 1
    IJDIM = max(IDIM, JDIM)

    # initialize domain sizes
    for i in range(1, NI + 1): XU[i] = DX0 * (i - 1)
    for i in range(1, NI):     XP[i] = (XU[i] + XU[i + 1]) / 2.0
    for j in range(1, NJ + 1): YV[j] = DY0 * (j - 1)
    for j in range(1, NJ):     YP[j] = (YV[j] + YV[j + 1]) / 2.0

    XP[0] = XU[1]
    XP[NI] = XU[NI]
    YP[0] = YV[1]
    YP[NJ] = YV[NJ]

    # initialize control volume sizes
    for i in range(1, NI):     DXP[i] = XU[i + 1] - XU[i]
    for i in range(1, NI + 1): DXU[i] = XP[i] - XP[i - 1]
    for j in range(1, NJ):     DYP[j] = YV[j + 1] - YV[j]
    for j in range(1, NJ + 1): DYV[j] = YP[j] - YP[j - 1]

    # initialize face dimensions
    if IR == 0:
        for i in range(1, NI + 1):
            RU[i] = 1.0
        for i in range(0, NI + 1):
            RP[i] = 1.0
    else:
        for i in range(1, NI + 1):
            RU[i] = np.fabs(XU[i])
        for i in range(0, NI + 1):
            RP[i] = np.fabs(XP[i])
