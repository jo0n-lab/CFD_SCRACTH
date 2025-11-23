from .init import init, IDIM, JDIM, NI, NJ, NIM, NJM, DX0, DY0, IR
from .init import XP, XU, YP, YV, DXP, DXU, DYP, DYV, RP, RU
from .init import T, TSOR, BT, U, USOR, BU, V, VSOR, BV
from .init import TKKP, ROCP, AW, AE, AN, AS, AP, APC, CORR, RES

# TDMA 관련 함수
from .tdma import TDMA, LLTDMA

__all__ = [
    "init",
    "IDIM", "JDIM", "NI", "NJ", "NIM", "NJM", "DX0", "DY0", "IR",
    "XP", "XU", "YP", "YV", "DXP", "DXU", "DYP", "DYV", "RP", "RU",
    "T", "TSOR", "BT", "U", "USOR", "BU", "V", "VSOR", "BV",
    "TKKP", "ROCP", "AW", "AE", "AN", "AS", "AP", "APC", "CORR", "RES",
    "TDMA", "LLTDMA",
]


@property
def Nx(self): return self.NI
@property
def Ny(self): return self.NJ