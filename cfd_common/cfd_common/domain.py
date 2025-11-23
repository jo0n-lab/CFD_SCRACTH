from dataclasses import dataclass

@dataclass
class DomainInfo:
    NI: int
    NJ: int
    IDIM: int
    JDIM: int
    IJDIM: int
    NIM: int
    NJM: int
    DX0: float
    DY0: float

# 안전성 확보: 초기값 None
domain: DomainInfo | None = None


def set_domain(NI, NJ, DX0, DY0):
    global domain

    IDIM = NI + 2
    JDIM = NJ + 2
    IJDIM = max(IDIM, JDIM)
    NIM   = NI - 1
    NJM   = NJ - 1

    domain = DomainInfo(
        NI=NI, NJ=NJ,
        IDIM=IDIM, JDIM=JDIM,
        IJDIM=IJDIM, NIM=NIM, NJM=NJM,
        DX0=DX0, DY0=DY0
    )
