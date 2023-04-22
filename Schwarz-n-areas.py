import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import diags, linalg, coo_matrix
import warnings
import time
import threading
import math
from multiprocessing import shared_memory, Process, Pipe

warnings.filterwarnings("ignore")

alpha = 1
beta = 1
# u_t = alpha*u_xx+beta*u*(1-u)


def fisher(x, t): return (1+np.exp((beta/(6*alpha))**0.5*x-5*beta*t/6))**(-2)
def FFF(u): return beta*u*(1-u)
def u_x_ta(x, t=0): return (1+np.exp((beta/(6*alpha))**0.5*x-5*beta*t/6))**(-2)
def u_xa_t(t, x=0): return (1+np.exp((beta/(6*alpha))**0.5*x-5*beta*t/6))**(-2)
def u_xb_t(t, x=1): return (1+np.exp((beta/(6*alpha))**0.5*x-5*beta*t/6))**(-2)

# def u_truth(x, t): return x * (x ** 2 + np.exp(t))
# def q_fn(x, t): return x * np.exp(t) - 6 * x
# def u_x_ta(x): return x ** 3 + x
# def u_xa_t(t): return np.zeros_like(t)
# def u_xb_t(t): return 1 + np.exp(t)


class subPro():
    def __init__(self, pipe, Nx=100, Nt=1000, xa=0, xb=1, ta=0, tb=1):
        self.Nx = Nx
        self.Nt = Nt
        self.xa = xa
        self.xb = xb
        self.ta = ta
        self.tb = tb
        self.h = (self.xb-self.xa)/Nx
        self.tau = (self.tb-self.ta)/Nt
        self.r = alpha*self.tau/self.h**2
        self.pipe = pipe

    def run(self, solName, AName, BName):
        t1 = time.time()
        solshm = shared_memory.SharedMemory(name=solName)
        Ashm = shared_memory.SharedMemory(name=AName)
        Bshm = shared_memory.SharedMemory(name=BName)

        sol = np.ndarray(shape=(self.Nx+1, self.Nt+1),
                         dtype=np.float64, buffer=solshm.buf)
        A = coo_matrix(np.ndarray(shape=(self.Nx-1, self.Nx-1),
                       dtype=np.float64, buffer=Ashm.buf))
        B = coo_matrix(np.ndarray(shape=(self.Nx-1, self.Nx-1),
                       dtype=np.float64, buffer=Bshm.buf))

        self.xx = np.linspace(self.xa, self.xb, self.Nx + 1, endpoint=True)
        self.tt = np.linspace(self.ta, self.tb, self.Nt + 1, endpoint=True)

        # print(solName,self.xa,self.xb,t2-t1)
        while(1):
            tmpstr=self.pipe.recv()
            if tmpstr=="break":
                break
            t2 = time.time()
            for t_idx in range(self.Nt):
                _f = sol[1:-1, t_idx]
                cc = np.zeros(self.Nx - 1)
                cc[0] += self.r / 2 * (sol[0, t_idx + 1] + sol[0, t_idx])
                cc[-1] += self.r / 2 * (sol[-1, t_idx + 1] + sol[-1, t_idx])
                dd = B * _f + self.tau / 2 * FFF(u=_f) + cc
                __f = _f.copy()
                for iii in range(2):
                    delta = linalg.inv(A - diags(self.tau / 2 * (1 - 2 * __f))) * \
                        (A * __f - self.tau / 2 * FFF(u=__f) - dd)
                    __f = __f - delta
                sol[1:-1, t_idx+1] = __f
            self.pipe.send("complete")

            print(solName,self.xa,self.xb,time.time()-t2)


class subProMatrix():
    def __init__(self, pipe, Nx=100, Nt=1000, xa=0, xb=1, ta=0, tb=1):
        self.Nx = Nx
        self.Nt = Nt
        self.xa = xa
        self.xb = xb
        self.ta = ta
        self.tb = tb
        self.h = (self.xb-self.xa)/Nx
        self.tau = (self.tb-self.ta)/Nt
        self.r = alpha*self.tau/self.h**2
        self.pipe = pipe
        self.xx = np.linspace(self.xa, self.xb, Nx + 1, endpoint=True)
        self.tt = np.linspace(self.ta, self.tb, Nt + 1, endpoint=True)
        A = diags([[-self.r / 2] * (Nx - 2), [1 + self.r] * (Nx - 1),
                   [-self.r / 2] * (Nx - 2)], [-1, 0, 1]).toarray()
        B = diags([[self.r / 2] * (Nx - 2), [1 - self.r] * (Nx - 1),
                   [self.r / 2] * (Nx - 2)], [-1, 0, 1]).toarray()
        xx_mesh, tt_mesh = np.meshgrid(self.xx, self.tt)
        self.gg = fisher(x=xx_mesh, t=tt_mesh).T
        # print(self.xx[-1],self.tt,self.gg[-1,:])

        self.solshm = shared_memory.SharedMemory(
            create=True, size=(self.Nt+1)*(self.Nx+1)*8+10000)
        self.Ashm = shared_memory.SharedMemory(
            create=True, size=A.size*8+10000)
        self.Bshm = shared_memory.SharedMemory(
            create=True, size=B.size*8+10000)
        self.sol = np.ndarray(shape=(Nx+1, Nt+1),
                              dtype=np.float64, buffer=self.solshm.buf)
        self.A = np.ndarray(shape=A.shape, dtype=A.dtype, buffer=self.Ashm.buf)
        self.B = np.ndarray(shape=B.shape, dtype=B.dtype, buffer=self.Bshm.buf)
        self.sol[:, 0] = self.gg[:, 0]
        for i in range(A.shape[0]):
            self.A[i, :] = A[i, :]
            self.B[i, :] = B[i, :]


if __name__ == '__main__':
    TA = 0
    TB = 25
    NT = 250
    pipes = [Pipe() for i in range(4)]
    proArgList = [[pipes[i][1], 20, NT, i, i+2, TA, TB]for i in range(len(pipes))]
    # print(proArgList)
    # proArgList = [[20, NT, 0, 2, TA, TB], [20, NT, 1, 3, TA, TB], [20, NT, 2, 4, TA, TB], 
    #               [20, NT, 3, 5, TA, TB], [20, NT, 4, 6, TA, TB], [20, NT, 5, 7, TA, TB],
    #               [20, NT, 6, 8, TA, TB], [20, NT, 7, 9, TA, TB], [20, NT, 8, 10, TA, TB]]
    # proArgList = [[30, NT, 0, 3, TA, TB], [30, NT, 2, 5, TA, TB], [30, NT, 4, 7, TA, TB],[30, NT, 6, 9, TA, TB], [20, NT, 8, 10, TA, TB]]
    # proArgList = [[40, NT, 0, 4, 0, 10], [40, NT, 3, 7, 0, 10], [40, NT, 6, 10, 0, 10]]
    # proArgList=[[70,NT,0,7,0,10],[70,NT,3,10,0,10]]
    XN = len(proArgList)
    subPros = []
    subProMatrixs = []
    for i in range(XN):
        ttt = proArgList[i]
        subPros.append(subPro(ttt[0], ttt[1], ttt[2], ttt[3], ttt[4], ttt[5], ttt[6]))
        subProMatrixs.append(subProMatrix(ttt[0], ttt[1], ttt[2], ttt[3], ttt[4], ttt[5], ttt[6]))

    for xn in range(XN):
        subProo = subProMatrixs[xn]
        if xn == 0:
            subProo.sol[0, :] = u_xa_t(t=subProo.tt, x=subProo.xa)
            subProo.sol[-1, 1:] = 0
        elif xn == XN-1:
            subProo.sol[0, 1:] = 0
            subProo.sol[-1, :] = u_xb_t(t=subProo.tt, x=subProo.xb)
        else:
            subProo.sol[0, 1:] = 0
            subProo.sol[-1, 1:] = 0
    processList = list(range(XN))
    for xn in range(XN):
        tmpMatrix = subProMatrixs[xn]
        processList[xn] = Process(target=subPros[xn].run, args=(
            tmpMatrix.solshm.name, tmpMatrix.Ashm.name, tmpMatrix.Bshm.name))
        processList[xn].start()
    for i in range(100):
        t1 = time.time()
        for xn in range(XN):
            pipes[xn][0].send("{} start".format(xn))
        for xn in range(XN):
            pipes[xn][0].recv()
        tmperr = []
        for xn in range(XN):
            tmperr.append(
                np.max(np.abs(subProMatrixs[xn].gg - subProMatrixs[xn].sol)))
            print(f"{i} r{xn}={subProMatrixs[xn].r:.4f} err:{tmperr[-1]:.6e}")
        print(time.time()-t1,"max_err:{}".format(max(tmperr)))
        if(max(tmperr)<4e-4):
            for xn in range(XN):
                pipes[xn][0].send("break")
            break

        # ssubpro=subProMatrixs[0]
        # xx=np.linspace(ssubpro.xa, ssubpro.xb, ssubpro.Nx + 1, endpoint=True)
        # tt=np.linspace(ssubpro.ta, ssubpro.tb, ssubpro.Nt + 1, endpoint=True)
        # xx_mesh, tt_mesh = np.meshgrid(xx,tt)
        # fig = plt.figure(figsize=(15, 5))
        # fig.set_tight_layout(True)
        # ax = fig.add_subplot(131, projection='3d')
        # for kllk in range(ssubpro.gg.T.shape[0]):
        #     ax.plot3D(xx_mesh[kllk],tt_mesh[kllk],ssubpro.gg.T[kllk,:])
        # # ax.plot_surface(xx_mesh, tt_mesh, ssubpro.gg.T, cmap="rainbow")
        # ax = fig.add_subplot(132, projection='3d')
        # # ax.plot_surface(xx_mesh, tt_mesh, ssubpro.sol.T, cmap="rainbow")
        # for kllk in range(ssubpro.sol.T.shape[0]):
        #     ax.plot3D(xx_mesh[kllk],tt_mesh[kllk],ssubpro.sol.T[kllk,:])
        # ax = fig.add_subplot(133, projection='3d')
        # ax.plot_surface(xx_mesh, tt_mesh, ssubpro.gg.T-ssubpro.sol.T, cmap="rainbow")
        # plt.show()

        for xn in range(XN-1):
            curPro = subProMatrixs[xn]
            nextPro = subProMatrixs[xn+1]
            if curPro.Nt > nextPro.Nt:
                curPro.sol[-1, 1:] = nextPro.sol[round(
                    (curPro.xb-nextPro.xa)/nextPro.h), 1:].repeat(curPro.Nt/nextPro.Nt)
                nextPro.sol[0, 1:] = curPro.sol[-1 -
                                                round((curPro.xb-nextPro.xa)/curPro.h), 1::int(curPro.Nt/nextPro.Nt)]
            else:
                curPro.sol[-1, 1:] = nextPro.sol[round(
                    (curPro.xb-nextPro.xa)/nextPro.h), 1::int(nextPro.Nt/curPro.Nt)]
                nextPro.sol[0, 1:] = curPro.sol[-1 -
                                                round((curPro.xb-nextPro.xa)/curPro.h), 1:].repeat(nextPro.Nt/curPro.Nt)
