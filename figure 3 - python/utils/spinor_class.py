import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from matplotlib import cm
import sys
from os import path
from os import makedirs


class Spinor:
  def __init__(self, time_resolution=2**16):
      # 1) initialize
      self.t = np.linspace(0, 2 * np.pi, time_resolution)                               # time vector
      self.I, self.li, self.lj, self.lk, self.ri, self.rj, self.rk = self.cayley()      # cayley matrices
      # 2) unitary
      self.U_4vector, self.U_4vector_dt, self.U_4vector_dt2 = self.unitary()
      # 3) S2 sphere
      self.surface = S2sphere()
      # 3) Hamiltonian
      self.H, self.H_dt = self.get_hamiltonian()
      # 4) Pre calculations for phase grids
      self.phase_pre_calcs = self.phase_grids_preamble()


  def get_hamiltonian(self):
      # quaternion shorthand
      a = self.U_4vector[:, 0]
      b = self.U_4vector[:, 1]
      c = self.U_4vector[:, 2]
      d = self.U_4vector[:, 3]
      # gradients
      a_dt = self.U_4vector_dt[:, 0]
      b_dt = self.U_4vector_dt[:, 1]
      c_dt = self.U_4vector_dt[:, 2]
      d_dt = self.U_4vector_dt[:, 3]
      # gradients - second order
      a_dt2 = self.U_4vector_dt2[:, 0]
      b_dt2 = self.U_4vector_dt2[:, 1]
      c_dt2 = self.U_4vector_dt2[:, 2]
      d_dt2 = self.U_4vector_dt2[:, 3]
      # electric and magentic
      Ei, Bi = 2 * (a * b_dt - a_dt * b), 2 * (c * d_dt - c_dt * d)
      Ej, Bj = 2 * (a * c_dt - a_dt * c), 2 * (b_dt * d - b * d_dt)
      Ek, Bk = 2 * (a * d_dt - a_dt * d), 2 * (b * c_dt - b_dt * c)
      # electric and magentic
      Ei_dt, Bi_dt = 2 * (a * b_dt2 - a_dt2 * b), 2 * (c * d_dt2 - c_dt2 * d)
      Ej_dt, Bj_dt = 2 * (a * c_dt2 - a_dt2 * c), 2 * (b_dt2 * d - b * d_dt2)
      Ek_dt, Bk_dt = 2 * (a * d_dt2 - a_dt2 * d), 2 * (b * c_dt2 - b_dt2 * c)
      # hamiltonian
      H = np.array([Bi + Ei, Bj + Ej, Bk + Ek]).T
      H_dt = np.array([Bi_dt + Ei_dt, Bj_dt + Ej_dt, Bk_dt + Ek_dt]).T
      return H, H_dt


  def phase_grids_preamble(self, N_t=2 ** 12):
      t = np.linspace(0, 2 * np.pi, N_t)
      dt = t[1] - t[0]
      # hamiltonian shorthand
      H = np.zeros((N_t, 3))
      H[:, 0] = np.interp(t, self.t, self.H[:, 0])
      H[:, 1] = np.interp(t, self.t, self.H[:, 1])
      H[:, 2] = np.interp(t, self.t, self.H[:, 2])
      Hi, Hj, Hk = H[:, 0].reshape(N_t, 1), H[:, 1].reshape(N_t, 1), H[:, 2].reshape(N_t, 1)
      # hamiltonian_dt shorthand
      H_dt = np.zeros((N_t, 3))
      H_dt[:, 0] = np.interp(t, self.t, self.H_dt[:, 0])
      H_dt[:, 1] = np.interp(t, self.t, self.H_dt[:, 1])
      H_dt[:, 2] = np.interp(t, self.t, self.H_dt[:, 2])
      Hi_dt, Hj_dt, Hk_dt = H_dt[:, 0].reshape(N_t, 1), H_dt[:, 1].reshape(N_t, 1), H_dt[:, 2].reshape(N_t, 1)
      # quaternion shorthand
      a = np.interp(t, self.t, self.U_4vector[:, 0])
      b = np.interp(t, self.t, self.U_4vector[:, 1])
      c = np.interp(t, self.t, self.U_4vector[:, 2])
      d = np.interp(t, self.t, self.U_4vector[:, 3])
      # SO(3) rotor
      rotor = self.quaternion2rotor(np.array([a, b, c, d]).T)
      return rotor, H, Hi, Hj, Hk, H_dt, Hi_dt, Hj_dt, Hk_dt, dt


  def phase_grids(self, principle_axis=None):
      dim = (2 ** 6 + 2 ** 5 + 2 ** 4, 2 ** 7)  # phase grid dimensions
      # extract
      h = self.phase_pre_calcs
      rotor, H, Hi, Hj, Hk, H_dt, Hi_dt, Hj_dt, Hk_dt, dt = h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9]
      # GLOBAL PHASE GRID:
      ν0, μ0 = np.meshgrid(np.linspace(0, 2 * np.pi, dim[1]), np.linspace(0, np.pi, dim[0]))
      ν0, μ0 = ν0.flatten(), μ0.flatten()
      # bloch vector
      np.seterr(divide='ignore', invalid='ignore')  # suppress divide by 0 warning
      if principle_axis == 'i':  # Ri, Rj, Rk = np.cos(μ0), np.sin(μ0) * np.sin(ν0), np.sin(μ0) * np.cos(ν0)
          Ri, Rj, Rk = self.evolve_R(μ0, ν0, rotor, principle_axis)
          # DYNAMIC PHASE:
          dξdt = Hi * Ri + Hj * Rj + Hk * Rk
          # TANGENT FRAME: Geometric phase σi
          dγdt = - Hi * Ri + (Hj * Rj + Hk * Rk) * Ri ** 2 / (Rj ** 2 + Rk ** 2)
          # DARBOUX FRAME: Geometric phase σi
          denom = (Hi * Hi + Hj * Hj + Hk * Hk - dξdt ** 2)
          H_cross_R_i = Hj * Rk - Hk * Rj
          H_cross_R_j = Hk * Ri - Hi * Rk
          H_cross_R_k = Hi * Rj - Hj * Ri
          Hdt_Rdt = Hi_dt * H_cross_R_i + Hj_dt * H_cross_R_j + Hk_dt * H_cross_R_k
          Rdt_Rdt = H_cross_R_i * H_cross_R_i + H_cross_R_j * H_cross_R_j + H_cross_R_k * H_cross_R_k
          DARBOUX_dγdt = (Hdt_Rdt - Rdt_Rdt * dξdt) / denom
      elif principle_axis == 'j':  # Ri, Rj, Rk = np.sin(μ0) * np.sin(ν0), np.cos(μ0), np.sin(μ0) * np.cos(ν0)
          Ri, Rj, Rk = self.evolve_R(μ0, ν0, rotor, principle_axis)
          # DYNAMIC PHASE:
          dξdt = Hi * Ri + Hj * Rj + Hk * Rk
          # TANGENT FRAME: Geometric phase σj
          dγdt = - Hj * Rj + (Hi * Ri + Hk * Rk) * Rj ** 2 / (Ri ** 2 + Rk ** 2)
          # DARBOUX FRAME: Geometric phase σj
          denom = (Hi * Hi + Hj * Hj + Hk * Hk - dξdt ** 2)
          H_cross_R_i = Hj * Rk - Hk * Rj
          H_cross_R_j = Hk * Ri - Hi * Rk
          H_cross_R_k = Hi * Rj - Hj * Ri
          Hdt_Rdt = Hi_dt * H_cross_R_i + Hj_dt * H_cross_R_j + Hk_dt * H_cross_R_k
          Rdt_Rdt = H_cross_R_i * H_cross_R_i + H_cross_R_j * H_cross_R_j + H_cross_R_k * H_cross_R_k
          DARBOUX_dγdt = (Hdt_Rdt - Rdt_Rdt * dξdt) / denom
      elif principle_axis == 'k':  # Ri, Rj, Rk = np.sin(μ0) * np.cos(ν0), np.sin(μ0) * np.sin(ν0), np.cos(μ0)
          Ri, Rj, Rk = self.evolve_R(μ0, ν0, rotor, principle_axis)
          # DYNAMIC PHASE:
          dξdt = Hi * Ri + Hj * Rj + Hk * Rk
          # TANGENT FRAME: Geometric phase σk
          dγdt = - Hk * Rk + (Hi * Ri + Hj * Rj) * Rk ** 2 / (Ri ** 2 + Rj ** 2)
          # DARBOUX FRAME: Geometric phase σk
          denom = (Hi * Hi + Hj * Hj + Hk * Hk - dξdt ** 2)
          H_cross_R_i = Hj * Rk - Hk * Rj
          H_cross_R_j = Hk * Ri - Hi * Rk
          H_cross_R_k = Hi * Rj - Hj * Ri
          Hdt_Rdt = Hi_dt * H_cross_R_i + Hj_dt * H_cross_R_j + Hk_dt * H_cross_R_k
          Rdt_Rdt = H_cross_R_i * H_cross_R_i + H_cross_R_j * H_cross_R_j + H_cross_R_k * H_cross_R_k
          DARBOUX_dγdt = (Hdt_Rdt - Rdt_Rdt * dξdt) / denom
      else:
          return 0
      # integrate
      ξ = (np.sum(dξdt * dt, axis=0) - dξdt[0, :] * dt).reshape(dim[0], dim[1])
      γ = (np.sum(dγdt * dt, axis=0) - dγdt[0, :] * dt).reshape(dim[0], dim[1])
      DARBOUX_γ = (np.sum(DARBOUX_dγdt * dt, axis=0) - DARBOUX_dγdt[0, :] * dt).reshape(dim[0], dim[1])
      # edge nans
      ξ[:, 0], ξ[:, -1], ξ[0, :], ξ[-1, :] = np.nan, np.nan, np.nan, np.nan
      γ[:, 0], γ[:, -1], γ[0, :], γ[-1, :] = np.nan, np.nan, np.nan, np.nan
      DARBOUX_γ[:, 0], DARBOUX_γ[:, -1], DARBOUX_γ[0, :], DARBOUX_γ[-1, :] = np.nan, np.nan, np.nan, np.nan
      # scale to pi
      ξ, γ, DARBOUX_γ = ξ / np.pi, γ / np.pi, DARBOUX_γ / np.pi
      return ξ, γ, DARBOUX_γ


  def unitary(self):
      # first term
      M1 = self.expmatrix(-self.t, self.li)
      dM1 = np.matmul(-self.li, self.expmatrix(-self.t, self.li))
      d2M1 = np.matmul(np.matmul(self.li, self.li), self.expmatrix(-self.t, self.li))
      # second term
      M2 = self.expmatrix(self.t / 2, self.lk)
      dM2 = np.matmul(0.5 * self.lk, self.expmatrix(self.t / 2, self.lk))
      d2M2 = np.matmul(0.25 * np.matmul(self.lk, self.lk), self.expmatrix(self.t / 2, self.lk))
      # third term
      M3 = self.expmatrix(self.t / 2, self.li)
      dM3 = np.matmul(0.5 * self.li, self.expmatrix(self.t / 2, self.li))
      d2M3 = np.matmul(0.25 * np.matmul(self.li, self.li), self.expmatrix(self.t / 2, self.li))
      # fourth term
      M4 = self.expmatrix(-self.t, self.lj)
      dM4 = np.matmul(-self.lj, self.expmatrix(-self.t, self.lj))
      d2M4 = np.matmul(np.matmul(self.lj, self.lj), self.expmatrix(-self.t, self.lj))
      # fifth term
      M5 = self.expmatrix(-self.t, self.li)
      dM5 = np.matmul(-self.li, self.expmatrix(-self.t, self.li))
      d2M5 = np.matmul(np.matmul(self.li, self.li), self.expmatrix(-self.t, self.li))
      # unitary
      U = np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(M4, M5))))
      # first derivative
      dU = np.matmul(dM1, np.matmul(M2, np.matmul(M3, np.matmul(M4, M5)))) \
           + np.matmul(M1, np.matmul(dM2, np.matmul(M3, np.matmul(M4, M5)))) \
           + np.matmul(M1, np.matmul(M2, np.matmul(dM3, np.matmul(M4, M5)))) \
           + np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(dM4, M5)))) \
           + np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(M4, dM5))))
      # second derivative
      d2Ua = np.matmul(d2M1, np.matmul(M2, np.matmul(M3, np.matmul(M4, M5)))) \
             + np.matmul(dM1, np.matmul(dM2, np.matmul(M3, np.matmul(M4, M5)))) \
             + np.matmul(dM1, np.matmul(M2, np.matmul(dM3, np.matmul(M4, M5)))) \
             + np.matmul(dM1, np.matmul(M2, np.matmul(M3, np.matmul(dM4, M5)))) \
             + np.matmul(dM1, np.matmul(M2, np.matmul(M3, np.matmul(M4, dM5))))
      d2Ub = np.matmul(dM1, np.matmul(dM2, np.matmul(M3, np.matmul(M4, M5)))) \
             + np.matmul(M1, np.matmul(d2M2, np.matmul(M3, np.matmul(M4, M5)))) \
             + np.matmul(M1, np.matmul(dM2, np.matmul(dM3, np.matmul(M4, M5)))) \
             + np.matmul(M1, np.matmul(dM2, np.matmul(M3, np.matmul(dM4, M5)))) \
             + np.matmul(M1, np.matmul(dM2, np.matmul(M3, np.matmul(M4, dM5))))
      d2Uc = np.matmul(dM1, np.matmul(M2, np.matmul(dM3, np.matmul(M4, M5)))) \
             + np.matmul(M1, np.matmul(dM2, np.matmul(dM3, np.matmul(M4, M5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(d2M3, np.matmul(M4, M5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(dM3, np.matmul(dM4, M5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(dM3, np.matmul(M4, dM5))))
      d2Ud = np.matmul(dM1, np.matmul(M2, np.matmul(M3, np.matmul(dM4, M5)))) \
             + np.matmul(M1, np.matmul(dM2, np.matmul(M3, np.matmul(dM4, M5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(dM3, np.matmul(dM4, M5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(d2M4, M5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(dM4, dM5))))
      d2Ue = np.matmul(dM1, np.matmul(M2, np.matmul(M3, np.matmul(M4, dM5)))) \
             + np.matmul(M1, np.matmul(dM2, np.matmul(M3, np.matmul(M4, dM5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(dM3, np.matmul(M4, dM5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(dM4, dM5)))) \
             + np.matmul(M1, np.matmul(M2, np.matmul(M3, np.matmul(M4, d2M5))))
      d2U = d2Ua + d2Ub + d2Uc + d2Ud + d2Ue
      # quaternion 4 vector
      U_4vector = U[:, :, 0]
      U_4vector_dt = dU[:, :, 0]
      U_4vector_dt2 = d2U[:, :, 0]
      return U_4vector, U_4vector_dt, U_4vector_dt2


  @staticmethod
  def evolve_R(μ0, ν0, rotor, principle_axis=None):
      if principle_axis == 'i':
          Ri0, Rj0, Rk0 = np.cos(μ0), np.sin(μ0) * np.sin(ν0), np.sin(μ0) * np.cos(ν0)
      elif principle_axis == 'j':
          Ri0, Rj0, Rk0 = np.sin(μ0) * np.sin(ν0), np.cos(μ0), np.sin(μ0) * np.cos(ν0)
      elif principle_axis == 'k':
          Ri0, Rj0, Rk0 = np.sin(μ0) * np.cos(ν0), np.sin(μ0) * np.sin(ν0), np.cos(μ0)
      else:
          return 0
      R0 = np.array([Ri0, Rj0, Rk0])
      # evolve the state
      R1 = np.matmul(rotor, R0)
      Ri, Rj, Rk = R1[:, 0, :], R1[:, 1, :], R1[:, 2, :]
      return Ri, Rj, Rk


  @staticmethod
  def quaternion2rotor(Q):
      # quaternion shorthand
      a = Q[:, 0]
      b = Q[:, 1]
      c = Q[:, 2]
      d = Q[:, 3]
      # rotation matrix
      rotor = np.zeros((len(a), 3, 3))
      rotor[:, 0, 0] = a ** 2 + b ** 2 - c ** 2 - d ** 2
      rotor[:, 1, 1] = a ** 2 - b ** 2 + c ** 2 - d ** 2
      rotor[:, 2, 2] = a ** 2 - b ** 2 - c ** 2 + d ** 2
      rotor[:, 0, 1] = 2 * (b * c - a * d)
      rotor[:, 1, 0] = 2 * (b * c + a * d)
      rotor[:, 0, 2] = 2 * (b * d + a * c)
      rotor[:, 2, 0] = 2 * (b * d - a * c)
      rotor[:, 1, 2] = 2 * (c * d - a * b)
      rotor[:, 2, 1] = 2 * (c * d + a * b)
      return rotor


  @staticmethod
  def expmatrix(theta, operator):
      theta = theta.reshape(len(theta), 1, 1)
      return np.cos(theta) * np.identity(4) + np.sin(theta) * operator


  @staticmethod
  def cayley():
      I = np.array([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]])
      # left cayley matrices
      li = np.array([[0, -1, 0, 0],
                     [1, 0, 0, 0],
                     [0, 0, 0, -1],
                     [0, 0, 1, 0]])
      lj = np.array([[0, 0, -1, 0],
                     [0, 0, 0, 1],
                     [1, 0, 0, 0],
                     [0, -1, 0, 0]])
      lk = np.array([[0, 0, 0, -1],
                     [0, 0, -1, 0],
                     [0, 1, 0, 0],
                     [1, 0, 0, 0]])
      # right cayley matrices
      ri = np.array([[0, 1, 0, 0],
                     [-1, 0, 0, 0],
                     [0, 0, 0, -1],
                     [0, 0, 1, 0]])
      rj = np.array([[0, 0, 1, 0],
                     [0, 0, 0, 1],
                     [-1, 0, 0, 0],
                     [0, -1, 0, 0]])
      rk = np.array([[0, 0, 0, 1],
                     [0, 0, -1, 0],
                     [0, 1, 0, 0],
                     [-1, 0, 0, 0]])
      return I, li, lj, lk, ri, rj, rk


  # PHASE GRIDS PLOT ###################################################################################################

  def plot_phase_grids(self, i, j, k, clip_geometricphase=False, progress=0, figsize=(17.8, 20.2)):
      if clip_geometricphase:
          string1, string2 = 'bounded geometric phase', ' [-π, π]'
          # clip geometric phase between [0, 2]
          i[1], i[2] = np.mod(i[1], 2), np.mod(i[2], 2)
          j[1], j[2] = np.mod(j[1], 2), np.mod(j[2], 2)
          k[1], k[2] = np.mod(k[1], 2), np.mod(k[2], 2)
          # shift geometric phase between [-1, 1]
          i[1][i[1] > 1], i[2][i[2] > 1] = i[1][i[1] > 1] - 2, i[2][i[2] > 1] - 2
          j[1][j[1] > 1], j[2][j[2] > 1] = j[1][j[1] > 1] - 2, j[2][j[2] > 1] - 2
          k[1][k[1] > 1], k[2][k[2] > 1] = k[1][k[1] > 1] - 2, k[2][k[2] > 1] - 2
      else:
          string1, string2 = 'unbounded geometric phase', ''
      # extract dimensions
      dim = i[0].shape
      # coordinates
      ν0, μ0 = np.meshgrid(np.linspace(0, 2, dim[1]), np.linspace(0, 1, dim[0]))
      # sphere grids
      ν, μ = np.meshgrid(np.linspace(0, 2 * np.pi, dim[1]), np.linspace(0, np.pi, dim[0]))
      x_i, y_i, z_i = self.mesh_i(μ[1:-1, 1:-1], ν[1:-1, 1:-1])
      x_j, y_j, z_j = self.mesh_j(μ[1:-1, 1:-1], ν[1:-1, 1:-1])
      x_k, y_k, z_k = self.mesh_k(μ[1:-1, 1:-1], ν[1:-1, 1:-1])
      # DARBOUX FRAME ##################################################################################################
      # declare figure
      fig1, ax1 = plt.figure('Darboux frame: phase grids', figsize=figsize), []
      # one
      for idx, global_phase, axis in zip([0, 1, 2],
                                         [i[0] + i[2], j[0] + j[2], k[0] + k[2]],
                                         ['=  \u03C3i  )', '=  \u03C3j  )', '=  \u03C3k  )']):
          ax1.append(fig1.add_subplot(4, 3, idx + 1))
          # SUBPLOT (3,3,1): global phase
          ax1[idx].set_title('Global phase (principle axis ' + axis)
          plt.pcolormesh(ν0, μ0, global_phase)
          plt.colorbar(orientation='vertical')  # , shrink=0.5
          ax1[idx].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
          ax1[idx].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
          ax1[idx].set_yticklabels(["0", " ", "π/2", " ", "π "])
          ax1[idx].set_xticklabels(["0", " ", " ", " ", " π ", " ", " ", " ", "2π "])
          ax1[idx].set_xlabel("                        ν\u2080", fontsize=14)
          ax1[idx].set_ylabel("                        μ\u2080", fontsize=14)
      for idx, geometric_phase in zip([3, 4, 5], [i[2], j[2], k[2]]):
          ax1.append(fig1.add_subplot(4, 3, idx + 1))
          # SUBPLOT (3,3,1): global phase
          ax1[idx].set_title('Geometric phase')
          plt.pcolormesh(ν0, μ0, geometric_phase)
          plt.colorbar(orientation='vertical')  # , shrink=0.5
          ax1[idx].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
          ax1[idx].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
          ax1[idx].set_yticklabels(["0", " ", "π/2", " ", "π "])
          ax1[idx].set_xticklabels(["0", " ", " ", " ", " π ", " ", " ", " ", "2π "])
          ax1[idx].set_xlabel("                        ν\u2080", fontsize=14)
          ax1[idx].set_ylabel("                        μ\u2080", fontsize=14)
      for idx, dynamic_phase in zip([6, 7, 8], [i[0], j[0], k[0]]):
          ax1.append(fig1.add_subplot(4, 3, idx + 1))
          # SUBPLOT (3,3,1): global phase
          ax1[idx].set_title('Dynamic phase')
          plt.pcolormesh(ν0, μ0, dynamic_phase)
          plt.colorbar(orientation='vertical')  # , shrink=0.5
          ax1[idx].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
          ax1[idx].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
          ax1[idx].set_yticklabels(["0", " ", "π/2", " ", "π "])
          ax1[idx].set_xticklabels(["0", " ", " ", " ", " π ", " ", " ", " ", "2π "])
          ax1[idx].set_xlabel("                        ν\u2080", fontsize=14)
          ax1[idx].set_ylabel("                        μ\u2080", fontsize=14)
      # sphere plotting
      Di, Dj, Dk = (i[0] + i[2])[1:-1, 1:-1], (j[0] + j[2])[1:-1, 1:-1], (k[0] + k[2])[1:-1, 1:-1]
      Dmap_i = (Di - Di.min()) / (Di.max() - Di.min())
      Dmap_j = (Dj - Dj.min()) / (Dj.max() - Dj.min())
      Dmap_k = (Dk - Dk.min()) / (Dk.max() - Dk.min())
      ax1.append(fig1.add_subplot(4, 3, 10, projection='3d'))
      ax1.append(fig1.add_subplot(4, 3, 11, projection='3d'))
      ax1.append(fig1.add_subplot(4, 3, 12, projection='3d'))
      ax1[9].plot_surface(x_i, y_i, z_i, cstride=1, rstride=1, facecolors=cm.summer(Dmap_i))
      ax1[10].plot_surface(x_j, y_j, z_j, cstride=1, rstride=1, facecolors=cm.summer(Dmap_j))
      ax1[11].plot_surface(x_k, y_k, z_k, cstride=1, rstride=1, facecolors=cm.summer(Dmap_k))
      # xyz ticks
      for idx in [9, 10, 11]:
          ax1[idx].set_xticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
          ax1[idx].set_yticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
          ax1[idx].set_zticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
          ax1[idx].set_xticklabels(["-1.0", "", "-0.5", "", "0", "", "0.5", "", "1.0"])
          ax1[idx].set_yticklabels(["-1.0", "", "-0.5", "", "0", "", "0.5", "", "1.0"])
          ax1[idx].set_zticklabels(["-1.0", "", "-0.5", "", "0", "", "0.5", "", "1.0"])
          ax1[idx].set_xlabel("i")
          ax1[idx].set_ylabel("j")
          ax1[idx].set_zlabel("k")
          ax1[idx].set_title('Global phase')
      fig1.suptitle('Darboux frame [phase grids]: ' + string1 + string2, fontsize=16)
      # save figure
      fig1.savefig("figures/Darboux_" + string1 + ".png")
      plt.close()
      self.drawProgressBar(progress + 25)
      # TANGENT FRAME ##################################################################################################
      # declare figure
      fig2, ax2 = plt.figure('Tangent frame: phase grids', figsize=figsize), []
      # one
      for idx, global_phase, axis in zip([0, 1, 2],
                                         [i[0] + i[1], j[0] + j[1], k[0] + k[1]],
                                         ['=  \u03C3i  )', '=  \u03C3j  )', '=  \u03C3k  )']):
          ax2.append(fig2.add_subplot(4, 3, idx + 1))
          # SUBPLOT (3,3,1): global phase
          ax2[idx].set_title('Global phase (principle axis ' + axis)
          plt.pcolormesh(ν0, μ0, global_phase)
          plt.colorbar(orientation='vertical')  # , shrink=0.5
          ax2[idx].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
          ax2[idx].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
          ax2[idx].set_yticklabels(["0", " ", "π/2", " ", "π "])
          ax2[idx].set_xticklabels(["0", " ", " ", " ", " π ", " ", " ", " ", "2π "])
          ax2[idx].set_xlabel("                        ν\u2080", fontsize=14)
          ax2[idx].set_ylabel("                        μ\u2080", fontsize=14)
      for idx, geometric_phase in zip([3, 4, 5], [i[1], j[1], k[1]]):
          ax2.append(fig2.add_subplot(4, 3, idx + 1))
          # SUBPLOT (3,3,1): global phase
          ax2[idx].set_title('Geometric phase')
          plt.pcolormesh(ν0, μ0, geometric_phase)
          plt.colorbar(orientation='vertical')  # , shrink=0.5
          ax2[idx].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
          ax2[idx].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
          ax2[idx].set_yticklabels(["0", " ", "π/2", " ", "π "])
          ax2[idx].set_xticklabels(["0", " ", " ", " ", " π ", " ", " ", " ", "2π "])
          ax2[idx].set_xlabel("                        ν\u2080", fontsize=14)
          ax2[idx].set_ylabel("                        μ\u2080", fontsize=14)
      for idx, dynamic_phase in zip([6, 7, 8], [i[0], j[0], k[0]]):
          ax2.append(fig2.add_subplot(4, 3, idx + 1))
          # SUBPLOT (3,3,1): global phase
          ax2[idx].set_title('Dynamic phase')
          plt.pcolormesh(ν0, μ0, dynamic_phase)
          plt.colorbar(orientation='vertical')  # , shrink=0.5
          ax2[idx].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
          ax2[idx].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
          ax2[idx].set_yticklabels(["0", " ", "π/2", " ", "π "])
          ax2[idx].set_xticklabels(["0", " ", " ", " ", " π ", " ", " ", " ", "2π "])
          ax2[idx].set_xlabel("                        ν\u2080", fontsize=14)
          ax2[idx].set_ylabel("                        μ\u2080", fontsize=14)
      # sphere plotting
      Di, Dj, Dk = (i[0] + i[1])[1:-1, 1:-1], (j[0] + j[1])[1:-1, 1:-1], (k[0] + k[1])[1:-1, 1:-1]
      Dmap_i = (Di - Di.min()) / (Di.max() - Di.min())
      Dmap_j = (Dj - Dj.min()) / (Dj.max() - Dj.min())
      Dmap_k = (Dk - Dk.min()) / (Dk.max() - Dk.min())
      ax2.append(fig2.add_subplot(4, 3, 10, projection='3d'))
      ax2.append(fig2.add_subplot(4, 3, 11, projection='3d'))
      ax2.append(fig2.add_subplot(4, 3, 12, projection='3d'))
      ax2[9].plot_surface(x_i, y_i, z_i, cstride=1, rstride=1, facecolors=cm.summer(Dmap_i))
      ax2[10].plot_surface(x_j, y_j, z_j, cstride=1, rstride=1, facecolors=cm.summer(Dmap_j))
      ax2[11].plot_surface(x_k, y_k, z_k, cstride=1, rstride=1, facecolors=cm.summer(Dmap_k))
      # xyz ticks
      for idx in [9, 10, 11]:
          ax2[idx].set_xticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
          ax2[idx].set_yticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
          ax2[idx].set_zticks([-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0])
          ax2[idx].set_xticklabels(["-1.0", "", "-0.5", "", "0", "", "0.5", "", "1.0"])
          ax2[idx].set_yticklabels(["-1.0", "", "-0.5", "", "0", "", "0.5", "", "1.0"])
          ax2[idx].set_zticklabels(["-1.0", "", "-0.5", "", "0", "", "0.5", "", "1.0"])
          ax2[idx].set_xlabel("i")
          ax2[idx].set_ylabel("j")
          ax2[idx].set_zlabel("k")
          ax2[idx].set_title('Global phase')
      fig2.suptitle('Tangent frame [phase grids]: ' + string1 + string2, fontsize=16)
      # save figure
      fig2.savefig("figures/Tangent_" + string1 + ".png")
      plt.close()
      self.drawProgressBar(progress + 50)


  @staticmethod
  def mesh_i(μ, ν):
      # surface # np.cos(μ0), np.sin(μ0) * np.sin(ν0), np.sin(μ0) * np.cos(ν0)
      x = np.cos(μ)
      y = np.sin(μ) * np.sin(ν)
      z = np.sin(μ) * np.cos(ν)
      return x, y, z


  @staticmethod
  def mesh_j(μ, ν):
      # surface # np.sin(μ0) * np.sin(ν0), np.cos(μ0), np.sin(μ0) * np.cos(ν0)
      x = np.sin(μ) * np.sin(ν)
      y = np.cos(μ)
      z = np.sin(μ) * np.cos(ν)
      return x, y, z


  @staticmethod
  def mesh_k(μ, ν):
      # surface # np.sin(μ0) * np.cos(ν0), np.sin(μ0) * np.sin(ν0), np.cos(μ0)
      x = np.sin(μ) * np.cos(ν)
      y = np.sin(μ) * np.sin(ν)
      z = np.cos(μ)
      return x, y, z

  @staticmethod
  # progress bar for data loading
  def drawProgressBar(percent, barLen=40):
      sys.stdout.write("\r")
      progress = ""
      for i in range(barLen):
          if i < int(barLen * percent / 100):
              progress += "="
          else:
              progress += " "
      sys.stdout.write("[ %s ] %.1f%% " % (progress, percent))
      sys.stdout.flush()

# BLOCH SPHERE #########################################################################################################

class S2sphere:
  def __init__(self, vertices=44):
    # meshgrids μ + ν
    μ, ν = np.meshgrid(np.linspace(0, np.pi, vertices).astype('float64'),
                                 np.linspace(0, 2 * np.pi, vertices).astype('float64'))
    self.x, self.y, self.z, self.tri = self.delunay_triangulation(μ, ν)

  # delunay triangulation
  @staticmethod
  def delunay_triangulation(μ, ν):
      # surface
      i = np.ravel(np.sin(μ) * np.cos(ν))
      j = np.ravel(np.sin(μ) * np.sin(ν))
      k = np.ravel(np.cos(μ))
      # delunay triangulation
      return i, j, k, Triangulation(np.ravel(μ), np.ravel(ν))

########################################################################################################################