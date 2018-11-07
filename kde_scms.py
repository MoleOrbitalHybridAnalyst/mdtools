# TODO data-dependent kernel function
# Ozertem, U. Locally Defined Principal Curves and Surfaces; 2011; Vol. 12.
import numpy as np

class KDE:
    """
    KDE class:
        initializer:
            define kernel function from given covariance
            compute inverse of the kernel covariance S_K^{-1} (\Sigma_i^{-1})
            read and store the training points {x_i}
        estimate(x):
            compute c_i = K(x, x_i)
            compute p(x) = mean({c_i})
        gradient(x):
            (assuming c_i have been updated)
            compute u_i = S_K^{-1} \cdot (x - x_i)
            compute g(x) = - {c_i} \cdot {u_i} / N
        hessian(x):
            (assuming u_i have been updated)
            compute H(x) = (c_i * {u_i})^T\cdot{u_i} / N - p(x) S_K^{-1}
        mean_shift(x):
            (assuming c_i have been updated)
            compute m(x) = (sum_i c_i S_K^{-1})^{-1} \cdot \sum_i c_i S_K^{-1} \cdot x_i
                         = p(x)^{-1} S_K \cdot S_K^{-1} \cdot ({c_i} \cdot {x_i}) / N
                         = p(x)^{-1} ({c_i} \cdot {x_i}) / N
    """

    def __init__(self, data = None, ker_cov = None):

        if data is None:
            raise Exception("no data provided")
        if ker_cov is None:
            raise Exception("no kernel covariance provided")
        self.X = np.array(data)
        self.SK = np.array(ker_cov)

        xshape = self.X.shape
        if len(xshape) != 2:
            raise Exception("Nxd input data expected")
        covshape = self.SK.shape
        if len(covshape) != 2 or covshape[0] != covshape[1]:
            raise Exception("dxd kernel covariance matrix expected")
        if xshape[1] != covshape[0]:
            raise Exception("inconsistent dimension between data and ker_cov")

        self.N = xshape[0]
        self.ISK = np.linalg.inv(self.SK)
        # x_: Nxd
        self.K = lambda x_: \
                np.exp(- np.sum(np.dot(x_, self.ISK) * x_, axis = 1) / 2)

        self.c = self.p = None
        self.u = self.g = None
        self.H = self.m = None

    def estimate(self, x):
        self.c = self.K(x - self.X)
        self.p = np.mean(self.c)
        return self.p

    def gradient(self, x):
        self.u = np.dot(x - self.X, self.ISK)
        self.g = - np.dot(self.c, self.u) / self.N
        return self.g

    def hessian(self, x):
        self.H = np.dot(self.c * self.u.transpose(), self.u) / self.N
        self.H -= self.p * self.ISK
        return self.H

    def mean_shift(self, x): 
        self.m = np.dot(self.c, self.X) / self.p / self.N
        return self.m

class SCMS:
    """
    SCMS class:
        Initializer:
            store KDE object using data and ker_cov
            store the targetted dimension d_
        covariance_inverse(x):
            compute IS (\Sigma^{-1}) = (g \cross g / p - H) / p
        normal_vectors():
            (assuming IS is ready)
            compute V = [v1 ... v_(d-d_)]
        project(x):
            (assuming V is ready)
            compute xp = V \cdot V^T \cdot m
        measure_orthogonality():
            compute |g^T \cdot V^T g| / (|g|\dot|V^Tg|)
    """
    _update_ = False

    def __init__(self, dim, data = None, ker_cov = None):
        self.dim = dim
        self.kde = KDE(data = data, ker_cov = ker_cov)
        if dim >= len(ker_cov):
            raise \
                Exception("target dimension not lower than original dimension")

        self.m = self.IS = self.V = None
        self.xp = self.cosTheta = None

    def covariance_inverse(self, x):
        if self._update_:
            self.kde.estimate(x)
            x = kde.mean_shift(x)
            self.m = x

        self.kde.estimate(x)
        self.kde.gradient(x)
        self.kde.hessian(x)
        self.IS = (
            np.outer(self.kde.g, self.kde.g) / self.kde.p - self.kde.H
            ) / self.kde.p
        return self.IS

    def normal_vectors(self):
        eigsys = np.linalg.eigh(self.IS)
        self.V = np.transpose(eigsys[1][:-self.dim])
        return self.V

    def project(self, x):
        if not self._update_:
            self.m = self.kde.mean_shift(x)
        # else self.m has been computed in covariance_inverse
        self.xp = np.dot(self.V, np.dot(self.V.transpose(), self.m))
        return self.xp
