from kde_scms import *
import matplotlib.pyplot as plt

def plot():
    plt.scatter(xs, ys)
    plt.show()

def shift2converge(x, epsilon = 0.001):
    dev = 1.0
    while dev > epsilon:
        x_pre = x
        test_kde.estimate(x)
        x = test_kde.mean_shift(x)
        dev = np.sum((x - x_pre)**2)
    return x

def project2converge(x, epsilon = 0.001):
    dev = 1.0
    while dev > epsilon:
        x_pre = x
        test_scms.covariance_inverse(x)
        test_scms.normal_vectors()
        x = test_scms.project(x)
        dev = np.sum((x - x_pre)**2)
    return x


np.random.seed(12345)

xs = np.arange(0, 3 * np.pi, 0.02)
ys = np.sin(xs)
ys += 0.4 * np.random.sample(len(xs))
xs += 0.2 * np.random.sample(len(xs))

xs2 = np.ones(20) * np.pi * 1.5 + 0.2 * np.random.sample(20)
ys2 = np.linspace(-1, -1.5, 20) + 0.4 * np.random.sample(20)

xs = np.append(xs, xs2)
ys = np.append(ys, ys2)

#plot()

test_kde = KDE(data = np.transpose([xs, ys]), ker_cov = [[1,0],[0,0.5]])
# test mean-shift
x0 = [1.5 * np.pi, -1.4]
print(x0)
print(test_kde.estimate(x0))
x0_ = shift2converge(x0)
print(x0_)
print(test_kde.estimate(x0_))
x0 = [xs[-1], ys[-1]]
print(x0)
print(test_kde.estimate(x0))
x0_ = shift2converge(x0)
print(x0_)
print(test_kde.estimate(x0_))

test_scms = SCMS(1, data = np.transpose([xs, ys]), ker_cov = [[1,0],[0,0.5]])
# test project
x0 = [xs[-1], ys[-1]]
print(x0)
#test_scms.covariance_inverse(x0)
#test_scms.normal_vectors()
#x0_ = test_scms.project(x0)
x0_ = project2converge(x0)
print(x0_)

curve = []
for i in range(len(xs)):
    nbatch = "
