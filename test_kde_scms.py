from kde_scms import *
import matplotlib.pyplot as plt

def plot():
    plt.scatter(xs, ys)
    plt.show()

def shift2converge(x, epsilon = 0.0001):
    dev = 1.0
    while dev > epsilon:
        x_pre = x
        test_kde.estimate(x)
        x = test_kde.mean_shift(x)
        dev = np.sum((x - x_pre)**2)
    return x

def project2converge(x, epsilon = 0.0001):
    dev = 1.0
    while dev > epsilon:
        x_pre = x
        test_scms.covariance_inverse(x)
        test_scms.normal_vectors()
        x = test_scms.project(x)
        dev = np.sum((x - x_pre)**2)
        #print(x)
    return x

def shift(x): 
    test_kde.estimate(x)
    return test_kde.mean_shift(x)

def project(x):
    test_scms.covariance_inverse(x)
    test_scms.normal_vectors()
    return test_scms.project(x)

def track_shift(n):
    for i in range(n):
        if i == 0:
            shifted = np.array([shift([x,y]) for x,y in zip(xs,ys)])
        else:
            shifted = np.array([shift(xy) for xy in shifted])
        plt.scatter(xs, ys)
        plt.scatter(shifted[:,0], shifted[:,1])
        plt.show()

def track_project(n):
    for i in range(n):
        if i == 0:
            projected = np.array([project([x,y]) for x,y in zip(xs,ys)])
        else:
            projected = np.array([project(xy) for xy in projected])
        plt.scatter(xs, ys)
        plt.scatter(projected[:,0], projected[:,1])
        plt.show()

def get_curve():
    curve = []
    nbatch = int(len(xs) / 100)
    for i in range(len(xs)):
        print(i, len(xs))
        curve.append(project2converge([xs[i], ys[i]]))
    return np.array(curve)

np.random.seed(12345)

xs = np.arange(0, 3 * np.pi, 0.02)
ys = np.sin(xs)
ys += 0.4 * np.random.sample(len(xs))
xs += 0.2 * np.random.sample(len(xs))

nextra = 20
xs2 = np.ones(nextra) * np.pi * 1.5 + 0.2 * np.random.sample(nextra)
ys2 = np.linspace(-1.1, -1.5, nextra) + 0.4 * np.random.sample(nextra)

xs = np.append(xs, xs2)
ys = np.append(ys, ys2)

#plot()

# kernel covariance
sigma_k = [[0.08,0], [0,0.08]]

test_kde = KDE(data = np.transpose([xs, ys]), ker_cov = sigma_k)
# test estimate
#xs_ = np.arange(0, 10, 0.01)
#ys_ = np.arange(-1.5, 1.5, 0.01)
#zs_ = [[test_kde.estimate([x, y]) for x in xs_] for y in ys_]
# test mean-shift
x0 = [1.5 * np.pi, -1.4]
print("initial", x0)
print("p", test_kde.estimate(x0))
x0_ = shift2converge(x0)
print("shifted", x0_)
print("p",test_kde.estimate(x0_))
x0 = [xs[-1], ys[-1]]
print("initial", x0)
print("p", test_kde.estimate(x0))
x0_ = shift2converge(x0)
print("projected", x0_)
print("p", test_kde.estimate(x0_))
print("g", test_kde.gradient(x0_))
print("H", test_kde.hessian(x0_))

test_scms = SCMS(1, data = np.transpose([xs, ys]), ker_cov = sigma_k)
# test project
x0 = [xs[-1], ys[-1]]
print(x0)
#test_scms.covariance_inverse(x0)
#test_scms.normal_vectors()
#x0_ = test_scms.project(x0)
x0_ = project2converge(x0)
print(x0_)

