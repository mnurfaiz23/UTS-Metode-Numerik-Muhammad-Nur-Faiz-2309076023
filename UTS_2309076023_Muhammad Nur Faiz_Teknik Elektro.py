# Import pustaka numpy untuk melakukan perhitungan numerik
import numpy as np

# Import pustaka matplotlib.pyplot untuk membuat visualisasi grafik
import matplotlib.pyplot as plt

# Diketahui
L = 0.5  # Induktansi dalam satuan Henry
C = 10e-6  # Kapasitansi dalam satuan Farad
fd = 1000  # Frekuensi yang diharapkan dalam satuan Hertz (Hz)

# Definisikan fungsi untuk menghitung frekuensi f(R) sebagai fungsi dari R
# Fungsi ini menghitung frekuensi berdasarkan nilai-nilai induktansi (L) dan kapasitansi (C)
def f_R(R):
    # Formula f(R) berasal dari rumus resonansi dalam sirkuit LC
    return (1 / (2 * np.pi)) * np.sqrt((1 / (L * C)) - (R**2 / (4 * L**2)))

# Definisikan fungsi F(R) yang merupakan perbedaan antara f(R) dan frekuensi yang diharapkan (fd)
def F_R(R):
    return f_R(R) - fd

# Mengecek nilai-nilai F_R pada dua titik awal yang dipilih yaitu R=10 dan R=100
# Hal ini berguna untuk memeriksa perubahan tanda dari F_R, yang dibutuhkan oleh metode biseksi
print("F_R(10):", F_R(10))
print("F_R(100):", F_R(100))

# Visualisasi grafik F(R) untuk menemukan interval yang tepat di mana akar berada
R_values = np.linspace(0, 200, 1000)  # Mengambil 1000 nilai dari 0 hingga 200 untuk R
F_R_values = F_R(R_values)  # Menghitung nilai F(R) untuk setiap nilai R

# Membuat plot dari F(R) terhadap R untuk membantu dalam analisis
plt.plot(R_values, F_R_values)
plt.axhline(0, color='red', lw=0.5)  # Menambahkan garis horizontal di y=0 sebagai acuan
plt.xlabel('Nilai R')  # Menambahkan label sumbu-x
plt.ylabel('F(R)')  # Menambahkan label sumbu-y
plt.title('Plot dari F(R)')  # Menambahkan judul grafik
plt.grid(True)
plt.show()

# Implementasi metode biseksi untuk mencari akar F(R)
# Metode biseksi mencari akar fungsi dengan membagi dua interval secara berulang
def bisection_method(func, a, b, tol=1e-6):
    # Memeriksa apakah ada perubahan tanda dalam interval [a, b]
    if func(a) * func(b) > 0:
        raise ValueError("Fungsi tidak berubah tanda pada interval ini.")

    # Daftar untuk menyimpan nilai R hasil iterasi
    R_values_bisection = []
    while (b - a) / 2.0 > tol:
        # Menentukan nilai tengah dari interval [a, b]
        mid = (a + b) / 2.0
        R_values_bisection.append(mid)  # Menyimpan nilai tengah

        # Memeriksa apakah akar telah ditemukan
        if func(mid) == 0:
            return mid, R_values_bisection
        elif func(a) * func(mid) < 0:
            b = mid
        else:
            a = mid

    return mid, R_values_bisection  # Mengembalikan nilai tengah dan daftar nilai R

# Menjalankan metode biseksi pada F(R) dengan batas awal 0 dan 200
R_bisection, R_values_bisection = bisection_method(F_R, 0, 200)

# Mendefinisikan turunan dari F(R) untuk metode Newton-Raphson
def dF_R(R):
    return -(R / (2 * np.pi * L**2 * np.sqrt((1 / (L * C)) - (R**2 / (4 * L**2)))))

# Implementasi metode Newton-Raphson untuk mencari akar F(R)
# Metode ini menggunakan turunan untuk memperbaiki perkiraan akar setiap iterasi
def newton_raphson(func, dfunc, x0, tol=1e-6, max_iter=100):
    x = x0  # Memulai dari nilai awal x0
    R_values_newton = []
    for i in range(max_iter):
        x_new = x - func(x) / dfunc(x)  # Memperbaiki x menggunakan formula Newton-Raphson
        R_values_newton.append(x_new)

        # Memeriksa apakah solusi sudah mencapai toleransi yang ditentukan
        if abs(x_new - x) < tol:
            return x_new, R_values_newton
        x = x_new

    raise ValueError("Metode Newton-Raphson gagal konvergen.")

# Menjalankan metode Newton-Raphson pada F(R) dengan titik awal 50
R_newton, R_values_newton = newton_raphson(F_R, dF_R, 50)

# Visualisasi konvergensi dari kedua metode untuk membandingkan hasil iterasi
def plot_comparison():
    plt.plot(R_values_bisection, label="Metode Bisection")
    plt.plot(R_values_newton, label="Metode Newton-Raphson")
    plt.xlabel('Iterasi')  # Menambahkan label untuk sumbu x
    plt.ylabel('Nilai R')  # Menambahkan label untuk sumbu y
    plt.title('Perbandingan Konvergensi')  # Menambahkan judul plot
    plt.legend()
    plt.grid(True)
    plt.show()

# Memanggil fungsi untuk membuat plot perbandingan konvergensi
plot_comparison()

# Implementasi eliminasi Gauss untuk menyelesaikan sistem persamaan linear
def gauss_elimination(A, B):
    n = len(B)  # Menentukan jumlah persamaan dalam sistem

    # Proses eliminasi maju untuk mengubah matriks ke bentuk segitiga atas
    for i in range(n):
        for j in range(i+1, n):
            # Menentukan rasio untuk eliminasi
            ratio = A[j][i] / A[i][i]
            # Melakukan eliminasi pada setiap elemen baris
            for k in range(n):
                A[j][k] = A[j][k] - ratio * A[i][k]
            # Memodifikasi elemen vektor B
            B[j] = B[j] - ratio * B[i]

    # Array untuk menyimpan solusi
    X = [0 for i in range(n)]
    # Melakukan substitusi mundur untuk menemukan nilai-nilai solusi
    X[n-1] = B[n-1] / A[n-1][n-1]
    for i in range(n-2, -1, -1):
        X[i] = B[i]
        for j in range(i+1, n):
            X[i] = X[i] - A[i][j] * X[j]
        X[i] = X[i] / A[i][i]

    return X

# Mendefinisikan sistem persamaan linear
A = [[4, -1, -1], [-1, 3, -1], [-1, -1, 5]]
B = [5, 3, 4]

# Menyelesaikan sistem persamaan menggunakan eliminasi Gauss dan menampilkan hasil
solutions = gauss_elimination(A, B)
print(solutions)
