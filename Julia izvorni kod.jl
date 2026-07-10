using LinearAlgebra

function svediNaHessenberga(A::Matrix{Float64})
    H = copy(A)
    n = size(H, 1)
    
    for k in 1:n-2
        x = H[k+1:n, k]
        norm_x = norm(x)
        if norm_x == 0 continue 
        end
        alpha = x[1] >= 0 ? -norm_x : norm_x
        v = copy(x)
        v[1] -= alpha
        beta = 2.0 / dot(v, v)
        for j in k:n
            s = beta * dot(v, H[k+1:n, j])
            H[k+1:n, j] -= v * s 
        end
        for i in 1:n
            s = beta * dot(H[i, k+1:n], v)
            H[i, k+1:n] -= (s * v')' 
        end
       for i in k+2:n
            H[i, k] = 0.0
        end
    end
    return H
end


function izracunajSvojstveneVrijednosti(Matrica_S::Matrix{Float64}; tolerancija=1e-12)
    n = size(Matrica_S, 1)
    kraj = false
    A_iter = svediNaHessenberga(Matrica_S)
    while !kraj
        kraj = true
        F = qr(A_iter)
        Q = Matrix(F.Q)
        R = F.R
        A_iter = R * Q
        for i in 2:n
            for j in 1:i-1
                if abs(A_iter[i, j]) > tolerancija
                    kraj = false
                    break
                end
            end
            if !kraj break end
        end
    end
    lambda = [A_iter[i, i] for i in 1:n]
    return lambda
end

A_test = [ 5.0  2.0  1.0  0.0 ;
           2.0  8.0  3.0  1.0 ;
           1.0  3.0  4.0  2.0 ;
           0.0  1.0  2.0  6.0 ]

println("Testna matrica A_test:")
display(A_test)

referentne_vrijednosti = sort(eigvals(A_test), rev=true)
moji_rezultati = izracunajSvojstveneVrijednosti(A_test, tolerancija=1e-12)

println("Referentne vrijednosti (eigvals): ", referentne_vrijednosti)
println("Rezultati tvoje QR metode:       ", moji_rezultati)

razlika = norm(referentne_vrijednosti - moji_rezultati)
println("Numerička razlika: ", razlika)


A_test2 = [ 10.0   1.0   2.0   3.0 ;
             4.0  12.0   1.0   0.0 ;
             2.0   1.0   8.0   2.0 ;
             1.0   0.0   5.0  15.0 ]

println("Matrica A_test2:")
display(A_test2)

referentne_vrijednosti2 = sort(eigvals(A_test2), rev=true)
moji_rezultati2 = izracunajSvojstveneVrijednosti(A_test2, tolerancija=1e-14)
println("Referentne vrijednosti: ", referentne_vrijednosti2)
println("Rezultati algoritma: ", moji_rezultati2)

razlika2 = norm(referentne_vrijednosti2 - moji_rezultati2)
println("Numerička razlika : ", razlika2)

A_test3 = [ 1.0   1.0   1.0   1.0 ;
            1.0   2.0   3.0   4.0 ;
            1.0   3.0   6.0  10.0 ;
            1.0   4.0  10.0  20.0 ]

println("Matrica A_test3:")
display(A_test3)

referentne_vrijednosti3 = sort(eigvals(A_test3), rev=true)
moji_rezultati3 = izracunajSvojstveneVrijednosti(A_test3, tolerancija=1e-14)
println("Referentne vrijednosti: ", referentne_vrijednosti3)
println("Rezultati algoritma: ", moji_rezultati3)

razlika3 = norm(referentne_vrijednosti3 - moji_rezultati3)
println("Numerička razlika: ", razlika3)

A_test4 = [ 1.0       0.5       0.33333333 0.25 ;
            0.5       0.33333333 0.25       0.2 ;
            0.33333333 0.25       0.2        0.16666667 ;
            0.25       0.2        0.16666667 0.14285714 ]

println("Matrica A_test4:")
display(A_test4)

referentne_vrijednosti4 = sort(eigvals(A_test4), rev=true)
moji_rezultati4 = izracunajSvojstveneVrijednosti(A_test4, tolerancija=1e-14)
println("Referentne vrijednosti: ", referentne_vrijednosti4)
println("Rezultati algoritma: ", moji_rezultati4)

razlika4 = norm(referentne_vrijednosti4 - moji_rezultati4)
println("Numerička razlika: ", razlika4)