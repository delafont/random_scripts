#!/usr/bin/env julia

function generate_key(N=20)
    key = rand(0:1, N)
    while countnz(key) == 0     # on ne veut pas une cle avec 0 partout
        key = rand(0:1, N)
    end
    return key
end

function encode_1bit(bit, key)
    last1 = findin(key,1)[end]  # dernier bit non-zero
    c = Float64[]
    total = 0
    #noise = 0
    noise = (rand()*2 - 1)/5    # "bruit", inferieur a 1/4
    # noise = 2^(-sqrt(length(key)))
    for i=1:length(key)
        a = rand()*2 - 1
        if i == last1
            push!(c, key[i]-total-(1-bit)+noise)
        else
            push!(c, a)
        end
        if key[i] != 0
            total = total + a
        end
    end
    return c
end

function encode(message::String, key)
    message = map(int, split(message,""))  # string -> array de 0/1
    encrypted = zeros(length(key), length(message))
    for (i,bit) in enumerate(message)
        c = encode_1bit(bit, key)
        encrypted[:,i] = c
    end
    return encrypted
end


function decode_1bit(encrypted_bit, key)
    return abs(round((key' * encrypted_bit)[1] % 2))
end

function decode(encrypted, key)
    M = size(encrypted,2)
    decrypted = zeros(Int, M)
    for i in 1:M
        decrypted[i] = decode_1bit(encrypted[:,i], key)
    end
    return join(decrypted)
end



function example()
    K = generate_key(20)
    message = bin(char(rand(65:90)))  # lettre aleatoire
    C = encode(message, K)
    decrypted = decode(C, K)
    #println("Encrypted message:", C)
    println("original: ", string(join(message)), " == ",
            string(join(decrypted)), ": decrypted",
            " (key=", string(join(K)), ")")
    @assert message == decrypted
end


# Now suppose we don't have the key
# We generate N encryptions of 1 and solve the system
function hack()
    K = generate_key(10)
    N = length(K)
    encrypt = zeros(N, N)
    for i=1:N
        encrypt[:,i] = encode_1bit(1, K)
    end
    O = ones(N)
    guess = int(encrypt' \ O)
    @assert guess == K
end

for _ = 1:10 example() end
for _ = 1:10 hack() end


function addbin()
    K = generate_key(20)
    c1 = encode(bits(9), K)
    c2 = encode(bits(6), K)
    cxor = (c1 .+ c2)
    cand = (c1 .* c2)
    rxor = decode(cxor, K)
    rand = decode(cand, K)
    r = parseint(decode(c, K), 2)
end
