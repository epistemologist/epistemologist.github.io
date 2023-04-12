---
layout: post
title: "Cracking a 9 Year Old Safe"

categories:
---

Recently, I came across [an old StackExchange thread](https://codegolf.stackexchange.com/q/36822/) wherein users submitted "lockers" or functions which only accepted a certain arithmetic sequence - then other users attempted to "crack" these lockers by providing said sequence. As a fun weekend challenge, I decided to crack one of these lockers that had remained uncracked for 9 years.

# The Locker

The provided function is as follows:
```py
def a(b):
    c=1
    for d in b:
        c=(c<<32)+d
    return pow(7,c,0xf494eca63dcab7b47ac21158799ffcabca8f2c6b3)==0xa3742a4abcb812e0c3664551dd3d6d2207aecb9be
```

Note that the function returns true only if $$ 7^c=\underbrace{\texttt{0xa3}\cdots\texttt{be}}_{\equiv H} $$
modulo a prime $$ p=\texttt{0xf4}\cdots\texttt{b3} $$. 
This is an example of the [discrete logarithm problem](https://en.wikipedia.org/wiki/Discrete_logarithm#Properties), an important primitive that has various applications in cryptography. 

# Breaking the Discrete Log

Note that the underlying group $$ G = \mathbb{Z}^{\times}_p $$ that we are working in has order 

$$ \begin{align*} |G| &= \phi(p) = p-1 \\ &= 2 \cdot 3 \cdot 23 \cdot 1057807 \cdot 2132567 \cdot \underbrace{717\cdots661}_{\text{35 digits}}  \end{align*}$$

Since the order of the group factors into (relatively) small factors, we can use the [Pohlig-Hellman algorithm](https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm) to solve the discrete logarithm in this case. Briefly, this algorithm solves the discrete logarithm by solving the discrete logarithm modulo each of the factors and then using the Chinese Remainder Theorem to reconstruct the final solution.

# A Slight Problem
To solve the discrete logarithm modulo each of the factors, the Pohlig-Hellman algorithm defaults to using the [baby-step giant-step algorithm (or BSGS)](https://en.wikipedia.org/wiki/Baby-step_giant-step) which calculates the discrete logarithm in $$ \mathbb{Z}_n $$ in time $$ O(\sqrt{n}) $$.

Note that the largest prime factor of $$ |G| $$ is $$ p'=71765404858975364469794424368755661 $$ 
and therefore, the application of BSGS at this step will take approximately $$ \sqrt{p'} \approx 10^{17} $$.
Therefore, a simple application of BSGS would be intractable - we instead opt for an application of the [number field sieve](https://en.wikipedia.org/wiki/General_number_field_sieve) as implemented by [CADO-NFS](https://cado-nfs.gitlabpages.inria.fr)

```
$ ./cado-nfs.py -dlp -ell 71765404858975364469794424368755661 target=14930497059761124952555442014579540362281335175614 22341037975022202529922405433894518512705789413043 -t all --workdir //home/epistemologist/cado-nfs/_work_dir/
[...]
Info:Complete Factorization / Discrete logarithm: Total cpu/elapsed time for entire Discrete logarithm: 63.44/47.7757
Info:root: log(target) = 55489157388397020258158645620533604 mod ell
Info:root: logbase = 3313713443676474566743726004211879966369713471586                                                                        
Info:root: target = 14930497059761124952555442014579540362281335175614                                                                        
Info:root: log(target) = 55489157388397020258158645620533604 mod ell                                                                          
55489157388397020258158645620533604 
$ ./cado-nfs.py //home/epistemologist/cado-nfs/_work_dir/p50.parameters_snapshot.0 target=7
Info:root: logbase = 3313713443676474566743726004211879966369713471586
Info:root: target = 7
Info:root: log(target) = 25109992438216120939941563608435562 mod ell
25109992438216120939941563608435562
```

We therefore have that 

$$c = \frac{55489157388397020258158645620533604}{25109992438216120939941563608435562}  \mod{p'}$$

# Solving the Rest of the Discrete Logarithm

We now apply Pohlig-Hellman for the rest of the factors to find the value of $c$.

```py
from sage.groups.generic import bsgs

p = 0xf494eca63dcab7b47ac21158799ffcabca8f2c6b3
g = 7
h = 0xa3742a4abcb812e0c3664551dd3d6d2207aecb9be

exponents, moduli = [], []

for p_, e_ in list(factor(p-1)):
    g_i = pow(g, (p-1) // (p_^e_), p)
    h_i = pow(h, (p-1) // (p_^e_), p)
    if p_ < 10^10:
        exponents.append(bsgs(g_i, h_i, (0, p_^e_-1)))
    else:
        exponents.append(int(Integers(71765404858975364469794424368755661)(55489157388397020258158645620533604/25109992438216120939941563608435562)))
    moduli.append(p_^e_)
    
c = crt(exponents, moduli)
```

Running this code, we get the value of $$ c = 1068574207815876554047411521461868356487653669046 $$. We can verify this:

```py
>>> c = 1068574207815876554047411521461868356487653669046
>>> pow(7,c,0xf494eca63dcab7b47ac21158799ffcabca8f2c6b3)==0xa3742a4abcb812e0c3664551dd3d6d2207aecb9be
True
```

# Getting the Sequence

To get the original sequence, it is possible to reverse the code by hand as I did originally - however, it is also possible to pass this into Z3 and let it spit out the answer:

```py
from z3 import *

a0, d = Ints("a0 d")
b = [a0+k*d for k in range(5)]

c = 1
for d in b:
    c = c * (1<<32) + d

c_actual = 1068574207815876554047411521461868356487653669046

s = Solver()
s.add(c == c_actual)
for i in b:
    s.add(-(1<<32) < i)
    s.add(i < (1<<32))

s.check()
print(s.model()) # [d = 316804249, a0 = -1154709934]
```

With this, we get that our final sequence is `[-1154709934, -837905685, -521101436, -204297187, 112507062]`.
