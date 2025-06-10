---
layout: post
title:  "Some Notes on Integer Factorization"
categories: 
---

# Some Notes on Integer Factorization

## A Factorization From "The Art of Computer Programming"
We aim to factor the integer $$N = 2^{214} + 1$$. Using the identity $$4x^4+1 = (2x^2+2x+1)(2x^2-2x+1)$$, we can break $$N$$ up into two factors as 

$$ N = \underbrace{ (2^{107} - 2^{54} + 1)}_{A}\underbrace{(2^{107} + 2^{54} + 1)}_{B}$$

### Factoring A
By applying trial division for all primes $$p < 1000$$, we have that 

$$ A = 162259276829213345377179500806145 = 5 \cdot 857 \cdot p_{29} $$

where $$p_{29} = 37866809061660057264219253397$$ is a 29-digit prime [confirmed by the Miller-Rabin test](https://tio.run/##TY7PboQgGMTvPsVcmmiLUWARNfHgpU/QW9M0trJd0g0a6mbdp7cf9M/ulxCGGfgN82U5TE7m82W@bFuPDrwSQjVCV7VoBJdyp6TWXDeqLOuy4jtVFGl9XyudJZ5hpBc9Q5mcD/Zo4HEHga5D2Sag8SiKDiLqEQ9ET/aTxwDr8CwYJINi0Ayc06ITD7phEKQF7ZJ8qV9aRMRKbfN0Tgc2sj6Llt0Hl8AgblR9zlu8T26x7mTinVD5Gir94D5M6nOe/Xzvlrky8ce85QbavxnmzZvhMzrm@GWu2eytW1I8DuTiyjGrXdIs@U2f/InCbfsG).

### Factoring B
We first try to divide out small primes from $$B$$ - we have that
$$B = 162259276829213381405976519770113 = 843589 \cdot c_{27}$$
where $$c_{27} = 192343993140277293096491917$$. We know that $$c_{27}$$ is composite as we have that $$3^{c_{27}-1} \not\equiv 1 \bmod(c_{27}) $$ which contradicts Fermat's little theorem. [The Pollard Rho algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm) [discovers](https://tio.run/##VVBbbsMgEPznFNtIlSAPyRgrFpZ6BX/lAm7ANpJtEMEqnN4B0tTqfg0z7OzsmuBGvbCLCSZsW2/1DHPnRlCz0dbBcBcok8pJ67SeHm/lrtfFIdTCF1BesopxzmhVlHVdclbwa8UppzVCQvZgR41b0iCI5WNDmVGIyGckkklGQ0RTN3@LDnwD2B/9iZLPNmspQpSL/PgZ1SThZlf5sn1bD9iTPyJkYsCB7FyaFbfCOFx8dD5Du2uqj/JHzLJ7pjJWLQ4f0vzmcM45yL8PVrrVLiB2oxT1lJZCr@bfEwDZtic) that $$c_{27}$$ factors into two primes:
$$c_{27} = 8174912477117 \cdot 23528569104401$$


## A Modern Sized Factorization A La Knuth
In a similar approach, we attempt to factor $$N = 2^{432}+1$$, a composite with 131 digits. 

We can factor the polynomial $$p(x) = x^{432} + 1$$ to get an algebraic factorization of $$N$$. We have that:

$$N = \underbrace{ (2^{16} + 1) }_{65537}\cdot \underbrace{ (2^{32} - 2^{16} + 1)}_{C_1} \; \cdot \; \underbrace{  (2^{96} - 2^{48} + 1) }_{C_2}  \cdot \underbrace{ (2^{288} - 2^{144} + 1)}_{C_3}$$

where $$C_1 = 4294901761 = 193 \cdot 22253377$$ factors via trial division. 

$C_2 = 79228162514264056118567239681$ is a 29-digit composite with relatively small factors; it can be factored with the [GNU `factor` utility](https://www.gnu.org/software/coreutils/manual/html_node/factor-invocation.html#factor-invocation):

```sh
$ factor 79228162514264056118567239681
79228162514264056118567239681: 1153 6337 38941695937 278452876033
```
and therefore, we have that $$C_2 = 1153\cdot 6337 \cdot  38941695937 \cdot 278452876033$$


With this luck, we hope that $$C_3$$ also has a small factor - we use the utility [`gmp-ecm`](https://www.rieselprime.de/ziki/GMP-ECM), an implementation of the [elliptic curve factoring method](https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization) which is useful for finding small to medium sized prime factors)

```txt
$ python3 -c "print(2**288 - 2**144 + 1)" | ecm 100000
GMP-ECM 7.0.5 [configured with GMP 6.3.0, --enable-asm-redc] [ECM]
Input number is 497323236409786642155382248146820840100456128496602518909840835357441224364171869552641 (87 digits)
Using B1=100000, B2=40868532, polynomial x^2, sigma=1:2036516883
Step 1 took 115ms
********** Factor found in step 1: 68016300334849
Found prime factor of 14 digits: 68016300334849
Prime cofactor 7311824282729722035859309520826138827918372038863677678267962656996811009 has 73 digits
```


While we broke up the factorization into stages, the entire factorization is small enough that it can be done in software packages such as `gp` and Sagemath:

```txt
$ echo  "factorint(2^432+1)" | gp -q 

[                                                                      193 1]

[                                                                     1153 1]

[                                                                     6337 1]

[                                                                    65537 1]

[                                                                 22253377 1]

[                                                              38941695937 1]

[                                                             278452876033 1]

[                                                           68016300334849 1]

[7311824282729722035859309520826138827918372038863677678267962656996811009 1]
```


## An Unbreakable Cipher

In the August 1977 issue of Scientific American, [Martin Gardener](https://www.jstor.org/stable/24954008?seq=1) presented an instance oF [RSA](https://en.wikipedia.org/wiki/RSA_cryptosystem) that was prediced to "take millions of years to break" that involved factoring a 129 digit number.

### Installing CADO-NFS
CADO-NFS is an open source implementation of the [ number field sieve algorithm ](https://en.wikipedia.org/wiki/General_number_field_sieve); we can install and build it from source:
```sh
sudo apt install make cmake g++ gcc libgmp3-dev python3 python3-flask git # Install necessary libraries
git clone https://gitlab.inria.fr/cado-nfs/cado-nfs.git && cd cado-nfs # Get source code for CADO-NFS
make -j8 # Change to however many cores you have
./cado-nfs.py $(python3 -c "print(10**57+1)") # Run a test factorization
```
### Breaking the Cipher

After installing, we can use CADO-NFS to factor the 129-digit composite on a 16 core AWS machine in around 2.5 hours:

```
$ time ./cado-nfs.py 114381625757888867669235779976146612010218296721242362562561842935706935245733897830597123563958705058989075147599290026879543541
[...]
Info:Complete Factorization / Discrete logarithm: Total cpu/elapsed time for entire Complete Factorization 243582/9509.24 [02:38:29]
3490529510847650949147849619903898133417764638493387843990820577 32769132993266709549961988190834461413177642967992942539798288533

real    149m9.555s
user    346m20.767s
sys     21m6.648s
```

We can then use these primes to decrypt the original ciphertext given in the original issue of Scientific American:
```python
N = 114381625757888867669235779976146612010218296721242362562561842935706935245733897830597123563958705058989075147599290026879543541
p = 3490529510847650949147849619903898133417764638493387843990820577
assert N % p == 0
q = N // p
c = int("""
9686 9613 7546 2206
1477 1409 2225 4355
8829 0575 9991 1245
7431 9874 6951 2093
0816 2982 2514 5708
3569 3147 6622 8839
8962 8013 3919 9055
1829 9451 5781 5154
""".replace("\n", "").replace(" ", ""))
e = 9007

from string import ascii_lowercase
ALPHABET = " " + ascii_lowercase
pt = str( pow(c, pow(e, -1, (p-1)*(q-1)), N) )
print( "".join( [ ALPHABET[ int( pt[i:i+2] ) ] for i in range(0, len(pt), 2) ] ) )
```

```
$ python3 decrypt_message.py 
the magic words are squeamish ossifrage
```

### A Historical Note
The original factorization was completed in 1993-1994 via a team lead by [Atkins et al](http://web.mit.edu/warlord/www/rsa129.ps)
 - the factoriztaion used the quadratic sieve algorithm instead of the number field sieve used above and took 825 mips years
 - the majority of this time was spent on the sieving step (~220 days over 1600 machines), as in our calculation
 - linear algebra took around 100 hours
 
#### Time to Factor 512-Bit Number
The technical paper mentioned above mentions that factoring a 512-bit number would take "a couple months" at 500,000 mips - we therefore have a total instruction count of 

$$500\cdot10^9 \text{ ops per second} \cdot 4 \text{ months} \approx 5.26 \cdot 10^{18} \text{ ops} $$

We can compare this instruction count to a [factorization of CADO-NFS on a 512-bit modulus on modern hardware](https://yurichev.com/news/20220210_RSA/):

$$\frac{3.6 \cdot 10^9 \text{ ops/second}}{\text{thread}} \cdot 12 \text{ threads } \cdot 4 \text{ days} \approx  1.5 \cdot 10^{16} \text{ ops} $$
