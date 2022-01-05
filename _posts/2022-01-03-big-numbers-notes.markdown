---
layout: post
title:  "Notes on Littlewood's 'Large Numbers' and Related Thoughts"
date:   2022-01-03 08:28:12 -0600
categories: 
---

In his essay titled ["Large Numbers"](https://www.jstor.org/stable/3609933?seq=1#metadata_info_tab_contents)(published as part of his larger collection of essays titled *A Mathematical Miscellany*), John Littlewood writes on the various applications of large numbers. In this post, we give a more modern take on some of his examples as well as give some new ones.


## Human Scale Numbers

### Scales of Measurement

In §7, Littlewood discusses the size of various scales of measurement (e.g. the range of sound from just perceptible to just tolerable - Littlewood mentions a figure of $$ 10^{12} $$, however the range between [0 decibels and 120 decibels](https://en.wikipedia.org/wiki/Sound_pressure#Examples) represents a range of $$ 10^6 $$).

Other examples include
 - **luminance**: $$ 10^{-6} \frac{\text{cd}}{\text{m}^2} $$ (threshold of vision) to $$  10^{6} \frac{\text{cd}}{\text{m}^2} $$ (incandescent lamp) - *range of $$ 10^{12} $$*
 - **length**: $$ 10^{-3} m$$ to $$ 10^{3} m$$ - *range of $$ 10^6 $$*
 - **time**: $$ 10^{-2} s $$ (human reaction time) to $$ 10^9 s$$ (human lifetime) - *range of $$ 10^{11} $$* 
 - **mass**: $$ 10^{-6} \text{kg} $$ (milligram) to $$ 10^3 \text{kg} $$ (metric ton) - *range of $$ 10^9 $$*

### Coincidences / Probabilities

#### Will West Case

Littlewood writes:
>Dorothy Sayers in Unpopular Opinions, cites the case of two negroes, each named Will West, confined simultaneously in Leavenworth Penitentiary, U.S.A. (in 1903), and with the same Bertillon measurements. (Is this really credible ?) 

We estimate this probability, let $$P$$ be the probability that two people share the same facial features, Bertillon measurements and name. We have:

$$ P \approx \underbrace{\text{Prob(same body measurements)}}_{:=P_b} $$
$$ \cdot \;  \underbrace{\text{Prob(same name)}}_{:=P_n}$$

##### Estimating $$ P_n $$

Let $$n_1, n_2$$ be the names of two arbitrary people. We have the crude approximation:

$$ P(n_1 = n_2) = \sum_{\text{names }n} P(n_1 = n) P(n_2 = n) \\ = \sum_{\text{names } n} P(n)^2 $$

Using data from [U.S. Social Security Adminbistration](https://www.ssa.gov/OACT/babynames/) and using the top 1000 names from the years 1880 and 2010 (assuming that 50% of births are male and 50% of births are female), we [get estimates](https://gist.github.com/epistemologist/12cff945e6b8c918055f5585fc373c7f) of probabilities of 1.8% and 0.21% respectively. Assuming that choosing first names and last names are independent, we estimate that $$  P_n \approx 10^{-2} \cdot 10^{-2} = 10^{-4} $$

##### Estimating $$P_b$$

Littlewood mentions that the two Wests' shared Bermillon measurements - these are a system of 11 different biometric measurements that were taken of prisoners as means of identification. In the original case, the measurements of the two Wests' were [not exactly the same, but remarkably similar](https://web.archive.org/web/20120505160629/opinionator.blogs.nytimes.com/2012/05/01/whats-in-a-name-part-2/). 

Bertillon [himself predicted](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1740-9713.2014.00739.x) that the probability of two people sharing all 11 measurements was $$P_b \approx 1/4^{32} = 1/268435456$$. However, [more recent studies](https://link.springer.com/article/10.1007%2Fs00414-015-1158-6) in forensic anthropometry place the probability of two people sharing 8 different body measurement traits at $$ P_b \approx 10^{-20} $$.

We therefore make the estimate $$ P_b \approx 10^{-12} $$ and therefore, $$ P \approx 10^{-12} \cdot 10^{-4} = 10^{-16} $$.

##### Wrapping Things Up

We now calculate the probability that for any one person being incarcerated that they share the same features to someone in the same prison that they are being held similarly to the two Will Wests.

Assuming an average prison houses an average of [250 inmates](https://www.americanjail.org/jail-statistics), the probability that a person who is imprisoned meets their doppelganger in the same prison is $$p = 1 - (1-P)^{250} $$. 

Therefore, if $n$ people are incarcerated, we have that the probability that any one person will find their doppelganger is $$ p' = 1 - (1-p)^n = 1 - (1-P)^{250n} $$.
Considering that the U.S. made [10 million arrests in 2019](https://www.ojjdp.gov/ojstatbb/crime/ucr.asp?table_in=2), we estimate $$ N \approx 5 \cdot 10^6 $$. We therefore have that:

$$ p' = 1 - (1-P)^{N} \\ = 1 - \left( 1 - NP + \frac{N(N-1)}{2} P^2 - \cdots \right) \\ = NP + O(P^2) $$ 

Plugging in $$N = 250 \cdot 5 \cdot 10^6, P = 10^{-16}$$, we get that $$p' \approx 10^{-7}$$. 

This probability is reasonable as [various doppelgangers have been found in real life](https://en.wikipedia.org/wiki/Doppelg%C3%A4nger#Examples_in_real_life) and so we expect it to be within the realm of possibility.

#### A Perfect Bridge Deal

Littlewood writes:
> A report of holding 13 of a suit at Bridge used to be an annual event. The chance of this in a given deal is 2.4 * 10^-9 ; if we suppose that 2 * 10^6 people in England each play an average of 30 hands a week the probability is of the right order. I confess that I used to suppose that Bridge hands were not random, on account of inadequate shuffling ; Borel's book on Bridge, however, shows that since the distribution within the separate hands is irrelevant the usual procedure of shuffling is adequate. (There is a marked difference where games of Patience are concerned : to destroy all organisation far more shuffling is necessary than one would naturally suppose ; I learned this from experience during a period of addiction, and have since compared notes with others.)

Note that Littlewood gives a probability for holding 13 of a suit to be $$2.4 \cdot 10^{-9} $$; however, the [actual probability](https://stats.stackexchange.com/a/288235) seems to be $$ \approx 2.5 \cdot 10^{-11} $$. Using Littlewood's assumptions, we have that the expected rate of perfect bridge suits is $$ 2.5 \cdot 10^{-11} \cdot \frac{\text{30 hands}}{\text{week}} \cdot 2 \cdot 10^6 \; \text{players} = \frac{0.0015 \text{ hands}}{\text{week}} \approx \frac{7.8 \text{ hands}}{\text{century}}$$

Meanwhile, the probability of all 4 players receiving perfect deals is $$ \approx 4.5 \cdot 10^{-28} $$, a much rarer occurance - this is discussed in [this video by Matt Parker](https://www.youtube.com/watch?v=s9-b-QJZdVA).


### Computation

#### Factorizations

In §13, Littlewood discusses the scope of factorizations of large numbers, quoting [D.H. Lehmer](https://en.wikipedia.org/wiki/D._H._Lehmer) 

> Professor Lehmer further tells me that numbers up to 2.7 * 10^9 can be completely factorised in 40 minutes; up to 10^15 in a day; up to 10^20 in a week; finally up to 10^100, with some luck, in a year

The cutting edge of the largest integers we can factor now has well evolved since the time of Littlewood due to both advances in better hardware and the development of better algorithms.

For reference, I can factor numbers of the magnitude $$ 10^{40} $$ nearly instantly on my nearly 10-year old laptop:

```gp
? x = random(10^40); x
%1 = 8986749531375339305571207626880980158008
? factor(x)
%2 = 
[                      2 3]

[                      3 3]

[                 319519 1]

[             2078212777 1]

[62655931019629323359651 1]

? ##
  ***   last result computed in 12 ms.
```

In a 1977 issue of Scientific American, Martin Gardener [encoded a message with RSA](https://en.wikipedia.org/wiki/The_Magic_Words_are_Squeamish_Ossifrage) using a modulus on the order of $$ 10^{129} $$. The message was eventually decrypted by a distributed computation effort lead by [Atkins et al.](https://link.springer.com/chapter/10.1007%2FBFb0000440). 

As of writing, the largest number that has been factored has [250 decimal digits](https://en.wikipedia.org/wiki/RSA_numbers#RSA-250) and was factored using 2700 CPU core-years with the Number Field Sieve algorithm.

#### Training AI Models

With the rise of large AI models like GPT-2/3 with billions of parameters, the amount of computation to train these models has risen exponentially in recent years. [This blog post from OpenAI](https://openai.com/blog/ai-and-compute/) shows a *3.4 month doubling time* in the increase of AI training. Note that models like AlphaGo and AlphaGoZero exceed $$ 10^2 $$ petaflop/s-days, or around $$ 10^{22} $$  operations $$ \approx 0.1 \text{mol} $$

#### Bitcoin Mining

Looking at the [hashrate for Bitcoin]( https://www.blockchain.com/charts/hash-rate ), as of January 2022 the Bitcoin network is performing $$ \approx 177 \cdot 10^{18} $$ SHA 256 hashes **per second**. Assuming that each hash is equivalent [to around 1000 arithmetic operations](https://bitcoin.stackexchange.com/questions/1293/how-many-integer-operations-on-a-gpu-are-necessary-for-one-hash), we have that the Bitcoin network can perform $$ 10^{23} \frac{\text{operations}}{s}
\approx 1 \frac{\text{mol}}{s}$$. Throughout its entire existence, the Bitcoin network has performed $$ 10^{28} $$ hashes or around $$ 10^{31} \text{operations} \approx 10^8 \text{mol} $$.

#### Breaking Cryptographic Schemes

##### Distributed.net

[Distributed.net](https://en.wikipedia.org/wiki/Distributed.net) is a platform designed to organize large-scale distributed computation projects. One of these projects was to break messages encrypted under [RC5](https://en.wikipedia.org/wiki/RC5) with varying key lengths. The projects' statuses are as follows:

 - [cracking a key of 56 bits](https://stats.distributed.net/projects.php?project_id=3): $$ 26 \cdot 10^{15} $$ keys attempted over 193 days, completed
 - [cracking a key of 64 bits](https://stats.distributed.net/projects.php?project_id=5): $$ 15 \cdot 10^{18} $$ keys attempted over 1726 days, completed
 - [cracking a key of 72 bits](https://stats.distributed.net/projects.php?project_id=8): $$ 390 \cdot 10^{18} $$ keys attempted over 6972 days, incomplete

## Nonhuman Scale Numbers

### A Mouse Surviving in Hell

We expand on some of the ommited details in Littlewood's calculation on the probability that a mouse can survive a week in Hell.

Following Littlewood's notation, let
 - $$ T_0 $$ be average room temperature (note all temperature units are measured [absolutely](https://en.wikipedia.org/wiki/Thermodynamic_temperature) - i.e. absolute zero is at $$T = 0$$)
 - $$ T_H $$ be the temperature of Hell (we assume that $$ T_h \gg T_0 $$)
 - $$ \mu = k n_0 $$ be the number of particles in the mouse (here, $$ n_0 $$ is [Avogadro's number](https://en.wikipedia.org/wiki/Avogadro%27s_number), so $$ k $$ represents the number of moles)

We assume, as Littlewood does, that "...we should treat the problem as classical, and suppose that the molecules and densities are terrestial." - let the particles in Hell (or somewhere extremely hot if you prefer) follow a [Maxwell-Boltzmann distribution](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution)


Littlewood also defines two other variables, $$ c_0, c_H $$ or the average speed of particles at temperatures of $$ T_0 $$ and $$ T_H $$ respectively. We can describe this by calculating the [mean velocity](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Typical_speeds) at various temperatures, giving $$ c_0 = \sqrt{\frac{8 k_B T_0}{ \pi m}}, c_H = \sqrt{\frac{8 k_B T_H}{\pi m}} $$.

Next, for any given particle at temperature $$ T =  T_H $$, the probability $$ p $$ that its speed does not exceed $$ c_0 $$ can be found by integrating [the distribution function of the particle's speed](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_speed) as follows:

$$ p = \int_0^{c_0} f(v) = \sqrt{\frac{2}{\pi}} \left(\frac{m}{k T_H}\right)^{3/2} \int_0^{c_0}  v^2 e^{-\frac{mv^2}{2k T_H}} dv \\ = \sqrt{\frac{2}{\pi}} \left(\frac{m}{k T_H}\right)^{3/2} \int_0^{c_0} \left( v^2 + O(1/T_H) \right) dv \\ \sim T_H^{-3/2} c_0^3 \sim (T_0 / T_H)^{3/2}  $$

Littlewood then states that the probability that most of the particles in the mouse have $$ c \le c_0 $$ is of the order of $$ p^\mu $$. He then breaks up the time interval of a week into time intervals of length $$ \tau $$ similar to the [mean free time](https://en.wikipedia.org/wiki/Mean_free_time) or average time between molecular collisions.

Here, we simply take Littlewood's word [TODO: actually derive this] and say that $$ \tau_H \sim \frac{n_0^{-1/3}}{c_0} \sqrt{\frac{T_0}{T_H}} $$ and so the total number of time periods in a week is $$ \nu = w / \tau_H $$ 

We therefore have that the total probability of survival is $$ p = 1/C$$, where

$$ C = (p^{-\mu})^\nu = \left(\sqrt{\frac{T_H}{T_0}}\right)^{\frac{3}{2} c_0 k w n_0^{4/3} \sqrt{\frac{T_H}{T_0}}} $$ 

Let's now substitute actual numbers into this equation (here, we depart from the values that Littlewood uses as his values are just to get a nice final result):
 - we set room temperature at $$ T_0 \approx 300 \text{K} $$ or around $$ 80^{\circ}\text{F} $$
 - the temperature of Hell is difficult to estimate - while Biblical accounts such as one in [Revelations 21:8](https://www.biblegateway.com/passage/?search=Revelation+21) describe a "fiery lake of burning sulfur" and therefore a temperature less than the boiling point of sulfur at [around 700 K](https://en.wikipedia.org/wiki/Sulfur), we take this as metaphor and instead take Littlewood's value of $$ T_H   \approx 10^{12} \text{K} $$. For comparison, [similar temperatures](https://en.wikipedia.org/wiki/Orders_of_magnitude_(temperature)) can be seen in the formation of stars and microseconds after the Big Bang.
 - we use a value of $$ k = 10^2 $$; for reference, $$ 10^2 \text{mol} $$ of water weighs around $$ 1.8 \text{kg} $$ 
 - we have that $$ w = 7 \text{days} \approx 6 \cdot 10^5 \text{s} $$
 - we [can calculate](https://www.wolframalpha.com/input/?i=sqrt%28+%288+*+Boltzmann%27s+constant+*+300+Kelvin%29+%2F+%28pi+*+18+grams+%2F+Avogadro%27s+number%29+%29+) $$ c_0 $$ from our earlier definition: $$ c_0 = \sqrt{\frac{8 \cdot (1.4 \cdot 10^{-23} \frac{J}{K} ) \cdot 300K }{\pi \left( \frac{18 g}{6 \cdot 10^{23}} \right) }} \approx 600 \frac{m}{s}$$

We can now explicitly calculate C:
```py
In [1]: from math import sqrt, log 
   ...: T_0 = 300 
   ...: T_H = 1e12 
   ...: k = 1e2 
   ...: w = 6e5 
   ...: c_0 = 600 
   ...: n_0 = 6e23 
   ...: exp = 3 / 2 * c_0 * k * w * n_0**(4/3) * sqrt(T_H / T_0) 
   ...: base = sqrt(T_H / T_0) 
   ...: print(f"C = 10^({exp*log(base) / log(10)})")                                                                          
C = 10^(7.512302344442364e+47)
```

#### Other Examples
This is not the only example of large numbers arising from the study of physical phenomena:
 - [Here](http://www.fpx.de/fp/Fun/Googolplex/GetAGoogol.html), it is described that a black hole of mass $$ 6.14 \cdot 10^{41} \text{kg} \approx 3 \cdot 10^{11} $$ times the mass of the sun has a dimensionless entropy of $$ 10^{100} $$ and therefore has $$ e^{10^{100}} \approx 10^{10^{99.6}}$$ different macroscopic states
 - [In this paper](https://arxiv.org/pdf/1712.08465.pdf), the probability that a human teleports due to random quantum fluctuations is calculated to be approximately $$ 10^{-4.5 \cdot 10^{29}} $$ 
