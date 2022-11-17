---
title: 'Numerical Methods for Conservation Laws (Part 1)'
subtitle: 'An introduction to solving nonlinear conservation laws.'
summary: An introduction to solving nonlinear conservation laws.
authors:
- admin
tags:
- Numerical PDE
- Conservation Laws
- Finite Volume Method
categories: []
date: "2019-05-15T00:00:00Z"
lastmod: "2019-05-25T00:00:00Z"
featured: false
draft: false

# Featured image
# To use, add an image named `featured.jpg/png` to your page's folder.
# Focal point options: Smart, Center, TopLeft, Top, TopRight, Left, Right, BottomLeft, Bottom, BottomRight
image:
  caption: 'Image credit: **Justin Dong**'
  focal_point: ""
  preview_only: false

# Projects (optional).
#   Associate this post with one or more of your projects.
#   Simply enter your project's folder or file name without extension.
#   E.g. `projects = ["internal-project"]` references `content/project/deep-learning/index.md`.
#   Otherwise, set `projects = []`.
projects: []

# Set captions for image gallery.
gallery_item:
- album: gallery
  caption: Default
  image: theme-default.png
- album: gallery
  caption: Ocean
  image: theme-ocean.png
- album: gallery
  caption: Forest
  image: theme-forest.png
- album: gallery
  caption: Dark
  image: theme-dark.png
- album: gallery
  caption: Apogee
  image: theme-apogee.png
- album: gallery
  caption: 1950s
  image: theme-1950s.png
- album: gallery
  caption: Coffee theme with Playfair font
  image: theme-coffee-playfair.png
- album: gallery
  caption: Cupcake
  image: theme-cupcake.png
---

## **1. Conservation Laws**

In this post, we'll take a look at conservation laws, the contexts in which they arise in nature, and some of the numerical methods used for solving them. In one spatial-dimension, conservation laws take the general form 
{{< math >}}
$$
\begin{equation}
  \begin{cases}
    u_{t} + f(u)_{x} = 0 &\text{in}\;\mathbb{R} \times [0,\infty)\\\\\\
    u(x,0) = u_{0}(x) &\text{on}\;\mathbb{R} \times \{t=0\}.
  \end{cases}
  \label{eq:conslaw}
\end{equation}
$$
{{< /math >}}

Here, we consider the spatial domain to be the entire real line (thus, we have no boundary conditions). The function $f(u)$ -- typically called the *flux function* -- is sufficiently smooth. As their name suggests, conservation laws preserve mass. We often think of the quantity $u$ as the density of some fluid, in which case $\int\_{\mathbb{R}} u(x,t)dx$ may be viewed as the mass of fluid. If we assume that $u$ is compactly supported in $x$ (i.e. it is nonzero outside of some compact set in $x$, a reasonable assumption as the fluid should have finite mass) and $f(0) = 0$, then taking derivatives over time yields
$$
\begin{align}
  \frac{\partial}{\partial t} \int\_{-\infty}^{\infty} u(x,t)\;dx &= \int\_{-\infty}^{\infty} \frac{\partial u}{\partial t}\;dx = -\int\_{-\infty}^{\infty} \frac{\partial f(u)}{\partial x}\;dx\notag\\\\\\
  &= f(u(-\infty)) - f(u(\infty)) = 0.\notag
\end{align}
$$

The above computation implies that $\int\_{\mathbb{R}} u(x,t)\;dx = \int\_{\mathbb{R}} u(x,0)\;dx$: the mass of $u$ is conserved over time. 

To gain some intuition for how solutions of this PDE might behave, we'll first consider the simplest conservation law there is: the linear advection equation. Taking $f(u) = u$, it is evident that our solution is given by $u(x,t) = u\_{0}(x-t)$. The PDE "transports" the initial profile in the $x$ direction. In particular, we note that the solution to the *linear* conservation law is exactly as smooth as our initial data. That is, if $u_{0} \in C^{1}\_{x}(\mathbb{R})$ then $u \in C\_{x}^{1}(\mathbb{R})$ as well. We haven't said anything particularly illuminating thus far and indeed, the linear advection equation isn't particularly interesting to begin with. 

Next, let's consider a simple nonlinear conservation law:
$$
\begin{cases}
  u\_{t} + uu\_{x} = 0 &\text{in}\;\mathbb{R} \times [0,\infty)\\\\\\
  u(x,0) = \sin{x} &\text{on}\;\mathbb{R} \times \\{t=0\\}.
\end{cases}
$$

The initial condition is as smooth as possible now: $\sin{x} \in C\_{x}^{\infty}(\mathbb{R})$. But what can we say about $u(x,t)$? Should we still expect it to be smooth in $x$? It turns out that the answer is a resounding *no*. In fact, solutions to the above nonlinear PDE (known as *Burgers equation*) are not even $C\_{x}^{0}(\mathbb{R})$! The exact mathematical cause of this phenomena is due to the crossing of something called characteristic curves, but we will first consider a more intuitive explanation.

We may view Burgers equation as an advection equation in which the speed of propagation is equal to $u$ itself. Our initial sine profile thus moves with variable speed in $x$. For $0 \leqslant x \leqslant \pi$, the sine wave is positive and the solution travels forward. For $\pi < x \leqslant 2\pi$, the sine wave is negative and the solution travels backwards. The animation below shows how the solution of Burgers equation evolves over time, and we see the development of a genuine discontinuity at $x=\pi$. 
{{< figure src="gallery/burgers.gif" caption="**Figure 1:** Solution of Burgers equation with sine initial data.">}}

Shocks such as the one above occur naturally in many settings, for instance high-speed compressible flows and aeroacoustics. However, a natural question we might ask next is how discontinuities and shocks fit with our notion of the classical solution of a PDE. For instance, the linear advection equation contains a first-order time derivative and a first-order spatial derivative. Naturally, we should require that $u \in C\_{x,t}^{1}(\mathbb{R} \times [0,\infty))$ -- our solution should be once differentiable in space and time. Consider, though, the linear advection equation with $u(x,0) = 𝟙\_{\\{x \leqslant 0\\}}$. You'll probably agree with me that the only possible solution to this problem is $u(x,t) = 𝟙\_{\\{x \leqslant t\\}}$, and yet this "solution" is not even continuous let alone differentiable. Before moving on to numerical solvers, we will spend some time developing alternative notions of PDE solutions that are capable of admitting discontinuities and shocks.


## **1.2 Weak Solutions of Conservation Laws**
**Definition 1.** *We call $u(x,t)$ a* **weak solution** *of the conservation law {{< math >}}$\eqref{eq:conslaw}${{< /math >}} if*
$$
\begin{equation} 
  \int\_{0}^{\infty}\int\_{-\infty}^{\infty} u\varphi\_{t} + f(u)\varphi\_{x}\;dxdt = -\int\_{-\infty}^{\infty} u(x,0)\varphi(x,0)\;dx
  \label{eq:weaksol}
\end{equation}
$$

*for all $\varphi \in C\_{c}^{\infty}(\mathbb{R} \times \mathbb{R})$.*

We refer to $\varphi$ as a test function. It is smooth in both space and time and compactly supported, i.e. $\varphi$ is zero outside of some compact set in space and time. Let's break down what's actually going on in the definition of the weak solution. The most notable feature is that if $u$ is a weak solution, it does *not* have to be differentiable in space or time since the definition only contains derivatives of $\varphi$!

So how do we arrive at this definition? First, we multiply the entire conservation law by $\varphi$:
$$
u\_{t}\varphi + f(u)\_{x}\varphi = 0.
$$

Then, we integrate in time over $[0,\infty)$ and in space over the entire real line:
$$
\int\_{0}^{\infty}\int\_{-\infty}^{\infty} u\_{t}\varphi + f(u)\_{x}\varphi\;dxdt = 0
$$

The last step is to integrate the first term by parts in time and the second term by parts in space:
$$
\begin{align}
  \int\_{0}^{\infty} \int\_{-\infty}^{\infty} u\_{t}\varphi dxdt &= -\int\_{0}^{\infty} \int\_{-\infty}^{\infty} u\varphi\_{t}dxdt + \int\_{-\infty}^{\infty} u(x,0)\varphi(x,0)dx\\\\\\
  \int\_{0}^{\infty} \int\_{-\infty}^{\infty} f(u)\_{x}\varphi\;dxdt &= -\int\_{0}^{\infty} \int\_{-\infty}^{\infty} f(u)\varphi\_{x}\;dxdt
\end{align}
$$

Most of the boundary terms vanish because $\varphi$ has compact support: we have $\varphi(\pm \infty,t) = \varphi(x,\pm\infty) = 0$. Note that in integrating by parts, we pass the derivatives from $u$ onto the test function. Altogether, we obtain
$$
\int\_{0}^{\infty}\int\_{-\infty}^{\infty} u\varphi\_{t} + f(u)\varphi\_{x}\;dxdt = -\int\_{-\infty}^{\infty} u(x,0)\varphi(x,0)\;dx
$$

But *why* do we do this? The idea is that rather than consider pointwise values of $u(x,t)$, we consider averaged values of $u(x,t)$ against test functions. If you've taken a functional analysis course, you can view the integration against a test function as the evaluation of a distribution induced by $u\_{t}$ and $f(u)\_{x}$. We note that all classical solutions ($C\_{x,t}^{1}$ solutions that satisfy the PDE in a pointwise sense) are weak solutions but the converse is certainly not true.

It is straightforward to verify that $𝟙\_{\\{x\leqslant t\\}}$ is a weak solution of the linear advection equation with $u(x,0) = 𝟙\_{\\{x \leqslant 0\\}}$. However, the definition $\eqref{eq:weaksol}$ of quite cumbersome to work with and we would like to develop a more convenient way to verify whether we have a weak solution or not. Notice that $u(x,t) = 𝟙\_{\\{x\leqslant t\\}}$ as well as the solution of Burgers equation in Figure 1 are piecewise smooth and only discontinuous at a single point that may (in the case of linear advection) or may not (in the case of Burgers equation) change with time. In this case, we can say much more about the structure of the weak solution.

**Theorem 1. (Rankine-Hugoniot)** *Suppose the solution of $\eqref{eq:conslaw}$ is piecewise smooth and contains a discontinuity along the curve $x(t)$. Then $u$ is a weak solution of $\eqref{eq:conslaw}$ if and only if
$$
\begin{equation}
  x'(t) = \frac{f(u^{-}) - f(u^{+})}{u^{-} - u^{+}},
  \label{eq:RH}
\end{equation}
$$
where $x'(t)$ is the speed of the discontinuity and $u^{\pm}$ are the values of $u$ along each side of the shock.*

Condition $\eqref{eq:RH}$ is referred to as the **Rankine-Hugoniot jump condition**. 

**Example 1.** *Let's return to the linear advection equation and use the jump condition to verify that $𝟙\_{\\{x\leqslant t\\}}$ is a weak solution. We have $x'(t) = 1$, $u^{-} = 1$, and $u^{+} = 0$. Then
$$
\frac{f(u^{-}) - f(u^{+})}{u^{-} - u^{+}} = \frac{u^{-} - u^{+}}{u^{-} - u^{+}} = 1
$$
and the jump condition is satisfied.*

Thus, the Rankine-Hugoniot condition is a handy way for us to verify whether something is a weak solution. It turns out that solutions to conservation laws often have the structure required by Theorem 1 -- that is, they are piecewise smooth except along finitely many curves of discontinuity. 

By relaxing the definition of the solution, we are able to admit much less smooth solutions to $\eqref{eq:conslaw}$. However, we will see weak solutions need not be unique! Consider Burgers equation with the initial condition $u(x,0) = -𝟙\_{\\{x\leqslant 0\\}} + 𝟙\_{\\{x>0\\}}$. It is straightforward to verify that $u(x,t) = -𝟙\_{\\{x\leqslant 0\\}} + 𝟙\_{\\{x>0\\}}$ as well as
$$
u(x,t) = \begin{cases}
  -1 &x \leqslant -t\\\\\\
  \frac{x}{t} &-t < x \leqslant t\\\\\\
  1 &x > t
\end{cases}
$$
are both weak solutions. The first weak solution does not change at all with time, so our intuition might tell us that this solution does not make much physical sense. But what about the second weak solution? At the very least, we need to establish more stringent criteria for our weak solutions in order to pick out the physically relevant solution. 


## **1.3 Entropy Solutions of Conservation Laws**

Consider the slightly modified PDE given by
$$
\begin{equation}
  \begin{cases}
    u^{(\varepsilon)}\_{t} + f(u^{(\varepsilon)})\_{x} = \varepsilon u^{(\varepsilon)}\_{xx} &\text{in}\;\mathbb{R} \times (0,\infty)\\\\\\
    u^{(\varepsilon)}(x,0) = u\_{0}(x) &\text{on}\;\mathbb{R} \times \\{t=0\\}.
    \label{eq:viscosityPDE}
  \end{cases}
\end{equation}
$$

This PDE is *parabolic* and it turns out that the solution $u^{(\varepsilon)}$, if it exists, is smooth and unique. 

**Definition 2.** *The entropy solution of $\eqref{eq:conslaw}$ is defined as 
$$
u(x,t) = \lim\_{\varepsilon \to 0} u^{(\varepsilon)}(x,t).
$$
This limit, if it exists, is unique and is a weak solution of  $\eqref{eq:conslaw}$.*

Again, this definition is quite difficult to work with, even moreso than the original definition of the weak solution, and we need to establish a more efficient way of verifying whether we have the entropy solution or not. There are alternative definitions based on the notions of *entropy flux* and *entropy flux functions*, but we'll skip right to the chase here. Thanks to the work of Oleinik and Lax, we have a quick way to identify entropy solutions. 

**Theorem 2. (Oleinik Entropy Condition)** *Suppose the solution of $\eqref{eq:conslaw}$ is piecewise smooth and contains a discontinuity along the curve $x(t)$. Then $u$ is the entropy solution if and only if
$$
\begin{equation}
  \frac{f(u) - f(u^{-})}{u - u^{-}} \geqslant \frac{f(u^{+}) - f(u^{-})}{u^{+} - u^{-}} \geqslant \frac{f(u^{+}) - f(u)}{u^{+} - u}
  \label{eq:oleinik}
\end{equation}
$$
for all $u$ between $u^{-}$ and $u^{+}$.*

Again, $u^{-}$ and $u^{+}$ denote the values of $u$ on each side of the shock $x(t)$. However, we can simplify $\eqref{eq:oleinik}$ even further if we have a *convex* conservation law, which is simply the case when $f(u)$ is convex. For Burgers equation, we have $f(u) = u^{2}/2$, which is indeed convex. The Lax Entropy condition tells us the following.

**Theorem 3. (Lax Entropy Condition)** *Suppose the solution of $\eqref{eq:conslaw}$ is piecewise smooth and contains a discontinuity along the curve $x(t)$. Suppose also that $f$ is convex. Then $u$ is the entropy solution if and only if
$$
\begin{equation}
  f'(u^{-}) \geqslant \frac{f(u^{+}) - f(u^{-})}{u^{+} - u^{-}} \geqslant f'(u^{+}).
  \label{eq:lax}
\end{equation}
$$
In particular, $\eqref{eq:lax}$ is equivalent to requiring that $u^{-} > u^{+}$.*

The Lax Entropy condition tells us that for convex conservation laws, the value of the weak solution on the left side of the shock ($u^{-}$) must be larger than the value of the solution on the right side of the shock ($u^{+}$) in order for the weak solution to be an entropy solution. The physical intuition behind the entropy solution is related to the entropy of a fluid in gas dynamics: in smooth flows, the entropy remains constant along particle paths, and if the particle crosses a shock, the entropy may only jump to a *higher* value. Entropy is inversely proportional to density, and so the density $u$ can only jump to lower values along a shock. 

We've spent a lot of time now developing the main ideas of weak solutions for conservation laws. We can concretely define what it means to have discontinuous solutions of $\eqref{eq:conslaw}$ and impose conditions to guarantee the uniqueness of this solution. Great. We can finally move on to numerical methods for solving conservation laws. In developing such schemes, our goal will be to ensure that our scheme converges (in some sense) to the entropy solution. 


## **2. Monotone Schemes**
First, a cautionary tale in constructing numerical methods for conservation laws. We consider Burgers equation with the initial condition $u(x,0) = 𝟙\_{\\{x \geqslant 0\\}}$ and the following finite difference scheme:
$$
u\_{j}^{n+1} = u\_{j}^{n} - \frac{\Delta t}{\Delta x} u\_{j}^{n}(u\_{j}^{n} - u\_{j-1}^{n}).
$$

Initially, we have $u\_{j}^{0} = 𝟙\_{\\{x\_{j} \geqslant 0\\}}$. However, this scheme returns $u\_{j}^{n} = u\_{j}^{0}$ for all $n$ and $j$, so the scheme converges to $u(x,t) = u(x,0)$. But in the preceding sections, we established that this can't even be a weak solution (indeed, the Rankine-Hugoniot condition is not satisfied). 

Clearly, we must exercise caution in constructing schemes to solve $\eqref{eq:conslaw}$. We have just seen that a perfectly reasonable-looking finite difference scheme (at least, at first glance) fails to converge to a weak solution, let alone the correct entropy solution. In other cases, it may be possible that a scheme converges to a weak solution but not an entropy solution. We will circle back to this idea shortly.

## **2.1 The Finite Volume Method**
We begin by introducing the **finite volume method**, which discretizes the spatial domain into cells (intervals in 1D, rectangles or other polygons in 2D), and computes an approximation to the average of the solution in each cell. The precise formulation is as follows: consider a bounded, open, connected domain $\Omega = (a,b) \subset \mathbb{R}$ and a finite stopping time $T$, and consider the conservation law given by
$$
\begin{equation}
  \begin{cases}
    u\_{t} + f(u)\_{x} = 0 &\text{in}\;\Omega \times (0,T)\\\\\\
    u(x,0) = u\_{0}(x) &\text{on}\;\Omega \times \\{t=0\\}.
  \end{cases}
  \label{eq:conslaw2}
\end{equation}
$$

We discretize $\Omega$ into equally-sized cells $I\_{j} = (x\_{j-1/ 2}, x\_{j+1/ 2})$:
$$
a = x\_{1/ 2} < x\_{3/ 2} < \dots < x\_{N+1/ 2} = b, \;\;\;\Delta x := x\_{j+1/ 2} - x\_{j-1/ 2}.
$$

Next, we integrate $\eqref{eq:conslaw2}$ over each cell $I\_{j}$:
$$
\begin{align}
  \int_{x\_{j-1/ 2}}^{x\_{j+1/ 2}} (u\_{t} + f(u)\_{x})\;dx &= 0\notag\\\\\\
  \Delta x \frac{d\bar{u}\_{j}}{dt} + f(u\_{j+1/ 2}^{-}) - f(u\_{j-1/ 2}^{+}) &= 0\notag\\\\\\
  \frac{d\bar{u}\_{j}}{dt} + \frac{f(u\_{j+1/ 2}) - f(u\_{j-1/ 2})}{\Delta x} &= 0
\end{align}
$$

Here, $\bar{u}\_{j}$ denotes the cell average of $u$ in cell $I\_{j}$: $\bar{u}\_{j} := \frac{1}{\Delta x}\int\_{I\_{j}} u\;dx$. If our goal is to solve for the cell averages $\bar{u}\_{j}$, you might notice that the second term poses a problem as it is formulated in terms of pointwise values of $u$ at $x\_{j-1/ 2}$ and $x\_{j+1/ 2}$ rather than $\bar{u}\_{j}$. In the first finite volume scheme we consider, we will approximate $u$ in $I\_{j} = (x\_{j-1/ 2}, x\_{j+1/ 2})$ by the cell average $\bar{u}\_{j}$ in $I\_{j}$. 

This approximation is not so well defined as the points $x\_{j+1/ 2}$, which lie at the intersection between the two cells $I\_{j}$ and $I\_{j+1}$. The finite volume scheme remedies this by introducing a *numerical flux function* which takes as inputs cell averages and returns an approximation of the flux function at the cell endpoints $x\_{j+1/ 2}$, i.e. $f(u\_{j+1/ 2}) \approx \hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1})$ and $f(u\_{j-1/ 2}) \approx \hat{f}(\bar{u}\_{j-1}, \bar{u}\_{j})$:
$$
\begin{equation}
  \frac{d\bar{u}\_{j}}{dt} + \frac{1}{\Delta x}\left(\hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1}) - \hat{f}(\bar{u}\_{j-1}, \bar{u}\_{j})\right) = 0.
  \label{eq:finitevol}
\end{equation}
$$

Thus, $\hat{f}$ can be viewed in this context as a function which takes in the dual values at the cell interfaces and returns a single, physically relevant value. This will be important when we try to develop higher order finite volume schemes later. 

For now, the only thing left to do to turn $\eqref{eq:finitevol}$ into a usable scheme is to discretize the time derivative. The simplest option is to use a first-order finite difference. In doing so, we arrive at the explicit first-order finite volume scheme.

**Definition 3.** *The explicit first-order finite volume scheme in space and time is given by
$$
\begin{equation}
  \bar{u}\_{j}^{n+1} = \bar{u}\_{j}^{n} - \frac{\Delta t}{\Delta x}\left( \hat{f}(\bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n}) - \hat{f}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}) \right),
  \label{eq:finitevol2}
\end{equation}
$$
where $\bar{u}\_{j}^{n}$ denotes the cell averages at time $t^{n}$ and $\bar{u}\_{j}^{n+1}$ the cell averages at time $t^{n+1} = t^{n} + \Delta t$.*

Notice that we haven't specified what, exactly, the flux function $\hat{f}$ is yet. The choice of $\hat{f}$ is immensely important as it will determine many properties of our scheme, such as convergence to the entropy solution. Before stating some possible choices of $\hat{f}$, we will first cover the properties $\hat{f}$ must satisfy in order to converge to the entropy solution.

**Definition 4.** *A* **monotone scheme** *is one that can be written in the form
$$
\begin{equation}
  \bar{u}\_{j}^{n+1} = \mathcal{H}(\bar{u}\_{j-p}^{n}, \dots, \bar{u}\_{j+q}^{n}),
\end{equation}
$$
where $\mathcal{H}$ is a nondecreasing function in all arguments. In particular, the monotone scheme on the three-point stencil $\\{\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n}\\}$ is given by $\bar{u}\_{j}^{n+1} = \mathcal{H}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n})$.*

This is all well and good, but why should we care about monotone schemes? To begin with, they have a number of remarkable properties that guarantee our numerical solution behaves "nicely." 

**Theorem 4.** *Monotone schemes satisfy the following properties:*

1. **Local maximum principle.**
$$
\min\_{j-p \leqslant i \leqslant j+q} \bar{u}\_{i}^{n} \leqslant \mathcal{H}(\bar{u})\_{j} \leqslant \max\_{j-p \leqslant i \leqslant j+q} \bar{u}\_{i}^{n} \;\;\;\forall j.
$$

2. **Total variation diminishing property.**
$$
||\mathcal{H}(\bar{u})||\_{BV} \leqslant ||\bar{u}||\_{BV},
$$
*where the bounded variation seminorm is given by $||u||\_{BV} = \sum\_{j} |u\_{j} - u\_{j-1}|$.*

3. **Entropy solution.** *Monotone schemes converge to the entropy solution with a rate in $\ell^{1}$ of half-order. This bound is sharp.*

Roughly speaking, the local maximum principle tells us that our cell averages at $t^{n+1}$ cannot be larger or smaller than our averages at $t^{n}$ around the stencil $\\{\bar{u}\_{j-p}^{n}, \dots, \bar{u}\_{j+q}^{n}\\}$. The total variation diminishing property tells us that our cell averages can only decrease in bounded variation. In particular, this means that our cell averages at $t^{n+1}$ cannot be more oscillatory than our cell averages at $t^{n}$. We will see soon some schemes which violate this property and tend to produce spurious oscillations in the numerical solution. The last property is of course self-explanatory: by using a monotone scheme, we are guaranteed to converge to the entropy solution!

We would like for the finite volume scheme in $\eqref{eq:finitevol2}$ to be a monotone scheme, and our goal now is to choose a function $\hat{f}$ in order to achieve this. We have
$$
\mathcal{H}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n}) = \bar{u}\_{j}^{n} - \frac{\Delta t}{\Delta x}\left( \hat{f}(\bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n}) - \hat{f}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}) \right).
$$

$\mathcal{H}$ must be nondecreasing in each argument, so we compute its partial derivatives:
$$
\begin{align}
  \mathcal{H}\_{1} &= \frac{\Delta t}{\Delta x} \hat{f}\_{1}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n})\label{eq:monotone1}\\\\\\
  \mathcal{H}\_{2} &= 1 - \frac{\Delta t}{\Delta x} \left( \hat{f}\_{1}(\bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n}) - \hat{f}\_{2}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}) \right)\label{eq:monotone2}\\\\\\
  \mathcal{H}\_{3} &= -\frac{\Delta t}{\Delta x} \hat{f}\_{2}(\bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n})\label{eq:monotone3}.
\end{align}
$$

Here, $\hat{f}\_{1}$ denotes the partial derivative of $\hat{f}$ with respect to the first argument, and so on. If $\hat{f}$ is increasing in the first argument and decreasing in the second argument, then $\eqref{eq:monotone1}$ and $\eqref{eq:monotone3}$ will be nonnegative. Furthermore, if
$$
\frac{\Delta t}{\Delta x} \left( \hat{f}\_{1}(\bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n}) - \hat{f}\_{2}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n}) \right) \leqslant 1,
$$

then $\eqref{eq:monotone2}$ will also be nonnegative. This last condition is merely a restriction on our time step. By choosing $\Delta t$ small enough, we can always satisfy this condition. The first two conditions tell us that we merely need to choose our numerical flux function $\hat{f}$ such that it is increasing in its first argument and decreasing in its second argument. 

**Proposition 1.** *If $\hat{f}$ satifies the following conditions:*

1. *$\hat{f}$ is Lipschitz continuous in all arguments.*

2. *$\hat{f}(u,u) = f(u)$. This is known as the consistency condition.*

3. *$\hat{f}$ is increasing in the first argument and decreasing in the second argument,*

*then the scheme $\eqref{eq:finitevol2}$ is monotone.*

We will not discuss the first two conditions here, but suffice it to say that these properties are very easy to check. 

### **2.1.1 Numerical Flux Functions**
We will briefly state some common choices of $\hat{f}$ that satisfy the conditions of Proposition 1 and lead to monotone schemes. 

1. **Lax-Friedrichs flux.**
$$
\begin{equation}
  \hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1}) = \frac{1}{2}\left( f(u\_{j}) + f(u\_{j+1}) - \alpha(u\_{j+1} - u\_{j}) \right),
\end{equation}
$$
where $\alpha = \max\_{u} |f'(u)|$. We can take $\alpha = \max\_{[\bar{u}\_{j}, \bar{u}\_{j+1}]} |f'(u)|$. 

2. **Godunov flux.**
$$
\begin{equation}
  \hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1}) = \begin{cases}
    \min\_{\bar{u}\_{j} \leqslant u \leqslant \bar{u}\_{j+1}} f(u), &\bar{u}\_{j} \leqslant \bar{u}\_{j+1}\\\\\\
    \max\_{\bar{u}\_{j} \geqslant u \geqslant \bar{u}\_{j+1}} f(u), &\bar{u}\_{j} > \bar{u}\_{j+1}.
  \end{cases}
\end{equation}
$$

3. **Engquist-Osher flux.**
$$
\begin{align}
  \hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1}) = &\int\_{0}^{\bar{u}\_{j}} \max\\{f'(u),0\\}\;du + \int\_{0}^{\bar{u}\_{j+1}} \min\\{f'(u),0\\}\;du.
\end{align}
$$

There are two more popular choices of flux that do not lead to monotone schemes but are nevertheless still popular.

4. **Roe flux.**
$$
\begin{equation}
  \hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1}) = \begin{cases}
    f(\bar{u}\_{j}), &\frac{f(\bar{u}\_{j+1}) - f(\bar{u}\_{j})}{\bar{u}\_{j+1} - \bar{u}\_{j}} \geqslant 0\\\\\\
    f(\bar{u}\_{j+1}), &\frac{f(\bar{u}\_{j+1}) - f(\bar{u}\_{j})}{\bar{u}\_{j+1} - \bar{u}\_{j}} < 0.
  \end{cases}
\end{equation}
$$

5. **Lax-Wendroff flux.**
$$
\begin{align}
  \hat{f}(\bar{u}\_{j}, \bar{u}\_{j+1}) = &\frac{1}{2}(f(\bar{u}\_{j}) + f(\bar{u}\_{j+1})) -\notag\\\\\\
  &\frac{\Delta t}{2\Delta x}f'\left( \frac{\bar{u}\_{j} + \bar{u}\_{j+1}}{2}\right) \left( f(\bar{u}\_{j+1}) - f(\bar{u}\_{j}) \right)
\end{align}
$$


## **2.2 Implementation**
Without further ado, we will implement the first-order finite volume method using all five of the aforementioned numerical fluxes for Burgers equation with $u(x,0) = \sin{x}$. All of the source code is publicly available on [GitHub](https://github.com/jdongg/numCL) in both a MATLAB and Python implementation. For convenience, we will explain the Python implementation here. First, we initialize the parameters of the scheme.
```python
# specify domain
a = 0
b = 2*np.pi

# initial condition: u(x,0) = alpha + beta*sin(x)
alpha = 0.0
beta  = 1.0

# number of grid points in spatial discretization
N  = 80

# setup grid points
x = np.linspace(a,b,N)     
dx = (b-a)/(N-1);  

# stopping time
T = 1.5
``` 

Next, we compute the cell averages of the initial condition, $\bar{u}\_{j}^{0}$. This is accomplished by integrating the initial condiiton over each cell $I\_{j}$:
```python
# setup array to store cell averages; due to periodicity, we omit the last cell
u = np.zeros((len(x)-1,1)); 

# returns the initial condition of the PDE
def initial_condition(z, alpha, beta):
    return alpha + beta*np.sin(z)

# compute cell averages at t=0
for i in range(0,N-1):
    u[i] = (1.0/dx)*integrate.quad(initial_condition, x[i], x[i+1], args=(alpha,beta))[0]
```

Next, we must choose an appropriate time step for the scheme. Since we are using an explicit time-stepping scheme, we must choose an appropriately small time step (recall that explicit time-stepping is conditionally stable). For the purposes of example, we choose 
$$
\frac{\Delta t}{\Delta x} \max\_{u} |f'(u)| = \frac{1}{2}.
$$

The last step is to perform the time-stepping, whereby for each time step we must compute the numerical fluxes $\hat{f}(\bar{u}\_{j}^{n}, \bar{u}\_{j+1}^{n})$ and $\hat{f}(\bar{u}\_{j-1}^{n}, \bar{u}\_{j}^{n})$:
```python
t = 0.0
while t < T:
    # alpha for the Lax-Friedrichs flux
    A  = np.amax(np.amax(u));

    # compute numerical fluxes fhat_{j+1/2}
    um = u
    up = np.roll(u,-1)
    fR = lf_flux(um, up, A)

    # compute numerical fluxes fhat_{j-1/2} (assuming periodic BCs)
    fL = np.roll(fR,1)

    # first order explicit time-stepping
    u -= dt/dx*(fR - fL)

    # increment time step
    t = t+dt
```

Each of the numerical fluxes in Section 2.1.1 have been implemented and can be viewed on [GitHub](https://github.com/jdongg/numCL). Figure 2 shows the finite volume solutions at $T=1.5$. 

{{< figure src="gallery/finiteVolumeBurgers.png" caption="**Figure 2:** Solution of Burgers equation with sine initial data using the first-order finite volume method.">}}

The first observation we make is that there is little difference between the Roe, Godunov, and Engquist-Osher fluxes, at least for this example problem. The Godunov and Engquist-Osher fluxes lead to monotone schemes which converge to the entropy solution. The Roe flux does not lead to a monotone scheme, and in very rare cases may provide solutions with the incorrect shock speed. Nevertheless, in most cases the Roe scheme works well and is a popular choice for the numerical flux. 

Next, the Lax-Friedrichs flux is noticeably more dissipative than any of the other results. It is known to "smear" discontinuities and add numerical diffusion. However, the Lax-Friedrichs flux is a very popular choice for the numerical flux since it leads to a monotone scheme and is straightforward to implement. 

Lastly, the Lax-Wendroff flux is...problematic, to say the least. Notice that the numerical solution develops spurious (nonphysical) oscillations in the vicinity of the shock at $x=\pi$. If the Lax-Wendroff scheme converges, it converges to a weak solution, but the scheme is not monotone and does not exhibit the total variation-diminishing (TVD) property. Lax-Wendroff is actually a *second-order method* in disguise, but we'll discuss this much later and how the order of a numerical method is related to the TVD property for conservation laws. 


## **3. Key Takeaways**
If I haven't bored you to death by this point, you should have a basic understanding of weak solutions for conservation laws and the mathematical foundations for discontinuous solutions. The study of such solutions is important because they occur frequently in fluid dynamics and other related fields. You should also be able to implement a basic first-order finite volume scheme for a variety of one-dimensional conservation laws. 

We've only covered first-order schemes in this post, and you might have already guessed that these methods converge too slowly to be of much use for practical problems. Next time, we'll cover higher order finite volume methods for conservation laws and touch on some of the theory regarding higher-order methods (remember how I said the Lax-Wendroff flux actually yields a second-order method? We'll circle back to that soon). Unfortunately, there's no free lunch and higher order methods come with their own set of issues. 

Below are two texts by Randy LeVeque which are helpful in fleshing out much of the details of the Theorems given in this post. I've found them to be immensely helpful in my Ph.D. studies. 


### **References**
[1] LeVeque, Randall J. [Finite volume methods for hyperbolic problems](http://staff.washington.edu/rjl/book2/sample.pdf). Vol. 31. Cambridge university press, 2002. 

[2] LeVeque, Randall J. [Numerical methods for conservation laws](https://pdfs.semanticscholar.org/1470/c6f43c769572c4cfc94ffc9c5710484ff1e5.pdf). Vol. 132. Basel: Birkhäuser, 1992.


### **Codes**
[NumCL repository on GitHub](https://github.com/jdongg/numCL)