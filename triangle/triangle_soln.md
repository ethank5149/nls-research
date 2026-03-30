## Problem Statement
Given:
$$\begin{aligned}
\sin^2\theta+\mu\cos^2\theta+\sqrt{3}(\mu-1)\sin\theta\cos\theta&\geq0\quad(\text{I})\\
\sin^2\theta+\mu\cos^2\theta-\sqrt{3}(\mu-1)\sin\theta\cos\theta&\geq0\quad(\text{II})\\
(\mu+1)-2(\mu-1)\cos(2\theta)&\geq0\quad(\text{III})\\
\end{aligned}$$
With
$$
\theta\in\left[0,\;\frac{\pi}{2}\right],\quad\mu\in\left(0,\;\infty\right)$$
## Solution
Let
$$C=\cos\theta,\;S=\sin\theta$$
Since
$$C^2=\tfrac{1+\cos2\theta}{2},\;S^2=\tfrac{1-\cos2\theta}{2},\;CS=\tfrac{\sin2\theta}{2}$$
(I) & (II) are equivalent to
$$(\mu+1)+(\mu-1)\cos(2\theta)\pm\sqrt3(\mu-1)\sin(2\theta)\ge0$$
$$\Updownarrow$$
$$(\mu+1)+2(\mu-1)\cos\bigl(2\theta\mp\tfrac{\pi}{3}\bigr)\ge0\quad(\text{IV}),(\text{V})$$
(III), (IV), and (V) require:
$$(\mu+1)+2(\mu-1)\cos\alpha \ge 0.\qquad\alpha\in{2\theta-\tfrac{\pi}{3},\,2\theta,\,2\theta+\tfrac{\pi}{3}}$$
Let
$$\nu(\mu)=\frac{\mu+1}{2|\mu-1|},\quad\varphi(\mu)=\arccos(\nu(\mu))$$
s.t. 
$$\varphi(\mu)=
\begin{cases}
\frac{\pi}{3},&\nu\leq\frac{1}{2}\\
0<\varphi<\frac{\pi}{3},&\frac{1}{2}<\nu<1\\
0,&1\leq\nu\\
\end{cases}$$

## Results
### $0<\mu<1/3$
Since
$$\mu-1<0\Rightarrow\nu=\frac{\mu+1}{2(1-\mu)}\in\left(\frac{1}{2},\;1\right).$$
We have
$$\varphi=\arccos\left(\frac{\mu+1}{2(1-\mu)}\right)\in\left(0,\;\frac{\pi}{3}\right)$$
s.t.
$$\theta\in\left[0,\;\frac{\pi/3-\varphi}{2}\right]\;\cup\;
\left[\frac{\pi/3+\varphi}{2},\;\frac{\pi-\varphi}{2}\right].$$
### $1/3\le\mu\le3$:
All $\theta\in\left[0,\frac{\pi}{2}\right]$
### $\mu>3$:
Since 
$$\mu-1>0\Rightarrow\nu=\frac{\mu+1}{2(\mu-1)}\in\left(\frac{1}{2},\;1\right).$$
We have
$$\varphi=\arccos\left(\frac{\mu+1}{2(\mu-1)}\right)\in\left(0,\;\frac{\pi}{3}\right)$$
s.t.
$$\theta\in\left[\frac{\varphi}{2},\;\frac{2\pi/3-\varphi}{2}\right]\;\cup\; \left[\frac{2\pi/3+\varphi}{2},\;\frac{\pi}{2}\right].$$

# Summary
| $$\mu-\text{range}$$ | $$\text{Allowed}\;\;\theta-\text{range}$$ |
| :-: | :-: |
| $$0<\mu<\dfrac{1}{3}$$ | $$\left[0,\;\frac{\pi}{6}-\frac{\varphi}{2}\right]\;\cup\;\left[\frac{\pi}{6}+\frac{\varphi}{2},\;\frac{\pi}{2}-\frac{\varphi}{2}\right]$$ |
| $$\frac{1}{3}\le\mu\le3$$ | $$\left[0,\;\frac{\pi}{2}\right]$$ |
| $$3<\mu$$ | $$\left[\frac{\varphi}{2},\;\frac{\pi}{3}-\frac{\varphi}{2}\right]\;\cup\;\left[\frac{\pi}{3}+\frac{\varphi}{2},\;\frac{\pi}{2}\right]$$ |