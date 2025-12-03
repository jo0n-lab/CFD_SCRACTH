# PLIC VOF

## Governing Equation

### Volume Fraction

$$
\begin{equation*}
\begin{align*}
\frac{\partial C}{\partial t} + \nabla \cdot (\mathbf{u} C) 
&= \frac{\partial C}{\partial t} + \frac{\partial C}{\partial x} \cdot u_x + \frac{\partial C}{\partial y} \cdot u_y + C\left(\frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y}\right)\\
&=
\end{align*}
\end{equation*}
$$

$$
\begin{equation*}
\begin{align*}
\frac{C^* - C^n}{\Delta t}\Delta V + F_{\text{out},x} - F_{\text{in},x}  - \int_{\Omega} C (\nabla \cdot \mathbf{u}) \, dV
\end{align*}
\end{equation*}
$$

