# CFD_SCRACTH
Study of Computational Fluid Dynamics from SCRATCH with Python &amp; Julia

### TODO LIST

1. Transient Solver
2. K-Epsilon Model
3. Multigrid Solver
4. Testcases & Benchmarks




### **Conservation vs Non-Conservation expressions**

Whoever experienced Computational Fluid Dynamics, one will find that "Continuity Equation" is widely used to evaluate the accuracy of the driven simulation.

Since "Continuity" stands for all physical phenomena, which is "Mass Conservation", it is useful how CFD handles this equation in every form of physical criteria.

It starts with very familiar equation: **"Material Derivative"**

$$
\frac{D\phi}{Dt}\equiv \frac{\partial \phi}{\partial t} + (u \cdot \nabla) \phi \qquad \text{(Material Derivative)}
$$

But this form ***actively*** hides **Continuity**. That is, it assumes mass conservation is always conserved.

However in ***simulation world***, thing don't always happen as it is. CFD conceive to fit simulation into real-world. Then, we cannot assume mass convervation happens naturally but it's our objective to fit!

Therefore another form raised:

$$
\frac{\partial(\rho \phi)}{\partial t} + \nabla \cdot (\rho \phi u) \qquad \text{(Conservative Expression)}
$$

The relationship between Material Derivative & Conservative Expression is that, former one is assuming mass conservation be conserved.

The following is derivation of Material Derivative from Conservative form.

$$
\frac{\partial(\rho \phi)}{\partial t} + \nabla \cdot (\rho \phi u)=

\rho\frac{\partial\phi}{\partial t}+ \phi\frac{\partial\rho}{\partial t} + \left(\phi \nabla \cdot (\rho u) + \rho u \cdot \nabla \phi \right) \\[1em]
\qquad\qquad\qquad\qquad= \underbrace{\phi\left( \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) \right)}_\text{Continuity Equation} + \left( \frac{\partial \phi}{\partial t} + (u\cdot \nabla) \phi \right)
$$

If $\phi$ is given as a vector, for example, $\Phi=\mathbf{u}$
(quite more complicated):

$$
\frac{\partial(\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \mathbf{u})= \frac{\partial(\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u})\\
$$

Since


$$
\qquad\qquad\qquad\qquad= \underbrace{\phi\left( \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) \right)}_\text{Continuity Equation} + \left( \frac{\partial \phi}{\partial t} + (u\cdot \nabla) \phi \right)
$$


### Matrix Multiplication

1. Dyadic Product(=Outer Product): result in tensor
$$ 
\mathbf{A} \otimes \mathbf{B} \equiv 
\begin{bmatrix}
a_1b_1 & a_1b_2 & a_1b_3 \\
a_2b_1 & a_2b_2 & a_2b_3 \\
a_3b_1 & a_3b_2 & a_3b_3 \\
\end{bmatrix} \\[1em]

\mathbf{A} \otimes \mathbf{B} = \mathbf{A}\mathbf{B}^T=(a \otimes b)_{i,j}=a_ib_j
$$  


2. Dot Product: result in scalar
$$
\mathbf{A}\cdot\mathbf{B} \equiv a_1b_1 + a_2b_2 + a_3b_3\\
\mathbf{A} \cdot \mathbf{B} = \sum a_ib_i
$$

3. Cross Product: result in vector
$$
\mathbf{A} \times \mathbf{B} \equiv 
\left<a_2b_3-a_3b_2,a_3b_1-a_1b_3,a_1b_2-a_2b_1\right>\\[1em]
\mathbf{A} \times \mathbf{B}=
\begin{vmatrix}
e_1 & e_2 & e_3 \\
a_1 & a_2 & a_3 \\
b_1 & b_2 & b_3 \\
\end{vmatrix} \\[1em]
$$

#### Vector vs Tensor
||order|component|meaning|
|--|--|--|--|
|Scalar|0|single|magnitude|
|Vector|1|multiple|magnitude + direction|
|2nd Tensor|2|square|linearity btw vectors|
|3rd Tensor|3|cubic|high-order btw vectors|



 

