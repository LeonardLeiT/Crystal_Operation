# Curvature calculation for the Schwarz D surface

The theoretical curvature associated with each box geometry was calculated using a trigonometric approximation to the Schwarz D surface. The three box vectors, $\mathbf a$, $\mathbf b$ and $\mathbf c$, form the columns of the matrix $\mathbf B=[\mathbf a\ \mathbf b\ \mathbf c]$. The surface was described in fractional coordinates $\mathbf u=(u,v,w)^{\mathsf T}$ along these vectors. With $X=\pi u$, $Y=\pi v$ and $Z=\pi w$, the theoretical surface was defined by the zero level set of the following field. [Gandy et al., 2001](https://doi.org/10.1016/S0009-2614(00)01418-4)

$$
\begin{aligned}
f(\mathbf u)={}&\sin X\sin Y\sin Z+\sin X\cos Y\cos Z\\
&+\cos X\sin Y\cos Z+\cos X\cos Y\sin Z=0.
\end{aligned}
\tag{1}
$$

The field was differentiated analytically with respect to fractional coordinates. Its derivatives were then transformed into Cartesian coordinates to obtain the gradient $\mathbf g$ and Hessian matrix $\mathbf Q$,

$$
\mathbf g=\mathbf B^{-\mathsf T}\nabla_{\mathbf u}f,
\qquad
\mathbf Q=\mathbf B^{-\mathsf T}
\left(\nabla_{\mathbf u}^{2}f\right)\mathbf B^{-1}.
\tag{2}
$$

This transformation accounts for both box edge lengths and inter-vector angles and applies to orthogonal and non-orthogonal boxes. Taking the unit normal as $\mathbf n=\mathbf g/\|\mathbf g\|$ and using the sign convention $H=\tfrac12\nabla\cdot\mathbf n$, the local mean curvature was calculated as

$$
H=\frac{\kappa_1+\kappa_2}{2}
=\frac{\operatorname{tr}(\mathbf Q)-\mathbf n^{\mathsf T}\mathbf Q\mathbf n}
{2\|\mathbf g\|}.
\tag{3}
$$

Here, $\kappa_1$ and $\kappa_2$ are the principal curvatures. They are the two eigenvalues of $\mathbf T^{\mathsf T}\mathbf Q\mathbf T/\|\mathbf g\|$, where the columns of $\mathbf T$ are two orthonormal tangent vectors. Equation (3) evaluates their mean directly at regular surface points, where $\|\mathbf g\|>0$. [Goldman, 2005](https://doi.org/10.1016/j.cagd.2005.06.005)

The fractional domain $[0,1]^3$ was sampled using 128 uniform intervals along each direction. The zero level set was extracted using the Marching Cubes algorithm, and the resulting triangular mesh was mapped into Cartesian space using the box matrix. Local mean curvature $H_i$ was evaluated at each triangle centroid, and the corresponding triangle area $A_i$ was calculated in Cartesian space. The dimensionless, area-weighted mean absolute curvature was then obtained as follows. [Lewiner et al., 2003](https://doi.org/10.1080/10867651.2003.10487582)

$$
H^*=V^{1/3}\frac{\sum_i A_i|H_i|}{\sum_i A_i},
\qquad V=|\det\mathbf B|.
\tag{4}
$$

Taking the absolute value prevents cancellation between positive and negative local curvatures. Normalization by $V^{1/3}$ removes the dependence on overall size, making $H^*$ invariant under uniform scaling of the box. This metric was used to compare surface curvature across boxes with different aspect ratios and inter-vector angles.
