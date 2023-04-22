# 重叠型Schwarz算法利用Crank-Nicolson格式解Fisher-kpp方程

Fisher-kpp方程有如下的形式和解析解：

![Fisher-kpp Equation](Fisher%E6%96%B9%E7%A8%8B%E8%A7%A3%E6%9E%90%E8%A7%A3.png)

这份代码解的是如下的Fisher-kpp方程：
$$ \frac{\partial u}{\partial t} = \alpha\frac{\partial^2 u}{\partial x^2} + \beta u(1-u) $$

这份代码使用了如下的Crank-Nicolson格式：
$$ 
\frac{\partial u}{\partial t}(x_j,t_{n+1/2})=\alpha\frac{\partial^2 u}{\partial x^2}(x_j,t_{n+1/2}) + f(u(x_j,t_{n+1/2}))
$$
其中
$$ f(u) = \beta u(1-u) $$
对上式差分得
$$
\begin{align*} 
\frac{U^{n+1}_j-U^n_j}{\tau} = \alpha(\frac{U^{n+1}_{j+1}-2U^{n+1}_j+U^{n+1}_{j-1}}{2h^2}+\frac{U^{n}_{j+1}-2U^{n}_j+U^{n}_{j-1}}{2h^2})  \\ 
+ \frac{1}{2}[f(U^{n+1_j})+f(U^n_j)] , j=1,2,\dots,J-1 
\end{align*}
$$
易知，在每一时间步长需要求解如下方程组：
$$ 
\begin{bmatrix} 
-\frac{r}{2} & 1+r & -\frac{r}{2} \\
             & -\frac{r}{2} & 1+r & -\frac{r}{2} \\
             &              & \ddots & \ddots & \ddots \\
             &              &        & -\frac{r}{2} & 1+r & -\frac{r}{2}
\end{bmatrix}
\begin{bmatrix} 
U^{n+1}_0 \\
U^{n+1}_1 \\
\vdots    \\
U^{n+1}_{J}
\end{bmatrix}
-\frac{\tau}{2}
\begin{bmatrix} 
f(U^{n+1}_1) \\
f(U^{n+1}_2) \\
\vdots    \\
f(U^{n+1}_{J-1})
\end{bmatrix}
=
\begin{bmatrix} 
\frac{r}{2} & 1-r & \frac{r}{2} \\
            & \frac{r}{2} & 1-r & \frac{r}{2} \\
            &              & \ddots & \ddots & \ddots \\
            &              &        & \frac{r}{2} & 1-r & \frac{r}{2}
\end{bmatrix}
\begin{bmatrix} 
U^{n}_0 \\
U^{n}_1 \\
\vdots    \\
U^{n}_{J}
\end{bmatrix}
+\frac{\tau}{2}
\begin{bmatrix} 
f(U^{n}_1) \\
f(U^{n}_2) \\
\vdots    \\
f(U^{n}_{J-1})
\end{bmatrix}
$$
其中
$$ r = \alpha\frac{\tau}{h^2} $$
由于
$U^n_0,U^n_{J}$已知且$U^n_0,U^n_{J} \neq 0$
简化上式得
$$
A\overrightarrow{U}^{n+1}-\frac{\tau}{2}F(\overrightarrow{U}^{n+1})=B\overrightarrow{U}^{n}+\frac{\tau}{2}F(\overrightarrow{U}^{n})+C
$$
其中
$$
A=
\begin{bmatrix} 
 1+r & -\frac{r}{2} \\
-\frac{r}{2} & 1+r & -\frac{r}{2} \\
             & \ddots & \ddots & \ddots \\
             &        & -\frac{r}{2} & 1+r
\end{bmatrix} \\
B=
\begin{bmatrix} 
1-r & \frac{r}{2} \\
\frac{r}{2} & 1-r & \frac{r}{2} \\
             & \ddots & \ddots & \ddots \\
             &        & \frac{r}{2} & 1-r
\end{bmatrix} \\
F(\overrightarrow{U^{n}})=
\begin{bmatrix} 
f(U^{n}_1) \\
f(U^{n}_2) \\
\vdots    \\
f(U^{n}_{J-1})
\end{bmatrix} \\
C=
\frac{r}{2}
\begin{bmatrix} 
U^n_0+U^{n+1}_0 \\
0 \\
\vdots    \\
0 \\
U^n_J+U^{n+1}_J
\end{bmatrix}
$$
为简化表达，由于$\overrightarrow{U}^{n}$已知，不妨令
$$ D=B\overrightarrow{U}^{n}+\frac{\tau}{2}F(\overrightarrow{U}^{n})+C $$
则$D$已知，只需求解以下方程组：
$$ A\overrightarrow{U}^{n+1}-\frac{\tau}{2}F(\overrightarrow{U}^{n+1})=D $$
考虑Newton迭代法：
$$
f(x)=0 \\
x_{n+1}=x_n-\frac{f(x)}{f'(x)}
$$
而矩阵形式的Newton迭代法为：
$$
\overrightarrow{x}_{n+1}=\overrightarrow{x}_n - J^{-1}_n F(\overrightarrow{x}_n) \\
$$
其中
$$
J_n=
\begin{bmatrix}
\frac{\partial f_1(x)}{x_1} & \frac{\partial f_1(x)}{x_2} & \dots & \frac{\partial f_1(x)}{x_{J-1}} \\
\frac{\partial f_2(x)}{x_1} & \frac{\partial f_2(x)}{x_2} & \dots & \frac{\partial f_2(x)}{x_{J-1}} \\
\vdots                      & \vdots                      &       & \vdots                          \\
\frac{\partial f_{J-1}(x)}{x_1} & \frac{\partial f_{J-1}(x)}{x_2} & \dots & \frac{\partial f_{J-1}(x)}{x_{J-1}}
\end{bmatrix}
$$
这里，使
$$
G(\overrightarrow{U})=A\overrightarrow{U}^{n+1}-\frac{\tau}{2}F(\overrightarrow{U}^{n+1})-D=0
$$
取$\overrightarrow{U}_0=\overrightarrow{U}^n$，则Newton迭代法过程如下：
$$
\overrightarrow{U}_{n+1}=\overrightarrow{U}_n-J^{-1}_nG(\overrightarrow{U}_n) \\
$$
其中
$$
J_n=A-\frac{\tau}{2}
\begin{bmatrix}
\beta(1-2U^{n}_1) &         &                   \\
                  & \ddots  &                   \\
                  &         & \beta(1-2U^{n}_{J-1}) \\
\end{bmatrix}
$$
迭代2次就可取得不错的结果，最后令$\overrightarrow{U}^{n+1}=\overrightarrow{U}_2$，进行下一个时间步长的运算。