# 重叠型Schwarz算法利用Crank-Nicolson格式解Fisher-kpp方程

Fisher-kpp方程有如下的形式和解析解：

![Fisher-kpp Equation](Fisher%E6%96%B9%E7%A8%8B%E8%A7%A3%E6%9E%90%E8%A7%A3.png)

这份代码解的是如下的Fisher-kpp方程：
$$\frac{\partial u}{\partial t} = \alpha\frac{\partial^2 u}{\partial x^2} + \beta u(1-u)$$

这份代码使用了如下的Crank-Nicolson格式：
<div align="center"><img style="background: white;" src="svg\pL9z30FgOA.svg"></div>

其中
$$f(u) = \beta u(1-u)$$
对上式差分得
<div align="center"><img style="background: white;" src="svg\eq3oWOtCor.svg"></div>

易知，在每一时间步长需要求解如下方程组：
<div align="center"><img style="background: white;" src="svg\ysBH8x27mD.svg"></div>

其中
$$r = \alpha\frac{\tau}{h^2}$$
由于 $U^n_0,U^n_{J}$ 已知且$U^n_0,U^n_{J} \neq 0$
简化上式得
<div align="center"><img style="background: white;" src="svg\ZUxt6uCYqb.svg"></div>

其中
<div align="center"><img style="background: white;" src="svg\uApcqarQKL.svg"></div>

为简化表达，由于$\overrightarrow{U}^{n}$已知，不妨令
$$D=B\overrightarrow{U}^{n}+\frac{\tau}{2}F(\overrightarrow{U}^{n})+C$$
则$D$已知，只需求解以下方程组：
$$A\overrightarrow{U}^{n+1}-\frac{\tau}{2}F(\overrightarrow{U}^{n+1})=D$$
考虑Newton迭代法：
$$f(x)=0 \\ x_{n+1}=x_n-\frac{f(x)}{f'(x)}$$
而矩阵形式的Newton迭代法为：
$$\overrightarrow{x}_{n+1}=\overrightarrow{x}_n - J^{-1}_n F(\overrightarrow{x}_n)$$
其中
<div align="center"><img style="background: white;" src="svg\WoPDhVt2u9.svg"></div>

这里，使
$$G(\overrightarrow{U})=A\overrightarrow{U}^{n+1}-\frac{\tau}{2}F(\overrightarrow{U}^{n+1})-D=0$$
取 $\overrightarrow{U}_0=\overrightarrow{U}^n$ ，则Newton迭代法过程如下：

$$\overrightarrow{U}_{n+1}=\overrightarrow{U}_n-J^{-1}_nG(\overrightarrow{U}_n)$$

其中
<div align="center"><img style="background: white;" src="svg\44eI3Khvmu.svg"></div>

迭代2次就可取得不错的结果，最后令 $\overrightarrow{U}^{n+1}=\overrightarrow{U}_2$ ，进行下一个时间步长的运算。
