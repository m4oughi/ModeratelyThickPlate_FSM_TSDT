# Buckling analysis of moderately Thick Plate using Finite Strip Method based on Third-order Shear Deformation Theory

This repository is dedicated to a little part of my master’s thesis. The focus of this code is on understanding how the mechanical buckling of a moderately thick plate behaves. The analysis is conducted using the Finite Strip Method (FSM) based on Third-order shear deformation theory. I used a Python package called Sympy to develop the formulation and get accurate matrix computations. You can verify this code with the research presented by Foroughi and Azhari [1]. For example, if you try (E=2000000, noo=0.3, G=770000, b=1, a=1, h=0.1, m=1, n=10, ax=1, ay=0,axy=0), you will get 36.8289 that you can find it in [1] of Tabel 2 as 36.8315. <br><br>

Feel free to explore the repository, and if you have any questions or want to chat about this, I’m all ears! <br><br>

Reference:
[1]. Foroughi, H., & Azhari, M. (2014). Mechanical buckling and free vibration of thick functionally graded plates resting on elastic foundation using the higher order B-spline finite strip method. Meccanica, 49, 981-993.
