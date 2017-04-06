using PDMP, Base.Test, Polynomials

# the tests have crude accuracy (because the benchmark is drawn from visual
# inspections with Grapher). This doesn't matter much.

# CASE A # p = (s+1)(s^2+1)
p = Poly([1.0,1.0])*Poly([1.0,0.0,1.0])

@test abs(pmf_caseA(2.0, p) - 0.9782) <= 1e-3

# CASE B # p = (s-1)(s^2+1)
p = Poly([-1.0,1.0])*Poly([1.0,0.0,1.0])

@test abs(pmf_caseB(2.0, p, 1.0) - 2.0165) <= 1e-3

# CASE C # p (s+1)(s-1)(s-2)
p = Poly([1.0,1.0])*Poly([-1.0,1.0])*Poly([-2.0,1.0])

@test abs(pmf_caseC(0.5014, p, 1.0, 2.0) - 0.2761) <= 1e-3
@test abs(pmf_caseC(1.4090, p, 1.0, 2.0) - 2.3968) <= 1e-3

# CASE D # p = (s-1)(s-2)(s-3)
p = Poly([-1.0,1.0])*Poly([-2.0,1.0])*Poly([-3.0,1.0])

@test abs(pmf_caseD(0.1398, p, 1.0, 2.0, 3.0) - 1.4978) <= 1e-3
@test abs(pmf_caseD(0.8405, p, 1.0, 2.0, 3.0) - 3.5928) <= 1e-3
