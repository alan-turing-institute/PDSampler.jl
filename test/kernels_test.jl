using PDMP, Base.Test

srand(1240)

d      = 10
normal = randn(d)
v1     = randn(d)
v2     = copy(v1)
I      = eye(d)

@test norm(PDMP.reflect_bps!(normal,v1)-PDMP.reflect_bps!(normal,v2,I)) <= 1e-12
