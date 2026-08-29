import math

RHO = 1.055
VISCOSITY = 0.045
A0 = 0.126
BETA = 606060.0
ALPHA = 4.0 / 3.0
LENGTH = 4.0

KINEMATIC_VISCOSITY = VISCOSITY / RHO
FRICTION_KR = 8.0 * math.pi * KINEMATIC_VISCOSITY

C0 = math.sqrt(BETA * math.sqrt(A0) / (2.0 * RHO * A0))

PULSE_AMPLITUDE = 5.0
PULSE_HALF_WIDTH = 4.0e-4
PULSE_CENTER = 1.0e-3
PULSE_WINDOW = 2.0 * PULSE_CENTER
PULSE_CSV_DT = 2.0e-7

def pulse_value(t):
    if abs(t-PULSE_CENTER) > PULSE_HALF_WIDTH:
        return 0.0
    return PULSE_AMPLITUDE * math.sin( \
        math.pi * (t - PULSE_CENTER + PULSE_HALF_WIDTH) / (2 * PULSE_HALF_WIDTH))**2

def write_pulse_csv(path):
    n = round(PULSE_WINDOW / PULSE_CSV_DT)
    with open(path, "w") as f:
        f.write("time, flow_rate\n")
        for i in range(n+1):
            t = i*PULSE_CSV_DT
            f.write(f"{t:.9f}, {pulse_value(t):.8f}\n")