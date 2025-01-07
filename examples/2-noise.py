from minke.noise import AdvancedLIGO

noise = AdvancedLIGO()

data = noise.time_series(duration=4, sample_rate=16384)

f = data.plot()
f.savefig("noise.png")
