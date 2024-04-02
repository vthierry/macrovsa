with(Statistics):
with(LinearAlgebra):

# Draw a unary random vector of dimension d

draw := d -> convert(Sample(RandomVariable(Normal(0, 1 / evalf(sqrt(d)))), d), Vector):

# Verifies that vectors are almost unary and orthogonal

ok := map(proc(d)
  local U := map(proc() Norm(draw(d), Euclidean) end, [$1..1000]),
        Z := map(proc() DotProduct(draw(d), draw(d)) end, [$1..1000]):
  [d,
   Mean(U) - 1,                        # Magnitude bias
   StandardDeviation(U) - 1 / sqrt(d), # Magnitude standard-deviation bias
   Mean(Z),                            # Orthogonality bias
   StandardDeviation(Z) - 1 / sqrt(d)] # Orthogonality standard-deviation bias
  end, [seq(100..1000, 100)]):

#plots[listplot](map(t->[op(1,t), op(2,t)], ok), title = "Magnitude bias");
#plots[listplot](map(t->[op(1,t), op(3,t)], ok), title = "Magnitude standard-deviation bias");
#plots[listplot](map(t->[op(1,t), op(4,t)], ok), title = "Orthogonality bias");
#plots[listplot](map(t->[op(1,t), op(5,t)], ok), title = "Orthogonality standard-deviation bias");

print("Magnitude bias"): DataSummary(map(t -> evalf(op(2,t)), ok));
print("Magnitude standard-deviation bias"): DataSummary(map(t -> evalf(op(3,t)), ok));
printf("Orthogonality bias"): DataSummary(map(t -> evalf(op(4,t)), ok));
printf("Orthogonality standard-deviation "): DataSummary(map(t -> evalf(op(5,t)), ok));

# Calculates the z-score between two vectors, with a common component of angle theta
zscore := proc(d, theta)
  local u, v: u := draw(d): v := draw(d):
  DotProduct(u, cos(theta) * v + u * sin(theta)) * evalf(sqrt(d));
end:

# Draws some z-score distribution of independent vectors
zscores := (d, k) -> map('zscore(d, 0)',[$1..k]);
zscores_2 := zscores(100,1000):
zscores_3 := zscores(1000,1000):
zscores_4 := zscores(10000,1000):

H0 := DensityPlot(Normal(0,1), color = black):
H2 := Histogram(zscores_2, color = gray):
H3 := Histogram(zscores_3, color = gray):

DataSummary(zscores_2);
DataSummary(zscores_3);
DataSummary(zscores_4);

plotsetup(jpeg, plotoutput="z_score_2.jpg", plotoptions="width=600,height=600");
plots[display](H2, H0, title="Z-score distribution, d = 100");
plotsetup(jpeg, plotoutput="z_score_3.jpg", plotoptions="width=600,height=600");
plots[display](H3, H0, title="Z-score distribution, d = 1000");
plotsetup(x11);

# Calculates z-score as a function of theta

T0 := plot(cos(Pi/2*a/90.0),a=1..90,color=black,thickness=4):

plotsetup(jpeg, plotoutput="z_score_2a.jpg", plotoptions="width=600,height=600");
zscores_theta := map(k-> map(a->zscore(100, Pi/180.0*(90-a))/sqrt(100),[$0..90]),[$1..10]):
T2 := ErrorPlot(zscores_theta, color = map(k->black,[$1..10])):
plots[display](T0, T2, title="Z-score as a function of angular dependency, d = 100");
plotsetup(jpeg, plotoutput="z_score_3a.jpg", plotoptions="width=600,height=600");
zscores_theta := map(k-> map(a->zscore(1000, Pi/180.0*(90-a))/sqrt(1000),[$0..90]),[$1..10]):
T3 := ErrorPlot(zscores_theta, color = map(k->black,[$1..10])):
plots[display](T0, T3, title="Z-score as a function of angular dependency, d = 1000");
plotsetup(x11);

#plot(zscore(100, Pi/180*(90-a)), a = 0..90);



