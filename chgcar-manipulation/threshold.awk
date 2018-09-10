# attempt to compute the threshold value for the charge density FT
($1 == 0){rhomin = 1e-4 * $2; next}
{if ($2>rhomin) {
  print $0;
  G[n] = log($1);
  rho[n] = log($2);

  G_sum += G[n];
  G2_sum += G[n]*G[n];
  rho_sum += rho[n];
  rhoG_sum += rho[n]*G[n];
  
  ++n;
 }
} 

END {
G_sum = G_sum/n;
G2_sum = G2_sum/n;
rho_sum = rho_sum/n;
rhoG_sum = rhoG_sum/n;
b = (rhoG_sum - rho_sum*G_sum)/(G2_sum - G_sum*G_sum);
a = 0;
for (i=0; i<n; ++i) {
  aguess = exp(rho[i] -b*G[i]);
  if (aguess>a) a=aguess;
}
printf("#rho_trunc= %.8le |G|^ %.8lf\n", a, b);
}
