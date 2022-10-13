#!/usr/bin/gawk -f
#
{
  i=$1;
  n[i]++;
  c[i"-"n[i]]=$3;
#
} END {
#
  for(i in n) {
    for(j=1;j<=n[i];j++) s[i]+=c[i"-"j];
    a[i]=s[i]/n[i];
    for(j=1;j<=n[i];j++) ss[i]+=(c[i"-"j]-a[i])**2
  }
#
#  for(i in n) print i,n[i],s[i]/n[i],sqrt(ss[i]/n[i])
#
   printf "1 1 %6.2f %6.2f\n",s["C11"]/n["C11"], sqrt(ss["C11"]/n["C11"]);
   printf "2 2 %6.2f %6.2f\n",s["C22"]/n["C22"], sqrt(ss["C22"]/n["C22"]);
   printf "3 3 %6.2f %6.2f\n",s["C33"]/n["C33"], sqrt(ss["C33"]/n["C33"]);
   printf "1 2 %6.2f %6.2f\n",s["C12"]/n["C12"], sqrt(ss["C12"]/n["C12"]);
   printf "1 3 %6.2f %6.2f\n",s["C13"]/n["C13"], sqrt(ss["C13"]/n["C13"]);
   printf "2 3 %6.2f %6.2f\n",s["C23"]/n["C23"], sqrt(ss["C23"]/n["C23"]);
   printf "4 4 %6.2f %6.2f\n",s["C44"]/n["C44"], sqrt(ss["C44"]/n["C44"]);
   printf "5 5 %6.2f %6.2f\n",s["C55"]/n["C55"], sqrt(ss["C55"]/n["C55"]);
   printf "6 6 %6.2f %6.2f\n",s["C66"]/n["C66"], sqrt(ss["C66"]/n["C66"]);
   printf "1 4 %6.2f %6.2f\n",s["C14"]/n["C14"], sqrt(ss["C14"]/n["C14"]);
   printf "1 5 %6.2f %6.2f\n",s["C15"]/n["C15"], sqrt(ss["C15"]/n["C15"]);
   printf "1 6 %6.2f %6.2f\n",s["C16"]/n["C16"], sqrt(ss["C16"]/n["C16"]);
   printf "2 4 %6.2f %6.2f\n",s["C24"]/n["C24"], sqrt(ss["C24"]/n["C24"]);
   printf "2 5 %6.2f %6.2f\n",s["C25"]/n["C25"], sqrt(ss["C25"]/n["C25"]);
   printf "2 6 %6.2f %6.2f\n",s["C26"]/n["C26"], sqrt(ss["C26"]/n["C26"]);
   printf "3 4 %6.2f %6.2f\n",s["C34"]/n["C34"], sqrt(ss["C34"]/n["C34"]);
   printf "3 5 %6.2f %6.2f\n",s["C35"]/n["C35"], sqrt(ss["C35"]/n["C35"]);
   printf "3 6 %6.2f %6.2f\n",s["C36"]/n["C36"], sqrt(ss["C36"]/n["C36"]);
   printf "4 5 %6.2f %6.2f\n",s["C45"]/n["C45"], sqrt(ss["C45"]/n["C45"]);
   printf "4 6 %6.2f %6.2f\n",s["C46"]/n["C46"], sqrt(ss["C46"]/n["C46"]);
   printf "5 6 %6.2f %6.2f\n",s["C56"]/n["C56"], sqrt(ss["C56"]/n["C56"]);
#
#   for(i=1;i<=6;i++) {
#     for(j=1;j<=6;j++) {
#       ii="C" i j;
#       if(i>j) printf "                 ";
#       else printf " %6.2f (%6.2f) ", s[ii]/n[ii],sqrt(ss[ii]/n[ii]);
#     }
#     printf "\n";
#  }
#
}
