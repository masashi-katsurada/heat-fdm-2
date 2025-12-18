#include<stdio.h>
#include<math.h>

int main()
{
  int d;
  double deg, pi;

  pi = 4 * atan(1.0);
  /*３６０度が２πラジアンなので、1度はπ／１８０ラジアン*/
 deg = pi / 180 ;
 for(d= 0; d <= 90; d++){
  printf("%2d  %1.7f %1.7f\n", d,sin(d * deg),cos(d * deg));
 }
  return 0;
}
