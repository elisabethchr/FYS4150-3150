//definitions for functions within Initialize
#include "initialize.h"

mat Initialize::Initialvelocity(mat vel, int dim, int N)
{
// JPL values are AU/day so multiply with 365.25
    double vx0 = vel(0,0) = -5.428888690270241E-03*365.25;      // AU/yr
    double vy0 = vel(1,0) = 1.636353485946535E-02*365.25;       // AU/yr
    double vz0 = vel(2,0) = -4.491683144318728E-07*365.25;      // AU/yr
    cout <<"Initial velocity in x, y, z direction:" << endl;
    cout << vx0 << " " << vy0 <<" " << vz0 << endl;
    return vel;
}

mat Initialize::Initialposition(mat pos, int dim, int N) //(, string obj)
{   /*pos_sun(0, 0) = obj(0, 0)
*/
    double x0 = pos(0,0) = 9.528047055398201E-01;       // AU
    double y0 = pos(1,0) = 3.053612869840809E-01;       // AU
    double z0 = pos(2,0) = -9.272902073041313E-05;      // AU
    cout <<"Initial position in x, y, z direction:" << endl;
    cout << x0 << " " << y0 <<" " << z0 << endl;
    return pos;
}
