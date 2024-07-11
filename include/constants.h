#ifndef _CONST_H
#define _CONST_H

#include <iostream>

#define _H_Theta(x) (x>0 ? x : 0 )
#define _sign(x) 	(x>0 ? 1 : -1)
#define _plus(x) 	(x>0 ? 1 : 0 )
#define _if(x) 		(x 	 ? 1 : 0 )



class ofstream_ps  : public std::ofstream {};  ///< output PostScript file
class ofstream_txt : public std::ofstream {};  ///< output text file


const double  gamma_0 = 1;		///< acid capacity number




#endif
