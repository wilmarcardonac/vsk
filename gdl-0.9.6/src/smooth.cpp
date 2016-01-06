/***************************************************************************
                          smooth.cpp  -  mathematical GDL library function
                             -------------------
    begin                : 05 September 2014
    copyright            : (C) 2014 by Levan Loria  (with Alain Coulais)
    email                : alaingdl@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "includefirst.hpp"

#include "datatypes.hpp"
#include "envt.hpp"

#include "smooth.hpp"

namespace lib {

using namespace std;
using namespace antlr;

BaseGDL* smooth3_fun( EnvT* e ) {

  if ( e->NParam( ) < 2 ) e->Throw( "Incorrect number of arguments." );

  BaseGDL *p0 = e->GetPar( 0 );
  BaseGDL *p1 = e->GetPar( 1 );

  //---------------------------------------------------------
  //keywords
  //keyword NAN
  bool NaN_B = false;
  static int nan_kw = e->KeywordIx( "NAN" );
  if ( e->KeywordSet( nan_kw ) ) {
    NaN_B = true;
  }

  //---------------------------------------------------------


  //first argument
  if ( p0 == NULL ) e->Throw( "Expression must be an array in this context " + e->GetParString( 0 ) );

  SizeT dimArray = p0->Rank( );

  //cout<<dimArray<<endl;


  //---------------------------------------------------------
  //second argument
  if ( p1 == NULL ) e->Throw( "Variable is undefined: " + e->GetParString( 1 ) );
  if ( dimArray < p1->N_Elements( ) ) e->Throw( "Number of Array dimensions does not match number of Width dimensions" );
  DLong64GDL *W_tmp = static_cast<DLong64GDL*>(p1->Convert2( GDL_LONG64, BaseGDL::COPY ));
#ifdef _MSC_VER
  long *W = (long*) alloca( sizeof (long)*dimArray );
#else
  long W[dimArray];
#endif
  if ( p1->Rank( ) != 0 ) {
    if ( dimArray != p1->N_Elements( ) ) e->Throw( "Number of Array dimensions does not match number of Width dimensions" );
    for ( int i = 0; i < dimArray; i++ ) {
      if ( (*W_tmp)[i] % 2 == 0 ) W[i] = (*W_tmp)[i] + 1;
      else W[i] = (*W_tmp)[i];
      if ( W[i] < 0 || W[i] >= p0->Dim( i ) ) e->Throw( "Width must be nonnegative and smaller than array dimensions" );
    }
  } else {
    if ( (*W_tmp)[0] % 2 == 0 )
      for ( int i = 0; i < dimArray; i++ ) {
        W[i] = (*W_tmp)[0] + 1;
        if ( W[i] < 0 || W[i] >= p0->Dim( i ) ) e->Throw( "Width must be nonnegative and smaller than array dimensions" );
      } else
      for ( int i = 0; i < dimArray; i++ ) {
        W[i] = (*W_tmp)[0];
        if ( W[i] < 0 || W[i] >= p0->Dim( i ) ) e->Throw( "Width must be nonnegative and smaller than array dimensions" );
      }
  }

  //passed array
  DDoubleGDL *A = static_cast<DDoubleGDL *>(p0->Convert2( GDL_DOUBLE, BaseGDL::COPY ));

  //copy of our array , result - that will be returned
  DDoubleGDL *R;




  //---------------------------------------------------------
  //1D SMOOTH

  if ( dimArray == 1 ) {

    R = new DDoubleGDL( *A );

    //initialization of variables
    long ind1 = (W[0] - 1) / 2;
    long ind2 = p0->Dim( 0 ) - (W[0] + 1) / 2;
    long double s = 0;
    long double w_inv = (1.0) / W[0];
    long w2 = W[0] / 2;
    //temporary variables
    long index1, index2;
    //variables for NaN keyword
    long double sub, add, dvalue;

    if ( !NaN_B ) {
      // 1D without NaN keyword

      index1 = ind1 - w2;


      for ( long k = 0; k < W[0]; k++ ) {
        s += (*A)[k + index1];
      }

      (*R)[ind1] = w_inv * s;

      index1 = w2 + 1;
      index2 = W[0] - 1 - w2;


      for ( long ii = ind1 + 1; ii <= ind2; ii++ ) {
        s = s - (*A)[ii - index1] + (*A)[ii + index2];
        (*R)[ii] = w_inv * s;
      }

      //end of if(!NaN_B)
    } else {
      // 1D with NaN keyword

      long n = W[0];

      index1 = ind1 - w2;
      for ( long k = 0; k < W[0]; k++ ) {
        dvalue = (*A)[k + index1];
        if ( isnan( dvalue ) || isinf( dvalue ) ) {
          n--;
        } else {
          s += dvalue;
        }
      }

      (*R)[ind1] = 1.0 / n * s;
      index1 = w2 + 1;
      index2 = W[0] - 1 - w2;
      for ( long ii = ind1 + 1; ii <= ind2; ii++ ) {
        sub = (*A)[ii - index1];
        add = (*A)[ii + index2];


        if ( isnan( sub ) || isinf( sub ) ) {
          sub = 0;
        } else {
          n--;
        }


        if ( isnan( add ) || isinf( add ) ) {
          add = 0;
        } else {
          n++;

        }

        s = s - sub + add;

        (*R)[ii] = 1.0 / n * s;
      }

      //end of if(!NaN_B){...}else{...}
    }
    //end of if(dimArray=1)
  }

  //---------------------------------------------------------
  //2D SMOOTH

  if ( dimArray == 2 ) {


    //initialization of variables
    //main variables
    //indexes of width and height,W/2,...
    long I1 = (W[0] - 1) / 2; //ind1 - i
    long I2 = p0->Dim( 0 ) - (W[0] + 1) / 2; //ind2 - i
    long J1 = (W[1] - 1) / 2; //ind1 - j
    long J2 = p0->Dim( 1 ) - (W[1] + 1) / 2; //ind2 - j
    long w2_i = W[0] / 2; //w2 for i
    long w2_j = W[1] / 2; //w2 for j
    long double w_inv; //1/W
    long double s;
    long double sub, add;

    long N = p0->Dim( 0 );
    long M = p0->Dim( 1 );

    //other variables
    long num;

    //tmp variables
    long index, index1, index2;
    double dvalue;

    if ( !NaN_B ) {

      // 2D Smooth without NaN keyword

      DDoubleGDL *intermed_A = new DDoubleGDL( *A );

      //cout<<"first dimension"<<endl;
      if ( W[0] > 1 )
        for ( long i = 0; i < M; i++ ) {

          w_inv = (1.0) / W[0];
          s = 0.0;


          index = i * N + I1 - w2_i;
          for ( long k = 0; k < W[0]; k++ ) {
            s += (*A)[k + index];
            //s += (*A)[i*N + I1+k-w2_i];
          }

          (*intermed_A)[i * N + I1] = w_inv * s;

          index1 = i * N - 1 - w2_i;
          index2 = i * N + W[0] - 1 - w2_i;
          for ( long ii = I1 + 1; ii <= I2; ii++ ) {
            s = s - (*A)[ii + index1] + (*A)[ii + index2];
            //s = s - (*A)[i*N + ii-1 - w2_i] + (*A)[i*N + ii + W[0]-1 - w2_i];
            (*intermed_A)[i * N + ii] = w_inv * s;
          }
        }

      //second dimension
      R = new DDoubleGDL( *intermed_A );
      if ( W[1] < 2 ) return R;

      for ( long i = 0; i < N; i++ ) {

        w_inv = (1.0) / W[1];
        s = 0.0;

        index = i + (J1 - w2_j) * N;
        for ( long k = 0; k < W[1]; k++ ) {
          s += (*intermed_A)[k * N + index];
          //s += (*intermed_A)[i+(k+J1-w2_j)*N];
        }

        (*R)[i + N * J1] = w_inv * s;


        index1 = i - (1 + w2_j) * N;
        index2 = i + (W[1] - 1 - w2_j) * N;
        for ( long ii = J1 + 1; ii <= J2; ii++ ) {
          s = s - (*intermed_A)[N * ii + index1] + (*intermed_A)[N * ii + index2];
          //s = s - (*intermed_A)[i+ N * (ii-1 - w2_j)] + (*intermed_A)[i+ N * (ii + W[1]-1 - w2_j)];
          (*R)[i + N * ii] = w_inv * s;
        }
      }


      //end of if(!NaN_B){...}
    } else {

      //2D Smooth with NaN keyword
      DDoubleGDL *intermed_A = new DDoubleGDL( *A );

      //cout<<"first dimension"<<endl;
      if ( W[0] > 1 )
        for ( long i = 0; i < M; i++ ) {

          s = 0.0;
          num = W[0];

          for ( long k = 0; k < W[0]; k++ ) {
            dvalue = (*A)[i * N + I1 + k - w2_i];

            if ( isnan( dvalue ) || isinf( dvalue ) ) {
              num--;
            } else {
              s += dvalue;
            }
          }

          (*intermed_A)[i * N + I1] = 1.0 / num * s;

          for ( long ii = I1 + 1; ii <= I2; ii++ ) {

            sub = (*A)[i * N + ii - 1 - w2_i];
            add = (*A)[i * N + ii + W[0] - 1 - w2_i];

            if ( isnan( sub ) || isinf( sub ) ) {
              sub = 0;
            } else {
              num--;
            }


            if ( isnan( add ) || isinf( add ) ) {
              add = 0;
            } else {
              num++;
            }

            s = s - sub + add;

            (*intermed_A)[i * N + ii] = 1.0 / num * s;

          }
        }

      //second dimension
      R = new DDoubleGDL( *intermed_A );
      if ( W[1] < 2 ) return R;

      for ( long i = 0; i < N; i++ ) {

        s = 0.0;
        num = W[1];


        for ( long k = 0; k < W[1]; k++ ) {
          dvalue = (*intermed_A)[i + (k + J1 - w2_j) * N];
#ifdef _MSC_VER
          if ( isnan( dvalue ) ) {
#else
          if ( isnan( dvalue ) ) {
#endif
            num--;
          } else {
            s += dvalue;
          }
        }


        (*R)[i + N * J1] = 1.0 / num * s;

        for ( long ii = J1 + 1; ii <= J2; ii++ ) {
          sub = (*intermed_A)[i + N * (ii - 1 - w2_j)];
          add = (*intermed_A)[i + N * (ii + W[1] - 1 - w2_j)];

          if ( isnan( sub ) || isinf( sub ) ) {
            sub = 0;
          } else {
            num--;
          }


#ifdef _MSC_VER
          if ( isnan( (long double) (add || isinf( add )) ) ) {
#else
          if ( isnan( add || isinf( add ) ) ) {
#endif
            add = 0;
          } else {
            num++;
          }

          s = s - sub + add;

          (*R)[i + N * ii] = 1.0 / num * s;

        }
      }


    }// end of if(!NaN_B){...}else{...}   ---- SMOOTH 2D

  }//end of if(dimArray == 2)


  //---------------------------------------------------------
  //type conversion

  //BYTE
  if ( p0->Type( ) == GDL_BYTE ) {
    return R->Convert2( GDL_BYTE, BaseGDL::CONVERT );
  } else


    //INT
    if ( p0->Type( ) == GDL_INT ) {
    return R->Convert2( GDL_INT, BaseGDL::CONVERT );
  } else

    //LONG
    if ( p0->Type( ) == GDL_LONG ) {
    return R->Convert2( GDL_LONG, BaseGDL::CONVERT );
  } else

    //UINT
    if ( p0->Type( ) == GDL_UINT ) {
    return R->Convert2( GDL_UINT, BaseGDL::CONVERT );
  } else


    //ULONG
    if ( p0->Type( ) == GDL_ULONG ) {
    return R->Convert2( GDL_ULONG, BaseGDL::CONVERT );
  } else


    //LONG64
    if ( p0->Type( ) == GDL_LONG64 ) {
    return R->Convert2( GDL_LONG64, BaseGDL::CONVERT );
  } else


    //ULONG64
    if ( p0->Type( ) == GDL_ULONG64 ) {
    return R->Convert2( GDL_ULONG64, BaseGDL::CONVERT );
  } else


    //FLOAT
    if ( p0->Type( ) == GDL_FLOAT ) {
    return R->Convert2( GDL_FLOAT, BaseGDL::CONVERT );
  } else


    //DOUBLE
    if ( p0->Type( ) == GDL_DOUBLE ) {
    return R;
  } else
 {
    e->Throw( "SMOOTH DOESN'T SUPPORT THIS TYPE YET : " + p0->TypeStr( ) );
  }

  //---------------------------------------------------------

}

BaseGDL* smooth2_fun( EnvT* e ) {

  //cout<<"Not ready"<<endl;

  BaseGDL *p0 = e->GetPar( 0 );
  BaseGDL *p1 = e->GetPar( 1 );

  if ( p0 == NULL || p1 == NULL ) e->Throw( "Incorrect number of arguments." );

  bool NaN_B = false;
  static int nan_kw = e->KeywordIx( "NAN" );
  if ( e->KeywordSet( nan_kw ) ) {
    NaN_B = true;
  }

  SizeT dimArray = p0->Rank( );
  SizeT N = p0->N_Elements( );


  DLongGDL * W = static_cast<DLongGDL*>(p1->Convert2( GDL_LONG, BaseGDL::COPY ));
  DDoubleGDL *A = static_cast<DDoubleGDL *>(p0->Convert2( GDL_DOUBLE, BaseGDL::COPY ));

  DDoubleGDL *R = new DDoubleGDL( *A );

  //cout<<N<<endl;
  //cout<<(*W)[0]<<endl;

  if ( p1->N_Elements( ) == 0 ) e->Throw( "Variable is undefined: <UNDEFINED>." );

  if ( (*W)[0] < 0 || (*W)[0] >= N ) e->Throw( "Width must be nonnegative and smaller than array dimensions" );
  if ( (*W)[0] <= 1 && dimArray == 0 ) return R;


  //cout<<ind1<<" "<<ind2<<endl;
  if ( dimArray == 1 )


    if ( !NaN_B ) {
      if ( (*W)[0] % 2 == 0 ) {
        (*W)[0] += 1;
      }

      long ind1 = ((*W)[0] - 1) / 2;
      long ind2 = N - ((*W)[0] + 1) / 2;
      long double s = 0.0;
      long double w_inv = (1.0) / (*W)[0];
      long w2 = (*W)[0] / 2;

      long index1 = ind1 - w2;
      if ( p1->N_Elements( ) > 1 ) e->Throw( "Number of Array dimensions does not match number of Width dimensions." );
      for ( long k = 0; k < (*W)[0]; k++ ) {
        s += (*A)[k + index1];
      }
      (*R)[ind1] = w_inv * s;

      index1 = w2 + 1;
      long index2 = (*W)[0] - 1 - w2;
      for ( ind1++; ind1 <= ind2; ind1++ ) {
        s = s - (*A)[ind1 - index1] + (*A)[ind1 + index2];
        (*R)[ind1] = w_inv * s;
      }
    } else {
//      double missing_data = 0.0; //unused


      if ( (*W)[0] % 2 == 0 ) {
        (*W)[0] += 1;
      }

      long n = (*W)[0];

      long ind1 = ((*W)[0] - 1) / 2;
      long ind2 = N - ((*W)[0] + 1) / 2;
      long double s = 0.0;
      long double w_inv = (1.0) / (*W)[0];
      long w2 = (*W)[0] / 2;

      for ( long k = 0; k < (*W)[0]; k++ ) {
        long double dvalue = (*A)[ind1 + k - w2];
        if ( isnan( dvalue ) ) {
          n--;
        } else {
          s += dvalue;
        }
      }

      (*R)[ind1] = 1.0 / n * s;
      for ( ind1++; ind1 <= ind2; ind1++ ) {
        long double sub = (*A)[ind1 - 1 - w2];
        long double add = (*A)[ind1 + (*W)[0] - 1 - w2];

        if ( isnan( sub ) ) {
          sub = 0;
        } else {
          n--;
        }


        if ( isnan( add ) ) {
          add = 0;
        } else {
          n++;
        }

        s = s - sub + add;

        (*R)[ind1] = 1.0 / n * s;
      }
    }

  if ( dimArray == 2 ) {

    long n1, n2;
    long n_elem = p1->N_Elements( );
    long N = p0->Dim( 0 );
    long M = p0->Dim( 1 );
    //cout<<N<<" "<<M<<endl;
    if ( p1->Rank( ) > 1 ) e->Throw( "Number of Array dimensions does not match number of Width dimensions." );
    if ( p1->Rank( ) == 1 ) {
      if ( n_elem != 2 ) e->Throw( "Number of Array dimensions does not match number of Width dimensions." );
      if ( (*W)[0] >= N || (*W)[1] >= N || (*W)[0] >= M || (*W)[1] >= M || (*W)[0] < 0 || (*W)[1] < 0 ) e->Throw( "Width must be nonnegative and smaller than array dimensions" );
      if ( (*W)[0] < 2 && (*W)[1] < 2 ) return R;
      n1 = (*W)[0];
      n2 = (*W)[1];
      if ( n1 % 2 == 0 ) n1++;
      if ( n2 % 2 == 0 ) n2++;

    }

    if ( p1->Rank( ) == 0 ) {

      n1 = (*W)[0];
      if ( n1 % 2 == 0 ) n1++;
      n2 = n1;

      if ( (*W)[0] < 0 || (*W)[0] >= N ) e->Throw( "Width must be nonnegative and smaller than array dimensions" );
      if ( (*W)[0] < 2 ) return R;

    }
    if ( !NaN_B ) {



      DDoubleGDL *intermed_A = new DDoubleGDL( *A );

      //cout<<"first dimension"<<endl;
      if ( n1 > 1 )
        for ( long i = 0; i < M; i++ ) {
          int ind1 = (n1 - 1) / 2;

          int ind2 = N - (n1 + 1) / 2;
          int w2 = n1 / 2;
          long double w_inv = (1.0) / n1;
          long double s = 0.0;


          for ( long k = 0; k < n1; k++ ) {
            s += (*A)[i * N + ind1 + k - w2];
          }
          (*intermed_A)[i * N + ind1] = w_inv * s;

          for ( ind1++; ind1 <= ind2; ind1++ ) {
            s = s - (*A)[i * N + ind1 - 1 - w2] + (*A)[i * N + ind1 + n1 - 1 - w2];
            (*intermed_A)[i * N + ind1] = w_inv * s;

          }
        }

      //second dimension
      R = new DDoubleGDL( *intermed_A );
      if ( n2 < 2 ) return R;
      for ( long i = 0; i < N; i++ ) {
        int ind1 = (n2 - 1) / 2;
        int ind2 = M - (n2 + 1) / 2;
        int w2 = n2 / 2;
        long double w_inv = (1.0) / n2;
        long double s = 0.0;


        for ( long k = 0; k < n2; k++ ) {
          s += (*intermed_A)[i + (k + ind1 - w2) * N];
        }
        (*R)[i + N * ind1] = w_inv * s;

        for ( ind1++; ind1 <= ind2; ind1++ ) {
          s = s - (*intermed_A)[i + N * (ind1 - 1 - w2)] + (*intermed_A)[i + N * (ind1 + n2 - 1 - w2)];
          //scout<<s<<"->"<<i*N + ind1<<endl;
          (*R)[i + N * ind1] = w_inv * s;

        }
      }

    } else {

      DDoubleGDL *intermed_A = new DDoubleGDL( *A );

      //cout<<"first dimension"<<endl;
      if ( n1 > 1 )
        for ( long i = 0; i < M; i++ ) {
          long ind1 = (n1 - 1) / 2;
          long ind2 = N - (n1 + 1) / 2;
          long w2 = n1 / 2;
          long double w_inv = (1.0) / n1;
          long double s = 0.0;

          long nn1 = n1;
          //long nn2 = n2;

          for ( long k = 0; k < n1; k++ ) {
            long double dvalue = (*A)[i * N + ind1 + k - w2];
            if ( isnan( dvalue ) ) {
              nn1--;
            } else {
              s += dvalue;
            }
          }
          (*intermed_A)[i * N + ind1] = 1.0 / nn1 * s;

          for ( ind1++; ind1 <= ind2; ind1++ ) {

            long double sub = (*A)[i * N + ind1 - 1 - w2];
            long double add = (*A)[i * N + ind1 + n1 - 1 - w2];

            if ( isnan( sub ) ) {
              sub = 0;
            } else {
              nn1--;
            }


            if ( isnan( add ) ) {
              add = 0;
            } else {
              nn1++;
            }

            s = s - sub + add;

            (*intermed_A)[i * N + ind1] = 1.0 / nn1 * s;

          }
        }

      //second dimension
      R = new DDoubleGDL( *intermed_A );
      if ( n2 < 2 ) return R;
      for ( long i = 0; i < N; i++ ) {
        long ind1 = (n2 - 1) / 2;
        long ind2 = M - (n2 + 1) / 2;
        long w2 = n2 / 2;
        long double w_inv = (1.0) / n2;
        long double s = 0.0;


        //long nn1 = n1;
        long nn2 = n2;


        for ( long k = 0; k < n2; k++ ) {
          long double dvalue = (*intermed_A)[i + (k + ind1 - w2) * N];
          if ( isnan( dvalue ) ) {
            nn2--;
          } else {
            s += dvalue;
          }
        }
        (*R)[i + N * ind1] = 1.0 / nn2 * s;

        for ( ind1++; ind1 <= ind2; ind1++ ) {
          long double sub = (*intermed_A)[i + N * (ind1 - 1 - w2)];
          long double add = (*intermed_A)[i + N * (ind1 + n2 - 1 - w2)];

          if ( isnan( sub ) ) {
            sub = 0;
          } else {
            nn2--;
          }


          if ( isnan( add ) ) {
            add = 0;
          } else {
            nn2++;
          }

          s = s - sub + add;

          (*R)[i + N * ind1] = 1.0 / nn2 * s;

        }
      }

    }




  }



  return R;
}
} // namespace

  