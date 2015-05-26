#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <cassert>
#include "common.h"

template<size_t N, typename T, typename Real>
class ExpansionBaseT
{
  protected:
    using vector_type = std::array<Real,N>;
    using matrix_type = std::array<vector_type,N>;
  public:
    using storage     = std::array<T,N>;
    static constexpr auto size() { return N; }
};

template<size_t N, typename T, typename Real>
class ExpansionT;

#if 1
#define PIF_
#else
#define DGF_
#endif

#if defined PIF_  /* not PIFt_ */
template<typename T, typename Real>
class ExpansionT<1,T,Real> : public ExpansionBaseT<1,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<1,T,Real>;
    static constexpr size_t N = 1;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    {
      return Real{1};
    }

    static constexpr auto weight(const size_t i) 
    {
      return Real{1};
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    { 
      return Real{0.5};
    }
    static constexpr auto maxAbsMEV() 
    {
      return Real{2};
    }
    static constexpr auto maxAbsPEV() 
    { 
      return Real{0.5};
    }
    static constexpr auto oneVec(const size_t i)
    {
      return Real{1};
    }
};
template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  /* PIF */
  protected:
    using base_type = ExpansionBaseT<2,T,Real>;
    static constexpr size_t N = 2;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {3.0,  0.4641016151377546},
        {-6.464101615137754,  3.0}
      };
      return matrix[i][j]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] = {0.5,0.5};
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.250000000000000000000,  0},
        {0.53867513459481288225, 0.250000000000000000000}
      };
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV = 3.5;
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV = 0.25;
      return maxAbsPEV; 
    }
    static constexpr auto oneVec(const size_t i)
    {
      constexpr Real vec[N] = 
      {-0.366025403784438646764,1.366025403784438646764};
      return vec[i];
    }
};
template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  /* PIF */
  protected:
    using base_type  = ExpansionBaseT<3,T,Real>;
    static constexpr size_t N = 3;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {5.0,  1.1639777949432226, -0.1639777949432225},
        {-5.727486121839514,  2.0, 0.7274861218395141},
        {10.163977794943223,  -9.163977794943223, 5.0}
      };
      return matrix[i][j]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.2777777777777778,
        0.4444444444444444,
        0.2777777777777778
      };
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)  
    { 
      constexpr Real weight[] =
      {0.30026319498086457,0.2222222222222222,-0.022485417203086805};
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.13888888888888888889, 0, 0},
        {0.30026319498086459244, 0.22222222222222222222,  0},
        {0.26798833376246945173, 0.48042111196938334790,  0.13888888888888888889}
      };
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{5.1};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.22222};
      return maxAbsPEV; 
    }

    static constexpr auto prolongateMatrix0(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {1.228830557701236147530,-0.312163891034569480863,0.083333333333333333333},
        {0.531081945517284740432,0.583333333333333333333,-0.1144152788506180737649},
        {0.0833333333333333333333,0.978830557701236147530,-0.0621638910345694808632}
      };
      return matrix[i][j];
    }
    static constexpr auto prolongateMatrix1(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {-0.0621638910345694808632,0.978830557701236147530,0.0833333333333333333333},
        {-0.1144152788506180737649,0.583333333333333333333,0.531081945517284740432},
        {0.083333333333333333333,-0.312163891034569480863,1.228830557701236147530}
      };
      return matrix[i][j];
    }
    static constexpr auto midVec(const size_t i)
    {
      constexpr Real vec[N] =
      {
        0., 1., 0.
      };
      return vec[i];
    }
    static constexpr auto oneVec(const size_t i)
    {
      constexpr Real vec[N] = 
      {0.187836108965430519137,-0.666666666666666666666,1.478830557701236147530};
      return vec[i];
    }
};
template<typename T, typename Real>
class ExpansionT<5,T,Real> : public ExpansionBaseT<5,T,Real>
{
  /* PIF */
  protected:
    static constexpr size_t N = 5;
    using base_type  = ExpansionBaseT<N,T,Real>;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {11.18330013267037774,3.13131216201181084,-0.758731795980807391,0.239101223353686050,-0.0543147605653389237},
        {-9.44759960151614989,2.81669986732962226,2.21788633227481812,-0.557122620293797301,0.118357949604667387},
        {6.42011650355933662,-6.22012045466975167,2.00000000000000000,1.86599528883177949,-0.315991337721364445},
        {-8.01592078481097030,6.19052235495304182,-7.39311627638238017,2.81669986732962226,1.55003676630984697},
        {22.4209150259060944,-16.1933902399992350,15.4154432215698509,-19.0856011786573598,11.18330013267037774}
      };
      return matrix[i][j]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.11846344252809503,
        0.2393143352496833,
        0.28444444444444444,
        0.23931433524968349,
        0.1184634425280951
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.0592317212640472719,0,0,0,0},
        {0.128151005670045283,0.119657167624841617,0,0,0},
        {0.1137762880042246025,0.260004651680641519,0.142222222222222222,0,0},
        {0.121232436926864147,0.228996054578999877,0.309036559064086645,0.119657167624841617,0},
        {0.1168753295602285452,0.244908128910495419,0.273190043625801489,0.258884699608759272,0.0592317212640472719}
      };
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{8.6};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.15};
      return maxAbsPEV; 
    }
    static constexpr auto nodeMatrix(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {0.213260005841349104302,1.71754726228747258129,9.3031177500426426728,29.2274326029943329589,54.3373967987902398687},
        {-0.74002965114900785695,-5.7585639443127711100,-30.1637945111551133825,-92.666703922546651957,-170.315662745074820114},
        {1.42879537949198108031,10.1954385103677075047,49.5333333333333333333,145.324986458118124091,261.184112985355520657},
        {-2.17509816702138543022,-12.4770986424139843478,-52.3901873373854831267,-142.564880053487944785,-247.782697697696973869},
        {2.27307243283706310256,7.3226768140715753718,24.7175307651646205030,61.6791649149221396919,103.576850658626033458}
      };
      return matrix[i][j];
    }

    static constexpr auto prolongateMatrix0(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {1.25614688083180737647,-0.40249692547716509839,0.226542265055292140728,-0.111885279656663511252,0.031693059246729092442},
        {0.445793658412247956549,0.74923312241633579261,-0.290802660952930342018,0.132214204855823862172,-0.036438324731477269318},
        {-0.032195241140255091699,0.96272564592129313619,0.095833333333333333333,-0.035663443357677330307,0.0092997052433059524822},
        {-0.071182755775665506236,0.442505803988867987154,0.76339820471556343181,-0.177008133472200604271,0.0422868805434346915484},
        {-0.0147464014247542751011,0.073003992669698522501,0.98975108007096365837,-0.061303468206632007628,0.0132947968907241018546}
      };
      return matrix[i][j];
    }
    static constexpr auto prolongateMatrix1(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {0.0132947968907241018546,-0.061303468206632007628,0.98975108007096365837,0.073003992669698522501,-0.0147464014247542751011},
        {0.0422868805434346915484,-0.177008133472200604271,0.76339820471556343181,0.442505803988867987154,-0.071182755775665506236},
        {0.0092997052433059524822,-0.035663443357677330307,0.095833333333333333333,0.96272564592129313619,-0.032195241140255091699},
        {-0.036438324731477269318,0.132214204855823862172,-0.290802660952930342018,0.74923312241633579261,0.445793658412247956549},
        {0.031693059246729092442,-0.111885279656663511252,0.226542265055292140728,-0.40249692547716509839,1.25614688083180737647}
      };
      return matrix[i][j];
    }
    static constexpr auto weight_half(const size_t i)  
    { 
      constexpr Real weight[] =
      {0.1137762880042246,0.2600046516806415,0.14222222222222236,-0.020690316430958255,0.004687154523869946};
      return weight[i];
    }
    static constexpr auto midVec(const size_t i)
    {
      constexpr Real vec[N] =
      {
        0., 0., 1., 0., 0.
      };
      return vec[i];
    }
    static constexpr auto oneVec(const size_t i)
    {
      constexpr Real vec[N] = 
      {0.076358661795812900484,-0.267941652223387509304,0.53333333333333333333,-0.89315839200007173733,1.55140804909431301281};
      return vec[i];
    }
};
template<typename T, typename Real>
class ExpansionT<7,T,Real> : public ExpansionBaseT<7,T,Real>
{
  /* PIF */
  protected:
    static constexpr size_t N = 7;
    using base_type  = ExpansionBaseT<N,T,Real>;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {20.162475126453963,5.939649595901784,-1.5713450065769194,0.6117568698401676,-0.2662682556100125,0.10823443910941748,-0.02751051928120109},
        {-15.629383503722168,4.443146050695284,4.125140303481284,-1.2701515428956040,0.5100682216357629,-0.20014625792437366,0.05010533623170008},
        {8.625187972409475,-8.605076196977838,2.3943788228507531,3.351834730017060,-1.0413621688017852,0.3736448491420073,-0.09029602574207080},
        {-7.258551263394934,5.727243007023915,-7.245304978834095,2.0000000000000000,3.062096190825881,-0.8500070065068292,0.18952405088606242},
        {8.182614584790622,-5.956895021669846,5.830119814503291,-7.930862828666295,2.3943788228507531,3.021826024450000,-0.5328694133609241},
        {-12.929826694060736,9.086438359314942,-8.131847279980262,8.558125387465308,-11.746919361825783,4.443146050695284,2.749662145893132},
        {40.35246077218913,-27.93020953980242,24.129196096015008,-23.429578355405607,25.43427284698192,-33.76162469659479,20.162475126453963}
      };
      return matrix[i][j]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.06474248308443512,
        0.1398526957446387,
        0.19091502525255918,
        0.20897959183673512,
        0.1909150252525592,
        0.13985269574463888,
        0.06474248308443513
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {{0.03237124154221742,0,0,0,0,0,0},{0.07004354137872608,0.06992634787231917,0,0,0,0,0},{0.06215393578734986,0.15200552205783099,0.09545751262627974,0,0,0,0},{0.06633292861768470,0.13359576922388229,0.20770188076597078,0.10448979591836735,0,0,0},{0.06366576746880693,0.14380627590344877,0.18220246265408429,0.22735483605218653,0.09545751262627974,0,0},{0.06548633327460677,0.13720685187790830,0.19631211718445561,0.19962996905329136,0.20750503183140724,0.06992634787231917,0},{0.06428437356068539,0.14145951478168444,0.18773996647887383,0.21411332539996004,0.18328182138013593,0.15130371302782220,0.03237124154221742}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{12.1};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.11};
      return maxAbsPEV; 
    }
    static constexpr auto nodeMatrix(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {0.119754539292883770824,1.22728019371609711879,10.4575118331634874210,59.758130571637751289,221.822382804114009194,540.70851913149474530,881.03018324838731946},
        {-0.41835867464874851868,-4.2412796572008277475,-35.651064814775981503,-201.262030327699632498,-740.30470969359352825,-1793.78716313596831775,-2913.54990834164420749},
        {0.81929749536172830115,8.1119477463172599059,66.265956470439013599,365.023534153547473868,1318.94049172829736763,3159.48999676330457420,5101.1851673164638649},
        {-1.30010385070212712168,-12.2809301491938185975,-95.170703437054145911,-502.657142857142857143,-1764.56772094090628785,-4152.27192312391146419,-6643.24378333339699150},
        {1.85025477095430382384,15.8344513639069835156,111.529774812171960644,550.88130274216226410,1854.22609525154686867,4257.46607507546489964,6729.2444447210924098},
        {-2.4239515078379212814,-16.4092510233253538662,-97.662893729409545089,-438.44243600366021949,-1399.98676532939864340,-3123.93628492915727242,-4870.8890047931245865},
        {2.35310722757988102594,8.7577815257796596709,41.231418865465210839,167.698641721155219876,510.87022617994021401,1113.33078021877283523,1717.22290118222219128}
      };
      return matrix[i][j];
    }
    static constexpr auto oneVec(const size_t i)
    {
      constexpr Real vec[N] = 
      {0.041115148862905928075,-0.144070103612068846929,0.28405414676522996668,-0.45714285714285714286,0.67210786192236178693,-0.97072669650612219065,1.57466249971055049874};
      return vec[i];
    }
    static constexpr auto prolongateMatrix0(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {1.26456599983916812566,-0.43234534056732839735,0.28195038322772238642,-0.18835192783047360857,0.11616036592589388227,-0.058707057505859117182,0.016727576910876728749},{0.42109278607973913381,0.79921466637368048938,-0.35358713521420274638,0.21611611605426510718,-0.128770565134950087797,0.064061427783256400240,-0.0181272959417882964250},{-0.045835148885397256768,0.91504833621469057431,0.18927447399255919105,-0.091573017089823911897,0.050713212378911497386,-0.024458203653981845121,0.0068303470430417510342},{-0.039817436193590984945,0.231802523127790484963,0.94641623008848683648,-0.204017857142857142857,0.098371864929076236398,-0.045095555299286550475,0.0123402304903811204335},{0.028074588343309993459,-0.128950219890518757533,0.83865263648751860443,0.35150041090653656033,-0.129770111565080330169,0.055181913285116808573,-0.0146892175668828790900},{0.0236298249953582552825,-0.099063291911618084239,0.34901105360735571686,0.85515214425013626248,-0.180422155391203671529,0.069658396009964208262,-0.0179659715599926871149},{0.0048012954939119962954,-0.019389837104433031921,0.058094207358141521295,0.99418078765695046114,-0.051239127713317722200,0.0181032418161813167151,-0.0045505675074345413282}
      };
      return matrix[i][j];
    }
    static constexpr auto prolongateMatrix1(const size_t i, const size_t j)
    {
      constexpr Real matrix[N][N] = 
      {
        {-0.0045505675074345413282,0.0181032418161813167151,-0.051239127713317722200,0.99418078765695046114,0.058094207358141521295,-0.019389837104433031921,0.0048012954939119962954},{-0.0179659715599926871149,0.069658396009964208262,-0.180422155391203671529,0.85515214425013626248,0.34901105360735571686,-0.099063291911618084239,0.0236298249953582552825},{-0.0146892175668828790900,0.055181913285116808573,-0.129770111565080330169,0.35150041090653656033,0.83865263648751860443,-0.128950219890518757533,0.028074588343309993459},{0.0123402304903811204335,-0.045095555299286550475,0.098371864929076236398,-0.204017857142857142857,0.94641623008848683648,0.231802523127790484963,-0.039817436193590984945},{0.0068303470430417510342,-0.024458203653981845121,0.050713212378911497386,-0.091573017089823911897,0.18927447399255919105,0.91504833621469057431,-0.045835148885397256768},{-0.0181272959417882964250,0.064061427783256400240,-0.128770565134950087797,0.21611611605426510718,-0.35358713521420274638,0.79921466637368048938,0.42109278607973913381},{0.016727576910876728749,-0.058707057505859117182,0.11616036592589388227,-0.18835192783047360857,0.28195038322772238642,-0.43234534056732839735,1.26456599983916812566}
      };
      return matrix[i][j];
    }
    static constexpr auto weight_half(const size_t i)  
    { 
      constexpr Real weight[N] =
      {0.06633292861768417,0.13359576922388383,0.20770188076596852,0.10448979591836771,-0.016786855513409416,0.006256926520754824,-0.0015904455332495793};
      return weight[i];
    }
    static constexpr auto midVec(const size_t i)
    {
      constexpr Real vec[N] =
      {
        0., 0., 0., 1., 0., 0., 0.
      };
      return vec[i];
    }
};
#elif defined DGF_  /* not PIF_ */
template<typename T, typename Real>
class ExpansionT<1,T,Real> : public ExpansionBaseT<1,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<1,T,Real>;
    static constexpr size_t N = 1;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    {
      return Real{1};
    }

    static constexpr auto weight(const size_t i) 
    {
      return Real{1};
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    { 
      return Real{1};
    }
    static constexpr auto maxAbsMEV() 
    {
      return Real{1};
    }
    static constexpr auto maxAbsPEV() 
    { 
      return Real{1};
    }
};
template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  /* DGF */
  protected:
    using base_type = ExpansionBaseT<2,T,Real>;
    static constexpr size_t N = 2;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {2.0,  0.7320508075688785},
        {-2.732050807568876,  1.9999999999999958}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] = {0.5,0.5};
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.3333333333333329,  0},
        {0.4553418012614798,  0.33333333333333365}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV = 2.5;
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV = 0.33;
      return maxAbsPEV; 
    }
};
template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  /* DGF */
  protected:
    using base_type  = ExpansionBaseT<3,T,Real>;
    static constexpr size_t N = 3;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {4.000000000000003, 1.6147844564602516, -0.2909944487358082},
        {-3.509240285287669,  1.0000000000000007, 1.009240285287668},
        {2.290994448735803, -5.6147844564602405, 3.9999999999999774}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.2777777777777778,
        0.4444444444444444,
        0.2777777777777778
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.16111111111111154, 0,0},
        {0.2724854172030868,  0.2777777777777769, 0},
        {0.29021055598469203, 0.43597666752493847,  0.16111111111111140}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{4.1};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.28};
      return maxAbsPEV; 
    }
};
template<typename T, typename Real>
class ExpansionT<4,T,Real> : public ExpansionBaseT<4,T,Real>
{
  /* DGF */
  protected:
    static constexpr size_t N = 4;
    using base_type  = ExpansionBaseT<N,T,Real>;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {6.738612787525834, 2.5779942727073806, -0.6995574654955786,  0.1612563383245323},
        {-5.324832600342981,  1.2613872124741736, 1.941340462561435,  -0.3731446166705182},
        {2.5339048478804855,  -3.9413404625614348,  1.2613872124741656, 1.3751045941403923},
        {-2.161256338324527,  4.750469319393827,  -9.982795494470142, 6.738612787525786}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.17392742256872748, 
        0.32607257743127305, 
        0.3260725774312732, 
        0.17392742256872748
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.0950400941860567,  0,0,0},
        {0.17720653136163217, 0.19067419152822931,  0,0},
        {0.17810350811242565, 0.32631510322115087,  0.19067419152822807,  0},
        {0.1694061893528294,  0.3339017452341196, 0.33222012702401943,  0.09504009418605698}

      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{5.8};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.20};
      return maxAbsPEV; 
    }
};
#endif


template<size_t ORDER, typename PDE>
class ODESolverT
{
  public:
    using Real           = typename PDE::Real;
    using Vector         = typename PDE::Vector;
  private:
    using Expansion      = ExpansionT<ORDER,Vector,Real>;
    using range_iterator = make_range_iterator<size_t>;


    PDE _pde;
    Real _time;
    bool _verbose;
    typename Expansion::storage _x, _rhs;
    Vector _y0;

    static constexpr Real omegaCFL = 0.8;

  public:
      auto expansionRange() const 
      {
        return make_range_iteratorT<0,Expansion::size()>{};
      }

    auto time() const {return _time;}

    ODESolverT(const PDE &pde) : _pde{pde}, _time{0}
    {
      for (auto k : expansionRange())
      {
        _x  [k].resize(_pde.resolution());
        _rhs[k].resize(_pde.resolution());
      }
    };

    PDE& pde() { return _pde; }
    const PDE& pde() const { return _pde; }

    void rhs(const Vector &u0)
    {
      using std::get;
      auto x = _x;

      /* compute RHS */
      for (auto k : expansionRange())
      {

        for (auto v : make_zip_iterator(x[k], u0))
          get<0>(v) += get<1>(v);
        _pde.compute_rhs(_rhs[k], x[k]);

        assert(_x[k].size() == u0.size());
        for (auto l : expansionRange())
          for (auto v : make_zip_iterator(_rhs[k], _x[l]))
            get<0>(v) -= Expansion::matrix(k,l) * get<1>(v);
      }

      /* precondition RHS */
      for (auto i : range_iterator{0,u0.size()})
      {
        std::array<Real,Expansion::size()> tmp;
        for (auto k : expansionRange())
        {
          tmp[k] = 0;
          for (auto l : expansionRange())
           tmp[k] += Expansion::preconditioner(k,l)*_rhs[l][i];
        }
        for (auto k : expansionRange())
          _rhs[k][i] = tmp[k];
      }
    }

    void iterate(const Vector &u0, int n)
    {
      using std::get;

      const Real omega = omegaCFL/(Expansion::maxAbsPEV()*(Expansion::maxAbsMEV() + _pde.AbsEV()));

      static decltype(_x) y0, y1,tmp;
      y0  = _x;

      rhs(u0);
      static auto res = _x;

#if 0
#define WP   /* prefered for high resolution & large time step, use scaleing 1x for nstate with cfl */
#elif 1
#define OPT
#endif

#ifdef WP
      auto scale = Real{1}/(2*n+1);
#elif defined OPT
      auto frac = Real{1}/(1+n)/(2+n)/(3+4*n);
      auto omega0 = frac*6;
      auto omega1 = frac*3*n*(3+n);
      auto omegak = [n,frac](const int k) 
      {
        return frac*3*(1-k+n)*(2+k+n);
      };
#endif

      for (auto k : expansionRange())
      {
        for (auto v : make_zip_iterator(_x[k], _rhs[k],res[k]))
        {
          auto&   x = get<0>(v);
          auto& rhs = get<1>(v);
#ifdef OPT
          auto& r = get<2>(v);
          r = x*omega0 + (3*x + 4.0*omega*rhs)*omega1;
#elif defined WP
          auto& r = get<2>(v);
          r = (3*x + 4.0*omega*rhs)*scale;
#endif
          x = x + 2.0*omega*rhs;
        }
        y1[k] = _x[k];
      }



      for (int i = 2; i <= n; i++)
      {
        rhs(u0);
        for (auto k : expansionRange())
        {
          for (auto v : make_zip_iterator(_x[k], y1[k], y0[k], _rhs[k],res[k]))
          {
            auto&   x = get<0>(v);
            auto&  y1 = get<1>(v);
            auto&  y0 = get<2>(v);
            auto& rhs = get<3>(v);
#ifdef OPT
            x = 2*y1 - y0 + 4*omega*rhs;
            auto& r = get<4>(v);
            r += 2*omegak(i)*x;
#elif defined WP
            x = 2*y1 - y0 + 4*omega*rhs;
            auto& r = get<4>(v);
            r += 2*scale*x;
#else
            const auto a = (Real{2}*i-1)/i;
            const auto b = Real{-1}*(i-1)/i;
            x = a*y1 + b*y0 + 2*omega*a*rhs;
#endif
          }
          y0[k] = y1[k];
          y1[k] = _x[k];
        }
      }
#if defined OPT || defined WP
      for (auto k : expansionRange())
        _x[k] = res[k];
#endif

    }


    void iterate(const int iter, const int niter, const Vector &u0, bool verbose)
    {
      const int nstage = static_cast<int>(1+2*std::sqrt(_pde.cfl()));  /* stiffff */
      iterate(u0, nstage);
      if (verbose)
      {
        printf(std::cerr, " nstage= % \n", nstage);
      }
    }

    void solve_system(const Vector& u0)
    {
      using std::get;
      size_t  niter = 5; //8*2*2; // * 32; //*2; //16 ;//1; //32; //50;
      niter = 31;
      constexpr Real tol = 1.0e-7;
      constexpr Real atol = tol;
      constexpr Real rtol = tol;

      bool verbose = _verbose;

      for (auto iter : range_iterator{0,niter})
      {
        iterate(iter,niter, u0, verbose);
        verbose = false;

        auto err = Real{0};
        for (auto i : range_iterator{0,u0.size()})
        {
          const auto y0 = _y0[i];
          auto y1 = Real{0};
          for (auto k : expansionRange())
            y1 += Expansion::oneVec(k)*_x[k][i];
          _y0[i] = y1;

          const auto aerr = std::abs(y1-y0);
          const auto ym  = std::max(std::abs(u0[i]+y0), std::abs(u0[i]+y1));
          err += square(aerr/(atol + rtol*ym));
        }
        err = std::sqrt(err/(u0.size()));

        if (_verbose)
        {
          printf(std::cerr, " >>  iter= %  err= % \n", iter, err);
          if (err < 1)
            break;
        }
        if (iter == niter - 1 && _verbose)
          printf(std::cerr, "   ** iter= %  err= % \n ", iter, err);
      }
    }

    void update(const bool verbose = true)
    {
      _verbose = verbose;
      using std::get;
      auto du = _pde.state();
      static auto u0 = _pde.state();
      _y0.resize(du.size());

      std::fill(_y0.begin(), _y0.end(), 0);


      static bool firstRun = true;
#if 0
      if (!firstRun)
      {
        auto xtmp = _x;
        const auto& u1 = _pde.state();
        for (auto i : range_iterator{0,u0.size()})
        {
          const auto du = u0[i] - u1[i];
          for (auto k : expansionRange())
          {
            _x[k][i] = du;
            for (auto l : expansionRange())
              _x[k][i] += Expansion::nodeMatrix(l,k)*xtmp[l][i];
          }
        }
      }
      else
#endif
        for (auto k : expansionRange())
          for (auto& x : _x[k])
            x = 0;

      firstRun = false;
        
      u0 = _pde.state();

      const auto cfl0 = _pde.get_cfl();
      _pde.set_cfl(1.0*cfl0);
      solve_system(_pde.state());
      auto x_coarse = _x;

      for (auto i : range_iterator{0,u0.size()})
        for (auto k : expansionRange())
        {
          _x[k][i] = 0;
          for (auto l : expansionRange())
            _x[k][i] += Expansion::prolongateMatrix0(k,l)*x_coarse[l][i];
          _y0[i] = 0;
          for (auto k : expansionRange())
            _y0[i] += Expansion::oneVec(k)*_x[k][i];
        }

      _pde.set_cfl(0.5*cfl0);
      solve_system(_pde.state());
      auto x_fine = _x;

      {
        auto du_err = u0;
        std::fill(du_err.begin(), du_err.end(), 0);
        for (auto k : expansionRange())
        {
          const auto wone = Expansion::oneVec(k);
          const auto wmid = Expansion::midVec(k);
          for (auto v : make_zip_iterator(du_err, x_coarse[k],x_fine[k]))
            get<0>(v) += wmid*get<1>(v) - wone*get<2>(v);
        }

        auto err = Real{0};
        for (auto& x : du_err)
          err += square(x);
        err = std::sqrt(err/du_err.size());
        printf(std::cerr, " -- err_half= % \n" ,err);
      }


      auto x = x_coarse;
      _pde.set_cfl(cfl0);
      for (auto k : expansionRange())
      {
        for (auto v : make_zip_iterator(x[k], _pde.state()))
          get<0>(v) += get<1>(v);
        _pde.compute_rhs(_rhs[k], x[k]);
      }

      std::fill(du.begin(), du.end(), 0);
      for (auto k : expansionRange())
      {
        const auto w = Expansion::weight(k);
        for (auto v : make_zip_iterator(du,_rhs[k]))
        {
          get<0>(v) += w*get<1>(v);
        }
      }

#if 0
      for (auto u : du)
        if (std::abs(u) > 1.0e-10)
          printf(std::cerr,  "u= % ", u);
#endif
      _pde.update(du);
      _time += _pde.dt();
    }

#if 0
    void integrate(Real dt)
    {
    }
#endif

};

template<typename real_type>
class PDEDiffusion
{
  public:
    using Real   = real_type;
    using Vector = std::vector<Real>;

  private:
    using range_iterator = make_range_iterator<size_t>;

    Vector _f;

    Real _cfl;
    Real _dx;
    Real _diff;

    size_t n_rhs_calls;

  public:
    void set_dx(const Real dx) { _dx = dx;}
    void set_diff(const Real diff) { _diff = diff;}
    void set_cfl(const Real cfl) { _cfl = cfl;}
    auto get_cfl() const { return _cfl; }
    Real dt() const {
      return  _cfl * 0.5*square(_dx)/_diff;
    }

    auto cost() const { return n_rhs_calls; }
    Real AbsEV() const
    {
      return dt() * 4.0*_diff/square(_dx);  /* 2.0 * cfl */
    }

    auto dx() const { return _dx; }

    PDEDiffusion(const size_t n) : _f{Vector(n+2)}, n_rhs_calls{0}
    {
    }

    static void periodic_bc(Vector &f) 
    {
      const auto n = f.size();
      f[0  ] = f[n-2];
      f[n-1] = f[1  ];
    }

    static void free_bc(Vector &f) 
    {
      const auto n = f.size();
      f[0  ] = f[  1];
      f[n-1] = f[n-2];
    }

    auto cfl() const { return _cfl; }
    auto resolution() const { return _f.size(); }

    void apply_bc(Vector &f) const
    {
      periodic_bc(f);
      //free_bc(f);
    }

    const Vector& state() const { return _f; }

    void update(const Vector &df) 
    {
      using std::get;
      for (auto v: make_zip_iterator(_f,df))
      {
        get<0>(v) += get<1>(v);
      }
    }

    template<typename Func>
      void compute_rhs(Vector &res, Vector &x, Func func)
      {
        apply_bc(x);
        n_rhs_calls++;
        const auto c = dt() * _diff/square(_dx);
        for (auto i : range_iterator{1,x.size() - 1})
        {
          res[i] = c * (x[i+1] - Real{2.0} * x[i] + x[i-1]);
          res[i] = func(res[i]);
        }
        apply_bc(res);
      }
    void compute_rhs(Vector &res, Vector &x)
    {
      compute_rhs(res, x, [](const auto x) { return x; });
    }

    void set_ic()
    {
      using std::max;

      /* set  a profile of delta function */

      auto &f = _f;

      const int n = f.size();
      const auto dx = _dx;
      const auto L = dx*(n-2);
      const auto dL = L * 0.1;
      const auto ic = n>>1;

      const auto ampl = Real{10.0};
      const auto slope = ampl/dL;


      const auto fmin = Real{1};
      std::fill(f.begin(), f.end(), fmin);
      const int m = static_cast<int>(dL/dx + 0.5);
      for (int i = -m; i <= m; i++)
      {
        const auto x = L/2 + dx*i;
        f[ic - ic/2 + i] = std::max(ampl - slope*(std::abs(L/2-x)),fmin);
      }
      for (int i = -m*2; i <= m*2; i++)
      {
        const auto x = L/2 + dx*i;
        f[ic + ic/2 + i] = std::max(1.5*ampl - slope*(std::abs(L/2-x)),fmin);
      }
    }
};

template<typename Solver>
void dump2file(const Solver &solver, std::string fileName = std::string{})
{
  const auto &f = solver.pde().state();
  const auto n = f.size();
  const auto dx = solver.pde().dx();
  std::ofstream foutFn(fileName);
  std::ostream &fout = fileName.empty() ? std::cout : foutFn;
  fout << "# time= " << solver.time() << std::endl;
  fout << "# n= " << n-2 << std::endl;
  for (size_t i = 1; i < n-1; i++)
  {
    const auto x = (i-1)*dx;
    fout << x << " " << std::setprecision(16) << f[i] << std::endl;
  }
}

  template<typename Solver>
auto compute_mass(const Solver &solver)
{
  const auto &f = solver.pde().state();
  const auto n = f.size();
  const auto dx = solver.pde().dx();
  typename Solver::Real sum = 0;
  for (int i = 1; i < n-1; i++)
    sum += f[i]*dx;
  return sum;
}

int main(int argc, char * argv[])
{
  using Real = double;

  const size_t ncell = argc > 1 ? atoi(argv[1]) : 128;
  printf(std::cerr, "ncell= %\n", ncell);

  const Real tau = argc > 2 ? atof(argv[2]) : 0.005;
  printf(std::cerr, "tau= %\n", tau);
  


  constexpr auto ORDER = 7;
  using PDE = PDEDiffusion<Real>;
  using Solver = ODESolverT<ORDER,PDE>;


  Solver solver(PDE{ncell});

  solver.pde().set_dx(1.0/ncell);
  solver.pde().set_diff(1);
  solver.pde().set_cfl(0.8*64); //*64); //*64); //*64); //*64);//*64); //*64);//*64); //*4); //*64/4); //*64); //*64); //*64/4); //*64*4);//*64); //*64); //*64); //*4*4*4);  /* stable for cfl <= 0.5 */

  const auto dt = solver.pde().dt();
  const size_t nstep = 1 + std::max(size_t{0}, static_cast<size_t>(tau/dt));

  printf(std::cerr, "dt= %  tau= %   nstep= %\n",
      dt, tau, nstep);

  solver.pde().set_ic();
  const auto mass0 = compute_mass(solver);
  dump2file(solver, "ic.txt");
  for (size_t step = 1; step <= nstep; step++)
  {
    auto verbose_step = (step-1)%1 == 0;
    auto verbose_iter = (step-1)%1  == 0;
    if (step == nstep || step == 1)
      verbose_step = verbose_iter = true;
    solver.update(verbose_iter);
    if (verbose_step)
    {
      const auto mass = compute_mass(solver);
      printf(std::cerr, "step= % : time= % ORDER= % cost= % -- mass_err= %\n", step, solver.time(), ORDER,
          solver.pde().cost(), (mass-mass0)/mass0);
    }
  }
  printf(std::cerr, " Writing output ... \n");
  dump2file(solver);
  printf(std::cerr, "cost= %\n", solver.pde().cost());



  return 0;  

}




