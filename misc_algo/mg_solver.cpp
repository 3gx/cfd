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


    static constexpr auto matrix(const size_t i, const size_t j);
    static constexpr auto maxAbsMEV();
    static constexpr auto preconditioner(const size_t i, const size_t j);
    static constexpr auto maxAbsPEV();

    static constexpr auto weight(const size_t i);
    static constexpr auto weight_half(const size_t i);
    static constexpr auto node(const size_t i);
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j);
};

template<size_t N, typename T, typename Real>
class ExpansionT;

template<typename T, typename Real>
class ExpansionT<1,T,Real> : public ExpansionBaseT<1,T,Real>
{
  protected:
    static constexpr size_t N = 1;
    using base_type   = ExpansionBaseT<1,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    {
      return Real{2};
    }
    static constexpr auto maxAbsMEV() 
    {
      return Real{2};
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    { 
      return Real{0.5};
    }
    static constexpr auto maxAbsPEV() 
    { 
      return Real{0.5};
    }

    static constexpr auto weight(const size_t i) 
    {
      return Real{1};
    }
    static constexpr auto weight_half(const size_t i)
    {
      assert(0);
      return Real{1};
    }
    static constexpr auto node(const size_t i) 
    {
      return Real{0.5};
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(2 == ORDER || 3 == ORDER, " Pronlogation order is not supported");
        if (N == ORDER)
        {
          return static_cast<Real>(i==j);
        }
        else if (2 == ORDER)
        {
          constexpr Real matrix[2][N]={{1},{1}};
          return matrix[i][j];
        }
        else if (3 == ORDER)
        {
          constexpr Real matrix[3][N]={{1},{1},{1}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  protected:
    static constexpr size_t N = 2;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{3.0,  0.4641016151377546},{-6.464101615137754,  3.0}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV = 3.5;
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = { {0.250000000000000000000,  0},{0.53867513459481288225, 0.250000000000000000000} };
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV = 0.25;
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] = {0.5,0.5};
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)
    {
      assert(0);
      return Real{1};
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.211324865405187118,0.788675134594812882};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(3 == ORDER || 4 == ORDER, " Pronlogation order is not supported");
        if (N == ORDER)
        {
          return static_cast<Real>(i==j);
        }
        else if (3 == ORDER)
        {
          constexpr Real matrix[3][N]={{1.170820393249936908923,-0.170820393249936908923},{0.500000000000000000000,0.500000000000000000000},{-0.170820393249936908923,1.170820393249936908923}};
          return matrix[i][j];
        }
        else if (4 == ORDER)
        {
          constexpr Real matrix[4][N]={{1.245765921961681556807,-0.245765921961681556807},{0.794432220549629981178,0.205567779450370018822},{0.205567779450370018822,0.794432220549629981178},{-0.245765921961681556807,1.245765921961681556807}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  protected:
    static constexpr size_t N = 3;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{5.0,  1.1639777949432226, -0.1639777949432225},{-5.727486121839514,  2.0, 0.7274861218395141},{10.163977794943223,  -9.163977794943223, 5.0}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{5.1};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = { {0.13888888888888888889, 0, 0}, {0.30026319498086459244, 0.22222222222222222222,  0}, {0.26798833376246945173, 0.48042111196938334790,  0.13888888888888888889}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.22222};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] = {0.2777777777777778,0.4444444444444444,0.2777777777777778};
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)
    {
      assert(0);
      return Real{1};
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.112701665379258311,0.500000000000000000,0.887298334620741689};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(4 == ORDER || 5 == ORDER, " Pronlogation order is not supported");
        if (N == ORDER)
        {
          return static_cast<Real>(i==j);
        }
        else if (4 == ORDER)
        {
          constexpr Real matrix[4][N]={{1.173824221557882097734,-0.235926245243015346149,0.062102023685133248415},{0.315779411635934322717,0.807354816671586774721,-0.123134228307521097438},{-0.123134228307521097438,0.807354816671586774721,0.315779411635934322717},{0.062102023685133248415,-0.235926245243015346149,1.173824221557882097734}};
          return matrix[i][j];
        }
        else if (5 == ORDER)
        {
          constexpr Real matrix[5][N]={{1.269238169652725404510,-0.368603188642368014803,0.099365018989642610293},{0.589204776685259874958,0.516751336790516162951,-0.105956113475776037910},{0,1.000000000000000000000,0},{-0.105956113475776037910,0.516751336790516162951,0.589204776685259874958},{0.099365018989642610293,-0.368603188642368014803,1.269238169652725404510}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<4,T,Real> : public ExpansionBaseT<4,T,Real>
{
  protected:
    static constexpr size_t N = 4;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{7.73861278752583057,2.045089650303908785,-0.4370708023957989036,0.0866440235032616747},{-7.201340999706890565,2.261387212474169433,1.448782034533681252,-0.2331339812124165654},{6.343622218624971333,-5.971556459482020118,2.261387212474169433,1.090852762294335798},{-15.56386959855492281,11.89278387805684136,-13.50080272596495125,7.73861278752583057}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{6.8};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = {{0.0869637112843634643,0,0,0},{0.1881181174998680717,0.1630362887156365357,0,0},{0.1671919219741887732,0.3539530060337439665,0.1630362887156365357,0},{0.1774825722545226118,0.3134451147418683468,0.3526767575162718646,0.0869637112843634643}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.17};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] = {0.17392742256872748,0.32607257743127305,0.3260725774312732,0.17392742256872748};
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)  
    { 
      constexpr Real weight[] =
      {0.1761059379386620,0.3049159394838622,0.02115663794741092,-0.002178515369935093};
      return weight[i];
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.0694318442029737124,0.330009478207571868,0.669990521792428132,0.930568155797026288};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(N==ORDER||5==ORDER || 6==ORDER, " Pronlogation order is not supported");
        if (N == ORDER)
        {
          constexpr Real matrix[N][N] = 
          {
            {1.24706002059021446916,-0.37134646293290967052,0.172612312282314955982,-0.048325869939619754622},{0.473360762073713481246,0.69446086758875954880,-0.226915961206226974307,0.059094331543753944264},{-0.007381463003215998440,0.99585142579386540219,0.014821400117421237263,-0.0032913629080706410102},{-0.095609594648018673259,0.70865608285499706932,0.468295800089416859158,-0.081342288296395255222}
          };
          return matrix[i][j];
          return static_cast<Real>(i==j);
        }
        else if (5==ORDER)
        {
          constexpr Real matrix[5][N]={{1.15665233447960049451,-0.23306848457062838863,0.105895713739470568226,-0.029479563648442674104},{0.226361866637914228271,0.93205208186650974176,-0.210599724002141635027,0.052185775497717664992},{-0.092326598440728820911,0.592326598440728820911,0.592326598440728820911,-0.092326598440728820911},{0.052185775497717664992,-0.210599724002141635027,0.93205208186650974176,0.226361866637914228271},{-0.029479563648442674104,0.105895713739470568226,-0.23306848457062838863,1.15665233447960049451}};
          return matrix[i][j];
        }
        else if (6==ORDER)
        {
          constexpr Real matrix[6][N]={{1.25427669651204691779,-0.38249201436315504170,0.178098950270956183941,-0.049883632419848060027},{0.454139585055516588660,0.71591916408312051758,-0.229700084831735121197,0.059641335693098014958},{-0.059826670151450920280,0.93065512687185849208,0.163036458535962670788,-0.033864915256370242588},{-0.033864915256370242588,0.163036458535962670788,0.93065512687185849208,-0.059826670151450920280},{0.059641335693098014958,-0.229700084831735121197,0.71591916408312051758,0.454139585055516588660},{-0.049883632419848060027,0.178098950270956183941,-0.38249201436315504170,1.25427669651204691779}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<5,T,Real> : public ExpansionBaseT<5,T,Real>
{
  protected:
    static constexpr size_t N = 5;
    using base_type = ExpansionBaseT<N,T,Real>;

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
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{8.6};
      return maxAbsMEV;
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
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.15};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] =
      {
        0.11846344252809503,
        0.2393143352496833,
        0.28444444444444444,
        0.23931433524968349,
        0.1184634425280951
      };
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)
    {
      assert(0);
      return Real{1};
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.0469100770306680036,0.230765344947158454,0.500000000000000000,0.769234655052841546,0.953089922969331996};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(6==ORDER || 7==ORDER, " Pronlogation order is not supported");
        if (N == ORDER)
        {
          return static_cast<Real>(i==j);
        }
        else if (6==ORDER)
        {
          constexpr Real matrix[6][N]={{1.13893085681645019777,-0.21522509329330665173,0.11765878765365128281,-0.057649393219534367543,0.016284842042739538701},{0.174920193735091690440,0.98872192203863736418,-0.237460403881069401988,0.101156921896816548177,-0.027338633789476200814},{-0.072959537243763693175,0.46001867145844468644,0.74790078978113712745,-0.177504436580858772501,0.042544512585040651789},{0.042544512585040651789,-0.177504436580858772501,0.74790078978113712745,0.46001867145844468644,-0.072959537243763693175},{-0.027338633789476200814,0.101156921896816548177,-0.237460403881069401988,0.98872192203863736418,0.174920193735091690440},{0.016284842042739538701,-0.057649393219534367543,0.11765878765365128281,-0.21522509329330665173,1.13893085681645019777}};
          return matrix[i][j];
        }
        else if (7==ORDER)
        {
          constexpr Real matrix[7][N]={{1.23293139187926606368,-0.36502898488802090147,0.20433426000800412366,-0.100764511494549110558,0.028527844495299824674},{0.364017066156523892597,0.83590791919108568654,-0.296160457335118513964,0.132610128183963825203,-0.036374656196454890377},{-0.076438912307116058844,0.81669399401294315353,0.345295804642875552215,-0.114700513771969879026,0.0291496274232672321267},{0,0,1.00000000000000000000,0,0},{0.0291496274232672321267,-0.114700513771969879026,0.345295804642875552215,0.81669399401294315353,-0.076438912307116058844},{-0.036374656196454890377,0.132610128183963825203,-0.296160457335118513964,0.83590791919108568654,0.364017066156523892597},{0.028527844495299824674,-0.100764511494549110558,0.20433426000800412366,-0.36502898488802090147,1.23293139187926606368}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<6,T,Real> : public ExpansionBaseT<6,T,Real>
{
  protected:
    static constexpr size_t N = 6;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{15.325599438771348,4.4287845932100730,-1.13579253120090091,0.41365582265249371,-0.15373639381799382,0.037475943927231130},{-12.2744915101690641,3.5536467118620921,3.1045943359124462,-0.89624888663104690,0.30843807418378992,-0.073008911144639749},{7.3152536048277402,-7.2146669386387277,2.1207538493665598,2.5760765596341566,-0.69100775654166641,0.151458139198923446},{-7.0508466361970577,5.5120329300063352,-6.8175842583672762,2.1207538493665598,2.3936417651740589,-0.41586510782960595},{10.2442850874922152,-7.4157314979079741,7.1492010470565008,-9.3575464963379001,3.5536467118620921,2.1032153338214886},{-30.688674821469927,21.571605738296697,-19.256962888352936,19.979099596901344,-25.846653937688777,15.325599438771348}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{10.3};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = {{0.042831123094792586,0,0,0,0,0},{0.092673491430378863,0.090190393262034652,0,0,0,0},{0.082247922612843874,0.19603216233324501,0.116978483643172762,0,0,0},{0.087737871974451507,0.17239079462440697,0.25443949503200162,0.116978483643172762,0,0},{0.084306685134100111,0.18526797945210698,0.22359381104609910,0.25425706957958511,0.090190393262034652,0},{0.086475026360849935,0.17752635320896997,0.23962582533582904,0.22463191657986777,0.19514451252126672,0.042831123094792586}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.12};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] = {0.08566224618958558,0.18038078652406941,0.23395696728634574,0.2339569672863455,0.1803807865240694,0.0856622461895856};
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)
    {
      assert(0);
      return Real{1};
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.0337652428984239861,0.169395306766867743,0.380690406958401546,0.619309593041598454,0.830604693233132257,0.966234757101576014};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(7 == ORDER || 8==ORDER , " Pronlogation order is not supported");
        if (N == ORDER)
        {
          return static_cast<Real>(i==j);
        }
        else if (8==ORDER)
        {
          constexpr Real matrix[8][N]={{1.2123234054623749459,-0.3398325104596944063,0.20762458113598606564,-0.12497742976428808870,0.06268104496578014909,-0.017819091340158665590},{0.30114498455171639951,0.90981981467510605944,-0.32557494616630780944,0.175493801603982138991,-0.084535006330491990878,0.023651351665995202375},{-0.078567615249208809356,0.71012814196288668962,0.49505923000084064385,-0.185878090697565371951,0.081187027842311063063,-0.0219286938592642152157},{0.016010396527596768376,-0.075640250077693814488,0.96543610387910536064,0.126233070636759027634,-0.042786073038554527321,0.0107467520727871851641},{0.0107467520727871851641,-0.042786073038554527321,0.126233070636759027634,0.96543610387910536064,-0.075640250077693814488,0.016010396527596768376},{-0.021928693859264215216,0.08118702784231106306,-0.18587809069756537195,0.49505923000084064384,0.71012814196288668962,-0.07856761524920880936},{0.023651351665995202376,-0.08453500633049199088,0.17549380160398213899,-0.32557494616630780944,0.9098198146751060594,0.30114498455171639951},{-0.017819091340158665590,0.062681044965780149092,-0.12497742976428808870,0.20762458113598606563,-0.33983251045969440634,1.21232340546237494591}};
          return matrix[i][j];
        }
        else if (7==ORDER)
        {
          constexpr Real matrix[7][N]={{1.12367272389532395323,-0.19569689243191671704,0.11690444536532644914,-0.06993129195399033189,0.034987419492685603320,-0.009936404367428956767},{0.141920168224500613243,1.01666020765414094417,-0.23937682454980414246,0.122823477856278604453,-0.058214596919792920440,0.0161875677346769010342},{-0.059992991312292368686,0.37283219741336188404,0.83933361670899686975,-0.217790754428853929500,0.089225063511153017508,-0.023607131892365473110},{0.0354336891832800019050,-0.150585800696685203166,0.61515211151340520126,0.61515211151340520126,-0.150585800696685203166,0.0354336891832800019050},{-0.023607131892365473110,0.089225063511153017508,-0.217790754428853929500,0.83933361670899686975,0.37283219741336188404,-0.059992991312292368686},{0.0161875677346769010342,-0.058214596919792920440,0.122823477856278604453,-0.23937682454980414246,1.01666020765414094417,0.141920168224500613243},{-0.009936404367428956767,0.034987419492685603320,-0.06993129195399033189,0.11690444536532644914,-0.19569689243191671704,1.12367272389532395323}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<7,T,Real> : public ExpansionBaseT<7,T,Real>
{
  protected:
    static constexpr size_t N = 7;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{20.162475126453963,5.939649595901784,-1.5713450065769194,0.6117568698401676,-0.2662682556100125,0.10823443910941748,-0.02751051928120109},{-15.629383503722168,4.443146050695284,4.125140303481284,-1.2701515428956040,0.5100682216357629,-0.20014625792437366,0.05010533623170008},{8.625187972409475,-8.605076196977838,2.3943788228507531,3.351834730017060,-1.0413621688017852,0.3736448491420073,-0.09029602574207080},{-7.258551263394934,5.727243007023915,-7.245304978834095,2.0000000000000000,3.062096190825881,-0.8500070065068292,0.18952405088606242},{8.182614584790622,-5.956895021669846,5.830119814503291,-7.930862828666295,2.3943788228507531,3.021826024450000,-0.5328694133609241},{-12.929826694060736,9.086438359314942,-8.131847279980262,8.558125387465308,-11.746919361825783,4.443146050695284,2.749662145893132},{40.35246077218913,-27.93020953980242,24.129196096015008,-23.429578355405607,25.43427284698192,-33.76162469659479,20.162475126453963}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{12.1};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = {{0.03237124154221742,0,0,0,0,0,0},{0.07004354137872608,0.06992634787231917,0,0,0,0,0},{0.06215393578734986,0.15200552205783099,0.09545751262627974,0,0,0,0},{0.06633292861768470,0.13359576922388229,0.20770188076597078,0.10448979591836735,0,0,0},{0.06366576746880693,0.14380627590344877,0.18220246265408429,0.22735483605218653,0.09545751262627974,0,0},{0.06548633327460677,0.13720685187790830,0.19631211718445561,0.19962996905329136,0.20750503183140724,0.06992634787231917,0},{0.06428437356068539,0.14145951478168444,0.18773996647887383,0.21411332539996004,0.18328182138013593,0.15130371302782220,0.03237124154221742}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.11};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] = {0.06474248308443512, 0.1398526957446387, 0.19091502525255918, 0.20897959183673512, 0.1909150252525592, 0.13985269574463888, 0.06474248308443513};
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)
    {
      assert(0);
      return Real{1};
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.0254460438286207377,0.129234407200302780,0.297077424311301417,0.500000000000000000,0.702922575688698583,0.870765592799697220,0.974553956171379262};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(8==ORDER , " Pronlogation order is not supported");
        if (N == ORDER)
        {
          return static_cast<Real>(i==j);
        }
        else if (8==ORDER)
        {
          constexpr Real matrix[8][N]={{1.1109764149740766555,-0.1777969288609955349,0.1116515579714447392,-0.07379620075760972588,0.04531368772400935322,-0.02285470615610564788,0.006506175105180160662},{0.119110236177449917140,1.03107393717281695193,-0.23151265182871470898,0.130013891533531084042,-0.075242583452459995331,0.036957904549612781124,-0.0104007341522360299254},{-0.050802380548726304847,0.31191120843947840758,0.89591714590041723546,-0.23357673335979253852,0.115130395585113660308,-0.053172105475361454533,0.0145924694588709945482},{0.030201104248416676361,-0.129725199899249985583,0.51809795025873322081,0.71911630754113613762,-0.195544510376441069161,0.078272281438027707195,-0.0204179332106226872421},{-0.0204179332106226872421,0.078272281438027707195,-0.195544510376441069161,0.71911630754113613762,0.51809795025873322081,-0.129725199899249985583,0.030201104248416676361},{0.014592469458870994548,-0.05317210547536145453,0.11513039558511366031,-0.23357673335979253852,0.8959171459004172355,0.31191120843947840758,-0.05080238054872630485},{-0.010400734152236029925,0.03695790454961278112,-0.07524258345245999533,0.13001389153353108404,-0.23151265182871470898,1.0310739371728169519,0.11911023617744991714},{0.006506175105180160662,-0.022854706156105647880,0.045313687724009353217,-0.07379620075760972588,0.11165155797144473924,-0.17779692886099553488,1.11097641497407665552}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<8,T,Real> : public ExpansionBaseT<8,T,Real>
{
  protected:
    static constexpr size_t N = 8;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{25.6926112963098,7.66481813729873,-2.066766305200644,0.834967836032650,-0.391322735056362,0.1880860790263934,-0.0807844483958189,0.02109496564232081},{-19.49249456998672,5.47461453891780,5.28490391970636,-1.686816817555106,0.728229725420268,-0.337061781204938,0.1420577647099491,-0.0367733313510094},{10.23940120319752,-10.29567243097482,2.76313437016009,4.22037639163604,-1.405147920784982,0.591814510044655,-0.238937072465555,0.0606926207226390},{-7.93798136827535,6.30583618258298,-8.09856974148034,2.069639794612294,3.76153680354670,-1.215395508135759,0.446517013113076,-0.1092219434734880},{7.81415574865183,-5.71807461470441,5.66351220252970,-7.90081639277129,2.069639794612294,3.65045304708640,-1.034278580991645,0.233047563097002},{-9.63317863204001,6.78823557067877,-6.11808325036483,6.54772240184721,-9.36295087269827,2.76313437016009,3.74637393276161,-0.666915191880153},{16.04021892340721,-11.09128684254554,9.57597223855476,-9.32567356743535,10.28426065957019,-14.52381437705618,5.47461453891780,3.48904897793052},{-51.4063175582620,35.2374993036919,-29.8531646168542,27.9967257720613,-28.4403708730376,31.7318448430285,-42.8215329925948,25.6926112963098}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{14};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = {{0.0253071340725941,0,0,0,0,0,0,0},{0.0547593217675543,0.0555952586133436,0,0,0,0,0,0},{0.0485875359989129,0.1208595249971732,0.0784266614694718,0,0,0,0,0},{0.0518655209705812,0.1061934901483484,0.1706711342745536,0.0906709458445905,0,0,0,0},{0.0497550315609233,0.1143883315370480,0.1496121163777214,0.1973629330102060,0.0906709458445905,0,0,0},{0.0512330738404386,0.1089648024514023,0.1614967888010481,0.1729701589689548,0.1973169950510594,0.0784266614694718,0,0},{0.0501714658408459,0.1127554521376362,0.1537135699534800,0.1865572437783281,0.1731921828308204,0.1704931191747253,0.0555952586133436,0},{0.0508917764723050,0.1102177595626280,0.1587709981935806,0.1782634003208542,0.1858249073022357,0.1505724917919132,0.1202964605326573,0.0253071340725941}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.1};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[N] = {0.050614268145188344,0.11119051722668769,0.15685332293894347,0.18134189168918116,0.18134189168918102,0.15685332293894338,0.11119051722668757,0.05061426814518831};
      return weight[i];
    }
    static constexpr auto weight_half(const size_t i)  
    { 
      constexpr Real weight[] =
      {0.05073647053793456,0.1106371236279321,0.1587799218056337,0.1697302482589022,0.01161164343027877,-0.001926598866690093,0.000553393598755149,-0.0001222023927464280};
      return weight[i];
    }
    static constexpr auto node(const size_t i)  
    { 
      constexpr Real node[N] = {0.0337652428984239861,0.169395306766867743,0.380690406958401546,0.619309593041598454,0.830604693233132257,0.966234757101576014};
      return node[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(N==ORDER , " Pronlogation order is not supported");
        if (N==ORDER)
        {
          constexpr Real matrix[N][N]={{1.2123234054623749459,-0.3398325104596944063,0.20762458113598606564,-0.12497742976428808870,0.06268104496578014909,-0.017819091340158665590},{0.30114498455171639951,0.90981981467510605944,-0.32557494616630780944,0.175493801603982138991,-0.084535006330491990878,0.023651351665995202375},{-0.078567615249208809356,0.71012814196288668962,0.49505923000084064385,-0.185878090697565371951,0.081187027842311063063,-0.0219286938592642152157},{0.016010396527596768376,-0.075640250077693814488,0.96543610387910536064,0.126233070636759027634,-0.042786073038554527321,0.0107467520727871851641},{0.0107467520727871851641,-0.042786073038554527321,0.126233070636759027634,0.96543610387910536064,-0.075640250077693814488,0.016010396527596768376},{-0.021928693859264215216,0.08118702784231106306,-0.18587809069756537195,0.49505923000084064384,0.71012814196288668962,-0.07856761524920880936},{0.023651351665995202376,-0.08453500633049199088,0.17549380160398213899,-0.32557494616630780944,0.9098198146751060594,0.30114498455171639951},{-0.017819091340158665590,0.062681044965780149092,-0.12497742976428808870,0.20762458113598606563,-0.33983251045969440634,1.21232340546237494591}};
          return matrix[i][j];
        }
      }
};

template<size_t ORDER_MAX, typename PDE>
class ODESolverT
{
  public:
    using Real           = typename PDE::Real;
    using Vector         = typename PDE::Vector;
  private:
    using range_iterator = make_range_iterator<size_t>;

    template<size_t ORDER>
      using Expansion = ExpansionT<ORDER,Vector,Real>;

    PDE _pde;
    Real _time;
    bool _verbose;

    Real _err, _err_pre, _cfl_pre, _cfl;
    Real _atol, _rtol;

    typename Expansion<ORDER_MAX>::storage _x, _rhs, _rhs_pde;

    template<size_t K>  
      auto expansionRange() const 
      {
        return make_range_iteratorT<0,K>{};
      }
    
    static constexpr Real omegaCFL = 0.5;

  public:
    auto time() const {return _time;}

    ODESolverT(const PDE &pde) : _pde{pde}, _time{0}
    {
      for (auto k : expansionRange<ORDER_MAX>())
      {
        _x  [k].resize(_pde.resolution());
        _rhs[k].resize(_pde.resolution());
        _rhs_pde[k].resize(_pde.resolution());
        std::fill(_x[k].begin(), _x[k].end(), 0);
        std::fill(_rhs_pde[k].begin(), _rhs_pde[k].end(), 0);
      }
      _err = _err_pre = -1;
      _cfl = _cfl_pre = -1;

      const Real tol = 1.0e-12;
      _atol = tol;
      _rtol = tol;
    };
    PDE& pde() { return _pde; }
    const PDE& pde() const { return _pde; }

    template<size_t ORDER>
      void rhs(const Vector &u0)
      {
        using std::get;
        auto x = _x;

        /* compute RHS */
        const auto h = _pde.dt();
        for (auto k : expansionRange<ORDER>())
        {
          for (auto v : make_zip_iterator(x[k], u0))
            get<0>(v) += get<1>(v);
          _pde.compute_rhs(_rhs_pde[k], x[k]);

          assert(_x[k].size() == u0.size());
          for (auto i : range_iterator{u0.size()})
          {
            Real r = 0;
            for (auto l : expansionRange<ORDER>())
              r += Expansion<ORDER>::matrix(k,l) * _x[l][i];
            _rhs[k][i] = h*_rhs_pde[k][i]  - r;
          }
        }

        /* precondition RHS */
        for (auto i : range_iterator{0,u0.size()})
        {
          std::array<Real,Expansion<ORDER>::size()> tmp;
          for (auto k : expansionRange<ORDER>())
          {
            tmp[k] = 0;
            for (auto l : expansionRange<ORDER>())
              tmp[k] += Expansion<ORDER>::preconditioner(k,l)*_rhs[l][i];
          }
          for (auto k : expansionRange<ORDER>())
            _rhs[k][i] = tmp[k];
        }
      }

    template<size_t ORDER>
      void smoother_OPT(const int n_smooth_iter, const Vector &u0)
      {
        const auto n = n_smooth_iter;

        using std::get;

        const Real omega = omegaCFL/(Expansion<ORDER>::maxAbsPEV()*(Expansion<ORDER>::maxAbsMEV() + _pde.AbsEV()));

        static decltype(_x) y0, y1,tmp;
        y0  = _x;

        rhs<ORDER>(u0);
        static auto res = _x;


        auto frac = Real{1}/(1+n)/(2+n)/(3+4*n);
        auto omega0 = frac*6;
        auto omega1 = frac*3*n*(3+n);
        auto omegak = [n,frac](const int k) 
        {
          return frac*3*(1-k+n)*(2+k+n);
        };

        for (auto k : expansionRange<ORDER>())
        {
          for (auto v : make_zip_iterator(_x[k], _rhs[k],res[k]))
          {
            auto&   x = get<0>(v);
            auto& rhs = get<1>(v);
            auto& r = get<2>(v);
            r = x*omega0 + (3*x + 4.0*omega*rhs)*omega1;
            x = x + 2.0*omega*rhs;
          }
          y1[k] = _x[k];
        }



        for (int i = 2; i <= n; i++)
        {
          rhs<ORDER>(u0);
          for (auto k : expansionRange<ORDER>())
          {
            for (auto v : make_zip_iterator(_x[k], y1[k], y0[k], _rhs[k],res[k]))
            {
              auto&   x = get<0>(v);
              auto&  y1 = get<1>(v);
              auto&  y0 = get<2>(v);
              auto& rhs = get<3>(v);
              x = 2*y1 - y0 + 4*omega*rhs;
              auto& r = get<4>(v);
              r += 2*omegak(i)*x;
            }
            y0[k] = y1[k];
            y1[k] = _x[k];
          }
        }
        for (auto k : expansionRange<ORDER>())
          _x[k] = res[k];
      }

    template<size_t ORDER>
      void smoother_WP(const int n_smooth_iter, const Vector &u0)
      {
        const auto n = n_smooth_iter;
        using std::get;

        const Real omega = omegaCFL/(Expansion<ORDER>::maxAbsPEV()*(Expansion<ORDER>::maxAbsMEV() + _pde.AbsEV()));

        static decltype(_x) y0, y1,tmp;
        y0  = _x;

        rhs<ORDER>(u0);
        static auto res = _x;

        auto scale = Real{1}/(2*n+1);

        for (auto k : expansionRange<ORDER>())
        {
          for (auto v : make_zip_iterator(_x[k], _rhs[k],res[k]))
          {
            auto&   x = get<0>(v);
            auto& rhs = get<1>(v);
            auto& r = get<2>(v);
            r = (3*x + 4.0*omega*rhs)*scale;
            x = x + 2.0*omega*rhs;
          }
          y1[k] = _x[k];
        }



        for (int i = 2; i <= n; i++)
        {
          rhs<ORDER>(u0);
          for (auto k : expansionRange<ORDER>())
          {
            for (auto v : make_zip_iterator(_x[k], y1[k], y0[k], _rhs[k],res[k]))
            {
              auto&   x = get<0>(v);
              auto&  y1 = get<1>(v);
              auto&  y0 = get<2>(v);
              auto& rhs = get<3>(v);
              x = 2*y1 - y0 + 4*omega*rhs;
              auto& r = get<4>(v);
              r += 2*scale*x;
            }
            y0[k] = y1[k];
            y1[k] = _x[k];
          }
        }
        for (auto k : expansionRange<ORDER>())
          _x[k] = res[k];
      }


    template<size_t ORDER_OLD,size_t ORDER_NEW, bool TYPE = true>
      void solve_system_mg(const size_t n_smooth_iter, const Vector &u0)
      {
        if (_verbose)
          printf(std::cerr, " ------ % -> % \n", ORDER_OLD, ORDER_NEW);
        const auto n = u0.size();
        const auto h = _pde.dt();

        Vector y0;
        y0.resize(n);

        using ExpansionNew = Expansion<ORDER_NEW>;
        using ExpansionOld = Expansion<ORDER_OLD>;

        /* interpolate from old to new order */
        for (auto i : range_iterator{n})
        {
          std::array<Real,ExpansionNew::size()> x;
          for (auto k : expansionRange<ExpansionNew::size()>())
          {
            x[k] = 0;
            for (auto l : expansionRange<ExpansionOld::size()>())
            {
              const auto weight = ExpansionOld::template prolongate<ExpansionNew::size()>(k,l);
              x[k] += weight * _x[l][i];
            }
          }
          for (auto k : expansionRange<ExpansionNew::size()>())
            _x[k][i] = x[k];

          y0[i] = 0;
          if (ORDER_OLD != ORDER_NEW)
          {
            for (auto l : expansionRange<ExpansionOld::size()>())
              y0[i] += ExpansionOld::weight(l)*h*_rhs_pde[l][i];
          }
          else
          {
            for (auto k : expansionRange<ExpansionNew::size()>())
              y0[i] += ExpansionNew::weight_half(k)*2.0*h*_rhs_pde[k][i];
          }

        }

        size_t n_iter = 16;

        const Real atol = _atol;
        const Real rtol = _rtol;

        Real err = 0;
        for (auto iter : range_iterator{n_iter})
        {
          if (TYPE)
            smoother_OPT<ORDER_NEW>(n_smooth_iter, u0);
          else
            smoother_WP<ORDER_NEW>(n_smooth_iter, u0);

          {
            using std::get;

            /* compute RHS */
            auto xtmp = _x;
            for (auto k : expansionRange<ExpansionNew::size()>())
            {
              for (auto v : make_zip_iterator(xtmp[k], u0))
                get<0>(v) += get<1>(v);
              _pde.compute_rhs(_rhs_pde[k], xtmp[k]);
            }
          }
          for (auto i : range_iterator{n})
          {
            auto du0 = y0[i];
            decltype(du0) du1 = 0;
            for (auto k : expansionRange<ExpansionNew::size()>())
              du1 += ExpansionNew::weight(k)*h*_rhs_pde[k][i];
            y0[i] = du1;

            const auto aerr = std::abs(du1-du0);
            const auto ym  = std::max(std::abs(u0[i]+du0), std::abs(u0[i]+du1));
            err += square(aerr/(atol + rtol*ym));
          }
          err = std::sqrt(err/n);

          if (_verbose)
          {
            printf(std::cerr, " >>  iter= %  err= % \n", iter, err);
          }
          if (err < 1)
            break;

          if (n_iter-1 == iter && _verbose)
            printf(std::cerr, "   ** iter= %  err= % \n ", iter, err);
        }
      }

    void update(const bool verbose = true)
    {
      _verbose = verbose;
      using std::get;
      const auto u0 = _pde.state();
      const auto n = u0.size();

      const auto h = _pde.dt();

      auto n_smooth_iter = static_cast<size_t>(1+4*std::sqrt(_pde.cfl()));
      if (_verbose)
        printf(std::cerr, " n_smooth_iter= % \n", n_smooth_iter);

      Vector du_ctrl, du_solv;
      du_ctrl.resize(n);  
      du_solv.resize(n);

      static_assert(4 == ORDER_MAX || 8 == ORDER_MAX, " Order mismatch ");
      if (4 == ORDER_MAX)
      {
        solve_system_mg<1,2>(n_smooth_iter, u0);
        solve_system_mg<2,3>(n_smooth_iter, u0);
        for (auto i : range_iterator{n})
        {
          du_ctrl[i] = 0;
          for (auto k : expansionRange<ORDER_MAX-1>())
            du_ctrl[i] += Expansion<ORDER_MAX-1>::weight(k)*h*_rhs_pde[k][i];
        }
        solve_system_mg<3,4>(n_smooth_iter, u0);
      }
      else if (8 == ORDER_MAX)
      {
        solve_system_mg<1,2>(n_smooth_iter, u0);
        solve_system_mg<2,3>(n_smooth_iter, u0);
        solve_system_mg<3,5>(n_smooth_iter, u0);
        solve_system_mg<5,6>(n_smooth_iter, u0);
        solve_system_mg<6,7>(n_smooth_iter, u0);
        for (auto i : range_iterator{n})
        {
          du_ctrl[i] = 0;
          for (auto k : expansionRange<ORDER_MAX-1>())
            du_ctrl[i] += Expansion<ORDER_MAX-1>::weight(k)*h*_rhs_pde[k][i];
        }
        solve_system_mg<7,8>(n_smooth_iter, u0);
      }

      for (auto i : range_iterator{n})
      {
        du_solv[i] = 0;
        for (auto k : expansionRange<ORDER_MAX>())
          du_solv[i] += Expansion<ORDER_MAX>::weight(k)*h*_rhs_pde[k][i];
      }
      auto du0 = du_solv;

      if (1)
      {
        for (auto i : range_iterator{n})
        {
          du_solv[i] = 0;
          for (auto k : expansionRange<ORDER_MAX>())
            du_solv[i] += Expansion<ORDER_MAX>::weight_half(k)*h*_rhs_pde[k][i];
        }

        const auto cfl0 = _pde.get_cfl();
        _pde.set_cfl(0.5*cfl0);
        solve_system_mg<ORDER_MAX,ORDER_MAX>(n_smooth_iter,u0);
        const auto h = _pde.dt();
        for (auto i : range_iterator{n})
        {
          du_ctrl[i] = 0;
          for (auto k : expansionRange<ORDER_MAX>())
            du_ctrl[i] += Expansion<ORDER_MAX>::weight(k)*h*_rhs_pde[k][i];
        }

        _pde.set_cfl(cfl0);

      }


      Real err = 0;
      for (auto i : range_iterator{n})
      {
        const auto u_solv = u0[i] + du_solv[i];
        const auto u_ctrl = u0[i] + du_ctrl[i]; 

        const auto um = std::max(std::abs(u0[i]), std::abs(u_solv));
        const auto sc1 = _atol + _rtol*um;
        const auto du_err = std::abs(u_solv - u_ctrl);
        err += square(du_err/sc1);
      }
      err = std::sqrt(err/n);
      _err_pre  = _err;
      _err      = err;
      _cfl_pre  = _cfl;
      _cfl      = _pde.get_cfl();

      _pde.update(du0);
      _time += _pde.dt();

      if (_verbose)
        printf(std::cerr, "# -- err_pre= %   err= %  cfl_pre= %  cfl= %\n", _err_pre, _err,
            _cfl_pre, _cfl);

      Real cfl_scale = 1;
      if (_err > 0 && _err_pre > 0)
      {
        const auto p = Real{1}/(Expansion<ORDER_MAX>::size());
        cfl_scale = 0.8*std::pow(1/_err,p)*_cfl/_cfl_pre*std::pow(_err_pre/_err,p);
      }
      else if (_err > 1)
      {
        const auto p = Real{1}/(Expansion<ORDER_MAX>::size());
        cfl_scale = 0.8*std::pow(1/_err,p);
      }
      if (_verbose)
        printf(std::cerr,"cfl_scale= % \n", cfl_scale);

      const auto cfl0 = _pde.get_cfl();
      const auto cfl1 = cfl0*cfl_scale;
      _pde.set_cfl(cfl1);
      
      for (auto i : range_iterator{n})
      {
        constexpr size_t ORDER0 = 1;
        for (auto k : expansionRange<ORDER0>())
        {
          const auto x = cfl1/cfl0 * Expansion<ORDER0>::node(k);
          _x[k][i] = du_solv[i]*x;
        }
        for (auto k : expansionRange<ORDER0>())
          _rhs_pde[k][i] = 0;
      }

    }
};

template<typename real_type>
class PDEBurger
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

    PDEBurger(const size_t n) : _f{Vector(n+2)}, n_rhs_calls{0}
    {
    }
    auto cfl() const { return _cfl; }
    auto resolution() const { return _f.size(); }


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
        n_rhs_calls++;
//        const auto c = dt() * _diff/square(_dx);
        const auto c =  _diff/square(_dx);
        const auto n = x.size();
        res[0] = c * (x[n-1] - Real{2.0} * x[0] + x[1]);
        res[0] = func(res[0]);
        for (auto i : range_iterator{1,x.size()-1})
        {
          res[i] = c * (x[i-1] - Real{2.0} * x[i] + x[i+1]);
          res[i] = func(res[i]);
        }
        res[n-1] = c * (x[n-2] - Real{2.0} * x[n-1] + x[0]);
        res[n-1] = func(res[n-1]);
      }
    void compute_rhs(Vector &res, Vector &x)
    {
      compute_rhs(res, x, [](const auto x) { return x; });
    }

    void wedge_ic()
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
    void burger_ic()
    {
      using std::max;

      /* set  a profile of delta function */

      auto &f = _f;

      const int n = f.size();
      const auto dx = _dx;
      for (int i = 0; i < n; i++)
      {
        const auto x = (i+1)*dx;
        f[i] = 1.5*x*square(1-x);
      }
    }
    void set_ic()
    {
    //  wedge_ic();
        burger_ic();
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
  for (int i = 0; i < n; i++)
    sum += f[i]*dx;
  return sum;
}

int main(int argc, char * argv[])
{
  using Real = double;

  const size_t ncell = argc > 1 ? atoi(argv[1]) : 128;
  printf(std::cerr, "ncell= %\n", ncell);

  const Real tend = argc > 2 ? atof(argv[2]) : 0.005;
  printf(std::cerr, "tend= %\n", tend);
  


  constexpr auto ORDER = 4;
  using PDE = PDEBurger<Real>;
  using Solver = ODESolverT<ORDER,PDE>;


  Solver solver(PDE{ncell});

  solver.pde().set_dx(1.0/ncell);
  solver.pde().set_diff(1);
  solver.pde().set_cfl(0.8);//*2*2*2); //*2*2); //*8); //*8); //*64); //*16); //*64); //*64); //*64); //*64); //*64);//*64); //*64);//*64); //*4); //*64/4); //*64); //*64); //*64/4); //*64*4);//*64); //*64); //*64); //*4*4*4);  /* stable for cfl <= 0.5 */

  const auto dt = solver.pde().dt();

  solver.pde().set_ic();
  const auto mass0 = compute_mass(solver);
  dump2file(solver, "ic.txt");
  {
    size_t kstep = 0;
    bool keep_looping = true;
    while (keep_looping)
    {
      auto verbose_step = (kstep%1) == 0;
      auto verbose_iter = (kstep%1) == 0;
      const auto dt = solver.pde().dt();
      bool break_loop = true;

      if (solver.time() + dt >= tend)
      {
        const auto dt_new = tend - solver.time();
        assert(dt_new >= 0);
        solver.pde().set_cfl(solver.pde().get_cfl() * dt_new/dt);
        keep_looping = false;
        verbose_step = verbose_iter = true;
      }

      if (verbose_step)
      {
        const auto mass = compute_mass(solver);
        printf(std::cerr, "step= % : time= % dt= % (cfl= %) ORDER= % cost= % -- mass_err= %  Tend= % \n", 
            kstep, solver.time(), solver.pde().dt(), solver.pde().cfl(), ORDER,
            solver.pde().cost(), (mass-mass0)/mass0, tend);
      }
//      if (kstep > 45)
//        solver.pde().set_cfl(1000);
      solver.update(verbose_iter);
      kstep++;
    }
#if 0
    const auto mass = compute_mass(solver);
    printf(std::cerr, "step= % : time= % dt= % (cfl= %) ORDER= % cost= % -- mass_err= %  Tend= % \n", 
        kstep, solver.time(), solver.pde().dt(), solver.pde().cfl(), ORDER,
        solver.pde().cost(), (mass-mass0)/mass0, tend);
#endif
  }
  printf(std::cerr, " Writing output ... \n");
//  dump2file(solver);
  dump2file(solver, "burger");
  printf(std::cerr, "cost= %\n", solver.pde().cost());



  return 0;  

}




