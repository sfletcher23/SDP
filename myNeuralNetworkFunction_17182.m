function [y1] = myNeuralNetworkFunction_17182(x1, adjustOutput)
%MYNEURALNETWORKFUNCTION_17182 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 22-Aug-2017 09:12:11.
% 
% [y1] = myNeuralNetworkFunction_17182(x1) takes these arguments:
%   x = 3xQ matrix, input #1
% and returns:
%   y = 108xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [0.200226052645337;0.0208150562914954;3];
x1_step1.gain = [0.356566960924732;7.3765183675583;0.000182748538011696];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.1716572442209876037;-0.23389926814490932805;-4.0462725780248156582;4.6108710937988526268;-12.50752864908271178];
IW1_1 = [-0.12012742107080752463 -2.3352457994157669319 0.39931130203831211034;-0.17865644756770843182 -0.1031181371162936955 0.1447872270133041539;-3.539650960409249425 -0.28908630665144940952 0.18034220841699188931;2.0284696369260437621 0.22307389432918475758 1.9630503117218152109;-12.240997175850766965 -1.4530250480110442624 1.0862323176457240326];

% Layer 2
b2 = [-0.063927238555105375672;0.38019848817577273081;0.11879378659911742078;-0.091309945997094091319;-0.34095298209633978503;-0.30140167780473292591;-0.95938753498684326004;-0.89308458649113864336;-0.88674781329709473354;0.13155366257879880232;-0.41232218113161661632;-0.45480408685609702291;-0.43811532500611261032;-0.48959890034531478031;-0.45237173870240759976;-0.37737673755980372103;-0.39827677290266466636;-0.4602271965840967205;-0.53013623017679956106;-0.53409995548759681139;-0.45804328824235207573;-0.48477680467000799913;-0.56971185623050135671;-0.60833686804021325223;-0.48419589200313012389;-0.46407897434047390428;-0.68978039385118050664;-0.31282269817656632727;-0.71613888248304224948;-0.56038490281514863334;-0.75243685372827551205;-0.99076779897569067224;-0.15091657144910156085;0.038065580703287228703;-0.27772703730931647126;0.36749065925651408504;-0.052881005505897846208;-0.25993665885985828279;-0.82100279138224696585;-0.49906584918956481411;-0.4058838293606493175;-0.35462920521717766764;-0.53249116295433884805;-0.54129797918893829412;-0.51879099810563245043;-0.70595858194365679328;-0.33324624770008731689;-0.35948394225289348469;-0.53096139774033368486;-0.52734352137396278692;-0.52021291146573467312;-0.43009204694578601069;-0.61685584481979716109;-0.66858558482849073101;-0.70135654564818972201;-0.74109088320064553379;-0.088889333135602183189;-0.48942541990797194806;-0.61358623297792114215;-0.66816204965991776721;-0.60805650268812549974;-0.83771919145231699133;-0.52424999710559982979;0.05640749658910716724;-0.68059308781906513452;-0.72317917428689304504;-0.77188579903117082726;-0.44318480530049592092;-0.8863676021570420227;-0.82690703723937386549;-0.68889860317802376066;-0.70313854685790433408;-0.48553777750580356942;-0.62268928874709061905;-0.15006579040164802485;-0.1721286880288361576;0.87965463041269587041;-0.0018914483607551770342;0.18895463684756941269;-0.36218601011717649341;-0.20674022416911588818;-0.66729532352147780117;-0.39677914391597901567;-0.8033230713342111029;-0.76380892772299402438;-0.88190050494562532979;-0.66913787873441421539;-0.68573829291876786662;-0.82939472456986595983;-0.98342571939510503576;-0.63819492299674551727;-0.29278662440267883005;-0.13898060923818092438;-0.94492699718200490189;-0.94432929027835799118;-0.25317998676438918659;-0.55880450062390174537;-0.68613434888133029155;-0.84710126097416682533;-0.043478298663850382944;-1.3619422277157460233;0.39986505076340450593;-0.22791835675514726667;-0.012013776398931169664;-0.91078686625863147786;-0.72345928386579871905;-0.26243946889287195212;-1.0722979100151901477];
LW2_1 = [-0.56553582827871007321 0.50096050475824915882 -0.082443304007463996674 -0.58397579388852693949 -0.77376896903824510243;-1.034611265427556992 -0.73949563976638865714 -0.16290748445231895136 -1.0739633265387436012 -0.11107363916003085147;0.42241748872102774515 0.42044016129681166527 -0.087267006437210969949 -0.30632784317676980956 -0.1760960730543572883;0.55457377412538533257 0.12778680952598239773 -0.17644523953861959709 -0.35483654188262764295 -0.43322961939100090856;0.85723050733895667097 0.020533004596643551071 -0.87918268566481272419 0.37763546643197543062 -0.079885362915251820204;1.0154508612903920817 -0.21130070646130327816 -0.21565203115339157125 0.58201759575407718827 0.38415828107744731668;0.0092395362096292649179 0.59506494950659294396 -0.29885327486138041353 -0.19641797523615084353 -0.16335091416488300986;0.71791589622879736776 0.54058323728709989364 0.032705615860695744224 0.75657317946315694179 -0.13423656621248231113;-0.075712777697939881438 -0.52096072245037683413 0.7257809495466998051 -0.10249050523724934081 -0.5981208747603388165;1.1361105972769292016 0.10249738523222677589 -0.43195984166192363807 -0.12735644677823781956 -0.43557143843713447273;0.13181642119251735235 -0.40701136367466189681 -0.62592760899091004489 -0.26519702960885266352 -0.4228598013117589316;-0.69692290680967483585 -0.73299078631658631799 -1.5389202962018020404 -1.2319979975736432554 -0.040635652237116785812;-1.07913059916587728 -0.95083586782526230419 -1.0105593754558681585 -1.1929799898691111348 -0.10623236282129064345;-1.1013015662128338157 -0.87190260043256717637 -1.105810669815436631 -1.2150529476233937753 -0.073674082423767106587;-0.91696809292814629444 -0.74240036637223938421 -1.3800851559783406941 -1.325689015272013993 -0.065090540480647093857;-1.0847182741962055896 -0.93932700242565159154 -0.99555629815655677373 -1.2393979516912119188 -0.11245061809994416624;-1.2221483926882257265 -1.0561084897946866334 -0.70644970346713575893 -1.0932579648663267591 -0.12982839414239608078;-1.581919655018895865 -1.2044628123730287683 -0.010251900680432397481 -0.64909011898514357153 -0.11537663849456594733;-1.3614892094196531236 -0.91771825879945756466 -0.7023811991654631548 -0.99537633722237850442 -0.072730513077826450252;-1.1471921961832372006 -0.84445450084997508711 -1.087736458310297305 -1.1781167359819151308 -0.064540825093058418194;-1.2173350052374822816 -0.97147505258175637088 -0.82252797234064700671 -1.1146123916847896496 -0.11326955132315799646;-1.2490707056491348048 -0.91270556869426211488 -0.81864225147302205343 -1.0776618841187650677 -0.095449579620198121699;-1.2765953886750938118 -0.85260649564143775603 -0.83471024049632125408 -0.99188981873264414624 -0.066369774634961797855;-1.1532698422427694052 -0.76694469631476924754 -1.1880679546912769062 -1.1736497159821830394 -0.042174781994540125329;-1.087780154092214735 -0.84692795252672314898 -1.1512977010738210115 -1.2467920513935704108 -0.069957315437235587097;-1.0682690859325454369 -0.86587370509630756654 -1.2189726200929928979 -1.3288624361655925199 -0.076037748744149744051;-1.2269068300394485593 -0.69795289993779130455 -1.2270415836720778291 -1.2190373735187800541 -0.061370353916564007735;-0.15007560495190885863 -0.84134931324873674097 -1.34140682989912019 -0.7269295913604129078 -0.11066090689932250868;-1.2597261956046590381 -0.98972275208731641527 -1.1112492105113878083 -1.2258518365968142927 -0.091539222422821586722;-0.78566002741151619304 -0.64328307919316451979 -1.615238228351744354 -1.2409464835010524109 -0.022841296189723057225;-0.88461633681805329488 -0.62280419656585817023 -1.5915411021917986734 -1.1387822291102893679 -0.027648729403380822162;-1.6216403281408999781 -0.74135722812334581278 -0.93983562658913755161 -0.92826463884786059033 0.00061904536944518509002;-0.36618067957480998631 -0.15063583652290962744 -1.7933353421312474563 -1.2598574033535114847 0.021549109000376355516;-0.40189828168531904762 -0.57329149245046506778 -0.86349313782959624675 -1.0475792100756353697 -0.24965422310581064647;-0.36734709885316202405 -0.38866583601559712147 -1.0360594594918590694 -0.82596710688332586336 -0.27154297131090199269;0.075824453877012115699 -0.54893133488366507855 -1.3826273572015228552 -1.2233108048173197968 -0.12848845628977242916;-0.38015303481192069812 -0.45493419779451366169 -1.2425306938466305073 -1.1940899134092877176 -0.17413720501443116784;-0.46695342488601127018 -0.24837630300124646787 -1.1559700040010012501 -0.85099512859844439472 -0.11908206929528938889;-1.2673873350419493367 -0.70193491945864305848 -1.3596430226779112616 -1.21587072963637155 -0.024680457613379697124;-1.1902933210880533021 -0.91896862390220113159 -1.0049929548139358015 -1.1962373892183888469 -0.073397738052385377849;-1.1079015365854369968 -1.0091466322842115311 -0.96111203208814743881 -1.219363788358288847 -0.10812137818639445952;-1.1637095888485931017 -1.0715807194056519158 -0.78892716409783336662 -1.1668837179454243369 -0.11785835291933749458;-1.3251505489492185763 -0.924459228832765878 -0.73418137740146638404 -0.99769365110999996826 -0.071810627285874728543;-1.2668005731053333385 -0.91248115595440359993 -0.8602633856476468166 -1.0707718090370765829 -0.075667096192046048864;-1.1907714825276181081 -0.94029011268394790957 -1.0333210529395073785 -1.2168157663370253729 -0.081667664471256506964;-1.2389948014342890836 -0.83125107270870046694 -1.2175406970762663317 -1.2495257503391012932 -0.066083694670206083455;-1.1094201152527534138 -1.0193240560134642081 -0.89228876641007637271 -1.2303779816337994113 -0.11529604197049289627;-1.1252260953276669309 -0.98349974499924430749 -0.94228651900950799458 -1.2516861727631218848 -0.10106234610813072394;-1.2240318832021808149 -0.92715006085048723428 -0.90741357033067526761 -1.0921037616786029023 -0.076181034761230920216;-1.3123637928717475365 -0.92298242899923166327 -0.73787713261067622295 -1.0004867321151502679 -0.075646590244059425112;-1.2429777255698202243 -0.94939782896672586521 -0.88658299752763547463 -1.1014252188232600549 -0.071193757718160777248;-1.1654797689885376677 -1.0201297927978287827 -0.93423986725718077828 -1.2091046103391160926 -0.08987242927598454445;-1.0340609636859301457 -1.0207243680334743274 -1.3009789971324168256 -1.2131916542383289759 -0.038659223102722004017;-1.0785197935810491021 -0.90975199724747080054 -1.2518990231879505437 -1.1058314813775336205 -0.031662072065034602142;-1.1084407072049784482 -0.90899085042525840095 -1.3388776860012805425 -1.213455884130224316 -0.027366591676074927819;-1.1181169271024082867 -0.81510249573220860952 -1.4107670646006886805 -1.2576440500464376449 -0.038785898556538390969;-0.37932274520403674867 -0.43499438327777101287 -1.3843609552318201228 -1.2763084962744790651 -0.29800831213588746538;-0.81378987121341694877 -0.64435207970053876014 -1.4978654701126534743 -1.5828479514464661726 -0.38931631781751868671;-0.93099758056727277555 -0.59788745116111097211 -1.6938794167886792597 -1.7308167638131088939 -0.35159783842116981711;-0.98869685109131588696 -0.5779159582906587822 -1.6526643964719847091 -1.6855960921640660288 -0.35199758395083263762;-1.06265361662237523 -0.41414376014564679984 -1.0068998180221166994 -1.0423639672174938031 -0.30421093230796986129;-1.6800753449739573853 -0.96814238935909413541 -0.90602655728552439207 -1.2896413067265231067 -0.10519468187701139228;-1.459076896932991918 -0.92472387788198429082 -0.49385305027016274115 -0.9517711983980110535 -0.14286735267392725035;-1.1771256337919480117 -1.3327852526621060125 -0.58842073537084782497 -1.3375751175773353818 0.043607341663187269842;-1.7490602716035268394 -1.1360968640606725621 -0.43196216008306714818 -1.0826820568232415898 -0.13830328193551014171;-1.746573591868940456 -1.0309178510533316686 -0.32980137468653514254 -0.85926695305281663995 -0.10925722877483304285;-1.295635418507935066 -0.45381892974092929993 -0.91020741797925097671 -0.97611666210209557804 -0.26076135350676310498;-1.2048031366900611161 -0.47889603772901223611 -0.74327032959260430633 -0.84996832929238574916 -0.08754992731118417526;-2.059491932845358253 -0.86346624733304422339 -0.15504362297560858597 -0.72455523048984260548 -0.061966515076046196908;-1.6975491352018712288 -0.68305532696593895814 -0.60760544299704521976 -0.8502011482038346335 -0.058576949178369433491;-1.1539390205841202341 -0.48890544943968400338 -1.2014752966579260729 -1.0734120830487550169 -0.099453431847043527458;-1.4250962133391273934 -0.67635284841889597995 -0.93466713113895472098 -1.0673693658143921237 -0.063815244596548414968;-1.147105904927150144 -0.63373766580615276123 -1.2978836539669862393 -1.4093009614487537196 -0.077569622366531823521;-1.206146803308540516 -0.77018529912254551917 -1.286731575977215547 -1.3353829084044859243 -0.050392565905801919857;0.08882625012896092076 -0.33898163687965748858 -1.1976185988843390806 -0.60182184924932680126 -0.26833134908008854369;-1.1307634299863189575 -1.4930118865286226004 -0.71432943073387999977 -1.4193404737899006918 -0.13676495221987380879;-0.50259874330565756928 -2.400903445435958794 -0.10518609670274563228 -1.5930050974242233597 -0.11814356868165530601;-1.2015567697485955989 -2.1774201755225757182 -4.5069697554699832409 -4.2274157743589801939 1.3793215956861815119;-1.1919588962588916203 -2.0288078853283564307 -0.17053829888783325619 -1.4742473376687603714 -0.1205586409537567627;-1.5017158468327929111 -1.3173934730002143656 -0.5627918145208982148 -1.3601701886634320893 -0.1304985114635736354;-0.95747766283217250471 -0.56758685729680946608 -0.81216383481992870941 -1.0661215113942594801 -0.20234587511004997751;-1.5740291360410638699 -0.79528502555043856237 -0.74588796379100208433 -1.0925383727334532225 -0.079848206135031754571;-1.1249837919064316072 -0.56405191251958930465 -1.0040826651939667169 -1.14114751215979493 -0.10648626272002473248;-1.2622620227028051865 -0.49901109799451898175 -1.2060706511841803312 -1.3266845246972485128 -0.3432005899877621502;-1.2822750233579265355 -0.47439822587966601342 -0.920345511440747166 -0.99321992829823657623 -0.26877320391748327211;-2.3182934363456020321 -1.2838303319423636228 0.27034261828353589596 -0.62971293495916957905 -0.034716828599310280146;-1.629993316156614247 -1.0088671844344598849 -0.40588161973616132538 -0.84534437714692745036 -0.07048074835852781439;-1.5384349952299973818 -0.86888303196663185801 -0.67176448247637643441 -1.0225133666316037484 -0.11374881970735438441;-1.4407683638873953758 -0.66328579078794547286 -1.1156645444000599099 -1.1305739837683239646 -0.038086281923200314414;-1.3491394765800448496 -0.68363224739941741426 -1.3390432583858138393 -1.1143869811062880881 -0.01616054377680282983;-1.0860420907194741069 -0.61807414696071072413 -1.2498372384338931607 -1.1271222000584433243 -0.070102565606412051147;-1.8213339492935676223 -1.8503204628216409766 -0.34570273630576686674 -1.4569785576148550632 0.12291739951868521641;-1.2732946437899663739 -1.4954304690229429653 -0.7001389487367228881 -1.4984221464774025723 -0.026144366658615327687;-1.7254062744547422614 -1.6506516039224692349 -4.2650382406739213437 -3.6791860393205637614 1.0572923407488437064;-1.7244921413179552516 -1.6514974730322651286 -4.2654561561220543098 -3.6789783456944813445 1.057692528156316536;-2.0039034290314283737 -2.1765384261837699142 0.065029836358701381682 -1.389731291095099186 0.10803945916785294423;-1.3133055115276375258 -0.94679779293581900479 -0.64311872910372858136 -1.0647339997595688388 -0.24353642364530725395;-1.569156992906742154 -1.1026534176940632204 -0.7441015184560048068 -1.2612218684885245779 -0.15827305013650450061;-1.5321421849397018899 -0.88070078007257879271 -0.955509831720449343 -1.1892636191483347829 -0.11886520752907027421;0.32613344446973951962 -0.17466006905483971723 -0.93972286283638151172 -0.51153857885801101713 -0.076780282492857773735;0.42784081793607564581 0.067757656243531957929 -0.21051309315970384572 0.73516321378602134118 -0.22992496785472046894;-0.13118315714592584942 -0.61495251493364733175 -1.4807931913520042677 -1.4691962503017772512 -0.089519552360230153498;-1.4508034509976461823 -1.0604845427883320408 -0.79651231495795304927 -1.5586304882186516885 -0.08498994388731290206;-1.3113334982453748179 -1.1986846689598493754 -0.68087355939185900677 -1.5574102674737893004 -0.10012389775699381822;-1.3411796476576469672 -0.49596044787455079161 -1.0939673097325295359 -0.88801977249119978808 -0.01582996936328873519;-1.1852224210922750469 -0.26075064807429787495 -1.4821348400376432863 -1.3018546502355010652 -0.0080343491799572935158;-0.37350618124378776042 -0.73836262894285165093 -1.1419067582229573965 -0.86698774880571161106 -0.16195856807849590142;-1.5571041564496572729 -1.8371440058594763478 -4.9715224122203904855 -3.7224388202245575208 1.5377491609221707325];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [0.235903084493303;0.235903084493303;0.421148686641497;0.421148686641497;0.160813687419867;0.394434012229766;0.405317566585545;0.412186470099877;0.39063235758692;0.250441853933121;0.201888079897602;0.0370026059803452;0.0291338183635166;0.0305933840219871;0.0347067723644359;0.0264898574314871;0.0246786245176426;0.029555321748299;0.0349884174542744;0.0325114037459352;0.0250763259907073;0.0247509762531537;0.0296876737671831;0.0342150484921149;0.0312745092658894;0.0355614882582679;0.0426030394709984;0.0654925456959308;0.0746842194031279;0.0439650968008979;0.0412085390008489;0.0395129731830173;0.0307028711560215;0.160499995101929;0.162144174613758;0.14223302551192;0.175170229400405;0.137133721909569;0.0450212548465722;0.0354777362805002;0.0342947867591228;0.029900408771839;0.0312571541948383;0.0334673504895039;0.0381338192385048;0.0429021536834015;0.0305902424053595;0.0311815601723322;0.0342320625240766;0.0314671480976434;0.0339506452648159;0.0348642755818404;0.0441494381954348;0.0409297631374822;0.047298058019798;0.0474363001770465;0.010006086618699;0.010006086618699;0.010006086618699;0.010006086618699;0.0100055236793017;0.0160046976679129;0.0100038665457093;0.0100033175980834;0.0130430020179924;0.0125969204810024;0.0100055236793017;0.0100044147899738;0.0173467431612796;0.016639162598394;0.0149782004058805;0.0210821591713311;0.0217505751401872;0.0365313500924767;0.152765993156113;0.0325594081516763;0.0241878767937407;0.0100075847890697;0.0212969839157174;0.0174233079414018;0.0100030420017562;0.0158461624298142;0.0100044147899738;0.0100055236793017;0.0100055236793017;0.0180540031940787;0.0132496924929461;0.0147149998818395;0.0203064999948874;0.0376507032074279;0.0166542178065776;0.0178651249309555;0.0100044147899738;0.0100100763662396;0.0100100763662396;0.0165593427404429;0.0100049776045959;0.0143915015565714;0.0179628109568049;0.255767211873568;0.29400008523562;0.0349802843589247;0.0720223311427127;0.072573406989022;0.0836527090184886;0.0897088595985715;0.0936035346989086;0.0100139032754369];
y1_step1.xoffset = [191.52165222168;191.52165222168;198.334030151367;198.334030151367;192.177291870117;199.930374145508;199.999984741211;199.999984741211;200;196.498336791992;192.524353027344;145.915267944336;131.316772460938;134.591903686523;142.339874267578;124.464218139648;118.926544189453;132.316040039062;142.815155029297;138.450714111328;120.209014892578;119.15991973877;132.601104736328;141.512451171875;136.016418457031;143.730010986328;153.02766418457;169.440673828125;173.206329345703;154.480743408203;151.423309326172;149.351287841797;134.815292358398;187.531799316406;187.658157348633;185.931427001953;188.575393676758;185.408554077148;155.549362182617;143.59944152832;141.654922485352;133.084106445312;135.98747253418;140.213088989258;147.525939941406;153.355117797852;134.592498779297;135.83235168457;141.548065185547;136.414474487305;141.063781738281;142.607513427734;154.672134399414;151.108627319336;157.687789916992;157.811019897461;0;0;0;0;0;74.9703598022461;0;0;46.6081657409668;41.1719856262207;0;0;84.6716232299805;79.7464981079102;66.3843536376953;105.077911376953;108.004295349121;145.224960327148;186.900848388672;138.551788330078;117.28466796875;0;106.060707092285;85.1819534301758;0;73.7423629760742;0;0;0;89.1992416381836;48.9867706298828;64.017936706543;101.44303894043;146.835968017578;79.8220443725586;88.0058670043945;0;0;0;79.1780853271484;0;60.9627647399902;88.5925216674805;194.632080078125;198.723526000977;142.810348510742;172.227600097656;172.438461303711;176.076995849609;177.691024780273;178.61865234375;0];

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);

% Adjusts neural net estimates in order to get more realistic behavior: only works for a single time series ie will erase samples if many samples 
% compatible with SDP
if adjustOutput
    indexNegatives = y1 < 0;
    y1(indexNegatives) = 0;  % Remove negative values
    indexAfterFirstZero = ~cumprod(y1 ~= 0, 2);
    y1(indexAfterFirstZero) = 0;    % Once drawdown gets to zero stays there
    indexAboveMax = y1 > 200;
    y1(indexAboveMax) = 200;    % Can't be higher than starting head
end

end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
