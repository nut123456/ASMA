
% Developed in MATLAB R2018b
%
% ################# Design Optimization of StructuresDesign Optimization of Structures #############
% #                 
% # ---------------------------------------------------------------------------------------
% # Author: Arnut Sutha, Ph.D Cadidate
% # Advisor: Assoc.Prof.Sawekchai Tangaramvong
% # Applied Machanics and Structures Research Unit Department of Civil Engineering, 
% #         Chulalongkorn University, Bangkok, Thailand 
% # Date: 25st Jun 2024
% # Version: 0.0.1

function [c,ceq] = Cstr1410(x)
%x=[5.35,	3.4,	26.5938,	8.0668,	3.6831,	2.3176,	21.3495,	10.3539,	6.3859,	4.5076,	13.8374,	10.876,	3.0903,	5.599,	14.2288,	9.8232,	2.9863,	7.6236,	8.6671,	12.1283,	5.1126,	6.4542,	1.2039,	4.2053,	2.2692,	3.1514,	3.9958,	7.6868,	3.3384,	3.4703,	4.498,	8.4914,	6.07,	4.7849,	4.1645,	1.0546,	7.2711,	4.9319,	4.3437,	1.0317,	6.1891,	7.2491,	4.8191,	1,	7.0979,	4.3028,	1.0716]
E = 2e11;     % Young's elastic modulus (kg/m^2)
rho = 7850;   % density of material (kg/m^3)
A = zeros(1,47);       % area of bar (m^2)
A1 = x(1); A2 = x(2); A3 = x(3); A4 = x(4); A5 = x(5); A6 = x(6); A7 = x(7); A8 = x(8); A9 = x(9); A10 = x(10);
A11 = x(11); A12 = x(12); A13 = x(13); A14 = x(14); A15 = x(15); A16 = x(16); A17 = x(17); A18 = x(18); A19 = x(19); A20 = x(20);
A21 = x(21); A22 = x(22); A23 = x(23); A24 = x(24); A25 = x(25); A26 = x(26); A27 = x(27); A28 = x(28); A29 = x(29); A30 = x(30);
A31 = x(31); A32 = x(32); A33 = x(33); A34 = x(34); A35 = x(35); A36 = x(36); A37 = x(37); A38 = x(38); A39 = x(39); A40 = x(40);
A41 = x(41); A42 = x(42); A43 = x(43); A44 = x(44); A45 = x(45); A46 = x(46); A47 = x(47);
A([1;48;95;142;189;236;283;330;377;424;471;518;565;612;659;706;753;800;847;894;941;988;1035;1082;1129;1176;1223;1270;1317;1364])= A1;
A([2;49;96;143;190;237;284;331;378;425;472;519;566;613;660;707;754;801;848;895;942;989;1036;1083;1130;1177;1224;1271;1318;1365])= A2;
A([3;50;97;144;191;238;285;332;379;426;473;520;567;614;661;708;755;802;849;896;943;990;1037;1084;1131;1178;1225;1272;1319;1366])= A3;
A([4;51;98;145;192;239;286;333;380;427;474;521;568;615;662;709;756;803;850;897;944;991;1038;1085;1132;1179;1226;1273;1320;1367])= A4;
A([5;52;99;146;193;240;287;334;381;428;475;522;569;616;663;710;757;804;851;898;945;992;1039;1086;1133;1180;1227;1274;1321;1368])= A5;
A([6;53;100;147;194;241;288;335;382;429;476;523;570;617;664;711;758;805;852;899;946;993;1040;1087;1134;1181;1228;1275;1322;1369])= A6;
A([7;54;101;148;195;242;289;336;383;430;477;524;571;618;665;712;759;806;853;900;947;994;1041;1088;1135;1182;1229;1276;1323;1370])= A7;
A([8;55;102;149;196;243;290;337;384;431;478;525;572;619;666;713;760;807;854;901;948;995;1042;1089;1136;1183;1230;1277;1324;1371])= A8;
A([9;56;103;150;197;244;291;338;385;432;479;526;573;620;667;714;761;808;855;902;949;996;1043;1090;1137;1184;1231;1278;1325;1372])= A9;
A([10;57;104;151;198;245;292;339;386;433;480;527;574;621;668;715;762;809;856;903;950;997;1044;1091;1138;1185;1232;1279;1326;1373])= A10;
A([11;58;105;152;199;246;293;340;387;434;481;528;575;622;669;716;763;810;857;904;951;998;1045;1092;1139;1186;1233;1280;1327;1374])= A11;
A([12;59;106;153;200;247;294;341;388;435;482;529;576;623;670;717;764;811;858;905;952;999;1046;1093;1140;1187;1234;1281;1328;1375])= A12;
A([13;60;107;154;201;248;295;342;389;436;483;530;577;624;671;718;765;812;859;906;953;1000;1047;1094;1141;1188;1235;1282;1329;1376])= A13;
A([14;61;108;155;202;249;296;343;390;437;484;531;578;625;672;719;766;813;860;907;954;1001;1048;1095;1142;1189;1236;1283;1330;1377])= A14;
A([15;62;109;156;203;250;297;344;391;438;485;532;579;626;673;720;767;814;861;908;955;1002;1049;1096;1143;1190;1237;1284;1331;1378])= A15;
A([16;63;110;157;204;251;298;345;392;439;486;533;580;627;674;721;768;815;862;909;956;1003;1050;1097;1144;1191;1238;1285;1332;1379])= A16;
A([17;64;111;158;205;252;299;346;393;440;487;534;581;628;675;722;769;816;863;910;957;1004;1051;1098;1145;1192;1239;1286;1333;1380])= A17;
A([18;65;112;159;206;253;300;347;394;441;488;535;582;629;676;723;770;817;864;911;958;1005;1052;1099;1146;1193;1240;1287;1334;1381])= A18;
A([19;66;113;160;207;254;301;348;395;442;489;536;583;630;677;724;771;818;865;912;959;1006;1053;1100;1147;1194;1241;1288;1335;1382])= A19;
A([20;67;114;161;208;255;302;349;396;443;490;537;584;631;678;725;772;819;866;913;960;1007;1054;1101;1148;1195;1242;1289;1336;1383])= A20;
A([21;68;115;162;209;256;303;350;397;444;491;538;585;632;679;726;773;820;867;914;961;1008;1055;1102;1149;1196;1243;1290;1337;1384])= A21;
A([22;69;116;163;210;257;304;351;398;445;492;539;586;633;680;727;774;821;868;915;962;1009;1056;1103;1150;1197;1244;1291;1338;1385])= A22;
A([23;70;117;164;211;258;305;352;399;446;493;540;587;634;681;728;775;822;869;916;963;1010;1057;1104;1151;1198;1245;1292;1339;1386])= A23;
A([24;71;118;165;212;259;306;353;400;447;494;541;588;635;682;729;776;823;870;917;964;1011;1058;1105;1152;1199;1246;1293;1340;1387])= A24;
A([25;72;119;166;213;260;307;354;401;448;495;542;589;636;683;730;777;824;871;918;965;1012;1059;1106;1153;1200;1247;1294;1341;1388])= A25;
A([26;73;120;167;214;261;308;355;402;449;496;543;590;637;684;731;778;825;872;919;966;1013;1060;1107;1154;1201;1248;1295;1342;1389])= A26;
A([27;74;121;168;215;262;309;356;403;450;497;544;591;638;685;732;779;826;873;920;967;1014;1061;1108;1155;1202;1249;1296;1343;1390])= A27;
A([28;75;122;169;216;263;310;357;404;451;498;545;592;639;686;733;780;827;874;921;968;1015;1062;1109;1156;1203;1250;1297;1344;1391])= A28;
A([29;76;123;170;217;264;311;358;405;452;499;546;593;640;687;734;781;828;875;922;969;1016;1063;1110;1157;1204;1251;1298;1345;1392])= A29;
A([30;77;124;171;218;265;312;359;406;453;500;547;594;641;688;735;782;829;876;923;970;1017;1064;1111;1158;1205;1252;1299;1346;1393])= A30;
A([31;78;125;172;219;266;313;360;407;454;501;548;595;642;689;736;783;830;877;924;971;1018;1065;1112;1159;1206;1253;1300;1347;1394])= A31;
A([32;79;126;173;220;267;314;361;408;455;502;549;596;643;690;737;784;831;878;925;972;1019;1066;1113;1160;1207;1254;1301;1348;1395])= A32;
A([33;80;127;174;221;268;315;362;409;456;503;550;597;644;691;738;785;832;879;926;973;1020;1067;1114;1161;1208;1255;1302;1349;1396])= A33;
A([34;81;128;175;222;269;316;363;410;457;504;551;598;645;692;739;786;833;880;927;974;1021;1068;1115;1162;1209;1256;1303;1350;1397])= A34;
A([35;82;129;176;223;270;317;364;411;458;505;552;599;646;693;740;787;834;881;928;975;1022;1069;1116;1163;1210;1257;1304;1351;1398])= A35;
A([36;83;130;177;224;271;318;365;412;459;506;553;600;647;694;741;788;835;882;929;976;1023;1070;1117;1164;1211;1258;1305;1352;1399])= A36;
A([37;84;131;178;225;272;319;366;413;460;507;554;601;648;695;742;789;836;883;930;977;1024;1071;1118;1165;1212;1259;1306;1353;1400])= A37;
A([38;85;132;179;226;273;320;367;414;461;508;555;602;649;696;743;790;837;884;931;978;1025;1072;1119;1166;1213;1260;1307;1354;1401])= A38;
A([39;86;133;180;227;274;321;368;415;462;509;556;603;650;697;744;791;838;885;932;979;1026;1073;1120;1167;1214;1261;1308;1355;1402])= A39;
A([40;87;134;181;228;275;322;369;416;463;510;557;604;651;698;745;792;839;886;933;980;1027;1074;1121;1168;1215;1262;1309;1356;1403])= A40;
A([41;88;135;182;229;276;323;370;417;464;511;558;605;652;699;746;793;840;887;934;981;1028;1075;1122;1169;1216;1263;1310;1357;1404])= A41;
A([42;89;136;183;230;277;324;371;418;465;512;559;606;653;700;747;794;841;888;935;982;1029;1076;1123;1170;1217;1264;1311;1358;1405])= A42;
A([43;90;137;184;231;278;325;372;419;466;513;560;607;654;701;748;795;842;889;936;983;1030;1077;1124;1171;1218;1265;1312;1359;1406])= A43;
A([44;91;138;185;232;279;326;373;420;467;514;561;608;655;702;749;796;843;890;937;984;1031;1078;1125;1172;1219;1266;1313;1360;1407])= A44;
A([45;92;139;186;233;280;327;374;421;468;515;562;609;656;703;750;797;844;891;938;985;1032;1079;1126;1173;1220;1267;1314;1361;1408])= A45;
A([46;93;140;187;234;281;328;375;422;469;516;563;610;657;704;751;798;845;892;939;986;1033;1080;1127;1174;1221;1268;1315;1362;1409])= A46;
A([47;94;141;188;235;282;329;376;423;470;517;564;611;658;705;752;799;846;893;940;987;1034;1081;1128;1175;1222;1269;1316;1363;1410])= A47;

% number of elements
num_ele = 1410;
% elements nodes
ele_nod = [1 2;1 8;1 14;2 3;2 8;2 9;2 15;3 4;3 9;3 10;3 16;4 5;4 10;4 11;4 17;5 6;5 11;5 12;5 18;6 7;6 12;6 13;6 19;7 13;8 9;8 14;8 15;8 21;9 10;9 15;9 16;9 22;10 11;10 16;10 17;10 23;11 12;11 17;11 18;11 24;12 13;12 18;12 19;12 25;13 19;13 20;13 26;...
           14 15;14 21;14 27;15 16;15 21;15 22;15 28;16 17;16 22;16 23;16 29;17 18;17 23;17 24;17 30;18 19;18 24;18 25;18 31;19 20;19 25;19 26;19 32;20 26;21 22;21 27;21 28;21 34;22 23;22 28;22 29;22 35;23 24;23 29;23 30;23 36;24 25;24 30;24 31;24 37;25 26;25 31;25 32;25 38;26 32;26 33;26 39;...
           27 28;27 34;27 40;28 29;28 34;28 35;28 41;29 30;29 35;29 36;29 42;30 31;30 36;30 37;30 43;31 32;31 37;31 38;31 44;32 33;32 38;32 39;32 45;33 39;34 35;34 40;34 41;34 47;35 36;35 41;35 42;35 48;36 37;36 42;36 43;36 49;37 38;37 43;37 44;37 50;38 39;38 44;38 45;38 51;39 45;39 46;39 52;...
           40 41;40 47;40 53;41 42;41 47;41 48;41 54;42 43;42 48;42 49;42 55;43 44;43 49;43 50;43 56;44 45;44 50;44 51;44 57;45 46;45 51;45 52;45 58;46 52;47 48;47 53;47 54;47 60;48 49;48 54;48 55;48 61;49 50;49 55;49 56;49 62;50 51;50 56;50 57;50 63;51 52;51 57;51 58;51 64;52 58;52 59;52 65;...
           53 54;53 60;53 66;54 55;54 60;54 61;54 67;55 56;55 61;55 62;55 68;56 57;56 62;56 63;56 69;57 58;57 63;57 64;57 70;58 59;58 64;58 65;58 71;59 65;60 61;60 66;60 67;60 73;61 62;61 67;61 68;61 74;62 63;62 68;62 69;62 75;63 64;63 69;63 70;63 76;64 65;64 70;64 71;64 77;65 71;65 72;65 78;...
           66 67;66 73;66 79;67 68;67 73;67 74;67 80;68 69;68 74;68 75;68 81;69 70;69 75;69 76;69 82;70 71;70 76;70 77;70 83;71 72;71 77;71 78;71 84;72 78;73 74;73 79;73 80;73 86;74 75;74 80;74 81;74 87;75 76;75 81;75 82;75 88;76 77;76 82;76 83;76 89;77 78;77 83;77 84;77 90;78 84;78 85;78 91;...
           79 80;79 86;79 92;80 81;80 86;80 87;80 93;81 82;81 87;81 88;81 94;82 83;82 88;82 89;82 95;83 84;83 89;83 90;83 96;84 85;84 90;84 91;84 97;85 91;86 87;86 92;86 93;86 99;87 88;87 93;87 94;87 100;88 89;88 94;88 95;88 101;89 90;89 95;89 96;89 102;90 91;90 96;90 97;90 103;91 97;91 98;91 104;...
           92 93;92 99;92 105;93 94;93 99;93 100;93 106;94 95;94 100;94 101;94 107;95 96;95 101;95 102;95 108;96 97;96 102;96 103;96 109;97 98;97 103;97 104;97 110;98 104;99 100;99 105;99 106;99 112;100 101;100 106;100 107;100 113;101 102;101 107;101 108;101 114;102 103;102 108;102 109;102 115;103 104;103 109;103 110;103 116;104 110;104 111;104 117;...
           105 106;105 112;105 118;106 107;106 112;106 113;106 119;107 108;107 113;107 114;107 120;108 109;108 114;108 115;108 121;109 110;109 115;109 116;109 122;110 111;110 116;110 117;110 123;111 117;112 113;112 118;112 119;112 125;113 114;113 119;113 120;113 126;114 115;114 120;114 121;114 127;115 116;115 121;115 122;115 128;116 117;116 122;116 123;116 129;117 123;117 124;117 130;...
           118 119;118 125;118 131;119 120;119 125;119 126;119 132;120 121;120 126;120 127;120 133;121 122;121 127;121 128;121 134;122 123;122 128;122 129;122 135;123 124;123 129;123 130;123 136;124 130;125 126;125 131;125 132;125 138;126 127;126 132;126 133;126 139;127 128;127 133;127 134;127 140;128 129;128 134;128 135;128 141;129 130;129 135;129 136;129 142;130 136;130 137;130 143;...
           131 132;131 138;131 144;132 133;132 138;132 139;132 145;133 134;133 139;133 140;133 146;134 135;134 140;134 141;134 147;135 136;135 141;135 142;135 148;136 137;136 142;136 143;136 149;137 143;138 139;138 144;138 145;138 151;139 140;139 145;139 146;139 152;140 141;140 146;140 147;140 153;141 142;141 147;141 148;141 154;142 143;142 148;142 149;142 155;143 149;143 150;143 156;...
           144 145;144 151;144 157;145 146;145 151;145 152;145 158;146 147;146 152;146 153;146 159;147 148;147 153;147 154;147 160;148 149;148 154;148 155;148 161;149 150;149 155;149 156;149 162;150 156;151 152;151 157;151 158;151 164;152 153;152 158;152 159;152 165;153 154;153 159;153 160;153 166;154 155;154 160;154 161;154 167;155 156;155 161;155 162;155 168;156 162;156 163;156 169;...
           157 158;157 164;157 170;158 159;158 164;158 165;158 171;159 160;159 165;159 166;159 172;160 161;160 166;160 167;160 173;161 162;161 167;161 168;161 174;162 163;162 168;162 169;162 175;163 169;164 165;164 170;164 171;164 177;165 166;165 171;165 172;165 178;166 167;166 172;166 173;166 179;167 168;167 173;167 174;167 180;168 169;168 174;168 175;168 181;169 175;169 176;169 182;...
           170 171;170 177;170 183;171 172;171 177;171 178;171 184;172 173;172 178;172 179;172 185;173 174;173 179;173 180;173 186;174 175;174 180;174 181;174 187;175 176;175 181;175 182;175 188;176 182;177 178;177 183;177 184;177 190;178 179;178 184;178 185;178 191;179 180;179 185;179 186;179 192;180 181;180 186;180 187;180 193;181 182;181 187;181 188;181 194;182 188;182 189;182 195;...
           183 184;183 190;183 196;184 185;184 190;184 191;184 197;185 186;185 191;185 192;185 198;186 187;186 192;186 193;186 199;187 188;187 193;187 194;187 200;188 189;188 194;188 195;188 201;189 195;190 191;190 196;190 197;190 203;191 192;191 197;191 198;191 204;192 193;192 198;192 199;192 205;193 194;193 199;193 200;193 206;194 195;194 200;194 201;194 207;195 201;195 202;195 208;...
           196 197;196 203;196 209;197 198;197 203;197 204;197 210;198 199;198 204;198 205;198 211;199 200;199 205;199 206;199 212;200 201;200 206;200 207;200 213;201 202;201 207;201 208;201 214;202 208;203 204;203 209;203 210;203 216;204 205;204 210;204 211;204 217;205 206;205 211;205 212;205 218;206 207;206 212;206 213;206 219;207 208;207 213;207 214;207 220;208 214;208 215;208 221;...
           209 210;209 216;209 222;210 211;210 216;210 217;210 223;211 212;211 217;211 218;211 224;212 213;212 218;212 219;212 225;213 214;213 219;213 220;213 226;214 215;214 220;214 221;214 227;215 221;216 217;216 222;216 223;216 229;217 218;217 223;217 224;217 230;218 219;218 224;218 225;218 231;219 220;219 225;219 226;219 232;220 221;220 226;220 227;220 233;221 227;221 228;221 234;...
           222 223;222 229;222 235;223 224;223 229;223 230;223 236;224 225;224 230;224 231;224 237;225 226;225 231;225 232;225 238;226 227;226 232;226 233;226 239;227 228;227 233;227 234;227 240;228 234;229 230;229 235;229 236;229 242;230 231;230 236;230 237;230 243;231 232;231 237;231 238;231 244;232 233;232 238;232 239;232 245;233 234;233 239;233 240;233 246;234 240;234 241;234 247;...
           235 236;235 242;235 248;236 237;236 242;236 243;236 249;237 238;237 243;237 244;237 250;238 239;238 244;238 245;238 251;239 240;239 245;239 246;239 252;240 241;240 246;240 247;240 253;241 247;242 243;242 248;242 249;242 255;243 244;243 249;243 250;243 256;244 245;244 250;244 251;244 257;245 246;245 251;245 252;245 258;246 247;246 252;246 253;246 259;247 253;247 254;247 260;...
           248 249;248 255;248 261;249 250;249 255;249 256;249 262;250 251;250 256;250 257;250 263;251 252;251 257;251 258;251 264;252 253;252 258;252 259;252 265;253 254;253 259;253 260;253 266;254 260;255 256;255 261;255 262;255 268;256 257;256 262;256 263;256 269;257 258;257 263;257 264;257 270;258 259;258 264;258 265;258 271;259 260;259 265;259 266;259 272;260 266;260 267;260 273;...
           261 262;261 268;261 274;262 263;262 268;262 269;262 275;263 264;263 269;263 270;263 276;264 265;264 270;264 271;264 277;265 266;265 271;265 272;265 278;266 267;266 272;266 273;266 279;267 273;268 269;268 274;268 275;268 281;269 270;269 275;269 276;269 282;270 271;270 276;270 277;270 283;271 272;271 277;271 278;271 284;272 273;272 278;272 279;272 285;273 279;273 280;273 285;...
           274 275;274 281;274 287;275 276;275 281;275 282;275 288;276 277;276 282;276 283;276 289;277 278;277 283;277 284;277 290;278 279;278 284;278 285;278 291;279 280;279 285;279 286;279 292;280 286;281 282;281 287;281 288;281 294;282 283;282 288;282 289;282 295;283 284;283 289;283 290;283 296;284 285;284 290;284 291;284 297;285 286;285 291;285 292;285 298;286 292;286 293;286 299;...
           287 288;287 294;287 300;288 289;288 294;288 295;288 301;289 290;289 295;289 296;289 302;290 291;290 296;290 297;290 303;291 292;291 297;291 298;291 304;292 293;292 298;292 299;292 305;293 299;294 295;294 300;294 301;294 307;295 296;295 301;295 302;295 308;296 297;296 302;296 303;296 309;297 298;297 303;297 304;297 310;298 299;298 304;298 305;298 311;299 305;299 306;299 312;...
           300 301;300 307;300 313;301 302;301 307;301 308;301 314;302 303;302 308;302 309;302 315;303 304;303 309;303 310;303 316;304 305;304 310;304 311;304 317;305 306;305 311;305 312;305 318;306 312;307 308;307 313;307 314;307 320;308 309;308 314;308 315;308 321;309 310;309 315;309 316;309 322;310 311;310 316;310 317;310 323;311 312;311 317;311 318;311 324;312 318;312 319;312 325;...
           313 314;313 320;313 326;314 315;314 320;314 321;314 327;315 316;315 321;315 322;315 328;316 317;316 322;316 323;316 329;317 318;317 323;317 324;317 330;318 319;318 324;318 325;318 331;319 325;320 321;320 326;320 327;320 333;321 322;321 327;321 328;321 334;322 323;322 328;322 329;322 335;323 324;323 329;323 330;323 336;324 325;324 330;324 331;324 337;325 331;325 332;325 338;...
           326 327;326 333;326 339;327 328;327 333;327 334;327 340;328 329;328 334;328 335;328 341;329 330;329 335;329 336;329 342;330 331;330 336;330 337;330 343;331 332;331 337;331 338;331 344;332 338;333 334;333 339;333 340;333 346;334 335;334 340;334 341;334 347;335 336;335 341;335 342;335 348;336 337;336 342;336 343;336 349;337 338;337 343;337 344;337 350;338 344;338 345;338 351;...
           339 340;339 346;339 352;340 341;340 346;340 347;340 353;341 342;341 347;341 348;341 354;342 343;342 348;342 349;342 355;343 344;343 349;343 350;343 356;344 345;344 350;344 351;344 357;345 351;346 347;346 352;346 353;346 359;347 348;347 353;347 354;347 360;348 349;348 354;348 355;348 361;349 350;349 355;349 356;349 362;350 351;350 356;350 357;350 363;351 357;351 358;351 364;...
           352 353;352 359;352 365;353 354;353 359;353 360;353 366;354 355;354 360;354 361;354 367;355 356;355 361;355 362;355 368;356 357;356 362;356 363;356 369;357 358;357 363;357 364;357 370;358 364;359 360;359 365;359 366;359 372;360 361;360 366;360 367;360 373;361 362;361 367;361 368;361 374;362 363;362 368;362 369;362 375;363 364;363 369;363 370;363 376;364 370;364 371;364 377;...
           365 366;365 372;365 378;366 367;366 372;366 373;366 379;367 368;367 373;367 374;367 380;368 369;368 374;368 375;368 381;369 370;369 375;369 376;369 382;370 371;370 376;370 377;370 383;371 377;372 373;372 378;372 379;372 385;373 374;373 379;373 380;373 386;374 375;374 380;374 381;374 387;375 376;375 381;375 382;375 388;376 377;376 382;376 383;376 389;377 383;377 384;377 390;...
           378 379;378 385;378 1;379 380;379 385;379 386;379 2;380 381;380 386;380 387;380 3;381 382;381 387;381 388;381 4;382 383;382 388;382 389;382 5;383 384;383 389;383 390;383 6;384 390;385 386;385 1;385 2;385 8;386 387;386 2;386 3;386 9;387 388;387 3;387 4;387 10;388 389;388 4;388 5;388 11;389 390;389 5;389 6;389 12;390 6;390 7;390 13];
% number of nodes
num_nod = 390;
% nodes coordinates
nod_coor = [1 0 4;3 0 3.75;5 0 3.25;7 0 2.75;9 0 2;11 0 1.25;13 0 0;1.989 0.209 3;3.978 0.418 2.75;5.967 0.627 2.25;7.956 0.836 1.75;9.945 1.0453 1;11.934 1.2543 -0.5;...
            0.9782 0.2079 4;2.9344 0.6237 3.75;4.8907 1.0396 3.25;6.847 1.4554 2.75;8.8033 1.8712 2;10.7596 2.287 1.25;12.7159 2.7029 0;1.9021 0.618 3;3.8042 1.2359 2.75;5.7063 1.8539 2.25;7.6083 2.4719 1.75;9.5104 3.0901 1;11.4124 3.7081 -0.5;...
            0.9136 0.4067 4;2.7406 1.2202 3.75;4.5677 2.0337 3.25;6.3948 2.8472 2.75;8.2219 3.6606 2;10.049 4.4741 1.25;11.8761 5.2876 0;1.732 0.9999 3;3.4641 1.9999 2.75;5.1961 2.9998 2.25;6.9281 3.9997 1.75;8.6601 4.9999 1;10.3921 5.9999 -0.5;...
            0.809 0.5878 4;2.4271 1.7634 3.75;4.0451 2.9389 3.25;5.6631 4.1145 2.75;7.2812 5.2901 2;8.8992 6.4656 1.25;10.5172 7.6412 0;1.4863 1.3382 3;2.9726 2.6764 2.75;4.4589 4.0146 2.25;5.9452 5.3528 1.75;7.4313	6.6912 1;8.9176 8.0294 -0.5;...
            0.6691 0.7431 4;2.0074 2.2294 3.75;3.3457 3.7157 3.25;4.6839 5.2020 2.75;6.0222 6.6883 2;7.3604	8.1746 1.25;8.6987 9.6609 0;1.1756 1.6180 3;2.3512 3.2359 2.75;3.5268 4.8539 2.25;4.7023 6.4719	1.75;5.8777	8.09 1;7.0533 9.708 -0.5;...
            0.5 0.866 4;1.5 2.5981 3.75;2.5 4.3301 3.25;3.5 6.0622 2.75;4.5 7.7942 2;5.5 9.5263 1.25;6.5 11.2583 0;0.8135 1.827 3;1.627 3.6541 2.75;2.4405 5.4811 2.25;3.2540 7.3081 1.75;4.0672 9.1353 1;4.8807 10.9623 -0.5;...
            0.309 0.9511 4;0.9271 2.8532 3.75;1.5451 4.7553 3.25;2.1631 6.6574 2.75;2.7812 8.5595 2;3.3992 10.4616 1.25;4.0172 12.3637 0;0.4159 1.9562 3;0.8317 3.9125 2.75;1.2476 5.8687 2.25;1.6635 7.8249 1.75;2.079 9.7813 1;2.4949 11.7375 -0.5;...
            0.1045 0.9945 4;0.3136 2.9836 3.75;0.5226 4.9726 3.25;0.7317 6.9617 2.75;0.9408 8.9507 2;1.1498 10.9397 1.25;1.3589 12.9288 0;0.0001 2 3;0.0001 3.9999 2.75;0.0002 5.9999 2.25;0.0002 7.9998 1.75;0 9.9998 1;0 11.9997 -0.5;...
            -0.1045 0.9945 4;-0.3136 2.9836 3.75;-0.5226 4.9726 3.25;-0.7317 6.9617 2.75;-0.9408 8.9507 2;-1.1498 10.9397 1.25;-1.3589 12.9288 0;-0.4158 1.9563 3;-0.8315 3.9125 2.75;-1.2473 5.8688 2.25;-1.6631 7.825 1.75;-2.0791 9.7813 1;-2.4949 11.7375 -0.5;...
            -0.309 0.9511 4;-0.9271 2.8532 3.75;-1.5451 4.7553 3.25;-2.1631 6.6574 2.75;-2.7812 8.5595 2;-3.3992 10.4616 1.25;-4.0172 12.3637 0;-0.8134 1.8271 3;-1.6268 3.6541 2.75;-2.4402 5.4812 2.25;-3.2536 7.3083 1.75;-4.0673 9.1352 1;-4.8807 10.9623 -0.5;...
            -0.5 0.866 4;-1.5 2.5981 3.75;-2.5 4.3301 3.25;-3.5 6.0622 2.75;-4.5 7.7942 2;-5.5 9.5263 1.25;-6.5 11.2583 0;-1.1755 1.6180 3;-2.3510 3.2361 2.75;-3.5265 4.8541 2.25;-4.702 6.4721 1.75;-5.8778 8.09 1;-7.0533 9.708 -0.5;...
            -0.6691 0.7431 4;-2.0074 2.2294 3.75;-3.3457 3.7157 3.25;-4.6839 5.2020 2.75;-6.0222 6.6883 2;-7.3604 8.1746 1.25;-8.6987 9.6609 0;-1.4862 1.3383 3;-2.9724 2.6765 2.75;-4.4587 4.0148 2.25;-5.9449 5.3531 1.75;-7.4313 6.6911 1;-8.9175 8.0294 -0.5;...
            -0.809 0.5878 4;-2.4271 1.7634 3.75;-4.0451 2.9389 3.25;-5.6631 4.1145 2.75;-7.2812 5.2901 2;-8.8992 6.4656 1.25;-10.5172 7.6412 0;-1.7320 1 3;-3.464 2 2.75;-5.196 3.0001 2.25;-6.9279 4.0001 1.75;-8.6601 4.9999 1;-10.3921 5.9999 -0.5;...
            -0.9136 0.4067 4;-2.7406 1.2202 3.75;-4.5677 2.0337 3.25;-6.3948 2.8472 2.75;-8.2219 3.6606 2;-10.049 4.4741 1.25;-11.8761 5.2876 0;-1.9021 0.6181 3;-3.8041 1.2361 2.75;-5.7062 1.8542 2.25;-7.6082 2.4723 1.75;-9.5104 3.0901 1;-11.4124 3.7081 -0.5;...
            -0.9782 0.2079 4;-2.9344 0.6237 3.75;-4.8907 1.0396 3.25;-6.847 1.4554 2.75;-8.8033 1.8712 2;-10.7596 2.287 1.25;-12.7159 2.7029 0;-1.9890 0.2091 3;-3.978 0.4182 2.75;-5.967 0.6273 2.25;-7.956 0.8364 1.75;-9.945 1.0452 1;-11.934 1.2543 -0.5;...
            -1 0 4;-3 0 3.75;-5 0 3.25;-7 0 2.75;-9 0 2;-11 0 1.25;-13 0 0;-1.989 -0.209 3;-3.978 -0.418 2.75;-5.967 -0.627 2.25;-7.956 -0.836 1.75;-9.945 -1.0453 1;-11.934 -1.2543 -0.5;...
            -0.9782 -0.2079 4;-2.9344 -0.6237 3.75;-4.8907 -1.0396 3.25;-6.847 -1.4554 2.75;-8.8033 -1.8712 2;-10.7596 -2.287 1.25;-12.7159 -2.7029 0;-1.9021 -0.618 3;-3.8042 -1.2359 2.75;-5.7063 -1.8539 2.25;-7.6083 -2.4719 1.75;-9.5104 -3.0901 1;-11.4124 -3.7081 -0.5;...
            -0.9136 -0.4067 4;-2.7406 -1.2202 3.75;-4.5677 -2.0337 3.25;-6.3948 -2.8472 2.75;-8.2219 -3.6606 2;-10.049 -4.4741 1.25;-11.8761 -5.2876 0;-1.732 -0.9999 3;-3.4641 -1.9999 2.75;-5.1961 -2.9998 2.25;-6.9281 -3.9997 1.75;-8.6601 -4.9999 1;-10.3921 -5.9999 -0.5;...
            -0.809 -0.5878 4;-2.4271 -1.7634 3.75;-4.0451 -2.9389 3.25;-5.6631 -4.1145 2.75;-7.2812 -5.2901 2;-8.8992 -6.4656 1.25;-10.5172 -7.6412 0;-1.4863 -1.3382 3;-2.9726 -2.6764 2.75;-4.4589 -4.0146 2.25;-5.9452 -5.3528 1.75;-7.4313	-6.6912 1;-8.9176 -8.0294 -0.5;...
            -0.6691 -0.7431 4;-2.0074 -2.2294 3.75;-3.3457 -3.7157 3.25;-4.6839 -5.2020 2.75;-6.0222 -6.6883 2;-7.3604	-8.1746 1.25;-8.6987 -9.6609 0;-1.1756 -1.6180 3;-2.3512 -3.2359 2.75;-3.5268 -4.8539 2.25;-4.7023 -6.4719 1.75;-5.8777	-8.09 1;-7.0533 -9.708 -0.5;...
            -0.5 -0.866 4;-1.5 -2.5981 3.75;-2.5 -4.3301 3.25;-3.5 -6.0622 2.75;-4.5 -7.7942 2;-5.5 -9.5263 1.25;-6.5 -11.2583 0;-0.8135 -1.827 3;-1.627 -3.6541 2.75;-2.4405 -5.4811 2.25;-3.2540 -7.3081 1.75;-4.0672 -9.1353 1;-4.8807 -10.9623 -0.5;...
            -0.309 -0.9511 4;-0.9271 -2.8532 3.75;-1.5451 -4.7553 3.25;-2.1631 -6.6574 2.75;-2.7812 -8.5595 2;-3.3992 -10.4616 1.25;-4.0172 -12.3637 0;-0.4159 -1.9562 3;-0.8317 -3.9125 2.75;-1.2476 -5.8687 2.25;-1.6635 -7.8249 1.75;-2.079 -9.7813 1;-2.4949 -11.7375 -0.5;...
            -0.1045 -0.9945 4;-0.3136 -2.9836 3.75;-0.5226 -4.9726 3.25;-0.7317 -6.9617 2.75;-0.9408 -8.9507 2;-1.1498 -10.9397 1.25;-1.3589 -12.9288 0;-0.0001 -2 3;-0.0001 -3.9999 2.75;-0.0002 -5.9999 2.25;-0.0002 -7.9998 1.75;0 -9.9998 1;0 -11.9997 -0.5;...
            0.1045 -0.9945 4;0.3136 -2.9836 3.75;0.5226 -4.9726 3.25;0.7317 -6.9617 2.75;0.9408 -8.9507 2;1.1498 -10.9397 1.25;1.3589 -12.9288 0;0.4158 -1.9563 3;0.8315 -3.9125 2.75;1.2473 -5.8688 2.25;1.6631 -7.825 1.75;2.0791 -9.7813 1;2.4949 -11.7375 -0.5;...
            0.309 -0.9511 4;0.9271 -2.8532 3.75;1.5451 -4.7553 3.25;2.1631 -6.6574 2.75;2.7812 -8.5595 2;3.3992 -10.4616 1.25;4.0172 -12.3637 0;0.8134 -1.8271 3;1.6268 -3.6541 2.75;2.4402 -5.4812 2.25;3.2536 -7.3083 1.75;4.0673 -9.1352 1;4.8807 -10.9623 -0.5;...
            0.5 -0.866 4;1.5 -2.5981 3.75;2.5 -4.3301 3.25;3.5 -6.0622 2.75;4.5 -7.7942 2;5.5 -9.5263 1.25;6.5 -11.2583 0;1.1755 -1.6180 3;2.351 -3.2361 2.75;3.5265 -4.8541 2.25;4.702 -6.4721 1.75;5.8778 -8.09 1;7.0533 -9.708 -0.5;...
            0.6691 -0.7431 4;2.0074 -2.2294 3.75;3.3457 -3.7157 3.25;4.6839 -5.2020 2.75;6.0222 -6.6883 2;7.3604 -8.1746 1.25;8.6987 -9.6609 0;1.4862 -1.3383 3;2.9724 -2.6765 2.75;4.4587 -4.0148 2.25;5.9449 -5.3531 1.75;7.4313 -6.6911 1;8.9175 -8.0294 -0.5;...
            0.809 -0.5878 4;2.4271 -1.7634 3.75;4.0451 -2.9389 3.25;5.6631 -4.1145 2.75;7.2812 -5.2901 2;8.8992 -6.4656 1.25;10.5172 -7.6412 0;1.7320 -1 3;3.464 -2 2.75;5.196 -3.0001 2.25;6.9279 -4.0001 1.75;8.6601 -4.9999 1;10.3921 -5.9999 -0.5;...
            0.9136 -0.4067 4;2.7406 -1.2202 3.75;4.5677 -2.0337 3.25;6.3948 -2.8472 2.75;8.2219 -3.6606 2;10.049 -4.4741 1.25;11.8761 -5.2876 0;1.9021 -0.6181 3;3.8041 -1.2361 2.75;5.7062 -1.8542 2.25;7.6082 -2.4723 1.75;9.5104 -3.0901 1;11.4124 -3.7081 -0.5;...
            0.9782 -0.2079 4;2.9344 -0.6237 3.75;4.8907 -1.0396 3.25;6.847 -1.4554 2.75;8.8033 -1.8712 2;10.7596 -2.287 1.25;12.7159 -2.7029 0;1.9890 -0.2091 3;3.978 -0.4182 2.75;5.967 -0.6273 2.25;7.956 -0.8364 1.75;9.945 -1.0452 1;11.934 -1.2543 -0.5];

dofPerNode = 3;
num_dof = num_nod*dofPerNode;
dof = 1:num_dof;

% calculate K and M
dofPerNode = 3;  % number of degree of freedom of one node
num_dof    = num_nod*dofPerNode; % total dgree of freedom of system
K = zeros(num_dof);
M = zeros(num_dof);
ele_dof = zeros(num_ele,6);
L = zeros(num_ele,1);
C = zeros(num_ele,3);
for ii=1:num_ele
   index1 = ele_nod(ii,1); 
   index2 = ele_nod(ii,2);
   dx = nod_coor(index2,1)-nod_coor(index1,1);
   dy = nod_coor(index2,2)-nod_coor(index1,2);
   dz = nod_coor(index2,3)-nod_coor(index1,3);
   % compute length of each bar
   L(ii) = sqrt(dx^2+dy^2+dz^2);
   C(ii,1) = dx/L(ii);
   C(ii,2) = dy/L(ii);
   C(ii,3) = dz/L(ii);
   ele_dof(ii,:) = [3*index1-2 3*index1-1 3*index1...
                    3*index2-2 3*index2-1 3*index2];
end
% Construct global stiffness matrix & Mass matrix
for ii = 1:num_ele
    index = ele_dof(ii,:);
    K(index,index) = K(index,index) + A(ii)*E/L(ii)*rotate3(C(ii,1),C(ii,2),C(ii,3)); % rotate() calls rotation matrix;
    M(index,index) = M(index,index) + rho*L(ii)*A(ii)*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0;0 1 0 0 2 0;0 0 1 0 0 2]/6;% lumped mass matrix
end

% boundary and loading conditions
dof_constrnt = [19 20 21 58 59 60 97 98 99 136 137 138 175 176 177 214 215 216 253 254 255 292 293 294 331 332 333 370 371 372 409 410 411 448 449 450 487 488 489 526 527 528 565 566 567 604 605 606 643 644 645 682 683 684 721 722 723 760 761 762 799 800 801 838 839 840 877 878 879 916 917 918 955 956 957 994 995 996 1033 1034 1035 1072 1073 1074 1111 1112 1113 1150 1151 1152]; % fixed constraints
% add non-structural mass 
dof_active = setdiff(dof,dof_constrnt);
             %  1  2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17 18 192021  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39 
add_M = diag([100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%1
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%2
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%3
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%4
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%5
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%6
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%7
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%8
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%9
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%10
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%11
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%12
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%13
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%14
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%15
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%16
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%17
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%18
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%19
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%20
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%21
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%22
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%23
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%24
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%25
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%26
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%27
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%28
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;...%29
              100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;0;0;0;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100]);%30
for dof = 1:num_dof
   M(dof,dof) = M(dof,dof) + add_M(dof,dof);
end

% Calculate frequency 
[omega,~]=eigens2s(K,M,dof_constrnt); %eigens2s(K,M,dof_constrnt); 
f=sqrt(omega)/(2*pi);
c1 = 7/f(1) -1;
c2 = 9/f(3)-1;
c = [c1 c2];
ceq = [];


function [L, X] = eigens2s(K, M, dof_constrnt)
% PURPOSE
%  Solve the generalized eigenvalue problem
%  [K - L*M]X = 0, considering boundary conditions.
%
% INPUT:
%    K : global stiffness matrix, dim(K) = ele_nod x ele_nod
%    M : global mass matrix, dim(M) = ele_nod x ele_nod
%    dof_constrnt : boundary condition matrix
%        dim(dof_active) = dof_active x 1
% OUTPUT:
%    L : eigenvalues stored in a vector with length (ele_nod - dof_active) 
%    X : eigenvectors dim(X) = ele_nod x num_dof, num_dof : number of dof's

[ele_nod, ~] = size(K);
dof = (1:ele_nod)';

if nargin == 3
    dof_constrnt = dof_constrnt(:);
    dof(dof_constrnt) = [];
    if nargout == 2
        % Use eigs to solve the generalized eigenvalue problem
        num_eigenvalues = min(length(dof), 6); % Change this number as needed
        [X1, D] = eigs(K(dof, dof), M(dof, dof), num_eigenvalues, 'smallestabs');
        d = diag(D);
        [L, ii] = sort(d);
        X2 = X1(:, ii);
        X = zeros(ele_nod, num_eigenvalues);
        X(dof, :) = X2;
        % Normalize eigenvectors
        for jj = 1:num_eigenvalues
            mnorm = sqrt(X(:, jj)' * M * X(:, jj));
            X(:, jj) = X(:, jj) / mnorm;
        end
    else
        d = eigs(K(dof, dof), M(dof, dof), 'smallestabs');
        L = sort(d);
    end
else
    if nargout == 2
        % Use eigs to solve the generalized eigenvalue problem
        num_eigenvalues = min(ele_nod, 6); % Change this number as needed
        [X1, D] = eigs(K, M, num_eigenvalues, 'smallestabs');
        d = diag(D);
        [L, ii] = sort(d);
        X = X1(:, ii);
        % Normalize eigenvectors
        for jj = 1:ele_nod
            mnorm = sqrt(X(:, jj)' * M * X(:, jj));
            X(:, jj) = X(:, jj) / mnorm;
        end
    else
        d = eigs(K, M, 'smallestabs');
        L = sort(d);
    end
end
%end


