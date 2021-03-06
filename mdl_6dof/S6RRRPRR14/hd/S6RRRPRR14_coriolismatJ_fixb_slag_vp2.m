% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRR14_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:34
% EndTime: 2019-03-09 20:13:46
% DurationCPUTime: 41.97s
% Computational Cost: add. (58292->1380), mult. (133801->1841), div. (0->0), fcn. (144181->10), ass. (0->678)
t833 = cos(qJ(6));
t1199 = t833 / 0.2e1;
t1023 = pkin(5) * t1199;
t829 = sin(qJ(6));
t1207 = t829 / 0.2e1;
t1024 = pkin(5) * t1207;
t1297 = m(7) * pkin(5);
t1029 = t1297 / 0.2e1;
t835 = cos(qJ(3));
t1139 = Ifges(6,3) * t835;
t1093 = qJ(4) * t835;
t1294 = pkin(3) + pkin(10);
t831 = sin(qJ(3));
t762 = t1294 * t831 - t1093;
t826 = t835 * pkin(9);
t794 = t835 * pkin(4) + t826;
t834 = cos(qJ(5));
t768 = t834 * t794;
t830 = sin(qJ(5));
t605 = -t762 * t830 + t768;
t1263 = t605 / 0.2e1;
t1296 = -mrSges(6,2) / 0.2e1;
t1137 = Ifges(7,3) * t835;
t1377 = mrSges(7,2) / 0.2e1;
t1380 = -mrSges(7,1) / 0.2e1;
t1047 = t833 * t834;
t1053 = t830 * t831;
t711 = t1047 * t831 - t1053 * t829;
t1142 = Ifges(7,6) * t711;
t900 = t829 * t834 + t833 * t830;
t714 = t900 * t831;
t1154 = Ifges(7,5) * t714;
t894 = -t1154 / 0.2e1 - t1142 / 0.2e1;
t519 = pkin(5) * t835 + t768 + (-pkin(11) * t831 - t762) * t830;
t1050 = t831 * t834;
t606 = t834 * t762 + t830 * t794;
t560 = pkin(11) * t1050 + t606;
t904 = t519 * t829 + t560 * t833;
t905 = t519 * t833 - t560 * t829;
t1314 = t894 - t1137 / 0.2e1 + t904 * t1377 + t905 * t1380;
t1145 = Ifges(6,6) * t834;
t1158 = Ifges(6,5) * t830;
t893 = t1145 / 0.2e1 + t1158 / 0.2e1;
t934 = -t835 * mrSges(7,2) + t711 * mrSges(7,3);
t936 = mrSges(7,1) * t835 - mrSges(7,3) * t714;
t1382 = t1139 / 0.2e1 + mrSges(6,1) * t1263 + t606 * t1296 + t934 * t1024 + t936 * t1023 + (t829 ^ 2 + t833 ^ 2) * t519 * t1029 + t893 * t831 - t1314;
t828 = sin(pkin(6));
t832 = sin(qJ(2));
t1060 = t828 * t832;
t1095 = cos(pkin(6));
t737 = t1060 * t835 + t1095 * t831;
t764 = t829 * t830 - t1047;
t1337 = t764 * t737;
t1378 = -mrSges(7,2) / 0.2e1;
t1379 = mrSges(7,1) / 0.2e1;
t882 = t900 * t737;
t1334 = t1337 * t1378 + t1379 * t882;
t1381 = -t1334 - (t1337 * t829 + t833 * t882) * t1297 / 0.2e1;
t1197 = t834 / 0.2e1;
t836 = cos(qJ(2));
t1059 = t828 * t836;
t736 = t1060 * t831 - t1095 * t835;
t615 = t1059 * t830 + t736 * t834;
t1040 = t834 * t836;
t616 = t1040 * t828 - t736 * t830;
t957 = t833 * t615 + t616 * t829;
t365 = Ifges(7,4) * t957;
t377 = t615 * t829 - t616 * t833;
t1354 = -Ifges(7,2) * t377 + t365;
t1376 = t1354 / 0.4e1;
t1375 = Ifges(5,4) - Ifges(4,5);
t1343 = -Ifges(4,6) + Ifges(5,5);
t989 = pkin(1) * t1095;
t950 = t832 * t989;
t740 = pkin(8) * t1059 + t950;
t686 = pkin(9) * t1095 + t740;
t687 = (-pkin(2) * t836 - pkin(9) * t832 - pkin(1)) * t828;
t487 = t686 * t831 - t835 * t687;
t388 = -pkin(4) * t737 - t487;
t813 = pkin(3) * t1059;
t326 = pkin(10) * t1059 - t388 + t813;
t1071 = t737 * qJ(4);
t738 = -pkin(8) * t1060 + t836 * t989;
t685 = -pkin(2) * t1095 - t738;
t864 = t685 - t1071;
t334 = t1294 * t736 + t864;
t159 = t326 * t830 + t334 * t834;
t132 = pkin(11) * t615 + t159;
t1086 = t132 * t829;
t158 = t834 * t326 - t334 * t830;
t131 = pkin(11) * t616 + t158;
t126 = pkin(5) * t737 + t131;
t86 = t126 * t833 - t1086;
t94 = t131 * t833 - t1086;
t1374 = -t86 + t94;
t1041 = t834 * t835;
t961 = -qJ(4) * t831 - pkin(2);
t746 = -t1294 * t835 + t961;
t793 = (pkin(4) + pkin(9)) * t831;
t597 = t746 * t834 + t793 * t830;
t551 = -pkin(11) * t1041 + t597;
t1077 = t551 * t829;
t1052 = t830 * t835;
t596 = -t746 * t830 + t834 * t793;
t550 = pkin(11) * t1052 + t596;
t517 = pkin(5) * t831 + t550;
t290 = t517 * t833 - t1077;
t306 = t550 * t833 - t1077;
t1373 = t290 - t306;
t959 = -mrSges(7,1) * t900 + t764 * mrSges(7,2);
t1372 = qJD(6) * t959;
t609 = -Ifges(7,5) * t900 + Ifges(7,6) * t764;
t1051 = t830 * t1294;
t780 = -t830 * pkin(11) - t1051;
t781 = (-pkin(11) - t1294) * t834;
t619 = t780 * t833 + t781 * t829;
t956 = -t780 * t829 + t833 * t781;
t114 = -t619 * mrSges(7,1) - t956 * mrSges(7,2) + t609;
t1371 = t114 * qJD(6);
t1353 = mrSges(7,1) * t377 + mrSges(7,2) * t957;
t1184 = pkin(5) * t615;
t1185 = pkin(4) * t736;
t1012 = qJ(4) * t1059;
t488 = t835 * t686 + t831 * t687;
t435 = t1012 - t488;
t341 = -t435 - t1185;
t258 = t341 - t1184;
t1370 = t258 * t1353;
t177 = t377 * Ifges(7,1) + t737 * Ifges(7,5) + t365;
t1369 = t1354 + t177;
t1125 = t377 * Ifges(7,4);
t1368 = Ifges(7,1) * t957 - t1125;
t1340 = Ifges(7,5) * t957;
t1361 = Ifges(7,6) * t377;
t199 = t1340 - t1361;
t1367 = Ifges(4,1) + Ifges(6,3) + Ifges(7,3);
t966 = t1361 / 0.2e1 - t1340 / 0.2e1;
t735 = pkin(5) * t1041 + t794;
t1230 = t735 / 0.2e1;
t1285 = t258 / 0.2e1;
t712 = t764 * t835;
t713 = t900 * t835;
t533 = -mrSges(7,1) * t713 + mrSges(7,2) * t712;
t691 = Ifges(7,4) * t712;
t502 = -t713 * Ifges(7,1) + Ifges(7,5) * t831 + t691;
t1268 = t502 / 0.4e1;
t536 = Ifges(7,2) * t713 + t691;
t965 = t536 / 0.4e1 + t1268;
t1366 = t1353 * t1230 + t533 * t1285 + t965 * t957;
t1365 = t377 / 0.4e1;
t1364 = t619 / 0.2e1;
t1224 = t737 / 0.4e1;
t1220 = t764 / 0.2e1;
t1219 = t764 / 0.4e1;
t1363 = -t828 / 0.2e1;
t1216 = -t900 / 0.4e1;
t1218 = t900 / 0.2e1;
t1362 = t957 / 0.4e1;
t819 = pkin(5) * t830 + qJ(4);
t1186 = m(7) * t819;
t1076 = t551 * t833;
t291 = t517 * t829 + t1076;
t305 = -t550 * t829 - t1076;
t1336 = t305 + t291;
t1246 = -t619 / 0.2e1;
t1358 = t377 * t1246;
t1161 = Ifges(7,4) * t764;
t612 = -Ifges(7,2) * t900 - t1161;
t1357 = -t377 * t612 / 0.4e1;
t389 = t488 - t1185;
t360 = t834 * t389;
t1094 = qJ(4) * t736;
t439 = t1294 * t737 + t1094;
t157 = -pkin(5) * t736 + t360 + (-pkin(11) * t737 - t439) * t830;
t1069 = t737 * t834;
t214 = t830 * t389 + t834 * t439;
t184 = pkin(11) * t1069 + t214;
t103 = t157 * t833 - t184 * t829;
t104 = t157 * t829 + t184 * t833;
t1138 = Ifges(7,3) * t736;
t1143 = Ifges(7,6) * t1337;
t1155 = Ifges(7,5) * t882;
t895 = t1155 / 0.2e1 - t1143 / 0.2e1;
t915 = t104 * t1378 + t103 * t1379 - t1138 / 0.2e1 + t895;
t1355 = t764 * t829 + t900 * t833;
t1233 = t713 / 0.2e1;
t1236 = -t712 / 0.2e1;
t1348 = -t957 / 0.2e1;
t1085 = t132 * t833;
t87 = t126 * t829 + t1085;
t1352 = t1233 * t87 + t1236 * t86 + t1348 * t290;
t1163 = Ifges(6,4) * t834;
t790 = -Ifges(6,2) * t830 + t1163;
t1170 = Ifges(6,1) * t830;
t930 = t1163 + t1170;
t1351 = (t790 + t930) * t1197;
t1038 = t835 * t836;
t1013 = t828 * t1038;
t739 = (pkin(2) * t832 - pkin(9) * t836) * t828;
t553 = -t831 * t738 + t739 * t835;
t398 = (pkin(4) * t1038 - t1294 * t832) * t828 - t553;
t1049 = t831 * t836;
t1014 = t828 * t1049;
t898 = pkin(3) * t1014 + t950;
t960 = pkin(8) - t1093;
t494 = (pkin(10) * t831 + t960) * t1059 + t898;
t225 = t834 * t398 - t494 * t830;
t672 = (t1049 * t830 + t832 * t834) * t828;
t183 = pkin(5) * t1013 - pkin(11) * t672 + t225;
t226 = t830 * t398 + t834 * t494;
t671 = (t1040 * t831 - t830 * t832) * t828;
t196 = pkin(11) * t671 + t226;
t110 = t183 * t833 - t196 * t829;
t111 = t183 * t829 + t196 * t833;
t455 = t671 * t833 - t672 * t829;
t1144 = Ifges(7,6) * t455;
t456 = t671 * t829 + t672 * t833;
t1156 = Ifges(7,5) * t456;
t1313 = t111 * t1377 + t110 * t1380 - t1144 / 0.2e1 - t1156 / 0.2e1;
t1349 = -t341 / 0.2e1;
t1249 = -t956 / 0.2e1;
t1248 = t956 / 0.2e1;
t1347 = -mrSges(4,1) + mrSges(5,2);
t1346 = mrSges(4,2) - mrSges(5,3);
t1345 = Ifges(5,1) + Ifges(4,3);
t1339 = t740 * mrSges(4,2);
t1335 = t435 + t488;
t1162 = Ifges(7,4) * t713;
t537 = Ifges(7,1) * t712 + t1162;
t1265 = t537 / 0.4e1;
t501 = Ifges(7,2) * t712 + Ifges(7,6) * t831 - t1162;
t964 = t1265 - t501 / 0.4e1;
t1333 = -t964 * t764 - t900 * t965;
t1194 = t835 / 0.2e1;
t1304 = t834 ^ 2;
t1332 = (t830 ^ 2 + t1304) * mrSges(6,3) * t1194;
t1020 = -t826 / 0.2e1;
t1181 = pkin(9) * t831;
t1021 = -t1181 / 0.2e1;
t1331 = t737 * t1020 + t736 * t1021;
t1330 = t619 * t1233 + t956 * t1236;
t1099 = t834 * mrSges(6,2);
t1104 = t830 * mrSges(6,1);
t1329 = t1104 / 0.2e1 + t1099 / 0.2e1;
t1327 = t1220 * t957;
t782 = -pkin(3) * t835 + t961;
t786 = mrSges(5,2) * t835 - mrSges(5,3) * t831;
t1326 = -m(5) * t782 - t786;
t1190 = -t1294 / 0.2e1;
t570 = -mrSges(6,2) * t1013 + mrSges(6,3) * t671;
t1325 = t570 * t1190 - t226 * mrSges(6,3) / 0.2e1;
t1323 = t1343 * t737 + t1375 * t736;
t1322 = Ifges(6,5) * t615 + Ifges(6,6) * t616 + t199;
t902 = t605 * t834 + t606 * t830;
t1321 = t103 * t764 - t104 * t900;
t176 = Ifges(7,2) * t957 + t737 * Ifges(7,6) + t1125;
t1320 = -t87 * mrSges(7,3) - t176 / 0.2e1;
t1164 = Ifges(6,4) * t830;
t792 = Ifges(6,1) * t834 - t1164;
t925 = Ifges(6,2) * t834 + t1164;
t1319 = t925 / 0.4e1 - t792 / 0.4e1;
t175 = t377 * Ifges(7,5) + Ifges(7,6) * t957 + t737 * Ifges(7,3);
t320 = -t616 * Ifges(6,5) + t615 * Ifges(6,6) + t737 * Ifges(6,3);
t720 = Ifges(4,4) * t736;
t512 = t737 * Ifges(4,1) - Ifges(4,5) * t1059 - t720;
t1318 = t512 + t320 + t175;
t1166 = Ifges(4,4) * t835;
t921 = Ifges(7,5) * t713 - Ifges(7,6) * t712;
t922 = t1145 + t1158;
t1317 = t1367 * t831 - t835 * t922 + t1166 - t921;
t1172 = mrSges(7,3) * t712;
t645 = -mrSges(7,2) * t831 + t1172;
t646 = mrSges(7,1) * t831 + mrSges(7,3) * t713;
t1316 = t1218 * t646 + t1220 * t645;
t1201 = t831 / 0.4e1;
t1209 = t819 / 0.2e1;
t1244 = t645 / 0.2e1;
t1106 = t900 * mrSges(7,2);
t1107 = t764 * mrSges(7,1);
t607 = -t1106 - t1107;
t1315 = t1201 * t609 + t1209 * t533 + t1230 * t607 + t956 * t1244;
t1312 = -t740 * mrSges(3,1) - t738 * mrSges(3,2);
t979 = t1059 / 0.2e1;
t946 = t835 * t979;
t1311 = Ifges(7,3) * t946 - t1313;
t1149 = Ifges(5,6) * t835;
t918 = t831 * Ifges(5,3) - t1149;
t927 = -t831 * Ifges(4,2) + t1166;
t1308 = t740 * mrSges(4,1) + (t832 * Ifges(4,6) + t836 * t927) * t1363 + (t832 * Ifges(5,5) + t836 * t918) * t828 / 0.2e1;
t1088 = t111 * t900;
t1090 = t110 * t764;
t1212 = t790 / 0.4e1;
t1221 = -t764 / 0.4e1;
t750 = Ifges(7,4) * t900;
t614 = -Ifges(7,1) * t764 - t750;
t1257 = t614 / 0.4e1;
t1259 = t612 / 0.4e1;
t495 = -pkin(3) * t1060 - t553;
t1272 = t495 / 0.2e1;
t1298 = m(7) / 0.2e1;
t1300 = m(6) / 0.2e1;
t1302 = m(5) / 0.2e1;
t228 = Ifges(7,4) * t456 + Ifges(7,2) * t455 + Ifges(7,6) * t1013;
t229 = Ifges(7,1) * t456 + Ifges(7,4) * t455 + Ifges(7,5) * t1013;
t253 = -mrSges(7,1) * t455 + mrSges(7,2) * t456;
t554 = t835 * t738 + t831 * t739;
t493 = -qJ(4) * t1060 - t554;
t432 = -pkin(4) * t1014 - t493;
t327 = -pkin(5) * t671 + t432;
t395 = -mrSges(7,2) * t1013 + mrSges(7,3) * t455;
t396 = mrSges(7,1) * t1013 - mrSges(7,3) * t456;
t489 = -mrSges(6,1) * t671 + mrSges(6,2) * t672;
t708 = (mrSges(5,1) * t1049 - mrSges(5,3) * t832) * t828;
t783 = mrSges(5,1) * t1013;
t709 = mrSges(5,2) * t1060 + t783;
t787 = t1099 + t1104;
t907 = t225 * t834 + t226 * t830;
t1307 = (-pkin(3) * t495 - qJ(4) * t493) * t1302 + (qJ(4) * t432 - t1294 * t907) * t1300 + (t110 * t956 + t111 * t619 + t327 * t819) * t1298 - pkin(3) * t709 / 0.2e1 - t327 * t959 / 0.2e1 + t432 * t787 / 0.2e1 + t455 * t1259 + t456 * t1257 - t493 * mrSges(5,3) / 0.2e1 + mrSges(5,2) * t1272 + t553 * mrSges(4,1) / 0.2e1 - t554 * mrSges(4,2) / 0.2e1 + t396 * t1248 + t395 * t1364 + t671 * t1212 + t672 * t792 / 0.4e1 + t229 * t1221 + t228 * t1216 + t253 * t1209 + (-t708 / 0.2e1 + t489 / 0.2e1) * qJ(4) + (-t1088 / 0.2e1 + t1090 / 0.2e1) * mrSges(7,3);
t611 = Ifges(7,2) * t764 - t750;
t613 = -Ifges(7,1) * t900 + t1161;
t1306 = t1216 * t177 + t1219 * t176 + t1224 * t609 + t1285 * t607 + t1362 * t611 + t1365 * t613;
t1303 = 2 * qJD(3);
t1301 = -m(6) / 0.2e1;
t1299 = -m(7) / 0.2e1;
t1295 = t86 / 0.2e1;
t1293 = t158 / 0.2e1;
t1292 = t159 / 0.2e1;
t1290 = t176 / 0.2e1;
t1289 = t176 / 0.4e1;
t1288 = t177 / 0.2e1;
t1287 = t177 / 0.4e1;
t1286 = -t225 / 0.2e1;
t1284 = t291 / 0.2e1;
t298 = -mrSges(7,2) * t737 + mrSges(7,3) * t957;
t1283 = t298 / 0.2e1;
t1282 = t305 / 0.2e1;
t1281 = t306 / 0.2e1;
t1279 = t957 / 0.2e1;
t1277 = -t377 / 0.2e1;
t1276 = t377 / 0.2e1;
t1108 = t737 * mrSges(6,2);
t1119 = t615 * mrSges(6,3);
t465 = -t1108 + t1119;
t1275 = t465 / 0.2e1;
t1109 = t737 * mrSges(6,1);
t1117 = t616 * mrSges(6,3);
t466 = t1109 + t1117;
t1274 = -t466 / 0.2e1;
t1271 = t501 / 0.2e1;
t1270 = t501 / 0.4e1;
t1269 = t502 / 0.2e1;
t1111 = t713 * mrSges(7,2);
t1112 = t712 * mrSges(7,1);
t534 = -t1111 - t1112;
t1267 = t534 / 0.2e1;
t1264 = t596 / 0.2e1;
t1262 = t959 / 0.2e1;
t1260 = t612 / 0.2e1;
t1258 = t614 / 0.2e1;
t1256 = -t615 / 0.2e1;
t1255 = t615 / 0.2e1;
t1254 = t615 / 0.4e1;
t1253 = -t616 / 0.2e1;
t1252 = -t616 / 0.4e1;
t1251 = t616 / 0.2e1;
t1250 = t616 / 0.4e1;
t1243 = -t646 / 0.2e1;
t1242 = t646 / 0.2e1;
t1241 = t672 / 0.2e1;
t1146 = Ifges(6,6) * t831;
t699 = -t835 * t925 + t1146;
t1240 = t699 / 0.2e1;
t1238 = -t711 / 0.2e1;
t1237 = t711 / 0.2e1;
t1235 = t712 / 0.4e1;
t1234 = -t713 / 0.4e1;
t1232 = t713 / 0.4e1;
t1231 = t714 / 0.2e1;
t1227 = -t737 / 0.2e1;
t1226 = -t737 / 0.4e1;
t1225 = t737 / 0.2e1;
t742 = t787 * t835;
t1223 = -t742 / 0.2e1;
t1222 = -t764 / 0.2e1;
t1217 = -t900 / 0.2e1;
t1102 = t831 * mrSges(6,1);
t769 = mrSges(6,3) * t1052 + t1102;
t1215 = t769 / 0.2e1;
t1101 = t831 * mrSges(6,2);
t770 = -mrSges(6,3) * t1041 - t1101;
t1214 = t770 / 0.2e1;
t1210 = t794 / 0.2e1;
t1208 = -t829 / 0.2e1;
t1206 = -t830 / 0.2e1;
t1205 = -t830 / 0.4e1;
t1204 = t830 / 0.2e1;
t1203 = -t831 / 0.4e1;
t1202 = t831 / 0.2e1;
t1200 = -t833 / 0.2e1;
t1198 = -t834 / 0.2e1;
t1196 = t834 / 0.4e1;
t1192 = -t836 / 0.2e1;
t1191 = t1294 / 0.2e1;
t1187 = m(7) * t735;
t1183 = pkin(5) * t616;
t1182 = pkin(5) * t834;
t1180 = t86 * mrSges(7,2);
t1179 = t86 * mrSges(7,3);
t1178 = t87 * mrSges(7,1);
t93 = -t131 * t829 - t1085;
t1176 = t93 * mrSges(7,1);
t1175 = t94 * mrSges(7,2);
t1174 = Ifges(4,4) + Ifges(5,6);
t1173 = m(7) * qJD(4);
t1171 = Ifges(6,1) * t615;
t1169 = Ifges(3,4) * t832;
t1168 = Ifges(4,4) * t737;
t1167 = Ifges(4,4) * t831;
t1165 = Ifges(6,4) * t616;
t1160 = Ifges(6,5) * t672;
t1159 = Ifges(6,5) * t737;
t1157 = Ifges(6,5) * t831;
t1153 = Ifges(7,5) * t764;
t1152 = Ifges(6,2) * t616;
t1151 = Ifges(5,6) * t737;
t1150 = Ifges(5,6) * t831;
t1148 = Ifges(6,6) * t671;
t1147 = Ifges(6,6) * t737;
t1141 = Ifges(7,6) * t900;
t1140 = Ifges(6,3) * t736;
t1136 = pkin(5) * qJD(5);
t1131 = t290 * mrSges(7,2);
t1130 = t291 * mrSges(7,1);
t1129 = t305 * mrSges(7,1);
t1128 = t306 * mrSges(7,2);
t1127 = t957 * mrSges(7,1);
t1126 = t377 * mrSges(7,2);
t198 = t1126 - t1127;
t227 = Ifges(7,3) * t1013 + t1144 + t1156;
t299 = mrSges(7,1) * t737 - mrSges(7,3) * t377;
t321 = t615 * Ifges(6,2) + t1147 - t1165;
t600 = Ifges(6,4) * t615;
t322 = -t616 * Ifges(6,1) + t1159 + t600;
t1118 = t616 * mrSges(6,2);
t1120 = t615 * mrSges(6,1);
t397 = -t1118 - t1120;
t417 = Ifges(6,3) * t1013 + t1148 + t1160;
t418 = Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t1013;
t419 = Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t1013;
t431 = t736 * pkin(3) + t864;
t436 = t487 + t813;
t509 = -Ifges(5,5) * t1059 + t736 * Ifges(5,3) - t1151;
t715 = Ifges(5,6) * t736;
t510 = -Ifges(5,4) * t1059 - t737 * Ifges(5,2) + t715;
t511 = -t736 * Ifges(4,2) - Ifges(4,6) * t1059 + t1168;
t564 = -mrSges(5,2) * t736 - mrSges(5,3) * t737;
t565 = t1059 * t960 + t898;
t571 = mrSges(6,1) * t1013 - mrSges(6,3) * t672;
t933 = Ifges(4,1) * t835 - t1167;
t624 = (t832 * Ifges(4,5) + t836 * t933) * t828;
t920 = -Ifges(5,2) * t835 + t1150;
t626 = (t832 * Ifges(5,4) + t836 * t920) * t828;
t1110 = t736 * mrSges(5,1);
t809 = mrSges(5,3) * t1059;
t641 = t809 + t1110;
t642 = t737 * mrSges(5,1) - mrSges(5,2) * t1059;
t643 = mrSges(4,2) * t1059 - t736 * mrSges(4,3);
t644 = -mrSges(4,1) * t1059 - t737 * mrSges(4,3);
t935 = -mrSges(5,2) * t831 - mrSges(5,3) * t835;
t680 = t935 * t1059;
t940 = mrSges(4,1) * t831 + mrSges(4,2) * t835;
t681 = t940 * t1059;
t706 = (-mrSges(4,2) * t832 - mrSges(4,3) * t1049) * t828;
t707 = (mrSges(4,1) * t832 - mrSges(4,3) * t1038) * t828;
t810 = Ifges(3,5) * t1059;
t886 = t836 * (Ifges(4,5) * t835 - Ifges(4,6) * t831);
t887 = t836 * (-Ifges(5,4) * t835 + Ifges(5,5) * t831);
t980 = -t1059 / 0.2e1;
t947 = t835 * t980;
t948 = t831 * t979;
t949 = t831 * t980;
t981 = t1060 / 0.2e1;
t5 = (t624 + t417 + t227) * t1225 + 0.2e1 * Ifges(3,4) * t1059 * t979 + t1318 * t946 + (-Ifges(3,6) * t1060 + t810 / 0.2e1 + Ifges(3,5) * t979 + t1312) * t1095 + (-pkin(1) * (mrSges(3,1) * t832 + mrSges(3,2) * t836) + t832 * (Ifges(3,1) * t836 - t1169) / 0.2e1 + (t1345 * t832 + t886 + t887) * t1192) * t828 ^ 2 + t737 * t1339 + (t1343 * t981 + t1308) * t736 + m(4) * (-t487 * t553 + t488 * t554 + t685 * t740) + m(5) * (t431 * t565 + t435 * t493 + t436 * t495) + m(6) * (t158 * t225 + t159 * t226 + t341 * t432) + m(7) * (t110 * t86 + t111 * t87 + t258 * t327) + t510 * t947 + t509 * t948 + t511 * t949 + ((Ifges(3,2) * t836 + t1169) * t1363 + (Ifges(3,1) - Ifges(3,2)) * t979) * t1060 - t487 * t707 + t435 * t708 + t436 * t709 + t488 * t706 + t685 * t681 + t431 * t680 + t671 * t321 / 0.2e1 + t493 * t641 + t495 * t642 + t554 * t643 + t553 * t644 + t565 * t564 + t159 * t570 + t158 * t571 + t341 * t489 + t226 * t465 + t225 * t466 + t432 * t397 + t87 * t395 + t86 * t396 + t327 * t198 + t110 * t299 + t111 * t298 + t258 * t253 + t228 * t1279 + t456 * t1288 + t455 * t1290 + t418 * t1255 + t229 * t1276 + t322 * t1241 + t419 * t1253 + t626 * t1227 + (-t1059 * t1345 - t1375 * t737) * t981;
t1122 = t5 * qJD(1);
t1017 = -Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1;
t1046 = t834 * t321;
t1057 = t830 * t322;
t213 = -t439 * t830 + t360;
t240 = Ifges(7,4) * t882 - Ifges(7,2) * t1337 - t736 * Ifges(7,6);
t241 = Ifges(7,1) * t882 - Ifges(7,4) * t1337 - t736 * Ifges(7,5);
t275 = mrSges(7,1) * t1337 + mrSges(7,2) * t882;
t990 = -pkin(4) - t1182;
t303 = t737 * t990 - t487;
t380 = mrSges(7,2) * t736 - mrSges(7,3) * t1337;
t381 = -mrSges(7,1) * t736 - mrSges(7,3) * t882;
t408 = -t736 * Ifges(6,6) + t737 * t925;
t409 = -t736 * Ifges(6,5) + t737 * t930;
t1100 = t834 * mrSges(6,1);
t1103 = t830 * mrSges(6,2);
t785 = t1100 - t1103;
t518 = t785 * t737;
t546 = -mrSges(6,3) * t737 * t830 - mrSges(6,1) * t736;
t547 = mrSges(6,2) * t736 + mrSges(6,3) * t1069;
t563 = pkin(3) * t737 + t1094;
t6 = t882 * t1288 - t1337 * t1290 + t563 * t564 + t158 * t546 + t159 * t547 - t341 * t518 + (t720 / 0.2e1 + t715 / 0.2e1 - t685 * mrSges(4,2) + t431 * mrSges(5,3) - t512 / 0.2e1 - t320 / 0.2e1 - t175 / 0.2e1 + t510 / 0.2e1 - t487 * mrSges(4,3) - t436 * mrSges(5,1)) * t736 + t214 * t465 + t213 * t466 + t388 * t397 + t87 * t380 + t86 * t381 + t103 * t299 + t303 * t198 + t104 * t298 + t258 * t275 + m(5) * (t431 * t563 + t435 * t487 + t436 * t488) + m(6) * (t158 * t213 + t159 * t214 + t341 * t388) + m(7) * (t103 * t86 + t104 * t87 + t258 * t303) + t240 * t1279 + t408 * t1255 + t241 * t1276 + t409 * t1253 + (t641 - t643) * t487 + (t642 - t644) * t488 + (-t511 / 0.2e1 + t685 * mrSges(4,1) - t431 * mrSges(5,2) + t509 / 0.2e1 + t1046 / 0.2e1 + t1057 / 0.2e1 - t488 * mrSges(4,3) + t435 * mrSges(5,1) + (-Ifges(4,4) / 0.2e1 - Ifges(5,6) / 0.2e1 + t893) * t737 + (-Ifges(4,1) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t1017) * t736 + t895) * t737 + t1323 * t828 * t1192;
t1121 = t6 * qJD(1);
t80 = t764 * t87;
t1105 = t900 * mrSges(7,3);
t1097 = t1294 * mrSges(6,3);
t931 = t1165 + t1171;
t939 = -mrSges(6,1) * t616 + mrSges(6,2) * t615;
t958 = t600 + t1152;
t9 = t198 * t1183 - t1370 - t93 * t299 + t957 * t1179 - t94 * t298 + t1368 * t1277 - m(7) * (-t1183 * t258 + t86 * t93 + t87 * t94) - t341 * t939 + t321 * t1253 + t931 * t1251 - t1320 * t377 + (t466 - t1117) * t159 + (t1119 - t465) * t158 + (t958 + t322) * t1256 + t1322 * t1227 + t1369 * t1348;
t1096 = t9 * qJD(1);
t1089 = t110 * t833;
t1087 = t111 * t829;
t14 = t1370 + t199 * t1225 - t87 * t299 + t86 * t298 + (t1368 / 0.2e1 + t1320) * t377 + (-t1179 + t1288 + t1354 / 0.2e1) * t957;
t1084 = t14 * qJD(1);
t1045 = t834 * t465;
t1056 = t830 * t466;
t909 = t158 * t830 - t159 * t834;
t23 = -t882 * t299 - t1337 * t298 - m(7) * (t1337 * t87 + t86 * t882) + (-m(5) * t435 + m(6) * t341 + m(7) * t258 + t198 + t397 - t641) * t1059 + (m(5) * t431 - m(6) * t909 + t1045 - t1056 + t564) * t737;
t1083 = t23 * qJD(1);
t1079 = t882 * t764;
t1078 = t1337 * t900;
t1073 = t616 * t959;
t1072 = t616 * t735;
t1066 = t764 * t714;
t1062 = t900 * t711;
t1058 = t829 * t713;
t700 = -t835 * t930 + t1157;
t1055 = t830 * t700;
t1054 = t830 * t769;
t1048 = t833 * t712;
t1044 = t834 * t699;
t1043 = t834 * t770;
t1039 = t834 * t1294;
t1037 = t1294 * t466;
t535 = Ifges(7,5) * t712 + Ifges(7,6) * t713;
t1022 = t1182 / 0.2e1;
t1019 = Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t1018 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t1016 = -t87 / 0.2e1 - t93 / 0.2e1;
t1015 = -t94 / 0.2e1 + t1295;
t999 = mrSges(7,3) * t1249;
t998 = mrSges(7,3) * t1364;
t997 = mrSges(7,3) * t1220;
t996 = t1105 / 0.2e1;
t995 = -t1105 / 0.2e1;
t992 = mrSges(6,3) * t1198;
t991 = -t1097 / 0.2e1;
t985 = t615 * t1198;
t984 = t534 * t1253;
t978 = -t1052 / 0.2e1;
t977 = -t1052 / 0.4e1;
t976 = -t1041 / 0.4e1;
t975 = -t1039 / 0.2e1;
t972 = t571 * t1190;
t971 = t769 * t1191;
t970 = t770 * t1190;
t969 = -t290 / 0.2e1 + t1281;
t968 = t1284 + t1282;
t967 = -t322 / 0.4e1 - t600 / 0.4e1;
t963 = -t613 / 0.4e1 + t1259;
t962 = t1257 + t611 / 0.4e1;
t953 = t830 * t991;
t743 = -Ifges(6,5) * t1041 + Ifges(6,6) * t1052;
t944 = t321 / 0.4e1 - t1165 / 0.4e1;
t741 = t785 * t835;
t942 = (-t534 / 0.2e1 - t741 / 0.2e1) * t836;
t937 = -mrSges(7,1) * t711 + mrSges(7,2) * t714;
t928 = Ifges(7,1) * t714 + Ifges(7,4) * t711;
t926 = Ifges(4,2) * t835 + t1167;
t923 = Ifges(7,4) * t714 + Ifges(7,2) * t711;
t919 = -Ifges(5,2) * t831 - t1149;
t917 = -Ifges(5,3) * t835 - t1150;
t916 = pkin(3) * t831 - t1093;
t610 = -t1141 - t1153;
t788 = Ifges(6,5) * t834 - Ifges(6,6) * t830;
t734 = (-pkin(9) + t990) * t831;
t885 = t785 * t831;
t896 = -mrSges(6,2) * t835 + mrSges(6,3) * t1050;
t897 = mrSges(6,1) * t835 - mrSges(6,3) * t1053;
t838 = (Ifges(7,5) * t835 + t928) * t1365 + (Ifges(7,6) * t835 + t923) * t1362 + (t916 * t431 + t782 * t563 + ((t436 - t487) * t835 + t1335 * t831) * pkin(9)) * t1302 + t1331 * mrSges(4,3) - (t927 + t1317) * t736 / 0.4e1 + (-Ifges(4,2) * t737 + t1318 - t720) * t835 / 0.4e1 + (t926 + t920) * t1226 + (t831 * t922 + t1044 + t1055 + t1137 + t1139 + t1142 + t1154 + t917 + t933) * t1224 + (Ifges(5,2) * t736 + t1151 + t511) * t1203 + (-Ifges(4,1) * t736 + t737 * t922 + t1046 - t1138 - t1140 - t1143 + t1155 - t1168 + t509) * t1201 - (Ifges(5,3) * t737 + t510 + t715) * t835 / 0.4e1 + (t919 + t918) * t736 / 0.4e1 + t882 * t1268 - t1337 * t1270 + t885 * t1349 + t916 * t564 / 0.2e1 + (-t887 / 0.4e1 - t886 / 0.4e1) * t828 + t644 * t1020 + t643 * t1021 - t793 * t397 / 0.2e1 + t563 * t786 / 0.2e1 + t782 * (-mrSges(5,2) * t737 + mrSges(5,3) * t736) / 0.2e1 + t388 * t741 / 0.2e1 - pkin(2) * (mrSges(4,1) * t737 - mrSges(4,2) * t736) / 0.2e1 + t734 * t198 / 0.2e1 + t597 * t547 / 0.2e1 + t290 * t381 / 0.2e1 + t408 * t976 + t409 * t977 + ((-t487 / 0.2e1 + t436 / 0.2e1) * t835 + t1331 + t1335 * t1202) * mrSges(5,1) + t87 * t934 / 0.2e1 + t431 * t935 / 0.2e1 + (t290 * t103 + t291 * t104 + t734 * t258 + t735 * t303 + t86 * t905 + t87 * t904) * t1298 + (t158 * t605 + t159 * t606 + t213 * t596 + t214 * t597 - t341 * t793 + t388 * t794) * t1300 + t904 * t1283 + t380 * t1284 + t937 * t1285 + t714 * t1287 + t711 * t1289 + t896 * t1292 + t897 * t1293 + t936 * t1295 + (Ifges(6,6) * t835 + t831 * t925) * t1254 + t466 * t1263 + t546 * t1264 + t303 * t1267 + t606 * t1275 + t103 * t1242 + t104 * t1244 + (Ifges(6,5) * t835 + t831 * t930) * t1252 + t275 * t1230 + t241 * t1234 + t240 * t1235 - t518 * t1210 + t214 * t1214 + t213 * t1215 + t905 * t299 / 0.2e1 + t642 * t826 / 0.2e1 + t641 * t1181 / 0.2e1 + t322 * t1053 / 0.4e1 + t685 * t940 / 0.2e1;
t2 = Ifges(5,4) * t947 + Ifges(4,5) * t946 + Ifges(5,5) * t948 + Ifges(4,6) * t949 + t419 * t1196 + t418 * t1205 + t225 * t992 + t834 * t972 - t838 + t1345 * t981 + (t788 + t610) * t1013 / 0.4e1 + t1325 * t830 + t1307;
t30 = -t714 * t502 / 0.2e1 + t501 * t1238 - t782 * t935 - t734 * t534 - t735 * t937 - t291 * t934 + t923 * t1236 + t793 * t741 + t794 * t885 + t928 * t1233 - t904 * t645 - t290 * t936 - t905 * t646 - t596 * t897 - t597 * t896 - t605 * t769 - t606 * t770 - m(7) * (t290 * t905 + t291 * t904 + t735 * t734) - m(6) * (t596 * t605 + t597 * t606 - t793 * t794) + (pkin(2) * mrSges(4,1) - t1044 / 0.2e1 - t1055 / 0.2e1 + (-t893 + t1174) * t831 + t894) * t831 + (pkin(2) * mrSges(4,2) + (t922 - t1174) * t835 + (Ifges(4,2) - Ifges(5,2) + Ifges(5,3) + Ifges(6,2) * t1304 / 0.2e1 + (t1163 + t1170 / 0.2e1) * t830 - t1367) * t831 + t921) * t835 + t1326 * t916;
t914 = -t2 * qJD(1) - t30 * qJD(2);
t744 = t835 * t790;
t745 = t835 * t792;
t839 = t1322 * t1201 + (-t744 + t700) * t1254 + t1366 + (t743 / 0.4e1 + t535 / 0.4e1) * t737 - t159 * t769 / 0.2e1 + t1369 * t1235 + t298 * t1281 + t299 * t1282 + t465 * t1264 + t597 * t1274 + t93 * t1242 + t94 * t1244 + t699 * t1250 - t745 * t1252 + t176 * t1232 + t1368 * t1234 + t341 * t1223 + t939 * t1210 + t158 * t1214 + (t1251 * t597 + t1256 * t596) * mrSges(6,3) + t1352 * mrSges(7,3) + (-mrSges(7,3) * t1284 + t1265 - t1270) * t377;
t850 = (mrSges(6,3) * t1293 - t1152 / 0.4e1 + t967) * t834 + (mrSges(6,3) * t1292 - t1171 / 0.4e1 + (t258 * t1299 - t198 / 0.2e1) * pkin(5) + t944) * t830;
t863 = -t1160 / 0.2e1 - t1148 / 0.2e1 + mrSges(6,1) * t1286 + t226 * mrSges(6,2) / 0.2e1;
t868 = t290 * t93 + t291 * t94 + t305 * t86 + t306 * t87;
t3 = (t1017 * t1059 + t850) * t835 + (t396 * t1200 + t395 * t1208 + t984 + 0.2e1 * (-t1072 / 0.4e1 - t1089 / 0.4e1 - t1087 / 0.4e1) * m(7)) * pkin(5) + t868 * t1298 + t839 + t863 + t1313;
t857 = t735 * t533 + (t536 / 0.2e1 + t1269) * t712 + (t291 * mrSges(7,3) + t1271 - t537 / 0.2e1) * t713 - t290 * t1172;
t39 = t306 * t645 + t305 * t646 + m(7) * (t290 * t305 + t291 * t306) - t597 * t769 + t596 * t770 - t794 * t742 + (t535 / 0.2e1 + t743 / 0.2e1) * t831 + ((t596 * mrSges(6,3) + t744 / 0.2e1 - t700 / 0.2e1) * t834 + (t597 * mrSges(6,3) + t1240 + t745 / 0.2e1 + (-t534 - t1187) * pkin(5)) * t830) * t835 + t857;
t913 = t3 * qJD(1) + t39 * qJD(2);
t842 = (t1376 + t1287) * t712 + (-t1368 / 0.4e1 + t1289) * t713 + t964 * t377 + (t1277 * t291 + t1352) * mrSges(7,3) + t290 * t1283 - t291 * t299 / 0.2e1 + t535 * t1224 + t199 * t1201 + t86 * t1244 + t87 * t1243 + t1366;
t11 = Ifges(7,3) * t947 + t1313 + t842;
t40 = t1202 * t535 + t290 * t645 - t291 * t646 + t857;
t912 = t11 * qJD(1) + t40 * qJD(2);
t888 = t1054 / 0.2e1 - t1043 / 0.2e1;
t889 = t1056 / 0.2e1 - t1045 / 0.2e1;
t903 = t596 * t830 - t597 * t834;
t840 = (-t564 / 0.2e1 + t889) * t831 + (-t786 / 0.2e1 + t888) * t737 + (-pkin(9) * t1013 - t431 * t831 - t782 * t737) * t1302 + (-t1059 * t794 + t737 * t903 + t831 * t909) * t1300 + (-t1059 * t735 + t1337 * t291 + t290 * t882 - t87 * t711 + t86 * t714) * t1298 + t882 * t1242 + t1337 * t1244 + t298 * t1238 + t299 * t1231;
t848 = m(5) * t1272 + t907 * t1300 + (t1088 - t1090) * t1298 + t396 * t1222 + t395 * t1218 + t570 * t1204 + t571 * t1197;
t15 = (-t832 * mrSges(5,2) / 0.2e1 + t942) * t828 + t840 - t783 - t848;
t67 = -t711 * t645 + t714 * t646 + m(7) * (t290 * t714 - t291 * t711) + (m(6) * t903 - t1043 + t1054 + t1326) * t831;
t911 = t15 * qJD(1) + t67 * qJD(2);
t890 = t1217 * t299 + t1222 * t298;
t853 = (t1217 * t377 + t1327) * mrSges(7,3) + t890;
t25 = t853 - t1334;
t892 = mrSges(7,1) * t1231 + mrSges(7,2) * t1237;
t59 = (t1217 * t713 + t1222 * t712) * mrSges(7,3) + t892 + t1316;
t910 = -t25 * qJD(1) + t59 * qJD(2);
t908 = t213 * t834 + t214 * t830;
t906 = t258 * t834 - t616 * t819;
t901 = t619 * t711 - t714 * t956;
t899 = -t788 / 0.2e1 - t610 / 0.2e1 + pkin(3) * mrSges(5,1);
t881 = t87 * t997 + t1306;
t847 = (t957 * t1249 + t1358 - t80 / 0.2e1) * mrSges(7,3) + t957 * t1257 + t1357 + t298 * t1248 + t299 * t1246 + t1368 * t1221 + t1354 * t1216 + t1353 * t1209 + t881;
t13 = t847 - t915;
t845 = (t999 + t962) * t712 + (t998 + t963) * t713 + t619 * t1243 + t1315 + t1333;
t37 = t845 + t1314;
t78 = t819 * t607 - (t611 / 0.2e1 + t1258) * t900 + (t1260 - t613 / 0.2e1) * t764;
t876 = t13 * qJD(1) + t37 * qJD(2) + t78 * qJD(3);
t852 = (-t1218 * t377 + t1327) * mrSges(7,3) + (t1374 * t900 - t764 * t93 - t80) * t1298 + t890;
t20 = (t1275 - t1119 / 0.2e1 - t1108 / 0.2e1) * t834 + (t1117 / 0.2e1 + t1274 - t1109 / 0.2e1) * t830 + t852 + t1381;
t869 = m(7) * (-t1336 * t764 - t1373 * t900);
t873 = (-t711 * t829 + t714 * t833) * t1029;
t45 = -t869 / 0.2e1 + t873 + (t1101 / 0.2e1 - t770 / 0.2e1 + t835 * t992) * t834 + (t1102 / 0.2e1 + t1215 + mrSges(6,3) * t978) * t830 + t59;
t875 = -t20 * qJD(1) + t45 * qJD(2);
t807 = -0.2e1 * t1012;
t843 = (t1079 / 0.2e1 - t1078 / 0.2e1) * mrSges(7,3) + (t1262 - t787 / 0.2e1) * t1059 - t809 + (t807 + t488) * t1302 + (t389 + t807) * t1300 + (-t1184 + t882 * t956 + t1337 * t619 + (-qJ(4) - t819) * t1059 + t389) * t1298 - t1127 / 0.2e1 + t1126 / 0.2e1 - t1120 / 0.2e1 - t1118 / 0.2e1;
t849 = -m(5) * t488 / 0.2e1 + t908 * t1301 - t1321 * t1299 + t381 * t1220 + t380 * t1217 + t547 * t1206 + t546 * t1198;
t29 = t843 + t849;
t872 = -t787 + t959;
t430 = t1186 + mrSges(5,3) + (m(6) + m(5)) * qJ(4) - t872;
t860 = -t764 * t905 + t900 * t904;
t851 = t1299 * t860 + t1301 * t902;
t856 = m(6) * t1210 + (-t901 + t735) * t1298 - t1112 / 0.2e1 - t1111 / 0.2e1;
t55 = (t1106 / 0.2e1 + t1107 / 0.2e1) * t835 + t851 + t856;
t874 = -qJD(1) * t29 - qJD(2) * t55 - qJD(3) * t430;
t867 = (t87 + t93) * t956 + t1374 * t619;
t60 = t712 * t997 + t713 * t996 - t1316 + t892;
t854 = qJ(4) * t1223 + t1210 * t785 - t1242 * t619 + t1315;
t865 = t1336 * t956 - t1373 * t619;
t22 = -pkin(5) * t959 * t978 + t1336 * t997 + t1330 * mrSges(7,3) + (t792 - t925) * t976 - t930 * t977 - t1332 * t1294 - t1382 + t306 * t995 + t290 * t996 + t534 * t1022 - t1055 / 0.4e1 + ((-t1052 * t819 + t735 * t834) * pkin(5) + t865) * t1298 + t612 * t1232 + t613 * t1234 + t501 * t1219 + t537 * t1221 + t1052 * t1212 + t834 * t970 + t830 * t971 + t922 * t1203 - t744 * t1205 - t745 * t1196 + t854 - t1044 / 0.4e1 + (t536 + t502) * t1216 + (t614 + t611) * t1235;
t56 = qJ(4) * t785 + (-t792 / 0.2e1 + t925 / 0.2e1) * t830 + t78 - t1351 + (-t959 + t1186) * t1182;
t844 = -qJ(4) * t939 / 0.2e1 + t785 * t1349 - t1357 - t957 * t614 / 0.4e1 + t790 * t1252 - t930 * t1250 + t299 * t1364 + t298 * t1249 + t1368 * t1219 + t900 * t1376 - t819 * t1353 / 0.2e1 + t1319 * t615;
t861 = -t1140 / 0.2e1 + t213 * mrSges(6,1) / 0.2e1 + t214 * t1296 + t915;
t862 = t381 * t1199 + (t103 * t833 + t104 * t829) * t1298 + t380 * t1207;
t7 = (t198 * t1198 - t1073 / 0.2e1 + t906 * t1299 + t862) * pkin(5) + t867 * t1299 + t861 + t844 + (0.3e1 / 0.4e1 * t1159 - t1037 / 0.2e1 + (Ifges(6,2) / 0.4e1 + t1097 / 0.2e1) * t616 - t967) * t830 + (-t1015 * t900 + t1016 * t764 + t1248 * t957 - t1358) * mrSges(7,3) + (t465 * t1191 + 0.3e1 / 0.4e1 * t1147 + (-Ifges(6,1) / 0.4e1 + t991) * t615 + t944) * t834 - t1306;
t866 = -t7 * qJD(1) + t22 * qJD(2) + t56 * qJD(3);
t117 = (t1249 + t1248) * mrSges(7,2) + (t1246 + t1364) * mrSges(7,1);
t855 = (t299 * t1208 + t298 * t1199 + (t1200 * t957 + t1208 * t377) * mrSges(7,3)) * pkin(5) - t966;
t18 = mrSges(7,1) * t1016 - mrSges(7,2) * t1015 + t855 + t966;
t858 = (t645 * t1199 + t646 * t1208 + (t1058 / 0.2e1 - t1048 / 0.2e1) * mrSges(7,3)) * pkin(5);
t51 = -mrSges(7,1) * t968 + mrSges(7,2) * t969 + t858;
t777 = (mrSges(7,1) * t829 + mrSges(7,2) * t833) * pkin(5);
t859 = -t18 * qJD(1) - t51 * qJD(2) - t117 * qJD(3) + t777 * qJD(5);
t763 = t777 * qJD(6);
t53 = t934 * t1218 + t936 * t1222 + t896 * t1204 + t897 * t1197 + (t1066 / 0.2e1 + t1062 / 0.2e1) * mrSges(7,3) + (mrSges(5,1) - t1103 / 0.2e1 + t1100 / 0.2e1 + m(5) * pkin(9)) * t835 - t851 + t856;
t49 = -t1130 / 0.2e1 - t1131 / 0.2e1 + t1129 / 0.2e1 - t1128 / 0.2e1 + t858 + t535;
t46 = t869 / 0.2e1 + t1329 * t831 + t873 + t60 - t888 + t1332;
t36 = t845 - t1314;
t27 = t843 - t849 - t1110;
t24 = t853 + t1334;
t21 = ((mrSges(6,3) * t975 + t1319) * t834 + (t930 / 0.4e1 + t1212 + t953 + (t1262 - t1186 / 0.2e1) * pkin(5)) * t830) * t835 + (-t1146 / 0.4e1 - t745 / 0.4e1 - t699 / 0.4e1 + t970 + (t1267 + t1187 / 0.2e1) * pkin(5)) * t834 + t865 * t1298 + (t764 * t968 - t900 * t969 + t1330) * mrSges(7,3) + (-t1157 / 0.4e1 + t744 / 0.4e1 - t700 / 0.4e1 + t971) * t830 + t963 * t713 + t962 * t712 + t854 + t1333 + t1382;
t19 = (t1204 * t616 + t985) * mrSges(6,3) + t852 - t889 + t1329 * t737 - t1381;
t17 = -t1178 / 0.2e1 - t1180 / 0.2e1 + t1176 / 0.2e1 - t1175 / 0.2e1 + t855 - t966;
t16 = mrSges(5,2) * t981 + t828 * t942 + t840 + t848;
t12 = t847 + t915;
t10 = t842 + t1311;
t8 = t867 * t1298 + t931 * t1196 - t1046 / 0.4e1 + t922 * t1226 + t958 * t1205 - t1057 / 0.4e1 + t861 + t893 * t737 - t844 - t985 * t1097 + t616 * t953 + t465 * t975 - t1037 * t1206 + t198 * t1022 + t94 * t995 + t86 * t996 + t93 * t997 - t377 * t998 + t957 * t999 + t881 + (t906 * t1298 + t862 + t1073 / 0.2e1) * pkin(5);
t4 = (t1087 + t1089) * t1029 + (-pkin(5) * t1072 + t868) * t1298 + t850 * t835 + t839 + Ifges(6,3) * t946 + t396 * t1023 + t395 * t1024 + pkin(5) * t984 - t863 + t1311;
t1 = (t419 / 0.4e1 + t972 + mrSges(6,3) * t1286) * t834 + (-t418 / 0.4e1 + t1325) * t830 + t838 + ((Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t832 + (-t1018 * t831 + (t788 / 0.4e1 + t610 / 0.4e1 + t1019) * t835) * t836) * t828 + t1307;
t26 = [qJD(2) * t5 + qJD(3) * t6 - qJD(4) * t23 - qJD(5) * t9 + qJD(6) * t14, t1 * qJD(3) + t16 * qJD(4) + t4 * qJD(5) + t10 * qJD(6) + t1122 + (((t1018 * t835 + t1019 * t831 - Ifges(3,6)) * t832 + (t917 * t1202 - t831 * t926 / 0.2e1 - t835 * t919 / 0.2e1 + t1317 * t1194) * t836) * t828 + t794 * t489 + t565 * t786 + t782 * t680 + t225 * t769 + t226 * t770 + t432 * t741 + t735 * t253 - t713 * t229 / 0.2e1 + t712 * t228 / 0.2e1 - pkin(2) * t681 + t111 * t645 + t110 * t646 + t597 * t570 + t596 * t571 + t327 * t534 + t291 * t395 + t290 * t396 + 0.2e1 * (t565 * t782 + (-t493 * t835 + t495 * t831) * pkin(9)) * t1302 + 0.2e1 * (t110 * t290 + t111 * t291 + t327 * t735) * t1298 + 0.2e1 * (t225 * t596 + t226 * t597 + t432 * t794) * t1300 + t456 * t1269 + t455 * t1271 + t671 * t1240 + t700 * t1241 + t810 + m(4) * (-pkin(2) * t740 + (-t553 * t831 + t554 * t835) * pkin(9)) + (t418 * t1198 + t419 * t1206 + t554 * mrSges(4,3) - t493 * mrSges(5,1) + (t706 - t708) * pkin(9) - t1308) * t835 + (t1339 + t624 / 0.2e1 - t626 / 0.2e1 + t417 / 0.2e1 + t227 / 0.2e1 - t553 * mrSges(4,3) + t495 * mrSges(5,1) + (-t707 + t709) * pkin(9)) * t831 + t1312) * qJD(2), t1121 + t1 * qJD(2) + t27 * qJD(4) + t8 * qJD(5) + t12 * qJD(6) + ((qJ(4) * t388 - t1294 * t908) * t1300 + (t103 * t956 + t104 * t619 + t303 * t819) * t1298 + (-pkin(3) * t488 - qJ(4) * t487) * t1302) * t1303 + (t819 * t275 + t388 * t787 + t240 * t1217 + t241 * t1222 + t619 * t380 + t956 * t381 - t303 * t959 - t1337 * t1260 + t882 * t1258 - qJ(4) * t518 - mrSges(5,1) * t1071 + t899 * t736 + (t409 / 0.2e1 - t1294 * t546 + t790 * t1225 - t213 * mrSges(6,3)) * t834 + (-t408 / 0.2e1 - t1294 * t547 + t792 * t1225 - t214 * mrSges(6,3)) * t830 + t1347 * t488 + t1346 * t487 + t1321 * mrSges(7,3) + t1323) * qJD(3), -t1083 + t16 * qJD(2) + t27 * qJD(3) + (t1078 - t1079) * t1173 + t19 * qJD(5) + t24 * qJD(6), -t1096 + t4 * qJD(2) + t8 * qJD(3) + t19 * qJD(4) + (-t159 * mrSges(6,1) - t158 * mrSges(6,2) - t1175 + t1176 + t1322) * qJD(5) + t17 * qJD(6) + (m(7) * (t829 * t94 + t833 * t93) + (-t377 * t829 - t833 * t957) * mrSges(7,3)) * t1136, t1084 + t10 * qJD(2) + t12 * qJD(3) + t24 * qJD(4) + t17 * qJD(5) + (-t1178 + t199 - t1180) * qJD(6); -qJD(3) * t2 + qJD(4) * t15 + qJD(5) * t3 + qJD(6) * t11 - t1122, -qJD(3) * t30 + qJD(4) * t67 + qJD(5) * t39 + qJD(6) * t40, t53 * qJD(4) + t21 * qJD(5) + t36 * qJD(6) + ((-qJ(4) * t793 - t1294 * t902) * t1300 + (t619 * t904 + t819 * t734 + t905 * t956) * t1298) * t1303 + t914 + (t1217 * t923 + t1222 * t928 + t1231 * t614 + t1237 * t612 - t734 * t959 - t793 * t787 + t819 * t937 + (-t860 + t901) * mrSges(7,3) + (-t1141 / 0.2e1 - t1153 / 0.2e1 + t956 * mrSges(7,1) - t619 * mrSges(7,2) + (Ifges(6,5) / 0.2e1 - t1294 * mrSges(6,1)) * t834 + (-Ifges(6,6) / 0.2e1 + t1294 * mrSges(6,2)) * t830 + (-m(5) * pkin(3) + t1347) * pkin(9) - t899 - t1375) * t835 + (t925 * t1206 + t792 * t1204 + (-mrSges(5,1) - t785) * qJ(4) + (-m(5) * qJ(4) + t1346) * pkin(9) + t1343 + t1351) * t831 - t902 * mrSges(6,3)) * qJD(3), t53 * qJD(3) + (-t1062 - t1066) * t1173 + t46 * qJD(5) + t60 * qJD(6) + t911, t21 * qJD(3) + t46 * qJD(4) + (-t597 * mrSges(6,1) - t596 * mrSges(6,2) - t1128 + t1129 + t535 + t743) * qJD(5) + t49 * qJD(6) + (m(7) * (t305 * t833 + t306 * t829) + (-t1048 + t1058) * mrSges(7,3)) * t1136 + t913, t36 * qJD(3) + t60 * qJD(4) + t49 * qJD(5) + (-t1130 + t535 - t1131) * qJD(6) + t912; qJD(2) * t2 + qJD(4) * t29 - qJD(5) * t7 + qJD(6) * t13 - t1121, qJD(4) * t55 + qJD(5) * t22 + qJD(6) * t37 - t914, qJD(4) * t430 + qJD(5) * t56 + qJD(6) * t78, -t874 (mrSges(6,1) * t1051 + mrSges(6,2) * t1039 + t114 - t922) * qJD(5) + t1371 + (m(7) * (-t619 * t833 + t829 * t956) + t1355 * mrSges(7,3)) * t1136 + t866, t114 * qJD(5) + t1371 + t876; -qJD(2) * t15 - qJD(3) * t29 + qJD(5) * t20 + qJD(6) * t25 + t1083, -qJD(3) * t55 - qJD(5) * t45 - qJD(6) * t59 - t911, t874, 0 (-t1297 * t1355 + t872) * qJD(5) + t1372 - t875, qJD(5) * t959 + t1372 - t910; -qJD(2) * t3 + qJD(3) * t7 - qJD(4) * t20 + qJD(6) * t18 + t1096, -qJD(3) * t22 + qJD(4) * t45 + qJD(6) * t51 - t913, qJD(6) * t117 - t866, t875, -t763, -t763 - t859; -qJD(2) * t11 - qJD(3) * t13 - qJD(4) * t25 - qJD(5) * t18 - t1084, -qJD(3) * t37 + qJD(4) * t59 - qJD(5) * t51 - t912, -qJD(5) * t117 - t876, t910, t859, 0;];
Cq  = t26;
