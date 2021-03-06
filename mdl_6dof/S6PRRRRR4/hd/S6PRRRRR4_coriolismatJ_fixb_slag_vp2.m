% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:55:14
% EndTime: 2019-03-09 00:56:18
% DurationCPUTime: 41.81s
% Computational Cost: add. (59693->1198), mult. (152547->1660), div. (0->0), fcn. (174937->14), ass. (0->635)
t742 = cos(qJ(6));
t1127 = -t742 / 0.2e1;
t738 = sin(qJ(6));
t1131 = -t738 / 0.2e1;
t1100 = Ifges(7,4) * t742;
t712 = Ifges(7,1) * t738 + t1100;
t970 = t742 * t712;
t1101 = Ifges(7,4) * t738;
t710 = Ifges(7,2) * t742 + t1101;
t979 = t738 * t710;
t1232 = -t979 / 0.2e1 + t970 / 0.2e1;
t1095 = Ifges(7,2) * t738;
t860 = -t1095 + t1100;
t1107 = Ifges(7,1) * t742;
t862 = -t1101 + t1107;
t769 = t1127 * t860 + t1131 * t862 - t1232;
t1120 = cos(qJ(5));
t960 = t1120 * pkin(4);
t727 = -t960 - pkin(5);
t1057 = t742 * mrSges(7,2);
t1059 = t738 * mrSges(7,1);
t866 = t1057 + t1059;
t817 = t727 * t866;
t828 = pkin(5) * t866;
t1339 = t769 + t828 / 0.2e1 - t817 / 0.2e1;
t1338 = qJD(4) + qJD(5);
t1051 = cos(pkin(6));
t1121 = cos(qJ(2));
t737 = sin(pkin(7));
t1049 = sin(pkin(6));
t1050 = cos(pkin(7));
t869 = t1050 * t1049;
t1231 = t1051 * t737 + t1121 * t869;
t740 = sin(qJ(3));
t744 = cos(qJ(3));
t741 = sin(qJ(2));
t894 = t741 * t1049;
t612 = -t1231 * t744 + t740 * t894;
t743 = cos(qJ(4));
t1119 = sin(qJ(5));
t739 = sin(qJ(4));
t926 = t1119 * t739;
t809 = t1120 * t743 - t926;
t399 = t809 * t612;
t613 = t1231 * t740 + t744 * t894;
t297 = -t399 * t742 + t613 * t738;
t1028 = t297 * t742;
t296 = t399 * t738 + t613 * t742;
t1029 = t296 * t738;
t707 = -mrSges(7,1) * t742 + mrSges(7,2) * t738;
t1136 = t707 / 0.2e1;
t1312 = t1136 - mrSges(6,1) / 0.2e1;
t928 = t1120 * t739;
t810 = t1119 * t743 + t928;
t398 = t810 * t612;
t1337 = t399 * mrSges(6,2) / 0.2e1 + (t1028 / 0.2e1 - t1029 / 0.2e1) * mrSges(7,3) - t1312 * t398;
t1058 = t738 * mrSges(7,3);
t936 = t1058 / 0.2e1;
t937 = -t1058 / 0.2e1;
t1279 = t936 + t937;
t870 = t1121 * t1049;
t786 = t1050 * t1051 - t737 * t870;
t470 = t613 * t743 + t739 * t786;
t776 = -t613 * t739 + t743 * t786;
t765 = t1119 * t776 + t1120 * t470;
t187 = t612 * t738 + t742 * t765;
t276 = t1119 * t470 - t1120 * t776;
t808 = t276 * (t1057 / 0.2e1 + t1059 / 0.2e1);
t819 = t276 * t866;
t1318 = t808 + t819 / 0.2e1 + t1279 * t187;
t1043 = t187 * t742;
t186 = t612 * t742 - t738 * t765;
t1045 = t186 * t738;
t1326 = m(7) * (t765 - t1043 + t1045) * t276;
t1334 = t1326 * qJD(1);
t1336 = t1318 * qJD(6) + t1334;
t1202 = m(7) / 0.2e1;
t841 = t1028 - t1029;
t1335 = (pkin(5) * t398 + pkin(12) * t841) * t1202 + t1337;
t595 = t866 * t809;
t1333 = t595 / 0.2e1;
t989 = t737 * t740;
t684 = t1050 * t739 + t743 * t989;
t807 = -t1050 * t743 + t739 * t989;
t567 = t1119 * t684 + t1120 * t807;
t780 = -t1119 * t807 + t1120 * t684;
t1267 = mrSges(6,1) * t780 - mrSges(6,2) * t567;
t1332 = t1267 / 0.2e1;
t1126 = -t742 / 0.4e1;
t932 = pkin(2) * t1050;
t988 = t737 * t744;
t688 = pkin(9) * t988 + t740 * t932;
t663 = pkin(10) * t1050 + t688;
t664 = (-pkin(3) * t744 - pkin(10) * t740 - pkin(2)) * t737;
t514 = -t739 * t663 + t743 * t664;
t444 = -t684 * pkin(11) + t514;
t406 = -pkin(4) * t988 + t444;
t515 = t663 * t743 + t664 * t739;
t445 = -pkin(11) * t807 + t515;
t927 = t1119 * t445;
t224 = t1120 * t406 - t927;
t207 = pkin(5) * t988 - t224;
t1187 = -t207 / 0.2e1;
t1092 = Ifges(7,6) * t738;
t729 = Ifges(7,5) * t742;
t1223 = t729 - t1092;
t1277 = t567 * t1223;
t490 = -t738 * t780 - t742 * t988;
t1297 = -t490 / 0.4e1;
t1098 = Ifges(7,5) * t567;
t476 = Ifges(7,4) * t490;
t491 = -t738 * t988 + t742 * t780;
t211 = t491 * Ifges(7,1) + t1098 + t476;
t977 = t742 * t211;
t1072 = t567 * Ifges(7,6);
t1102 = Ifges(7,4) * t491;
t210 = Ifges(7,2) * t490 + t1072 + t1102;
t987 = t738 * t210;
t1307 = t977 / 0.4e1 - t987 / 0.4e1;
t1293 = t710 / 0.4e1;
t825 = t1293 * t491 + t1297 * t712;
t890 = -Ifges(7,2) * t491 + t476;
t1330 = t860 * t1297 + t866 * t1187 + t890 * t1126 - t1277 / 0.4e1 + t825 - t1307;
t1193 = -pkin(11) - pkin(10);
t732 = t743 * pkin(10);
t964 = t743 * pkin(11) + t732;
t1219 = t1120 * t964 + t1193 * t926;
t1278 = t1219 * t707;
t1284 = t1219 * mrSges(6,1);
t644 = t1119 * t964 - t1193 * t928;
t1322 = t644 * mrSges(6,2);
t1125 = t742 / 0.2e1;
t1130 = t738 / 0.2e1;
t1137 = t810 / 0.2e1;
t502 = Ifges(7,6) * t810 + t809 * t860;
t504 = Ifges(7,5) * t810 + t809 * t862;
t616 = Ifges(6,5) * t809 - Ifges(6,6) * t810;
t709 = Ifges(7,5) * t738 + Ifges(7,6) * t742;
t812 = t502 * t1125 + t504 * t1130 + t709 * t1137 + t1232 * t809 + t616;
t1329 = t1278 + t812 - t1284 + t1322;
t1328 = -t1278 / 0.2e1 + t1284 / 0.2e1 - t1322 / 0.2e1;
t1290 = t810 / 0.4e1;
t1291 = t809 / 0.4e1;
t1294 = -t567 / 0.4e1;
t1296 = t491 / 0.4e1;
t209 = Ifges(7,5) * t491 + Ifges(7,6) * t490 + Ifges(7,3) * t567;
t1103 = Ifges(6,4) * t780;
t356 = -Ifges(6,2) * t567 - Ifges(6,6) * t988 + t1103;
t549 = Ifges(6,4) * t567;
t357 = Ifges(6,1) * t780 - Ifges(6,5) * t988 - t549;
t1088 = Ifges(7,3) * t810;
t500 = t1223 * t809 + t1088;
t686 = -pkin(9) * t989 + t744 * t932;
t662 = -pkin(3) * t1050 - t686;
t582 = pkin(4) * t807 + t662;
t998 = t809 * t738;
t605 = -mrSges(7,2) * t810 - mrSges(7,3) * t998;
t1056 = t742 * mrSges(7,3);
t607 = mrSges(7,1) * t810 - t1056 * t809;
t614 = mrSges(6,1) * t810 + mrSges(6,2) * t809;
t691 = Ifges(6,4) * t809;
t617 = -Ifges(6,2) * t810 + t691;
t1060 = t810 * Ifges(6,4);
t619 = Ifges(6,1) * t809 - t1060;
t929 = t1120 * t445;
t225 = t1119 * t406 + t929;
t208 = -pkin(12) * t988 + t225;
t272 = t567 * pkin(5) - pkin(12) * t780 + t582;
t97 = -t208 * t738 + t272 * t742;
t98 = t208 * t742 + t272 * t738;
t1327 = t1307 * t809 + t582 * t614 / 0.2e1 + t780 * t619 / 0.4e1 + t617 * t1294 + t567 * t500 / 0.4e1 + t504 * t1296 + t490 * t502 / 0.4e1 + t207 * t1333 + t98 * t605 / 0.2e1 + t97 * t607 / 0.2e1 - t616 * t988 / 0.4e1 - t810 * t356 / 0.4e1 + t209 * t1290 + t357 * t1291;
t1239 = t780 * t709;
t1265 = Ifges(7,6) * t780 - t567 * t860;
t1315 = t742 * t1265;
t1264 = Ifges(7,5) * t780 - t567 * t862;
t1316 = t738 * t1264;
t1325 = t1239 / 0.2e1 + t1315 / 0.2e1 + t1316 / 0.2e1;
t1247 = Ifges(6,6) * t780;
t1287 = Ifges(6,5) * t567;
t1324 = -t1239 / 0.4e1 + t1247 / 0.2e1 - (t979 / 0.4e1 - t970 / 0.4e1) * t567 + t1287 / 0.2e1 - t1315 / 0.4e1 - t1316 / 0.4e1;
t1189 = -t186 / 0.2e1;
t1323 = -t187 / 0.2e1;
t1151 = -t612 / 0.2e1;
t1282 = t866 * t567;
t1178 = -t1282 / 0.2e1;
t1321 = t225 * t644;
t1320 = t644 * t738;
t1319 = t644 * t742;
t1242 = t644 * t765;
t1317 = t1119 * t644;
t1002 = t1219 * t644;
t1314 = -Ifges(6,2) * t780 - t549;
t1146 = t644 / 0.2e1;
t1313 = t1282 * t1146;
t1240 = t765 * t707;
t1244 = t765 * mrSges(6,1);
t1285 = t276 * mrSges(6,2);
t1311 = t1240 - t1244 + t1285;
t967 = -t1287 - t1247;
t1221 = t614 * t1151 + t1189 * t607 + t605 * t1323;
t1310 = t1333 * t276 - t1221;
t1308 = -t1240 / 0.2e1 + t1244 / 0.2e1 - t1285 / 0.2e1;
t733 = t738 ^ 2;
t735 = t742 ^ 2;
t963 = t733 + t735;
t1281 = t963 * t276;
t1304 = -pkin(5) * t765 - pkin(12) * t1281;
t728 = -pkin(4) * t743 - pkin(3);
t599 = -pkin(5) * t809 - pkin(12) * t810 + t728;
t428 = -t1219 * t738 + t599 * t742;
t429 = t1219 * t742 + t599 * t738;
t501 = -Ifges(7,3) * t809 + t1223 * t810;
t596 = t866 * t810;
t618 = Ifges(6,2) * t809 + t1060;
t620 = Ifges(6,1) * t810 + t691;
t1063 = t809 * Ifges(7,5);
t505 = t810 * t862 - t1063;
t972 = t742 * t505;
t1062 = t809 * Ifges(7,6);
t503 = t810 * t860 - t1062;
t981 = t738 * t503;
t1303 = -(-t972 / 0.2e1 + t981 / 0.2e1 + t500 / 0.2e1 - t617 / 0.2e1 - t620 / 0.2e1) * t809 + t1219 * t596 - (t504 * t1127 + t502 * t1130 - t501 / 0.2e1 + t618 / 0.2e1 - t619 / 0.2e1) * t810 + t428 * t607 + t429 * t605 + t644 * t595 + t728 * t614;
t1071 = t780 * mrSges(6,3);
t1262 = -mrSges(7,2) * t780 + t1058 * t567;
t1263 = mrSges(7,1) * t780 + t1056 * t567;
t307 = -mrSges(7,1) * t490 + mrSges(7,2) * t491;
t1183 = t307 / 0.2e1;
t507 = -mrSges(6,1) * t988 - t1071;
t895 = -t507 / 0.2e1 + t1183;
t1302 = t1267 * t1151 + t1263 * t1189 + t1262 * t1323 + (t1071 / 0.2e1 - t895) * t765;
t959 = t1119 * pkin(4);
t726 = t959 + pkin(12);
t1301 = -t726 * t1281 + t727 * t765;
t1145 = -t1219 / 0.2e1;
t1176 = t428 / 0.2e1;
t1300 = t728 * t1332 + t1262 * t429 / 0.2e1 + t1263 * t1176 + (mrSges(6,3) * t1145 - t618 / 0.4e1 + t501 / 0.4e1) * t780;
t1169 = t491 / 0.2e1;
t1171 = t490 / 0.2e1;
t1299 = t1264 * t1169 + t1265 * t1171 + t98 * t1262 + t97 * t1263 + t582 * t1267 - t207 * t1282;
t1298 = -mrSges(6,2) / 0.2e1;
t1295 = -t567 / 0.2e1;
t1292 = t712 / 0.4e1;
t1245 = Ifges(7,3) * t780;
t1288 = t1245 / 0.2e1;
t1286 = pkin(5) * t1219;
t946 = t729 / 0.2e1;
t1283 = (t946 - t1092 / 0.2e1) * t567;
t1276 = t727 * t1219;
t1204 = m(6) / 0.2e1;
t1275 = t728 * t1204;
t1065 = t684 * mrSges(5,3);
t642 = -mrSges(5,1) * t988 - t1065;
t1273 = t642 + t1065;
t232 = t1120 * t444 - t927;
t1118 = pkin(4) * t684;
t382 = pkin(5) * t780 + pkin(12) * t567;
t329 = t1118 + t382;
t110 = -t232 * t738 + t329 * t742;
t111 = t232 * t742 + t329 * t738;
t846 = -t110 * t738 + t111 * t742;
t995 = t810 * t738;
t958 = mrSges(7,3) * t995;
t606 = mrSges(7,2) * t809 - t958;
t994 = t810 * t742;
t608 = -mrSges(7,1) * t809 - mrSges(7,3) * t994;
t1229 = t1127 * t606 + t1130 * t608;
t341 = -mrSges(7,2) * t567 + mrSges(7,3) * t490;
t342 = mrSges(7,1) * t567 - mrSges(7,3) * t491;
t1270 = t1125 * t342 + t1130 * t341;
t231 = t1119 * t444 + t929;
t1040 = t231 * t276;
t1140 = -t684 / 0.2e1;
t1164 = -t780 / 0.2e1;
t1188 = t187 / 0.2e1;
t1243 = t207 * t765;
t1257 = t776 / 0.2e1;
t799 = t807 * mrSges(5,3);
t641 = mrSges(5,2) * t988 - t799;
t1073 = t567 * mrSges(6,3);
t506 = mrSges(6,2) * t988 - t1073;
t1166 = t506 / 0.2e1;
t974 = t742 * t341;
t983 = t738 * t342;
t801 = t974 / 0.2e1 - t983 / 0.2e1 + t1166;
t1052 = t98 * t742;
t1054 = t97 * t738;
t856 = t1052 - t1054;
t1268 = (t1118 * t612 + t1040) * t1204 + (t110 * t186 + t111 * t187 + t1040 + t1243) * t1202 + t186 * t1263 / 0.2e1 + t1262 * t1188 - t470 * t642 / 0.2e1 + t641 * t1257 + (t470 * t1140 + t807 * t1257) * mrSges(5,3) + ((-t224 + t232) * t1204 + t1164 * mrSges(6,3) + t895) * t765 + (mrSges(6,3) * t1295 - t1202 * t856 - t1204 * t225 + t1178 - t801) * t276;
t1266 = -Ifges(6,1) * t567 - t1103;
t1261 = t1245 - t1277;
t1258 = Ifges(6,6) / 0.2e1;
t1256 = t981 / 0.4e1;
t1255 = t987 / 0.2e1;
t911 = -t988 / 0.2e1;
t1235 = t307 - t507;
t592 = mrSges(5,1) * t807 + t684 * mrSges(5,2);
t957 = mrSges(4,3) * t989;
t696 = mrSges(4,1) * t1050 - t957;
t1234 = t592 - t696;
t615 = -mrSges(6,1) * t809 + mrSges(6,2) * t810;
t708 = -mrSges(5,1) * t743 + mrSges(5,2) * t739;
t1233 = t615 + t708;
t1230 = t890 + t211;
t655 = t809 * t988;
t587 = -t655 * t738 + t742 * t989;
t654 = t810 * t988;
t442 = -mrSges(7,2) * t654 + mrSges(7,3) * t587;
t588 = t655 * t742 + t738 * t989;
t443 = mrSges(7,1) * t654 - mrSges(7,3) * t588;
t1228 = t442 * t1125 + t443 * t1131;
t1225 = -Ifges(5,5) * t807 - Ifges(5,6) * t684 + t967;
t1155 = t596 / 0.2e1;
t1061 = t810 * mrSges(6,3);
t939 = t1061 / 0.2e1;
t1224 = t1155 + t939;
t730 = Ifges(5,4) * t743;
t891 = -Ifges(5,2) * t739 + t730;
t1117 = pkin(4) * t739;
t621 = pkin(5) * t810 - pkin(12) * t809;
t603 = t621 + t1117;
t436 = t603 * t742 + t1320;
t437 = t603 * t738 - t1319;
t838 = -t436 * t738 + t437 * t742;
t1167 = -t506 / 0.2e1;
t1222 = t1167 + t1178;
t1138 = -t809 / 0.4e1;
t940 = t1063 / 0.2e1;
t1218 = -Ifges(7,6) * t998 / 0.2e1 + t742 * t940 + t1088 / 0.2e1 + t1223 * t1138;
t935 = -t1056 / 0.2e1;
t1217 = t490 * t936 + t491 * t935 - t1270;
t938 = -t1061 / 0.2e1;
t1214 = t765 * t938 + t1310;
t836 = t741 * t869;
t660 = -t740 * t836 + t744 * t870;
t873 = t737 * t894;
t581 = t660 * t743 + t739 * t873;
t800 = t660 * t739 - t743 * t873;
t380 = t1119 * t581 + t1120 * t800;
t381 = -t1119 * t800 + t1120 * t581;
t1213 = t381 * t1298 + t1312 * t380;
t1211 = t1202 * t207 - t1204 * t224 + t895;
t1099 = Ifges(6,5) * t655;
t687 = (pkin(3) * t740 - pkin(10) * t744) * t737;
t579 = -t739 * t686 + t743 * t687;
t968 = t743 * t744;
t473 = (pkin(4) * t740 - pkin(11) * t968) * t737 + t579;
t580 = t743 * t686 + t739 * t687;
t978 = t739 * t744;
t951 = t737 * t978;
t513 = -pkin(11) * t951 + t580;
t304 = -t1119 * t513 + t1120 * t473;
t282 = -pkin(5) * t989 - t304;
t305 = t1119 * t473 + t1120 * t513;
t1209 = t588 * t1292 + t587 * t1293 + t282 * t1136 + t305 * t1298 + t304 * mrSges(6,1) / 0.2e1 + t1099 / 0.2e1;
t1207 = t737 ^ 2;
t1206 = 2 * qJD(3);
t1205 = m(5) / 0.2e1;
t1203 = -m(7) / 0.2e1;
t1201 = m(6) * pkin(4);
t1200 = m(7) * pkin(4);
t1198 = -mrSges(7,1) / 0.2e1;
t1197 = mrSges(7,1) / 0.2e1;
t1196 = -mrSges(5,2) / 0.2e1;
t1195 = -mrSges(7,2) / 0.2e1;
t1194 = Ifges(7,3) / 0.2e1;
t1192 = pkin(12) * mrSges(7,3);
t283 = pkin(12) * t989 + t305;
t647 = pkin(4) * t951 + t688;
t393 = t654 * pkin(5) - t655 * pkin(12) + t647;
t142 = -t283 * t738 + t393 * t742;
t1191 = t142 / 0.2e1;
t143 = t283 * t742 + t393 * t738;
t1190 = -t143 / 0.2e1;
t1186 = t209 / 0.2e1;
t1185 = -t276 / 0.2e1;
t1184 = t276 / 0.2e1;
t1182 = t341 / 0.2e1;
t1181 = -t342 / 0.2e1;
t1180 = t342 / 0.2e1;
t385 = -mrSges(7,1) * t587 + mrSges(7,2) * t588;
t1177 = t385 / 0.2e1;
t1174 = t437 / 0.2e1;
t1173 = t443 / 0.2e1;
t457 = t621 * t738 - t1319;
t1172 = t457 / 0.2e1;
t1170 = -t491 / 0.2e1;
t1162 = t567 / 0.2e1;
t1159 = t780 / 0.2e1;
t1158 = t587 / 0.2e1;
t1157 = t588 / 0.2e1;
t594 = t707 * t810;
t1156 = -t594 / 0.2e1;
t1154 = t606 / 0.2e1;
t1153 = -t608 / 0.2e1;
t1152 = t608 / 0.2e1;
t1150 = t612 / 0.2e1;
t1147 = -t641 / 0.2e1;
t1144 = -t654 / 0.2e1;
t1143 = t654 / 0.2e1;
t1142 = t655 / 0.2e1;
t1139 = t684 / 0.2e1;
t1135 = t709 / 0.4e1;
t1134 = t727 / 0.2e1;
t1129 = t738 / 0.4e1;
t1128 = t739 / 0.2e1;
t1124 = t742 / 0.4e1;
t1123 = t743 / 0.2e1;
t306 = mrSges(7,1) * t491 + mrSges(7,2) * t490;
t1116 = pkin(5) * t306;
t1115 = pkin(5) * t595;
t1113 = pkin(12) * t738;
t1112 = pkin(12) * t742;
t1109 = mrSges(7,3) * t733;
t1108 = mrSges(7,3) * t735;
t1106 = Ifges(4,4) * t740;
t1105 = Ifges(5,4) * t684;
t1104 = Ifges(5,4) * t739;
t1097 = Ifges(7,5) * t588;
t1094 = Ifges(6,6) * t654;
t1093 = Ifges(7,6) * t587;
t1089 = Ifges(7,3) * t654;
t1087 = t224 * mrSges(6,2);
t1086 = t225 * mrSges(6,1);
t1085 = t225 * mrSges(6,3);
t1084 = t231 * mrSges(6,1);
t1083 = t232 * mrSges(6,2);
t1074 = t429 * mrSges(7,3);
t1068 = t644 * mrSges(6,3);
t1064 = t809 * mrSges(6,3);
t1041 = t225 * t707;
t1039 = t231 * t644;
t1038 = t231 * t707;
t1035 = t276 * t380;
t1034 = t276 * t398;
t1032 = t765 * t596;
t1027 = t380 * t644;
t1025 = t398 * t644;
t659 = t740 * t870 + t744 * t836;
t302 = -t381 * t738 + t659 * t742;
t303 = t381 * t742 + t659 * t738;
t487 = t612 * t659;
t782 = t737 * t786;
t42 = m(7) * (t186 * t302 + t187 * t303 + t1035) + m(6) * (t381 * t765 + t1035 + t487) + m(5) * (t470 * t581 - t776 * t800 + t487) + m(4) * (t613 * t660 + t782 * t894 + t487);
t1024 = t42 * qJD(1);
t455 = t612 * t613;
t43 = m(7) * (t186 * t296 + t187 * t297 - t1034) + m(6) * (-t399 * t765 - t1034 + t455) + m(5) * (t455 + (-t470 * t743 + t739 * t776) * t612);
t1023 = t43 * qJD(1);
t1014 = t491 * t742;
t1004 = t612 * t739;
t1003 = t612 * t743;
t993 = t726 * t738;
t992 = t726 * t742;
t991 = t727 * t306;
t990 = t727 * t595;
t982 = t738 * t490;
t969 = t743 * t581;
t308 = Ifges(7,5) * t490 - Ifges(7,6) * t491;
t962 = -t739 ^ 2 - t743 ^ 2;
t961 = t1201 / 0.2e1;
t956 = mrSges(4,3) * t988;
t955 = t1113 / 0.2e1;
t954 = -t1112 / 0.2e1;
t952 = -Ifges(7,1) / 0.4e1 + Ifges(7,2) / 0.4e1;
t949 = -t1109 / 0.2e1;
t948 = -t1108 / 0.2e1;
t947 = t1102 / 0.2e1;
t942 = t1064 / 0.2e1;
t934 = t1056 / 0.2e1;
t931 = t738 * t1120;
t930 = t742 * t1120;
t916 = t1004 / 0.2e1;
t915 = -t993 / 0.2e1;
t914 = t992 / 0.2e1;
t913 = -t989 / 0.2e1;
t912 = t989 / 0.2e1;
t910 = t988 / 0.2e1;
t680 = Ifges(5,4) * t807;
t892 = -Ifges(5,2) * t684 - t680;
t886 = mrSges(6,3) * t960;
t885 = mrSges(6,3) * t959;
t883 = t960 / 0.2e1;
t881 = t959 / 0.2e1;
t880 = t765 * t939;
t876 = t739 * t911;
t875 = t743 * t910;
t874 = -t931 / 0.2e1;
t872 = t1207 * t894;
t871 = mrSges(7,3) * (t735 / 0.2e1 + t733 / 0.2e1);
t868 = mrSges(4,1) * t740 + mrSges(4,2) * t744;
t867 = mrSges(5,1) * t739 + mrSges(5,2) * t743;
t865 = Ifges(5,1) * t743 - t1104;
t863 = Ifges(7,1) * t490 - t1102;
t859 = Ifges(5,5) * t743 - Ifges(5,6) * t739;
t375 = mrSges(6,1) * t567 + mrSges(6,2) * t780;
t528 = -Ifges(5,2) * t807 - Ifges(5,6) * t988 + t1105;
t529 = Ifges(5,1) * t684 - Ifges(5,5) * t988 - t680;
t591 = t684 * mrSges(5,1) - mrSges(5,2) * t807;
t784 = -Ifges(5,1) * t807 - t1105;
t7 = t1299 + t1225 * t911 + t567 * t1255 + (t224 * t567 - t225 * t780) * mrSges(6,3) + t356 * t1164 + t662 * t591 + t514 * t641 + t784 * t1139 + t528 * t1140 + t375 * t1118 + t232 * t506 + t111 * t341 + t110 * t342 + m(6) * (t1118 * t582 - t224 * t231 + t225 * t232) + t1261 * t1162 + t1266 * t1159 + m(7) * (t110 * t97 + t111 * t98 + t207 * t231) + t780 * t1186 + (-t529 / 0.2e1 + t514 * mrSges(5,3) - t892 / 0.2e1) * t807 - t1273 * t515 + t1235 * t231 + (t977 + t357 + t1314) * t1295;
t840 = -t302 * t738 + t303 * t742;
t759 = (t380 * t727 + t726 * t840) * t1202 + t581 * t1196 - t800 * mrSges(5,1) / 0.2e1 + (t1119 * t381 - t1120 * t380) * t961;
t9 = t303 * t935 + t302 * t936 - t759 + (t591 + t1267) * t1150 - t1213 + t1268;
t857 = t9 * qJD(1) + t7 * qJD(2);
t853 = pkin(4) * t874;
t852 = t742 * t883;
t851 = t302 * t937 + t303 * t934 + t1213;
t318 = t1089 + t1093 + t1097;
t319 = Ifges(7,4) * t588 + Ifges(7,2) * t587 + Ifges(7,6) * t654;
t320 = Ifges(7,1) * t588 + Ifges(7,4) * t587 + Ifges(7,5) * t654;
t463 = Ifges(6,4) * t655 - Ifges(6,2) * t654 + Ifges(6,6) * t989;
t464 = Ifges(6,1) * t655 - Ifges(6,4) * t654 + Ifges(6,5) * t989;
t494 = mrSges(6,1) * t654 + mrSges(6,2) * t655;
t601 = -mrSges(6,2) * t989 - mrSges(6,3) * t654;
t602 = mrSges(6,1) * t989 - mrSges(6,3) * t655;
t626 = (Ifges(5,6) * t740 + t744 * t891) * t737;
t627 = (Ifges(5,5) * t740 + t744 * t865) * t737;
t661 = t867 * t988;
t676 = (-mrSges(5,2) * t740 - mrSges(5,3) * t978) * t737;
t677 = (mrSges(5,1) * t740 - mrSges(5,3) * t968) * t737;
t697 = -mrSges(4,2) * t1050 + t956;
t720 = Ifges(4,5) * t988;
t814 = t744 * t859;
t10 = (t697 - t956) * t686 + m(5) * (t514 * t579 + t515 * t580) - t807 * t626 / 0.2e1 + t320 * t1169 + t319 * t1171 + t211 * t1157 + t210 * t1158 + t464 * t1159 + t318 * t1162 + t357 * t1142 + t209 * t1143 + t356 * t1144 + m(6) * (t224 * t304 + t225 * t305 + t582 * t647) + m(7) * (t142 * t97 + t143 * t98 + t207 * t282) + t515 * t676 + t514 * t677 + t662 * t661 + t647 * t375 + t580 * t641 + t579 * t642 + t224 * t602 + t225 * t601 + t582 * t494 + t627 * t1139 + t305 * t506 + t304 * t507 + t98 * t442 + t97 * t443 + t207 * t385 + t143 * t341 + t142 * t342 + t282 * t307 + (Ifges(4,6) * t1050 + (Ifges(4,2) * t744 + t1106) * t737) * t913 + (Ifges(6,3) * t989 - t1094 + t1099) * t911 + t1050 * (-Ifges(4,6) * t989 + t720) / 0.2e1 + t463 * t1295 + t528 * t876 + t529 * t875 + (Ifges(5,5) * t684 + Ifges(6,5) * t780 - Ifges(5,6) * t807 - Ifges(6,6) * t567 + (-Ifges(5,3) - Ifges(6,3)) * t988) * t912 + (Ifges(4,5) * t1050 + 0.2e1 * Ifges(4,4) * t988 + (Ifges(4,1) - Ifges(4,2)) * t989) * t910 + (t740 * (Ifges(4,1) * t744 - t1106) / 0.2e1 - pkin(2) * t868 - t744 * (Ifges(5,3) * t740 + t814) / 0.2e1) * t1207 + (m(5) * t662 + t1234 - t957) * t688;
t746 = t677 * t1257 + (-t225 * t399 - t276 * t304 + t305 * t765 + t582 * t613 + t612 * t647) * t1204 + (t580 * t470 + t579 * t776 + t662 * t613 + (t514 * t739 - t515 * t743 + t688) * t612) * t1205 + (t142 * t186 + t143 * t187 + t276 * t282 + t296 * t97 + t297 * t98) * t1202 + t297 * t1182 + t602 * t1185 + t442 * t1188 + t186 * t1173 + t276 * t1177 + t296 * t1180 - t399 * t1166 + t1003 * t1147 + t697 * t1151 - t613 * t696 / 0.2e1 + t470 * t676 / 0.2e1 + t765 * t601 / 0.2e1 + t868 * t782 / 0.2e1 + t642 * t916 - t1211 * t398 + (t661 + t494) * t1150 + (t592 + t375) * t613 / 0.2e1 + (t612 * t910 + t613 * t913) * mrSges(4,3);
t794 = t739 * t800;
t749 = (t794 + t969) * pkin(10) * t1205 + (t1219 * t381 + t1027) * t1204 + (t302 * t428 + t303 * t429 + t1027) * t1202 + t302 * t1152 + t303 * t1154 - t660 * mrSges(4,2) / 0.2e1 + t381 * t942 + t1224 * t380 + (t969 / 0.2e1 + t794 / 0.2e1) * mrSges(5,3) + (-pkin(3) * t1205 + t1275 - mrSges(4,1) / 0.2e1 + t1233 / 0.2e1) * t659;
t6 = -t746 + t749;
t850 = -t6 * qJD(1) + t10 * qJD(2);
t120 = -t224 * t738 + t382 * t742;
t121 = t224 * t742 + t382 * t738;
t11 = t120 * t342 + m(7) * (t120 * t97 + t121 * t98) + t121 * t341 + t224 * t506 + t967 * t911 + (t1186 - t356 / 0.2e1 + t1266 / 0.2e1 - t1085) * t780 + (t1261 / 0.2e1 - t1314 / 0.2e1 - t357 / 0.2e1 - t977 / 0.2e1 + t1255 + t224 * mrSges(6,3)) * t567 + (m(7) * t207 + t1235) * t225 + t1299;
t774 = (-pkin(5) * t380 + pkin(12) * t840) * t1202 + t851;
t803 = t120 * t186 + t121 * t187 + t1243;
t835 = t225 - t856;
t14 = t803 * t1203 + (t835 * t1203 + t1282 / 0.2e1 + t1073 / 0.2e1 + t801) * t276 + t774 + t1302;
t848 = -t14 * qJD(1) + t11 * qJD(2);
t22 = t207 * t306 + t308 * t1162 + t863 * t1169 - t98 * t342 + t210 * t1170 + t97 * t341 + (-t490 * t97 - t491 * t98) * mrSges(7,3) + t1230 * t1171;
t772 = (t1170 * t187 + t1189 * t490) * mrSges(7,3) + t186 * t1182 + t187 * t1181 + t306 * t1184;
t832 = t1195 * t303 + t1197 * t302;
t27 = t772 - t832;
t847 = t27 * qJD(1) + t22 * qJD(2);
t845 = -t120 * t738 + t121 * t742;
t844 = -t142 * t738 + t143 * t742;
t842 = t1219 * t276 + t1242;
t839 = -t428 * t738 + t429 * t742;
t456 = t621 * t742 + t1320;
t837 = -t456 * t738 + t457 * t742;
t833 = t1195 * t297 + t1197 * t296;
t831 = mrSges(7,2) * t1174 + t1198 * t436;
t830 = mrSges(7,2) * t1172 + t1198 * t456;
t826 = t1219 - t839;
t824 = t1256 - t972 / 0.4e1;
t823 = t1127 * t608 + t1131 * t606;
t711 = t743 * Ifges(5,2) + t1104;
t713 = Ifges(5,1) * t739 + t730;
t821 = t713 * t1123 - t739 * t711 / 0.2e1;
t811 = t963 * t1120;
t754 = -m(6) * (pkin(4) * t1004 - t1242 + t842) / 0.2e1 + (t186 * t436 + t187 * t437 + t842) * t1203 + t867 * t1151 + (-t596 / 0.2e1 + t938) * t765 + (-m(6) * t1145 - t839 * t1203 - t1229) * t276;
t756 = (-t398 * t727 + t726 * t841) * t1202 + (-t1119 * t399 + t1120 * t398) * t961 + mrSges(5,1) * t916 + mrSges(5,2) * t1003 / 0.2e1 + t1337;
t18 = t754 + t880 + t756 - t1310;
t745 = t1300 + (t972 + t620) * t1294 + t567 * t1256 + t1211 * t1219 + t341 * t1174 + t436 * t1180 + t110 * t1152 + t111 * t1154 + (t1219 * t232 + t1039 - t1321) * t1204 - t1313 + (t784 / 0.4e1 + (t1147 - t799 / 0.2e1) * pkin(10) - t528 / 0.4e1 + t582 * pkin(4) * t1204) * t739 - pkin(3) * t591 / 0.2e1 + t232 * t942 + t375 * t1117 / 0.2e1 + t615 * t1118 / 0.2e1 + t1261 * t1138 + t1266 * t1290 + t1068 * t1295 - t644 * t1166 - t1265 * t995 / 0.4e1 + t1264 * t994 / 0.4e1 - t224 * t1064 / 0.2e1 + t1314 * t1291 + t662 * t867 / 0.2e1 + (-t713 / 0.4e1 - t891 / 0.4e1) * t807 + (t110 * t428 + t111 * t429 + t436 * t97 + t437 * t98 + t1039) * t1202 + (t892 + t529) * t743 / 0.4e1 + t1224 * t231 - t1273 * t732 / 0.2e1 + (pkin(4) * t1275 - t711 / 0.4e1 + t865 / 0.4e1) * t684 - t737 * t814 / 0.4e1;
t783 = t142 * t937 + t143 * t934 + t320 * t1129 + t319 * t1124 + t654 * t1135 - t1094 / 0.2e1 + Ifges(6,3) * t912 + t1209;
t750 = (t1119 * t305 + t1120 * t304) * t961 + Ifges(5,3) * t912 + t601 * t881 + t602 * t883 + t282 * t727 * t1202 + t385 * t1134 + t579 * mrSges(5,1) / 0.2e1 + t580 * t1196 + t783 + Ifges(5,5) * t875 + Ifges(5,6) * t876 + (t1202 * t844 + t1228) * t726;
t2 = t225 * t939 - t1327 - t745 + t750;
t36 = -pkin(3) * t867 + t865 * t1128 + t891 * t1123 + t437 * t606 + t436 * t608 + m(7) * (t428 * t436 + t429 * t437 + t1002) + t821 + t1303 + (m(6) * t728 + t615) * t1117;
t806 = -t18 * qJD(1) - t2 * qJD(2) + t36 * qJD(3);
t802 = t186 * t456 + t187 * t457 + t1242;
t20 = -t1032 / 0.2e1 + t802 * t1203 + (t826 * t1203 - t595 / 0.2e1 - t1229) * t276 + t1221 + t1335;
t777 = t225 * t938 + t1327;
t748 = -(t1264 * t1126 + t1265 * t1129 - t1085 / 0.2e1 - t1266 / 0.4e1) * t810 - (t1261 / 0.4e1 - t1314 / 0.4e1) * t809 + (t120 * t428 + t121 * t429 + t1219 * t207 + t456 * t97 + t457 * t98 + t1321) * t1202 + (-t1068 / 0.2e1 - t620 / 0.4e1 + t824) * t567 + t777 + t121 * t1154 + t120 * t1152 + t225 * t1155 + t456 * t1180 + t341 * t1172 + t1300;
t781 = (-pkin(5) * t282 + pkin(12) * t844) * t1203 + pkin(5) * t1177;
t3 = (-t709 / 0.4e1 + t1258) * t654 + t748 + t1222 * t644 + t895 * t1219 + (pkin(12) * t1173 + mrSges(7,3) * t1191 - t320 / 0.4e1) * t738 + (mrSges(7,3) * t1190 - pkin(12) * t442 / 0.2e1 - t319 / 0.4e1) * t742 + Ifges(6,3) * t913 + t781 - t1209;
t41 = t456 * t608 + m(7) * (t428 * t456 + t429 * t457 + t1002) + t457 * t606 + t1303;
t805 = -t20 * qJD(1) + t3 * qJD(2) + t41 * qJD(3);
t597 = t710 * t810;
t598 = t712 * t810;
t751 = (-t428 * mrSges(7,3) / 0.2e1 - t597 / 0.4e1 + t505 / 0.4e1) * t490 + (-t1074 / 0.2e1 - t598 / 0.4e1 - t503 / 0.4e1) * t491 - (t863 * t1126 + t210 * t1124 + t567 * t1135 + (t1052 / 0.2e1 - t1054 / 0.2e1) * mrSges(7,3) + t1230 * t1129) * t810 + t207 * t1156 + t341 * t1176 + t429 * t1181 + t306 * t1146 + t308 * t1138 + t97 * t1154 + t98 * t1153;
t778 = t1097 / 0.2e1 + t1093 / 0.2e1 + t1089 / 0.2e1 + mrSges(7,1) * t1191 + mrSges(7,2) * t1190;
t17 = t751 - t778;
t771 = -(t1043 / 0.2e1 - t1045 / 0.2e1) * t810 * mrSges(7,3) + t186 * t1154 + t187 * t1153 + t276 * t1156;
t33 = t771 - t833;
t60 = t429 * t608 + t644 * t594 - ((-t505 / 0.2e1 + t597 / 0.2e1 + t940) * t738 + (-t598 / 0.2e1 - t503 / 0.2e1 + t1062 / 0.2e1 - t1074) * t742) * t810 + (-t958 - t606) * t428;
t804 = t33 * qJD(1) + t17 * qJD(2) - t60 * qJD(3);
t798 = t121 * t1195 + t120 * t1197 + t1288;
t797 = (-0.3e1 / 0.4e1 * t1092 + t729 / 0.4e1 + t946) * t809;
t747 = (t727 * t225 + t845 * t726 + (t1119 * t207 + t930 * t98 - t931 * t97) * pkin(4)) * t1202 - t1087 / 0.2e1 - t1086 / 0.2e1 + t1041 / 0.2e1 - t1282 * t1134 + t120 * t937 + t121 * t934 + t1263 * t915 + t1262 * t914 + t307 * t881 + t342 * t853 + t341 * t852 + (t506 + t1073) * t883 - (t507 + t1071) * t959 / 0.2e1 - t1324;
t755 = (-pkin(5) * t231 + pkin(12) * t846) * t1203 + pkin(5) * t1178 + t1084 / 0.2e1 - t1038 / 0.2e1 + t1083 / 0.2e1 + t110 * t936 + t111 * t935 + t1263 * t955 + t1262 * t954 + t1324;
t12 = t747 + t755;
t757 = ((t1119 * t276 - t186 * t931 + t187 * t930) * pkin(4) + t1301) * t1202 + t276 * t949 + t276 * t948 - t1308;
t762 = t1304 * t1203 - (t949 + t948) * t276 + t1308;
t30 = t757 + t762;
t764 = (-mrSges(6,1) + t707) * t959 + (mrSges(7,3) * t963 - mrSges(6,2)) * t960;
t383 = (t1119 * t727 + t726 * t811) * t1200 + t764;
t752 = (t1276 + t837 * t726 + (-t428 * t931 + t429 * t930 + t1317) * pkin(4)) * t1202 + t990 / 0.2e1 + t456 * t937 + t457 * t934 + t607 * t915 + t605 * t914 + t596 * t881 + t608 * t853 + t606 * t852 - t1328;
t758 = (pkin(12) * t838 - t1286) * t1203 + t1115 / 0.2e1 + t436 * t936 + t437 * t935 + t607 * t955 + t605 * t954 + t1328;
t45 = t752 + t758;
t792 = t30 * qJD(1) + t12 * qJD(2) + t45 * qJD(3) + t383 * qJD(4);
t773 = t110 * t1197 + t111 * t1195 - t1283 + t1288;
t23 = -t991 / 0.2e1 + t738 * t947 + (-t982 / 0.4e1 - t1014 / 0.4e1) * Ifges(7,1) + ((-t982 / 0.2e1 + t1014 / 0.2e1) * mrSges(7,3) + t1270) * t726 + t773 + t1330;
t413 = t817 - t769;
t53 = -t819 / 0.2e1 + t808;
t770 = -t1124 * t597 - t1129 * t598 + t1146 * t866 - t824;
t760 = -t1134 * t594 + t726 * t823 + t770;
t775 = (t1293 - t1107 / 0.4e1) * t742 + (t1292 + t1100 / 0.2e1 - t1095 / 0.4e1) * t738;
t767 = t726 * t871 + t775;
t59 = -t797 - (t1194 + t767) * t810 + t760 + t831;
t791 = -t53 * qJD(1) - t23 * qJD(2) + t59 * qJD(3) + t413 * qJD(4);
t25 = t1116 / 0.2e1 + t729 * t1294 + (-t1098 / 0.2e1 + mrSges(7,2) * t1187 - t211 / 0.4e1 + pkin(12) * t1180 - t476 / 0.2e1 + (t1192 / 0.2e1 + t952) * t491) * t742 + (0.3e1 / 0.4e1 * t1072 + mrSges(7,1) * t1187 + t947 + t210 / 0.4e1 + pkin(12) * t1182 + (-t1192 / 0.2e1 + t952) * t490) * t738 + t798 + t825;
t787 = (mrSges(7,1) * t874 + t1195 * t930) * pkin(4);
t327 = t787 + t1339;
t438 = t828 + t769;
t761 = t823 * pkin(12) + pkin(5) * t594 / 0.2e1 + t770;
t768 = pkin(12) * t871 + t775;
t62 = -t797 - (t1194 + t768) * t810 + t761 + t830;
t790 = t25 * qJD(2) - t62 * qJD(3) + t327 * qJD(4) + t438 * qJD(5);
t785 = t863 * t1129 + t1279 * t98 + t1296 * t862 - t1330;
t328 = t787 - t1339;
t61 = -t768 * t810 + t1218 + t761 - t830;
t58 = -t767 * t810 + t1218 + t760 - t831;
t38 = t752 - t758 + t812;
t34 = t771 + t833;
t29 = t757 - t762;
t28 = t772 + t832;
t26 = -t1116 / 0.2e1 - t1283 + t785 + t798 + t1217 * pkin(12);
t24 = t773 + t991 / 0.2e1 + t785 + t1217 * t726;
t21 = t802 * t1202 + t1032 / 0.2e1 + t880 + (t1202 * t826 + t1229) * t276 + t1214 + t1335;
t19 = -t754 + t756 + t1214;
t16 = t751 + t778;
t15 = t983 * t1184 + t803 * t1202 + t774 + (t1202 * t835 + t1222) * t276 + (t974 + t1073) * t1185 - t1302;
t13 = t747 - t755;
t8 = (t591 / 0.2e1 + t1332) * t612 + t759 + t851 + t1268;
t5 = t746 + t749;
t4 = pkin(12) * t1228 + t1145 * t507 + t1167 * t644 + t1183 * t1219 - t1313 + t748 - t781 + t783;
t1 = t745 + t750 + t777;
t31 = [t42 * qJD(2) + t43 * qJD(3) + t1326 * t1338, t5 * qJD(3) + t8 * qJD(4) + t15 * qJD(5) + t28 * qJD(6) + t1024 + (-mrSges(3,2) * t870 - mrSges(3,1) * t894 + t303 * t341 + t302 * t342 + t380 * t307 + t381 * t506 - t380 * t507 + t581 * t641 - t800 * t642 + t660 * t697 + (-mrSges(4,1) * t744 + mrSges(4,2) * t740) * t872 + (t375 + t1234) * t659 + 0.2e1 * (t207 * t380 + t302 * t97 + t303 * t98) * t1202 + 0.2e1 * (-t224 * t380 + t225 * t381 + t582 * t659) * t1204 + 0.2e1 * (-t514 * t800 + t515 * t581 + t662 * t659) * t1205 + m(4) * (-pkin(2) * t872 - t686 * t659 + t688 * t660)) * qJD(2), t1023 + t5 * qJD(2) + t19 * qJD(4) + t21 * qJD(5) + t34 * qJD(6) + ((t296 * t428 + t297 * t429 - t1025) * t1202 + (-t1219 * t399 + t613 * t728 - t1025) * t1204 + (pkin(10) * t612 * t962 - pkin(3) * t613) * t1205) * t1206 + (-t1064 * t399 + t296 * t608 + t297 * t606 + (-mrSges(4,1) + t1233) * t613 + (mrSges(5,3) * t962 + mrSges(4,2)) * t612 - (t1061 + t596) * t398) * qJD(3), t8 * qJD(2) + t19 * qJD(3) + (m(7) * t1301 - t1120 * t765 * t1201 - t776 * mrSges(5,2) - t470 * mrSges(5,1) - (t1119 * t1201 + t1108 + t1109) * t276 + t1311) * qJD(4) + t29 * qJD(5) + t1336, t15 * qJD(2) + t21 * qJD(3) + t29 * qJD(4) + (m(7) * t1304 - mrSges(7,3) * t1281 + t1311) * qJD(5) + t1336, t28 * qJD(2) + t34 * qJD(3) + (-mrSges(7,1) * t187 - mrSges(7,2) * t186) * qJD(6) + t1338 * t1318; -qJD(3) * t6 + qJD(4) * t9 - qJD(5) * t14 + qJD(6) * t27 - t1024, qJD(3) * t10 + qJD(4) * t7 + qJD(5) * t11 + qJD(6) * t22, t1 * qJD(4) + t4 * qJD(5) + t16 * qJD(6) + ((-pkin(3) * t688 + (-t579 * t739 + t580 * t743) * pkin(10)) * t1205 + (t1219 * t305 - t304 * t644 + t647 * t728) * t1204 + (t142 * t428 + t143 * t429 + t282 * t644) * t1202) * t1206 + t850 + ((pkin(10) * t676 + t580 * mrSges(5,3) + t626 / 0.2e1) * t743 + (-pkin(10) * t677 - t579 * mrSges(5,3) + t627 / 0.2e1) * t739 - (-t305 * mrSges(6,3) + t318 / 0.2e1 - t463 / 0.2e1) * t809 - (t320 * t1127 + t319 * t1130 + t304 * mrSges(6,3) - t464 / 0.2e1) * t810 + ((Ifges(5,5) * t1128 + Ifges(6,5) * t1137 + Ifges(5,6) * t1123 + t1258 * t809 - Ifges(4,6)) * t740 + t821 * t744) * t737 + t505 * t1157 + t503 * t1158 + t620 * t1142 + t501 * t1143 + t618 * t1144 + t720 + t728 * t494 - t686 * mrSges(4,2) - pkin(3) * t661 + t647 * t615 + t1219 * t601 + t143 * t606 + t142 * t608 + t282 * t596 + t429 * t442 + t428 * t443 + (t708 - mrSges(4,1)) * t688 + (-t602 + t385) * t644) * qJD(3), t1 * qJD(3) + t13 * qJD(5) + t24 * qJD(6) + t857 + ((t1119 * t232 - t1120 * t231) * t1201 - t780 * t885 + t1262 * t992 - t1263 * t993 - t727 * t1282 + t1038 - t514 * mrSges(5,2) - t515 * mrSges(5,1) - t1083 - t1084 + m(7) * (t231 * t727 + t726 * t846) + t1225 - (-t886 + t1232) * t567 + t846 * mrSges(7,3) + t1325) * qJD(4), t4 * qJD(3) + t13 * qJD(4) + (t1041 + pkin(5) * t1282 - t1263 * t1113 + m(7) * (-pkin(5) * t225 + pkin(12) * t845) + t1262 * t1112 - t1087 - t1086 - t1232 * t567 + t845 * mrSges(7,3) + t967 + t1325) * qJD(5) + t26 * qJD(6) + t848, t16 * qJD(3) + t24 * qJD(4) + t26 * qJD(5) + (-mrSges(7,1) * t98 - mrSges(7,2) * t97 + t308) * qJD(6) + t847; qJD(2) * t6 - qJD(4) * t18 - qJD(5) * t20 + qJD(6) * t33 - t1023, -qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t17 - t850, qJD(4) * t36 + qJD(5) * t41 - qJD(6) * t60 ((-t1120 * t1219 - t1317) * t1201 + t605 * t992 - t607 * t993 + m(7) * (t726 * t838 + t1276) + t990 - t809 * t886 - t810 * t885 + t859 + t708 * pkin(10) + t838 * mrSges(7,3) + t1329) * qJD(4) + t38 * qJD(5) + t58 * qJD(6) + t806, t38 * qJD(4) + (m(7) * (pkin(12) * t837 - t1286) + t605 * t1112 - t607 * t1113 - t1115 + t837 * mrSges(7,3) + t1329) * qJD(5) + t61 * qJD(6) + t805, t58 * qJD(4) + t61 * qJD(5) + (-mrSges(7,1) * t429 - mrSges(7,2) * t428 - t709 * t810) * qJD(6) + t804; -qJD(2) * t9 + qJD(3) * t18 + qJD(5) * t30 - qJD(6) * t53 - t1334, qJD(3) * t2 + qJD(5) * t12 - qJD(6) * t23 - t857, qJD(5) * t45 + qJD(6) * t59 - t806, qJD(5) * t383 + qJD(6) * t413 ((-pkin(5) * t1119 + pkin(12) * t811) * t1200 + t764) * qJD(5) + t328 * qJD(6) + t792, t328 * qJD(5) + (t707 * t726 + t1223) * qJD(6) + t791; qJD(2) * t14 + qJD(3) * t20 - qJD(4) * t30 - t1334, -qJD(3) * t3 - qJD(4) * t12 - qJD(6) * t25 - t848, -qJD(4) * t45 + qJD(6) * t62 - t805, -qJD(6) * t327 - t792, -t438 * qJD(6) (pkin(12) * t707 + t1223) * qJD(6) - t790; -t27 * qJD(2) - t33 * qJD(3) + t53 * qJD(4), -qJD(3) * t17 + qJD(4) * t23 + qJD(5) * t25 - t847, -qJD(4) * t59 - qJD(5) * t62 - t804, qJD(5) * t327 - t791, t790, 0;];
Cq  = t31;
