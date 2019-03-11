% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:45:08
% EndTime: 2019-03-09 14:46:10
% DurationCPUTime: 37.45s
% Computational Cost: add. (56284->1368), mult. (130423->1838), div. (0->0), fcn. (138794->10), ass. (0->669)
t788 = sin(qJ(6));
t792 = cos(qJ(6));
t793 = cos(qJ(5));
t794 = cos(qJ(4));
t972 = t793 * t794;
t789 = sin(qJ(5));
t986 = t789 * t794;
t679 = t788 * t986 - t792 * t972;
t1215 = pkin(3) + pkin(8);
t787 = cos(pkin(6));
t1109 = pkin(1) * t787;
t795 = cos(qJ(2));
t771 = t795 * t1109;
t786 = sin(pkin(6));
t791 = sin(qJ(2));
t995 = t786 * t791;
t622 = -t1215 * t995 + t771;
t796 = -pkin(2) - pkin(9);
t548 = t787 * t796 - t622;
t895 = -qJ(3) * t791 - pkin(1);
t587 = (t795 * t796 + t895) * t786;
t790 = sin(qJ(4));
t356 = t548 * t790 + t587 * t794;
t309 = pkin(10) * t995 + t356;
t994 = t786 * t795;
t702 = pkin(8) * t994 + t791 * t1109;
t768 = pkin(3) * t994;
t623 = t768 + t702;
t781 = t787 * qJ(3);
t586 = t781 + t623;
t695 = t787 * t790 + t794 * t994;
t696 = t787 * t794 - t790 * t994;
t375 = pkin(4) * t695 - pkin(10) * t696 + t586;
t163 = t309 * t793 + t375 * t789;
t980 = t791 * t793;
t570 = -t696 * t789 + t786 * t980;
t117 = pkin(11) * t570 + t163;
t1017 = t117 * t792;
t162 = -t309 * t789 + t793 * t375;
t571 = t696 * t793 + t789 * t995;
t116 = -pkin(11) * t571 + t162;
t113 = pkin(5) * t695 + t116;
t82 = t113 * t788 + t1017;
t1023 = t82 * t679;
t853 = t788 * t793 + t792 * t789;
t681 = t853 * t794;
t1018 = t117 * t788;
t81 = t113 * t792 - t1018;
t1024 = t81 * t681;
t1105 = pkin(5) * t789;
t893 = -t796 + t1105;
t716 = t893 * t794;
t1141 = t716 / 0.2e1;
t1148 = t695 / 0.4e1;
t1102 = pkin(10) * t794;
t1108 = pkin(4) * t790;
t738 = qJ(3) - t1102 + t1108;
t981 = t790 * t796;
t618 = t789 * t738 + t793 * t981;
t561 = -pkin(11) * t986 + t618;
t1009 = t561 * t788;
t714 = t793 * t738;
t894 = -t789 * t796 + pkin(5);
t958 = pkin(11) * t972;
t533 = t790 * t894 + t714 - t958;
t322 = t533 * t792 - t1009;
t1201 = -t322 / 0.2e1;
t355 = t548 * t794 - t790 * t587;
t308 = -pkin(4) * t995 - t355;
t238 = -pkin(5) * t570 + t308;
t1206 = t238 / 0.2e1;
t361 = t570 * t788 + t571 * t792;
t890 = t792 * t570 - t571 * t788;
t1257 = mrSges(7,1) * t361 + mrSges(7,2) * t890;
t1269 = -t361 / 0.2e1;
t1008 = t561 * t792;
t323 = t533 * t788 + t1008;
t499 = -mrSges(7,1) * t679 - mrSges(7,2) * t681;
t502 = -Ifges(7,5) * t681 + Ifges(7,6) * t679;
t1090 = Ifges(7,4) * t679;
t480 = -Ifges(7,2) * t681 + t790 * Ifges(7,6) - t1090;
t1185 = -t480 / 0.4e1;
t504 = -Ifges(7,1) * t681 + t1090;
t905 = t504 / 0.4e1 + t1185;
t651 = Ifges(7,4) * t681;
t482 = -Ifges(7,1) * t679 + t790 * Ifges(7,5) - t651;
t503 = Ifges(7,2) * t679 - t651;
t907 = t482 / 0.4e1 + t503 / 0.4e1;
t345 = Ifges(7,4) * t890;
t143 = Ifges(7,1) * t361 + t695 * Ifges(7,5) + t345;
t1211 = t143 / 0.4e1;
t1258 = -Ifges(7,2) * t361 + t345;
t914 = t1258 / 0.4e1 + t1211;
t1091 = Ifges(7,4) * t361;
t142 = Ifges(7,2) * t890 + t695 * Ifges(7,6) + t1091;
t1213 = t142 / 0.4e1;
t1275 = Ifges(7,1) * t890 - t1091;
t915 = t1213 - t1275 / 0.4e1;
t1293 = t907 * t890 + t1148 * t502 + t1206 * t499 + t1257 * t1141 + (t1269 * t323 + t1201 * t890 + t1024 / 0.2e1 + t1023 / 0.2e1) * mrSges(7,3) + t905 * t361 - t914 * t681 + t915 * t679;
t984 = t790 * t791;
t635 = (-t789 * t984 + t793 * t795) * t786;
t636 = (t789 * t795 + t790 * t980) * t786;
t431 = t635 * t792 - t636 * t788;
t1190 = t431 / 0.2e1;
t432 = t635 * t788 + t636 * t792;
t1189 = t432 / 0.2e1;
t1218 = mrSges(6,2) / 0.2e1;
t1219 = -mrSges(6,1) / 0.2e1;
t1286 = mrSges(7,1) / 0.2e1;
t979 = t791 * t794;
t948 = t786 * t979;
t675 = t787 * t793 + t789 * t948;
t676 = t787 * t789 - t793 * t948;
t476 = t675 * t792 - t676 * t788;
t477 = t675 * t788 + t676 * t792;
t1251 = t476 * t1286 - t477 * mrSges(7,2) / 0.2e1;
t1221 = m(7) * pkin(5);
t960 = -t1221 / 0.2e1;
t1291 = t1251 - t675 * t1219 - t676 * t1218 - (t476 * t792 + t477 * t788) * t960;
t1029 = t793 * mrSges(6,1);
t1032 = t789 * mrSges(6,2);
t1035 = t853 * mrSges(7,2);
t977 = t792 * t793;
t852 = t788 * t789 - t977;
t1037 = t852 * mrSges(7,1);
t1250 = -t1037 / 0.2e1 - t1035 / 0.2e1;
t854 = t788 * t853 - t792 * t852;
t1290 = t1032 / 0.2e1 - t1029 / 0.2e1 + t854 * t960 - t1250;
t767 = pkin(2) * t995;
t621 = t767 + (pkin(9) * t791 - qJ(3) * t795) * t786;
t424 = t794 * t621 + t790 * t623;
t393 = pkin(10) * t994 + t424;
t753 = pkin(4) * t794 + pkin(10) * t790;
t489 = t771 + (-t753 - t1215) * t995;
t227 = -t393 * t789 + t793 * t489;
t166 = -pkin(5) * t948 - pkin(11) * t636 + t227;
t228 = t793 * t393 + t789 * t489;
t189 = pkin(11) * t635 + t228;
t109 = t166 * t792 - t189 * t788;
t110 = t166 * t788 + t189 * t792;
t1262 = Ifges(7,5) * t1189 + Ifges(7,6) * t1190;
t1285 = mrSges(7,2) / 0.2e1;
t1287 = -mrSges(7,1) / 0.2e1;
t1236 = t109 * t1287 + t110 * t1285 - t1262;
t1137 = t853 / 0.4e1;
t1139 = -t852 / 0.4e1;
t711 = Ifges(7,4) * t852;
t566 = -Ifges(7,2) * t853 - t711;
t1178 = t566 / 0.4e1;
t562 = mrSges(7,1) * t853 - mrSges(7,2) * t852;
t1268 = t562 / 0.2e1;
t564 = -Ifges(7,5) * t852 - Ifges(7,6) * t853;
t1241 = t1148 * t564 + t1268 * t238;
t1089 = Ifges(7,4) * t853;
t568 = -Ifges(7,1) * t852 - t1089;
t1034 = t853 * mrSges(7,3);
t933 = -t1034 / 0.2e1;
t1288 = t1275 * t1137 - t853 * t1213 + t82 * t933 + t890 * t1178 + t361 * t568 / 0.4e1 + t1241 + (t1258 + t143) * t1139;
t1284 = -t1257 / 0.2e1;
t90 = t116 * t792 - t1018;
t1283 = -t81 + t90;
t617 = -t789 * t981 + t714;
t560 = t617 - t958;
t352 = t560 * t792 - t1009;
t1280 = -t322 + t352;
t351 = -t560 * t788 - t1008;
t1279 = t323 + t351;
t680 = t853 * t790;
t987 = t789 * t790;
t682 = -t788 * t987 + t790 * t977;
t892 = -t682 * mrSges(7,1) + t680 * mrSges(7,2);
t1278 = qJD(6) * t892;
t1214 = -pkin(11) - pkin(10);
t752 = t1214 * t789;
t754 = t1214 * t793;
t597 = t752 * t788 - t754 * t792;
t889 = t792 * t752 + t754 * t788;
t111 = -t597 * mrSges(7,1) - t889 * mrSges(7,2) + t564;
t1277 = t111 * qJD(6);
t1276 = t238 * t1257;
t1254 = Ifges(7,5) * t890;
t1266 = Ifges(7,6) * t361;
t1274 = t1254 - t1266;
t1071 = Ifges(7,3) * t794;
t723 = t793 * t753;
t983 = t790 * t793;
t545 = pkin(11) * t983 + t794 * t894 + t723;
t970 = t794 * t796;
t634 = t789 * t753 + t793 * t970;
t578 = pkin(11) * t987 + t634;
t333 = t545 * t792 - t578 * t788;
t334 = t545 * t788 + t578 * t792;
t1075 = Ifges(7,6) * t680;
t1083 = Ifges(7,5) * t682;
t849 = -t1083 / 0.2e1 + t1075 / 0.2e1;
t1242 = t334 * t1285 + t1287 * t333 - t849;
t1233 = t1071 / 0.2e1 - t1242;
t912 = t1266 / 0.2e1 - t1254 / 0.2e1;
t1271 = -Ifges(6,5) / 0.2e1;
t526 = pkin(4) * t696 + pkin(10) * t695;
t225 = -t355 * t789 + t793 * t526;
t996 = t695 * t793;
t157 = pkin(5) * t696 + pkin(11) * t996 + t225;
t226 = t793 * t355 + t789 * t526;
t997 = t695 * t789;
t192 = pkin(11) * t997 + t226;
t106 = t157 * t788 + t192 * t792;
t1270 = -t106 / 0.2e1;
t89 = -t116 * t788 - t1017;
t1267 = t82 + t89;
t885 = -t948 / 0.2e1;
t1265 = Ifges(6,3) * t885;
t1103 = pkin(5) * t793;
t775 = -pkin(4) - t1103;
t1127 = t775 / 0.2e1;
t1078 = Ifges(6,6) * t789;
t782 = Ifges(6,5) * t793;
t1245 = t782 - t1078;
t1130 = t1245 / 0.4e1;
t1028 = t793 * mrSges(6,2);
t1033 = t789 * mrSges(6,1);
t740 = t1028 + t1033;
t1132 = t740 / 0.2e1;
t1136 = -t853 / 0.2e1;
t1168 = t597 / 0.2e1;
t1169 = t889 / 0.2e1;
t1170 = -t889 / 0.2e1;
t569 = Ifges(7,1) * t853 - t711;
t1174 = t569 / 0.4e1;
t567 = -Ifges(7,2) * t852 + t1089;
t1176 = t567 / 0.4e1;
t378 = mrSges(6,1) * t571 + mrSges(6,2) * t570;
t1192 = -t378 / 0.2e1;
t282 = mrSges(7,1) * t695 - mrSges(7,3) * t361;
t1204 = -t282 / 0.2e1;
t1223 = m(7) / 0.2e1;
t281 = -mrSges(7,2) * t695 + mrSges(7,3) * t890;
t783 = Ifges(6,4) * t793;
t1244 = -Ifges(6,2) * t789 + t783;
t747 = Ifges(6,1) * t789 + t783;
t896 = -t747 / 0.4e1 - t1244 / 0.4e1;
t1092 = Ifges(6,4) * t789;
t744 = Ifges(6,2) * t793 + t1092;
t748 = Ifges(6,1) * t793 - t1092;
t897 = -t744 / 0.4e1 + t748 / 0.4e1;
t1217 = -t81 / 0.2e1;
t950 = t1217 + t90 / 0.2e1;
t1260 = t1267 * t889 * t1223 + t1257 * t1127 + t308 * t1132 + t695 * t1130 + t281 * t1169 - t361 * t1176 + t890 * t1174 + pkin(4) * t1192 - t896 * t570 + (t1136 * t89 - t1168 * t361 + t1170 * t890 - t852 * t950) * mrSges(7,3) + (t1223 * t1283 + t1204) * t597 + t897 * t571 + t1288;
t876 = t1029 - t1032;
t1021 = -t876 - mrSges(5,1);
t1072 = Ifges(7,3) * t696;
t1074 = Ifges(6,3) * t696;
t474 = t853 * t695;
t1076 = Ifges(7,6) * t474;
t475 = t852 * t695;
t1084 = Ifges(7,5) * t475;
t1249 = -t1245 * t695 + t1072 + t1074 + t1076 + t1084;
t1155 = t680 / 0.2e1;
t1246 = t1155 * t890;
t633 = -t789 * t970 + t723;
t857 = -t633 * t789 + t634 * t793;
t1030 = t790 * Ifges(6,6);
t668 = t1244 * t794 + t1030;
t707 = t747 * t794;
t901 = -t668 / 0.4e1 - t707 / 0.4e1;
t1096 = Ifges(5,4) * t696;
t141 = Ifges(7,5) * t361 + Ifges(7,6) * t890 + t695 * Ifges(7,3);
t300 = Ifges(6,5) * t571 + Ifges(6,6) * t570 + t695 * Ifges(6,3);
t1243 = -Ifges(5,1) * t695 - t1096 + t141 + t300;
t1044 = t681 * mrSges(7,1);
t1046 = t679 * mrSges(7,2);
t965 = -t1044 / 0.2e1 + t1046 / 0.2e1;
t1240 = -Ifges(6,6) * t635 / 0.2e1 + t636 * t1271;
t1043 = t681 * mrSges(7,3);
t602 = -mrSges(7,2) * t790 - t1043;
t1167 = t602 / 0.2e1;
t1239 = t1127 * t499 + t1141 * t562 + t889 * t1167 - t852 * t907 + t853 * t905;
t1124 = t789 / 0.2e1;
t1172 = t571 / 0.2e1;
t182 = -mrSges(7,1) * t890 + mrSges(7,2) * t361;
t563 = t1035 + t1037;
t1235 = (t238 * t789 + t571 * t775) * t1223 + t182 * t1124 + t563 * t1172;
t1234 = Ifges(7,3) * t885 - t1236;
t1039 = t696 * mrSges(5,3);
t1041 = t695 * mrSges(5,3);
t1114 = t794 / 0.2e1;
t1115 = -t794 / 0.2e1;
t1120 = t790 / 0.2e1;
t1122 = -t790 / 0.2e1;
t379 = -mrSges(6,1) * t570 + mrSges(6,2) * t571;
t490 = t740 * t695;
t593 = mrSges(5,1) * t995 - t1039;
t592 = -mrSges(5,2) * t995 - t1041;
t971 = t794 * t592;
t1232 = t1041 * t1114 + t971 / 0.2e1 - t490 * t1115 + t379 * t1120 + (t1039 + t593) * t1122;
t105 = t157 * t792 - t192 * t788;
t1231 = mrSges(7,2) * t1270 + t105 * t1286 + t1076 / 0.2e1 + t1084 / 0.2e1;
t1230 = 0.2e1 * m(7);
t1229 = 2 * qJD(4);
t1228 = m(4) / 0.2e1;
t1227 = m(5) / 0.2e1;
t1226 = -m(6) / 0.2e1;
t1225 = m(6) / 0.2e1;
t1224 = -m(7) / 0.2e1;
t1220 = mrSges(5,1) / 0.2e1;
t1216 = -t82 / 0.2e1;
t1212 = t143 / 0.2e1;
t1210 = -t163 / 0.2e1;
t1209 = -t182 / 0.2e1;
t1208 = -t227 / 0.2e1;
t1207 = t228 / 0.2e1;
t1205 = t281 / 0.2e1;
t1079 = Ifges(6,6) * t695;
t1093 = Ifges(6,4) * t571;
t301 = Ifges(6,2) * t570 + t1079 + t1093;
t1203 = -t301 / 0.4e1;
t1086 = Ifges(6,5) * t695;
t551 = Ifges(6,4) * t570;
t302 = Ifges(6,1) * t571 + t1086 + t551;
t1202 = t302 / 0.2e1;
t1200 = -t323 / 0.2e1;
t1199 = t351 / 0.2e1;
t1198 = t352 / 0.2e1;
t1196 = t890 / 0.2e1;
t1194 = -t361 / 0.4e1;
t1193 = t361 / 0.2e1;
t1191 = -t424 / 0.2e1;
t1054 = t570 * mrSges(6,3);
t439 = -mrSges(6,2) * t695 + t1054;
t1188 = -t439 / 0.2e1;
t1187 = t439 / 0.2e1;
t1053 = t571 * mrSges(6,3);
t440 = mrSges(6,1) * t695 - t1053;
t1186 = -t440 / 0.2e1;
t500 = -mrSges(7,1) * t680 - mrSges(7,2) * t682;
t1183 = -t500 / 0.2e1;
t501 = t1044 - t1046;
t1182 = -t501 / 0.2e1;
t1179 = t563 / 0.2e1;
t1177 = t567 / 0.2e1;
t1175 = t569 / 0.2e1;
t1173 = t570 / 0.2e1;
t1171 = t586 / 0.2e1;
t1045 = t679 * mrSges(7,3);
t604 = mrSges(7,1) * t790 + t1045;
t1166 = -t604 / 0.2e1;
t1165 = t604 / 0.2e1;
t1164 = -t617 / 0.2e1;
t1163 = -t618 / 0.2e1;
t1162 = -t633 / 0.2e1;
t1160 = t668 / 0.2e1;
t1031 = t790 * Ifges(6,5);
t670 = t748 * t794 + t1031;
t1159 = t670 / 0.2e1;
t1158 = -t679 / 0.2e1;
t1157 = t679 / 0.2e1;
t1156 = -t680 / 0.2e1;
t1154 = -t681 / 0.2e1;
t1153 = -t682 / 0.2e1;
t1152 = t682 / 0.2e1;
t1151 = -t695 / 0.2e1;
t1150 = -t695 / 0.4e1;
t1149 = t695 / 0.2e1;
t1147 = t696 / 0.2e1;
t703 = t876 * t794;
t1145 = -t703 / 0.2e1;
t704 = t740 * t790;
t1144 = t704 / 0.2e1;
t705 = t794 * t740;
t1143 = -t705 / 0.2e1;
t1140 = -t852 / 0.2e1;
t1138 = t853 / 0.2e1;
t957 = mrSges(6,3) * t986;
t725 = -mrSges(6,2) * t790 - t957;
t1135 = -t725 / 0.2e1;
t1134 = t725 / 0.2e1;
t727 = mrSges(6,1) * t790 - mrSges(6,3) * t972;
t1133 = -t727 / 0.2e1;
t742 = Ifges(6,5) * t789 + Ifges(6,6) * t793;
t1131 = t742 / 0.2e1;
t1126 = -t788 / 0.2e1;
t1125 = -t789 / 0.2e1;
t1123 = t789 / 0.4e1;
t1119 = t790 / 0.4e1;
t1118 = -t792 / 0.2e1;
t1117 = -t793 / 0.2e1;
t1116 = t793 / 0.2e1;
t1112 = -t796 / 0.2e1;
t1110 = m(7) * t716;
t1107 = pkin(5) * t571;
t1106 = pkin(5) * t788;
t1104 = pkin(5) * t792;
t1101 = t81 * mrSges(7,2);
t1100 = t82 * mrSges(7,1);
t1099 = t89 * mrSges(7,1);
t1098 = t90 * mrSges(7,2);
t1095 = Ifges(5,4) * t790;
t1094 = Ifges(5,4) * t794;
t1088 = Ifges(5,5) * t790;
t1081 = Ifges(5,6) * t794;
t1073 = Ifges(6,3) * t794;
t1070 = pkin(5) * qJD(5);
t1069 = qJ(3) * mrSges(5,1);
t1064 = t162 * mrSges(6,3);
t1063 = t322 * mrSges(7,2);
t1062 = t323 * mrSges(7,1);
t1059 = t351 * mrSges(7,1);
t1058 = t352 * mrSges(7,2);
t701 = pkin(8) * t995 - t771;
t1038 = t701 * mrSges(3,2);
t204 = Ifges(7,4) * t432 + Ifges(7,2) * t431 - Ifges(7,6) * t948;
t205 = Ifges(7,1) * t432 + Ifges(7,4) * t431 - Ifges(7,5) * t948;
t237 = -mrSges(7,1) * t431 + mrSges(7,2) * t432;
t423 = -t790 * t621 + t623 * t794;
t392 = -pkin(4) * t994 - t423;
t293 = -pkin(5) * t635 + t392;
t376 = mrSges(7,2) * t948 + mrSges(7,3) * t431;
t377 = -mrSges(7,1) * t948 - mrSges(7,3) * t432;
t400 = Ifges(6,4) * t636 + Ifges(6,2) * t635 - Ifges(6,6) * t948;
t401 = Ifges(6,1) * t636 + Ifges(6,4) * t635 - Ifges(6,5) * t948;
t451 = -mrSges(6,1) * t635 + mrSges(6,2) * t636;
t487 = -t695 * Ifges(5,2) + Ifges(5,6) * t995 + t1096;
t1040 = t696 * mrSges(5,2);
t1042 = t695 * mrSges(5,1);
t525 = t1040 + t1042;
t534 = mrSges(6,2) * t948 + mrSges(6,3) * t635;
t535 = -mrSges(6,1) * t948 - mrSges(6,3) * t636;
t874 = Ifges(5,1) * t790 + t1094;
t580 = (Ifges(5,5) * t795 + t791 * t874) * t786;
t877 = mrSges(5,1) * t794 - mrSges(5,2) * t790;
t641 = t877 * t995;
t642 = -t781 - t702;
t643 = (-pkin(2) * t795 + t895) * t786;
t647 = -pkin(2) * t787 + t701;
t677 = (mrSges(5,1) * t795 - mrSges(5,3) * t984) * t786;
t678 = (-mrSges(5,2) * t795 + mrSges(5,3) * t979) * t786;
t699 = (mrSges(4,2) * t795 - mrSges(4,3) * t791) * t786;
t700 = -qJ(3) * t994 + t767;
t718 = -mrSges(4,1) * t994 - t787 * mrSges(4,3);
t764 = Ifges(4,5) * t995;
t765 = Ifges(3,5) * t994;
t873 = Ifges(5,2) * t794 + t1095;
t888 = Ifges(7,3) * t948;
t879 = -(Ifges(5,6) * t795 + t791 * t873) * t786 / 0.2e1 + t1265 - t888 / 0.2e1 - t1240 + t1262;
t947 = t1088 / 0.2e1;
t688 = Ifges(5,4) * t695;
t956 = Ifges(5,5) * t995;
t488 = t696 * Ifges(5,1) - t688 + t956;
t985 = t790 * t488;
t5 = t701 * t718 + t700 * t699 + t356 * t678 + t355 * t677 - t586 * t641 + t635 * t301 / 0.2e1 + t622 * t525 + t424 * t592 + t423 * t593 + t163 * t534 + t162 * t535 + t227 * t440 + t308 * t451 + t228 * t439 + t392 * t379 + t82 * t376 + t81 * t377 + t293 * t182 + t205 * t1193 + t204 * t1196 + t636 * t1202 + t401 * t1172 + t400 * t1173 + t143 * t1189 + t142 * t1190 + t110 * t281 + t109 * t282 + t238 * t237 + t580 * t1147 + m(5) * (t355 * t423 + t356 * t424 + t586 * t622) + m(6) * (t162 * t227 + t163 * t228 + t308 * t392) + m(7) * (t109 * t81 + t110 * t82 + t238 * t293) + m(4) * (t642 * t701 + t643 * t700 + t647 * t702) + t879 * t695 + ((-t643 * mrSges(4,3) + t647 * mrSges(4,1) + Ifges(5,5) * t1147 + Ifges(5,6) * t1151 + Ifges(3,4) * t994 + (-Ifges(4,4) + Ifges(3,5) / 0.2e1) * t787 + (-mrSges(3,2) * pkin(1) + Ifges(4,6) * t795) * t786) * t795 + (-t643 * mrSges(4,2) + t985 / 0.2e1 + (-Ifges(3,6) + Ifges(4,5) / 0.2e1) * t787 + (t487 / 0.2e1 - t300 / 0.2e1 - t141 / 0.2e1) * t794 + (-pkin(1) * mrSges(3,1) + (-Ifges(3,4) - Ifges(4,6) + t947 + t1081 / 0.2e1) * t791) * t786 + (t702 + t642) * mrSges(4,1) + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(4,3) + Ifges(5,3)) * t994) * t791) * t786 + (t1038 + t765 / 0.2e1 + t764 / 0.2e1 + (-mrSges(3,1) + mrSges(4,2)) * t702) * t787;
t1055 = t5 * qJD(1);
t1048 = t597 * mrSges(7,3);
t1047 = t618 * mrSges(6,3);
t1036 = t852 * mrSges(7,3);
t1027 = t794 * mrSges(5,2);
t1026 = t796 * mrSges(5,3);
t223 = Ifges(7,4) * t475 + Ifges(7,2) * t474 + Ifges(7,6) * t696;
t224 = Ifges(7,1) * t475 + Ifges(7,4) * t474 + Ifges(7,5) * t696;
t254 = -mrSges(7,1) * t474 + mrSges(7,2) * t475;
t280 = -pkin(5) * t997 + t356;
t364 = -mrSges(7,2) * t696 + mrSges(7,3) * t474;
t365 = mrSges(7,1) * t696 - mrSges(7,3) * t475;
t394 = Ifges(6,6) * t696 - t1244 * t695;
t395 = Ifges(6,5) * t696 - t695 * t748;
t517 = -mrSges(6,2) * t696 + mrSges(6,3) * t997;
t518 = mrSges(6,1) * t696 + mrSges(6,3) * t996;
t878 = mrSges(5,1) * t696 - mrSges(5,2) * t695;
t891 = -Ifges(5,2) * t696 - t688;
t925 = t995 / 0.2e1;
t926 = -t996 / 0.2e1;
t927 = t997 / 0.2e1;
t964 = -Ifges(5,5) * t695 - Ifges(5,6) * t696;
t968 = t593 - t379;
t8 = t302 * t926 + t301 * t927 + (t891 + t488) * t1151 - t696 * t487 / 0.2e1 + t163 * t517 + t162 * t518 - t308 * t490 + t474 * t142 / 0.2e1 + t225 * t440 + t226 * t439 + t82 * t364 + t81 * t365 + t475 * t1212 + t224 * t1193 + t223 * t1196 + t395 * t1172 + t394 * t1173 + t280 * t182 + t106 * t281 + t105 * t282 + t238 * t254 + m(6) * (t162 * t225 + t163 * t226 + t308 * t356) + m(7) * (t105 * t81 + t106 * t82 + t238 * t280) + t586 * t878 + t1243 * t1147 + t1249 * t1149 + t964 * t925 + (-t968 - t1039) * t356 + (t592 + t1041) * t355;
t1025 = t8 * qJD(1);
t380 = Ifges(6,5) * t570 - Ifges(6,6) * t571;
t381 = -Ifges(6,2) * t571 + t551;
t382 = Ifges(6,1) * t570 - t1093;
t9 = t90 * t281 + t89 * t282 + t142 * t1269 + t1258 * t1196 + t890 * t1212 + t1275 * t1193 + m(7) * (t1107 * t238 + t81 * t89 + t82 * t90) + t182 * t1107 + t1276 + t382 * t1172 - t571 * t301 / 0.2e1 + t308 * t378 + t162 * t439 + (t1274 / 0.2e1 + t380 / 0.2e1) * t695 + (-t440 - t1053) * t163 + (-t361 * t82 - t81 * t890) * mrSges(7,3) + (-t1064 + t1202 + t381 / 0.2e1) * t570;
t1022 = t9 * qJD(1);
t1020 = t109 * t792;
t1019 = t110 * t788;
t14 = t81 * t281 - t82 * t282 + t1274 * t1149 + t1276 + (-t82 * mrSges(7,3) - t142 / 0.2e1 + t1275 / 0.2e1) * t361 + (-t81 * mrSges(7,3) + t1258 / 0.2e1 + t1212) * t890;
t1016 = t14 * qJD(1);
t1010 = t356 * t794;
t949 = t786 * t984;
t25 = t477 * t281 + t476 * t282 + t676 * t439 + t675 * t440 + (t525 - t718) * t787 + (-t971 - t699 + (-t182 + t968) * t790) * t995 + m(7) * (-t238 * t949 + t476 * t81 + t477 * t82) + m(6) * (t162 * t675 + t163 * t676 - t308 * t949) + m(5) * (t586 * t787 + (t355 * t790 - t1010) * t995) + m(4) * (-t642 * t787 - t643 * t995);
t1014 = t25 * qJD(1);
t1013 = t323 * t679;
t1012 = t333 * t792;
t1011 = t334 * t788;
t1007 = t571 * t716;
t1006 = t889 * t681;
t1005 = t597 * t679;
t1002 = t675 * t789;
t1001 = t676 * t793;
t1000 = t680 * t681;
t999 = t682 * t679;
t993 = t788 * t604;
t992 = t788 * t679;
t991 = t789 * t668;
t990 = t789 * t725;
t726 = mrSges(6,1) * t794 + mrSges(6,3) * t983;
t989 = t789 * t726;
t988 = t789 * t744;
t982 = t790 * t794;
t978 = t792 * t681;
t976 = t793 * t302;
t975 = t793 * t670;
t724 = -mrSges(6,2) * t794 + mrSges(6,3) * t987;
t974 = t793 * t724;
t973 = t793 * t727;
t784 = t789 ^ 2;
t785 = t793 ^ 2;
t963 = t784 + t785;
t962 = qJD(4) * t790;
t961 = qJD(4) * t794;
t959 = t1221 / 0.2e1;
t953 = t1106 / 0.2e1;
t952 = t1104 / 0.2e1;
t951 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t944 = Ifges(5,3) * t795 / 0.2e1;
t937 = mrSges(7,3) * t1170;
t934 = t1034 / 0.2e1;
t932 = t1030 / 0.2e1;
t931 = t562 * t1115 + (t933 + t934) * t682;
t930 = 0.2e1 * t781 + t702;
t929 = t501 * t1172;
t923 = t789 * t1203;
t922 = -t987 / 0.2e1;
t919 = t983 / 0.2e1;
t917 = t499 * t1115;
t913 = t1203 + t382 / 0.4e1;
t911 = t1199 + t323 / 0.2e1;
t910 = t1198 + t1201;
t909 = -t379 / 0.2e1 + t593 / 0.2e1;
t908 = t381 / 0.4e1 + t302 / 0.4e1;
t906 = -t490 / 0.2e1 - t592 / 0.2e1;
t565 = Ifges(7,5) * t853 - Ifges(7,6) * t852;
t904 = -t565 / 0.4e1 - t742 / 0.4e1;
t903 = t1176 - t568 / 0.4e1;
t902 = t1174 + t1178;
t706 = t794 * t744;
t900 = t670 / 0.4e1 - t706 / 0.4e1;
t899 = t1143 + t1182;
t898 = t1130 + t564 / 0.4e1;
t728 = t790 * t948;
t884 = mrSges(6,3) * (-t784 / 0.2e1 - t785 / 0.2e1);
t883 = -t1053 / 0.2e1 + t1186;
t882 = t1220 + t876 / 0.2e1 - t563 / 0.2e1;
t881 = t791 * t899;
t478 = -Ifges(7,5) * t679 - Ifges(7,6) * t681 + Ifges(7,3) * t790;
t666 = t790 * Ifges(6,3) + t1245 * t794;
t746 = -Ifges(5,2) * t790 + t1094;
t880 = -t478 / 0.2e1 - t666 / 0.2e1 + t746 / 0.2e1;
t741 = t790 * mrSges(5,1) + t1027;
t870 = t1075 - t1083;
t869 = t1072 / 0.2e1 + t1231;
t479 = -Ifges(7,4) * t682 + Ifges(7,2) * t680 + Ifges(7,6) * t794;
t481 = -Ifges(7,1) * t682 + Ifges(7,4) * t680 + Ifges(7,5) * t794;
t601 = -mrSges(7,2) * t794 + mrSges(7,3) * t680;
t603 = mrSges(7,1) * t794 + mrSges(7,3) * t682;
t667 = Ifges(6,6) * t794 - t1244 * t790;
t669 = Ifges(6,5) * t794 - t748 * t790;
t715 = t893 * t790;
t749 = t794 * Ifges(5,1) - t1095;
t797 = -t890 * t479 / 0.4e1 + t602 * t1270 - t162 * t726 / 0.2e1 - t716 * t254 / 0.2e1 + t681 * t223 / 0.4e1 + t679 * t224 / 0.4e1 - t680 * t142 / 0.4e1 - t571 * t669 / 0.4e1 - t570 * t667 / 0.4e1 - t475 * t482 / 0.4e1 - t334 * t281 / 0.2e1 - t715 * t1209 + t724 * t1210 + t682 * t1211 + t601 * t1216 + t603 * t1217 + t481 * t1194 + t364 * t1200 + t365 * t1201 + t333 * t1204 + t280 * t1182 + t238 * t1183 + t474 * t1185 + t634 * t1188 + t440 * t1162 + t517 * t1163 + t518 * t1164 + t105 * t1166 + t356 * t1143 + t308 * t1144 + t749 * t1148 + t225 * t1133 + t226 * t1135 + (t105 * t322 + t106 * t323 - t238 * t715 + t280 * t716 + t333 * t81 + t334 * t82) * t1224;
t860 = -t227 * t789 + t228 * t793;
t802 = (t109 * t1136 + t110 * t1140) * mrSges(7,3) + (-pkin(4) * t392 + pkin(10) * t860) * t1225 + (t109 * t889 + t110 * t597 + t293 * t775) * t1223 - pkin(4) * t451 / 0.2e1 + t293 * t1179 - t392 * t876 / 0.2e1 + t423 * t1220 + t431 * t1176 + t432 * t1174 + t377 * t1169 + t376 * t1168 + t635 * t744 / 0.4e1 + t636 * t747 / 0.4e1 + t204 * t1139 + t205 * t1137 + t237 * t1127;
t824 = t633 * t162 + t634 * t163 + t617 * t225 + t618 * t226;
t833 = t401 / 0.4e1 - pkin(10) * t535 / 0.2e1 + mrSges(6,3) * t1208;
t834 = t400 / 0.4e1 + pkin(10) * t534 / 0.2e1 + mrSges(6,3) * t1207;
t850 = t782 / 0.2e1 - t1078 / 0.2e1;
t851 = Ifges(5,1) / 0.4e1 - Ifges(5,2) / 0.4e1 - Ifges(6,3) / 0.4e1 - Ifges(7,3) / 0.4e1;
t2 = (qJ(3) * t1149 + t1191) * mrSges(5,2) + t802 + t797 + t870 * t1150 + (mrSges(5,2) * t1171 + t488 / 0.4e1 - t1084 / 0.4e1 - t1076 / 0.4e1 - t688 / 0.4e1 + 0.3e1 / 0.4e1 * t956 + t923 + t976 / 0.4e1 + (t1226 * t308 + t909) * t796 + (t1026 / 0.2e1 + t851) * t696 + (-Ifges(5,4) / 0.4e1 + t850) * t695) * t790 + (t1096 / 0.2e1 - t586 * mrSges(5,1) / 0.2e1 + t487 / 0.4e1 - t300 / 0.4e1 - t141 / 0.4e1 - t793 * t395 / 0.4e1 + t394 * t1123 + (t1225 * t356 + t906) * t796 + (0.3e1 / 0.4e1 * Ifges(5,6) + t904) * t995 + (-t1026 / 0.2e1 + t851) * t695) * t794 + (t746 / 0.4e1 - t1069 / 0.2e1 - t666 / 0.4e1 - t478 / 0.4e1) * t696 + t824 * t1226 + t786 * t944 + (t1150 * t668 + t833) * t789 + (t670 * t1148 + t834) * t793;
t32 = t634 * t725 + t617 * t726 + t633 * t727 + t618 * t724 - t715 * t501 + t716 * t500 + t482 * t1153 + t479 * t1154 + t481 * t1158 + t480 * t1155 + t334 * t602 + t322 * t603 + t333 * t604 + t323 * t601 + (t1069 + t667 * t1125 + t669 * t1116 + t796 * t704 - t1094 / 0.2e1 - t880) * t794 + (-qJ(3) * mrSges(5,2) - t749 / 0.2e1 + t991 / 0.2e1 - t975 / 0.2e1 + t796 * t705 + (Ifges(5,4) / 0.2e1 - t850) * t790 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + t951) * t794 + t849) * t790 + m(7) * (t322 * t333 + t323 * t334 - t715 * t716) + m(6) * (-t796 ^ 2 * t982 + t617 * t633 + t618 * t634);
t868 = -t2 * qJD(1) + t32 * qJD(2);
t798 = (t380 / 0.4e1 + t1274 / 0.4e1) * t790 + (mrSges(6,3) * t1164 + t900) * t570 + t162 * t1134 + t163 * t1133 + t308 * t703 / 0.2e1 + t282 * t1199 + t281 * t1198 + t617 * t1187 + t89 * t1165 + t90 * t1167 + t1293;
t811 = t378 * t1112 + (-t1086 / 0.4e1 + t1064 / 0.2e1 - t908) * t789 + (-t1079 / 0.4e1 + mrSges(6,3) * t1210 + (t182 / 0.2e1 + m(7) * t1206) * pkin(5) + t913) * t793;
t823 = mrSges(6,1) * t1208 + mrSges(6,2) * t1207 + t1240;
t829 = t322 * t89 + t323 * t90 + t351 * t81 + t352 * t82;
t3 = t798 + t829 * t1223 + (t929 + t376 * t1126 + t377 * t1118 + (t1007 / 0.4e1 - t1020 / 0.4e1 - t1019 / 0.4e1) * t1230) * pkin(5) + (t951 * t995 + t811) * t794 + t440 * t1163 + (-t1047 / 0.2e1 + t901) * t571 + t823 + t1236;
t866 = mrSges(7,3) * t1013 + t1120 * t502 + t716 * t499;
t39 = -t351 * t604 - t352 * t602 - m(7) * (t322 * t351 + t323 * t352) - pkin(5) * t501 * t972 + t618 * t727 + (t504 / 0.2e1 - t480 / 0.2e1) * t679 - (t322 * mrSges(7,3) - t482 / 0.2e1 - t503 / 0.2e1) * t681 + (t796 * t703 + (t1159 - t706 / 0.2e1 + t1031 / 0.2e1) * t789 + (t707 / 0.2e1 + t1160 - pkin(5) * t1110 + t1047 + t932) * t793) * t794 - t866 + (-t725 - t957) * t617;
t867 = t3 * qJD(1) - t39 * qJD(2);
t803 = t1119 * t1274 + t82 * t1166 + t81 * t1167 + t282 * t1200 + t322 * t1205 + t1293;
t10 = t803 + t888 / 0.2e1 + t1236;
t41 = -t323 * t604 + t480 * t1157 + t504 * t1158 + (t602 + t1043) * t322 + t866 + (t482 + t503) * t1154;
t865 = t10 * qJD(1) + t41 * qJD(2);
t839 = t1117 * t440 + t1125 * t439;
t799 = (mrSges(4,3) + t741 / 0.2e1) * t787 + t930 * t1228 + (t768 + t930) * t1227 + (t793 * t162 + t789 * t163 + t617 * t675 + t618 * t676 + t728 * t796) * t1225 + (t322 * t476 + t323 * t477 - t716 * t949 - t81 * t852 + t82 * t853) * t1223 + t476 * t1165 + t477 * t1167 + t675 * t727 / 0.2e1 + t676 * t1134 + t1042 / 0.2e1 + t1040 / 0.2e1 + t282 * t1140 + t281 * t1138 - t839;
t858 = t423 * t794 + t424 * t790;
t810 = -m(4) * t702 / 0.2e1 - m(5) * t858 / 0.2e1 + (-t392 * t794 + t790 * t860) * t1226 + (-t109 * t680 + t110 * t682 - t293 * t794) * t1224 + t377 * t1155 + t376 * t1153;
t18 = (-t678 / 0.2e1 + t535 * t1124 + t534 * t1117 + t786 * t881) * t790 + (-t677 / 0.2e1 + t451 / 0.2e1 + t237 / 0.2e1) * t794 + t799 + t810;
t91 = t853 * t602 - t852 * t604 + t990 + t973 + mrSges(4,3) + (m(5) + m(4)) * qJ(3) + m(7) * (-t322 * t852 + t323 * t853) + m(6) * (t617 * t793 + t618 * t789) + t741;
t864 = t18 * qJD(1) + t91 * qJD(2);
t819 = (t1117 * t571 + t1124 * t570) * mrSges(6,3) + t839;
t840 = t1153 * t282 + t1156 * t281;
t806 = (-t1152 * t361 + t1246) * mrSges(7,3) + t819 * t790 + (-t1107 * t794 - t1267 * t680 + t1283 * t682) * t1223 + t840;
t19 = (t1284 + t1192) * t794 + t806 - t1291;
t818 = (t999 / 0.2e1 - t1000 / 0.2e1) * mrSges(7,3) + t602 * t1156 + t604 * t1153;
t809 = (-t973 / 0.2e1 - t990 / 0.2e1 + t794 * t884) * t790 + (-t1103 * t794 ^ 2 - t1279 * t680 + t1280 * t682) * t1223 + t818;
t45 = (-t499 / 0.2e1 + t1145) * t794 + t809 + t1290;
t863 = -t19 * qJD(1) - t45 * qJD(2);
t813 = (t1153 * t361 + t1246) * mrSges(7,3) + t1257 * t1115 + t840;
t26 = t813 - t1251;
t812 = t917 + t818;
t53 = t812 - t1250;
t862 = -t26 * qJD(1) - t53 * qJD(2);
t861 = -t225 * t789 + t226 * t793;
t856 = t1001 - t1002;
t855 = t978 + t992;
t848 = t1218 * t226 + t1219 * t225;
t847 = mrSges(6,1) * t1162 + t1218 * t634;
t846 = -t1028 / 0.2e1 - t1033 / 0.2e1;
t844 = pkin(10) * t1135 + t901;
t843 = t1048 / 0.2e1 + t903;
t842 = t937 + t902;
t838 = m(7) * (t105 * t792 + t106 * t788);
t837 = m(7) * (t597 * t853 - t852 * t889);
t805 = ((-t162 * t789 + t163 * t793 - t356) * t794 + (t308 + t861) * t790) * t1226 + (-t105 * t680 + t106 * t682 + t238 * t790 - t280 * t794 - t1023 - t1024) * t1224 + t281 * t1157 + t365 * t1155 - t681 * t1204 + t364 * t1153;
t808 = (t1136 * t476 + t1140 * t477) * mrSges(7,3) + (-t1002 / 0.2e1 + t1001 / 0.2e1) * mrSges(6,3) + (pkin(4) * t949 + pkin(10) * t856) * t1225 + (t476 * t889 + t477 * t597 - t775 * t949) * t1223;
t15 = (t254 / 0.2e1 + mrSges(5,2) * t925 - t1041 / 0.2e1 + t439 * t1117 + t440 * t1124 + t906) * t794 + (t1209 + t1039 / 0.2e1 + t518 * t1124 + t517 * t1117 + t882 * t995 + t909) * t790 + t805 + t808;
t211 = m(6) * (-0.1e1 + t963) * t982 + (-t999 + t1000 - t982) * m(7);
t801 = (t725 * t1116 + t727 * t1125 + t1144 + t1183) * t794 + (-t989 / 0.2e1 + t974 / 0.2e1 - t899) * t790 + ((-t617 * t789 + t618 * t793) * t794 + (t857 - 0.2e1 * t970) * t790) * t1225 + (-t322 * t681 - t333 * t680 + t334 * t682 + t715 * t794 + t716 * t790 - t1013) * t1223 + t602 * t1158 + t603 * t1156 + t604 * t1154 + t601 * t1152;
t43 = -t837 / 0.2e1 + t801;
t832 = -t15 * qJD(1) + t43 * qJD(2) + t211 * qJD(3);
t830 = t282 * t1168 + t281 * t1170 + t1284 * t775;
t800 = t903 * t679 - t902 * t681 + (pkin(10) * t1133 + t900) * t793 + (t1006 / 0.2e1 + t1005 / 0.2e1 - t911 * t853 - t910 * t852) * mrSges(7,3) + (t1279 * t889 + t1280 * t597) * t1223 + pkin(4) * t1145 - t597 * t1165 + t1239;
t814 = t740 * t1112 + t896 * t789 + pkin(10) * t884 + ((m(7) * t1127 + t1179) * pkin(5) + t897) * t793;
t24 = t800 + t844 * t789 + (t850 + t898) * t790 + (t601 * t1126 + t501 * t1124 + t603 * t1118 + (t716 * t1123 - t1012 / 0.4e1 - t1011 / 0.4e1) * t1230) * pkin(5) + (t814 - t951) * t794 + t847 + t1242;
t62 = t775 * t562 - (t1177 - t568 / 0.2e1) * t853 - (t566 / 0.2e1 + t1175) * t852;
t54 = -pkin(4) * t740 + t748 * t1124 - t988 / 0.2e1 + (t747 / 0.2e1 + t1244 / 0.2e1) * t793 + t62 + (m(7) * t775 + t563) * t1105;
t6 = (t364 * t1126 + t365 * t1118 - t838 / 0.2e1 + t1235) * pkin(5) + (t1086 / 0.2e1 + t883 * pkin(10) + t908) * t793 + (-t1079 / 0.2e1 + (t1188 + t1054 / 0.2e1) * pkin(10) + t913) * t789 - t951 * t696 + t848 - t1231 + t1260;
t822 = t986 * t1221;
t826 = t855 * t960 + t965;
t60 = (t1132 + t1268 + t846) * t794 + t822 / 0.2e1 + t826;
t827 = t6 * qJD(1) + t24 * qJD(2) - t60 * qJD(3) + t54 * qJD(4);
t12 = t843 * t361 - t842 * t890 + t852 * t914 + t853 * t915 - t1241 + t830 + t869;
t807 = t1119 * t564 + t1166 * t597 + t679 * t843 - t681 * t842 + t1239;
t38 = t807 - t1233;
t92 = t931 - t965;
t825 = -t12 * qJD(1) + t38 * qJD(2) + t92 * qJD(3) + t62 * qJD(4);
t112 = (t1170 + t1169) * mrSges(7,2);
t815 = (t792 * t1205 + t282 * t1126 + (t1118 * t890 + t1126 * t361) * mrSges(7,3)) * pkin(5) - t912;
t21 = t950 * mrSges(7,2) + (t1216 - t89 / 0.2e1) * mrSges(7,1) + t815 + t912;
t58 = -t910 * mrSges(7,2) + t911 * mrSges(7,1) + (t602 * t1118 + t993 / 0.2e1 + (-t992 / 0.2e1 - t978 / 0.2e1) * mrSges(7,3)) * pkin(5);
t735 = (mrSges(7,1) * t788 + mrSges(7,2) * t792) * pkin(5);
t817 = -t21 * qJD(1) + t58 * qJD(2) - t112 * qJD(4) + t735 * qJD(5);
t717 = t735 * qJD(6);
t93 = t931 + t965;
t61 = t740 * t1115 - t822 / 0.2e1 + t846 * t794 + t826 + t931;
t55 = -t1062 / 0.2e1 - t1063 / 0.2e1 + t602 * t952 + t953 * t1045 + t1059 / 0.2e1 - t1058 / 0.2e1 + t502 - (-mrSges(7,3) * t978 + t993) * pkin(5) / 0.2e1;
t52 = t812 + t1250;
t44 = t1115 * t703 - t1290 + t809 + t917;
t42 = t837 / 0.2e1 + t801;
t37 = t807 + t1233;
t27 = t813 + t1251;
t23 = (t1011 + t1012) * t959 + t800 + t1073 / 0.2e1 + t814 * t794 + t898 * t790 + t601 * t953 + t603 * t952 + t983 * t1271 - t847 + ((t1110 / 0.2e1 + t501 / 0.2e1) * pkin(5) + t844 + t932) * t789 + t1233;
t22 = -t1101 / 0.2e1 - t1100 / 0.2e1 - t1098 / 0.2e1 + t1099 / 0.2e1 + t815 - t912;
t20 = t806 + (t1257 + t378) * t1115 + t1291;
t17 = t799 + t677 * t1114 + t678 * t1120 + t535 * t922 + t534 * t919 + (t795 * mrSges(4,1) + t790 * t881) * t786 - t810 + (t451 + t237) * t1115;
t16 = -t805 + (t1027 / 0.2e1 + t882 * t790) * t995 + t254 * t1115 + t182 * t1120 + t518 * t922 + t517 * t919 + t972 * t1187 + t986 * t1186 + t808 + t1232;
t13 = t1048 * t1269 + t567 * t1194 + t82 * t934 - t830 + t869 + (t1174 + t937) * t890 + t1288;
t11 = t803 + t1234;
t7 = t382 * t1123 + t923 + t1074 / 0.2e1 - t848 + t869 + t364 * t953 + Ifges(6,6) * t927 + t365 * t952 + Ifges(6,5) * t926 + t908 * t793 + t819 * pkin(10) + (t838 / 0.2e1 + t1235) * pkin(5) + t1260;
t4 = (t1019 + t1020) * t959 + t883 * t618 + t798 + (pkin(5) * t1007 + t829) * t1223 + pkin(5) * t929 - t823 + t1265 + t376 * t953 + t377 * t952 + t811 * t794 + t901 * t571 + t1234;
t1 = t802 - t797 - t985 / 0.4e1 - t794 * t487 / 0.4e1 + t877 * t1171 + mrSges(5,2) * t1191 + t975 * t1150 + qJ(3) * t878 / 0.2e1 + t833 * t789 + t834 * t793 + (t944 + (t947 + (Ifges(5,6) / 0.2e1 + t904) * t794) * t791) * t786 + t824 * t1225 + t395 * t972 / 0.4e1 - t394 * t986 / 0.4e1 + t301 * t987 / 0.4e1 + (-t1081 - t1088) * t995 / 0.4e1 + (-t746 / 0.4e1 - t874 / 0.4e1) * t696 + (t666 + t478) * t696 / 0.4e1 - (t976 + t891) * t790 / 0.4e1 + t1249 * t1119 + t1243 * t794 / 0.4e1 + (-t1245 * t790 + t1071 + t1073 + t870 + t873 + t991) * t1148 + ((t308 * t790 - t1010) * t1225 + t1232) * t796;
t28 = [qJD(2) * t5 + qJD(3) * t25 + qJD(4) * t8 + qJD(5) * t9 + qJD(6) * t14, t17 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + t11 * qJD(6) + t1055 + (t622 * t741 + t228 * t725 + t227 * t727 + t716 * t237 + t702 * mrSges(4,2) - t702 * mrSges(3,1) + t392 * t705 - t701 * mrSges(4,3) - qJ(3) * t641 + t618 * t534 + t617 * t535 + t110 * t602 + t109 * t604 + t293 * t501 + 0.2e1 * (qJ(3) * t622 + t796 * t858) * t1227 + 0.2e1 * (-pkin(2) * t702 - qJ(3) * t701) * t1228 + 0.2e1 * (t617 * t227 + t618 * t228 - t392 * t970) * t1225 + t323 * t376 + t322 * t377 + 0.2e1 * (t109 * t322 + t110 * t323 + t293 * t716) * t1223 + t482 * t1189 + t480 * t1190 + t205 * t1158 + t636 * t1159 + t635 * t1160 + t204 * t1154 + (-t424 * mrSges(5,3) + t796 * t678 + t879) * t790 + (t580 / 0.2e1 - t423 * mrSges(5,3) + t400 * t1125 + t401 * t1116 + (t677 - t451) * t796) * t794 + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t1114 + Ifges(5,6) * t1122 - Ifges(4,4)) * t795 + (-qJ(3) * mrSges(4,1) + t1120 * t749 + t794 * t880 - Ifges(3,6)) * t791) * t786 + t764 + t765 + t1038) * qJD(2), t1014 + t17 * qJD(2) + 0.2e1 * ((-t476 * t680 + t477 * t682 + t728) * t1223 + (t790 * t856 + t728) * t1225) * qJD(3) + t16 * qJD(4) + t20 * qJD(5) + t27 * qJD(6), t1025 + t1 * qJD(2) + t16 * qJD(3) + t7 * qJD(5) + t13 * qJD(6) + ((t105 * t889 + t106 * t597 + t280 * t775) * t1223 + (-pkin(4) * t356 + pkin(10) * t861) * t1225) * t1229 + (t775 * t254 + t696 * t1131 + t224 * t1138 + t223 * t1140 + t565 * t1147 + t597 * t364 + t889 * t365 + t280 * t563 + t474 * t1177 + t475 * t1175 + pkin(4) * t490 - t355 * mrSges(5,2) - t106 * t1036 - t105 * t1034 + t964 + (t394 / 0.2e1 + pkin(10) * t517 + t226 * mrSges(6,3) + t747 * t1151) * t793 + (t395 / 0.2e1 - pkin(10) * t518 - t225 * mrSges(6,3) + t744 * t1149) * t789 + t1021 * t356) * qJD(4), t1022 + t4 * qJD(2) + t20 * qJD(3) + t7 * qJD(4) + (-t163 * mrSges(6,1) - t162 * mrSges(6,2) - t1098 + t1099 + t1274 + t380) * qJD(5) + t22 * qJD(6) + (m(7) * (t788 * t90 + t792 * t89) + (-t361 * t788 - t792 * t890) * mrSges(7,3)) * t1070, t1016 + t11 * qJD(2) + t27 * qJD(3) + t13 * qJD(4) + t22 * qJD(5) + (-t1100 + t1274 - t1101) * qJD(6); qJD(3) * t18 - qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t10 - t1055, qJD(3) * t91 + qJD(4) * t32 - qJD(5) * t39 + qJD(6) * t41, m(7) * (t680 * t852 + t682 * t853) * qJD(3) + t42 * qJD(4) + t44 * qJD(5) + t52 * qJD(6) + t864, t42 * qJD(3) + t23 * qJD(5) + t37 * qJD(6) + (-Ifges(5,5) + t988 / 0.2e1 + t747 * t1117 + (-m(6) * pkin(4) + t1021) * t796) * t962 + (t1131 + t565 / 0.2e1 - Ifges(5,6) - t796 * mrSges(5,2)) * t961 + t868 + (t667 * t1116 + t669 * t1124 + t775 * t500 + t481 * t1138 + t479 * t1140 - t715 * t563 + pkin(4) * t704 + t569 * t1153 + t567 * t1155 + t889 * t603 + t597 * t601 - t334 * t1036 - t333 * t1034 + m(7) * (t333 * t889 + t334 * t597 - t715 * t775) + (m(6) * t857 + t974 - t989) * pkin(10) + t857 * mrSges(6,3)) * qJD(4), t44 * qJD(3) + t23 * qJD(4) + (-t618 * mrSges(6,1) - t617 * mrSges(6,2) - Ifges(6,5) * t986 - Ifges(6,6) * t972 - t1058 + t1059 + t502) * qJD(5) + t55 * qJD(6) + (m(7) * (t351 * t792 + t352 * t788) + t855 * mrSges(7,3)) * t1070 + t867, t52 * qJD(3) + t37 * qJD(4) + t55 * qJD(5) + (-t1062 + t502 - t1063) * qJD(6) + t865; -qJD(2) * t18 - qJD(4) * t15 + qJD(5) * t19 + qJD(6) * t26 - t1014, qJD(4) * t43 + qJD(5) * t45 + qJD(6) * t53 - t864, t211 * qJD(4), t61 * qJD(5) + t93 * qJD(6) + (t679 * t852 + t681 * t853) * qJD(4) * mrSges(7,3) + (mrSges(6,3) * t963 - mrSges(5,2)) * t961 + (t563 + t1021) * t962 + ((t1102 * t963 - t1108) * t1225 + (t775 * t790 - t1005 - t1006) * t1223) * t1229 + t832, t61 * qJD(4) + (m(7) * (-t682 * t1104 - t1106 * t680) - t876 * t790 + t892) * qJD(5) + t1278 - t863, t93 * qJD(4) + qJD(5) * t892 + t1278 - t862; qJD(2) * t2 + qJD(3) * t15 + qJD(5) * t6 - qJD(6) * t12 - t1025, -qJD(3) * t43 + qJD(5) * t24 + qJD(6) * t38 - t868, -qJD(5) * t60 + qJD(6) * t92 - t832, qJD(5) * t54 + qJD(6) * t62 (-pkin(10) * t876 + t111 + t1245) * qJD(5) + t1277 + (m(7) * (-t597 * t792 + t788 * t889) - t854 * mrSges(7,3)) * t1070 + t827, t111 * qJD(5) + t1277 + t825; -qJD(2) * t3 - qJD(3) * t19 - qJD(4) * t6 + qJD(6) * t21 - t1022, -qJD(3) * t45 - qJD(4) * t24 - qJD(6) * t58 - t867, t60 * qJD(4) + t863, qJD(6) * t112 - t827, -t717, -t717 - t817; -qJD(2) * t10 - qJD(3) * t26 + qJD(4) * t12 - qJD(5) * t21 - t1016, -qJD(3) * t53 - qJD(4) * t38 + qJD(5) * t58 - t865, -t92 * qJD(4) + t862, -qJD(5) * t112 - t825, t817, 0;];
Cq  = t28;
