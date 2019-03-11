% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:44
% EndTime: 2019-03-10 01:27:46
% DurationCPUTime: 34.91s
% Computational Cost: add. (58400->1093), mult. (122249->1406), div. (0->0), fcn. (130182->8), ass. (0->567)
t1037 = Ifges(6,1) + Ifges(7,1);
t1173 = Ifges(7,4) + Ifges(6,5);
t1051 = cos(qJ(5));
t649 = sin(qJ(4));
t651 = sin(qJ(2));
t1052 = cos(qJ(4));
t652 = cos(qJ(3));
t871 = t1052 * t652;
t650 = sin(qJ(3));
t924 = t650 * t651;
t569 = t649 * t924 - t651 * t871;
t648 = sin(qJ(5));
t925 = t649 * t652;
t600 = -t1052 * t650 - t925;
t718 = t651 * t600;
t451 = t1051 * t718 + t569 * t648;
t1221 = t451 / 0.2e1;
t1108 = -pkin(9) - pkin(8);
t805 = t1108 * t1052;
t767 = t650 * t805;
t524 = t1108 * t925 + t767;
t598 = t600 * pkin(10);
t1150 = t598 + t524;
t890 = t649 * t1108;
t816 = t650 * t890;
t526 = -t1108 * t871 + t816;
t753 = t649 * t650 - t871;
t597 = t753 * pkin(10);
t702 = -t597 + t526;
t290 = t1051 * t1150 - t648 * t702;
t1262 = t1221 * t290;
t694 = -t1051 * t600 - t648 * t753;
t1156 = qJ(6) * t694;
t497 = -t1051 * t753 + t600 * t648;
t634 = -t652 * pkin(3) - pkin(2);
t547 = pkin(4) * t753 + t634;
t284 = -pkin(5) * t497 - t1156 + t547;
t1048 = m(7) * t284;
t314 = -mrSges(7,1) * t497 - mrSges(7,3) * t694;
t1261 = -t314 - t1048;
t1017 = Ifges(7,5) * t694;
t1025 = Ifges(6,4) * t694;
t776 = -Ifges(7,3) * t497 + t1017;
t835 = t1037 * t497;
t1260 = t1017 - t1025 + t776 + t835;
t1181 = t1051 * t702 + t1150 * t648;
t1238 = t1181 * mrSges(7,1);
t1239 = t1181 * mrSges(6,1);
t1255 = t290 * mrSges(7,3);
t1256 = t290 * mrSges(6,2);
t1259 = -t1238 - t1239 + t1255 - t1256;
t688 = -t1051 * t569 + t648 * t718;
t1081 = -t688 / 0.2e1;
t1085 = -t451 / 0.2e1;
t1027 = Ifges(6,4) * t688;
t432 = Ifges(7,5) * t688;
t653 = cos(qJ(2));
t250 = -Ifges(7,6) * t653 - Ifges(7,3) * t451 + t432;
t1133 = t1037 * t451 - t1027 + t250 + t432;
t1019 = Ifges(7,5) * t451;
t437 = Ifges(6,4) * t451;
t1153 = t1037 * t688 - t1173 * t653 - t1019 + t437;
t1135 = -Ifges(6,2) * t688 + t1153 + t437;
t1179 = t688 / 0.2e1;
t1184 = mrSges(6,1) * t688 + mrSges(6,2) * t451;
t1230 = Ifges(7,3) * t688 + t1019;
t252 = Ifges(6,2) * t451 - Ifges(6,6) * t653 + t1027;
t646 = t651 * pkin(7);
t609 = pkin(3) * t924 + t646;
t500 = -pkin(4) * t718 + t609;
t1258 = t252 * t1081 + t1230 * t1085 + t1133 * t1179 + t1135 * t1221 + t500 * t1184;
t695 = -t1256 / 0.2e1 + t1255 / 0.2e1 - t1238 / 0.2e1 - t1239 / 0.2e1;
t1202 = t497 * mrSges(7,2);
t1242 = t1202 / 0.2e1;
t1257 = 0.2e1 * t1242;
t1104 = t1181 / 0.2e1;
t1162 = Ifges(7,6) * t694;
t1164 = Ifges(6,6) * t694;
t1205 = Ifges(6,5) * t497;
t1207 = Ifges(7,4) * t497;
t1254 = t1207 + t1162 + t1205 - t1164;
t1026 = Ifges(6,4) * t497;
t1018 = Ifges(7,5) * t497;
t1149 = t1037 * t694 - t1018 + t1026;
t1253 = -Ifges(6,2) * t694 + t1026 + t1149;
t1250 = -pkin(5) * t1181 + qJ(6) * t290;
t1249 = -t1051 * t1181 + t290 * t648;
t1041 = t648 * pkin(4);
t631 = qJ(6) + t1041;
t905 = t1051 * pkin(4);
t632 = -t905 - pkin(5);
t1248 = t632 * t1181 + t290 * t631;
t1157 = t694 * mrSges(7,2);
t1210 = -t1157 / 0.2e1;
t1247 = -t631 * t1210 - t632 * t1242;
t1028 = Ifges(5,4) * t600;
t1038 = mrSges(7,2) + mrSges(6,3);
t1039 = t651 * pkin(8);
t611 = -pkin(2) * t653 - pkin(1) - t1039;
t599 = t652 * t611;
t922 = t651 * t652;
t804 = -pkin(9) * t922 + t599;
t499 = (-pkin(7) * t650 - pkin(3)) * t653 + t804;
t920 = t652 * t653;
t546 = pkin(7) * t920 + t650 * t611;
t516 = -pkin(9) * t924 + t546;
t503 = t649 * t516;
t324 = t1052 * t499 - t503;
t558 = t569 * pkin(10);
t294 = t324 + t558;
t262 = -pkin(4) * t653 + t294;
t872 = t1052 * t516;
t325 = t649 * t499 + t872;
t559 = pkin(10) * t718;
t295 = t559 + t325;
t928 = t648 * t295;
t117 = t1051 * t262 - t928;
t870 = t1051 * t295;
t929 = t648 * t262;
t118 = t870 + t929;
t1209 = mrSges(6,3) * t497;
t1217 = t497 / 0.4e1;
t1219 = -t497 / 0.4e1;
t1220 = t451 / 0.4e1;
t1222 = -t451 / 0.4e1;
t1182 = mrSges(7,1) * t688 - mrSges(7,3) * t451;
t1183 = mrSges(7,1) * t694 - mrSges(7,3) * t497;
t1185 = mrSges(6,1) * t694 + mrSges(6,2) * t497;
t1243 = -t694 / 0.4e1;
t1244 = -t688 / 0.4e1;
t775 = -pkin(5) * t451 - qJ(6) * t688;
t201 = t500 + t775;
t781 = Ifges(6,2) * t497 + t1025;
t1223 = t1243 * t252 + t1244 * t781 + t201 * t1183 / 0.2e1 + t284 * t1182 / 0.2e1 + t500 * t1185 / 0.2e1 + t547 * t1184 / 0.2e1;
t1231 = Ifges(7,3) * t694 + t1018;
t715 = Ifges(5,4) * t718;
t684 = -Ifges(5,1) * t569 - Ifges(5,5) * t653 + t715;
t1237 = Ifges(5,2) * t569 + t684 + t715;
t1029 = Ifges(5,4) * t569;
t683 = Ifges(5,2) * t718 - Ifges(5,6) * t653 - t1029;
t697 = -t569 * mrSges(5,1) + mrSges(5,2) * t718;
t699 = Ifges(5,1) * t718 + t1029;
t722 = Ifges(5,2) * t753;
t704 = -t722 - t1028;
t723 = Ifges(5,4) * t753;
t705 = -Ifges(5,1) * t600 - t723;
t714 = mrSges(5,3) * t718;
t721 = t753 * mrSges(5,2);
t724 = Ifges(5,1) * t753;
t759 = -Ifges(5,5) * t753 + Ifges(5,6) * t600 + t1254;
t960 = t569 * mrSges(5,3);
t965 = t694 * mrSges(6,3);
t1246 = t1038 * (t1104 * t688 + t1262) + t118 * t965 / 0.2e1 - t1223 + t524 * t714 / 0.2e1 + t1237 * t753 / 0.4e1 - (Ifges(5,2) * t600 + t705 - t723) * t718 / 0.4e1 + t1133 * t1243 + t1260 * t1244 + t759 * t653 / 0.4e1 + t1230 * t1217 - t526 * t960 / 0.2e1 + t1231 * t1220 + t1253 * t1222 + t1135 * t1219 - t609 * (-t600 * mrSges(5,1) - t721) / 0.2e1 + t569 * (-t724 + t1028) / 0.4e1 - t569 * t704 / 0.4e1 - t634 * t697 / 0.2e1 + t600 * t699 / 0.4e1 + t117 * t1209 / 0.2e1 - t600 * t683 / 0.4e1;
t1112 = m(7) / 0.2e1;
t879 = -t1202 / 0.2e1;
t1204 = t451 * mrSges(7,2);
t1241 = t1204 / 0.2e1;
t1036 = Ifges(7,5) - Ifges(6,4);
t1203 = t451 * mrSges(6,3);
t919 = t653 * qJ(6);
t103 = t118 - t919;
t1236 = t103 * t1210;
t1057 = t652 / 0.2e1;
t643 = Ifges(4,4) * t652;
t613 = Ifges(4,1) * t650 + t643;
t1235 = t613 * t1057;
t813 = t905 / 0.2e1;
t1159 = t688 * mrSges(7,2);
t899 = t631 * t1159;
t1228 = t813 * t1203 + t899 / 0.2e1;
t1163 = Ifges(7,6) * t688;
t1165 = Ifges(6,6) * t688;
t1206 = Ifges(6,5) * t451;
t1208 = Ifges(7,4) * t451;
t1142 = t1206 + t1208 + t1163 - t1165;
t1123 = t1205 / 0.2e1 + t1207 / 0.2e1 + t1162 / 0.2e1 - t1164 / 0.2e1;
t1126 = -t1208 / 0.2e1 - t1206 / 0.2e1 - t1163 / 0.2e1 + t1165 / 0.2e1;
t1224 = t284 * t1183 + t547 * t1185;
t1218 = t497 / 0.2e1;
t876 = t1157 / 0.2e1;
t525 = t652 * t890 + t767;
t464 = t598 + t525;
t523 = t652 * t805 - t816;
t703 = t597 + t523;
t288 = -t1051 * t703 + t464 * t648;
t291 = t1051 * t464 + t648 * t703;
t822 = mrSges(6,3) * t905;
t802 = t497 * t822;
t985 = t291 * mrSges(7,3);
t986 = t291 * mrSges(6,2);
t991 = t288 * mrSges(7,1);
t992 = t288 * mrSges(6,1);
t806 = -t991 / 0.2e1 - t992 / 0.2e1 + t985 / 0.2e1 - t986 / 0.2e1;
t894 = -t1041 / 0.2e1;
t820 = mrSges(6,3) * t894;
t1110 = m(6) * pkin(4);
t908 = t1110 / 0.2e1;
t962 = t525 * mrSges(5,2);
t964 = t523 * mrSges(5,1);
t1200 = -(t288 * t632 + t291 * t631) * t1112 - t964 / 0.2e1 + t962 / 0.2e1 - (-t1051 * t288 + t291 * t648) * t908 - t694 * t820 + t802 / 0.2e1 - t806 + t1247;
t1113 = -m(7) / 0.2e1;
t923 = t650 * t653;
t907 = pkin(7) * t923;
t515 = t804 - t907;
t346 = t1052 * t515 - t503;
t302 = t558 + t346;
t345 = -t649 * t515 - t872;
t729 = -t559 + t345;
t143 = -t1051 * t729 + t302 * t648;
t144 = t1051 * t302 + t648 * t729;
t995 = t144 * mrSges(7,3);
t996 = t144 * mrSges(6,2);
t997 = t143 * mrSges(7,1);
t998 = t143 * mrSges(6,1);
t807 = -t997 / 0.2e1 - t998 / 0.2e1 + t995 / 0.2e1 - t996 / 0.2e1;
t909 = mrSges(6,3) * t1041;
t827 = t688 * t909;
t881 = -t1204 / 0.2e1;
t979 = t346 * mrSges(5,2);
t980 = t345 * mrSges(5,1);
t1199 = (t143 * t632 + t144 * t631) * t1113 - t980 / 0.2e1 + t979 / 0.2e1 - (-t1051 * t143 + t144 * t648) * t1110 / 0.2e1 + t632 * t881 + t827 / 0.2e1 - t807 + t1228;
t716 = t1123 + t695;
t814 = -t905 / 0.2e1;
t1198 = -(-pkin(4) * t1249 + t1248) * t1112 - t876 * t1041 + t814 * t1202 - t716 + t1247;
t641 = m(7) * qJ(6) + mrSges(7,3);
t1194 = qJD(5) * t641;
t1193 = t103 * t1112;
t1191 = t641 * qJD(6);
t133 = t294 * t648 + t870;
t134 = t1051 * t294 - t928;
t1000 = t134 * mrSges(6,2);
t1001 = t133 * mrSges(7,1);
t1002 = t133 * mrSges(6,1);
t999 = t134 * mrSges(7,3);
t808 = -t1001 / 0.2e1 - t1002 / 0.2e1 + t999 / 0.2e1 - t1000 / 0.2e1;
t1189 = t808 + (-pkin(5) * t133 + qJ(6) * t134) * t1112;
t774 = pkin(5) * t694 - qJ(6) * t497;
t263 = pkin(5) * t688 - qJ(6) * t451;
t1114 = m(6) / 0.2e1;
t869 = t1051 * t649;
t819 = pkin(3) * t869;
t906 = t1052 * pkin(3);
t633 = t906 + pkin(4);
t927 = t648 * t633;
t585 = t819 + t927;
t568 = qJ(6) + t585;
t926 = t648 * t649;
t584 = -pkin(3) * t926 + t1051 * t633;
t575 = -pkin(5) - t584;
t592 = (t1052 * t648 + t869) * pkin(3);
t1061 = t592 / 0.2e1;
t593 = (t1051 * t1052 - t926) * pkin(3);
t754 = t694 * t1061 + t1218 * t593;
t823 = t1181 * t593 - t290 * t592;
t930 = t585 * t694;
t931 = t584 * t497;
t932 = t575 * t497;
t933 = t568 * t694;
t961 = t526 * mrSges(5,1);
t963 = t524 * mrSges(5,2);
t1180 = (-t1181 * t584 + t290 * t585 + t823) * t1114 + (t1181 * t575 + t290 * t568 + t823) * t1112 - t963 / 0.2e1 - t961 / 0.2e1 + (t754 + t932 / 0.2e1 - t933 / 0.2e1) * mrSges(7,2) + (t754 - t930 / 0.2e1 - t931 / 0.2e1) * mrSges(6,3);
t1178 = m(7) * t103;
t893 = t1041 / 0.2e1;
t1175 = mrSges(6,1) + mrSges(7,1);
t1174 = -mrSges(7,3) + mrSges(6,2);
t1172 = Ifges(6,6) - Ifges(7,6);
t725 = mrSges(5,3) * t753;
t1158 = t688 * mrSges(6,3);
t104 = t653 * pkin(5) - t117;
t1154 = t117 + t104;
t570 = t600 * t653;
t571 = t753 * t653;
t450 = -t1051 * t570 - t571 * t648;
t453 = -t1051 * t571 + t648 * t570;
t1152 = t1036 * t450 + t1037 * t453 + t1173 * t651;
t642 = t653 * mrSges(7,3);
t387 = -t642 + t1204;
t388 = mrSges(6,2) * t653 + t1203;
t1151 = t387 + t388;
t1040 = t651 * pkin(2);
t616 = -pkin(8) * t653 + t1040;
t548 = pkin(7) * t924 + t652 * t616;
t549 = -pkin(7) * t922 + t650 * t616;
t1141 = -t548 * t650 + t549 * t652;
t1139 = -Ifges(4,2) * t650 + t643;
t1138 = -mrSges(4,1) * t652 + mrSges(4,2) * t650;
t390 = -mrSges(6,1) * t653 - t1158;
t1095 = -t390 / 0.2e1;
t391 = mrSges(7,1) * t653 + t1159;
t1137 = t391 / 0.2e1 + t1095;
t1136 = t388 / 0.2e1 + t387 / 0.2e1;
t386 = -mrSges(7,2) * t450 + mrSges(7,3) * t651;
t1100 = t386 / 0.2e1;
t502 = t651 * pkin(3) - pkin(9) * t920 + t548;
t519 = -pkin(9) * t923 + t549;
t332 = t1052 * t502 - t519 * t649;
t283 = pkin(4) * t651 + pkin(10) * t571 + t332;
t333 = t1052 * t519 + t649 * t502;
t297 = pkin(10) * t570 + t333;
t131 = t1051 * t297 + t648 * t283;
t114 = qJ(6) * t651 + t131;
t130 = t1051 * t283 - t648 * t297;
t115 = -t651 * pkin(5) - t130;
t949 = t651 * mrSges(7,1);
t970 = t453 * mrSges(7,2);
t393 = -t949 + t970;
t1134 = (-pkin(5) * t115 + qJ(6) * t114) * t1112 + qJ(6) * t1100 - pkin(5) * t393 / 0.2e1;
t761 = Ifges(5,5) * t718 + Ifges(5,6) * t569 + t1142;
t1132 = (-t649 * mrSges(5,1) - t1052 * mrSges(5,2)) * pkin(3);
t265 = -mrSges(6,1) * t451 + mrSges(6,2) * t688;
t1131 = m(6) * t500 + t265;
t264 = -mrSges(7,1) * t451 - mrSges(7,3) * t688;
t1130 = m(7) * t201 + t264;
t1003 = t118 * mrSges(7,1);
t1004 = t118 * mrSges(6,1);
t1005 = t117 * mrSges(7,3);
t1006 = t117 * mrSges(6,2);
t1129 = -t1003 / 0.2e1 - t1004 / 0.2e1 + t1005 / 0.2e1 - t1006 / 0.2e1;
t1128 = m(7) * t104 - t390 + t391;
t1127 = t1151 - t1203;
t883 = -t1159 / 0.2e1;
t1125 = pkin(5) * t881 + qJ(6) * t883 - t1126;
t1093 = t393 / 0.2e1;
t389 = -mrSges(6,2) * t651 - mrSges(6,3) * t450;
t392 = mrSges(6,1) * t651 - mrSges(6,3) * t453;
t1124 = (t1051 * t130 + t131 * t648) * t908 + (t114 * t631 + t115 * t632) * t1112 + t392 * t813 + t389 * t893 + t632 * t1093 + t631 * t1100;
t1122 = m(6) * t118 + t1151 + t1178;
t1121 = -m(6) * t117 + t1128;
t1009 = Ifges(4,3) * t651;
t1046 = pkin(3) * t649;
t1062 = t585 / 0.2e1;
t1064 = t584 / 0.2e1;
t1111 = m(5) * pkin(3);
t528 = -mrSges(5,2) * t651 + mrSges(5,3) * t570;
t530 = mrSges(5,1) * t651 + mrSges(5,3) * t571;
t815 = t906 / 0.2e1;
t839 = t920 / 0.2e1;
t842 = -t923 / 0.2e1;
t1120 = (t1052 * t332 + t333 * t649) * t1111 / 0.2e1 + (t114 * t568 + t115 * t575) * t1112 + (t130 * t584 + t131 * t585) * t1114 + t530 * t815 + Ifges(4,6) * t842 + Ifges(4,5) * t839 + t528 * t1046 / 0.2e1 + t389 * t1062 + t392 * t1064 + t575 * t1093 + t568 * t1100 - t549 * mrSges(4,2) / 0.2e1 + t548 * mrSges(4,1) / 0.2e1 + t1009 / 0.2e1;
t1118 = 0.2e1 * m(7);
t1116 = 0.2e1 * qJ(6);
t1115 = m(5) / 0.2e1;
t1109 = m(7) * pkin(4);
t1107 = -t115 / 0.2e1;
t1106 = t264 / 0.2e1;
t1105 = t265 / 0.2e1;
t1103 = -t1181 / 0.4e1;
t1102 = t314 / 0.2e1;
t315 = -mrSges(6,1) * t497 + mrSges(6,2) * t694;
t1101 = t315 / 0.2e1;
t1099 = -t387 / 0.2e1;
t1097 = -t388 / 0.2e1;
t1084 = -t450 / 0.2e1;
t1083 = t450 / 0.2e1;
t1078 = t453 / 0.2e1;
t1075 = -t497 / 0.2e1;
t1073 = -t694 / 0.2e1;
t1072 = t694 / 0.2e1;
t529 = -mrSges(5,1) * t653 + t960;
t1070 = t529 / 0.2e1;
t1069 = -t568 / 0.2e1;
t1068 = -t569 / 0.2e1;
t1067 = t570 / 0.2e1;
t1066 = -t571 / 0.2e1;
t1065 = -t575 / 0.2e1;
t1059 = t650 / 0.2e1;
t1058 = t651 / 0.2e1;
t1055 = -t653 / 0.2e1;
t1053 = t653 / 0.2e1;
t1050 = m(7) * t133;
t1049 = m(7) * t143;
t1047 = m(7) * t1181;
t1045 = pkin(4) * t569;
t1044 = pkin(4) * t600;
t1043 = pkin(8) * t650;
t1042 = pkin(8) * t652;
t645 = t650 * pkin(3);
t647 = t653 * pkin(7);
t1035 = Ifges(6,2) + Ifges(7,3);
t1031 = Ifges(3,4) * t651;
t1030 = Ifges(4,4) * t650;
t1024 = Ifges(7,4) * t453;
t1022 = Ifges(5,5) * t571;
t1021 = Ifges(6,5) * t453;
t1015 = Ifges(7,2) * t651;
t1014 = Ifges(5,6) * t570;
t1013 = Ifges(6,6) * t450;
t1011 = Ifges(7,6) * t450;
t1008 = Ifges(5,3) * t651;
t1007 = Ifges(6,3) * t651;
t982 = t324 * mrSges(5,2);
t981 = t325 * mrSges(5,1);
t959 = t584 * mrSges(6,2);
t958 = t584 * mrSges(7,3);
t957 = t585 * mrSges(6,1);
t956 = t585 * mrSges(7,1);
t955 = t592 * mrSges(6,1);
t954 = t592 * mrSges(7,1);
t953 = t593 * mrSges(6,2);
t952 = t593 * mrSges(7,3);
t951 = t600 * mrSges(5,3);
t610 = pkin(3) * t923 + t647;
t501 = -pkin(4) * t570 + t610;
t202 = pkin(5) * t450 - qJ(6) * t453 + t501;
t251 = Ifges(7,5) * t453 + Ifges(7,6) * t651 + Ifges(7,3) * t450;
t253 = Ifges(6,4) * t453 - Ifges(6,2) * t450 + Ifges(6,6) * t651;
t266 = mrSges(7,1) * t450 - mrSges(7,3) * t453;
t267 = mrSges(6,1) * t450 + mrSges(6,2) * t453;
t441 = -Ifges(5,4) * t571 + Ifges(5,2) * t570 + Ifges(5,6) * t651;
t442 = -Ifges(5,1) * t571 + Ifges(5,4) * t570 + Ifges(5,5) * t651;
t461 = -mrSges(5,1) * t718 - t569 * mrSges(5,2);
t462 = -mrSges(5,1) * t570 - mrSges(5,2) * t571;
t527 = t653 * mrSges(5,2) + t714;
t545 = t599 - t907;
t563 = -t653 * Ifges(4,6) + t1139 * t651;
t564 = Ifges(4,6) * t651 + t1139 * t653;
t790 = Ifges(4,1) * t652 - t1030;
t733 = t790 * t651;
t565 = -t653 * Ifges(4,5) + t733;
t566 = Ifges(4,5) * t651 + t653 * t790;
t799 = mrSges(4,1) * t650 + mrSges(4,2) * t652;
t591 = t799 * t653;
t898 = mrSges(4,3) * t924;
t604 = mrSges(4,2) * t653 - t898;
t605 = -t651 * mrSges(4,2) - mrSges(4,3) * t923;
t897 = mrSges(4,3) * t922;
t606 = -mrSges(4,1) * t653 - t897;
t607 = t651 * mrSges(4,1) - mrSges(4,3) * t920;
t713 = t718 / 0.2e1;
t731 = t799 * t646;
t780 = Ifges(4,5) * t652 - Ifges(4,6) * t650;
t734 = t653 * t780;
t840 = t922 / 0.2e1;
t844 = -t924 / 0.2e1;
t921 = t651 * t653;
t7 = t565 * t839 + t566 * t840 + t563 * t842 + t564 * t844 + t441 * t713 + (0.2e1 * Ifges(3,4) * t653 + (Ifges(3,1) - Ifges(3,2)) * t651) * t1053 + t201 * t266 + t202 * t264 + t1152 * t1179 + (-Ifges(5,5) * t569 + Ifges(5,6) * t718 + t651 * t780 - t1031 + t1173 * t688 + t1172 * t451 + (Ifges(3,1) - Ifges(7,2) - Ifges(4,3) - Ifges(5,3) - Ifges(6,3)) * t653) * t1058 + t253 * t1221 - t651 * (Ifges(3,2) * t653 + t1031) / 0.2e1 + m(4) * (pkin(7) ^ 2 * t921 + t545 * t548 + t546 * t549) + m(7) * (t103 * t114 + t104 * t115 + t201 * t202) + m(6) * (t117 * t130 + t118 * t131 + t500 * t501) + m(5) * (t324 * t332 + t325 * t333 + t609 * t610) + t114 * t387 + t131 * t388 + t118 * t389 + t130 * t390 + t115 * t391 + t117 * t392 + t252 * t1084 + t251 * t1085 + t653 * t731 + t103 * t386 + (t1011 + t1015 + t1024 + t1008 + t1014 - t1022 + t1007 - t1013 + t1021 + t734 + t1009) * t1055 + t104 * t393 + t1153 * t1078 + t591 * t646 + t684 * t1066 + t683 * t1067 + t442 * t1068 + t250 * t1083 + t500 * t267 + t501 * t265 + t333 * t527 + t325 * t528 + t332 * t529 + t324 * t530 + t549 * t604 + t546 * t605 + t548 * t606 + t545 * t607 + t609 * t462 + t610 * t461 - pkin(1) * (t651 * mrSges(3,1) + mrSges(3,2) * t653);
t948 = t7 * qJD(1);
t224 = -t1045 + t263;
t625 = pkin(3) * t922;
t203 = t224 + t625;
t517 = t625 - t1045;
t612 = Ifges(4,2) * t652 + t1030;
t657 = t201 * t1182 + t609 * t697 + t104 * t1204 + t325 * t960 - t324 * t714 - t103 * t1159 - t117 * t1203 + t699 * t1068 + t569 * t683 / 0.2e1 - t118 * t1158 + t1237 * t713 + t761 * t1055 + t1258;
t779 = Ifges(4,5) * t650 + Ifges(4,6) * t652;
t841 = -t922 / 0.2e1;
t8 = t657 + m(5) * (t324 * t345 + t325 * t346 + t609 * t625) + t779 * t921 / 0.2e1 + t461 * t625 + t346 * t527 + t345 * t529 + t565 * t844 + t563 * t841 + (-t897 - t606) * t546 + (t898 + t604) * t545 + t1131 * t517 + t1130 * t203 + t1122 * t144 + t1121 * t143 + (-pkin(7) * t1138 + t1059 * t612 - t1235) * t651 ^ 2;
t947 = t8 * qJD(1);
t9 = -t1045 * t1131 + t1121 * t133 + t1122 * t134 + t1130 * t224 + t324 * t527 - t325 * t529 + t657;
t946 = t9 * qJD(1);
t39 = m(7) * (-t103 * t653 - t201 * t688) - t653 * t387 - t688 * t264;
t945 = qJD(1) * t39;
t10 = t263 * t264 + (-t103 * t688 + t104 * t451) * mrSges(7,2) + (m(7) * t263 + t1182) * t201 + (t1128 - t1158) * t118 + (t1127 + t1178) * t117 + t1142 * t1055 + t1258;
t944 = t10 * qJD(1);
t910 = mrSges(5,3) * t1046;
t901 = t584 * t1203;
t900 = t585 * t1158;
t896 = m(7) * t1061;
t892 = -mrSges(7,1) / 0.2e1 - mrSges(6,1) / 0.2e1;
t891 = mrSges(6,2) / 0.2e1 - mrSges(7,3) / 0.2e1;
t873 = t951 / 0.2e1;
t867 = qJ(6) * t1179;
t838 = t497 * t1053;
t834 = t1035 * t694;
t826 = t569 * t910;
t824 = t1181 * t291 - t288 * t290;
t810 = -t325 * t951 / 0.2e1;
t809 = -t642 + t1241;
t803 = (t1218 + t1075) * mrSges(7,2);
t771 = t1047 / 0.2e1 + t1257;
t696 = t114 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t1107 + t130 * mrSges(6,1) / 0.2e1 - t131 * mrSges(6,2) / 0.2e1;
t668 = t1024 / 0.2e1 + t1021 / 0.2e1 + t1015 / 0.2e1 - t1013 / 0.2e1 + t1011 / 0.2e1 + t1007 / 0.2e1 + t696;
t663 = t668 + t1014 / 0.2e1 - t1022 / 0.2e1 + t1008 / 0.2e1 + t332 * mrSges(5,1) / 0.2e1 - t333 * mrSges(5,2) / 0.2e1;
t654 = t103 * t876 + t104 * t879 + t1246 + t663;
t298 = -t1044 + t774;
t286 = t298 + t645;
t510 = mrSges(5,1) * t753 - t600 * mrSges(5,2);
t567 = t645 - t1044;
t764 = t1181 * t144 - t143 * t290;
t655 = (t1114 * t118 + t1136 + t1193) * t291 - (t650 ^ 2 + t652 ^ 2) * mrSges(4,3) * t1039 / 0.2e1 - (t346 / 0.2e1 - t324 / 0.2e1) * t725 + t612 * t841 + (t733 + t565) * t652 / 0.4e1 + t1138 * t1040 / 0.2e1 + t1038 * (t143 * t1072 + t1218 * t144) - (t1139 + 0.2e1 * t613) * t924 / 0.4e1 - t604 * t1043 / 0.2e1 - t606 * t1042 / 0.2e1 + t461 * t645 / 0.2e1 + (t324 * t523 + t325 * t525 + t345 * t524 + t346 * t526 + (t609 * t650 + t634 * t922) * pkin(3)) * t1115 + (t500 * t567 + t517 * t547 + t764) * t1114 + (t201 * t286 + t203 * t284 + t764) * t1112 + pkin(3) * t510 * t840 + t567 * t1105 + t286 * t1106 + t517 * t1101 + t203 * t1102 + t731 / 0.2e1 - t734 / 0.4e1 + (t104 * t1112 - t1114 * t117 + t1137) * t288 + t345 * t873 + t523 * t1070 + t525 * t527 / 0.2e1 - t650 * t563 / 0.4e1;
t1 = t654 + t1120 + t810 - t655;
t662 = -t634 * t721 + t1036 * (-t497 ^ 2 + t694 ^ 2) + t1038 * (-t1181 * t694 - t290 * t497) + t1224;
t678 = -t634 * mrSges(5,1) - t722 + t724;
t732 = t753 ^ 2;
t752 = -t650 * t612 / 0.2e1 + t1235;
t17 = m(7) * (t284 * t286 + t824) + m(6) * (t547 * t567 + t824) + m(5) * (t523 * t524 + t525 * t526 + t634 * t645) + (-t1028 + (t526 + t523) * mrSges(5,3) + t678) * t600 - (-t1038 * t291 + t834) * t497 + (t1038 * t288 + t835) * t694 + t286 * t314 + Ifges(5,4) * t732 + t662 + t567 * t315 - pkin(2) * t799 + t790 * t1059 + t1139 * t1057 + t510 * t645 + t752 + (t524 - t525) * t725;
t770 = -t1 * qJD(1) + t17 * qJD(2);
t18 = -(-t1038 * t290 + t834) * t497 + (t1038 * t1181 + t835) * t694 + (-pkin(4) * t315 + t678) * t600 - m(6) * t1044 * t547 + (-t600 ^ 2 + t732) * Ifges(5,4) + t662 - t1261 * t298;
t765 = t1181 * t134 - t133 * t290;
t660 = -m(6) * (-t117 * t1181 + t118 * t290 + (-t500 * t600 - t547 * t569) * pkin(4) + t765) / 0.2e1 + (t103 * t290 + t104 * t1181 + t201 * t298 + t224 * t284 + t765) * t1113 - t224 * t314 / 0.2e1 + t390 * t1104 - t1181 * t391 / 0.2e1 + t290 * t1099 + t290 * t1097 - t298 * t264 / 0.2e1 - t524 * t527 / 0.2e1 + t526 * t1070 + t1045 * t1101 + t1044 * t1105 + t1038 * (t133 * t1073 + t134 * t1075);
t3 = t654 + t660 + t1124;
t769 = -t3 * qJD(1) + t18 * qJD(2);
t21 = t1181 * t965 + t290 * t1202 - t774 * t314 + (-mrSges(6,3) * t1181 - pkin(5) * t1048 - t1036 * t694) * t694 - (t290 * mrSges(7,2) - qJ(6) * t1048 - t1036 * t497 + (-t1035 + t1037) * t694) * t497 - t1224;
t656 = t1133 * t694 / 0.4e1 + t1260 * t688 / 0.4e1 + t1223 - ((-t103 + t118) * t1112 + t1097 + t1099) * t290 + t1038 * (t1081 * t1181 - t1262) + t1236 - t1254 * t653 / 0.4e1 + (t1112 * t1154 + t1137) * t1181 + t1135 * t1217 + t1230 * t1219 + t1253 * t1220 + t1231 * t1222 + t774 * t1106 + t263 * t1102 + (t774 * t201 + t284 * t263) * t1112 + t1154 * t1242 + t118 * t876;
t5 = -t656 + t668 + t1134;
t768 = -t5 * qJD(1) - t21 * qJD(2);
t681 = (-t1181 * t653 - t201 * t694 - t284 * t688) * t1112 + t314 * t1081 + t264 * t1073;
t762 = m(7) * t1107 + t949 / 0.2e1;
t33 = (-t838 - t453 / 0.2e1) * mrSges(7,2) + t681 + t762;
t74 = t1261 * t694;
t766 = qJD(1) * t33 + qJD(2) * t74;
t760 = pkin(5) * t879 + qJ(6) * t1210 + t1123;
t757 = m(7) * (-pkin(5) * t143 + qJ(6) * t144);
t756 = m(7) * (-pkin(5) * t288 + qJ(6) * t291);
t755 = -t1065 * t451 + t1069 * t688;
t701 = t714 * t906;
t658 = t808 - t701 / 0.2e1 - t982 / 0.2e1 - t981 / 0.2e1 + t826 / 0.2e1 + (-t117 * t592 + t118 * t593 - t133 * t584 + t134 * t585) * t1114 + (t103 * t593 + t104 * t592 + t133 * t575 + t134 * t568) * t1112 + t527 * t815 + t568 * t883 + t575 * t1241 - t901 / 0.2e1 - t900 / 0.2e1 + t391 * t1061 + t592 * t1095 - t529 * t1046 / 0.2e1 + t1151 * t593 / 0.2e1;
t12 = t658 + t1199;
t20 = t695 + t1180 + t1200;
t91 = -t1174 * t593 - t1175 * t592 + t1132 + m(6) * (-t584 * t592 + t585 * t593) + m(7) * (t568 * t593 + t575 * t592);
t720 = t12 * qJD(1) + t20 * qJD(2) + t91 * qJD(3);
t677 = (t103 * t584 + t104 * t585 + t117 * t568 + t118 * t575) * t1112 - t1126;
t14 = (-t1158 / 0.2e1 + t1137) * t585 + (-t1203 / 0.2e1 + t1136) * t584 - t757 / 0.2e1 + (pkin(5) * t1221 + t755 + t867) * mrSges(7,2) + t677 + t1175 * (t143 / 0.2e1 - t118 / 0.2e1) + t1174 * (t144 / 0.2e1 - t117 / 0.2e1) + t1126;
t671 = (-(-t568 + t585) * t290 + (t575 + t584) * t1181) * t1112 + t716;
t707 = (t1069 + t1062) * t694 - (t1065 - t584 / 0.2e1) * t497;
t23 = (t1156 / 0.2e1 + pkin(5) * t1218 + t707) * mrSges(7,2) - t892 * t288 + t891 * t291 - t756 / 0.2e1 + t671 - t1123;
t709 = -t956 - t957 + t958 - t959;
t79 = m(7) * (t568 * t584 + t575 * t585) + t709;
t719 = t14 * qJD(1) + t23 * qJD(2) + t79 * qJD(3);
t680 = ((-qJ(6) - t568) * t653 + t118) * t1112 + t809;
t42 = t881 - t1049 / 0.2e1 + t680;
t533 = m(7) * t568 + mrSges(7,3);
t69 = t803 + (t288 / 0.4e1 + t1103) * t1118;
t717 = -qJD(1) * t42 + qJD(2) * t69 - qJD(3) * t533;
t661 = (t631 * t117 + t632 * t118 + (t103 * t1051 + t104 * t648) * pkin(4)) * t1112 + t632 * t1241 + t390 * t894 + t391 * t893 + t688 * t820 - t1126 + t1151 * t813 + t1129 - t1228;
t16 = mrSges(7,2) * t867 + pkin(5) * t1241 + t1126 - t1189 + t661;
t670 = t1112 * t1250 + t695 + t760;
t25 = t670 + t1198;
t640 = mrSges(7,3) * t905;
t700 = mrSges(6,2) * t905 + t1041 * t1175 - t640;
t406 = -(t1051 * t631 + t632 * t648) * t1109 + t700;
t664 = t640 / 0.2e1 + (t631 * t584 + t632 * t585 + (t1051 * t568 + t575 * t648) * pkin(4)) * t1112 - t959 / 0.2e1 + t958 / 0.2e1 - t957 / 0.2e1 - t956 / 0.2e1 + mrSges(6,2) * t814 + t1175 * t894;
t673 = (-pkin(5) * t592 + qJ(6) * t593) * t1112 - t955 / 0.2e1 - t954 / 0.2e1 - t953 / 0.2e1 + t952 / 0.2e1;
t66 = t664 - t673;
t708 = t16 * qJD(1) - t25 * qJD(2) + t66 * qJD(3) - t406 * qJD(4);
t690 = -mrSges(7,3) + (t819 + t1116 + (pkin(4) + t633) * t648) * t1113;
t407 = t896 + t690;
t679 = ((-qJ(6) - t631) * t653 + t118) * t1112 + t809;
t44 = t881 - t1050 / 0.2e1 + t679;
t608 = m(7) * t631 + mrSges(7,3);
t68 = t803 + (t1181 / 0.4e1 + t1103) * t1118;
t706 = -qJD(1) * t44 + qJD(2) * t68 + qJD(3) * t407 - qJD(4) * t608;
t309 = t1242 + t879;
t438 = -mrSges(7,3) + (t585 / 0.4e1 - t819 / 0.4e1 - t927 / 0.4e1 - qJ(6) / 0.2e1) * t1118;
t46 = -t642 + (t870 / 0.4e1 - t919 / 0.2e1 + t929 / 0.4e1 - t118 / 0.4e1) * t1118;
t685 = qJD(1) * t46 + qJD(2) * t309 - qJD(3) * t438 + qJD(4) * t641 + t1194;
t659 = t104 * t1242 + t325 * t873 + t1236 - t1246 + t663;
t587 = mrSges(7,3) + (qJ(6) + 0.2e1 * t893) * m(7);
t439 = mrSges(7,3) + (t1116 + t585) * t1112 + m(7) * t1062;
t408 = t896 - t690;
t145 = t1257 + t1047;
t71 = t1112 * t288 + t771;
t70 = m(7) * t1104 + t771;
t67 = t664 + t673;
t45 = t387 + 0.2e1 * t1193;
t43 = t1241 + t1050 / 0.2e1 + t679;
t41 = t1241 + t1049 / 0.2e1 + t680;
t32 = -mrSges(7,2) * t838 + t970 / 0.2e1 + t681 - t762;
t24 = t670 - t1198;
t22 = t756 / 0.2e1 + t707 * mrSges(7,2) + t671 + t760 + t806;
t19 = t1181 * t892 - t891 * t290 + t1180 - t1200 + t759;
t15 = t661 + t1125 + t1189;
t13 = t807 + t677 + t755 * mrSges(7,2) + t757 / 0.2e1 + t391 * t1062 - (t390 + t1158) * t585 / 0.2e1 + t1127 * t1064 + t1125 + t1129;
t11 = t658 + t761 - t1199;
t6 = t656 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t651 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t450 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t453 + t696 + t1134;
t4 = -t660 + t659 + t810 + t1124;
t2 = t1120 + t659 + t655;
t26 = [qJD(2) * t7 + qJD(3) * t8 + qJD(4) * t9 + qJD(5) * t10 + qJD(6) * t39, t2 * qJD(3) + t4 * qJD(4) + t6 * qJD(5) + t32 * qJD(6) + t948 + (t202 * t314 + t284 * t266 + (-Ifges(5,5) * t600 - Ifges(5,6) * t753 + t1172 * t497 + t1173 * t694 + t779) * t1058 + t115 * t1157 + (Ifges(3,5) + (-mrSges(3,1) + t1138) * pkin(7) + t752) * t653 + t114 * t1202 + t253 * t1218 + t1181 * t389 + t1181 * t386 - t607 * t1043 - t130 * t965 + t131 * t1209 + t332 * t951 - t333 * t725 - t753 * t441 / 0.2e1 + 0.2e1 * (t332 * t524 + t333 * t526 + t610 * t634) * t1115 + 0.2e1 * (t1181 * t131 + t130 * t290 + t501 * t547) * t1114 + 0.2e1 * (t114 * t1181 - t115 * t290 + t202 * t284) * t1112 + t290 * t392 - t290 * t393 + m(4) * (-pkin(2) * t647 + pkin(8) * t1141) + t1141 * mrSges(4,3) + t1149 * t1078 + t1152 * t1072 + mrSges(3,2) * t646 + t605 * t1042 + t564 * t1057 + t566 * t1059 + t705 * t1066 + t704 * t1067 + t251 * t1075 + t776 * t1083 + t781 * t1084 + t501 * t315 + t526 * t528 + t524 * t530 + t547 * t267 - pkin(2) * t591 - t600 * t442 / 0.2e1 + t610 * t510 + t634 * t462 - Ifges(3,6) * t651) * qJD(2), t947 + t2 * qJD(2) + (-t701 - t998 - t996 + t995 - t997 + m(7) * (t143 * t575 + t144 * t568) + m(6) * (-t143 * t584 + t144 * t585) + t826 + t761 + (t1052 * t345 + t346 * t649) * t1111 - t568 * t1159 + t575 * t1204 - t901 - t900 - t545 * mrSges(4,2) - t546 * mrSges(4,1) + t980 - t979 - Ifges(4,5) * t924 - Ifges(4,6) * t922) * qJD(3) + t11 * qJD(4) + t13 * qJD(5) + t41 * qJD(6), t946 + t4 * qJD(2) + t11 * qJD(3) + (-t982 - t981 - t1000 + t999 - t1001 - t1002 + m(7) * (t133 * t632 + t134 * t631) - t827 - t451 * t822 + t761 + (-t1051 * t133 + t134 * t648) * t1110 - t899 + t632 * t1204) * qJD(4) + t15 * qJD(5) + t43 * qJD(6), t944 + t6 * qJD(2) + t13 * qJD(3) + t15 * qJD(4) + (m(7) * (-pkin(5) * t118 + qJ(6) * t117) + t1005 - t1004 - t1006 - t1003 + t775 * mrSges(7,2) + t1142) * qJD(5) + t45 * qJD(6), qJD(2) * t32 + qJD(3) * t41 + qJD(4) * t43 + qJD(5) * t45 + t945; -qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t5 + qJD(6) * t33 - t948, qJD(3) * t17 + qJD(4) * t18 - qJD(5) * t21 + qJD(6) * t74 (-t986 + t985 - t991 - t992 + t725 * t906 + m(7) * (t288 * t575 + t291 * t568) + m(6) * (-t288 * t584 + t291 * t585) + t600 * t910 + t759 + (t1052 * t523 + t525 * t649) * t1111 + t964 - t962 + t780 + t1138 * pkin(8) + (-t931 - t930) * mrSges(6,3) + (-t933 + t932) * mrSges(7,2)) * qJD(3) + t19 * qJD(4) + t22 * qJD(5) + t71 * qJD(6) + t770, t19 * qJD(3) + (m(7) * t1248 + t1110 * t1249 - t631 * t1157 + t632 * t1202 - t694 * t909 + t1259 + t759 - t802 - t961 - t963) * qJD(4) + t24 * qJD(5) + t70 * qJD(6) + t769, t22 * qJD(3) + t24 * qJD(4) + t145 * qJD(6) + t768 + (m(7) * t1250 + (-mrSges(7,2) * qJ(6) - t1172) * t694 - (mrSges(7,2) * pkin(5) - t1173) * t497 + t1259) * qJD(5), qJD(3) * t71 + qJD(4) * t70 + qJD(5) * t145 + t766; qJD(2) * t1 + qJD(4) * t12 + qJD(5) * t14 + qJD(6) * t42 - t947, qJD(4) * t20 + qJD(5) * t23 - qJD(6) * t69 - t770, qJD(4) * t91 + qJD(5) * t79 + qJD(6) * t533 ((-t1051 * t592 + t593 * t648) * t1110 - t953 - t955 - t954 + t952 + m(7) * (t592 * t632 + t593 * t631) + t1132) * qJD(4) + t67 * qJD(5) + t408 * qJD(6) + t720, t67 * qJD(4) + (m(7) * (-pkin(5) * t585 + qJ(6) * t584) + t709) * qJD(5) + t439 * qJD(6) + t719, qJD(4) * t408 + qJD(5) * t439 - t717; qJD(2) * t3 - qJD(3) * t12 + qJD(5) * t16 + qJD(6) * t44 - t946, -qJD(3) * t20 - qJD(5) * t25 - qJD(6) * t68 - t769, qJD(5) * t66 - qJD(6) * t407 - t720, -qJD(5) * t406 + qJD(6) * t608 ((-pkin(5) * t648 + qJ(6) * t1051) * t1109 - t700) * qJD(5) + t587 * qJD(6) + t708, qJD(5) * t587 - t706; qJD(2) * t5 - qJD(3) * t14 - qJD(4) * t16 + qJD(6) * t46 - t944, -qJD(3) * t23 + qJD(4) * t25 + qJD(6) * t309 - t768, -qJD(4) * t66 - qJD(6) * t438 - t719, -t708 + t1191, t1191, t685; -qJD(2) * t33 - qJD(3) * t42 - qJD(4) * t44 - qJD(5) * t46 - t945, qJD(3) * t69 + qJD(4) * t68 - qJD(5) * t309 - t766, qJD(4) * t407 + qJD(5) * t438 + t717, t706 - t1194, -t685, 0;];
Cq  = t26;
