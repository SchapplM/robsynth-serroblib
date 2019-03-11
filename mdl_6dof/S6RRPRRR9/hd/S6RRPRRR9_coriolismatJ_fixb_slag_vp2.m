% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:41
% EndTime: 2019-03-09 14:09:45
% DurationCPUTime: 38.87s
% Computational Cost: add. (96591->1203), mult. (226435->1644), div. (0->0), fcn. (266608->12), ass. (0->653)
t1047 = cos(qJ(5));
t1046 = sin(qJ(4));
t1048 = cos(qJ(4));
t678 = sin(pkin(12));
t680 = cos(pkin(12));
t681 = cos(pkin(6));
t679 = sin(pkin(6));
t684 = sin(qJ(2));
t909 = t679 * t684;
t631 = -t678 * t909 + t680 * t681;
t906 = t680 * t684;
t910 = t678 * t681;
t741 = t679 * t906 + t910;
t559 = t1046 * t631 + t1048 * t741;
t683 = sin(qJ(5));
t710 = t1046 * t741 - t1048 * t631;
t428 = t1047 * t710 + t559 * t683;
t1272 = t428 / 0.2e1;
t685 = cos(qJ(6));
t1051 = -t685 / 0.2e1;
t1271 = pkin(11) * t428;
t1248 = t1047 * t559 - t683 * t710;
t682 = sin(qJ(6));
t1015 = Ifges(7,4) * t682;
t656 = Ifges(7,1) * t685 - t1015;
t1217 = Ifges(7,5) * t1248 - t428 * t656;
t1270 = t1217 / 0.4e1;
t1114 = -t1248 / 0.2e1;
t1269 = mrSges(6,1) * t1114;
t644 = -t1046 * t678 + t1048 * t680;
t645 = -t1046 * t680 - t1048 * t678;
t592 = t1047 * t644 + t645 * t683;
t738 = -t1047 * t645 + t683 * t644;
t916 = t738 * t682;
t475 = mrSges(7,2) * t592 - mrSges(7,3) * t916;
t1268 = t475 * t1051;
t1105 = t1248 / 0.2e1;
t673 = Ifges(7,4) * t685;
t655 = Ifges(7,1) * t682 + t673;
t888 = t685 * t655;
t653 = Ifges(7,2) * t685 + t1015;
t897 = t682 * t653;
t1183 = t897 / 0.4e1 - t888 / 0.4e1;
t1267 = t1183 * t428;
t652 = Ifges(7,5) * t682 + Ifges(7,6) * t685;
t1190 = t738 * t652;
t1031 = pkin(9) + qJ(3);
t815 = t1031 * t680;
t816 = t1031 * t678;
t602 = -t1046 * t815 - t1048 * t816;
t552 = t645 * pkin(10) + t602;
t603 = -t1046 * t816 + t1048 * t815;
t711 = t644 * pkin(10) + t603;
t405 = -t1047 * t552 + t683 * t711;
t1257 = t405 * mrSges(6,2);
t1266 = t1190 / 0.2e1 + t1257;
t1194 = t1248 * mrSges(6,1);
t1056 = -t682 / 0.2e1;
t947 = t685 * mrSges(7,3);
t1215 = mrSges(7,1) * t1248 + t428 * t947;
t952 = t682 * mrSges(7,3);
t1216 = -mrSges(7,2) * t1248 + t428 * t952;
t1237 = t428 * mrSges(6,2);
t719 = t1237 / 0.2e1 + t1216 * t1056 + t1215 * t1051;
t1265 = t1194 / 0.2e1 - t719;
t1238 = Ifges(6,5) * t592;
t1264 = t1238 / 0.2e1 + t1257 / 0.2e1;
t1199 = Ifges(6,6) * t738;
t1263 = -t1183 * t592 + t1190 / 0.4e1 - t1199 / 0.2e1;
t669 = -pkin(3) * t680 - pkin(2);
t619 = -pkin(4) * t644 + t669;
t773 = -t1237 + t1194;
t915 = t738 * t685;
t1262 = t915 * t1270 + t619 * t773 / 0.2e1;
t1130 = -t1215 / 0.2e1;
t1170 = t1047 * t711 + t552 * t683;
t408 = -pkin(5) * t592 - pkin(11) * t738 + t619;
t252 = t1170 * t685 + t408 * t682;
t1134 = -t252 / 0.2e1;
t1178 = -Ifges(7,2) * t682 + t673;
t1214 = Ifges(7,6) * t738 + t1178 * t592;
t1218 = Ifges(7,5) * t738 + t592 * t656;
t251 = -t1170 * t682 + t408 * t685;
t686 = cos(qJ(2));
t908 = t679 * t686;
t375 = -t1248 * t682 - t685 * t908;
t376 = t1248 * t685 - t682 * t908;
t1032 = t631 * pkin(3);
t1039 = pkin(1) * t686;
t664 = pkin(8) * t909;
t635 = t1039 * t681 - t664;
t454 = -t681 * pkin(2) + pkin(4) * t710 - t1032 - t635;
t1236 = t592 * mrSges(6,2);
t963 = t738 * mrSges(6,1);
t802 = t1236 + t963;
t1261 = t251 * t1130 - t1218 * t376 / 0.4e1 - t1214 * t375 / 0.4e1 + t1216 * t1134 - t802 * t454 / 0.2e1;
t1155 = m(7) / 0.4e1;
t1260 = 0.2e1 * t1155;
t1118 = -t405 / 0.2e1;
t1091 = -t559 / 0.2e1;
t1054 = t682 / 0.2e1;
t1050 = t685 / 0.2e1;
t1259 = -t710 / 0.2e1;
t1089 = -t1236 / 0.2e1;
t676 = t682 ^ 2;
t677 = t685 ^ 2;
t821 = t677 / 0.2e1 + t676 / 0.2e1;
t1193 = t821 * mrSges(7,3);
t904 = t681 * t684;
t636 = pkin(1) * t904 + pkin(8) * t908;
t620 = qJ(3) * t681 + t636;
t621 = (-pkin(2) * t686 - qJ(3) * t684 - pkin(1)) * t679;
t534 = -t678 * t620 + t680 * t621;
t483 = -pkin(3) * t908 - pkin(9) * t741 + t534;
t535 = t680 * t620 + t678 * t621;
t499 = pkin(9) * t631 + t535;
t309 = -t1046 * t499 + t1048 * t483;
t259 = -t559 * pkin(10) + t309;
t238 = -pkin(4) * t908 + t259;
t310 = t1046 * t483 + t1048 * t499;
t260 = -pkin(10) * t710 + t310;
t850 = t1047 * t260;
t121 = t683 * t238 + t850;
t936 = t121 * t405;
t1256 = t405 * t682;
t1255 = t405 * t683;
t1254 = t405 * t685;
t425 = Ifges(6,4) * t428;
t768 = -Ifges(6,2) * t1248 - t425;
t1049 = t685 / 0.4e1;
t1053 = t682 / 0.4e1;
t1066 = t655 / 0.4e1;
t1070 = -t653 / 0.4e1;
t1010 = Ifges(7,5) * t428;
t370 = Ifges(7,4) * t375;
t162 = Ifges(7,1) * t376 + t1010 + t370;
t894 = t685 * t162;
t1016 = Ifges(7,4) * t376;
t974 = t428 * Ifges(7,6);
t161 = Ifges(7,2) * t375 + t1016 + t974;
t903 = t682 * t161;
t1184 = -t894 / 0.4e1 + t903 / 0.4e1;
t672 = Ifges(7,5) * t685;
t999 = Ifges(7,6) * t682;
t1188 = t672 - t999;
t1226 = t428 * t1188;
t1024 = Ifges(7,1) * t375;
t769 = -t1016 + t1024;
t1007 = Ifges(7,2) * t376;
t799 = t370 - t1007;
t1251 = t769 * t1053 + t799 * t1049 + t1226 / 0.4e1 + (t656 / 0.4e1 + t1070) * t376 + (t1066 + t1178 / 0.4e1) * t375 - t1184;
t764 = -t1199 + t1238;
t1025 = Ifges(6,1) * t592;
t1204 = Ifges(6,4) * t738;
t1250 = -t1204 + t1025;
t1200 = Ifges(6,6) * t1248;
t1239 = Ifges(6,5) * t428;
t887 = -t1239 - t1200;
t1191 = t1248 * t652;
t1249 = t1191 / 0.4e1 - t1200 / 0.2e1 - t1239 / 0.2e1;
t342 = -Ifges(7,3) * t592 + t1188 * t738;
t487 = Ifges(6,2) * t592 + t1204;
t1174 = t342 / 0.2e1 - t487 / 0.2e1;
t958 = t619 * mrSges(6,1);
t1247 = t1218 * t1050 + t1214 * t1056 + t958 + t1174 - t1204 / 0.2e1;
t1120 = t376 / 0.2e1;
t1122 = t375 / 0.2e1;
t896 = t683 * t260;
t120 = t1047 * t238 - t896;
t113 = pkin(5) * t908 - t120;
t1213 = Ifges(7,6) * t1248 - t1178 * t428;
t948 = t685 * mrSges(7,2);
t953 = t682 * mrSges(7,1);
t651 = t948 + t953;
t1235 = t651 * t428;
t1246 = t1217 * t1120 + t1213 * t1122 - t113 * t1235;
t975 = t428 * mrSges(6,3);
t385 = mrSges(6,2) * t908 - t975;
t1244 = t385 / 0.2e1;
t1112 = -t428 / 0.2e1;
t1243 = t592 / 0.2e1;
t1241 = t1191 / 0.2e1;
t1196 = Ifges(7,3) * t1248;
t1240 = t1196 / 0.2e1;
t1234 = t651 * t592;
t746 = -t999 / 0.2e1 + t672 / 0.2e1;
t729 = t746 * t428;
t1233 = t746 * t592;
t1232 = t680 * t631 + t678 * t741;
t983 = t375 * mrSges(7,3);
t242 = -mrSges(7,2) * t428 + t983;
t982 = t376 * mrSges(7,3);
t243 = mrSges(7,1) * t428 - t982;
t1231 = t242 * t1050 + t243 * t1056;
t1037 = pkin(4) * t645;
t1034 = pkin(5) * t738;
t489 = -pkin(11) * t592 + t1034;
t432 = -t1037 + t489;
t269 = t432 * t685 + t1256;
t270 = t432 * t682 - t1254;
t1230 = -t269 * t682 + t270 * t685;
t127 = t1047 * t259 - t896;
t1038 = pkin(4) * t559;
t280 = pkin(5) * t1248 + t1271;
t192 = t1038 + t280;
t82 = -t127 * t682 + t192 * t685;
t83 = t127 * t685 + t192 * t682;
t1229 = -t682 * t82 + t685 * t83;
t966 = t559 * mrSges(5,3);
t524 = -mrSges(5,1) * t908 - t966;
t1228 = t966 + t524;
t1187 = t676 + t677;
t1227 = t1187 * t1260;
t876 = t1047 * pkin(4);
t671 = -t876 - pkin(5);
t1157 = 0.2e1 * t671;
t1225 = t1170 * t1157;
t992 = t120 * mrSges(6,3);
t1224 = -t992 - t903 / 0.2e1;
t587 = Ifges(6,4) * t592;
t488 = Ifges(6,1) * t738 + t587;
t1222 = -t488 / 0.4e1 + mrSges(6,3) * t1118;
t1221 = t1050 * t243 + t1054 * t242;
t1211 = t113 / 0.2e1;
t1088 = -t738 / 0.2e1;
t1087 = t738 / 0.2e1;
t840 = -t908 / 0.2e1;
t1208 = mrSges(7,1) * t738;
t1206 = mrSges(7,2) * t738;
t1205 = mrSges(6,3) * t121;
t1019 = Ifges(6,4) * t1248;
t1195 = Ifges(7,3) * t738;
t772 = mrSges(7,1) * t685 - mrSges(7,2) * t682;
t1192 = -mrSges(6,1) - t772;
t1189 = -t897 / 0.2e1 + t888 / 0.2e1;
t1186 = Ifges(5,5) * t644 + Ifges(5,6) * t645 + t764;
t1185 = t799 + t162;
t1181 = -Ifges(5,5) * t710 - Ifges(5,6) * t559 + t887;
t478 = -mrSges(7,1) * t592 - mrSges(7,3) * t915;
t1180 = m(7) * t251 + t478;
t486 = -mrSges(6,1) * t592 + mrSges(6,2) * t738;
t1179 = m(6) * t619 + t486;
t1177 = t1155 * t1157;
t1176 = mrSges(4,1) * t678 + mrSges(4,2) * t680;
t1072 = -t772 / 0.2e1;
t1175 = t1072 - mrSges(6,1) / 0.2e1;
t160 = Ifges(7,5) * t376 + Ifges(7,6) * t375 + t428 * Ifges(7,3);
t261 = -Ifges(6,2) * t428 - Ifges(6,6) * t908 + t1019;
t826 = -t261 / 0.2e1 + t160 / 0.2e1;
t209 = -mrSges(7,1) * t375 + mrSges(7,2) * t376;
t971 = t1248 * mrSges(6,3);
t386 = -mrSges(6,1) * t908 - t971;
t1173 = -t209 / 0.2e1 + t386 / 0.2e1;
t1008 = Ifges(7,5) * t592;
t348 = t656 * t738 - t1008;
t890 = t685 * t348;
t1000 = Ifges(7,6) * t592;
t345 = t1178 * t738 - t1000;
t900 = t682 * t345;
t1172 = t488 / 0.2e1 - t900 / 0.2e1 + t890 / 0.2e1 + t587 / 0.2e1;
t1115 = t1170 / 0.2e1;
t1171 = mrSges(6,3) * t1115 + t487 / 0.4e1 - t342 / 0.4e1;
t857 = -t947 / 0.2e1;
t858 = t952 / 0.2e1;
t1169 = t375 * t858 + t376 * t857 - t1221;
t1168 = -t890 / 0.4e1 + t900 / 0.4e1 + t1222;
t1167 = qJD(4) * (-mrSges(6,3) * t876 + t1189);
t611 = t645 * t908;
t612 = t644 * t908;
t517 = t1047 * t612 + t683 * t611;
t496 = -t517 * t682 + t685 * t909;
t497 = t517 * t685 + t682 * t909;
t516 = -t1047 * t611 + t612 * t683;
t228 = Ifges(7,4) * t497 + Ifges(7,2) * t496 + Ifges(7,6) * t516;
t229 = Ifges(7,1) * t497 + Ifges(7,4) * t496 + Ifges(7,5) * t516;
t1013 = Ifges(6,5) * t517;
t1068 = t653 / 0.4e1;
t1148 = -mrSges(6,2) / 0.2e1;
t634 = (pkin(2) * t684 - qJ(3) * t686) * t679;
t572 = t680 * t634 - t678 * t635;
t905 = t680 * t686;
t518 = (pkin(3) * t684 - pkin(9) * t905) * t679 + t572;
t573 = t678 * t634 + t680 * t635;
t870 = t678 * t908;
t525 = -pkin(9) * t870 + t573;
t355 = -t1046 * t525 + t1048 * t518;
t312 = pkin(4) * t909 - t612 * pkin(10) + t355;
t356 = t1046 * t518 + t1048 * t525;
t317 = pkin(10) * t611 + t356;
t168 = t1047 * t312 - t683 * t317;
t158 = -pkin(5) * t909 - t168;
t169 = t1047 * t317 + t683 * t312;
t703 = t1013 / 0.2e1 + t158 * t1072 + t168 * mrSges(6,1) / 0.2e1 + t169 * t1148 + t496 * t1068 + t497 * t1066;
t841 = t909 / 0.2e1;
t159 = pkin(11) * t909 + t169;
t610 = pkin(3) * t870 + t636;
t522 = -pkin(4) * t611 + t610;
t293 = pkin(5) * t516 - pkin(11) * t517 + t522;
t100 = t159 * t685 + t293 * t682;
t940 = t100 * t685;
t99 = -t159 * t682 + t293 * t685;
t942 = t99 * t682;
t1166 = Ifges(6,3) * t841 + t1049 * t228 + t1053 * t229 + t703 + (-t942 / 0.2e1 + t940 / 0.2e1) * mrSges(7,3);
t1165 = Ifges(6,6) * t516 / 0.2e1 - Ifges(5,6) * t611 / 0.2e1 - Ifges(5,5) * t612 / 0.2e1;
t881 = 0.2e1 * m(6);
t806 = t881 / 0.4e1;
t788 = pkin(4) * t806;
t1164 = -t1047 * t788 + t1177;
t1036 = pkin(4) * t683;
t670 = pkin(11) + t1036;
t1163 = t670 * t1227 + t683 * t788 + t1193;
t1161 = -0.2e1 * pkin(5);
t1159 = t680 ^ 2;
t1158 = -0.2e1 * t1170;
t1156 = m(7) / 0.2e1;
t1154 = -pkin(5) / 0.2e1;
t1153 = m(7) * pkin(4);
t1152 = m(7) * pkin(5);
t1150 = mrSges(7,1) / 0.2e1;
t1149 = -mrSges(5,2) / 0.2e1;
t1147 = -mrSges(7,2) / 0.2e1;
t1146 = -mrSges(7,3) / 0.2e1;
t1145 = -Ifges(7,5) / 0.2e1;
t114 = -pkin(11) * t908 + t121;
t181 = t428 * pkin(5) - pkin(11) * t1248 + t454;
t73 = t114 * t685 + t181 * t682;
t1144 = -t73 / 0.2e1;
t1143 = m(7) * t82;
t1142 = m(7) * t83;
t1141 = m(7) * t99;
t1140 = t100 / 0.2e1;
t1138 = t1213 / 0.4e1;
t1136 = -t242 / 0.2e1;
t1135 = t242 / 0.2e1;
t1132 = -t269 / 0.2e1;
t1131 = t270 / 0.2e1;
t287 = t489 * t685 + t1256;
t1129 = -t287 / 0.2e1;
t313 = -mrSges(7,1) * t496 + mrSges(7,2) * t497;
t1128 = -t313 / 0.2e1;
t1125 = -t345 / 0.4e1;
t1124 = t348 / 0.4e1;
t1117 = t405 / 0.2e1;
t1116 = -t1170 / 0.2e1;
t1111 = -t428 / 0.4e1;
t1109 = t428 / 0.4e1;
t1106 = -t1248 / 0.4e1;
t1102 = -t1234 / 0.2e1;
t468 = t651 * t738;
t1101 = -t468 / 0.2e1;
t1100 = -t475 / 0.2e1;
t477 = -t592 * t947 + t1208;
t1099 = -t477 / 0.2e1;
t1098 = -t478 / 0.2e1;
t1094 = t496 / 0.2e1;
t1093 = t497 / 0.2e1;
t1092 = t517 / 0.2e1;
t1090 = t559 / 0.2e1;
t1085 = t592 / 0.4e1;
t1084 = -t592 / 0.2e1;
t1083 = -t592 / 0.4e1;
t1081 = -t738 / 0.4e1;
t1080 = t738 / 0.4e1;
t1078 = t611 / 0.2e1;
t1077 = t612 / 0.2e1;
t1075 = t644 / 0.2e1;
t1074 = -t645 / 0.2e1;
t1073 = t645 / 0.2e1;
t1071 = t652 / 0.4e1;
t1069 = t653 / 0.2e1;
t1067 = -t655 / 0.2e1;
t1065 = -t670 / 0.2e1;
t1064 = t671 / 0.2e1;
t1061 = -t678 / 0.2e1;
t1060 = t678 / 0.2e1;
t1058 = t680 / 0.2e1;
t1057 = t681 / 0.2e1;
t1052 = t684 / 0.2e1;
t1045 = m(6) * t1170;
t1043 = m(7) * t158;
t1041 = m(7) * t269;
t1040 = m(7) * t270;
t208 = mrSges(7,1) * t376 + mrSges(7,2) * t375;
t1035 = pkin(5) * t208;
t1033 = pkin(5) * t651;
t1028 = mrSges(4,3) * t678;
t1027 = mrSges(5,3) * t644;
t1026 = mrSges(5,3) * t645;
t1022 = Ifges(4,4) * t678;
t1021 = Ifges(5,4) * t559;
t1020 = Ifges(5,4) * t645;
t1009 = Ifges(7,5) * t497;
t1001 = Ifges(7,6) * t496;
t996 = Ifges(7,3) * t516;
t993 = t120 * mrSges(6,2);
t126 = t259 * t683 + t850;
t991 = t126 * mrSges(6,1);
t990 = t127 * mrSges(6,2);
t989 = t251 * mrSges(7,3);
t89 = -t120 * t682 + t280 * t685;
t988 = t251 * t89;
t987 = t252 * mrSges(7,3);
t90 = t120 * t685 + t280 * t682;
t986 = t252 * t90;
t72 = -t114 * t682 + t181 * t685;
t985 = t287 * t72;
t288 = t489 * t682 - t1254;
t984 = t288 * t73;
t981 = t1170 * mrSges(6,1);
t262 = Ifges(6,1) * t1248 - Ifges(6,5) * t908 - t425;
t279 = mrSges(6,1) * t428 + mrSges(6,2) * t1248;
t402 = -Ifges(5,2) * t710 - Ifges(5,6) * t908 + t1021;
t555 = Ifges(5,4) * t710;
t403 = Ifges(5,1) * t559 - Ifges(5,5) * t908 - t555;
t708 = t710 * mrSges(5,3);
t523 = mrSges(5,2) * t908 - t708;
t622 = t664 + (-pkin(2) - t1039) * t681;
t575 = t622 - t1032;
t709 = t710 * mrSges(5,2);
t967 = t559 * mrSges(5,1);
t704 = -t709 + t967;
t705 = -Ifges(5,1) * t710 - t1021;
t717 = t1196 - t1226;
t783 = m(7) * t113 + t209 - t386;
t733 = -m(6) * t120 + t783;
t770 = -Ifges(6,1) * t428 - t1019;
t801 = -Ifges(5,2) * t559 - t555;
t812 = m(6) * t121 + t385;
t5 = t1181 * t840 + t812 * t127 - t1228 * t310 - t1224 * t428 + t402 * t1091 + t309 * t523 + t83 * t242 + t82 * t243 + t770 * t1105 + t717 * t1272 + t705 * t1090 + (t826 - t1205) * t1248 + (m(6) * t454 + t279) * t1038 + t454 * t773 + (t1216 + t1142) * t73 + (t1215 + t1143) * t72 + t575 * t704 + t1246 + (t768 + t894 + t262) * t1112 + t733 * t126 + (-t801 / 0.2e1 + t309 * mrSges(5,3) - t403 / 0.2e1) * t710;
t970 = t5 * qJD(1);
t969 = t516 * mrSges(6,1);
t968 = t517 * mrSges(6,2);
t964 = t592 * mrSges(6,3);
t962 = t738 * mrSges(6,3);
t227 = t1001 + t996 + t1009;
t325 = -mrSges(7,2) * t516 + mrSges(7,3) * t496;
t326 = mrSges(7,1) * t516 - mrSges(7,3) * t497;
t329 = Ifges(6,4) * t517 - Ifges(6,2) * t516 + Ifges(6,6) * t909;
t330 = Ifges(6,1) * t517 - Ifges(6,4) * t516 + Ifges(6,5) * t909;
t337 = t968 + t969;
t500 = -mrSges(6,2) * t909 - mrSges(6,3) * t516;
t505 = Ifges(5,4) * t612 + Ifges(5,2) * t611 + Ifges(5,6) * t909;
t506 = Ifges(5,1) * t612 + Ifges(5,4) * t611 + Ifges(5,5) * t909;
t959 = t612 * mrSges(5,2);
t960 = t611 * mrSges(5,1);
t520 = t959 - t960;
t588 = -mrSges(5,2) * t909 + mrSges(5,3) * t611;
t589 = mrSges(5,1) * t909 - mrSges(5,3) * t612;
t954 = t680 * Ifges(4,4);
t598 = (Ifges(4,6) * t684 + (-t678 * Ifges(4,2) + t954) * t686) * t679;
t599 = (Ifges(4,5) * t684 + (Ifges(4,1) * t680 - t1022) * t686) * t679;
t608 = mrSges(4,2) * t908 + t631 * mrSges(4,3);
t609 = -mrSges(4,1) * t908 - mrSges(4,3) * t741;
t618 = t1176 * t908;
t662 = Ifges(3,5) * t908;
t810 = t313 + t1043;
t501 = mrSges(6,1) * t909 - mrSges(6,3) * t517;
t811 = m(6) * t168 + t501;
t813 = m(4) * t573 + (-mrSges(4,2) * t684 - t1028 * t686) * t679;
t814 = m(4) * t572 + (mrSges(4,1) * t684 - mrSges(4,3) * t905) * t679;
t819 = m(7) * t73 + t242;
t956 = t635 * mrSges(3,2);
t6 = (Ifges(3,5) * t1057 - t1013 / 0.2e1 + (Ifges(4,4) * t910 + Ifges(4,2) * t631) * t1061 + (Ifges(4,1) * t910 + Ifges(4,4) * t631) * t1058 + t1165) * t908 + t505 * t1259 + t810 * t113 + t811 * t120 + t813 * t535 + t814 * t534 + (-Ifges(3,6) * t904 + (-t906 * t1022 / 0.2e1 + Ifges(4,1) * t1159 * t1052 - pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t680 + Ifges(4,6) * t678) * t686) * t908 + ((-t684 ^ 2 + t686 ^ 2) * Ifges(3,4) + ((Ifges(3,1) - Ifges(6,3) - Ifges(5,3) - Ifges(4,3) - Ifges(3,2)) * t686 - pkin(1) * mrSges(3,1)) * t684) * t679) * t679 + (Ifges(4,5) * t741 + Ifges(5,5) * t559 + Ifges(6,5) * t1248 + Ifges(4,6) * t631 - Ifges(5,6) * t710 - Ifges(6,6) * t428) * t841 + t631 * t598 / 0.2e1 + t622 * t618 + t573 * t608 + t572 * t609 + t356 * t523 + t355 * t524 + t522 * t279 + t169 * t385 + t168 * t386 + t73 * t325 + t99 * t243 + t158 * t209 + t229 * t1120 + t228 * t1122 + t819 * t100 + t330 * t1105 + t329 * t1112 + t506 * t1090 + t262 * t1092 + t162 * t1093 + t161 * t1094 + t403 * t1077 + t402 * t1078 + t662 * t1057 + (t326 + t1141) * t72 + (m(5) * t610 + t520) * t575 - t681 * t956 + t610 * (mrSges(5,1) * t710 + t559 * mrSges(5,2)) + (m(6) * t522 + t337) * t454 + (m(5) * t356 + t588) * t310 + (m(5) * t355 + t589) * t309 + (m(6) * t169 + t500) * t121 + t826 * t516 + t741 * t599 / 0.2e1 + (m(4) * t622 - mrSges(3,1) * t681 - t631 * mrSges(4,1) + mrSges(4,2) * t741) * t636 + t227 * t1272;
t961 = t6 * qJD(1);
t957 = t619 * mrSges(6,2);
t951 = t682 * t72;
t949 = t682 * t90;
t946 = t685 * t73;
t944 = t685 * t89;
t871 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t784 = -Ifges(6,1) / 0.2e1 + t871;
t817 = m(7) * t89 + t1215;
t818 = m(7) * t90 + t1216;
t7 = t90 * t242 + t89 * t243 + t120 * t385 + t887 * t840 + t818 * t73 + (t454 * mrSges(6,1) - t1019 / 0.2e1 + t826) * t1248 + (-t454 * mrSges(6,2) + t425 / 0.2e1 - t262 / 0.2e1 - t894 / 0.2e1 - t729 + t784 * t1248 - t1224) * t428 + t817 * t72 + (t783 - t971) * t121 + t1246;
t943 = t7 * qJD(1);
t20 = -(m(5) * t310 + t523) * t710 + (m(4) * t535 + t608) * t631 - (m(5) * t309 + t524) * t559 + t733 * t1248 - (t819 * t685 + (-m(7) * t72 - t243) * t682 + t812) * t428 + (-m(4) * t534 - t609) * t741;
t941 = qJD(1) * t20;
t939 = t113 * t1170;
t938 = t113 * t651;
t210 = Ifges(7,5) * t375 - Ifges(7,6) * t376;
t12 = -t376 * t161 / 0.2e1 + t210 * t1272 + t769 * t1120 + t113 * t208 + t72 * t242 - t73 * t243 + (-t375 * t72 - t376 * t73) * mrSges(7,3) + t1185 * t1122;
t937 = t12 * qJD(1);
t935 = t126 * t405;
t934 = t126 * t772;
t931 = t287 * t685;
t930 = t288 * t682;
t925 = t1170 * t772;
t924 = t405 * t1248;
t919 = t592 * t682;
t918 = t592 * t685;
t914 = t670 * t682;
t913 = t670 * t685;
t912 = t671 * t208;
t911 = t671 * t651;
t902 = t682 * t1217;
t899 = t682 * t1218;
t898 = t682 * t478;
t893 = t685 * t1213;
t891 = t685 * t1214;
t884 = qJD(5) * t682;
t883 = qJD(5) * t685;
t882 = 0.2e1 * m(5);
t880 = 0.2e1 * m(7);
t879 = 0.2e1 * pkin(4);
t878 = 0.2e1 * t670;
t877 = mrSges(6,3) * t1036;
t875 = -t1038 / 0.2e1;
t874 = t1036 / 0.2e1;
t873 = pkin(11) * t1056;
t872 = pkin(11) * t1050;
t866 = -t1008 / 0.2e1;
t864 = t1000 / 0.4e1;
t860 = -t971 / 0.2e1;
t859 = -t952 / 0.2e1;
t856 = t947 / 0.2e1;
t854 = t1071 - Ifges(6,6) / 0.2e1;
t853 = t682 * t1047;
t852 = t685 * t1047;
t846 = -t938 / 0.2e1;
t845 = t1248 * t1072;
t844 = t738 * t1072;
t843 = -t914 / 0.2e1;
t842 = t913 / 0.2e1;
t471 = t655 * t738;
t825 = t1125 - t471 / 0.4e1;
t470 = t738 * t653;
t824 = t1124 - t470 / 0.4e1;
t823 = t1068 - t656 / 0.4e1;
t822 = -t1178 / 0.4e1 - t655 / 0.4e1;
t809 = m(7) * t252 + t475;
t476 = -mrSges(7,3) * t918 + t1208;
t808 = m(7) * t287 + t476;
t473 = -mrSges(7,3) * t919 - t1206;
t807 = m(7) * t288 + t473;
t805 = -t880 / 0.4e1;
t804 = t880 / 0.4e1;
t803 = -t645 * mrSges(5,1) + t644 * mrSges(5,2);
t638 = Ifges(5,4) * t644;
t800 = Ifges(5,2) * t645 + t638;
t795 = t876 / 0.2e1;
t792 = t1050 * t1178 + t656 * t1054 + t1189;
t790 = t1047 * t1147;
t789 = m(6) * t879 / 0.2e1;
t787 = -t853 / 0.2e1;
t786 = -t1141 / 0.2e1 - t326 / 0.2e1;
t785 = Ifges(6,1) / 0.4e1 - Ifges(6,2) / 0.4e1 - Ifges(7,3) / 0.4e1;
t782 = m(7) * t1140 + t325 / 0.2e1;
t780 = t854 * t516;
t778 = t1135 - t983 / 0.2e1;
t777 = t468 + (m(6) + m(7)) * t405;
t776 = t982 / 0.2e1 + t243 / 0.2e1;
t775 = t592 * t821;
t774 = t964 + t1045;
t771 = Ifges(5,1) * t644 + t1020;
t760 = pkin(4) * t787;
t759 = t685 * t795;
t474 = -t592 * t952 - t1206;
t596 = t644 * Ifges(5,2) - t1020;
t597 = -Ifges(5,1) * t645 + t638;
t687 = t1087 * t1205 - t575 * t803 / 0.2e1 + (t1101 - t962 / 0.2e1) * t126 + (t800 / 0.4e1 + t597 / 0.4e1) * t710 + t486 * t875 + t1228 * t603 / 0.2e1 + t1171 * t1248 + (-t523 / 0.2e1 - t708 / 0.2e1) * t602 + (-t160 / 0.4e1 + t261 / 0.4e1 - Ifges(6,2) * t1109) * t738 - (t801 + t403) * t644 / 0.4e1 - t645 * t402 / 0.4e1 + t559 * t596 / 0.4e1 + t243 * t1132 + t270 * t1136 + t916 * t1138 + t474 * t1144 + t83 * t1100 + t113 * t1102 + t717 * t1085 + t82 * t1098 + t72 * t1099 + t770 * t1081 + t768 * t1083 + t1186 * t908 / 0.4e1 + t279 * t1037 / 0.2e1 - t127 * t964 / 0.2e1 - t559 * t771 / 0.4e1 + (t251 * t82 + t252 * t83 + t269 * t72 + t270 * t73 + t935) * t805 + t992 * t1243 - t1235 * t1118 + t405 * t1244 - t669 * t704 / 0.2e1 + t645 * t705 / 0.4e1 - m(6) * (t120 * t1158 - 0.2e1 * t936 + 0.2e1 * t935 + 0.2e1 * t127 * t1170 + (-t454 * t645 + t559 * t619) * t879) / 0.4e1 + t1250 * t1106 + (-t262 / 0.4e1 + Ifges(6,4) * t1109 + t1188 * t1111 + t1184) * t592 + t1195 * t1111 + (t113 * t805 + t1173) * t1170 - t1168 * t428 + t1261 - t1262;
t690 = t501 * t795 + t313 * t1064 + (t1047 * t168 + t169 * t683) * t788 + t516 * t1071 + t355 * mrSges(5,1) / 0.2e1 + t356 * t1149 + (t158 * t1157 + (t940 - t942) * t878) * t1155 + t326 * t843 + t500 * t874 + Ifges(5,3) * t841 + t325 * t842 - t1165 + t1166;
t1 = t690 + t687;
t19 = t669 * t803 + t596 * t1073 + t771 * t1074 + t619 * t1236 + t270 * t475 + t269 * t478 + (t474 + t1040) * t252 + (t477 + t1041) * t251 - t1179 * t1037 + t777 * t1170 + (t800 + t597) * t1075 + (t1025 / 0.2e1 + t1247) * t738 + (t1084 * t1188 - t871 * t738 + t1172) * t592 + (-t1045 + t1234) * t405;
t758 = -t1 * qJD(1) + t19 * qJD(2);
t23 = t287 * t478 + t288 * t475 + t1170 * t468 + t1247 * t738 - (t784 * t738 - t1172 + t1233 - t957) * t592 + (m(7) * t1170 + t1234) * t405 + t807 * t252 + t808 * t251;
t694 = -t1234 * t1211 + t121 * t1101 + t243 * t1129 + t288 * t1136 + t262 * t1083 + t160 * t1081 + t261 * t1080 - t72 * t476 / 0.2e1 + t473 * t1144 + t89 * t1098 + t90 * t1100 + t1261;
t734 = t686 * t764;
t3 = t1173 * t1170 + (t1235 / 0.2e1 + t1244) * t405 + t425 * t1085 + (-t939 / 0.2e1 - t936 / 0.2e1 - t988 / 0.2e1 - t986 / 0.2e1 - t985 / 0.2e1 - t984 / 0.2e1 + t158 * t1154) * m(7) + (t957 / 0.2e1 + t587 / 0.4e1 + t785 * t738 - t1222) * t428 + (Ifges(6,3) * t1052 + t734 / 0.4e1) * t679 + t694 + pkin(5) * t1128 + (t228 / 0.4e1 + t1217 * t1081 + t162 * t1083 + mrSges(7,3) * t1140 + (t866 + t1124) * t428 + t782 * pkin(11)) * t685 + t780 + t703 + (t229 / 0.4e1 + t1213 * t1080 + t161 * t1085 + t99 * t1146 + (t1000 / 0.2e1 + t1125) * t428 + t786 * pkin(11)) * t682 + (-t958 / 0.2e1 + t1204 / 0.2e1 - t785 * t592 + t1171) * t1248;
t757 = -t3 * qJD(1) + t23 * qJD(2);
t465 = t772 * t738;
t692 = (-t989 / 0.2e1 + t824) * t375 + (-t987 / 0.2e1 + t825) * t376 + (t769 * t1049 - t685 * t161 / 0.4e1 + (t951 / 0.2e1 - t946 / 0.2e1) * mrSges(7,3) - t1185 * t682 / 0.4e1) * t738 + t465 * t1211 + t251 * t1135 + t243 * t1134 + t208 * t1117 - t1190 * t1109 + t210 * t1083 + t72 * t475 / 0.2e1 + t73 * t1098;
t707 = t1009 / 0.2e1 + t1001 / 0.2e1 + t996 / 0.2e1 + t100 * t1147 + t99 * t1150;
t11 = t692 - t707;
t31 = t405 * t465 - t1190 * t1084 + t251 * t475 - t252 * t478 + ((-t987 - t345 / 0.2e1 - t471 / 0.2e1) * t685 + (t989 + t470 / 0.2e1 - t348 / 0.2e1) * t682) * t738;
t755 = t11 * qJD(1) + t31 * qJD(2);
t688 = t209 * t1087 + t386 * t1088 + t385 * t1243 + t523 * t1075 + t524 * t1073 + t609 * t1061 + t608 * t1058 + (t309 * t645 + t310 * t644 - t559 * t602 - t603 * t710) * t882 / 0.4e1 + (-t1170 * t428 - t120 * t738 + t924) * t806 + t1026 * t1091 + t1027 * t1259 + t428 * t1268 + (t121 * t806 + t1231) * t592 + t1232 * mrSges(4,3) / 0.2e1 + (t468 + t962) * t1105 + (-t964 + t898) * t1272 + m(4) * (qJ(3) * t1232 - t678 * t534 + t680 * t535) / 0.2e1 + (t113 * t738 + t924 + (t946 - t951) * t592 - (-t251 * t682 + t252 * t685) * t428) * t1260;
t693 = -t969 / 0.2e1 - t968 / 0.2e1 + t960 / 0.2e1 - t959 / 0.2e1 + t325 * t1056 + t326 * t1051 - m(4) * t636 / 0.2e1 - t610 * t882 / 0.4e1 - t522 * t881 / 0.4e1 + (t100 * t682 + t685 * t99) * t805 + t1176 * t840;
t14 = t693 + t688;
t43 = (m(5) * t602 + t1026) * t645 + (m(5) * t603 + t1027) * t644 + (t777 + t962) * t738 + (-t1180 * t682 + t809 * t685 + t774) * t592 + (t678 ^ 2 + t1159) * (m(4) * qJ(3) + mrSges(4,3));
t754 = -qJD(1) * t14 - qJD(2) * t43;
t695 = mrSges(5,1) * t1091 + mrSges(6,2) * t1272 - t1149 * t710 - t1163 * t428 + t1164 * t1248 + t1269 + t845;
t696 = -t967 / 0.2e1 + t709 / 0.2e1 + m(6) * t875 + (t682 * t83 + t685 * t82) * t805 - t1265;
t28 = t695 + t696;
t697 = t1163 * t592 + t1164 * t738 + t844;
t701 = t474 * t1054 + t477 * t1050 - m(6) * t1037 / 0.2e1 + (t269 * t685 + t270 * t682) * t804;
t52 = t697 - t701 - t802 - t803;
t753 = qJD(1) * t28 + qJD(2) * t52;
t32 = t845 + 0.2e1 * t1269 - (t1148 + t1193) * t428 + (pkin(5) * t1114 - t949 / 0.2e1 - t944 / 0.2e1 - t821 * t1271) * m(7) + t719;
t718 = t476 * t1051 + t473 * t1056 + t1089;
t55 = t844 + t1089 + mrSges(7,3) * t775 + 0.2e1 * t1088 * mrSges(6,1) + (-t1034 / 0.2e1 - t931 / 0.2e1 - t930 / 0.2e1 + pkin(11) * t775) * m(7) + t718;
t752 = qJD(1) * t32 + qJD(2) * t55;
t38 = (mrSges(7,2) * t1272 - t778) * t685 + (mrSges(7,1) * t1272 + t776) * t682;
t744 = -t948 / 0.2e1 - t953 / 0.2e1;
t728 = t744 * t592;
t739 = t898 / 0.2e1 + t1268;
t76 = t728 + t739;
t751 = qJD(1) * t38 + qJD(2) * t76;
t750 = pkin(11) * t1193;
t749 = qJD(5) * (-t1152 + t1192);
t748 = -t1152 / 0.2e1 + t1175;
t745 = mrSges(7,1) * t1129 + t288 * mrSges(7,2) / 0.2e1;
t743 = pkin(11) * t1100 + t825;
t742 = pkin(11) * t1098 + t824;
t740 = t1083 * t672 + t1117 * t651;
t737 = t1065 * t475 + t825;
t736 = t1065 * t478 + t824;
t732 = 0.2e1 * t1229;
t731 = 0.2e1 * t1230;
t730 = t1187 * t1047;
t691 = (t1225 + (-t287 * t682 + t288 * t685) * t878 + (-t251 * t853 + t252 * t852 + t1255) * t879) * t1155 + mrSges(6,1) * t1116 + t1170 * t1072 + t1234 * t1064 + t1218 * t1053 + t1214 * t1049 + t287 * t859 + t288 * t856 + t476 * t843 + t468 * t874 + t473 * t842 + t478 * t760 + t475 * t759 + t1263 + t1264;
t722 = pkin(5) * t1102 + t1264;
t22 = -t691 + t891 / 0.4e1 + t899 / 0.4e1 - t925 / 0.2e1 - t981 / 0.2e1 + (pkin(5) * t1158 + pkin(11) * t731) * t1155 + t477 * t873 + t269 * t859 + t474 * t872 + t270 * t856 + t722 + t1263;
t700 = t1192 * t1036 + (mrSges(7,3) * t1187 - mrSges(6,2)) * t876;
t464 = (t670 * t730 + t671 * t683) * t1153 + t700;
t689 = ((-t682 * t89 + t685 * t90) * t878 + (t113 * t683 - t72 * t853 + t73 * t852) * t879) * t1155 - t993 / 0.2e1 - t1235 * t1064 + t1217 * t1053 + t1213 * t1049 + t1215 * t843 + t209 * t874 + t1216 * t842 + t89 * t859 + t90 * t856 + t243 * t760 + t242 * t759 + (t385 + t975) * t795 + t1267 + (t1175 + t1177) * t121 + (-t386 / 0.2e1 + t860) * t1036 + t1249;
t712 = -t1235 * t1154 - t990 / 0.2e1 + t1249;
t8 = -t689 + (pkin(11) * t732 + t1161 * t126) * t1155 + t893 / 0.4e1 + t902 / 0.4e1 - t934 / 0.2e1 - t991 / 0.2e1 + t1215 * t873 + t82 * t859 + t1216 * t872 + t83 * t856 + t712 + t1267;
t726 = -t8 * qJD(1) - t22 * qJD(2) + t464 * qJD(4);
t706 = t83 * t1147 + t82 * t1150 + t1240 - t729;
t15 = t846 - t912 / 0.2e1 + ((t1050 * t376 + t1056 * t375) * mrSges(7,3) + t1221) * t670 + t706 - t1251;
t716 = t682 * t822 - t685 * t823;
t699 = (-t1193 * t670 + t716) * t738 + t465 * t1064 + t740;
t721 = -t1195 / 0.2e1 + mrSges(7,1) * t1132 + mrSges(7,2) * t1131;
t34 = ((t1085 + t1243) * Ifges(7,6) + t737) * t682 + (t1145 * t592 + t736) * t685 + t699 + t721;
t495 = t792 + t911;
t724 = -t15 * qJD(1) + t34 * qJD(2) + t495 * qJD(4);
t723 = t90 * t1147 + t89 * t1150 + t1240;
t720 = t1154 * t465 + t740;
t17 = t672 * t1111 + t846 + t1035 / 0.2e1 + t823 * t376 + t822 * t375 + (-t1010 / 0.2e1 - t162 / 0.4e1 + t1007 / 0.4e1 - t370 / 0.4e1 + t776 * pkin(11)) * t685 + (0.3e1 / 0.4e1 * t974 - t1024 / 0.4e1 + t1016 / 0.4e1 + t161 / 0.4e1 + t778 * pkin(11)) * t682 + t723;
t36 = (-Ifges(7,3) / 0.2e1 - t750) * t738 + (-t738 * t823 + t742 + t866) * t685 + (0.3e1 / 0.4e1 * t1000 + t822 * t738 + t743) * t682 + t720 + t745;
t379 = (-t671 / 0.2e1 + pkin(5) / 0.2e1) * t651 + (pkin(4) * t790 + t1067 - t1178 / 0.2e1) * t685 + (-mrSges(7,1) * t876 / 0.2e1 - t656 / 0.2e1 + t1069) * t682;
t504 = t792 - t1033;
t715 = t17 * qJD(1) - t36 * qJD(2) + t379 * qJD(4) - t504 * qJD(5);
t714 = t938 / 0.2e1 + (t858 + t859) * t73 + t1251;
t713 = pkin(11) * t1227 + t1193;
t380 = t911 / 0.2e1 - t1033 / 0.2e1 + (mrSges(7,1) * t787 + t685 * t790) * pkin(4) + t792;
t77 = t728 - t739;
t57 = t697 + t701;
t56 = (t930 + t931) * t804 + t963 / 0.2e1 + t1089 + t748 * t738 + t713 * t592 - t718;
t39 = t375 * t857 + t376 * t859 - t428 * t744 + t1231;
t37 = -Ifges(7,6) * t919 / 0.2e1 - t918 * t1145 + t1195 / 0.2e1 + (t864 + t743) * t682 + t742 * t685 + (-t750 + t716) * t738 + t720 - t745;
t35 = (t864 + t737) * t682 + t736 * t685 + t699 - t721 + t1233;
t33 = (t944 + t949) * t804 + t748 * t1248 - (t1148 + t713) * t428 + t1265;
t29 = t695 - t696;
t21 = t691 + (t1214 / 0.4e1 + t592 * t1066 + mrSges(7,3) * t1131 + (t1040 / 0.2e1 + t474 / 0.2e1) * pkin(11)) * t685 + t854 * t738 + t748 * t1170 + (t1218 / 0.4e1 + t592 * t1070 + mrSges(7,3) * t1132 + (-t1041 / 0.2e1 + t1099) * pkin(11)) * t682 + t722;
t18 = -t729 - t1035 / 0.2e1 + t714 + t723 + t1169 * pkin(11);
t16 = t912 / 0.2e1 + t706 + t714 + t1169 * t670;
t13 = -t693 + t688;
t10 = t692 + t707;
t9 = t748 * t126 + t689 + (t1270 - t428 * t1070 + t82 * t1146 + (-t1143 / 0.2e1 + t1130) * pkin(11)) * t682 + (t1138 - t428 * t1066 + t83 * mrSges(7,3) / 0.2e1 + (t1142 / 0.2e1 + t1216 / 0.2e1) * pkin(11)) * t685 + t712;
t4 = t780 + t1166 + (t342 + t1250) * t1248 / 0.4e1 - t694 + t209 * t1115 + t386 * t1116 + t385 * t1118 + t487 * t1106 + (t682 * t786 + t685 * t782) * pkin(11) + (-t1043 / 0.2e1 + t1128) * pkin(5) + (t936 + t939 + t984 + t985 + t986 + t988) * t804 - t1235 * t1117 - t1213 * t916 / 0.4e1 + t768 * t1085 - t1019 * t1080 - t679 * t734 / 0.4e1 + (-Ifges(6,2) * t738 + t587) * t1111 - (-t1109 * t1188 + t1184) * t592 + (-Ifges(6,1) * t1080 - t1083 * t1188 + t1168) * t428 + t1195 * t1109 + t1196 * t1083 + t1170 * t860 + t1262;
t2 = t690 - t687;
t24 = [qJD(2) * t6 + qJD(3) * t20 + qJD(4) * t5 + qJD(5) * t7 + qJD(6) * t12, t13 * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t10 * qJD(6) + t961 + ((-t168 * mrSges(6,3) + t330 / 0.2e1 + t228 * t1056 + t229 * t1050) * t738 + t669 * t520 + t610 * (-mrSges(5,1) * t644 - mrSges(5,2) * t645) + t619 * t337 + (-m(4) * pkin(2) - mrSges(4,1) * t680 + mrSges(4,2) * t678 - mrSges(3,1)) * t636 + t774 * t169 + t809 * t100 + (t810 - t811) * t405 + (-t678 * t814 + t680 * t813) * qJ(3) - pkin(2) * t618 + t602 * t589 + t603 * t588 + t158 * t468 + t252 * t325 + t251 * t326 + t227 * t1084 + t488 * t1092 + t348 * t1093 + t345 * t1094 + t506 * t1074 + t505 * t1075 + t597 * t1077 + t596 * t1078 + t598 * t1058 + t599 * t1060 + t1179 * t522 + t1180 * t99 + t1174 * t516 + t355 * t1026 + t356 * t1027 - t572 * t1028 + (t355 * t602 + t356 * t603 + t610 * t669) * t882 / 0.2e1 + ((Ifges(4,5) * t1060 + Ifges(5,5) * t1074 + Ifges(6,5) * t1087 + Ifges(4,6) * t1058 + Ifges(5,6) * t1075 + Ifges(6,6) * t1243 - Ifges(3,6)) * t684 + ((Ifges(4,2) * t680 + t1022) * t1061 + (Ifges(4,1) * t678 + t954) * t1058) * t686) * t679 + t329 * t1243 - t956 + t1170 * t500 + t573 * t680 * mrSges(4,3) + t662) * qJD(2), qJD(2) * t13 + qJD(4) * t29 + qJD(5) * t33 + qJD(6) * t39 + t941, t970 + t2 * qJD(2) + t29 * qJD(3) + ((-t1047 * t126 + t127 * t683) * t789 + t893 / 0.2e1 + t902 / 0.2e1 - t671 * t1235 - t934 + t1241 - t310 * mrSges(5,1) - t309 * mrSges(5,2) - t990 - t991 - t1248 * t877 + (t1157 * t126 + t670 * t732) * t1156 - t1215 * t914 + t1216 * t913 + t1181 + t1229 * mrSges(7,3)) * qJD(4) + t9 * qJD(5) + t16 * qJD(6) - t428 * t1167, t943 + t4 * qJD(2) + t33 * qJD(3) + t9 * qJD(4) + (pkin(5) * t1235 + t1241 + t887 - t993) * qJD(5) + t18 * qJD(6) + t121 * t749 + (t428 * t1067 + t1213 / 0.2e1 + t90 * mrSges(7,3) + t818 * pkin(11)) * t883 + (t428 * t1069 + t1217 / 0.2e1 - t89 * mrSges(7,3) - t817 * pkin(11)) * t884, t937 + t10 * qJD(2) + t39 * qJD(3) + t16 * qJD(4) + t18 * qJD(5) + (-mrSges(7,1) * t73 - mrSges(7,2) * t72 + t210) * qJD(6); qJD(3) * t14 - qJD(4) * t1 - qJD(5) * t3 + qJD(6) * t11 - t961, qJD(3) * t43 + qJD(4) * t19 + qJD(5) * t23 + qJD(6) * t31, qJD(4) * t57 + qJD(5) * t56 + qJD(6) * t77 - t754, t57 * qJD(3) + (t891 / 0.2e1 + t899 / 0.2e1 + t671 * t1234 + (-t1047 * t1170 - t1255) * t789 - t925 - t602 * mrSges(5,2) - t603 * mrSges(5,1) - t981 + (t670 * t731 + t1225) * t1156 - t738 * t877 - t477 * t914 + t474 * t913 + t1186 + t1230 * mrSges(7,3) + t1266) * qJD(4) + t21 * qJD(5) + t35 * qJD(6) + t758 + t592 * t1167, t56 * qJD(3) + t21 * qJD(4) + (-pkin(5) * t1234 + t1266 + t764) * qJD(5) + t37 * qJD(6) + t1170 * t749 + (-t592 * t1067 + t1214 / 0.2e1 + t288 * mrSges(7,3) + t807 * pkin(11)) * t883 + (-t592 * t1069 + t1218 / 0.2e1 - t287 * mrSges(7,3) - t808 * pkin(11)) * t884 + t757, t77 * qJD(3) + t35 * qJD(4) + t37 * qJD(5) + (-mrSges(7,1) * t252 - mrSges(7,2) * t251 - t1190) * qJD(6) + t755; -qJD(2) * t14 - qJD(4) * t28 - qJD(5) * t32 - qJD(6) * t38 - t941, -qJD(4) * t52 - qJD(5) * t55 - qJD(6) * t76 + t754, 0, -t753, -t752, -qJD(6) * t651 - t751; qJD(2) * t1 + qJD(3) * t28 - qJD(5) * t8 - qJD(6) * t15 - t970, qJD(3) * t52 - qJD(5) * t22 + qJD(6) * t34 - t758, t753, qJD(5) * t464 + qJD(6) * t495 ((0.2e1 * pkin(11) * t730 + t1161 * t683) * t1153 / 0.2e1 + t700) * qJD(5) + t380 * qJD(6) + t726, t380 * qJD(5) + (-t670 * t772 + t1188) * qJD(6) + t724; qJD(2) * t3 + qJD(3) * t32 + qJD(4) * t8 - qJD(6) * t17 - t943, qJD(3) * t55 + qJD(4) * t22 + qJD(6) * t36 - t757, t752, -qJD(6) * t379 - t726, t504 * qJD(6) (-pkin(11) * t772 + t1188) * qJD(6) - t715; -qJD(2) * t11 + qJD(3) * t38 + qJD(4) * t15 + qJD(5) * t17 - t937, qJD(3) * t76 - qJD(4) * t34 - qJD(5) * t36 - t755, t751, qJD(5) * t379 - t724, t715, 0;];
Cq  = t24;
