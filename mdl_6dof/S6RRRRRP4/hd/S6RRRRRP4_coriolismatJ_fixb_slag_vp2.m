% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRRP4
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:15
% EndTime: 2019-03-10 01:13:03
% DurationCPUTime: 28.80s
% Computational Cost: add. (54136->970), mult. (107522->1214), div. (0->0), fcn. (118714->8), ass. (0->549)
t1212 = Ifges(6,1) + Ifges(7,1);
t1165 = Ifges(7,4) + Ifges(6,5);
t1200 = Ifges(6,6) - Ifges(7,6);
t607 = sin(qJ(4));
t1035 = pkin(4) * t607;
t1046 = cos(qJ(5));
t858 = t1046 * t607;
t606 = sin(qJ(5));
t608 = cos(qJ(4));
t911 = t606 * t608;
t759 = t858 + t911;
t1056 = t759 / 0.2e1;
t912 = t606 * t607;
t1124 = t1046 * t608 - t912;
t1057 = t1124 / 0.2e1;
t1044 = sin(qJ(3));
t1047 = cos(qJ(3));
t1045 = sin(qJ(2));
t890 = t1045 * pkin(7);
t583 = -pkin(8) * t1045 - t890;
t1048 = cos(qJ(2));
t894 = t1048 * pkin(7);
t584 = pkin(8) * t1048 + t894;
t1125 = t1044 * t583 + t1047 * t584;
t569 = t1044 * t1045 - t1047 * t1048;
t570 = -t1044 * t1048 - t1045 * t1047;
t598 = -pkin(2) * t1048 - pkin(1);
t458 = t569 * pkin(3) + t570 * pkin(9) + t598;
t282 = -t1125 * t607 + t608 * t458;
t920 = t570 * t608;
t221 = pkin(10) * t920 + t282;
t926 = t1125 * t608;
t222 = t926 + (pkin(10) * t570 + t458) * t607;
t860 = t1046 * t222;
t122 = t221 * t606 + t860;
t913 = t606 * t222;
t123 = t1046 * t221 - t913;
t1152 = t1056 * t122 + t1057 * t123;
t1171 = -t920 / 0.2e1;
t477 = -mrSges(7,1) * t1124 - mrSges(7,3) * t759;
t1176 = -t477 / 0.2e1;
t1132 = t570 * t759;
t434 = t1124 * t570;
t239 = -mrSges(7,1) * t1132 + t434 * mrSges(7,3);
t1178 = -t239 / 0.2e1;
t1026 = Ifges(5,4) * t608;
t1210 = -Ifges(5,1) * t607 / 0.2e1 - t1026 / 0.2e1;
t240 = -mrSges(6,1) * t1132 - t434 * mrSges(6,2);
t473 = pkin(5) * t759 - qJ(6) * t1124;
t459 = t473 + t1035;
t478 = -mrSges(6,1) * t1124 + mrSges(6,2) * t759;
t505 = t1044 * t584 - t1047 * t583;
t979 = t607 * mrSges(5,1);
t580 = mrSges(5,2) * t608 + t979;
t669 = -t434 * pkin(5) - qJ(6) * t1132;
t896 = pkin(4) * t920;
t643 = t669 - t896;
t1027 = Ifges(5,4) * t607;
t788 = Ifges(5,1) * t608 - t1027;
t697 = Ifges(5,5) * t569 - t570 * t788;
t680 = t608 * t697;
t782 = -Ifges(5,2) * t607 + t1026;
t696 = Ifges(5,6) * t569 - t570 * t782;
t681 = t607 * t696;
t778 = Ifges(5,5) * t608 - Ifges(5,6) * t607;
t746 = t569 * t778;
t781 = Ifges(5,2) * t608 + t1027;
t921 = t570 * t607;
t1211 = -t1152 * mrSges(7,2) + t643 * t1176 - t680 / 0.4e1 + t681 / 0.4e1 + t781 * t1171 + t788 * t920 / 0.4e1 - t746 / 0.4e1 + t478 * t896 / 0.2e1 - t240 * t1035 / 0.2e1 - t505 * t580 / 0.2e1 + t459 * t1178 + (-t782 / 0.4e1 + t1210) * t921;
t1050 = t607 / 0.2e1;
t1136 = mrSges(7,1) + mrSges(6,1);
t1135 = mrSges(6,3) + mrSges(7,2);
t1081 = -pkin(10) - pkin(9);
t602 = t608 * pkin(10);
t603 = t608 * pkin(9);
t902 = t603 + t602;
t503 = -t1081 * t858 + t606 * t902;
t1189 = t503 * t1124;
t416 = mrSges(7,2) * t1189;
t1160 = Ifges(7,5) * t1132;
t706 = Ifges(6,4) * t1132;
t1127 = t1165 * t569 - t1212 * t434 - t1160 + t706;
t1209 = Ifges(6,2) * t434 + t1127 + t706;
t792 = t1047 * t1046;
t522 = (t1047 * t911 + t607 * t792) * pkin(2);
t861 = t607 * t1047;
t523 = (-t606 * t861 + t608 * t792) * pkin(2);
t1208 = t1135 * (t1124 * t523 + t522 * t759);
t1062 = t503 / 0.2e1;
t1063 = -t503 / 0.2e1;
t1207 = t1062 + t1063;
t923 = t569 * t607;
t1142 = -pkin(4) * t923 + t1125;
t435 = t759 * t569;
t436 = t569 * t1124;
t1129 = -t435 * pkin(5) + t436 * qJ(6) + t1142;
t1145 = t1129 * t477;
t579 = -mrSges(5,1) * t608 + t607 * mrSges(5,2);
t1149 = t1125 * t579;
t1159 = t1125 * mrSges(4,1);
t1163 = t505 * mrSges(4,2);
t1192 = t1142 * t478;
t1049 = t608 / 0.2e1;
t1052 = -t570 / 0.2e1;
t1071 = -t436 / 0.2e1;
t1072 = -t435 / 0.2e1;
t1073 = t435 / 0.2e1;
t1153 = t1050 * t781 + t608 * t1210;
t1156 = -t1165 * t570 - t1212 * t436 + (-Ifges(7,5) + Ifges(6,4)) * t435;
t1175 = -t1124 / 0.2e1;
t1018 = Ifges(7,5) * t1124;
t1023 = Ifges(6,4) * t1124;
t1195 = t1212 * t759 - t1018 + t1023;
t191 = -Ifges(7,5) * t436 - Ifges(7,6) * t570 - Ifges(7,3) * t435;
t192 = -Ifges(6,4) * t436 + Ifges(6,2) * t435 - Ifges(6,6) * t570;
t347 = -Ifges(5,6) * t570 - t569 * t782;
t348 = -Ifges(5,5) * t570 - t569 * t788;
t745 = t570 * (Ifges(5,5) * t607 + Ifges(5,6) * t608);
t1017 = Ifges(7,5) * t759;
t777 = -Ifges(7,3) * t1124 + t1017;
t1022 = Ifges(6,4) * t759;
t780 = Ifges(6,2) * t1124 + t1022;
t636 = t191 * t1175 + t192 * t1057 + t777 * t1072 + t780 * t1073 + t348 * t1050 + t347 * t1049 - t745 / 0.2e1 + Ifges(4,6) * t570 + t1195 * t1071 + t1156 * t1056 + (t1124 * t1200 + t1165 * t759) * t1052 + (-Ifges(4,5) + t1153) * t569;
t1206 = t1145 + t1149 + t1163 + t636 - t1159 + t1192;
t1205 = t1145 / 0.2e1 + t1163 / 0.2e1 + t1192 / 0.2e1;
t984 = t1124 * mrSges(7,2);
t866 = -t984 / 0.2e1;
t867 = t984 / 0.2e1;
t1123 = t867 + t866;
t1204 = (qJD(2) + qJD(3)) * t1123;
t1083 = -mrSges(6,2) / 0.2e1;
t1082 = mrSges(7,3) / 0.2e1;
t670 = -t434 * mrSges(7,1) - mrSges(7,3) * t1132;
t1202 = -t670 / 0.2e1;
t671 = -t434 * mrSges(6,1) + mrSges(6,2) * t1132;
t1201 = -t671 / 0.2e1;
t794 = t1124 * t1165 - t1200 * t759;
t605 = t608 ^ 2;
t1128 = t607 ^ 2 + t605;
t893 = t1047 * pkin(2);
t1194 = t1128 * t893;
t389 = -pkin(4) * t921 + t505;
t1193 = t1142 * t389;
t1033 = t608 * pkin(4);
t596 = -t893 - pkin(3);
t577 = t596 - t1033;
t1191 = t1142 * t577;
t597 = -pkin(3) - t1033;
t1190 = t1142 * t597;
t1091 = m(7) / 0.2e1;
t1111 = (t1082 + t1083) * t523 - t1136 * t522 / 0.2e1;
t1036 = pkin(4) * t606;
t588 = qJ(6) + t1036;
t892 = t1046 * pkin(4);
t595 = -t892 - pkin(5);
t806 = -t893 / 0.2e1;
t774 = t608 * t806;
t1089 = m(6) * pkin(4);
t899 = t1089 / 0.2e1;
t1188 = t1111 + (t522 * t595 + t523 * t588) * t1091 + (-t1046 * t522 + t523 * t606) * t899 + t806 * t979 + mrSges(5,2) * t774;
t461 = t1123 * qJD(6);
t1019 = Ifges(7,5) * t434;
t1024 = Ifges(6,4) * t434;
t648 = Ifges(7,6) * t569 - Ifges(7,3) * t1132 - t1019;
t1187 = t1132 * t1212 - t1019 + t1024 + t648;
t1186 = t1124 * t1212 + t1017 - t1022 + t777;
t1185 = -Ifges(6,2) * t759 + t1023 + t1195;
t1183 = t1149 / 0.2e1 - t1159 / 0.2e1;
t1182 = (qJD(4) - qJD(5)) * t1123;
t1012 = Ifges(6,3) * t570;
t1013 = Ifges(5,3) * t570;
t1014 = Ifges(7,6) * t435;
t1015 = Ifges(6,6) * t435;
t1016 = Ifges(7,2) * t570;
t1020 = Ifges(6,5) * t436;
t1021 = Ifges(7,4) * t436;
t1028 = Ifges(4,4) * t570;
t1053 = t569 / 0.2e1;
t1054 = -t569 / 0.2e1;
t1075 = -t434 / 0.2e1;
t180 = t569 * pkin(4) + t221;
t109 = t1046 * t180 - t913;
t110 = t606 * t180 + t860;
t1167 = -t1132 / 0.2e1;
t951 = t434 * qJ(6);
t150 = -pkin(5) * t1132 + t389 + t951;
t237 = -mrSges(7,1) * t435 + mrSges(7,3) * t436;
t238 = -mrSges(6,1) * t435 - mrSges(6,2) * t436;
t283 = t607 * t458 + t926;
t323 = mrSges(6,2) * t570 + mrSges(6,3) * t435;
t325 = -mrSges(6,1) * t570 + mrSges(6,3) * t436;
t1001 = t436 * mrSges(7,2);
t980 = t570 * mrSges(7,1);
t326 = t980 - t1001;
t330 = mrSges(7,2) * t435 - mrSges(7,3) * t570;
t455 = t580 * t569;
t456 = t580 * t570;
t465 = mrSges(5,2) * t570 + mrSges(5,3) * t923;
t922 = t569 * t608;
t467 = -t570 * mrSges(5,1) + mrSges(5,3) * t922;
t649 = Ifges(6,2) * t1132 + Ifges(6,6) * t569 - t1024;
t700 = t1132 / 0.2e1;
t829 = t921 / 0.2e1;
t552 = t569 * qJ(6);
t85 = t552 + t110;
t86 = -t569 * pkin(5) - t109;
t1181 = t348 * t1171 + t347 * t829 + t192 * t700 - t505 * t455 + t283 * t465 + t282 * t467 + t389 * t238 + t85 * t330 + t110 * t323 + t109 * t325 + t86 * t326 + t150 * t237 + t1156 * t1075 - t1125 * t456 + t1129 * t239 + t648 * t1072 + t649 * t1073 + (-t1012 - t1013 - t1014 + t1015 - t1016 - t1020 - t1021 + t681 - t746) * t1053 + t191 * t1167 + t1142 * t240 + (-0.2e1 * Ifges(4,4) * t569 + t680 + (-Ifges(4,1) + Ifges(4,2)) * t570) * t1054 + t1127 * t1071 + (-t570 * t778 + t1028 - t1165 * t434 + (-Ifges(4,1) + Ifges(7,2) + Ifges(5,3) + Ifges(6,3)) * t569 + t1200 * t1132) * t1052 + t598 * (-mrSges(4,1) * t570 - mrSges(4,2) * t569) + t570 * (-Ifges(4,2) * t569 - t1028) / 0.2e1;
t1055 = -t759 / 0.2e1;
t1113 = t1046 * t902 + t1081 * t912;
t1170 = -t1113 / 0.2e1;
t1169 = t1113 / 0.2e1;
t704 = mrSges(6,3) * t1132;
t324 = -t569 * mrSges(6,2) + t704;
t1161 = mrSges(7,2) * t1132;
t981 = t569 * mrSges(7,3);
t329 = t1161 + t981;
t1131 = t324 + t329;
t1168 = -t1131 / 0.2e1;
t1134 = mrSges(7,3) - mrSges(6,2);
t1164 = pkin(3) * t1125;
t1162 = t110 - t85;
t1158 = t607 * t505;
t1157 = t608 * t505;
t318 = t473 * t477;
t776 = Ifges(7,3) * t759 + t1018;
t714 = t780 * t1055 + t1056 * t1186 + t1057 * t1185 + t776 * t1175;
t1155 = t318 + t714;
t1154 = t595 + t892;
t600 = m(7) * qJ(6) + mrSges(7,3);
t1150 = qJD(5) * t600;
t1148 = t1129 * t150;
t775 = -pkin(5) * t1124 - qJ(6) * t759;
t457 = t597 + t775;
t425 = -t893 + t457;
t1147 = t1129 * t425;
t1146 = t1129 * t457;
t933 = t505 * t1125;
t1144 = t596 * t1125;
t1143 = t600 * qJD(6);
t1141 = t1201 * t597 + t1202 * t457;
t1140 = t1201 * t577 + t1202 * t425;
t1084 = -mrSges(5,2) / 0.2e1;
t1087 = mrSges(5,1) / 0.2e1;
t479 = -t570 * pkin(3) + t569 * pkin(9);
t304 = t608 * t479 + t1158;
t791 = -t570 * pkin(4) + pkin(10) * t922;
t195 = t304 + t791;
t305 = t607 * t479 - t1157;
t898 = pkin(10) * t923;
t242 = t898 + t305;
t116 = t1046 * t195 - t606 * t242;
t117 = t1046 * t242 + t606 * t195;
t553 = t570 * qJ(6);
t97 = -t553 + t117;
t1034 = t570 * pkin(5);
t98 = -t116 + t1034;
t1139 = (t1046 * t116 + t117 * t606) * t899 + (t588 * t97 + t595 * t98) * t1091 + t305 * t1084 + t304 * t1087;
t891 = t1045 * pkin(2);
t460 = t891 + t479;
t292 = t608 * t460 + t1158;
t181 = t292 + t791;
t293 = t607 * t460 - t1157;
t223 = t898 + t293;
t111 = t1046 * t181 - t606 * t223;
t112 = t1046 * t223 + t606 * t181;
t89 = -t553 + t112;
t90 = -t111 + t1034;
t1138 = (t1046 * t111 + t112 * t606) * t899 + (t588 * t89 + t595 * t90) * t1091 + t293 * t1084 + t292 * t1087;
t880 = t1036 / 0.2e1;
t1133 = t109 + t86;
t982 = t759 * mrSges(6,3);
t417 = t1113 * t982;
t1130 = -t416 + t417;
t795 = t1132 * t1165 + t1200 * t434;
t475 = mrSges(6,1) * t759 + mrSges(6,2) * t1124;
t918 = t577 * t475;
t474 = mrSges(7,1) * t759 - mrSges(7,3) * t1124;
t952 = t425 * t474;
t1122 = t918 + t952;
t960 = t293 * t608;
t961 = t292 * t607;
t1121 = t960 - t961;
t1120 = Ifges(5,4) * t605 + (-t1027 + pkin(4) * t478 + (Ifges(5,1) - Ifges(5,2)) * t608) * t607;
t1118 = t1176 * t669 + t1178 * t473;
t1117 = m(7) * t150 + t239;
t1116 = mrSges(5,3) * t1052 * t1128;
t1043 = m(7) * t110;
t666 = t981 + t1161 / 0.2e1;
t632 = (0.2e1 * t552 + t110) * t1091 + t666;
t694 = mrSges(7,2) * t1167;
t50 = t694 - t1043 / 0.2e1 + t632;
t1115 = qJD(1) * t50 + qJD(4) * t600 + t1150 + t1204;
t889 = t1044 * pkin(2);
t793 = t889 + pkin(9);
t744 = t607 * (-pkin(10) - t793);
t586 = t608 * t793;
t903 = t586 + t602;
t452 = -t1046 * t744 + t606 * t903;
t1114 = t1046 * t903 + t606 * t744;
t1112 = -m(7) * pkin(5) - t1136;
t1110 = mrSges(6,2) - t600;
t1109 = m(7) * t85 + t1131;
t1002 = t434 * mrSges(6,3);
t327 = mrSges(6,1) * t569 + t1002;
t1003 = t434 * mrSges(7,2);
t328 = -mrSges(7,1) * t569 - t1003;
t1108 = m(7) * t86 - t327 + t328;
t1107 = -t1047 * mrSges(4,2) + (-mrSges(4,1) + t477 + t478 + t579) * t1044;
t938 = t1113 * t759;
t1106 = -t1189 / 0.2e1 + t938 / 0.2e1 + (t1055 + t1056) * t1114;
t1085 = -mrSges(7,1) / 0.2e1;
t1086 = mrSges(6,1) / 0.2e1;
t668 = t1082 * t89 + t1083 * t112 + t1085 * t90 + t1086 * t111;
t950 = t434 * t1114;
t1105 = -t1140 + t668 + t1135 * (t452 * t700 + t950 / 0.2e1);
t667 = t1082 * t97 + t1083 * t117 + t1085 * t98 + t1086 * t116;
t928 = t1113 * t434;
t1104 = -t1141 + t667 + t1135 * (t503 * t700 + t928 / 0.2e1);
t1103 = t1141 + t667 + t1135 * (t503 * t1167 - t928 / 0.2e1);
t1102 = t1140 + t668 + t1135 * (t452 * t1167 - t950 / 0.2e1);
t1051 = t588 / 0.2e1;
t1078 = t326 / 0.2e1;
t805 = t892 / 0.2e1;
t1101 = t325 * t805 + Ifges(5,6) * t923 / 0.2e1 - Ifges(5,5) * t922 / 0.2e1 + t323 * t880 + t595 * t1078 + t330 * t1051 - t1013 / 0.2e1;
t1038 = pkin(3) * t455;
t1093 = m(6) / 0.2e1;
t1095 = m(5) / 0.2e1;
t909 = t608 * t465;
t910 = t607 * t467;
t914 = t597 * t238;
t945 = t457 * t237;
t969 = t112 * t1124;
t970 = t111 * t759;
t975 = t90 * t759;
t976 = t89 * t1124;
t1100 = (pkin(9) * t1121 - t1164) * t1095 + (-t111 * t503 + t1113 * t112 + t1190) * t1093 + (t1113 * t89 + t503 * t90 + t1146) * t1091 + t1038 / 0.2e1 + t945 / 0.2e1 + t914 / 0.2e1 + (t909 / 0.2e1 - t910 / 0.2e1) * pkin(9) + (t960 / 0.2e1 - t961 / 0.2e1) * mrSges(5,3) + (-t970 / 0.2e1 + t969 / 0.2e1) * mrSges(6,3) + (t975 / 0.2e1 + t976 / 0.2e1) * mrSges(7,2) + t1205;
t672 = -Ifges(7,3) * t434 + t1160;
t1098 = -t389 * t475 / 0.2e1 - t150 * t474 / 0.2e1 - t794 * t569 / 0.4e1 + (t649 / 0.4e1 - t1187 / 0.4e1) * t759 + (-t780 / 0.4e1 + t1186 / 0.4e1) * t434 + (t776 / 0.4e1 - t1185 / 0.4e1) * t1132 + (t672 / 0.4e1 - t1209 / 0.4e1) * t1124;
t1094 = -m(6) / 0.2e1;
t1092 = -m(7) / 0.2e1;
t1090 = -pkin(5) / 0.2e1;
t1088 = m(7) * pkin(4);
t1080 = m(7) * t90;
t1079 = m(7) * t98;
t1077 = -t328 / 0.2e1;
t1076 = t330 / 0.2e1;
t1074 = t434 / 0.2e1;
t1070 = t1114 / 0.2e1;
t1067 = -t1114 / 0.2e1;
t887 = mrSges(5,3) * t920;
t468 = t569 * mrSges(5,1) + t887;
t1066 = t468 / 0.2e1;
t1059 = t522 / 0.2e1;
t1042 = m(7) * t122;
t1041 = m(7) * t1114;
t1040 = m(7) * t1113;
t1039 = m(7) * t759;
t1037 = pkin(3) * t580;
t1029 = mrSges(5,3) * t607;
t1011 = t109 * mrSges(6,2);
t1010 = t109 * mrSges(7,3);
t1009 = t110 * mrSges(6,1);
t1008 = t110 * mrSges(7,1);
t1007 = t122 * mrSges(6,1);
t1006 = t122 * mrSges(7,1);
t1005 = t123 * mrSges(6,2);
t1004 = t123 * mrSges(7,3);
t1000 = t1114 * mrSges(6,1);
t999 = t1114 * mrSges(7,1);
t998 = t452 * mrSges(6,2);
t997 = t452 * mrSges(7,3);
t996 = t1113 * mrSges(6,1);
t995 = t1113 * mrSges(7,1);
t992 = t503 * mrSges(6,2);
t991 = t503 * mrSges(7,3);
t983 = t759 * mrSges(7,2);
t888 = mrSges(5,3) * t921;
t466 = -mrSges(5,2) * t569 + t888;
t7 = m(5) * (t282 * t292 + t283 * t293 + t933) + m(6) * (t109 * t111 + t110 * t112 + t1193) + m(7) * (t85 * t89 + t86 * t90 + t1148) + t292 * t468 + t293 * t466 + t89 * t329 + t112 * t324 + t111 * t327 + t90 * t328 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t1045 + (Ifges(3,1) - Ifges(3,2)) * t1048) * t1045 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t1048) * t1048 + (m(4) * t598 + mrSges(4,1) * t569 - mrSges(4,2) * t570) * t891 + t1181;
t977 = t7 * qJD(1);
t46 = m(7) * (t150 * t434 + t569 * t85) + t569 * t329 + t434 * t239;
t974 = qJD(1) * t46;
t10 = m(5) * (t282 * t304 + t283 * t305 + t933) + m(6) * (t109 * t116 + t110 * t117 + t1193) + m(7) * (t85 * t97 + t86 * t98 + t1148) + t304 * t468 + t305 * t466 + t97 * t329 + t117 * t324 + t116 * t327 + t98 * t328 + t1181;
t973 = t10 * qJD(1);
t972 = t109 * t1124;
t971 = t110 * t759;
t613 = t110 * t1002 + t85 * t1003 + t795 * t1053 + t649 * t1074 + t1187 * t1075 - t109 * t704 + t86 * t1161 + t672 * t1167 + t1209 * t700 + t150 * t670 + t389 * t671;
t741 = t579 * t570;
t13 = t613 + t696 * t920 / 0.2e1 + t697 * t829 + t745 * t1053 + t505 * t741 + (-m(6) * t389 - t240) * t896 + t1117 * t643 + (-t468 + t887) * t283 + (t466 - t888) * t282 + (m(6) * t110 + t1109) * t123 + (-m(6) * t109 + t1108) * t122 + t1153 * t570 ^ 2;
t966 = t13 * qJD(1);
t14 = t109 * t1109 + t110 * t1108 + t1117 * t669 + t613;
t965 = t14 * qJD(1);
t959 = t304 * t607;
t958 = t305 * t608;
t955 = t389 * t607;
t953 = t425 * t237;
t949 = t452 * t325;
t948 = t452 * t326;
t947 = t1114 * t323;
t946 = t1114 * t330;
t944 = t457 * t473;
t943 = t459 * t150;
t941 = t459 * t477;
t940 = t473 * t150;
t937 = t503 * t325;
t936 = t503 * t326;
t930 = t1113 * t323;
t929 = t1113 * t330;
t919 = t577 * t238;
t917 = t588 * t759;
t916 = t596 * t455;
t915 = t596 * t580;
t906 = 0.2e1 * t867;
t900 = mrSges(6,3) * t1036;
t886 = mrSges(7,2) * t951;
t885 = t117 * t1124 * mrSges(6,3);
t884 = mrSges(5,3) * t958;
t883 = t588 * t1003;
t882 = m(7) * t1059;
t881 = -t1036 / 0.2e1;
t879 = t1086 + mrSges(7,1) / 0.2e1;
t878 = -mrSges(7,3) / 0.2e1 + mrSges(6,2) / 0.2e1;
t873 = t1029 / 0.2e1;
t865 = -t983 / 0.2e1;
t864 = t982 / 0.2e1;
t859 = t1046 * t1124;
t852 = t972 / 0.2e1;
t851 = -t972 / 0.2e1;
t850 = t971 / 0.2e1;
t849 = -t971 / 0.2e1;
t830 = t1124 * t1054;
t823 = t1067 + t1070;
t820 = -t588 + t1036;
t817 = t434 * t900;
t814 = t1114 * t523 + t452 * t522;
t812 = t1113 * t523 + t503 * t522;
t83 = mrSges(6,3) * t850;
t299 = t457 * t474;
t450 = t597 * t475;
t796 = -t299 - t450 + t1130;
t790 = t1170 - t823;
t635 = -t1021 / 0.2e1 - t1020 / 0.2e1 - t1016 / 0.2e1 + t1015 / 0.2e1 - t1014 / 0.2e1 - t1012 / 0.2e1;
t612 = t1098 + t86 * t866 + t635 + t85 * t983 / 0.2e1;
t609 = t83 + (t852 - t1152) * mrSges(6,3) + t612 + t1101 + t1211;
t765 = t607 * t793;
t738 = t765 / 0.2e1;
t769 = t1114 * t123 + t452 * t122;
t615 = (-t1114 * t109 - t452 * t110 + (-t577 * t920 + t955) * pkin(4) + t769) * t1094 + (t1114 * t86 + t425 * t643 - t452 * t85 + t769 + t943) * t1092 + t327 * t1070 + t1114 * t1077 - t596 * t741 / 0.2e1 + t466 * t738 + t586 * t1066 + t793 * t1116 + t1131 * t452 / 0.2e1;
t1 = t609 + t615 + t1102 + t1138;
t254 = t425 * t459;
t565 = t577 * t1035;
t630 = t714 + t941 + t1120;
t21 = m(6) * t565 + m(7) * t254 + t1122 + t630 + t915;
t772 = -t1 * qJD(1) + t21 * qJD(2);
t278 = t425 * t473;
t28 = m(7) * t278 + t1122 + t1155;
t764 = qJ(6) * t1076 + t1090 * t326;
t610 = t612 + t764 + (t849 + t851) * mrSges(7,2) + t1118;
t624 = (t425 * t669 + t940) * t1091 + t327 * t1067 + (t1133 * t1091 + t328 / 0.2e1) * t1114 + (t1091 * t1162 + t1168) * t452;
t737 = (-pkin(5) * t90 + qJ(6) * t89) * t1091;
t8 = t610 + t737 - t624 + t1102;
t771 = -t8 * qJD(1) + t28 * qJD(2);
t726 = t1128 * t1047;
t45 = (mrSges(5,3) * t726 + t1107) * pkin(2) + m(7) * (t425 * t889 + t814) + m(6) * (t577 * t889 + t814) + m(5) * (t1194 * t793 + t596 * t889) + t1208;
t739 = t465 * t586;
t614 = t523 * t1168 - t1205 - t884 / 0.2e1 + t116 * t864 + t98 * t865 + t97 * t866 + t304 * t873 + t467 * t738 - t947 / 0.2e1 + (t1114 * t97 + t150 * t889 + t452 * t98 + t522 * t86 + t523 * t85 + t1147) * t1092 + t466 * t774 - t946 / 0.2e1 + pkin(2) * t861 * t1066 + t522 * t1077 + t327 * t1059 - m(5) * (t505 * t889 + t1144 + (t283 * t893 + t305 * t793) * t608 + (-t282 * t893 - t304 * t793) * t607) / 0.2e1 - t739 / 0.2e1 - t919 / 0.2e1 - t1183 - t948 / 0.2e1 + t949 / 0.2e1 + t916 / 0.2e1 - (-t456 + t240 + t239) * t889 / 0.2e1 - t953 / 0.2e1 - t885 / 0.2e1 + (-t522 * t109 + t523 * t110 + t1114 * t117 - t452 * t116 + t389 * t889 + t1191) * t1094;
t6 = t929 / 0.2e1 + t930 / 0.2e1 + t936 / 0.2e1 - t937 / 0.2e1 + t614 + t1100 + t1183;
t770 = -t6 * qJD(1) + t45 * qJD(2);
t768 = t1113 * t123 + t503 * t122;
t381 = t759 * t477;
t182 = -t1039 * t425 - t381;
t760 = t1055 * t239 + t1074 * t477;
t638 = (-t830 + t436 / 0.2e1) * mrSges(7,2) - t980 / 0.2e1 + t760;
t134 = t759 * t150;
t710 = (t1114 * t569 + t425 * t434 - t134) * t1091;
t39 = t710 - t1080 / 0.2e1 + t638;
t766 = -qJD(1) * t39 - qJD(2) * t182;
t763 = m(7) * (-pkin(5) * t1114 - qJ(6) * t452);
t762 = m(7) * (-pkin(5) * t1113 - qJ(6) * t503);
t761 = m(7) * (-pkin(5) * t522 + qJ(6) * t523);
t740 = t381 + (t425 + t457) * t1039 / 0.2e1;
t736 = (-pkin(5) * t98 + qJ(6) * t97) * t1091;
t279 = t457 * t459;
t578 = t597 * t1035;
t633 = t299 / 0.2e1 + t416 / 0.2e1 - t417 / 0.2e1 + t450 / 0.2e1 + t952 / 0.2e1 + t918 / 0.2e1 + t714;
t622 = t633 + (t565 + t578) * t1093 + (t254 + t279) * t1091 + t941;
t17 = t1113 * t865 + t782 * t1049 + t788 * t1050 + t478 * t1035 - t1037 / 0.2e1 + t915 / 0.2e1 + t622 + t1106 * mrSges(7,2) + (t1189 / 0.2e1 + t1106) * mrSges(6,3) - t1153 - t1188;
t22 = m(6) * t578 + m(7) * t279 + t938 * mrSges(6,3) - t1037 - t416 + t630 - t796;
t616 = (-t1113 * t109 - t503 * t110 + (-t597 * t920 + t955) * pkin(4) + t768) * t1094 + (t1113 * t86 + t457 * t643 - t503 * t85 + t768 + t943) * t1092 + t327 * t1169 + t328 * t1170 + pkin(3) * t741 / 0.2e1 + t603 * t1066 + t1131 * t1062 + (t1050 * t466 + t1116) * pkin(9);
t3 = t609 + t616 + t1103 + t1139;
t722 = -t3 * qJD(1) + t17 * qJD(2) + t22 * qJD(3);
t623 = (t1133 * t1113 + t1162 * t503 + t457 * t669 + t940) * t1091 + t327 * t1170 + t328 * t1169 + t1131 * t1063;
t11 = t610 + t736 - t623 + t1103;
t619 = t633 + (t278 + t944) * t1091 + t318 + t1113 * t864 + t503 * t866;
t19 = t878 * t523 + t879 * t522 - t761 / 0.2e1 + t619;
t31 = -m(7) * t944 - t1130 - t1155 + t796;
t720 = -t11 * qJD(1) + t19 * qJD(2) - t31 * qJD(3);
t521 = t595 * t984;
t715 = t521 / 0.2e1;
t126 = t882 + t740;
t196 = -t1039 * t457 - t381;
t709 = (t1113 * t569 + t457 * t434 - t134) * t1091;
t41 = t709 - t1079 / 0.2e1 + t638;
t713 = -qJD(1) * t41 + qJD(2) * t126 - qJD(3) * t196;
t405 = mrSges(7,2) * t700;
t677 = t704 * t892;
t618 = (t588 * t109 + t595 * t110 + (t1046 * t85 + t606 * t86) * pkin(4)) * t1091 - t1011 / 0.2e1 + t1010 / 0.2e1 - t1009 / 0.2e1 - t1008 / 0.2e1 + t883 / 0.2e1 + t327 * t881 + t328 * t880 + t595 * t405 + t817 / 0.2e1 - t677 / 0.2e1 + t1131 * t805;
t625 = (-pkin(5) * t122 + qJ(6) * t123) * t1092 + t1007 / 0.2e1 + t1006 / 0.2e1 + t1005 / 0.2e1 - t1004 / 0.2e1 - t886 / 0.2e1 + pkin(5) * t405;
t16 = t618 + t625;
t676 = -t1136 * t1036 + t1134 * t892;
t352 = -(t1046 * t588 + t595 * t606) * t1088 - t676;
t628 = (-(t1090 - t892 / 0.2e1) * t1124 - (-qJ(6) / 0.2e1 + t1051 + t881) * t759) * mrSges(7,2) + t715;
t679 = (t1114 * t1154 + t452 * t820) * t1091;
t37 = -t763 / 0.2e1 + t679 + t628 + t1136 * t823;
t678 = (t1113 * t1154 + t503 * t820) * t1091;
t43 = t678 - t762 / 0.2e1 + t628 + t1134 * t1207 + t1136 * (t1170 + t1169);
t688 = t16 * qJD(1) + t37 * qJD(2) + t43 * qJD(3) - t352 * qJD(4);
t629 = (t588 * t569 + t85) * t1091 + t666;
t48 = t694 - t1042 / 0.2e1 + t629;
t576 = m(7) * t588 + mrSges(7,3);
t687 = -qJD(1) * t48 - qJD(4) * t576 + t1204;
t682 = mrSges(7,2) * t775 + t794;
t665 = -mrSges(7,2) * t830 - t1001 / 0.2e1 + t980 / 0.2e1 + t760;
t634 = -mrSges(6,3) * pkin(4) * t859 - mrSges(7,2) * t917 - t759 * t900 + t521 + t778 + t794;
t631 = mrSges(6,3) * t849 + t85 * t865 + t86 * t867 - t1098;
t627 = qJ(6) * t865 + pkin(5) * t866 + (-t917 / 0.2e1 + (t859 / 0.2e1 + t606 * t1056) * pkin(4)) * mrSges(7,2) + t715 + t794;
t621 = t83 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t570 + t631 - (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t435 - (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t436 + t764 - t1118 + (t850 + t852) * mrSges(7,2);
t617 = t635 + t631 + t1101 + (t873 - t1029 / 0.2e1) * t283 + (t851 + t1152) * mrSges(6,3) - t1211;
t544 = mrSges(7,3) + (qJ(6) + 0.2e1 * t880) * m(7);
t310 = t906 + t1040;
t261 = t906 + t1041;
t212 = t1040 / 0.2e1 + m(7) * t1169 + t906;
t164 = t1041 / 0.2e1 + m(7) * t1070 + t906;
t127 = t882 - t740;
t49 = t405 + t1043 / 0.2e1 + t632;
t47 = t405 + t1042 / 0.2e1 + t629;
t42 = t678 + t992 / 0.2e1 - t991 / 0.2e1 - t996 / 0.2e1 - t995 / 0.2e1 + t627 - t879 * t1113 + t878 * t503 + t762 / 0.2e1;
t40 = t709 + t1079 / 0.2e1 + t665;
t38 = t710 + t1080 / 0.2e1 + t665;
t35 = t679 + t998 / 0.2e1 - t997 / 0.2e1 - t999 / 0.2e1 - t1000 / 0.2e1 + t627 + t763 / 0.2e1 + t878 * t452 - t879 * t1114;
t20 = t761 / 0.2e1 + t619 + t1111;
t18 = (-t1062 * t1124 - (t1169 + t790) * t759) * mrSges(7,2) + (-t1124 * t1207 - t790 * t759) * mrSges(6,3) + (t596 / 0.2e1 - pkin(3) / 0.2e1) * t580 + t622 + t1120 + t1188;
t15 = t618 - t625 + t795;
t12 = t736 + t621 + t623 + t1104;
t9 = t737 + t621 + t624 + t1105;
t5 = -t614 + t636 + (t1076 + t323 / 0.2e1) * t1113 + (t579 / 0.2e1 - mrSges(4,1) / 0.2e1) * t1125 + (t1078 - t325 / 0.2e1) * t503 + t1100;
t4 = -t616 + t617 + t1104 + t1139;
t2 = -t615 + t617 + t1105 + t1138;
t23 = [qJD(2) * t7 + qJD(3) * t10 + qJD(4) * t13 + qJD(5) * t14 + qJD(6) * t46, t977 + (m(4) * (-t1044 * t505 - t1047 * t1125) * pkin(2) + t1206 + t947 + m(7) * (t1114 * t89 + t452 * t90 + t1147) - t467 * t765 + t946 + (-t970 + t969) * mrSges(6,3) + (t975 + t976) * mrSges(7,2) + m(5) * (-t292 * t765 + t293 * t586 + t1144) + t739 + t919 + t948 - t949 - t916 + (t569 * t893 + t570 * t889) * mrSges(4,3) + t953 + t1121 * mrSges(5,3) + mrSges(3,2) * t890 + m(6) * (-t111 * t452 + t1114 * t112 + t1191) - mrSges(3,1) * t894 - Ifges(3,6) * t1045 + Ifges(3,5) * t1048) * qJD(2) + t5 * qJD(3) + t2 * qJD(4) + t9 * qJD(5) + t38 * qJD(6), t973 + t5 * qJD(2) + (t884 + (t909 - t910) * pkin(9) + t945 + t930 + t1038 + t929 + t914 + t936 - t937 + t885 - mrSges(5,3) * t959 - t116 * t982 + t97 * t984 + t98 * t983 + t1206) * qJD(3) + t4 * qJD(4) + t12 * qJD(5) + t40 * qJD(6) + 0.2e1 * ((-t1164 + (t958 - t959) * pkin(9)) * t1095 + (t1113 * t117 - t116 * t503 + t1190) * t1093 + (t1113 * t97 + t503 * t98 + t1146) * t1091) * qJD(3), t966 + t2 * qJD(2) + t4 * qJD(3) + (Ifges(5,5) * t921 + Ifges(5,6) * t920 + t595 * t1161 + m(7) * (t122 * t595 + t123 * t588) + t883 + t1004 - t1006 + (-t1046 * t122 + t123 * t606) * t1089 - t1005 - t1007 - t283 * mrSges(5,1) - t282 * mrSges(5,2) + t817 - t677 + t795) * qJD(4) + t15 * qJD(5) + t47 * qJD(6), t965 + t9 * qJD(2) + t12 * qJD(3) + t15 * qJD(4) + (-pkin(5) * t1161 + m(7) * (-pkin(5) * t110 + qJ(6) * t109) + t886 + t1010 - t1008 - t1011 - t1009 + t795) * qJD(5) + t49 * qJD(6), qJD(2) * t38 + qJD(3) * t40 + qJD(4) * t47 + qJD(5) * t49 + t974; -qJD(3) * t6 - qJD(4) * t1 - qJD(5) * t8 + qJD(6) * t39 - t977, qJD(3) * t45 + qJD(4) * t21 + qJD(5) * t28 + qJD(6) * t182, t18 * qJD(4) + t20 * qJD(5) + t127 * qJD(6) + t770 + (m(7) * (t457 * t889 + t812) + m(6) * (t597 * t889 + t812) + (m(5) * (-pkin(3) * t1044 + pkin(9) * t726) + t1107) * pkin(2) + mrSges(5,3) * t1194 + t1208) * qJD(3), t18 * qJD(3) + ((-t1046 * t1114 - t452 * t606) * t1089 + m(7) * (t1114 * t595 - t452 * t588) + t998 - t997 - t999 - t1000 + t634 - mrSges(5,1) * t586 + mrSges(5,2) * t765) * qJD(4) + t35 * qJD(5) + t164 * qJD(6) + t772, t20 * qJD(3) + t35 * qJD(4) + (t1110 * t452 + t1112 * t1114 + t682) * qJD(5) + t261 * qJD(6) + t771, qJD(3) * t127 + qJD(4) * t164 + qJD(5) * t261 - t766; qJD(2) * t6 - qJD(4) * t3 - qJD(5) * t11 + qJD(6) * t41 - t973, qJD(4) * t17 + qJD(5) * t19 - qJD(6) * t126 - t770, qJD(4) * t22 - qJD(5) * t31 + qJD(6) * t196 ((-t1046 * t1113 - t503 * t606) * t1089 + m(7) * (t1113 * t595 - t503 * t588) + t992 - t991 - t996 - t995 + t634 + t579 * pkin(9)) * qJD(4) + t42 * qJD(5) + t212 * qJD(6) + t722, t42 * qJD(4) + (t1110 * t503 + t1112 * t1113 + t682) * qJD(5) + t310 * qJD(6) + t720, qJD(4) * t212 + qJD(5) * t310 - t713; qJD(2) * t1 + qJD(3) * t3 + qJD(5) * t16 + qJD(6) * t48 - t966, -qJD(3) * t17 + qJD(5) * t37 - t461 - t772, qJD(5) * t43 - t461 - t722, -qJD(5) * t352 + qJD(6) * t576 ((-pkin(5) * t606 + qJ(6) * t1046) * t1088 + t676) * qJD(5) + t544 * qJD(6) + t688, qJD(5) * t544 - t687; qJD(2) * t8 + qJD(3) * t11 - qJD(4) * t16 + qJD(6) * t50 - t965, -qJD(3) * t19 - qJD(4) * t37 + t461 - t771, -qJD(4) * t43 + t461 - t720, -t688 + t1143, t1143, t1115; -qJD(2) * t39 - qJD(3) * t41 - qJD(4) * t48 - qJD(5) * t50 - t974, qJD(3) * t126 + t1182 + t766, t1182 + t713, t687 - t1150, -t1115, 0;];
Cq  = t23;
