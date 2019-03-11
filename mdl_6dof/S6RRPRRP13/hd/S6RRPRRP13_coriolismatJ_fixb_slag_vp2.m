% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:43
% EndTime: 2019-03-09 12:56:29
% DurationCPUTime: 25.83s
% Computational Cost: add. (28752->1262), mult. (68210->1662), div. (0->0), fcn. (68394->8), ass. (0->592)
t730 = cos(pkin(6));
t732 = sin(qJ(4));
t735 = cos(qJ(4));
t729 = sin(pkin(6));
t736 = cos(qJ(2));
t962 = t729 * t736;
t615 = t730 * t735 - t732 * t962;
t731 = sin(qJ(5));
t733 = sin(qJ(2));
t734 = cos(qJ(5));
t955 = t733 * t734;
t448 = -t615 * t731 + t729 * t955;
t614 = t730 * t732 + t735 * t962;
t963 = t729 * t733;
t449 = t615 * t734 + t731 * t963;
t997 = t449 * Ifges(7,4);
t181 = t448 * Ifges(7,2) + Ifges(7,6) * t614 + t997;
t998 = t449 * Ifges(6,4);
t182 = t448 * Ifges(6,2) + Ifges(6,6) * t614 + t998;
t1175 = t182 + t181;
t1196 = t1175 / 0.4e1;
t1185 = Ifges(6,5) + Ifges(7,5);
t1184 = Ifges(6,6) + Ifges(7,6);
t727 = t731 ^ 2;
t728 = t734 ^ 2;
t947 = t727 + t728;
t1188 = mrSges(6,3) + mrSges(7,3);
t251 = Ifges(7,1) * t448 - t997;
t252 = Ifges(6,1) * t448 - t998;
t1195 = t252 + t251;
t1133 = pkin(2) + pkin(9);
t1132 = pkin(3) + pkin(8);
t1041 = pkin(1) * t730;
t707 = t736 * t1041;
t518 = -t1132 * t963 + t707;
t433 = -t1133 * t730 - t518;
t856 = -qJ(3) * t733 - pkin(1);
t469 = (-t1133 * t736 + t856) * t729;
t234 = t433 * t735 - t732 * t469;
t196 = -pkin(4) * t963 - t234;
t135 = -pkin(5) * t448 + t196;
t1045 = m(7) * t135;
t247 = -mrSges(7,1) * t448 + mrSges(7,2) * t449;
t1194 = t247 + t1045;
t235 = t433 * t732 + t469 * t735;
t197 = pkin(10) * t963 + t235;
t704 = pkin(3) * t962;
t618 = pkin(8) * t962 + t733 * t1041;
t719 = t730 * qJ(3);
t911 = t719 + t618;
t849 = t704 + t911;
t244 = pkin(4) * t614 - pkin(10) * t615 + t849;
t103 = -t197 * t731 + t734 * t244;
t104 = t197 * t734 + t244 * t731;
t1147 = -m(7) / 0.2e1;
t1149 = -m(6) / 0.2e1;
t1034 = pkin(10) * t735;
t1038 = pkin(4) * t732;
t664 = qJ(3) - t1034 + t1038;
t640 = t734 * t664;
t855 = t1133 * t731 + pkin(5);
t953 = t734 * t735;
t915 = qJ(6) * t953;
t412 = t732 * t855 + t640 - t915;
t956 = t732 * t1133;
t509 = t731 * t664 - t734 * t956;
t960 = t731 * t735;
t446 = -qJ(6) * t960 + t509;
t508 = t731 * t956 + t640;
t954 = t733 * t735;
t916 = t729 * t954;
t587 = t730 * t734 + t731 * t916;
t588 = t730 * t731 - t734 * t916;
t1037 = pkin(5) * t731;
t854 = t1133 + t1037;
t642 = t854 * t735;
t720 = t730 * mrSges(4,3);
t75 = -qJ(6) * t449 + t103;
t62 = pkin(5) * t614 + t75;
t76 = qJ(6) * t448 + t104;
t820 = -t62 * t734 - t731 * t76;
t652 = mrSges(6,1) * t732 - mrSges(6,3) * t953;
t1071 = t652 / 0.2e1;
t984 = t732 * mrSges(7,1);
t651 = -mrSges(7,3) * t953 + t984;
t1073 = t651 / 0.2e1;
t866 = t1071 + t1073;
t932 = mrSges(7,3) * t960;
t983 = t732 * mrSges(7,2);
t647 = -t932 - t983;
t1078 = t647 / 0.2e1;
t933 = mrSges(6,3) * t960;
t648 = -mrSges(6,2) * t732 - t933;
t867 = t648 / 0.2e1 + t1078;
t910 = 0.2e1 * t719 + t618;
t959 = t732 * t733;
t917 = t729 * t959;
t988 = t615 * mrSges(5,2);
t992 = t614 * mrSges(5,1);
t1193 = -t720 - m(4) * t910 / 0.2e1 - m(5) * (t704 + t910) / 0.2e1 + (t734 * t103 + t731 * t104 + t508 * t587 + t509 * t588 - t916 * t956) * t1149 + (t412 * t587 + t446 * t588 - t642 * t917 - t820) * t1147 - t992 / 0.2e1 - t988 / 0.2e1 - t867 * t588 - t866 * t587;
t1192 = m(7) * t76;
t1190 = m(7) * (t412 * t734 + t446 * t731);
t1044 = m(7) * t642;
t721 = t731 * mrSges(7,1);
t722 = t734 * mrSges(7,2);
t1166 = t722 + t721;
t623 = t1166 * t735;
t1083 = t623 / 0.2e1;
t1189 = pkin(5) * (t1044 / 0.2e1 + t1083);
t1187 = -Ifges(4,4) + Ifges(3,5);
t1186 = Ifges(4,5) - Ifges(3,6);
t1183 = -Ifges(6,3) - Ifges(7,3);
t617 = pkin(8) * t963 - t707;
t1182 = t617 * mrSges(3,2);
t922 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t1181 = t614 * t922;
t1180 = t849 * (mrSges(5,1) * t735 - mrSges(5,2) * t732);
t921 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t1177 = t921 * t614;
t1176 = t922 * t732;
t1144 = m(7) * pkin(5);
t919 = -mrSges(7,1) - t1144;
t440 = Ifges(7,4) * t448;
t183 = t449 * Ifges(7,1) + Ifges(7,5) * t614 + t440;
t441 = Ifges(6,4) * t448;
t184 = t449 * Ifges(6,1) + Ifges(6,5) * t614 + t441;
t1174 = t184 + t183;
t1001 = t448 * mrSges(7,3);
t990 = t614 * mrSges(7,2);
t317 = -t990 + t1001;
t1002 = t448 * mrSges(6,3);
t318 = -mrSges(6,2) * t614 + t1002;
t1173 = t317 + t318;
t991 = t614 * mrSges(7,1);
t999 = t449 * mrSges(7,3);
t319 = t991 - t999;
t1000 = t449 * mrSges(6,3);
t320 = mrSges(6,1) * t614 - t1000;
t1172 = t320 + t319;
t725 = Ifges(7,4) * t734;
t1160 = -Ifges(7,2) * t731 + t725;
t576 = Ifges(7,6) * t732 + t1160 * t735;
t726 = Ifges(6,4) * t734;
t1159 = -Ifges(6,2) * t731 + t726;
t578 = Ifges(6,6) * t732 + t1159 * t735;
t1171 = t576 + t578;
t1019 = Ifges(7,4) * t731;
t681 = Ifges(7,1) * t734 - t1019;
t580 = Ifges(7,5) * t732 + t681 * t735;
t1020 = Ifges(6,4) * t731;
t683 = Ifges(6,1) * t734 - t1020;
t582 = Ifges(6,5) * t732 + t683 * t735;
t1170 = t580 + t582;
t670 = mrSges(6,1) * t731 + mrSges(6,2) * t734;
t624 = t735 * t670;
t1169 = t624 + t623;
t1168 = t647 + t648;
t1167 = t651 + t652;
t1012 = Ifges(7,6) * t731;
t723 = Ifges(7,5) * t734;
t822 = t723 - t1012;
t1013 = Ifges(6,6) * t731;
t724 = Ifges(6,5) * t734;
t823 = t724 - t1013;
t1165 = -t822 - t823;
t1162 = -t1184 * t449 + t1185 * t448;
t969 = t449 * t731;
t970 = t448 * t734;
t1161 = t969 + t970;
t680 = Ifges(7,1) * t731 + t725;
t682 = Ifges(6,1) * t731 + t726;
t1158 = t104 * mrSges(6,3) + t76 * mrSges(7,3);
t1157 = t103 * mrSges(6,3) + t62 * mrSges(7,3);
t671 = Ifges(7,5) * t731 + Ifges(7,6) * t734;
t673 = Ifges(6,5) * t731 + Ifges(6,6) * t734;
t1156 = t673 / 0.2e1 + t671 / 0.2e1;
t869 = t582 / 0.4e1 + t580 / 0.4e1;
t1155 = t578 / 0.2e1 + t576 / 0.2e1;
t871 = -t578 / 0.4e1 - t576 / 0.4e1;
t249 = -Ifges(7,2) * t449 + t440;
t250 = -Ifges(6,2) * t449 + t441;
t1154 = t250 + t249 + t1174;
t1028 = mrSges(5,3) * t615;
t1048 = t735 / 0.2e1;
t1055 = t732 / 0.2e1;
t1057 = -t732 / 0.2e1;
t248 = -mrSges(6,1) * t448 + mrSges(6,2) * t449;
t989 = t614 * mrSges(5,3);
t477 = -mrSges(5,2) * t963 - t989;
t936 = mrSges(5,1) * t963;
t799 = t936 - t1028;
t778 = t732 * t799;
t1153 = -t778 / 0.2e1 + t1028 * t1057 + t248 * t1055 + (t989 + t477) * t1048;
t1152 = t729 ^ 2;
t1151 = m(4) / 0.2e1;
t1150 = m(5) / 0.2e1;
t1148 = m(6) / 0.2e1;
t1146 = m(7) / 0.2e1;
t1145 = m(6) * pkin(4);
t1143 = mrSges(5,1) / 0.2e1;
t1142 = -mrSges(6,1) / 0.2e1;
t1141 = -mrSges(6,2) / 0.2e1;
t1140 = mrSges(6,2) / 0.2e1;
t1139 = -mrSges(7,2) / 0.2e1;
t1138 = mrSges(7,2) / 0.2e1;
t1137 = -mrSges(7,3) / 0.2e1;
t1136 = mrSges(7,3) / 0.2e1;
t1135 = Ifges(5,1) / 0.2e1;
t400 = pkin(4) * t615 + pkin(10) * t614;
t130 = -t234 * t731 + t734 * t400;
t964 = t614 * t734;
t97 = pkin(5) * t615 + qJ(6) * t964 + t130;
t1134 = -t97 / 0.2e1;
t1131 = mrSges(6,3) * pkin(10);
t977 = qJ(3) * t736;
t512 = (t1133 * t733 - t977) * t729;
t519 = t704 + t618;
t304 = t735 * t512 + t732 * t519;
t264 = pkin(10) * t962 + t304;
t687 = pkin(4) * t735 + pkin(10) * t732;
t349 = t707 + (-t687 - t1132) * t963;
t132 = -t264 * t731 + t734 * t349;
t540 = (t731 * t736 + t732 * t955) * t729;
t105 = -pkin(5) * t916 - qJ(6) * t540 + t132;
t1130 = t105 / 0.2e1;
t1129 = -t132 / 0.2e1;
t133 = t734 * t264 + t731 * t349;
t1128 = t133 / 0.2e1;
t1127 = -t135 / 0.2e1;
t1126 = t135 / 0.2e1;
t1125 = -t196 / 0.2e1;
t246 = mrSges(6,1) * t449 + mrSges(6,2) * t448;
t1124 = t246 / 0.2e1;
t1123 = -t247 / 0.2e1;
t1122 = -t248 / 0.2e1;
t1121 = -t304 / 0.2e1;
t1120 = -t317 / 0.2e1;
t1119 = -t318 / 0.2e1;
t1118 = t319 / 0.2e1;
t1117 = t320 / 0.2e1;
t539 = (-t731 * t959 + t734 * t736) * t729;
t326 = -mrSges(6,1) * t539 + mrSges(6,2) * t540;
t1116 = -t326 / 0.2e1;
t965 = t614 * t731;
t384 = -mrSges(6,2) * t615 + mrSges(6,3) * t965;
t1115 = -t384 / 0.2e1;
t385 = mrSges(7,1) * t615 + mrSges(7,3) * t964;
t1114 = t385 / 0.2e1;
t1113 = -t412 / 0.2e1;
t1112 = t412 / 0.2e1;
t414 = mrSges(6,2) * t916 + mrSges(6,3) * t539;
t1111 = t414 / 0.2e1;
t415 = -mrSges(7,1) * t916 - mrSges(7,3) * t540;
t1110 = t415 / 0.2e1;
t416 = -mrSges(6,1) * t916 - mrSges(6,3) * t540;
t1109 = -t416 / 0.2e1;
t644 = t734 * t687;
t958 = t732 * t734;
t417 = qJ(6) * t958 + t735 * t855 + t644;
t1108 = -t417 / 0.2e1;
t445 = t508 - t915;
t1107 = -t445 / 0.2e1;
t1106 = -t446 / 0.2e1;
t1105 = t446 / 0.2e1;
t1104 = -t448 / 0.2e1;
t1100 = -t508 / 0.2e1;
t1099 = t509 / 0.2e1;
t952 = t735 * t1133;
t536 = t731 * t952 + t644;
t1098 = -t536 / 0.2e1;
t1087 = t614 / 0.2e1;
t1086 = t614 / 0.4e1;
t982 = t734 * mrSges(6,1);
t986 = t731 * mrSges(6,2);
t832 = t982 - t986;
t620 = t832 * t735;
t1084 = -t620 / 0.2e1;
t1082 = -t642 / 0.2e1;
t961 = t731 * t732;
t645 = -mrSges(7,2) * t735 + mrSges(7,3) * t961;
t1081 = -t645 / 0.2e1;
t646 = -mrSges(6,2) * t735 + mrSges(6,3) * t961;
t1080 = -t646 / 0.2e1;
t1079 = -t647 / 0.2e1;
t1077 = -t648 / 0.2e1;
t649 = mrSges(7,1) * t735 + mrSges(7,3) * t958;
t1075 = -t649 / 0.2e1;
t1074 = -t651 / 0.2e1;
t1072 = -t652 / 0.2e1;
t1030 = -qJ(6) - pkin(10);
t665 = t1030 * t731;
t1070 = t665 / 0.2e1;
t985 = t731 * mrSges(7,2);
t666 = -t734 * mrSges(7,1) + t985;
t1069 = t666 / 0.2e1;
t668 = t1030 * t734;
t1068 = t668 / 0.2e1;
t1067 = -t668 / 0.2e1;
t1066 = t1166 / 0.2e1;
t1065 = t670 / 0.2e1;
t1036 = pkin(5) * t734;
t710 = -pkin(4) - t1036;
t1062 = t710 / 0.2e1;
t1061 = -t730 / 0.2e1;
t1060 = -t731 / 0.2e1;
t1059 = t731 / 0.2e1;
t1056 = -t732 / 0.4e1;
t1053 = -t734 / 0.2e1;
t1051 = t734 / 0.2e1;
t1049 = -t735 / 0.2e1;
t1047 = t735 / 0.4e1;
t1046 = t1133 / 0.2e1;
t1043 = m(7) * t710;
t1042 = m(7) * t732;
t1040 = pkin(2) * t730;
t303 = -t732 * t512 + t519 * t735;
t263 = -pkin(4) * t962 - t303;
t1039 = pkin(4) * t263;
t1035 = pkin(5) * t735 ^ 2;
t1029 = -t62 + t75;
t1025 = Ifges(3,4) * t733;
t1024 = Ifges(3,4) * t736;
t1023 = Ifges(5,4) * t615;
t1022 = Ifges(5,4) * t732;
t1021 = Ifges(5,4) * t735;
t1018 = Ifges(5,5) * t732;
t1015 = Ifges(4,6) * t733;
t1014 = Ifges(4,6) * t736;
t1011 = Ifges(5,3) * t736;
t1010 = Ifges(6,3) * t615;
t1009 = Ifges(6,3) * t735;
t1008 = Ifges(7,3) * t615;
t1007 = Ifges(7,3) * t735;
t1006 = qJ(3) * mrSges(5,1);
t109 = qJ(6) * t539 + t133;
t176 = -pkin(5) * t539 + t263;
t272 = Ifges(7,4) * t540 + Ifges(7,2) * t539 - Ifges(7,6) * t916;
t273 = Ifges(6,4) * t540 + Ifges(6,2) * t539 - Ifges(6,6) * t916;
t274 = Ifges(7,1) * t540 + Ifges(7,4) * t539 - Ifges(7,5) * t916;
t275 = Ifges(6,1) * t540 + Ifges(6,4) * t539 - Ifges(6,5) * t916;
t994 = t540 * mrSges(7,2);
t995 = t539 * mrSges(7,1);
t325 = t994 - t995;
t413 = mrSges(7,2) * t916 + mrSges(7,3) * t539;
t590 = (mrSges(5,1) * t736 - mrSges(5,3) * t959) * t729;
t591 = (-mrSges(5,2) * t736 + mrSges(5,3) * t954) * t729;
t765 = t735 * (-Ifges(5,2) * t614 + Ifges(5,6) * t963 + t1023);
t770 = t735 * (Ifges(7,5) * t449 + Ifges(7,6) * t448 + Ifges(7,3) * t614);
t771 = t735 * (Ifges(6,5) * t449 + Ifges(6,6) * t448 + Ifges(6,3) * t614);
t600 = Ifges(5,4) * t614;
t931 = Ifges(5,5) * t963;
t773 = t732 * (Ifges(5,1) * t615 - t600 + t931);
t802 = -pkin(2) * t736 + t856;
t779 = t1152 * t802;
t782 = t1152 * (mrSges(4,2) * t736 - mrSges(4,3) * t733);
t821 = pkin(2) * t733 - t977;
t824 = Ifges(5,6) * t735 + t1018;
t827 = Ifges(5,2) * t735 + t1022;
t830 = Ifges(5,1) * t732 + t1021;
t835 = t988 + t992;
t930 = mrSges(4,1) * t962;
t838 = -t720 - t930;
t901 = t962 / 0.2e1;
t902 = -t962 / 0.2e1;
t903 = t963 / 0.2e1;
t904 = -t963 / 0.2e1;
t3 = -t303 * t799 - t518 * t835 - ((Ifges(3,1) * t736 + t733 * t824 + t1011 - t1025) * t733 + t736 * (-Ifges(3,2) * t733 + t1024)) * t1152 / 0.2e1 + (t1186 * t733 + t1187 * t736) * t729 * t1061 + (t1180 + (t911 - t618) * mrSges(4,1) + Ifges(5,3) * t902) * t963 + ((Ifges(5,6) * t736 + t733 * t827) * t1087 - t615 * (Ifges(5,5) * t736 + t733 * t830) / 0.2e1 + (Ifges(3,1) * t733 + t1024) * t902 + (Ifges(3,2) * t736 + t1025) * t903 + (-Ifges(4,2) * t733 - t1014) * t901 + (-Ifges(4,3) * t736 - t1015) * t904 + (-pkin(8) * t733 * t618 - t729 * t802 * t821) * m(4)) * t729 + (t771 + t770) * t903 + (Ifges(5,5) * t615 - Ifges(5,6) * t614) * t902 - m(5) * (t234 * t303 + t235 * t304 + t518 * t849) - m(4) * (-t911 * t617 + (-t707 - t1040) * t618) + (t773 + t765) * t904 - m(6) * (t103 * t132 + t104 * t133 + t196 * t263) - m(7) * (t105 * t62 + t109 * t76 + t135 * t176) + t1152 * pkin(1) * (mrSges(3,1) * t733 + mrSges(3,2) * t736) - (-mrSges(4,2) * t733 - mrSges(4,3) * t736) * t779 - t234 * t590 - t235 * t591 + (t273 + t272) * t1104 - t304 * t477 - t103 * t416 - t76 * t413 - t104 * t414 - t62 * t415 - t135 * t325 - t196 * t326 - t105 * t319 - t132 * t320 - t109 * t317 - t133 * t318 + (t736 * (Ifges(4,3) * t733 - t1014) + t733 * (-Ifges(4,2) * t736 + t1015)) * t1152 / 0.2e1 - (t275 + t274) * t449 / 0.2e1 - t263 * t248 - t176 * t247 - t821 * t782 - (t617 - t1040) * t930 - t1174 * t540 / 0.2e1 - t1175 * t539 / 0.2e1 - t617 * t838 + (-t1182 + Ifges(4,4) * t901 + Ifges(3,5) * t902 + Ifges(4,5) * t904 + Ifges(3,6) * t903 + (mrSges(3,1) - mrSges(4,2)) * t618) * t730 - (t1183 * t916 + t1184 * t539 + t1185 * t540) * t614 / 0.2e1;
t1003 = t3 * qJD(1);
t996 = t508 * mrSges(6,3);
t131 = t734 * t234 + t731 * t400;
t110 = qJ(6) * t965 + t131;
t162 = -pkin(5) * t965 + t235;
t353 = t1166 * t614;
t354 = t670 * t614;
t383 = -mrSges(7,2) * t615 + mrSges(7,3) * t965;
t386 = mrSges(6,1) * t615 + mrSges(6,3) * t964;
t805 = t849 * mrSges(5,2);
t806 = t849 * mrSges(5,1);
t267 = Ifges(7,5) * t615 - t614 * t681;
t268 = Ifges(6,5) * t615 - t614 * t683;
t876 = t267 / 0.2e1 + t268 / 0.2e1;
t265 = t615 * Ifges(7,6) - t1160 * t614;
t266 = t615 * Ifges(6,6) - t1159 * t614;
t877 = t265 / 0.2e1 + t266 / 0.2e1;
t879 = -t183 / 0.2e1 - t184 / 0.2e1;
t881 = t181 / 0.2e1 + t182 / 0.2e1;
t948 = -Ifges(5,5) * t614 - Ifges(5,6) * t615;
t6 = t234 * t477 + t76 * t383 + t104 * t384 + t62 * t385 + t103 * t386 - t196 * t354 - t135 * t353 + t97 * t319 + t130 * t320 + t110 * t317 + t131 * t318 + t162 * t247 + t948 * t903 + t876 * t449 + t877 * t448 + (t248 - t936) * t235 + m(6) * (t103 * t130 + t104 * t131 + t196 * t235) + m(7) * (t110 * t76 + t135 * t162 + t62 * t97) + (Ifges(5,6) * t904 + t448 * t921 - t449 * t922 - t1023 + t806) * t615 + (Ifges(5,5) * t904 + t600 - t805 + t234 * mrSges(5,3) + (t879 + t1181) * t734 + (t881 + t1177) * t731 + (-Ifges(5,1) + Ifges(5,2) - t1183) * t615) * t614;
t993 = t6 * qJD(1);
t987 = t62 * t731;
t981 = t735 * mrSges(5,2);
t980 = t76 * t734;
t435 = t448 * mrSges(7,2);
t245 = mrSges(7,1) * t449 + t435;
t918 = m(7) * t1029;
t9 = t196 * t246 + t135 * t245 + t75 * t317 - t104 * t320 + t103 * t318 + (t918 - t319) * t76 + (t251 / 0.2e1 + t252 / 0.2e1 + t1194 * pkin(5) - t881 - t1158) * t449 + (t250 / 0.2e1 + t249 / 0.2e1 - t879 - t1157) * t448 + t1162 * t1087;
t979 = t9 * qJD(1);
t978 = -t832 - mrSges(5,1);
t35 = m(7) * (t448 * t76 - t449 * t62) - t449 * t319 + t448 * t317;
t976 = qJD(1) * t35;
t975 = t103 * t731;
t974 = t104 * t734;
t972 = t235 * t735;
t14 = t477 * t916 + (-m(5) * (t234 * t732 - t972) - t778) * t963 + (m(4) * t779 + t782) * t733 + (-m(4) * t911 - m(5) * t849 - t835 + t838) * t730 + (m(6) * t196 + t1194 + t248) * t917 + (-m(6) * t104 - t1173 - t1192) * t588 + (-m(6) * t103 - m(7) * t62 - t1172) * t587;
t973 = t14 * qJD(1);
t971 = t446 * t734;
t968 = t509 * t734;
t967 = t587 * t731;
t966 = t588 * t734;
t957 = t732 * t735;
t951 = -t412 + t445;
t537 = t731 * t687 - t734 * t952;
t163 = (t904 - t970 / 0.2e1 - t969 / 0.2e1) * t1042;
t946 = qJD(1) * t163;
t945 = qJD(4) * t731;
t944 = qJD(4) * t732;
t943 = qJD(4) * t734;
t942 = qJD(4) * t735;
t939 = t1144 / 0.2e1;
t938 = pkin(5) * t953;
t929 = t162 * t1146;
t641 = t854 * t732;
t928 = t641 * t1146;
t927 = t1042 / 0.2e1;
t926 = pkin(5) * t1083;
t925 = -mrSges(7,1) / 0.2e1 + t1142;
t924 = t1141 + t1139;
t923 = t1136 + mrSges(6,3) / 0.2e1;
t920 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t914 = -t1001 / 0.2e1;
t913 = -t1000 / 0.2e1;
t912 = -t999 / 0.2e1;
t900 = -t961 / 0.2e1;
t899 = t961 / 0.2e1;
t897 = -t960 / 0.2e1;
t896 = -t960 / 0.4e1;
t893 = -t958 / 0.2e1;
t891 = t958 / 0.2e1;
t888 = t319 * t1053;
t887 = -t953 / 0.2e1;
t885 = t953 / 0.2e1;
t884 = t953 / 0.4e1;
t882 = t952 / 0.2e1;
t880 = -t182 / 0.4e1 - t181 / 0.4e1;
t878 = t184 / 0.4e1 + t183 / 0.4e1;
t875 = t1119 + t1120;
t874 = t1117 + t1118;
t873 = -t354 / 0.2e1 - t477 / 0.2e1;
t575 = t735 * Ifges(7,6) - t1160 * t732;
t577 = t735 * Ifges(6,6) - t1159 * t732;
t872 = -t575 / 0.2e1 - t577 / 0.2e1;
t579 = t735 * Ifges(7,5) - t681 * t732;
t581 = t735 * Ifges(6,5) - t683 * t732;
t870 = t581 / 0.2e1 + t579 / 0.2e1;
t868 = t624 / 0.2e1 + t1083;
t865 = -t666 / 0.2e1 + t832 / 0.2e1;
t864 = -t671 / 0.4e1 - t673 / 0.4e1;
t863 = t723 / 0.4e1 - t1012 / 0.4e1 + t724 / 0.4e1 - t1013 / 0.4e1;
t675 = Ifges(7,2) * t734 + t1019;
t677 = Ifges(6,2) * t734 + t1020;
t862 = t675 / 0.2e1 + t677 / 0.2e1;
t861 = t675 / 0.4e1 + t677 / 0.4e1;
t860 = -t680 / 0.2e1 - t682 / 0.2e1;
t859 = t680 / 0.4e1 + t682 / 0.4e1;
t858 = -t728 / 0.2e1 - t727 / 0.2e1;
t857 = m(7) * t951;
t853 = t947 * mrSges(7,3);
t852 = t978 - t1145;
t844 = t732 * t904;
t843 = Ifges(5,2) / 0.2e1 + t920;
t842 = t924 * t734;
t841 = t921 * t731;
t840 = t247 / 0.2e1 + t1045 / 0.2e1;
t708 = mrSges(7,1) * t953;
t837 = -mrSges(7,2) * t960 + t708;
t836 = t1116 - t325 / 0.2e1 + t590 / 0.2e1;
t833 = t732 * mrSges(5,1) + t981;
t819 = t980 - t987;
t450 = qJ(6) * t961 + t537;
t621 = t1166 * t732;
t622 = t670 * t732;
t650 = mrSges(6,1) * t735 + mrSges(6,3) * t958;
t573 = t732 * Ifges(7,3) + t735 * t822;
t574 = t732 * Ifges(6,3) + t735 * t823;
t679 = -t732 * Ifges(5,2) + t1021;
t763 = -t573 / 0.2e1 - t574 / 0.2e1 + t679 / 0.2e1 + t1021 / 0.2e1 - t1006;
t684 = Ifges(5,1) * t735 - t1022;
t780 = qJ(3) * mrSges(5,2) + t684 / 0.2e1 - t1022 / 0.2e1;
t800 = t1135 - t843;
t19 = t412 * t649 + t417 * t651 + t446 * t645 + t450 * t647 + t508 * t650 + t509 * t646 + t536 * t652 + t537 * t648 - t642 * t621 - t641 * t623 + m(7) * (t412 * t417 + t446 * t450 - t641 * t642) + m(6) * (t508 * t536 + t509 * t537) + (-t1133 * t622 + t731 * t872 + t734 * t870 - t763) * t735 + (-t1133 * t624 + (-t580 / 0.2e1 - t582 / 0.2e1 + t1176) * t734 + (t732 * t921 + t1155) * t731 + (-m(6) * t1133 ^ 2 - t800) * t735 - t780) * t732;
t738 = (t110 * t446 - t135 * t641 + t162 * t642 + t412 * t97 + t417 * t62 + t450 * t76) * t1147 - t103 * t650 / 0.2e1 + t104 * t1080 + t110 * t1079 + t130 * t1072 + t131 * t1077 - t621 * t1127 - t162 * t623 / 0.2e1 - t622 * t1125 - t235 * t624 / 0.2e1 + t385 * t1113 + t319 * t1108 + t383 * t1106 + t450 * t1120 + t386 * t1100 + t509 * t1115 + t320 * t1098 + t537 * t1119 + t684 * t1086 + t62 * t1075 - t641 * t1123 - t353 * t1082 + t76 * t1081 + t97 * t1074;
t742 = t861 * t539 + t859 * t540 + (t105 * t665 - t109 * t668 + t176 * t710) * t1146 + pkin(4) * t1116 + t176 * t1069 - t263 * t832 / 0.2e1 + t303 * t1143 + t415 * t1070 + t413 * t1067 + t325 * t1062;
t760 = t536 * t103 + t537 * t104 + t508 * t130 + t509 * t131;
t767 = t274 / 0.4e1 + t275 / 0.4e1 + t105 * t1137 + mrSges(6,3) * t1129;
t768 = t272 / 0.4e1 + t273 / 0.4e1 + t109 * t1136 + mrSges(6,3) * t1128;
t814 = -t132 * t731 + t133 * t734;
t2 = (-t579 / 0.4e1 - t581 / 0.4e1) * t449 + (-t575 / 0.4e1 - t577 / 0.4e1) * t448 + t738 + (pkin(10) * t814 - t1039) * t1148 + (pkin(10) * t1109 + t614 * t871 + t767) * t731 + (pkin(10) * t1111 + t614 * t869 + t768) * t734 + t760 * t1149 + (t679 / 0.4e1 - t1006 / 0.2e1 - t574 / 0.4e1 - t573 / 0.4e1) * t615 + (qJ(3) * t1087 + t1121) * mrSges(5,2) + (t805 / 0.2e1 - t600 / 0.2e1 + t931 + t878 * t734 + t880 * t731 - (m(6) * t1125 + mrSges(5,1) * t903 + t1122) * t1133 + (-Ifges(5,2) / 0.4e1 + t1135 - Ifges(6,3) / 0.4e1 - Ifges(7,3) / 0.4e1) * t615 + (-Ifges(5,4) / 0.4e1 - t922 * t734 - t841) * t614) * t732 + (-t806 / 0.2e1 + 0.3e1 / 0.4e1 * t1023 + (-t268 / 0.4e1 - t267 / 0.4e1) * t734 + (t266 / 0.4e1 + t265 / 0.4e1) * t731 + (-Ifges(6,5) / 0.4e1 - Ifges(7,5) / 0.4e1) * t449 + (-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t448 - (t1148 * t235 + t873) * t1133 + (Ifges(5,6) + t864) * t963 + (mrSges(5,3) * t1046 + Ifges(5,1) / 0.4e1 - t843) * t614) * t735 + t742 + Ifges(5,3) * t901;
t818 = -t2 * qJD(1) + t19 * qJD(2);
t24 = t509 * t652 - t623 * t938 - t642 * t837 - t445 * t647 - t412 * t932 + (-t1133 * t620 + mrSges(6,3) * t968 + mrSges(7,3) * t971 - t1036 * t1044 + ((t675 + t677) * t1060 + (-t680 - t682) * t1053) * t735 + t1170 * t1059 + (t671 + t673) * t1055 + t1171 * t1051) * t735 + (-t933 - t648) * t508 + (t651 - t857) * t446;
t743 = t103 * t1077 + t104 * t1071 + t1056 * t1162 + t1073 * t76 + t1079 * t75 + t1082 * t245 + t1084 * t196 + t1099 * t320 + t1100 * t318 + t1105 * t319 + t1107 * t317;
t764 = t539 * t921 - t540 * t922;
t745 = mrSges(7,1) * t1130 + t109 * t1139 + t132 * mrSges(6,1) / 0.2e1 + t133 * t1141 + t764;
t783 = t614 * t671;
t784 = t614 * t673;
t785 = t449 * t680;
t786 = t449 * t682;
t787 = t448 * t675;
t788 = t448 * t677;
t804 = t920 * t963;
t4 = t708 * t1127 + (t996 / 0.2e1 + mrSges(7,3) * t1112 - t869) * t448 + pkin(5) * t1110 + (pkin(5) * t1130 + t1105 * t62 + t1106 * t75 + (t1107 + t1112) * t76) * m(7) + (mrSges(6,3) * t1099 + mrSges(7,3) * t1105 - t1189 - t871) * t449 + (t786 / 0.4e1 + t785 / 0.4e1 + t788 / 0.4e1 + t787 / 0.4e1 + t985 * t1126 + t784 / 0.4e1 + t783 / 0.4e1 - t1133 * t1124 - t804 - t840 * t1036 + (t980 / 0.2e1 - t987 / 0.2e1) * mrSges(7,3) + (t974 / 0.2e1 - t975 / 0.2e1) * mrSges(6,3) + (-t1195 / 0.4e1 + t1196) * t734) * t735 + t743 + t745 + t1154 * t960 / 0.4e1;
t817 = -t4 * qJD(1) - t24 * qJD(2);
t813 = t303 * t735 + t304 * t732;
t744 = t618 * t1151 + t813 * t1150 + (-t263 * t735 + t732 * t814) * t1148 + (-t176 * t735 + (-t105 * t731 + t109 * t734) * t732) * t1146;
t755 = t591 / 0.2e1 + (t1111 + t413 / 0.2e1) * t734 + (t1109 - t415 / 0.2e1) * t731;
t13 = (mrSges(5,1) * t1061 + t868 * t963 + t755) * t732 - t874 * t734 + t875 * t731 + (mrSges(5,2) * t1061 + t836) * t735 + t744 + t1193;
t55 = mrSges(4,3) + t1167 * t734 + t1168 * t731 + (m(5) + m(4)) * qJ(3) + m(6) * (t508 * t734 + t509 * t731) + t1190 + t833;
t816 = qJD(1) * t13 - qJD(2) * t55;
t803 = t939 - t925;
t753 = t587 * t803 + t588 * t924;
t801 = -t319 / 0.2e1 + t918 / 0.2e1;
t15 = (t1124 + t245 / 0.2e1 + t449 * t939) * t735 + ((-t448 * t923 - t875) * t731 + (t449 * t923 + t1117 - t801) * t734) * t732 + t753;
t762 = (t620 + t837) * t1048;
t43 = -t985 / 0.2e1 - t986 / 0.2e1 + (-t1188 * t735 * t858 + t867 * t731) * t732 + (t866 * t732 + (pkin(5) / 0.2e1 + t1035 / 0.2e1 + t412 * t1055 + t445 * t1057) * m(7) - t925) * t734 + t762;
t815 = -t15 * qJD(1) - t43 * qJD(2);
t811 = -t508 * t731 + t968;
t810 = t587 * t665 - t588 * t668;
t809 = t966 - t967;
t808 = -t665 * t731 - t668 * t734;
t118 = (t731 * t647 + t734 * t651 + t1190) * t735;
t746 = (-t412 * t449 + t446 * t448 + t735 * t820) * t1146 + t448 * t1078 + t449 * t1074;
t769 = t176 * t1147 + t995 / 0.2e1 - t994 / 0.2e1;
t789 = t1060 * t317 + t888;
t26 = t735 * t789 + t746 + t769;
t807 = qJD(1) * t26 - qJD(2) * t118;
t798 = -t162 + t819;
t797 = -t722 / 0.2e1 - t721 / 0.2e1;
t796 = pkin(10) * t1077 + t871;
t795 = t110 * t734 - t731 * t97 + t135;
t794 = -t235 + t974 - t975;
t793 = -t130 * t731 + t131 * t734 + t196;
t792 = -t412 * t731 + t641 + t971;
t791 = -t417 * t731 + t450 * t734 + t642;
t781 = t809 * pkin(10);
t777 = (t665 * t734 - t668 * t731) * t1146;
t776 = -t536 * t731 + t537 * t734 + 0.2e1 * t952;
t748 = t1188 * (-t967 / 0.2e1 + t966 / 0.2e1);
t10 = t810 * t1146 + t781 * t1148 + (-t353 / 0.2e1 + mrSges(5,2) * t903 - t989 / 0.2e1 + t875 * t734 + t874 * t731 + t794 * t1149 + t798 * t1147 + t873) * t735 + t748 + (t1122 + t1123 + (t1115 - t383 / 0.2e1) * t734 + (t386 / 0.2e1 + t1114) * t731 + t793 * t1149 + t795 * t1147 + (mrSges(5,1) + t865 - t1043 / 0.2e1 + t1145 / 0.2e1) * t963) * t732;
t270 = 0.4e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * (-0.1e1 + t947) * t957;
t36 = t777 + (-t622 / 0.2e1 - t621 / 0.2e1 - t867 * t734 + t866 * t731 + t811 * t1149 + t792 * t1147) * t735 + ((t1080 + t1081) * t734 + (t650 / 0.2e1 + t649 / 0.2e1) * t731 + t776 * t1149 + t791 * t1147 - t868) * t732;
t775 = -t10 * qJD(1) - t36 * qJD(2) + t270 * qJD(3);
t143 = t449 * t919 - t435;
t514 = -m(7) * t938 - t837;
t625 = -m(7) * t1037 - t1166;
t774 = qJD(1) * t143 + qJD(2) * t514 + qJD(4) * t625;
t766 = t1160 / 0.4e1 + t1159 / 0.4e1 + t665 * t1137 + t859;
t128 = (t731 * t925 + t1065 + t1066 + t842) * t735;
t749 = t681 / 0.4e1 + t683 / 0.4e1 + mrSges(7,3) * t1068 + (t1069 + t1043 / 0.2e1) * pkin(5) - t861;
t739 = t670 * t1046 + t858 * t1131 + (t710 * t1139 + (-Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1) * t731 - t766) * t731 + ((-Ifges(6,2) / 0.4e1 - Ifges(7,2) / 0.4e1) * t734 + (-Ifges(7,4) / 0.2e1 - Ifges(6,4) / 0.2e1) * t731 + t749) * t734;
t750 = mrSges(6,1) * t1098 + mrSges(7,1) * t1108 + pkin(5) * t1075 + t1138 * t450 + t1140 * t537;
t754 = pkin(4) * t1084 + t1062 * t708 + t1066 * t642 + t1070 * t647;
t758 = pkin(10) * t1072 + (t445 / 0.2e1 + t1113) * mrSges(7,3) + t869;
t18 = t651 * t1068 + (t926 + t796) * t731 + (pkin(5) * t1108 + t642 * t1037 / 0.2e1 + t412 * t1068 + t445 * t1067) * m(7) + (-t841 + t863) * t732 + (t758 - t1176) * t734 + (t739 - t920) * t735 + t750 + t754;
t53 = pkin(4) * t670 - t710 * t1166 - t666 * t1037 + (-t1159 / 0.2e1 - t1160 / 0.2e1 + t860) * t734 + (-t681 / 0.2e1 - pkin(5) * t1043 - t683 / 0.2e1 + t862) * t731;
t741 = -t801 * t668 + t863 * t614 + t749 * t449 - pkin(4) * t246 / 0.2e1 + t135 * t1066 + t196 * t1065 + t317 * t1070 + t245 * t1062;
t747 = t250 / 0.4e1 + t249 / 0.4e1 + (t913 - t320 / 0.2e1) * pkin(10) + (t75 / 0.2e1 - t62 / 0.2e1) * mrSges(7,3) + t878;
t752 = t840 * pkin(5) + t251 / 0.4e1 + t252 / 0.4e1 + t880;
t756 = mrSges(7,1) * t1134 + t110 * t1138 + t1140 * t131 + t1142 * t130;
t8 = (t747 - t1181) * t734 + (m(7) * t1134 - t385 / 0.2e1) * pkin(5) + t766 * t448 - t920 * t615 + (-t1177 + (t1002 / 0.2e1 + t1119) * pkin(10) + t752) * t731 + t741 + t756;
t761 = t8 * qJD(1) + t18 * qJD(2) - t128 * qJD(3) - t53 * qJD(4);
t757 = m(7) * (-t448 * t668 - t449 * t665 + t819);
t29 = (-t990 / 0.2e1 + t1120 + t914) * t734 + (-t991 / 0.2e1 + t1118 + t912) * t731 + t929 - t757 / 0.2e1;
t371 = m(7) * t808 + t853;
t516 = (-0.1e1 / 0.2e1 - t858) * t1042;
t751 = m(7) * ((-t665 * t735 + t446) * t734 + (t668 * t735 - t412) * t731);
t77 = (-t983 / 0.2e1 + t1079) * t734 + (-t984 / 0.2e1 + t1073) * t731 - t928 - t751 / 0.2e1;
t759 = -qJD(1) * t29 - qJD(2) * t77 + qJD(3) * t516 + qJD(4) * t371;
t515 = (t947 + 0.1e1) * t927;
t164 = m(7) * t844 + t1161 * t927;
t129 = t897 * t1144 + (-t731 * t803 + t842) * t735 + (t670 + t1166) * t1049;
t78 = t647 * t1051 + t751 / 0.2e1 + t651 * t1060 - t928 + t797 * t732;
t44 = m(7) * (t732 * t951 - t1035) * t1051 + t924 * t731 + t803 * t734 - t762 + t1168 * t900 + t1167 * t893 - t1188 * t947 * t957 / 0.2e1;
t37 = (t732 * t776 + t735 * t811) * t1148 + (t732 * t791 + t735 * t792) * t1146 + t777 + (t650 + t649) * t900 + t1167 * t897 + (t646 + t645) * t891 + t1168 * t885 + t1169 * t1055 + (-t622 - t621) * t1049;
t30 = t757 / 0.2e1 + t319 * t1060 + t317 * t1051 + t929 + t797 * t614 + t1161 * t1136;
t25 = t317 * t897 + t319 * t887 + t746 - t769;
t17 = t758 * t734 - t750 + t417 * t939 + t739 * t735 + t1009 / 0.2e1 + t1007 / 0.2e1 + t863 * t732 - (t1074 + t857 / 0.2e1) * t668 + (t796 + t1189) * t731 + t754 + t1184 * t899 + t1185 * t893;
t16 = (-pkin(5) * t449 * t735 + t1029 * t958) * t1146 + t732 * t888 + t320 * t893 + t753 + t1173 * t900 + (t246 + t245) * t1049 + t1188 * (t448 * t899 + t449 * t893);
t12 = t730 * t833 / 0.2e1 + t320 * t1051 + t318 * t1059 + t755 * t732 + t836 * t735 + t930 + t744 - t789 + t1169 * t844 - t1193;
t11 = (t981 / 0.2e1 + (t1143 + t865) * t732) * t963 + t748 + t247 * t1055 + (t386 + t385) * t900 + t1172 * t897 + (t384 + t383) * t891 + t1173 * t885 + (pkin(4) * t917 + t732 * t793 + t735 * t794 + t781) * t1148 + (-t710 * t917 + t732 * t795 + t735 * t798 + t810) * t1146 + (-t354 - t353) * t1049 + t1153;
t7 = t747 * t734 + (t1059 * t1131 + t766) * t448 + t1010 / 0.2e1 + t1008 / 0.2e1 + t97 * t939 + pkin(5) * t1114 + (pkin(10) * t1119 + t752) * t731 + t741 - t756 + t1184 * t965 / 0.2e1 - t1185 * t964 / 0.2e1;
t5 = t745 + t446 * t912 + t509 * t913 + t412 * t914 + (t1029 * t446 + t951 * t76) * t1146 + t837 * t1126 + t996 * t1104 - t735 * t804 - t743 + t246 * t882 + t1157 * t960 / 0.2e1 + t1158 * t887 - t1175 * t953 / 0.4e1 + t1195 * t884 + t869 * t448 + (t926 + t871) * t449 + (t247 * t885 + (t135 * t953 + t449 * t642) * t1146 + m(7) * t1130 + t1110) * pkin(5) + t1154 * t896 + (-t788 - t787 - t786 - t785 - t784 - t783) * t1047;
t1 = t771 / 0.4e1 + t770 / 0.4e1 - t765 / 0.4e1 + (t266 + t265) * t896 - (t830 + t679) * t615 / 0.4e1 - ((t196 * t732 - t972) * t1148 + t1153) * t1133 - t738 + (t268 + t267) * t884 - t773 / 0.4e1 + t961 * t1196 + t1039 * t1149 + mrSges(5,2) * t1121 + qJ(3) * (mrSges(5,1) * t615 - mrSges(5,2) * t614) / 0.2e1 + (-Ifges(5,2) * t615 - t600) * t1056 + (-Ifges(5,1) * t614 - t1023) * t1047 + t1180 / 0.2e1 + (t574 + t573) * t615 / 0.4e1 + (t581 + t579) * t449 / 0.4e1 + (t577 + t575) * t448 / 0.4e1 + ((m(6) * t1128 + t1111) * pkin(10) + t768) * t734 + ((m(6) * t1129 + t1109) * pkin(10) + t767) * t731 + (t1011 / 0.2e1 + (t1018 / 0.2e1 + (Ifges(5,6) / 0.2e1 + t864) * t735) * t733) * t729 - t824 * t963 / 0.4e1 + t760 * t1148 + t742 + (t1165 * t614 + t1008 + t1010) * t732 / 0.4e1 + (t1165 * t732 + t1007 + t1009 + t827) * t1086 - t1170 * t964 / 0.4e1 + t1171 * t965 / 0.4e1 - t1174 * t958 / 0.4e1 - t354 * t882;
t20 = [-qJD(2) * t3 - qJD(3) * t14 + qJD(4) * t6 + qJD(5) * t9 + qJD(6) * t35, t12 * qJD(3) + t1 * qJD(4) + t5 * qJD(5) + t25 * qJD(6) - t1003 + ((t518 * mrSges(5,1) - t304 * mrSges(5,3) - t1133 * t591 + t764) * t732 + (t518 * mrSges(5,2) - t303 * mrSges(5,3) - (t590 - t326) * t1133 + (t274 / 0.2e1 + t275 / 0.2e1) * t734 + (-t272 / 0.2e1 - t273 / 0.2e1) * t731) * t735 + ((-mrSges(4,1) * pkin(2) + Ifges(5,5) * t735 - Ifges(5,6) * t732 + t1187) * t736 + (-qJ(3) * mrSges(4,1) + t780 * t732 + (t732 * t800 + t763) * t735 + t1186) * t733) * t729 + 0.2e1 * (t508 * t132 + t509 * t133 + t263 * t952) * t1148 + 0.2e1 * (qJ(3) * t518 - t1133 * t813) * t1150 + 0.2e1 * (-pkin(2) * t618 - qJ(3) * t617) * t1151 + 0.2e1 * (t105 * t412 + t109 * t446 + t176 * t642) * t1146 + t105 * t651 + t132 * t652 + t642 * t325 + t109 * t647 + t133 * t648 + t1182 + t618 * mrSges(4,2) - t618 * mrSges(3,1) + t176 * t623 + t263 * t624 - t617 * mrSges(4,3) + t509 * t414 + t508 * t416 + t446 * t413 + t412 * t415 + t1155 * t539 + t1170 * t540 / 0.2e1) * qJD(2), -t973 + t12 * qJD(2) + 0.2e1 * (t1146 + t1148) * (t809 + t916) * t732 * qJD(3) + t11 * qJD(4) + t16 * qJD(5) + t164 * qJD(6), t993 + t1 * qJD(2) + t11 * qJD(3) + t7 * qJD(5) + t30 * qJD(6) + (t131 * mrSges(6,3) + t110 * mrSges(7,3) + t860 * t614 + (m(6) * t131 + t384) * pkin(10) + t877) * t943 + (-t130 * mrSges(6,3) - t97 * mrSges(7,3) + t862 * t614 + (-m(6) * t130 - t386) * pkin(10) + t876) * t945 + (-t710 * t353 + t665 * t385 + t162 * t666 - t668 * t383 + pkin(4) * t354 - t234 * mrSges(5,2) + m(7) * (-t110 * t668 + t162 * t710 + t665 * t97) + t948 + t852 * t235 + t1156 * t615) * qJD(4), t979 + t5 * qJD(2) + t16 * qJD(3) + t7 * qJD(4) + (-mrSges(6,1) * t104 - mrSges(7,1) * t76 - mrSges(6,2) * t103 - mrSges(7,2) * t75 + (-t1001 - t1192) * pkin(5) + t1162) * qJD(5), qJD(2) * t25 + qJD(3) * t164 + qJD(4) * t30 + t976; -qJD(3) * t13 - qJD(4) * t2 - qJD(5) * t4 + qJD(6) * t26 + t1003, qJD(3) * t55 + qJD(4) * t19 - qJD(5) * t24 - qJD(6) * t118, qJD(4) * t37 + qJD(5) * t44 - t816, t37 * qJD(3) + (-t710 * t621 + t665 * t649 - t641 * t666 - t668 * t645 + pkin(4) * t622 + m(7) * (t417 * t665 - t450 * t668 - t641 * t710)) * qJD(4) + t17 * qJD(5) + t78 * qJD(6) + (mrSges(5,2) * t1133 - Ifges(5,6) + t1156) * t942 + (-t1133 * t852 - Ifges(5,5)) * t944 + (t537 * mrSges(6,3) + t450 * mrSges(7,3) + t860 * t732 + (m(6) * t537 + t646) * pkin(10) - t872) * t943 + (-t536 * mrSges(6,3) - t417 * mrSges(7,3) + t862 * t732 + (-m(6) * t536 - t650) * pkin(10) + t870) * t945 + t818, t44 * qJD(3) + t17 * qJD(4) + t817 + (-mrSges(6,1) * t509 - mrSges(6,2) * t508 - mrSges(7,2) * t445 + (-t1184 * t734 + (mrSges(7,3) * pkin(5) - t1185) * t731) * t735 + t919 * t446) * qJD(5), qJD(4) * t78 + t807; qJD(2) * t13 - qJD(4) * t10 - qJD(5) * t15 - qJD(6) * t163 + t973, -qJD(4) * t36 - qJD(5) * t43 + t816, t270 * qJD(4), t129 * qJD(5) + t515 * qJD(6) + (t666 + t978) * t944 + (mrSges(6,3) * t947 - mrSges(5,2) + t853) * t942 + 0.2e1 * ((t1034 * t947 - t1038) * t1148 + (t710 * t732 + t735 * t808) * t1146) * qJD(4) + t775, t129 * qJD(4) + ((mrSges(6,2) + mrSges(7,2)) * t731 + (-mrSges(6,1) + t919) * t734) * qJD(5) * t732 + t815, qJD(4) * t515 - t946; qJD(2) * t2 + qJD(3) * t10 + qJD(5) * t8 - qJD(6) * t29 - t993, qJD(3) * t36 + qJD(5) * t18 - qJD(6) * t77 - t818, -qJD(5) * t128 + qJD(6) * t516 - t775, -qJD(5) * t53 + qJD(6) * t371, t761 + (-mrSges(7,2) * t665 - mrSges(7,3) * t1036 - pkin(10) * t982 + t723 + t724 + (mrSges(6,2) * pkin(10) - t1184) * t731 - t919 * t668) * qJD(5), t759; qJD(2) * t4 + qJD(3) * t15 - qJD(4) * t8 + qJD(6) * t143 - t979, qJD(3) * t43 - qJD(4) * t18 + qJD(6) * t514 - t817, qJD(4) * t128 - t815, qJD(6) * t625 - t761, 0, t774; -qJD(2) * t26 + qJD(3) * t163 + qJD(4) * t29 - qJD(5) * t143 - t976, qJD(4) * t77 - qJD(5) * t514 - t807, -qJD(4) * t516 + t946, -qJD(5) * t625 - t759, -t774, 0;];
Cq  = t20;
