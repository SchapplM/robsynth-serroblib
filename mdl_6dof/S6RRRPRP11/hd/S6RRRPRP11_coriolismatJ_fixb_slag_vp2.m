% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:42
% EndTime: 2019-03-09 17:42:37
% DurationCPUTime: 30.08s
% Computational Cost: add. (29651->1257), mult. (70333->1659), div. (0->0), fcn. (71596->8), ass. (0->578)
t738 = sin(qJ(5));
t735 = t738 ^ 2;
t741 = cos(qJ(5));
t736 = t741 ^ 2;
t918 = t735 + t736;
t1135 = Ifges(6,4) + Ifges(7,4);
t1134 = Ifges(5,4) - Ifges(4,5);
t1124 = Ifges(5,5) - Ifges(4,6);
t1123 = Ifges(7,5) + Ifges(6,5);
t1122 = Ifges(7,6) + Ifges(6,6);
t742 = cos(qJ(3));
t1133 = t918 * t742 / 0.2e1;
t995 = Ifges(7,4) * t738;
t681 = t741 * Ifges(7,1) - t995;
t997 = Ifges(6,4) * t738;
t682 = t741 * Ifges(6,1) - t997;
t1131 = -t681 / 0.4e1 - t682 / 0.4e1;
t1132 = pkin(4) + pkin(9);
t737 = sin(pkin(6));
t739 = sin(qJ(3));
t740 = sin(qJ(2));
t743 = cos(qJ(2));
t924 = t741 * t743;
t533 = (-t738 * t740 + t739 * t924) * t737;
t929 = t739 * t743;
t534 = (t738 * t929 + t740 * t741) * t737;
t1130 = t1135 * t534 + (Ifges(6,2) + Ifges(7,2)) * t533;
t1129 = (Ifges(6,1) + Ifges(7,1)) * t534 + t1135 * t533;
t728 = t741 * mrSges(7,1);
t1128 = -t728 / 0.2e1;
t934 = t737 * t743;
t875 = -t934 / 0.2e1;
t1127 = -mrSges(4,1) + mrSges(5,2);
t1126 = mrSges(4,2) - mrSges(5,3);
t1121 = Ifges(4,3) + Ifges(5,1);
t1120 = Ifges(6,3) + Ifges(7,3);
t935 = t737 * t740;
t953 = cos(pkin(6));
t605 = t739 * t935 - t742 * t953;
t885 = pkin(1) * t953;
t607 = -pkin(8) * t935 + t743 * t885;
t542 = -pkin(2) * t953 - t607;
t606 = t739 * t953 + t742 * t935;
t767 = -t606 * qJ(4) + t542;
t281 = t605 * pkin(3) + t767;
t669 = -mrSges(5,2) * t739 - mrSges(5,3) * t742;
t1119 = t281 * t669;
t670 = mrSges(4,1) * t739 + mrSges(4,2) * t742;
t1118 = t542 * t670;
t902 = Ifges(6,4) / 0.2e1 + Ifges(7,4) / 0.2e1;
t1117 = t902 * t741;
t1088 = m(7) * pkin(5);
t893 = mrSges(7,1) + t1088;
t1078 = pkin(3) + pkin(10);
t1116 = -qJ(6) - t1078;
t474 = t605 * t741 + t738 * t934;
t475 = -t605 * t738 + t737 * t924;
t971 = t475 * Ifges(7,4);
t179 = t474 * Ifges(7,2) + t606 * Ifges(7,6) - t971;
t972 = t475 * Ifges(6,4);
t180 = t474 * Ifges(6,2) + t606 * Ifges(6,6) - t972;
t1115 = t180 + t179;
t469 = Ifges(7,4) * t474;
t181 = -Ifges(7,1) * t475 + Ifges(7,5) * t606 + t469;
t470 = Ifges(6,4) * t474;
t182 = -Ifges(6,1) * t475 + Ifges(6,5) * t606 + t470;
t1114 = t182 + t181;
t967 = t606 * mrSges(7,1);
t973 = t475 * mrSges(7,3);
t319 = t967 + t973;
t320 = mrSges(6,1) * t606 + t475 * mrSges(6,3);
t1113 = t319 + t320;
t931 = t738 * t742;
t959 = t739 * mrSges(7,1);
t646 = mrSges(7,3) * t931 + t959;
t647 = mrSges(6,1) * t739 + mrSges(6,3) * t931;
t1112 = t646 + t647;
t925 = t741 * t742;
t911 = mrSges(7,3) * t925;
t958 = t739 * mrSges(7,2);
t650 = -t911 - t958;
t651 = -t739 * mrSges(6,2) - mrSges(6,3) * t925;
t1111 = t651 + t650;
t994 = Ifges(7,4) * t741;
t677 = -t738 * Ifges(7,2) + t994;
t996 = Ifges(6,4) * t741;
t678 = -t738 * Ifges(6,2) + t996;
t1110 = t677 + t678;
t962 = t738 * mrSges(7,2);
t1108 = t728 - t962;
t729 = Ifges(5,5) * t739;
t730 = Ifges(4,5) * t742;
t1107 = t729 + t730;
t731 = Ifges(4,4) * t742;
t680 = -Ifges(4,2) * t739 + t731;
t923 = t742 * t743;
t890 = t737 * t923;
t1106 = t1122 * t890 + t1130;
t1105 = t1123 * t890 + t1129;
t1020 = t741 / 0.2e1;
t957 = t741 * mrSges(7,2);
t666 = t738 * mrSges(7,1) + t957;
t1032 = -t666 / 0.2e1;
t248 = -t474 * mrSges(7,1) - t475 * mrSges(7,2);
t1104 = t248 * t1020 + t475 * t1032;
t1022 = -t741 / 0.2e1;
t1023 = t738 / 0.2e1;
t1103 = t474 * t1022 + t475 * t1023;
t1102 = t1124 * t606 + t1134 * t605;
t608 = (pkin(2) * t740 - pkin(9) * t743) * t737;
t406 = -t739 * t607 + t608 * t742;
t256 = (pkin(4) * t923 - t1078 * t740) * t737 - t406;
t950 = qJ(4) * t742;
t811 = pkin(10) * t739 - t950;
t609 = pkin(8) * t934 + t740 * t885;
t891 = t737 * t929;
t887 = pkin(3) * t891 + t609;
t347 = t811 * t934 + t887;
t123 = t741 * t256 - t347 * t738;
t124 = t738 * t256 + t741 * t347;
t803 = t123 * t741 + t124 * t738;
t843 = t677 / 0.4e1 + t678 / 0.4e1;
t819 = t738 * Ifges(7,1) + t994;
t784 = t819 * t742;
t563 = t739 * Ifges(7,5) - t784;
t821 = t738 * Ifges(6,1) + t996;
t785 = t821 * t742;
t565 = t739 * Ifges(6,5) - t785;
t848 = t565 / 0.2e1 + t563 / 0.2e1;
t815 = t741 * Ifges(7,2) + t995;
t782 = t815 * t742;
t559 = t739 * Ifges(7,6) - t782;
t817 = t741 * Ifges(6,2) + t997;
t783 = t817 * t742;
t561 = t739 * Ifges(6,6) - t783;
t851 = t559 / 0.2e1 + t561 / 0.2e1;
t1025 = -t738 / 0.2e1;
t1101 = mrSges(6,3) * t1133 + t647 * t1025;
t955 = t741 * Ifges(7,6);
t960 = t738 * Ifges(7,5);
t813 = t955 + t960;
t956 = t741 * Ifges(6,6);
t961 = t738 * Ifges(6,5);
t814 = t956 + t961;
t1100 = -t813 / 0.4e1 - t814 / 0.4e1;
t903 = Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1;
t1099 = -t903 * t738 - t843;
t1098 = t609 * mrSges(3,1) + t607 * mrSges(3,2);
t318 = -t606 * mrSges(6,2) + t474 * mrSges(6,3);
t927 = t741 * t318;
t1097 = t1103 * mrSges(6,3) + t320 * t1025 + t927 / 0.2e1;
t1095 = 2 * qJD(3);
t1094 = -m(5) / 0.2e1;
t1093 = m(5) / 0.2e1;
t1092 = -m(6) / 0.2e1;
t1091 = m(6) / 0.2e1;
t1090 = -m(7) / 0.2e1;
t1089 = m(7) / 0.2e1;
t1087 = -mrSges(6,1) / 0.2e1;
t1086 = -mrSges(6,2) / 0.2e1;
t1085 = mrSges(6,2) / 0.2e1;
t1084 = -mrSges(7,2) / 0.2e1;
t1083 = mrSges(7,2) / 0.2e1;
t1082 = -mrSges(6,3) / 0.2e1;
t1081 = -mrSges(7,3) / 0.2e1;
t1080 = mrSges(7,3) / 0.2e1;
t1010 = pkin(4) * t605;
t543 = pkin(9) * t953 + t609;
t544 = (-pkin(2) * t743 - pkin(9) * t740 - pkin(1)) * t737;
t335 = t742 * t543 + t739 * t544;
t237 = t335 - t1010;
t217 = t741 * t237;
t952 = qJ(4) * t605;
t297 = t1078 * t606 + t952;
t92 = -pkin(5) * t605 + t217 + (-qJ(6) * t606 - t297) * t738;
t1079 = -t92 / 0.2e1;
t1077 = pkin(3) * mrSges(5,1);
t1076 = -qJ(4) / 0.2e1;
t101 = pkin(5) * t890 - qJ(6) * t534 + t123;
t1075 = -t101 / 0.2e1;
t108 = qJ(6) * t533 + t124;
t1074 = t108 / 0.2e1;
t1073 = t124 / 0.2e1;
t1009 = pkin(5) * t474;
t889 = qJ(4) * t934;
t285 = t889 - t335;
t204 = -t285 - t1010;
t131 = t204 - t1009;
t1072 = t131 / 0.2e1;
t1071 = t204 / 0.2e1;
t966 = t606 * mrSges(7,2);
t974 = t474 * mrSges(7,3);
t317 = -t966 + t974;
t1070 = t317 / 0.2e1;
t1069 = t318 / 0.2e1;
t1068 = -t319 / 0.2e1;
t1067 = -t320 / 0.2e1;
t732 = t739 * pkin(3);
t638 = t732 + t811;
t686 = t1132 * t742;
t642 = t741 * t686;
t368 = pkin(5) * t742 + t642 + (-qJ(6) * t739 - t638) * t738;
t1066 = t368 / 0.2e1;
t937 = t606 * t738;
t396 = -mrSges(7,1) * t605 - mrSges(7,3) * t937;
t1065 = -t396 / 0.2e1;
t1064 = t396 / 0.2e1;
t951 = qJ(4) * t739;
t623 = -t1078 * t742 - pkin(2) - t951;
t685 = t1132 * t739;
t462 = -t623 * t738 + t741 * t685;
t404 = qJ(6) * t931 + t462;
t1063 = t404 / 0.2e1;
t421 = -mrSges(5,2) * t605 - mrSges(5,3) * t606;
t1062 = t421 / 0.2e1;
t1061 = t462 / 0.2e1;
t471 = -t638 * t738 + t642;
t1060 = t471 / 0.2e1;
t1059 = -t474 / 0.2e1;
t1058 = t475 / 0.2e1;
t1008 = pkin(5) * t741;
t886 = -pkin(4) - t1008;
t603 = (-pkin(9) + t886) * t739;
t1048 = t603 / 0.2e1;
t604 = pkin(5) * t925 + t686;
t1047 = t604 / 0.2e1;
t613 = t1108 * t742;
t1045 = t613 / 0.2e1;
t825 = mrSges(6,1) * t741 - mrSges(6,2) * t738;
t614 = t825 * t742;
t1044 = t614 / 0.2e1;
t932 = t738 * t739;
t644 = mrSges(7,1) * t742 - mrSges(7,3) * t932;
t1043 = t644 / 0.2e1;
t1042 = -t646 / 0.2e1;
t1041 = t646 / 0.2e1;
t1040 = -t647 / 0.2e1;
t1039 = t647 / 0.2e1;
t1038 = t650 / 0.2e1;
t1037 = t651 / 0.2e1;
t659 = t1116 * t738;
t1036 = -t659 / 0.2e1;
t660 = t1116 * t741;
t1035 = -t660 / 0.2e1;
t1034 = t660 / 0.2e1;
t665 = t742 * mrSges(5,2) - t739 * mrSges(5,3);
t1033 = t665 / 0.2e1;
t667 = t738 * mrSges(6,1) + t741 * mrSges(6,2);
t1031 = -t667 / 0.2e1;
t1028 = -t685 / 0.2e1;
t1027 = t686 / 0.2e1;
t723 = pkin(5) * t738 + qJ(4);
t1026 = -t723 / 0.2e1;
t1024 = -t738 / 0.4e1;
t1021 = -t741 / 0.4e1;
t1019 = t741 / 0.4e1;
t1016 = t1078 / 0.2e1;
t1015 = -t1078 / 0.2e1;
t1014 = m(7) * t131;
t334 = t543 * t739 - t742 * t544;
t165 = t606 * t886 - t334;
t1013 = m(7) * t165;
t1012 = m(7) * t604;
t1011 = m(7) * t723;
t236 = -pkin(4) * t606 - t334;
t715 = pkin(3) * t934;
t185 = pkin(10) * t934 - t236 + t715;
t198 = t1078 * t605 + t767;
t94 = t741 * t185 - t198 * t738;
t74 = qJ(6) * t475 + t94;
t61 = pkin(5) * t606 + t74;
t1007 = t61 * mrSges(7,3);
t95 = t185 * t738 + t198 * t741;
t75 = qJ(6) * t474 + t95;
t1006 = t75 * mrSges(7,3);
t1005 = t94 * mrSges(6,3);
t1004 = t95 * mrSges(6,3);
t1003 = mrSges(7,2) + mrSges(6,2);
t1002 = -t61 + t74;
t1001 = mrSges(5,2) * t740;
t1000 = Ifges(3,4) * t740;
t999 = Ifges(3,4) * t743;
t998 = Ifges(4,4) * t739;
t993 = Ifges(6,5) * t534;
t992 = Ifges(7,5) * t534;
t990 = Ifges(5,6) * t739;
t989 = Ifges(5,6) * t742;
t988 = Ifges(6,6) * t533;
t987 = Ifges(7,6) * t533;
t986 = Ifges(6,3) * t605;
t985 = Ifges(7,3) * t605;
t177 = -t475 * Ifges(7,5) + t474 * Ifges(7,6) + t606 * Ifges(7,3);
t178 = -t475 * Ifges(6,5) + t474 * Ifges(6,6) + t606 * Ifges(6,3);
t407 = t742 * t607 + t739 * t608;
t346 = -qJ(4) * t935 - t407;
t282 = -pkin(4) * t891 - t346;
t188 = -pkin(5) * t533 + t282;
t249 = -t474 * mrSges(6,1) - t475 * mrSges(6,2);
t286 = t334 + t715;
t348 = -pkin(3) * t935 - t406;
t964 = t606 * Ifges(5,6);
t353 = -Ifges(5,5) * t934 + t605 * Ifges(5,3) - t964;
t573 = Ifges(5,6) * t605;
t354 = -Ifges(5,4) * t934 - t606 * Ifges(5,2) + t573;
t965 = t606 * Ifges(4,4);
t355 = -t605 * Ifges(4,2) - Ifges(4,6) * t934 + t965;
t578 = Ifges(4,4) * t605;
t356 = t606 * Ifges(4,1) - Ifges(4,5) * t934 - t578;
t430 = -t742 * t889 + t887;
t437 = -mrSges(7,2) * t890 + mrSges(7,3) * t533;
t438 = -mrSges(6,2) * t890 + mrSges(6,3) * t533;
t439 = mrSges(7,1) * t890 - mrSges(7,3) * t534;
t440 = mrSges(6,1) * t890 - mrSges(6,3) * t534;
t712 = mrSges(5,3) * t934;
t968 = t605 * mrSges(5,1);
t492 = t712 + t968;
t493 = t606 * mrSges(5,1) - mrSges(5,2) * t934;
t494 = mrSges(4,2) * t934 - t605 * mrSges(4,3);
t495 = -mrSges(4,1) * t934 - t606 * mrSges(4,3);
t672 = Ifges(5,3) * t739 - t989;
t674 = -Ifges(5,2) * t742 + t990;
t684 = Ifges(4,1) * t742 - t998;
t780 = (mrSges(5,1) * t929 - mrSges(5,3) * t740) * t737;
t781 = (mrSges(5,1) * t923 + t1001) * t737;
t969 = t534 * mrSges(7,2);
t970 = t533 * mrSges(7,1);
t824 = t969 - t970;
t826 = -mrSges(6,1) * t533 + mrSges(6,2) * t534;
t874 = t934 / 0.2e1;
t830 = t742 * t874;
t831 = t742 * t875;
t832 = t739 * t874;
t833 = t739 * t875;
t876 = t935 / 0.2e1;
t938 = t606 * t737;
t939 = t605 * t737;
t3 = -(-t1121 * t934 + t1124 * t605 - t1134 * t606) * t935 / 0.2e1 + (-t743 * (-Ifges(3,2) * t740 + t999) / 0.2e1 - t740 * (Ifges(3,1) * t743 - t1000) / 0.2e1 + pkin(1) * (mrSges(3,1) * t740 + mrSges(3,2) * t743) + ((-Ifges(5,4) * t742 - Ifges(4,6) * t739 + t1107) * t743 + t1121 * t740) * t743 / 0.2e1) * t737 ^ 2 + ((Ifges(3,2) * t743 + t1000) * t876 + (Ifges(3,1) * t740 + t999) * t875 + t334 * (mrSges(4,1) * t740 - mrSges(4,3) * t923) - t335 * (-mrSges(4,2) * t740 - mrSges(4,3) * t929)) * t737 + (t356 + t178 + t177) * t831 - m(4) * (-t334 * t406 + t335 * t407 + t542 * t609) - m(5) * (t281 * t430 + t285 * t346 + t286 * t348) - m(6) * (t123 * t94 + t124 * t95 + t204 * t282) - m(7) * (t101 * t61 + t108 * t75 + t131 * t188) - t285 * t780 - t286 * t781 - (t1120 * t890 + t987 + t988 + t992 + t993) * t606 / 0.2e1 + (-t1118 - t1119) * t934 - t1114 * t534 / 0.2e1 - t1115 * t533 / 0.2e1 + t1105 * t1058 + t1106 * t1059 + (Ifges(4,6) * t740 + t680 * t743) * t939 / 0.2e1 + (Ifges(3,5) * t875 + Ifges(3,6) * t876 - (Ifges(3,5) * t743 - Ifges(3,6) * t740) * t737 / 0.2e1 + t1098) * t953 - t609 * (mrSges(4,1) * t605 + mrSges(4,2) * t606) - t346 * t492 - t348 * t493 - t407 * t494 - t406 * t495 - t430 * t421 - t75 * t437 - t95 * t438 - t61 * t439 - t94 * t440 - t124 * t318 - t101 * t319 - t123 * t320 - t108 * t317 - t282 * t249 - t188 * t248 - (Ifges(4,5) * t740 + t684 * t743) * t938 / 0.2e1 + (Ifges(5,4) * t740 + t674 * t743) * t938 / 0.2e1 - (Ifges(5,5) * t740 + t672 * t743) * t939 / 0.2e1 - t204 * t826 - t131 * t824 + t354 * t830 + t355 * t832 + t353 * t833;
t984 = t3 * qJD(1);
t366 = pkin(5) * t739 + t404;
t983 = t366 * mrSges(7,3);
t982 = t366 * t75;
t114 = t738 * t237 + t741 * t297;
t936 = t606 * t741;
t102 = qJ(6) * t936 + t114;
t113 = -t297 * t738 + t217;
t265 = t813 * t606 - t985;
t266 = t814 * t606 - t986;
t369 = t1108 * t606;
t370 = t825 * t606;
t397 = -mrSges(6,1) * t605 - mrSges(6,3) * t937;
t398 = mrSges(7,2) * t605 + mrSges(7,3) * t936;
t399 = mrSges(6,2) * t605 + mrSges(6,3) * t936;
t418 = -mrSges(5,2) * t606 + mrSges(5,3) * t605;
t419 = pkin(3) * t606 + t952;
t420 = mrSges(4,1) * t606 - mrSges(4,2) * t605;
t422 = Ifges(5,3) * t606 + t573;
t423 = Ifges(5,2) * t605 + t964;
t424 = -Ifges(4,2) * t606 - t578;
t425 = -Ifges(4,1) * t605 - t965;
t269 = -Ifges(7,5) * t605 + t606 * t819;
t270 = -Ifges(6,5) * t605 + t606 * t821;
t857 = -t270 / 0.2e1 - t269 / 0.2e1;
t267 = -Ifges(7,6) * t605 + t606 * t815;
t268 = -Ifges(6,6) * t605 + t606 * t817;
t858 = t267 / 0.2e1 + t268 / 0.2e1;
t860 = t181 / 0.2e1 + t182 / 0.2e1;
t862 = t179 / 0.2e1 + t180 / 0.2e1;
t4 = t1102 * t875 + t857 * t475 + t858 * t474 + (-t423 / 0.2e1 - t355 / 0.2e1 + t425 / 0.2e1 + t353 / 0.2e1 + t265 / 0.2e1 + t266 / 0.2e1 - t335 * mrSges(4,3) + t285 * mrSges(5,1) + t862 * t741 + t860 * t738) * t606 + (t492 - t494) * t334 + t542 * t420 + t281 * t418 + t419 * t421 + (-t424 / 0.2e1 - t356 / 0.2e1 - t177 / 0.2e1 - t178 / 0.2e1 + t422 / 0.2e1 + t354 / 0.2e1 - t286 * mrSges(5,1) - t334 * mrSges(4,3)) * t605 + t61 * t396 + t94 * t397 + t75 * t398 + t95 * t399 - t131 * t369 - t204 * t370 + t114 * t318 + t92 * t319 + t113 * t320 + t102 * t317 + t165 * t248 + t236 * t249 + m(5) * (t281 * t419 + t285 * t334 + t286 * t335) + m(6) * (t113 * t94 + t114 * t95 + t204 * t236) + m(7) * (t102 * t75 + t131 * t165 + t61 * t92) + (t493 - t495) * t335;
t981 = t4 * qJD(1);
t980 = t404 * t75;
t463 = t623 * t741 + t685 * t738;
t405 = -qJ(6) * t925 + t463;
t979 = t405 * mrSges(7,3);
t978 = t405 * t61;
t977 = t405 * t74;
t976 = t462 * mrSges(6,3);
t975 = t463 * mrSges(6,3);
t963 = t659 * mrSges(7,3);
t464 = t474 * mrSges(7,2);
t246 = -mrSges(7,1) * t475 + t464;
t247 = -mrSges(6,1) * t475 + mrSges(6,2) * t474;
t250 = Ifges(7,5) * t474 + Ifges(7,6) * t475;
t251 = Ifges(6,5) * t474 + Ifges(6,6) * t475;
t252 = Ifges(7,2) * t475 + t469;
t253 = Ifges(6,2) * t475 + t470;
t254 = Ifges(7,1) * t474 + t971;
t255 = Ifges(6,1) * t474 + t972;
t892 = m(7) * t1002;
t9 = t131 * t246 + t204 * t247 + t74 * t317 + t94 * t318 - t95 * t320 + (t251 / 0.2e1 + t250 / 0.2e1) * t606 + (t892 - t319) * t75 + (t1006 + t1004 - t255 / 0.2e1 - t254 / 0.2e1 + (-t248 - t1014) * pkin(5) + t862) * t475 + (-t1007 - t1005 + t252 / 0.2e1 + t253 / 0.2e1 + t860) * t474;
t954 = t9 * qJD(1);
t809 = t738 * t94 - t741 * t95;
t810 = t61 * t738 - t741 * t75;
t12 = (-t248 - t249 + t492) * t934 + (-t421 + (-t317 - t318) * t741 + t1113 * t738) * t606 + m(7) * (-t131 * t934 + t606 * t810) + m(6) * (-t204 * t934 + t606 * t809) + m(5) * (-t281 * t606 + t285 * t934);
t949 = qJD(1) * t12;
t33 = m(7) * (t474 * t75 + t475 * t61) + t474 * t317 + t475 * t319;
t948 = qJD(1) * t33;
t947 = t101 * t741;
t944 = t131 * t741;
t943 = t366 * t738;
t942 = t474 * t738;
t941 = t475 * t723;
t940 = t475 * t741;
t933 = t738 * t646;
t930 = t739 * t741;
t928 = t741 * t317;
t926 = t741 * t650;
t922 = t1078 * t438;
t921 = -t366 + t404;
t472 = t741 * t638 + t738 * t686;
t917 = qJD(3) * t738;
t916 = qJD(3) * t741;
t915 = qJD(5) * t738;
t914 = qJD(5) * t741;
t912 = t1088 / 0.2e1;
t908 = m(7) * t1048;
t907 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t906 = t1083 + t1085;
t905 = -mrSges(4,3) / 0.2e1 - mrSges(5,1) / 0.2e1;
t904 = t1081 + t1082;
t901 = Ifges(6,4) / 0.4e1 + Ifges(7,4) / 0.4e1;
t900 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t899 = 0.3e1 / 0.4e1 * Ifges(7,5) + 0.3e1 / 0.4e1 * Ifges(6,5);
t898 = Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.4e1;
t897 = -Ifges(5,6) / 0.2e1 - Ifges(4,4) / 0.2e1;
t896 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t895 = 0.3e1 / 0.4e1 * Ifges(6,6) + 0.3e1 / 0.4e1 * Ifges(7,6);
t894 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t888 = mrSges(7,3) * t1025;
t882 = t937 / 0.2e1;
t881 = -t936 / 0.2e1;
t873 = t934 / 0.4e1;
t869 = t931 / 0.2e1;
t866 = -t928 / 0.2e1;
t864 = t741 * t1015;
t861 = t180 / 0.4e1 + t179 / 0.4e1;
t859 = t182 / 0.4e1 + t181 / 0.4e1;
t856 = t1070 + t1069;
t855 = t335 / 0.2e1 + t285 / 0.2e1;
t854 = t437 / 0.2e1 + t438 / 0.2e1;
t853 = t440 / 0.2e1 + t439 / 0.2e1;
t558 = Ifges(7,6) * t742 + t739 * t815;
t560 = Ifges(6,6) * t742 + t739 * t817;
t852 = -t558 / 0.2e1 - t560 / 0.2e1;
t850 = t561 / 0.4e1 + t559 / 0.4e1;
t562 = Ifges(7,5) * t742 + t739 * t819;
t564 = Ifges(6,5) * t742 + t739 * t821;
t849 = -t562 / 0.2e1 - t564 / 0.2e1;
t847 = t565 / 0.4e1 + t563 / 0.4e1;
t846 = t1038 + t1037;
t671 = -Ifges(5,3) * t742 - t990;
t679 = Ifges(4,2) * t742 + t998;
t845 = t671 / 0.2e1 - t679 / 0.2e1;
t844 = t677 / 0.2e1 + t678 / 0.2e1;
t842 = t681 / 0.2e1 + t682 / 0.2e1;
t841 = -t736 / 0.2e1 - t735 / 0.2e1;
t840 = -m(6) * t1078 - mrSges(6,3);
t839 = m(7) * t921;
t837 = mrSges(7,3) * pkin(5) - t1123;
t834 = t892 / 0.2e1;
t828 = t841 * mrSges(6,3);
t827 = (-t730 / 0.4e1 - t729 / 0.4e1) * t743;
t812 = -pkin(3) * t742 - t951;
t675 = t741 * Ifges(7,5) - t738 * Ifges(7,6);
t676 = t741 * Ifges(6,5) - t738 * Ifges(6,6);
t408 = qJ(6) * t930 + t472;
t554 = Ifges(7,3) * t742 + t813 * t739;
t555 = t739 * Ifges(7,3) - t742 * t813;
t556 = Ifges(6,3) * t742 + t814 * t739;
t557 = t739 * Ifges(6,3) - t742 * t814;
t611 = t1108 * t739;
t612 = t825 * t739;
t645 = mrSges(6,1) * t742 - mrSges(6,3) * t932;
t648 = -mrSges(7,2) * t742 + mrSges(7,3) * t930;
t649 = -mrSges(6,2) * t742 + mrSges(6,3) * t930;
t663 = -pkin(2) + t812;
t668 = t732 - t950;
t673 = -t739 * Ifges(5,2) - t989;
t683 = t739 * Ifges(4,1) + t731;
t746 = (-t679 / 0.4e1 + t684 / 0.4e1 + t671 / 0.4e1 - t674 / 0.4e1 + t556 / 0.4e1 + t554 / 0.4e1 + t850 * t741 + t847 * t738) * t606 + (-t680 / 0.4e1 - t683 / 0.4e1 + t672 / 0.4e1 + t673 / 0.4e1 - t557 / 0.4e1 - t555 / 0.4e1) * t605 + t1119 / 0.2e1 + (t281 * t668 + t419 * t663) * t1093 + (t102 * t405 + t131 * t603 + t165 * t604 + t366 * t92 + t368 * t61 + t408 * t75) * t1089 + (t113 * t462 + t114 * t463 - t204 * t685 + t236 * t686 + t471 * t94 + t472 * t95) * t1091 + t472 * t1069 + t408 * t1070 - t612 * t1071 - t611 * t1072 + t248 * t1048 + t320 * t1060 + t397 * t1061 + t668 * t1062 + t366 * t1064 + t319 * t1066 + t663 * t418 / 0.2e1 + t75 * t648 / 0.2e1 + t95 * t649 / 0.2e1 + t94 * t645 / 0.2e1 + t463 * t399 / 0.2e1 - pkin(2) * t420 / 0.2e1 + t419 * t1033 + t114 * t1037 + t102 * t1038 + t113 * t1039 + t92 * t1041 + t61 * t1043 + t236 * t1044 + t165 * t1045 - t369 * t1047 - t370 * t1027 + t249 * t1028 + t405 * t398 / 0.2e1 + (t560 / 0.4e1 + t558 / 0.4e1) * t474 + t1118 / 0.2e1 + (-t562 / 0.4e1 - t564 / 0.4e1) * t475;
t748 = (-pkin(3) * t348 - qJ(4) * t346) * t1094 + (qJ(4) * t282 - t1078 * t803) * t1092 + (t101 * t660 + t108 * t659 + t188 * t723) * t1090 + t826 * t1076 + t188 * t1032 + t282 * t1031 + t346 * mrSges(5,3) / 0.2e1 - t348 * mrSges(5,2) / 0.2e1 - t406 * mrSges(4,1) / 0.2e1 + t407 * mrSges(4,2) / 0.2e1 + t437 * t1036 + t439 * t1035 + t824 * t1026 - t1110 * t533 / 0.4e1 + t534 * t1131;
t752 = (-t268 / 0.4e1 - t267 / 0.4e1) * t741 + (-t269 / 0.4e1 - t270 / 0.4e1) * t738 + (-t334 / 0.2e1 + t286 / 0.2e1) * mrSges(5,1) + ((t286 - t334) * t1093 - t495 / 0.2e1 + t493 / 0.2e1 + t905 * t606) * pkin(9) + t177 / 0.4e1 + t178 / 0.4e1 - t354 / 0.4e1 + t356 / 0.4e1 - t422 / 0.4e1 + t424 / 0.4e1;
t754 = t861 * t741 + t859 * t738 + ((t285 + t335) * t1093 - t494 / 0.2e1 + t492 / 0.2e1 + t905 * t605) * pkin(9) + t265 / 0.4e1 + t266 / 0.4e1 + t353 / 0.4e1 - t355 / 0.4e1 - t423 / 0.4e1 + t425 / 0.4e1;
t1 = (t922 / 0.2e1 + mrSges(6,3) * t1073 + mrSges(7,3) * t1074 + t901 * t534 + t898 * t533) * t738 + (t440 * t1016 + t101 * t1080 + t123 * mrSges(6,3) / 0.2e1 - t903 * t534 - t901 * t533) * t741 + t748 + t746 + (t827 + (-Ifges(5,1) / 0.2e1 - Ifges(4,3) / 0.2e1 + pkin(3) * mrSges(5,2) / 0.2e1 + mrSges(5,3) * t1076) * t740) * t737 + ((-Ifges(5,5) / 0.2e1 + 0.3e1 / 0.4e1 * Ifges(4,6)) * t934 + (qJ(4) * t874 + t855) * mrSges(5,1) + t754) * t739 + ((0.3e1 / 0.4e1 * Ifges(5,4) + t1077 / 0.2e1 - t676 / 0.4e1 - t675 / 0.4e1 - Ifges(4,5) / 0.2e1 + (-Ifges(6,5) / 0.4e1 - Ifges(7,5) / 0.4e1) * t741 + (Ifges(6,6) / 0.4e1 + Ifges(7,6) / 0.4e1) * t738) * t934 + t752) * t742;
t788 = t555 / 0.2e1 + t557 / 0.2e1 - t673 / 0.2e1 + t683 / 0.2e1;
t17 = -pkin(2) * t670 + t366 * t644 + t368 * t646 + t405 * t648 + t408 * t650 + t462 * t645 + t463 * t649 + t471 * t647 + t472 * t651 + t603 * t613 - t604 * t611 - t686 * t612 - t685 * t614 + t668 * t665 + (m(5) * t668 + t669) * t663 + m(6) * (t462 * t471 + t463 * t472 - t685 * t686) + m(7) * (t366 * t368 + t405 * t408 + t603 * t604) + (t680 / 0.2e1 - t672 / 0.2e1 + t852 * t741 + t849 * t738 + t788) * t742 + (t684 / 0.2e1 - t674 / 0.2e1 + t554 / 0.2e1 + t556 / 0.2e1 + t851 * t741 + t848 * t738 + t845) * t739;
t808 = t1 * qJD(1) + t17 * qJD(2);
t615 = t666 * t742;
t616 = t667 * t742;
t719 = Ifges(7,6) * t931;
t617 = -Ifges(7,5) * t925 + t719;
t720 = Ifges(6,6) * t931;
t618 = -Ifges(6,5) * t925 + t720;
t619 = t742 * t677;
t620 = t742 * t678;
t621 = t742 * t681;
t622 = t742 * t682;
t22 = t404 * t650 + t462 * t651 - t463 * t647 - t604 * t615 - t686 * t616 + (t617 / 0.2e1 + t618 / 0.2e1) * t739 + (t839 - t646) * t405 + ((t983 + t976 + t619 / 0.2e1 + t620 / 0.2e1 - t848) * t741 + (t979 + t975 + t621 / 0.2e1 + t622 / 0.2e1 + (-t613 - t1012) * pkin(5) + t851) * t738) * t742;
t764 = (-t1012 / 0.2e1 - t613 / 0.2e1) * pkin(5) + t621 / 0.4e1 + t622 / 0.4e1 + t850;
t787 = -t619 / 0.4e1 - t620 / 0.4e1 + t847;
t747 = (t251 / 0.4e1 + t250 / 0.4e1) * t739 + (t617 / 0.4e1 + t618 / 0.4e1) * t606 + (-t983 / 0.2e1 - t976 / 0.2e1 + t787) * t474 + (t975 / 0.2e1 + t979 / 0.2e1 + t764) * t475 - t615 * t1072 - t616 * t1071 + t317 * t1063 + t405 * t1068 + t318 * t1061 + t463 * t1067 + t246 * t1047 + t247 * t1027 + t74 * t1038 + t75 * t1042 + t94 * t1037 + t95 * t1040;
t789 = -t252 / 0.4e1 - t253 / 0.4e1 - t859;
t790 = t254 / 0.4e1 + t255 / 0.4e1 - t861;
t753 = (t1005 / 0.2e1 + t1007 / 0.2e1 + t789) * t741 + (t1004 / 0.2e1 + t1006 / 0.2e1 + (-t248 / 0.2e1 - t1014 / 0.2e1) * pkin(5) - t790) * t738;
t761 = -pkin(5) * t439 / 0.2e1 + mrSges(7,1) * t1075 + mrSges(7,2) * t1074 + t123 * t1087 + mrSges(6,2) * t1073;
t5 = t747 - t900 * t534 - t896 * t533 + (-t894 * t934 + t753) * t742 + (-t982 / 0.2e1 + t980 / 0.2e1 - t978 / 0.2e1 + t977 / 0.2e1 + pkin(5) * t1075) * m(7) + t761;
t807 = t5 * qJD(1) + t22 * qJD(2);
t801 = t462 * t738 - t463 * t741;
t802 = -t405 * t741 + t943;
t749 = (-pkin(9) * t890 - t281 * t739 - t606 * t663) * t1094 + (t606 * t801 - t686 * t934 + t739 * t809) * t1092 + (-t604 * t934 + t606 * t802 + t739 * t810) * t1090 + t606 * t1033 + t739 * t1062;
t757 = t348 * t1093 + t803 * t1091 + (t108 * t738 + t947) * t1089;
t10 = (t1001 / 0.2e1 + (t742 * mrSges(5,1) + t1044 + t1045) * t743) * t737 + (t606 * t846 + t739 * t856 + t853) * t741 + ((t1067 + t1068) * t739 + (t1042 + t1040) * t606 + t854) * t738 + t749 + t757;
t786 = m(7) * t802;
t50 = (-m(5) * t663 + m(6) * t801 - t1111 * t741 + t1112 * t738 - t665 + t786) * t739;
t806 = qJD(1) * t10 - qJD(2) * t50;
t795 = t1068 + t834;
t797 = t912 + t907;
t13 = (-t474 * t904 + t606 * t906 - t856) * t741 + (t320 / 0.2e1 + t904 * t475 + t797 * t606 - t795) * t738;
t794 = -t839 / 0.2e1 + t1041;
t39 = (mrSges(7,3) * t841 + t828) * t742 + (t739 * t906 - t846) * t741 + (t739 * t797 + t1039 + t794) * t738;
t805 = t13 * qJD(1) + t39 * qJD(2);
t804 = t113 * t741 + t114 * t738;
t800 = -t659 * t741 + t660 * t738;
t799 = -t675 / 0.2e1 - t676 / 0.2e1 + t1077;
t109 = (t786 - t926 + t933) * t742;
t755 = (t366 * t475 + t405 * t474 + t742 * t810) * t1090 + t650 * t1059 + t475 * t1042;
t774 = t188 * t1089 - t970 / 0.2e1 + t969 / 0.2e1;
t791 = t928 / 0.2e1 + t319 * t1025;
t24 = t742 * t791 + t755 + t774;
t798 = -qJD(1) * t24 + qJD(2) * t109;
t793 = -t942 / 0.2e1 - t940 / 0.2e1;
t186 = (t875 + t793) * m(7);
t637 = (-0.1e1 / 0.2e1 + t841) * m(7);
t796 = qJD(1) * t186 - qJD(3) * t637;
t756 = (m(7) * t1066 + t1043) * pkin(5) + mrSges(7,1) * t1066 + t408 * t1084 + mrSges(6,1) * t1060 + t472 * t1086;
t769 = -t1026 * t615 + t1035 * t650 - t1076 * t616;
t776 = -t898 * t741 - t1131;
t15 = t604 * t1128 + t794 * t659 + (-t1078 * t828 + t894) * t742 + (t686 * t1087 + t651 * t1016 + t895 * t739 + (mrSges(7,3) * t1035 + t776) * t742 + t764) * t741 + (mrSges(6,2) * t1027 + mrSges(7,2) * t1047 + t647 * t1015 + t899 * t739 + (t1063 - t366 / 0.2e1) * mrSges(7,3) + (-t963 / 0.2e1 - t1117 + (t666 / 0.2e1 + t1011 / 0.2e1) * pkin(5) + t1099) * t742 + t787) * t738 + t756 + t769;
t53 = -t666 * t1008 - t723 * t1108 - qJ(4) * t825 + (-t738 * t902 + t842) * t738 + (-pkin(5) * t1011 + t1117 + (-Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t738 + t844) * t741;
t750 = (t659 * t1058 + t474 * t1035 + (-t74 / 0.2e1 + t61 / 0.2e1) * t738) * mrSges(7,3) - t1097 * t1078 + (-t738 * t901 + t776) * t474 + qJ(4) * t247 / 0.2e1 + t317 * t1034 + t723 * t246 / 0.2e1;
t768 = mrSges(7,1) * t1079 + t102 * t1083 + t1085 * t114 + t1087 * t113;
t7 = (t741 * t901 - t1099) * t475 + (mrSges(6,1) * t1071 + t790) * t741 + (t1084 * t131 + t1086 * t204 + t789) * t738 + t894 * t605 + t750 + (-t738 * t899 - t741 * t895) * t606 + t728 * t1072 + t795 * t659 + (t1065 + (t944 / 0.2e1 - t941 / 0.2e1 + t1079) * m(7) + t1104) * pkin(5) + t768;
t779 = t7 * qJD(1) - t15 * qJD(2) - t53 * qJD(3);
t710 = -0.2e1 * t889;
t751 = -t907 * t474 - t906 * t475 + (t1032 + t1031) * t934 - t712 + (t710 + t335) * t1093 + (t237 + t710) * t1091 + (-t1009 + (-qJ(4) - t723) * t934 + t800 * t606 + t237) * t1089;
t758 = t335 * t1094 + t804 * t1092 + (t102 * t738 + t741 * t92) * t1090;
t19 = (-t397 / 0.2e1 + t1065) * t741 + (-t398 / 0.2e1 - t399 / 0.2e1) * t738 + t751 + t758;
t429 = t1011 + mrSges(5,3) + t1003 * t741 + (mrSges(7,1) + mrSges(6,1)) * t738 + (m(6) + m(5)) * qJ(4);
t759 = m(6) * t1027 + (t739 * t800 + t604) * t1089;
t760 = (t471 * t741 + t472 * t738) * t1092 + (t368 * t741 + t408 * t738) * t1090;
t45 = (-t644 / 0.2e1 - t645 / 0.2e1 + t907 * t742) * t741 + (-t648 / 0.2e1 - t649 / 0.2e1 - t906 * t742) * t738 + t759 + t760;
t778 = qJD(1) * t19 + qJD(2) * t45 + qJD(3) * t429;
t766 = (t474 * t659 + t475 * t660 - t61 * t741 - t738 * t75) * t1089;
t29 = (t1068 - t973 / 0.2e1 + t967 / 0.2e1) * t741 + (-t317 / 0.2e1 - t974 / 0.2e1 - t966 / 0.2e1) * t738 + t766 - t1013 / 0.2e1;
t361 = m(7) * (-t659 * t738 - t660 * t741) + t918 * mrSges(7,3);
t763 = m(7) * ((-t659 * t742 - t366) * t741 + (t660 * t742 - t405) * t738);
t70 = (-t959 / 0.2e1 + t1041) * t741 + (t958 / 0.2e1 + t1038) * t738 + t908 - t763 / 0.2e1;
t777 = qJD(1) * t29 - qJD(2) * t70 + qJD(3) * t361;
t135 = t475 * t893 - t464;
t518 = (t738 * t893 + t957) * t742;
t624 = -m(7) * t1008 - t1108;
t775 = qJD(1) * t135 + qJD(2) * t518 + qJD(3) * t624;
t765 = t738 * t797 + t741 * t906;
t762 = -qJ(4) * mrSges(5,1) + t738 * t842 + t741 * t844;
t636 = t1090 * t918 + t1089;
t187 = (t940 + t942) * t1089 + m(7) * t875;
t71 = t763 / 0.2e1 + t650 * t1025 + t646 * t1022 + t908 + (t962 / 0.2e1 + t1128) * t739;
t42 = (m(5) * pkin(9) - t738 * t906 + t741 * t907 + mrSges(5,1)) * t742 + t759 - t760 + (t648 + t649) * t1023 + (t644 + t645) * t1020;
t40 = t839 * t1023 + t926 / 0.2e1 - t933 / 0.2e1 + t651 * t1020 + t765 * t739 + t1133 * mrSges(7,3) + t1101;
t28 = t766 + t317 * t1025 + t319 * t1022 + t1013 / 0.2e1 + mrSges(7,2) * t882 + mrSges(7,1) * t881 + t793 * mrSges(7,3);
t25 = t319 * t869 + t742 * t866 - t755 + t774;
t18 = t751 - t758 - t968 + (t398 + t399) * t1023 + (t396 + t397) * t1020;
t16 = pkin(5) * t931 * t1032 + t404 * t888 - t769 + t651 * t864 + t894 * t742 + (t921 * t659 + (t604 * t741 - t723 * t931) * pkin(5)) * t1089 + t943 * t1080 + t911 * t1034 + t646 * t1036 + t1008 * t1045 + t1108 * t1047 + t825 * t1027 + t869 * t963 + t756 + t1110 * t931 / 0.4e1 + t925 * t1131 + (-t622 - t621) * t1019 - t1101 * t1078 + (t738 * t900 + t741 * t896 + t1100) * t739 + (-t783 - t782 + t561 + t559) * t1021 + (-t785 - t784 - t620 - t619 + t565 + t563) * t1024;
t14 = mrSges(7,3) * t1103 + t606 * t765 + t738 * t834 + t1097 + t791;
t11 = mrSges(5,2) * t876 + t854 * t738 + t853 * t741 - t749 + t757 + t1112 * t882 + t1111 * t881 + (t614 + t613) * t875 + t1113 * t932 / 0.2e1 + (-t927 / 0.2e1 + t866) * t739;
t8 = t92 * t912 - t985 / 0.2e1 - t986 / 0.2e1 + t750 + t1002 * t659 * t1089 + t825 * t1071 + t1108 * t1072 + t319 * t1036 - t768 + t1123 * t882 + t1122 * t936 / 0.2e1 + t1100 * t606 + t843 * t475 + (t821 + t819) * t475 / 0.4e1 + t1115 * t1021 + (t255 + t254) * t1019 + ((-t941 + t944) * t1089 + t1064 + t1104) * pkin(5) + (t253 + t252 + t1114) * t1024;
t6 = t747 + t988 / 0.2e1 + t987 / 0.2e1 + t993 / 0.2e1 + t992 / 0.2e1 + t753 * t742 + (t977 - t978 + t980 - t982) * t1089 + t101 * t912 - t761 + t1120 * t830;
t2 = t108 * t888 + t440 * t864 - pkin(3) * t781 / 0.2e1 - t748 + (mrSges(5,1) * t855 + Ifges(4,6) * t873 + t754) * t739 + t780 * t1076 + t947 * t1081 - t922 * t1023 + t746 + t737 * t827 + Ifges(4,5) * t830 + Ifges(5,4) * t831 + Ifges(5,5) * t832 + Ifges(4,6) * t833 + t1121 * t876 + t803 * t1082 + t1106 * t1024 + t1105 * t1019 + (t752 + (Ifges(5,4) + t676 + t675) * t873) * t742;
t20 = [-qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t12 + qJD(5) * t9 + qJD(6) * t33, t2 * qJD(3) + t11 * qJD(4) + t6 * qJD(5) + t25 * qJD(6) - t984 + ((-t609 * mrSges(4,1) - t346 * mrSges(5,1) + t407 * mrSges(4,3) + t1130 * t1022 + t1129 * t1025) * t742 - t1098 + m(4) * (-pkin(2) * t609 + (-t406 * t739 + t407 * t742) * pkin(9)) + (t348 * mrSges(5,1) + t609 * mrSges(4,2) - t406 * mrSges(4,3)) * t739 + (t739 * t900 + t848) * t534 + (t739 * t896 + t851) * t533 + 0.2e1 * (t430 * t663 + (-t346 * t742 + t348 * t739) * pkin(9)) * t1093 + 0.2e1 * (t101 * t366 + t108 * t405 + t188 * t604) * t1089 + 0.2e1 * (t123 * t462 + t124 * t463 + t282 * t686) * t1091 + t430 * t665 + t108 * t650 + t124 * t651 + t101 * t646 + t123 * t647 + t188 * t613 + t282 * t614 + t462 * t440 + t463 * t438 + t405 * t437 + t366 * t439 + ((-Ifges(3,6) - t1124 * t742 - t1134 * t739 + (-t1126 * t742 + t1127 * t739) * pkin(9)) * t740 + (Ifges(3,5) + (-pkin(2) * mrSges(4,1) - t663 * mrSges(5,2) + t739 * t897 + t845) * t739 + (-pkin(2) * mrSges(4,2) - t663 * mrSges(5,3) + (-t956 / 0.2e1 - t955 / 0.2e1 - t961 / 0.2e1 - t960 / 0.2e1 - t897) * t742 + (-Ifges(5,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + t894) * t739 + t788) * t742) * t743) * t737 + t686 * t826 + t604 * t824) * qJD(2), t981 + t2 * qJD(2) + t18 * qJD(4) + t8 * qJD(5) + t28 * qJD(6) + (-t113 * mrSges(6,3) - t92 * mrSges(7,3) - t1078 * t397 - t857) * t916 + (-t114 * mrSges(6,3) - t102 * mrSges(7,3) - t1078 * t399 - t858) * t917 + ((qJ(4) * t236 - t1078 * t804) * t1091 + (t102 * t659 + t165 * t723 + t660 * t92) * t1089 + (-pkin(3) * t335 - qJ(4) * t334) * t1093) * t1095 + (-qJ(4) * t370 + t1126 * t334 + t1127 * t335 + t165 * t666 + t236 * t667 - t723 * t369 + t660 * t396 + t659 * t398 + t799 * t605 + t762 * t606 + t1102) * qJD(3), qJD(2) * t11 + qJD(3) * t18 + qJD(5) * t14 + qJD(6) * t187 + t949, t954 + t6 * qJD(2) + t8 * qJD(3) + t14 * qJD(4) + (-mrSges(6,1) * t95 - mrSges(7,1) * t75 - mrSges(6,2) * t94 - mrSges(7,2) * t74 + (-m(7) * t75 - t974) * pkin(5) + t251 + t250) * qJD(5), qJD(2) * t25 + qJD(3) * t28 + qJD(4) * t187 + t948; qJD(3) * t1 - qJD(4) * t10 + qJD(5) * t5 - qJD(6) * t24 + t984, qJD(3) * t17 + qJD(4) * t50 + qJD(5) * t22 + qJD(6) * t109, t42 * qJD(4) + t16 * qJD(5) + t71 * qJD(6) + (m(6) * qJ(4) * t1028 + (t368 * t660 + t408 * t659 + t603 * t723) * t1089) * t1095 + (-t368 * mrSges(7,3) - t1078 * t645 + t471 * t840 - t849) * t916 + (-t408 * mrSges(7,3) - t1078 * t649 + t472 * t840 + t852) * t917 + t808 + (-qJ(4) * t612 + t603 * t666 - t723 * t611 + t660 * t644 + t659 * t648 - t685 * t667 + (-Ifges(5,4) - t799) * t742 + (-Ifges(4,6) + t762) * t739 + (m(5) * t812 - mrSges(4,1) * t742 + mrSges(4,2) * t739 + t665) * pkin(9) + t1107) * qJD(3), qJD(3) * t42 + qJD(5) * t40 - t806, t16 * qJD(3) + t40 * qJD(4) + (-mrSges(6,1) * t463 - mrSges(6,2) * t462 - mrSges(7,2) * t404 - t405 * t893 + t719 + t720) * qJD(5) + t837 * t742 * t914 + t807, qJD(3) * t71 + t798; -qJD(2) * t1 + qJD(4) * t19 + qJD(5) * t7 + qJD(6) * t29 - t981, qJD(4) * t45 - qJD(5) * t15 - qJD(6) * t70 - t808, qJD(4) * t429 - qJD(5) * t53 + qJD(6) * t361, qJD(6) * t636 + t778 (-mrSges(7,2) * t660 - t659 * t893) * qJD(5) + (mrSges(6,2) * t1078 - t1122) * t914 + (mrSges(6,1) * t1078 + t837) * t915 + t779, qJD(4) * t636 + t777; qJD(2) * t10 - qJD(3) * t19 - qJD(5) * t13 - qJD(6) * t186 - t949, -qJD(3) * t45 - qJD(5) * t39 + t806, qJD(6) * t637 - t778, 0, -t1003 * t914 + (-mrSges(6,1) - t893) * t915 - t805, -t796; -qJD(2) * t5 - qJD(3) * t7 + qJD(4) * t13 + qJD(6) * t135 - t954, qJD(3) * t15 + qJD(4) * t39 + qJD(6) * t518 - t807, t624 * qJD(6) - t779, t805, 0, t775; qJD(2) * t24 - qJD(3) * t29 + qJD(4) * t186 - qJD(5) * t135 - t948, qJD(3) * t70 - qJD(5) * t518 - t798, -qJD(4) * t637 - qJD(5) * t624 - t777, t796, -t775, 0;];
Cq  = t20;
