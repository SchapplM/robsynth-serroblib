% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:10
% EndTime: 2019-12-31 21:20:00
% DurationCPUTime: 38.84s
% Computational Cost: add. (88843->898), mult. (116182->1163), div. (0->0), fcn. (122924->8), ass. (0->528)
t1222 = -Icges(5,4) + Icges(4,5);
t1221 = Icges(5,5) - Icges(4,6);
t756 = qJ(2) + qJ(3);
t742 = sin(t756);
t743 = cos(t756);
t1225 = t1221 * t742 + t1222 * t743;
t1129 = m(6) / 0.2e1;
t1130 = m(5) / 0.2e1;
t759 = sin(qJ(1));
t758 = sin(qJ(2));
t1057 = pkin(2) * t758;
t1053 = pkin(8) * t742;
t757 = sin(qJ(5));
t760 = cos(qJ(5));
t823 = rSges(6,1) * t757 + rSges(6,2) * t760;
t588 = rSges(6,3) * t742 - t743 * t823;
t1055 = pkin(3) * t742;
t695 = -qJ(4) * t743 + t1055;
t828 = t588 + t695 + t1053;
t790 = t828 + t1057;
t455 = t790 * t759;
t762 = cos(qJ(1));
t457 = t790 * t762;
t1048 = rSges(5,2) * t742;
t822 = rSges(5,3) * t743 + t1048;
t895 = t695 - t822;
t831 = t895 + t1057;
t506 = t831 * t759;
t508 = t831 * t762;
t1024 = qJ(4) * t742;
t995 = t742 * t759;
t851 = rSges(5,1) * t762 - rSges(5,3) * t995;
t1128 = -pkin(7) - pkin(6);
t761 = cos(qJ(2));
t1056 = pkin(2) * t761;
t739 = pkin(1) + t1056;
t891 = -t762 * t1128 - t759 * t739;
t467 = (-t1024 + (rSges(5,2) - pkin(3)) * t743) * t759 + t851 + t891;
t1042 = rSges(5,3) + qJ(4);
t1054 = pkin(3) * t743;
t740 = t759 * t1128;
t987 = t743 * t762;
t850 = t759 * rSges(5,1) - rSges(5,2) * t987;
t468 = -t740 + (t1042 * t742 + t1054 + t739) * t762 + t850;
t988 = t743 * t759;
t930 = t467 * t987 + t468 * t988;
t1108 = -rSges(6,3) - pkin(3);
t972 = t759 * t760;
t986 = t757 * t762;
t669 = t742 * t972 + t986;
t970 = t760 * t762;
t973 = t759 * t757;
t670 = -t742 * t973 + t970;
t1147 = t670 * rSges(6,1) - t669 * rSges(6,2);
t732 = pkin(8) * t988;
t684 = pkin(4) * t762 - t732;
t397 = (t1108 * t743 - t1024) * t759 + t684 + t891 + t1147;
t1052 = t759 * pkin(4);
t667 = t742 * t970 - t973;
t668 = t742 * t986 + t972;
t825 = t668 * rSges(6,1) + t667 * rSges(6,2);
t885 = pkin(8) - t1108;
t398 = t1052 - t740 + (t743 * t885 + t1024 + t739) * t762 + t825;
t936 = t397 * t987 + t398 * t988;
t994 = t742 * t762;
t943 = (-t455 * t994 + t457 * t995 + t936) * t1129 + (-t506 * t994 + t508 * t995 + t930) * t1130;
t776 = t1053 + (-qJ(4) - t823) * t743;
t724 = rSges(6,3) * t995;
t731 = pkin(3) * t995;
t890 = t724 + t731;
t431 = (t776 + t1057) * t759 + t890;
t835 = t885 * t742;
t715 = qJ(4) * t987;
t892 = (rSges(6,1) * t986 + rSges(6,2) * t970) * t743;
t868 = t715 + t892;
t432 = (-t835 - t1057) * t762 + t868;
t791 = -t1042 * t743 - t1048;
t499 = t731 + (t791 + t1057) * t759;
t889 = rSges(5,2) * t994 + rSges(5,3) * t987;
t500 = t715 + (-t1055 - t1057) * t762 + t889;
t944 = ((t431 * t762 + t432 * t759) * t742 + t936) * t1129 + ((t499 * t762 + t500 * t759) * t742 + t930) * t1130;
t24 = t944 - t943;
t1224 = t24 * qJD(1);
t1223 = Icges(5,1) + Icges(4,3);
t492 = -rSges(6,3) * t988 + t1147;
t783 = rSges(6,3) * t987 + t825;
t698 = t1024 + t1054;
t754 = t759 ^ 2;
t755 = t762 ^ 2;
t887 = t754 + t755;
t902 = t887 * t698;
t321 = t902 + (pkin(8) * t987 + t1052 + t783) * t762 + (-t492 - t684) * t759;
t1148 = t887 * t742;
t543 = t823 * t988 - t724;
t544 = -rSges(6,3) * t994 + t892;
t832 = -pkin(3) * t994 + t715;
t903 = t759 * (qJ(4) * t988 - t731) + t762 * t832;
t350 = -pkin(8) * t1148 + t759 * t543 + t762 * t544 + t903;
t469 = t828 * t759;
t1047 = rSges(6,3) * t743;
t587 = t742 * t823 + t1047;
t912 = -t587 - t698;
t470 = t759 * t912 - t732;
t471 = t828 * t762;
t829 = -pkin(8) * t743 + t912;
t472 = t829 * t762;
t136 = t321 * t350 - t469 * t470 - t471 * t472;
t390 = -t759 * (rSges(5,2) * t988 + t851) + t762 * (rSges(5,3) * t994 + t850) + t902;
t410 = t754 * t822 + t762 * t889 + t903;
t533 = t895 * t759;
t699 = -rSges(5,2) * t743 + rSges(5,3) * t742;
t894 = -t698 - t699;
t534 = t894 * t759;
t535 = t895 * t762;
t536 = t894 * t762;
t210 = t390 * t410 - t533 * t534 - t535 * t536;
t606 = rSges(4,1) * t988 - rSges(4,2) * t995 - t762 * rSges(4,3);
t849 = -rSges(4,2) * t994 + t759 * rSges(4,3);
t473 = t759 * t606 + t762 * (rSges(4,1) * t987 + t849);
t697 = rSges(4,1) * t742 + rSges(4,2) * t743;
t663 = t697 * t759;
t666 = t697 * t762;
t498 = -t759 * t663 - t762 * t666;
t1049 = rSges(4,1) * t743;
t700 = -rSges(4,2) * t742 + t1049;
t333 = t887 * t697 * t700 + t473 * t498;
t1220 = m(4) * t333 + m(5) * t210 + m(6) * t136;
t753 = t762 * pkin(6);
t917 = -t759 * (pkin(1) * t759 - t753 + t891) + t762 * (-t759 * pkin(6) - t740 + (-pkin(1) + t739) * t762);
t264 = t321 + t917;
t106 = t264 * t350 - t455 * t470 - t457 * t472;
t340 = t390 + t917;
t174 = t340 * t410 - t506 * t534 - t508 * t536;
t787 = t697 + t1057;
t1150 = t787 * t762;
t1151 = t787 * t759;
t379 = t473 + t917;
t258 = t379 * t498 + (t1150 * t762 + t1151 * t759) * t700;
t1219 = m(4) * t258 + m(5) * t174 + m(6) * t106;
t1114 = t742 / 0.2e1;
t1184 = t743 / 0.2e1;
t483 = Icges(6,5) * t668 + Icges(6,6) * t667 + Icges(6,3) * t987;
t1033 = Icges(6,4) * t757;
t817 = Icges(6,2) * t760 + t1033;
t779 = -Icges(6,6) * t742 + t743 * t817;
t540 = t779 * t762;
t1032 = Icges(6,4) * t760;
t819 = Icges(6,1) * t757 + t1032;
t780 = -Icges(6,5) * t742 + t743 * t819;
t542 = t780 * t762;
t813 = Icges(6,5) * t757 + Icges(6,6) * t760;
t778 = -Icges(6,3) * t742 + t743 * t813;
t1034 = Icges(6,4) * t668;
t486 = Icges(6,2) * t667 + Icges(6,6) * t987 + t1034;
t643 = Icges(6,4) * t667;
t489 = Icges(6,1) * t668 + Icges(6,5) * t987 + t643;
t808 = -t486 * t760 - t489 * t757;
t794 = t778 * t762 - t808;
t239 = (-t540 * t760 - t542 * t757 + t483) * t743 + t794 * t742;
t583 = Icges(6,6) * t743 + t742 * t817;
t585 = Icges(6,5) * t743 + t742 * t819;
t1014 = t778 * t742;
t581 = Icges(6,3) * t743 + t742 * t813;
t805 = t780 * t757 + t779 * t760;
t792 = t581 - t805;
t773 = t743 * t792 + t1014;
t271 = t667 * t583 + t668 * t585 + t762 * t773;
t1188 = -t583 * t760 - t585 * t757 - t778;
t289 = t1188 * t743 + t792 * t742;
t1016 = t483 * t742;
t330 = t743 * t808 + t1016;
t485 = -Icges(6,5) * t670 + Icges(6,6) * t669 + Icges(6,3) * t988;
t1015 = t485 * t742;
t645 = Icges(6,4) * t670;
t488 = Icges(6,2) * t669 + Icges(6,6) * t988 - t645;
t644 = Icges(6,4) * t669;
t490 = Icges(6,1) * t670 - Icges(6,5) * t988 - t644;
t1187 = t488 * t760 - t490 * t757;
t331 = t1187 * t743 - t1015;
t360 = -t667 * t779 - t668 * t780 - t778 * t987;
t362 = -t669 * t779 + t670 * t780 - t778 * t988;
t389 = t743 * t805 - t1014;
t1218 = t289 * t1114 + t389 * t1184 - (t330 + t360) * t994 / 0.4e1 - (-t331 + t362) * t995 / 0.4e1 + (t239 + t271) * t987 / 0.4e1;
t650 = (-Icges(6,5) * t760 + Icges(6,6) * t757) * t743;
t1001 = t742 * t650;
t516 = rSges(6,1) * t667 - rSges(6,2) * t668;
t517 = rSges(6,1) * t669 + rSges(6,2) * t670;
t653 = (Icges(6,2) * t757 - t1032) * t743;
t658 = (-Icges(6,1) * t760 + t1033) * t743;
t163 = (-(t658 / 0.2e1 + t779 / 0.2e1) * t757 - (-t780 / 0.2e1 + t653 / 0.2e1) * t760) * t743 + m(6) * (-t397 * t517 + t398 * t516) + t1001 / 0.2e1;
t1217 = t163 * qJD(1);
t1216 = t1225 * t762;
t1215 = t1223 * t759 + t1216;
t1208 = t1221 * t995 + t1222 * t988 - t1223 * t762;
t1027 = Icges(5,6) * t743;
t738 = Icges(4,4) * t743;
t1214 = t1027 + t738 + (Icges(4,1) + Icges(5,2)) * t742;
t941 = (-t469 * t994 + t471 * t995 + t936) * t1129 + (-t533 * t994 + t535 * t995 + t930) * t1130;
t461 = t759 * t776 + t890;
t462 = -t762 * t835 + t868;
t520 = t759 * t791 + t731;
t521 = t832 + t889;
t942 = ((t461 * t762 + t462 * t759) * t742 + t936) * t1129 + ((t520 * t762 + t521 * t759) * t742 + t930) * t1130;
t1163 = t941 - t942;
t1213 = t1163 * qJD(1);
t1035 = Icges(4,4) * t742;
t694 = Icges(4,1) * t743 - t1035;
t599 = Icges(4,5) * t759 + t694 * t762;
t686 = Icges(5,3) * t742 - t1027;
t600 = Icges(5,5) * t759 + t686 * t762;
t1212 = -t599 * t988 - t600 * t995;
t596 = Icges(4,4) * t988 - Icges(4,2) * t995 - Icges(4,6) * t762;
t601 = Icges(5,5) * t762 + Icges(5,6) * t988 - Icges(5,3) * t995;
t1211 = t596 + t601;
t722 = Icges(4,4) * t995;
t598 = Icges(4,1) * t988 - Icges(4,5) * t762 - t722;
t717 = Icges(5,6) * t995;
t603 = Icges(5,4) * t762 + Icges(5,2) * t988 - t717;
t1210 = t598 + t603;
t309 = t483 * t988 + t669 * t486 - t670 * t489;
t310 = t485 * t988 + t488 * t669 + t490 * t670;
t810 = t309 * t762 + t310 * t759;
t145 = t362 * t742 + t743 * t810;
t181 = t309 * t759 - t310 * t762;
t307 = t483 * t987 + t667 * t486 + t668 * t489;
t308 = t485 * t987 + t667 * t488 - t668 * t490;
t1196 = t307 * t759 - t308 * t762;
t861 = t988 / 0.4e1;
t863 = -t988 / 0.4e1;
t1209 = (t861 + t863) * t1196;
t1113 = -t759 / 0.2e1;
t1110 = -t762 / 0.2e1;
t1028 = Icges(5,6) * t742;
t1207 = -Icges(5,3) * t743 - t1028 + t694;
t691 = Icges(4,2) * t743 + t1035;
t1206 = -Icges(5,2) * t743 + t1028 + t691;
t1205 = -t1221 * t743 + t1222 * t742;
t1204 = -t1215 * t762 - t1212;
t1203 = t1210 - t717 - t722 + (-Icges(4,2) - Icges(5,3)) * t988;
t718 = Icges(5,6) * t994;
t602 = Icges(5,4) * t759 - Icges(5,2) * t987 + t718;
t1202 = Icges(5,3) * t987 + t691 * t762 - t599 + t602 + t718;
t1201 = t1214 * t759 + t1211;
t692 = -Icges(4,2) * t742 + t738;
t597 = Icges(4,6) * t759 + t692 * t762;
t1200 = -t1214 * t762 - t597 + t600;
t1199 = -t1208 * t759 - t598 * t987 + t601 * t994;
t1190 = t1215 * t759 + t599 * t987 + t600 * t994;
t1198 = t692 + t1214;
t1109 = t762 / 0.2e1;
t1112 = t759 / 0.2e1;
t1171 = t1110 * t1196 + t181 * t1113;
t1189 = t596 * t742 - t603 * t743;
t1191 = t597 * t742 + t602 * t743 - t1208;
t1193 = -t597 * t994 - t602 * t987 + t1190;
t1194 = t596 * t994 - t603 * t987 + t1199;
t1195 = -t597 * t995 - t602 * t988 + t1204;
t771 = (t1193 * t759 + t1194 * t762) * t1109 + (t1196 + t1190 * t759 + ((t1189 + t1215) * t762 + t1195 + t1199 + t1212) * t762) * t1110 + (-t181 + (t1191 * t759 - t1194 + t1195 - t1204) * t759 + ((t1191 + t1208) * t762 + (-t598 * t743 + t601 * t742 + t1189) * t759 - t1190 + t1193) * t762) * t1112 - t1171;
t811 = t307 * t762 + t759 * t308;
t1197 = t360 * t742 + t743 * t811;
t413 = t492 * t742 + t588 * t988;
t1131 = m(4) / 0.2e1;
t1186 = -t742 / 0.2e1;
t1185 = -t743 / 0.2e1;
t1183 = t759 / 0.4e1;
t1182 = -t762 / 0.4e1;
t539 = t779 * t759;
t541 = t780 * t759;
t793 = t778 * t759 + t1187;
t240 = (-t539 * t760 - t541 * t757 + t485) * t743 + t793 * t742;
t270 = t583 * t669 - t585 * t670 + t759 * t773;
t1169 = t240 + t270;
t1168 = t1205 * t759;
t1167 = t1205 * t762;
t1166 = (-t1206 + t1207) * t743 + (t686 - t1198) * t742;
t1165 = t1210 * t742 + t1211 * t743;
t518 = -t606 + t891;
t519 = -t740 + (t739 + t1049) * t762 + t849;
t788 = (-t518 * t762 - t519 * t759) * t700;
t933 = t536 * t467 + t534 * t468;
t937 = t472 * t397 + t470 * t398;
t873 = (-t431 * t471 - t432 * t469 + t937) * t1129 + (-t499 * t535 - t500 * t533 + t933) * t1130 + (t788 + (t1150 * t759 - t1151 * t762) * t697) * t1131;
t874 = (-t455 * t462 - t457 * t461 + t937) * t1129 + (-t506 * t521 - t508 * t520 + t933) * t1130 + (-t1150 * t663 + t1151 * t666 + t788) * t1131;
t1164 = t873 - t874;
t1036 = Icges(3,4) * t758;
t707 = Icges(3,2) * t761 + t1036;
t710 = Icges(3,1) * t761 - t1036;
t1161 = (t710 / 0.2e1 - t707 / 0.2e1) * t758;
t931 = -t455 * t988 - t457 * t987;
t175 = t1148 * t264 + t931;
t929 = -t469 * t988 - t471 * t987;
t211 = t1148 * t321 + t929;
t928 = -t506 * t988 - t508 * t987;
t241 = t1148 * t340 + t928;
t923 = -t533 * t988 - t535 * t987;
t284 = t1148 * t390 + t923;
t1041 = (t284 + t241) * t1130 + (t211 + t175) * t1129;
t827 = -t350 * t743 + t470 * t995 + t472 * t994;
t1077 = m(6) * (t321 * t742 + t827 + t929);
t826 = -t410 * t743 + t534 * t995 + t536 * t994;
t1094 = m(5) * (t390 * t742 + t826 + t923);
t950 = t1077 / 0.2e1 + t1094 / 0.2e1;
t1160 = t950 - t1041;
t1159 = (t1200 * t759 + t1201 * t762) * t743 + (t1202 * t759 + t1203 * t762) * t742;
t996 = t742 * t743;
t893 = t887 * t996;
t1149 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (t893 - t996);
t750 = Icges(3,4) * t761;
t708 = -Icges(3,2) * t758 + t750;
t709 = Icges(3,1) * t758 + t750;
t510 = Icges(6,5) * t667 - Icges(6,6) * t668;
t925 = -Icges(6,2) * t668 + t489 + t643;
t927 = -Icges(6,1) * t667 + t1034 + t486;
t229 = t510 * t987 + t667 * t925 - t668 * t927;
t511 = Icges(6,5) * t669 + Icges(6,6) * t670;
t924 = Icges(6,2) * t670 - t490 + t644;
t926 = -Icges(6,1) * t669 + t488 - t645;
t230 = t511 * t987 + t667 * t924 - t668 * t926;
t113 = t229 * t759 - t230 * t762;
t231 = t510 * t988 + t669 * t925 + t670 * t927;
t232 = t511 * t988 + t669 * t924 + t670 * t926;
t114 = t231 * t759 - t232 * t762;
t951 = t114 * t1110 + t113 * t1112;
t869 = t742 * t340 + t928;
t870 = t742 * t264 + t931;
t952 = (t827 + t870) * t1129 + (t826 + t869) * t1130;
t1144 = t762 * t492 + t759 * t783;
t254 = t511 * t742 + (t757 * t926 - t760 * t924) * t743;
t1020 = t254 * t762;
t253 = t510 * t742 + (t757 * t927 - t760 * t925) * t743;
t1021 = t253 * t759;
t913 = -t780 + t653;
t914 = -t779 - t658;
t315 = t650 * t987 + t667 * t913 - t668 * t914;
t316 = t650 * t988 + t669 * t913 + t670 * t914;
t834 = t1021 / 0.4e1 - t1020 / 0.4e1 + t315 * t1183 + t316 * t1182;
t965 = t762 * t1197;
t980 = t759 * t145;
t1141 = t1197 * t1182 - t145 * t1183 + t965 / 0.4e1 + t980 / 0.4e1;
t641 = Icges(3,5) * t759 + t710 * t762;
t897 = -t707 * t762 + t641;
t985 = t758 * t759;
t735 = Icges(3,4) * t985;
t971 = t759 * t761;
t640 = Icges(3,1) * t971 - Icges(3,5) * t762 - t735;
t898 = -Icges(3,2) * t971 + t640 - t735;
t639 = Icges(3,6) * t759 + t708 * t762;
t899 = -t709 * t762 - t639;
t638 = Icges(3,4) * t971 - Icges(3,2) * t985 - Icges(3,6) * t762;
t900 = t709 * t759 + t638;
t1140 = (-t897 * t759 + t762 * t898) * t758 + (t899 * t759 + t762 * t900) * t761;
t769 = t1198 * t1184 + (t581 + t1207) * t1114 + (t805 + t1206) * t1186 + (t686 - t1188) * t1185;
t1136 = 0.4e1 * qJD(1);
t1135 = 2 * qJD(2);
t1133 = 2 * qJD(3);
t368 = t1144 * t743;
t414 = t588 * t987 - t742 * t783;
t871 = -t368 * t350 + t413 * t472 - t414 * t470;
t292 = (t543 * t762 - t544 * t759) * t743 + t1144 * t742;
t334 = (t587 * t759 + t492) * t743 + (-t588 * t759 - t543) * t742;
t335 = ((-t587 + t1047) * t762 + t825) * t743 + (t588 * t762 + t544) * t742;
t872 = t292 * t264 - t334 * t457 - t335 * t455;
t1126 = m(6) * (t871 + t872);
t41 = t292 * t321 - t334 * t471 - t335 * t469 + t871;
t1125 = m(6) * t41;
t415 = t516 * t762 + t759 * t517;
t214 = t264 * t415;
t243 = t321 * t415;
t665 = (-rSges(6,1) * t760 + rSges(6,2) * t757) * t743;
t1121 = m(6) * (t214 + t243 + ((t457 + t471) * t762 + (t455 + t469) * t759) * t665);
t935 = t413 * t987 - t414 * t988;
t1120 = m(6) * (-t292 * t743 + (t334 * t762 + t335 * t759 - t368) * t742 + t935);
t939 = t334 * t397 + t335 * t398;
t1119 = m(6) * (t413 * t431 - t414 * t432 + t939);
t1118 = m(6) * (t413 * t461 - t414 * t462 + t939);
t1117 = m(6) * (-t292 * t368 + t334 * t413 - t335 * t414);
t1050 = rSges(3,1) * t761;
t855 = pkin(1) + t1050;
t888 = rSges(3,2) * t985 + t762 * rSges(3,3);
t558 = -t759 * t855 + t753 + t888;
t984 = t758 * t762;
t737 = rSges(3,2) * t984;
t559 = -t737 + t855 * t762 + (rSges(3,3) + pkin(6)) * t759;
t711 = rSges(3,1) * t758 + rSges(3,2) * t761;
t682 = t711 * t759;
t683 = t711 * t762;
t1107 = m(3) * (t558 * t682 - t559 * t683);
t1100 = m(4) * (-t1150 * t519 + t1151 * t518);
t1099 = m(4) * (t518 * t663 - t519 * t666);
t1084 = m(5) * (t467 * t499 + t468 * t500);
t1082 = m(5) * (t467 * t520 + t468 * t521);
t1081 = m(5) * (-t467 * t995 + t468 * t994);
t789 = (-t397 * t762 - t398 * t759) * t665;
t1072 = m(6) * (-t455 * t516 + t457 * t517 + t789);
t1071 = m(6) * (-t469 * t516 + t471 * t517 + t789);
t1065 = m(6) * (t397 * t431 + t398 * t432);
t1064 = m(6) * (-t1148 * t368 + t935);
t1062 = m(6) * (t397 * t461 + t398 * t462);
t1061 = m(6) * (-t397 * t995 + t398 * t994);
t1060 = m(6) * (-t413 * t995 - t414 * t994);
t1059 = m(6) * (-t1148 * t665 - t415 * t743);
t806 = t516 * t759 - t517 * t762;
t1058 = m(6) * t806 * t742;
t1018 = t331 * t759;
t809 = t330 * t762 - t1018;
t156 = t389 * t742 + t743 * t809;
t775 = t743 * t794 - t1016;
t216 = t540 * t669 - t542 * t670 + t759 * t775;
t774 = t743 * t793 - t1015;
t217 = t539 * t669 - t541 * t670 + t759 * t774;
t37 = (t216 * t762 + t217 * t759 + t362) * t743 + (t270 - t810) * t742;
t218 = t667 * t540 + t668 * t542 + t762 * t775;
t219 = t667 * t539 + t668 * t541 + t762 * t774;
t38 = (t218 * t762 + t219 * t759 + t360) * t743 + (t271 - t811) * t742;
t60 = (t239 * t762 + t240 * t759 + t389) * t743 + (t289 - t809) * t742;
t14 = t1117 + (t38 * t1109 + t37 * t1112 + t156 / 0.2e1) * t743 + (-t965 / 0.2e1 - t980 / 0.2e1 + t60 / 0.2e1) * t742;
t212 = t1064 / 0.2e1;
t77 = t1120 / 0.2e1;
t48 = t212 + t77 - t1059 / 0.2e1;
t1051 = t48 * qJD(4) + t14 * qJD(5);
t343 = t1059 / 0.2e1;
t47 = t343 + t212 - t1120 / 0.2e1;
t1040 = t47 * qJD(5) + (-0.4e1 * t1149 + 0.2e1 * (t1129 + t1130) * (-t1148 * t743 + t893)) * qJD(4);
t417 = 0.4e1 * t1149;
t46 = t343 + t77 - t1064 / 0.2e1;
t1039 = t417 * qJD(4) + t46 * qJD(5);
t1023 = t239 * t759;
t1022 = t240 * t762;
t1007 = t638 * t758;
t969 = t761 * t762;
t636 = Icges(3,5) * t971 - Icges(3,6) * t985 - Icges(3,3) * t762;
t916 = -t759 * t636 - t640 * t969;
t816 = Icges(3,5) * t761 - Icges(3,6) * t758;
t637 = Icges(3,3) * t759 + t762 * t816;
t915 = t759 * t637 + t641 * t969;
t886 = qJD(2) + qJD(3);
t881 = t1121 / 0.2e1 + t951;
t862 = t988 / 0.2e1;
t859 = t987 / 0.2e1;
t853 = -t698 - t1056;
t852 = -t700 - t1056;
t576 = t641 * t971;
t839 = t762 * t637 - t576;
t836 = t639 * t758 - t636;
t833 = t887 * t1057;
t830 = -t699 + t853;
t97 = m(5) * t241 + m(6) * t175;
t137 = m(5) * t284 + m(6) * t211;
t815 = -Icges(3,5) * t758 - Icges(3,6) * t761;
t104 = t216 * t759 - t217 * t762;
t105 = t218 * t759 - t219 * t762;
t797 = (t105 + (t1168 * t759 + t1159) * t762 - t1167 * t754) * t1112 + (t104 + (t1167 * t762 + t1159) * t759 - t1168 * t755) * t1110;
t786 = -t368 * t415 + (-t413 * t762 + t414 * t759) * t665;
t782 = t1141 + t1209;
t781 = t104 * t862 + t105 * t859 + t37 * t1110 + t38 * t1112 + (-t1022 + t1023) * t1114 + (t330 * t759 + t331 * t762) * t1184 - t951 + t1171 * t742;
t777 = t1169 * t861 + t1218;
t88 = t315 * t742 + (t229 * t762 + t230 * t759) * t743;
t89 = t316 * t742 + (t231 * t762 + t232 * t759) * t743;
t772 = t114 * t862 + t113 * t859 + t156 * t1185 - t37 * t988 / 0.2e1 - t38 * t987 / 0.2e1 + t60 * t1186 + t88 * t1112 + t89 * t1110 - t1117 + (t980 + t965 - t1020 + t1021) * t1114;
t767 = t777 + t782 - t834;
t766 = t1169 * t863 - t1218 + t782 + t834;
t765 = -t1141 + t1209 + t777 + t834;
t764 = t1023 / 0.2e1 - t1022 / 0.2e1 + (t1166 * t762 + t1200 * t742 - t1202 * t743 + t1225 * t759 + t271) * t1112 + (t1166 * t759 - t1201 * t742 + t1203 * t743 - t1216 + t270) * t1110 - t771;
t763 = -t769 - t1018 / 0.2e1 + t1165 * t1113 + (t331 + t1165) * t1112;
t713 = -rSges(3,2) * t758 + t1050;
t677 = t815 * t762;
t676 = t815 * t759;
t592 = t852 * t762;
t590 = t852 * t759;
t509 = t830 * t762;
t507 = t830 * t759;
t463 = -t833 + t498;
t458 = (t829 - t1056) * t762;
t456 = -t732 + (-t587 + t853) * t759;
t425 = t742 * t516 - t665 * t987;
t424 = -t517 * t742 + t665 * t988;
t423 = -t639 * t984 + t915;
t422 = -t638 * t984 - t916;
t421 = -t639 * t985 - t839;
t396 = t806 * t743;
t387 = t1058 / 0.2e1;
t383 = -t833 + t410;
t346 = -t833 + t350;
t342 = (t1001 + (t757 * t914 - t760 * t913) * t743) * t742;
t337 = -t422 * t762 + t423 * t759;
t336 = -(-t759 * (-t640 * t761 + t1007) - t762 * t636) * t762 + t421 * t759;
t291 = t1060 / 0.2e1;
t202 = t1061 + t1081;
t164 = t243 + (t469 * t759 + t471 * t762) * t665;
t162 = t1071 / 0.2e1;
t158 = t1072 / 0.2e1;
t154 = (t421 - t576 + (t637 + t1007) * t762 + t916) * t762 + t915 * t759;
t153 = (t762 * t836 + t423 - t915) * t762 + (t759 * t836 + t422 + t839) * t759;
t149 = t214 + (t455 * t759 + t457 * t762) * t665;
t93 = t291 - t1058 / 0.2e1;
t92 = t387 + t291;
t91 = t387 - t1060 / 0.2e1;
t83 = t1118 / 0.2e1;
t79 = t1119 / 0.2e1;
t74 = t1062 + t1082 + t769 + t1099;
t61 = t1107 + t1100 + t1084 + t1065 + t769 + (t709 / 0.2e1 + t708 / 0.2e1) * t761 + t1161;
t44 = t47 * qJD(4);
t40 = t1125 / 0.2e1;
t36 = t1126 / 0.2e1;
t30 = t941 + t942;
t28 = m(6) * t164 + t951;
t27 = m(6) * t149 + t951;
t25 = t943 + t944;
t23 = t950 + t1041 - t952;
t22 = t952 - t1160;
t21 = t952 + t1160;
t20 = t797 + t1220;
t19 = t20 * qJD(3);
t18 = t797 + t1219;
t16 = t40 - t1126 / 0.2e1 + t881;
t15 = t36 - t1125 / 0.2e1 + t881;
t11 = t771 + (t153 / 0.2e1 + t336 / 0.2e1) * t759 + (t337 / 0.2e1 - t154 / 0.2e1) * t762;
t10 = t36 + t40 - t1121 / 0.2e1 + t781;
t9 = t771 + t1164;
t8 = t771 - t1164;
t7 = t162 + t83 + t765;
t6 = t767 + t83 - t1071 / 0.2e1;
t5 = t162 + t766 - t1118 / 0.2e1;
t4 = t79 + t158 + t765;
t3 = t79 + t767 - t1072 / 0.2e1;
t2 = t766 + t158 - t1119 / 0.2e1;
t1 = t764 + t873 + t874;
t12 = [t61 * qJD(2) + t74 * qJD(3) + t202 * qJD(4) + t163 * qJD(5), t61 * qJD(1) + t1 * qJD(3) + t25 * qJD(4) + t4 * qJD(5) + (m(3) * ((-t558 * t762 - t559 * t759) * t713 + (-t682 * t762 + t683 * t759) * t711) / 0.2e1 + (t518 * t592 + t519 * t590) * t1131 + (t467 * t509 + t468 * t507 - t499 * t508 - t500 * t506) * t1130 + (t397 * t458 + t398 * t456 - t431 * t457 - t432 * t455) * t1129) * t1135 + (t764 + (t758 * t899 + t761 * t897) * t1112 + t154 * t1109 + (t153 + t336) * t1113 + (-t758 * t900 + t761 * t898 + t337) * t1110 + (t754 / 0.2e1 + t755 / 0.2e1) * t816) * qJD(2), t74 * qJD(1) + t1 * qJD(2) + t764 * qJD(3) + t30 * qJD(4) + t7 * qJD(5) + ((t788 + (-t663 * t762 + t666 * t759) * t697) * t1131 + (-t520 * t535 - t521 * t533 + t933) * t1130 + (-t461 * t471 - t462 * t469 + t937) * t1129) * t1133, qJD(1) * t202 + qJD(2) * t25 + qJD(3) * t30 + qJD(5) * t92, t1217 + t4 * qJD(2) + t7 * qJD(3) + t92 * qJD(4) + (t342 + m(6) * (t397 * t424 + t398 * t425 - t413 * t517 - t414 * t516) + ((t315 / 0.2e1 + t253 / 0.2e1) * t762 + (t316 / 0.2e1 + t254 / 0.2e1) * t759) * t743) * qJD(5); (t763 - t1161 - (t709 + t708) * t761 / 0.2e1) * qJD(1) + t11 * qJD(2) + t8 * qJD(3) - t24 * qJD(4) + t2 * qJD(5) + (-t1107 / 0.4e1 - t1100 / 0.4e1 - t1084 / 0.4e1 - t1065 / 0.4e1) * t1136, t11 * qJD(1) + (m(6) * (t264 * t346 - t455 * t456 - t457 * t458) + m(5) * (t340 * t383 - t506 * t507 - t508 * t509) + m(4) * (-t1150 * t592 - t1151 * t590 + t379 * t463) + (t754 * t677 + (-t759 * t676 + t1140) * t762) * t1112 + (t755 * t676 + (-t762 * t677 + t1140) * t759) * t1110 + m(3) * ((t759 * (rSges(3,1) * t971 - t888) + t762 * (rSges(3,1) * t969 + t759 * rSges(3,3) - t737)) * (-t759 * t682 - t683 * t762) + t887 * t713 * t711) + t797) * qJD(2) + t18 * qJD(3) + t97 * qJD(4) + t27 * qJD(5), t8 * qJD(1) + t18 * qJD(2) + t22 * qJD(4) + t15 * qJD(5) + ((t136 + t106) * t1129 + (t210 + t174) * t1130 + (t258 + t333) * t1131) * t1133 + (t797 - t1220) * qJD(3), qJD(2) * t97 + qJD(3) * t22 + t1040 - t1224, t2 * qJD(1) + t27 * qJD(2) + t15 * qJD(3) + t44 + (t772 + m(6) * (-t396 * t264 - t424 * t457 - t425 * t455 + t786)) * qJD(5); t763 * qJD(1) + t9 * qJD(2) + t771 * qJD(3) + t1163 * qJD(4) + t5 * qJD(5) + (-t1099 / 0.4e1 - t1082 / 0.4e1 - t1062 / 0.4e1) * t1136, t9 * qJD(1) + t19 + t23 * qJD(4) + t16 * qJD(5) + ((t321 * t346 - t456 * t469 - t458 * t471 + t106) * t1129 + (t383 * t390 - t507 * t533 - t509 * t535 + t174) * t1130 + (t473 * t463 + (-t590 * t759 - t592 * t762) * t697 + t258) * t1131) * t1135 + (t797 - t1219) * qJD(2), qJD(1) * t771 + qJD(2) * t20 + qJD(4) * t137 + qJD(5) * t28 + t19, qJD(2) * t23 + qJD(3) * t137 + t1040 + t1213, t5 * qJD(1) + t16 * qJD(2) + t28 * qJD(3) + t44 + (t772 + m(6) * (-t396 * t321 - t424 * t471 - t425 * t469 + t786)) * qJD(5); t24 * qJD(2) - t1163 * qJD(3) + t91 * qJD(5) + (-t1061 / 0.4e1 - t1081 / 0.4e1) * t1136, t1224 + (m(6) * (-t743 * t346 + (t456 * t759 + t458 * t762) * t742 + t870) + m(5) * (-t743 * t383 + (t507 * t759 + t509 * t762) * t742 + t869) - t97) * qJD(2) + t21 * qJD(3) + t1039, -t1213 + t21 * qJD(2) + (t1077 - t137 + t1094) * qJD(3) + t1039, t886 * t417, t91 * qJD(1) + m(6) * (t396 * t743 + (t424 * t762 + t425 * t759) * t742) * qJD(5) + t886 * t46; t3 * qJD(2) + t6 * qJD(3) + t93 * qJD(4) - t1217, t3 * qJD(1) + ((-t346 * t368 + t413 * t458 - t414 * t456 - t149 + t872) * m(6) + t781) * qJD(2) + t10 * qJD(3) + t1051, t6 * qJD(1) + t10 * qJD(2) + ((t41 - t164) * m(6) + t781) * qJD(3) + t1051, qJD(1) * t93 + t48 * t886, (m(6) * (t368 * t396 + t413 * t424 - t414 * t425) + t342 * t1114 + (t88 * t1109 + t89 * t1112 + (t253 * t762 + t254 * t759) * t1114) * t743) * qJD(5) + t886 * t14;];
Cq = t12;
