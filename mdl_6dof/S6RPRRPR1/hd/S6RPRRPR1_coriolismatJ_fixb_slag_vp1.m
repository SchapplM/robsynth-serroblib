% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:33
% EndTime: 2019-03-09 04:58:16
% DurationCPUTime: 35.34s
% Computational Cost: add. (164096->880), mult. (108833->1173), div. (0->0), fcn. (115277->12), ass. (0->550)
t751 = qJ(3) + qJ(4);
t745 = pkin(11) + t751;
t738 = sin(t745);
t739 = cos(t745);
t746 = sin(t751);
t747 = cos(t751);
t1226 = Icges(5,5) * t747 + Icges(6,5) * t739 - Icges(5,6) * t746 - Icges(6,6) * t738;
t1225 = Icges(5,3) + Icges(6,3);
t750 = qJ(1) + pkin(10);
t743 = sin(t750);
t1006 = t738 * t743;
t744 = cos(t750);
t997 = t739 * t743;
t569 = rSges(6,1) * t997 - rSges(6,2) * t1006 - t744 * rSges(6,3);
t1005 = t738 * t744;
t847 = -rSges(6,2) * t1005 + rSges(6,3) * t743;
t1066 = pkin(4) * t747;
t756 = cos(qJ(3));
t749 = t756 * pkin(3);
t740 = t749 + pkin(2);
t698 = t740 + t1066;
t758 = -pkin(8) - pkin(7);
t889 = -qJ(5) + t758;
t719 = t743 * t889;
t722 = t743 * t758;
t901 = -t743 * t698 - t744 * t889;
t962 = t744 * t758;
t927 = -t743 * (t743 * t740 + t901 + t962) + t744 * (-t719 + t722 + (t698 - t740) * t744);
t996 = t739 * t744;
t318 = t743 * t569 + t744 * (rSges(6,1) * t996 + t847) + t927;
t1063 = pkin(2) - t740;
t736 = t744 * pkin(7);
t921 = -t743 * (t1063 * t743 - t736 - t962) + t744 * (-pkin(7) * t743 - t1063 * t744 - t722);
t260 = t318 + t921;
t1067 = pkin(4) * t746;
t741 = t743 ^ 2;
t742 = t744 ^ 2;
t895 = t741 + t742;
t836 = t895 * t1067;
t675 = rSges(6,1) * t738 + rSges(6,2) * t739;
t898 = rSges(6,1) * t1006 + rSges(6,2) * t997;
t913 = -t742 * t675 - t743 * t898;
t420 = -t836 + t913;
t753 = sin(qJ(3));
t1068 = pkin(3) * t753;
t831 = t1067 + t1068;
t801 = t675 + t831;
t514 = t801 * t743;
t516 = t801 * t744;
t1058 = rSges(6,1) * t739;
t676 = -rSges(6,2) * t738 + t1058;
t850 = -t676 - t1066;
t559 = t850 * t743;
t561 = t850 * t744;
t934 = -t514 * t559 - t516 * t561;
t166 = t260 * t420 + t934;
t686 = rSges(5,1) * t746 + rSges(5,2) * t747;
t786 = t686 + t1068;
t1158 = t786 * t744;
t1159 = t786 * t743;
t968 = t744 * t746;
t848 = -rSges(5,2) * t968 + rSges(5,3) * t743;
t986 = t743 * t746;
t897 = rSges(5,2) * t986 + t744 * rSges(5,3);
t967 = t744 * t747;
t985 = t743 * t747;
t448 = t743 * (rSges(5,1) * t985 - t897) + t744 * (rSges(5,1) * t967 + t848);
t350 = t448 + t921;
t644 = t686 * t743;
t645 = t686 * t744;
t480 = -t743 * t644 - t744 * t645;
t1059 = rSges(5,1) * t747;
t687 = -rSges(5,2) * t746 + t1059;
t236 = t350 * t480 + (t1158 * t744 + t1159 * t743) * t687;
t755 = cos(qJ(6));
t964 = t744 * t755;
t752 = sin(qJ(6));
t984 = t743 * t752;
t640 = t739 * t984 + t964;
t966 = t744 * t752;
t982 = t743 * t755;
t641 = t739 * t982 - t966;
t1153 = -t641 * rSges(7,1) + t640 * rSges(7,2);
t477 = -rSges(7,3) * t1006 + t1153;
t642 = -t739 * t966 + t982;
t643 = t739 * t964 + t984;
t826 = rSges(7,1) * t643 + rSges(7,2) * t642;
t478 = rSges(7,3) * t1005 + t826;
t1065 = pkin(5) * t739;
t678 = pkin(9) * t738 + t1065;
t234 = t744 * t478 + t742 * t678 + t927 + (t678 * t743 - t477) * t743;
t212 = t234 + t921;
t1056 = rSges(7,3) * t739;
t1003 = t738 * t755;
t886 = rSges(7,1) * t1003;
t679 = t743 * t886;
t1004 = t738 * t752;
t885 = rSges(7,2) * t1004;
t510 = -t679 + (t885 + t1056) * t743;
t899 = rSges(7,3) * t996 + t744 * t885;
t511 = -t744 * t886 + t899;
t696 = pkin(5) * t1006;
t697 = pkin(9) * t996;
t837 = (-pkin(5) * t1005 + t511 + t697) * t744 + (pkin(9) * t997 + t510 - t696) * t743;
t316 = -t836 + t837;
t1057 = rSges(7,1) * t755;
t825 = -rSges(7,2) * t752 + t1057;
t599 = t738 * t825 - t1056;
t1155 = pkin(5) * t738 - pkin(9) * t739 + t599;
t793 = t831 + t1155;
t434 = t793 * t743;
t436 = t793 * t744;
t600 = rSges(7,3) * t738 + t739 * t825;
t1154 = -t600 - t678;
t832 = -t1066 + t1154;
t456 = t832 * t743;
t458 = t832 * t744;
t939 = -t434 * t456 - t436 * t458;
t82 = t212 * t316 + t939;
t1224 = m(5) * t236 + m(6) * t166 + m(7) * t82;
t1119 = -t739 / 0.2e1;
t1191 = t738 / 0.2e1;
t465 = Icges(7,5) * t643 + Icges(7,6) * t642 + Icges(7,3) * t1005;
t1040 = Icges(7,4) * t755;
t819 = -Icges(7,2) * t752 + t1040;
t594 = -Icges(7,6) * t739 + t738 * t819;
t507 = t594 * t744;
t1041 = Icges(7,4) * t752;
t820 = Icges(7,1) * t755 - t1041;
t596 = -Icges(7,5) * t739 + t738 * t820;
t509 = t596 * t744;
t814 = Icges(7,5) * t755 - Icges(7,6) * t752;
t592 = -Icges(7,3) * t739 + t738 * t814;
t1042 = Icges(7,4) * t643;
t468 = Icges(7,2) * t642 + Icges(7,6) * t1005 + t1042;
t623 = Icges(7,4) * t642;
t471 = Icges(7,1) * t643 + Icges(7,5) * t1005 + t623;
t809 = -t468 * t752 + t471 * t755;
t795 = -t592 * t744 - t809;
t215 = -t795 * t739 + (t507 * t752 - t509 * t755 + t465) * t738;
t595 = Icges(7,6) * t738 + t739 * t819;
t597 = Icges(7,5) * t738 + t739 * t820;
t1017 = t592 * t739;
t593 = Icges(7,3) * t738 + t739 * t814;
t1015 = t596 * t755;
t1016 = t594 * t752;
t803 = t1015 - t1016;
t794 = t593 - t803;
t769 = t738 * t794 + t1017;
t246 = t595 * t642 + t597 * t643 + t744 * t769;
t290 = -t794 * t739 + (-t595 * t752 + t597 * t755 + t592) * t738;
t463 = Icges(7,5) * t641 - Icges(7,6) * t640 + Icges(7,3) * t1006;
t1023 = t463 * t739;
t622 = Icges(7,4) * t641;
t466 = -Icges(7,2) * t640 + Icges(7,6) * t1006 + t622;
t621 = Icges(7,4) * t640;
t470 = -Icges(7,1) * t641 - Icges(7,5) * t1006 + t621;
t1194 = t466 * t752 + t470 * t755;
t300 = t1194 * t738 + t1023;
t1022 = t465 * t739;
t301 = t738 * t809 - t1022;
t345 = t1006 * t592 - t594 * t640 + t596 * t641;
t347 = t1005 * t592 + t594 * t642 + t596 * t643;
t388 = t738 * t803 - t1017;
t1223 = t290 * t1119 + t388 * t1191 + (-t300 + t345) * t997 / 0.4e1 + (t301 + t347) * t996 / 0.4e1 + (t215 + t246) * t1005 / 0.4e1;
t646 = (-Icges(7,5) * t752 - Icges(7,6) * t755) * t738;
t1000 = t739 * t646;
t1113 = rSges(7,3) + pkin(9);
t1148 = -t1113 * t738 - t1065;
t1064 = sin(qJ(1)) * pkin(1);
t834 = t901 - t1064;
t1167 = t1148 * t743 + t1153 + t834;
t1069 = cos(qJ(1)) * pkin(1);
t853 = -t719 + t1069;
t379 = (-t1148 + t698) * t744 + t826 + t853;
t494 = -rSges(7,1) * t640 - rSges(7,2) * t641;
t495 = rSges(7,1) * t642 - rSges(7,2) * t643;
t647 = (-Icges(7,2) * t755 - t1041) * t738;
t648 = (-Icges(7,1) * t752 - t1040) * t738;
t158 = (-(t596 / 0.2e1 + t647 / 0.2e1) * t752 + (t648 / 0.2e1 - t594 / 0.2e1) * t755) * t738 - t1000 / 0.2e1 + m(7) * (-t1167 * t494 + t379 * t495);
t1222 = t158 * qJD(1);
t1221 = t1226 * t744;
t1220 = t1225 * t743 + t1221;
t1214 = Icges(5,5) * t985 + Icges(6,5) * t997 - Icges(5,6) * t986 - Icges(6,6) * t1006 - t1225 * t744;
t1043 = Icges(6,4) * t738;
t673 = Icges(6,1) * t739 - t1043;
t567 = Icges(6,5) * t743 + t673 * t744;
t1044 = Icges(5,4) * t746;
t685 = Icges(5,1) * t747 - t1044;
t591 = Icges(5,5) * t743 + t685 * t744;
t1219 = -t567 * t997 - t591 * t985;
t267 = t1006 * t463 - t466 * t640 - t470 * t641;
t268 = t465 * t1006 - t640 * t468 + t641 * t471;
t813 = t267 * t743 + t268 * t744;
t128 = -t345 * t739 + t738 * t813;
t159 = -t267 * t744 + t268 * t743;
t1207 = t599 * t1006 - t477 * t739;
t393 = t1005 * t599 + t478 * t739;
t1156 = -t1207 * t744 + t393 * t743;
t269 = t463 * t1005 + t642 * t466 - t470 * t643;
t270 = t465 * t1005 + t642 * t468 + t643 * t471;
t1205 = -t269 * t744 + t270 * t743;
t865 = t1006 / 0.4e1;
t867 = -t1006 / 0.4e1;
t1217 = (t865 + t867) * t1205;
t1117 = t743 / 0.2e1;
t1115 = t744 / 0.2e1;
t1135 = m(6) / 0.2e1;
t1152 = t743 * t434 + t436 * t744;
t1193 = -m(7) / 0.2e1;
t942 = t1152 * t1193 + (-t743 * t514 - t516 * t744) * t1135;
t1134 = m(7) / 0.2e1;
t674 = t744 * t831;
t829 = (-pkin(5) - t1057) * t738;
t868 = t697 + t899;
t422 = t744 * t829 - t674 + t868;
t787 = -t1113 * t739 - t885;
t900 = t679 + t696;
t423 = (t831 + t787) * t743 + t900;
t518 = t743 * t831 + t898;
t519 = -t675 * t744 - t674;
t947 = (-t422 * t744 + t423 * t743) * t1134 + (t518 * t743 - t519 * t744) * t1135;
t83 = t947 - t942;
t1215 = t83 * qJD(1);
t1213 = -Icges(5,5) * t746 - Icges(6,5) * t738 - Icges(5,6) * t747 - Icges(6,6) * t739;
t1212 = -t1220 * t744 - t1219;
t690 = Icges(6,4) * t1006;
t566 = Icges(6,1) * t997 - Icges(6,5) * t744 - t690;
t701 = Icges(5,4) * t986;
t590 = Icges(5,1) * t985 - Icges(5,5) * t744 - t701;
t1211 = -t1214 * t743 - t566 * t996 - t590 * t967;
t1198 = t1220 * t743 + t567 * t996 + t591 * t967;
t833 = t1067 + t1155;
t455 = t833 * t743;
t457 = t833 * t744;
t1151 = t743 * t455 + t457 * t744;
t851 = t675 + t1067;
t828 = t851 * t744;
t1160 = t744 * t828;
t558 = t851 * t743;
t940 = t1151 * t1193 + (-t743 * t558 - t1160) * t1135;
t427 = (t787 + t1067) * t743 + t900;
t428 = (t829 - t1067) * t744 + t868;
t547 = pkin(4) * t986 + t898;
t944 = (t427 * t743 - t428 * t744) * t1134 + (t547 * t743 + t1160) * t1135;
t106 = t944 - t940;
t1210 = qJD(1) * t106;
t808 = -t744 * t477 - t478 * t743;
t259 = t808 * t739 + (t510 * t744 - t511 * t743) * t738;
t389 = t494 * t743 + t495 * t744;
t209 = (t259 + t389) * t1134;
t1209 = t209 * qJD(6);
t208 = 0.2e1 * (t259 / 0.4e1 - t389 / 0.4e1) * m(7);
t1208 = t208 * qJD(6);
t890 = qJD(3) + qJD(4);
t812 = t269 * t743 + t270 * t744;
t1206 = -t347 * t739 + t738 * t812;
t723 = Icges(6,4) * t739;
t671 = -Icges(6,2) * t738 + t723;
t565 = Icges(6,6) * t743 + t671 * t744;
t737 = Icges(5,4) * t747;
t683 = -Icges(5,2) * t746 + t737;
t589 = Icges(5,6) * t743 + t683 * t744;
t1204 = -t1006 * t565 - t589 * t986 + t1212;
t564 = Icges(6,4) * t997 - Icges(6,2) * t1006 - Icges(6,6) * t744;
t588 = Icges(5,4) * t985 - Icges(5,2) * t986 - Icges(5,6) * t744;
t1203 = t1005 * t564 + t588 * t968 + t1211;
t1202 = -t1005 * t565 - t589 * t968 + t1198;
t672 = Icges(6,1) * t738 + t723;
t1201 = -t671 - t672;
t684 = Icges(5,1) * t746 + t737;
t1200 = t683 + t684;
t1199 = t565 * t738 + t589 * t746 - t1214;
t1195 = t564 * t738 + t588 * t746;
t1116 = -t744 / 0.2e1;
t1179 = t1115 * t1205 + t1117 * t159;
t767 = (t1202 * t743 + t1203 * t744) * t1115 + (t1205 + t1198 * t743 + ((t1195 + t1220) * t744 + t1204 + t1211 + t1219) * t744) * t1116 + (-t159 + (t1199 * t743 - t1203 + t1204 - t1212) * t743 + ((t1199 + t1214) * t744 + (-t566 * t739 - t590 * t747 + t1195) * t743 - t1198 + t1202) * t744) * t1117 + t1179;
t1157 = -t1167 * t744 - t379 * t743;
t1136 = m(5) / 0.2e1;
t1192 = -t738 / 0.2e1;
t1190 = t739 / 0.2e1;
t1189 = t743 / 0.4e1;
t1188 = -t744 / 0.4e1;
t506 = t594 * t743;
t508 = t596 * t743;
t796 = -t592 * t743 + t1194;
t214 = -t796 * t739 + (t506 * t752 - t508 * t755 + t463) * t738;
t245 = -t595 * t640 + t597 * t641 + t743 * t769;
t1178 = t214 + t245;
t1176 = t1213 * t743;
t1175 = t1213 * t744;
t670 = Icges(6,2) * t739 + t1043;
t682 = Icges(5,2) * t747 + t1044;
t1174 = (-t682 + t685) * t747 - t1200 * t746 + (-t670 + t673) * t739 + t1201 * t738;
t1173 = t564 * t739 + t566 * t738 + t588 * t747 + t590 * t746;
t849 = t740 + t1059;
t481 = -t743 * t849 - t1064 + t897 - t962;
t482 = t744 * t849 + t1069 - t722 + t848;
t790 = (-t481 * t744 - t482 * t743) * t687;
t472 = -t569 + t834;
t473 = (t698 + t1058) * t744 + t847 + t853;
t935 = t561 * t472 + t559 * t473;
t946 = t1167 * t458 + t456 * t379;
t872 = (-t422 * t455 - t423 * t457 + t946) * t1134 + (-t518 * t828 - t519 * t558 + t935) * t1135 + (t790 + (t1158 * t743 - t1159 * t744) * t686) * t1136;
t873 = (-t427 * t436 - t428 * t434 + t946) * t1134 + (t514 * t828 - t516 * t547 + t935) * t1135 + (-t1158 * t644 + t1159 * t645 + t790) * t1136;
t1172 = t872 - t873;
t1045 = Icges(4,4) * t753;
t706 = Icges(4,2) * t756 + t1045;
t709 = Icges(4,1) * t756 - t1045;
t1170 = (t709 / 0.2e1 - t706 / 0.2e1) * t753;
t748 = Icges(4,4) * t756;
t707 = -Icges(4,2) * t753 + t748;
t708 = Icges(4,1) * t753 + t748;
t488 = -Icges(7,5) * t640 - Icges(7,6) * t641;
t930 = -Icges(7,2) * t641 - t470 - t621;
t932 = -Icges(7,1) * t640 - t466 - t622;
t199 = t1006 * t488 - t640 * t930 + t641 * t932;
t489 = Icges(7,5) * t642 - Icges(7,6) * t643;
t929 = -Icges(7,2) * t643 + t471 + t623;
t931 = Icges(7,1) * t642 - t1042 - t468;
t200 = t1006 * t489 - t640 * t929 + t641 * t931;
t102 = -t199 * t744 + t200 * t743;
t201 = t1005 * t488 + t642 * t930 + t643 * t932;
t202 = t1005 * t489 + t642 * t929 + t643 * t931;
t103 = -t201 * t744 + t202 * t743;
t1049 = t102 * t1116 + t103 * t1117;
t230 = -t489 * t739 + (-t752 * t929 + t755 * t931) * t738;
t1030 = t230 * t743;
t229 = -t488 * t739 + (-t752 * t930 + t755 * t932) * t738;
t1031 = t229 * t744;
t907 = t596 + t647;
t908 = -t594 + t648;
t288 = t1006 * t646 - t640 * t907 + t641 * t908;
t289 = t1005 * t646 + t642 * t907 + t643 * t908;
t838 = -t1031 / 0.4e1 + t1030 / 0.4e1 + t289 * t1189 + t288 * t1188;
t978 = t744 * t1206;
t993 = t743 * t128;
t1145 = t1206 * t1188 - t128 * t1189 + t978 / 0.4e1 + t993 / 0.4e1;
t915 = -t670 * t744 + t567;
t916 = -Icges(6,2) * t997 + t566 - t690;
t917 = -t672 * t744 - t565;
t918 = t672 * t743 + t564;
t1144 = (-t915 * t743 + t744 * t916) * t738 + (t917 * t743 + t744 * t918) * t739;
t765 = -t595 * t1004 / 0.2e1 + t597 * t1003 / 0.2e1 + t670 * t1192 + (-t682 / 0.2e1 + t685 / 0.2e1) * t746 + (t592 + t673) * t1191 + t1200 * t747 / 0.2e1 + (t1016 + t593) * t1119 + (t1015 - t1201) * t1190;
t1142 = 0.4e1 * qJD(1);
t1141 = 2 * qJD(3);
t1139 = 2 * qJD(4);
t1138 = 4 * qJD(4);
t1137 = m(4) / 0.2e1;
t344 = t808 * t738;
t224 = t344 * t316;
t306 = (t599 * t743 + t510) * t739 + (t600 * t743 + t477) * t738;
t307 = (-t599 * t744 - t511) * t739 + (-t600 * t744 + t478) * t738;
t871 = t259 * t212 - t306 * t436 - t307 * t434;
t945 = t1207 * t458 - t393 * t456;
t1131 = m(7) * (t224 + t871 + t945);
t965 = t744 * t753;
t914 = t741 * (-t831 + t1068) + t744 * (pkin(3) * t965 - t674);
t296 = t837 + t914;
t824 = t259 * t234 - t306 * t457 - t307 * t455 + t945;
t1130 = m(7) * (t296 * t344 + t824);
t167 = t212 * t389;
t188 = t234 * t389;
t649 = (-rSges(7,1) * t752 - rSges(7,2) * t755) * t738;
t1127 = m(7) * (t167 + t188 + ((t436 + t457) * t744 + (t434 + t455) * t743) * t649);
t1125 = m(7) * (t1207 * t306 + t259 * t344 - t307 * t393);
t949 = t1167 * t306 + t307 * t379;
t1124 = m(7) * (t1207 * t423 - t393 * t422 + t949);
t1123 = m(7) * (t1207 * t427 - t393 * t428 + t949);
t938 = -t455 * t456 - t457 * t458;
t1120 = m(7) * (t234 * t296 + t938);
t1118 = -t743 / 0.2e1;
t1060 = rSges(4,1) * t756;
t855 = pkin(2) + t1060;
t983 = t743 * t753;
t896 = rSges(4,2) * t983 + t744 * rSges(4,3);
t512 = -t743 * t855 - t1064 + t736 + t896;
t718 = rSges(4,2) * t965;
t513 = t1069 - t718 + t855 * t744 + (rSges(4,3) + pkin(7)) * t743;
t710 = rSges(4,1) * t753 + rSges(4,2) * t756;
t664 = t710 * t743;
t665 = t710 * t744;
t1112 = m(4) * (t512 * t664 - t513 * t665);
t310 = t686 * t687 * t895 + t448 * t480;
t309 = m(5) * t310;
t1106 = m(5) * (-t1158 * t482 + t1159 * t481);
t1105 = m(5) * (t481 * t644 - t482 * t645);
t368 = t913 + t914;
t933 = -t558 * t559 - t561 * t828;
t1100 = m(6) * (t318 * t368 + t933);
t1097 = m(6) * (t472 * t518 + t473 * t519);
t1096 = m(6) * (t472 * t547 - t473 * t828);
t1095 = m(6) * (t472 * t744 + t473 * t743);
t1091 = m(6) * t420;
t791 = t1157 * t649;
t1085 = m(7) * (-t434 * t495 + t436 * t494 + t791);
t1084 = m(7) * (-t455 * t495 + t457 * t494 + t791);
t187 = t306 * t743 - t307 * t744;
t1082 = m(7) * t187;
t1081 = m(7) * (t1167 * t423 + t379 * t422);
t1080 = m(7) * (t1167 * t427 + t379 * t428);
t1079 = m(7) * t1157;
t1078 = m(7) * t1156;
t1077 = m(7) * t316;
t1070 = m(7) * t389;
t1061 = m(7) * qJD(6);
t1033 = t214 * t744;
t1032 = t215 * t743;
t1028 = t300 * t743;
t981 = t743 * t756;
t615 = Icges(4,4) * t981 - Icges(4,2) * t983 - Icges(4,6) * t744;
t1013 = t615 * t753;
t963 = t744 * t756;
t286 = m(7) * (-t456 * t744 + t458 * t743) + m(6) * (-t559 * t744 + t561 * t743);
t952 = t286 * qJD(4) + t187 * t1061 / 0.2e1;
t883 = -t1082 / 0.2e1;
t892 = t208 * qJD(2);
t951 = qJD(5) * t883 - t892;
t882 = t1082 / 0.2e1;
t950 = t890 * t882;
t613 = Icges(4,5) * t981 - Icges(4,6) * t983 - Icges(4,3) * t744;
t716 = Icges(4,4) * t983;
t617 = Icges(4,1) * t981 - Icges(4,5) * t744 - t716;
t920 = -t743 * t613 - t617 * t963;
t818 = Icges(4,5) * t756 - Icges(4,6) * t753;
t614 = Icges(4,3) * t743 + t744 * t818;
t618 = Icges(4,5) * t743 + t709 * t744;
t919 = t743 * t614 + t618 * t963;
t912 = t684 * t743 + t588;
t911 = -t684 * t744 - t589;
t910 = -Icges(5,2) * t985 + t590 - t701;
t909 = -t682 * t744 + t591;
t906 = t708 * t743 + t615;
t616 = Icges(4,6) * t743 + t707 * t744;
t905 = -t708 * t744 - t616;
t904 = -Icges(4,2) * t981 + t617 - t716;
t903 = -t706 * t744 + t618;
t132 = 0.2e1 * (t316 / 0.4e1 - t296 / 0.4e1) * m(7) + 0.2e1 * (t420 / 0.4e1 - t368 / 0.4e1) * m(6);
t894 = t132 * qJD(2);
t887 = t1127 / 0.2e1 + t1049;
t811 = t301 * t744 - t1028;
t137 = -t388 * t739 + t738 * t811;
t771 = t738 * t796 + t1023;
t183 = t506 * t640 - t508 * t641 + t743 * t771;
t770 = t738 * t795 + t1022;
t184 = t507 * t640 - t509 * t641 + t743 * t770;
t32 = (-t245 + t813) * t739 + (t183 * t743 + t184 * t744 + t345) * t738;
t185 = -t506 * t642 - t508 * t643 + t744 * t771;
t186 = -t507 * t642 - t509 * t643 + t744 * t770;
t33 = (-t246 + t812) * t739 + (t185 * t743 + t186 * t744 + t347) * t738;
t45 = (-t290 + t811) * t739 + (t214 * t743 + t215 * t744 + t388) * t738;
t14 = t1125 + (t978 / 0.2e1 + t993 / 0.2e1 - t45 / 0.2e1) * t739 + (t33 * t1115 + t32 * t1117 + t137 / 0.2e1) * t738;
t881 = qJD(5) * t882 + t14 * qJD(6) + t892;
t870 = t234 * t316 + t938;
t869 = t318 * t420 + t933;
t866 = t1006 / 0.2e1;
t863 = t1005 / 0.2e1;
t852 = -t687 - t749;
t555 = t618 * t981;
t842 = t614 * t744 - t555;
t839 = t616 * t753 - t613;
t835 = t895 * t1068;
t830 = -t749 - t1066;
t817 = -Icges(4,5) * t753 - Icges(4,6) * t756;
t800 = -t676 + t830;
t783 = -t746 * t909 + t747 * t911;
t784 = t746 * t910 + t747 * t912;
t80 = -t183 * t744 + t184 * t743;
t81 = -t185 * t744 + t186 * t743;
t799 = (t81 + (t784 * t744 + t1144 + (t783 - t1176) * t743) * t744 + t1175 * t741) * t1117 + (t80 + (t783 * t743 + t1144 + (t784 - t1175) * t744) * t743 + t1176 * t742) * t1116;
t792 = t830 + t1154;
t789 = t309 + t799;
t788 = -t835 + t914;
t785 = t1156 * t649 + t344 * t389;
t782 = t753 * t904 + t756 * t906;
t781 = -t753 * t903 + t756 * t905;
t778 = t1145 + t1217;
t776 = t32 * t1116 + t33 * t1117 + t80 * t866 + t81 * t863 + (t1032 - t1033) * t1119 + (t300 * t744 + t301 * t743) * t1191 - t1049 + t1179 * t739;
t772 = t1178 * t865 + t1223;
t69 = -t288 * t739 + (t199 * t743 + t200 * t744) * t738;
t70 = -t289 * t739 + (t201 * t743 + t202 * t744) * t738;
t768 = t137 * t1192 - t32 * t1006 / 0.2e1 - t33 * t1005 / 0.2e1 + t45 * t1190 - t1125 + t70 * t1117 + t69 * t1116 + t102 * t866 + t103 * t863 + (t993 + t978 + t1030 - t1031) * t1119;
t763 = t772 + t778 - t838;
t762 = t1178 * t867 - t1223 + t778 + t838;
t761 = -t1145 + t1217 + t772 + t838;
t760 = t1032 / 0.2e1 - t1033 / 0.2e1 + (t1174 * t744 + t1226 * t743 + t738 * t917 + t739 * t915 + t746 * t911 + t747 * t909 + t246) * t1117 + (t1174 * t743 - t738 * t918 + t739 * t916 - t746 * t912 + t747 * t910 - t1221 + t245) * t1116 - t767;
t759 = -t1028 / 0.2e1 - t765 + t1173 * t1118 + (t300 + t1173) * t1117;
t712 = -rSges(4,2) * t753 + t1060;
t659 = t817 * t744;
t658 = t817 * t743;
t604 = t852 * t744;
t602 = t852 * t743;
t517 = t800 * t744;
t515 = t800 * t743;
t496 = -t664 * t743 - t665 * t744;
t479 = m(5) * t480;
t440 = -t835 + t480;
t437 = t792 * t744;
t435 = t792 * t743;
t407 = -t1005 * t649 - t495 * t739;
t406 = t1006 * t649 + t494 * t739;
t400 = -t616 * t965 + t919;
t399 = -t615 * t965 - t920;
t398 = -t616 * t983 - t842;
t381 = -t1070 / 0.2e1;
t370 = (t494 * t744 - t495 * t743) * t738;
t354 = t788 + t913;
t327 = -t1000 + (-t752 * t907 + t755 * t908) * t738;
t303 = -t399 * t744 + t400 * t743;
t302 = -(-t743 * (-t617 * t756 + t1013) - t613 * t744) * t744 + t398 * t743;
t282 = -t1078 / 0.2e1;
t277 = t788 + t837;
t194 = -t1079 + t1095;
t178 = qJD(6) * t883;
t145 = t1084 / 0.2e1;
t142 = t1085 / 0.2e1;
t135 = t1151 * t649 + t188;
t134 = (t398 - t555 + (t614 + t1013) * t744 + t920) * t744 + t919 * t743;
t133 = (t744 * t839 + t400 - t919) * t744 + (t743 * t839 + t399 + t842) * t743;
t116 = t479 + t368 * t1135 + t296 * t1134 + t1091 / 0.2e1 + t1077 / 0.2e1;
t115 = t1152 * t649 + t167;
t107 = t940 + t944;
t88 = t282 + t1070 / 0.2e1;
t87 = t381 + t282;
t86 = t381 + t1078 / 0.2e1;
t84 = t942 + t947;
t71 = t1123 / 0.2e1;
t64 = t1124 / 0.2e1;
t61 = t1080 + t1096 + t765 + t1105;
t56 = t1170 + t765 + (t708 / 0.2e1 + t707 / 0.2e1) * t756 + t1081 + t1097 + t1106 + t1112;
t26 = t1130 / 0.2e1;
t24 = t1131 / 0.2e1;
t21 = m(7) * t135 + t1049;
t20 = m(7) * t115 + t1049;
t19 = t1100 + t789 + t1120;
t18 = t799 + t1224;
t16 = t26 - t1131 / 0.2e1 + t887;
t15 = t24 - t1130 / 0.2e1 + t887;
t11 = (t133 / 0.2e1 + t302 / 0.2e1) * t743 + (t303 / 0.2e1 - t134 / 0.2e1) * t744 + t767;
t10 = t24 + t26 - t1127 / 0.2e1 + t776;
t9 = t767 + t1172;
t8 = t767 - t1172;
t7 = t145 + t71 + t761;
t6 = t145 + t762 - t1123 / 0.2e1;
t5 = t71 + t763 - t1084 / 0.2e1;
t4 = t64 + t763 - t1085 / 0.2e1;
t3 = t64 + t761 + t142;
t2 = t762 - t1124 / 0.2e1 + t142;
t1 = t760 + t872 + t873;
t12 = [t56 * qJD(3) + t61 * qJD(4) + t194 * qJD(5) + t158 * qJD(6), 0, t56 * qJD(1) + t1 * qJD(4) + t84 * qJD(5) + t3 * qJD(6) + (((-t512 * t744 - t513 * t743) * t712 + (-t664 * t744 + t665 * t743) * t710) * t1137 + (t481 * t604 + t482 * t602) * t1136 + (t472 * t517 + t473 * t515 - t514 * t519 - t516 * t518) * t1135 + (t1167 * t437 + t379 * t435 - t422 * t434 - t423 * t436) * t1134) * t1141 + (t134 * t1115 + t760 + (t753 * t905 + t756 * t903) * t1117 + (t133 + t302) * t1118 + (-t753 * t906 + t756 * t904 + t303) * t1116 + (t742 / 0.2e1 + t741 / 0.2e1) * t818) * qJD(3), t61 * qJD(1) + t1 * qJD(3) + t760 * qJD(4) + t107 * qJD(5) + t7 * qJD(6) + ((t790 + (-t644 * t744 + t645 * t743) * t686) * t1136 + (t935 + (-t547 + t558) * t828) * t1135 + (-t427 * t457 - t428 * t455 + t946) * t1134) * t1139, qJD(1) * t194 + qJD(3) * t84 + qJD(4) * t107 + qJD(6) * t87, t1222 + t3 * qJD(3) + t7 * qJD(4) + t87 * qJD(5) + (-t327 * t739 + (t1167 * t406 - t1207 * t494 + t379 * t407 - t393 * t495) * m(7) + ((t289 / 0.2e1 + t230 / 0.2e1) * t744 + (t288 / 0.2e1 + t229 / 0.2e1) * t743) * t738) * qJD(6); 0, 0, t116 * qJD(4) + (t1134 * t277 + t1135 * t354 + t1136 * t440 + t1137 * t496) * t1141 + t1209, t116 * qJD(3) + (t1077 + t479 + t1091) * qJD(4) + t1209, 0, t1061 * t370 + t209 * t890; (t759 - t1170 - (t708 + t707) * t756 / 0.2e1) * qJD(1) + t11 * qJD(3) + t8 * qJD(4) - t83 * qJD(5) + t2 * qJD(6) + (-t1081 / 0.4e1 - t1097 / 0.4e1 - t1106 / 0.4e1 - t1112 / 0.4e1) * t1142, qJD(4) * t132 - t1208, t11 * qJD(1) + (m(7) * (t212 * t277 - t434 * t435 - t436 * t437) + m(6) * (t260 * t354 - t514 * t515 - t516 * t517) + m(5) * (-t1158 * t604 - t1159 * t602 + t350 * t440) + (t659 * t741 + (t782 * t744 + (-t658 + t781) * t743) * t744) * t1117 + (t658 * t742 + (t781 * t743 + (-t659 + t782) * t744) * t743) * t1116 + m(4) * (t710 * t712 * t895 + (t743 * (rSges(4,1) * t981 - t896) + t744 * (rSges(4,1) * t963 + t743 * rSges(4,3) - t718)) * t496) + t799) * qJD(3) + t18 * qJD(4) + t20 * qJD(6), t8 * qJD(1) + t894 + t18 * qJD(3) + t799 * qJD(4) + t15 * qJD(6) + (-t309 / 0.4e1 - t1120 / 0.4e1 - t1100 / 0.4e1) * t1138 + ((t870 + t82) * t1134 + (t869 + t166) * t1135 + (t236 + t310) * t1136) * t1139, t178 - t1215, t2 * qJD(1) + t20 * qJD(3) + t15 * qJD(4) + (t768 + m(7) * (t212 * t370 - t406 * t436 - t407 * t434 + t785)) * qJD(6) + t951; t759 * qJD(1) + t9 * qJD(3) + t767 * qJD(4) - t106 * qJD(5) + t6 * qJD(6) + (-t1080 / 0.4e1 - t1096 / 0.4e1 - t1105 / 0.4e1) * t1142, -qJD(3) * t132 - t1208, t9 * qJD(1) - t894 + t19 * qJD(4) + t16 * qJD(6) + ((t212 * t296 + t234 * t277 - t435 * t455 - t437 * t457 + t939) * t1134 + (t260 * t368 + t318 * t354 - t515 * t558 - t517 * t828 + t934) * t1135 + (t440 * t448 + (-t602 * t743 - t604 * t744) * t686 + t236) * t1136) * t1141 + (t799 - t1224) * qJD(3), t767 * qJD(1) + t19 * qJD(3) + t789 * qJD(4) + t21 * qJD(6) + (m(7) * t870 / 0.4e1 + m(6) * t869 / 0.4e1) * t1138, t178 - t1210, t6 * qJD(1) + t16 * qJD(3) + t21 * qJD(4) + (t768 + m(7) * (t234 * t370 - t406 * t457 - t407 * t455 + t785)) * qJD(6) + t951; t83 * qJD(3) + t106 * qJD(4) + t86 * qJD(6) + (t1079 / 0.4e1 - t1095 / 0.4e1) * t1142, 0, t1215 + ((-t435 * t744 + t437 * t743) * t1134 + (-t515 * t744 + t517 * t743) * t1135) * t1141 + t952, qJD(3) * t286 + t1210 + t952, 0, t86 * qJD(1) + (t406 * t743 - t407 * t744) * t1061 + t950; t4 * qJD(3) + t5 * qJD(4) + t88 * qJD(5) - t1222, t208 * t890, t4 * qJD(1) + ((t1207 * t437 + t277 * t344 - t393 * t435 - t115 + t871) * m(7) + t776) * qJD(3) + t10 * qJD(4) + t881, t5 * qJD(1) + t10 * qJD(3) + ((t224 - t135 + t824) * m(7) + t776) * qJD(4) + t881, qJD(1) * t88 + t950 (t739 ^ 2 * t327 / 0.2e1 + m(7) * (t1207 * t406 + t344 * t370 - t393 * t407) + (t70 * t1115 + t69 * t1117 + (t229 * t743 + t230 * t744) * t1119) * t738) * qJD(6) + t890 * t14;];
Cq  = t12;
