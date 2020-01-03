% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR15_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:14
% EndTime: 2019-12-31 20:42:19
% DurationCPUTime: 51.99s
% Computational Cost: add. (93913->1172), mult. (147608->1584), div. (0->0), fcn. (161275->8), ass. (0->660)
t1180 = m(6) / 0.2e1;
t1182 = m(5) / 0.2e1;
t1183 = m(4) / 0.2e1;
t1177 = -pkin(8) - pkin(7);
t783 = sin(qJ(2));
t784 = sin(qJ(1));
t1041 = t783 * t784;
t768 = pkin(2) * t1041;
t786 = cos(qJ(2));
t782 = sin(qJ(4));
t1100 = pkin(4) * t782;
t781 = qJ(4) + qJ(5);
t770 = sin(t781);
t771 = cos(t781);
t872 = rSges(6,1) * t770 + rSges(6,2) * t771;
t821 = -t872 - t1100;
t816 = qJ(3) - t821;
t488 = t768 + ((rSges(6,3) - t1177) * t783 - t816 * t786) * t784;
t787 = cos(qJ(1));
t1039 = t783 * t787;
t1156 = rSges(6,3) + pkin(2);
t1024 = t786 * t787;
t754 = qJ(3) * t1024;
t930 = t782 * t1024;
t931 = t787 * t1177;
t947 = pkin(4) * t930 + t783 * t931;
t951 = t872 * t1024;
t489 = -t1039 * t1156 + t754 + t947 + t951;
t1157 = rSges(5,3) + pkin(7);
t785 = cos(qJ(4));
t1094 = rSges(5,2) * t785;
t875 = rSges(5,1) * t782 + t1094;
t526 = t768 + (t1157 * t783 + (-qJ(3) - t875) * t786) * t784;
t937 = pkin(2) + t1157;
t950 = rSges(5,1) * t930 + t1024 * t1094;
t527 = -t1039 * t937 + t754 + t950;
t1086 = rSges(4,3) + qJ(3);
t1095 = rSges(4,2) * t783;
t597 = t768 + (-t1086 * t786 - t1095) * t784;
t882 = -pkin(2) * t1039 + t754;
t943 = rSges(4,2) * t1039 + rSges(4,3) * t1024;
t598 = t882 + t943;
t1030 = t784 * t786;
t777 = t787 * pkin(6);
t896 = rSges(4,1) * t787 - rSges(4,3) * t1041;
t1048 = t783 * qJ(3);
t899 = -pkin(1) - t1048;
t563 = t777 + ((rSges(4,2) - pkin(2)) * t786 + t899) * t784 + t896;
t1101 = pkin(2) * t786;
t767 = rSges(4,2) * t1024;
t564 = -t767 + (rSges(4,1) + pkin(6)) * t784 + (t1086 * t783 + pkin(1) + t1101) * t787;
t978 = t563 * t1024 + t564 * t1030;
t1031 = t784 * t785;
t695 = t1031 * t783 + t782 * t787;
t1027 = t785 * t787;
t1032 = t784 * t782;
t928 = t783 * t1032;
t696 = -t928 + t1027;
t1203 = t696 * rSges(5,1) - t695 * rSges(5,2);
t769 = pkin(7) * t1030;
t721 = pkin(3) * t787 - t769;
t457 = t777 + ((-rSges(5,3) - pkin(2)) * t786 + t899) * t784 + t721 + t1203;
t1178 = pkin(3) + pkin(6);
t693 = t783 * t1027 - t1032;
t694 = t1039 * t782 + t1031;
t877 = t694 * rSges(5,1) + t693 * rSges(5,2);
t459 = t1178 * t784 + (t786 * t937 - t899) * t787 + t877;
t990 = t457 * t1024 + t459 * t1030;
t1198 = -t1156 * t786 - pkin(1);
t1033 = t784 * t771;
t672 = t1033 * t783 + t770 * t787;
t1034 = t784 * t770;
t673 = -t1034 * t783 + t771 * t787;
t1204 = t673 * rSges(6,1) - t672 * rSges(6,2);
t1099 = pkin(4) * t785;
t946 = t787 * (pkin(3) + t1099) + t1177 * t1030;
t421 = t777 + ((-qJ(3) - t1100) * t783 + t1198) * t784 + t946 + t1204;
t753 = t786 * t931;
t873 = -rSges(6,1) * t771 + rSges(6,2) * t770;
t820 = t873 - t1099;
t423 = -t753 + (-t820 + t1178) * t784 + (t783 * t816 - t1198) * t787;
t994 = t421 * t1024 + t423 * t1030;
t922 = ((t488 * t787 + t489 * t784) * t783 + t994) * t1180 + ((t526 * t787 + t527 * t784) * t783 + t990) * t1182 + ((t597 * t787 + t598 * t784) * t783 + t978) * t1183;
t1026 = t786 * qJ(3);
t737 = pkin(2) * t783 - t1026;
t898 = t783 * pkin(7) + t737;
t1196 = -t783 * rSges(6,3) + t786 * t872;
t1098 = pkin(4) * t786;
t841 = -pkin(8) * t783 + t1098 * t782;
t962 = -t1196 - t841;
t844 = t898 + t962;
t484 = t844 * t784;
t486 = t844 * t787;
t1195 = -rSges(5,3) * t783 + t786 * t875;
t880 = -t1195 + t898;
t539 = t880 * t784;
t541 = t880 * t787;
t871 = rSges(4,3) * t786 + t1095;
t949 = t737 - t871;
t613 = t949 * t784;
t615 = t949 * t787;
t924 = (-t484 * t1039 + t1041 * t486 + t994) * t1180 + (-t539 * t1039 + t1041 * t541 + t990) * t1182 + (-t613 * t1039 + t1041 * t615 + t978) * t1183;
t27 = t924 - t922;
t1276 = t27 * qJD(1);
t1275 = Icges(4,1) + Icges(3,3);
t835 = t1039 * t770 + t1033;
t836 = t1039 * t771 - t1034;
t549 = t836 * rSges(6,1) - t835 * rSges(6,2);
t550 = t672 * rSges(6,1) + t673 * rSges(6,2);
t1274 = m(6) * (-t421 * t550 + t423 * t549);
t517 = -rSges(6,3) * t1030 + t1204;
t565 = pkin(4) * t928 + t721 - t946;
t1206 = -t517 + t565;
t595 = t1196 * t1030;
t969 = -t1030 * t841 - t595;
t377 = -t1206 * t783 + t969;
t506 = t783 * t517;
t997 = t565 * t783 + t377 - t506 - t969;
t1273 = m(6) * t997;
t908 = t1030 / 0.4e1;
t910 = -t1030 / 0.4e1;
t1272 = t910 + t908;
t443 = -t595 + t506;
t729 = Icges(3,5) * t786 - Icges(3,6) * t783;
t730 = -Icges(4,4) * t786 + Icges(4,5) * t783;
t1271 = (t729 + t730) * t787 + t1275 * t784;
t1260 = t1275 * t787 + (-Icges(4,5) + Icges(3,6)) * t1041 + (Icges(4,4) - Icges(3,5)) * t1030;
t861 = Icges(6,5) * t770 + Icges(6,6) * t771;
t805 = -Icges(6,3) * t783 + t786 * t861;
t1079 = Icges(6,4) * t770;
t864 = Icges(6,2) * t771 + t1079;
t807 = -Icges(6,6) * t783 + t786 * t864;
t1078 = Icges(6,4) * t771;
t867 = Icges(6,1) * t770 + t1078;
t809 = -Icges(6,5) * t783 + t786 * t867;
t405 = -t1024 * t805 - t807 * t836 - t809 * t835;
t508 = Icges(6,5) * t835 + Icges(6,6) * t836 + Icges(6,3) * t1024;
t1219 = Icges(6,4) * t835;
t511 = Icges(6,2) * t836 + Icges(6,6) * t1024 + t1219;
t643 = Icges(6,4) * t836;
t514 = Icges(6,1) * t835 + Icges(6,5) * t1024 + t643;
t310 = t508 * t1024 + t836 * t511 + t835 * t514;
t510 = -Icges(6,5) * t673 + Icges(6,6) * t672 + Icges(6,3) * t1030;
t645 = Icges(6,4) * t673;
t513 = Icges(6,2) * t672 + Icges(6,6) * t1030 - t645;
t644 = Icges(6,4) * t672;
t515 = Icges(6,1) * t673 - Icges(6,5) * t1030 - t644;
t311 = t510 * t1024 + t836 * t513 - t515 * t835;
t859 = t310 * t787 + t784 * t311;
t1267 = t405 * t783 + t786 * t859;
t911 = -t1030 / 0.2e1;
t1270 = t1267 * t911;
t862 = Icges(5,5) * t782 + Icges(5,6) * t785;
t806 = -Icges(5,3) * t783 + t786 * t862;
t1081 = Icges(5,4) * t782;
t865 = Icges(5,2) * t785 + t1081;
t808 = -Icges(5,6) * t783 + t786 * t865;
t1080 = Icges(5,4) * t785;
t868 = Icges(5,1) * t782 + t1080;
t810 = -Icges(5,5) * t783 + t786 * t868;
t426 = -t1030 * t806 - t695 * t808 + t696 * t810;
t551 = Icges(5,5) * t694 + Icges(5,6) * t693 + Icges(5,3) * t1024;
t1082 = Icges(5,4) * t694;
t554 = Icges(5,2) * t693 + Icges(5,6) * t1024 + t1082;
t686 = Icges(5,4) * t693;
t557 = Icges(5,1) * t694 + Icges(5,5) * t1024 + t686;
t345 = t551 * t1030 + t695 * t554 - t696 * t557;
t553 = -Icges(5,5) * t696 + Icges(5,6) * t695 + Icges(5,3) * t1030;
t688 = Icges(5,4) * t696;
t556 = Icges(5,2) * t695 + Icges(5,6) * t1030 - t688;
t687 = Icges(5,4) * t695;
t558 = Icges(5,1) * t696 - Icges(5,5) * t1030 - t687;
t346 = t1030 * t553 + t556 * t695 + t558 * t696;
t855 = t345 * t787 + t346 * t784;
t183 = t426 * t783 + t786 * t855;
t407 = -t1030 * t805 - t672 * t807 + t673 * t809;
t312 = t508 * t1030 + t672 * t511 - t673 * t514;
t313 = t1030 * t510 + t513 * t672 + t515 * t673;
t858 = t312 * t787 + t313 * t784;
t161 = t407 * t783 + t786 * t858;
t1061 = t553 * t783;
t1255 = t556 * t785 - t558 * t782;
t385 = t1255 * t786 - t1061;
t1063 = t510 * t783;
t1256 = t513 * t771 - t515 * t770;
t341 = t1256 * t786 - t1063;
t1264 = t310 * t784 - t311 * t787;
t1269 = t1272 * t1264;
t1083 = Icges(3,4) * t783;
t734 = Icges(3,1) * t786 - t1083;
t661 = Icges(3,5) * t784 + t734 * t787;
t1073 = Icges(4,6) * t786;
t726 = Icges(4,3) * t783 - t1073;
t662 = Icges(4,5) * t784 + t726 * t787;
t1268 = -t661 * t1030 - t662 * t1041;
t424 = -t1024 * t806 - t693 * t808 - t694 * t810;
t343 = t551 * t1024 + t693 * t554 + t694 * t557;
t344 = t553 * t1024 + t693 * t556 - t694 * t558;
t856 = t343 * t787 + t784 * t344;
t1266 = t424 * t783 + t786 * t856;
t561 = -rSges(5,3) * t1030 + t1203;
t476 = -t1195 * t1030 + t561 * t783;
t1265 = t343 * t784 - t344 * t787;
t1262 = -m(6) / 0.2e1;
t1254 = -t786 / 0.2e1;
t916 = -t1041 / 0.4e1;
t963 = -t809 + (Icges(6,2) * t770 - t1078) * t786;
t1261 = t771 * t963;
t778 = t784 ^ 2;
t780 = t787 ^ 2;
t942 = t778 + t780;
t1259 = t1271 * t787 + t1268;
t761 = Icges(3,4) * t1041;
t660 = Icges(3,1) * t1030 - Icges(3,5) * t787 - t761;
t663 = Icges(4,5) * t787 + Icges(4,6) * t1030 - Icges(4,3) * t1041;
t1258 = -t660 * t1024 + t663 * t1039 + t1260 * t784;
t1236 = t661 * t1024 + t662 * t1039 + t1271 * t784;
t611 = t1195 * t784;
t592 = t1196 * t784;
t656 = Icges(3,4) * t1030 - Icges(3,2) * t1041 - Icges(3,6) * t787;
t756 = Icges(4,6) * t1041;
t665 = Icges(4,4) * t787 + Icges(4,2) * t1030 - t756;
t1257 = t656 * t783 - t665 * t786;
t1163 = t783 / 0.2e1;
t1162 = t784 / 0.2e1;
t1161 = t784 / 0.4e1;
t1160 = t786 / 0.2e1;
t1159 = -t787 / 0.2e1;
t1158 = -t787 / 0.4e1;
t689 = t693 * pkin(4);
t977 = -t549 - t689;
t1244 = t784 * t977;
t853 = -t511 * t771 - t514 * t770;
t832 = t805 * t787 - t853;
t589 = t807 * t787;
t591 = t809 * t787;
t834 = t589 * t771 + t591 * t770 - t508;
t250 = t783 * t832 - t786 * t834;
t632 = Icges(6,6) * t786 + t783 * t864;
t634 = Icges(6,5) * t786 + t783 * t867;
t630 = Icges(6,3) * t786 + t783 * t861;
t1052 = t771 * t807;
t1055 = t770 * t809;
t848 = t1052 + t1055;
t828 = t630 - t848;
t817 = t828 * t786;
t1053 = t771 * t632;
t1056 = t770 * t634;
t827 = t1053 + t805 + t1056;
t302 = (-t632 * t770 + t634 * t771) * t784 + (t783 * t827 + t817) * t787;
t1243 = t250 + t302;
t831 = t805 * t784 + t1256;
t588 = t807 * t784;
t590 = t809 * t784;
t833 = t588 * t771 + t590 * t770 - t510;
t251 = t783 * t831 - t786 * t833;
t1060 = t805 * t783;
t301 = t632 * t672 - t634 * t673 + (t817 + t1060) * t784;
t1242 = t251 + t301;
t1064 = t508 * t783;
t340 = t786 * t853 + t1064;
t1241 = t340 + t405;
t1240 = -t341 + t407;
t1075 = Icges(3,2) * t783;
t775 = Icges(3,4) * t786;
t657 = Icges(3,6) * t784 + (t775 - t1075) * t787;
t757 = Icges(4,6) * t1039;
t664 = Icges(4,4) * t784 - Icges(4,2) * t1024 + t757;
t1239 = -t1024 * t664 - t1039 * t657 + t1236;
t690 = t695 * pkin(4);
t493 = -t550 - t690;
t1238 = t783 * t657 + t664 * t786 + t1260;
t793 = rSges(6,1) * t835 + rSges(6,2) * t836 + rSges(6,3) * t1024;
t505 = t783 * t793;
t935 = pkin(7) * t1024;
t799 = pkin(4) * t694 - t753 - t935;
t795 = t783 * t799;
t379 = t1024 * t962 - t505 - t795;
t370 = t379 * t1039;
t815 = rSges(5,3) * t1024 + t877;
t477 = -t1024 * t1195 - t783 * t815;
t1002 = (-t1041 * t377 - t370) * t1180 + (-t477 * t1039 - t1041 * t476) * t1182;
t1013 = (-t370 + (t379 * t787 - t784 * t997) * t783) * t1180;
t1231 = t1002 - t1013;
t1207 = t942 * t783;
t1018 = t787 * t565;
t504 = t517 * t1024;
t307 = -t504 + (t1018 + (t753 + t820 * t784 + ((-rSges(6,3) + pkin(7)) * t786 + t821 * t783) * t787) * t784) * t786;
t1197 = t561 * t787 + t784 * t815;
t430 = t1197 * t786;
t989 = t476 * t1024 - t477 * t1030;
t998 = t377 * t1024 - t379 * t1030;
t1011 = (t1207 * t307 + t998) * t1180 + (-t1207 * t430 + t989) * t1182;
t593 = -rSges(6,3) * t1039 + t951;
t322 = t592 * t1024 - t1030 * t593 + t1039 * t517 + t784 * t505;
t623 = t841 * t784;
t624 = pkin(7) * t1039 + t947;
t222 = -t1018 * t783 + t1024 * t623 - t1030 * t624 + t784 * t795 + t322;
t637 = rSges(6,3) * t786 + t783 * t872;
t594 = t637 * t1030;
t685 = pkin(8) * t786 + t1100 * t783;
t970 = t592 + t623;
t272 = t594 + (t685 * t784 - t1206) * t786 + (-t784 * t962 - t970) * t783;
t507 = t786 * t793;
t381 = -t1024 * t637 - t1039 * t1196 + t783 * t593 + t507;
t273 = -t1024 * t685 - t1039 * t841 + t783 * t624 + t786 * t799 + t381;
t612 = -rSges(5,3) * t1039 + t950;
t350 = (t611 * t787 - t612 * t784) * t786 + t1197 * t783;
t1091 = rSges(5,3) * t786;
t668 = t783 * t875 + t1091;
t395 = (t668 * t784 + t561) * t786;
t396 = ((-t668 + t1091) * t787 + t877) * t786 + (-t1195 * t787 + t612) * t783;
t1085 = (-t350 * t786 + (t395 * t787 + t396 * t784 - t430) * t783 + t989) * t1182 + (-t222 * t786 + (t272 * t787 + t273 * t784 + t307) * t783 + t998) * t1180;
t1230 = t1011 - t1085;
t573 = rSges(5,1) * t693 - rSges(5,2) * t694;
t574 = rSges(5,1) * t695 + rSges(5,2) * t696;
t677 = t873 * t786;
t881 = t1098 * t785 - t677;
t584 = t881 * t784;
t585 = t881 * t787;
t718 = (-rSges(5,1) * t785 + rSges(5,2) * t782) * t786;
t1229 = (t421 * t585 + t423 * t584 + t484 * t977 - t486 * t493) * t1180 + (-t539 * t573 + t541 * t574 + (-t457 * t787 - t459 * t784) * t718) * t1182;
t1228 = -(t395 * t457 + t396 * t459 + t476 * t526 - t477 * t527) * m(5) / 0.2e1 + (t272 * t421 + t273 * t423 + t377 * t488 - t379 * t489) * t1262;
t1227 = t1024 * t665 - t1030 * t664 - t1039 * t656 - t1041 * t657 - t1258 - t1259;
t1226 = t484 * t997 * t1262;
t674 = (-Icges(6,5) * t771 + Icges(6,6) * t770) * t786;
t1043 = t783 * t674;
t1225 = t1043 / 0.2e1 + t1261 * t1254;
t1223 = t787 / 0.2e1;
t1220 = m(6) * t783;
t522 = t550 * t1024;
t428 = -t1030 * t549 + t522;
t626 = t677 * t1030;
t473 = -t550 * t783 + t626;
t531 = t783 * t549;
t474 = -t1024 * t677 + t531;
t543 = Icges(6,5) * t836 - Icges(6,6) * t835;
t983 = -Icges(6,2) * t835 + t514 + t643;
t985 = -Icges(6,1) * t836 + t1219 + t511;
t235 = t1030 * t543 + t672 * t983 + t673 * t985;
t544 = Icges(6,5) * t672 + Icges(6,6) * t673;
t982 = Icges(6,2) * t673 - t515 + t644;
t984 = -Icges(6,1) * t672 + t513 - t645;
t236 = t1030 * t544 + t672 * t982 + t673 * t984;
t676 = (-Icges(6,1) * t771 + t1079) * t786;
t964 = -t807 - t676;
t319 = t1030 * t674 + t672 * t963 + t673 * t964;
t100 = t319 * t783 + (t235 * t787 + t236 * t784) * t786;
t268 = t543 * t783 + (t770 * t985 - t771 * t983) * t786;
t269 = t544 * t783 + (t770 * t984 - t771 * t982) * t786;
t388 = (t1043 + (t770 * t964 - t1261) * t786) * t783;
t905 = t1024 / 0.2e1;
t909 = t1030 / 0.2e1;
t233 = t543 * t1024 - t835 * t985 + t836 * t983;
t234 = t544 * t1024 - t835 * t984 + t836 * t982;
t318 = t674 * t1024 - t835 * t964 + t836 * t963;
t99 = t318 * t783 + (t233 * t787 + t234 * t784) * t786;
t934 = t100 * t909 + t99 * t905 + (t388 + (t268 * t787 + t269 * t784) * t786) * t1163;
t22 = t934 + m(6) * (t307 * t428 + t377 * t473 - t379 * t474);
t1218 = t22 * qJD(5);
t1028 = t785 * t808;
t1049 = t782 * t810;
t650 = Icges(5,3) * t786 + t783 * t862;
t731 = Icges(3,2) * t786 + t1083;
t1212 = t783 * (t1028 / 0.2e1 + t1055 / 0.2e1 + t1052 / 0.2e1 + t1049 / 0.2e1 - t734 / 0.2e1 + t731 / 0.2e1 + Icges(4,2) * t1254 + Icges(4,6) * t783 + Icges(4,3) * t1160 - t650 / 0.2e1 - t630 / 0.2e1);
t1040 = t783 * t786;
t1179 = m(6) / 0.4e1;
t1181 = m(5) / 0.4e1;
t945 = t942 * t1040;
t1208 = (m(4) / 0.4e1 + t1181 + t1179) * (t945 - t1040);
t1205 = -t637 - t685;
t863 = -Icges(3,5) * t783 - Icges(3,6) * t786;
t866 = Icges(4,4) * t783 + Icges(4,5) * t786;
t1202 = (t863 + t866) * t787;
t326 = t783 * t828 - t786 * t827;
t435 = t786 * t848 - t1060;
t857 = t340 * t787 - t341 * t784;
t1200 = (t435 * t783 + t786 * t857) * t1160 + ((t250 * t787 + t251 * t784 + t435) * t786 + (t326 - t857) * t783) * t1163;
t134 = t233 * t784 - t234 * t787;
t135 = t235 * t784 - t236 * t787;
t1012 = t135 * t1159 + t134 * t1162;
t1199 = t435 * t1160 + t326 * t1163;
t849 = t573 * t784 - t574 * t787;
t995 = (t493 * t787 - t1244) * t1220 / 0.2e1 + m(5) * t849 * t1163;
t434 = t787 * t549 + t784 * t550;
t392 = t689 * t787 + t784 * t690 + t434;
t453 = t573 * t787 + t784 * t574;
t1001 = (-t786 * t392 + (t584 * t784 + t585 * t787) * t783) * t1180 + (-t1207 * t718 - t453 * t786) * t1182;
t704 = (Icges(5,2) * t782 - t1080) * t786;
t709 = (-Icges(5,1) * t785 + t1081) * t786;
t1194 = -t782 * (t709 / 0.2e1 + t808 / 0.2e1) - t785 * (-t810 / 0.2e1 + t704 / 0.2e1);
t567 = Icges(5,5) * t693 - Icges(5,6) * t694;
t974 = -Icges(5,2) * t694 + t557 + t686;
t976 = -Icges(5,1) * t693 + t1082 + t554;
t299 = t567 * t783 + (t782 * t976 - t785 * t974) * t786;
t568 = Icges(5,5) * t695 + Icges(5,6) * t696;
t973 = Icges(5,2) * t696 - t558 + t687;
t975 = -Icges(5,1) * t695 + t556 - t688;
t300 = t568 * t783 + (t782 * t975 - t785 * t973) * t786;
t701 = (-Icges(5,5) * t785 + Icges(5,6) * t782) * t786;
t958 = -t810 + t704;
t960 = -t808 - t709;
t363 = t1024 * t701 + t693 * t958 - t694 * t960;
t364 = t1030 * t701 + t695 * t958 + t696 * t960;
t1193 = t1229 + (t299 + t363) * t1161 + (t300 + t364) * t1158;
t1067 = t269 * t787;
t1068 = t268 * t784;
t885 = t1068 / 0.4e1 - t1067 / 0.4e1 + t318 * t1161 + t319 * t1158;
t1022 = t787 * t1267;
t1037 = t784 * t161;
t1192 = t1022 / 0.4e1 + t1037 / 0.4e1 + t1267 * t1158 - t161 * t1161;
t1021 = t787 * t1266;
t1036 = t784 * t183;
t1191 = t1036 / 0.4e1 + t1021 / 0.4e1 - t183 * t1161 + t1266 * t1158 - t1226 + t1272 * t1265;
t608 = t808 * t787;
t610 = t810 * t787;
t851 = -t554 * t785 - t557 * t782;
t830 = t806 * t787 - t851;
t287 = (-t608 * t785 - t610 * t782 + t551) * t786 + t830 * t783;
t607 = t808 * t784;
t609 = t810 * t784;
t829 = t806 * t784 + t1255;
t288 = (-t607 * t785 - t609 * t782 + t553) * t786 + t829 * t783;
t654 = Icges(5,6) * t786 + t783 * t865;
t658 = Icges(5,5) * t786 + t783 * t868;
t1059 = t806 * t783;
t847 = t1028 + t1049;
t826 = t650 - t847;
t800 = t786 * t826 + t1059;
t327 = t654 * t695 - t658 * t696 + t784 * t800;
t328 = t693 * t654 + t694 * t658 + t787 * t800;
t1029 = t785 * t654;
t1050 = t782 * t658;
t376 = (-t1029 - t806 - t1050) * t786 + t826 * t783;
t1062 = t551 * t783;
t384 = t786 * t851 + t1062;
t462 = t786 * t847 - t1059;
t904 = t1024 / 0.4e1;
t913 = -t1039 / 0.4e1;
t1190 = t462 * t1160 + t376 * t1163 - t1228 + (-t385 + t426) * t916 + (t384 + t424) * t913 + (t288 + t327) * t908 + (t287 + t328) * t904;
t1188 = 0.4e1 * qJD(1);
t1187 = 2 * qJD(2);
t1186 = 4 * qJD(2);
t1185 = 2 * qJD(4);
t1184 = 4 * qJD(4);
t380 = t517 * t786 + t594;
t410 = -t507 * t784 - t504;
t444 = -t1024 * t1196 - t505;
t1175 = m(6) * (t222 * t410 + t272 * t443 - t273 * t444 + t307 * t322 + t377 * t380 - t379 * t381);
t1173 = t444 * t1273;
t1171 = m(6) * (t222 * t307 + t272 * t377 - t273 * t379);
t741 = t1048 + t1101;
t953 = t942 * t741;
t884 = -t784 * t721 + t787 * (t784 * pkin(3) + t935) + t953;
t281 = t884 + (t793 + t799) * t787 + t1206 * t784;
t923 = t428 * t281 - t473 * t486 - t474 * t484;
t1169 = m(6) * (t307 * t434 + (-t377 * t787 + t379 * t784) * t677 + t923);
t1167 = m(6) * (t392 * t410 + t443 * t585 - t444 * t584 + t923);
t1166 = t379 * t1273;
t1165 = t183 / 0.2e1;
t1164 = t299 / 0.2e1;
t1097 = rSges(3,1) * t786;
t900 = pkin(1) + t1097;
t944 = rSges(3,2) * t1041 + t787 * rSges(3,3);
t603 = -t784 * t900 + t777 + t944;
t765 = rSges(3,2) * t1039;
t604 = -t765 + t900 * t787 + (rSges(3,3) + pkin(6)) * t784;
t739 = rSges(3,1) * t783 + rSges(3,2) * t786;
t716 = t739 * t784;
t719 = t739 * t787;
t1155 = m(3) * (t603 * t716 - t604 * t719);
t447 = -t784 * (rSges(4,2) * t1030 + t896) + t787 * (t784 * rSges(4,1) + rSges(4,3) * t1039 - t767) + t953;
t971 = -t615 * t1024 - t613 * t1030;
t1151 = m(4) * (t1207 * t447 + t971);
t1149 = m(4) * (t563 * t597 + t564 * t598);
t1148 = m(4) * (t564 * t1039 - t1041 * t563);
t1146 = m(5) * (-t350 * t430 + t395 * t476 - t396 * t477);
t356 = -t561 * t784 + t787 * t815 + t884;
t1141 = m(5) * (t356 * t453 + (t539 * t784 + t541 * t787) * t718);
t979 = -t541 * t1024 - t539 * t1030;
t1138 = m(5) * (t1207 * t356 + t979);
t1134 = m(5) * (t457 * t526 + t459 * t527);
t1133 = m(5) * (-t457 * t574 + t459 * t573);
t1132 = m(5) * (t459 * t1039 - t1041 * t457);
t992 = t443 * t1024 - t444 * t1030;
t1126 = m(6) * (-t322 * t786 + (t380 * t787 + t381 * t784 + t410) * t783 + t992);
t1125 = m(6) * (t380 * t421 + t381 * t423 + t443 * t488 - t444 * t489);
t1000 = t473 * t421 + t474 * t423;
t1123 = m(6) * (-t377 * t550 - t379 * t549 + t1000);
t1122 = m(6) * (t443 * t493 + t444 * t977 + t1000);
t1121 = m(6) * (t281 * t392 - t484 * t584 - t486 * t585);
t1117 = m(6) * (-t484 * t549 + t486 * t550 + (-t421 * t787 - t423 * t784) * t677);
t987 = -t486 * t1024 - t484 * t1030;
t1116 = m(6) * (t1207 * t281 + t987);
t1112 = m(6) * (t421 * t488 + t423 * t489);
t1111 = m(6) * (t421 * t493 - t423 * t977);
t1109 = m(6) * (t1207 * t410 + t992);
t1108 = m(6) * (-t428 * t786 + (t473 * t787 + t474 * t784) * t783);
t1106 = m(6) * (t423 * t1039 - t1041 * t421);
t1105 = m(6) * (-t444 * t1039 - t1041 * t443);
t1104 = m(6) * (-t1207 * t677 - t434 * t786);
t1102 = (t549 * t784 - t550 * t787) * t1220;
t1084 = Icges(3,1) * t783;
t1054 = t770 * t786;
t1042 = t783 * t701;
t869 = -t775 - t1084;
t710 = t869 * t784;
t959 = t656 - t710;
t707 = -Icges(3,2) * t1030 - t761;
t957 = t660 + t707;
t860 = Icges(4,2) * t783 + t1073;
t699 = t860 * t784;
t956 = t663 + t699;
t697 = Icges(4,3) * t1030 + t756;
t955 = t665 - t697;
t954 = t784 * (t1026 * t784 - t768) + t787 * t882;
t948 = rSges(4,2) * t786 - rSges(4,3) * t783 - t741;
t941 = qJD(1) * t786;
t940 = qJD(2) * t783;
t939 = qJD(2) * t786;
t938 = qJD(4) * t786;
t936 = t786 ^ 2 * t1099;
t926 = -t183 / 0.2e1 + t1165;
t921 = -t1054 / 0.2e1;
t920 = t1054 / 0.2e1;
t917 = -t1041 / 0.2e1;
t914 = -t1039 / 0.2e1;
t901 = t730 / 0.2e1 + t729 / 0.2e1;
t897 = -t786 * pkin(7) - t741;
t711 = t869 * t787;
t895 = (-t657 + t711) * t784;
t708 = t731 * t787;
t894 = (-t661 + t708) * t784;
t700 = t860 * t787;
t893 = (t662 - t700) * t784;
t698 = Icges(4,3) * t1024 + t757;
t892 = (t664 + t698) * t784;
t878 = t1267 * t909 + t1270;
t870 = t676 * t921 - t807 * t920 + t1225;
t854 = t384 * t787 - t385 * t784;
t819 = t832 * t786;
t228 = t589 * t672 - t591 * t673 + (t819 - t1064) * t784;
t818 = t831 * t786;
t229 = t588 * t672 - t590 * t673 + (t818 - t1063) * t784;
t46 = (t228 * t787 + t229 * t784 + t407) * t786 + (t301 - t858) * t783;
t230 = (-t589 * t770 + t591 * t771) * t784 + (t783 * t834 + t819) * t787;
t231 = (-t588 * t770 + t590 * t771) * t784 + (t783 * t833 + t818) * t787;
t47 = (t230 * t787 + t231 * t784 + t405) * t786 + (t302 - t859) * t783;
t843 = t1267 * t914 + t161 * t917 + t46 * t909 + t47 * t905 + t1200;
t842 = t1173 / 0.2e1 + t878;
t276 = t1024 * t567 + t693 * t974 - t694 * t976;
t277 = t1024 * t568 + t693 * t973 - t694 * t975;
t157 = t276 * t784 - t277 * t787;
t278 = t1030 * t567 + t695 * t974 + t696 * t976;
t279 = t1030 * t568 + t695 * t973 + t696 * t975;
t158 = t278 * t784 - t279 * t787;
t825 = t1159 * t158 + t1162 * t157;
t823 = t1175 / 0.2e1 + t843;
t822 = -pkin(7) * t1207 + t954;
t814 = t1192 + t1269;
t813 = t676 * t920 - t807 * t921 - t1225;
t128 = t228 * t784 - t229 * t787;
t129 = t230 * t784 - t231 * t787;
t812 = t46 * t1159 + t47 * t1162 + t128 * t909 + t129 * t905 + (t250 * t784 - t251 * t787) * t1163 + (t312 * t784 - t313 * t787) * t917 + t1264 * t914 + (t340 * t784 + t341 * t787) * t1160 - t1012;
t811 = t1270 + t388 + (t268 + t318) * t905 + (t1267 + t269 + t319) * t909;
t804 = t1240 * t916 + t1241 * t913 + t1242 * t908 + t1243 * t904 + t1199;
t803 = t100 * t1159 + t99 * t1162 + t134 * t905 + t135 * t909 + t46 * t911 - t47 * t1024 / 0.2e1 - t1200 + (-t1067 + t1068 + t1037 + t1022) * t1163;
t802 = t786 * t830 - t1062;
t801 = t786 * t829 - t1061;
t798 = -t1236 * t784 / 0.2e1 + t1239 * t1162 + ((t1257 + t1271) * t787 + t1227 + t1258 + t1268) * t1159;
t797 = (t1238 * t787 - t1236 + t1239) * t1223 + (t1260 * t787 + (t660 * t786 - t663 * t783 - t1257) * t784) * t1159 + (t1238 * t784 + t1227 + t1259) * t1162;
t792 = t1050 / 0.2e1 + t1029 / 0.2e1 + t1056 / 0.2e1 + t1053 / 0.2e1 - t775 - t1084 / 0.2e1 + t1075 / 0.2e1 - t860 / 0.2e1 + t726 / 0.2e1 + t806 / 0.2e1 + t805 / 0.2e1;
t790 = t804 + t814 - t885;
t789 = -t1199 + t814 + t885 + t1240 * t1041 / 0.4e1 + t1241 * t1039 / 0.4e1 + t1242 * t910 - t1243 * t1024 / 0.4e1;
t788 = -t1192 + t1269 + t804 + t885;
t743 = -rSges(3,2) * t783 + t1097;
t705 = t866 * t784;
t702 = t863 * t784;
t616 = t948 * t787;
t614 = t948 * t784;
t542 = (-t668 + t897) * t787;
t540 = -t769 + (-t668 - t741) * t784;
t491 = -t1024 * t718 + t783 * t573;
t490 = t1030 * t718 - t574 * t783;
t487 = (t897 + t1205) * t787;
t485 = -t769 + (-t741 + t1205) * t784;
t472 = t778 * t871 + t787 * t943 + t954;
t446 = t849 * t786;
t441 = 0.4e1 * t1208;
t419 = t783 * t689 + t531 + (-t677 * t786 + t936) * t787;
t418 = t493 * t783 - t784 * t936 + t626;
t417 = t1102 / 0.2e1;
t411 = t784 * t611 + t612 * t787 + t822;
t408 = (t1042 + (t782 * t960 - t785 * t958) * t786) * t783;
t383 = t522 + (t690 * t787 + t1244) * t786;
t373 = t1104 / 0.2e1;
t336 = (t593 + t624) * t787 + t970 * t784 + t822;
t325 = t1105 / 0.2e1;
t274 = qJD(5) * t1108;
t259 = t693 * t607 + t694 * t609 + t787 * t801;
t258 = t693 * t608 + t694 * t610 + t787 * t802;
t257 = t607 * t695 - t609 * t696 + t784 * t801;
t256 = t608 * t695 - t610 * t696 + t784 * t802;
t239 = t1109 / 0.2e1;
t202 = t462 * t783 + t786 * t854;
t184 = t870 + t1274;
t171 = t1117 / 0.2e1;
t167 = t1106 + t1132 + t1148;
t156 = t281 * t434 + (t484 * t784 + t486 * t787) * t677;
t142 = t258 * t784 - t259 * t787;
t141 = t256 * t784 - t257 * t787;
t136 = t1122 / 0.2e1;
t117 = t364 * t783 + (t278 * t787 + t279 * t784) * t786;
t116 = t363 * t783 + (t276 * t787 + t277 * t784) * t786;
t114 = t1123 / 0.2e1;
t112 = t322 * t410 + t380 * t443 - t381 * t444;
t111 = t325 - t1102 / 0.2e1;
t110 = t417 + t325;
t109 = t417 - t1105 / 0.2e1;
t105 = t1125 / 0.2e1;
t104 = t1126 / 0.2e1;
t103 = t1042 / 0.2e1 + t1133 + t1111 + t1194 * t786 + t870;
t94 = t1116 + t1138 + t1151;
t86 = (t287 * t787 + t288 * t784 + t462) * t786 + (t376 - t854) * t783;
t84 = t1167 / 0.2e1;
t77 = (t258 * t787 + t259 * t784 + t424) * t786 + (t328 - t856) * t783;
t76 = (t256 * t787 + t257 * t784 + t426) * t786 + (t327 - t855) * t783;
t66 = t239 + t104 - t1104 / 0.2e1;
t65 = t373 + t239 - t1126 / 0.2e1;
t64 = t373 + t104 - t1109 / 0.2e1;
t63 = -t786 * t792 + t1112 + t1134 + t1149 + t1155 - t1212;
t61 = t1169 / 0.2e1;
t33 = m(6) * t156 + t1012;
t32 = t1002 + t1013 - t995;
t31 = t995 - t1231;
t30 = t995 + t1231;
t26 = t922 + t924;
t24 = m(6) * (t410 * t428 + t443 * t473 - t444 * t474) + t934;
t23 = t24 * qJD(5);
t21 = t1011 + t1085 - t1001;
t20 = t1001 - t1230;
t19 = t1001 + t1230;
t18 = t1012 + t1121 + t825 + t1141;
t16 = m(6) * t112 + t843;
t15 = t136 - t1123 / 0.2e1 + t842;
t14 = t114 - t1122 / 0.2e1 + t842;
t13 = t114 + t136 - t1173 / 0.2e1 + t811;
t12 = t926 * t1024 + t1166 + t878;
t11 = t84 - t1169 / 0.2e1 + t823;
t10 = t61 - t1167 / 0.2e1 + t823;
t9 = t784 * t797 + t787 * t798;
t8 = t1146 + t1171 + (t202 / 0.2e1 + t77 * t1223 + t76 * t1162) * t786 + (t86 / 0.2e1 - t1021 / 0.2e1 - t1036 / 0.2e1) * t783 + t843;
t7 = t105 + t788 + t171;
t6 = -t1125 / 0.2e1 + t789 + t171;
t5 = t105 - t1117 / 0.2e1 + t790;
t4 = t84 + t803 + t61 - t1175 / 0.2e1;
t3 = t1190 + t788 + t1193 + t1226;
t2 = t1190 + t1191 + t790 + (-t363 / 0.4e1 - t299 / 0.4e1) * t784 + (t364 / 0.4e1 + t300 / 0.4e1) * t787 - t1229;
t1 = t1191 + (-t376 / 0.2e1 + (t384 / 0.4e1 + t424 / 0.4e1) * t787 + (t426 / 0.4e1 - t385 / 0.4e1) * t784) * t783 + (-t462 / 0.2e1 + (-t328 / 0.4e1 - t287 / 0.4e1) * t787 + (-t327 / 0.4e1 - t288 / 0.4e1) * t784) * t786 + t789 + t1193 + t1228;
t17 = [t63 * qJD(2) + t167 * qJD(3) + t103 * qJD(4) + t184 * qJD(5), t63 * qJD(1) + t26 * qJD(3) + t3 * qJD(4) + t7 * qJD(5) + ((t563 * t616 + t564 * t614 - t597 * t615 - t598 * t613) * t1183 + (t457 * t542 + t459 * t540 - t526 * t541 - t527 * t539) * t1182 + (t421 * t487 + t423 * t485 - t484 * t489 - t486 * t488) * t1180) * t1187 + ((m(3) * (-t603 * t743 - t716 * t739) - t288 / 0.2e1 - t301 / 0.2e1 - t327 / 0.2e1 - t251 / 0.2e1 + t901 * t787 - t798) * qJD(2) + (-t660 / 0.2e1 - t707 / 0.2e1 - t665 / 0.2e1 + t697 / 0.2e1) * t939 + (t656 / 0.2e1 - t710 / 0.2e1 + t663 / 0.2e1 + t699 / 0.2e1) * t940) * t787 + ((t328 / 0.2e1 + t287 / 0.2e1 + t250 / 0.2e1 + m(3) * (-t604 * t743 + t719 * t739) + t302 / 0.2e1 + t901 * t784 - t797) * qJD(2) + (-t664 / 0.2e1 - t698 / 0.2e1 + t661 / 0.2e1 - t708 / 0.2e1) * t939 + (t662 / 0.2e1 - t700 / 0.2e1 + t711 / 0.2e1 - t657 / 0.2e1) * t940) * t784, qJD(1) * t167 + qJD(2) * t26 + qJD(4) * t30 + qJD(5) * t110, t103 * qJD(1) + t3 * qJD(2) + t30 * qJD(3) + (t408 + t811) * qJD(4) + t13 * qJD(5) - t1166 * t1184 / 0.4e1 + ((t377 * t493 + t379 * t977 + t418 * t421 + t419 * t423) * t1180 + (t457 * t490 + t459 * t491 - t476 * t574 - t477 * t573) * t1182) * t1185 + ((t363 / 0.2e1 + t1164 - t926) * t787 + (t300 / 0.2e1 + t364 / 0.2e1) * t784) * t938, t184 * qJD(1) + t7 * qJD(2) + t110 * qJD(3) + t13 * qJD(4) + ((-t443 * t550 - t444 * t549 + t1000) * m(6) + t811) * qJD(5); t9 * qJD(2) + t27 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + (-t1155 / 0.4e1 - t1149 / 0.4e1 - t1134 / 0.4e1 - t1112 / 0.4e1) * t1188 + t792 * t941 + t1212 * qJD(1), t9 * qJD(1) + t94 * qJD(3) + t18 * qJD(4) + t33 * qJD(5) + (m(6) * (t281 * t336 - t484 * t485 - t486 * t487) + m(5) * (t356 * t411 - t539 * t540 - t541 * t542) + m(4) * (t447 * t472 - t613 * t614 - t615 * t616) + m(3) * ((t784 * (rSges(3,1) * t1030 - t944) + t787 * (rSges(3,1) * t1024 + t784 * rSges(3,3) - t765)) * (-t784 * t716 - t719 * t787) + t942 * t743 * t739) + (t129 + t142 + (-t784 * t705 + (t787 * t956 + t893) * t786 + (t787 * t955 + t892) * t783) * t787 + (-t784 * t702 + (t787 * t959 + t895) * t786 + (t787 * t957 + t894) * t783) * t787 + t1202 * t778) * t1162 + (t128 + t141 + (-t1202 * t787 + (t893 + t895 + (t956 + t959) * t787) * t786 + (t892 + t894 + (t955 + t957) * t787) * t783) * t784 + (t705 + t702) * t780) * t1159) * qJD(2), t1276 + t94 * qJD(2) + t19 * qJD(4) + t65 * qJD(5) + (-0.4e1 * t1208 + 0.2e1 * (t1180 + t1182 + t1183) * (-t1207 * t786 + t945)) * qJD(3), t1 * qJD(1) + t18 * qJD(2) + t19 * qJD(3) + t4 * qJD(5) + (-t1146 / 0.4e1 - t1171 / 0.4e1) * t1184 + ((-t446 * t356 - t430 * t453 - t490 * t541 - t491 * t539 + (-t476 * t787 + t477 * t784) * t718) * t1182 + (t281 * t383 + t307 * t392 + t377 * t585 - t379 * t584 - t418 * t486 - t419 * t484) * t1180) * t1185 + (-t202 / 0.2e1 + (-t77 / 0.2e1 + t157 / 0.2e1) * t787 + (t158 / 0.2e1 - t76 / 0.2e1) * t784) * t938 + (t117 * t1159 + t116 * t1162 + t803 + (-t86 / 0.2e1 + (t1266 / 0.2e1 - t300 / 0.2e1) * t787 + (t1165 + t1164) * t784) * t783) * qJD(4), t6 * qJD(1) + t33 * qJD(2) + t65 * qJD(3) + t4 * qJD(4) + ((t410 * t434 + (-t443 * t787 + t444 * t784) * t677 - t112 + t923) * m(6) + t803) * qJD(5); -t27 * qJD(2) + t31 * qJD(4) + t109 * qJD(5) + (-t1106 / 0.4e1 - t1132 / 0.4e1 - t1148 / 0.4e1) * t1188, -t1276 + t441 * qJD(3) + t20 * qJD(4) + t64 * qJD(5) + (-t1116 / 0.4e1 - t1138 / 0.4e1 - t1151 / 0.4e1) * t1186 + ((-t786 * t336 + t987) * t1180 + (-t786 * t411 + t979) * t1182 + (-t786 * t472 + t971) * t1183 + ((t485 * t784 + t487 * t787 + t281) * t1180 + (t540 * t784 + t542 * t787 + t356) * t1182 + (t614 * t784 + t616 * t787 + t447) * t1183) * t783) * t1187, t441 * qJD(2), t31 * qJD(1) + t20 * qJD(2) + ((t446 * t786 + (t490 * t787 + t491 * t784) * t783) * t1182 + (-t383 * t786 + (t418 * t787 + t419 * t784) * t783) * t1180) * t1185 + t274, t109 * qJD(1) + t64 * qJD(2) + qJD(4) * t1108 + t274; t2 * qJD(2) + t32 * qJD(3) + t12 * qJD(4) + t14 * qJD(5) + (-t1111 / 0.4e1 - t1133 / 0.4e1) * t1188 - t1194 * t941 + (-t1042 / 0.2e1 + t813 - 0.2e1 * t423 * t997 * t1180) * qJD(1), t2 * qJD(1) + (t812 + t142 * t905 + t1265 * t914 + t141 * t909 + (t345 * t784 - t346 * t787) * t917 + (t287 * t784 - t288 * t787) * t1163 + t76 * t1159 + (t384 * t784 + t385 * t787) * t1160 + t77 * t1162 - t825) * qJD(2) + t21 * qJD(3) + t8 * qJD(4) + t10 * qJD(5) + (-t1121 / 0.4e1 - t1141 / 0.4e1) * t1186 + ((t350 * t356 - t395 * t541 - t396 * t539 - t411 * t430 + t476 * t542 - t477 * t540) * t1182 + (t222 * t281 - t272 * t486 - t273 * t484 + t307 * t336 + t377 * t487 - t379 * t485) * t1180) * t1187, qJD(1) * t32 + qJD(2) * t21, t12 * qJD(1) + t8 * qJD(2) + (t116 * t905 + t117 * t909 + (t408 + (t299 * t787 + t300 * t784) * t786) * t1163 + t934) * qJD(4) + t1218 + ((t307 * t383 + t377 * t418 - t379 * t419) * t1179 + (t430 * t446 + t476 * t490 - t477 * t491) * t1181) * t1184, t14 * qJD(1) + t10 * qJD(2) + t22 * qJD(4) + t1218; (t813 - t1274) * qJD(1) + t5 * qJD(2) + t111 * qJD(3) + t15 * qJD(4) + t878 * qJD(5), t5 * qJD(1) + ((t281 * t322 + t336 * t410 - t380 * t486 - t381 * t484 + t443 * t487 - t444 * t485 - t156) * m(6) + t812) * qJD(2) + t66 * qJD(3) + t11 * qJD(4) + t16 * qJD(5), qJD(1) * t111 + qJD(2) * t66, t15 * qJD(1) + t11 * qJD(2) + ((t383 * t410 + t418 * t443 - t419 * t444) * m(6) + t934) * qJD(4) + t23, qJD(1) * t878 + qJD(2) * t16 + qJD(4) * t24 + t23;];
Cq = t17;
