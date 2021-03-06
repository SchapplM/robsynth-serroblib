% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:14
% EndTime: 2019-03-09 02:06:06
% DurationCPUTime: 47.15s
% Computational Cost: add. (73764->1126), mult. (165902->1599), div. (0->0), fcn. (211618->8), ass. (0->687)
t608 = cos(qJ(4));
t1011 = -t608 / 0.2e1;
t897 = sin(pkin(9));
t898 = cos(pkin(9));
t942 = sin(qJ(1));
t943 = cos(qJ(1));
t581 = -t897 * t942 - t898 * t943;
t582 = t897 * t943 - t898 * t942;
t605 = sin(qJ(5));
t607 = cos(qJ(5));
t846 = t607 * t608;
t517 = -t581 * t846 + t582 * t605;
t606 = sin(qJ(4));
t866 = t581 * t606;
t780 = rSges(7,2) * t866;
t842 = t608 * t605;
t516 = -t581 * t842 - t582 * t607;
t903 = t516 * rSges(7,3);
t634 = t608 * (t517 * rSges(7,1) - t780 + t903);
t499 = t516 * qJ(6);
t704 = t608 * (t517 * pkin(5) + t499);
t613 = -t704 - t634;
t699 = rSges(7,1) * t607 + rSges(7,3) * t605;
t649 = t699 * t606;
t1030 = -rSges(7,2) * t608 + t649;
t698 = pkin(5) * t607 + qJ(6) * t605;
t648 = t698 * t606;
t800 = -t648 - t1030;
t712 = t800 * t606;
t267 = -t581 * t712 + t613;
t1014 = -0.2e1 * t267;
t1050 = 2 * m(6);
t1067 = -t1050 / 0.4e1;
t917 = pkin(4) * t608;
t596 = -t606 * pkin(8) - t917;
t528 = t596 * t582;
t487 = t581 * t528;
t862 = t582 * t608;
t863 = t582 * t606;
t526 = pkin(4) * t862 + pkin(8) * t863;
t512 = -t581 * t607 + t582 * t842;
t513 = t581 * t605 + t582 * t846;
t554 = rSges(7,2) * t863;
t385 = -t513 * rSges(7,1) - t512 * rSges(7,3) - t554;
t431 = -t513 * pkin(5) - t512 * qJ(6);
t812 = t431 + t385;
t900 = rSges(7,3) + qJ(6);
t946 = rSges(7,1) + pkin(5);
t1051 = t512 * t900 + t513 * t946;
t817 = t1051 + t554;
t709 = t812 + t817;
t697 = t526 + t709;
t177 = t581 * t697 + t487;
t198 = t709 * t863;
t224 = t526 * t581 + t487;
t553 = rSges(6,3) * t863;
t702 = t513 * rSges(6,1) - t512 * rSges(6,2);
t386 = -t702 - t553;
t701 = t517 * rSges(6,1) - t516 * rSges(6,2);
t277 = (t582 * t701 + (-t386 - t553) * t581) * t606;
t906 = rSges(6,1) * t607;
t700 = -rSges(6,2) * t605 + t906;
t650 = t700 * t606;
t1025 = -rSges(6,3) * t608 + t650;
t640 = -rSges(6,3) * t866 + t701;
t635 = t608 * t640;
t321 = t1025 * t866 - t635;
t595 = -t606 * pkin(4) + pkin(8) * t608;
t799 = -t1025 + t595;
t466 = t799 * t581;
t532 = t581 * t595;
t469 = t1025 * t581 - t532;
t808 = -t466 - t469;
t756 = t595 + t800;
t397 = t756 * t581;
t400 = -t581 * t800 - t532;
t814 = -t397 - t400;
t705 = t582 * t712;
t265 = t608 * t817 - t705;
t713 = t812 * t608;
t266 = -t713 - t705;
t828 = t265 - t266;
t398 = t756 * t582;
t983 = -0.2e1 * t398;
t645 = t946 * t517 + t499 + t903;
t199 = (t645 * t582 + (-t554 - t812) * t581) * t606;
t986 = 0.2e1 * t199;
t489 = t582 * t528;
t867 = t581 * t596;
t178 = t489 + t812 * t582 + (t867 + t645 - t780) * t581;
t988 = 0.2e1 * t178;
t1085 = -m(7) * (t1014 * t814 + t177 * t986 + t198 * t988 + t828 * t983) / 0.4e1 + (t224 * t277 - t321 * t808) * t1067;
t500 = Icges(7,5) * t513;
t366 = -Icges(7,6) * t863 - Icges(7,3) * t512 - t500;
t888 = Icges(7,5) * t512;
t378 = -Icges(7,1) * t513 - Icges(7,4) * t863 - t888;
t1074 = t366 * t605 + t378 * t607;
t371 = Icges(7,4) * t513 + Icges(7,2) * t863 + Icges(7,6) * t512;
t876 = t371 * t608;
t227 = t1074 * t606 + t876;
t416 = Icges(6,5) * t512 + Icges(6,6) * t513;
t418 = Icges(7,4) * t512 - Icges(7,6) * t513;
t1083 = -t416 - t418;
t417 = -Icges(6,5) * t516 - Icges(6,6) * t517;
t419 = -Icges(7,4) * t516 + Icges(7,6) * t517;
t1082 = -t417 - t419;
t845 = t608 * t386;
t1081 = t1025 * t863 - t845;
t503 = Icges(6,4) * t516;
t382 = Icges(6,1) * t517 - Icges(6,5) * t866 - t503;
t818 = -Icges(6,2) * t517 + t382 - t503;
t887 = Icges(7,5) * t516;
t379 = Icges(7,1) * t517 - Icges(7,4) * t866 + t887;
t820 = -Icges(7,3) * t517 + t379 + t887;
t1080 = t818 + t820;
t502 = Icges(6,4) * t512;
t381 = -Icges(6,1) * t513 - Icges(6,5) * t863 + t502;
t819 = Icges(6,2) * t513 + t381 + t502;
t821 = Icges(7,3) * t513 + t378 - t888;
t1079 = t819 + t821;
t893 = Icges(6,4) * t517;
t376 = -Icges(6,2) * t516 - Icges(6,6) * t866 + t893;
t822 = Icges(6,1) * t516 + t376 + t893;
t501 = Icges(7,5) * t517;
t367 = -Icges(7,6) * t866 + Icges(7,3) * t516 + t501;
t824 = -Icges(7,1) * t516 + t367 + t501;
t1078 = t822 - t824;
t894 = Icges(6,4) * t513;
t375 = Icges(6,2) * t512 - Icges(6,6) * t863 - t894;
t823 = -Icges(6,1) * t512 + t375 - t894;
t825 = Icges(7,1) * t512 + t366 - t500;
t1077 = t823 - t825;
t675 = t375 * t605 - t381 * t607;
t368 = Icges(6,5) * t513 - Icges(6,6) * t512 + Icges(6,3) * t863;
t878 = t368 * t608;
t231 = t606 * t675 - t878;
t1076 = -t366 * t512 - t378 * t513;
t1075 = -t516 * t366 - t517 * t378;
t1073 = t375 * t512 - t381 * t513;
t1072 = t516 * t375 - t517 * t381;
t911 = m(7) * qJD(1);
t860 = t605 * t606;
t599 = Icges(6,4) * t860;
t847 = t606 * t607;
t889 = Icges(6,5) * t608;
t545 = -Icges(6,1) * t847 + t599 + t889;
t892 = Icges(6,4) * t607;
t690 = -Icges(6,2) * t605 + t892;
t630 = -Icges(6,6) * t608 + t606 * t690;
t686 = Icges(6,5) * t607 - Icges(6,6) * t605;
t628 = -Icges(6,3) * t608 + t606 * t686;
t849 = t606 * t628;
t309 = t512 * t630 + t513 * t545 - t582 * t849;
t1071 = t309 * t608;
t885 = Icges(7,5) * t607;
t685 = Icges(7,3) * t605 + t885;
t629 = -Icges(7,6) * t608 + t606 * t685;
t886 = Icges(7,5) * t605;
t693 = Icges(7,1) * t607 + t886;
t633 = -Icges(7,4) * t608 + t606 * t693;
t689 = Icges(7,4) * t607 + Icges(7,6) * t605;
t631 = -Icges(7,2) * t608 + t606 * t689;
t848 = t606 * t631;
t310 = t512 * t629 + t513 * t633 + t582 * t848;
t1070 = t608 * t310;
t1069 = t227 - t231;
t676 = -t367 * t605 - t379 * t607;
t373 = Icges(7,4) * t517 - Icges(7,2) * t866 + Icges(7,6) * t516;
t875 = t373 * t608;
t229 = t606 * t676 + t875;
t674 = t376 * t605 - t382 * t607;
t370 = Icges(6,5) * t517 - Icges(6,6) * t516 - Icges(6,3) * t866;
t877 = t370 * t608;
t232 = t606 * t674 + t877;
t1068 = -t229 - t232;
t688 = Icges(5,5) * t608 - Icges(5,6) * t606;
t477 = Icges(5,3) * t582 - t581 * t688;
t1065 = t582 * t477;
t1064 = t1077 * t513 + t1079 * t512 + t1083 * t863;
t1063 = t1078 * t513 + t1080 * t512 + t1082 * t863;
t1062 = -t1077 * t517 - t1079 * t516 + t1083 * t866;
t1061 = -t1078 * t517 - t1080 * t516 + t1082 * t866;
t568 = (Icges(7,4) * t605 - Icges(7,6) * t607) * t606;
t566 = (-Icges(7,3) * t607 + t886) * t606;
t802 = -t633 - t566;
t570 = (Icges(7,1) * t605 - t885) * t606;
t804 = -t629 + t570;
t253 = t512 * t802 - t513 * t804 - t568 * t863;
t567 = (Icges(6,5) * t605 + Icges(6,6) * t607) * t606;
t569 = Icges(6,2) * t847 + t599;
t801 = t545 + t569;
t571 = (Icges(6,1) * t605 + t892) * t606;
t803 = -t630 - t571;
t254 = t512 * t801 + t513 * t803 - t567 * t863;
t1060 = t253 + t254;
t255 = -t516 * t802 + t517 * t804 - t568 * t866;
t256 = -t516 * t801 - t517 * t803 - t567 * t866;
t1059 = t256 + t255;
t859 = t606 * t368;
t1015 = t582 * t859 + t1073;
t1058 = t1015 * t581;
t1057 = t1015 * t582;
t857 = t606 * t371;
t1016 = t582 * t857 + t1076;
t1056 = t1016 * t581;
t1055 = t1016 * t582;
t1049 = 0.2e1 * t512;
t873 = t516 * t267;
t1054 = t1049 * t266 - 0.2e1 * t873;
t872 = t516 * t398;
t1053 = t1049 * t400 - 0.2e1 * t872;
t895 = Icges(5,4) * t608;
t692 = -Icges(5,2) * t606 + t895;
t478 = Icges(5,6) * t581 + t582 * t692;
t896 = Icges(5,4) * t606;
t696 = Icges(5,1) * t608 - t896;
t481 = Icges(5,5) * t581 + t582 * t696;
t1052 = t478 * t606 - t481 * t608;
t977 = m(5) / 0.4e1;
t976 = m(6) / 0.4e1;
t974 = m(7) / 0.4e1;
t783 = 0.2e1 * m(7);
t1010 = -t783 / 0.4e1;
t717 = t783 / 0.4e1;
t720 = t1050 / 0.4e1;
t785 = 0.2e1 * m(5);
t1048 = t785 / 0.4e1;
t204 = t581 * t857 - t1075;
t1046 = t204 * t581;
t1045 = t204 * t582;
t206 = t581 * t859 - t1072;
t1044 = t206 * t581;
t1043 = t206 * t582;
t706 = -pkin(1) * t942 + t943 * qJ(2);
t643 = -pkin(2) * t942 + t706;
t639 = t581 * pkin(7) + t643;
t755 = pkin(3) + t917;
t945 = rSges(7,2) + pkin(8);
t995 = t606 * t945 + t755;
t259 = t582 * t995 + t1051 + t639;
t1042 = t259 * t512;
t703 = t606 * rSges(5,1) + rSges(5,2) * t608;
t524 = t703 * t582;
t1041 = t524 * t942;
t525 = t703 * t581;
t707 = t525 * t943;
t448 = t629 * t582;
t456 = t633 * t582;
t660 = -t631 * t582 + t1074;
t180 = -t660 * t608 + (-t448 * t605 - t456 * t607 + t371) * t606;
t454 = t630 * t582;
t694 = Icges(6,1) * t607 - Icges(6,4) * t605;
t632 = t606 * t694 - t889;
t458 = t632 * t582;
t658 = t628 * t582 + t675;
t182 = t658 * t608 + (t454 * t605 - t458 * t607 + t368) * t606;
t1040 = -t180 - t182;
t449 = t629 * t581;
t457 = t633 * t581;
t659 = -t631 * t581 - t676;
t181 = -t659 * t608 + (-t449 * t605 - t457 * t607 - t373) * t606;
t455 = t630 * t581;
t459 = t632 * t581;
t657 = t628 * t581 + t674;
t183 = t657 * t608 + (t455 * t605 - t459 * t607 - t370) * t606;
t1039 = -t181 - t183;
t189 = t418 * t608 + (t605 * t821 - t607 * t825) * t606;
t191 = t416 * t608 + (t605 * t819 + t607 * t823) * t606;
t1038 = -t189 - t191;
t190 = t419 * t608 + (t605 * t820 - t607 * t824) * t606;
t192 = t417 * t608 + (t605 * t818 + t607 * t822) * t606;
t1037 = -t190 - t192;
t536 = -Icges(7,6) * t606 - t608 * t685;
t544 = -Icges(7,4) * t606 - t608 * t693;
t540 = -Icges(7,2) * t606 - t608 * t689;
t670 = t605 * t629 + t607 * t633;
t656 = -t540 - t670;
t869 = t631 * t608;
t620 = t606 * t656 + t869;
t248 = -t512 * t536 - t513 * t544 + t582 * t620;
t542 = -Icges(6,6) * t606 - t608 * t690;
t546 = -Icges(6,5) * t606 - t608 * t694;
t538 = -Icges(6,3) * t606 - t608 * t686;
t669 = -t607 * t545 - t605 * t630;
t655 = t538 + t669;
t870 = t628 * t608;
t621 = -t606 * t655 + t870;
t249 = t512 * t542 - t513 * t546 + t582 * t621;
t1036 = t248 + t249;
t250 = t516 * t536 + t517 * t544 + t581 * t620;
t251 = -t516 * t542 + t517 * t546 + t581 * t621;
t1035 = t250 + t251;
t1034 = t310 - t309;
t312 = -t516 * t629 - t517 * t633 + t581 * t848;
t313 = t516 * t630 + t517 * t545 + t581 * t849;
t1033 = t312 + t313;
t314 = -t656 * t608 + (-t605 * t536 - t607 * t544 + t631) * t606;
t315 = t655 * t608 + (t605 * t542 - t607 * t546 + t628) * t606;
t1032 = t314 + t315;
t435 = t606 * t670 - t869;
t436 = t606 * t669 - t870;
t1031 = t435 + t436;
t1029 = t1068 * t581 + t1069 * t582;
t865 = t581 * t608;
t558 = rSges(7,2) * t865;
t562 = pkin(8) * t865;
t993 = -t605 * t900 - t607 * t946;
t322 = -t558 - t562 + (pkin(4) - t993) * t866;
t561 = pkin(4) * t863;
t323 = -t561 + (t606 * t993 + t945 * t608) * t582;
t944 = rSges(6,3) + pkin(8);
t390 = -t561 + (t608 * t944 - t650) * t582;
t764 = t581 * t860;
t805 = -rSges(6,2) * t764 - rSges(6,3) * t865;
t391 = -t562 + (pkin(4) + t906) * t866 + t805;
t759 = (-t322 * t943 + t323 * t942) * t717 + (t390 * t942 - t391 * t943) * t720 + (-t707 - t1041) * t1048;
t749 = t398 * t942;
t786 = 0.2e1 * t943;
t291 = t400 * t786 - 0.2e1 * t749;
t467 = t799 * t582;
t746 = t467 * t942;
t342 = t469 * t786 - 0.2e1 * t746;
t531 = 0.2e1 * t707;
t787 = -0.2e1 * t943;
t760 = (-t397 * t787 + t291 + 0.2e1 * t749) * t974 + (-t466 * t787 + t342 + 0.2e1 * t746) * t976 + (t531 - 0.2e1 * t707) * t977;
t1028 = t759 - t760;
t753 = t267 * t942;
t179 = t266 * t786 - 0.2e1 * t753;
t751 = t321 * t942;
t240 = t1081 * t786 - 0.2e1 * t751;
t833 = t179 * t974 + t240 * t976;
t899 = (t1081 * t787 + t240 + 0.2e1 * t751) * t976 + (t265 * t787 + t179 + 0.2e1 * t753) * t974;
t1027 = t833 - t899;
t460 = t1030 * t582;
t496 = t582 * t648;
t809 = t460 + t496;
t1024 = t606 * t809 + t713;
t646 = pkin(1) * t943 + qJ(2) * t942;
t626 = pkin(2) * t943 + t646;
t1019 = -m(3) * ((rSges(3,3) * t943 + t706) * t943 + (rSges(3,3) * t942 + t646) * t942) - m(4) * ((-t581 * rSges(4,1) - t582 * rSges(4,2) + t626) * t942 + t943 * (t582 * rSges(4,1) - t581 * rSges(4,2) + t643));
t1018 = t367 * t512 + t379 * t513;
t1017 = -t376 * t512 + t382 * t513;
t957 = -t581 / 0.2e1;
t1013 = t581 / 0.2e1;
t1012 = -t582 / 0.2e1;
t953 = t582 / 0.2e1;
t950 = -t606 / 0.2e1;
t1008 = t1060 * t608 + (-t1063 * t581 - t1064 * t582) * t606;
t1007 = t1059 * t608 + (-t1061 * t581 - t1062 * t582) * t606;
t428 = pkin(5) * t512 - qJ(6) * t513;
t429 = rSges(7,1) * t512 - rSges(7,3) * t513;
t813 = t428 + t429;
t1006 = t582 * t813;
t871 = t516 * t582;
t874 = t512 * t581;
t360 = (-t871 - t874) * t606;
t1004 = t800 * t608;
t1003 = ((t568 + t567) * t608 + ((t803 - t804) * t607 + (t801 + t802) * t605) * t606) * t608;
t1002 = t183 / 0.2e1 + t181 / 0.2e1;
t1001 = -t182 / 0.2e1 - t180 / 0.2e1;
t302 = -t512 * t946 + t513 * t900;
t432 = -t516 * pkin(5) + qJ(6) * t517;
t433 = -t516 * rSges(7,1) + t517 * rSges(7,3);
t303 = t432 + t433;
t430 = rSges(6,1) * t512 + rSges(6,2) * t513;
t434 = -rSges(6,1) * t516 - rSges(6,2) * t517;
t832 = (t302 * t942 - t303 * t943) * t717 + (-t430 * t942 - t434 * t943) * t720;
t618 = t582 * pkin(7) + t626;
t260 = -t581 * t995 + t618 + t645;
t996 = t606 * t944 + t755;
t293 = t582 * t996 + t639 + t702;
t294 = -t581 * t996 + t618 + t701;
t797 = (t605 * t946 - t607 * t900) * t606;
t464 = t797 * t582;
t465 = t797 * t581;
t578 = (rSges(6,1) * t605 + rSges(6,2) * t607) * t606;
t999 = (-t430 * t469 - t434 * t467 + (-t293 * t581 - t294 * t582) * t578) * t1067 + (-t259 * t465 - t260 * t464 + t302 * t400 - t303 * t398) * t1010;
t902 = t606 * rSges(7,2);
t798 = t902 + (t698 + t699) * t608;
t195 = -t809 * t608 + t812 * t606 + (t606 * t798 - t1004) * t582;
t462 = t581 * t649 - t558;
t497 = t581 * t648;
t806 = t497 + t462;
t196 = t806 * t608 - t645 * t606 + (t1004 + (-t798 + t902) * t606) * t581;
t461 = t1025 * t582;
t901 = t606 * rSges(6,3);
t550 = -t608 * t700 - t901;
t843 = t608 * t1025;
t271 = t606 * t386 - t461 * t608 + (-t550 * t606 + t843) * t582;
t463 = rSges(6,1) * t581 * t847 + t805;
t272 = -t606 * t701 + t608 * t463 + (-t843 + (t550 + t901) * t606) * t581;
t998 = (t195 * t259 + t196 * t260 + t266 * t323 - t267 * t322) * t1010 + (t1081 * t390 + t271 * t293 + t272 * t294 - t321 * t391) * t1067;
t475 = Icges(5,3) * t581 + t582 * t688;
t758 = t291 * t974 + t342 * t976 + (t531 + 0.2e1 * t1041) * t977;
t724 = t545 / 0.2e1 - t633 / 0.2e1;
t726 = -t630 / 0.2e1 + t629 / 0.2e1;
t994 = t605 * (t569 / 0.2e1 - t566 / 0.2e1 + t724) + t607 * (-t571 / 0.2e1 - t570 / 0.2e1 + t726);
t695 = Icges(5,1) * t606 + t895;
t992 = t605 * t726 - t607 * t724 + t538 / 0.2e1 + t540 / 0.2e1 + t692 / 0.2e1 + t695 / 0.2e1;
t691 = Icges(5,2) * t608 + t896;
t991 = t605 * (t542 / 0.2e1 - t536 / 0.2e1) - t607 * (t546 / 0.2e1 + t544 / 0.2e1) + t628 / 0.2e1 + t631 / 0.2e1 - t691 / 0.2e1 + t696 / 0.2e1;
t990 = t581 ^ 2;
t989 = t582 ^ 2;
t987 = 0.4e1 * t178;
t225 = t582 * t386 + t489 + (t867 + t640) * t581;
t985 = 0.4e1 * t225;
t984 = 0.2e1 * t260;
t861 = t606 ^ 2 * t605;
t446 = t512 * t608 + t582 * t861;
t982 = 0.2e1 * t446;
t979 = -0.2e1 * t512;
t978 = 0.2e1 * t516;
t975 = m(7) / 0.2e1;
t447 = t516 * t608 - t581 * t861;
t95 = -0.4e1 * t199 * t360 + 0.4e1 * t266 * t446 - 0.4e1 * t267 * t447;
t973 = t95 / 0.4e1;
t143 = -t1024 * t581 + (t645 * t608 + (-t558 + t806) * t606) * t582;
t710 = 0.2e1 * t764;
t243 = t266 * t710;
t782 = 0.2e1 * t606;
t969 = m(7) * (t195 * t978 + t196 * t979 + t243 + (-0.2e1 * t199 * t608 + (-t267 * t582 - t143) * t782) * t605);
t864 = t582 * t512;
t868 = t581 * t516;
t362 = -t864 + t868;
t763 = t582 * t860;
t967 = m(7) * (t1014 * t763 - t360 * t988 + t362 * t986 + t400 * t982 + t447 * t983 + t243);
t781 = -0.2e1 * t860;
t966 = m(7) * (t198 * t781 + t265 * t979 + t1054 + 0.2e1 * t873);
t212 = t303 * t581 + t1006;
t965 = m(7) * (-t513 * t983 + 0.2e1 * t400 * t517 - t464 * t979 - t465 * t978 + (-t178 * t607 - t212 * t605) * t782);
t964 = m(7) * (t259 * t982 + t447 * t984 + t1054);
t958 = t478 / 0.2e1;
t956 = -t581 / 0.4e1;
t952 = t582 / 0.4e1;
t951 = -t688 / 0.2e1;
t948 = t608 / 0.2e1;
t929 = m(7) * (t177 * t781 - t397 * t979 + t1053 + 0.2e1 * t872);
t928 = (t259 * t517 - t260 * t513 + t302 * t516 - t303 * t512) * t783;
t665 = t259 * t710 + t763 * t984;
t927 = m(7) * (t322 * t979 + t323 * t978 + t665);
t926 = m(7) * (t1053 + t665);
t741 = t943 * t516;
t744 = t512 * t942;
t920 = (t744 - t741) * t783;
t498 = -0.2e1 * t744;
t919 = m(7) * (0.2e1 * t741 + t498);
t918 = m(7) * (t1049 * t942 + t498);
t915 = m(5) * qJD(1);
t914 = m(6) * qJD(1);
t913 = m(6) * qJD(4);
t912 = m(6) * qJD(5);
t910 = m(7) * qJD(4);
t909 = m(7) * qJD(5);
t908 = m(7) * qJD(6);
t907 = rSges(5,1) * t608;
t572 = t581 * rSges(5,3);
t858 = t606 * t370;
t856 = t606 * t373;
t855 = t606 * t461;
t854 = t606 * t462;
t853 = t606 * t463;
t480 = Icges(5,6) * t582 - t581 * t692;
t851 = t606 * t480;
t850 = t606 * t497;
t840 = (-t371 * t581 + t373 * t582) * t606 + t204 + t1018 + t1075;
t838 = (-t368 * t581 + t370 * t582) * t606 + t206 + t1017 + t1072;
t205 = t516 * t367 + t517 * t379 - t581 * t856;
t836 = -t371 * t863 + t1016 - t1076 + t205;
t207 = -t516 * t376 + t517 * t382 - t581 * t858;
t834 = -t368 * t863 + t1015 - t1073 + t207;
t619 = t582 * pkin(3) + t639;
t614 = t619 + t526;
t261 = t614 + t817;
t829 = t259 - t261;
t292 = t614 - t386;
t827 = t292 - t293;
t593 = t606 * rSges(5,2) - t907;
t392 = t572 + (pkin(3) - t593) * t582 + t639;
t666 = rSges(5,1) * t862 - rSges(5,2) * t863 + t572;
t394 = t619 + t666;
t815 = t392 - t394;
t483 = Icges(5,5) * t582 - t581 * t696;
t811 = -t581 * t477 - t483 * t862;
t810 = -t483 * t865 + t1065;
t796 = qJD(1) * t606;
t795 = qJD(1) * t608;
t794 = qJD(4) * t581;
t793 = qJD(4) * t582;
t792 = qJD(5) * t606;
t484 = t582 * t593 - t572;
t647 = -t484 - t666;
t109 = (-m(6) / 0.2e1 - m(7) / 0.2e1) * t487 + (t697 * t1010 + t1048 * t647 - t526 * t720) * t581;
t791 = t109 * qJD(3);
t148 = t1010 * t198;
t790 = t148 * qJD(3);
t735 = t842 / 0.2e1;
t350 = (t735 - t864 / 0.2e1 + t868 / 0.2e1) * m(7);
t789 = t350 * qJD(3);
t623 = t606 * t660 + t876;
t149 = -t448 * t512 - t456 * t513 + t582 * t623;
t622 = t606 * t659 - t875;
t150 = -t449 * t512 - t457 * t513 + t582 * t622;
t625 = -t606 * t658 + t878;
t151 = t454 * t512 - t458 * t513 + t582 * t625;
t624 = -t606 * t657 - t877;
t152 = t455 * t512 - t459 * t513 + t582 * t624;
t203 = -t582 * t858 - t1017;
t682 = -t203 * t581 - t1057;
t201 = -t582 * t856 - t1018;
t683 = -t201 * t581 - t1055;
t779 = ((-t149 - t151) * t582 + (-t150 - t152) * t581 - t1034) * t950 + (t683 + t682 + t1036) * t1011;
t153 = t516 * t448 + t517 * t456 + t581 * t623;
t154 = t516 * t449 + t517 * t457 + t581 * t622;
t155 = -t516 * t454 + t517 * t458 + t581 * t625;
t156 = -t516 * t455 + t517 * t459 + t581 * t624;
t680 = -t207 * t581 - t1043;
t681 = -t205 * t581 - t1045;
t778 = ((-t153 - t155) * t582 + (-t154 - t156) * t581 - t1033) * t950 + (t681 + t680 + t1035) * t1011;
t777 = (t1039 * t581 + t1040 * t582 - t1031) * t950 + (t1029 + t1032) * t1011;
t776 = t1012 * t1063 + t1013 * t1064;
t775 = t1061 * t953 + t1062 * t957;
t761 = t908 / 0.4e1;
t754 = t259 * t943;
t752 = t293 * t943;
t750 = t392 * t943;
t748 = t446 * t942;
t747 = t447 * t943;
t743 = t513 * t943;
t742 = t517 * t942;
t740 = -t866 / 0.4e1;
t738 = -t863 / 0.4e1;
t111 = t606 * t683 + t1070;
t112 = t606 * t682 - t1071;
t734 = t111 / 0.2e1 + t112 / 0.2e1;
t304 = t312 * t608;
t113 = t606 * t681 + t304;
t305 = t313 * t608;
t114 = t606 * t680 + t305;
t733 = -t114 / 0.2e1 - t113 / 0.2e1;
t126 = t201 * t582 - t1056;
t127 = t203 * t582 - t1058;
t732 = -t126 / 0.2e1 - t127 / 0.2e1;
t128 = t205 * t582 - t1046;
t129 = t207 * t582 - t1044;
t731 = -t128 / 0.2e1 - t129 / 0.2e1;
t730 = t1011 * t1031 + t1029 * t950;
t729 = -t189 / 0.2e1 - t191 / 0.2e1;
t728 = t190 / 0.2e1 + t192 / 0.2e1;
t722 = t567 / 0.2e1 + t568 / 0.2e1;
t721 = t1050 / 0.2e1;
t718 = t783 / 0.2e1;
t715 = rSges(5,2) * t866 + t582 * rSges(5,3);
t711 = 0.4e1 * t860;
t708 = 0.2e1 * t303;
t687 = Icges(5,5) * t606 + Icges(5,6) * t608;
t55 = ((t432 / 0.2e1 + t433 / 0.2e1 + (t385 / 0.2e1 + t431 / 0.2e1) * t608 + (t460 / 0.2e1 + t496 / 0.2e1) * t606) * t581 + (-t634 / 0.2e1 - t854 / 0.2e1 - t704 / 0.2e1 - t850 / 0.2e1 + t428 / 0.2e1 + t429 / 0.2e1) * t582) * m(7) + ((t845 / 0.2e1 + t855 / 0.2e1 + t434 / 0.2e1) * t581 + (-t635 / 0.2e1 - t853 / 0.2e1 + t430 / 0.2e1) * t582) * m(6);
t644 = t581 * t942 - t582 * t943;
t636 = 0.2e1 * t644;
t609 = (t464 * t943 - t465 * t942) * t717 - t578 * t636 * t976;
t610 = (t271 * t942 - t272 * t943) * t1067 + (t195 * t942 - t196 * t943) * t1010;
t64 = t609 + t610;
t684 = -t64 * qJD(2) + t55 * qJD(3);
t673 = t845 + t855;
t672 = -t398 * t582 + t400 * t581;
t289 = t430 * t582 + t434 * t581;
t279 = (t748 / 0.2e1 - t747 / 0.2e1 - t742 / 0.2e1 - t743 / 0.2e1) * m(7);
t344 = (t607 / 0.2e1 + t874 / 0.2e1 + t871 / 0.2e1) * t606 * m(7);
t668 = t279 * qJD(2) - t344 * qJD(3);
t667 = qJD(4) * t109 + qJD(5) * t148;
t664 = t776 - t779;
t663 = -t775 - t778;
t30 = t304 + (-t581 * t836 - t1045) * t606;
t31 = t305 + (-t581 * t834 - t1043) * t606;
t662 = -t30 / 0.2e1 - t31 / 0.2e1 - t733;
t28 = -t1070 + (-t581 * t840 + t1055) * t606;
t29 = t1071 + (-t581 * t838 + t1057) * t606;
t661 = -t28 / 0.2e1 - t29 / 0.2e1 - t734;
t520 = t691 * t582;
t522 = t695 * t582;
t642 = (t478 + t522) * t608 + (t481 - t520) * t606;
t521 = t691 * t581;
t523 = t695 * t581;
t641 = (t480 - t523) * t608 + (t483 + t521) * t606;
t281 = t582 * t851 + t811;
t283 = t581 * t851 + t810;
t43 = t582 * t836 - t1046;
t44 = t582 * t834 - t1044;
t638 = -t43 / 0.2e1 + t283 * t953 + (t281 - t811) * t957 + t810 * t1012 - t44 / 0.2e1 - t731;
t41 = t582 * t840 + t1056;
t42 = t582 * t838 + t1058;
t637 = (t283 + (t475 - t851) * t581 + t1065 - t810) * t1013 + (-t1052 * t582 + t581 * t475) * t957 + t42 / 0.2e1 + t41 / 0.2e1 - t732 + (-t1052 * t581 + t281 + (t483 * t608 - t851) * t582) * t953;
t627 = -t999 + (-t1038 + t1060) * t956 + (-t1037 + t1059) * t952;
t615 = t635 + t853;
t612 = (t30 + t31) * t956 + (t113 + t114) * t581 / 0.4e1 + (t43 + t44) * t738 + (t128 + t129) * t863 / 0.4e1 + (t111 + t112 + t28 + t29) * t952 + (t126 + t127 + t41 + t42) * t740 - t1085;
t611 = -t998 + t1031 * t950 + t1032 * t948 + (t1035 - t1039) * t740 - (t1033 - t1068) * t865 / 0.4e1 + (t1036 - t1040) * t738 - (t1034 - t1069) * t862 / 0.4e1;
t529 = pkin(4) * t866 - t562;
t519 = t687 * t581;
t518 = t687 * t582;
t488 = t582 * (-pkin(8) * t862 + t561);
t470 = -t550 * t581 - t867;
t468 = (-t550 - t596) * t582;
t404 = t918 / 0.4e1;
t403 = t919 / 0.4e1;
t402 = t920 / 0.4e1;
t401 = t581 * t798 - t867;
t399 = (-t596 + t798) * t582;
t393 = (-pkin(3) - t907) * t581 + t618 + t715;
t388 = t644 * t718 * t860;
t351 = m(7) * t735 + t1010 * t362;
t349 = t434 * t608 + t578 * t866;
t348 = -t430 * t608 - t578 * t863;
t345 = t360 * t717 + t847 * t975;
t343 = (t362 + t842) * t711;
t337 = t582 * t484 + t581 * (-rSges(5,1) * t865 + t715);
t301 = 0.4e1 * t512 * t513 + 0.4e1 * t516 * t517 + 0.4e1 * t607 * t861;
t287 = (-t430 * t581 + t434 * t582) * t606;
t286 = 0.4e1 * t393 * t942 + 0.4e1 * t750;
t278 = (t748 - t747 + t742 + t743) * t717;
t275 = t303 * t608 + t797 * t866;
t274 = -t464 * t606 - t608 * t813;
t268 = t461 * t582 + t488 + (t463 + t529) * t581;
t264 = t403 + t404 - t920 / 0.4e1;
t263 = t402 + t404 - t919 / 0.4e1;
t262 = t402 + t403 - t918 / 0.4e1;
t252 = -0.4e1 * t392 * t524 + 0.4e1 * t393 * t525;
t217 = -0.4e1 * t337 * t647 * t581;
t216 = 0.4e1 * t294 * t942 + 0.4e1 * t752;
t214 = t488 + t809 * t582 + (t529 + t806) * t581;
t210 = -t581 * t673 + t582 * t615;
t209 = (t303 * t582 - t581 * t813) * t606;
t176 = -0.4e1 * t293 * t430 + 0.4e1 * t294 * t434;
t173 = 0.4e1 * t260 * t942 + 0.4e1 * t754;
t168 = 0.4e1 * t293 * t390 + 0.4e1 * t294 * t391;
t165 = 0.4e1 * t260 * t516 + 0.4e1 * t1042;
t147 = t148 * qJD(1);
t145 = t289 * t985 + 0.4e1 * (t467 * t582 - t469 * t581) * t578;
t131 = 0.4e1 * t259 * t323 + 0.4e1 * t260 * t322;
t125 = t926 / 0.4e1;
t120 = 0.4e1 * t259 * t302 + 0.4e1 * t260 * t303;
t119 = t362 * t987 + t672 * t711;
t117 = t224 * t985 - 0.4e1 * t467 * t808;
t115 = t927 / 0.4e1;
t108 = t109 * qJD(1);
t107 = t928 / 0.4e1;
t99 = t929 / 0.4e1;
t94 = 0.4e1 * t178 * t212 + 0.4e1 * t398 * t464 - 0.4e1 * t400 * t465;
t92 = t964 / 0.4e1;
t91 = t173 * t974 + t216 * t976 + t286 * t977 - t1019;
t86 = -t155 * t581 + t156 * t582;
t85 = -t153 * t581 + t154 * t582;
t84 = -t151 * t581 + t152 * t582;
t83 = -t149 * t581 + t150 * t582;
t71 = t965 / 0.4e1;
t70 = 0.4e1 * t1081 * t271 + 0.4e1 * t210 * t277 - 0.4e1 * t272 * t321;
t69 = t177 * t987 - 0.4e1 * t398 * t814;
t63 = t609 - t610;
t62 = t120 * t974 + t176 * t976 + t606 * t994 + t722 * t608;
t60 = t966 / 0.4e1;
t59 = t758 - t1028;
t58 = t759 + t760 - t758;
t57 = t758 + t1028;
t56 = t289 * t1067 + (-t581 * t708 - 0.2e1 * t1006) * t974 + (t1024 * t717 + t673 * t720) * t581 + (t615 * t1067 + (t613 - t850 - t854) * t717) * t582;
t54 = 0.4e1 * t198 * t199 - 0.4e1 * t267 * t828;
t49 = t131 * t974 + t168 * t976 + t252 * t977 + t606 * t991 + t608 * t992;
t47 = t967 / 0.4e1;
t36 = 0.4e1 * t143 * t199 + 0.4e1 * t195 * t266 - 0.4e1 * t196 * t267;
t22 = t969 / 0.4e1;
t21 = t115 + t99 - t926 / 0.4e1;
t20 = t125 + t115 - t929 / 0.4e1;
t19 = t125 + t99 - t927 / 0.4e1;
t18 = t833 + t899 - t832;
t17 = t832 - t1027;
t16 = t832 + t1027;
t13 = t107 + t60 - t964 / 0.4e1;
t12 = t92 + t107 - t966 / 0.4e1;
t11 = t92 + t60 - t928 / 0.4e1;
t10 = t71 + t22 - t967 / 0.4e1;
t9 = t47 + t71 - t969 / 0.4e1;
t8 = t47 + t22 - t965 / 0.4e1;
t7 = t145 * t976 + t581 * t776 + t582 * t775 + t94 * t974;
t6 = t117 * t976 + t217 * t977 + t581 * t638 + t582 * t637 + t69 * t974;
t5 = t54 * t974 + (t581 * t661 + t582 * t662) * t606;
t4 = t36 * t974 + t70 * t976 + (t581 * t733 - t582 * t734 - t777) * t608 + (t581 * t778 + t582 * t779 + t730) * t606;
t3 = (-t112 / 0.4e1 - t111 / 0.4e1 - t29 / 0.4e1 - t28 / 0.4e1 + (t44 / 0.4e1 + t43 / 0.4e1 - t129 / 0.4e1 - t128 / 0.4e1) * t606) * t582 + t611 + (-t114 / 0.4e1 - t113 / 0.4e1 + t31 / 0.4e1 + t30 / 0.4e1 + (t42 / 0.4e1 + t41 / 0.4e1 + t127 / 0.4e1 + t126 / 0.4e1) * t606) * t581 + t627 + t1085;
t2 = t612 + t611 + (t254 / 0.4e1 + t253 / 0.4e1 + t191 / 0.4e1 + t189 / 0.4e1) * t581 + (-t256 / 0.4e1 - t255 / 0.4e1 - t192 / 0.4e1 - t190 / 0.4e1) * t582 + t999;
t1 = (t436 / 0.2e1 + t435 / 0.2e1 + (t182 / 0.4e1 + t180 / 0.4e1 + t249 / 0.4e1 + t248 / 0.4e1) * t582 + (t183 / 0.4e1 + t181 / 0.4e1 + t251 / 0.4e1 + t250 / 0.4e1) * t581) * t606 + t612 + (-t315 / 0.2e1 - t314 / 0.2e1 + (t231 / 0.4e1 - t227 / 0.4e1 - t309 / 0.4e1 + t310 / 0.4e1) * t582 + (t232 / 0.4e1 + t229 / 0.4e1 + t313 / 0.4e1 + t312 / 0.4e1) * t581) * t608 + t627 + t998;
t14 = [(-m(5) * t815 * t393 + m(6) * t294 * t827 - m(7) * t260 * t829) * qJD(1) + t91 * qJD(2) + t49 * qJD(4) + t62 * qJD(5) + t165 * t761, qJD(1) * t91 + qJD(4) * t57 + qJD(5) * t16 + qJD(6) * t262, -t667, t49 * qJD(1) + t57 * qJD(2) - t791 + t3 * qJD(5) + t20 * qJD(6) + (t250 / 0.2e1 + t251 / 0.2e1 + t582 * t951 + (-t483 / 0.2e1 - t521 / 0.2e1) * t608 + (t480 / 0.2e1 - t523 / 0.2e1) * t606 - t637 + t1002) * t793 + (-t248 / 0.2e1 - t249 / 0.2e1 + t581 * t951 + (-t481 / 0.2e1 + t520 / 0.2e1) * t608 + (t958 + t522 / 0.2e1) * t606 - t638 + t1001) * t794 + (t259 * t401 + t260 * t399 - t322 * t398 + t323 * t400 - t69 / 0.4e1) * t910 + (t293 * t470 + t294 * t468 + t390 * t469 - t391 * t467 - t117 / 0.4e1) * t913 + (-t217 / 0.4e1 + (-t393 * t593 + t525 * t703) * t582 + (-t392 * t593 - t524 * t703) * t581) * qJD(4) * m(5), t62 * qJD(1) + t16 * qJD(2) - t790 + t3 * qJD(4) + t1003 * qJD(5) + t12 * qJD(6) + (-t54 / 0.4e1 + t259 * t274 + t260 * t275 + t266 * t302 - t267 * t303) * t909 + (-t1081 * t430 + t293 * t348 + t294 * t349 - t321 * t434) * t912 + ((-t254 / 0.2e1 - t253 / 0.2e1 - t662 + t729) * t582 + (-t256 / 0.2e1 - t255 / 0.2e1 - t661 - t728) * t581) * t792, t165 * t911 / 0.4e1 + t262 * qJD(2) + t20 * qJD(4) + t12 * qJD(5); t58 * qJD(4) + t17 * qJD(5) + t263 * qJD(6) + (t754 - t943 * t261 - t173 / 0.4e1) * t911 + (-t943 * t292 + t752 - t216 / 0.4e1) * t914 + (t750 - t943 * t394 - t286 / 0.4e1) * t915 + t1019 * qJD(1), 0, 0, t58 * qJD(1) + (-m(5) * t593 * t636 / 0.2e1 + (-t468 * t943 + t470 * t942) * t721 + (-t399 * t943 + t401 * t942) * t718) * qJD(4) + t63 * qJD(5) + t388 * qJD(6), t17 * qJD(1) + t63 * qJD(4) + ((t348 * t942 - t349 * t943) * t721 + (t274 * t942 - t275 * t943) * t718) * qJD(5) + t278 * qJD(6), qJD(1) * t263 + qJD(4) * t388 + qJD(5) * t278; t667, 0, 0, t108 + ((-m(5) * t524 - m(6) * t461 - t718 * t809) * t582 + (-m(5) * t525 - m(6) * t463 - t718 * t806) * t581 + 0.2e1 * (m(6) / 0.2e1 + t975) * (-t581 * t529 - t488)) * qJD(4) + t56 * qJD(5) + t351 * qJD(6), t147 + t56 * qJD(4) + ((-m(6) * t434 - t708 * t975) * t582 + (m(6) * t430 + t718 * t813) * t581) * t792 + t345 * qJD(6), qJD(4) * t351 + qJD(5) * t345; -t252 * t915 / 0.4e1 + t59 * qJD(2) + t791 + t6 * qJD(4) + t1 * qJD(5) + t19 * qJD(6) + (-t131 / 0.4e1 + t829 * t398 + t814 * t260) * t911 + (-t168 / 0.4e1 - t827 * t467 + t808 * t294) * t914 - t992 * t795 - t991 * t796 + ((t958 - t478 / 0.2e1) * t608 - t815 * t703 * t785 / 0.2e1) * qJD(1) * t582, qJD(1) * t59 + qJD(5) * t64, -qJD(5) * t55 - qJD(6) * t350 + t108, t6 * qJD(1) + t7 * qJD(5) + t119 * t761 + ((t178 * t214 - t398 * t399 + t400 * t401) * m(7) + (t225 * t268 - t467 * t468 + t469 * t470) * m(6) + 0.4e1 * (t337 * (t524 * t582 + t525 * t581) - (t989 + t990) * t593 * t703) * t977 + (t84 + t83 + t990 * t518 + (t641 * t582 + (-t519 + t642) * t581) * t582) * t957 + (t86 + t85 + t989 * t519 + (t642 * t581 + (-t518 + t641) * t582) * t581) * t953) * qJD(4), t1 * qJD(1) + t7 * qJD(4) + t9 * qJD(6) + (t209 * t178 + t199 * t212 - t266 * t465 + t267 * t464 + t274 * t400 - t275 * t398 - t36 / 0.4e1) * t909 + (-t70 / 0.4e1 + t287 * t225 + t277 * t289 + t348 * t469 - t349 * t467 + (-t1081 * t581 + t321 * t582) * t578) * t912 + (t581 * t663 + t582 * t664 - t730) * t792 - t684 + (((t728 + t734) * t582 + (t729 - t733) * t581 + t777) * t608 + t1008 * t957 + t1007 * t953) * qJD(5), t19 * qJD(1) - t789 + t119 * t910 / 0.4e1 + t9 * qJD(5) - t343 * t908 / 0.4e1; t18 * qJD(2) + t790 + t2 * qJD(4) + t5 * qJD(5) + t11 * qJD(6) - t722 * t795 + (-t120 / 0.4e1 + t829 * t267 + t828 * t260) * t911 + (-t176 / 0.4e1 - t827 * t321) * t914 - t994 * t796, qJD(1) * t18 - qJD(4) * t64 + qJD(6) * t279, qJD(4) * t55 - qJD(6) * t344 + t147, t2 * qJD(1) + t4 * qJD(5) + t8 * qJD(6) + ((t732 + t1002) * t608 + (-t84 / 0.2e1 - t83 / 0.2e1 - t232 / 0.2e1 - t229 / 0.2e1) * t606 + t663) * t793 + ((t731 + t1001) * t608 + (-t86 / 0.2e1 - t85 / 0.2e1 - t227 / 0.2e1 + t231 / 0.2e1) * t606 - t664) * t794 + (-t94 / 0.4e1 + t143 * t178 + t195 * t400 - t196 * t398 + t199 * t214 + t266 * t401 - t267 * t399) * t910 + (-t145 / 0.4e1 + t210 * t225 + t268 * t277 + t271 * t469 - t272 * t467 + t1081 * t470 - t321 * t468) * t913 + t684, t5 * qJD(1) + t4 * qJD(4) + ((t1081 * t348 + t277 * t287 - t321 * t349) * m(6) + t1003 * t948) * qJD(5) + (t1007 * t957 + t1008 * t1012 + (t1037 * t581 + t1038 * t582) * t948) * t792 + ((t199 * t209 + t266 * t274 - t267 * t275) * qJD(5) + qJD(6) * t973) * m(7), t11 * qJD(1) + t8 * qJD(4) + t909 * t973 + (t360 * t860 + t516 * t446 - t512 * t447 - t301 / 0.4e1) * t908 + t668; (t1042 - t261 * t512 - t165 / 0.4e1) * t911 + t264 * qJD(2) + t21 * qJD(4) + t13 * qJD(5), qJD(1) * t264 - qJD(5) * t279, qJD(4) * t350 + qJD(5) * t344, t21 * qJD(1) + t789 + (-t512 * t399 + t516 * t401 - t119 / 0.4e1 + (-t178 * t608 + (-t214 + t672) * t606) * t605) * t910 + t10 * qJD(5) + t343 * t761, t13 * qJD(1) + t10 * qJD(4) + (t517 * t266 + t513 * t267 + t516 * t274 - t512 * t275 - t95 / 0.4e1 + (-t199 * t607 - t209 * t605) * t606) * t909 + t301 * t761 - t668 (t343 * qJD(4) / 0.4e1 + t301 * qJD(5) / 0.4e1) * m(7);];
Cq  = t14;
