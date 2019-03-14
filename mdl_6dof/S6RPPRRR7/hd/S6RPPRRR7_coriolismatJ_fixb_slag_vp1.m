% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:09
% EndTime: 2019-03-09 02:33:42
% DurationCPUTime: 26.12s
% Computational Cost: add. (106590->826), mult. (97576->1112), div. (0->0), fcn. (104407->10), ass. (0->515)
t672 = sin(qJ(1));
t674 = cos(qJ(1));
t666 = pkin(10) + qJ(4);
t653 = sin(t666);
t899 = pkin(4) * t653;
t669 = sin(pkin(10));
t900 = pkin(3) * t669;
t629 = t899 + t900;
t896 = -pkin(7) - qJ(3);
t771 = -pkin(8) + t896;
t777 = -t674 * t629 - t672 * t771;
t483 = t674 * (t672 * t896 + t674 * t900 + t777);
t648 = t674 * t771;
t651 = t674 * t896;
t502 = -t648 + t651 + (t629 - t900) * t672;
t655 = qJ(5) + t666;
t649 = sin(t655);
t673 = cos(qJ(6));
t831 = t673 * t674;
t671 = sin(qJ(6));
t833 = t672 * t671;
t596 = -t649 * t833 + t831;
t832 = t672 * t673;
t841 = t671 * t674;
t597 = t649 * t832 + t841;
t650 = cos(t655);
t846 = t650 * t672;
t447 = t597 * rSges(7,1) + t596 * rSges(7,2) - rSges(7,3) * t846;
t1005 = -pkin(9) * t846 + t447;
t852 = t649 * t672;
t794 = -pkin(5) * t852 - t1005;
t598 = t649 * t841 + t832;
t599 = t649 * t831 - t833;
t1002 = t599 * rSges(7,1) - t598 * rSges(7,2);
t844 = t650 * t674;
t448 = -rSges(7,3) * t844 + t1002;
t639 = pkin(9) * t844;
t851 = t649 * t674;
t799 = (-pkin(5) * t851 - t448 + t639) * t674;
t236 = t483 + (-t502 + t794) * t672 + t799;
t847 = t650 * t671;
t767 = rSges(7,2) * t847;
t845 = t650 * t673;
t768 = rSges(7,1) * t845;
t776 = rSges(7,3) * t851 + t674 * t768;
t496 = t674 * t767 - t776;
t775 = -pkin(5) * t844 - pkin(9) * t851;
t1024 = t496 + t775;
t495 = rSges(7,3) * t852 + (-t767 + t768) * t672;
t792 = -pkin(5) * t846 - pkin(9) * t852 - t495;
t310 = t1024 * t674 + t672 * t792;
t726 = rSges(7,1) * t673 - rSges(7,2) * t671;
t530 = rSges(7,3) * t649 + t650 * t726;
t511 = t672 * t530;
t610 = pkin(5) * t650 + pkin(9) * t649;
t435 = t672 * t610 + t511;
t654 = cos(t666);
t842 = t654 * t672;
t644 = pkin(4) * t842;
t414 = t644 + t435;
t789 = t530 + t610;
t898 = pkin(4) * t654;
t416 = (t789 + t898) * t674;
t893 = rSges(7,3) * t650;
t529 = -t649 * t726 + t893;
t510 = t672 * t529;
t897 = pkin(5) * t649;
t609 = pkin(9) * t650 - t897;
t434 = t672 * t609 + t510;
t436 = (-t529 - t609) * t674;
t106 = t236 * t310 + t414 * t434 - t416 * t436;
t574 = rSges(6,1) * t846 - rSges(6,2) * t852;
t608 = rSges(6,1) * t650 - rSges(6,2) * t649;
t576 = t608 * t674;
t1003 = t672 * t574 + t674 * t576;
t728 = rSges(6,1) * t649 + rSges(6,2) * t650;
t858 = t674 * t728;
t687 = -t672 * rSges(6,3) + t858;
t535 = t674 * t687;
t859 = t672 * t728;
t545 = rSges(6,3) * t674 + t859;
t304 = t483 - t535 + (-t502 - t545) * t672;
t532 = t608 * t672 + t644;
t731 = (t608 + t898) * t674;
t208 = -t1003 * t304 - t532 * t859 - t731 * t858;
t1054 = -m(6) * t208 - m(7) * t106;
t656 = t674 * qJ(2);
t1021 = -t672 * pkin(1) + t656;
t699 = t1021 - t777;
t1022 = (-t893 + t897) * t674 - t639 + t699 + t1002;
t824 = qJ(2) + t629;
t901 = pkin(1) * t674;
t352 = t901 - t648 + (t824 + t897) * t672 + t1005;
t470 = rSges(7,1) * t596 - rSges(7,2) * t597;
t471 = rSges(7,1) * t598 + rSges(7,2) * t599;
t880 = Icges(7,4) * t673;
t720 = -Icges(7,2) * t671 + t880;
t523 = Icges(7,6) * t649 + t650 * t720;
t881 = Icges(7,4) * t671;
t723 = Icges(7,1) * t673 - t881;
t525 = Icges(7,5) * t649 + t650 * t723;
t568 = (-Icges(7,2) * t673 - t881) * t650;
t571 = (-Icges(7,1) * t671 - t880) * t650;
t565 = (-Icges(7,5) * t671 - Icges(7,6) * t673) * t650;
t855 = t649 * t565;
t131 = t650 * (-(t525 / 0.2e1 + t568 / 0.2e1) * t671 + (t571 / 0.2e1 - t523 / 0.2e1) * t673) + m(7) * (-t1022 * t471 + t352 * t470) + t855 / 0.2e1;
t1053 = t131 * qJD(1);
t438 = Icges(7,5) * t597 + Icges(7,6) * t596 - Icges(7,3) * t846;
t882 = Icges(7,4) * t597;
t441 = Icges(7,2) * t596 - Icges(7,6) * t846 + t882;
t581 = Icges(7,4) * t596;
t444 = Icges(7,1) * t597 - Icges(7,5) * t846 + t581;
t258 = t438 * t844 + t598 * t441 - t599 * t444;
t1052 = t258 * t674;
t1051 = t672 * t258;
t583 = Icges(7,4) * t599;
t443 = Icges(7,2) * t598 + Icges(7,6) * t844 - t583;
t582 = Icges(7,4) * t598;
t445 = Icges(7,1) * t599 - Icges(7,5) * t844 - t582;
t1042 = t443 * t671 + t445 * t673;
t440 = -Icges(7,5) * t599 + Icges(7,6) * t598 + Icges(7,3) * t844;
t867 = t440 * t649;
t272 = t1042 * t650 - t867;
t963 = t672 / 0.2e1;
t962 = -t674 / 0.2e1;
t1009 = t672 * t731;
t518 = t574 + t644;
t402 = -t518 * t674 + t1009;
t404 = -t792 + t644;
t405 = (-t767 + t898) * t674 - t775 + t776;
t618 = rSges(5,1) * t654 - rSges(5,2) * t653;
t594 = t618 * t672;
t595 = t618 * t674;
t458 = -t594 * t674 + t672 * t595;
t979 = m(7) / 0.2e1;
t980 = m(6) / 0.2e1;
t981 = m(5) / 0.2e1;
t754 = (-t404 * t674 + t672 * t405) * t979 + t402 * t980 + t458 * t981;
t808 = (t414 * t674 - t416 * t672) * t979 + (t532 * t674 - t1009) * t980;
t89 = t808 - t754;
t1050 = qJD(1) * t89;
t1008 = t674 * t731;
t1041 = -m(5) / 0.2e1;
t667 = t672 ^ 2;
t668 = t674 ^ 2;
t646 = t667 + t668;
t752 = (-t414 * t672 - t416 * t674) * t979 + (-t532 * t672 - t1008) * t980 + t646 * t618 * t1041;
t459 = t672 * t594 + t595 * t674;
t753 = (t672 * t404 + t405 * t674) * t979 + (t672 * t518 + t1008) * t980 + t459 * t981;
t75 = t753 - t752;
t1049 = t75 * qJD(1);
t1016 = m(6) * t646;
t437 = t789 * t674;
t803 = (-t435 * t672 - t437 * t674) * t979 - t608 * t1016 / 0.2e1;
t804 = t1003 * t980 - t310 * t979;
t135 = t804 - t803;
t1048 = t135 * qJD(1);
t710 = -t598 * t443 - t599 * t445;
t800 = t596 * t441 + t597 * t444;
t1047 = t800 + (-t438 * t672 - t440 * t674) * t650 + t710;
t884 = Icges(6,4) * t649;
t721 = Icges(6,2) * t650 + t884;
t542 = -Icges(6,6) * t672 + t674 * t721;
t883 = Icges(6,4) * t650;
t724 = Icges(6,1) * t649 + t883;
t544 = -Icges(6,5) * t672 + t674 * t724;
t1006 = (t542 * t650 + t544 * t649) * t674;
t257 = -t440 * t846 + t596 * t443 - t445 * t597;
t1013 = t257 * t672;
t256 = -t438 * t846 + t800;
t156 = t256 * t674 + t1013;
t259 = t440 * t844 - t710;
t157 = t259 * t672 + t1052;
t1028 = t156 * t963 + t157 * t962;
t541 = Icges(6,6) * t674 + t672 * t721;
t630 = Icges(6,4) * t846;
t543 = Icges(6,1) * t852 + Icges(6,5) * t674 + t630;
t705 = -t541 * t650 - t543 * t649;
t1044 = t674 * t705;
t716 = Icges(6,5) * t649 + Icges(6,6) * t650;
t1032 = t672 * t716;
t354 = t674 * (Icges(6,3) * t674 + t1032) + t541 * t846 + t543 * t852;
t1030 = t674 * t716;
t540 = -Icges(6,3) * t672 + t1030;
t355 = -t674 * t540 - t542 * t846 - t544 * t852;
t357 = -t672 * t540 + t1006;
t819 = t259 + t1047;
t53 = t674 * t819 + t1013;
t816 = -t256 + t1047;
t54 = t672 * t816 - t1052;
t961 = t674 / 0.2e1;
t964 = -t672 / 0.2e1;
t688 = (t354 * t674 + t355 * t672) * t964 + ((t355 - t1044) * t672 + (t357 - t1006 + (t540 + t705) * t672 + t354) * t674 + t53) * t963 + (t357 * t672 + t667 * t540 + t54 + (t355 + (t540 - t705) * t674 + t1044) * t674) * t961 - t1028;
t1046 = t448 * t649 + t530 * t844;
t604 = -Icges(6,2) * t649 + t883;
t1043 = t604 + t724;
t1040 = -t649 / 0.2e1;
t966 = t649 / 0.2e1;
t1039 = -t650 / 0.2e1;
t1038 = t650 / 0.2e1;
t1037 = t672 / 0.4e1;
t1036 = t674 / 0.4e1;
t1035 = t157 + t54;
t715 = Icges(7,5) * t673 - Icges(7,6) * t671;
t521 = Icges(7,3) * t649 + t650 * t715;
t302 = t521 * t844 + t598 * t523 - t599 * t525;
t1034 = t302 * t649;
t885 = Icges(5,4) * t654;
t614 = -Icges(5,2) * t653 + t885;
t725 = Icges(5,1) * t653 + t885;
t1029 = (t725 / 0.2e1 + t614 / 0.2e1) * t654;
t491 = t523 * t672;
t493 = t525 * t672;
t711 = -t441 * t671 + t444 * t673;
t693 = t521 * t672 - t711;
t200 = (-t491 * t671 + t493 * t673 + t438) * t650 + t693 * t649;
t522 = Icges(7,6) * t650 - t649 * t720;
t524 = Icges(7,5) * t650 - t649 * t723;
t520 = Icges(7,3) * t650 - t649 * t715;
t863 = t525 * t673;
t864 = t523 * t671;
t707 = t863 - t864;
t691 = t520 - t707;
t865 = t521 * t649;
t997 = t650 * t691 - t865;
t223 = t522 * t596 + t524 * t597 - t672 * t997;
t1027 = t200 + t223;
t492 = t523 * t674;
t494 = t525 * t674;
t692 = -t521 * t674 + t1042;
t201 = (t492 * t671 - t494 * t673 + t440) * t650 + t692 * t649;
t224 = t598 * t522 - t599 * t524 + t674 * t997;
t1026 = t201 + t224;
t1025 = -t1022 * t672 + t674 * t352;
t425 = -t648 + (rSges(6,3) + pkin(1)) * t674 + (t728 + t824) * t672;
t419 = t674 * t425;
t424 = t687 + t699;
t732 = t419 * t728 - t424 * t859;
t809 = t1022 * t434 + t436 * t352;
t822 = (-t404 * t437 + t405 * t435 + t809) * t979 + (t402 * t608 + t732) * t980;
t823 = (-t1024 * t414 + t416 * t792 + t809) * t979 + (t532 * t576 - t574 * t731 + t732) * t980;
t1023 = t822 - t823;
t772 = qJD(4) + qJD(5);
t785 = t604 * t674 + t544;
t786 = -Icges(6,2) * t852 + t543 + t630;
t606 = Icges(6,1) * t650 - t884;
t787 = -t606 * t674 + t542;
t788 = -t606 * t672 + t541;
t1020 = -(t672 * t787 - t674 * t788) * t649 + (t672 * t785 - t674 * t786) * t650;
t563 = -Icges(5,5) * t672 + t674 * t725;
t781 = t614 * t674 + t563;
t643 = Icges(5,4) * t842;
t843 = t653 * t672;
t562 = Icges(5,1) * t843 + Icges(5,5) * t674 + t643;
t782 = -Icges(5,2) * t843 + t562 + t643;
t886 = Icges(5,4) * t653;
t722 = Icges(5,2) * t654 + t886;
t561 = -Icges(5,6) * t672 + t674 * t722;
t616 = Icges(5,1) * t654 - t886;
t783 = -t616 * t674 + t561;
t560 = Icges(5,6) * t674 + t672 * t722;
t784 = -t616 * t672 + t560;
t1019 = -(t672 * t783 - t674 * t784) * t653 + (t672 * t781 - t674 * t782) * t654;
t237 = (-t522 * t671 + t524 * t673 + t521) * t650 + t691 * t649;
t869 = t438 * t649;
t271 = t650 * t711 + t869;
t300 = -t521 * t846 + t523 * t596 + t525 * t597;
t336 = t650 * t707 + t865;
t1017 = t336 * t1038 + t237 * t966 + (t271 + t300) * t852 / 0.4e1 - (-t272 + t302) * t851 / 0.4e1;
t1012 = t257 * t674;
t1007 = (t654 * t561 + t653 * t563) * t674;
t664 = t672 * rSges(5,3);
t729 = rSges(5,1) * t653 + rSges(5,2) * t654;
t689 = t729 + t900;
t460 = t656 - t664 + (-pkin(1) + t896) * t672 + t689 * t674;
t461 = -t651 + (rSges(5,3) + pkin(1)) * t674 + (qJ(2) + t689) * t672;
t1004 = -t460 * t672 + t674 * t461;
t282 = t672 * t794 + t799;
t127 = t282 * t310 + t435 * t434 - t437 * t436;
t421 = -t545 * t672 - t535;
t274 = -t608 * t646 * t728 - t1003 * t421;
t1001 = m(6) * t274 + m(7) * t127;
t464 = Icges(7,5) * t596 - Icges(7,6) * t597;
t796 = -Icges(7,2) * t597 + t444 + t581;
t798 = -Icges(7,1) * t596 + t441 + t882;
t202 = -t464 * t846 + t596 * t796 - t597 * t798;
t465 = Icges(7,5) * t598 + Icges(7,6) * t599;
t795 = Icges(7,2) * t599 - t445 + t582;
t797 = -Icges(7,1) * t598 + t443 - t583;
t203 = -t465 * t846 + t596 * t795 - t597 * t797;
t97 = t202 * t674 + t203 * t672;
t204 = t464 * t844 + t598 * t796 + t599 * t798;
t205 = t465 * t844 + t598 * t795 + t599 * t797;
t98 = t204 * t674 + t205 * t672;
t895 = t961 * t97 + t963 * t98;
t450 = -t574 * t674 + t672 * t576;
t805 = (-t1024 * t672 + t674 * t792) * t979 + t450 * t980;
t996 = t650 * t692 - t867;
t995 = t650 * t693 - t869;
t994 = t649 * t1043 - t650 * (-t721 + t606);
t790 = t525 + t568;
t791 = t523 - t571;
t250 = -t565 * t846 + t596 * t790 - t597 * t791;
t251 = t565 * t844 + t598 * t790 + t599 * t791;
t216 = t465 * t649 + (-t671 * t795 - t673 * t797) * t650;
t873 = t216 * t672;
t215 = t464 * t649 + (-t671 * t796 - t673 * t798) * t650;
t874 = t215 * t674;
t735 = t874 / 0.4e1 + t873 / 0.4e1 + t251 * t1037 + t250 * t1036;
t295 = t300 * t649;
t36 = t295 + (-t672 * t819 + t1012) * t650;
t37 = -t1034 + (t674 * t816 + t1051) * t650;
t713 = t259 * t674 - t1051;
t114 = t650 * t713 + t1034;
t829 = t674 * t114;
t714 = -t256 * t672 + t1012;
t113 = t650 * t714 + t295;
t839 = t672 * t113;
t988 = t829 / 0.4e1 - t839 / 0.4e1 + t37 * t1036 + t36 * t1037;
t682 = -t522 * t847 / 0.2e1 + t524 * t845 / 0.2e1 + t521 * t1038 + (t863 + t606) * t1040 + t1043 * t1039 + (t864 + t520 + t721) * t966;
t986 = 0.2e1 * t646;
t985 = 0.4e1 * qJD(1);
t984 = 2 * qJD(4);
t982 = 2 * qJD(5);
t708 = t447 * t674 - t448 * t672;
t324 = t708 * t650;
t361 = t447 * t649 + t511 * t650;
t755 = t1046 * t434 - t324 * t310 + t361 * t436;
t242 = (-t495 * t674 - t496 * t672) * t650 + t708 * t649;
t275 = (t447 + t510) * t650 + (t495 - t511) * t649;
t276 = (t529 * t674 + t448) * t650 + (-t530 * t674 - t496) * t649;
t756 = t242 * t236 - t275 * t416 + t276 * t414;
t978 = m(7) * (t755 + t756);
t40 = t242 * t282 - t275 * t437 + t276 * t435 + t755;
t977 = m(7) * t40;
t376 = t672 * t470 - t471 * t674;
t575 = (-rSges(7,1) * t671 - rSges(7,2) * t673) * t650;
t862 = t575 * t672;
t811 = -t282 * t376 + t435 * t862;
t813 = -t236 * t376 + t414 * t862;
t861 = t575 * t674;
t973 = m(7) * ((t416 + t437) * t861 + t811 + t813);
t812 = t1022 * t276 + t275 * t352;
t972 = m(7) * (t1046 * t405 + t361 * t404 + t812);
t971 = m(7) * (t1046 * t276 - t242 * t324 + t275 * t361);
t970 = m(7) * (-t1024 * t1046 - t361 * t792 + t812);
t968 = t646 / 0.2e1;
t960 = m(3) * ((rSges(3,3) * t674 + t1021) * t674 + (t901 + (rSges(3,3) + qJ(2)) * t672) * t672);
t730 = rSges(4,1) * t669 + rSges(4,2) * cos(pkin(10));
t770 = rSges(4,3) + pkin(1) + qJ(3);
t508 = -t672 * t770 + t674 * t730 + t656;
t509 = t770 * t674 + (qJ(2) + t730) * t672;
t959 = m(4) * (-t508 * t672 + t674 * t509);
t958 = m(4) * (t508 * t674 + t672 * t509);
t957 = m(5) * (t460 * t595 + t461 * t594);
t956 = m(5) * t1004;
t955 = m(5) * (t460 * t674 + t672 * t461);
t945 = m(6) * (t424 * t731 + t425 * t518);
t944 = m(6) * (t424 * t576 + t425 * t574);
t943 = m(6) * (-t424 * t672 + t419);
t942 = m(6) * (t424 * t674 + t672 * t425);
t733 = t1025 * t575;
t927 = m(7) * (-t414 * t471 - t416 * t470 - t733);
t926 = m(7) * (-t435 * t471 - t437 * t470 - t733);
t923 = m(7) * (t275 * t672 + t276 * t674);
t922 = m(7) * (t1022 * t405 + t352 * t404);
t921 = m(7) * (-t1022 * t1024 - t352 * t792);
t920 = m(7) * t1025;
t919 = m(7) * (t1022 * t674 + t672 * t352);
t918 = m(7) * (-t1046 * t672 + t361 * t674);
t917 = m(7) * (t1046 * t674 + t361 * t672);
t908 = m(7) * (t434 * t674 + t436 * t672);
t906 = m(7) * (t435 * t674 - t437 * t672);
t374 = -t470 * t674 - t672 * t471;
t903 = m(7) * t374;
t902 = m(7) * t376;
t894 = m(7) * qJD(6);
t876 = t200 * t674;
t875 = t201 * t672;
t871 = t272 * t674;
t690 = m(7) * (-t275 * t674 + t276 * t672);
t698 = m(7) * t575 * t968;
t150 = t690 / 0.2e1 + t698;
t280 = m(7) * (t434 * t672 - t436 * t674) - t728 * t1016;
t821 = t280 * qJD(5) + t150 * qJD(6);
t765 = -t923 / 0.2e1;
t151 = -t690 / 0.2e1 + t698;
t774 = t151 * qJD(2);
t820 = qJD(3) * t765 + t774;
t764 = t923 / 0.2e1;
t815 = qJD(5) * t908 + qJD(6) * t764;
t814 = t772 * t764;
t369 = (-m(7) / 0.2e1 - m(6) / 0.2e1 + t1041 - m(4) / 0.2e1) * t986;
t773 = t369 * qJD(1);
t769 = t973 / 0.2e1 + t895;
t712 = -t271 * t672 - t871;
t126 = t336 * t649 + t650 * t712;
t185 = t491 * t596 + t493 * t597 - t672 * t995;
t186 = -t492 * t596 - t494 * t597 - t672 * t996;
t27 = (-t185 * t672 + t186 * t674 + t300) * t650 + (t223 - t714) * t649;
t187 = t598 * t491 - t599 * t493 + t674 * t995;
t188 = -t598 * t492 + t599 * t494 + t674 * t996;
t28 = (-t187 * t672 + t188 * t674 + t302) * t650 + (t224 - t713) * t649;
t38 = (-t200 * t672 + t201 * t674 + t336) * t650 + (t237 - t712) * t649;
t13 = t971 + (t27 * t964 + t28 * t961 + t126 / 0.2e1) * t650 + (t839 / 0.2e1 - t829 / 0.2e1 + t38 / 0.2e1) * t649;
t763 = qJD(3) * t764 + t13 * qJD(6) - t774;
t758 = -t113 / 0.2e1 + t36 / 0.2e1;
t757 = -t37 / 0.2e1 - t114 / 0.2e1;
t718 = Icges(5,5) * t653 + Icges(5,6) * t654;
t558 = Icges(5,3) * t674 + t672 * t718;
t370 = t674 * t558 + t560 * t842 + t562 * t843;
t559 = -Icges(5,3) * t672 + t674 * t718;
t371 = -t674 * t559 - t561 * t842 - t563 * t843;
t747 = -t846 / 0.2e1;
t746 = -t846 / 0.4e1;
t745 = t846 / 0.4e1;
t744 = -t844 / 0.4e1;
t743 = t844 / 0.2e1;
t742 = t844 / 0.4e1;
t737 = t646 * t729;
t717 = Icges(6,5) * t650 - Icges(6,6) * t649;
t566 = t672 * t717;
t567 = t717 * t674;
t85 = t185 * t674 + t186 * t672;
t86 = t187 * t674 + t188 * t672;
t736 = (t86 - t667 * t567 + (t672 * t566 + t1020) * t674) * t963 + (t85 + t668 * t566 + (-t674 * t567 - t1020) * t672) * t961;
t734 = t646 * t898;
t719 = Icges(5,5) * t654 - Icges(5,6) * t653;
t531 = (-t728 - t899) * t672;
t645 = t674 * t899;
t533 = t645 + t858;
t706 = t531 * t672 - t533 * t674;
t703 = t542 * t649 - t544 * t650;
t702 = -t654 * t560 - t653 * t562;
t696 = t1046 * t862 + t324 * t376 - t361 * t861;
t686 = t1035 * t746 + t156 * t744 + t53 * t742 + t988;
t685 = t27 * t961 + t28 * t963 + t86 * t743 + t85 * t747 + (t271 * t674 - t272 * t672) * t1038 + (t875 + t876) * t966 - t895 + t1028 * t649;
t683 = t1026 * t742 + t1027 * t746 + t1017;
t68 = t250 * t649 + (-t202 * t672 + t203 * t674) * t650;
t69 = t251 * t649 + (-t204 * t672 + t205 * t674) * t650;
t681 = t126 * t1039 + t27 * t846 / 0.2e1 - t28 * t844 / 0.2e1 - t971 + t69 * t963 + t68 * t961 + t97 * t747 + t98 * t743 + (t839 + t38) * t1040 + (t829 + t873 + t874) * t966;
t679 = t876 / 0.2e1 + t875 / 0.2e1 + (t649 * t785 + t650 * t787 + t674 * t994 - t1032 + t224) * t963 + (-t649 * t786 - t650 * t788 - t672 * t994 - t1030 + t223) * t961 - t688;
t678 = -t871 / 0.2e1 - t682 + t703 * t962 + (t703 + t272) * t961;
t677 = t1035 * t745 + t156 * t742 + t53 * t744 + t683 + t735 - t988;
t676 = t1026 * t744 + t1027 * t745 - t1017 + t686 + t735;
t675 = t683 + t686 - t735;
t587 = t719 * t674;
t586 = t672 * t719;
t536 = t672 * t558;
t415 = t645 + t436;
t413 = -pkin(4) * t843 + t434;
t401 = -t1003 - t734;
t379 = -t649 * t471 + t575 * t844;
t378 = t470 * t649 + t575 * t846;
t373 = -t672 * t559 + t1007;
t372 = t674 * t702 + t536;
t368 = (-m(7) / 0.4e1 - m(6) / 0.4e1 - m(5) / 0.4e1 - m(4) / 0.4e1) * t986 + (m(4) + m(5) + m(6) + m(7)) * t968;
t367 = t902 / 0.2e1;
t366 = t903 / 0.2e1;
t353 = t374 * t650;
t333 = t906 / 0.2e1;
t292 = -t734 + t310;
t277 = (t855 + (-t671 * t790 - t673 * t791) * t650) * t649;
t267 = t372 * t674 + t373 * t672;
t266 = t370 * t674 + t371 * t672;
t255 = t917 / 0.2e1;
t254 = t918 / 0.2e1;
t176 = qJD(6) * t765;
t167 = t437 * t861 + t811;
t160 = t333 - t805;
t159 = -t906 / 0.2e1 + t805;
t158 = t333 + t805;
t141 = t151 * qJD(6);
t138 = t926 / 0.2e1;
t134 = t803 + t804;
t132 = t927 / 0.2e1;
t128 = t416 * t861 + t813;
t124 = t920 + t943 + t956 + t959;
t117 = t667 * t559 + (t371 - t536 + (t559 - t702) * t674) * t674;
t116 = (-t372 + t536 + t371) * t672 + (t373 - t1007 + (t559 + t702) * t672 + t370) * t674;
t115 = t919 + t942 + t955 + t958 + t960;
t87 = t754 + t808;
t84 = t255 - t903 / 0.2e1;
t83 = t254 - t902 / 0.2e1;
t82 = t367 + t254;
t81 = t367 - t918 / 0.2e1;
t80 = t366 + t255;
t79 = t366 - t917 / 0.2e1;
t74 = t752 + t753;
t70 = t682 + t921 + t944;
t64 = t970 / 0.2e1;
t60 = t972 / 0.2e1;
t57 = -t1029 + t957 + t945 + t922 + t682 + (-t616 / 0.2e1 + t722 / 0.2e1) * t653;
t39 = t977 / 0.2e1;
t29 = t978 / 0.2e1;
t22 = m(7) * t167 + t895;
t21 = m(7) * t128 + t895;
t20 = t736 + t1001;
t19 = t20 * qJD(5);
t18 = t736 - t1054;
t17 = (t672 * t757 + t674 * t758) * t650;
t16 = t39 - t978 / 0.2e1 + t769;
t15 = t29 - t977 / 0.2e1 + t769;
t11 = (t117 / 0.2e1 + t267 / 0.2e1) * t674 + (-t266 / 0.2e1 + t116 / 0.2e1) * t672 + t688;
t10 = t688 - t1023;
t9 = t688 + t1023;
t8 = t679 + t822 + t823;
t7 = t29 + t39 - t973 / 0.2e1 + t685;
t6 = t677 + t138 + t64;
t5 = t675 - t926 / 0.2e1 + t64;
t4 = t138 + t676 - t970 / 0.2e1;
t3 = t677 + t132 + t60;
t2 = t675 - t927 / 0.2e1 + t60;
t1 = t676 - t972 / 0.2e1 + t132;
t12 = [t115 * qJD(2) + t124 * qJD(3) + t57 * qJD(4) + t70 * qJD(5) + t131 * qJD(6), qJD(1) * t115 + qJD(3) * t368 + qJD(4) * t87 + qJD(5) * t158 + qJD(6) * t80, qJD(1) * t124 + qJD(2) * t368 + qJD(4) * t74 + qJD(5) * t134 + qJD(6) * t82, t57 * qJD(1) + t87 * qJD(2) + t74 * qJD(3) + t8 * qJD(5) + t3 * qJD(6) + ((t1004 * t729 + t458 * t618) * t981 + (t424 * t531 + t425 * t533 + (-t518 + t532) * t731) * t980 + (t1022 * t413 + t352 * t415 - t404 * t416 + t405 * t414) * t979) * t984 + (t679 + t116 * t964 + (-t653 * t782 - t654 * t784) * t961 + (t653 * t781 + t654 * t783 + t266) * t963 + (t117 + t267) * t962 - (t667 / 0.2e1 + t668 / 0.2e1) * t718) * qJD(4), t70 * qJD(1) + t158 * qJD(2) + t134 * qJD(3) + t8 * qJD(4) + t679 * qJD(5) + t6 * qJD(6) + ((-t1024 * t435 + t437 * t792 + t809) * t979 + (t450 * t608 + t732) * t980) * t982, t1053 + t80 * qJD(2) + t82 * qJD(3) + t3 * qJD(4) + t6 * qJD(5) + (t277 + m(7) * (t1022 * t379 - t1046 * t471 + t352 * t378 + t361 * t470) + ((t216 / 0.2e1 + t251 / 0.2e1 - t758) * t674 + (-t250 / 0.2e1 - t215 / 0.2e1 - t757) * t672) * t650) * qJD(6); t369 * qJD(3) - t89 * qJD(4) + t159 * qJD(5) + t79 * qJD(6) + (-t919 / 0.4e1 - t942 / 0.4e1 - t955 / 0.4e1 - t960 / 0.4e1 - t958 / 0.4e1) * t985, 0, t773, -t1050 + (t706 * t980 + (t413 * t672 - t415 * t674) * t979 - t737 * t981) * t984 + t821, qJD(1) * t159 + qJD(4) * t280 + t821, t79 * qJD(1) + (-t378 * t674 + t379 * t672) * t894 + t772 * t150; -t369 * qJD(2) + t75 * qJD(4) + t135 * qJD(5) + t81 * qJD(6) + (-t920 / 0.4e1 - t943 / 0.4e1 - t956 / 0.4e1 - t959 / 0.4e1) * t985, -t773, 0, t1049 + ((t531 * t674 + t533 * t672) * t980 + (t413 * t674 + t415 * t672) * t979) * t984 + t815, qJD(4) * t908 + t1048 + t815, t81 * qJD(1) + (t378 * t672 + t379 * t674) * t894 + t814; (t678 + (-t722 + t616) * t653 / 0.2e1 + t1029) * qJD(1) + t89 * qJD(2) - t75 * qJD(3) + t11 * qJD(4) + t10 * qJD(5) + t1 * qJD(6) + (-t957 / 0.4e1 - t945 / 0.4e1 - t922 / 0.4e1) * t985, t141 + t1050, t176 - t1049, t11 * qJD(1) + ((-t667 * t587 + (t672 * t586 + t1019) * t674) * t963 + m(7) * (t236 * t292 + t413 * t414 - t415 * t416) + m(6) * (t304 * t401 + t531 * t532 - t533 * t731) + m(5) * (-(-t674 * (t674 * t729 - t664) + (-t674 * rSges(5,3) - t672 * t729) * t672) * t459 - t618 * t737) + (t668 * t586 + (-t674 * t587 - t1019) * t672) * t961 + t736) * qJD(4) + t18 * qJD(5) + t21 * qJD(6), t10 * qJD(1) + t18 * qJD(4) + t15 * qJD(6) + ((t127 + t106) * t979 + (t208 + t274) * t980) * t982 + (t736 - t1001) * qJD(5), t1 * qJD(1) + t21 * qJD(4) + t15 * qJD(5) + (t681 + m(7) * (t353 * t236 - t378 * t416 + t379 * t414 + t696)) * qJD(6) + t820; t678 * qJD(1) + t160 * qJD(2) - t135 * qJD(3) + t9 * qJD(4) + t688 * qJD(5) + t4 * qJD(6) + (-t921 / 0.4e1 - t944 / 0.4e1) * t985, qJD(1) * t160 + t141, t176 - t1048, t9 * qJD(1) + t19 + t16 * qJD(6) + ((t282 * t292 + t413 * t435 - t415 * t437 + t106) * t979 + (t421 * t401 + t608 * t706 + t208) * t980) * t984 + (t736 + t1054) * qJD(4), qJD(1) * t688 + qJD(4) * t20 + qJD(6) * t22 + t19, t4 * qJD(1) + t16 * qJD(4) + t22 * qJD(5) + (t681 + m(7) * (t353 * t282 - t378 * t437 + t379 * t435 + t696)) * qJD(6) + t820; t84 * qJD(2) + t83 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t17 * qJD(6) - t1053, qJD(1) * t84 - t151 * t772, qJD(1) * t83 + t814, t2 * qJD(1) + ((t1046 * t413 - t292 * t324 + t361 * t415 - t128 + t756) * m(7) + t685) * qJD(4) + t7 * qJD(5) + t763, t5 * qJD(1) + t7 * qJD(4) + ((t40 - t167) * m(7) + t685) * qJD(5) + t763, t17 * qJD(1) + (t277 * t966 + m(7) * (t1046 * t379 - t324 * t353 + t361 * t378) + (t68 * t964 + t69 * t961 + (-t215 * t672 + t216 * t674) * t966) * t650) * qJD(6) + t772 * t13;];
Cq  = t12;