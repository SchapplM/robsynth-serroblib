% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:05
% EndTime: 2019-12-31 19:43:54
% DurationCPUTime: 39.95s
% Computational Cost: add. (39473->922), mult. (100870->1281), div. (0->0), fcn. (119855->8), ass. (0->513)
t565 = sin(qJ(1));
t567 = cos(qJ(2));
t564 = sin(qJ(2));
t561 = sin(pkin(8));
t562 = cos(pkin(8));
t563 = sin(qJ(5));
t566 = cos(qJ(5));
t602 = t561 * t566 - t562 * t563;
t887 = t602 * t564;
t601 = t561 * t563 + t562 * t566;
t888 = t601 * t564;
t370 = rSges(6,1) * t888 + rSges(6,2) * t887 + rSges(6,3) * t567;
t791 = pkin(4) * t562;
t676 = t564 * t791;
t519 = pkin(7) * t567 + t676;
t529 = pkin(2) * t564 - qJ(3) * t567;
t757 = qJ(4) * t561;
t626 = pkin(3) * t562 + t757;
t886 = t626 * t564;
t688 = t529 + t886;
t638 = t370 + t519 + t688;
t272 = t638 * t565;
t568 = cos(qJ(1));
t274 = t638 * t568;
t628 = rSges(5,1) * t562 + rSges(5,3) * t561;
t885 = t628 * t564;
t910 = rSges(5,2) * t567 - t885;
t659 = t688 - t910;
t348 = t659 * t565;
t350 = t659 * t568;
t787 = rSges(4,1) * t562;
t629 = -rSges(4,2) * t561 + t787;
t595 = t629 * t564;
t911 = rSges(4,3) * t567 - t595;
t697 = t529 - t911;
t422 = t697 * t565;
t424 = t697 * t568;
t557 = t568 * pkin(6);
t792 = pkin(2) * t567;
t656 = pkin(1) + t792;
t777 = rSges(4,3) + qJ(3);
t869 = t777 * t564 + t656;
t735 = t567 * t561;
t503 = t562 * t568 + t565 * t735;
t732 = t568 * t561;
t737 = t565 * t567;
t504 = t562 * t737 - t732;
t876 = -t504 * rSges(4,1) + t503 * rSges(4,2);
t318 = -t565 * t869 + t557 + t876;
t738 = t565 * t562;
t505 = t567 * t732 - t738;
t734 = t567 * t568;
t739 = t565 * t561;
t506 = t562 * t734 + t739;
t630 = t506 * rSges(4,1) - t505 * rSges(4,2);
t790 = t565 * pkin(6);
t319 = t568 * t869 + t630 + t790;
t713 = t318 * t734 + t319 * t737;
t495 = t503 * qJ(4);
t498 = t503 * rSges(5,3);
t778 = rSges(5,2) + qJ(3);
t868 = t778 * t564 + t656;
t575 = -t565 * t868 + t557;
t826 = rSges(5,1) + pkin(3);
t268 = -t826 * t504 - t495 - t498 + t575;
t496 = t505 * qJ(4);
t645 = t496 + t790;
t782 = t505 * rSges(5,3);
t269 = t826 * t506 + t568 * t868 + t645 + t782;
t722 = t268 * t734 + t269 * t737;
t855 = pkin(3) + pkin(4);
t532 = qJ(3) * t564 + t792;
t599 = pkin(1) + t532;
t416 = t503 * t566 - t504 * t563;
t417 = t503 * t563 + t504 * t566;
t742 = t564 * t565;
t287 = -rSges(6,1) * t417 - rSges(6,2) * t416 + rSges(6,3) * t742;
t948 = -pkin(7) * t742 - t287;
t935 = -t565 * t599 + t557 - t948;
t221 = -t855 * t504 - t495 + t935;
t420 = t505 * t566 - t506 * t563;
t421 = t505 * t563 + t506 * t566;
t740 = t564 * t568;
t288 = t421 * rSges(6,1) + t420 * rSges(6,2) - rSges(6,3) * t740;
t882 = -pkin(7) * t740 + t288;
t222 = t855 * t506 + t599 * t568 + t645 + t882;
t730 = t221 * t734 + t222 * t737;
t857 = m(6) / 0.2e1;
t858 = m(5) / 0.2e1;
t860 = m(4) / 0.2e1;
t660 = (-t272 * t740 + t274 * t742 + t730) * t857 + (-t348 * t740 + t350 * t742 + t722) * t858 + (-t422 * t740 + t424 * t742 + t713) * t860;
t537 = qJ(3) * t734;
t453 = t568 * t887;
t454 = t568 * t888;
t700 = -t454 * rSges(6,1) - t453 * rSges(6,2);
t825 = -rSges(6,3) - pkin(7);
t260 = t537 + (t825 * t567 + (-t855 * t562 - pkin(2) - t757) * t564) * t568 + t700;
t451 = t565 * t887;
t452 = t565 * t888;
t627 = rSges(6,1) * t452 + rSges(6,2) * t451;
t548 = pkin(2) * t742;
t685 = (-pkin(3) * t738 - qJ(4) * t739) * t564;
t657 = t548 - t685;
t261 = (t676 + (-qJ(3) - t825) * t567) * t565 + t627 + t657;
t546 = rSges(5,2) * t734;
t321 = t537 + t546 + (-pkin(2) - t826 * t562 + (-rSges(5,3) - qJ(4)) * t561) * t740;
t322 = (-t778 * t567 + t885) * t565 + t657;
t401 = t548 + (-t777 * t567 + t595) * t565;
t686 = t564 * rSges(4,2) * t732 + rSges(4,3) * t734;
t402 = t537 + (-pkin(2) - t787) * t740 + t686;
t669 = ((t321 * t565 + t322 * t568) * t564 + t722) * t858 + ((t401 * t568 + t402 * t565) * t564 + t713) * t860 + ((t260 * t565 + t261 * t568) * t564 + t730) * t857;
t7 = t669 - t660;
t956 = t7 * qJD(1);
t743 = t561 * t564;
t585 = (-t268 * t568 - t269 * t565) * t743;
t613 = -t221 * t568 - t222 * t565;
t586 = t613 * t743;
t774 = (-t505 * t348 + t503 * t350 + t585) * t858 + (-t505 * t272 + t503 * t274 + t586) * t857;
t775 = (t503 * t321 + t505 * t322 + t585) * t858 + (t503 * t260 + t505 * t261 + t586) * t857;
t12 = t775 - t774;
t955 = t12 * qJD(1);
t364 = Icges(6,5) * t888 + Icges(6,6) * t887 + Icges(6,3) * t567;
t487 = t602 * t567;
t488 = t601 * t567;
t365 = Icges(6,5) * t488 + Icges(6,6) * t487 - Icges(6,3) * t564;
t768 = Icges(3,4) * t564;
t524 = Icges(3,2) * t567 + t768;
t527 = Icges(3,1) * t567 - t768;
t555 = Icges(3,4) * t567;
t618 = Icges(4,5) * t562 - Icges(4,6) * t561;
t620 = Icges(5,4) * t562 + Icges(5,6) * t561;
t622 = Icges(5,1) * t562 + Icges(5,5) * t561;
t623 = Icges(4,1) * t562 - Icges(4,4) * t561;
t901 = -t567 / 0.2e1;
t903 = -t564 / 0.2e1;
t945 = Icges(5,4) + Icges(4,5);
t647 = t945 * t903 + (t622 + t623) * t901;
t617 = Icges(5,5) * t562 + Icges(5,3) * t561;
t621 = Icges(4,4) * t562 - Icges(4,2) * t561;
t829 = t567 / 0.2e1;
t902 = t564 / 0.2e1;
t649 = Icges(4,6) * t902 + Icges(5,6) * t903 + t617 * t901 + t621 * t829;
t461 = Icges(6,4) * t887;
t368 = Icges(6,1) * t888 + Icges(6,5) * t567 + t461;
t744 = t488 * t368;
t764 = Icges(6,4) * t888;
t366 = Icges(6,2) * t887 + Icges(6,6) * t567 + t764;
t745 = t487 * t366;
t369 = Icges(6,1) * t488 + Icges(6,4) * t487 - Icges(6,5) * t564;
t746 = t888 * t369;
t367 = Icges(6,4) * t488 + Icges(6,2) * t487 - Icges(6,6) * t564;
t747 = t887 * t367;
t758 = Icges(4,3) * t567;
t759 = Icges(5,6) * t567;
t760 = Icges(4,6) * t567;
t761 = Icges(5,2) * t567;
t762 = Icges(3,2) * t564;
t763 = Icges(4,5) * t567;
t767 = Icges(5,4) * t567;
t769 = Icges(3,1) * t564;
t889 = t564 * t623;
t890 = t564 * t622;
t891 = t564 * t621;
t892 = t564 * t620;
t893 = t564 * t618;
t894 = t564 * t617;
t954 = t564 * (t561 * t649 + t562 * t647 + t364 / 0.2e1 - t893 / 0.2e1 + t758 / 0.2e1 - t892 / 0.2e1 + t761 / 0.2e1 + t524 / 0.2e1 - t527 / 0.2e1) + (t561 * (t891 / 0.2e1 - t760 / 0.2e1 - t894 / 0.2e1 + t759 / 0.2e1) - t562 * (t890 / 0.2e1 - t767 / 0.2e1 + t889 / 0.2e1 - t763 / 0.2e1) - t365 / 0.2e1 - t555 + t762 / 0.2e1 - t769 / 0.2e1 + (Icges(4,3) + Icges(5,2)) * t902 + (t618 + t620) * t829) * t567 - t747 / 0.2e1 - t746 / 0.2e1 - t745 / 0.2e1 - t744 / 0.2e1;
t501 = t504 * pkin(4);
t642 = -t504 * pkin(3) - t495;
t223 = -t501 + t642 + t935;
t729 = t221 - t223;
t952 = m(6) * t729;
t951 = Icges(4,1) + Icges(5,1);
t946 = -Icges(4,4) + Icges(5,5);
t950 = Icges(4,2) + Icges(5,3);
t949 = -Icges(4,6) + Icges(5,6);
t277 = Icges(6,5) * t417 + Icges(6,6) * t416 - Icges(6,3) * t742;
t279 = Icges(6,5) * t421 + Icges(6,6) * t420 - Icges(6,3) * t740;
t765 = Icges(6,4) * t421;
t282 = Icges(6,2) * t420 - Icges(6,6) * t740 + t765;
t408 = Icges(6,4) * t420;
t285 = Icges(6,1) * t421 - Icges(6,5) * t740 + t408;
t724 = t420 * t282 + t421 * t285;
t766 = Icges(6,4) * t417;
t281 = -Icges(6,2) * t416 + Icges(6,6) * t742 - t766;
t407 = Icges(6,4) * t416;
t284 = -Icges(6,1) * t417 + Icges(6,5) * t742 - t407;
t905 = -t416 * t281 - t417 * t284;
t944 = -t905 + t724 + (t277 * t565 - t279 * t568) * t564;
t725 = -t420 * t281 - t421 * t284;
t726 = t416 * t282 + t417 * t285;
t943 = t564 * (t277 * t568 + t279 * t565) - t726 - t725;
t914 = t946 * t503 + t951 * t504 + t945 * t742;
t915 = t950 * t503 + t946 * t504 + t949 * t742;
t924 = t915 * t503 + t914 * t504;
t942 = t287 * t567 - t370 * t742;
t754 = t277 * t567;
t173 = t281 * t887 + t284 * t888 - t754;
t197 = -t364 * t742 + t366 * t416 + t368 * t417;
t938 = t197 * t567;
t929 = t950 * t505 + t946 * t506 + t949 * t740;
t381 = Icges(4,5) * t504 - Icges(4,6) * t503 + Icges(4,3) * t742;
t384 = Icges(5,4) * t504 + Icges(5,2) * t742 + Icges(5,6) * t503;
t937 = t381 + t384;
t383 = Icges(4,5) * t506 - Icges(4,6) * t505 + Icges(4,3) * t740;
t386 = Icges(5,4) * t506 + Icges(5,2) * t740 + Icges(5,6) * t505;
t936 = t383 + t386;
t927 = t946 * t505 + t951 * t506 + t945 * t740;
t923 = m(5) + m(6);
t932 = t923 / 0.2e1;
t375 = t503 * t740 - t505 * t742;
t677 = m(5) / 0.4e1 + m(6) / 0.4e1;
t904 = -0.2e1 * t677;
t931 = t375 * t904;
t930 = t923 * t375;
t559 = t565 ^ 2;
t560 = t568 ^ 2;
t682 = t559 + t560;
t926 = t759 - t894 - t760 + t891;
t925 = t763 - t889 + t767 - t890;
t689 = -t504 * rSges(5,1) - t498;
t270 = t575 + t642 + t689;
t720 = t268 - t270;
t920 = t720 * t858;
t919 = t729 * t857;
t917 = t929 * t505 + t927 * t506 + t936 * t740;
t480 = Icges(3,4) * t737 - Icges(3,2) * t742 - Icges(3,6) * t568;
t481 = Icges(3,6) * t565 + (t555 - t762) * t568;
t483 = Icges(3,5) * t565 + t527 * t568;
t442 = t483 * t737;
t523 = Icges(3,5) * t567 - Icges(3,6) * t564;
t895 = t523 * t568;
t479 = Icges(3,3) * t565 + t895;
t641 = t479 * t568 - t442;
t478 = Icges(3,5) * t737 - Icges(3,6) * t742 - Icges(3,3) * t568;
t540 = Icges(3,4) * t742;
t482 = Icges(3,1) * t737 - Icges(3,5) * t568 - t540;
t702 = -t565 * t478 - t482 * t734;
t916 = -t480 * t740 - t481 * t742 - t641 - t702;
t298 = rSges(6,1) * t416 - rSges(6,2) * t417;
t299 = rSges(6,1) * t420 - rSges(6,2) * t421;
t405 = Icges(6,1) * t887 - t764;
t403 = Icges(6,5) * t887 - Icges(6,6) * t888;
t736 = t567 * t403;
t906 = (t405 / 0.2e1 - t366 / 0.2e1) * t888 + m(6) * (-t221 * t298 + t222 * t299) + t736 / 0.2e1;
t832 = -t565 / 0.2e1;
t831 = t565 / 0.2e1;
t828 = -t568 / 0.2e1;
t404 = -Icges(6,2) * t888 + t461;
t704 = t368 + t404;
t896 = t887 * t704;
t741 = t564 * t567;
t684 = t682 * t741;
t884 = (m(4) / 0.4e1 + t677) * (t684 - t741);
t415 = t565 * t503 + t505 * t568;
t883 = t677 * (-t415 + t735) * t743;
t518 = t682 * t564;
t880 = t926 * t565;
t879 = t926 * t568;
t878 = t925 * t565;
t877 = t925 * t568;
t576 = t758 - t893;
t430 = t576 * t565;
t431 = t576 * t568;
t579 = t761 - t892;
t432 = t579 * t565;
t433 = t579 * t568;
t511 = -Icges(3,2) * t737 - t540;
t512 = t524 * t568;
t624 = -t555 - t769;
t513 = t624 * t565;
t514 = t624 * t568;
t865 = ((t561 * t915 + t562 * t914 - t430 - t432 + t482 + t511) * t568 + (-t561 * t929 - t562 * t927 + t431 + t433 - t483 + t512) * t565) * t564 + ((t480 - t513 - t937) * t568 + (-t481 + t514 + t936) * t565) * t567;
t864 = 0.2e1 * qJD(1);
t863 = 0.4e1 * qJD(1);
t862 = 2 * qJD(2);
t199 = -t364 * t740 + t420 * t366 + t421 * t368;
t196 = t199 * t567;
t152 = -t277 * t740 + t725;
t153 = -t279 * t740 + t724;
t615 = -t565 * t152 - t153 * t568;
t64 = t564 * t615 + t196;
t856 = t64 / 0.2e1;
t346 = -rSges(6,3) * t737 - t627;
t347 = -rSges(6,3) * t734 + t700;
t608 = t287 * t568 + t288 * t565;
t169 = t608 * t567 + (-t346 * t568 + t347 * t565) * t564;
t371 = rSges(6,1) * t488 + rSges(6,2) * t487 - rSges(6,3) * t564;
t185 = (-t370 * t565 - t346) * t567 + (-t371 * t565 - t287) * t564;
t186 = (t370 * t568 + t347) * t567 + (t371 * t568 - t288) * t564;
t209 = t608 * t564;
t242 = t567 * t288 + t370 * t740;
t612 = -t242 * t565 - t568 * t942;
t852 = m(6) * (t185 * t505 + t186 * t503 + (t209 * t567 + (t169 + t612) * t564) * t561);
t851 = m(6) * (t185 * t221 + t186 * t222 + t242 * t260 + t261 * t942);
t728 = t242 * t737 + t734 * t942;
t850 = m(6) * (-t169 * t567 + (t185 * t568 + t186 * t565 + t209) * t564 + t728);
t849 = m(6) * (t169 * t209 + t185 * t942 + t186 * t242);
t406 = rSges(6,1) * t887 - rSges(6,2) * t888;
t844 = m(6) * (-t272 * t299 + t274 * t298 + t406 * t613);
t340 = -Icges(6,5) * t452 - Icges(6,6) * t451 - Icges(6,3) * t737;
t342 = -Icges(6,4) * t452 - Icges(6,2) * t451 - Icges(6,6) * t737;
t344 = -Icges(6,1) * t452 - Icges(6,4) * t451 - Icges(6,5) * t737;
t115 = -t277 * t564 - t281 * t487 - t284 * t488 + t340 * t567 + t342 * t887 + t344 * t888;
t842 = -t115 / 0.2e1;
t341 = -Icges(6,5) * t454 - Icges(6,6) * t453 - Icges(6,3) * t734;
t343 = -Icges(6,4) * t454 - Icges(6,2) * t453 - Icges(6,6) * t734;
t345 = -Icges(6,1) * t454 - Icges(6,4) * t453 - Icges(6,5) * t734;
t116 = -t279 * t564 + t282 * t487 + t285 * t488 + t341 * t567 + t343 * t887 + t345 * t888;
t841 = t116 / 0.2e1;
t292 = Icges(6,5) * t416 - Icges(6,6) * t417;
t717 = -Icges(6,2) * t417 - t284 + t407;
t719 = Icges(6,1) * t416 + t281 - t766;
t121 = t292 * t567 + t717 * t887 + t719 * t888;
t840 = -t121 / 0.2e1;
t830 = t565 / 0.4e1;
t827 = -t568 / 0.4e1;
t788 = rSges(3,1) * t567;
t646 = pkin(1) + t788;
t683 = rSges(3,2) * t742 + t568 * rSges(3,3);
t440 = -t565 * t646 + t557 + t683;
t544 = rSges(3,2) * t740;
t441 = -t544 + t646 * t568 + (rSges(3,3) + pkin(6)) * t565;
t530 = rSges(3,1) * t564 + rSges(3,2) * t567;
t515 = t530 * t565;
t517 = t530 * t568;
t824 = m(3) * (t440 * t515 - t441 * t517);
t692 = t682 * t532;
t253 = t565 * (rSges(4,3) * t742 - t876) + t568 * (rSges(4,3) * t740 + t630) + t692;
t703 = -t422 * t737 - t424 * t734;
t820 = m(4) * (t253 * t518 + t703);
t818 = m(4) * (t318 * t401 + t319 * t402);
t311 = t319 * t740;
t817 = m(4) * (-t318 * t742 + t311);
t636 = -t565 * t642 + t568 * (t506 * pkin(3) + t496) + t692;
t204 = t565 * (rSges(5,2) * t742 - t689) + t568 * (t506 * rSges(5,1) + rSges(5,2) * t740 + t782) + t636;
t709 = -t348 * t737 - t350 * t734;
t812 = m(5) * (t204 * t518 + t709);
t809 = m(5) * (t268 * t322 + t269 * t321);
t264 = t269 * t740;
t808 = m(5) * (-t268 * t742 + t264);
t148 = (t506 * pkin(4) + t882) * t568 + (t501 + t948) * t565 + t636;
t721 = -t272 * t737 - t274 * t734;
t805 = m(6) * (t148 * t518 + t721);
t803 = m(6) * (t221 * t261 + t222 * t260);
t801 = m(6) * (t209 * t415 + t612 * t743);
t800 = m(6) * (t209 * t518 + t728);
t219 = t222 * t740;
t799 = m(6) * (-t221 * t742 + t219);
t798 = m(6) * (t242 * t505 - t503 * t942);
t797 = m(6) * (t242 * t740 - t742 * t942);
t220 = t565 * t298 + t299 * t568;
t796 = m(6) * (t220 * t743 - t406 * t415);
t795 = m(6) * (-t220 * t567 - t406 * t518);
t794 = m(6) * (-t298 * t505 + t299 * t503);
t633 = -t298 * t740 + t299 * t742;
t793 = m(6) * t633;
t789 = m(6) * qJD(5);
t150 = -t277 * t742 + t905;
t151 = -t279 * t742 + t726;
t616 = -t150 * t565 - t151 * t568;
t63 = t564 * t616 + t938;
t781 = t565 * t63;
t779 = t568 * t64;
t773 = t150 + t944;
t772 = t151 + t943;
t771 = t152 + t943;
t770 = t153 - t944;
t211 = t222 * t505;
t257 = t269 * t505;
t752 = t279 * t567;
t750 = t364 * t567;
t748 = t480 * t564;
t142 = -t221 * t503 + t211;
t190 = -t268 * t503 + t257;
t718 = Icges(6,1) * t420 - t282 - t765;
t716 = -Icges(6,2) * t421 + t285 + t408;
t558 = t564 ^ 2;
t698 = t503 * t737 + t505 * t734;
t302 = (-t567 ^ 2 + (0.1e1 - t682) * t558) * t561 + t698;
t715 = t302 * t932;
t320 = t518 * t743 + t698;
t712 = t320 * t932;
t326 = -t558 * t561 * t682 - t415 * t567;
t711 = t326 * t932;
t708 = t930 / 0.2e1;
t707 = -t930 / 0.2e1;
t705 = -t366 + t405;
t701 = t565 * t479 + t483 * t734;
t696 = -rSges(4,3) * t564 - t567 * t629 - t532;
t693 = t565 * (qJ(3) * t737 - t548) + t568 * (-pkin(2) * t740 + t537);
t687 = -t626 * t567 - t532;
t680 = qJD(2) * t564;
t679 = qJD(2) * t567;
t678 = qJD(5) * t564;
t18 = -t938 + (t565 * t770 - t568 * t771) * t564;
t673 = -t18 / 0.2e1 - t63 / 0.2e1;
t609 = -t341 * t564 - t752;
t100 = -t282 * t451 - t285 * t452 + t343 * t416 + t345 * t417 + t565 * t609;
t606 = -t365 * t564 - t750;
t133 = -t366 * t451 + t367 * t416 - t368 * t452 + t369 * t417 + t565 * t606;
t610 = -t340 * t564 - t754;
t99 = t281 * t451 + t284 * t452 + t342 * t416 + t344 * t417 + t565 * t610;
t14 = (t133 + t616) * t567 + (-t100 * t568 - t565 * t99 - t197) * t564;
t93 = -t292 * t742 + t416 * t717 + t417 * t719;
t293 = Icges(6,5) * t420 - Icges(6,6) * t421;
t94 = -t293 * t742 + t416 * t716 + t417 * t718;
t49 = t94 * t565 - t568 * t93;
t672 = -t49 / 0.2e1 + t14 / 0.2e1;
t101 = t281 * t453 + t284 * t454 + t420 * t342 + t421 * t344 + t568 * t610;
t102 = -t453 * t282 - t454 * t285 + t420 * t343 + t421 * t345 + t568 * t609;
t134 = -t453 * t366 + t420 * t367 - t454 * t368 + t421 * t369 + t568 * t606;
t15 = (t134 + t615) * t567 + (-t101 * t565 - t102 * t568 - t199) * t564;
t95 = -t292 * t740 + t420 * t717 + t421 * t719;
t96 = -t293 * t740 + t420 * t716 + t421 * t718;
t50 = t96 * t565 - t568 * t95;
t671 = -t50 / 0.2e1 + t15 / 0.2e1;
t19 = t196 + (t565 * t772 - t568 * t773) * t564;
t670 = t856 - t19 / 0.2e1;
t658 = -rSges(5,2) * t564 - t567 * t628 + t687;
t655 = -t742 / 0.4e1;
t654 = -t740 / 0.4e1;
t640 = t481 * t564 - t478;
t637 = pkin(7) * t564 - t567 * t791 - t371 + t687;
t634 = -t560 * t886 + t565 * t685 + t693;
t607 = t348 * t565 + t350 * t568;
t611 = t272 * t565 + t274 * t568;
t59 = m(5) * (t204 * t415 + t607 * t743) + m(6) * (t148 * t415 + t611 * t743);
t89 = m(5) * t190 + m(6) * t142;
t122 = t293 * t567 + t716 * t887 + t718 * t888;
t138 = -t403 * t742 + t416 * t704 + t417 * t705;
t139 = -t403 * t740 + t420 * t704 + t421 * t705;
t632 = t844 / 0.2e1 + (t122 + t139) * t830 + (t121 + t138) * t827;
t619 = -Icges(3,5) * t564 - Icges(3,6) * t567;
t174 = t282 * t887 + t285 * t888 + t752;
t614 = t173 * t565 - t174 * t568;
t26 = t565 * t771 + t568 * t770;
t27 = t565 * t773 + t568 * t772;
t79 = -t150 * t568 + t151 * t565;
t80 = -t152 * t568 + t153 * t565;
t582 = t18 * t830 + t19 * t827 + t27 * t655 + t781 / 0.4e1 + t779 / 0.4e1 + t80 * t742 / 0.4e1 + (t26 + t79) * t654;
t147 = -t564 * t364 + t567 * t365 + t744 + t745 + t746 + t747;
t214 = t366 * t887 + t368 * t888 + t750;
t572 = t147 * t829 + t214 * t903 + t851 / 0.2e1 + (t115 + t133) * t655 + (t116 + t134) * t654 - (-t173 + t197) * t737 / 0.4e1 - (t174 + t199) * t734 / 0.4e1;
t306 = -t481 * t740 + t701;
t571 = t80 / 0.2e1 - t27 / 0.2e1 + (-t442 + (t479 + t748) * t568 + t702 + t916) * t828 + (t701 + t917) * t832 + (t306 + t917) * t831;
t570 = t26 / 0.2e1 + t79 / 0.2e1 + (t568 * t640 + t306 - t701 + t924) * t568 / 0.2e1 + (-t565 * (-t482 * t567 + t748) - t478 * t568 + t924 + t937 * t742) * t828 + (t565 * t640 + t937 * t740 + t641 + t916) * t831;
t533 = -rSges(3,2) * t564 + t788;
t510 = t619 * t568;
t509 = t619 * t565;
t425 = t696 * t568;
t423 = t696 * t565;
t351 = t658 * t568;
t349 = t658 * t565;
t291 = 0.4e1 * t884;
t275 = t637 * t568;
t273 = t637 * t565;
t271 = t568 * (-t740 * t787 + t686) + t911 * t559 + t693;
t255 = 0.4e1 * t883;
t246 = t567 * t299 + t406 * t740;
t245 = -t298 * t567 - t406 * t742;
t239 = (-t568 * t885 + t546) * t568 + t910 * t559 + t634;
t213 = t793 / 0.2e1;
t206 = t794 / 0.2e1;
t187 = (-t519 * t568 + t347) * t568 + (-t519 * t565 + t346) * t565 + t634;
t183 = t795 / 0.2e1;
t182 = t707 + t931;
t181 = t707 + t708;
t180 = t708 - t931;
t178 = (t705 * t888 + t736 + t896) * t567;
t176 = t796 / 0.2e1;
t175 = t797 / 0.2e1;
t166 = t798 / 0.2e1;
t132 = t326 * t904 + t712 + t715;
t131 = t320 * t904 + t711 + t715;
t130 = t302 * t904 + t711 + t712;
t119 = t800 / 0.2e1;
t111 = t801 / 0.2e1;
t78 = t799 + t808 + t817;
t71 = t148 * t220 + t406 * t611;
t70 = (t368 / 0.2e1 + t404 / 0.2e1) * t887 + t906;
t67 = t214 * t567 + t564 * t614;
t53 = t805 + t812 + t820;
t52 = -t101 * t568 + t102 * t565;
t51 = t100 * t565 - t568 * t99;
t48 = t175 - t793 / 0.2e1;
t47 = t213 + t175;
t46 = t213 - t797 / 0.2e1;
t43 = t850 / 0.2e1;
t40 = t206 + t166;
t39 = t206 - t798 / 0.2e1;
t38 = t166 - t794 / 0.2e1;
t36 = t852 / 0.2e1;
t35 = t139 * t567 + (-t565 * t95 - t568 * t96) * t564;
t34 = t138 * t567 + (-t565 * t93 - t568 * t94) * t564;
t33 = t803 + t809 + t818 + t824 - t954;
t30 = t119 + t43 - t795 / 0.2e1;
t29 = t183 + t119 - t850 / 0.2e1;
t28 = t183 + t43 - t800 / 0.2e1;
t23 = t176 + t111 - t852 / 0.2e1;
t22 = t176 + t36 - t801 / 0.2e1;
t21 = t111 + t36 - t796 / 0.2e1;
t20 = (t147 + t614) * t567 + (-t115 * t565 - t116 * t568 - t214) * t564;
t13 = m(6) * t71 + t49 * t828 + t50 * t831;
t11 = t774 + t775;
t8 = t660 + t669;
t6 = (t565 * t670 + t568 * t673) * t564;
t5 = t849 + (-t779 / 0.2e1 - t781 / 0.2e1 + t20 / 0.2e1) * t567 + (t15 * t828 + t14 * t832 - t67 / 0.2e1) * t564;
t4 = t565 * t570 + t568 * t571;
t3 = ((t79 / 0.4e1 + t26 / 0.4e1) * t568 + (-t80 / 0.4e1 + t27 / 0.4e1) * t565) * t564 + (-t63 / 0.4e1 - t18 / 0.4e1) * t565 + (-t64 / 0.4e1 + t19 / 0.4e1) * t568 + t572 + t632;
t2 = -t851 / 0.2e1 + (t214 / 0.2e1 + (t134 / 0.4e1 + t116 / 0.4e1) * t568 + (t115 / 0.4e1 + t133 / 0.4e1) * t565) * t564 + (-t147 / 0.2e1 + (t199 / 0.4e1 + t174 / 0.4e1) * t568 + (t197 / 0.4e1 - t173 / 0.4e1) * t565) * t567 + t582 + t632;
t1 = -t844 / 0.2e1 + (t138 / 0.4e1 + t121 / 0.4e1) * t568 + t582 + (-t139 / 0.4e1 - t122 / 0.4e1) * t565 + t572;
t9 = [(-m(5) * t720 * t269 / 0.4e1 - t222 * t952 / 0.4e1) * t863 + t33 * qJD(2) + t78 * qJD(3) + t89 * qJD(4) + t70 * qJD(5), t33 * qJD(1) + t8 * qJD(3) + t11 * qJD(4) + t3 * qJD(5) + ((t318 * t425 + t319 * t423 - t401 * t424 - t402 * t422) * t860 + (t268 * t351 + t269 * t349 - t321 * t348 - t322 * t350) * t858 + (t221 * t275 + t222 * t273 - t260 * t272 - t261 * t274) * t857) * t862 + ((t432 / 0.2e1 + t430 / 0.2e1 - t482 / 0.2e1 - t511 / 0.2e1) * t679 + (t480 / 0.2e1 - t513 / 0.2e1 - t381 / 0.2e1 - t384 / 0.2e1) * t680) * t568 + ((-t431 / 0.2e1 + t483 / 0.2e1 - t512 / 0.2e1 - t433 / 0.2e1) * t679 + (t383 / 0.2e1 + t514 / 0.2e1 + t386 / 0.2e1 - t481 / 0.2e1) * t680) * t565 + ((m(3) * (-t440 * t533 - t515 * t530) + t842 - t133 / 0.2e1 + t895 / 0.2e1 + t647 * t504 + t649 * t503 - t571) * t568 + (t134 / 0.2e1 + t841 + m(3) * (-t441 * t533 + t517 * t530) + t523 * t831 - t647 * t506 - t649 * t505 - t570) * t565 + ((t877 * t564 + t567 * t927) * t831 + (t878 * t564 + t914 * t567) * t828) * t562 + ((t879 * t564 + t567 * t929) * t831 + (t880 * t564 + t915 * t567) * t828) * t561) * qJD(2), qJD(1) * t78 + qJD(2) * t8 + qJD(4) * t181 + qJD(5) * t47, t89 * qJD(1) + t11 * qJD(2) + t181 * qJD(3) + t40 * qJD(5), t70 * qJD(1) + t3 * qJD(2) + t47 * qJD(3) + t40 * qJD(4) + (m(6) * (t221 * t245 + t222 * t246 + t242 * t299 - t298 * t942) + t178) * qJD(5) + ((-t139 / 0.2e1 - t122 / 0.2e1 - t673) * t568 + (-t138 / 0.2e1 + t840 - t670) * t565) * t678; t4 * qJD(2) - t7 * qJD(3) - t12 * qJD(4) + t2 * qJD(5) + (t272 * t919 + t348 * t920) * t864 + (-t824 / 0.4e1 - t818 / 0.4e1 - t809 / 0.4e1 - t803 / 0.4e1) * t863 + t954 * qJD(1), t4 * qJD(1) + t53 * qJD(3) + t59 * qJD(4) + t13 * qJD(5) + (m(5) * (t204 * t239 - t348 * t349 - t350 * t351) + m(4) * (t253 * t271 - t422 * t423 - t424 * t425) + m(3) * ((t565 * (rSges(3,1) * t737 - t683) + t568 * (rSges(3,1) * t734 + t565 * rSges(3,3) - t544)) * (-t565 * t515 - t517 * t568) + t682 * t533 * t530) + m(6) * (t148 * t187 - t272 * t273 - t274 * t275) + (t559 * t510 + t52 + (-t505 * t880 - t506 * t878 + t865) * t568 + (t505 * t879 + t506 * t877 - t509 * t568) * t565) * t831 + (t560 * t509 + t51 + (-t503 * t880 - t504 * t878) * t568 + (t503 * t879 + t504 * t877 - t568 * t510 + t865) * t565) * t828) * qJD(2), -t956 + t53 * qJD(2) + t130 * qJD(4) + t29 * qJD(5) + (-0.4e1 * t884 + 0.2e1 * (t857 + t858 + t860) * (-t518 * t567 + t684)) * qJD(3), t59 * qJD(2) + t130 * qJD(3) - 0.4e1 * qJD(4) * t883 + t23 * qJD(5) - t955, t2 * qJD(1) + t13 * qJD(2) + t29 * qJD(3) + t23 * qJD(4) + (t67 / 0.2e1 + t671 * t568 + t672 * t565) * t678 + (m(6) * (t148 * t633 + t209 * t220 - t245 * t274 - t246 * t272 + t406 * t612) + t35 * t831 + t34 * t828 - t849 + (-t20 / 0.2e1 + (t840 + t856) * t568 + (t122 / 0.2e1 + t63 / 0.2e1) * t565) * t567) * qJD(5); t7 * qJD(2) + t182 * qJD(4) + t46 * qJD(5) + (-t799 / 0.4e1 - t808 / 0.4e1 - t817 / 0.4e1) * t863 + (t219 * t857 + t264 * t858 + t311 * t860 + ((-t222 * t857 - t269 * t858 - t319 * t860) * t568 + (-t919 - t920) * t565) * t564) * t864, t956 + t291 * qJD(3) + t131 * qJD(4) + t28 * qJD(5) + 0.4e1 * (-t805 / 0.4e1 - t812 / 0.4e1 - t820 / 0.4e1) * qJD(2) + ((-t567 * t187 + t721) * t857 + (-t567 * t239 + t709) * t858 + (-t567 * t271 + t703) * t860 + ((t273 * t565 + t275 * t568 + t148) * t857 + (t349 * t565 + t351 * t568 + t204) * t858 + (t423 * t565 + t425 * t568 + t253) * t860) * t564) * t862, t291 * qJD(2), qJD(1) * t182 + qJD(2) * t131, t46 * qJD(1) + t28 * qJD(2) + (-t633 * t567 + (t245 * t568 + t246 * t565) * t564) * t789; (m(6) * (t223 * t503 + t142 - t211) + m(5) * (t270 * t503 + t190 - t257) - t89) * qJD(1) + t12 * qJD(2) + t180 * qJD(3) + t39 * qJD(5), t955 + (m(6) * (t503 * t273 + t505 * t275) + m(5) * (t503 * t349 + t505 * t351) + 0.2e1 * ((t148 * t857 + t204 * t858) * t567 + ((t187 + t611) * t857 + (t239 + t607) * t858) * t564) * t561 - t59) * qJD(2) + t132 * qJD(3) + t255 * qJD(4) + t22 * qJD(5), qJD(1) * t180 + qJD(2) * t132, qJD(2) * t255, t39 * qJD(1) + t22 * qJD(2) + (t245 * t505 + t246 * t503 + t633 * t743) * t789; (-t242 * t952 - t896 / 0.2e1 - t906) * qJD(1) + t1 * qJD(2) + t48 * qJD(3) + t38 * qJD(4) + t6 * qJD(5), t1 * qJD(1) + (((-t80 / 0.2e1 + t842) * t567 - t672) * t568 + ((-t79 / 0.2e1 + t841) * t567 + t671) * t565 + ((-t52 / 0.2e1 - t173 / 0.2e1) * t568 + (-t51 / 0.2e1 - t174 / 0.2e1) * t565) * t564 + (t148 * t169 - t185 * t274 - t186 * t272 + t187 * t209 + t242 * t273 + t275 * t942 - t71) * m(6)) * qJD(2) + t30 * qJD(3) + t21 * qJD(4) + t5 * qJD(5), qJD(1) * t48 + qJD(2) * t30, qJD(1) * t38 + qJD(2) * t21, t6 * qJD(1) + t5 * qJD(2) + (m(6) * (t209 * t633 + t242 * t246 + t245 * t942) + t178 * t829 + (t35 * t828 + t34 * t832 + (-t121 * t565 - t122 * t568) * t829) * t564) * qJD(5);];
Cq = t9;
