% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:17:06
% DurationCPUTime: 46.31s
% Computational Cost: add. (37150->1265), mult. (43102->1641), div. (0->0), fcn. (39641->10), ass. (0->620)
t556 = qJ(2) + qJ(3);
t536 = sin(t556);
t561 = sin(qJ(1));
t836 = t536 * t561;
t498 = Icges(4,4) * t836;
t537 = cos(t556);
t832 = t537 * t561;
t563 = cos(qJ(1));
t871 = Icges(4,5) * t563;
t381 = Icges(4,1) * t832 - t498 - t871;
t558 = cos(pkin(9));
t827 = t558 * t563;
t557 = sin(pkin(9));
t829 = t557 * t561;
t432 = t537 * t829 + t827;
t820 = t563 * t557;
t828 = t558 * t561;
t433 = t537 * t828 - t820;
t238 = Icges(5,4) * t433 - Icges(5,2) * t432 + Icges(5,6) * t836;
t241 = Icges(5,1) * t433 - Icges(5,4) * t432 + Icges(5,5) * t836;
t953 = t238 * t557 - t241 * t558;
t991 = t381 - t953;
t665 = Icges(5,5) * t558 - Icges(5,6) * t557;
t354 = -Icges(5,3) * t537 + t536 * t665;
t668 = Icges(5,4) * t558 - Icges(5,2) * t557;
t356 = -Icges(5,6) * t537 + t536 * t668;
t673 = Icges(5,1) * t558 - Icges(5,4) * t557;
t358 = -Icges(5,5) * t537 + t536 * t673;
t876 = Icges(4,4) * t536;
t449 = Icges(4,2) * t537 + t876;
t520 = Icges(4,4) * t537;
t451 = Icges(4,1) * t536 + t520;
t649 = t449 * t536 - t451 * t537;
t447 = Icges(4,5) * t536 + Icges(4,6) * t537;
t841 = t447 * t563;
t990 = -t354 * t836 + t356 * t432 - t358 * t433 + t561 * t649 + t841;
t434 = -t537 * t820 + t828;
t435 = t537 * t827 + t829;
t835 = t536 * t563;
t842 = t447 * t561;
t989 = t354 * t835 + t356 * t434 + t358 * t435 - t563 * t649 + t842;
t552 = pkin(9) + qJ(5);
t532 = sin(t552);
t533 = cos(t552);
t396 = t532 * t832 + t533 * t563;
t821 = t563 * t532;
t397 = t533 * t832 - t821;
t201 = Icges(6,5) * t397 - Icges(6,6) * t396 + Icges(6,3) * t836;
t375 = Icges(6,4) * t397;
t204 = -Icges(6,2) * t396 + Icges(6,6) * t836 + t375;
t374 = Icges(6,4) * t396;
t208 = -Icges(6,1) * t397 - Icges(6,5) * t836 + t374;
t984 = t204 * t532 + t208 * t533;
t93 = -t201 * t537 - t984 * t536;
t669 = -Icges(4,2) * t536 + t520;
t380 = Icges(4,6) * t561 + t563 * t669;
t452 = Icges(4,1) * t537 - t876;
t382 = Icges(4,5) * t561 + t452 * t563;
t323 = t382 * t832;
t448 = Icges(4,5) * t537 - Icges(4,6) * t536;
t378 = Icges(4,3) * t561 + t448 * t563;
t708 = t378 * t563 - t323;
t149 = -t380 * t836 - t708;
t237 = Icges(5,5) * t435 + Icges(5,6) * t434 + Icges(5,3) * t835;
t240 = Icges(5,4) * t435 + Icges(5,2) * t434 + Icges(5,6) * t835;
t243 = Icges(5,1) * t435 + Icges(5,4) * t434 + Icges(5,5) * t835;
t98 = t237 * t836 - t432 * t240 + t433 * t243;
t972 = t149 + t98;
t100 = t237 * t835 + t434 * t240 + t435 * t243;
t831 = t537 * t563;
t809 = t561 * t378 + t382 * t831;
t151 = -t380 * t835 + t809;
t966 = t100 + t151;
t355 = Icges(5,3) * t536 + t537 * t665;
t535 = qJD(2) * t561;
t475 = qJD(3) * t561 + t535;
t553 = qJD(2) + qJD(3);
t476 = t553 * t563;
t654 = -t356 * t557 + t358 * t558;
t657 = -t240 * t557 + t243 * t558;
t957 = t451 + t669;
t959 = -t449 * t563 + t382;
t988 = (-Icges(4,2) * t832 + t354 * t561 - t498 + t991) * t476 + (-t354 * t563 - t657 - t959) * t475 + (t355 - t654 - t957) * qJD(1);
t235 = Icges(5,5) * t433 - Icges(5,6) * t432 + Icges(5,3) * t836;
t99 = t235 * t835 + t434 * t238 + t435 * t241;
t986 = t989 * qJD(1) - t476 * t99;
t97 = t235 * t836 - t238 * t432 + t241 * t433;
t985 = t990 * qJD(1) + t476 * t97;
t664 = Icges(6,5) * t533 - Icges(6,6) * t532;
t328 = -Icges(6,3) * t537 + t536 * t664;
t873 = Icges(6,4) * t533;
t667 = -Icges(6,2) * t532 + t873;
t330 = -Icges(6,6) * t537 + t536 * t667;
t874 = Icges(6,4) * t532;
t672 = Icges(6,1) * t533 - t874;
t332 = -Icges(6,5) * t537 + t536 * t672;
t120 = t328 * t836 - t330 * t396 + t332 * t397;
t765 = qJD(5) * t563;
t415 = t536 * t765 + t475;
t766 = qJD(5) * t561;
t416 = -t536 * t766 + t476;
t767 = qJD(5) * t537;
t505 = qJD(1) - t767;
t85 = t201 * t836 - t204 * t396 - t208 * t397;
t398 = t533 * t561 - t537 * t821;
t399 = t532 * t561 + t533 * t831;
t203 = Icges(6,5) * t399 + Icges(6,6) * t398 + Icges(6,3) * t835;
t875 = Icges(6,4) * t399;
t206 = Icges(6,2) * t398 + Icges(6,6) * t835 + t875;
t376 = Icges(6,4) * t398;
t209 = Icges(6,1) * t399 + Icges(6,5) * t835 + t376;
t86 = t203 * t836 - t396 * t206 + t397 * t209;
t37 = t120 * t505 + t415 * t86 - t416 * t85;
t121 = t328 * t835 + t330 * t398 + t332 * t399;
t87 = t201 * t835 + t398 * t204 - t208 * t399;
t88 = t203 * t835 + t398 * t206 + t399 * t209;
t38 = t121 * t505 + t415 * t88 - t416 * t87;
t830 = t553 * t561;
t757 = t536 * t830;
t309 = qJD(1) * t434 + t557 * t757;
t310 = qJD(1) * t435 - t558 * t757;
t772 = qJD(1) * t563;
t740 = t536 * t772;
t612 = t537 * t830 + t740;
t139 = Icges(5,5) * t310 + Icges(5,6) * t309 + Icges(5,3) * t612;
t141 = Icges(5,4) * t310 + Icges(5,2) * t309 + Icges(5,6) * t612;
t143 = Icges(5,1) * t310 + Icges(5,4) * t309 + Icges(5,5) * t612;
t756 = t536 * t476;
t307 = qJD(1) * t432 + t557 * t756;
t308 = -qJD(1) * t433 - t558 * t756;
t861 = Icges(4,3) * t563;
t377 = Icges(4,5) * t832 - Icges(4,6) * t836 - t861;
t865 = Icges(4,6) * t563;
t379 = Icges(4,4) * t832 - Icges(4,2) * t836 - t865;
t961 = -t451 * t561 - t379;
t709 = qJD(1) * t382 + t553 * t961;
t711 = qJD(1) * t380 + t381 * t553 - t449 * t830;
t580 = qJD(1) * t377 - t536 * t711 + t537 * t709;
t653 = t379 * t536 - t381 * t537;
t777 = qJD(1) * t378;
t935 = qJD(1) * t653 - t553 * t842 + t777;
t773 = qJD(1) * t561;
t741 = t536 * t773;
t755 = t476 * t537;
t943 = -t741 + t755;
t983 = -t139 * t835 - t141 * t434 - t143 * t435 - t235 * t943 - t238 * t307 - t241 * t308 - t561 * t935 - t580 * t563;
t138 = Icges(5,5) * t308 + Icges(5,6) * t307 + Icges(5,3) * t943;
t140 = Icges(5,4) * t308 + Icges(5,2) * t307 + Icges(5,6) * t943;
t142 = Icges(5,1) * t308 + Icges(5,4) * t307 + Icges(5,5) * t943;
t960 = -t451 * t563 - t380;
t710 = (-t452 * t561 + t871) * qJD(1) + t960 * t553;
t712 = (-t561 * t669 + t865) * qJD(1) + t959 * t553;
t579 = -t536 * t712 + t537 * t710 + t777;
t847 = t380 * t536;
t936 = -t553 * t841 + (-t382 * t537 - t448 * t561 + t847 + t861) * qJD(1);
t982 = t138 * t835 + t140 * t434 + t142 * t435 + t237 * t943 + t240 * t307 + t243 * t308 + t561 * t936 + t579 * t563;
t981 = -t139 * t836 + t141 * t432 - t143 * t433 - t235 * t612 - t238 * t309 - t241 * t310 - t580 * t561 + t563 * t935;
t980 = t138 * t836 - t140 * t432 + t142 * t433 + t237 * t612 + t240 * t309 + t243 * t310 + t579 * t561 - t563 * t936;
t621 = t653 * t561;
t848 = t377 * t563;
t148 = -t621 - t848;
t979 = -t148 * t476 + t972 * t475 - t985;
t810 = -t561 * t377 - t381 * t831;
t150 = -t379 * t835 - t810;
t978 = -t150 * t476 + t966 * t475 + t986;
t299 = t355 * t553;
t357 = Icges(5,6) * t536 + t537 * t668;
t300 = t357 * t553;
t359 = Icges(5,5) * t536 + t537 * t673;
t301 = t359 * t553;
t958 = -t449 + t452;
t704 = t958 * t553;
t705 = t957 * t553;
t578 = qJD(1) * t447 - t536 * t705 + t537 * t704;
t932 = qJD(1) * t649 + t448 * t553;
t977 = t299 * t835 + t300 * t434 + t301 * t435 + t307 * t356 + t308 * t358 + t354 * t943 + t561 * t932 + t578 * t563;
t976 = t299 * t836 - t300 * t432 + t301 * t433 + t309 * t356 + t310 * t358 + t354 * t612 + t578 * t561 - t563 * t932;
t975 = (t553 * t953 + t139 - t711) * t537 + (t141 * t557 - t143 * t558 - t235 * t553 - t709) * t536;
t974 = (t553 * t657 - t138 + t712) * t537 + (-t140 * t557 + t142 * t558 + t237 * t553 + t710) * t536;
t973 = t148 + t97;
t971 = t150 + t99;
t965 = (-t235 + t379) * t537 + t991 * t536;
t964 = (-t237 + t380) * t537 + (t382 + t657) * t536;
t685 = rSges(6,1) * t397 - rSges(6,2) * t396;
t210 = rSges(6,3) * t836 + t685;
t522 = pkin(4) * t558 + pkin(3);
t444 = t522 * t832;
t760 = pkin(4) * t820;
t905 = pkin(3) * t537;
t559 = -pkin(8) - qJ(4);
t819 = qJ(4) + t559;
t948 = t536 * t819;
t232 = t760 - t444 + (t905 + t948) * t561;
t963 = t210 - t232;
t713 = t819 * t537;
t900 = pkin(3) - t522;
t726 = t900 * t536;
t320 = t713 - t726;
t684 = rSges(6,1) * t533 - rSges(6,2) * t532;
t334 = -rSges(6,3) * t537 + t536 * t684;
t962 = t320 + t334;
t583 = qJD(1) * t958 + t475 * t960 - t476 * t961;
t956 = (qJD(1) * t354 - t235 * t476 + t237 * t475 + t583) * t537 + t988 * t536;
t955 = t37 * t561 + t38 * t563;
t954 = -t210 * t505 - t416 * t334;
t952 = 0.2e1 * qJD(2);
t560 = sin(qJ(2));
t951 = rSges(3,2) * t560;
t645 = t563 * t505;
t699 = qJD(1) * t537 - qJD(5);
t938 = t561 * t699 + t756;
t181 = t532 * t938 + t533 * t645;
t182 = t532 * t645 - t533 * t938;
t754 = t182 * rSges(6,1) + t181 * rSges(6,2) + rSges(6,3) * t755;
t117 = -rSges(6,3) * t741 + t754;
t464 = qJ(4) * t755;
t787 = qJD(1) * t760 + t559 * t741;
t833 = t537 * t559;
t858 = qJ(4) * t536;
t947 = t537 * t900;
t134 = -t464 + (t726 - t833) * t476 + (t858 + t947) * t773 + t787;
t628 = t684 * t537;
t683 = -rSges(6,1) * t532 - rSges(6,2) * t533;
t177 = t553 * t628 + (rSges(6,3) * t553 + qJD(5) * t683) * t536;
t212 = t399 * rSges(6,1) + t398 * rSges(6,2) + rSges(6,3) * t835;
t762 = qJD(1) * qJD(2);
t526 = t563 * t762;
t761 = qJD(1) * qJD(3);
t458 = t563 * t761 + t526;
t303 = qJD(5) * t943 + t458;
t769 = qJD(4) * t563;
t491 = t536 * t769;
t613 = -t537 * t773 - t756;
t213 = pkin(3) * t613 - qJ(4) * t741 + t464 + t491;
t770 = qJD(4) * t561;
t490 = t537 * t770;
t564 = -pkin(7) - pkin(6);
t530 = t563 * t564;
t531 = pkin(6) * t772;
t771 = qJD(2) * t563;
t735 = t560 * t771;
t700 = pkin(2) * t735;
t562 = cos(qJ(2));
t527 = pkin(2) * t562 + pkin(1);
t901 = pkin(1) - t527;
t284 = -t700 - t531 + (t561 * t901 - t530) * qJD(1);
t446 = qJD(1) * (-pkin(1) * t773 + t531);
t822 = t562 * qJD(2) ^ 2;
t587 = qJD(1) * t284 + t446 + (-t526 * t560 - t561 * t822) * pkin(2);
t574 = t553 * t490 + t587 + (t213 + t491) * qJD(1);
t768 = qJD(5) * t536;
t732 = t553 * t768;
t453 = pkin(3) * t536 - qJ(4) * t537;
t811 = -t320 - t453;
t321 = -t947 - t948;
t259 = t321 * t553;
t838 = t536 * t553;
t348 = qJ(4) * t838 + (pkin(3) * t553 - qJD(4)) * t537;
t814 = -t259 - t348;
t29 = qJD(1) * t134 + t117 * t505 - t177 * t415 + t212 * t732 - t303 * t334 + t458 * t811 + t475 * t814 + t574;
t950 = t29 * t561;
t949 = -rSges(5,3) - qJ(4);
t550 = t563 * pkin(6);
t493 = pkin(1) * t561 - t550;
t782 = -t561 * t527 - t530;
t372 = t493 + t782;
t472 = qJD(1) * t493;
t944 = qJD(1) * t372 - t472;
t544 = Icges(3,4) * t562;
t670 = -Icges(3,2) * t560 + t544;
t482 = Icges(3,1) * t560 + t544;
t549 = t561 * pkin(6);
t494 = t563 * pkin(1) + t549;
t294 = t334 * t561;
t894 = rSges(6,3) * t536;
t335 = t628 + t894;
t431 = t453 * t773;
t731 = t537 * t766;
t455 = t858 + t905;
t492 = t537 * t769;
t783 = -t476 * t455 + t492;
t942 = t210 * t768 - t294 * t505 + t476 * t321 - t334 * t731 + t416 * t335 + t773 * t962 + t431 - t783;
t390 = rSges(4,1) * t832 - rSges(4,2) * t836 - t563 * rSges(4,3);
t545 = t561 * rSges(4,3);
t391 = rSges(4,1) * t831 - rSges(4,2) * t835 + t545;
t500 = t563 * t527;
t702 = -t561 * t564 + t500;
t373 = t702 - t494;
t804 = -t372 * t535 + t373 * t771;
t128 = t390 * t475 + t391 * t476 + t804;
t454 = rSges(4,1) * t536 + rSges(4,2) * t537;
t614 = -t454 * t476 - t700;
t798 = t372 - t493;
t146 = (-t390 + t798) * qJD(1) + t614;
t892 = pkin(2) * qJD(2);
t759 = t560 * t892;
t509 = t561 * t759;
t797 = t373 + t494;
t147 = -t454 * t475 - t509 + (t391 + t797) * qJD(1);
t426 = t454 * t561;
t429 = t454 * t563;
t898 = rSges(4,1) * t537;
t456 = -rSges(4,2) * t536 + t898;
t939 = -t146 * (qJD(1) * t426 - t476 * t456) - t128 * (-t475 * t426 - t429 * t476) - t147 * (-qJD(1) * t429 - t456 * t475);
t687 = rSges(5,1) * t558 - rSges(5,2) * t557;
t361 = -rSges(5,3) * t537 + t536 * t687;
t824 = t561 * t562;
t826 = t560 * t561;
t862 = Icges(3,3) * t563;
t407 = Icges(3,5) * t824 - Icges(3,6) * t826 - t862;
t516 = Icges(3,4) * t826;
t872 = Icges(3,5) * t563;
t411 = Icges(3,1) * t824 - t516 - t872;
t866 = Icges(3,6) * t563;
t409 = Icges(3,4) * t824 - Icges(3,2) * t826 - t866;
t845 = t409 * t560;
t651 = -t411 * t562 + t845;
t160 = -t407 * t563 - t561 * t651;
t937 = t563 * t699 - t757;
t479 = Icges(3,5) * t562 - Icges(3,6) * t560;
t478 = Icges(3,5) * t560 + Icges(3,6) * t562;
t617 = qJD(2) * t478;
t877 = Icges(3,4) * t560;
t483 = Icges(3,1) * t562 - t877;
t412 = Icges(3,5) * t561 + t483 * t563;
t410 = Icges(3,6) * t561 + t563 * t670;
t844 = t410 * t560;
t650 = -t412 * t562 + t844;
t934 = -t563 * t617 + (-t479 * t561 + t650 + t862) * qJD(1);
t408 = Icges(3,3) * t561 + t479 * t563;
t776 = qJD(1) * t408;
t933 = qJD(1) * t651 - t561 * t617 + t776;
t480 = Icges(3,2) * t562 + t877;
t648 = t480 * t560 - t482 * t562;
t931 = t648 * qJD(1) + t479 * qJD(2);
t930 = t561 * (-t480 * t563 + t412) - t563 * (-Icges(3,2) * t824 + t411 - t516);
t644 = t505 * t561;
t183 = -t532 * t937 + t533 * t644;
t184 = t532 * t644 + t533 * t937;
t686 = rSges(6,1) * t184 + rSges(6,2) * t183;
t118 = rSges(6,3) * t612 + t686;
t469 = pkin(3) * t757;
t512 = pkin(4) * t829;
t135 = t469 + (-t522 * t536 - t713) * t830 + (t321 * t563 + t512) * qJD(1);
t525 = t561 * t762;
t457 = t561 * t761 + t525;
t302 = qJD(5) * t612 + t457;
t647 = -pkin(2) * t563 * t822 + qJD(1) * t509;
t603 = -t476 * t348 + t457 * t453 + t553 * t492 + t647;
t506 = pkin(3) * t831;
t733 = t536 * t770;
t214 = qJ(4) * t612 + qJD(1) * t506 - t469 + t733;
t780 = t564 * t773 + t509;
t285 = (-t563 * t901 - t549) * qJD(1) - t780;
t471 = t494 * qJD(1);
t813 = -t285 - t471;
t620 = -t214 - t733 + t813;
t30 = -t210 * t732 - t118 * t505 - t177 * t416 - t259 * t476 + t302 * t334 + t320 * t457 + (-t135 + t620) * qJD(1) + t603;
t430 = qJ(4) * t835 + t506;
t789 = t522 * t831 + t512;
t233 = -t559 * t835 - t430 + t789;
t641 = -t475 * t453 - t509 + t733;
t743 = t430 + t797;
t75 = t212 * t505 - t320 * t475 - t334 * t415 + (t233 + t743) * qJD(1) + t641;
t929 = (qJD(1) * t75 + t30) * t563 + t950;
t622 = t664 * t537;
t655 = -t330 * t532 + t332 * t533;
t660 = -t206 * t532 + t209 * t533;
t928 = t415 * (-t328 * t563 - t660) - t416 * (-t328 * t561 + t984) + t505 * (Icges(6,3) * t536 + t622 - t655);
t666 = -Icges(6,2) * t533 - t874;
t927 = t415 * (-Icges(6,2) * t399 + t209 + t376) - t416 * (-Icges(6,2) * t397 - t208 - t374) + t505 * (t666 * t536 + t332);
t924 = m(5) / 0.2e1;
t923 = m(6) / 0.2e1;
t922 = t302 / 0.2e1;
t921 = t303 / 0.2e1;
t920 = -t415 / 0.2e1;
t919 = t415 / 0.2e1;
t918 = -t416 / 0.2e1;
t917 = t416 / 0.2e1;
t916 = t457 / 0.2e1;
t915 = t458 / 0.2e1;
t914 = -t475 / 0.2e1;
t913 = t475 / 0.2e1;
t912 = -t476 / 0.2e1;
t911 = t476 / 0.2e1;
t910 = -t505 / 0.2e1;
t909 = t505 / 0.2e1;
t908 = t561 / 0.2e1;
t907 = -t563 / 0.2e1;
t906 = pkin(2) * t560;
t904 = pkin(4) * t557;
t903 = -qJD(1) / 0.2e1;
t902 = qJD(1) / 0.2e1;
t899 = rSges(3,1) * t562;
t897 = rSges(3,2) * t562;
t896 = rSges(5,3) * t536;
t112 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t612;
t114 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t612;
t116 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t612;
t27 = (-t553 * t984 - t112) * t537 + (-t114 * t532 + t116 * t533 + t201 * t553 + (-t204 * t533 + t208 * t532) * qJD(5)) * t536;
t891 = t27 * t416;
t111 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t943;
t113 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t943;
t115 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t943;
t28 = (t553 * t660 - t111) * t537 + (-t113 * t532 + t115 * t533 + t203 * t553 + (-t206 * t533 - t209 * t532) * qJD(5)) * t536;
t890 = t28 * t415;
t546 = t561 * rSges(3,3);
t883 = t93 * t302;
t94 = -t203 * t537 + t536 * t660;
t882 = t94 * t303;
t881 = -rSges(6,3) + t559;
t133 = -t328 * t537 + t536 * t655;
t663 = -Icges(6,5) * t532 - Icges(6,6) * t533;
t173 = t553 * t622 + (Icges(6,3) * t553 + qJD(5) * t663) * t536;
t623 = t667 * t537;
t174 = t553 * t623 + (Icges(6,6) * t553 + qJD(5) * t666) * t536;
t624 = t672 * t537;
t671 = -Icges(6,1) * t532 - t873;
t175 = t553 * t624 + (Icges(6,5) * t553 + qJD(5) * t671) * t536;
t57 = (t553 * t655 - t173) * t537 + (-t174 * t532 + t175 * t533 + t328 * t553 + (-t330 * t533 - t332 * t532) * qJD(5)) * t536;
t880 = t133 * t732 + t57 * t505;
t427 = t455 * t561;
t640 = -qJD(4) * t537 + t475 * t427 + t804;
t816 = t233 + t430;
t65 = t210 * t415 + t212 * t416 - t232 * t475 + t476 * t816 + t640;
t857 = qJD(1) * t65;
t688 = rSges(5,1) * t433 - rSges(5,2) * t432;
t244 = rSges(5,3) * t836 + t688;
t246 = t435 * rSges(5,1) + t434 * rSges(5,2) + rSges(5,3) * t835;
t815 = t246 + t430;
t84 = t244 * t475 + t476 * t815 + t640;
t855 = qJD(1) * t84;
t102 = -t361 * t475 + (t246 + t743) * qJD(1) + t641;
t854 = t102 * t561;
t852 = t146 * t561;
t779 = rSges(3,2) * t826 + t563 * rSges(3,3);
t417 = rSges(3,1) * t824 - t779;
t485 = rSges(3,1) * t560 + t897;
t736 = t485 * t771;
t257 = -t736 + (-t417 - t493) * qJD(1);
t851 = t257 * t561;
t850 = t257 * t563;
t737 = t485 * t535;
t823 = t562 * t563;
t825 = t560 * t563;
t418 = rSges(3,1) * t823 - rSges(3,2) * t825 + t546;
t790 = t418 + t494;
t258 = qJD(1) * t790 - t737;
t443 = t485 * t563;
t849 = t258 * t443;
t840 = t478 * t561;
t839 = t478 * t563;
t837 = t536 * t559;
t834 = t537 * t553;
t362 = t537 * t687 + t896;
t305 = t362 * t553;
t812 = -t305 - t348;
t808 = t361 * t773 + t431;
t425 = t453 * t561;
t519 = qJD(4) * t536;
t806 = -t475 * t425 + t519;
t803 = -t561 * t372 + t563 * t373;
t802 = t561 * t390 + t563 * t391;
t801 = -t361 - t453;
t800 = -t561 * t407 - t411 * t823;
t799 = t561 * t408 + t412 * t823;
t794 = t561 * t427 + t563 * t430;
t428 = t453 * t563;
t793 = -qJD(1) * t428 + t490;
t786 = rSges(4,2) * t741 + rSges(4,3) * t772;
t785 = -t480 + t483;
t784 = t482 + t670;
t781 = rSges(3,3) * t772 + t773 * t951;
t775 = qJD(1) * t425;
t774 = qJD(1) * t479;
t263 = -t561 * t648 - t839;
t763 = t263 * qJD(1);
t758 = t562 * t892;
t753 = -t177 + t814;
t752 = t563 * t213 + t561 * t214 + t427 * t772;
t751 = -t212 - t816;
t230 = rSges(4,1) * t613 - rSges(4,2) * t755 + t786;
t231 = -t553 * t426 + (t456 * t563 + t545) * qJD(1);
t750 = t563 * t230 + t561 * t231 + t390 * t772;
t749 = t563 * t284 + t561 * t285 - t372 * t772;
t747 = t308 * rSges(5,1) + t307 * rSges(5,2) + rSges(5,3) * t755;
t744 = -t427 + t798;
t742 = -t476 * t362 + t783;
t738 = t560 * t772;
t730 = t537 * t765;
t729 = t838 / 0.2e1;
t727 = -pkin(1) - t899;
t725 = t773 / 0.2e1;
t724 = t772 / 0.2e1;
t723 = -t535 / 0.2e1;
t720 = t771 / 0.2e1;
t718 = -t453 - t906;
t717 = -t454 - t906;
t715 = t102 * t801;
t714 = t560 * (-t561 ^ 2 - t563 ^ 2);
t363 = t412 * t824;
t707 = t408 * t563 - t363;
t706 = -t377 + t847;
t703 = -t407 + t844;
t698 = t561 * t244 + t563 * t246 + t794;
t697 = qJD(5) * t729;
t696 = -t361 + t718;
t693 = -t348 - t758;
t406 = t456 * t553;
t692 = -t406 - t758;
t690 = t899 - t951;
t689 = rSges(5,1) * t310 + rSges(5,2) * t309;
t643 = t491 - t700;
t627 = -t476 * t453 + t643;
t74 = -t320 * t476 + (t232 + t744) * qJD(1) + t627 + t954;
t682 = t561 * t75 + t563 * t74;
t681 = t561 * t86 - t563 * t85;
t680 = t561 * t85 + t563 * t86;
t679 = t561 * t88 - t563 * t87;
t678 = t561 * t87 + t563 * t88;
t677 = t561 * t94 - t563 * t93;
t676 = t561 * t93 + t563 * t94;
t662 = -t146 * t563 - t147 * t561;
t659 = t210 * t563 - t212 * t561;
t656 = -t258 * t561 - t850;
t261 = t409 * t562 + t411 * t560;
t262 = t410 * t562 + t412 * t560;
t642 = -t305 + t693;
t144 = -rSges(5,3) * t741 + t747;
t145 = rSges(5,3) * t612 + t689;
t639 = t563 * t144 + t561 * t145 + t244 * t772 + t752;
t638 = t794 + (t212 + t233) * t563 + t963 * t561;
t637 = -t522 * t537 - t527 - t894;
t442 = t485 * t561;
t626 = -t177 - t259 + t693;
t619 = qJD(2) * t482;
t618 = qJD(2) * t480;
t161 = -t410 * t826 - t707;
t616 = (-t160 * t563 + t161 * t561) * qJD(2);
t162 = -t409 * t825 - t800;
t163 = -t410 * t825 + t799;
t615 = (-t162 * t563 + t163 * t561) * qJD(2);
t247 = (t417 * t561 + t418 * t563) * qJD(2);
t610 = -t455 - t896;
t609 = -t201 * t416 + t203 * t415 + t328 * t505;
t608 = (-Icges(6,5) * t396 - Icges(6,6) * t397) * t416 - (Icges(6,5) * t398 - Icges(6,6) * t399) * t415 - t663 * t536 * t505;
t607 = qJD(1) * t448 - t475 * t841 + t476 * t842;
t606 = t284 * t771 + t285 * t535 - t372 * t526 - t373 * t525;
t604 = t409 * t563 - t410 * t561;
t602 = -qJD(1) * t427 + t643 + t944;
t601 = t752 + t963 * t772 + (t117 + t134) * t563 + (t118 + t135) * t561;
t600 = t536 * t608;
t591 = (-t560 * t784 + t562 * t785) * qJD(1);
t590 = t475 * t214 + t458 * t427 + t553 * t519 + t606;
t586 = (Icges(6,1) * t398 - t206 - t875) * t415 - (-Icges(6,1) * t396 - t204 - t375) * t416 + (t671 * t536 - t330) * t505;
t277 = qJD(1) * t410 - t561 * t618;
t279 = qJD(1) * t412 - t561 * t619;
t577 = qJD(1) * t407 - qJD(2) * t261 - t277 * t560 + t279 * t562;
t276 = -t563 * t618 + (-t561 * t670 + t866) * qJD(1);
t278 = -t563 * t619 + (-t483 * t561 + t872) * qJD(1);
t576 = -qJD(2) * t262 - t276 * t560 + t278 * t562 + t776;
t462 = t670 * qJD(2);
t463 = t483 * qJD(2);
t575 = qJD(1) * t478 - t462 * t560 + t463 * t562 + (-t480 * t562 - t482 * t560) * qJD(2);
t318 = t361 * t561;
t319 = t361 * t563;
t573 = t102 * (-qJD(1) * t319 + (-t362 - t455) * t475 + t793) + t84 * (-t475 * t318 + (-t319 - t428) * t476 + t806);
t572 = -t560 * t930 + t604 * t562;
t571 = t928 * t536;
t21 = t112 * t835 + t114 * t398 + t116 * t399 + t181 * t204 - t182 * t208 + t201 * t943;
t22 = t111 * t835 + t113 * t398 + t115 * t399 + t181 * t206 + t182 * t209 + t203 * t943;
t23 = t112 * t836 - t114 * t396 + t116 * t397 + t183 * t204 - t184 * t208 + t201 * t612;
t24 = t111 * t836 - t113 * t396 + t115 * t397 + t183 * t206 + t184 * t209 + t203 * t612;
t288 = t330 * t561;
t289 = t330 * t563;
t290 = t332 * t561;
t291 = t332 * t563;
t49 = t173 * t835 + t174 * t398 + t175 * t399 + t181 * t330 + t182 * t332 + t328 * t943;
t3 = t121 * t732 - t21 * t416 + t22 * t415 + t302 * t87 + t303 * t88 + t49 * t505;
t314 = t356 * t561;
t315 = t356 * t563;
t316 = t358 * t561;
t317 = t358 * t563;
t331 = Icges(6,6) * t536 + t623;
t333 = Icges(6,5) * t536 + t624;
t50 = t173 * t836 - t174 * t396 + t175 * t397 + t183 * t330 + t184 * t332 + t328 * t612;
t4 = t120 * t732 - t23 * t416 + t24 * t415 + t302 * t85 + t303 * t86 + t50 * t505;
t44 = t133 * t505 + t415 * t94 - t416 * t93;
t569 = t677 * t697 + (((t289 * t532 - t291 * t533 + t203) * t415 - (t288 * t532 - t290 * t533 + t201) * t416 + (-t331 * t532 + t333 * t533 + t328) * t505 + t133 * qJD(5)) * t536 + (qJD(5) * t676 - t928) * t537) * t910 - t44 * t768 / 0.2e1 + (qJD(1) * t678 - t21 * t563 + t22 * t561) * t919 + ((-t289 * t398 - t291 * t399) * t415 - (-t288 * t398 - t290 * t399) * t416 + (t331 * t398 + t333 * t399) * t505 + (t121 * t536 + t832 * t87) * qJD(5) + ((qJD(5) * t88 + t609) * t537 + t571) * t563) * t920 + t679 * t921 + t681 * t922 + ((t289 * t396 - t291 * t397) * t415 - (t288 * t396 - t290 * t397) * t416 + (-t331 * t396 + t333 * t397) * t505 + (t120 * t536 + t831 * t86) * qJD(5) + ((qJD(5) * t85 + t609) * t537 + t571) * t561) * t917 + (qJD(1) * t680 - t23 * t563 + t24 * t561) * t918 + (qJD(1) * t676 - t27 * t563 + t28 * t561) * t909 + (t561 * t972 - t563 * t973) * t916 + (t561 * t966 - t563 * t971) * t915 + ((-t315 * t434 - t317 * t435) * t475 - (-t314 * t434 - t316 * t435) * t476 + (t357 * t434 + t359 * t435) * qJD(1) + t561 * t607 + t956 * t563) * t914 + (t983 * t563 + t982 * t561 + (t561 * t971 + t563 * t966) * qJD(1)) * t913 + (t981 * t563 + t980 * t561 + (t561 * t973 + t563 * t972) * qJD(1)) * t912 + ((t315 * t432 - t317 * t433) * t475 - (t314 * t432 - t316 * t433) * t476 + (-t357 * t432 + t359 * t433) * qJD(1) - t563 * t607 + t956 * t561) * t911 + (-t988 * t537 + ((t315 * t557 - t317 * t558 + t237) * t475 - (t314 * t557 - t316 * t558 + t235) * t476 + (-t357 * t557 + t359 * t558 + t354) * qJD(1) + t583) * t536) * t903 + (t975 * t563 + t974 * t561 + (t561 * t965 + t563 * t964) * qJD(1)) * t902 - t955 * t767 / 0.2e1 + (t977 * qJD(1) + t971 * t457 + t966 * t458 + t982 * t475 + t476 * t983 + t3) * t908 + (qJD(1) * t976 + t457 * t973 + t458 * t972 + t475 * t980 + t476 * t981 + t4) * t907 + (t37 + t979) * t725 + (t38 + t978) * t724;
t268 = t320 * t561;
t269 = t320 * t563;
t295 = t334 * t563;
t568 = t65 * (-t212 * t731 - t295 * t416 + t210 * t730 - t415 * t294 - t475 * t268 + (-t269 - t428) * t476 + t806) + t75 * (-t334 * t730 - t335 * t415 + t212 * t768 - t505 * t295 - qJD(1) * t269 + (-t321 - t455) * t475 + t793);
t467 = t690 * qJD(2);
t389 = t683 * t536;
t366 = t475 * t537 + t740;
t296 = (t475 * t561 + t476 * t563) * t536;
t281 = -qJD(2) * t442 + (t563 * t690 + t546) * qJD(1);
t280 = -t771 * t897 + (-t562 * t773 - t735) * rSges(3,1) + t781;
t264 = -t563 * t648 + t840;
t256 = t264 * qJD(1);
t255 = rSges(6,1) * t398 - rSges(6,2) * t399;
t254 = -rSges(6,1) * t396 - rSges(6,2) * t397;
t153 = -t467 * t771 + (-t281 - t471 + t737) * qJD(1);
t152 = -t467 * t535 + t446 + (t280 - t736) * qJD(1);
t125 = t575 * t561 - t563 * t931;
t124 = t561 * t931 + t575 * t563;
t123 = -qJD(2) * t650 + t276 * t562 + t278 * t560;
t122 = -t651 * qJD(2) + t277 * t562 + t279 * t560;
t104 = -t406 * t476 + t454 * t457 + (-t231 + t813) * qJD(1) + t647;
t103 = qJD(1) * t230 - t406 * t475 - t454 * t458 + t587;
t101 = -t361 * t476 + (-t244 + t744) * qJD(1) + t627;
t90 = t256 + t615;
t89 = t616 + t763;
t68 = t230 * t476 + t231 * t475 + t390 * t458 - t391 * t457 + t606;
t64 = -t305 * t476 + t361 * t457 + (-t145 + t620) * qJD(1) + t603;
t63 = qJD(1) * t144 + t458 * t801 + t475 * t812 + t574;
t32 = t145 * t475 + t244 * t458 + (t144 + t213) * t476 - t815 * t457 + t590;
t16 = t117 * t416 + t118 * t415 + t135 * t475 + t210 * t303 - t212 * t302 - t232 * t458 + (t134 + t213) * t476 - t816 * t457 + t590;
t1 = [((t262 + t264) * t563 + (t261 + t263) * t561) * t762 / 0.2e1 + (-(-qJD(1) * t244 - t101 + t602) * t102 - t715 * t476 + t64 * (-t688 + t782) + t101 * (t469 - t689 + t780) + t63 * (t500 + t815) + t102 * (-pkin(3) * t756 + t464 + t643 + t747) + (t64 * t610 + t101 * (t834 * t949 - t519) - t63 * t564) * t561 + ((t536 * t949 - t527 - t905) * t854 + (t101 * (-t527 + t610) - t102 * t564) * t563) * qJD(1)) * m(5) + (-(-qJD(1) * t417 - t257 - t472 - t736) * t258 + t153 * (t561 * t727 + t550 + t779) + t152 * t790 + t258 * (t531 + t781) + (t485 * t851 - t849) * qJD(2) + ((-pkin(1) - t690) * t850 + (t257 * (-rSges(3,3) - pkin(6)) + t258 * t727) * t561) * qJD(1)) * m(3) + t890 / 0.2e1 - t891 / 0.2e1 + t882 / 0.2e1 + t883 / 0.2e1 - (t122 + t125 + t90) * t771 / 0.2e1 + (t123 + t124) * t535 / 0.2e1 + t880 + (t964 + t989) * t915 + (t965 - t990) * t916 + ((t149 + (t379 * t563 + t380 * t561) * t536 + t708 + t810) * t476 + (-t381 * t832 + t848 + t148 + (t379 * t561 - t380 * t563) * t536 + t809 + t100) * t475 + t986) * t911 + t38 * t917 + (t38 + t50) * t918 + (-t975 + t976 + t978) * t912 + t49 * t919 + t121 * t921 + t120 * t922 + (t256 + ((t161 - t363 + (t408 + t845) * t563 + t800) * t563 + t799 * t561) * qJD(2)) * t720 + (t30 * (-t444 - t685 + t782) + t74 * (-t686 + t780) + t29 * (t500 + t212 + t789) + t75 * (t491 + t754 + t787) + (t30 * t904 - t29 * t837 + t75 * (-t522 * t838 - t553 * t833 - t759)) * t563 + (-t29 * t564 + t74 * t881 * t834 + (t30 * t881 + t74 * (t522 * t553 - qJD(4))) * t536) * t561 + ((t637 * t75 - t74 * t904) * t561 + (t74 * (t637 + t837) - t75 * t564) * t563) * qJD(1) - (qJD(1) * t232 + t476 * t811 + t602 - t74 + t954) * t75) * m(6) + ((t706 * t563 + t151 - t621 - t809) * t476 + (t706 * t561 + t150 - t323 + (t378 + t653) * t563 - t98) * t475 + t979 + t985) * t914 + (-qJD(2) * t648 + t462 * t562 + t463 * t560 + (t553 * t654 - t299 + t705) * t537 + (-t300 * t557 + t301 * t558 + t354 * t553 + t704) * t536) * qJD(1) + (-t763 + ((t563 * t703 + t163 - t799) * t563 + (t561 * t703 + t162 + t707) * t561) * qJD(2) + t89) * t723 + (t974 + t977) * t913 + (-(-qJD(1) * t390 - t146 + t614 + t944) * t147 + t104 * (-t390 + t782) + t146 * t780 + t103 * (t391 + t702) + t147 * (-t700 + t786) + (-t147 * t429 + t454 * t852) * t553 + ((-t146 * rSges(4,3) + t147 * (-t527 - t898)) * t561 + (t146 * (-t456 - t527) - t147 * t564) * t563) * qJD(1)) * m(4); (-t122 * t563 + t123 * t561 + (t261 * t561 + t262 * t563) * qJD(1)) * t902 + ((-t771 * t840 - t774) * t563 + (t591 + (t563 * t839 + t572) * qJD(2)) * t561) * t720 + ((t560 * t785 + t562 * t784) * qJD(1) + (t604 * t560 + t562 * t930) * qJD(2)) * t903 + t569 + ((-t535 * t839 + t774) * t561 + (t591 + (t561 * t840 + t572) * qJD(2)) * t563) * t723 + (qJD(1) * t124 + (t561 * (t561 * t934 + t576 * t563) - t563 * (t561 * t933 + t577 * t563) + (t162 * t561 + t163 * t563) * qJD(1)) * t952) * t908 + (qJD(1) * t125 + (t561 * (t576 * t561 - t563 * t934) - t563 * (t577 * t561 - t563 * t933) + (t160 * t561 + t161 * t563) * qJD(1)) * t952) * t907 + (t616 + t89) * t725 + (t615 + t90) * t724 + (t16 * (t638 + t803) + t65 * (t601 + t749) + (t75 * t626 + (-t373 + t751) * t857) * t561 - (-t75 * t738 + (-t562 * t682 + t65 * t714) * qJD(2)) * pkin(2) - t568 + t929 * (t718 - t962) + (-qJD(1) * t268 + t563 * t626 - t775 + t942) * t74) * m(6) + (-t101 * (qJD(1) * t318 + t742 + t775) - (-t102 * t738 + ((-t101 * t563 - t854) * t562 + t84 * t714) * qJD(2)) * pkin(2) - t573 + t101 * t808 + t32 * (t698 + t803) + t84 * (t639 + t749) + (t101 * t642 + (qJD(1) * t102 + t64) * t696) * t563 + (t63 * t696 + t102 * t642 + (-t373 - t815) * t855) * t561) * m(5) + (-(-t147 * t738 + (t128 * t714 + t562 * t662) * qJD(2)) * pkin(2) + t68 * (t802 + t803) + t128 * (t749 + t750) + (t146 * t692 + (qJD(1) * t147 + t104) * t717) * t563 + (t103 * t717 + t147 * t692 + (t146 * t454 + t128 * (-t373 - t391)) * qJD(1)) * t561 + t939) * m(4) + (-(t257 * t442 - t849) * qJD(1) - (t247 * (-t442 * t561 - t443 * t563) + t656 * t690) * qJD(2) + 0.2e1 * t247 * (t280 * t563 + t281 * t561 + (t417 * t563 - t418 * t561) * qJD(1)) + t656 * t467 + (-t152 * t561 - t153 * t563 + (-t258 * t563 + t851) * qJD(1)) * t485) * m(3); t569 + (t16 * t638 + t65 * t601 + (t75 * t753 + t751 * t857) * t561 - t568 + t929 * (-t334 + t811) + (t753 * t563 - (t268 + t425) * qJD(1) + t942) * t74) * m(6) + (-t573 + t32 * t698 + t84 * t639 + (qJD(1) * t715 + t64 * t801) * t563 + (t102 * t812 + t63 * t801 - t815 * t855) * t561 + (-(t318 + t425) * qJD(1) - t742 + t808 + t812 * t563) * t101) * m(5) + (t68 * t802 + t128 * (-t391 * t773 + t750) + t662 * t406 + (-t103 * t561 - t104 * t563 + (-t147 * t563 + t852) * qJD(1)) * t454 + t939) * m(4); -m(5) * (t101 * t943 + t102 * t366 + t296 * t84) - m(6) * (t296 * t65 + t366 * t75 + t74 * t943) + 0.2e1 * ((t101 * t476 + t102 * t830 - t32) * t924 + (t476 * t74 + t75 * t830 - t16) * t923) * t537 + 0.2e1 * ((-t101 * t773 + t102 * t772 + t553 * t84 + t561 * t63 + t563 * t64) * t924 + (t30 * t563 + t553 * t65 - t74 * t773 + t75 * t772 + t950) * t923) * t536; -t38 * t741 / 0.2e1 + t3 * t835 / 0.2e1 + (-t121 * t537 + t536 * t678) * t921 + ((t553 * t678 - t49) * t537 + (-qJD(1) * t679 + t121 * t553 + t21 * t561 + t22 * t563) * t536) * t919 + t536 * t37 * t724 + t4 * t836 / 0.2e1 + (-t120 * t537 + t536 * t680) * t922 + ((t553 * t680 - t50) * t537 + (-qJD(1) * t681 + t120 * t553 + t23 * t561 + t24 * t563) * t536) * t918 + t44 * t729 - t537 * (t880 + t882 + t883 + t890 - t891) / 0.2e1 + (-t133 * t537 + t536 * t676) * t697 + ((t553 * t676 - t57) * t537 + (-qJD(1) * t677 + t133 * t553 + t27 * t561 + t28 * t563) * t536) * t909 + (t398 * t927 + t586 * t399 - t563 * t600) * t920 + (-t396 * t927 + t397 * t586 - t561 * t600) * t917 + (t608 * t537 + (-t532 * t927 + t533 * t586) * t536) * t910 + t955 * t834 / 0.2e1 + ((-t75 * t117 + t74 * t118 + t30 * t210 - t29 * t212 + (t65 * t659 + (t561 * t74 - t563 * t75) * t334) * t553) * t537 + (t74 * (t177 * t561 - t210 * t553) + t75 * (-t177 * t563 + t212 * t553) + t16 * t659 + t65 * (-t117 * t561 + t118 * t563 - t210 * t773 - t212 * t772) + (qJD(1) * t682 - t29 * t563 + t30 * t561) * t334) * t536 - t74 * (-t254 * t505 - t389 * t416) - t75 * (t255 * t505 - t389 * t415) - t65 * (t254 * t415 + t255 * t416)) * m(6);];
tauc = t1(:);
