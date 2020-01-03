% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:54
% EndTime: 2019-12-31 20:18:00
% DurationCPUTime: 56.41s
% Computational Cost: add. (32636->1112), mult. (34675->1471), div. (0->0), fcn. (31008->10), ass. (0->581)
t991 = Icges(3,3) + Icges(4,3);
t522 = qJ(2) + pkin(9);
t493 = sin(t522);
t494 = cos(t522);
t527 = sin(qJ(2));
t530 = cos(qJ(2));
t985 = Icges(3,5) * t530 + Icges(4,5) * t494 - Icges(3,6) * t527 - Icges(4,6) * t493;
t531 = cos(qJ(1));
t990 = t991 * t531;
t528 = sin(qJ(1));
t799 = t528 * t530;
t802 = t527 * t528;
t806 = t494 * t528;
t808 = t493 * t528;
t977 = -Icges(3,5) * t799 - Icges(4,5) * t806 + Icges(3,6) * t802 + Icges(4,6) * t808 + t990;
t986 = t991 * t528 + t985 * t531;
t987 = Icges(3,5) * t527 + Icges(4,5) * t493 + Icges(3,6) * t530 + Icges(4,6) * t494;
t841 = Icges(4,6) * t531;
t319 = Icges(4,4) * t806 - Icges(4,2) * t808 - t841;
t842 = Icges(3,6) * t531;
t362 = Icges(3,4) * t799 - Icges(3,2) * t802 - t842;
t989 = t319 * t493 + t362 * t527;
t854 = Icges(4,4) * t493;
t410 = Icges(4,1) * t494 - t854;
t322 = Icges(4,5) * t528 + t410 * t531;
t855 = Icges(3,4) * t527;
t450 = Icges(3,1) * t530 - t855;
t365 = Icges(3,5) * t528 + t450 * t531;
t988 = -t322 * t806 - t365 * t799;
t407 = Icges(4,2) * t494 + t854;
t477 = Icges(4,4) * t494;
t409 = Icges(4,1) * t493 + t477;
t447 = Icges(3,2) * t530 + t855;
t509 = Icges(3,4) * t530;
t449 = Icges(3,1) * t527 + t509;
t984 = t407 * t493 - t409 * t494 + t447 * t527 - t449 * t530;
t458 = Icges(4,4) * t808;
t848 = Icges(4,5) * t531;
t321 = Icges(4,1) * t806 - t458 - t848;
t473 = Icges(3,4) * t802;
t849 = Icges(3,5) * t531;
t364 = Icges(3,1) * t799 - t473 - t849;
t973 = -t321 * t494 - t364 * t530 + t989;
t983 = t986 * t531 + t988;
t797 = t530 * t531;
t805 = t494 * t531;
t982 = -t321 * t805 - t364 * t797 + t977 * t528;
t929 = t322 * t805 + t365 * t797 + t986 * t528;
t927 = t987 * t531;
t926 = t987 * t528;
t638 = -Icges(4,2) * t493 + t477;
t320 = Icges(4,6) * t528 + t531 * t638;
t639 = -Icges(3,2) * t527 + t509;
t363 = Icges(3,6) * t528 + t531 * t639;
t981 = t320 * t493 + t363 * t527;
t980 = -t973 * t528 + t977 * t531;
t968 = -t320 * t808 - t363 * t802 - t983;
t801 = t527 * t531;
t807 = t493 * t531;
t967 = -t319 * t807 - t362 * t801 - t982;
t966 = -t320 * t807 - t363 * t801 + t929;
t979 = -t984 * t528 - t927;
t978 = -t984 * t531 + t926;
t944 = t319 * t494 + t321 * t493 + t362 * t530 + t364 * t527;
t943 = t320 * t494 + t322 * t493 + t363 * t530 + t365 * t527;
t975 = t987 * qJD(2);
t974 = t322 * t494 + t365 * t530 - t981;
t495 = qJ(4) + t522;
t479 = cos(t495);
t529 = cos(qJ(5));
t798 = t529 * t531;
t526 = sin(qJ(5));
t803 = t526 * t528;
t376 = t479 * t803 + t798;
t795 = t531 * t526;
t800 = t528 * t529;
t377 = t479 * t800 - t795;
t478 = sin(t495);
t812 = t478 * t528;
t186 = Icges(6,5) * t377 - Icges(6,6) * t376 + Icges(6,3) * t812;
t356 = Icges(6,4) * t377;
t189 = -Icges(6,2) * t376 + Icges(6,6) * t812 + t356;
t355 = Icges(6,4) * t376;
t193 = -Icges(6,1) * t377 - Icges(6,5) * t812 + t355;
t957 = t189 * t526 + t193 * t529;
t84 = -t186 * t479 - t957 * t478;
t972 = rSges(3,2) * t527;
t382 = t638 * qJD(2);
t383 = t410 * qJD(2);
t420 = t639 * qJD(2);
t421 = t450 * qJD(2);
t971 = -t382 * t493 + t383 * t494 - t420 * t527 + t421 * t530 + (-t407 * t494 - t409 * t493 - t447 * t530 - t449 * t527) * qJD(2) + t987 * qJD(1);
t970 = t986 * qJD(1);
t969 = t984 * qJD(1) + t985 * qJD(2);
t965 = t978 * qJD(1);
t593 = qJD(2) * t407;
t203 = -t531 * t593 + (-t528 * t638 + t841) * qJD(1);
t595 = qJD(2) * t409;
t205 = -t531 * t595 + (-t410 * t528 + t848) * qJD(1);
t594 = qJD(2) * t447;
t243 = -t531 * t594 + (-t528 * t639 + t842) * qJD(1);
t596 = qJD(2) * t449;
t245 = -t531 * t596 + (-t450 * t528 + t849) * qJD(1);
t964 = -t943 * qJD(2) - t203 * t493 + t205 * t494 - t243 * t527 + t245 * t530 + t970;
t204 = qJD(1) * t320 - t528 * t593;
t206 = qJD(1) * t322 - t528 * t595;
t244 = qJD(1) * t363 - t528 * t594;
t246 = qJD(1) * t365 - t528 * t596;
t963 = t977 * qJD(1) + t944 * qJD(2) + t204 * t493 - t206 * t494 + t244 * t527 - t246 * t530;
t962 = (t966 * t528 - t967 * t531) * qJD(2);
t961 = (t968 * t528 - t980 * t531) * qJD(2);
t960 = t979 * qJD(1);
t959 = t973 * qJD(1) - t975 * t528 + t970;
t958 = -t975 * t531 + (-t528 * t985 - t974 + t990) * qJD(1);
t634 = Icges(6,5) * t529 - Icges(6,6) * t526;
t274 = -Icges(6,3) * t479 + t478 * t634;
t850 = Icges(6,4) * t529;
t636 = -Icges(6,2) * t526 + t850;
t276 = -Icges(6,6) * t479 + t478 * t636;
t851 = Icges(6,4) * t526;
t641 = Icges(6,1) * t529 - t851;
t278 = -Icges(6,5) * t479 + t478 * t641;
t105 = t274 * t812 - t276 * t376 + t278 * t377;
t499 = qJD(2) * t528;
t443 = qJD(4) * t528 + t499;
t736 = qJD(5) * t531;
t352 = t478 * t736 + t443;
t521 = qJD(2) + qJD(4);
t444 = t521 * t531;
t737 = qJD(5) * t528;
t353 = -t478 * t737 + t444;
t738 = qJD(5) * t479;
t440 = qJD(1) - t738;
t76 = t186 * t812 - t189 * t376 - t193 * t377;
t378 = -t479 * t795 + t800;
t379 = t479 * t798 + t803;
t811 = t478 * t531;
t188 = Icges(6,5) * t379 + Icges(6,6) * t378 + Icges(6,3) * t811;
t852 = Icges(6,4) * t379;
t191 = Icges(6,2) * t378 + Icges(6,6) * t811 + t852;
t357 = Icges(6,4) * t378;
t194 = Icges(6,1) * t379 + Icges(6,5) * t811 + t357;
t77 = t188 * t812 - t376 * t191 + t377 * t194;
t29 = t105 * t440 + t352 * t77 - t353 * t76;
t106 = t274 * t811 + t276 * t378 + t278 * t379;
t78 = t186 * t811 + t378 * t189 - t193 * t379;
t79 = t188 * t811 + t378 * t191 + t379 * t194;
t30 = t106 * t440 + t352 * t79 - t353 * t78;
t956 = 0.2e1 * qJD(2);
t955 = t960 + t961;
t954 = t962 + t965;
t953 = t973 * qJD(2) - t204 * t494 - t206 * t493 - t244 * t530 - t246 * t527;
t952 = t974 * qJD(2) + t203 * t494 + t205 * t493 + t243 * t530 + t245 * t527;
t951 = t969 * t528 + t971 * t531;
t950 = t971 * t528 - t969 * t531;
t741 = qJD(1) * t531;
t804 = t521 * t528;
t581 = t478 * t741 + t479 * t804;
t616 = t440 * t529;
t667 = qJD(1) * t479 - qJD(5);
t727 = t478 * t804;
t179 = t528 * t616 + (-t531 * t667 + t727) * t526;
t617 = t440 * t526;
t813 = t478 * t521;
t180 = t667 * t798 + (-t529 * t813 + t617) * t528;
t656 = rSges(6,1) * t180 + rSges(6,2) * t179;
t104 = rSges(6,3) * t581 + t656;
t654 = rSges(6,1) * t529 - rSges(6,2) * t526;
t603 = t654 * t479;
t653 = -rSges(6,1) * t526 - rSges(6,2) * t529;
t154 = t521 * t603 + (rSges(6,3) * t521 + qJD(5) * t653) * t478;
t655 = rSges(6,1) * t377 - rSges(6,2) * t376;
t195 = rSges(6,3) * t812 + t655;
t199 = t581 * pkin(8) + (t479 * t741 - t727) * pkin(4);
t732 = qJD(1) * qJD(2);
t484 = t528 * t732;
t730 = qJD(1) * qJD(4);
t417 = t528 * t730 + t484;
t251 = qJD(5) * t581 + t417;
t282 = -rSges(6,3) * t479 + t478 * t654;
t880 = t479 * pkin(4);
t394 = pkin(8) * t478 + t880;
t331 = t394 * t521;
t881 = t478 * pkin(4);
t393 = -pkin(8) * t479 + t881;
t883 = pkin(3) * t493;
t435 = t499 * t883;
t517 = t530 * pkin(2);
t532 = qJD(2) ^ 2;
t882 = pkin(3) * t494;
t609 = (-t517 - t882) * t532;
t884 = pkin(2) * t527;
t467 = t499 * t884;
t731 = qJD(1) * qJD(3);
t753 = qJD(1) * t467 + t531 * t731;
t563 = qJD(1) * t435 + t531 * t609 + t753;
t739 = qJD(5) * t478;
t704 = t521 * t739;
t525 = -qJ(3) - pkin(6);
t742 = qJD(1) * t528;
t470 = t525 * t742;
t428 = -t883 - t884;
t411 = t428 * qJD(2);
t520 = -pkin(7) + t525;
t466 = t520 * t742;
t672 = t411 * t528 - t466;
t488 = t517 + pkin(1);
t422 = t488 + t882;
t759 = t422 - t488;
t147 = t741 * t759 + t467 + t470 + t672;
t516 = t528 * pkin(6);
t498 = qJD(3) * t531;
t750 = t467 + t498;
t715 = t470 + t750;
t877 = pkin(1) - t488;
t230 = (-t531 * t877 - t516) * qJD(1) - t715;
t462 = t531 * pkin(1) + t516;
t425 = t462 * qJD(1);
t790 = -t230 - t425;
t724 = -t147 + t790;
t23 = -t195 * t704 - t104 * t440 - t154 * t353 + t251 * t282 - t331 * t444 + t393 * t417 + (-t199 + t724) * qJD(1) + t563;
t197 = t379 * rSges(6,1) + t378 * rSges(6,2) + rSges(6,3) * t811;
t809 = t479 * t531;
t351 = pkin(4) * t809 + pkin(8) * t811;
t716 = -t435 - t750;
t396 = t531 * t422;
t464 = t531 * t488;
t748 = -t520 + t525;
t254 = t528 * t748 + t396 - t464;
t671 = -t525 * t528 + t464;
t315 = t671 - t462;
t770 = t315 + t462;
t717 = t254 + t770;
t66 = t197 * t440 - t282 * t352 - t393 * t443 + (t351 + t717) * qJD(1) + t716;
t833 = qJD(1) * t66;
t946 = t531 * (t23 + t833);
t327 = rSges(4,1) * t806 - rSges(4,2) * t808 - t531 * rSges(4,3);
t511 = t528 * rSges(4,3);
t328 = rSges(4,1) * t805 - rSges(4,2) * t807 + t511;
t623 = t327 * t528 + t328 * t531;
t518 = t531 * pkin(6);
t461 = pkin(1) * t528 - t518;
t486 = t531 * t525;
t752 = t528 * t488 + t486;
t314 = t461 - t752;
t740 = qJD(2) * t531;
t779 = -t314 * t499 + t315 * t740;
t113 = qJD(2) * t623 + t779;
t778 = -t528 * t314 + t531 * t315;
t584 = t623 + t778;
t945 = qJD(2) * t584 + t113;
t349 = t394 * t528;
t942 = t195 + t349;
t941 = t197 + t351;
t940 = t282 + t393;
t810 = t479 * t528;
t840 = Icges(5,6) * t531;
t301 = Icges(5,4) * t810 - Icges(5,2) * t812 - t840;
t469 = Icges(5,4) * t479;
t389 = Icges(5,1) * t478 + t469;
t939 = -t389 * t528 - t301;
t637 = -Icges(5,2) * t478 + t469;
t302 = Icges(5,6) * t528 + t531 * t637;
t938 = -t389 * t531 - t302;
t853 = Icges(5,4) * t478;
t390 = Icges(5,1) * t479 - t853;
t304 = Icges(5,5) * t528 + t390 * t531;
t387 = Icges(5,2) * t479 + t853;
t937 = -t387 * t531 + t304;
t936 = -t387 + t390;
t935 = t389 + t637;
t573 = t362 * t531 - t363 * t528;
t574 = t319 * t531 - t320 * t528;
t905 = t528 * (-t407 * t531 + t322) - t531 * (-Icges(4,2) * t806 + t321 - t458);
t906 = t528 * (-t447 * t531 + t365) - t531 * (-Icges(3,2) * t799 + t364 - t473);
t934 = -t493 * t905 + t574 * t494 - t527 * t906 + t573 * t530;
t754 = t449 + t639;
t755 = -t447 + t450;
t760 = t409 + t638;
t761 = -t407 + t410;
t933 = (-t493 * t760 + t494 * t761 - t527 * t754 + t530 * t755) * qJD(1);
t932 = t977 + t981;
t931 = t985 * qJD(1);
t689 = t428 + t884;
t334 = t689 * t528;
t713 = t493 * t742;
t930 = pkin(3) * t713 + qJD(1) * t334;
t261 = t282 * t528;
t872 = rSges(6,3) * t478;
t283 = t603 + t872;
t348 = t393 * t528;
t703 = t479 * t737;
t928 = -qJD(1) * t348 + t195 * t739 - t261 * t440 - t282 * t703 + t353 * t283 + t444 * t394 + t940 * t742;
t925 = t29 * t528 + t30 * t531;
t497 = qJD(3) * t528;
t796 = t531 * t411;
t920 = t796 + t497;
t924 = -t195 * t440 - t353 * t282 - t444 * t393 + t920;
t429 = qJD(1) * t461;
t921 = qJD(1) * t314 - t429;
t762 = -t528 * t422 - t531 * t520;
t253 = t752 + t762;
t917 = qJD(1) * t253 + t921;
t726 = t478 * t444;
t916 = t528 * t667 + t726;
t386 = Icges(5,5) * t479 - Icges(5,6) * t478;
t385 = Icges(5,5) * t478 + Icges(5,6) * t479;
t818 = t385 * t531;
t826 = t302 * t478;
t836 = Icges(5,3) * t531;
t915 = -t521 * t818 + (-t304 * t479 - t386 * t528 + t826 + t836) * qJD(1);
t434 = Icges(5,4) * t812;
t847 = Icges(5,5) * t531;
t303 = Icges(5,1) * t810 - t434 - t847;
t627 = t301 * t478 - t303 * t479;
t300 = Icges(5,3) * t528 + t386 * t531;
t747 = qJD(1) * t300;
t819 = t385 * t528;
t914 = qJD(1) * t627 - t521 * t819 + t747;
t620 = t387 * t478 - t389 * t479;
t909 = qJD(1) * t620 + t386 * t521;
t598 = t634 * t479;
t628 = -t276 * t526 + t278 * t529;
t631 = -t191 * t526 + t194 * t529;
t904 = t352 * (-t274 * t531 - t631) - t353 * (-t274 * t528 + t957) + t440 * (Icges(6,3) * t478 + t598 - t628);
t635 = -Icges(6,2) * t529 - t851;
t903 = t352 * (-Icges(6,2) * t379 + t194 + t357) - t353 * (-Icges(6,2) * t377 - t193 - t355) + t440 * (t635 * t478 + t278);
t902 = qJD(1) * t935 + t443 * t937 - t444 * (-Icges(5,2) * t810 + t303 - t434);
t901 = t251 / 0.2e1;
t485 = t531 * t732;
t418 = t531 * t730 + t485;
t714 = t478 * t742;
t725 = t479 * t444;
t580 = -t714 + t725;
t252 = qJD(5) * t580 + t418;
t900 = t252 / 0.2e1;
t899 = -t352 / 0.2e1;
t898 = t352 / 0.2e1;
t897 = -t353 / 0.2e1;
t896 = t353 / 0.2e1;
t895 = t417 / 0.2e1;
t894 = t418 / 0.2e1;
t893 = -t440 / 0.2e1;
t892 = t440 / 0.2e1;
t891 = -t443 / 0.2e1;
t890 = t443 / 0.2e1;
t889 = -t444 / 0.2e1;
t888 = t444 / 0.2e1;
t887 = t528 / 0.2e1;
t886 = -t531 / 0.2e1;
t885 = -rSges(6,3) - pkin(8);
t879 = -qJD(1) / 0.2e1;
t878 = qJD(1) / 0.2e1;
t876 = rSges(3,1) * t530;
t875 = rSges(4,1) * t494;
t874 = rSges(5,1) * t479;
t873 = rSges(4,2) * t494;
t100 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t581;
t102 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t581;
t98 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t581;
t20 = (-t521 * t957 - t98) * t479 + (-t100 * t526 + t102 * t529 + t186 * t521 + (-t189 * t529 + t193 * t526) * qJD(5)) * t478;
t870 = t20 * t353;
t177 = t526 * t916 + t531 * t616;
t178 = -t529 * t916 + t531 * t617;
t101 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t580;
t97 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t580;
t99 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t580;
t21 = (t521 * t631 - t97) * t479 + (t101 * t529 + t188 * t521 - t526 * t99 + (-t191 * t529 - t194 * t526) * qJD(5)) * t478;
t869 = t21 * t352;
t722 = t178 * rSges(6,1) + t177 * rSges(6,2) + rSges(6,3) * t725;
t103 = -rSges(6,3) * t714 + t722;
t416 = pkin(8) * t725;
t582 = -t479 * t742 - t726;
t198 = pkin(4) * t582 - pkin(8) * t714 + t416;
t708 = t527 * t740;
t668 = pkin(2) * t708;
t146 = t668 + t796 + (-t528 * t759 + t531 * t748) * qJD(1);
t492 = pkin(6) * t741;
t229 = -t668 - t492 + t497 + (t528 * t877 - t486) * qJD(1);
t415 = qJD(1) * (-pkin(1) * t742 + t492);
t719 = qJD(1) * t229 + t528 * t731 + t415;
t561 = qJD(1) * t146 + t528 * t609 + t719;
t22 = t197 * t704 + t103 * t440 - t154 * t352 - t252 * t282 - t331 * t443 - t393 * t418 + (t198 + t796) * qJD(1) + t561;
t868 = t22 * t531;
t867 = t23 * t528;
t512 = t528 * rSges(3,3);
t510 = t528 * rSges(5,3);
t862 = t84 * t251;
t85 = -t188 * t479 + t478 * t631;
t861 = t85 * t252;
t391 = rSges(5,1) * t478 + rSges(5,2) * t479;
t306 = rSges(5,1) * t810 - rSges(5,2) * t812 - t531 * rSges(5,3);
t562 = -t444 * t391 + t920;
t772 = t314 - t461;
t718 = t253 + t772;
t92 = (-t306 + t718) * qJD(1) + t562;
t860 = t92 * t391;
t112 = -t274 * t479 + t478 * t628;
t633 = -Icges(6,5) * t526 - Icges(6,6) * t529;
t151 = t521 * t598 + (Icges(6,3) * t521 + qJD(5) * t633) * t478;
t599 = t636 * t479;
t152 = t521 * t599 + (Icges(6,6) * t521 + qJD(5) * t635) * t478;
t600 = t641 * t479;
t640 = -Icges(6,1) * t526 - t850;
t153 = t521 * t600 + (Icges(6,5) * t521 + qJD(5) * t640) * t478;
t41 = (t521 * t628 - t151) * t479 + (-t152 * t526 + t153 * t529 + t274 * t521 + (-t276 * t529 - t278 * t526) * qJD(5)) * t478;
t859 = t112 * t704 + t41 * t440;
t665 = -t253 * t499 + t254 * t740 + t779;
t48 = t195 * t352 + t197 * t353 + t349 * t443 + t351 * t444 + t665;
t834 = qJD(1) * t48;
t412 = rSges(4,1) * t493 + t873;
t691 = -t412 - t884;
t659 = t531 * t691;
t612 = qJD(2) * t659;
t579 = t497 + t612;
t120 = (-t327 + t772) * qJD(1) + t579;
t831 = t120 * t412;
t749 = rSges(3,2) * t802 + t531 * rSges(3,3);
t374 = rSges(3,1) * t799 - t749;
t452 = rSges(3,1) * t527 + rSges(3,2) * t530;
t709 = t452 * t740;
t214 = -t709 + (-t374 - t461) * qJD(1);
t830 = t214 * t528;
t829 = t214 * t531;
t710 = t452 * t499;
t375 = rSges(3,1) * t797 - rSges(3,2) * t801 + t512;
t765 = t375 + t462;
t215 = qJD(1) * t765 - t710;
t404 = t452 * t531;
t828 = t215 * t404;
t299 = Icges(5,5) * t810 - Icges(5,6) * t812 - t836;
t827 = t299 * t531;
t794 = -t154 - t331;
t789 = -t254 - t315;
t788 = -t528 * t299 - t303 * t809;
t787 = t528 * t300 + t304 * t809;
t307 = rSges(5,1) * t809 - rSges(5,2) * t811 + t510;
t781 = t528 * t306 + t531 * t307;
t335 = t689 * t531;
t777 = t334 * t499 + t335 * t740;
t758 = rSges(5,2) * t714 + rSges(5,3) * t741;
t757 = rSges(4,2) * t713 + rSges(4,3) * t741;
t751 = rSges(3,3) * t741 + t742 * t972;
t135 = -t528 * t620 - t818;
t735 = t135 * qJD(1);
t729 = pkin(2) * t801;
t728 = qJD(2) * t517;
t175 = rSges(5,1) * t582 - rSges(5,2) * t725 + t758;
t346 = t391 * t528;
t392 = -rSges(5,2) * t478 + t874;
t176 = -t521 * t346 + (t392 * t531 + t510) * qJD(1);
t723 = t531 * t175 + t528 * t176 + t306 * t741;
t721 = t229 * t740 + t230 * t499 - t314 * t485;
t720 = t531 * t229 + t528 * t230 - t314 * t741;
t711 = t412 * t499;
t707 = t530 * t740;
t702 = t479 * t736;
t701 = t813 / 0.2e1;
t699 = -pkin(1) - t876;
t698 = t742 / 0.2e1;
t697 = t741 / 0.2e1;
t696 = -t499 / 0.2e1;
t693 = t740 / 0.2e1;
t413 = -rSges(4,2) * t493 + t875;
t690 = -t413 - t517;
t687 = t527 * (-t528 ^ 2 - t531 ^ 2);
t686 = (-t528 * t637 + t840) * qJD(1) + t937 * t521;
t685 = qJD(1) * t302 + t303 * t521 - t387 * t804;
t684 = (-t390 * t528 + t847) * qJD(1) + t938 * t521;
t683 = qJD(1) * t304 + t521 * t939;
t263 = t304 * t810;
t682 = t300 * t531 - t263;
t347 = t391 * t531;
t681 = -t443 * t346 - t347 * t444;
t679 = -t299 + t826;
t676 = t935 * t521;
t675 = t936 * t521;
t674 = -qJD(1) * t347 - t392 * t443;
t670 = qJD(1) * t346 - t444 * t392;
t664 = -t528 * t253 + t531 * t254 + t778;
t663 = t528 * t942 + t531 * t941;
t662 = qJD(5) * t701;
t657 = t876 - t972;
t65 = (-t349 + t718) * qJD(1) + t924;
t652 = t528 * t66 + t531 * t65;
t651 = t528 * t77 - t531 * t76;
t650 = t528 * t76 + t531 * t77;
t649 = t528 * t79 - t531 * t78;
t648 = t528 * t78 + t531 * t79;
t647 = t528 * t85 - t531 * t84;
t646 = t528 * t84 + t531 * t85;
t93 = -t391 * t443 + (t307 + t717) * qJD(1) + t716;
t645 = -t528 * t93 - t531 * t92;
t630 = t195 * t531 - t197 * t528;
t629 = -t215 * t528 - t829;
t148 = t301 * t479 + t303 * t478;
t615 = (t335 - t729) * qJD(1);
t614 = -t391 + t428;
t613 = t942 * t741 + (t103 + t198) * t531 + (t104 + t199) * t528;
t611 = t531 * t146 + t528 * t147 - t253 * t741 + t720;
t384 = t413 * qJD(2);
t610 = -qJD(2) * t384 - t517 * t532;
t403 = t452 * t528;
t372 = t412 * t528;
t601 = t428 - t940;
t597 = t627 * t528;
t210 = (t374 * t528 + t375 * t531) * qJD(2);
t586 = -qJD(2) * t882 - t728;
t585 = -t394 - t872;
t578 = -t186 * t353 + t188 * t352 + t274 * t440;
t577 = (-Icges(6,5) * t376 - Icges(6,6) * t377) * t353 - (Icges(6,5) * t378 - Icges(6,6) * t379) * t352 - t633 * t478 * t440;
t329 = t392 * t521;
t576 = -t329 + t586;
t575 = qJD(1) * t386 - t443 * t818 + t444 * t819;
t571 = t586 + t794;
t570 = t478 * t577;
t559 = (Icges(6,1) * t378 - t191 - t852) * t352 - (-Icges(6,1) * t376 - t189 - t356) * t353 + (t640 * t478 - t276) * t440;
t557 = t146 * t740 + t147 * t499 - t253 * t485 + t484 * t789 + t721;
t555 = qJD(1) * t936 + t443 * t938 - t444 * t939;
t262 = t282 * t531;
t350 = t393 * t531;
t554 = t195 * t702 - t197 * t703 - t352 * t261 - t262 * t353 - t443 * t348 - t350 * t444;
t553 = -qJD(1) * t350 + t197 * t739 - t440 * t262 - t282 * t702 - t283 * t352 - t394 * t443;
t547 = qJD(1) * t299 - t478 * t685 + t479 * t683;
t546 = -t478 * t686 + t479 * t684 + t747;
t545 = qJD(1) * t385 - t478 * t676 + t479 * t675;
t114 = -t597 - t827;
t115 = -t302 * t812 - t682;
t116 = -t301 * t811 - t788;
t117 = -t302 * t811 + t787;
t149 = t302 * t479 + t304 * t478;
t16 = t100 * t378 + t102 * t379 + t177 * t189 - t178 * t193 + t186 * t580 + t811 * t98;
t17 = t101 * t379 + t177 * t191 + t178 * t194 + t188 * t580 + t378 * t99 + t811 * t97;
t18 = -t100 * t376 + t102 * t377 + t179 * t189 - t180 * t193 + t186 * t581 + t812 * t98;
t19 = t101 * t377 + t179 * t191 + t180 * t194 + t188 * t581 - t376 * t99 + t812 * t97;
t257 = t276 * t528;
t258 = t276 * t531;
t259 = t278 * t528;
t260 = t278 * t531;
t277 = Icges(6,6) * t478 + t599;
t279 = Icges(6,5) * t478 + t600;
t34 = t151 * t811 + t152 * t378 + t153 * t379 + t177 * t276 + t178 * t278 + t274 * t580;
t3 = t106 * t704 - t16 * t353 + t17 * t352 + t251 * t78 + t252 * t79 + t34 * t440;
t33 = t112 * t440 + t352 * t85 - t353 * t84;
t35 = t151 * t812 - t152 * t376 + t153 * t377 + t179 * t276 + t180 * t278 + t274 * t581;
t4 = t105 * t704 - t18 * t353 + t19 * t352 + t251 * t76 + t252 * t77 + t35 * t440;
t43 = t528 * t914 + t547 * t531;
t44 = t528 * t915 + t546 * t531;
t45 = t547 * t528 - t531 * t914;
t46 = t546 * t528 - t531 * t915;
t534 = -t478 * t902 + t555 * t479;
t535 = t904 * t478;
t55 = -t114 * t444 + t115 * t443 + t735;
t136 = -t531 * t620 + t819;
t128 = t136 * qJD(1);
t56 = -t116 * t444 + t117 * t443 + t128;
t69 = t528 * t909 + t545 * t531;
t70 = t545 * t528 - t531 * t909;
t80 = t478 * t683 + t479 * t685;
t81 = t478 * t684 + t479 * t686;
t536 = t649 * t900 + t651 * t901 + (t528 * t534 - t531 * t575) * t888 + (-t45 * t531 + t46 * t528 + (t114 * t528 + t115 * t531) * qJD(1)) * t889 + (-t43 * t531 + t44 * t528 + (t116 * t528 + t117 * t531) * qJD(1)) * t890 + (t528 * t575 + t531 * t534) * t891 + (qJD(1) * t646 - t20 * t531 + t21 * t528) * t892 + (-t116 * t531 + t117 * t528) * t894 + (-t114 * t531 + t115 * t528) * t895 + ((t258 * t376 - t260 * t377) * t352 - (t257 * t376 - t259 * t377) * t353 + (-t277 * t376 + t279 * t377) * t440 + (t105 * t478 + t77 * t809) * qJD(5) + ((qJD(5) * t76 + t578) * t479 + t535) * t528) * t896 + (qJD(1) * t650 - t18 * t531 + t19 * t528) * t897 + (qJD(1) * t648 - t16 * t531 + t17 * t528) * t898 + ((-t258 * t378 - t260 * t379) * t352 - (-t257 * t378 - t259 * t379) * t353 + (t277 * t378 + t279 * t379) * t440 + (t106 * t478 + t78 * t810) * qJD(5) + ((qJD(5) * t79 + t578) * t479 + t535) * t531) * t899 + (t528 * t81 - t531 * t80 + (t148 * t528 + t149 * t531) * qJD(1)) * t878 - t33 * t739 / 0.2e1 + (t555 * t478 + t479 * t902) * t879 + (((t258 * t526 - t260 * t529 + t188) * t352 - (t257 * t526 - t259 * t529 + t186) * t353 + (-t277 * t526 + t279 * t529 + t274) * t440 + t112 * qJD(5)) * t478 + (qJD(5) * t646 - t904) * t479) * t893 + t647 * t662 + (qJD(1) * t69 + t116 * t417 + t117 * t418 - t43 * t444 + t44 * t443 + t3) * t887 + (qJD(1) * t70 + t114 * t417 + t115 * t418 + t443 * t46 - t444 * t45 + t4) * t886 + (t55 + t29) * t698 + (t56 + t30) * t697 - t925 * t738 / 0.2e1;
t423 = t657 * qJD(2);
t373 = t412 * t531;
t345 = t653 * t478;
t250 = -qJD(2) * t403 + (t531 * t657 + t512) * qJD(1);
t249 = -rSges(3,2) * t707 + (-t530 * t742 - t708) * rSges(3,1) + t751;
t240 = rSges(6,1) * t378 - rSges(6,2) * t379;
t239 = -rSges(6,1) * t376 - rSges(6,2) * t377;
t208 = -qJD(2) * t372 + (t413 * t531 + t511) * qJD(1);
t207 = -t740 * t873 + (-t493 * t740 - t494 * t742) * rSges(4,1) + t757;
t127 = -t423 * t740 + (-t250 - t425 + t710) * qJD(1);
t126 = -t423 * t499 + t415 + (t249 - t709) * qJD(1);
t121 = -t711 + (t328 + t770) * qJD(1) - t750;
t87 = t610 * t531 + (-t208 + t711 + t790) * qJD(1) + t753;
t86 = t610 * t528 + (t207 + t612) * qJD(1) + t719;
t73 = t306 * t443 + t307 * t444 + t665;
t58 = -t329 * t444 + t391 * t417 + (-t176 + t724) * qJD(1) + t563;
t57 = -t329 * t443 - t391 * t418 + (t175 + t796) * qJD(1) + t561;
t28 = t175 * t444 + t176 * t443 + t306 * t418 - t307 * t417 + t557;
t13 = t103 * t353 + t104 * t352 + t195 * t252 - t197 * t251 + t198 * t444 + t199 * t443 + t349 * t418 - t351 * t417 + t557;
t1 = [(t23 * (-t655 + t762) + t22 * (t396 + t941) + (t23 * t585 - t22 * t520 + (t478 * t885 - t422 - t880) * t833) * t528 + (t466 + t498 - t656 + (-t422 + t585) * t741 + (-t411 + (t479 * t885 + t881) * t521) * t528) * t65 + (qJD(1) * t349 + t65 - t917 + t416 + t722 + t920 + (-pkin(4) * t813 - t520 * qJD(1)) * t531 - t924) * t66) * m(6) - (t950 - t953 + t954) * t740 / 0.2e1 + (t135 + t148) * t895 + (t136 + t149) * t894 + (t951 + t952) * t499 / 0.2e1 + (t55 - t735 + (t117 - t597 - t787) * t444 + (t528 * t679 + t116 - t263) * t443 + ((t300 + t627) * t443 + t679 * t444) * t531) * t891 + (t80 + t70 + t56) * t889 + t106 * t900 + t105 * t901 + (t128 + (t115 + (t301 * t531 + t302 * t528) * t478 + t682 + t788) * t444 + (-t303 * t810 + t827 + t114 + (t301 * t528 - t302 * t531) * t478 + t787) * t443) * t888 + t34 * t898 + (t35 + t30) * t897 + (t81 + t69) * t890 + ((t943 + t978) * t531 + (t944 + t979) * t528) * t732 / 0.2e1 + t30 * t896 + t862 / 0.2e1 + (t87 * (-t327 - t752) + t120 * t715 + t86 * (t328 + t671) + t121 * (t497 + t757) + (t121 * t659 + t528 * t831) * qJD(2) + ((-t120 * rSges(4,3) + t121 * (-t488 - t875)) * t528 + (t120 * (-t413 - t488) - t121 * t525) * t531) * qJD(1) - (-qJD(1) * t327 - t120 + t579 + t921) * t121) * m(4) + (t58 * (-t306 + t762) + t92 * (t498 - t672) + t57 * (-t520 * t528 + t307 + t396) + t93 * (t758 + t920) + (-t347 * t93 + t528 * t860) * t521 + ((-t92 * rSges(5,3) + t93 * (-t422 - t874)) * t528 + (t92 * (-t392 - t422) - t93 * t520) * t531) * qJD(1) - (-qJD(1) * t306 + t562 + t917 - t92) * t93) * m(5) + (((t531 * t932 - t929 + t966) * t531 + (t528 * t932 + t967 + t983) * t528) * qJD(2) + t955 - t960) * t696 + (-t984 * qJD(2) + t382 * t494 + t383 * t493 + t420 * t530 + t421 * t527 + t478 * t675 + t479 * t676) * qJD(1) + (-(-qJD(1) * t374 - t214 - t429 - t709) * t215 + t127 * (t528 * t699 + t518 + t749) + t126 * t765 + t215 * (t492 + t751) + (t452 * t830 - t828) * qJD(2) + ((-pkin(1) - t657) * t829 + (t214 * (-rSges(3,3) - pkin(6)) + t215 * t699) * t528) * qJD(1)) * m(3) + t859 + t861 / 0.2e1 - t870 / 0.2e1 + t869 / 0.2e1 + ((t929 * t528 + ((t986 + t989) * t531 + t968 + t982 + t988) * t531) * qJD(2) + t965) * t693; t536 + ((t574 * t493 + t494 * t905 + t573 * t527 + t530 * t906) * qJD(2) + (t493 * t761 + t494 * t760 + t527 * t755 + t530 * t754) * qJD(1)) * t879 + (t953 * t531 + t952 * t528 + (t528 * t944 + t531 * t943) * qJD(1)) * t878 + ((-t499 * t927 + t931) * t528 + ((t528 * t926 + t934) * qJD(2) + t933) * t531) * t696 + ((-t740 * t926 - t931) * t531 + ((t531 * t927 + t934) * qJD(2) + t933) * t528) * t693 + (-t66 * (t553 + t615) - t48 * (t554 + t777) - (-t652 * t882 + (t48 * t687 - t530 * t652) * pkin(2)) * qJD(2) + t13 * (t663 + t664) + t48 * (t611 + t613) + t601 * t946 + (t22 * t601 + t66 * t571 + (t789 - t941) * t834) * t528 + (t531 * t571 + t928 + t930) * t65) * m(6) + (t28 * (t664 + t781) + t73 * (t611 + t723) + (qJD(1) * t93 + t58) * t614 * t531 + (t57 * t614 + t93 * t576 + (t860 + t73 * (-t307 + t789)) * qJD(1)) * t528 - t93 * (t615 + t674) - t73 * (t681 + t777) - (t645 * t882 + (t530 * t645 + t687 * t73) * pkin(2)) * qJD(2) + (t531 * t576 - t670 + t930) * t92) * m(5) + (t87 * t659 - t120 * pkin(2) * t707 + (t207 * t740 + t208 * t499 + t721) * t584 + t113 * t720 + (-t120 * t384 + t113 * t207 + (t121 * t691 + t327 * t945) * qJD(1)) * t531 + (t86 * t691 + t121 * (-t384 - t728) + t113 * t208 + (t831 + t945 * (-t315 - t328)) * qJD(1)) * t528 - (t120 * t372 + t121 * (-t373 - t729)) * qJD(1) - (t113 * pkin(2) * t687 + (-t113 * t373 + t120 * t690) * t531 + (-t113 * t372 + t121 * t690) * t528) * qJD(2)) * m(4) + (0.2e1 * t210 * (t249 * t531 + t250 * t528 + (t374 * t531 - t375 * t528) * qJD(1)) + t629 * t423 + (-t126 * t528 - t127 * t531 + (-t215 * t531 + t830) * qJD(1)) * t452 - (t214 * t403 - t828) * qJD(1) - (t210 * (-t403 * t528 - t404 * t531) + t629 * t657) * qJD(2)) * m(3) + (t951 * qJD(1) + ((t966 * qJD(1) + t963 * t531) * t531 + (t958 * t528 + t967 * qJD(1) + (-t959 + t964) * t531) * t528) * t956) * t887 + (t950 * qJD(1) + ((t968 * qJD(1) + t959 * t531) * t531 + (t964 * t528 + t980 * qJD(1) + (-t958 + t963) * t531) * t528) * t956) * t886 + (t955 + t961) * t698 + (t954 + t962) * t697; 0.2e1 * (-t868 / 0.2e1 + t867 / 0.2e1) * m(6) + 0.2e1 * (t57 * t886 + t58 * t887) * m(5) + 0.2e1 * (t86 * t886 + t87 * t887) * m(4); t536 + (-t553 * t66 + t13 * t663 + (t66 * t794 - t834 * t941) * t528 - (t22 * t528 + t946) * t940 + (t531 * t794 + t928) * t65 + (-t554 + t613) * t48) * m(6) + (-t670 * t92 - t674 * t93 + t28 * t781 + t645 * t329 + (-t57 * t528 - t58 * t531 + (t528 * t92 - t531 * t93) * qJD(1)) * t391 + (-t307 * t742 - t681 + t723) * t73) * m(5); -t30 * t714 / 0.2e1 + t3 * t811 / 0.2e1 + (-t106 * t479 + t478 * t648) * t900 + ((t521 * t648 - t34) * t479 + (-qJD(1) * t649 + t106 * t521 + t16 * t528 + t17 * t531) * t478) * t898 + t478 * t29 * t697 + t4 * t812 / 0.2e1 + (-t105 * t479 + t478 * t650) * t901 + ((t521 * t650 - t35) * t479 + (-qJD(1) * t651 + t105 * t521 + t18 * t528 + t19 * t531) * t478) * t897 + t33 * t701 - t479 * (t859 + t861 + t862 + t869 - t870) / 0.2e1 + (-t112 * t479 + t478 * t646) * t662 + ((t521 * t646 - t41) * t479 + (-qJD(1) * t647 + t112 * t521 + t20 * t528 + t21 * t531) * t478) * t892 + (t378 * t903 + t559 * t379 - t531 * t570) * t899 + (-t376 * t903 + t377 * t559 - t528 * t570) * t896 + (t577 * t479 + (-t526 * t903 + t529 * t559) * t478) * t893 + t925 * t479 * t521 / 0.2e1 + ((-t66 * t103 + t65 * t104 + t23 * t195 - t22 * t197 + (t48 * t630 + (t528 * t65 - t531 * t66) * t282) * t521) * t479 + (t65 * (t154 * t528 - t195 * t521) + t66 * (-t154 * t531 + t197 * t521) + t13 * t630 + t48 * (-t103 * t528 + t104 * t531 - t195 * t742 - t197 * t741) + (qJD(1) * t652 + t867 - t868) * t282) * t478 - t65 * (-t239 * t440 - t345 * t353) - t66 * (t240 * t440 - t345 * t352) - t48 * (t239 * t352 + t240 * t353)) * m(6);];
tauc = t1(:);
