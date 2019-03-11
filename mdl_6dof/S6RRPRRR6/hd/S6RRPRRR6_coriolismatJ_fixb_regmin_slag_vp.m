% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x35]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:30
% EndTime: 2019-03-09 13:54:06
% DurationCPUTime: 23.41s
% Computational Cost: add. (11858->683), mult. (22580->861), div. (0->0), fcn. (25945->8), ass. (0->532)
t471 = sin(qJ(5));
t472 = sin(qJ(2));
t844 = pkin(7) - pkin(8);
t437 = t844 * t472;
t474 = cos(qJ(2));
t438 = t844 * t474;
t800 = sin(qJ(4));
t802 = cos(qJ(4));
t276 = t437 * t800 + t438 * t802;
t419 = -t472 * t800 - t474 * t802;
t488 = t419 * pkin(9) + t276;
t842 = t471 * t488;
t420 = t472 * t802 - t474 * t800;
t422 = t800 * t438;
t618 = t802 * t437;
t831 = t422 - t618;
t541 = -pkin(9) * t420 - t831;
t801 = cos(qJ(5));
t884 = t801 * t541;
t114 = t884 - t842;
t470 = sin(qJ(6));
t963 = t470 * t114;
t105 = -t963 / 0.2e1;
t106 = t963 / 0.2e1;
t958 = t105 + t106;
t473 = cos(qJ(6));
t962 = t473 * t114;
t107 = -t962 / 0.2e1;
t108 = t962 / 0.2e1;
t957 = t107 + t108;
t837 = t801 * t488;
t885 = t471 * t541;
t855 = t837 + t885;
t960 = t855 * t470;
t911 = t960 / 0.2e1;
t968 = -t960 / 0.2e1;
t971 = t911 + t968;
t959 = t855 * t473;
t914 = -t959 / 0.2e1;
t969 = t959 / 0.2e1;
t970 = t914 + t969;
t464 = qJD(6) * t473;
t448 = t470 * t464;
t396 = t801 * t419;
t735 = t471 * t420;
t835 = -t396 + t735;
t861 = t473 * t835;
t894 = -t861 / 0.2e1;
t907 = t470 * t894;
t893 = t861 / 0.2e1;
t908 = t470 * t893;
t933 = t908 + t907;
t950 = qJD(1) * t933;
t967 = t448 - t950;
t966 = t448 + t950;
t628 = qJD(4) + qJD(5);
t732 = t472 * qJ(3);
t823 = pkin(2) + pkin(3);
t411 = t474 * t823 + pkin(1) + t732;
t362 = -t419 * pkin(4) + t411;
t397 = t801 * t420;
t736 = t471 * t419;
t834 = t397 + t736;
t538 = pkin(5) * t835 - pkin(10) * t834;
t480 = t362 + t538;
t43 = -t473 * t480 + t960;
t961 = t834 * t855;
t965 = -t43 * t834 + t470 * t961;
t44 = t470 * t480 + t959;
t964 = -t44 * t834 + t473 * t961;
t915 = 0.2e1 * t893;
t922 = t915 * qJD(1);
t956 = t464 + t922;
t862 = t473 * t834;
t892 = -t862 / 0.2e1;
t916 = 0.2e1 * t892;
t863 = t470 * t835;
t889 = t863 / 0.2e1;
t890 = -t863 / 0.2e1;
t901 = t890 + t889;
t927 = qJD(6) * t901;
t348 = t834 ^ 2;
t881 = -t835 ^ 2 + t348;
t905 = t881 * t470;
t931 = qJD(1) * t905;
t955 = qJD(2) * t916 + t927 - t931;
t864 = t470 * t834;
t888 = -t864 / 0.2e1;
t918 = 0.2e1 * t888;
t899 = t893 + t894;
t928 = qJD(6) * t899;
t904 = t881 * t473;
t932 = qJD(1) * t904;
t954 = qJD(2) * t918 + t928 + t932;
t670 = qJD(6) * t470;
t917 = 0.2e1 * t890;
t921 = t917 * qJD(1);
t953 = t921 - t670;
t952 = qJD(6) * t917 + t931;
t951 = qJD(6) * t915 - t932;
t949 = qJD(2) * t933;
t948 = qJD(4) * t933;
t947 = qJD(5) * t933;
t945 = qJD(2) - t628;
t944 = qJD(2) + qJD(4);
t942 = qJD(2) + t628;
t467 = t470 ^ 2;
t469 = t473 ^ 2;
t569 = t467 / 0.2e1 - t469 / 0.2e1;
t138 = t569 * t834;
t668 = t138 * qJD(6);
t938 = 0.2e1 * t908;
t941 = t938 * qJD(2) - t668;
t940 = 0.2e1 * t105;
t939 = 0.2e1 * t107;
t937 = 0.2e1 * t968;
t936 = 0.2e1 * t969;
t546 = -t884 / 0.2e1;
t851 = t842 / 0.2e1 + t546;
t871 = t884 / 0.2e1;
t935 = -t842 / 0.2e1 + t871;
t891 = t862 / 0.2e1;
t900 = t891 + t892;
t930 = qJD(5) * t900;
t887 = t864 / 0.2e1;
t902 = t887 + t888;
t929 = qJD(5) * t902;
t924 = t900 * qJD(1);
t923 = t902 * qJD(1);
t896 = t885 / 0.2e1 + t837 / 0.2e1;
t421 = t471 * t802 + t800 * t801;
t910 = t421 * t887;
t535 = t471 * t800 - t801 * t802;
t873 = -t835 / 0.2e1;
t502 = t535 * t873;
t872 = t835 / 0.2e1;
t503 = t535 * t872;
t903 = -t502 - t503;
t880 = 0.2e1 * t896;
t898 = qJD(2) * t880;
t897 = t881 * qJD(1);
t549 = t834 * t448;
t262 = -0.2e1 * t549;
t703 = t467 * t873 + t469 * t872;
t845 = t467 * t872 + t469 * t873;
t879 = -t845 + t703;
t895 = t879 * qJD(2) + t262;
t794 = t471 * pkin(4);
t454 = pkin(10) + t794;
t886 = t454 * t834;
t796 = t834 * pkin(5);
t797 = t835 * pkin(10);
t846 = -t797 - t796;
t874 = t834 / 0.2e1;
t875 = -t834 / 0.2e1;
t573 = t874 + t875;
t630 = qJD(2) - qJD(4);
t878 = t421 * t573;
t693 = qJD(1) * t834;
t856 = t835 * t693;
t877 = t469 * t856 - t668;
t866 = t362 * t834;
t865 = t362 * t835;
t742 = t470 * t473;
t608 = qJD(1) * t742;
t547 = t834 * t608;
t536 = 0.2e1 * t547;
t860 = t536 * t835;
t694 = qJD(1) * t835;
t859 = t694 * t834;
t810 = t421 / 0.2e1;
t858 = t810 * t834;
t857 = -t886 / 0.2e1;
t636 = t421 * qJD(2);
t839 = t628 * t421;
t225 = t636 - t839;
t840 = t945 * t742;
t86 = t138 * qJD(1) + t840;
t384 = t397 / 0.2e1;
t585 = t736 / 0.2e1;
t830 = t585 + t384;
t854 = qJD(6) * t830;
t264 = t419 ^ 2 - t420 ^ 2;
t853 = t264 * qJD(1);
t852 = t830 * qJD(1);
t805 = t470 / 0.2e1;
t501 = t535 * t805;
t345 = 0.2e1 * t501;
t848 = t345 * qJD(6) - t473 * t839;
t828 = t535 * qJD(2);
t224 = t535 * t628 - t828;
t444 = t469 - t467;
t165 = t444 * t348;
t847 = -t165 * qJD(1) + 0.2e1 * t834 * t840;
t809 = t618 / 0.2e1;
t843 = t444 * t945;
t803 = t473 / 0.2e1;
t500 = t535 * t803;
t426 = qJ(3) * t802 - t800 * t823;
t412 = t801 * t426;
t456 = t802 * t823;
t612 = t800 * qJ(3);
t425 = t612 + t456;
t520 = -pkin(4) - t425;
t505 = t471 * t520;
t359 = t412 + t505;
t629 = qJD(2) - qJD(5);
t556 = t629 * t359;
t838 = t630 * t419;
t558 = t630 * t420;
t376 = -t736 / 0.2e1;
t811 = -t397 / 0.2e1;
t337 = t376 + t811;
t734 = t471 * t425;
t360 = t412 - t734;
t409 = t421 * qJD(3);
t833 = -qJD(2) * t360 - t409;
t832 = -qJD(2) * t359 - t409;
t583 = t735 / 0.2e1;
t584 = -t735 / 0.2e1;
t340 = t584 + t583;
t829 = t630 * t276;
t405 = t535 * qJD(3);
t347 = 0.2e1 * t500;
t672 = qJD(5) * t470;
t677 = qJD(4) * t470;
t827 = t347 * qJD(6) + (t672 + t677) * t421;
t826 = t138 * t945 - t348 * t608;
t825 = t878 + t903;
t404 = t801 * t520;
t733 = t471 * t426;
t358 = -t404 + t733;
t356 = pkin(5) + t358;
t814 = -t356 / 0.2e1;
t357 = -pkin(10) + t359;
t813 = t357 / 0.2e1;
t812 = -t396 / 0.2e1;
t626 = t801 * pkin(4);
t455 = -t626 - pkin(5);
t807 = -t455 / 0.2e1;
t806 = -t470 / 0.2e1;
t804 = -t473 / 0.2e1;
t799 = pkin(5) * t470;
t798 = pkin(5) * t473;
t795 = t420 * pkin(4);
t39 = t825 * t473;
t793 = t39 * qJD(4);
t36 = t825 * t470;
t792 = t36 * qJD(4);
t465 = t474 * qJ(3);
t417 = -t472 * t823 + t465;
t363 = t417 - t795;
t88 = t363 + t846;
t1 = t861 * t88 - t965;
t791 = t1 * qJD(1);
t2 = -t863 * t88 - t964;
t790 = t2 * qJD(1);
t5 = -t846 * t861 + (-t43 + t960) * t834;
t789 = t5 * qJD(1);
t6 = t846 * t863 + (-t44 + t959) * t834;
t788 = t6 * qJD(1);
t613 = t801 * t425;
t361 = -t613 - t733;
t512 = t455 * t872 - t857;
t476 = t360 * t875 + t361 * t872 + t813 * t834 - t814 * t835 + t512;
t7 = t470 * t476 + t970;
t787 = t7 * qJD(1);
t9 = t473 * t476 + t971;
t786 = t9 * qJD(1);
t27 = t114 * t864 + t43 * t835;
t785 = qJD(1) * t27;
t28 = -t114 * t862 - t44 * t835;
t784 = qJD(1) * t28;
t783 = qJD(1) * t39;
t57 = t569 * t835 - t845;
t782 = qJD(1) * t57;
t61 = t845 + t703;
t780 = qJD(1) * t61;
t69 = t363 * t835 - t866;
t775 = qJD(1) * t69;
t70 = t363 * t834 + t865;
t774 = qJD(1) * t70;
t71 = -t795 * t835 - t866;
t773 = qJD(1) * t71;
t72 = -t795 * t834 + t865;
t772 = qJD(1) * t72;
t571 = t358 / 0.2e1 + t814;
t484 = (-t357 / 0.2e1 + t359 / 0.2e1) * t834 + t571 * t835;
t477 = pkin(5) * t872 + pkin(10) * t875 + t484;
t11 = t477 * t470 + t970;
t769 = t11 * qJD(1);
t14 = t477 * t473 + t971;
t762 = t14 * qJD(1);
t761 = t36 * qJD(1);
t760 = t360 * t470;
t759 = t411 * t419;
t757 = t455 * t470;
t756 = t455 * t473;
t51 = t935 + t851;
t721 = t51 * qJD(1);
t53 = t871 + t546;
t720 = t53 * qJD(1);
t485 = t535 * t875 + t810 * t835;
t481 = t472 / 0.2e1 + t485;
t74 = t481 * t470;
t717 = t74 * qJD(1);
t77 = t481 * t473;
t716 = t77 * qJD(1);
t714 = t470 * t503 + t910;
t692 = qJD(1) * t362;
t691 = qJD(1) * t420;
t690 = qJD(1) * t472;
t689 = qJD(1) * t474;
t688 = qJD(2) * qJ(3);
t687 = qJD(2) * t834;
t686 = qJD(2) * t835;
t683 = qJD(2) * t470;
t682 = qJD(2) * t473;
t681 = qJD(3) * t472;
t680 = qJD(4) * t834;
t679 = qJD(4) * t835;
t678 = qJD(4) * t411;
t676 = qJD(4) * t473;
t675 = qJD(5) * t835;
t674 = qJD(5) * t834;
t673 = qJD(5) * t362;
t671 = qJD(5) * t473;
t139 = t573 * t470;
t667 = t139 * qJD(1);
t143 = 0.2e1 * t889;
t664 = t143 * qJD(1);
t151 = t573 * t473;
t661 = t151 * qJD(1);
t160 = 0.2e1 * t894;
t656 = t160 * qJD(1);
t185 = t830 + t337;
t654 = t185 * qJD(1);
t341 = t812 + t396 / 0.2e1;
t187 = t341 + t340;
t653 = t187 * qJD(1);
t199 = -t411 * t420 - t417 * t419;
t652 = t199 * qJD(1);
t200 = t417 * t420 - t759;
t651 = t200 * qJD(1);
t646 = t337 * qJD(1);
t338 = t585 + t376;
t645 = t338 * qJD(1);
t643 = t340 * qJD(1);
t642 = t341 * qJD(1);
t343 = t811 + t384;
t640 = t343 * qJD(1);
t534 = -pkin(2) * t474 - t732;
t431 = -pkin(1) + t534;
t434 = pkin(2) * t472 - t465;
t364 = t431 * t474 + t434 * t472;
t639 = t364 * qJD(1);
t365 = -t431 * t472 + t434 * t474;
t638 = t365 * qJD(1);
t367 = t809 - t618 / 0.2e1;
t637 = t367 * qJD(1);
t635 = t444 * qJD(6);
t468 = t472 ^ 2;
t445 = t474 ^ 2 - t468;
t634 = t445 * qJD(1);
t633 = t468 * qJD(1);
t632 = t472 * qJD(2);
t463 = t474 * qJD(2);
t631 = t474 * qJD(3);
t625 = pkin(1) * t690;
t624 = pkin(1) * t689;
t623 = qJD(4) * t794;
t622 = pkin(7) * t463;
t621 = qJD(5) * t794;
t620 = pkin(7) * t632;
t619 = t794 / 0.2e1;
t611 = t835 * t692;
t610 = t834 * t692;
t609 = t469 * t693;
t607 = t835 * t681;
t606 = qJD(6) * t835 * t834;
t602 = qJD(1) * t759;
t601 = t411 * t691;
t600 = t419 * t691;
t599 = t431 * t434 * qJD(1);
t598 = t431 * t690;
t597 = t835 * t690;
t596 = t834 * t690;
t595 = t419 * t690;
t594 = t420 * t690;
t586 = t361 * t806;
t574 = t361 * t804;
t570 = -t361 / 0.2e1 + t814;
t568 = t802 * qJD(2);
t567 = t802 * qJD(3);
t566 = t801 * qJD(4);
t565 = t801 * qJD(5);
t564 = t800 * qJD(2);
t563 = t800 * qJD(3);
t37 = (t502 - t858) * t470 + t714;
t393 = t473 * t636;
t560 = t37 * qJD(1) + t393;
t392 = t470 * t636;
t42 = (t878 - t903) * t473;
t559 = qJD(1) * t42 - t392;
t557 = t630 * t360;
t553 = -t626 / 0.2e1;
t552 = -t392 + t827;
t551 = t470 * t597;
t550 = t473 * t597;
t542 = t613 / 0.2e1;
t540 = t619 - t359 / 0.2e1;
t539 = t628 * t794;
t537 = -0.2e1 * t547;
t533 = t470 * t553;
t532 = t473 * t553;
t117 = t795 - t846;
t3 = t117 * t861 + t965;
t531 = t3 * qJD(1) - t36 * qJD(3);
t4 = -t117 * t863 + t964;
t530 = t4 * qJD(1) - t39 * qJD(3);
t528 = t356 * t835 + t357 * t834;
t526 = -t455 * t835 - t886;
t525 = 0.2e1 * t812;
t514 = t360 / 0.2e1 + t540;
t130 = t514 * t473;
t478 = t857 + t835 * t807 + (t471 * t874 + t801 * t873) * pkin(4);
t475 = pkin(5) * t873 + pkin(10) * t874 + t478;
t15 = t470 * t475 + t970;
t524 = t15 * qJD(1) + t130 * qJD(2);
t189 = (t425 / 0.2e1 - pkin(4) - t612 / 0.2e1 - t456 / 0.2e1) * t471;
t523 = t189 * qJD(2);
t518 = t404 / 0.2e1 + t553;
t191 = t542 + t518;
t522 = t191 * qJD(2);
t521 = t834 * (-qJD(6) - t694);
t519 = -qJD(4) * t360 - qJD(5) * t359;
t517 = t797 / 0.2e1 + t796 / 0.2e1;
t513 = t813 * t835 + t814 * t834;
t511 = t454 * t872 + t807 * t834;
t497 = t88 / 0.2e1 + t513;
t19 = t497 * t470 + t957;
t510 = -qJD(1) * t19 + t356 * t682;
t21 = -t497 * t473 + t958;
t509 = -qJD(1) * t21 + t356 * t683;
t508 = qJD(6) * t337 - t859;
t507 = t854 + t859;
t506 = t856 + t854;
t504 = -t846 / 0.2e1 + t517;
t498 = t393 + t848;
t494 = t117 / 0.2e1 + t511;
t493 = -t675 - t679 + t686;
t492 = -t519 + t409;
t427 = t757 / 0.2e1;
t118 = t470 * t570 + t427;
t25 = -t494 * t473 + t958;
t491 = -qJD(1) * t25 + qJD(2) * t118 - t455 * t677;
t429 = t756 / 0.2e1;
t119 = t473 * t570 + t429;
t23 = t494 * t470 + t957;
t490 = -qJD(1) * t23 + qJD(2) * t119 - t455 * t676;
t489 = t493 * t834;
t129 = t514 * t470;
t18 = t473 * t475 + t971;
t487 = -qJD(1) * t18 + qJD(2) * t129 - t470 * t623;
t486 = qJD(2) * t534 + t631;
t460 = -t799 / 0.2e1;
t134 = t470 * t571 + t460;
t31 = -t473 * t504 + t958;
t428 = -t757 / 0.2e1;
t459 = t799 / 0.2e1;
t368 = t533 + t459 + t428;
t483 = pkin(5) * t672 - qJD(1) * t31 + qJD(2) * t134 + qJD(4) * t368;
t462 = -t798 / 0.2e1;
t135 = t473 * t571 + t462;
t29 = t470 * t504 + t957;
t430 = -t756 / 0.2e1;
t461 = t798 / 0.2e1;
t369 = t532 + t461 + t430;
t482 = pkin(5) * t671 - qJD(1) * t29 + qJD(2) * t135 + qJD(4) * t369;
t449 = t472 * t689;
t443 = t470 * t621;
t433 = -qJD(4) * t802 + t568;
t432 = -qJD(4) * t800 + t564;
t391 = t470 * t409;
t371 = t462 + t429 + t532;
t370 = t460 + t427 + t533;
t328 = t356 * t803;
t327 = t356 * t805;
t275 = -t422 + 0.2e1 * t809;
t261 = 0.2e1 * t549;
t204 = t397 + 0.2e1 * t585;
t203 = 0.2e1 * t811 - t736;
t202 = t525 + t735;
t201 = 0.2e1 * t584 + t396;
t192 = -t518 + t542 + t733;
t190 = t619 - t412 - t505 / 0.2e1 + t734 / 0.2e1;
t188 = 0.2e1 * t583 + t525;
t186 = 0.2e1 * t830;
t167 = t537 - t843;
t166 = t536 + t843;
t161 = 0.2e1 * t891;
t149 = 0.2e1 * t887;
t137 = t358 * t803 + t328 + t461;
t136 = t358 * t805 + t327 + t459;
t131 = t360 * t804 + t473 * t540;
t128 = t760 / 0.2e1 - t540 * t470;
t122 = 0.2e1 * t907;
t121 = t430 + t328 + t574;
t120 = t428 + t327 + t586;
t76 = t472 * t803 - t473 * t485;
t75 = t470 * t485 + t472 * t806;
t73 = qJD(2) * t337 + t628 * t830;
t60 = 0.2e1 * t845;
t56 = -0.2e1 * t896;
t54 = -t851 + t935;
t52 = 0.2e1 * t935;
t50 = 0.2e1 * t851;
t41 = t421 * t891 + t835 * t500 + (t503 + t858) * t473;
t38 = t501 * t835 + t714 + t910;
t32 = -t473 * t517 - t803 * t846 + t940;
t30 = t470 * t517 - t806 * t846 + t939;
t26 = t117 * t803 - t511 * t473 + t940;
t24 = t117 * t806 + t511 * t470 + t939;
t22 = -t513 * t473 + t88 * t803 + 0.2e1 * t106;
t20 = t513 * t470 + t88 * t806 + 0.2e1 * t108;
t17 = pkin(5) * t893 + pkin(10) * t892 + t473 * t478 + 0.2e1 * t911;
t16 = pkin(5) * t889 + pkin(10) * t888 + t470 * t478 + 0.2e1 * t914;
t13 = pkin(5) * t894 + pkin(10) * t891 + t473 * t484 + t937;
t12 = pkin(5) * t890 + pkin(10) * t887 + t470 * t484 + t936;
t10 = t356 * t894 + t357 * t892 + t360 * t891 + t473 * t512 + t574 * t835 + t937;
t8 = t356 * t890 + t357 * t888 + t360 * t887 + t470 * t512 + t586 * t835 + t936;
t33 = [0, 0, 0, t472 * t463, t445 * qJD(2), 0, 0, 0, -pkin(1) * t632, -pkin(1) * t463, -t365 * qJD(2) + t472 * t631, 0, -qJD(2) * t364 + qJD(3) * t468 (qJD(2) * t434 - t681) * t431, -t419 * t558, -t630 * t264, 0, 0, 0, qJD(2) * t199 - t419 * t681 + t420 * t678, qJD(2) * t200 + t419 * t678 + t420 * t681, t489, t945 * t881, 0, 0, 0, qJD(2) * t69 - qJD(4) * t71 + t673 * t834 + t607, qJD(2) * t70 - qJD(4) * t72 - t673 * t835 + t681 * t834, -t348 * t448 + t469 * t489, -0.2e1 * t470 * t493 * t862 - t165 * qJD(6), -t470 * t606 - t904 * t945, -t473 * t606 + t905 * t945 (t674 + t680 - t687) * t835, qJD(2) * t1 + qJD(4) * t3 + qJD(5) * t5 + qJD(6) * t28 + t473 * t607, qJD(2) * t2 + qJD(4) * t4 + qJD(5) * t6 + qJD(6) * t27 - t470 * t607; 0, 0, 0, t449, t634, t463, -t632, 0, -t622 - t625, t620 - t624, -t622 - t638, t486, -t620 - t639, pkin(7) * t486 + t599, -t600, -t853, t838, -t558, 0, t652 - t829, qJD(2) * t831 + t275 * qJD(4) + t651, t856, t897, qJD(4) * t188 + qJD(5) * t202 - t686, qJD(4) * t186 + qJD(5) * t204 - t687, 0, -qJD(2) * t855 + t628 * t880 + t775, -qJD(2) * t114 + qJD(4) * t54 + qJD(5) * t52 + t774, t668 + (-t470 * t682 + t609) * t835 + t628 * t938, t261 + (-qJD(2) * t444 + t537) * t835 + t628 * t879, t628 * t918 + t683 * t834 + t928 - t932, t628 * t916 + t682 * t834 + t927 + t931, t508, t791 + (t470 * t528 - t959) * qJD(2) + t38 * qJD(3) + t8 * qJD(4) + t12 * qJD(5) + t22 * qJD(6), t790 + (t473 * t528 + t960) * qJD(2) + t41 * qJD(3) + t10 * qJD(4) + t13 * qJD(5) + t20 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t449, t463, t633, -t598 + t622, 0, 0, 0, 0, 0, -t595, t594, 0, 0, 0, 0, 0, t597, t596, 0, 0, 0, 0, 0, qJD(2) * t38 + qJD(6) * t76 + t550 - t792, qJD(2) * t41 + qJD(6) * t75 - t551 - t793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t600, t853, -t838, t558, 0, t601 + t829, t275 * qJD(2) + qJD(4) * t831 + t602, -t856, -t897, qJD(2) * t188 + qJD(5) * t201 - t679, qJD(2) * t186 + qJD(5) * t203 - t680, 0, -qJD(4) * t855 + qJD(5) * t56 - t773 + t898, qJD(2) * t54 - qJD(4) * t114 + qJD(5) * t50 - t772, t122 * qJD(5) - (t470 * t676 + t609) * t835 + t941, t60 * qJD(5) - (qJD(4) * t444 + t537) * t835 + t895, qJD(5) * t149 + t677 * t834 + t954, qJD(5) * t161 + t676 * t834 + t955, t507, t8 * qJD(2) + (t470 * t526 - t959) * qJD(4) + t16 * qJD(5) + t26 * qJD(6) + t531, t10 * qJD(2) + (t473 * t526 + t960) * qJD(4) + t17 * qJD(5) + t24 * qJD(6) + t530; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t856, -t897, qJD(2) * t202 + qJD(4) * t201 - t675, qJD(2) * t204 + qJD(4) * t203 - t674, 0, qJD(4) * t56 - qJD(5) * t855 + t610 + t898, qJD(2) * t52 + qJD(4) * t50 - qJD(5) * t114 - t611, t122 * qJD(4) + (-t470 * t671 - t609) * t835 + t941, t60 * qJD(4) + (-qJD(5) * t444 + t536) * t835 + t895, qJD(4) * t149 + t672 * t834 + t954, qJD(4) * t161 + t671 * t834 + t955, t506, t789 + t12 * qJD(2) + t16 * qJD(4) + (t470 * t538 - t959) * qJD(5) + t32 * qJD(6), t788 + t13 * qJD(2) + t17 * qJD(4) + (t473 * t538 + t960) * qJD(5) + t30 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t826, t847, t470 * t521 + t899 * t942, t473 * t521 + t901 * t942, t73, qJD(2) * t22 + qJD(3) * t76 + qJD(4) * t26 + qJD(5) * t32 - qJD(6) * t44 + t784, qJD(2) * t20 + qJD(3) * t75 + qJD(4) * t24 + qJD(5) * t30 + qJD(6) * t43 + t785; 0, 0, 0, -t449, -t634, 0, 0, 0, t625, t624, t638, 0, t639, -t599, t600, t853, 0, 0, 0, -t652, qJD(4) * t367 - t651, -t856, -t897, -qJD(4) * t187 - qJD(5) * t341, -qJD(4) * t185 - qJD(5) * t338, 0, -t775, qJD(4) * t53 + qJD(5) * t51 - t774, -t877 + t947 - t948, -qJD(4) * t61 + qJD(5) * t57 + t261 + t860, -qJD(4) * t139 + qJD(6) * t160 - t929 + t932, -qJD(4) * t151 + qJD(6) * t143 - t930 - t931, -t508, qJD(3) * t37 - qJD(4) * t7 + qJD(5) * t11 + qJD(6) * t21 - t791, qJD(3) * t42 - qJD(4) * t9 + qJD(5) * t14 + qJD(6) * t19 - t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, 0, 0, 0, 0, t426 * qJD(4) + t563, -t425 * qJD(4) + t567, 0, 0, 0, 0, 0, t492, t361 * qJD(4) - t358 * qJD(5) - t405, t448, t635, 0, 0, 0, -t356 * t670 + t473 * t492, -t356 * t464 + t470 * t519 - t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t688, 0, 0, 0, 0, 0, t564, t568, 0, 0, 0, 0, 0, t636, -t828, 0, 0, 0, 0, 0, t560, t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t630 * t426, -t425 * t630 + t637, 0, 0, -t653, -t654, 0, t190 * qJD(5) + t557, t192 * qJD(5) + t361 * t630 + t720, -t966, -t635 - t780, -t667, -t661, 0, t131 * qJD(5) + t120 * qJD(6) + t473 * t557 - t787, t128 * qJD(5) + t121 * qJD(6) - t630 * t760 - t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t642, -t645, 0, t190 * qJD(4) + t556, t192 * qJD(4) - t358 * t629 + t721, -t967, -t635 + t782, -t923, -t924, 0, t131 * qJD(4) + t136 * qJD(6) + t473 * t556 + t769, t128 * qJD(4) + t137 * qJD(6) - t470 * t556 + t762; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t166, -t464 + t656, t664 + t670, -t646, qJD(4) * t120 + qJD(5) * t136 - t357 * t464 - t509, qJD(4) * t121 + qJD(5) * t137 + t357 * t670 - t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t449, 0, -t633, t598, 0, 0, 0, 0, 0, t595, -t594, 0, 0, 0, 0, 0, -t597, -t596, 0, 0, 0, 0, 0, -qJD(2) * t37 - qJD(6) * t77 - t550 - t792, -qJD(2) * t42 + qJD(6) * t74 + t551 - t793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t688, 0, 0, 0, 0, 0, -t432, -t433, 0, 0, 0, 0, 0, -t225, -t224, 0, 0, 0, 0, 0, -t560 - t848, -t559 - t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t432, t433, 0, 0, 0, 0, 0, t225, t224, 0, 0, 0, 0, 0, t498 - t761, t552 - t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, t224, 0, 0, 0, 0, 0, t498, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345 * t945 - t421 * t464 - t716, -t347 * t945 + t421 * t670 + t717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t600, -t853, 0, 0, 0, -t601, -qJD(2) * t367 - t602, t856, t897, qJD(2) * t187 + qJD(5) * t340, qJD(2) * t185 + qJD(5) * t343, 0, t773, -qJD(2) * t53 + t772, t877 + t947 + t949, qJD(2) * t61 + t262 - t860, qJD(2) * t139 - t929 + t951, qJD(2) * t151 - t930 + t952, -t507, qJD(2) * t7 + qJD(5) * t15 + qJD(6) * t25 - t531, qJD(2) * t9 + qJD(5) * t18 + qJD(6) * t23 - t530; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426 * qJD(2) - t563, t425 * qJD(2) - t567 - t637, 0, 0, t653, t654, 0, -qJD(5) * t189 + t833, -qJD(2) * t361 - qJD(5) * t191 + t405 - t720, -t967, -t635 + t780, t667, t661, 0, t130 * qJD(5) - t118 * qJD(6) + t473 * t833 + t787, -qJD(5) * t129 - qJD(6) * t119 + t360 * t683 + t391 + t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t564, -t568, 0, 0, 0, 0, 0, -t636, t828, 0, 0, 0, 0, 0, -t393 + t761, t392 + t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t621, -pkin(4) * t565, t448, t635, 0, 0, 0, t455 * t670 - t473 * t621, t455 * t464 + t443; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t643, t640, 0, -t539 - t523 (-t566 - t565) * pkin(4) - t522, t966, t635, -t923, -t924, 0, t370 * qJD(6) - t473 * t539 + t524, qJD(6) * t371 + t443 - t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t167, t956, t953, -t852, qJD(5) * t370 - t454 * t464 - t491, qJD(5) * t371 + t454 * t670 - t490; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t856, t897, qJD(2) * t341 - qJD(4) * t340, qJD(2) * t338 - qJD(4) * t343, 0, -t610, -qJD(2) * t51 + t611, t877 - t948 - t949, -qJD(2) * t57 + t537 * t835 + t262, t902 * t944 + t951, t900 * t944 + t952, -t506, -qJD(2) * t11 - qJD(4) * t15 + qJD(6) * t31 - t789, -qJD(2) * t14 - qJD(4) * t18 + qJD(6) * t29 - t788; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t642, t645, 0, qJD(4) * t189 + t832, qJD(2) * t358 + qJD(4) * t191 + t405 - t721, -t966, -t635 - t782, t923, t924, 0, -t130 * qJD(4) - t134 * qJD(6) + t473 * t832 - t769, qJD(4) * t129 - qJD(6) * t135 + t359 * t683 + t391 - t762; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t636, t828, 0, 0, 0, 0, 0, -t393, t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t643, -t640, 0, t523 + t623, pkin(4) * t566 + t522, t967, t635, t923, t924, 0, -qJD(6) * t368 + t473 * t623 - t524, -qJD(6) * t369 + t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t448, t635, 0, 0, 0, -pkin(5) * t670, -pkin(5) * t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t167, t956, t953, -t852, -pkin(10) * t464 - t483, pkin(10) * t670 - t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t826, -t847, -qJD(2) * t160 + t470 * t856 - t628 * t915, -qJD(2) * t143 + t473 * t856 - t628 * t917, t73, -qJD(2) * t21 + qJD(3) * t77 - qJD(4) * t25 - qJD(5) * t31 - t784, -qJD(2) * t19 - qJD(3) * t74 - qJD(4) * t23 - qJD(5) * t29 - t785; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t167, -t656, -t664, t646, qJD(3) * t345 + qJD(4) * t118 + qJD(5) * t134 + t509, qJD(3) * t347 + qJD(4) * t119 + qJD(5) * t135 + t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345 * qJD(2) + t716, t347 * qJD(2) - t717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t166, -t922, -t921, t852, qJD(5) * t368 + t491, qJD(5) * t369 + t490; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t166, -t922, -t921, t852, t483, t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t33;
