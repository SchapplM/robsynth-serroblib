% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x35]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:55
% EndTime: 2019-03-09 07:10:33
% DurationCPUTime: 22.64s
% Computational Cost: add. (20243->702), mult. (40916->899), div. (0->0), fcn. (50796->10), ass. (0->557)
t645 = qJD(3) + qJD(4);
t644 = qJD(5) + qJD(6);
t499 = sin(pkin(11));
t500 = cos(pkin(11));
t504 = sin(qJ(3));
t507 = cos(qJ(3));
t460 = -t504 * t499 + t500 * t507;
t461 = t499 * t507 + t504 * t500;
t503 = sin(qJ(4));
t832 = cos(qJ(4));
t414 = -t832 * t460 + t503 * t461;
t502 = sin(qJ(5));
t505 = cos(qJ(6));
t730 = t505 * t502;
t501 = sin(qJ(6));
t506 = cos(qJ(5));
t742 = t501 * t506;
t464 = t730 + t742;
t887 = t464 * t414;
t312 = -t887 / 0.2e1;
t545 = t742 / 0.2e1 + t730 / 0.2e1;
t529 = t414 * t545;
t907 = t312 + t529;
t721 = t644 * t907;
t447 = t832 * t461;
t733 = t503 * t460;
t858 = t447 + t733;
t883 = t858 * t464;
t772 = t883 * t858;
t773 = t887 * t414;
t921 = -t772 + t773;
t942 = qJD(1) * t921;
t950 = t721 + t942;
t885 = t506 * t414;
t913 = t885 / 0.2e1;
t914 = -t885 / 0.2e1;
t920 = t914 + t913;
t938 = qJD(5) * t920;
t761 = t414 ^ 2;
t896 = t858 ^ 2;
t906 = t896 - t761;
t925 = t906 * t506;
t940 = qJD(1) * t925;
t949 = t938 + t940;
t948 = 0.2e1 * t502;
t895 = -t414 / 0.2e1;
t934 = -t883 / 0.2e1;
t820 = pkin(7) + qJ(2);
t472 = t820 * t499;
t473 = t820 * t500;
t558 = t504 * t472 - t473 * t507;
t332 = t460 * pkin(8) - t558;
t588 = t507 * t472 + t473 * t504;
t551 = -pkin(8) * t461 - t588;
t881 = t832 * t332 + t503 * t551;
t886 = t502 * t414;
t904 = -pkin(5) * t886 + t881;
t927 = t904 * t464;
t947 = t927 / 0.2e1;
t729 = t505 * t506;
t743 = t501 * t502;
t462 = -t729 + t743;
t928 = t904 * t462;
t946 = t928 / 0.2e1;
t734 = t503 * t332;
t861 = t832 * t551;
t209 = t861 - t734;
t931 = t506 * t209;
t617 = -t931 / 0.2e1;
t726 = t506 * t858;
t882 = t858 * t502;
t233 = t501 * t882 - t505 * t726;
t71 = t464 * t233 + t462 * t883;
t818 = t644 * t71;
t888 = t462 * t414;
t899 = -t233 * t887 - t883 * t888;
t919 = t899 * qJD(1);
t945 = t818 - t919;
t576 = t447 / 0.2e1;
t857 = t576 + t733 / 0.2e1;
t874 = t858 * qJD(1);
t279 = t414 * t874;
t924 = t857 * qJD(5) + t279;
t944 = qJD(6) * t857 + t924;
t891 = pkin(5) * t858;
t932 = t502 * t209;
t943 = pkin(10) * t885 + t891 - t932;
t770 = t888 * t414;
t774 = t233 * t858;
t922 = t770 - t774;
t941 = qJD(1) * t922;
t926 = t906 * t502;
t939 = qJD(1) * t926;
t390 = t882 / 0.2e1;
t614 = -t882 / 0.2e1;
t923 = t614 + t390;
t937 = t923 * qJD(2);
t936 = t818 + t919 + t645 * (-t888 * t462 + t464 * t887);
t935 = -t233 / 0.2e1;
t844 = t414 / 0.2e1;
t933 = t209 * t414;
t918 = t906 * qJD(1);
t821 = t506 * pkin(5);
t491 = -pkin(4) - t821;
t848 = pkin(9) + pkin(10);
t476 = t848 * t502;
t477 = t848 * t506;
t557 = -t505 * t476 - t501 * t477;
t905 = t491 * t934 + t557 * t895;
t149 = pkin(5) * t882 - t209;
t827 = pkin(5) * t414;
t485 = -pkin(2) * t500 - pkin(1);
t430 = -pkin(3) * t460 + t485;
t824 = pkin(9) * t858;
t571 = pkin(4) * t414 - t824;
t520 = t430 + t571;
t910 = t881 * t502;
t105 = -t506 * t520 + t910;
t97 = -pkin(10) * t726 - t105;
t75 = t97 + t827;
t909 = t881 * t506;
t106 = t502 * t520 + t909;
t98 = -pkin(10) * t882 + t106;
t802 = t505 * t98;
t44 = t501 * t75 + t802;
t917 = t149 * t888 - t904 * t233 - t44 * t858;
t74 = t505 * t75;
t806 = t501 * t98;
t43 = -t74 + t806;
t916 = -t149 * t887 - t43 * t858 + t904 * t883;
t915 = t645 * t881;
t615 = t883 / 0.2e1;
t304 = -t888 / 0.2e1;
t637 = t891 / 0.2e1;
t912 = pkin(10) * t886;
t823 = pkin(9) * t414;
t829 = pkin(4) * t858;
t271 = t823 + t829;
t269 = t502 * t271;
t720 = -t931 - t269;
t104 = -t720 + t912;
t731 = t505 * t104;
t604 = -t731 / 0.2e1;
t270 = t506 * t271;
t83 = t270 + t943;
t808 = t501 * t83;
t547 = -t808 / 0.2e1 + t604;
t781 = t149 * t462;
t619 = t781 / 0.2e1;
t26 = t619 + t547 - t905;
t423 = -t501 * t476 + t505 * t477;
t871 = t423 * t895 + t491 * t935;
t641 = t832 * pkin(3);
t490 = -t641 - pkin(4);
t474 = t490 - t821;
t822 = t503 * pkin(3);
t489 = pkin(9) + t822;
t819 = pkin(10) + t489;
t453 = t819 * t502;
t454 = t819 * t506;
t589 = t505 * t453 + t501 * t454;
t873 = t474 * t934 + t589 * t844;
t409 = -t501 * t453 + t505 * t454;
t872 = t409 * t895 + t474 * t935;
t700 = qJD(1) * t233;
t840 = t464 / 0.2e1;
t841 = t462 / 0.2e1;
t111 = t233 * t841 - t840 * t883;
t723 = t644 * t111;
t900 = t700 * t888 + t723;
t833 = -t506 / 0.2e1;
t835 = t502 / 0.2e1;
t898 = (-t883 * t835 + (t462 * t833 + t505 / 0.2e1) * t858) * pkin(5);
t897 = (t233 * t835 + (t464 * t833 - t501 / 0.2e1) * t858) * pkin(5);
t894 = t887 / 0.2e1;
t893 = -t858 / 0.2e1;
t892 = t858 / 0.2e1;
t847 = -t861 / 0.2e1;
t451 = t464 * qJD(6);
t714 = -t464 * qJD(5) - t451;
t880 = t645 * t414;
t497 = t502 ^ 2;
t498 = t506 ^ 2;
t593 = t497 / 0.2e1 - t498 / 0.2e1;
t246 = t593 * t858;
t737 = t502 * t506;
t862 = t645 * t737;
t202 = -qJD(1) * t246 + t862;
t879 = qJD(1) * t414;
t878 = qJD(2) * t414;
t876 = t857 * qJD(1);
t744 = t501 * t104;
t611 = -t744 / 0.2e1;
t804 = t505 * t83;
t549 = t611 + t804 / 0.2e1;
t780 = t149 * t464;
t618 = -t780 / 0.2e1;
t25 = t618 + t549 - t871;
t830 = pkin(3) * t461;
t220 = t830 + t271;
t214 = t502 * t220;
t719 = t931 + t214;
t101 = t719 + t912;
t745 = t501 * t101;
t612 = -t745 / 0.2e1;
t215 = t506 * t220;
t78 = t215 + t943;
t805 = t505 * t78;
t550 = t612 + t805 / 0.2e1;
t21 = t618 + t550 - t872;
t732 = t505 * t101;
t605 = -t732 / 0.2e1;
t809 = t501 * t78;
t548 = -t809 / 0.2e1 + t605;
t22 = t619 + t548 - t873;
t870 = t644 * t409;
t869 = t644 * t423;
t868 = t644 * t557;
t867 = -0.2e1 * t858;
t866 = t414 * t644;
t482 = t498 - t497;
t865 = t482 * t645;
t263 = 0.2e1 * t914;
t306 = t462 * t893;
t178 = t858 * t462;
t616 = -t178 / 0.2e1;
t579 = t306 + t616;
t854 = qJD(4) * t579;
t84 = -t233 ^ 2 + t883 ^ 2;
t19 = qJD(1) * t84 + t645 * t71;
t39 = t111 * t645 + t700 * t883;
t624 = qJD(1) * t737;
t853 = t246 * t645 + t624 * t896;
t355 = t462 ^ 2 - t464 ^ 2;
t59 = qJD(1) * t71 + t355 * t645;
t82 = t462 * t464 * t645 - qJD(1) * t111;
t852 = t644 * t589;
t850 = -t74 / 0.2e1;
t849 = -t75 / 0.2e1;
t613 = t888 / 0.2e1;
t839 = t474 / 0.2e1;
t838 = -t490 / 0.2e1;
t837 = t491 / 0.2e1;
t831 = pkin(3) * t414;
t825 = pkin(5) * t502;
t817 = pkin(5) * qJD(5);
t816 = pkin(5) * qJD(6);
t638 = -t827 / 0.2e1;
t574 = t638 + t97 / 0.2e1;
t5 = (t849 + t574) * t501;
t815 = qJD(1) * t5;
t7 = t505 * t574 + t850;
t814 = qJD(1) * t7;
t1 = (-t745 + t805) * t414 + t916;
t813 = t1 * qJD(1);
t2 = -(t732 + t809) * t414 + t917;
t812 = t2 * qJD(1);
t3 = (-t744 + t804) * t414 + t916;
t811 = t3 * qJD(1);
t4 = -(t731 + t808) * t414 + t917;
t810 = t4 * qJD(1);
t807 = t501 * t97;
t803 = t505 * t97;
t45 = -t802 - t807;
t642 = pkin(5) * t726;
t783 = t149 * t233;
t17 = t414 * t45 + t642 * t883 - t783;
t801 = qJD(1) * t17;
t46 = t803 - t806;
t782 = t149 * t883;
t18 = -t233 * t642 - t414 * t46 - t782;
t800 = qJD(1) * t18;
t29 = t414 * t43 - t782;
t799 = qJD(1) * t29;
t30 = -t414 * t44 - t783;
t798 = qJD(1) * t30;
t511 = t489 * t893 + t414 * t838 + (t503 * t892 + t832 * t895) * pkin(3);
t510 = t824 / 0.2e1 + pkin(4) * t895 + t511;
t47 = t502 * t510;
t797 = qJD(1) * t47;
t597 = t613 + t304;
t54 = (t894 + t312) * t464 - t597 * t462;
t796 = qJD(1) * t54;
t61 = t105 * t414 + t209 * t882;
t795 = qJD(1) * t61;
t62 = -t106 * t414 - t209 * t726;
t794 = qJD(1) * t62;
t93 = t772 + t773;
t789 = qJD(1) * t93;
t96 = -t770 - t774;
t786 = qJD(1) * t96;
t562 = t858 * t881 + t933;
t35 = -t105 * t858 + t215 * t414 + (t562 - t933) * t502;
t768 = t35 * qJD(1);
t36 = -t106 * t858 - t414 * t719 + t506 * t562;
t767 = t36 * qJD(1);
t37 = t270 * t414 + (-t105 + t910) * t858;
t766 = t37 * qJD(1);
t38 = (-t106 + t909) * t858 + (t720 + t931) * t414;
t765 = t38 * qJD(1);
t751 = t474 * t462;
t750 = t474 * t464;
t747 = t491 * t462;
t746 = t491 * t464;
t152 = t894 + t529;
t722 = t644 * t152;
t225 = t644 * t355;
t610 = t743 / 0.2e1;
t718 = -t414 * t610 + t505 * t913;
t607 = t886 / 0.2e1;
t717 = t501 * t607 + t505 * t914;
t595 = 0.2e1 * t893;
t258 = t595 * t506;
t688 = qJD(3) * t506;
t716 = -t258 * qJD(4) + t688 * t858;
t296 = t714 * t462;
t449 = t462 * qJD(6);
t715 = -t462 * qJD(5) - t449;
t481 = t499 ^ 2 + t500 ^ 2;
t110 = t597 * t464;
t713 = qJD(1) * t110;
t561 = t761 + t896;
t130 = t561 * t502;
t708 = qJD(1) * t130;
t132 = t561 * t506;
t706 = qJD(1) * t132;
t163 = -t414 * t830 - t430 * t858;
t703 = qJD(1) * t163;
t164 = t414 * t430 - t830 * t858;
t702 = qJD(1) * t164;
t228 = (t895 + t844) * t737;
t701 = qJD(1) * t228;
t268 = t482 * t896;
t698 = qJD(1) * t268;
t695 = qJD(1) * t430;
t694 = qJD(2) * t858;
t693 = qJD(3) * t858;
t691 = qJD(3) * t464;
t690 = qJD(3) * t474;
t689 = qJD(3) * t502;
t686 = qJD(4) * t858;
t685 = qJD(4) * t430;
t684 = qJD(4) * t464;
t683 = qJD(4) * t491;
t682 = qJD(4) * t502;
t681 = qJD(4) * t506;
t680 = qJD(5) * t502;
t494 = qJD(5) * t506;
t113 = t847 + t861 / 0.2e1;
t679 = t113 * qJD(1);
t151 = (t840 + t545) * t414;
t135 = t151 * qJD(1);
t137 = t152 * qJD(1);
t544 = -t729 / 0.2e1 + t610;
t153 = (t841 + t544) * t414;
t138 = t153 * qJD(1);
t154 = t613 + t717;
t139 = t154 * qJD(1);
t155 = t304 + t718;
t141 = t155 * qJD(1);
t596 = t892 + t893;
t176 = t596 * t462;
t676 = t176 * qJD(1);
t675 = t178 * qJD(1);
t674 = t579 * qJD(1);
t673 = t883 * qJD(1);
t183 = t596 * t464;
t672 = t183 * qJD(1);
t184 = t595 * t464;
t671 = t184 * qJD(1);
t670 = t882 * qJD(1);
t669 = t886 * qJD(1);
t251 = t895 * t948;
t243 = t251 * qJD(1);
t254 = t596 * t502;
t668 = t254 * qJD(1);
t256 = t595 * t502;
t667 = t256 * qJD(1);
t257 = t596 * t506;
t666 = t257 * qJD(1);
t665 = t258 * qJD(1);
t664 = t885 * qJD(1);
t262 = 0.2e1 * t913;
t663 = t262 * qJD(1);
t662 = t263 * qJD(1);
t302 = 0.2e1 * t576 + t733;
t660 = t302 * qJD(1);
t337 = t460 ^ 2 - t461 ^ 2;
t659 = t337 * qJD(1);
t411 = t576 - t447 / 0.2e1;
t654 = t411 * qJD(1);
t653 = t411 * qJD(4);
t650 = t460 * qJD(1);
t448 = t460 * qJD(3);
t649 = t461 * qJD(1);
t648 = t461 * qJD(3);
t469 = t481 * qJ(2);
t647 = t469 * qJD(1);
t646 = t481 * qJD(1);
t640 = qJD(3) * t822;
t639 = qJD(4) * t822;
t636 = -t825 / 0.2e1;
t635 = -t822 / 0.2e1;
t125 = -t781 / 0.2e1;
t585 = t506 * t637;
t634 = t233 * t636 + t464 * t585 + t125;
t126 = t780 / 0.2e1;
t633 = t462 * t585 - t636 * t883 + t126;
t632 = t502 * t832;
t631 = t506 * t832;
t627 = t414 * t695;
t626 = t858 * t695;
t625 = t498 * t874;
t623 = qJD(5) * t414 * t858;
t622 = t858 * t879;
t621 = t460 * t649;
t483 = t502 * t494;
t620 = t506 * t874;
t600 = t931 / 0.2e1 + t269 / 0.2e1;
t599 = t617 - t214 / 0.2e1;
t594 = t837 + t839;
t592 = t832 * qJD(3);
t591 = t832 * qJD(4);
t590 = pkin(5) * t644;
t584 = -qJD(1) * t485 - qJD(2);
t583 = -qJD(5) - t879;
t582 = t501 * t632;
t580 = t858 * t624;
t578 = -t632 / 0.2e1;
t577 = -t631 / 0.2e1;
t575 = t637 + t78 / 0.2e1;
t573 = t637 + t83 / 0.2e1;
t572 = t645 * t822;
t569 = -0.2e1 * t580;
t568 = 0.2e1 * t580;
t565 = qJD(6) - t583;
t445 = t464 * t825;
t354 = t445 - t751;
t9 = t897 + t22;
t563 = -t9 * qJD(1) + t354 * qJD(3);
t559 = -t490 * t414 - t489 * t858;
t10 = t898 + t21;
t444 = t462 * t825;
t353 = t444 + t750;
t556 = -t10 * qJD(1) + t353 * qJD(3);
t555 = t858 * t583;
t554 = qJD(4) * t302 + t693;
t552 = t823 / 0.2e1 + t829 / 0.2e1;
t546 = t489 * t844 + t858 * t838;
t543 = qJD(1) * t21 - t464 * t690;
t542 = qJD(1) * t22 + t462 * t690;
t530 = t546 * t506;
t57 = -t215 / 0.2e1 - t530;
t541 = -qJD(1) * t57 - t490 * t689;
t518 = t502 * t546 + t617;
t55 = t518 - t599;
t540 = -qJD(1) * t55 - t490 * t688;
t536 = -t641 / 0.2e1 + pkin(4) / 0.2e1 + t838;
t535 = t880 * t858;
t509 = t946 + t589 * t893 - t887 * t839 - t883 * t635 + (-t501 * t631 - t505 * t632) * t831 / 0.2e1;
t522 = -t557 * t893 - t887 * t837 + t946;
t32 = t509 - t522;
t534 = -qJD(1) * t32 - t462 * t640;
t508 = t947 + t409 * t893 + t888 * t839 + t233 * t635 - (t505 * t631 - t582) * t831 / 0.2e1;
t521 = t423 * t893 + t888 * t837 + t947;
t34 = t508 - t521;
t533 = -qJD(1) * t34 - t464 * t640;
t50 = t506 * t510;
t532 = -qJD(1) * t50 - t502 * t640;
t531 = t552 * t506;
t13 = t897 + t26;
t516 = (t505 * t577 + t582 / 0.2e1) * pkin(3);
t281 = t462 * t594 + t516;
t230 = -t445 + t281;
t405 = t445 - t747;
t528 = t13 * qJD(1) + t230 * qJD(3) - t405 * qJD(4);
t14 = t898 + t25;
t517 = (t501 * t577 + t505 * t578) * pkin(3);
t280 = -t464 * t594 + t517;
t229 = -t444 + t280;
t404 = t444 + t746;
t527 = t14 * qJD(1) + t229 * qJD(3) - t404 * qJD(4);
t426 = t536 * t502;
t65 = -t270 / 0.2e1 - t531;
t526 = pkin(4) * t682 - qJD(1) * t65 + qJD(3) * t426;
t427 = t536 * t506;
t519 = t502 * t552 + t617;
t63 = t519 + t600;
t525 = pkin(4) * t681 - qJD(1) * t63 + qJD(3) * t427;
t524 = qJD(1) * t25 + qJD(3) * t280 - t464 * t683;
t523 = qJD(1) * t26 + qJD(3) * t281 + t462 * t683;
t282 = -t747 / 0.2e1 - t751 / 0.2e1 + t516;
t283 = t746 / 0.2e1 + t750 / 0.2e1 + t517;
t480 = t502 * t639;
t475 = t482 * qJD(5);
t438 = t464 * t639;
t437 = t462 * t639;
t429 = pkin(4) * t833 + t490 * t506 / 0.2e1 + pkin(3) * t577;
t428 = -pkin(4) * t502 / 0.2e1 + t490 * t835 + pkin(3) * t578;
t338 = t483 * t867;
t307 = t178 / 0.2e1;
t291 = t569 + t865;
t290 = t568 - t865;
t267 = t645 * t857;
t253 = 0.2e1 * t390;
t252 = -t886 / 0.2e1 + t607;
t245 = t257 * qJD(4);
t242 = t252 * qJD(5);
t241 = t251 * qJD(5);
t240 = t246 * qJD(5);
t232 = t445 + t282;
t231 = t444 + t283;
t227 = t263 * t502;
t226 = t243 - t680;
t187 = t615 + t934;
t186 = 0.2e1 * t615;
t180 = t306 + t307;
t177 = t616 + t307;
t158 = t613 + t718;
t157 = t304 + t717;
t156 = t414 * t544 + t304;
t147 = t907 * qJD(2);
t142 = t158 * qJD(2);
t140 = t155 * qJD(2);
t136 = t152 * qJD(2);
t133 = (t686 + t693) * t414;
t122 = t593 * t414 + t497 * t844 + t498 * t895;
t120 = -t137 + t714;
t119 = -t135 + t714;
t118 = t462 * t644 - t141;
t117 = -t139 + t715;
t116 = -t138 + t715;
t114 = t734 + 0.2e1 * t847;
t109 = 0.2e1 * t888 * t840;
t66 = -t932 + t270 / 0.2e1 - t531;
t64 = t519 - t600;
t58 = -t932 + t215 / 0.2e1 - t530;
t56 = t518 + t599;
t52 = qJD(3) * t154 + qJD(4) * t153 + t879 * t883;
t51 = qJD(3) * t152 + qJD(4) * t151 - t233 * t879;
t49 = pkin(4) * t913 + t511 * t506 + t824 * t833 + t910;
t48 = pkin(4) * t607 + pkin(9) * t614 + t511 * t502 - t909;
t42 = qJD(3) * t157 + qJD(4) * t156 - t565 * t883;
t41 = t233 * t565 + t645 * t907;
t33 = t508 + t521;
t31 = t509 + t522;
t28 = t126 + t549 + t871;
t27 = t125 + t547 + t905;
t24 = t126 + t550 + t872;
t23 = t125 + t548 + t873;
t16 = -t573 * t501 + t604 + t634 + t905;
t15 = t573 * t505 + t611 + t633 + t871;
t12 = -t575 * t501 + t605 + t634 + t873;
t11 = t575 * t505 + t612 + t633 + t872;
t8 = t505 * t638 + t806 + t850 - t803 / 0.2e1;
t6 = -t802 - t807 / 0.2e1 + (t638 + t849) * t501;
t20 = [0, 0, 0, 0, 0, t481 * qJD(2), t469 * qJD(2), t460 * t648, t337 * qJD(3), 0, 0, 0, t485 * t648, t485 * t448, -t535, -t645 * t906, 0, 0, 0, -qJD(3) * t163 + t685 * t858, -qJD(3) * t164 - t414 * t685, -t483 * t896 - t498 * t535, t726 * t880 * t948 - qJD(5) * t268, -t502 * t623 + t645 * t925, -t506 * t623 - t645 * t926, t133, qJD(2) * t130 + qJD(3) * t35 + qJD(4) * t37 + qJD(5) * t62, qJD(2) * t132 + qJD(3) * t36 + qJD(4) * t38 + qJD(5) * t61 (t644 * t883 - t645 * t888) * t233, t644 * t84 + t645 * t899, t645 * t922 - t866 * t883, t233 * t866 + t645 * t921, t133, qJD(2) * t93 + qJD(3) * t1 + qJD(4) * t3 + qJD(5) * t17 + qJD(6) * t30, qJD(2) * t96 + qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t18 + qJD(6) * t29; 0, 0, 0, 0, 0, t646, t647, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t653, 0, 0, 0, 0, 0, 0, t242 - t245 + t708, t645 * t923 + t706 + t938, 0, 0, 0, 0, 0, qJD(3) * t177 + qJD(4) * t180 + t721 + t789, qJD(4) * t187 + t158 * t644 + t786; 0, 0, 0, 0, 0, 0, 0, t621, t659, t448, -t648, 0, qJD(3) * t558 + t485 * t649, qJD(3) * t588 + t485 * t650, -t622, -t918, -t880, -t554, 0, -t703 - t915, -qJD(3) * t209 + qJD(4) * t114 - t702, qJD(4) * t227 - t240 - (t502 * t688 + t625) * t414, t122 * qJD(4) + t338 - (qJD(3) * t482 + t569) * t414, qJD(4) * t253 + t689 * t858 + t949, t242 - t939 + t716, t924, t768 + (t502 * t559 - t909) * qJD(3) + t48 * qJD(4) + t58 * qJD(5), t767 + t937 + (t506 * t559 + t910) * qJD(3) + t49 * qJD(4) + t56 * qJD(5), qJD(4) * t109 + (t691 - t700) * t888 + t723, t936, qJD(4) * t186 + t157 * t644 + t691 * t858 + t941, -t462 * t693 + t854 + t950, t944, t813 + t177 * qJD(2) + (-t474 * t887 - t589 * t858 + t928) * qJD(3) + t31 * qJD(4) + t11 * qJD(5) + t24 * qJD(6), t812 + (-t409 * t858 + t474 * t888 + t927) * qJD(3) + t33 * qJD(4) + t12 * qJD(5) + t23 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, -t918, -t880, -qJD(3) * t302 - t686, 0, qJD(2) * t411 + t626 - t915, qJD(3) * t114 - qJD(4) * t209 - t627, qJD(3) * t227 - t240 + (-t502 * t681 - t625) * t414, t122 * qJD(3) + t338 + (-qJD(4) * t482 + t568) * t414, qJD(3) * t253 + t682 * t858 + t949, -qJD(3) * t258 + t681 * t858 - t939, t924, t766 - t257 * qJD(2) + t48 * qJD(3) + (t502 * t571 - t909) * qJD(4) + t66 * qJD(5), t765 + t937 + t49 * qJD(3) + (t506 * t571 + t910) * qJD(4) + t64 * qJD(5), qJD(3) * t109 + (t684 - t700) * t888 + t723, t936, qJD(3) * t186 + t156 * t644 + t684 * t858 + t941, qJD(3) * t579 - t462 * t686 + t950, t944, t811 + t180 * qJD(2) + t31 * qJD(3) + (-t491 * t887 + t557 * t858 + t928) * qJD(4) + t15 * qJD(5) + t28 * qJD(6), t810 + t187 * qJD(2) + t33 * qJD(3) + (-t423 * t858 + t491 * t888 + t927) * qJD(4) + t16 * qJD(5) + t27 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t853, t862 * t867 - t698, t502 * t555 + t645 * t920, qJD(3) * t252 + t506 * t555, t267, qJD(2) * t252 + qJD(3) * t58 + qJD(4) * t66 - qJD(5) * t106 + t794, qJD(2) * t920 + qJD(3) * t56 + qJD(4) * t64 + qJD(5) * t105 + t795, t39, t19, t42, t41, t267, qJD(3) * t11 + qJD(4) * t15 + qJD(5) * t45 + qJD(6) * t6 + t147 + t801, qJD(3) * t12 + qJD(4) * t16 - qJD(5) * t46 + qJD(6) * t8 + t142 + t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t19, t42, t41, t267, qJD(3) * t24 + qJD(4) * t28 + qJD(5) * t6 - qJD(6) * t44 + t147 + t798, qJD(3) * t23 + qJD(4) * t27 + qJD(5) * t8 + qJD(6) * t43 + t142 + t799; 0, 0, 0, 0, 0, -t646, -t647, 0, 0, 0, 0, 0, t648, t448, 0, 0, 0, 0, 0, t554, -t880, 0, 0, 0, 0, 0, t241 - t708 + t716, -qJD(3) * t882 + qJD(4) * t256 + qJD(5) * t263 - t706, 0, 0, 0, 0, 0, -qJD(3) * t178 - t722 - t789 + t854, -qJD(3) * t883 + qJD(4) * t184 - t155 * t644 - t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t649, t650, 0, 0, 0, 0, 0, t874, -t879, 0, 0, 0, 0, 0, t620, -t670, 0, 0, 0, 0, 0, -t675, -t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t660, -t879, 0, 0, 0, 0, 0, -t665, t667, 0, 0, 0, 0, 0, t674, t671; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t494 + t662, 0, 0, 0, 0, 0, t120, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t118; 0, 0, 0, 0, 0, 0, 0, -t621, -t659, 0, 0, 0, t584 * t461, t584 * t460, t622, t918, 0, -t653, 0, -t694 + t703, qJD(4) * t113 + t702 + t878, qJD(4) * t228 + t498 * t622 - t240, -t414 * t568 + t338, -qJD(4) * t254 + qJD(5) * t262 - t940, t241 - t245 + t939, -t924, qJD(4) * t47 + qJD(5) * t57 - t506 * t694 - t768, qJD(2) * t882 + qJD(4) * t50 + qJD(5) * t55 - t767, qJD(4) * t110 + t900, qJD(4) * t54 + t945, -qJD(4) * t183 - t154 * t644 - t941, qJD(4) * t176 - t722 - t942, -t944, qJD(2) * t178 + qJD(4) * t32 - qJD(5) * t10 - qJD(6) * t21 - t813, qJD(2) * t883 + qJD(4) * t34 - qJD(5) * t9 - qJD(6) * t22 - t812; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t649, -t650, 0, 0, 0, 0, 0, -t874, t879, 0, 0, 0, 0, 0, -t620, t670, 0, 0, 0, 0, 0, t675, t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t639, -pkin(3) * t591, t483, t475, 0, 0, 0, t490 * t680 - t506 * t639, t490 * t494 + t480, t296, t225, 0, 0, 0, qJD(5) * t353 + t451 * t474 + t437, qJD(5) * t354 - t449 * t474 + t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t654, 0, -t572, t679 + (-t592 - t591) * pkin(3), t483 + t701, t475, -t668, -t666, 0, qJD(5) * t428 - t506 * t572 + t797, qJD(5) * t429 + t480 - t532, t713 + t296, t225 + t796, -t672, t676, 0, qJD(5) * t231 + qJD(6) * t283 + t437 - t534, qJD(5) * t232 + qJD(6) * t282 + t438 - t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t291, t494 + t663, t226, -t876, qJD(4) * t428 - t489 * t494 - t541, qJD(4) * t429 + t489 * t680 - t540, -t82, t59, t117, t120, -t876, t231 * qJD(4) + t556 - t870, t232 * qJD(4) + t563 + t852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t59, t117, t120, -t876, qJD(4) * t283 - t543 - t870, qJD(4) * t282 - t542 + t852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, t918, 0, t411 * qJD(3), 0, -qJD(2) * t302 - t626, -qJD(3) * t113 + t627 + t878, -qJD(3) * t228 + t279 * t498 - t240, t414 * t569 + t338, qJD(3) * t254 + qJD(5) * t885 - t940, qJD(3) * t257 - qJD(5) * t886 + t939, -t924, qJD(2) * t258 - qJD(3) * t47 + qJD(5) * t65 - t766, -qJD(2) * t256 - qJD(3) * t50 + qJD(5) * t63 - t765, -qJD(3) * t110 + t900, -qJD(3) * t54 + t945, qJD(3) * t183 - t153 * t644 - t941, -qJD(3) * t176 - t151 * t644 - t942, -t944, -qJD(2) * t579 - qJD(3) * t32 - qJD(5) * t14 - qJD(6) * t25 - t811, -qJD(2) * t184 - qJD(3) * t34 - qJD(5) * t13 - qJD(6) * t26 - t810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t660, t879, 0, 0, 0, 0, 0, t665, -t667, 0, 0, 0, 0, 0, -t674, -t671; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t654, 0, t640, pkin(3) * t592 - t679, t483 - t701, t475, t668, t666, 0, -qJD(5) * t426 + t506 * t640 - t797, -qJD(5) * t427 + t532, -t713 + t296, t225 - t796, t672, -t676, 0, -qJD(5) * t229 - qJD(6) * t280 + t534, -qJD(5) * t230 - qJD(6) * t281 + t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, t475, 0, 0, 0, -pkin(4) * t680, -pkin(4) * t494, t296, t225, 0, 0, 0, qJD(5) * t404 + t451 * t491, qJD(5) * t405 - t449 * t491; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t291, t494 + t664, -t669 - t680, -t876, -pkin(9) * t494 - t526, pkin(9) * t680 - t525, -t82, t59, t116, t119, -t876, -t527 - t869, -t528 - t868; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t59, t116, t119, -t876, -t524 - t869, -t523 - t868; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t853, 0.2e1 * t858 * t862 + t698, -qJD(3) * t262 - qJD(4) * t885 + t279 * t502, -qJD(3) * t251 + qJD(4) * t886 + t279 * t506, t267, -qJD(2) * t251 - qJD(3) * t57 - qJD(4) * t65 - t794, -qJD(2) * t263 - qJD(3) * t55 - qJD(4) * t63 - t795, -t39, -t19, t52, t51, t267, qJD(3) * t10 + qJD(4) * t14 + qJD(6) * t5 + t136 - t801, qJD(3) * t9 + qJD(4) * t13 + qJD(6) * t7 + t140 - t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, -t662, 0, 0, 0, 0, 0, t137, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t290, -t663, -t243, t876, qJD(4) * t426 + t541, qJD(4) * t427 + t540, t82, -t59, t139, t137, t876, qJD(4) * t229 - t556, qJD(4) * t230 - t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t290, -t664, t669, t876, t526, t525, t82, -t59, t138, t135, t876, t527, t528; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501 * t816, -t505 * t816; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501 * t590 + t815, -t505 * t590 + t814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t19, t52, t51, t267, qJD(3) * t21 + qJD(4) * t25 - qJD(5) * t5 + t136 - t798, qJD(3) * t22 + qJD(4) * t26 - qJD(5) * t7 + t140 - t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t59, t139, t137, t876, qJD(4) * t280 + t543, qJD(4) * t281 + t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t59, t138, t135, t876, t524, t523; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501 * t817 - t815, t505 * t817 - t814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t20;
