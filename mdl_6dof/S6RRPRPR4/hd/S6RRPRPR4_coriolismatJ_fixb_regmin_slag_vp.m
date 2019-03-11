% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:21
% EndTime: 2019-03-09 10:26:59
% DurationCPUTime: 21.77s
% Computational Cost: add. (25705->792), mult. (65075->1098), div. (0->0), fcn. (75179->12), ass. (0->607)
t541 = cos(pkin(6));
t546 = cos(qJ(4));
t540 = sin(pkin(6));
t544 = sin(qJ(2));
t847 = cos(pkin(11));
t653 = t847 * t544;
t539 = sin(pkin(11));
t547 = cos(qJ(2));
t798 = t539 * t547;
t485 = (-t653 - t798) * t540;
t543 = sin(qJ(4));
t809 = t485 * t543;
t438 = -t541 * t546 - t809;
t440 = -t485 * t546 + t541 * t543;
t538 = sin(pkin(12));
t846 = cos(pkin(12));
t321 = t846 * t438 + t440 * t538;
t545 = cos(qJ(6));
t903 = t545 * t321;
t907 = t903 / 0.2e1;
t908 = -t903 / 0.2e1;
t909 = t907 + t908;
t910 = qJD(6) * t909;
t529 = pkin(2) * t539 + pkin(9);
t771 = qJ(5) + t529;
t508 = t771 * t546;
t649 = t771 * t543;
t573 = t508 * t846 - t538 * t649;
t898 = t573 / 0.2e1;
t899 = -t573 / 0.2e1;
t906 = t899 + t898;
t884 = -t321 / 0.2e1;
t883 = t321 / 0.2e1;
t401 = t508 * t538 + t846 * t649;
t905 = -t401 / 0.2e1;
t650 = t846 * t546;
t511 = t538 * t543 - t650;
t872 = -t511 / 0.2e1;
t651 = t846 * t543;
t800 = t538 * t546;
t514 = t651 + t800;
t870 = -t514 / 0.2e1;
t428 = t846 * t440;
t801 = t538 * t438;
t892 = t428 - t801;
t904 = t892 / 0.2e1;
t542 = sin(qJ(6));
t413 = t542 * t511;
t664 = t413 / 0.2e1;
t642 = t321 * t664;
t197 = t542 * t892;
t525 = t541 * t547 * pkin(1);
t797 = t540 * t544;
t856 = pkin(8) + qJ(3);
t476 = -t797 * t856 + t525;
t448 = t541 * pkin(2) + t476;
t859 = pkin(1) * t544;
t702 = t541 * t859;
t796 = t540 * t547;
t477 = t796 * t856 + t702;
t654 = t847 * t477;
t340 = t539 * t448 + t654;
t334 = pkin(9) * t541 + t340;
t483 = t539 * t797 - t796 * t847;
t690 = -pkin(2) * t547 - pkin(1);
t360 = t483 * pkin(3) + t485 * pkin(9) + t540 * t690;
t211 = t334 * t543 - t546 * t360;
t183 = -qJ(5) * t440 - t211;
t166 = pkin(4) * t483 + t183;
t212 = t334 * t546 + t360 * t543;
t184 = -qJ(5) * t438 + t212;
t802 = t538 * t184;
t78 = t166 * t846 - t802;
t179 = t846 * t184;
t79 = t538 * t166 + t179;
t855 = t78 * t870 + t79 * t872;
t877 = t401 / 0.2e1;
t902 = t321 * t898 - t877 * t892 - t855;
t381 = t543 * t483;
t371 = t381 * t538 - t483 * t650;
t776 = t545 * t371;
t810 = t485 * t542;
t894 = t776 / 0.2e1 - t810 / 0.2e1;
t901 = t79 / 0.2e1;
t536 = t545 ^ 2;
t900 = -t536 / 0.2e1;
t864 = t538 / 0.2e1;
t482 = t483 ^ 2;
t897 = t482 * t546;
t818 = t573 * t542;
t817 = t573 * t545;
t657 = 0.2e1 * t883;
t639 = t657 * t511;
t790 = t542 * t321;
t666 = t790 / 0.2e1;
t273 = t511 * t666;
t415 = t542 * t514;
t662 = t415 / 0.2e1;
t896 = t662 * t892 + t273;
t663 = -t415 / 0.2e1;
t269 = t892 * t663;
t895 = t269 - t642;
t406 = t413 * qJD(6);
t730 = qJD(4) * t545;
t502 = t514 * t730;
t891 = t502 - t406;
t509 = t511 ^ 2;
t510 = t514 ^ 2;
t890 = -t510 - t509;
t534 = t542 ^ 2;
t523 = t536 - t534;
t889 = qJD(4) * t523;
t479 = t485 * t545;
t788 = t542 * t371;
t888 = t479 / 0.2e1 + t788 / 0.2e1;
t737 = qJD(2) * t545;
t684 = t542 * t737;
t643 = t514 * t684;
t887 = 0.2e1 * t643 - t889;
t533 = t540 ^ 2;
t102 = t183 * t846 - t802;
t886 = t102 / 0.2e1;
t453 = t539 * t477;
t339 = t448 * t847 - t453;
t333 = -t541 * pkin(3) - t339;
t881 = t333 / 0.2e1;
t873 = t440 / 0.2e1;
t871 = t511 / 0.2e1;
t869 = t514 / 0.2e1;
t531 = -pkin(2) * t847 - pkin(3);
t518 = -t546 * pkin(4) + t531;
t868 = -t518 / 0.2e1;
t528 = pkin(4) * t538 + pkin(10);
t867 = t528 / 0.2e1;
t530 = -pkin(4) * t846 - pkin(5);
t866 = -t530 / 0.2e1;
t865 = t530 / 0.2e1;
t863 = -t542 / 0.2e1;
t862 = t542 / 0.2e1;
t861 = -t543 / 0.2e1;
t860 = t545 / 0.2e1;
t858 = t440 * pkin(4);
t857 = t543 * pkin(4);
t854 = pkin(2) * qJD(2);
t853 = pkin(4) * qJD(4);
t703 = pkin(2) * t797;
t372 = -pkin(3) * t485 + pkin(9) * t483 + t703;
t369 = t546 * t372;
t384 = t546 * t483;
t363 = t476 * t847 - t453;
t785 = t543 * t363;
t185 = -pkin(4) * t485 + qJ(5) * t384 + t369 - t785;
t355 = t546 * t363;
t368 = t543 * t372;
t768 = -t355 - t368;
t205 = qJ(5) * t381 - t768;
t116 = t185 * t846 - t538 * t205;
t103 = t485 * pkin(5) - t116;
t247 = -t483 * t545 + t197;
t297 = t479 + t788;
t370 = t514 * t483;
t254 = t438 * pkin(4) + t333;
t550 = t321 * pkin(5) - pkin(10) * t892 + t254;
t70 = pkin(10) * t483 + t79;
t39 = t542 * t70 - t545 * t550;
t799 = t539 * t476;
t362 = t654 + t799;
t474 = pkin(4) * t381;
t301 = t362 - t474;
t161 = -pkin(5) * t370 - pkin(10) * t371 + t301;
t782 = t545 * t161;
t117 = t538 * t185 + t846 * t205;
t104 = -pkin(10) * t485 + t117;
t794 = t542 * t104;
t43 = t782 - t794;
t69 = -t483 * pkin(5) - t78;
t7 = t103 * t247 + t297 * t69 + t321 * t43 + t370 * t39;
t852 = t7 * qJD(1);
t249 = t483 * t542 + t545 * t892;
t298 = t776 - t810;
t40 = t542 * t550 + t545 * t70;
t783 = t545 * t104;
t793 = t542 * t161;
t44 = t783 + t793;
t8 = t103 * t249 + t298 * t69 - t321 * t44 + t370 * t40;
t849 = t8 * qJD(1);
t101 = t538 * t183 + t179;
t169 = pkin(5) * t892 + pkin(10) * t321 + t858;
t781 = t545 * t169;
t795 = t542 * t102;
t631 = t781 - t795;
t9 = t101 * t247 + t321 * t631 - t39 * t892 - t69 * t790;
t848 = t9 * qJD(1);
t23 = -t247 * t69 + t321 * t39;
t845 = qJD(1) * t23;
t24 = t249 * t69 - t321 * t40;
t844 = qJD(1) * t24;
t564 = -t370 * t877 + t371 * t898 + t485 * t868;
t606 = t116 * t871 + t117 * t870;
t33 = t564 + t606;
t843 = qJD(1) * t33;
t34 = -t321 * t79 - t78 * t892;
t842 = qJD(1) * t34;
t777 = t545 * t370;
t352 = -t777 / 0.2e1;
t806 = t514 * t247;
t591 = t806 / 0.2e1 + t642;
t579 = t352 + t591;
t51 = t579 + t895;
t841 = qJD(1) * t51;
t696 = t514 * t777;
t825 = t298 * t511;
t93 = t696 + t825;
t840 = qJD(1) * t93;
t784 = t545 * t102;
t792 = t542 * t169;
t632 = t784 + t792;
t10 = t101 * t249 - t321 * t632 - t40 * t892 - t69 * t903;
t839 = t10 * qJD(1);
t838 = t101 * t542;
t837 = t101 * t545;
t18 = -t116 * t892 - t117 * t321 + t370 * t79 - t371 * t78;
t836 = t18 * qJD(1);
t21 = t116 * t78 + t117 * t79 + t254 * t301;
t835 = t21 * qJD(1);
t22 = -t101 * t78 + t102 * t79 + t254 * t858;
t834 = t22 * qJD(1);
t833 = t247 * t370;
t832 = t249 * t892;
t831 = t249 * t370;
t830 = t249 * t511;
t829 = t249 * t542;
t828 = t249 * t545;
t827 = t297 * t321;
t826 = t298 * t321;
t824 = t298 * t542;
t31 = -t254 * t485 + t370 * t78 + t371 * t79;
t823 = t31 * qJD(1);
t822 = t892 * t247;
t821 = t892 * t511;
t820 = t321 * t514;
t819 = t333 * t546;
t816 = t440 * t485;
t780 = t545 * t247;
t601 = t780 / 0.2e1 + t829 / 0.2e1;
t565 = t415 * t903 + t511 * t601;
t605 = t297 * t863 + t298 * t860;
t46 = t565 - t605;
t815 = t46 * qJD(1);
t47 = t211 * t485 + t362 * t438 + (t369 + (-t333 - t363) * t543) * t483;
t814 = t47 * qJD(1);
t48 = t212 * t485 + t362 * t440 + (t768 - t819) * t483;
t813 = t48 * qJD(1);
t812 = t482 * t543;
t811 = t485 * t438;
t808 = t511 * t297;
t807 = t511 * t538;
t805 = t514 * t892;
t804 = t514 * t545;
t803 = t533 * t547;
t791 = t542 * t247;
t789 = t542 * t370;
t787 = t401 * t542;
t418 = pkin(5) * t514 + pkin(10) * t511 + t857;
t786 = t542 * t418;
t775 = t401 * t545;
t774 = t545 * t418;
t77 = -t247 * t298 - t249 * t297;
t773 = t77 * qJD(1);
t119 = -t780 - t829;
t90 = t119 * t321;
t772 = t90 * qJD(1);
t655 = t905 + t877;
t151 = t511 * t906 + t655 * t514;
t770 = t151 * qJD(4);
t524 = -t543 ^ 2 + t546 ^ 2;
t699 = t321 * t790;
t112 = -t699 - t822;
t766 = qJD(1) * t112;
t113 = -t699 + t822;
t765 = qJD(1) * t113;
t698 = t321 * t903;
t114 = -t698 + t832;
t764 = qJD(1) * t114;
t115 = t698 + t832;
t763 = qJD(1) * t115;
t120 = -t827 - t833;
t762 = qJD(1) * t120;
t123 = -t826 - t831;
t761 = qJD(1) * t123;
t676 = t247 * t872;
t592 = t321 * t662 + t676;
t129 = t592 + t894;
t760 = qJD(1) * t129;
t641 = t514 * t908;
t578 = t641 + t888;
t674 = t830 / 0.2e1;
t131 = t674 + t578;
t759 = qJD(1) * t131;
t137 = -t514 * t370 - t371 * t511;
t758 = qJD(1) * t137;
t144 = -t321 * t371 - t370 * t892;
t757 = qJD(1) * t144;
t145 = t211 * t483 - t333 * t438;
t756 = qJD(1) * t145;
t146 = -t212 * t483 + t333 * t440;
t755 = qJD(1) * t146;
t481 = (t653 / 0.2e1 + t798 / 0.2e1) * t540;
t603 = -t821 / 0.2e1 + t820 / 0.2e1;
t159 = t481 + t603;
t754 = qJD(1) * t159;
t182 = t339 * t485 - t340 * t483;
t753 = qJD(1) * t182;
t752 = qJD(1) * t249;
t255 = t811 - t812;
t751 = qJD(1) * t255;
t256 = -t811 - t812;
t750 = qJD(1) * t256;
t257 = -t816 - t897;
t749 = qJD(1) * t257;
t304 = -t816 + t897;
t748 = qJD(1) * t304;
t747 = qJD(1) * t321;
t507 = t651 / 0.2e1 + t800 / 0.2e1;
t361 = t507 * t483;
t746 = qJD(1) * t361;
t745 = qJD(1) * t440;
t744 = qJD(1) * t483;
t743 = qJD(1) * t541;
t742 = qJD(1) * t546;
t741 = qJD(2) * t483;
t740 = qJD(2) * t511;
t739 = qJD(2) * t514;
t738 = qJD(2) * t543;
t736 = qJD(2) * t546;
t735 = qJD(3) * t546;
t417 = t545 * t511;
t734 = qJD(4) * t417;
t733 = qJD(4) * t483;
t732 = qJD(4) * t542;
t731 = qJD(4) * t543;
t729 = qJD(4) * t546;
t728 = qJD(5) * t545;
t727 = qJD(6) * t321;
t726 = qJD(6) * t542;
t725 = qJD(6) * t545;
t109 = -(-t340 + t362) * t485 + (t339 - t363) * t483;
t724 = t109 * qJD(1);
t121 = -t827 + t833;
t723 = t121 * qJD(1);
t122 = t826 - t831;
t722 = t122 * qJD(1);
t675 = -t830 / 0.2e1;
t580 = t545 * t675 + t820 * t900;
t126 = -t824 / 0.2e1 + t580;
t721 = t126 * qJD(1);
t593 = t321 * t663 + t676;
t128 = t593 - t894;
t720 = t128 * qJD(1);
t133 = t675 + t578;
t719 = t133 * qJD(1);
t318 = t370 * t415;
t138 = t318 + t808;
t718 = t138 * qJD(1);
t143 = pkin(2) * t533 * t544 * t690 - t339 * t362 + t340 * t363;
t717 = t143 * qJD(1);
t202 = 0.2e1 * t908;
t716 = t202 * qJD(1);
t629 = t438 * t546 + t440 * t543;
t234 = t629 * t483;
t715 = t234 * qJD(1);
t586 = -t539 * t483 / 0.2e1 + t847 * t485 / 0.2e1;
t346 = (-t797 / 0.2e1 + t586) * pkin(2);
t714 = t346 * qJD(1);
t713 = t361 * qJD(6);
t712 = t381 * qJD(1);
t711 = t384 * qJD(1);
t710 = t809 * qJD(1);
t395 = t485 ^ 2 + t482;
t709 = t395 * qJD(1);
t594 = -pkin(8) * t796 - t702;
t451 = t533 * t859 - t541 * t594;
t708 = t451 * qJD(1);
t506 = pkin(8) * t797 - t525;
t452 = pkin(1) * t803 - t506 * t541;
t707 = t452 * qJD(1);
t706 = t481 * qJD(1);
t513 = (-t544 ^ 2 + t547 ^ 2) * t533;
t705 = t513 * qJD(1);
t704 = t510 - t509;
t701 = t858 / 0.2e1;
t700 = t857 / 0.2e1;
t697 = t544 * t803;
t695 = t69 * t860;
t694 = -t846 / 0.2e1;
t693 = t886 - t78 / 0.2e1;
t692 = t69 / 0.2e1 + t886;
t691 = t901 - t101 / 0.2e1;
t689 = t540 * t743;
t688 = t511 * t739;
t687 = t536 * t739;
t686 = t514 * t737;
t685 = qJD(2) * t540 * t541;
t683 = t542 * t730;
t682 = t514 * t726;
t681 = t514 * t725;
t680 = t483 * t742;
t679 = t485 * t742;
t678 = t511 * t514 * qJD(4);
t677 = t542 * t725;
t673 = t249 * t869;
t672 = t805 / 0.2e1;
t671 = t804 / 0.2e1;
t670 = t321 * t867;
t669 = t247 * t866;
t668 = t249 * t865;
t667 = -t794 / 0.2e1;
t665 = t789 / 0.2e1;
t349 = -t789 / 0.2e1;
t661 = -t783 / 0.2e1;
t658 = t777 / 0.2e1;
t656 = -t355 / 0.2e1 - t368 / 0.2e1;
t652 = t846 * t514;
t648 = qJD(2) + t743;
t647 = -qJD(4) - t744;
t646 = -qJD(6) - t747;
t645 = -qJD(6) - t740;
t644 = qJD(1) * t697;
t640 = t161 / 0.2e1 + t69 * t870;
t638 = (t883 + t884) * t511;
t637 = 0.2e1 * t542 * t502;
t635 = t670 + t169 / 0.2e1;
t556 = t511 * pkin(5) - t514 * pkin(10) + t518;
t282 = -t545 * t556 + t818;
t549 = t282 * t904 + t247 * t899 + t631 * t872 + (t774 + t787) * t884 + t401 * t666 + t69 * t664;
t559 = -t103 * t545 / 0.2e1 + t297 * t865 + t528 * t665;
t1 = (t39 / 0.2e1 - t838 / 0.2e1) * t514 + t549 + t559;
t110 = (-t282 + t818) * t514 + t774 * t511;
t634 = -t1 * qJD(1) + t110 * qJD(2);
t17 = (t101 - t79) * t892 + (-t102 + t78) * t321;
t91 = (-t892 / 0.2e1 + t904) * t514 + t638;
t633 = t17 * qJD(1) + t91 * qJD(3);
t630 = -t321 * t530 - t528 * t892;
t628 = -t483 * t531 + t485 * t529;
t627 = -t511 * t530 - t514 * t528;
t283 = t542 * t556 + t817;
t567 = t249 * t905 + t283 * t883 + t40 * t871;
t13 = t545 * t640 + t567 + t667;
t191 = -t283 * t511 + t401 * t804;
t626 = qJD(1) * t13 - qJD(2) * t191;
t135 = t321 ^ 2 + t892 ^ 2;
t92 = 0.2e1 * t514 * t904 + t639;
t625 = qJD(1) * t135 + qJD(2) * t92;
t568 = t247 * t877 + t282 * t884 + t39 * t872;
t14 = -t542 * t640 + t568 + t661;
t190 = t282 * t511 - t401 * t415;
t624 = -qJD(1) * t14 + qJD(2) * t190;
t577 = (-t370 * t694 + t371 * t864) * pkin(4);
t19 = t511 * t691 - t514 * t693 + t577;
t623 = t19 * qJD(1) - t151 * qJD(2);
t235 = t401 * t514 - t511 * t573;
t600 = t654 / 0.2e1 - t474 / 0.2e1 + t799 / 0.2e1;
t29 = t600 + t902;
t622 = -qJD(1) * t29 + qJD(2) * t235;
t553 = (-t805 / 0.2e1 + t638) * t545 + t673;
t58 = t665 + t553;
t621 = qJD(1) * t58;
t620 = qJD(1) * t91;
t378 = t704 * t542;
t566 = t269 + t273 - t806 / 0.2e1 + t642;
t50 = t658 + t566;
t619 = -qJD(1) * t50 + qJD(2) * t378;
t379 = t890 * t542;
t557 = t658 + t591;
t54 = t557 + t896;
t618 = qJD(1) * t54 - qJD(2) * t379;
t380 = t704 * t545;
t552 = (t672 - t639) * t545 + t673;
t56 = t665 + t552;
t617 = -qJD(1) * t56 - qJD(2) * t380;
t426 = t890 * t545;
t551 = (t672 + t639) * t545 + t673;
t60 = t349 + t551;
t616 = qJD(1) * t60 - qJD(2) * t426;
t615 = qJD(1) * t92 - qJD(2) * t890;
t614 = t645 * t545;
t588 = -t321 * t864 + t694 * t892;
t164 = (-t440 / 0.2e1 + t588) * pkin(4);
t587 = -t807 / 0.2e1 - t652 / 0.2e1;
t390 = (t861 + t587) * pkin(4);
t613 = qJD(1) * t164 + qJD(2) * t390;
t612 = qJD(1) * t197 + qJD(2) * t415;
t198 = t657 * t542;
t162 = qJD(1) * t198 + qJD(2) * t413;
t201 = 0.2e1 * t907;
t611 = -qJD(1) * t201 - qJD(2) * t417;
t260 = t438 ^ 2 - t440 ^ 2;
t610 = qJD(1) * t260 - qJD(2) * t629;
t609 = -qJD(1) * t629 + qJD(2) * t524;
t317 = -t801 / 0.2e1 + t428 / 0.2e1;
t608 = qJD(1) * t317 + qJD(2) * t507;
t607 = t641 - t888;
t602 = t511 * t867 + t514 * t866;
t558 = t819 / 0.2e1 - t531 * t438 / 0.2e1 + t529 * t381 / 0.2e1;
t139 = t558 - t656;
t599 = -qJD(1) * t139 - t531 * t736;
t590 = t531 * t873 - t529 * t384 / 0.2e1;
t141 = -t369 / 0.2e1 + (t881 + t363 / 0.2e1) * t543 + t590;
t598 = -qJD(1) * t141 - t531 * t738;
t319 = t438 * t861 + t546 * t873;
t597 = -qJD(2) * t319 + t438 * t745;
t596 = qJD(1) * t319 + t543 * t736;
t595 = -qJD(4) * t481 + t485 * t744;
t589 = t117 * t864 + t116 * t846 / 0.2e1;
t585 = t418 / 0.2e1 + t602;
t111 = (-t283 + t817) * t514 - t786 * t511;
t548 = t283 * t904 + t249 * t899 + t632 * t871 + (-t775 + t786) * t883 + t401 * t907 + t511 * t695;
t560 = t103 * t862 + t298 * t865 + t528 * t658;
t2 = (t40 / 0.2e1 - t837 / 0.2e1) * t514 + t548 + t560;
t584 = -t2 * qJD(1) + t111 * qJD(2);
t175 = t518 * t857;
t555 = t101 * t905 + t102 * t899 + t401 * t901 + t78 * t898;
t5 = (t254 * t861 + t440 * t868 + t589) * pkin(4) + t555;
t583 = -t5 * qJD(1) + t175 * qJD(2) + t151 * qJD(3);
t554 = t906 * t892;
t576 = (t370 * t864 + t371 * t694) * pkin(4);
t11 = t511 * t693 + t514 * t691 + t554 + t576;
t582 = t11 * qJD(1);
t124 = t247 ^ 2 - t249 ^ 2;
t87 = (-t791 + t828) * t514;
t581 = qJD(1) * t124 - qJD(2) * t87 + qJD(4) * t119;
t170 = t542 * t585 + t545 * t655;
t25 = t542 * t635 + t545 * t692 + t669;
t575 = -qJD(1) * t25 - qJD(2) * t170 - t530 * t730;
t172 = t542 * t655 - t545 * t585;
t27 = t542 * t692 - t545 * t635 + t668;
t574 = -qJD(1) * t27 - qJD(2) * t172 - t530 * t732;
t149 = t601 * t514;
t154 = -t791 / 0.2e1 + t828 / 0.2e1;
t572 = qJD(2) * t149 - qJD(4) * t154 + t247 * t752;
t411 = (t534 / 0.2e1 + t900) * t514;
t571 = qJD(1) * t154 - qJD(2) * t411 + t683;
t176 = t821 / 0.2e1 + t321 * t869;
t570 = qJD(2) * t176 + qJD(6) * t317 + t747 * t892;
t569 = qJD(1) * t176 + qJD(6) * t507 + t688;
t425 = t523 * t510;
t563 = qJD(1) * t87 + qJD(2) * t425 + t637;
t562 = -qJD(1) * t119 + t887;
t561 = qJD(1) * t149 + qJD(4) * t411 + t510 * t684;
t503 = t507 * qJD(4);
t501 = t514 * t732;
t480 = -0.2e1 * t514 * t677;
t475 = t481 * qJD(2);
t473 = t485 * t736;
t410 = t417 * qJD(6);
t405 = t413 * qJD(4);
t404 = t411 * qJD(6);
t389 = pkin(4) * t587 + t700;
t375 = t381 * qJD(4);
t357 = -t712 - t731;
t345 = t703 / 0.2e1 + t586 * pkin(2);
t314 = t319 * qJD(4);
t313 = t892 * t730;
t261 = (qJD(1) * t892 + t739) * t545;
t236 = t629 * qJD(4);
t210 = -qJD(2) * t361 + qJD(4) * t317;
t199 = t321 * t863 + t666;
t195 = t199 * qJD(6);
t194 = t198 * qJD(6);
t174 = t176 * qJD(4);
t173 = t401 * t862 + t787 / 0.2e1 + t774 / 0.2e1 - t602 * t545;
t171 = t401 * t860 + t775 / 0.2e1 - t786 / 0.2e1 + t602 * t542;
t165 = pkin(4) * t588 + t701;
t160 = t481 - t603;
t155 = -t162 - t726;
t153 = t154 * qJD(6);
t148 = t149 * qJD(6);
t142 = t543 * t881 - t785 / 0.2e1 + t369 / 0.2e1 + t590;
t140 = t558 + t656;
t134 = t675 + t607;
t132 = t674 + t607;
t130 = t592 - t894;
t127 = t593 + t894;
t125 = t824 / 0.2e1 + t580;
t118 = t119 * qJD(6);
t89 = t92 * qJD(5);
t88 = t91 * qJD(4);
t82 = t87 * qJD(6);
t59 = t665 + t551;
t57 = t349 + t553;
t55 = t349 + t552;
t53 = t579 + t896;
t52 = t557 + t895;
t49 = t352 + t566;
t45 = t565 + t605;
t32 = t564 - t606;
t30 = t600 - t902;
t28 = t528 * t908 + t668 + t69 * t862 - t795 / 0.2e1 + t781 / 0.2e1;
t26 = t542 * t670 + t669 + t695 - t784 / 0.2e1 - t792 / 0.2e1;
t20 = t101 * t871 + t102 * t869 + t577 + t855;
t16 = t69 * t671 + t667 + t782 / 0.2e1 - t567;
t15 = t69 * t663 + t661 - t793 / 0.2e1 - t568;
t12 = t101 * t869 + t102 * t872 + t78 * t871 + t79 * t870 - t554 + t576;
t6 = pkin(4) * t589 + t254 * t700 + t518 * t701 - t555;
t4 = t101 * t671 + t40 * t870 - t548 + t560;
t3 = t101 * t662 + t39 * t870 - t549 + t559;
t35 = [0, 0, 0, qJD(2) * t697, t513 * qJD(2), t547 * t685, -t544 * t685, 0, -t451 * qJD(2), -t452 * qJD(2), qJD(2) * t109 + qJD(3) * t395, qJD(2) * t143 + qJD(3) * t182 (-qJD(4) * t438 - t483 * t736) * t440, qJD(2) * t234 + qJD(4) * t260, qJD(2) * t257 - t438 * t733, -qJD(2) * t256 - t440 * t733, -t485 * t741, qJD(2) * t47 - qJD(3) * t255 + qJD(4) * t146, qJD(2) * t48 + qJD(3) * t304 + qJD(4) * t145, qJD(2) * t18 + qJD(3) * t144 + qJD(4) * t17 + qJD(5) * t135, qJD(2) * t21 + qJD(3) * t31 + qJD(4) * t22 + qJD(5) * t34 (qJD(2) * t298 - qJD(6) * t247 - t321 * t730) * t249, qJD(2) * t77 - qJD(4) * t90 + qJD(6) * t124, qJD(2) * t122 + qJD(4) * t114 - t247 * t727, qJD(2) * t121 - qJD(4) * t113 - t249 * t727 (-qJD(2) * t370 + qJD(4) * t892) * t321, qJD(2) * t7 + qJD(3) * t120 + qJD(4) * t9 - qJD(5) * t112 + qJD(6) * t24, qJD(2) * t8 + qJD(3) * t123 + qJD(4) * t10 + qJD(5) * t115 + qJD(6) * t23; 0, 0, 0, t644, t705, t648 * t796, -t648 * t797, 0, qJD(2) * t594 - t708, qJD(2) * t506 - t707, t724 + (t483 * t847 + t485 * t539) * t854, t717 + (-t362 * t847 + t363 * t539) * t854 + t345 * qJD(3), t314 + (-t738 - t745) * t384, -t524 * t741 - t236 + t715, -t485 * t738 + t749, -t473 - t750, -t595, t814 + (-t362 * t546 + t543 * t628) * qJD(2) + t142 * qJD(4), t813 + (t362 * t543 + t546 * t628) * qJD(2) + t140 * qJD(4), t836 + (-t116 * t514 - t117 * t511 + t370 * t573 + t371 * t401) * qJD(2) + t12 * qJD(4) + t89, t835 + (-t116 * t401 + t117 * t573 + t301 * t518) * qJD(2) + t32 * qJD(3) + t6 * qJD(4) + t30 * qJD(5), t125 * qJD(4) - t148 + (t686 + t752) * t298, t773 + t45 * qJD(4) - t82 + (-t297 * t545 - t824) * t739, t722 + (-t696 + t825) * qJD(2) + t55 * qJD(4) + t127 * qJD(6), t723 + (t318 - t808) * qJD(2) + t49 * qJD(4) + t134 * qJD(6), -t713 + t174 - (t740 + t747) * t370, t852 + (t103 * t415 + t282 * t370 + t297 * t401 + t43 * t511) * qJD(2) + t3 * qJD(4) + t53 * qJD(5) + t16 * qJD(6), t849 + (t103 * t804 + t283 * t370 + t298 * t401 - t44 * t511) * qJD(2) + t4 * qJD(4) + t59 * qJD(5) + t15 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t709, qJD(2) * t345 + t753, 0, 0, 0, 0, 0, -t751, t748, t88 + t757, t823 + t32 * qJD(2) + (-t370 * t511 + t371 * t514) * qJD(3) + t20 * qJD(4) + t160 * qJD(5), 0, 0, 0, 0, 0, qJD(4) * t52 + qJD(6) * t132 + t762, qJD(4) * t57 + qJD(6) * t130 + t761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t597, t610, t647 * t438, t647 * t440, t475, qJD(2) * t142 - qJD(4) * t212 + t755, qJD(2) * t140 + qJD(4) * t211 + t756, t12 * qJD(2) + (t321 * t846 - t538 * t892) * t853 + t633, t834 + t6 * qJD(2) + t20 * qJD(3) + (-t101 * t846 + t102 * t538) * t853 + t165 * qJD(5), t125 * qJD(2) + t153 - (t732 + t752) * t903, t45 * qJD(2) - t321 * t889 + t118 - t772, qJD(2) * t55 + t732 * t892 + t764 + t910, qJD(2) * t49 + t195 + t313 - t765, t570, t848 + t3 * qJD(2) + t52 * qJD(3) + (t542 * t630 - t837) * qJD(4) + t28 * qJD(6), t839 + t4 * qJD(2) + t57 * qJD(3) + (t545 * t630 + t838) * qJD(4) + t26 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t625, qJD(2) * t30 + qJD(3) * t160 + qJD(4) * t165 + t842, 0, 0, 0, 0, 0, qJD(2) * t53 + t195 - t766, qJD(2) * t59 + t763 + t910; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t572, t581, qJD(2) * t127 + qJD(4) * t909 + t247 * t646, qJD(2) * t134 + qJD(4) * t199 + t249 * t646, t210, qJD(2) * t16 + qJD(3) * t132 + qJD(4) * t28 + qJD(5) * t199 - qJD(6) * t40 + t844, qJD(2) * t15 + qJD(3) * t130 + qJD(4) * t26 + qJD(5) * t909 + qJD(6) * t39 + t845; 0, 0, 0, -t644, -t705, -t547 * t689, t544 * t689, 0, t708, t707, -t724, qJD(3) * t346 - t717, t440 * t680 + t314, -t236 - t715, qJD(4) * t384 - t749, -t375 + t750, t595, qJD(4) * t141 + t485 * t735 - t814, -qJD(3) * t809 + qJD(4) * t139 - t813, qJD(3) * t137 - qJD(4) * t11 - t836 + t89, qJD(3) * t33 - qJD(4) * t5 - qJD(5) * t29 - t835, qJD(4) * t126 - t298 * t752 - t148, qJD(4) * t46 - t773 - t82, qJD(4) * t56 + qJD(6) * t128 - t722, qJD(4) * t50 + qJD(6) * t133 - t723, t370 * t747 + t174 + t713, -qJD(3) * t138 - qJD(4) * t1 + qJD(5) * t54 - qJD(6) * t13 - t852, -qJD(3) * t93 - qJD(4) * t2 + qJD(5) * t60 - qJD(6) * t14 - t849; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t543 * t729, t524 * qJD(4), 0, 0, 0, t531 * t731, t531 * t729, -qJD(5) * t890, qJD(4) * t175 + qJD(5) * t235, -t510 * t677 - t536 * t678, -qJD(6) * t425 + t511 * t637, qJD(4) * t380 - t511 * t682, -qJD(4) * t378 - t511 * t681, t678, qJD(4) * t110 - qJD(5) * t379 + qJD(6) * t191, qJD(4) * t111 - qJD(5) * t426 + qJD(6) * t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, 0, 0, 0, 0, 0, t679, -t710, t758, t770 + t843, 0, 0, 0, 0, 0, -t718, -t840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, t609, t711 + t729, t357, -t706, -t529 * t729 - t598, t529 * t731 - t599 (t511 * t846 - t514 * t538) * t853 - t582 (-t401 * t538 - t573 * t846) * t853 + t389 * qJD(5) + t583, t721 - t404 + (-t683 - t687) * t511, t511 * t887 + t480 + t815, t501 - t617, t502 - t619, t569 (t542 * t627 - t817) * qJD(4) + t173 * qJD(6) + t634 (t545 * t627 + t818) * qJD(4) + t171 * qJD(6) + t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t615, qJD(4) * t389 + t622, 0, 0, 0, 0, 0, t618, t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t561, -t563, t415 * t645 + t720, t514 * t614 + t719, t503 + t746, qJD(4) * t173 - qJD(6) * t283 - t626, qJD(4) * t171 + qJD(6) * t282 + t624; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t709, -qJD(2) * t346 - t753, 0, 0, 0, 0, 0, -t375 - t473 + t751, qJD(2) * t809 - t483 * t729 - t748, -qJD(2) * t137 - t757 + t88, -qJD(2) * t33 - qJD(4) * t19 - qJD(5) * t159 - t823, 0, 0, 0, 0, 0, qJD(2) * t138 + qJD(4) * t51 + qJD(6) * t131 - t762, qJD(2) * t93 + qJD(4) * t58 + qJD(6) * t129 - t761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t714, 0, 0, 0, 0, 0, -t679, t710, -t758, t770 - t843, 0, 0, 0, 0, 0, t718, t840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t357, t647 * t546, t620 (-t652 - t807) * t853 - t623, 0, 0, 0, 0, 0, t841 - t891, t410 + t501 + t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t754, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t405 - t681 + t759, t682 + t734 + t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t597, -t610, -qJD(2) * t384 + t438 * t744, qJD(2) * t381 + t440 * t744, t475, -qJD(2) * t141 + qJD(3) * t381 - t755, -qJD(2) * t139 + t483 * t735 - t756, qJD(2) * t11 - t633, qJD(2) * t5 + qJD(3) * t19 + qJD(5) * t164 - t834, -qJD(2) * t126 + t752 * t903 + t153, -qJD(2) * t46 + t118 + t772, -qJD(2) * t56 + qJD(6) * t201 - t764, -qJD(2) * t50 - t194 + t765, -t570, qJD(2) * t1 - qJD(3) * t51 + qJD(6) * t27 - t728 * t892 - t848, qJD(2) * t2 - qJD(3) * t58 + qJD(5) * t197 + qJD(6) * t25 - t839; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t596, -t609, -t711, t712, t706, t598, t599, t582, qJD(5) * t390 - t583, t511 * t687 - t404 - t721, -0.2e1 * t511 * t643 + t480 - t815, t410 + t617, -t406 + t619, -t569, qJD(6) * t172 - t514 * t728 - t634, qJD(5) * t415 + qJD(6) * t170 - t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t712, t680, -t620, t623, 0, 0, 0, 0, 0, -t841, -t621; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t677, t523 * qJD(6), 0, 0, 0, t530 * t726, t530 * t725; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t613, 0, 0, 0, 0, 0, -t261, t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t571, -t562, -t611 + t725, t155, -t608, -t528 * t725 - t574, t528 * t726 - t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t625, qJD(2) * t29 + qJD(3) * t159 - qJD(4) * t164 - t842, 0, 0, 0, 0, 0, -qJD(2) * t54 - t194 + t313 + t766, -qJD(2) * t60 - qJD(4) * t197 + qJD(6) * t202 - t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t615, -qJD(4) * t390 - t622, 0, 0, 0, 0, 0, -t618 + t891, -qJD(4) * t415 - t511 * t725 - t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t754, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t613, 0, 0, 0, 0, 0, t261, -t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t614 + t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t572, -t581, -qJD(2) * t128 - qJD(4) * t201 + t247 * t747, -qJD(2) * t133 + qJD(4) * t198 + t249 * t747, t210, qJD(2) * t13 - qJD(3) * t131 - qJD(4) * t27 + qJD(5) * t198 - t844, qJD(2) * t14 - qJD(3) * t129 - qJD(4) * t25 - qJD(5) * t202 - t845; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t561, t563, t542 * t688 - t720 - t734, t511 * t686 + t405 - t719, t503 - t746, -qJD(4) * t172 + qJD(5) * t413 + t626, -qJD(4) * t170 + t511 * t728 - t624; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t759, -t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t571, t562, t611, t162, t608, t574, t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t511 * t737 - t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t35;
