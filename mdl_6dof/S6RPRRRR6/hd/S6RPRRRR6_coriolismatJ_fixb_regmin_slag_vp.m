% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:38
% EndTime: 2019-03-09 07:15:24
% DurationCPUTime: 30.01s
% Computational Cost: add. (20781->696), mult. (43116->891), div. (0->0), fcn. (51651->10), ass. (0->566)
t947 = qJD(5) + qJD(6);
t956 = qJD(4) + t947;
t542 = cos(qJ(4));
t878 = pkin(8) + pkin(9);
t515 = t878 * t542;
t541 = cos(qJ(5));
t508 = t541 * t515;
t538 = sin(qJ(4));
t514 = t878 * t538;
t537 = sin(qJ(5));
t769 = t537 * t514;
t432 = t508 - t769;
t751 = t541 * t542;
t504 = t537 * t538 - t751;
t362 = -pkin(10) * t504 + t432;
t536 = sin(qJ(6));
t540 = cos(qJ(6));
t752 = t541 * t538;
t768 = t537 * t542;
t506 = t752 + t768;
t628 = t541 * t514 + t537 * t515;
t888 = -t506 * pkin(10) - t628;
t259 = t540 * t362 + t536 * t888;
t953 = t956 * t259;
t950 = t536 * t362 - t540 * t888;
t952 = t956 * t950;
t534 = sin(pkin(11));
t535 = cos(pkin(11));
t539 = sin(qJ(3));
t858 = cos(qJ(3));
t501 = t534 * t539 - t535 * t858;
t570 = t950 * t501;
t948 = -t570 / 0.2e1;
t866 = -t501 / 0.2e1;
t951 = t259 * t866;
t629 = t540 * t504 + t536 * t506;
t679 = qJD(4) + qJD(5);
t884 = t629 * qJD(6);
t918 = t629 * t679 + t884;
t369 = t504 * t501;
t755 = t540 * t369;
t366 = t506 * t501;
t777 = t536 * t366;
t736 = t777 / 0.2e1 + t755 / 0.2e1;
t904 = t629 * t501;
t922 = t904 / 0.2e1;
t929 = t922 + t736;
t935 = t929 * qJD(1);
t83 = -t935 - t918;
t923 = -t904 / 0.2e1;
t928 = t923 - t736;
t936 = t928 * qJD(1);
t84 = -t936 + t918;
t934 = qJD(6) + t679;
t927 = t923 + t736;
t937 = t927 * qJD(3);
t150 = t928 * qJD(2);
t930 = t922 - t736;
t156 = t930 * qJD(2);
t662 = t858 * t534;
t765 = t539 * t535;
t503 = t662 + t765;
t385 = t538 * t503;
t844 = pkin(7) + qJ(2);
t512 = t844 * t534;
t513 = t844 * t535;
t427 = t858 * t512 + t513 * t539;
t322 = pkin(4) * t385 + t427;
t527 = -pkin(4) * t542 - pkin(3);
t863 = -t506 / 0.2e1;
t365 = t385 * t537 - t503 * t751;
t871 = t365 / 0.2e1;
t933 = -t322 * t863 + t432 * t866 - t527 * t871;
t298 = t503 * t506;
t630 = -t298 * t540 + t536 * t365;
t684 = t501 * qJD(1);
t71 = qJD(3) * t929 - t630 * t684;
t179 = t298 * t536 + t540 * t365;
t911 = -t179 ^ 2 + t630 ^ 2;
t926 = qJD(1) * t911;
t350 = t540 * t366;
t775 = t536 * t369;
t610 = t350 / 0.2e1 - t775 / 0.2e1;
t486 = t540 * t506;
t774 = t536 * t504;
t889 = t486 - t774;
t900 = t889 * t501;
t913 = -t900 / 0.2e1 + t610;
t925 = qJD(3) * t913;
t912 = t900 / 0.2e1 + t610;
t924 = t912 * qJD(1);
t149 = t912 * qJD(2);
t155 = t913 * qJD(2);
t388 = t542 * t501;
t849 = t503 * pkin(3);
t850 = t501 * pkin(8);
t413 = t849 + t850;
t396 = t542 * t413;
t798 = t427 * t538;
t848 = t503 * pkin(4);
t238 = pkin(9) * t388 + t396 + t798 + t848;
t231 = t541 * t238;
t383 = t538 * t501;
t395 = t538 * t413;
t412 = t427 * t542;
t734 = t412 - t395;
t270 = pkin(9) * t383 - t734;
t770 = t537 * t270;
t611 = t231 / 0.2e1 - t770 / 0.2e1;
t68 = t611 + t933;
t66 = t611 - t933;
t70 = qJD(3) * t912 - t179 * t684;
t261 = pkin(5) * t298 + t322;
t441 = t504 * pkin(5) + t527;
t868 = t889 / 0.2e1;
t874 = -t179 / 0.2e1;
t891 = t261 * t868 + t441 * t874;
t870 = -t629 / 0.2e1;
t875 = t630 / 0.2e1;
t890 = t261 * t870 + t441 * t875;
t921 = t179 * t679;
t820 = t261 * t179;
t488 = t662 / 0.2e1 + t765 / 0.2e1;
t683 = t503 * qJD(1);
t419 = t501 * t683;
t732 = t488 * qJD(4) + t419;
t917 = t488 * qJD(5) + t732;
t898 = t629 ^ 2 - t889 ^ 2;
t916 = qJD(3) * t898;
t915 = t179 * qJD(6);
t914 = t527 * t298 / 0.2e1 + t322 * t504 / 0.2e1 + t628 * t866;
t819 = t261 * t630;
t907 = t179 * t889;
t716 = qJD(3) * t889;
t903 = t629 * t716;
t722 = qJD(1) * t179;
t901 = t630 * t722;
t230 = t537 * t238;
t260 = t541 * t270;
t638 = -t230 / 0.2e1 - t260 / 0.2e1;
t67 = t638 + t914;
t69 = t638 - t914;
t897 = t501 * t679;
t896 = t504 * t679;
t892 = t630 * t679;
t492 = t501 ^ 2;
t493 = t503 ^ 2;
t887 = -t493 - t492;
t677 = t493 - t492;
t525 = -pkin(2) * t535 - pkin(1);
t616 = pkin(3) * t501 - pkin(8) * t503;
t381 = t525 + t616;
t428 = -t539 * t512 + t858 * t513;
t767 = t538 * t428;
t273 = -t542 * t381 + t767;
t789 = t503 * t542;
t252 = -pkin(9) * t789 - t273;
t753 = t541 * t252;
t750 = t542 * t428;
t274 = t538 * t381 + t750;
t253 = -pkin(9) * t385 + t274;
t771 = t537 * t253;
t127 = t753 - t771;
t852 = t365 * pkin(10);
t112 = t127 + t852;
t851 = t501 * pkin(4);
t219 = t252 + t851;
t251 = t541 * t253;
t123 = t537 * t219 + t251;
t855 = pkin(10) * t298;
t103 = t123 - t855;
t772 = t537 * t252;
t126 = -t251 - t772;
t111 = t126 + t855;
t639 = t111 / 0.2e1 + t103 / 0.2e1;
t614 = t639 * t540;
t845 = t541 * pkin(4);
t663 = pkin(5) + t845;
t617 = t536 * t663;
t754 = t540 * t537;
t482 = pkin(4) * t754 + t617;
t792 = t482 * t501;
t216 = t541 * t219;
t122 = -t216 + t771;
t588 = -t122 + t852;
t86 = t501 * pkin(5) + t588;
t15 = t792 / 0.2e1 + t614 + (-t112 / 0.2e1 + t86 / 0.2e1) * t536;
t642 = -t486 / 0.2e1;
t393 = t642 + t486 / 0.2e1;
t691 = t393 * qJD(2);
t886 = qJD(1) * t15 + qJD(4) * t482 - t691;
t433 = t536 * pkin(5);
t565 = t536 * t588;
t553 = t565 / 0.2e1;
t669 = -t433 / 0.2e1;
t764 = t540 * t103;
t645 = -t764 / 0.2e1;
t833 = t536 * t86;
t843 = -t833 / 0.2e1 + t645;
t591 = t501 * t669 + t843;
t7 = t553 + t764 / 0.2e1 + t591;
t885 = t7 * qJD(1) - t433 * qJD(4) + t691;
t883 = qJD(6) * t488 + t917;
t532 = t538 ^ 2;
t533 = t542 ^ 2;
t522 = t533 - t532;
t766 = t538 * t542;
t632 = 0.2e1 * t503 * t766;
t561 = qJD(1) * t632 - t522 * qJD(3);
t879 = t679 * t628;
t877 = -t216 / 0.2e1;
t876 = -t219 / 0.2e1;
t641 = -t251 / 0.2e1;
t869 = -t889 / 0.2e1;
t856 = pkin(5) * t506;
t857 = pkin(4) * t538;
t443 = t856 + t857;
t867 = -t443 / 0.2e1;
t864 = -t503 / 0.2e1;
t862 = t506 / 0.2e1;
t640 = -t508 / 0.2e1;
t861 = t538 / 0.2e1;
t860 = t541 / 0.2e1;
t859 = -t542 / 0.2e1;
t853 = t365 * pkin(5);
t847 = t503 * pkin(5);
t846 = t537 * pkin(4);
t786 = t536 * t103;
t648 = t786 / 0.2e1;
t85 = t540 * t86;
t842 = -t85 / 0.2e1 + t648;
t841 = -t565 / 0.2e1 + t645;
t564 = t540 * t588;
t840 = t648 - t564 / 0.2e1;
t839 = pkin(4) * qJD(4);
t838 = pkin(4) * qJD(5);
t837 = pkin(5) * qJD(5);
t836 = pkin(5) * qJD(6);
t323 = -pkin(4) * t383 + t428;
t262 = -pkin(5) * t366 + t323;
t264 = -t350 + t775;
t41 = -t85 + t786;
t743 = t260 + t230;
t110 = pkin(10) * t366 + t743;
t785 = t536 * t110;
t631 = t231 - t770;
t97 = -pkin(10) * t369 + t631 + t847;
t831 = t540 * t97;
t1 = (-t785 + t831) * t501 - t41 * t503 - t262 * t630 + t261 * t264;
t835 = t1 * qJD(1);
t268 = t755 + t777;
t42 = t764 + t833;
t763 = t540 * t110;
t832 = t536 * t97;
t2 = -(t763 + t832) * t501 - t42 * t503 - t262 * t179 + t261 * t268;
t834 = t2 * qJD(1);
t43 = -t565 - t764;
t29 = t43 * t501 + t630 * t853 - t820;
t830 = qJD(1) * t29;
t44 = t564 - t786;
t30 = t179 * t853 - t44 * t501 + t819;
t829 = qJD(1) * t30;
t676 = pkin(4) * t789;
t294 = t676 - t853;
t762 = t540 * t111;
t783 = t536 * t112;
t48 = t762 - t783;
t31 = -t294 * t630 + t48 * t501 - t820;
t828 = qJD(1) * t31;
t761 = t540 * t112;
t784 = t536 * t111;
t49 = t761 + t784;
t32 = -t179 * t294 - t49 * t501 + t819;
t827 = qJD(1) * t32;
t33 = t41 * t501 + t819;
t826 = qJD(1) * t33;
t34 = -t42 * t501 - t820;
t825 = qJD(1) * t34;
t811 = t322 * t298;
t64 = -t127 * t501 - t365 * t676 - t811;
t824 = qJD(1) * t64;
t812 = t322 * t365;
t65 = t126 * t501 + t298 * t676 - t812;
t823 = qJD(1) * t65;
t74 = t122 * t501 - t811;
t822 = qJD(1) * t74;
t75 = -t123 * t501 - t812;
t821 = qJD(1) * t75;
t816 = t264 * t501;
t815 = t630 * t503;
t814 = t268 * t501;
t813 = t179 * t503;
t808 = t365 * t503;
t807 = t366 * t501;
t806 = t298 * t503;
t805 = t369 * t501;
t37 = -t123 * t503 + t322 * t369 - t323 * t365 - t501 * t743;
t804 = t37 * qJD(1);
t38 = -t122 * t503 + t298 * t323 - t322 * t366 + t501 * t631;
t803 = t38 * qJD(1);
t516 = t540 * t663;
t523 = t536 * t846;
t481 = t523 - t516;
t793 = t481 * t501;
t773 = t536 * t541;
t489 = (t754 + t773) * pkin(4);
t791 = t489 * t501;
t490 = t540 * t845 - t523;
t790 = t490 * t501;
t671 = -t851 / 0.2e1;
t618 = t671 + t252 / 0.2e1;
t58 = t541 * t618 + t877;
t749 = t58 * qJD(1);
t63 = t179 * t264 + t268 * t630;
t748 = t63 * qJD(1);
t87 = t396 * t501 + (-t273 + t767) * t503;
t747 = t87 * qJD(1);
t88 = (-t274 + t750) * t503 + (t734 - t412) * t501;
t746 = t88 * qJD(1);
t145 = t298 * t504 + t506 * t365;
t745 = t679 * t145;
t213 = -t298 * t862 + t504 * t871;
t744 = t679 * t213;
t579 = t768 / 0.2e1 + t752 / 0.2e1;
t286 = (t862 + t579) * t501;
t742 = t679 * t286;
t291 = (t579 + t863) * t501;
t741 = t679 * t291;
t666 = t537 * t383;
t731 = -t666 / 0.2e1 + t388 * t860;
t651 = t541 * t866;
t730 = t666 / 0.2e1 + t542 * t651;
t519 = t534 ^ 2 + t535 ^ 2;
t118 = t815 - t816;
t729 = qJD(1) * t118;
t119 = -t815 - t816;
t728 = qJD(1) * t119;
t121 = -t813 - t814;
t727 = qJD(1) * t121;
t177 = -t274 * t501 + t427 * t789;
t726 = qJD(1) * t177;
t209 = -t806 + t807;
t725 = qJD(1) * t209;
t210 = t806 + t807;
t724 = qJD(1) * t210;
t212 = -t805 - t808;
t723 = qJD(1) * t212;
t317 = t677 * t538;
t721 = qJD(1) * t317;
t318 = t887 * t538;
t720 = qJD(1) * t318;
t319 = t677 * t542;
t719 = qJD(1) * t319;
t718 = qJD(1) * t365;
t717 = qJD(2) * t542;
t715 = qJD(3) * t441;
t714 = qJD(3) * t506;
t713 = qJD(3) * t527;
t712 = qJD(3) * t542;
t711 = qJD(4) * t538;
t710 = qJD(4) * t542;
t709 = qJD(5) * t527;
t120 = -t813 + t814;
t708 = t120 * qJD(1);
t146 = -t298 * t369 - t365 * t366;
t707 = t146 * qJD(1);
t176 = t273 * t501 - t385 * t427;
t702 = t176 * qJD(1);
t211 = t805 - t808;
t701 = t211 * qJD(1);
t240 = t503 * t629;
t700 = t240 * qJD(1);
t242 = t503 * t889;
t699 = t242 * qJD(1);
t276 = t286 * qJD(1);
t649 = t369 / 0.2e1;
t287 = t649 + t730;
t277 = t287 * qJD(1);
t650 = -t369 / 0.2e1;
t288 = t650 + t731;
t279 = t288 * qJD(1);
t296 = t503 * t504;
t698 = t296 * qJD(1);
t697 = t298 * qJD(1);
t695 = t677 * qJD(1);
t694 = t383 * qJD(1);
t693 = t385 * qJD(1);
t692 = t388 * qJD(1);
t394 = t887 * t542;
t690 = t394 * qJD(1);
t688 = t889 * qJD(6);
t687 = t488 * qJD(1);
t487 = t501 * qJD(3);
t682 = t503 * qJD(3);
t509 = t519 * qJ(2);
t681 = t509 * qJD(1);
t680 = t519 * qJD(1);
t47 = -t629 * t875 - t179 * t869 + t630 * t870 + t907 / 0.2e1;
t61 = -t629 * t630 + t907;
t678 = t61 * qJD(6) + t47 * t679;
t101 = -t629 * t874 + t630 * t868;
t102 = -t179 * t870 - t630 * t869;
t675 = t102 * qJD(6) + t101 * t679;
t674 = -t857 / 0.2e1;
t673 = t856 / 0.2e1;
t672 = -t853 / 0.2e1;
t670 = t848 / 0.2e1;
t668 = -t540 * pkin(5) / 0.2e1;
t667 = -t845 / 0.2e1;
t665 = t934 * t912;
t664 = t934 * t913;
t661 = t538 * t712;
t660 = t501 * t710;
t418 = t501 * t682;
t659 = t538 * t710;
t658 = t542 * t683;
t653 = t481 * t864;
t652 = t482 * t864;
t647 = -t785 / 0.2e1;
t646 = -t783 / 0.2e1;
t644 = -t763 / 0.2e1;
t643 = -t761 / 0.2e1;
t635 = t395 / 0.2e1 - t412 / 0.2e1;
t634 = pkin(4) * t679;
t633 = pkin(5) * t947;
t625 = t679 * t506;
t623 = t542 * t670;
t622 = qJD(1) * t525 + qJD(2);
t621 = -qJD(4) - t684;
t620 = -qJD(6) - t684;
t619 = t847 / 0.2e1 + t97 / 0.2e1;
t613 = qJD(3) * t632;
t609 = qJD(5) - t621;
t303 = t441 * t889;
t217 = t443 * t629 + t303;
t544 = t294 * t870 - t630 * t867 - t951;
t587 = t647 + t831 / 0.2e1;
t546 = t587 - t891;
t3 = t653 + t544 + t546;
t605 = -t3 * qJD(1) + t217 * qJD(3);
t304 = t441 * t629;
t218 = t443 * t889 - t304;
t543 = -t179 * t867 + t294 * t869 + t948;
t586 = -t832 / 0.2e1 + t644;
t545 = t586 - t890;
t4 = t652 + t543 + t545;
t604 = -t4 * qJD(1) + t218 * qJD(3);
t23 = qJD(3) * t47 + t926;
t603 = qJD(3) * t61 + t926;
t11 = t948 + (-t179 * t863 + t365 * t868 + t536 * t864) * pkin(5) + t545;
t233 = -t856 * t889 + t304;
t598 = -t11 * qJD(1) - t233 * qJD(3);
t569 = t259 * t501;
t12 = t569 / 0.2e1 + (t540 * t503 / 0.2e1 - t630 * t863 + t629 * t871) * pkin(5) + t546;
t232 = -t629 * t856 - t303;
t597 = -t12 * qJD(1) - t232 * qJD(3);
t35 = qJD(1) * t47 + t916;
t596 = qJD(1) * t61 + t916;
t375 = t504 * t857 + t506 * t527;
t52 = (-t298 * t861 + (t504 * t859 + t860) * t503) * pkin(4) + t66;
t595 = -t52 * qJD(1) + t375 * qJD(3);
t376 = -t504 * t527 + t506 * t857;
t50 = (t365 * t861 + (-t537 / 0.2e1 + t506 * t859) * t503) * pkin(4) + t67;
t594 = t50 * qJD(1) - t376 * qJD(3);
t430 = t640 + t508 / 0.2e1;
t56 = t641 + t251 / 0.2e1 + (t876 + t618) * t537;
t593 = t56 * qJD(1) + t430 * qJD(3);
t592 = t621 * t542;
t200 = t298 ^ 2 - t365 ^ 2;
t76 = qJD(1) * t200 + qJD(3) * t145;
t357 = t504 ^ 2 - t506 ^ 2;
t116 = qJD(1) * t145 + qJD(3) * t357;
t182 = qJD(3) * t393;
t590 = t501 * t668 + t842;
t589 = t850 / 0.2e1 + t849 / 0.2e1;
t548 = t589 * t538 + t412 / 0.2e1;
t159 = t548 + t635;
t585 = pkin(3) * t712 - t159 * qJD(1);
t566 = t589 * t542;
t161 = -t396 / 0.2e1 - t566;
t584 = pkin(3) * t538 * qJD(3) - t161 * qJD(1);
t583 = -t784 / 0.2e1 + t643;
t582 = t646 + t762 / 0.2e1;
t550 = t891 + t951;
t25 = -t550 + t587;
t578 = qJD(1) * t25 - t715 * t889;
t551 = t890 - t948;
t26 = -t551 + t586;
t577 = qJD(1) * t26 + t629 * t715;
t576 = qJD(1) * t67 + t504 * t713;
t575 = qJD(1) * t66 - t506 * t713;
t574 = t503 * t592;
t54 = qJD(3) * t101 - t901;
t72 = qJD(1) * t101 - t903;
t573 = -qJD(3) * t102 + t901;
t572 = -qJD(1) * t102 + t903;
t125 = -qJD(3) * t213 - t298 * t718;
t148 = -qJD(1) * t213 + t504 * t714;
t382 = (t532 / 0.2e1 - t533 / 0.2e1) * t503;
t571 = -t382 * qJD(1) + t661;
t567 = qJD(6) * t630 + t892;
t563 = qJD(1) * t493 * t766 + t382 * qJD(3);
t392 = t522 * t493;
t562 = t392 * qJD(1) + t613;
t434 = t516 / 0.2e1 + (t667 + pkin(5) / 0.2e1) * t540;
t552 = t564 / 0.2e1;
t9 = t552 - t786 / 0.2e1 + t590;
t559 = t9 * qJD(1) - t434 * qJD(4);
t19 = t646 + t791 / 0.2e1 + t553 + t614;
t558 = qJD(1) * t19 + qJD(4) * t489;
t554 = -t536 * t639 + t643;
t16 = -t793 / 0.2e1 + t85 / 0.2e1 + t554;
t556 = qJD(1) * t16 - qJD(4) * t481;
t20 = t790 / 0.2e1 + t552 + t554;
t555 = qJD(1) * t20 + qJD(4) * t490;
t484 = t490 * qJD(5);
t483 = t489 * qJD(5);
t476 = t488 * qJD(3);
t466 = t542 * t682;
t457 = t482 * qJD(6);
t456 = t481 * qJD(6);
t415 = t668 + t523 - t516 / 0.2e1 + t540 * t667;
t414 = t669 - t617 / 0.2e1 + (-t754 - t773 / 0.2e1) * pkin(4);
t372 = t383 * qJD(4);
t371 = t382 * qJD(4);
t364 = 0.2e1 * t640 + t769;
t348 = -t694 - t711;
t306 = 0.2e1 * t642 + t774;
t290 = t649 + t731;
t289 = t650 + t730;
t285 = t291 * qJD(2);
t280 = t290 * qJD(2);
t278 = t288 * qJD(2);
t275 = t286 * qJD(2);
t236 = -t625 - t276;
t235 = -t279 + t896;
t234 = -t277 - t896;
t185 = qJD(3) * t287 + t298 * t684;
t184 = qJD(3) * t286 - t365 * t684;
t162 = t798 + t396 / 0.2e1 - t566;
t160 = t548 - t635;
t140 = qJD(6) * t393 + t924;
t115 = t289 * qJD(3) - t298 * t609;
t114 = t291 * qJD(3) + t365 * t609;
t113 = -t393 * t679 + t924;
t82 = t306 * qJD(6) - t679 * t889 - t924;
t78 = t306 * t679 - t688 - t924;
t59 = pkin(4) * t651 + t771 + t877 - t753 / 0.2e1;
t57 = 0.2e1 * t641 - t772 / 0.2e1 + (t671 + t876) * t537;
t53 = -t298 * t674 + t504 * t623 + t541 * t670 + t68;
t51 = t365 * t674 + t506 * t623 + t846 * t864 + t69;
t40 = t937 + (qJD(6) + t609) * t630;
t39 = t179 * t609 + t915 + t925;
t28 = t550 + t587;
t27 = t551 + t586;
t22 = -t790 / 0.2e1 + t583 + t840;
t21 = -t791 / 0.2e1 + t582 + t841;
t18 = -t792 / 0.2e1 + t582 + t843;
t17 = t793 / 0.2e1 + t583 + t842;
t14 = t570 / 0.2e1 - t179 * t673 + t889 * t672 + t644 - t619 * t536 + t890;
t13 = -t569 / 0.2e1 - t630 * t673 + t629 * t672 + t647 + t619 * t540 + t891;
t10 = t590 + t840;
t8 = t591 + t841;
t6 = t652 - t543 + t586 + t890;
t5 = t653 - t544 + t587 + t891;
t24 = [0, 0, 0, 0, 0, t519 * qJD(2), t509 * qJD(2), -t418, -t677 * qJD(3), 0, 0, 0, t525 * t682, -t525 * t487, -t418 * t533 - t493 * t659, -t392 * qJD(4) + t501 * t613, -t501 * t503 * t711 + qJD(3) * t319, -t317 * qJD(3) - t503 * t660, t418, -qJD(2) * t318 + qJD(3) * t87 + qJD(4) * t177, -qJD(2) * t394 + qJD(3) * t88 + qJD(4) * t176 (-qJD(3) * t369 + t298 * t679) * t365, qJD(3) * t146 + t200 * t679, t211 * qJD(3) - t298 * t897, t209 * qJD(3) + t365 * t897, t418, qJD(2) * t210 + qJD(3) * t38 + qJD(4) * t65 + qJD(5) * t75, qJD(2) * t212 + qJD(3) * t37 + qJD(4) * t64 + qJD(5) * t74 -(qJD(3) * t268 + t567) * t179, qJD(3) * t63 + t911 * t934, t120 * qJD(3) + t501 * t567, t118 * qJD(3) + (t915 + t921) * t501, t418, qJD(2) * t119 + qJD(3) * t1 + qJD(4) * t31 + qJD(5) * t29 + qJD(6) * t34, qJD(2) * t121 + qJD(3) * t2 + qJD(4) * t32 + qJD(5) * t30 + qJD(6) * t33; 0, 0, 0, 0, 0, t680, t681, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t720, -t690, 0, 0, 0, 0, 0, t724 + t741, t290 * t679 + t723, 0, 0, 0, 0, 0, t664 + t728, t930 * t934 + t727; 0, 0, 0, 0, 0, 0, 0, -t419, -t695, -t487, -t682, 0, -qJD(3) * t428 + t525 * t683, qJD(3) * t427 - t525 * t684, -t371 + (-t533 * t683 - t661) * t501, t501 * t561 - 0.2e1 * t503 * t659, t538 * t682 + t719, t466 - t721, t732, t747 + (t538 * t616 - t750) * qJD(3) + t162 * qJD(4), t746 + (t542 * t616 + t767) * qJD(3) + t160 * qJD(4) (t714 - t718) * t369 + t744, t707 + (t366 * t506 - t369 * t504) * qJD(3) + t745, t289 * t679 + t506 * t682 + t701, -t504 * t682 + t725 + t741, t917, t803 + (t323 * t504 - t366 * t527 - t503 * t628) * qJD(3) + t53 * qJD(4) + t68 * qJD(5), t804 + (t323 * t506 + t369 * t527 - t432 * t503) * qJD(3) + t51 * qJD(4) + t69 * qJD(5) (t716 - t722) * t268 + t675, t748 + (-t264 * t889 - t268 * t629) * qJD(3) + t678, t682 * t889 + t927 * t934 + t708, -t629 * t682 + t664 + t729, t883, t835 + (t262 * t629 + t264 * t441 - t503 * t950) * qJD(3) + t5 * qJD(4) + t13 * qJD(5) + t28 * qJD(6), t834 + (-t259 * t503 + t262 * t889 + t268 * t441) * qJD(3) + t6 * qJD(4) + t14 * qJD(5) + t27 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t563, -t562, t621 * t385, t574, t476, qJD(3) * t162 - qJD(4) * t274 + t726, qJD(3) * t160 + qJD(4) * t273 + t702, -t125, t76, t115, t114, t476, qJD(3) * t53 + qJD(4) * t126 + qJD(5) * t57 + t285 + t823, qJD(3) * t51 - qJD(4) * t127 + qJD(5) * t59 + t280 + t824, t54, t23, t40, t39, t476, qJD(3) * t5 + qJD(4) * t48 + qJD(5) * t21 + qJD(6) * t18 + t155 + t828, qJD(3) * t6 - qJD(4) * t49 + qJD(5) * t22 + qJD(6) * t17 + t156 + t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, t76, t115, t114, t476, qJD(3) * t68 + qJD(4) * t57 - qJD(5) * t123 + t285 + t821, qJD(3) * t69 + qJD(4) * t59 + qJD(5) * t122 + t280 + t822, t54, t23, t40, t39, t476, qJD(3) * t13 + qJD(4) * t21 + qJD(5) * t43 + qJD(6) * t8 + t155 + t830, qJD(3) * t14 + qJD(4) * t22 - qJD(5) * t44 + qJD(6) * t10 + t156 + t829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t573, t603, -t620 * t630 + t892 + t937, -t179 * t620 + t921 + t925, t476, qJD(3) * t28 + qJD(4) * t18 + qJD(5) * t8 - qJD(6) * t42 + t155 + t825, qJD(3) * t27 + qJD(4) * t17 + qJD(5) * t10 + qJD(6) * t41 + t156 + t826; 0, 0, 0, 0, 0, -t680, -t681, 0, 0, 0, 0, 0, t682, -t487, 0, 0, 0, 0, 0, -t372 + t466 + t720, -t385 * qJD(3) - t660 + t690, 0, 0, 0, 0, 0, -qJD(3) * t296 - t724 - t742, -qJD(3) * t298 - t288 * t679 - t723, 0, 0, 0, 0, 0, -qJD(3) * t240 - t665 - t728, -qJD(3) * t242 - t928 * t934 - t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t683, -t684, 0, 0, 0, 0, 0, t658, -t693, 0, 0, 0, 0, 0, -t698, -t697, 0, 0, 0, 0, 0, -t700, -t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, t592, 0, 0, 0, 0, 0, t236, t235, 0, 0, 0, 0, 0, t82, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, t235, 0, 0, 0, 0, 0, t82, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t84; 0, 0, 0, 0, 0, 0, 0, t419, t695, 0, 0, 0, -t622 * t503, t622 * t501, t419 * t533 - t371, 0.2e1 * t538 * t574, qJD(4) * t388 - t719, -t372 + t721, -t732, t161 * qJD(4) - t503 * t717 - t747, qJD(2) * t385 + qJD(4) * t159 - t746, t369 * t718 + t744, -t707 + t745, -t287 * t679 - t701, -t725 - t742, -t917, qJD(2) * t296 - qJD(4) * t52 - qJD(5) * t66 - t803, qJD(2) * t298 - qJD(4) * t50 - qJD(5) * t67 - t804, t268 * t722 + t675, t678 - t748, -t929 * t934 - t708, -t665 - t729, -t883, qJD(2) * t240 - qJD(4) * t3 - qJD(5) * t12 - qJD(6) * t25 - t835, qJD(2) * t242 - qJD(4) * t4 - qJD(5) * t11 - qJD(6) * t26 - t834; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t683, t684, 0, 0, 0, 0, 0, -t658, t693, 0, 0, 0, 0, 0, t698, t697, 0, 0, 0, 0, 0, t700, t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t659, t522 * qJD(4), 0, 0, 0, -pkin(3) * t711, -pkin(3) * t710, -t504 * t625, t679 * t357, 0, 0, 0, qJD(4) * t375 + t506 * t709, qJD(4) * t376 - t504 * t709, -t918 * t889, t934 * t898, 0, 0, 0, qJD(4) * t217 - qJD(5) * t232 + t441 * t688, qJD(4) * t218 - qJD(5) * t233 - t441 * t884; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t571, -t561, t692 + t710, t348, -t687, -pkin(8) * t710 - t584, pkin(8) * t711 - t585, -t148, t116, t234, t236, -t687, -qJD(4) * t432 + t364 * qJD(5) + t595, -t594 + t879, t72, t35, t83, t82, -t687, t605 - t953, t604 + t952; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t116, t234, t236, -t687, qJD(4) * t364 - qJD(5) * t432 - t575, -t576 + t879, t72, t35, t83, t82, -t687, t597 - t953, t598 + t952; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t572, t596, t83, t78, -t687, -t578 - t953, -t577 + t952; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t563, t562, -qJD(3) * t388 + t419 * t538, t383 * qJD(3) + t501 * t658, t476, qJD(2) * t383 - qJD(3) * t161 - t726, -t159 * qJD(3) + t501 * t717 - t702, t125, -t76, t185, t184, t476, qJD(3) * t52 + qJD(5) * t56 + t275 - t823, qJD(3) * t50 + qJD(5) * t58 + t278 - t824, -t54, -t23, t71, t70, t476, qJD(3) * t3 - qJD(5) * t19 - qJD(6) * t15 + t149 - t828, qJD(3) * t4 - qJD(5) * t20 - qJD(6) * t16 + t150 - t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t694, t542 * t684, 0, 0, 0, 0, 0, t276, t279, 0, 0, 0, 0, 0, t140, t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t571, t561, -t692, t694, t687, t584, t585, t148, -t116, t277, t276, t687, qJD(5) * t430 - t595, t594, -t72, -t35, t935, t140, t687, -t605, -t604; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t537 * t838, -t541 * t838, 0, 0, 0, 0, 0, -t483 - t457, -t484 + t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t537 * t634 + t593, -t541 * t634 + t749, 0, 0, 0, 0, 0, qJD(6) * t414 - t483 - t558, qJD(6) * t415 - t484 - t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, 0, qJD(5) * t414 - t457 - t886, qJD(5) * t415 + t456 - t556; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t76, t185, t184, t476, qJD(3) * t66 - qJD(4) * t56 + t275 - t821, qJD(3) * t67 - qJD(4) * t58 + t278 - t822, -t54, -t23, t71, t70, t476, qJD(3) * t12 + qJD(4) * t19 + qJD(6) * t7 + t149 - t830, qJD(3) * t11 + qJD(4) * t20 + qJD(6) * t9 + t150 - t829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, t279, 0, 0, 0, 0, 0, t140, t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t116, t277, t276, t687, -qJD(4) * t430 + t575, t576, -t72, -t35, t935, t140, t687, -t597, -t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t537 * t839 - t593, t541 * t839 - t749, 0, 0, 0, 0, 0, -qJD(6) * t433 + t558, -qJD(6) * t434 + t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t536 * t836, -t540 * t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, 0, -t536 * t633 + t885, -t540 * t633 + t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t573, -t603, t71, t70, t476, qJD(3) * t25 + qJD(4) * t15 - qJD(5) * t7 + t149 - t825, qJD(3) * t26 + qJD(4) * t16 - qJD(5) * t9 + t150 - t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t936; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t572, -t596, t935, t113, t687, t578, t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, 0, qJD(5) * t433 + t886, qJD(5) * t434 + t556; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, 0, t536 * t837 - t885, t540 * t837 - t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t24;
