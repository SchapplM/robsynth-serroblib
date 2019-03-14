% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x33]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:23
% EndTime: 2019-03-09 13:26:13
% DurationCPUTime: 32.18s
% Computational Cost: add. (21839->709), mult. (44317->926), div. (0->0), fcn. (52819->10), ass. (0->584)
t969 = qJD(5) + qJD(6);
t978 = qJD(4) + t969;
t549 = cos(qJ(4));
t542 = sin(pkin(11));
t534 = pkin(2) * t542 + pkin(8);
t869 = pkin(9) + t534;
t507 = t869 * t549;
t548 = cos(qJ(5));
t482 = t548 * t507;
t545 = sin(qJ(4));
t506 = t869 * t545;
t544 = sin(qJ(5));
t787 = t544 * t506;
t389 = t482 - t787;
t770 = t548 * t549;
t517 = t544 * t545 - t770;
t332 = -pkin(10) * t517 + t389;
t543 = sin(qJ(6));
t547 = cos(qJ(6));
t771 = t548 * t545;
t786 = t544 * t549;
t519 = t771 + t786;
t632 = t548 * t506 + t544 * t507;
t909 = -t519 * pkin(10) - t632;
t218 = t547 * t332 + t543 * t909;
t975 = t978 * t218;
t972 = t543 * t332 - t547 * t909;
t974 = t978 * t972;
t546 = sin(qJ(2));
t550 = cos(qJ(2));
t853 = cos(pkin(11));
t511 = t542 * t546 - t550 * t853;
t888 = -t511 / 0.2e1;
t973 = t218 * t888;
t887 = t511 / 0.2e1;
t970 = t887 * t972;
t631 = t547 * t517 + t543 * t519;
t687 = qJD(4) + qJD(5);
t956 = qJD(6) + t687;
t901 = t956 * t631;
t377 = t517 * t511;
t774 = t547 * t377;
t374 = t519 * t511;
t795 = t543 * t374;
t753 = t795 / 0.2e1 + t774 / 0.2e1;
t926 = t631 * t511;
t945 = t926 / 0.2e1;
t952 = t945 + t753;
t957 = t952 * qJD(1);
t81 = -t957 - t901;
t946 = -t926 / 0.2e1;
t951 = t946 - t753;
t958 = t951 * qJD(1);
t82 = -t958 + t901;
t950 = t946 + t753;
t959 = t950 * qJD(2);
t155 = t951 * qJD(3);
t953 = t945 - t753;
t161 = t953 * qJD(3);
t638 = t853 * t546;
t805 = t542 * t550;
t513 = t638 + t805;
t303 = t519 * t513;
t399 = t545 * t513;
t373 = t399 * t544 - t513 * t770;
t633 = -t303 * t547 + t543 * t373;
t732 = qJD(1) * t511;
t71 = qJD(2) * t952 - t633 * t732;
t885 = t517 / 0.2e1;
t178 = t303 * t543 + t547 * t373;
t935 = -t178 ^ 2 + t633 ^ 2;
t949 = qJD(1) * t935;
t366 = t547 * t374;
t793 = t543 * t377;
t618 = t366 / 0.2e1 - t793 / 0.2e1;
t502 = t547 * t519;
t792 = t543 * t517;
t910 = t502 - t792;
t922 = t910 * t511;
t937 = -t922 / 0.2e1 + t618;
t948 = qJD(2) * t937;
t936 = t922 / 0.2e1 + t618;
t947 = t936 * qJD(1);
t154 = t936 * qJD(3);
t160 = t937 * qJD(3);
t70 = qJD(2) * t936 - t178 * t732;
t868 = -qJ(3) - pkin(7);
t523 = t868 * t546;
t524 = t868 * t550;
t438 = -t853 * t523 - t524 * t542;
t351 = pkin(4) * t399 + t438;
t272 = pkin(5) * t303 + t351;
t535 = -pkin(2) * t853 - pkin(3);
t522 = -t549 * pkin(4) + t535;
t441 = t517 * pkin(5) + t522;
t891 = t910 / 0.2e1;
t896 = -t178 / 0.2e1;
t912 = t272 * t891 + t441 * t896;
t893 = -t631 / 0.2e1;
t897 = t633 / 0.2e1;
t911 = t272 * t893 + t441 * t897;
t929 = t438 * t549;
t944 = t929 / 0.2e1;
t842 = t178 * t272;
t943 = t178 * t687;
t503 = t638 / 0.2e1 + t805 / 0.2e1;
t731 = qJD(1) * t513;
t429 = t511 * t731;
t749 = t503 * qJD(4) + t429;
t941 = t503 * qJD(5) + t749;
t919 = t631 ^ 2 - t910 ^ 2;
t940 = qJD(2) * t919;
t939 = t178 * qJD(6);
t933 = -t522 / 0.2e1;
t938 = -t303 * t933 + t351 * t885 + t632 * t888;
t883 = t519 / 0.2e1;
t930 = t178 * t910;
t728 = qJD(2) * t910;
t925 = t631 * t728;
t839 = t633 * t272;
t739 = qJD(1) * t178;
t923 = t633 * t739;
t907 = t542 * t523 - t853 * t524;
t818 = t907 * t545;
t817 = t907 * t549;
t758 = t351 * t883 + t373 * t933;
t920 = t389 * t888 + t758;
t871 = t546 * pkin(2);
t402 = pkin(3) * t513 + pkin(8) * t511 + t871;
t380 = t549 * t402;
t401 = t549 * t511;
t784 = t438 * t545;
t874 = t513 * pkin(4);
t236 = pkin(9) * t401 + t380 + t784 + t874;
t233 = t544 * t236;
t396 = t545 * t511;
t379 = t545 * t402;
t751 = -t379 + t929;
t263 = pkin(9) * t396 - t751;
t258 = t548 * t263;
t642 = -t233 / 0.2e1 - t258 / 0.2e1;
t63 = t642 + t938;
t64 = t642 - t938;
t918 = t511 * t687;
t917 = t517 * t687;
t913 = t633 * t687;
t509 = t511 ^ 2;
t510 = t513 ^ 2;
t908 = -t510 - t509;
t670 = -pkin(2) * t550 - pkin(1);
t394 = t511 * pkin(3) - t513 * pkin(8) + t670;
t287 = -t549 * t394 + t818;
t809 = t513 * t549;
t259 = -pkin(9) * t809 - t287;
t772 = t548 * t259;
t288 = t394 * t545 + t817;
t260 = -pkin(9) * t399 + t288;
t789 = t544 * t260;
t139 = t772 - t789;
t876 = t373 * pkin(10);
t112 = t139 + t876;
t875 = t511 * pkin(4);
t235 = t259 + t875;
t257 = t548 * t260;
t135 = t544 * t235 + t257;
t878 = pkin(10) * t303;
t105 = t135 - t878;
t790 = t544 * t259;
t138 = -t257 - t790;
t111 = t138 + t878;
t643 = t111 / 0.2e1 + t105 / 0.2e1;
t621 = t643 * t547;
t870 = t548 * pkin(4);
t669 = pkin(5) + t870;
t623 = t543 * t669;
t773 = t547 * t544;
t498 = pkin(4) * t773 + t623;
t812 = t498 * t511;
t230 = t548 * t235;
t134 = -t230 + t789;
t596 = -t134 + t876;
t90 = t511 * pkin(5) + t596;
t15 = t812 / 0.2e1 + t621 + (-t112 / 0.2e1 + t90 / 0.2e1) * t543;
t647 = -t502 / 0.2e1;
t405 = t647 + t502 / 0.2e1;
t693 = t405 * qJD(3);
t906 = qJD(1) * t15 + qJD(4) * t498 - t693;
t574 = t543 * t596;
t560 = t574 / 0.2e1;
t446 = t543 * pkin(5);
t676 = -t446 / 0.2e1;
t783 = t547 * t105;
t650 = -t783 / 0.2e1;
t856 = t543 * t90;
t867 = -t856 / 0.2e1 + t650;
t598 = t511 * t676 + t867;
t11 = t560 + t783 / 0.2e1 + t598;
t905 = t11 * qJD(1) - t446 * qJD(4) + t693;
t903 = qJD(6) * t503 + t941;
t234 = t548 * t236;
t788 = t544 * t263;
t762 = t234 / 0.2e1 - t788 / 0.2e1;
t540 = t545 ^ 2;
t541 = t549 ^ 2;
t530 = t541 - t540;
t635 = 0.2e1 * t549 * t399;
t570 = qJD(1) * t635 - qJD(2) * t530;
t900 = t687 * t632;
t899 = -t230 / 0.2e1;
t898 = -t235 / 0.2e1;
t646 = -t257 / 0.2e1;
t892 = -t910 / 0.2e1;
t890 = t438 / 0.2e1;
t879 = pkin(5) * t519;
t881 = pkin(4) * t545;
t457 = t879 + t881;
t889 = -t457 / 0.2e1;
t645 = -t482 / 0.2e1;
t886 = -t513 / 0.2e1;
t884 = -t519 / 0.2e1;
t882 = -t548 / 0.2e1;
t880 = pkin(5) * t373;
t873 = t513 * pkin(5);
t872 = t544 * pkin(4);
t804 = t543 * t105;
t654 = t804 / 0.2e1;
t85 = t547 * t90;
t866 = -t85 / 0.2e1 + t654;
t865 = -t574 / 0.2e1 + t650;
t573 = t547 * t596;
t864 = t654 - t573 / 0.2e1;
t863 = pkin(4) * qJD(4);
t862 = pkin(4) * qJD(5);
t861 = pkin(5) * qJD(5);
t860 = pkin(5) * qJD(6);
t859 = qJD(2) * pkin(2);
t265 = -t366 + t793;
t350 = -pkin(4) * t396 + t907;
t271 = -pkin(5) * t374 + t350;
t41 = -t85 + t804;
t761 = t258 + t233;
t106 = pkin(10) * t374 + t761;
t803 = t543 * t106;
t634 = t234 - t788;
t91 = -pkin(10) * t377 + t634 + t873;
t854 = t547 * t91;
t1 = (-t803 + t854) * t511 - t41 * t513 - t271 * t633 + t272 * t265;
t858 = t1 * qJD(1);
t269 = t774 + t795;
t42 = t783 + t856;
t782 = t547 * t106;
t855 = t543 * t91;
t2 = -(t782 + t855) * t511 - t42 * t513 - t271 * t178 + t272 * t269;
t857 = t2 * qJD(1);
t43 = -t574 - t783;
t29 = t43 * t511 + t633 * t880 - t842;
t852 = qJD(1) * t29;
t44 = t573 - t804;
t30 = t178 * t880 - t44 * t511 + t839;
t851 = qJD(1) * t30;
t684 = pkin(4) * t809;
t297 = t684 - t880;
t781 = t547 * t111;
t801 = t543 * t112;
t48 = t781 - t801;
t31 = -t297 * t633 + t48 * t511 - t842;
t850 = qJD(1) * t31;
t780 = t547 * t112;
t802 = t543 * t111;
t49 = t780 + t802;
t32 = -t178 * t297 - t49 * t511 + t839;
t849 = qJD(1) * t32;
t33 = t41 * t511 + t839;
t848 = qJD(1) * t33;
t34 = -t42 * t511 - t842;
t847 = qJD(1) * t34;
t834 = t351 * t373;
t68 = t138 * t511 + t303 * t684 - t834;
t846 = qJD(1) * t68;
t833 = t351 * t303;
t69 = -t139 * t511 - t373 * t684 - t833;
t845 = qJD(1) * t69;
t74 = t134 * t511 - t833;
t844 = qJD(1) * t74;
t75 = -t135 * t511 - t834;
t843 = qJD(1) * t75;
t841 = t265 * t511;
t840 = t633 * t513;
t838 = t269 * t511;
t837 = t178 * t513;
t37 = -t134 * t513 + t303 * t350 - t351 * t374 + t511 * t634;
t830 = t37 * qJD(1);
t829 = t373 * t513;
t828 = t374 * t511;
t827 = t303 * t513;
t826 = t377 * t511;
t38 = -t135 * t513 - t350 * t373 + t351 * t377 - t511 * t761;
t825 = t38 * qJD(1);
t525 = t547 * t669;
t532 = t543 * t872;
t497 = t532 - t525;
t813 = t497 * t511;
t791 = t543 * t548;
t504 = (t773 + t791) * pkin(4);
t811 = t504 * t511;
t505 = t547 * t870 - t532;
t810 = t505 * t511;
t785 = t545 * t373;
t769 = t549 * t519;
t678 = -t875 / 0.2e1;
t624 = t678 + t259 / 0.2e1;
t58 = t548 * t624 + t899;
t768 = t58 * qJD(1);
t67 = t178 * t265 + t269 * t633;
t767 = t67 * qJD(1);
t88 = (-t287 + t818) * t513 + t380 * t511;
t766 = t88 * qJD(1);
t89 = (-t288 + t817) * t513 + (t751 - t929) * t511;
t765 = t89 * qJD(1);
t146 = t303 * t517 + t519 * t373;
t764 = t687 * t146;
t223 = -t303 * t883 + t373 * t885;
t763 = t687 * t223;
t588 = t786 / 0.2e1 + t771 / 0.2e1;
t289 = (t883 + t588) * t511;
t760 = t687 * t289;
t294 = (t588 + t884) * t511;
t759 = t687 * t294;
t673 = t544 * t396;
t747 = -t673 / 0.2e1 + t548 * t401 / 0.2e1;
t657 = t511 * t882;
t746 = t673 / 0.2e1 + t549 * t657;
t130 = t840 - t841;
t745 = qJD(1) * t130;
t131 = -t840 - t841;
t744 = qJD(1) * t131;
t179 = t287 * t511 - t399 * t438;
t743 = qJD(1) * t179;
t180 = -t288 * t511 + t438 * t809;
t742 = qJD(1) * t180;
t219 = -t827 + t828;
t741 = qJD(1) * t219;
t222 = t827 + t828;
t740 = qJD(1) * t222;
t273 = t438 * t513 - t511 * t907;
t738 = qJD(1) * t273;
t685 = t510 - t509;
t345 = t685 * t545;
t737 = qJD(1) * t345;
t346 = t908 * t545;
t736 = qJD(1) * t346;
t347 = t685 * t549;
t735 = qJD(1) * t347;
t734 = qJD(1) * t373;
t406 = t908 * t549;
t733 = qJD(1) * t406;
t730 = qJD(1) * t549;
t729 = qJD(1) * t550;
t727 = qJD(2) * t441;
t726 = qJD(2) * t513;
t725 = qJD(2) * t519;
t724 = qJD(2) * t522;
t723 = qJD(2) * t545;
t722 = qJD(2) * t546;
t721 = qJD(2) * t549;
t720 = qJD(2) * t550;
t719 = qJD(3) * t549;
t718 = qJD(4) * t545;
t717 = qJD(4) * t549;
t716 = qJD(5) * t522;
t713 = qJD(6) * t441;
t132 = -t837 + t838;
t712 = t132 * qJD(1);
t133 = -t837 - t838;
t711 = t133 * qJD(1);
t147 = -t303 * t377 - t373 * t374;
t709 = t147 * qJD(1);
t203 = t670 * t871;
t704 = t203 * qJD(1);
t220 = t826 - t829;
t703 = t220 * qJD(1);
t221 = -t826 - t829;
t702 = t221 * qJD(1);
t243 = t513 * t631;
t701 = t243 * qJD(1);
t246 = t910 * t513;
t700 = t246 * qJD(1);
t277 = t289 * qJD(1);
t655 = t377 / 0.2e1;
t290 = t655 + t746;
t278 = t290 * qJD(1);
t656 = -t377 / 0.2e1;
t291 = t656 + t747;
t280 = t291 * qJD(1);
t300 = t513 * t517;
t699 = t300 * qJD(1);
t698 = t303 * qJD(1);
t569 = t542 * t888 + t853 * t886;
t364 = (-t546 / 0.2e1 + t569) * pkin(2);
t697 = t364 * qJD(1);
t696 = t396 * qJD(1);
t695 = t399 * qJD(1);
t694 = t401 * qJD(1);
t692 = t908 * qJD(1);
t691 = t503 * qJD(1);
t531 = -t546 ^ 2 + t550 ^ 2;
t688 = t531 * qJD(1);
t47 = -t631 * t897 - t178 * t892 + t633 * t893 + t930 / 0.2e1;
t61 = -t631 * t633 + t930;
t686 = t61 * qJD(6) + t47 * t687;
t103 = -t631 * t896 + t633 * t891;
t104 = -t178 * t893 - t633 * t892;
t683 = t104 * qJD(6) + t103 * t687;
t682 = pkin(1) * t546 * qJD(1);
t681 = pkin(1) * t729;
t680 = -t880 / 0.2e1;
t679 = t879 / 0.2e1;
t677 = t874 / 0.2e1;
t675 = -t547 * pkin(5) / 0.2e1;
t674 = -t870 / 0.2e1;
t672 = t956 * t936;
t671 = t956 * t937;
t668 = t545 * t721;
t667 = t511 * t717;
t428 = t511 * t726;
t666 = t545 * t717;
t665 = t546 * t729;
t664 = t513 * t730;
t659 = t497 * t886;
t658 = t498 * t886;
t653 = -t803 / 0.2e1;
t652 = -t801 / 0.2e1;
t651 = t545 * t303 / 0.2e1;
t649 = -t782 / 0.2e1;
t648 = -t780 / 0.2e1;
t644 = t549 * t885;
t639 = -t379 / 0.2e1 + t944;
t637 = pkin(4) * t687;
t636 = pkin(5) * t969;
t629 = t687 * t519;
t627 = -qJD(4) - t732;
t626 = -qJD(6) - t732;
t625 = t873 / 0.2e1 + t91 / 0.2e1;
t620 = qJD(2) * t635;
t617 = qJD(5) - t627;
t615 = -t389 * t887 + t758;
t298 = t441 * t910;
t213 = t457 * t631 + t298;
t552 = t297 * t893 - t633 * t889 - t973;
t595 = t653 + t854 / 0.2e1;
t554 = t595 - t912;
t3 = t659 + t552 + t554;
t613 = -t3 * qJD(1) + t213 * qJD(2);
t299 = t441 * t631;
t214 = t457 * t910 - t299;
t551 = -t178 * t889 + t297 * t892 - t970;
t594 = -t855 / 0.2e1 + t649;
t553 = t594 - t911;
t4 = t658 + t551 + t553;
t612 = -t4 * qJD(1) + t214 * qJD(2);
t231 = -t631 * t879 - t298;
t576 = t218 * t511;
t8 = t576 / 0.2e1 + (t547 * t513 / 0.2e1 - t633 * t884 + t373 * t631 / 0.2e1) * pkin(5) + t554;
t611 = -t8 * qJD(1) - t231 * qJD(2);
t232 = -t879 * t910 + t299;
t577 = t972 * t511;
t7 = -t577 / 0.2e1 + (-t178 * t884 + t373 * t891 + t543 * t886) * pkin(5) + t553;
t610 = -t7 * qJD(1) - t232 * qJD(2);
t27 = qJD(2) * t47 + t949;
t609 = qJD(2) * t61 + t949;
t604 = -t511 * t535 - t513 * t534;
t35 = qJD(1) * t47 + t940;
t603 = qJD(1) * t61 + t940;
t368 = t517 * t881 + t519 * t522;
t52 = (t651 + (t644 + t882) * t513) * pkin(4) + t615 - t762;
t602 = -qJD(1) * t52 - qJD(2) * t368;
t369 = -t517 * t522 + t519 * t881;
t50 = (t785 / 0.2e1 + (-t544 / 0.2e1 - t769 / 0.2e1) * t513) * pkin(4) + t63;
t601 = t50 * qJD(1) - t369 * qJD(2);
t382 = t645 + t482 / 0.2e1;
t56 = t646 + t257 / 0.2e1 + (t898 + t624) * t544;
t600 = t56 * qJD(1) + t382 * qJD(2);
t599 = t627 * t549;
t208 = t303 ^ 2 - t373 ^ 2;
t76 = qJD(1) * t208 + qJD(2) * t146;
t370 = t517 ^ 2 - t519 ^ 2;
t128 = qJD(1) * t146 + qJD(2) * t370;
t183 = qJD(2) * t405;
t597 = t511 * t675 + t866;
t593 = t534 * t887 + t535 * t886;
t592 = -t802 / 0.2e1 + t648;
t591 = t652 + t781 / 0.2e1;
t557 = t912 + t973;
t23 = -t557 + t595;
t587 = qJD(1) * t23 - t727 * t910;
t558 = t911 + t970;
t24 = -t558 + t594;
t586 = qJD(1) * t24 + t631 * t727;
t585 = qJD(1) * t63 + t517 * t724;
t62 = t762 - t920;
t584 = qJD(1) * t62 - t519 * t724;
t583 = t513 * t599;
t54 = qJD(2) * t103 - t923;
t72 = qJD(1) * t103 - t925;
t582 = -qJD(2) * t104 + t923;
t581 = -qJD(1) * t104 + t925;
t555 = t593 * t545 + t944;
t148 = t555 - t639;
t580 = -qJD(1) * t148 - t535 * t721;
t568 = t593 * t549;
t150 = -t380 / 0.2e1 - t568 + (t890 - t438 / 0.2e1) * t545;
t579 = -qJD(1) * t150 - t535 * t723;
t137 = -qJD(2) * t223 - t303 * t734;
t153 = -qJD(1) * t223 + t517 * t725;
t395 = (t540 / 0.2e1 - t541 / 0.2e1) * t513;
t578 = -qJD(1) * t395 + t668;
t575 = qJD(6) * t633 + t913;
t572 = t510 * t545 * t730 + qJD(2) * t395;
t404 = t530 * t510;
t571 = qJD(1) * t404 + t620;
t19 = t652 + t811 / 0.2e1 + t560 + t621;
t566 = qJD(1) * t19 + qJD(4) * t504;
t559 = t573 / 0.2e1;
t13 = t559 - t804 / 0.2e1 + t597;
t447 = t525 / 0.2e1 + (t674 + pkin(5) / 0.2e1) * t547;
t564 = t13 * qJD(1) - t447 * qJD(4);
t561 = -t543 * t643 + t648;
t16 = -t813 / 0.2e1 + t85 / 0.2e1 + t561;
t563 = qJD(1) * t16 - qJD(4) * t497;
t20 = t810 / 0.2e1 + t559 + t561;
t562 = qJD(1) * t20 + qJD(4) * t505;
t500 = t505 * qJD(5);
t499 = t504 * qJD(5);
t492 = t503 * qJD(2);
t480 = t513 * t721;
t471 = t498 * qJD(6);
t470 = t497 * qJD(6);
t424 = t675 + t532 - t525 / 0.2e1 + t547 * t674;
t423 = t676 - t623 / 0.2e1 + (-t773 - t791 / 0.2e1) * pkin(4);
t384 = t396 * qJD(4);
t383 = t395 * qJD(4);
t363 = t871 / 0.2e1 + t569 * pkin(2);
t362 = -t696 - t718;
t334 = 0.2e1 * t647 + t792;
t305 = 0.2e1 * t645 + t787;
t293 = t655 + t747;
t292 = t656 + t746;
t286 = t294 * qJD(3);
t281 = t293 * qJD(3);
t279 = t291 * qJD(3);
t276 = t289 * qJD(3);
t241 = -t629 - t277;
t240 = -t280 + t917;
t239 = -t278 - t917;
t186 = qJD(2) * t290 + t303 * t732;
t185 = qJD(2) * t289 - t373 * t732;
t151 = t545 * t890 + t784 / 0.2e1 + t380 / 0.2e1 - t568;
t149 = t555 + t639;
t140 = qJD(6) * t405 + t947;
t127 = t292 * qJD(2) - t303 * t617;
t126 = t294 * qJD(2) + t373 * t617;
t125 = -t405 * t687 + t947;
t80 = qJD(6) * t334 - t687 * t910 - t947;
t78 = -qJD(6) * t910 + t334 * t687 - t947;
t65 = t762 + t920;
t59 = pkin(4) * t657 + t789 + t899 - t772 / 0.2e1;
t57 = 0.2e1 * t646 - t790 / 0.2e1 + (t678 + t898) * t544;
t53 = -pkin(4) * t785 / 0.2e1 + t677 * t769 + t872 * t886 + t64;
t51 = t548 * t677 + (t513 * t644 + t651) * pkin(4) + t615 + t762;
t40 = t959 + (qJD(6) + t617) * t633;
t39 = t178 * t617 + t939 + t948;
t26 = t557 + t595;
t25 = t558 + t594;
t22 = -t810 / 0.2e1 + t592 + t864;
t21 = -t811 / 0.2e1 + t591 + t865;
t18 = -t812 / 0.2e1 + t591 + t867;
t17 = t813 / 0.2e1 + t592 + t866;
t14 = t597 + t864;
t12 = t598 + t865;
t10 = t577 / 0.2e1 - t178 * t679 + t910 * t680 + t649 - t625 * t543 + t911;
t9 = -t576 / 0.2e1 - t633 * t679 + t631 * t680 + t653 + t625 * t547 + t912;
t6 = t658 - t551 + t594 + t911;
t5 = t659 - t552 + t595 + t912;
t28 = [0, 0, 0, t546 * t720, t531 * qJD(2), 0, 0, 0, -pkin(1) * t722, -pkin(1) * t720, -qJD(3) * t908, qJD(2) * t203 + qJD(3) * t273, -t428 * t541 - t510 * t666, -qJD(4) * t404 + t511 * t620, -t511 * t513 * t718 + qJD(2) * t347, -qJD(2) * t345 - t513 * t667, t428, qJD(2) * t88 - qJD(3) * t346 + qJD(4) * t180, qJD(2) * t89 - qJD(3) * t406 + qJD(4) * t179 (-qJD(2) * t377 + t303 * t687) * t373, qJD(2) * t147 + t208 * t687, t220 * qJD(2) - t303 * t918, t219 * qJD(2) + t373 * t918, t428, qJD(2) * t37 + qJD(3) * t222 + qJD(4) * t68 + qJD(5) * t75, qJD(2) * t38 + qJD(3) * t221 + qJD(4) * t69 + qJD(5) * t74 -(qJD(2) * t269 + t575) * t178, qJD(2) * t67 + t935 * t956, t132 * qJD(2) + t511 * t575, t130 * qJD(2) + (t939 + t943) * t511, t428, qJD(2) * t1 + qJD(3) * t131 + qJD(4) * t31 + qJD(5) * t29 + qJD(6) * t34, qJD(2) * t2 + qJD(3) * t133 + qJD(4) * t32 + qJD(5) * t30 + qJD(6) * t33; 0, 0, 0, t665, t688, t720, -t722, 0, -pkin(7) * t720 - t682, pkin(7) * t722 - t681 (t511 * t853 - t513 * t542) * t859, t704 + (-t438 * t542 - t853 * t907) * t859 + t363 * qJD(3), -t383 + (-t541 * t731 - t668) * t511, t511 * t570 - 0.2e1 * t513 * t666, t513 * t723 + t735, t480 - t737, t749, t766 + (t545 * t604 - t817) * qJD(2) + t151 * qJD(4), t765 + (t549 * t604 + t818) * qJD(2) + t149 * qJD(4) (t725 - t734) * t377 + t763, t709 + (t374 * t519 - t377 * t517) * qJD(2) + t764, t292 * t687 + t513 * t725 + t703, -t517 * t726 + t741 + t759, t941, t830 + (t350 * t517 - t374 * t522 - t513 * t632) * qJD(2) + t51 * qJD(4) + t65 * qJD(5), t825 + (t350 * t519 + t377 * t522 - t389 * t513) * qJD(2) + t53 * qJD(4) + t64 * qJD(5) (t728 - t739) * t269 + t683, t767 + (-t265 * t910 - t269 * t631) * qJD(2) + t686, t726 * t910 + t950 * t956 + t712, -t631 * t726 + t671 + t745, t903, t858 + (t265 * t441 + t271 * t631 - t513 * t972) * qJD(2) + t5 * qJD(4) + t9 * qJD(5) + t26 * qJD(6), t857 + (-t218 * t513 + t269 * t441 + t271 * t910) * qJD(2) + t6 * qJD(4) + t10 * qJD(5) + t25 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t692, qJD(2) * t363 + t738, 0, 0, 0, 0, 0, -t736, -t733, 0, 0, 0, 0, 0, t740 + t759, t293 * t687 + t702, 0, 0, 0, 0, 0, t671 + t744, t953 * t956 + t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t572, -t571, t627 * t399, t583, t492, qJD(2) * t151 - qJD(4) * t288 + t742, qJD(2) * t149 + qJD(4) * t287 + t743, -t137, t76, t127, t126, t492, qJD(2) * t51 + qJD(4) * t138 + qJD(5) * t57 + t286 + t846, qJD(2) * t53 - qJD(4) * t139 + qJD(5) * t59 + t281 + t845, t54, t27, t40, t39, t492, qJD(2) * t5 + qJD(4) * t48 + qJD(5) * t21 + qJD(6) * t18 + t160 + t850, qJD(2) * t6 - qJD(4) * t49 + qJD(5) * t22 + qJD(6) * t17 + t161 + t849; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t76, t127, t126, t492, qJD(2) * t65 + qJD(4) * t57 - qJD(5) * t135 + t286 + t843, qJD(2) * t64 + qJD(4) * t59 + qJD(5) * t134 + t281 + t844, t54, t27, t40, t39, t492, qJD(2) * t9 + qJD(4) * t21 + qJD(5) * t43 + qJD(6) * t12 + t160 + t852, qJD(2) * t10 + qJD(4) * t22 - qJD(5) * t44 + qJD(6) * t14 + t161 + t851; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t582, t609, -t626 * t633 + t913 + t959, -t178 * t626 + t943 + t948, t492, qJD(2) * t26 + qJD(4) * t18 + qJD(5) * t12 - qJD(6) * t42 + t160 + t847, qJD(2) * t25 + qJD(4) * t17 + qJD(5) * t14 + qJD(6) * t41 + t161 + t848; 0, 0, 0, -t665, -t688, 0, 0, 0, t682, t681, 0, qJD(3) * t364 - t704, t429 * t541 - t383, 0.2e1 * t545 * t583, qJD(4) * t401 - t735, -t384 + t737, -t749, qJD(4) * t150 - t513 * t719 - t766, qJD(3) * t399 + qJD(4) * t148 - t765, t377 * t734 + t763, -t709 + t764, -t290 * t687 - t703, -t741 - t760, -t941, qJD(3) * t300 + qJD(4) * t52 - qJD(5) * t62 - t830, qJD(3) * t303 - qJD(4) * t50 - qJD(5) * t63 - t825, t269 * t739 + t683, t686 - t767, -t952 * t956 - t712, -t672 - t745, -t903, qJD(3) * t243 - qJD(4) * t3 - qJD(5) * t8 - qJD(6) * t23 - t858, qJD(3) * t246 - qJD(4) * t4 - qJD(5) * t7 - qJD(6) * t24 - t857; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666, t530 * qJD(4), 0, 0, 0, t535 * t718, t535 * t717, -t517 * t629, t687 * t370, 0, 0, 0, qJD(4) * t368 + t519 * t716, qJD(4) * t369 - t517 * t716, -t901 * t910, t956 * t919, 0, 0, 0, qJD(4) * t213 - qJD(5) * t231 + t713 * t910, qJD(4) * t214 - qJD(5) * t232 - t631 * t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t697, 0, 0, 0, 0, 0, -t664, t695, 0, 0, 0, 0, 0, t699, t698, 0, 0, 0, 0, 0, t701, t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578, -t570, t694 + t717, t362, -t691, -t534 * t717 - t579, t534 * t718 - t580, -t153, t128, t239, t241, -t691, -qJD(4) * t389 + qJD(5) * t305 - t602, -t601 + t900, t72, t35, t81, t80, -t691, t613 - t975, t612 + t974; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t128, t239, t241, -t691, qJD(4) * t305 - qJD(5) * t389 - t584, -t585 + t900, t72, t35, t81, t80, -t691, t611 - t975, t610 + t974; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t581, t603, t81, t78, -t691, -t587 - t975, -t586 + t974; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t692, -qJD(2) * t364 - t738, 0, 0, 0, 0, 0, -t384 + t480 + t736, -qJD(2) * t399 - t667 + t733, 0, 0, 0, 0, 0, -qJD(2) * t300 - t740 - t760, -qJD(2) * t303 - t291 * t687 - t702, 0, 0, 0, 0, 0, -qJD(2) * t243 - t672 - t744, -qJD(2) * t246 - t951 * t956 - t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t697, 0, 0, 0, 0, 0, t664, -t695, 0, 0, 0, 0, 0, -t699, -t698, 0, 0, 0, 0, 0, -t701, -t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t362, t599, 0, 0, 0, 0, 0, t241, t240, 0, 0, 0, 0, 0, t80, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t240, 0, 0, 0, 0, 0, t80, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t572, t571, -qJD(2) * t401 + t429 * t545, qJD(2) * t396 + t511 * t664, t492, -qJD(2) * t150 + qJD(3) * t396 - t742, -qJD(2) * t148 + t511 * t719 - t743, t137, -t76, t186, t185, t492, -qJD(2) * t52 + qJD(5) * t56 + t276 - t846, qJD(2) * t50 + qJD(5) * t58 + t279 - t845, -t54, -t27, t71, t70, t492, qJD(2) * t3 - qJD(5) * t19 - qJD(6) * t15 + t154 - t850, qJD(2) * t4 - qJD(5) * t20 - qJD(6) * t16 + t155 - t849; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t578, t570, -t694, t696, t691, t579, t580, t153, -t128, t278, t277, t691, qJD(5) * t382 + t602, t601, -t72, -t35, t957, t140, t691, -t613, -t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t696, t511 * t730, 0, 0, 0, 0, 0, t277, t280, 0, 0, 0, 0, 0, t140, t958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t544 * t862, -t548 * t862, 0, 0, 0, 0, 0, -t499 - t471, -t500 + t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t544 * t637 + t600, -t548 * t637 + t768, 0, 0, 0, 0, 0, qJD(6) * t423 - t499 - t566, qJD(6) * t424 - t500 - t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, 0, qJD(5) * t423 - t471 - t906, qJD(5) * t424 + t470 - t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t76, t186, t185, t492, qJD(2) * t62 - qJD(4) * t56 + t276 - t843, qJD(2) * t63 - qJD(4) * t58 + t279 - t844, -t54, -t27, t71, t70, t492, qJD(2) * t8 + qJD(4) * t19 + qJD(6) * t11 + t154 - t852, qJD(2) * t7 + qJD(4) * t20 + qJD(6) * t13 + t155 - t851; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t128, t278, t277, t691, -qJD(4) * t382 + t584, t585, -t72, -t35, t957, t140, t691, -t611, -t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, t280, 0, 0, 0, 0, 0, t140, t958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t544 * t863 - t600, t548 * t863 - t768, 0, 0, 0, 0, 0, -qJD(6) * t446 + t566, -qJD(6) * t447 + t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t543 * t860, -t547 * t860; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, 0, -t543 * t636 + t905, -t547 * t636 + t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t582, -t609, t71, t70, t492, qJD(2) * t23 + qJD(4) * t15 - qJD(5) * t11 + t154 - t847, qJD(2) * t24 + qJD(4) * t16 - qJD(5) * t13 + t155 - t848; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t581, -t603, t957, t125, t691, t587, t586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t958; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, 0, qJD(5) * t446 + t906, qJD(5) * t447 + t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, 0, t543 * t861 - t905, t547 * t861 - t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t28;