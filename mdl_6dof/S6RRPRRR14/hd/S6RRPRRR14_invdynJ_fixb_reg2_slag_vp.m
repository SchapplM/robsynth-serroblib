% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR14_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynJ_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:08:52
% EndTime: 2019-03-09 15:10:28
% DurationCPUTime: 56.64s
% Computational Cost: add. (107846->1203), mult. (325568->1688), div. (0->0), fcn. (284745->18), ass. (0->522)
t494 = sin(pkin(6));
t495 = cos(pkin(7));
t499 = sin(qJ(2));
t776 = sin(pkin(14));
t677 = t776 * t499;
t502 = cos(qJ(2));
t778 = cos(pkin(14));
t682 = t778 * t502;
t421 = (-t495 * t677 + t682) * t494;
t413 = qJD(1) * t421;
t498 = sin(qJ(4));
t676 = t776 * t502;
t683 = t778 * t499;
t812 = t494 * (t495 * t683 + t676);
t410 = qJD(1) * t812;
t493 = sin(pkin(7));
t777 = sin(pkin(8));
t696 = t493 * t777;
t757 = t494 * t499;
t625 = t696 * t757;
t779 = cos(pkin(8));
t564 = qJD(1) * t625 - t410 * t779;
t795 = cos(qJ(4));
t655 = t777 * t795;
t620 = qJD(4) * t655;
t640 = t779 * t778;
t614 = t493 * t640;
t695 = t493 * t776;
t801 = t498 * t695 - t795 * t614;
t745 = -qJD(4) * t801 - t413 * t795 + t495 * t620 - t498 * t564;
t679 = t777 * t495;
t386 = t498 * t679 + (t498 * t640 + t776 * t795) * t493;
t624 = t493 * t655;
t599 = t624 * t757;
t656 = t779 * t795;
t744 = qJD(1) * t599 + t386 * qJD(4) - t410 * t656 - t413 * t498;
t503 = cos(qJ(1));
t780 = cos(pkin(6));
t688 = t503 * t780;
t794 = sin(qJ(1));
t443 = t499 * t794 - t502 * t688;
t444 = t499 * t688 + t502 * t794;
t694 = t494 * t776;
t661 = t493 * t694;
t692 = t495 * t776;
t315 = -t443 * t692 + t444 * t778 - t503 * t661;
t697 = t493 * t778;
t662 = t494 * t697;
t693 = t495 * t778;
t316 = t443 * t693 + t444 * t776 + t503 * t662;
t755 = t494 * t503;
t398 = -t443 * t493 + t495 * t755;
t802 = t779 * t316 + t777 * t398;
t199 = -t315 * t795 + t498 * t802;
t271 = t316 * t777 - t398 * t779;
t497 = sin(qJ(5));
t501 = cos(qJ(5));
t141 = t199 * t501 - t271 * t497;
t496 = sin(qJ(6));
t845 = t141 * t496;
t500 = cos(qJ(6));
t844 = t141 * t500;
t717 = pkin(1) * t780;
t483 = t502 * t717;
t473 = qJD(1) * t483;
t774 = qJ(3) * t495;
t704 = -pkin(10) - t774;
t641 = t704 * t499;
t616 = t494 * t641;
t393 = qJD(1) * t616 + t473;
t482 = t499 * t717;
t756 = t494 * t502;
t565 = t704 * t756 - t482;
t394 = t565 * qJD(1);
t775 = qJ(3) * t493;
t613 = pkin(2) * t499 - t502 * t775;
t737 = qJD(1) * t494;
t423 = t613 * t737;
t262 = t778 * t393 + t394 * t692 + t423 * t695;
t206 = pkin(11) * t564 + t262;
t261 = -t393 * t776 + t394 * t693 + t423 * t697;
t713 = t499 * t737;
t667 = t493 * t713;
t716 = pkin(11) * t779;
t208 = pkin(3) * t667 - t413 * t716 + t261;
t309 = -t493 * t394 + t495 * t423;
t715 = pkin(11) * t777;
t248 = t410 * pkin(3) - t413 * t715 + t309;
t440 = pkin(2) * t692 + qJ(3) * t697;
t376 = (t614 + t679) * pkin(11) + t440;
t477 = pkin(2) * t693;
t639 = t779 * t776;
t678 = t776 * qJ(3);
t389 = t495 * pkin(3) + t477 + (-pkin(11) * t639 - t678) * t493;
t638 = t777 * t776;
t416 = (-pkin(3) * t778 - pkin(11) * t638 - pkin(2)) * t493;
t252 = -t498 * t376 + t389 * t656 + t416 * t655;
t690 = t498 * t777;
t691 = t498 * t779;
t735 = qJD(3) * t493;
t810 = -t795 * t206 - t208 * t691 - t248 * t690 + (-t498 * t639 + t778 * t795) * t735 + t252 * qJD(4);
t588 = t208 * t777 - t779 * t248 + t638 * t735;
t698 = t493 * t779;
t626 = t698 * t757;
t338 = -qJD(1) * t626 - t410 * t777;
t843 = pkin(12) * t338 + t810;
t842 = t744 * pkin(4) - t745 * pkin(12) + t588;
t253 = t795 * t376 + (t779 * t389 + t777 * t416) * t498;
t809 = -t498 * t206 + t208 * t656 + t248 * t655 + (t498 * t778 + t639 * t795) * t735 + t253 * qJD(4);
t436 = -t495 * t779 + t696 * t778;
t730 = qJD(5) * t501;
t731 = qJD(5) * t497;
t749 = t338 * t497 - t386 * t731 - t436 * t730 + t501 * t745;
t314 = t386 * t501 - t436 * t497;
t748 = qJD(5) * t314 - t501 * t338 + t497 * t745;
t140 = t199 * t497 + t271 * t501;
t302 = -t389 * t777 + t779 * t416;
t385 = -t495 * t655 + t801;
t222 = t385 * pkin(4) - t386 * pkin(12) + t302;
t229 = -pkin(12) * t436 + t253;
t819 = t222 * t730 - t229 * t731 + t842 * t497 + t501 * t843;
t815 = -t338 * pkin(4) + t809;
t674 = t780 * qJD(1);
t629 = t674 + qJD(2);
t607 = t629 * t493;
t712 = t502 * t737;
t419 = t495 * t712 + t607;
t796 = t713 * (-t493 ^ 2 - t495 ^ 2);
t319 = -t776 * t419 + t778 * t796;
t320 = t778 * t419 + t776 * t796;
t242 = t319 * t691 + t320 * t795;
t836 = -t242 + t620;
t591 = t629 * t778;
t654 = qJD(1) * t694;
t460 = t499 * t654;
t623 = t693 * t737;
t592 = -t502 * t623 + t460;
t545 = -t493 * t591 + t592;
t462 = t493 * t712;
t590 = -t495 * t629 + t462;
t841 = t779 * t545 + t777 * t590;
t580 = t495 * t676 + t683;
t569 = t580 * t494;
t578 = t776 * t607;
t533 = qJD(1) * t569 + t578;
t243 = t498 * t533 + t795 * t841;
t523 = qJD(5) + t243;
t840 = -t744 * pkin(13) - t819;
t839 = t748 * pkin(5) - t749 * pkin(13) + t815;
t687 = t780 * t493;
t615 = t778 * t687;
t824 = t494 * (t495 * t682 - t677);
t547 = t824 + t615;
t598 = t493 * t756 - t495 * t780;
t803 = t547 * t779 - t598 * t777;
t489 = t494 ^ 2;
t835 = 0.2e1 * t489;
t619 = qJD(2) * t654;
t736 = qJD(2) * t499;
t573 = qJDD(1) * t824 - t502 * t619 - t623 * t736;
t670 = t780 * qJDD(1);
t622 = t670 + qJDD(2);
t589 = t778 * t622;
t538 = -t493 * t589 - t573;
t723 = qJDD(1) * t502;
t705 = t494 * t723;
t459 = t493 * t705;
t587 = t495 * t622 - t459;
t726 = qJD(1) * qJD(2);
t707 = t499 * t726;
t664 = t494 * t707;
t555 = t493 * t664 + t587;
t251 = -t538 * t777 - t555 * t779 - qJDD(4);
t633 = qJD(2) * t704;
t754 = t495 * t502;
t574 = qJD(3) * t754 + t499 * t633;
t645 = qJD(2) * t674;
t627 = pkin(1) * t645;
t663 = pkin(1) * t670;
t718 = pkin(10) * t705 + t499 * t663 + t502 * t627;
t258 = (qJ(3) * t622 + qJD(3) * t629) * t493 + (qJD(1) * t574 + t723 * t774) * t494 + t718;
t472 = t502 * t663;
t594 = -t499 * t627 + t472;
t734 = qJD(3) * t499;
t710 = t495 * t734;
t277 = t622 * pkin(2) + (qJDD(1) * t641 + (t502 * t633 - t710) * qJD(1)) * t494 + t594;
t560 = qJD(2) * t613 - t493 * t734;
t609 = -pkin(2) * t502 - t499 * t775 - pkin(1);
t301 = (qJD(1) * t560 + qJDD(1) * t609) * t494;
t132 = t778 * t258 + t277 * t692 + t301 * t695;
t828 = t538 * t779 - t555 * t777;
t109 = -pkin(11) * t828 + t132;
t131 = -t258 * t776 + t277 * t693 + t301 * t697;
t566 = t495 * t499 * t619 - t622 * t695;
t568 = t580 * qJDD(1);
t653 = qJD(2) * t682;
t823 = -(qJD(1) * t653 + t568) * t494 + t566;
t516 = pkin(11) * t823;
t110 = pkin(3) * t555 + t516 * t779 + t131;
t191 = -t493 * t277 + t495 * t301 + qJDD(3);
t128 = pkin(3) * t538 + t516 * t777 + t191;
t740 = pkin(10) * t756 + t482;
t429 = t740 * qJD(1);
t346 = qJ(3) * t419 + t429;
t576 = pkin(2) * t780 + t616;
t354 = qJD(2) * pkin(2) + qJD(1) * t576 + t473;
t407 = t609 * t737;
t227 = t778 * t346 + t354 * t692 + t407 * t695;
t162 = -pkin(11) * t841 + t227;
t226 = -t346 * t776 + t354 * t693 + t407 * t697;
t531 = pkin(11) * t533;
t166 = -pkin(3) * t590 - t531 * t779 + t226;
t283 = -t493 * t354 + t495 * t407 + qJD(3);
t189 = pkin(3) * t545 - t531 * t777 + t283;
t650 = qJD(4) * t690;
t652 = qJD(4) * t691;
t709 = qJD(4) * t795;
t30 = -t498 * t109 + t110 * t656 + t128 * t655 - t162 * t709 - t166 * t652 - t189 * t650;
t22 = pkin(4) * t251 - t30;
t122 = qJD(4) * t243 + t498 * t828 + t795 * t823;
t527 = t795 * t533;
t245 = -t498 * t841 + t527;
t556 = -t777 * t545 + t590 * t779 - qJD(4);
t553 = qJD(5) * t556;
t72 = t245 * t731 + t497 * t251 + (t122 + t553) * t501;
t673 = -t497 * t122 + t501 * t251;
t165 = t501 * t245 - t497 * t556;
t732 = qJD(5) * t165;
t73 = t673 + t732;
t10 = pkin(5) * t73 + pkin(13) * t72 + t22;
t512 = -t498 * t823 + t795 * t828;
t733 = qJD(4) * t498;
t123 = qJD(4) * t527 - t733 * t841 + t512;
t506 = qJDD(5) + t123;
t621 = qJD(4) * t656;
t583 = -t795 * t109 - t110 * t691 - t128 * t690 + t162 * t733 - t166 * t621 - t189 * t620;
t21 = -pkin(12) * t251 - t583;
t67 = -t110 * t777 + t779 * t128;
t37 = t123 * pkin(4) + t122 * pkin(12) + t67;
t78 = t795 * t162 + (t166 * t779 + t189 * t777) * t498;
t69 = -pkin(12) * t556 + t78;
t111 = -t166 * t777 + t779 * t189;
t75 = t243 * pkin(4) - t245 * pkin(12) + t111;
t612 = -t501 * t21 - t497 * t37 + t69 * t731 - t75 * t730;
t5 = pkin(13) * t506 - t612;
t41 = t497 * t75 + t501 * t69;
t34 = pkin(13) * t523 + t41;
t163 = t245 * t497 + t501 * t556;
t77 = -t498 * t162 + t166 * t656 + t189 * t655;
t68 = pkin(4) * t556 - t77;
t49 = t163 * pkin(5) - t165 * pkin(13) + t68;
t636 = t34 * t496 - t49 * t500;
t1 = -t636 * qJD(6) + t496 * t10 + t500 * t5;
t157 = qJD(6) + t163;
t834 = t636 * t157 + t1;
t808 = t497 * t222 + t501 * t229;
t818 = qJD(5) * t808 + t497 * t843 - t842 * t501;
t833 = qJD(6) * t523 - t72;
t832 = t315 * t498;
t441 = t497 * t690 - t501 * t779;
t701 = t319 * t777;
t743 = -qJD(5) * t441 + t497 * t701 + t501 * t836;
t657 = t780 * t794;
t445 = -t499 * t657 + t502 * t503;
t586 = t503 * t499 + t502 * t657;
t562 = t586 * t778;
t528 = t445 * t776 + t495 * t562 - t662 * t794;
t714 = t494 * t794;
t827 = t586 * t493 + t495 * t714;
t829 = t528 * t777 + t779 * t827;
t40 = -t497 * t69 + t501 * t75;
t826 = -t523 * t40 - t612;
t241 = -t319 * t656 + t320 * t498;
t797 = t650 - t241;
t353 = t812 * t777 + t626;
t150 = pkin(4) * t245 + pkin(12) * t243;
t62 = t497 * t150 + t501 * t77;
t825 = pkin(12) * t731 + pkin(13) * t245 + t62;
t822 = -pkin(12) * qJD(6) * t501 - t78 + t523 * (pkin(5) * t497 - pkin(13) * t501);
t16 = t34 * t500 + t49 * t496;
t2 = -qJD(6) * t16 + t500 * t10 - t496 * t5;
t821 = -t16 * t157 - t2;
t130 = pkin(13) * t385 + t808;
t228 = t436 * pkin(4) - t252;
t313 = t386 * t497 + t501 * t436;
t147 = t313 * pkin(5) - t314 * pkin(13) + t228;
t79 = -t130 * t496 + t147 * t500;
t820 = qJD(6) * t79 + t839 * t496 - t840 * t500;
t80 = t130 * t500 + t147 * t496;
t792 = -qJD(6) * t80 + t840 * t496 + t839 * t500;
t790 = -t744 * pkin(5) + t818;
t816 = t523 * t41;
t814 = t163 * t523;
t813 = t165 * t523;
t504 = qJD(1) ^ 2;
t758 = t489 * t504;
t811 = t501 * t506;
t689 = t501 * t777;
t442 = t497 * t779 + t498 * t689;
t742 = qJD(5) * t442 - t319 * t689 + t836 * t497;
t521 = qJD(5) * t523;
t807 = -t497 * t506 - t501 * t521;
t806 = qJD(4) * t245 + qJDD(5) + t512;
t805 = t528 * t779 - t777 * t827;
t804 = t803 * t795;
t411 = qJD(2) * t812;
t800 = qJD(2) * t625 - t779 * t411;
t798 = (qJDD(2) + 0.2e1 * t670) * t494;
t431 = t740 * qJD(2);
t339 = -qJD(2) * t626 - t411 * t777;
t474 = qJD(2) * t483;
t322 = qJD(3) * t687 + t494 * t574 + t474;
t349 = qJD(2) * t565 - t494 * t710;
t383 = t560 * t494;
t212 = t778 * t322 + t349 * t692 + t383 * t695;
t170 = pkin(11) * t800 + t212;
t211 = -t322 * t776 + t349 * t693 + t383 * t697;
t412 = qJD(2) * t421;
t711 = t494 * t736;
t666 = t493 * t711;
t174 = pkin(3) * t666 - t412 * t716 + t211;
t381 = (t494 * t754 + t687) * qJ(3) + t740;
t391 = t483 + t576;
t720 = t493 * t757;
t741 = pkin(2) * t756 + qJ(3) * t720;
t417 = -pkin(1) * t494 - t741;
t255 = t778 * t381 + t391 * t692 + t417 * t695;
t186 = pkin(11) * t803 + t255;
t254 = -t381 * t776 + t391 * t693 + t417 * t697;
t384 = t687 * t776 + t569;
t190 = -pkin(3) * t598 - t384 * t716 + t254;
t278 = -t493 * t349 + t495 * t383;
t215 = t411 * pkin(3) - t412 * t715 + t278;
t303 = -t493 * t391 + t495 * t417;
t221 = -pkin(3) * t547 - t384 * t715 + t303;
t50 = t795 * t170 + t174 * t691 - t186 * t733 + t190 * t621 + t215 * t690 + t221 * t620;
t45 = -pkin(12) * t339 + t50;
t113 = -t174 * t777 + t779 * t215;
t175 = -qJD(2) * t599 + t384 * t709 + t411 * t656 + t412 * t498 + t733 * t803;
t176 = qJD(4) * t804 - t384 * t733 + t412 * t795 + t498 * t800;
t66 = t175 * pkin(4) - t176 * pkin(12) + t113;
t310 = t547 * t777 + t598 * t779;
t92 = t795 * t186 + t190 * t691 + t221 * t690;
t85 = -pkin(12) * t310 + t92;
t127 = -t190 * t777 + t779 * t221;
t265 = t384 * t498 - t804;
t260 = t265 * pkin(4);
t266 = t384 * t795 + t498 * t803;
t88 = -t266 * pkin(12) + t127 + t260;
t791 = t497 * t88 + t501 * t85;
t14 = -qJD(5) * t791 - t45 * t497 + t501 * t66;
t729 = qJD(6) * t496;
t38 = t165 * t729 - t496 * t506 - t500 * t833;
t787 = t38 * t496;
t728 = qJD(6) * t500;
t39 = t165 * t728 + t496 * t833 - t500 * t506;
t786 = t39 * t500;
t71 = qJDD(6) + t73;
t785 = t496 * t71;
t784 = t500 * t71;
t783 = t68 * t243;
t466 = -pkin(5) * t501 - pkin(13) * t497 - pkin(4);
t782 = t466 * t728 + t496 * t822 - t500 * t825;
t781 = -t466 * t729 + t496 * t825 + t500 * t822;
t114 = t165 * t496 - t500 * t523;
t773 = t114 * t157;
t772 = t114 * t496;
t116 = t500 * t165 + t496 * t523;
t771 = t116 * t114;
t770 = t116 * t157;
t769 = t165 * t163;
t196 = t795 * t802 + t832;
t768 = t196 * t497;
t767 = t196 * t501;
t561 = t586 * t776;
t318 = t445 * t778 - t495 * t561 + t661 * t794;
t200 = t318 * t498 + t795 * t805;
t766 = t200 * t497;
t765 = t200 * t501;
t764 = t243 * t245;
t763 = t265 * t497;
t762 = t265 * t501;
t760 = t444 * t493;
t759 = t445 * t493;
t753 = t496 * t501;
t752 = t500 * t501;
t751 = t314 * t729 - t385 * t728 - t496 * t744 - t500 * t749;
t270 = t314 * t500 + t385 * t496;
t750 = qJD(6) * t270 + t496 * t749 - t500 * t744;
t401 = -t496 * t442 - t500 * t655;
t747 = qJD(6) * t401 + t496 * t797 + t500 * t743;
t585 = -t500 * t442 + t496 * t655;
t746 = qJD(6) * t585 - t496 * t743 + t500 * t797;
t739 = t503 * pkin(1) + pkin(10) * t714;
t491 = t499 ^ 2;
t492 = t502 ^ 2;
t738 = t491 - t492;
t724 = qJDD(1) * t499;
t721 = t502 * t758;
t708 = pkin(1) * t835;
t706 = t502 * t726;
t703 = -t497 * t21 + t501 * t37;
t675 = qJD(3) * t776;
t672 = t523 * t497;
t671 = t157 * t500;
t669 = t499 * t721;
t665 = t499 * t706;
t658 = -pkin(1) * t794 + pkin(10) * t755;
t651 = t494 * t504 * t780;
t201 = t318 * t795 - t498 * t805;
t142 = t201 * t497 - t501 * t829;
t648 = g(1) * t140 + g(2) * t142;
t198 = -t316 * t656 - t398 * t655 - t832;
t647 = g(1) * t198 + g(2) * t200;
t646 = -g(1) * t445 - g(2) * t444;
t144 = -t243 * t753 - t500 * t245;
t644 = -t496 * t730 + t144;
t145 = -t243 * t752 + t245 * t496;
t643 = t500 * t730 - t145;
t637 = -t16 * t496 + t500 * t636;
t43 = pkin(13) * t265 + t791;
t187 = t266 * t497 + t310 * t501;
t188 = t266 * t501 - t310 * t497;
t91 = -t498 * t186 + t190 * t656 + t221 * t655;
t84 = t310 * pkin(4) - t91;
t58 = t187 * pkin(5) - t188 * pkin(13) + t84;
t18 = t43 * t500 + t496 * t58;
t17 = -t43 * t496 + t500 * t58;
t47 = -t497 * t85 + t501 * t88;
t61 = t150 * t501 - t497 * t77;
t134 = t188 * t500 + t265 * t496;
t133 = t188 * t496 - t265 * t500;
t136 = t222 * t501 - t229 * t497;
t628 = 0.2e1 * t674 + qJD(2);
t617 = g(1) * t503 + g(2) * t794;
t13 = t501 * t45 + t497 * t66 + t88 * t730 - t731 * t85;
t33 = -pkin(5) * t523 - t40;
t611 = -pkin(13) * t71 + t157 * t33;
t608 = -t706 - t724;
t606 = g(1) * t142 - g(2) * t140 + g(3) * t187;
t143 = t201 * t501 + t497 * t829;
t605 = -g(1) * t143 + g(2) * t141 - g(3) * t188;
t341 = t443 * t776 - t444 * t693;
t342 = -t443 * t778 - t444 * t692;
t237 = t342 * t795 + (t341 * t779 + t444 * t696) * t498;
t298 = -t341 * t777 + t444 * t698;
t158 = t237 * t497 - t298 * t501;
t343 = -t445 * t693 + t561;
t344 = -t445 * t692 - t562;
t239 = t344 * t795 + (t343 * t779 + t445 * t696) * t498;
t299 = -t343 * t777 + t445 * t698;
t160 = t239 * t497 - t299 * t501;
t295 = t421 * t795 + (-t779 * t812 + t625) * t498;
t246 = t295 * t497 - t353 * t501;
t604 = -g(1) * t160 - g(2) * t158 - g(3) * t246;
t603 = g(1) * t200 + g(2) * t196 + g(3) * t265;
t602 = g(1) * t201 - g(2) * t199 + g(3) * t266;
t236 = -t341 * t656 + t342 * t498 - t444 * t624;
t238 = -t343 * t656 + t344 * t498 - t445 * t624;
t294 = t421 * t498 + t656 * t812 - t599;
t601 = -g(1) * t238 - g(2) * t236 - g(3) * t294;
t572 = t69 * t730 + t731 * t75 - t703;
t6 = -pkin(5) * t506 + t572;
t597 = -t6 + t606;
t596 = t608 * pkin(10);
t584 = t421 * pkin(3) + pkin(11) * t353 + t741;
t577 = -t444 * pkin(2) + qJ(3) * t398 + t658;
t567 = pkin(13) * qJD(6) * t157 - t597;
t559 = t295 * pkin(4) + t294 * pkin(12) + t584;
t552 = t445 * pkin(2) + qJ(3) * t827 + t739;
t51 = -t498 * t170 + t174 * t656 - t186 * t709 - t190 * t652 + t215 * t655 - t221 * t650;
t432 = t443 * pkin(2);
t551 = t342 * pkin(3) + pkin(11) * t298 + qJ(3) * t760 - t432;
t434 = t586 * pkin(2);
t550 = t344 * pkin(3) + pkin(11) * t299 + qJ(3) * t759 - t434;
t549 = -t315 * pkin(3) - pkin(11) * t271 + t577;
t46 = t339 * pkin(4) - t51;
t542 = t493 * t545;
t541 = t237 * pkin(4) + t236 * pkin(12) + t551;
t540 = t239 * pkin(4) + t238 * pkin(12) + t550;
t539 = t199 * pkin(4) + t198 * pkin(12) + t549;
t535 = t493 * t538;
t534 = t580 * t737 + t578;
t526 = t533 * t720;
t522 = t243 * t523;
t515 = t493 * t823;
t513 = -g(1) * t827 + g(2) * t398 + g(3) * t598;
t510 = t318 * pkin(3) + pkin(11) * t829 + t552;
t508 = t201 * pkin(4) + t200 * pkin(12) + t510;
t446 = -pkin(10) * t757 + t483;
t439 = -t493 * t678 + t477;
t430 = -pkin(10) * t711 + t474;
t428 = -pkin(10) * t713 + t473;
t425 = pkin(12) * t752 + t466 * t496;
t424 = -pkin(12) * t753 + t466 * t500;
t358 = t494 * t596 + t594;
t357 = -pkin(10) * t664 + t718;
t269 = t314 * t496 - t500 * t385;
t247 = t295 * t501 + t353 * t497;
t194 = t200 * pkin(4);
t192 = t196 * pkin(4);
t161 = t239 * t501 + t299 * t497;
t159 = t237 * t501 + t298 * t497;
t129 = -pkin(5) * t385 - t136;
t104 = -qJD(5) * t187 + t176 * t501 - t339 * t497;
t103 = qJD(5) * t188 + t176 * t497 + t339 * t501;
t102 = pkin(5) * t165 + pkin(13) * t163;
t99 = t143 * t500 + t200 * t496;
t98 = -t143 * t496 + t200 * t500;
t60 = -qJD(6) * t133 + t104 * t500 + t175 * t496;
t59 = qJD(6) * t134 + t104 * t496 - t175 * t500;
t54 = -pkin(5) * t245 - t61;
t42 = -pkin(5) * t265 - t47;
t28 = t102 * t496 + t40 * t500;
t27 = t102 * t500 - t40 * t496;
t19 = t103 * pkin(5) - t104 * pkin(13) + t46;
t12 = -pkin(5) * t175 - t14;
t11 = pkin(13) * t175 + t13;
t8 = -qJD(5) * t41 + t703;
t4 = -qJD(6) * t18 - t11 * t496 + t19 * t500;
t3 = qJD(6) * t17 + t11 * t500 + t19 * t496;
t7 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t794 - g(2) * t503, t617, 0, 0 (qJDD(1) * t491 + 0.2e1 * t665) * t489 (t499 * t723 - t726 * t738) * t835, qJD(2) * t628 * t756 + t499 * t798 (qJDD(1) * t492 - 0.2e1 * t665) * t489, t502 * t798 - t628 * t711, t622 * t780, -t431 * t629 + t446 * t622 + t358 * t780 + g(1) * t444 - g(2) * t445 + (-t707 + t723) * t708, -g(1) * t443 + g(2) * t586 - t357 * t780 - t430 * t629 + t608 * t708 - t622 * t740 ((-t428 * qJD(2) + qJDD(1) * t740 + t357 + (-qJD(2) * t446 + t430) * qJD(1)) * t502 + (-t429 * qJD(2) - qJDD(1) * t446 - t358) * t499 - t617) * t494, t489 * qJDD(1) * pkin(1) ^ 2 - g(1) * t658 - g(2) * t739 + t357 * t740 + t358 * t446 - t428 * t431 + t429 * t430, -t566 * t384 + t412 * t578 + (t384 * t568 + (t384 * t653 + t412 * t580) * qJD(1)) * t494, -t384 * t538 - t411 * t533 - t412 * t545 - t547 * t823, qJD(2) * t526 + t384 * t555 - t412 * t590 + t598 * t823, t411 * t545 - t538 * t547, t411 * t590 + t538 * t598 - t542 * t711 + t547 * t555, -t555 * t598 - t590 * t666, g(1) * t315 - g(2) * t318 - t131 * t598 - t191 * t547 - t211 * t590 + t226 * t666 + t254 * t555 + t278 * t545 + t283 * t411 + t303 * t538, -g(1) * t316 + g(2) * t528 + t132 * t598 + t191 * t384 + t212 * t590 - t227 * t666 - t255 * t555 + t278 * t533 + t283 * t412 - t303 * t823, -g(1) * t398 - g(2) * t827 - t131 * t384 + t132 * t547 - t211 * t534 - t212 * t545 - t226 * t412 - t227 * t411 + t254 * t823 - t255 * t538, -g(1) * t577 - g(2) * t552 + t131 * t254 + t132 * t255 + t191 * t303 + t226 * t211 + t227 * t212 + t283 * t278, -t122 * t266 + t176 * t245, t122 * t265 - t123 * t266 - t175 * t245 - t176 * t243, t310 * t122 - t176 * t556 - t339 * t245 - t251 * t266, t123 * t265 + t175 * t243, t310 * t123 + t175 * t556 + t339 * t243 + t251 * t265, t251 * t310 + t339 * t556, -g(1) * t199 - g(2) * t201 + t111 * t175 + t113 * t243 + t127 * t123 - t91 * t251 + t67 * t265 - t30 * t310 - t77 * t339 - t51 * t556, t111 * t176 + t113 * t245 - t127 * t122 + t92 * t251 + t67 * t266 - t310 * t583 + t78 * t339 + t50 * t556 + t647, g(1) * t271 - g(2) * t829 + t91 * t122 - t92 * t123 - t78 * t175 - t77 * t176 - t50 * t243 - t51 * t245 + t265 * t583 - t30 * t266, -g(1) * t549 - g(2) * t510 + t111 * t113 + t67 * t127 + t30 * t91 + t78 * t50 + t77 * t51 - t583 * t92, t104 * t165 - t188 * t72, -t103 * t165 - t104 * t163 + t187 * t72 - t188 * t73, t104 * t523 + t165 * t175 + t188 * t506 - t72 * t265, t103 * t163 + t187 * t73, -t103 * t523 - t163 * t175 - t187 * t506 - t73 * t265, t175 * t523 + t265 * t806, -g(1) * t141 - g(2) * t143 + t68 * t103 + t14 * t523 + t46 * t163 + t40 * t175 + t22 * t187 + t8 * t265 + t47 * t506 + t84 * t73, t68 * t104 - t13 * t523 + t46 * t165 - t41 * t175 + t22 * t188 + t265 * t612 - t506 * t791 - t84 * t72 + t648, -t103 * t41 - t104 * t40 - t13 * t163 - t14 * t165 + t187 * t612 - t188 * t8 + t47 * t72 - t73 * t791 - t647, -g(1) * t539 - g(2) * t508 + t41 * t13 + t40 * t14 + t22 * t84 + t68 * t46 + t8 * t47 - t612 * t791, t116 * t60 - t134 * t38, -t114 * t60 - t116 * t59 + t133 * t38 - t134 * t39, t103 * t116 + t134 * t71 + t157 * t60 - t187 * t38, t114 * t59 + t133 * t39, -t103 * t114 - t133 * t71 - t157 * t59 - t187 * t39, t103 * t157 + t187 * t71, t4 * t157 + t17 * t71 + t2 * t187 - t636 * t103 + t12 * t114 + t42 * t39 + t6 * t133 + t33 * t59 - g(1) * (t198 * t496 + t844) - g(2) * t99, -t3 * t157 - t18 * t71 - t1 * t187 - t16 * t103 + t12 * t116 - t42 * t38 + t6 * t134 + t33 * t60 - g(1) * (t198 * t500 - t845) - g(2) * t98, -t1 * t133 - t114 * t3 - t116 * t4 - t134 * t2 - t16 * t59 + t17 * t38 - t18 * t39 + t60 * t636 - t648, t1 * t18 + t16 * t3 + t2 * t17 - t636 * t4 + t6 * t42 + t33 * t12 - g(1) * (t141 * pkin(5) + t140 * pkin(13) + t539) - g(2) * (t143 * pkin(5) + t142 * pkin(13) + t508); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t669, t738 * t758, t494 * t724 - t502 * t651, t669, t499 * t651 + t705, t622, t472 + t429 * t629 + g(1) * t586 + g(2) * t443 + (-t645 + t758) * t499 * pkin(1) + (-g(3) * t502 + t596) * t494, t428 * t629 + pkin(1) * t721 + (pkin(10) * t726 + g(3)) * t757 - t646 - t718, 0, 0, -t413 * t533 - t515 * t776, t410 * t533 + t413 * t545 - t515 * t778 - t535 * t776, -qJD(1) * t526 + t413 * t590 - t495 * t823 + t555 * t695, -t410 * t545 - t535 * t778, -t459 * t697 - t410 * t462 - (qJD(1) * t615 - t460) * t667 + (t410 * t629 + (-t499 * t682 * t758 + 0.2e1 * t589) * t493 + t573) * t495, t587 * t495 - (t495 * t674 - t462) * t667, t439 * t587 + t131 * t495 - t309 * t592 - t283 * t410 + t261 * t590 - g(1) * t344 - g(2) * t342 - g(3) * t421 + (-t778 * t191 + t590 * t675 + t309 * t591 + (qJD(2) * t439 - t226) * t713 - t538 * pkin(2)) * t493, pkin(2) * t515 - g(1) * t343 - g(2) * t341 + g(3) * t812 - t132 * t495 + t191 * t695 + t227 * t667 - t283 * t413 - t309 * t533 - t440 * t555 + (qJD(3) * t697 - t262) * t590, -qJD(3) * t542 * t778 - g(1) * t759 - g(2) * t760 - g(3) * t720 - t131 * t695 + t132 * t697 + t226 * t413 + t227 * t410 + t262 * t545 + t439 * t823 - t440 * t538 + (t493 * t675 + t261) * t534, t132 * t440 + t131 * t439 - t227 * t262 - t226 * t261 - t283 * t309 + g(1) * t434 + g(2) * t432 - g(3) * t741 + (-t191 * pkin(2) + (-t226 * t776 + t227 * t778) * qJD(3) + t646 * qJ(3)) * t493, -t122 * t386 + t245 * t745, t122 * t385 - t123 * t386 - t243 * t745 - t245 * t744, t436 * t122 + t338 * t245 - t251 * t386 - t556 * t745, t123 * t385 + t243 * t744, t436 * t123 - t338 * t243 + t251 * t385 + t556 * t744, t251 * t436 - t338 * t556, -g(1) * t239 - g(2) * t237 - g(3) * t295 + t111 * t744 + t302 * t123 + t243 * t588 - t252 * t251 - t30 * t436 + t77 * t338 + t67 * t385 + t556 * t809, t111 * t745 - t302 * t122 + t245 * t588 + t253 * t251 - t78 * t338 + t67 * t386 - t436 * t583 + t556 * t810 - t601, -g(1) * t299 - g(2) * t298 - g(3) * t353 + t122 * t252 - t123 * t253 - t243 * t810 + t245 * t809 - t30 * t386 + t385 * t583 - t744 * t78 - t745 * t77, -g(1) * t550 - g(2) * t551 - g(3) * t584 + t111 * t588 + t30 * t252 - t253 * t583 + t67 * t302 - t77 * t809 + t78 * t810, t165 * t749 - t314 * t72, -t163 * t749 - t165 * t748 + t313 * t72 - t314 * t73, t165 * t744 + t314 * t506 - t72 * t385 + t523 * t749, t163 * t748 + t313 * t73, -t163 * t744 - t313 * t506 - t73 * t385 - t523 * t748, t385 * t806 + t523 * t744, -g(1) * t161 - g(2) * t159 - g(3) * t247 + t136 * t506 + t163 * t815 + t22 * t313 + t228 * t73 + t8 * t385 + t40 * t744 - t523 * t818 + t68 * t748, t165 * t815 + t22 * t314 - t228 * t72 + t385 * t612 - t41 * t744 - t506 * t808 - t523 * t819 + t68 * t749 - t604, t136 * t72 - t163 * t819 + t165 * t818 + t313 * t612 - t314 * t8 - t749 * t40 - t748 * t41 - t73 * t808 + t601, -g(1) * t540 - g(2) * t541 - g(3) * t559 + t8 * t136 + t22 * t228 - t40 * t818 + t41 * t819 - t612 * t808 + t68 * t815, -t116 * t751 - t270 * t38, t114 * t751 - t116 * t750 + t269 * t38 - t270 * t39, t116 * t748 - t157 * t751 + t270 * t71 - t313 * t38, t114 * t750 + t269 * t39, -t114 * t748 - t157 * t750 - t269 * t71 - t313 * t39, t157 * t748 + t313 * t71, t79 * t71 + t2 * t313 + t129 * t39 + t6 * t269 - g(1) * (t161 * t500 + t238 * t496) - g(2) * (t159 * t500 + t236 * t496) - g(3) * (t247 * t500 + t294 * t496) + t750 * t33 + t792 * t157 - t748 * t636 + t790 * t114, -t80 * t71 - t1 * t313 - t129 * t38 + t6 * t270 - g(1) * (-t161 * t496 + t238 * t500) - g(2) * (-t159 * t496 + t236 * t500) - g(3) * (-t247 * t496 + t294 * t500) - t751 * t33 - t748 * t16 - t820 * t157 + t790 * t116, -t1 * t269 - t114 * t820 - t116 * t792 - t16 * t750 - t2 * t270 + t38 * t79 - t39 * t80 - t636 * t751 + t604, t1 * t80 + t2 * t79 + t6 * t129 - g(1) * (t161 * pkin(5) + t160 * pkin(13) + t540) - g(2) * (t159 * pkin(5) + t158 * pkin(13) + t541) - g(3) * (t247 * pkin(5) + t246 * pkin(13) + t559) + t790 * t33 + t820 * t16 - t792 * t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t319 * t590 + t538, -t320 * t590 - t823, t319 * t534 + t320 * t545, -t226 * t319 - t227 * t320 + t191 + t513, 0, 0, 0, 0, 0, 0, t123 * t779 + t243 * t701 - t251 * t655 + t556 * t797, -t122 * t779 + t245 * t701 + t251 * t690 + t836 * t556, t122 * t655 - t123 * t690 - t241 * t245 + t242 * t243 + (-t243 * t655 + t245 * t690) * qJD(4), -t583 * t690 + t30 * t655 + t67 * t779 - t78 * t242 + t77 * t241 + t111 * t701 + (t655 * t78 - t690 * t77) * qJD(4) + t513, 0, 0, 0, 0, 0, 0, t163 * t797 - t441 * t506 - t523 * t742 - t655 * t73, t165 * t797 - t442 * t506 - t523 * t743 + t655 * t72, -t163 * t743 + t165 * t742 - t441 * t72 - t442 * t73, -t22 * t655 - t742 * t40 + t743 * t41 - t8 * t441 - t442 * t612 + t68 * t797 + t513, 0, 0, 0, 0, 0, 0, t114 * t742 + t157 * t746 + t39 * t441 + t401 * t71, t116 * t742 - t157 * t747 - t38 * t441 + t585 * t71, -t114 * t747 - t116 * t746 + t38 * t401 + t39 * t585, -t1 * t585 + t16 * t747 + t2 * t401 + t33 * t742 + t6 * t441 - t636 * t746 + t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t764, -t243 ^ 2 + t245 ^ 2, -t243 * t556 - t122, -t764, -t245 * t556 - t123, -t251, -t111 * t245 - t556 * t78 + t30 + t603, t111 * t243 - t556 * t77 + t583 + t602, 0, 0, -t497 * t72 + t501 * t813 (-t72 - t814) * t501 + (-t73 - t813) * t497, -t165 * t245 + t501 * t522 - t807, t163 * t672 - t501 * t73, t163 * t245 + t811 + (-t521 - t522) * t497, -t523 * t245, -pkin(4) * t73 + pkin(12) * t807 + g(1) * t765 + g(2) * t767 + g(3) * t762 - t78 * t163 - t22 * t501 - t40 * t245 + t497 * t783 - t523 * t61 + t68 * t731, pkin(4) * t72 - g(1) * t766 - g(2) * t768 - g(3) * t763 - t78 * t165 + t22 * t497 + t41 * t245 + t501 * t783 + t523 * t62 + t68 * t730 + (t523 * t731 - t811) * pkin(12), t163 * t62 + t165 * t61 + ((-t73 + t732) * pkin(12) + t826) * t501 + (-t8 - t816 + (qJD(5) * t163 - t72) * pkin(12)) * t497 - t602, -t22 * pkin(4) + g(1) * t194 + g(2) * t192 + g(3) * t260 - t40 * t61 - t41 * t62 - t68 * t78 + (-t8 * t497 - t612 * t501 + (-t40 * t501 - t41 * t497) * qJD(5) - t602) * pkin(12), -t38 * t497 * t500 + (-t497 * t729 + t643) * t116, t114 * t145 + t116 * t144 + (-t114 * t500 - t116 * t496) * t730 + (t787 - t786 + (-t116 * t500 + t772) * qJD(6)) * t497, t38 * t501 + t643 * t157 + (t116 * t523 - t157 * t729 + t784) * t497, t39 * t496 * t497 + (t497 * t728 - t644) * t114, t39 * t501 + t644 * t157 + (-t114 * t523 - t157 * t728 - t785) * t497, t157 * t672 - t501 * t71, -t54 * t114 - t33 * t144 + t424 * t71 + t781 * t157 - t602 * t496 + (-t2 + (pkin(12) * t114 + t33 * t496) * qJD(5) + t603 * t500) * t501 + (pkin(12) * t39 + t33 * t728 + t6 * t496 - t523 * t636) * t497, -t54 * t116 - t33 * t145 - t425 * t71 - t782 * t157 - t602 * t500 + (t1 + (pkin(12) * t116 + t33 * t500) * qJD(5) - t603 * t496) * t501 + (-pkin(12) * t38 - t16 * t523 - t33 * t729 + t6 * t500) * t497, t144 * t16 - t145 * t636 + t38 * t424 - t39 * t425 - t781 * t116 - t782 * t114 + t637 * t730 + (-t1 * t496 - t2 * t500 + (-t16 * t500 - t496 * t636) * qJD(6) + t603) * t497, t1 * t425 + t2 * t424 - t33 * t54 - g(1) * (-pkin(5) * t765 - pkin(13) * t766 - t194) - g(2) * (-pkin(5) * t767 - pkin(13) * t768 - t192) - g(3) * (-pkin(5) * t762 - pkin(13) * t763 - t260) + t782 * t16 - t781 * t636 + (t33 * t730 + t497 * t6 - t602) * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t769, -t163 ^ 2 + t165 ^ 2, -t72 + t814, -t769, -t245 * t730 + t497 * t553 - t673 + t813, t506, -t68 * t165 - t572 + t606 + t816, t68 * t163 - t605 - t826, 0, 0, t116 * t671 - t787 (-t38 - t773) * t500 + (-t39 - t770) * t496, -t116 * t165 + t157 * t671 + t785, t157 * t772 - t786, -t157 ^ 2 * t496 + t114 * t165 + t784, -t157 * t165, -pkin(5) * t39 - t114 * t41 - t157 * t27 + t165 * t636 + t496 * t611 - t500 * t567, pkin(5) * t38 - t116 * t41 + t157 * t28 + t16 * t165 + t496 * t567 + t500 * t611, t114 * t28 + t116 * t27 + ((qJD(6) * t116 - t39) * pkin(13) + t834) * t500 + ((qJD(6) * t114 - t38) * pkin(13) + t821) * t496 + t605, t636 * t27 - t16 * t28 - t33 * t41 + t597 * pkin(5) + (qJD(6) * t637 + t1 * t500 - t2 * t496 + t605) * pkin(13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t771, -t114 ^ 2 + t116 ^ 2, -t38 + t773, -t771, -t39 + t770, t71, -t33 * t116 - g(1) * t98 - g(2) * (t196 * t500 + t845) + g(3) * t133 - t821, t33 * t114 + g(1) * t99 - g(2) * (-t196 * t496 + t844) + g(3) * t134 - t834, 0, 0;];
tau_reg  = t7;
