% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:26
% EndTime: 2019-03-10 05:41:20
% DurationCPUTime: 65.90s
% Computational Cost: add. (85802->1249), mult. (241707->1700), div. (0->0), fcn. (204516->18), ass. (0->496)
t529 = cos(qJ(2));
t767 = cos(pkin(6));
t701 = pkin(1) * t767;
t508 = t529 * t701;
t493 = qJD(1) * t508;
t526 = sin(qJ(2));
t521 = sin(pkin(6));
t766 = cos(pkin(7));
t608 = t521 * (-pkin(10) * t766 - pkin(9));
t587 = t526 * t608;
t374 = qJD(1) * t587 + t493;
t507 = t526 * t701;
t552 = t529 * t608 - t507;
t375 = t552 * qJD(1);
t765 = sin(pkin(7));
t674 = t529 * t765;
t585 = pkin(2) * t526 - pkin(10) * t674;
t722 = qJD(1) * t521;
t416 = t585 * t722;
t794 = cos(qJ(3));
t648 = t766 * t794;
t525 = sin(qJ(3));
t679 = t525 * t765;
t443 = pkin(2) * t648 - pkin(10) * t679;
t680 = t525 * t766;
t730 = t443 * qJD(3) - t794 * t374 - t375 * t680 - t416 * t679;
t572 = -t526 * t680 + t529 * t794;
t412 = t572 * t521;
t396 = qJD(1) * t412;
t647 = t765 * t794;
t612 = qJD(3) * t647;
t582 = t612 - t396;
t289 = -t375 * t765 + t766 * t416;
t570 = t525 * t529 + t526 * t648;
t411 = t570 * t521;
t395 = qJD(1) * t411;
t852 = -t395 * pkin(3) + t396 * pkin(11) - t289 + (pkin(3) * t679 - pkin(11) * t647) * qJD(3);
t682 = t521 * t765;
t646 = qJD(1) * t682;
t617 = t526 * t646;
t851 = pkin(11) * t617 - t730;
t524 = sin(qJ(4));
t528 = cos(qJ(4));
t437 = t524 * t679 - t528 * t766;
t732 = qJD(4) * t437 + t524 * t617 - t582 * t528;
t676 = t528 * t765;
t438 = t524 * t766 + t525 * t676;
t731 = qJD(4) * t438 + t582 * t524 + t528 * t617;
t445 = pkin(2) * t680 + pkin(10) * t647;
t728 = t445 * qJD(3) - t525 * t374 + t375 * t648 + t416 * t647;
t643 = qJD(3) * t679;
t605 = t643 - t395;
t422 = pkin(11) * t766 + t445;
t699 = t765 * pkin(2);
t423 = -pkin(3) * t647 - pkin(11) * t679 - t699;
t717 = qJD(4) * t528;
t718 = qJD(4) * t524;
t738 = -t422 * t718 + t423 * t717 + t852 * t524 - t851 * t528;
t729 = pkin(3) * t617 + t728;
t698 = t526 * t722;
t660 = t525 * t698;
t673 = t767 * qJD(1);
t628 = t673 + qJD(2);
t575 = t765 * t628;
t621 = t529 * t648;
t599 = t521 * t621;
t821 = -qJD(1) * t599 - t794 * t575;
t338 = t660 + t821;
t334 = qJD(4) + t338;
t850 = -pkin(12) * t605 - t738;
t849 = t731 * pkin(4) + t732 * pkin(12) + t729;
t795 = cos(qJ(1));
t650 = t767 * t795;
t792 = sin(qJ(1));
t439 = t526 * t650 + t529 * t792;
t583 = t792 * t526 - t529 * t650;
t838 = t583 * t766 + t795 * t682;
t309 = -t439 * t794 + t525 * t838;
t683 = t521 * t766;
t819 = t583 * t765 - t795 * t683;
t236 = t309 * t528 - t524 * t819;
t306 = t439 * t525 + t794 * t838;
t520 = qJ(5) + qJ(6);
t514 = sin(t520);
t515 = cos(t520);
t848 = t236 * t514 + t306 * t515;
t847 = t236 * t515 - t306 * t514;
t523 = sin(qJ(5));
t527 = cos(qJ(5));
t846 = t236 * t523 + t306 * t527;
t845 = t236 * t527 - t306 * t523;
t750 = t521 * t529;
t725 = pkin(9) * t750 + t507;
t426 = t725 * qJD(1);
t675 = t529 * t766;
t652 = t521 * t675;
t331 = t426 + (qJD(1) * t652 + t575) * pkin(10);
t555 = pkin(2) * t767 + t587;
t337 = qJD(2) * pkin(2) + qJD(1) * t555 + t493;
t678 = t526 * t765;
t580 = -pkin(2) * t529 - pkin(10) * t678 - pkin(1);
t389 = t580 * t722;
t197 = t794 * t331 + (t337 * t766 + t389 * t765) * t525;
t844 = -t197 + t334 * (pkin(4) * t524 - pkin(12) * t528);
t571 = t525 * t675 + t526 * t794;
t561 = t571 * t521;
t566 = t525 * t575;
t340 = qJD(1) * t561 + t566;
t747 = t523 * t528;
t240 = -t338 * t747 - t527 * t340;
t822 = -t523 * t717 + t240;
t381 = t523 * t438 + t527 * t647;
t736 = qJD(5) * t381 - t523 * t605 + t527 * t732;
t622 = t523 * t647;
t715 = qJD(5) * t527;
t735 = -qJD(5) * t622 + t438 * t715 - t523 * t732 - t527 * t605;
t421 = -pkin(3) * t766 - t443;
t302 = t437 * pkin(4) - t438 * pkin(12) + t421;
t316 = t528 * t422 + t524 * t423;
t305 = -pkin(12) * t647 + t316;
t716 = qJD(5) * t523;
t781 = t302 * t715 - t305 * t716 + t849 * t523 - t527 * t850;
t204 = t523 * t302 + t527 * t305;
t776 = -qJD(5) * t204 + t523 * t850 + t849 * t527;
t839 = -t524 * t715 + t822;
t837 = t309 * t524 + t528 * t819;
t836 = pkin(5) * t731 + pkin(13) * t736 + t776;
t835 = pkin(13) * t735 - t781;
t196 = -t525 * t331 + t337 * t648 + t389 * t647;
t246 = pkin(3) * t340 + pkin(11) * t338;
t146 = t528 * t196 + t524 * t246;
t119 = pkin(12) * t340 + t146;
t803 = -pkin(4) * t528 - pkin(12) * t524;
t481 = -pkin(3) + t803;
t770 = -t527 * t119 + t481 * t715 + (-t527 * t718 - t528 * t716) * pkin(11) + t844 * t523;
t788 = pkin(11) * t523;
t834 = t119 * t523 + t527 * t844 + t718 * t788;
t737 = -t422 * t717 - t423 * t718 + t851 * t524 + t852 * t528;
t517 = t521 ^ 2;
t831 = 0.2e1 * t517;
t744 = t527 * t528;
t241 = -t338 * t744 + t340 * t523;
t509 = pkin(11) * t744;
t790 = pkin(5) * t524;
t830 = pkin(13) * t241 + t338 * t790 + (-pkin(13) * t744 + t790) * qJD(4) + (-t509 + (pkin(13) * t524 - t481) * t523) * qJD(5) + t834;
t829 = pkin(13) * t839 + t770;
t530 = -pkin(13) - pkin(12);
t700 = qJD(5) * t530;
t712 = t529 * t646 - qJD(3);
t558 = -t628 * t766 + t712;
t270 = t528 * t340 - t524 * t558;
t670 = t524 * t340 + t528 * t558;
t183 = pkin(4) * t270 + pkin(12) * t670;
t262 = -t337 * t765 + t766 * t389;
t172 = t338 * pkin(3) - t340 * pkin(11) + t262;
t179 = -pkin(11) * t558 + t197;
t97 = t172 * t528 - t524 * t179;
t73 = t523 * t183 + t527 * t97;
t757 = t670 * t523;
t828 = -pkin(13) * t757 + t523 * t700 - t73;
t72 = t527 * t183 - t523 * t97;
t827 = -pkin(5) * t270 - t72 + (-pkin(13) * t670 + t700) * t527;
t205 = t523 * t270 - t334 * t527;
t207 = t270 * t527 + t334 * t523;
t522 = sin(qJ(6));
t793 = cos(qJ(6));
t121 = t793 * t205 + t207 * t522;
t604 = -t522 * t205 + t207 * t793;
t764 = t121 * t604;
t801 = qJD(5) + t670;
t826 = t205 * t801;
t825 = t207 * t801;
t692 = t793 * qJD(6);
t824 = t527 * (t793 * qJD(5) + t692);
t739 = -pkin(4) * t605 - t737;
t649 = t767 * t792;
t565 = t526 * t795 + t529 * t649;
t820 = t565 * t765 + t792 * t683;
t635 = t527 * t717 - t241;
t818 = -t524 * t716 + t635;
t664 = t767 * qJDD(1);
t615 = t664 + qJDD(2);
t654 = t521 * t674;
t554 = -qJDD(1) * t654 + t615 * t766 + qJDD(3);
t655 = t521 * t678;
t616 = qJD(2) * t655;
t540 = qJD(1) * t616 + t554;
t662 = pkin(9) * t698;
t637 = qJD(2) * t673;
t624 = pkin(1) * t637;
t656 = pkin(1) * t664;
t709 = qJDD(1) * t529;
t688 = t521 * t709;
t703 = pkin(9) * t688 + t526 * t656 + t529 * t624;
t335 = -qJD(2) * t662 + t703;
t574 = t615 * t765;
t669 = qJDD(1) * t766;
t638 = t529 * t669;
t672 = t766 * qJD(2);
t639 = qJD(1) * t672;
t256 = (t574 + (-t526 * t639 + t638) * t521) * pkin(10) + t335;
t492 = t529 * t656;
t577 = -t526 * t624 + t492;
t711 = qJD(1) * qJD(2);
t690 = t529 * t711;
t710 = qJDD(1) * t526;
t595 = -t690 - t710;
t578 = t595 * pkin(9);
t264 = t615 * pkin(2) + ((-t526 * t669 - t529 * t639) * pkin(10) + t578) * t521 + t577;
t573 = qJD(2) * t585;
t312 = (qJD(1) * t573 + qJDD(1) * t580) * t521;
t613 = qJD(3) * t648;
t720 = qJD(3) * t525;
t569 = -t794 * t256 - t264 * t680 - t312 * t679 + t331 * t720 - t337 * t613 - t389 * t612;
t83 = pkin(11) * t540 - t569;
t186 = -t264 * t765 + t766 * t312;
t689 = t521 * t710;
t216 = (t672 + qJD(3)) * t660 - t525 * t574 - t794 * t689 + (-t525 * t638 - t690 * t794) * t521 + t821 * qJD(3);
t533 = qJD(2) * t570 + qJD(3) * t571;
t217 = qJD(3) * t566 - qJDD(1) * t599 + t521 * (qJD(1) * t533 + t525 * t710) - t794 * t574;
t89 = t217 * pkin(3) + t216 * pkin(11) + t186;
t663 = t172 * t718 + t179 * t717 + t524 * t83 - t528 * t89;
t98 = t172 * t524 + t179 * t528;
t817 = -t334 * t98 + t663;
t816 = -t121 ^ 2 + t604 ^ 2;
t178 = pkin(3) * t558 - t196;
t108 = pkin(4) * t670 - t270 * pkin(12) + t178;
t88 = pkin(12) * t334 + t98;
t53 = t108 * t523 + t527 * t88;
t213 = qJDD(4) + t217;
t34 = t172 * t717 - t179 * t718 + t524 * t89 + t528 * t83;
t24 = pkin(12) * t213 + t34;
t802 = qJD(4) * t670;
t137 = t528 * t216 - t524 * t540 + t802;
t719 = qJD(4) * t270;
t138 = -t524 * t216 - t528 * t540 + t719;
t645 = qJD(3) * t680;
t694 = qJD(3) * t794;
t607 = t525 * t256 - t264 * t648 - t312 * t647 + t331 * t694 + t337 * t645 + t389 * t643;
t84 = -pkin(3) * t540 + t607;
t40 = t138 * pkin(4) + t137 * pkin(12) + t84;
t9 = -t53 * qJD(5) - t523 * t24 + t527 * t40;
t815 = -t801 * t53 - t9;
t705 = -t523 * t137 + t270 * t715 + t334 * t716;
t614 = t527 * t213 - t705;
t714 = qJD(6) * t522;
t74 = t527 * t137 - t523 * t213 + t270 * t716 - t334 * t715;
t21 = t205 * t692 + t207 * t714 - t522 * t614 + t793 * t74;
t265 = qJD(6) + t801;
t814 = t121 * t265 - t21;
t52 = t527 * t108 - t523 * t88;
t44 = -pkin(13) * t207 + t52;
t37 = pkin(5) * t801 + t44;
t45 = -pkin(13) * t205 + t53;
t136 = qJDD(5) + t138;
t6 = pkin(5) * t136 + pkin(13) * t74 + t9;
t602 = -t108 * t715 - t527 * t24 - t523 * t40 + t716 * t88;
t7 = pkin(13) * t614 - t602;
t1 = t37 * t692 - t45 * t714 + t522 * t6 + t7 * t793;
t440 = -t526 * t649 + t529 * t795;
t806 = t565 * t766 - t792 * t682;
t311 = t440 * t794 - t525 * t806;
t238 = t311 * t528 + t524 * t820;
t310 = t440 * t525 + t794 * t806;
t168 = t238 * t515 + t310 * t514;
t633 = t765 * t767;
t609 = t525 * t633;
t367 = t609 + t561;
t634 = t767 * t766;
t560 = t654 - t634;
t301 = t367 * t528 - t524 * t560;
t586 = t794 * t633;
t745 = t525 * t526;
t366 = t521 * t745 - t586 - t599;
t87 = -pkin(4) * t334 - t97;
t66 = pkin(5) * t205 + t87;
t813 = t66 * t121 + g(1) * t168 - g(2) * t847 - g(3) * (-t301 * t515 - t366 * t514) - t1;
t812 = -t52 * t801 - t602;
t203 = t527 * t302 - t305 * t523;
t382 = t527 * t438 - t622;
t162 = pkin(5) * t437 - pkin(13) * t382 + t203;
t180 = -pkin(13) * t381 + t204;
t94 = t162 * t793 - t522 * t180;
t810 = qJD(6) * t94 + t522 * t836 - t835 * t793;
t95 = t522 * t162 + t180 * t793;
t782 = -qJD(6) * t95 + t835 * t522 + t793 * t836;
t773 = pkin(5) * t735 + t739;
t809 = -t334 * t97 + t34;
t808 = t524 * t334;
t373 = t508 + t555;
t726 = pkin(2) * t750 + pkin(10) * t655;
t405 = -pkin(1) * t521 - t726;
t283 = -t373 * t765 + t766 * t405;
t363 = t366 * pkin(3);
t192 = -t367 * pkin(11) + t283 + t363;
t355 = (t652 + t633) * pkin(10) + t725;
t220 = t794 * t355 + t373 * t680 + t405 * t679;
t200 = -pkin(11) * t560 + t220;
t114 = t524 * t192 + t528 * t200;
t111 = pkin(12) * t366 + t114;
t219 = -t525 * t355 + t373 * t648 + t405 * t647;
t199 = pkin(3) * t560 - t219;
t300 = t367 * t524 + t528 * t560;
t128 = t300 * pkin(4) - t301 * pkin(12) + t199;
t62 = t527 * t111 + t523 * t128;
t419 = t523 * t481 + t509;
t805 = (qJDD(2) + 0.2e1 * t664) * t521;
t430 = t725 * qJD(2);
t800 = qJD(5) + qJD(6);
t174 = -t238 * t523 + t310 * t527;
t229 = t301 * t523 - t366 * t527;
t799 = -g(1) * t174 - g(2) * t846 + g(3) * t229;
t167 = -t238 * t514 + t310 * t515;
t707 = t793 * t45;
t16 = t522 * t37 + t707;
t2 = -qJD(6) * t16 - t522 * t7 + t793 * t6;
t798 = -t66 * t604 - g(1) * t167 - g(2) * t848 - g(3) * (-t301 * t514 + t366 * t515) + t2;
t22 = qJD(6) * t604 - t522 * t74 - t793 * t614;
t797 = t265 * t604 - t22;
t531 = qJD(1) ^ 2;
t789 = pkin(9) * t521;
t784 = t213 * pkin(4);
t778 = t522 * t45;
t777 = t74 * t523;
t459 = t527 * t481;
t746 = t524 * t527;
t365 = -pkin(13) * t746 + t459 + (-pkin(5) - t788) * t528;
t748 = t523 * t524;
t384 = -pkin(13) * t748 + t419;
t275 = t365 * t793 - t522 * t384;
t775 = qJD(6) * t275 + t522 * t830 + t793 * t829;
t276 = t522 * t365 + t384 * t793;
t774 = -qJD(6) * t276 - t522 * t829 + t793 * t830;
t487 = t530 * t523;
t488 = t530 * t527;
t391 = t487 * t793 + t522 * t488;
t772 = qJD(6) * t391 + t522 * t827 + t793 * t828;
t392 = t522 * t487 - t488 * t793;
t771 = -qJD(6) * t392 - t522 * t828 + t793 * t827;
t769 = -qJD(5) * t419 + t834;
t145 = -t524 * t196 + t246 * t528;
t118 = -pkin(4) * t340 - t145;
t768 = -pkin(5) * t839 + pkin(11) * t717 - t118;
t763 = t136 * t523;
t762 = t136 * t527;
t761 = t205 * t523;
t760 = t207 * t205;
t759 = t670 * t270;
t758 = t670 * t334;
t756 = t270 * t334;
t755 = t340 * t338;
t754 = t514 * t528;
t753 = t515 * t528;
t752 = t517 * t531;
t751 = t521 * t526;
t749 = t522 * t523;
t743 = t381 * t692 + t382 * t714 + t522 * t735 + t736 * t793;
t278 = -t522 * t381 + t382 * t793;
t742 = qJD(6) * t278 - t522 * t736 + t735 * t793;
t461 = t522 * t527 + t523 * t793;
t372 = t800 * t461;
t741 = -t461 * t670 - t372;
t603 = t527 * t793 - t749;
t740 = -t603 * t670 + t749 * t800 - t824;
t659 = t793 * t717;
t734 = -t241 * t793 - t372 * t524 + t522 * t822 + t527 * t659;
t733 = -t240 * t793 + t522 * t818 + t523 * t659 + t524 * t824 - t714 * t748;
t724 = t795 * pkin(1) + t792 * t789;
t518 = t526 ^ 2;
t519 = t529 ^ 2;
t723 = t518 - t519;
t721 = qJD(2) * t521;
t706 = t529 * t752;
t704 = t412 * pkin(3) + t726;
t702 = pkin(5) * t523 + pkin(11);
t697 = t526 * t721;
t691 = pkin(1) * t831;
t685 = t439 * t765;
t684 = t440 * t765;
t681 = t524 * t765;
t61 = -t111 * t523 + t527 * t128;
t113 = t192 * t528 - t524 * t200;
t315 = -t524 * t422 + t423 * t528;
t667 = t527 * t801;
t665 = t334 * t528;
t661 = t526 * t706;
t658 = t526 * t690;
t657 = -t98 + (t716 + t757) * pkin(5);
t651 = -pkin(1) * t792 + t795 * t789;
t644 = t521 * t531 * t767;
t237 = t311 * t524 - t528 * t820;
t641 = -g(1) * t837 - g(2) * t237;
t640 = -g(1) * t306 + g(2) * t310;
t304 = pkin(4) * t647 - t315;
t632 = -t52 * t527 - t523 * t53;
t377 = t552 * qJD(2);
t417 = t521 * t573;
t290 = -t377 * t765 + t766 * t417;
t230 = t301 * t527 + t366 * t523;
t513 = pkin(5) * t527 + pkin(4);
t631 = -t513 * t528 + t524 * t530;
t627 = 0.2e1 * t673 + qJD(2);
t625 = pkin(11) * t411 + t704;
t494 = qJD(2) * t508;
t376 = qJD(2) * t587 + t494;
t150 = -t355 * t720 + t373 * t613 + t794 * t376 + t377 * t680 + t405 * t612 + t417 * t679;
t141 = pkin(11) * t616 + t150;
t287 = qJD(3) * t609 + t521 * t533;
t288 = qJD(3) * t586 + ((t621 - t745) * qJD(3) + t572 * qJD(2)) * t521;
t158 = t287 * pkin(3) - t288 * pkin(11) + t290;
t51 = -t524 * t141 + t158 * t528 - t192 * t718 - t200 * t717;
t611 = -t583 * pkin(2) + pkin(10) * t685;
t610 = -t565 * pkin(2) + pkin(10) * t684;
t110 = -pkin(4) * t366 - t113;
t49 = pkin(5) * t300 - pkin(13) * t230 + t61;
t54 = -pkin(13) * t229 + t62;
t19 = t49 * t793 - t522 * t54;
t20 = t522 * t49 + t54 * t793;
t25 = t663 - t784;
t156 = -t522 * t229 + t230 * t793;
t600 = -pkin(12) * t136 + t801 * t87;
t50 = t528 * t141 + t524 * t158 + t192 * t717 - t200 * t718;
t47 = pkin(12) * t287 + t50;
t151 = -t355 * t694 - t373 * t645 - t525 * t376 + t377 * t648 - t405 * t643 + t417 * t647;
t142 = -pkin(3) * t616 - t151;
t184 = qJD(4) * t301 + t288 * t524 - t528 * t616;
t185 = -qJD(4) * t300 + t288 * t528 + t524 * t616;
t65 = t184 * pkin(4) - t185 * pkin(12) + t142;
t12 = -t111 * t716 + t128 * t715 + t527 * t47 + t523 * t65;
t598 = -pkin(11) * t213 + t178 * t334;
t328 = -t439 * t680 - t583 * t794;
t597 = t328 * pkin(3) + t611;
t330 = -t440 * t680 - t565 * t794;
t596 = t330 * pkin(3) + t610;
t594 = g(1) * t795 + g(2) * t792;
t593 = g(1) * t237 - g(2) * t837 + g(3) * t300;
t592 = -g(1) * t238 + g(2) * t236 - g(3) * t301;
t271 = t328 * t524 - t439 * t676;
t273 = t330 * t524 - t440 * t676;
t341 = t412 * t524 - t528 * t655;
t591 = -g(1) * t273 - g(2) * t271 - g(3) * t341;
t590 = g(1) * t310 + g(2) * t306 + g(3) * t366;
t589 = g(1) * t311 - g(2) * t309 + g(3) * t367;
t327 = t439 * t648 - t525 * t583;
t329 = t440 * t648 - t525 * t565;
t588 = -g(1) * t329 - g(2) * t327 - g(3) * t411;
t584 = t614 * t527;
t579 = -t25 + t593;
t48 = -pkin(4) * t287 - t51;
t568 = t327 * pkin(11) + t597;
t567 = t329 * pkin(11) + t596;
t13 = -qJD(5) * t62 - t47 * t523 + t527 * t65;
t557 = pkin(11) * qJD(4) * t334 - t590 + t84;
t556 = pkin(12) * qJD(5) * t801 - t579;
t548 = t558 * t765;
t547 = qJD(3) * t548;
t544 = -t439 * pkin(2) - pkin(10) * t819 + t651;
t543 = t440 * pkin(2) + pkin(10) * t820 + t724;
t539 = t309 * pkin(3) + t544;
t538 = t311 * pkin(3) + t543;
t536 = -pkin(11) * t306 + t539;
t535 = t310 * pkin(11) + t538;
t534 = t540 * t765;
t469 = t702 * t524;
t444 = -pkin(9) * t751 + t508;
t432 = t603 * t524;
t431 = t461 * t524;
t429 = -pkin(9) * t697 + t494;
t424 = t493 - t662;
t418 = -pkin(11) * t747 + t459;
t342 = t412 * t528 + t524 * t655;
t336 = t521 * t578 + t577;
t298 = t310 * pkin(3);
t296 = t306 * pkin(3);
t277 = t381 * t793 + t382 * t522;
t274 = t330 * t528 + t440 * t681;
t272 = t328 * t528 + t439 * t681;
t242 = pkin(5) * t381 + t304;
t175 = t238 * t527 + t310 * t523;
t155 = t229 * t793 + t230 * t522;
t133 = qJDD(6) + t136;
t105 = -qJD(5) * t229 + t185 * t527 + t287 * t523;
t104 = qJD(5) * t230 + t185 * t523 - t287 * t527;
t82 = pkin(5) * t229 + t110;
t42 = qJD(6) * t156 + t104 * t793 + t522 * t105;
t41 = t522 * t104 - t105 * t793 + t229 * t692 + t230 * t714;
t36 = pkin(5) * t104 + t48;
t18 = t44 * t793 - t778;
t17 = -t522 * t44 - t707;
t15 = t37 * t793 - t778;
t14 = -pkin(5) * t614 + t25;
t11 = -pkin(13) * t104 + t12;
t10 = pkin(5) * t184 - pkin(13) * t105 + t13;
t4 = -qJD(6) * t20 + t10 * t793 - t522 * t11;
t3 = qJD(6) * t19 + t522 * t10 + t11 * t793;
t5 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t792 - g(2) * t795, t594, 0, 0 (qJDD(1) * t518 + 0.2e1 * t658) * t517 (t526 * t709 - t711 * t723) * t831, t529 * t627 * t721 + t526 * t805 (qJDD(1) * t519 - 0.2e1 * t658) * t517, t529 * t805 - t627 * t697, t615 * t767, -t430 * t628 + t444 * t615 + t336 * t767 + g(1) * t439 - g(2) * t440 + (-t526 * t711 + t709) * t691, -g(1) * t583 + g(2) * t565 - t335 * t767 - t429 * t628 + t595 * t691 - t615 * t725 ((-t424 * qJD(2) + qJDD(1) * t725 + t335 + (-qJD(2) * t444 + t429) * qJD(1)) * t529 + (-t426 * qJD(2) - qJDD(1) * t444 - t336) * t526 - t594) * t521, t517 * qJDD(1) * pkin(1) ^ 2 - g(1) * t651 - g(2) * t724 + t335 * t725 + t336 * t444 - t424 * t430 + t426 * t429, -t216 * t367 + t288 * t340, t216 * t366 - t217 * t367 - t287 * t340 - t288 * t338, t216 * t560 - t288 * t558 + t340 * t616 + t367 * t540, t217 * t366 + t287 * t338, t217 * t560 + t287 * t558 - t338 * t616 - t366 * t540, -t540 * t560 - t548 * t697, -g(1) * t309 - g(2) * t311 - t151 * t558 + t186 * t366 + t196 * t616 + t283 * t217 + t219 * t540 + t262 * t287 + t290 * t338 + t560 * t607, t150 * t558 + t186 * t367 - t197 * t616 - t283 * t216 - t220 * t540 + t262 * t288 + t290 * t340 - t560 * t569 + t640, g(1) * t819 - g(2) * t820 - t150 * t338 - t151 * t340 - t196 * t288 - t197 * t287 + t219 * t216 - t220 * t217 + t366 * t569 + t367 * t607, -g(1) * t544 - g(2) * t543 + t197 * t150 + t196 * t151 + t186 * t283 - t219 * t607 - t220 * t569 + t262 * t290, -t137 * t301 + t185 * t270, t137 * t300 - t138 * t301 - t184 * t270 - t185 * t670, -t137 * t366 + t185 * t334 + t213 * t301 + t270 * t287, t138 * t300 + t184 * t670, -t138 * t366 - t184 * t334 - t213 * t300 - t287 * t670, t213 * t366 + t287 * t334, -g(1) * t236 - g(2) * t238 + t113 * t213 + t138 * t199 + t142 * t670 + t178 * t184 + t287 * t97 + t300 * t84 + t334 * t51 - t366 * t663, -t114 * t213 - t137 * t199 + t142 * t270 + t178 * t185 - t287 * t98 + t301 * t84 - t334 * t50 - t34 * t366 - t641, t113 * t137 - t114 * t138 - t184 * t98 - t185 * t97 - t270 * t51 - t300 * t34 + t301 * t663 - t50 * t670 - t640, -g(1) * t536 - g(2) * t535 - t113 * t663 + t34 * t114 + t178 * t142 + t84 * t199 + t98 * t50 + t97 * t51, t105 * t207 - t230 * t74, -t207 * t104 - t105 * t205 + t74 * t229 + t230 * t614, t105 * t801 + t136 * t230 + t184 * t207 - t300 * t74, t205 * t104 - t229 * t614, -t104 * t801 - t136 * t229 - t184 * t205 + t300 * t614, t136 * t300 + t184 * t801, -g(1) * t845 - g(2) * t175 + t87 * t104 - t110 * t614 + t13 * t801 + t61 * t136 + t52 * t184 + t48 * t205 + t25 * t229 + t9 * t300, g(1) * t846 - g(2) * t174 + t87 * t105 - t110 * t74 - t12 * t801 - t62 * t136 - t53 * t184 + t48 * t207 + t25 * t230 + t602 * t300, -t53 * t104 - t52 * t105 - t12 * t205 - t13 * t207 + t229 * t602 - t9 * t230 + t61 * t74 + t614 * t62 + t641, -t602 * t62 + t53 * t12 + t9 * t61 + t52 * t13 + t25 * t110 + t87 * t48 - g(1) * (t236 * pkin(4) + pkin(12) * t837 + t536) - g(2) * (t238 * pkin(4) + t237 * pkin(12) + t535) -t156 * t21 - t41 * t604, t121 * t41 + t155 * t21 - t156 * t22 - t42 * t604, t133 * t156 + t184 * t604 - t21 * t300 - t265 * t41, t121 * t42 + t155 * t22, -t121 * t184 - t133 * t155 - t22 * t300 - t265 * t42, t133 * t300 + t184 * t265, -g(1) * t847 - g(2) * t168 + t36 * t121 + t19 * t133 + t14 * t155 + t15 * t184 + t2 * t300 + t82 * t22 + t4 * t265 + t66 * t42, g(1) * t848 - g(2) * t167 - t1 * t300 - t20 * t133 + t14 * t156 - t16 * t184 - t82 * t21 - t3 * t265 + t36 * t604 - t66 * t41, -t1 * t155 - t121 * t3 + t15 * t41 - t156 * t2 - t16 * t42 + t19 * t21 - t20 * t22 - t4 * t604 + t641, t1 * t20 + t16 * t3 + t2 * t19 + t15 * t4 + t14 * t82 + t66 * t36 - g(1) * (t236 * t513 - t306 * t702 - t530 * t837 + t539) - g(2) * (-t237 * t530 + t238 * t513 + t310 * t702 + t538); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t661, t723 * t752, -t529 * t644 + t689, t661, t526 * t644 + t688, t615, t492 + t426 * t628 + g(1) * t565 + g(2) * t583 + (-t637 + t752) * t526 * pkin(1) + (-g(3) * t529 + t578) * t521, t424 * t628 + pkin(1) * t706 + g(1) * t440 + g(2) * t439 + (pkin(9) * t711 + g(3)) * t751 - t703, 0, 0, -t216 * t679 + t340 * t582, -t216 * t647 - t217 * t679 + t396 * t338 + t340 * t395 + (-t338 * t647 - t340 * t679) * qJD(3), -t216 * t766 - t340 * t617 + t396 * t558 + t525 * t534 - t547 * t794, -t217 * t647 + t338 * t605, -t217 * t766 + t338 * t617 - t395 * t558 + t525 * t547 + t534 * t794, t554 * t766 - (qJD(1) * t634 - t712) * t617, -g(1) * t330 - g(2) * t328 - g(3) * t412 - t186 * t647 - t196 * t617 - t217 * t699 + t262 * t605 - t289 * t338 + t443 * t540 + t558 * t728 - t607 * t766, t186 * t679 + t197 * t617 + t216 * t699 + t262 * t582 - t289 * t340 - t445 * t540 + t558 * t730 + t569 * t766 - t588, -g(3) * t655 - t569 * t647 - g(1) * t684 - g(2) * t685 + t607 * t679 + t196 * t396 + t197 * t395 + t443 * t216 - t445 * t217 + t728 * t340 - t730 * t338 + (-t196 * t647 - t197 * t679) * qJD(3), -g(1) * t610 - g(2) * t611 - g(3) * t726 - t186 * t699 - t196 * t728 + t197 * t730 - t262 * t289 - t443 * t607 - t445 * t569, -t137 * t438 - t270 * t732, t137 * t437 - t138 * t438 - t270 * t731 + t670 * t732, t137 * t647 + t213 * t438 + t270 * t605 - t334 * t732, t138 * t437 + t670 * t731, t138 * t647 - t437 * t213 - t334 * t731 - t605 * t670, -t213 * t647 + t334 * t605, -g(1) * t274 - g(2) * t272 - g(3) * t342 + t421 * t138 + t178 * t731 + t315 * t213 + t334 * t737 + t84 * t437 + t605 * t97 + t647 * t663 + t670 * t729, -t421 * t137 - t178 * t732 - t316 * t213 + t270 * t729 - t334 * t738 + t34 * t647 + t84 * t438 - t605 * t98 - t591, t137 * t315 - t138 * t316 - t270 * t737 - t34 * t437 + t438 * t663 - t670 * t738 - t731 * t98 + t732 * t97 + t588, -g(1) * t567 - g(2) * t568 - g(3) * t625 + t178 * t729 - t315 * t663 + t34 * t316 + t84 * t421 + t737 * t97 + t738 * t98, -t207 * t736 - t382 * t74, t205 * t736 - t207 * t735 + t74 * t381 + t382 * t614, t136 * t382 + t207 * t731 - t437 * t74 - t736 * t801, t205 * t735 - t381 * t614, -t136 * t381 - t205 * t731 + t437 * t614 - t735 * t801, t136 * t437 + t731 * t801, t203 * t136 + t9 * t437 - t304 * t614 + t25 * t381 - g(1) * (t274 * t527 + t329 * t523) - g(2) * (t272 * t527 + t327 * t523) - g(3) * (t342 * t527 + t411 * t523) + t735 * t87 + t731 * t52 + t776 * t801 + t739 * t205, -t204 * t136 + t602 * t437 - t304 * t74 + t25 * t382 - g(1) * (-t274 * t523 + t329 * t527) - g(2) * (-t272 * t523 + t327 * t527) - g(3) * (-t342 * t523 + t411 * t527) - t736 * t87 - t731 * t53 - t781 * t801 + t739 * t207, t203 * t74 + t204 * t614 - t205 * t781 - t207 * t776 + t381 * t602 - t9 * t382 + t52 * t736 - t53 * t735 + t591, -t602 * t204 + t9 * t203 + t25 * t304 - g(1) * (t274 * pkin(4) + t273 * pkin(12) + t567) - g(2) * (t272 * pkin(4) + t271 * pkin(12) + t568) - g(3) * (pkin(4) * t342 + pkin(12) * t341 + t625) + t739 * t87 + t781 * t53 + t776 * t52, -t21 * t278 - t604 * t743, t121 * t743 + t21 * t277 - t22 * t278 - t604 * t742, t133 * t278 - t21 * t437 - t265 * t743 + t604 * t731, t121 * t742 + t22 * t277, -t121 * t731 - t133 * t277 - t22 * t437 - t265 * t742, t133 * t437 + t265 * t731, t94 * t133 + t2 * t437 + t242 * t22 + t14 * t277 - g(1) * (t274 * t515 + t329 * t514) - g(2) * (t272 * t515 + t327 * t514) - g(3) * (t342 * t515 + t411 * t514) + t742 * t66 + t782 * t265 + t731 * t15 + t773 * t121, -t95 * t133 - t1 * t437 - t242 * t21 + t14 * t278 - g(1) * (-t274 * t514 + t329 * t515) - g(2) * (-t272 * t514 + t327 * t515) - g(3) * (-t342 * t514 + t411 * t515) - t743 * t66 - t810 * t265 - t731 * t16 + t773 * t604, -t1 * t277 - t121 * t810 + t15 * t743 - t16 * t742 - t2 * t278 + t21 * t94 - t22 * t95 - t604 * t782 + t591, t1 * t95 + t2 * t94 + t14 * t242 - g(1) * (-t273 * t530 + t274 * t513 + t329 * t702 + t596) - g(2) * (-t271 * t530 + t272 * t513 + t327 * t702 + t597) - g(3) * (-t341 * t530 + t342 * t513 + t411 * t702 + t704) + t773 * t66 + t810 * t16 + t782 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t755, -t338 ^ 2 + t340 ^ 2, -t338 * t558 - t216, -t755, -t340 * t558 - t217, t540, -t197 * t558 - t262 * t340 + t590 - t607, -t196 * t558 + t262 * t338 + t569 + t589, 0, 0, -t137 * t524 + t270 * t665 (-t137 - t758) * t528 + (-t138 - t756) * t524, t213 * t524 - t270 * t340 + t334 * t665, -t138 * t528 + t670 * t808, t213 * t528 - t334 * t808 + t340 * t670, -t334 * t340, -pkin(3) * t138 - t145 * t334 - t197 * t670 - t340 * t97 + t524 * t598 - t528 * t557, pkin(3) * t137 + t146 * t334 - t197 * t270 + t340 * t98 + t524 * t557 + t528 * t598, t145 * t270 + t146 * t670 + ((-t138 + t719) * pkin(11) + t809) * t528 + ((-t137 + t802) * pkin(11) + t817) * t524 - t589, -t84 * pkin(3) + g(1) * t298 + g(2) * t296 + g(3) * t363 - t97 * t145 - t98 * t146 - t178 * t197 + (t34 * t528 + t663 * t524 + (-t524 * t98 - t528 * t97) * qJD(4) - t589) * pkin(11), t207 * t818 - t74 * t746, t241 * t205 + t207 * t240 + (-t205 * t527 - t207 * t523) * t717 + (t584 + t777 + (-t207 * t527 + t761) * qJD(5)) * t524, t528 * t74 + t635 * t801 + (t207 * t334 - t716 * t801 + t762) * t524, -t205 * t839 - t614 * t748, -t528 * t614 + t822 * t801 + (-t205 * t334 - t715 * t801 - t763) * t524, -t136 * t528 + t801 * t808, -t118 * t205 + t418 * t136 - t87 * t240 + t769 * t801 - t589 * t523 + (-t9 + (pkin(11) * t205 + t523 * t87) * qJD(4) + t590 * t527) * t528 + (-pkin(11) * t614 + t25 * t523 + t334 * t52 + t715 * t87) * t524, -t118 * t207 - t419 * t136 - t87 * t241 - t770 * t801 - t589 * t527 + (-t602 + (pkin(11) * t207 + t527 * t87) * qJD(4) - t590 * t523) * t528 + (-pkin(11) * t74 + t25 * t527 - t334 * t53 - t716 * t87) * t524, t419 * t614 + t418 * t74 + t53 * t240 + t52 * t241 - t769 * t207 - t770 * t205 + t632 * t717 + (t602 * t523 - t9 * t527 + (t52 * t523 - t527 * t53) * qJD(5) + t590) * t524, -t602 * t419 + t9 * t418 - t87 * t118 - g(1) * (t310 * t803 - t298) - g(2) * (t306 * t803 - t296) - g(3) * (t366 * t803 - t363) + t770 * t53 + t769 * t52 + (t25 * t524 + t717 * t87 - t589) * pkin(11), -t21 * t432 + t604 * t734, -t121 * t734 + t21 * t431 - t22 * t432 - t604 * t733, t133 * t432 + t21 * t528 + t265 * t734 + t604 * t808, t121 * t733 + t22 * t431, -t121 * t808 - t133 * t431 + t22 * t528 - t265 * t733, -t133 * t528 + t265 * t808, t275 * t133 - t2 * t528 + t469 * t22 + t14 * t431 - g(1) * (-t310 * t753 + t311 * t514) - g(2) * (-t306 * t753 - t309 * t514) - g(3) * (-t366 * t753 + t367 * t514) + t733 * t66 + t774 * t265 + t15 * t808 + t768 * t121, -t276 * t133 + t1 * t528 - t469 * t21 + t14 * t432 - g(1) * (t310 * t754 + t311 * t515) - g(2) * (t306 * t754 - t309 * t515) - g(3) * (t366 * t754 + t367 * t515) + t734 * t66 - t775 * t265 - t16 * t808 + t768 * t604, -t1 * t431 - t121 * t775 - t15 * t734 - t16 * t733 - t2 * t432 + t21 * t275 - t22 * t276 + t524 * t590 - t604 * t774, t1 * t276 + t2 * t275 + t14 * t469 - g(1) * (t310 * t631 + t311 * t702 - t298) - g(2) * (t306 * t631 - t309 * t702 - t296) - g(3) * (t366 * t631 + t367 * t702 - t363) + t768 * t66 + t775 * t16 + t774 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t759, t270 ^ 2 - t670 ^ 2, -t137 + t758, -t759, t756 - t138, t213, -t178 * t270 + t593 - t817, t178 * t670 - t592 - t809, 0, 0, t207 * t667 - t777 (-t74 - t826) * t527 + (t614 - t825) * t523, -t207 * t270 + t667 * t801 + t763, t761 * t801 + t584, -t523 * t801 ^ 2 + t205 * t270 + t762, -t801 * t270, -pkin(4) * t705 - t72 * t801 - t52 * t270 - t98 * t205 + t600 * t523 + (-t556 + t784) * t527, pkin(4) * t74 - t207 * t98 + t270 * t53 + t523 * t556 + t527 * t600 + t73 * t801, t73 * t205 + t72 * t207 + ((qJD(5) * t207 + t614) * pkin(12) + t812) * t527 + ((qJD(5) * t205 - t74) * pkin(12) + t815) * t523 + t592, -t52 * t72 - t53 * t73 - t87 * t98 + t579 * pkin(4) + (qJD(5) * t632 - t9 * t523 - t527 * t602 + t592) * pkin(12), -t21 * t461 - t604 * t740, t121 * t740 - t21 * t603 - t22 * t461 + t604 * t741, t133 * t461 - t265 * t740 - t270 * t604, -t121 * t741 - t22 * t603, t121 * t270 + t133 * t603 + t265 * t741, -t265 * t270, t121 * t657 + t133 * t391 - t14 * t603 - t15 * t270 - t22 * t513 + t265 * t771 + t515 * t593 - t66 * t741, -t133 * t392 + t14 * t461 + t16 * t270 + t21 * t513 - t265 * t772 - t514 * t593 + t604 * t657 - t66 * t740, t1 * t603 - t121 * t772 + t15 * t740 + t16 * t741 - t2 * t461 + t21 * t391 - t22 * t392 - t604 * t771 + t592, t1 * t392 + t2 * t391 - t14 * t513 - g(1) * (-t237 * t513 - t238 * t530) - g(2) * (t236 * t530 + t513 * t837) - g(3) * (-t300 * t513 - t301 * t530) + t657 * t66 + t772 * t16 + t771 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t760, -t205 ^ 2 + t207 ^ 2, -t74 + t826, -t760, t614 + t825, t136, -t87 * t207 + t799 - t815, g(1) * t175 - g(2) * t845 + g(3) * t230 + t87 * t205 - t812, 0, 0, t764, t816, t814, -t764, t797, t133, -t17 * t265 + (-t121 * t207 + t133 * t793 - t265 * t714) * pkin(5) + t798, t18 * t265 + (-t133 * t522 - t207 * t604 - t265 * t692) * pkin(5) + t813, t16 * t604 + t18 * t121 - t15 * t121 + t17 * t604 + (t793 * t21 - t22 * t522 + (-t121 * t793 + t522 * t604) * qJD(6)) * pkin(5), -t15 * t17 - t16 * t18 + (t1 * t522 + t2 * t793 - t66 * t207 + (-t15 * t522 + t16 * t793) * qJD(6) + t799) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t764, t816, t814, -t764, t797, t133, t16 * t265 + t798, t15 * t265 + t813, 0, 0;];
tau_reg  = t5;
