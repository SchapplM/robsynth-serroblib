% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:14
% EndTime: 2019-12-31 22:02:40
% DurationCPUTime: 17.43s
% Computational Cost: add. (12658->689), mult. (27192->891), div. (0->0), fcn. (27264->6), ass. (0->497)
t496 = cos(qJ(3));
t770 = -pkin(8) - pkin(7);
t452 = t770 * t496;
t743 = cos(qJ(4));
t444 = t743 * t452;
t493 = sin(qJ(4));
t494 = sin(qJ(3));
t672 = t493 * t494;
t569 = t770 * t672;
t362 = -t444 + t569;
t598 = t743 * t496;
t435 = -t598 + t672;
t683 = t435 * qJ(5);
t264 = t362 - t683;
t768 = -t264 / 0.2e1;
t497 = cos(qJ(2));
t745 = -t497 / 0.2e1;
t495 = sin(qJ(2));
t669 = t494 * t495;
t446 = pkin(3) * t669 + pkin(6) * t495;
t599 = t743 * t494;
t671 = t493 * t496;
t439 = t599 + t671;
t405 = t439 * t495;
t739 = t405 * pkin(4);
t321 = t446 + t739;
t483 = -pkin(3) * t496 - pkin(2);
t735 = t435 * pkin(4);
t393 = t483 + t735;
t753 = t439 / 0.2e1;
t408 = -t493 * t669 + t495 * t598;
t757 = t408 / 0.2e1;
t657 = t321 * t753 + t393 * t757;
t801 = t264 * t745 - t657;
t779 = -t452 * t493 - t599 * t770;
t789 = -t779 / 0.2e1;
t807 = t408 * t768;
t728 = t497 * pkin(2);
t559 = -pkin(7) * t495 - t728;
t542 = -pkin(1) + t559;
t665 = t496 * t497;
t613 = pkin(6) * t665;
t384 = t494 * t542 + t613;
t342 = -pkin(8) * t669 + t384;
t326 = t743 * t342;
t514 = (t495 * t770 - pkin(1) - t728) * t496;
t742 = pkin(6) * t494;
t504 = (-pkin(3) - t742) * t497 + t514;
t501 = t493 * t504;
t184 = t326 + t501;
t689 = t405 * qJ(5);
t158 = t184 - t689;
t785 = qJ(5) * t439 + t779;
t767 = t785 / 0.2e1;
t799 = -t785 / 0.2e1;
t802 = t799 + t767;
t806 = t802 * t158;
t668 = t494 * t497;
t616 = pkin(6) * t668;
t341 = t514 - t616;
t676 = t493 * t341;
t198 = -t326 - t676;
t163 = t198 + t689;
t805 = t158 + t163;
t790 = -t483 / 0.2e1;
t791 = -t446 / 0.2e1;
t804 = t362 * t745 + t408 * t790 + t439 * t791;
t295 = t743 * t504;
t675 = t493 * t342;
t183 = -t295 + t675;
t600 = t743 * t341;
t199 = t600 - t675;
t796 = t183 + t199;
t617 = qJD(3) + qJD(4);
t800 = t405 * t790 + t435 * t791 + t497 * t789;
t727 = t497 * pkin(7);
t731 = t495 * pkin(2);
t453 = -t727 + t731;
t445 = t496 * t453;
t474 = pkin(6) * t669;
t394 = t474 + t445;
t730 = t495 * pkin(3);
t325 = -pkin(8) * t665 + t394 + t730;
t297 = t743 * t325;
t442 = t494 * t453;
t666 = t496 * t495;
t615 = pkin(6) * t666;
t395 = t442 - t615;
t347 = -pkin(8) * t668 + t395;
t674 = t493 * t347;
t655 = t297 / 0.2e1 - t674 / 0.2e1;
t77 = t655 + t804;
t80 = t655 - t804;
t762 = -t362 / 0.2e1;
t756 = -t435 / 0.2e1;
t195 = -t405 * t753 + t408 * t756;
t797 = t617 * t195;
t795 = t184 + t198;
t296 = t493 * t325;
t334 = t743 * t347;
t581 = -t296 / 0.2e1 - t334 / 0.2e1;
t79 = t581 + t800;
t78 = t581 - t800;
t754 = -t439 / 0.2e1;
t760 = t405 / 0.2e1;
t793 = t158 * t754 - t760 * t785 + t807;
t792 = t617 * t785;
t759 = -t405 / 0.2e1;
t788 = t497 * t767;
t780 = t617 * t439;
t787 = t435 * t780;
t518 = -t671 / 0.2e1 - t599 / 0.2e1;
t281 = (t754 - t518) * t497;
t619 = t497 * qJD(1);
t377 = t408 * t619;
t576 = t617 * t408;
t124 = -t281 * qJD(2) + t377 - t576;
t645 = qJD(2) * t439;
t784 = -qJD(1) * t195 + t435 * t645;
t646 = qJD(1) * t408;
t783 = qJD(2) * t195 - t405 * t646;
t782 = t617 * t779;
t781 = t184 / 0.2e1;
t398 = t408 * qJ(5);
t157 = t183 + t398;
t151 = -pkin(4) * t497 - t157;
t660 = t151 + t157;
t455 = t493 * t668;
t410 = t497 * t598 - t455;
t685 = t410 * qJ(5);
t601 = -t685 / 0.2e1 + t655;
t776 = t321 * t756 + t393 * t759;
t489 = t494 ^ 2;
t491 = t496 ^ 2;
t466 = t491 - t489;
t621 = t495 * qJD(1);
t596 = t496 * t621;
t774 = qJD(2) * t466 - 0.2e1 * t494 * t596;
t773 = t408 ^ 2;
t772 = t439 ^ 2;
t771 = -pkin(4) / 0.2e1;
t191 = t297 - t674;
t729 = t495 * pkin(4);
t156 = t191 - t685 + t729;
t769 = t156 / 0.2e1;
t763 = t264 / 0.2e1;
t315 = -t326 / 0.2e1;
t761 = t779 / 0.2e1;
t755 = t435 / 0.2e1;
t752 = t444 / 0.2e1;
t611 = t743 * pkin(3);
t482 = t611 + pkin(4);
t751 = -t482 / 0.2e1;
t750 = t482 / 0.2e1;
t749 = -t493 / 0.2e1;
t748 = t493 / 0.2e1;
t747 = -t494 / 0.2e1;
t746 = -t496 / 0.2e1;
t741 = t163 * pkin(4);
t740 = t264 * pkin(4);
t407 = t439 * t497;
t738 = t407 * pkin(4);
t737 = t408 * pkin(4);
t736 = t410 * pkin(4);
t734 = t439 * pkin(4);
t733 = t493 * pkin(3);
t732 = t494 * pkin(3);
t488 = t497 * pkin(6);
t726 = pkin(3) * qJD(3);
t725 = pkin(4) * qJD(4);
t18 = t660 * t405;
t723 = qJD(1) * t18;
t231 = t321 * t408;
t21 = pkin(4) * t231 - t158 * t660;
t722 = qJD(1) * t21;
t47 = -t151 * t408 - t158 * t405;
t721 = qJD(1) * t47;
t614 = pkin(3) * t666;
t344 = t614 + t737;
t63 = t163 * t497 - t344 * t405 - t231;
t720 = qJD(1) * t63;
t164 = -t398 + t199;
t704 = t321 * t405;
t64 = t164 * t497 + t344 * t408 - t704;
t719 = qJD(1) * t64;
t65 = pkin(4) * t773 - t157 * t497 - t704;
t718 = qJD(1) * t65;
t66 = -t158 * t497 - t405 * t737 - t231;
t717 = qJD(1) * t66;
t192 = t334 + t296;
t687 = t407 * qJ(5);
t161 = t192 - t687;
t17 = -t151 * t410 - t156 * t408 - t158 * t407 - t161 * t405;
t714 = t17 * qJD(1);
t19 = -t805 * t408 + (t151 - t164) * t405;
t711 = t19 * qJD(1);
t472 = pkin(3) * t668;
t447 = t488 + t472;
t322 = t447 + t738;
t20 = t151 * t156 + t158 * t161 + t321 * t322;
t710 = t20 * qJD(1);
t22 = t151 * t163 + t158 * t164 + t321 * t344;
t709 = t22 * qJD(1);
t28 = t183 * t410 - t184 * t407 - t191 * t408 - t192 * t405;
t706 = t28 * qJD(1);
t30 = -t405 * t796 - t408 * t795;
t705 = t30 * qJD(1);
t35 = -t183 * t191 + t184 * t192 + t446 * t447;
t701 = t35 * qJD(1);
t36 = -t151 * t495 + t156 * t497 - t321 * t407 - t322 * t405;
t699 = t36 * qJD(1);
t383 = -t496 * t542 + t616;
t693 = t383 * t497;
t692 = t384 * t497;
t303 = t393 * t439;
t688 = t405 * t435;
t41 = -t183 * t198 + t184 * t199 + t446 * t614;
t686 = t41 * qJD(1);
t42 = -t158 * t495 + t161 * t497 + t321 * t410 + t322 * t408;
t684 = t42 * qJD(1);
t682 = t439 * t408;
t681 = t446 * t405;
t670 = t494 * t405;
t490 = t495 ^ 2;
t667 = t496 * t490;
t53 = t183 * t495 + t191 * t497 - t405 * t447 - t407 * t446;
t664 = t53 * qJD(1);
t54 = -t184 * t495 + t192 * t497 + t408 * t447 + t410 * t446;
t663 = t54 * qJD(1);
t612 = t743 / 0.2e1;
t562 = t497 * t612;
t525 = -t295 / 0.2e1 + pkin(3) * t562;
t91 = t600 / 0.2e1 + t525;
t662 = t91 * qJD(1);
t133 = -t682 + t688;
t661 = t617 * t133;
t421 = t435 * qJD(4);
t654 = -qJD(3) * t435 - t421;
t492 = t497 ^ 2;
t467 = t492 - t490;
t332 = t446 * t408;
t105 = t198 * t497 - t405 * t614 - t332;
t653 = qJD(1) * t105;
t106 = t199 * t497 + t408 * t614 - t681;
t652 = qJD(1) * t106;
t125 = -t183 * t497 - t681;
t651 = qJD(1) * t125;
t126 = -t184 * t497 - t332;
t650 = qJD(1) * t126;
t316 = t490 * t742 + t693;
t649 = qJD(1) * t316;
t317 = -pkin(6) * t667 - t692;
t648 = qJD(1) * t317;
t647 = qJD(1) * t405;
t644 = qJD(2) * t483;
t643 = qJD(2) * t494;
t642 = qJD(2) * t496;
t641 = qJD(3) * t494;
t640 = qJD(3) * t496;
t639 = qJD(3) * t497;
t638 = qJD(4) * t158;
t637 = qJD(4) * t264;
t636 = qJD(4) * t439;
t111 = (t394 * t495 - t693) * t496 + (t395 * t495 + t692) * t494;
t635 = t111 * qJD(1);
t152 = pkin(6) ^ 2 * t495 * t497 - t383 * t394 + t384 * t395;
t634 = t152 * qJD(1);
t173 = -t405 * t410 - t407 * t408;
t633 = t173 * qJD(1);
t215 = t395 * t497 + (-t384 + 0.2e1 * t613) * t495;
t632 = t215 * qJD(1);
t216 = t383 * t495 + (t394 - 0.2e1 * t474) * t497;
t631 = t216 * qJD(1);
t254 = t405 * t495 - t407 * t497;
t630 = t254 * qJD(1);
t255 = -t408 * t495 + t410 * t497;
t629 = t255 * qJD(1);
t628 = t281 * qJD(1);
t282 = (t754 + t518) * t497;
t272 = t282 * qJD(1);
t283 = t435 * t745 + t496 * t562 - t455 / 0.2e1;
t273 = t283 * qJD(1);
t284 = t455 / 0.2e1 + (t756 - t598 / 0.2e1) * t497;
t627 = t284 * qJD(1);
t626 = t405 * qJD(5);
t397 = t408 * qJD(5);
t425 = (t489 / 0.2e1 - t491 / 0.2e1) * t495;
t625 = t425 * qJD(3);
t438 = t467 * t494;
t624 = t438 * qJD(1);
t423 = t439 * qJD(5);
t441 = t492 * t496 - t667;
t623 = t441 * qJD(1);
t622 = t467 * qJD(1);
t620 = t495 * qJD(2);
t618 = t497 * qJD(2);
t610 = pkin(1) * t621;
t609 = pkin(1) * t619;
t608 = qJD(4) * t733;
t607 = t739 / 0.2e1;
t606 = t735 / 0.2e1;
t605 = t733 / 0.2e1;
t604 = t732 / 0.2e1;
t603 = t730 / 0.2e1;
t602 = t751 + t771;
t594 = t494 * t642;
t593 = t496 * t620;
t592 = t494 * t639;
t591 = t496 * t639;
t590 = t494 * t640;
t589 = t495 * t618;
t588 = t495 * t619;
t587 = t405 * t748;
t586 = t407 * t749;
t585 = t435 * t749;
t584 = -t666 / 0.2e1;
t583 = t157 / 0.2e1 + t151 / 0.2e1;
t580 = t743 * qJD(3);
t579 = t743 * qJD(4);
t202 = qJD(2) * t282 - t377;
t574 = t617 * t497;
t573 = t730 * t749;
t572 = t496 * t603;
t571 = -qJD(3) + t619;
t570 = pkin(3) * t579;
t568 = t151 * t755 + t793;
t565 = t494 * t593;
t564 = t490 * t590;
t560 = t619 - qJD(3) / 0.2e1;
t399 = t732 + t734;
t499 = t151 * t768 + t164 * t763 + t321 * t399 / 0.2e1 + t344 * t393 / 0.2e1 + t805 * t799;
t523 = t156 * t750 + t161 * t605;
t1 = -t499 + t523;
t52 = t393 * t399;
t558 = -qJD(1) * t1 + qJD(2) * t52;
t3 = t583 * t264 + t806 + (t769 - t657) * pkin(4);
t59 = pkin(4) * t303;
t557 = -qJD(1) * t3 + qJD(2) * t59;
t500 = t157 * t755 + t568 - t793;
t8 = t736 / 0.2e1 + t500;
t556 = qJD(1) * t8;
t522 = pkin(3) * t586 + t410 * t751;
t5 = (t158 / 0.2e1 + t163 / 0.2e1) * t439 + (t164 / 0.2e1 - t151 / 0.2e1) * t435 + (t763 + t768) * t408 + t802 * t405 + t522;
t555 = t5 * qJD(1);
t511 = (t586 - t743 * t410 / 0.2e1) * pkin(3);
t9 = (t781 + t198 / 0.2e1) * t439 + (t199 / 0.2e1 + t183 / 0.2e1) * t435 + (t362 / 0.2e1 + t762) * t408 + (t789 + t761) * t405 + t511;
t554 = t9 * qJD(1);
t552 = -t394 * t494 + t395 * t496;
t107 = t483 * t732;
t506 = t198 * t761 + t762 * t796 + t779 * t781;
t519 = t191 * t612 + t192 * t748;
t11 = (t446 * t747 + t483 * t584 + t519) * pkin(3) + t506;
t551 = -t11 * qJD(1) + t107 * qJD(2);
t110 = -t264 * t435 + t439 * t785;
t507 = t151 * t753 + t158 * t755 + t264 * t760 - t757 * t785;
t534 = t472 / 0.2e1 + t488 / 0.2e1 + t738 / 0.2e1;
t25 = t507 + t534;
t550 = -qJD(1) * t25 + qJD(2) * t110;
t189 = t399 * t435 + t303;
t508 = t344 * t755 + t399 * t760 - t801;
t32 = t495 * t602 + t508 - t601;
t549 = -qJD(1) * t32 - qJD(2) * t189;
t302 = t393 * t435;
t190 = t399 * t439 - t302;
t530 = t687 / 0.2e1 + t581;
t505 = t530 - t776;
t510 = t344 * t753 + t399 * t757 - t788;
t33 = t573 + t505 - t510;
t548 = qJD(1) * t33 - qJD(2) * t190;
t203 = -t435 * t734 - t303;
t37 = (t495 + t195) * pkin(4) + t601 + t801;
t547 = qJD(1) * t37 + qJD(2) * t203;
t204 = -pkin(4) * t772 + t302;
t526 = -pkin(4) * t682 + t788;
t39 = t505 + t526;
t546 = qJD(1) * t39 + qJD(2) * t204;
t277 = t435 * t732 + t439 * t483;
t56 = (-t670 / 0.2e1 + (t435 * t746 + t612) * t495) * pkin(3) + t77;
t545 = qJD(1) * t56 - qJD(2) * t277;
t278 = -t435 * t483 + t439 * t732;
t55 = (t408 * t747 + (t439 * t746 + t749) * t495) * pkin(3) + t78;
t544 = qJD(1) * t55 - qJD(2) * t278;
t346 = t752 - t444 / 0.2e1;
t498 = -t501 / 0.2e1 + t315 + t497 * t605;
t89 = t676 / 0.2e1 + t326 / 0.2e1 + t498;
t543 = qJD(1) * t89 + qJD(2) * t346;
t541 = t571 * t495;
t401 = t405 ^ 2;
t217 = t401 - t773;
t75 = qJD(1) * t217 + qJD(2) * t133;
t434 = t435 ^ 2;
t265 = t434 - t772;
t87 = qJD(1) * t133 + qJD(2) * t265;
t172 = t602 * t408 + (-t587 + t584) * pkin(3);
t214 = t602 * t439 + (t585 + t747) * pkin(3);
t540 = qJD(1) * t172 + qJD(2) * t214;
t536 = -t611 / 0.2e1 + t750;
t524 = t771 + t536;
t180 = t524 * t405;
t225 = t524 * t435;
t539 = qJD(1) * t180 + qJD(2) * t225;
t197 = t682 + t688;
t274 = t401 + t773;
t538 = qJD(1) * t274 + qJD(2) * t197;
t343 = t434 + t772;
t537 = qJD(1) * t197 + qJD(2) * t343;
t292 = -qJD(2) * t435 - t647;
t293 = t645 + t646;
t535 = t731 / 0.2e1 - t727 / 0.2e1;
t520 = t535 * t494;
t339 = t442 / 0.2e1 + t520;
t533 = pkin(2) * t642 - qJD(1) * t339;
t521 = t535 * t496;
t340 = -t445 / 0.2e1 - t521;
t532 = pkin(2) * t643 - qJD(1) * t340;
t529 = qJD(1) * t78 + t435 * t644;
t528 = qJD(1) * t77 - t439 * t644;
t527 = t496 * t541;
t367 = -qJD(1) * t425 + t594;
t348 = qJD(1) * t494 * t667 + qJD(2) * t425;
t437 = t466 * t490;
t517 = qJD(1) * t437 + 0.2e1 * t565;
t503 = (t158 * t612 - t493 * t583) * pkin(3) + t158 * t751;
t16 = -t741 / 0.2e1 + t503;
t411 = (t611 - t482) * t733;
t502 = (t264 * t612 - t493 * t802) * pkin(3) + t264 * t751;
t46 = t740 / 0.2e1 + t502;
t515 = -qJD(1) * t16 - qJD(2) * t46 - qJD(3) * t411;
t512 = t530 + t776;
t269 = 0.2e1 * t752 - t569;
t92 = t675 - t600 / 0.2e1 + t525;
t90 = -t676 / 0.2e1 + t315 + t498;
t479 = -t621 / 0.2e1;
t478 = t621 / 0.2e1;
t477 = t620 / 0.2e1;
t433 = t560 * t495;
t420 = t435 * qJD(5);
t415 = t439 * t733;
t400 = (-qJD(4) / 0.2e1 + t560) * t495;
t378 = t408 * t733;
t338 = t346 * qJD(3);
t337 = t346 * qJD(4);
t286 = t474 + t445 / 0.2e1 - t521;
t285 = t615 - t442 / 0.2e1 + t520;
t276 = t293 * pkin(4);
t224 = t435 * t536 + t606;
t223 = t269 + t683;
t213 = pkin(3) * t585 + t439 * t751 + t604 + t734 / 0.2e1;
t201 = qJD(2) * t283 - t405 * t619;
t193 = t197 * qJD(5);
t182 = -t780 - t272;
t181 = -t273 + t654;
t179 = t405 * t536 + t607;
t171 = -pkin(3) * t587 + t408 * t751 + t572 + t737 / 0.2e1;
t162 = t617 * t265;
t150 = -t255 * qJD(2) + t405 * t574;
t149 = -t254 * qJD(2) + t408 * t574;
t123 = -t284 * qJD(2) + (-qJD(4) + t571) * t405;
t109 = (qJD(2) * t410 - t405 * t617) * t408;
t108 = (qJD(2) * t407 + t576) * t405;
t104 = -t283 * t617 + t629;
t103 = -t282 * t617 + t630;
t94 = -t284 * t617 + t439 * t620 - t629;
t93 = -t281 * t617 - t435 * t620 - t630;
t86 = t91 * qJD(3);
t85 = t89 * qJD(3);
t84 = t91 * qJD(4);
t83 = t89 * qJD(4);
t72 = t398 + t92;
t71 = t90 + t689;
t70 = -t410 * t646 + t797;
t69 = -t407 * t647 - t797;
t68 = pkin(3) * t580 - t662;
t67 = t493 * t726 - t543;
t62 = t662 + (-t580 - t579) * pkin(3);
t61 = -t617 * t733 + t543;
t60 = t173 * qJD(2) + t217 * t617;
t58 = t408 * t604 + t439 * t572 + t573 + t79;
t57 = pkin(3) * t670 / 0.2e1 + t435 * t572 + t743 * t603 + t80;
t51 = t293 * t410 + t797;
t50 = -t292 * t407 - t797;
t45 = -t740 / 0.2e1 + t502;
t43 = -t633 + t661;
t40 = t512 - t526;
t38 = t408 * t606 + t439 * t607 + t601 + t729 - t801;
t34 = t573 + t510 + t512;
t31 = t495 * t750 + t729 / 0.2e1 + t508 + t601;
t29 = t633 + (-t407 * t439 - t410 * t435) * qJD(2) + t661;
t26 = -t507 + t534;
t15 = t741 / 0.2e1 + t503;
t12 = pkin(3) * t519 + t446 * t604 + t483 * t572 - t506;
t10 = t405 * t789 + t754 * t795 + t756 * t796 - t779 * t759 + t511;
t7 = -t736 / 0.2e1 + t500;
t6 = t163 * t754 + t164 * t756 - t759 * t785 + t522 + t568 - t807;
t4 = t660 * t768 + t806 + (t769 + t657) * pkin(4);
t2 = t499 + t523;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t589, t467 * qJD(2), 0, -t589, 0, 0, -pkin(1) * t620, -pkin(1) * t618, 0, 0, t491 * t589 - t564, -qJD(3) * t437 - 0.2e1 * t497 * t565, -qJD(2) * t441 + t495 * t592, t489 * t589 + t564, qJD(2) * t438 + t495 * t591, -t589, -qJD(2) * t216 - qJD(3) * t317, qJD(2) * t215 - qJD(3) * t316, -qJD(2) * t111, qJD(2) * t152, t109, t60, t150, t108, t149, -t589, -qJD(2) * t53 - qJD(3) * t105 - qJD(4) * t126, qJD(2) * t54 + qJD(3) * t106 + qJD(4) * t125, qJD(2) * t28 + qJD(3) * t30, qJD(2) * t35 + qJD(3) * t41, t109, t60, t150, t108, t149, -t589, -qJD(2) * t36 - qJD(3) * t63 - qJD(4) * t66 + t397 * t497, qJD(2) * t42 + qJD(3) * t64 + qJD(4) * t65 - t497 * t626, qJD(2) * t17 + qJD(3) * t19 + qJD(4) * t18 + qJD(5) * t274, qJD(2) * t20 + qJD(3) * t22 + qJD(4) * t21 + qJD(5) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t588, t622, t618, -t588, -t620, 0, -pkin(6) * t618 - t610, pkin(6) * t620 - t609, 0, 0, -t625 + (t491 * t621 + t594) * t497, -0.2e1 * t495 * t590 + t497 * t774, t494 * t620 - t623, t625 + (t489 * t621 - t594) * t497, t593 + t624, -t433, -t631 + (t494 * t559 - t613) * qJD(2) + t286 * qJD(3), t632 + (t496 * t559 + t616) * qJD(2) + t285 * qJD(3), qJD(2) * t552 - t635, t634 + (-pkin(2) * t488 + pkin(7) * t552) * qJD(2), t51, t29, t94, t50, t93, -t400, -t664 + (t407 * t483 + t435 * t447 - t495 * t779) * qJD(2) + t57 * qJD(3) + t80 * qJD(4), t663 + (-t362 * t495 + t410 * t483 + t439 * t447) * qJD(2) + t58 * qJD(3) + t79 * qJD(4), t706 + (-t191 * t439 - t192 * t435 - t362 * t407 + t410 * t779) * qJD(2) + t10 * qJD(3), t701 + (-t191 * t779 + t192 * t362 + t447 * t483) * qJD(2) + t12 * qJD(3), t51, t29, t94, t50, t93, -t400, -t699 + (t322 * t435 + t393 * t407 - t495 * t785) * qJD(2) + t31 * qJD(3) + t38 * qJD(4) - t282 * qJD(5), t684 + (-t264 * t495 + t322 * t439 + t393 * t410) * qJD(2) + t34 * qJD(3) + t40 * qJD(4) + t283 * qJD(5), t714 + (-t156 * t439 - t161 * t435 - t264 * t407 + t410 * t785) * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + t193, t710 + (-t156 * t785 + t161 * t264 + t322 * t393) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t26 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t348, -t517, t494 * t541, t348, t527, t477, qJD(2) * t286 - qJD(3) * t384 - t648, qJD(2) * t285 + qJD(3) * t383 - t649, 0, 0, t783, t75, t123, -t783, t124, t477, qJD(2) * t57 + qJD(3) * t198 + qJD(4) * t90 - t653, qJD(2) * t58 - qJD(3) * t199 + qJD(4) * t92 + t652, t705 + t10 * qJD(2) + (t405 * t611 - t378) * qJD(3), t686 + t12 * qJD(2) + (t198 * t743 + t199 * t493) * t726, t783, t75, t123, -t783, t124, t477, qJD(2) * t31 + qJD(3) * t163 + qJD(4) * t71 - t720, qJD(2) * t34 - qJD(3) * t164 + qJD(4) * t72 + t719, t711 + t6 * qJD(2) + (t482 * t405 - t378) * qJD(3) + t179 * qJD(4), t709 + t2 * qJD(2) + (t163 * t482 + t164 * t733) * qJD(3) + t15 * qJD(4) + t171 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t783, t75, t123, -t783, t124, t477, qJD(2) * t80 + qJD(3) * t90 - qJD(4) * t184 - t650, qJD(2) * t79 + qJD(3) * t92 + qJD(4) * t183 + t651, 0, 0, t783, t75, t123, -t783, t124, t477, qJD(2) * t38 + qJD(3) * t71 - t638 - t717, qJD(2) * t40 + qJD(3) * t72 + qJD(4) * t157 + t718, qJD(2) * t7 + qJD(3) * t179 + t405 * t725 + t723, -pkin(4) * t638 + qJD(2) * t4 + qJD(3) * t15 + t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t201, t538, qJD(2) * t26 + qJD(3) * t171 + t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, -t622, 0, t588, 0, 0, t610, t609, 0, 0, -t491 * t588 - t625, 0.2e1 * t494 * t527, -t591 + t623, -t489 * t588 + t625, t592 - t624, t433, qJD(3) * t340 + t631, qJD(3) * t339 - t632, t635, -t634, t70, t43, t104, t69, t103, t400, -qJD(3) * t56 - qJD(4) * t77 + t664, -qJD(3) * t55 - qJD(4) * t78 - t663, -qJD(3) * t9 - t706, -qJD(3) * t11 - t701, t70, t43, t104, t69, t103, t400, qJD(3) * t32 - qJD(4) * t37 - qJD(5) * t281 + t699, -qJD(3) * t33 - qJD(4) * t39 + qJD(5) * t284 - t684, -qJD(3) * t5 + qJD(4) * t8 + t193 - t714, -qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t25 - t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t590, t466 * qJD(3), 0, -t590, 0, 0, -pkin(2) * t641, -pkin(2) * t640, 0, 0, -t787, t162, 0, t787, 0, 0, qJD(3) * t277 + t483 * t636, qJD(3) * t278 - t421 * t483, 0, qJD(3) * t107, -t787, t162, 0, t787, 0, 0, qJD(3) * t189 - qJD(4) * t203, qJD(3) * t190 - qJD(4) * t204, qJD(5) * t343, qJD(3) * t52 + qJD(4) * t59 + qJD(5) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t367, t774, -t571 * t496, -t367, t571 * t494, t479, -pkin(7) * t640 - t532, pkin(7) * t641 - t533, 0, 0, -t784, t87, t181, t784, t182, t479, -qJD(3) * t362 + qJD(4) * t269 - t545, -t544 + t782, (t435 * t611 - t415) * qJD(3) - t554, (-t362 * t743 - t493 * t779) * t726 + t551, -t784, t87, t181, t784, t182, t479, -qJD(3) * t264 + qJD(4) * t223 - t549, -t548 + t792, (t482 * t435 - t415) * qJD(3) + t224 * qJD(4) - t555, (-t264 * t482 - t733 * t785) * qJD(3) + t45 * qJD(4) + t213 * qJD(5) + t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t784, t87, t181, t784, t182, t479, qJD(3) * t269 - qJD(4) * t362 - t528, -t529 + t782, 0, 0, -t784, t87, t181, t784, t182, t479, qJD(3) * t223 - t547 - t637, -t546 + t792, pkin(4) * t421 + qJD(3) * t224 + t556, -pkin(4) * t637 + qJD(3) * t45 + t557; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t628, t627, t537, qJD(3) * t213 + t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, t517, (-t494 * t621 + t642) * t497, -t348, (-t596 - t643) * t497, t477, -qJD(2) * t340 + t648, -qJD(2) * t339 + t649, 0, 0, -t783, -t75, t201, t783, t202, t477, qJD(2) * t56 + t653 + t83, qJD(2) * t55 - t652 + t84, qJD(2) * t9 - t705, qJD(2) * t11 - t686, -t783, -t75, t201, t783, t202, t477, -qJD(2) * t32 - t397 + t720 + t83, qJD(2) * t33 + t626 - t719 + t84, qJD(2) * t5 + qJD(4) * t180 - t711, qJD(2) * t1 + qJD(4) * t16 + qJD(5) * t172 - t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t367, -t774, t496 * t619, t367, -t494 * t619, t478, t532, t533, 0, 0, t784, -t87, t273, -t784, t272, t478, t337 + t545, t544, t554, -t551, t784, -t87, t273, -t784, t272, t478, t337 - t423 + t549, t420 + t548, qJD(4) * t225 + t555, qJD(4) * t46 + qJD(5) * t214 - t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t608, -t570, 0, 0, 0, 0, 0, 0, 0, 0, -t608, -t570, 0, t411 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, t539, -pkin(4) * t608 - t515; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, -t292, 0, t540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t783, -t75, t201, t783, t202, t477, qJD(2) * t77 + t650 - t85, qJD(2) * t78 - t651 - t86, 0, 0, -t783, -t75, t201, t783, t202, t477, qJD(2) * t37 - t397 + t717 - t85, qJD(2) * t39 + t626 - t718 - t86, -qJD(2) * t8 - qJD(3) * t180 - t723, -pkin(4) * t397 + qJD(2) * t3 - qJD(3) * t16 - t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t784, -t87, t273, -t784, t272, t478, -t338 + t528, t529, 0, 0, t784, -t87, t273, -t784, t272, t478, -t338 - t423 + t547, t420 + t546, -qJD(3) * t225 - t556, -pkin(4) * t423 - qJD(3) * t46 - t557; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t68, 0, 0, 0, 0, 0, 0, 0, 0, t67, t68, -t539, t515; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, -t292, 0, -t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t123, -t538, qJD(2) * t25 - qJD(3) * t172 + t408 * t725 - t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t780 + t628, -t627 + t654, -t537, pkin(4) * t636 - qJD(3) * t214 - t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, t292, 0, -t540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, t292, 0, t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;
