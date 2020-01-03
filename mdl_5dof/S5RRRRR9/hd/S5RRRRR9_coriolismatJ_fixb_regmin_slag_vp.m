% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x31]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:56
% EndTime: 2019-12-31 22:30:39
% DurationCPUTime: 23.53s
% Computational Cost: add. (10158->615), mult. (22758->843), div. (0->0), fcn. (24461->8), ass. (0->499)
t803 = qJD(4) + qJD(5);
t812 = qJD(3) + t803;
t450 = cos(qJ(3));
t752 = pkin(7) + pkin(8);
t400 = t752 * t450;
t449 = cos(qJ(4));
t388 = t449 * t400;
t446 = sin(qJ(3));
t399 = t752 * t446;
t445 = sin(qJ(4));
t670 = t445 * t399;
t316 = t388 - t670;
t649 = t449 * t450;
t669 = t445 * t446;
t380 = -t649 + t669;
t245 = -pkin(9) * t380 + t316;
t444 = sin(qJ(5));
t448 = cos(qJ(5));
t650 = t449 * t446;
t668 = t445 * t450;
t516 = t650 + t668;
t548 = t399 * t449 + t400 * t445;
t756 = -pkin(9) * t516 - t548;
t141 = t448 * t245 + t444 * t756;
t809 = t812 * t141;
t806 = t444 * t245 - t448 * t756;
t808 = t812 * t806;
t451 = cos(qJ(2));
t735 = t451 / 0.2e1;
t807 = t141 * t735;
t736 = -t451 / 0.2e1;
t804 = t736 * t806;
t214 = t380 * t444 - t448 * t516;
t596 = qJD(3) + qJD(4);
t785 = qJD(5) + t596;
t353 = t516 * t451;
t338 = t448 * t353;
t354 = t451 * t380;
t676 = t444 * t354;
t528 = -t338 / 0.2e1 + t676 / 0.2e1;
t780 = t214 * t451;
t783 = t780 / 0.2e1 + t528;
t793 = qJD(1) * t783;
t48 = t214 * t785 - t793;
t549 = -t380 * t448 - t444 * t516;
t762 = t785 * t549;
t653 = t448 * t354;
t677 = t444 * t353;
t494 = -t677 / 0.2e1 - t653 / 0.2e1;
t766 = t549 * t451;
t777 = t766 / 0.2e1 + t494;
t788 = qJD(1) * t777;
t49 = -t788 + t762;
t784 = -t780 / 0.2e1 + t528;
t792 = qJD(2) * t784;
t447 = sin(qJ(2));
t667 = t446 * t447;
t350 = t445 * t667 - t447 * t649;
t351 = t516 * t447;
t190 = t350 * t448 + t351 * t444;
t598 = t451 * qJD(1);
t61 = qJD(2) * t783 + t190 * t598;
t550 = t350 * t444 - t351 * t448;
t776 = -t190 ^ 2 + t550 ^ 2;
t789 = qJD(1) * t776;
t775 = -t214 ^ 2 + t549 ^ 2;
t787 = qJD(2) * t775;
t778 = -t766 / 0.2e1 + t494;
t786 = t778 * qJD(2);
t62 = qJD(2) * t777 + t550 * t598;
t782 = t190 * t596;
t734 = pkin(3) * t446;
t577 = pkin(6) + t734;
t534 = t577 * t447;
t731 = pkin(4) * t351;
t475 = -t534 - t731;
t172 = t475 * t190;
t779 = t190 * qJD(5);
t750 = -t190 / 0.2e1;
t460 = t475 * t549;
t771 = -t460 / 0.2e1;
t462 = t475 * t550;
t619 = qJD(2) * t214;
t765 = t549 * t619;
t623 = qJD(1) * t190;
t763 = t550 * t623;
t738 = -t516 / 0.2e1;
t759 = t550 * t596;
t461 = t475 * t214;
t733 = pkin(3) * t450;
t436 = -pkin(2) - t733;
t343 = pkin(4) * t380 + t436;
t697 = t343 * t190;
t758 = t461 / 0.2e1 - t697 / 0.2e1;
t696 = t343 * t550;
t757 = t771 + t696 / 0.2e1;
t755 = -qJD(5) * t550 - t759;
t440 = t446 ^ 2;
t442 = t450 ^ 2;
t418 = t442 - t440;
t600 = t447 * qJD(1);
t576 = t450 * t600;
t754 = qJD(2) * t418 - 0.2e1 * t446 * t576;
t753 = t596 * t548;
t751 = t550 / 0.2e1;
t533 = -pkin(2) * t451 - pkin(7) * t447;
t394 = -pkin(1) + t533;
t379 = t450 * t394;
t665 = t447 * t450;
t539 = -pkin(8) * t665 + t379;
t727 = pkin(6) * t446;
t286 = (-pkin(3) - t727) * t451 + t539;
t266 = t449 * t286;
t747 = -t266 / 0.2e1;
t746 = t549 / 0.2e1;
t745 = t214 / 0.2e1;
t744 = -t286 / 0.2e1;
t647 = t450 * t451;
t591 = pkin(6) * t647;
t340 = t394 * t446 + t591;
t307 = -pkin(8) * t667 + t340;
t293 = t449 * t307;
t560 = -t293 / 0.2e1;
t730 = pkin(4) * t516;
t344 = t730 + t734;
t743 = -t344 / 0.2e1;
t742 = t350 / 0.2e1;
t741 = -t351 / 0.2e1;
t740 = -t380 / 0.2e1;
t739 = t380 / 0.2e1;
t559 = -t388 / 0.2e1;
t737 = -t447 / 0.2e1;
t732 = pkin(4) * t350;
t729 = pkin(4) * t447;
t728 = pkin(4) * t451;
t726 = pkin(9) * t351;
t724 = t350 * pkin(9);
t322 = t444 * pkin(4);
t723 = t445 * pkin(3);
t722 = t447 * pkin(2);
t721 = t447 * pkin(3);
t720 = t449 * pkin(3);
t719 = t451 * pkin(7);
t718 = pkin(3) * qJD(3);
t717 = pkin(3) * qJD(4);
t716 = pkin(4) * qJD(4);
t715 = pkin(4) * qJD(5);
t234 = t338 - t676;
t666 = t446 * t451;
t390 = pkin(3) * t666 + pkin(6) * t451;
t291 = pkin(4) * t353 + t390;
t672 = t445 * t307;
t181 = -t266 + t672;
t504 = -t181 + t724;
t136 = t504 - t728;
t664 = t448 * t136;
t182 = t445 * t286 + t293;
t146 = t182 - t726;
t687 = t444 * t146;
t57 = -t664 + t687;
t401 = -t719 + t722;
t389 = t450 * t401;
t426 = pkin(6) * t667;
t292 = -pkin(8) * t647 + t389 + t426 + t721;
t268 = t449 * t292;
t386 = t446 * t401;
t590 = pkin(6) * t665;
t310 = -pkin(8) * t666 + t386 - t590;
t671 = t445 * t310;
t551 = t268 - t671;
t139 = pkin(9) * t354 + t551 + t729;
t663 = t448 * t139;
t267 = t445 * t292;
t301 = t449 * t310;
t633 = t301 + t267;
t147 = -pkin(9) * t353 + t633;
t686 = t444 * t147;
t5 = (t663 - t686) * t451 + t57 * t447 + t291 * t550 + t475 * t234;
t714 = t5 * qJD(1);
t239 = -t653 - t677;
t662 = t448 * t146;
t689 = t444 * t136;
t58 = t662 + t689;
t661 = t448 * t147;
t688 = t444 * t139;
t6 = (t661 + t688) * t451 - t291 * t190 + t239 * t731 + (t239 * t577 - t58) * t447;
t713 = t6 * qJD(1);
t477 = t444 * t504;
t59 = -t477 - t662;
t31 = t451 * t59 - t550 * t732 - t172;
t712 = qJD(1) * t31;
t476 = t448 * t504;
t60 = t476 - t687;
t32 = t190 * t732 + t60 * t451 - t462;
t711 = qJD(1) * t32;
t593 = pkin(3) * t665;
t510 = t593 - t732;
t592 = pkin(6) * t666;
t306 = t539 - t592;
t673 = t445 * t306;
t198 = -t293 - t673;
t152 = t198 + t726;
t660 = t448 * t152;
t651 = t449 * t306;
t199 = t651 - t672;
t153 = t199 + t724;
t684 = t444 * t153;
t77 = t660 - t684;
t33 = t451 * t77 + t510 * t550 - t172;
t710 = qJD(1) * t33;
t659 = t448 * t153;
t685 = t444 * t152;
t78 = t659 + t685;
t34 = -t190 * t510 + t78 * t451 - t462;
t709 = qJD(1) * t34;
t40 = -t57 * t451 - t462;
t708 = qJD(1) * t40;
t41 = -t451 * t58 - t172;
t707 = qJD(1) * t41;
t579 = pkin(3) * t735;
t537 = t579 + t306 / 0.2e1;
t94 = t449 * t537 + t747;
t706 = qJD(1) * t94;
t299 = t350 * t534;
t98 = t198 * t451 - t351 * t593 + t299;
t705 = qJD(1) * t98;
t99 = t199 * t451 + (-t350 * t733 - t351 * t577) * t447;
t704 = qJD(1) * t99;
t699 = t548 * t451;
t698 = t316 * t451;
t578 = pkin(4) + t720;
t402 = t448 * t578;
t422 = t444 * t723;
t361 = t422 - t402;
t695 = t361 * t451;
t535 = t444 * t578;
t652 = t448 * t445;
t362 = pkin(3) * t652 + t535;
t694 = t362 * t451;
t674 = t444 * t449;
t368 = (t652 + t674) * pkin(3);
t693 = t368 * t451;
t369 = t448 * t720 - t422;
t692 = t369 * t451;
t691 = t436 * t350;
t690 = t436 * t351;
t441 = t447 ^ 2;
t648 = t450 * t441;
t50 = t181 * t447 - t351 * t390 - t353 * t534 + t451 * t551;
t646 = t50 * qJD(1);
t51 = t633 * t451 - t390 * t350 + (-t354 * t577 - t182) * t447;
t645 = t51 * qJD(1);
t76 = t190 * t234 + t239 * t550;
t644 = t76 * qJD(1);
t115 = t350 * t516 + t351 * t380;
t643 = t596 * t115;
t563 = -t662 / 0.2e1;
t642 = -t689 / 0.2e1 + t563;
t566 = t687 / 0.2e1;
t641 = -t664 / 0.2e1 + t566;
t640 = -t477 / 0.2e1 + t563;
t639 = t566 - t476 / 0.2e1;
t195 = t350 * t739 + t351 * t738;
t638 = t596 * t195;
t443 = t451 ^ 2;
t419 = t443 - t441;
t108 = -t181 * t451 - t351 * t534;
t632 = qJD(1) * t108;
t109 = -t182 * t451 + t299;
t631 = qJD(1) * t109;
t120 = -t234 * t451 - t447 * t550;
t630 = qJD(1) * t120;
t121 = t190 * t447 + t239 * t451;
t629 = qJD(1) * t121;
t235 = -t351 * t447 + t353 * t451;
t625 = qJD(1) * t235;
t236 = -t350 * t447 + t354 * t451;
t624 = qJD(1) * t236;
t339 = -t379 + t592;
t287 = -t339 * t451 - t441 * t727;
t622 = qJD(1) * t287;
t288 = -pkin(6) * t648 - t340 * t451;
t621 = qJD(1) * t288;
t620 = qJD(1) * t350;
t618 = qJD(2) * t343;
t617 = qJD(2) * t516;
t616 = qJD(2) * t436;
t615 = qJD(2) * t446;
t614 = qJD(2) * t450;
t613 = qJD(3) * t446;
t612 = qJD(3) * t450;
t611 = qJD(3) * t451;
t610 = qJD(4) * t436;
t608 = qJD(5) * t343;
t176 = t350 * t353 + t351 * t354;
t607 = t176 * qJD(1);
t206 = t339 * t447 + (-t426 + t389) * t451;
t606 = t206 * qJD(1);
t207 = t386 * t451 + (-t340 + t591) * t447;
t605 = t207 * qJD(1);
t493 = -t668 / 0.2e1 - t650 / 0.2e1;
t258 = (t738 + t493) * t451;
t248 = t258 * qJD(1);
t492 = t649 / 0.2e1 - t669 / 0.2e1;
t259 = (t740 + t492) * t451;
t249 = t259 * qJD(1);
t384 = t419 * t446;
t603 = t384 * qJD(1);
t385 = t443 * t450 - t648;
t602 = t385 * qJD(1);
t601 = t419 * qJD(1);
t599 = t447 * qJD(2);
t597 = t451 * qJD(2);
t37 = -t190 * t745 + t214 * t750 + t549 * t751 + t550 * t746;
t45 = -t190 * t214 + t549 * t550;
t595 = t45 * qJD(5) + t37 * t596;
t84 = -t550 * t214 / 0.2e1 + t549 * t750;
t85 = -t190 * t746 - t550 * t745;
t594 = t85 * qJD(5) + t596 * t84;
t589 = pkin(1) * t600;
t588 = pkin(1) * t598;
t587 = t734 / 0.2e1;
t586 = t733 / 0.2e1;
t585 = -t732 / 0.2e1;
t584 = t730 / 0.2e1;
t583 = t728 / 0.2e1;
t581 = -t720 / 0.2e1;
t580 = t720 / 0.2e1;
t575 = t446 * t614;
t574 = t450 * t599;
t573 = t446 * t611;
t572 = t450 * t611;
t571 = t446 * t612;
t570 = t447 * t597;
t569 = t447 * t598;
t568 = t697 / 0.2e1;
t567 = -t696 / 0.2e1;
t565 = -t686 / 0.2e1;
t564 = -t684 / 0.2e1;
t562 = -t661 / 0.2e1;
t561 = -t659 / 0.2e1;
t558 = t152 / 0.2e1 + t146 / 0.2e1;
t557 = -t153 / 0.2e1 + t136 / 0.2e1;
t554 = -t267 / 0.2e1 - t301 / 0.2e1;
t553 = pkin(3) * t596;
t552 = pkin(4) * t803;
t545 = t596 * t351;
t544 = t596 * t516;
t543 = -qJD(3) + t598;
t542 = -qJD(5) + t598;
t541 = t446 * t574;
t538 = t729 / 0.2e1 + t139 / 0.2e1;
t536 = t598 - qJD(3) / 0.2e1;
t531 = t558 * t444;
t530 = t558 * t448;
t529 = t268 / 0.2e1 - t671 / 0.2e1;
t527 = -qJD(4) + t543;
t509 = -t534 / 0.2e1;
t526 = t380 * t509 + t554 - t690 / 0.2e1;
t453 = t510 * t746 - t550 * t743 - t807;
t499 = t565 + t663 / 0.2e1;
t459 = t361 * t737 + t499;
t1 = t568 - t461 / 0.2e1 + t453 + t459;
t211 = t343 * t214;
t102 = -t344 * t549 - t211;
t525 = -qJD(1) * t1 + qJD(2) * t102;
t212 = t343 * t549;
t103 = -t214 * t344 + t212;
t452 = -t190 * t743 + t510 * t745 - t804;
t500 = -t688 / 0.2e1 + t562;
t458 = t362 * t737 + t500;
t2 = t567 + t460 / 0.2e1 + t452 + t458;
t524 = -qJD(1) * t2 + qJD(2) * t103;
t110 = t549 * t730 + t211;
t481 = t141 * t451;
t8 = -t481 / 0.2e1 + t568 - t214 * t509 + (t448 * t447 / 0.2e1 - t550 * t738 - t549 * t742 - t214 * t741) * pkin(4) + t499;
t523 = -qJD(1) * t8 - qJD(2) * t110;
t111 = t214 * t730 - t212;
t482 = t806 * t451;
t7 = t482 / 0.2e1 + t567 + t549 * t509 + (-t190 * t738 - t214 * t742 + t444 * t737 + t549 * t741) * pkin(4) + t500;
t522 = -qJD(1) * t7 - qJD(2) * t111;
t11 = qJD(2) * t37 + t789;
t13 = qJD(1) * t37 + t787;
t521 = qJD(2) * t45 + t789;
t515 = qJD(1) * t45 + t787;
t254 = t380 * t734 + t436 * t516;
t456 = -t316 * t735 - t351 * t587;
t471 = t691 / 0.2e1 + t529;
t53 = (t577 * t738 + t733 * t740 + t580) * t447 + t456 + t471;
t514 = -qJD(1) * t53 + qJD(2) * t254;
t255 = -t380 * t436 + t516 * t734;
t457 = t350 * t587 + t548 * t735;
t491 = t690 / 0.2e1 + t554;
t52 = (-t723 / 0.2e1 - t516 * t586 + t577 * t739) * t447 + t457 + t491;
t513 = qJD(1) * t52 - qJD(2) * t255;
t309 = t559 + t388 / 0.2e1;
t92 = t560 + t293 / 0.2e1 + (t744 + t537) * t445;
t512 = qJD(1) * t92 + qJD(2) * t309;
t511 = t543 * t447;
t208 = -t350 ^ 2 + t351 ^ 2;
t82 = qJD(1) * t208 + qJD(2) * t115;
t240 = t380 ^ 2 - t516 ^ 2;
t90 = qJD(1) * t115 + qJD(2) * t240;
t508 = t534 / 0.2e1;
t507 = t719 / 0.2e1 - t722 / 0.2e1;
t506 = t444 * t583 + t642;
t505 = t448 * t583 + t641;
t503 = -qJD(4) / 0.2e1 + t536;
t479 = t507 * t446;
t304 = t386 / 0.2e1 - t479;
t502 = pkin(2) * t614 - qJD(1) * t304;
t478 = t507 * t450;
t305 = -t389 / 0.2e1 + t478;
t501 = pkin(2) * t615 - qJD(1) * t305;
t498 = -t685 / 0.2e1 + t561;
t497 = t564 + t660 / 0.2e1;
t454 = t343 * t750 + t475 * t745 + t807;
t15 = -t454 + t499;
t490 = qJD(1) * t15 + t214 * t618;
t455 = t343 * t751 + t771 + t804;
t16 = -t455 + t500;
t489 = qJD(1) * t16 - t549 * t618;
t488 = -qJD(2) * t85 + t763;
t42 = qJD(2) * t84 - t763;
t487 = -qJD(1) * t85 + t765;
t46 = qJD(1) * t84 - t765;
t87 = t699 / 0.2e1 + t380 * t508 + t491;
t486 = qJD(1) * t87 + t380 * t616;
t86 = -t698 / 0.2e1 - t516 * t508 + t471;
t485 = qJD(1) * t86 - t516 * t616;
t484 = t450 * t511;
t105 = -qJD(2) * t195 - t351 * t620;
t118 = qJD(1) * t195 - t380 * t617;
t371 = (t440 / 0.2e1 - t442 / 0.2e1) * t447;
t483 = -qJD(1) * t371 + t575;
t480 = -t516 * t509 + t529 - t691 / 0.2e1;
t474 = qJD(1) * t446 * t648 + qJD(2) * t371;
t383 = t418 * t441;
t473 = qJD(1) * t383 + 0.2e1 * t541;
t19 = -t694 / 0.2e1 + t530 + t557 * t444;
t470 = qJD(1) * t19 + qJD(3) * t362;
t20 = t695 / 0.2e1 + t557 * t448 - t531;
t469 = qJD(1) * t20 - qJD(3) * t361;
t464 = t477 / 0.2e1;
t23 = t464 + t662 / 0.2e1 + t506;
t468 = qJD(1) * t23 - qJD(3) * t322;
t463 = t476 / 0.2e1;
t25 = t463 - t687 / 0.2e1 + t505;
t323 = t402 / 0.2e1 + (t581 + pkin(4) / 0.2e1) * t448;
t467 = qJD(1) * t25 - qJD(3) * t323;
t27 = t564 - t693 / 0.2e1 + t464 + t530;
t466 = qJD(1) * t27 + qJD(3) * t368;
t28 = t561 - t692 / 0.2e1 + t463 - t531;
t465 = qJD(1) * t28 + qJD(3) * t369;
t431 = -t600 / 0.2e1;
t430 = t600 / 0.2e1;
t429 = t599 / 0.2e1;
t378 = t536 * t447;
t365 = t371 * qJD(3);
t364 = t369 * qJD(4);
t363 = t368 * qJD(4);
t349 = t362 * qJD(5);
t348 = t361 * qJD(5);
t347 = t503 * t447;
t318 = (-qJD(5) / 0.2e1 + t503) * t447;
t303 = t422 - t402 / 0.2e1 + (-pkin(4) / 0.2e1 + t581) * t448;
t302 = -t322 / 0.2e1 - t535 / 0.2e1 + (-t652 - t674 / 0.2e1) * pkin(3);
t263 = t426 + t389 / 0.2e1 + t478;
t262 = t590 - t386 / 0.2e1 - t479;
t261 = t451 * t493 - t516 * t736;
t260 = t380 * t735 + t451 * t492;
t247 = 0.2e1 * t559 + t670;
t203 = qJD(2) * t259 - t351 * t598;
t202 = qJD(2) * t258 + t350 * t598;
t180 = -t544 - t248;
t179 = -t380 * t596 - t249;
t107 = t260 * qJD(2) + t351 * t527;
t106 = t261 * qJD(2) - t350 * t527;
t95 = t449 * t579 + t672 + t747 - t651 / 0.2e1;
t93 = 0.2e1 * t560 - t673 / 0.2e1 + (t579 + t744) * t445;
t89 = t698 / 0.2e1 + t480;
t88 = -t699 / 0.2e1 + t526;
t55 = -t457 + t526 - (-t450 * t516 + t445) * t721 / 0.2e1;
t54 = -t456 + t480 + (t380 * t586 + t580) * t447;
t39 = t786 + (qJD(5) - t527) * t550;
t38 = -t190 * t527 + t779 + t792;
t30 = t692 / 0.2e1 + t498 + t639;
t29 = t693 / 0.2e1 + t497 + t640;
t26 = t505 + t639;
t24 = t506 + t640;
t22 = t694 / 0.2e1 + t497 + t642;
t21 = -t695 / 0.2e1 + t498 + t641;
t18 = t454 + t499;
t17 = t455 + t500;
t10 = -t482 / 0.2e1 - t190 * t584 - t214 * t585 + t562 - t538 * t444 + t757;
t9 = t481 / 0.2e1 - t550 * t584 - t549 * t585 + t565 + t538 * t448 + t758;
t4 = -t452 + t458 + t757;
t3 = -t453 + t459 + t758;
t12 = [0, 0, 0, t570, t419 * qJD(2), 0, 0, 0, -pkin(1) * t599, -pkin(1) * t597, -t441 * t571 + t442 * t570, -qJD(3) * t383 - 0.2e1 * t451 * t541, -qJD(2) * t385 + t447 * t573, qJD(2) * t384 + t447 * t572, -t570, -qJD(2) * t206 - qJD(3) * t288, qJD(2) * t207 + qJD(3) * t287, (qJD(2) * t354 + t545) * t350, qJD(2) * t176 + t208 * t596, t236 * qJD(2) + t451 * t545, -t350 * t451 * t596 + t235 * qJD(2), -t570, -qJD(2) * t50 - qJD(3) * t98 - qJD(4) * t109, qJD(2) * t51 + qJD(3) * t99 + qJD(4) * t108, -(qJD(2) * t239 - t755) * t190, qJD(2) * t76 + t776 * t785, -t121 * qJD(2) + t451 * t755, -t120 * qJD(2) + (-t779 - t782) * t451, -t570, -qJD(2) * t5 - qJD(3) * t33 - qJD(4) * t31 - qJD(5) * t41, qJD(2) * t6 + qJD(3) * t34 + qJD(4) * t32 + qJD(5) * t40; 0, 0, 0, t569, t601, t597, -t599, 0, -pkin(6) * t597 - t589, pkin(6) * t599 - t588, -t365 + (t442 * t600 + t575) * t451, -0.2e1 * t447 * t571 + t451 * t754, t446 * t599 - t602, t574 + t603, -t378, -t606 + (t446 * t533 - t591) * qJD(2) + t263 * qJD(3), t605 + (t450 * t533 + t592) * qJD(2) + t262 * qJD(3), -(t617 - t620) * t354 + t638, t607 + (-t353 * t516 + t354 * t380) * qJD(2) + t643, t260 * t596 + t516 * t599 + t624, t261 * t596 - t380 * t599 + t625, -t347, -t646 + (t353 * t436 + t380 * t390 - t447 * t548) * qJD(2) + t54 * qJD(3) + t89 * qJD(4), t645 + (-t316 * t447 - t354 * t436 + t390 * t516) * qJD(2) + t55 * qJD(3) + t88 * qJD(4), (-t619 - t623) * t239 + t594, t644 + (t214 * t234 + t239 * t549) * qJD(2) + t595, -t214 * t599 + t778 * t785 - t629, t549 * t599 + t784 * t785 - t630, -t318, -t714 + (t234 * t343 - t291 * t549 - t447 * t806) * qJD(2) + t3 * qJD(3) + t9 * qJD(4) + t18 * qJD(5), t713 + (-t141 * t447 - t214 * t291 + t239 * t343) * qJD(2) + t4 * qJD(3) + t10 * qJD(4) + t17 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t474, -t473, t446 * t511, t484, t429, qJD(2) * t263 - qJD(3) * t340 - t621, qJD(2) * t262 + qJD(3) * t339 + t622, -t105, t82, t107, t106, t429, qJD(2) * t54 + qJD(3) * t198 + qJD(4) * t93 - t705, qJD(2) * t55 - qJD(3) * t199 + qJD(4) * t95 + t704, t42, t11, t39, t38, t429, qJD(2) * t3 + qJD(3) * t77 + qJD(4) * t29 + qJD(5) * t22 - t710, qJD(2) * t4 - qJD(3) * t78 + qJD(4) * t30 + qJD(5) * t21 + t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t82, t107, t106, t429, qJD(2) * t89 + qJD(3) * t93 - qJD(4) * t182 - t631, qJD(2) * t88 + qJD(3) * t95 + qJD(4) * t181 + t632, t42, t11, t39, t38, t429, qJD(2) * t9 + qJD(3) * t29 + qJD(4) * t59 + qJD(5) * t24 - t712, qJD(2) * t10 + qJD(3) * t30 - qJD(4) * t60 + qJD(5) * t26 + t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t488, t521, -t542 * t550 + t759 + t786, -t190 * t542 + t782 + t792, t429, qJD(2) * t18 + qJD(3) * t22 + qJD(4) * t24 - qJD(5) * t58 - t707, qJD(2) * t17 + qJD(3) * t21 + qJD(4) * t26 + qJD(5) * t57 + t708; 0, 0, 0, -t569, -t601, 0, 0, 0, t589, t588, -t442 * t569 - t365, 0.2e1 * t446 * t484, -t572 + t602, t573 - t603, t378, qJD(3) * t305 + t606, qJD(3) * t304 - t605, -t354 * t620 + t638, -t607 + t643, -t259 * t596 - t624, -t258 * t596 - t625, t347, -qJD(3) * t53 - qJD(4) * t86 + t646, -qJD(3) * t52 - qJD(4) * t87 - t645, t239 * t623 + t594, t595 - t644, -t777 * t785 + t629, -t783 * t785 + t630, t318, -qJD(3) * t1 - qJD(4) * t8 - qJD(5) * t15 + t714, -qJD(3) * t2 - qJD(4) * t7 - qJD(5) * t16 - t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t571, t418 * qJD(3), 0, 0, 0, -pkin(2) * t613, -pkin(2) * t612, -t380 * t544, t596 * t240, 0, 0, 0, qJD(3) * t254 + t516 * t610, qJD(3) * t255 - t380 * t610, -t762 * t214, t785 * t775, 0, 0, 0, qJD(3) * t102 - qJD(4) * t110 - t214 * t608, qJD(3) * t103 - qJD(4) * t111 + t549 * t608; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, t754, -t543 * t450, t543 * t446, t431, -pkin(7) * t612 - t501, pkin(7) * t613 - t502, t118, t90, t179, t180, t431, -qJD(3) * t316 + qJD(4) * t247 + t514, -t513 + t753, t46, t13, t49, t48, t431, t525 - t809, t524 + t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t90, t179, t180, t431, qJD(3) * t247 - qJD(4) * t316 - t485, -t486 + t753, t46, t13, t49, t48, t431, t523 - t809, t522 + t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t487, t515, t49, t48, t431, -t490 - t809, -t489 + t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474, t473, (-t446 * t600 + t614) * t451, (-t576 - t615) * t451, t429, -qJD(2) * t305 + t621, -qJD(2) * t304 - t622, t105, -t82, t203, t202, t429, qJD(2) * t53 + qJD(4) * t92 + t705, qJD(2) * t52 + qJD(4) * t94 - t704, -t42, -t11, t62, t61, t429, qJD(2) * t1 - qJD(4) * t27 - qJD(5) * t19 + t710, qJD(2) * t2 - qJD(4) * t28 - qJD(5) * t20 - t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t483, -t754, t450 * t598, -t446 * t598, t430, t501, t502, -t118, -t90, t249, t248, t430, qJD(4) * t309 - t514, t513, -t46, -t13, t788, t793, t430, -t525, -t524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t445 * t717, -t449 * t717, 0, 0, 0, 0, 0, -t363 - t349, -t364 + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t445 * t553 + t512, -t449 * t553 + t706, 0, 0, 0, 0, 0, qJD(5) * t302 - t363 - t466, qJD(5) * t303 - t364 - t465; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t302 - t349 - t470, qJD(4) * t303 + t348 - t469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t82, t203, t202, t429, qJD(2) * t86 - qJD(3) * t92 + t631, qJD(2) * t87 - qJD(3) * t94 - t632, -t42, -t11, t62, t61, t429, qJD(2) * t8 + qJD(3) * t27 + qJD(5) * t23 + t712, qJD(2) * t7 + qJD(3) * t28 + qJD(5) * t25 - t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, -t90, t249, t248, t430, -qJD(3) * t309 + t485, t486, -t46, -t13, t788, t793, t430, -t523, -t522; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t445 * t718 - t512, t449 * t718 - t706, 0, 0, 0, 0, 0, -qJD(5) * t322 + t466, -qJD(5) * t323 + t465; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444 * t715, -t448 * t715; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444 * t552 + t468, -t448 * t552 + t467; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488, -t521, t62, t61, t429, qJD(2) * t15 + qJD(3) * t19 - qJD(4) * t23 + t707, qJD(2) * t16 + qJD(3) * t20 - qJD(4) * t25 - t708; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, -t515, t788, t793, t430, t490, t489; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t322 + t470, qJD(4) * t323 + t469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t444 * t716 - t468, t448 * t716 - t467; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t12;
