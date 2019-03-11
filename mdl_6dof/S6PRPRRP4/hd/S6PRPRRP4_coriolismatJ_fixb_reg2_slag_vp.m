% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRRP4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:26
% EndTime: 2019-03-08 20:12:46
% DurationCPUTime: 14.95s
% Computational Cost: add. (10653->575), mult. (24541->733), div. (0->0), fcn. (28179->10), ass. (0->459)
t465 = sin(pkin(11));
t720 = cos(qJ(4));
t598 = t720 * t465;
t467 = cos(pkin(11));
t469 = sin(qJ(4));
t663 = t469 * t467;
t435 = t598 + t663;
t466 = sin(pkin(6));
t471 = cos(qJ(2));
t669 = t466 * t471;
t376 = t435 * t669;
t779 = -t376 / 0.2e1;
t468 = sin(qJ(5));
t668 = t468 * qJ(6);
t470 = cos(qJ(5));
t710 = t470 * pkin(5);
t552 = t668 + t710;
t442 = -pkin(4) - t552;
t719 = sin(qJ(2));
t599 = t466 * t719;
t450 = t470 * t599;
t597 = t720 * t467;
t664 = t469 * t465;
t508 = t597 - t664;
t377 = t508 * t669;
t665 = t468 * t377;
t309 = -t450 + t665;
t572 = t468 * t599;
t659 = t470 * t377;
t310 = t572 + t659;
t722 = t468 / 0.2e1;
t769 = t310 * t470 / 0.2e1 + t309 * t722;
t776 = t769 * pkin(9);
t778 = t442 * t779 - t776;
t777 = pkin(4) * t779 + t776;
t429 = t508 ^ 2;
t430 = t435 ^ 2;
t606 = -t430 - t429;
t770 = t606 * t470;
t775 = qJD(2) * t770;
t771 = t606 * t468;
t774 = qJD(2) * t771;
t773 = qJD(3) * t770;
t772 = qJD(3) * t771;
t463 = t468 ^ 2;
t464 = t470 ^ 2;
t768 = -t463 - t464;
t708 = cos(pkin(6));
t415 = t465 * t708 + t467 * t599;
t493 = t465 * t599 - t467 * t708;
t305 = t415 * t469 + t493 * t720;
t728 = t435 / 0.2e1;
t589 = t305 * t728;
t306 = t415 * t720 - t469 * t493;
t256 = t306 * t470 - t468 * t669;
t686 = t256 * t470;
t667 = t468 * t306;
t255 = t470 * t669 + t667;
t688 = t255 * t468;
t482 = -(-t686 / 0.2e1 - t688 / 0.2e1) * t508 + t589;
t682 = t310 * t468;
t683 = t309 * t470;
t747 = -t682 / 0.2e1 + t683 / 0.2e1;
t477 = t482 - t747;
t685 = t305 * t376;
t498 = t255 * t309 + t256 * t310 + t685;
t749 = t498 * qJD(1);
t767 = qJD(3) * t477 + t749;
t475 = t482 + t747;
t766 = qJD(3) * t475 - t749;
t765 = t493 * t465;
t605 = t430 - t429;
t328 = t468 * t508;
t619 = t328 * qJD(2);
t628 = qJD(5) * t468;
t291 = -t619 + t628;
t314 = t328 * qJD(5);
t267 = t605 * t468;
t624 = t267 * qJD(2);
t764 = t624 + t314;
t610 = t435 * qJD(4);
t410 = t470 * t610;
t763 = t624 - t410;
t762 = qJD(1) * t475;
t721 = -t470 / 0.2e1;
t497 = (t255 * t721 + t256 * t722) * t508;
t478 = -t497 - t769;
t761 = qJD(1) * t478;
t733 = -t255 / 0.2e1;
t515 = (t733 + t667 / 0.2e1) * t435;
t673 = t376 * t470;
t485 = t673 / 0.2e1 + t515;
t760 = qJD(1) * t485;
t333 = t470 * t435;
t583 = -t333 / 0.2e1;
t687 = t256 * t508;
t505 = -t687 / 0.2e1 + t305 * t583;
t648 = t665 / 0.2e1 - t450 / 0.2e1;
t487 = t505 - t648;
t759 = qJD(1) * t487;
t758 = qJD(2) * t475;
t757 = qJD(2) * t477;
t756 = qJD(2) * t478;
t755 = qJD(2) * t485;
t754 = qJD(2) * t487;
t753 = qJD(2) * (t310 * t508 + t333 * t376);
t751 = qJD(4) * t478;
t499 = (t306 - t686 - t688) * t305;
t748 = t499 * qJD(1);
t612 = t508 * qJD(2);
t573 = qJD(5) - t612;
t416 = t463 * t508;
t417 = t464 * t508;
t746 = 0.2e1 * t468 * t333 * (-t612 - qJD(5)) - (t416 - t417) * qJD(4);
t481 = -t497 + t769;
t633 = qJD(4) * t305;
t745 = qJD(2) * t481 + t633 * t768;
t150 = t305 * t468;
t495 = -t673 / 0.2e1 + t515;
t631 = qJD(4) * t470;
t744 = qJD(2) * t495 + qJD(5) * t150 - t306 * t631;
t582 = t333 / 0.2e1;
t566 = t305 * t582 + t687 / 0.2e1;
t517 = t566 - t648;
t653 = qJD(4) * t150 - qJD(5) * t256;
t743 = qJD(2) * t517 + t653;
t611 = t435 * qJD(2);
t742 = qJD(4) * t481 + (-t682 + t683) * t611;
t741 = qJD(2) * t498 + qJD(4) * t499;
t740 = qJD(4) * t485 - qJD(5) * t487;
t329 = t468 * t435;
t739 = qJD(4) * t495 + qJD(5) * t517 + (t309 * t508 + t329 * t376) * qJD(2);
t738 = -qJ(6) / 0.2e1;
t459 = -pkin(3) * t467 - pkin(2);
t716 = pkin(9) * t435;
t559 = -pkin(4) * t508 - t716;
t500 = t459 + t559;
t709 = pkin(8) + qJ(3);
t447 = t709 * t467;
t579 = t709 * t465;
t368 = t447 * t720 - t469 * t579;
t660 = t470 * t368;
t188 = t468 * t500 + t660;
t671 = t508 * qJ(6);
t155 = t188 - t671;
t737 = -t155 / 0.2e1;
t666 = t468 * t368;
t187 = -t470 * t500 + t666;
t714 = t508 * pkin(5);
t156 = t187 + t714;
t736 = t156 / 0.2e1;
t712 = t435 * pkin(4);
t713 = t508 * pkin(9);
t354 = t712 - t713;
t661 = t470 * t354;
t367 = t447 * t469 + t579 * t720;
t677 = t367 * t468;
t212 = t661 + t677;
t711 = t435 * pkin(5);
t164 = -t212 - t711;
t735 = t164 / 0.2e1;
t734 = t188 / 0.2e1;
t732 = t256 / 0.2e1;
t731 = t306 / 0.2e1;
t351 = t367 * t470;
t346 = t351 / 0.2e1;
t730 = -t354 / 0.2e1;
t729 = t508 / 0.2e1;
t727 = t442 / 0.2e1;
t662 = t470 * qJ(6);
t718 = pkin(5) * t468;
t448 = t662 - t718;
t726 = -t448 / 0.2e1;
t725 = t448 / 0.2e1;
t724 = -t463 / 0.2e1;
t723 = -t464 / 0.2e1;
t580 = -t187 / 0.2e1 + t736;
t520 = -t714 / 0.2e1 - t580;
t556 = -t671 / 0.2e1 + t155 / 0.2e1;
t16 = (t734 - t556) * t470 + t520 * t468;
t707 = qJD(2) * t16;
t14 = ((t155 - t188) * t470 + (t156 - t187) * t468) * t435;
t706 = t14 * qJD(2);
t697 = t188 * t468;
t15 = -t697 / 0.2e1 + t556 * t468 + t520 * t470;
t705 = t15 * qJD(2);
t704 = t155 * t508;
t703 = t155 * t470;
t702 = t156 * t468;
t701 = t187 * t508;
t700 = t187 * t468;
t699 = t187 * t470;
t698 = t188 * t508;
t696 = t188 * t470;
t695 = t212 * t470;
t343 = t468 * t354;
t213 = -t351 + t343;
t694 = t213 * t468;
t558 = pkin(5) * t329 - qJ(6) * t333;
t214 = t367 + t558;
t693 = t214 * t468;
t692 = t214 * t470;
t215 = -t448 * t508 + t368;
t691 = t215 * t468;
t690 = t215 * t470;
t689 = t255 * t508;
t678 = t367 * t435;
t676 = t376 * t367;
t674 = t376 * t468;
t672 = t415 * t467;
t670 = t448 * t468;
t332 = t470 * t508;
t528 = (t723 + t724) * t713;
t587 = t435 * t727;
t490 = -t528 + t587;
t427 = t435 * qJ(6);
t163 = t427 + t213;
t513 = t163 * t722 + t164 * t721;
t51 = t490 - t513;
t656 = t51 * qJD(2);
t602 = -t712 / 0.2e1;
t491 = -t528 + t602;
t511 = t695 / 0.2e1 + t694 / 0.2e1;
t66 = t491 - t511;
t655 = t66 * qJD(2);
t449 = t466 ^ 2 * t719 * t471;
t94 = t306 * t377 - t449 + t685;
t654 = t94 * qJD(1);
t652 = t693 / 0.2e1 + t442 * t582;
t651 = t768 * pkin(9) * t305;
t650 = t343 / 0.2e1 - t351 / 0.2e1;
t649 = -t343 / 0.2e1 + t346;
t563 = t599 / 0.2e1;
t646 = t659 / 0.2e1 + t468 * t563;
t645 = -t659 / 0.2e1 - t572 / 0.2e1;
t644 = t435 * t724 + t464 * t728;
t451 = t465 ^ 2 + t467 ^ 2;
t269 = t605 * t470;
t641 = qJD(2) * t269;
t637 = qJD(2) * t466;
t636 = qJD(3) * t328;
t635 = qJD(3) * t468;
t634 = qJD(3) * t470;
t632 = qJD(4) * t468;
t630 = qJD(5) * t187;
t629 = qJD(5) * t508;
t627 = qJD(5) * t470;
t626 = qJD(6) * t468;
t216 = -t449 + (t672 + t765) * t669;
t625 = t216 * qJD(1);
t428 = t598 / 0.2e1 + t663 / 0.2e1;
t272 = (t728 - t428) * t669;
t623 = t272 * qJD(1);
t503 = -t597 / 0.2e1 + t664 / 0.2e1;
t273 = (t729 + t503) * t669;
t622 = t273 * qJD(1);
t621 = t605 * qJD(2);
t326 = (t463 / 0.2e1 + t723) * t435;
t620 = t326 * qJD(5);
t618 = t329 * qJD(2);
t319 = t332 * qJD(2);
t617 = t333 * qJD(2);
t339 = -t416 - t417;
t615 = t339 * qJD(2);
t614 = t606 * qJD(2);
t613 = t428 * qJD(2);
t425 = t508 * qJD(4);
t423 = t508 * qJD(6);
t609 = t451 * qJD(2);
t455 = t464 - t463;
t608 = t455 * qJD(5);
t607 = t470 * qJD(6);
t604 = pkin(9) * t628;
t603 = pkin(9) * t627;
t419 = t711 / 0.2e1;
t584 = -t332 / 0.2e1;
t601 = -t693 / 0.2e1 + t442 * t583 + pkin(9) * t584;
t337 = t661 / 0.2e1;
t344 = t677 / 0.2e1;
t600 = t337 + t344 + t419;
t596 = t508 * t628;
t595 = t435 * t627;
t594 = t468 * t607;
t593 = t508 * t611;
t363 = t508 * t610;
t457 = t468 * t627;
t592 = t468 * t611;
t591 = t435 * t626;
t456 = t468 * t631;
t590 = t470 * t611;
t586 = -t329 / 0.2e1;
t585 = t332 / 0.2e1;
t581 = t737 + t734;
t578 = t451 * t471;
t262 = t428 + t644;
t577 = qJD(2) * t262 + t456;
t279 = qJD(2) * t326 - t456;
t383 = t470 * t430 * t468 * qJD(2);
t221 = qJD(4) * t326 + t383;
t350 = t508 * t590;
t576 = -qJD(4) * t328 - t350;
t575 = qJD(2) * t459 + qJD(3);
t571 = t468 * t590;
t570 = t468 * t410;
t569 = t430 * t457;
t230 = t305 * t585;
t568 = t468 * t589 + t689 / 0.2e1;
t567 = t305 * t586 - t689 / 0.2e1;
t565 = -t711 / 0.2e1 - t677 / 0.2e1;
t564 = qJD(2) * t599;
t321 = t552 * t435;
t562 = -t321 / 0.2e1 + t713 / 0.2e1;
t557 = 0.2e1 * t570;
t554 = -t350 + t595;
t551 = -t442 * t508 + t716;
t11 = t155 * t163 + t156 * t164 + t214 * t215;
t514 = -t703 / 0.2e1 - t702 / 0.2e1;
t473 = (t215 / 0.2e1 + t514) * t305 + t255 * t735 + t163 * t732 + t214 * t731;
t2 = t473 + t778;
t550 = t2 * qJD(1) + t11 * qJD(2);
t12 = -t155 * t187 + t156 * t188 + t214 * t321;
t480 = t581 * t255 + t580 * t256 + t305 * t321 / 0.2e1;
t522 = t309 * pkin(5) / 0.2e1 + t310 * t738;
t6 = t480 + t522;
t549 = qJD(1) * t6 + qJD(2) * t12;
t22 = -t187 * t212 + t188 * t213 + t367 * t368;
t512 = -t700 / 0.2e1 - t696 / 0.2e1;
t472 = (t368 / 0.2e1 + t512) * t305 + t212 * t733 + t213 * t732 + t367 * t731;
t4 = t472 - t777;
t548 = t4 * qJD(1) + t22 * qJD(2);
t13 = -t156 * t332 - t164 * t333 + (t163 * t435 + t704) * t468;
t547 = -qJD(2) * t13 + t761;
t19 = (t155 - t690) * t435 - (t163 + t692) * t508;
t359 = t674 / 0.2e1;
t504 = t256 * t728 + t230;
t488 = t305 * t584 + t306 * t583 + t504;
t37 = t359 + t488;
t546 = t37 * qJD(1) + t19 * qJD(2);
t20 = (-t156 + t691) * t435 - (-t164 - t693) * t508;
t545 = t20 * qJD(2) + t760;
t21 = (t694 + t695) * t435 - (-t697 + t699) * t508;
t544 = -qJD(2) * t21 + t761;
t24 = t214 * t435 - (-t702 - t703) * t508;
t543 = qJD(2) * t24 + t762;
t61 = t678 - (-t696 - t700) * t508;
t542 = -qJD(2) * t61 - t762;
t47 = (-t187 + t666) * t435 - (t212 - t677) * t508;
t541 = t47 * qJD(2) + t760;
t360 = -t674 / 0.2e1;
t489 = t306 * t582 + t230 - t504;
t40 = t360 + t489;
t48 = (-t188 + t660) * t435 - (-t213 - t351) * t508;
t540 = t40 * qJD(1) + t48 * qJD(2);
t59 = t698 + (t321 * t468 + t692) * t435;
t539 = qJD(2) * t59 - t759;
t60 = t701 + (-t321 * t470 + t693) * t435;
t74 = t567 + t646;
t538 = qJD(1) * t74 - qJD(2) * t60;
t68 = -t214 * t333 - t704;
t71 = t566 + t648;
t537 = -qJD(1) * t71 + qJD(2) * t68;
t536 = t163 * t470 + t164 * t468;
t535 = -t212 * t468 + t213 * t470;
t103 = -t329 * t367 - t701;
t75 = t568 + t645;
t533 = -qJD(1) * t75 + qJD(2) * t103;
t104 = t333 * t367 + t698;
t532 = qJD(2) * t104 - t759;
t31 = (-t670 / 0.2e1 - pkin(5) / 0.2e1) * t435 + (t730 + t562) * t470 + t565 + t652;
t373 = -t442 * t468 - t448 * t470;
t531 = -qJD(2) * t31 + qJD(4) * t373;
t372 = t442 * t470 - t670;
t476 = (t435 * t725 - t214 / 0.2e1) * t470 + (t587 + t562) * t468;
t44 = -t427 + t476 + t649;
t530 = -qJD(2) * t44 + qJD(4) * t372;
t529 = t573 * t468;
t510 = t306 * t729 + t589;
t113 = t563 - t510;
t159 = t368 * t508 + t678;
t527 = qJD(1) * t113 - qJD(2) * t159;
t479 = t672 / 0.2e1 + t765 / 0.2e1;
t260 = t563 - t479;
t441 = t451 * qJ(3);
t526 = qJD(1) * t260 - qJD(2) * t441;
t152 = t305 * t470;
t525 = -qJD(4) * t152 - qJD(5) * t255;
t524 = -qJD(5) * t448 - t626;
t523 = pkin(5) * t735 + t163 * t738;
t521 = -t662 / 0.2e1 + t718 / 0.2e1;
t394 = pkin(9) * t585;
t99 = t394 + (t602 + t730) * t470;
t519 = pkin(4) * t632 - qJD(2) * t99;
t483 = (-t713 / 0.2e1 + t712 / 0.2e1) * t468 + t346;
t97 = t483 + t650;
t518 = pkin(4) * t631 - qJD(2) * t97;
t65 = t600 + t601;
t507 = -qJD(2) * t65 + t442 * t632;
t153 = t435 * t529;
t276 = qJD(5) * t428 - t593;
t506 = -qJD(5) * t152 - t306 * t632;
t418 = t464 * t430;
t340 = t430 * t463 - t418;
t245 = -qJD(2) * t340 + t557;
t352 = -qJD(4) * t455 + 0.2e1 * t571;
t502 = qJD(4) * t267 - t508 * t595;
t494 = qJD(5) * t340 - t508 * t557;
t45 = (t725 + t521) * t305;
t474 = (t468 * t581 + t470 * t580) * pkin(9) + t214 * t726 + t321 * t727;
t8 = t474 + t523;
t492 = t442 * t448 * qJD(4) + t45 * qJD(1) - t8 * qJD(2);
t486 = -qJD(5) * t552 + t607;
t349 = t418 + t429;
t484 = qJD(2) * t349 + t570 - t629;
t412 = t428 * qJD(4);
t409 = t468 * t610;
t384 = t470 * t591;
t382 = t573 * qJ(6);
t374 = qJD(4) * t463 + t571;
t324 = t339 * qJD(3);
t322 = t339 * qJD(4);
t318 = t332 * qJD(5);
t293 = -t319 + t627;
t275 = (-t435 / 0.2e1 - t428) * t669;
t274 = (-t508 / 0.2e1 + t503) * t669;
t263 = -t428 + t644;
t261 = t563 + t479;
t258 = t363 * t464 - t569;
t257 = t363 * t463 + t569;
t254 = t470 * t153;
t199 = -t464 * t593 - t620;
t198 = qJD(4) * t332 - t508 * t592;
t197 = -t463 * t593 + t620;
t176 = t188 * qJD(5);
t162 = qJD(4) * t269 + t435 * t596;
t154 = -t318 - t641;
t148 = -t620 - (-t464 * t611 - t456) * t508;
t147 = t620 - (-t463 * t611 + t456) * t508;
t129 = t409 + t641;
t114 = t563 + t510;
t100 = pkin(4) * t583 + t337 + 0.2e1 * t344 + t394;
t98 = t483 + t649;
t78 = t505 + t648;
t77 = t567 + t645;
t76 = t568 + t646;
t67 = t491 + t511;
t64 = -t661 / 0.2e1 + t565 + t601;
t52 = t490 + t513;
t46 = (t521 + t726) * t305;
t43 = t427 + t476 + t650;
t39 = t359 + t489;
t38 = t360 + t488;
t32 = t448 * t586 + t470 * t562 + t419 + t600 + t652;
t18 = -t699 / 0.2e1 + t468 * t737 + t697 / 0.2e1 + t470 * t736 - (t668 / 0.2e1 + t710 / 0.2e1) * t508;
t17 = -t508 * t521 + t512 - t514;
t7 = t474 - t523;
t5 = t480 - t522;
t3 = t472 + t777;
t1 = t473 - t778;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t94, 0, 0, 0, 0, 0, 0, 0, 0, 0, t741, 0, 0, 0, 0, 0, 0, 0, 0, 0, t741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t564, -t471 * t637, 0, 0, 0, 0, 0, 0, 0, 0, -t467 * t564, t465 * t564, t578 * t637, t625 + t261 * qJD(3) + (-pkin(2) * t719 + qJ(3) * t578) * t637, 0, 0, 0, 0, 0, 0, qJD(4) * t275 - t508 * t564, qJD(4) * t274 + t435 * t564 (t376 * t435 + t377 * t508) * qJD(2), t654 + (t368 * t377 + t459 * t599 + t676) * qJD(2) + t114 * qJD(3), 0, 0, 0, 0, 0, 0, t739, qJD(4) * t39 + qJD(5) * t77 + t753, t742 (t187 * t309 + t188 * t310 + t676) * qJD(2) + t3 * qJD(4) + t767, 0, 0, 0, 0, 0, 0, t739, t742, qJD(4) * t38 + qJD(5) * t76 - t753 (t155 * t310 + t156 * t309 + t214 * t376) * qJD(2) + t1 * qJD(4) + t5 * qJD(5) + t78 * qJD(6) + t767; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t757, 0, 0, 0, 0, 0, 0, 0, 0, 0, t757; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t275 - qJD(4) * t306, qJD(2) * t274 + t633, 0, 0, 0, 0, 0, 0, 0, 0, t744, qJD(2) * t39 - t506, t745, t748 + t3 * qJD(2) + (-pkin(4) * t306 + t651) * qJD(4), 0, 0, 0, 0, 0, 0, t744, t745, qJD(2) * t38 + t506, t748 + t1 * qJD(2) + (t306 * t442 + t651) * qJD(4) + t46 * qJD(5) - t150 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t743, qJD(2) * t77 - t525, 0, 0, 0, 0, 0, 0, 0, 0, t743, 0, qJD(2) * t76 + t525, t5 * qJD(2) + t46 * qJD(4) + (-pkin(5) * t256 - qJ(6) * t255) * qJD(5) + t256 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t78 - t653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t260 - t625, 0, 0, 0, 0, 0, 0, -t272 * qJD(4), -t273 * qJD(4), 0, -qJD(3) * t113 - t654, 0, 0, 0, 0, 0, 0, t740, qJD(4) * t40 - qJD(5) * t75, t751, qJD(4) * t4 + t766, 0, 0, 0, 0, 0, 0, t740, t751, qJD(4) * t37 - qJD(5) * t74, qJD(4) * t2 + qJD(5) * t6 - qJD(6) * t71 + t766; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451 * qJD(3), t441 * qJD(3), t363, -t605 * qJD(4), 0, -t363, 0, 0, t459 * t610, t459 * t425, -qJD(3) * t606, qJD(3) * t159, t258, t494, t162, t257, -t502, -t363, qJD(4) * t47 + qJD(5) * t104 - t772, qJD(4) * t48 + qJD(5) * t103 - t773, -qJD(4) * t21, qJD(3) * t61 + qJD(4) * t22, t258, t162, -t494, -t363, t502, t257, qJD(4) * t20 + qJD(5) * t59 - t430 * t594 - t772, -qJD(4) * t13 - qJD(5) * t14 + t508 * t591, qJD(4) * t19 + qJD(5) * t60 + qJD(6) * t349 + t773, qJD(3) * t24 + qJD(4) * t11 + qJD(5) * t12 + qJD(6) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t609, -t526, 0, 0, 0, 0, 0, 0, 0, 0, -t614, -t527, 0, 0, 0, 0, 0, 0, -t774, -t775, 0, qJD(4) * t67 - t542, 0, 0, 0, 0, 0, 0, -t774, 0, t775, qJD(4) * t52 + qJD(5) * t17 + t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t593, -t621, t425, -t593, -t610, 0, -qJD(4) * t368 + t459 * t611 - t623, qJD(4) * t367 + t459 * t612 - t622, 0, 0, t148, t746, t129, t147, -t763, t276 (t468 * t559 - t660) * qJD(4) + t100 * qJD(5) + t541 (t470 * t559 + t666) * qJD(4) + t98 * qJD(5) + t540, qJD(4) * t535 + t544, t67 * qJD(3) + (-t368 * pkin(4) + pkin(9) * t535) * qJD(4) + t548, t148, t129, -t746, t276, t763, t147 (-t468 * t551 - t690) * qJD(4) + t32 * qJD(5) + t263 * qJD(6) + t545, qJD(4) * t536 + qJD(5) * t18 + t547 (t470 * t551 - t691) * qJD(4) + t43 * qJD(5) + t384 + t546, t52 * qJD(3) + (pkin(9) * t536 + t215 * t442) * qJD(4) + t7 * qJD(5) + t64 * qJD(6) + t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, -t245, -t153, t221, -t554, t412, qJD(4) * t100 - t176 + t532, qJD(4) * t98 + t533 + t630, 0, 0, -t221, -t153, t245, t412, t554, t221, qJD(4) * t32 - t176 + t539, t18 * qJD(4) + qJD(5) * t558 - t591 - t706, qJD(4) * t43 - t423 - t538 - t630, t17 * qJD(3) + t7 * qJD(4) + (-pkin(5) * t188 - qJ(6) * t187) * qJD(5) + t155 * qJD(6) + t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t263 - t383, -t153, t484, qJD(4) * t64 + qJD(5) * t155 + t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t758, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t609, t526, 0, 0, 0, 0, 0, 0, t610, t425, t614, t527, 0, 0, 0, 0, 0, 0, t314 + t410 + t774, -qJD(4) * t329 + t508 * t627 + t775, t322, -qJD(4) * t66 + t542, 0, 0, 0, 0, 0, 0, qJD(4) * t333 + t596 + t774, t322, -t318 + t409 - t775, -qJD(4) * t51 - qJD(5) * t16 - qJD(6) * t328 - t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t611, t612, 0, 0, 0, 0, 0, 0, 0, 0, t590, -t618, t615, -t655, 0, 0, 0, 0, 0, 0, t617, t615, t592, -t656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, -t573 * t470, 0, 0, 0, 0, 0, 0, 0, 0, -t529, 0, t293, -t524 - t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272 * qJD(2), t273 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t755, -qJD(2) * t40, -t756, -qJD(2) * t4 - t748, 0, 0, 0, 0, 0, 0, -t755, -t756, -qJD(2) * t37, -qJD(2) * t2 - qJD(5) * t45 - t748; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t593, t621, 0, t593, 0, 0, -t435 * t575 + t623, -t508 * t575 + t622, 0, 0, t199, -0.2e1 * t254, t154, t197, t764, -t276, qJD(5) * t99 - t435 * t634 - t541, qJD(3) * t329 + qJD(5) * t97 - t540, -t324 - t544, qJD(3) * t66 - t548, t199, t154, 0.2e1 * t254, -t276, -t764, t197, -qJD(3) * t333 + qJD(5) * t31 + qJD(6) * t262 - t545, -qJD(5) * t15 - qJD(6) * t332 - t324 - t547, qJD(5) * t44 - t435 * t635 + t384 - t546, qJD(3) * t51 + qJD(5) * t8 + qJD(6) * t65 - t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t611, -t612, 0, 0, 0, 0, 0, 0, 0, 0, -t590, t618, -t615, t655, 0, 0, 0, 0, 0, 0, -t617, -t615, -t592, t656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t457, t608, 0, -t457, 0, 0, -pkin(4) * t628, -pkin(4) * t627, 0, 0, t457, 0, -t608, 0, 0, -t457, -qJD(5) * t373 + t594, 0, -qJD(5) * t372 + qJD(6) * t463, t524 * t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, -t352, t293, t279, -t291, -t613, -t519 - t603, -t518 + t604, 0, 0, -t279, t293, t352, -t613, t291, t279, -t531 - t603, t486 - t705, -t530 - t604, pkin(9) * t486 - t492; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, t293, t374, -t507 + t603; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t754, qJD(2) * t75, 0, 0, 0, 0, 0, 0, 0, 0, t754, 0, qJD(2) * t74, -qJD(2) * t6 + qJD(4) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t245, t198, -t221, t576, t412, -qJD(4) * t99 - t532 - t636, -qJD(4) * t97 - t508 * t634 - t533, 0, 0, t221, t198, -t245, t412, -t576, -t221, -qJD(4) * t31 - t508 * t635 - t539, qJD(4) * t15 + t706, qJD(3) * t332 - qJD(4) * t44 - t423 + t538, -qJ(6) * t423 + qJD(3) * t16 - qJD(4) * t8 - t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t619, -t470 * t612, 0, 0, 0, 0, 0, 0, 0, 0, -t468 * t612, 0, t319, t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, t352, t319, -t279, -t619, t613, t519, t518, 0, 0, t279, t319, -t352, t613, t619, -t279, t531, t705, t530, t492; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t573, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t262 + t383, t198, -t484, qJ(6) * t629 - qJD(4) * t65 - t537 + t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t577, t319, -t374, t507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t573, -t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t9;
