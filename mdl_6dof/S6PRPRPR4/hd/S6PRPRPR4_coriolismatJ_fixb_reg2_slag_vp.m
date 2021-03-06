% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRPR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:24
% EndTime: 2019-03-08 19:41:45
% DurationCPUTime: 14.04s
% Computational Cost: add. (16425->592), mult. (37085->839), div. (0->0), fcn. (44806->12), ass. (0->447)
t480 = sin(pkin(11));
t737 = cos(qJ(4));
t619 = t737 * t480;
t483 = cos(pkin(11));
t485 = sin(qJ(4));
t681 = t485 * t483;
t451 = t619 + t681;
t481 = sin(pkin(6));
t486 = cos(qJ(2));
t690 = t481 * t486;
t394 = t451 * t690;
t586 = -t394 / 0.2e1;
t600 = t690 / 0.2e1;
t675 = t451 * t600;
t731 = pkin(8) + qJ(3);
t597 = t731 * t480;
t584 = t485 * t597;
t462 = t731 * t483;
t620 = t737 * t462;
t392 = t620 - t584;
t482 = cos(pkin(12));
t370 = t482 * t392;
t479 = sin(pkin(12));
t474 = -t483 * pkin(3) - pkin(2);
t618 = t737 * t483;
t682 = t485 * t480;
t527 = t618 - t682;
t733 = t527 * pkin(4);
t553 = t474 - t733;
t726 = qJ(5) * t451;
t510 = t553 - t726;
t232 = t479 * t510 + t370;
t362 = t479 * t451;
t200 = -pkin(9) * t362 + t232;
t484 = sin(qJ(6));
t691 = t479 * t392;
t730 = pkin(9) + qJ(5);
t489 = -t527 * pkin(5) - t691 + (-t451 * t730 + t553) * t482;
t736 = cos(qJ(6));
t107 = t200 * t484 - t489 * t736;
t108 = t200 * t736 + t484 * t489;
t613 = t736 * t479;
t683 = t484 * t482;
t449 = t613 + t683;
t280 = t451 * t449;
t461 = t730 * t482;
t596 = t730 * t479;
t389 = t461 * t484 + t596 * t736;
t741 = t449 / 0.2e1;
t612 = t736 * t482;
t684 = t484 * t479;
t445 = -t612 + t684;
t746 = -t445 / 0.2e1;
t391 = t461 * t736 - t484 * t596;
t748 = t391 / 0.2e1;
t693 = t451 * t445;
t756 = -t693 / 0.2e1;
t777 = -t107 * t741 - t108 * t746 + t280 * t748 - t389 * t756;
t395 = t527 * t690;
t735 = sin(qJ(2));
t622 = t481 * t735;
t358 = -t479 * t395 + t482 * t622;
t359 = t482 * t395 + t479 * t622;
t762 = t359 * t482 / 0.2e1 - t358 * t479 / 0.2e1;
t776 = pkin(4) * t586 + t762 * qJ(5);
t419 = t451 * t684;
t354 = t451 * t612 - t419;
t181 = -t280 * t741 + t354 * t746;
t775 = t181 * qJD(6);
t774 = t280 * qJD(2);
t656 = qJD(4) * t449;
t773 = -qJD(2) * t181 + t445 * t656;
t664 = qJD(2) * t354;
t772 = qJD(4) * t181 - t280 * t664;
t771 = t280 ^ 2;
t728 = cos(pkin(6));
t508 = t480 * t622 - t483 * t728;
t500 = t485 * t508;
t431 = t480 * t728 + t483 * t622;
t621 = t737 * t431;
t347 = t621 - t500;
t689 = t482 * t347;
t308 = -t479 * t690 + t689;
t545 = -t479 * t347 - t482 * t690;
t173 = t308 * t736 + t484 * t545;
t770 = t173 / 0.2e1;
t769 = t280 / 0.2e1;
t768 = t508 * t480;
t629 = t527 * qJD(2);
t767 = t280 * t629;
t441 = t527 ^ 2;
t442 = t451 ^ 2;
t766 = -t442 - t441;
t624 = t442 - t441;
t517 = t683 / 0.2e1 + t613 / 0.2e1;
t765 = t517 + t741;
t657 = qJD(4) * t445;
t764 = t657 + t774;
t763 = qJD(6) * t280;
t439 = t619 / 0.2e1 + t681 / 0.2e1;
t172 = t308 * t484 - t545 * t736;
t676 = t172 * t741 + t173 * t746;
t440 = t445 ^ 2;
t761 = t449 ^ 2;
t346 = t431 * t485 + t508 * t737;
t198 = t346 * t449;
t760 = -t198 / 0.2e1;
t759 = -t308 / 0.2e1;
t758 = t308 / 0.2e1;
t757 = t347 / 0.2e1;
t349 = t449 * t527;
t755 = t349 / 0.2e1;
t753 = -t280 / 0.2e1;
t418 = t527 * t684;
t353 = t527 * t612 - t418;
t751 = t353 / 0.2e1;
t750 = t354 / 0.2e1;
t749 = t389 / 0.2e1;
t747 = -t418 / 0.2e1;
t745 = t445 / 0.2e1;
t744 = t527 / 0.2e1;
t743 = -t527 / 0.2e1;
t742 = -t449 / 0.2e1;
t740 = -t451 / 0.2e1;
t473 = -pkin(5) * t482 - pkin(4);
t739 = t473 / 0.2e1;
t92 = (t750 + t756) * t449 + 0.2e1 * t769 * t445;
t738 = t92 * qJD(5);
t732 = t451 * pkin(4);
t199 = t346 * t445;
t96 = t198 * t746 + t199 * t741;
t729 = t96 * qJD(4);
t727 = qJ(5) * t527;
t694 = t449 * t349;
t707 = t353 * t445;
t94 = t694 - t707;
t725 = qJD(2) * t94;
t231 = t482 * t510 - t691;
t718 = t231 * t479;
t717 = t232 * t482;
t371 = -t727 + t732;
t390 = t462 * t485 + t737 * t597;
t698 = t390 * t479;
t260 = t482 * t371 + t698;
t716 = t260 * t482;
t369 = t482 * t390;
t261 = t479 * t371 - t369;
t715 = t261 * t479;
t714 = t308 * t482;
t713 = t346 * t394;
t712 = t346 * t451;
t711 = t349 * t354;
t710 = t349 * t527;
t709 = t280 * t353;
t708 = t280 * t451;
t706 = t353 * t527;
t705 = t354 * t451;
t703 = t358 * t482;
t702 = t359 * t479;
t615 = t736 * t358;
t685 = t484 * t359;
t211 = t615 - t685;
t614 = t736 * t359;
t686 = t484 * t358;
t212 = t614 + t686;
t38 = -t172 * t211 + t173 * t212 + t713;
t700 = t38 * qJD(1);
t699 = t390 * t451;
t697 = t394 * t390;
t696 = t394 * t479;
t695 = t431 * t483;
t692 = t451 * t482;
t360 = t479 * t527;
t193 = -pkin(9) * t482 * t527 + pkin(5) * t451 + t260;
t688 = t484 * t193;
t208 = -pkin(9) * t360 + t261;
t687 = t484 * t208;
t514 = t545 * t479;
t71 = (t347 + t514 - t714) * t346;
t680 = t71 * qJD(1);
t79 = t308 * t359 + t358 * t545 + t713;
t679 = t79 * qJD(1);
t475 = t479 ^ 2;
t477 = t482 ^ 2;
t464 = t475 + t477;
t465 = t480 ^ 2 + t483 ^ 2;
t516 = -t618 / 0.2e1 + t682 / 0.2e1;
t323 = (t744 + t516) * t690;
t672 = qJD(1) * t323;
t126 = -t709 + t711;
t671 = qJD(2) * t126;
t530 = -t280 * t742 + t693 * t745;
t163 = t530 + t439;
t670 = qJD(2) * t163;
t186 = -t708 + t710;
t669 = qJD(2) * t186;
t187 = t708 + t710;
t668 = qJD(2) * t187;
t188 = t705 - t706;
t667 = qJD(2) * t188;
t189 = t705 + t706;
t666 = qJD(2) * t189;
t365 = t766 * t482;
t663 = qJD(2) * t365;
t662 = qJD(2) * t481;
t507 = t742 + t517;
t254 = t507 * t527;
t659 = qJD(4) * t254;
t255 = t765 * t527;
t658 = qJD(4) * t255;
t655 = qJD(4) * t473;
t654 = qJD(4) * t479;
t653 = qJD(4) * t482;
t652 = qJD(5) * t527;
t651 = qJD(6) * t527;
t650 = qJD(6) * t449;
t598 = -t475 / 0.2e1 - t477 / 0.2e1;
t506 = -t598 * t727 - t732 / 0.2e1;
t533 = t716 / 0.2e1 + t715 / 0.2e1;
t109 = t506 - t533;
t649 = t109 * qJD(2);
t127 = -t709 - t711;
t648 = t127 * qJD(2);
t463 = t481 ^ 2 * t735 * t486;
t148 = t347 * t395 - t463 + t713;
t647 = t148 * qJD(1);
t646 = t254 * qJD(2);
t235 = t255 * qJD(2);
t581 = -t612 / 0.2e1;
t256 = t747 - (t581 + t745) * t527;
t645 = t256 * qJD(2);
t257 = t747 - (t581 + t746) * t527;
t644 = t257 * qJD(2);
t522 = -t527 * t581 + t747;
t601 = t445 * t744;
t258 = t601 - t522;
t643 = t258 * qJD(2);
t262 = -t463 + (t695 + t768) * t690;
t642 = t262 * qJD(1);
t263 = t464 * t442;
t641 = t263 * qJD(2);
t640 = t693 * qJD(2);
t579 = t598 * t451;
t310 = t579 - t439;
t638 = t310 * qJD(2);
t319 = t624 * t479;
t637 = t319 * qJD(2);
t320 = t766 * t479;
t636 = t320 * qJD(2);
t321 = t624 * t482;
t635 = t321 * qJD(2);
t634 = t624 * qJD(2);
t633 = t362 * qJD(2);
t428 = t475 * t527;
t429 = t477 * t527;
t364 = -t428 - t429;
t632 = t364 * qJD(2);
t631 = t766 * qJD(2);
t630 = t439 * qJD(2);
t435 = t445 * qJD(6);
t436 = t527 * qJD(4);
t628 = t451 * qJD(2);
t627 = t451 * qJD(4);
t626 = t464 * qJD(4);
t625 = t465 * qJD(2);
t623 = pkin(5) * t360;
t617 = t736 * t193;
t616 = t736 * t208;
t611 = t693 * t629;
t607 = t479 * t653;
t606 = t451 * t652;
t388 = t527 * t628;
t387 = t527 * t627;
t605 = t449 * t435;
t604 = t482 * t628;
t427 = t482 * t627;
t283 = t712 / 0.2e1;
t603 = t394 * t745;
t602 = t394 * t741;
t595 = t464 * t346;
t594 = t465 * t486;
t593 = qJD(6) * t439 - t388;
t592 = qJD(2) * t474 + qJD(3);
t591 = -qJD(6) + t629;
t590 = qJD(5) + t655;
t588 = t527 * t604;
t587 = t479 * t388;
t585 = qJD(2) * t622;
t583 = t622 / 0.2e1;
t580 = t612 / 0.2e1;
t578 = -t726 - t733;
t111 = t617 - t687;
t112 = t616 + t688;
t326 = pkin(5) * t362 + t390;
t327 = t392 + t623;
t11 = -t107 * t111 + t108 * t112 + t326 * t327;
t488 = -t172 * t111 / 0.2e1 + t112 * t770 + t107 * t760 + t199 * t108 / 0.2e1 + t346 * t327 / 0.2e1 + t326 * t757;
t505 = t211 * t749 - t212 * t391 / 0.2e1 + t473 * t586;
t2 = t488 + t505;
t577 = t2 * qJD(1) + t11 * qJD(2);
t497 = t172 * t751 + t199 * t753 - t349 * t770 + t354 * t760;
t535 = t211 * t742 + t212 * t746;
t16 = t497 - t535;
t8 = t107 * t353 - t108 * t349 - t111 * t354 - t112 * t280;
t576 = t16 * qJD(1) + t8 * qJD(2);
t574 = t482 * t587;
t534 = t718 / 0.2e1 - t717 / 0.2e1;
t487 = (t392 / 0.2e1 + t534) * t346 + t261 * t758 + t390 * t757 + t545 * t260 / 0.2e1;
t18 = t487 - t776;
t51 = t231 * t260 + t232 * t261 + t390 * t392;
t572 = t18 * qJD(1) + t51 * qJD(2);
t23 = -t107 * t451 - t111 * t527 + t280 * t327 + t326 * t349;
t496 = t172 * t740 + t198 * t743 + t346 * t755 + t347 * t769;
t40 = t603 - t496;
t571 = -t40 * qJD(1) + t23 * qJD(2);
t24 = -t108 * t451 + t112 * t527 + t326 * t353 + t327 * t354;
t495 = t173 * t740 + t199 * t744 + t346 * t751 + t347 * t750;
t41 = t602 - t495;
t570 = -t41 * qJD(1) + t24 * qJD(2);
t502 = t349 * t749 + t353 * t748 + t451 * t739;
t542 = t111 * t746 + t112 * t741;
t27 = t502 - t542;
t569 = -qJD(1) * t96 + qJD(2) * t27;
t536 = t211 * t746 + t212 * t741;
t540 = t172 * t755 + t173 * t751;
t29 = -t712 / 0.2e1 + t536 - t540;
t31 = t107 * t349 + t108 * t353 + t326 * t451;
t568 = -t29 * qJD(1) + t31 * qJD(2);
t37 = -t107 * t693 - t108 * t280;
t541 = t172 * t693 / 0.2e1 + t173 * t769;
t49 = t541 + t675;
t567 = qJD(1) * t49 - qJD(2) * t37;
t559 = -t231 * t482 - t232 * t479;
t44 = (t715 + t716) * t451 - t559 * t527;
t513 = t545 * t482;
t499 = t479 * t758 + t513 / 0.2e1;
t494 = t499 * t527;
t85 = -t494 - t762;
t566 = t85 * qJD(1) - t44 * qJD(2);
t60 = -t107 * t527 - t280 * t326;
t519 = -t686 / 0.2e1 - t614 / 0.2e1;
t539 = t172 * t743 + t346 * t753;
t65 = t519 - t539;
t565 = qJD(1) * t65 - qJD(2) * t60;
t61 = t108 * t527 + t326 * t354;
t518 = -t685 / 0.2e1 + t615 / 0.2e1;
t537 = t173 * t744 + t346 * t750;
t64 = t518 - t537;
t564 = qJD(1) * t64 - qJD(2) * t61;
t498 = t714 / 0.2e1 - t514 / 0.2e1;
t491 = t498 * t527 + t283;
t529 = -t703 / 0.2e1 - t702 / 0.2e1;
t69 = t491 + t529;
t91 = t699 - (-t717 + t718) * t527;
t563 = -qJD(1) * t69 - qJD(2) * t91;
t73 = (t586 + t394 / 0.2e1) * t482;
t77 = (t231 + t691) * t451 - (t260 - t698) * t527;
t562 = t73 * qJD(1) + t77 * qJD(2);
t544 = (t759 + t689 / 0.2e1) * t451;
t76 = -t696 / 0.2e1 + t544;
t78 = (-t232 + t370) * t451 - (-t261 - t369) * t527;
t561 = t76 * qJD(1) + t78 * qJD(2);
t36 = -t172 * t198 + t173 * t199 + t346 * t347;
t560 = t36 * qJD(1) + t96 * qJD(3);
t558 = -t479 * t260 + t261 * t482;
t140 = t507 * t346;
t503 = t326 * t741 + t354 * t739 + t391 * t744;
t520 = -t687 / 0.2e1 + t617 / 0.2e1;
t45 = -t503 + t520;
t557 = qJD(1) * t140 + qJD(2) * t45;
t515 = t580 - t684 / 0.2e1;
t141 = (t745 + t515) * t346;
t504 = -t280 * t739 + t326 * t746 + t389 * t743;
t521 = -t688 / 0.2e1 - t616 / 0.2e1;
t46 = -t504 + t521;
t556 = qJD(1) * t141 + qJD(2) * t46;
t125 = -t354 * t693 + t771;
t555 = qJD(2) * t125 + qJD(4) * t92;
t374 = t440 + t761;
t554 = qJD(2) * t92 + qJD(4) * t374;
t115 = t559 * t451;
t118 = t451 * t499 + t675;
t552 = qJD(1) * t118 - qJD(2) * t115;
t532 = t347 * t744 + t283;
t168 = t583 - t532;
t210 = t392 * t527 + t699;
t551 = qJD(1) * t168 - qJD(2) * t210;
t492 = t695 / 0.2e1 + t768 / 0.2e1;
t311 = t583 - t492;
t460 = t465 * qJ(3);
t550 = qJD(1) * t311 - qJD(2) * t460;
t124 = t280 * t445 - t449 * t354;
t167 = -t354 ^ 2 + t771;
t549 = qJD(2) * t167 + qJD(4) * t124;
t344 = t440 - t761;
t548 = qJD(2) * t124 + qJD(4) * t344;
t219 = t419 / 0.2e1 + (t684 / 0.2e1 - t612) * t451;
t547 = qJD(2) * t219 - t656;
t209 = t449 * t389 - t391 * t445;
t501 = t620 / 0.2e1 - t584 / 0.2e1;
t493 = t623 / 0.2e1 + t501;
t21 = t493 + t777;
t490 = -t500 / 0.2e1 + t621 / 0.2e1;
t58 = t490 - t676;
t511 = qJD(1) * t58 + qJD(2) * t21 - qJD(4) * t209;
t113 = t501 + t534;
t128 = t490 - t498;
t459 = t464 * qJ(5);
t509 = qJD(1) * t128 + qJD(2) * t113 - qJD(4) * t459;
t430 = t439 * qJD(4);
t325 = t586 - t675;
t324 = t516 * t690 - t527 * t600;
t312 = t583 + t492;
t309 = t579 + t439;
t259 = t601 + t522;
t237 = t254 * qJD(6);
t236 = t255 * qJD(6);
t222 = -t419 / 0.2e1 + (t580 - t515) * t451;
t206 = t235 - t650;
t169 = t583 + t532;
t164 = -t530 + t439;
t143 = t765 * t346;
t142 = (t515 + t746) * t346;
t129 = t490 + t498;
t122 = t124 * qJD(6);
t119 = t362 * t759 + t513 * t740 + t675;
t114 = t501 - t534;
t110 = t506 + t533;
t84 = -t494 + t762;
t75 = t696 / 0.2e1 + t544;
t74 = 0.2e1 * t482 * t586;
t68 = t491 - t529;
t67 = t518 + t537;
t66 = t519 + t539;
t59 = t490 + t676;
t50 = -t541 + t675;
t48 = t503 + t520;
t47 = t504 + t521;
t43 = t602 + t495;
t42 = t603 + t496;
t30 = t283 + t536 + t540;
t28 = t502 + t542;
t22 = t493 - t777;
t17 = t487 + t776;
t15 = t497 + t535;
t1 = t488 - t505;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t148, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t79 + qJD(4) * t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t38 + qJD(4) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t585, -t486 * t662, 0, 0, 0, 0, 0, 0, 0, 0, -t483 * t585, t480 * t585, t594 * t662, t642 + t312 * qJD(3) + (-pkin(2) * t735 + qJ(3) * t594) * t662, 0, 0, 0, 0, 0, 0, t325 * qJD(4) - t527 * t585, t324 * qJD(4) + t451 * t585 (t394 * t451 + t395 * t527) * qJD(2), t647 + (t395 * t392 + t474 * t622 + t697) * qJD(2) + t169 * qJD(3), 0, 0, 0, 0, 0, 0 (-t358 * t527 + t362 * t394) * qJD(2) + t74 * qJD(4) (t359 * t527 + t394 * t692) * qJD(2) + t75 * qJD(4), t84 * qJD(4) + (-t702 - t703) * t628, t679 + (t231 * t358 + t232 * t359 + t697) * qJD(2) + t68 * qJD(3) + t17 * qJD(4) + t119 * qJD(5), 0, 0, 0, 0, 0, 0 (-t211 * t527 + t280 * t394) * qJD(2) + t42 * qJD(4) + t67 * qJD(6) (t212 * t527 + t354 * t394) * qJD(2) + t43 * qJD(4) + t66 * qJD(6) (-t211 * t354 - t212 * t280) * qJD(2) + t15 * qJD(4), t700 + (-t107 * t211 + t108 * t212 + t326 * t394) * qJD(2) + t30 * qJD(3) + t1 * qJD(4) + t50 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t30 + t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t325 - qJD(4) * t347, qJD(2) * t324 + qJD(4) * t346, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t74 - t347 * t653, qJD(2) * t75 + t347 * t654, t84 * qJD(2) - qJD(4) * t595, t680 + t17 * qJD(2) + (-pkin(4) * t347 - qJ(5) * t595) * qJD(4) + t129 * qJD(5), 0, 0, 0, 0, 0, 0, qJD(2) * t42 + qJD(6) * t143 + t347 * t657, qJD(2) * t43 + qJD(6) * t142 + t347 * t656, t15 * qJD(2) + (-t198 * t449 - t199 * t445) * qJD(4), t1 * qJD(2) + (-t198 * t389 + t199 * t391 + t347 * t473) * qJD(4) + t59 * qJD(5) + t560; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t119 + qJD(4) * t129, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t50 + qJD(4) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t67 + qJD(4) * t143 - qJD(6) * t173, qJD(2) * t66 + qJD(4) * t142 + qJD(6) * t172, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t311 - t642, 0, 0, 0, 0, 0, 0, 0, -t323 * qJD(4), 0, -qJD(3) * t168 - t647, 0, 0, 0, 0, 0, 0, t73 * qJD(4), t76 * qJD(4), t85 * qJD(4), qJD(3) * t69 + qJD(4) * t18 - qJD(5) * t118 - t679, 0, 0, 0, 0, 0, 0, -qJD(4) * t40 - qJD(6) * t64, -qJD(4) * t41 - qJD(6) * t65, qJD(4) * t16, -qJD(3) * t29 + qJD(4) * t2 - qJD(5) * t49 - t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t465 * qJD(3), t460 * qJD(3), t387, -t624 * qJD(4), 0, -t387, 0, 0, t474 * t627, t474 * t436, -qJD(3) * t766, qJD(3) * t210, t477 * t387, -0.2e1 * t427 * t360, t321 * qJD(4), t475 * t387, -t319 * qJD(4), -t387, -qJD(3) * t320 + qJD(4) * t77 + t482 * t606, -qJD(3) * t365 + qJD(4) * t78 - t479 * t606, -qJD(4) * t44 + qJD(5) * t263, qJD(3) * t91 + qJD(4) * t51 + qJD(5) * t115 (qJD(4) * t353 - t763) * t354, qJD(4) * t127 + qJD(6) * t167, qJD(4) * t188 + t280 * t651 (qJD(4) * t349 + qJD(6) * t354) * t280, qJD(4) * t186 + t354 * t651, -t387, qJD(3) * t187 + qJD(4) * t23 + qJD(6) * t61 - t652 * t693, qJD(3) * t189 + qJD(4) * t24 + qJD(6) * t60 - t280 * t652, qJD(3) * t126 + qJD(4) * t8 + qJD(5) * t125, qJD(3) * t31 + qJD(4) * t11 + qJD(5) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t625, -t550, 0, 0, 0, 0, 0, 0, 0, 0, -t631, -t551, 0, 0, 0, 0, 0, 0, -t636, -t663, 0, qJD(4) * t110 + qJD(5) * t309 - t563, 0, 0, 0, 0, 0, 0, -t237 + t668, -qJD(6) * t257 + t666, t671 (t349 * t445 + t353 * t449) * qJD(3) + t28 * qJD(4) + t164 * qJD(5) + t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t388, -t634, t436, -t388, -t627, 0, -qJD(4) * t392 + t474 * t628, qJD(4) * t390 + t474 * t629 - t672, 0, 0 -(-t477 * t628 - t607) * t527, -0.2e1 * t574 + (-t428 + t429) * qJD(4), t479 * t627 + t635 -(-t475 * t628 + t607) * t527, t427 - t637, -t388 (t479 * t578 - t370) * qJD(4) + t360 * qJD(5) + t562, t392 * t654 + (qJD(4) * t578 + t652) * t482 + t561, qJD(4) * t558 + t566, t110 * qJD(3) + (-pkin(4) * t392 + qJ(5) * t558) * qJD(4) + t114 * qJD(5) + t572, t775 + (t656 + t664) * t353, t648 + (-t694 - t707) * qJD(4) + t122, qJD(6) * t259 + t449 * t627 + t667, t349 * t764 - t775, -t445 * t627 - t237 + t669, t593 (t327 * t445 + t349 * t473 - t389 * t451) * qJD(4) + t255 * qJD(5) + t48 * qJD(6) + t571 (t327 * t449 + t353 * t473 - t391 * t451) * qJD(4) - t258 * qJD(5) + t47 * qJD(6) + t570 (-t111 * t449 - t112 * t445 - t349 * t391 + t353 * t389) * qJD(4) + t576 + t738, t28 * qJD(3) + (-t111 * t389 + t112 * t391 + t327 * t473) * qJD(4) + t22 * qJD(5) + t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t360 + t588 -(t479 * t628 - t653) * t527, t641, qJD(3) * t309 + qJD(4) * t114 - t552, 0, 0, 0, 0, 0, 0, qJD(6) * t222 - t611 + t658, -qJD(4) * t258 - t767, t555, qJD(3) * t164 + qJD(4) * t22 - t567; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t772, t549, qJD(4) * t259 + t280 * t591, -t772, t354 * t591 - t659, t430, -qJD(3) * t254 + qJD(4) * t48 + qJD(5) * t222 - qJD(6) * t108 - t564, -qJD(3) * t257 + qJD(4) * t47 + qJD(6) * t107 - t565, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t168 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t29 + t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t625, t550, 0, 0, 0, 0, 0, 0, t627, t436, t631, t551, 0, 0, 0, 0, 0, 0, t427 + t636, -qJD(4) * t362 + t663, t364 * qJD(4), -qJD(4) * t109 + qJD(5) * t310 + t563, 0, 0, 0, 0, 0, 0, -qJD(4) * t693 + t236 - t668, -qJD(4) * t280 - qJD(6) * t258 - t666, -qJD(4) * t94 - t671, -qJD(4) * t27 - qJD(5) * t163 - t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t628, t629, 0, 0, 0, 0, 0, 0, 0, 0, t604, -t633, t632, -t649, 0, 0, 0, 0, 0, 0, -t640, -t774, -t725, -t569; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t638, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t670; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, t435 - t643, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t73 * qJD(2), -t76 * qJD(2), -t85 * qJD(2), -qJD(2) * t18 - qJD(5) * t128 - t680, 0, 0, 0, 0, 0, 0, qJD(2) * t40 - qJD(6) * t140, qJD(2) * t41 - qJD(6) * t141, -qJD(2) * t16, -qJD(2) * t2 - qJD(5) * t58 - t560; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t388, t634, 0, t388, 0, 0, -t451 * t592, -t527 * t592 + t672, 0, 0, -t477 * t388, 0.2e1 * t574, -t635, -t475 * t388, t637, t388, -qJD(3) * t692 - t562, qJD(3) * t362 - t561, -qJD(3) * t364 - t566, qJD(3) * t109 - qJD(5) * t113 - t572, -t353 * t664 + t775, t122 - t648, -qJD(6) * t256 - t667, -t349 * t774 - t775, t236 - t669, -t593, qJD(3) * t693 - qJD(5) * t254 - qJD(6) * t45 - t571, qJD(3) * t280 - qJD(5) * t257 - qJD(6) * t46 - t570, qJD(3) * t94 - t576 + t738, qJD(3) * t27 - qJD(5) * t21 - t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t628, -t629, 0, 0, 0, 0, 0, 0, 0, 0, -t604, t633, -t632, t649, 0, 0, 0, 0, 0, 0, t640, t774, t725, t569; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t464 * qJD(5), t459 * qJD(5), -t605, t344 * qJD(6), 0, t605, 0, 0, t473 * t650, -t473 * t435, qJD(5) * t374, qJD(5) * t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t626, -t509, 0, 0, 0, 0, 0, 0, -t646, -t644, t554, -t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, t548, -t435 - t645, t773, t206, -t630, -qJD(6) * t391 + t449 * t655 - t557, qJD(6) * t389 - t445 * t655 - t556, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t118 + qJD(4) * t128, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t49 + qJD(4) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, t587, -t641, -qJD(3) * t310 + qJD(4) * t113 + t552, 0, 0, 0, 0, 0, 0, -qJD(6) * t219 + t611 + t659, qJD(4) * t257 - t763 + t767, -t555, qJD(3) * t163 + qJD(4) * t21 + t567; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t638, 0, 0, 0, 0, 0, 0, 0, 0, 0, t670; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t626, t509, 0, 0, 0, 0, 0, 0, t646 + t650, -t435 + t644, -t554, t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t547, -t764, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t64 + qJD(4) * t140, qJD(2) * t65 + qJD(4) * t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t772, -t549, qJD(4) * t256 - t767, t772, -t354 * t629 - t658, t430, -qJD(3) * t255 + qJD(4) * t45 + qJD(5) * t219 + t564, qJD(3) * t258 + qJD(4) * t46 + qJD(5) * t280 + t565, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t643, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, -t548, t645, -t773, -t235, t630, -t449 * t590 + t557, t445 * t590 + t556, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547, t764, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t3;
