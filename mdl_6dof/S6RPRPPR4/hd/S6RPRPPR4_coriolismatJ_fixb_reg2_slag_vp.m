% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPPR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:29
% EndTime: 2019-03-09 02:48:50
% DurationCPUTime: 13.56s
% Computational Cost: add. (13554->583), mult. (27054->752), div. (0->0), fcn. (31329->8), ass. (0->456)
t437 = sin(pkin(9));
t439 = cos(pkin(9));
t441 = sin(qJ(3));
t694 = cos(qJ(3));
t408 = t437 * t441 - t439 * t694;
t438 = cos(pkin(10));
t693 = cos(qJ(6));
t549 = t693 * t438;
t436 = sin(pkin(10));
t440 = sin(qJ(6));
t623 = t440 * t436;
t736 = (t623 / 0.2e1 + t549 / 0.2e1) * t408;
t553 = t694 * t437;
t621 = t441 * t439;
t411 = t553 + t621;
t403 = t411 ^ 2;
t712 = t408 ^ 2;
t720 = t403 + t712;
t730 = t720 * t438;
t735 = t730 * qJD(1);
t731 = t720 * t436;
t734 = t731 * qJD(1);
t733 = t731 * qJD(2);
t595 = qJD(4) * t408;
t543 = t411 * t595;
t732 = qJD(2) * t730 + t436 * t543;
t283 = t438 * t411;
t574 = t283 * qJD(1);
t597 = qJD(3) * t436;
t727 = t574 + t597;
t550 = t693 * t436;
t622 = t440 * t438;
t410 = t550 - t622;
t598 = qJD(3) * t410;
t406 = t549 + t623;
t631 = t411 * t406;
t715 = t631 * qJD(1);
t484 = t598 + t715;
t360 = t411 * t550;
t261 = t411 * t622 - t360;
t703 = -t410 / 0.2e1;
t705 = -t406 / 0.2e1;
t721 = t631 * t705;
t108 = t261 * t703 + t721;
t726 = t108 * qJD(6);
t725 = -qJD(1) * t108 + t406 * t598;
t724 = -qJD(3) * t108 + t261 * t715;
t723 = t631 ^ 2;
t722 = t631 * t410;
t569 = t408 * qJD(1);
t548 = t631 * t569;
t555 = t410 * t693;
t384 = t555 / 0.2e1;
t696 = t440 / 0.2e1;
t462 = t406 * t696 + t384;
t432 = t436 ^ 2;
t434 = t438 ^ 2;
t424 = t432 + t434;
t418 = t424 * qJ(4);
t565 = t418 * qJD(3);
t689 = pkin(7) + qJ(2);
t522 = t689 * t437;
t508 = t441 * t522;
t420 = t689 * t439;
t554 = t694 * t420;
t446 = t554 / 0.2e1 - t508 / 0.2e1;
t337 = t554 - t508;
t312 = t438 * t337;
t682 = qJ(4) * t411;
t501 = pkin(3) * t408 - t682;
t557 = t439 * pkin(2) + pkin(1);
t455 = t501 - t557;
t154 = t436 * t455 + t312;
t659 = t154 * t438;
t628 = t436 * t337;
t153 = t438 * t455 - t628;
t660 = t153 * t436;
t477 = -t660 / 0.2e1 + t659 / 0.2e1;
t79 = t446 - t477;
t719 = -qJD(1) * t79 + t565;
t627 = t438 * qJ(5);
t530 = t627 / 0.2e1;
t690 = t436 * pkin(4);
t443 = (t530 - t690 / 0.2e1) * t408 + t446;
t135 = -t408 * pkin(4) - t153;
t661 = t135 * t436;
t134 = t408 * qJ(5) + t154;
t662 = t134 * t438;
t479 = t662 / 0.2e1 + t661 / 0.2e1;
t43 = t443 - t479;
t718 = -qJD(1) * t43 + t565;
t567 = t411 * qJD(3);
t375 = t438 * t567;
t258 = t712 - t403;
t235 = t258 * t436;
t582 = t235 * qJD(1);
t717 = -t582 - t375;
t716 = qJD(6) * t631;
t613 = -t621 / 0.2e1 - t553 / 0.2e1;
t710 = pkin(4) + pkin(5);
t714 = t436 * t710 - t627;
t713 = t406 ^ 2;
t711 = t410 ^ 2;
t688 = -pkin(8) + qJ(4);
t419 = t688 * t438;
t521 = t688 * t436;
t334 = t419 * t440 - t521 * t693;
t709 = -t334 / 0.2e1;
t336 = t419 * t693 + t440 * t521;
t708 = t336 / 0.2e1;
t359 = t408 * t550;
t707 = -t359 / 0.2e1;
t520 = t436 * qJ(5) + pkin(3);
t390 = t438 * t710 + t520;
t706 = -t390 / 0.2e1;
t704 = t406 / 0.2e1;
t702 = t410 / 0.2e1;
t701 = -t432 / 0.2e1;
t700 = -t436 / 0.2e1;
t699 = t436 / 0.2e1;
t698 = t438 / 0.2e1;
t697 = -t440 / 0.2e1;
t198 = t410 * t411;
t71 = -t722 + (-t261 / 0.2e1 + t198 / 0.2e1) * t406;
t695 = t71 * qJD(4);
t692 = t411 * pkin(3);
t691 = t411 * pkin(4);
t279 = t436 * t411;
t115 = pkin(8) * t279 + t134;
t444 = t628 + (t411 * t688 + t557) * t438 + (-pkin(3) * t438 - t710) * t408;
t60 = t440 * t115 - t444 * t693;
t684 = t60 * t410;
t47 = -t684 / 0.2e1;
t558 = t684 / 0.2e1;
t9 = t47 + t558;
t687 = qJD(2) * t9;
t335 = t441 * t420 + t694 * t522;
t147 = -t411 * t714 - t335;
t148 = t408 * t714 - t337;
t61 = t115 * t693 + t440 * t444;
t310 = t436 * t335;
t635 = t408 * qJ(4);
t314 = t635 + t692;
t101 = -t310 + (pkin(8) * t408 - t314) * t438 - t710 * t411;
t552 = t693 * t101;
t311 = t438 * t335;
t181 = t436 * t314 - t311;
t146 = t411 * qJ(5) + t181;
t278 = t436 * t408;
t119 = -pkin(8) * t278 + t146;
t624 = t440 * t119;
t66 = t552 - t624;
t551 = t693 * t119;
t625 = t440 * t101;
t67 = t551 + t625;
t4 = t147 * t148 - t60 * t66 + t61 * t67;
t686 = t4 * qJD(1);
t260 = -t408 * t622 + t359;
t636 = t406 * t408;
t5 = -t260 * t61 - t261 * t67 - t60 * t636 - t631 * t66;
t685 = t5 * qJD(1);
t447 = t260 * t334 / 0.2e1 - t636 * t708 + t411 * t706;
t480 = t66 * t704 + t67 * t703;
t18 = t447 + t480;
t681 = qJD(1) * t18;
t20 = t147 * t283 + (t440 * t60 + t61 * t693) * t408;
t680 = qJD(1) * t20;
t21 = -t198 * t61 - t60 * t631;
t679 = qJD(1) * t21;
t31 = t147 * t631 + t408 * t61;
t678 = qJD(1) * t31;
t500 = -t627 + t690;
t184 = t411 * t500 + t335;
t36 = t184 * t411 + (-t661 - t662) * t408;
t677 = qJD(1) * t36;
t639 = t335 * t411;
t68 = t639 + (-t659 + t660) * t408;
t675 = qJD(1) * t68;
t633 = t410 * t260;
t646 = t636 * t406;
t70 = t633 + t646;
t674 = qJD(1) * t70;
t73 = (-t134 * t436 + t135 * t438) * t411;
t673 = qJD(1) * t73;
t663 = t134 * t408;
t78 = -t184 * t283 + t663;
t672 = qJD(1) * t78;
t490 = -t153 * t438 - t154 * t436;
t81 = t490 * t411;
t670 = qJD(1) * t81;
t648 = t261 * t636;
t650 = t260 * t631;
t89 = t648 + t650;
t669 = qJD(1) * t89;
t560 = -t693 / 0.2e1;
t453 = (t406 * t560 + t410 * t696) * t408;
t459 = t260 * t697 - t560 * t636;
t92 = t453 - t459;
t668 = qJD(1) * t92;
t454 = t462 * t408;
t559 = t693 / 0.2e1;
t463 = t260 * t559 - t636 * t697;
t94 = t454 + t463;
t667 = qJD(1) * t94;
t475 = -t198 * t702 + t721;
t95 = -t613 - t475;
t666 = qJD(1) * t95;
t11 = t147 * t260 + t148 * t261 - t408 * t66 + t411 * t60;
t665 = t11 * qJD(1);
t12 = -t147 * t636 + t148 * t631 + t408 * t67 + t411 * t61;
t664 = t12 * qJD(1);
t17 = -t147 * t411 + t260 * t60 - t61 * t636;
t658 = t17 * qJD(1);
t626 = t438 * t314;
t180 = t310 + t626;
t657 = t180 * t438;
t656 = t181 * t436;
t655 = t184 * t436;
t185 = -t408 * t500 + t337;
t654 = t185 * t436;
t149 = -t180 - t691;
t24 = t134 * t146 + t135 * t149 + t184 * t185;
t653 = t24 * qJD(1);
t448 = t147 * t702 + t408 * t708 + t390 * t631 / 0.2e1;
t460 = -t624 / 0.2e1 + t552 / 0.2e1;
t25 = -t448 + t460;
t652 = t25 * qJD(1);
t449 = t147 * t705 + t261 * t706 + t408 * t709;
t461 = -t625 / 0.2e1 - t551 / 0.2e1;
t26 = -t449 + t461;
t651 = t26 * qJD(1);
t649 = t260 * t408;
t647 = t261 * t411;
t645 = t636 * t408;
t644 = t631 * t411;
t282 = t438 * t408;
t29 = t135 * t282 - t149 * t283 + (t146 * t411 - t663) * t436;
t643 = t29 * qJD(1);
t30 = -t147 * t261 - t408 * t60;
t642 = t30 * qJD(1);
t32 = (-t185 * t438 + t134) * t411 + (t184 * t438 + t146) * t408;
t641 = t32 * qJD(1);
t33 = (-t135 + t654) * t411 + (-t149 - t655) * t408;
t640 = t33 * qJD(1);
t34 = (t656 + t657) * t411 + t490 * t408;
t638 = t34 * qJD(1);
t35 = t153 * t180 + t154 * t181 + t335 * t337;
t637 = t35 * qJD(1);
t41 = (t153 + t628) * t411 + (t180 - t310) * t408;
t634 = t41 * qJD(1);
t632 = t410 * t408;
t416 = -pkin(4) * t438 - t520;
t630 = t411 * t416;
t42 = (-t154 + t312) * t411 + (-t181 - t311) * t408;
t629 = t42 * qJD(1);
t524 = t701 - t434 / 0.2e1;
t482 = t524 * t635;
t450 = t482 + t630 / 0.2e1;
t478 = t146 * t699 - t149 * t438 / 0.2e1;
t58 = t450 - t478;
t620 = t58 * qJD(1);
t474 = t635 / 0.2e1 - t630 / 0.2e1;
t504 = t310 / 0.2e1 + t691 / 0.2e1;
t536 = -t655 / 0.2e1;
t75 = t536 + (t314 / 0.2e1 + t474) * t438 + t504;
t619 = t75 * qJD(1);
t452 = t482 - t692 / 0.2e1;
t476 = t657 / 0.2e1 + t656 / 0.2e1;
t76 = t452 - t476;
t618 = t76 * qJD(1);
t90 = t648 - t650;
t617 = t90 * qJD(1);
t425 = t437 ^ 2 + t439 ^ 2;
t102 = (-t261 * t693 + t440 * t631) * t408;
t612 = qJD(1) * t102;
t109 = t647 + t649;
t611 = qJD(1) * t109;
t110 = -t647 + t649;
t610 = qJD(1) * t110;
t111 = -t644 - t645;
t609 = qJD(1) * t111;
t112 = -t644 + t645;
t608 = qJD(1) * t112;
t132 = -t261 * t283 - t440 * t712;
t607 = qJD(1) * t132;
t606 = qJD(1) * t261;
t289 = t384 - t555 / 0.2e1;
t602 = qJD(2) * t289;
t601 = qJD(2) * t411;
t385 = t432 * t408;
t386 = t434 * t408;
t600 = qJD(3) * (-t385 + t386);
t599 = qJD(3) * t406;
t596 = qJD(3) * t438;
t594 = qJD(5) * t436;
t593 = qJD(6) * t408;
t141 = -t337 * t408 + t639;
t592 = t141 * qJD(1);
t152 = t283 * t631 + t693 * t712;
t591 = t152 * qJD(1);
t525 = t622 / 0.2e1;
t168 = t707 + (t525 + t702) * t408;
t590 = t168 * qJD(1);
t169 = t707 + (t525 + t703) * t408;
t157 = t169 * qJD(1);
t535 = -t636 / 0.2e1;
t170 = t535 - t736;
t589 = t170 * qJD(1);
t497 = t535 + t736;
t588 = t497 * qJD(1);
t534 = t636 / 0.2e1;
t496 = t534 + t736;
t587 = t496 * qJD(1);
t189 = t424 * t403;
t586 = t189 * qJD(1);
t584 = t198 * qJD(1);
t503 = t524 * t411;
t230 = t503 + t613;
t583 = t230 * qJD(1);
t581 = t235 * qJD(3);
t236 = t258 * t438;
t234 = t236 * qJD(1);
t578 = t258 * qJD(1);
t577 = t278 * qJD(1);
t576 = t279 * qJD(1);
t575 = t282 * qJD(1);
t286 = t440 * t408;
t573 = t286 * qJD(1);
t293 = t385 + t386;
t572 = t293 * qJD(1);
t571 = t720 * qJD(1);
t570 = t613 * qJD(1);
t393 = t406 * qJD(6);
t394 = t408 * qJD(3);
t395 = t410 * qJD(6);
t568 = t411 * qJD(1);
t417 = t425 * qJ(2);
t566 = t417 * qJD(1);
t564 = t424 * qJD(3);
t563 = t425 * qJD(1);
t556 = t408 * t693;
t546 = t198 * t569;
t544 = t436 * t596;
t542 = t408 * t568;
t541 = t408 * t567;
t540 = t406 * t395;
t539 = t438 * t594;
t538 = t436 * t568;
t537 = t438 * t568;
t532 = t261 * t699;
t531 = t631 * t699;
t529 = -t283 / 0.2e1;
t528 = t283 / 0.2e1;
t523 = qJD(6) * t693;
t519 = -qJD(6) * t613 - t542;
t518 = -qJD(1) * t557 + qJD(2);
t517 = -qJD(6) + t569;
t516 = qJD(3) * t390 - qJD(4);
t515 = qJD(1) * t403 * t436 * t438;
t514 = t436 * t375;
t512 = t432 * t542;
t511 = t434 * t542;
t510 = t408 * t538;
t509 = t408 * t537;
t507 = qJD(1) * t556;
t502 = (t701 + t434 / 0.2e1) * t411;
t499 = t436 * t509;
t498 = t408 * t514;
t495 = t408 * t416 + t682;
t88 = t198 * t261 - t723;
t494 = qJD(1) * t88 + qJD(3) * t71;
t87 = t406 * t261 - t722;
t98 = t261 ^ 2 - t723;
t493 = qJD(1) * t98 + qJD(3) * t87;
t492 = qJD(1) * t9 + qJD(5) * t289;
t491 = t146 * t438 + t149 * t436;
t489 = -t180 * t436 + t181 * t438;
t442 = (t530 + (-pkin(4) / 0.2e1 - pkin(5) / 0.2e1) * t436) * t408 + t446;
t451 = -t198 * t708 + t61 * t704 + t631 * t709;
t13 = t558 + t442 - t451;
t140 = -t410 * t334 + t336 * t406;
t488 = -qJD(1) * t13 + qJD(3) * t140;
t199 = -t711 - t713;
t487 = qJD(1) * t71 + qJD(3) * t199;
t257 = -t711 + t713;
t486 = qJD(1) * t87 + qJD(3) * t257;
t464 = -t198 * t697 + t560 * t631;
t114 = t529 + t464;
t268 = t700 - t462;
t485 = qJD(1) * t114 + qJD(3) * t268;
t143 = -t360 / 0.2e1 + (t622 - t550 / 0.2e1) * t411;
t483 = qJD(1) * t143 + t599;
t481 = t408 * t525 + t707;
t445 = (t334 * t696 + t336 * t559) * t408 + t147 * t699 + t390 * t528;
t465 = t560 * t66 + t67 * t697;
t16 = t445 + t465;
t473 = qJD(1) * t16 + t390 * t597;
t124 = t532 + (t406 * t698 + t559) * t411;
t468 = qJD(1) * t124 + t406 * t597;
t126 = t531 + (t410 * t698 + t697) * t411;
t467 = qJD(1) * t126 + t410 * t597;
t227 = t502 - t613;
t466 = qJD(1) * t227 + t544;
t458 = -qJD(3) * t282 + t510;
t301 = t403 * t434 + t712;
t457 = qJD(1) * t301 + t514;
t421 = t424 * qJD(4);
t413 = t418 * qJD(4);
t387 = t613 * qJD(3);
t374 = t436 * t567;
t340 = t411 * t539;
t338 = qJD(3) * t432 + t436 * t537;
t305 = t434 * t541;
t304 = t432 * t541;
t291 = 0.2e1 * t499;
t290 = -0.2e1 * t499;
t275 = t293 * qJD(2);
t273 = t293 * qJD(3);
t272 = t289 * qJD(6);
t267 = t700 + t462;
t233 = t236 * qJD(3);
t229 = t503 - t613;
t228 = t502 + t613;
t225 = t230 * qJD(2);
t223 = t230 * qJD(4);
t222 = t229 * qJD(2);
t221 = t229 * qJD(4);
t220 = (-t434 * t568 - t544) * t408;
t219 = (-t432 * t568 + t544) * t408;
t187 = t189 * qJD(4);
t178 = t374 - t234;
t177 = t632 / 0.2e1 + t481;
t176 = -t632 / 0.2e1 + t481;
t175 = t534 - t736;
t159 = t177 * qJD(6);
t158 = t169 * qJD(6);
t133 = -t395 - t157;
t125 = t410 * t528 + t411 * t696 + t531;
t123 = t406 * t528 + t411 * t560 + t532;
t113 = t529 - t464;
t96 = -t613 + t475;
t93 = t454 - t463;
t91 = t453 + t459;
t86 = t87 * qJD(6);
t80 = t446 + t477;
t77 = t452 + t476;
t74 = t536 - t626 / 0.2e1 + t474 * t438 - t504;
t59 = t450 + t478;
t44 = t443 + t479;
t28 = t448 + t460;
t27 = t449 + t461;
t19 = t447 - t480;
t15 = t445 - t465;
t14 = t47 + t442 + t451;
t8 = t9 * qJD(6);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425 * qJD(2), t417 * qJD(2), -t541, t258 * qJD(3), 0, t541, 0, 0, -t557 * t567, t557 * t394, qJD(2) * t720, qJD(2) * t141, -t305, 0.2e1 * t498, -t233, -t304, t581, t541, qJD(3) * t41 - t438 * t543 + t733, qJD(3) * t42 + t732, -qJD(3) * t34 + t187, qJD(2) * t68 + qJD(3) * t35 + qJD(4) * t81, -t305, -t233, -0.2e1 * t498, t541, -t581, -t304, t733 + t33 * qJD(3) + (-t403 * t594 - t543) * t438, -t408 * t411 * t594 - qJD(3) * t29 + t187, qJD(3) * t32 + qJD(5) * t301 - t732, qJD(2) * t36 + qJD(3) * t24 + qJD(4) * t73 + qJD(5) * t78 (-qJD(3) * t636 - qJD(6) * t261) * t631, qJD(3) * t90 + qJD(6) * t98, qJD(3) * t112 + t261 * t593 (qJD(3) * t260 + t716) * t261, qJD(3) * t109 + t593 * t631, t541, qJD(2) * t110 + qJD(3) * t11 - qJD(5) * t132 + qJD(6) * t31 - t595 * t631, qJD(2) * t111 + qJD(3) * t12 + qJD(5) * t152 + qJD(6) * t30 - t198 * t595, qJD(2) * t89 + qJD(3) * t5 + qJD(4) * t88 + qJD(5) * t102, qJD(2) * t17 + qJD(3) * t4 + qJD(4) * t21 + qJD(5) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t563, t566, 0, 0, 0, 0, 0, 0, 0, 0, t571, t592, 0, 0, 0, 0, 0, 0, t734, t735, 0, qJD(3) * t77 + t221 + t675, 0, 0, 0, 0, 0, 0, t734, 0, -t735, qJD(3) * t59 + t221 + t677, 0, 0, 0, 0, 0, 0, t159 + t610, qJD(6) * t497 + t609, t669, t658 + (t260 * t406 - t410 * t636) * qJD(2) + t19 * qJD(3) + t96 * qJD(4) + t93 * qJD(5) + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t542, t578, -t394, t542, -t567, 0, -qJD(3) * t337 - t557 * t568, qJD(3) * t335 + t557 * t569, 0, 0, t220, t291 - t600, t178, t219, -t717, t542, t634 + (t436 * t501 - t312) * qJD(3) - t278 * qJD(4), t337 * t597 + t629 + (qJD(3) * t501 - t595) * t438, qJD(3) * t489 - t638, t637 + t77 * qJD(2) + (-t337 * pkin(3) + qJ(4) * t489) * qJD(3) + t80 * qJD(4), t220, t178, t290 + t600, t542, t717, t219, -t185 * t596 + t640 + t228 * qJD(5) + (-qJD(3) * t495 - t595) * t436, qJD(3) * t491 - t643, t641 + (t438 * t495 - t654) * qJD(3) + t282 * qJD(4) + t340, t653 + t59 * qJD(2) + (qJ(4) * t491 + t185 * t416) * qJD(3) + t44 * qJD(4) + t74 * qJD(5), -t484 * t636 + t726, t617 + (-t633 + t646) * qJD(3) + t86, qJD(6) * t175 - t410 * t567 + t608, -t726 + (t599 + t606) * t260, t406 * t567 + t159 + t611, -t519, t665 + (t148 * t406 + t260 * t390 + t334 * t411) * qJD(3) + t176 * qJD(4) + t123 * qJD(5) + t28 * qJD(6), t664 + (t148 * t410 + t336 * t411 - t390 * t636) * qJD(3) + t496 * qJD(4) + t125 * qJD(5) + t27 * qJD(6), t685 + (-t260 * t336 - t334 * t636 - t406 * t67 - t410 * t66) * qJD(3) + t91 * qJD(5) + t695, t686 + t19 * qJD(2) + (t148 * t390 - t334 * t66 + t336 * t67) * qJD(3) + t14 * qJD(4) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t278 - t509 (t538 - t596) * t408, t586, qJD(3) * t80 + t222 + t670, 0, 0, 0, 0, 0, 0 (-t537 - t597) * t408, t586, -t458, qJD(3) * t44 + t222 + t673, 0, 0, 0, 0, 0, 0, qJD(3) * t176 - t548, qJD(3) * t496 - t546, t494, qJD(2) * t96 + qJD(3) * t14 + qJD(5) * t113 + t679; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t228 - t515, -t510, t457, qJD(3) * t74 + t672, 0, 0, 0, 0, 0, 0, qJD(3) * t123 - t607, qJD(3) * t125 + t591, qJD(3) * t91 + t612, qJD(2) * t93 + qJD(3) * t15 + qJD(4) * t113 + t680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t724, t493, t175 * qJD(3) + t261 * t517, t724, t177 * qJD(3) + t517 * t631, t387, qJD(2) * t177 + qJD(3) * t28 - qJD(6) * t61 + t678, qJD(2) * t497 + qJD(3) * t27 + qJD(6) * t60 + t642, 0, t687; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t563, -t566, 0, 0, 0, 0, 0, 0, t567, -t394, -t571, -t592, 0, 0, 0, 0, 0, 0, t375 - t734, -qJD(3) * t279 - t735, t273, -qJD(3) * t76 + t223 - t675, 0, 0, 0, 0, 0, 0, qJD(3) * t283 - t734, t273, t374 + t735, -qJD(3) * t58 + qJD(5) * t278 + t223 - t677, 0, 0, 0, 0, 0, 0, qJD(3) * t631 - t158 - t610, qJD(3) * t198 - qJD(6) * t496 - t609, -qJD(3) * t70 - t669, -qJD(3) * t18 - qJD(4) * t95 + qJD(5) * t94 - t658 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, -t569, 0, 0, 0, 0, 0, 0, 0, 0, t537, -t576, t572, -t618, 0, 0, 0, 0, 0, 0, t574, t572, t538, -t620, 0, 0, 0, 0, 0, 0, t715, t584, -t674, -t681; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583, 0, 0, 0, 0, 0, 0, 0, 0, 0, t583, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272 + t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t393 - t587, 0, t492; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t542, -t578, 0, -t542, 0, 0, -t518 * t411, t518 * t408, 0, 0, t511, t290, t234, t512, -t582, -t542, -t438 * t601 - t634, qJD(2) * t279 - t629, -t275 + t638, qJD(2) * t76 - qJD(4) * t79 - t637, t511, t234, t291, -t542, t582, t512, -qJD(2) * t283 + qJD(5) * t227 - t640, qJD(5) * t282 - t275 + t643, -t436 * t601 + t340 - t641, qJD(2) * t58 - qJD(4) * t43 + qJD(5) * t75 - t653, t636 * t715 + t726, t86 - t617, -qJD(6) * t170 - t608, -t260 * t606 - t726, -t158 - t611, t519, -qJD(2) * t631 - qJD(4) * t168 + qJD(5) * t124 - qJD(6) * t25 - t665, -qJD(2) * t198 - qJD(4) * t497 + qJD(5) * t126 - qJD(6) * t26 - t664, qJD(2) * t70 + qJD(5) * t92 - t685 + t695, qJD(2) * t18 - qJD(4) * t13 + qJD(5) * t16 - t686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t568, t569, 0, 0, 0, 0, 0, 0, 0, 0, -t537, t576, -t572, t618, 0, 0, 0, 0, 0, 0, -t574, -t572, -t538, t620, 0, 0, 0, 0, 0, 0, -t715, -t584, t674, t681; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, t413, 0, 0, 0, 0, 0, 0, t539, t421, t432 * qJD(5), -t416 * t594 + t413, -t540, t257 * qJD(6), 0, t540, 0, 0, t390 * t395 + t406 * t594, -t390 * t393 + t410 * t594, qJD(4) * t199, qJD(4) * t140 + t390 * t594; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t564, t719, 0, 0, 0, 0, 0, 0, 0, t564, 0, t718, 0, 0, 0, 0, 0, 0, -t590, -t588, t487, qJD(5) * t267 + t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t466, t575, t338, -t416 * t597 + t619, 0, 0, 0, 0, 0, 0, t468, t467, t668, qJD(4) * t267 + t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t725, t486, -t393 - t589, t725, t133, -t570, -qJD(6) * t336 + t390 * t598 - t652, qJD(6) * t334 - t390 * t599 - t651, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t509, -t510, -t586, qJD(3) * t79 - t225 - t670, 0, 0, 0, 0, 0, 0, t509, -t586, t510, qJD(3) * t43 - qJD(5) * t283 - t225 - t673, 0, 0, 0, 0, 0, 0, qJD(3) * t168 + t548 - t716, qJD(3) * t497 + qJD(6) * t143 + t546, -t494, qJD(2) * t95 + qJD(3) * t13 + qJD(5) * t114 - t679; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t583, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t583, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t564, -t719, 0, 0, 0, 0, 0, 0, 0, -t564, 0, -t594 - t718, 0, 0, 0, 0, 0, 0, -t395 + t590, t393 + t588, -t487, qJD(5) * t268 - t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t727, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, t483, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t227 + t515, t458, -t457, -qJD(2) * t278 - qJD(3) * t75 + qJD(4) * t283 - t672, 0, 0, 0, 0, 0, 0, -qJD(3) * t124 + qJD(6) * t286 + t607, -t126 * qJD(3) + t408 * t523 - t591, -qJD(3) * t92 - t612, -qJD(2) * t94 - qJD(3) * t16 - qJD(4) * t114 - t680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t577, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272 - t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t466, -t575, -t338, -t619 + (qJD(3) * t416 + qJD(4)) * t436, 0, 0, 0, 0, 0, 0, -t468, -t467, -t668, -qJD(4) * t268 - t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t727, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t440 + t573, t507 - t523, 0, t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t724, -t493, qJD(3) * t170 - t261 * t569, -t724, qJD(3) * t169 - t548, t387, qJD(2) * t169 + qJD(3) * t25 + qJD(4) * t631 - qJD(5) * t286 - t678, qJD(2) * t496 + t26 * qJD(3) - t143 * qJD(4) - qJD(5) * t556 - t642, 0, -t687; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t587, 0, -t492; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t725, -t486, t589, -t725, t157, t570, -t410 * t516 + t652, t406 * t516 + t651, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t484, -t483, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t573, -t507, 0, -t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
