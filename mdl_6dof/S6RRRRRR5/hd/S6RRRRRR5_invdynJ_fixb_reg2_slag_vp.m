% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:58
% EndTime: 2019-03-10 04:06:05
% DurationCPUTime: 41.05s
% Computational Cost: add. (58367->1042), mult. (141641->1371), div. (0->0), fcn. (115465->18), ass. (0->457)
t428 = cos(qJ(2));
t424 = sin(qJ(2));
t418 = sin(pkin(6));
t588 = qJD(1) * t418;
t545 = t424 * t588;
t419 = cos(pkin(6));
t587 = qJD(1) * t419;
t561 = pkin(1) * t587;
t317 = -pkin(8) * t545 + t428 * t561;
t496 = pkin(2) * t424 - pkin(9) * t428;
t318 = t496 * t588;
t423 = sin(qJ(3));
t427 = cos(qJ(3));
t224 = -t317 * t423 + t427 * t318;
t429 = -pkin(10) - pkin(9);
t551 = qJD(3) * t429;
t598 = t427 * t428;
t733 = -(pkin(3) * t424 - pkin(10) * t598) * t588 - t224 + t427 * t551;
t225 = t427 * t317 + t423 * t318;
t544 = t428 * t588;
t505 = t423 * t544;
t732 = -pkin(10) * t505 - t423 * t551 + t225;
t422 = sin(qJ(4));
t661 = cos(qJ(4));
t347 = t422 * t427 + t423 * t661;
t275 = t347 * t544;
t677 = qJD(3) + qJD(4);
t437 = t677 * t347;
t682 = t275 - t437;
t476 = t422 * t423 - t427 * t661;
t276 = t476 * t544;
t439 = t677 * t476;
t680 = t276 - t439;
t377 = t429 * t423;
t378 = t429 * t427;
t538 = qJD(4) * t661;
t582 = qJD(4) * t422;
t686 = t377 * t538 + t378 * t582 + t422 * t733 - t732 * t661;
t731 = t732 * t422 + t661 * t733;
t584 = qJD(3) * t423;
t730 = t505 - t584;
t426 = cos(qJ(6));
t577 = qJD(6) * t426;
t421 = sin(qJ(5));
t513 = qJD(2) + t587;
t293 = t423 * t545 - t427 * t513;
t295 = t423 * t513 + t427 * t545;
t474 = -t422 * t293 + t295 * t661;
t475 = -t293 * t661 - t422 * t295;
t660 = cos(qJ(5));
t147 = t421 * t474 - t475 * t660;
t722 = t147 * t426;
t729 = t577 + t722;
t728 = -pkin(4) * t545 - pkin(11) * t680 - t377 * t582 + t378 * t538 + t731;
t727 = -pkin(11) * t682 - t686;
t420 = sin(qJ(6));
t376 = -qJD(3) + t544;
t362 = -qJD(4) + t376;
t464 = -qJD(5) + t362;
t451 = t426 * t464;
t663 = t421 * t475 + t474 * t660;
t133 = t420 * t663 + t451;
t135 = -t420 * t464 + t426 * t663;
t572 = qJD(1) * qJD(2);
t535 = t424 * t572;
t501 = t418 * t535;
t569 = qJDD(1) * t428;
t392 = t418 * t569;
t568 = t392 - qJDD(3);
t314 = t501 - t568;
t302 = qJDD(4) + t314;
t448 = qJDD(5) + t302;
t578 = qJD(6) * t420;
t534 = t428 * t572;
t570 = qJDD(1) * t424;
t468 = t534 + t570;
t481 = qJD(3) * t513;
t442 = t418 * t468 + t481;
t571 = qJDD(1) * t419;
t393 = qJDD(2) + t571;
t502 = qJD(3) * t545;
t518 = (t393 - t502) * t423;
t433 = t442 * t427 + t518;
t533 = t418 * t570;
t708 = t418 * t534 + t533;
t703 = t481 + t708;
t507 = t423 * t703 + t427 * t502;
t469 = t393 * t427 - t507;
t115 = t293 * t538 + t295 * t582 - t422 * t469 - t661 * t433;
t706 = t422 * t433 - t661 * t469;
t116 = qJD(4) * t474 + t706;
t537 = qJD(5) * t660;
t581 = qJD(5) * t421;
t61 = t660 * t115 + t421 * t116 + t474 * t581 - t475 * t537;
t41 = qJD(6) * t451 - t420 * t448 + t426 * t61 + t578 * t663;
t580 = qJD(6) * t135;
t42 = -t420 * t61 - t426 * t448 + t580;
t719 = qJD(6) + t147;
t726 = t719 * t420;
t724 = -t133 * t729 - t135 * t726 - t41 * t426 - t420 * t42;
t524 = t421 * t115 - t660 * t116;
t62 = qJD(5) * t663 - t524;
t60 = qJDD(6) + t62;
t58 = t426 * t60;
t723 = t133 * t663 - t719 * t726 + t58;
t449 = t660 * t476;
t155 = qJD(5) * t449 + t347 * t581 + t421 * t437 + t439 * t660;
t194 = -t421 * t275 - t276 * t660;
t597 = t155 + t194;
t458 = t421 * t476;
t596 = -qJD(5) * t458 + t347 * t537 + t421 * t680 - t660 * t682;
t320 = pkin(8) * t544 + t424 * t561;
t684 = pkin(3) * t730 + t320;
t40 = t42 * t426;
t629 = t133 * t420;
t725 = t629 * t719 - t40;
t38 = t41 * t420;
t715 = t135 * t729 - t38;
t714 = -t135 * t663 + t420 * t60 + t719 * t729;
t281 = t661 * t377 + t422 * t378;
t228 = -pkin(11) * t347 + t281;
t282 = t422 * t377 - t661 * t378;
t229 = -pkin(11) * t476 + t282;
t640 = t228 * t537 - t229 * t581 + t421 * t728 - t727 * t660;
t274 = pkin(9) * t513 + t320;
t482 = -pkin(2) * t428 - pkin(9) * t424 - pkin(1);
t309 = t482 * t418;
t286 = qJD(1) * t309;
t197 = -t274 * t423 + t427 * t286;
t175 = -pkin(10) * t295 + t197;
t163 = -pkin(3) * t376 + t175;
t198 = t274 * t427 + t286 * t423;
t176 = -pkin(10) * t293 + t198;
t547 = t661 * t176;
t111 = t422 * t163 + t547;
t694 = pkin(11) * t475;
t96 = t111 + t694;
t634 = t421 * t96;
t173 = t422 * t176;
t110 = t661 * t163 - t173;
t695 = pkin(11) * t474;
t95 = t110 - t695;
t89 = -pkin(4) * t362 + t95;
t53 = t660 * t89 - t634;
t51 = pkin(5) * t464 - t53;
t712 = t147 * t51;
t559 = t660 * t96;
t54 = t421 * t89 + t559;
t52 = -pkin(12) * t464 + t54;
t273 = -pkin(2) * t513 - t317;
t214 = t293 * pkin(3) + t273;
t157 = -pkin(4) * t475 + t214;
t75 = t147 * pkin(5) - pkin(12) * t663 + t157;
t486 = t420 * t52 - t426 * t75;
t721 = t147 * t486;
t623 = t147 * t663;
t595 = -pkin(4) * t682 - t684;
t704 = -t147 ^ 2 + t663 ^ 2;
t21 = t420 * t75 + t426 * t52;
t509 = qJD(2) * t561;
t558 = pkin(1) * t571;
t506 = pkin(8) * t708 + t424 * t509 - t428 * t558;
t644 = t393 * pkin(2);
t219 = t506 - t644;
t154 = -pkin(3) * t469 + t219;
t82 = t116 * pkin(4) + t154;
t16 = t62 * pkin(5) + t61 * pkin(12) + t82;
t552 = -pkin(8) * t392 - t424 * t558 - t428 * t509;
t234 = -pkin(8) * t501 - t552;
t218 = pkin(9) * t393 + t234;
t477 = t496 * qJD(2);
t222 = (qJD(1) * t477 + qJDD(1) * t482) * t418;
t583 = qJD(3) * t427;
t467 = -t427 * t218 - t423 * t222 + t274 * t584 - t286 * t583;
t101 = pkin(10) * t469 - t467;
t521 = -t423 * t218 + t427 * t222;
t92 = t314 * pkin(3) - pkin(10) * t433 - t274 * t583 - t286 * t584 + t521;
t36 = -qJD(4) * t111 - t422 * t101 + t661 * t92;
t19 = t302 * pkin(4) + t115 * pkin(11) + t36;
t514 = -t661 * t101 - t163 * t538 + t176 * t582 - t422 * t92;
t23 = -pkin(11) * t116 - t514;
t9 = t421 * t19 + t660 * t23 + t89 * t537 - t96 * t581;
t7 = pkin(12) * t448 + t9;
t3 = -qJD(6) * t21 + t426 * t16 - t420 * t7;
t718 = -t147 * t21 - t3;
t100 = pkin(5) * t663 + pkin(12) * t147;
t702 = -t147 * t464 - t61;
t717 = -t53 * t147 + t54 * t663;
t662 = cos(qJ(1));
t549 = t662 * t424;
t425 = sin(qJ(1));
t600 = t425 * t428;
t338 = t419 * t549 + t600;
t417 = qJ(3) + qJ(4);
t411 = qJ(5) + t417;
t402 = sin(t411);
t403 = cos(t411);
t550 = t418 * t662;
t243 = t338 * t403 - t402 * t550;
t548 = t662 * t428;
t601 = t424 * t425;
t340 = -t419 * t601 + t548;
t606 = t418 * t425;
t247 = t340 * t403 + t402 * t606;
t607 = t418 * t424;
t292 = t402 * t419 + t403 * t607;
t465 = -g(1) * t247 - g(2) * t243 - g(3) * t292;
t701 = t157 * t147 - t465 - t9;
t246 = t340 * t402 - t403 * t606;
t291 = -t402 * t607 + t403 * t419;
t520 = -t338 * t402 - t403 * t550;
t466 = g(1) * t246 - g(2) * t520 - g(3) * t291;
t528 = -t660 * t19 + t421 * t23;
t10 = -qJD(5) * t54 - t528;
t8 = -pkin(5) * t448 - t10;
t459 = t466 - t8;
t713 = pkin(12) * t545 - t640;
t710 = pkin(5) * t596 + pkin(12) * t597 + t595;
t690 = t719 * t663;
t409 = sin(t417);
t410 = cos(t417);
t255 = -t340 * t409 + t410 * t606;
t707 = -t409 * t607 + t410 * t419;
t394 = pkin(8) * t607;
t658 = pkin(1) * t428;
t343 = t419 * t658 - t394;
t321 = qJD(2) * t343;
t671 = t663 * t21 - t459 * t420 + t51 * t577;
t673 = t663 * t486 + t51 * t578;
t668 = -t157 * t663 + t466 - t528;
t337 = -t419 * t548 + t601;
t700 = t243 * t420 - t337 * t426;
t699 = t243 * t426 + t337 * t420;
t664 = -t663 * t362 + t524;
t413 = t418 ^ 2;
t697 = 0.2e1 * t413;
t656 = pkin(4) * t474;
t2 = -t486 * qJD(6) + t420 * t16 + t426 * t7;
t1 = t2 * t426;
t693 = -t3 * t420 + t1;
t692 = -t21 * t719 - t3;
t169 = t421 * t228 + t229 * t660;
t639 = -qJD(5) * t169 + t727 * t421 + t660 * t728;
t563 = t661 * pkin(3);
t406 = t563 + pkin(4);
t602 = t421 * t422;
t240 = t406 * t537 + (-t422 * t581 + (t660 * t661 - t602) * qJD(4)) * pkin(3);
t117 = -t422 * t175 - t547;
t452 = t117 - t694;
t118 = t661 * t175 - t173;
t98 = t118 - t695;
t64 = t421 * t452 + t660 * t98;
t631 = t240 - t64;
t546 = t660 * t422;
t630 = -t421 * t98 + t452 * t660 + t406 * t581 + (t422 * t537 + (t421 * t661 + t546) * qJD(4)) * pkin(3);
t55 = t421 * t95 + t559;
t691 = -pkin(4) * t581 + t55;
t339 = t419 * t600 + t549;
t604 = t418 * t428;
t457 = -g(1) * t339 - g(2) * t337 + g(3) * t604;
t688 = t457 * t402;
t622 = t475 * t474;
t659 = pkin(1) * t424;
t344 = pkin(8) * t604 + t419 * t659;
t308 = pkin(9) * t419 + t344;
t220 = -t308 * t423 + t427 * t309;
t605 = t418 * t427;
t334 = t419 * t423 + t424 * t605;
t184 = -pkin(3) * t604 - pkin(10) * t334 + t220;
t221 = t427 * t308 + t423 * t309;
t333 = -t419 * t427 + t423 * t607;
t188 = -pkin(10) * t333 + t221;
t126 = t661 * t184 - t188 * t422;
t231 = -t422 * t333 + t334 * t661;
t104 = -pkin(4) * t604 - pkin(11) * t231 + t126;
t127 = t422 * t184 + t661 * t188;
t230 = t333 * t661 + t334 * t422;
t106 = -pkin(11) * t230 + t127;
t687 = t421 * t104 + t660 * t106;
t685 = -qJD(4) * t282 + t731;
t450 = -t219 - t457;
t504 = t427 * t544;
t679 = t504 - t583;
t678 = -t314 * t427 - t376 * t584;
t322 = qJD(2) * t344;
t56 = t660 * t95 - t634;
t85 = t100 + t656;
t28 = t420 * t85 + t426 * t56;
t655 = pkin(4) * t421;
t404 = pkin(12) + t655;
t510 = pkin(4) * t537;
t676 = -t404 * t578 + t426 * t510 - t28;
t263 = -t340 * t423 + t425 * t605;
t463 = t338 * t423 + t427 * t550;
t675 = -g(1) * t263 + g(2) * t463 + g(3) * t333;
t674 = t474 ^ 2 - t475 ^ 2;
t670 = t362 * t475 - t115;
t614 = t338 * t409;
t253 = -t410 * t550 - t614;
t667 = -g(1) * t255 - g(2) * t253 - g(3) * t707 - t214 * t474 + t36;
t254 = -t338 * t410 + t409 * t550;
t256 = t340 * t410 + t409 * t606;
t666 = -t214 * t475 + g(1) * t256 - g(2) * t254 - g(3) * (-t409 * t419 - t410 * t607) + t514;
t665 = t293 * t582 - t295 * t538 - t362 * t474 - t706;
t585 = qJD(2) * t428;
t542 = t418 * t585;
t260 = qJD(3) * t334 + t423 * t542;
t261 = -qJD(3) * t333 + t427 * t542;
t152 = t422 * t260 - t261 * t661 + t333 * t538 + t334 * t582;
t586 = qJD(2) * t424;
t543 = t418 * t586;
t319 = t418 * t477;
t165 = -qJD(3) * t221 + t427 * t319 - t321 * t423;
t132 = pkin(3) * t543 - pkin(10) * t261 + t165;
t164 = -t308 * t584 + t309 * t583 + t423 * t319 + t427 * t321;
t139 = -pkin(10) * t260 + t164;
t68 = -qJD(4) * t127 + t661 * t132 - t422 * t139;
t44 = pkin(4) * t543 + t152 * pkin(11) + t68;
t153 = qJD(4) * t231 + t260 * t661 + t422 * t261;
t67 = t422 * t132 + t661 * t139 + t184 * t538 - t188 * t582;
t48 = -pkin(11) * t153 + t67;
t14 = -qJD(5) * t687 - t421 * t48 + t44 * t660;
t430 = qJD(1) ^ 2;
t657 = pkin(3) * t295;
t646 = g(3) * t418;
t643 = t427 * pkin(3);
t251 = t347 * t421 + t449;
t252 = t347 * t660 - t458;
t407 = pkin(2) + t643;
t304 = pkin(4) * t476 - t407;
t167 = t251 * pkin(5) - t252 * pkin(12) + t304;
t107 = t167 * t426 - t169 * t420;
t642 = qJD(6) * t107 + t710 * t420 - t426 * t713;
t108 = t167 * t420 + t169 * t426;
t641 = -qJD(6) * t108 + t420 * t713 + t710 * t426;
t638 = pkin(5) * t545 - t639;
t637 = t486 * t719;
t636 = t486 * t420;
t628 = t135 * t133;
t621 = t252 * t420;
t620 = t252 * t426;
t619 = t295 * t293;
t618 = t295 * t376;
t612 = t376 * t423;
t611 = t403 * t420;
t610 = t403 * t426;
t608 = t413 * t430;
t603 = t420 * t428;
t599 = t426 * t428;
t365 = pkin(4) * t410 + t643;
t354 = pkin(2) + t365;
t414 = -pkin(11) + t429;
t594 = -t337 * t354 - t338 * t414;
t593 = -t339 * t354 - t340 * t414;
t364 = pkin(3) * t423 + pkin(4) * t409;
t592 = -t340 * t364 + t365 * t606;
t591 = -t364 * t607 + t419 * t365;
t332 = pkin(3) * t546 + t421 * t406;
t590 = t662 * pkin(1) + pkin(8) * t606;
t415 = t424 ^ 2;
t416 = t428 ^ 2;
t589 = t415 - t416;
t579 = qJD(6) * t719;
t576 = t293 * qJD(3);
t575 = t295 * qJD(3);
t564 = g(1) * t662;
t562 = t660 * pkin(4);
t555 = t428 * t608;
t554 = t423 * t606;
t553 = t418 * t599;
t382 = t418 * t603;
t539 = t418 * t419 * t430;
t536 = pkin(1) * t697;
t532 = -t425 * pkin(1) + pkin(8) * t550;
t238 = t520 * pkin(5);
t531 = t243 * pkin(12) + t238;
t239 = t246 * pkin(5);
t530 = pkin(12) * t247 - t239;
t287 = t291 * pkin(5);
t529 = pkin(12) * t292 + t287;
t527 = qJD(5) + t677;
t179 = t657 + t656;
t81 = t100 + t179;
t26 = t420 * t81 + t426 * t64;
t526 = t240 * t426 - t26;
t180 = t194 * t420 - t426 * t545;
t523 = t155 * t420 + t180;
t181 = t194 * t426 + t420 * t545;
t522 = t155 * t426 + t181;
t383 = t423 * t550;
t519 = t338 * t427 - t383;
t512 = qJD(2) + 0.2e1 * t587;
t511 = t393 + t571;
t508 = t424 * t555;
t503 = t424 * t534;
t498 = t255 * pkin(4);
t495 = pkin(5) * t403 + pkin(12) * t402;
t494 = g(1) * t520 + g(2) * t246;
t493 = g(1) * t337 - g(2) * t339;
t492 = g(1) * t340 + g(2) * t338;
t490 = qJDD(4) + qJDD(5) - t568;
t489 = -t404 * t60 + t712;
t324 = pkin(12) + t332;
t488 = -t324 * t60 + t712;
t487 = t21 * t420 - t426 * t486;
t70 = -pkin(12) * t604 + t687;
t170 = t230 * t660 + t231 * t421;
t171 = -t421 * t230 + t231 * t660;
t307 = t394 + (-pkin(2) - t658) * t419;
t237 = pkin(3) * t333 + t307;
t177 = pkin(4) * t230 + t237;
t86 = pkin(5) * t170 - pkin(12) * t171 + t177;
t34 = t420 * t86 + t426 * t70;
t33 = -t420 * t70 + t426 * t86;
t484 = -t339 * t414 + t340 * t354 + t364 * t606 + t590;
t483 = t707 * pkin(4);
t480 = -t338 * t364 - t365 * t550;
t479 = g(2) * t425 + t564;
t158 = t171 * t420 + t553;
t71 = t104 * t660 - t421 * t106;
t470 = t337 * t414 - t338 * t354 + t364 * t550 + t532;
t13 = t104 * t537 - t106 * t581 + t421 * t44 + t660 * t48;
t331 = -pkin(3) * t602 + t406 * t660;
t462 = t252 * t577 - t523;
t461 = -t252 * t578 - t522;
t460 = t1 + t465;
t456 = -g(3) * t607 - t492;
t455 = t469 * t427;
t453 = t253 * pkin(4);
t122 = -qJD(3) * t198 + t521;
t215 = pkin(3) * t260 + t322;
t444 = -qJD(6) * t487 + t693;
t443 = -t122 * t423 - t427 * t467 + t456;
t123 = pkin(4) * t153 + t215;
t405 = -t562 - pkin(5);
t323 = -pkin(5) - t331;
t316 = t354 * t604;
t264 = t340 * t427 + t554;
t200 = t247 * t426 + t339 * t420;
t199 = -t247 * t420 + t339 * t426;
t168 = -t228 * t660 + t229 * t421;
t159 = t171 * t426 - t382;
t80 = qJD(5) * t171 - t421 * t152 + t153 * t660;
t79 = t152 * t660 + t421 * t153 + t230 * t537 + t231 * t581;
t69 = pkin(5) * t604 - t71;
t66 = -qJD(6) * t382 + t171 * t577 - t420 * t79 - t426 * t543;
t65 = qJD(6) * t158 - t420 * t543 + t426 * t79;
t32 = t100 * t420 + t426 * t53;
t31 = t100 * t426 - t420 * t53;
t27 = -t420 * t56 + t426 * t85;
t25 = -t420 * t64 + t426 * t81;
t24 = pkin(5) * t80 + pkin(12) * t79 + t123;
t12 = -pkin(5) * t543 - t14;
t11 = pkin(12) * t543 + t13;
t5 = -qJD(6) * t34 - t11 * t420 + t24 * t426;
t4 = qJD(6) * t33 + t11 * t426 + t24 * t420;
t6 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t425 - g(2) * t662, t479, 0, 0 (qJDD(1) * t415 + 0.2e1 * t503) * t413 (t424 * t569 - t572 * t589) * t697 (t424 * t511 + t512 * t585) * t418 (qJDD(1) * t416 - 0.2e1 * t503) * t413 (t428 * t511 - t512 * t586) * t418, t393 * t419, -t322 * t513 + t343 * t393 - t506 * t419 + g(1) * t338 - g(2) * t340 + (-t535 + t569) * t536, -t234 * t419 - t321 * t513 - t344 * t393 - t468 * t536 - t493 ((-t317 * qJD(2) + qJDD(1) * t344 + t234) * t428 + (-qJD(2) * t320 - qJDD(1) * t343 + t506) * t424 - t479) * t418, t413 * qJDD(1) * pkin(1) ^ 2 - g(1) * t532 - g(2) * t590 + t234 * t344 - t317 * t322 + t320 * t321 - t343 * t506, t295 * t261 + t334 * t433, -t261 * t293 - t334 * t507 - t518 * t333 - t295 * t260 + (-t333 * t703 + t334 * t393) * t427, -t261 * t376 + t334 * t314 + (t295 * t586 - t428 * t433) * t418, t293 * t260 - t333 * t469, t260 * t376 - t333 * t314 + (-t293 * t586 - t428 * t469) * t418 (-t314 * t428 - t376 * t586) * t418, -t165 * t376 + t220 * t314 + t322 * t293 - t307 * t469 + t219 * t333 + t273 * t260 + g(1) * t519 - g(2) * t264 + (-t122 * t428 + t197 * t586) * t418, -g(1) * t463 - g(2) * t263 + t164 * t376 - t198 * t543 + t219 * t334 - t221 * t314 + t273 * t261 + t322 * t295 + t307 * t433 - t467 * t604, -t122 * t334 - t164 * t293 - t165 * t295 - t197 * t261 - t198 * t260 - t220 * t433 + t221 * t469 + t333 * t467 + t493, -t467 * t221 + t198 * t164 + t122 * t220 + t197 * t165 + t219 * t307 + t273 * t322 - g(1) * (-pkin(2) * t338 - pkin(9) * t337 + t532) - g(2) * (pkin(2) * t340 + pkin(9) * t339 + t590) -t115 * t231 - t152 * t474, t115 * t230 - t116 * t231 - t152 * t475 - t153 * t474, t152 * t362 + t231 * t302 + (t115 * t428 + t474 * t586) * t418, t116 * t230 - t153 * t475, t153 * t362 - t230 * t302 + (t116 * t428 + t475 * t586) * t418 (-t302 * t428 - t362 * t586) * t418, -t68 * t362 + t126 * t302 - t215 * t475 + t237 * t116 + t154 * t230 + t214 * t153 - g(1) * t254 - g(2) * t256 + (t110 * t586 - t36 * t428) * t418, -g(1) * t614 - g(2) * t255 - t237 * t115 - t127 * t302 - t214 * t152 + t154 * t231 + t215 * t474 + t67 * t362 + (-t111 * t586 - t410 * t564 - t428 * t514) * t418, t110 * t152 - t111 * t153 + t115 * t126 - t116 * t127 + t230 * t514 - t231 * t36 - t474 * t68 + t475 * t67 + t493, -t514 * t127 + t111 * t67 + t36 * t126 + t110 * t68 + t154 * t237 + t214 * t215 - g(1) * (pkin(3) * t383 + t337 * t429 - t338 * t407 + t532) - g(2) * (pkin(3) * t554 - t339 * t429 + t340 * t407 + t590) -t171 * t61 - t663 * t79, t147 * t79 + t170 * t61 - t171 * t62 - t663 * t80, -t79 * t527 + t171 * t490 + (t663 * t586 + t61 * t428 + (t171 * t586 + t428 * t79) * qJD(1)) * t418, t147 * t80 + t170 * t62, -t80 * t527 - t170 * t490 + (-t147 * t586 + t62 * t428 + (-t170 * t586 + t428 * t80) * qJD(1)) * t418 (-t490 * t428 + (t527 - 0.2e1 * t544) * t586) * t418, t14 * t527 + t71 * t490 + t123 * t147 + t177 * t62 + t82 * t170 + t157 * t80 + g(1) * t243 - g(2) * t247 + (t53 * t586 - t10 * t428 + (-t14 * t428 + t586 * t71) * qJD(1)) * t418, -t13 * t527 - t687 * t490 + t123 * t663 - t177 * t61 + t82 * t171 - t157 * t79 + (-t54 * t586 + t9 * t428 + (t13 * t428 - t586 * t687) * qJD(1)) * t418 + t494, -t10 * t171 - t13 * t147 - t14 * t663 - t170 * t9 + t53 * t79 - t54 * t80 + t61 * t71 - t62 * t687 + t493, -g(1) * t470 - g(2) * t484 + t10 * t71 + t157 * t123 + t54 * t13 + t53 * t14 + t82 * t177 + t687 * t9, -t135 * t65 - t159 * t41, t133 * t65 - t135 * t66 + t158 * t41 - t159 * t42, t135 * t80 + t159 * t60 - t170 * t41 - t65 * t719, t133 * t66 + t158 * t42, -t133 * t80 - t158 * t60 - t170 * t42 - t66 * t719, t170 * t60 + t719 * t80, g(1) * t699 - g(2) * t200 + t12 * t133 + t8 * t158 + t3 * t170 + t33 * t60 + t69 * t42 - t486 * t80 + t5 * t719 + t51 * t66, -g(1) * t700 - g(2) * t199 + t12 * t135 + t8 * t159 - t2 * t170 - t21 * t80 - t34 * t60 - t4 * t719 - t69 * t41 - t51 * t65, -t133 * t4 - t135 * t5 - t158 * t2 - t159 * t3 - t21 * t66 + t33 * t41 - t34 * t42 - t486 * t65 - t494, t2 * t34 + t21 * t4 + t3 * t33 - t486 * t5 + t8 * t69 + t51 * t12 - g(1) * (-pkin(5) * t243 + pkin(12) * t520 + t470) - g(2) * (pkin(5) * t247 + pkin(12) * t246 + t484); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t508, t589 * t608, -t428 * t539 + t533, t508, t424 * t539 + t392, t393, t320 * t513 + t608 * t659 - t457 - t506, pkin(1) * t555 + t317 * t513 + (pkin(8) * t572 + g(3)) * t607 + t492 + t552, 0, 0, t518 * t423 + (t423 * t442 - t618) * t427, t293 * t504 + t295 * t505 + (t433 - t576) * t427 + (t469 - t575) * t423, -t376 * t583 + t314 * t423 + (-t295 * t424 + t376 * t598) * t588, -t293 * t612 + t455 (t293 * t424 - t428 * t612) * t588 - t678, t376 * t545, -pkin(2) * t507 + t224 * t376 - t197 * t545 - t320 * t293 + (-pkin(9) * t314 - t273 * t376) * t423 + (pkin(9) * qJD(3) * t376 + t450 + t644) * t427, -pkin(2) * t433 + pkin(9) * t678 + t198 * t545 - t225 * t376 - t273 * t679 - t320 * t295 - t423 * t450, t224 * t295 + t225 * t293 + t443 + t730 * t198 + t679 * t197 + ((t433 + t576) * t423 + t427 * t575 + t455) * pkin(9), -t197 * t224 - t198 * t225 - t273 * t320 + t450 * pkin(2) + ((-t197 * t427 - t198 * t423) * qJD(3) + t443) * pkin(9), -t115 * t347 + t474 * t680, t115 * t476 - t347 * t116 + t474 * t682 + t475 * t680, t347 * t302 - t362 * t680 - t474 * t545, t116 * t476 + t475 * t682, -t476 * t302 - t362 * t682 - t475 * t545, t362 * t545, -t110 * t545 - t407 * t116 + t154 * t476 - t214 * t682 + t281 * t302 - t362 * t685 - t410 * t457 + t475 * t684, t111 * t545 + t407 * t115 + t154 * t347 + t214 * t680 - t282 * t302 + t362 * t686 + t409 * t457 - t474 * t684, -t110 * t680 + t111 * t682 + t281 * t115 - t282 * t116 - t36 * t347 - t474 * t685 + t475 * t686 + t476 * t514 + t456, -t514 * t282 + t36 * t281 - t154 * t407 - g(1) * (-t339 * t407 - t340 * t429) - g(2) * (-t337 * t407 - t338 * t429) - (t407 * t428 - t424 * t429) * t646 - t684 * t214 + t686 * t111 + t685 * t110, -t252 * t61 - t597 * t663, t147 * t597 + t251 * t61 - t252 * t62 - t596 * t663, t252 * t448 + t464 * t597 - t545 * t663, t147 * t596 + t251 * t62, t147 * t545 - t251 * t448 + t464 * t596, t464 * t545, t147 * t595 + t157 * t596 - t168 * t448 + t82 * t251 + t304 * t62 - t403 * t457 - t464 * t639 - t53 * t545, -t157 * t597 - t169 * t448 + t82 * t252 - t304 * t61 + t464 * t640 + t54 * t545 + t595 * t663 + t688, -t10 * t252 - t147 * t640 - t168 * t61 - t169 * t62 - t251 * t9 + t53 * t597 - t54 * t596 - t639 * t663 + t456, t9 * t169 - t10 * t168 + t82 * t304 - g(1) * t593 - g(2) * t594 - g(3) * (-t414 * t607 + t316) + t640 * t54 + t639 * t53 + t595 * t157, t135 * t461 - t41 * t620, t523 * t135 + t522 * t133 + (t38 - t40 + (-t135 * t426 + t629) * qJD(6)) * t252, t135 * t596 - t251 * t41 + t461 * t719 + t60 * t620, t133 * t462 + t42 * t621, -t133 * t596 - t251 * t42 - t462 * t719 - t60 * t621, t251 * t60 + t596 * t719, t107 * t60 + t3 * t251 + t168 * t42 + t8 * t621 - g(1) * (-t339 * t610 + t340 * t420) - g(2) * (-t337 * t610 + t338 * t420) - (t403 * t599 + t420 * t424) * t646 - t596 * t486 + t641 * t719 + t638 * t133 + t462 * t51, -t108 * t60 - t2 * t251 - t168 * t41 + t8 * t620 - g(1) * (t339 * t611 + t340 * t426) - g(2) * (t337 * t611 + t338 * t426) - (-t403 * t603 + t424 * t426) * t646 - t596 * t21 - t642 * t719 + t638 * t135 + t461 * t51, t107 * t41 - t108 * t42 + t180 * t21 - t181 * t486 + t487 * t155 - t641 * t135 - t642 * t133 - t688 + (-t2 * t420 - t3 * t426 + (-t21 * t426 - t636) * qJD(6)) * t252, t2 * t108 + t3 * t107 + t8 * t168 - g(1) * (-t339 * t495 + t593) - g(2) * (-t337 * t495 + t594) - g(3) * t316 + t638 * t51 - (-t414 * t424 + t428 * t495) * t646 + t642 * t21 - t641 * t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t619, -t293 ^ 2 + t295 ^ 2, -t293 * t376 + t433, -t619, t469 - t618, t314, -t198 * t376 - t273 * t295 + t122 + t675, g(1) * t264 + g(2) * t519 + g(3) * t334 - t197 * t376 + t273 * t293 + t467, 0, 0, -t622, t674, t670, t622, t665, t302, t117 * t362 + (t295 * t475 + t302 * t661 + t362 * t582) * pkin(3) + t667, -t118 * t362 + (-t295 * t474 - t422 * t302 + t362 * t538) * pkin(3) + t666, t110 * t475 + t111 * t474 + t117 * t474 - t118 * t475 + (t115 * t661 - t116 * t422 + (t422 * t474 + t475 * t661) * qJD(4)) * pkin(3), -t110 * t117 - t111 * t118 - t214 * t657 + t36 * t563 + (-t110 * t582 + t111 * t538 - t422 * t514 + t675) * pkin(3), t623, t704, t702, -t623, t664, t448, t331 * t448 - t179 * t147 + (-t54 - t630) * qJD(5) + t630 * t362 + t668, -t179 * t663 - t332 * t448 + t464 * t631 + t701, -t147 * t631 + t331 * t61 - t332 * t62 + t630 * t663 + t717, -g(1) * t592 - g(2) * t480 - g(3) * t591 + t10 * t331 - t157 * t179 + t9 * t332 - t53 * t630 + t54 * t631, t715, t724, t714, t725, t723, -t690, -t719 * t25 + t323 * t42 + t630 * t133 + (-t240 * t719 + t488) * t420 + (-t324 * t579 + t459) * t426 + t673, -t323 * t41 + t488 * t426 + t630 * t135 + (t324 * t578 - t526) * t719 + t671, t133 * t26 + t135 * t25 + (-t133 * t240 + t721 - t324 * t42 + (t135 * t324 + t486) * qJD(6)) * t426 + (t135 * t240 - t324 * t41 + (t133 * t324 - t21) * qJD(6) + t718) * t420 + t460, t8 * t323 - g(1) * (t530 + t592) - g(2) * (t480 + t531) - g(3) * (t529 + t591) + t630 * t51 + t526 * t21 - (-t240 * t420 - t25) * t486 + t444 * t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t622, t674, t670, t622, t665, t302, -t111 * t362 + t667, -t110 * t362 + t666, 0, 0, t623, t704, t702, -t623, t664, t448, -t55 * t362 + (t55 - t54) * qJD(5) + (-t147 * t474 + t448 * t660 + t464 * t581) * pkin(4) + t668, -t56 * t464 + (-t421 * t448 + t464 * t537 - t474 * t663) * pkin(4) + t701, t56 * t147 - t55 * t663 + (t660 * t61 - t421 * t62 + (-t147 * t660 + t421 * t663) * qJD(5)) * pkin(4) + t717, -g(1) * t498 - g(2) * t453 - g(3) * t483 + t10 * t562 - t157 * t656 + t655 * t9 + (t510 - t56) * t54 + t691 * t53, t715, t724, t714, t725, t723, -t690, -t27 * t719 + t405 * t42 - t691 * t133 + (-t510 * t719 + t489) * t420 + (-t404 * t579 + t459) * t426 + t673, -t135 * t691 - t405 * t41 + t489 * t426 - t676 * t719 + t671, t28 * t133 + t27 * t135 + (-t133 * t510 + t721 - t404 * t42 + (t135 * t404 + t486) * qJD(6)) * t426 + (t135 * t510 - t404 * t41 + (t133 * t404 - t21) * qJD(6) + t718) * t420 + t460, t510 * t636 + t8 * t405 + t486 * t27 - g(1) * (t498 + t530) - g(2) * (t453 + t531) - g(3) * (t483 + t529) - t691 * t51 + (t486 * t577 + t693) * t404 + t676 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623, t704, t702, -t623, t664, t448, -t362 * t54 + t668, -t464 * t53 + t701, 0, 0, t715, t724, t714, t133 * t726 - t40, t723, -t690, -pkin(5) * t42 - t133 * t54 - t719 * t31 + (-pkin(12) * t60 + t712) * t420 + (-pkin(12) * t579 + t459) * t426 + t673, t51 * t722 + pkin(5) * t41 - t135 * t54 + t719 * t32 + (t578 * t719 - t58) * pkin(12) + t671, t133 * t32 + t135 * t31 + (t637 + (-t42 + t580) * pkin(12)) * t426 + ((qJD(6) * t133 - t41) * pkin(12) + t692) * t420 + t460, -t8 * pkin(5) + g(1) * t239 - g(2) * t238 - g(3) * t287 + t486 * t31 - t21 * t32 - t51 * t54 + (t444 + t465) * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t628, -t133 ^ 2 + t135 ^ 2, t133 * t719 - t41, -t628, t135 * t719 - t42, t60, -t51 * t135 - g(1) * t199 + g(2) * t700 - g(3) * (-t292 * t420 - t553) - t692, -t637 + t51 * t133 + g(1) * t200 + g(2) * t699 - g(3) * (-t292 * t426 + t382) - t2, 0, 0;];
tau_reg  = t6;