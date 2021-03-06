% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRP8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:58
% EndTime: 2019-03-10 01:56:47
% DurationCPUTime: 26.44s
% Computational Cost: add. (35215->985), mult. (85192->1232), div. (0->0), fcn. (68465->14), ass. (0->437)
t409 = cos(qJ(2));
t406 = sin(qJ(2));
t401 = sin(pkin(6));
t580 = qJD(1) * t401;
t544 = t406 * t580;
t402 = cos(pkin(6));
t579 = qJD(1) * t402;
t561 = pkin(1) * t579;
t315 = -pkin(8) * t544 + t409 * t561;
t499 = pkin(2) * t406 - pkin(9) * t409;
t316 = t499 * t580;
t405 = sin(qJ(3));
t660 = cos(qJ(3));
t219 = -t315 * t405 + t316 * t660;
t696 = pkin(10) + pkin(9);
t367 = t696 * t660;
t548 = t409 * t660;
t724 = (pkin(3) * t406 - pkin(10) * t548) * t580 + t219 + qJD(3) * t367;
t220 = t315 * t660 + t316 * t405;
t366 = t696 * t405;
t543 = t409 * t580;
t511 = t405 * t543;
t723 = pkin(10) * t511 - qJD(3) * t366 - t220;
t404 = sin(qJ(4));
t659 = cos(qJ(4));
t451 = -t404 * t405 + t659 * t660;
t566 = qJD(3) + qJD(4);
t250 = t451 * t566;
t263 = t451 * t543;
t587 = -t250 + t263;
t335 = t404 * t660 + t405 * t659;
t586 = (-t543 + t566) * t335;
t408 = cos(qJ(5));
t572 = qJD(5) * t408;
t550 = t401 * t660;
t506 = qJD(1) * t550;
t520 = qJD(2) + t579;
t291 = t405 * t520 + t406 * t506;
t290 = -t405 * t544 + t520 * t660;
t687 = t659 * t290;
t202 = t291 * t404 - t687;
t708 = t202 * t408;
t722 = t572 + t708;
t600 = t401 * t406;
t705 = t402 * t660 - t405 * t600;
t472 = t705 * pkin(3);
t534 = qJD(4) * t659;
t576 = qJD(4) * t404;
t590 = -t366 * t534 - t367 * t576 - t404 * t724 + t659 * t723;
t598 = t401 * t409;
t384 = pkin(8) * t598;
t657 = pkin(1) * t406;
t583 = t402 * t657 + t384;
t318 = t583 * qJD(1);
t577 = qJD(3) * t405;
t501 = -t318 + (-t511 + t577) * pkin(3);
t403 = sin(qJ(5));
t362 = -qJD(3) + t543;
t466 = -qJD(4) + t362;
t449 = t408 * t466;
t464 = t290 * t404 + t291 * t659;
t166 = t403 * t464 + t449;
t168 = -t403 * t466 + t408 * t464;
t261 = qJD(2) * pkin(9) + (t384 + (pkin(9) + t657) * t402) * qJD(1);
t471 = -pkin(2) * t409 - pkin(9) * t406 - pkin(1);
t303 = t471 * t401;
t275 = qJD(1) * t303;
t183 = -t261 * t405 + t275 * t660;
t152 = -pkin(10) * t291 + t183;
t141 = -pkin(3) * t362 + t152;
t184 = t261 * t660 + t275 * t405;
t153 = pkin(10) * t290 + t184;
t148 = t404 * t153;
t95 = t141 * t659 - t148;
t91 = pkin(4) * t466 - t95;
t50 = pkin(5) * t166 - qJ(6) * t168 + t91;
t703 = qJD(5) + t202;
t721 = t50 * t703;
t720 = pkin(11) * t544 - t590;
t719 = pkin(4) * t586 + pkin(11) * t587 + t501;
t573 = qJD(5) * t403;
t474 = t409 * t506;
t569 = qJDD(1) * t402;
t509 = qJDD(2) + t569;
t513 = t406 * t550;
t717 = qJD(3) * t290;
t187 = -qJD(2) * t474 - qJDD(1) * t513 - t405 * t509 - t717;
t568 = qJDD(1) * t406;
t531 = t401 * t568;
t570 = qJD(1) * qJD(2);
t532 = t409 * t570;
t706 = t401 * t532 + t531;
t421 = t291 * qJD(3) + t405 * t706 - t660 * t509;
t414 = t187 * t659 + t291 * t576 + t404 * t421;
t413 = -t290 * t534 + t414;
t533 = t406 * t570;
t504 = t401 * t533;
t567 = qJDD(1) * t409;
t378 = t401 * t567;
t565 = t378 - qJDD(3);
t304 = t504 - t565;
t438 = qJDD(4) + t304;
t68 = qJD(5) * t449 - t403 * t438 + t408 * t413 + t464 * t573;
t69 = qJD(5) * t168 - t403 * t413 - t408 * t438;
t718 = t166 * t572 + t168 * t573 + t403 * t69 + t68 * t408 + (t166 * t408 + t168 * t403) * t202;
t489 = pkin(5) * t408 + qJ(6) * t403;
t260 = -pkin(2) * t520 - t315;
t206 = -pkin(3) * t290 + t260;
t111 = pkin(4) * t202 - pkin(11) * t464 + t206;
t149 = t659 * t153;
t96 = t141 * t404 + t149;
t92 = -pkin(11) * t466 + t96;
t44 = t111 * t408 - t403 * t92;
t591 = qJD(6) - t44;
t36 = -pkin(5) * t703 + t591;
t518 = qJD(2) * t561;
t558 = pkin(1) * t569;
t552 = pkin(8) * t378 + t406 * t558 + t409 * t518;
t226 = -pkin(8) * t504 + t552;
t210 = pkin(9) * t509 + t226;
t465 = t499 * qJD(2);
t214 = (qJD(1) * t465 + qJDD(1) * t471) * t401;
t108 = -qJD(3) * t184 - t405 * t210 + t214 * t660;
t78 = t304 * pkin(3) + t187 * pkin(10) + t108;
t535 = qJD(3) * t660;
t458 = -t210 * t660 - t214 * t405 + t261 * t577 - t275 * t535;
t85 = -pkin(10) * t421 - t458;
t21 = t141 * t534 - t153 * t576 + t404 * t78 + t659 * t85;
t19 = pkin(11) * t438 + t21;
t527 = t404 * t187 - t421 * t659;
t103 = qJD(4) * t464 - t527;
t514 = pkin(8) * t706 + t406 * t518 - t409 * t558;
t211 = -pkin(2) * t509 + t514;
t137 = pkin(3) * t421 + t211;
t31 = t103 * pkin(4) + pkin(11) * t413 + t137;
t6 = t111 * t572 + t19 * t408 + t31 * t403 - t573 * t92;
t100 = qJDD(5) + t103;
t628 = qJ(6) * t100;
t2 = qJD(6) * t703 + t6 + t628;
t652 = pkin(5) * t100;
t7 = -t111 * t573 - t19 * t403 + t31 * t408 - t572 * t92;
t4 = qJDD(6) - t7 - t652;
t661 = cos(qJ(1));
t547 = t661 * t406;
t407 = sin(qJ(1));
t593 = t407 * t409;
t326 = t402 * t547 + t593;
t400 = qJ(3) + qJ(4);
t394 = sin(t400);
t395 = cos(t400);
t551 = t401 * t661;
t244 = t326 * t395 - t394 * t551;
t546 = t661 * t409;
t594 = t406 * t407;
t328 = -t402 * t594 + t546;
t599 = t401 * t407;
t248 = t328 * t395 + t394 * t599;
t300 = t394 * t402 + t395 * t600;
t454 = g(1) * t248 + g(2) * t244 + g(3) * t300;
t676 = t2 * t408 + t4 * t403 - t454;
t716 = t36 * t722 + t676;
t625 = t166 * t403;
t65 = t69 * t408;
t715 = t625 * t703 - t65;
t709 = t202 * t403;
t714 = -t65 + (t573 + t709) * t166;
t63 = t68 * t403;
t713 = t168 * t722 - t63;
t97 = t403 * t100;
t712 = -t168 * t464 + t703 * t722 + t97;
t98 = t408 * t100;
t629 = -t573 * t703 + t98;
t442 = -t166 * t464 + t703 * t709 - t629;
t711 = t202 * t44;
t710 = t202 * t50;
t636 = t202 * t91;
t615 = t202 * t464;
t105 = t152 * t404 + t149;
t686 = -pkin(3) * t576 + t105;
t266 = -t366 * t404 + t367 * t659;
t588 = qJD(4) * t266 + t404 * t723 + t659 * t724;
t256 = -t328 * t405 + t407 * t550;
t704 = -pkin(4) * t395 - pkin(11) * t394;
t381 = pkin(8) * t600;
t656 = pkin(1) * t409;
t330 = t402 * t656 - t381;
t319 = qJD(2) * t330;
t702 = -t202 ^ 2 + t464 ^ 2;
t247 = t328 * t394 - t395 * t599;
t299 = -t394 * t600 + t395 * t402;
t524 = -t326 * t394 - t395 * t551;
t455 = g(1) * t247 - g(2) * t524 - g(3) * t299;
t138 = pkin(4) * t464 + pkin(11) * t202;
t701 = t202 * t206 - t21 + t454;
t488 = pkin(5) * t403 - qJ(6) * t408;
t700 = pkin(5) * t573 - qJ(6) * t572 - qJD(6) * t403 + t202 * t488;
t699 = -t202 * t466 - t413;
t325 = -t402 * t546 + t594;
t190 = t244 * t403 - t325 * t408;
t191 = t244 * t408 + t325 * t403;
t222 = t263 * t408 + t403 * t544;
t525 = -t250 * t408 + t222;
t221 = t263 * t403 - t408 * t544;
t526 = -t250 * t403 + t221;
t698 = t525 * t166 + t526 * t168 + t335 * (qJD(5) * (-t168 * t408 + t625) - t65 + t63);
t694 = pkin(5) * t464;
t447 = t6 * t408 - t454;
t392 = pkin(3) * t660 + pkin(2);
t236 = -pkin(4) * t451 - pkin(11) * t335 - t392;
t641 = t236 * t572 - t266 * t573 + t403 * t719 - t408 * t720;
t519 = pkin(3) * t534;
t482 = t408 * t519;
t654 = pkin(3) * t404;
t390 = pkin(11) + t654;
t604 = t390 * t408;
t693 = -t166 * t482 - t604 * t69;
t631 = t700 - t686;
t692 = -t168 * t703 + t69;
t691 = qJ(6) * t464;
t690 = t703 * t464;
t327 = t402 * t593 + t547;
t677 = g(1) * t327 + g(2) * t325;
t446 = g(3) * t598 - t677;
t432 = t446 * t394;
t685 = t184 * t362 - t108;
t302 = pkin(9) * t402 + t583;
t212 = -t302 * t405 + t303 * t660;
t323 = t402 * t405 + t513;
t164 = -pkin(3) * t598 - pkin(10) * t323 + t212;
t213 = t302 * t660 + t303 * t405;
t176 = pkin(10) * t705 + t213;
t116 = t164 * t404 + t176 * t659;
t563 = pkin(11) * t598;
t114 = -t563 + t116;
t223 = t323 * t404 - t659 * t705;
t224 = t323 * t659 + t404 * t705;
t301 = t381 + (-pkin(2) - t656) * t402;
t229 = t301 - t472;
t132 = pkin(4) * t223 - pkin(11) * t224 + t229;
t684 = t114 * t408 + t132 * t403;
t589 = pkin(4) * t544 + t588;
t683 = t489 * t524;
t682 = t489 * t247;
t681 = t236 * t403 + t266 * t408;
t680 = t489 * t299;
t540 = t390 * t572;
t679 = -t403 * t519 - t540;
t539 = t390 * t573;
t678 = t539 - t482;
t320 = t583 * qJD(2);
t106 = t152 * t659 - t148;
t655 = pkin(3) * t291;
t124 = t138 + t655;
t54 = t106 * t408 + t124 * t403;
t675 = -t54 - t678;
t674 = t36 * t464 + t50 * t573;
t673 = -t44 * t464 + t573 * t91;
t45 = t111 * t403 + t408 * t92;
t37 = qJ(6) * t703 + t45;
t553 = t455 * t403;
t523 = t141 * t576 + t153 * t534 + t404 * t85 - t659 * t78;
t20 = -pkin(4) * t438 + t523;
t8 = pkin(5) * t69 + qJ(6) * t68 - qJD(6) * t168 + t20;
t672 = -t37 * t464 - t8 * t403 + t553;
t671 = t20 * t403 + t45 * t464 + t572 * t91 - t553;
t668 = -t206 * t464 + t455 - t523;
t667 = -t362 * t464 + t527;
t317 = t401 * t465;
t143 = -qJD(3) * t213 + t317 * t660 - t319 * t405;
t541 = qJD(2) * t598;
t253 = qJD(3) * t705 + t541 * t660;
t578 = qJD(2) * t406;
t542 = t401 * t578;
t123 = pkin(3) * t542 - pkin(10) * t253 + t143;
t142 = -t302 * t577 + t303 * t535 + t317 * t405 + t319 * t660;
t252 = qJD(3) * t323 + t405 * t541;
t130 = -pkin(10) * t252 + t142;
t38 = t123 * t404 + t130 * t659 + t164 * t534 - t176 * t576;
t33 = pkin(11) * t542 + t38;
t135 = qJD(4) * t223 + t252 * t404 - t253 * t659;
t136 = qJD(4) * t224 + t252 * t659 + t253 * t404;
t207 = pkin(3) * t252 + t320;
t71 = pkin(4) * t136 + pkin(11) * t135 + t207;
t12 = -qJD(5) * t684 - t33 * t403 + t408 * t71;
t640 = -qJD(5) * t681 + t403 * t720 + t408 * t719;
t453 = t335 * t572 - t526;
t664 = t166 * t586 + t335 * t97 - t451 * t69 + t453 * t703;
t663 = t168 ^ 2;
t411 = qJD(1) ^ 2;
t662 = qJD(2) ^ 2;
t397 = t401 ^ 2;
t658 = pkin(1) * t397;
t651 = pkin(11) * t100;
t643 = qJ(6) * t586 - qJD(6) * t451 + t641;
t642 = -pkin(5) * t586 - t640;
t639 = -t221 * pkin(5) + t222 * qJ(6) + t488 * t250 + (qJD(5) * t489 - qJD(6) * t408) * t335 + t589;
t638 = t703 * t37;
t637 = t703 * t45;
t630 = t700 - t96;
t56 = t138 * t403 + t408 * t95;
t627 = t100 * t390;
t626 = t166 * t390;
t624 = t168 * t166;
t614 = t290 * t362;
t613 = t291 * t290;
t612 = t291 * t362;
t609 = t326 * t405;
t607 = t335 * t403;
t606 = t335 * t408;
t605 = t390 * t403;
t603 = t395 * t403;
t602 = t395 * t408;
t601 = t397 * t411;
t595 = t405 * t409;
t592 = t408 * t409;
t585 = -t325 * t392 + t326 * t696;
t584 = -t327 * t392 + t328 * t696;
t582 = pkin(1) * t661 + pkin(8) * t599;
t398 = t406 ^ 2;
t399 = t409 ^ 2;
t581 = t398 - t399;
t574 = qJD(5) * t703;
t383 = pkin(4) * t598;
t562 = t659 * pkin(3);
t557 = t409 * t601;
t555 = t401 * t592;
t370 = t403 * t598;
t554 = t405 * t599;
t545 = t660 * t304;
t537 = t401 * t402 * t411;
t536 = t166 ^ 2 - t663;
t530 = -pkin(1) * t407 + pkin(8) * t551;
t528 = qJDD(4) - t565;
t265 = t366 * t659 + t367 * t404;
t517 = t406 * t557;
t516 = t325 * t704 + t585;
t515 = t327 * t704 + t584;
t371 = t405 * t551;
t508 = t661 * t660;
t505 = t406 * t532;
t502 = t256 * pkin(3);
t194 = t248 * t403 - t327 * t408;
t498 = -g(1) * t190 + g(2) * t194;
t195 = t248 * t408 + t327 * t403;
t497 = g(1) * t191 - g(2) * t195;
t496 = g(1) * t524 + g(2) * t247;
t495 = g(1) * t325 - g(2) * t327;
t494 = g(1) * t328 + g(2) * t326;
t493 = t392 * t598 + t600 * t696;
t115 = t164 * t659 - t176 * t404;
t492 = t326 * t660 - t371;
t490 = (qJD(5) * t166 - t68) * pkin(11);
t486 = t36 * t408 - t37 * t403;
t485 = t403 * t45 + t408 * t44;
t481 = pkin(3) * t554 + t327 * t696 + t328 * t392 + t582;
t102 = -qJD(5) * t370 - t135 * t403 + t224 * t572 - t408 * t542;
t197 = t224 * t403 + t555;
t480 = t102 * t166 + t197 * t69;
t55 = t138 * t408 - t403 * t95;
t53 = -t106 * t403 + t124 * t408;
t66 = -t114 * t403 + t132 * t408;
t173 = t236 * t408 - t266 * t403;
t353 = -pkin(4) - t489;
t468 = g(1) * t661 + g(2) * t407;
t113 = t383 - t115;
t463 = pkin(3) * t371 - t325 * t696 - t326 * t392 + t530;
t11 = -t114 * t573 + t132 * t572 + t33 * t408 + t403 * t71;
t460 = t383 * t395 + t394 * t563 + t493;
t39 = t123 * t659 - t130 * t404 - t164 * t576 - t176 * t534;
t459 = t532 + t568;
t215 = -t325 * t603 - t326 * t408;
t217 = -t327 * t603 - t328 * t408;
t276 = t370 * t395 - t408 * t600;
t457 = g(1) * t217 + g(2) * t215 + g(3) * t276;
t216 = -t325 * t602 + t326 * t403;
t218 = -t327 * t602 + t328 * t403;
t277 = (t395 * t592 + t403 * t406) * t401;
t456 = -g(1) * t218 - g(2) * t216 - g(3) * t277;
t452 = t335 * t573 + t525;
t450 = -t304 * t405 + t362 * t535;
t239 = t247 * pkin(4);
t448 = pkin(11) * t248 - t239 + t502;
t445 = -g(3) * t600 - t494;
t444 = t168 * t519 - t390 * t68;
t441 = -t519 * t703 - t627;
t440 = pkin(4) * t248 + pkin(11) * t247 + t481;
t439 = -t401 * t508 - t609;
t292 = t299 * pkin(4);
t436 = pkin(11) * t300 + t292 + t472;
t435 = -0.2e1 * t533 + 0.2e1 * t567;
t434 = -t211 - t446;
t431 = -pkin(4) * t244 + pkin(11) * t524 + t463;
t430 = t439 * pkin(3);
t429 = t474 - t535;
t101 = qJD(5) * t197 + t135 * t408 - t403 * t542;
t198 = t224 * t408 - t370;
t428 = t101 * t166 - t102 * t168 + t197 * t68 - t198 * t69;
t427 = t100 * t197 + t102 * t703 + t136 * t166 + t223 * t69;
t426 = -pkin(11) * t574 + t455;
t241 = t300 * t403 + t555;
t425 = g(1) * t194 + g(2) * t190 + g(3) * t241 + t7;
t424 = -t390 * t574 + t455;
t34 = -pkin(4) * t542 - t39;
t420 = -t458 * t660 + t445;
t419 = t166 * t453 + t607 * t69;
t242 = t300 * t408 - t370;
t418 = -g(1) * t195 - g(2) * t191 - g(3) * t242 + t6;
t237 = t524 * pkin(4);
t417 = pkin(11) * t244 + t237 + t430;
t416 = t421 * t660;
t415 = t168 * t50 + qJDD(6) - t425;
t391 = -t562 - pkin(4);
t332 = -t562 + t353;
t257 = t328 * t660 + t554;
t180 = t335 * t488 + t265;
t151 = pkin(5) * t451 - t173;
t150 = -qJ(6) * t451 + t681;
t117 = pkin(5) * t168 + qJ(6) * t166;
t72 = pkin(5) * t197 - qJ(6) * t198 + t113;
t60 = pkin(11) * t65;
t52 = -pkin(5) * t223 - t66;
t51 = qJ(6) * t223 + t684;
t48 = -t55 - t694;
t47 = t56 + t691;
t46 = -t100 * t451 + t586 * t703;
t43 = -t53 - t694;
t42 = t54 + t691;
t41 = t100 * t223 + t136 * t703;
t40 = t166 * t703 - t68;
t24 = -t101 * t168 - t198 * t68;
t23 = -t168 * t452 - t606 * t68;
t15 = pkin(5) * t102 + qJ(6) * t101 - qJD(6) * t198 + t34;
t14 = t168 * t586 + t335 * t98 + t451 * t68 - t452 * t703;
t13 = t100 * t198 - t101 * t703 + t136 * t168 - t223 * t68;
t10 = -pkin(5) * t136 - t12;
t9 = qJ(6) * t136 + qJD(6) * t223 + t11;
t1 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t407 - g(2) * t661, t468, 0, 0 (qJDD(1) * t398 + 0.2e1 * t505) * t397, 0.2e1 * (t406 * t567 - t570 * t581) * t397 (qJDD(2) * t406 + 0.2e1 * t402 * t459 + t409 * t662) * t401 (qJDD(1) * t399 - 0.2e1 * t505) * t397 (qJDD(2) * t409 + t402 * t435 - t406 * t662) * t401, t509 * t402, g(1) * t326 - g(2) * t328 - t320 * t520 + t330 * t509 - t402 * t514 + t435 * t658, -t226 * t402 - t319 * t520 - 0.2e1 * t459 * t658 - t509 * t583 - t495 ((-t315 * qJD(2) + t583 * qJDD(1) + t226) * t409 + (-qJD(2) * t318 - t330 * qJDD(1) + t514) * t406 - t468) * t401, pkin(1) ^ 2 * qJDD(1) * t397 - g(1) * t530 - g(2) * t582 + t226 * t583 - t315 * t320 + t318 * t319 - t330 * t514, -t187 * t323 + t253 * t291, -t187 * t705 - t252 * t291 + t253 * t290 - t323 * t421, -t253 * t362 + t304 * t323 + (t187 * t409 + t291 * t578) * t401, -t252 * t290 - t421 * t705, t252 * t362 + t705 * t304 + (t290 * t578 + t409 * t421) * t401 (-t304 * t409 - t362 * t578) * t401, -t143 * t362 + t212 * t304 - t320 * t290 + t301 * t421 - t211 * t705 + t260 * t252 + g(1) * t492 - g(2) * t257 + (-t108 * t409 + t183 * t578) * t401, -g(1) * t609 - g(2) * t256 + t142 * t362 - t301 * t187 + t211 * t323 - t213 * t304 + t260 * t253 + t320 * t291 + (-g(1) * t508 - t184 * t578 - t409 * t458) * t401, -t108 * t323 + t142 * t290 - t143 * t291 - t183 * t253 - t184 * t252 + t187 * t212 - t213 * t421 - t458 * t705 + t495, -t458 * t213 + t184 * t142 + t108 * t212 + t183 * t143 + t211 * t301 + t260 * t320 - g(1) * (-pkin(2) * t326 - pkin(9) * t325 + t530) - g(2) * (pkin(2) * t328 + pkin(9) * t327 + t582) -t135 * t464 - t224 * t413, -t103 * t224 + t135 * t202 - t136 * t464 + t223 * t413, -t135 * t566 + t224 * t528 + (t413 * t409 + t464 * t578 + (t135 * t409 + t224 * t578) * qJD(1)) * t401, t103 * t223 + t136 * t202, -t136 * t566 - t223 * t528 + (-t202 * t578 + t103 * t409 + (t136 * t409 - t223 * t578) * qJD(1)) * t401 (-t528 * t409 + (-0.2e1 * t543 + t566) * t578) * t401, t39 * t566 + t115 * t528 + t207 * t202 + t229 * t103 + t137 * t223 + t206 * t136 + g(1) * t244 - g(2) * t248 + (t95 * t578 + t523 * t409 + (t115 * t578 - t39 * t409) * qJD(1)) * t401, -t38 * t566 - t116 * t528 + t207 * t464 - t229 * t413 + t137 * t224 - t206 * t135 + (-t96 * t578 + t21 * t409 + (-t116 * t578 + t38 * t409) * qJD(1)) * t401 + t496, -t103 * t116 + t115 * t413 + t135 * t95 - t136 * t96 - t202 * t38 - t21 * t223 + t224 * t523 - t39 * t464 + t495, -g(1) * t463 - g(2) * t481 - t115 * t523 + t116 * t21 + t137 * t229 + t206 * t207 + t38 * t96 + t39 * t95, t24, t428, t13, t480, -t427, t41, t100 * t66 + t102 * t91 + t113 * t69 + t12 * t703 + t136 * t44 + t166 * t34 + t197 * t20 + t223 * t7 + t497, -t100 * t684 - t101 * t91 - t11 * t703 - t113 * t68 - t136 * t45 + t168 * t34 + t198 * t20 - t223 * t6 + t498, t101 * t44 - t102 * t45 - t11 * t166 - t12 * t168 - t197 * t6 - t198 * t7 + t66 * t68 - t684 * t69 - t496, -g(1) * t431 - g(2) * t440 + t45 * t11 + t20 * t113 + t44 * t12 + t91 * t34 + t6 * t684 + t7 * t66, t24, t13, -t428, t41, t427, t480, -t10 * t703 - t100 * t52 + t102 * t50 - t136 * t36 + t15 * t166 + t197 * t8 - t223 * t4 + t69 * t72 + t497, t10 * t168 - t101 * t36 - t102 * t37 - t166 * t9 - t197 * t2 + t198 * t4 - t51 * t69 - t52 * t68 - t496, t100 * t51 + t101 * t50 + t136 * t37 - t15 * t168 - t198 * t8 + t2 * t223 + t68 * t72 + t703 * t9 - t498, t2 * t51 + t37 * t9 + t8 * t72 + t50 * t15 + t4 * t52 + t36 * t10 - g(1) * (-pkin(5) * t191 - qJ(6) * t190 + t431) - g(2) * (pkin(5) * t195 + qJ(6) * t194 + t440); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t517, t581 * t601, -t409 * t537 + t531, t517, t406 * t537 + t378, t509, t318 * t520 + t601 * t657 - t446 - t514, pkin(1) * t557 + t315 * t520 + (pkin(8) * t570 + g(3)) * t600 + t494 - t552, 0, 0, -t187 * t405 - t291 * t429, -t187 * t660 - t429 * t290 + (-t421 + t612) * t405 (-t291 * t406 + t362 * t548) * t580 - t450, t405 * t614 - t416, t362 * t577 + t545 + (-t290 * t406 - t362 * t595) * t580, t362 * t544, -pkin(2) * t421 + t260 * t577 + t219 * t362 + t318 * t290 + t450 * pkin(9) + (-g(3) * t548 + (-t183 * t406 - t260 * t595) * qJD(1)) * t401 + (-t211 + t677) * t660, -pkin(9) * t545 + t260 * t535 + pkin(2) * t187 - t220 * t362 - t318 * t291 + (t184 * t406 - t260 * t548) * t580 + (-pkin(9) * qJD(3) * t362 - t434) * t405, t219 * t291 - t220 * t290 + t429 * t183 + (t291 * t535 - t416) * pkin(9) + ((-t187 - t717) * pkin(9) + t685) * t405 + t420, -t183 * t219 - t184 * t220 - t260 * t318 + t434 * pkin(2) + (-t108 * t405 + (-t183 * t660 - t184 * t405) * qJD(3) + t420) * pkin(9), -t335 * t413 - t464 * t587, -t103 * t335 + t202 * t587 - t413 * t451 - t464 * t586, t335 * t438 - t464 * t544 + t466 * t587, -t103 * t451 + t202 * t586, t202 * t544 + t438 * t451 + t466 * t586, t466 * t544, -t392 * t103 - t137 * t451 + t202 * t501 + t206 * t586 - t265 * t438 - t395 * t446 + t466 * t588 - t544 * t95, t137 * t335 - t206 * t587 - t266 * t438 + t392 * t413 + t464 * t501 + t466 * t590 + t544 * t96 + t432, -t266 * t103 - t202 * t590 + t21 * t451 - t265 * t413 + t335 * t523 + t464 * t588 - t586 * t96 + t587 * t95 + t445, -g(1) * t584 - g(2) * t585 - g(3) * t493 - t137 * t392 + t206 * t501 + t21 * t266 + t265 * t523 - t588 * t95 + t590 * t96, t23, t698, t14, t419, -t664, t46, t100 * t173 + t166 * t589 + t20 * t607 + t265 * t69 + t44 * t586 - t451 * t7 + t453 * t91 + t640 * t703 + t456, -t100 * t681 + t168 * t589 + t20 * t606 - t265 * t68 - t45 * t586 + t451 * t6 - t452 * t91 - t641 * t703 + t457, t173 * t68 - t681 * t69 + t221 * t45 + t222 * t44 - t485 * t250 - t640 * t168 - t641 * t166 - t432 + (-t403 * t6 - t408 * t7 + (t403 * t44 - t408 * t45) * qJD(5)) * t335, -g(1) * t515 - g(2) * t516 - g(3) * t460 + t7 * t173 + t20 * t265 + t44 * t640 + t45 * t641 + t589 * t91 + t6 * t681, t23, t14, -t698, t46, t664, t419, -t100 * t151 + t166 * t639 + t180 * t69 - t36 * t586 + t4 * t451 + t453 * t50 + t607 * t8 - t642 * t703 + t456, -t150 * t69 - t151 * t68 + t221 * t37 - t222 * t36 + t486 * t250 + t642 * t168 - t643 * t166 - t432 + (-t2 * t403 + t4 * t408 + (-t36 * t403 - t37 * t408) * qJD(5)) * t335, t100 * t150 - t168 * t639 + t180 * t68 - t2 * t451 + t37 * t586 + t452 * t50 - t606 * t8 + t643 * t703 - t457, t2 * t150 + t8 * t180 + t4 * t151 - g(1) * (pkin(5) * t218 + qJ(6) * t217 + t515) - g(2) * (pkin(5) * t216 + qJ(6) * t215 + t516) - g(3) * (pkin(5) * t277 + qJ(6) * t276 + t460) + t639 * t50 + t643 * t37 + t642 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t613, -t290 ^ 2 + t291 ^ 2, -t187 + t614, t613, -t421 - t612, t304, -g(1) * t256 - g(2) * t439 - g(3) * t705 - t260 * t291 - t685, g(1) * t257 + g(2) * t492 + g(3) * t323 - t183 * t362 - t260 * t290 + t458, 0, 0, t615, t702, t699, -t615, t667, t438, -t105 * t466 + (-t291 * t202 + t438 * t659 + t466 * t576) * pkin(3) + t668, -t106 * t466 + (-t291 * t464 - t404 * t438 + t466 * t534) * pkin(3) + t701, -t105 * t464 + t106 * t202 + t96 * t464 - t95 * t202 + (-t404 * t103 + t659 * t414 + (t404 * t464 + (-t202 - t687) * t659) * qJD(4)) * pkin(3), -g(1) * t502 - g(2) * t430 - g(3) * t472 - t206 * t655 + t21 * t654 - t523 * t562 + (-t106 + t519) * t96 + t686 * t95, t713, -t718, t712, t715, -t442, -t690, -t53 * t703 + t391 * t69 - t686 * t166 + (t441 + t636) * t403 + (-t20 + t424) * t408 + t673, -t391 * t68 + (-t627 + t636) * t408 - t686 * t168 - t675 * t703 + t671, t54 * t166 + t53 * t168 + (-t711 + (t168 * t390 - t44) * qJD(5)) * t408 + (-t202 * t45 - t7 + (-t45 + t626) * qJD(5) + t444) * t403 + t447 + t693, -g(1) * t448 - g(2) * t417 - g(3) * t436 + t20 * t391 + t6 * t604 - t605 * t7 - t686 * t91 + t675 * t45 + (-t53 + t679) * t44, t713, t712, t718, -t690, t442, t714, t43 * t703 + t332 * t69 + t631 * t166 + (t441 + t710) * t403 + (t424 - t8) * t408 + t674, t42 * t166 + (-t43 + t540) * t168 + (-t202 * t37 + (-t37 + t626) * qJD(5) + t444) * t403 + t693 + t716, t332 * t68 + (-t42 - t539) * t703 - t631 * t168 + (-t441 - t721) * t408 + t672, t2 * t604 + t8 * t332 + t4 * t605 - g(1) * (t448 - t682) - g(2) * (t417 + t683) - g(3) * (t436 + t680) + t631 * t50 + (-t42 - t678) * t37 + (-t43 - t679) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t615, t702, t699, -t615, t667, t438, -t466 * t96 + t668, -t466 * t95 + t701, 0, 0, t713, -t718, t712, t715, -t442, -t690, -pkin(4) * t69 - t166 * t96 - t703 * t55 + (-t651 + t636) * t403 + (-t20 + t426) * t408 + t673, pkin(4) * t68 - pkin(11) * t629 - t168 * t96 + t56 * t703 + t708 * t91 + t671, t166 * t56 + t168 * t55 - t60 + (-t711 + (pkin(11) * t168 - t44) * qJD(5)) * t408 + (t490 - t7 - t637) * t403 + t447, -t20 * pkin(4) + g(1) * t239 - g(2) * t237 - g(3) * t292 - t44 * t55 - t45 * t56 - t91 * t96 + (-qJD(5) * t485 - t7 * t403 + t447) * pkin(11), t713, t712, t718, -t690, t442, t714, t703 * t48 + t353 * t69 + (-t651 + t710) * t403 + t630 * t166 + (t426 - t8) * t408 + t674, t166 * t47 - t60 + (pkin(11) * t572 - t48) * t168 + (t490 - t638) * t403 + t716, t353 * t68 + (-pkin(11) * t573 - t47) * t703 - t630 * t168 + (t651 - t721) * t408 + t672, t8 * t353 - t37 * t47 - t36 * t48 - g(1) * (-t239 - t682) - g(2) * (t237 + t683) - g(3) * (t292 + t680) + t630 * t50 + (qJD(5) * t486 + t676) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t624, -t536, t40, -t624, -t692, t100, -t168 * t91 + t425 + t637, t166 * t91 + t44 * t703 - t418, 0, 0, t624, t40, t536, t100, t692, -t624, -t117 * t166 - t415 + t637 + 0.2e1 * t652, pkin(5) * t68 - qJ(6) * t69 + (t37 - t45) * t168 + (t36 - t591) * t166, 0.2e1 * t628 + t117 * t168 - t166 * t50 + (0.2e1 * qJD(6) - t44) * t703 + t418, t2 * qJ(6) - t4 * pkin(5) - t50 * t117 - t36 * t45 - g(1) * (-pkin(5) * t194 + qJ(6) * t195) - g(2) * (-pkin(5) * t190 + qJ(6) * t191) - g(3) * (-pkin(5) * t241 + qJ(6) * t242) + t591 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t624 - t100, t40, -t703 ^ 2 - t663, t415 - t638 - t652;];
tau_reg  = t1;
