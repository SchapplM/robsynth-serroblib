% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRP10
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:35
% EndTime: 2019-03-10 02:28:26
% DurationCPUTime: 27.78s
% Computational Cost: add. (33863->1006), mult. (81753->1283), div. (0->0), fcn. (65449->14), ass. (0->431)
t407 = sin(qJ(2));
t634 = cos(pkin(6));
t549 = pkin(1) * t634;
t386 = t407 * t549;
t406 = sin(qJ(3));
t409 = cos(qJ(3));
t494 = pkin(3) * t406 - pkin(10) * t409;
t403 = sin(pkin(6));
t410 = cos(qJ(2));
t603 = t403 * t410;
t702 = -(t386 + (pkin(8) + t494) * t603) * qJD(1) + t494 * qJD(3);
t405 = sin(qJ(4));
t408 = cos(qJ(4));
t574 = qJD(1) * t403;
t542 = t410 * t574;
t509 = t409 * t542;
t543 = t407 * t574;
t273 = t405 * t509 - t408 * t543;
t570 = qJD(3) * t409;
t701 = -t405 * t570 + t273;
t567 = qJD(4) * t408;
t700 = -t406 * t567 + t701;
t524 = t634 * qJD(1);
t477 = t524 + qJD(2);
t450 = qJD(3) * t477;
t516 = t634 * qJDD(1);
t472 = t516 + qJDD(2);
t572 = qJD(2) * t410;
t539 = t406 * t572;
t564 = qJDD(1) * t407;
t166 = t403 * (qJD(1) * (t407 * t570 + t539) + t406 * t564) + t406 * t450 - t409 * t472;
t159 = qJDD(4) + t166;
t155 = qJDD(5) + t159;
t284 = t406 * t477 + t409 * t543;
t563 = qJDD(1) * t410;
t375 = t403 * t563;
t565 = qJD(1) * qJD(2);
t529 = t407 * t565;
t504 = t403 * t529;
t307 = qJDD(3) - t375 + t504;
t356 = -qJD(3) + t542;
t568 = qJD(4) * t405;
t528 = t410 * t565;
t456 = t528 + t564;
t418 = t403 * t456 + t450;
t511 = t406 * t543;
t438 = qJD(3) * t511 - t406 * t472;
t672 = t418 * t409 - t438;
t105 = t284 * t568 - t405 * t307 - (-qJD(4) * t356 + t672) * t408;
t596 = t406 * t408;
t508 = pkin(1) * t524;
t309 = -pkin(8) * t543 + t410 * t508;
t496 = pkin(2) * t407 - pkin(9) * t410;
t310 = t496 * t574;
t211 = t409 * t309 + t406 * t310;
t185 = pkin(10) * t543 + t211;
t654 = pkin(3) * t409;
t495 = pkin(10) * t406 + t654;
t349 = -pkin(2) - t495;
t571 = qJD(3) * t406;
t585 = -t408 * t185 + t349 * t567 + (-t408 * t571 - t409 * t568) * pkin(9) + t702 * t405;
t591 = t409 * t410;
t597 = t405 * t407;
t274 = (t408 * t591 + t597) * t574;
t488 = t408 * t570 - t274;
t699 = -t406 * t568 + t488;
t559 = pkin(9) * t571;
t698 = t702 * t408 + (t185 + t559) * t405;
t593 = t408 * t409;
t388 = pkin(9) * t593;
t510 = t406 * t542;
t696 = -pkin(4) * t510 + pkin(11) * t274 + (pkin(4) * t406 - pkin(11) * t593) * qJD(3) + (-t388 + (pkin(11) * t406 - t349) * t405) * qJD(4) + t698;
t282 = -t409 * t477 + t511;
t276 = qJD(4) + t282;
t269 = qJD(5) + t276;
t404 = sin(qJ(5));
t215 = -t405 * t284 - t356 * t408;
t216 = t284 * t408 - t356 * t405;
t656 = cos(qJ(5));
t465 = t404 * t215 + t216 * t656;
t551 = t284 * t567 - t356 * t568 + t405 * t672;
t471 = t307 * t408 - t551;
t38 = qJD(5) * t465 - t404 * t105 - t656 * t471;
t691 = -t269 * t465 + t38;
t695 = -pkin(11) * t700 - t585;
t604 = t403 * t407;
t322 = t406 * t634 + t409 * t604;
t468 = -t322 * t408 + t405 * t603;
t694 = pkin(11) * t468;
t693 = pkin(4) * t405 + pkin(9);
t121 = -t656 * t215 + t216 * t404;
t577 = pkin(8) * t603 + t386;
t305 = t634 * pkin(9) + t577;
t257 = qJD(2) * pkin(9) + qJD(1) * t305;
t475 = -pkin(2) * t410 - pkin(9) * t407 - pkin(1);
t275 = t475 * t574;
t164 = -t406 * t257 + t275 * t409;
t143 = pkin(3) * t356 - t164;
t111 = -pkin(4) * t215 + t143;
t51 = pkin(5) * t121 - qJ(6) * t465 + t111;
t692 = t121 * t51;
t690 = t111 * t121;
t689 = t121 * t465;
t601 = t404 * t408;
t338 = t405 * t656 + t601;
t667 = qJD(4) + qJD(5);
t239 = t667 * t338;
t507 = t656 * t570;
t589 = t239 * t406 + t274 * t656 - t701 * t404 - t408 * t507;
t531 = t656 * qJD(5);
t532 = t656 * qJD(4);
t598 = t405 * t406;
t553 = t404 * t598;
t588 = t405 * t507 - qJD(5) * t553 + (t532 + t531) * t596 - t273 * t656 + t699 * t404;
t688 = t338 * t282 + t239;
t544 = t656 * t408;
t602 = t404 * t405;
t337 = -t544 + t602;
t506 = t408 * t532;
t586 = t337 * t282 - t408 * t531 + t602 * t667 - t506;
t256 = -pkin(2) * t477 - t309;
t138 = t282 * pkin(3) - t284 * pkin(10) + t256;
t165 = t257 * t409 + t275 * t406;
t144 = -pkin(10) * t356 + t165;
t473 = qJD(2) * t508;
t502 = pkin(1) * t516;
t550 = pkin(8) * t375 + t407 * t502 + t410 * t473;
t217 = -pkin(8) * t504 + t550;
t195 = pkin(9) * t472 + t217;
t467 = t496 * qJD(2);
t206 = (qJD(1) * t467 + qJDD(1) * t475) * t403;
t88 = t409 * t195 + t406 * t206 - t257 * t571 + t275 * t570;
t75 = pkin(10) * t307 + t88;
t527 = t403 * t564;
t513 = t407 * t473 - t410 * t502 + (t403 * t528 + t527) * pkin(8);
t196 = -pkin(2) * t472 + t513;
t84 = t166 * pkin(3) - pkin(10) * t672 + t196;
t460 = -t138 * t567 + t144 * t568 - t405 * t84 - t408 * t75;
t91 = t408 * t138 - t144 * t405;
t686 = -t91 * t276 - t460;
t669 = t510 - t571;
t393 = pkin(4) * t408 + pkin(3);
t411 = -pkin(11) - pkin(10);
t684 = -t393 * t409 + t406 * t411;
t330 = -pkin(8) * t604 + t410 * t549;
t313 = qJD(2) * t330;
t658 = t465 ^ 2;
t533 = t121 ^ 2 - t658;
t566 = qJD(5) * t404;
t37 = t656 * t105 - t215 * t531 + t216 * t566 - t404 * t471;
t683 = t121 * t269 - t37;
t78 = pkin(5) * t465 + qJ(6) * t121;
t315 = t338 * t406;
t461 = t356 * t406;
t682 = t121 * t461 - t155 * t315 - t269 * t588 + t38 * t409;
t657 = cos(qJ(1));
t499 = t634 * t657;
t655 = sin(qJ(1));
t324 = t407 * t499 + t410 * t655;
t547 = t403 * t657;
t246 = t324 * t409 - t406 * t547;
t323 = t407 * t655 - t410 * t499;
t402 = qJ(4) + qJ(5);
t395 = sin(t402);
t396 = cos(t402);
t174 = t246 * t395 - t323 * t396;
t175 = t246 * t396 + t323 * t395;
t612 = t323 * t408;
t681 = t246 * t405 - t612;
t614 = t323 * t405;
t680 = t246 * t408 + t614;
t399 = t403 ^ 2;
t679 = 0.2e1 * t399;
t336 = t408 * t349;
t652 = pkin(9) * t405;
t231 = -pkin(11) * t596 + t336 + (-pkin(4) - t652) * t409;
t292 = t405 * t349 + t388;
t253 = -pkin(11) * t598 + t292;
t648 = t231 * t531 - t253 * t566 + t696 * t404 - t695 * t656;
t358 = t411 * t405;
t359 = t411 * t408;
t263 = t404 * t358 - t359 * t656;
t548 = qJD(4) * t411;
t340 = t405 * t548;
t199 = pkin(3) * t284 + pkin(10) * t282;
t109 = -t164 * t405 + t408 * t199;
t87 = pkin(11) * t282 * t408 + pkin(4) * t284 + t109;
t110 = t408 * t164 + t405 * t199;
t620 = t282 * t405;
t97 = pkin(11) * t620 + t110;
t636 = qJD(5) * t263 - t411 * t506 + t656 * t87 + (t340 - t97) * t404;
t92 = t138 * t405 + t144 * t408;
t23 = -qJD(4) * t92 - t405 * t75 + t408 * t84;
t677 = -t92 * t276 - t23;
t500 = -t165 + (t568 + t620) * pkin(4);
t151 = t155 * qJ(6);
t258 = t269 * qJD(6);
t674 = t151 + t258;
t304 = -pkin(2) * t634 - t330;
t321 = t406 * t604 - t409 * t634;
t181 = t321 * pkin(3) - t322 * pkin(10) + t304;
t578 = pkin(2) * t603 + pkin(9) * t604;
t306 = -pkin(1) * t403 - t578;
t201 = t409 * t305 + t406 * t306;
t183 = -pkin(10) * t603 + t201;
t108 = t405 * t181 + t408 * t183;
t673 = t404 * t231 + t656 * t253;
t210 = -t406 * t309 + t310 * t409;
t184 = -pkin(3) * t543 - t210;
t394 = pkin(9) * t570;
t581 = -t700 * pkin(4) - t184 + t394;
t671 = (qJDD(2) + 0.2e1 * t516) * t403;
t670 = t509 - t570;
t668 = -t307 * t409 - t356 * t571;
t314 = t577 * qJD(2);
t152 = t155 * pkin(5);
t666 = t152 - qJDD(6);
t498 = t634 * t655;
t326 = -t407 * t498 + t410 * t657;
t546 = t403 * t655;
t250 = t326 * t409 + t406 * t546;
t325 = t407 * t657 + t410 * t498;
t178 = t250 * t395 - t325 * t396;
t228 = t322 * t395 + t396 * t603;
t16 = pkin(4) * t159 + pkin(11) * t105 + t23;
t20 = pkin(11) * t471 - t460;
t72 = -pkin(11) * t216 + t91;
t63 = pkin(4) * t276 + t72;
t73 = pkin(11) * t215 + t92;
t525 = -t656 * t16 + t404 * t20 + t73 * t531 + t63 * t566;
t429 = g(1) * t178 + g(2) * t174 + g(3) * t228 - t525;
t420 = t465 * t51 - t429 - t666;
t665 = -t111 * t465 + t429;
t664 = t121 * t284 - t155 * t337 - t269 * t688;
t249 = t326 * t406 - t409 * t546;
t519 = -t324 * t406 - t409 * t547;
t452 = g(1) * t249 - g(2) * t519 + g(3) * t321;
t662 = -t155 * t263 - t395 * t452;
t573 = qJD(2) * t407;
t538 = t408 * t573;
t540 = t403 * t572;
t242 = -qJD(3) * t321 + t409 * t540;
t623 = t242 * t405;
t448 = t403 * t538 - t623;
t419 = qJD(4) * t468 + t448;
t555 = t408 * t603;
t243 = t322 * t405 + t555;
t541 = t403 * t573;
t137 = -qJD(4) * t243 + t242 * t408 + t405 * t541;
t241 = qJD(3) * t322 + t403 * t539;
t311 = t403 * t467;
t117 = -t305 * t571 + t306 * t570 + t406 * t311 + t409 * t313;
t113 = pkin(10) * t541 + t117;
t130 = t241 * pkin(3) - t242 * pkin(10) + t314;
t50 = -t108 * qJD(4) - t113 * t405 + t408 * t130;
t29 = pkin(4) * t241 - pkin(11) * t137 + t50;
t552 = t408 * t113 + t405 * t130 + t181 * t567;
t599 = t405 * t183;
t35 = t448 * pkin(11) + (-t599 + t694) * qJD(4) + t552;
t107 = t408 * t181 - t599;
t86 = pkin(4) * t321 + t107 + t694;
t93 = -pkin(11) * t243 + t108;
t645 = t404 * t86 + t656 * t93;
t8 = -qJD(5) * t645 + t29 * t656 - t404 * t35;
t647 = -qJD(5) * t673 + t695 * t404 + t696 * t656;
t660 = t121 * t586 + t337 * t37 - t338 * t38 - t465 * t688;
t316 = t406 * t544 - t553;
t659 = t121 * t589 + t315 * t37 - t316 * t38 - t465 * t588;
t412 = qJD(1) ^ 2;
t651 = t307 * pkin(3);
t397 = t406 * pkin(9);
t650 = -t669 * qJ(6) - qJD(6) * t409 + t648;
t649 = t669 * pkin(5) - t647;
t646 = t588 * pkin(5) + t589 * qJ(6) - qJD(6) * t316 + t581;
t48 = t404 * t87 + t656 * t97;
t558 = t656 * t73;
t31 = t404 * t63 + t558;
t644 = t269 * t31;
t643 = t404 * t73;
t640 = t688 * pkin(5) + t586 * qJ(6) - qJD(6) * t338 + t500;
t463 = t358 * t656 + t404 * t359;
t160 = qJD(5) * t463 + t340 * t656 + t548 * t601;
t41 = qJ(6) * t284 + t48;
t639 = t160 - t41;
t638 = t160 - t48;
t637 = t284 * pkin(5) + t636;
t33 = t656 * t72 - t643;
t635 = pkin(4) * t531 + qJD(6) - t33;
t633 = t105 * t405;
t627 = t159 * t408;
t626 = t215 * t276;
t625 = t216 * t215;
t624 = t216 * t276;
t621 = t269 * t284;
t619 = t284 * t282;
t618 = t284 * t356;
t613 = t323 * t406;
t611 = t325 * t405;
t610 = t325 * t406;
t609 = t325 * t408;
t607 = t395 * t409;
t606 = t396 * t409;
t605 = t399 * t412;
t600 = t405 * t159;
t595 = t406 * t410;
t592 = t409 * t166;
t30 = t63 * t656 - t643;
t590 = qJD(6) - t30;
t584 = -qJD(4) * t292 + t698;
t583 = -t246 * t411 + t393 * t519;
t582 = -t249 * t393 - t250 * t411;
t580 = -t321 * t393 - t322 * t411;
t342 = pkin(4) * t598 + t397;
t576 = t657 * pkin(1) + pkin(8) * t546;
t400 = t407 ^ 2;
t401 = t410 ^ 2;
t575 = t400 - t401;
t561 = g(3) * t604;
t557 = t410 * t605;
t556 = t403 * t595;
t554 = t395 * t603;
t530 = pkin(1) * t679;
t3 = t404 * t16 + t656 * t20 + t63 * t531 - t73 * t566;
t523 = -t174 * pkin(5) + qJ(6) * t175;
t179 = t250 * t396 + t325 * t395;
t522 = -t178 * pkin(5) + qJ(6) * t179;
t200 = -t406 * t305 + t306 * t409;
t518 = t276 * t405;
t517 = t276 * t408;
t515 = t407 * t557;
t514 = t406 * t195 - t409 * t206 + t257 * t570 + t275 * t571;
t505 = t407 * t528;
t32 = t404 * t72 + t558;
t503 = pkin(4) * t566 - t32;
t501 = -pkin(1) * t655 + pkin(8) * t547;
t497 = t403 * t412 * t634;
t493 = g(1) * t174 - g(2) * t178;
t492 = g(1) * t175 - g(2) * t179;
t491 = g(1) * t519 + g(2) * t249;
t490 = -g(1) * t323 + g(2) * t325;
t489 = g(1) * t326 + g(2) * t324;
t182 = pkin(3) * t603 - t200;
t229 = t322 * t396 - t554;
t487 = g(3) * (-t228 * pkin(5) + qJ(6) * t229);
t486 = pkin(5) * t396 + qJ(6) * t395;
t148 = t243 * t656 - t404 * t468;
t149 = -t404 * t243 - t468 * t656;
t57 = qJD(5) * t149 + t404 * t137 - t419 * t656;
t485 = t121 * t57 + t148 * t38;
t482 = -t405 * t92 - t408 * t91;
t476 = 0.2e1 * t524 + qJD(2);
t474 = t326 * pkin(2) + pkin(9) * t325 + t576;
t118 = -t305 * t570 - t306 * t571 + t409 * t311 - t406 * t313;
t44 = -t404 * t93 + t656 * t86;
t7 = t404 * t29 + t656 * t35 + t86 * t531 - t566 * t93;
t145 = t231 * t656 - t404 * t253;
t76 = t514 - t651;
t459 = -pkin(10) * t159 + t143 * t276;
t458 = t121 * t588 + t315 * t38;
t457 = t121 * t688 + t337 * t38;
t455 = g(1) * t657 + g(2) * t655;
t202 = -t323 * t607 - t324 * t396;
t204 = -t325 * t607 - t326 * t396;
t264 = -t396 * t604 + t409 * t554;
t454 = -g(1) * t204 - g(2) * t202 - g(3) * t264;
t203 = -t323 * t606 + t324 * t395;
t205 = -t325 * t606 + t326 * t395;
t265 = (t395 * t407 + t396 * t591) * t403;
t453 = -g(1) * t205 - g(2) * t203 - g(3) * t265;
t451 = g(1) * t250 + g(2) * t246 + g(3) * t322;
t126 = pkin(4) * t243 + t182;
t449 = t471 * t408;
t445 = -t324 * pkin(2) - t323 * pkin(9) + t501;
t444 = t452 - t76;
t443 = -g(1) * t325 - g(2) * t323 + g(3) * t603;
t442 = -t489 - t561;
t440 = pkin(4) * t611 - t249 * t411 + t250 * t393 + t474;
t439 = -t411 * t556 + t578 + (pkin(4) * t597 + t393 * t591) * t403;
t317 = t323 * pkin(2);
t437 = t684 * t323 + t693 * t324 - t317;
t319 = t325 * pkin(2);
t436 = t684 * t325 + t693 * t326 - t319;
t435 = t443 * t406;
t434 = t406 * t514 + t88 * t409 - t489;
t56 = -t137 * t656 + t243 * t531 - t404 * t419 - t468 * t566;
t432 = t121 * t56 + t148 * t37 - t149 * t38 - t465 * t57;
t431 = t121 * t241 + t148 * t155 + t269 * t57 + t321 * t38;
t430 = g(1) * t179 + g(2) * t175 + g(3) * t229 - t3;
t428 = -pkin(4) * t614 - t246 * t393 - t411 * t519 + t445;
t427 = -t160 * t121 - t263 * t38 + t37 * t463 - t451;
t426 = pkin(10) * qJD(4) * t276 - t444;
t425 = t155 * t463 + t396 * t452;
t422 = t269 * t30 + t430;
t421 = g(3) * t555 + t405 * t451;
t46 = -pkin(4) * t471 + t76;
t77 = -(-t322 * t567 - t623) * pkin(4) + (-pkin(3) * t573 - (t410 * t568 + t538) * pkin(4)) * t403 - t118;
t392 = -pkin(4) * t656 - pkin(5);
t390 = pkin(4) * t404 + qJ(6);
t312 = t577 * qJD(1);
t302 = pkin(4) * t609;
t298 = pkin(4) * t612;
t291 = -t409 * t652 + t336;
t225 = pkin(5) * t337 - qJ(6) * t338 - t393;
t188 = pkin(5) * t315 - qJ(6) * t316 + t342;
t187 = t250 * t408 + t611;
t186 = -t250 * t405 + t609;
t141 = t409 * pkin(5) - t145;
t140 = -qJ(6) * t409 + t673;
t114 = -pkin(3) * t541 - t118;
t100 = -t155 * t409 - t269 * t461;
t99 = t155 * t321 + t241 * t269;
t64 = pkin(4) * t216 + t78;
t61 = pkin(5) * t148 - qJ(6) * t149 + t126;
t49 = -t183 * t568 + t552;
t43 = t155 * t338 - t269 * t586 - t284 * t465;
t40 = -t321 * pkin(5) - t44;
t39 = qJ(6) * t321 + t645;
t27 = t269 * qJ(6) + t31;
t26 = -t269 * pkin(5) + t590;
t21 = -t338 * t37 - t465 * t586;
t17 = -t316 * t37 - t465 * t589;
t15 = t57 * pkin(5) + t56 * qJ(6) - t149 * qJD(6) + t77;
t12 = -t149 * t37 - t465 * t56;
t11 = t155 * t316 - t269 * t589 + t37 * t409 - t461 * t465;
t10 = t149 * t155 + t241 * t465 - t269 * t56 - t321 * t37;
t9 = t38 * pkin(5) + t37 * qJ(6) - qJD(6) * t465 + t46;
t6 = -t241 * pkin(5) - t8;
t5 = qJ(6) * t241 + qJD(6) * t321 + t7;
t2 = t525 - t666;
t1 = t3 + t674;
t4 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t655 - g(2) * t657, t455, 0, 0 (qJDD(1) * t400 + 0.2e1 * t505) * t399 (t407 * t563 - t565 * t575) * t679, t407 * t671 + t476 * t540 (qJDD(1) * t401 - 0.2e1 * t505) * t399, t410 * t671 - t476 * t541, t472 * t634, -t314 * t477 + t330 * t472 - t513 * t634 + g(1) * t324 - g(2) * t326 + (-t529 + t563) * t530, -t217 * t634 - t313 * t477 - t456 * t530 - t472 * t577 + t490 ((-t309 * qJD(2) + qJDD(1) * t577 + t217) * t410 + (-qJD(2) * t312 - qJDD(1) * t330 + t513) * t407 - t455) * t403, t399 * qJDD(1) * pkin(1) ^ 2 - g(1) * t501 - g(2) * t576 + t217 * t577 - t309 * t314 + t312 * t313 - t330 * t513, t284 * t242 + t322 * t672, -t322 * t166 - t284 * t241 - t242 * t282 - t321 * t672, -t242 * t356 + t322 * t307 + (t284 * t573 - t410 * t672) * t403, t166 * t321 + t241 * t282, t241 * t356 - t307 * t321 + (t166 * t410 - t282 * t573) * t403 (-t307 * t410 - t356 * t573) * t403, g(1) * t246 - g(2) * t250 - t118 * t356 + t166 * t304 + t196 * t321 + t200 * t307 + t241 * t256 + t282 * t314 + (t164 * t573 + t410 * t514) * t403, t117 * t356 - t165 * t541 + t196 * t322 - t201 * t307 + t256 * t242 + t314 * t284 + t304 * t672 + t603 * t88 + t491, -t117 * t282 - t118 * t284 - t164 * t242 - t165 * t241 - t201 * t166 - t200 * t672 - t88 * t321 + t322 * t514 - t490, -g(1) * t445 - g(2) * t474 + t165 * t117 + t164 * t118 + t196 * t304 - t200 * t514 + t88 * t201 + t256 * t314, t105 * t468 + t137 * t216, t105 * t243 + t137 * t215 + t216 * t419 - t468 * t471, -t105 * t321 + t137 * t276 - t159 * t468 + t216 * t241, t215 * t419 - t243 * t471, -t243 * t159 + t215 * t241 + t276 * t419 + t321 * t471, t159 * t321 + t241 * t276, g(1) * t680 - g(2) * t187 + t107 * t159 - t114 * t215 - t143 * t419 - t182 * t471 + t23 * t321 + t91 * t241 + t76 * t243 + t50 * t276, -g(1) * t681 - g(2) * t186 - t182 * t105 - t108 * t159 + t114 * t216 + t143 * t137 - t92 * t241 - t49 * t276 + t460 * t321 - t76 * t468, t107 * t105 + t108 * t471 - t91 * t137 + t49 * t215 - t50 * t216 + t23 * t468 + t243 * t460 + t419 * t92 - t491, -t460 * t108 + t92 * t49 + t23 * t107 + t91 * t50 + t76 * t182 + t143 * t114 - g(1) * (-pkin(3) * t246 + pkin(10) * t519 + t445) - g(2) * (pkin(3) * t250 + pkin(10) * t249 + t474) t12, t432, t10, t485, -t431, t99, t111 * t57 + t121 * t77 + t126 * t38 + t148 * t46 + t155 * t44 + t241 * t30 + t269 * t8 - t321 * t525 + t492, -t111 * t56 - t126 * t37 + t149 * t46 - t155 * t645 - t241 * t31 - t269 * t7 - t3 * t321 + t465 * t77 - t493, -t121 * t7 - t148 * t3 + t149 * t525 + t30 * t56 - t31 * t57 + t37 * t44 - t38 * t645 - t465 * t8 - t491, -g(1) * t428 - g(2) * t440 + t111 * t77 + t46 * t126 + t3 * t645 + t30 * t8 + t31 * t7 - t44 * t525, t12, t10, -t432, t99, t431, t485, t121 * t15 + t148 * t9 - t155 * t40 - t2 * t321 - t241 * t26 - t269 * t6 + t38 * t61 + t51 * t57 + t492, -t1 * t148 - t121 * t5 + t149 * t2 - t26 * t56 - t27 * t57 - t37 * t40 - t38 * t39 + t465 * t6 - t491, t1 * t321 - t149 * t9 - t15 * t465 + t155 * t39 + t241 * t27 + t269 * t5 + t37 * t61 + t51 * t56 + t493, t1 * t39 + t27 * t5 + t9 * t61 + t51 * t15 + t2 * t40 + t26 * t6 - g(1) * (-pkin(5) * t175 - qJ(6) * t174 + t428) - g(2) * (pkin(5) * t179 + qJ(6) * t178 + t440); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t515, t575 * t605, -t410 * t497 + t527, t515, t407 * t497 + t375, t472, pkin(1) * t407 * t605 + t312 * t477 - t443 - t513, pkin(1) * t557 + t309 * t477 + (pkin(8) * t565 + g(3)) * t604 + t489 - t550, 0, 0, -t438 * t406 + (t406 * t418 - t618) * t409, -t406 * t166 + t282 * t670 + t284 * t669 + t409 * t672, -t356 * t570 + t307 * t406 + (-t284 * t407 + t356 * t591) * t574, -t282 * t461 - t592 (t282 * t407 - t356 * t595) * t574 - t668, t356 * t543, -t164 * t543 - pkin(2) * t166 + t210 * t356 - t282 * t312 + (-pkin(9) * t307 - t256 * t356) * t406 + (qJD(3) * pkin(9) * t356 - t196 - t443) * t409, -pkin(2) * t672 + pkin(9) * t668 - g(1) * t610 - g(2) * t613 + g(3) * t556 + t165 * t543 + t196 * t406 - t211 * t356 - t256 * t670 - t312 * t284, -pkin(9) * t592 + t397 * t672 + t434 - t561 + (t210 + t394) * t284 + (t211 + t559) * t282 + t669 * t165 + t670 * t164, -t196 * pkin(2) - t165 * t211 - t164 * t210 - t256 * t312 + g(1) * t319 + g(2) * t317 - g(3) * t578 + ((-t164 * t409 - t165 * t406) * qJD(3) + t434) * pkin(9), -t105 * t596 + t699 * t216, -t274 * t215 + t216 * t273 + (t215 * t408 - t216 * t405) * t570 + (t449 + t633 + (-t215 * t405 - t216 * t408) * qJD(4)) * t406, t105 * t409 + t488 * t276 + (-t216 * t356 - t276 * t568 + t627) * t406, t700 * t215 - t471 * t598, -t471 * t409 + t701 * t276 + (-t215 * t356 - t276 * t567 - t600) * t406, -t159 * t409 - t276 * t461, -t143 * t273 + t291 * t159 + t184 * t215 + t584 * t276 + t442 * t405 + (-t23 + (-pkin(9) * t215 + t143 * t405) * qJD(3) - t443 * t408) * t409 + (-pkin(9) * t471 + t143 * t567 - t356 * t91 + t76 * t405) * t406, -t143 * t274 - t292 * t159 - t184 * t216 - t585 * t276 + t442 * t408 + (-t460 + (pkin(9) * t216 + t143 * t408) * qJD(3) + t443 * t405) * t409 + (-pkin(9) * t105 - t143 * t568 + t356 * t92 + t76 * t408) * t406, t292 * t471 + t291 * t105 + t92 * t273 + t91 * t274 - t584 * t216 + t585 * t215 + t482 * t570 + (t460 * t405 - t23 * t408 + (t405 * t91 - t408 * t92) * qJD(4) - t443) * t406, -t460 * t292 + t23 * t291 - t143 * t184 - g(1) * (-pkin(10) * t610 - t325 * t654 - t319) - g(2) * (-pkin(10) * t613 - t323 * t654 - t317) - g(3) * (t495 * t603 + t578) + t585 * t92 + t584 * t91 + (t143 * t570 + t406 * t76 - t489) * pkin(9), t17, t659, t11, t458, t682, t100, t111 * t588 + t121 * t581 + t145 * t155 + t269 * t647 - t30 * t461 + t315 * t46 + t342 * t38 + t409 * t525 + t453, -t111 * t589 - t155 * t673 - t269 * t648 + t3 * t409 + t31 * t461 + t316 * t46 - t342 * t37 + t465 * t581 - t454, -t121 * t648 + t145 * t37 - t3 * t315 + t30 * t589 - t31 * t588 + t316 * t525 - t38 * t673 - t465 * t647 - t435, -g(1) * t436 - g(2) * t437 - g(3) * t439 + t111 * t581 - t145 * t525 + t3 * t673 + t30 * t647 + t31 * t648 + t46 * t342, t17, t11, -t659, t100, -t682, t458, t121 * t646 - t141 * t155 + t188 * t38 + t2 * t409 + t26 * t461 - t269 * t649 + t315 * t9 + t51 * t588 + t453, -t1 * t315 - t121 * t650 - t140 * t38 - t141 * t37 + t2 * t316 - t26 * t589 - t27 * t588 + t465 * t649 - t435, -t1 * t409 + t140 * t155 + t188 * t37 + t269 * t650 - t27 * t461 - t316 * t9 - t465 * t646 + t51 * t589 + t454, t1 * t140 + t9 * t188 + t2 * t141 - g(1) * (pkin(5) * t205 + qJ(6) * t204 + t436) - g(2) * (pkin(5) * t203 + qJ(6) * t202 + t437) - g(3) * (pkin(5) * t265 + qJ(6) * t264 + t439) + t646 * t51 + t650 * t27 + t649 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t619, -t282 ^ 2 + t284 ^ 2, -t282 * t356 + t672, -t619, -t166 - t618, t307, -t165 * t356 - t256 * t284 + t452 - t514, -t164 * t356 + t256 * t282 + t451 - t88, 0, 0, t216 * t517 - t633 (-t105 + t626) * t408 + (t471 - t624) * t405, -t216 * t284 + t276 * t517 + t600, -t215 * t518 + t449, -t215 * t284 - t276 * t518 + t627, -t276 * t284, -pkin(3) * t551 - t109 * t276 - t91 * t284 + t165 * t215 + t459 * t405 + (-t426 + t651) * t408, pkin(3) * t105 + t110 * t276 - t165 * t216 + t284 * t92 + t405 * t426 + t408 * t459, t109 * t216 - t110 * t215 + ((qJD(4) * t216 + t471) * pkin(10) + t686) * t408 + ((-qJD(4) * t215 - t105) * pkin(10) + t677) * t405 - t451, -t91 * t109 - t92 * t110 - t143 * t165 + t444 * pkin(3) + (qJD(4) * t482 - t23 * t405 - t408 * t460 - t451) * pkin(10), t21, t660, t43, t457, t664, -t621, t111 * t688 + t121 * t500 - t269 * t636 - t284 * t30 + t337 * t46 - t38 * t393 + t425, -t111 * t586 - t269 * t638 + t284 * t31 + t338 * t46 + t37 * t393 + t465 * t500 + t662, t121 * t48 - t3 * t337 + t30 * t586 - t31 * t688 + t338 * t525 + t465 * t636 + t427, -g(1) * t582 - g(2) * t583 - g(3) * t580 + t111 * t500 + t3 * t263 - t30 * t636 + t31 * t638 - t46 * t393 - t463 * t525, t21, t43, -t660, -t621, -t664, t457, t121 * t640 + t225 * t38 + t26 * t284 - t269 * t637 + t337 * t9 + t51 * t688 + t425, -t1 * t337 + t121 * t41 + t2 * t338 - t26 * t586 - t27 * t688 + t465 * t637 + t427, t225 * t37 + t269 * t639 - t27 * t284 - t338 * t9 - t465 * t640 + t51 * t586 - t662, t1 * t263 + t9 * t225 - t2 * t463 - g(1) * (-t249 * t486 + t582) - g(2) * (t486 * t519 + t583) - g(3) * (-t321 * t486 + t580) + t640 * t51 + t639 * t27 + t637 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t625, -t215 ^ 2 + t216 ^ 2, -t105 - t626, t625, t471 + t624, t159, -g(1) * t186 + g(2) * t681 + g(3) * t243 - t143 * t216 - t677, g(1) * t187 + g(2) * t680 - g(3) * t468 - t143 * t215 - t686, 0, 0, t689, -t533, t683, -t689, -t691, t155, t32 * t269 + (-t121 * t216 + t155 * t656 - t269 * t566) * pkin(4) + t665, t690 + t33 * t269 + (-t155 * t404 - t216 * t465 - t269 * t531) * pkin(4) + t430, t31 * t465 + t33 * t121 - t121 * t30 - t32 * t465 + (t656 * t37 - t38 * t404 + (-t121 * t656 + t404 * t465) * qJD(5)) * pkin(4), -g(1) * t302 - g(2) * t298 + t30 * t32 - t31 * t33 + (-t525 * t656 - t111 * t216 + t3 * t404 + (-t30 * t404 + t31 * t656) * qJD(5) + t421) * pkin(4), t689, t683, t533, t155, t691, -t689, -t121 * t64 - t155 * t392 - t269 * t503 - t420, -t37 * t392 - t38 * t390 + (t27 + t503) * t465 + (-t635 + t26) * t121, t155 * t390 + t269 * t635 + t465 * t64 - t430 + t674 - t692, t1 * t390 + t2 * t392 - t51 * t64 - t26 * t32 - g(1) * (t302 + t522) - g(2) * (t298 + t523) - t487 + t635 * t27 + (t26 * t566 + t421) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t689, -t533, t683, -t689, -t691, t155, t644 + t665, t422 + t690, 0, 0, t689, t683, t533, t155, t691, -t689, -t121 * t78 + t152 - t420 + t644, pkin(5) * t37 - qJ(6) * t38 + (t27 - t31) * t465 + (t26 - t590) * t121, t465 * t78 + 0.2e1 * t151 + 0.2e1 * t258 - t422 - t692, -t2 * pkin(5) - g(1) * t522 - g(2) * t523 + t1 * qJ(6) - t26 * t31 + t27 * t590 - t51 * t78 - t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t689 - t155, t683, -t269 ^ 2 - t658, -t269 * t27 + t420;];
tau_reg  = t4;
