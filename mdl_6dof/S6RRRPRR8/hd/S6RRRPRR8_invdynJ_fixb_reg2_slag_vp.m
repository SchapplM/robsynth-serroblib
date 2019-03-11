% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:31
% EndTime: 2019-03-09 18:53:30
% DurationCPUTime: 36.64s
% Computational Cost: add. (42173->985), mult. (102632->1311), div. (0->0), fcn. (83917->18), ass. (0->416)
t401 = sin(pkin(6));
t411 = cos(qJ(2));
t534 = qJD(1) * t411;
t499 = t401 * t534;
t669 = qJD(3) - t499;
t407 = sin(qJ(2));
t536 = qJD(1) * t401;
t500 = t407 * t536;
t402 = cos(pkin(6));
t535 = qJD(1) * t402;
t519 = pkin(1) * t535;
t310 = -pkin(8) * t500 + t411 * t519;
t461 = pkin(2) * t407 - pkin(9) * t411;
t311 = t461 * t536;
t406 = sin(qJ(3));
t410 = cos(qJ(3));
t230 = -t310 * t406 + t410 * t311;
t403 = -qJ(4) - pkin(9);
t486 = qJD(3) * t403;
t551 = t410 * t411;
t668 = -(pkin(3) * t407 - qJ(4) * t551) * t536 - t230 - qJD(4) * t406 + t410 * t486;
t231 = t410 * t310 + t406 * t311;
t472 = t406 * t499;
t667 = -qJ(4) * t472 - qJD(4) * t410 - t406 * t486 + t231;
t400 = sin(pkin(12));
t564 = t400 * t406;
t603 = cos(pkin(12));
t437 = t410 * t603 - t564;
t269 = t437 * t499;
t320 = t437 * qJD(3);
t543 = t269 - t320;
t483 = t603 * t406;
t334 = t400 * t410 + t483;
t540 = t669 * t334;
t544 = t400 * t668 - t667 * t603;
t560 = t401 * t411;
t624 = pkin(1) * t407;
t539 = pkin(8) * t560 + t402 * t624;
t313 = t539 * qJD(1);
t531 = qJD(3) * t406;
t462 = -t313 + (-t472 + t531) * pkin(3);
t666 = pkin(10) * t500 - t544;
t665 = t540 * pkin(4) + t543 * pkin(10) + t462;
t405 = sin(qJ(5));
t409 = cos(qJ(5));
t235 = t269 * t405 - t409 * t500;
t664 = -t320 * t405 + t235;
t528 = qJD(5) * t409;
t429 = t334 * t528 - t664;
t236 = t269 * t409 + t405 * t500;
t583 = t320 * t409;
t480 = t236 - t583;
t529 = qJD(5) * t405;
t428 = -t334 * t529 - t480;
t524 = qJD(1) * qJD(2);
t491 = t411 * t524;
t522 = qJDD(1) * t407;
t663 = t491 + t522;
t389 = pkin(3) * t410 + pkin(2);
t251 = -pkin(4) * t437 - pkin(10) * t334 - t389;
t361 = t403 * t410;
t272 = -t361 * t603 + t403 * t564;
t612 = t251 * t528 - t272 * t529 + t665 * t405 - t409 * t666;
t621 = pkin(3) * t400;
t385 = pkin(10) + t621;
t616 = pkin(11) + t385;
t485 = qJD(5) * t616;
t374 = qJD(2) + t535;
t470 = t406 * t500;
t288 = t374 * t410 - t470;
t289 = t374 * t406 + t410 * t500;
t212 = -t603 * t288 + t289 * t400;
t438 = t400 * t288 + t289 * t603;
t622 = pkin(3) * t289;
t127 = pkin(4) * t438 + pkin(10) * t212 + t622;
t274 = pkin(9) * t374 + t313;
t449 = -pkin(2) * t411 - pkin(9) * t407 - pkin(1);
t303 = t449 * t401;
t278 = qJD(1) * t303;
t196 = t274 * t410 + t278 * t406;
t166 = qJ(4) * t288 + t196;
t160 = t400 * t166;
t195 = -t274 * t406 + t410 * t278;
t165 = -qJ(4) * t289 + t195;
t93 = t165 * t603 - t160;
t53 = t405 * t127 + t409 * t93;
t588 = t212 * t405;
t662 = pkin(11) * t588 + t405 * t485 + t53;
t52 = t409 * t127 - t405 * t93;
t661 = -pkin(5) * t438 - t52 + (-pkin(11) * t212 - t485) * t409;
t660 = t405 * t666 + t665 * t409;
t659 = t529 + t588;
t273 = -pkin(2) * t374 - t310;
t210 = -pkin(3) * t288 + qJD(4) + t273;
t100 = pkin(4) * t212 - pkin(10) * t438 + t210;
t151 = pkin(3) * t669 + t165;
t484 = t603 * t166;
t88 = t400 * t151 + t484;
t83 = pkin(10) * t669 + t88;
t48 = t100 * t405 + t409 * t83;
t649 = qJD(5) + t212;
t521 = qJDD(1) * t411;
t371 = t401 * t521;
t492 = t407 * t524;
t468 = t401 * t492;
t305 = qJDD(3) - t371 + t468;
t523 = qJDD(1) * t402;
t469 = qJDD(2) + t523;
t530 = qJD(3) * t410;
t652 = t663 * t401;
t198 = qJD(3) * t470 - t374 * t530 - t406 * t469 - t652 * t410;
t475 = qJD(2) * t519;
t515 = pkin(1) * t523;
t508 = pkin(8) * t371 + t407 * t515 + t411 * t475;
t244 = -pkin(8) * t468 + t508;
t219 = pkin(9) * t469 + t244;
t443 = t461 * qJD(2);
t228 = (qJD(1) * t443 + qJDD(1) * t449) * t401;
t99 = -qJD(3) * t196 - t406 * t219 + t410 * t228;
t65 = pkin(3) * t305 + qJ(4) * t198 - qJD(4) * t289 + t99;
t532 = qJD(2) * t411;
t497 = t406 * t532;
t199 = t374 * t531 + t401 * (qJD(1) * (t407 * t530 + t497) + t406 * t522) - t410 * t469;
t434 = -t410 * t219 - t406 * t228 + t274 * t531 - t278 * t530;
t68 = -qJ(4) * t199 + qJD(4) * t288 - t434;
t34 = t400 * t65 + t603 * t68;
t32 = pkin(10) * t305 + t34;
t141 = -t198 * t603 - t400 * t199;
t473 = t652 * pkin(8) + t407 * t475 - t411 * t515;
t220 = -pkin(2) * t469 + t473;
t147 = t199 * pkin(3) + qJDD(4) + t220;
t482 = t198 * t400 - t603 * t199;
t49 = -pkin(4) * t482 - t141 * pkin(10) + t147;
t9 = -qJD(5) * t48 - t405 * t32 + t409 * t49;
t658 = t48 * t649 + t9;
t444 = -t100 * t528 - t409 * t32 - t405 * t49 + t529 * t83;
t47 = t409 * t100 - t405 * t83;
t453 = -t47 * t649 - t444;
t563 = t401 * t407;
t321 = -t402 * t410 + t406 * t563;
t657 = pkin(3) * t321;
t259 = t409 * t272;
t656 = pkin(11) * t236 - pkin(11) * t583 + (-t259 + (pkin(11) * t334 - t251) * t405) * qJD(5) + t660 + t540 * pkin(5);
t655 = t429 * pkin(11) - t612;
t178 = t405 * t438 - t409 * t669;
t180 = t405 * t669 + t409 * t438;
t404 = sin(qJ(6));
t625 = cos(qJ(6));
t111 = t625 * t178 + t180 * t404;
t442 = -t404 * t178 + t180 * t625;
t601 = t111 * t442;
t654 = t212 * t438;
t494 = t625 * qJD(6);
t653 = t409 * (t625 * qJD(5) + t494);
t503 = t625 * t405;
t337 = t404 * t409 + t503;
t632 = qJD(5) + qJD(6);
t261 = t632 * t337;
t548 = -t337 * t212 - t261;
t557 = t404 * t405;
t441 = t625 * t409 - t557;
t547 = -t441 * t212 + t557 * t632 - t653;
t546 = t667 * t400 + t603 * t668;
t626 = cos(qJ(1));
t504 = t626 * t411;
t408 = sin(qJ(1));
t554 = t407 * t408;
t326 = -t402 * t554 + t504;
t561 = t401 * t410;
t265 = -t326 * t406 + t408 * t561;
t476 = t649 * t409;
t137 = qJDD(5) - t482;
t556 = t405 * t137;
t651 = -t476 * t649 - t556;
t650 = -t195 * t669 - t434;
t375 = pkin(8) * t563;
t623 = pkin(1) * t411;
t327 = t402 * t623 - t375;
t314 = qJD(2) * t327;
t648 = -t111 ^ 2 + t442 ^ 2;
t206 = qJD(6) + t649;
t509 = t405 * t141 + t438 * t528 + t529 * t669;
t447 = t305 * t409 - t509;
t527 = qJD(6) * t404;
t77 = -t409 * t141 - t405 * t305 + t438 * t529 - t528 * t669;
t26 = t178 * t494 + t180 * t527 - t404 * t447 + t625 * t77;
t647 = t111 * t206 - t26;
t38 = -pkin(11) * t180 + t47;
t35 = pkin(5) * t649 + t38;
t39 = -pkin(11) * t178 + t48;
t6 = pkin(5) * t137 + pkin(11) * t77 + t9;
t7 = pkin(11) * t447 - t444;
t1 = t35 * t494 - t39 * t527 + t404 * t6 + t625 * t7;
t396 = qJ(3) + pkin(12);
t390 = sin(t396);
t391 = cos(t396);
t562 = t401 * t408;
t257 = t326 * t391 + t390 * t562;
t505 = t626 * t407;
t553 = t408 * t411;
t325 = t402 * t553 + t505;
t399 = qJ(5) + qJ(6);
t392 = sin(t399);
t393 = cos(t399);
t192 = t257 * t393 + t325 * t392;
t300 = t390 * t402 + t391 * t563;
t87 = t151 * t603 - t160;
t82 = -pkin(4) * t669 - t87;
t64 = t178 * pkin(5) + t82;
t324 = t402 * t505 + t553;
t507 = t401 * t626;
t253 = t324 * t391 - t390 * t507;
t323 = -t402 * t504 + t554;
t644 = t253 * t393 + t323 * t392;
t646 = t64 * t111 + g(1) * t192 + g(2) * t644 - g(3) * (-t300 * t393 + t392 * t560) - t1;
t645 = t253 * t392 - t323 * t393;
t643 = t253 * t405 - t323 * t409;
t580 = t323 * t405;
t642 = t253 * t409 + t580;
t395 = t401 ^ 2;
t640 = 0.2e1 * t395;
t331 = t616 * t405;
t332 = t616 * t409;
t249 = -t331 * t625 - t404 * t332;
t639 = qJD(6) * t249 + t661 * t404 - t662 * t625;
t250 = -t404 * t331 + t332 * t625;
t605 = -qJD(6) * t250 + t662 * t404 + t661 * t625;
t92 = t165 * t400 + t484;
t464 = pkin(5) * t659 - t92;
t638 = -t196 * t669 - t99;
t477 = t649 * t405;
t637 = t180 * t477;
t635 = -g(1) * t325 - g(2) * t323 + g(3) * t560;
t420 = t635 * t390;
t302 = pkin(9) * t402 + t539;
t221 = -t302 * t406 + t410 * t303;
t322 = t402 * t406 + t407 * t561;
t172 = -pkin(3) * t560 - qJ(4) * t322 + t221;
t222 = t410 * t302 + t406 * t303;
t186 = -qJ(4) * t321 + t222;
t109 = t400 * t172 + t603 * t186;
t104 = -pkin(10) * t560 + t109;
t238 = t321 * t603 + t322 * t400;
t239 = -t400 * t321 + t322 * t603;
t301 = t375 + (-pkin(2) - t623) * t402;
t248 = t301 + t657;
t142 = pkin(4) * t238 - pkin(10) * t239 + t248;
t60 = t409 * t104 + t405 * t142;
t545 = pkin(4) * t500 - t546;
t176 = t405 * t251 + t259;
t201 = -t257 * t405 + t325 * t409;
t552 = t409 * t411;
t510 = t401 * t552;
t634 = g(2) * t643 - g(3) * (-t300 * t405 - t510) - g(1) * t201;
t315 = t539 * qJD(2);
t191 = -t257 * t392 + t325 * t393;
t516 = t625 * t39;
t13 = t404 * t35 + t516;
t2 = -qJD(6) * t13 - t404 * t7 + t625 * t6;
t631 = -t64 * t442 - g(1) * t191 + g(2) * t645 - g(3) * (-t300 * t392 - t393 * t560) + t2;
t27 = qJD(6) * t442 - t404 * t77 - t625 * t447;
t630 = t206 * t442 - t27;
t133 = qJDD(6) + t137;
t629 = t133 * t337 - t206 * t547;
t628 = t26 * t441 - t442 * t548;
t338 = t389 * t560;
t619 = g(3) * t338;
t618 = g(3) * t401;
t617 = t409 * pkin(5);
t175 = t409 * t251 - t272 * t405;
t572 = t334 * t409;
t149 = -pkin(5) * t437 - pkin(11) * t572 + t175;
t573 = t334 * t405;
t158 = -pkin(11) * t573 + t176;
t80 = t149 * t625 - t404 * t158;
t615 = qJD(6) * t80 + t656 * t404 - t655 * t625;
t81 = t404 * t149 + t158 * t625;
t614 = -qJD(6) * t81 + t655 * t404 + t656 * t625;
t613 = t400 * t68 - t603 * t65;
t611 = -qJD(5) * t176 + t660;
t609 = t404 * t39;
t608 = t77 * t405;
t607 = t429 * pkin(5) + t545;
t73 = t405 * t447;
t604 = -t178 * t528 + t73;
t602 = t111 * t438;
t600 = t442 * t438;
t598 = t178 * t438;
t597 = t178 * t212;
t596 = t180 * t178;
t595 = t180 * t438;
t592 = t438 ^ 2;
t591 = t438 * t669;
t590 = t212 ^ 2;
t589 = t212 * t669;
t587 = t288 * t669;
t586 = t289 * t288;
t585 = t289 * t669;
t578 = t324 * t405;
t577 = t324 * t406;
t576 = t325 * t405;
t575 = t326 * t405;
t571 = t669 * t406;
t570 = t391 * t392;
t569 = t391 * t393;
t568 = t391 * t405;
t567 = t391 * t409;
t566 = t391 * t411;
t565 = t395 * qJD(1) ^ 2;
t558 = t403 * t407;
t555 = t405 * t411;
t128 = t409 * t137;
t312 = t401 * t443;
t153 = -qJD(3) * t222 + t410 * t312 - t314 * t406;
t263 = -qJD(3) * t321 + t532 * t561;
t533 = qJD(2) * t407;
t498 = t401 * t533;
t105 = pkin(3) * t498 - qJ(4) * t263 - qJD(4) * t322 + t153;
t152 = -t302 * t531 + t303 * t530 + t406 * t312 + t410 * t314;
t262 = qJD(3) * t322 + t401 * t497;
t115 = -qJ(4) * t262 - qJD(4) * t321 + t152;
t58 = t400 * t105 + t603 * t115;
t550 = -t404 * t235 + t625 * t236 + t261 * t334 - t320 * t441;
t549 = -t625 * t235 + t320 * t503 + t334 * t653 + t428 * t404 - t527 * t573;
t542 = -t323 * t389 - t324 * t403;
t541 = -t325 * t389 - t326 * t403;
t538 = pkin(1) * t626 + pkin(8) * t562;
t397 = t407 ^ 2;
t398 = t411 ^ 2;
t537 = t397 - t398;
t526 = qJD(2) - t374;
t514 = t411 * t565;
t512 = t406 * t562;
t363 = t401 * t555;
t506 = t410 * t626;
t501 = t603 * pkin(3);
t493 = pkin(1) * t640;
t488 = -pkin(1) * t408 + pkin(8) * t507;
t59 = -t104 * t405 + t409 * t142;
t479 = -t324 * t390 - t391 * t507;
t364 = t406 * t507;
t478 = t324 * t410 - t364;
t271 = -t361 * t400 - t403 * t483;
t474 = t407 * t514;
t31 = -t305 * pkin(4) + t613;
t471 = t669 * t500;
t466 = t407 * t491;
t465 = t265 * pkin(3);
t463 = t547 * t111 - t337 * t27;
t386 = -t501 - pkin(4);
t460 = pkin(4) * t391 + pkin(10) * t390;
t256 = t326 * t390 - t391 * t562;
t459 = g(1) * t479 + g(2) * t256;
t457 = g(1) * t323 - g(2) * t325;
t456 = g(1) * t326 + g(2) * t324;
t455 = t441 * t133 + t548 * t206;
t57 = t105 * t603 - t400 * t115;
t108 = t172 * t603 - t400 * t186;
t452 = pkin(3) * t512 - t325 * t403 + t326 * t389 + t538;
t388 = pkin(4) + t617;
t412 = -pkin(11) - pkin(10);
t451 = t388 * t391 - t390 * t412;
t446 = g(1) * t626 + g(2) * t408;
t445 = -t649 * t659 + t128;
t208 = t239 * t409 - t363;
t43 = pkin(5) * t238 - pkin(11) * t208 + t59;
t207 = t239 * t405 + t510;
t50 = -pkin(11) * t207 + t60;
t21 = -t404 * t50 + t43 * t625;
t22 = t404 * t43 + t50 * t625;
t103 = pkin(4) * t560 - t108;
t146 = -t404 * t207 + t208 * t625;
t439 = pkin(3) * t364 + t323 * t403 - t324 * t389 + t488;
t55 = pkin(10) * t498 + t58;
t184 = t262 * t603 + t263 * t400;
t185 = -t400 * t262 + t263 * t603;
t216 = pkin(3) * t262 + t315;
t86 = pkin(4) * t184 - pkin(10) * t185 + t216;
t18 = -t104 * t529 + t142 * t528 + t405 * t86 + t409 * t55;
t435 = -t385 * t137 + t649 * t82;
t299 = -t390 * t563 + t391 * t402;
t433 = g(1) * t256 - g(2) * t479 - g(3) * t299;
t432 = g(1) * t257 + g(2) * t253 + g(3) * t300;
t431 = t447 * t409;
t430 = -t401 * t506 - t577;
t424 = -g(3) * t563 - t456;
t423 = t430 * pkin(3);
t422 = -t220 - t635;
t421 = -pkin(9) * t305 + t273 * t669;
t54 = -pkin(4) * t498 - t57;
t19 = -qJD(5) * t60 - t405 * t55 + t409 * t86;
t419 = qJD(5) * t385 * t649 + t31 - t433;
t417 = -pkin(9) * qJD(3) * t669 + t422;
t348 = t386 - t617;
t266 = t326 * t410 + t512;
t241 = t441 * t334;
t240 = t337 * t334;
t229 = (-t305 * t411 + t533 * t669) * t401;
t227 = pkin(5) * t573 + t271;
t202 = t257 * t409 + t576;
t145 = t207 * t625 + t208 * t404;
t122 = -qJD(5) * t363 + t185 * t405 + t239 * t528 - t409 * t498;
t121 = qJD(5) * t207 - t409 * t185 - t405 * t498;
t78 = t207 * pkin(5) + t103;
t42 = qJD(6) * t146 - t404 * t121 + t122 * t625;
t41 = t121 * t625 + t404 * t122 + t207 * t494 + t208 * t527;
t36 = t122 * pkin(5) + t54;
t20 = -pkin(5) * t447 + t31;
t15 = t38 * t625 - t609;
t14 = -t404 * t38 - t516;
t12 = t35 * t625 - t609;
t11 = -pkin(11) * t122 + t18;
t10 = pkin(5) * t184 + pkin(11) * t121 + t19;
t4 = -qJD(6) * t22 + t10 * t625 - t404 * t11;
t3 = qJD(6) * t21 + t404 * t10 + t11 * t625;
t5 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t408 - g(2) * t626, t446, 0, 0 (qJDD(1) * t397 + 0.2e1 * t466) * t395 (t407 * t521 - t524 * t537) * t640 (t374 * t532 + t407 * qJDD(2) + (t491 + 0.2e1 * t522) * t402) * t401 (qJDD(1) * t398 - 0.2e1 * t466) * t395 (-t374 * t533 + t411 * qJDD(2) + (-t492 + 0.2e1 * t521) * t402) * t401, t469 * t402, -t315 * t374 + t327 * t469 - t473 * t402 + g(1) * t324 - g(2) * t326 + (-t492 + t521) * t493, -t244 * t402 - t314 * t374 - t539 * t469 - t493 * t663 - t457 ((-qJD(2) * t310 + qJDD(1) * t539 + t244) * t411 + (-qJD(2) * t313 - qJDD(1) * t327 + t473) * t407 - t446) * t401, t395 * qJDD(1) * pkin(1) ^ 2 - g(1) * t488 - g(2) * t538 + t244 * t539 - t310 * t315 + t313 * t314 - t327 * t473, -t198 * t322 + t263 * t289, t198 * t321 - t199 * t322 - t262 * t289 + t263 * t288, t263 * t669 + t305 * t322 + (t198 * t411 + t289 * t533) * t401, t199 * t321 - t262 * t288, -t262 * t669 - t305 * t321 + (t199 * t411 + t288 * t533) * t401, t229, t153 * t669 + t221 * t305 - t315 * t288 + t301 * t199 + t220 * t321 + t273 * t262 + g(1) * t478 - g(2) * t266 + (t195 * t533 - t411 * t99) * t401, -g(1) * t577 - g(2) * t265 - t152 * t669 - t301 * t198 + t220 * t322 - t222 * t305 + t273 * t263 + t315 * t289 + (-g(1) * t506 - t196 * t533 - t411 * t434) * t401, t152 * t288 - t153 * t289 - t195 * t263 - t196 * t262 + t198 * t221 - t199 * t222 + t321 * t434 - t322 * t99 + t457, -t434 * t222 + t196 * t152 + t99 * t221 + t195 * t153 + t220 * t301 + t273 * t315 - g(1) * (-pkin(2) * t324 - pkin(9) * t323 + t488) - g(2) * (pkin(2) * t326 + pkin(9) * t325 + t538) t141 * t239 + t185 * t438, -t141 * t238 - t184 * t438 - t185 * t212 + t239 * t482, t185 * t669 + t239 * t305 + (-t141 * t411 + t438 * t533) * t401, t184 * t212 - t238 * t482, -t184 * t669 - t238 * t305 + (-t212 * t533 - t411 * t482) * t401, t229, g(1) * t253 - g(2) * t257 + t108 * t305 - t482 * t248 + t147 * t238 + t184 * t210 + t212 * t216 + t669 * t57 + (t411 * t613 + t533 * t87) * t401, -t109 * t305 + t141 * t248 + t147 * t239 + t185 * t210 + t438 * t216 - t669 * t58 + (t34 * t411 - t533 * t88) * t401 + t459, -t108 * t141 + t109 * t482 - t184 * t88 - t185 * t87 - t212 * t58 - t238 * t34 + t239 * t613 - t438 * t57 + t457, -g(1) * t439 - g(2) * t452 - t108 * t613 + t34 * t109 + t147 * t248 + t210 * t216 + t87 * t57 + t88 * t58, -t121 * t180 - t208 * t77, t121 * t178 - t180 * t122 + t77 * t207 + t208 * t447, -t121 * t649 + t137 * t208 + t180 * t184 - t238 * t77, t178 * t122 - t207 * t447, -t122 * t649 - t207 * t137 - t178 * t184 + t238 * t447, t137 * t238 + t184 * t649, g(1) * t642 - g(2) * t202 - t103 * t447 + t82 * t122 + t59 * t137 + t54 * t178 + t47 * t184 + t19 * t649 + t31 * t207 + t9 * t238, -g(1) * t643 - g(2) * t201 - t103 * t77 - t82 * t121 - t60 * t137 - t18 * t649 + t54 * t180 - t48 * t184 + t31 * t208 + t444 * t238, t47 * t121 - t48 * t122 - t18 * t178 - t19 * t180 + t207 * t444 - t9 * t208 + t447 * t60 + t59 * t77 - t459, -t444 * t60 + t48 * t18 + t9 * t59 + t47 * t19 + t31 * t103 + t82 * t54 - g(1) * (-pkin(4) * t253 + pkin(10) * t479 + t439) - g(2) * (pkin(4) * t257 + pkin(10) * t256 + t452) -t146 * t26 - t41 * t442, t111 * t41 + t145 * t26 - t146 * t27 - t42 * t442, t133 * t146 + t184 * t442 - t206 * t41 - t238 * t26, t111 * t42 + t145 * t27, -t111 * t184 - t133 * t145 - t206 * t42 - t238 * t27, t133 * t238 + t184 * t206, g(1) * t644 - g(2) * t192 + t36 * t111 + t12 * t184 + t21 * t133 + t20 * t145 + t2 * t238 + t4 * t206 + t78 * t27 + t64 * t42, -g(1) * t645 - g(2) * t191 - t1 * t238 - t13 * t184 - t22 * t133 + t20 * t146 - t3 * t206 - t78 * t26 + t36 * t442 - t64 * t41, -t1 * t145 - t111 * t3 + t12 * t41 - t13 * t42 - t146 * t2 + t21 * t26 - t22 * t27 - t4 * t442 - t459, t1 * t22 + t13 * t3 + t2 * t21 + t12 * t4 + t20 * t78 + t64 * t36 - g(1) * (-pkin(5) * t580 - t253 * t388 - t412 * t479 + t439) - g(2) * (pkin(5) * t576 - t256 * t412 + t257 * t388 + t452); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t474, t537 * t565 (t526 * t534 + t522) * t401, t474, -t500 * t526 + t371, t469, t313 * t374 + t565 * t624 - t473 - t635, pkin(1) * t514 + t310 * t374 + (pkin(8) * t524 + g(3)) * t563 + t456 - t508, 0, 0, -t198 * t406 + t410 * t585 (-t198 + t587) * t410 + (-t199 - t585) * t406, t669 * t530 + t305 * t406 + (-t289 * t407 - t551 * t669) * t536, -t199 * t410 - t288 * t571, -t669 * t531 + t305 * t410 + (-t288 * t407 + t411 * t571) * t536, -t471, -pkin(2) * t199 - t195 * t500 - t230 * t669 + t288 * t313 + t406 * t421 + t410 * t417, pkin(2) * t198 + t196 * t500 + t231 * t669 - t289 * t313 - t406 * t417 + t410 * t421, t230 * t289 - t231 * t288 + ((qJD(3) * t289 - t199) * pkin(9) + t650) * t410 + ((-qJD(3) * t288 - t198) * pkin(9) + t638) * t406 + t424, -t195 * t230 - t196 * t231 - t273 * t313 + t422 * pkin(2) + (-t99 * t406 - t434 * t410 + (-t195 * t410 - t196 * t406) * qJD(3) + t424) * pkin(9), t141 * t334 - t438 * t543, t141 * t437 + t212 * t543 + t334 * t482 - t438 * t540, t305 * t334 - t438 * t500 - t543 * t669, t212 * t540 + t437 * t482, t212 * t500 + t305 * t437 - t540 * t669, -t471, -t147 * t437 + t210 * t540 + t212 * t462 - t271 * t305 + t389 * t482 - t391 * t635 - t500 * t87 + t546 * t669, -t141 * t389 + t147 * t334 - t210 * t543 - t272 * t305 + t438 * t462 + t500 * t88 - t544 * t669 + t420, t141 * t271 - t212 * t544 + t272 * t482 + t334 * t613 + t34 * t437 - t438 * t546 - t540 * t88 + t543 * t87 + t424, t34 * t272 + t613 * t271 - t147 * t389 - g(1) * t541 - g(2) * t542 - g(3) * (-t401 * t558 + t338) + t544 * t88 + t546 * t87 + t462 * t210, t180 * t428 - t572 * t77, t664 * t180 + t480 * t178 + (t431 + t608 + (t178 * t405 - t180 * t409) * qJD(5)) * t334, t128 * t334 + t180 * t540 + t428 * t649 + t437 * t77, t178 * t429 - t334 * t73, -t178 * t540 - t334 * t556 - t429 * t649 - t437 * t447, -t137 * t437 + t540 * t649, t175 * t137 - t9 * t437 - t271 * t447 + t31 * t573 - g(1) * (-t325 * t567 + t575) - g(2) * (-t323 * t567 + t578) + t540 * t47 - (t391 * t552 + t405 * t407) * t618 + t611 * t649 + t545 * t178 + t429 * t82, -t176 * t137 - t444 * t437 - t271 * t77 + t31 * t572 - g(1) * (t325 * t568 + t326 * t409) - g(2) * (t323 * t568 + t324 * t409) - t540 * t48 - (-t391 * t555 + t407 * t409) * t618 - t612 * t649 + t545 * t180 + t428 * t82, t176 * t447 + t175 * t77 + t48 * t235 + t47 * t236 + (-t405 * t48 - t409 * t47) * t320 - t611 * t180 - t612 * t178 - t420 + (t444 * t405 - t9 * t409 + (t405 * t47 - t409 * t48) * qJD(5)) * t334, -t444 * t176 + t9 * t175 + t31 * t271 - g(1) * (-t325 * t460 + t541) - g(2) * (-t323 * t460 + t542) - t619 + t545 * t82 + t612 * t48 + t611 * t47 - (t411 * t460 - t558) * t618, -t241 * t26 - t442 * t550, t111 * t550 + t240 * t26 - t241 * t27 - t442 * t549, t133 * t241 - t206 * t550 + t26 * t437 + t442 * t540, t111 * t549 + t240 * t27, -t111 * t540 - t133 * t240 - t206 * t549 + t27 * t437, -t133 * t437 + t206 * t540, t80 * t133 - t2 * t437 + t227 * t27 + t20 * t240 - g(1) * (-t325 * t569 + t326 * t392) - g(2) * (-t323 * t569 + t324 * t392) + t549 * t64 - (t392 * t407 + t393 * t566) * t618 + t614 * t206 + t540 * t12 + t607 * t111, -t81 * t133 + t1 * t437 - t227 * t26 + t20 * t241 - g(1) * (t325 * t570 + t326 * t393) - g(2) * (t323 * t570 + t324 * t393) - t550 * t64 - (-t392 * t566 + t393 * t407) * t618 - t615 * t206 - t540 * t13 + t607 * t442, -t1 * t240 - t111 * t615 + t12 * t550 - t13 * t549 - t2 * t241 + t26 * t80 - t27 * t81 - t442 * t614 - t420, t1 * t81 + t2 * t80 + t20 * t227 - g(1) * (pkin(5) * t575 - t325 * t451 + t541) - g(2) * (pkin(5) * t578 - t323 * t451 + t542) - t619 + t607 * t64 + t615 * t13 + t614 * t12 - (t451 * t411 + (pkin(5) * t405 - t403) * t407) * t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t586, -t288 ^ 2 + t289 ^ 2, -t198 - t587, t586, -t199 + t585, t305, -g(1) * t265 - g(2) * t430 + g(3) * t321 - t273 * t289 - t638, g(1) * t266 + g(2) * t478 + g(3) * t322 - t273 * t288 - t650, 0, 0, t654, -t590 + t592, t141 + t589, -t654, t482 + t591, t305, -t210 * t438 + t92 * t669 + (-t212 * t289 + t305 * t603) * pkin(3) + t433 - t613, t210 * t212 + t669 * t93 + (-t289 * t438 - t305 * t400) * pkin(3) + t432 - t34 (-t141 * t603 + t400 * t482) * pkin(3) + (t88 - t92) * t438 + (t93 - t87) * t212, -g(1) * t465 - g(2) * t423 + g(3) * t657 - t210 * t622 + t34 * t621 - t501 * t613 + t87 * t92 - t88 * t93, t180 * t476 - t608 (-t77 - t597) * t409 - t637 + t604, -t595 - t651, t178 * t477 + t431, t445 + t598, -t649 * t438, t386 * t509 - t52 * t649 - t47 * t438 - t92 * t178 + t435 * t405 + (-t305 * t386 - t419) * t409, -t180 * t92 - t386 * t77 + t405 * t419 + t409 * t435 + t438 * t48 + t53 * t649, t53 * t178 + t52 * t180 + ((qJD(5) * t180 + t447) * t385 + t453) * t409 + (-t48 * t212 - t385 * t77 - t9 + (t178 * t385 - t48) * qJD(5)) * t405 - t432, t31 * t386 - t48 * t53 - t47 * t52 - t82 * t92 - g(1) * (-pkin(4) * t256 + pkin(10) * t257 + t465) - g(2) * (pkin(4) * t479 + t253 * pkin(10) + t423) - g(3) * (pkin(4) * t299 + pkin(10) * t300 - t657) + (-t9 * t405 - t409 * t444 - t47 * t528 - t48 * t529) * t385, -t26 * t337 - t442 * t547, t463 - t628, -t600 + t629, -t111 * t548 - t27 * t441, t455 + t602, -t206 * t438, t111 * t464 - t12 * t438 + t133 * t249 - t20 * t441 + t206 * t605 + t27 * t348 + t393 * t433 - t548 * t64, t13 * t438 - t133 * t250 + t20 * t337 - t206 * t639 - t26 * t348 - t392 * t433 + t442 * t464 - t547 * t64, t1 * t441 - t111 * t639 + t12 * t547 + t13 * t548 - t2 * t337 + t249 * t26 - t250 * t27 - t442 * t605 - t432, t1 * t250 + t2 * t249 + t20 * t348 - g(1) * (-t256 * t388 - t257 * t412 + t465) - g(2) * (-t253 * t412 + t388 * t479 + t423) - g(3) * (t299 * t388 - t300 * t412 - t657) + t464 * t64 + t639 * t13 + t605 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t482 + t591, t141 - t589, -t590 - t592, t212 * t88 + t438 * t87 + t147 + t635, 0, 0, 0, 0, 0, 0, t445 - t598, -t595 + t651 (t77 - t597) * t409 + t637 + t604, t453 * t405 + t409 * t658 - t438 * t82 + t635, 0, 0, 0, 0, 0, 0, t455 - t602, -t600 - t629, t463 + t628, t1 * t337 + t12 * t548 - t13 * t547 + t2 * t441 - t438 * t64 + t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, -t178 ^ 2 + t180 ^ 2, t178 * t649 - t77, -t596, t180 * t649 + t447, t137, -t82 * t180 + t634 + t658, t82 * t178 + g(1) * t202 + g(2) * t642 - g(3) * (-t300 * t409 + t363) - t453, 0, 0, t601, t648, t647, -t601, t630, t133, -t14 * t206 + (-t111 * t180 + t133 * t625 - t206 * t527) * pkin(5) + t631, t15 * t206 + (-t133 * t404 - t180 * t442 - t206 * t494) * pkin(5) + t646, t13 * t442 + t15 * t111 - t12 * t111 + t14 * t442 + (t625 * t26 - t27 * t404 + (-t111 * t625 + t404 * t442) * qJD(6)) * pkin(5), -t12 * t14 - t13 * t15 + (t1 * t404 + t2 * t625 - t64 * t180 + (-t12 * t404 + t13 * t625) * qJD(6) + t634) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t601, t648, t647, -t601, t630, t133, t13 * t206 + t631, t12 * t206 + t646, 0, 0;];
tau_reg  = t5;
