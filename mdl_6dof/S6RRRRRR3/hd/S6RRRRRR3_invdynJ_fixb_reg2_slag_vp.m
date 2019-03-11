% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:47
% EndTime: 2019-03-10 03:43:25
% DurationCPUTime: 20.37s
% Computational Cost: add. (40246->867), mult. (88130->1099), div. (0->0), fcn. (66420->18), ass. (0->407)
t447 = cos(qJ(2));
t630 = cos(qJ(3));
t529 = t630 * t447;
t404 = qJD(1) * t529;
t442 = sin(qJ(3));
t443 = sin(qJ(2));
t549 = qJD(1) * t443;
t527 = t442 * t549;
t337 = -t404 + t527;
t358 = t442 * t447 + t443 * t630;
t339 = t358 * qJD(1);
t257 = pkin(3) * t339 + pkin(9) * t337;
t238 = pkin(2) * t549 + t257;
t450 = -pkin(8) - pkin(7);
t387 = t450 * t447;
t366 = qJD(1) * t387;
t341 = t442 * t366;
t385 = t450 * t443;
t364 = qJD(1) * t385;
t272 = t364 * t630 + t341;
t441 = sin(qJ(4));
t446 = cos(qJ(4));
t186 = t446 * t238 - t272 * t441;
t524 = qJD(3) * t630;
t507 = pkin(2) * t524;
t674 = -t441 * t507 - t186;
t187 = t441 * t238 + t446 * t272;
t673 = -t446 * t507 + t187;
t578 = t337 * t446;
t503 = t339 * pkin(4) + pkin(10) * t578;
t420 = pkin(2) * t442 + pkin(9);
t610 = -pkin(10) - t420;
t516 = qJD(4) * t610;
t672 = t446 * t516 - t503 + t674;
t605 = qJD(2) * pkin(2);
t346 = t364 + t605;
t261 = t346 * t630 + t341;
t188 = t446 * t257 - t261 * t441;
t449 = -pkin(10) - pkin(9);
t530 = qJD(4) * t449;
t671 = t446 * t530 - t188 - t503;
t579 = t337 * t441;
t537 = pkin(10) * t579;
t670 = -t441 * t516 + t537 + t673;
t189 = t441 * t257 + t446 * t261;
t669 = -t441 * t530 + t189 + t537;
t445 = cos(qJ(5));
t538 = qJD(4) + qJD(5);
t544 = qJD(5) * t445;
t546 = qJD(4) * t446;
t440 = sin(qJ(5));
t567 = t440 * t441;
t638 = t445 * t446 - t567;
t555 = -t638 * t337 - t445 * t546 - t446 * t544 + t538 * t567;
t566 = t440 * t446;
t357 = t441 * t445 + t566;
t277 = t538 * t357;
t554 = t357 * t337 + t277;
t539 = qJD(2) + qJD(3);
t508 = t446 * t539;
t289 = t339 * t441 - t508;
t291 = t446 * t339 + t441 * t539;
t204 = t289 * t445 + t291 * t440;
t439 = sin(qJ(6));
t489 = -t289 * t440 + t445 * t291;
t629 = cos(qJ(6));
t118 = t629 * t204 + t439 * t489;
t649 = -t439 * t204 + t489 * t629;
t596 = t118 * t649;
t656 = -t118 ^ 2 + t649 ^ 2;
t329 = qJD(4) + t337;
t316 = qJD(5) + t329;
t646 = pkin(11) * t489;
t627 = pkin(2) * t447;
t424 = pkin(1) + t627;
t383 = t424 * qJD(1);
t234 = pkin(3) * t337 - pkin(9) * t339 - t383;
t344 = t630 * t366;
t262 = t442 * t346 - t344;
t242 = pkin(9) * t539 + t262;
t177 = t446 * t234 - t242 * t441;
t134 = -pkin(10) * t291 + t177;
t114 = pkin(4) * t329 + t134;
t178 = t234 * t441 + t242 * t446;
t135 = -pkin(10) * t289 + t178;
t129 = t440 * t135;
t75 = t445 * t114 - t129;
t60 = t75 - t646;
t54 = pkin(5) * t316 + t60;
t663 = pkin(11) * t204;
t131 = t445 * t135;
t76 = t114 * t440 + t131;
t61 = t76 - t663;
t604 = t439 * t61;
t24 = t54 * t629 - t604;
t534 = t629 * t61;
t25 = t439 * t54 + t534;
t668 = -t24 * t118 + t25 * t649;
t509 = t442 * t539;
t496 = t443 * t509;
t515 = qJDD(1) * t630;
t540 = t447 * qJDD(1);
t505 = -t404 * t539 - t442 * t540 - t443 * t515;
t217 = qJD(1) * t496 + t505;
t433 = qJDD(2) + qJDD(3);
t547 = qJD(4) * t441;
t174 = -qJD(4) * t508 + t446 * t217 + t339 * t547 - t441 * t433;
t175 = qJD(4) * t291 - t441 * t217 - t446 * t433;
t545 = qJD(5) * t440;
t483 = -t440 * t174 + t175 * t445 - t289 * t545 + t291 * t544;
t523 = qJD(6) * t629;
t543 = qJD(6) * t439;
t85 = t445 * t174 + t440 * t175 + t289 * t544 + t291 * t545;
t28 = t204 * t523 + t439 * t483 + t489 * t543 + t629 * t85;
t311 = qJD(6) + t316;
t653 = t118 * t311 - t28;
t241 = -pkin(3) * t539 - t261;
t195 = t289 * pkin(4) + t241;
t115 = t204 * pkin(5) + t195;
t279 = t539 * t358;
t541 = t443 * qJDD(1);
t494 = t442 * t541 - t447 * t515;
t218 = qJD(1) * t279 + t494;
t216 = qJDD(4) + t218;
t542 = qJD(1) * qJD(2);
t522 = t443 * t542;
t331 = pkin(2) * t522 - qJDD(1) * t424;
t128 = pkin(3) * t218 + pkin(9) * t217 + t331;
t126 = t446 * t128;
t521 = t447 * t542;
t283 = qJDD(2) * pkin(2) - t450 * (-t521 - t541);
t288 = t450 * (-t522 + t540);
t548 = qJD(3) * t442;
t506 = -t442 * t283 + t630 * t288 - t346 * t524 - t366 * t548;
t147 = pkin(9) * t433 - t506;
t59 = -qJD(4) * t178 - t147 * t441 + t126;
t42 = pkin(4) * t216 + pkin(10) * t174 + t59;
t58 = t441 * t128 + t446 * t147 + t234 * t546 - t242 * t547;
t49 = -pkin(10) * t175 + t58;
t13 = t114 * t544 - t135 * t545 + t440 * t42 + t445 * t49;
t10 = -pkin(11) * t483 + t13;
t14 = -qJD(5) * t76 + t445 * t42 - t440 * t49;
t214 = qJDD(5) + t216;
t9 = pkin(5) * t214 + pkin(11) * t85 + t14;
t3 = t10 * t629 + t439 * t9 + t523 * t54 - t61 * t543;
t437 = qJ(4) + qJ(5);
t431 = qJ(6) + t437;
t418 = sin(t431);
t419 = cos(t431);
t448 = cos(qJ(1));
t438 = qJ(2) + qJ(3);
t430 = cos(t438);
t444 = sin(qJ(1));
t570 = t430 * t444;
t300 = t418 * t448 - t419 * t570;
t569 = t430 * t448;
t302 = t418 * t444 + t419 * t569;
t428 = sin(t438);
t414 = g(3) * t428;
t651 = g(1) * t302 - g(2) * t300 + t115 * t118 + t419 * t414 - t3;
t350 = t610 * t441;
t432 = t446 * pkin(10);
t351 = t420 * t446 + t432;
t603 = t350 * t544 - t351 * t545 + t440 * t672 - t670 * t445;
t256 = t440 * t350 + t445 * t351;
t602 = -qJD(5) * t256 + t670 * t440 + t445 * t672;
t384 = t449 * t441;
t386 = pkin(9) * t446 + t432;
t601 = t384 * t544 - t386 * t545 + t440 * t671 - t445 * t669;
t294 = t440 * t384 + t445 * t386;
t600 = -qJD(5) * t294 + t440 * t669 + t445 * t671;
t665 = t554 * pkin(11);
t664 = -t339 * pkin(5) + pkin(11) * t555;
t29 = qJD(6) * t649 - t439 * t85 + t629 * t483;
t632 = t311 * t649 - t29;
t299 = t418 * t570 + t419 * t448;
t301 = -t418 * t569 + t419 * t444;
t4 = -qJD(6) * t25 - t439 * t10 + t629 * t9;
t634 = -g(1) * t301 + g(2) * t299 - t115 * t649 + t418 * t414 + t4;
t662 = -t665 + t603;
t661 = t664 + t602;
t660 = -t665 + t601;
t659 = t664 + t600;
t591 = t204 * t489;
t559 = t357 * t543 + t439 * t554 - t523 * t638 + t555 * t629;
t269 = t357 * t629 + t439 * t638;
t558 = qJD(6) * t269 - t439 * t555 + t554 * t629;
t271 = t364 * t442 - t344;
t501 = pkin(2) * t548 - t271;
t572 = t428 * t448;
t573 = t428 * t444;
t658 = g(1) * t572 + g(2) * t573;
t498 = g(1) * t448 + g(2) * t444;
t480 = t498 * t428;
t616 = g(3) * t430;
t657 = -t480 + t616;
t278 = -qJD(2) * t529 - t447 * t524 + t496;
t589 = t278 * t441;
t477 = t358 * t546 - t589;
t655 = -t204 ^ 2 + t489 ^ 2;
t654 = -t616 + t658;
t652 = t204 * t316 - t85;
t427 = sin(t437);
t429 = cos(t437);
t305 = t427 * t448 - t429 * t570;
t307 = t427 * t444 + t429 * t569;
t650 = g(1) * t307 - g(2) * t305 + t195 * t204 + t429 * t414 - t13;
t356 = t442 * t443 - t529;
t260 = pkin(3) * t356 - pkin(9) * t358 - t424;
t295 = t442 * t385 - t387 * t630;
t196 = t446 * t260 - t295 * t441;
t575 = t358 * t446;
t160 = pkin(4) * t356 - pkin(10) * t575 + t196;
t284 = t446 * t295;
t197 = t441 * t260 + t284;
t576 = t358 * t441;
t182 = -pkin(10) * t576 + t197;
t102 = t440 * t160 + t445 * t182;
t425 = pkin(4) * t547;
t643 = pkin(5) * t554 + t425;
t611 = t446 * pkin(4);
t377 = pkin(5) * t429 + t611;
t363 = pkin(3) + t377;
t434 = -pkin(11) + t449;
t642 = t430 * t363 - t428 * t434;
t641 = t630 * t385 + t442 * t387;
t422 = pkin(3) + t611;
t640 = t430 * t422 - t428 * t449;
t499 = t430 * pkin(3) + t428 * pkin(9);
t312 = pkin(4) * t579;
t639 = t312 + t501;
t637 = -t177 * t441 + t178 * t446;
t560 = t446 * t448;
t565 = t441 * t444;
t322 = t430 * t565 + t560;
t562 = t444 * t446;
t564 = t441 * t448;
t324 = -t430 * t564 + t562;
t636 = -g(1) * t324 + g(2) * t322;
t304 = t427 * t570 + t429 * t448;
t306 = -t427 * t569 + t429 * t444;
t635 = -g(1) * t306 + g(2) * t304 + t427 * t414;
t633 = -t195 * t489 + t14 + t635;
t631 = t316 * t489 - t483;
t628 = pkin(2) * t443;
t626 = pkin(4) * t441;
t624 = pkin(11) * t357;
t388 = t448 * t424;
t618 = g(2) * t388;
t615 = g(3) * t441;
t614 = g(3) * t447;
t612 = t638 * pkin(5);
t255 = t445 * t350 - t351 * t440;
t220 = t255 - t624;
t349 = t638 * pkin(11);
t221 = t349 + t256;
t145 = t220 * t629 - t439 * t221;
t609 = qJD(6) * t145 + t439 * t661 + t629 * t662;
t146 = t439 * t220 + t221 * t629;
t608 = -qJD(6) * t146 - t439 * t662 + t629 * t661;
t292 = t445 * t384 - t386 * t440;
t239 = t292 - t624;
t240 = t349 + t294;
t180 = t239 * t629 - t439 * t240;
t607 = qJD(6) * t180 + t439 * t659 + t629 * t660;
t181 = t439 * t239 + t240 * t629;
t606 = -qJD(6) * t181 - t439 * t660 + t629 * t659;
t56 = t58 * t446;
t421 = pkin(4) * t445 + pkin(5);
t568 = t439 * t440;
t87 = -t134 * t440 - t131;
t62 = t87 + t663;
t88 = t445 * t134 - t129;
t63 = t88 - t646;
t599 = -t439 * t62 - t629 * t63 + t421 * t523 + (-t440 * t543 + (t445 * t629 - t568) * qJD(5)) * pkin(4);
t528 = t629 * t440;
t598 = t439 * t63 - t629 * t62 - t421 * t543 + (-t440 * t523 + (-t439 * t445 - t528) * qJD(5)) * pkin(4);
t597 = pkin(7) * qJDD(1);
t595 = t174 * t441;
t594 = t175 * t446;
t590 = t241 * t337;
t588 = t278 * t446;
t587 = t289 * t329;
t586 = t289 * t441;
t585 = t291 * t289;
t584 = t291 * t329;
t583 = t291 * t446;
t582 = t311 * t339;
t581 = t316 * t339;
t580 = t329 * t339;
t577 = t339 * t337;
t557 = t639 + t643;
t215 = -t312 + t262;
t556 = -t215 + t643;
t553 = t425 + t639;
t376 = pkin(5) * t427 + t626;
t552 = t376 - t450;
t435 = t443 ^ 2;
t436 = t447 ^ 2;
t551 = t435 - t436;
t550 = t435 + t436;
t536 = t443 * t605;
t535 = pkin(9) * qJD(4) * t329;
t452 = qJD(1) ^ 2;
t533 = t443 * t452 * t447;
t532 = g(1) * t569 + g(2) * t570 + t414;
t531 = qJD(2) * t450;
t526 = t358 * t547;
t520 = -t450 + t626;
t156 = t630 * t283 + t442 * t288 - t346 * t548 + t366 * t524;
t148 = -t433 * pkin(3) - t156;
t518 = -t148 - t616;
t101 = t445 * t160 - t182 * t440;
t194 = pkin(3) * t279 + pkin(9) * t278 + t536;
t365 = t443 * t531;
t367 = t447 * t531;
t210 = qJD(3) * t641 + t630 * t365 + t442 * t367;
t514 = t446 * t194 - t210 * t441;
t510 = t329 * t446;
t423 = -pkin(2) * t630 - pkin(3);
t504 = t443 * t521;
t502 = -g(1) * t573 + g(2) * t572;
t500 = -t215 + t425;
t497 = g(1) * t444 - g(2) * t448;
t236 = pkin(4) * t576 - t641;
t495 = -pkin(9) * t216 + t590;
t491 = t177 * t446 + t178 * t441;
t490 = -t216 * t420 + t590;
t487 = t363 * t428 + t430 * t434;
t485 = t422 * t428 + t430 * t449;
t484 = t148 * t441 + t178 * t339 + t241 * t546 + t430 * t615;
t482 = -t177 * t339 + t241 * t547 + t446 * t658;
t481 = -t177 * t578 - t178 * t579 - t532 + t56;
t248 = t638 * t358;
t86 = pkin(5) * t356 - pkin(11) * t248 + t101;
t247 = t357 * t358;
t91 = -pkin(11) * t247 + t102;
t38 = -t439 * t91 + t629 * t86;
t39 = t439 * t86 + t629 * t91;
t479 = -0.2e1 * pkin(1) * t542 - pkin(7) * qJDD(2);
t184 = -t439 * t247 + t248 * t629;
t476 = -t526 - t588;
t382 = t423 - t611;
t68 = pkin(10) * t588 + pkin(4) * t279 + (-t284 + (pkin(10) * t358 - t260) * t441) * qJD(4) + t514;
t98 = t441 * t194 + t446 * t210 + t260 * t546 - t295 * t547;
t79 = -pkin(10) * t477 + t98;
t22 = t160 * t544 - t182 * t545 + t440 * t68 + t445 * t79;
t100 = t175 * pkin(4) + t148;
t451 = qJD(2) ^ 2;
t470 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t451 + t497;
t469 = pkin(1) * t452 + t498 - t597;
t468 = -t383 * t337 + t506 + t532;
t467 = -qJD(4) * t491 - t59 * t441;
t23 = -qJD(5) * t102 - t440 * t79 + t445 * t68;
t466 = t383 * t339 + t156 + t654;
t268 = t357 * t439 - t629 * t638;
t465 = t24 * t559 - t25 * t558 - t3 * t268 - t4 * t269 - t532;
t464 = t13 * t638 - t14 * t357 - t554 * t76 + t555 * t75 - t532;
t47 = pkin(5) * t483 + t100;
t462 = t115 * t558 - t24 * t339 + t47 * t268 + t419 * t654;
t211 = qJD(3) * t295 + t442 * t365 - t630 * t367;
t461 = -t100 * t638 + t195 * t554 - t75 * t339 + t429 * t654;
t460 = t467 + t56;
t140 = pkin(4) * t477 + t211;
t459 = -t115 * t559 + t25 * t339 + t47 * t269 + t418 * t657;
t458 = t100 * t357 - t195 * t555 + t339 * t76 + t427 * t657;
t457 = -g(1) * (-pkin(3) * t572 + pkin(9) * t569) - g(2) * (-pkin(3) * t573 + pkin(9) * t570) - g(3) * t499;
t333 = pkin(4) * t528 + t439 * t421;
t332 = -pkin(4) * t568 + t421 * t629;
t325 = t430 * t560 + t565;
t323 = -t430 * t562 + t564;
t308 = -t422 - t612;
t298 = t382 - t612;
t222 = -t337 ^ 2 + t339 ^ 2;
t212 = qJDD(6) + t214;
t192 = -t505 + (t337 - t527) * t539;
t185 = pkin(5) * t247 + t236;
t183 = t247 * t629 + t248 * t439;
t176 = pkin(4) * t291 + pkin(5) * t489;
t111 = -t278 * t566 - t440 * t526 - t545 * t576 + (t538 * t575 - t589) * t445;
t110 = t277 * t358 + t278 * t638;
t109 = t216 * t441 - t291 * t339 + t329 * t510;
t108 = -t329 ^ 2 * t441 + t216 * t446 + t289 * t339;
t104 = t329 * t586 - t594;
t103 = t291 * t510 - t595;
t99 = -qJD(4) * t197 + t514;
t81 = t204 * t339 + t214 * t638 - t316 * t554;
t80 = t214 * t357 - t316 * t555 - t339 * t489;
t74 = t111 * pkin(5) + t140;
t55 = (-t174 - t587) * t446 + (-t175 - t584) * t441;
t51 = qJD(6) * t184 - t439 * t110 + t111 * t629;
t50 = t110 * t629 + t439 * t111 + t247 * t523 + t248 * t543;
t46 = t204 * t554 - t483 * t638;
t45 = -t357 * t85 - t489 * t555;
t44 = t118 * t339 - t212 * t268 - t311 * t558;
t43 = t212 * t269 - t311 * t559 - t339 * t649;
t27 = t60 * t629 - t604;
t26 = -t439 * t60 - t534;
t19 = -pkin(11) * t111 + t22;
t18 = t204 * t555 - t357 * t483 - t489 * t554 - t638 * t85;
t17 = pkin(5) * t279 + pkin(11) * t110 + t23;
t16 = t118 * t558 + t268 * t29;
t15 = -t269 * t28 - t559 * t649;
t7 = t118 * t559 + t268 * t28 - t269 * t29 - t558 * t649;
t6 = -qJD(6) * t39 + t17 * t629 - t439 * t19;
t5 = qJD(6) * t38 + t439 * t17 + t19 * t629;
t1 = [0, 0, 0, 0, 0, qJDD(1), t497, t498, 0, 0, qJDD(1) * t435 + 0.2e1 * t504, 0.2e1 * t443 * t540 - 0.2e1 * t542 * t551, qJDD(2) * t443 + t447 * t451, qJDD(1) * t436 - 0.2e1 * t504, qJDD(2) * t447 - t443 * t451, 0, t443 * t479 + t447 * t470, -t443 * t470 + t447 * t479, 0.2e1 * t550 * t597 - t498, -g(1) * (-pkin(1) * t444 + pkin(7) * t448) - g(2) * (pkin(1) * t448 + pkin(7) * t444) + (pkin(7) ^ 2 * t550 + pkin(1) ^ 2) * qJDD(1), -t217 * t358 - t278 * t339, t217 * t356 - t218 * t358 + t278 * t337 - t279 * t339, -t278 * t539 + t358 * t433, t218 * t356 + t279 * t337, -t279 * t539 - t356 * t433, 0, -t211 * t539 - t424 * t218 - t383 * t279 + t331 * t356 + t337 * t536 + t430 * t497 + t433 * t641, -t210 * t539 + t424 * t217 + t383 * t278 - t295 * t433 + t331 * t358 + t339 * t536 + t502, -t156 * t358 - t210 * t337 + t211 * t339 + t217 * t641 - t218 * t295 + t261 * t278 - t262 * t279 + t356 * t506 - t498, -t506 * t295 + t262 * t210 + t156 * t641 - t261 * t211 - t331 * t424 - t383 * t536 - g(1) * (-t424 * t444 - t448 * t450) - g(2) * (-t444 * t450 + t388) -t174 * t575 + t291 * t476 (t289 * t446 + t291 * t441) * t278 + (t595 - t594 + (-t583 + t586) * qJD(4)) * t358, -t174 * t356 + t216 * t575 + t279 * t291 + t329 * t476, t175 * t576 + t289 * t477, -t175 * t356 - t216 * t576 - t279 * t289 - t329 * t477, t216 * t356 + t279 * t329, -g(1) * t323 - g(2) * t325 + t148 * t576 - t175 * t641 + t177 * t279 + t196 * t216 + t211 * t289 + t241 * t477 + t329 * t99 + t356 * t59, -g(1) * t322 - g(2) * t324 + t148 * t575 + t174 * t641 - t178 * t279 - t197 * t216 + t211 * t291 + t241 * t476 - t329 * t98 - t356 * t58, t174 * t196 - t175 * t197 - t289 * t98 - t291 * t99 + t491 * t278 + (-qJD(4) * t637 - t441 * t58 - t446 * t59) * t358 - t502, -t618 - t148 * t641 + t177 * t99 + t178 * t98 + t59 * t196 + t58 * t197 + t241 * t211 + (g(1) * t450 - g(2) * t499) * t448 + (-g(1) * (-t424 - t499) + g(2) * t450) * t444, -t110 * t489 - t248 * t85, t110 * t204 - t111 * t489 + t85 * t247 - t248 * t483, -t110 * t316 + t214 * t248 + t279 * t489 - t356 * t85, t204 * t111 + t247 * t483, -t111 * t316 - t204 * t279 - t247 * t214 - t356 * t483, t214 * t356 + t279 * t316, -g(1) * t305 - g(2) * t307 + t100 * t247 + t101 * t214 + t195 * t111 + t14 * t356 + t140 * t204 + t23 * t316 + t236 * t483 + t75 * t279, -g(1) * t304 - g(2) * t306 + t100 * t248 - t102 * t214 - t110 * t195 - t13 * t356 + t140 * t489 - t22 * t316 - t236 * t85 - t279 * t76, t101 * t85 - t102 * t483 + t75 * t110 - t76 * t111 - t13 * t247 - t14 * t248 - t22 * t204 - t23 * t489 - t502, -t618 + t100 * t236 + t14 * t101 + t13 * t102 + t195 * t140 + t76 * t22 + t75 * t23 + (-g(1) * t520 - g(2) * t640) * t448 + (-g(1) * (-t424 - t640) - g(2) * t520) * t444, -t184 * t28 - t50 * t649, t118 * t50 + t183 * t28 - t184 * t29 - t51 * t649, t184 * t212 + t279 * t649 - t28 * t356 - t311 * t50, t118 * t51 + t183 * t29, -t118 * t279 - t183 * t212 - t29 * t356 - t311 * t51, t212 * t356 + t279 * t311, -g(1) * t300 - g(2) * t302 + t115 * t51 + t118 * t74 + t183 * t47 + t185 * t29 + t212 * t38 + t24 * t279 + t311 * t6 + t356 * t4, -g(1) * t299 - g(2) * t301 - t115 * t50 + t184 * t47 - t185 * t28 - t212 * t39 - t25 * t279 - t3 * t356 - t311 * t5 + t649 * t74, -t118 * t5 - t183 * t3 - t184 * t4 + t24 * t50 - t25 * t51 + t28 * t38 - t29 * t39 - t6 * t649 - t502, -t618 + t115 * t74 + t47 * t185 + t24 * t6 + t25 * t5 + t3 * t39 + t4 * t38 + (-g(1) * t552 - g(2) * t642) * t448 + (-g(1) * (-t424 - t642) - g(2) * t552) * t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t533, t551 * t452, t541, t533, t540, qJDD(2), t443 * t469 - t614, g(3) * t443 + t447 * t469, 0, 0, t577, t222, t192, -t577, -t494, t433, t271 * t539 + (-qJD(3) * t509 - t337 * t549 + t433 * t630) * pkin(2) + t466, t272 * t539 + (-t339 * t549 - t442 * t433 - t524 * t539) * pkin(2) + t468 (t262 - t271) * t339 + (-t261 + t272) * t337 + (t630 * t217 - t218 * t442 + (-t337 * t630 + t339 * t442) * qJD(3)) * pkin(2), t261 * t271 - t262 * t272 + (t630 * t156 - t614 - t506 * t442 + (-t261 * t442 + t262 * t630) * qJD(3) + (qJD(1) * t383 + t498) * t443) * pkin(2), t103, t55, t109, t104, t108, -t580, t423 * t175 + t518 * t446 + t490 * t441 + t501 * t289 + (-t420 * t546 + t674) * t329 + t482, -t423 * t174 + t490 * t446 - t441 * t480 + t501 * t291 + (t420 * t547 + t673) * t329 + t484, t186 * t291 + t187 * t289 + (-t289 * t507 - t175 * t420 + (t291 * t420 - t177) * qJD(4)) * t446 + (t291 * t507 - t174 * t420 - t59 + (t289 * t420 - t178) * qJD(4)) * t441 + t481, t148 * t423 - t178 * t187 - t177 * t186 - t241 * t271 + (-t614 + t498 * t443 + (t241 * t442 + t630 * t637) * qJD(3)) * pkin(2) + t460 * t420 + t457, t45, t18, t80, t46, t81, -t581, t204 * t553 + t255 * t214 + t316 * t602 + t382 * t483 + t461, -t214 * t256 - t316 * t603 - t382 * t85 + t489 * t553 + t458, -t204 * t603 + t255 * t85 - t256 * t483 - t489 * t602 + t464, t13 * t256 + t14 * t255 + t100 * t382 - g(3) * (t640 + t627) + t603 * t76 + t602 * t75 + t553 * t195 + t498 * (t485 + t628) t15, t7, t43, t16, t44, -t582, t118 * t557 + t145 * t212 + t29 * t298 + t311 * t608 + t462, -t146 * t212 - t28 * t298 - t311 * t609 + t557 * t649 + t459, -t118 * t609 + t145 * t28 - t146 * t29 - t608 * t649 + t465, t3 * t146 + t4 * t145 + t47 * t298 - g(3) * (t642 + t627) + t609 * t25 + t608 * t24 + t557 * t115 + t498 * (t487 + t628); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, t222, t192, -t577, -t494, t433, t262 * t539 + t466, t261 * t539 + t468, 0, 0, t103, t55, t109, t104, t108, -t580, -pkin(3) * t175 - t188 * t329 - t262 * t289 + t495 * t441 + (t518 - t535) * t446 + t482, pkin(3) * t174 + t189 * t329 - t262 * t291 + t495 * t446 + (-t480 + t535) * t441 + t484, t188 * t291 + t189 * t289 + (-t595 - t594 + (t583 + t586) * qJD(4)) * pkin(9) + t467 + t481, -t148 * pkin(3) + pkin(9) * t460 - t177 * t188 - t178 * t189 - t241 * t262 + t457, t45, t18, t80, t46, t81, -t581, t204 * t500 + t292 * t214 + t316 * t600 - t422 * t483 + t461, -t214 * t294 - t316 * t601 + t422 * t85 + t489 * t500 + t458, -t204 * t601 + t292 * t85 - t294 * t483 - t489 * t600 + t464, -g(3) * t640 - t100 * t422 + t13 * t294 + t14 * t292 + t195 * t500 + t485 * t498 + t600 * t75 + t601 * t76, t15, t7, t43, t16, t44, -t582, t118 * t556 + t180 * t212 + t29 * t308 + t311 * t606 + t462, -t181 * t212 - t28 * t308 - t311 * t607 + t556 * t649 + t459, -t118 * t607 + t180 * t28 - t181 * t29 - t606 * t649 + t465, -g(3) * t642 + t115 * t556 + t4 * t180 + t3 * t181 + t24 * t606 + t25 * t607 + t47 * t308 + t487 * t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t585, -t289 ^ 2 + t291 ^ 2, -t174 + t587, -t585, -t175 + t584, t216, -t242 * t546 + t178 * t329 - t241 * t291 + t126 + (-qJD(4) * t234 - t147 + t414) * t441 + t636, g(1) * t325 - g(2) * t323 + t177 * t329 + t241 * t289 + t414 * t446 - t58, 0, 0, t591, t655, t652, -t591, t631, t214, -t316 * t87 + (-t204 * t291 + t214 * t445 - t316 * t545) * pkin(4) + t633, t316 * t88 + (-t214 * t440 - t291 * t489 - t316 * t544) * pkin(4) + t650, t76 * t489 + t88 * t204 - t75 * t204 + t87 * t489 + (-t440 * t483 + t445 * t85 + (-t204 * t445 + t440 * t489) * qJD(5)) * pkin(4), -t75 * t87 - t76 * t88 + (t13 * t440 + t14 * t445 - t195 * t291 + t428 * t615 + (-t440 * t75 + t445 * t76) * qJD(5) + t636) * pkin(4), t596, t656, t653, -t596, t632, t212, -t176 * t118 + t332 * t212 + t311 * t598 + t634, -t176 * t649 - t333 * t212 - t311 * t599 + t651, -t118 * t599 + t28 * t332 - t29 * t333 - t598 * t649 + t668, t3 * t333 + t4 * t332 - t115 * t176 - g(1) * (-t376 * t569 + t377 * t444) - g(2) * (-t376 * t570 - t377 * t448) + t376 * t414 + t599 * t25 + t598 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t591, t655, t652, -t591, t631, t214, t316 * t76 + t633, t316 * t75 + t650, 0, 0, t596, t656, t653, -t596, t632, t212, -t26 * t311 + (-t118 * t489 + t212 * t629 - t311 * t543) * pkin(5) + t634, t27 * t311 + (-t212 * t439 - t311 * t523 - t489 * t649) * pkin(5) + t651, t27 * t118 + t26 * t649 + (t629 * t28 - t29 * t439 + (-t118 * t629 + t439 * t649) * qJD(6)) * pkin(5) + t668, -t24 * t26 - t25 * t27 + (t3 * t439 + t4 * t629 - t115 * t489 + (-t24 * t439 + t25 * t629) * qJD(6) + t635) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, t656, t653, -t596, t632, t212, t25 * t311 + t634, t24 * t311 + t651, 0, 0;];
tau_reg  = t1;
