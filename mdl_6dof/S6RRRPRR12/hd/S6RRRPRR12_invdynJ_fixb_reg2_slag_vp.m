% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRPRR12
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:40
% EndTime: 2019-03-09 19:47:38
% DurationCPUTime: 33.43s
% Computational Cost: add. (43511->1018), mult. (105265->1371), div. (0->0), fcn. (86202->18), ass. (0->388)
t432 = sin(qJ(2));
t590 = cos(pkin(6));
t522 = pkin(1) * t590;
t408 = t432 * t522;
t431 = sin(qJ(3));
t434 = cos(qJ(3));
t476 = pkin(3) * t431 - qJ(4) * t434;
t426 = sin(pkin(6));
t435 = cos(qJ(2));
t563 = t426 * t435;
t647 = (t408 + (pkin(8) + t476) * t563) * qJD(1) - qJD(3) * t476 + qJD(4) * t431;
t499 = t590 * qJD(1);
t470 = t499 + qJD(2);
t542 = qJD(1) * t426;
t518 = t432 * t542;
t313 = t431 * t470 + t434 * t518;
t541 = qJD(1) * t435;
t517 = t426 * t541;
t384 = -qJD(3) + t517;
t425 = sin(pkin(12));
t427 = cos(pkin(12));
t253 = t313 * t425 + t427 * t384;
t255 = t313 * t427 - t384 * t425;
t430 = sin(qJ(5));
t433 = cos(qJ(5));
t157 = t253 * t433 + t430 * t255;
t159 = t253 * t430 - t255 * t433;
t429 = sin(qJ(6));
t610 = cos(qJ(6));
t87 = t610 * t157 - t429 * t159;
t646 = t87 ^ 2;
t89 = t429 * t157 + t159 * t610;
t645 = t89 ^ 2;
t494 = t431 * t518;
t311 = -t434 * t470 + t494;
t303 = qJD(5) + t311;
t297 = qJD(6) + t303;
t644 = t297 * t87;
t643 = t89 * t297;
t492 = pkin(1) * t499;
t333 = -pkin(8) * t518 + t435 * t492;
t483 = pkin(2) * t432 - pkin(9) * t435;
t334 = t483 * t542;
t245 = t434 * t333 + t431 * t334;
t210 = qJ(4) * t518 + t245;
t642 = t427 * t210 + t647 * t425;
t538 = qJD(3) * t431;
t528 = pkin(9) * t538;
t551 = -t647 * t427 + (t210 + t528) * t425;
t565 = t425 * t434;
t295 = -t427 * t518 + t517 * t565;
t537 = qJD(3) * t434;
t479 = t425 * t537 - t295;
t559 = t434 * t435;
t566 = t425 * t432;
t296 = (t427 * t559 + t566) * t542;
t641 = -t427 * t537 + t296;
t516 = t431 * t541;
t493 = t426 * t516;
t561 = t427 * t434;
t640 = -pkin(4) * t493 + pkin(10) * t296 + (pkin(4) * t431 - pkin(10) * t561) * qJD(3) + t551;
t562 = t427 * t431;
t639 = -pkin(10) * t295 - (-pkin(9) * t562 - pkin(10) * t565) * qJD(3) + t642;
t638 = t157 ^ 2;
t637 = t159 ^ 2;
t636 = t157 * t303;
t635 = t159 * t303;
t364 = t425 * t433 + t427 * t430;
t346 = t364 * qJD(5);
t553 = t346 * t431 + t479 * t430 + t641 * t433;
t535 = qJD(5) * t433;
t536 = qJD(5) * t430;
t617 = -t425 * t536 + t427 * t535;
t552 = -t433 * t295 - t296 * t430 + t364 * t537 + t431 * t617;
t548 = t364 * t311 + t346;
t560 = t433 * t427;
t363 = t425 * t430 - t560;
t549 = t363 * t311 - t617;
t477 = pkin(3) * t434 + qJ(4) * t431;
t376 = -pkin(2) - t477;
t358 = t427 * t376;
t268 = -pkin(10) * t562 + t358 + (-pkin(9) * t425 - pkin(4)) * t434;
t309 = pkin(9) * t561 + t425 * t376;
t567 = t425 * t431;
t280 = -pkin(10) * t567 + t309;
t594 = t268 * t535 - t280 * t536 + t430 * t640 - t639 * t433;
t177 = t430 * t268 + t433 * t280;
t593 = -qJD(5) * t177 + t639 * t430 + t433 * t640;
t545 = pkin(8) * t563 + t408;
t324 = t590 * pkin(9) + t545;
t293 = qJD(2) * pkin(9) + qJD(1) * t324;
t468 = -pkin(2) * t435 - pkin(9) * t432 - pkin(1);
t301 = t468 * t542;
t195 = -t431 * t293 + t301 * t434;
t225 = pkin(3) * t313 + qJ(4) * t311;
t134 = -t195 * t425 + t427 * t225;
t578 = t311 * t427;
t104 = pkin(4) * t313 + pkin(10) * t578 + t134;
t135 = t427 * t195 + t425 * t225;
t579 = t311 * t425;
t121 = pkin(10) * t579 + t135;
t600 = pkin(10) + qJ(4);
t379 = t600 * t425;
t380 = t600 * t427;
t592 = qJD(4) * t560 - t433 * t121 - t379 * t535 + (-qJD(4) * t425 - qJD(5) * t380 - t104) * t430;
t291 = -t430 * t379 + t433 * t380;
t591 = -t364 * qJD(4) - qJD(5) * t291 - t433 * t104 + t121 * t430;
t634 = t593 + t553 * pkin(11) + (-t493 + t538) * pkin(5);
t633 = pkin(11) * t552 - t594;
t602 = t87 * t89;
t632 = -pkin(11) * t548 + t592;
t631 = -pkin(5) * t313 + pkin(11) * t549 + t591;
t531 = qJDD(1) * t432;
t503 = t426 * t531;
t532 = qJD(1) * qJD(2);
t504 = t426 * t532;
t630 = t435 * t504 + t503;
t564 = t426 * t432;
t353 = -pkin(8) * t564 + t435 * t522;
t337 = qJD(2) * t353;
t629 = t645 - t646;
t454 = qJD(3) * t470;
t498 = t590 * qJDD(1);
t464 = t498 + qJDD(2);
t198 = qJD(3) * t494 - t431 * t464 + (-t454 - t630) * t434;
t530 = qJDD(1) * t435;
t396 = t426 * t530;
t490 = t432 * t504;
t329 = qJDD(3) - t396 + t490;
t163 = -t198 * t425 - t427 * t329;
t164 = -t198 * t427 + t329 * t425;
t463 = t163 * t433 + t430 * t164 - t253 * t536 + t255 * t535;
t507 = qJD(6) * t610;
t534 = qJD(6) * t429;
t64 = t430 * t163 - t433 * t164 + t253 * t535 + t255 * t536;
t22 = t157 * t507 - t159 * t534 + t429 * t463 + t610 * t64;
t628 = -t22 + t644;
t292 = -pkin(2) * t470 - t333;
t171 = t311 * pkin(3) - t313 * qJ(4) + t292;
t196 = t293 * t434 + t301 * t431;
t178 = -qJ(4) * t384 + t196;
t111 = t427 * t171 - t178 * t425;
t71 = pkin(4) * t311 - pkin(10) * t255 + t111;
t112 = t425 * t171 + t427 * t178;
t81 = -pkin(10) * t253 + t112;
t40 = -t430 * t81 + t433 * t71;
t34 = pkin(11) * t159 + t40;
t32 = pkin(5) * t303 + t34;
t41 = t430 * t71 + t433 * t81;
t35 = -pkin(11) * t157 + t41;
t539 = qJD(2) * t435;
t513 = t431 * t539;
t199 = t426 * (qJD(1) * (t432 * t537 + t513) + t431 * t531) + t431 * t454 - t434 * t464;
t187 = qJDD(5) + t199;
t466 = qJD(2) * t492;
t488 = pkin(1) * t498;
t524 = pkin(8) * t396 + t432 * t488 + t435 * t466;
t258 = -pkin(8) * t490 + t524;
t223 = pkin(9) * t464 + t258;
t461 = t483 * qJD(2);
t235 = (qJD(1) * t461 + qJDD(1) * t468) * t426;
t107 = t434 * t223 + t431 * t235 - t293 * t538 + t301 * t537;
t80 = qJ(4) * t329 - qJD(4) * t384 + t107;
t496 = pkin(8) * t630 + t432 * t466 - t435 * t488;
t224 = -pkin(2) * t464 + t496;
t84 = t199 * pkin(3) + t198 * qJ(4) - t313 * qJD(4) + t224;
t45 = -t425 * t80 + t427 * t84;
t30 = pkin(4) * t199 - pkin(10) * t164 + t45;
t46 = t425 * t84 + t427 * t80;
t39 = -pkin(10) * t163 + t46;
t9 = -qJD(5) * t41 + t433 * t30 - t430 * t39;
t6 = pkin(5) * t187 + pkin(11) * t64 + t9;
t462 = -t430 * t30 - t433 * t39 - t71 * t535 + t536 * t81;
t7 = -pkin(11) * t463 - t462;
t1 = t32 * t507 - t35 * t534 + t429 * t6 + t610 * t7;
t609 = sin(qJ(1));
t485 = t590 * t609;
t611 = cos(qJ(1));
t352 = -t432 * t485 + t435 * t611;
t520 = t426 * t609;
t286 = t352 * t434 + t431 * t520;
t351 = t432 * t611 + t435 * t485;
t422 = pkin(12) + qJ(5);
t417 = qJ(6) + t422;
t411 = sin(t417);
t412 = cos(t417);
t194 = t286 * t412 + t351 * t411;
t348 = t431 * t590 + t434 * t564;
t486 = t590 * t611;
t350 = t432 * t486 + t435 * t609;
t521 = t426 * t611;
t282 = t350 * t434 - t431 * t521;
t349 = t432 * t609 - t435 * t486;
t625 = t282 * t412 + t349 * t411;
t174 = pkin(3) * t384 + qJD(4) - t195;
t137 = pkin(4) * t253 + t174;
t76 = pkin(5) * t157 + t137;
t627 = t76 * t87 + g(1) * t194 + g(2) * t625 - g(3) * (-t348 * t412 + t411 * t563) - t1;
t626 = t282 * t411 - t349 * t412;
t415 = sin(t422);
t416 = cos(t422);
t624 = t282 * t415 - t349 * t416;
t623 = t282 * t416 + t349 * t415;
t421 = t426 ^ 2;
t621 = 0.2e1 * t421;
t323 = -pkin(2) * t590 - t353;
t347 = t431 * t564 - t434 * t590;
t207 = t347 * pkin(3) - t348 * qJ(4) + t323;
t546 = pkin(2) * t563 + pkin(9) * t564;
t325 = -pkin(1) * t426 - t546;
t227 = t434 * t324 + t431 * t325;
t208 = -qJ(4) * t563 + t227;
t133 = t425 * t207 + t427 * t208;
t276 = t348 * t425 + t427 * t563;
t113 = -pkin(10) * t276 + t133;
t132 = t427 * t207 - t208 * t425;
t277 = t348 * t427 - t425 * t563;
t99 = pkin(4) * t347 - pkin(10) * t277 + t132;
t51 = t433 * t113 + t430 * t99;
t620 = t195 * t384 + t107;
t244 = -t431 * t333 + t334 * t434;
t211 = -pkin(3) * t518 - t244;
t547 = pkin(4) * t479 + pkin(9) * t537 - t211;
t204 = -t286 * t415 + t351 * t416;
t619 = g(2) * t624 - g(3) * (-t348 * t415 - t416 * t563) - g(1) * t204;
t618 = (qJDD(2) + 0.2e1 * t498) * t426;
t338 = t545 * qJD(2);
t616 = -t329 * pkin(3) + qJDD(4);
t193 = -t286 * t411 + t351 * t412;
t527 = t610 * t35;
t13 = t429 * t32 + t527;
t2 = -qJD(6) * t13 - t429 * t7 + t610 * t6;
t615 = t76 * t89 - g(1) * t193 + g(2) * t626 - g(3) * (-t348 * t411 - t412 * t563) + t2;
t23 = -qJD(6) * t89 - t429 * t64 + t610 * t463;
t614 = -t23 - t643;
t612 = t311 ^ 2;
t436 = qJD(1) ^ 2;
t607 = pkin(4) * t425;
t605 = g(3) * t426;
t604 = t163 * pkin(4);
t368 = pkin(5) * t415 + t607;
t601 = pkin(9) + t368;
t176 = t433 * t268 - t280 * t430;
t331 = t363 * t431;
t145 = -pkin(5) * t434 + pkin(11) * t331 + t176;
t330 = t364 * t431;
t151 = -pkin(11) * t330 + t177;
t82 = t145 * t610 - t429 * t151;
t599 = qJD(6) * t82 + t429 * t634 - t633 * t610;
t83 = t429 * t145 + t151 * t610;
t598 = -qJD(6) * t83 + t633 * t429 + t610 * t634;
t290 = -t433 * t379 - t380 * t430;
t247 = -pkin(11) * t364 + t290;
t248 = -pkin(11) * t363 + t291;
t146 = t247 * t610 - t429 * t248;
t597 = qJD(6) * t146 + t429 * t631 + t610 * t632;
t147 = t429 * t247 + t248 * t610;
t596 = -qJD(6) * t147 - t429 * t632 + t610 * t631;
t595 = t429 * t35;
t588 = t159 * t157;
t587 = t163 * t427;
t586 = t164 * t425;
t584 = t196 * t384;
t583 = t199 * t425;
t582 = t199 * t427;
t581 = t311 * t313;
t580 = t311 * t384;
t577 = t313 * t384;
t458 = t384 * t431;
t572 = t411 * t434;
t571 = t412 * t434;
t570 = t415 * t434;
t569 = t416 * t434;
t568 = t421 * t436;
t558 = t330 * t507 - t331 * t534 + t429 * t552 + t553 * t610;
t234 = -t429 * t330 - t331 * t610;
t557 = qJD(6) * t234 - t429 * t553 + t552 * t610;
t267 = -t429 * t363 + t364 * t610;
t556 = -qJD(6) * t267 + t429 * t549 - t548 * t610;
t555 = t363 * t507 + t364 * t534 + t429 * t548 + t549 * t610;
t335 = t426 * t461;
t148 = -t324 * t538 + t325 * t537 + t431 * t335 + t434 * t337;
t540 = qJD(2) * t432;
t138 = (qJ(4) * t540 - qJD(4) * t435) * t426 + t148;
t278 = qJD(3) * t348 + t426 * t513;
t514 = t426 * t539;
t279 = -qJD(3) * t347 + t434 * t514;
t144 = t278 * pkin(3) - t279 * qJ(4) - t348 * qJD(4) + t338;
t73 = t427 * t138 + t425 * t144;
t554 = pkin(5) * t552 + t547;
t550 = -t427 * t528 - t642;
t366 = pkin(4) * t567 + t431 * pkin(9);
t544 = t611 * pkin(1) + pkin(8) * t520;
t423 = t432 ^ 2;
t424 = t435 ^ 2;
t543 = t423 - t424;
t533 = -qJD(4) + t174;
t526 = t435 * t568;
t525 = t352 * pkin(2) + t544;
t413 = pkin(4) * t427 + pkin(3);
t523 = pkin(9) + t607;
t515 = t426 * t540;
t508 = g(3) * t546;
t506 = pkin(1) * t621;
t505 = t435 * t532;
t165 = -pkin(4) * t579 + t196;
t501 = pkin(5) * t548 - t165;
t50 = -t113 * t430 + t433 * t99;
t72 = -t138 * t425 + t427 * t144;
t226 = -t431 * t324 + t325 * t434;
t497 = t432 * t526;
t108 = -t431 * t223 + t434 * t235 - t293 * t537 - t301 * t538;
t491 = t432 * t505;
t487 = -pkin(1) * t609 + pkin(8) * t521;
t484 = t426 * t436 * t590;
t281 = t350 * t431 + t434 * t521;
t285 = t352 * t431 - t434 * t520;
t482 = -g(1) * t281 + g(2) * t285;
t481 = -g(1) * t349 + g(2) * t351;
t480 = g(1) * t352 + g(2) * t350;
t209 = pkin(3) * t563 - t226;
t475 = t199 * t347 + t278 * t311;
t474 = -t253 * t427 - t255 * t425;
t180 = -t276 * t430 + t277 * t433;
t365 = pkin(5) * t416 + t413;
t420 = -pkin(11) - t600;
t473 = t365 * t434 - t420 * t431;
t472 = t413 * t434 + t431 * t600;
t469 = 0.2e1 * t499 + qJD(2);
t467 = pkin(9) * t351 + t525;
t465 = -t350 * pkin(2) + t487;
t149 = -t324 * t537 - t325 * t538 + t335 * t434 - t431 * t337;
t42 = pkin(5) * t347 - pkin(11) * t180 + t50;
t179 = t433 * t276 + t277 * t430;
t44 = -pkin(11) * t179 + t51;
t18 = t42 * t610 - t429 * t44;
t19 = t429 * t42 + t44 * t610;
t120 = -t429 * t179 + t180 * t610;
t240 = t279 * t427 + t425 * t515;
t57 = pkin(4) * t278 - pkin(10) * t240 + t72;
t239 = t279 * t425 - t427 * t515;
t66 = -pkin(10) * t239 + t73;
t16 = -t113 * t536 + t430 * t57 + t433 * t66 + t99 * t535;
t457 = g(1) * t611 + g(2) * t609;
t456 = g(1) * t285 + g(2) * t281 + g(3) * t347;
t455 = -g(1) * t286 - g(2) * t282 - g(3) * t348;
t162 = pkin(4) * t276 + t209;
t85 = -t108 + t616;
t452 = -t349 * pkin(9) + t465;
t451 = t456 - t85;
t450 = -g(1) * t351 - g(2) * t349 + g(3) * t563;
t449 = -g(3) * t564 - t480;
t447 = -pkin(9) * t329 - t292 * t384;
t446 = t450 * t431;
t63 = t85 + t604;
t143 = -pkin(3) * t515 - t149;
t444 = -t199 * t434 - t311 * t458;
t17 = -qJD(5) * t51 - t430 * t66 + t433 * t57;
t443 = t456 + t108;
t441 = pkin(9) * qJD(3) * t384 - t224 - t450;
t114 = pkin(4) * t239 + t143;
t440 = -t443 + t616;
t31 = pkin(5) * t463 + t63;
t341 = t351 * pkin(2);
t339 = t349 * pkin(2);
t336 = t545 * qJD(1);
t320 = pkin(5) * t363 - t413;
t308 = -pkin(9) * t565 + t358;
t271 = pkin(5) * t330 + t366;
t266 = t363 * t610 + t364 * t429;
t233 = t330 * t610 - t331 * t429;
t205 = t286 * t416 + t351 * t415;
t184 = qJDD(6) + t187;
t119 = t179 * t610 + t180 * t429;
t109 = pkin(5) * t179 + t162;
t103 = qJD(5) * t180 + t433 * t239 + t240 * t430;
t102 = t430 * t239 - t433 * t240 + t276 * t535 + t277 * t536;
t52 = pkin(5) * t103 + t114;
t38 = qJD(6) * t120 - t429 * t102 + t103 * t610;
t37 = t102 * t610 + t429 * t103 + t179 * t507 + t180 * t534;
t15 = t34 * t610 - t595;
t14 = -t429 * t34 - t527;
t12 = t32 * t610 - t595;
t11 = -pkin(11) * t103 + t16;
t10 = pkin(5) * t278 + pkin(11) * t102 + t17;
t4 = -qJD(6) * t19 + t10 * t610 - t429 * t11;
t3 = qJD(6) * t18 + t429 * t10 + t11 * t610;
t5 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t609 - g(2) * t611, t457, 0, 0 (qJDD(1) * t423 + 0.2e1 * t491) * t421 (t432 * t530 - t532 * t543) * t621, t432 * t618 + t469 * t514 (qJDD(1) * t424 - 0.2e1 * t491) * t421, t435 * t618 - t469 * t515, t464 * t590, -t338 * t470 + t353 * t464 - t496 * t590 + g(1) * t350 - g(2) * t352 + (-t432 * t532 + t530) * t506, -t337 * t470 - t545 * t464 - t258 * t590 + (-t505 - t531) * t506 + t481 ((-t333 * qJD(2) + qJDD(1) * t545 + t258) * t435 + (-qJD(2) * t336 - qJDD(1) * t353 + t496) * t432 - t457) * t426, t421 * qJDD(1) * pkin(1) ^ 2 - g(1) * t487 - g(2) * t544 + t258 * t545 - t333 * t338 + t336 * t337 - t353 * t496, -t198 * t348 + t279 * t313, t198 * t347 - t199 * t348 - t278 * t313 - t279 * t311, -t279 * t384 + t329 * t348 + (t198 * t435 + t313 * t540) * t426, t475, t278 * t384 - t329 * t347 + (t199 * t435 - t311 * t540) * t426 (-t329 * t435 - t384 * t540) * t426, g(1) * t282 - g(2) * t286 - t149 * t384 + t199 * t323 + t224 * t347 + t226 * t329 + t278 * t292 + t311 * t338 + (-t108 * t435 + t195 * t540) * t426, t148 * t384 - t198 * t323 + t224 * t348 - t227 * t329 + t279 * t292 + t313 * t338 + (t107 * t435 - t196 * t540) * t426 + t482, -t107 * t347 - t108 * t348 - t148 * t311 - t149 * t313 - t195 * t279 - t196 * t278 + t198 * t226 - t199 * t227 - t481, -g(1) * t452 - g(2) * t467 + t107 * t227 + t108 * t226 + t196 * t148 + t195 * t149 + t224 * t323 + t292 * t338, t164 * t277 + t240 * t255, -t163 * t277 - t164 * t276 - t239 * t255 - t240 * t253, t164 * t347 + t199 * t277 + t240 * t311 + t255 * t278, t163 * t276 + t239 * t253, -t163 * t347 - t199 * t276 - t239 * t311 - t253 * t278, t475, t72 * t311 + t132 * t199 + t45 * t347 + t111 * t278 + t143 * t253 + t209 * t163 + t85 * t276 + t174 * t239 - g(1) * (-t282 * t427 - t349 * t425) - g(2) * (t286 * t427 + t351 * t425) -t73 * t311 - t133 * t199 - t46 * t347 - t112 * t278 + t143 * t255 + t209 * t164 + t85 * t277 + t174 * t240 - g(1) * (t282 * t425 - t349 * t427) - g(2) * (-t286 * t425 + t351 * t427) -t111 * t240 - t112 * t239 - t132 * t164 - t133 * t163 - t253 * t73 - t255 * t72 - t276 * t46 - t277 * t45 - t482, t46 * t133 + t112 * t73 + t45 * t132 + t111 * t72 + t85 * t209 + t174 * t143 - g(1) * (-pkin(3) * t282 - qJ(4) * t281 + t452) - g(2) * (pkin(3) * t286 + qJ(4) * t285 + t467) t102 * t159 - t180 * t64, t102 * t157 + t103 * t159 + t64 * t179 - t180 * t463, -t102 * t303 - t159 * t278 + t180 * t187 - t347 * t64, t157 * t103 + t179 * t463, -t103 * t303 - t157 * t278 - t179 * t187 - t347 * t463, t187 * t347 + t278 * t303, g(1) * t623 - g(2) * t205 + t137 * t103 + t114 * t157 + t162 * t463 + t17 * t303 + t63 * t179 + t50 * t187 + t40 * t278 + t9 * t347, -g(1) * t624 - g(2) * t204 - t137 * t102 - t114 * t159 - t16 * t303 - t162 * t64 + t63 * t180 - t51 * t187 - t41 * t278 + t462 * t347, t40 * t102 - t41 * t103 - t16 * t157 + t159 * t17 + t179 * t462 - t9 * t180 - t463 * t51 + t50 * t64 - t482, -t462 * t51 + t41 * t16 + t9 * t50 + t40 * t17 + t63 * t162 + t137 * t114 - g(1) * (-t281 * t600 - t282 * t413 - t349 * t523 + t465) - g(2) * (t285 * t600 + t286 * t413 + t351 * t523 + t525) -t120 * t22 + t37 * t89, t119 * t22 - t120 * t23 + t37 * t87 + t38 * t89, t120 * t184 - t22 * t347 - t278 * t89 - t297 * t37, t119 * t23 + t38 * t87, -t119 * t184 - t23 * t347 - t278 * t87 - t297 * t38, t184 * t347 + t278 * t297, g(1) * t625 - g(2) * t194 + t109 * t23 + t31 * t119 + t12 * t278 + t18 * t184 + t2 * t347 + t4 * t297 + t76 * t38 + t52 * t87, -g(1) * t626 - g(2) * t193 - t1 * t347 - t109 * t22 + t31 * t120 - t13 * t278 - t19 * t184 - t3 * t297 - t76 * t37 - t52 * t89, -t1 * t119 + t12 * t37 - t120 * t2 - t13 * t38 + t18 * t22 - t19 * t23 - t3 * t87 + t4 * t89 - t482, t1 * t19 + t13 * t3 + t2 * t18 + t12 * t4 + t31 * t109 + t76 * t52 - g(1) * (t281 * t420 - t282 * t365 - t349 * t601 + t465) - g(2) * (-t285 * t420 + t286 * t365 + t351 * t601 + t525); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t497, t543 * t568, -t435 * t484 + t503, t497, t432 * t484 + t396, t464, pkin(1) * t432 * t568 + t336 * t470 - t450 - t496, pkin(1) * t526 + t333 * t470 + (pkin(8) * t532 + g(3)) * t564 + t480 - t524, 0, 0, -t198 * t431 - t434 * t577 (-t198 + t580) * t434 + (-t199 + t577) * t431, -t384 * t537 + t329 * t431 + (-t313 * t432 + t384 * t559) * t542, t444, t384 * t538 + t329 * t434 + (t311 * t432 - t435 * t458) * t542, t384 * t518, -pkin(2) * t199 - t195 * t518 + t244 * t384 - t311 * t336 + t431 * t447 + t434 * t441, pkin(2) * t198 + t196 * t518 - t245 * t384 - t313 * t336 - t431 * t441 + t434 * t447, t244 * t313 + t245 * t311 + ((qJD(3) * t313 - t199) * pkin(9) + t620) * t434 + (-t108 + t584 + (qJD(3) * t311 - t198) * pkin(9)) * t431 + t449, -t224 * pkin(2) - t196 * t245 - t195 * t244 - t292 * t336 + g(1) * t341 + g(2) * t339 - t508 + (t107 * t434 - t108 * t431 + (-t195 * t434 - t196 * t431) * qJD(3) - t480) * pkin(9), t164 * t562 - t255 * t641, t253 * t296 + t255 * t295 + (-t586 - t587) * t431 + t474 * t537, -t164 * t434 - t641 * t311 + (-t255 * t384 + t582) * t431, t163 * t567 + t253 * t479, t163 * t434 - t479 * t311 + (t253 * t384 - t583) * t431, t444, -t174 * t295 + t308 * t199 - t211 * t253 + t551 * t311 + t449 * t425 + (-t45 + (pkin(9) * t253 + t174 * t425) * qJD(3) - t450 * t427) * t434 + (pkin(9) * t163 - t111 * t384 + t85 * t425) * t431, -t174 * t296 - t309 * t199 - t211 * t255 - t550 * t311 + t449 * t427 + (t46 + (pkin(9) * t255 + t174 * t427) * qJD(3) + t450 * t425) * t434 + (pkin(9) * t164 + t112 * t384 + t85 * t427) * t431, t111 * t296 + t112 * t295 - t163 * t309 - t164 * t308 - t551 * t255 - t550 * t253 + (-t111 * t427 - t112 * t425) * t537 + (-t425 * t46 - t427 * t45 - t450) * t431, t46 * t309 + t45 * t308 - t174 * t211 - g(1) * (-t351 * t477 - t341) - g(2) * (-t349 * t477 - t339) - g(3) * (t477 * t563 + t546) + t550 * t112 + t551 * t111 + (t174 * t537 + t431 * t85 - t480) * pkin(9), t159 * t553 + t331 * t64, t157 * t553 + t159 * t552 + t64 * t330 + t331 * t463, t159 * t458 - t187 * t331 - t303 * t553 + t434 * t64, t157 * t552 + t330 * t463, t157 * t458 - t330 * t187 - t303 * t552 + t434 * t463, -t187 * t434 - t303 * t458, t176 * t187 - t9 * t434 + t40 * t538 + t366 * t463 + t63 * t330 - g(1) * (-t351 * t569 + t352 * t415) - g(2) * (-t349 * t569 + t350 * t415) + t593 * t303 + t547 * t157 + t552 * t137 + (-t40 * t516 - g(3) * (t415 * t432 + t416 * t559)) * t426, -t177 * t187 - t462 * t434 - t41 * t538 - t366 * t64 - t63 * t331 - g(1) * (t351 * t570 + t352 * t416) - g(2) * (t349 * t570 + t350 * t416) - t594 * t303 - t547 * t159 - t553 * t137 + (t41 * t516 - g(3) * (-t415 * t559 + t416 * t432)) * t426, -t157 * t594 + t159 * t593 + t176 * t64 - t177 * t463 + t330 * t462 + t9 * t331 + t40 * t553 - t41 * t552 - t446, -t462 * t177 + t9 * t176 + t63 * t366 - g(1) * (-t351 * t472 + t352 * t523 - t341) - g(2) * (-t349 * t472 + t350 * t523 - t339) - t508 - (pkin(4) * t566 + t435 * t472) * t605 + t594 * t41 + t593 * t40 + t547 * t137, -t22 * t234 + t558 * t89, t22 * t233 - t23 * t234 + t557 * t89 + t558 * t87, t184 * t234 + t22 * t434 - t297 * t558 + t458 * t89, t23 * t233 + t557 * t87, -t184 * t233 + t23 * t434 - t297 * t557 + t458 * t87, -t184 * t434 - t297 * t458, t82 * t184 - t2 * t434 + t12 * t538 + t271 * t23 + t31 * t233 - g(1) * (-t351 * t571 + t352 * t411) - g(2) * (-t349 * t571 + t350 * t411) + t554 * t87 + t557 * t76 + t598 * t297 + (-t12 * t516 - g(3) * (t411 * t432 + t412 * t559)) * t426, -t83 * t184 + t1 * t434 - t13 * t538 - t271 * t22 + t31 * t234 - g(1) * (t351 * t572 + t352 * t412) - g(2) * (t349 * t572 + t350 * t412) - t554 * t89 - t558 * t76 - t599 * t297 + (t13 * t516 - g(3) * (-t411 * t559 + t412 * t432)) * t426, -t1 * t233 + t12 * t558 - t13 * t557 - t2 * t234 + t22 * t82 - t23 * t83 + t598 * t89 - t599 * t87 - t446, t1 * t83 + t2 * t82 + t31 * t271 - g(1) * (-t351 * t473 + t352 * t601 - t341) - g(2) * (-t349 * t473 + t350 * t601 - t339) - t508 + t554 * t76 - (t368 * t432 + t435 * t473) * t605 + t599 * t13 + t598 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t581, t313 ^ 2 - t612, -t198 - t580, -t581, -t199 - t577, t329, -t292 * t313 + t443 - t584, t292 * t311 - t455 - t620, 0, 0, t255 * t578 + t586, -t163 * t425 + t164 * t427 + t311 * t474, -t255 * t313 + t427 * t612 + t583, t253 * t579 - t587, t253 * t313 - t425 * t612 + t582, -t581, -qJ(4) * t583 - pkin(3) * t163 - t111 * t313 - t196 * t253 + (t425 * t533 - t134) * t311 + t451 * t427, -qJ(4) * t582 - pkin(3) * t164 + t112 * t313 - t196 * t255 + (t427 * t533 + t135) * t311 - t451 * t425, t134 * t255 + t135 * t253 + (-qJ(4) * t163 - qJD(4) * t253 - t111 * t311 + t46) * t427 + (qJ(4) * t164 + qJD(4) * t255 - t112 * t311 - t45) * t425 + t455, -t111 * t134 - t112 * t135 - t174 * t196 + (-t111 * t425 + t112 * t427) * qJD(4) + t451 * pkin(3) + (-t425 * t45 + t427 * t46 + t455) * qJ(4), t159 * t549 - t364 * t64, t157 * t549 + t159 * t548 + t64 * t363 - t364 * t463, t159 * t313 + t187 * t364 - t303 * t549, t157 * t548 + t363 * t463, t157 * t313 - t187 * t363 - t303 * t548, -t303 * t313, t137 * t548 - t165 * t157 + t290 * t187 + t303 * t591 - t40 * t313 + t63 * t363 - t413 * t463 + t416 * t456, -t137 * t549 + t159 * t165 - t187 * t291 - t303 * t592 + t313 * t41 + t364 * t63 + t413 * t64 - t415 * t456, -t157 * t592 + t159 * t591 + t290 * t64 - t291 * t463 + t363 * t462 - t9 * t364 + t40 * t549 - t41 * t548 + t455, -t462 * t291 + t9 * t290 - t63 * t413 - t137 * t165 - g(1) * (-t285 * t413 + t286 * t600) - g(2) * (-t281 * t413 + t282 * t600) - g(3) * (-t347 * t413 + t348 * t600) + t592 * t41 + t591 * t40, -t22 * t267 + t555 * t89, t22 * t266 - t23 * t267 + t555 * t87 - t556 * t89, t184 * t267 - t297 * t555 + t313 * t89, t23 * t266 - t556 * t87, -t184 * t266 + t297 * t556 + t313 * t87, -t297 * t313, -t12 * t313 + t146 * t184 + t23 * t320 + t266 * t31 + t297 * t596 + t412 * t456 + t501 * t87 - t556 * t76, t13 * t313 - t147 * t184 - t22 * t320 + t267 * t31 - t297 * t597 - t411 * t456 - t501 * t89 - t555 * t76, -t1 * t266 + t12 * t555 + t13 * t556 + t146 * t22 - t147 * t23 - t2 * t267 + t596 * t89 - t597 * t87 + t455, t1 * t147 + t2 * t146 + t31 * t320 - g(1) * (-t285 * t365 - t286 * t420) - g(2) * (-t281 * t365 - t282 * t420) - g(3) * (-t347 * t365 - t348 * t420) + t501 * t76 + t597 * t13 + t596 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255 * t311 + t163, -t253 * t311 + t164, -t253 ^ 2 - t255 ^ 2, t111 * t255 + t112 * t253 + t440, 0, 0, 0, 0, 0, 0, t463 - t635, -t64 - t636, -t637 - t638, t157 * t41 - t159 * t40 + t440 + t604, 0, 0, 0, 0, 0, 0, t23 - t643, -t22 - t644, -t645 - t646, -t12 * t89 + t13 * t87 + t31 - t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, t637 - t638, -t64 + t636, t588, -t463 - t635, t187, t137 * t159 + t41 * t303 + t619 + t9, t40 * t303 + t137 * t157 + g(1) * t205 + g(2) * t623 - g(3) * (-t348 * t416 + t415 * t563) + t462, 0, 0, -t602, t629, t628, t602, t614, t184, -t14 * t297 + (t159 * t87 + t184 * t610 - t297 * t534) * pkin(5) + t615, t15 * t297 + (-t159 * t89 - t184 * t429 - t297 * t507) * pkin(5) + t627, -t12 * t87 - t13 * t89 - t14 * t89 + t15 * t87 + (t610 * t22 - t23 * t429 + (-t429 * t89 - t610 * t87) * qJD(6)) * pkin(5), -t12 * t14 - t13 * t15 + (t1 * t429 + t2 * t610 + t76 * t159 + (-t12 * t429 + t13 * t610) * qJD(6) + t619) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t602, t629, t628, t602, t614, t184, t13 * t297 + t615, t12 * t297 + t627, 0, 0;];
tau_reg  = t5;