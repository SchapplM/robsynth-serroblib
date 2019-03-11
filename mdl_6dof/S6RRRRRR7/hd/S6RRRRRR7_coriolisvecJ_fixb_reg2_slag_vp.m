% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:22
% EndTime: 2019-03-10 04:42:12
% DurationCPUTime: 27.15s
% Computational Cost: add. (48209->853), mult. (120749->1179), div. (0->0), fcn. (96850->12), ass. (0->336)
t407 = sin(qJ(2));
t528 = cos(pkin(6));
t465 = pkin(1) * t528;
t390 = t407 * t465;
t406 = sin(qJ(3));
t410 = cos(qJ(3));
t433 = pkin(3) * t406 - pkin(10) * t410;
t402 = sin(pkin(6));
t411 = cos(qJ(2));
t510 = t402 * t411;
t590 = (t390 + (pkin(8) + t433) * t510) * qJD(1) - t433 * qJD(3);
t482 = qJD(1) * t411;
t384 = t402 * t482;
t409 = cos(qJ(4));
t483 = qJD(1) * t407;
t461 = t402 * t483;
t405 = sin(qJ(4));
t505 = t405 * t410;
t298 = t384 * t505 - t409 * t461;
t477 = qJD(3) * t410;
t587 = -t405 * t477 + t298;
t475 = qJD(4) * t405;
t484 = qJD(1) * t402;
t501 = t410 * t411;
t299 = (t405 * t407 + t409 * t501) * t484;
t586 = -t409 * t477 + t299;
t589 = -t406 * t475 - t586;
t474 = qJD(4) * t409;
t581 = t406 * t474 - t587;
t441 = t411 * t465;
t330 = -pkin(8) * t461 + qJD(1) * t441;
t423 = (pkin(2) * t407 - pkin(9) * t411) * t402;
t331 = qJD(1) * t423;
t248 = t410 * t330 + t406 * t331;
t232 = pkin(10) * t461 + t248;
t478 = qJD(3) * t406;
t540 = pkin(9) * t405;
t588 = -t232 * t405 + t590 * t409 - t478 * t540;
t374 = -pkin(3) * t410 - pkin(10) * t406 - pkin(2);
t493 = t409 * t232 - t374 * t474 - (-t409 * t478 - t410 * t475) * pkin(9) + t590 * t405;
t486 = pkin(8) * t510 + t390;
t328 = t528 * pkin(9) + t486;
t289 = qJD(2) * pkin(9) + qJD(1) * t328;
t329 = (-pkin(2) * t411 - pkin(9) * t407 - pkin(1)) * t402;
t304 = qJD(1) * t329;
t219 = -t406 * t289 + t410 * t304;
t449 = t528 * qJD(1);
t426 = t449 + qJD(2);
t311 = t406 * t461 - t410 * t426;
t313 = t406 * t426 + t410 * t461;
t237 = pkin(3) * t313 + pkin(10) * t311;
t141 = -t219 * t405 + t409 * t237;
t542 = -pkin(11) - pkin(10);
t464 = qJD(4) * t542;
t585 = -pkin(4) * t313 - t141 + (-pkin(11) * t311 + t464) * t409;
t142 = t409 * t219 + t405 * t237;
t514 = t311 * t405;
t584 = pkin(11) * t514 - t405 * t464 + t142;
t502 = t409 * t410;
t392 = pkin(9) * t502;
t439 = t406 * t384;
t583 = pkin(4) * t439 - pkin(11) * t299 - (pkin(4) * t406 - pkin(11) * t502) * qJD(3) - (-t392 + (pkin(11) * t406 - t374) * t405) * qJD(4) + t588;
t582 = t581 * pkin(11) + t493;
t431 = t384 - qJD(3);
t446 = t409 * t431;
t255 = t313 * t405 + t446;
t257 = t409 * t313 - t405 * t431;
t404 = sin(qJ(5));
t408 = cos(qJ(5));
t175 = t255 * t408 + t404 * t257;
t403 = sin(qJ(6));
t429 = -t255 * t404 + t408 * t257;
t541 = cos(qJ(6));
t100 = t541 * t175 + t403 * t429;
t562 = -t403 * t175 + t429 * t541;
t527 = t100 * t562;
t358 = t404 * t409 + t405 * t408;
t469 = qJD(4) + qJD(5);
t278 = t469 * t358;
t492 = t278 * t406 - t587 * t404 + t586 * t408;
t473 = qJD(5) * t404;
t504 = t406 * t409;
t506 = t405 * t406;
t491 = -t473 * t506 + (t469 * t504 - t587) * t408 + t589 * t404;
t508 = t404 * t405;
t357 = -t408 * t409 + t508;
t472 = qJD(5) * t408;
t490 = t357 * t311 - t408 * t474 - t409 * t472 + t469 * t508;
t489 = t358 * t311 + t278;
t568 = -t100 ^ 2 + t562 ^ 2;
t305 = qJD(4) + t311;
t294 = qJD(5) + t305;
t288 = -pkin(2) * t426 - t330;
t192 = t311 * pkin(3) - t313 * pkin(10) + t288;
t220 = t410 * t289 + t406 * t304;
t198 = -pkin(10) * t431 + t220;
t117 = t409 * t192 - t198 * t405;
t94 = -pkin(11) * t257 + t117;
t85 = pkin(4) * t305 + t94;
t118 = t192 * t405 + t198 * t409;
t95 = -pkin(11) * t255 + t118;
t91 = t404 * t95;
t47 = t408 * t85 - t91;
t557 = pkin(12) * t429;
t37 = t47 - t557;
t32 = pkin(5) * t294 + t37;
t93 = t408 * t95;
t48 = t404 * t85 + t93;
t576 = pkin(12) * t175;
t38 = t48 - t576;
t533 = t403 * t38;
t14 = t541 * t32 - t533;
t467 = t541 * t38;
t15 = t403 * t32 + t467;
t580 = -t14 * t100 + t15 * t562;
t470 = qJD(1) * qJD(2);
t452 = t402 * t470;
t437 = t411 * t452;
t577 = qJD(3) * t311;
t260 = -t410 * t437 + t577;
t438 = t407 * t452;
t159 = qJD(4) * t446 + t409 * t260 + t313 * t475 - t405 * t438;
t476 = qJD(4) * t257;
t160 = -t405 * t260 - t409 * t438 + t476;
t424 = -t404 * t159 + t160 * t408 - t255 * t473 + t257 * t472;
t454 = qJD(6) * t541;
t471 = qJD(6) * t403;
t71 = t408 * t159 + t404 * t160 + t255 * t472 + t257 * t473;
t24 = t175 * t454 + t403 * t424 + t429 * t471 + t541 * t71;
t290 = qJD(6) + t294;
t566 = t100 * t290 - t24;
t427 = t406 * t437;
t479 = qJD(3) * t313;
t261 = t427 + t479;
t332 = qJD(2) * t423;
t323 = qJD(1) * t332;
t511 = t402 * t407;
t349 = -pkin(8) * t511 + t441;
t334 = t349 * qJD(2);
t324 = qJD(1) * t334;
t418 = t289 * t478 - t304 * t477 - t406 * t323 - t410 * t324;
t135 = pkin(10) * t438 - t418;
t335 = t486 * qJD(2);
t325 = qJD(1) * t335;
t155 = t261 * pkin(3) + t260 * pkin(10) + t325;
t58 = -qJD(4) * t118 - t405 * t135 + t409 * t155;
t35 = t261 * pkin(4) + t159 * pkin(11) + t58;
t57 = t409 * t135 + t405 * t155 + t192 * t474 - t198 * t475;
t43 = -pkin(11) * t160 + t57;
t9 = -qJD(5) * t48 + t408 * t35 - t404 * t43;
t6 = t261 * pkin(5) + t71 * pkin(12) + t9;
t450 = -t404 * t35 - t408 * t43 - t85 * t472 + t95 * t473;
t7 = -pkin(12) * t424 - t450;
t414 = -t32 * t454 + t38 * t471 - t403 * t6 - t541 * t7;
t197 = pkin(3) * t431 - t219;
t143 = t255 * pkin(4) + t197;
t89 = t175 * pkin(5) + t143;
t565 = t89 * t100 + t414;
t355 = t409 * t374;
t273 = -pkin(11) * t504 + t355 + (-pkin(4) - t540) * t410;
t322 = t405 * t374 + t392;
t286 = -pkin(11) * t506 + t322;
t535 = -t273 * t472 + t286 * t473 + t583 * t404 + t582 * t408;
t201 = t404 * t273 + t408 * t286;
t534 = -qJD(5) * t201 + t582 * t404 - t583 * t408;
t378 = t542 * t405;
t379 = t542 * t409;
t530 = t378 * t472 + t379 * t473 + t585 * t404 - t584 * t408;
t293 = t404 * t378 - t408 * t379;
t529 = -qJD(5) * t293 + t584 * t404 + t585 * t408;
t2 = -qJD(6) * t15 - t403 * t7 + t541 * t6;
t547 = -t562 * t89 + t2;
t25 = qJD(6) * t562 - t403 * t71 + t541 * t424;
t544 = t290 * t562 - t25;
t575 = t534 + t492 * pkin(12) + (-t439 + t478) * pkin(5);
t574 = t491 * pkin(12) + t535;
t573 = -t489 * pkin(12) + t530;
t572 = pkin(5) * t313 - t490 * pkin(12) - t529;
t520 = t175 * t429;
t443 = qJD(3) * t431;
t571 = -pkin(9) * t443 + t325;
t569 = t431 * t219 - t418;
t567 = -t175 ^ 2 + t429 ^ 2;
t564 = t175 * t294 - t71;
t563 = t143 * t175 + t450;
t399 = t402 ^ 2;
t453 = t399 * t470;
t559 = -0.2e1 * t453;
t556 = -t117 * t305 + t57;
t555 = -t118 * t305 - t58;
t554 = t311 * t431;
t553 = t313 * t431;
t551 = t407 * t411;
t548 = t431 * t406;
t327 = -pkin(2) * t528 - t349;
t341 = t406 * t511 - t410 * t528;
t342 = t406 * t528 + t410 * t511;
t227 = t341 * pkin(3) - t342 * pkin(10) + t327;
t239 = t410 * t328 + t406 * t329;
t229 = -pkin(10) * t510 + t239;
t139 = t409 * t227 - t229 * t405;
t466 = t405 * t510;
t283 = t342 * t409 - t466;
t113 = pkin(4) * t341 - pkin(11) * t283 + t139;
t140 = t405 * t227 + t409 * t229;
t282 = t342 * t405 + t409 * t510;
t119 = -pkin(11) * t282 + t140;
t63 = t404 * t113 + t408 * t119;
t436 = -t220 + (t475 + t514) * pkin(4);
t246 = -t406 * t330 + t331 * t410;
t231 = -pkin(3) * t461 - t246;
t488 = t581 * pkin(4) + pkin(9) * t477 - t231;
t546 = -t143 * t429 + t9;
t145 = -t289 * t477 - t304 * t478 + t410 * t323 - t406 * t324;
t545 = t431 * t220 - t145;
t543 = t294 * t429 - t424;
t412 = qJD(1) ^ 2;
t292 = t408 * t378 + t379 * t404;
t252 = -pkin(12) * t358 + t292;
t253 = -pkin(12) * t357 + t293;
t171 = t403 * t252 + t541 * t253;
t539 = t171 * qJD(6) + t573 * t403 + t572 * t541;
t170 = t541 * t252 - t403 * t253;
t538 = -t170 * qJD(6) + t572 * t403 - t573 * t541;
t200 = t408 * t273 - t286 * t404;
t337 = t357 * t406;
t165 = -pkin(5) * t410 + pkin(12) * t337 + t200;
t336 = t358 * t406;
t172 = -pkin(12) * t336 + t201;
t96 = t541 * t165 - t403 * t172;
t537 = t96 * qJD(6) + t575 * t403 - t574 * t541;
t97 = t403 * t165 + t541 * t172;
t536 = -t97 * qJD(6) + t574 * t403 + t575 * t541;
t52 = t408 * t94 - t91;
t395 = pkin(4) * t408 + pkin(5);
t51 = -t404 * t94 - t93;
t40 = t51 + t576;
t41 = t52 - t557;
t462 = t541 * t404;
t532 = t541 * t40 - t403 * t41 + t395 * t471 - (-t404 * t454 + (-t403 * t408 - t462) * qJD(5)) * pkin(4);
t509 = t403 * t404;
t531 = t403 * t40 + t541 * t41 - t395 * t454 - (-t404 * t471 + (t541 * t408 - t509) * qJD(5)) * pkin(4);
t136 = -pkin(3) * t438 - t145;
t524 = t136 * t405;
t523 = t136 * t409;
t522 = t159 * t405;
t521 = t160 * t409;
t519 = t255 * t305;
t518 = t255 * t405;
t517 = t257 * t255;
t516 = t257 * t305;
t230 = t261 * t341;
t515 = t261 * t410;
t513 = t313 * t311;
t512 = t399 * t412;
t507 = t405 * t261;
t503 = t409 * t261;
t500 = t336 * t454 - t337 * t471 + t491 * t403 + t492 * t541;
t247 = -t403 * t336 - t541 * t337;
t499 = t247 * qJD(6) - t492 * t403 + t491 * t541;
t498 = t489 * pkin(5) + t436;
t497 = t491 * pkin(5) + t488;
t271 = -t403 * t357 + t541 * t358;
t496 = -t271 * qJD(6) + t490 * t403 - t489 * t541;
t495 = t357 * t454 + t358 * t471 + t489 * t403 + t490 * t541;
t494 = qJD(4) * t322 + t588;
t366 = pkin(4) * t506 + t406 * pkin(9);
t485 = t407 ^ 2 - t411 ^ 2;
t481 = qJD(2) * t407;
t480 = qJD(2) * t410;
t396 = -pkin(4) * t409 - pkin(3);
t460 = t402 * t481;
t459 = qJD(2) * t510;
t62 = t408 * t113 - t119 * t404;
t238 = -t406 * t328 + t329 * t410;
t447 = t305 * t409;
t445 = t411 * t431;
t444 = t431 * t402;
t442 = t512 * t551;
t435 = t402 * t412 * t528;
t434 = pkin(1) * t559;
t228 = pkin(3) * t510 - t238;
t430 = -t117 * t409 - t118 * t405;
t204 = -t282 * t404 + t283 * t408;
t428 = t453 * t551;
t425 = 0.2e1 * t449 + qJD(2);
t162 = -t328 * t477 - t329 * t478 + t332 * t410 - t406 * t334;
t49 = pkin(5) * t341 - pkin(12) * t204 + t62;
t203 = t408 * t282 + t283 * t404;
t54 = -pkin(12) * t203 + t63;
t20 = -t403 * t54 + t541 * t49;
t21 = t403 * t49 + t541 * t54;
t126 = -t403 * t203 + t541 * t204;
t281 = -qJD(3) * t341 + t410 * t459;
t191 = -qJD(4) * t282 + t281 * t409 + t405 * t460;
t280 = qJD(3) * t342 + t406 * t459;
t161 = -t328 * t478 + t329 * t477 + t406 * t332 + t410 * t334;
t149 = pkin(10) * t460 + t161;
t183 = t280 * pkin(3) - t281 * pkin(10) + t335;
t67 = -qJD(4) * t140 - t405 * t149 + t409 * t183;
t46 = t280 * pkin(4) - t191 * pkin(11) + t67;
t190 = -qJD(4) * t466 + t281 * t405 + t342 * t474 - t409 * t460;
t66 = t409 * t149 + t405 * t183 + t227 * t474 - t229 * t475;
t56 = -pkin(11) * t190 + t66;
t12 = t113 * t472 - t119 * t473 + t404 * t46 + t408 * t56;
t419 = -pkin(10) * t261 + t197 * t305;
t179 = pkin(4) * t282 + t228;
t415 = pkin(1) * (-qJD(2) * t449 + t512);
t150 = -pkin(3) * t460 - t162;
t13 = -qJD(5) * t63 - t404 * t56 + t408 * t46;
t86 = t160 * pkin(4) + t136;
t105 = pkin(4) * t190 + t150;
t340 = pkin(4) * t462 + t403 * t395;
t339 = -pkin(4) * t509 + t541 * t395;
t333 = t486 * qJD(1);
t326 = pkin(5) * t357 + t396;
t321 = -pkin(9) * t505 + t355;
t275 = pkin(5) * t336 + t366;
t270 = t541 * t357 + t358 * t403;
t245 = t541 * t336 - t337 * t403;
t130 = pkin(4) * t257 + pkin(5) * t429;
t125 = t541 * t203 + t204 * t403;
t115 = pkin(5) * t203 + t179;
t82 = qJD(5) * t204 + t408 * t190 + t404 * t191;
t81 = t404 * t190 - t408 * t191 + t282 * t472 + t283 * t473;
t50 = pkin(5) * t82 + t105;
t39 = pkin(5) * t424 + t86;
t31 = t126 * qJD(6) - t403 * t81 + t541 * t82;
t30 = t203 * t454 + t204 * t471 + t403 * t82 + t541 * t81;
t17 = t541 * t37 - t533;
t16 = -t403 * t37 - t467;
t11 = -pkin(12) * t82 + t12;
t10 = t280 * pkin(5) + t81 * pkin(12) + t13;
t4 = -t21 * qJD(6) + t541 * t10 - t403 * t11;
t3 = t20 * qJD(6) + t403 * t10 + t541 * t11;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t428, t485 * t559, t425 * t459, -0.2e1 * t428, -t425 * t460, 0, -t325 * t528 - t335 * t426 + t407 * t434, -t324 * t528 - t334 * t426 + t411 * t434 (t324 * t411 + t325 * t407 + (-t330 * t411 - t333 * t407) * qJD(2) + (t334 * t411 + t335 * t407 + (-t349 * t411 - t407 * t486) * qJD(2)) * qJD(1)) * t402, t324 * t486 - t325 * t349 - t330 * t335 + t333 * t334, -t260 * t342 + t281 * t313, t260 * t341 - t261 * t342 - t280 * t313 - t281 * t311, -t281 * t431 + (t260 * t411 + (qJD(1) * t342 + t313) * t481) * t402, t280 * t311 + t230, t280 * t431 + (t261 * t411 + (-qJD(1) * t341 - t311) * t481) * t402 (-t399 * t482 - t444) * t481, -t162 * t431 + t335 * t311 + t327 * t261 + t325 * t341 + t288 * t280 + (-t145 * t411 + (qJD(1) * t238 + t219) * t481) * t402, t161 * t431 + t335 * t313 - t327 * t260 + t325 * t342 + t288 * t281 + (-t418 * t411 + (-qJD(1) * t239 - t220) * t481) * t402, -t145 * t342 - t161 * t311 - t162 * t313 - t219 * t281 - t220 * t280 + t238 * t260 - t239 * t261 + t341 * t418, t145 * t238 + t161 * t220 + t162 * t219 - t239 * t418 + t288 * t335 + t325 * t327, -t159 * t283 + t191 * t257, t159 * t282 - t160 * t283 - t190 * t257 - t191 * t255, -t159 * t341 + t191 * t305 + t257 * t280 + t261 * t283, t160 * t282 + t190 * t255, -t160 * t341 - t190 * t305 - t255 * t280 - t261 * t282, t280 * t305 + t230, t117 * t280 + t136 * t282 + t139 * t261 + t150 * t255 + t160 * t228 + t190 * t197 + t305 * t67 + t341 * t58, -t118 * t280 + t136 * t283 - t140 * t261 + t150 * t257 - t159 * t228 + t191 * t197 - t305 * t66 - t341 * t57, -t117 * t191 - t118 * t190 + t139 * t159 - t140 * t160 - t255 * t66 - t257 * t67 - t282 * t57 - t283 * t58, t117 * t67 + t118 * t66 + t136 * t228 + t139 * t58 + t140 * t57 + t150 * t197, -t204 * t71 - t429 * t81, t81 * t175 + t71 * t203 - t204 * t424 - t429 * t82, t204 * t261 + t280 * t429 - t294 * t81 - t341 * t71, t175 * t82 + t203 * t424, -t175 * t280 - t203 * t261 - t82 * t294 - t341 * t424, t280 * t294 + t230, t105 * t175 + t13 * t294 + t143 * t82 + t179 * t424 + t86 * t203 + t62 * t261 + t47 * t280 + t9 * t341, t105 * t429 - t12 * t294 - t143 * t81 - t179 * t71 + t204 * t86 - t261 * t63 - t280 * t48 + t341 * t450, -t12 * t175 - t13 * t429 + t203 * t450 - t9 * t204 - t424 * t63 + t47 * t81 - t48 * t82 + t62 * t71, t105 * t143 + t12 * t48 + t13 * t47 + t179 * t86 - t450 * t63 + t62 * t9, -t126 * t24 - t30 * t562, t100 * t30 + t125 * t24 - t126 * t25 - t31 * t562, t126 * t261 - t24 * t341 + t280 * t562 - t290 * t30, t100 * t31 + t125 * t25, -t100 * t280 - t125 * t261 - t25 * t341 - t290 * t31, t280 * t290 + t230, t100 * t50 + t115 * t25 + t125 * t39 + t14 * t280 + t2 * t341 + t20 * t261 + t290 * t4 + t31 * t89, -t115 * t24 + t126 * t39 - t15 * t280 - t21 * t261 - t290 * t3 - t30 * t89 + t341 * t414 + t50 * t562, -t100 * t3 + t125 * t414 - t126 * t2 + t14 * t30 - t15 * t31 + t20 * t24 - t21 * t25 - t4 * t562, t115 * t39 + t14 * t4 + t15 * t3 + t2 * t20 - t21 * t414 + t50 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, t485 * t512, -t411 * t435, t442, t407 * t435, 0, -pkin(8) * t437 + t333 * t426 + t407 * t415, pkin(8) * t438 + t330 * t426 + t411 * t415, 0, 0, -t260 * t406 - t410 * t553 (-t260 + t554) * t410 + (-t261 + t553) * t406, -t410 * t443 + (t410 * t445 + (t406 * qJD(2) - t313) * t407) * t484, -t311 * t548 - t515, t406 * t443 + (-t406 * t445 + (t311 + t480) * t407) * t484, t444 * t483, -pkin(2) * t261 + t288 * t478 + t246 * t431 - t333 * t311 - t571 * t410 + (-t219 * t407 + (-pkin(9) * t481 - t288 * t411) * t406) * t484, pkin(2) * t260 + t288 * t477 - t248 * t431 - t333 * t313 + t571 * t406 + (-t288 * t501 + (-pkin(9) * t480 + t220) * t407) * t484, t246 * t313 + t248 * t311 + ((-t261 + t479) * pkin(9) + t569) * t410 + ((-t260 + t577) * pkin(9) + t545) * t406, -pkin(2) * t325 - t219 * t246 - t220 * t248 - t288 * t333 + (-t418 * t410 - t145 * t406 + (-t219 * t410 - t220 * t406) * qJD(3)) * pkin(9), -t159 * t504 + t589 * t257, t255 * t299 + t257 * t298 + (-t255 * t409 - t257 * t405) * t477 + (t522 - t521 + (-t257 * t409 + t518) * qJD(4)) * t406, t159 * t410 - t586 * t305 + (-t257 * t431 - t305 * t475 + t503) * t406, t160 * t506 + t581 * t255, t160 * t410 + t587 * t305 + (t255 * t431 - t305 * t474 - t507) * t406, -t305 * t548 - t515, -t197 * t298 - t231 * t255 + t261 * t321 - t494 * t305 + (-t58 + (pkin(9) * t255 + t197 * t405) * qJD(3)) * t410 + (pkin(9) * t160 - t117 * t431 + t197 * t474 + t524) * t406, -t197 * t299 - t231 * t257 - t261 * t322 + t493 * t305 + (t57 + (pkin(9) * t257 + t197 * t409) * qJD(3)) * t410 + (-pkin(9) * t159 + t118 * t431 - t197 * t475 + t523) * t406, t117 * t299 + t118 * t298 + t159 * t321 - t160 * t322 + t494 * t257 + t493 * t255 + t430 * t477 + (-t405 * t57 - t409 * t58 + (t117 * t405 - t118 * t409) * qJD(4)) * t406, -t197 * t231 + t321 * t58 + t322 * t57 - t493 * t118 - t494 * t117 + (t136 * t406 + t197 * t477) * pkin(9), t71 * t337 - t429 * t492, t175 * t492 + t71 * t336 + t337 * t424 - t429 * t491, -t337 * t261 - t294 * t492 + t71 * t410 - t429 * t548, t175 * t491 + t336 * t424, t175 * t548 - t336 * t261 - t294 * t491 + t410 * t424, -t294 * t548 - t515, t143 * t491 + t175 * t488 + t200 * t261 + t294 * t534 + t86 * t336 + t366 * t424 - t9 * t410 - t47 * t548, -t143 * t492 - t201 * t261 + t294 * t535 - t337 * t86 - t366 * t71 - t410 * t450 + t429 * t488 + t48 * t548, t175 * t535 + t200 * t71 - t201 * t424 + t336 * t450 + t9 * t337 - t429 * t534 + t47 * t492 - t48 * t491, t143 * t488 + t200 * t9 - t201 * t450 + t366 * t86 + t47 * t534 - t48 * t535, -t24 * t247 - t500 * t562, t100 * t500 + t24 * t245 - t247 * t25 - t499 * t562, t24 * t410 + t247 * t261 - t290 * t500 - t548 * t562, t100 * t499 + t25 * t245, t100 * t548 - t245 * t261 + t25 * t410 - t290 * t499, -t290 * t548 - t515, t100 * t497 - t14 * t548 - t2 * t410 + t245 * t39 + t25 * t275 + t261 * t96 + t290 * t536 + t499 * t89, t15 * t548 - t24 * t275 + t247 * t39 - t261 * t97 - t290 * t537 - t410 * t414 + t497 * t562 - t500 * t89, -t100 * t537 + t14 * t500 - t15 * t499 - t2 * t247 + t24 * t96 + t245 * t414 - t25 * t97 - t536 * t562, t14 * t536 + t15 * t537 + t2 * t96 + t275 * t39 - t414 * t97 + t497 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t513, -t311 ^ 2 + t313 ^ 2, -t260 - t554, -t513, -t313 * t384 - t427, t438, -t288 * t313 - t545, t288 * t311 - t569, 0, 0, t257 * t447 - t522 (-t159 - t519) * t409 + (-t160 - t516) * t405, -t257 * t313 + t305 * t447 + t507, t305 * t518 - t521, -t305 ^ 2 * t405 + t255 * t313 + t503, -t305 * t313, -pkin(3) * t160 - t117 * t313 - t523 - t220 * t255 + (-pkin(10) * t474 - t141) * t305 + t419 * t405, pkin(3) * t159 + t118 * t313 + t524 - t220 * t257 + (pkin(10) * t475 + t142) * t305 + t419 * t409, t141 * t257 + t142 * t255 + ((-t160 + t476) * pkin(10) + t556) * t409 + ((qJD(4) * t255 - t159) * pkin(10) + t555) * t405, -pkin(3) * t136 - t117 * t141 - t118 * t142 - t197 * t220 + (qJD(4) * t430 - t405 * t58 + t409 * t57) * pkin(10), -t71 * t358 - t429 * t490, t175 * t490 + t71 * t357 - t358 * t424 - t429 * t489, t358 * t261 - t294 * t490 - t313 * t429, t175 * t489 + t357 * t424, t175 * t313 - t357 * t261 - t294 * t489, -t294 * t313, t143 * t489 + t175 * t436 + t292 * t261 + t294 * t529 - t47 * t313 + t86 * t357 + t396 * t424, -t143 * t490 - t261 * t293 - t294 * t530 + t313 * t48 + t358 * t86 - t396 * t71 + t429 * t436, -t175 * t530 + t292 * t71 - t293 * t424 + t357 * t450 - t9 * t358 - t429 * t529 + t47 * t490 - t48 * t489, t143 * t436 + t292 * t9 - t293 * t450 + t396 * t86 + t47 * t529 + t48 * t530, -t24 * t271 - t495 * t562, t100 * t495 + t24 * t270 - t271 * t25 + t496 * t562, t271 * t261 - t290 * t495 - t313 * t562, -t100 * t496 + t25 * t270, t100 * t313 - t270 * t261 + t290 * t496, -t290 * t313, t498 * t100 - t14 * t313 + t170 * t261 + t25 * t326 + t270 * t39 - t290 * t539 - t496 * t89, t15 * t313 - t171 * t261 - t24 * t326 + t271 * t39 + t290 * t538 - t495 * t89 + t498 * t562, t538 * t100 + t495 * t14 + t496 * t15 + t170 * t24 - t171 * t25 - t2 * t271 + t270 * t414 + t539 * t562, -t14 * t539 - t538 * t15 + t170 * t2 - t171 * t414 + t326 * t39 + t498 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t517, -t255 ^ 2 + t257 ^ 2, -t159 + t519, -t517, t516 - t160, t261, -t197 * t257 - t555, t197 * t255 - t556, 0, 0, t520, t567, t564, -t520, t543, t261, -t294 * t51 + (-t175 * t257 + t261 * t408 - t294 * t473) * pkin(4) + t546, t294 * t52 + (-t257 * t429 - t261 * t404 - t294 * t472) * pkin(4) + t563, t48 * t429 + t52 * t175 - t47 * t175 + t51 * t429 + (-t404 * t424 + t408 * t71 + (-t175 * t408 + t404 * t429) * qJD(5)) * pkin(4), -t47 * t51 - t48 * t52 + (-t143 * t257 - t404 * t450 + t408 * t9 + (-t404 * t47 + t408 * t48) * qJD(5)) * pkin(4), t527, t568, t566, -t527, t544, t261, -t130 * t100 + t339 * t261 - t290 * t532 + t547, -t130 * t562 - t340 * t261 + t290 * t531 + t565, t100 * t531 + t24 * t339 - t25 * t340 + t532 * t562 + t580, -t89 * t130 - t14 * t532 - t15 * t531 + t2 * t339 - t340 * t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t520, t567, t564, -t520, t543, t261, t48 * t294 + t546, t294 * t47 + t563, 0, 0, t527, t568, t566, -t527, t544, t261, -t16 * t290 + (-t100 * t429 + t541 * t261 - t290 * t471) * pkin(5) + t547, t17 * t290 + (-t261 * t403 - t290 * t454 - t429 * t562) * pkin(5) + t565, t17 * t100 + t16 * t562 + (t541 * t24 - t25 * t403 + (-t541 * t100 + t403 * t562) * qJD(6)) * pkin(5) + t580, -t14 * t16 - t15 * t17 + (t541 * t2 - t414 * t403 - t429 * t89 + (-t14 * t403 + t541 * t15) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t527, t568, t566, -t527, t544, t261, t15 * t290 + t547, t14 * t290 + t565, 0, 0;];
tauc_reg  = t1;
