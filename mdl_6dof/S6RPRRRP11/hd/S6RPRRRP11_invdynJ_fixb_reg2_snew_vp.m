% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 02:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:08:59
% EndTime: 2019-05-06 02:09:38
% DurationCPUTime: 25.96s
% Computational Cost: add. (182753->782), mult. (569683->1201), div. (0->0), fcn. (489844->14), ass. (0->492)
t338 = sin(pkin(12));
t339 = sin(pkin(7));
t340 = sin(pkin(6));
t341 = cos(pkin(12));
t342 = cos(pkin(7));
t559 = sin(qJ(1));
t560 = cos(qJ(1));
t393 = g(1) * t559 - g(2) * t560;
t375 = qJDD(1) * pkin(1) + t393;
t343 = cos(pkin(6));
t450 = t343 * g(3) - qJDD(2);
t470 = qJDD(1) * t340;
t453 = t341 * t470;
t352 = -pkin(2) * t453 - t340 * t375 - t450;
t364 = t343 * t375;
t469 = qJDD(1) * t343;
t461 = pkin(2) * t469;
t494 = t340 * t341;
t353 = -g(3) * t494 + t341 * t364 + t461;
t394 = g(1) * t560 + g(2) * t559;
t366 = qJ(2) * t470 - t394;
t337 = t343 ^ 2;
t547 = pkin(9) * t339;
t387 = (t338 * pkin(1) + t337 * t547) * qJD(1);
t502 = t338 * t339;
t557 = pkin(2) * t341;
t405 = -pkin(9) * t502 - t557;
t480 = qJD(1) * t343;
t481 = qJD(1) * t340;
t546 = pkin(9) * t342;
t558 = pkin(2) * t338;
t334 = t338 ^ 2;
t336 = t341 ^ 2;
t569 = -t334 - t336;
t586 = -2 * qJD(2);
t584 = t338 * t586;
t230 = qJD(1) * (t339 * t387 + t340 * (t480 * ((qJ(2) + 0.2e1 * t546) * t339 * t341 - t342 * t558) - (t405 * t502 + t342 * (t546 * t569 - qJ(2))) * t481 + t339 * t584)) + t339 * (-t338 * t366 + t353) - t342 * t352;
t349 = cos(qJ(3));
t346 = sin(qJ(3));
t491 = t342 * t346;
t497 = t339 * t346;
t361 = t343 * t497 + (t338 * t349 + t341 * t491) * t340;
t308 = t361 * qJD(1);
t464 = t339 * t494;
t315 = -qJD(1) * t464 + t342 * t480 + qJD(3);
t345 = sin(qJ(4));
t348 = cos(qJ(4));
t298 = t308 * t348 + t315 * t345;
t344 = sin(qJ(5));
t347 = cos(qJ(5));
t500 = t338 * t346;
t493 = t341 * t342;
t462 = t340 * t493;
t496 = t339 * t349;
t578 = -t343 * t496 - t349 * t462;
t306 = qJD(1) * t578 + t481 * t500;
t382 = qJD(4) + t306;
t274 = t298 * t344 - t347 * t382;
t276 = t347 * t298 + t344 * t382;
t234 = t276 * t274;
t359 = t361 * qJDD(1);
t288 = -t306 * qJD(3) + t359;
t314 = -t339 * t453 + t342 * t469 + qJDD(3);
t448 = t345 * t288 - t348 * t314;
t251 = -t298 * qJD(4) - t448;
t250 = qJDD(5) - t251;
t572 = -t234 + t250;
t583 = pkin(5) * t572;
t296 = t308 * t345 - t348 * t315;
t409 = -t348 * t288 - t345 * t314;
t252 = -t296 * qJD(4) - t409;
t454 = t338 * t470;
t391 = qJDD(1) * t578 + t346 * t454;
t374 = -qJD(3) * t308 - t391;
t365 = qJDD(4) - t374;
t203 = -t274 * qJD(5) + t347 * t252 + t344 * t365;
t291 = qJD(5) + t296;
t247 = t291 * t274;
t183 = t247 + t203;
t582 = qJ(6) * t183;
t265 = t298 * t296;
t570 = -t265 + t365;
t581 = t345 * t570;
t580 = t348 * t570;
t515 = t572 * t344;
t514 = t572 * t347;
t272 = t276 ^ 2;
t290 = t291 ^ 2;
t231 = -t272 - t290;
t271 = t274 ^ 2;
t449 = t252 * t344 - t347 * t365;
t202 = -qJD(5) * t276 - t449;
t243 = pkin(5) * t291 - qJ(6) * t276;
t537 = t340 * g(3);
t358 = t364 - t537;
t388 = t405 * t340;
t452 = qJ(2) + t546;
t443 = t343 * t452;
t563 = 2 * qJD(2);
t498 = t339 * t343;
t398 = t462 + t498;
t575 = pkin(9) * t398;
t256 = t341 * t366 + t358 * t338 + qJDD(1) * t575 + ((-pkin(1) * t341 - pkin(2) * t337) * qJD(1) + (t341 * t563 + (t338 * t443 + t341 * t388) * qJD(1)) * t340) * qJD(1);
t492 = t342 * t343;
t495 = t340 * qJ(2);
t564 = qJD(1) ^ 2;
t568 = -t464 + t492;
t190 = t349 * t256 + (t342 * (t564 * pkin(1) + t394) + ((-t339 ^ 2 * pkin(9) - t342 * t452) * qJDD(1) + (pkin(2) * qJD(1) * t398 + t342 * t586) * qJD(1)) * t340) * t500 + (t339 * t352 + t342 * t353 + ((-t339 * t340 + t341 * t492) * t495 + t568 * t575) * t564) * t346;
t287 = pkin(3) * t306 - pkin(10) * t308;
t565 = t315 ^ 2;
t168 = -pkin(3) * t565 + t314 * pkin(10) - t306 * t287 + t190;
t473 = qJD(3) + t315;
t258 = t308 * t473 + t391;
t301 = t315 * t306;
t262 = -t301 + t288;
t350 = t258 * pkin(3) - t262 * pkin(10) - t230;
t130 = t348 * t168 + t345 * t350;
t257 = pkin(4) * t296 - pkin(11) * t298;
t377 = t382 ^ 2;
t108 = -pkin(4) * t377 + pkin(11) * t365 - t296 * t257 + t130;
t397 = -pkin(1) + t405;
t476 = t340 * t564;
t490 = t342 * t349;
t189 = t346 * t256 - (t461 + t358 * t341 + (-t452 * t470 + t394) * t338 + (t387 + (t584 + (-t338 * t388 + t341 * t443) * qJD(1)) * t340) * qJD(1)) * t490 - ((qJDD(1) * t397 - t393) * t340 + (t343 * t558 - t495 + (t340 * t342 * t569 - t341 * t498) * pkin(9)) * t476 - t450) * t496;
t167 = -t314 * pkin(3) - t565 * pkin(10) + t308 * t287 + t189;
t279 = t382 * t296;
t226 = t252 - t279;
t376 = t382 * t298;
t127 = -t226 * pkin(11) + (-t251 + t376) * pkin(4) + t167;
t84 = t347 * t108 + t344 * t127;
t396 = t202 * qJ(6) - 0.2e1 * qJD(6) * t274 - t243 * t291 + t84;
t579 = -t396 + pkin(5) * (t231 + t271);
t506 = t308 * t306;
t373 = t314 - t506;
t574 = t346 * t373;
t573 = t349 * t373;
t571 = -t247 + t203;
t180 = (qJD(5) - t291) * t276 + t449;
t294 = t296 ^ 2;
t295 = t298 ^ 2;
t304 = t306 ^ 2;
t305 = t308 ^ 2;
t477 = qJD(6) * t276;
t268 = -0.2e1 * t477;
t83 = t108 * t344 - t347 * t127;
t383 = -t582 - t83 + t583;
t67 = t268 + t383;
t562 = pkin(5) * t67;
t129 = t168 * t345 - t348 * t350;
t90 = t129 * t345 + t348 * t130;
t79 = t346 * t167 + t349 * t90;
t561 = pkin(9) * t79;
t136 = -t180 * t347 + t183 * t344;
t211 = -t271 - t272;
t104 = t136 * t345 - t211 * t348;
t556 = pkin(3) * t104;
t218 = -t290 - t271;
t155 = t218 * t347 - t515;
t179 = (qJD(5) + t291) * t276 + t449;
t114 = t155 * t345 - t179 * t348;
t555 = pkin(3) * t114;
t194 = t234 + t250;
t516 = t194 * t347;
t160 = -t231 * t344 - t516;
t119 = t160 * t345 - t348 * t571;
t554 = pkin(3) * t119;
t553 = pkin(3) * t349;
t154 = t218 * t344 + t514;
t552 = pkin(4) * t154;
t517 = t194 * t344;
t159 = t231 * t347 - t517;
t551 = pkin(4) * t159;
t550 = pkin(4) * t345;
t549 = pkin(4) * t348;
t548 = pkin(5) * t183;
t545 = pkin(10) * t104;
t544 = pkin(10) * t114;
t543 = pkin(10) * t119;
t134 = -t180 * t344 - t183 * t347;
t542 = pkin(11) * t134;
t541 = pkin(11) * t154;
t540 = pkin(11) * t159;
t529 = t344 * t67;
t68 = -pkin(5) * t271 + t396;
t37 = t347 * t68 - t529;
t107 = -t365 * pkin(4) - t377 * pkin(11) + t257 * t298 + t129;
t91 = -t202 * pkin(5) - t271 * qJ(6) + t243 * t276 + qJDD(6) + t107;
t26 = t345 * t37 - t348 * t91;
t27 = t345 * t91 + t348 * t37;
t527 = t347 * t67;
t36 = t344 * t68 + t527;
t437 = t27 * t346 - t349 * t36;
t5 = -t339 * t26 + t342 * t437;
t539 = t338 * t5;
t51 = t344 * t83 + t347 * t84;
t40 = -t107 * t348 + t345 * t51;
t41 = t107 * t345 + t348 * t51;
t50 = t344 * t84 - t347 * t83;
t434 = t346 * t41 - t349 * t50;
t8 = -t339 * t40 + t342 * t434;
t538 = t338 * t8;
t105 = t136 * t348 + t211 * t345;
t431 = t105 * t346 - t134 * t349;
t62 = t342 * t104 + t339 * t431;
t63 = -t339 * t104 + t342 * t431;
t87 = t105 * t349 + t346 * t134;
t536 = pkin(1) * (-t340 * t62 + (t338 * t87 + t341 * t63) * t343) + (-t338 * t63 + t341 * t87) * t495;
t115 = t155 * t348 + t179 * t345;
t429 = t115 * t346 - t154 * t349;
t72 = t342 * t114 + t339 * t429;
t73 = -t339 * t114 + t342 * t429;
t96 = t115 * t349 + t346 * t154;
t535 = pkin(1) * (-t340 * t72 + (t338 * t96 + t341 * t73) * t343) + (-t338 * t73 + t341 * t96) * t495;
t120 = t160 * t348 + t345 * t571;
t428 = t120 * t346 - t159 * t349;
t75 = t342 * t119 + t339 * t428;
t76 = -t339 * t119 + t342 * t428;
t99 = t120 * t349 + t346 * t159;
t534 = pkin(1) * (-t340 * t75 + (t338 * t99 + t341 * t76) * t343) + (-t338 * t76 + t341 * t99) * t495;
t533 = pkin(2) * t63 + t87 * t547;
t532 = pkin(2) * t73 + t96 * t547;
t531 = pkin(2) * t76 + t99 * t547;
t432 = -t167 * t349 + t346 * t90;
t89 = -t129 * t348 + t130 * t345;
t53 = -t339 * t89 + t342 * t432;
t530 = t338 * t53;
t528 = t346 * t89;
t526 = t349 * t89;
t525 = t107 * t344;
t524 = t107 * t347;
t223 = -t298 * t306 + t448;
t227 = t252 + t279;
t170 = -t223 * t345 - t227 * t348;
t172 = -t223 * t348 + t227 * t345;
t241 = t294 + t295;
t422 = t172 * t346 + t241 * t349;
t110 = -t339 * t170 + t342 * t422;
t523 = t110 * t338;
t253 = -t377 - t294;
t196 = t253 * t345 + t580;
t197 = t253 * t348 - t581;
t381 = 0.2e1 * qJD(4) + t306;
t224 = -t298 * t381 - t448;
t419 = t197 * t346 + t224 * t349;
t140 = -t339 * t196 + t342 * t419;
t522 = t140 * t338;
t255 = -t295 - t377;
t238 = t265 + t365;
t511 = t238 * t345;
t206 = t255 * t348 - t511;
t510 = t238 * t348;
t207 = -t255 * t345 - t510;
t228 = t296 * t381 + t409;
t418 = t207 * t346 + t228 * t349;
t142 = -t339 * t206 + t342 * t418;
t521 = t142 * t338;
t520 = t167 * t345;
t519 = t167 * t348;
t273 = -t304 - t305;
t263 = t301 + t288;
t362 = (-qJD(3) + t315) * t308 - t391;
t414 = -t263 * t349 + t346 * t362;
t192 = -t339 * t273 + t342 * t414;
t518 = t192 * t338;
t513 = t230 * t346;
t512 = t230 * t349;
t281 = -t314 - t506;
t509 = t281 * t349;
t508 = t291 * t344;
t507 = t291 * t347;
t505 = t308 * t346;
t335 = t340 ^ 2;
t504 = t334 * t335;
t503 = t335 * t336;
t501 = t338 * t340;
t489 = t346 * t281;
t488 = -pkin(3) * t134 + pkin(10) * t105;
t487 = -pkin(3) * t154 + pkin(10) * t115;
t486 = -pkin(3) * t159 + pkin(10) * t120;
t485 = -pkin(4) * t211 + pkin(11) * t136;
t484 = pkin(4) * t179 - pkin(11) * t155;
t483 = pkin(4) * t571 - pkin(11) * t160;
t475 = t341 * t564;
t474 = t343 * t564;
t471 = qJDD(1) * t335;
t467 = t345 * t234;
t466 = t348 * t234;
t465 = t349 * t265;
t460 = t346 * t265;
t459 = -pkin(2) * t62 + t87 * t546;
t458 = -pkin(2) * t72 + t96 * t546;
t457 = -pkin(2) * t75 + t99 * t546;
t455 = t340 * t474;
t446 = -pkin(4) * t107 + pkin(11) * t51;
t445 = -t484 - t524;
t444 = -t483 + t525;
t280 = -t565 - t304;
t240 = t280 * t349 - t574;
t442 = pkin(9) * t240 + t512;
t284 = -t305 - t565;
t242 = -t346 * t284 + t509;
t441 = pkin(9) * t242 - t513;
t100 = -pkin(4) * t134 + t548;
t57 = -qJ(6) * t180 + (-t211 - t271) * pkin(5) + t396;
t269 = 0.2e1 * t477;
t58 = t269 - t383 + t582;
t25 = -t344 * t57 + t347 * t58 - t542;
t14 = -t100 * t345 + t25 * t348 - t545;
t390 = t344 * t58 + t347 * t57 + t485;
t16 = -t390 - t556;
t440 = t14 * t346 + t16 * t349;
t372 = t383 + t583;
t55 = t269 - t372 - t552;
t77 = -pkin(5) * t179 + qJ(6) * t218 - t91;
t59 = -qJ(6) * t514 - t344 * t77 - t541;
t21 = -t345 * t55 + t348 * t59 - t544;
t379 = -qJ(6) * t515 + t347 * t77 - t484;
t44 = -t379 - t555;
t439 = t21 * t346 + t349 * t44;
t56 = -t551 - t579;
t144 = -pkin(5) * t571 - qJ(6) * t194;
t88 = -qJ(6) * t231 + t91;
t60 = -t144 * t344 + t347 * t88 - t540;
t24 = -t345 * t56 + t348 * t60 - t543;
t389 = t144 * t347 + t344 * t88 - t483;
t45 = -t389 - t554;
t438 = t24 * t346 + t349 * t45;
t43 = -t50 - t542;
t28 = t134 * t550 + t348 * t43 - t545;
t406 = t51 + t485;
t31 = -t406 - t556;
t436 = t28 * t346 + t31 * t349;
t69 = t83 - t552;
t92 = t525 - t541;
t38 = -t345 * t69 + t348 * t92 - t544;
t64 = -t445 - t555;
t435 = t346 * t38 + t349 * t64;
t70 = t84 - t551;
t93 = t524 - t540;
t42 = -t345 * t70 + t348 * t93 - t543;
t66 = -t444 - t554;
t433 = t346 * t42 + t349 * t66;
t135 = -t179 * t347 - t344 * t571;
t232 = t272 - t271;
t112 = t135 * t348 + t232 * t345;
t133 = -t179 * t344 + t347 * t571;
t430 = t112 * t346 - t133 * t349;
t245 = -t272 + t290;
t164 = -t245 * t344 + t514;
t123 = t164 * t348 + t183 * t345;
t162 = t245 * t347 + t515;
t427 = t123 * t346 - t162 * t349;
t244 = t271 - t290;
t165 = t244 * t347 - t517;
t124 = t165 * t348 - t180 * t345;
t163 = t244 * t344 + t516;
t426 = t124 * t346 - t163 * t349;
t174 = -t202 * t344 + t274 * t507;
t149 = t174 * t348 - t467;
t173 = -t347 * t202 - t274 * t508;
t425 = t149 * t346 + t173 * t349;
t176 = t203 * t347 - t276 * t508;
t150 = t176 * t348 + t467;
t175 = t203 * t344 + t276 * t507;
t424 = t150 * t346 - t175 * t349;
t171 = t224 * t348 - t226 * t345;
t264 = t295 - t294;
t423 = t171 * t346 - t264 * t349;
t205 = (-t274 * t347 + t276 * t344) * t291;
t186 = t205 * t348 + t250 * t345;
t204 = (-t274 * t344 - t276 * t347) * t291;
t421 = t186 * t346 - t204 * t349;
t420 = t189 * t349 - t346 * t190;
t146 = t346 * t189 + t190 * t349;
t278 = -t295 + t377;
t216 = -t278 * t345 + t580;
t417 = t216 * t346 - t227 * t349;
t277 = t294 - t377;
t217 = t277 * t348 - t511;
t416 = t217 * t346 + t223 * t349;
t415 = -t258 * t346 + t262 * t349;
t413 = t280 * t346 + t573;
t412 = t284 * t349 + t489;
t299 = t304 - t565;
t411 = t299 * t346 - t509;
t300 = -t305 + t565;
t410 = t300 * t349 + t574;
t355 = (-pkin(1) * qJD(1) + t340 * t563) * qJD(1) + t366;
t360 = qJ(2) * t476 + t375;
t356 = t343 * t360 - t537;
t292 = t338 * t355 - t341 * t356;
t293 = t338 * t356 + t341 * t355;
t408 = t338 * t292 + t341 * t293;
t407 = -t306 * t346 - t308 * t349;
t403 = (-t339 * t62 - t342 * t63) * pkin(9);
t402 = (-t339 * t72 - t342 * t73) * pkin(9);
t401 = (-t339 * t75 - t342 * t76) * pkin(9);
t368 = t348 * t279;
t220 = -t345 * t251 + t368;
t400 = t220 * t346 + t465;
t369 = t345 * t376;
t222 = t348 * t252 - t369;
t399 = t222 * t346 - t465;
t10 = t27 * t349 + t346 * t36;
t48 = -pkin(5) * t91 + qJ(6) * t68;
t11 = -pkin(11) * t36 - qJ(6) * t527 - t344 * t48;
t18 = -pkin(4) * t36 - t562;
t2 = -pkin(10) * t26 + t11 * t348 - t18 * t345;
t371 = -pkin(4) * t91 + pkin(11) * t37 - qJ(6) * t529 + t347 * t48;
t3 = -pkin(3) * t26 - t371;
t395 = pkin(9) * t10 + t2 * t346 + t3 * t349;
t13 = -pkin(3) * t40 - t446;
t15 = t346 * t50 + t349 * t41;
t9 = -pkin(10) * t40 + (-pkin(11) * t348 + t550) * t50;
t392 = pkin(9) * t15 + t13 * t349 + t346 * t9;
t101 = -pkin(3) * t196 + t129;
t143 = -pkin(10) * t196 + t520;
t156 = t197 * t349 - t346 * t224;
t386 = pkin(9) * t156 + t101 * t349 + t143 * t346;
t102 = -pkin(3) * t206 + t130;
t145 = -pkin(10) * t206 + t519;
t161 = t207 * t349 - t346 * t228;
t385 = pkin(9) * t161 + t102 * t349 + t145 * t346;
t229 = t346 * t263 + t349 * t362;
t384 = pkin(9) * t229 + t146;
t151 = t172 * t349 - t346 * t241;
t81 = -pkin(10) * t170 - t89;
t380 = pkin(9) * t151 - t170 * t553 + t346 * t81;
t370 = t345 * t279;
t367 = t348 * t376;
t363 = t349 * t374;
t236 = -t368 + t369;
t357 = t346 * t236 - t349 * t365;
t330 = t341 * t455;
t329 = t338 * t455;
t328 = t335 * t338 * t475;
t324 = (-t337 - t503) * t564;
t323 = (-t337 - t504) * t564;
t322 = t569 * t335 * t564;
t321 = -t328 + t469;
t320 = t328 + t469;
t319 = -t329 + t453;
t318 = t329 + t453;
t317 = -t330 + t454;
t316 = t330 + t454;
t309 = t340 * t360 + t450;
t289 = t305 - t304;
t261 = -t306 * t473 + t359;
t235 = -t370 - t367;
t221 = t345 * t252 + t367;
t219 = t348 * t251 + t370;
t215 = t277 * t345 + t510;
t214 = t278 * t348 + t581;
t213 = -t339 * t261 + t342 * t412;
t212 = t342 * t261 + t339 * t412;
t209 = -t339 * t258 + t342 * t413;
t208 = t342 * t258 + t339 * t413;
t191 = t342 * t273 + t339 * t414;
t185 = t205 * t345 - t250 * t348;
t169 = t224 * t345 + t226 * t348;
t148 = t176 * t345 - t466;
t147 = t174 * t345 + t466;
t141 = t342 * t206 + t339 * t418;
t139 = t342 * t196 + t339 * t419;
t138 = t339 * t230 - t342 * t420;
t137 = -t342 * t230 - t339 * t420;
t122 = t165 * t345 + t180 * t348;
t121 = t164 * t345 - t183 * t348;
t117 = pkin(3) * t228 + pkin(10) * t207 + t520;
t116 = pkin(3) * t224 + pkin(10) * t197 - t519;
t111 = t135 * t345 - t232 * t348;
t109 = t342 * t170 + t339 * t422;
t80 = -pkin(3) * t167 + pkin(10) * t90;
t78 = pkin(3) * t241 + pkin(10) * t172 + t90;
t65 = t343 * (t342 * t185 + t339 * t421) + (t338 * (t186 * t349 + t346 * t204) + t341 * (-t339 * t185 + t342 * t421)) * t340;
t52 = t339 * t432 + t342 * t89;
t47 = t343 * (t342 * t148 + t339 * t424) + (t338 * (t150 * t349 + t346 * t175) + t341 * (-t339 * t148 + t342 * t424)) * t340;
t46 = t343 * (t342 * t147 + t339 * t425) + (t338 * (t149 * t349 - t346 * t173) + t341 * (-t339 * t147 + t342 * t425)) * t340;
t35 = t345 * t93 + t348 * t70 + t486;
t34 = t345 * t92 + t348 * t69 + t487;
t33 = t343 * (t342 * t122 + t339 * t426) + (t338 * (t124 * t349 + t346 * t163) + t341 * (-t339 * t122 + t342 * t426)) * t340;
t32 = t343 * (t342 * t121 + t339 * t427) + (t338 * (t123 * t349 + t346 * t162) + t341 * (-t339 * t121 + t342 * t427)) * t340;
t23 = -t134 * t549 + t345 * t43 + t488;
t22 = t343 * (t342 * t111 + t339 * t430) + (t338 * (t112 * t349 + t346 * t133) + t341 * (-t339 * t111 + t342 * t430)) * t340;
t20 = t345 * t60 + t348 * t56 + t486;
t17 = t345 * t59 + t348 * t55 + t487;
t12 = t100 * t348 + t25 * t345 + t488;
t7 = t339 * t434 + t342 * t40;
t6 = pkin(10) * t41 + (-pkin(11) * t345 - pkin(3) - t549) * t50;
t4 = t342 * t26 + t339 * t437;
t1 = -pkin(3) * t36 + pkin(10) * t27 + t11 * t345 + t18 * t348;
t19 = [0, 0, 0, 0, 0, qJDD(1), t393, t394, 0, 0, t334 * t471, (t316 * t341 + t319 * t338) * t340 + (t334 - t336) * t335 * t474, t320 * t501 + t343 * t317 + t340 * (t337 - t504) * t475, t336 * t471, t321 * t494 + t343 * t318 + t338 * (-t337 + t503) * t476, t337 * qJDD(1), (-t292 + pkin(1) * (t320 * t341 + t324 * t338)) * t343 + (t341 * t309 + pkin(1) * t319 + qJ(2) * (-t320 * t338 + t324 * t341)) * t340, (-t293 + pkin(1) * (-t321 * t338 + t323 * t341)) * t343 + (-t338 * t309 - pkin(1) * t316 + qJ(2) * (-t321 * t341 - t323 * t338)) * t340, pkin(1) * (-t317 * t341 + t318 * t338) * t343 + (-pkin(1) * t322 + qJ(2) * (t317 * t338 + t318 * t341) + t408) * t340, pkin(1) * (t309 * t340 + (-t292 * t341 + t293 * t338) * t343) + t408 * t495, t343 * (t288 * t497 + (t306 * t342 + t315 * t496) * t308) + (t338 * (t288 * t349 - t315 * t505) + t341 * (t288 * t491 + (-t306 * t339 + t315 * t490) * t308)) * t340, t343 * (t342 * t289 + t339 * t415) + (t338 * (-t258 * t349 - t346 * t262) + t341 * (-t339 * t289 + t342 * t415)) * t340, t343 * (t342 * t263 + t339 * t410) + (t338 * (-t346 * t300 + t573) + t341 * (-t339 * t263 + t342 * t410)) * t340, (t301 * t349 - t346 * t374) * t501 + (t342 * t363 + (t308 * t339 + t315 * t491) * t306) * t494 + t343 * (t339 * t363 + (-t308 * t342 + t315 * t497) * t306), t343 * (t339 * t411 + t342 * t362) + (t338 * (t299 * t349 + t489) + t341 * (-t339 * t362 + t342 * t411)) * t340, t568 * t314 + ((t338 * (-t306 * t349 + t505) + t407 * t493) * t340 + t407 * t498) * t315, (pkin(2) * t209 - t342 * t189 + pkin(1) * (t209 * t341 + t240 * t338) + t442 * t339) * t343 + (t338 * (-t513 + (-t208 * t339 - t209 * t342) * pkin(9)) + t341 * (-pkin(2) * t208 + t339 * t189 + t342 * t442) - pkin(1) * t208 + qJ(2) * (-t209 * t338 + t240 * t341)) * t340, (pkin(2) * t213 - t190 * t342 + pkin(1) * (t213 * t341 + t242 * t338) + t441 * t339) * t343 + (t338 * (-t512 + (-t212 * t339 - t213 * t342) * pkin(9)) + t341 * (-pkin(2) * t212 + t190 * t339 + t342 * t441) - pkin(1) * t212 + qJ(2) * (-t213 * t338 + t242 * t341)) * t340, (pkin(2) * t192 + pkin(1) * (t192 * t341 + t229 * t338) + t384 * t339) * t343 + (t338 * t420 + qJ(2) * (t229 * t341 - t518) + t397 * t191 + (-pkin(9) * t518 + t341 * t384) * t342) * t340, (t146 * t547 + pkin(2) * t138 + pkin(1) * (t138 * t341 + t146 * t338)) * t343 + (qJ(2) * (-t138 * t338 + t146 * t341) + (-pkin(1) - t557) * t137 + (t338 * (-t137 * t339 - t138 * t342) + t146 * t493) * pkin(9)) * t340, t343 * (t342 * t221 + t339 * t399) + (t338 * (t222 * t349 + t460) + t341 * (-t339 * t221 + t342 * t399)) * t340, t343 * (t342 * t169 + t339 * t423) + (t338 * (t171 * t349 + t346 * t264) + t341 * (-t339 * t169 + t342 * t423)) * t340, t343 * (t342 * t214 + t339 * t417) + (t338 * (t216 * t349 + t346 * t227) + t341 * (-t339 * t214 + t342 * t417)) * t340, t343 * (t342 * t219 + t339 * t400) + (t338 * (t220 * t349 - t460) + t341 * (-t339 * t219 + t342 * t400)) * t340, t343 * (t342 * t215 + t339 * t416) + (t338 * (t217 * t349 - t346 * t223) + t341 * (-t339 * t215 + t342 * t416)) * t340, (t349 * t236 + t346 * t365) * t501 + (-t339 * t235 + t342 * t357) * t494 + t343 * (t342 * t235 + t339 * t357), (pkin(2) * t140 + t342 * t116 + pkin(1) * (t140 * t341 + t156 * t338) + t386 * t339) * t343 + (t338 * (-t346 * t101 - t139 * t547 + t143 * t349) + t341 * (-pkin(2) * t139 - t339 * t116) - pkin(1) * t139 + qJ(2) * (t156 * t341 - t522) + (-pkin(9) * t522 + t341 * t386) * t342) * t340, (pkin(2) * t142 + t342 * t117 + pkin(1) * (t142 * t341 + t161 * t338) + t385 * t339) * t343 + (t338 * (-t346 * t102 - t141 * t547 + t145 * t349) + t341 * (-pkin(2) * t141 - t339 * t117) - pkin(1) * t141 + qJ(2) * (t161 * t341 - t521) + (-pkin(9) * t521 + t341 * t385) * t342) * t340, (pkin(2) * t110 + t342 * t78 + pkin(1) * (t110 * t341 + t151 * t338) + t380 * t339) * t343 + (t338 * (t346 * pkin(3) * t170 - t109 * t547 + t349 * t81) + t341 * (-pkin(2) * t109 - t339 * t78) - pkin(1) * t109 + qJ(2) * (t151 * t341 - t523) + (-pkin(9) * t523 + t341 * t380) * t342) * t340, (pkin(2) * t53 + t342 * t80 + pkin(1) * (t338 * t79 + t341 * t53) + (t561 + (-pkin(10) * t346 - t553) * t89) * t339) * t343 + (t338 * (pkin(3) * t528 - pkin(10) * t526 - t52 * t547) + t341 * (-pkin(2) * t52 - t339 * t80) - pkin(1) * t52 + qJ(2) * (t341 * t79 - t530) + (-pkin(9) * t530 + t341 * (-pkin(3) * t526 - pkin(10) * t528 + t561)) * t342) * t340, t47, t22, t32, t46, t33, t65, t343 * (t339 * t435 + t342 * t34 + t532) + (t338 * (-t346 * t64 + t349 * t38 + t402) + t341 * (-t339 * t34 + t342 * t435 + t458)) * t340 + t535, t343 * (t339 * t433 + t342 * t35 + t531) + (t338 * (-t346 * t66 + t349 * t42 + t401) + t341 * (-t339 * t35 + t342 * t433 + t457)) * t340 + t534, t343 * (t342 * t23 + t339 * t436 + t533) + (t338 * (t28 * t349 - t346 * t31 + t403) + t341 * (-t339 * t23 + t342 * t436 + t459)) * t340 + t536, (pkin(2) * t8 + t342 * t6 + pkin(1) * (t15 * t338 + t341 * t8) + t392 * t339) * t343 + (t338 * (-t346 * t13 + t349 * t9 - t547 * t7) + t341 * (-pkin(2) * t7 - t339 * t6) - pkin(1) * t7 + qJ(2) * (t15 * t341 - t538) + (-pkin(9) * t538 + t341 * t392) * t342) * t340, t47, t22, t32, t46, t33, t65, t343 * (t342 * t17 + t339 * t439 + t532) + (t338 * (t21 * t349 - t346 * t44 + t402) + t341 * (-t339 * t17 + t342 * t439 + t458)) * t340 + t535, t343 * (t342 * t20 + t339 * t438 + t531) + (t338 * (t24 * t349 - t346 * t45 + t401) + t341 * (-t339 * t20 + t342 * t438 + t457)) * t340 + t534, t343 * (t342 * t12 + t339 * t440 + t533) + (t338 * (t14 * t349 - t346 * t16 + t403) + t341 * (-t339 * t12 + t342 * t440 + t459)) * t340 + t536, (pkin(2) * t5 + t342 * t1 + pkin(1) * (t10 * t338 + t341 * t5) + t395 * t339) * t343 + (t338 * (t2 * t349 - t346 * t3 - t4 * t547) + t341 * (-pkin(2) * t4 - t339 * t1) - pkin(1) * t4 + qJ(2) * (t10 * t341 - t539) + (-pkin(9) * t539 + t341 * t395) * t342) * t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319, t316, t322, -t309, 0, 0, 0, 0, 0, 0, t208, t212, t191, t137, 0, 0, 0, 0, 0, 0, t139, t141, t109, t52, 0, 0, 0, 0, 0, 0, t72, t75, t62, t7, 0, 0, 0, 0, 0, 0, t72, t75, t62, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t506, t289, t263, -t506, t362, t314, -t189, -t190, 0, 0, t221, t169, t214, t219, t215, t235, t116, t117, t78, t80, t148, t111, t121, t147, t122, t185, t34, t35, t23, t6, t148, t111, t121, t147, t122, t185, t17, t20, t12, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, t264, t227, -t265, -t223, t365, -t129, -t130, 0, 0, t175, t133, t162, -t173, t163, t204, t445, t444, t406, t446, t175, t133, t162, -t173, t163, t204, t379, t389, t390, t371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t232, t183, -t234, -t180, t250, -t83, -t84, 0, 0, t234, t232, t183, -t234, -t180, t250, t268 + t372, t579, -t548, t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t571, t211, t91;];
tauJ_reg  = t19;
