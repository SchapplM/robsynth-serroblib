% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 23:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRRPR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:59:14
% EndTime: 2019-05-07 23:00:11
% DurationCPUTime: 21.38s
% Computational Cost: add. (101884->664), mult. (219149->927), div. (0->0), fcn. (178334->12), ass. (0->419)
t384 = cos(pkin(6));
t378 = qJD(1) * t384 + qJD(2);
t387 = sin(qJ(3));
t391 = cos(qJ(3));
t383 = sin(pkin(6));
t388 = sin(qJ(2));
t459 = qJD(1) * t388;
t446 = t383 * t459;
t349 = t378 * t391 - t387 * t446;
t350 = t378 * t387 + t391 * t446;
t386 = sin(qJ(4));
t390 = cos(qJ(4));
t327 = -t390 * t349 + t350 * t386;
t324 = t327 ^ 2;
t392 = cos(qJ(2));
t460 = qJD(1) * t383;
t445 = t392 * t460;
t371 = -qJD(3) + t445;
t366 = -qJD(4) + t371;
t365 = t366 ^ 2;
t266 = -t365 - t324;
t329 = t349 * t386 + t350 * t390;
t274 = t329 * t327;
t444 = qJD(2) * t460;
t455 = qJDD(1) * t383;
t357 = -t388 * t444 + t392 * t455;
t352 = -qJDD(3) + t357;
t351 = -qJDD(4) + t352;
t409 = t351 + t274;
t533 = t390 * t409;
t187 = t266 * t386 - t533;
t535 = t386 * t409;
t189 = t266 * t390 + t535;
t127 = t187 * t387 - t189 * t391;
t437 = qJDD(1) * t384 + qJDD(2);
t462 = t388 * t455 + t392 * t444;
t407 = -t387 * t437 - t391 * t462;
t315 = t349 * qJD(3) - t407;
t408 = -t387 * t462 + t391 * t437;
t403 = -t350 * qJD(3) + t408;
t243 = qJD(4) * t329 + t315 * t386 - t390 * t403;
t477 = t366 * t329;
t531 = t243 - t477;
t586 = pkin(1) * (t127 * t388 + t392 * t531);
t514 = t329 ^ 2;
t304 = -t514 + t365;
t221 = -t304 * t390 + t535;
t225 = t304 * t386 + t533;
t163 = t221 * t391 + t225 * t387;
t585 = t384 * t163;
t124 = t187 * t391 + t189 * t387;
t584 = -pkin(8) * (t127 * t392 - t388 * t531) - pkin(1) * t124;
t524 = -t514 - t365;
t528 = -t274 + t351;
t541 = t528 * t386;
t215 = t390 * t524 + t541;
t540 = t528 * t390;
t217 = -t386 * t524 + t540;
t150 = t215 * t391 + t217 * t387;
t582 = pkin(1) * t150;
t581 = pkin(2) * t124;
t580 = pkin(2) * t150;
t579 = pkin(9) * t124;
t578 = pkin(9) * t150;
t153 = t215 * t387 - t217 * t391;
t577 = pkin(9) * t153;
t576 = t153 * t388;
t575 = t153 * t392;
t574 = t388 * (t221 * t387 - t225 * t391);
t573 = -pkin(2) * t531 - pkin(9) * t127;
t301 = t365 - t324;
t223 = t301 * t386 + t540;
t227 = t301 * t390 - t541;
t165 = t223 * t391 + t227 * t387;
t530 = t243 + t477;
t570 = t384 * t165 - (t388 * (t223 * t387 - t227 * t391) + t392 * t530) * t383;
t567 = pkin(3) * t187;
t566 = pkin(3) * t215;
t565 = pkin(10) * t187;
t564 = pkin(10) * t189;
t563 = pkin(10) * t215;
t562 = pkin(10) * t217;
t511 = -2 * qJD(5);
t398 = t390 * t315 + t386 * t403;
t244 = -t327 * qJD(4) + t398;
t309 = t327 * t366;
t527 = t309 + t244;
t435 = -pkin(2) * t392 - pkin(9) * t388;
t356 = t435 * t460;
t393 = qJD(1) ^ 2;
t507 = sin(qJ(1));
t508 = cos(qJ(1));
t418 = g(1) * t508 + g(2) * t507;
t353 = -t393 * pkin(1) + pkin(8) * t455 - t418;
t417 = g(1) * t507 - g(2) * t508;
t469 = t383 * t393;
t404 = qJDD(1) * pkin(1) + pkin(8) * t469 + t417;
t401 = t384 * t404;
t440 = t388 * t353 - t392 * t401;
t512 = t378 ^ 2;
t278 = t383 * (g(3) * t392 + t356 * t459) - t437 * pkin(2) - t512 * pkin(9) + t440;
t333 = -pkin(3) * t371 - pkin(10) * t350;
t513 = t349 ^ 2;
t201 = -t403 * pkin(3) - t513 * pkin(10) + t350 * t333 + t278;
t396 = t243 * pkin(4) - qJ(5) * t527 + t201;
t558 = t329 * t511 + t396;
t510 = -pkin(4) - pkin(11);
t526 = -t324 - t514;
t557 = pkin(2) * t526;
t555 = pkin(3) * t526;
t385 = sin(qJ(6));
t240 = qJDD(6) + t244;
t389 = cos(qJ(6));
t296 = -t389 * t327 - t366 * t385;
t298 = t327 * t385 - t366 * t389;
t248 = t298 * t296;
t529 = -t248 + t240;
t554 = t385 * t529;
t552 = t386 * t530;
t551 = t386 * t531;
t550 = t388 * t526;
t548 = t389 * t529;
t546 = t390 * t530;
t523 = t514 - t324;
t545 = t392 * t523;
t544 = t392 * t526;
t493 = t531 * t390;
t476 = t366 * t386;
t432 = t390 * t244 + t329 * t476;
t475 = t366 * t390;
t434 = t386 * t244 - t329 * t475;
t449 = t392 * t274;
t518 = t387 * t432 + t391 * t434;
t539 = t384 * t518 + (t388 * (-t387 * t434 + t391 * t432) - t449) * t383;
t422 = t243 * t386 - t327 * t475;
t433 = -t390 * t243 - t327 * t476;
t519 = t387 * t422 + t391 * t433;
t538 = t384 * t519 + (t388 * (-t387 * t433 + t391 * t422) + t449) * t383;
t420 = (t327 * t386 + t329 * t390) * t366;
t421 = (t327 * t390 - t329 * t386) * t366;
t470 = t383 * t392;
t471 = t383 * t388;
t520 = t387 * t421 + t391 * t420;
t537 = t384 * t520 + t351 * t470 + (-t387 * t420 + t391 * t421) * t471;
t207 = -t309 + t244;
t399 = -g(3) * t471 + t388 * t401;
t279 = t437 * pkin(9) - t512 * pkin(2) + (t356 * t460 + t353) * t392 + t399;
t461 = qJD(1) * t378;
t502 = t384 * g(3);
t280 = -t502 - t462 * pkin(9) - t357 * pkin(2) + ((pkin(2) * t388 - pkin(9) * t392) * t461 - t404) * t383;
t219 = t387 * t279 - t391 * t280;
t339 = t349 * t371;
t287 = t339 + t315;
t479 = t349 * t350;
t412 = -t352 + t479;
t174 = t412 * pkin(3) - pkin(10) * t287 - t219;
t220 = t391 * t279 + t387 * t280;
t175 = -pkin(3) * t513 + pkin(10) * t403 + t371 * t333 + t220;
t111 = -t390 * t174 + t386 * t175;
t269 = pkin(4) * t327 - qJ(5) * t329;
t97 = t351 * pkin(4) - t365 * qJ(5) + t329 * t269 + qJDD(5) + t111;
t69 = pkin(5) * t207 + t409 * pkin(11) + t97;
t300 = pkin(5) * t329 + pkin(11) * t366;
t442 = -pkin(4) * t366 + t511;
t73 = t243 * pkin(11) + (-t300 + t442) * t329 + t396 - t324 * pkin(5);
t43 = t385 * t73 - t389 * t69;
t44 = t385 * t69 + t389 * t73;
t25 = t385 * t44 - t389 * t43;
t380 = t383 ^ 2;
t536 = t380 * (-t384 * t393 + t461);
t534 = t387 * t412;
t532 = t391 * t412;
t317 = g(3) * t470 + t440;
t318 = t392 * t353 + t399;
t525 = t388 * t317 + t392 * t318;
t283 = (qJD(3) + t371) * t350 - t408;
t294 = t296 ^ 2;
t295 = t298 ^ 2;
t323 = qJD(6) + t329;
t321 = t323 ^ 2;
t348 = t350 ^ 2;
t368 = t371 ^ 2;
t112 = t386 * t174 + t390 * t175;
t66 = -t111 * t390 + t112 * t386;
t509 = pkin(3) * t66;
t138 = -t207 * t390 - t552;
t506 = pkin(3) * t138;
t505 = pkin(4) * t386;
t504 = pkin(4) * t390;
t503 = pkin(8) * t383;
t406 = -t365 * pkin(4) - t351 * qJ(5) - t327 * t269 + t112;
t93 = t366 * t511 + t406;
t501 = -pkin(4) * t97 + qJ(5) * t93;
t72 = -t243 * pkin(5) - t324 * pkin(11) + (t511 - t300) * t366 + t406;
t499 = t385 * t72;
t498 = t387 * t66;
t70 = t389 * t72;
t497 = t391 * t66;
t179 = t248 + t240;
t496 = t179 * t389;
t495 = t201 * t386;
t494 = t201 * t390;
t486 = t278 * t387;
t485 = t278 * t391;
t484 = t296 * t323;
t306 = t352 + t479;
t483 = t306 * t387;
t482 = t306 * t391;
t481 = t323 * t385;
t480 = t323 * t389;
t474 = t371 * t387;
t473 = t371 * t391;
t472 = t380 * t393;
t468 = t385 * t179;
t370 = t392 * t388 * t472;
t355 = t370 + t437;
t466 = t388 * t355;
t354 = -t370 + t437;
t464 = t392 * t354;
t209 = (-qJD(4) - t366) * t327 + t398;
t463 = -pkin(4) * t209 - qJ(5) * t530;
t457 = -qJD(3) + t371;
t454 = -t295 - t321;
t453 = t386 * t248;
t452 = t390 * t248;
t381 = t388 ^ 2;
t451 = t381 * t472;
t382 = t392 ^ 2;
t450 = t382 * t472;
t448 = t392 * t479;
t182 = -t296 * qJD(6) + t385 * t243 - t389 * t351;
t361 = t378 * t445;
t447 = t361 + t462;
t443 = qJ(5) * t386 + pkin(3);
t67 = t111 * t386 + t390 * t112;
t155 = t219 * t387 + t391 * t220;
t441 = t389 * t243 + t385 * t351;
t55 = t386 * t93 - t390 * t97;
t439 = pkin(3) * t55 + t501;
t438 = qJ(5) * t72 + t25 * t510;
t135 = -t209 * t390 - t552;
t436 = pkin(3) * t135 + t463;
t26 = t385 * t43 + t389 * t44;
t431 = t219 * t391 - t220 * t387;
t128 = t389 * t454 - t468;
t427 = t182 - t484;
t430 = qJ(5) * t427 + t128 * t510 + t70;
t428 = -pkin(1) + t435;
t425 = -t111 + t567;
t21 = -t25 * t390 + t386 * t72;
t423 = pkin(3) * t21 + t438;
t229 = -t321 - t294;
t121 = t385 * t229 + t548;
t167 = (qJD(6) + t323) * t298 - t441;
t419 = qJ(5) * t167 + t121 * t510 + t499;
t94 = -t128 * t390 + t386 * t427;
t416 = pkin(3) * t94 + t430;
t415 = -t112 + t566;
t171 = t182 + t484;
t413 = (-qJD(6) + t323) * t298 + t441;
t105 = -t171 * t389 + t385 * t413;
t213 = -t294 - t295;
t414 = qJ(5) * t213 + t105 * t510 - t25;
t90 = -t121 * t390 + t167 * t386;
t411 = pkin(3) * t90 + t419;
t78 = -t105 * t390 + t213 * t386;
t410 = pkin(3) * t78 + t414;
t405 = pkin(4) * t409 - qJ(5) * t266 + t97;
t402 = t405 - t567;
t400 = -pkin(4) * t524 - qJ(5) * t528 + t93;
t397 = t400 - t566;
t395 = (-qJD(4) + t366) * t327 + t398;
t360 = t378 * t446;
t359 = (t381 - t382) * t472;
t358 = -t450 - t512;
t345 = -t451 - t512;
t340 = t383 * t404 + t502;
t338 = t357 - t360;
t337 = t357 + t360;
t336 = -t361 + t462;
t335 = -t348 + t368;
t334 = -t368 + t513;
t332 = -t348 - t368;
t331 = t348 - t513;
t316 = -t368 - t513;
t299 = t348 + t513;
t293 = (-t349 * t387 + t350 * t391) * t371;
t288 = t349 * t457 + t407;
t286 = -t339 + t315;
t284 = t350 * t457 + t408;
t282 = t315 * t387 - t350 * t473;
t281 = t349 * t474 + t391 * t403;
t271 = t334 * t387 - t482;
t270 = t335 * t391 + t534;
t268 = -t332 * t387 + t482;
t267 = t332 * t391 + t483;
t258 = -t295 + t321;
t257 = t294 - t321;
t254 = t316 * t391 - t534;
t253 = t316 * t387 + t532;
t247 = t295 - t294;
t232 = -t283 * t391 + t287 * t387;
t230 = t284 * t387 + t286 * t391;
t186 = (-t296 * t389 + t298 * t385) * t323;
t185 = (t296 * t385 + t298 * t389) * t323;
t181 = -qJD(6) * t298 + t441;
t180 = pkin(2) * t288 + pkin(9) * t268 + t486;
t176 = pkin(2) * t284 + pkin(9) * t254 - t485;
t159 = t182 * t389 - t298 * t481;
t158 = -t182 * t385 - t298 * t480;
t157 = t181 * t385 - t296 * t480;
t156 = -t181 * t389 - t296 * t481;
t149 = -t185 * t386 + t240 * t390;
t148 = t185 * t390 + t240 * t386;
t147 = t257 * t389 - t468;
t146 = -t258 * t385 + t548;
t145 = -t257 * t385 - t496;
t144 = -t258 * t389 - t554;
t143 = t494 - t563;
t142 = t207 * t386 - t546;
t141 = -t386 * t395 - t493;
t140 = -t386 * t527 - t493;
t139 = t209 * t386 - t546;
t137 = t390 * t395 - t551;
t136 = t390 * t527 - t551;
t134 = t495 - t565;
t129 = -t385 * t454 - t496;
t122 = t389 * t229 - t554;
t120 = -pkin(2) * t278 + pkin(9) * t155;
t118 = -t158 * t386 + t452;
t117 = -t156 * t386 - t452;
t116 = t158 * t390 + t453;
t115 = t156 * t390 - t453;
t114 = pkin(2) * t299 + pkin(9) * t232 + t155;
t113 = -pkin(3) * t527 + t495 + t562;
t109 = -pkin(3) * t531 - t494 + t564;
t108 = t329 * t442 + t396;
t107 = t385 * t171 + t389 * t413;
t106 = -t167 * t389 - t385 * t427;
t104 = t167 * t385 - t389 * t427;
t102 = -t144 * t386 + t171 * t390;
t101 = -t145 * t386 + t390 * t413;
t100 = t144 * t390 + t171 * t386;
t99 = t145 * t390 + t386 * t413;
t98 = t148 * t391 + t149 * t387;
t95 = t128 * t386 + t390 * t427;
t91 = t121 * t386 + t167 * t390;
t89 = -t104 * t386 + t247 * t390;
t88 = t104 * t390 + t247 * t386;
t87 = -t138 * t387 + t142 * t391;
t86 = -t135 * t387 + t139 * t391;
t85 = t138 * t391 + t142 * t387;
t84 = t137 * t391 + t141 * t387;
t83 = t136 * t391 + t140 * t387;
t82 = t135 * t391 + t139 * t387;
t81 = (t531 - t477) * pkin(4) + t558;
t80 = pkin(4) * t477 + qJ(5) * t395 - t558;
t79 = t105 * t386 + t213 * t390;
t77 = -qJ(5) * t526 + t97;
t76 = -pkin(4) * t526 + t93;
t75 = t116 * t391 + t118 * t387;
t74 = t115 * t391 + t117 * t387;
t65 = t390 * t80 - t395 * t505 + t563;
t64 = qJ(5) * t493 - t386 * t81 + t565;
t63 = pkin(5) * t105 - qJ(5) * t107;
t62 = -pkin(3) * t201 + pkin(10) * t67;
t61 = -t562 + t386 * t80 + (pkin(3) + t504) * t395;
t60 = t100 * t391 + t102 * t387;
t59 = t101 * t387 + t391 * t99;
t58 = t390 * t81 + t443 * t531 - t564;
t57 = -pkin(2) * t527 + t113 * t391 + t143 * t387 - t577;
t56 = t386 * t97 + t390 * t93;
t54 = -t387 * t94 + t391 * t95;
t53 = t387 * t95 + t391 * t94;
t52 = -t387 * t90 + t391 * t91;
t51 = t387 * t91 + t391 * t90;
t50 = -pkin(10) * t138 - t66;
t49 = t387 * t89 + t391 * t88;
t48 = t109 * t391 + t134 * t387 + t573;
t47 = -t387 * t78 + t391 * t79;
t46 = t387 * t79 + t391 * t78;
t45 = pkin(10) * t142 - t555 + t67;
t40 = -pkin(10) * t135 - t386 * t76 + t390 * t77;
t39 = pkin(5) * t427 + t129 * t510 - t499;
t38 = pkin(5) * t167 + t122 * t510 + t70;
t37 = pkin(10) * t139 + t386 * t77 + t390 * t76 - t555;
t36 = t391 * t67 - t498;
t35 = t387 * t67 + t497;
t34 = -pkin(10) * t55 + (-qJ(5) * t390 + t505) * t108;
t33 = pkin(5) * t128 - qJ(5) * t129 - t44;
t32 = pkin(5) * t121 - qJ(5) * t122 - t43;
t31 = -t387 * t55 + t391 * t56;
t30 = t387 * t56 + t391 * t55;
t29 = pkin(2) * t395 + t387 * t65 + t391 * t61 + t577;
t28 = t387 * t64 + t391 * t58 - t573;
t27 = pkin(10) * t56 + (-t443 - t504) * t108;
t23 = pkin(9) * t87 + t387 * t50 + t391 * t45 - t557;
t22 = t25 * t386 + t390 * t72;
t20 = -pkin(2) * t201 + pkin(9) * t36 - pkin(10) * t498 + t391 * t62;
t19 = pkin(9) * t86 + t37 * t391 + t387 * t40 - t557;
t18 = pkin(5) * t213 + t107 * t510 - t26;
t17 = -pkin(10) * t94 + t33 * t390 - t386 * t39;
t16 = -pkin(10) * t90 + t32 * t390 - t38 * t386;
t15 = -pkin(3) * t129 + pkin(10) * t95 + t33 * t386 + t39 * t390;
t14 = -pkin(3) * t122 + pkin(10) * t91 + t32 * t386 + t38 * t390;
t13 = -pkin(10) * t78 - t18 * t386 + t390 * t63;
t12 = pkin(5) * t25 - qJ(5) * t26;
t11 = -pkin(3) * t107 + pkin(10) * t79 + t18 * t390 + t386 * t63;
t10 = pkin(5) * t72 + t26 * t510;
t9 = -pkin(2) * t108 + pkin(9) * t31 + t27 * t391 + t34 * t387;
t8 = -t21 * t387 + t22 * t391;
t7 = t21 * t391 + t22 * t387;
t6 = -pkin(2) * t129 + pkin(9) * t54 + t15 * t391 + t17 * t387;
t5 = -pkin(2) * t122 + pkin(9) * t52 + t14 * t391 + t16 * t387;
t4 = -pkin(2) * t107 + pkin(9) * t47 + t11 * t391 + t13 * t387;
t3 = -pkin(10) * t21 - t10 * t386 + t12 * t390;
t2 = -pkin(3) * t26 + pkin(10) * t22 + t10 * t390 + t12 * t386;
t1 = -pkin(2) * t26 + pkin(9) * t8 + t2 * t391 + t3 * t387;
t24 = [0, 0, 0, 0, 0, qJDD(1), t417, t418, 0, 0, (t383 * t462 + t392 * t536) * t388, t384 * t359 + (t388 * t338 + t392 * t447) * t383, t384 * t336 + (t466 + t392 * (-t451 + t512)) * t383, (t357 * t383 - t388 * t536) * t392, t384 * t337 + (t388 * (t450 - t512) + t464) * t383, t384 * t437, (-t317 + pkin(1) * (t355 * t392 + t358 * t388)) * t384 + (t392 * t340 + pkin(1) * t338 + pkin(8) * (t358 * t392 - t466)) * t383, -t340 * t471 - t384 * t318 + pkin(1) * (-t383 * t447 + (t345 * t392 - t354 * t388) * t384) + (-t345 * t388 - t464) * t503, pkin(1) * ((-t336 * t392 + t337 * t388) * t384 - (-t381 - t382) * t380 * t469) + (t336 * t388 + t337 * t392) * t503 + t525 * t383, pkin(1) * (t383 * t340 + (-t317 * t392 + t318 * t388) * t384) + t525 * t503, t384 * t282 + (t388 * (t315 * t391 + t350 * t474) + t448) * t383, t384 * t230 + (t388 * (t284 * t391 - t286 * t387) - t392 * t331) * t383, t384 * t270 + (t388 * (-t335 * t387 + t532) - t392 * t287) * t383, t384 * t281 + (t388 * (t349 * t473 - t387 * t403) - t448) * t383, t384 * t271 + (t388 * (t334 * t391 + t483) + t392 * t283) * t383, t352 * t470 + t384 * t293 + (-t349 * t391 - t350 * t387) * t371 * t471, (t176 + pkin(1) * (t254 * t388 + t284 * t392)) * t384 + (t388 * (-pkin(9) * t253 + t486) + t392 * (-pkin(2) * t253 + t219) - pkin(1) * t253 + pkin(8) * (t254 * t392 - t284 * t388)) * t383, (t180 + pkin(1) * (t268 * t388 + t288 * t392)) * t384 + (t388 * (-pkin(9) * t267 + t485) + t392 * (-pkin(2) * t267 + t220) - pkin(1) * t267 + pkin(8) * (t268 * t392 - t288 * t388)) * t383, (t114 + pkin(1) * (t232 * t388 + t299 * t392)) * t384 + (t388 * t431 + pkin(8) * (t232 * t392 - t299 * t388) + t428 * (-t283 * t387 - t287 * t391)) * t383, (t120 + pkin(1) * (t155 * t388 - t278 * t392)) * t384 + (pkin(8) * (t155 * t392 + t278 * t388) - t428 * t431) * t383, t539, t384 * t83 + (t388 * (-t136 * t387 + t140 * t391) - t545) * t383, -t585 + (-t392 * t207 + t574) * t383, t538, -t570, t537, (t48 - t586) * t384 + (t388 * (-t109 * t387 + t134 * t391 - t579) + t392 * (-t425 - t581) + t584) * t383, (t57 + pkin(1) * (-t392 * t527 - t576)) * t384 + (t388 * (-t113 * t387 + t143 * t391 - t578) + t392 * (-t415 - t580) - t582 + pkin(8) * (t388 * t527 - t575)) * t383, (t23 + pkin(1) * (t388 * t87 - t544)) * t384 + (t388 * (-pkin(9) * t85 - t387 * t45 + t391 * t50) + t392 * (-pkin(2) * t85 - t506) - pkin(1) * t85 + pkin(8) * (t392 * t87 + t550)) * t383, (t20 + pkin(1) * (-t201 * t392 + t36 * t388)) * t384 + (t388 * (-pkin(9) * t35 - pkin(10) * t497 - t387 * t62) + t392 * (-pkin(2) * t35 - t509) - pkin(1) * t35 + pkin(8) * (t201 * t388 + t36 * t392)) * t383, t537, t585 + (t392 * t209 - t574) * t383, t570, t539, t384 * t84 + (t388 * (-t137 * t387 + t141 * t391) - t545) * t383, t538, (t19 + pkin(1) * (t388 * t86 - t544)) * t384 + (t388 * (-pkin(9) * t82 - t37 * t387 + t391 * t40) + t392 * (-pkin(2) * t82 - t436) - pkin(1) * t82 + pkin(8) * (t392 * t86 + t550)) * t383, (t28 + t586) * t384 + (t388 * (-t387 * t58 + t391 * t64 + t579) + t392 * (-t402 + t581) - t584) * t383, (t29 + pkin(1) * (t392 * t395 + t576)) * t384 + (t388 * (-t387 * t61 + t391 * t65 + t578) + t392 * (-t397 + t580) + t582 + pkin(8) * (-t388 * t395 + t575)) * t383, (t9 + pkin(1) * (-t108 * t392 + t31 * t388)) * t384 + (t388 * (-pkin(9) * t30 - t27 * t387 + t34 * t391) + t392 * (-pkin(2) * t30 - t439) - pkin(1) * t30 + pkin(8) * (t108 * t388 + t31 * t392)) * t383, t384 * t75 + (t388 * (-t116 * t387 + t118 * t391) - t392 * t159) * t383, t384 * t49 + (t388 * (-t387 * t88 + t391 * t89) - t392 * t106) * t383, t384 * t60 + (t388 * (-t100 * t387 + t102 * t391) - t392 * t146) * t383, t384 * t74 + (t388 * (-t115 * t387 + t117 * t391) + t392 * t157) * t383, t384 * t59 + (t388 * (t101 * t391 - t387 * t99) - t392 * t147) * t383, t384 * t98 + (t388 * (-t148 * t387 + t149 * t391) - t392 * t186) * t383, (t5 + pkin(1) * (-t122 * t392 + t388 * t52)) * t384 + (t388 * (-pkin(9) * t51 - t14 * t387 + t16 * t391) + t392 * (-pkin(2) * t51 - t411) - pkin(1) * t51 + pkin(8) * (t122 * t388 + t392 * t52)) * t383, (t6 + pkin(1) * (-t129 * t392 + t388 * t54)) * t384 + (t388 * (-pkin(9) * t53 - t15 * t387 + t17 * t391) + t392 * (-pkin(2) * t53 - t416) - pkin(1) * t53 + pkin(8) * (t129 * t388 + t392 * t54)) * t383, (t4 + pkin(1) * (-t107 * t392 + t388 * t47)) * t384 + (t388 * (-pkin(9) * t46 - t11 * t387 + t13 * t391) + t392 * (-pkin(2) * t46 - t410) - pkin(1) * t46 + pkin(8) * (t107 * t388 + t392 * t47)) * t383, (t1 + pkin(1) * (-t26 * t392 + t388 * t8)) * t384 + (t388 * (-pkin(9) * t7 - t2 * t387 + t3 * t391) + t392 * (-pkin(2) * t7 - t423) - pkin(1) * t7 + pkin(8) * (t26 * t388 + t392 * t8)) * t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t370, t359, t336, t370, t337, t437, -t317, -t318, 0, 0, t282, t230, t270, t281, t271, t293, t176, t180, t114, t120, t518, t83, -t163, t519, -t165, t520, t48, t57, t23, t20, t520, t163, t165, t518, t84, t519, t19, t28, t29, t9, t75, t49, t60, t74, t59, t98, t5, t6, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t479, t331, t287, t479, -t283, -t352, -t219, -t220, 0, 0, t274, t523, t207, -t274, -t530, -t351, t425, t415, t506, t509, -t351, -t209, t530, t274, t523, -t274, t436, t402, t397, t439, t159, t106, t146, -t157, t147, t186, t411, t416, t410, t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, t523, t207, -t274, -t530, -t351, -t111, -t112, 0, 0, -t351, -t209, t530, t274, t523, -t274, t463, t405, t400, t501, t159, t106, t146, -t157, t147, t186, t419, t430, t414, t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -t409, t524, t97, 0, 0, 0, 0, 0, 0, t121, t128, t105, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t247, t171, -t248, t413, t240, -t43, -t44, 0, 0;];
tauJ_reg  = t24;