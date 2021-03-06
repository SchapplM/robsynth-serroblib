% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRRRR2
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 08:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:28:22
% EndTime: 2019-05-08 08:29:13
% DurationCPUTime: 18.08s
% Computational Cost: add. (162769->694), mult. (351410->998), div. (0->0), fcn. (271938->12), ass. (0->434)
t413 = cos(qJ(2));
t395 = t413 * qJDD(1);
t407 = sin(qJ(2));
t475 = qJD(1) * qJD(2);
t463 = t407 * t475;
t382 = t395 - t463;
t402 = t413 ^ 2;
t415 = qJD(1) ^ 2;
t408 = sin(qJ(1));
t518 = cos(qJ(1));
t461 = g(1) * t408 - t518 * g(2);
t431 = qJDD(1) * pkin(1) + t461;
t479 = qJD(1) * t407;
t432 = qJD(2) * pkin(2) - pkin(8) * t479;
t342 = pkin(2) * t382 - t432 * t479 + t431 + (pkin(8) * t402 + pkin(7)) * t415;
t403 = sin(qJ(6));
t406 = sin(qJ(3));
t412 = cos(qJ(3));
t372 = qJD(1) * t412 * t413 - t406 * t479;
t373 = (t406 * t413 + t407 * t412) * qJD(1);
t405 = sin(qJ(4));
t411 = cos(qJ(4));
t355 = t372 * t405 + t373 * t411;
t400 = qJD(2) + qJD(3);
t393 = qJD(4) + t400;
t404 = sin(qJ(5));
t410 = cos(qJ(5));
t335 = t355 * t404 - t410 * t393;
t337 = t355 * t410 + t393 * t404;
t409 = cos(qJ(6));
t301 = t409 * t335 + t337 * t403;
t303 = -t335 * t403 + t337 * t409;
t252 = t303 * t301;
t394 = t407 * qJDD(1);
t462 = t413 * t475;
t381 = t394 + t462;
t455 = t406 * t381 - t412 * t382;
t339 = -qJD(3) * t373 - t455;
t340 = t372 * qJD(3) + t412 * t381 + t406 * t382;
t457 = -t411 * t339 + t340 * t405;
t282 = -qJD(4) * t355 - t457;
t281 = qJDD(5) - t282;
t280 = qJDD(6) + t281;
t526 = -t252 + t280;
t534 = t403 * t526;
t305 = t337 * t335;
t524 = t281 - t305;
t533 = t404 * t524;
t353 = -t411 * t372 + t373 * t405;
t320 = t355 * t353;
t473 = qJDD(2) + qJDD(3);
t392 = qJDD(4) + t473;
t523 = -t320 + t392;
t532 = t405 * t523;
t361 = t372 * t373;
t522 = t361 + t473;
t531 = t406 * t522;
t530 = t409 * t526;
t529 = t410 * t524;
t528 = t411 * t523;
t527 = t412 * t522;
t370 = t372 ^ 2;
t444 = pkin(3) * t400 - pkin(9) * t373;
t274 = pkin(3) * t339 + pkin(9) * t370 - t373 * t444 + t342;
t433 = g(1) * t518 + t408 * g(2);
t511 = qJDD(1) * pkin(7);
t375 = -t415 * pkin(1) - t433 + t511;
t397 = t402 * t415;
t515 = t407 * g(3);
t419 = -pkin(2) * t397 + t382 * pkin(8) - qJD(2) * t432 - t515;
t481 = t407 * t415;
t425 = pkin(2) * t481 + pkin(8) * t475 - g(3);
t487 = t375 * t407;
t426 = qJDD(2) * pkin(2) - pkin(8) * t381 - t487;
t298 = t412 * t419 + t406 * t426 + (t412 * t375 + t406 * t425) * t413;
t255 = -t370 * pkin(3) + t339 * pkin(9) - t400 * t444 + t298;
t297 = t406 * t419 - t412 * t426 + (t406 * t375 - t412 * t425) * t413;
t488 = t372 * t400;
t325 = -t340 + t488;
t417 = pkin(3) * t522 + pkin(9) * t325 - t297;
t197 = t411 * t255 + t405 * t417;
t316 = pkin(4) * t353 - pkin(10) * t355;
t520 = t393 ^ 2;
t181 = -pkin(4) * t520 + t392 * pkin(10) - t353 * t316 + t197;
t283 = -t353 * qJD(4) + t405 * t339 + t411 * t340;
t491 = t353 * t393;
t435 = t283 - t491;
t188 = -t435 * pkin(10) + (t355 * t393 - t282) * pkin(4) - t274;
t122 = t181 * t404 - t410 * t188;
t123 = t410 * t181 + t404 * t188;
t71 = t404 * t122 + t410 * t123;
t437 = -t283 * t410 - t392 * t404;
t256 = -qJD(5) * t335 - t437;
t349 = qJD(5) + t353;
t311 = t349 * t335;
t239 = t256 + t311;
t458 = t283 * t404 - t410 * t392;
t429 = qJD(5) * t337 + t458;
t187 = -t301 * qJD(6) + t409 * t256 - t403 * t429;
t346 = qJD(6) + t349;
t278 = t346 * t301;
t525 = -t278 + t187;
t414 = qJD(2) ^ 2;
t449 = t397 + t414;
t459 = t256 * t403 + t409 * t429;
t158 = (qJD(6) - t346) * t303 + t459;
t235 = (qJD(5) - t349) * t337 + t458;
t299 = t301 ^ 2;
t300 = t303 ^ 2;
t521 = t335 ^ 2;
t334 = t337 ^ 2;
t345 = t346 ^ 2;
t348 = t349 ^ 2;
t351 = t353 ^ 2;
t352 = t355 ^ 2;
t371 = t373 ^ 2;
t399 = t400 ^ 2;
t86 = pkin(5) * t524 - pkin(11) * t239 - t122;
t306 = pkin(5) * t349 - pkin(11) * t337;
t98 = -pkin(5) * t521 - pkin(11) * t429 - t349 * t306 + t123;
t49 = t403 * t98 - t409 * t86;
t50 = t403 * t86 + t409 * t98;
t28 = t403 * t50 - t409 * t49;
t519 = pkin(5) * t28;
t517 = pkin(4) * t405;
t161 = t278 + t187;
t110 = -t158 * t403 - t161 * t409;
t516 = pkin(5) * t110;
t514 = t28 * t404;
t513 = t28 * t410;
t196 = t255 * t405 - t411 * t417;
t180 = -t392 * pkin(4) - t520 * pkin(10) + t316 * t355 + t196;
t512 = -pkin(4) * t180 + pkin(10) * t71;
t130 = pkin(5) * t429 - pkin(11) * t521 + t306 * t337 + t180;
t510 = t130 * t403;
t509 = t130 * t409;
t134 = -t196 * t411 + t197 * t405;
t508 = t134 * t406;
t507 = t134 * t412;
t214 = t252 + t280;
t506 = t214 * t403;
t505 = t214 * t409;
t246 = t281 + t305;
t504 = t246 * t404;
t503 = t246 * t410;
t248 = -t297 * t412 + t298 * t406;
t502 = t248 * t407;
t501 = t274 * t405;
t500 = t274 * t411;
t314 = t320 + t392;
t499 = t314 * t405;
t498 = t314 * t411;
t497 = t342 * t406;
t496 = t342 * t412;
t495 = t346 * t403;
t494 = t346 * t409;
t493 = t349 * t404;
t492 = t349 * t410;
t358 = -t361 + t473;
t490 = t358 * t406;
t489 = t358 * t412;
t486 = t393 * t405;
t485 = t393 * t411;
t484 = t400 * t406;
t483 = t400 * t412;
t176 = t404 * t180;
t389 = t413 * t481;
t474 = qJDD(2) + t389;
t482 = t407 * t474;
t177 = t410 * t180;
t480 = t413 * (qJDD(2) - t389);
t477 = qJD(5) + t349;
t59 = -t180 * t411 + t405 * t71;
t472 = pkin(3) * t59 + t512;
t471 = t405 * t252;
t470 = t411 * t252;
t469 = t405 * t305;
t468 = t411 * t305;
t294 = -t334 - t348;
t207 = -t294 * t404 - t503;
t240 = t335 * t477 + t437;
t467 = pkin(4) * t240 + pkin(10) * t207 + t176;
t286 = -t348 - t521;
t204 = t286 * t410 - t533;
t236 = -t337 * t477 - t458;
t466 = pkin(4) * t236 + pkin(10) * t204 - t177;
t465 = -pkin(4) * t411 - pkin(3);
t29 = t403 * t49 + t409 * t50;
t112 = -t158 * t409 + t161 * t403;
t226 = -t299 - t300;
t22 = -pkin(5) * t226 + pkin(11) * t112 + t29;
t24 = -pkin(11) * t110 - t28;
t65 = -t110 * t404 + t112 * t410;
t460 = -pkin(4) * t226 + pkin(10) * t65 + t410 * t22 + t404 * t24;
t135 = t196 * t405 + t411 * t197;
t249 = t297 * t406 + t412 * t298;
t364 = g(3) * t413 + t487;
t442 = -t413 * t375 + t515;
t456 = t407 * t364 - t413 * t442;
t244 = -t345 - t299;
t167 = t244 * t403 + t530;
t168 = t244 * t409 - t534;
t116 = -t167 * t404 + t168 * t410;
t157 = (qJD(6) + t346) * t303 + t459;
t67 = -pkin(5) * t157 + pkin(11) * t168 - t509;
t88 = -pkin(11) * t167 + t510;
t454 = -pkin(4) * t157 + pkin(10) * t116 + t404 * t88 + t410 * t67;
t261 = -t300 - t345;
t182 = t261 * t409 - t506;
t183 = -t261 * t403 - t505;
t119 = -t182 * t404 + t183 * t410;
t72 = -pkin(5) * t525 + pkin(11) * t183 + t510;
t94 = -pkin(11) * t182 + t509;
t453 = -pkin(4) * t525 + pkin(10) * t119 + t404 * t94 + t410 * t72;
t174 = -t235 * t410 + t239 * t404;
t275 = t334 + t521;
t452 = pkin(4) * t275 + pkin(10) * t174 + t71;
t149 = t207 * t405 + t240 * t411;
t451 = pkin(3) * t149 + t467;
t146 = t204 * t405 + t236 * t411;
t450 = pkin(3) * t146 + t466;
t53 = -t226 * t411 + t405 * t65;
t448 = pkin(3) * t53 + t460;
t74 = t116 * t405 - t157 * t411;
t446 = pkin(3) * t74 + t454;
t77 = t119 * t405 - t411 * t525;
t445 = pkin(3) * t77 + t453;
t312 = -t520 - t351;
t271 = t312 * t405 + t528;
t443 = pkin(3) * t271 - t196;
t441 = -t381 + t462;
t440 = t382 + t463;
t141 = t174 * t405 + t275 * t411;
t438 = pkin(3) * t141 + t452;
t70 = -t122 * t410 + t123 * t404;
t436 = pkin(5) * t167 - t49;
t434 = t340 + t488;
t17 = t29 * t410 - t514;
t26 = -pkin(5) * t130 + pkin(11) * t29;
t430 = -pkin(4) * t130 + pkin(10) * t17 - pkin(11) * t514 + t410 * t26;
t11 = -t130 * t411 + t17 * t405;
t428 = pkin(3) * t11 + t430;
t427 = pkin(5) * t182 - t50;
t424 = (-qJD(4) + t393) * t355 - t457;
t423 = (-qJD(3) + t400) * t373 - t455;
t338 = -t352 - t520;
t288 = t338 * t411 - t499;
t416 = pkin(3) * t288 - t197;
t401 = t407 ^ 2;
t396 = t401 * t415;
t383 = t395 - 0.2e1 * t463;
t380 = t394 + 0.2e1 * t462;
t374 = pkin(7) * t415 + t431;
t367 = -t371 + t399;
t366 = t370 - t399;
t363 = -t371 - t399;
t360 = t371 - t370;
t356 = -t399 - t370;
t344 = -t352 + t520;
t343 = t351 - t520;
t341 = -t370 - t371;
t328 = -t363 * t406 - t489;
t327 = t363 * t412 - t490;
t321 = (qJD(3) + t400) * t373 + t455;
t319 = t352 - t351;
t318 = t356 * t412 - t531;
t317 = t356 * t406 + t527;
t310 = -t334 + t348;
t309 = -t348 + t521;
t308 = (-t353 * t411 + t355 * t405) * t393;
t307 = (-t353 * t405 - t355 * t411) * t393;
t304 = t334 - t521;
t296 = -t351 - t352;
t293 = t343 * t411 - t499;
t292 = -t344 * t405 + t528;
t291 = t343 * t405 + t498;
t290 = t344 * t411 + t532;
t289 = -t338 * t405 - t498;
t285 = -t325 * t406 + t412 * t423;
t284 = t325 * t412 + t406 * t423;
t277 = -t300 + t345;
t276 = t299 - t345;
t272 = t312 * t411 - t532;
t269 = (-t335 * t404 - t337 * t410) * t349;
t268 = (-t335 * t410 + t337 * t404) * t349;
t266 = -t283 - t491;
t262 = (qJD(4) + t393) * t355 + t457;
t260 = t283 * t411 - t355 * t486;
t259 = t283 * t405 + t355 * t485;
t258 = -t282 * t405 + t353 * t485;
t257 = t282 * t411 + t353 * t486;
t251 = t300 - t299;
t243 = -t288 * t406 + t289 * t412;
t242 = t288 * t412 + t289 * t406;
t241 = -pkin(9) * t288 - t500;
t238 = t256 - t311;
t232 = (-t301 * t409 + t303 * t403) * t346;
t231 = (-t301 * t403 - t303 * t409) * t346;
t230 = t256 * t410 - t337 * t493;
t229 = t335 * t493 - t410 * t429;
t228 = t256 * t404 + t337 * t492;
t227 = t335 * t492 + t404 * t429;
t224 = -pkin(9) * t271 - t501;
t223 = t268 * t411 + t281 * t405;
t222 = t268 * t405 - t281 * t411;
t221 = t310 * t410 + t533;
t220 = t309 * t410 - t504;
t219 = t309 * t404 + t503;
t218 = -t310 * t404 + t529;
t217 = -t271 * t406 + t272 * t412;
t216 = t271 * t412 + t272 * t406;
t212 = -t266 * t405 + t411 * t424;
t211 = -t262 * t411 - t405 * t435;
t210 = t266 * t411 + t405 * t424;
t209 = -t262 * t405 + t411 * t435;
t208 = pkin(3) * t210;
t206 = t294 * t410 - t504;
t203 = t286 * t404 + t529;
t201 = t230 * t411 + t469;
t200 = t227 * t411 - t469;
t199 = t230 * t405 - t468;
t198 = t227 * t405 + t468;
t194 = -pkin(3) * t435 + pkin(9) * t289 - t501;
t193 = -pkin(3) * t262 + pkin(9) * t272 + t500;
t192 = t276 * t409 - t506;
t191 = -t277 * t403 + t530;
t190 = t276 * t403 + t505;
t189 = t277 * t409 + t534;
t186 = -qJD(6) * t303 - t459;
t175 = t236 * t404 + t238 * t410;
t173 = t236 * t410 - t238 * t404;
t172 = -t235 * t404 - t239 * t410;
t170 = -t231 * t404 + t232 * t410;
t169 = t231 * t410 + t232 * t404;
t166 = t220 * t411 - t235 * t405;
t165 = t218 * t411 + t239 * t405;
t164 = t220 * t405 + t235 * t411;
t163 = t218 * t405 - t239 * t411;
t154 = t187 * t409 - t303 * t495;
t153 = t187 * t403 + t303 * t494;
t152 = -t186 * t403 + t301 * t494;
t151 = t186 * t409 + t301 * t495;
t150 = t207 * t411 - t240 * t405;
t147 = t204 * t411 - t236 * t405;
t144 = t173 * t411 + t304 * t405;
t143 = t173 * t405 - t304 * t411;
t142 = t174 * t411 - t275 * t405;
t139 = t170 * t411 + t280 * t405;
t138 = t170 * t405 - t280 * t411;
t137 = -t210 * t406 + t212 * t412;
t136 = t210 * t412 + t212 * t406;
t133 = pkin(3) * t134;
t132 = -pkin(10) * t206 + t177;
t131 = -pkin(10) * t203 + t176;
t128 = pkin(3) * t274 + pkin(9) * t135;
t127 = -t190 * t404 + t192 * t410;
t126 = -t189 * t404 + t191 * t410;
t125 = t190 * t410 + t192 * t404;
t124 = t189 * t410 + t191 * t404;
t118 = t182 * t410 + t183 * t404;
t115 = t167 * t410 + t168 * t404;
t113 = -pkin(9) * t210 - t134;
t111 = -t157 * t409 - t403 * t525;
t109 = -t157 * t403 + t409 * t525;
t108 = -t153 * t404 + t154 * t410;
t107 = t153 * t410 + t154 * t404;
t106 = -t151 * t404 + t152 * t410;
t105 = t151 * t410 + t152 * t404;
t104 = -pkin(3) * t296 + pkin(9) * t212 + t135;
t103 = -t149 * t406 + t150 * t412;
t102 = t149 * t412 + t150 * t406;
t101 = -t146 * t406 + t147 * t412;
t100 = t146 * t412 + t147 * t406;
t99 = -pkin(4) * t206 + t123;
t97 = -pkin(4) * t203 + t122;
t96 = -t141 * t406 + t142 * t412;
t95 = t141 * t412 + t142 * t406;
t92 = t108 * t411 + t471;
t91 = t106 * t411 - t471;
t90 = t108 * t405 - t470;
t89 = t106 * t405 + t470;
t84 = t127 * t411 - t158 * t405;
t83 = t126 * t411 + t161 * t405;
t82 = t127 * t405 + t158 * t411;
t81 = t126 * t405 - t161 * t411;
t80 = t135 * t412 - t508;
t79 = t135 * t406 + t507;
t78 = t119 * t411 + t405 * t525;
t75 = t116 * t411 + t157 * t405;
t64 = -t109 * t404 + t111 * t410;
t63 = t109 * t410 + t111 * t404;
t62 = t110 * t410 + t112 * t404;
t60 = t180 * t405 + t411 * t71;
t57 = -pkin(10) * t172 - t70;
t56 = t251 * t405 + t411 * t64;
t55 = -t251 * t411 + t405 * t64;
t54 = t226 * t405 + t411 * t65;
t51 = -pkin(9) * t149 + t132 * t411 - t405 * t99;
t47 = -pkin(9) * t146 + t131 * t411 - t405 * t97;
t46 = -pkin(3) * t206 + pkin(9) * t150 + t132 * t405 + t411 * t99;
t45 = -pkin(3) * t203 + pkin(9) * t147 + t131 * t405 + t411 * t97;
t44 = -t406 * t77 + t412 * t78;
t43 = t406 * t78 + t412 * t77;
t42 = -t406 * t74 + t412 * t75;
t41 = t406 * t75 + t412 * t74;
t40 = -pkin(4) * t62 - t516;
t39 = -pkin(9) * t141 + t172 * t517 + t411 * t57;
t38 = pkin(9) * t142 + t172 * t465 + t405 * t57;
t37 = -pkin(4) * t118 - t427;
t36 = -pkin(10) * t118 - t404 * t72 + t410 * t94;
t35 = -pkin(4) * t115 - t436;
t34 = -pkin(10) * t115 - t404 * t67 + t410 * t88;
t33 = -t406 * t59 + t412 * t60;
t32 = t406 * t60 + t412 * t59;
t31 = -t406 * t53 + t412 * t54;
t30 = t406 * t54 + t412 * t53;
t27 = -pkin(9) * t59 + (-pkin(10) * t411 + t517) * t70;
t20 = pkin(9) * t60 + (-pkin(10) * t405 + t465) * t70;
t19 = -pkin(9) * t77 + t36 * t411 - t37 * t405;
t18 = -pkin(9) * t74 + t34 * t411 - t35 * t405;
t16 = t29 * t404 + t513;
t14 = -pkin(3) * t118 + pkin(9) * t78 + t36 * t405 + t37 * t411;
t13 = -pkin(3) * t115 + pkin(9) * t75 + t34 * t405 + t35 * t411;
t12 = t130 * t405 + t17 * t411;
t9 = -pkin(10) * t62 - t22 * t404 + t24 * t410;
t8 = -pkin(4) * t16 - t519;
t7 = -pkin(9) * t53 - t40 * t405 + t411 * t9;
t6 = -pkin(10) * t16 - pkin(11) * t513 - t26 * t404;
t5 = -t11 * t406 + t12 * t412;
t4 = t11 * t412 + t12 * t406;
t3 = -pkin(3) * t62 + pkin(9) * t54 + t40 * t411 + t405 * t9;
t2 = -pkin(9) * t11 - t405 * t8 + t411 * t6;
t1 = -pkin(3) * t16 + pkin(9) * t12 + t405 * t6 + t411 * t8;
t10 = [0, 0, 0, 0, 0, qJDD(1), t461, t433, 0, 0, (t381 + t462) * t407, t380 * t413 + t383 * t407, t482 + t413 * (-t396 + t414), (t382 - t463) * t413, t407 * (t397 - t414) + t480, 0, t413 * t374 + pkin(1) * t383 + pkin(7) * (-t413 * t449 - t482), -t407 * t374 - pkin(1) * t380 + pkin(7) * (-t480 - t407 * (-t396 - t414)), pkin(1) * (t396 + t397) + (t401 + t402) * t511 + t456, pkin(1) * t374 + pkin(7) * t456, t407 * (t340 * t412 - t373 * t484) + t413 * (t340 * t406 + t373 * t483), t407 * (-t321 * t412 - t406 * t434) + t413 * (-t321 * t406 + t412 * t434), t407 * (-t367 * t406 + t527) + t413 * (t367 * t412 + t531), t407 * (-t339 * t406 - t372 * t483) + t413 * (t339 * t412 - t372 * t484), t407 * (t366 * t412 - t490) + t413 * (t366 * t406 + t489), (t407 * (t372 * t412 + t373 * t406) + t413 * (t372 * t406 - t373 * t412)) * t400, t407 * (-pkin(8) * t317 - t497) + t413 * (-pkin(2) * t321 + pkin(8) * t318 + t496) - pkin(1) * t321 + pkin(7) * (-t317 * t407 + t318 * t413), t407 * (-pkin(8) * t327 - t496) + t413 * (-pkin(2) * t434 + pkin(8) * t328 - t497) - pkin(1) * t434 + pkin(7) * (-t327 * t407 + t328 * t413), t407 * (-pkin(8) * t284 - t248) + t413 * (-pkin(2) * t341 + pkin(8) * t285 + t249) - pkin(1) * t341 + pkin(7) * (-t284 * t407 + t285 * t413), -pkin(8) * t502 + t413 * (pkin(2) * t342 + pkin(8) * t249) + pkin(1) * t342 + pkin(7) * (t249 * t413 - t502), t407 * (-t259 * t406 + t260 * t412) + t413 * (t259 * t412 + t260 * t406), t407 * (-t209 * t406 + t211 * t412) + t413 * (t209 * t412 + t211 * t406), t407 * (-t290 * t406 + t292 * t412) + t413 * (t290 * t412 + t292 * t406), t407 * (-t257 * t406 + t258 * t412) + t413 * (t257 * t412 + t258 * t406), t407 * (-t291 * t406 + t293 * t412) + t413 * (t291 * t412 + t293 * t406), t407 * (-t307 * t406 + t308 * t412) + t413 * (t307 * t412 + t308 * t406), t407 * (-pkin(8) * t216 - t193 * t406 + t224 * t412) + t413 * (-pkin(2) * t262 + pkin(8) * t217 + t193 * t412 + t224 * t406) - pkin(1) * t262 + pkin(7) * (-t216 * t407 + t217 * t413), t407 * (-pkin(8) * t242 - t194 * t406 + t241 * t412) + t413 * (-pkin(2) * t435 + pkin(8) * t243 + t194 * t412 + t241 * t406) - pkin(1) * t435 + pkin(7) * (-t242 * t407 + t243 * t413), t407 * (-pkin(8) * t136 - t104 * t406 + t113 * t412) + t413 * (-pkin(2) * t296 + pkin(8) * t137 + t104 * t412 + t113 * t406) - pkin(1) * t296 + pkin(7) * (-t136 * t407 + t137 * t413), t407 * (-pkin(8) * t79 - pkin(9) * t507 - t128 * t406) + t413 * (pkin(2) * t274 + pkin(8) * t80 - pkin(9) * t508 + t128 * t412) + pkin(1) * t274 + pkin(7) * (-t407 * t79 + t413 * t80), t407 * (-t199 * t406 + t201 * t412) + t413 * (t199 * t412 + t201 * t406), t407 * (-t143 * t406 + t144 * t412) + t413 * (t143 * t412 + t144 * t406), t407 * (-t163 * t406 + t165 * t412) + t413 * (t163 * t412 + t165 * t406), t407 * (-t198 * t406 + t200 * t412) + t413 * (t198 * t412 + t200 * t406), t407 * (-t164 * t406 + t166 * t412) + t413 * (t164 * t412 + t166 * t406), t407 * (-t222 * t406 + t223 * t412) + t413 * (t222 * t412 + t223 * t406), t407 * (-pkin(8) * t100 - t406 * t45 + t412 * t47) + t413 * (-pkin(2) * t203 + pkin(8) * t101 + t406 * t47 + t412 * t45) - pkin(1) * t203 + pkin(7) * (-t100 * t407 + t101 * t413), t407 * (-pkin(8) * t102 - t406 * t46 + t412 * t51) + t413 * (-pkin(2) * t206 + pkin(8) * t103 + t406 * t51 + t412 * t46) - pkin(1) * t206 + pkin(7) * (-t102 * t407 + t103 * t413), t407 * (-pkin(8) * t95 - t38 * t406 + t39 * t412) + t413 * (-pkin(2) * t172 + pkin(8) * t96 + t38 * t412 + t39 * t406) - pkin(1) * t172 + pkin(7) * (-t407 * t95 + t413 * t96), t407 * (-pkin(8) * t32 - t20 * t406 + t27 * t412) + t413 * (-pkin(2) * t70 + pkin(8) * t33 + t20 * t412 + t27 * t406) - pkin(1) * t70 + pkin(7) * (-t32 * t407 + t33 * t413), t407 * (-t406 * t90 + t412 * t92) + t413 * (t406 * t92 + t412 * t90), t407 * (-t406 * t55 + t412 * t56) + t413 * (t406 * t56 + t412 * t55), t407 * (-t406 * t81 + t412 * t83) + t413 * (t406 * t83 + t412 * t81), t407 * (-t406 * t89 + t412 * t91) + t413 * (t406 * t91 + t412 * t89), t407 * (-t406 * t82 + t412 * t84) + t413 * (t406 * t84 + t412 * t82), t407 * (-t138 * t406 + t139 * t412) + t413 * (t138 * t412 + t139 * t406), t407 * (-pkin(8) * t41 - t13 * t406 + t18 * t412) + t413 * (-pkin(2) * t115 + pkin(8) * t42 + t13 * t412 + t18 * t406) - pkin(1) * t115 + pkin(7) * (-t407 * t41 + t413 * t42), t407 * (-pkin(8) * t43 - t14 * t406 + t19 * t412) + t413 * (-pkin(2) * t118 + pkin(8) * t44 + t14 * t412 + t19 * t406) - pkin(1) * t118 + pkin(7) * (-t407 * t43 + t413 * t44), t407 * (-pkin(8) * t30 - t3 * t406 + t412 * t7) + t413 * (-pkin(2) * t62 + pkin(8) * t31 + t3 * t412 + t406 * t7) - pkin(1) * t62 + pkin(7) * (-t30 * t407 + t31 * t413), t407 * (-pkin(8) * t4 - t1 * t406 + t2 * t412) + t413 * (-pkin(2) * t16 + pkin(8) * t5 + t1 * t412 + t2 * t406) - pkin(1) * t16 + pkin(7) * (-t4 * t407 + t413 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t389, -t397 + t396, t394, t389, t395, qJDD(2), -t364, t442, 0, 0, -t361, t360, -t325, t361, t423, t473, t406 * t442 - t412 * t364 + (-t406 * t440 + t412 * t441) * pkin(8) + (t406 * t449 + t412 * t474 + t317) * pkin(2), t412 * t442 + t406 * t364 + (-t406 * t441 - t412 * t440) * pkin(8) + (-t406 * t474 + t412 * t449 + t327) * pkin(2), pkin(2) * t284, pkin(2) * t248, t320, t319, -t266, -t320, t424, t392, pkin(2) * t216 + t443, pkin(2) * t242 + t416, pkin(2) * t136 + t208, pkin(2) * t79 + t133, t228, t175, t221, t229, t219, t269, pkin(2) * t100 + t450, pkin(2) * t102 + t451, pkin(2) * t95 + t438, pkin(2) * t32 + t472, t107, t63, t124, t105, t125, t169, pkin(2) * t41 + t446, pkin(2) * t43 + t445, pkin(2) * t30 + t448, pkin(2) * t4 + t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361, t360, -t325, t361, t423, t473, -t297, -t298, 0, 0, t320, t319, -t266, -t320, t424, t392, t443, t416, t208, t133, t228, t175, t221, t229, t219, t269, t450, t451, t438, t472, t107, t63, t124, t105, t125, t169, t446, t445, t448, t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, t319, -t266, -t320, t424, t392, -t196, -t197, 0, 0, t228, t175, t221, t229, t219, t269, t466, t467, t452, t512, t107, t63, t124, t105, t125, t169, t454, t453, t460, t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, t304, t239, -t305, -t235, t281, -t122, -t123, 0, 0, t252, t251, t161, -t252, -t158, t280, t436, t427, t516, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t251, t161, -t252, -t158, t280, -t49, -t50, 0, 0;];
tauJ_reg  = t10;
