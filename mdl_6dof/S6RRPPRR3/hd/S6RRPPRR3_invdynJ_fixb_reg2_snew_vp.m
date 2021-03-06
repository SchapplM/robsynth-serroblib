% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 10:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:04:20
% EndTime: 2019-05-06 10:04:53
% DurationCPUTime: 19.05s
% Computational Cost: add. (152373->775), mult. (413150->1158), div. (0->0), fcn. (333656->14), ass. (0->470)
t388 = sin(pkin(12));
t389 = sin(pkin(11));
t392 = cos(pkin(11));
t390 = sin(pkin(6));
t399 = cos(qJ(2));
t470 = qJD(1) * t399;
t448 = t390 * t470;
t396 = sin(qJ(2));
t471 = qJD(1) * t396;
t449 = t390 * t471;
t359 = t389 * t448 + t392 * t449;
t393 = cos(pkin(6));
t472 = qJD(1) * t393;
t383 = qJD(2) + t472;
t391 = cos(pkin(12));
t339 = -t388 * t359 + t383 * t391;
t340 = t359 * t391 + t383 * t388;
t301 = t339 * t340;
t445 = qJD(1) * qJD(2) * t390;
t464 = qJDD(1) * t390;
t365 = t396 * t464 + t399 * t445;
t366 = -t396 * t445 + t399 * t464;
t440 = t365 * t389 - t392 * t366;
t530 = t301 + t440;
t542 = t388 * t530;
t357 = t389 * t449 - t392 * t448;
t331 = t359 * t357;
t438 = qJDD(1) * t393 + qJDD(2);
t528 = -t331 + t438;
t541 = t389 * t528;
t540 = t391 * t530;
t539 = t392 * t528;
t394 = sin(qJ(6));
t395 = sin(qJ(5));
t398 = cos(qJ(5));
t299 = t339 * t395 + t340 * t398;
t397 = cos(qJ(6));
t428 = qJD(5) + t357;
t280 = t299 * t394 - t397 * t428;
t282 = t397 * t299 + t394 * t428;
t232 = t282 * t280;
t332 = t392 * t365 + t389 * t366;
t317 = t391 * t332 + t388 * t438;
t441 = t332 * t388 - t391 * t438;
t442 = t395 * t317 + t398 * t441;
t243 = -t299 * qJD(5) - t442;
t237 = qJDD(6) - t243;
t532 = -t232 + t237;
t538 = t394 * t532;
t297 = -t398 * t339 + t340 * t395;
t255 = t299 * t297;
t329 = qJDD(5) + t440;
t531 = -t255 + t329;
t537 = t395 * t531;
t536 = t397 * t532;
t535 = t398 * t531;
t244 = -qJD(5) * t297 + t317 * t398 - t395 * t441;
t286 = t428 * t297;
t534 = t244 - t286;
t425 = t428 ^ 2;
t385 = t390 ^ 2;
t533 = t385 * (-t383 + t472);
t314 = t357 * t339;
t270 = t314 - t317;
t529 = t314 + t317;
t346 = t383 * t357;
t527 = -t346 + t332;
t369 = t383 ^ 2;
t479 = t390 * t396;
t384 = g(3) * t479;
t517 = sin(qJ(1));
t518 = cos(qJ(1));
t429 = g(1) * t517 - g(2) * t518;
t418 = qJDD(1) * pkin(1) + t429;
t416 = t393 * t418;
t430 = g(1) * t518 + g(2) * t517;
t417 = pkin(8) * t464 - t430;
t482 = t383 * t390;
t454 = qJ(3) * t482;
t512 = pkin(8) * t390;
t463 = t393 * t512;
t387 = t399 ^ 2;
t480 = t385 * t387;
t283 = t399 * t417 + t396 * t416 - t384 - pkin(2) * t369 + t366 * qJ(3) + (t396 * t454 + (-pkin(1) * t399 - pkin(2) * t480 + t396 * t463) * qJD(1)) * qJD(1);
t521 = qJD(1) ^ 2;
t415 = -pkin(1) * t521 + t417;
t412 = t396 * t415;
t410 = pkin(2) * t438 - t365 * qJ(3) - t412;
t526 = t392 * t283 + t389 * t410;
t468 = t390 * t521;
t414 = pkin(8) * t468 + t418;
t413 = t393 * t414;
t509 = t390 * g(3);
t318 = t412 + (-t413 + t509) * t399;
t319 = t396 * t413 + t399 * t415 - t384;
t525 = t396 * t318 + t399 * t319;
t524 = qJDD(3) + (pkin(2) * t383 - qJ(3) * t449) * t449 - t366 * pkin(2);
t523 = t527 * qJ(4);
t293 = qJD(6) + t297;
t444 = t244 * t394 - t397 * t329;
t171 = (qJD(6) - t293) * t282 + t444;
t508 = t393 * g(3);
t294 = t390 * t418 + (pkin(8) * t385 + qJ(3) * t480) * t521 + t508 - t524;
t278 = t280 ^ 2;
t279 = t282 ^ 2;
t291 = t293 ^ 2;
t295 = t297 ^ 2;
t296 = t299 ^ 2;
t337 = t339 ^ 2;
t338 = t340 ^ 2;
t522 = t357 ^ 2;
t356 = t359 ^ 2;
t520 = 2 * qJD(3);
t519 = 2 * qJD(4);
t408 = (t416 - t509 + (t454 + (pkin(2) * t385 * t396 + t463) * qJD(1)) * qJD(1)) * t399;
t443 = t389 * t283 - t392 * (t408 + t410);
t223 = t359 * t520 + t443;
t407 = t389 * t408;
t402 = t407 + t526;
t462 = -0.2e1 * qJD(3) * t357;
t224 = t402 + t462;
t180 = -t223 * t392 + t224 * t389;
t516 = pkin(2) * t180;
t307 = t346 + t332;
t484 = t359 * t383;
t431 = -t440 + t484;
t256 = -t307 * t392 + t389 * t431;
t515 = pkin(2) * t256;
t514 = pkin(3) * t389;
t513 = pkin(5) * t395;
t510 = t369 * pkin(3);
t433 = t438 * qJ(4);
t324 = pkin(3) * t357 - qJ(4) * t359;
t467 = t520 + t324;
t401 = -t357 * t467 + t402 + t433 - t510;
t302 = t440 + t484;
t404 = t302 * pkin(3) - t294 - t523;
t164 = t339 * t519 + t388 * t404 + t391 * t401;
t309 = pkin(4) * t357 - pkin(9) * t340;
t139 = -pkin(4) * t337 - pkin(9) * t441 - t309 * t357 + t164;
t461 = t340 * t519;
t403 = -t391 * t404 + t461;
t400 = pkin(4) * t530 + pkin(9) * t270 - t388 * t401 - t403;
t81 = t139 * t395 - t398 * t400;
t82 = t398 * t139 + t395 * t400;
t50 = t395 * t82 - t398 * t81;
t507 = t388 * t50;
t506 = t391 * t50;
t253 = pkin(5) * t297 - pkin(10) * t299;
t76 = -t329 * pkin(5) - pkin(10) * t425 + t253 * t299 + t81;
t505 = t394 * t76;
t504 = t397 * t76;
t202 = -t438 * pkin(3) - t369 * qJ(4) + t359 * t467 + qJDD(4) + t443;
t179 = pkin(4) * t441 - t337 * pkin(9) + t340 * t309 + t202;
t503 = t179 * t395;
t502 = t179 * t398;
t185 = t232 + t237;
t501 = t185 * t394;
t500 = t185 * t397;
t499 = t202 * t388;
t498 = t202 * t391;
t241 = t255 + t329;
t497 = t241 * t395;
t496 = t241 * t398;
t272 = -t301 + t440;
t495 = t272 * t388;
t494 = t272 * t391;
t493 = t293 * t394;
t492 = t293 * t397;
t491 = t294 * t389;
t490 = t294 * t392;
t322 = t331 + t438;
t489 = t322 * t389;
t488 = t322 * t392;
t487 = t340 * t357;
t486 = t357 * t388;
t485 = t357 * t391;
t483 = t383 * t389;
t481 = t383 * t392;
t477 = t396 * t180;
t469 = t385 * t521;
t377 = t399 * t396 * t469;
t363 = t377 + t438;
t475 = t396 * t363;
t364 = -t377 + t438;
t473 = t399 * t364;
t465 = qJD(6) + t293;
t460 = t395 * t232;
t459 = t398 * t232;
t458 = t389 * t255;
t457 = t392 * t255;
t456 = t389 * t301;
t455 = t392 * t301;
t453 = t393 * t331;
t371 = t383 * t448;
t452 = t371 + t365;
t451 = -pkin(3) * t392 - pkin(2);
t450 = -pkin(5) * t398 - pkin(4);
t386 = t396 ^ 2;
t447 = t386 * t469;
t446 = t387 * t469;
t423 = t428 * t299;
t115 = -t534 * pkin(10) + (-t243 + t423) * pkin(5) + t179;
t77 = -pkin(5) * t425 + t329 * pkin(10) - t297 * t253 + t82;
t56 = -t397 * t115 + t394 * t77;
t57 = t115 * t394 + t397 * t77;
t38 = t394 * t56 + t397 * t57;
t51 = t395 * t81 + t398 * t82;
t405 = -t357 * t324 + t433 + t462 + t526;
t406 = t388 * t407;
t163 = t388 * (t405 - t510) + t406 + t403;
t106 = t163 * t388 + t391 * t164;
t181 = t223 * t389 + t392 * t224;
t437 = -pkin(5) * t76 + pkin(10) * t38;
t37 = t394 * t57 - t397 * t56;
t435 = t389 * t440;
t105 = -t163 * t391 + t164 * t388;
t434 = -t244 * t397 - t329 * t394;
t266 = t441 - t487;
t220 = -t291 - t278;
t141 = t220 * t397 - t538;
t172 = -t282 * t465 - t444;
t427 = pkin(5) * t172 + pkin(10) * t141 - t504;
t221 = -t279 - t291;
t143 = -t221 * t394 - t500;
t176 = t280 * t465 + t434;
t426 = pkin(5) * t176 + pkin(10) * t143 + t505;
t192 = -qJD(6) * t280 - t434;
t247 = t293 * t280;
t175 = t192 + t247;
t119 = -t171 * t397 + t175 * t394;
t203 = t278 + t279;
t424 = pkin(5) * t203 + pkin(10) * t119 + t38;
t422 = t395 * t286;
t421 = t395 * t423;
t420 = t398 * t286;
t419 = t398 * t423;
t215 = -t299 * t357 + t442;
t347 = t390 * t414 + t508;
t373 = t393 * t438;
t370 = t383 * t449;
t368 = (t386 - t387) * t469;
t367 = -t446 - t369;
t353 = -t369 - t447;
t345 = -t356 + t369;
t344 = t522 - t369;
t343 = t366 - t370;
t342 = t366 + t370;
t341 = -t371 + t365;
t336 = -t356 - t369;
t328 = t356 - t522;
t326 = t392 * t440;
t320 = -t522 - t369;
t312 = -t338 + t522;
t311 = t337 - t522;
t308 = -t522 - t356;
t300 = -t338 + t337;
t292 = -t338 - t522;
t290 = -t522 - t337;
t289 = -t336 * t389 - t488;
t288 = t336 * t392 - t489;
t285 = -t296 + t425;
t284 = t295 - t425;
t276 = t320 * t392 - t541;
t275 = t320 * t389 + t539;
t274 = -t337 - t338;
t265 = t441 + t487;
t264 = t317 * t391 - t340 * t486;
t263 = t317 * t388 + t340 * t485;
t262 = -t339 * t485 + t388 * t441;
t261 = -t339 * t486 - t391 * t441;
t260 = (t339 * t391 + t340 * t388) * t357;
t259 = (t339 * t388 - t340 * t391) * t357;
t258 = -t296 - t425;
t257 = t307 * t389 + t392 * t431;
t254 = t296 - t295;
t252 = t311 * t391 - t495;
t251 = -t312 * t388 + t540;
t250 = t311 * t388 + t494;
t249 = t312 * t391 + t542;
t248 = -t425 - t295;
t246 = -t279 + t291;
t245 = t278 - t291;
t239 = -t292 * t388 - t494;
t238 = t292 * t391 - t495;
t236 = -t420 + t421;
t235 = -t422 - t419;
t234 = t290 * t391 - t542;
t233 = t290 * t388 + t540;
t231 = t279 - t278;
t229 = -t266 * t391 - t270 * t388;
t228 = -t265 * t391 - t388 * t529;
t227 = -t266 * t388 + t270 * t391;
t226 = -t265 * t388 + t391 * t529;
t225 = -t295 - t296;
t218 = t244 + t286;
t214 = (0.2e1 * qJD(5) + t357) * t299 + t442;
t213 = t398 * t244 - t421;
t212 = t395 * t244 + t419;
t211 = -t395 * t243 + t420;
t210 = t398 * t243 + t422;
t209 = t284 * t398 - t497;
t208 = -t285 * t395 + t535;
t207 = t284 * t395 + t496;
t206 = t285 * t398 + t537;
t205 = t239 * t392 + t389 * t529;
t204 = t239 * t389 - t392 * t529;
t200 = -t258 * t395 - t496;
t199 = t258 * t398 - t497;
t198 = pkin(2) * t288 - t224;
t197 = t234 * t392 + t265 * t389;
t196 = t234 * t389 - t265 * t392;
t195 = pkin(2) * t275 - t223;
t194 = t229 * t392 + t274 * t389;
t193 = t229 * t389 - t274 * t392;
t191 = -qJD(6) * t282 - t444;
t190 = (-t280 * t397 + t282 * t394) * t293;
t189 = (-t280 * t394 - t282 * t397) * t293;
t188 = t248 * t398 - t537;
t187 = t248 * t395 + t535;
t183 = -t235 * t388 + t236 * t391;
t182 = t235 * t391 + t236 * t388;
t178 = -qJ(4) * t238 + t498;
t177 = -qJ(4) * t233 + t499;
t174 = t192 - t247;
t170 = t192 * t397 - t282 * t493;
t169 = t192 * t394 + t282 * t492;
t168 = -t191 * t394 + t280 * t492;
t167 = -t191 * t397 - t280 * t493;
t166 = t190 * t398 + t237 * t395;
t165 = t190 * t395 - t237 * t398;
t161 = -t215 * t398 + t218 * t395;
t160 = -t214 * t398 - t395 * t534;
t159 = -t215 * t395 - t218 * t398;
t158 = -t214 * t395 + t398 * t534;
t157 = t245 * t397 - t501;
t156 = -t246 * t394 + t536;
t155 = t245 * t394 + t500;
t154 = t246 * t397 + t538;
t153 = -t212 * t388 + t213 * t391;
t152 = -t210 * t388 + t211 * t391;
t151 = t212 * t391 + t213 * t388;
t150 = t210 * t391 + t211 * t388;
t149 = -t207 * t388 + t209 * t391;
t148 = -t206 * t388 + t208 * t391;
t147 = t207 * t391 + t209 * t388;
t146 = t206 * t391 + t208 * t388;
t145 = -t199 * t388 + t200 * t391;
t144 = t199 * t391 + t200 * t388;
t142 = t221 * t397 - t501;
t140 = t220 * t394 + t536;
t137 = t170 * t398 + t460;
t136 = t168 * t398 - t460;
t135 = t170 * t395 - t459;
t134 = t168 * t395 + t459;
t133 = -t187 * t388 + t188 * t391;
t132 = t187 * t391 + t188 * t388;
t131 = -pkin(3) * t238 + t164;
t130 = t388 * t405 - t391 * (-qJ(3) * t446 - t347 - t523 + t524) + t461 + t406 + (-t302 * t391 - t369 * t388 - t233) * pkin(3);
t127 = -pkin(9) * t199 + t502;
t126 = -pkin(9) * t187 + t503;
t125 = pkin(2) * t204 - pkin(3) * t529 + qJ(4) * t239 + t499;
t124 = t145 * t392 + t389 * t534;
t123 = t145 * t389 - t392 * t534;
t122 = pkin(2) * t196 - pkin(3) * t265 + qJ(4) * t234 - t498;
t121 = t133 * t392 + t214 * t389;
t120 = t133 * t389 - t214 * t392;
t118 = t172 * t397 - t174 * t394;
t117 = -t171 * t394 - t175 * t397;
t116 = t172 * t394 + t174 * t397;
t113 = -pkin(4) * t534 + pkin(9) * t200 + t503;
t112 = t157 * t398 - t171 * t395;
t111 = t156 * t398 + t175 * t395;
t110 = t157 * t395 + t171 * t398;
t109 = t156 * t395 - t175 * t398;
t108 = -t165 * t388 + t166 * t391;
t107 = t165 * t391 + t166 * t388;
t104 = -pkin(4) * t214 + pkin(9) * t188 - t502;
t103 = -t159 * t388 + t161 * t391;
t102 = -t158 * t388 + t160 * t391;
t101 = t159 * t391 + t161 * t388;
t100 = t158 * t391 + t160 * t388;
t99 = t143 * t398 - t176 * t395;
t98 = t143 * t395 + t176 * t398;
t97 = t141 * t398 - t172 * t395;
t96 = t141 * t395 + t172 * t398;
t95 = t118 * t398 + t231 * t395;
t94 = t118 * t395 - t231 * t398;
t93 = -qJ(4) * t227 - t105;
t92 = t119 * t398 - t203 * t395;
t91 = t119 * t395 + t203 * t398;
t90 = t103 * t392 + t225 * t389;
t89 = t103 * t389 - t225 * t392;
t88 = t106 * t392 + t202 * t389;
t87 = t106 * t389 - t202 * t392;
t86 = -t135 * t388 + t137 * t391;
t85 = -t134 * t388 + t136 * t391;
t84 = t135 * t391 + t137 * t388;
t83 = t134 * t391 + t136 * t388;
t79 = -pkin(3) * t101 - pkin(4) * t159;
t78 = pkin(2) * t193 - pkin(3) * t274 + qJ(4) * t229 + t106;
t75 = -t110 * t388 + t112 * t391;
t74 = -t109 * t388 + t111 * t391;
t73 = t110 * t391 + t112 * t388;
t72 = t109 * t391 + t111 * t388;
t71 = -t388 * t98 + t391 * t99;
t70 = t388 * t99 + t391 * t98;
t69 = -t388 * t96 + t391 * t97;
t68 = t388 * t97 + t391 * t96;
t67 = -pkin(10) * t142 + t504;
t66 = -pkin(3) * t144 - pkin(4) * t199 + t82;
t65 = -pkin(10) * t140 + t505;
t64 = -t388 * t94 + t391 * t95;
t63 = t388 * t95 + t391 * t94;
t62 = -qJ(4) * t144 - t113 * t388 + t127 * t391;
t61 = -pkin(3) * t132 - pkin(4) * t187 + t81;
t60 = -t388 * t91 + t391 * t92;
t59 = t388 * t92 + t391 * t91;
t58 = -qJ(4) * t132 - t104 * t388 + t126 * t391;
t54 = pkin(2) * t87 - pkin(3) * t202 + qJ(4) * t106;
t53 = t142 * t389 + t392 * t71;
t52 = -t142 * t392 + t389 * t71;
t49 = t140 * t389 + t392 * t69;
t48 = -t140 * t392 + t389 * t69;
t47 = -pkin(4) * t179 + pkin(9) * t51;
t46 = -pkin(5) * t142 + t57;
t45 = -pkin(5) * t140 + t56;
t44 = -pkin(9) * t159 - t50;
t43 = pkin(2) * t123 - pkin(3) * t534 + qJ(4) * t145 + t113 * t391 + t127 * t388;
t42 = t117 * t389 + t392 * t60;
t41 = -t117 * t392 + t389 * t60;
t40 = -pkin(4) * t225 + pkin(9) * t161 + t51;
t39 = pkin(2) * t120 - pkin(3) * t214 + qJ(4) * t133 + t104 * t391 + t126 * t388;
t36 = t391 * t51 - t507;
t35 = t388 * t51 + t506;
t34 = -pkin(3) * t70 - pkin(4) * t98 - t426;
t33 = t179 * t389 + t36 * t392;
t32 = -t179 * t392 + t36 * t389;
t31 = -pkin(3) * t68 - pkin(4) * t96 - t427;
t30 = -pkin(10) * t117 - t37;
t29 = -pkin(9) * t98 - t395 * t46 + t398 * t67;
t28 = -pkin(9) * t96 - t395 * t45 + t398 * t65;
t27 = t38 * t398 + t395 * t76;
t26 = t38 * t395 - t398 * t76;
t25 = -pkin(4) * t142 + pkin(9) * t99 + t395 * t67 + t398 * t46;
t24 = -pkin(4) * t140 + pkin(9) * t97 + t395 * t65 + t398 * t45;
t23 = -qJ(4) * t101 - t388 * t40 + t391 * t44;
t22 = -pkin(3) * t35 - pkin(4) * t50;
t21 = -pkin(9) * t91 + t117 * t513 + t30 * t398;
t20 = pkin(2) * t89 - pkin(3) * t225 + qJ(4) * t103 + t388 * t44 + t391 * t40;
t19 = pkin(9) * t92 + t117 * t450 + t30 * t395;
t18 = -pkin(3) * t59 - pkin(4) * t91 - t424;
t17 = -pkin(9) * t506 - qJ(4) * t35 - t388 * t47;
t16 = -t26 * t388 + t27 * t391;
t15 = t26 * t391 + t27 * t388;
t14 = -qJ(4) * t70 - t25 * t388 + t29 * t391;
t13 = -qJ(4) * t68 - t24 * t388 + t28 * t391;
t12 = -pkin(9) * t26 + (-pkin(10) * t398 + t513) * t37;
t11 = pkin(2) * t32 - pkin(3) * t179 - pkin(9) * t507 + qJ(4) * t36 + t391 * t47;
t10 = t16 * t392 + t37 * t389;
t9 = t16 * t389 - t37 * t392;
t8 = pkin(2) * t52 - pkin(3) * t142 + qJ(4) * t71 + t25 * t391 + t29 * t388;
t7 = pkin(2) * t48 - pkin(3) * t140 + qJ(4) * t69 + t24 * t391 + t28 * t388;
t6 = -qJ(4) * t59 - t19 * t388 + t21 * t391;
t5 = pkin(9) * t27 + (-pkin(10) * t395 + t450) * t37;
t4 = pkin(2) * t41 - pkin(3) * t117 + qJ(4) * t60 + t19 * t391 + t21 * t388;
t3 = -pkin(3) * t15 - pkin(4) * t26 - t437;
t2 = -qJ(4) * t15 + t12 * t391 - t388 * t5;
t1 = pkin(2) * t9 - pkin(3) * t37 + qJ(4) * t16 + t12 * t388 + t391 * t5;
t55 = [0, 0, 0, 0, 0, qJDD(1), t429, t430, 0, 0, (t365 * t390 - t470 * t533) * t396, t393 * t368 + (t396 * t343 + t399 * t452) * t390, t393 * t341 + (t475 + t399 * (-t447 + t369)) * t390, (t366 * t390 + t471 * t533) * t399, t393 * t342 + (t396 * (t446 - t369) + t473) * t390, t373, (-t318 + pkin(1) * (t363 * t399 + t367 * t396)) * t393 + (t399 * t347 + pkin(1) * t343 + pkin(8) * (t367 * t399 - t475)) * t390, -t347 * t479 - t393 * t319 + pkin(1) * (-t390 * t452 + (t353 * t399 - t364 * t396) * t393) + (-t396 * t353 - t473) * t512, pkin(1) * ((-t341 * t399 + t342 * t396) * t393 - (-t386 - t387) * t385 * t468) + (t396 * t341 + t342 * t399) * t512 + t525 * t390, pkin(1) * (t390 * t347 + (-t318 * t399 + t319 * t396) * t393) + t525 * t512, t453 + (t396 * (t332 * t392 - t359 * t483) + t399 * (t332 * t389 + t359 * t481)) * t390, t393 * t328 + (t396 * (-t302 * t392 - t389 * t527) + t399 * (-t302 * t389 + t392 * t527)) * t390, t393 * t307 + (t396 * (-t345 * t389 + t539) + t399 * (t345 * t392 + t541)) * t390, -t453 + (t396 * (t357 * t481 + t435) + t399 * (t357 * t483 - t326)) * t390, t393 * t431 + (t396 * (t344 * t392 - t489) + t399 * (t344 * t389 + t488)) * t390, t373 + (t396 * (-t357 * t392 + t359 * t389) + t399 * (-t357 * t389 - t359 * t392)) * t482, (t195 + pkin(1) * (t275 * t399 + t276 * t396)) * t393 + (t396 * (-qJ(3) * t275 - t491) + t399 * (-pkin(2) * t302 + qJ(3) * t276 + t490) - pkin(1) * t302 + pkin(8) * (-t396 * t275 + t276 * t399)) * t390, (t198 + pkin(1) * (t288 * t399 + t289 * t396)) * t393 + (t396 * (-qJ(3) * t288 - t490) + t399 * (-pkin(2) * t527 + qJ(3) * t289 - t491) - pkin(1) * t527 + pkin(8) * (-t396 * t288 + t289 * t399)) * t390, (t515 + pkin(1) * (t256 * t399 + t257 * t396)) * t393 + (t396 * (-qJ(3) * t256 - t180) + t399 * (-pkin(2) * t308 + qJ(3) * t257 + t181) - pkin(1) * t308 + pkin(8) * (-t396 * t256 + t257 * t399)) * t390, (t516 + pkin(1) * (t180 * t399 + t181 * t396)) * t393 + (-qJ(3) * t477 + t399 * (pkin(2) * t294 + qJ(3) * t181) + pkin(1) * t294 + pkin(8) * (t181 * t399 - t477)) * t390, t393 * t263 + (t396 * (t264 * t392 - t456) + t399 * (t264 * t389 + t455)) * t390, t393 * t226 + (t396 * (t228 * t392 - t300 * t389) + t399 * (t228 * t389 + t300 * t392)) * t390, t393 * t249 + (t396 * (t251 * t392 - t270 * t389) + t399 * (t251 * t389 + t270 * t392)) * t390, t393 * t261 + (t396 * (t262 * t392 + t456) + t399 * (t262 * t389 - t455)) * t390, t393 * t250 + (t396 * (t252 * t392 - t266 * t389) + t399 * (t252 * t389 + t266 * t392)) * t390, t393 * t259 + (t396 * (t392 * t260 + t435) + t399 * (t260 * t389 - t326)) * t390, (t122 + pkin(1) * (t196 * t399 + t197 * t396)) * t393 + (t396 * (-qJ(3) * t196 - t130 * t389 + t177 * t392) + t399 * (-pkin(2) * t233 + qJ(3) * t197 + t130 * t392 + t177 * t389) - pkin(1) * t233 + pkin(8) * (-t396 * t196 + t197 * t399)) * t390, (t125 + pkin(1) * (t204 * t399 + t205 * t396)) * t393 + (t396 * (-qJ(3) * t204 - t131 * t389 + t178 * t392) + t399 * (-pkin(2) * t238 + qJ(3) * t205 + t131 * t392 + t178 * t389) - pkin(1) * t238 + pkin(8) * (-t396 * t204 + t205 * t399)) * t390, (t78 + pkin(1) * (t193 * t399 + t194 * t396)) * t393 + (t396 * (-qJ(3) * t193 + t392 * t93) + t399 * (qJ(3) * t194 + t389 * t93) + pkin(8) * (-t396 * t193 + t194 * t399) + (t396 * t514 + t399 * t451 - pkin(1)) * t227) * t390, (t54 + pkin(1) * (t396 * t88 + t399 * t87)) * t393 + ((t396 * (-qJ(4) * t392 + t514) + t399 * (-qJ(4) * t389 + t451) - pkin(1)) * t105 + (pkin(8) + qJ(3)) * (-t396 * t87 + t399 * t88)) * t390, t393 * t151 + (t396 * (t153 * t392 + t458) + t399 * (t153 * t389 - t457)) * t390, t393 * t100 + (t396 * (t102 * t392 + t254 * t389) + t399 * (t102 * t389 - t254 * t392)) * t390, t393 * t146 + (t396 * (t148 * t392 + t218 * t389) + t399 * (t148 * t389 - t218 * t392)) * t390, t393 * t150 + (t396 * (t152 * t392 - t458) + t399 * (t152 * t389 + t457)) * t390, t393 * t147 + (t396 * (t149 * t392 - t215 * t389) + t399 * (t149 * t389 + t215 * t392)) * t390, t393 * t182 + (t396 * (t183 * t392 + t329 * t389) + t399 * (t183 * t389 - t329 * t392)) * t390, (t39 + pkin(1) * (t120 * t399 + t121 * t396)) * t393 + (t396 * (-qJ(3) * t120 - t389 * t61 + t392 * t58) + t399 * (-pkin(2) * t132 + qJ(3) * t121 + t389 * t58 + t392 * t61) - pkin(1) * t132 + pkin(8) * (-t396 * t120 + t121 * t399)) * t390, (t43 + pkin(1) * (t123 * t399 + t124 * t396)) * t393 + (t396 * (-qJ(3) * t123 - t389 * t66 + t392 * t62) + t399 * (-pkin(2) * t144 + qJ(3) * t124 + t389 * t62 + t392 * t66) - pkin(1) * t144 + pkin(8) * (-t396 * t123 + t124 * t399)) * t390, (t20 + pkin(1) * (t396 * t90 + t399 * t89)) * t393 + (t396 * (-qJ(3) * t89 + t23 * t392 - t389 * t79) + t399 * (-pkin(2) * t101 + qJ(3) * t90 + t23 * t389 + t392 * t79) - pkin(1) * t101 + pkin(8) * (-t396 * t89 + t399 * t90)) * t390, (t11 + pkin(1) * (t32 * t399 + t33 * t396)) * t393 + (t396 * (-qJ(3) * t32 + t17 * t392 - t22 * t389) + t399 * (-pkin(2) * t35 + qJ(3) * t33 + t17 * t389 + t22 * t392) - pkin(1) * t35 + pkin(8) * (-t396 * t32 + t33 * t399)) * t390, t393 * t84 + (t396 * (t169 * t389 + t392 * t86) + t399 * (-t169 * t392 + t389 * t86)) * t390, t393 * t63 + (t396 * (t116 * t389 + t392 * t64) + t399 * (-t116 * t392 + t389 * t64)) * t390, t393 * t72 + (t396 * (t154 * t389 + t392 * t74) + t399 * (-t154 * t392 + t389 * t74)) * t390, t393 * t83 + (t396 * (-t167 * t389 + t392 * t85) + t399 * (t167 * t392 + t389 * t85)) * t390, t393 * t73 + (t396 * (t155 * t389 + t392 * t75) + t399 * (-t155 * t392 + t389 * t75)) * t390, t393 * t107 + (t396 * (t108 * t392 + t189 * t389) + t399 * (t108 * t389 - t189 * t392)) * t390, (t7 + pkin(1) * (t396 * t49 + t399 * t48)) * t393 + (t396 * (-qJ(3) * t48 + t13 * t392 - t31 * t389) + t399 * (-pkin(2) * t68 + qJ(3) * t49 + t13 * t389 + t31 * t392) - pkin(1) * t68 + pkin(8) * (-t396 * t48 + t399 * t49)) * t390, (t8 + pkin(1) * (t396 * t53 + t399 * t52)) * t393 + (t396 * (-qJ(3) * t52 + t14 * t392 - t34 * t389) + t399 * (-pkin(2) * t70 + qJ(3) * t53 + t14 * t389 + t34 * t392) - pkin(1) * t70 + pkin(8) * (-t396 * t52 + t399 * t53)) * t390, (t4 + pkin(1) * (t396 * t42 + t399 * t41)) * t393 + (t396 * (-qJ(3) * t41 - t18 * t389 + t392 * t6) + t399 * (-pkin(2) * t59 + qJ(3) * t42 + t18 * t392 + t389 * t6) - pkin(1) * t59 + pkin(8) * (-t396 * t41 + t399 * t42)) * t390, (t1 + pkin(1) * (t10 * t396 + t399 * t9)) * t393 + (t396 * (-qJ(3) * t9 + t2 * t392 - t3 * t389) + t399 * (-pkin(2) * t15 + qJ(3) * t10 + t2 * t389 + t3 * t392) - pkin(1) * t15 + pkin(8) * (t10 * t399 - t396 * t9)) * t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377, t368, t341, t377, t342, t438, -t318, -t319, 0, 0, t331, t328, t307, -t331, t431, t438, t195, t198, t515, t516, t263, t226, t249, t261, t250, t259, t122, t125, t78, t54, t151, t100, t146, t150, t147, t182, t39, t43, t20, t11, t84, t63, t72, t83, t73, t107, t7, t8, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t527, t308, -t294, 0, 0, 0, 0, 0, 0, t233, t238, t227, t105, 0, 0, 0, 0, 0, 0, t132, t144, t101, t35, 0, 0, 0, 0, 0, 0, t68, t70, t59, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, t529, t274, t202, 0, 0, 0, 0, 0, 0, t214, t534, t225, t179, 0, 0, 0, 0, 0, 0, t140, t142, t117, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, t254, t218, -t255, -t215, t329, -t81, -t82, 0, 0, t169, t116, t154, -t167, t155, t189, t427, t426, t424, t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t231, t175, -t232, -t171, t237, -t56, -t57, 0, 0;];
tauJ_reg  = t55;
