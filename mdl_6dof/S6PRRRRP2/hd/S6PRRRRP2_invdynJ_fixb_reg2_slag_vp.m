% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:27
% EndTime: 2019-03-09 00:04:46
% DurationCPUTime: 10.12s
% Computational Cost: add. (11411->693), mult. (26015->867), div. (0->0), fcn. (20307->14), ass. (0->353)
t266 = sin(qJ(2));
t261 = sin(pkin(6));
t409 = qJD(1) * t261;
t378 = t266 * t409;
t265 = sin(qJ(3));
t404 = qJD(3) * t265;
t307 = pkin(3) * t404 - t378;
t267 = cos(qJ(5));
t400 = qJD(5) * t267;
t264 = sin(qJ(4));
t406 = qJD(2) * t265;
t373 = t264 * t406;
t268 = cos(qJ(3));
t479 = cos(qJ(4));
t379 = t479 * t268;
t206 = -qJD(2) * t379 + t373;
t436 = t206 * t267;
t511 = t400 + t436;
t369 = t479 * qJD(3);
t347 = t268 * t369;
t393 = qJD(3) + qJD(4);
t352 = t264 * t393;
t368 = t479 * qJD(4);
t153 = t265 * t352 - t268 * t368 - t347;
t418 = t264 * t268;
t212 = t479 * t265 + t418;
t154 = t393 * t212;
t510 = pkin(4) * t154 + pkin(10) * t153 + t307;
t478 = pkin(3) * t264;
t250 = pkin(10) + t478;
t350 = pkin(3) * t368;
t330 = t267 * t350;
t263 = sin(qJ(5));
t401 = qJD(5) * t263;
t488 = -t250 * t401 + t330;
t208 = t212 * qJD(2);
t196 = t208 * qJ(6);
t137 = pkin(4) * t208 + pkin(10) * t206;
t390 = pkin(3) * t406;
t121 = t137 + t390;
t218 = qJD(2) * pkin(8) + t378;
t357 = pkin(9) * qJD(2) + t218;
t262 = cos(pkin(6));
t408 = qJD(1) * t262;
t376 = t265 * t408;
t159 = t268 * t357 + t376;
t150 = t264 * t159;
t239 = t268 * t408;
t158 = -t265 * t357 + t239;
t93 = t479 * t158 - t150;
t58 = t263 * t121 + t267 * t93;
t51 = t196 + t58;
t509 = -t51 + t488;
t508 = t58 - t488;
t270 = -pkin(9) - pkin(8);
t380 = qJD(3) * t270;
t215 = t265 * t380;
t319 = -t264 * t265 + t379;
t269 = cos(qJ(2));
t377 = t269 * t409;
t228 = t270 * t265;
t229 = t270 * t268;
t489 = t479 * t228 + t264 * t229;
t496 = qJD(4) * t489 + t479 * t215 - t319 * t377 + t380 * t418;
t339 = t267 * pkin(5) + t263 * qJ(6);
t403 = qJD(4) * t264;
t151 = t479 * t159;
t92 = t158 * t264 + t151;
t345 = pkin(3) * t403 - t92;
t437 = t206 * t263;
t507 = -qJD(6) * t263 - t511 * qJ(6) + (t401 + t437) * pkin(5);
t424 = t261 * t266;
t204 = t262 * t268 - t265 * t424;
t454 = cos(pkin(11));
t359 = t454 * t266;
t260 = sin(pkin(11));
t425 = t260 * t269;
t201 = t262 * t359 + t425;
t259 = qJ(3) + qJ(4);
t254 = sin(t259);
t255 = cos(t259);
t360 = t261 * t454;
t141 = -t201 * t254 - t255 * t360;
t358 = t454 * t269;
t426 = t260 * t266;
t203 = -t262 * t426 + t358;
t427 = t260 * t261;
t143 = -t203 * t254 + t255 * t427;
t186 = -t254 * t424 + t255 * t262;
t506 = g(1) * t143 + g(2) * t141 + g(3) * t186;
t351 = t267 * t393;
t162 = t208 * t263 - t351;
t164 = t267 * t208 + t263 * t393;
t441 = t164 * t267;
t444 = t162 * t263;
t303 = t369 + t368;
t361 = qJDD(2) * t479;
t398 = qJD(2) * qJD(3);
t365 = t265 * t398;
t412 = -qJD(4) * t373 - t264 * t365;
t275 = t265 * t361 + (qJD(2) * t303 + qJDD(2) * t264) * t268 + t412;
t392 = qJDD(3) + qJDD(4);
t68 = qJD(5) * t164 + t263 * t275 - t267 * t392;
t463 = t267 * t68;
t67 = -qJD(5) * t351 + t208 * t401 - t263 * t392 - t267 * t275;
t464 = t263 * t67;
t505 = t212 * (qJD(5) * (-t441 + t444) - t463 + t464) + (t162 * t267 + t164 * t263) * t153;
t199 = qJD(5) + t206;
t442 = t164 * t199;
t446 = t162 * t199;
t9 = t263 * (t68 + t442) + t267 * (t67 + t446);
t395 = t265 * qJDD(2);
t336 = t264 * t395 - t268 * t361;
t112 = qJD(2) * t154 + t336;
t108 = qJDD(5) + t112;
t449 = t108 * t267;
t504 = pkin(10) * (t199 * t401 - t449);
t142 = t201 * t255 - t254 * t360;
t144 = t203 * t255 + t254 * t427;
t187 = t254 * t262 + t255 * t424;
t309 = g(1) * t144 + g(2) * t142 + g(3) * t187;
t252 = t268 * pkin(3) + pkin(2);
t197 = -qJD(2) * t252 - t377;
t111 = pkin(4) * t206 - pkin(10) * t208 + t197;
t152 = qJD(3) * pkin(3) + t158;
t399 = qJD(1) * qJD(2);
t366 = t269 * t399;
t451 = qJDD(2) * pkin(8);
t185 = t451 + (qJDD(1) * t266 + t366) * t261;
t397 = qJDD(1) * t262;
t238 = t268 * t397;
t78 = qJDD(3) * pkin(3) + t238 + (-pkin(9) * qJDD(2) - t185) * t265 - t159 * qJD(3);
t394 = t268 * qJDD(2);
t381 = qJD(3) * t239 + t268 * t185 + t265 * t397;
t98 = -t218 * t404 + t381;
t83 = (-t365 + t394) * pkin(9) + t98;
t19 = t152 * t368 - t159 * t403 + t264 * t78 + t479 * t83;
t17 = pkin(10) * t392 + t19;
t367 = t266 * t399;
t422 = t261 * t269;
t337 = -qJDD(1) * t422 + t261 * t367;
t452 = qJDD(2) * pkin(2);
t184 = t337 - t452;
t247 = pkin(3) * t365;
t44 = -pkin(3) * t394 + t112 * pkin(4) - pkin(10) * t275 + t184 + t247;
t89 = t264 * t152 + t151;
t81 = pkin(10) * t393 + t89;
t6 = t111 * t400 + t267 * t17 + t263 * t44 - t401 * t81;
t503 = t6 * t267 - t309;
t353 = t152 * t403 + t159 * t368 + t264 * t83 - t479 * t78;
t18 = -pkin(4) * t392 + t353;
t88 = t479 * t152 - t150;
t80 = -pkin(4) * t393 - t88;
t502 = t18 * t263 + t80 * t400;
t430 = t250 * t267;
t501 = -t162 * t330 - t68 * t430;
t460 = t507 + t345;
t139 = -pkin(4) * t319 - pkin(10) * t212 - t252;
t167 = t264 * t228 - t479 * t229;
t458 = t139 * t400 - t167 * t401 + t263 * t510 + t496 * t267;
t500 = t350 - t93;
t498 = t68 - t442;
t455 = qJD(4) * t167 - t212 * t377 + t264 * t215 - t270 * t347;
t495 = t208 * t393;
t493 = t339 * t141;
t492 = t339 * t143;
t491 = t263 * t139 + t267 * t167;
t490 = t339 * t186;
t372 = t250 * t400;
t487 = -t263 * t350 - t372;
t453 = qJ(6) * t108;
t2 = qJD(6) * t199 + t453 + t6;
t356 = t111 * t401 + t263 * t17 - t267 * t44 + t81 * t400;
t476 = pkin(5) * t108;
t4 = qJDD(6) + t356 - t476;
t486 = t2 * t267 + t4 * t263 - t309;
t448 = t153 * t263;
t317 = t212 * t400 - t448;
t434 = t212 * t263;
t485 = t108 * t434 + t154 * t162 + t199 * t317 - t319 * t68;
t271 = qJD(3) ^ 2;
t200 = -t262 * t358 + t426;
t202 = t262 * t425 + t359;
t343 = g(1) * t202 + g(2) * t200;
t322 = -t184 + t343;
t484 = -pkin(8) * t271 + t261 * (-g(3) * t269 + t367) + t322 + t452;
t457 = -qJD(5) * t491 - t496 * t263 + t267 * t510;
t481 = t164 ^ 2;
t480 = t199 ^ 2;
t477 = pkin(4) * t255;
t475 = pkin(5) * t208;
t474 = pkin(10) * t108;
t8 = t68 * pkin(5) + t67 * qJ(6) - t164 * qJD(6) + t18;
t469 = t8 * t263;
t468 = qJD(2) * pkin(2);
t49 = t111 * t263 + t267 * t81;
t467 = t199 * t49;
t53 = t162 * pkin(5) - t164 * qJ(6) + t80;
t466 = t206 * t53;
t465 = t206 * t80;
t462 = qJ(6) * t154 - qJD(6) * t319 + t458;
t461 = -pkin(5) * t154 - t457;
t459 = -t89 + t507;
t338 = pkin(5) * t263 - qJ(6) * t267;
t456 = -t338 * t153 + (qJD(5) * t339 - qJD(6) * t267) * t212 + t455;
t60 = t263 * t137 + t267 * t88;
t450 = t108 * t250;
t447 = t153 * t267;
t445 = t162 * t250;
t443 = t164 * t162;
t440 = t199 * t208;
t439 = t200 * t254;
t438 = t202 * t254;
t435 = t208 * t206;
t433 = t212 * t267;
t432 = t218 * t265;
t431 = t250 * t263;
t429 = t255 * t263;
t428 = t255 * t267;
t423 = t261 * t268;
t417 = t267 * t269;
t48 = t111 * t267 - t263 * t81;
t416 = qJD(6) - t48;
t415 = qJDD(1) - g(3);
t414 = -t200 * t252 - t201 * t270;
t413 = -t202 * t252 - t203 * t270;
t257 = t265 ^ 2;
t258 = t268 ^ 2;
t411 = t257 - t258;
t410 = t257 + t258;
t219 = -t377 - t468;
t407 = qJD(2) * t219;
t405 = qJD(2) * t266;
t402 = qJD(5) * t199;
t391 = t479 * pkin(3);
t387 = t254 * t422;
t385 = t261 * t417;
t233 = t263 * t422;
t272 = qJD(2) ^ 2;
t384 = t265 * t272 * t268;
t47 = t53 * t401;
t383 = t53 * t400;
t71 = t80 * t401;
t382 = t506 * t263;
t375 = t261 * t405;
t374 = qJD(2) * t422;
t370 = t162 ^ 2 - t481;
t364 = t269 * t398;
t36 = -pkin(5) * t199 + t416;
t363 = t36 * t208 + t47;
t362 = -t48 * t208 + t71;
t355 = t199 * t263;
t354 = t199 * t267;
t349 = -pkin(10) * t439 - t200 * t477 + t414;
t348 = -pkin(10) * t438 - t202 * t477 + t413;
t346 = t268 * t365;
t344 = (-t203 * t265 + t260 * t423) * pkin(3);
t342 = g(1) * t203 + g(2) * t201;
t340 = t252 * t422 - t270 * t424;
t37 = qJ(6) * t199 + t49;
t335 = t263 * t37 - t267 * t36;
t334 = t263 * t49 + t267 * t48;
t57 = t121 * t267 - t263 * t93;
t59 = t137 * t267 - t263 * t88;
t329 = -t162 * t208 - t449;
t90 = t139 * t267 - t167 * t263;
t325 = t441 + t444;
t324 = t204 * pkin(3);
t323 = qJDD(2) * t269 - t266 * t272;
t221 = -pkin(4) - t339;
t205 = t262 * t265 + t266 * t423;
t123 = t264 * t204 + t479 * t205;
t109 = t123 * t263 + t385;
t320 = t479 * t204 - t264 * t205;
t175 = t218 * t268 + t376;
t316 = -t212 * t401 - t447;
t110 = t123 * t267 - t233;
t155 = qJD(3) * t204 + t268 * t374;
t156 = -qJD(3) * t205 - t265 * t374;
t62 = t320 * qJD(4) + t479 * t155 + t264 * t156;
t26 = -qJD(5) * t109 + t263 * t375 + t267 * t62;
t27 = -qJD(5) * t233 + t123 * t400 + t263 * t62 - t267 * t375;
t315 = -t109 * t67 - t110 * t68 - t26 * t162 + t164 * t27;
t314 = pkin(10) * t387 + t422 * t477 + t340;
t63 = t123 * qJD(4) + t264 * t155 - t479 * t156;
t313 = -t108 * t109 + t63 * t162 - t199 * t27 - t320 * t68;
t114 = -t200 * t429 - t201 * t267;
t116 = -t202 * t429 - t203 * t267;
t172 = t233 * t255 - t267 * t424;
t312 = -g(1) * t116 - g(2) * t114 - g(3) * t172;
t115 = -t200 * t428 + t201 * t263;
t117 = -t202 * t428 + t203 * t263;
t173 = (t255 * t417 + t263 * t266) * t261;
t311 = -g(1) * t117 - g(2) * t115 - g(3) * t173;
t306 = t49 * t208 + t382 + t502;
t136 = t143 * pkin(4);
t305 = pkin(10) * t144 + t136 + t344;
t304 = t162 * t355 - t463;
t302 = -g(3) * t422 + t343;
t301 = t164 * t350 - t250 * t67;
t300 = -t199 * t350 - t450;
t183 = t186 * pkin(4);
t298 = pkin(10) * t187 + t183 + t324;
t297 = -t37 * t208 - t53 * t436 - t382 - t469;
t296 = (-t201 * t265 - t268 * t360) * pkin(3);
t294 = -t48 * t436 - t49 * t437 + t503;
t293 = t302 * t254;
t292 = t108 * t110 - t164 * t63 + t199 * t26 - t320 * t67;
t291 = -pkin(10) * t402 - t506;
t100 = t142 * t263 - t200 * t267;
t102 = t144 * t263 - t202 * t267;
t147 = t187 * t263 + t385;
t289 = g(1) * t102 + g(2) * t100 + g(3) * t147 - t356;
t288 = -t250 * t402 - t506;
t287 = t162 * t317 + t434 * t68;
t286 = -pkin(8) * qJDD(3) + (t219 + t377 - t468) * qJD(3);
t135 = t141 * pkin(4);
t285 = t142 * pkin(10) + t135 + t296;
t101 = t142 * t267 + t200 * t263;
t103 = t144 * t267 + t202 * t263;
t148 = t187 * t267 - t233;
t284 = -g(1) * t103 - g(2) * t101 - g(3) * t148 + t6;
t283 = -t197 * t208 - t353 - t506;
t282 = t197 * t206 - t19 + t309;
t281 = t511 * t36 - t37 * t437 + t486;
t280 = t164 * t53 + qJDD(6) - t289;
t174 = t239 - t432;
t99 = -t175 * qJD(3) - t265 * t185 + t238;
t277 = -t99 * t265 + t98 * t268 + (-t174 * t268 - t175 * t265) * qJD(3) - t342;
t274 = -t336 - t495;
t251 = -t391 - pkin(4);
t210 = -t391 + t221;
t140 = -qJDD(2) * t252 + t247 + t337;
t113 = -t206 ^ 2 + t208 ^ 2;
t104 = pkin(5) * t164 + qJ(6) * t162;
t95 = t212 * t338 - t489;
t85 = t274 + t495;
t84 = t206 * t393 + t275;
t70 = pkin(5) * t319 - t90;
t69 = -qJ(6) * t319 + t491;
t66 = pkin(10) * t463;
t56 = -t108 * t319 + t154 * t199;
t55 = -t59 - t475;
t54 = t196 + t60;
t52 = -t57 - t475;
t45 = -t67 + t446;
t34 = t108 * t263 - t164 * t208 + t199 * t354;
t33 = -t480 * t263 - t329;
t32 = t199 * t355 + t329;
t28 = t164 * t354 - t464;
t21 = t164 * t316 - t433 * t67;
t11 = t108 * t433 + t154 * t164 + t199 * t316 + t319 * t67;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t415, 0, 0, 0, 0, 0, 0, t323 * t261 (-qJDD(2) * t266 - t269 * t272) * t261, 0, -g(3) + (t262 ^ 2 + (t266 ^ 2 + t269 ^ 2) * t261 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, qJD(3) * t156 + qJDD(3) * t204 + (-t265 * t364 + t268 * t323) * t261, -qJD(3) * t155 - qJDD(3) * t205 + (-t265 * t323 - t268 * t364) * t261 (-t204 * t265 + t205 * t268) * qJDD(2) + (t155 * t268 - t156 * t265 + (-t204 * t268 - t205 * t265) * qJD(3)) * qJD(2), t155 * t175 + t156 * t174 + t204 * t99 + t205 * t98 - g(3) + (-t184 * t269 + t219 * t405) * t261, 0, 0, 0, 0, 0, 0, -t63 * t393 + t320 * t392 + (-t112 * t269 + t206 * t405) * t261, -t62 * t393 - t123 * t392 + (-t269 * (qJDD(2) * t212 + t412) + (-t268 * t269 * t303 + t266 * t208) * qJD(2)) * t261, -t123 * t112 - t62 * t206 + t63 * t208 - t275 * t320, -t320 * t353 + t123 * t19 + t62 * t89 - t63 * t88 - g(3) + (-t140 * t269 + t197 * t405) * t261, 0, 0, 0, 0, 0, 0, t313, -t292, t315, t109 * t356 + t110 * t6 - t18 * t320 + t26 * t49 - t27 * t48 + t63 * t80 - g(3), 0, 0, 0, 0, 0, 0, t313, t315, t292, t109 * t4 + t110 * t2 + t26 * t37 + t27 * t36 - t320 * t8 + t53 * t63 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t415 * t422 + t343, -t415 * t424 + t342, 0, 0, qJDD(2) * t257 + 0.2e1 * t346, 0.2e1 * t265 * t394 - 0.2e1 * t398 * t411, qJDD(3) * t265 + t268 * t271, qJDD(2) * t258 - 0.2e1 * t346, qJDD(3) * t268 - t265 * t271, 0, t286 * t265 + t268 * t484, -t265 * t484 + t286 * t268, t410 * t451 + (-g(3) * t266 - t366 * t410) * t261 + t277, t322 * pkin(2) + t277 * pkin(8) + (-g(3) * (pkin(2) * t269 + pkin(8) * t266) + (-t219 * t266 + (t174 * t265 - t175 * t268) * t269) * qJD(1)) * t261, -t208 * t153 + t212 * t275, -t212 * t112 + t153 * t206 - t208 * t154 + t275 * t319, -t153 * t393 + t212 * t392, -t112 * t319 + t154 * t206, -t154 * t393 + t319 * t392, 0, -t252 * t112 - t140 * t319 + t197 * t154 + t206 * t307 + t255 * t302 + t392 * t489 - t393 * t455, -g(1) * t438 - g(2) * t439 + g(3) * t387 + t140 * t212 - t197 * t153 - t167 * t392 + t208 * t307 - t252 * t275 - t393 * t496, -g(3) * t424 - t167 * t112 + t88 * t153 - t89 * t154 + t19 * t319 - t206 * t496 + t208 * t455 + t212 * t353 - t275 * t489 - t342, -g(1) * t413 - g(2) * t414 - g(3) * t340 - t140 * t252 + t19 * t167 + t307 * t197 - t353 * t489 - t455 * t88 + t496 * t89, t21, t505, t11, t287, -t485, t56, t108 * t90 + t154 * t48 + t455 * t162 + t457 * t199 + t502 * t212 + t319 * t356 - t80 * t448 - t489 * t68 + t311, -t80 * t447 - t108 * t491 - t154 * t49 + t489 * t67 + t319 * t6 + (t18 * t267 - t71) * t212 - t458 * t199 + t455 * t164 - t312, t67 * t90 - t68 * t491 - t457 * t164 - t458 * t162 + t334 * t153 + t293 + (-t263 * t6 + t267 * t356 + (t263 * t48 - t267 * t49) * qJD(5)) * t212, -g(1) * t348 - g(2) * t349 - g(3) * t314 - t18 * t489 - t356 * t90 + t455 * t80 + t457 * t48 + t458 * t49 + t491 * t6, t21, t11, -t505, t56, t485, t287, -t53 * t448 - t108 * t70 - t154 * t36 + t319 * t4 + t68 * t95 + (t383 + t469) * t212 - t461 * t199 + t456 * t162 + t311, -t67 * t70 - t68 * t69 + t461 * t164 - t462 * t162 + t335 * t153 + t293 + (-t2 * t263 + t267 * t4 + (-t263 * t36 - t267 * t37) * qJD(5)) * t212, t53 * t447 + t108 * t69 + t154 * t37 - t2 * t319 + t67 * t95 + (-t8 * t267 + t47) * t212 + t462 * t199 - t456 * t164 + t312, t2 * t69 + t8 * t95 + t4 * t70 - g(1) * (pkin(5) * t117 + qJ(6) * t116 + t348) - g(2) * (pkin(5) * t115 + qJ(6) * t114 + t349) - g(3) * (pkin(5) * t173 + qJ(6) * t172 + t314) + t456 * t53 + t462 * t37 + t461 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, t411 * t272, t395, t384, t394, qJDD(3), -g(3) * t204 + t238 + (-g(1) * t260 + g(2) * t454) * t423 + (-t185 + t342 - t407) * t265, -t268 * t407 - g(1) * (-t203 * t268 - t265 * t427) - g(2) * (-t201 * t268 + t265 * t360) + g(3) * t205 + (t174 + t432) * qJD(3) - t381, 0, 0, t435, t113, t84, -t435, t85, t392, t92 * t393 + (-qJD(4) * t352 - t206 * t406 + t479 * t392) * pkin(3) + t283, t93 * t393 + (-t208 * t406 - t264 * t392 - t368 * t393) * pkin(3) + t282, -t112 * t478 - t275 * t391 + (t89 + t345) * t208 + (-t88 - t500) * t206, -g(1) * t344 - g(2) * t296 - g(3) * t324 + t19 * t478 - t197 * t390 - t345 * t88 - t353 * t391 + t500 * t89, t28, -t9, t34, t304, t33, -t440, -t57 * t199 + t251 * t68 + t345 * t162 + (t300 + t465) * t263 + (-t18 + t288) * t267 + t362, -t251 * t67 + (-t450 + t465) * t267 + t345 * t164 + t508 * t199 + t306, t58 * t162 + t57 * t164 + (t164 * t250 - t48) * t400 + (t356 + (-t49 + t445) * qJD(5) + t301) * t263 + t294 + t501, -g(1) * t305 - g(2) * t285 - g(3) * t298 + t18 * t251 + t430 * t6 + t431 * t356 + t345 * t80 - t508 * t49 + (-t57 + t487) * t48, t28, t34, t9, -t440, t32, t304, t52 * t199 + t210 * t68 + t460 * t162 + (t300 + t466) * t263 + (t288 - t8) * t267 + t363, t51 * t162 + (-t52 + t372) * t164 + ((-t37 + t445) * qJD(5) + t301) * t263 + t281 + t501, t210 * t67 + (-qJD(5) * t53 + t450) * t267 - t460 * t164 + t509 * t199 + t297, t2 * t430 + t8 * t210 + t4 * t431 - g(1) * (t305 + t492) - g(2) * (t285 + t493) - g(3) * (t298 + t490) + t460 * t53 + t509 * t37 + (-t52 - t487) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t435, t113, t84, -t435, t85, t392, t393 * t89 + t283, t393 * t88 + t282, 0, 0, t28, -t9, t34, t304, t33, -t440, -pkin(4) * t68 - t162 * t89 - t199 * t59 + (t465 - t474) * t263 + (-t18 + t291) * t267 + t362, pkin(4) * t67 - t164 * t89 + t199 * t60 + t436 * t80 + t306 + t504, t162 * t60 + t164 * t59 - t66 + (-pkin(10) * t67 + t356) * t263 + (pkin(10) * t325 - t334) * qJD(5) + t294, -t18 * pkin(4) - g(1) * t136 - g(2) * t135 - g(3) * t183 - t48 * t59 - t49 * t60 - t80 * t89 + (-qJD(5) * t334 + t263 * t356 + t503) * pkin(10), t28, t34, t9, -t440, t32, t304, t199 * t55 + t221 * t68 + (t466 - t474) * t263 + t459 * t162 + (t291 - t8) * t267 + t363, -t37 * t401 + t162 * t54 - t164 * t55 - t66 + (qJD(5) * t325 - t464) * pkin(10) + t281, -t164 * t459 - t199 * t54 + t221 * t67 + t297 - t383 - t504, t8 * t221 - t37 * t54 - t36 * t55 - g(1) * (t136 + t492) - g(2) * (t135 + t493) - g(3) * (t183 + t490) + t459 * t53 + (-qJD(5) * t335 + t486) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t443, -t370, t45, -t443, -t498, t108, -t164 * t80 + t289 + t467, t162 * t80 + t199 * t48 - t284, 0, 0, t443, t45, t370, t108, t498, -t443, -t104 * t162 - t280 + t467 + 0.2e1 * t476, pkin(5) * t67 - qJ(6) * t68 + (t37 - t49) * t164 + (t36 - t416) * t162, 0.2e1 * t453 + t104 * t164 - t162 * t53 + (0.2e1 * qJD(6) - t48) * t199 + t284, t2 * qJ(6) - t4 * pkin(5) - t53 * t104 - t36 * t49 - g(1) * (-pkin(5) * t102 + qJ(6) * t103) - g(2) * (-pkin(5) * t100 + qJ(6) * t101) - g(3) * (-pkin(5) * t147 + qJ(6) * t148) + t416 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) + t274 + t443, t45, -t480 - t481, -t199 * t37 + t280 - t476;];
tau_reg  = t1;