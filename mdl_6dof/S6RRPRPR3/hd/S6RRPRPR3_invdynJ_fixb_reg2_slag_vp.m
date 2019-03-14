% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:16
% EndTime: 2019-03-09 10:19:39
% DurationCPUTime: 13.79s
% Computational Cost: add. (22736->739), mult. (53541->928), div. (0->0), fcn. (41175->18), ass. (0->330)
t323 = cos(qJ(2));
t437 = cos(pkin(10));
t375 = t437 * t323;
t284 = qJD(1) * t375;
t314 = sin(pkin(10));
t320 = sin(qJ(2));
t397 = qJD(1) * t320;
t251 = t314 * t397 - t284;
t267 = t314 * t323 + t320 * t437;
t255 = t267 * qJD(1);
t166 = pkin(2) * t397 + pkin(3) * t255 + pkin(8) * t251;
t454 = qJ(3) + pkin(7);
t279 = t454 * t323;
t273 = qJD(1) * t279;
t258 = t314 * t273;
t278 = t454 * t320;
t272 = qJD(1) * t278;
t197 = -t272 * t437 - t258;
t319 = sin(qJ(4));
t322 = cos(qJ(4));
t118 = t319 * t166 + t322 * t197;
t470 = pkin(2) * t314;
t292 = pkin(8) + t470;
t403 = qJ(5) + t292;
t373 = qJD(4) * t403;
t418 = t251 * t319;
t500 = -qJ(5) * t418 + t322 * qJD(5) - t319 * t373 - t118;
t117 = t322 * t166 - t197 * t319;
t499 = -pkin(4) * t255 - t319 * qJD(5) - t117 + (-qJ(5) * t251 - t373) * t322;
t208 = -t322 * qJD(2) + t255 * t319;
t210 = qJD(2) * t319 + t255 * t322;
t313 = sin(pkin(11));
t315 = cos(pkin(11));
t141 = t208 * t315 + t313 * t210;
t318 = sin(qJ(6));
t358 = -t208 * t313 + t315 * t210;
t471 = cos(qJ(6));
t382 = qJD(6) * t471;
t394 = qJD(6) * t318;
t393 = qJD(1) * qJD(2);
t381 = t320 * t393;
t336 = qJDD(1) * t267 - t314 * t381;
t331 = qJD(2) * t284 + t336;
t392 = qJD(2) * qJD(4);
t396 = qJD(4) * t319;
t130 = -t319 * qJDD(2) + t255 * t396 + (-t331 - t392) * t322;
t328 = -t322 * qJDD(2) + t319 * t331;
t131 = qJD(4) * t210 + t328;
t68 = -t130 * t313 + t315 * t131;
t69 = -t130 * t315 - t131 * t313;
t23 = t141 * t382 + t318 * t68 + t358 * t394 - t471 * t69;
t239 = qJD(4) + t251;
t231 = qJD(6) + t239;
t77 = t141 * t471 + t318 * t358;
t449 = t231 * t77;
t498 = -t23 + t449;
t456 = t77 ^ 2;
t487 = -t318 * t141 + t358 * t471;
t457 = t487 ^ 2;
t497 = -t456 + t457;
t455 = t77 * t487;
t266 = t313 * t322 + t315 * t319;
t154 = t266 * t251;
t479 = t266 * qJD(4) + t154;
t264 = t313 * t319 - t315 * t322;
t254 = t264 * qJD(4);
t401 = -t264 * t251 - t254;
t377 = qJD(2) * t454;
t248 = -t320 * qJD(3) - t323 * t377;
t185 = qJDD(2) * pkin(2) + qJD(1) * t248 - qJDD(1) * t278;
t247 = t323 * qJD(3) - t320 * t377;
t195 = qJD(1) * t247 + qJDD(1) * t279;
t122 = t437 * t185 - t314 * t195;
t120 = -qJDD(2) * pkin(3) - t122;
t310 = qJ(2) + pkin(10);
t299 = sin(t310);
t301 = cos(t310);
t321 = sin(qJ(1));
t324 = cos(qJ(1));
t363 = g(1) * t324 + g(2) * t321;
t341 = -g(3) * t301 + t363 * t299;
t496 = -qJD(4) * t239 * t292 - t120 + t341;
t458 = t323 * pkin(2);
t297 = pkin(1) + t458;
t276 = -qJD(1) * t297 + qJD(3);
t153 = t251 * pkin(3) - t255 * pkin(8) + t276;
t451 = qJD(2) * pkin(2);
t261 = -t272 + t451;
t376 = t437 * t273;
t190 = t314 * t261 + t376;
t177 = qJD(2) * pkin(8) + t190;
t101 = t322 * t153 - t177 * t319;
t86 = -qJ(5) * t210 + t101;
t72 = pkin(4) * t239 + t86;
t102 = t153 * t319 + t177 * t322;
t87 = -qJ(5) * t208 + t102;
t81 = t313 * t87;
t42 = t315 * t72 - t81;
t486 = pkin(9) * t358;
t30 = pkin(5) * t239 + t42 - t486;
t447 = t315 * t87;
t43 = t313 * t72 + t447;
t492 = pkin(9) * t141;
t32 = t43 - t492;
t253 = t267 * qJD(2);
t391 = t320 * qJDD(1);
t360 = -qJDD(1) * t375 + t314 * t391;
t191 = qJD(1) * t253 + t360;
t186 = qJDD(4) + t191;
t389 = pkin(2) * t381 + qJDD(3);
t390 = t323 * qJDD(1);
t435 = qJDD(1) * pkin(1);
t107 = -pkin(2) * t390 + t191 * pkin(3) - pkin(8) * t331 + t389 - t435;
t100 = t322 * t107;
t123 = t314 * t185 + t437 * t195;
t121 = qJDD(2) * pkin(8) + t123;
t39 = -qJD(4) * t102 - t319 * t121 + t100;
t27 = t186 * pkin(4) + t130 * qJ(5) - t210 * qJD(5) + t39;
t395 = qJD(4) * t322;
t38 = t319 * t107 + t322 * t121 + t153 * t395 - t177 * t396;
t29 = -qJ(5) * t131 - qJD(5) * t208 + t38;
t8 = t315 * t27 - t313 * t29;
t6 = pkin(5) * t186 - pkin(9) * t69 + t8;
t9 = t313 * t27 + t315 * t29;
t7 = -pkin(9) * t68 + t9;
t1 = t30 * t382 + t318 * t6 - t32 * t394 + t471 * t7;
t309 = qJ(4) + pkin(11);
t304 = qJ(6) + t309;
t290 = sin(t304);
t291 = cos(t304);
t412 = t301 * t321;
t212 = t290 * t324 - t291 * t412;
t411 = t301 * t324;
t214 = t290 * t321 + t291 * t411;
t461 = g(3) * t299;
t189 = t261 * t437 - t258;
t176 = -qJD(2) * pkin(3) - t189;
t136 = t208 * pkin(4) + qJD(5) + t176;
t73 = t141 * pkin(5) + t136;
t495 = g(1) * t214 - g(2) * t212 + t291 * t461 + t73 * t77 - t1;
t442 = -t313 * t500 + t499 * t315;
t441 = t499 * t313 + t315 * t500;
t24 = qJD(6) * t487 + t318 * t69 + t471 * t68;
t446 = t487 * t231;
t494 = -t24 + t446;
t11 = t318 * t30 + t32 * t471;
t2 = -qJD(6) * t11 - t318 * t7 + t471 * t6;
t211 = t290 * t412 + t291 * t324;
t213 = -t290 * t411 + t291 * t321;
t493 = -g(1) * t213 + g(2) * t211 + t290 * t461 - t487 * t73 + t2;
t444 = t264 * t382 + t266 * t394 + t318 * t479 - t401 * t471;
t194 = -t318 * t264 + t266 * t471;
t443 = qJD(6) * t194 + t318 * t401 + t471 * t479;
t491 = -pkin(5) * t255 - pkin(9) * t401 + t442;
t490 = -pkin(9) * t479 + t441;
t489 = t141 * t358;
t488 = -g(1) * t321 + g(2) * t324;
t347 = -t314 * t320 + t375;
t257 = t347 * qJD(2);
t383 = t267 * t395;
t349 = t257 * t319 + t383;
t485 = -t101 * t239 + t38;
t372 = t239 * t319;
t484 = t210 * t372;
t188 = -pkin(3) * t347 - pkin(8) * t267 - t297;
t203 = -t314 * t278 + t279 * t437;
t198 = t322 * t203;
t135 = t319 * t188 + t198;
t404 = t322 * t324;
t407 = t319 * t321;
t240 = t301 * t407 + t404;
t405 = t321 * t322;
t406 = t319 * t324;
t242 = -t301 * t406 + t405;
t477 = -g(1) * t242 + g(2) * t240;
t193 = t264 * t471 + t266 * t318;
t476 = t193 * t23 - t443 * t487;
t179 = qJDD(6) + t186;
t475 = t194 * t179 - t231 * t444;
t474 = -t266 * t186 - t239 * t401;
t339 = -t363 * t301 - t461;
t473 = t255 ^ 2;
t472 = t68 * pkin(5);
t469 = pkin(2) * t320;
t468 = pkin(4) * t313;
t282 = t324 * t297;
t463 = g(2) * t282;
t460 = g(3) * t323;
t459 = t319 * pkin(4);
t316 = -qJ(5) - pkin(8);
t262 = t403 * t319;
t263 = t403 * t322;
t180 = -t315 * t262 - t263 * t313;
t147 = -pkin(9) * t266 + t180;
t181 = -t313 * t262 + t315 * t263;
t148 = -pkin(9) * t264 + t181;
t88 = t147 * t471 - t318 * t148;
t453 = qJD(6) * t88 + t318 * t491 + t471 * t490;
t89 = t318 * t147 + t148 * t471;
t452 = -qJD(6) * t89 - t318 * t490 + t471 * t491;
t355 = -qJ(5) * t257 - qJD(5) * t267;
t165 = t247 * t437 + t314 * t248;
t388 = t320 * t451;
t167 = pkin(3) * t253 - pkin(8) * t257 + t388;
t374 = -t319 * t165 + t322 * t167;
t52 = t253 * pkin(4) + t355 * t322 + (-t198 + (qJ(5) * t267 - t188) * t319) * qJD(4) + t374;
t386 = t322 * t165 + t319 * t167 + t188 * t395;
t56 = -qJ(5) * t383 + (-qJD(4) * t203 + t355) * t319 + t386;
t21 = t313 * t52 + t315 * t56;
t48 = t315 * t86 - t81;
t448 = t255 * t77;
t445 = t487 * t255;
t414 = t267 * t319;
t108 = -qJ(5) * t414 + t135;
t134 = t322 * t188 - t203 * t319;
t413 = t267 * t322;
t93 = -pkin(4) * t347 - qJ(5) * t413 + t134;
t59 = t315 * t108 + t313 * t93;
t196 = -t272 * t314 + t376;
t144 = -pkin(4) * t418 + t196;
t440 = pkin(5) * t154 - t144 + (pkin(5) * t266 + t459) * qJD(4);
t293 = pkin(4) * t315 + pkin(5);
t245 = t293 * t471 - t318 * t468;
t47 = -t313 * t86 - t447;
t34 = t47 + t492;
t35 = t48 - t486;
t439 = t245 * qJD(6) - t318 * t34 - t35 * t471;
t246 = t318 * t293 + t468 * t471;
t438 = -t246 * qJD(6) + t318 * t35 - t34 * t471;
t436 = pkin(7) * qJDD(1);
t433 = t102 * t239;
t432 = t130 * t319;
t431 = t131 * t322;
t430 = t358 ^ 2;
t429 = t358 * t239;
t428 = t141 ^ 2;
t427 = t141 * t255;
t426 = t141 * t239;
t425 = t358 * t255;
t423 = t208 * t251;
t422 = t208 * t255;
t421 = t210 * t208;
t420 = t210 * t255;
t419 = t239 * t255;
t417 = t255 * t251;
t128 = t319 * t131;
t408 = t319 * t186;
t173 = t322 * t186;
t402 = -t208 * t395 - t128;
t298 = sin(t309);
t274 = pkin(5) * t298 + t459;
t400 = t274 + t454;
t300 = cos(t309);
t306 = t322 * pkin(4);
t275 = pkin(5) * t300 + t306;
t311 = t320 ^ 2;
t312 = t323 ^ 2;
t399 = t311 - t312;
t398 = t311 + t312;
t326 = qJD(1) ^ 2;
t387 = t320 * t326 * t323;
t385 = t437 * pkin(2);
t380 = t454 + t459;
t20 = -t313 * t56 + t315 * t52;
t58 = -t108 * t313 + t315 * t93;
t164 = t247 * t314 - t437 * t248;
t202 = t437 * t278 + t279 * t314;
t371 = t239 * t322;
t370 = -t194 * t24 + t444 * t77;
t369 = t323 * t381;
t368 = t488 * t299;
t367 = -t193 * t179 - t231 * t443;
t294 = -t385 - pkin(3);
t366 = pkin(4) * t396 - t144;
t365 = -t141 * t401 - t266 * t68;
t364 = pkin(3) * t301 + pkin(8) * t299;
t361 = -t264 * t186 - t239 * t479;
t162 = pkin(4) * t414 + t202;
t359 = -t101 * t322 - t102 * t319;
t271 = pkin(3) + t275;
t308 = -pkin(9) + t316;
t357 = t271 * t301 - t299 * t308;
t296 = t306 + pkin(3);
t356 = t296 * t301 - t299 * t316;
t277 = -t306 + t294;
t119 = pkin(4) * t349 + t164;
t353 = t173 + (-t396 - t418) * t239;
t170 = t264 * t267;
t46 = -pkin(5) * t347 + pkin(9) * t170 + t58;
t169 = t266 * t267;
t51 = -pkin(9) * t169 + t59;
t18 = -t318 * t51 + t46 * t471;
t19 = t318 * t46 + t471 * t51;
t350 = -0.2e1 * pkin(1) * t393 - pkin(7) * qJDD(2);
t112 = -t318 * t169 - t170 * t471;
t348 = t257 * t322 - t267 * t396;
t345 = t176 * t239 - t186 * t292;
t237 = -qJDD(1) * t297 + t389;
t64 = t131 * pkin(4) + qJDD(5) + t120;
t325 = qJD(2) ^ 2;
t338 = -pkin(7) * t325 + 0.2e1 * t435 - t488;
t337 = pkin(1) * t326 + t363 - t436;
t334 = t69 * t264 + t358 * t479;
t333 = t64 - t341;
t249 = t251 ^ 2;
t243 = t301 * t404 + t407;
t241 = -t301 * t405 + t406;
t222 = t298 * t321 + t300 * t411;
t221 = -t298 * t411 + t300 * t321;
t220 = t298 * t324 - t300 * t412;
t219 = t298 * t412 + t300 * t324;
t207 = t264 * pkin(5) + t277;
t115 = -t186 * t347 + t239 * t253;
t114 = pkin(5) * t169 + t162;
t113 = pkin(4) * t210 + pkin(5) * t358;
t111 = t169 * t471 - t170 * t318;
t110 = qJD(4) * t169 + t257 * t264;
t109 = t254 * t267 - t257 * t266;
t63 = -qJD(4) * t135 + t374;
t62 = -t203 * t396 + t386;
t60 = -pkin(5) * t109 + t119;
t41 = qJD(6) * t112 - t109 * t471 - t318 * t110;
t40 = -t318 * t109 + t110 * t471 + t169 * t382 - t170 * t394;
t33 = t64 + t472;
t17 = pkin(9) * t109 + t21;
t16 = pkin(5) * t253 + pkin(9) * t110 + t20;
t10 = t30 * t471 - t318 * t32;
t4 = -qJD(6) * t19 + t16 * t471 - t318 * t17;
t3 = qJD(6) * t18 + t318 * t16 + t17 * t471;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t488, t363, 0, 0, qJDD(1) * t311 + 0.2e1 * t369, 0.2e1 * t320 * t390 - 0.2e1 * t393 * t399, qJDD(2) * t320 + t323 * t325, qJDD(1) * t312 - 0.2e1 * t369, qJDD(2) * t323 - t320 * t325, 0, t320 * t350 + t323 * t338, -t320 * t338 + t323 * t350, 0.2e1 * t398 * t436 - t363, -g(1) * (-pkin(1) * t321 + pkin(7) * t324) - g(2) * (pkin(1) * t324 + pkin(7) * t321) + (pkin(7) ^ 2 * t398 + pkin(1) ^ 2) * qJDD(1), t255 * t257 + t267 * t331, -t267 * t191 - t257 * t251 - t255 * t253 + t331 * t347, qJD(2) * t257 + qJDD(2) * t267, -t191 * t347 + t251 * t253, -qJD(2) * t253 + qJDD(2) * t347, 0, -t202 * qJDD(2) - t297 * t191 - t237 * t347 + t276 * t253 - t488 * t301 + (t251 * t469 - t164) * qJD(2), -t165 * qJD(2) - t203 * qJDD(2) + t237 * t267 + t255 * t388 + t276 * t257 - t297 * t331 + t368, -t122 * t267 + t123 * t347 + t164 * t255 - t165 * t251 - t189 * t257 - t190 * t253 - t203 * t191 + t202 * t331 - t363, t123 * t203 + t190 * t165 - t122 * t202 - t189 * t164 - t237 * t297 + t276 * t388 - g(1) * (-t321 * t297 + t324 * t454) - g(2) * (t321 * t454 + t282) -t130 * t413 + t210 * t348 (-t208 * t322 - t210 * t319) * t257 + (t432 - t431 + (t208 * t319 - t210 * t322) * qJD(4)) * t267, t130 * t347 + t173 * t267 + t210 * t253 + t239 * t348, t128 * t267 + t208 * t349, t131 * t347 - t208 * t253 - t239 * t349 - t267 * t408, t115, -g(1) * t241 - g(2) * t243 + t101 * t253 + t120 * t414 + t202 * t131 + t134 * t186 + t164 * t208 + t176 * t349 + t63 * t239 - t347 * t39, -g(1) * t240 - g(2) * t242 - t102 * t253 + t120 * t413 - t202 * t130 - t135 * t186 + t164 * t210 + t176 * t348 - t62 * t239 + t347 * t38, t134 * t130 - t135 * t131 - t62 * t208 - t63 * t210 + t359 * t257 + (-t38 * t319 - t39 * t322 + (t101 * t319 - t102 * t322) * qJD(4)) * t267 - t368, -t463 + t101 * t63 + t102 * t62 + t120 * t202 + t39 * t134 + t38 * t135 + t176 * t164 + (-g(1) * t454 - g(2) * t364) * t324 + (-g(1) * (-t297 - t364) - g(2) * t454) * t321, -t110 * t358 - t170 * t69, t109 * t358 + t110 * t141 - t169 * t69 + t170 * t68, -t110 * t239 - t170 * t186 + t253 * t358 - t347 * t69, -t109 * t141 + t169 * t68, t109 * t239 - t141 * t253 - t169 * t186 + t347 * t68, t115, -g(1) * t220 - g(2) * t222 - t109 * t136 + t119 * t141 + t162 * t68 + t169 * t64 + t186 * t58 + t20 * t239 + t253 * t42 - t347 * t8, -g(1) * t219 - g(2) * t221 - t110 * t136 + t119 * t358 + t162 * t69 - t170 * t64 - t186 * t59 - t21 * t239 - t253 * t43 + t347 * t9, t109 * t43 + t110 * t42 - t141 * t21 - t169 * t9 + t170 * t8 - t20 * t358 - t58 * t69 - t59 * t68 - t368, -t463 + t136 * t119 + t64 * t162 + t42 * t20 + t43 * t21 + t8 * t58 + t9 * t59 + (-g(1) * t380 - g(2) * t356) * t324 + (-g(1) * (-t297 - t356) - g(2) * t380) * t321, -t112 * t23 - t40 * t487, t111 * t23 - t112 * t24 + t40 * t77 - t41 * t487, t112 * t179 + t23 * t347 - t231 * t40 + t253 * t487, t111 * t24 + t41 * t77, -t111 * t179 - t231 * t41 + t24 * t347 - t253 * t77, -t179 * t347 + t231 * t253, -g(1) * t212 - g(2) * t214 + t10 * t253 + t111 * t33 + t114 * t24 + t179 * t18 - t2 * t347 + t231 * t4 + t41 * t73 + t60 * t77, -g(1) * t211 - g(2) * t213 + t1 * t347 - t11 * t253 + t112 * t33 - t114 * t23 - t179 * t19 - t231 * t3 - t40 * t73 + t487 * t60, -t1 * t111 + t10 * t40 - t11 * t41 - t112 * t2 + t18 * t23 - t19 * t24 - t3 * t77 - t4 * t487 - t368, -t463 + t1 * t19 + t10 * t4 + t11 * t3 + t33 * t114 + t2 * t18 + t73 * t60 + (-g(1) * t400 - g(2) * t357) * t324 + (-g(1) * (-t297 - t357) - g(2) * t400) * t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t387, t399 * t326, t391, t387, t390, qJDD(2), t320 * t337 - t460, g(3) * t320 + t323 * t337, 0, 0, t417, -t249 + t473 (t284 + t251) * qJD(2) + t336, -t417, -t360, qJDD(2), t196 * qJD(2) - t276 * t255 + (qJDD(2) * t437 - t251 * t397) * pkin(2) + t122 + t341, t197 * qJD(2) + t276 * t251 + (-qJDD(2) * t314 - t255 * t397) * pkin(2) - t123 - t339, -t191 * t470 - t331 * t385 - (-t190 + t196) * t255 + (-t189 + t197) * t251, t189 * t196 - t190 * t197 + (t437 * t122 - t460 + t123 * t314 + (-qJD(1) * t276 + t363) * t320) * pkin(2), t210 * t371 - t432 (-t130 - t423) * t322 - t484 + t402, t239 * t371 + t408 - t420, t208 * t372 - t431, t353 + t422, -t419, -t101 * t255 - t117 * t239 + t294 * t131 - t196 * t208 + t345 * t319 + t322 * t496, t102 * t255 + t118 * t239 - t294 * t130 - t196 * t210 - t319 * t496 + t345 * t322, t117 * t210 + t118 * t208 + (-t101 * t251 - t131 * t292 + t38 + (t210 * t292 - t101) * qJD(4)) * t322 + (-t102 * t251 - t130 * t292 - t39 + (t208 * t292 - t102) * qJD(4)) * t319 + t339, t120 * t294 - t102 * t118 - t101 * t117 - t176 * t196 - g(3) * (t364 + t458) + (qJD(4) * t359 - t39 * t319 + t38 * t322) * t292 + t363 * (pkin(3) * t299 - pkin(8) * t301 + t469) t69 * t266 + t358 * t401, -t334 + t365, -t425 - t474, t141 * t479 + t68 * t264, t361 + t427, -t419, t136 * t154 - t144 * t141 + t180 * t186 - t42 * t255 + t64 * t264 + t277 * t68 + t442 * t239 + t341 * t300 + (t136 * t266 + t141 * t459) * qJD(4), t136 * t401 - t181 * t186 - t239 * t441 + t43 * t255 + t64 * t266 + t277 * t69 - t298 * t341 + t358 * t366, -t141 * t441 - t180 * t69 - t181 * t68 - t9 * t264 - t8 * t266 - t358 * t442 - t401 * t42 - t43 * t479 + t339, t9 * t181 + t8 * t180 + t64 * t277 - g(3) * (t356 + t458) + t441 * t43 + t442 * t42 + t366 * t136 + t363 * (t296 * t299 + t301 * t316 + t469) -t194 * t23 - t444 * t487, t370 + t476, -t445 + t475, t193 * t24 + t443 * t77, t367 + t448, -t231 * t255, -t10 * t255 + t88 * t179 + t33 * t193 + t207 * t24 + t231 * t452 + t291 * t341 + t440 * t77 + t443 * t73, t11 * t255 - t89 * t179 + t33 * t194 - t207 * t23 - t231 * t453 - t290 * t341 + t440 * t487 - t444 * t73, -t1 * t193 + t10 * t444 - t11 * t443 - t194 * t2 + t23 * t88 - t24 * t89 - t452 * t487 - t453 * t77 + t339, t1 * t89 + t2 * t88 + t33 * t207 - g(3) * (t357 + t458) + t440 * t73 + t453 * t11 + t452 * t10 + t363 * (t271 * t299 + t301 * t308 + t469); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t255 * qJD(2) + t360 (t284 - t251) * qJD(2) + t336, -t249 - t473, t189 * t255 + t190 * t251 + t237 + t488, 0, 0, 0, 0, 0, 0, t353 - t422, -t239 ^ 2 * t322 - t408 - t420 (t130 - t423) * t322 + t484 + t402, -t176 * t255 + (t39 + t433) * t322 + t485 * t319 + t488, 0, 0, 0, 0, 0, 0, t361 - t427, -t425 + t474, t334 + t365, -t136 * t255 - t8 * t264 + t9 * t266 + t401 * t43 - t42 * t479 + t488, 0, 0, 0, 0, 0, 0, t367 - t448, -t445 - t475, t370 - t476, t1 * t194 - t10 * t443 - t11 * t444 - t2 * t193 - t73 * t255 + t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, -t208 ^ 2 + t210 ^ 2, t208 * t239 - t130, -t421, t210 * t239 - t255 * t395 - t319 * t392 - t328, t186, -t177 * t395 + t433 - t176 * t210 + t100 + (-qJD(4) * t153 - t121 + t461) * t319 + t477, g(1) * t243 - g(2) * t241 + t176 * t208 + t322 * t461 - t485, 0, 0, t489, -t428 + t430, t69 + t426, -t489, -t68 + t429, t186, t298 * t461 - g(1) * t221 + g(2) * t219 - t136 * t358 - t47 * t239 + (-t141 * t210 + t186 * t315) * pkin(4) + t8, t300 * t461 + g(1) * t222 - g(2) * t220 + t136 * t141 + t48 * t239 + (-t186 * t313 - t210 * t358) * pkin(4) - t9 (-t313 * t68 - t315 * t69) * pkin(4) + (t43 + t47) * t358 + (t48 - t42) * t141, -t42 * t47 - t43 * t48 + (-t136 * t210 + t9 * t313 + t8 * t315 + t319 * t461 + t477) * pkin(4), t455, t497, t498, -t455, t494, t179, -t113 * t77 + t245 * t179 + t231 * t438 + t493, -t113 * t487 - t246 * t179 - t231 * t439 + t495, t245 * t23 - t246 * t24 + (t11 - t438) * t487 + (-t10 - t439) * t77, t1 * t246 + t2 * t245 - t73 * t113 - g(1) * (-t274 * t411 + t275 * t321) - g(2) * (-t274 * t412 - t275 * t324) + t274 * t461 + t439 * t11 + t438 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 + t429, t69 - t426, -t428 - t430, t141 * t43 + t358 * t42 + t333, 0, 0, 0, 0, 0, 0, t24 + t446, -t23 - t449, -t456 - t457, t10 * t487 + t11 * t77 + t333 + t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t455, t497, t498, -t455, t494, t179, t11 * t231 + t493, t10 * t231 + t495, 0, 0;];
tau_reg  = t5;