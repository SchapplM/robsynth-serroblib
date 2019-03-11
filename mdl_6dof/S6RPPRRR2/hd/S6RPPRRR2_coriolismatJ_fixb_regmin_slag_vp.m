% Calculate minimal parameter regressor of coriolis matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPPRRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:30
% EndTime: 2019-03-09 02:21:40
% DurationCPUTime: 5.21s
% Computational Cost: add. (6195->338), mult. (12872->469), div. (0->0), fcn. (14861->10), ass. (0->299)
t287 = sin(pkin(11));
t288 = cos(pkin(11));
t292 = sin(qJ(4));
t449 = cos(qJ(4));
t266 = t292 * t287 - t288 * t449;
t475 = -t266 / 0.2e1;
t279 = sin(pkin(10)) * pkin(1) + qJ(3);
t443 = pkin(7) + t279;
t255 = t443 * t287;
t256 = t443 * t288;
t176 = t449 * t255 + t256 * t292;
t344 = t449 * t287;
t391 = t292 * t288;
t268 = t344 + t391;
t291 = sin(qJ(5));
t191 = t291 * t268;
t124 = pkin(5) * t191 + t176;
t293 = cos(qJ(6));
t390 = t293 * t291;
t290 = sin(qJ(6));
t294 = cos(qJ(5));
t396 = t290 * t294;
t271 = t390 + t396;
t129 = t268 * t271;
t282 = -pkin(5) * t294 - pkin(4);
t454 = pkin(8) + pkin(9);
t273 = t454 * t291;
t274 = t454 * t294;
t317 = -t293 * t273 - t290 * t274;
t389 = t293 * t294;
t397 = t290 * t291;
t269 = -t389 + t397;
t453 = -t269 / 0.2e1;
t474 = -t282 * t129 / 0.2e1 + t124 * t453 + t317 * t475;
t451 = t271 / 0.2e1;
t351 = qJD(5) + qJD(6);
t175 = -t191 * t290 + t268 * t389;
t222 = -t290 * t273 + t293 * t274;
t472 = t282 * t175 / 0.2e1 + t124 * t451 + t222 * t475;
t190 = t291 * t266;
t160 = t176 * t294;
t445 = t268 * pkin(4);
t446 = t266 * pkin(8);
t207 = t445 + t446;
t199 = t291 * t207;
t380 = t160 - t199;
t91 = pkin(9) * t190 - t380;
t428 = t293 * t91;
t345 = -t428 / 0.2e1;
t194 = t294 * t266;
t200 = t294 * t207;
t407 = t176 * t291;
t444 = t268 * pkin(5);
t71 = pkin(9) * t194 + t200 + t407 + t444;
t435 = t290 * t71;
t314 = -t435 / 0.2e1 + t345;
t16 = t314 - t474;
t353 = t268 * qJD(1);
t209 = t266 * t353;
t259 = t344 / 0.2e1 + t391 / 0.2e1;
t464 = t259 * qJD(5) + t209;
t432 = t290 * t91;
t346 = -t432 / 0.2e1;
t431 = t293 * t71;
t313 = t346 + t431 / 0.2e1;
t15 = t313 - t472;
t471 = t351 * t222;
t470 = t351 * t317;
t469 = pkin(5) / 0.2e1;
t466 = t129 * t351;
t465 = t269 * t351;
t181 = t190 * qJD(5);
t258 = t268 * qJD(4);
t248 = t294 * t258;
t463 = -t248 + t181;
t260 = t266 ^ 2;
t261 = t268 ^ 2;
t462 = -t261 - t260;
t350 = t261 - t260;
t310 = t396 / 0.2e1 + t390 / 0.2e1;
t115 = (-t271 / 0.2e1 + t310) * t266;
t366 = t115 * qJD(2);
t371 = qJD(4) * t282;
t461 = qJD(1) * t15 - t271 * t371 + t366;
t448 = pkin(5) * t291;
t186 = t269 * t448 + t271 * t282;
t395 = t291 * t129;
t450 = -t294 / 0.2e1;
t8 = (-t395 / 0.2e1 + (t293 / 0.2e1 + t269 * t450) * t268) * pkin(5) + t15;
t460 = t8 * qJD(1) - t186 * qJD(4) + t366;
t459 = qJD(6) * t259 + t464;
t309 = t389 / 0.2e1 - t397 / 0.2e1;
t452 = t269 / 0.2e1;
t118 = (t452 + t309) * t266;
t364 = t118 * qJD(2);
t458 = qJD(1) * t16 + t269 * t371 + t364;
t187 = -t269 * t282 + t271 * t448;
t394 = t291 * t175;
t7 = (-t394 / 0.2e1 + (-t290 / 0.2e1 + t271 * t450) * t268) * pkin(5) + t16;
t457 = t7 * qJD(1) - t187 * qJD(4) + t364;
t285 = t291 ^ 2;
t286 = t294 ^ 2;
t278 = t286 - t285;
t392 = t291 * t294;
t332 = 0.2e1 * t268 * t392;
t299 = qJD(1) * t332 - t278 * qJD(4);
t447 = t266 * pkin(5);
t401 = t268 * t294;
t272 = -cos(pkin(10)) * pkin(1) - pkin(3) * t288 - pkin(2);
t323 = pkin(4) * t266 - pkin(8) * t268;
t167 = t272 + t323;
t177 = -t292 * t255 + t256 * t449;
t393 = t291 * t177;
t92 = -t294 * t167 + t393;
t83 = -pkin(9) * t401 - t92;
t60 = t83 + t447;
t58 = t293 * t60;
t456 = -t58 / 0.2e1;
t455 = -t60 / 0.2e1;
t44 = t129 * t269 - t271 * t175;
t442 = t351 * t44;
t116 = (t451 + t310) * t266;
t441 = t116 * qJD(3);
t82 = -t129 * t451 + t175 * t453;
t440 = t351 * t82;
t439 = pkin(5) * qJD(5);
t438 = pkin(5) * qJD(6);
t125 = -pkin(5) * t190 + t177;
t388 = t294 * t177;
t93 = t291 * t167 + t388;
t84 = -pkin(9) * t191 + t93;
t433 = t290 * t84;
t27 = -t58 + t433;
t400 = t271 * t266;
t1 = (t431 - t432) * t266 - t27 * t268 + t125 * t129 - t124 * t400;
t437 = t1 * qJD(1);
t174 = t269 * t266;
t429 = t293 * t84;
t28 = t290 * t60 + t429;
t2 = -(t428 + t435) * t266 - t28 * t268 + t125 * t175 + t124 * t174;
t436 = t2 * qJD(1);
t434 = t290 * t83;
t430 = t293 * t83;
t348 = -t447 / 0.2e1;
t325 = t348 + t83 / 0.2e1;
t3 = (t455 + t325) * t290;
t427 = t3 * qJD(1);
t5 = t293 * t325 + t456;
t426 = t5 * qJD(1);
t119 = (t453 + t309) * t266;
t425 = t119 * qJD(3);
t335 = t266 * t452;
t120 = t266 * t309 + t335;
t424 = t120 * qJD(3);
t298 = t310 * t266;
t122 = -t400 / 0.2e1 + t298;
t423 = t122 * qJD(3);
t31 = -t429 - t434;
t349 = pkin(5) * t401;
t412 = t124 * t175;
t13 = t129 * t349 + t31 * t266 + t412;
t422 = qJD(1) * t13;
t32 = t430 - t433;
t413 = t124 * t129;
t14 = t175 * t349 - t32 * t266 - t413;
t421 = qJD(1) * t14;
t19 = t266 * t27 - t413;
t420 = qJD(1) * t19;
t20 = -t266 * t28 + t412;
t419 = qJD(1) * t20;
t51 = t176 * t401 - t93 * t266;
t418 = qJD(1) * t51;
t402 = t268 * t129;
t403 = t266 * t400;
t78 = -t402 + t403;
t417 = qJD(1) * t78;
t79 = t402 + t403;
t416 = qJD(1) * t79;
t408 = t175 * t268;
t409 = t174 * t266;
t80 = t408 + t409;
t415 = qJD(1) * t80;
t21 = t200 * t266 + (-t92 + t393) * t268;
t406 = t21 * qJD(1);
t45 = -t129 * t174 + t175 * t400;
t387 = t45 * qJD(1);
t50 = -t176 * t191 + t266 * t92;
t386 = t50 * qJD(1);
t81 = t408 - t409;
t385 = t81 * qJD(1);
t384 = t351 * t118;
t336 = t293 * t475;
t117 = t294 * t336 + t290 * t190 / 0.2e1 + t335;
t383 = t351 * t117;
t382 = t351 * t116;
t381 = t351 * t122;
t379 = t351 * t175;
t275 = t287 ^ 2 + t288 ^ 2;
t146 = t350 * t291;
t378 = qJD(1) * t146;
t147 = t462 * t291;
t377 = qJD(1) * t147;
t148 = t350 * t294;
t376 = qJD(1) * t148;
t375 = qJD(1) * t175;
t374 = qJD(3) * t294;
t373 = qJD(4) * t194;
t372 = qJD(4) * t271;
t370 = qJD(4) * t294;
t369 = qJD(5) * t291;
t368 = qJD(5) * t294;
t367 = qJD(6) * t282;
t98 = t116 * qJD(1);
t99 = t117 * qJD(1);
t365 = t117 * qJD(4);
t39 = t118 * qJD(4);
t102 = t119 * qJD(1);
t127 = t268 * t269;
t363 = t127 * qJD(1);
t362 = t129 * qJD(1);
t361 = t350 * qJD(1);
t180 = t190 * qJD(1);
t360 = t191 * qJD(1);
t359 = t194 * qJD(1);
t198 = t462 * t294;
t358 = t198 * qJD(1);
t223 = t275 * t279;
t357 = t223 * qJD(1);
t356 = t259 * qJD(1);
t354 = t266 * qJD(1);
t257 = t266 * qJD(4);
t352 = t275 * qJD(1);
t347 = t444 / 0.2e1;
t343 = t175 * t354;
t342 = t291 * t370;
t341 = t268 * t369;
t340 = t268 * t368;
t208 = t266 * t258;
t339 = t269 * t258;
t338 = t291 * t368;
t337 = t294 * t353;
t334 = -t160 / 0.2e1 + t199 / 0.2e1;
t333 = pkin(5) * t351;
t330 = t351 * t266;
t329 = t351 * t271;
t328 = t294 * t347;
t327 = qJD(1) * t272 + qJD(3);
t326 = -qJD(5) - t354;
t324 = t347 + t71 / 0.2e1;
t322 = qJD(4) * t332;
t22 = (-t93 + t388) * t268 + (t380 - t160) * t266;
t320 = t22 * qJD(1);
t59 = t129 ^ 2 - t175 ^ 2;
t23 = qJD(1) * t59 + qJD(4) * t44;
t168 = t269 ^ 2 - t271 ^ 2;
t35 = qJD(1) * t44 + qJD(4) * t168;
t316 = t326 * t294;
t315 = t446 / 0.2e1 + t445 / 0.2e1;
t302 = t315 * t294;
t48 = -t200 / 0.2e1 - t302;
t312 = pkin(4) * t291 * qJD(4) - t48 * qJD(1);
t297 = t315 * t291 + t160 / 0.2e1;
t46 = t297 + t334;
t311 = pkin(4) * t370 - t46 * qJD(1);
t37 = qJD(4) * t82 - t129 * t375;
t53 = -qJD(1) * t82 + t269 * t372;
t306 = t268 * t316;
t188 = (t285 / 0.2e1 - t286 / 0.2e1) * t268;
t305 = -t188 * qJD(1) + t342;
t301 = qJD(1) * t261 * t392 + t188 * qJD(4);
t197 = t278 * t261;
t300 = t197 * qJD(1) + t322;
t251 = t259 * qJD(4);
t247 = t291 * t258;
t218 = t271 * t258;
t185 = t194 * qJD(5);
t179 = t190 * qJD(4);
t178 = t188 * qJD(5);
t154 = -t180 - t369;
t121 = t400 / 0.2e1 + t298;
t87 = -t329 - t98;
t86 = -t102 + t465;
t85 = -t99 - t465;
t55 = qJD(4) * t116 + t343;
t54 = t129 * t354 + t365;
t49 = t407 + t200 / 0.2e1 - t302;
t47 = t297 - t334;
t40 = t115 * qJD(4);
t34 = qJD(4) * t122 - t343 - t379;
t33 = -t39 + (-qJD(6) + t326) * t129;
t26 = qJD(4) * t121 - t379;
t25 = -t365 + t466;
t18 = t313 + t472;
t17 = t314 + t474;
t10 = t271 * t328 - t324 * t290 + t394 * t469 + t345 + t474;
t9 = t269 * t328 + t324 * t293 + t395 * t469 + t346 + t472;
t6 = pkin(5) * t336 + t433 + t456 - t430 / 0.2e1;
t4 = -t429 - t434 / 0.2e1 + (t348 + t455) * t290;
t11 = [0, 0, 0, 0, 0, 0, t275 * qJD(3), t223 * qJD(3), -t208, -t350 * qJD(4), 0, 0, 0, t272 * t258, -t272 * t257, -t208 * t286 - t261 * t338, -t197 * qJD(5) + t266 * t322, qJD(4) * t148 - t266 * t341, -t146 * qJD(4) - t266 * t340, t208, -qJD(3) * t147 + qJD(4) * t21 + qJD(5) * t51, -qJD(3) * t198 + qJD(4) * t22 + qJD(5) * t50 (qJD(4) * t174 - t466) * t175, t45 * qJD(4) + t351 * t59, t80 * qJD(4) - t129 * t330, t78 * qJD(4) - t175 * t330, t208, qJD(3) * t79 + qJD(4) * t1 + qJD(5) * t13 + qJD(6) * t20, qJD(3) * t81 + qJD(4) * t2 + qJD(5) * t14 + qJD(6) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t352, t357, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377, -t358, 0, 0, 0, 0, 0, t381 + t416, t120 * t351 + t385; 0, 0, 0, 0, 0, 0, 0, 0, -t209, -t361, -t257, -t258, 0, -qJD(4) * t177 + t272 * t353, qJD(4) * t176 - t272 * t354, -t178 + (-t286 * t353 - t342) * t266, t266 * t299 - 0.2e1 * t268 * t338, t247 + t376, t248 - t378, t464, t406 + (t291 * t323 - t388) * qJD(4) + t49 * qJD(5) (t294 * t323 + t393) * qJD(4) + t47 * qJD(5) + t320 (t372 + t375) * t174 + t440, t387 + (-t174 * t269 + t271 * t400) * qJD(4) + t442, t218 - t384 + t415, -t339 + t381 + t417, t459, t437 + (t125 * t269 + t268 * t317 - t282 * t400) * qJD(4) + t9 * qJD(5) + t18 * qJD(6), t436 + (t125 * t271 + t174 * t282 - t222 * t268) * qJD(4) + t10 * qJD(5) + t17 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, -t300, t326 * t191, t306, t251, qJD(4) * t49 - qJD(5) * t93 + t418, qJD(4) * t47 + qJD(5) * t92 + t386, t37, t23, t33, t34, t251, qJD(4) * t9 + qJD(5) * t31 + qJD(6) * t4 + t422 + t423, qJD(4) * t10 - qJD(5) * t32 + qJD(6) * t6 + t421 + t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t23, t33, t34, t251, qJD(4) * t18 + qJD(5) * t4 - qJD(6) * t28 + t419 + t423, qJD(4) * t17 + qJD(5) * t6 + qJD(6) * t27 + t420 + t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, t257, 0, 0, 0, 0, 0, t463, t185 + t247, 0, 0, 0, 0, 0, t121 * t351 + t339, t218 - t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179 - t340, t341 + t373, 0, 0, 0, 0, 0, t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t25; 0, 0, 0, 0, 0, 0, -t352, -t357, 0, 0, 0, 0, 0, t258, -t257, 0, 0, 0, 0, 0, t377 - t463, -t191 * qJD(4) - t266 * t368 + t358, 0, 0, 0, 0, 0, -qJD(4) * t127 - t382 - t416, -t129 * qJD(4) - t119 * t351 - t385; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, -t354, 0, 0, 0, 0, 0, t337, -t360, 0, 0, 0, 0, 0, -t363, -t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t316, 0, 0, 0, 0, 0, t87, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t86; 0, 0, 0, 0, 0, 0, 0, 0, t209, t361, 0, 0, 0, -t327 * t268, t327 * t266, t209 * t286 - t178, 0.2e1 * t291 * t306, t185 - t376, -t181 + t378, -t464, t48 * qJD(5) - t268 * t374 - t406, qJD(3) * t191 + qJD(5) * t46 - t320, -t174 * t375 + t440, -t387 + t442, -t383 - t415, -t382 - t417, -t459, qJD(3) * t127 - qJD(5) * t8 - qJD(6) * t15 - t437, qJD(3) * t129 - qJD(5) * t7 - qJD(6) * t16 - t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t351 * t115, -t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, t354, 0, 0, 0, 0, 0, -t337, t360, 0, 0, 0, 0, 0, t363, t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t278 * qJD(5), 0, 0, 0, -pkin(4) * t369, -pkin(4) * t368, -t269 * t329, t351 * t168, 0, 0, 0, qJD(5) * t186 + t271 * t367, qJD(5) * t187 - t269 * t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, -t299, t359 + t368, t154, -t356, -pkin(8) * t368 - t312, pkin(8) * t369 - t311, -t53, t35, t85, t87, -t356, -t460 - t471, -t457 - t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t35, t85, t87, -t356, -t461 - t471, -t458 - t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, t300, t209 * t291 - t373, t266 * t337 + t179, t251, qJD(3) * t190 - qJD(4) * t48 - t418, -t46 * qJD(4) + t266 * t374 - t386, -t37, -t23, t54, t55, t251, qJD(4) * t8 + qJD(6) * t3 - t422 + t441, qJD(4) * t7 + qJD(6) * t5 - t421 + t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t294 * t354, 0, 0, 0, 0, 0, t98, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t305, t299, -t359, t180, t356, t312, t311, t53, -t35, t99, t98, t356, t460, t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290 * t438, -t293 * t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290 * t333 + t427, -t293 * t333 + t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t23, t54, t55, t251, qJD(4) * t15 - qJD(5) * t3 - t419 + t441, qJD(4) * t16 - qJD(5) * t5 - t420 + t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t35, t99, t98, t356, t461, t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290 * t439 - t427, t293 * t439 - t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t11;
