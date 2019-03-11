% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:30
% EndTime: 2019-03-09 14:05:58
% DurationCPUTime: 12.87s
% Computational Cost: add. (27589->614), mult. (67623->854), div. (0->0), fcn. (51553->10), ass. (0->280)
t324 = sin(pkin(11));
t325 = cos(pkin(11));
t328 = sin(qJ(4));
t331 = cos(qJ(4));
t284 = t324 * t331 + t325 * t328;
t266 = t284 * qJD(4);
t332 = cos(qJ(2));
t341 = t284 * t332;
t389 = qJD(1) * t341 - t266;
t398 = t331 * t325;
t283 = t324 * t328 - t398;
t383 = qJD(1) * t332;
t377 = qJD(4) * t331;
t378 = qJD(4) * t328;
t428 = -t324 * t378 + t325 * t377;
t458 = t283 * t383 + t428;
t329 = sin(qJ(2));
t384 = qJD(1) * t329;
t368 = t324 * t384;
t381 = qJD(2) * t325;
t276 = t368 - t381;
t366 = t325 * t384;
t382 = qJD(2) * t324;
t278 = t366 + t382;
t216 = t276 * t328 - t278 * t331;
t217 = t331 * t276 + t278 * t328;
t327 = sin(qJ(5));
t330 = cos(qJ(5));
t137 = t216 * t330 + t217 * t327;
t141 = t216 * t327 - t330 * t217;
t326 = sin(qJ(6));
t421 = cos(qJ(6));
t79 = t137 * t421 - t326 * t141;
t462 = t79 ^ 2;
t84 = t326 * t137 + t141 * t421;
t461 = t84 ^ 2;
t311 = -qJD(4) + t383;
t304 = -qJD(5) + t311;
t297 = -qJD(6) + t304;
t460 = t297 * t84;
t459 = t79 * t297;
t349 = pkin(2) * t329 - qJ(3) * t332;
t286 = t349 * qJD(1);
t235 = pkin(7) * t368 + t325 * t286;
t400 = t325 * t332;
t345 = pkin(3) * t329 - pkin(8) * t400;
t200 = qJD(1) * t345 + t235;
t267 = t324 * t286;
t401 = t325 * t329;
t402 = t324 * t332;
t342 = -pkin(7) * t401 - pkin(8) * t402;
t223 = qJD(1) * t342 + t267;
t418 = pkin(8) + qJ(3);
t294 = t418 * t324;
t295 = t418 * t325;
t233 = -t328 * t294 + t331 * t295;
t393 = t284 * qJD(3) + qJD(4) * t233 + t331 * t200 - t223 * t328;
t392 = -qJD(3) * t398 + t331 * t223 + t294 * t377 + (qJD(3) * t324 + qJD(4) * t295 + t200) * t328;
t457 = -pkin(4) * t384 - pkin(9) * t458 - t393;
t456 = -pkin(9) * t389 + t392;
t455 = t137 ^ 2;
t454 = t141 ^ 2;
t419 = t84 * t79;
t453 = t137 * t304;
t452 = t141 * t304;
t375 = qJD(5) * t330;
t376 = qJD(5) * t327;
t391 = t283 * t375 + t284 * t376 - t327 * t389 - t330 * t458;
t222 = -t283 * t327 + t284 * t330;
t390 = qJD(5) * t222 + t327 * t458 - t330 * t389;
t441 = -t461 + t462;
t291 = -pkin(2) * t332 - qJ(3) * t329 - pkin(1);
t269 = t291 * qJD(1);
t317 = pkin(7) * t383;
t298 = qJD(2) * qJ(3) + t317;
t225 = t325 * t269 - t298 * t324;
t169 = -pkin(3) * t383 - pkin(8) * t278 + t225;
t226 = t324 * t269 + t325 * t298;
t175 = -pkin(8) * t276 + t226;
t107 = t331 * t169 - t175 * t328;
t93 = pkin(9) * t216 + t107;
t88 = -pkin(4) * t311 + t93;
t108 = t169 * t328 + t175 * t331;
t94 = -pkin(9) * t217 + t108;
t90 = t327 * t94;
t38 = t330 * t88 - t90;
t431 = pkin(10) * t137;
t30 = t38 + t431;
t28 = -pkin(5) * t304 + t30;
t92 = t330 * t94;
t39 = t327 * t88 + t92;
t446 = pkin(10) * t141;
t31 = t39 + t446;
t412 = t326 * t31;
t12 = t28 * t421 - t412;
t370 = t421 * t31;
t13 = t326 * t28 + t370;
t451 = t12 * t84 - t13 * t79;
t362 = qJD(6) * t421;
t374 = qJD(6) * t326;
t373 = qJD(1) * qJD(2);
t361 = t332 * t373;
t352 = t324 * t361;
t379 = qJD(2) * t332;
t353 = t379 * t398;
t161 = t328 * (qJD(4) * t278 + t352) - qJD(1) * t353 + t276 * t377;
t338 = qJD(2) * t341;
t426 = qJD(4) * t216;
t162 = qJD(1) * t338 - t426;
t69 = t330 * t161 + t327 * t162 - t216 * t376 + t217 * t375;
t70 = -t327 * t161 + t330 * t162 - t216 * t375 - t217 * t376;
t22 = -t137 * t374 - t141 * t362 + t326 * t70 + t421 * t69;
t439 = -t22 + t460;
t314 = t329 * t373;
t260 = qJD(2) * t349 - t329 * qJD(3);
t246 = t260 * qJD(1);
t316 = pkin(7) * t384;
t289 = (qJD(3) - t316) * qJD(2);
t198 = t325 * t246 - t289 * t324;
t340 = t345 * qJD(2);
t165 = qJD(1) * t340 + t198;
t199 = t324 * t246 + t325 * t289;
t171 = -pkin(8) * t352 + t199;
t59 = -qJD(4) * t108 + t331 * t165 - t328 * t171;
t44 = pkin(4) * t314 + pkin(9) * t161 + t59;
t58 = t328 * t165 + t169 * t377 + t331 * t171 - t175 * t378;
t46 = -pkin(9) * t162 + t58;
t9 = -qJD(5) * t39 - t327 * t46 + t330 * t44;
t6 = pkin(5) * t314 + pkin(10) * t69 + t9;
t358 = -t327 * t44 - t330 * t46 - t88 * t375 + t94 * t376;
t7 = -pkin(10) * t70 - t358;
t337 = -t28 * t362 + t31 * t374 - t326 * t6 - t421 * t7;
t413 = qJD(2) * pkin(2);
t356 = qJD(3) - t413;
t290 = t316 + t356;
t234 = pkin(3) * t276 + t290;
t158 = pkin(4) * t217 + t234;
t98 = -pkin(5) * t141 + t158;
t438 = -t84 * t98 + t337;
t232 = -t331 * t294 - t295 * t328;
t193 = -pkin(9) * t284 + t232;
t194 = -pkin(9) * t283 + t233;
t415 = t193 * t375 - t194 * t376 + t327 * t457 - t456 * t330;
t122 = t327 * t193 + t330 * t194;
t414 = -qJD(5) * t122 + t456 * t327 + t330 * t457;
t2 = -qJD(6) * t13 - t326 * t7 + t421 * t6;
t425 = t79 * t98 + t2;
t23 = -qJD(6) * t79 - t326 * t69 + t421 * t70;
t422 = -t23 + t459;
t448 = t216 ^ 2;
t447 = t217 ^ 2;
t445 = -pkin(10) * t390 + t415;
t444 = -pkin(5) * t384 + pkin(10) * t391 + t414;
t407 = t141 * t137;
t443 = t216 * t311;
t442 = t217 * t311;
t440 = -t454 + t455;
t437 = -t69 + t452;
t436 = -t141 * t158 + t358;
t432 = -0.2e1 * t373;
t275 = t325 * t291;
t224 = -pkin(8) * t401 + t275 + (-pkin(7) * t324 - pkin(3)) * t332;
t309 = pkin(7) * t400;
t242 = t324 * t291 + t309;
t403 = t324 * t329;
t231 = -pkin(8) * t403 + t242;
t147 = t331 * t224 - t231 * t328;
t256 = t283 * t329;
t119 = -pkin(4) * t332 + pkin(9) * t256 + t147;
t148 = t328 * t224 + t331 * t231;
t255 = t284 * t329;
t123 = -pkin(9) * t255 + t148;
t74 = t327 * t119 + t330 * t123;
t367 = t324 * t383;
t270 = pkin(3) * t367 + t317;
t360 = -pkin(4) * t389 - t270;
t322 = t329 ^ 2;
t323 = t332 ^ 2;
t427 = qJD(1) * (t322 - 0.2e1 * t323);
t424 = t137 * t158 + t9;
t423 = -t70 + t453;
t121 = t330 * t193 - t194 * t327;
t101 = -pkin(10) * t222 + t121;
t221 = t330 * t283 + t284 * t327;
t102 = -pkin(10) * t221 + t122;
t50 = t101 * t421 - t326 * t102;
t417 = t50 * qJD(6) + t326 * t444 + t421 * t445;
t51 = t326 * t101 + t102 * t421;
t416 = -t51 * qJD(6) - t326 * t445 + t421 * t444;
t43 = t330 * t93 - t90;
t146 = -t326 * t221 + t222 * t421;
t411 = -t146 * qJD(6) + t326 * t391 - t390 * t421;
t410 = t221 * t362 + t222 * t374 + t326 * t390 + t391 * t421;
t315 = pkin(4) * t330 + pkin(5);
t42 = -t327 * t93 - t92;
t34 = t42 - t446;
t35 = t43 + t431;
t369 = t421 * t327;
t409 = -t326 * t35 + t421 * t34 + t315 * t374 - (-t327 * t362 + (-t326 * t330 - t369) * qJD(5)) * pkin(4);
t399 = t326 * t327;
t408 = t326 * t34 + t421 * t35 - t315 * t362 - (-t327 * t374 + (t330 * t421 - t399) * qJD(5)) * pkin(4);
t406 = t216 * t217;
t405 = t276 * t325;
t334 = qJD(1) ^ 2;
t404 = t323 * t334;
t397 = t332 * t334;
t333 = qJD(2) ^ 2;
t396 = t333 * t329;
t395 = t333 * t332;
t394 = pkin(5) * t390 + t360;
t380 = qJD(2) * t329;
t371 = pkin(7) * t380;
t229 = t325 * t260 + t324 * t371;
t310 = pkin(7) * t361;
t259 = pkin(3) * t352 + t310;
t318 = pkin(7) * t379;
t365 = t324 * t379;
t271 = pkin(3) * t365 + t318;
t287 = pkin(3) * t403 + t329 * pkin(7);
t385 = t322 - t323;
t372 = pkin(7) * t402;
t313 = -pkin(3) * t325 - pkin(2);
t357 = pkin(1) * t432;
t73 = t330 * t119 - t123 * t327;
t355 = t276 + t381;
t354 = -t278 + t382;
t351 = t332 * t314;
t128 = pkin(4) * t162 + t259;
t196 = t329 * t428 + t338;
t163 = pkin(4) * t196 + t271;
t228 = pkin(4) * t255 + t287;
t350 = -t290 + t356;
t189 = -t255 * t327 - t256 * t330;
t247 = pkin(4) * t283 + t313;
t347 = qJD(1) * t355;
t346 = qJD(1) * t354;
t52 = -pkin(5) * t332 - pkin(10) * t189 + t73;
t188 = t330 * t255 - t256 * t327;
t53 = -pkin(10) * t188 + t74;
t26 = -t326 * t53 + t421 * t52;
t27 = t326 * t52 + t421 * t53;
t117 = -t326 * t188 + t189 * t421;
t343 = t332 * t346;
t47 = pkin(5) * t70 + t128;
t195 = t266 * t329 + t328 * t365 - t353;
t190 = t340 + t229;
t248 = t324 * t260;
t205 = qJD(2) * t342 + t248;
t76 = -qJD(4) * t148 + t331 * t190 - t328 * t205;
t64 = pkin(4) * t380 + pkin(9) * t195 + t76;
t75 = t328 * t190 + t331 * t205 + t224 * t377 - t231 * t378;
t71 = -pkin(9) * t196 + t75;
t20 = t119 * t375 - t123 * t376 + t327 * t64 + t330 * t71;
t21 = -qJD(5) * t74 - t327 * t71 + t330 * t64;
t321 = t325 ^ 2;
t320 = t324 ^ 2;
t307 = t329 * t397;
t299 = -0.2e1 * t351;
t263 = pkin(4) * t369 + t326 * t315;
t262 = -pkin(4) * t399 + t315 * t421;
t241 = t275 - t372;
t236 = -pkin(7) * t366 + t267;
t230 = -t325 * t371 + t248;
t167 = pkin(5) * t221 + t247;
t145 = t221 * t421 + t222 * t326;
t135 = pkin(5) * t188 + t228;
t116 = t188 * t421 + t189 * t326;
t103 = -pkin(4) * t216 - pkin(5) * t137;
t100 = qJD(5) * t189 - t327 * t195 + t330 * t196;
t99 = t330 * t195 + t327 * t196 + t255 * t375 - t256 * t376;
t72 = pkin(5) * t100 + t163;
t33 = t117 * qJD(6) + t100 * t421 - t326 * t99;
t32 = t326 * t100 + t188 * t362 + t189 * t374 + t421 * t99;
t17 = -pkin(10) * t100 + t20;
t16 = pkin(5) * t380 + pkin(10) * t99 + t21;
t15 = t30 * t421 - t412;
t14 = -t326 * t30 - t370;
t4 = -t27 * qJD(6) + t16 * t421 - t326 * t17;
t3 = t26 * qJD(6) + t326 * t16 + t17 * t421;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t351, t385 * t432, t395, t299, -t396, 0, -pkin(7) * t395 + t329 * t357, pkin(7) * t396 + t332 * t357, 0, 0 (t278 * t325 + t321 * t384) * t379 (-t405 + (-t278 - 0.2e1 * t366) * t324) * t379 (t278 * t329 + t325 * t427) * qJD(2) (t276 * t324 + t320 * t384) * t379 (-t276 * t329 - t324 * t427) * qJD(2), t299 (-qJD(1) * t229 - t198) * t332 + ((pkin(7) * t276 + t290 * t324) * t332 + (t225 + (t241 + 0.2e1 * t372) * qJD(1)) * t329) * qJD(2) (qJD(1) * t230 + t199) * t332 + ((pkin(7) * t278 + t290 * t325) * t332 + (-t226 + (-t242 + 0.2e1 * t309) * qJD(1)) * t329) * qJD(2), -t229 * t278 - t230 * t276 + (-t198 * t325 - t199 * t324) * t329 + (-t225 * t325 - t226 * t324 + (-t241 * t325 - t242 * t324) * qJD(1)) * t379, t198 * t241 + t199 * t242 + t225 * t229 + t226 * t230 + (t290 + t316) * t318, t161 * t256 + t195 * t216, t161 * t255 + t162 * t256 + t195 * t217 + t196 * t216, t161 * t332 + t195 * t311 + (-qJD(1) * t256 - t216) * t380, t162 * t255 + t196 * t217, t162 * t332 + t196 * t311 + (-qJD(1) * t255 - t217) * t380 (-t311 - t383) * t380, t162 * t287 + t196 * t234 + t217 * t271 + t255 * t259 - t311 * t76 - t332 * t59 + (qJD(1) * t147 + t107) * t380, -t161 * t287 - t195 * t234 - t216 * t271 - t256 * t259 + t311 * t75 + t332 * t58 + (-qJD(1) * t148 - t108) * t380, t107 * t195 - t108 * t196 + t147 * t161 - t148 * t162 + t216 * t76 - t217 * t75 - t255 * t58 + t256 * t59, t107 * t76 + t108 * t75 + t147 * t59 + t148 * t58 + t234 * t271 + t259 * t287, t137 * t99 - t189 * t69, t100 * t137 - t141 * t99 + t188 * t69 - t189 * t70, t304 * t99 + t332 * t69 + (qJD(1) * t189 - t137) * t380, -t100 * t141 + t188 * t70, t100 * t304 + t70 * t332 + (-qJD(1) * t188 + t141) * t380 (-t304 - t383) * t380, t100 * t158 + t128 * t188 - t141 * t163 - t21 * t304 + t228 * t70 - t332 * t9 + (qJD(1) * t73 + t38) * t380, t128 * t189 - t137 * t163 - t158 * t99 + t20 * t304 - t228 * t69 - t332 * t358 + (-qJD(1) * t74 - t39) * t380, -t100 * t39 + t137 * t21 + t141 * t20 + t188 * t358 - t189 * t9 + t38 * t99 + t69 * t73 - t70 * t74, t128 * t228 + t158 * t163 + t20 * t39 + t21 * t38 - t358 * t74 + t73 * t9, -t117 * t22 + t32 * t79, t116 * t22 - t117 * t23 - t32 * t84 + t33 * t79, t22 * t332 + t297 * t32 + (qJD(1) * t117 - t79) * t380, t116 * t23 - t33 * t84, t23 * t332 + t33 * t297 + (-qJD(1) * t116 + t84) * t380 (-t297 - t383) * t380, t116 * t47 + t135 * t23 - t2 * t332 - t297 * t4 + t33 * t98 - t72 * t84 + (qJD(1) * t26 + t12) * t380, -t337 * t332 + t117 * t47 - t135 * t22 + t297 * t3 - t32 * t98 - t72 * t79 + (-qJD(1) * t27 - t13) * t380, t116 * t337 - t117 * t2 + t12 * t32 - t13 * t33 + t22 * t26 - t23 * t27 + t3 * t84 + t4 * t79, t12 * t4 + t13 * t3 + t135 * t47 + t2 * t26 - t27 * t337 + t72 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t385 * t334, 0, t307, 0, 0, t334 * pkin(1) * t329, pkin(1) * t397, 0, 0, t325 * t343 (t405 + t278 * t324 + (-t320 + t321) * qJD(2)) * t383, t325 * t404 + t329 * t346, -t355 * t367, -t324 * t404 + t329 * t347, t307 ((-qJ(3) * t382 - t225) * t329 + (-pkin(7) * t355 + t324 * t350 + t235) * t332) * qJD(1) ((-qJ(3) * t381 + t226) * t329 + (pkin(7) * t354 + t325 * t350 - t236) * t332) * qJD(1), t235 * t278 + t236 * t276 + (-qJD(3) * t276 + t225 * t383 + t199) * t325 + (qJD(3) * t278 + t226 * t383 - t198) * t324, -t225 * t235 - t226 * t236 + (-t225 * t324 + t226 * t325) * qJD(3) + (-t198 * t324 + t199 * t325) * qJ(3) + (-t290 - t413) * t317, -t161 * t284 - t216 * t458, t161 * t283 - t162 * t284 - t216 * t389 - t217 * t458, -t458 * t311 + (qJD(2) * t284 + t216) * t384, t162 * t283 - t217 * t389, -t389 * t311 + (-qJD(2) * t283 + t217) * t384, t311 * t384, t162 * t313 - t217 * t270 + t259 * t283 + t393 * t311 - t389 * t234 + (qJD(2) * t232 - t107) * t384, -t161 * t313 + t216 * t270 + t259 * t284 - t392 * t311 + t458 * t234 + (-qJD(2) * t233 + t108) * t384, -t107 * t458 + t108 * t389 + t161 * t232 - t162 * t233 - t216 * t393 + t217 * t392 - t283 * t58 - t284 * t59, -t107 * t393 - t108 * t392 + t232 * t59 + t233 * t58 - t234 * t270 + t259 * t313, t137 * t391 - t69 * t222, t137 * t390 - t141 * t391 + t221 * t69 - t222 * t70, t391 * t304 + (qJD(2) * t222 + t137) * t384, -t141 * t390 + t70 * t221, t390 * t304 + (-qJD(2) * t221 - t141) * t384, t304 * t384, t128 * t221 + t247 * t70 - t414 * t304 + t390 * t158 - t360 * t141 + (qJD(2) * t121 - t38) * t384, t128 * t222 - t247 * t69 + t415 * t304 - t391 * t158 - t360 * t137 + (-qJD(2) * t122 + t39) * t384, t121 * t69 - t122 * t70 + t137 * t414 + t141 * t415 + t221 * t358 - t222 * t9 + t38 * t391 - t39 * t390, t121 * t9 - t122 * t358 + t128 * t247 + t158 * t360 + t38 * t414 + t39 * t415, -t22 * t146 + t410 * t79, t145 * t22 - t146 * t23 - t410 * t84 - t411 * t79, t410 * t297 + (qJD(2) * t146 + t79) * t384, t23 * t145 + t411 * t84, -t411 * t297 + (-qJD(2) * t145 - t84) * t384, t297 * t384, t145 * t47 + t167 * t23 - t411 * t98 - t394 * t84 - t416 * t297 + (qJD(2) * t50 - t12) * t384, t146 * t47 - t167 * t22 - t410 * t98 - t394 * t79 + t417 * t297 + (-qJD(2) * t51 + t13) * t384, t12 * t410 + t13 * t411 + t145 * t337 - t146 * t2 + t22 * t50 - t23 * t51 + t416 * t79 + t417 * t84, t12 * t416 + t13 * t417 + t167 * t47 + t2 * t50 - t337 * t51 + t394 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, t332 * t347, -t276 ^ 2 - t278 ^ 2, t225 * t278 + t226 * t276 + t310, 0, 0, 0, 0, 0, 0, t162 + t443, -t161 + t442, -t447 - t448, -t107 * t216 + t108 * t217 + t259, 0, 0, 0, 0, 0, 0, t70 + t453, -t69 - t452, -t454 - t455, -t137 * t38 - t141 * t39 + t128, 0, 0, 0, 0, 0, 0, t23 + t459, -t22 - t460, -t461 - t462, -t12 * t79 - t13 * t84 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t406, -t447 + t448, -t161 - t442, t406, -t284 * t361 + t426 + t443, t314, -t108 * t311 + t216 * t234 + t59, -t107 * t311 + t217 * t234 - t58, 0, 0, t407, t440, t437, -t407, t423, t314, t304 * t42 + (-t141 * t216 + t304 * t376 + t314 * t330) * pkin(4) + t424, -t304 * t43 + (-t137 * t216 + t304 * t375 - t314 * t327) * pkin(4) + t436, -t137 * t39 - t141 * t43 + t141 * t38 - t137 * t42 + (-t327 * t70 + t330 * t69 + (-t137 * t327 + t141 * t330) * qJD(5)) * pkin(4), -t38 * t42 - t39 * t43 + (t158 * t216 - t327 * t358 + t330 * t9 + (-t327 * t38 + t330 * t39) * qJD(5)) * pkin(4), t419, t441, t439, -t419, t422, t314, t103 * t84 + t262 * t314 + t297 * t409 + t425, t103 * t79 - t263 * t314 - t297 * t408 + t438, t22 * t262 - t23 * t263 - t408 * t84 - t409 * t79 + t451, -t103 * t98 - t12 * t409 - t13 * t408 + t2 * t262 - t263 * t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t407, t440, t437, -t407, t423, t314, -t304 * t39 + t424, -t304 * t38 + t436, 0, 0, t419, t441, t439, -t419, t422, t314, t14 * t297 + (-t137 * t84 + t297 * t374 + t314 * t421) * pkin(5) + t425, -t15 * t297 + (-t137 * t79 + t297 * t362 - t314 * t326) * pkin(5) + t438, -t14 * t79 - t15 * t84 + (t421 * t22 - t23 * t326 + (-t326 * t79 + t421 * t84) * qJD(6)) * pkin(5) + t451, -t12 * t14 - t13 * t15 + (t421 * t2 - t337 * t326 + t137 * t98 + (-t12 * t326 + t13 * t421) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, t441, t439, -t419, t422, t314, -t13 * t297 + t425, -t12 * t297 + t438, 0, 0;];
tauc_reg  = t1;
