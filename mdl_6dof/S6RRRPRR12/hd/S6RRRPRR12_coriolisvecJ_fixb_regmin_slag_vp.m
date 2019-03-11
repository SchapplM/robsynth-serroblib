% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:30
% EndTime: 2019-03-09 19:47:01
% DurationCPUTime: 13.41s
% Computational Cost: add. (15928->617), mult. (41490->872), div. (0->0), fcn. (34005->12), ass. (0->269)
t312 = sin(pkin(6));
t321 = cos(qJ(2));
t389 = qJD(1) * t321;
t373 = t312 * t389;
t292 = -qJD(3) + t373;
t317 = sin(qJ(2));
t421 = cos(pkin(6));
t363 = t421 * qJD(1);
t352 = pkin(1) * t363;
t262 = pkin(8) * t373 + t317 * t352;
t316 = sin(qJ(3));
t320 = cos(qJ(3));
t455 = qJD(4) * t316 + t262 + t292 * (pkin(3) * t316 - qJ(4) * t320);
t337 = t363 + qJD(2);
t390 = qJD(1) * t312;
t374 = t317 * t390;
t442 = -t316 * t374 + t320 * t337;
t231 = qJD(5) - t442;
t222 = qJD(6) + t231;
t240 = t316 * t337 + t320 * t374;
t311 = sin(pkin(12));
t313 = cos(pkin(12));
t187 = t240 * t311 + t313 * t292;
t189 = t240 * t313 - t292 * t311;
t315 = sin(qJ(5));
t319 = cos(qJ(5));
t119 = t319 * t187 + t189 * t315;
t314 = sin(qJ(6));
t318 = cos(qJ(6));
t437 = t187 * t315 - t189 * t319;
t57 = t318 * t119 - t314 * t437;
t454 = t222 * t57;
t453 = t119 * t231;
t340 = t119 * t314 + t318 * t437;
t452 = t222 * t340;
t259 = -pkin(8) * t374 + t321 * t352;
t333 = (pkin(2) * t317 - pkin(9) * t321) * t312;
t260 = qJD(1) * t333;
t397 = t320 * t259 + t316 * t260;
t163 = qJ(4) * t374 + t397;
t386 = qJD(3) * t316;
t376 = pkin(9) * t386;
t403 = t455 * t313 + (-t163 - t376) * t311;
t451 = -t313 * t163 - t311 * t455;
t404 = t320 * t321;
t221 = (t311 * t317 + t313 * t404) * t390;
t353 = t316 * t373;
t405 = t313 * t320;
t450 = -pkin(4) * t353 + pkin(10) * t221 + (pkin(4) * t316 - pkin(10) * t405) * qJD(3) - t403;
t409 = t311 * t320;
t220 = -t313 * t374 + t373 * t409;
t406 = t313 * t316;
t443 = pkin(10) * t220 + (-pkin(9) * t406 - pkin(10) * t409) * qJD(3) + t451;
t218 = pkin(9) * t337 + t262;
t255 = (-pkin(2) * t321 - pkin(9) * t317 - pkin(1)) * t312;
t230 = qJD(1) * t255;
t146 = -t316 * t218 + t230 * t320;
t170 = pkin(3) * t240 - qJ(4) * t442;
t104 = t313 * t146 + t311 * t170;
t415 = t442 * t311;
t84 = -pkin(10) * t415 + t104;
t449 = -qJD(4) * t313 + t84;
t448 = t231 * t437;
t133 = pkin(3) * t292 + qJD(4) - t146;
t106 = pkin(4) * t187 + t133;
t54 = pkin(5) * t119 + t106;
t445 = t54 * t57;
t444 = t340 * t57;
t411 = t311 * t315;
t279 = -t319 * t313 + t411;
t280 = t311 * t319 + t313 * t315;
t382 = qJD(5) * t316;
t385 = qJD(3) * t320;
t401 = t220 * t315 - t221 * t319 - t279 * t385 - t280 * t382;
t381 = qJD(5) * t319;
t400 = -t319 * t220 - t221 * t315 + t280 * t385 + t381 * t406 - t382 * t411;
t396 = t231 * t279;
t395 = t231 * t280;
t441 = t340 ^ 2 - t57 ^ 2;
t217 = -pkin(2) * t337 - t259;
t131 = -pkin(3) * t442 - t240 * qJ(4) + t217;
t147 = t320 * t218 + t316 * t230;
t135 = -qJ(4) * t292 + t147;
t74 = t313 * t131 - t135 * t311;
t44 = -pkin(4) * t442 - pkin(10) * t189 + t74;
t75 = t311 * t131 + t313 * t135;
t55 = -pkin(10) * t187 + t75;
t20 = -t315 * t55 + t319 * t44;
t16 = pkin(11) * t437 + t20;
t14 = pkin(5) * t231 + t16;
t21 = t315 * t44 + t319 * t55;
t17 = -pkin(11) * t119 + t21;
t425 = t17 * t318;
t11 = t14 * t314 + t425;
t377 = qJD(1) * qJD(2);
t367 = t312 * t377;
t350 = t321 * t367;
t196 = qJD(3) * t240 + t316 * t350;
t195 = qJD(3) * t442 + t320 * t350;
t351 = t317 * t367;
t161 = t195 * t311 - t313 * t351;
t162 = t195 * t313 + t311 * t351;
t383 = qJD(5) * t315;
t52 = -t315 * t161 + t319 * t162 - t187 * t381 - t189 * t383;
t261 = qJD(2) * t333;
t248 = qJD(1) * t261;
t375 = pkin(1) * t421;
t408 = t312 * t317;
t433 = -pkin(8) * t408 + t321 * t375;
t263 = t433 * qJD(2);
t249 = qJD(1) * t263;
t329 = t218 * t386 - t230 * t385 - t316 * t248 - t320 * t249;
t90 = qJ(4) * t351 - qJD(4) * t292 - t329;
t407 = t312 * t321;
t326 = pkin(8) * t407 + t317 * t375;
t264 = t326 * qJD(2);
t250 = qJD(1) * t264;
t99 = t196 * pkin(3) - t195 * qJ(4) - t240 * qJD(4) + t250;
t40 = -t311 * t90 + t313 * t99;
t29 = pkin(4) * t196 - pkin(10) * t162 + t40;
t41 = t311 * t99 + t313 * t90;
t33 = -pkin(10) * t161 + t41;
t7 = -qJD(5) * t21 + t319 * t29 - t315 * t33;
t4 = pkin(5) * t196 - pkin(11) * t52 + t7;
t53 = -qJD(5) * t437 + t319 * t161 + t162 * t315;
t6 = t315 * t29 + t319 * t33 + t44 * t381 - t383 * t55;
t5 = -pkin(11) * t53 + t6;
t2 = -qJD(6) * t11 - t314 * t5 + t318 * t4;
t440 = t54 * t340 + t2;
t379 = qJD(6) * t318;
t380 = qJD(6) * t314;
t12 = -t119 * t379 - t314 * t53 + t318 * t52 + t380 * t437;
t439 = t12 + t454;
t323 = qJD(6) * t340 - t314 * t52 - t318 * t53;
t438 = t323 - t452;
t436 = t450 * t319;
t435 = t317 * t321;
t288 = -pkin(3) * t320 - qJ(4) * t316 - pkin(2);
t274 = t313 * t288;
t199 = -pkin(10) * t406 + t274 + (-pkin(9) * t311 - pkin(4)) * t320;
t236 = pkin(9) * t405 + t311 * t288;
t410 = t311 * t316;
t212 = -pkin(10) * t410 + t236;
t399 = t315 * t199 + t319 * t212;
t359 = -t316 * t259 + t260 * t320;
t164 = -pkin(3) * t374 - t359;
t394 = -pkin(4) * t220 - t164 + (pkin(4) * t311 + pkin(9)) * t385;
t431 = pkin(10) + qJ(4);
t290 = t431 * t311;
t291 = t431 * t313;
t392 = -t315 * t290 + t319 * t291;
t434 = -t199 * t381 + t212 * t383 - t315 * t450 - t443 * t319;
t15 = t17 * t380;
t366 = qJD(6) * t14 + t5;
t1 = t314 * t4 + t318 * t366 - t15;
t335 = qJD(4) * t311 + qJD(5) * t291;
t103 = -t146 * t311 + t313 * t170;
t71 = -pkin(10) * t313 * t442 + pkin(4) * t240 + t103;
t432 = t290 * t381 + t449 * t319 + (t335 + t71) * t315;
t322 = qJD(1) ^ 2;
t253 = -pkin(2) * t421 - t433;
t269 = t316 * t408 - t320 * t421;
t270 = t316 * t421 + t320 * t408;
t157 = t269 * pkin(3) - t270 * qJ(4) + t253;
t254 = pkin(9) * t421 + t326;
t398 = t320 * t254 + t316 * t255;
t158 = -qJ(4) * t407 + t398;
t101 = t313 * t157 - t158 * t311;
t209 = t270 * t313 - t311 * t407;
t66 = pkin(4) * t269 - pkin(10) * t209 + t101;
t102 = t311 * t157 + t313 * t158;
t208 = t270 * t311 + t313 * t407;
t76 = -pkin(10) * t208 + t102;
t430 = t315 * t66 + t319 * t76;
t257 = t280 * t316;
t258 = t279 * t316;
t172 = t318 * t257 - t258 * t314;
t428 = -qJD(6) * t172 - t314 * t400 + t318 * t401;
t173 = -t257 * t314 - t258 * t318;
t427 = qJD(6) * t173 + t314 * t401 + t318 * t400;
t197 = t318 * t279 + t280 * t314;
t424 = -qJD(6) * t197 - t314 * t395 - t318 * t396;
t198 = -t279 * t314 + t280 * t318;
t423 = qJD(6) * t198 - t314 * t396 + t318 * t395;
t422 = pkin(5) * t400 + t394;
t420 = qJ(4) * t196;
t356 = t218 * t385 + t230 * t386 - t320 * t248 + t316 * t249;
t100 = -pkin(3) * t351 + t356;
t419 = t100 * t311;
t418 = t100 * t313;
t417 = t196 * t320;
t416 = t442 * t292;
t414 = t240 * t292;
t330 = t292 * t316;
t413 = t292 * t320;
t308 = t312 ^ 2;
t412 = t308 * t322;
t328 = -t254 * t386 + t255 * t385 + t316 * t261 + t320 * t263;
t388 = qJD(2) * t317;
t107 = (qJ(4) * t388 - qJD(4) * t321) * t312 + t328;
t371 = qJD(2) * t407;
t210 = qJD(3) * t270 + t316 * t371;
t211 = -qJD(3) * t269 + t320 * t371;
t113 = t210 * pkin(3) - t211 * qJ(4) - t270 * qJD(4) + t264;
t46 = t313 * t107 + t311 * t113;
t402 = t313 * t376 - t451;
t281 = pkin(4) * t410 + t316 * pkin(9);
t391 = t317 ^ 2 - t321 ^ 2;
t387 = qJD(2) * t320;
t378 = -qJD(4) + t133;
t125 = pkin(4) * t415 + t147;
t305 = -pkin(4) * t313 - pkin(3);
t372 = t312 * t388;
t370 = pkin(5) * t395 - t125;
t368 = t308 * t377;
t364 = -t315 * t76 + t319 * t66;
t45 = -t107 * t311 + t313 * t113;
t361 = t319 * t199 - t212 * t315;
t360 = -t316 * t254 + t255 * t320;
t358 = -t319 * t290 - t291 * t315;
t357 = 0.2e1 * t368;
t114 = -pkin(5) * t320 + pkin(11) * t258 + t361;
t349 = pkin(11) * t400 - qJD(6) * t114 + t434;
t116 = -pkin(11) * t257 + t399;
t348 = t399 * qJD(5) + qJD(6) * t116 - t436 + t443 * t315 + t401 * pkin(11) + (t353 - t386) * pkin(5);
t347 = t312 * t322 * t421;
t346 = -0.2e1 * pkin(1) * t368;
t181 = -pkin(11) * t280 + t358;
t345 = pkin(11) * t395 - qJD(6) * t181 + t432;
t182 = -pkin(11) * t279 + t392;
t68 = t319 * t71;
t344 = pkin(5) * t240 - pkin(11) * t396 + t280 * qJD(4) + t392 * qJD(5) + qJD(6) * t182 - t315 * t84 + t68;
t159 = pkin(3) * t407 - t360;
t137 = -t208 * t315 + t209 * t319;
t22 = pkin(5) * t269 - pkin(11) * t137 + t364;
t136 = t319 * t208 + t209 * t315;
t24 = -pkin(11) * t136 + t430;
t342 = t22 * t318 - t24 * t314;
t341 = t22 * t314 + t24 * t318;
t82 = t318 * t136 + t137 * t314;
t83 = -t136 * t314 + t137 * t318;
t336 = 0.2e1 * t363 + qJD(2);
t334 = -t254 * t385 - t255 * t386 + t261 * t320 - t316 * t263;
t178 = t211 * t313 + t311 * t372;
t37 = pkin(4) * t210 - pkin(10) * t178 + t45;
t177 = t211 * t311 - t313 * t372;
t39 = -pkin(10) * t177 + t46;
t332 = t315 * t37 + t319 * t39 + t66 * t381 - t383 * t76;
t124 = pkin(4) * t208 + t159;
t325 = pkin(1) * (-qJD(2) * t363 + t412);
t112 = -pkin(3) * t372 - t334;
t324 = -qJD(5) * t430 - t315 * t39 + t319 * t37;
t61 = pkin(4) * t161 + t100;
t77 = pkin(4) * t177 + t112;
t247 = pkin(5) * t279 + t305;
t235 = -pkin(9) * t409 + t274;
t202 = pkin(5) * t257 + t281;
t160 = t196 * t269;
t72 = pkin(5) * t136 + t124;
t70 = qJD(5) * t137 + t319 * t177 + t178 * t315;
t69 = -qJD(5) * t136 - t177 * t315 + t178 * t319;
t34 = pkin(5) * t70 + t77;
t25 = pkin(5) * t53 + t61;
t19 = qJD(6) * t83 + t314 * t69 + t318 * t70;
t18 = -qJD(6) * t82 - t314 * t70 + t318 * t69;
t10 = t14 * t318 - t17 * t314;
t9 = -pkin(11) * t70 + t332;
t8 = pkin(5) * t210 - pkin(11) * t69 + t324;
t3 = [0, 0, 0, t357 * t435, -t391 * t357, t336 * t371, -t336 * t372, 0, -t250 * t421 - t264 * t337 + t317 * t346, -t249 * t421 - t263 * t337 + t321 * t346, t195 * t270 + t211 * t240, -t195 * t269 - t196 * t270 - t210 * t240 + t211 * t442, -t211 * t292 + (-t195 * t321 + (qJD(1) * t270 + t240) * t388) * t312, t210 * t292 + (t196 * t321 + (-qJD(1) * t269 + t442) * t388) * t312 (-t292 * t312 - t308 * t389) * t388, -t334 * t292 - t264 * t442 + t253 * t196 + t250 * t269 + t217 * t210 + (t356 * t321 + (qJD(1) * t360 + t146) * t388) * t312, t328 * t292 + t264 * t240 + t253 * t195 + t250 * t270 + t217 * t211 + (-t329 * t321 + (-qJD(1) * t398 - t147) * t388) * t312, t100 * t208 + t101 * t196 + t112 * t187 + t133 * t177 + t159 * t161 + t210 * t74 + t269 * t40 - t442 * t45, t100 * t209 - t102 * t196 + t112 * t189 + t133 * t178 + t159 * t162 - t210 * t75 - t269 * t41 + t442 * t46, -t101 * t162 - t102 * t161 - t177 * t75 - t178 * t74 - t187 * t46 - t189 * t45 - t208 * t41 - t209 * t40, t100 * t159 + t101 * t40 + t102 * t41 + t112 * t133 + t45 * t74 + t46 * t75, t137 * t52 - t437 * t69, -t119 * t69 - t136 * t52 - t137 * t53 + t437 * t70, t137 * t196 - t210 * t437 + t231 * t69 + t269 * t52, -t119 * t210 - t136 * t196 - t231 * t70 - t269 * t53, t210 * t231 + t160, t106 * t70 + t77 * t119 + t124 * t53 + t61 * t136 + t196 * t364 + t20 * t210 + t231 * t324 + t7 * t269, t106 * t69 + t124 * t52 + t61 * t137 - t196 * t430 - t21 * t210 - t231 * t332 - t6 * t269 - t437 * t77, t12 * t83 - t18 * t340, -t12 * t82 - t18 * t57 + t19 * t340 + t323 * t83, t12 * t269 + t18 * t222 + t196 * t83 - t210 * t340, -t19 * t222 - t196 * t82 - t210 * t57 + t269 * t323, t210 * t222 + t160 (-qJD(6) * t341 - t314 * t9 + t318 * t8) * t222 + t342 * t196 + t2 * t269 + t10 * t210 + t34 * t57 - t72 * t323 + t25 * t82 + t54 * t19 -(qJD(6) * t342 + t314 * t8 + t318 * t9) * t222 - t341 * t196 - t1 * t269 - t11 * t210 - t34 * t340 + t72 * t12 + t25 * t83 + t54 * t18; 0, 0, 0, -t412 * t435, t391 * t412, -t321 * t347, t317 * t347, 0, -pkin(8) * t350 + t262 * t337 + t317 * t325, pkin(8) * t351 + t259 * t337 + t321 * t325, t195 * t316 - t240 * t413 (t195 - t416) * t320 + (-t196 + t414) * t316, -t292 * t385 + (t292 * t404 + (qJD(2) * t316 - t240) * t317) * t390, t292 * t386 + (-t321 * t330 + (-t442 + t387) * t317) * t390, t292 * t374, -pkin(2) * t196 - t250 * t320 + t359 * t292 + t262 * t442 + (pkin(9) * t413 + t217 * t316) * qJD(3) + (-t146 * t317 + (-pkin(9) * t388 - t217 * t321) * t316) * t390, -pkin(2) * t195 + t250 * t316 - t397 * t292 - t262 * t240 + (-pkin(9) * t330 + t217 * t320) * qJD(3) + (-t217 * t404 + (-pkin(9) * t387 + t147) * t317) * t390, -t133 * t220 - t164 * t187 + t196 * t235 + t403 * t442 + (-t40 + (pkin(9) * t187 + t133 * t311) * qJD(3)) * t320 + (pkin(9) * t161 - t292 * t74 + t419) * t316, -t133 * t221 - t164 * t189 - t196 * t236 - t402 * t442 + (t41 + (pkin(9) * t189 + t133 * t313) * qJD(3)) * t320 + (pkin(9) * t162 + t292 * t75 + t418) * t316, -t161 * t236 - t162 * t235 + t220 * t75 + t221 * t74 + (-t311 * t41 - t313 * t40) * t316 + t403 * t189 + t402 * t187 + (-t311 * t75 - t313 * t74) * t385, -t133 * t164 + t235 * t40 + t236 * t41 - t402 * t75 - t403 * t74 + (t100 * t316 + t133 * t385) * pkin(9), -t258 * t52 - t401 * t437, -t119 * t401 - t257 * t52 + t258 * t53 + t400 * t437, -t196 * t258 + t231 * t401 - t320 * t52 + t330 * t437, t119 * t330 - t196 * t257 - t231 * t400 + t320 * t53, -t231 * t330 - t417, t361 * t196 - t7 * t320 + t281 * t53 + t61 * t257 + (-t212 * t381 + (-qJD(5) * t199 - t443) * t315 + t436) * t231 - t20 * t330 + t394 * t119 + t400 * t106, t401 * t106 - t399 * t196 + t21 * t330 + t231 * t434 - t61 * t258 + t281 * t52 + t6 * t320 - t394 * t437, t12 * t173 - t340 * t428, -t12 * t172 + t173 * t323 + t340 * t427 - t428 * t57, -t12 * t320 + t173 * t196 + t222 * t428 + t330 * t340, -t172 * t196 - t222 * t427 - t320 * t323 + t330 * t57, -t222 * t330 - t417 (t114 * t318 - t116 * t314) * t196 - t2 * t320 - t202 * t323 + t25 * t172 + t422 * t57 + t427 * t54 + (t314 * t349 - t318 * t348) * t222 - t10 * t330 -(t114 * t314 + t116 * t318) * t196 + t1 * t320 + t202 * t12 + t25 * t173 - t422 * t340 + t428 * t54 + (t314 * t348 + t318 * t349) * t222 + t11 * t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240 * t442, t240 ^ 2 - t442 ^ 2, t195 + t416, -t196 - t414, t351, -t147 * t292 - t217 * t240 - t356, -t146 * t292 - t217 * t442 + t329, -t311 * t420 - pkin(3) * t161 - t418 - t147 * t187 - t240 * t74 - (t311 * t378 - t103) * t442, -t313 * t420 - pkin(3) * t162 + t419 - t147 * t189 + t240 * t75 - (t313 * t378 + t104) * t442, t103 * t189 + t104 * t187 + (-qJ(4) * t161 - qJD(4) * t187 + t442 * t74 + t41) * t313 + (qJ(4) * t162 + qJD(4) * t189 + t442 * t75 - t40) * t311, -pkin(3) * t100 - t103 * t74 - t104 * t75 - t133 * t147 + (-t311 * t74 + t313 * t75) * qJD(4) + (-t311 * t40 + t313 * t41) * qJ(4), t280 * t52 + t396 * t437, t119 * t396 - t279 * t52 - t280 * t53 + t395 * t437, t196 * t280 - t231 * t396 + t240 * t437, t119 * t240 - t196 * t279 - t231 * t395, -t231 * t240, t358 * t196 + t305 * t53 + t61 * t279 - t20 * t240 - t125 * t119 + (-t68 - t335 * t319 + (qJD(5) * t290 + t449) * t315) * t231 + t395 * t106, -t396 * t106 + t125 * t437 - t392 * t196 + t21 * t240 + t231 * t432 + t61 * t280 + t305 * t52, t12 * t198 - t340 * t424, -t12 * t197 + t198 * t323 + t340 * t423 - t424 * t57, t196 * t198 + t222 * t424 + t240 * t340, -t196 * t197 - t222 * t423 + t240 * t57, -t222 * t240 (t181 * t318 - t182 * t314) * t196 - t247 * t323 + t25 * t197 - t10 * t240 + t370 * t57 + t423 * t54 + (t314 * t345 - t318 * t344) * t222 -(t181 * t314 + t182 * t318) * t196 + t247 * t12 + t25 * t198 + t11 * t240 - t370 * t340 + t424 * t54 + (t314 * t344 + t318 * t345) * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189 * t442 + t161, t187 * t442 + t162, -t187 ^ 2 - t189 ^ 2, t187 * t75 + t189 * t74 + t100, 0, 0, 0, 0, 0, t53 - t448, t52 - t453, 0, 0, 0, 0, 0, -t323 - t452, t12 - t454; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t437 * t119, -t119 ^ 2 + t437 ^ 2, t52 + t453, -t53 - t448, t196, t106 * t437 + t21 * t231 + t7, t106 * t119 + t20 * t231 - t6, -t444, t441, t439, t438, t196 -(-t16 * t314 - t425) * t222 + (t196 * t318 - t222 * t380 + t437 * t57) * pkin(5) + t440, t445 + t15 + (-t17 * t222 - t4) * t314 + (t16 * t222 - t366) * t318 + (-t196 * t314 - t222 * t379 - t340 * t437) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444, t441, t439, t438, t196, t11 * t222 + t440, t10 * t222 - t1 + t445;];
tauc_reg  = t3;
