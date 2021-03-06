% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRRR2
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:10
% EndTime: 2019-03-09 13:19:32
% DurationCPUTime: 10.79s
% Computational Cost: add. (24828->570), mult. (62918->746), div. (0->0), fcn. (49144->10), ass. (0->298)
t289 = cos(qJ(5));
t285 = sin(qJ(6));
t286 = sin(qJ(5));
t378 = t285 * t286;
t426 = cos(qJ(6));
t310 = t426 * t289 - t378;
t282 = qJD(2) + qJD(4);
t288 = sin(qJ(2));
t290 = cos(qJ(2));
t401 = cos(pkin(11));
t331 = qJD(1) * t401;
t320 = t290 * t331;
t400 = sin(pkin(11));
t330 = qJD(1) * t400;
t239 = -t288 * t330 + t320;
t305 = -t288 * t401 - t290 * t400;
t240 = t305 * qJD(1);
t287 = sin(qJ(4));
t427 = cos(qJ(4));
t312 = -t287 * t239 + t240 * t427;
t165 = -t289 * t282 - t286 * t312;
t168 = t282 * t286 - t289 * t312;
t311 = -t285 * t165 + t168 * t426;
t339 = t426 * qJD(6);
t317 = qJD(2) * t330;
t263 = t288 * t317;
t318 = qJD(2) * t331;
t227 = t290 * t318 - t263;
t341 = qJD(4) * t427;
t360 = qJD(4) * t287;
t365 = t288 * t318 + t290 * t317;
t324 = -t427 * t227 - t239 * t341 - t240 * t360 + t287 * t365;
t358 = qJD(5) * t289;
t359 = qJD(5) * t286;
t350 = t282 * t359 - t286 * t324 - t312 * t358;
t357 = qJD(6) * t285;
t86 = -t282 * t358 + t289 * t324 - t312 * t359;
t36 = t165 * t339 + t168 * t357 + t285 * t350 + t426 * t86;
t184 = t427 * t239 + t240 * t287;
t347 = t426 * t286;
t252 = t285 * t289 + t347;
t432 = qJD(5) + qJD(6);
t193 = t432 * t252;
t370 = t252 * t184 - t193;
t490 = -t310 * t36 + t311 * t370;
t452 = qJD(5) - t184;
t173 = qJD(6) + t452;
t295 = qJD(4) * t312 - t287 * t227 - t427 * t365;
t435 = t426 * qJD(5) + t339;
t369 = t310 * t184 - t289 * t435 + t378 * t432;
t489 = -t173 * t369 - t252 * t295;
t80 = t286 * t350;
t402 = -t165 * t358 - t80;
t480 = t452 * t286;
t437 = t168 * t480;
t455 = t184 * t289;
t486 = t165 * t455 - t289 * t86 + t402 - t437;
t107 = t426 * t165 + t168 * t285;
t37 = qJD(6) * t311 - t285 * t86 + t426 * t350;
t321 = t107 * t369 - t252 * t37;
t485 = t321 + t490;
t456 = t184 * t286;
t484 = pkin(10) * t456;
t278 = -t290 * pkin(2) - pkin(1);
t363 = qJD(1) * t278;
t257 = qJD(3) + t363;
t194 = -t239 * pkin(3) + t257;
t103 = -pkin(4) * t184 + pkin(9) * t312 + t194;
t421 = -qJ(3) - pkin(7);
t259 = t421 * t290;
t254 = qJD(1) * t259;
t242 = t400 * t254;
t258 = t421 * t288;
t253 = qJD(1) * t258;
t412 = qJD(2) * pkin(2);
t246 = t253 + t412;
t186 = t401 * t246 + t242;
t424 = t240 * pkin(8);
t154 = qJD(2) * pkin(3) + t186 + t424;
t332 = t401 * t254;
t187 = t400 * t246 - t332;
t425 = t239 * pkin(8);
t161 = t187 + t425;
t99 = t287 * t154 + t161 * t427;
t95 = t282 * pkin(9) + t99;
t54 = t289 * t103 - t286 * t95;
t483 = t452 * t54;
t275 = pkin(2) * t401 + pkin(3);
t345 = t400 * pkin(2);
t233 = t427 * t275 - t287 * t345;
t220 = t233 * qJD(4);
t191 = t401 * t253 + t242;
t169 = t191 + t424;
t190 = -t253 * t400 + t332;
t307 = t190 - t425;
t112 = t169 * t427 + t287 * t307;
t135 = -pkin(4) * t312 - pkin(9) * t184;
t362 = qJD(1) * t288;
t279 = pkin(2) * t362;
t206 = -pkin(3) * t240 + t279;
t113 = t135 + t206;
t58 = -t112 * t286 + t289 * t113;
t482 = -t286 * t220 - t58;
t59 = t289 * t112 + t286 * t113;
t481 = -t289 * t220 + t59;
t389 = t312 * t282;
t479 = t295 - t389;
t319 = t173 * t370 - t295 * t310;
t399 = t107 * t312;
t478 = t319 - t399;
t385 = t184 * t282;
t477 = -t324 - t385;
t388 = t184 ^ 2;
t390 = t312 ^ 2;
t476 = -t388 + t390;
t475 = (t359 - t456) * pkin(5);
t474 = pkin(5) * t312 + pkin(10) * t455;
t334 = qJD(2) * t421;
t235 = t290 * qJD(3) + t288 * t334;
t212 = t235 * qJD(1);
t236 = -t288 * qJD(3) + t290 * t334;
t213 = t236 * qJD(1);
t164 = -t212 * t400 + t401 * t213;
t144 = -t227 * pkin(8) + t164;
t167 = t401 * t212 + t400 * t213;
t145 = -pkin(8) * t365 + t167;
t294 = -t287 * t144 - t145 * t427 - t154 * t341 + t161 * t360;
t55 = t103 * t286 + t289 * t95;
t355 = qJD(1) * qJD(2);
t338 = t288 * t355;
t273 = pkin(2) * t338;
t195 = pkin(3) * t365 + t273;
t66 = -pkin(4) * t295 + pkin(9) * t324 + t195;
t13 = -qJD(5) * t55 + t286 * t294 + t289 * t66;
t438 = t452 * t55 + t13;
t46 = -pkin(10) * t165 + t55;
t352 = t426 * t46;
t45 = -pkin(10) * t168 + t54;
t39 = pkin(5) * t452 + t45;
t18 = t285 * t39 + t352;
t50 = qJD(4) * t99 - t427 * t144 + t287 * t145;
t35 = pkin(5) * t350 + t50;
t98 = t154 * t427 - t287 * t161;
t94 = -t282 * pkin(4) - t98;
t76 = t165 * pkin(5) + t94;
t473 = -t18 * t312 + t35 * t252 - t369 * t76;
t409 = t285 * t46;
t17 = t39 * t426 - t409;
t472 = t17 * t312 - t310 * t35 - t370 * t76;
t471 = -t107 * t370 - t310 * t37;
t83 = t86 * t286;
t469 = -t83 + (t358 - t455) * t168;
t468 = -t36 * t252 - t311 * t369;
t397 = t311 * t312;
t466 = t397 + t489;
t391 = t168 * t312;
t128 = t286 * t295;
t436 = -t358 * t452 + t128;
t465 = -t452 * t455 + t391 - t436;
t8 = -pkin(5) * t295 + t86 * pkin(10) + t13;
t12 = t103 * t358 + t286 * t66 - t289 * t294 - t359 * t95;
t9 = -pkin(10) * t350 + t12;
t299 = -t285 * t8 - t339 * t39 + t46 * t357 - t426 * t9;
t4 = -qJD(6) * t18 - t285 * t9 + t426 * t8;
t464 = t17 * t369 + t18 * t370 - t4 * t252 - t299 * t310;
t234 = t287 * t275 + t427 * t345;
t231 = pkin(9) + t234;
t422 = -pkin(10) - t231;
t333 = qJD(5) * t422;
t462 = t286 * t333 - t481 + t484;
t461 = -t289 * t333 - t474 - t482;
t428 = -pkin(10) - pkin(9);
t349 = qJD(5) * t428;
t65 = t286 * t135 + t289 * t98;
t460 = t286 * t349 + t484 - t65;
t64 = t289 * t135 - t286 * t98;
t459 = t289 * t349 + t474 - t64;
t398 = t107 * t311;
t394 = t165 * t312;
t458 = t173 * t312;
t457 = t452 * t312;
t382 = t184 * t312;
t451 = -t107 ^ 2 + t311 ^ 2;
t450 = t107 * t173 - t36;
t449 = t76 * t107 + t299;
t335 = -t50 * t289 + t94 * t359;
t446 = t54 * t312 + t335;
t414 = t50 * t286 + t94 * t358;
t445 = -t55 * t312 + t414;
t444 = -t194 * t184 + t294;
t443 = t194 * t312 - t50;
t441 = -0.2e1 * t355;
t439 = t12 - t483;
t196 = t401 * t258 + t259 * t400;
t176 = pkin(8) * t305 + t196;
t197 = t400 * t258 - t401 * t259;
t306 = t288 * t400 - t290 * t401;
t177 = -pkin(8) * t306 + t197;
t126 = t287 * t176 + t177 * t427;
t119 = t289 * t126;
t300 = t427 * t306;
t188 = -t287 * t305 + t300;
t304 = t287 * t306;
t189 = -t305 * t427 - t304;
t216 = pkin(3) * t306 + t278;
t120 = t188 * pkin(4) - t189 * pkin(9) + t216;
t70 = t286 * t120 + t119;
t366 = t234 * qJD(4) - t169 * t287 + t427 * t307;
t130 = t289 * t295;
t434 = t359 * t452 + t130;
t433 = qJD(2) * t306;
t431 = -t76 * t311 + t4;
t430 = t173 * t311 - t37;
t429 = t240 ^ 2;
t291 = qJD(2) ^ 2;
t423 = t289 * pkin(5);
t198 = t422 * t286;
t281 = t289 * pkin(10);
t199 = t231 * t289 + t281;
t148 = t285 * t198 + t199 * t426;
t420 = qJD(6) * t148 + t285 * t462 + t426 * t461;
t147 = t198 * t426 - t285 * t199;
t419 = -qJD(6) * t147 + t285 * t461 - t426 * t462;
t11 = t12 * t289;
t125 = -t427 * t176 + t177 * t287;
t411 = t125 * t50;
t261 = t428 * t286;
t262 = pkin(9) * t289 + t281;
t201 = t261 * t426 - t285 * t262;
t406 = qJD(6) * t201 + t285 * t459 + t426 * t460;
t202 = t285 * t261 + t262 * t426;
t405 = -qJD(6) * t202 - t285 * t460 + t426 * t459;
t403 = t366 + t475;
t96 = t295 * t188;
t301 = qJD(2) * t305;
t139 = t282 * t300 - t287 * t301 - t305 * t360;
t396 = t139 * t286;
t395 = t139 * t289;
t393 = t165 * t286;
t392 = t168 * t165;
t381 = t189 * t286;
t380 = t189 * t289;
t379 = t240 * t239;
t292 = qJD(1) ^ 2;
t375 = t290 * t292;
t374 = t291 * t288;
t373 = t291 * t290;
t175 = t401 * t235 + t400 * t236;
t367 = -t220 + t112;
t364 = t288 ^ 2 - t290 ^ 2;
t361 = qJD(2) * t240;
t280 = t288 * t412;
t351 = t288 * t375;
t343 = t189 * t359;
t149 = pkin(8) * t301 + t175;
t174 = -t235 * t400 + t401 * t236;
t293 = pkin(8) * t433 + t174;
t60 = t149 * t427 + t176 * t341 - t177 * t360 + t287 * t293;
t140 = -qJD(4) * t304 - t287 * t433 - t301 * t427 - t305 * t341;
t207 = -pkin(3) * t301 + t280;
t73 = t140 * pkin(4) + t139 * pkin(9) + t207;
t336 = -t286 * t60 + t289 * t73;
t329 = pkin(1) * t441;
t69 = t289 * t120 - t126 * t286;
t323 = t290 * t338;
t322 = t475 - t99;
t85 = t350 * t289;
t316 = t286 * t55 + t289 * t54;
t315 = t286 * t54 - t289 * t55;
t314 = -t184 * t94 + t231 * t295;
t230 = -pkin(4) - t233;
t313 = t452 * t456 - t434;
t52 = pkin(5) * t188 - pkin(10) * t380 + t69;
t56 = -pkin(10) * t381 + t70;
t27 = -t285 * t56 + t426 * t52;
t28 = t285 * t52 + t426 * t56;
t309 = t189 * t358 - t396;
t308 = -t343 - t395;
t19 = t120 * t358 - t126 * t359 + t286 * t73 + t289 * t60;
t297 = -qJD(5) * t316 - t13 * t286 + t11;
t61 = qJD(4) * t126 + t287 * t149 - t427 * t293;
t277 = -pkin(4) - t423;
t237 = t239 ^ 2;
t208 = t230 - t423;
t138 = t310 * t189;
t137 = t252 * t189;
t91 = pkin(5) * t381 + t125;
t41 = -t139 * t347 - t285 * t343 - t357 * t381 + (-t139 * t285 + t189 * t435) * t289;
t40 = t139 * t310 + t189 * t193;
t38 = pkin(5) * t309 + t61;
t22 = t426 * t45 - t409;
t21 = -t285 * t45 - t352;
t20 = -qJD(5) * t70 + t336;
t14 = -pkin(10) * t309 + t19;
t10 = pkin(10) * t395 + t140 * pkin(5) + (-t119 + (pkin(10) * t189 - t120) * t286) * qJD(5) + t336;
t6 = -qJD(6) * t28 + t10 * t426 - t285 * t14;
t5 = qJD(6) * t27 + t285 * t10 + t14 * t426;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t323, t364 * t441, t373, -0.2e1 * t323, -t374, 0, -pkin(7) * t373 + t288 * t329, pkin(7) * t374 + t290 * t329, 0, 0, -t227 * t305 + t240 * t433, -t227 * t306 - t239 * t433 - t240 * t301 + t305 * t365, -t306 * t291, t239 * t301 + t306 * t365, t305 * t291, 0, -t239 * t280 + t278 * t365 + t279 * t433 + (-t257 * t305 + t174) * qJD(2), -t175 * qJD(2) + t278 * t227 - t240 * t280 - t257 * t433 - t273 * t305, t164 * t305 - t167 * t306 + t174 * t240 + t175 * t239 + t186 * t433 + t187 * t301 - t196 * t227 - t197 * t365, t164 * t196 + t167 * t197 + t186 * t174 + t187 * t175 + (t257 + t363) * t280, t139 * t312 - t189 * t324, -t139 * t184 + t140 * t312 + t188 * t324 + t189 * t295, -t139 * t282, -t140 * t184 - t96, -t140 * t282, 0, t140 * t194 - t184 * t207 + t188 * t195 - t216 * t295 - t282 * t61, -t194 * t139 + t195 * t189 - t207 * t312 - t216 * t324 - t60 * t282, -t125 * t324 + t126 * t295 + t98 * t139 - t99 * t140 + t184 * t60 + t188 * t294 + t50 * t189 - t312 * t61, -t126 * t294 + t194 * t207 + t195 * t216 + t60 * t99 - t61 * t98 + t411, t168 * t308 - t380 * t86 (t165 * t289 + t168 * t286) * t139 + (-t85 + t83 + (-t168 * t289 + t393) * qJD(5)) * t189, -t130 * t189 + t168 * t140 - t86 * t188 + t308 * t452, t165 * t309 + t189 * t80, t128 * t189 - t165 * t140 - t188 * t350 - t309 * t452, t140 * t452 - t96, t125 * t350 + t13 * t188 + t54 * t140 + t61 * t165 + t189 * t414 + t20 * t452 - t295 * t69 - t396 * t94, -t12 * t188 - t125 * t86 - t55 * t140 + t61 * t168 - t189 * t335 - t19 * t452 + t295 * t70 - t94 * t395, -t19 * t165 - t70 * t350 - t20 * t168 + t69 * t86 + t316 * t139 + (qJD(5) * t315 - t12 * t286 - t13 * t289) * t189, t12 * t70 + t13 * t69 + t19 * t55 + t20 * t54 + t61 * t94 + t411, -t138 * t36 - t311 * t40, t107 * t40 + t137 * t36 - t138 * t37 - t311 * t41, -t138 * t295 + t140 * t311 - t173 * t40 - t188 * t36, t107 * t41 + t137 * t37, -t107 * t140 + t137 * t295 - t173 * t41 - t188 * t37, t140 * t173 - t96, t107 * t38 + t137 * t35 + t140 * t17 + t173 * t6 + t188 * t4 - t27 * t295 + t37 * t91 + t41 * t76, t138 * t35 - t140 * t18 - t173 * t5 + t188 * t299 + t28 * t295 + t311 * t38 - t36 * t91 - t40 * t76, -t107 * t5 + t137 * t299 - t138 * t4 + t17 * t40 - t18 * t41 + t27 * t36 - t28 * t37 - t311 * t6, t17 * t6 + t18 * t5 + t27 * t4 - t28 * t299 + t35 * t91 + t38 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t351, t364 * t292, 0, t351, 0, 0, t292 * pkin(1) * t288, pkin(1) * t375, 0, 0, t379, -t237 + t429, -t263 + (t320 - t239) * qJD(2), -t379, -t361 - t365, 0, -t190 * qJD(2) + t239 * t279 + t257 * t240 + t164, qJD(2) * t191 - t239 * t257 + t240 * t279 - t167 -(t187 + t190) * t240 + (-t191 + t186) * t239 + (-t227 * t401 - t365 * t400) * pkin(2), -t186 * t190 - t187 * t191 + (t164 * t401 + t167 * t400 - t257 * t362) * pkin(2), t382, t476, t477, -t382, t479, 0, t184 * t206 - t282 * t366 + t443, t206 * t312 + t282 * t367 + t444, t234 * t295 + t233 * t324 + (-t366 - t99) * t312 + (-t367 + t98) * t184, -t194 * t206 - t50 * t233 - t234 * t294 - t366 * t98 - t367 * t99, t469, t486, t465, t165 * t480 - t85, t313 - t394, t457, t230 * t350 + t314 * t286 + t366 * t165 + (-t231 * t358 + t482) * t452 + t446, -t230 * t86 + t314 * t289 + t366 * t168 + (t231 * t359 + t481) * t452 + t445, t59 * t165 + t58 * t168 + t11 + (-t220 * t165 - t231 * t350 + t54 * t184 + (t168 * t231 - t54) * qJD(5)) * t289 + (t220 * t168 + t55 * t184 - t231 * t86 - t13 + (t165 * t231 - t55) * qJD(5)) * t286, -t220 * t315 + t50 * t230 + t231 * t297 + t366 * t94 - t54 * t58 - t55 * t59, t468, t485, t466, t471, t478, t458, t107 * t403 - t147 * t295 - t173 * t420 + t208 * t37 + t472, t148 * t295 + t173 * t419 - t208 * t36 + t311 * t403 + t473, t107 * t419 + t147 * t36 - t148 * t37 + t311 * t420 + t464, t4 * t147 - t148 * t299 - t17 * t420 - t18 * t419 + t35 * t208 + t403 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361 + t365, -t263 + (t320 + t239) * qJD(2), -t237 - t429, -t186 * t240 - t187 * t239 + t273, 0, 0, 0, 0, 0, 0, -t295 - t389, -t324 + t385, -t388 - t390, -t99 * t184 - t312 * t98 + t195, 0, 0, 0, 0, 0, 0, t313 + t394, -t289 * t452 ^ 2 + t128 + t391 (t165 * t184 + t86) * t289 + t437 + t402, t286 * t439 + t289 * t438 + t312 * t94, 0, 0, 0, 0, 0, 0, t319 + t399, t397 - t489, t321 - t490, t17 * t370 - t18 * t369 - t252 * t299 + t310 * t4 + t312 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t382, t476, t477, -t382, t479, 0, t99 * t282 + t443, t98 * t282 + t444, 0, 0, t469, t486, t465, t393 * t452 - t85, -t452 * t480 - t130 - t394, t457, -pkin(4) * t350 + pkin(9) * t436 - t99 * t165 - t452 * t64 - t456 * t94 + t446, pkin(4) * t86 + pkin(9) * t434 - t99 * t168 + t452 * t65 - t455 * t94 + t445, t65 * t165 + t64 * t168 + t11 + (-t483 + (qJD(5) * t168 - t350) * pkin(9)) * t289 + ((qJD(5) * t165 - t86) * pkin(9) - t438) * t286, -t50 * pkin(4) + pkin(9) * t297 - t54 * t64 - t55 * t65 - t94 * t99, t468, t485, t466, t471, t478, t458, t107 * t322 + t173 * t405 - t201 * t295 + t277 * t37 + t472, -t173 * t406 + t202 * t295 - t277 * t36 + t311 * t322 + t473, -t107 * t406 + t201 * t36 - t202 * t37 - t311 * t405 + t464, t17 * t405 + t18 * t406 + t4 * t201 - t202 * t299 + t35 * t277 + t322 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, -t165 ^ 2 + t168 ^ 2, t165 * t452 - t86, -t392, t168 * t452 - t350, -t295, -t94 * t168 + t438, t165 * t94 - t439, 0, 0, t398, t451, t450, -t398, t430, -t295, -t21 * t173 + (-t107 * t168 - t173 * t357 - t295 * t426) * pkin(5) + t431, t22 * t173 + (-t168 * t311 - t173 * t339 + t285 * t295) * pkin(5) + t449, t18 * t311 + t22 * t107 - t17 * t107 + t21 * t311 + (t426 * t36 - t285 * t37 + (-t107 * t426 + t285 * t311) * qJD(6)) * pkin(5), -t17 * t21 - t18 * t22 + (t426 * t4 - t168 * t76 - t285 * t299 + (-t17 * t285 + t18 * t426) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t398, t451, t450, -t398, t430, -t295, t18 * t173 + t431, t17 * t173 + t449, 0, 0;];
tauc_reg  = t1;
