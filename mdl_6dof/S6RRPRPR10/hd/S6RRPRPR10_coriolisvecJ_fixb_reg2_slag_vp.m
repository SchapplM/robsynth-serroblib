% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:25
% EndTime: 2019-03-09 11:09:49
% DurationCPUTime: 10.14s
% Computational Cost: add. (16318->604), mult. (43318->828), div. (0->0), fcn. (34963->10), ass. (0->294)
t240 = sin(pkin(11));
t242 = cos(pkin(11));
t245 = sin(qJ(4));
t406 = cos(qJ(4));
t336 = t406 * t242;
t241 = sin(pkin(6));
t248 = cos(qJ(2));
t367 = t241 * t248;
t297 = t336 * t367;
t281 = qJD(1) * t297;
t355 = qJD(1) * t248;
t334 = t241 * t355;
t311 = t240 * t334;
t330 = qJD(4) * t406;
t350 = qJD(4) * t245;
t361 = t240 * t350 - t242 * t330 - t245 * t311 + t281;
t205 = t240 * t406 + t245 * t242;
t263 = t205 * t367;
t165 = qJD(1) * t263;
t198 = t205 * qJD(4);
t362 = t165 - t198;
t243 = cos(pkin(6));
t356 = qJD(1) * t243;
t230 = qJD(2) + t356;
t246 = sin(qJ(2));
t357 = qJD(1) * t241;
t335 = t246 * t357;
t178 = t230 * t242 - t240 * t335;
t179 = t230 * t240 + t242 * t335;
t127 = -t406 * t178 + t179 * t245;
t222 = -qJD(4) + t334;
t428 = t127 * t222;
t344 = qJD(1) * qJD(2);
t328 = t248 * t344;
t308 = t241 * t328;
t292 = t240 * t308;
t82 = t245 * (qJD(4) * t179 + t292) - qJD(2) * t281 - t178 * t330;
t438 = t82 - t428;
t437 = t82 + t428;
t244 = sin(qJ(6));
t247 = cos(qJ(6));
t103 = t127 * t244 - t222 * t247;
t280 = t245 * t178 + t179 * t406;
t415 = qJD(6) + t280;
t427 = t247 * t415;
t436 = t103 * t427;
t302 = pkin(2) * t246 - qJ(3) * t248;
t191 = t302 * t357;
t343 = pkin(1) * t356;
t192 = -pkin(8) * t335 + t248 * t343;
t137 = t242 * t191 - t192 * t240;
t365 = t242 * t248;
t273 = t241 * (pkin(3) * t246 - pkin(9) * t365);
t109 = qJD(1) * t273 + t137;
t138 = t240 * t191 + t242 * t192;
t119 = -pkin(9) * t311 + t138;
t402 = pkin(9) + qJ(3);
t217 = t402 * t240;
t218 = t402 * t242;
t164 = -t245 * t217 + t218 * t406;
t384 = qJD(3) * t205 + qJD(4) * t164 + t109 * t406 - t245 * t119;
t193 = pkin(8) * t334 + t246 * t343;
t159 = pkin(3) * t311 + t193;
t435 = t361 * qJ(5) - qJD(5) * t205 - t159;
t268 = t244 * t82 - t415 * t427;
t434 = t127 ^ 2;
t407 = pkin(4) + pkin(10);
t433 = t362 * t407 - t435;
t368 = t241 * t246;
t313 = t407 * t368;
t432 = pkin(5) * t361 - qJD(1) * t313 - t384;
t386 = -qJD(3) * t336 + t406 * t119 + t217 * t330 + (qJD(3) * t240 + qJD(4) * t218 + t109) * t245;
t101 = -t247 * t127 - t222 * t244;
t320 = t415 * t101;
t352 = qJD(2) * t246;
t333 = t241 * t352;
t223 = qJD(1) * t333;
t348 = qJD(6) * t247;
t349 = qJD(6) * t244;
t261 = qJD(2) * t263;
t416 = qJD(4) * t280;
t83 = qJD(1) * t261 + t416;
t45 = -t127 * t348 - t222 * t349 - t247 * t223 - t244 * t83;
t431 = t45 - t320;
t430 = t101 * t127;
t429 = t103 * t127;
t167 = qJ(3) * t230 + t193;
t188 = (-pkin(2) * t248 - qJ(3) * t246 - pkin(1)) * t241;
t173 = qJD(1) * t188;
t112 = -t167 * t240 + t173 * t242;
t259 = -pkin(3) * t334 - pkin(9) * t179 + t112;
t70 = t406 * t259;
t113 = t242 * t167 + t240 * t173;
t84 = pkin(9) * t178 + t113;
t41 = t245 * t84 - t70;
t364 = -qJD(5) - t41;
t342 = pkin(1) * qJD(2) * t243;
t312 = qJD(1) * t342;
t186 = pkin(8) * t308 + t246 * t312;
t150 = pkin(3) * t292 + t186;
t267 = qJ(5) * t82 - qJD(5) * t280 + t150;
t20 = t407 * t83 + t267;
t363 = pkin(5) * t280 - t364;
t28 = t222 * t407 + t363;
t161 = -pkin(2) * t230 + qJD(3) - t192;
t124 = -pkin(3) * t178 + t161;
t255 = -qJ(5) * t280 + t124;
t36 = t127 * t407 + t255;
t6 = t244 * t28 + t247 * t36;
t296 = qJD(2) * t313;
t256 = t245 * t259;
t170 = (qJD(2) * t302 - qJD(3) * t246) * t241;
t154 = qJD(1) * t170;
t185 = -pkin(8) * t223 + t248 * t312;
t155 = qJD(3) * t230 + t185;
t105 = t242 * t154 - t155 * t240;
t266 = qJD(2) * t273;
t73 = qJD(1) * t266 + t105;
t106 = t240 * t154 + t242 * t155;
t85 = -pkin(9) * t292 + t106;
t326 = qJD(4) * t256 + t245 * t85 + t84 * t330 - t406 * t73;
t9 = -pkin(5) * t82 - qJD(1) * t296 + t326;
t2 = -qJD(6) * t6 - t20 * t244 + t247 * t9;
t426 = t415 * t6 + t2;
t299 = t244 * t36 - t247 * t28;
t1 = -qJD(6) * t299 + t20 * t247 + t244 * t9;
t425 = t299 * t415 + t1;
t338 = t242 * t368;
t196 = t240 * t243 + t338;
t339 = t240 * t368;
t366 = t242 * t243;
t282 = t339 - t366;
t264 = t406 * t282;
t144 = t196 * t245 + t264;
t145 = t196 * t406 - t245 * t282;
t97 = qJD(4) * t145 + t261;
t423 = -t222 * t97 + t241 * ((qJD(1) * t144 + t127) * t352 - t248 * t83);
t351 = qJD(2) * t248;
t332 = t241 * t351;
t310 = t240 * t332;
t96 = qJD(4) * t264 - qJD(2) * t297 + (qJD(4) * t196 + t310) * t245;
t422 = t222 * t96 + t241 * ((qJD(1) * t145 + t280) * t352 + t248 * t82);
t409 = t280 ^ 2;
t237 = t241 ^ 2;
t421 = -0.2e1 * t237 * t344;
t387 = qJ(5) * t335 + t386;
t61 = t82 * t145;
t420 = -t280 * t96 - t61;
t67 = t82 * t205;
t419 = -t361 * t280 - t67;
t417 = t280 * t222;
t42 = t406 * t84 + t256;
t40 = t222 * qJ(5) - t42;
t403 = t127 * pkin(5);
t30 = -t40 - t403;
t414 = t30 * t415 + t407 * t82;
t194 = -pkin(8) * t333 + t248 * t342;
t177 = qJD(3) * t243 + t194;
t117 = t240 * t170 + t242 * t177;
t104 = -pkin(9) * t310 + t117;
t405 = pkin(1) * t246;
t200 = pkin(8) * t367 + t243 * t405;
t187 = qJ(3) * t243 + t200;
t134 = t242 * t187 + t240 * t188;
t108 = -pkin(9) * t282 + t134;
t116 = t242 * t170 - t177 * t240;
t89 = t116 + t266;
t133 = -t187 * t240 + t242 * t188;
t91 = -pkin(3) * t367 - pkin(9) * t196 + t133;
t279 = -t406 * t104 + t108 * t350 - t245 * t89 - t91 * t330;
t21 = -t241 * (qJ(5) * t352 - qJD(5) * t248) + t279;
t411 = t222 * t361 + (qJD(2) * t205 - t280) * t335;
t204 = t240 * t245 - t336;
t410 = t222 * t362 + (qJD(2) * t204 - t127) * t335;
t408 = t242 ^ 2;
t404 = pkin(1) * t248;
t236 = -t242 * pkin(3) - pkin(2);
t288 = -qJ(5) * t205 + t236;
t132 = t204 * t407 + t288;
t163 = t406 * t217 + t218 * t245;
t140 = pkin(5) * t205 + t163;
t66 = t132 * t247 + t140 * t244;
t401 = qJD(6) * t66 - t244 * t433 + t247 * t432;
t65 = -t132 * t244 + t140 * t247;
t400 = -qJD(6) * t65 + t244 * t432 + t247 * t433;
t50 = pkin(4) * t127 + t255;
t399 = t280 * t50;
t274 = t223 * t244 - t247 * t83;
t46 = qJD(6) * t103 + t274;
t395 = t244 * t46;
t393 = t247 * t45;
t78 = t247 * t82;
t54 = t406 * t108 + t245 * t91;
t389 = pkin(5) * t362 - t387;
t388 = -pkin(4) * t362 + t435;
t385 = pkin(4) * t335 + t384;
t383 = qJ(5) * t127;
t382 = t103 * t101;
t381 = t112 * t246;
t380 = t113 * t246;
t378 = t127 * t280;
t374 = t186 * t246;
t373 = t204 * t244;
t372 = t204 * t247;
t239 = t248 ^ 2;
t371 = t237 * t239;
t370 = t237 * qJD(1) ^ 2;
t369 = t240 * t248;
t358 = t246 ^ 2 - t239;
t354 = qJD(2) * t163;
t353 = qJD(2) * t164;
t347 = qJD(6) * t407;
t341 = t239 * t370;
t340 = t248 * t370;
t337 = t247 * t367;
t331 = t409 - t434;
t327 = -qJD(4) * t70 - t245 * t73 + t84 * t350 - t406 * t85;
t142 = -t247 * t165 + t244 * t335;
t324 = t198 * t247 + t142;
t143 = t165 * t244 + t247 * t335;
t323 = -t198 * t244 + t143;
t321 = t415 * t244;
t317 = t230 + t356;
t314 = qJD(2) * t242 - t178;
t309 = t222 * t335;
t53 = -t245 * t108 + t406 * t91;
t307 = -qJD(5) * t222 - t327;
t306 = pkin(1) * t421;
t305 = t244 * t299 + t247 * t6;
t301 = t127 * t97 + t144 * t83;
t300 = -t163 * t82 - t164 * t83;
t52 = pkin(4) * t367 - t53;
t37 = t145 * pkin(5) + pkin(10) * t367 + t52;
t231 = pkin(8) * t368;
t146 = pkin(3) * t339 + t231 + (t236 - t404) * t243;
t253 = -t145 * qJ(5) + t146;
t48 = t144 * t407 + t253;
t15 = -t244 * t48 + t247 * t37;
t16 = t244 * t37 + t247 * t48;
t298 = pkin(4) * t223;
t295 = -t409 - t434;
t294 = -t242 * t178 + t179 * t240;
t293 = qJ(5) * t223;
t291 = t237 * t246 * t328;
t51 = qJ(5) * t367 - t54;
t285 = -t222 * t42 - t326;
t284 = -t321 * t415 - t78;
t283 = t245 * t104 + t108 * t330 + t91 * t350 - t406 * t89;
t121 = t144 * t247 + t244 * t367;
t278 = t314 * t357;
t277 = (qJD(2) * t240 - t179) * t357;
t275 = -t127 * t362 + t204 * t83;
t272 = t204 * t349 - t324;
t271 = t204 * t348 - t323;
t270 = t248 * t277;
t269 = t241 * t282;
t195 = t200 * qJD(2);
t17 = -t298 + t326;
t262 = t127 * t96 + t144 * t82 - t145 * t83 - t280 * t97;
t160 = pkin(3) * t310 + t195;
t29 = pkin(4) * t83 + t267;
t14 = -t293 - t307;
t260 = -qJ(3) * t352 + (-pkin(2) * qJD(2) + qJD(3) - t161) * t248;
t257 = t127 * t361 + t204 * t82 - t205 * t83 + t280 * t362;
t254 = qJ(5) * t96 - qJD(5) * t145 + t160;
t251 = -t205 * t308 - t416;
t215 = t246 * t340;
t209 = -0.2e1 * t291;
t199 = t243 * t404 - t231;
t190 = t231 + (-pkin(2) - t404) * t243;
t156 = (-t222 * t241 - t237 * t355) * t352;
t147 = pkin(4) * t204 + t288;
t141 = -t204 * pkin(5) + t164;
t122 = -t244 * t144 + t337;
t64 = pkin(4) * t280 + t383;
t62 = t144 * pkin(4) + t253;
t56 = qJD(6) * t121 + t244 * t97 + t247 * t333;
t55 = -qJD(6) * t337 - t247 * t97 + (qJD(6) * t144 + t333) * t244;
t49 = t280 * t407 + t383;
t43 = t247 * t46;
t39 = pkin(4) * t222 - t364;
t38 = -pkin(5) * t144 - t51;
t35 = pkin(4) * t97 + t254;
t34 = t42 - t403;
t27 = t407 * t97 + t254;
t22 = -pkin(4) * t333 + t283;
t13 = -pkin(5) * t97 - t21;
t12 = -t96 * pkin(5) + t283 - t296;
t11 = t244 * t34 + t247 * t49;
t10 = -t244 * t49 + t247 * t34;
t8 = -pkin(5) * t83 - t14;
t4 = -qJD(6) * t16 + t12 * t247 - t244 * t27;
t3 = qJD(6) * t15 + t12 * t244 + t247 * t27;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t291, t358 * t421, t317 * t332, t209, -t317 * t333, 0, -t186 * t243 - t195 * t230 + t246 * t306, -t185 * t243 - t194 * t230 + t248 * t306 (t185 * t248 + t374 + (-t192 * t248 - t193 * t246) * qJD(2) + (t194 * t248 + t195 * t246 + (-t199 * t248 - t200 * t246) * qJD(2)) * qJD(1)) * t241, t185 * t200 - t186 * t199 - t192 * t195 + t193 * t194 (qJD(1) * t196 + t179) * t242 * t332 ((t408 * t243 + (-t196 - t338) * t240) * qJD(1) - t294) * t332 (t179 * t368 + (t196 * t368 - 0.2e1 * t242 * t371) * qJD(1)) * qJD(2) (qJD(1) * t269 - t178 * t241) * t240 * t351 (t178 * t368 + (0.2e1 * t240 * t371 - t246 * t269) * qJD(1)) * qJD(2), t209, -t186 * t366 - t195 * t178 + (t240 * t374 + (-qJD(1) * t116 - t105) * t248 + (t161 * t369 + t381 + (t133 * t246 + t190 * t369) * qJD(1)) * qJD(2)) * t241, t179 * t195 + t186 * t196 + ((qJD(1) * t117 + t106) * t248 + (t161 * t365 - t380 + (-t134 * t246 + t190 * t365) * qJD(1)) * qJD(2)) * t241, t117 * t178 - t106 * t282 - t116 * t179 - t105 * t196 + (-t112 * t242 - t113 * t240 + (-t133 * t242 - t134 * t240) * qJD(1)) * t332, t105 * t133 + t106 * t134 + t112 * t116 + t113 * t117 + t161 * t195 + t186 * t190, t420, t262, t422, t301, -t423, t156, t124 * t97 + t127 * t160 + t144 * t150 + t146 * t83 + t222 * t283 + (t326 * t248 + (qJD(1) * t53 - t41) * t352) * t241, -t124 * t96 + t280 * t160 + t145 * t150 - t146 * t82 - t222 * t279 + (-t327 * t248 + (-qJD(1) * t54 - t42) * t352) * t241, t127 * t279 + t144 * t327 + t145 * t326 + t280 * t283 - t41 * t96 - t42 * t97 + t53 * t82 - t54 * t83, t124 * t160 + t146 * t150 - t279 * t42 + t283 * t41 - t326 * t53 - t327 * t54, t156, -t422, t423, t420, t262, t301, t127 * t21 + t14 * t144 + t145 * t17 + t22 * t280 - t39 * t96 + t40 * t97 + t51 * t83 - t52 * t82, -t127 * t35 - t144 * t29 - t22 * t222 - t50 * t97 - t62 * t83 + (-t17 * t248 + (qJD(1) * t52 + t39) * t352) * t241, -t280 * t35 - t145 * t29 + t21 * t222 + t50 * t96 + t62 * t82 + (t14 * t248 + (-qJD(1) * t51 - t40) * t352) * t241, t14 * t51 + t17 * t52 + t21 * t40 + t22 * t39 + t29 * t62 + t35 * t50, t103 * t56 + t122 * t45, -t101 * t56 - t103 * t55 - t121 * t45 + t122 * t46, -t103 * t96 + t122 * t82 - t145 * t45 + t415 * t56, t101 * t55 - t121 * t46, t101 * t96 - t121 * t82 - t145 * t46 - t415 * t55, -t415 * t96 - t61, t101 * t13 - t121 * t8 + t145 * t2 - t15 * t82 + t299 * t96 + t30 * t55 + t38 * t46 + t4 * t415, -t1 * t145 + t103 * t13 - t122 * t8 + t16 * t82 - t3 * t415 + t30 * t56 - t38 * t45 + t6 * t96, t1 * t121 - t101 * t3 - t103 * t4 + t122 * t2 + t15 * t45 - t16 * t46 + t299 * t56 - t55 * t6, t1 * t16 + t13 * t30 + t15 * t2 - t299 * t4 + t3 * t6 + t38 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t358 * t370 (qJD(2) - t230) * t334, t215, t230 * t335 - t223, 0, t193 * t230 + t370 * t405 - t186, pkin(1) * t340 + t192 * t230 - t185, 0, 0, t242 * t270 ((-t240 ^ 2 + t408) * qJD(2) + t294) * t334, t242 * t341 + t246 * t277, -t314 * t311, -t240 * t341 + t246 * t278, t215, t178 * t193 - t186 * t242 + (t137 * t248 + t240 * t260 - t381) * t357, -t179 * t193 + t186 * t240 + (-t138 * t248 + t242 * t260 + t380) * t357, t137 * t179 - t138 * t178 + (qJD(3) * t178 + t112 * t334 + t106) * t242 + (qJD(3) * t179 + t113 * t334 - t105) * t240, -pkin(2) * t186 - t112 * t137 - t113 * t138 - t161 * t193 + (-t112 * t240 + t113 * t242) * qJD(3) + (-t105 * t240 + t106 * t242) * qJ(3), t419, t257, t411, t275, -t410, t309, -t127 * t159 + t150 * t204 + t236 * t83 + t384 * t222 - t362 * t124 + (t41 - t354) * t335, -t280 * t159 + t150 * t205 - t236 * t82 - t386 * t222 - t361 * t124 + (t42 - t353) * t335, t127 * t386 + t204 * t327 + t205 * t326 + t280 * t384 - t361 * t41 + t362 * t42 + t300, -t124 * t159 + t150 * t236 + t163 * t326 - t164 * t327 + t384 * t41 - t386 * t42, t309, -t411, t410, t419, t257, t275, t127 * t387 + t14 * t204 + t17 * t205 + t280 * t385 - t361 * t39 - t362 * t40 + t300, -t147 * t83 - t204 * t29 + t362 * t50 - t385 * t222 - t388 * t127 + (-t39 + t354) * t335, t147 * t82 - t205 * t29 + t361 * t50 + t387 * t222 - t388 * t280 + (t40 + t353) * t335, -t14 * t164 + t147 * t29 + t163 * t17 + t385 * t39 + t387 * t40 + t388 * t50, t103 * t271 - t373 * t45, t324 * t103 + t323 * t101 + (-t395 - t393 + (-t101 * t247 - t103 * t244) * qJD(6)) * t204, -t103 * t361 - t205 * t45 + t271 * t415 - t373 * t82, t101 * t272 - t372 * t46, t101 * t361 - t205 * t46 - t272 * t415 - t372 * t82, -t361 * t415 - t67, t101 * t389 + t141 * t46 + t2 * t205 + t272 * t30 + t299 * t361 - t372 * t8 - t401 * t415 - t65 * t82, -t1 * t205 + t103 * t389 - t141 * t45 + t271 * t30 + t361 * t6 + t373 * t8 + t400 * t415 + t66 * t82, t142 * t6 - t143 * t299 + t45 * t65 - t46 * t66 + t305 * t198 + t401 * t103 + t400 * t101 + (t1 * t247 - t2 * t244 + (-t244 * t6 + t247 * t299) * qJD(6)) * t204, t1 * t66 + t141 * t8 + t2 * t65 + t299 * t401 + t30 * t389 - t400 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, t248 * t278, -t178 ^ 2 - t179 ^ 2, t112 * t179 - t113 * t178 + t186, 0, 0, 0, 0, 0, 0, t83 - t417, -t438, t295, t127 * t42 - t280 * t41 + t150, 0, 0, 0, 0, 0, 0, t295, t251 + t417, t438, -t127 * t40 - t280 * t39 + t29, 0, 0, 0, 0, 0, 0, t268 + t430, t244 * t415 ^ 2 + t429 + t78, -t244 * t431 - t43 + t436, t127 * t30 - t244 * t426 + t247 * t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378, t331, -t437, -t378, t251 - t417, t223, -t124 * t280 + t285, t124 * t127 + t222 * t41 + t327, 0, 0, t223, t437, t83 + t417, t378, t331, -t378, pkin(4) * t82 - qJ(5) * t83 + (-t40 - t42) * t280 + (t39 + t364) * t127, t127 * t64 - t285 - 0.2e1 * t298 + t399, -t127 * t50 + t222 * t364 + t280 * t64 + 0.2e1 * t293 + t307, -pkin(4) * t17 - qJ(5) * t14 + t364 * t40 - t39 * t42 - t50 * t64, -t103 * t321 - t393, -t43 - t436 + (t45 + t320) * t244, t284 + t429, t247 * t320 + t395, t268 - t430, t415 * t127, qJ(5) * t46 - t127 * t299 + t244 * t8 + (t244 * t347 - t10) * t415 + t363 * t101 + t414 * t247, -qJ(5) * t45 - t127 * t6 + t247 * t8 + (t247 * t347 + t11) * t415 + t363 * t103 - t414 * t244, t10 * t103 + t101 * t11 + (-t280 * t6 - t407 * t45 - t2 + (t101 * t407 - t6) * qJD(6)) * t247 + (-t280 * t299 + t407 * t46 - t1 + (-t103 * t407 - t299) * qJD(6)) * t244, qJ(5) * t8 + t10 * t299 - t11 * t6 + t363 * t30 - (qJD(6) * t305 + t1 * t244 + t2 * t247) * t407; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t437, t223 - t378, -t222 ^ 2 - t409, -t222 * t40 + t17 + t399, 0, 0, 0, 0, 0, 0, t101 * t222 + t284, t103 * t222 + t268, t431 * t247 + (t103 * t415 - t46) * t244, t222 * t30 + t244 * t425 + t247 * t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t382, -t101 ^ 2 + t103 ^ 2, -t431, -t382, -t274 + (-qJD(6) + t415) * t103, -t82, -t103 * t30 + t426, t101 * t30 - t425, 0, 0;];
tauc_reg  = t5;
