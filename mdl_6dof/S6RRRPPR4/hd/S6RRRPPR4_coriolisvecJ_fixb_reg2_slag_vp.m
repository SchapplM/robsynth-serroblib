% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:01
% EndTime: 2019-03-09 15:36:22
% DurationCPUTime: 9.50s
% Computational Cost: add. (11587->596), mult. (28375->771), div. (0->0), fcn. (19811->8), ass. (0->286)
t265 = cos(qJ(3));
t266 = cos(qJ(2));
t339 = qJD(1) * qJD(2);
t321 = t266 * t339;
t264 = sin(qJ(2));
t263 = sin(qJ(3));
t344 = qJD(3) * t263;
t328 = t264 * t344;
t338 = qJD(2) * qJD(3);
t163 = qJD(1) * t328 + (-t321 - t338) * t265;
t345 = qJD(2) * t266;
t330 = t263 * t345;
t343 = qJD(3) * t265;
t278 = t264 * t343 + t330;
t164 = qJD(1) * t278 + t263 * t338;
t261 = sin(pkin(10));
t384 = cos(pkin(10));
t105 = -t163 * t261 + t384 * t164;
t106 = -t163 * t384 - t261 * t164;
t262 = sin(qJ(6));
t409 = cos(qJ(6));
t351 = qJD(1) * t264;
t333 = t263 * t351;
t341 = t265 * qJD(2);
t211 = t333 - t341;
t331 = t265 * t351;
t347 = qJD(2) * t263;
t213 = t331 + t347;
t145 = t384 * t211 + t213 * t261;
t281 = -t261 * t211 + t213 * t384;
t88 = t262 * t145 + t281 * t409;
t26 = qJD(6) * t88 - t409 * t105 + t262 * t106;
t350 = qJD(1) * t266;
t241 = -qJD(3) + t350;
t231 = qJD(6) + t241;
t395 = t88 * t231;
t434 = t26 - t395;
t423 = -t145 * t409 + t262 * t281;
t405 = t423 ^ 2;
t406 = t88 ^ 2;
t433 = t405 - t406;
t404 = t88 * t423;
t317 = t384 * t263;
t205 = t261 * t265 + t317;
t194 = t205 * qJD(3);
t356 = t205 * t350 - t194;
t316 = t384 * t265;
t195 = qJD(3) * t316 - t261 * t344;
t332 = t263 * t350;
t355 = -t261 * t332 + t316 * t350 - t195;
t248 = t264 * t339;
t222 = -pkin(2) * t266 - pkin(8) * t264 - pkin(1);
t200 = t222 * qJD(1);
t254 = pkin(7) * t350;
t226 = qJD(2) * pkin(8) + t254;
t154 = t200 * t263 + t226 * t265;
t299 = pkin(2) * t264 - pkin(8) * t266;
t215 = t299 * qJD(2);
t201 = qJD(1) * t215;
t308 = pkin(7) * t248;
t96 = -qJD(3) * t154 + t265 * t201 + t263 * t308;
t54 = pkin(3) * t248 + t163 * qJ(4) - t213 * qJD(4) + t96;
t95 = t200 * t343 + t263 * t201 - t226 * t344 - t265 * t308;
t60 = -qJ(4) * t164 - qJD(4) * t211 + t95;
t20 = t261 * t54 + t384 * t60;
t228 = t241 * qJD(5);
t238 = qJ(5) * t248;
t15 = t238 - t228 + t20;
t10 = pkin(9) * t105 + t15;
t410 = pkin(4) + pkin(5);
t335 = t264 * t410;
t307 = qJD(2) * t335;
t400 = t261 * t60 - t384 * t54;
t11 = -t106 * pkin(9) - qJD(1) * t307 + t400;
t153 = t265 * t200 - t226 * t263;
t119 = -qJ(4) * t213 + t153;
t107 = -pkin(3) * t241 + t119;
t120 = -qJ(4) * t211 + t154;
t366 = t261 * t120;
t55 = t107 * t384 - t366;
t284 = qJD(5) - t55;
t421 = pkin(9) * t281;
t36 = t241 * t410 + t284 - t421;
t428 = pkin(9) * t145;
t318 = t384 * t120;
t56 = t261 * t107 + t318;
t49 = -t241 * qJ(5) + t56;
t40 = t49 + t428;
t6 = t262 * t36 + t40 * t409;
t2 = -qJD(6) * t6 - t262 * t10 + t409 * t11;
t253 = pkin(7) * t351;
t398 = qJD(2) * pkin(2);
t225 = t253 - t398;
t279 = -t211 * pkin(3) - qJD(4) - t225;
t271 = qJ(5) * t281 + t279;
t45 = -t145 * t410 + t271;
t432 = -t45 * t88 + t2;
t431 = t423 * t45;
t403 = -qJ(4) - pkin(8);
t319 = qJD(3) * t403;
t340 = t265 * qJD(4);
t190 = t263 * t319 + t340;
t276 = -t263 * qJD(4) + t265 * t319;
t129 = t384 * t190 + t261 * t276;
t214 = t299 * qJD(1);
t165 = pkin(7) * t333 + t265 * t214;
t361 = t265 * t266;
t285 = pkin(3) * t264 - qJ(4) * t361;
t133 = qJD(1) * t285 + t165;
t196 = t263 * t214;
t362 = t264 * t265;
t363 = t263 * t266;
t149 = t196 + (-pkin(7) * t362 - qJ(4) * t363) * qJD(1);
t76 = t261 * t133 + t384 * t149;
t70 = qJ(5) * t351 + t76;
t387 = t129 - t70;
t322 = qJD(6) * t409;
t342 = qJD(6) * t262;
t25 = -t262 * t105 - t409 * t106 - t145 * t322 + t281 * t342;
t397 = t231 * t423;
t430 = t25 - t397;
t62 = t384 * t119 - t366;
t357 = qJD(5) - t62;
t379 = t145 * t241;
t63 = t106 - t379;
t388 = (t133 - t276) * t384 + (-t149 + t190) * t261;
t427 = -t355 * pkin(9) - qJD(1) * t335 - t388;
t426 = -t356 * pkin(9) + t387;
t425 = t145 * t281;
t424 = -t421 + t357;
t302 = -t254 + (-t332 + t344) * pkin(3);
t383 = t281 ^ 2;
t422 = -0.2e1 * t339;
t1 = t409 * t10 + t262 * t11 + t36 * t322 - t342 * t40;
t5 = -t262 * t40 + t36 * t409;
t420 = -t231 * t5 + t1;
t65 = pkin(4) * t145 - t271;
t419 = t65 * t281;
t418 = t153 * t241 + t95;
t417 = t154 * t241 - t96;
t244 = pkin(7) * t361;
t174 = t263 * t222 + t244;
t414 = qJ(5) * t355 - t205 * qJD(5) + t302;
t143 = t164 * pkin(3) + pkin(7) * t321;
t274 = -t106 * qJ(5) - qJD(5) * t281 + t143;
t14 = -t105 * t410 - t274;
t314 = qJD(2) * t384;
t301 = t266 * t314;
t329 = t266 * t341;
t126 = -t195 * t264 - t261 * t329 - t263 * t301;
t184 = t205 * t264;
t346 = qJD(2) * t264;
t413 = -t266 * t105 + t126 * t241 + (qJD(1) * t184 + t145) * t346;
t365 = t261 * t263;
t204 = -t316 + t365;
t412 = t241 * t356 + (qJD(2) * t204 - t145) * t351;
t411 = -t205 * t105 - t106 * t204 + t145 * t355 + t281 * t356;
t408 = pkin(7) * t263;
t224 = t403 * t265;
t159 = -t224 * t261 - t403 * t317;
t131 = -pkin(9) * t205 + t159;
t160 = -t384 * t224 + t403 * t365;
t132 = pkin(9) * t204 + t160;
t68 = t262 * t131 + t132 * t409;
t402 = qJD(6) * t68 + t426 * t262 + t427 * t409;
t67 = t131 * t409 - t262 * t132;
t401 = -qJD(6) * t67 + t427 * t262 - t426 * t409;
t399 = t356 * t410 - t414;
t353 = t265 * t215 + t346 * t408;
t77 = -t264 * t340 + t285 * qJD(2) + (-t244 + (qJ(4) * t264 - t222) * t263) * qJD(3) + t353;
t354 = t263 * t215 + t222 * t343;
t92 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t362 + (-qJD(4) * t264 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t266) * t263 + t354;
t38 = t261 * t77 + t384 * t92;
t61 = t119 * t261 + t318;
t396 = t61 * t281;
t393 = -t356 * pkin(4) + t414;
t142 = t262 * t204 + t205 * t409;
t392 = -qJD(6) * t142 + t355 * t262 - t356 * t409;
t391 = -t204 * t322 + t205 * t342 + t356 * t262 + t355 * t409;
t247 = -pkin(3) * t384 - pkin(4);
t239 = -pkin(5) + t247;
t245 = pkin(3) * t261 + qJ(5);
t176 = t262 * t239 + t245 * t409;
t43 = t61 + t428;
t390 = qJD(6) * t176 + t424 * t262 + t409 * t43;
t389 = pkin(4) * t351 + t388;
t386 = t129 - t76;
t175 = t239 * t409 - t262 * t245;
t385 = -qJD(6) * t175 + t262 * t43 - t424 * t409;
t382 = t281 * t241;
t380 = t145 ^ 2;
t375 = t163 * t263;
t374 = t164 * t265;
t373 = t211 * t241;
t372 = t213 * t211;
t371 = t213 * t241;
t370 = t225 * t263;
t369 = t225 * t265;
t368 = t241 * t263;
t367 = t241 * t265;
t364 = t263 * t264;
t268 = qJD(1) ^ 2;
t360 = t266 * t268;
t267 = qJD(2) ^ 2;
t359 = t267 * t264;
t358 = t267 * t266;
t207 = t265 * t222;
t150 = -qJ(4) * t362 + t207 + (-pkin(3) - t408) * t266;
t157 = -qJ(4) * t364 + t174;
t94 = t261 * t150 + t384 * t157;
t216 = pkin(3) * t364 + t264 * pkin(7);
t259 = t264 ^ 2;
t352 = -t266 ^ 2 + t259;
t349 = qJD(2) * t159;
t348 = qJD(2) * t160;
t337 = pkin(7) * t363;
t255 = pkin(7) * t345;
t334 = t264 * t360;
t168 = t278 * pkin(3) + t255;
t252 = -t265 * pkin(3) - pkin(2);
t326 = t241 * t343;
t325 = t241 * t351;
t313 = pkin(1) * t422;
t310 = t211 + t341;
t309 = -t213 + t347;
t306 = -t62 * t241 - t20;
t305 = -t61 * t241 - t400;
t304 = t266 * t248;
t303 = -t160 * t105 + t106 * t159 - t129 * t145;
t89 = -qJ(5) * t266 + t94;
t185 = -t261 * t364 + t264 * t316;
t297 = t185 * qJ(5) - t216;
t296 = t261 * t92 - t384 * t77;
t295 = -t213 * pkin(3) - qJ(5) * t145;
t93 = t150 * t384 - t261 * t157;
t292 = t105 * t184 - t126 * t145;
t290 = -t380 - t383;
t289 = -t380 + t383;
t288 = -t153 * t265 - t154 * t263;
t287 = qJD(1) * t259 - t241 * t266;
t33 = qJ(5) * t346 - qJD(5) * t266 + t38;
t286 = t205 * qJ(5) - t252;
t91 = t266 * pkin(4) - t93;
t64 = t266 * pkin(5) - t185 * pkin(9) + t91;
t66 = pkin(9) * t184 + t89;
t27 = -t262 * t66 + t409 * t64;
t28 = t262 * t64 + t409 * t66;
t283 = t105 - t382;
t282 = t105 + t382;
t125 = t262 * t184 + t185 * t409;
t18 = -pkin(4) * t248 + t400;
t280 = t105 * t204 - t145 * t356;
t277 = t241 * t409 + t322;
t127 = t194 * t264 + t261 * t330 - t265 * t301;
t275 = -t127 * qJ(5) + t185 * qJD(5) - t168;
t272 = -t106 - t379;
t270 = t105 * t185 + t106 * t184 - t126 * t281 - t127 * t145;
t29 = pkin(4) * t105 + t274;
t173 = t207 - t337;
t172 = (-t241 - t350) * t346;
t166 = -pkin(7) * t331 + t196;
t141 = -t204 * t409 + t205 * t262;
t139 = pkin(4) * t204 - t286;
t124 = -t184 * t409 + t185 * t262;
t118 = -t174 * qJD(3) + t353;
t117 = (-t264 * t341 - t266 * t344) * pkin(7) + t354;
t109 = -t204 * t410 + t286;
t108 = pkin(4) * t184 - t297;
t90 = -t184 * t410 + t297;
t69 = pkin(4) * t281 - t295;
t53 = t355 * t241 + (qJD(2) * t205 - t281) * t351;
t48 = t241 * pkin(4) + t284;
t47 = -t281 * t410 + t295;
t46 = -pkin(4) * t126 - t275;
t42 = qJD(6) * t125 + t126 * t409 - t262 * t127;
t41 = t262 * t126 + t127 * t409 - t184 * t322 + t185 * t342;
t39 = t106 * t185 - t127 * t281;
t35 = -pkin(4) * t346 + t296;
t32 = -t106 * t266 + t127 * t241 + (qJD(1) * t185 + t281) * t346;
t31 = t126 * t410 + t275;
t30 = t106 * t205 - t281 * t355;
t24 = -pkin(9) * t126 + t33;
t23 = t127 * pkin(9) + t296 - t307;
t4 = -qJD(6) * t28 + t23 * t409 - t262 * t24;
t3 = qJD(6) * t27 + t262 * t23 + t24 * t409;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t304, t352 * t422, t358, -0.2e1 * t304, -t359, 0, -pkin(7) * t358 + t264 * t313, pkin(7) * t359 + t266 * t313, 0, 0, -t163 * t362 + (-t328 + t329) * t213 (-t211 * t265 - t213 * t263) * t345 + (t375 - t374 + (t211 * t263 - t213 * t265) * qJD(3)) * t264, t241 * t328 + t163 * t266 + (t213 * t264 + t265 * t287) * qJD(2), t164 * t364 + t211 * t278, t264 * t326 + t164 * t266 + (-t211 * t264 - t263 * t287) * qJD(2), t172, -t118 * t241 - t96 * t266 + (pkin(7) * t164 + t225 * t343) * t264 + ((pkin(7) * t211 + t370) * t266 + (t153 + (t173 + t337) * qJD(1)) * t264) * qJD(2), t117 * t241 + t95 * t266 + (-pkin(7) * t163 - t225 * t344) * t264 + ((pkin(7) * t213 + t369) * t266 + (-t154 + (-t174 + t244) * qJD(1)) * t264) * qJD(2), -t117 * t211 - t118 * t213 + t163 * t173 - t164 * t174 + t288 * t345 + (-t263 * t95 - t265 * t96 + (t153 * t263 - t154 * t265) * qJD(3)) * t264, t154 * t117 + t153 * t118 + t96 * t173 + t95 * t174 + (t225 + t253) * t255, t39, -t270, t32, t292, -t413, t172, t216 * t105 + t279 * t126 + t143 * t184 + t168 * t145 + t400 * t266 + t296 * t241 + (qJD(1) * t93 + t55) * t346, t216 * t106 + t279 * t127 + t143 * t185 + t168 * t281 + t20 * t266 + t38 * t241 + (-qJD(1) * t94 - t56) * t346, -t105 * t94 - t106 * t93 + t126 * t56 + t127 * t55 - t145 * t38 - t184 * t20 + t185 * t400 + t281 * t296, t143 * t216 - t168 * t279 + t20 * t94 - t296 * t55 + t38 * t56 - t400 * t93, t39, t32, t270, t172, t413, t292, t108 * t105 - t65 * t126 + t46 * t145 + t18 * t266 + t29 * t184 + t35 * t241 + (-qJD(1) * t91 - t48) * t346, -t105 * t89 + t106 * t91 + t126 * t49 - t127 * t48 - t145 * t33 - t15 * t184 + t18 * t185 + t281 * t35, -t108 * t106 + t65 * t127 - t46 * t281 - t15 * t266 - t29 * t185 - t33 * t241 + (qJD(1) * t89 + t49) * t346, t108 * t29 + t15 * t89 + t18 * t91 + t33 * t49 + t35 * t48 + t46 * t65, -t125 * t25 - t41 * t88, t124 * t25 - t125 * t26 + t41 * t423 - t42 * t88, -t41 * t231 - t25 * t266 + (-qJD(1) * t125 - t88) * t346, t124 * t26 + t42 * t423, -t42 * t231 - t26 * t266 + (qJD(1) * t124 + t423) * t346 (-t231 - t350) * t346, t14 * t124 + t2 * t266 + t4 * t231 + t90 * t26 + t31 * t423 + t45 * t42 + (-qJD(1) * t27 - t5) * t346, -t1 * t266 + t14 * t125 - t3 * t231 - t90 * t25 + t31 * t88 - t45 * t41 + (qJD(1) * t28 + t6) * t346, -t1 * t124 - t125 * t2 + t25 * t27 - t26 * t28 - t3 * t423 - t4 * t88 + t41 * t5 - t42 * t6, t1 * t28 + t14 * t90 + t2 * t27 + t3 * t6 + t31 * t45 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t334, t352 * t268, 0, t334, 0, 0, t268 * pkin(1) * t264, pkin(1) * t360, 0, 0, -t213 * t367 - t375 (-t163 + t373) * t265 + (-t164 + t371) * t263, -t326 + (t241 * t361 + t264 * t309) * qJD(1), -t211 * t368 - t374, t241 * t344 + (-t241 * t363 + t264 * t310) * qJD(1), t325, -pkin(2) * t164 + t165 * t241 + (pkin(8) * t367 + t370) * qJD(3) + ((-pkin(8) * t347 - t153) * t264 + (-pkin(7) * t310 - t370) * t266) * qJD(1), pkin(2) * t163 - t166 * t241 + (-pkin(8) * t368 + t369) * qJD(3) + ((-pkin(8) * t341 + t154) * t264 + (pkin(7) * t309 - t369) * t266) * qJD(1), t165 * t213 + t166 * t211 + ((qJD(3) * t213 - t164) * pkin(8) + t418) * t265 + ((qJD(3) * t211 - t163) * pkin(8) + t417) * t263, -t153 * t165 - t154 * t166 + (-t225 - t398) * t254 + (qJD(3) * t288 - t96 * t263 + t95 * t265) * pkin(8), t30, t411, t53, t280, -t412, t325, t105 * t252 + t143 * t204 + t388 * t241 + t356 * t279 + t302 * t145 + (-t55 - t349) * t351, t106 * t252 + t143 * t205 + t386 * t241 + t355 * t279 + t302 * t281 + (t56 - t348) * t351, t145 * t76 - t20 * t204 + t205 * t400 + t281 * t388 + t355 * t55 + t356 * t56 + t303, t143 * t252 + t159 * t400 + t20 * t160 - t279 * t302 + t386 * t56 - t388 * t55, t30, t53, -t411, t325, t412, t280, t139 * t105 + t29 * t204 - t356 * t65 + t389 * t241 + t393 * t145 + (t48 - t349) * t351, t145 * t70 - t15 * t204 + t18 * t205 + t281 * t389 - t355 * t48 + t356 * t49 + t303, -t139 * t106 - t29 * t205 + t355 * t65 - t387 * t241 - t393 * t281 + (-t49 + t348) * t351, t29 * t139 + t15 * t160 + t18 * t159 + t387 * t49 + t389 * t48 + t393 * t65, -t25 * t142 - t391 * t88, t25 * t141 - t142 * t26 + t391 * t423 + t392 * t88, -t391 * t231 + (-qJD(2) * t142 + t88) * t351, t26 * t141 - t392 * t423, t392 * t231 + (qJD(2) * t141 - t423) * t351, t231 * t351, t109 * t26 + t14 * t141 + t399 * t423 - t392 * t45 - t402 * t231 + (-qJD(2) * t67 + t5) * t351, -t109 * t25 + t14 * t142 + t399 * t88 - t391 * t45 + t401 * t231 + (qJD(2) * t68 - t6) * t351, -t1 * t141 - t142 * t2 + t25 * t67 - t26 * t68 + t391 * t5 + t392 * t6 + t401 * t423 + t402 * t88, t1 * t68 + t14 * t109 + t2 * t67 + t399 * t45 - t401 * t6 - t402 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, -t211 ^ 2 + t213 ^ 2, -t163 - t373, -t372, -t164 - t371, t248, -t225 * t213 - t417, t211 * t225 - t418, 0, 0, t425, t289, t63, -t425, -t282, t248, t279 * t281 + (-t145 * t213 + t314 * t351) * pkin(3) + t305, -t279 * t145 + (-t213 * t281 - t248 * t261) * pkin(3) + t306, t56 * t281 - t396 + (-t105 * t261 - t106 * t384) * pkin(3) + (t62 - t55) * t145, t55 * t61 - t56 * t62 + (t20 * t261 + t213 * t279 - t384 * t400) * pkin(3), t425, t63, -t289, t248, t282, -t425, -t419 - t69 * t145 + (pkin(4) - t247) * t248 + t305, -t245 * t105 + t247 * t106 + t281 * t49 - t396 + (-t357 + t48) * t145, -t145 * t65 + t245 * t248 + t281 * t69 - 0.2e1 * t228 + t238 - t306, t15 * t245 + t18 * t247 + t357 * t49 - t48 * t61 - t65 * t69, -t404, t433, t430, t404, t434, t248, -t175 * t248 - t231 * t390 - t423 * t47 - t432, t176 * t248 + t231 * t385 - t47 * t88 + t1 - t431, t175 * t25 - t176 * t26 + (t385 + t5) * t423 + (t390 - t6) * t88, t1 * t176 + t2 * t175 - t385 * t6 - t390 * t5 - t45 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, -t272, t290, t145 * t56 + t281 * t55 + t143, 0, 0, 0, 0, 0, 0, t283, t290, t272, t145 * t49 - t281 * t48 + t29, 0, 0, 0, 0, 0, 0, -t26 - t395, t25 + t397, t405 + t406, -t423 * t6 - t5 * t88 - t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248 + t425, t63, -t241 ^ 2 - t383, t241 * t49 + t18 + t419, 0, 0, 0, 0, 0, 0, -t231 ^ 2 * t262 - t248 * t409 - t281 * t423, -t231 * t277 + t248 * t262 - t281 * t88, t409 * t25 - t262 * t434 - t277 * t423, t2 * t409 + t262 * t420 + t277 * t6 - t45 * t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, -t433, -t430, -t404, -t434, -t248, t6 * t231 + t432, -t420 + t431, 0, 0;];
tauc_reg  = t7;
