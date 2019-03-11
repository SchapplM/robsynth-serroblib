% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPPR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:05
% EndTime: 2019-03-09 15:36:23
% DurationCPUTime: 8.31s
% Computational Cost: add. (7496->591), mult. (16824->760), div. (0->0), fcn. (12047->12), ass. (0->290)
t263 = sin(qJ(6));
t267 = cos(qJ(6));
t264 = sin(qJ(3));
t268 = cos(qJ(3));
t269 = cos(qJ(2));
t353 = qJD(1) * qJD(2);
t332 = t269 * t353;
t265 = sin(qJ(2));
t352 = t265 * qJDD(1);
t355 = t268 * qJD(2);
t361 = qJD(3) * t265;
t444 = -qJD(1) * t361 + qJDD(2);
t117 = qJD(3) * t355 + (t332 + t352) * t268 + t444 * t264;
t368 = qJD(1) * t269;
t118 = t264 * (qJD(2) * (qJD(3) + t368) + t352) - t444 * t268;
t260 = sin(pkin(10));
t261 = cos(pkin(10));
t58 = t117 * t260 + t261 * t118;
t59 = t117 * t261 - t118 * t260;
t369 = qJD(1) * t265;
t342 = t264 * t369;
t192 = t342 - t355;
t365 = qJD(2) * t264;
t194 = t268 * t369 + t365;
t126 = t261 * t192 + t194 * t260;
t304 = -t192 * t260 + t261 * t194;
t74 = t126 * t263 + t267 * t304;
t15 = qJD(6) * t74 + t263 * t59 - t267 * t58;
t226 = -qJD(3) + t368;
t354 = -qJD(6) - t226;
t407 = t354 * t74;
t457 = t15 + t407;
t356 = qJD(6) * t267;
t357 = qJD(6) * t263;
t326 = -t126 * t356 - t263 * t58 - t267 * t59 + t304 * t357;
t443 = -t267 * t126 + t263 * t304;
t406 = t354 * t443;
t456 = t326 + t406;
t241 = pkin(7) * t352;
t317 = qJDD(2) * pkin(2) - pkin(7) * t332 - t241;
t256 = g(3) * t269;
t266 = sin(qJ(1));
t270 = cos(qJ(1));
t318 = g(1) * t270 + g(2) * t266;
t437 = -t318 * t265 + t256;
t455 = -qJD(3) * pkin(8) * t226 - t317 + t437;
t454 = -t443 ^ 2 + t74 ^ 2;
t257 = qJ(3) + pkin(10);
t245 = sin(t257);
t246 = cos(t257);
t300 = t245 * t267 - t246 * t263;
t249 = t269 * qJDD(1);
t434 = -t265 * t353 + t249;
t184 = qJDD(3) - t434;
t319 = pkin(2) * t265 - pkin(8) * t269;
t196 = t319 * qJD(2);
t203 = -pkin(2) * t269 - pkin(8) * t265 - pkin(1);
t138 = qJD(1) * t196 + qJDD(1) * t203;
t132 = t268 * t138;
t182 = t203 * qJD(1);
t243 = pkin(7) * t368;
t207 = qJD(2) * pkin(8) + t243;
t137 = t182 * t264 + t207 * t268;
t165 = pkin(7) * t434 + qJDD(2) * pkin(8);
t27 = pkin(3) * t184 - qJ(4) * t117 - qJD(3) * t137 - qJD(4) * t194 - t264 * t165 + t132;
t360 = qJD(3) * t268;
t362 = qJD(3) * t264;
t288 = t264 * t138 + t268 * t165 + t182 * t360 - t207 * t362;
t31 = -qJ(4) * t118 - qJD(4) * t192 + t288;
t11 = -t260 * t31 + t261 * t27;
t345 = -qJDD(5) + t11;
t425 = pkin(4) + pkin(5);
t3 = -pkin(9) * t59 - t184 * t425 - t345;
t210 = t226 * qJD(5);
t12 = t260 * t27 + t261 * t31;
t350 = t184 * qJ(5) + t12;
t7 = -t210 + t350;
t5 = pkin(9) * t58 + t7;
t343 = -t263 * t5 + t267 * t3;
t206 = -qJD(2) * pkin(2) + pkin(7) * t369;
t286 = -pkin(3) * t192 - qJD(4) - t206;
t278 = qJ(5) * t304 + t286;
t36 = -t126 * t425 + t278;
t414 = g(3) * t265;
t385 = t266 * t269;
t154 = t245 * t385 + t246 * t270;
t381 = t270 * t245;
t155 = t246 * t385 - t381;
t432 = t154 * t267 - t155 * t263;
t156 = -t266 * t246 + t269 * t381;
t382 = t269 * t270;
t157 = t245 * t266 + t246 * t382;
t97 = t156 * t267 - t157 * t263;
t453 = -g(1) * t97 - g(2) * t432 - t300 * t414 - t36 * t74 + t343;
t451 = t74 * t443;
t334 = t260 * t362;
t341 = t264 * t368;
t394 = t261 * t268;
t377 = -t260 * t341 - t261 * t360 + t368 * t394 + t334;
t299 = t245 * t263 + t246 * t267;
t307 = t154 * t263 + t155 * t267;
t98 = t156 * t263 + t157 * t267;
t449 = -g(1) * t98 - g(2) * t307 - t299 * t414 - t36 * t443;
t448 = pkin(9) * t126;
t447 = t126 * t226;
t186 = t260 * t268 + t261 * t264;
t170 = t186 * qJD(3);
t378 = t186 * t368 - t170;
t363 = qJD(2) * t269;
t339 = t264 * t363;
t446 = t265 * t360 + t339;
t445 = -t243 + (-t341 + t362) * pkin(3);
t442 = t304 ^ 2;
t441 = pkin(9) * t304;
t384 = t268 * t269;
t296 = pkin(3) * t265 - qJ(4) * t384;
t195 = t319 * qJD(1);
t374 = pkin(7) * t342 + t268 * t195;
t112 = qJD(1) * t296 + t374;
t176 = t264 * t195;
t388 = t265 * t268;
t391 = t264 * t269;
t130 = t176 + (-pkin(7) * t388 - qJ(4) * t391) * qJD(1);
t262 = -qJ(4) - pkin(8);
t329 = qJD(3) * t262;
t359 = qJD(4) * t268;
t163 = t264 * t329 + t359;
t164 = -qJD(4) * t264 + t268 * t329;
t403 = (t112 - t164) * t261 + (-t130 + t163) * t260;
t107 = t261 * t163 + t260 * t164;
t62 = t260 * t112 + t261 * t130;
t52 = qJ(5) * t369 + t62;
t402 = t107 - t52;
t100 = -qJ(4) * t192 + t137;
t396 = t260 * t100;
t136 = t268 * t182 - t207 * t264;
t99 = -qJ(4) * t194 + t136;
t47 = t261 * t99 - t396;
t379 = qJD(5) - t47;
t435 = qJ(5) * t377 - qJD(5) * t186 + t445;
t433 = qJ(5) * t59 + qJD(5) * t304;
t49 = pkin(4) * t126 - t278;
t431 = g(1) * t156 + g(2) * t154 + t245 * t414 - t304 * t49 + t345;
t430 = t269 * t318 + t414;
t426 = -0.2e1 * pkin(1);
t424 = pkin(4) * t58;
t421 = pkin(4) * t184;
t420 = pkin(7) * t264;
t419 = g(1) * t266;
t416 = g(2) * t270;
t185 = t260 * t264 - t394;
t305 = t267 * t185 - t186 * t263;
t413 = qJD(6) * t305 - t263 * t378 - t267 * t377;
t124 = t185 * t263 + t186 * t267;
t412 = qJD(6) * t124 - t263 * t377 + t267 * t378;
t231 = pkin(7) * t384;
t364 = qJD(2) * t265;
t375 = t268 * t196 + t364 * t420;
t63 = -t265 * t359 + t296 * qJD(2) + (-t231 + (qJ(4) * t265 - t203) * t264) * qJD(3) + t375;
t376 = t264 * t196 + t203 * t360;
t78 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t388 + (-qJD(4) * t265 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t269) * t264 + t376;
t26 = t260 * t63 + t261 * t78;
t411 = t378 * t425 - t435;
t410 = -pkin(4) * t378 + t435;
t395 = t261 * t100;
t87 = -pkin(3) * t226 + t99;
t44 = t260 * t87 + t395;
t46 = t260 * t99 + t395;
t408 = t304 * t46;
t40 = -t226 * qJ(5) + t44;
t28 = t40 + t448;
t405 = t263 * t28;
t404 = pkin(4) * t369 + t403;
t401 = t117 * t264;
t399 = t192 * t226;
t398 = t194 * t226;
t397 = t194 * t268;
t393 = t264 * t265;
t392 = t264 * t266;
t390 = t264 * t270;
t389 = t265 * t266;
t387 = t265 * t270;
t386 = t266 * t268;
t383 = t268 * t270;
t240 = pkin(3) * t268 + pkin(2);
t208 = t269 * t240;
t380 = -t379 + t441;
t188 = t268 * t203;
t133 = -qJ(4) * t388 + t188 + (-pkin(3) - t420) * t269;
t372 = t264 * t203 + t231;
t139 = -qJ(4) * t393 + t372;
t80 = t260 * t133 + t261 * t139;
t204 = t262 * t264;
t205 = t262 * t268;
t141 = t260 * t204 - t261 * t205;
t228 = pkin(3) * t393;
t371 = t265 * pkin(7) + t228;
t258 = t265 ^ 2;
t370 = -t269 ^ 2 + t258;
t367 = qJD(2) * t192;
t366 = qJD(2) * t194;
t43 = t261 * t87 - t396;
t322 = qJD(5) - t43;
t22 = t226 * t425 + t322 - t441;
t351 = t22 * t356 + t263 * t3 + t267 * t5;
t347 = t425 * t265;
t346 = t264 * t382;
t344 = pkin(3) * t446 + pkin(7) * t363;
t236 = -pkin(3) * t261 - pkin(4);
t340 = t226 * t355;
t338 = t269 * t355;
t337 = t226 * t362;
t336 = t226 * t360;
t25 = -t260 * t78 + t261 * t63;
t79 = t133 * t261 - t260 * t139;
t140 = -t261 * t204 - t205 * t260;
t327 = t354 ^ 2;
t325 = -qJD(3) * t182 - t165;
t233 = g(1) * t389;
t323 = -g(2) * t387 + t233;
t109 = -pkin(9) * t186 + t140;
t321 = pkin(9) * t378 - qJD(6) * t109 - t402;
t110 = pkin(9) * t185 + t141;
t320 = -pkin(9) * t377 - qJD(1) * t347 + qJD(6) * t110 - t403;
t75 = -qJ(5) * t269 + t80;
t161 = -t260 * t393 + t261 * t388;
t315 = qJ(5) * t161 - t371;
t77 = t269 * pkin(4) - t79;
t314 = t207 * t360 - t132;
t313 = -pkin(3) * t194 - qJ(5) * t126;
t312 = pkin(4) * t246 + qJ(5) * t245;
t9 = t263 * t22 + t267 * t28;
t48 = pkin(5) * t269 - pkin(9) * t161 + t77;
t160 = t186 * t265;
t50 = pkin(9) * t160 + t75;
t311 = -t263 * t50 + t267 * t48;
t310 = t263 * t48 + t267 * t50;
t309 = -pkin(8) * t184 + qJD(3) * t206;
t308 = -t126 ^ 2 - t442;
t306 = t267 * t160 - t161 * t263;
t103 = t160 * t263 + t161 * t267;
t224 = -pkin(5) + t236;
t234 = pkin(3) * t260 + qJ(5);
t303 = t224 * t267 - t234 * t263;
t302 = t224 * t263 + t234 * t267;
t19 = qJ(5) * t364 - qJD(5) * t269 + t26;
t297 = qJ(5) * t186 + t240;
t293 = t28 * t357 - t351;
t291 = -pkin(7) * qJDD(2) + t353 * t426;
t172 = t264 * t385 + t383;
t290 = t184 * t264 - t336;
t289 = t184 * t268 + t337;
t272 = qJD(1) ^ 2;
t287 = pkin(1) * t272 + t318;
t271 = qJD(2) ^ 2;
t285 = pkin(7) * t271 + qJDD(1) * t426 + t416;
t284 = t270 * pkin(1) + pkin(3) * t392 + t266 * pkin(7) + t240 * t382 - t262 * t387;
t283 = -pkin(3) * t118 - qJDD(4) + t317;
t105 = t170 * t265 + t260 * t339 - t261 * t338;
t282 = -qJ(5) * t105 + qJD(5) * t161 - t344;
t281 = t262 * t389 + pkin(3) * t390 + t270 * pkin(7) + (-pkin(1) - t208) * t266;
t275 = -t107 * t126 + t140 * t59 - t141 * t58 - t430;
t274 = t283 + t433;
t273 = -t283 + t437;
t230 = pkin(3) * t386;
t181 = -qJDD(6) + t184;
t175 = t268 * t382 + t392;
t174 = -t346 + t386;
t173 = -t266 * t384 + t390;
t121 = pkin(4) * t185 - t297;
t104 = -t260 * t338 - t261 * t446 + t265 * t334;
t89 = -t185 * t425 + t297;
t88 = pkin(4) * t160 - t315;
t76 = -t160 * t425 + t315;
t51 = pkin(4) * t304 - t313;
t39 = pkin(4) * t226 + t322;
t38 = -t304 * t425 + t313;
t37 = -pkin(4) * t104 - t282;
t34 = t46 + t448;
t33 = qJD(6) * t103 + t267 * t104 - t105 * t263;
t32 = qJD(6) * t306 - t104 * t263 - t105 * t267;
t21 = -pkin(4) * t364 - t25;
t18 = t104 * t425 + t282;
t17 = -pkin(9) * t104 + t19;
t16 = pkin(9) * t105 - qJD(2) * t347 - t25;
t13 = -t274 + t424;
t10 = -t345 - t421;
t8 = t22 * t267 - t405;
t6 = -t425 * t58 + t274;
t1 = [qJDD(1), -t416 + t419, t318, qJDD(1) * t258 + 0.2e1 * t265 * t332, 0.2e1 * t249 * t265 - 0.2e1 * t353 * t370, qJDD(2) * t265 + t269 * t271, qJDD(2) * t269 - t265 * t271, 0, t291 * t265 + (-t285 + t419) * t269, t265 * t285 + t269 * t291 - t233, t117 * t388 + (-t264 * t361 + t338) * t194 (-t192 * t268 - t194 * t264) * t363 + (-t401 - t118 * t268 + (t192 * t264 - t397) * qJD(3)) * t265 (-t117 - t340) * t269 + (t289 + t366) * t265 (t226 * t365 + t118) * t269 + (-t290 - t367) * t265, -t184 * t269 - t226 * t364 -(-t203 * t362 + t375) * t226 + t188 * t184 - g(1) * t173 - g(2) * t175 + ((t336 + t367) * pkin(7) + (-pkin(7) * t184 + qJD(2) * t206 - t325) * t264 + t314) * t269 + (pkin(7) * t118 + qJD(2) * t136 + t206 * t360 - t264 * t317) * t265, t376 * t226 - t372 * t184 - g(1) * t172 - g(2) * t174 + (t206 * t355 + (-t337 + t366) * pkin(7) + t288) * t269 + (-t206 * t362 - t137 * qJD(2) - t317 * t268 + (t117 - t340) * pkin(7)) * t265, t104 * t44 + t105 * t43 - t11 * t161 - t12 * t160 - t126 * t26 - t25 * t304 - t58 * t80 - t59 * t79 + t323, -g(1) * t281 - g(2) * t284 + t11 * t79 + t12 * t80 + t43 * t25 + t44 * t26 - t283 * t371 - t286 * t344, g(1) * t155 - g(2) * t157 + t10 * t269 - t104 * t49 + t126 * t37 + t13 * t160 - t184 * t77 + t21 * t226 - t364 * t39 + t58 * t88, t10 * t161 + t104 * t40 - t105 * t39 - t126 * t19 - t160 * t7 + t21 * t304 - t58 * t75 + t59 * t77 + t323, g(1) * t154 - g(2) * t156 + t105 * t49 - t13 * t161 + t184 * t75 - t19 * t226 - t269 * t7 - t304 * t37 + t364 * t40 - t59 * t88, t7 * t75 + t40 * t19 + t13 * t88 + t49 * t37 + t10 * t77 + t39 * t21 - g(1) * (-pkin(4) * t155 - qJ(5) * t154 + t281) - g(2) * (pkin(4) * t157 + qJ(5) * t156 + t284) -t103 * t326 + t32 * t74, -t103 * t15 - t306 * t326 - t32 * t443 - t33 * t74, -t103 * t181 - t269 * t326 - t32 * t354 - t364 * t74, -t15 * t269 - t181 * t306 + t33 * t354 + t364 * t443, -t181 * t269 + t354 * t364 -(t16 * t267 - t17 * t263) * t354 - t311 * t181 + t343 * t269 - t8 * t364 + t18 * t443 + t76 * t15 - t6 * t306 + t36 * t33 + g(1) * t307 - g(2) * t98 + (-t269 * t9 + t310 * t354) * qJD(6) (qJD(6) * t311 + t16 * t263 + t17 * t267) * t354 + t310 * t181 + t293 * t269 + t9 * t364 + t18 * t74 - t76 * t326 + t6 * t103 + t36 * t32 + g(1) * t432 - g(2) * t97; 0, 0, 0, -t265 * t272 * t269, t370 * t272, t352, t249, qJDD(2), t265 * t287 - t241 - t256, t414 + (-pkin(7) * qJDD(1) + t287) * t269, -t226 * t397 + t401 (t117 + t399) * t268 + (-t118 + t398) * t264 (-t194 * t265 + t226 * t384) * qJD(1) + t290 (t192 * t265 - t226 * t391) * qJD(1) + t289, t226 * t369, -pkin(2) * t118 + t374 * t226 + t309 * t264 + (-t136 * t265 + (-pkin(7) * t192 - t206 * t264) * t269) * qJD(1) - t455 * t268, -pkin(2) * t117 - t176 * t226 + t309 * t268 + (-t206 * t384 + t137 * t265 + (-t194 * t269 + t226 * t388) * pkin(7)) * qJD(1) + t455 * t264, -t11 * t186 - t12 * t185 + t126 * t62 + t304 * t403 + t377 * t43 + t378 * t44 + t275, t12 * t141 - t11 * t140 + t283 * t240 - g(3) * (-t262 * t265 + t208) + (t107 - t62) * t44 - t403 * t43 - t445 * t286 + t318 * (t240 * t265 + t262 * t269) t121 * t58 + t126 * t410 + t13 * t185 - t140 * t184 + t226 * t404 - t246 * t437 + t369 * t39 - t378 * t49, t10 * t186 + t126 * t52 - t185 * t7 + t304 * t404 - t377 * t39 + t378 * t40 + t275, -t121 * t59 - t13 * t186 + t141 * t184 - t226 * t402 - t245 * t437 - t304 * t410 - t369 * t40 + t377 * t49, -g(3) * t208 + t10 * t140 + t13 * t121 + t7 * t141 + t410 * t49 + t402 * t40 + t404 * t39 + (-g(3) * t312 + t262 * t318) * t269 + (g(3) * t262 + t318 * (t240 + t312)) * t265, -t124 * t326 + t413 * t74, -t124 * t15 - t305 * t326 - t412 * t74 - t413 * t443, -t124 * t181 - t354 * t413 + t369 * t74, -t181 * t305 + t354 * t412 - t369 * t443, -t354 * t369 -(t109 * t267 - t110 * t263) * t181 + t89 * t15 - t6 * t305 + t411 * t443 + t412 * t36 - (t263 * t321 - t267 * t320) * t354 + t8 * t369 - t437 * t299 (t109 * t263 + t110 * t267) * t181 - t89 * t326 + t6 * t124 + t411 * t74 + t413 * t36 - (t263 * t320 + t267 * t321) * t354 - t9 * t369 - t437 * t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194 * t192, -t192 ^ 2 + t194 ^ 2, t117 - t399, -t118 - t398, t184, -g(1) * t174 + g(2) * t172 - t137 * t226 - t194 * t206 + (t325 + t414) * t264 - t314, g(1) * t175 - g(2) * t173 + g(3) * t388 - t136 * t226 + t192 * t206 - t288, t304 * t44 - t408 + (-t260 * t58 - t261 * t59) * pkin(3) + (t47 - t43) * t126, -g(1) * t230 + t43 * t46 - t44 * t47 + (g(2) * t383 + t11 * t261 + t12 * t260 + t194 * t286 + t264 * t430) * pkin(3), -t126 * t51 - t226 * t46 + (pkin(4) - t236) * t184 + t431, -t234 * t58 + t236 * t59 + t304 * t40 - t408 + (-t379 + t39) * t126, -g(1) * t157 - g(2) * t155 - t126 * t49 + t184 * t234 + t226 * t47 - t246 * t414 + t304 * t51 - 0.2e1 * t210 + t350, t7 * t234 + t10 * t236 - t49 * t51 - t39 * t46 - g(1) * (-pkin(3) * t346 - pkin(4) * t156 + qJ(5) * t157 + t230) - g(2) * (-pkin(3) * t172 - pkin(4) * t154 + qJ(5) * t155) - g(3) * (-t228 + (-pkin(4) * t245 + qJ(5) * t246) * t265) + t379 * t40, -t451, -t454, t456, t457, t181, -t303 * t181 - t38 * t443 - (t263 * t380 - t267 * t34) * t354 + (t302 * t354 + t9) * qJD(6) - t453, t302 * t181 - t38 * t74 - (t263 * t34 + t267 * t380) * t354 + (t303 * t354 - t405) * qJD(6) + t351 + t449; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t308, t126 * t44 + t304 * t43 + t273, -t226 * t304 + t58, t308, -t59 - t447, t126 * t40 - t304 * t39 + t273 + t424 - t433, 0, 0, 0, 0, 0, -t15 + t407, t326 - t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126 * t304 - t184, t59 - t447, -t226 ^ 2 - t442, t226 * t40 - t421 - t431, 0, 0, 0, 0, 0, -t267 * t181 - t263 * t327 - t304 * t443, t263 * t181 - t267 * t327 - t304 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451, t454, -t456, -t457, -t181 (-t354 - qJD(6)) * t9 + t453, -t354 * t8 + t293 - t449;];
tau_reg  = t1;
