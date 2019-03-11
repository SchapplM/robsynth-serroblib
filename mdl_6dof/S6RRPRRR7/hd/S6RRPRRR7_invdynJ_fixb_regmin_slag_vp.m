% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:50
% EndTime: 2019-03-09 14:00:05
% DurationCPUTime: 6.33s
% Computational Cost: add. (6633->545), mult. (13943->698), div. (0->0), fcn. (10197->12), ass. (0->292)
t228 = sin(qJ(4));
t233 = cos(qJ(4));
t234 = cos(qJ(2));
t334 = qJD(1) * t234;
t229 = sin(qJ(2));
t335 = qJD(1) * t229;
t145 = -t228 * t334 + t233 * t335;
t217 = qJD(2) - qJD(4);
t227 = sin(qJ(5));
t232 = cos(qJ(5));
t107 = t145 * t227 + t232 * t217;
t150 = t228 * t229 + t233 * t234;
t394 = t150 * qJD(1);
t411 = qJD(5) + t394;
t426 = t107 * t411;
t328 = qJD(4) * t233;
t329 = qJD(4) * t228;
t331 = qJD(2) * t234;
t425 = t228 * t331 + t229 * t328 - t234 * t329;
t203 = pkin(7) * t335;
t156 = pkin(8) * t335 - t203;
t424 = qJD(3) - t156;
t216 = qJDD(2) - qJDD(4);
t236 = -pkin(2) - pkin(3);
t311 = t236 * qJD(2);
t124 = t311 + t424;
t204 = pkin(7) * t334;
t158 = -pkin(8) * t334 + t204;
t220 = qJD(2) * qJ(3);
t146 = t158 + t220;
t320 = qJD(1) * qJD(2);
t305 = t234 * t320;
t186 = pkin(7) * t305;
t207 = t229 * qJDD(1);
t198 = pkin(7) * t207;
t304 = qJDD(3) + t186 + t198;
t413 = -t305 - t207;
t89 = pkin(8) * t413 + t236 * qJDD(2) + t304;
t319 = t234 * qJDD(1);
t199 = pkin(7) * t319;
t218 = qJDD(2) * qJ(3);
t219 = qJD(2) * qJD(3);
t306 = t229 * t320;
t118 = -pkin(7) * t306 + t199 + t218 + t219;
t412 = t306 - t319;
t90 = pkin(8) * t412 + t118;
t267 = t124 * t329 + t146 * t328 + t228 * t90 - t233 * t89;
t23 = pkin(4) * t216 + t267;
t109 = t145 * t232 - t217 * t227;
t251 = t228 * t412 - t233 * t413;
t254 = t150 * qJD(4);
t59 = -qJD(1) * t254 + t251;
t299 = t232 * t216 + t227 * t59;
t36 = qJD(5) * t109 + t299;
t10 = pkin(5) * t36 + t23;
t231 = cos(qJ(6));
t226 = sin(qJ(6));
t349 = t226 * t227;
t149 = -t231 * t232 + t349;
t223 = qJ(5) + qJ(6);
t210 = cos(t223);
t230 = sin(qJ(1));
t342 = t230 * t234;
t345 = t229 * t233;
t134 = t228 * t342 - t230 * t345;
t235 = cos(qJ(1));
t344 = t229 * t235;
t346 = t228 * t234;
t136 = -t233 * t344 + t235 * t346;
t252 = g(1) * t136 + g(2) * t134 + g(3) * t150;
t348 = t226 * t232;
t152 = t227 * t231 + t348;
t318 = qJD(5) + qJD(6);
t373 = (t318 + t394) * t152;
t77 = t124 * t233 - t228 * t146;
t64 = pkin(4) * t217 - t77;
t40 = pkin(5) * t107 + t64;
t147 = -qJD(1) * pkin(1) - pkin(2) * t334 - qJ(3) * t335;
t117 = pkin(3) * t334 - t147;
t54 = pkin(4) * t394 - pkin(9) * t145 + t117;
t78 = t228 * t124 + t233 * t146;
t65 = -pkin(9) * t217 + t78;
t27 = -t227 * t65 + t232 * t54;
t18 = -pkin(10) * t109 + t27;
t13 = pkin(5) * t411 + t18;
t28 = t227 * t54 + t232 * t65;
t19 = -pkin(10) * t107 + t28;
t5 = t13 * t231 - t19 * t226;
t423 = -t10 * t149 + t5 * t145 - t252 * t210 - t373 * t40;
t209 = sin(t223);
t325 = qJD(5) * t232;
t396 = -qJD(6) * t232 - t325;
t249 = t396 * t231;
t375 = t149 * t394 + t318 * t349 + t249;
t370 = t19 * t231;
t6 = t13 * t226 + t370;
t422 = -t10 * t152 - t6 * t145 + t252 * t209 + t375 * t40;
t132 = qJD(6) + t411;
t272 = t107 * t226 - t231 * t109;
t60 = qJD(1) * t425 + qJDD(1) * t150 - t233 * t306;
t57 = qJDD(5) + t60;
t56 = qJDD(6) + t57;
t371 = t152 * t56;
t421 = t132 * t375 - t145 * t272 - t371;
t372 = t149 * t56;
t44 = t231 * t107 + t109 * t226;
t420 = t132 * t373 - t44 * t145 + t372;
t326 = qJD(5) * t227;
t35 = -t145 * t326 - t227 * t216 - t217 * t325 + t232 * t59;
t399 = t411 * t109;
t419 = t227 * (t36 + t399) + (-t35 + t426) * t232;
t323 = qJD(6) * t231;
t324 = qJD(6) * t226;
t8 = -t107 * t323 - t109 * t324 - t226 * t36 + t231 * t35;
t418 = t8 * t152 + t272 * t375;
t246 = qJD(6) * t272 - t226 * t35 - t231 * t36;
t417 = -t8 * t149 + t152 * t246 + t272 * t373 + t375 * t44;
t416 = t272 * t44;
t297 = -t228 * qJ(3) + t233 * t236;
t127 = t233 * qJD(3) + qJD(4) * t297;
t96 = t156 * t233 + t158 * t228;
t362 = t127 - t96;
t410 = t272 ^ 2 - t44 ^ 2;
t409 = t132 * t44 + t8;
t17 = t19 * t324;
t151 = t345 - t346;
t381 = g(3) * t151;
t135 = t150 * t230;
t92 = -t135 * t210 - t209 * t235;
t137 = t150 * t235;
t94 = t137 * t210 - t230 * t209;
t408 = g(1) * t94 - g(2) * t92 + t210 * t381 + t40 * t44 + t17;
t271 = t135 * t209 - t210 * t235;
t208 = t229 * qJD(3);
t224 = qJDD(1) * pkin(1);
t275 = pkin(2) * t319 - qJ(3) * t413 + qJD(1) * t208 + t224;
t290 = t229 * t311;
t62 = pkin(3) * t319 + qJD(1) * t290 + t275;
t14 = t60 * pkin(4) - t59 * pkin(9) + t62;
t12 = t232 * t14;
t257 = t124 * t328 - t146 * t329 + t228 * t89 + t233 * t90;
t22 = -pkin(9) * t216 + t257;
t2 = t57 * pkin(5) - t35 * pkin(10) - qJD(5) * t28 - t227 * t22 + t12;
t262 = -t227 * t14 - t232 * t22 - t54 * t325 + t326 * t65;
t3 = -pkin(10) * t36 - t262;
t313 = t231 * t2 - t226 * t3;
t93 = -t137 * t209 - t210 * t230;
t407 = -g(1) * t93 + g(2) * t271 - qJD(6) * t6 + t209 * t381 + t40 * t272 + t313;
t406 = -t132 * t272 + t246;
t368 = t232 * t57;
t400 = t411 ^ 2;
t405 = -t107 * t145 + t227 * t400 - t368;
t369 = t227 * t57;
t404 = -t109 * t145 + t232 * t400 + t369;
t365 = t35 * t227;
t403 = t232 * t399 + t365;
t337 = t234 * pkin(2) + t229 * qJ(3);
t163 = -pkin(1) - t337;
t402 = t411 * t64;
t338 = t233 * qJ(3) + t228 * t236;
t361 = qJD(4) * t338 + t158 * t233 + t228 * t424;
t401 = t145 * t217 + t60;
t81 = t152 * t151;
t398 = -t226 * t326 - t227 * t324;
t354 = t394 * t227;
t397 = (t326 + t354) * pkin(5);
t382 = g(2) * t235;
t385 = g(1) * t230;
t395 = -t382 + t385;
t383 = g(2) * t230;
t384 = g(1) * t235;
t285 = t383 + t384;
t360 = pkin(7) * qJDD(2);
t390 = (qJD(1) * t163 + t147) * qJD(2) - t360;
t389 = pkin(7) - pkin(8);
t388 = pkin(9) + pkin(10);
t387 = pkin(5) * t232;
t386 = g(1) * t135;
t155 = -pkin(9) + t338;
t379 = pkin(10) - t155;
t87 = t145 * pkin(4) + pkin(9) * t394;
t378 = t227 * t87 + t232 * t77;
t194 = qJ(3) * t334;
t133 = t236 * t335 + t194;
t58 = t133 - t87;
t377 = t227 * t58 + t232 * t96;
t367 = t27 * t145;
t366 = t28 * t145;
t166 = t389 * t229;
t168 = t389 * t234;
t111 = t166 * t228 + t168 * t233;
t104 = t232 * t111;
t148 = t234 * pkin(3) - t163;
t75 = pkin(4) * t150 - pkin(9) * t151 + t148;
t364 = t227 * t75 + t104;
t363 = -t397 + t361;
t359 = qJDD(2) * pkin(2);
t101 = qJD(2) * t150 - t254;
t358 = t101 * t232;
t357 = t132 * t145;
t356 = t411 * t145;
t355 = t394 * t217;
t353 = t145 * t394;
t351 = t151 * t227;
t350 = t151 * t232;
t347 = t227 * t101;
t238 = qJD(1) ^ 2;
t343 = t229 * t238;
t339 = qJ(3) * t331 + t208;
t221 = t229 ^ 2;
t222 = t234 ^ 2;
t336 = t221 - t222;
t333 = qJD(2) * t229;
t332 = qJD(2) * t233;
t330 = qJD(4) * t132;
t327 = qJD(5) * t411;
t316 = pkin(10) * t354;
t314 = t234 * t343;
t312 = qJD(5) * t388;
t307 = -t145 ^ 2 + t394 ^ 2;
t303 = qJD(6) * t13 + t3;
t301 = qJD(5) * t379;
t300 = -qJD(5) * t54 - t22;
t298 = -qJD(2) * pkin(2) + qJD(3);
t291 = t217 ^ 2;
t289 = -t78 + t397;
t154 = pkin(4) - t297;
t119 = t379 * t227;
t288 = -qJD(6) * t119 - t232 * t127 - t227 * t301 - t316 + t377;
t120 = t379 * t232;
t266 = -pkin(10) * t232 * t394 - pkin(5) * t145;
t53 = t232 * t58;
t287 = -qJD(6) * t120 + t227 * t362 - t232 * t301 + t266 + t53;
t237 = qJD(2) ^ 2;
t286 = pkin(7) * t237 + t382;
t284 = g(2) * t135 + t381;
t165 = t388 * t227;
t282 = qJD(6) * t165 + t227 * t312 + t316 + t378;
t167 = t388 * t232;
t80 = t232 * t87;
t281 = qJD(6) * t167 - t227 * t77 + t232 * t312 - t266 + t80;
t278 = pkin(2) * t229 - qJ(3) * t234;
t162 = t203 + t298;
t164 = t204 + t220;
t270 = t162 * t234 - t164 * t229;
t269 = t166 * t233 - t168 * t228;
t268 = g(1) * t344 - g(3) * t234 + t229 * t383 - t198;
t264 = -0.2e1 * pkin(1) * t320 - t360;
t261 = -pkin(9) * t57 + t402;
t260 = t151 * t325 + t347;
t259 = -t151 * t326 + t358;
t100 = -t229 * t332 + t425;
t115 = t290 + t339;
t39 = t100 * pkin(4) - t101 * pkin(9) + t115;
t157 = t389 * t333;
t159 = qJD(2) * t168;
t49 = qJD(4) * t269 - t233 * t157 + t228 * t159;
t258 = -t111 * t326 + t227 * t39 + t232 * t49 + t75 * t325;
t256 = -t155 * t57 - t402;
t255 = -qJDD(3) + t268;
t250 = -t286 + 0.2e1 * t224;
t247 = -t23 + t252;
t138 = pkin(2) * t333 - t339;
t86 = pkin(2) * t306 - t275;
t245 = -qJD(1) * t138 - qJDD(1) * t163 - t286 - t86;
t244 = pkin(9) * t327 - t247;
t243 = t155 * t327 + t247;
t50 = qJD(4) * t111 - t228 * t157 - t233 * t159;
t242 = t117 * t145 - t252 + t267;
t240 = -g(1) * t137 - t117 * t394 + t257 - t284;
t131 = t304 - t359;
t239 = qJD(2) * t270 + t118 * t234 + t131 * t229 - t285;
t197 = -pkin(4) - t387;
t188 = g(1) * t342;
t153 = pkin(2) * t335 - t194;
t144 = t227 * t335 + t232 * t332;
t141 = -t227 * t332 + t232 * t335;
t140 = t154 + t387;
t106 = t137 * t232 - t227 * t230;
t105 = -t137 * t227 - t230 * t232;
t82 = t149 * t151;
t76 = pkin(5) * t351 - t269;
t67 = t232 * t75;
t38 = t232 * t39;
t34 = -pkin(10) * t351 + t364;
t29 = pkin(5) * t150 - pkin(10) * t350 - t111 * t227 + t67;
t25 = pkin(5) * t260 + t50;
t16 = t101 * t348 + (t318 * t350 + t347) * t231 + t398 * t151;
t15 = -t149 * t101 - t318 * t81;
t7 = -pkin(10) * t260 + t258;
t4 = -pkin(10) * t358 + t100 * pkin(5) - t227 * t49 + t38 + (-t104 + (pkin(10) * t151 - t75) * t227) * qJD(5);
t1 = [qJDD(1), t395, t285, qJDD(1) * t221 + 0.2e1 * t229 * t305, 0.2e1 * t229 * t319 - 0.2e1 * t320 * t336, qJDD(2) * t229 + t234 * t237, qJDD(2) * t234 - t229 * t237, 0, t229 * t264 + t234 * t250 + t188, t264 * t234 + (-t250 - t385) * t229, t229 * t390 + t245 * t234 + t188 (t221 + t222) * qJDD(1) * pkin(7) + t239, -t390 * t234 + (t245 + t385) * t229, pkin(7) * t239 + t147 * t138 + (-t395 + t86) * t163, t101 * t145 + t151 * t59, -t100 * t145 - t101 * t394 - t150 * t59 - t151 * t60, -t101 * t217 - t151 * t216, t100 * t217 + t150 * t216, 0, -g(2) * t137 + t100 * t117 + t115 * t394 + t148 * t60 + t150 * t62 - t216 * t269 + t217 * t50 + t386, -g(1) * t134 + g(2) * t136 + t101 * t117 + t111 * t216 + t115 * t145 + t148 * t59 + t151 * t62 + t217 * t49, t109 * t259 + t35 * t350 (-t107 * t232 - t109 * t227) * t101 + (-t365 - t232 * t36 + (t107 * t227 - t109 * t232) * qJD(5)) * t151, t109 * t100 + t35 * t150 + t259 * t411 + t350 * t57, -t107 * t100 - t36 * t150 - t260 * t411 - t351 * t57, t100 * t411 + t150 * t57, -g(2) * t106 + t27 * t100 + t50 * t107 - t269 * t36 + t12 * t150 + t38 * t411 + t67 * t57 + (t386 + (-t111 * t411 - t150 * t65 + t151 * t64) * qJD(5)) * t232 + ((-qJD(5) * t75 - t49) * t411 - t111 * t57 + t300 * t150 + t23 * t151 + t64 * t101 + t384) * t227, -t258 * t411 - t364 * t57 + t262 * t150 - t28 * t100 + t50 * t109 - t269 * t35 + t64 * t358 - g(1) * (t135 * t227 - t232 * t235) - g(2) * t105 + (t23 * t232 - t326 * t64) * t151, -t15 * t272 - t8 * t82, -t15 * t44 + t16 * t272 - t246 * t82 - t8 * t81, -t100 * t272 + t132 * t15 + t150 * t8 - t56 * t82, -t100 * t44 - t132 * t16 + t150 * t246 - t56 * t81, t100 * t132 + t150 * t56 (-t226 * t7 + t231 * t4) * t132 + (-t226 * t34 + t231 * t29) * t56 + t313 * t150 + t5 * t100 + t25 * t44 - t76 * t246 + t10 * t81 + t40 * t16 - g(1) * t92 - g(2) * t94 + ((-t226 * t29 - t231 * t34) * t132 - t6 * t150) * qJD(6), t17 * t150 - t6 * t100 - t25 * t272 + t76 * t8 - t10 * t82 + t40 * t15 - g(1) * t271 - g(2) * t93 + (-(-qJD(6) * t34 + t4) * t132 - t29 * t56 - t2 * t150) * t226 + (-(qJD(6) * t29 + t7) * t132 - t34 * t56 - t303 * t150) * t231; 0, 0, 0, -t314, t336 * t238, t207, t319, qJDD(2), pkin(1) * t343 + t268, g(3) * t229 - t199 + (pkin(1) * t238 + t285) * t234, 0.2e1 * t359 + (-t147 * t229 + t153 * t234) * qJD(1) + t255, -t278 * qJDD(1) + ((t164 - t220) * t229 + (-t162 + t298) * t234) * qJD(1), t199 + 0.2e1 * t218 + 0.2e1 * t219 + (qJD(1) * t153 - g(3)) * t229 + (qJD(1) * t147 - t285) * t234, -t270 * qJD(1) * pkin(7) - t131 * pkin(2) - g(3) * t337 + t118 * qJ(3) + t164 * qJD(3) - t147 * t153 + t278 * t285, -t353, t307, qJD(4) * t394 - t251 + t355, t401, t216, -t133 * t394 - t216 * t297 + t217 * t361 + t242, -t133 * t145 + t216 * t338 + t217 * t362 + t240, -t403, t419, -t404, t405, t356, -t53 * t411 + t367 + t154 * t36 + t361 * t107 + (-t362 * t411 + t256) * t227 - t243 * t232, t154 * t35 + t377 * t411 - t366 + t361 * t109 + (-t127 * t411 + t256) * t232 + t243 * t227, -t418, -t417, t421, t420, t357 (t119 * t231 + t120 * t226) * t56 - t140 * t246 + t363 * t44 + (t226 * t288 - t231 * t287) * t132 + t423 -(t119 * t226 - t120 * t231) * t56 + t140 * t8 - t363 * t272 + (t226 * t287 + t231 * t288) * t132 + t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t314, t207, -t221 * t238 - t237, -qJD(2) * t164 + t147 * t335 + t186 - t255 - t359, 0, 0, 0, 0, 0, -t233 * t216 - t228 * t291 - t335 * t394, -t145 * t335 + t228 * t216 - t233 * t291, 0, 0, 0, 0, 0, -t233 * t36 + (-t227 * t328 - t141) * t411 + (-t107 * t217 - t325 * t411 - t369) * t228, -t233 * t35 + (-t232 * t328 + t144) * t411 + (-t109 * t217 + t326 * t411 - t368) * t228, 0, 0, 0, 0, 0 -(t141 * t231 - t144 * t226) * t132 + (-t152 * t330 + t246) * t233 + ((t249 - t398) * t132 - t371 - t217 * t44) * t228 (t141 * t226 + t144 * t231) * t132 + (t149 * t330 - t8) * t233 + (-(t226 * t396 - t227 * t323 - t231 * t326) * t132 + t372 + t217 * t272) * t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, -t307, t59 - t355, -t401, -t216, -t217 * t78 - t242, -t217 * t77 - t240, t403, -t419, t404, -t405, -t356, -pkin(4) * t36 - t78 * t107 - t80 * t411 - t367 + (t411 * t77 + t261) * t227 - t244 * t232, -pkin(4) * t35 - t78 * t109 + t227 * t244 + t232 * t261 + t378 * t411 + t366, t418, t417, -t421, -t420, -t357 (-t165 * t231 - t167 * t226) * t56 - t197 * t246 + t289 * t44 + (t226 * t282 - t231 * t281) * t132 - t423 -(-t165 * t226 + t167 * t231) * t56 + t197 * t8 - t289 * t272 + (t226 * t281 + t231 * t282) * t132 - t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t107, -t107 ^ 2 + t109 ^ 2, t35 + t426, -t299 + (-qJD(5) + t411) * t109, t57, -g(1) * t105 - t64 * t109 + t28 * t411 + t12 + (-qJD(5) * t65 - t382) * t232 + (t284 + t300) * t227, t27 * t411 + t64 * t107 + g(1) * t106 - g(2) * (-t135 * t232 - t227 * t235) + g(3) * t350 + t262, -t416, t410, t409, t406, t56 -(-t18 * t226 - t370) * t132 + (-t109 * t44 - t132 * t324 + t231 * t56) * pkin(5) + t407 (-t132 * t19 - t2) * t226 + (t132 * t18 - t303) * t231 + (t109 * t272 - t132 * t323 - t226 * t56) * pkin(5) + t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t416, t410, t409, t406, t56, t6 * t132 + t407, t5 * t132 - t226 * t2 - t231 * t303 + t408;];
tau_reg  = t1;
