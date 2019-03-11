% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR1
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
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:10:01
% EndTime: 2019-03-09 10:10:16
% DurationCPUTime: 7.51s
% Computational Cost: add. (11209->492), mult. (26937->648), div. (0->0), fcn. (21228->18), ass. (0->272)
t263 = sin(pkin(10));
t265 = cos(pkin(10));
t269 = sin(qJ(2));
t272 = cos(qJ(2));
t203 = -t263 * t269 + t265 * t272;
t189 = t203 * qJD(1);
t205 = t263 * t272 + t265 * t269;
t191 = t205 * qJD(1);
t268 = sin(qJ(4));
t382 = cos(qJ(4));
t145 = t382 * t189 - t191 * t268;
t139 = qJD(6) - t145;
t262 = sin(pkin(11));
t264 = cos(pkin(11));
t267 = sin(qJ(6));
t271 = cos(qJ(6));
t206 = t262 * t271 + t264 * t267;
t204 = t262 * t267 - t271 * t264;
t417 = t139 * t204;
t332 = qJD(1) * qJD(2);
t323 = t272 * t332;
t324 = t269 * t332;
t154 = qJDD(1) * t205 - t263 * t324 + t265 * t323;
t190 = t205 * qJD(2);
t286 = qJD(1) * t190;
t276 = qJDD(1) * t203 - t286;
t294 = -t268 * t189 - t382 * t191;
t278 = qJD(4) * t294 - t268 * t154 + t382 * t276;
t74 = qJDD(6) - t278;
t419 = t139 * t417 - t206 * t74;
t195 = t206 * qJD(6);
t418 = -t206 * t145 + t195;
t400 = t145 * t262;
t416 = pkin(5) * t400;
t415 = pkin(9) * t400;
t253 = t272 * pkin(2);
t246 = t253 + pkin(1);
t174 = -pkin(3) * t203 - t246;
t258 = qJD(2) + qJD(4);
t357 = t294 * t258;
t414 = t278 - t357;
t314 = -t139 * t418 - t204 * t74;
t128 = -t264 * t258 - t262 * t294;
t130 = t258 * t262 - t264 * t294;
t88 = t271 * t128 + t130 * t267;
t366 = t294 * t88;
t413 = t314 - t366;
t412 = qJD(6) - t139;
t369 = qJ(3) + pkin(7);
t229 = t369 * t272;
t215 = qJD(1) * t229;
t196 = t263 * t215;
t228 = t369 * t269;
t214 = qJD(1) * t228;
t367 = qJD(2) * pkin(2);
t200 = -t214 + t367;
t152 = t265 * t200 - t196;
t377 = pkin(8) * t191;
t123 = qJD(2) * pkin(3) + t152 - t377;
t341 = t265 * t215;
t153 = t263 * t200 + t341;
t378 = pkin(8) * t189;
t127 = t153 + t378;
t325 = qJD(4) * t382;
t335 = qJD(4) * t268;
t319 = qJD(2) * t369;
t187 = -qJD(3) * t269 - t272 * t319;
t151 = qJDD(2) * pkin(2) + qJD(1) * t187 - qJDD(1) * t228;
t186 = qJD(3) * t272 - t269 * t319;
t157 = qJD(1) * t186 + qJDD(1) * t229;
t109 = t265 * t151 - t157 * t263;
t78 = qJDD(2) * pkin(3) - pkin(8) * t154 + t109;
t110 = t263 * t151 + t265 * t157;
t83 = pkin(8) * t276 + t110;
t315 = -t123 * t335 - t127 * t325 - t268 * t83 + t382 * t78;
t255 = qJDD(2) + qJDD(4);
t383 = -pkin(4) * t255 + qJDD(5);
t24 = -t315 + t383;
t75 = t382 * t154 + t189 * t325 - t191 * t335 + t268 * t276;
t63 = -t264 * t255 + t262 * t75;
t13 = pkin(5) * t63 + t24;
t259 = qJ(2) + pkin(10);
t251 = qJ(4) + t259;
t243 = cos(t251);
t257 = pkin(11) + qJ(6);
t250 = cos(t257);
t80 = t268 * t123 + t382 * t127;
t71 = qJ(5) * t258 + t80;
t218 = -qJD(1) * t246 + qJD(3);
t160 = -pkin(3) * t189 + t218;
t85 = -pkin(4) * t145 + qJ(5) * t294 + t160;
t38 = -t262 * t71 + t264 * t85;
t28 = -pkin(5) * t145 - pkin(9) * t130 + t38;
t39 = t262 * t85 + t264 * t71;
t31 = -pkin(9) * t128 + t39;
t306 = t267 * t31 - t271 * t28;
t373 = g(3) * t250;
t242 = sin(t251);
t273 = cos(qJ(1));
t349 = t242 * t273;
t270 = sin(qJ(1));
t350 = t242 * t270;
t396 = g(1) * t349 + g(2) * t350;
t79 = t382 * t123 - t268 * t127;
t70 = -t258 * pkin(4) + qJD(5) - t79;
t55 = t128 * pkin(5) + t70;
t411 = t13 * t204 - t243 * t373 + t396 * t250 - t306 * t294 + t418 * t55;
t249 = sin(t257);
t313 = g(1) * t273 + g(2) * t270;
t300 = t242 * t313;
t374 = g(3) * t249;
t8 = t267 * t28 + t271 * t31;
t410 = t13 * t206 + t243 * t374 - t249 * t300 - t8 * t294 - t417 * t55;
t333 = qJD(6) * t271;
t334 = qJD(6) * t267;
t64 = t255 * t262 + t264 * t75;
t19 = -t128 * t333 - t130 * t334 - t267 * t63 + t271 * t64;
t303 = t128 * t267 - t130 * t271;
t409 = t19 * t206 + t303 * t417;
t365 = t294 * t303;
t408 = -t365 - t419;
t20 = -qJD(6) * t303 + t267 * t64 + t271 * t63;
t407 = -t19 * t204 - t206 * t20 + t303 * t418 + t417 * t88;
t406 = pkin(5) * t294;
t405 = t139 * t88;
t404 = t145 * t38;
t355 = t145 * t258;
t403 = t75 - t355;
t402 = t139 * t294;
t401 = t139 * t303;
t397 = t294 * t145;
t339 = t243 * pkin(4) + t242 * qJ(5);
t395 = pkin(3) * cos(t259) + t253 + t339;
t394 = -g(3) * t243 + t396;
t393 = -t145 ^ 2 + t294 ^ 2;
t279 = t123 * t325 - t127 * t335 + t268 * t78 + t382 * t83;
t22 = t255 * qJ(5) + t258 * qJD(5) + t279;
t329 = pkin(2) * t324 + qJDD(3);
t125 = pkin(3) * t286 + qJDD(1) * t174 + t329;
t27 = -pkin(4) * t278 - t75 * qJ(5) + qJD(5) * t294 + t125;
t5 = -t22 * t262 + t264 * t27;
t392 = t145 * t39 - t5;
t348 = t243 * t262;
t390 = -t294 * t39 + g(3) * t348 + (t24 - t300) * t262;
t389 = t294 * t38 + (-t24 + t394) * t264;
t346 = t243 * t273;
t347 = t243 * t270;
t328 = -g(1) * t346 - g(2) * t347 - g(3) * t242;
t388 = -t160 * t145 - t279 - t328;
t285 = t315 + t394;
t387 = t160 * t294 + t285;
t106 = -pkin(4) * t294 - qJ(5) * t145;
t158 = t214 * t263 - t341;
t132 = t158 - t378;
t159 = -t265 * t214 - t196;
t133 = t159 - t377;
t244 = pkin(2) * t265 + pkin(3);
t381 = pkin(2) * t263;
t338 = t268 * t244 + t382 * t381;
t360 = t338 * qJD(4) + t382 * t132 - t268 * t133;
t384 = g(1) * t270 - g(2) * t273;
t385 = t384 * t242;
t372 = g(3) * t272;
t371 = t264 * pkin(5);
t252 = t264 * pkin(9);
t370 = t269 * pkin(2);
t6 = t264 * t22 + t262 * t27;
t4 = t6 * t264;
t134 = -t186 * t263 + t265 * t187;
t193 = t203 * qJD(2);
t114 = -pkin(8) * t193 + t134;
t135 = t265 * t186 + t263 * t187;
t115 = -pkin(8) * t190 + t135;
t163 = -t265 * t228 - t229 * t263;
t136 = -pkin(8) * t205 + t163;
t164 = -t263 * t228 + t265 * t229;
t137 = pkin(8) * t203 + t164;
t295 = t382 * t136 - t268 * t137;
t44 = t295 * qJD(4) + t268 * t114 + t382 * t115;
t293 = t382 * t203 - t268 * t205;
t111 = t293 * qJD(4) - t268 * t190 + t382 * t193;
t156 = t268 * t203 + t382 * t205;
t112 = t156 * qJD(4) + t382 * t190 + t268 * t193;
t248 = t269 * t367;
t166 = pkin(3) * t190 + t248;
t48 = pkin(4) * t112 - qJ(5) * t111 - qJD(5) * t156 + t166;
t15 = t262 * t48 + t264 * t44;
t92 = t268 * t132 + t382 * t133;
t165 = pkin(3) * t191 + qJD(1) * t370;
t93 = t106 + t165;
t43 = t262 * t93 + t264 * t92;
t101 = -pkin(4) * t293 - qJ(5) * t156 + t174;
t103 = t268 * t136 + t382 * t137;
t54 = t262 * t101 + t264 * t103;
t364 = t262 * t278;
t361 = t360 - t416;
t50 = t262 * t106 + t264 * t79;
t359 = qJ(5) * t264;
t358 = t111 * t262;
t353 = t145 * t264;
t352 = t156 * t262;
t351 = t156 * t264;
t345 = t249 * t270;
t344 = t249 * t273;
t343 = t250 * t270;
t342 = t250 * t273;
t340 = -qJD(5) + t70;
t260 = t269 ^ 2;
t336 = -t272 ^ 2 + t260;
t331 = t269 * qJDD(1);
t330 = t272 * qJDD(1);
t236 = t268 * t381;
t2 = -pkin(5) * t278 - pkin(9) * t64 + t5;
t3 = -pkin(9) * t63 + t6;
t326 = t271 * t2 - t267 * t3;
t322 = -pkin(4) * t242 - pkin(3) * sin(t259) - t370;
t14 = -t262 * t44 + t264 * t48;
t42 = -t262 * t92 + t264 * t93;
t53 = t264 * t101 - t103 * t262;
t49 = t264 * t106 - t262 * t79;
t316 = t4 + t328;
t311 = t2 * t267 + t271 * t3;
t310 = -t262 * t6 - t264 * t5;
t309 = -t262 * t5 + t4;
t308 = t382 * t244 - t236;
t307 = t262 * t38 - t264 * t39;
t34 = -pkin(5) * t293 - pkin(9) * t351 + t53;
t40 = -pkin(9) * t352 + t54;
t305 = -t267 * t40 + t271 * t34;
t304 = t267 * t34 + t271 * t40;
t185 = -pkin(4) - t308;
t302 = -qJD(4) * t236 + t244 * t325;
t301 = pkin(1) + t395;
t299 = t384 * t243;
t298 = -0.2e1 * pkin(1) * t332 - pkin(7) * qJDD(2);
t184 = qJ(5) + t338;
t161 = (-pkin(9) - t184) * t262;
t172 = qJD(5) + t302;
t297 = -qJD(6) * t161 - t172 * t264 - t415 + t43;
t162 = t184 * t264 + t252;
t296 = -pkin(9) * t353 + qJD(6) * t162 + t172 * t262 - t406 + t42;
t226 = (-pkin(9) - qJ(5)) * t262;
t292 = -qJD(5) * t264 - qJD(6) * t226 - t415 + t50;
t227 = t252 + t359;
t291 = qJD(5) * t262 + qJD(6) * t227 - t145 * t252 - t406 + t49;
t288 = t184 * t278 + (t172 - t70) * t145;
t287 = -qJDD(1) * t246 + t329;
t274 = qJD(2) ^ 2;
t284 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t274 + t384;
t275 = qJD(1) ^ 2;
t283 = pkin(1) * t275 - pkin(7) * qJDD(1) + t313;
t282 = t111 * t70 + t156 * t24 - t313;
t45 = t103 * qJD(4) - t382 * t114 + t268 * t115;
t256 = -pkin(8) - t369;
t245 = -pkin(4) - t371;
t220 = qJ(5) * t346;
t219 = qJ(5) * t347;
t171 = t185 - t371;
t170 = t243 * t342 + t345;
t169 = -t243 * t344 + t343;
t168 = -t243 * t343 + t344;
t167 = t243 * t345 + t342;
t108 = t204 * t156;
t107 = t206 * t156;
t65 = pkin(5) * t352 - t295;
t56 = t80 + t416;
t37 = t111 * t206 + t333 * t351 - t334 * t352;
t36 = -t111 * t204 - t156 * t195;
t30 = pkin(5) * t358 + t45;
t10 = -pkin(9) * t358 + t15;
t9 = pkin(5) * t112 - t111 * t252 + t14;
t1 = [qJDD(1), t384, t313, qJDD(1) * t260 + 0.2e1 * t269 * t323, 0.2e1 * t269 * t330 - 0.2e1 * t332 * t336, qJDD(2) * t269 + t272 * t274, qJDD(2) * t272 - t269 * t274, 0, t269 * t298 + t272 * t284, -t269 * t284 + t272 * t298, -t109 * t205 + t110 * t203 - t134 * t191 + t135 * t189 - t152 * t193 - t153 * t190 - t163 * t154 + t164 * t276 - t313, t110 * t164 + t153 * t135 + t109 * t163 + t152 * t134 - t287 * t246 + t218 * t248 - g(1) * (-t246 * t270 + t273 * t369) - g(2) * (t246 * t273 + t270 * t369) -t111 * t294 + t156 * t75, t111 * t145 + t112 * t294 + t156 * t278 + t293 * t75, t111 * t258 + t156 * t255, -t112 * t258 + t255 * t293, 0, t112 * t160 - t125 * t293 - t145 * t166 - t174 * t278 + t255 * t295 - t258 * t45 + t299, -t103 * t255 + t111 * t160 + t125 * t156 - t166 * t294 + t174 * t75 - t258 * t44 - t385, t38 * t112 + t45 * t128 - t14 * t145 + t262 * t282 + t264 * t299 - t278 * t53 - t293 * t5 - t295 * t63, -t39 * t112 + t45 * t130 + t145 * t15 + t264 * t282 + t278 * t54 + t293 * t6 - t295 * t64 - t348 * t384, -t128 * t15 - t130 * t14 - t53 * t64 - t54 * t63 + t385 + t310 * t156 + (-t262 * t39 - t264 * t38) * t111, -t24 * t295 + t38 * t14 + t39 * t15 + t70 * t45 + t5 * t53 + t6 * t54 + (g(1) * t256 - g(2) * t301) * t273 + (g(1) * t301 + g(2) * t256) * t270, -t108 * t19 - t303 * t36, -t107 * t19 + t108 * t20 + t303 * t37 - t36 * t88, -t108 * t74 - t112 * t303 + t139 * t36 - t19 * t293, -t107 * t74 - t112 * t88 - t139 * t37 + t20 * t293, t112 * t139 - t293 * t74 (-t10 * t267 + t271 * t9) * t139 + t305 * t74 - t326 * t293 - t306 * t112 + t30 * t88 + t65 * t20 + t13 * t107 + t55 * t37 - g(1) * t168 - g(2) * t170 + (-t139 * t304 + t293 * t8) * qJD(6) -(t10 * t271 + t267 * t9) * t139 - t304 * t74 + t311 * t293 - t8 * t112 - t30 * t303 + t65 * t19 - t13 * t108 + t55 * t36 - g(1) * t167 - g(2) * t169 + (-t139 * t305 - t293 * t306) * qJD(6); 0, 0, 0, -t269 * t275 * t272, t336 * t275, t331, t330, qJDD(2), t269 * t283 - t372, g(3) * t269 + t272 * t283 (t153 + t158) * t191 + (-t159 + t152) * t189 + (-t265 * t154 + ((-t323 - t331) * t263 + (-t324 + t330) * t265) * t263) * pkin(2), -t152 * t158 - t153 * t159 + (-t372 + t109 * t265 + t110 * t263 + (-qJD(1) * t218 + t313) * t269) * pkin(2), t397, t393, t403, t414, t255, t145 * t165 + t255 * t308 - t258 * t360 + t387, -t338 * t255 + t165 * t294 + (-t302 + t92) * t258 + t388, t128 * t360 + t145 * t42 + t185 * t63 + t262 * t288 + t389, t130 * t360 - t145 * t43 + t185 * t64 + t264 * t288 + t390, t128 * t43 + t130 * t42 + (-t128 * t172 - t184 * t63 + t404) * t264 + (t130 * t172 + t184 * t64 + t392) * t262 + t316, t24 * t185 - t39 * t43 - t38 * t42 - g(1) * (t273 * t322 + t220) - g(2) * (t270 * t322 + t219) - g(3) * t395 + t360 * t70 + t309 * t184 - t307 * t172, t409, t407, t408, t413, t402 (t161 * t271 - t162 * t267) * t74 + t171 * t20 + t361 * t88 + (t267 * t297 - t271 * t296) * t139 + t411 -(t161 * t267 + t162 * t271) * t74 + t171 * t19 - t361 * t303 + (t267 * t296 + t271 * t297) * t139 + t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189 ^ 2 - t191 ^ 2, t152 * t191 - t153 * t189 + t287 - t384, 0, 0, 0, 0, 0, -t278 - t357, t75 + t355, t128 * t294 - t145 * t400 - t264 * t278, t130 * t294 - t145 * t353 + t364, -t262 * t63 - t264 * t64 + (t128 * t264 - t130 * t262) * t145, t145 * t307 + t294 * t70 - t310 - t384, 0, 0, 0, 0, 0, t314 + t366, -t365 + t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t397, t393, t403, t414, t255, t258 * t80 + t387, t79 * t258 + t388, qJ(5) * t364 - pkin(4) * t63 - t128 * t80 - (t262 * t340 - t49) * t145 + t389, t278 * t359 - pkin(4) * t64 - t130 * t80 - (t264 * t340 + t50) * t145 + t390, t128 * t50 + t130 * t49 + (-qJ(5) * t63 - qJD(5) * t128 + t404) * t264 + (qJ(5) * t64 + qJD(5) * t130 + t392) * t262 + t316, -t24 * pkin(4) - t39 * t50 - t38 * t49 - t70 * t80 - g(1) * (-pkin(4) * t349 + t220) - g(2) * (-pkin(4) * t350 + t219) - g(3) * t339 - t307 * qJD(5) + t309 * qJ(5), t409, t407, t408, t413, t402 (t226 * t271 - t227 * t267) * t74 + t245 * t20 - t56 * t88 + (t267 * t292 - t271 * t291) * t139 + t411 -(t226 * t267 + t227 * t271) * t74 + t245 * t19 + t56 * t303 + (t267 * t291 + t271 * t292) * t139 + t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130 * t145 + t63, t128 * t145 + t64, -t128 ^ 2 - t130 ^ 2, t128 * t39 + t130 * t38 - t285 + t383, 0, 0, 0, 0, 0, t20 - t401, t19 - t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303 * t88, t303 ^ 2 - t88 ^ 2, t19 + t405, -t20 - t401, t74, -g(1) * t169 + g(2) * t167 + t242 * t374 + t55 * t303 - t412 * t8 + t326, g(1) * t170 - g(2) * t168 + t242 * t373 + t306 * t412 + t55 * t88 - t311;];
tau_reg  = t1;
