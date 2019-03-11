% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x34]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:16
% EndTime: 2019-03-09 07:21:29
% DurationCPUTime: 4.95s
% Computational Cost: add. (6584->469), mult. (12963->623), div. (0->0), fcn. (9548->14), ass. (0->269)
t361 = cos(qJ(4));
t290 = t361 * qJD(4);
t233 = -pkin(1) - pkin(7);
t181 = t233 * qJD(1) + qJD(2);
t227 = sin(qJ(3));
t319 = qJD(1) * t227;
t136 = -pkin(8) * t319 + t181 * t227;
t226 = sin(qJ(4));
t128 = t226 * t136;
t231 = cos(qJ(3));
t318 = qJD(1) * t231;
t137 = -pkin(8) * t318 + t231 * t181;
t90 = t361 * t137 - t128;
t384 = pkin(3) * t290 - t90;
t296 = t361 * t231;
t145 = -qJD(1) * t296 + t226 * t319;
t216 = qJD(3) + qJD(4);
t225 = sin(qJ(5));
t230 = cos(qJ(5));
t115 = -t145 * t225 - t230 * t216;
t254 = -t226 * t231 - t361 * t227;
t146 = t254 * qJD(1);
t379 = qJD(5) - t146;
t390 = t115 * t379;
t261 = t145 * t230 - t216 * t225;
t389 = t261 * t379;
t232 = cos(qJ(1));
t214 = g(2) * t232;
t228 = sin(qJ(1));
t360 = g(1) * t228;
t372 = -t214 + t360;
t223 = qJ(3) + qJ(4);
t211 = cos(t223);
t177 = t211 * t214;
t209 = sin(t223);
t258 = -g(3) * t209 - t177;
t302 = t211 * t360;
t382 = t258 + t302;
t102 = -pkin(4) * t145 - pkin(9) * t146;
t93 = pkin(3) * t318 + t102;
t388 = -t225 * t384 - t230 * t93;
t224 = sin(qJ(6));
t229 = cos(qJ(6));
t262 = t115 * t224 + t229 * t261;
t62 = t229 * t115 - t224 * t261;
t387 = t262 * t62;
t154 = t224 * t225 - t229 * t230;
t366 = qJD(5) + qJD(6);
t346 = (t146 - t366) * t154;
t329 = t224 * t230;
t157 = t225 * t229 + t329;
t108 = t366 * t157;
t386 = -t157 * t146 + t108;
t215 = qJDD(3) + qJDD(4);
t130 = qJD(3) * pkin(3) + t137;
t315 = qJD(4) * t226;
t179 = t233 * qJDD(1) + qJDD(2);
t161 = t231 * t179;
t308 = qJD(1) * qJD(3);
t289 = t227 * t308;
t304 = t231 * qJDD(1);
t317 = qJD(3) * t227;
t92 = -t181 * t317 + qJDD(3) * pkin(3) + t161 + (t289 - t304) * pkin(8);
t316 = qJD(3) * t231;
t288 = t231 * t308;
t305 = t227 * qJDD(1);
t381 = t288 + t305;
t95 = -pkin(8) * t381 + t179 * t227 + t181 * t316;
t279 = -t130 * t315 - t136 * t290 - t226 * t95 + t361 * t92;
t27 = -pkin(4) * t215 - t279;
t385 = t27 + t302;
t383 = qJDD(2) - t372;
t217 = qJDD(1) * qJ(2);
t218 = qJD(1) * qJD(2);
t270 = g(1) * t232 + g(2) * t228;
t245 = -t270 + 0.2e1 * t218;
t380 = 0.2e1 * t217 + t245;
t378 = t262 ^ 2 - t62 ^ 2;
t311 = qJD(6) * t229;
t312 = qJD(6) * t224;
t313 = qJD(5) * t230;
t314 = qJD(5) * t225;
t284 = qJDD(1) * t361;
t373 = t146 * t216;
t71 = -t226 * t305 + t231 * t284 + t373;
t47 = t145 * t314 + t225 * t215 + t216 * t313 + t230 * t71;
t48 = -t261 * qJD(5) - t230 * t215 + t225 * t71;
t13 = -t115 * t311 - t224 * t48 + t229 * t47 + t261 * t312;
t135 = qJD(6) + t379;
t377 = t135 * t62 + t13;
t222 = qJ(5) + qJ(6);
t208 = sin(t222);
t210 = cos(t222);
t331 = t210 * t228;
t125 = t208 * t232 + t209 * t331;
t330 = t210 * t232;
t127 = -t208 * t228 + t209 * t330;
t197 = g(3) * t211;
t129 = t361 * t136;
t81 = t226 * t130 + t129;
t74 = pkin(9) * t216 + t81;
t171 = pkin(3) * t319 + qJD(1) * qJ(2);
t84 = -pkin(4) * t146 + pkin(9) * t145 + t171;
t42 = t225 * t84 + t230 * t74;
t30 = -pkin(10) * t115 + t42;
t28 = t30 * t312;
t80 = t361 * t130 - t128;
t73 = -t216 * pkin(4) - t80;
t50 = t115 * pkin(5) + t73;
t376 = g(1) * t125 - g(2) * t127 + t210 * t197 + t50 * t62 + t28;
t333 = t209 * t228;
t124 = -t208 * t333 + t330;
t332 = t209 * t232;
t126 = t208 * t332 + t331;
t241 = t130 * t290 - t136 * t315 + t226 * t92 + t361 * t95;
t26 = t215 * pkin(9) + t241;
t131 = pkin(3) * t381 + t217 + t218;
t244 = t361 * qJD(3) + t290;
t295 = t227 * t315;
t72 = -qJD(1) * t295 - t226 * t289 + t227 * t284 + t231 * (qJD(1) * t244 + qJDD(1) * t226);
t34 = pkin(4) * t72 - pkin(9) * t71 + t131;
t32 = t230 * t34;
t70 = qJDD(5) + t72;
t3 = pkin(5) * t70 - pkin(10) * t47 - qJD(5) * t42 - t225 * t26 + t32;
t256 = t225 * t34 + t230 * t26 + t84 * t313 - t314 * t74;
t4 = -pkin(10) * t48 + t256;
t298 = -t224 * t4 + t229 * t3;
t41 = -t225 * t74 + t230 * t84;
t29 = pkin(10) * t261 + t41;
t23 = pkin(5) * t379 + t29;
t349 = t229 * t30;
t9 = t224 * t23 + t349;
t375 = -g(1) * t124 - g(2) * t126 - qJD(6) * t9 + t208 * t197 + t50 * t262 + t298;
t243 = qJD(6) * t262 - t224 * t47 - t229 * t48;
t374 = -t135 * t262 + t243;
t89 = t226 * t137 + t129;
t275 = pkin(3) * t315 - t89;
t343 = pkin(1) * qJDD(1);
t371 = t343 - t383;
t370 = t224 * t314 + t225 * t312;
t336 = t146 * t225;
t369 = (t314 - t336) * pkin(5);
t368 = -qJD(6) * t230 - t313;
t367 = t225 * t93 - t230 * t384;
t362 = -pkin(9) - pkin(10);
t359 = g(3) * t230;
t358 = t230 * pkin(5);
t213 = t230 * pkin(10);
t357 = pkin(8) - t233;
t200 = pkin(3) * t226 + pkin(9);
t356 = -pkin(10) - t200;
t355 = t225 * t102 + t230 * t80;
t69 = qJDD(6) + t70;
t353 = t154 * t69;
t352 = t157 * t69;
t351 = t225 * t47;
t350 = t225 * t70;
t348 = t230 * t70;
t347 = t73 * t146;
t344 = t369 + t275;
t235 = qJD(1) ^ 2;
t342 = qJ(2) * t235;
t109 = -t226 * t316 - t227 * t244 - t231 * t315;
t341 = t109 * t230;
t110 = -t226 * t317 + t231 * t244 - t295;
t340 = t110 * t135;
t339 = t135 * t145;
t338 = t379 * t145;
t337 = t145 * t146;
t155 = t226 * t227 - t296;
t335 = t155 * t225;
t334 = t155 * t230;
t328 = t225 * t109;
t327 = t225 * t228;
t326 = t225 * t232;
t325 = t228 * t230;
t165 = t357 * t227;
t166 = t357 * t231;
t114 = -t361 * t165 - t226 * t166;
t105 = t230 * t114;
t324 = t230 * t232;
t195 = t227 * pkin(3) + qJ(2);
t323 = t109 * t216 - t155 * t215;
t103 = -pkin(4) * t254 + pkin(9) * t155 + t195;
t322 = t225 * t103 + t105;
t221 = t231 ^ 2;
t321 = t227 ^ 2 - t221;
t234 = qJD(3) ^ 2;
t320 = -t234 - t235;
t182 = pkin(3) * t316 + qJD(2);
t306 = qJDD(3) * t227;
t303 = pkin(10) * t336;
t300 = qJD(5) * pkin(9) * t379;
t66 = t73 * t314;
t297 = qJD(5) * t362;
t293 = t155 * t313;
t287 = qJD(6) * t23 + t4;
t285 = qJD(5) * t356;
t283 = -qJD(5) * t84 - t26;
t280 = t230 * t379;
t278 = -qJD(5) * t254 + qJD(1);
t201 = -t361 * pkin(3) - pkin(4);
t274 = -t81 + t369;
t273 = -t145 * pkin(5) - t146 * t213;
t272 = -t313 * t74 + t32;
t271 = -pkin(9) * t70 - t347;
t148 = t356 * t225;
t269 = -qJD(6) * t148 - t225 * t285 - t303 + t367;
t149 = t200 * t230 + t213;
t268 = qJD(6) * t149 - t230 * t285 + t273 - t388;
t172 = t362 * t225;
t267 = -qJD(6) * t172 - t225 * t297 - t303 + t355;
t173 = pkin(9) * t230 + t213;
t97 = t230 * t102;
t266 = qJD(6) * t173 - t225 * t80 - t230 * t297 + t273 + t97;
t265 = -t200 * t70 - t347;
t263 = -t110 * t216 + t215 * t254;
t260 = -t42 * t145 + t225 * t385 + t73 * t313;
t259 = t41 * t145 + t230 * t177 + t209 * t359 + t66;
t255 = t226 * t165 - t361 * t166;
t253 = t293 - t328;
t252 = t155 * t314 + t341;
t53 = pkin(4) * t110 - pkin(9) * t109 + t182;
t150 = t357 * t317;
t151 = qJD(3) * t166;
t58 = t255 * qJD(4) + t226 * t150 - t361 * t151;
t251 = t103 * t313 - t114 * t314 + t225 * t53 + t230 * t58;
t249 = t135 * t154;
t248 = 0.2e1 * qJ(2) * t308 + qJDD(3) * t233;
t246 = -t342 - t372;
t242 = t171 * t145 + t279 - t382;
t17 = pkin(5) * t48 + t27;
t8 = -t224 * t30 + t229 * t23;
t240 = t8 * t145 + t17 * t154 - t210 * t382 + t386 * t50;
t239 = -t233 * t234 + t380;
t59 = t114 * qJD(4) - t361 * t150 - t226 * t151;
t238 = -t9 * t145 + t17 * t157 + t208 * t382 + t346 * t50;
t236 = g(1) * t333 - g(2) * t332 - t171 * t146 + t197 - t241;
t207 = qJDD(3) * t231;
t202 = -pkin(4) - t358;
t170 = t201 - t358;
t141 = t209 * t324 - t327;
t140 = t209 * t326 + t325;
t139 = t209 * t325 + t326;
t138 = -t209 * t327 + t324;
t101 = t230 * t103;
t99 = t154 * t155;
t98 = t157 * t155;
t82 = -pkin(5) * t335 - t255;
t76 = t145 ^ 2 - t146 ^ 2;
t55 = -t145 * t216 - t72;
t54 = t71 - t373;
t52 = t230 * t53;
t49 = pkin(10) * t335 + t322;
t43 = -pkin(5) * t254 + pkin(10) * t334 - t114 * t225 + t101;
t36 = -pkin(5) * t253 + t59;
t22 = t109 * t329 + (-t334 * t366 + t328) * t229 + t370 * t155;
t21 = t108 * t155 - t109 * t154;
t20 = -t145 * t261 + t280 * t379 + t350;
t19 = -t225 * t379 ^ 2 - t115 * t145 + t348;
t18 = -t261 * t280 + t351;
t12 = -t135 * t386 - t145 * t62 - t353;
t11 = t346 * t135 - t145 * t262 + t352;
t10 = pkin(10) * t253 + t251;
t7 = -pkin(10) * t341 + pkin(5) * t110 - t225 * t58 + t52 + (-t105 + (-pkin(10) * t155 - t103) * t225) * qJD(5);
t6 = (t47 - t390) * t230 + (-t48 + t389) * t225;
t5 = t13 * t157 - t262 * t346;
t1 = -t13 * t154 + t157 * t243 + t262 * t386 - t346 * t62;
t2 = [qJDD(1), t372, t270, -0.2e1 * t343 + t383, t380, t371 * pkin(1) + (t245 + t217) * qJ(2), qJDD(1) * t221 - 0.2e1 * t227 * t288, -0.2e1 * t227 * t304 + 0.2e1 * t308 * t321, -t227 * t234 + t207, -t231 * t234 - t306, 0, t227 * t239 + t231 * t248, -t227 * t248 + t231 * t239, -t109 * t145 - t155 * t71, t109 * t146 + t110 * t145 + t155 * t72 + t254 * t71, t323, t263, 0, t110 * t171 - t131 * t254 - t146 * t182 + t195 * t72 - t209 * t270 + t215 * t255 - t216 * t59, t109 * t171 - t114 * t215 - t131 * t155 - t145 * t182 + t195 * t71 - t211 * t270 - t216 * t58, -t252 * t261 - t47 * t334 (-t115 * t230 + t225 * t261) * t109 + (t351 + t230 * t48 + (-t115 * t225 - t230 * t261) * qJD(5)) * t155, -t110 * t261 + t252 * t379 - t254 * t47 - t70 * t334, -t110 * t115 + t253 * t379 + t254 * t48 + t70 * t335, t110 * t379 - t254 * t70 (-t114 * t313 + t52) * t379 + t101 * t70 - t272 * t254 + t41 * t110 + t59 * t115 - t255 * t48 - t73 * t293 - g(1) * t141 - g(2) * t139 + ((-qJD(5) * t103 - t58) * t379 - t114 * t70 - t283 * t254 - t27 * t155 + t73 * t109) * t225, -t251 * t379 - t322 * t70 + t256 * t254 - t42 * t110 - t59 * t261 - t255 * t47 + t73 * t341 + g(1) * t140 - g(2) * t138 + (-t27 * t230 + t66) * t155, t13 * t99 - t21 * t262, t13 * t98 - t21 * t62 + t22 * t262 + t243 * t99, -t110 * t262 - t13 * t254 + t135 * t21 + t69 * t99, -t110 * t62 - t135 * t22 - t243 * t254 + t69 * t98, -t254 * t69 + t340 (-t10 * t224 + t229 * t7) * t135 + (-t224 * t49 + t229 * t43) * t69 - t298 * t254 + t8 * t110 + t36 * t62 - t82 * t243 - t17 * t98 + t50 * t22 - g(1) * t127 - g(2) * t125 + ((-t224 * t43 - t229 * t49) * t135 + t9 * t254) * qJD(6), g(1) * t126 - g(2) * t124 - t9 * t110 + t82 * t13 - t28 * t254 + t17 * t99 + t50 * t21 - t36 * t262 + (-(-qJD(6) * t49 + t7) * t135 - t43 * t69 + t3 * t254) * t224 + (-(qJD(6) * t43 + t10) * t135 - t49 * t69 + t287 * t254) * t229; 0, 0, 0, qJDD(1), -t235, -t342 - t371, 0, 0, 0, 0, 0, t227 * t320 + t207, t231 * t320 - t306, 0, 0, 0, 0, 0, qJD(1) * t146 + t323, qJD(1) * t145 + t263, 0, 0, 0, 0, 0, t254 * t350 - t109 * t115 + t155 * t48 + (-t110 * t225 - t230 * t278) * t379, t254 * t348 + t109 * t261 + t155 * t47 + (-t110 * t230 + t225 * t278) * t379, 0, 0, 0, 0, 0, -t109 * t62 - t155 * t243 - t157 * t340 + qJD(1) * t249 - ((t229 * t368 + t370) * t135 - t352) * t254, t109 * t262 + t155 * t13 + t110 * t249 + t157 * t135 * qJD(1) - (-(t224 * t368 - t225 * t311 - t229 * t314) * t135 + t353) * t254; 0, 0, 0, 0, 0, 0, t231 * t235 * t227, -t321 * t235, t304, -t305, qJDD(3), g(3) * t227 + t231 * t246 + t161, g(3) * t231 + (-t179 - t246) * t227, t337, t76, t54, t55, t215, t89 * t216 + (t146 * t318 + t361 * t215 - t216 * t315) * pkin(3) + t242, t90 * t216 + (t145 * t318 - t215 * t226 - t216 * t290) * pkin(3) + t236, t18, t6, t20, t19, t338, t201 * t48 - t385 * t230 + t265 * t225 + t275 * t115 + (-t200 * t313 + t388) * t379 + t259, t201 * t47 + t265 * t230 + t258 * t225 - t275 * t261 + (t200 * t314 + t367) * t379 + t260, t5, t1, t11, t12, t339 (t148 * t229 - t149 * t224) * t69 - t170 * t243 + t344 * t62 + (t224 * t269 - t229 * t268) * t135 + t240 -(t148 * t224 + t149 * t229) * t69 + t170 * t13 - t344 * t262 + (t224 * t268 + t229 * t269) * t135 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, t76, t54, t55, t215, t216 * t81 + t242, t80 * t216 + t236, t18, t6, t20, t19, t338, -pkin(4) * t48 - t81 * t115 - t97 * t379 + (t379 * t80 + t271) * t225 + (-t385 - t300) * t230 + t259, -pkin(4) * t47 + t355 * t379 + t81 * t261 + t271 * t230 + (t258 + t300) * t225 + t260, t5, t1, t11, t12, t339 (t172 * t229 - t173 * t224) * t69 - t202 * t243 + t274 * t62 + (t224 * t267 - t229 * t266) * t135 + t240 -(t172 * t224 + t173 * t229) * t69 + t202 * t13 - t274 * t262 + (t224 * t266 + t229 * t267) * t135 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261 * t115, -t115 ^ 2 + t261 ^ 2, t47 + t390, -t48 - t389, t70, -g(1) * t138 - g(2) * t140 + t261 * t73 + t379 * t42 + (t283 + t197) * t225 + t272, g(1) * t139 - g(2) * t141 + t115 * t73 + t211 * t359 + t379 * t41 - t256, -t387, t378, t377, t374, t69 -(-t224 * t29 - t349) * t135 + (-t135 * t312 + t229 * t69 + t261 * t62) * pkin(5) + t375 (-t30 * t135 - t3) * t224 + (t29 * t135 - t287) * t229 + (-t135 * t311 - t224 * t69 - t261 * t262) * pkin(5) + t376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t387, t378, t377, t374, t69, t135 * t9 + t375, t135 * t8 - t224 * t3 - t229 * t287 + t376;];
tau_reg  = t2;
