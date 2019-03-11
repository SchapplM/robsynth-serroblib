% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:59
% EndTime: 2019-03-09 04:30:13
% DurationCPUTime: 7.77s
% Computational Cost: add. (7986->590), mult. (16992->721), div. (0->0), fcn. (11391->14), ass. (0->290)
t237 = cos(qJ(3));
t337 = qJD(1) * t237;
t198 = -qJD(4) + t337;
t224 = g(3) * t237;
t234 = sin(qJ(3));
t226 = qJ(1) + pkin(9);
t216 = cos(t226);
t214 = sin(t226);
t394 = g(2) * t214;
t285 = g(1) * t216 + t394;
t252 = -t285 * t234 + t224;
t230 = sin(pkin(9));
t206 = pkin(1) * t230 + pkin(7);
t181 = t206 * qJD(1);
t333 = qJD(3) * t237;
t179 = t206 * qJDD(1);
t416 = -qJD(2) * qJD(3) - t179;
t85 = t237 * qJDD(2) - t181 * t333 + t234 * t416;
t76 = -qJDD(3) * pkin(3) - t85;
t248 = -t76 - t252;
t419 = pkin(8) * qJD(4) * t198 + t248;
t236 = cos(qJ(4));
t390 = qJ(5) + pkin(8);
t183 = t390 * t236;
t229 = sin(pkin(10));
t233 = sin(qJ(4));
t352 = t229 * t233;
t370 = cos(pkin(10));
t110 = t183 * t370 - t352 * t390;
t221 = t237 * qJDD(1);
t325 = qJD(1) * qJD(3);
t155 = t234 * t325 + qJDD(4) - t221;
t225 = qJ(4) + pkin(10);
t213 = sin(t225);
t418 = t110 * t155 - t213 * t252;
t327 = t236 * qJD(3);
t338 = qJD(1) * t234;
t162 = t233 * t338 - t327;
t335 = qJD(3) * t233;
t164 = t236 * t338 + t335;
t99 = t370 * t162 + t164 * t229;
t368 = t99 * t198;
t307 = t237 * t325;
t323 = t234 * qJDD(1);
t330 = qJD(4) * t234;
t410 = qJD(1) * t330 - qJDD(3);
t94 = -qJD(4) * t327 + (-t307 - t323) * t236 + t410 * t233;
t95 = t233 * (qJD(3) * (qJD(4) + t337) + t323) + t410 * t236;
t49 = -t229 * t95 - t370 * t94;
t28 = t49 - t368;
t264 = -t229 * t162 + t164 * t370;
t417 = t99 * t264;
t301 = t370 * t233;
t157 = t229 * t236 + t301;
t142 = t157 * qJD(4);
t343 = t157 * t337 - t142;
t300 = t370 * t236;
t331 = qJD(4) * t233;
t143 = qJD(4) * t300 - t229 * t331;
t315 = t233 * t337;
t342 = -t229 * t315 + t300 * t337 - t143;
t212 = pkin(4) * t236 + pkin(3);
t297 = t237 * t212 + t234 * t390;
t131 = t157 * t234;
t291 = t370 * t333;
t311 = t237 * t327;
t77 = -t143 * t234 - t229 * t311 - t233 * t291;
t279 = t49 * t131 - t264 * t77;
t350 = t233 * t234;
t132 = -t229 * t350 + t234 * t300;
t48 = -t229 * t94 + t370 * t95;
t312 = t233 * t333;
t78 = t142 * t234 + t229 * t312 - t236 * t291;
t389 = -t132 * t48 + t78 * t99;
t415 = t279 - t389;
t329 = qJD(4) * t236;
t310 = t234 * t329;
t414 = t310 + t312;
t334 = qJD(3) * t234;
t303 = -t237 * t49 + t264 * t334;
t372 = t132 * t155 + t78 * t198;
t413 = t372 - t303;
t304 = -t237 * t48 + t99 * t334;
t409 = -t155 * t131 - t77 * t198;
t412 = t409 - t304;
t411 = t409 + t304;
t382 = t264 ^ 2;
t302 = qJD(4) * t390;
t326 = t236 * qJD(5);
t138 = -t233 * t302 + t326;
t260 = -t233 * qJD(5) - t236 * t302;
t345 = t236 * t237;
t270 = pkin(4) * t234 - qJ(5) * t345;
t161 = t234 * t181;
t133 = qJD(2) * t237 - t161;
t288 = pkin(3) * t234 - pkin(8) * t237;
t166 = t288 * qJD(1);
t91 = -t233 * t133 + t236 * t166;
t65 = qJD(1) * t270 + t91;
t92 = t236 * t133 + t233 * t166;
t71 = -qJ(5) * t315 + t92;
t387 = (t260 - t65) * t370 + (-t138 + t71) * t229;
t134 = qJD(2) * t234 + t181 * t237;
t120 = qJD(3) * pkin(8) + t134;
t289 = pkin(3) * t237 + pkin(8) * t234;
t273 = -pkin(2) - t289;
t231 = cos(pkin(9));
t398 = pkin(1) * t231;
t149 = t273 - t398;
t123 = t149 * qJD(1);
t317 = -t234 * qJDD(2) + t237 * t416;
t84 = -t181 * t334 - t317;
t75 = qJDD(3) * pkin(8) + t84;
t167 = t288 * qJD(3);
t97 = qJD(1) * t167 + qJDD(1) * t149;
t23 = -t120 * t331 + t123 * t329 + t233 * t97 + t236 * t75;
t67 = -t120 * t233 + t236 * t123;
t408 = t198 * t67 + t23;
t369 = pkin(1) * qJDD(1);
t292 = -t134 + (-t315 + t331) * pkin(4);
t165 = t206 * t345;
t106 = t233 * t149 + t165;
t407 = t279 + t389;
t406 = t48 * pkin(5) - t49 * qJ(6) - t264 * qJD(6);
t261 = -t233 * t330 + t311;
t346 = t236 * t155;
t405 = -t198 * t261 + t234 * t346;
t156 = -t300 + t352;
t404 = t155 * t156 + t198 * t343 - t338 * t99;
t392 = g(3) * t234;
t403 = t237 * t285 + t392;
t119 = -qJD(3) * pkin(3) - t133;
t93 = pkin(4) * t162 + qJD(5) + t119;
t33 = pkin(5) * t99 - qJ(6) * t264 + t93;
t401 = -t33 * t264 - qJDD(6);
t400 = -t49 * t156 - t157 * t48 + t264 * t343 + t342 * t99;
t15 = -qJ(5) * t95 - qJD(5) * t162 + t23;
t68 = t120 * t236 + t123 * t233;
t87 = t236 * t97;
t24 = -qJD(4) * t68 - t233 * t75 + t87;
t9 = t155 * pkin(4) + t94 * qJ(5) - t164 * qJD(5) + t24;
t3 = -t229 * t15 + t370 * t9;
t4 = t370 * t15 + t229 * t9;
t397 = pkin(5) * t155;
t396 = g(1) * t214;
t393 = g(2) * t216;
t362 = t206 * t233;
t340 = t236 * t167 + t334 * t362;
t40 = -t234 * t326 + t270 * qJD(3) + (-t165 + (qJ(5) * t234 - t149) * t233) * qJD(4) + t340;
t341 = t149 * t329 + t233 * t167;
t347 = t234 * t236;
t45 = (-qJ(5) * qJD(4) - qJD(3) * t206) * t347 + (-qJD(5) * t234 + (-qJ(5) * qJD(3) - qJD(4) * t206) * t237) * t233 + t341;
t17 = t229 * t40 + t370 * t45;
t58 = -qJ(5) * t164 + t67;
t55 = -pkin(4) * t198 + t58;
t59 = -qJ(5) * t162 + t68;
t56 = t370 * t59;
t22 = t229 * t55 + t56;
t388 = -pkin(5) * t343 + qJ(6) * t342 - qJD(6) * t157 + t292;
t35 = t229 * t65 + t370 * t71;
t386 = pkin(5) * t338 - t387;
t30 = qJ(6) * t338 + t35;
t81 = t138 * t370 + t229 * t260;
t385 = t81 - t30;
t384 = t81 - t35;
t136 = t236 * t149;
t79 = -qJ(5) * t347 + t136 + (-pkin(4) - t362) * t237;
t90 = -qJ(5) * t350 + t106;
t44 = t229 * t79 + t370 * t90;
t383 = t99 ^ 2;
t380 = t198 * t264;
t379 = t229 * t59;
t25 = t229 * t58 + t56;
t376 = t25 * t264;
t374 = t68 * t198;
t371 = -t162 * t311 - t95 * t347;
t365 = t162 * t198;
t364 = t164 * t162;
t363 = t164 * t198;
t361 = t213 * t214;
t360 = t214 * t233;
t359 = t214 * t236;
t215 = cos(t225);
t358 = t215 * t234;
t357 = t215 * t237;
t356 = t216 * t233;
t355 = t216 * t234;
t354 = t216 * t236;
t353 = t216 * t237;
t351 = t390 * t237;
t349 = t233 * t237;
t26 = t370 * t58 - t379;
t344 = qJD(6) - t26;
t200 = pkin(4) * t350;
t139 = t234 * t206 + t200;
t227 = t234 ^ 2;
t228 = t237 ^ 2;
t339 = t227 - t228;
t208 = -pkin(2) - t398;
t182 = qJD(1) * t208;
t336 = qJD(3) * t162;
t332 = qJD(4) * t162;
t180 = qJDD(1) * t208;
t320 = t216 * t349;
t240 = qJD(1) ^ 2;
t318 = t237 * t240 * t234;
t108 = pkin(4) * t414 + t206 * t333;
t238 = cos(qJ(1));
t316 = t238 * pkin(1) + t216 * pkin(2) + t214 * pkin(7);
t314 = t164 * t333;
t313 = t198 * t335;
t309 = t198 * t338;
t235 = sin(qJ(1));
t305 = -t235 * pkin(1) + t216 * pkin(7);
t299 = t164 * t334 + t94 * t237;
t298 = -t94 + t332;
t295 = t164 * t310;
t294 = t234 * t307;
t189 = t234 * t396;
t293 = -g(2) * t355 + t189;
t113 = t215 * t216 + t237 * t361;
t115 = t213 * t353 - t214 * t215;
t287 = g(1) * t113 - g(2) * t115;
t114 = -t216 * t213 + t214 * t357;
t116 = t215 * t353 + t361;
t286 = g(1) * t114 - g(2) * t116;
t284 = g(1) * t235 - g(2) * t238;
t283 = t131 * t48 - t77 * t99;
t282 = pkin(5) * t215 + qJ(6) * t213;
t281 = -t382 - t383;
t280 = t382 - t383;
t277 = -t233 * t68 - t236 * t67;
t276 = t233 * t67 - t236 * t68;
t272 = t48 - t380;
t271 = t48 + t380;
t126 = t214 * t349 + t354;
t16 = -t229 * t45 + t370 * t40;
t21 = t370 * t55 - t379;
t43 = -t229 * t90 + t370 * t79;
t266 = -t233 * t155 + t198 * t329;
t265 = t48 * t156 - t343 * t99;
t263 = g(1) * t115 + g(2) * t113 + t213 * t392 + t3;
t262 = -qJD(1) * t182 + t285;
t109 = t229 * t183 + t301 * t390;
t259 = g(1) * t215 * t355 - g(3) * t357 - t109 * t155 + t358 * t394;
t258 = -pkin(8) * t155 - t119 * t198;
t257 = pkin(4) * t360 + t216 * t297 + t316;
t239 = qJD(3) ^ 2;
t256 = t206 * t239 + 0.2e1 * t180 + t393;
t255 = 0.2e1 * qJD(3) * t182 - qJDD(3) * t206;
t254 = -t49 - t368;
t250 = g(1) * t116 + g(2) * t114 + g(3) * t358 - t4;
t249 = -t25 * t198 + t263;
t47 = t95 * pkin(4) + qJDD(5) + t76;
t246 = pkin(4) * t356 + t305 + (-pkin(2) - t297) * t214;
t245 = qJD(4) * t277 + t23 * t236 - t24 * t233;
t244 = -t85 * t234 + t84 * t237 + (-t133 * t237 - t134 * t234) * qJD(3);
t243 = t109 * t49 - t110 * t48 - t81 * t99 - t403;
t242 = t47 + t252;
t207 = -pkin(4) * t370 - pkin(5);
t201 = pkin(4) * t229 + qJ(6);
t186 = pkin(4) * t359;
t173 = qJDD(3) * t237 - t234 * t239;
t172 = qJDD(3) * t234 + t237 * t239;
t148 = t155 * qJ(6);
t129 = t216 * t345 + t360;
t128 = -t320 + t359;
t127 = -t214 * t345 + t356;
t107 = -t155 * t237 - t198 * t334;
t105 = -t206 * t349 + t136;
t96 = pkin(5) * t156 - qJ(6) * t157 - t212;
t62 = pkin(5) * t131 - qJ(6) * t132 + t139;
t61 = -qJD(4) * t106 + t340;
t60 = (-t234 * t327 - t237 * t331) * t206 + t341;
t46 = pkin(4) * t164 + pkin(5) * t264 + qJ(6) * t99;
t41 = t237 * pkin(5) - t43;
t39 = -qJ(6) * t237 + t44;
t29 = t157 * t155 + t198 * t342 - t264 * t338;
t27 = -pkin(5) * t77 + qJ(6) * t78 - qJD(6) * t132 + t108;
t20 = -qJ(6) * t198 + t22;
t19 = t198 * pkin(5) + qJD(6) - t21;
t18 = t132 * t49 - t264 * t78;
t14 = -pkin(5) * t334 - t16;
t11 = qJ(6) * t334 - qJD(6) * t237 + t17;
t10 = t49 * t157 - t264 * t342;
t8 = t303 + t372;
t5 = t47 + t406;
t2 = qJDD(6) - t397 - t3;
t1 = -qJD(6) * t198 + t148 + t4;
t6 = [0, 0, 0, 0, 0, qJDD(1), t284, g(1) * t238 + g(2) * t235, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t231 * t369 - t393 + t396, -0.2e1 * t230 * t369 + t285, 0 (t284 + (t230 ^ 2 + t231 ^ 2) * t369) * pkin(1), qJDD(1) * t227 + 0.2e1 * t294, 0.2e1 * t221 * t234 - 0.2e1 * t325 * t339, t172, qJDD(1) * t228 - 0.2e1 * t294, t173, 0, t255 * t234 + (-t256 + t396) * t237, t234 * t256 + t237 * t255 - t189 (t227 + t228) * t179 + t244 - t285, t180 * t208 - g(1) * (-pkin(2) * t214 + t305) - g(2) * t316 + t244 * t206, t164 * t261 - t347 * t94, -t295 + (-t314 + (t94 + t332) * t234) * t233 + t371, t299 + t405, t162 * t414 + t95 * t350 (t95 + t313) * t237 + (t266 - t336) * t234, t107, -g(1) * t127 - g(2) * t129 + t105 * t155 - t61 * t198 + (-t24 + (t119 * t233 + t162 * t206) * qJD(3)) * t237 + (qJD(3) * t67 + t119 * t329 + t206 * t95 + t76 * t233) * t234, -g(1) * t126 - g(2) * t128 - t106 * t155 + t60 * t198 + (t23 + (t119 * t236 + t164 * t206) * qJD(3)) * t237 + (-qJD(3) * t68 - t119 * t331 - t206 * t94 + t76 * t236) * t234, t105 * t94 - t106 * t95 - t60 * t162 - t61 * t164 + t189 + t277 * t333 + (qJD(4) * t276 - t23 * t233 - t236 * t24 - t393) * t234, t23 * t106 + t68 * t60 + t24 * t105 + t67 * t61 - g(1) * t305 - g(2) * (t216 * t289 + t316) - t273 * t396 + (t119 * t333 + t234 * t76) * t206, t18, -t415, t8, t283, t412, t107, t108 * t99 + t131 * t47 + t139 * t48 + t155 * t43 - t16 * t198 + t21 * t334 - t237 * t3 - t77 * t93 + t286, t108 * t264 + t132 * t47 + t139 * t49 - t155 * t44 + t17 * t198 - t22 * t334 + t237 * t4 - t78 * t93 - t287, -t131 * t4 - t132 * t3 - t16 * t264 - t17 * t99 + t21 * t78 + t22 * t77 - t43 * t49 - t44 * t48 + t293, -g(1) * t246 - g(2) * t257 + t93 * t108 + t47 * t139 + t21 * t16 + t22 * t17 + t3 * t43 + t4 * t44, t18, t8, t415, t107, -t412, t283, t131 * t5 + t14 * t198 - t155 * t41 - t19 * t334 + t2 * t237 + t27 * t99 - t33 * t77 + t48 * t62 + t286, -t1 * t131 - t11 * t99 + t132 * t2 + t14 * t264 - t19 * t78 + t20 * t77 - t39 * t48 + t41 * t49 + t293, -t1 * t237 - t11 * t198 - t132 * t5 + t155 * t39 + t20 * t334 - t264 * t27 + t33 * t78 - t49 * t62 + t287, t1 * t39 + t20 * t11 + t5 * t62 + t33 * t27 + t2 * t41 + t19 * t14 - g(1) * (-t114 * pkin(5) - t113 * qJ(6) + t246) - g(2) * (pkin(5) * t116 + qJ(6) * t115 + t257); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t173, -t172, 0, t84 * t234 + t85 * t237 - g(3) + (-t133 * t234 + t134 * t237) * qJD(3), 0, 0, 0, 0, 0, 0 (-t95 + t313) * t237 + (t266 + t336) * t234, t299 - t405, t295 + (t234 * t298 + t314) * t233 + t371, -g(3) + (-qJD(3) * t276 - t76) * t237 + (qJD(3) * t119 + t245) * t234, 0, 0, 0, 0, 0, 0, t411, -t413, t407, -t131 * t3 + t132 * t4 + t21 * t77 - t22 * t78 - t237 * t47 + t334 * t93 - g(3), 0, 0, 0, 0, 0, 0, t411, t407, t413, t1 * t132 + t131 * t2 - t19 * t77 - t20 * t78 - t237 * t5 + t33 * t334 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, t339 * t240, t323, t318, t221, qJDD(3), t134 * qJD(3) + t234 * t262 - t224 + t85, t392 + (t133 + t161) * qJD(3) + t262 * t237 + t317, 0, 0, -t94 * t233 - t236 * t363 (-t94 + t365) * t236 + (-t95 + t363) * t233 (-t164 * t234 + t198 * t345) * qJD(1) - t266, -t233 * t365 - t95 * t236, t198 * t331 + t346 + (t162 * t234 - t198 * t349) * qJD(1), t309, -pkin(3) * t95 - t134 * t162 + t91 * t198 + t258 * t233 + t236 * t419 - t67 * t338, pkin(3) * t94 - t134 * t164 - t92 * t198 - t233 * t419 + t258 * t236 + t68 * t338, t92 * t162 + t91 * t164 + ((qJD(4) * t164 - t95) * pkin(8) + t408) * t236 + (pkin(8) * t298 - t24 + t374) * t233 - t403, -t119 * t134 - t67 * t91 - t68 * t92 + t248 * pkin(3) + (t245 - t403) * pkin(8), t10, t400, t29, t265, -t404, t309, t47 * t156 - t198 * t387 - t21 * t338 - t212 * t48 + t292 * t99 - t343 * t93 + t259, t47 * t157 + t198 * t384 - t212 * t49 + t22 * t338 + t264 * t292 - t342 * t93 - t418, -t4 * t156 - t3 * t157 + t21 * t342 + t22 * t343 - t264 * t387 + t35 * t99 + t243, -g(3) * t297 - t3 * t109 + t4 * t110 + t21 * t387 - t47 * t212 + t22 * t384 + t292 * t93 + t285 * (t212 * t234 - t351) t10, t29, -t400, t309, t404, t265, t5 * t156 + t19 * t338 + t198 * t386 - t33 * t343 + t388 * t99 + t96 * t48 + t259, -t1 * t156 + t2 * t157 - t19 * t342 + t20 * t343 + t264 * t386 + t30 * t99 + t243, -t5 * t157 - t198 * t385 - t20 * t338 - t264 * t388 + t33 * t342 - t96 * t49 + t418, t1 * t110 + t5 * t96 + t2 * t109 - g(3) * (t237 * t282 + t297) + t388 * t33 + t385 * t20 + t386 * t19 + t285 * (-t351 - (-t212 - t282) * t234); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t364, -t162 ^ 2 + t164 ^ 2, -t94 - t365, -t364, -t363 - t95, t155, -t120 * t329 - g(1) * t128 + g(2) * t126 - t119 * t164 - t374 + t87 + (-qJD(4) * t123 + t392 - t75) * t233, g(1) * t129 - g(2) * t127 + g(3) * t347 + t119 * t162 - t408, 0, 0, t417, t280, t28, -t417, -t271, t155, -t93 * t264 + (t155 * t370 - t164 * t99) * pkin(4) + t249, t93 * t99 - t26 * t198 + (-t155 * t229 - t164 * t264) * pkin(4) + t250, -t376 + t22 * t264 + (-t229 * t48 - t370 * t49) * pkin(4) + (-t21 + t26) * t99, -g(1) * t186 + t21 * t25 - t22 * t26 + (g(2) * t354 - t93 * t164 + t4 * t229 + t233 * t403 + t3 * t370) * pkin(4), t417, t28, -t280, t155, t271, -t417, -t46 * t99 + (pkin(5) - t207) * t155 + t249 + t401, t20 * t264 - t201 * t48 + t207 * t49 - t376 + (t19 - t344) * t99, -t33 * t99 + t46 * t264 + t201 * t155 + t148 + (-0.2e1 * qJD(6) + t26) * t198 - t250, t1 * t201 + t2 * t207 - t33 * t46 - t19 * t25 - g(1) * (-pkin(4) * t320 - pkin(5) * t115 + qJ(6) * t116 + t186) - g(2) * (-pkin(4) * t126 - t113 * pkin(5) + t114 * qJ(6)) - g(3) * (-t200 + (-pkin(5) * t213 + qJ(6) * t215) * t234) + t344 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, -t254, t281, t21 * t264 + t22 * t99 + t242, 0, 0, 0, 0, 0, 0, t272, t281, t254, -t19 * t264 + t20 * t99 + t242 + t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155 + t417, t28, -t198 ^ 2 - t382, t198 * t20 - t263 - t397 - t401;];
tau_reg  = t6;
