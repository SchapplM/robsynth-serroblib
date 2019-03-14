% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:26
% EndTime: 2019-03-09 03:35:38
% DurationCPUTime: 7.04s
% Computational Cost: add. (12931->558), mult. (28669->678), div. (0->0), fcn. (21177->18), ass. (0->300)
t229 = sin(qJ(3));
t342 = sin(pkin(11));
t287 = t342 * t229;
t232 = cos(qJ(3));
t343 = cos(pkin(11));
t288 = t343 * t232;
t253 = t287 - t288;
t259 = t253 * qJDD(1);
t252 = -t229 * t343 - t232 * t342;
t379 = t252 * qJD(1);
t390 = t379 * qJD(3);
t110 = t259 - t390;
t275 = qJD(1) * t287;
t243 = -qJD(3) * t275 - qJDD(1) * t252;
t276 = qJD(1) * t288;
t111 = qJD(3) * t276 + t243;
t160 = -t275 + t276;
t228 = sin(qJ(5));
t375 = cos(qJ(5));
t293 = qJD(5) * t375;
t317 = qJD(5) * t228;
t257 = -t228 * t110 + t375 * t111 + t160 * t293 + t317 * t379;
t284 = t375 * t110 + t228 * t111;
t260 = -t228 * t160 + t375 * t379;
t391 = qJD(5) * t260;
t56 = t284 - t391;
t313 = qJD(1) * qJD(3);
t292 = t229 * t313;
t225 = cos(pkin(10));
t204 = -t225 * pkin(1) - pkin(2);
t217 = t232 * pkin(3);
t380 = t204 - t217;
t127 = pkin(3) * t292 + qJDD(1) * t380 + qJDD(4);
t371 = t110 * pkin(4);
t79 = t127 + t371;
t12 = t56 * pkin(5) - pkin(9) * t257 + t79;
t231 = cos(qJ(6));
t11 = t231 * t12;
t227 = sin(qJ(6));
t309 = qJD(3) + qJD(5);
t367 = t379 * pkin(8);
t314 = t229 * qJD(2);
t224 = sin(pkin(10));
t202 = pkin(1) * t224 + pkin(7);
t183 = t202 * qJD(1);
t281 = qJ(4) * qJD(1) + t183;
t382 = t232 * t281;
t129 = t314 + t382;
t120 = t342 * t129;
t216 = t232 * qJD(2);
t128 = -t229 * t281 + t216;
t123 = qJD(3) * pkin(3) + t128;
t80 = t343 * t123 - t120;
t67 = qJD(3) * pkin(4) + t367 + t80;
t368 = t160 * pkin(8);
t122 = t343 * t129;
t81 = t342 * t123 + t122;
t69 = t81 + t368;
t42 = t228 * t67 + t375 * t69;
t37 = pkin(9) * t309 + t42;
t101 = t375 * t160 + t228 * t379;
t158 = qJD(1) * t380 + qJD(4);
t114 = -t160 * pkin(4) + t158;
t57 = -pkin(5) * t101 + pkin(9) * t260 + t114;
t14 = t227 * t57 + t231 * t37;
t308 = qJDD(3) + qJDD(5);
t214 = t232 * qJDD(2);
t181 = t202 * qJDD(1);
t279 = -qJD(2) * qJD(3) - t181;
t312 = qJD(1) * qJD(4);
t78 = qJDD(3) * pkin(3) + t214 - qJD(3) * t382 + (-qJ(4) * qJDD(1) + t279 - t312) * t229;
t297 = -t229 * qJDD(2) + t232 * t279;
t318 = qJD(3) * t229;
t103 = -t183 * t318 - t297;
t310 = t232 * qJDD(1);
t82 = t232 * t312 + (-t292 + t310) * qJ(4) + t103;
t47 = -t342 * t82 + t343 * t78;
t39 = qJDD(3) * pkin(4) - t111 * pkin(8) + t47;
t48 = t342 * t78 + t343 * t82;
t40 = -pkin(8) * t110 + t48;
t9 = t228 * t39 + t67 * t293 - t69 * t317 + t375 * t40;
t7 = pkin(9) * t308 + t9;
t3 = -qJD(6) * t14 - t227 * t7 + t11;
t389 = qJD(6) - t101;
t369 = t14 * t389;
t403 = t3 + t369;
t13 = -t227 * t37 + t231 * t57;
t370 = t13 * t389;
t203 = pkin(3) * t343 + pkin(4);
t294 = t342 * pkin(3);
t156 = t228 * t203 + t375 * t294;
t83 = -t128 * t342 - t122;
t246 = t83 - t368;
t84 = t343 * t128 - t120;
t72 = t84 + t367;
t344 = qJD(5) * t156 - t228 * t72 + t375 * t246;
t402 = -t42 - t344;
t401 = t227 * t389;
t383 = t101 * t309;
t400 = t257 - t383;
t355 = t260 ^ 2;
t356 = t101 ^ 2;
t399 = t355 - t356;
t280 = t231 * t309;
t316 = qJD(6) * t227;
t29 = -qJD(6) * t280 - t227 * t308 - t231 * t257 - t260 * t316;
t26 = t29 * t227;
t315 = qJD(6) * t231;
t394 = t101 * t231;
t94 = t227 * t309 - t231 * t260;
t398 = -t26 + (t315 - t394) * t94;
t53 = qJDD(6) + t56;
t50 = t227 * t53;
t359 = t315 * t389 + t50;
t363 = t94 * t260;
t397 = -t389 * t394 + t359 + t363;
t41 = -t228 * t69 + t375 * t67;
t36 = -pkin(5) * t309 - t41;
t396 = t101 * t36;
t92 = -t227 * t260 - t280;
t365 = t92 * t260;
t395 = t389 * t260;
t354 = t101 * t260;
t220 = qJ(3) + pkin(11);
t215 = qJ(5) + t220;
t200 = sin(t215);
t221 = qJ(1) + pkin(10);
t212 = cos(t221);
t333 = t200 * t212;
t210 = sin(t221);
t334 = t200 * t210;
t393 = g(1) * t333 + g(2) * t334;
t326 = t260 * qJD(3);
t392 = -t326 - t284;
t32 = t36 * t316;
t388 = t13 * t260 + t231 * t393 + t32;
t201 = cos(t215);
t373 = g(3) * t201;
t289 = t228 * t40 - t375 * t39;
t10 = -qJD(5) * t42 - t289;
t8 = -pkin(5) * t308 - t10;
t384 = t8 * t227 + t36 * t315;
t387 = -t14 * t260 + t227 * t373 + t384;
t192 = g(3) * t200;
t331 = t201 * t212;
t332 = t201 * t210;
t296 = -g(1) * t331 - g(2) * t332 - t192;
t386 = -t114 * t101 - t296 - t9;
t385 = t114 * t260 - t289 - t373 + t393;
t65 = -pkin(5) * t260 - pkin(9) * t101;
t155 = t203 * t375 - t228 * t294;
t139 = t155 * qJD(5);
t44 = t228 * t246 + t375 * t72;
t345 = t139 - t44;
t51 = t231 * t53;
t262 = t316 * t389 - t51;
t341 = pkin(1) * qJDD(1);
t290 = -g(1) * t210 + g(2) * t212;
t381 = t290 * t200;
t322 = t201 * pkin(5) + t200 * pkin(9);
t340 = qJD(6) * t94;
t30 = t227 * t257 - t231 * t308 + t340;
t250 = t228 * t253;
t116 = -t252 * t375 - t250;
t245 = t375 * t253;
t248 = qJD(3) * t252;
t70 = -t228 * t248 + t245 * t309 - t252 * t317;
t350 = t231 * t70;
t378 = -t116 * t262 - t350 * t389;
t377 = t116 * t308 - t309 * t70;
t376 = t379 ^ 2;
t234 = qJD(3) ^ 2;
t372 = g(3) * t232;
t2 = qJD(6) * t13 + t227 * t12 + t231 * t7;
t1 = t2 * t231;
t230 = sin(qJ(1));
t366 = t230 * pkin(1);
t364 = t94 * t92;
t226 = -qJ(4) - pkin(7);
t362 = -t227 * t30 - t92 * t315;
t28 = t30 * t231;
t349 = t231 * t92;
t361 = -t116 * t28 + t70 * t349;
t115 = -t228 * t252 + t245;
t249 = qJD(3) * t253;
t71 = -qJD(5) * t250 - t228 * t249 - t248 * t375 - t252 * t293;
t360 = -t29 * t115 + t94 * t71;
t358 = -t101 * t70 - t116 * t56;
t353 = t227 * t92;
t352 = t227 * t94;
t27 = t231 * t29;
t348 = t231 * t94;
t346 = t110 * t252 - t160 * t249;
t339 = qJD(6) * t389;
t337 = t379 * t160;
t336 = t183 * t229;
t335 = t183 * t232;
t330 = t210 * t227;
t329 = t210 * t231;
t328 = t212 * t227;
t327 = t212 * t231;
t325 = qJ(4) + t202;
t283 = qJD(3) * t325;
t136 = t232 * qJD(4) - t229 * t283;
t137 = -t229 * qJD(4) - t232 * t283;
t89 = t343 * t136 + t342 * t137;
t165 = t325 * t229;
t166 = t325 * t232;
t109 = -t342 * t165 + t343 * t166;
t211 = cos(t220);
t321 = pkin(4) * t211 + t217;
t175 = pkin(2) + t321;
t233 = cos(qJ(1));
t218 = t233 * pkin(1);
t323 = t212 * t175 + t218;
t222 = t229 ^ 2;
t223 = t232 ^ 2;
t320 = t222 - t223;
t184 = qJD(1) * t204;
t319 = qJD(1) * t229;
t182 = qJDD(1) * t204;
t306 = -t27 + t362;
t305 = pkin(9) * t339;
t303 = t70 * t352;
t208 = pkin(3) * t318;
t152 = pkin(9) + t156;
t301 = t152 * t339;
t235 = qJD(1) ^ 2;
t299 = t229 * t235 * t232;
t298 = t94 * t316;
t295 = -t8 - t373;
t130 = pkin(3) * t319 - pkin(4) * t379;
t209 = sin(t220);
t176 = -pkin(3) * t229 - pkin(4) * t209;
t291 = -pkin(5) * t200 + t176;
t282 = t1 + t296;
t278 = t232 * t292;
t277 = -pkin(9) * t53 - t396;
t274 = g(1) * t212 + g(2) * t210;
t272 = g(1) * t230 - g(2) * t233;
t271 = -t115 * t30 - t71 * t92;
t219 = -pkin(8) + t226;
t270 = -t212 * t219 - t366;
t269 = -t115 * t257 + t260 * t71;
t268 = t13 * t231 + t14 * t227;
t267 = t13 * t227 - t14 * t231;
t108 = -t343 * t165 - t166 * t342;
t90 = pkin(8) * t252 + t108;
t91 = -pkin(8) * t253 + t109;
t60 = t228 * t90 + t375 * t91;
t126 = pkin(4) * t253 + t380;
t62 = t115 * pkin(5) - t116 * pkin(9) + t126;
t24 = t227 * t62 + t231 * t60;
t23 = -t227 * t60 + t231 * t62;
t147 = t314 + t335;
t265 = t101 * t401 - t262;
t264 = -qJD(6) * t57 + t192 - t7;
t263 = -t298 - t27;
t261 = t274 * t200;
t256 = -qJD(1) * t184 + t274;
t255 = -t139 * t389 - t152 * t53 - t396;
t88 = -t136 * t342 + t137 * t343;
t251 = 0.2e1 * qJD(3) * t184 - qJDD(3) * t202;
t244 = -t202 * t234 - 0.2e1 * t182 - t290;
t242 = -t116 * t359 + t401 * t70;
t241 = -qJD(6) * t268 - t3 * t227 + t1;
t131 = -pkin(4) * t248 + t208;
t104 = -t147 * qJD(3) - t229 * t181 + t214;
t146 = t216 - t336;
t240 = t103 * t232 - t104 * t229 + (-t146 * t232 - t147 * t229) * qJD(3);
t239 = t127 + t290;
t238 = pkin(8) * t249 + t88;
t237 = -t111 * t253 - t248 * t379;
t206 = t217 + pkin(2);
t179 = qJDD(3) * t232 - t229 * t234;
t178 = qJDD(3) * t229 + t232 * t234;
t170 = pkin(9) * t331;
t169 = pkin(9) * t332;
t157 = t160 ^ 2;
t151 = -pkin(5) - t155;
t135 = t201 * t327 + t330;
t134 = -t201 * t328 + t329;
t133 = -t201 * t329 + t328;
t132 = t201 * t330 + t327;
t113 = -qJDD(3) * t252 - t234 * t253;
t112 = -qJDD(3) * t253 + t234 * t252;
t74 = pkin(8) * t248 + t89;
t61 = t130 + t65;
t59 = t228 * t91 - t375 * t90;
t54 = -t115 * t308 - t309 * t71;
t31 = t71 * pkin(5) + t70 * pkin(9) + t131;
t20 = t227 * t65 + t231 * t41;
t19 = -t227 * t41 + t231 * t65;
t18 = qJD(5) * t60 + t228 * t74 - t238 * t375;
t17 = t228 * t238 + t293 * t90 - t317 * t91 + t375 * t74;
t16 = t227 * t61 + t231 * t44;
t15 = -t227 * t44 + t231 * t61;
t5 = -qJD(6) * t24 - t227 * t17 + t231 * t31;
t4 = qJD(6) * t23 + t231 * t17 + t227 * t31;
t6 = [0, 0, 0, 0, 0, qJDD(1), t272, g(1) * t233 + g(2) * t230, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t225 * t341 - t290, -0.2e1 * t224 * t341 + t274, 0 (t272 + (t224 ^ 2 + t225 ^ 2) * t341) * pkin(1), qJDD(1) * t222 + 0.2e1 * t278, 0.2e1 * t229 * t310 - 0.2e1 * t313 * t320, t178, qJDD(1) * t223 - 0.2e1 * t278, t179, 0, t229 * t251 + t232 * t244, -t229 * t244 + t232 * t251 (t222 + t223) * t181 + t240 - t274, t182 * t204 - g(1) * (-pkin(2) * t210 + pkin(7) * t212 - t366) - g(2) * (pkin(2) * t212 + pkin(7) * t210 + t218) + t240 * t202, -t111 * t252 + t249 * t379, t237 + t346, t113, t110 * t253 + t160 * t248, t112, 0, t108 * qJDD(3) + t380 * t110 + t127 * t253 - t160 * t208 - t290 * t211 + (-t158 * t252 + t88) * qJD(3), -t89 * qJD(3) - t109 * qJDD(3) + t111 * t380 - t127 * t252 - t158 * t249 - t208 * t379 + t209 * t290, -t108 * t111 - t109 * t110 + t89 * t160 + t88 * t379 + t47 * t252 + t248 * t81 - t274 + (qJD(3) * t80 - t48) * t253, t48 * t109 + t81 * t89 + t47 * t108 + t80 * t88 + t127 * t380 + t158 * t208 - g(1) * (-t206 * t210 - t212 * t226 - t366) - g(2) * (t206 * t212 - t210 * t226 + t218) t116 * t257 + t260 * t70, t269 + t358, t377, -t101 * t71 + t115 * t56, t54, 0, -t101 * t131 + t114 * t71 + t79 * t115 + t126 * t56 - t18 * t309 - t201 * t290 - t308 * t59, -t114 * t70 + t79 * t116 + t126 * t257 - t131 * t260 - t17 * t309 - t308 * t60 + t381, -t10 * t116 + t101 * t17 - t115 * t9 - t18 * t260 + t257 * t59 + t41 * t70 - t42 * t71 - t56 * t60 - t274, t9 * t60 + t42 * t17 - t10 * t59 - t41 * t18 + t79 * t126 + t114 * t131 - g(1) * (-t175 * t210 + t270) - g(2) * (-t210 * t219 + t323) t116 * t263 - t348 * t70, t303 + (t26 + (-t348 + t353) * qJD(6)) * t116 + t361, t360 + t378, -t116 * t362 - t353 * t70, t242 + t271, t115 * t53 + t389 * t71, -t36 * t227 * t70 - g(1) * t133 - g(2) * t135 + t3 * t115 + t116 * t384 + t13 * t71 + t18 * t92 + t23 * t53 + t59 * t30 + t389 * t5, -t36 * t350 - g(1) * t132 - g(2) * t134 - t2 * t115 - t14 * t71 + t18 * t94 - t24 * t53 - t59 * t29 - t4 * t389 + (t8 * t231 - t32) * t116, t23 * t29 - t24 * t30 - t4 * t92 - t5 * t94 + t268 * t70 - t381 + (qJD(6) * t267 - t2 * t227 - t3 * t231) * t116, t2 * t24 + t14 * t4 + t3 * t23 + t13 * t5 + t8 * t59 + t36 * t18 - g(1) * t270 - g(2) * (pkin(5) * t331 + pkin(9) * t333 + t323) + (-g(1) * (-t175 - t322) + g(2) * t219) * t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t179, -t178, 0, t103 * t229 + t104 * t232 - g(3) + (-t146 * t229 + t147 * t232) * qJD(3), 0, 0, 0, 0, 0, 0, t112, -t113, -t237 + t346, t248 * t80 - t249 * t81 - t252 * t48 - t253 * t47 - g(3), 0, 0, 0, 0, 0, 0, t54, -t377, -t269 + t358, -t10 * t115 + t116 * t9 - t41 * t71 - t42 * t70 - g(3), 0, 0, 0, 0, 0, 0, t242 - t271, t360 - t378, -t303 + (-t26 + (t348 + t353) * qJD(6)) * t116 + t361, t8 * t115 + t116 * t241 + t267 * t70 + t36 * t71 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t320 * t235, t229 * qJDD(1), t299, t310, qJDD(3), -t372 + t214 + (t147 - t335) * qJD(3) + (t256 + t279) * t229, g(3) * t229 + (t146 + t336) * qJD(3) + t256 * t232 + t297, 0, 0, t337, -t157 + t376 (t276 - t160) * qJD(3) + t243, -t337, -t259, qJDD(3), -g(3) * t211 - t83 * qJD(3) + t158 * t379 + t274 * t209 + (qJDD(3) * t343 + t160 * t319) * pkin(3) + t47, g(3) * t209 + t84 * qJD(3) - t158 * t160 + t274 * t211 + (-qJDD(3) * t342 + t319 * t379) * pkin(3) - t48 -(t81 + t83) * t379 + (t80 - t84) * t160 + (-t110 * t342 - t111 * t343) * pkin(3), -t80 * t83 - t81 * t84 + (t342 * t48 + t47 * t343 - t372 + (-qJD(1) * t158 + t274) * t229) * pkin(3), t354, t399, t400, -t354, t392, t308, -t344 * qJD(3) + qJD(5) * t402 + t130 * t101 + t155 * t308 + t385, t130 * t260 - t156 * t308 - t309 * t345 + t386, -t155 * t257 - t156 * t56 + t402 * t260 + (t345 + t41) * t101, -g(3) * t321 + t10 * t155 - t114 * t130 + t9 * t156 - t176 * t274 - t344 * t41 + t345 * t42, t398, -t298 + (t349 + t352) * t101 + t306, t397, t401 * t92 - t28, t265 - t365, t395, -t15 * t389 + t151 * t30 + t344 * t92 + (t295 - t301) * t231 + t255 * t227 + t388, -t151 * t29 + t16 * t389 + t344 * t94 + t255 * t231 + (-t261 + t301) * t227 + t387, t15 * t94 + t16 * t92 + (t101 * t13 - t139 * t92 - t152 * t30 + (t152 * t94 - t13) * qJD(6)) * t231 + (t101 * t14 + t139 * t94 - t152 * t29 - t3 + (t152 * t92 - t14) * qJD(6)) * t227 + t282, t8 * t151 - g(1) * (t212 * t291 + t170) - g(2) * (t210 * t291 + t169) - g(3) * (t321 + t322) + t344 * t36 + (t139 * t231 - t16) * t14 + (-t139 * t227 - t15) * t13 + t241 * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259 - 0.2e1 * t390 (t276 + t160) * qJD(3) + t243, -t157 - t376, -t81 * t160 - t379 * t80 + t239, 0, 0, 0, 0, 0, 0, t284 - t326 - 0.2e1 * t391, t257 + t383, -t355 - t356, -t42 * t101 - t260 * t41 + t239 + t371, 0, 0, 0, 0, 0, 0, t265 + t365, -t231 * t389 ^ 2 + t363 - t50 (t349 - t352) * t101 - t263 + t362, t36 * t260 + t403 * t231 + (t2 - t370) * t227 + t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, t399, t400, -t354, t392, t308, t42 * qJD(3) + t385, t309 * t41 + t386, 0, 0, t398, t394 * t92 - t401 * t94 + t306, t397, t353 * t389 - t28, -t389 * t401 - t365 + t51, t395, -pkin(5) * t30 - t19 * t389 - t42 * t92 + t277 * t227 + (t295 - t305) * t231 + t388, pkin(5) * t29 + t20 * t389 - t42 * t94 + t277 * t231 + (-t261 + t305) * t227 + t387, t19 * t94 + t20 * t92 + (-t370 + (-t30 + t340) * pkin(9)) * t231 + ((qJD(6) * t92 - t29) * pkin(9) - t403) * t227 + t282, -t8 * pkin(5) - t14 * t20 - t13 * t19 - t36 * t42 - g(1) * (-pkin(5) * t333 + t170) - g(2) * (-pkin(5) * t334 + t169) - g(3) * t322 + t241 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t364, -t92 ^ 2 + t94 ^ 2, t389 * t92 - t29, -t364, t389 * t94 - t30, t53, -g(1) * t134 + g(2) * t132 + t227 * t264 - t315 * t37 - t36 * t94 + t11 + t369, g(1) * t135 - g(2) * t133 + t370 + t36 * t92 + (qJD(6) * t37 - t12) * t227 + t264 * t231, 0, 0;];
tau_reg  = t6;