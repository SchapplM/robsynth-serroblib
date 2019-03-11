% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:32
% EndTime: 2019-03-09 18:13:43
% DurationCPUTime: 6.03s
% Computational Cost: add. (7376->477), mult. (16330->579), div. (0->0), fcn. (12254->12), ass. (0->253)
t225 = qJDD(2) + qJDD(3);
t217 = -qJDD(5) + t225;
t232 = sin(qJ(6));
t237 = cos(qJ(6));
t239 = cos(qJ(2));
t370 = cos(qJ(3));
t311 = t370 * t239;
t289 = qJD(1) * t311;
t234 = sin(qJ(3));
t235 = sin(qJ(2));
t326 = qJD(1) * t235;
t310 = t234 * t326;
t149 = -t289 + t310;
t165 = t234 * t239 + t370 * t235;
t151 = t165 * qJD(1);
t233 = sin(qJ(5));
t238 = cos(qJ(5));
t323 = qJD(5) * t238;
t324 = qJD(5) * t233;
t226 = qJD(2) + qJD(3);
t334 = t234 * t235;
t282 = t226 * t334;
t302 = qJDD(1) * t370;
t317 = t239 * qJDD(1);
t290 = -t226 * t289 - t234 * t317 - t235 * t302;
t73 = qJD(1) * t282 + t290;
t112 = t226 * t165;
t318 = t235 * qJDD(1);
t280 = t234 * t318 - t239 * t302;
t74 = qJD(1) * t112 + t280;
t28 = t149 * t323 - t151 * t324 + t233 * t74 - t238 * t73;
t320 = qJD(5) - t226;
t321 = qJD(6) * t237;
t322 = qJD(6) * t232;
t377 = t149 * t233 + t238 * t151;
t15 = -t232 * t217 + t237 * t28 + t320 * t321 - t322 * t377;
t278 = -t232 * t320 - t237 * t377;
t16 = -qJD(6) * t278 + t237 * t217 + t232 * t28;
t379 = -t238 * t149 + t151 * t233;
t387 = qJD(6) + t379;
t360 = t278 * t387;
t75 = t232 * t377 - t237 * t320;
t361 = t75 * t387;
t399 = (t15 - t361) * t237 + (-t16 + t360) * t232;
t141 = t151 * pkin(9);
t371 = pkin(8) + pkin(7);
t180 = t371 * t235;
t170 = qJD(1) * t180;
t355 = qJD(2) * pkin(2);
t159 = -t170 + t355;
t181 = t371 * t239;
t172 = qJD(1) * t181;
t335 = t234 * t172;
t98 = t370 * t159 - t335;
t309 = qJD(4) - t98;
t389 = -t141 + t309;
t353 = t15 * t232;
t394 = t237 * t387;
t9 = t278 * t394 - t353;
t29 = qJD(5) * t377 - t233 * t73 - t238 * t74;
t27 = qJDD(6) + t29;
t350 = t232 * t27;
t396 = t278 * t377 + t387 * t394 + t350;
t158 = t370 * t172;
t109 = -t234 * t170 + t158;
t325 = qJD(3) * t234;
t288 = pkin(2) * t325 - t109;
t110 = -t370 * t170 - t335;
t308 = qJD(3) * t370;
t330 = pkin(2) * t308 + qJD(4) - t110;
t393 = pkin(5) * t377;
t224 = t239 * pkin(2);
t211 = t224 + pkin(1);
t179 = t211 * qJD(1);
t82 = t149 * pkin(3) - t151 * qJ(4) - t179;
t60 = -pkin(4) * t149 - t82;
t30 = pkin(5) * t379 - pkin(10) * t377 + t60;
t242 = -pkin(3) - pkin(4);
t58 = t242 * t226 + t389;
t219 = t226 * qJ(4);
t362 = t149 * pkin(9);
t99 = t234 * t159 + t158;
t72 = t99 + t362;
t64 = t219 + t72;
t37 = t233 * t58 + t238 * t64;
t34 = pkin(10) * t320 + t37;
t12 = t232 * t30 + t237 * t34;
t392 = t12 * t377;
t359 = t387 * t377;
t391 = -t141 + t330;
t390 = -t362 + t288;
t358 = t377 * t379;
t319 = qJD(1) * qJD(2);
t388 = t235 * t319 - t317;
t236 = sin(qJ(1));
t229 = qJ(2) + qJ(3);
t222 = cos(t229);
t337 = t222 * t236;
t221 = sin(t229);
t339 = t221 * t238;
t121 = t233 * t337 - t236 * t339;
t240 = cos(qJ(1));
t336 = t222 * t240;
t338 = t221 * t240;
t123 = t233 * t336 - t238 * t338;
t272 = t221 * t233 + t222 * t238;
t263 = g(1) * t123 + g(2) * t121 + g(3) * t272;
t386 = t377 ^ 2 - t379 ^ 2;
t11 = -t232 * t34 + t237 * t30;
t306 = t239 * t319;
t115 = qJDD(2) * pkin(2) + t371 * (-t306 - t318);
t116 = t371 * t388;
t292 = -t370 * t115 - t234 * t116 + t159 * t325 + t172 * t308;
t218 = t225 * pkin(3);
t376 = qJDD(4) - t218;
t43 = t292 + t376;
t23 = -pkin(4) * t225 + pkin(9) * t73 + t43;
t291 = -t234 * t115 + t370 * t116 - t159 * t308 + t172 * t325;
t214 = t225 * qJ(4);
t216 = t226 * qJD(4);
t378 = t214 + t216;
t40 = -t291 + t378;
t24 = pkin(9) * t74 + t40;
t277 = -t238 * t23 + t233 * t24 + t64 * t323 + t58 * t324;
t6 = pkin(5) * t217 + t277;
t313 = t11 * t377 + t6 * t237;
t124 = t272 * t240;
t266 = t233 * t23 + t238 * t24 + t58 * t323 - t324 * t64;
t122 = t272 * t236;
t283 = g(2) * t122 + g(3) * (-t222 * t233 + t339);
t245 = -g(1) * t124 - t379 * t60 + t266 - t283;
t246 = t377 * t60 - t263 + t277;
t385 = t320 * t379 + t28;
t384 = t320 * t377 - t29;
t348 = t237 * t27;
t383 = t377 * t75 + t348;
t328 = t238 * qJ(4) + t233 * t242;
t169 = -pkin(10) + t328;
t369 = pkin(4) * t151;
t96 = t151 * pkin(3) + t149 * qJ(4);
t70 = -t96 - t369;
t35 = -pkin(10) * t379 - t393 + t70;
t381 = (qJD(6) * t169 + t35) * t387;
t210 = -t370 * pkin(2) - pkin(3);
t204 = -pkin(4) + t210;
t206 = pkin(2) * t234 + qJ(4);
t331 = t233 * t204 + t238 * t206;
t126 = -pkin(10) + t331;
t316 = pkin(2) * t326;
t380 = (qJD(6) * t126 - t316 + t35) * t387;
t120 = -t234 * t180 + t370 * t181;
t329 = t222 * pkin(3) + t221 * qJ(4);
t375 = -(t387 * pkin(10) + t393) * t387 + t263;
t364 = g(2) * t240;
t285 = g(1) * t236 - t364;
t312 = qJD(2) * t371;
t171 = t235 * t312;
t173 = t239 * t312;
t61 = -t370 * t171 - t234 * t173 - t180 * t308 - t181 * t325;
t374 = -t120 * t225 - t221 * t285 - t61 * t226;
t164 = -t311 + t334;
t104 = t164 * t233 + t165 * t238;
t48 = pkin(9) * t112 + t61;
t111 = -qJD(2) * t311 - t239 * t308 + t282;
t62 = t120 * qJD(3) - t234 * t171 + t370 * t173;
t49 = t111 * pkin(9) + t62;
t119 = t370 * t180 + t234 * t181;
t87 = -t165 * pkin(9) + t119;
t88 = pkin(9) * t164 + t120;
t50 = t233 * t88 - t238 * t87;
t13 = -qJD(5) * t50 + t233 * t49 + t238 * t48;
t274 = t238 * t164 - t165 * t233;
t304 = -pkin(10) * t217 + qJD(6) * t30 + t266;
t36 = -t233 * t64 + t238 * t58;
t33 = -pkin(5) * t320 - t36;
t366 = g(1) * t240;
t97 = t164 * pkin(3) - t165 * qJ(4) - t211;
t78 = -pkin(4) * t164 - t97;
t38 = -pkin(5) * t274 - pkin(10) * t104 + t78;
t46 = qJD(5) * t274 - t238 * t111 + t233 * t112;
t51 = t233 * t87 + t238 * t88;
t373 = t274 * t304 + t6 * t104 - t51 * t27 - (qJD(6) * t38 + t13) * t387 + t33 * t46 + t366;
t372 = t151 ^ 2;
t368 = g(1) * t122;
t273 = t204 * t238 - t206 * t233;
t357 = qJD(5) * t273 + t390 * t233 + t391 * t238;
t356 = t331 * qJD(5) + t391 * t233 - t390 * t238;
t354 = t104 * t33;
t352 = t226 * t98;
t351 = t226 * t99;
t349 = t232 * t387;
t347 = t237 * t278;
t276 = -qJ(4) * t233 + t238 * t242;
t346 = qJD(5) * t276 - t233 * t72 + t389 * t238;
t345 = t328 * qJD(5) + t389 * t233 + t238 * t72;
t344 = qJD(6) * t34;
t230 = qJDD(1) * pkin(1);
t342 = t151 * t149;
t340 = t221 * t236;
t227 = t235 ^ 2;
t327 = -t239 ^ 2 + t227;
t315 = t235 * t355;
t143 = t388 * pkin(2) - t230;
t301 = t387 * t33;
t296 = t320 * t387;
t293 = t320 ^ 2;
t287 = -pkin(2) * t235 - pkin(3) * t221;
t286 = g(2) * t236 + t366;
t284 = t38 * t27 + t368;
t281 = -t344 - t364;
t271 = t263 * t232 - t392;
t86 = t316 + t96;
t31 = t74 * pkin(3) + t73 * qJ(4) - t151 * qJD(4) + t143;
t270 = t211 + t329;
t269 = -t321 * t387 - t350;
t268 = -t322 * t387 + t348;
t267 = -0.2e1 * pkin(1) * t319 - pkin(7) * qJDD(2);
t52 = t112 * pkin(3) + t111 * qJ(4) - t165 * qJD(4) + t315;
t262 = g(1) * t336 + g(2) * t337 + g(3) * t221 + t291;
t261 = t283 - t304;
t19 = -pkin(4) * t74 - t31;
t260 = g(1) * t338 + g(2) * t340 - g(3) * t222 - t292;
t259 = -pkin(10) * t27 + (t33 + t36) * t387;
t39 = -pkin(4) * t112 - t52;
t257 = g(1) * t337 - g(2) * t336 - t119 * t225 - t226 * t62;
t243 = qJD(2) ^ 2;
t256 = -pkin(7) * t243 + 0.2e1 * t230 + t285;
t244 = qJD(1) ^ 2;
t255 = pkin(1) * t244 - pkin(7) * qJDD(1) + t286;
t254 = -t82 * t149 - t262;
t253 = -t179 * t149 + t262;
t252 = -t126 * t27 - t357 * t387 - t301;
t251 = -t169 * t27 - t346 * t387 - t301;
t250 = t179 * t151 + t260;
t247 = t82 * t151 - t260 + t376;
t183 = qJ(4) * t336;
t182 = qJ(4) * t337;
t168 = pkin(5) - t276;
t125 = pkin(5) - t273;
t107 = t124 * t237 - t232 * t236;
t106 = -t124 * t232 - t236 * t237;
t90 = t219 + t99;
t85 = -pkin(3) * t226 + t309;
t79 = -t149 ^ 2 + t372;
t63 = -t86 - t369;
t54 = -t290 + (t149 - t310) * t226;
t47 = qJD(5) * t104 - t233 * t111 - t238 * t112;
t14 = qJD(5) * t51 + t233 * t48 - t238 * t49;
t10 = pkin(5) * t47 - pkin(10) * t46 + t39;
t7 = t349 * t387 - t383;
t3 = pkin(5) * t29 - pkin(10) * t28 + t19;
t2 = t237 * t3;
t1 = [qJDD(1), t285, t286, qJDD(1) * t227 + 0.2e1 * t235 * t306, 0.2e1 * t235 * t317 - 0.2e1 * t319 * t327, qJDD(2) * t235 + t239 * t243, qJDD(2) * t239 - t235 * t243, 0, t235 * t267 + t239 * t256, -t235 * t256 + t239 * t267, -t111 * t151 - t165 * t73, t111 * t149 - t112 * t151 + t164 * t73 - t165 * t74, -t111 * t226 + t165 * t225, -t112 * t226 - t164 * t225, 0, -t112 * t179 + t143 * t164 + t149 * t315 - t211 * t74 + t257, t179 * t111 + t143 * t165 + t151 * t315 + t211 * t73 + t374, t112 * t82 + t149 * t52 + t164 * t31 + t74 * t97 + t257, -t111 * t85 - t112 * t90 - t119 * t73 - t120 * t74 - t149 * t61 + t151 * t62 - t164 * t40 + t165 * t43 - t286, t82 * t111 - t52 * t151 - t31 * t165 + t97 * t73 - t374, t43 * t119 + t40 * t120 + t31 * t97 + t82 * t52 + t90 * t61 + t85 * t62 + (-g(1) * t371 - g(2) * t270) * t240 + (g(1) * t270 - g(2) * t371) * t236, t104 * t28 + t377 * t46, -t104 * t29 + t274 * t28 - t377 * t47 - t379 * t46, -t104 * t217 + t320 * t46, -t217 * t274 - t320 * t47, 0, -g(2) * t124 - t14 * t320 - t19 * t274 + t217 * t50 + t29 * t78 + t379 * t39 + t47 * t60 + t368, -g(1) * t121 + g(2) * t123 + t104 * t19 - t13 * t320 + t217 * t51 + t28 * t78 + t377 * t39 + t46 * t60, -t46 * t347 + (t15 * t237 + t278 * t322) * t104 (t232 * t278 - t237 * t75) * t46 + (-t353 - t16 * t237 + (t232 * t75 + t347) * qJD(6)) * t104, t104 * t268 - t15 * t274 - t278 * t47 + t394 * t46, t104 * t269 + t16 * t274 - t349 * t46 - t75 * t47, -t27 * t274 + t387 * t47, -g(2) * t107 - t2 * t274 + t11 * t47 + t14 * t75 + t50 * t16 + (t10 * t387 + (t274 * t34 - t387 * t51 + t354) * qJD(6) + t284) * t237 + t373 * t232, -g(2) * t106 - t12 * t47 - t14 * t278 + t50 * t15 + (-(-qJD(6) * t51 + t10) * t387 + (t3 - t344) * t274 - qJD(6) * t354 - t284) * t232 + t373 * t237; 0, 0, 0, -t235 * t244 * t239, t327 * t244, t318, t317, qJDD(2), -g(3) * t239 + t235 * t255, g(3) * t235 + t239 * t255, t342, t79, t54, -t280, t225, t109 * t226 + (-t149 * t326 + t370 * t225 - t226 * t325) * pkin(2) + t250, t110 * t226 + (-t151 * t326 - t225 * t234 - t226 * t308) * pkin(2) + t253, -t86 * t149 - t210 * t225 - t226 * t288 - t247, -t206 * t74 - t210 * t73 + (t288 + t90) * t151 + (t85 - t330) * t149, t86 * t151 + t206 * t225 + t226 * t330 + t254 + t378, t40 * t206 + t43 * t210 - t82 * t86 - g(1) * (t240 * t287 + t183) - g(2) * (t236 * t287 + t182) - g(3) * (t224 + t329) + t330 * t90 + t288 * t85, -t358, -t386, -t385, -t384, t217, -t217 * t273 - t320 * t356 - t379 * t63 + t246, t217 * t331 - t320 * t357 - t377 * t63 + t245, t9, -t399, -t396, t7, t359, t125 * t16 + t356 * t75 + (-t263 - t380) * t237 + t252 * t232 + t313, t125 * t15 - t356 * t278 + (-t6 + t380) * t232 + t252 * t237 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, t79, t54, -t280, t225, t250 + t351, t253 + t352, -t149 * t96 + t218 - t247 + t351, pkin(3) * t73 - t74 * qJ(4) + (t90 - t99) * t151 + (t85 - t309) * t149, t151 * t96 + 0.2e1 * t214 + 0.2e1 * t216 + t254 - t352, t40 * qJ(4) - t43 * pkin(3) - t82 * t96 - t85 * t99 - g(1) * (-pkin(3) * t338 + t183) - g(2) * (-pkin(3) * t340 + t182) - g(3) * t329 + t309 * t90, -t358, -t386, -t385, -t384, t217, -t217 * t276 - t320 * t345 - t379 * t70 + t246, t217 * t328 - t320 * t346 - t377 * t70 + t245, t9, -t399, -t396, t7, t359, t168 * t16 + t345 * t75 + (-t263 - t381) * t237 + t251 * t232 + t313, t168 * t15 - t345 * t278 + (-t6 + t381) * t232 + t251 * t237 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225 + t342, t54, -t226 ^ 2 - t372, -t226 * t90 + t247, 0, 0, 0, 0, 0, -t151 * t379 - t238 * t217 - t233 * t293, -t151 * t377 + t233 * t217 - t238 * t293, 0, 0, 0, 0, 0, -t151 * t394 + (-t232 * t296 - t16) * t238 + (t320 * t75 + t269) * t233, t151 * t349 + (-t237 * t296 - t15) * t238 + (-t278 * t320 - t268) * t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t358, t386, t385, t384, -t217, t320 * t37 - t246, t320 * t36 - t245, -t9, t399, t396, -t232 * t387 ^ 2 + t383, -t359, -pkin(5) * t16 + t259 * t232 + t237 * t375 - t37 * t75 - t313, -pkin(5) * t15 + t392 + t37 * t278 + t259 * t237 + (-t375 + t6) * t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278 * t75, t278 ^ 2 - t75 ^ 2, t15 + t361, -t16 - t360, t27, -g(1) * t106 + t12 * t387 + t232 * t261 + t237 * t281 + t278 * t33 + t2, g(1) * t107 + t11 * t387 + t33 * t75 + (-t281 - t3) * t232 + t261 * t237;];
tau_reg  = t1;
