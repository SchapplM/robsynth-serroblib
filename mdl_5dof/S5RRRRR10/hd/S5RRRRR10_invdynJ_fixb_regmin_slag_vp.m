% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:36:02
% EndTime: 2019-12-31 22:36:21
% DurationCPUTime: 8.12s
% Computational Cost: add. (8296->502), mult. (20314->719), div. (0->0), fcn. (16565->14), ass. (0->260)
t230 = cos(qJ(5));
t319 = qJD(5) * t230;
t224 = cos(pkin(5));
t328 = qJD(1) * t224;
t208 = qJD(2) + t328;
t232 = cos(qJ(3));
t227 = sin(qJ(3));
t228 = sin(qJ(2));
t223 = sin(pkin(5));
t329 = qJD(1) * t223;
t302 = t228 * t329;
t279 = t227 * t302;
t143 = t208 * t232 - t279;
t144 = t208 * t227 + t232 * t302;
t226 = sin(qJ(4));
t231 = cos(qJ(4));
t92 = -t231 * t143 + t144 * t226;
t396 = t230 * t92;
t405 = -t319 - t396;
t233 = cos(qJ(2));
t317 = qJD(1) * qJD(2);
t296 = t233 * t317;
t315 = qJDD(1) * t228;
t404 = t296 + t315;
t327 = qJD(1) * t233;
t301 = t223 * t327;
t280 = t227 * t301;
t370 = pkin(8) + pkin(9);
t303 = qJD(3) * t370;
t311 = pkin(1) * t328;
t159 = -pkin(7) * t302 + t233 * t311;
t272 = pkin(2) * t228 - pkin(8) * t233;
t160 = t272 * t329;
t332 = t232 * t159 + t227 * t160;
t403 = pkin(9) * t280 - t227 * t303 - t332;
t148 = t232 * t160;
t339 = t232 * t233;
t402 = t232 * t303 - t159 * t227 + t148 + (pkin(3) * t228 - pkin(9) * t339) * t329;
t314 = qJDD(1) * t233;
t206 = t223 * t314;
t297 = t228 * t317;
t278 = t223 * t297;
t157 = qJDD(3) - t206 + t278;
t151 = qJDD(4) + t157;
t192 = -qJD(3) + t301;
t179 = -qJD(4) + t192;
t225 = sin(qJ(5));
t260 = t143 * t226 + t231 * t144;
t320 = qJD(5) * t225;
t321 = qJD(4) * t231;
t322 = qJD(4) * t226;
t316 = qJDD(1) * t224;
t207 = qJDD(2) + t316;
t323 = qJD(3) * t232;
t393 = t404 * t223;
t81 = -qJD(3) * t279 + t227 * t207 + t208 * t323 + t393 * t232;
t325 = qJD(2) * t233;
t298 = t227 * t325;
t324 = qJD(3) * t227;
t82 = -t232 * t207 + t208 * t324 + t223 * (qJD(1) * (t228 * t323 + t298) + t227 * t315);
t35 = t143 * t321 - t144 * t322 - t226 * t82 + t231 * t81;
t15 = t225 * t151 - t179 * t319 + t230 * t35 - t260 * t320;
t13 = t15 * t225;
t261 = t179 * t225 - t230 * t260;
t401 = t405 * t261 + t13;
t36 = qJD(4) * t260 + t226 * t81 + t231 * t82;
t32 = qJDD(5) + t36;
t29 = t225 * t32;
t336 = -qJD(5) - t92;
t400 = t260 * t261 + t405 * t336 + t29;
t30 = t230 * t32;
t394 = t225 * t336;
t66 = t230 * t179 + t225 * t260;
t399 = t260 * t66 - t336 * t394 + t30;
t16 = -qJD(5) * t261 - t230 * t151 + t225 * t35;
t398 = t15 * t230 - t225 * t16 - t261 * t394 + t405 * t66;
t162 = pkin(7) * t301 + t228 * t311;
t127 = pkin(8) * t208 + t162;
t258 = -pkin(2) * t233 - pkin(8) * t228 - pkin(1);
t156 = t258 * t223;
t139 = qJD(1) * t156;
t79 = t127 * t232 + t139 * t227;
t59 = pkin(9) * t143 + t79;
t360 = t226 * t59;
t78 = -t127 * t227 + t232 * t139;
t58 = -pkin(9) * t144 + t78;
t55 = -pkin(3) * t192 + t58;
t27 = t231 * t55 - t360;
t25 = pkin(4) * t179 - t27;
t397 = t25 * t92;
t234 = cos(qJ(1));
t337 = t233 * t234;
t229 = sin(qJ(1));
t342 = t228 * t229;
t171 = -t224 * t342 + t337;
t222 = qJ(3) + qJ(4);
t217 = sin(t222);
t218 = cos(t222);
t347 = t223 * t229;
t114 = -t171 * t217 + t218 * t347;
t344 = t223 * t234;
t348 = t223 * t228;
t340 = t229 * t233;
t341 = t228 * t234;
t169 = t224 * t341 + t340;
t352 = t169 * t217;
t252 = -g(3) * (-t217 * t348 + t218 * t224) - g(2) * (-t218 * t344 - t352) - g(1) * t114;
t310 = pkin(1) * qJD(2) * t224;
t282 = qJD(1) * t310;
t308 = pkin(1) * t316;
t304 = -pkin(7) * t206 - t228 * t308 - t233 * t282;
t247 = -pkin(7) * t278 - t304;
t100 = pkin(8) * t207 + t247;
t255 = t272 * qJD(2);
t102 = (qJD(1) * t255 + qJDD(1) * t258) * t223;
t241 = -qJD(3) * t79 - t227 * t100 + t232 * t102;
t20 = pkin(3) * t157 - pkin(9) * t81 + t241;
t254 = -t232 * t100 - t227 * t102 + t127 * t324 - t139 * t323;
t22 = -pkin(9) * t82 - t254;
t357 = t231 * t59;
t28 = t226 * t55 + t357;
t243 = -qJD(4) * t28 + t231 * t20 - t22 * t226;
t5 = -pkin(4) * t151 - t243;
t249 = t252 - t5;
t395 = t260 * t92;
t173 = t226 * t227 - t231 * t232;
t380 = qJD(3) + qJD(4);
t117 = t380 * t173;
t129 = t173 * t301;
t335 = -t117 + t129;
t174 = t226 * t232 + t227 * t231;
t334 = (-t301 + t380) * t174;
t216 = -pkin(3) * t232 - pkin(2);
t110 = pkin(4) * t173 - pkin(10) * t174 + t216;
t193 = t370 * t227;
t194 = t370 * t232;
t132 = -t193 * t226 + t194 * t231;
t168 = -t224 * t337 + t342;
t170 = t224 * t340 + t341;
t345 = t223 * t233;
t248 = -g(1) * t170 - g(2) * t168 + g(3) * t345;
t245 = t248 * t218;
t275 = -t162 + (-t280 + t324) * pkin(3);
t392 = -t25 * qJD(5) * t174 - t110 * t32 - (-t334 * pkin(4) + t335 * pkin(10) + qJD(5) * t132 - t275) * t336 + t245;
t391 = t260 ^ 2 - t92 ^ 2;
t53 = pkin(4) * t260 + pkin(10) * t92;
t390 = -t179 * t92 + t35;
t112 = t169 * t218 - t217 * t344;
t115 = t171 * t218 + t217 * t347;
t153 = t217 * t224 + t218 * t348;
t291 = -t226 * t20 - t231 * t22 - t55 * t321 + t59 * t322;
t126 = -pkin(2) * t208 - t159;
t96 = -pkin(3) * t143 + t126;
t389 = g(1) * t115 + g(2) * t112 + g(3) * t153 + t92 * t96 + t291;
t388 = t112 * t225 - t168 * t230;
t387 = t112 * t230 + t168 * t225;
t219 = t223 ^ 2;
t313 = 0.2e1 * t219;
t214 = pkin(3) * t226 + pkin(10);
t384 = (pkin(3) * t144 + qJD(5) * t214 + t53) * t336;
t383 = t336 * t260;
t382 = qJD(4) * t132 + t403 * t226 + t402 * t231;
t259 = -t193 * t231 - t194 * t226;
t381 = qJD(4) * t259 - t402 * t226 + t403 * t231;
t369 = pkin(1) * t228;
t331 = pkin(7) * t345 + t224 * t369;
t155 = pkin(8) * t224 + t331;
t333 = t232 * t155 + t227 * t156;
t26 = -pkin(10) * t179 + t28;
t39 = pkin(4) * t92 - pkin(10) * t260 + t96;
t268 = t225 * t26 - t230 * t39;
t379 = t25 * t320 + t260 * t268;
t11 = t225 * t39 + t230 * t26;
t377 = t11 * t260 - t249 * t225 + t25 * t319;
t376 = -t260 * t96 + t243 + t252;
t373 = -t179 * t260 - t36;
t346 = t223 * t232;
t166 = t224 * t227 + t228 * t346;
t286 = -t155 * t227 + t232 * t156;
t64 = -pkin(3) * t345 - pkin(9) * t166 + t286;
t165 = -t224 * t232 + t227 * t348;
t71 = -pkin(9) * t165 + t333;
t263 = t226 * t64 + t231 * t71;
t299 = t223 * t325;
t120 = -qJD(3) * t165 + t232 * t299;
t161 = t223 * t255;
t209 = pkin(7) * t348;
t343 = t224 * t233;
t163 = (pkin(1) * t343 - t209) * qJD(2);
t240 = -qJD(3) * t333 + t232 * t161 - t163 * t227;
t326 = qJD(2) * t228;
t300 = t223 * t326;
t45 = pkin(3) * t300 - pkin(9) * t120 + t240;
t119 = qJD(3) * t166 + t223 * t298;
t253 = -t155 * t324 + t156 * t323 + t227 * t161 + t232 * t163;
t48 = -pkin(9) * t119 + t253;
t371 = -t263 * qJD(4) - t226 * t48 + t231 * t45;
t4 = pkin(10) * t151 - t291;
t281 = t393 * pkin(7) + t228 * t282 - t233 * t308;
t101 = -pkin(2) * t207 + t281;
t52 = pkin(3) * t82 + t101;
t7 = pkin(4) * t36 - pkin(10) * t35 + t52;
t1 = -t268 * qJD(5) + t225 * t7 + t230 * t4;
t367 = g(1) * t234;
t361 = pkin(4) * t302 + t382;
t356 = t143 * t192;
t355 = t144 * t192;
t351 = t169 * t227;
t350 = t174 * t230;
t349 = t219 * qJD(1) ^ 2;
t338 = t232 * t234;
t164 = pkin(7) * t299 + t228 * t310;
t220 = t228 ^ 2;
t330 = -t233 ^ 2 + t220;
t318 = qJD(2) - t208;
t307 = t233 * t349;
t306 = t225 * t345;
t305 = t230 * t345;
t285 = t169 * t232 - t227 * t344;
t284 = t208 + t328;
t283 = t207 + t316;
t37 = t226 * t58 + t357;
t276 = pkin(3) * t322 - t37;
t97 = pkin(3) * t119 + t164;
t271 = -g(1) * t171 - g(2) * t169;
t41 = -pkin(10) * t345 + t263;
t105 = t231 * t165 + t166 * t226;
t106 = -t165 * t226 + t166 * t231;
t154 = t209 + (-pkin(1) * t233 - pkin(2)) * t224;
t107 = pkin(3) * t165 + t154;
t49 = pkin(4) * t105 - pkin(10) * t106 + t107;
t267 = t225 * t49 + t230 * t41;
t266 = -t225 * t41 + t230 * t49;
t264 = -t226 * t71 + t231 * t64;
t87 = t106 * t225 + t305;
t256 = t226 * t45 + t231 * t48 + t64 * t321 - t71 * t322;
t104 = -t129 * t230 + t225 * t302;
t251 = -t230 * t117 - t174 * t320 - t104;
t246 = -pkin(8) * t157 - t192 * t126;
t2 = -t11 * qJD(5) - t225 * t4 + t230 * t7;
t38 = t231 * t58 - t360;
t239 = -t214 * t32 + t397 - (-pkin(3) * t321 + t38) * t336;
t238 = pkin(8) * qJD(3) * t192 - t101 - t248;
t236 = -g(3) * t348 - t25 * t117 - t132 * t32 + t5 * t174 - (pkin(10) * t302 - qJD(5) * t110 - t381) * t336 + t271;
t215 = -pkin(3) * t231 - pkin(4);
t123 = t171 * t232 + t227 * t347;
t122 = -t171 * t227 + t229 * t346;
t103 = -t129 * t225 - t230 * t302;
t88 = t106 * t230 - t306;
t85 = t115 * t230 + t170 * t225;
t84 = -t115 * t225 + t170 * t230;
t51 = qJD(4) * t106 + t231 * t119 + t120 * t226;
t50 = -qJD(4) * t105 - t119 * t226 + t120 * t231;
t40 = pkin(4) * t345 - t264;
t34 = -qJD(5) * t306 + t106 * t319 + t225 * t50 - t230 * t300;
t33 = -qJD(5) * t87 + t225 * t300 + t230 * t50;
t17 = pkin(4) * t51 - pkin(10) * t50 + t97;
t9 = -pkin(4) * t300 - t371;
t8 = pkin(10) * t300 + t256;
t3 = [qJDD(1), g(1) * t229 - g(2) * t234, g(2) * t229 + t367, (qJDD(1) * t220 + 0.2e1 * t228 * t296) * t219, (t228 * t314 - t330 * t317) * t313, (t283 * t228 + t284 * t325) * t223, (t283 * t233 - t284 * t326) * t223, t207 * t224, -t164 * t208 - t209 * t207 - t281 * t224 + g(1) * t169 - g(2) * t171 + (t207 * t343 + (-t297 + t314) * t313) * pkin(1), -t404 * pkin(1) * t313 - g(1) * t168 + g(2) * t170 - t163 * t208 - t331 * t207 - t247 * t224, t120 * t144 + t166 * t81, -t119 * t144 + t120 * t143 - t165 * t81 - t166 * t82, -t120 * t192 + t157 * t166 + (t144 * t326 - t233 * t81) * t223, t119 * t192 - t157 * t165 + (t143 * t326 + t233 * t82) * t223, (-t157 * t233 - t192 * t326) * t223, -t240 * t192 + t286 * t157 - t164 * t143 + t154 * t82 + t101 * t165 + t126 * t119 + g(1) * t285 - g(2) * t123 + (-t233 * t241 + t326 * t78) * t223, t253 * t192 - t333 * t157 + t164 * t144 + t154 * t81 + t101 * t166 + t126 * t120 - g(1) * t351 - g(2) * t122 + (-g(1) * t338 - t233 * t254 - t326 * t79) * t223, t106 * t35 + t260 * t50, -t105 * t35 - t106 * t36 - t260 * t51 - t50 * t92, t106 * t151 - t179 * t50 + (-t233 * t35 + t260 * t326) * t223, -t105 * t151 + t179 * t51 + (t233 * t36 - t326 * t92) * t223, (-t151 * t233 - t179 * t326) * t223, -t371 * t179 + t264 * t151 + t97 * t92 + t107 * t36 + t52 * t105 + t96 * t51 + g(1) * t112 - g(2) * t115 + (-t233 * t243 + t27 * t326) * t223, t256 * t179 - t263 * t151 + t97 * t260 + t107 * t35 + t52 * t106 + t96 * t50 - g(1) * t352 - g(2) * t114 + (-t218 * t367 - t233 * t291 - t28 * t326) * t223, t15 * t88 - t261 * t33, -t15 * t87 - t16 * t88 + t261 * t34 - t33 * t66, t105 * t15 - t261 * t51 + t32 * t88 - t33 * t336, -t105 * t16 - t32 * t87 + t336 * t34 - t51 * t66, t105 * t32 - t336 * t51, -(-qJD(5) * t267 + t17 * t230 - t225 * t8) * t336 + t266 * t32 + t2 * t105 - t268 * t51 + t9 * t66 + t40 * t16 + t5 * t87 + t25 * t34 + g(1) * t387 - g(2) * t85, (qJD(5) * t266 + t17 * t225 + t230 * t8) * t336 - t267 * t32 - t1 * t105 - t11 * t51 - t9 * t261 + t40 * t15 + t5 * t88 + t25 * t33 - g(1) * t388 - g(2) * t84; 0, 0, 0, -t228 * t307, t330 * t349, (t318 * t327 + t315) * t223, -t318 * t302 + t206, t207, t162 * t208 + t349 * t369 - t248 - t281, pkin(1) * t307 + t159 * t208 + (pkin(7) * t317 + g(3)) * t348 - t271 + t304, t227 * t81 - t232 * t355, (t81 - t356) * t232 + (-t82 + t355) * t227, -t192 * t323 + t157 * t227 + (-t144 * t228 + t192 * t339) * t329, t192 * t324 + t157 * t232 + (-t192 * t227 * t233 - t143 * t228) * t329, t192 * t302, -t78 * t302 - pkin(2) * t82 + t162 * t143 + t148 * t192 + (-t159 * t192 + t246) * t227 + t238 * t232, -pkin(2) * t81 - t162 * t144 - t192 * t332 - t227 * t238 + t232 * t246 + t302 * t79, t174 * t35 + t260 * t335, -t173 * t35 - t174 * t36 - t260 * t334 - t335 * t92, t151 * t174 - t179 * t335 - t260 * t302, -t151 * t173 + t179 * t334 + t302 * t92, t179 * t302, t151 * t259 + t52 * t173 + t179 * t382 + t216 * t36 - t27 * t302 + t275 * t92 + t334 * t96 - t245, -t132 * t151 + t52 * t174 + t179 * t381 + t216 * t35 + t248 * t217 + t275 * t260 + t28 * t302 + t335 * t96, t15 * t350 - t251 * t261, -t103 * t261 + t104 * t66 - (t225 * t261 - t230 * t66) * t117 + (-t13 - t16 * t230 + (t225 * t66 + t230 * t261) * qJD(5)) * t174, t15 * t173 - t251 * t336 - t261 * t334 + t32 * t350, -t174 * t29 - t16 * t173 - t334 * t66 - (t225 * t117 - t174 * t319 + t103) * t336, t173 * t32 - t334 * t336, -t25 * t103 - t259 * t16 + t2 * t173 + t236 * t225 - t230 * t392 - t334 * t268 + t361 * t66, -t1 * t173 - t25 * t104 - t334 * t11 - t259 * t15 + t225 * t392 + t236 * t230 - t361 * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144 * t143, -t143 ^ 2 + t144 ^ 2, t81 + t356, -t355 - t82, t157, -t79 * t192 - t126 * t144 - g(1) * t122 - g(2) * (-t223 * t338 - t351) + g(3) * t165 + t241, g(1) * t123 + g(2) * t285 + g(3) * t166 - t126 * t143 - t78 * t192 + t254, t395, t391, t390, t373, t151, -t179 * t37 + (-t144 * t92 + t151 * t231 + t179 * t322) * pkin(3) + t376, -t179 * t38 + (-t144 * t260 - t151 * t226 + t179 * t321) * pkin(3) + t389, t401, t398, t400, t399, t383, t215 * t16 + t276 * t66 + t239 * t225 + (t249 + t384) * t230 + t379, t215 * t15 - t225 * t384 + t230 * t239 - t261 * t276 + t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t395, t391, t390, t373, t151, -t179 * t28 + t376, -t179 * t27 + t389, t401, t398, t400, t399, t383, -pkin(4) * t16 - t28 * t66 + (-pkin(10) * t32 - t27 * t336 + t397) * t225 + (-(-pkin(10) * qJD(5) - t53) * t336 + t249) * t230 + t379, -pkin(4) * t15 - (t225 * t53 + t230 * t27) * t336 + t28 * t261 + t25 * t396 + (-t320 * t336 - t30) * pkin(10) + t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261 * t66, t261 ^ 2 - t66 ^ 2, -t336 * t66 + t15, t261 * t336 - t16, t32, -t11 * t336 + t25 * t261 - g(1) * t84 + g(2) * t388 - g(3) * (-t153 * t225 - t305) + t2, t268 * t336 + t25 * t66 + g(1) * t85 + g(2) * t387 - g(3) * (-t153 * t230 + t306) - t1;];
tau_reg = t3;
