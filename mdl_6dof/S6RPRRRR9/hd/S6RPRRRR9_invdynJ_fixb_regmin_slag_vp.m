% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR9
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:23
% EndTime: 2019-03-09 07:25:39
% DurationCPUTime: 6.54s
% Computational Cost: add. (6871->517), mult. (14024->708), div. (0->0), fcn. (10117->14), ass. (0->255)
t223 = sin(qJ(4));
t228 = cos(qJ(4));
t304 = t228 * qJD(3);
t229 = cos(qJ(3));
t319 = qJD(1) * t229;
t168 = t223 * t319 - t304;
t316 = qJD(3) * t223;
t170 = t228 * t319 + t316;
t222 = sin(qJ(5));
t227 = cos(qJ(5));
t102 = t227 * t168 + t170 * t222;
t221 = sin(qJ(6));
t226 = cos(qJ(6));
t260 = t168 * t222 - t227 * t170;
t261 = t102 * t221 + t226 * t260;
t45 = t226 * t102 - t221 * t260;
t382 = t261 * t45;
t376 = t261 ^ 2 - t45 ^ 2;
t224 = sin(qJ(3));
t320 = qJD(1) * t224;
t202 = qJD(4) + t320;
t197 = qJD(5) + t202;
t189 = qJD(6) + t197;
t309 = qJD(4) * t229;
t249 = -t223 * t309 - t224 * t304;
t299 = t229 * qJDD(1);
t91 = qJD(1) * t249 + qJD(4) * t304 + t223 * qJDD(3) + t228 * t299;
t315 = qJD(3) * t224;
t290 = t223 * t315;
t92 = -qJD(1) * t290 + qJD(4) * t170 - t228 * qJDD(3) + t223 * t299;
t238 = qJD(5) * t260 - t222 * t91 - t227 * t92;
t307 = qJD(5) * t227;
t308 = qJD(5) * t222;
t30 = -t168 * t307 - t170 * t308 - t222 * t92 + t227 * t91;
t305 = qJD(6) * t226;
t306 = qJD(6) * t221;
t6 = -t102 * t305 + t221 * t238 + t226 * t30 + t260 * t306;
t374 = t189 * t45 + t6;
t220 = qJ(4) + qJ(5);
t216 = qJ(6) + t220;
t206 = sin(t216);
t207 = cos(t216);
t230 = cos(qJ(1));
t225 = sin(qJ(1));
t334 = t224 * t225;
t120 = t206 * t230 + t207 * t334;
t332 = t224 * t230;
t122 = -t206 * t225 + t207 * t332;
t181 = pkin(3) * t224 - pkin(8) * t229 + qJ(2);
t147 = t181 * qJD(1);
t231 = -pkin(1) - pkin(7);
t199 = qJD(1) * t231 + qJD(2);
t180 = t224 * t199;
t157 = qJD(3) * pkin(8) + t180;
t86 = t228 * t147 - t157 * t223;
t66 = -pkin(9) * t170 + t86;
t61 = pkin(4) * t202 + t66;
t87 = t147 * t223 + t157 * t228;
t67 = -pkin(9) * t168 + t87;
t65 = t227 * t67;
t25 = t222 * t61 + t65;
t384 = pkin(10) * t102;
t19 = t25 - t384;
t17 = t19 * t306;
t357 = g(3) * t229;
t340 = t199 * t229;
t158 = -qJD(3) * pkin(3) - t340;
t114 = pkin(4) * t168 + t158;
t54 = pkin(5) * t102 + t114;
t373 = g(1) * t120 - g(2) * t122 + t207 * t357 + t45 * t54 + t17;
t119 = -t206 * t334 + t207 * t230;
t121 = t206 * t332 + t207 * t225;
t303 = qJD(1) * qJD(3);
t285 = t229 * t303;
t300 = t224 * qJDD(1);
t165 = qJDD(4) + t285 + t300;
t160 = qJDD(5) + t165;
t195 = qJDD(1) * t231 + qJDD(2);
t314 = qJD(3) * t229;
t116 = qJDD(3) * pkin(8) + t195 * t224 + t199 * t314;
t267 = pkin(3) * t229 + pkin(8) * t224;
t166 = qJD(3) * t267 + qJD(2);
t100 = qJD(1) * t166 + qJDD(1) * t181;
t94 = t228 * t100;
t16 = pkin(4) * t165 - pkin(9) * t91 - qJD(4) * t87 - t116 * t223 + t94;
t310 = qJD(4) * t228;
t295 = -t223 * t100 - t228 * t116 - t147 * t310;
t312 = qJD(4) * t223;
t253 = -t157 * t312 - t295;
t21 = -pkin(9) * t92 + t253;
t282 = t227 * t16 - t222 * t21;
t240 = -qJD(5) * t25 + t282;
t2 = pkin(5) * t160 - pkin(10) * t30 + t240;
t277 = -t222 * t16 - t227 * t21 - t61 * t307 + t67 * t308;
t3 = pkin(10) * t238 - t277;
t293 = t226 * t2 - t221 * t3;
t388 = -g(1) * t119 - g(2) * t121 + t206 * t357 + t54 * t261 + t293;
t239 = qJD(6) * t261 - t221 * t30 + t226 * t238;
t368 = -t189 * t261 + t239;
t291 = t223 * t320;
t360 = pkin(8) + pkin(9);
t292 = qJD(4) * t360;
t177 = t267 * qJD(1);
t329 = t228 * t229;
t325 = t223 * t177 + t199 * t329;
t387 = pkin(9) * t291 + t223 * t292 + t325;
t156 = t228 * t177;
t333 = t224 * t228;
t298 = pkin(9) * t333;
t336 = t223 * t229;
t386 = t228 * t292 - t199 * t336 + t156 + (pkin(4) * t229 + t298) * qJD(1);
t265 = g(1) * t225 - g(2) * t230;
t115 = -qJDD(3) * pkin(3) - t195 * t229 + t199 * t315;
t358 = g(3) * t224;
t244 = t229 * t265 - t358;
t385 = qJD(4) * pkin(8) * t202 + t115 + t244;
t383 = pkin(10) * t260;
t380 = t260 * t102;
t173 = t222 * t228 + t223 * t227;
t361 = qJD(4) + qJD(5);
t379 = t361 * t173;
t172 = t222 * t223 - t227 * t228;
t365 = t224 * t172;
t327 = -qJD(1) * t365 - t361 * t172;
t146 = t173 * qJD(1);
t326 = t224 * t146 + t379;
t378 = qJDD(2) - t265;
t377 = t228 * t309 - t290;
t375 = -t102 ^ 2 + t260 ^ 2;
t372 = t102 * t197 + t30;
t63 = t222 * t67;
t24 = t227 * t61 - t63;
t18 = t24 + t383;
t15 = pkin(5) * t197 + t18;
t351 = t19 * t226;
t5 = t15 * t221 + t351;
t371 = -qJD(6) * t5 + t388;
t214 = sin(t220);
t215 = cos(t220);
t130 = t214 * t230 + t215 * t334;
t132 = -t214 * t225 + t215 * t332;
t370 = g(1) * t130 - g(2) * t132 + t102 * t114 + t215 * t357 + t277;
t129 = -t214 * t334 + t215 * t230;
t131 = t214 * t332 + t215 * t225;
t369 = -g(1) * t129 - g(2) * t131 + t114 * t260 + t214 * t357 + t240;
t367 = -t197 * t260 + t238;
t196 = t231 * t333;
t323 = t223 * t181 + t196;
t113 = -pkin(9) * t336 + t323;
t164 = t228 * t181;
t335 = t223 * t231;
t284 = pkin(4) - t335;
t98 = -pkin(9) * t329 + t224 * t284 + t164;
t349 = t227 * t113 + t222 * t98;
t366 = t386 * t227;
t268 = -t180 + (t291 + t312) * pkin(4);
t190 = t360 * t223;
t191 = t360 * t228;
t324 = -t222 * t190 + t227 * t191;
t363 = t190 * t307 + t191 * t308 + t386 * t222 + t227 * t387;
t346 = pkin(1) * qJDD(1);
t362 = t346 - t378;
t108 = t226 * t172 + t173 * t221;
t356 = -qJD(6) * t108 - t221 * t326 + t226 * t327;
t109 = -t172 * t221 + t173 * t226;
t355 = qJD(6) * t109 + t221 * t327 + t226 * t326;
t354 = t227 * t66 - t63;
t352 = pkin(5) * t326 + t268;
t350 = t223 * t91;
t348 = t172 * qJD(1) - t173 * t314 + t361 * t365;
t137 = t172 * t229;
t347 = qJD(3) * t137 + t224 * t379 + t146;
t233 = qJD(1) ^ 2;
t345 = qJ(2) * t233;
t141 = qJDD(6) + t160;
t344 = t141 * t222;
t343 = t168 * t202;
t342 = t170 * t202;
t341 = t170 * t228;
t208 = pkin(4) * t227 + pkin(5);
t339 = t208 * t141;
t338 = t222 * t226;
t337 = t223 * t165;
t331 = t225 * t228;
t330 = t228 * t165;
t328 = t228 * t230;
t219 = t229 ^ 2;
t322 = t224 ^ 2 - t219;
t232 = qJD(3) ^ 2;
t321 = -t232 - t233;
t318 = qJD(3) * t168;
t317 = qJD(3) * t170;
t313 = qJD(3) * t231;
t311 = qJD(4) * t224;
t302 = qJDD(1) * qJ(2);
t301 = qJDD(3) * t224;
t288 = t229 * t313;
t294 = t223 * t166 + t181 * t310 + t228 * t288;
t209 = -pkin(4) * t228 - pkin(3);
t283 = qJD(6) * t15 + t3;
t140 = t228 * t166;
t42 = t140 + (-t196 + (pkin(9) * t229 - t181) * t223) * qJD(4) + (t229 * t284 + t298) * qJD(3);
t51 = -pkin(9) * t377 - t311 * t335 + t294;
t280 = -t222 * t51 + t227 * t42;
t279 = -t222 * t66 - t65;
t276 = -t113 * t222 + t227 * t98;
t274 = t202 * t231 + t157;
t273 = -t227 * t190 - t191 * t222;
t167 = pkin(4) * t336 - t229 * t231;
t272 = -qJD(4) * t147 - t116;
t271 = qJD(1) + t311;
t82 = -pkin(10) * t172 + t324;
t270 = pkin(5) * t319 + pkin(10) * t327 + t324 * qJD(5) + qJD(6) * t82 - t222 * t387 + t366;
t81 = -pkin(10) * t173 + t273;
t269 = pkin(10) * t326 - qJD(6) * t81 + t363;
t266 = g(1) * t230 + g(2) * t225;
t134 = t173 * t224;
t264 = qJD(6) * t134 + t347;
t263 = -qJD(6) * t365 - t348;
t135 = t173 * t229;
t74 = t226 * t135 - t137 * t221;
t75 = -t135 * t221 - t137 * t226;
t257 = t202 * t310 + t337;
t256 = -t202 * t312 + t330;
t255 = -t113 * t308 + t222 * t42 + t227 * t51 + t98 * t307;
t117 = pkin(4) * t377 + t224 * t313;
t251 = 0.2e1 * qJ(2) * t303 + qJDD(3) * t231;
t248 = 0.2e1 * qJD(1) * qJD(2) - t266;
t247 = -t195 + t265 + t345;
t246 = -pkin(8) * t165 + t158 * t202;
t52 = pkin(4) * t92 + t115;
t241 = t248 + 0.2e1 * t302;
t237 = -t231 * t232 + t241;
t212 = qJDD(3) * t229;
t151 = -t223 * t225 + t224 * t328;
t150 = t223 * t332 + t331;
t149 = t223 * t230 + t224 * t331;
t148 = -t223 * t334 + t328;
t127 = pkin(5) * t172 + t209;
t107 = pkin(5) * t135 + t167;
t68 = pkin(4) * t170 - pkin(5) * t260;
t60 = -t308 * t336 + (t329 * t361 - t290) * t227 + t249 * t222;
t58 = qJD(3) * t365 - t229 * t379;
t39 = pkin(5) * t60 + t117;
t38 = -pkin(10) * t135 + t349;
t35 = pkin(5) * t224 + pkin(10) * t137 + t276;
t23 = t354 + t383;
t22 = t279 + t384;
t12 = qJD(6) * t75 + t221 * t58 + t226 * t60;
t11 = -qJD(6) * t74 - t221 * t60 + t226 * t58;
t10 = -pkin(5) * t238 + t52;
t9 = -pkin(10) * t60 + t255;
t8 = pkin(5) * t314 - pkin(10) * t58 - qJD(5) * t349 + t280;
t4 = t15 * t226 - t19 * t221;
t1 = [qJDD(1), t265, t266, -0.2e1 * t346 + t378, t241, t362 * pkin(1) + (t248 + t302) * qJ(2), qJDD(1) * t219 - 0.2e1 * t224 * t285, -0.2e1 * t224 * t299 + 0.2e1 * t303 * t322, -t224 * t232 + t212, -t229 * t232 - t301, 0, t224 * t237 + t229 * t251, -t224 * t251 + t229 * t237, t170 * t249 + t329 * t91 (t168 * t228 + t170 * t223) * t315 + (-t350 - t228 * t92 + (t168 * t223 - t341) * qJD(4)) * t229 (-t202 * t304 + t91) * t224 + (t256 + t317) * t229 (t202 * t316 - t92) * t224 + (-t257 - t318) * t229, t165 * t224 + t202 * t314, -g(1) * t151 - g(2) * t149 + t140 * t202 + t164 * t165 + (t168 * t313 - t274 * t310 + t94) * t224 + (qJD(3) * t86 + t158 * t310 - t231 * t92) * t229 + ((-qJD(4) * t181 - t288) * t202 + t115 * t229 + (-qJD(3) * t158 - t165 * t231 + t272) * t224) * t223, -t294 * t202 - t323 * t165 + g(1) * t150 - g(2) * t148 + (t274 * t312 + (-t158 * t228 + t170 * t231) * qJD(3) + t295) * t224 + (-qJD(3) * t87 + t115 * t228 - t158 * t312 - t231 * t91) * t229, -t137 * t30 - t260 * t58, -t102 * t58 - t135 * t30 - t137 * t238 + t260 * t60, -t137 * t160 + t197 * t58 + t224 * t30 - t260 * t314, -t102 * t314 - t135 * t160 - t197 * t60 + t224 * t238, t160 * t224 + t197 * t314, t280 * t197 + t276 * t160 + t282 * t224 + t24 * t314 + t117 * t102 - t167 * t238 + t52 * t135 + t114 * t60 - g(1) * t132 - g(2) * t130 + (-t197 * t349 - t224 * t25) * qJD(5), g(1) * t131 - g(2) * t129 + t114 * t58 - t117 * t260 - t52 * t137 - t160 * t349 + t167 * t30 - t197 * t255 + t224 * t277 - t25 * t314, -t11 * t261 + t6 * t75, -t11 * t45 + t12 * t261 + t239 * t75 - t6 * t74, t11 * t189 + t141 * t75 + t224 * t6 - t261 * t314, -t12 * t189 - t141 * t74 + t224 * t239 - t314 * t45, t141 * t224 + t189 * t314 (-t221 * t9 + t226 * t8) * t189 + (-t221 * t38 + t226 * t35) * t141 + t293 * t224 + t4 * t314 + t39 * t45 - t107 * t239 + t10 * t74 + t54 * t12 - g(1) * t122 - g(2) * t120 + ((-t221 * t35 - t226 * t38) * t189 - t5 * t224) * qJD(6), -t5 * t314 + g(1) * t121 - g(2) * t119 + t10 * t75 + t107 * t6 + t54 * t11 + t17 * t224 - t39 * t261 + (-(-qJD(6) * t38 + t8) * t189 - t35 * t141 - t2 * t224) * t221 + (-(qJD(6) * t35 + t9) * t189 - t38 * t141 - t283 * t224) * t226; 0, 0, 0, qJDD(1), -t233, -t345 - t362, 0, 0, 0, 0, 0, t224 * t321 + t212, t229 * t321 - t301, 0, 0, 0, 0, 0, -t229 * t92 + (t318 - t337) * t224 + (-t223 * t314 - t228 * t271) * t202, -t229 * t91 + (t317 - t330) * t224 + (t223 * t271 - t229 * t304) * t202, 0, 0, 0, 0, 0, t102 * t315 - t134 * t160 + t197 * t348 + t229 * t238, t160 * t365 + t197 * t347 - t229 * t30 - t260 * t315, 0, 0, 0, 0, 0 (-t134 * t226 + t221 * t365) * t141 + t45 * t315 + t229 * t239 + (t221 * t264 - t226 * t263) * t189 -(-t134 * t221 - t226 * t365) * t141 - t261 * t315 - t229 * t6 + (t221 * t263 + t226 * t264) * t189; 0, 0, 0, 0, 0, 0, t229 * t233 * t224, -t322 * t233, t299, -t300, qJDD(3), -t229 * t247 + t358, t224 * t247 + t357, t202 * t341 + t350 (t91 - t343) * t228 + (-t92 - t342) * t223 (-t170 * t229 + t202 * t333) * qJD(1) + t257 (-t202 * t223 * t224 + t168 * t229) * qJD(1) + t256, -t202 * t319, -t86 * t319 - t168 * t180 - pkin(3) * t92 - t156 * t202 + (t202 * t340 + t246) * t223 - t385 * t228, -pkin(3) * t91 - t170 * t180 + t325 * t202 + t223 * t385 + t246 * t228 + t87 * t319, t173 * t30 - t260 * t327, -t102 * t327 - t172 * t30 + t173 * t238 + t260 * t326, t160 * t173 + t197 * t327 + t260 * t319, t102 * t319 - t160 * t172 - t197 * t326, -t197 * t319, t273 * t160 - t209 * t238 + t52 * t172 - t24 * t319 + (-t191 * t307 + (qJD(5) * t190 + t387) * t222 - t366) * t197 + t326 * t114 + t268 * t102 - t244 * t215, t327 * t114 - t324 * t160 + t52 * t173 + t197 * t363 + t209 * t30 + t244 * t214 + t25 * t319 - t260 * t268, t109 * t6 - t261 * t356, -t108 * t6 + t109 * t239 + t261 * t355 - t356 * t45, t109 * t141 + t189 * t356 + t261 * t319, -t108 * t141 - t189 * t355 + t319 * t45, -t189 * t319 (-t221 * t82 + t226 * t81) * t141 - t127 * t239 + t10 * t108 - t4 * t319 + t355 * t54 + t352 * t45 + (t221 * t269 - t226 * t270) * t189 - t244 * t207 -(t221 * t81 + t226 * t82) * t141 + t127 * t6 + t10 * t109 + t5 * t319 + t356 * t54 - t352 * t261 + (t221 * t270 + t226 * t269) * t189 + t244 * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170 * t168, -t168 ^ 2 + t170 ^ 2, t91 + t343, t342 - t92, t165, -t157 * t310 - g(1) * t148 - g(2) * t150 - t158 * t170 + t202 * t87 + t94 + (t272 + t357) * t223, g(1) * t149 - g(2) * t151 + g(3) * t329 + t158 * t168 + t202 * t86 - t253, -t380, t375, t372, t367, t160, -t279 * t197 + (-t102 * t170 + t160 * t227 - t197 * t308) * pkin(4) + t369, t354 * t197 + (-t160 * t222 + t170 * t260 - t197 * t307) * pkin(4) + t370, -t382, t376, t374, t368, t141, t226 * t339 - (t22 * t226 - t221 * t23) * t189 - t68 * t45 + (-t221 * t344 + (-t221 * t227 - t338) * t189 * qJD(5)) * pkin(4) + ((-pkin(4) * t338 - t208 * t221) * t189 - t5) * qJD(6) + t388, t68 * t261 + (-t339 - t2 + (t22 - (-qJD(5) - qJD(6)) * t222 * pkin(4)) * t189) * t221 + (-pkin(4) * t344 + (-pkin(4) * t307 - qJD(6) * t208 + t23) * t189 - t283) * t226 + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, t375, t372, t367, t160, t197 * t25 + t369, t197 * t24 + t370, -t382, t376, t374, t368, t141 -(-t18 * t221 - t351) * t189 + (t141 * t226 - t189 * t306 + t260 * t45) * pkin(5) + t371 (-t189 * t19 - t2) * t221 + (t18 * t189 - t283) * t226 + (-t141 * t221 - t189 * t305 - t260 * t261) * pkin(5) + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t382, t376, t374, t368, t141, t189 * t5 + t371, t189 * t4 - t221 * t2 - t226 * t283 + t373;];
tau_reg  = t1;
