% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR3
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
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:11
% EndTime: 2019-03-09 10:19:27
% DurationCPUTime: 6.72s
% Computational Cost: add. (9856->492), mult. (23367->660), div. (0->0), fcn. (18124->16), ass. (0->260)
t248 = sin(pkin(10));
t255 = sin(qJ(2));
t321 = qJD(1) * t255;
t250 = cos(pkin(10));
t259 = cos(qJ(2));
t334 = t250 * t259;
t199 = qJD(1) * t334 - t248 * t321;
t189 = qJD(4) - t199;
t214 = t248 * t259 + t250 * t255;
t202 = t214 * qJD(1);
t254 = sin(qJ(4));
t258 = cos(qJ(4));
t316 = t258 * qJD(2);
t168 = t202 * t254 - t316;
t170 = qJD(2) * t254 + t202 * t258;
t247 = sin(pkin(11));
t249 = cos(pkin(11));
t108 = t168 * t247 - t170 * t249;
t253 = sin(qJ(6));
t257 = cos(qJ(6));
t317 = qJD(6) * t253;
t285 = -t168 * t249 - t170 * t247;
t329 = t257 * t285;
t315 = qJD(1) * qJD(2);
t305 = t259 * t315;
t306 = t255 * t315;
t155 = qJDD(1) * t214 - t248 * t306 + t250 * t305;
t319 = qJD(4) * t254;
t100 = qJD(4) * t316 + t254 * qJDD(2) + t258 * t155 - t202 * t319;
t101 = qJD(4) * t170 - t258 * qJDD(2) + t155 * t254;
t49 = -t100 * t247 - t101 * t249;
t50 = t100 * t249 - t101 * t247;
t12 = qJD(6) * t329 + t108 * t317 + t253 * t49 + t257 * t50;
t183 = qJD(6) + t189;
t60 = t108 * t253 + t329;
t357 = t183 * t60;
t395 = t12 - t357;
t380 = -t257 * t108 + t253 * t285;
t394 = t380 * t60;
t213 = t247 * t258 + t249 * t254;
t383 = t189 * t213;
t282 = t247 * t254 - t249 * t258;
t382 = t189 * t282;
t393 = t380 ^ 2 - t60 ^ 2;
t241 = qJ(4) + pkin(11) + qJ(6);
t229 = sin(t241);
t230 = cos(t241);
t260 = cos(qJ(1));
t244 = qJ(2) + pkin(10);
t238 = cos(t244);
t256 = sin(qJ(1));
t339 = t238 * t256;
t172 = t229 * t260 - t230 * t339;
t338 = t238 * t260;
t174 = t229 * t256 + t230 * t338;
t200 = t214 * qJD(2);
t313 = t259 * qJDD(1);
t314 = t255 * qJDD(1);
t287 = -t248 * t314 + t250 * t313;
t149 = qJD(1) * t200 + qJDD(4) - t287;
t154 = -qJD(2) * t202 + t287;
t371 = pkin(2) * t259;
t236 = pkin(1) + t371;
t270 = pkin(2) * t306 - qJDD(1) * t236 + qJDD(3);
t84 = -pkin(3) * t154 - pkin(8) * t155 + t270;
t79 = t258 * t84;
t220 = -qJD(1) * t236 + qJD(3);
t119 = -pkin(3) * t199 - pkin(8) * t202 + t220;
t252 = -qJ(3) - pkin(7);
t221 = t252 * t255;
t218 = qJD(1) * t221;
t359 = qJD(2) * pkin(2);
t208 = t218 + t359;
t222 = t252 * t259;
t219 = qJD(1) * t222;
t335 = t250 * t219;
t153 = t248 * t208 - t335;
t140 = qJD(2) * pkin(8) + t153;
t81 = t119 * t254 + t140 * t258;
t300 = qJD(2) * t252;
t196 = -qJD(3) * t255 + t259 * t300;
t148 = qJDD(2) * pkin(2) + qJD(1) * t196 + qJDD(1) * t221;
t195 = qJD(3) * t259 + t255 * t300;
t158 = qJD(1) * t195 - qJDD(1) * t222;
t96 = t248 * t148 + t250 * t158;
t94 = qJDD(2) * pkin(8) + t96;
t16 = pkin(4) * t149 - qJ(5) * t100 - qJD(4) * t81 - qJD(5) * t170 - t254 * t94 + t79;
t318 = qJD(4) * t258;
t273 = t119 * t318 - t140 * t319 + t254 * t84 + t258 * t94;
t18 = -qJ(5) * t101 - qJD(5) * t168 + t273;
t4 = t249 * t16 - t18 * t247;
t2 = pkin(5) * t149 - pkin(9) * t50 + t4;
t68 = -qJ(5) * t168 + t81;
t354 = t249 * t68;
t80 = t258 * t119 - t140 * t254;
t67 = -qJ(5) * t170 + t80;
t53 = pkin(4) * t189 + t67;
t30 = t247 * t53 + t354;
t379 = pkin(9) * t285;
t21 = t30 + t379;
t20 = t21 * t317;
t237 = sin(t244);
t364 = g(3) * t237;
t205 = t248 * t219;
t152 = t208 * t250 + t205;
t139 = -qJD(2) * pkin(3) - t152;
t104 = pkin(4) * t168 + qJD(5) + t139;
t54 = -pkin(5) * t285 + t104;
t392 = g(1) * t174 - g(2) * t172 - t253 * t2 + t230 * t364 - t54 * t60 + t20;
t13 = qJD(6) * t380 + t253 * t50 - t257 * t49;
t358 = t183 * t380;
t390 = -t13 + t358;
t231 = pkin(2) * t248 + pkin(8);
t327 = qJ(5) + t231;
t297 = qJD(4) * t327;
t130 = pkin(2) * t321 + pkin(3) * t202 - pkin(8) * t199;
t160 = t218 * t250 + t205;
t324 = t254 * t130 + t258 * t160;
t343 = t199 * t254;
t389 = qJ(5) * t343 + qJD(5) * t258 - t254 * t297 - t324;
t125 = t258 * t130;
t388 = -pkin(4) * t202 - t125 + (qJ(5) * t199 - t297) * t258 + (-qJD(5) + t160) * t254;
t294 = g(1) * t260 + g(2) * t256;
t267 = -g(3) * t238 + t237 * t294;
t95 = t148 * t250 - t248 * t158;
t93 = -qJDD(2) * pkin(3) - t95;
t387 = -qJD(4) * t231 * t189 + t267 - t93;
t386 = t319 - t343;
t171 = t229 * t339 + t230 * t260;
t173 = -t229 * t338 + t230 * t256;
t5 = t247 * t16 + t249 * t18;
t3 = pkin(9) * t49 + t5;
t309 = t257 * t2 - t253 * t3;
t385 = -g(1) * t173 + g(2) * t171 + t229 * t364 - t54 * t380 + t309;
t384 = pkin(9) * t108;
t157 = t213 * t257 - t253 * t282;
t360 = qJD(6) * t157 - t253 * t382 + t257 * t383;
t212 = t248 * t255 - t334;
t204 = t212 * qJD(2);
t308 = t214 * t318;
t381 = -t204 * t254 + t308;
t352 = -t247 * t389 + t249 * t388;
t351 = t247 * t388 + t249 * t389;
t376 = g(1) * t256 - g(2) * t260;
t159 = t218 * t248 - t335;
t375 = pkin(4) * t386 - t159;
t328 = t258 * t260;
t332 = t254 * t256;
t190 = t238 * t332 + t328;
t330 = t256 * t258;
t331 = t254 * t260;
t192 = -t238 * t331 + t330;
t374 = -g(1) * t192 + g(2) * t190;
t142 = qJDD(6) + t149;
t284 = -t213 * t253 - t257 * t282;
t361 = qJD(6) * t284 - t253 * t383 - t257 * t382;
t373 = -t142 * t157 - t183 * t361;
t372 = pkin(2) * t250;
t370 = pkin(4) * t247;
t362 = g(3) * t259;
t312 = t255 * t359;
t131 = pkin(3) * t200 + pkin(8) * t204 + t312;
t126 = t258 * t131;
t129 = t195 * t250 + t196 * t248;
t151 = pkin(3) * t212 - pkin(8) * t214 - t236;
t163 = t221 * t248 - t222 * t250;
t161 = t258 * t163;
t281 = qJ(5) * t204 - qJD(5) * t214;
t37 = pkin(4) * t200 - t129 * t254 + t126 + t281 * t258 + (-t161 + (qJ(5) * t214 - t151) * t254) * qJD(4);
t310 = t258 * t129 + t254 * t131 + t151 * t318;
t41 = -qJ(5) * t308 + (-qJD(4) * t163 + t281) * t254 + t310;
t11 = t247 * t37 + t249 * t41;
t62 = t247 * t68;
t33 = t249 * t67 - t62;
t138 = t258 * t151;
t340 = t214 * t258;
t72 = pkin(4) * t212 - qJ(5) * t340 - t163 * t254 + t138;
t323 = t254 * t151 + t161;
t341 = t214 * t254;
t85 = -qJ(5) * t341 + t323;
t43 = t247 * t72 + t249 * t85;
t356 = t202 * t60;
t355 = t202 * t380;
t29 = t249 * t53 - t62;
t19 = pkin(5) * t189 + t29 + t384;
t353 = t257 * t19;
t350 = pkin(5) * t383 + t375;
t349 = t100 * t254;
t347 = t168 * t189;
t346 = t168 * t202;
t345 = t170 * t189;
t344 = t170 * t202;
t333 = t254 * t149;
t136 = t258 * t149;
t209 = t327 * t254;
t210 = t327 * t258;
t144 = -t247 * t209 + t249 * t210;
t245 = t255 ^ 2;
t322 = -t259 ^ 2 + t245;
t320 = qJD(4) * t214;
t235 = pkin(4) * t258 + pkin(3);
t304 = pkin(4) * t254 - t252;
t302 = qJD(6) * t19 + t3;
t10 = -t247 * t41 + t249 * t37;
t32 = -t247 * t67 - t354;
t42 = -t247 * t85 + t249 * t72;
t299 = -qJD(4) * t119 - t94;
t128 = t195 * t248 - t250 * t196;
t143 = -t249 * t209 - t210 * t247;
t162 = -t250 * t221 - t222 * t248;
t296 = t189 * t258;
t295 = t284 * t142 - t183 * t360;
t292 = -t140 * t318 + t79;
t115 = -pkin(9) * t282 + t144;
t291 = pkin(5) * t202 - pkin(9) * t382 + qJD(6) * t115 - t352;
t114 = -pkin(9) * t213 + t143;
t290 = pkin(9) * t383 - qJD(6) * t114 - t351;
t288 = pkin(4) * t341 + t162;
t7 = t253 * t19 + t257 * t21;
t133 = t213 * t214;
t134 = t282 * t214;
t286 = -t257 * t133 + t134 * t253;
t89 = -t133 * t253 - t134 * t257;
t251 = -qJ(5) - pkin(8);
t283 = t235 * t238 - t237 * t251;
t280 = -t235 - t372;
t279 = pkin(4) * t381 + t128;
t278 = -t189 * t386 + t136;
t232 = pkin(4) * t249 + pkin(5);
t277 = t232 * t253 + t257 * t370;
t276 = t232 * t257 - t253 * t370;
t275 = -0.2e1 * pkin(1) * t315 - pkin(7) * qJDD(2);
t274 = -t204 * t258 - t214 * t319;
t272 = t139 * t189 - t231 * t149;
t261 = qJD(2) ^ 2;
t265 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t261 + t376;
t262 = qJD(1) ^ 2;
t264 = pkin(1) * t262 - pkin(7) * qJDD(1) + t294;
t45 = pkin(4) * t101 + qJDD(5) + t93;
t233 = -pkin(3) - t372;
t224 = t260 * t236;
t193 = t238 * t328 + t332;
t191 = -t238 * t330 + t331;
t167 = pkin(5) * t282 + t280;
t91 = pkin(5) * t133 + t288;
t90 = pkin(4) * t170 - pkin(5) * t108;
t87 = -t204 * t282 + t213 * t320;
t86 = t204 * t213 + t282 * t320;
t44 = -pkin(5) * t86 + t279;
t36 = -pkin(9) * t133 + t43;
t31 = pkin(5) * t212 + pkin(9) * t134 + t42;
t28 = qJD(6) * t89 - t253 * t87 - t257 * t86;
t27 = qJD(6) * t286 + t253 * t86 - t257 * t87;
t24 = t33 + t384;
t23 = t32 - t379;
t22 = -pkin(5) * t49 + t45;
t9 = pkin(9) * t86 + t11;
t8 = pkin(5) * t200 + pkin(9) * t87 + t10;
t6 = -t21 * t253 + t353;
t1 = [qJDD(1), t376, t294, qJDD(1) * t245 + 0.2e1 * t255 * t305, 0.2e1 * t255 * t313 - 0.2e1 * t322 * t315, qJDD(2) * t255 + t259 * t261, qJDD(2) * t259 - t255 * t261, 0, t255 * t275 + t259 * t265, -t255 * t265 + t259 * t275, t128 * t202 + t129 * t199 + t152 * t204 - t153 * t200 + t154 * t163 + t155 * t162 - t212 * t96 - t214 * t95 - t294, t96 * t163 + t153 * t129 - t95 * t162 - t152 * t128 - t270 * t236 + t220 * t312 - g(1) * (-t236 * t256 - t252 * t260) - g(2) * (-t252 * t256 + t224) t100 * t340 + t170 * t274 -(-t168 * t258 - t170 * t254) * t204 + (-t349 - t101 * t258 + (t168 * t254 - t170 * t258) * qJD(4)) * t214, t100 * t212 + t136 * t214 + t170 * t200 + t189 * t274, -t101 * t212 - t168 * t200 - t189 * t381 - t214 * t333, t149 * t212 + t189 * t200 (-t163 * t318 + t126) * t189 + t138 * t149 + t292 * t212 + t80 * t200 + t128 * t168 + t162 * t101 + t139 * t308 - g(1) * t191 - g(2) * t193 + ((-qJD(4) * t151 - t129) * t189 - t163 * t149 + t299 * t212 + t93 * t214 - t139 * t204) * t254 -(-t163 * t319 + t310) * t189 - t323 * t149 - t273 * t212 - t81 * t200 + t128 * t170 + t162 * t100 + t93 * t340 - g(1) * t190 - g(2) * t192 + t274 * t139, t10 * t108 + t11 * t285 - t133 * t5 + t134 * t4 + t237 * t376 + t29 * t87 + t30 * t86 - t42 * t50 + t43 * t49, t5 * t43 + t30 * t11 + t4 * t42 + t29 * t10 + t45 * t288 + t104 * t279 - g(2) * t224 + (-g(1) * t304 - g(2) * t283) * t260 + (-g(1) * (-t236 - t283) - g(2) * t304) * t256, t12 * t89 + t27 * t380, t12 * t286 - t13 * t89 + t27 * t60 - t28 * t380, t12 * t212 + t142 * t89 + t183 * t27 + t200 * t380, -t13 * t212 + t142 * t286 - t183 * t28 + t200 * t60, t142 * t212 + t183 * t200 (-t253 * t9 + t257 * t8) * t183 + (-t253 * t36 + t257 * t31) * t142 + t309 * t212 + t6 * t200 - t44 * t60 + t91 * t13 - t22 * t286 + t54 * t28 - g(1) * t172 - g(2) * t174 + ((-t253 * t31 - t257 * t36) * t183 - t7 * t212) * qJD(6), -g(1) * t171 - g(2) * t173 + t91 * t12 + t20 * t212 - t7 * t200 + t22 * t89 + t54 * t27 + t44 * t380 + (-(-qJD(6) * t36 + t8) * t183 - t31 * t142 - t2 * t212) * t253 + (-(qJD(6) * t31 + t9) * t183 - t36 * t142 - t302 * t212) * t257; 0, 0, 0, -t255 * t262 * t259, t322 * t262, t314, t313, qJDD(2), t255 * t264 - t362, g(3) * t255 + t259 * t264 (t153 - t159) * t202 + (t152 - t160) * t199 + (t154 * t248 - t155 * t250) * pkin(2), t152 * t159 - t153 * t160 + (-t362 + t248 * t96 + t250 * t95 + (-qJD(1) * t220 + t294) * t255) * pkin(2), t170 * t296 + t349 (t100 - t347) * t258 + (-t101 - t345) * t254, t189 * t296 + t333 - t344, t278 + t346, -t189 * t202, t233 * t101 - t125 * t189 - t159 * t168 - t80 * t202 + (t160 * t189 + t272) * t254 + t387 * t258, t233 * t100 - t159 * t170 + t324 * t189 + t81 * t202 - t254 * t387 + t272 * t258, t108 * t352 - t143 * t50 + t144 * t49 - t213 * t4 - t238 * t294 - t282 * t5 + t285 * t351 + t29 * t382 - t30 * t383 - t364, t5 * t144 + t4 * t143 + t45 * t280 - g(3) * (t283 + t371) + t351 * t30 + t352 * t29 + t375 * t104 + t294 * (pkin(2) * t255 + t235 * t237 + t238 * t251) t12 * t157 + t361 * t380, t12 * t284 - t13 * t157 - t360 * t380 + t361 * t60, -t355 - t373, t295 - t356, -t183 * t202 (t114 * t257 - t115 * t253) * t142 + t167 * t13 - t22 * t284 - t6 * t202 - t350 * t60 + t360 * t54 + (t253 * t290 - t257 * t291) * t183 + t267 * t230 -(t114 * t253 + t115 * t257) * t142 + t167 * t12 + t22 * t157 + t7 * t202 + t350 * t380 + t361 * t54 + (t253 * t291 + t257 * t290) * t183 - t267 * t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199 ^ 2 - t202 ^ 2, t152 * t202 - t153 * t199 + t270 - t376, 0, 0, 0, 0, 0, t278 - t346, -t189 ^ 2 * t258 - t333 - t344, -t108 * t383 + t213 * t49 + t282 * t50 - t285 * t382, -t104 * t202 + t213 * t5 - t282 * t4 - t29 * t383 - t30 * t382 - t376, 0, 0, 0, 0, 0, t295 + t356, -t355 + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170 * t168, -t168 ^ 2 + t170 ^ 2, t100 + t347, -t101 + t345, t149, -t139 * t170 + t189 * t81 + (t299 + t364) * t254 + t292 + t374, g(1) * t193 - g(2) * t191 + t139 * t168 + t189 * t80 + t258 * t364 - t273 (t247 * t49 - t249 * t50) * pkin(4) + (-t33 + t29) * t285 + (-t30 - t32) * t108, -t29 * t32 - t30 * t33 + (-t104 * t170 + t5 * t247 + t4 * t249 + t254 * t364 + t374) * pkin(4), -t394, t393, t395, t390, t142, t276 * t142 - (t23 * t257 - t24 * t253) * t183 + t90 * t60 + (-t183 * t277 - t7) * qJD(6) + t385, -t277 * t142 - t257 * t3 + (t23 * t253 + t24 * t257) * t183 - t90 * t380 + (-t183 * t276 - t353) * qJD(6) + t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 ^ 2 - t285 ^ 2, -t108 * t29 - t285 * t30 - t267 + t45, 0, 0, 0, 0, 0, t13 + t358, t12 + t357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t394, t393, t395, t390, t142 (-qJD(6) + t183) * t7 + t385, t183 * t6 - t257 * t302 + t392;];
tau_reg  = t1;
