% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:11
% EndTime: 2019-03-09 15:30:21
% DurationCPUTime: 3.75s
% Computational Cost: add. (4384->455), mult. (9538->513), div. (0->0), fcn. (6650->10), ass. (0->236)
t187 = sin(qJ(3));
t188 = sin(qJ(2));
t191 = cos(qJ(2));
t316 = cos(qJ(3));
t116 = t187 * t191 + t316 * t188;
t103 = t116 * qJD(1);
t329 = qJD(6) + t103;
t259 = t316 * t191;
t238 = qJD(1) * t259;
t272 = qJD(1) * t188;
t258 = t187 * t272;
t101 = -t238 + t258;
t177 = qJD(2) + qJD(3);
t186 = sin(qJ(6));
t190 = cos(qJ(6));
t76 = t101 * t186 + t190 * t177;
t332 = t329 * t76;
t242 = t190 * t329;
t172 = t191 * pkin(2);
t163 = t172 + pkin(1);
t265 = t191 * qJDD(1);
t267 = qJD(1) * qJD(2);
t331 = t188 * t267 - t265;
t193 = -pkin(8) - pkin(7);
t125 = t193 * t191;
t119 = qJD(1) * t125;
t105 = t187 * t119;
t124 = t193 * t188;
t117 = qJD(1) * t124;
t308 = qJD(2) * pkin(2);
t109 = t117 + t308;
t295 = -t316 * t109 - t105;
t330 = qJD(4) + t295;
t182 = qJ(2) + qJ(3);
t170 = sin(t182);
t159 = g(3) * t170;
t171 = cos(t182);
t192 = cos(qJ(1));
t283 = t171 * t192;
t189 = sin(qJ(1));
t284 = t171 * t189;
t328 = g(1) * t283 + g(2) * t284 + t159;
t317 = -pkin(4) - pkin(9);
t179 = -pkin(3) + t317;
t59 = t103 * pkin(3) + t101 * qJ(4);
t26 = -pkin(5) * t101 + t317 * t103 - t59;
t327 = (qJD(6) * t179 + t26) * t329;
t162 = -t316 * pkin(2) - pkin(3);
t154 = -pkin(4) + t162;
t142 = -pkin(9) + t154;
t264 = pkin(2) * t272;
t326 = (qJD(6) * t142 + t26 - t264) * t329;
t256 = qJD(3) * t316;
t69 = t316 * t117 + t105;
t293 = pkin(2) * t256 + qJD(4) - t69;
t271 = qJD(3) * t187;
t108 = t316 * t119;
t68 = t187 * t117 - t108;
t236 = pkin(2) * t271 - t68;
t100 = t103 ^ 2;
t325 = -t177 ^ 2 - t100;
t275 = t171 * pkin(3) + t170 * qJ(4);
t176 = qJDD(2) + qJDD(3);
t165 = t176 * qJ(4);
t167 = t177 * qJD(4);
t324 = -t165 - t167;
t253 = g(1) * t189 - g(2) * t192;
t175 = g(1) * t192;
t274 = g(2) * t189 + t175;
t168 = t176 * pkin(3);
t322 = qJDD(4) - t168;
t280 = t187 * t188;
t231 = t177 * t280;
t250 = qJDD(1) * t316;
t240 = -t177 * t238 - t187 * t265 - t188 * t250;
t47 = qJD(1) * t231 + t240;
t42 = -qJDD(6) + t47;
t299 = t190 * t42;
t300 = t186 * t329;
t214 = t300 * t329 + t299;
t115 = -t259 + t280;
t260 = qJD(2) * t193;
t118 = t188 * t260;
t120 = t191 * t260;
t80 = t187 * t124 - t316 * t125;
t36 = t80 * qJD(3) + t187 * t118 - t316 * t120;
t70 = -qJD(2) * t259 - t191 * t256 + t231;
t21 = t70 * qJ(5) - t116 * qJD(5) + t36;
t254 = t191 * t267;
t266 = t188 * qJDD(1);
t74 = qJDD(2) * pkin(2) - t193 * (-t254 - t266);
t75 = t193 * t331;
t249 = t109 * t271 - t119 * t256 + t187 * t75 - t316 * t74;
t19 = t249 + t322;
t309 = qJ(5) * t47;
t201 = -qJD(5) * t103 + t19 + t309;
t123 = t163 * qJD(1);
t53 = t101 * pkin(3) - t103 * qJ(4) - t123;
t230 = qJD(5) - t53;
t22 = pkin(5) * t103 + t317 * t101 + t230;
t252 = qJD(6) * t22 + t317 * t176 + t201;
t241 = -t115 * pkin(3) + t116 * qJ(4) + t163;
t28 = pkin(5) * t116 + t317 * t115 + t241;
t169 = t177 * qJ(4);
t292 = t101 * qJ(5);
t62 = t187 * t109 - t108;
t46 = t292 + t62;
t38 = -t169 - t46;
t34 = pkin(5) * t177 - t38;
t79 = -t316 * t124 - t187 * t125;
t56 = -t116 * qJ(5) + t79;
t248 = -t109 * t256 - t119 * t271 - t187 * t74 - t316 * t75;
t227 = t248 + t324;
t229 = t187 * t266 - t191 * t250;
t71 = t177 * t116;
t48 = qJD(1) * t71 + t229;
t9 = -t48 * qJ(5) - t101 * qJD(5) + t227;
t7 = pkin(5) * t176 - t9;
t320 = t7 * t115 - t116 * t252 - (qJD(6) * t28 + t21) * t329 + t34 * t71 + t56 * t42;
t318 = pkin(3) + pkin(4);
t315 = pkin(2) * t188;
t314 = pkin(4) * t103;
t313 = pkin(4) * t176;
t160 = g(3) * t171;
t157 = t171 * pkin(4);
t311 = t28 * t42;
t310 = qJ(4) * t48;
t307 = t101 * t76;
t78 = t101 * t190 - t177 * t186;
t306 = t101 * t78;
t156 = pkin(2) * t187 + qJ(4);
t305 = t156 * t48;
t304 = t177 * t295;
t303 = t177 * t62;
t269 = qJD(6) * t190;
t270 = qJD(6) * t186;
t24 = -t101 * t270 - t186 * t176 - t177 * t269 + t190 * t48;
t302 = t186 * t24;
t301 = t186 * t42;
t298 = t190 * t78;
t297 = t34 * t115;
t296 = t329 * t101;
t290 = t103 * qJ(5);
t294 = -t290 + t293;
t183 = qJDD(1) * pkin(1);
t291 = t101 * t177;
t289 = t103 * t101;
t287 = t156 * t176;
t286 = t170 * t189;
t285 = t170 * t192;
t282 = t186 * t189;
t281 = t186 * t192;
t279 = t189 * t190;
t278 = t190 * t192;
t277 = qJ(5) + t193;
t224 = -t290 + t330;
t180 = t188 ^ 2;
t273 = -t191 ^ 2 + t180;
t263 = t188 * t308;
t96 = t331 * pkin(2) - t183;
t261 = t172 + t275;
t31 = t177 * t179 + t224;
t11 = -t186 * t31 + t190 * t22;
t257 = t11 * t101 + t7 * t190;
t16 = t48 * pkin(3) + t47 * qJ(4) - t103 * qJD(4) + t96;
t217 = qJDD(5) - t16;
t3 = -pkin(5) * t47 + t317 * t48 + t217;
t251 = qJD(6) * t31 - t3;
t246 = t190 * t176 + t186 * t48;
t245 = t329 * t34;
t239 = g(2) * (pkin(3) * t283 + qJ(4) * t285 + t192 * t163);
t237 = -t292 + t236;
t235 = g(1) * t286 - g(2) * t285;
t234 = -g(1) * t284 + g(2) * t283;
t233 = -pkin(3) * t170 - t315;
t228 = t176 - t289;
t55 = t264 + t59;
t12 = t186 * t22 + t190 * t31;
t226 = -t12 * t101 + t328 * t186;
t225 = -t252 - t160;
t223 = -t163 - t275;
t222 = -0.2e1 * pkin(1) * t267 - pkin(7) * qJDD(2);
t27 = t71 * pkin(3) + t70 * qJ(4) - t116 * qJD(4) + t263;
t221 = -g(1) * t285 - g(2) * t286 + t160 + t249;
t220 = -t248 - t328;
t35 = t316 * t118 + t187 * t120 + t124 * t256 + t125 * t271;
t215 = -t242 * t329 + t301;
t213 = -t221 - t322;
t212 = t179 * t42 + t329 * t46 - t245;
t211 = -t101 * t53 + t220;
t210 = -t171 * t274 - t159;
t4 = -pkin(4) * t48 + t217;
t209 = -t123 * t101 - t220;
t208 = t123 * t103 - t221;
t207 = t176 * t80 + t177 * t35 + t235;
t206 = -t176 * t79 - t177 * t36 - t234;
t195 = qJD(2) ^ 2;
t205 = -pkin(7) * t195 + 0.2e1 * t183 + t253;
t196 = qJD(1) ^ 2;
t204 = pkin(1) * t196 - pkin(7) * qJDD(1) + t274;
t203 = t274 * t318 * t170;
t202 = -t177 * t258 - t240;
t200 = -t103 * t53 + t213;
t33 = -pkin(4) * t101 + t230;
t199 = t101 * t33 - t328 - t9;
t198 = t142 * t42 - t237 * t329 - t245;
t197 = t309 + (-qJD(5) - t33) * t103 - t213;
t185 = qJ(4) + pkin(5);
t151 = pkin(5) + t156;
t128 = qJ(4) * t283;
t126 = qJ(4) * t284;
t99 = t101 ^ 2;
t92 = t170 * t278 - t282;
t91 = -t170 * t281 - t279;
t90 = -t170 * t279 - t281;
t89 = t170 * t282 - t278;
t58 = t169 + t62;
t57 = t115 * qJ(5) + t80;
t54 = -pkin(3) * t177 + t330;
t52 = -t99 + t100;
t49 = -pkin(4) * t115 + t241;
t44 = -t59 - t314;
t37 = -t55 - t314;
t32 = -t318 * t177 + t224;
t29 = t202 + t291;
t25 = qJD(6) * t78 + t246;
t20 = -t71 * qJ(5) - t115 * qJD(5) - t35;
t17 = -pkin(4) * t71 - t27;
t15 = t215 + t306;
t14 = t214 - t307;
t13 = -t242 * t78 - t302;
t10 = -pkin(5) * t70 + t317 * t71 - t27;
t8 = t201 - t313;
t2 = t190 * t3;
t1 = (-t24 + t332) * t190 + (t329 * t78 + t25) * t186;
t5 = [qJDD(1), t253, t274, qJDD(1) * t180 + 0.2e1 * t188 * t254, 0.2e1 * t188 * t265 - 0.2e1 * t267 * t273, qJDD(2) * t188 + t191 * t195, qJDD(2) * t191 - t188 * t195, 0, t188 * t222 + t191 * t205, -t188 * t205 + t191 * t222, -t103 * t70 - t116 * t47, t101 * t70 - t103 * t71 + t115 * t47 - t116 * t48, t116 * t176 - t177 * t70, -t115 * t176 - t177 * t71, 0, t101 * t263 + t115 * t96 - t123 * t71 - t163 * t48 + t206, t103 * t263 + t116 * t96 + t123 * t70 + t163 * t47 - t207, t101 * t27 + t115 * t16 - t241 * t48 + t53 * t71 + t206, -t101 * t35 + t103 * t36 + t115 * t227 + t116 * t19 - t47 * t79 - t48 * t80 - t54 * t70 - t58 * t71 - t274, -t103 * t27 - t116 * t16 - t241 * t47 + t53 * t70 + t207, -t227 * t80 + t58 * t35 - t16 * t241 + t53 * t27 + t19 * t79 + t54 * t36 + t193 * t175 - t239 + (-g(1) * t223 + g(2) * t193) * t189, t103 * t17 + t116 * t4 + t176 * t57 - t177 * t20 - t33 * t70 - t47 * t49 + t235, t101 * t17 + t115 * t4 + t176 * t56 + t177 * t21 + t33 * t71 + t48 * t49 + t234, -t101 * t20 - t103 * t21 - t115 * t9 - t116 * t8 + t32 * t70 - t38 * t71 + t47 * t56 + t48 * t57 + t274, t8 * t56 + t32 * t21 - t9 * t57 + t38 * t20 + t4 * t49 + t33 * t17 - t239 + (g(1) * t277 - g(2) * t157) * t192 + (-g(1) * (t223 - t157) + g(2) * t277) * t189, t71 * t298 + (t190 * t24 - t270 * t78) * t115 (-t186 * t78 - t190 * t76) * t71 + (-t302 - t190 * t25 + (t186 * t76 - t298) * qJD(6)) * t115, t71 * t242 + t116 * t24 - t70 * t78 + (-t270 * t329 - t299) * t115, -t71 * t300 - t116 * t25 + t70 * t76 + (-t269 * t329 + t301) * t115, -t116 * t42 - t329 * t70, -g(1) * t90 - g(2) * t92 - t11 * t70 + t2 * t116 - t20 * t76 + t57 * t25 + (t10 * t329 - t311 + (-t31 * t116 - t329 * t56 + t297) * qJD(6)) * t190 + t320 * t186, -g(1) * t89 - g(2) * t91 + t12 * t70 - t20 * t78 + t57 * t24 + (-(-qJD(6) * t56 + t10) * t329 + t311 + t251 * t116 - qJD(6) * t297) * t186 + t320 * t190; 0, 0, 0, -t188 * t196 * t191, t273 * t196, t266, t265, qJDD(2), -g(3) * t191 + t188 * t204, g(3) * t188 + t191 * t204, t289, t52, t29, -t229, t176, t68 * t177 + (-t101 * t272 + t316 * t176 - t177 * t271) * pkin(2) + t208, t69 * t177 + (-t103 * t272 - t176 * t187 - t177 * t256) * pkin(2) + t209, -t101 * t55 - t162 * t176 - t177 * t236 + t200, -t305 - t162 * t47 + (t236 + t58) * t103 + (t54 - t293) * t101, t103 * t55 + t177 * t293 + t211 + t287 - t324, -t227 * t156 + t19 * t162 - t53 * t55 - g(1) * (t192 * t233 + t128) - g(2) * (t189 * t233 + t126) - g(3) * t261 + t293 * t58 + t236 * t54, -t103 * t37 + t177 * t294 + t199 + t287 (-pkin(4) + t154) * t176 - t101 * t37 + t237 * t177 + t197, t154 * t47 + t305 + (-t237 + t38) * t103 + (-t32 + t294) * t101, t8 * t154 - t9 * t156 - t33 * t37 - g(1) * (-t192 * t315 + t128) - g(2) * (-t189 * t315 + t126) - g(3) * (t157 + t261) - t294 * t38 + t237 * t32 + t203, t13, t1, t15, t14, t296, t151 * t25 + t294 * t76 + (t210 - t326) * t190 + t198 * t186 + t257, t151 * t24 + t294 * t78 + (-t7 + t326) * t186 + t198 * t190 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, t52, t29, -t229, t176, t208 + t303, t209 - t304, -t101 * t59 + t168 + t200 + t303, pkin(3) * t47 - t310 + (t58 - t62) * t103 + (t54 - t330) * t101, t103 * t59 + 0.2e1 * t165 + 0.2e1 * t167 + t211 + t304, -t227 * qJ(4) - t19 * pkin(3) - t53 * t59 - t54 * t62 - g(1) * (-pkin(3) * t285 + t128) - g(2) * (-pkin(3) * t286 + t126) - g(3) * t275 + t330 * t58, -t103 * t44 + t177 * t224 + t165 + t199 (-pkin(4) - t318) * t176 + t197 - t101 * t44 - t177 * t46, t310 - t318 * t47 + (t38 + t46) * t103 + (-t32 + t224) * t101, -t8 * t318 - t9 * qJ(4) - t32 * t46 - t33 * t44 - g(1) * t128 - g(2) * t126 - g(3) * (t157 + t275) - t224 * t38 + t203, t13, t1, t15, t14, t296, t185 * t25 + t224 * t76 + t212 * t186 + (t210 - t327) * t190 + t257, t185 * t24 + t224 * t78 + (-t7 + t327) * t186 + t212 * t190 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, t29, t325, -t177 * t58 - t200, t325, t228, t47 - t291, t177 * t38 + t197 - t313, 0, 0, 0, 0, 0, -t177 * t76 + t215, -t177 * t78 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202 - t291, t103 * t177 + t48, -t100 - t99, t101 * t38 + t103 * t32 + t253 + t4, 0, 0, 0, 0, 0, -t214 - t307, t215 - t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t24 + t332, -t246 + (-qJD(6) + t329) * t78, -t42, -g(1) * t91 + g(2) * t89 + t12 * t329 + t186 * t225 - t269 * t31 - t34 * t78 + t2, g(1) * t92 - g(2) * t90 + t11 * t329 + t186 * t251 + t190 * t225 + t34 * t76;];
tau_reg  = t5;
