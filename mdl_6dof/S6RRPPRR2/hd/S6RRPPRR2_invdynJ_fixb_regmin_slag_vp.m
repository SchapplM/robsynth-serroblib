% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:40
% EndTime: 2019-03-09 08:52:55
% DurationCPUTime: 6.73s
% Computational Cost: add. (8871->499), mult. (21147->670), div. (0->0), fcn. (16901->18), ass. (0->239)
t252 = cos(qJ(2));
t326 = cos(pkin(10));
t286 = t326 * t252;
t220 = qJD(1) * t286;
t243 = sin(pkin(10));
t248 = sin(qJ(2));
t308 = qJD(1) * t248;
t190 = t243 * t308 - t220;
t185 = qJD(5) + t190;
t179 = qJD(6) + t185;
t206 = t243 * t252 + t326 * t248;
t193 = t206 * qJD(1);
t242 = sin(pkin(11));
t244 = cos(pkin(11));
t163 = qJD(2) * t242 + t193 * t244;
t164 = t244 * qJD(2) - t193 * t242;
t247 = sin(qJ(5));
t251 = cos(qJ(5));
t101 = t163 * t247 - t164 * t251;
t246 = sin(qJ(6));
t250 = cos(qJ(6));
t350 = -t163 * t251 - t164 * t247;
t47 = t250 * t101 - t246 * t350;
t363 = t179 * t47;
t362 = t101 * t185;
t273 = t246 * t101 + t250 * t350;
t361 = t179 * t273;
t205 = t242 * t247 - t251 * t244;
t311 = t185 * t205;
t207 = t242 * t251 + t244 * t247;
t197 = t207 * qJD(5);
t310 = t207 * t190 + t197;
t319 = t190 * t242;
t128 = pkin(2) * t308 + pkin(3) * t193 + qJ(4) * t190;
t245 = -qJ(3) - pkin(7);
t216 = t245 * t252;
t212 = qJD(1) * t216;
t198 = t243 * t212;
t215 = t245 * t248;
t211 = qJD(1) * t215;
t154 = t326 * t211 + t198;
t83 = t242 * t128 + t244 * t154;
t65 = pkin(8) * t319 + t83;
t360 = -t244 * qJD(4) + t65;
t239 = qJ(2) + pkin(10);
t233 = sin(t239);
t235 = cos(t239);
t249 = sin(qJ(1));
t253 = cos(qJ(1));
t280 = g(1) * t253 + g(2) * t249;
t263 = -g(3) * t235 + t280 * t233;
t288 = qJD(2) * t245;
t187 = -qJD(3) * t248 + t252 * t288;
t141 = qJDD(2) * pkin(2) + t187 * qJD(1) + qJDD(1) * t215;
t186 = qJD(3) * t252 + t248 * t288;
t152 = t186 * qJD(1) - qJDD(1) * t216;
t85 = t326 * t141 - t243 * t152;
t84 = -qJDD(2) * pkin(3) + qJDD(4) - t85;
t349 = t84 - t263;
t359 = t185 * t350;
t151 = -t205 * t246 + t207 * t250;
t331 = qJD(6) * t151 - t311 * t246 + t310 * t250;
t356 = t273 * t47;
t355 = t273 ^ 2 - t47 ^ 2;
t304 = qJD(6) * t250;
t305 = qJD(6) * t246;
t301 = qJD(1) * qJD(2);
t295 = t248 * t301;
t148 = qJD(2) * t220 + t206 * qJDD(1) - t243 * t295;
t124 = -t244 * qJDD(2) + t148 * t242;
t125 = qJDD(2) * t242 + t148 * t244;
t306 = qJD(5) * t251;
t307 = qJD(5) * t247;
t33 = -t247 * t124 + t251 * t125 - t163 * t307 + t164 * t306;
t34 = -qJD(5) * t350 + t251 * t124 + t125 * t247;
t8 = -t101 * t304 - t246 * t34 + t250 * t33 + t305 * t350;
t354 = t8 + t363;
t343 = pkin(2) * t252;
t229 = pkin(1) + t343;
t214 = -t229 * qJD(1) + qJD(3);
t117 = pkin(3) * t190 - qJ(4) * t193 + t214;
t330 = qJD(2) * pkin(2);
t203 = t211 + t330;
t287 = t326 * t212;
t146 = t243 * t203 - t287;
t137 = qJD(2) * qJ(4) + t146;
t70 = t244 * t117 - t137 * t242;
t41 = pkin(4) * t190 - pkin(8) * t163 + t70;
t71 = t242 * t117 + t244 * t137;
t54 = pkin(8) * t164 + t71;
t23 = t247 * t41 + t251 * t54;
t13 = -pkin(9) * t101 + t23;
t11 = t13 * t305;
t238 = pkin(11) + qJ(5);
t236 = qJ(6) + t238;
t226 = sin(t236);
t227 = cos(t236);
t315 = t235 * t249;
t166 = t226 * t253 - t227 * t315;
t314 = t235 * t253;
t168 = t226 * t249 + t227 * t314;
t338 = g(3) * t233;
t145 = t326 * t203 + t198;
t131 = -qJD(2) * pkin(3) + qJD(4) - t145;
t94 = -pkin(4) * t164 + t131;
t45 = t101 * pkin(5) + t94;
t353 = g(1) * t168 - g(2) * t166 + t227 * t338 + t45 * t47 + t11;
t165 = t226 * t315 + t227 * t253;
t167 = -t226 * t314 + t227 * t249;
t192 = t206 * qJD(2);
t300 = t248 * qJDD(1);
t147 = qJD(1) * t192 - qJDD(1) * t286 + t243 * t300;
t142 = qJDD(5) + t147;
t264 = pkin(2) * t295 - t229 * qJDD(1) + qJDD(3);
t57 = pkin(3) * t147 - qJ(4) * t148 - qJD(4) * t193 + t264;
t86 = t243 * t141 + t326 * t152;
t81 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t86;
t28 = -t242 * t81 + t244 * t57;
t18 = pkin(4) * t147 - pkin(8) * t125 + t28;
t29 = t242 * t57 + t244 * t81;
t26 = -pkin(8) * t124 + t29;
t292 = t251 * t18 - t247 * t26;
t258 = -t23 * qJD(5) + t292;
t2 = pkin(5) * t142 - pkin(9) * t33 + t258;
t269 = t247 * t18 + t251 * t26 + t41 * t306 - t54 * t307;
t3 = -pkin(9) * t34 + t269;
t297 = t250 * t2 - t246 * t3;
t22 = -t247 * t54 + t251 * t41;
t12 = pkin(9) * t350 + t22;
t10 = pkin(5) * t185 + t12;
t329 = t13 * t250;
t5 = t10 * t246 + t329;
t352 = -g(1) * t167 + g(2) * t165 - t5 * qJD(6) + t226 * t338 + t45 * t273 + t297;
t257 = t273 * qJD(6) - t246 * t33 - t250 * t34;
t351 = t257 - t361;
t347 = g(1) * t249 - g(2) * t253;
t348 = t235 * t347;
t223 = pkin(2) * t243 + qJ(4);
t335 = pkin(8) + t223;
t201 = t335 * t242;
t202 = t335 * t244;
t312 = -t247 * t201 + t251 * t202;
t138 = qJDD(6) + t142;
t150 = t250 * t205 + t207 * t246;
t332 = -t150 * qJD(6) - t246 * t310 - t250 * t311;
t346 = -t138 * t151 - t179 * t332;
t345 = -t142 * t207 + t185 * t311;
t271 = t242 * qJD(4) + qJD(5) * t202;
t342 = pkin(8) * t244;
t82 = t244 * t128 - t154 * t242;
t53 = pkin(4) * t193 + t190 * t342 + t82;
t344 = t201 * t306 + t360 * t251 + (t271 + t53) * t247;
t188 = t190 ^ 2;
t336 = g(3) * t252;
t266 = -t243 * t248 + t286;
t316 = t206 * t244;
t144 = -pkin(3) * t266 - qJ(4) * t206 - t229;
t159 = t243 * t215 - t326 * t216;
t92 = t244 * t144 - t159 * t242;
t64 = -pkin(4) * t266 - pkin(8) * t316 + t92;
t317 = t206 * t242;
t93 = t242 * t144 + t244 * t159;
t74 = -pkin(8) * t317 + t93;
t333 = t247 * t64 + t251 * t74;
t328 = t193 * t47;
t327 = t193 * t273;
t195 = t266 * qJD(2);
t298 = t248 * t330;
t106 = pkin(3) * t192 - qJ(4) * t195 - qJD(4) * t206 + t298;
t130 = t326 * t186 + t243 * t187;
t60 = t242 * t106 + t244 * t130;
t325 = t101 * t193;
t324 = t350 * t193;
t321 = t147 * t242;
t320 = t147 * t244;
t318 = t195 * t242;
t240 = t248 ^ 2;
t309 = -t252 ^ 2 + t240;
t302 = -qJD(4) + t131;
t299 = t252 * qJDD(1);
t153 = t211 * t243 - t287;
t107 = -pkin(4) * t319 + t153;
t296 = pkin(5) * t310 - t107;
t293 = qJD(6) * t10 + t3;
t59 = t244 * t106 - t130 * t242;
t38 = pkin(4) * t192 - t195 * t342 + t59;
t43 = -pkin(8) * t318 + t60;
t290 = -t247 * t43 + t251 * t38;
t289 = -t247 * t74 + t251 * t64;
t129 = t186 * t243 - t326 * t187;
t284 = -t251 * t201 - t202 * t247;
t158 = -t326 * t215 - t216 * t243;
t283 = -t150 * t138 - t179 * t331;
t228 = -t326 * pkin(2) - pkin(3);
t108 = -pkin(9) * t207 + t284;
t282 = pkin(9) * t310 - qJD(6) * t108 + t344;
t109 = -pkin(9) * t205 + t312;
t52 = t251 * t53;
t281 = pkin(5) * t193 - pkin(9) * t311 + t207 * qJD(4) + t312 * qJD(5) + qJD(6) * t109 - t247 * t65 + t52;
t278 = -t205 * t142 - t185 * t310;
t96 = pkin(4) * t318 + t129;
t123 = pkin(4) * t317 + t158;
t277 = pkin(3) * t235 + qJ(4) * t233;
t276 = -t242 * t29 - t244 * t28;
t275 = -t70 * t242 + t71 * t244;
t132 = t207 * t206;
t133 = t205 * t206;
t78 = t250 * t132 - t133 * t246;
t79 = -t132 * t246 - t133 * t250;
t270 = -0.2e1 * pkin(1) * t301 - pkin(7) * qJDD(2);
t268 = t247 * t38 + t251 * t43 + t64 * t306 - t74 * t307;
t213 = -t244 * pkin(4) + t228;
t254 = qJD(2) ^ 2;
t261 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t254 + t347;
t255 = qJD(1) ^ 2;
t260 = pkin(1) * t255 - pkin(7) * qJDD(1) + t280;
t259 = t131 * t195 + t84 * t206 - t280;
t44 = t124 * pkin(4) + t84;
t234 = cos(t238);
t232 = sin(t238);
t218 = t253 * t229;
t172 = t232 * t249 + t234 * t314;
t171 = -t232 * t314 + t234 * t249;
t170 = t232 * t253 - t234 * t315;
t169 = t232 * t315 + t234 * t253;
t160 = t205 * pkin(5) + t213;
t80 = pkin(5) * t132 + t123;
t77 = t195 * t207 + t306 * t316 - t307 * t317;
t76 = -t195 * t205 - t197 * t206;
t35 = pkin(5) * t77 + t96;
t27 = -pkin(9) * t132 + t333;
t24 = -pkin(5) * t266 + pkin(9) * t133 + t289;
t21 = t79 * qJD(6) + t246 * t76 + t250 * t77;
t20 = -t78 * qJD(6) - t246 * t77 + t250 * t76;
t17 = t34 * pkin(5) + t44;
t7 = -pkin(9) * t77 + t268;
t6 = pkin(5) * t192 - pkin(9) * t76 - qJD(5) * t333 + t290;
t4 = t10 * t250 - t13 * t246;
t1 = [qJDD(1), t347, t280, qJDD(1) * t240 + 0.2e1 * t252 * t295, 0.2e1 * t248 * t299 - 0.2e1 * t309 * t301, qJDD(2) * t248 + t252 * t254, qJDD(2) * t252 - t248 * t254, 0, t248 * t270 + t252 * t261, -t248 * t261 + t252 * t270, t129 * t193 - t130 * t190 - t145 * t195 - t146 * t192 - t147 * t159 + t148 * t158 - t206 * t85 + t266 * t86 - t280, t86 * t159 + t146 * t130 - t85 * t158 - t145 * t129 - t264 * t229 + t214 * t298 - g(1) * (-t229 * t249 - t245 * t253) - g(2) * (-t245 * t249 + t218) t158 * t124 - t129 * t164 + t92 * t147 + t59 * t190 + t70 * t192 + t259 * t242 + t244 * t348 - t266 * t28, t158 * t125 + t129 * t163 - t93 * t147 - t60 * t190 - t71 * t192 - t242 * t348 + t259 * t244 + t266 * t29, -t124 * t93 - t125 * t92 + t164 * t60 - t163 * t59 + t347 * t233 + t276 * t206 + (-t242 * t71 - t244 * t70) * t195, -g(2) * t218 + t131 * t129 + t84 * t158 + t28 * t92 + t29 * t93 + t70 * t59 + t71 * t60 + (g(1) * t245 - g(2) * t277) * t253 + (-g(1) * (-t229 - t277) + g(2) * t245) * t249, -t133 * t33 - t350 * t76, -t101 * t76 - t132 * t33 + t133 * t34 + t350 * t77, -t133 * t142 + t185 * t76 - t192 * t350 - t266 * t33, -t101 * t192 - t132 * t142 - t185 * t77 + t266 * t34, -t142 * t266 + t185 * t192, t290 * t185 + t289 * t142 - t292 * t266 + t22 * t192 + t96 * t101 + t123 * t34 + t44 * t132 + t94 * t77 - g(1) * t170 - g(2) * t172 + (-t185 * t333 + t23 * t266) * qJD(5), -g(1) * t169 - g(2) * t171 + t123 * t33 - t44 * t133 - t142 * t333 - t185 * t268 - t23 * t192 + t266 * t269 - t350 * t96 + t94 * t76, -t20 * t273 + t79 * t8, -t20 * t47 + t21 * t273 + t257 * t79 - t78 * t8, t138 * t79 + t179 * t20 - t192 * t273 - t266 * t8, -t138 * t78 - t179 * t21 - t192 * t47 - t257 * t266, -t138 * t266 + t179 * t192 (-t246 * t7 + t250 * t6) * t179 + (t24 * t250 - t246 * t27) * t138 - t297 * t266 + t4 * t192 + t35 * t47 - t80 * t257 + t17 * t78 + t45 * t21 - g(1) * t166 - g(2) * t168 + ((-t24 * t246 - t250 * t27) * t179 + t5 * t266) * qJD(6), -g(1) * t165 - g(2) * t167 - t11 * t266 + t17 * t79 - t5 * t192 + t45 * t20 - t35 * t273 + t80 * t8 + (-(-qJD(6) * t27 + t6) * t179 - t24 * t138 + t2 * t266) * t246 + (-(qJD(6) * t24 + t7) * t179 - t27 * t138 + t293 * t266) * t250; 0, 0, 0, -t248 * t255 * t252, t309 * t255, t300, t299, qJDD(2), t248 * t260 - t336, g(3) * t248 + t252 * t260 (t146 - t153) * t193 + (-t145 + t154) * t190 + (-t147 * t243 - t326 * t148) * pkin(2), t145 * t153 - t146 * t154 + (t326 * t85 - t336 + t243 * t86 + (-qJD(1) * t214 + t280) * t248) * pkin(2), -t223 * t321 + t124 * t228 + t153 * t164 - t193 * t70 + (t302 * t242 - t82) * t190 - t349 * t244, -t223 * t320 + t125 * t228 - t153 * t163 + t193 * t71 + (t302 * t244 + t83) * t190 + t349 * t242, -t338 - t164 * t83 + t163 * t82 - t280 * t235 + (qJD(4) * t164 - t124 * t223 - t190 * t70 + t29) * t244 + (qJD(4) * t163 + t125 * t223 - t190 * t71 - t28) * t242, t84 * t228 - t71 * t83 - t70 * t82 - t131 * t153 - g(3) * (t277 + t343) + (-t28 * t242 + t29 * t244) * t223 + t275 * qJD(4) + t280 * (pkin(2) * t248 + pkin(3) * t233 - qJ(4) * t235) t207 * t33 + t311 * t350, t101 * t311 - t205 * t33 - t207 * t34 + t310 * t350, t324 - t345, t278 + t325, -t185 * t193, t284 * t142 + t213 * t34 + t44 * t205 - t22 * t193 - t107 * t101 + t310 * t94 + (-t52 - t271 * t251 + (qJD(5) * t201 + t360) * t247) * t185 + t263 * t234, t107 * t350 - t312 * t142 + t185 * t344 + t23 * t193 + t44 * t207 + t213 * t33 - t232 * t263 - t311 * t94, t151 * t8 - t273 * t332, -t150 * t8 + t151 * t257 + t273 * t331 - t332 * t47, t327 - t346, t283 + t328, -t179 * t193 (t108 * t250 - t109 * t246) * t138 - t160 * t257 + t17 * t150 - t4 * t193 + t296 * t47 + t331 * t45 + (t246 * t282 - t250 * t281) * t179 + t263 * t227 -(t108 * t246 + t109 * t250) * t138 + t160 * t8 + t17 * t151 + t5 * t193 - t296 * t273 + t332 * t45 + (t246 * t281 + t250 * t282) * t179 - t263 * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193 ^ 2 - t188, t145 * t193 + t146 * t190 + t264 - t347, t164 * t193 - t188 * t242 + t320, -t163 * t193 - t188 * t244 - t321, -t124 * t242 - t125 * t244 + (t163 * t242 + t164 * t244) * t190, -t131 * t193 + t190 * t275 - t276 - t347, 0, 0, 0, 0, 0, t278 - t325, t324 + t345, 0, 0, 0, 0, 0, t283 - t328, t327 + t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163 * t190 + t124, t164 * t190 + t125, -t163 ^ 2 - t164 ^ 2, t163 * t70 - t71 * t164 + t349, 0, 0, 0, 0, 0, t34 - t359, t33 - t362, 0, 0, 0, 0, 0, -t257 - t361, t8 - t363; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350 * t101, -t101 ^ 2 + t350 ^ 2, t33 + t362, -t34 - t359, t142, -g(1) * t171 + g(2) * t169 + t185 * t23 + t232 * t338 + t350 * t94 + t258, g(1) * t172 - g(2) * t170 + t101 * t94 + t185 * t22 + t234 * t338 - t269, -t356, t355, t354, t351, t138 -(-t12 * t246 - t329) * t179 + (t250 * t138 - t179 * t305 + t350 * t47) * pkin(5) + t352 (-t13 * t179 - t2) * t246 + (t12 * t179 - t293) * t250 + (-t246 * t138 - t179 * t304 - t273 * t350) * pkin(5) + t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t356, t355, t354, t351, t138, t179 * t5 + t352, t179 * t4 - t246 * t2 - t250 * t293 + t353;];
tau_reg  = t1;
