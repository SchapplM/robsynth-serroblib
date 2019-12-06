% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:56
% EndTime: 2019-12-05 18:39:06
% DurationCPUTime: 4.92s
% Computational Cost: add. (9243->434), mult. (22732->555), div. (0->0), fcn. (16954->16), ass. (0->230)
t249 = cos(qJ(2));
t345 = cos(qJ(3));
t292 = t345 * t249;
t280 = qJD(1) * t292;
t245 = sin(qJ(3));
t246 = sin(qJ(2));
t303 = qJD(1) * t246;
t291 = t245 * t303;
t164 = -t280 + t291;
t177 = t245 * t249 + t246 * t345;
t166 = t177 * qJD(1);
t242 = sin(pkin(9));
t243 = cos(pkin(9));
t121 = t164 * t243 + t166 * t242;
t244 = sin(qJ(5));
t248 = cos(qJ(5));
t270 = -t164 * t242 + t166 * t243;
t300 = qJD(5) * t248;
t301 = qJD(5) * t244;
t238 = qJD(2) + qJD(3);
t307 = t245 * t246;
t274 = t238 * t307;
t283 = qJDD(1) * t345;
t297 = t249 * qJDD(1);
t281 = -t238 * t280 - t245 * t297 - t246 * t283;
t104 = qJD(1) * t274 + t281;
t137 = t238 * t177;
t298 = t246 * qJDD(1);
t273 = t245 * t298 - t249 * t283;
t105 = qJD(1) * t137 + t273;
t54 = -t104 * t242 + t105 * t243;
t55 = -t104 * t243 - t105 * t242;
t15 = t121 * t300 + t244 * t54 - t248 * t55 + t270 * t301;
t229 = qJD(5) + t238;
t71 = t121 * t248 + t244 * t270;
t324 = t229 * t71;
t8 = -t15 + t324;
t271 = t121 * t244 - t248 * t270;
t261 = qJD(5) * t271 - t244 * t55 - t248 * t54;
t325 = t229 * t271;
t9 = t261 - t325;
t332 = t71 ^ 2;
t334 = t271 ^ 2;
t12 = -t332 + t334;
t355 = t271 * t71;
t117 = t270 * pkin(8);
t251 = -pkin(7) - pkin(6);
t198 = t251 * t249;
t184 = qJD(1) * t198;
t167 = t245 * t184;
t197 = t251 * t246;
t182 = qJD(1) * t197;
t326 = qJD(2) * pkin(2);
t173 = t182 + t326;
t127 = t173 * t345 + t167;
t160 = t166 * qJ(4);
t102 = t127 - t160;
t92 = pkin(3) * t238 + t102;
t171 = t345 * t184;
t128 = t173 * t245 - t171;
t312 = t164 * qJ(4);
t103 = t128 - t312;
t93 = t242 * t103;
t47 = t243 * t92 - t93;
t35 = pkin(4) * t238 - t117 + t47;
t340 = pkin(8) * t121;
t309 = t243 * t103;
t48 = t242 * t92 + t309;
t38 = t48 - t340;
t236 = qJDD(2) + qJDD(3);
t299 = qJD(1) * qJD(2);
t287 = t249 * t299;
t139 = qJDD(2) * pkin(2) - t251 * (-t287 - t298);
t288 = t246 * t299;
t141 = t251 * (-t288 + t297);
t65 = -qJD(3) * t128 + t139 * t345 + t245 * t141;
t32 = t236 * pkin(3) + t104 * qJ(4) - t166 * qJD(4) + t65;
t289 = qJD(3) * t345;
t302 = qJD(3) * t245;
t282 = -t139 * t245 + t141 * t345 - t173 * t289 - t184 * t302;
t34 = -qJ(4) * t105 - qJD(4) * t164 - t282;
t10 = -t242 * t34 + t243 * t32;
t6 = pkin(4) * t236 - pkin(8) * t55 + t10;
t11 = t242 * t32 + t243 * t34;
t7 = -pkin(8) * t54 + t11;
t1 = (qJD(5) * t35 + t7) * t248 + t244 * t6 - t38 * t301;
t241 = qJ(2) + qJ(3);
t230 = pkin(9) + t241;
t222 = qJ(5) + t230;
t211 = sin(t222);
t212 = cos(t222);
t247 = sin(qJ(1));
t250 = cos(qJ(1));
t277 = g(1) * t250 + g(2) * t247;
t234 = t249 * pkin(2);
t225 = t234 + pkin(1);
t196 = t225 * qJD(1);
t140 = pkin(3) * t164 + qJD(4) - t196;
t85 = pkin(4) * t121 + t140;
t258 = g(3) * t211 + t212 * t277 + t71 * t85 - t1;
t14 = t244 * t35 + t248 * t38;
t2 = -qJD(5) * t14 - t244 * t7 + t248 * t6;
t256 = -g(3) * t212 + t211 * t277 + t271 * t85 + t2;
t134 = -t182 * t245 + t171;
t106 = t134 + t312;
t135 = t182 * t345 + t167;
t107 = -t160 + t135;
t308 = t243 * t245;
t327 = pkin(2) * qJD(3);
t321 = -t243 * t106 + t107 * t242 + (-t242 * t345 - t308) * t327;
t310 = t242 * t245;
t320 = -t242 * t106 - t243 * t107 + (t243 * t345 - t310) * t327;
t231 = sin(t241);
t232 = cos(t241);
t354 = -g(3) * t232 + t231 * t277;
t353 = -t340 + t321;
t352 = t117 + t320;
t351 = t121 * t270;
t143 = t197 * t245 - t198 * t345;
t347 = g(1) * t247 - g(2) * t250;
t346 = pkin(4) * t54;
t344 = pkin(2) * t246;
t343 = pkin(3) * t166;
t342 = pkin(3) * t231;
t341 = pkin(3) * t242;
t335 = g(3) * t249;
t176 = -t292 + t307;
t293 = qJD(2) * t251;
t183 = t246 * t293;
t185 = t249 * t293;
t89 = t183 * t345 + t185 * t245 + t197 * t289 + t198 * t302;
t60 = -qJ(4) * t137 - qJD(4) * t176 + t89;
t136 = -qJD(2) * t292 - t249 * t289 + t274;
t90 = -qJD(3) * t143 - t245 * t183 + t185 * t345;
t61 = t136 * qJ(4) - t177 * qJD(4) + t90;
t28 = t242 * t61 + t243 * t60;
t224 = pkin(2) * t345 + pkin(3);
t156 = -pkin(2) * t310 + t224 * t243;
t152 = pkin(4) + t156;
t158 = pkin(2) * t308 + t224 * t242;
t113 = t152 * t248 - t158 * t244;
t329 = qJD(5) * t113 + t244 * t353 + t248 * t352;
t114 = t152 * t244 + t158 * t248;
t328 = -qJD(5) * t114 - t244 * t352 + t248 * t353;
t53 = t102 * t243 - t93;
t219 = pkin(3) * t243 + pkin(4);
t157 = t219 * t248 - t244 * t341;
t52 = -t102 * t242 - t309;
t41 = t52 + t340;
t42 = -t117 + t53;
t323 = qJD(5) * t157 - t244 * t41 - t248 * t42;
t159 = t219 * t244 + t248 * t341;
t322 = -qJD(5) * t159 + t244 * t42 - t248 * t41;
t319 = pkin(6) * qJDD(1);
t318 = t270 ^ 2;
t317 = t270 * t238;
t315 = t121 ^ 2;
t314 = t121 * t238;
t311 = t166 * t164;
t142 = t197 * t345 + t198 * t245;
t118 = -qJ(4) * t177 + t142;
t119 = -qJ(4) * t176 + t143;
t67 = t118 * t242 + t119 * t243;
t221 = pkin(3) * t232;
t306 = t221 + t234;
t239 = t246 ^ 2;
t240 = t249 ^ 2;
t305 = t239 - t240;
t304 = t239 + t240;
t237 = -qJ(4) + t251;
t227 = t246 * t326;
t253 = qJD(1) ^ 2;
t296 = t246 * t253 * t249;
t218 = cos(t230);
t207 = pkin(4) * t218;
t295 = t207 + t306;
t124 = pkin(3) * t137 + t227;
t27 = -t242 * t60 + t243 * t61;
t66 = t118 * t243 - t119 * t242;
t147 = pkin(3) * t176 - t225;
t279 = t246 * t287;
t91 = pkin(4) * t270 + t343;
t217 = sin(t230);
t278 = -pkin(4) * t217 - t342;
t13 = -t244 * t38 + t248 * t35;
t275 = -t13 * t71 - t14 * t271;
t272 = -t121 * t47 + t270 * t48;
t131 = -t176 * t242 + t177 * t243;
t45 = -pkin(8) * t131 + t66;
t130 = t176 * t243 + t177 * t242;
t46 = -pkin(8) * t130 + t67;
t23 = -t244 * t46 + t248 * t45;
t24 = t244 * t45 + t248 * t46;
t76 = -t130 * t244 + t131 * t248;
t269 = -0.2e1 * pkin(1) * t299 - pkin(6) * qJDD(2);
t161 = pkin(2) * t288 - qJDD(1) * t225;
t266 = g(3) * t217 + t121 * t140 + t218 * t277 - t11;
t252 = qJD(2) ^ 2;
t265 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t252 + t347;
t264 = pkin(1) * t253 + t277 - t319;
t263 = g(3) * t231 - t196 * t164 + t232 * t277 + t282;
t80 = pkin(3) * t105 + qJDD(4) + t161;
t262 = -g(3) * t218 - t140 * t270 + t217 * t277 + t10;
t259 = -t347 + t80;
t255 = t196 * t166 + t354 + t65;
t233 = -pkin(8) + t237;
t228 = qJDD(5) + t236;
t226 = pkin(2) * t303;
t181 = pkin(1) + t306;
t148 = pkin(1) + t295;
t144 = t226 + t343;
t108 = -t164 ^ 2 + t166 ^ 2;
t97 = pkin(4) * t130 + t147;
t86 = t226 + t91;
t84 = -t136 * t243 - t137 * t242;
t83 = -t136 * t242 + t137 * t243;
t81 = -t281 + (t164 - t291) * t238;
t75 = t130 * t248 + t131 * t244;
t59 = pkin(4) * t83 + t124;
t40 = t55 + t314;
t39 = -t54 + t317;
t37 = -t315 + t318;
t29 = t80 + t346;
t26 = qJD(5) * t76 + t244 * t84 + t248 * t83;
t25 = t130 * t300 + t131 * t301 + t244 * t83 - t248 * t84;
t22 = -pkin(8) * t83 + t28;
t21 = -pkin(8) * t84 + t27;
t4 = -qJD(5) * t24 + t21 * t248 - t22 * t244;
t3 = qJD(5) * t23 + t21 * t244 + t22 * t248;
t5 = [0, 0, 0, 0, 0, qJDD(1), t347, t277, 0, 0, qJDD(1) * t239 + 0.2e1 * t279, 0.2e1 * t246 * t297 - 0.2e1 * t299 * t305, qJDD(2) * t246 + t249 * t252, qJDD(1) * t240 - 0.2e1 * t279, qJDD(2) * t249 - t246 * t252, 0, t246 * t269 + t249 * t265, -t246 * t265 + t249 * t269, 0.2e1 * t304 * t319 - t277, -g(1) * (-pkin(1) * t247 + pkin(6) * t250) - g(2) * (pkin(1) * t250 + pkin(6) * t247) + (pkin(6) ^ 2 * t304 + pkin(1) ^ 2) * qJDD(1), -t104 * t177 - t136 * t166, t104 * t176 - t105 * t177 + t136 * t164 - t137 * t166, -t136 * t238 + t177 * t236, t105 * t176 + t137 * t164, -t137 * t238 - t176 * t236, 0, -t105 * t225 - t137 * t196 + t142 * t236 + t161 * t176 + t164 * t227 + t232 * t347 + t238 * t90, t104 * t225 + t136 * t196 - t143 * t236 + t161 * t177 + t166 * t227 - t231 * t347 - t238 * t89, t104 * t142 - t105 * t143 + t127 * t136 - t128 * t137 - t164 * t89 - t166 * t90 + t176 * t282 - t177 * t65 - t277, -t282 * t143 + t128 * t89 + t65 * t142 + t127 * t90 - t161 * t225 - t196 * t227 - g(1) * (-t225 * t247 - t250 * t251) - g(2) * (t225 * t250 - t247 * t251), t131 * t55 + t270 * t84, -t121 * t84 - t130 * t55 - t131 * t54 - t270 * t83, t131 * t236 + t238 * t84, t121 * t83 + t130 * t54, -t130 * t236 - t238 * t83, 0, t121 * t124 + t130 * t80 + t140 * t83 + t147 * t54 + t218 * t347 + t236 * t66 + t238 * t27, t124 * t270 + t131 * t80 + t140 * t84 + t147 * t55 - t217 * t347 - t236 * t67 - t238 * t28, -t10 * t131 - t11 * t130 - t121 * t28 - t27 * t270 - t47 * t84 - t48 * t83 - t54 * t67 - t55 * t66 - t277, t11 * t67 + t48 * t28 + t10 * t66 + t47 * t27 + t80 * t147 + t140 * t124 - g(1) * (-t181 * t247 - t237 * t250) - g(2) * (t181 * t250 - t237 * t247), -t15 * t76 + t25 * t271, t15 * t75 + t25 * t71 + t26 * t271 + t261 * t76, t228 * t76 - t229 * t25, t26 * t71 - t261 * t75, -t228 * t75 - t229 * t26, 0, t212 * t347 + t228 * t23 + t229 * t4 + t26 * t85 - t261 * t97 + t29 * t75 + t59 * t71, -t15 * t97 - t211 * t347 - t228 * t24 - t229 * t3 - t25 * t85 - t271 * t59 + t29 * t76, -t1 * t75 + t13 * t25 - t14 * t26 + t15 * t23 - t2 * t76 + t24 * t261 + t271 * t4 - t3 * t71 - t277, t1 * t24 + t14 * t3 + t2 * t23 + t13 * t4 + t29 * t97 + t85 * t59 - g(1) * (-t148 * t247 - t233 * t250) - g(2) * (t148 * t250 - t233 * t247); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t296, t305 * t253, t298, t296, t297, qJDD(2), t246 * t264 - t335, g(3) * t246 + t249 * t264, 0, 0, t311, t108, t81, -t311, -t273, t236, -t134 * t238 + (-t164 * t303 + t236 * t345 - t238 * t302) * pkin(2) + t255, t135 * t238 + (-t166 * t303 - t236 * t245 - t238 * t289) * pkin(2) + t263, (t128 + t134) * t166 + (-t127 + t135) * t164 + (t345 * t104 - t105 * t245 + (-t164 * t345 + t166 * t245) * qJD(3)) * pkin(2), -t127 * t134 - t128 * t135 + (t345 * t65 - t335 - t245 * t282 + (-t127 * t245 + t128 * t345) * qJD(3) + (qJD(1) * t196 + t277) * t246) * pkin(2), t351, t37, t40, -t351, t39, t236, -t121 * t144 + t156 * t236 + t238 * t321 + t262, -t144 * t270 - t158 * t236 - t238 * t320 + t266, -t121 * t320 - t156 * t55 - t158 * t54 - t270 * t321 + t272, t11 * t158 + t10 * t156 - t140 * t144 - g(3) * t306 + t320 * t48 + t321 * t47 - t277 * (-t342 - t344), -t355, t12, t8, t355, t9, t228, t113 * t228 + t229 * t328 - t71 * t86 + t256, -t114 * t228 - t229 * t329 + t271 * t86 + t258, t113 * t15 + t114 * t261 + t271 * t328 - t329 * t71 + t275, t1 * t114 + t2 * t113 - t85 * t86 - g(3) * t295 - t277 * (t278 - t344) + t329 * t14 + t328 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, t108, t81, -t311, -t273, t236, t128 * t238 + t255, t127 * t238 + t263, 0, 0, t351, t37, t40, -t351, t39, t236, -t238 * t52 + (-t121 * t166 + t236 * t243) * pkin(3) + t262, t238 * t53 + (-t166 * t270 - t236 * t242) * pkin(3) + t266, t121 * t53 + t270 * t52 + (-t242 * t54 - t243 * t55) * pkin(3) + t272, -t47 * t52 - t48 * t53 + (t10 * t243 + t11 * t242 - t140 * t166 + t354) * pkin(3), -t355, t12, t8, t355, t9, t228, t157 * t228 + t229 * t322 - t71 * t91 + t256, -t159 * t228 - t229 * t323 + t271 * t91 + t258, t15 * t157 + t159 * t261 + t271 * t322 - t323 * t71 + t275, t1 * t159 + t2 * t157 - t85 * t91 - g(3) * (t207 + t221) - t277 * t278 + t323 * t14 + t322 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 + t317, t55 - t314, -t315 - t318, t121 * t48 + t270 * t47 + t259, 0, 0, 0, 0, 0, 0, -t261 - t325, -t15 - t324, -t332 - t334, -t13 * t271 + t14 * t71 + t259 + t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355, t12, t8, t355, t9, t228, t14 * t229 + t256, t13 * t229 + t258, 0, 0;];
tau_reg = t5;
