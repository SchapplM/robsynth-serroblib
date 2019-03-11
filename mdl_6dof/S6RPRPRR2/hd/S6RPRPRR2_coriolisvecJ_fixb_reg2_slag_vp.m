% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:58
% EndTime: 2019-03-09 03:39:10
% DurationCPUTime: 4.65s
% Computational Cost: add. (9944->440), mult. (23376->587), div. (0->0), fcn. (16745->10), ass. (0->234)
t194 = sin(qJ(5));
t196 = cos(qJ(5));
t190 = sin(pkin(11));
t195 = sin(qJ(3));
t197 = cos(qJ(3));
t287 = cos(pkin(11));
t165 = t190 * t197 + t287 * t195;
t264 = qJD(1) * t165;
t127 = -t196 * qJD(3) + t194 * t264;
t129 = qJD(3) * t194 + t196 * t264;
t193 = sin(qJ(6));
t312 = cos(qJ(6));
t212 = -t193 * t127 + t312 * t129;
t65 = t312 * t127 + t129 * t193;
t310 = t65 * t212;
t231 = t287 * t197;
t173 = qJD(1) * t231;
t262 = qJD(1) * t195;
t154 = t190 * t262 - t173;
t270 = t193 * t194;
t211 = t312 * t196 - t270;
t320 = qJD(5) + qJD(6);
t239 = t312 * qJD(6);
t321 = t312 * qJD(5) + t239;
t290 = -t211 * t154 - t321 * t196 + t320 * t270;
t243 = t312 * t194;
t167 = t193 * t196 + t243;
t117 = t320 * t167;
t289 = t167 * t154 + t117;
t178 = pkin(3) * t190 + pkin(8);
t309 = pkin(9) + t178;
t234 = qJD(5) * t309;
t278 = t154 * t194;
t179 = sin(pkin(10)) * pkin(1) + pkin(7);
t170 = t179 * qJD(1);
t227 = qJ(4) * qJD(1) + t170;
t253 = t195 * qJD(2);
t131 = t197 * t227 + t253;
t122 = t190 * t131;
t187 = t197 * qJD(2);
t130 = -t195 * t227 + t187;
t70 = t287 * t130 - t122;
t250 = pkin(3) * t262;
t97 = pkin(4) * t264 + pkin(8) * t154 + t250;
t42 = t194 * t97 + t196 * t70;
t333 = pkin(9) * t278 + t194 * t234 + t42;
t41 = -t194 * t70 + t196 * t97;
t332 = pkin(5) * t264 + t41 + (pkin(9) * t154 + t234) * t196;
t252 = qJD(1) * qJD(3);
t238 = t195 * t252;
t172 = t190 * t238;
t144 = qJD(3) * t173 - t172;
t331 = qJD(3) * qJD(5) + t144;
t256 = qJD(5) * t194;
t330 = t256 + t278;
t329 = t212 ^ 2 - t65 ^ 2;
t149 = qJD(5) + t154;
t142 = qJD(6) + t149;
t255 = qJD(5) * t196;
t244 = t331 * t194 + t264 * t255;
t254 = qJD(6) * t193;
t84 = -t331 * t196 + t264 * t256;
t27 = t127 * t239 + t129 * t254 + t193 * t244 + t312 * t84;
t328 = t142 * t65 - t27;
t125 = qJD(3) * pkin(3) + t130;
t232 = t287 * t131;
t62 = t190 * t125 + t232;
t59 = qJD(3) * pkin(8) + t62;
t181 = -cos(pkin(10)) * pkin(1) - pkin(2);
t169 = -pkin(3) * t197 + t181;
t263 = qJD(1) * t169;
t153 = qJD(4) + t263;
t81 = pkin(4) * t154 - pkin(8) * t264 + t153;
t38 = -t194 * t59 + t196 * t81;
t30 = -pkin(9) * t129 + t38;
t26 = pkin(5) * t149 + t30;
t39 = t194 * t81 + t196 * t59;
t31 = -pkin(9) * t127 + t39;
t156 = t165 * qJD(3);
t143 = qJD(1) * t156;
t184 = qJD(3) * t187;
t259 = qJD(3) * t195;
t137 = -t170 * t259 + t184;
t257 = qJD(4) * t197;
t110 = (-qJ(4) * t259 + t257) * qJD(1) + t137;
t233 = t287 * t110;
t258 = qJD(4) * t195;
t322 = -qJD(1) * t258 - t131 * qJD(3);
t53 = t322 * t190 + t233;
t176 = pkin(3) * t238;
t86 = pkin(4) * t143 - pkin(8) * t144 + t176;
t15 = -qJD(5) * t39 - t194 * t53 + t196 * t86;
t6 = pkin(5) * t143 + pkin(9) * t84 + t15;
t14 = t194 * t86 + t196 * t53 + t81 * t255 - t59 * t256;
t9 = -t244 * pkin(9) + t14;
t204 = -t193 * t6 - t26 * t239 + t31 * t254 - t312 * t9;
t61 = t287 * t125 - t122;
t58 = -qJD(3) * pkin(4) - t61;
t45 = t127 * pkin(5) + t58;
t327 = t45 * t65 + t204;
t325 = -t149 * t38 + t14;
t324 = t149 * t39 + t15;
t229 = t149 * t194;
t323 = t129 * t229;
t248 = t312 * t31;
t8 = t193 * t26 + t248;
t2 = -qJD(6) * t8 - t193 * t9 + t312 * t6;
t319 = -t45 * t212 + t2;
t28 = qJD(6) * t212 - t193 * t84 + t312 * t244;
t318 = t142 * t212 - t28;
t317 = -t142 * t290 + t143 * t167;
t134 = t196 * t143;
t241 = t165 * t256;
t208 = -t190 * t195 + t231;
t159 = t208 * qJD(3);
t274 = t159 * t196;
t209 = t241 - t274;
t316 = t165 * t134 - t149 * t209;
t315 = t156 * t127 - t208 * t244;
t314 = -t211 * t27 - t212 * t289;
t313 = t264 ^ 2;
t311 = pkin(3) * t195;
t162 = t309 * t194;
t163 = t309 * t196;
t109 = -t193 * t162 + t312 * t163;
t308 = t109 * qJD(6) - t333 * t193 + t332 * t312;
t108 = -t312 * t162 - t193 * t163;
t307 = -t108 * qJD(6) + t332 * t193 + t333 * t312;
t103 = t211 * t165;
t35 = t117 * t165 - t211 * t159;
t306 = -t103 * t28 + t35 * t65;
t102 = t167 * t165;
t273 = t165 * t194;
t36 = t159 * t243 - t193 * t241 - t254 * t273 + (t159 * t193 + t321 * t165) * t196;
t305 = -t102 * t143 - t36 * t142;
t304 = t156 * t212 + t208 * t27;
t221 = t244 * t196;
t303 = -t127 * t274 - t165 * t221;
t302 = t129 * t156 + t208 * t84;
t267 = qJ(4) + t179;
t160 = t267 * t195;
t161 = t267 * t197;
t105 = t287 * t160 + t161 * t190;
t52 = t190 * t110 - t287 * t322;
t301 = t105 * t52;
t298 = t264 * t65;
t297 = t264 * t212;
t296 = t208 * t52;
t294 = t193 * t31;
t293 = t194 * t52;
t292 = t52 * t196;
t291 = t84 * t194;
t106 = -t190 * t160 + t287 * t161;
t100 = t196 * t106;
t101 = -pkin(4) * t208 - pkin(8) * t165 + t169;
t48 = t194 * t101 + t100;
t76 = t194 * t244;
t288 = -t127 * t255 - t76;
t286 = t127 * t154;
t285 = t127 * t264;
t284 = t127 * t194;
t283 = t129 * t127;
t282 = t129 * t264;
t281 = t129 * t196;
t111 = t143 * t208;
t279 = t143 * t194;
t276 = t264 * t154;
t275 = t159 * t194;
t272 = t165 * t196;
t271 = t170 * t195;
t198 = qJD(3) ^ 2;
t269 = t198 * t195;
t268 = t198 * t197;
t266 = -t165 * t143 - t159 * t154;
t265 = t195 ^ 2 - t197 ^ 2;
t171 = qJD(1) * t181;
t261 = qJD(1) * t197;
t260 = qJD(3) * t159;
t249 = pkin(3) * t259;
t247 = t129 * t275;
t199 = qJD(1) ^ 2;
t245 = t195 * t199 * t197;
t230 = qJD(3) * t267;
t132 = -t195 * t230 + t257;
t206 = -t197 * t230 - t258;
t80 = t287 * t132 + t190 * t206;
t98 = pkin(4) * t156 - pkin(8) * t159 + t249;
t235 = -t194 * t80 + t196 * t98;
t47 = t196 * t101 - t106 * t194;
t69 = t130 * t190 + t232;
t79 = t132 * t190 - t287 * t206;
t228 = t149 * t196;
t226 = 0.2e1 * t264;
t225 = -t167 * t28 + t290 * t65;
t224 = t197 * t238;
t223 = t330 * pkin(5) - t69;
t222 = -t289 * t142 + t211 * t143;
t180 = -t287 * pkin(3) - pkin(4);
t220 = -t102 * t27 + t212 * t36;
t219 = -t156 * t65 + t208 * t28;
t218 = -t194 * t39 - t196 * t38;
t217 = t194 * t38 - t196 * t39;
t216 = -t103 * t143 + t142 * t35;
t215 = t144 * t208 - t156 * t264;
t146 = t170 * t197 + t253;
t214 = 0.2e1 * qJD(3) * t171;
t213 = -t330 * t149 + t134;
t40 = -pkin(5) * t208 - pkin(9) * t272 + t47;
t43 = -pkin(9) * t273 + t48;
t18 = -t193 * t43 + t312 * t40;
t19 = t193 * t40 + t312 * t43;
t210 = t165 * t255 + t275;
t21 = t101 * t255 - t106 * t256 + t194 * t98 + t196 * t80;
t207 = -t178 * t143 + t149 * t58;
t202 = -t143 * t273 - t149 * t210;
t201 = qJD(5) * t218 + t14 * t196 - t15 * t194;
t138 = t146 * qJD(3);
t145 = t187 - t271;
t200 = t137 * t197 + t138 * t195 + (-t145 * t197 - t146 * t195) * qJD(3);
t168 = -t196 * pkin(5) + t180;
t152 = t154 ^ 2;
t148 = qJD(3) * t156;
t74 = pkin(5) * t273 + t105;
t44 = pkin(5) * t210 + t79;
t34 = t244 * pkin(5) + t52;
t22 = -t48 * qJD(5) + t235;
t17 = -pkin(9) * t210 + t21;
t16 = -pkin(9) * t274 + pkin(5) * t156 + (-t100 + (pkin(9) * t165 - t101) * t194) * qJD(5) + t235;
t11 = t312 * t30 - t294;
t10 = -t193 * t30 - t248;
t7 = t312 * t26 - t294;
t4 = -t19 * qJD(6) + t312 * t16 - t193 * t17;
t3 = t18 * qJD(6) + t193 * t16 + t312 * t17;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t224, -0.2e1 * t265 * t252, t268, -0.2e1 * t224, -t269, 0, -t179 * t268 + t195 * t214, t179 * t269 + t197 * t214, t200, t200 * t179, t144 * t165 + t159 * t264, t215 + t266, t260, t154 * t156 - t111, -t148, 0, t143 * t169 + t153 * t156 + (-t79 + (-qJD(1) * t208 + t154) * t311) * qJD(3), t144 * t169 + t153 * t159 + (t226 * t311 - t80) * qJD(3), t105 * t144 - t106 * t143 - t154 * t80 - t156 * t62 - t159 * t61 + t165 * t52 + t208 * t53 + t264 * t79, t301 + t106 * t53 - t61 * t79 + t62 * t80 + (t153 + t263) * t249, -t129 * t209 - t272 * t84, -t247 + (t291 + (-t281 + t284) * qJD(5)) * t165 + t303, t302 + t316, t127 * t210 + t165 * t76, t202 - t315, t149 * t156 - t111, t22 * t149 + t47 * t143 - t15 * t208 + t38 * t156 + t79 * t127 + t105 * t244 + t58 * t275 + (t255 * t58 + t293) * t165, t58 * t274 - t105 * t84 + t129 * t79 + t14 * t208 - t143 * t48 - t149 * t21 - t156 * t39 + (-t256 * t58 + t292) * t165, -t21 * t127 - t48 * t244 - t22 * t129 + t47 * t84 + t218 * t159 + (qJD(5) * t217 - t14 * t194 - t15 * t196) * t165, t14 * t48 + t15 * t47 + t21 * t39 + t22 * t38 + t58 * t79 + t301, -t103 * t27 - t212 * t35, -t220 + t306, -t216 + t304, t102 * t28 + t36 * t65, t219 + t305, t142 * t156 - t111, t102 * t34 + t142 * t4 + t143 * t18 + t156 * t7 - t2 * t208 + t28 * t74 + t36 * t45 + t44 * t65, t103 * t34 - t142 * t3 - t143 * t19 - t156 * t8 - t204 * t208 + t212 * t44 - t27 * t74 - t35 * t45, t102 * t204 - t103 * t2 + t18 * t27 - t19 * t28 - t212 * t4 - t3 * t65 + t35 * t7 - t36 * t8, t18 * t2 - t19 * t204 + t3 * t8 + t34 * t74 + t4 * t7 + t44 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269, -t268, 0, t137 * t195 - t138 * t197 + (-t145 * t195 + t146 * t197) * qJD(3), 0, 0, 0, 0, 0, 0, -t148, -t260, -t215 + t266, -t156 * t61 + t159 * t62 + t165 * t53 - t296, 0, 0, 0, 0, 0, 0, t202 + t315, t302 - t316, t247 + (-t291 + (t281 + t284) * qJD(5)) * t165 + t303, t156 * t58 - t159 * t217 + t165 * t201 - t296, 0, 0, 0, 0, 0, 0, -t219 + t305, t216 + t304, t220 + t306, -t102 * t2 - t103 * t204 + t156 * t45 - t208 * t34 - t35 * t8 - t36 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t265 * t199, 0, t245, 0, 0, -t171 * t262, -t171 * t261 - t184 + (t145 + t271) * qJD(3), 0, 0, t276, -t152 + t313, -t172 + (t173 + t154) * qJD(3), -t276, 0, 0, qJD(3) * t69 - t153 * t264 - t154 * t250 - t52, -t233 + t153 * t154 + (-pkin(3) * t264 + qJD(4) * t190) * t262 + (-t190 * (-qJ(4) * t261 - t146) + t70) * qJD(3) (t62 - t69) * t264 + (-t61 + t70) * t154 + (-t143 * t190 - t144 * t287) * pkin(3), t61 * t69 - t62 * t70 + (-t153 * t262 + t190 * t53 - t287 * t52) * pkin(3), t129 * t228 - t291 (-t84 - t286) * t196 - t323 + t288, t149 * t228 + t279 - t282, t127 * t229 - t221, t213 + t285, -t149 * t264, t180 * t244 - t292 - t38 * t264 - t69 * t127 + (-t178 * t255 - t41) * t149 + t207 * t194, -t129 * t69 + t264 * t39 - t180 * t84 + t293 + (t178 * t256 + t42) * t149 + t207 * t196, t42 * t127 + t41 * t129 + (-t178 * t244 + t14 - t38 * t154 + (t129 * t178 - t38) * qJD(5)) * t196 + (-t39 * t154 - t178 * t84 - t15 + (t127 * t178 - t39) * qJD(5)) * t194, t178 * t201 + t180 * t52 - t38 * t41 - t39 * t42 - t58 * t69, -t167 * t27 - t212 * t290, t225 + t314, -t297 + t317, -t211 * t28 + t289 * t65, t222 + t298, -t142 * t264, t108 * t143 - t308 * t142 + t168 * t28 - t211 * t34 + t223 * t65 - t264 * t7 + t289 * t45, -t109 * t143 + t307 * t142 + t167 * t34 - t168 * t27 + t212 * t223 + t264 * t8 - t290 * t45, t108 * t27 - t109 * t28 - t167 * t2 - t204 * t211 + t212 * t308 - t289 * t8 + t290 * t7 + t307 * t65, t108 * t2 - t109 * t204 + t168 * t34 + t223 * t45 - t307 * t8 - t308 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226 * qJD(3), -t172 + (t173 - t154) * qJD(3), -t152 - t313, t154 * t62 + t264 * t61 + t176, 0, 0, 0, 0, 0, 0, t213 - t285, -t149 ^ 2 * t196 - t279 - t282 (t84 - t286) * t196 + t323 + t288, t325 * t194 + t324 * t196 - t264 * t58, 0, 0, 0, 0, 0, 0, t222 - t298, -t297 - t317, t225 - t314, -t167 * t204 + t2 * t211 - t264 * t45 - t289 * t7 - t290 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, -t127 ^ 2 + t129 ^ 2, t127 * t149 - t84, -t283, t129 * t149 - t244, t143, -t129 * t58 + t324, t127 * t58 - t325, 0, 0, t310, t329, t328, -t310, t318, t143, -t10 * t142 + (-t129 * t65 - t142 * t254 + t312 * t143) * pkin(5) + t319, t11 * t142 + (-t129 * t212 - t142 * t239 - t143 * t193) * pkin(5) + t327, t10 * t212 + t11 * t65 + t8 * t212 - t7 * t65 + (t312 * t27 - t193 * t28 + (t193 * t212 - t312 * t65) * qJD(6)) * pkin(5), -t7 * t10 - t8 * t11 + (t312 * t2 - t204 * t193 - t129 * t45 + (-t193 * t7 + t312 * t8) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t329, t328, -t310, t318, t143, t8 * t142 + t319, t7 * t142 + t327, 0, 0;];
tauc_reg  = t1;
