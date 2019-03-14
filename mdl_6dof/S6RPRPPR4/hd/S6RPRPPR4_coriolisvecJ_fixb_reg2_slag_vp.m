% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:23
% EndTime: 2019-03-09 02:48:34
% DurationCPUTime: 4.23s
% Computational Cost: add. (6887->417), mult. (18210->540), div. (0->0), fcn. (13943->8), ass. (0->211)
t188 = sin(pkin(9));
t190 = cos(pkin(9));
t192 = sin(qJ(3));
t285 = cos(qJ(3));
t164 = t188 * t285 + t190 * t192;
t151 = t164 * qJD(1);
t187 = sin(pkin(10));
t189 = cos(pkin(10));
t119 = qJD(3) * t189 - t151 * t187;
t191 = sin(qJ(6));
t212 = qJD(3) * t187 + t151 * t189;
t284 = cos(qJ(6));
t305 = t119 * t284 + t191 * t212;
t314 = t305 ^ 2;
t306 = -t119 * t191 + t212 * t284;
t313 = t306 ^ 2;
t238 = qJD(1) * t285;
t177 = t190 * t238;
t255 = t192 * t188;
t239 = qJD(1) * t255;
t149 = -t177 + t239;
t245 = qJD(6) - t149;
t309 = t245 * t306;
t176 = qJD(3) * t177;
t136 = qJD(3) * t239 - t176;
t163 = t187 * t284 - t189 * t191;
t34 = qJD(6) * t306 + t136 * t163;
t312 = -t34 + t309;
t311 = t245 * t305;
t126 = t189 * t136;
t270 = t119 * t149;
t310 = t126 - t270;
t298 = t187 * t191 + t189 * t284;
t153 = t298 * qJD(6);
t274 = -t149 * t298 + t153;
t236 = qJD(6) * t284;
t246 = qJD(6) * t191;
t154 = t187 * t236 - t189 * t246;
t273 = -t149 * t163 + t154;
t281 = pkin(7) + qJ(2);
t170 = t281 * t188;
t172 = t281 * t190;
t114 = -t170 * t192 + t172 * t285;
t84 = qJD(2) * t164 + qJD(3) * t114;
t308 = t189 * t164 * qJD(5) - t84;
t156 = t164 * qJD(3);
t137 = qJD(1) * t156;
t288 = pkin(4) + pkin(5);
t166 = qJD(1) * t172;
t147 = t192 * t166;
t244 = qJD(1) * qJD(2);
t235 = t192 * t244;
t165 = qJD(1) * t170;
t225 = qJD(2) * t238;
t237 = qJD(3) * t285;
t252 = -t165 * t237 + t190 * t225;
t203 = -t188 * t235 + t252;
t72 = (qJD(4) - t147) * qJD(3) + t203;
t64 = t187 * t72;
t66 = pkin(3) * t137 + qJ(4) * t136 - qJD(4) * t151;
t12 = t64 + (pkin(8) * t136 - t66) * t189 - t288 * t137;
t124 = t187 * t136;
t27 = t187 * t66 + t189 * t72;
t14 = qJ(5) * t137 + qJD(5) * t149 + t27;
t13 = -pkin(8) * t124 + t14;
t181 = -pkin(2) * t190 - pkin(1);
t167 = qJD(1) * t181 + qJD(2);
t80 = t149 * pkin(3) - t151 * qJ(4) + t167;
t103 = -t192 * t165 + t166 * t285;
t96 = qJD(3) * qJ(4) + t103;
t46 = -t187 * t96 + t189 * t80;
t223 = qJD(5) - t46;
t17 = -pkin(8) * t212 - t149 * t288 + t223;
t47 = t187 * t80 + t189 * t96;
t36 = qJ(5) * t149 + t47;
t22 = -pkin(8) * t119 + t36;
t208 = -t17 * t284 + t191 * t22;
t1 = -qJD(6) * t208 + t191 * t12 + t13 * t284;
t307 = t208 * t245 + t1;
t304 = t212 ^ 2;
t26 = t189 * t66 - t64;
t20 = -t137 * pkin(4) - t26;
t303 = -t149 * t36 + t20;
t206 = t190 * t285 - t255;
t104 = t137 * t206;
t61 = t149 * t156 - t104;
t300 = -t170 * t285 - t192 * t172;
t146 = t149 ^ 2;
t229 = t137 * t187 + t146 * t189;
t268 = t212 * t151;
t299 = t229 + t268;
t253 = t137 * t189 - t146 * t187;
t263 = t151 * t119;
t297 = t253 + t263;
t296 = qJD(3) * t151;
t33 = t119 * t236 + t136 * t298 + t212 * t246;
t295 = -t273 * t306 + t298 * t33;
t231 = t149 * t212 - t124;
t249 = qJD(3) * t192;
t75 = -t165 * t249 + t166 * t237 + t188 * t225 + t190 * t235;
t198 = pkin(4) * t124 + qJD(5) * t212 - t75;
t271 = qJ(5) * t189;
t18 = t136 * (pkin(5) * t187 - t271) + t198;
t294 = -t163 * t137 - t245 * t274;
t234 = t187 * qJ(5) + pkin(3);
t168 = -pkin(4) * t189 - t234;
t272 = qJ(4) * t137;
t102 = -t165 * t285 - t147;
t210 = qJD(3) * pkin(3) - qJD(4) + t102;
t197 = qJ(5) * t212 + t210;
t45 = -pkin(4) * t119 - t197;
t293 = t136 * t168 + (qJD(4) - t45) * t149 + t272;
t155 = t188 * t249 - t190 * t237;
t267 = t212 * t187;
t269 = t119 * t189;
t213 = t267 - t269;
t259 = t164 * t187;
t292 = 0.2e1 * t126 * t259 + t155 * t213;
t199 = -t136 * t206 - t137 * t164 + t149 * t155;
t291 = t119 * t156 + t187 * t199;
t207 = t187 * t288 - t271;
t183 = t187 ^ 2;
t185 = t189 ^ 2;
t290 = (t183 - t185) * t136 - t213 * t149;
t289 = t151 ^ 2;
t280 = -pkin(8) + qJ(4);
t169 = t280 * t187;
t171 = t280 * t189;
t113 = t169 * t191 + t171 * t284;
t283 = pkin(8) * t149;
t94 = t187 * t102;
t97 = pkin(3) * t151 + qJ(4) * t149;
t23 = t94 + (-t97 + t283) * t189 - t288 * t151;
t55 = t102 * t189 + t187 * t97;
t42 = qJ(5) * t151 + t55;
t32 = -t187 * t283 + t42;
t287 = -qJD(4) * t163 + qJD(6) * t113 - t191 * t32 + t23 * t284;
t111 = t169 * t284 - t171 * t191;
t286 = -qJD(4) * t298 - qJD(6) * t111 + t191 * t23 + t284 * t32;
t282 = t306 * t305;
t76 = pkin(3) * t156 + qJ(4) * t155 - qJD(4) * t164;
t83 = qJD(2) * t206 + qJD(3) * t300;
t39 = t187 * t76 + t189 * t83;
t279 = t151 * t305;
t276 = t306 * t151;
t275 = t75 * t300;
t98 = -pkin(3) * t206 - qJ(4) * t164 + t181;
t59 = t114 * t189 + t187 * t98;
t266 = t136 * t164;
t265 = t149 * t151;
t262 = t155 * t187;
t250 = t188 ^ 2 + t190 ^ 2;
t248 = qJD(4) * t212;
t247 = qJD(5) * t187;
t50 = -qJ(5) * t206 + t59;
t78 = t187 * t83;
t38 = t189 * t76 - t78;
t54 = t189 * t97 - t94;
t106 = t187 * t114;
t58 = t189 * t98 - t106;
t233 = t250 * qJD(1) ^ 2;
t232 = -t149 * t207 + t103 + t247;
t21 = qJ(5) * t156 - qJD(5) * t206 + t39;
t227 = -t163 * t34 + t274 * t305;
t224 = t119 * t262 - t183 * t266;
t222 = t137 * t298 - t245 * t273;
t220 = -t253 + t263;
t219 = pkin(4) * t187 - t271;
t217 = -t187 * t46 + t189 * t47;
t216 = t136 * t300 + t164 * t75;
t214 = -t119 ^ 2 - t304;
t211 = 0.2e1 * t250 * t244;
t30 = t106 + (-pkin(8) * t164 - t98) * t189 + t288 * t206;
t40 = pkin(8) * t259 + t50;
t9 = -t191 * t40 + t284 * t30;
t6 = t17 * t191 + t22 * t284;
t10 = t191 * t30 + t284 * t40;
t204 = t310 * t187;
t24 = qJ(5) * t126 - t198;
t63 = t164 * t219 - t300;
t202 = t136 * t63 + t155 * t45 - t164 * t24;
t201 = t155 * t210 + t216;
t200 = -t149 * t284 + t236;
t196 = pkin(3) * t136 - t272 + (-qJD(4) - t210) * t149;
t195 = (t267 + t269) * t149 + (t183 + t185) * t136;
t2 = -qJD(6) * t6 + t12 * t284 - t191 * t13;
t157 = t189 * t288 + t234;
t108 = qJD(4) * t269;
t91 = t298 * t164;
t90 = t163 * t164;
t74 = (-qJD(3) * t166 - t188 * t244) * t192 + t252;
t60 = t231 * t189;
t57 = -t149 * t219 + t103;
t53 = -t164 * t207 + t300;
t52 = -t155 * t189 * t212 - t185 * t266;
t51 = pkin(4) * t206 - t58;
t49 = t153 * t164 + t155 * t163;
t48 = -t154 * t164 + t155 * t298;
t44 = -t151 * pkin(4) - t54;
t41 = t229 - t268;
t37 = -t155 * t219 - t308;
t35 = -pkin(4) * t149 + t223;
t31 = t119 * t288 + t197;
t29 = -t156 * pkin(4) - t38;
t25 = t155 * t207 + t308;
t19 = t156 * t212 - t189 * t199;
t16 = -pkin(8) * t262 + t21;
t15 = t78 + (pkin(8) * t155 - t76) * t189 - t288 * t156;
t4 = -qJD(6) * t10 + t15 * t284 - t191 * t16;
t3 = qJD(6) * t9 + t15 * t191 + t16 * t284;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, qJ(2) * t211, -t151 * t155 - t266, -t151 * t156 + t199, -t155 * qJD(3), t61, -t156 * qJD(3), 0, -qJD(3) * t84 + t137 * t181 + t156 * t167, -qJD(3) * t83 - t136 * t181 - t155 * t167, t102 * t155 - t103 * t156 - t114 * t137 - t149 * t83 + t151 * t84 + t206 * t74 + t216, -t102 * t84 + t103 * t83 + t114 * t74 - t275, t52, t292, t19, t224, t291, t61, -t119 * t84 + t58 * t137 + t38 * t149 + t46 * t156 + t187 * t201 - t206 * t26, -t59 * t137 - t39 * t149 - t47 * t156 + t189 * t201 + t206 * t27 + t212 * t84, t39 * t119 - t38 * t212 + (t136 * t58 + t155 * t46 - t164 * t26) * t189 + (t136 * t59 + t155 * t47 - t164 * t27) * t187, -t210 * t84 + t26 * t58 + t27 * t59 + t38 * t46 + t39 * t47 - t275, t52, t19, -t292, t61, -t291, t224, -t119 * t37 - t51 * t137 - t29 * t149 - t35 * t156 - t187 * t202 + t20 * t206, t21 * t119 + t29 * t212 + (-t136 * t51 - t155 * t35 + t164 * t20) * t189 + (t136 * t50 - t14 * t164 + t155 * t36) * t187, t50 * t137 - t14 * t206 + t21 * t149 + t36 * t156 + t189 * t202 - t212 * t37, t14 * t50 + t20 * t51 + t21 * t36 + t24 * t63 + t29 * t35 + t37 * t45, -t306 * t48 - t33 * t91, t305 * t48 - t306 * t49 - t33 * t90 - t34 * t91, -t137 * t91 - t156 * t306 - t206 * t33 - t245 * t48, t305 * t49 - t34 * t90, -t137 * t90 + t156 * t305 - t206 * t34 - t245 * t49, -t156 * t245 - t104, -t137 * t9 + t156 * t208 - t18 * t90 + t2 * t206 + t245 * t4 + t25 * t305 + t31 * t49 + t34 * t53, -t1 * t206 + t10 * t137 + t156 * t6 + t18 * t91 - t245 * t3 + t25 * t306 - t31 * t48 - t33 * t53, t1 * t90 - t10 * t34 - t2 * t91 - t208 * t48 - t3 * t305 - t306 * t4 + t33 * t9 - t49 * t6, t1 * t10 + t18 * t53 + t2 * t9 - t208 * t4 + t25 * t31 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233, -qJ(2) * t233, 0, 0, 0, 0, 0, 0, 0.2e1 * t296, t176 + (-t149 - t239) * qJD(3), -t146 - t289, t102 * t151 + t103 * t149, 0, 0, 0, 0, 0, 0, t297, -t299, t195, t149 * t217 + t151 * t210 + t27 * t187 + t26 * t189, 0, 0, 0, 0, 0, 0, t297, t195, t299, t14 * t187 - t45 * t151 - t20 * t189 + (t187 * t35 + t189 * t36) * t149, 0, 0, 0, 0, 0, 0, t222 + t279, t276 - t294, t227 - t295, t1 * t163 + t31 * t151 - t2 * t298 + t208 * t273 - t274 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, -t146 + t289, t176 + (t149 - t239) * qJD(3), -t265, 0, 0, qJD(3) * t103 - t151 * t167 - t75, t167 * t149 + (t102 + t147) * qJD(3) - t203, 0, 0, t60, t290, t41, t204, -t220, -t265, t103 * t119 - t54 * t149 - t46 * t151 + t187 * t196 - t75 * t189, -t103 * t212 + t55 * t149 + t47 * t151 + t75 * t187 + t189 * t196, -t55 * t119 + t54 * t212 + t108 + (-t149 * t46 + t27) * t189 + (-t149 * t47 + t248 - t26) * t187, -t75 * pkin(3) + t210 * t103 - t46 * t54 - t47 * t55 + t217 * qJD(4) + (-t26 * t187 + t27 * t189) * qJ(4), t60, t41, -t290, -t265, t220, t204, t57 * t119 + t44 * t149 + t35 * t151 - t24 * t189 + (qJD(5) * t119 - t293) * t187, -t42 * t119 - t44 * t212 + t108 + (t149 * t35 + t14) * t189 + (t248 + t303) * t187, -t42 * t149 - t36 * t151 - t24 * t187 + (t57 + t247) * t212 + t293 * t189, t24 * t168 - t35 * t44 - t36 * t42 - t45 * t57 + (qJ(4) * t14 + qJD(4) * t36) * t189 + (qJ(4) * t20 + qJD(4) * t35 - qJD(5) * t45) * t187, -t33 * t163 - t274 * t306, t227 + t295, t276 + t294, t273 * t305 + t298 * t34, t222 - t279, t245 * t151, -t111 * t137 - t151 * t208 + t157 * t34 + t18 * t298 + t232 * t305 - t245 * t287 + t273 * t31, t113 * t137 - t6 * t151 - t157 * t33 + t18 * t163 + t232 * t306 + t245 * t286 - t274 * t31, -t1 * t298 + t111 * t33 - t113 * t34 - t2 * t163 - t208 * t274 - t273 * t6 + t286 * t305 + t287 * t306, t1 * t113 + t2 * t111 + t18 * t157 + t208 * t287 + t232 * t31 - t286 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, -t310, t214, -t119 * t47 + t212 * t46 + t75, 0, 0, 0, 0, 0, 0, t231, t214, t310, -t119 * t36 - t212 * t35 + t24, 0, 0, 0, 0, 0, 0, -t34 - t309, t33 + t311, t313 + t314, t208 * t306 - t305 * t6 - t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 * t212 - t296, -t126 - t270, -t146 - t304, t212 * t45 + t303, 0, 0, 0, 0, 0, 0, -t191 * t245 ^ 2 - t137 * t284 - t212 * t305, t191 * t137 - t200 * t245 - t212 * t306, t191 * t312 - t200 * t305 + t284 * t33, t191 * t307 + t2 * t284 + t200 * t6 - t212 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, t313 - t314, -t33 + t311, -t282, t312, -t137, t245 * t6 - t306 * t31 + t2, t305 * t31 - t307, 0, 0;];
tauc_reg  = t5;