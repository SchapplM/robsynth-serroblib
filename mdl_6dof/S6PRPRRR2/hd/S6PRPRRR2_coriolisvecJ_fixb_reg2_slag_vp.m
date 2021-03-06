% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:55
% EndTime: 2019-03-08 20:30:09
% DurationCPUTime: 5.43s
% Computational Cost: add. (6913->428), mult. (17082->613), div. (0->0), fcn. (13157->12), ass. (0->218)
t165 = sin(pkin(6));
t164 = sin(pkin(12));
t166 = cos(pkin(12));
t171 = sin(qJ(2));
t175 = cos(qJ(2));
t192 = t164 * t175 + t166 * t171;
t110 = t192 * t165;
t102 = qJD(1) * t110;
t170 = sin(qJ(4));
t174 = cos(qJ(4));
t199 = pkin(4) * t170 - pkin(9) * t174;
t137 = t199 * qJD(4);
t303 = t137 - t102;
t240 = qJD(1) * t165;
t221 = t171 * t240;
t139 = t164 * t221;
t217 = t175 * t240;
t105 = t166 * t217 - t139;
t189 = -t174 * pkin(4) - t170 * pkin(9) - pkin(3);
t281 = t166 * pkin(2);
t126 = t189 - t281;
t154 = t164 * pkin(2) + pkin(8);
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t229 = t173 * qJD(4);
t234 = qJD(5) * t173;
t236 = qJD(5) * t169;
t246 = t173 * t174;
t275 = t126 * t234 + (-t170 * t229 - t174 * t236) * t154 - t105 * t246 + t303 * t169;
t230 = t169 * qJD(4);
t251 = t169 * t174;
t302 = -t170 * t154 * t230 - t105 * t251 - t173 * t303;
t133 = t154 * t246;
t188 = pkin(5) * t170 - pkin(10) * t246;
t301 = -t188 * qJD(4) - (-t133 + (pkin(10) * t170 - t126) * t169) * qJD(5) + t302;
t216 = t174 * t230;
t218 = t170 * t234;
t182 = t216 + t218;
t300 = -pkin(10) * t182 + t275;
t239 = qJD(2) * t170;
t128 = t169 * t239 - t229;
t130 = t173 * t239 + t230;
t168 = sin(qJ(6));
t172 = cos(qJ(6));
t193 = t168 * t128 - t172 * t130;
t76 = t172 * t128 + t168 * t130;
t279 = t76 * t193;
t285 = -pkin(10) - pkin(9);
t222 = qJD(5) * t285;
t228 = t174 * qJD(2);
t134 = t199 * qJD(2);
t167 = cos(pkin(6));
t151 = t167 * qJD(1) + qJD(3);
t138 = qJD(2) * pkin(2) + t217;
t95 = t164 * t138 + t166 * t221;
t89 = qJD(2) * pkin(8) + t95;
t289 = t174 * t151 - t170 * t89;
t46 = t169 * t134 + t173 * t289;
t299 = -t46 + (pkin(10) * t228 + t222) * t169;
t45 = t173 * t134 - t169 * t289;
t298 = qJD(2) * t188 - t173 * t222 + t45;
t227 = qJD(2) * qJD(4);
t297 = qJD(4) * qJD(5) + t174 * t227;
t296 = t193 ^ 2 - t76 ^ 2;
t152 = -qJD(5) + t228;
t248 = t170 * t151;
t66 = t174 * t89 + t248;
t63 = qJD(4) * pkin(9) + t66;
t94 = t166 * t138 - t139;
t71 = qJD(2) * t189 - t94;
t30 = -t169 * t63 + t173 * t71;
t27 = -t130 * pkin(10) + t30;
t22 = -t152 * pkin(5) + t27;
t31 = t169 * t71 + t173 * t63;
t28 = -t128 * pkin(10) + t31;
t272 = t172 * t28;
t11 = t168 * t22 + t272;
t235 = qJD(5) * t170;
t212 = qJD(2) * t235;
t100 = t169 * t212 - t297 * t173;
t109 = (t164 * t171 - t166 * t175) * t165;
t104 = qJD(2) * t109;
t99 = qJD(1) * t104;
t41 = qJD(4) * t289 - t174 * t99;
t72 = (t137 + t102) * qJD(2);
t15 = -qJD(5) * t31 - t169 * t41 + t173 * t72;
t156 = t170 * t227;
t8 = pkin(5) * t156 + t100 * pkin(10) + t15;
t14 = t169 * t72 + t173 * t41 + t71 * t234 - t63 * t236;
t223 = t297 * t169 + t173 * t212;
t9 = -t223 * pkin(10) + t14;
t2 = -qJD(6) * t11 - t168 * t9 + t172 * t8;
t62 = -qJD(4) * pkin(4) - t289;
t49 = t128 * pkin(5) + t62;
t295 = t49 * t193 + t2;
t148 = -qJD(6) + t152;
t232 = qJD(6) * t172;
t233 = qJD(6) * t168;
t33 = t172 * t100 + t128 * t232 + t130 * t233 + t168 * t223;
t294 = -t76 * t148 - t33;
t1 = (qJD(6) * t22 + t9) * t172 + t168 * t8 - t28 * t233;
t293 = t49 * t76 - t1;
t180 = qJD(6) * t193 + t168 * t100 - t172 * t223;
t292 = t148 * t193 + t180;
t291 = t30 * t152 + t14;
t290 = t31 * t152 - t15;
t288 = t174 * t223;
t82 = t169 * t126 + t133;
t215 = t174 * t229;
t219 = t169 * t235;
t287 = t215 - t219;
t225 = qJD(5) + qJD(6);
t162 = t170 ^ 2;
t190 = qJD(2) * t162 - t152 * t174;
t286 = -t152 * t219 - t190 * t229;
t115 = t173 * t126;
t247 = t170 * t173;
t70 = -pkin(10) * t247 + t115 + (-t154 * t169 - pkin(5)) * t174;
t252 = t169 * t170;
t74 = -pkin(10) * t252 + t82;
t36 = t168 * t70 + t172 * t74;
t284 = qJD(6) * t36 + t300 * t168 + t301 * t172;
t35 = -t168 * t74 + t172 * t70;
t283 = -qJD(6) * t35 + t301 * t168 - t300 * t172;
t282 = pkin(5) * t169;
t42 = t66 * qJD(4) - t170 * t99;
t86 = t110 * t170 - t167 * t174;
t280 = t42 * t86;
t143 = t285 * t169;
t144 = t285 * t173;
t97 = t168 * t143 - t172 * t144;
t278 = qJD(6) * t97 + t299 * t168 + t298 * t172;
t96 = t172 * t143 + t168 * t144;
t277 = -qJD(6) * t96 + t298 * t168 - t299 * t172;
t253 = t168 * t169;
t131 = -t172 * t173 + t253;
t113 = t131 * t170;
t132 = t168 * t173 + t172 * t169;
t85 = t225 * t132;
t50 = t168 * t216 + t170 * t85 - t172 * t215;
t276 = -t113 * t180 + t50 * t76;
t274 = -t82 * qJD(5) - t302;
t273 = t168 * t28;
t271 = t173 * t62;
t268 = t42 * t169;
t267 = t42 * t170;
t266 = t42 * t173;
t265 = t62 * t169;
t103 = qJD(2) * t110;
t98 = qJD(1) * t103;
t264 = t98 * t109;
t263 = t132 * t228 - t85;
t262 = -t131 * t228 - t172 * t234 - t173 * t232 + t225 * t253;
t112 = t132 * t170;
t51 = -t233 * t252 + (t225 * t247 + t216) * t172 + t287 * t168;
t261 = -t112 * t156 + t51 * t148;
t200 = t223 * t173;
t260 = -t128 * t215 - t170 * t200;
t259 = t128 * t152;
t258 = t130 * t128;
t257 = t130 * t152;
t256 = t152 * t169;
t255 = t152 * t173;
t177 = qJD(2) ^ 2;
t254 = t165 * t177;
t250 = t170 * t105;
t249 = t170 * t128;
t176 = qJD(4) ^ 2;
t244 = t176 * t170;
t243 = t176 * t174;
t163 = t174 ^ 2;
t241 = t162 - t163;
t238 = qJD(4) * t170;
t237 = qJD(4) * t174;
t231 = t128 * qJD(5);
t224 = t170 * t177 * t174;
t220 = t130 * t237;
t210 = t33 * t174 - t193 * t238;
t88 = -qJD(2) * pkin(3) - t94;
t209 = -qJD(2) * t88 + t99;
t207 = t100 * t174 + t130 * t238;
t206 = -t100 + t231;
t204 = t152 * t218;
t203 = t130 * t218;
t202 = t174 * t156;
t201 = pkin(5) * t236 - t248 - (qJD(2) * t282 + t89) * t174;
t198 = -t112 * t33 - t193 * t51;
t87 = t110 * t174 + t167 * t170;
t56 = t109 * t173 - t87 * t169;
t57 = t109 * t169 + t87 * t173;
t20 = -t168 * t57 + t172 * t56;
t21 = t168 * t56 + t172 * t57;
t197 = -t169 * t31 - t173 * t30;
t196 = t169 * t30 - t173 * t31;
t195 = t170 * t289 - t174 * t66;
t187 = -t174 * t180 - t76 * t238;
t186 = t102 * qJD(2) - t154 * t176 - t98;
t155 = -pkin(3) - t281;
t185 = qJD(4) * (qJD(2) * t155 + t105 + t88);
t184 = t190 * t169;
t183 = t113 * t156 - t50 * t148;
t179 = qJD(5) * t197 + t14 * t173 - t15 * t169;
t178 = t267 + t41 * t174 + (-t170 * t66 - t174 * t289) * qJD(4);
t160 = -t173 * pkin(5) - pkin(4);
t118 = (t154 + t282) * t170;
t90 = pkin(5) * t182 + t154 * t237;
t81 = -t154 * t251 + t115;
t55 = -qJD(4) * t86 - t104 * t174;
t54 = qJD(4) * t87 - t104 * t170;
t29 = t223 * pkin(5) + t42;
t19 = qJD(5) * t56 + t103 * t169 + t55 * t173;
t18 = -qJD(5) * t57 + t103 * t173 - t55 * t169;
t13 = t172 * t27 - t273;
t12 = -t168 * t27 - t272;
t10 = t172 * t22 - t273;
t4 = -qJD(6) * t21 - t168 * t19 + t172 * t18;
t3 = qJD(6) * t20 + t168 * t18 + t172 * t19;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171 * t254, -t175 * t254, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * qJD(2), t104 * qJD(2), 0, -t94 * t103 - t95 * t104 - t99 * t110 + t264, 0, 0, 0, 0, 0, 0, -t54 * qJD(4) + (-t103 * t174 + t109 * t238) * qJD(2), -t55 * qJD(4) + (t103 * t170 + t109 * t237) * qJD(2) (t170 * t54 + t174 * t55 + (-t170 * t87 + t174 * t86) * qJD(4)) * qJD(2), t88 * t103 - t289 * t54 + t41 * t87 + t66 * t55 + t264 + t280, 0, 0, 0, 0, 0, 0, t54 * t128 - t18 * t152 + t156 * t56 + t223 * t86, -t86 * t100 + t54 * t130 + t19 * t152 - t156 * t57, t56 * t100 - t19 * t128 - t18 * t130 - t223 * t57, t14 * t57 + t15 * t56 + t30 * t18 + t31 * t19 + t62 * t54 + t280, 0, 0, 0, 0, 0, 0, -t4 * t148 + t156 * t20 - t180 * t86 + t54 * t76, t3 * t148 - t156 * t21 - t193 * t54 - t86 * t33, t180 * t21 + t193 * t4 + t20 * t33 - t3 * t76, t1 * t21 + t10 * t4 + t11 * t3 + t2 * t20 + t29 * t86 + t49 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t192 * t240 + t102) * qJD(2) (qJD(1) * t109 + t105) * qJD(2), 0, t94 * t102 - t95 * t105 + (-t164 * t99 - t166 * t98) * pkin(2), 0.2e1 * t202, -0.2e1 * t241 * t227, t243, -0.2e1 * t202, -t244, 0, t170 * t185 + t174 * t186, -t170 * t186 + t174 * t185 (-t162 - t163) * t105 * qJD(2) + t178, -t88 * t102 + t105 * t195 + t154 * t178 + t98 * t155, -t100 * t247 + t287 * t130, -t203 + (-t220 + (t100 + t231) * t170) * t169 + t260, t207 - t286, t182 * t128 + t223 * t252, t204 + t288 + (-t184 - t249) * qJD(4) (-t152 - t228) * t238, -t274 * t152 + (-t15 + (t128 * t154 + t265) * qJD(4)) * t174 + (t154 * t223 + t268 + t62 * t234 - t105 * t128 + (t81 * qJD(2) + t30) * qJD(4)) * t170, t275 * t152 + (t14 + (t130 * t154 + t271) * qJD(4)) * t174 + (-t62 * t236 - t100 * t154 - t105 * t130 + t266 + (-qJD(2) * t82 - t31) * qJD(4)) * t170, -t82 * t223 + t81 * t100 - t274 * t130 - t275 * t128 + t197 * t237 + (qJD(5) * t196 - t14 * t169 - t15 * t173) * t170, -t62 * t250 + t14 * t82 + t15 * t81 + t275 * t31 + t274 * t30 + (t237 * t62 + t267) * t154, t33 * t113 + t193 * t50, -t198 + t276, -t183 + t210, -t112 * t180 + t76 * t51, t187 + t261 (-t148 - t228) * t238, t29 * t112 - t118 * t180 - t2 * t174 + t49 * t51 + t90 * t76 + t284 * t148 + (-t105 * t76 + (qJD(2) * t35 + t10) * qJD(4)) * t170, t1 * t174 - t29 * t113 - t118 * t33 - t49 * t50 - t90 * t193 - t283 * t148 + (t105 * t193 + (-qJD(2) * t36 - t11) * qJD(4)) * t170, -t1 * t112 + t10 * t50 - t11 * t51 + t2 * t113 + t180 * t36 - t193 * t284 + t283 * t76 + t35 * t33, t1 * t36 + t29 * t118 + t2 * t35 + (t90 - t250) * t49 - t283 * t11 - t284 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, -t243, 0, -qJD(4) * t195 + t41 * t170 - t42 * t174, 0, 0, 0, 0, 0, 0, t204 - t288 + (-t184 + t249) * qJD(4), t207 + t286, t203 + (t170 * t206 + t220) * t169 + t260 (-qJD(4) * t196 - t42) * t174 + (qJD(4) * t62 + t179) * t170, 0, 0, 0, 0, 0, 0, -t187 + t261, t183 + t210, t198 + t276, -t1 * t113 - t10 * t51 - t11 * t50 - t2 * t112 - t29 * t174 + t238 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, t241 * t177, 0, t224, 0, 0, t209 * t170, t209 * t174, 0, 0, -t100 * t169 - t130 * t255 (-t100 + t259) * t173 + (-t223 + t257) * t169, -t152 * t234 + (t152 * t246 + (-t130 + t230) * t170) * qJD(2), -t128 * t256 - t200, t152 * t236 + (-t152 * t251 + (t128 + t229) * t170) * qJD(2), t152 * t239, -pkin(4) * t223 - t266 + t45 * t152 - t66 * t128 + (pkin(9) * t255 + t265) * qJD(5) + (-t30 * t170 + (-pkin(9) * t238 - t174 * t62) * t169) * qJD(2), pkin(4) * t100 - t66 * t130 - t46 * t152 + t268 + (-pkin(9) * t256 + t271) * qJD(5) + (-t62 * t246 + (-pkin(9) * t229 + t31) * t170) * qJD(2), t46 * t128 + t45 * t130 + ((t130 * qJD(5) - t223) * pkin(9) + t291) * t173 + (pkin(9) * t206 + t290) * t169, -t42 * pkin(4) + pkin(9) * t179 - t30 * t45 - t31 * t46 - t62 * t66, -t33 * t132 + t193 * t262, t33 * t131 + t132 * t180 - t193 * t263 + t262 * t76, t262 * t148 + (qJD(4) * t132 + t193) * t239, -t131 * t180 - t263 * t76, -t263 * t148 + (-qJD(4) * t131 + t76) * t239, t148 * t239, t29 * t131 - t160 * t180 + t201 * t76 - t263 * t49 + t278 * t148 + (qJD(4) * t96 - t10) * t239, t29 * t132 - t160 * t33 - t201 * t193 - t262 * t49 - t277 * t148 + (-qJD(4) * t97 + t11) * t239, -t1 * t131 + t262 * t10 + t263 * t11 - t2 * t132 + t180 * t97 - t193 * t278 + t277 * t76 + t96 * t33, t1 * t97 - t278 * t10 - t277 * t11 + t29 * t160 + t2 * t96 + t201 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, -t128 ^ 2 + t130 ^ 2, -t100 - t259, -t258, -t223 - t257, t156, -t62 * t130 - t290, t62 * t128 - t291, 0, 0, -t279, t296, t294, t279, t292, t156, t12 * t148 + (-t130 * t76 + t148 * t233 + t156 * t172) * pkin(5) + t295, -t13 * t148 + (t130 * t193 + t148 * t232 - t156 * t168) * pkin(5) + t293, -t10 * t76 - t11 * t193 - t12 * t193 + t13 * t76 + (t168 * t180 + t172 * t33 + (-t168 * t193 - t172 * t76) * qJD(6)) * pkin(5), -t10 * t12 - t11 * t13 + (t1 * t168 - t130 * t49 + t172 * t2 + (-t10 * t168 + t11 * t172) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, t296, t294, t279, t292, t156, -t11 * t148 + t295, -t10 * t148 + t293, 0, 0;];
tauc_reg  = t5;
