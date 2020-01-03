% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:43
% EndTime: 2019-12-31 20:30:51
% DurationCPUTime: 4.44s
% Computational Cost: add. (5080->485), mult. (10786->592), div. (0->0), fcn. (6952->8), ass. (0->245)
t161 = sin(qJ(5));
t165 = cos(qJ(5));
t167 = cos(qJ(2));
t259 = qJD(1) * t167;
t163 = sin(qJ(2));
t260 = qJD(1) * t163;
t78 = -qJD(1) * pkin(1) - pkin(2) * t259 - qJ(3) * t260;
t57 = pkin(3) * t259 - t78;
t166 = cos(qJ(4));
t271 = t166 * t167;
t162 = sin(qJ(4));
t277 = t162 * t163;
t84 = t271 + t277;
t73 = t84 * qJD(1);
t76 = -t162 * t259 + t166 * t260;
t24 = pkin(4) * t73 - pkin(8) * t76 + t57;
t248 = qJD(2) - qJD(4);
t302 = pkin(2) + pkin(3);
t241 = t302 * qJD(2);
t138 = pkin(6) * t260;
t90 = pkin(7) * t260 - t138;
t320 = qJD(3) - t90;
t61 = -t241 + t320;
t156 = qJD(2) * qJ(3);
t139 = pkin(6) * t259;
t216 = -pkin(7) * t259 + t139;
t77 = t156 + t216;
t34 = t162 * t61 + t166 * t77;
t31 = -t248 * pkin(8) + t34;
t11 = -t161 * t31 + t165 * t24;
t319 = qJD(5) + t73;
t324 = t11 * t319;
t12 = t161 * t24 + t165 * t31;
t323 = t319 * t12;
t224 = t165 * t248;
t49 = t161 * t76 + t224;
t322 = t319 * t49;
t51 = -t161 * t248 + t165 * t76;
t321 = t319 * t51;
t142 = t163 * qJDD(1);
t251 = qJD(1) * qJD(2);
t239 = t167 * t251;
t318 = -t239 - t142;
t254 = qJD(4) * t166;
t256 = qJD(2) * t167;
t187 = t162 * t256 + t163 * t254;
t250 = qJD(1) * qJD(4);
t238 = t167 * t250;
t240 = t163 * t251;
t28 = t187 * qJD(1) + t84 * qJDD(1) - t162 * t238 - t166 * t240;
t26 = qJDD(5) + t28;
t285 = t165 * t26;
t310 = t161 * t319;
t317 = t310 * t319 - t49 * t76 - t285;
t288 = t161 * t26;
t309 = t165 * t319;
t316 = t309 * t319 - t51 * t76 + t288;
t249 = t167 * qJDD(1);
t191 = -t162 * t249 - t250 * t277 + (t142 - t238) * t166;
t192 = t84 * qJD(2);
t174 = qJD(1) * t192 + t191;
t247 = qJDD(2) - qJDD(4);
t280 = qJD(5) * t51;
t14 = t161 * t174 + t165 * t247 + t280;
t289 = t14 * t165;
t315 = t310 * t49 - t289;
t253 = qJD(5) * t161;
t13 = qJD(5) * t224 + t161 * t247 - t165 * t174 + t76 * t253;
t290 = t13 * t161;
t314 = t309 * t51 - t290;
t313 = (t13 + t322) * t165 + (t14 + t321) * t161;
t265 = t167 * pkin(2) + t163 * qJ(3);
t312 = -pkin(1) - t265;
t44 = t162 * t216 + t166 * t90;
t92 = -qJ(3) * t162 - t166 * t302;
t62 = qJD(3) * t166 + qJD(4) * t92;
t292 = t62 - t44;
t93 = t166 * qJ(3) - t162 * t302;
t291 = t93 * qJD(4) + t320 * t162 + t166 * t216;
t33 = -t162 * t77 + t166 * t61;
t30 = t248 * pkin(4) - t33;
t311 = t319 * t30;
t203 = -t76 * t248 - t28;
t168 = cos(qJ(1));
t153 = g(1) * t168;
t164 = sin(qJ(1));
t308 = g(2) * t164 + t153;
t307 = t167 * t162 - t163 * t166;
t282 = pkin(6) * qJDD(2);
t304 = qJD(2) * (qJD(1) * t312 + t78) - t282;
t303 = t248 ^ 2;
t301 = pkin(6) - pkin(7);
t67 = t84 * t164;
t300 = g(1) * t67;
t299 = g(1) * t164;
t298 = g(2) * t168;
t297 = t11 * t76;
t296 = t12 * t76;
t147 = t167 * pkin(3);
t295 = t51 * t49;
t294 = t319 * t76;
t293 = t76 * t73;
t287 = t161 * t49;
t286 = t161 * t51;
t284 = t165 * t49;
t283 = t165 * t51;
t281 = qJD(4) * t319;
t279 = qJD(5) * t319;
t159 = qJDD(1) * pkin(1);
t278 = qJDD(2) * pkin(2);
t276 = t163 * t164;
t274 = t163 * t168;
t171 = qJD(1) ^ 2;
t273 = t163 * t171;
t272 = t164 * t167;
t269 = t167 * t168;
t143 = t163 * qJD(3);
t266 = qJ(3) * t256 + t143;
t264 = t168 * pkin(1) + t164 * pkin(6);
t157 = t163 ^ 2;
t158 = t167 ^ 2;
t262 = -t157 + t158;
t261 = t157 + t158;
t258 = qJD(2) * t163;
t257 = qJD(2) * t166;
t255 = qJD(4) * t162;
t252 = qJD(5) * t165;
t246 = t73 ^ 2 - t76 ^ 2;
t101 = t301 * t167;
t245 = t307 * t253;
t244 = t307 * t252;
t243 = t147 + t265;
t242 = -g(1) * t274 - g(2) * t276 + g(3) * t167;
t133 = pkin(6) * t142;
t237 = pkin(6) * t239 + qJDD(3) + t133;
t149 = t168 * pkin(6);
t236 = -pkin(7) * t168 + t149;
t41 = t318 * pkin(7) - t302 * qJDD(2) + t237;
t134 = pkin(6) * t249;
t154 = qJDD(2) * qJ(3);
t155 = qJD(2) * qJD(3);
t58 = -pkin(6) * t240 + t134 + t154 + t155;
t42 = (t240 - t249) * pkin(7) + t58;
t235 = t162 * t42 - t166 * t41 + t77 * t254 + t61 * t255;
t229 = -qJD(2) * pkin(2) + qJD(3);
t228 = t33 * t248;
t227 = t34 * t248;
t226 = t73 * t248;
t223 = pkin(2) * t269 + qJ(3) * t274 + t264;
t222 = -t133 - t242;
t221 = qJD(2) * t101;
t220 = t163 * t241;
t219 = t163 * t239;
t40 = pkin(4) * t76 + pkin(8) * t73;
t66 = t307 * t164;
t68 = t162 * t269 - t166 * t274;
t218 = g(1) * t66 - g(2) * t68;
t217 = g(2) * t67 - g(3) * t307;
t215 = t261 * qJDD(1) * pkin(6);
t170 = qJD(2) ^ 2;
t214 = pkin(6) * t170 + t298;
t46 = -qJD(4) * t84 + t192;
t212 = -t26 * t307 + t319 * t46;
t209 = -qJD(5) * t31 - t298;
t208 = -t11 * t165 - t12 * t161;
t207 = t11 * t161 - t12 * t165;
t188 = pkin(4) * t84 + pkin(8) * t307 + t243;
t32 = pkin(1) + t188;
t100 = t301 * t163;
t53 = t100 * t162 + t101 * t166;
t20 = -t161 * t53 + t165 * t32;
t21 = t161 * t32 + t165 * t53;
t94 = t138 + t229;
t99 = t139 + t156;
t206 = t163 * t99 - t167 * t94;
t205 = pkin(2) * t249 - t318 * qJ(3) + qJD(1) * t143 + t159;
t204 = pkin(3) * t269 + t223;
t202 = t166 * t100 - t101 * t162;
t201 = -qJD(5) * t30 * t307 + t300;
t130 = qJ(3) * t259;
t65 = -t302 * t260 + t130;
t64 = t237 - t278;
t199 = pkin(3) * t249 + t205;
t198 = -pkin(8) * t26 + t311;
t197 = -0.2e1 * pkin(1) * t251 - t282;
t196 = -t162 * t41 - t166 * t42 - t61 * t254 + t77 * t255;
t195 = g(1) * t68 + g(2) * t66 + g(3) * t84;
t69 = t84 * t168;
t194 = g(1) * t69 + t217;
t56 = -t220 + t266;
t8 = t247 * pkin(4) + t235;
t193 = t30 * t46 - t307 * t8 + t153;
t190 = -t214 + 0.2e1 * t159;
t189 = t195 - t8;
t7 = -t247 * pkin(8) - t196;
t186 = -qJD(5) * t24 + t217 - t7;
t89 = -pkin(8) + t93;
t185 = -t26 * t89 - t319 * t62 - t311;
t184 = t308 * t302 * t163;
t183 = pkin(8) * t279 - t189;
t182 = t89 * t279 + t189;
t39 = pkin(2) * t240 - t205;
t70 = pkin(2) * t258 - t266;
t181 = -qJD(1) * t70 - qJDD(1) * t312 - t214 - t39;
t6 = -t191 * pkin(8) + t28 * pkin(4) + (-pkin(8) * t271 + (-pkin(8) * t162 - t302) * t163) * t251 + t199;
t1 = qJD(5) * t11 + t161 * t6 + t165 * t7;
t5 = t165 * t6;
t2 = -qJD(5) * t12 - t161 * t7 + t5;
t180 = t208 * qJD(5) + t1 * t165 - t2 * t161;
t179 = -qJD(2) * t206 + t64 * t163 + t58 * t167;
t178 = -t57 * t76 + t195 - t235;
t177 = (-g(1) * (t312 - t147) + g(2) * pkin(7)) * t164;
t176 = t57 * t73 + t194 + t196;
t175 = -qJD(2) * t73 - t191;
t124 = g(1) * t272;
t117 = qJ(3) * t269;
t115 = qJ(3) * t272;
t112 = t167 * t273;
t98 = t262 * t171;
t97 = qJDD(2) * t167 - t163 * t170;
t96 = qJDD(2) * t163 + t167 * t170;
t91 = t301 * t258;
t88 = pkin(4) - t92;
t87 = pkin(2) * t260 - t130;
t83 = qJDD(1) * t158 - 0.2e1 * t219;
t82 = qJDD(1) * t157 + 0.2e1 * t219;
t81 = pkin(1) + t243;
t75 = t161 * t260 + t165 * t257;
t72 = -t161 * t257 + t165 * t260;
t59 = t163 * t249 + t262 * t251;
t48 = -t161 * t164 + t165 * t69;
t47 = -t161 * t69 - t164 * t165;
t45 = -t163 * t257 - t167 * t255 + t187;
t29 = -qJD(1) * t220 + t199;
t27 = -t40 + t65;
t23 = t53 * qJD(4) - t162 * t91 - t166 * t221;
t22 = t202 * qJD(4) + t162 * t221 - t166 * t91;
t19 = t161 * t40 + t165 * t33;
t18 = -t161 * t33 + t165 * t40;
t17 = pkin(4) * t45 - pkin(8) * t46 + t56;
t16 = t161 * t27 + t165 * t44;
t15 = -t161 * t44 + t165 * t27;
t4 = -qJD(5) * t21 - t161 * t22 + t165 * t17;
t3 = qJD(5) * t20 + t161 * t17 + t165 * t22;
t9 = [0, 0, 0, 0, 0, qJDD(1), -t298 + t299, t308, 0, 0, t82, 0.2e1 * t59, t96, t83, t97, 0, t163 * t197 + t167 * t190 + t124, t197 * t167 + (-t190 - t299) * t163, 0.2e1 * t215 - t308, -g(1) * (-pkin(1) * t164 + t149) - g(2) * t264 + (pkin(6) ^ 2 * t261 + pkin(1) ^ 2) * qJDD(1), t82, t96, -0.2e1 * t59, 0, -t97, t83, t304 * t163 + t181 * t167 + t124, t215 + t179 - t308, -t304 * t167 + (t181 + t299) * t163, pkin(6) * t179 - g(1) * t149 - g(2) * t223 + t78 * t70 + (-t299 + t39) * t312, -t174 * t307 + t76 * t46, -t174 * t84 + t28 * t307 - t76 * t45 - t46 * t73, t247 * t307 - t248 * t46, t28 * t84 + t45 * t73, t247 * t84 + t248 * t45, 0, -g(2) * t69 - t202 * t247 + t23 * t248 + t81 * t28 + t29 * t84 + t57 * t45 + t56 * t73 + t300, t174 * t81 + t22 * t248 + t247 * t53 - t29 * t307 + t57 * t46 + t56 * t76 - t218, t175 * t202 + t196 * t84 - t22 * t73 + t23 * t76 - t235 * t307 - t53 * t28 - t33 * t46 - t34 * t45 + t308, -g(1) * t236 - g(2) * t204 - t196 * t53 - t202 * t235 + t34 * t22 - t33 * t23 + t29 * t81 + t57 * t56 + t177, t51 * t245 + (t13 * t307 + t46 * t51) * t165, (-t284 - t286) * t46 - (t290 - t289 + (-t283 + t287) * qJD(5)) * t307, -t13 * t84 + t165 * t212 + t245 * t319 + t45 * t51, -t49 * t244 + (-t14 * t307 + t46 * t49) * t161, -t14 * t84 - t161 * t212 + t244 * t319 - t45 * t49, t26 * t84 + t319 * t45, -g(2) * t48 + t11 * t45 - t14 * t202 + t161 * t193 + t165 * t201 + t2 * t84 + t20 * t26 + t23 * t49 + t319 * t4, -g(2) * t47 - t1 * t84 - t12 * t45 + t13 * t202 - t161 * t201 + t165 * t193 - t21 * t26 + t23 * t51 - t3 * t319, t13 * t20 - t14 * t21 - t3 * t49 - t4 * t51 + t208 * t46 - (qJD(5) * t207 - t1 * t161 - t165 * t2) * t307 + t218, t1 * t21 + t12 * t3 + t2 * t20 + t11 * t4 - t8 * t202 + t30 * t23 - g(1) * (-pkin(4) * t67 - pkin(8) * t66 + t236) - g(2) * (pkin(4) * t69 + pkin(8) * t68 + t204) + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t98, t142, t112, t249, qJDD(2), pkin(1) * t273 + t222, g(3) * t163 - t134 + (pkin(1) * t171 + t308) * t167, 0, 0, -t112, t142, t98, qJDD(2), -t249, t112, 0.2e1 * t278 - qJDD(3) + (-t163 * t78 + t167 * t87) * qJD(1) + t222, (-t163 * pkin(2) + qJ(3) * t167) * qJDD(1) + ((t99 - t156) * t163 + (t229 - t94) * t167) * qJD(1), t134 + 0.2e1 * t154 + 0.2e1 * t155 + (qJD(1) * t87 - g(3)) * t163 + (qJD(1) * t78 - t308) * t167, t58 * qJ(3) + t99 * qJD(3) - t64 * pkin(2) - t78 * t87 - g(1) * (-pkin(2) * t274 + t117) - g(2) * (-pkin(2) * t276 + t115) - g(3) * t265 + t206 * qJD(1) * pkin(6), -t293, t246, t226 + t175, t293, -t203, t247, -t247 * t92 + t291 * t248 - t65 * t73 - t178, t247 * t93 + t292 * t248 - t65 * t76 - t176, -t93 * t28 + t92 * t175 + (-t34 + t291) * t76 + (t33 - t292) * t73, -g(1) * t117 - g(2) * t115 - g(3) * t243 - t196 * t93 - t235 * t92 - t291 * t33 + t292 * t34 - t57 * t65 + t184, -t314, t313, -t316, -t315, t317, t294, t14 * t88 - t15 * t319 + t161 * t185 - t165 * t182 + t291 * t49 + t297, -t13 * t88 + t16 * t319 + t161 * t182 + t165 * t185 + t291 * t51 - t296, t15 * t51 + t16 * t49 + (t11 * t73 - t14 * t89 - t49 * t62 - t1 + (t51 * t89 + t11) * qJD(5)) * t165 + (t12 * t73 - t13 * t89 + t51 * t62 + t2 + (t49 * t89 + t12) * qJD(5)) * t161 + t194, t8 * t88 - g(1) * (pkin(4) * t68 - pkin(8) * t69 + t117) - g(2) * (pkin(4) * t66 - pkin(8) * t67 + t115) - g(3) * t188 + t291 * t30 + t184 + (t165 * t62 - t16) * t12 + (-t161 * t62 - t15) * t11 + t180 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t112, t142, -t157 * t171 - t170, -qJD(2) * t99 + t260 * t78 + t242 + t64, 0, 0, 0, 0, 0, 0, -t162 * t303 - t166 * t247 - t73 * t260, t162 * t247 - t166 * t303 - t76 * t260, t203 * t162 + (-qJD(4) * t73 - t191) * t166, -t57 * t260 + (-t235 - t227) * t166 + (-t196 + t228) * t162 + t242, 0, 0, 0, 0, 0, 0, -t319 * t72 + (-t161 * t281 - t14) * t166 + (-t248 * t49 - t252 * t319 - t288) * t162, t319 * t75 + (-t165 * t281 + t13) * t166 + (-t248 * t51 + t253 * t319 - t285) * t162, t49 * t75 + t51 * t72 + (-t284 + t286) * t254 + (-t290 - t289 + (t283 + t287) * qJD(5)) * t162, -t11 * t72 - t12 * t75 + (-qJD(4) * t207 - t8) * t166 + (-t248 * t30 + t180) * t162 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, -t246, -t226 + t174, -t293, t203, -t247, -t227 + t178, -t228 + t176, 0, 0, t314, -t313, t316, t315, -t317, -t294, -pkin(4) * t14 + t161 * t198 - t165 * t183 - t18 * t319 - t34 * t49 - t297, pkin(4) * t13 + t161 * t183 + t165 * t198 + t19 * t319 - t34 * t51 + t296, t18 * t51 + t19 * t49 + (t1 - t324 + (-t14 + t280) * pkin(8)) * t165 + (-t2 - t323 + (qJD(5) * t49 - t13) * pkin(8)) * t161 - t194, -t11 * t18 - t12 * t19 - t30 * t34 + t189 * pkin(4) + (t180 - t194) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t49 ^ 2 + t51 ^ 2, -t13 + t322, -t295, -t14 + t321, t26, -g(1) * t47 + t161 * t186 + t165 * t209 - t30 * t51 + t323 + t5, g(1) * t48 + t324 + t30 * t49 + (-t209 - t6) * t161 + t186 * t165, 0, 0;];
tau_reg = t9;
