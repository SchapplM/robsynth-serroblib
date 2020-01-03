% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:17
% EndTime: 2019-12-31 21:17:29
% DurationCPUTime: 4.28s
% Computational Cost: add. (7602->381), mult. (18923->521), div. (0->0), fcn. (13618->8), ass. (0->214)
t201 = sin(qJ(2));
t200 = sin(qJ(3));
t203 = cos(qJ(2));
t262 = t200 * t203;
t289 = cos(qJ(3));
t167 = t289 * t201 + t262;
t155 = qJD(1) * t167;
t194 = qJD(2) + qJD(3);
t197 = sin(pkin(9));
t198 = cos(pkin(9));
t136 = t198 * t155 + t197 * t194;
t199 = sin(qJ(5));
t202 = cos(qJ(5));
t144 = t197 * t155;
t295 = t198 * t194 - t144;
t218 = t202 * t295;
t84 = -t199 * t136 + t218;
t307 = t84 ^ 2;
t244 = t289 * t203;
t182 = qJD(1) * t244;
t253 = qJD(1) * t201;
t241 = t200 * t253;
t153 = -t182 + t241;
t149 = qJD(5) + t153;
t306 = t84 * t149;
t83 = t202 * t136 + t199 * t295;
t305 = t83 ^ 2;
t261 = t202 * t198;
t292 = -t199 * t197 + t261;
t165 = t202 * t197 + t199 * t198;
t152 = t165 * qJD(5);
t279 = t165 * t153 + t152;
t293 = t292 * qJD(5);
t303 = t292 * t153 + t293;
t263 = t200 * t201;
t227 = t194 * t263;
t255 = t194 * t182;
t119 = qJD(1) * t227 - t255;
t268 = t197 * t119;
t302 = t136 * t153 - t268;
t267 = t198 * t119;
t301 = t153 * t295 - t267;
t282 = qJD(2) * pkin(2);
t248 = t201 * t282;
t300 = 0.2e1 * t248;
t251 = qJD(1) * qJD(2);
t299 = -0.2e1 * t251;
t131 = t194 * t167;
t120 = t131 * qJD(1);
t239 = t201 * t251;
t50 = pkin(2) * t239 + t120 * pkin(3) + t119 * qJ(4) - t155 * qJD(4);
t290 = -pkin(7) - pkin(6);
t176 = t290 * t201;
t169 = qJD(1) * t176;
t161 = t169 + t282;
t245 = qJD(2) * t290;
t233 = qJD(1) * t245;
t162 = t201 * t233;
t163 = t203 * t233;
t177 = t290 * t203;
t171 = qJD(1) * t177;
t240 = t289 * qJD(3);
t252 = qJD(3) * t200;
t234 = -t161 * t240 - t289 * t162 - t200 * t163 - t171 * t252;
t68 = t194 * qJD(4) - t234;
t27 = -t197 * t68 + t198 * t50;
t28 = t197 * t50 + t198 * t68;
t224 = -t27 * t197 + t28 * t198;
t166 = -t244 + t263;
t91 = t120 * t166;
t298 = t153 * t131 + t91;
t297 = t119 * t165;
t157 = t289 * t171;
t127 = t200 * t169 - t157;
t296 = -pkin(2) * t252 + t127;
t294 = t289 * t176 + t200 * t177;
t33 = t199 * (qJD(5) * t136 - t268) - qJD(5) * t218 + t119 * t261;
t13 = t120 * pkin(4) + pkin(8) * t267 + t27;
t16 = pkin(8) * t268 + t28;
t190 = -t203 * pkin(2) - pkin(1);
t175 = qJD(1) * t190;
t101 = t153 * pkin(3) - t155 * qJ(4) + t175;
t125 = t200 * t161 - t157;
t108 = t194 * qJ(4) + t125;
t55 = t198 * t101 - t197 * t108;
t35 = t153 * pkin(4) - t136 * pkin(8) + t55;
t56 = t197 * t101 + t198 * t108;
t38 = pkin(8) * t295 + t56;
t222 = t199 * t38 - t202 * t35;
t3 = -qJD(5) * t222 + t199 * t13 + t202 * t16;
t291 = t153 ^ 2;
t288 = t198 * pkin(4);
t191 = t198 * pkin(8);
t287 = t83 * t84;
t186 = t200 * pkin(2) + qJ(4);
t159 = (-pkin(8) - t186) * t197;
t160 = t198 * t186 + t191;
t117 = t199 * t159 + t202 * t160;
t180 = pkin(2) * t240 + qJD(4);
t270 = t153 * t198;
t229 = t155 * pkin(4) + pkin(8) * t270;
t121 = t155 * pkin(3) + t153 * qJ(4);
t247 = pkin(2) * t253;
t106 = t121 + t247;
t156 = t200 * t171;
t128 = t289 * t169 + t156;
t66 = t198 * t106 - t197 * t128;
t41 = t229 + t66;
t271 = t153 * t197;
t250 = pkin(8) * t271;
t67 = t197 * t106 + t198 * t128;
t51 = t250 + t67;
t286 = qJD(5) * t117 + t165 * t180 - t199 * t51 + t202 * t41;
t116 = t202 * t159 - t199 * t160;
t285 = -qJD(5) * t116 - t180 * t292 + t199 * t41 + t202 * t51;
t173 = (-pkin(8) - qJ(4)) * t197;
t174 = t198 * qJ(4) + t191;
t134 = t199 * t173 + t202 * t174;
t124 = t289 * t161 + t156;
t69 = t198 * t121 - t197 * t124;
t44 = t229 + t69;
t70 = t197 * t121 + t198 * t124;
t53 = t250 + t70;
t284 = qJD(4) * t165 + qJD(5) * t134 - t199 * t53 + t202 * t44;
t133 = t202 * t173 - t199 * t174;
t283 = -qJD(4) * t292 - qJD(5) * t133 + t199 * t44 + t202 * t53;
t232 = qJD(2) * t244;
t130 = -t203 * t240 + t227 - t232;
t64 = t131 * pkin(3) + t130 * qJ(4) - t167 * qJD(4) + t248;
t170 = t201 * t245;
t85 = t294 * qJD(3) + t289 * t170 + t245 * t262;
t32 = t197 * t64 + t198 * t85;
t235 = t200 * t162 - t289 * t163;
t73 = qJD(3) * t125 + t235;
t280 = t73 * t294;
t278 = t119 * t167;
t277 = t130 * t197;
t276 = t130 * t198;
t275 = t136 * t197;
t274 = t149 * t155;
t272 = t153 * t155;
t269 = t167 * t197;
t205 = qJD(1) ^ 2;
t260 = t203 * t205;
t204 = qJD(2) ^ 2;
t259 = t204 * t201;
t258 = t204 * t203;
t123 = t166 * pkin(3) - t167 * qJ(4) + t190;
t139 = t200 * t176 - t289 * t177;
t78 = t197 * t123 + t198 * t139;
t254 = t201 ^ 2 - t203 ^ 2;
t246 = t201 * t260;
t238 = t56 * t155 + t73 * t197;
t31 = -t197 * t85 + t198 * t64;
t237 = pkin(1) * t299;
t77 = t198 * t123 - t197 * t139;
t189 = -t289 * pkin(2) - pkin(3);
t231 = t203 * t239;
t142 = pkin(4) * t271;
t230 = t142 - t296;
t225 = -t55 * t155 - t73 * t198;
t223 = -t197 * t55 + t198 * t56;
t15 = t199 * t35 + t202 * t38;
t52 = t166 * pkin(4) - t167 * t191 + t77;
t63 = -pkin(8) * t269 + t78;
t24 = -t199 * t63 + t202 * t52;
t25 = t199 * t52 + t202 * t63;
t221 = t119 * t294 + t73 * t167;
t220 = -t55 * t270 - t56 * t271 + t224;
t219 = t198 * t295;
t216 = -t175 * t155 - t235;
t45 = -pkin(4) * t268 + t73;
t105 = -t194 * pkin(3) + qJD(4) - t124;
t74 = -pkin(4) * t295 + t105;
t215 = t222 * t155 + t279 * t74 - t292 * t45;
t214 = t15 * t155 + t45 * t165 + t303 * t74;
t213 = t175 * t153 + t234;
t212 = -t105 * t130 + t221;
t211 = t119 * t166 - t167 * t120 + t130 * t153;
t4 = -qJD(5) * t15 + t202 * t13 - t199 * t16;
t210 = -t279 * t15 - t4 * t165 + t303 * t222 + t292 * t3;
t209 = t219 - t275;
t208 = pkin(3) * t119 - qJ(4) * t120 + (-qJD(4) + t105) * t153;
t207 = -t189 * t119 - t186 * t120 + (t105 - t180) * t153;
t86 = qJD(3) * t139 + t200 * t170 - t290 * t232;
t34 = qJD(5) * t83 - t297;
t193 = t198 ^ 2;
t192 = t197 ^ 2;
t187 = -pkin(3) - t288;
t172 = t189 - t288;
t112 = t292 * t167;
t111 = t165 * t167;
t103 = pkin(4) * t269 - t294;
t92 = t155 ^ 2 - t291;
t89 = -t142 + t125;
t87 = t255 + (t153 - t241) * t194;
t58 = t301 * t197;
t57 = t302 * t198;
t54 = -pkin(4) * t277 + t86;
t43 = -t130 * t165 + t293 * t167;
t42 = t292 * t130 + t152 * t167;
t40 = t198 * t120 - t155 * t295 - t197 * t291;
t39 = t197 * t120 - t136 * t155 + t198 * t291;
t30 = t209 * t153 + (t192 - t193) * t119;
t29 = pkin(8) * t277 + t32;
t23 = t131 * pkin(4) + pkin(8) * t276 + t31;
t20 = t120 * t292 - t279 * t149 - t155 * t84;
t19 = t165 * t120 + t149 * t303 - t83 * t155;
t9 = -t279 * t84 - t292 * t34;
t8 = -t33 * t165 + t303 * t83;
t7 = -qJD(5) * t25 - t199 * t29 + t202 * t23;
t6 = qJD(5) * t24 + t199 * t23 + t202 * t29;
t5 = -t165 * t34 - t279 * t83 - t292 * t33 + t303 * t84;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t231, t254 * t299, t258, -0.2e1 * t231, -t259, 0, -pkin(6) * t258 + t201 * t237, pkin(6) * t259 + t203 * t237, 0, 0, -t155 * t130 - t278, -t155 * t131 + t211, -t130 * t194, t298, -t131 * t194, 0, t190 * t120 + t175 * t131 - t86 * t194 + (qJD(1) * t166 + t153) * t248, -t190 * t119 - t175 * t130 + t155 * t300 - t85 * t194, -t139 * t120 + t124 * t130 - t125 * t131 - t85 * t153 + t86 * t155 + t166 * t234 + t221, -t124 * t86 + t125 * t85 - t139 * t234 + t175 * t300 - t280, -t136 * t276 - t193 * t278, -t130 * t209 + 0.2e1 * t267 * t269, t136 * t131 - t198 * t211, -t192 * t278 + t277 * t295, t131 * t295 + t197 * t211, t298, t77 * t120 + t55 * t131 + t31 * t153 + t27 * t166 + t197 * t212 - t295 * t86, -t78 * t120 - t56 * t131 + t86 * t136 - t32 * t153 - t28 * t166 + t198 * t212, -t31 * t136 - t32 * t144 + (t77 * t119 + t55 * t130 - t27 * t167 + t32 * t194) * t198 + (t78 * t119 + t56 * t130 - t28 * t167) * t197, t105 * t86 + t27 * t77 + t28 * t78 + t55 * t31 + t56 * t32 - t280, -t33 * t112 - t83 * t42, t33 * t111 - t112 * t34 - t42 * t84 - t83 * t43, t112 * t120 + t83 * t131 - t42 * t149 - t33 * t166, t34 * t111 - t43 * t84, -t111 * t120 + t131 * t84 - t43 * t149 - t34 * t166, t149 * t131 + t91, t103 * t34 + t45 * t111 + t24 * t120 - t131 * t222 + t7 * t149 + t4 * t166 + t74 * t43 - t54 * t84, -t103 * t33 + t45 * t112 - t25 * t120 - t15 * t131 - t6 * t149 - t3 * t166 - t74 * t42 + t54 * t83, -t3 * t111 - t4 * t112 - t15 * t43 - t222 * t42 + t24 * t33 - t25 * t34 + t6 * t84 - t7 * t83, t45 * t103 + t15 * t6 - t222 * t7 + t4 * t24 + t3 * t25 + t74 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, t254 * t205, 0, t246, 0, 0, t205 * pkin(1) * t201, pkin(1) * t260, 0, 0, t272, t92, t87, -t272, 0, 0, -t153 * t247 + t127 * t194 + (t157 + (-pkin(2) * t194 - t161) * t200) * qJD(3) + t216, t128 * t194 + (-t155 * t253 - t194 * t240) * pkin(2) + t213, (t125 - t127) * t155 + (-t124 + t128) * t153 + (t289 * t119 - t120 * t200 + (-t289 * t153 + t155 * t200) * qJD(3)) * pkin(2), t124 * t127 - t125 * t128 + (-t175 * t253 - t289 * t73 - t200 * t234 + (-t124 * t200 + t289 * t125) * qJD(3)) * pkin(2), t57, t30, t39, -t58, t40, -t272, -t66 * t153 + t197 * t207 + t295 * t296 + t225, -t136 * t296 + t67 * t153 + t198 * t207 + t238, t180 * t219 - t67 * t295 + (t197 * t180 + t66) * t136 + t220, -t105 * t296 + t180 * t223 + t186 * t224 + t73 * t189 - t55 * t66 - t56 * t67, t8, t5, t19, t9, t20, -t274, t116 * t120 - t286 * t149 + t172 * t34 - t230 * t84 + t215, -t117 * t120 + t285 * t149 - t172 * t33 + t230 * t83 + t214, t116 * t33 - t117 * t34 - t285 * t84 + t286 * t83 + t210, t4 * t116 + t3 * t117 - t285 * t15 + t45 * t172 + t222 * t286 + t230 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t92, t87, -t272, 0, 0, t216 + (-qJD(3) + t194) * t125, t124 * t194 + t213, 0, 0, t57, t30, t39, -t58, t40, -t272, t125 * t295 - t69 * t153 + t197 * t208 + t225, -t125 * t136 + t70 * t153 + t198 * t208 + t238, t69 * t136 - t70 * t295 + (t219 + t275) * qJD(4) + t220, -t73 * pkin(3) + qJ(4) * t224 + qJD(4) * t223 - t105 * t125 - t55 * t69 - t56 * t70, t8, t5, t19, t9, t20, -t274, t133 * t120 - t284 * t149 + t187 * t34 + t84 * t89 + t215, -t134 * t120 + t283 * t149 - t187 * t33 - t89 * t83 + t214, t133 * t33 - t134 * t34 - t283 * t84 + t284 * t83 + t210, t4 * t133 + t3 * t134 - t283 * t15 + t45 * t187 + t222 * t284 - t74 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t301, -t136 ^ 2 - t295 ^ 2, t136 * t55 - t295 * t56 + t73, 0, 0, 0, 0, 0, 0, t83 * t149 + t34, -t33 + t306, -t305 - t307, -t15 * t84 - t222 * t83 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t305 - t307, -t33 - t306, t287, t297 + (-qJD(5) + t149) * t83, t120, t15 * t149 - t74 * t83 + t4, -t149 * t222 - t74 * t84 - t3, 0, 0;];
tauc_reg = t1;
