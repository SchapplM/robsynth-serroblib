% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR13_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:19
% DurationCPUTime: 1.61s
% Computational Cost: add. (7481->314), mult. (13868->443), div. (0->0), fcn. (8881->10), ass. (0->220)
t161 = qJD(1) + qJD(2);
t173 = cos(qJ(2));
t264 = pkin(1) * qJD(1);
t231 = t173 * t264;
t136 = t161 * pkin(2) + t231;
t172 = cos(qJ(3));
t168 = sin(qJ(3));
t169 = sin(qJ(2));
t232 = t169 * t264;
t205 = t168 * t232;
t109 = t172 * t136 - t205;
t110 = t136 * t168 + t172 * t232;
t171 = cos(qJ(4));
t165 = cos(pkin(5));
t238 = qJD(4) * t171;
t220 = t165 * t238;
t167 = sin(qJ(4));
t250 = t165 * t167;
t287 = pkin(3) * t220 - t109 * t171 + t110 * t250;
t247 = t169 * t172;
t185 = -t168 * t173 - t247;
t120 = t185 * t264;
t248 = t168 * t169;
t184 = t172 * t173 - t248;
t121 = t184 * t264;
t249 = t165 * t171;
t286 = -t120 * t249 + t121 * t167 + (-t167 * t172 - t168 * t249) * qJD(3) * pkin(2);
t154 = pkin(2) * t172 + pkin(3);
t240 = qJD(3) * t172;
t285 = -t120 * t250 + t154 * t220 + (pkin(2) * t240 - t121) * t171;
t241 = qJD(3) * t168;
t230 = pkin(2) * t241;
t206 = t165 * t230;
t164 = sin(pkin(5));
t158 = t164 * pkin(9);
t142 = pkin(2) * t168 + t158;
t274 = pkin(10) * t164;
t217 = -t142 - t274;
t284 = (t217 * qJD(4) - t206) * t167 + t285;
t134 = t154 * t250;
t283 = (t217 * t171 - t134) * qJD(4) + t286;
t228 = t164 * (-pkin(9) - pkin(10));
t204 = t167 * t228;
t282 = qJD(4) * t204 + t287;
t151 = pkin(3) * t250;
t59 = -t109 * t167 - t110 * t249;
t281 = -(t171 * t228 - t151) * qJD(4) + t59;
t195 = t120 + t230;
t170 = cos(qJ(5));
t221 = t164 * t238;
t251 = t164 * t171;
t225 = t170 * t251;
t280 = -qJD(5) * t225 - t170 * t221;
t104 = t171 * t142 + t134;
t279 = qJD(4) * t104;
t166 = sin(qJ(5));
t186 = t166 * t171 + t167 * t170;
t114 = t186 * t164;
t157 = qJD(3) + t161;
t256 = t157 * t164;
t92 = pkin(9) * t256 + t110;
t197 = pkin(10) * t256 + t92;
t182 = t197 * t167;
t203 = qJD(2) * t231;
t244 = (qJD(2) + qJD(3)) * t205;
t76 = (qJD(3) * t136 + t203) * t172 - t244;
t175 = (t185 * qJD(2) - t169 * t240) * pkin(1);
t174 = qJD(1) * t175;
t77 = -t136 * t241 + t174;
t276 = pkin(3) * t157;
t94 = t109 + t276;
t234 = t171 * t76 + t94 * t220 + t77 * t250;
t14 = -qJD(4) * t182 + t234;
t255 = t157 * t165;
t137 = qJD(4) + t255;
t87 = t94 * t249;
t47 = t87 - t182;
t40 = t137 * pkin(4) + t47;
t278 = (qJD(5) * t40 + t14) * t170;
t155 = pkin(1) * t173 + pkin(2);
t196 = -pkin(1) * t248 + t172 * t155;
t118 = pkin(3) + t196;
t107 = t118 * t250;
t243 = pkin(1) * t247 + t168 * t155;
t112 = t158 + t243;
t70 = t171 * t112 + t107;
t127 = pkin(9) * t251 + t151;
t277 = qJD(4) * t127;
t236 = qJD(4) + qJD(5);
t229 = t94 * t250;
t48 = t197 * t171 + t229;
t275 = pkin(4) * t171;
t254 = t157 * t171;
t75 = (-pkin(4) * t254 - t94) * t164;
t99 = t157 * t114;
t273 = t75 * t99;
t252 = t164 * t167;
t226 = t166 * t252;
t202 = t157 * t226;
t97 = -t157 * t225 + t202;
t272 = t99 * t97;
t135 = t154 * t249;
t159 = t165 * pkin(4);
t86 = t217 * t167 + t135 + t159;
t149 = pkin(10) * t251;
t93 = t149 + t104;
t51 = -t166 * t93 + t170 * t86;
t271 = t51 * qJD(5) + t283 * t166 + t284 * t170;
t52 = t166 * t86 + t170 * t93;
t270 = -t52 * qJD(5) - t284 * t166 + t283 * t170;
t68 = t77 * t249;
t213 = -t167 * t76 + t68;
t54 = t171 * t92 + t229;
t23 = -qJD(4) * t54 + t213;
t160 = t164 ^ 2;
t253 = t160 * t171;
t269 = t23 * t165 + t77 * t253;
t152 = pkin(3) * t249;
t100 = t152 + t159 + t204;
t111 = t149 + t127;
t64 = t100 * t166 + t111 * t170;
t268 = t64 * qJD(5) + t282 * t166 + t281 * t170;
t63 = t100 * t170 - t111 * t166;
t267 = -t63 * qJD(5) + t281 * t166 - t282 * t170;
t266 = (-qJD(4) * t142 - t206) * t167 + t285;
t265 = -t279 + t286;
t263 = t166 * t48;
t262 = t167 * t23;
t261 = t170 * t48;
t239 = qJD(4) * t167;
t22 = -t92 * t239 + t234;
t260 = t22 * t165;
t222 = t164 * t239;
t259 = -pkin(9) * t222 + t287;
t258 = -t59 - t277;
t257 = t157 * t160;
t246 = t280 * t157;
t242 = t167 ^ 2 - t171 ^ 2;
t237 = qJD(4) - t137;
t113 = -t225 + t226;
t17 = t166 * t40 + t261;
t15 = -qJD(4) * t48 + t213;
t216 = -t166 * t14 + t170 * t15;
t5 = -t17 * qJD(5) + t216;
t223 = t157 * t239;
t55 = (pkin(4) * t223 - t77) * t164;
t83 = t236 * t114;
t235 = t55 * t113 + t5 * t165 + t75 * t83;
t88 = t155 * t240 + (t184 * qJD(2) - t169 * t241) * pkin(1);
t89 = -t155 * t241 + t175;
t233 = t118 * t220 + t171 * t88 + t89 * t250;
t227 = t157 * t252;
t16 = t170 * t40 - t263;
t215 = -qJD(5) * t263 + t166 * t15;
t4 = t215 + t278;
t181 = t236 * t226;
t82 = t181 + t280;
t224 = -t4 * t113 - t5 * t114 + t16 * t82 - t17 * t83;
t133 = qJD(5) + t137;
t219 = -pkin(4) * t133 - t40;
t218 = -t112 - t274;
t214 = -t89 * t157 - t77;
t212 = -t167 * t88 + t89 * t249;
t210 = -t110 * t157 - t77;
t146 = pkin(4) * t222;
t209 = t195 * t164 + t146;
t208 = t137 + t255;
t207 = pkin(4) * t227;
t156 = t157 ^ 2;
t201 = t167 * t156 * t253;
t200 = t55 * t114 - t4 * t165 - t75 * t82;
t194 = (-qJD(2) + t161) * t264;
t193 = pkin(1) * qJD(2) * (-qJD(1) - t161);
t192 = qJD(4) * (-t94 - t276);
t191 = t218 * t167;
t108 = t118 * t249;
t58 = t108 + t159 + t191;
t65 = t149 + t70;
t30 = -t166 * t65 + t170 * t58;
t31 = t166 * t58 + t170 * t65;
t53 = -t167 * t92 + t87;
t190 = -t167 * t54 - t171 * t53;
t189 = (-pkin(2) * t157 - t136) * qJD(3);
t188 = qJD(4) * (-t118 * t157 - t94);
t187 = t223 * t253;
t183 = t75 * t97 - t215;
t180 = t164 * (pkin(4) * t239 - t110);
t177 = t190 * qJD(4) - t262;
t132 = (-pkin(3) - t275) * t164;
t126 = -pkin(9) * t252 + t152;
t123 = -0.2e1 * t187;
t122 = 0.2e1 * t187;
t119 = (-t154 - t275) * t164;
t103 = -t142 * t167 + t135;
t96 = -0.2e1 * t242 * qJD(4) * t257;
t95 = (-t118 - t275) * t164;
t91 = t208 * t221;
t90 = t208 * t222;
t72 = -t164 * t89 + t146;
t69 = -t112 * t167 + t108;
t57 = t157 * t83;
t56 = t157 * t181 + t246;
t46 = -t97 ^ 2 + t99 ^ 2;
t39 = -t186 * t236 * t256 + t99 * t133;
t38 = t97 * t133 - t236 * t202 - t246;
t35 = -t133 * t83 - t165 * t57;
t34 = -t133 * t82 - t165 * t56;
t29 = -qJD(4) * t70 + t212;
t28 = -t112 * t239 + t233;
t27 = (t218 * t171 - t107) * qJD(4) + t212;
t26 = qJD(4) * t191 + t233;
t25 = t113 * t57 + t83 * t97;
t24 = -t114 * t56 - t82 * t99;
t21 = t170 * t47 - t263;
t20 = -t166 * t47 - t261;
t18 = t22 * t251;
t8 = t113 * t56 - t114 * t57 + t82 * t97 - t83 * t99;
t7 = -t31 * qJD(5) - t166 * t26 + t170 * t27;
t6 = t30 * qJD(5) + t166 * t27 + t170 * t26;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 * t193, t173 * t193, 0, 0, 0, 0, 0, 0, 0, 0, -t214, -t88 * t157 - t76, 0, t109 * t89 + t110 * t88 + t77 * t196 + t76 * t243, t122, t96, t91, t123, -t90, 0, t29 * t137 + (t167 * t188 + t89 * t254) * t160 + t269, -t28 * t137 - t260 + (t214 * t167 + t171 * t188) * t160, t18 + (-t262 + (-t167 * t29 + t171 * t28) * t157 + ((-t167 * t70 - t171 * t69) * t157 + t190) * qJD(4)) * t164, t22 * t70 + t23 * t69 + t54 * t28 + t53 * t29 + (t118 * t77 + t89 * t94) * t160, t24, t8, t34, t25, t35, 0, t133 * t7 + t57 * t95 + t72 * t97 + t235, -t133 * t6 - t56 * t95 + t72 * t99 + t200, t30 * t56 - t31 * t57 - t6 * t97 - t7 * t99 + t224, t16 * t7 + t17 * t6 + t30 * t5 + t31 * t4 + t55 * t95 + t72 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 * t194, t173 * t194, 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t157 + t168 * t189 + t174, t121 * t157 + (t189 - t203) * t172 + t244, 0, -t109 * t120 - t110 * t121 + (t168 * t76 + t172 * t77 + (-t109 * t168 + t110 * t172) * qJD(3)) * pkin(2), t122, t96, t91, t123, -t90, 0, t265 * t137 + (-t94 * t239 + (-t154 * t239 - t195 * t171) * t157) * t160 + t269, -t260 - t266 * t137 + (-t94 * t238 - t167 * t77 + (-t154 * t238 + t195 * t167) * t157) * t160, t18 + (((-qJD(4) * t103 + t266) * t171 + (-t265 - t279) * t167) * t157 + t177) * t164, t23 * t103 + t22 * t104 + t266 * t54 + t265 * t53 + (t154 * t77 - t195 * t94) * t160, t24, t8, t34, t25, t35, 0, t119 * t57 + t270 * t133 + t209 * t97 + t235, -t119 * t56 - t271 * t133 + t209 * t99 + t200, -t270 * t99 - t271 * t97 + t51 * t56 - t52 * t57 + t224, t55 * t119 + t270 * t16 + t271 * t17 + t209 * t75 + t4 * t52 + t5 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, t109 * t157 - t76, 0, 0, t122, t96, t91, t123, -t90, 0, t258 * t137 + (t110 * t254 + t167 * t192) * t160 + t269, -t260 - t259 * t137 + (t210 * t167 + t171 * t192) * t160, t18 + (((-qJD(4) * t126 + t259) * t171 + (-t258 - t277) * t167) * t157 + t177) * t164, t23 * t126 + t22 * t127 + t259 * t54 + t258 * t53 + (pkin(3) * t77 + t110 * t94) * t160, t24, t8, t34, t25, t35, 0, t132 * t57 - t268 * t133 + t97 * t180 + t235, -t132 * t56 + t267 * t133 + t99 * t180 + t200, t267 * t97 + t268 * t99 + t63 * t56 - t64 * t57 + t224, t55 * t132 - t268 * t16 - t267 * t17 + t75 * t180 + t4 * t64 + t5 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, t242 * t160 * t156, t237 * t157 * t251, t201, -t237 * t227, 0, -t92 * t238 + t54 * t137 + t68 + (-t76 + (-qJD(4) * t165 + t257) * t94) * t167, t157 * t94 * t253 + t137 * t53 - t22, 0, 0, t272, t46, t38, -t272, t39, 0, -t97 * t207 - t20 * t133 - t273 + (t219 * t166 - t261) * qJD(5) + t216, -t99 * t207 + t21 * t133 + (t219 * qJD(5) - t14) * t170 + t183, (t17 + t20) * t99 + (-t16 + t21) * t97 + (-t166 * t57 + t170 * t56 + (t166 * t99 - t170 * t97) * qJD(5)) * pkin(4), -t16 * t20 - t17 * t21 + (-t75 * t227 + t166 * t4 + t170 * t5 + (-t16 * t166 + t17 * t170) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t46, t38, -t272, t39, 0, t17 * t133 - t273 + t5, t16 * t133 + t183 - t278, 0, 0;];
tauc_reg = t1;
