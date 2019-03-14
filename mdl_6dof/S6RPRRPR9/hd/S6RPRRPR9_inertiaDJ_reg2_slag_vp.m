% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_inertiaDJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:35
% EndTime: 2019-03-09 05:31:59
% DurationCPUTime: 9.98s
% Computational Cost: add. (19773->470), mult. (57106->870), div. (0->0), fcn. (61678->14), ass. (0->210)
t115 = sin(pkin(12));
t117 = sin(pkin(6));
t118 = cos(pkin(12));
t258 = cos(pkin(7));
t272 = sin(qJ(3));
t205 = t258 * t272;
t116 = sin(pkin(7));
t259 = cos(pkin(6));
t222 = t116 * t259;
t273 = cos(qJ(3));
t277 = (t115 * t273 + t118 * t205) * t117 + t272 * t222;
t286 = t277 * qJD(3);
t285 = 0.2e1 * t286;
t221 = t117 * t258;
t281 = t118 * t221 + t222;
t156 = t281 * t273;
t252 = t115 * t117;
t213 = t272 * t252;
t80 = t213 - t156;
t56 = t80 * t285;
t231 = pkin(1) * t259;
t169 = pkin(2) * t259 + t118 * t231;
t153 = (-pkin(9) * t258 - qJ(2)) * t252 + t169;
t232 = -pkin(2) * t118 - pkin(1);
t253 = t115 * t116;
t171 = (-pkin(9) * t253 + t232) * t117;
t289 = t116 * t171 + t153 * t258;
t122 = cos(qJ(4));
t120 = sin(qJ(4));
t246 = qJD(4) * t120;
t288 = t122 * t286 - t246 * t80;
t98 = qJD(3) * t213;
t152 = qJD(3) * t156 - t98;
t251 = t117 * t118;
t177 = t116 * t251 - t258 * t259;
t287 = -qJD(4) * t177 + t152;
t148 = pkin(9) * t281 + qJ(2) * t251 + t115 * t231;
t49 = t273 * t148 + t272 * t289;
t284 = qJD(3) * t49;
t114 = sin(pkin(13));
t143 = pkin(4) * t286;
t134 = -pkin(10) * t177 + t49;
t133 = qJD(4) * t134;
t136 = t80 * pkin(3) - pkin(10) * t277 - t116 * t153 + t171 * t258;
t135 = qJD(4) * t136;
t126 = -t120 * t135 - t122 * t133;
t184 = qJD(2) * t115 * t221;
t247 = qJD(2) * t117;
t228 = t118 * t247;
t154 = -t184 * t272 + t228 * t273;
t278 = t272 * t148 - t273 * t289;
t132 = qJD(3) * t278 - t154;
t211 = t247 * t253;
t137 = pkin(3) * t286 - pkin(10) * t152 + t211;
t16 = t120 * t132 + t122 * t137 + t126;
t146 = qJD(4) * t277;
t51 = -t120 * t146 + t122 * t287;
t61 = -t120 * t177 + t122 * t277;
t279 = -t51 * qJ(5) - t61 * qJD(5);
t125 = t143 + t16 + t279;
t15 = (t132 - t135) * t122 + (t133 - t137) * t120;
t165 = t120 * t277 + t122 * t177;
t239 = t120 * t287 + t122 * t146;
t129 = qJ(5) * t239 + qJD(5) * t165 + t15;
t257 = cos(pkin(13));
t128 = t257 * t129;
t124 = t114 * t125 - t128;
t45 = pkin(3) * t177 + t278;
t31 = pkin(4) * t165 + t45;
t46 = t114 * t61 + t165 * t257;
t47 = -t114 * t165 + t257 * t61;
t131 = t46 * pkin(5) - t47 * pkin(11) + t31;
t283 = -pkin(11) * t286 - qJD(6) * t131 - t124;
t245 = qJD(4) * t122;
t282 = t120 * t286 + t245 * t80;
t229 = t116 * t272;
t280 = t120 * t229 - t122 * t258;
t119 = sin(qJ(6));
t112 = t119 ^ 2;
t121 = cos(qJ(6));
t113 = t121 ^ 2;
t248 = t112 - t113;
t216 = qJD(6) * t248;
t230 = t116 * t273;
t92 = t120 * t258 + t122 * t229;
t69 = -t114 * t280 + t257 * t92;
t175 = t119 * t69 + t121 * t230;
t176 = t119 * t230 - t121 * t69;
t210 = qJD(3) * t230;
t81 = qJD(4) * t280 - t122 * t210;
t82 = -qJD(4) * t92 - t120 * t210;
t180 = -t114 * t82 + t257 * t81;
t225 = qJD(3) * t272;
t209 = t116 * t225;
t38 = qJD(6) * t175 - t119 * t209 + t121 * t180;
t39 = qJD(6) * t176 + t119 * t180 + t121 * t209;
t276 = (-t119 * t175 + t121 * t176) * qJD(6) + t119 * t38 - t121 * t39;
t267 = -qJ(5) - pkin(10);
t223 = qJD(4) * t267;
t173 = t120 * qJD(5) - t122 * t223;
t88 = t122 * qJD(5) + t120 * t223;
t155 = -t114 * t173 + t257 * t88;
t109 = -pkin(4) * t122 - pkin(3);
t218 = t257 * t122;
t254 = t114 * t120;
t94 = -t218 + t254;
t219 = t257 * t120;
t95 = t114 * t122 + t219;
t174 = -pkin(5) * t94 + pkin(11) * t95 - t109;
t167 = t121 * t174;
t237 = pkin(4) * t246;
t89 = t95 * qJD(4);
t90 = qJD(4) * t218 - t114 * t246;
t172 = pkin(5) * t89 - pkin(11) * t90 + t237;
t244 = qJD(6) * t119;
t102 = t267 * t122;
t84 = -t102 * t257 + t254 * t267;
t33 = qJD(6) * t167 - t119 * t172 - t121 * t155 + t244 * t84;
t54 = -t119 * t174 + t121 * t84;
t34 = -qJD(6) * t54 - t119 * t155 + t121 * t172;
t53 = -t119 * t84 - t167;
t275 = t119 * t33 - t121 * t34 + (t119 * t53 - t121 * t54) * qJD(6);
t26 = -t120 * t134 + t122 * t136;
t22 = t80 * pkin(4) - t61 * qJ(5) + t26;
t27 = t120 * t136 + t122 * t134;
t25 = -qJ(5) * t165 + t27;
t13 = t114 * t22 + t25 * t257;
t10 = pkin(11) * t80 + t13;
t168 = -t114 * t239 + t257 * t51;
t40 = t273 * t184 + t272 * t228 + t284;
t28 = pkin(4) * t239 + t40;
t30 = t114 * t51 + t239 * t257;
t127 = t30 * pkin(5) - pkin(11) * t168 + t28;
t1 = t10 * t244 - t119 * t127 + t121 * t283;
t243 = qJD(6) * t121;
t2 = -t10 * t243 + t119 * t283 + t121 * t127;
t7 = -t10 * t119 + t121 * t131;
t8 = t121 * t10 + t119 * t131;
t274 = t1 * t119 - t121 * t2 + (t119 * t7 - t121 * t8) * qJD(6);
t271 = pkin(4) * t114;
t57 = -t114 * t81 - t257 * t82;
t68 = t114 * t92 + t257 * t280;
t270 = t68 * t57;
t66 = t114 * t88 + t173 * t257;
t83 = -t102 * t114 - t219 * t267;
t269 = t83 * t66;
t268 = t95 * t90;
t36 = t119 * t80 + t121 * t47;
t20 = qJD(6) * t36 + t119 * t168 - t121 * t286;
t35 = t119 * t47 - t121 * t80;
t266 = -t119 * t20 - t243 * t35;
t19 = -t119 * t286 - t121 * t168 - t243 * t80 + t244 * t47;
t265 = t119 * t19;
t264 = t119 * t35;
t263 = t121 * t20;
t262 = t121 * t36;
t260 = t51 * t120;
t107 = pkin(11) + t271;
t256 = t107 * t119;
t255 = t107 * t121;
t250 = t119 * t121;
t242 = 0.2e1 * t46 * t30;
t241 = 0.2e1 * t94 * t89;
t240 = -0.2e1 * pkin(3) * qJD(4);
t108 = -pkin(4) * t257 - pkin(5);
t238 = 0.2e1 * qJD(6) * t108;
t234 = t95 * t244;
t233 = t95 * t243;
t227 = t119 * t243;
t226 = t120 * t245;
t224 = qJD(4) * t273;
t217 = -0.4e1 * t95 * t250;
t93 = t95 ^ 2;
t215 = t93 * t227;
t214 = t239 * t122;
t208 = t120 * t224;
t5 = t114 * t129 + t125 * t257;
t4 = -pkin(5) * t286 - t5;
t12 = -t114 * t25 + t22 * t257;
t9 = -pkin(5) * t80 - t12;
t207 = t4 * t95 + t9 * t90;
t204 = -t1 * t121 - t119 * t2;
t202 = -t119 * t8 - t121 * t7;
t200 = t30 * t94 + t46 * t89;
t199 = t30 * t95 + t46 * t90;
t198 = t57 * t83 + t68 * t66;
t197 = t57 * t95 + t68 * t90;
t196 = t66 * t95 + t83 * t90;
t195 = t89 * t95 + t90 * t94;
t194 = -t107 * t89 + t108 * t90;
t193 = t107 * t94 - t108 * t95;
t190 = -t119 * t54 - t121 * t53;
t188 = t119 * t176 + t121 * t175;
t186 = -0.2e1 * t259 * t247;
t183 = -t121 * t19 - t244 * t36;
t182 = t119 * t30 + t243 * t46;
t181 = t119 * t89 + t243 * t94;
t179 = t116 ^ 2 * t273 * t225;
t178 = 0.2e1 * (t115 ^ 2 + t118 ^ 2) * t117 ^ 2 * qJD(2);
t163 = t165 * t120;
t160 = qJD(6) * t190 - t119 * t34 - t121 * t33;
t159 = qJD(6) * t188 - t119 * t39 - t121 * t38;
t158 = -t16 * t120 - t15 * t122 + (-t27 * t120 - t26 * t122) * qJD(4);
t157 = -t82 * t120 - t122 * t81 + (-t120 * t92 + t122 * t280) * qJD(4);
t73 = t121 * t89 - t244 * t94;
t60 = t216 * t95 - t250 * t90;
t58 = -t116 * t169 + (qJ(2) * t253 + t232 * t258) * t117;
t23 = t121 * t30 - t244 * t46;
t6 = -t128 + (-t120 * t154 + t122 * (t98 * pkin(10) + t211) + t126 + (t120 * t278 + t122 * (pkin(3) * t277 - pkin(10) * t156) + t277 * pkin(4)) * qJD(3) + t279) * t114;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t186, t118 * t186, t178, qJ(2) * t178, 0.2e1 * t277 * t152, 0.2e1 * t98 * t80 + 0.2e1 * (-t156 * t80 - t277 ^ 2) * qJD(3), -0.2e1 * t152 * t177, t56, t177 * t285, 0, 0.2e1 * t177 * t40 + 0.2e1 * t211 * t80 + 0.2e1 * t286 * t58, -0.2e1 * t132 * t177 + 0.2e1 * t152 * t58 + 0.2e1 * t211 * t277, 0.2e1 * t132 * t80 + 0.2e1 * t152 * t278 + 0.2e1 * t277 * t40 - 0.2e1 * t286 * t49, 0.2e1 * (t49 * (-t115 * t205 + t118 * t273) + t58 * t253) * t247 - 0.2e1 * (-t40 + t284) * t278, 0.2e1 * t61 * t51, -0.2e1 * t165 * t51 - 0.2e1 * t239 * t61, 0.2e1 * t286 * t61 + 0.2e1 * t51 * t80, 0.2e1 * t165 * t239, -0.2e1 * t165 * t286 - 0.2e1 * t239 * t80, t56, 0.2e1 * t16 * t80 + 0.2e1 * t165 * t40 + 0.2e1 * t239 * t45 + 0.2e1 * t26 * t286, 0.2e1 * t15 * t80 - 0.2e1 * t27 * t286 + 0.2e1 * t40 * t61 + 0.2e1 * t45 * t51, 0.2e1 * t15 * t165 - 0.2e1 * t16 * t61 - 0.2e1 * t239 * t27 - 0.2e1 * t26 * t51, -0.2e1 * t15 * t27 + 0.2e1 * t16 * t26 + 0.2e1 * t40 * t45, 0.2e1 * t47 * t168, -0.2e1 * t168 * t46 - 0.2e1 * t30 * t47, 0.2e1 * t168 * t80 + 0.2e1 * t286 * t47, t242, -0.2e1 * t286 * t46 - 0.2e1 * t30 * t80, t56, 0.2e1 * t12 * t286 + 0.2e1 * t28 * t46 + 0.2e1 * t31 * t30 + 0.2e1 * t5 * t80, -0.2e1 * t13 * t286 + 0.2e1 * t168 * t31 + 0.2e1 * t28 * t47 - 0.2e1 * t6 * t80, -0.2e1 * t12 * t168 - 0.2e1 * t13 * t30 - 0.2e1 * t46 * t6 - 0.2e1 * t47 * t5, 0.2e1 * t12 * t5 + 0.2e1 * t13 * t6 + 0.2e1 * t28 * t31, -0.2e1 * t36 * t19, 0.2e1 * t19 * t35 - 0.2e1 * t20 * t36, -0.2e1 * t19 * t46 + 0.2e1 * t30 * t36, 0.2e1 * t35 * t20, -0.2e1 * t20 * t46 - 0.2e1 * t30 * t35, t242, 0.2e1 * t2 * t46 + 0.2e1 * t20 * t9 + 0.2e1 * t30 * t7 + 0.2e1 * t35 * t4, 0.2e1 * t1 * t46 - 0.2e1 * t19 * t9 - 0.2e1 * t30 * t8 + 0.2e1 * t36 * t4, 0.2e1 * t1 * t35 + 0.2e1 * t19 * t7 - 0.2e1 * t2 * t36 - 0.2e1 * t20 * t8, -0.2e1 * t1 * t8 + 0.2e1 * t2 * t7 + 0.2e1 * t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t177 * t229 + t258 * t277) * qJD(3), t152 * t258 + t177 * t210, -t152 * t230 - t210 * t80, t116 * t184 - t132 * t229 + t209 * t278 + t210 * t49 - t230 * t40, 0, 0, 0, 0, 0, 0, t165 * t209 - t230 * t239 - t280 * t286 + t82 * t80, t209 * t61 - t230 * t51 - t286 * t92 + t80 * t81, t165 * t81 - t239 * t92 + t280 * t51 - t61 * t82, -t15 * t92 - t16 * t280 + t26 * t82 - t27 * t81 + (t225 * t45 - t273 * t40) * t116, 0, 0, 0, 0, 0, 0, t209 * t46 - t230 * t30 - t286 * t68 - t57 * t80, -t168 * t230 + t180 * t80 + t209 * t47 - t286 * t69, t168 * t68 + t180 * t46 - t69 * t30 + t57 * t47, t6 * t69 - t13 * t180 - t5 * t68 - t12 * t57 + (t225 * t31 - t273 * t28) * t116, 0, 0, 0, 0, 0, 0, -t175 * t30 + t20 * t68 + t35 * t57 + t39 * t46, t176 * t30 - t19 * t68 + t36 * t57 + t38 * t46, -t175 * t19 + t176 * t20 + t35 * t38 - t36 * t39, t1 * t176 - t175 * t2 - t38 * t8 + t39 * t7 + t4 * t68 + t57 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t280 * t82 - 0.2e1 * t81 * t92 - 0.2e1 * t179, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t180 * t69 - 0.2e1 * t179 + 0.2e1 * t270, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t175 * t39 + 0.2e1 * t176 * t38 + 0.2e1 * t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, -t286, 0, -t40, t132, 0, 0, t245 * t61 + t260, -t120 * t239 + t51 * t122 + (-t61 * t120 - t122 * t165) * qJD(4), t282, qJD(4) * t163 - t214, t288, 0, -pkin(3) * t239 - pkin(10) * t282 - t40 * t122 + t246 * t45, -pkin(3) * t51 - pkin(10) * t288 + t40 * t120 + t245 * t45 (-t214 + t260 + (t122 * t61 + t163) * qJD(4)) * pkin(10) + t158, -t40 * pkin(3) + pkin(10) * t158, t168 * t95 + t47 * t90, -t168 * t94 - t47 * t89 - t199, t286 * t95 + t90 * t80, t200, -t286 * t94 - t89 * t80, 0, t109 * t30 + t237 * t46 + t28 * t94 - t286 * t83 + t31 * t89 - t66 * t80, t109 * t168 - t155 * t80 + t237 * t47 + t28 * t95 - t286 * t84 + t31 * t90, -t12 * t90 - t13 * t89 - t155 * t46 + t168 * t83 - t84 * t30 + t66 * t47 - t5 * t95 - t6 * t94, t28 * t109 - t12 * t66 + t13 * t155 + t237 * t31 - t5 * t83 + t6 * t84, -t36 * t234 + (-t19 * t95 + t36 * t90) * t121 (-t119 * t36 - t121 * t35) * t90 + (t265 - t263 + (-t262 + t264) * qJD(6)) * t95, t121 * t199 - t19 * t94 - t234 * t46 + t36 * t89, t35 * t233 + (t20 * t95 + t35 * t90) * t119, -t119 * t199 - t20 * t94 - t233 * t46 - t35 * t89, t200, t119 * t207 + t2 * t94 + t20 * t83 + t233 * t9 + t30 * t53 + t34 * t46 + t35 * t66 + t7 * t89, t1 * t94 + t121 * t207 - t19 * t83 - t234 * t9 - t30 * t54 + t33 * t46 + t36 * t66 - t8 * t89, t19 * t53 - t20 * t54 + t202 * t90 + t274 * t95 + t33 * t35 - t34 * t36, -t1 * t54 + t2 * t53 - t33 * t8 + t34 * t7 + t4 * t83 + t66 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, -t210, 0, 0, 0, 0, 0, 0, 0, 0 (-t122 * t225 - t208) * t116 (t120 * t225 - t122 * t224) * t116, t157, -pkin(3) * t209 + pkin(10) * t157, 0, 0, 0, 0, 0, 0 (t225 * t94 - t273 * t89) * t116 (t225 * t95 - t273 * t90) * t116, t180 * t94 - t69 * t89 + t197, -t180 * t84 + t69 * t155 + (-pkin(4) * t208 + t109 * t225) * t116 + t198, 0, 0, 0, 0, 0, 0, t119 * t197 - t175 * t89 + t233 * t68 + t39 * t94, t121 * t197 + t176 * t89 - t234 * t68 + t38 * t94, t188 * t90 + t276 * t95, -t175 * t34 + t176 * t33 - t38 * t54 + t39 * t53 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t226, 0.2e1 * (-t120 ^ 2 + t122 ^ 2) * qJD(4), 0, -0.2e1 * t226, 0, 0, t120 * t240, t122 * t240, 0, 0, 0.2e1 * t268, -0.2e1 * t195, 0, t241, 0, 0, 0.2e1 * t109 * t89 + 0.2e1 * t237 * t94, 0.2e1 * t109 * t90 + 0.2e1 * t237 * t95, -0.2e1 * t155 * t94 - 0.2e1 * t84 * t89 + 0.2e1 * t196, 0.2e1 * t109 * t237 + 0.2e1 * t155 * t84 + 0.2e1 * t269, 0.2e1 * t113 * t268 - 0.2e1 * t215, 0.2e1 * t216 * t93 + t217 * t90, 0.2e1 * t121 * t195 - 0.2e1 * t234 * t94, 0.2e1 * t112 * t268 + 0.2e1 * t215, -0.2e1 * t119 * t195 - 0.2e1 * t233 * t94, t241, 0.2e1 * t119 * t196 + 0.2e1 * t233 * t83 + 0.2e1 * t34 * t94 + 0.2e1 * t53 * t89, 0.2e1 * t121 * t196 - 0.2e1 * t234 * t83 + 0.2e1 * t33 * t94 - 0.2e1 * t54 * t89, 0.2e1 * t190 * t90 + 0.2e1 * t275 * t95, -0.2e1 * t33 * t54 + 0.2e1 * t34 * t53 + 0.2e1 * t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t239, t286, t16, t15, 0, 0, 0, 0, t168, 0, -t30, t286, t143 * t257 + t5, -t271 * t286 - t124 (-t114 * t30 - t168 * t257) * pkin(4) (t114 * t6 + t257 * t5) * pkin(4), t243 * t36 - t265, t183 + t266, t182, t244 * t35 - t263, t23, 0, -t30 * t256 + t108 * t20 - t121 * t4 + (t119 * t9 - t255 * t46) * qJD(6), -t30 * t255 - t108 * t19 + t119 * t4 + (t121 * t9 + t256 * t46) * qJD(6) (-t263 - t265) * t107 + ((t262 + t264) * t107 + t202) * qJD(6) + t204, t108 * t4 + (qJD(6) * t202 + t204) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t81, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t180, 0 (-t114 * t180 - t257 * t57) * pkin(4), 0, 0, 0, 0, 0, 0, -t121 * t57 + t244 * t68, t119 * t57 + t243 * t68, t159, t107 * t159 + t108 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, 0, -t246, 0, -pkin(10) * t245, pkin(10) * t246, 0, 0, 0, 0, t90, 0, -t89, 0, -t66, -t155 (-t114 * t89 - t257 * t90) * pkin(4) (t114 * t155 - t257 * t66) * pkin(4), -t60, qJD(6) * t217 - t248 * t90, t181, t60, t73, 0, -t121 * t66 + t194 * t119 + (t119 * t83 - t121 * t193) * qJD(6), t119 * t66 + t194 * t121 + (t119 * t193 + t121 * t83) * qJD(6), t160, t107 * t160 + t108 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t227, -0.2e1 * t216, 0, -0.2e1 * t227, 0, 0, t119 * t238, t121 * t238, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t168, 0, t28, 0, 0, 0, 0, 0, 0, t23, -t182, -t183 + t266, -t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t90, 0, t237, 0, 0, 0, 0, 0, 0, t73, -t181 (-t112 - t113) * t90, -t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, -t20, t30, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t90 - t234, 0, -t119 * t90 - t233, t89, t34, t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, 0, -t244, 0, -t107 * t243, t107 * t244, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, -t243, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;