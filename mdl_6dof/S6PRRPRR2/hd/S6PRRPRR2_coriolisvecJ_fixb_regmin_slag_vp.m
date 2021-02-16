% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:37
% EndTime: 2021-01-16 03:46:53
% DurationCPUTime: 3.97s
% Computational Cost: add. (4062->377), mult. (10631->545), div. (0->0), fcn. (8480->12), ass. (0->210)
t165 = cos(qJ(3));
t254 = cos(pkin(12));
t204 = t254 * t165;
t144 = qJD(2) * t204;
t156 = sin(pkin(12));
t161 = sin(qJ(3));
t233 = qJD(2) * t161;
t123 = t156 * t233 - t144;
t275 = qJD(5) + qJD(6);
t288 = t123 + t275;
t162 = sin(qJ(2));
t157 = sin(pkin(6));
t235 = qJD(1) * t157;
t215 = t162 * t235;
t231 = qJD(3) * t161;
t287 = pkin(3) * t231 - t215;
t179 = -t156 * t161 + t204;
t166 = cos(qJ(2));
t214 = t166 * t235;
t102 = t179 * t214;
t268 = qJ(4) + pkin(8);
t206 = qJD(3) * t268;
t121 = t165 * qJD(4) - t161 * t206;
t177 = -t161 * qJD(4) - t165 * t206;
t70 = t254 * t121 + t156 * t177;
t257 = t102 - t70;
t133 = t156 * t165 + t254 * t161;
t125 = t133 * qJD(3);
t128 = t179 * qJD(3);
t286 = -t125 * pkin(4) + t128 * pkin(9) - t287;
t126 = t133 * qJD(2);
t160 = sin(qJ(5));
t164 = cos(qJ(5));
t226 = t164 * qJD(3);
t105 = t160 * t126 - t226;
t137 = qJD(2) * pkin(8) + t215;
t202 = qJ(4) * qJD(2) + t137;
t158 = cos(pkin(6));
t234 = qJD(1) * t158;
t213 = t161 * t234;
t98 = t202 * t165 + t213;
t86 = t156 * t98;
t145 = t165 * t234;
t97 = -t202 * t161 + t145;
t91 = qJD(3) * pkin(3) + t97;
t34 = t254 * t91 - t86;
t31 = -qJD(3) * pkin(4) - t34;
t27 = t105 * pkin(5) + t31;
t107 = t160 * qJD(3) + t164 * t126;
t159 = sin(qJ(6));
t163 = cos(qJ(6));
t45 = t163 * t105 + t159 * t107;
t285 = t27 * t45;
t136 = t159 * t164 + t163 * t160;
t264 = t288 * t136;
t189 = t159 * t105 - t163 * t107;
t284 = t189 * t45;
t230 = qJD(5) * t160;
t250 = t123 * t160;
t283 = t230 + t250;
t282 = t189 ^ 2 - t45 ^ 2;
t118 = qJD(5) + t123;
t114 = qJD(6) + t118;
t227 = qJD(6) * t163;
t228 = qJD(6) * t159;
t225 = qJD(2) * qJD(3);
t208 = t161 * t225;
t142 = t156 * t208;
t116 = qJD(3) * t144 - t142;
t58 = qJD(5) * t226 + t164 * t116 - t126 * t230;
t59 = qJD(5) * t107 + t160 * t116;
t9 = -t105 * t227 - t107 * t228 - t159 * t59 + t163 * t58;
t281 = t45 * t114 + t9;
t115 = qJD(2) * t125;
t209 = t254 * t98;
t35 = t156 * t91 + t209;
t32 = qJD(3) * pkin(9) + t35;
t152 = -t165 * pkin(3) - pkin(2);
t117 = t152 * qJD(2) + qJD(4) - t214;
t51 = t123 * pkin(4) - t126 * pkin(9) + t117;
t20 = t160 * t51 + t164 * t32;
t194 = qJD(4) + t214;
t57 = (-t161 * t137 + t145) * qJD(3) + (-qJ(4) * t231 + t194 * t165) * qJD(2);
t210 = t254 * t57;
t184 = -t165 * t137 - t213;
t218 = qJD(3) * t165 * qJ(4);
t276 = t184 * qJD(3) + (-t194 * t161 - t218) * qJD(2);
t25 = t276 * t156 + t210;
t122 = pkin(3) * t208 + qJD(2) * t215;
t43 = t115 * pkin(4) - t116 * pkin(9) + t122;
t41 = t164 * t43;
t172 = -qJD(5) * t20 - t160 * t25 + t41;
t2 = t115 * pkin(5) - t58 * pkin(10) + t172;
t229 = qJD(5) * t164;
t185 = t160 * t43 + t164 * t25 + t51 * t229 - t32 * t230;
t3 = -t59 * pkin(10) + t185;
t221 = -t159 * t3 + t163 * t2;
t14 = -t105 * pkin(10) + t20;
t261 = t163 * t14;
t19 = -t160 * t32 + t164 * t51;
t13 = -t107 * pkin(10) + t19;
t8 = t118 * pkin(5) + t13;
t5 = t159 * t8 + t261;
t280 = -t5 * qJD(6) + t27 * t189 + t221;
t171 = qJD(6) * t189 - t159 * t58 - t163 * t59;
t279 = -t114 * t189 + t171;
t278 = t160 * t102 - t286 * t164;
t140 = t268 * t161;
t141 = t268 * t165;
t104 = -t156 * t140 + t254 * t141;
t84 = -pkin(4) * t179 - t133 * pkin(9) + t152;
t277 = t104 * t230 + t286 * t160 + t257 * t164 - t84 * t229;
t255 = t156 * t121 - t133 * t214 - t254 * t177;
t78 = t136 * t133;
t239 = t164 * t128;
t182 = -t133 * t230 + t239;
t135 = t159 * t160 - t163 * t164;
t265 = t288 * t135;
t274 = t114 * t265 - t136 * t115;
t11 = t14 * t228;
t212 = qJD(6) * t8 + t3;
t273 = t159 * t2 + t163 * t212 - t11;
t92 = t164 * t104;
t272 = -pkin(10) * t239 + t125 * pkin(5) - t160 * t70 + (-t92 + (pkin(10) * t133 - t84) * t160) * qJD(5) + t278;
t216 = t133 * t229;
t240 = t160 * t128;
t183 = t216 + t240;
t271 = pkin(10) * t183 + t277;
t270 = pkin(3) * t161;
t148 = t156 * pkin(3) + pkin(9);
t269 = pkin(10) + t148;
t39 = t254 * t97 - t86;
t223 = pkin(3) * t233;
t71 = t126 * pkin(4) + t123 * pkin(9) + t223;
t267 = t160 * t71 + t164 * t39;
t266 = t160 * t84 + t92;
t263 = qJD(2) * pkin(2);
t262 = t126 * t45;
t24 = t156 * t57 - t254 * t276;
t260 = t24 * t164;
t259 = t189 * t126;
t258 = t58 * t160;
t256 = pkin(5) * t183 + t255;
t253 = t105 * t118;
t252 = t107 * t118;
t251 = t107 * t126;
t249 = t126 * t105;
t248 = t133 * t160;
t247 = t133 * t164;
t245 = t157 * t162;
t244 = t157 * t166;
t168 = qJD(2) ^ 2;
t243 = t157 * t168;
t242 = t160 * t115;
t108 = t164 * t115;
t167 = qJD(3) ^ 2;
t238 = t167 * t161;
t237 = t167 * t165;
t236 = t161 ^ 2 - t165 ^ 2;
t232 = qJD(2) * t162;
t222 = t162 * t243;
t220 = t157 * t232;
t219 = qJD(2) * t244;
t37 = t156 * t97 + t209;
t205 = qJD(5) * t269;
t103 = t254 * t140 + t156 * t141;
t203 = t118 * t164;
t201 = t161 * t219;
t200 = t165 * t219;
t199 = -t264 * t114 - t135 * t115;
t198 = t283 * pkin(5) - t37;
t149 = -t254 * pkin(3) - pkin(4);
t130 = t269 * t160;
t197 = pkin(10) * t250 + qJD(6) * t130 + t160 * t205 + t267;
t131 = t269 * t164;
t64 = t164 * t71;
t196 = t126 * pkin(5) + qJD(6) * t131 - t160 * t39 + t64 + (pkin(10) * t123 + t205) * t164;
t81 = t164 * t84;
t26 = -pkin(5) * t179 - pkin(10) * t247 - t160 * t104 + t81;
t28 = -pkin(10) * t248 + t266;
t193 = t159 * t26 + t163 * t28;
t129 = t158 * t161 + t165 * t245;
t186 = t158 * t165 - t161 * t245;
t77 = t254 * t129 + t156 * t186;
t187 = t160 * t244 - t164 * t77;
t55 = -t160 * t77 - t164 * t244;
t192 = t159 * t187 + t163 * t55;
t191 = t159 * t55 - t163 * t187;
t190 = -t104 * t115 + t24 * t133;
t188 = -t283 * t118 + t108;
t180 = t263 * qJD(2);
t178 = -t148 * t115 + t118 * t31;
t176 = t129 * qJD(3);
t174 = -0.2e1 * qJD(3) * t263;
t169 = -t176 - t201;
t139 = -t164 * pkin(5) + t149;
t96 = qJD(3) * t186 + t200;
t85 = t115 * t179;
t79 = t135 * t133;
t76 = t156 * t129 - t254 * t186;
t68 = pkin(5) * t248 + t103;
t38 = t156 * t169 + t254 * t96;
t36 = t156 * t96 - t254 * t169;
t22 = -t228 * t248 + (t275 * t247 + t240) * t163 + t182 * t159;
t21 = -t135 * t128 - t275 * t78;
t17 = qJD(5) * t187 - t160 * t38 + t164 * t220;
t16 = qJD(5) * t55 + t160 * t220 + t164 * t38;
t12 = t59 * pkin(5) + t24;
t4 = -t159 * t14 + t163 * t8;
t1 = [0, 0, -t222, -t166 * t243, 0, 0, 0, 0, 0, -t165 * t222 + (-t176 - 0.2e1 * t201) * qJD(3), t161 * t222 + (-t96 - t200) * qJD(3), -t36 * qJD(3) + (-t115 * t166 + t123 * t232) * t157, -t38 * qJD(3) + (-t116 * t166 + t126 * t232) * t157, -t77 * t115 + t76 * t116 - t38 * t123 + t36 * t126, t24 * t76 + t25 * t77 - t34 * t36 + t35 * t38 + (t117 * t232 - t122 * t166) * t157, 0, 0, 0, 0, 0, t36 * t105 + t55 * t115 + t17 * t118 + t76 * t59, t36 * t107 + t115 * t187 - t16 * t118 + t76 * t58, 0, 0, 0, 0, 0, (-qJD(6) * t191 - t159 * t16 + t163 * t17) * t114 + t192 * t115 + t36 * t45 - t76 * t171, -(qJD(6) * t192 + t159 * t17 + t163 * t16) * t114 - t191 * t115 - t36 * t189 + t76 * t9; 0, 0, 0, 0, 0.2e1 * t165 * t208, -0.2e1 * t236 * t225, t237, -t238, 0, -pkin(8) * t237 + t161 * t174, pkin(8) * t238 + t165 * t174, -t123 * t215 + t152 * t115 + t117 * t125 - t122 * t179 + (t123 * t270 - t255) * qJD(3), -t126 * t215 + t152 * t116 + t117 * t128 + t122 * t133 + (t126 * t270 + t257) * qJD(3), t103 * t116 + t257 * t123 - t35 * t125 + t255 * t126 - t34 * t128 + t179 * t25 + t190, t24 * t103 + t25 * t104 + t287 * t117 + t122 * t152 - t255 * t34 - t257 * t35, t107 * t182 + t58 * t247, (-t105 * t164 - t107 * t160) * t128 + (-t258 - t164 * t59 + (t105 * t160 - t107 * t164) * qJD(5)) * t133, t107 * t125 + t133 * t108 + t118 * t182 - t179 * t58, -t105 * t125 - t118 * t183 - t133 * t242 + t179 * t59, t118 * t125 - t85, t81 * t115 - (-t229 * t32 + t41) * t179 + t19 * t125 + t103 * t59 + t31 * t216 + (-t104 * t229 + t278) * t118 + t255 * t105 + ((-qJD(5) * t84 - t70) * t118 - (-qJD(5) * t51 - t25) * t179 + t31 * t128 + t190) * t160, -t266 * t115 + t185 * t179 - t20 * t125 + t103 * t58 + t31 * t239 + (-t230 * t31 + t260) * t133 + t277 * t118 + t255 * t107, -t189 * t21 - t9 * t79, -t171 * t79 + t189 * t22 - t21 * t45 - t9 * t78, t21 * t114 - t79 * t115 - t125 * t189 - t179 * t9, -t22 * t114 - t78 * t115 - t45 * t125 - t171 * t179, t114 * t125 - t85, (-t159 * t28 + t163 * t26) * t115 - t221 * t179 + t4 * t125 - t68 * t171 + t12 * t78 + t27 * t22 + t256 * t45 + (t271 * t159 + t272 * t163) * t114 + (-t114 * t193 + t179 * t5) * qJD(6), -t193 * t115 + t273 * t179 - t5 * t125 + t68 * t9 - t12 * t79 + t27 * t21 - t256 * t189 + ((-qJD(6) * t26 + t271) * t163 + (qJD(6) * t28 - t272) * t159) * t114; 0, 0, 0, 0, -t161 * t168 * t165, t236 * t168, 0, 0, 0, t161 * t180, t165 * t180, t37 * qJD(3) - t117 * t126 - t123 * t223 - t24, -t210 + t117 * t123 + (-t156 * t184 + t39) * qJD(3) + (t156 * t218 + (-pkin(3) * t126 + t156 * t194) * t161) * qJD(2), (t35 - t37) * t126 + (-t34 + t39) * t123 + (-t115 * t156 - t254 * t116) * pkin(3), t34 * t37 - t35 * t39 + (-t117 * t233 + t156 * t25 - t254 * t24) * pkin(3), t107 * t203 + t258, (t58 - t253) * t164 + (-t59 - t252) * t160, t118 * t203 + t242 - t251, t188 + t249, -t118 * t126, -t37 * t105 - t19 * t126 + t149 * t59 - t260 + (-t148 * t229 - t64) * t118 + (t39 * t118 + t178) * t160, -t37 * t107 + t20 * t126 + t149 * t58 + t24 * t160 + (t148 * t230 + t267) * t118 + t178 * t164, t9 * t136 + t189 * t265, -t9 * t135 + t136 * t171 + t189 * t264 + t265 * t45, t259 - t274, t199 + t262, -t114 * t126, (-t163 * t130 - t159 * t131) * t115 - t139 * t171 + t12 * t135 - t4 * t126 + t198 * t45 + t264 * t27 + (t159 * t197 - t163 * t196) * t114, -(-t159 * t130 + t163 * t131) * t115 + t139 * t9 + t12 * t136 + t5 * t126 - t198 * t189 - t265 * t27 + (t159 * t196 + t163 * t197) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t126 * qJD(3), -t142 + (t144 - t123) * qJD(3), -t123 ^ 2 - t126 ^ 2, t35 * t123 + t34 * t126 + t122, 0, 0, 0, 0, 0, t188 - t249, -t118 ^ 2 * t164 - t242 - t251, 0, 0, 0, 0, 0, t199 - t262, t259 + t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t105, -t105 ^ 2 + t107 ^ 2, t58 + t253, t252 - t59, t115, -t31 * t107 + t20 * t118 + t172, t31 * t105 + t19 * t118 - t185, -t284, t282, t281, t279, t115, -(-t159 * t13 - t261) * t114 + (-t107 * t45 - t114 * t228 + t163 * t115) * pkin(5) + t280, t285 + t11 + (-t14 * t114 - t2) * t159 + (t13 * t114 - t212) * t163 + (t107 * t189 - t114 * t227 - t159 * t115) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t284, t282, t281, t279, t115, t5 * t114 + t280, t4 * t114 - t273 + t285;];
tauc_reg = t1;
