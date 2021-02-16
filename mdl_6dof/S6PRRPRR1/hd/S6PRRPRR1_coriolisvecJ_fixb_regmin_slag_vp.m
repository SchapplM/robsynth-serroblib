% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR1
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
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:23
% EndTime: 2021-01-16 03:28:37
% DurationCPUTime: 3.57s
% Computational Cost: add. (4419->325), mult. (11670->475), div. (0->0), fcn. (9508->12), ass. (0->196)
t167 = cos(qJ(6));
t217 = qJD(6) * t167;
t159 = sin(pkin(12));
t161 = cos(pkin(12));
t169 = cos(qJ(3));
t228 = t161 * t169;
t206 = qJD(2) * t228;
t165 = sin(qJ(3));
t221 = qJD(2) * t165;
t128 = -t159 * t221 + t206;
t168 = cos(qJ(5));
t116 = t168 * t128;
t136 = t159 * t169 + t161 * t165;
t130 = t136 * qJD(2);
t164 = sin(qJ(5));
t83 = -t164 * t130 + t116;
t278 = t167 * t83;
t285 = t217 - t278;
t248 = qJ(4) + pkin(8);
t199 = qJD(3) * t248;
t124 = t169 * qJD(4) - t165 * t199;
t125 = -t165 * qJD(4) - t169 * t199;
t170 = cos(qJ(2));
t160 = sin(pkin(6));
t223 = qJD(1) * t160;
t207 = t170 * t223;
t236 = t159 * t124 - t161 * t125 - t136 * t207;
t135 = t159 * t165 - t228;
t235 = t161 * t124 + t159 * t125 + t135 * t207;
t156 = qJD(3) + qJD(5);
t237 = t83 * t156;
t129 = t136 * qJD(3);
t119 = qJD(2) * t129;
t215 = qJD(2) * qJD(3);
t204 = t165 * t215;
t144 = t159 * t204;
t203 = t169 * t215;
t120 = t161 * t203 - t144;
t219 = qJD(5) * t164;
t43 = qJD(5) * t116 - t164 * t119 + t168 * t120 - t130 * t219;
t284 = t43 - t237;
t225 = -qJD(6) + t83;
t283 = qJD(6) + t225;
t163 = sin(qJ(6));
t185 = t164 * t128 + t168 * t130;
t218 = qJD(6) * t163;
t24 = t156 * t217 + t167 * t43 - t185 * t218;
t71 = t163 * t156 + t167 * t185;
t25 = qJD(6) * t71 + t163 * t43;
t69 = -t167 * t156 + t163 * t185;
t282 = -t163 * t25 + t24 * t167 - t285 * t69;
t22 = t24 * t163;
t281 = t285 * t71 + t22;
t44 = qJD(5) * t185 + t168 * t119 + t164 * t120;
t37 = t163 * t44;
t72 = t225 * t217;
t246 = t37 - t72;
t251 = t71 * t185;
t280 = t225 * t278 + t246 - t251;
t254 = t130 * pkin(9);
t166 = sin(qJ(2));
t208 = t166 * t223;
t140 = qJD(2) * pkin(8) + t208;
t194 = qJ(4) * qJD(2) + t140;
t162 = cos(pkin(6));
t222 = qJD(1) * t162;
t205 = t165 * t222;
t99 = t169 * t194 + t205;
t91 = t159 * t99;
t148 = t169 * t222;
t98 = -t165 * t194 + t148;
t95 = qJD(3) * pkin(3) + t98;
t52 = t161 * t95 - t91;
t40 = qJD(3) * pkin(4) - t254 + t52;
t255 = t128 * pkin(9);
t244 = t161 * t99;
t53 = t159 * t95 + t244;
t45 = t53 + t255;
t16 = -t164 * t45 + t168 * t40;
t14 = -t156 * pkin(5) - t16;
t279 = t14 * t83;
t277 = t185 * t83;
t132 = t135 * qJD(3);
t276 = -t132 * pkin(9) + t236;
t275 = -t129 * pkin(9) + t235;
t274 = t163 * t225;
t238 = t185 * t156;
t273 = -t44 + t238;
t216 = t165 * qJD(3);
t272 = pkin(3) * t216 - t208;
t270 = t185 ^ 2 - t83 ^ 2;
t49 = pkin(5) * t185 - t83 * pkin(10);
t187 = qJD(4) + t207;
t67 = (-t165 * t140 + t148) * qJD(3) + (-qJ(4) * t216 + t169 * t187) * qJD(2);
t68 = (-t169 * t140 - t205) * qJD(3) + (-qJD(3) * t169 * qJ(4) - t165 * t187) * qJD(2);
t30 = -t159 * t67 + t161 * t68;
t27 = -t120 * pkin(9) + t30;
t31 = t159 * t68 + t161 * t67;
t28 = -t119 * pkin(9) + t31;
t2 = (qJD(5) * t40 + t28) * t168 + t164 * t27 - t45 * t219;
t153 = -t169 * pkin(3) - pkin(2);
t121 = t153 * qJD(2) + qJD(4) - t207;
t86 = -t128 * pkin(4) + t121;
t269 = -t86 * t83 - t2;
t250 = t185 * t69;
t266 = t225 * t185;
t39 = t167 * t44;
t265 = -t218 * t225 - t39;
t189 = t129 * pkin(4) + t272;
t17 = t164 * t40 + t168 * t45;
t15 = t156 * pkin(10) + t17;
t29 = -pkin(5) * t83 - pkin(10) * t185 + t86;
t186 = t163 * t15 - t167 * t29;
t264 = t14 * t218 + t185 * t186;
t3 = qJD(5) * t17 + t164 * t28 - t168 * t27;
t5 = t167 * t15 + t163 * t29;
t263 = t14 * t217 + t3 * t163 + t5 * t185;
t142 = t248 * t165;
t143 = t248 * t169;
t102 = -t161 * t142 - t159 * t143;
t75 = -t136 * pkin(9) + t102;
t103 = -t159 * t142 + t161 * t143;
t76 = -t135 * pkin(9) + t103;
t35 = t164 * t76 - t168 * t75;
t259 = qJD(5) * t35 + t276 * t164 - t275 * t168;
t110 = t135 * pkin(4) + t153;
t89 = t168 * t135 + t164 * t136;
t90 = -t164 * t135 + t168 * t136;
t33 = t89 * pkin(5) - t90 * pkin(10) + t110;
t36 = t164 * t75 + t168 * t76;
t50 = -qJD(5) * t89 - t164 * t129 - t168 * t132;
t262 = -(qJD(6) * t29 + t2) * t89 + t14 * t50 + t3 * t90 - (-qJD(6) * t33 + t259) * t225 - t36 * t44;
t261 = -t185 * t86 - t3;
t258 = qJD(5) * t36 + t275 * t164 + t276 * t168;
t257 = pkin(3) * t159;
t256 = pkin(3) * t165;
t253 = t14 * t90;
t252 = t33 * t44;
t249 = t90 * t44;
t57 = t161 * t98 - t91;
t245 = qJD(2) * pkin(2);
t242 = t163 * t71;
t152 = t161 * pkin(3) + pkin(4);
t182 = t168 * t152 - t164 * t257;
t55 = -t159 * t98 - t244;
t46 = t55 - t255;
t47 = t57 - t254;
t234 = -t182 * qJD(5) + t164 * t46 + t168 * t47;
t183 = t164 * t152 + t168 * t257;
t233 = t183 * qJD(5) - t164 * t47 + t168 * t46;
t231 = t160 * t166;
t230 = t160 * t170;
t172 = qJD(2) ^ 2;
t229 = t160 * t172;
t171 = qJD(3) ^ 2;
t227 = t171 * t165;
t226 = t171 * t169;
t126 = pkin(3) * t204 + qJD(2) * t208;
t224 = t165 ^ 2 - t169 ^ 2;
t220 = qJD(2) * t166;
t154 = pkin(3) * t221;
t212 = t90 * t218;
t211 = t166 * t229;
t210 = t160 * t220;
t209 = qJD(2) * t230;
t104 = t130 * pkin(4) + t154;
t123 = pkin(10) + t183;
t196 = qJD(6) * t123 + t104 + t49;
t193 = t165 * t209;
t192 = t169 * t209;
t85 = t119 * pkin(4) + t126;
t51 = qJD(5) * t90 + t168 * t129 - t164 * t132;
t191 = t51 * pkin(5) - t50 * pkin(10) + t189;
t190 = -t225 * t50 + t249;
t133 = t162 * t169 - t165 * t231;
t134 = t162 * t165 + t169 * t231;
t78 = t161 * t133 - t159 * t134;
t79 = t159 * t133 + t161 * t134;
t41 = t164 * t79 - t168 * t78;
t42 = t164 * t78 + t168 * t79;
t184 = -t274 * t83 - t265;
t181 = -t163 * t42 - t167 * t230;
t180 = t163 * t230 - t167 * t42;
t179 = t245 * qJD(2);
t178 = -t123 * t44 - t225 * t234 - t279;
t177 = -0.2e1 * qJD(3) * t245;
t122 = -pkin(5) - t182;
t97 = -qJD(3) * t134 - t193;
t96 = qJD(3) * t133 + t192;
t56 = t159 * t97 + t161 * t96;
t54 = -t159 * t96 + t161 * t97;
t11 = t44 * pkin(5) - t43 * pkin(10) + t85;
t10 = t167 * t11;
t7 = qJD(5) * t42 + t164 * t56 - t168 * t54;
t6 = -qJD(5) * t41 + t164 * t54 + t168 * t56;
t1 = [0, 0, -t211, -t170 * t229, 0, 0, 0, 0, 0, -t169 * t211 + (t97 - t193) * qJD(3), t165 * t211 + (-t96 - t192) * qJD(3), t54 * qJD(3) + (-t119 * t170 - t128 * t220) * t160, -t56 * qJD(3) + (-t120 * t170 + t130 * t220) * t160, -t79 * t119 - t78 * t120 + t56 * t128 - t54 * t130, t30 * t78 + t31 * t79 + t52 * t54 + t53 * t56 + (t121 * t220 - t126 * t170) * t160, 0, 0, 0, 0, 0, -t7 * t156 + (-t170 * t44 - t220 * t83) * t160, -t6 * t156 + (-t170 * t43 + t185 * t220) * t160, 0, 0, 0, 0, 0, -(qJD(6) * t180 - t163 * t6 + t167 * t210) * t225 + t181 * t44 + t7 * t69 + t41 * t25, (qJD(6) * t181 + t163 * t210 + t167 * t6) * t225 + t180 * t44 + t7 * t71 + t41 * t24; 0, 0, 0, 0, 0.2e1 * t165 * t203, -0.2e1 * t224 * t215, t226, -t227, 0, -pkin(8) * t226 + t165 * t177, pkin(8) * t227 + t169 * t177, t128 * t208 + t153 * t119 + t121 * t129 + t126 * t135 + (-t128 * t256 - t236) * qJD(3), -t130 * t208 + t153 * t120 - t121 * t132 + t126 * t136 + (t130 * t256 - t235) * qJD(3), -t102 * t120 - t103 * t119 + t128 * t235 - t53 * t129 + t130 * t236 + t52 * t132 - t31 * t135 - t30 * t136, t30 * t102 + t31 * t103 + t272 * t121 + t126 * t153 + t235 * t53 - t236 * t52, t185 * t50 + t43 * t90, -t185 * t51 - t43 * t89 + t50 * t83 - t249, t50 * t156, -t51 * t156, 0, t110 * t44 - t258 * t156 - t189 * t83 + t86 * t51 + t85 * t89, t110 * t43 + t259 * t156 + t185 * t189 + t86 * t50 + t85 * t90, -t71 * t212 + (t24 * t90 + t50 * t71) * t167, (-t167 * t69 - t242) * t50 + (-t22 - t167 * t25 + (t163 * t69 - t167 * t71) * qJD(6)) * t90, t167 * t190 + t212 * t225 + t24 * t89 + t71 * t51, -t163 * t190 - t25 * t89 - t69 * t51 + t72 * t90, -t225 * t51 + t44 * t89, t10 * t89 + t35 * t25 - t186 * t51 + t258 * t69 + (t252 - t191 * t225 + (-t15 * t89 + t225 * t36 + t253) * qJD(6)) * t167 + t262 * t163, t35 * t24 - t5 * t51 + t258 * t71 + (-t252 - (-qJD(6) * t15 + t11) * t89 - qJD(6) * t253 - (qJD(6) * t36 - t191) * t225) * t163 + t262 * t167; 0, 0, 0, 0, -t165 * t172 * t169, t224 * t172, 0, 0, 0, t165 * t179, t169 * t179, -t55 * qJD(3) - t121 * t130 + t128 * t154 + t30, t57 * qJD(3) - t121 * t128 - t130 * t154 - t31, (t53 + t55) * t130 + (t52 - t57) * t128 + (-t119 * t159 - t120 * t161) * pkin(3), -t52 * t55 - t53 * t57 + (-t121 * t221 + t159 * t31 + t161 * t30) * pkin(3), -t277, t270, t284, t273, 0, t104 * t83 - t156 * t233 + t261, -t104 * t185 + t156 * t234 + t269, t281, t225 * t242 + t282, t280, t184 + t250, t266, t122 * t25 + t233 * t69 + (t196 * t225 - t3) * t167 + t178 * t163 + t264, t122 * t24 + t167 * t178 - t196 * t274 + t233 * t71 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t130 * qJD(3), -t144 + (t128 + t206) * qJD(3), -t128 ^ 2 - t130 ^ 2, -t53 * t128 + t52 * t130 + t126, 0, 0, 0, 0, 0, t44 + t238, t43 + t237, 0, 0, 0, 0, 0, t184 - t250, -t167 * t225 ^ 2 - t251 - t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t277, t270, t284, t273, 0, t17 * t156 + t261, t16 * t156 + t269, t281, t274 * t71 + t282, t280, -t225 * t274 + t250 + t39, t266, -pkin(5) * t25 - t3 * t167 + (-t163 * t16 + t167 * t49) * t225 - t17 * t69 - t163 * t279 - t246 * pkin(10) + t264, -pkin(5) * t24 - (t167 * t16 + t163 * t49) * t225 - t17 * t71 - t14 * t278 + t265 * pkin(10) + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t69, -t69 ^ 2 + t71 ^ 2, -t225 * t69 + t24, -t225 * t71 - t25, t44, -t14 * t71 - t163 * t2 - t283 * t5 + t10, -t163 * t11 + t14 * t69 - t167 * t2 + t186 * t283;];
tauc_reg = t1;
