% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:39
% EndTime: 2019-03-09 02:18:50
% DurationCPUTime: 4.20s
% Computational Cost: add. (9977->361), mult. (24078->476), div. (0->0), fcn. (18454->10), ass. (0->202)
t148 = cos(qJ(6));
t146 = sin(qJ(6));
t208 = qJD(6) * t148;
t247 = cos(qJ(4));
t197 = qJD(4) * t247;
t144 = cos(pkin(11));
t211 = qJD(1) * t144;
t131 = t197 * t211;
t142 = sin(pkin(11));
t245 = sin(qJ(4));
t200 = t245 * t142;
t184 = qJD(1) * t200;
t111 = qJD(4) * t184 - t131;
t167 = -t247 * t142 - t245 * t144;
t256 = t167 * qJD(1);
t112 = qJD(4) * t256;
t168 = -t247 * t144 + t200;
t121 = t168 * qJD(1);
t147 = sin(qJ(5));
t246 = cos(qJ(5));
t195 = qJD(5) * t246;
t210 = qJD(5) * t147;
t170 = -t246 * t111 + t147 * t112 - t121 * t195 + t210 * t256;
t171 = t147 * t121 + t246 * t256;
t207 = qJD(4) + qJD(5);
t79 = t146 * t207 - t148 * t171;
t216 = qJD(6) * t79;
t37 = t146 * t170 + t216;
t188 = t148 * t207;
t77 = -t146 * t171 - t188;
t230 = -t146 * t37 - t77 * t208;
t85 = -t246 * t121 + t147 * t256;
t270 = qJD(6) - t85;
t281 = t146 * t270;
t260 = t79 * t281;
t273 = t148 * t85;
t209 = qJD(6) * t146;
t36 = -qJD(6) * t188 - t148 * t170 - t171 * t209;
t283 = -t148 * t36 + t273 * t77 + t230 - t260;
t35 = t37 * t148;
t282 = t281 * t77 - t35;
t135 = sin(pkin(10)) * pkin(1) + qJ(3);
t129 = t135 * qJD(1);
t139 = t144 * qJD(2);
t97 = t139 + (-pkin(7) * qJD(1) - t129) * t142;
t105 = t142 * qJD(2) + t144 * t129;
t98 = pkin(7) * t211 + t105;
t169 = -t245 * t97 - t247 * t98;
t64 = -t121 * pkin(8) - t169;
t205 = t246 * t64;
t204 = t245 * t98;
t68 = t247 * t97 - t204;
t63 = pkin(8) * t256 + t68;
t61 = qJD(4) * pkin(4) + t63;
t29 = t147 * t61 + t205;
t27 = t207 * pkin(9) + t29;
t128 = -cos(pkin(10)) * pkin(1) - t144 * pkin(3) - pkin(2);
t118 = t128 * qJD(1) + qJD(3);
t89 = pkin(4) * t121 + t118;
t41 = -pkin(5) * t85 + pkin(9) * t171 + t89;
t174 = t146 * t27 - t148 * t41;
t280 = t270 * t174;
t259 = t85 * t207;
t279 = t170 - t259;
t233 = t85 ^ 2;
t234 = t171 ^ 2;
t278 = -t233 + t234;
t10 = t146 * t41 + t148 * t27;
t160 = t167 * qJD(3);
t156 = qJD(1) * t160;
t59 = t169 * qJD(4) + t156;
t152 = t111 * pkin(8) + t59;
t198 = qJD(3) * t245;
t199 = qJD(3) * t247;
t269 = -t142 * t198 + t144 * t199;
t166 = t269 * qJD(1) + t97 * t197;
t196 = qJD(4) * t245;
t58 = -t98 * t196 + t166;
t47 = t112 * pkin(8) + t58;
t150 = -t147 * t152 - t61 * t195 + t64 * t210 - t246 * t47;
t244 = pkin(4) * t112;
t190 = t147 * t111 + t246 * t112;
t267 = t171 * qJD(5);
t53 = -t190 - t267;
t20 = pkin(5) * t53 - pkin(9) * t170 - t244;
t3 = -qJD(6) * t10 + t146 * t150 + t148 * t20;
t262 = t10 * t270 + t3;
t73 = t79 * t208;
t276 = -t36 * t146 - t273 * t79 + t73;
t49 = t146 * t53;
t226 = t208 * t270 + t49;
t235 = t79 * t171;
t275 = -t270 * t273 + t226 + t235;
t220 = t147 * t64;
t28 = t246 * t61 - t220;
t26 = -t207 * pkin(5) - t28;
t274 = t26 * t85;
t237 = t77 * t171;
t272 = t270 * t171;
t232 = t85 * t171;
t271 = t89 * t171;
t214 = t171 * qJD(4);
t268 = -t214 + t190;
t266 = -t89 * t85 + t150;
t25 = t26 * t208;
t45 = t246 * t152;
t194 = t147 * t47 - t45;
t8 = t29 * qJD(5) + t194;
t265 = -t10 * t171 + t8 * t146 + t25;
t24 = t26 * t209;
t264 = -t174 * t171 + t24;
t55 = -pkin(5) * t171 - pkin(9) * t85;
t2 = -t174 * qJD(6) + t146 * t20 - t148 * t150;
t263 = t2 + t280;
t158 = t246 * t168;
t163 = qJD(4) * t167;
t65 = -t147 * t163 + t207 * t158 - t167 * t210;
t261 = t65 * t207;
t51 = t148 * t53;
t258 = t209 * t270 - t51;
t257 = t10 * t148 + t146 * t174;
t165 = t147 * t168;
t91 = -t167 * t246 - t165;
t40 = t91 * t53;
t178 = -t270 * t65 + t40;
t202 = t91 * t209;
t255 = -t148 * t178 + t202 * t270;
t254 = t256 ^ 2;
t253 = qJD(4) ^ 2;
t231 = pkin(7) + t135;
t124 = t231 * t142;
t125 = t231 * t144;
t87 = -t247 * t124 - t245 * t125;
t75 = pkin(8) * t167 + t87;
t88 = -t245 * t124 + t247 * t125;
t76 = -t168 * pkin(8) + t88;
t43 = t147 * t76 - t246 * t75;
t252 = t8 * t43;
t90 = -t147 * t167 + t158;
t251 = t8 * t90;
t250 = t8 * t91;
t243 = pkin(4) * t256;
t1 = t2 * t148;
t240 = t53 * t90;
t239 = t65 * t77;
t238 = t65 * t79;
t236 = t79 * t77;
t229 = t148 * t239 - t91 * t35;
t164 = qJD(4) * t168;
t66 = -qJD(5) * t165 - t147 * t164 - t246 * t163 - t167 * t195;
t228 = -t36 * t90 + t79 * t66;
t227 = -t65 * t85 - t40;
t224 = -t112 * t167 + t121 * t164;
t223 = pkin(4) * qJD(5);
t217 = qJD(6) * t77;
t215 = t256 * t121;
t212 = t142 ^ 2 + t144 ^ 2;
t201 = t91 * t208;
t193 = -t36 + t217;
t189 = qJD(1) * t212;
t186 = t91 * t73;
t185 = pkin(4) * t195;
t30 = t147 * t63 + t205;
t182 = pkin(4) * t210 - t30;
t181 = -t26 * t65 + t250;
t180 = -t37 * t90 - t66 * t77;
t179 = t170 * t90 - t171 * t66;
t176 = t10 * t146 - t148 * t174;
t44 = t147 * t75 + t246 * t76;
t96 = t168 * pkin(4) + t128;
t48 = t90 * pkin(5) - t91 * pkin(9) + t96;
t18 = t146 * t48 + t148 * t44;
t17 = -t146 * t44 + t148 * t48;
t173 = (-t129 * t142 + t139) * t142 - t105 * t144;
t172 = t281 * t85 - t258;
t161 = t168 * t253;
t159 = pkin(4) * t163;
t136 = t147 * pkin(4) + pkin(9);
t157 = -t136 * t53 - t185 * t270 - t274;
t155 = -t178 * t146 - t201 * t270;
t154 = -t176 * qJD(6) - t3 * t146 + t1;
t70 = -t124 * t197 - t125 * t196 + t269;
t153 = t168 * t111 - t163 * t256;
t151 = pkin(8) * t164 + t124 * t196 - t125 * t197 - t142 * t199 - t144 * t198;
t137 = -t246 * pkin(4) - pkin(5);
t119 = t121 ^ 2;
t116 = t167 * t253;
t71 = -t88 * qJD(4) + t160;
t67 = pkin(8) * t163 + t70;
t57 = t66 * t207;
t46 = -t243 + t55;
t31 = t246 * t63 - t220;
t23 = t66 * pkin(5) + t65 * pkin(9) - t159;
t16 = t146 * t55 + t148 * t28;
t15 = -t146 * t28 + t148 * t55;
t14 = t146 * t46 + t148 * t31;
t13 = -t146 * t31 + t148 * t46;
t12 = t44 * qJD(5) + t147 * t67 - t246 * t151;
t11 = t147 * t151 + t75 * t195 - t76 * t210 + t246 * t67;
t5 = -t18 * qJD(6) - t11 * t146 + t148 * t23;
t4 = t17 * qJD(6) + t11 * t148 + t146 * t23;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t189 (t135 * t189 - t173) * qJD(3), t111 * t167 + t164 * t256, t153 + t224, -t161, -t112 * t168 - t121 * t163, t116, 0, -t128 * t112 + (-t118 * t167 + t71) * qJD(4), -t70 * qJD(4) - t128 * t111 - t118 * t164, t87 * t111 + t88 * t112 - t70 * t121 + t71 * t256 + t59 * t167 - t169 * t163 + (t68 * qJD(4) - t58) * t168, -t169 * t70 + t58 * t88 + t59 * t87 + t68 * t71, t170 * t91 + t171 * t65, -t179 + t227, -t261, -t66 * t85 + t240, -t57, 0, -t12 * t207 + t159 * t85 - t90 * t244 + t96 * t53 + t89 * t66, -t11 * t207 + t159 * t171 + t170 * t96 - t91 * t244 - t89 * t65, t11 * t85 - t12 * t171 + t150 * t90 + t170 * t43 + t28 * t65 - t29 * t66 - t44 * t53 + t250, t29 * t11 - t28 * t12 - t150 * t44 - t159 * t89 - t96 * t244 + t252, -t79 * t202 + (-t36 * t91 - t238) * t148, -t186 + (t238 + (t36 + t217) * t91) * t146 + t229, t228 - t255, t77 * t201 + (t37 * t91 - t239) * t146, t155 + t180, t270 * t66 + t240, t12 * t77 + t146 * t181 + t17 * t53 - t174 * t66 + t25 * t91 + t270 * t5 + t3 * t90 + t37 * t43, -t10 * t66 + t12 * t79 + t148 * t181 - t18 * t53 - t2 * t90 - t24 * t91 - t270 * t4 - t36 * t43, t17 * t36 - t18 * t37 - t4 * t77 - t5 * t79 + t176 * t65 + (-qJD(6) * t257 - t146 * t2 - t148 * t3) * t91, t10 * t4 + t12 * t26 + t17 * t3 - t174 * t5 + t18 * t2 + t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t161, -t153 + t224, t68 * t163 + t164 * t169 - t167 * t58 - t59 * t168, 0, 0, 0, 0, 0, 0, -t57, t261, t179 + t227, -t150 * t91 - t28 * t66 - t29 * t65 + t251, 0, 0, 0, 0, 0, 0, t155 - t180, t228 + t255, t186 + (t193 * t91 - t238) * t146 + t229, t154 * t91 - t257 * t65 + t26 * t66 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212 * qJD(1) ^ 2, t173 * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t112, t131 + (-t184 - t121) * qJD(4), -t119 - t254, -t121 * t169 - t256 * t68, 0, 0, 0, 0, 0, 0, -t190 - t214 - 0.2e1 * t267, t170 + t259, -t233 - t234, -t171 * t28 - t29 * t85 - t244, 0, 0, 0, 0, 0, 0, t172 + t237, -t148 * t270 ^ 2 + t235 - t49 (t77 * t85 + t36) * t148 + t260 + t230, t263 * t146 + t262 * t148 + t171 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, -t119 + t254, t131 + (-t184 + t121) * qJD(4), t215, 0, 0, t118 * t256 + t156, t118 * t121 + (t204 + t68) * qJD(4) - t166, 0, 0, t232, t278, t279, -t232, t268, 0, -t64 * t195 + t45 - t85 * t243 + t271 + t30 * t207 + (-qJD(5) * t61 - t207 * t223 - t47) * t147, -t171 * t243 + (-t185 + t31) * t207 + t266, t28 * t85 - t29 * t171 + t30 * t171 - t31 * t85 + (-t246 * t170 - t147 * t53 + (-t147 * t171 + t246 * t85) * qJD(5)) * pkin(4), t28 * t30 - t29 * t31 + (-t246 * t8 + t256 * t89 - t147 * t150 + (-t147 * t28 + t246 * t29) * qJD(5)) * pkin(4), t276, t283, t275, t282, t172 - t237, t272, -t13 * t270 + t137 * t37 + t182 * t77 + (-qJD(6) * t136 * t270 - t8) * t148 + t157 * t146 + t264, -t137 * t36 + (t136 * t209 + t14) * t270 + t182 * t79 + t157 * t148 + t265, t13 * t79 + t14 * t77 + t1 + (-t77 * t185 - t136 * t37 - t85 * t174 + (t136 * t79 + t174) * qJD(6)) * t148 + (t79 * t185 + t10 * t85 - t136 * t36 - t3 + (t136 * t77 - t10) * qJD(6)) * t146, -t10 * t14 + t174 * t13 + t8 * t137 - t26 * t30 + (t147 * t26 + t257 * t246) * t223 + t154 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t278, t279, -t232, t268, 0, t29 * qJD(4) - t194 + t271, t28 * t207 + t266, 0, 0, t276, t283, t275, t282, -t270 * t281 - t237 + t51, t272, -pkin(5) * t37 - pkin(9) * t226 - t146 * t274 - t8 * t148 - t15 * t270 - t29 * t77 + t264, pkin(5) * t36 + t258 * pkin(9) + t16 * t270 - t26 * t273 - t29 * t79 + t265, t15 * t79 + t16 * t77 + t1 + (t280 + (-t37 + t216) * pkin(9)) * t148 + (pkin(9) * t193 - t262) * t146, -pkin(5) * t8 + pkin(9) * t154 - t10 * t16 + t15 * t174 - t26 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t77 ^ 2 + t79 ^ 2, t270 * t77 - t36, -t236, t270 * t79 - t37, t53, -t26 * t79 + t262, t26 * t77 - t263, 0, 0;];
tauc_reg  = t6;
