% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:20
% EndTime: 2019-03-09 02:36:31
% DurationCPUTime: 4.46s
% Computational Cost: add. (9406->422), mult. (20618->561), div. (0->0), fcn. (14960->8), ass. (0->205)
t171 = sin(pkin(10));
t275 = sin(qJ(4));
t217 = t275 * t171;
t172 = cos(pkin(10));
t178 = cos(qJ(4));
t240 = t172 * t178;
t137 = -t217 + t240;
t175 = sin(qJ(5));
t177 = cos(qJ(5));
t228 = qJD(5) * t177;
t229 = qJD(5) * t175;
t136 = t178 * t171 + t275 * t172;
t185 = t136 * qJD(3);
t173 = -pkin(1) - qJ(3);
t283 = t173 * qJD(1);
t151 = qJD(2) + t283;
t207 = -pkin(7) * qJD(1) + t151;
t122 = t207 * t171;
t123 = t207 * t172;
t284 = -t275 * t122 + t178 * t123;
t48 = -qJD(1) * t185 + qJD(4) * t284;
t80 = t178 * t122 + t275 * t123;
t74 = qJD(4) * pkin(8) + t80;
t216 = qJD(1) * t240;
t129 = -qJD(1) * t217 + t216;
t170 = qJD(1) * qJ(2);
t164 = qJD(3) + t170;
t166 = t171 * pkin(3);
t144 = qJD(1) * t166 + t164;
t184 = qJD(1) * t136;
t75 = pkin(4) * t184 - pkin(8) * t129 + t144;
t212 = qJD(4) * t275;
t199 = qJD(1) * t212;
t148 = t171 * t199;
t230 = qJD(4) * t178;
t211 = qJD(1) * t230;
t121 = t172 * t211 - t148;
t169 = qJD(1) * qJD(2);
t233 = t171 * t211 + t172 * t199;
t76 = t121 * pkin(4) + t233 * pkin(8) + t169;
t15 = t175 * t76 + t177 * t48 + t75 * t228 - t74 * t229;
t291 = qJD(5) + t184;
t37 = -t175 * t74 + t177 * t75;
t195 = -t291 * t37 + t15;
t38 = t175 * t75 + t177 * t74;
t16 = -qJD(5) * t38 - t175 * t48 + t177 * t76;
t281 = t291 * t38 + t16;
t109 = -t177 * qJD(4) + t129 * t175;
t111 = qJD(4) * t175 + t129 * t177;
t174 = sin(qJ(6));
t176 = cos(qJ(6));
t194 = t109 * t174 - t176 * t111;
t54 = t176 * t109 + t111 * t174;
t273 = t54 * t194;
t235 = t176 * t177;
t239 = t174 * t175;
t138 = -t235 + t239;
t222 = qJD(5) + qJD(6);
t226 = qJD(6) * t176;
t260 = t138 * t184 - t176 * t228 - t177 * t226 + t222 * t239;
t238 = t174 * t177;
t139 = t175 * t176 + t238;
t102 = t222 * t139;
t259 = t139 * t184 + t102;
t276 = -pkin(9) - pkin(8);
t219 = qJD(5) * t276;
t247 = t184 * t175;
t95 = pkin(4) * t129 + pkin(8) * t184;
t42 = t175 * t95 + t177 * t284;
t297 = -pkin(9) * t247 + t175 * t219 - t42;
t41 = -t175 * t284 + t177 * t95;
t296 = pkin(5) * t129 + t41 + (pkin(9) * t184 - t219) * t177;
t295 = t137 * qJD(3);
t294 = qJD(4) * qJD(5) - t233;
t293 = t229 + t247;
t204 = t177 * t291;
t237 = t175 * t121;
t292 = -t204 * t291 - t237;
t290 = t194 ^ 2 - t54 ^ 2;
t65 = t129 * t229 - t177 * t294;
t6 = t121 * pkin(5) + t65 * pkin(9) + t16;
t27 = -pkin(9) * t111 + t37;
t25 = pkin(5) * t291 + t27;
t28 = -pkin(9) * t109 + t38;
t265 = t176 * t28;
t8 = t174 * t25 + t265;
t220 = t129 * t228 + t175 * t294;
t9 = -t220 * pkin(9) + t15;
t2 = -qJD(6) * t8 - t174 * t9 + t176 * t6;
t73 = -qJD(4) * pkin(4) - t284;
t46 = pkin(5) * t109 + t73;
t289 = t46 * t194 + t2;
t120 = qJD(6) + t291;
t227 = qJD(6) * t174;
t23 = t109 * t226 + t111 * t227 + t174 * t220 + t176 * t65;
t288 = t120 * t54 - t23;
t1 = (qJD(6) * t25 + t9) * t176 + t174 * t6 - t28 * t227;
t287 = t46 * t54 - t1;
t181 = qJD(6) * t194 + t174 * t65 - t176 * t220;
t286 = -t120 * t194 + t181;
t90 = t138 * t136;
t88 = t139 * t136;
t272 = -pkin(7) + t173;
t142 = t272 * t171;
t143 = t272 * t172;
t99 = t275 * t142 - t178 * t143;
t232 = t171 ^ 2 + t172 ^ 2;
t282 = t232 * qJD(3);
t280 = -t138 * t23 - t194 * t259;
t279 = -t120 * t260 + t139 * t121;
t131 = t136 * qJD(4);
t278 = -t129 * t131 - t137 * t233;
t277 = t129 ^ 2;
t221 = 0.2e1 * t169;
t190 = t295 * qJD(1);
t49 = t80 * qJD(4) + t190;
t274 = t49 * t99;
t145 = t276 * t175;
t146 = t276 * t177;
t108 = t145 * t174 - t146 * t176;
t271 = qJD(6) * t108 + t174 * t297 + t176 * t296;
t107 = t145 * t176 + t146 * t174;
t270 = -qJD(6) * t107 + t174 * t296 - t176 * t297;
t157 = qJ(2) + t166;
t94 = pkin(4) * t136 - pkin(8) * t137 + t157;
t100 = t178 * t142 + t275 * t143;
t96 = t177 * t100;
t45 = t175 * t94 + t96;
t269 = t129 * t54;
t267 = t174 * t28;
t266 = t175 * t49;
t264 = t49 * t137;
t263 = t49 * t177;
t262 = t194 * t129;
t261 = t65 * t175;
t132 = -t171 * t212 + t172 * t230;
t258 = -t138 * qJD(1) + t139 * t132 - t222 * t90;
t257 = t139 * qJD(1) + t138 * t132 + t222 * t88;
t59 = t175 * t220;
t256 = -t109 * t228 - t59;
t255 = t109 * t184;
t254 = t109 * t129;
t253 = t109 * t175;
t252 = t109 * t177;
t251 = t111 * t109;
t250 = t111 * t129;
t249 = t111 * t175;
t248 = t111 * t177;
t97 = t121 * t136;
t246 = t129 * t184;
t244 = t131 * t177;
t243 = t137 * t175;
t242 = t137 * t177;
t236 = t175 * t131;
t114 = t177 * t121;
t231 = qJD(4) * t184;
t225 = t132 * qJD(4);
t215 = t137 * t229;
t67 = -t99 * qJD(4) - t185;
t92 = pkin(4) * t132 + pkin(8) * t131 + qJD(2);
t208 = -t175 * t67 + t177 * t92;
t44 = -t100 * t175 + t177 * t94;
t206 = qJD(1) * t232;
t205 = t175 * t291;
t203 = qJD(5) * t136 + qJD(1);
t202 = t139 * t181 + t260 * t54;
t201 = -t120 * t259 - t138 * t121;
t200 = pkin(5) * t293 - t80;
t198 = t220 * t177;
t36 = pkin(5) * t136 - pkin(9) * t242 + t44;
t39 = -pkin(9) * t243 + t45;
t18 = -t174 * t39 + t176 * t36;
t19 = t174 * t36 + t176 * t39;
t197 = t175 * t38 + t177 * t37;
t196 = t175 * t37 - t177 * t38;
t193 = t248 + t253;
t192 = t132 * t184 + t97;
t191 = -t291 * t293 + t114;
t188 = t137 * t228 - t236;
t187 = -t215 - t244;
t21 = -t100 * t229 + t175 * t92 + t177 * t67 + t94 * t228;
t186 = -pkin(8) * t121 + t291 * t73;
t182 = t131 * t284 - t80 * t132 - t48 * t136 + t264;
t180 = -qJD(5) * t197 + t15 * t177 - t16 * t175;
t68 = qJD(4) * t100 + t295;
t179 = qJD(1) ^ 2;
t163 = -pkin(5) * t177 - pkin(4);
t126 = t184 ^ 2;
t124 = t131 * qJD(4);
t91 = t138 * t137;
t89 = t139 * t137;
t71 = pkin(5) * t243 + t99;
t40 = pkin(5) * t188 + t68;
t34 = t220 * pkin(5) + t49;
t32 = -t131 * t238 - t174 * t215 - t227 * t243 + (t222 * t242 - t236) * t176;
t30 = t102 * t137 + t131 * t235 - t174 * t236;
t22 = -t45 * qJD(5) + t208;
t17 = -pkin(9) * t188 + t21;
t14 = pkin(9) * t244 + t132 * pkin(5) + (-t96 + (pkin(9) * t137 - t94) * t175) * qJD(5) + t208;
t11 = t176 * t27 - t267;
t10 = -t174 * t27 - t265;
t7 = t176 * t25 - t267;
t4 = -qJD(6) * t19 + t176 * t14 - t174 * t17;
t3 = qJD(6) * t18 + t174 * t14 + t176 * t17;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, qJ(2) * t221, 0, 0, 0, 0, 0, 0, t171 * t221, t172 * t221, 0.2e1 * qJD(3) * t206 (t164 + t170) * qJD(2) + (-t151 - t283) * t282, t278, -t137 * t121 - t129 * t132 + t131 * t184 + t233 * t136, -t124, t192, -t225, 0, 0.2e1 * t184 * qJD(2) - t68 * qJD(4) + t157 * t121 + t144 * t132, -t157 * t233 - t144 * t131 - t67 * qJD(4) + (qJD(1) * t137 + t129) * qJD(2), -t100 * t121 + t68 * t129 - t184 * t67 - t99 * t233 + t182, t48 * t100 + t274 + t80 * t67 - t284 * t68 + (qJD(1) * t157 + t144) * qJD(2), t111 * t187 - t65 * t242 (t249 + t252) * t131 + (-t198 + t261 + (-t248 + t253) * qJD(5)) * t137, t111 * t132 + t137 * t114 - t65 * t136 + t187 * t291, t109 * t188 + t137 * t59, -t109 * t132 - t220 * t136 - t137 * t237 - t188 * t291, t132 * t291 + t97, t22 * t291 + t44 * t121 + t16 * t136 + t37 * t132 + t68 * t109 + t99 * t220 - t73 * t236 + (t73 * t228 + t266) * t137, -t73 * t244 + t111 * t68 - t121 * t45 - t291 * t21 - t132 * t38 - t136 * t15 - t65 * t99 + (-t73 * t229 + t263) * t137, -t21 * t109 - t45 * t220 - t22 * t111 + t44 * t65 + t197 * t131 + (qJD(5) * t196 - t15 * t175 - t16 * t177) * t137, t15 * t45 + t16 * t44 + t21 * t38 + t22 * t37 + t68 * t73 + t274, t194 * t30 + t23 * t91, -t181 * t91 + t194 * t32 + t23 * t89 + t30 * t54, -t120 * t30 - t121 * t91 - t132 * t194 - t136 * t23, -t181 * t89 + t32 * t54, -t120 * t32 - t121 * t89 - t132 * t54 + t136 * t181, t120 * t132 + t97, t120 * t4 + t121 * t18 + t132 * t7 + t136 * t2 - t181 * t71 + t32 * t46 + t34 * t89 + t40 * t54, -t1 * t136 - t120 * t3 - t121 * t19 - t132 * t8 - t194 * t40 - t23 * t71 - t30 * t46 - t34 * t91, -t1 * t89 + t18 * t23 + t181 * t19 + t194 * t4 + t2 * t91 - t3 * t54 + t30 * t7 - t32 * t8, t1 * t19 + t18 * t2 + t3 * t8 + t34 * t71 + t4 * t7 + t40 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -t179 * qJ(2), 0, 0, 0, 0, 0, 0, -t179 * t171, -t179 * t172, 0 (-t164 - t282) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t184 - t124, -qJD(1) * t129 - t225, -t192 - t278, -qJD(1) * t144 - t182, 0, 0, 0, 0, 0, 0, -t136 * t237 + t131 * t109 - t137 * t220 + (-t175 * t132 - t177 * t203) * t291, -t136 * t114 + t131 * t111 + t137 * t65 + (-t177 * t132 + t175 * t203) * t291 (t249 - t252) * t132 + t193 * qJD(1) + (qJD(5) * t193 - t198 - t261) * t136, -qJD(1) * t197 + t131 * t73 - t132 * t196 + t136 * t180 - t264, 0, 0, 0, 0, 0, 0, -t120 * t258 - t88 * t121 + t131 * t54 + t137 * t181, t120 * t257 + t90 * t121 - t131 * t194 + t137 * t23, -t181 * t90 - t194 * t258 - t88 * t23 + t257 * t54, -t1 * t90 + t131 * t46 - t137 * t34 - t2 * t88 - t257 * t8 - t258 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232 * t179, t151 * t206 + t169, 0, 0, 0, 0, 0, 0, -t148 + (t129 + t216) * qJD(4), -t231 - t233, -t126 - t277, t129 * t284 + t184 * t80 + t169, 0, 0, 0, 0, 0, 0, t191 - t254, -t250 + t292 (t65 - t255) * t177 + t111 * t205 + t256, -t73 * t129 + t195 * t175 + t281 * t177, 0, 0, 0, 0, 0, 0, t201 - t269, t262 - t279, t202 + t280, t1 * t139 - t129 * t46 - t138 * t2 - t259 * t7 - t260 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t126 + t277, t231 - t233, -t246, t148 + (t129 - t216) * qJD(4), 0, -t144 * t129 - t190 (qJD(3) + t144) * t184, 0, 0, t111 * t204 - t261 (-t65 - t255) * t177 - t291 * t249 + t256, -t250 - t292, t109 * t205 - t198, t191 + t254, -t291 * t129, -pkin(4) * t220 - t263 - t37 * t129 - t80 * t109 + (-pkin(8) * t228 - t41) * t291 + t186 * t175, pkin(4) * t65 - t111 * t80 + t129 * t38 + t266 + (pkin(8) * t229 + t42) * t291 + t186 * t177, t42 * t109 + t41 * t111 + ((qJD(5) * t111 - t220) * pkin(8) + t195) * t177 + ((qJD(5) * t109 - t65) * pkin(8) - t281) * t175, -pkin(4) * t49 + pkin(8) * t180 - t37 * t41 - t38 * t42 - t73 * t80, -t23 * t139 + t194 * t260, t202 - t280, t262 + t279, -t138 * t181 + t259 * t54, t201 + t269, -t120 * t129, t107 * t121 - t120 * t271 - t129 * t7 + t138 * t34 - t163 * t181 + t200 * t54 + t259 * t46, -t108 * t121 + t120 * t270 + t129 * t8 + t139 * t34 - t163 * t23 - t194 * t200 - t260 * t46, -t1 * t138 + t107 * t23 + t108 * t181 - t139 * t2 - t194 * t271 - t259 * t8 + t260 * t7 + t270 * t54, t1 * t108 + t107 * t2 + t163 * t34 + t200 * t46 - t270 * t8 - t271 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, -t109 ^ 2 + t111 ^ 2, t109 * t291 - t65, -t251, t111 * t291 - t220, t121, -t73 * t111 + t281, t109 * t73 - t195, 0, 0, -t273, t290, t288, t273, t286, t121, -t10 * t120 + (-t111 * t54 - t120 * t227 + t121 * t176) * pkin(5) + t289, t11 * t120 + (t111 * t194 - t120 * t226 - t121 * t174) * pkin(5) + t287, -t10 * t194 + t11 * t54 - t194 * t8 - t54 * t7 + (t174 * t181 + t176 * t23 + (-t174 * t194 - t176 * t54) * qJD(6)) * pkin(5), -t7 * t10 - t8 * t11 + (t1 * t174 - t111 * t46 + t176 * t2 + (-t174 * t7 + t176 * t8) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, t290, t288, t273, t286, t121, t8 * t120 + t289, t7 * t120 + t287, 0, 0;];
tauc_reg  = t5;
