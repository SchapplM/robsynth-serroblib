% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:16
% EndTime: 2021-01-15 23:41:34
% DurationCPUTime: 4.76s
% Computational Cost: add. (5874->394), mult. (16102->569), div. (0->0), fcn. (12543->10), ass. (0->190)
t179 = sin(pkin(5));
t187 = cos(qJ(2));
t244 = qJD(1) * t187;
t229 = t179 * t244;
t284 = qJD(3) - t229;
t181 = cos(pkin(5));
t245 = qJD(1) * t181;
t168 = qJD(2) + t245;
t186 = cos(qJ(3));
t183 = sin(qJ(3));
t184 = sin(qJ(2));
t246 = qJD(1) * t179;
t230 = t184 * t246;
t213 = t183 * t230;
t117 = -t186 * t168 + t213;
t119 = t168 * t183 + t186 * t230;
t178 = sin(pkin(10));
t180 = cos(pkin(10));
t79 = t180 * t117 + t119 * t178;
t279 = qJD(5) + t79;
t182 = sin(qJ(5));
t185 = cos(qJ(5));
t201 = -t117 * t178 + t180 * t119;
t60 = t182 * t201 - t185 * t284;
t283 = t279 * t60;
t212 = t183 * t229;
t269 = qJ(4) + pkin(8);
t221 = qJD(3) * t269;
t234 = pkin(1) * t245;
t134 = -pkin(7) * t230 + t187 * t234;
t199 = (pkin(2) * t184 - pkin(8) * t187) * t179;
t135 = qJD(1) * t199;
t249 = t186 * t134 + t183 * t135;
t282 = qJ(4) * t212 + qJD(4) * t186 - t183 * t221 - t249;
t216 = -t134 * t183 + t186 * t135;
t253 = t186 * t187;
t281 = qJD(4) * t183 + t186 * t221 + (pkin(3) * t184 - qJ(4) * t253) * t246 + t216;
t202 = -t182 * t284 - t185 * t201;
t280 = t202 * t279;
t146 = t178 * t183 - t180 * t186;
t101 = t146 * t229;
t142 = t146 * qJD(3);
t278 = t101 - t142;
t147 = t178 * t186 + t180 * t183;
t248 = t284 * t147;
t219 = t185 * t279;
t235 = qJD(1) * qJD(2);
t222 = t179 * t235;
t211 = t187 * t222;
t239 = qJD(3) * t186;
t92 = -qJD(3) * t213 + t168 * t239 + t186 * t211;
t241 = qJD(2) * t187;
t226 = t183 * t241;
t240 = qJD(3) * t183;
t93 = (t184 * t239 + t226) * t246 + t168 * t240;
t56 = t178 * t92 + t180 * t93;
t262 = t182 * t56;
t277 = -t219 * t279 - t262;
t268 = t178 * t282 + t281 * t180;
t267 = t281 * t178 - t180 * t282;
t165 = t184 * t234;
t137 = pkin(7) * t229 + t165;
t209 = -t137 + (-t212 + t240) * pkin(3);
t254 = t179 * t187;
t131 = pkin(7) * t254 + (pkin(1) * t184 + pkin(8)) * t181;
t132 = (-pkin(2) * t187 - pkin(8) * t184 - pkin(1)) * t179;
t250 = t186 * t131 + t183 * t132;
t172 = pkin(3) * t178 + pkin(9);
t210 = t184 * t222;
t136 = qJD(2) * t199;
t126 = qJD(1) * t136;
t255 = t179 * t184;
t169 = pkin(7) * t255;
t273 = pkin(1) * t187;
t138 = (t181 * t273 - t169) * qJD(2);
t127 = qJD(1) * t138;
t106 = pkin(8) * t168 + t137;
t113 = qJD(1) * t132;
t72 = t106 * t186 + t113 * t183;
t189 = -qJD(3) * t72 + t186 * t126 - t183 * t127;
t22 = pkin(3) * t210 - qJ(4) * t92 - qJD(4) * t119 + t189;
t198 = -t106 * t240 + t113 * t239 + t183 * t126 + t186 * t127;
t26 = -qJ(4) * t93 - qJD(4) * t117 + t198;
t5 = -t178 * t26 + t180 * t22;
t3 = -pkin(4) * t210 - t5;
t275 = (pkin(3) * t119 + pkin(4) * t201 + pkin(9) * t79 + qJD(5) * t172) * t279 + t3;
t162 = t269 * t186;
t224 = t269 * t183;
t104 = t180 * t162 - t178 * t224;
t174 = -pkin(3) * t186 - pkin(2);
t94 = pkin(4) * t146 - pkin(9) * t147 + t174;
t274 = (-t248 * pkin(4) + t278 * pkin(9) + qJD(5) * t104 - t209) * t279 - t94 * t56;
t57 = -t178 * t93 + t180 * t92;
t28 = -qJD(5) * t202 + t182 * t57 - t185 * t210;
t128 = pkin(7) * t211 + qJD(2) * t165;
t69 = pkin(3) * t93 + t128;
t14 = pkin(4) * t56 - pkin(9) * t57 + t69;
t71 = -t106 * t183 + t186 * t113;
t54 = -qJ(4) * t119 + t71;
t48 = pkin(3) * t284 + t54;
t55 = -qJ(4) * t117 + t72;
t50 = t180 * t55;
t21 = t178 * t48 + t50;
t16 = pkin(9) * t284 + t21;
t105 = -pkin(2) * t168 - t134;
t77 = pkin(3) * t117 + qJD(4) + t105;
t29 = pkin(4) * t79 - pkin(9) * t201 + t77;
t205 = t16 * t182 - t185 * t29;
t6 = t178 * t22 + t180 * t26;
t4 = pkin(9) * t210 + t6;
t1 = -qJD(5) * t205 + t14 * t182 + t185 * t4;
t272 = t60 * t201;
t271 = t202 * t201;
t144 = t181 * t183 + t186 * t255;
t190 = -t250 * qJD(3) + t186 * t136 - t138 * t183;
t243 = qJD(2) * t184;
t228 = t179 * t243;
t143 = -t181 * t186 + t183 * t255;
t227 = t179 * t241;
t97 = -qJD(3) * t143 + t186 * t227;
t33 = pkin(3) * t228 - qJ(4) * t97 - qJD(4) * t144 + t190;
t197 = -t131 * t240 + t132 * t239 + t183 * t136 + t186 * t138;
t96 = qJD(3) * t144 + t179 * t226;
t37 = -qJ(4) * t96 - qJD(4) * t143 + t197;
t12 = t178 * t33 + t180 * t37;
t217 = -t131 * t183 + t186 * t132;
t59 = -pkin(3) * t254 - qJ(4) * t144 + t217;
t67 = -qJ(4) * t143 + t250;
t36 = t178 * t59 + t180 * t67;
t266 = pkin(4) * t230 + t268;
t265 = t104 * t56;
t264 = t178 * t55;
t237 = qJD(5) * t185;
t238 = qJD(5) * t182;
t27 = t182 * t210 + t185 * t57 - t201 * t238 + t237 * t284;
t263 = t182 * t27;
t261 = t117 * t284;
t260 = t119 * t284;
t259 = t147 * t185;
t258 = t284 * t183;
t257 = t284 * t186;
t175 = t179 ^ 2;
t256 = t175 * qJD(1) ^ 2;
t139 = t181 * pkin(1) * t243 + pkin(7) * t227;
t247 = t184 ^ 2 - t187 ^ 2;
t242 = qJD(2) * t186;
t236 = qJD(2) - t168;
t232 = t184 * t256;
t231 = t182 * t254;
t223 = t175 * t235;
t215 = t168 + t245;
t214 = 0.2e1 * t223;
t82 = pkin(3) * t96 + t139;
t207 = -0.2e1 * pkin(1) * t223;
t8 = t16 * t185 + t182 * t29;
t11 = -t178 * t37 + t180 * t33;
t20 = t180 * t48 - t264;
t35 = -t178 * t67 + t180 * t59;
t32 = -pkin(9) * t254 + t36;
t88 = t180 * t143 + t144 * t178;
t89 = -t143 * t178 + t144 * t180;
t130 = t169 + (-pkin(2) - t273) * t181;
t91 = pkin(3) * t143 + t130;
t45 = pkin(4) * t88 - pkin(9) * t89 + t91;
t204 = t182 * t45 + t185 * t32;
t203 = -t182 * t32 + t185 * t45;
t200 = t185 * t56 + (-t182 * t79 - t238) * t279;
t74 = t182 * t89 + t185 * t254;
t86 = -t101 * t182 - t185 * t230;
t196 = -t142 * t182 + t147 * t237 - t86;
t87 = -t101 * t185 + t182 * t230;
t195 = -t142 * t185 - t147 * t238 - t87;
t15 = -pkin(4) * t284 - t20;
t25 = t180 * t54 - t264;
t192 = -t172 * t56 + (t15 + t25) * t279;
t2 = -qJD(5) * t8 + t185 * t14 - t182 * t4;
t191 = -t265 + t3 * t147 + (pkin(9) * t230 - qJD(5) * t94 + t267) * t279;
t173 = -pkin(3) * t180 - pkin(4);
t103 = t162 * t178 + t180 * t224;
t75 = t185 * t89 - t231;
t66 = -t178 * t96 + t180 * t97;
t65 = t178 * t97 + t180 * t96;
t39 = -qJD(5) * t231 + t182 * t66 - t185 * t228 + t89 * t237;
t38 = -qJD(5) * t74 + t182 * t228 + t185 * t66;
t31 = pkin(4) * t254 - t35;
t24 = t178 * t54 + t50;
t17 = pkin(4) * t65 - pkin(9) * t66 + t82;
t10 = pkin(9) * t228 + t12;
t9 = -pkin(4) * t228 - t11;
t7 = [0, 0, 0, t184 * t187 * t214, -t247 * t214, t215 * t227, -t215 * t228, 0, -t128 * t181 - t139 * t168 + t184 * t207, -t127 * t181 - t138 * t168 + t187 * t207, t119 * t97 + t144 * t92, -t117 * t97 - t119 * t96 - t143 * t92 - t144 * t93, t284 * t97 + (-t187 * t92 + (qJD(1) * t144 + t119) * t243) * t179, -t284 * t96 + (t187 * t93 + (-qJD(1) * t143 - t117) * t243) * t179, (-t175 * t244 + t179 * t284) * t243, t190 * t284 + t139 * t117 + t130 * t93 + t128 * t143 + t105 * t96 + (-t189 * t187 + (qJD(1) * t217 + t71) * t243) * t179, -t197 * t284 + t139 * t119 + t130 * t92 + t128 * t144 + t105 * t97 + (t198 * t187 + (-t250 * qJD(1) - t72) * t243) * t179, t11 * t284 + t56 * t91 + t65 * t77 + t69 * t88 + t79 * t82 + (-t187 * t5 + (qJD(1) * t35 + t20) * t243) * t179, -t12 * t284 + t57 * t91 + t66 * t77 + t69 * t89 + t201 * t82 + (t187 * t6 + (-qJD(1) * t36 - t21) * t243) * t179, -t11 * t201 - t12 * t79 - t20 * t66 - t21 * t65 - t35 * t57 - t36 * t56 - t5 * t89 - t6 * t88, t11 * t20 + t12 * t21 + t35 * t5 + t36 * t6 + t69 * t91 + t77 * t82, -t202 * t38 + t27 * t75, t202 * t39 - t27 * t74 - t28 * t75 - t38 * t60, -t202 * t65 + t27 * t88 + t279 * t38 + t56 * t75, -t279 * t39 - t28 * t88 - t56 * t74 - t60 * t65, t279 * t65 + t56 * t88, (-qJD(5) * t204 - t10 * t182 + t17 * t185) * t279 + t203 * t56 + t2 * t88 - t205 * t65 + t9 * t60 + t31 * t28 + t3 * t74 + t15 * t39, -(qJD(5) * t203 + t10 * t185 + t17 * t182) * t279 - t204 * t56 - t1 * t88 - t8 * t65 - t9 * t202 + t31 * t27 + t3 * t75 + t15 * t38; 0, 0, 0, -t187 * t232, t247 * t256, t236 * t229, -t236 * t230, 0, pkin(1) * t232 + t137 * t168 - t128, pkin(7) * t210 + t134 * t168 + (-t181 * t235 + t256) * t273, t119 * t257 + t183 * t92, (t92 - t261) * t186 + (-t260 - t93) * t183, t284 * t239 + (-t284 * t253 + (qJD(2) * t183 - t119) * t184) * t246, -t284 * t240 + (t187 * t258 + (t117 + t242) * t184) * t246, -t284 * t230, -pkin(2) * t93 - t128 * t186 - t216 * t284 - t137 * t117 + (-pkin(8) * t257 + t105 * t183) * qJD(3) + (-t71 * t184 + (-pkin(8) * t243 - t105 * t187) * t183) * t246, -pkin(2) * t92 + t128 * t183 + t249 * t284 - t137 * t119 + (pkin(8) * t258 + t105 * t186) * qJD(3) + (-t105 * t253 + (-pkin(8) * t242 + t72) * t184) * t246, t146 * t69 + t174 * t56 + t209 * t79 + t248 * t77 - t268 * t284 + (-qJD(2) * t103 - t20) * t230, t147 * t69 + t174 * t57 + t209 * t201 + t278 * t77 + t267 * t284 + (-qJD(2) * t104 + t21) * t230, t103 * t57 - t146 * t6 - t147 * t5 - t20 * t278 + t201 * t268 - t248 * t21 + t267 * t79 - t265, -t103 * t5 + t104 * t6 + t174 * t69 - t268 * t20 + t209 * t77 - t267 * t21, -t195 * t202 + t27 * t259, t60 * t87 - t202 * t86 - (t182 * t202 - t185 * t60) * t142 + (-t263 - t185 * t28 + (t182 * t60 + t185 * t202) * qJD(5)) * t147, t146 * t27 + t195 * t279 - t202 * t248 + t56 * t259, -t146 * t28 - t147 * t262 - t196 * t279 - t248 * t60, t146 * t56 + t248 * t279, t103 * t28 + t2 * t146 + t196 * t15 + t191 * t182 - t274 * t185 - t205 * t248 + t266 * t60, -t1 * t146 + t103 * t27 + t195 * t15 + t274 * t182 + t191 * t185 - t202 * t266 - t248 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t117, -t117 ^ 2 + t119 ^ 2, t92 + t261, -t93 + t260, t210, -t105 * t119 + t284 * t72 + t189, t105 * t117 + t284 * t71 - t198, t284 * t24 - t77 * t201 + (-t119 * t79 + t180 * t210) * pkin(3) + t5, t284 * t25 + t77 * t79 + (-t119 * t201 - t178 * t210) * pkin(3) - t6, (-t178 * t56 - t180 * t57) * pkin(3) + (t21 - t24) * t201 + (-t20 + t25) * t79, t20 * t24 - t21 * t25 + (-t119 * t77 + t178 * t6 + t180 * t5) * pkin(3), -t202 * t219 + t263, (t27 - t283) * t185 + (-t28 + t280) * t182, t271 - t277, t200 + t272, -t279 * t201, t173 * t28 + t192 * t182 - t275 * t185 + t201 * t205 - t24 * t60, t173 * t27 + t275 * t182 + t192 * t185 + t8 * t201 + t202 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201 * t284 + t56, -t284 * t79 + t57, -t201 ^ 2 - t79 ^ 2, t20 * t201 + t21 * t79 + t69, 0, 0, 0, 0, 0, t200 - t272, t271 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202 * t60, t202 ^ 2 - t60 ^ 2, t27 + t283, -t28 - t280, t56, t15 * t202 + t279 * t8 + t2, t15 * t60 - t205 * t279 - t1;];
tauc_reg = t7;
