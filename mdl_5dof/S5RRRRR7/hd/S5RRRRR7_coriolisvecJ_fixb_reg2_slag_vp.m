% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:27
% EndTime: 2019-12-31 22:22:40
% DurationCPUTime: 4.76s
% Computational Cost: add. (10626->425), mult. (27135->569), div. (0->0), fcn. (19721->8), ass. (0->218)
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t235 = qJD(5) * t177;
t178 = cos(qJ(2));
t283 = cos(qJ(3));
t227 = t283 * t178;
t209 = qJD(1) * t227;
t175 = sin(qJ(3));
t176 = sin(qJ(2));
t239 = qJD(1) * t176;
t225 = t175 * t239;
t135 = t209 - t225;
t145 = t175 * t178 + t283 * t176;
t136 = t145 * qJD(1);
t174 = sin(qJ(4));
t282 = cos(qJ(4));
t196 = t174 * t135 + t282 * t136;
t170 = qJD(2) + qJD(3);
t219 = qJD(4) + t170;
t85 = t173 * t219 + t177 * t196;
t250 = qJD(5) * t85;
t241 = t170 * t209;
t103 = t170 * t225 - t241;
t298 = t170 * t145;
t104 = t298 * qJD(1);
t222 = qJD(4) * t282;
t237 = qJD(4) * t174;
t47 = t282 * t103 + t174 * t104 - t135 * t222 + t136 * t237;
t36 = -t173 * t47 + t250;
t201 = t177 * t219;
t83 = t173 * t196 - t201;
t291 = t173 * t36 + t83 * t235;
t96 = -t282 * t135 + t174 * t136;
t300 = t177 * t96;
t297 = qJD(5) + t96;
t302 = t173 * t297;
t236 = qJD(5) * t173;
t35 = -qJD(5) * t201 + t177 * t47 + t196 * t236;
t33 = t35 * t177;
t6 = -t300 * t83 - t302 * t85 - t291 - t33;
t218 = t174 * t103 - t282 * t104;
t48 = qJD(4) * t196 - t218;
t44 = t177 * t48;
t10 = t196 * t83 - t297 * t302 + t44;
t32 = t35 * t173;
t12 = -t32 + (t235 + t300) * t85;
t270 = t173 * t48 + t235 * t297;
t11 = -t196 * t85 + t297 * t300 + t270;
t299 = t96 * t196;
t286 = -pkin(7) - pkin(6);
t155 = t286 * t176;
t149 = qJD(1) * t155;
t267 = qJD(2) * pkin(2);
t141 = t149 + t267;
t156 = t286 * t178;
t151 = qJD(1) * t156;
t228 = t283 * t151;
t106 = t175 * t141 - t228;
t229 = qJD(2) * t286;
t210 = qJD(1) * t229;
t142 = t176 * t210;
t143 = t178 * t210;
t214 = -t175 * t142 + t283 * t143;
t70 = -qJD(3) * t106 + t214;
t53 = t103 * pkin(8) + t70;
t51 = t282 * t53;
t223 = qJD(3) * t283;
t238 = qJD(3) * t175;
t211 = -t141 * t223 - t283 * t142 - t175 * t143 - t151 * t238;
t52 = -t104 * pkin(8) - t211;
t220 = t174 * t52 - t51;
t167 = -t178 * pkin(2) - pkin(1);
t154 = qJD(1) * t167;
t114 = -t135 * pkin(3) + t154;
t266 = t114 * t196;
t296 = -t220 - t266;
t39 = t196 ^ 2 - t96 ^ 2;
t68 = pkin(4) * t196 + t96 * pkin(9);
t37 = t219 * t96 - t47;
t137 = t175 * t151;
t105 = t283 * t141 + t137;
t130 = t136 * pkin(8);
t81 = t105 - t130;
t73 = t170 * pkin(3) + t81;
t281 = t135 * pkin(8);
t82 = t106 + t281;
t188 = -t174 * t53 - t73 * t222 + t82 * t237 - t282 * t52;
t186 = t114 * t96 + t188;
t234 = qJD(1) * qJD(2);
t294 = -0.2e1 * t234;
t231 = t282 * t82;
t50 = t174 * t73 + t231;
t46 = t219 * pkin(9) + t50;
t56 = t96 * pkin(4) - pkin(9) * t196 + t114;
t22 = t173 * t56 + t177 * t46;
t221 = t176 * t234;
t162 = pkin(2) * t221;
t88 = t104 * pkin(3) + t162;
t15 = t48 * pkin(4) + t47 * pkin(9) + t88;
t3 = -qJD(5) * t22 + t177 * t15 + t173 * t188;
t292 = -t22 * t297 - t3;
t274 = t297 * t196;
t166 = t283 * pkin(2) + pkin(3);
t246 = t174 * t175;
t101 = t166 * t222 + (-t175 * t237 + (t283 * t282 - t246) * qJD(3)) * pkin(2);
t112 = -t175 * t149 + t228;
t191 = t112 - t281;
t113 = t283 * t149 + t137;
t87 = -t130 + t113;
t58 = t174 * t191 + t282 * t87;
t252 = t101 - t58;
t226 = t282 * t175;
t251 = -t174 * t87 + t282 * t191 + t166 * t237 + (t175 * t222 + (t283 * t174 + t226) * qJD(3)) * pkin(2);
t116 = t175 * t155 - t283 * t156;
t202 = t173 * t46 - t177 * t56;
t289 = t173 * t202 + t177 * t22;
t259 = t174 * t82;
t49 = t282 * t73 - t259;
t45 = -t219 * pkin(4) - t49;
t40 = t45 * t236;
t224 = t196 * t202 + t40;
t9 = t50 * qJD(4) + t220;
t284 = t9 * t173 + t45 * t235;
t213 = t196 * t22 + t284;
t38 = t196 * t170 + t218;
t2 = -t202 * qJD(5) + t173 * t15 - t177 * t188;
t115 = t283 * t155 + t175 * t156;
t91 = -t145 * pkin(8) + t115;
t198 = t175 * t176 - t227;
t92 = -pkin(8) * t198 + t116;
t66 = t174 * t92 - t282 * t91;
t285 = t9 * t66;
t280 = t136 * pkin(3);
t1 = t2 * t177;
t279 = t202 * t297;
t277 = t45 * t96;
t276 = t85 * t83;
t275 = t9 * t177;
t268 = pkin(3) * qJD(4);
t264 = t173 * t22;
t263 = t173 * t45;
t261 = t173 * t83;
t184 = t170 * t198;
t190 = t282 * t198;
t60 = qJD(4) * t190 + t145 * t237 + t174 * t298 + t282 * t184;
t257 = t177 * t60;
t256 = t177 * t85;
t34 = t36 * t177;
t110 = t174 * t145 + t190;
t253 = t48 * t110;
t249 = qJD(5) * t297;
t248 = t136 * t135;
t247 = t154 * t136;
t180 = qJD(1) ^ 2;
t245 = t178 * t180;
t179 = qJD(2) ^ 2;
t244 = t179 * t176;
t243 = t179 * t178;
t132 = pkin(2) * t226 + t174 * t166;
t240 = t176 ^ 2 - t178 ^ 2;
t233 = t202 * t300 - t264 * t96 + t1;
t169 = t176 * t267;
t168 = pkin(2) * t239;
t230 = t176 * t245;
t215 = pkin(1) * t294;
t212 = pkin(3) * t222;
t208 = t178 * t221;
t54 = t174 * t81 + t231;
t207 = pkin(3) * t237 - t54;
t206 = t196 * t50 - t49 * t96;
t204 = -t177 * t202 + t264;
t193 = t174 * t198;
t111 = t282 * t145 - t193;
t119 = pkin(3) * t198 + t167;
t65 = t110 * pkin(4) - t111 * pkin(9) + t119;
t67 = t174 * t91 + t282 * t92;
t29 = -t173 * t67 + t177 * t65;
t30 = t173 * t65 + t177 * t67;
t199 = t236 * t297 - t44;
t150 = t176 * t229;
t152 = t178 * t229;
t197 = -t175 * t150 + t283 * t152;
t62 = t280 + t68;
t195 = -t154 * t135 + t211;
t71 = t283 * t150 + t175 * t152 + t155 * t223 + t156 * t238;
t131 = -pkin(2) * t246 + t282 * t166;
t128 = pkin(9) + t132;
t194 = -t101 * t297 - t128 * t48 + t277;
t192 = t198 * qJD(2);
t164 = t174 * pkin(3) + pkin(9);
t189 = -t164 * t48 - t212 * t297 + t277;
t187 = -t204 * qJD(5) - t3 * t173 + t1;
t100 = pkin(3) * t298 + t169;
t181 = t184 * pkin(8) - t155 * t238 + t156 * t223 + t197;
t165 = -t282 * pkin(3) - pkin(4);
t127 = -pkin(4) - t131;
t117 = t168 + t280;
t86 = -t135 ^ 2 + t136 ^ 2;
t76 = t136 * t170 - t104;
t75 = t241 + (-t135 - t225) * t170;
t72 = -t116 * qJD(3) + t197;
t64 = -pkin(8) * t298 + t71;
t61 = -qJD(4) * t193 + t145 * t222 - t174 * t184 + t282 * t298;
t59 = t168 + t62;
t55 = t282 * t81 - t259;
t28 = t173 * t68 + t177 * t49;
t27 = -t173 * t49 + t177 * t68;
t26 = t173 * t59 + t177 * t58;
t25 = -t173 * t58 + t177 * t59;
t24 = t173 * t62 + t177 * t55;
t23 = -t173 * t55 + t177 * t62;
t20 = t61 * pkin(4) + t60 * pkin(9) + t100;
t17 = t67 * qJD(4) + t174 * t64 - t282 * t181;
t16 = t174 * t181 + t91 * t222 - t92 * t237 + t282 * t64;
t13 = t261 * t297 - t34;
t5 = -t30 * qJD(5) - t173 * t16 + t177 * t20;
t4 = t29 * qJD(5) + t177 * t16 + t173 * t20;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t208, t240 * t294, t243, -0.2e1 * t208, -t244, 0, -pkin(6) * t243 + t176 * t215, pkin(6) * t244 + t178 * t215, 0, 0, -t103 * t145 - t136 * t184, t103 * t198 - t145 * t104 + t170 * (-t135 * t198 - t136 * t145), -t184 * t170, t104 * t198 - t135 * t298, -t298 * t170, 0, t167 * t104 - t135 * t169 + t154 * t298 + t192 * t168 + t72 * t170, -t167 * t103 + t136 * t169 + t145 * t162 - t154 * t184 - t71 * t170, t71 * t135 - t116 * t104 + t211 * t198 - t106 * t298 - t72 * t136 + t115 * t103 - t70 * t145 + t105 * (qJD(3) * t198 + t192), t105 * t72 + t106 * t71 + t70 * t115 - t116 * t211 + 0.2e1 * t154 * t169, -t47 * t111 - t196 * t60, t47 * t110 - t111 * t48 - t196 * t61 + t60 * t96, -t60 * t219, t96 * t61 + t253, -t61 * t219, 0, t100 * t96 + t88 * t110 + t114 * t61 + t119 * t48 - t17 * t219, t100 * t196 + t88 * t111 - t114 * t60 - t119 * t47 - t16 * t219, t110 * t188 + t9 * t111 - t16 * t96 + t17 * t196 - t66 * t47 - t67 * t48 + t49 * t60 - t50 * t61, t114 * t100 + t88 * t119 + t50 * t16 - t49 * t17 - t188 * t67 + t285, -t60 * t256 + (-t236 * t85 - t33) * t111, (t173 * t85 + t177 * t83) * t60 + (t32 - t34 + (-t256 + t261) * qJD(5)) * t111, -t35 * t110 - t111 * t199 - t257 * t297 + t85 * t61, t291 * t111 - t60 * t261, -t36 * t110 - t111 * t270 + t302 * t60 - t83 * t61, t297 * t61 + t253, t3 * t110 + t284 * t111 + t17 * t83 - t202 * t61 - t60 * t263 + t29 * t48 + t297 * t5 + t66 * t36, -t45 * t257 - t2 * t110 + t17 * t85 - t22 * t61 - t30 * t48 - t66 * t35 - t4 * t297 + (-t40 + t275) * t111, t29 * t35 - t30 * t36 - t4 * t83 - t5 * t85 + t204 * t60 + (-qJD(5) * t289 - t173 * t2 - t177 * t3) * t111, t17 * t45 + t2 * t30 - t202 * t5 + t22 * t4 + t29 * t3 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, t240 * t180, 0, t230, 0, 0, t180 * pkin(1) * t176, pkin(1) * t245, 0, 0, -t248, t86, t75, t248, t76, 0, t135 * t168 - t112 * t170 - t247 + (t228 + (-pkin(2) * t170 - t141) * t175) * qJD(3) + t214, t113 * t170 + (-t136 * t239 - t170 * t223) * pkin(2) + t195, (t106 + t112) * t136 + (t105 - t113) * t135 + (t103 * t283 - t104 * t175 + (t283 * t135 + t136 * t175) * qJD(3)) * pkin(2), -t105 * t112 - t106 * t113 + (-t154 * t239 + t283 * t70 - t175 * t211 + (-t105 * t175 + t283 * t106) * qJD(3)) * pkin(2), t299, t39, t37, -t299, t38, 0, -t117 * t96 + (-t50 - t251) * qJD(4) - t251 * t170 + t296, -t117 * t196 - t252 * t219 + t186, t131 * t47 - t132 * t48 + t196 * t251 - t252 * t96 + t206, -t114 * t117 - t9 * t131 - t132 * t188 - t251 * t49 + t252 * t50, t12, t6, t11, t13, t10, -t274, t127 * t36 - t25 * t297 + t251 * t83 + (-t128 * t249 - t9) * t177 + t194 * t173 + t224, -t127 * t35 + (t128 * t236 + t26) * t297 + t251 * t85 + t194 * t177 + t213, t25 * t85 + t26 * t83 + (-t101 * t83 - t128 * t36 + (t128 * t85 + t202) * qJD(5)) * t177 + (t101 * t85 - t128 * t35 - t3 + (t128 * t83 - t22) * qJD(5)) * t173 + t233, t101 * t289 + t9 * t127 + t128 * t187 + t202 * t25 - t22 * t26 + t251 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248, t86, t75, t248, t76, 0, t106 * t170 - t247 + t70, t105 * t170 + t195, 0, 0, t299, t39, t37, -t299, t38, 0, -t82 * t222 + t51 - t96 * t280 - t266 + t54 * t219 + (-qJD(4) * t73 - t219 * t268 - t52) * t174, t55 * t219 + (-t136 * t196 - t219 * t222) * pkin(3) + t186, -t54 * t196 + t55 * t96 + (t282 * t47 - t174 * t48 + (t174 * t196 - t282 * t96) * qJD(4)) * pkin(3) + t206, t49 * t54 - t50 * t55 + (-t282 * t9 - t114 * t136 - t174 * t188 + (-t174 * t49 + t282 * t50) * qJD(4)) * pkin(3), t12, t6, t11, t13, t10, -t274, t165 * t36 - t23 * t297 + t207 * t83 + (-t164 * t249 - t9) * t177 + t189 * t173 + t224, -t165 * t35 + (t164 * t236 + t24) * t297 + t207 * t85 + t189 * t177 + t213, t23 * t85 + t24 * t83 + (-t83 * t212 - t164 * t36 + (t164 * t85 + t202) * qJD(5)) * t177 + (t85 * t212 - t164 * t35 - t3 + (t164 * t83 - t22) * qJD(5)) * t173 + t233, t9 * t165 + t202 * t23 - t22 * t24 - t45 * t54 + (t174 * t45 + t282 * t289) * t268 + t187 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t39, t37, -t299, t38, 0, t170 * t50 + t296, t219 * t49 + t186, 0, 0, t12, t6, t11, t302 * t83 - t34, t10, -t274, -pkin(4) * t36 - pkin(9) * t270 + t96 * t263 - t27 * t297 - t50 * t83 + t224 - t275, pkin(4) * t35 + pkin(9) * t199 + t28 * t297 + t300 * t45 - t50 * t85 + t213, t27 * t85 + t28 * t83 + t1 + (t279 + (-t36 + t250) * pkin(9)) * t177 + ((qJD(5) * t83 - t35) * pkin(9) + t292) * t173, -t9 * pkin(4) + pkin(9) * t187 + t202 * t27 - t22 * t28 - t45 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, -t83 ^ 2 + t85 ^ 2, t297 * t83 - t35, -t276, t297 * t85 - t36, t48, -t45 * t85 - t292, t45 * t83 - t2 - t279, 0, 0;];
tauc_reg = t7;
