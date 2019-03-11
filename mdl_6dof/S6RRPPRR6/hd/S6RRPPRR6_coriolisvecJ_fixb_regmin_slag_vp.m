% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:16
% EndTime: 2019-03-09 09:15:24
% DurationCPUTime: 3.29s
% Computational Cost: add. (4029->339), mult. (9798->460), div. (0->0), fcn. (7248->8), ass. (0->194)
t164 = qJD(2) - qJD(5);
t175 = cos(qJ(6));
t170 = cos(pkin(10));
t174 = sin(qJ(2));
t215 = qJD(1) * qJD(2);
t210 = t174 * t215;
t142 = t170 * t210;
t169 = sin(pkin(10));
t177 = cos(qJ(2));
t209 = t177 * t215;
t100 = t169 * t209 - t142;
t124 = t174 * t169 + t177 * t170;
t116 = t124 * qJD(2);
t101 = qJD(1) * t116;
t223 = qJD(1) * t177;
t224 = qJD(1) * t174;
t107 = t169 * t224 + t170 * t223;
t211 = t169 * t223;
t110 = t170 * t224 - t211;
t173 = sin(qJ(5));
t176 = cos(qJ(5));
t219 = qJD(5) * t176;
t220 = qJD(5) * t173;
t185 = t173 * t100 - t176 * t101 + t107 * t219 + t110 * t220;
t217 = qJD(6) * t175;
t172 = sin(qJ(6));
t218 = qJD(6) * t172;
t61 = t173 * t107 - t176 * t110;
t10 = -t164 * t217 - t175 * t185 + t218 * t61;
t190 = t172 * t164 + t175 * t61;
t11 = -t190 * qJD(6) - t172 * t185;
t252 = t176 * t107 + t173 * t110;
t261 = qJD(6) + t252;
t243 = t190 * t261;
t44 = t175 * t164 - t172 * t61;
t244 = t44 * t261;
t280 = (t10 - t244) * t175 + (-t11 + t243) * t172;
t204 = t261 ^ 2;
t26 = -t61 * qJD(5) + t176 * t100 + t173 * t101;
t24 = t175 * t26;
t277 = -t172 * t204 + t24;
t237 = t10 * t172;
t273 = t175 * t261;
t276 = t190 * t273 - t237;
t235 = t172 * t26;
t266 = t190 * t61;
t275 = t261 * t273 + t235 - t266;
t272 = qJD(6) - t261;
t268 = t44 * t61;
t271 = t268 - t277;
t247 = t110 * pkin(8);
t157 = pkin(7) * t223;
t133 = -qJ(4) * t223 + t157;
t166 = qJD(2) * qJ(3);
t120 = t133 + t166;
t156 = pkin(7) * t224;
t131 = qJ(4) * t224 - t156;
t178 = -pkin(2) - pkin(3);
t212 = t178 * qJD(2);
t96 = qJD(3) + t212 - t131;
t52 = -t169 * t120 + t170 * t96;
t38 = -qJD(2) * pkin(4) - t247 + t52;
t248 = t107 * pkin(8);
t53 = t170 * t120 + t169 * t96;
t39 = t53 - t248;
t15 = t173 * t38 + t176 * t39;
t13 = -t164 * pkin(9) + t15;
t121 = -qJD(1) * pkin(1) - pkin(2) * t223 - qJ(3) * t224;
t86 = pkin(3) * t223 + qJD(4) - t121;
t60 = t107 * pkin(4) + t86;
t16 = pkin(5) * t252 + pkin(9) * t61 + t60;
t4 = t175 * t13 + t172 * t16;
t270 = t4 * t61;
t269 = t61 * pkin(5);
t267 = t261 * t61;
t192 = t172 * t13 - t175 * t16;
t265 = t192 * t61;
t264 = t61 * t164;
t263 = t61 * t252;
t262 = t164 * t252;
t260 = t252 ^ 2 - t61 ^ 2;
t222 = qJD(2) * t174;
t241 = pkin(7) - qJ(4);
t103 = -t177 * qJD(4) - t241 * t222;
t165 = qJD(2) * qJD(3);
t82 = t103 * qJD(1) + t165;
t149 = pkin(7) * t209;
t216 = t174 * qJD(4);
t221 = qJD(2) * t177;
t90 = t149 + (-qJ(4) * t221 - t216) * qJD(1);
t205 = t169 * t82 - t170 * t90;
t33 = -t101 * pkin(8) - t205;
t43 = t169 * t90 + t170 * t82;
t34 = -t100 * pkin(8) + t43;
t2 = t173 * t34 - t176 * t33 + t39 * t219 + t38 * t220;
t259 = t60 * t61 - t2;
t1 = t173 * t33 + t176 * t34 + t38 * t219 - t39 * t220;
t258 = t252 * t60 - t1;
t256 = -t26 + t264;
t254 = t185 + t262;
t253 = -0.2e1 * t215;
t154 = qJ(3) * t223;
t102 = t178 * t224 + t154;
t68 = -t110 * pkin(4) + t102;
t134 = -t169 * qJ(3) + t170 * t178;
t130 = -pkin(4) + t134;
t135 = t170 * qJ(3) + t169 * t178;
t188 = t173 * t130 + t176 * t135;
t72 = -pkin(9) + t188;
t251 = (-pkin(9) * t252 + qJD(6) * t72 + t269 + t68) * t261 - t2;
t250 = (t261 * pkin(9) - t269) * t261 + t2;
t14 = -t173 * t39 + t176 * t38;
t12 = t164 * pkin(5) - t14;
t125 = -t177 * t169 + t174 * t170;
t69 = t176 * t124 + t173 * t125;
t70 = -t173 * t124 + t176 * t125;
t138 = -t177 * pkin(2) - t174 * qJ(3) - pkin(1);
t122 = t177 * pkin(3) - t138;
t75 = t124 * pkin(4) + t122;
t20 = t69 * pkin(5) - t70 * pkin(9) + t75;
t140 = t241 * t174;
t141 = t241 * t177;
t78 = t170 * t140 - t169 * t141;
t56 = -t125 * pkin(8) + t78;
t79 = t169 * t140 + t170 * t141;
t57 = -t124 * pkin(8) + t79;
t22 = t173 * t56 + t176 * t57;
t115 = t169 * t221 - t170 * t222;
t31 = -t69 * qJD(5) - t173 * t115 + t176 * t116;
t21 = t173 * t57 - t176 * t56;
t104 = qJD(2) * t141 - t216;
t54 = -t169 * t103 + t170 * t104;
t40 = -t116 * pkin(8) + t54;
t55 = t170 * t103 + t169 * t104;
t41 = -t115 * pkin(8) + t55;
t7 = -t21 * qJD(5) + t173 * t40 + t176 * t41;
t249 = -(qJD(6) * t20 + t7) * t261 - t22 * t26 - (qJD(6) * t16 + t1) * t69 + t2 * t70 + t12 * t31;
t246 = t12 * t70;
t245 = t20 * t26;
t242 = t70 * t26;
t123 = t173 * t169 - t176 * t170;
t189 = t176 * t130 - t173 * t135;
t73 = -t169 * t131 + t170 * t133;
t47 = t73 - t248;
t74 = t170 * t131 + t169 * t133;
t48 = t74 + t247;
t240 = t123 * qJD(3) - t189 * qJD(5) + t173 * t47 + t176 * t48;
t126 = t176 * t169 + t173 * t170;
t239 = t126 * qJD(3) + t188 * qJD(5) - t173 * t48 + t176 * t47;
t238 = qJD(2) * pkin(2);
t180 = qJD(1) ^ 2;
t233 = t177 * t180;
t179 = qJD(2) ^ 2;
t232 = t179 * t174;
t231 = t179 * t177;
t229 = t164 * t123;
t228 = t164 * t126;
t160 = t174 * qJD(3);
t227 = qJ(3) * t209 + qJD(1) * t160;
t226 = qJ(3) * t221 + t160;
t167 = t174 ^ 2;
t225 = -t177 ^ 2 + t167;
t214 = t70 * t218;
t213 = t174 * t233;
t199 = pkin(1) * t253;
t198 = qJD(3) - t238;
t197 = qJD(3) * t169 + t73;
t196 = qJD(3) * t170 - t74;
t195 = qJD(1) * t138 + t121;
t194 = t174 * t212;
t193 = t261 * t31 + t242;
t191 = t169 * t52 - t170 * t53;
t187 = qJD(6) * t126 + t224;
t105 = pkin(2) * t222 - t226;
t91 = pkin(2) * t210 - t227;
t186 = -pkin(7) * t179 - qJD(1) * t105 - t91;
t85 = t194 + t226;
t77 = qJD(1) * t194 + t227;
t183 = -pkin(9) * t26 + (t12 + t14) * t261;
t59 = t115 * pkin(4) + t85;
t182 = -t72 * t26 + (-t12 + t240) * t261;
t51 = t100 * pkin(4) + t77;
t136 = -pkin(7) * t210 + t165;
t137 = t156 + t198;
t139 = t157 + t166;
t181 = t136 * t177 + (t137 * t177 + (-t139 + t157) * t174) * qJD(2);
t132 = pkin(2) * t224 - t154;
t71 = pkin(5) - t189;
t32 = t70 * qJD(5) + t176 * t115 + t173 * t116;
t9 = t32 * pkin(5) - t31 * pkin(9) + t59;
t8 = t22 * qJD(5) + t173 * t41 - t176 * t40;
t6 = t26 * pkin(5) + pkin(9) * t185 + t51;
t5 = t175 * t6;
t3 = [0, 0, 0, 0.2e1 * t174 * t209, t225 * t253, t231, -t232, 0, -pkin(7) * t231 + t174 * t199, pkin(7) * t232 + t177 * t199, t186 * t177 + t195 * t222, t181, t186 * t174 - t195 * t221, t181 * pkin(7) + t121 * t105 + t91 * t138, -t54 * qJD(2) + t122 * t100 + t85 * t107 + t86 * t115 + t77 * t124, t55 * qJD(2) + t122 * t101 + t85 * t110 + t86 * t116 + t77 * t125, -t79 * t100 - t78 * t101 - t55 * t107 - t54 * t110 - t53 * t115 - t52 * t116 - t43 * t124 + t125 * t205, t77 * t122 - t205 * t78 + t43 * t79 + t52 * t54 + t53 * t55 + t86 * t85, -t185 * t70 - t31 * t61, t185 * t69 - t252 * t31 + t32 * t61 - t242, -t31 * t164, t32 * t164, 0, t8 * t164 + t252 * t59 + t75 * t26 + t60 * t32 + t51 * t69, t7 * t164 - t185 * t75 + t60 * t31 + t51 * t70 - t59 * t61, t190 * t214 + (t10 * t70 - t190 * t31) * t175 (t172 * t190 - t175 * t44) * t31 + (-t237 - t11 * t175 + (t172 * t44 + t175 * t190) * qJD(6)) * t70, t10 * t69 + t175 * t193 - t190 * t32 - t214 * t261, -t217 * t261 * t70 - t11 * t69 - t172 * t193 - t44 * t32, t26 * t69 + t261 * t32, t21 * t11 - t192 * t32 + t8 * t44 + t5 * t69 + (t245 + t9 * t261 + (-t13 * t69 - t22 * t261 + t246) * qJD(6)) * t175 + t249 * t172, t21 * t10 - t4 * t32 - t8 * t190 + (-(-qJD(6) * t22 + t9) * t261 - t245 - (-qJD(6) * t13 + t6) * t69 - qJD(6) * t246) * t172 + t249 * t175; 0, 0, 0, -t213, t225 * t180, 0, 0, 0, t180 * pkin(1) * t174, pkin(1) * t233 (-t121 * t174 + t132 * t177) * qJD(1) ((t139 - t166) * t174 + (-t137 + t198) * t177) * qJD(1), 0.2e1 * t165 + (t121 * t177 + t132 * t174) * qJD(1), t136 * qJ(3) + t139 * qJD(3) - t121 * t132 + (t139 * t174 + (-t137 - t238) * t177) * qJD(1) * pkin(7), t197 * qJD(2) - t102 * t107 + t86 * t110 + t205, t196 * qJD(2) - t102 * t110 - t86 * t107 + t43, -t135 * t100 - t134 * t101 + (t197 - t53) * t110 + (-t196 + t52) * t107, -t191 * qJD(3) - t86 * t102 - t134 * t205 + t43 * t135 - t52 * t73 - t53 * t74, t263, t260, t254, -t256, 0, t239 * t164 - t252 * t68 - t259, -t240 * t164 + t61 * t68 - t258, t276, -t280, -t275, t271, -t267, t71 * t11 + t182 * t172 - t251 * t175 + t239 * t44 + t265, t71 * t10 + t251 * t172 + t182 * t175 - t190 * t239 + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, 0, -t167 * t180 - t179, -t139 * qJD(2) + t121 * t224 + t149, -t107 * t224 - t169 * t179, -t110 * t224 - t170 * t179, -t169 * t100 - t170 * t101 + (t107 * t170 - t110 * t169) * qJD(2), t191 * qJD(2) + t43 * t169 - t170 * t205 - t86 * t224, 0, 0, 0, 0, 0, -t228 * t164 - t224 * t252, t229 * t164 + t224 * t61, 0, 0, 0, 0, 0, -t126 * t235 + t123 * t11 - t228 * t44 + (-t229 * t172 - t187 * t175) * t261, -t126 * t24 + t123 * t10 + t228 * t190 + (t187 * t172 - t229 * t175) * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142 + (-t110 + t211) * qJD(2) (t124 * qJD(1) + t107) * qJD(2), -t107 ^ 2 - t110 ^ 2, t53 * t107 + t52 * t110 + t77, 0, 0, 0, 0, 0, t26 + t264, -t185 + t262, 0, 0, 0, 0, 0, t268 + t277, -t175 * t204 - t235 - t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, -t260, -t254, t256, 0, -t15 * t164 + t259, -t14 * t164 + t258, -t276, t280, t275, -t271, t267, -pkin(5) * t11 - t15 * t44 + t183 * t172 - t250 * t175 - t265, -pkin(5) * t10 + t15 * t190 + t250 * t172 + t183 * t175 - t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190 * t44, t190 ^ 2 - t44 ^ 2, t10 + t244, -t11 - t243, t26, -t172 * t1 + t12 * t190 - t272 * t4 + t5, -t175 * t1 + t12 * t44 - t172 * t6 + t192 * t272;];
tauc_reg  = t3;
