% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:32
% EndTime: 2022-01-20 10:48:39
% DurationCPUTime: 2.20s
% Computational Cost: add. (3706->310), mult. (6062->388), div. (0->0), fcn. (3824->16), ass. (0->209)
t156 = qJDD(1) + qJDD(2);
t163 = sin(pkin(9));
t164 = cos(pkin(9));
t171 = cos(qJ(2));
t262 = pkin(1) * qJD(2);
t223 = qJD(1) * t262;
t167 = sin(qJ(2));
t230 = qJDD(1) * t167;
t282 = pkin(1) * t230 + t171 * t223;
t276 = pkin(1) * t171;
t143 = qJDD(1) * t276;
t77 = pkin(2) * t156 - t167 * t223 + t143;
t53 = t163 * t77 + t282 * t164;
t49 = pkin(7) * t156 + t53;
t284 = qJD(3) * qJD(4) + t49;
t162 = qJ(1) + qJ(2);
t146 = pkin(9) + t162;
t131 = cos(t146);
t271 = g(2) * t131;
t225 = t282 * t163 - t164 * t77;
t48 = -pkin(3) * t156 + t225;
t283 = t48 + t271;
t169 = cos(qJ(5));
t170 = cos(qJ(4));
t233 = qJD(5) * t169;
t235 = qJD(4) * t170;
t281 = -t169 * t235 - t170 * t233;
t158 = qJD(1) + qJD(2);
t165 = sin(qJ(5));
t166 = sin(qJ(4));
t243 = t166 * t169;
t94 = t165 * t170 + t243;
t76 = t94 * t158;
t224 = t166 * qJDD(3) + t284 * t170;
t236 = qJD(4) * t166;
t263 = pkin(1) * qJD(1);
t227 = t171 * t263;
t102 = pkin(2) * t158 + t227;
t228 = t167 * t263;
t69 = t163 * t102 + t164 * t228;
t65 = pkin(7) * t158 + t69;
t20 = -t65 * t236 + t224;
t145 = t170 * qJDD(3);
t232 = t166 * qJD(3);
t254 = t170 * t65;
t55 = t232 + t254;
t21 = -qJD(4) * t55 - t166 * t49 + t145;
t147 = t170 * qJD(3);
t256 = t166 * t65;
t54 = t147 - t256;
t179 = -t21 * t166 + t20 * t170 + (-t166 * t55 - t170 * t54) * qJD(4);
t130 = sin(t146);
t280 = g(1) * t131 + g(2) * t130;
t149 = sin(t162);
t151 = cos(t162);
t279 = g(1) * t149 - g(2) * t151;
t157 = qJD(4) + qJD(5);
t215 = pkin(8) * t158 + t65;
t51 = t170 * t215 + t232;
t13 = qJDD(4) * pkin(4) + t145 + (-pkin(8) * t156 - t49) * t166 - t51 * qJD(4);
t218 = t158 * t236;
t241 = t170 * t156;
t188 = t218 - t241;
t14 = -pkin(8) * t188 + t20;
t234 = qJD(5) * t165;
t50 = -t166 * t215 + t147;
t45 = qJD(4) * pkin(4) + t50;
t3 = (qJD(5) * t45 + t14) * t169 + t13 * t165 - t51 * t234;
t142 = pkin(2) + t276;
t246 = t164 * t167;
t87 = pkin(1) * t246 + t163 * t142;
t81 = pkin(7) + t87;
t278 = -pkin(8) - t81;
t168 = sin(qJ(1));
t277 = pkin(1) * t168;
t275 = pkin(2) * t149;
t274 = pkin(2) * t163;
t273 = pkin(2) * t164;
t272 = pkin(4) * t170;
t124 = g(1) * t130;
t269 = g(3) * t170;
t245 = t165 * t166;
t221 = t158 * t245;
t242 = t169 * t170;
t74 = -t158 * t242 + t221;
t268 = t76 * t74;
t132 = pkin(7) + t274;
t267 = -pkin(8) - t132;
t91 = t267 * t166;
t152 = t170 * pkin(8);
t92 = t132 * t170 + t152;
t57 = -t165 * t92 + t169 * t91;
t119 = t163 * t228;
t84 = t164 * t227 - t119;
t211 = qJD(4) * t267;
t88 = t166 * t211;
t89 = t170 * t211;
t93 = -t242 + t245;
t266 = qJD(5) * t57 + t165 * t89 + t169 * t88 + t93 * t84;
t58 = t165 * t91 + t169 * t92;
t265 = -qJD(5) * t58 - t165 * t88 + t169 * t89 + t94 * t84;
t244 = t166 * t156;
t198 = t165 * t244 - t169 * t241;
t62 = t157 * t94;
t33 = t158 * t62 + t198;
t196 = t157 * t245;
t61 = t196 + t281;
t264 = -t94 * t33 + t61 * t74;
t187 = pkin(1) * (t163 * t171 + t246);
t82 = qJD(1) * t187;
t261 = t158 * t82;
t83 = qJD(2) * t187;
t260 = t158 * t83;
t259 = t158 * t84;
t247 = t163 * t167;
t85 = (t164 * t171 - t247) * t262;
t258 = t158 * t85;
t257 = t165 * t51;
t255 = t169 * t51;
t68 = t102 * t164 - t119;
t64 = -pkin(3) * t158 - t68;
t253 = t170 * t124 + t64 * t236;
t161 = qJ(4) + qJ(5);
t148 = sin(t161);
t252 = t130 * t148;
t150 = cos(t161);
t251 = t130 * t150;
t250 = t131 * t148;
t249 = t131 * t150;
t248 = t158 * t166;
t239 = g(1) * t151 + g(2) * t149;
t159 = t166 ^ 2;
t160 = t170 ^ 2;
t238 = t159 - t160;
t237 = t159 + t160;
t229 = t283 * t166 + t64 * t235;
t226 = pkin(4) * t236;
t154 = t158 ^ 2;
t220 = t166 * t154 * t170;
t139 = pkin(2) * t151;
t219 = t131 * pkin(3) + t130 * pkin(7) + t139;
t141 = pkin(3) + t272;
t213 = qJD(4) * t278;
t209 = t237 * t156;
t86 = -pkin(1) * t247 + t142 * t164;
t208 = qJD(1) * (-qJD(2) + t158);
t207 = qJD(2) * (-qJD(1) - t158);
t205 = -t156 * t243 + t281 * t158 - t165 * t241;
t204 = t170 * t218;
t203 = -t82 + t226;
t80 = -pkin(3) - t86;
t202 = t143 + t279;
t173 = -pkin(8) - pkin(7);
t201 = -t130 * t173 + t131 * t141 + t139;
t172 = cos(qJ(1));
t199 = g(1) * t168 - g(2) * t172;
t32 = t158 * t196 + t205;
t197 = -t32 * t93 + t62 * t76;
t195 = -t53 + t280;
t15 = t169 * t45 - t257;
t16 = t165 * t45 + t255;
t4 = -qJD(5) * t16 + t169 * t13 - t14 * t165;
t194 = t15 * t61 - t16 * t62 - t3 * t93 - t4 * t94 - t280;
t155 = qJDD(4) + qJDD(5);
t37 = t155 * t94 - t157 * t61;
t66 = t278 * t166;
t67 = t170 * t81 + t152;
t35 = -t165 * t67 + t169 * t66;
t36 = t165 * t66 + t169 * t67;
t193 = t166 * t54 - t170 * t55;
t192 = -pkin(3) * t130 + t131 * pkin(7) - t275;
t31 = pkin(4) * t188 + t48;
t56 = -t141 * t158 - t68;
t191 = -g(1) * t252 + g(2) * t250 + t31 * t94 - t56 * t61;
t190 = g(1) * t251 - g(2) * t249 + t31 * t93 + t56 * t62;
t189 = t124 - t225 - t271;
t186 = -t158 * t64 + t280;
t174 = qJD(4) ^ 2;
t185 = t156 * t80 + t174 * t81 + t260;
t184 = -t130 * t141 - t131 * t173 - t275;
t133 = -pkin(3) - t273;
t183 = t132 * t174 + t133 * t156 - t261;
t182 = -qJDD(4) * t81 + (t158 * t80 - t85) * qJD(4);
t181 = -qJDD(4) * t132 + (t133 * t158 + t84) * qJD(4);
t178 = -t280 + t179;
t177 = g(1) * t249 + g(2) * t251 + g(3) * t148 + t56 * t74 - t3;
t176 = g(1) * t250 + g(2) * t252 - g(3) * t150 - t56 * t76 + t4;
t153 = t172 * pkin(1);
t107 = qJDD(4) * t170 - t166 * t174;
t106 = qJDD(4) * t166 + t170 * t174;
t105 = -t141 - t273;
t79 = t156 * t160 - 0.2e1 * t204;
t78 = t156 * t159 + 0.2e1 * t204;
t73 = t80 - t272;
t70 = t83 + t226;
t63 = -0.2e1 * t238 * t158 * qJD(4) + 0.2e1 * t166 * t241;
t43 = -t166 * t85 + t170 * t213;
t42 = t166 * t213 + t170 * t85;
t38 = -t155 * t93 - t157 * t62;
t34 = -t74 ^ 2 + t76 ^ 2;
t22 = -t205 + (-t221 + t74) * t157;
t19 = t169 * t50 - t257;
t18 = -t165 * t50 - t255;
t11 = t33 * t93 + t62 * t74;
t10 = -t32 * t94 - t61 * t76;
t7 = -qJD(5) * t36 - t165 * t42 + t169 * t43;
t6 = qJD(5) * t35 + t165 * t43 + t169 * t42;
t5 = -t197 + t264;
t1 = [0, 0, 0, 0, 0, qJDD(1), t199, g(1) * t172 + g(2) * t168, 0, 0, 0, 0, 0, 0, 0, t156, (t156 * t171 + t167 * t207) * pkin(1) + t202, ((-qJDD(1) - t156) * t167 + t171 * t207) * pkin(1) + t239, 0, (t199 + (t167 ^ 2 + t171 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t156, t156 * t86 + t189 - t260, -t156 * t87 + t195 - t258, 0, t53 * t87 + t69 * t85 - t225 * t86 - t68 * t83 - g(1) * (-t275 - t277) - g(2) * (t139 + t153), t78, t63, t106, t79, t107, 0, t182 * t166 + (-t185 - t283) * t170 + t253, t182 * t170 + (t185 - t124) * t166 + t229, t81 * t209 + t237 * t258 + t178, t48 * t80 + t64 * t83 - g(1) * (t192 - t277) - g(2) * (t153 + t219) - t193 * t85 + t179 * t81, t10, t5, t37, t11, t38, 0, t155 * t35 + t157 * t7 + t33 * t73 + t70 * t74 + t190, -t155 * t36 - t157 * t6 - t32 * t73 + t70 * t76 + t191, t32 * t35 - t33 * t36 - t6 * t74 - t7 * t76 + t194, t3 * t36 + t16 * t6 + t4 * t35 + t15 * t7 + t31 * t73 + t56 * t70 - g(1) * (t184 - t277) - g(2) * (t153 + t201); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, pkin(1) * t167 * t208 + t202, (t171 * t208 - t230) * pkin(1) + t239, 0, 0, 0, 0, 0, 0, 0, t156, t156 * t273 + t189 + t261, -t156 * t274 + t195 + t259, 0, t68 * t82 - t69 * t84 + (t163 * t53 - t164 * t225 + t279) * pkin(2), t78, t63, t106, t79, t107, 0, t181 * t166 + (-t183 - t283) * t170 + t253, t181 * t170 + (t183 - t124) * t166 + t229, t132 * t209 - t237 * t259 + t178, -g(1) * t192 - g(2) * t219 + t132 * t179 + t48 * t133 + t193 * t84 - t64 * t82, t10, t5, t37, t11, t38, 0, t105 * t33 + t155 * t57 + t157 * t265 + t203 * t74 + t190, -t105 * t32 - t155 * t58 - t157 * t266 + t203 * t76 + t191, -t265 * t76 - t266 * t74 + t32 * t57 - t33 * t58 + t194, -g(1) * t184 - g(2) * t201 + t31 * t105 + t15 * t265 + t16 * t266 + t203 * t56 + t3 * t58 + t4 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0, 0, 0, 0, t107, -t106, 0, -qJD(4) * t193 + t166 * t20 + t170 * t21 - g(3), 0, 0, 0, 0, 0, 0, t38, -t37, t197 + t264, -t15 * t62 - t16 * t61 + t3 * t94 - t4 * t93 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t238 * t154, t244, t220, t241, qJDD(4), -t269 + t145 + (t55 - t254) * qJD(4) + (t186 - t284) * t166, g(3) * t166 + (t54 + t256) * qJD(4) + t186 * t170 - t224, 0, 0, t268, t34, t22, -t268, -t198, t155, -t157 * t18 + (t155 * t169 - t157 * t234 - t248 * t74) * pkin(4) + t176, t157 * t19 + (-t155 * t165 - t157 * t233 - t248 * t76) * pkin(4) + t177, (t16 + t18) * t76 + (-t15 + t19) * t74 + (-t165 * t33 + t169 * t32 + (t165 * t76 - t169 * t74) * qJD(5)) * pkin(4), -t15 * t18 - t16 * t19 + (-t269 + t165 * t3 + t169 * t4 + (-t15 * t165 + t16 * t169) * qJD(5) + (-t158 * t56 + t280) * t166) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t34, t22, -t268, -t198, t155, t157 * t16 + t176, t15 * t157 + t177, 0, 0;];
tau_reg = t1;
