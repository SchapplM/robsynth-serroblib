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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:32:13
% EndTime: 2019-12-05 18:32:17
% DurationCPUTime: 2.11s
% Computational Cost: add. (3706->308), mult. (6062->384), div. (0->0), fcn. (3824->16), ass. (0->205)
t151 = qJDD(1) + qJDD(2);
t158 = sin(pkin(9));
t159 = cos(pkin(9));
t166 = cos(qJ(2));
t254 = pkin(1) * qJD(2);
t219 = qJD(1) * t254;
t162 = sin(qJ(2));
t226 = qJDD(1) * t162;
t274 = pkin(1) * t226 + t166 * t219;
t263 = t166 * pkin(1);
t139 = qJDD(1) * t263;
t77 = t151 * pkin(2) - t162 * t219 + t139;
t53 = t158 * t77 + t274 * t159;
t49 = t151 * pkin(7) + t53;
t276 = qJD(3) * qJD(4) + t49;
t157 = qJ(1) + qJ(2);
t142 = pkin(9) + t157;
t127 = sin(t142);
t128 = cos(t142);
t275 = g(2) * t128 + g(3) * t127;
t235 = g(2) * t127 - g(3) * t128;
t164 = cos(qJ(5));
t165 = cos(qJ(4));
t229 = qJD(5) * t164;
t231 = qJD(4) * t165;
t273 = -t164 * t231 - t165 * t229;
t153 = qJD(1) + qJD(2);
t160 = sin(qJ(5));
t161 = sin(qJ(4));
t238 = t164 * t161;
t94 = t160 * t165 + t238;
t76 = t94 * t153;
t220 = t161 * qJDD(3) + t276 * t165;
t232 = qJD(4) * t161;
t255 = pkin(1) * qJD(1);
t225 = t162 * t255;
t223 = t166 * t255;
t99 = t153 * pkin(2) + t223;
t69 = t158 * t99 + t159 * t225;
t65 = t153 * pkin(7) + t69;
t20 = -t65 * t232 + t220;
t141 = t165 * qJDD(3);
t228 = t161 * qJD(3);
t250 = t165 * t65;
t55 = t228 + t250;
t21 = -t55 * qJD(4) - t161 * t49 + t141;
t143 = t165 * qJD(3);
t252 = t161 * t65;
t54 = t143 - t252;
t174 = -t21 * t161 + t20 * t165 + (-t161 * t55 - t165 * t54) * qJD(4);
t145 = sin(t157);
t147 = cos(t157);
t272 = g(2) * t147 + g(3) * t145;
t152 = qJD(4) + qJD(5);
t211 = pkin(8) * t153 + t65;
t51 = t211 * t165 + t228;
t13 = qJDD(4) * pkin(4) + t141 + (-pkin(8) * t151 - t49) * t161 - t51 * qJD(4);
t214 = t153 * t232;
t236 = t165 * t151;
t185 = t214 - t236;
t14 = -t185 * pkin(8) + t20;
t230 = qJD(5) * t160;
t50 = -t211 * t161 + t143;
t45 = qJD(4) * pkin(4) + t50;
t3 = (qJD(5) * t45 + t14) * t164 + t160 * t13 - t51 * t230;
t138 = pkin(2) + t263;
t241 = t159 * t162;
t87 = pkin(1) * t241 + t158 * t138;
t81 = pkin(7) + t87;
t271 = -pkin(8) - t81;
t270 = pkin(2) * t145;
t269 = pkin(2) * t147;
t268 = g(1) * t165;
t267 = t158 * pkin(2);
t266 = t159 * pkin(2);
t163 = sin(qJ(1));
t265 = t163 * pkin(1);
t264 = t165 * pkin(4);
t167 = cos(qJ(1));
t262 = t167 * pkin(1);
t240 = t160 * t161;
t218 = t153 * t240;
t237 = t164 * t165;
t74 = -t153 * t237 + t218;
t261 = t76 * t74;
t129 = pkin(7) + t267;
t260 = -pkin(8) - t129;
t91 = t260 * t161;
t148 = t165 * pkin(8);
t92 = t165 * t129 + t148;
t57 = -t160 * t92 + t164 * t91;
t117 = t158 * t225;
t84 = t159 * t223 - t117;
t207 = qJD(4) * t260;
t88 = t161 * t207;
t89 = t165 * t207;
t93 = -t237 + t240;
t259 = t57 * qJD(5) + t160 * t89 + t164 * t88 + t93 * t84;
t58 = t160 * t91 + t164 * t92;
t258 = -t58 * qJD(5) - t160 * t88 + t164 * t89 + t94 * t84;
t221 = t274 * t158 - t159 * t77;
t48 = -t151 * pkin(3) + t221;
t68 = t159 * t99 - t117;
t64 = -t153 * pkin(3) - t68;
t257 = t48 * t161 + t64 * t231;
t239 = t161 * t151;
t194 = t160 * t239 - t164 * t236;
t62 = t152 * t94;
t33 = t62 * t153 + t194;
t192 = t152 * t240;
t61 = t192 + t273;
t256 = -t94 * t33 + t61 * t74;
t253 = t160 * t51;
t251 = t164 * t51;
t184 = pkin(1) * (t158 * t166 + t241);
t82 = qJD(1) * t184;
t249 = t82 * t153;
t83 = qJD(2) * t184;
t248 = t83 * t153;
t247 = t84 * t153;
t242 = t158 * t162;
t85 = (t159 * t166 - t242) * t254;
t246 = t85 * t153;
t156 = qJ(4) + qJ(5);
t146 = cos(t156);
t245 = t127 * t146;
t244 = t128 * t146;
t243 = t153 * t161;
t154 = t161 ^ 2;
t155 = t165 ^ 2;
t234 = t154 - t155;
t233 = t154 + t155;
t224 = pkin(4) * t232;
t222 = t275 * t165 + t64 * t232;
t149 = t153 ^ 2;
t217 = t161 * t149 * t165;
t215 = t139 + t272;
t137 = pkin(3) + t264;
t210 = qJD(4) * t271;
t209 = -g(2) * t145 + g(3) * t147;
t31 = t185 * pkin(4) + t48;
t56 = -t137 * t153 - t68;
t206 = g(2) * t244 + g(3) * t245 + t31 * t93 + t56 * t62;
t204 = t233 * t151;
t86 = -pkin(1) * t242 + t159 * t138;
t203 = qJD(1) * (-qJD(2) + t153);
t202 = qJD(2) * (-qJD(1) - t153);
t200 = -t151 * t238 + t273 * t153 - t160 * t236;
t199 = t165 * t214;
t198 = -t82 + t224;
t80 = -pkin(3) - t86;
t195 = g(2) * t167 + g(3) * t163;
t32 = t153 * t192 + t200;
t193 = -t93 * t32 + t76 * t62;
t191 = -t221 + t275;
t190 = -t53 - t235;
t15 = t164 * t45 - t253;
t16 = t160 * t45 + t251;
t4 = -t16 * qJD(5) + t164 * t13 - t160 * t14;
t189 = t15 * t61 - t16 * t62 - t3 * t93 - t4 * t94 + t235;
t150 = qJDD(4) + qJDD(5);
t37 = t94 * t150 - t61 * t152;
t66 = t271 * t161;
t67 = t165 * t81 + t148;
t35 = -t160 * t67 + t164 * t66;
t36 = t160 * t66 + t164 * t67;
t188 = t54 * t161 - t55 * t165;
t187 = -t127 * pkin(3) + t128 * pkin(7) - t270;
t168 = -pkin(8) - pkin(7);
t186 = t127 * t168 - t128 * t137 - t269;
t183 = -t128 * pkin(3) - t127 * pkin(7) - t269;
t182 = -t153 * t64 - t235;
t169 = qJD(4) ^ 2;
t181 = -t151 * t80 - t169 * t81 - t248;
t180 = -t127 * t137 - t128 * t168 - t270;
t130 = -pkin(3) - t266;
t179 = -t129 * t169 - t130 * t151 + t249;
t178 = -qJDD(4) * t81 + (t153 * t80 - t85) * qJD(4);
t144 = sin(t156);
t177 = -t144 * t275 + t31 * t94 - t56 * t61;
t176 = -qJDD(4) * t129 + (t130 * t153 + t84) * qJD(4);
t173 = t235 + t174;
t172 = g(1) * t144 - g(2) * t245 + g(3) * t244 + t56 * t74 - t3;
t171 = -g(1) * t146 - t235 * t144 - t56 * t76 + t4;
t104 = qJDD(4) * t165 - t169 * t161;
t103 = qJDD(4) * t161 + t169 * t165;
t102 = -t137 - t266;
t79 = t155 * t151 - 0.2e1 * t199;
t78 = t154 * t151 + 0.2e1 * t199;
t73 = t80 - t264;
t70 = t83 + t224;
t63 = -0.2e1 * t234 * t153 * qJD(4) + 0.2e1 * t161 * t236;
t43 = -t161 * t85 + t165 * t210;
t42 = t161 * t210 + t165 * t85;
t38 = -t93 * t150 - t62 * t152;
t34 = -t74 ^ 2 + t76 ^ 2;
t22 = -t200 + (-t218 + t74) * t152;
t19 = t164 * t50 - t253;
t18 = -t160 * t50 - t251;
t11 = t33 * t93 + t74 * t62;
t10 = -t32 * t94 - t76 * t61;
t7 = -t36 * qJD(5) - t160 * t42 + t164 * t43;
t6 = t35 * qJD(5) + t160 * t43 + t164 * t42;
t5 = -t193 + t256;
t1 = [0, 0, 0, 0, 0, qJDD(1), t195, -g(2) * t163 + g(3) * t167, 0, 0, 0, 0, 0, 0, 0, t151, (t151 * t166 + t162 * t202) * pkin(1) + t215, ((-qJDD(1) - t151) * t162 + t166 * t202) * pkin(1) + t209, 0, (t195 + (t162 ^ 2 + t166 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t151, t86 * t151 + t191 - t248, -t87 * t151 + t190 - t246, 0, t53 * t87 + t69 * t85 - t221 * t86 - t68 * t83 - g(2) * (-t262 - t269) - g(3) * (-t265 - t270), t78, t63, t103, t79, t104, 0, t178 * t161 + (t181 - t48) * t165 + t222, t178 * t165 + (-t181 - t275) * t161 + t257, t204 * t81 + t233 * t246 + t173, t48 * t80 + t64 * t83 - g(2) * (t183 - t262) - g(3) * (t187 - t265) - t188 * t85 + t174 * t81, t10, t5, t37, t11, t38, 0, t35 * t150 + t7 * t152 + t73 * t33 + t70 * t74 + t206, -t36 * t150 - t6 * t152 - t73 * t32 + t70 * t76 + t177, t35 * t32 - t36 * t33 - t6 * t74 - t7 * t76 + t189, t3 * t36 + t16 * t6 + t4 * t35 + t15 * t7 + t31 * t73 + t56 * t70 - g(2) * (t186 - t262) - g(3) * (t180 - t265); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, pkin(1) * t162 * t203 + t215, (t166 * t203 - t226) * pkin(1) + t209, 0, 0, 0, 0, 0, 0, 0, t151, t151 * t266 + t191 + t249, -t151 * t267 + t190 + t247, 0, t68 * t82 - t69 * t84 + (t158 * t53 - t159 * t221 + t272) * pkin(2), t78, t63, t103, t79, t104, 0, t176 * t161 + (t179 - t48) * t165 + t222, t176 * t165 + (-t179 - t275) * t161 + t257, t129 * t204 - t233 * t247 + t173, -g(2) * t183 - g(3) * t187 + t129 * t174 + t48 * t130 + t188 * t84 - t64 * t82, t10, t5, t37, t11, t38, 0, t102 * t33 + t57 * t150 + t258 * t152 + t198 * t74 + t206, -t102 * t32 - t58 * t150 - t259 * t152 + t198 * t76 + t177, -t258 * t76 - t259 * t74 + t57 * t32 - t58 * t33 + t189, -g(2) * t186 - g(3) * t180 + t31 * t102 + t258 * t15 + t259 * t16 + t198 * t56 + t3 * t58 + t4 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, t104, -t103, 0, -qJD(4) * t188 + t20 * t161 + t21 * t165 - g(1), 0, 0, 0, 0, 0, 0, t38, -t37, t193 + t256, -t15 * t62 - t16 * t61 + t3 * t94 - t4 * t93 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t234 * t149, t239, t217, t236, qJDD(4), -t268 + t141 + (t55 - t250) * qJD(4) + (t182 - t276) * t161, g(1) * t161 + (t54 + t252) * qJD(4) + t182 * t165 - t220, 0, 0, t261, t34, t22, -t261, -t194, t150, -t18 * t152 + (t150 * t164 - t152 * t230 - t74 * t243) * pkin(4) + t171, t19 * t152 + (-t150 * t160 - t152 * t229 - t76 * t243) * pkin(4) + t172, (t16 + t18) * t76 + (-t15 + t19) * t74 + (-t160 * t33 + t164 * t32 + (t160 * t76 - t164 * t74) * qJD(5)) * pkin(4), -t15 * t18 - t16 * t19 + (-t268 + t160 * t3 + t164 * t4 + (-t15 * t160 + t16 * t164) * qJD(5) + (-t153 * t56 - t235) * t161) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t34, t22, -t261, -t194, t150, t16 * t152 + t171, t15 * t152 + t172, 0, 0;];
tau_reg = t1;
