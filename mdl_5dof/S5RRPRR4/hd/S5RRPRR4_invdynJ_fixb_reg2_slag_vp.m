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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:02:11
% EndTime: 2020-01-03 12:02:14
% DurationCPUTime: 1.76s
% Computational Cost: add. (3706->306), mult. (6062->383), div. (0->0), fcn. (3824->16), ass. (0->206)
t161 = qJ(1) + qJ(2);
t144 = pkin(9) + t161;
t129 = sin(t144);
t130 = cos(t144);
t200 = -g(2) * t130 - g(3) * t129;
t155 = qJDD(1) + qJDD(2);
t162 = sin(pkin(9));
t163 = cos(pkin(9));
t170 = cos(qJ(2));
t260 = pkin(1) * qJD(2);
t225 = qJD(1) * t260;
t166 = sin(qJ(2));
t232 = qJDD(1) * t166;
t277 = pkin(1) * t232 + t170 * t225;
t267 = t170 * pkin(1);
t141 = qJDD(1) * t267;
t77 = pkin(2) * t155 - t166 * t225 + t141;
t228 = t162 * t277 - t163 * t77;
t48 = -t155 * pkin(3) + t228;
t190 = t200 - t48;
t53 = t162 * t77 + t163 * t277;
t49 = pkin(7) * t155 + t53;
t278 = qJD(3) * qJD(4) + t49;
t275 = g(2) * t129 - g(3) * t130;
t168 = cos(qJ(5));
t169 = cos(qJ(4));
t235 = qJD(5) * t168;
t237 = qJD(4) * t169;
t276 = -t168 * t237 - t169 * t235;
t157 = qJD(1) + qJD(2);
t164 = sin(qJ(5));
t165 = sin(qJ(4));
t244 = t168 * t165;
t94 = t164 * t169 + t244;
t76 = t94 * t157;
t227 = t165 * qJDD(3) + t169 * t278;
t238 = qJD(4) * t165;
t261 = pkin(1) * qJD(1);
t229 = t170 * t261;
t101 = pkin(2) * t157 + t229;
t231 = t166 * t261;
t69 = t101 * t162 + t163 * t231;
t65 = pkin(7) * t157 + t69;
t20 = -t238 * t65 + t227;
t143 = t169 * qJDD(3);
t234 = t165 * qJD(3);
t256 = t169 * t65;
t55 = t234 + t256;
t21 = -t55 * qJD(4) - t165 * t49 + t143;
t145 = t169 * qJD(3);
t258 = t165 * t65;
t54 = t145 - t258;
t178 = t20 * t169 + (-t165 * t55 - t169 * t54) * qJD(4) - t21 * t165;
t156 = qJD(4) + qJD(5);
t217 = pkin(8) * t157 + t65;
t51 = t169 * t217 + t234;
t13 = qJDD(4) * pkin(4) + t143 + (-pkin(8) * t155 - t49) * t165 - t51 * qJD(4);
t220 = t157 * t238;
t242 = t169 * t155;
t188 = t220 - t242;
t14 = -pkin(8) * t188 + t20;
t236 = qJD(5) * t164;
t50 = -t165 * t217 + t145;
t45 = qJD(4) * pkin(4) + t50;
t3 = (qJD(5) * t45 + t14) * t168 + t164 * t13 - t51 * t236;
t140 = pkin(2) + t267;
t247 = t163 * t166;
t87 = pkin(1) * t247 + t140 * t162;
t81 = pkin(7) + t87;
t274 = -pkin(8) - t81;
t273 = g(1) * t169;
t270 = t162 * pkin(2);
t269 = t163 * pkin(2);
t268 = t169 * pkin(4);
t246 = t164 * t165;
t224 = t157 * t246;
t243 = t168 * t169;
t74 = -t157 * t243 + t224;
t266 = t76 * t74;
t131 = pkin(7) + t270;
t265 = -pkin(8) - t131;
t91 = t265 * t165;
t151 = t169 * pkin(8);
t92 = t131 * t169 + t151;
t57 = -t164 * t92 + t168 * t91;
t119 = t162 * t231;
t84 = t163 * t229 - t119;
t213 = qJD(4) * t265;
t88 = t165 * t213;
t89 = t169 * t213;
t93 = -t243 + t246;
t264 = qJD(5) * t57 + t164 * t89 + t168 * t88 + t84 * t93;
t58 = t164 * t91 + t168 * t92;
t263 = -qJD(5) * t58 - t164 * t88 + t168 * t89 + t84 * t94;
t245 = t165 * t155;
t196 = t164 * t245 - t168 * t242;
t62 = t156 * t94;
t33 = t157 * t62 + t196;
t194 = t156 * t246;
t61 = t194 + t276;
t262 = -t33 * t94 + t61 * t74;
t259 = t164 * t51;
t257 = t168 * t51;
t187 = pkin(1) * (t162 * t170 + t247);
t82 = qJD(1) * t187;
t255 = t82 * t157;
t83 = qJD(2) * t187;
t254 = t83 * t157;
t253 = t84 * t157;
t248 = t162 * t166;
t85 = (t163 * t170 - t248) * t260;
t252 = t85 * t157;
t160 = qJ(4) + qJ(5);
t146 = sin(t160);
t251 = t129 * t146;
t250 = t130 * t146;
t249 = t157 * t165;
t158 = t165 ^ 2;
t159 = t169 ^ 2;
t240 = t158 - t159;
t239 = t158 + t159;
t230 = pkin(4) * t238;
t147 = sin(t161);
t136 = pkin(2) * t147;
t139 = pkin(3) + t268;
t172 = -pkin(8) - pkin(7);
t226 = t129 * t139 + t130 * t172 + t136;
t153 = t157 ^ 2;
t223 = t165 * t153 * t169;
t149 = cos(t161);
t137 = pkin(2) * t149;
t221 = pkin(3) * t130 + pkin(7) * t129 + t137;
t216 = qJD(4) * t274;
t215 = g(2) * t147 - g(3) * t149;
t211 = t239 * t155;
t31 = pkin(4) * t188 + t48;
t68 = t163 * t101 - t119;
t56 = -t139 * t157 - t68;
t210 = g(2) * t250 + g(3) * t251 + t31 * t94 - t56 * t61;
t86 = -pkin(1) * t248 + t140 * t163;
t64 = -pkin(3) * t157 - t68;
t209 = -t165 * t190 + t64 * t237;
t208 = qJD(1) * (-qJD(2) + t157);
t207 = qJD(2) * (-qJD(1) - t157);
t205 = -t155 * t244 + t157 * t276 - t164 * t242;
t204 = t169 * t220;
t203 = -t82 + t230;
t80 = -pkin(3) - t86;
t202 = pkin(3) * t129 - pkin(7) * t130 + t136;
t201 = -t129 * t172 + t130 * t139 + t137;
t198 = -g(2) * t149 - g(3) * t147;
t167 = sin(qJ(1));
t171 = cos(qJ(1));
t197 = -g(2) * t171 - g(3) * t167;
t32 = t157 * t194 + t205;
t195 = -t32 * t93 + t62 * t76;
t193 = -t53 + t275;
t15 = t168 * t45 - t259;
t16 = t164 * t45 + t257;
t4 = -qJD(5) * t16 + t13 * t168 - t164 * t14;
t192 = t15 * t61 - t16 * t62 - t3 * t93 - t4 * t94 - t275;
t154 = qJDD(4) + qJDD(5);
t37 = t154 * t94 - t156 * t61;
t66 = t274 * t165;
t67 = t169 * t81 + t151;
t35 = -t164 * t67 + t168 * t66;
t36 = t164 * t66 + t168 * t67;
t191 = t165 * t54 - t169 * t55;
t189 = t141 + t198;
t186 = -t157 * t64 + t275;
t173 = qJD(4) ^ 2;
t185 = t155 * t80 + t173 * t81 + t254;
t184 = t200 - t228;
t132 = -pkin(3) - t269;
t183 = t131 * t173 + t132 * t155 - t255;
t182 = -qJDD(4) * t81 + (t157 * t80 - t85) * qJD(4);
t148 = cos(t160);
t181 = t148 * t200 + t31 * t93 + t56 * t62;
t180 = -qJDD(4) * t131 + (t132 * t157 + t84) * qJD(4);
t177 = -t275 + t178;
t176 = g(1) * t146 + t148 * t275 + t56 * t74 - t3;
t175 = -g(1) * t148 + g(2) * t251 - g(3) * t250 - t56 * t76 + t4;
t152 = t171 * pkin(1);
t150 = t167 * pkin(1);
t106 = qJDD(4) * t169 - t165 * t173;
t105 = qJDD(4) * t165 + t169 * t173;
t104 = -t139 - t269;
t79 = t155 * t159 - 0.2e1 * t204;
t78 = t155 * t158 + 0.2e1 * t204;
t73 = t80 - t268;
t70 = t83 + t230;
t63 = -0.2e1 * qJD(4) * t157 * t240 + 0.2e1 * t165 * t242;
t59 = t64 * t238;
t43 = -t165 * t85 + t169 * t216;
t42 = t165 * t216 + t169 * t85;
t38 = -t154 * t93 - t156 * t62;
t34 = -t74 ^ 2 + t76 ^ 2;
t22 = -t205 + (-t224 + t74) * t156;
t19 = t168 * t50 - t259;
t18 = -t164 * t50 - t257;
t11 = t33 * t93 + t62 * t74;
t10 = -t32 * t94 - t61 * t76;
t7 = -qJD(5) * t36 - t164 * t42 + t168 * t43;
t6 = qJD(5) * t35 + t164 * t43 + t168 * t42;
t5 = -t195 + t262;
t1 = [0, 0, 0, 0, 0, qJDD(1), t197, g(2) * t167 - g(3) * t171, 0, 0, 0, 0, 0, 0, 0, t155, (t155 * t170 + t166 * t207) * pkin(1) + t189, ((-qJDD(1) - t155) * t166 + t170 * t207) * pkin(1) + t215, 0, (t197 + (t166 ^ 2 + t170 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t155, t155 * t86 + t184 - t254, -t155 * t87 + t193 - t252, 0, t53 * t87 + t69 * t85 - t228 * t86 - t68 * t83 - g(2) * (t137 + t152) - g(3) * (t136 + t150), t78, t63, t105, t79, t106, 0, t59 + t182 * t165 + (-t185 + t190) * t169, t165 * t185 + t169 * t182 + t209, t211 * t81 + t239 * t252 + t177, t48 * t80 + t64 * t83 - g(2) * (t152 + t221) - g(3) * (t150 + t202) - t191 * t85 + t178 * t81, t10, t5, t37, t11, t38, 0, t35 * t154 + t7 * t156 + t73 * t33 + t70 * t74 + t181, -t154 * t36 - t156 * t6 - t32 * t73 + t70 * t76 + t210, t32 * t35 - t33 * t36 - t6 * t74 - t7 * t76 + t192, t3 * t36 + t16 * t6 + t4 * t35 + t15 * t7 + t31 * t73 + t56 * t70 - g(2) * (t152 + t201) - g(3) * (t150 + t226); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, pkin(1) * t166 * t208 + t189, (t170 * t208 - t232) * pkin(1) + t215, 0, 0, 0, 0, 0, 0, 0, t155, t155 * t269 + t184 + t255, -t155 * t270 + t193 + t253, 0, t68 * t82 - t69 * t84 + (t162 * t53 - t163 * t228 + t198) * pkin(2), t78, t63, t105, t79, t106, 0, t59 + t180 * t165 + (-t183 + t190) * t169, t165 * t183 + t169 * t180 + t209, t131 * t211 - t239 * t253 + t177, -g(2) * t221 - g(3) * t202 + t131 * t178 + t48 * t132 + t191 * t84 - t64 * t82, t10, t5, t37, t11, t38, 0, t104 * t33 + t57 * t154 + t156 * t263 + t203 * t74 + t181, -t104 * t32 - t58 * t154 - t156 * t264 + t203 * t76 + t210, -t263 * t76 - t264 * t74 + t57 * t32 - t58 * t33 + t192, -g(2) * t201 - g(3) * t226 + t31 * t104 + t15 * t263 + t16 * t264 + t203 * t56 + t3 * t58 + t4 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, t106, -t105, 0, -qJD(4) * t191 + t20 * t165 + t21 * t169 - g(1), 0, 0, 0, 0, 0, 0, t38, -t37, t195 + t262, -t15 * t62 - t16 * t61 + t3 * t94 - t4 * t93 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t240 * t153, t245, t223, t242, qJDD(4), -t273 + t143 + (t55 - t256) * qJD(4) + (t186 - t278) * t165, g(1) * t165 + (t54 + t258) * qJD(4) + t186 * t169 - t227, 0, 0, t266, t34, t22, -t266, -t196, t154, -t18 * t156 + (t154 * t168 - t156 * t236 - t249 * t74) * pkin(4) + t175, t19 * t156 + (-t154 * t164 - t156 * t235 - t249 * t76) * pkin(4) + t176, (t16 + t18) * t76 + (-t15 + t19) * t74 + (-t164 * t33 + t168 * t32 + (t164 * t76 - t168 * t74) * qJD(5)) * pkin(4), -t15 * t18 - t16 * t19 + (-t273 + t164 * t3 + t168 * t4 + (-t15 * t164 + t16 * t168) * qJD(5) + (-t157 * t56 + t275) * t165) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, t34, t22, -t266, -t196, t154, t16 * t156 + t175, t15 * t156 + t176, 0, 0;];
tau_reg = t1;
