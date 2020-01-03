% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:06:03
% DurationCPUTime: 3.29s
% Computational Cost: add. (2501->385), mult. (6291->469), div. (0->0), fcn. (3738->4), ass. (0->193)
t121 = sin(qJ(3));
t123 = cos(qJ(3));
t124 = cos(qJ(2));
t197 = qJD(1) * t124;
t107 = -qJD(3) + t197;
t122 = sin(qJ(2));
t198 = qJD(1) * t122;
t177 = t123 * t198;
t196 = qJD(2) * t121;
t86 = t177 + t196;
t213 = t86 * t107;
t192 = qJD(3) * t123;
t175 = t122 * t192;
t185 = qJD(2) * qJD(3);
t194 = qJD(2) * t124;
t54 = qJD(1) * (t121 * t194 + t175) + t121 * t185;
t264 = t54 - t213;
t178 = t121 * t198;
t188 = t123 * qJD(2);
t84 = t178 - t188;
t222 = t107 * t84;
t186 = qJD(1) * qJD(2);
t172 = t124 * t186;
t193 = qJD(3) * t121;
t176 = t122 * t193;
t53 = qJD(1) * t176 + (-t172 - t185) * t123;
t265 = t53 - t222;
t2 = t264 * t121 + t265 * t123;
t268 = t54 + t213;
t111 = t122 * t186;
t266 = qJ(4) * t111 - t107 * qJD(4);
t241 = pkin(3) + pkin(4);
t263 = t111 * t241;
t165 = pkin(6) * t178;
t92 = -pkin(2) * t124 - pkin(7) * t122 - pkin(1);
t73 = t92 * qJD(1);
t160 = pkin(2) * t122 - pkin(7) * t124;
t90 = t160 * qJD(2);
t76 = qJD(1) * t90;
t116 = pkin(6) * t197;
t98 = qJD(2) * pkin(7) + t116;
t169 = -qJD(2) * t165 - t123 * t76 + t98 * t192 + t73 * t193;
t41 = t121 * t73 + t123 * t98;
t142 = -t107 * t41 - t169;
t215 = t54 * t123;
t216 = t53 * t121;
t262 = (qJD(3) * (t121 * t84 - t123 * t86) - t215 + t216) * t122 - (t121 * t86 + t123 * t84) * t194;
t167 = -t86 + t196;
t174 = t107 * t192;
t206 = t123 * t124;
t24 = (t107 * t206 + t122 * t167) * qJD(1) - t174;
t119 = t122 ^ 2;
t148 = qJD(1) * t119 - t107 * t124;
t220 = t122 * t84;
t261 = (t121 * t148 + t220) * qJD(2) - t122 * t174 - t124 * t54;
t6 = (t122 * t86 + t123 * t148) * qJD(2) + t107 * t176 + t53 * t124;
t260 = -0.2e1 * t186;
t259 = t54 * qJ(5) + t84 * qJD(5);
t115 = pkin(6) * t198;
t228 = qJD(2) * pkin(2);
t97 = t115 - t228;
t144 = qJ(4) * t86 - t97;
t22 = -t241 * t84 + qJD(5) + t144;
t258 = (qJD(5) + t22) * t86;
t257 = 0.2e1 * t266;
t164 = pkin(3) * t111;
t13 = -t164 + t169;
t102 = t107 * qJ(4);
t35 = -t102 + t41;
t256 = t107 * t35 + t13;
t40 = -t121 * t98 + t123 * t73;
t201 = qJD(4) - t40;
t242 = t86 ^ 2;
t253 = -t107 ^ 2 - t242;
t195 = qJD(2) * t122;
t252 = qJ(4) * t195 - qJD(4) * t124;
t251 = -t121 * qJD(4) - t116;
t243 = t84 ^ 2;
t250 = t242 - t243;
t209 = t121 * qJ(4);
t246 = -t123 * t241 - t209;
t100 = t123 * t111;
t208 = t121 * t124;
t245 = qJD(1) * (t107 * t208 - t220) - t107 * t193 - t100;
t240 = pkin(7) * t86;
t239 = pkin(6) * t121;
t238 = pkin(7) * t107;
t37 = pkin(3) * t84 - t144;
t237 = t37 * t86;
t236 = t86 * t84;
t235 = pkin(7) - qJ(5);
t179 = -pkin(3) - t239;
t166 = -pkin(4) + t179;
t89 = t160 * qJD(1);
t218 = t123 * t89;
t96 = t235 * t123;
t234 = -t218 + (-qJ(5) * t206 + t122 * t166) * qJD(1) - qJD(3) * t96 + t121 * qJD(5);
t113 = qJ(4) * t198;
t187 = t123 * qJD(5);
t207 = t122 * t123;
t71 = t121 * t89;
t233 = t113 + t71 + (-pkin(6) * t207 + qJ(5) * t208) * qJD(1) + t193 * t235 + t187;
t211 = qJ(4) * t123;
t140 = -t121 * t241 + t211;
t232 = t107 * t140 + t251;
t156 = pkin(3) * t121 - t211;
t231 = t107 * t156 - t251;
t230 = t121 * t90 + t92 * t192;
t229 = qJ(4) * t54;
t10 = t54 * pkin(3) + pkin(6) * t172 + t53 * qJ(4) - t86 * qJD(4);
t227 = t10 * t121;
t226 = t10 * t123;
t29 = qJ(5) * t84 + t41;
t21 = -t102 + t29;
t225 = t107 * t21;
t221 = t121 * t97;
t217 = t123 * t97;
t214 = t84 * qJ(4);
t110 = pkin(6) * t206;
t212 = qJD(3) * t110 + t92 * t193;
t61 = t121 * t92 + t110;
t210 = qJ(5) * t122;
t127 = qJD(1) ^ 2;
t205 = t124 * t127;
t126 = qJD(2) ^ 2;
t204 = t126 * t122;
t203 = t126 * t124;
t28 = qJ(5) * t86 + t40;
t202 = qJD(4) - t28;
t199 = -t124 ^ 2 + t119;
t191 = qJD(4) * t123;
t184 = t121 * t238;
t183 = t123 * t238;
t109 = pkin(6) * t208;
t182 = pkin(7) * t188;
t180 = t122 * t205;
t173 = t107 * t198;
t170 = -t111 + t236;
t60 = t123 * t92 - t109;
t168 = pkin(1) * t260;
t163 = t124 * t111;
t162 = -t123 * t90 + t212;
t51 = -qJ(4) * t124 + t61;
t161 = t179 * t122;
t159 = (qJD(3) * t84 - t53) * pkin(7);
t157 = pkin(3) * t123 + t209;
t33 = pkin(3) * t107 + t201;
t155 = -t121 * t35 + t123 * t33;
t154 = -t121 * t41 - t123 * t40;
t56 = -pkin(6) * t177 + t71;
t147 = pkin(6) + t156;
t5 = -pkin(4) * t54 - t10;
t146 = -t121 * t5 - t192 * t22;
t145 = t123 * t5 - t193 * t22;
t143 = qJ(5) * t53 + t169;
t59 = (-t107 - t197) * t195;
t139 = -pkin(6) + t140;
t15 = -t121 * t222 - t215;
t135 = t53 + t222;
t134 = pkin(6) * t100 - t121 * t76 - t73 * t192 + t193 * t98;
t26 = (-t122 * t188 - t124 * t193) * pkin(6) + t230;
t133 = t143 - t263;
t12 = t84 * t175 + (t122 * t54 + t194 * t84) * t121;
t9 = -t134 + t266;
t131 = -t107 * t40 + t134;
t118 = t124 * pkin(3);
t95 = t235 * t121;
t91 = -pkin(2) - t157;
t79 = pkin(2) - t246;
t64 = t147 * t122;
t55 = t165 + t218;
t52 = t118 - t60;
t50 = t139 * t122;
t46 = pkin(7) * t215;
t45 = pkin(3) * t86 + t214;
t44 = qJD(1) * t161 - t218;
t43 = t113 + t56;
t39 = t121 * t210 + t51;
t38 = t124 * pkin(4) + t109 + t118 + (-t92 - t210) * t123;
t31 = -t241 * t86 - t214;
t27 = t195 * t239 - t162;
t25 = (qJD(3) * t157 - t191) * t122 + t147 * t194;
t23 = qJD(2) * t161 + t162;
t20 = t26 + t252;
t17 = (t246 * qJD(3) + t191) * t122 + t139 * t194;
t16 = t107 * t241 + t202;
t14 = -t123 * t213 - t216;
t11 = t86 * t124 * t188 + (-t123 * t53 - t193 * t86) * t122;
t8 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t207 + (qJD(5) * t122 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t124) * t121 + t230 + t252;
t7 = (-qJ(5) * t194 - t90) * t123 + (qJ(5) * t193 + qJD(2) * t166 - t187) * t122 + t212;
t4 = -qJD(5) * t86 + t133;
t3 = t9 + t259;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t163, t199 * t260, t203, -0.2e1 * t163, -t204, 0, -pkin(6) * t203 + t122 * t168, pkin(6) * t204 + t124 * t168, 0, 0, t11, t262, t6, t12, -t261, t59, -t107 * t27 + t124 * t169 + (pkin(6) * t54 + t192 * t97) * t122 + ((pkin(6) * t84 + t221) * t124 + (t40 + (t60 + t109) * qJD(1)) * t122) * qJD(2), t107 * t26 - t124 * t134 + (-pkin(6) * t53 - t193 * t97) * t122 + ((pkin(6) * t86 + t217) * t124 + (-t41 + (-t61 + t110) * qJD(1)) * t122) * qJD(2), -t26 * t84 - t27 * t86 + t53 * t60 - t54 * t61 + t154 * t194 + (t121 * t134 + t123 * t169 + (t121 * t40 - t123 * t41) * qJD(3)) * t122, -t134 * t61 - t169 * t60 + t26 * t41 + t27 * t40 + (t97 + t115) * pkin(6) * t194, t11, t6, -t262, t59, t261, t12, t107 * t23 + t25 * t84 + t54 * t64 + (t196 * t37 + t13) * t124 + (t37 * t192 + t227 + (-qJD(1) * t52 - t33) * qJD(2)) * t122, -t20 * t84 + t23 * t86 - t51 * t54 - t52 * t53 + t155 * t194 + (-t121 * t9 + t123 * t13 + (-t121 * t33 - t123 * t35) * qJD(3)) * t122, -t107 * t20 - t25 * t86 + t53 * t64 + (-t188 * t37 - t9) * t124 + (t37 * t193 - t226 + (qJD(1) * t51 + t35) * qJD(2)) * t122, t10 * t64 + t13 * t52 + t20 * t35 + t23 * t33 + t25 * t37 + t51 * t9, t11, -t262, -t6, t12, -t261, t59, t107 * t7 - t17 * t84 - t50 * t54 + (-t196 * t22 + t4) * t124 + ((-qJD(1) * t38 - t16) * qJD(2) + t146) * t122, -t107 * t8 + t17 * t86 - t50 * t53 + (t188 * t22 - t3) * t124 + ((qJD(1) * t39 + t21) * qJD(2) + t145) * t122, t38 * t53 + t39 * t54 - t7 * t86 + t8 * t84 + (t121 * t21 - t123 * t16) * t194 + (t121 * t3 - t123 * t4 + (t121 * t16 + t123 * t21) * qJD(3)) * t122, t16 * t7 + t17 * t22 + t21 * t8 + t3 * t39 + t38 * t4 + t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t199 * t127, 0, t180, 0, 0, t127 * pkin(1) * t122, pkin(1) * t205, 0, 0, t14, -t2, t24, t15, -t245, t173, -pkin(2) * t54 + t107 * t55 + (t183 + t221) * qJD(3) + ((-pkin(7) * t196 - t40) * t122 + (-t221 + (-t84 - t188) * pkin(6)) * t124) * qJD(1), pkin(2) * t53 - t107 * t56 + (-t184 + t217) * qJD(3) + ((t41 - t182) * t122 + (pkin(6) * t167 - t217) * t124) * qJD(1), t55 * t86 + t56 * t84 - t46 + (t40 * t197 - t134 + (-t40 + t240) * qJD(3)) * t123 + (t159 - t142) * t121, -t40 * t55 - t41 * t56 + (-t97 - t228) * t116 + (qJD(3) * t154 + t121 * t169 - t123 * t134) * pkin(7), t14, t24, t2, t173, t245, t15, -t226 - t107 * t44 + t54 * t91 - t231 * t84 + (t121 * t37 + t183) * qJD(3) + (t122 * t33 + (-pkin(7) * t195 - t124 * t37) * t121) * qJD(1), t43 * t84 - t44 * t86 - t46 + (-t33 * t197 + t9 + (t33 + t240) * qJD(3)) * t123 + (t159 + t256) * t121, -t227 + t107 * t43 + t53 * t91 + t231 * t86 + (-t123 * t37 + t184) * qJD(3) + (t37 * t206 + (-t35 + t182) * t122) * qJD(1), t10 * t91 - t33 * t44 - t35 * t43 - t231 * t37 + (qJD(3) * t155 + t121 * t13 + t123 * t9) * pkin(7), t14, t2, -t24, t15, -t245, t173, -t54 * t79 + t232 * t84 - t234 * t107 + (t22 * t208 + (-qJD(2) * t95 + t16) * t122) * qJD(1) + t145, -t53 * t79 - t232 * t86 + t233 * t107 + (-t22 * t206 + (qJD(2) * t96 - t21) * t122) * qJD(1) - t146, t53 * t95 + t54 * t96 + t234 * t86 - t233 * t84 + (t107 * t16 - t3) * t123 + (-t4 - t225) * t121, -t16 * t234 - t21 * t233 - t22 * t232 + t3 * t96 + t4 * t95 + t5 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, t250, -t135, -t236, -t268, t111, -t86 * t97 + t142, t84 * t97 + t131, 0, 0, t236, -t135, -t250, t111, t268, -t236, -t45 * t84 + t142 + 0.2e1 * t164 - t237, pkin(3) * t53 - t229 + (t35 - t41) * t86 + (t33 - t201) * t84, -t37 * t84 + t45 * t86 - t131 + t257, -pkin(3) * t13 + qJ(4) * t9 + t201 * t35 - t33 * t41 - t37 * t45, t236, -t250, t135, -t236, -t268, t111, -t107 * t29 + t31 * t84 - t143 + t258 + 0.2e1 * t263, t107 * t28 + t22 * t84 - t31 * t86 - t134 + t257 + t259, t229 - t241 * t53 + (-t21 + t29) * t86 + (-t16 + t202) * t84, qJ(4) * t3 - t16 * t29 + t202 * t21 - t22 * t31 - t241 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, -t135, t253, t237 + t256, 0, 0, 0, 0, 0, 0, t170, t253, t135, t133 + t225 - t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t264, -t265, -t242 - t243, t16 * t86 - t21 * t84 + t5;];
tauc_reg = t1;
