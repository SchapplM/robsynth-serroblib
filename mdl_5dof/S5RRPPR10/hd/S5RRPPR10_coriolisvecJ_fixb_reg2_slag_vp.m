% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:41
% EndTime: 2019-12-31 19:44:49
% DurationCPUTime: 2.96s
% Computational Cost: add. (2742->356), mult. (6984->498), div. (0->0), fcn. (4515->6), ass. (0->202)
t138 = sin(pkin(8));
t143 = cos(qJ(2));
t195 = qJD(1) * qJD(2);
t183 = t143 * t195;
t118 = t138 * t183;
t140 = sin(qJ(5));
t142 = cos(qJ(5));
t141 = sin(qJ(2));
t205 = qJD(1) * t141;
t188 = t138 * t205;
t139 = cos(pkin(8));
t196 = t139 * qJD(2);
t93 = t188 - t196;
t186 = t139 * t205;
t203 = qJD(2) * t138;
t95 = t186 + t203;
t160 = t140 * t93 + t142 * t95;
t173 = t139 * t183;
t16 = t160 * qJD(5) - t142 * t118 + t140 * t173;
t204 = qJD(1) * t143;
t126 = qJD(5) + t204;
t230 = t126 * t160;
t257 = -t16 + t230;
t136 = t141 ^ 2;
t137 = t143 ^ 2;
t174 = qJD(1) * (t136 - 0.2e1 * t137);
t256 = (t138 * t174 + t141 * t93) * qJD(2);
t134 = t138 ^ 2;
t135 = t139 ^ 2;
t228 = t139 * t93;
t255 = ((t134 - t135) * qJD(2) - t138 * t95 - t228) * t204;
t201 = qJD(2) * t143;
t254 = ((t95 + 0.2e1 * t186) * t138 + t228) * t201;
t252 = t160 ^ 2;
t161 = t140 * t95 - t142 * t93;
t251 = t161 ^ 2;
t250 = -0.2e1 * t195;
t214 = t139 * t143;
t194 = pkin(7) * t214;
t241 = pkin(3) + pkin(4);
t129 = pkin(6) * t205;
t106 = (qJD(3) - t129) * qJD(2);
t166 = pkin(2) * t141 - qJ(3) * t143;
t81 = t166 * qJD(2) - t141 * qJD(3);
t72 = t81 * qJD(1);
t32 = -t138 * t106 + t139 * t72;
t12 = (-t241 * t141 - t194) * t195 - t32;
t182 = t141 * t195;
t33 = t139 * t106 + t138 * t72;
t193 = qJ(4) * t182 + t33;
t13 = (pkin(7) * t203 - qJD(4)) * t204 + t193;
t130 = pkin(6) * t204;
t115 = qJD(2) * qJ(3) + t130;
t224 = qJ(3) * t141;
t109 = -pkin(2) * t143 - pkin(1) - t224;
t87 = t109 * qJD(1);
t50 = -t138 * t115 + t139 * t87;
t35 = pkin(3) * t204 + qJD(4) - t50;
t14 = pkin(4) * t204 - pkin(7) * t95 + t35;
t51 = t139 * t115 + t138 * t87;
t37 = -qJ(4) * t204 + t51;
t19 = pkin(7) * t93 + t37;
t5 = t14 * t142 - t140 * t19;
t1 = t5 * qJD(5) + t12 * t140 + t13 * t142;
t249 = -t126 * t5 + t1;
t6 = t14 * t140 + t142 * t19;
t2 = -t6 * qJD(5) + t142 * t12 - t13 * t140;
t248 = t126 * t6 + t2;
t91 = t95 ^ 2;
t247 = -t93 ^ 2 - t91;
t197 = qJD(5) * t142;
t198 = qJD(5) * t140;
t15 = -t140 * t118 - t142 * t173 - t93 * t197 + t95 * t198;
t231 = t126 * t161;
t246 = t15 - t231;
t199 = qJD(4) * t143;
t202 = qJD(2) * t141;
t244 = qJ(4) * t202 - t199;
t125 = pkin(6) * t183;
t209 = pkin(3) * t118 + t125;
t221 = qJD(4) * t95;
t223 = qJ(4) * t139;
t20 = -(pkin(4) * t138 - t223) * t183 - t209 + t221;
t189 = -pkin(6) * t138 - pkin(3);
t147 = -t194 + (-pkin(4) + t189) * t141;
t101 = t166 * qJD(1);
t220 = t101 * t139;
t26 = t147 * qJD(1) - t220;
t127 = qJ(4) * t205;
t215 = t139 * t141;
t216 = t138 * t143;
t154 = -pkin(6) * t215 + pkin(7) * t216;
t85 = t138 * t101;
t36 = t154 * qJD(1) + t127 + t85;
t235 = -pkin(7) + qJ(3);
t112 = t235 * t138;
t113 = t235 * t139;
t59 = t112 * t142 - t113 * t140;
t97 = t138 * t140 + t139 * t142;
t240 = t97 * qJD(3) + t59 * qJD(5) - t140 * t26 - t142 * t36;
t60 = t112 * t140 + t113 * t142;
t217 = t138 * t142;
t98 = -t139 * t140 + t217;
t239 = t98 * qJD(3) - t60 * qJD(5) + t140 * t36 - t142 * t26;
t236 = t160 * t161;
t187 = t138 * t204;
t213 = t140 * t143;
t190 = t139 * t213;
t234 = qJD(1) * t190 - t138 * t197 + t139 * t198 - t142 * t187;
t153 = t97 * t143;
t83 = t97 * qJD(5);
t233 = qJD(1) * t153 + t83;
t232 = qJD(2) * pkin(2);
t229 = t139 * t81;
t227 = t141 * t35;
t226 = t141 * t37;
t124 = pkin(6) * t214;
t69 = t138 * t109 + t124;
t222 = qJD(3) * t95;
t145 = qJD(1) ^ 2;
t219 = t137 * t145;
t218 = t138 * t141;
t212 = t143 * t145;
t144 = qJD(2) ^ 2;
t211 = t144 * t141;
t210 = t144 * t143;
t206 = t136 - t137;
t200 = qJD(4) * t138;
t123 = pkin(6) * t216;
t192 = pkin(6) * t202;
t191 = t138 * t219;
t185 = t138 * t201;
t184 = qJD(4) * t215;
t181 = qJ(4) * t138 + pkin(2);
t178 = -qJD(3) + t232;
t168 = -t129 + t178;
t180 = t168 - t232;
t179 = pkin(1) * t250;
t177 = -t95 + t203;
t176 = t93 + t196;
t155 = -t241 * t138 + t223;
t175 = -t155 * t204 + t130 + t200;
t68 = t109 * t139 - t123;
t172 = t143 * t182;
t62 = -qJ(4) * t143 + t69;
t171 = t189 * t141;
t170 = -t95 * t204 + t118;
t165 = pkin(3) * t138 - t223;
t133 = t143 * pkin(3);
t38 = pkin(4) * t143 + t123 + t133 + (-pkin(7) * t141 - t109) * t139;
t48 = pkin(7) * t218 + t62;
t9 = -t140 * t48 + t142 * t38;
t10 = t140 * t38 + t142 * t48;
t159 = qJD(1) * t177;
t158 = qJD(1) * t176;
t65 = -pkin(6) * t186 + t85;
t73 = t138 * t81;
t57 = -t139 * t192 + t73;
t157 = t134 * t172 + t93 * t185;
t156 = pkin(6) + t165;
t152 = -pkin(6) + t155;
t150 = qJ(4) * t95 + t168;
t149 = t158 * t216;
t148 = -qJ(4) * t173 + t209;
t122 = t141 * t212;
t117 = qJD(3) * t187;
t116 = -0.2e1 * t172;
t107 = -pkin(3) * t139 - t181;
t86 = t241 * t139 + t181;
t79 = t93 * t204;
t78 = qJD(3) * t228;
t77 = t97 * t141;
t76 = t140 * t215 - t141 * t217;
t74 = t156 * t141;
t67 = t165 * t204 + t130;
t64 = pkin(6) * t188 + t220;
t63 = t133 - t68;
t61 = t152 * t141;
t58 = t79 + t173;
t56 = t138 * t192 + t229;
t55 = qJD(1) * t171 - t220;
t54 = t127 + t65;
t53 = t156 * t201 - t184;
t49 = (t135 * t205 + t139 * t95) * t201;
t47 = t159 * t214;
t41 = qJD(2) * t171 - t229;
t40 = t139 * t219 + t141 * t159;
t39 = (t139 * t174 + t141 * t95) * qJD(2);
t34 = t152 * t201 + t184;
t31 = t57 + t244;
t30 = pkin(3) * t93 - t150;
t29 = t98 * t141 * qJD(5) + qJD(2) * t153;
t28 = qJD(2) * t190 + t141 * t83 - t142 * t185;
t27 = t148 - t221;
t23 = -pkin(3) * t182 - t32;
t22 = t154 * qJD(2) + t244 + t73;
t21 = t147 * qJD(2) - t229;
t18 = -qJD(1) * t199 + t193;
t17 = -t241 * t93 + t150;
t4 = -t10 * qJD(5) - t140 * t22 + t142 * t21;
t3 = t9 * qJD(5) + t140 * t21 + t142 * t22;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t172, t206 * t250, t210, t116, -t211, 0, -pkin(6) * t210 + t141 * t179, pkin(6) * t211 + t143 * t179, 0, 0, t49, -t254, t39, t157, -t256, t116, (-qJD(1) * t56 - t32) * t143 + ((pkin(6) * t93 - t138 * t168) * t143 + (t50 + (t68 + 0.2e1 * t123) * qJD(1)) * t141) * qJD(2), (qJD(1) * t57 + t33) * t143 + ((pkin(6) * t95 - t139 * t168) * t143 + (-t51 + (-t69 + 0.2e1 * t124) * qJD(1)) * t141) * qJD(2), -t56 * t95 - t57 * t93 + (-t138 * t33 - t139 * t32) * t141 + (-t138 * t51 - t139 * t50 + (-t138 * t69 - t139 * t68) * qJD(1)) * t201, t32 * t68 + t33 * t69 + t50 * t56 + t51 * t57 + (-t168 + t129) * pkin(6) * t201, t49, t39, t254, t116, t256, t157, t27 * t218 + t53 * t93 + (qJD(1) * t41 + t23) * t143 + (t30 * t216 - t227 + (-t141 * t63 + t74 * t216) * qJD(1)) * qJD(2), -t31 * t93 + t41 * t95 + (-t138 * t18 + t139 * t23) * t141 + (-t138 * t37 + t139 * t35 + (-t138 * t62 + t139 * t63) * qJD(1)) * t201, -t27 * t215 - t53 * t95 + (-qJD(1) * t31 - t18) * t143 + (-t30 * t214 + t226 + (t141 * t62 - t74 * t214) * qJD(1)) * qJD(2), t18 * t62 + t23 * t63 + t27 * t74 + t30 * t53 + t31 * t37 + t35 * t41, -t15 * t77 + t160 * t29, t15 * t76 - t16 * t77 - t160 * t28 - t161 * t29, t126 * t29 - t143 * t15 + (-qJD(1) * t77 - t160) * t202, t16 * t76 + t161 * t28, -t126 * t28 - t143 * t16 + (qJD(1) * t76 + t161) * t202, (-t126 - t204) * t202, t126 * t4 + t143 * t2 + t16 * t61 + t17 * t28 + t20 * t76 + t34 * t161 + (-qJD(1) * t9 - t5) * t202, -t1 * t143 - t126 * t3 - t15 * t61 + t17 * t29 + t20 * t77 + t34 * t160 + (qJD(1) * t10 + t6) * t202, -t1 * t76 - t10 * t16 + t15 * t9 - t160 * t4 - t161 * t3 - t2 * t77 - t28 * t6 - t29 * t5, t1 * t10 + t17 * t34 + t2 * t9 + t20 * t61 + t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t206 * t145, 0, t122, 0, 0, t145 * pkin(1) * t141, pkin(1) * t212, 0, 0, t47, -t255, t40, -t149, t176 * t205 - t191, t122, t117 + ((-qJ(3) * t203 - t50) * t141 + (-t176 * pkin(6) + t180 * t138 + t64) * t143) * qJD(1), ((-qJ(3) * t196 + t51) * t141 + (-t65 + t177 * pkin(6) + (t168 - t178) * t139) * t143) * qJD(1), t64 * t95 + t65 * t93 - t78 + (t50 * t204 + t33) * t139 + (t51 * t204 + t222 - t32) * t138, -t50 * t64 - t51 * t65 + (-t138 * t50 + t139 * t51) * qJD(3) + (-t138 * t32 + t139 * t33) * qJ(3) + t180 * t130, t47, t40, t255, t122, -t141 * t158 + t191, -t149, -t139 * t27 + t117 + (-t67 - t200) * t93 + (t227 - t143 * t55 + (-t143 * t30 + (t107 * t143 - t224) * qJD(2)) * t138) * qJD(1), t54 * t93 - t55 * t95 - t78 + (-t35 * t204 + t18) * t139 + (t37 * t204 + t222 + t23) * t138, t67 * t95 + (-t27 + t221) * t138 + (-t226 + t143 * t54 + (qJ(3) * t202 + (-qJD(2) * t107 - qJD(3) + t30) * t143) * t139) * qJD(1), t107 * t27 - t30 * t67 - t35 * t55 - t37 * t54 + (qJ(3) * t18 + qJD(3) * t37) * t139 + (qJ(3) * t23 + qJD(3) * t35 - qJD(4) * t30) * t138, -t15 * t98 - t160 * t233, t15 * t97 - t16 * t98 + t160 * t234 + t161 * t233, -t233 * t126 + (-qJD(2) * t98 + t160) * t205, t16 * t97 - t161 * t234, t234 * t126 + (qJD(2) * t97 - t161) * t205, t126 * t205, t16 * t86 + t20 * t97 + t175 * t161 - t234 * t17 + t239 * t126 + (-qJD(2) * t59 + t5) * t205, -t15 * t86 + t20 * t98 + t175 * t160 - t233 * t17 - t240 * t126 + (qJD(2) * t60 - t6) * t205, -t1 * t97 + t15 * t59 - t16 * t60 - t160 * t239 - t161 * t240 - t2 * t98 + t233 * t5 + t234 * t6, t1 * t60 + t175 * t17 + t2 * t59 + t20 * t86 + t239 * t5 + t240 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t58, t247, t50 * t95 + t51 * t93 + t125, 0, 0, 0, 0, 0, 0, t170, t247, -t58, t37 * t93 + (-qJD(4) - t35) * t95 + t148, 0, 0, 0, 0, 0, 0, -t16 - t230, t15 + t231, t251 + t252, -t160 * t5 - t161 * t6 - t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t95 - t182, -t79 + t173, -t91 - t219, t30 * t95 + (-pkin(3) * t202 + t143 * t37) * qJD(1) - t32, 0, 0, 0, 0, 0, 0, -t126 * t198 - t161 * t95 + (-t126 * t213 - t142 * t202) * qJD(1), -t126 * t197 - t160 * t95 + (-t126 * t142 * t143 + t140 * t202) * qJD(1), t257 * t140 + t246 * t142, t249 * t140 + t248 * t142 - t17 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t251 + t252, -t246, -t236, t257, -t182, -t160 * t17 + t248, t161 * t17 - t249, 0, 0;];
tauc_reg = t7;
