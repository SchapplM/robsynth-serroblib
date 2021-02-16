% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:11
% EndTime: 2021-01-15 21:35:20
% DurationCPUTime: 2.20s
% Computational Cost: add. (3207->258), mult. (8496->366), div. (0->0), fcn. (6557->8), ass. (0->162)
t143 = cos(qJ(5));
t177 = qJD(5) * t143;
t141 = sin(qJ(4));
t138 = sin(pkin(9));
t139 = cos(pkin(9));
t142 = sin(qJ(2));
t145 = cos(qJ(2));
t116 = t138 * t145 + t139 * t142;
t182 = qJD(1) * t116;
t188 = t139 * t145;
t172 = qJD(1) * t188;
t180 = qJD(1) * t142;
t105 = -t138 * t180 + t172;
t144 = cos(qJ(4));
t93 = t144 * t105;
t64 = -t141 * t182 + t93;
t229 = t143 * t64;
t236 = t177 - t229;
t135 = qJD(2) + qJD(4);
t191 = t64 * t135;
t179 = qJD(4) * t141;
t106 = t116 * qJD(2);
t96 = qJD(1) * t106;
t176 = qJD(1) * qJD(2);
t171 = t142 * t176;
t125 = t138 * t171;
t170 = t145 * t176;
t97 = t139 * t170 - t125;
t32 = qJD(4) * t93 - t141 * t96 + t144 * t97 - t179 * t182;
t235 = t32 - t191;
t184 = -qJD(5) + t64;
t234 = qJD(5) + t184;
t156 = t141 * t105 + t144 * t182;
t140 = sin(qJ(5));
t178 = qJD(5) * t140;
t14 = t135 * t177 + t143 * t32 - t156 * t178;
t52 = t140 * t135 + t143 * t156;
t15 = qJD(5) * t52 + t140 * t32;
t49 = -t143 * t135 + t140 * t156;
t233 = t14 * t143 - t140 * t15 - t236 * t49;
t12 = t14 * t140;
t232 = t236 * t52 + t12;
t33 = t156 * qJD(4) + t141 * t97 + t144 * t96;
t29 = t140 * t33;
t55 = t184 * t177;
t200 = t29 - t55;
t205 = t52 * t156;
t231 = t184 * t229 + t200 - t205;
t209 = t182 * pkin(7);
t203 = -qJ(3) - pkin(6);
t124 = t203 * t145;
t121 = qJD(1) * t124;
t110 = t138 * t121;
t123 = t203 * t142;
t120 = qJD(1) * t123;
t198 = qJD(2) * pkin(2);
t114 = t120 + t198;
t66 = t139 * t114 + t110;
t44 = qJD(2) * pkin(3) - t209 + t66;
t210 = t105 * pkin(7);
t189 = t139 * t121;
t67 = t138 * t114 - t189;
t47 = t67 + t210;
t20 = -t141 * t47 + t144 * t44;
t18 = -t135 * pkin(4) - t20;
t230 = t18 * t64;
t228 = t156 * t64;
t227 = t140 * t184;
t192 = t156 * t135;
t226 = -t33 + t192;
t224 = t156 ^ 2 - t64 ^ 2;
t35 = pkin(4) * t156 - t64 * pkin(8);
t164 = qJD(2) * t203;
t102 = t145 * qJD(3) + t142 * t164;
t84 = t102 * qJD(1);
t103 = -t142 * qJD(3) + t145 * t164;
t85 = t103 * qJD(1);
t48 = -t138 * t84 + t139 * t85;
t39 = -t97 * pkin(7) + t48;
t51 = t138 * t85 + t139 * t84;
t40 = -t96 * pkin(7) + t51;
t2 = (qJD(4) * t44 + t40) * t144 + t141 * t39 - t47 * t179;
t132 = -t145 * pkin(2) - pkin(1);
t181 = qJD(1) * t132;
t122 = qJD(3) + t181;
t72 = -t105 * pkin(3) + t122;
t223 = -t72 * t64 - t2;
t221 = -0.2e1 * t176;
t204 = t156 * t49;
t219 = t184 * t156;
t31 = t143 * t33;
t218 = -t178 * t184 - t31;
t21 = t141 * t44 + t144 * t47;
t19 = t135 * pkin(8) + t21;
t22 = -pkin(4) * t64 - pkin(8) * t156 + t72;
t157 = t140 * t19 - t143 * t22;
t217 = t156 * t157 + t18 * t178;
t3 = t21 * qJD(4) + t141 * t40 - t144 * t39;
t5 = t140 * t22 + t143 * t19;
t216 = t3 * t140 + t5 * t156 + t18 * t177;
t215 = -t156 * t72 - t3;
t115 = t138 * t142 - t188;
t68 = t144 * t115 + t141 * t116;
t69 = -t141 * t115 + t144 * t116;
t87 = t115 * pkin(3) + t132;
t26 = t68 * pkin(4) - t69 * pkin(8) + t87;
t74 = t139 * t123 + t138 * t124;
t58 = -t116 * pkin(7) + t74;
t75 = t138 * t123 - t139 * t124;
t59 = -t115 * pkin(7) + t75;
t28 = t141 * t58 + t144 * t59;
t109 = t115 * qJD(2);
t36 = -t68 * qJD(4) - t141 * t106 - t144 * t109;
t27 = t141 * t59 - t144 * t58;
t56 = -t138 * t102 + t139 * t103;
t41 = t109 * pkin(7) + t56;
t57 = t139 * t102 + t138 * t103;
t42 = -t106 * pkin(7) + t57;
t6 = -t27 * qJD(4) + t141 * t41 + t144 * t42;
t214 = (qJD(5) * t26 + t6) * t184 - (qJD(5) * t22 + t2) * t68 + t18 * t36 - t28 * t33 + t3 * t69;
t212 = pkin(2) * t138;
t211 = pkin(2) * t142;
t208 = t18 * t69;
t207 = t26 * t33;
t206 = t33 * t69;
t131 = t139 * pkin(2) + pkin(3);
t153 = t144 * t131 - t141 * t212;
t70 = -t138 * t120 + t189;
t53 = t70 - t210;
t71 = t139 * t120 + t110;
t54 = t71 - t209;
t201 = -t153 * qJD(4) + t141 * t53 + t144 * t54;
t154 = t141 * t131 + t144 * t212;
t199 = t154 * qJD(4) - t141 * t54 + t144 * t53;
t196 = t140 * t52;
t147 = qJD(1) ^ 2;
t187 = t145 * t147;
t146 = qJD(2) ^ 2;
t186 = t146 * t142;
t185 = t146 * t145;
t183 = t142 ^ 2 - t145 ^ 2;
t134 = t142 * t198;
t133 = pkin(2) * t180;
t173 = t69 * t178;
t128 = pkin(2) * t171;
t73 = t96 * pkin(3) + t128;
t80 = t106 * pkin(3) + t134;
t79 = pkin(3) * t182 + t133;
t161 = pkin(1) * t221;
t101 = pkin(8) + t154;
t160 = qJD(5) * t101 + t35 + t79;
t159 = 0.2e1 * t182;
t158 = -t184 * t36 + t206;
t155 = -t227 * t64 - t218;
t151 = -t101 * t33 - t184 * t201 - t230;
t100 = -pkin(4) - t153;
t37 = t69 * qJD(4) + t144 * t106 - t141 * t109;
t10 = t37 * pkin(4) - t36 * pkin(8) + t80;
t9 = t33 * pkin(4) - t32 * pkin(8) + t73;
t8 = t143 * t9;
t7 = t28 * qJD(4) + t141 * t42 - t144 * t41;
t1 = [0, 0, 0, 0.2e1 * t142 * t170, t183 * t221, t185, -t186, 0, -pkin(6) * t185 + t142 * t161, pkin(6) * t186 + t145 * t161, t122 * t106 + t132 * t96 + (t56 + (qJD(1) * t115 - t105) * t211) * qJD(2), -t122 * t109 + t132 * t97 + (t159 * t211 - t57) * qJD(2), t57 * t105 - t67 * t106 + t66 * t109 - t51 * t115 - t48 * t116 - t182 * t56 - t74 * t97 - t75 * t96, t48 * t74 + t51 * t75 + t66 * t56 + t67 * t57 + (t122 + t181) * t134, t156 * t36 + t32 * t69, -t156 * t37 - t32 * t68 + t36 * t64 - t206, t36 * t135, -t37 * t135, 0, -t7 * t135 + t87 * t33 + t72 * t37 - t64 * t80 + t73 * t68, -t6 * t135 + t156 * t80 + t87 * t32 + t72 * t36 + t73 * t69, -t52 * t173 + (t14 * t69 + t36 * t52) * t143, (-t143 * t49 - t196) * t36 + (-t12 - t143 * t15 + (t140 * t49 - t143 * t52) * qJD(5)) * t69, t14 * t68 + t143 * t158 + t173 * t184 + t52 * t37, -t140 * t158 - t15 * t68 - t49 * t37 + t55 * t69, -t184 * t37 + t33 * t68, t27 * t15 - t157 * t37 + t7 * t49 + t8 * t68 + (-t10 * t184 + t207 + (t184 * t28 - t19 * t68 + t208) * qJD(5)) * t143 + t214 * t140, t27 * t14 - t5 * t37 + t7 * t52 + ((-qJD(5) * t28 + t10) * t184 - t207 - (-qJD(5) * t19 + t9) * t68 - qJD(5) * t208) * t140 + t214 * t143; 0, 0, 0, -t142 * t187, t183 * t147, 0, 0, 0, t147 * pkin(1) * t142, pkin(1) * t187, -t70 * qJD(2) + t105 * t133 - t122 * t182 + t48, t71 * qJD(2) - t122 * t105 - t133 * t182 - t51, (t67 + t70) * t182 + (t66 - t71) * t105 + (-t138 * t96 - t139 * t97) * pkin(2), -t66 * t70 - t67 * t71 + (-t122 * t180 + t138 * t51 + t139 * t48) * pkin(2), -t228, t224, t235, t226, 0, -t199 * t135 + t64 * t79 + t215, t201 * t135 - t156 * t79 + t223, t232, t184 * t196 + t233, t231, t155 + t204, t219, t100 * t15 + t199 * t49 + (t160 * t184 - t3) * t143 + t151 * t140 + t217, t100 * t14 + t151 * t143 - t160 * t227 + t199 * t52 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159 * qJD(2), -t125 + (t105 + t172) * qJD(2), -t105 ^ 2 - t182 ^ 2, -t67 * t105 + t182 * t66 + t128, 0, 0, 0, 0, 0, t33 + t192, t32 + t191, 0, 0, 0, 0, 0, t155 - t204, -t143 * t184 ^ 2 - t205 - t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, t224, t235, t226, 0, t21 * t135 + t215, t20 * t135 + t223, t232, t227 * t52 + t233, t231, -t184 * t227 + t204 + t31, t219, -pkin(4) * t15 - t3 * t143 + (-t140 * t20 + t143 * t35) * t184 - t21 * t49 - t140 * t230 - t200 * pkin(8) + t217, -pkin(4) * t14 - (t140 * t35 + t143 * t20) * t184 - t21 * t52 - t18 * t229 + t218 * pkin(8) + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t49, -t49 ^ 2 + t52 ^ 2, -t184 * t49 + t14, -t184 * t52 - t15, t33, -t140 * t2 - t18 * t52 - t234 * t5 + t8, -t140 * t9 - t143 * t2 + t234 * t157 + t18 * t49;];
tauc_reg = t1;
