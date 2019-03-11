% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:03
% EndTime: 2019-03-08 19:23:08
% DurationCPUTime: 2.10s
% Computational Cost: add. (3003->290), mult. (6322->432), div. (0->0), fcn. (4459->10), ass. (0->165)
t82 = sin(pkin(6));
t154 = qJD(1) * t82;
t91 = cos(qJ(2));
t134 = t91 * t154;
t88 = sin(qJ(2));
t135 = t88 * t154;
t81 = sin(pkin(11));
t83 = cos(pkin(11));
t40 = t81 * t134 - t83 * t135;
t204 = qJD(3) * t81 - t40;
t109 = qJD(3) - t134;
t92 = -pkin(2) - pkin(3);
t54 = t92 * qJD(2) + t109;
t66 = qJD(2) * qJ(3) + t135;
t31 = t81 * t54 + t83 * t66;
t29 = -qJD(2) * pkin(8) + t31;
t84 = cos(pkin(6));
t72 = -qJD(1) * t84 + qJD(4);
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t22 = t29 * t90 + t72 * t87;
t18 = qJD(5) * pkin(9) + t22;
t116 = pkin(5) * t90 + pkin(9) * t87;
t30 = t54 * t83 - t81 * t66;
t28 = qJD(2) * pkin(4) - t30;
t23 = t116 * qJD(2) + t28;
t86 = sin(qJ(6));
t89 = cos(qJ(6));
t108 = t18 * t86 - t23 * t89;
t115 = -pkin(5) * t87 + pkin(9) * t90;
t100 = t115 * qJD(5);
t152 = qJD(2) * t82;
t133 = t88 * t152;
t118 = qJD(1) * t133;
t57 = (qJD(3) + t134) * qJD(2);
t32 = -t83 * t118 + t57 * t81;
t25 = qJD(2) * t100 + t32;
t105 = t29 * t87 - t72 * t90;
t33 = t81 * t118 + t83 * t57;
t9 = -t105 * qJD(5) + t33 * t90;
t1 = -t108 * qJD(6) + t25 * t86 + t89 * t9;
t150 = qJD(2) * t90;
t73 = qJD(6) + t150;
t203 = t108 * t73 + t1;
t202 = t204 * qJD(2) + t32;
t159 = t83 * qJ(3) + t81 * t92;
t62 = -pkin(8) + t159;
t201 = -qJD(6) * t62 * t90 + t100 + t204;
t148 = qJD(3) * t83;
t46 = (t81 * t88 + t83 * t91) * t82;
t43 = qJD(1) * t46;
t121 = -t43 + t148;
t200 = qJD(2) * t121;
t6 = t18 * t89 + t23 * t86;
t2 = -qJD(6) * t6 + t89 * t25 - t86 * t9;
t199 = -t6 * t73 - t2;
t79 = t87 ^ 2;
t198 = (qJD(2) * t79 - t73 * t90) * t86;
t93 = qJD(5) ^ 2;
t197 = -t62 * t93 + t202;
t10 = qJD(5) * t22 + t87 * t33;
t190 = t10 * t87;
t95 = -(-t105 * t90 + t22 * t87) * qJD(5) + t9 * t90 + t190;
t165 = t86 * t90;
t147 = qJD(5) * t87;
t124 = -t81 * qJ(3) + t83 * t92;
t61 = pkin(4) - t124;
t44 = t116 + t61;
t97 = qJD(6) * t44 - t62 * t147 + t90 * t148;
t194 = t43 * t165 + t201 * t89 - t97 * t86;
t164 = t89 * t90;
t193 = t43 * t164 - t201 * t86 - t97 * t89;
t168 = t81 * t91;
t47 = (t83 * t88 - t168) * t82;
t34 = t47 * t87 + t84 * t90;
t192 = t10 * t34;
t191 = t10 * t86;
t189 = t10 * t89;
t17 = -qJD(5) * pkin(5) + t105;
t188 = t17 * t86;
t187 = t17 * t89;
t186 = t32 * t46;
t140 = qJD(2) * qJD(5);
t127 = t89 * t140;
t146 = qJD(5) * t89;
t151 = qJD(2) * t87;
t59 = t86 * t151 + t146;
t144 = qJD(6) * t59;
t38 = -t90 * t127 + t144;
t185 = t38 * t86;
t184 = t38 * t90;
t141 = t86 * qJD(5);
t142 = qJD(6) * t89;
t129 = t87 * t142;
t99 = t90 * t141 + t129;
t39 = t99 * qJD(2) - qJD(6) * t141;
t183 = t39 * t89;
t182 = t39 * t90;
t60 = t89 * t151 - t141;
t181 = t59 * t60;
t180 = t59 * t73;
t179 = t59 * t86;
t178 = t59 * t87;
t177 = t59 * t89;
t176 = t60 * t73;
t175 = t60 * t86;
t174 = t60 * t87;
t173 = t60 * t89;
t172 = t66 * t91;
t171 = t73 * t86;
t170 = t73 * t89;
t169 = t81 * t87;
t94 = qJD(2) ^ 2;
t167 = t82 * t94;
t166 = t83 * t94;
t163 = t93 * t87;
t162 = t93 * t90;
t132 = t81 * t147;
t52 = -t81 * t165 - t83 * t89;
t161 = t52 * qJD(6) - t89 * t132 - (t83 * t164 + t81 * t86) * qJD(2);
t53 = t81 * t164 - t83 * t86;
t160 = -t53 * qJD(6) + t86 * t132 - (-t83 * t165 + t81 * t89) * qJD(2);
t80 = t90 ^ 2;
t158 = t79 - t80;
t157 = t79 + t80;
t156 = t93 + t94;
t155 = qJD(2) * pkin(2);
t153 = qJD(2) * t28;
t145 = qJD(5) * t90;
t143 = qJD(6) * t86;
t139 = t73 * t164;
t138 = t88 * t167;
t137 = t91 * t167;
t136 = t87 * t94 * t90;
t131 = t81 * t145;
t130 = t87 * t143;
t128 = 0.2e1 * t140;
t126 = t87 * t140;
t123 = -t33 + t153;
t120 = t73 * t129;
t119 = t83 * t128;
t117 = t90 * t126;
t113 = -t108 * t89 + t6 * t86;
t112 = t108 * t86 + t6 * t89;
t106 = -t105 * t87 - t22 * t90;
t104 = t30 * t81 - t31 * t83;
t35 = t47 * t90 - t84 * t87;
t16 = t35 * t89 + t46 * t86;
t15 = -t35 * t86 + t46 * t89;
t101 = t79 * t127 + t73 * t130;
t98 = qJD(5) * (-qJD(2) * t61 - t121 - t28);
t96 = -t113 * qJD(6) + t1 * t89 - t2 * t86;
t65 = t115 * qJD(2);
t58 = t109 - t155;
t42 = qJD(2) * t46;
t41 = -t83 * t133 + t152 * t168;
t27 = t62 * t164 + t44 * t86;
t26 = -t62 * t165 + t44 * t89;
t14 = -t34 * qJD(5) + t42 * t90;
t13 = t47 * t145 - t84 * t147 + t42 * t87;
t12 = -t105 * t89 + t65 * t86;
t11 = t105 * t86 + t65 * t89;
t4 = t15 * qJD(6) + t14 * t89 + t41 * t86;
t3 = -t16 * qJD(6) - t14 * t86 + t41 * t89;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, 0, 0, 0, 0, 0, 0, 0, 0, -t138, 0, t137 (t57 * t88 + (t172 + (t58 - t134) * t88) * qJD(2)) * t82, 0, 0, 0, 0, 0, 0, t41 * qJD(2), t42 * qJD(2), 0, -t30 * t41 + t31 * t42 + t33 * t47 + t186, 0, 0, 0, 0, 0, 0, -qJD(5) * t13 + (-t147 * t46 + t41 * t90) * qJD(2), -qJD(5) * t14 + (-t145 * t46 - t41 * t87) * qJD(2) (-t13 * t87 - t14 * t90 + (-t34 * t90 + t35 * t87) * qJD(5)) * qJD(2), t105 * t13 + t14 * t22 + t28 * t41 + t35 * t9 + t186 + t192, 0, 0, 0, 0, 0, 0, -t126 * t15 - t13 * t59 + t3 * t73 - t34 * t39, t126 * t16 - t13 * t60 + t34 * t38 - t4 * t73, -t15 * t38 + t16 * t39 + t3 * t60 + t4 * t59, t1 * t16 - t108 * t3 + t13 * t17 + t15 * t2 + t4 * t6 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), qJ(3) * t57 + qJD(3) * t66 + (-t172 + (-t58 - t155) * t88) * t154, 0, 0, 0, 0, 0, 0, t202, t33 + t200, 0, -qJD(3) * t104 - t124 * t32 + t159 * t33 + t30 * t40 - t31 * t43, 0.2e1 * t117, -t158 * t128, -t162, -0.2e1 * t117, t163, 0, t197 * t90 + t87 * t98, -t197 * t87 + t90 * t98, -t157 * t200 - t95, -t28 * t40 + t32 * t61 + t106 * t43 + (-t106 * t83 + t28 * t81) * qJD(3) + t95 * t62, -t38 * t87 * t89 + (t145 * t89 - t130) * t60 (-t175 - t177) * t145 + (t185 - t183 + (-t173 + t179) * qJD(6)) * t87, t184 + (-t139 + t174) * qJD(5) + t101, t39 * t86 * t87 + t59 * t99, t120 + t182 + (-t178 - t198) * qJD(5) (-t73 - t150) * t147, t194 * t73 + (t2 + (-t59 * t62 - t188) * qJD(5)) * t90 + (-t17 * t142 - t191 - t39 * t62 - t121 * t59 + (-qJD(2) * t26 + t108) * qJD(5)) * t87, t193 * t73 + (-t1 + (-t60 * t62 - t187) * qJD(5)) * t90 + (t17 * t143 - t189 + t38 * t62 - t121 * t60 + (qJD(2) * t27 + t6) * qJD(5)) * t87, -t26 * t38 + t27 * t39 + t194 * t60 - t193 * t59 + t113 * t145 + (qJD(6) * t112 + t1 * t86 + t2 * t89) * t87, t62 * t190 + t1 * t27 + t2 * t26 - t193 * t6 - t194 * t108 + (t121 * t87 + t145 * t62) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 (-t66 + t135) * qJD(2), 0, 0, 0, 0, 0, 0, -t81 * t94, -t166, 0, qJD(2) * t104 - t32 * t83 + t33 * t81, 0, 0, 0, 0, 0, 0, -t156 * t81 * t90 + t119 * t87, t119 * t90 + t156 * t169, t157 * t166 (qJD(2) * t106 - t32) * t83 + (t95 - t153) * t81, 0, 0, 0, 0, 0, 0, -t59 * t131 + t160 * t73 + (-t39 * t81 + (-qJD(5) * t52 + t59 * t83) * qJD(2)) * t87, -t60 * t131 - t161 * t73 + (t38 * t81 + (qJD(5) * t53 + t60 * t83) * qJD(2)) * t87, t160 * t60 + t161 * t59 - t38 * t52 + t39 * t53, t10 * t169 + t1 * t53 + t2 * t52 + t161 * t6 - t160 * t108 + (-t151 * t83 + t131) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t162, 0, -qJD(5) * t106 - t10 * t90 + t87 * t9, 0, 0, 0, 0, 0, 0, -t120 + t182 + (-t178 + t198) * qJD(5), -t184 + (-t139 - t174) * qJD(5) + t101 (-t175 + t177) * t145 + (t185 + t183 + (-t173 - t179) * qJD(6)) * t87 (qJD(5) * t112 - t10) * t90 + (qJD(5) * t17 + t96) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t158 * t94, 0, t136, 0, 0, t123 * t87, t123 * t90, 0, 0, -t170 * t60 + t185 (t38 + t180) * t89 + (t39 + t176) * t86, t73 * t142 + (t139 + (-t60 - t141) * t87) * qJD(2), -t171 * t59 + t183, -t73 * t143 + (-t73 * t165 + (t59 - t146) * t87) * qJD(2), t73 * t151, pkin(5) * t39 - t189 - t11 * t73 + t22 * t59 + (-pkin(9) * t170 + t188) * qJD(6) + (-t108 * t87 + (pkin(9) * t147 + t17 * t90) * t86) * qJD(2), -pkin(5) * t38 + t191 + t12 * t73 + t22 * t60 + (pkin(9) * t171 + t187) * qJD(6) + (t17 * t164 + (pkin(9) * t146 - t6) * t87) * qJD(2), -t11 * t60 - t12 * t59 + ((-qJD(6) * t60 + t39) * pkin(9) + t203) * t89 + ((t38 - t144) * pkin(9) + t199) * t86, -pkin(5) * t10 + pkin(9) * t96 + t108 * t11 - t12 * t6 - t17 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -t59 ^ 2 + t60 ^ 2, t38 - t180, -t181, t39 - t176, -t126, t17 * t60 - t199, -t17 * t59 - t203, 0, 0;];
tauc_reg  = t5;
