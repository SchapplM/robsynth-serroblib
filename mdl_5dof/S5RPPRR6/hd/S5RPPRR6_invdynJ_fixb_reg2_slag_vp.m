% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:13
% EndTime: 2019-12-31 17:58:16
% DurationCPUTime: 1.99s
% Computational Cost: add. (3338->308), mult. (7169->391), div. (0->0), fcn. (5152->14), ass. (0->173)
t115 = qJ(1) + pkin(8);
t105 = sin(t115);
t107 = cos(t115);
t162 = -g(1) * t105 + g(2) * t107;
t119 = cos(pkin(8));
t204 = t119 * pkin(1);
t98 = -pkin(2) - t204;
t174 = qJDD(1) * t98;
t85 = qJDD(3) + t174;
t221 = -t162 - t85;
t114 = pkin(9) + qJ(4);
t104 = sin(t114);
t106 = cos(t114);
t150 = g(1) * t107 + g(2) * t105;
t131 = -g(3) * t106 + t150 * t104;
t122 = sin(qJ(4));
t211 = cos(qJ(4));
t118 = cos(pkin(9));
t101 = t118 * qJDD(2);
t116 = sin(pkin(9));
t117 = sin(pkin(8));
t94 = t117 * pkin(1) + qJ(3);
t80 = qJD(1) * qJD(3) + t94 * qJDD(1);
t50 = t101 + (-pkin(6) * qJDD(1) - t80) * t116;
t170 = t118 * qJDD(1);
t59 = t116 * qJDD(2) + t118 * t80;
t51 = pkin(6) * t170 + t59;
t160 = t122 * t51 - t211 * t50;
t103 = t118 * qJD(2);
t191 = pkin(6) * qJD(1);
t87 = t94 * qJD(1);
t53 = t103 + (-t87 - t191) * t116;
t65 = t116 * qJD(2) + t118 * t87;
t54 = t118 * t191 + t65;
t25 = t122 * t53 + t211 * t54;
t8 = -t25 * qJD(4) - t160;
t6 = -qJDD(4) * pkin(4) - t8;
t130 = -t6 + t131;
t179 = t122 * t116;
t163 = qJD(1) * t179;
t164 = t211 * t118;
t92 = qJD(1) * t164;
t72 = -t92 + t163;
t68 = qJD(5) + t72;
t220 = -pkin(7) * qJD(5) * t68 + t130;
t121 = sin(qJ(5));
t124 = cos(qJ(5));
t23 = qJD(4) * pkin(7) + t25;
t97 = t118 * pkin(3) + pkin(2);
t86 = -t97 - t204;
t70 = qJD(1) * t86 + qJD(3);
t83 = t211 * t116 + t122 * t118;
t74 = t83 * qJD(1);
t28 = t72 * pkin(4) - t74 * pkin(7) + t70;
t10 = t121 * t28 + t124 * t23;
t159 = qJDD(1) * t211;
t168 = qJD(4) * t92 + t116 * t159 + t122 * t170;
t41 = qJD(4) * t163 - t168;
t171 = t116 * qJDD(1);
t144 = -t118 * t159 + t122 * t171;
t77 = t83 * qJD(4);
t42 = qJD(1) * t77 + t144;
t67 = qJDD(1) * t86 + qJDD(3);
t14 = t42 * pkin(4) + t41 * pkin(7) + t67;
t13 = t124 * t14;
t161 = qJD(4) * t211;
t169 = -t122 * t50 - t53 * t161 - t211 * t51;
t177 = qJD(4) * t122;
t7 = -t54 * t177 - t169;
t5 = qJDD(4) * pkin(7) + t7;
t2 = -qJD(5) * t10 - t121 * t5 + t13;
t207 = t10 * t68;
t219 = t2 + t207;
t156 = t121 * t68;
t57 = t121 * qJD(4) + t124 * t74;
t218 = t57 * t156;
t186 = pkin(1) * qJDD(1);
t217 = t162 * t104;
t184 = qJD(5) * t57;
t21 = -t124 * qJDD(4) - t121 * t41 + t184;
t216 = t106 * pkin(4) + t104 * pkin(7);
t209 = g(3) * t104;
t132 = -t150 * t106 - t209;
t38 = qJDD(5) + t42;
t76 = t116 * t177 - t118 * t161;
t146 = t38 * t83 - t68 * t76;
t176 = qJD(5) * t121;
t166 = t83 * t176;
t215 = -t124 * t146 + t68 * t166;
t214 = t74 ^ 2;
t9 = -t121 * t23 + t124 * t28;
t213 = t9 * t68;
t212 = pkin(6) + t94;
t123 = sin(qJ(1));
t203 = t123 * pkin(1);
t173 = t124 * qJD(4);
t55 = t121 * t74 - t173;
t202 = t55 * t72;
t201 = t55 * t76;
t200 = t57 * t55;
t199 = t57 * t74;
t198 = t57 * t76;
t197 = t74 * t55;
t196 = t74 * t72;
t175 = qJD(5) * t124;
t195 = -t121 * t21 - t55 * t175;
t188 = t21 * t124;
t194 = t124 * t201 - t83 * t188;
t134 = t164 - t179;
t20 = -qJD(5) * t173 - t121 * qJDD(4) + t124 * t41 + t74 * t176;
t193 = t134 * t20 + t57 * t77;
t192 = -t83 * t42 + t76 * t72;
t190 = t121 * t38;
t189 = t122 * t54;
t125 = cos(qJ(1));
t111 = t125 * pkin(1);
t187 = t107 * t97 + t111;
t185 = qJD(5) * t55;
t183 = t105 * t121;
t182 = t105 * t124;
t181 = t107 * t121;
t180 = t107 * t124;
t112 = t116 ^ 2;
t113 = t118 ^ 2;
t178 = t112 + t113;
t165 = t83 * t175;
t158 = -t20 + t185;
t155 = t124 * t68;
t153 = t57 * t165;
t24 = t211 * t53 - t189;
t22 = -qJD(4) * pkin(4) - t24;
t152 = -t22 * t76 + t6 * t83;
t1 = qJD(5) * t9 + t121 * t14 + t124 * t5;
t151 = t1 - t213;
t148 = g(1) * t123 - g(2) * t125;
t147 = t134 * t21 - t77 * t55;
t145 = t134 * t41 + t74 * t77;
t143 = -t10 * t124 + t121 * t9;
t142 = t10 * t121 + t124 * t9;
t120 = -pkin(6) - qJ(3);
t141 = -t107 * t120 - t203;
t58 = -t116 * t80 + t101;
t140 = -t58 * t116 + t59 * t118;
t139 = (-t116 * t87 + t103) * t116 - t65 * t118;
t32 = -pkin(4) * t134 - t83 * pkin(7) + t86;
t78 = t212 * t116;
t79 = t212 * t118;
t35 = -t122 * t78 + t211 * t79;
t16 = -t121 * t35 + t124 * t32;
t17 = t121 * t32 + t124 * t35;
t138 = t124 * t38 - t72 * t156 - t68 * t176;
t137 = -qJD(5) * t28 + t209 - t5;
t136 = -t122 * t79 - t211 * t78;
t135 = -pkin(7) * t38 + t68 * t22;
t133 = -t174 + t221;
t129 = -t146 * t121 - t68 * t165;
t128 = -qJD(5) * t142 + t1 * t124 - t2 * t121;
t71 = t72 ^ 2;
t63 = t106 * t180 + t183;
t62 = -t106 * t181 + t182;
t61 = -t106 * t182 + t181;
t60 = t106 * t183 + t180;
t44 = -t77 * qJD(4) + qJDD(4) * t134;
t43 = -t76 * qJD(4) + t83 * qJDD(4);
t40 = t77 * pkin(4) + t76 * pkin(7);
t39 = t74 * pkin(4) + t72 * pkin(7);
t27 = qJD(3) * t83 + qJD(4) * t35;
t26 = qJD(3) * t134 + qJD(4) * t136;
t12 = t121 * t39 + t124 * t24;
t11 = -t121 * t24 + t124 * t39;
t4 = -qJD(5) * t17 - t121 * t26 + t124 * t40;
t3 = qJD(5) * t16 + t121 * t40 + t124 * t26;
t15 = [0, 0, 0, 0, 0, qJDD(1), t148, g(1) * t125 + g(2) * t123, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t119 * t186 - t162, -0.2e1 * t117 * t186 + t150, 0, (t148 + (t117 ^ 2 + t119 ^ 2) * t186) * pkin(1), t112 * qJDD(1), 0.2e1 * t116 * t170, 0, t113 * qJDD(1), 0, 0, t133 * t118, -t133 * t116, t80 * t178 + t140 - t150, t85 * t98 - g(1) * (-t105 * pkin(2) + t107 * qJ(3) - t203) - g(2) * (t107 * pkin(2) + t105 * qJ(3) + t111) + t140 * t94 - t139 * qJD(3), -t41 * t83 - t74 * t76, -t145 + t192, t43, -t134 * t42 + t72 * t77, t44, 0, -t27 * qJD(4) + qJDD(4) * t136 - t106 * t162 - t134 * t67 + t86 * t42 + t70 * t77, -t26 * qJD(4) - t35 * qJDD(4) - t86 * t41 + t67 * t83 - t70 * t76 + t217, t134 * t7 + t136 * t41 + t24 * t76 - t25 * t77 - t26 * t72 + t27 * t74 - t35 * t42 - t8 * t83 - t150, t7 * t35 + t25 * t26 + t8 * t136 - t24 * t27 + t67 * t86 - g(1) * (-t105 * t97 + t141) - g(2) * (-t105 * t120 + t187), -t57 * t166 + (-t20 * t83 - t198) * t124, -t153 + (t198 + (t20 + t185) * t83) * t121 + t194, t193 - t215, t55 * t165 + (t21 * t83 - t201) * t121, t129 + t147, -t134 * t38 + t68 * t77, -g(1) * t61 - g(2) * t63 + t121 * t152 - t134 * t2 - t136 * t21 + t16 * t38 + t165 * t22 + t27 * t55 + t4 * t68 + t9 * t77, -g(1) * t60 - g(2) * t62 + t1 * t134 - t10 * t77 + t124 * t152 + t136 * t20 - t166 * t22 - t17 * t38 + t27 * t57 - t3 * t68, t16 * t20 - t17 * t21 - t3 * t55 - t4 * t57 + t142 * t76 - t217 + (t143 * qJD(5) - t1 * t121 - t2 * t124) * t83, t1 * t17 + t10 * t3 + t2 * t16 + t9 * t4 - t6 * t136 + t22 * t27 - g(1) * t141 - g(2) * (t216 * t107 + t187) + (-g(1) * (-t97 - t216) + g(2) * t120) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t116 + t58 * t118 - g(3), 0, 0, 0, 0, 0, 0, t44, -t43, t145 + t192, t134 * t8 - t24 * t77 - t25 * t76 + t7 * t83 - g(3), 0, 0, 0, 0, 0, 0, t129 - t147, t193 + t215, t153 + (t158 * t83 - t198) * t121 + t194, t128 * t83 - t134 * t6 + t143 * t76 + t22 * t77 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, t171, -t178 * qJD(1) ^ 2, qJD(1) * t139 - t221, 0, 0, 0, 0, 0, 0, 0.2e1 * t74 * qJD(4) + t144, (-t72 - t163) * qJD(4) + t168, -t71 - t214, t24 * t74 + t25 * t72 + t162 + t67, 0, 0, 0, 0, 0, 0, t138 - t197, -t124 * t68 ^ 2 - t190 - t199, (t20 - t202) * t124 + t218 + t195, t151 * t121 + t219 * t124 - t22 * t74 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, -t71 + t214, (t72 - t163) * qJD(4) + t168, -t196, -t144, qJDD(4), -t70 * t74 + t131 - t160, t70 * t72 + (t24 + t189) * qJD(4) + t169 - t132, 0, 0, -t20 * t121 + t155 * t57, (-t20 - t202) * t124 - t218 + t195, t155 * t68 + t190 - t199, t156 * t55 - t188, t138 + t197, -t68 * t74, -pkin(4) * t21 - t11 * t68 + t135 * t121 + t220 * t124 - t25 * t55 - t9 * t74, pkin(4) * t20 + t10 * t74 + t12 * t68 - t220 * t121 + t135 * t124 - t25 * t57, t11 * t57 + t12 * t55 + ((-t21 + t184) * pkin(7) + t151) * t124 + (pkin(7) * t158 - t219) * t121 + t132, -t10 * t12 - t9 * t11 - t22 * t25 + t130 * pkin(4) + (t128 + t132) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t55 ^ 2 + t57 ^ 2, t55 * t68 - t20, -t200, t57 * t68 - t21, t38, -g(1) * t62 + g(2) * t60 + t121 * t137 - t175 * t23 - t22 * t57 + t13 + t207, g(1) * t63 - g(2) * t61 + t22 * t55 + t213 + (qJD(5) * t23 - t14) * t121 + t137 * t124, 0, 0;];
tau_reg = t15;
