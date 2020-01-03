% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:20
% DurationCPUTime: 2.15s
% Computational Cost: add. (2777->298), mult. (7095->446), div. (0->0), fcn. (4573->6), ass. (0->155)
t122 = cos(qJ(3));
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t175 = t122 * t123;
t129 = pkin(3) * t120 - pkin(7) * t175;
t201 = -pkin(7) - pkin(6);
t152 = qJD(3) * t201;
t119 = sin(qJ(3));
t170 = qJD(1) * t120;
t151 = t119 * t170;
t133 = pkin(2) * t120 - pkin(6) * t123;
t84 = t133 * qJD(1);
t52 = pkin(5) * t151 + t122 * t84;
t212 = t129 * qJD(1) - t122 * t152 + t52;
t176 = t120 * t122;
t177 = t119 * t123;
t69 = t119 * t84;
t211 = t69 + (-pkin(5) * t176 - pkin(7) * t177) * qJD(1) - t119 * t152;
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t161 = t122 * qJD(2);
t79 = t151 - t161;
t150 = t122 * t170;
t169 = qJD(2) * t119;
t81 = t150 + t169;
t132 = t118 * t79 - t121 * t81;
t35 = t118 * t81 + t121 * t79;
t198 = t35 * t132;
t159 = qJD(1) * qJD(2);
t145 = t123 * t159;
t210 = qJD(2) * qJD(3) + t145;
t209 = t132 ^ 2 - t35 ^ 2;
t108 = t120 * t159;
t137 = pkin(5) * t108;
t89 = -pkin(2) * t123 - pkin(6) * t120 - pkin(1);
t73 = t89 * qJD(1);
t160 = t123 * qJD(1);
t114 = pkin(5) * t160;
t95 = qJD(2) * pkin(6) + t114;
t42 = t119 * t73 + t122 * t95;
t87 = t133 * qJD(2);
t74 = qJD(1) * t87;
t21 = -t42 * qJD(3) + t119 * t137 + t122 * t74;
t165 = qJD(3) * t120;
t144 = qJD(1) * t165;
t51 = t119 * t144 - t210 * t122;
t10 = pkin(3) * t108 + pkin(7) * t51 + t21;
t154 = t210 * t119 + t122 * t144;
t164 = qJD(3) * t122;
t166 = qJD(3) * t119;
t20 = t119 * t74 - t122 * t137 + t73 * t164 - t95 * t166;
t13 = -t154 * pkin(7) + t20;
t31 = -pkin(7) * t79 + t42;
t185 = t121 * t31;
t106 = -qJD(3) + t160;
t41 = -t119 * t95 + t122 * t73;
t30 = -pkin(7) * t81 + t41;
t26 = -pkin(3) * t106 + t30;
t6 = t118 * t26 + t185;
t2 = -t6 * qJD(4) + t121 * t10 - t118 * t13;
t113 = pkin(5) * t170;
t192 = qJD(2) * pkin(2);
t94 = t113 - t192;
t54 = pkin(3) * t79 + t94;
t208 = t132 * t54 + t2;
t101 = -qJD(4) + t106;
t162 = qJD(4) * t121;
t163 = qJD(4) * t118;
t11 = t118 * t154 + t121 * t51 + t79 * t162 + t81 * t163;
t207 = -t101 * t35 - t11;
t1 = (qJD(4) * t26 + t13) * t121 + t10 * t118 - t31 * t163;
t206 = t35 * t54 - t1;
t126 = t132 * qJD(4) + t118 * t51 - t121 * t154;
t205 = t101 * t132 + t126;
t204 = -0.2e1 * t159;
t107 = pkin(5) * t175;
t57 = t119 * t89 + t107;
t203 = t106 * t41 + t20;
t202 = t106 * t42 - t21;
t157 = qJD(3) + qJD(4);
t200 = pkin(3) * t119;
t199 = pkin(5) * t119;
t197 = t81 * t79;
t96 = t201 * t119;
t97 = t201 * t122;
t50 = t118 * t96 - t121 * t97;
t196 = t50 * qJD(4) - t211 * t118 + t212 * t121;
t49 = t118 * t97 + t121 * t96;
t195 = -t49 * qJD(4) + t212 * t118 + t211 * t121;
t179 = t118 * t119;
t82 = -t121 * t122 + t179;
t194 = -t121 * t164 - t122 * t162 + t157 * t179 - t82 * t160;
t83 = t118 * t122 + t119 * t121;
t44 = t157 * t83;
t193 = -t83 * t160 + t44;
t189 = t106 * t79;
t188 = t106 * t81;
t187 = t118 * t31;
t186 = t119 * t94;
t184 = t122 * t94;
t183 = t51 * t119;
t168 = qJD(2) * t120;
t182 = t122 * t87 + t168 * t199;
t181 = t106 * t119;
t180 = t106 * t122;
t178 = t119 * t120;
t125 = qJD(1) ^ 2;
t174 = t123 * t125;
t124 = qJD(2) ^ 2;
t173 = t124 * t120;
t172 = t124 * t123;
t116 = t120 ^ 2;
t171 = -t123 ^ 2 + t116;
t167 = qJD(2) * t123;
t156 = pkin(5) * t177;
t155 = pkin(5) * t167;
t153 = t120 * t174;
t149 = t119 * t167;
t148 = t123 * t161;
t147 = t119 * t165;
t146 = t120 * t164;
t140 = pkin(1) * t204;
t139 = -t81 + t169;
t138 = t79 + t161;
t136 = t123 * t108;
t135 = pkin(3) * t166 - t160 * t200 - t114;
t134 = t154 * t122;
t78 = t122 * t89;
t40 = -pkin(7) * t176 + t78 + (-pkin(3) - t199) * t123;
t45 = -pkin(7) * t178 + t57;
t17 = -t118 * t45 + t121 * t40;
t18 = t118 * t40 + t121 * t45;
t131 = -t119 * t42 - t41 * t122;
t130 = qJD(1) * t116 - t106 * t123;
t128 = t146 + t149;
t28 = t119 * t87 + t89 * t164 + (-t120 * t161 - t123 * t166) * pkin(5);
t112 = -pkin(3) * t122 - pkin(2);
t88 = (pkin(5) + t200) * t120;
t65 = t82 * t120;
t64 = t83 * t120;
t56 = t78 - t156;
t55 = t128 * pkin(3) + t155;
t53 = -pkin(5) * t150 + t69;
t33 = t154 * pkin(3) + pkin(5) * t145;
t29 = -t57 * qJD(3) + t182;
t23 = -t163 * t178 + (t157 * t176 + t149) * t121 + (t148 - t147) * t118;
t22 = t118 * t149 + t44 * t120 - t121 * t148;
t19 = -t128 * pkin(7) + t28;
t16 = t129 * qJD(2) + (-t107 + (pkin(7) * t120 - t89) * t119) * qJD(3) + t182;
t8 = t121 * t30 - t187;
t7 = -t118 * t30 - t185;
t5 = t121 * t26 - t187;
t4 = -t18 * qJD(4) - t118 * t19 + t121 * t16;
t3 = t17 * qJD(4) + t118 * t16 + t121 * t19;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, t171 * t204, t172, -0.2e1 * t136, -t173, 0, -pkin(5) * t172 + t120 * t140, pkin(5) * t173 + t123 * t140, 0, 0, t81 * t148 + (-t51 * t122 - t81 * t166) * t120, (-t119 * t81 - t122 * t79) * t167 + (-t134 + t183 + (t119 * t79 - t122 * t81) * qJD(3)) * t120, t106 * t147 + t123 * t51 + (t120 * t81 + t130 * t122) * qJD(2), t79 * t146 + (t154 * t120 + t79 * t167) * t119, t106 * t146 + t154 * t123 + (-t130 * t119 - t79 * t120) * qJD(2), (-t106 - t160) * t168, -t29 * t106 - t21 * t123 + (pkin(5) * t154 + t94 * t164) * t120 + ((pkin(5) * t79 + t186) * t123 + (t41 + (t56 + t156) * qJD(1)) * t120) * qJD(2), t106 * t28 + t123 * t20 + (-pkin(5) * t51 - t94 * t166) * t120 + ((pkin(5) * t81 + t184) * t123 + (-t42 + (-t57 + t107) * qJD(1)) * t120) * qJD(2), -t28 * t79 - t57 * t154 - t29 * t81 + t56 * t51 + t131 * t167 + (-t20 * t119 - t21 * t122 + (t119 * t41 - t122 * t42) * qJD(3)) * t120, t20 * t57 + t21 * t56 + t28 * t42 + t29 * t41 + (t94 + t113) * t155, t11 * t65 + t132 * t22, t11 * t64 - t126 * t65 + t132 * t23 + t22 * t35, t101 * t22 + t11 * t123 + (-qJD(1) * t65 - t132) * t168, -t126 * t64 + t23 * t35, t101 * t23 - t126 * t123 + (-qJD(1) * t64 - t35) * t168, (-t101 - t160) * t168, -t101 * t4 - t126 * t88 - t123 * t2 + t23 * t54 + t33 * t64 + t35 * t55 + (qJD(1) * t17 + t5) * t168, t1 * t123 + t101 * t3 - t11 * t88 - t22 * t54 - t33 * t65 - t132 * t55 + (-qJD(1) * t18 - t6) * t168, -t1 * t64 + t11 * t17 + t126 * t18 + t132 * t4 + t2 * t65 + t22 * t5 - t23 * t6 - t3 * t35, t1 * t18 + t17 * t2 + t3 * t6 + t33 * t88 + t4 * t5 + t54 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t171 * t125, 0, t153, 0, 0, t125 * pkin(1) * t120, pkin(1) * t174, 0, 0, -t81 * t180 - t183, (-t51 + t189) * t122 + (-t154 + t188) * t119, -t106 * t164 + (t106 * t175 + t139 * t120) * qJD(1), -t79 * t181 - t134, t106 * t166 + (-t106 * t177 + t138 * t120) * qJD(1), t106 * t170, -pkin(2) * t154 + t52 * t106 + (pkin(6) * t180 + t186) * qJD(3) + ((-pkin(6) * t169 - t41) * t120 + (-t138 * pkin(5) - t186) * t123) * qJD(1), pkin(2) * t51 - t106 * t53 + (-pkin(6) * t181 + t184) * qJD(3) + ((-pkin(6) * t161 + t42) * t120 + (t139 * pkin(5) - t184) * t123) * qJD(1), t52 * t81 + t53 * t79 + ((qJD(3) * t81 - t154) * pkin(6) + t203) * t122 + ((qJD(3) * t79 - t51) * pkin(6) + t202) * t119, -t41 * t52 - t42 * t53 + (-t94 - t192) * t114 + (t131 * qJD(3) - t21 * t119 + t20 * t122) * pkin(6), -t11 * t83 + t132 * t194, t11 * t82 + t126 * t83 + t132 * t193 + t194 * t35, t194 * t101 + (qJD(2) * t83 + t132) * t170, -t126 * t82 + t193 * t35, t193 * t101 + (-qJD(2) * t82 + t35) * t170, t101 * t170, -t112 * t126 + t33 * t82 + t193 * t54 + t135 * t35 + t196 * t101 + (qJD(2) * t49 - t5) * t170, -t11 * t112 + t33 * t83 - t194 * t54 - t135 * t132 - t195 * t101 + (-qJD(2) * t50 + t6) * t170, -t1 * t82 + t11 * t49 + t126 * t50 - t132 * t196 - t193 * t6 + t194 * t5 + t195 * t35 - t2 * t83, t1 * t50 + t112 * t33 + t135 * t54 - t195 * t6 - t196 * t5 + t2 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t79 ^ 2 + t81 ^ 2, -t51 - t189, -t197, -t154 - t188, t108, -t81 * t94 - t202, t79 * t94 - t203, 0, 0, -t198, t209, t207, t198, t205, t108, t101 * t7 + (t101 * t163 + t108 * t121 - t35 * t81) * pkin(3) + t208, -t101 * t8 + (t101 * t162 - t108 * t118 + t132 * t81) * pkin(3) + t206, -t132 * t6 + t35 * t8 - t35 * t5 - t132 * t7 + (t11 * t121 + t118 * t126 + (-t118 * t132 - t121 * t35) * qJD(4)) * pkin(3), -t5 * t7 - t6 * t8 + (t1 * t118 + t121 * t2 - t54 * t81 + (-t118 * t5 + t121 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, t209, t207, t198, t205, t108, -t101 * t6 + t208, -t101 * t5 + t206, 0, 0;];
tauc_reg = t9;
