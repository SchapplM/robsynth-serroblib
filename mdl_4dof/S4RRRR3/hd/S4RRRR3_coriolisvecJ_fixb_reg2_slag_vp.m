% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRR3
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:41
% DurationCPUTime: 1.44s
% Computational Cost: add. (2569->220), mult. (6880->315), div. (0->0), fcn. (4737->6), ass. (0->134)
t124 = cos(qJ(4));
t126 = cos(qJ(2));
t181 = -pkin(6) - pkin(5);
t108 = t181 * t126;
t104 = qJD(1) * t108;
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t155 = qJD(3) * t125;
t156 = qJD(3) * t122;
t123 = sin(qJ(2));
t107 = t181 * t123;
t102 = qJD(1) * t107;
t170 = qJD(2) * pkin(2);
t93 = t102 + t170;
t148 = qJD(2) * t181;
t136 = qJD(1) * t148;
t94 = t123 * t136;
t95 = t126 * t136;
t139 = -t104 * t156 - t122 * t95 - t125 * t94 - t93 * t155;
t118 = qJD(2) + qJD(3);
t98 = t122 * t126 + t123 * t125;
t64 = t118 * t98;
t53 = t64 * qJD(1);
t13 = -pkin(7) * t53 - t139;
t121 = sin(qJ(4));
t140 = -t122 * t94 + t125 * t95;
t91 = t125 * t104;
t55 = t122 * t93 - t91;
t26 = -t55 * qJD(3) + t140;
t165 = t122 * t123;
t133 = t118 * t165;
t152 = qJD(1) * qJD(2);
t144 = t126 * t152;
t163 = t125 * t126;
t146 = qJD(1) * t163;
t159 = -qJD(3) * t146 - t125 * t144;
t52 = qJD(1) * t133 + t159;
t14 = t52 * pkin(7) + t26;
t154 = qJD(4) * t121;
t157 = qJD(1) * t123;
t147 = t122 * t157;
t84 = -t146 + t147;
t179 = pkin(7) * t84;
t34 = t55 - t179;
t142 = t121 * t14 - t34 * t154;
t87 = t122 * t104;
t54 = t125 * t93 + t87;
t86 = qJD(1) * t98;
t79 = t86 * pkin(7);
t33 = t54 - t79;
t29 = pkin(3) * t118 + t33;
t1 = (qJD(4) * t29 + t13) * t124 + t142;
t45 = t121 * t86 + t124 * t84;
t116 = -pkin(2) * t126 - pkin(1);
t106 = qJD(1) * t116;
t65 = t84 * pkin(3) + t106;
t176 = t65 * t45;
t186 = t176 - t1;
t132 = t121 * t84 - t124 * t86;
t178 = t132 * t65;
t167 = t124 * t34;
t11 = t121 * t29 + t167;
t143 = -t121 * t13 + t124 * t14;
t2 = -t11 * qJD(4) + t143;
t185 = t178 + t2;
t184 = t45 * t132;
t7 = t132 ^ 2 - t45 ^ 2;
t117 = qJD(4) + t118;
t153 = qJD(4) * t124;
t8 = t121 * t53 + t124 * t52 + t84 * t153 + t86 * t154;
t5 = t117 * t45 - t8;
t129 = t132 * qJD(4) + t121 * t52 - t124 * t53;
t6 = -t117 * t132 + t129;
t150 = t123 * t170;
t183 = 0.2e1 * t150;
t182 = -0.2e1 * t152;
t180 = pkin(3) * t86;
t173 = t86 * t84;
t115 = pkin(2) * t125 + pkin(3);
t164 = t122 * t124;
t61 = -t102 * t122 + t91;
t36 = t61 + t179;
t62 = t125 * t102 + t87;
t37 = -t79 + t62;
t172 = -t121 * t37 + t124 * t36 + t115 * t154 - (-t122 * t153 + (-t121 * t125 - t164) * qJD(3)) * pkin(2);
t166 = t121 * t122;
t171 = t121 * t36 + t124 * t37 - t115 * t153 - (-t122 * t154 + (t124 * t125 - t166) * qJD(3)) * pkin(2);
t169 = t106 * t86;
t168 = t121 * t34;
t67 = t122 * t107 - t125 * t108;
t128 = qJD(1) ^ 2;
t162 = t126 * t128;
t127 = qJD(2) ^ 2;
t161 = t127 * t123;
t160 = t127 * t126;
t158 = t123 ^ 2 - t126 ^ 2;
t151 = pkin(2) * t157;
t149 = t123 * t162;
t145 = -pkin(3) * t117 - t29;
t137 = pkin(1) * t182;
t66 = t125 * t107 + t108 * t122;
t135 = t123 * t144;
t10 = t124 * t29 - t168;
t134 = -t10 * t45 - t11 * t132;
t39 = -pkin(7) * t98 + t66;
t97 = -t163 + t165;
t40 = -pkin(7) * t97 + t67;
t23 = -t121 * t40 + t124 * t39;
t24 = t121 * t39 + t124 * t40;
t60 = -t121 * t97 + t124 * t98;
t131 = t106 * t84 + t139;
t103 = t123 * t148;
t105 = t126 * t148;
t27 = t125 * t103 + t122 * t105 + t107 * t155 + t108 * t156;
t28 = -t67 * qJD(3) - t122 * t103 + t125 * t105;
t81 = pkin(2) * t164 + t115 * t121;
t80 = -pkin(2) * t166 + t115 * t124;
t70 = pkin(3) * t97 + t116;
t68 = t151 + t180;
t63 = -qJD(2) * t163 - t126 * t155 + t133;
t59 = t121 * t98 + t124 * t97;
t49 = pkin(3) * t64 + t150;
t38 = pkin(3) * t53 + qJD(1) * t150;
t35 = -t84 ^ 2 + t86 ^ 2;
t31 = -t159 + (-t147 + t84) * t118;
t22 = t63 * pkin(7) + t28;
t21 = -pkin(7) * t64 + t27;
t20 = t60 * qJD(4) - t121 * t63 + t124 * t64;
t19 = t121 * t64 + t124 * t63 + t97 * t153 + t98 * t154;
t16 = t124 * t33 - t168;
t15 = -t121 * t33 - t167;
t4 = -t24 * qJD(4) - t121 * t21 + t124 * t22;
t3 = t23 * qJD(4) + t121 * t22 + t124 * t21;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t135, t158 * t182, t160, -0.2e1 * t135, -t161, 0, -pkin(5) * t160 + t123 * t137, pkin(5) * t161 + t126 * t137, 0, 0, -t52 * t98 - t63 * t86, t52 * t97 - t53 * t98 + t63 * t84 - t64 * t86, -t63 * t118, t53 * t97 + t64 * t84, -t64 * t118, 0, t106 * t64 + t116 * t53 + t28 * t118 + (qJD(1) * t97 + t84) * t150, -t106 * t63 - t116 * t52 - t27 * t118 + t86 * t183, t139 * t97 - t26 * t98 - t27 * t84 - t28 * t86 + t52 * t66 - t53 * t67 + t54 * t63 - t55 * t64, t106 * t183 - t139 * t67 + t26 * t66 + t55 * t27 + t54 * t28, t132 * t19 - t60 * t8, t129 * t60 + t132 * t20 + t19 * t45 + t59 * t8, -t19 * t117, -t129 * t59 + t20 * t45, -t20 * t117, 0, t117 * t4 - t129 * t70 + t20 * t65 + t38 * t59 + t45 * t49, -t117 * t3 - t132 * t49 - t19 * t65 + t38 * t60 - t70 * t8, -t1 * t59 + t10 * t19 - t11 * t20 + t129 * t24 + t132 * t4 - t2 * t60 + t23 * t8 - t3 * t45, t1 * t24 + t10 * t4 + t11 * t3 + t2 * t23 + t38 * t70 + t49 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t158 * t128, 0, t149, 0, 0, t128 * pkin(1) * t123, pkin(1) * t162, 0, 0, t173, t35, t31, -t173, 0, 0, -t84 * t151 - t169 - t118 * t61 + (t91 + (-pkin(2) * t118 - t93) * t122) * qJD(3) + t140, t118 * t62 + (-t118 * t155 - t86 * t157) * pkin(2) + t131, (t55 + t61) * t86 + (-t54 + t62) * t84 + (-t122 * t53 + t125 * t52 + (t122 * t86 - t125 * t84) * qJD(3)) * pkin(2), -t54 * t61 - t55 * t62 + (-t106 * t157 - t122 * t139 + t125 * t26 + (-t122 * t54 + t125 * t55) * qJD(3)) * pkin(2), -t184, t7, t5, t184, t6, 0, -t172 * t117 - t45 * t68 + t185, t171 * t117 + t132 * t68 + t186, t129 * t81 - t132 * t172 + t171 * t45 + t8 * t80 + t134, t1 * t81 - t172 * t10 - t171 * t11 + t2 * t80 - t65 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t35, t31, -t173, 0, 0, t55 * t118 - t169 + t26, t118 * t54 + t131, 0, 0, -t184, t7, t5, t184, t6, 0, -t45 * t180 - t117 * t15 + t178 + (t121 * t145 - t167) * qJD(4) + t143, t132 * t180 + t117 * t16 + t176 + (qJD(4) * t145 - t13) * t124 - t142, -t15 * t132 + t16 * t45 + (t121 * t129 + t124 * t8 + (-t121 * t132 - t124 * t45) * qJD(4)) * pkin(3) + t134, -t10 * t15 - t11 * t16 + (t1 * t121 + t124 * t2 - t65 * t86 + (-t10 * t121 + t11 * t124) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t7, t5, t184, t6, 0, t11 * t117 + t185, t10 * t117 + t186, 0, 0;];
tauc_reg = t9;
