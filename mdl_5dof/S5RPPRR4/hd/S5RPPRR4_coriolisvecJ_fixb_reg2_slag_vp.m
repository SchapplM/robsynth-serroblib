% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:55
% EndTime: 2019-12-05 17:45:08
% DurationCPUTime: 2.28s
% Computational Cost: add. (3849->240), mult. (10406->356), div. (0->0), fcn. (7843->8), ass. (0->135)
t128 = cos(pkin(8));
t155 = t128 * qJD(1);
t190 = qJD(4) - t155;
t114 = -qJD(5) - t190;
t129 = sin(qJ(5));
t131 = cos(qJ(5));
t126 = sin(pkin(8));
t125 = sin(pkin(9));
t127 = cos(pkin(9));
t130 = sin(qJ(4));
t132 = cos(qJ(4));
t106 = t132 * t125 + t130 * t127;
t137 = qJD(1) * t106;
t80 = t126 * t137;
t163 = qJD(1) * t126;
t149 = t125 * t163;
t166 = t132 * t127;
t151 = t126 * t166;
t83 = qJD(1) * t151 - t130 * t149;
t37 = t129 * t83 + t131 * t80;
t173 = t37 * t114;
t157 = qJD(5) * t131;
t158 = qJD(5) * t129;
t95 = t106 * t126;
t86 = qJD(4) * t95;
t72 = qJD(1) * t86;
t73 = t83 * qJD(4);
t24 = t129 * t73 + t131 * t72 + t80 * t157 + t83 * t158;
t189 = -t24 - t173;
t138 = -t129 * t80 + t131 * t83;
t181 = t138 ^ 2;
t182 = t37 ^ 2;
t188 = t181 - t182;
t115 = qJ(2) * t163 + qJD(3);
t98 = pkin(3) * t149 + t115;
t48 = t80 * pkin(4) + t98;
t187 = t48 * t37;
t180 = t37 * t138;
t105 = -t130 * t125 + t166;
t170 = t190 * t105;
t169 = -t106 * qJD(4) + t128 * t137;
t174 = t138 * t114;
t25 = t138 * qJD(5) - t129 * t72 + t131 * t73;
t186 = -t25 - t174;
t185 = t48 * t138;
t159 = qJD(4) * t132;
t160 = qJD(4) * t130;
t136 = -t127 * t126 * pkin(6) + (-qJ(2) * t125 - pkin(3)) * t128;
t109 = -t128 * pkin(2) - t126 * qJ(3) - pkin(1);
t97 = t109 * qJD(1) + qJD(2);
t89 = t127 * t97;
t45 = t136 * qJD(1) + t89;
t150 = qJ(2) * t155;
t60 = t125 * t97 + t127 * t150;
t50 = -pkin(6) * t149 + t60;
t161 = qJD(3) * t126;
t162 = qJD(2) * t128;
t99 = -t125 * t162 - t127 * t161;
t93 = t99 * qJD(1);
t100 = -t125 * t161 + t127 * t162;
t94 = t100 * qJD(1);
t14 = t130 * t93 + t132 * t94 + t45 * t159 - t50 * t160;
t12 = -t73 * pkin(7) + t14;
t31 = t130 * t45 + t132 * t50;
t15 = -t31 * qJD(4) - t130 * t94 + t132 * t93;
t13 = t72 * pkin(7) + t15;
t21 = -t80 * pkin(7) + t31;
t145 = t129 * t13 - t21 * t158;
t30 = -t130 * t50 + t132 * t45;
t20 = -t83 * pkin(7) + t30;
t18 = pkin(4) * t190 + t20;
t1 = (qJD(5) * t18 + t12) * t131 + t145;
t184 = t83 ^ 2;
t183 = pkin(4) * t83;
t179 = t83 * t80;
t56 = t131 * t105 - t129 * t106;
t178 = t56 * qJD(5) + t169 * t129 + t170 * t131;
t57 = t129 * t105 + t131 * t106;
t177 = -t57 * qJD(5) - t170 * t129 + t169 * t131;
t104 = t127 * t109;
t53 = t104 + t136;
t168 = qJ(2) * t128;
t165 = t125 * t109 + t127 * t168;
t167 = t125 * t126;
t61 = -pkin(6) * t167 + t165;
t35 = t130 * t53 + t132 * t61;
t176 = t129 * t21;
t175 = t131 * t21;
t172 = t80 * t190;
t171 = t83 * t190;
t107 = pkin(3) * t167 + t126 * qJ(2);
t123 = t126 ^ 2;
t124 = t128 ^ 2;
t164 = t123 + t124;
t156 = t126 * qJD(2);
t154 = qJD(1) * qJD(2);
t153 = 0.2e1 * qJD(2) * t123;
t133 = qJD(1) ^ 2;
t152 = t126 * t128 * t133;
t148 = pkin(4) * t114 - t18;
t120 = t126 * t154;
t58 = t73 * pkin(4) + t120;
t147 = qJ(2) * t154;
t146 = -t129 * t12 + t131 * t13;
t34 = -t130 * t61 + t132 * t53;
t142 = t164 * t133;
t139 = -t94 * t125 - t93 * t127;
t6 = t129 * t18 + t175;
t96 = t105 * t126;
t28 = -t128 * pkin(4) - t96 * pkin(7) + t34;
t29 = -t95 * pkin(7) + t35;
t9 = -t129 * t29 + t131 * t28;
t10 = t129 * t28 + t131 * t29;
t47 = -t129 * t95 + t131 * t96;
t22 = t132 * t100 + t130 * t99 + t53 * t159 - t61 * t160;
t135 = qJD(4) * t80;
t2 = -t6 * qJD(5) + t146;
t23 = -t35 * qJD(4) - t130 * t100 + t132 * t99;
t117 = t123 * t147;
t87 = qJD(4) * t151 - t160 * t167;
t78 = t80 ^ 2;
t62 = t87 * pkin(4) + t156;
t59 = -t125 * t150 + t89;
t55 = t95 * pkin(4) + t107;
t46 = t129 * t96 + t131 * t95;
t27 = t47 * qJD(5) - t129 * t86 + t131 * t87;
t26 = t129 * t87 + t131 * t86 + t95 * t157 + t96 * t158;
t17 = t86 * pkin(7) + t23;
t16 = -t87 * pkin(7) + t22;
t8 = t131 * t20 - t176;
t7 = -t129 * t20 - t175;
t5 = t131 * t18 - t176;
t4 = -t10 * qJD(5) - t129 * t16 + t131 * t17;
t3 = t9 * qJD(5) + t129 * t17 + t131 * t16;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t164 * t154, 0.2e1 * t124 * t147 + 0.2e1 * t117, 0, 0, 0, 0, 0, 0, -t93 * t128 + (t125 * t153 - t128 * t99) * qJD(1), t94 * t128 + (t100 * t128 + t127 * t153) * qJD(1), ((-t100 * t125 - t127 * t99) * qJD(1) + t139) * t126, t94 * t165 + t60 * t100 + t93 * (-t125 * t168 + t104) + t59 * t99 + t117 + t115 * t156, -t72 * t96 - t83 * t86, t72 * t95 - t96 * t73 + t86 * t80 - t83 * t87, t72 * t128 - t190 * t86, t73 * t95 + t80 * t87, t73 * t128 - t190 * t87, 0, t107 * t73 + t23 * t190 - t15 * t128 + t98 * t87 + (qJD(1) * t95 + t80) * t156, -t107 * t72 - t22 * t190 + t14 * t128 - t98 * t86 + (qJD(1) * t96 + t83) * t156, -t14 * t95 - t15 * t96 - t22 * t80 - t23 * t83 + t30 * t86 - t31 * t87 + t34 * t72 - t35 * t73, t14 * t35 + t15 * t34 + t31 * t22 + t30 * t23 + (qJD(1) * t107 + t98) * t156, -t138 * t26 - t24 * t47, -t138 * t27 + t24 * t46 - t47 * t25 + t26 * t37, t26 * t114 + t24 * t128, t25 * t46 + t37 * t27, t27 * t114 + t25 * t128, 0, -t4 * t114 - t2 * t128 + t55 * t25 + t48 * t27 + t62 * t37 + t58 * t46, t1 * t128 + t3 * t114 + t138 * t62 - t55 * t24 - t48 * t26 + t58 * t47, -t1 * t46 - t10 * t25 - t138 * t4 - t2 * t47 + t9 * t24 + t5 * t26 - t6 * t27 - t3 * t37, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4 + t48 * t62 + t58 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -qJ(2) * t142, 0, 0, 0, 0, 0, 0, -t125 * t142, -t127 * t142, 0, (-t115 * t126 + (t125 * t59 - t127 * t60) * t128) * qJD(1) - t139, 0, 0, 0, 0, 0, 0, -t80 * t163 + t169 * t190, -t83 * t163 - t170 * t190, t105 * t72 - t106 * t73 - t169 * t83 - t170 * t80, t15 * t105 + t14 * t106 - t98 * t163 + t169 * t30 + t170 * t31, 0, 0, 0, 0, 0, 0, -t177 * t114 - t37 * t163, t178 * t114 - t138 * t163, -t138 * t177 - t178 * t37 + t56 * t24 - t57 * t25, t1 * t57 - t48 * t163 + t177 * t5 + t178 * t6 + t2 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 * t152, t125 * t152, (-t125 ^ 2 - t127 ^ 2) * t133 * t123, t120 + (t125 * t60 + t127 * t59) * t163, 0, 0, 0, 0, 0, 0, t73 + t171, -t135 - t172, -t78 - t184, t30 * t83 + t31 * t80 + t120, 0, 0, 0, 0, 0, 0, t25 - t174, -t24 + t173, -t181 - t182, t138 * t5 + t37 * t6 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t78 + t184, -t135 + t172, -t179, -t73 + t171, 0, t190 * t31 - t98 * t83 + t15, t190 * t30 + t98 * t80 - t14, 0, 0, t180, t188, t189, -t180, t186, 0, -t37 * t183 + t7 * t114 - t185 + (t148 * t129 - t175) * qJD(5) + t146, -t138 * t183 - t8 * t114 + t187 + (qJD(5) * t148 - t12) * t131 - t145, t6 * t138 + t8 * t37 - t5 * t37 + t7 * t138 + (-t129 * t25 + t131 * t24 + (t129 * t138 - t131 * t37) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t1 * t129 + t131 * t2 - t48 * t83 + (-t129 * t5 + t131 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t188, t189, -t180, t186, 0, -t6 * t114 - t185 + t2, -t5 * t114 - t1 + t187, 0, 0;];
tauc_reg = t11;
