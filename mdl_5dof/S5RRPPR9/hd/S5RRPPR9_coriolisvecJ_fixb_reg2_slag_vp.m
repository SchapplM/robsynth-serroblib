% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:48
% EndTime: 2019-12-31 19:41:52
% DurationCPUTime: 1.38s
% Computational Cost: add. (1611->258), mult. (3542->352), div. (0->0), fcn. (1791->4), ass. (0->145)
t86 = cos(qJ(5));
t135 = t86 * qJD(2);
t87 = cos(qJ(2));
t144 = qJD(1) * t87;
t84 = sin(qJ(5));
t37 = t84 * t144 - t135;
t85 = sin(qJ(2));
t137 = t85 * qJD(1);
t57 = qJD(5) + t137;
t101 = t37 * t57;
t139 = qJD(5) * t87;
t124 = t84 * t139;
t94 = t85 * t135 + t124;
t18 = t94 * qJD(1) - qJD(5) * t135;
t163 = t18 - t101;
t138 = t84 * qJD(2);
t38 = t86 * t144 + t138;
t154 = t38 * t57;
t129 = qJD(1) * qJD(2);
t119 = t85 * t129;
t19 = t38 * qJD(5) - t84 * t119;
t95 = t19 - t154;
t148 = pkin(6) - qJ(4);
t118 = t87 * t129;
t69 = t85 * qJD(3);
t147 = qJ(3) * t118 + qJD(1) * t69;
t88 = -pkin(2) - pkin(3);
t78 = -pkin(7) + t88;
t96 = pkin(4) * t87 + t78 * t85;
t93 = t96 * qJD(2);
t10 = qJD(1) * t93 + t147;
t130 = qJ(4) * qJD(2);
t121 = t87 * t130;
t136 = t85 * qJD(4);
t56 = pkin(6) * t118;
t24 = t56 + (-t121 - t136) * qJD(1);
t109 = t85 * pkin(4) + t87 * pkin(7);
t132 = qJ(3) * qJD(1);
t34 = -qJD(1) * pkin(1) - pkin(2) * t144 - t85 * t132;
t23 = pkin(3) * t144 + qJD(4) - t34;
t14 = t109 * qJD(1) + t23;
t126 = pkin(6) * t137;
t106 = qJD(3) + t126;
t131 = qJ(4) * qJD(1);
t60 = t85 * t131;
t100 = t106 - t60;
t21 = t78 * qJD(2) + t100;
t5 = t86 * t14 - t84 * t21;
t1 = qJD(5) * t5 + t84 * t10 + t86 * t24;
t165 = -t5 * t57 + t1;
t6 = t84 * t14 + t86 * t21;
t2 = -qJD(5) * t6 + t86 * t10 - t84 * t24;
t164 = -t6 * t57 - t2;
t65 = pkin(6) * t144;
t41 = -t87 * t131 + t65;
t79 = qJD(2) * qJ(3);
t32 = -t41 - t79;
t120 = qJD(2) * t88;
t53 = qJ(4) * t119;
t77 = qJD(2) * qJD(3);
t143 = qJD(2) * t85;
t97 = pkin(6) * t143 + t87 * qJD(4);
t17 = t97 * qJD(1) - t53 - t77;
t48 = t148 * t87;
t160 = t17 * t48;
t159 = t17 * t84;
t158 = t17 * t86;
t157 = t18 * t84;
t156 = t19 * t86;
t155 = t37 * t38;
t153 = t57 * t85;
t152 = t57 * t86;
t81 = t87 ^ 2;
t90 = qJD(1) ^ 2;
t151 = t81 * t90;
t150 = t87 * t90;
t89 = qJD(2) ^ 2;
t149 = t89 * t85;
t71 = t89 * t87;
t146 = t87 * t79 + t69;
t145 = qJD(2) * pkin(2);
t142 = qJD(2) * t87;
t141 = qJD(5) * t84;
t140 = qJD(5) * t86;
t39 = -t60 + t126;
t134 = qJD(3) + t39;
t133 = -qJD(4) - t23;
t44 = -t87 * pkin(2) - t85 * qJ(3) - pkin(1);
t128 = t84 * t153;
t127 = t85 * t152;
t125 = t57 * t141;
t123 = t57 * t140;
t122 = t86 * t139;
t35 = t87 * pkin(3) - t44;
t43 = t106 - t145;
t117 = -t43 - t145;
t116 = t133 * t85;
t111 = t85 * t120;
t16 = qJD(1) * t111 + t147;
t22 = t111 + t146;
t115 = qJD(1) * t22 + t16;
t114 = qJD(1) * t35 + t23;
t113 = qJD(1) * t44 + t34;
t112 = -0.2e1 * pkin(1) * t129;
t110 = t85 * t118;
t108 = -t5 * t86 - t6 * t84;
t107 = -t5 * t84 + t6 * t86;
t20 = t109 + t35;
t47 = t148 * t85;
t11 = t86 * t20 - t84 * t47;
t12 = t84 * t20 + t86 * t47;
t104 = qJD(1) * t81 - t153;
t25 = pkin(2) * t119 - t147;
t31 = pkin(2) * t143 - t146;
t99 = -pkin(6) * t89 - qJD(1) * t31 - t25;
t27 = qJD(2) * pkin(4) - t32;
t98 = -t78 * t142 - t27 * t85;
t92 = t107 * qJD(5) + t1 * t84 + t2 * t86;
t42 = -pkin(6) * t119 + t77;
t46 = t65 + t79;
t91 = t42 * t87 + (t43 * t87 + (-t46 + t65) * t85) * qJD(2);
t83 = qJ(3) + pkin(4);
t80 = t85 ^ 2;
t74 = 0.2e1 * t77;
t68 = t80 * t90;
t62 = t87 * t132;
t55 = t85 * t150;
t52 = -t68 - t89;
t50 = -0.2e1 * t110;
t49 = 0.2e1 * t110;
t45 = -t68 + t151;
t40 = pkin(2) * t137 - t62;
t36 = (-t80 + t81) * t129;
t33 = 0.2e1 * t36;
t30 = t148 * t142 - t136;
t29 = -t85 * t130 + t97;
t28 = t88 * t137 + t62;
t26 = t120 + t100;
t15 = t96 * qJD(1) + t62;
t13 = t93 + t146;
t8 = t84 * t15 + t86 * t41;
t7 = t86 * t15 - t84 * t41;
t4 = -t12 * qJD(5) + t86 * t13 - t84 * t30;
t3 = t11 * qJD(5) + t84 * t13 + t86 * t30;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t33, t71, t50, -t149, 0, -pkin(6) * t71 + t85 * t112, pkin(6) * t149 + t87 * t112, 0, 0, t49, t71, -0.2e1 * t36, 0, t149, t50, t113 * t143 + t99 * t87, t91, -t113 * t142 + t99 * t85, t91 * pkin(6) + t25 * t44 + t34 * t31, t50, t33, -t149, t49, t71, 0, t115 * t85 + (t114 * t87 - t29) * qJD(2), -t115 * t87 + (t114 * t85 + t30) * qJD(2), t17 * t87 - t24 * t85 + (-t26 * t87 - t32 * t85) * qJD(2) + (t29 * t87 - t30 * t85 + (-t47 * t87 + t48 * t85) * qJD(2)) * qJD(1), t16 * t35 + t23 * t22 + t24 * t47 + t26 * t30 + t32 * t29 - t160, -t18 * t86 * t87 - t94 * t38, (t37 * t86 + t38 * t84) * t143 + (t157 - t156 + (t37 * t84 - t38 * t86) * qJD(5)) * t87, t57 * t124 + t18 * t85 + (-t104 * t86 - t38 * t87) * qJD(2), t19 * t84 * t87 + (-t85 * t138 + t122) * t37, t57 * t122 + t19 * t85 + (t104 * t84 + t37 * t87) * qJD(2), (t57 + t137) * t142, -t48 * t19 + t29 * t37 + t4 * t57 + (t27 * t138 + t2) * t85 + (-t27 * t140 + t159 + (qJD(1) * t11 + t5) * qJD(2)) * t87, t48 * t18 + t29 * t38 - t3 * t57 + (t27 * t135 - t1) * t85 + (t27 * t141 + t158 + (-qJD(1) * t12 - t6) * qJD(2)) * t87, t108 * t143 - t11 * t18 + t12 * t19 + t3 * t37 + t4 * t38 + t87 * t92, t1 * t12 + t2 * t11 - t27 * t29 + t6 * t3 + t5 * t4 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t45, 0, t55, 0, 0, t90 * pkin(1) * t85, pkin(1) * t150, 0, 0, -t55, 0, t45, 0, 0, t55, (-t34 * t85 + t40 * t87) * qJD(1), ((t46 - t79) * t85 + (qJD(3) + t117) * t87) * qJD(1), t74 + (t34 * t87 + t40 * t85) * qJD(1), t42 * qJ(3) + t46 * qJD(3) - t34 * t40 + (t117 * t87 + t46 * t85) * qJD(1) * pkin(6), t55, -t45, 0, -t55, 0, 0, t39 * qJD(2) + t53 + t74 + (t133 * t87 + (-pkin(6) * qJD(2) - t28) * t85) * qJD(1), -t41 * qJD(2) + t56 + ((t28 - t130) * t87 + t116) * qJD(1), (-t134 + t26 - t120) * t144, -t17 * qJ(3) - t134 * t32 - t23 * t28 + t24 * t88 - t26 * t41, t38 * t152 - t157, (-t18 - t101) * t86 + (-t19 - t154) * t84, -t123 + (-t127 + (t38 - t138) * t87) * qJD(1), t84 * t101 - t156, t125 + (t128 + (-t37 - t135) * t87) * qJD(1), -t57 * t144, -t158 - t83 * t19 - t7 * t57 - t134 * t37 + (-t78 * t152 - t27 * t84) * qJD(5) + (-t5 * t87 + t98 * t84) * qJD(1), t159 + t83 * t18 + t8 * t57 - t134 * t38 + (t57 * t78 * t84 - t27 * t86) * qJD(5) + (t6 * t87 + t98 * t86) * qJD(1), -t8 * t37 - t7 * t38 + (t5 * t137 + t19 * t78 - t1 + (-t38 * t78 + t5) * qJD(5)) * t86 + (t6 * t137 + t18 * t78 + t2 + (-t37 * t78 + t6) * qJD(5)) * t84, -t17 * t83 - t5 * t7 - t6 * t8 + t134 * t27 + (qJD(5) * t108 + t1 * t86 - t2 * t84) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, t52, -t46 * qJD(2) + t34 * t137 + t56, 0, 0, 0, 0, 0, 0, t52, t55, 0, t32 * qJD(2) + t56 + (t116 - t121) * qJD(1), 0, 0, 0, 0, 0, 0, -t123 + qJD(2) * t37 + (-t87 * t138 - t127) * qJD(1), t125 + qJD(2) * t38 + (-t87 * t135 + t128) * qJD(1), t163 * t84 + t95 * t86, -t27 * qJD(2) + t164 * t84 + t165 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t118, 0.2e1 * t119, -t68 - t151, (-t32 * t87 + (t26 + t120) * t85) * qJD(1) + t147, 0, 0, 0, 0, 0, 0, -t125 + (-t128 + (-t37 + t135) * t87) * qJD(1), -t123 + (-t127 + (-t38 - t138) * t87) * qJD(1), -t163 * t86 + t95 * t84, (t107 * t85 + t27 * t87) * qJD(1) + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t37 ^ 2 + t38 ^ 2, t163, -t155, t95, t118, t27 * t38 - t164, -t27 * t37 - t165, 0, 0;];
tauc_reg = t9;
