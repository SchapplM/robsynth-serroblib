% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:42
% DurationCPUTime: 1.26s
% Computational Cost: add. (1528->213), mult. (3413->298), div. (0->0), fcn. (1868->6), ass. (0->137)
t65 = sin(pkin(8)) * pkin(1) + pkin(6);
t53 = t65 * qJD(1);
t81 = cos(qJ(3));
t71 = t81 * qJD(2);
t79 = sin(qJ(3));
t34 = t79 * t53 - t71;
t162 = -qJD(4) - t34;
t80 = cos(qJ(5));
t121 = t80 * qJD(3);
t132 = qJD(1) * t81;
t78 = sin(qJ(5));
t50 = -t78 * t132 + t121;
t114 = t79 * t121;
t127 = qJD(5) * t81;
t161 = t78 * t127 + t114;
t118 = qJD(1) * qJD(3);
t109 = t79 * t118;
t63 = pkin(3) * t109;
t122 = t79 * qJD(4);
t99 = pkin(7) * t79 - qJ(4) * t81;
t87 = t99 * qJD(3) - t122;
t14 = t87 * qJD(1) + t63;
t107 = pkin(4) * qJD(1) + t53;
t123 = t79 * qJD(2);
t15 = (t107 * t81 + t123) * qJD(3);
t120 = t107 * t79 + qJD(4) - t71;
t82 = -pkin(3) - pkin(7);
t13 = t82 * qJD(3) + t120;
t66 = -cos(pkin(8)) * pkin(1) - pkin(2);
t92 = -t79 * qJ(4) + t66;
t33 = t82 * t81 + t92;
t18 = t33 * qJD(1);
t5 = t80 * t13 - t78 * t18;
t1 = qJD(5) * t5 + t80 * t14 + t78 * t15;
t124 = t79 * qJD(1);
t64 = qJD(5) + t124;
t160 = -t5 * t64 + t1;
t6 = t78 * t13 + t80 * t18;
t2 = -qJD(5) * t6 - t78 * t14 + t80 * t15;
t159 = t6 * t64 + t2;
t125 = t78 * qJD(3);
t48 = t80 * t132 + t125;
t21 = t48 * qJD(5) - t78 * t109;
t93 = t48 * t64;
t158 = t21 - t93;
t25 = -qJD(3) * pkin(3) - t162;
t143 = t50 * t64;
t22 = t50 * qJD(5) - t80 * t109;
t157 = -t22 + t143;
t119 = qJD(3) * qJ(4);
t35 = t81 * t53 + t123;
t27 = -t35 - t119;
t111 = t80 * t127;
t75 = t81 ^ 2;
t133 = qJD(1) * t75;
t141 = t64 * t79;
t156 = -(t133 - t141) * t125 - t64 * t111;
t100 = t5 * t78 - t6 * t80;
t155 = -qJD(5) * t100 + t1 * t78 + t2 * t80;
t152 = pkin(4) + t65;
t131 = qJD(3) * t79;
t136 = -qJD(3) * t71 + t53 * t131;
t11 = (-pkin(4) * t124 + qJD(4)) * qJD(3) - t136;
t151 = t11 * t78;
t150 = t11 * t80;
t149 = t22 * t78;
t31 = t35 * qJD(3);
t148 = t31 * t79;
t147 = t31 * t81;
t146 = t48 * t80;
t145 = t48 * t81;
t144 = t50 * t48;
t142 = t50 * t78;
t140 = t64 * t82;
t139 = t79 * t22;
t138 = t80 * t21;
t83 = qJD(3) ^ 2;
t72 = t83 * t79;
t73 = t83 * t81;
t130 = qJD(3) * t81;
t137 = t50 * t130 - t21 * t79;
t135 = t161 * t64;
t74 = t79 ^ 2;
t134 = t74 - t75;
t42 = -t81 * pkin(3) + t92;
t28 = qJD(1) * t42;
t54 = qJD(1) * t66;
t129 = qJD(5) * t78;
t128 = qJD(5) * t80;
t69 = pkin(4) * t132;
t17 = -t27 + t69;
t126 = t17 * qJD(3);
t117 = t80 * t141;
t116 = t80 * t133;
t115 = t79 * t125;
t112 = t64 * t128;
t44 = t152 * t81;
t108 = t81 * t118;
t104 = t50 * t114;
t103 = t79 * t108;
t101 = t5 * t80 + t6 * t78;
t43 = t152 * t79;
t10 = t80 * t33 + t78 * t43;
t9 = -t78 * t33 + t80 * t43;
t96 = -0.2e1 * qJD(3) * t28;
t95 = 0.2e1 * qJD(3) * t54;
t94 = t64 * t78;
t91 = t82 * t130 + t17 * t79;
t88 = -t81 * t119 - t122;
t29 = t88 * qJD(1) + t63;
t68 = pkin(3) * t131;
t40 = t68 + t88;
t90 = qJD(1) * t40 + t65 * t83 + t29;
t20 = -qJD(3) * qJD(4) + t136;
t86 = -t20 * t81 + t148 + (t25 * t81 + t27 * t79) * qJD(3);
t85 = -t136 * t81 + t148 + (t34 * t81 - t35 * t79) * qJD(3);
t84 = qJD(1) ^ 2;
t70 = pkin(3) * t124;
t62 = t79 * t84 * t81;
t59 = t80 * t108;
t56 = -0.2e1 * t103;
t55 = 0.2e1 * t103;
t52 = t134 * t84;
t51 = -qJ(4) * t132 + t70;
t41 = -0.2e1 * t134 * t118;
t39 = qJD(3) * t44;
t38 = t152 * t131;
t36 = t99 * qJD(1) + t70;
t32 = t48 * t115;
t26 = t68 + t87;
t24 = t35 + t69;
t19 = t28 * t124;
t8 = t78 * t24 + t80 * t36;
t7 = t80 * t24 - t78 * t36;
t4 = -t10 * qJD(5) - t78 * t26 + t80 * t39;
t3 = t9 * qJD(5) + t80 * t26 + t78 * t39;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t41, t73, t56, -t72, 0, -t65 * t73 + t79 * t95, t65 * t72 + t81 * t95, t85, t85 * t65, 0, -t73, t72, t55, t41, t56, t86, t79 * t96 + t81 * t90, -t79 * t90 + t81 * t96, t28 * t40 + t29 * t42 + t65 * t86, t21 * t78 * t81 + (-t111 + t115) * t50, t104 - t32 + (t138 + t149 + (t142 + t146) * qJD(5)) * t81, t137 + t156, t22 * t80 * t81 - t161 * t48, -t139 + (-t116 - t145) * qJD(3) + t135, (t64 + t124) * t130, t44 * t22 - t38 * t48 + t4 * t64 + (-t17 * t121 + t2) * t79 + (-t17 * t129 + t150 + (qJD(1) * t9 + t5) * qJD(3)) * t81, -t44 * t21 - t3 * t64 - t38 * t50 + (t17 * t125 - t1) * t79 + (-t17 * t128 - t151 + (-qJD(1) * t10 - t6) * qJD(3)) * t81, -t10 * t22 + t9 * t21 - t3 * t48 - t4 * t50 - t100 * t131 + (qJD(5) * t101 - t1 * t80 + t2 * t78) * t81, t1 * t10 + t11 * t44 - t17 * t38 + t2 * t9 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t73, 0, -t136 * t79 - t147 + (t34 * t79 + t35 * t81) * qJD(3), 0, 0, 0, 0, 0, 0, 0, t72, t73, -t20 * t79 - t147 + (t25 * t79 - t27 * t81) * qJD(3), 0, 0, 0, 0, 0, 0, t139 + (-t116 + t145) * qJD(3) + t135, t137 - t156, -t104 - t32 + (-t138 + t149 + (-t142 + t146) * qJD(5)) * t81, (qJD(3) * t101 + t11) * t79 + (t126 - t155) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t52, 0, t62, 0, 0, -t54 * t124, -t34 * qJD(3) - t54 * t132 + t136, 0, 0, 0, 0, 0, -t62, t52, t62, 0, -t51 * t132 + t19, (0.2e1 * qJD(4) + t34) * qJD(3) + (t28 * t81 + t51 * t79) * qJD(1) - t136, -t31 * pkin(3) - t20 * qJ(4) + t162 * t27 - t25 * t35 - t28 * t51, -t50 * t94 - t138, (-t22 - t143) * t80 + (t21 + t93) * t78, -t64 * t129 + t59 + (-t78 * t141 - t50 * t81) * qJD(1), t80 * t93 + t149, -t112 + (-t117 + (t48 - t125) * t81) * qJD(1), -t64 * t132, qJ(4) * t22 + t151 - t7 * t64 + t120 * t48 + (-t78 * t140 + t17 * t80) * qJD(5) + (-t5 * t81 + t80 * t91) * qJD(1), -qJ(4) * t21 + t150 + t8 * t64 + t120 * t50 + (-t80 * t140 - t17 * t78) * qJD(5) + (t6 * t81 - t78 * t91) * qJD(1), t8 * t48 + t7 * t50 + (-t6 * t124 + t21 * t82 - t2 + (-t48 * t82 - t6) * qJD(5)) * t80 + (t5 * t124 - t22 * t82 - t1 + (t50 * t82 + t5) * qJD(5)) * t78, t11 * qJ(4) + t120 * t17 + t155 * t82 - t5 * t7 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t74 * t84 - t83, t19 + (t27 + t35) * qJD(3), 0, 0, 0, 0, 0, 0, -qJD(3) * t48 - t64 * t94 + t59, -t112 - qJD(3) * t50 + (-t81 * t125 - t117) * qJD(1), t157 * t78 + t158 * t80, t159 * t80 + t160 * t78 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t48 ^ 2 + t50 ^ 2, -t158, -t144, t157, t108, -t17 * t50 + t159, t17 * t48 - t160, 0, 0;];
tauc_reg = t12;
