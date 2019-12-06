% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:36
% EndTime: 2019-12-05 17:00:42
% DurationCPUTime: 1.49s
% Computational Cost: add. (907->171), mult. (2677->321), div. (0->0), fcn. (2300->8), ass. (0->118)
t69 = sin(qJ(2));
t128 = qJD(2) * t69;
t65 = sin(pkin(5));
t107 = t65 * t128;
t72 = cos(qJ(2));
t134 = t65 * t72;
t67 = sin(qJ(4));
t115 = t67 * t134;
t127 = qJD(2) * t72;
t106 = t65 * t127;
t135 = t65 * t69;
t66 = cos(pkin(5));
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t35 = t68 * t135 - t66 * t71;
t16 = -qJD(3) * t35 + t71 * t106;
t36 = t71 * t135 + t66 * t68;
t70 = cos(qJ(4));
t60 = qJD(4) * t70;
t5 = -qJD(4) * t115 - t70 * t107 + t16 * t67 + t36 * t60;
t18 = t70 * t134 + t36 * t67;
t6 = -t18 * qJD(4) + t67 * t107 + t16 * t70;
t19 = t36 * t70 - t115;
t90 = t18 * t70 - t19 * t67;
t79 = t90 * qJD(4) + t5 * t67 + t6 * t70;
t137 = t68 * pkin(8);
t95 = -t71 * pkin(3) - t137;
t49 = -pkin(2) + t95;
t139 = pkin(7) * t71;
t57 = t70 * t139;
t144 = t67 * t49 + t57;
t61 = t67 ^ 2;
t63 = t70 ^ 2;
t130 = t61 - t63;
t48 = t130 * qJD(4);
t123 = qJD(4) * t71;
t110 = t67 * t123;
t59 = t68 * qJD(3);
t38 = t70 * t59 + t110;
t138 = pkin(8) * t71;
t94 = pkin(3) * t68 - t138;
t83 = t94 * qJD(3);
t12 = t38 * pkin(7) - t49 * t60 - t67 * t83;
t91 = pkin(4) * t67 - qJ(5) * t70;
t85 = pkin(7) + t91;
t33 = t85 * t68;
t34 = t91 * qJD(4) - t67 * qJD(5);
t92 = t70 * pkin(4) + t67 * qJ(5);
t45 = -pkin(3) - t92;
t143 = qJD(3) * (-t45 * t71 + t137) - qJD(4) * t33 - t34 * t68;
t142 = t92 * qJD(4) - t70 * qJD(5);
t103 = t67 * t59;
t13 = pkin(7) * t103 - qJD(4) * t144 + t70 * t83;
t141 = 0.2e1 * qJD(5);
t140 = pkin(7) * t67;
t17 = t36 * qJD(3) + t68 * t106;
t136 = t17 * t68;
t14 = t35 * t17;
t133 = t70 * t49;
t62 = t68 ^ 2;
t129 = -t71 ^ 2 + t62;
t126 = qJD(3) * t70;
t125 = qJD(4) * t67;
t124 = qJD(4) * t68;
t121 = t71 * qJD(3);
t120 = t71 * qJD(5);
t119 = qJ(5) * qJD(3);
t118 = t67 * t139;
t117 = -0.2e1 * pkin(2) * qJD(3);
t116 = -0.2e1 * pkin(3) * qJD(4);
t114 = pkin(8) * t125;
t113 = pkin(8) * t60;
t112 = pkin(4) * t59;
t111 = pkin(7) * t121;
t109 = t70 * t123;
t108 = t35 * t60;
t105 = t67 * t121;
t104 = t67 * t60;
t102 = t68 * t121;
t101 = t70 * t121;
t100 = t68 * t119;
t99 = t129 * qJD(3);
t98 = 0.2e1 * t102;
t97 = t67 * t101;
t96 = t62 * t104;
t24 = -t71 * qJ(5) + t144;
t26 = -t133 + (pkin(4) + t140) * t71;
t89 = -t24 * t67 + t26 * t70;
t31 = -t118 + t133;
t88 = -t144 * t67 - t31 * t70;
t9 = t35 * t125 - t17 * t70;
t84 = 0.2e1 * t18 * t5 + 0.2e1 * t19 * t6 + 0.2e1 * t14;
t82 = t35 * t105 + t68 * t108 + t67 * t136 - t18 * t59 + t5 * t71;
t10 = -t112 - t13;
t7 = t100 - t12 - t120;
t77 = t89 * qJD(4) + t10 * t67 + t7 * t70;
t76 = t88 * qJD(4) - t12 * t70 - t13 * t67;
t75 = t16 * t71 + t136 + (t35 * t71 - t36 * t68) * qJD(3);
t74 = t79 * pkin(8);
t73 = t90 * t121 + (t5 * t70 - t6 * t67 + (-t18 * t67 - t19 * t70) * qJD(4)) * t68;
t55 = -0.2e1 * t102;
t54 = -0.2e1 * t104;
t53 = 0.2e1 * t104;
t52 = pkin(8) * t109;
t40 = t103 - t109;
t39 = t68 * t60 + t105;
t37 = -t67 * t124 + t101;
t30 = 0.2e1 * t63 * t102 - 0.2e1 * t96;
t29 = 0.2e1 * t61 * t102 + 0.2e1 * t96;
t25 = t130 * t124 - t97;
t23 = t68 * t109 - t67 * t99;
t22 = 0.4e1 * t68 * t104 + t130 * t121;
t21 = 0.2e1 * t68 * t110 + 0.2e1 * t129 * t126;
t20 = t62 * t48 - 0.2e1 * t68 * t97;
t11 = t85 * t121 + t142 * t68;
t8 = t17 * t67 + t108;
t1 = (t35 * t126 + t6) * t71 + (-qJD(3) * t19 - t9) * t68;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t65 ^ 2 * t69 * t127 + 0.2e1 * t36 * t16 + 0.2e1 * t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t106, 0, 0, 0, 0, 0, 0, 0, 0, (-t71 * t128 - t72 * t59) * t65, (-t72 * t121 + t68 * t128) * t65, t75, -pkin(2) * t107 + pkin(7) * t75, 0, 0, 0, 0, 0, 0, t82, t1, t73, -t19 * t12 - t18 * t13 - t5 * t31 + t6 * t144 + (t35 * t121 + t136) * pkin(7), 0, 0, 0, 0, 0, 0, t82, t73, -t1, t18 * t10 + t35 * t11 + t17 * t33 + t19 * t7 + t6 * t24 + t5 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -0.2e1 * t99, 0, t55, 0, 0, t68 * t117, t71 * t117, 0, 0, t30, 0.2e1 * t20, t21, t29, 0.2e1 * t23, t55, 0.2e1 * t31 * t59 - 0.2e1 * t13 * t71 + 0.2e1 * (t62 * t60 + t67 * t98) * pkin(7), -0.2e1 * t144 * t59 - 0.2e1 * t12 * t71 + 0.2e1 * (-t62 * t125 + t70 * t98) * pkin(7), 0.2e1 * t88 * t121 + 0.2e1 * (t12 * t67 - t13 * t70 + (-t144 * t70 + t31 * t67) * qJD(4)) * t68, 0.2e1 * pkin(7) ^ 2 * t102 - 0.2e1 * t12 * t144 + 0.2e1 * t31 * t13, t30, t21, -0.2e1 * t20, t55, -0.2e1 * t23, t29, 0.2e1 * (qJD(3) * t33 * t67 + t10) * t71 + 0.2e1 * (-qJD(3) * t26 + t11 * t67 + t33 * t60) * t68, 0.2e1 * t89 * t121 + 0.2e1 * (t10 * t70 - t67 * t7 + (-t24 * t70 - t26 * t67) * qJD(4)) * t68, 0.2e1 * (-t33 * t126 - t7) * t71 + 0.2e1 * (qJD(3) * t24 - t11 * t70 + t33 * t125) * t68, 0.2e1 * t26 * t10 + 0.2e1 * t33 * t11 + 0.2e1 * t24 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t16, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, t79, -t17 * pkin(3) + t74, 0, 0, 0, 0, 0, 0, t9, t79, -t8, t17 * t45 + t35 * t34 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, -t59, 0, -t111, pkin(7) * t59, 0, 0, -t25, -t22, t40, t25, t38, 0, t52 + (-pkin(3) * t70 + t140) * t124 + (t67 * t95 - t57) * qJD(3), (pkin(7) * t68 * t70 + t67 * t94) * qJD(4) + (t70 * t95 + t118) * qJD(3), t76, -pkin(3) * t111 + pkin(8) * t76, -t25, t40, t22, 0, -t38, t25, t52 + (t45 * t124 - t11) * t70 - t143 * t67, t77, (-t11 + (t45 * t68 + t138) * qJD(4)) * t67 + t143 * t70, pkin(8) * t77 + t11 * t45 + t33 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -0.2e1 * t48, 0, t54, 0, 0, t67 * t116, t70 * t116, 0, 0, t53, 0, 0.2e1 * t48, 0, 0, t54, 0.2e1 * t45 * t125 - 0.2e1 * t34 * t70, 0, -0.2e1 * t34 * t67 - 0.2e1 * t45 * t60, 0.2e1 * t45 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, t6, -t5 * pkin(4) + t6 * qJ(5) + t19 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t39, t59, t13, t12, 0, 0, 0, t37, 0, t59, t39, 0, t13 + 0.2e1 * t112, (-pkin(4) * t121 - qJ(5) * t124) * t70 + (-t71 * t119 + (pkin(4) * qJD(4) - qJD(5)) * t68) * t67, 0.2e1 * t100 - t12 - 0.2e1 * t120, -t10 * pkin(4) + t7 * qJ(5) + t24 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t125, 0, -t113, t114, 0, 0, 0, t60, 0, 0, t125, 0, -t113, -t142, -t114, -t142 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, qJ(5) * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t37, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
