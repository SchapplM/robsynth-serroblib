% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:52
% DurationCPUTime: 0.67s
% Computational Cost: add. (1089->165), mult. (2389->206), div. (0->0), fcn. (1471->4), ass. (0->106)
t80 = -pkin(1) - pkin(6);
t59 = t80 * qJD(1) + qJD(2);
t138 = -qJ(4) * qJD(1) + t59;
t77 = cos(pkin(7));
t102 = qJD(1) * qJD(3);
t79 = cos(qJ(3));
t97 = t79 * t102;
t58 = t77 * t97;
t78 = sin(qJ(3));
t115 = qJD(1) * t78;
t76 = sin(pkin(7));
t98 = t76 * t115;
t35 = qJD(3) * t98 - t58;
t52 = t76 * t79 + t77 * t78;
t45 = t52 * qJD(1);
t112 = qJD(3) * t79;
t113 = qJD(3) * t78;
t47 = -t77 * t112 + t76 * t113;
t92 = -t52 * t35 - t45 * t47;
t50 = t52 * qJD(3);
t39 = t50 * qJD(3);
t137 = -qJD(1) * t45 - t39;
t114 = qJD(1) * t79;
t48 = t77 * t114 - t98;
t132 = t48 ^ 2;
t42 = t45 ^ 2;
t136 = -t42 - t132;
t135 = -t42 + t132;
t134 = qJD(3) * (t48 + t98) - t58;
t108 = t78 * qJD(4);
t117 = qJ(4) - t80;
t96 = t117 * t79;
t37 = -qJD(3) * t96 - t108;
t107 = t79 * qJD(4);
t87 = t117 * t113 - t107;
t14 = t76 * t37 - t77 * t87;
t15 = t77 * t37 + t76 * t87;
t57 = t117 * t78;
t27 = -t76 * t57 + t77 * t96;
t28 = -t77 * t57 - t76 * t96;
t36 = qJD(1) * t50;
t133 = t14 * t48 - t15 * t45 - t27 * t36 + t28 * t35;
t72 = qJD(1) * qJD(2);
t131 = 0.2e1 * t72;
t103 = qJ(4) * qJD(3);
t29 = t59 * t112 + (-t79 * t103 - t108) * qJD(1);
t83 = -t59 * t113 + (t78 * t103 - t107) * qJD(1);
t4 = t76 * t29 - t77 * t83;
t130 = t4 * t27;
t53 = -t76 * t78 + t77 * t79;
t129 = t4 * t53;
t56 = pkin(3) * t115 + qJD(1) * qJ(2) + qJD(4);
t13 = t45 * pkin(4) - t48 * qJ(5) + t56;
t128 = t13 * t48;
t124 = t48 * t45;
t40 = t138 * t78;
t123 = t76 * t40;
t32 = t77 * t40;
t81 = qJD(3) ^ 2;
t122 = t81 * t78;
t121 = t81 * t79;
t5 = t77 * t29 + t76 * t83;
t41 = t138 * t79;
t34 = qJD(3) * pkin(3) + t41;
t12 = t76 * t34 + t32;
t55 = pkin(3) * t97 + t72;
t120 = t78 ^ 2 - t79 ^ 2;
t82 = qJD(1) ^ 2;
t119 = -t81 - t82;
t118 = t82 * qJ(2);
t67 = t78 * pkin(3) + qJ(2);
t111 = t14 * qJD(3);
t110 = t15 * qJD(3);
t38 = t47 * qJD(3);
t109 = t56 * qJD(1);
t60 = pkin(3) * t112 + qJD(2);
t17 = t77 * t41 - t123;
t106 = qJD(5) - t17;
t105 = qJ(2) * qJD(3);
t101 = t79 * t82 * t78;
t100 = 0.2e1 * qJD(1);
t99 = pkin(3) * t114;
t95 = qJD(1) * t48 - t38;
t94 = t78 * t97;
t11 = t77 * t34 - t123;
t2 = -t53 * t36 - t50 * t48;
t16 = t76 * t41 + t32;
t91 = t16 * qJD(3) - t4;
t90 = -t35 * pkin(4) + t36 * qJ(5) + t55;
t89 = -t2 - t92;
t88 = t58 + (t48 - t98) * qJD(3);
t3 = qJD(3) * qJD(5) + t5;
t8 = -qJD(3) * pkin(4) + qJD(5) - t11;
t9 = qJD(3) * qJ(5) + t12;
t86 = t3 * t52 - t9 * t47 + t8 * t50 - t129;
t85 = -t11 * t50 - t12 * t47 + t5 * t52 - t129;
t84 = t53 * t35 + t36 * t52 + t50 * t45 + t48 * t47;
t70 = qJ(2) * t131;
t66 = -t77 * pkin(3) - pkin(4);
t64 = t76 * pkin(3) + qJ(5);
t21 = t52 * pkin(4) - t53 * qJ(5) + t67;
t19 = 0.2e1 * t45 * qJD(3);
t18 = t48 * pkin(4) + t45 * qJ(5) + t99;
t7 = -t47 * pkin(4) + t50 * qJ(5) - t53 * qJD(5) + t60;
t1 = -t48 * qJD(5) + t90;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t70, -0.2e1 * t94, 0.2e1 * t120 * t102, -t122, 0.2e1 * t94, -t121, 0, -t80 * t122 + (qJD(2) * t78 + t79 * t105) * t100, -t80 * t121 + (qJD(2) * t79 - t78 * t105) * t100, 0, t70, t2, t84, -t39, t92, t38, 0, -t67 * t35 + t60 * t45 - t56 * t47 + t55 * t52 - t111, -t67 * t36 + t60 * t48 - t56 * t50 + t55 * t53 - t110, -t85 + t133, -t11 * t14 + t12 * t15 + t5 * t28 + t55 * t67 + t56 * t60 + t130, t2, -t39, -t84, 0, -t38, t92, t1 * t52 - t13 * t47 - t21 * t35 + t7 * t45 - t111, -t86 + t133, -t1 * t53 + t13 * t50 + t21 * t36 - t7 * t48 + t110, t1 * t21 + t13 * t7 + t8 * t14 + t9 * t15 + t3 * t28 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t118, 0, 0, 0, 0, 0, 0, t119 * t78, t119 * t79, 0, -t118, 0, 0, 0, 0, 0, 0, t137, -t95, t89, t85 - t109, 0, 0, 0, 0, 0, 0, t137, t89, t95, -t13 * qJD(1) + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t120 * t82, 0, -t101, 0, 0, -t79 * t118, t78 * t118, 0, 0, t124, t135, 0, -t124, t134, 0, -t45 * t99 - t56 * t48 + t91, t17 * qJD(3) + t56 * t45 - t48 * t99 - t5, (t12 - t16) * t48 + (-t11 + t17) * t45 + (t35 * t76 + t36 * t77) * pkin(3), t11 * t16 - t12 * t17 + (-t79 * t109 - t4 * t77 + t5 * t76) * pkin(3), t124, 0, -t135, 0, -t134, -t124, -t18 * t45 - t128 + t91, t64 * t35 - t66 * t36 + (-t16 + t9) * t48 + (t8 - t106) * t45, -t13 * t45 + t18 * t48 + (0.2e1 * qJD(5) - t17) * qJD(3) + t5, t106 * t9 - t13 * t18 - t8 * t16 + t3 * t64 + t4 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t19, t136, t11 * t48 + t12 * t45 + t55, 0, 0, 0, 0, 0, 0, t88, t136, t19, t9 * t45 + (-qJD(5) - t8) * t48 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, -t132 - t81, -t9 * qJD(3) + t128 + t4;];
tauc_reg = t6;
