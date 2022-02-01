% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:19
% EndTime: 2022-01-23 09:30:22
% DurationCPUTime: 0.77s
% Computational Cost: add. (1296->152), mult. (3107->202), div. (0->0), fcn. (2024->6), ass. (0->107)
t132 = cos(qJ(4));
t77 = sin(qJ(3));
t78 = cos(qJ(3));
t66 = sin(pkin(8)) * pkin(1) + pkin(6);
t133 = pkin(7) + t66;
t98 = t133 * qJD(1);
t37 = t77 * qJD(2) + t98 * t78;
t76 = sin(qJ(4));
t31 = t76 * t37;
t120 = qJD(3) * pkin(3);
t36 = t78 * qJD(2) - t98 * t77;
t34 = t36 + t120;
t103 = t132 * t34 - t31;
t55 = t132 * t77 + t76 * t78;
t116 = qJD(1) * t55;
t117 = t116 * qJ(5);
t6 = t103 - t117;
t71 = qJD(3) + qJD(4);
t136 = t116 ^ 2;
t5 = t71 * pkin(4) + t6;
t135 = t5 - t6;
t134 = pkin(3) * t71;
t100 = t132 * qJD(4);
t106 = t132 * t78;
t127 = t76 * t77;
t96 = t71 * t127;
t26 = -qJD(3) * t106 - t78 * t100 + t96;
t21 = t26 * t71;
t115 = qJD(1) * t77;
t107 = t76 * t115;
t97 = qJD(1) * t106;
t47 = -t97 + t107;
t131 = t47 * t71;
t130 = t116 * t47;
t67 = -cos(pkin(8)) * pkin(1) - pkin(2);
t56 = -t78 * pkin(3) + t67;
t51 = t56 * qJD(1);
t128 = t51 * t116;
t79 = qJD(3) ^ 2;
t126 = t79 * t77;
t125 = t79 * t78;
t27 = t71 * t55;
t20 = t27 * qJD(1);
t124 = -t55 * t20 + t26 * t47;
t123 = t132 * t36 - t31;
t122 = t71 * t97;
t121 = t77 ^ 2 - t78 ^ 2;
t19 = qJD(1) * t96 - t122;
t119 = t19 * qJ(5);
t118 = t47 * qJ(5);
t58 = qJD(1) * t67;
t114 = qJD(4) * t76;
t113 = t58 * qJD(1);
t99 = t47 * pkin(4) + qJD(5);
t23 = t51 + t99;
t111 = qJD(5) + t23;
t110 = qJD(1) * qJD(3);
t109 = t77 * t120;
t108 = pkin(3) * t115;
t33 = t132 * t37;
t105 = t77 * t110;
t65 = pkin(3) * t105;
t17 = t20 * pkin(4) + t65;
t29 = t36 * qJD(3);
t30 = t37 * qJD(3);
t104 = -t132 * t30 - t76 * t29;
t102 = -t76 * t36 - t33;
t101 = qJD(3) * t133;
t54 = -t106 + t127;
t95 = t116 * t27 - t54 * t19;
t93 = 0.2e1 * qJD(3) * t58;
t92 = -t76 * t34 - t33;
t52 = t133 * t77;
t53 = t133 * t78;
t91 = -t132 * t53 + t76 * t52;
t44 = t77 * t101;
t45 = t78 * t101;
t90 = -t52 * t100 - t53 * t114 - t132 * t44 - t76 * t45;
t89 = -t71 * t107 + t122;
t88 = t92 * qJD(4) + t104;
t87 = t91 * qJD(4) - t132 * t45 + t76 * t44;
t86 = t34 * t100 - t37 * t114 + t132 * t29 - t76 * t30;
t85 = t88 + t119;
t84 = t51 * t47 - t86;
t83 = -t20 * qJ(5) + t86;
t82 = (-t33 + (-t34 - t134) * t76) * qJD(4) + t104;
t81 = t111 * t47 - t83;
t80 = qJD(1) ^ 2;
t69 = t132 * pkin(3) + pkin(4);
t59 = t100 * t134;
t46 = t47 ^ 2;
t38 = pkin(4) * t116 + t108;
t35 = t54 * pkin(4) + t56;
t22 = t27 * t71;
t18 = t27 * pkin(4) + t109;
t15 = -t46 + t136;
t13 = -t54 * qJ(5) - t91;
t12 = -t55 * qJ(5) - t132 * t52 - t76 * t53;
t10 = t89 + t131;
t9 = -t117 + t123;
t8 = t102 + t118;
t7 = -t92 - t118;
t4 = t26 * qJ(5) - t55 * qJD(5) + t87;
t3 = -t27 * qJ(5) - t54 * qJD(5) + t90;
t2 = -qJD(5) * t116 + t85;
t1 = -t47 * qJD(5) + t83;
t11 = [0, 0, 0, 0, 0.2e1 * t78 * t105, -0.2e1 * t121 * t110, t125, -t126, 0, -t66 * t125 + t77 * t93, t66 * t126 + t78 * t93, -t116 * t26 - t19 * t55, -t95 + t124, -t21, -t22, 0, t47 * t109 + t56 * t20 + t51 * t27 + t54 * t65 + t87 * t71, 0.2e1 * t116 * t109 - t56 * t19 - t51 * t26 - t90 * t71, t17 * t54 + t18 * t47 + t35 * t20 + t23 * t27 + t4 * t71, t116 * t18 + t17 * t55 - t35 * t19 - t23 * t26 - t3 * t71, -t1 * t54 - t116 * t4 + t12 * t19 - t13 * t20 - t2 * t55 + t5 * t26 - t7 * t27 - t3 * t47, t1 * t13 + t2 * t12 + t17 * t35 + t23 * t18 + t7 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t125, 0, 0, 0, 0, 0, -t22, t21, -t22, t21, t95 + t124, t1 * t55 - t2 * t54 - t7 * t26 - t5 * t27; 0, 0, 0, 0, -t77 * t80 * t78, t121 * t80, 0, 0, 0, -t77 * t113, -t78 * t113, t130, t15, t10, 0, 0, -t102 * t71 - t47 * t108 - t128 + t82, -t108 * t116 + t123 * t71 - t59 + t84, -t111 * t116 - t38 * t47 - t8 * t71 + t119 + t82, -t116 * t38 + t9 * t71 - t59 + t81, t69 * t19 + (t7 + t8) * t116 + (-t5 + t9) * t47 + (-t20 * t76 + (t116 * t76 - t132 * t47) * qJD(4)) * pkin(3), t2 * t69 - t23 * t38 - t5 * t8 - t7 * t9 + (t1 * t76 + (t132 * t7 - t5 * t76) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t15, t10, 0, 0, -t92 * t71 - t128 + t88, t103 * t71 + t84, t7 * t71 + (-t23 - t99) * t116 + t85, -t136 * pkin(4) + t6 * t71 + t81, t19 * pkin(4) - t135 * t47, t135 * t7 + (-t116 * t23 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t116 + t20, t89 - t131, -t46 - t136, t116 * t5 + t7 * t47 + t17;];
tauc_reg = t11;
