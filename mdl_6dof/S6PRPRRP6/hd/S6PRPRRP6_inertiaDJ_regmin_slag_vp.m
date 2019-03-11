% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:55
% EndTime: 2019-03-08 20:20:59
% DurationCPUTime: 1.19s
% Computational Cost: add. (865->173), mult. (2154->294), div. (0->0), fcn. (1836->8), ass. (0->110)
t48 = sin(pkin(6));
t55 = cos(qJ(2));
t117 = t48 * t55;
t49 = cos(pkin(6));
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t24 = -t51 * t117 + t49 * t54;
t52 = sin(qJ(2));
t118 = t48 * t52;
t40 = qJD(2) * t118;
t13 = t24 * qJD(4) - t54 * t40;
t23 = t54 * t117 + t49 * t51;
t53 = cos(qJ(5));
t43 = qJD(5) * t53;
t50 = sin(qJ(5));
t6 = t13 * t50 + t23 * t43;
t107 = qJD(4) * t23;
t121 = t54 * pkin(9);
t74 = t51 * pkin(4) - t121;
t34 = qJ(3) + t74;
t32 = t50 * t34;
t56 = -pkin(2) - pkin(8);
t39 = t53 * t51 * t56;
t64 = t39 + t32;
t44 = t50 ^ 2;
t46 = t53 ^ 2;
t113 = t44 - t46;
t80 = t113 * qJD(5);
t71 = pkin(5) * t50 - qJ(6) * t53;
t20 = t71 * qJD(5) - t50 * qJD(6);
t115 = t54 * t20;
t72 = t53 * pkin(5) + t50 * qJ(6);
t35 = -pkin(4) - t72;
t119 = t35 * t51;
t63 = -t56 + t71;
t19 = t63 * t54;
t126 = (t119 + t121) * qJD(4) - qJD(5) * t19 - t115;
t122 = pkin(9) * t51;
t75 = pkin(4) * t54 + t122;
t33 = t75 * qJD(4) + qJD(3);
t61 = -t64 * qJD(5) + t53 * t33;
t125 = t72 * qJD(5) - t53 * qJD(6);
t124 = 0.2e1 * qJD(3);
t123 = 0.2e1 * qJD(6);
t116 = t50 * t56;
t114 = t54 * t56;
t112 = t44 + t46;
t45 = t51 ^ 2;
t47 = t54 ^ 2;
t111 = t45 - t47;
t110 = t45 + t47;
t109 = qJD(2) * t55;
t108 = qJD(4) * t19;
t106 = qJD(4) * t53;
t105 = qJD(5) * t50;
t104 = qJD(5) * t54;
t103 = qJD(5) * t56;
t102 = t51 * qJD(4);
t100 = t54 * qJD(4);
t99 = qJ(3) * qJD(4);
t98 = qJ(6) * qJD(4);
t97 = -0.2e1 * pkin(4) * qJD(5);
t84 = t56 * t100;
t96 = t50 * t33 + t34 * t43 + t53 * t84;
t95 = pkin(9) * t105;
t94 = pkin(9) * t43;
t93 = t50 * t104;
t92 = t50 * t103;
t91 = t53 * t104;
t89 = t47 * t103;
t88 = t48 * t109;
t87 = t50 * t43;
t86 = t53 * t102;
t85 = t51 * t100;
t83 = t54 * t98;
t82 = -pkin(5) + t116;
t81 = t112 * t54;
t79 = t111 * qJD(4);
t78 = 0.2e1 * t85;
t77 = t50 * t84;
t76 = t50 * t86;
t15 = t50 * t118 + t24 * t53;
t65 = t53 * t118 - t24 * t50;
t70 = -t15 * t50 - t53 * t65;
t69 = t15 * t53 - t50 * t65;
t16 = t51 * qJ(6) + t64;
t17 = -t53 * t34 + t82 * t51;
t68 = t16 * t53 + t17 * t50;
t67 = t16 * t50 - t17 * t53;
t7 = t23 * t105 - t13 * t53;
t10 = -t63 * t102 + t125 * t54;
t62 = -t10 + (t35 * t54 - t122) * qJD(5);
t12 = -t51 * t40 + t107;
t2 = t15 * qJD(5) - t12 * t50 - t53 * t88;
t3 = t65 * qJD(5) - t12 * t53 + t50 * t88;
t59 = t70 * qJD(5) + t2 * t50 + t3 * t53;
t4 = t83 + (qJD(6) - t92) * t51 + t96;
t5 = t82 * t100 - t61;
t58 = -t67 * qJD(5) + t4 * t53 + t5 * t50;
t57 = (-t50 * t107 - t2) * t51 + (qJD(4) * t65 + t6) * t54;
t30 = -t50 * t102 + t91;
t29 = t50 * t100 + t51 * t43;
t28 = t110 * t43;
t27 = -t86 - t93;
t26 = -t53 * t100 + t51 * t105;
t25 = t110 * t105;
t9 = t61 - t77;
t8 = t51 * t92 - t96;
t1 = (t23 * t106 + t3) * t51 + (qJD(4) * t15 + t7) * t54;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23 * t13 + 0.2e1 * t15 * t3 - 0.2e1 * t2 * t65; 0, 0, -t40, -t88, t40, t88 (qJD(3) * t52 + (-pkin(2) * t52 + qJ(3) * t55) * qJD(2)) * t48, 0, 0, 0, 0, 0 (t52 * t100 + t51 * t109) * t48 (-t52 * t102 + t54 * t109) * t48, 0, 0, 0, 0, 0, t57, -t1, t57, -t70 * t102 + (-qJD(5) * t69 + t2 * t53 - t3 * t50) * t54, t1, t23 * t10 + t13 * t19 + t15 * t4 + t3 * t16 + t2 * t17 - t5 * t65; 0, 0, 0, 0, 0, t124, qJ(3) * t124, -0.2e1 * t85, 0.2e1 * t79, 0, 0, 0, 0.2e1 * qJD(3) * t51 + 0.2e1 * t54 * t99, 0.2e1 * qJD(3) * t54 - 0.2e1 * t51 * t99, -0.2e1 * t46 * t85 - 0.2e1 * t47 * t87, 0.2e1 * t47 * t80 + 0.4e1 * t54 * t76, -0.2e1 * t111 * t106 - 0.2e1 * t51 * t93, 0.2e1 * t50 * t79 - 0.2e1 * t51 * t91, t78, 0.2e1 * (t100 * t34 - t89) * t53 + 0.2e1 * (t9 + t77) * t51, 0.2e1 * t50 * t89 + 0.2e1 * t8 * t51 + 0.2e1 * (-t32 + t39) * t100, 0.2e1 * (-t50 * t108 - t5) * t51 + 0.2e1 * (-qJD(4) * t17 + t10 * t50 + t19 * t43) * t54, 0.2e1 * t67 * t102 + 0.2e1 * (-qJD(5) * t68 - t4 * t50 + t5 * t53) * t54, 0.2e1 * (t19 * t106 + t4) * t51 + 0.2e1 * (qJD(4) * t16 - t10 * t53 + t19 * t105) * t54, 0.2e1 * t19 * t10 + 0.2e1 * t16 * t4 + 0.2e1 * t17 * t5; 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t69 * qJD(4) - t13) * t54 + (t59 + t107) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t25, -t28, 0, -t25 (t68 * qJD(4) - t10) * t54 + (t58 + t108) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t112) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, 0, 0, 0, 0, 0, t7, t6, t7, t59, -t6, pkin(9) * t59 + t13 * t35 + t23 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t100, 0, -t56 * t102, -t84, -t54 * t80 - t76, t113 * t102 - 0.4e1 * t54 * t87, t29, -t26, 0 (-t50 * t114 - t53 * t75) * qJD(5) + (t50 * t74 - t39) * qJD(4) (-t53 * t114 + t50 * t75) * qJD(5) + (-t53 * t121 + (pkin(4) * t53 + t116) * t51) * qJD(4), -t126 * t50 + t62 * t53, t58, t126 * t53 + t62 * t50, pkin(9) * t58 + t10 * t35 + t19 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t100, 0, 0, 0, 0, 0, t27, -t30, t27, qJD(4) * t81, t30, -t115 + (pkin(9) * t81 + t119) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t80, 0, 0, 0, t50 * t97, t53 * t97, 0.2e1 * t35 * t105 - 0.2e1 * t20 * t53, 0, -0.2e1 * t20 * t50 - 0.2e1 * t35 * t43, 0.2e1 * t35 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(5) + t3 * qJ(6) + t15 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t30, t100, t9, t8 (0.2e1 * pkin(5) - t116) * t100 + t61 (pkin(5) * t102 - qJ(6) * t104) * t53 + (t51 * t98 + (pkin(5) * qJD(5) - qJD(6)) * t54) * t50, 0.2e1 * t83 + (t123 - t92) * t51 + t96, -t5 * pkin(5) + t4 * qJ(6) + t16 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t26, -t29, 0, -t26, -t100 * t71 - t125 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t105, 0, -t94, t95, -t94, -t125, -t95, -t125 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, qJ(6) * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t27, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
