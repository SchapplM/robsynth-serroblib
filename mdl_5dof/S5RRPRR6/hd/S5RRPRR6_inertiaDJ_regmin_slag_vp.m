% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:54
% EndTime: 2020-01-03 12:05:56
% DurationCPUTime: 0.64s
% Computational Cost: add. (770->113), mult. (1955->182), div. (0->0), fcn. (1591->8), ass. (0->109)
t137 = qJD(4) + qJD(5);
t83 = sin(qJ(5));
t84 = sin(qJ(4));
t86 = cos(qJ(5));
t87 = cos(qJ(4));
t90 = t83 * t87 + t84 * t86;
t38 = t137 * t90;
t81 = sin(pkin(9));
t23 = t38 * t81;
t138 = 0.2e1 * t23;
t136 = 2 * qJD(4);
t135 = pkin(8) * t81;
t88 = cos(qJ(2));
t134 = t88 * pkin(1);
t82 = cos(pkin(9));
t133 = t38 * t82;
t118 = pkin(1) * qJD(2);
t106 = t88 * t118;
t69 = qJD(3) + t106;
t79 = t81 ^ 2;
t56 = t79 * t69;
t132 = t81 * t84;
t131 = t82 * t84;
t130 = t82 * t87;
t129 = t86 * t87;
t115 = qJD(4) * t84;
t100 = t81 * t115;
t113 = qJD(5) * t83;
t24 = t137 * t81 * t129 - t83 * t100 - t113 * t132;
t114 = qJD(4) * t87;
t99 = t81 * t114;
t66 = pkin(4) * t99;
t40 = t81 * t69 + t66;
t43 = t90 * t81;
t75 = pkin(4) * t132;
t85 = sin(qJ(2));
t76 = t85 * pkin(1) + qJ(3);
t45 = t81 * t76 + t75;
t128 = t45 * t24 + t40 * t43;
t48 = t81 * qJD(3) + t66;
t53 = t81 * qJ(3) + t75;
t127 = t53 * t24 + t48 * t43;
t116 = qJD(3) * t82;
t59 = -t82 * pkin(3) - t81 * pkin(7) - pkin(2);
t126 = t59 * t114 + t87 * t116;
t101 = t79 * t114;
t125 = t76 * t101 + t84 * t56;
t107 = t85 * t118;
t47 = t59 - t134;
t108 = t84 * t107 + t47 * t114 + t69 * t130;
t73 = t82 * t115;
t15 = t76 * t73 - t108;
t124 = -t15 * t82 + t87 * t56;
t80 = t82 ^ 2;
t123 = t80 * t69 + t56;
t77 = t79 * qJD(3);
t122 = qJ(3) * t101 + t84 * t77;
t98 = qJ(3) * t115;
t35 = t82 * t98 - t126;
t121 = -t35 * t82 + t87 * t77;
t120 = t80 * qJD(3) + t77;
t119 = t79 + t80;
t117 = qJ(3) * t84;
t112 = t87 * t135;
t109 = t76 * t130;
t16 = -t69 * t131 + t87 * t107 + (-t84 * t47 - t109) * qJD(4);
t65 = pkin(8) * t100;
t10 = t65 + t16;
t97 = t47 - t135;
t22 = t97 * t87 + (-t76 * t84 - pkin(4)) * t82;
t29 = t97 * t84 + t109;
t9 = (-t76 * t131 - t112) * qJD(4) + t108;
t2 = -t83 * t10 - t86 * t9 + (-t22 * t86 + t29 * t83) * qJD(5);
t89 = t83 * t84 - t129;
t44 = t89 * t81;
t111 = -t2 * t82 - t45 * t23 - t40 * t44;
t25 = (-t82 * t117 - t112) * qJD(4) + t126;
t103 = qJ(3) * t130;
t36 = -t84 * t116 + (-t84 * t59 - t103) * qJD(4);
t26 = t65 + t36;
t96 = t59 - t135;
t31 = t96 * t87 + (-pkin(4) - t117) * t82;
t39 = t96 * t84 + t103;
t5 = -t86 * t25 - t83 * t26 + (-t31 * t86 + t39 * t83) * qJD(5);
t110 = -t53 * t23 - t48 * t44 - t5 * t82;
t105 = pkin(4) * t113;
t104 = qJD(5) * t86 * pkin(4);
t102 = t79 * t115;
t95 = t119 * t69;
t94 = t119 * qJD(3);
t93 = t81 * t82 * t136;
t92 = t81 * t107;
t91 = t82 * t107;
t3 = -t83 * t9 + t86 * t10 + (-t22 * t83 - t29 * t86) * qJD(5);
t6 = -t83 * t25 + t86 * t26 + (-t31 * t83 - t39 * t86) * qJD(5);
t74 = t82 * t114;
t64 = t82 * t104;
t63 = t82 * t105;
t58 = -0.2e1 * t84 * t101;
t55 = t87 * t93;
t54 = t84 * t93;
t42 = (t84 ^ 2 - t87 ^ 2) * t79 * t136;
t37 = t137 * t89;
t32 = t37 * t82;
t20 = 0.2e1 * t24 * t82;
t19 = t82 * t138;
t11 = t44 * t138;
t7 = 0.2e1 * t23 * t43 + 0.2e1 * t44 * t24;
t1 = [0, 0, 0, 0, -0.2e1 * t107, -0.2e1 * t106, -0.2e1 * t91, 0.2e1 * t92, 0.2e1 * t123, 0.2e1 * (-pkin(2) - t134) * t107 + 0.2e1 * t76 * t95, t58, t42, t54, t55, 0, -0.2e1 * t16 * t82 + 0.2e1 * t125, -0.2e1 * t76 * t102 + 0.2e1 * t124, t11, t7, t19, t20, 0, -0.2e1 * t3 * t82 + 0.2e1 * t128, 0.2e1 * t111; 0, 0, 0, 0, -t107, -t106, -t91, t92, t120 + t123, -pkin(2) * t107 + qJ(3) * t95 + t76 * t94, t58, t42, t54, t55, 0, (-t16 - t36) * t82 + t122 + t125, (-qJ(3) - t76) * t102 + t121 + t124, t11, t7, t19, t20, 0, (-t3 - t6) * t82 + t127 + t128, t110 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t120, 0.2e1 * qJ(3) * t94, t58, t42, t54, t55, 0, -0.2e1 * t36 * t82 + 0.2e1 * t122, -0.2e1 * t79 * t98 + 0.2e1 * t121, t11, t7, t19, t20, 0, -0.2e1 * t6 * t82 + 0.2e1 * t127, 0.2e1 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, t73, t74, 0, 0, 0, 0, 0, t133, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t74, 0, 0, 0, 0, 0, t133, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t99, 0, t16, t15, 0, 0, -t23, -t24, 0, t63 + t3, t64 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t99, 0, t36, t35, 0, 0, -t23, -t24, 0, t63 + t6, t64 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t114, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t105, -0.2e1 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
