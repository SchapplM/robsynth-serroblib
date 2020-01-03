% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:31
% DurationCPUTime: 1.32s
% Computational Cost: add. (1901->168), mult. (4239->289), div. (0->0), fcn. (3752->6), ass. (0->113)
t61 = sin(qJ(5));
t59 = t61 ^ 2;
t64 = cos(qJ(5));
t60 = t64 ^ 2;
t115 = t59 - t60;
t132 = t115 * qJD(5);
t131 = qJD(2) + qJD(3);
t67 = 2 * qJD(4);
t130 = pkin(3) + pkin(8);
t129 = -pkin(7) - pkin(6);
t110 = qJD(5) * t61;
t62 = sin(qJ(3));
t63 = sin(qJ(2));
t120 = t62 * t63;
t65 = cos(qJ(2));
t124 = cos(qJ(3));
t95 = qJD(3) * t124;
t97 = t124 * t65;
t25 = -qJD(2) * t97 + t120 * t131 - t65 * t95;
t127 = t25 * pkin(4);
t43 = t129 * t65;
t91 = t129 * t124;
t84 = qJD(2) * t91;
t81 = -t43 * t95 - t65 * t84;
t96 = t129 * qJD(2);
t16 = (qJD(3) * t129 + t96) * t120 + t81;
t37 = -t97 + t120;
t119 = t62 * t65;
t38 = t124 * t63 + t119;
t54 = -t65 * pkin(2) - pkin(1);
t83 = -t38 * qJ(4) + t54;
t19 = t130 * t37 + t83;
t26 = t131 * t38;
t107 = t63 * qJD(2);
t101 = pkin(2) * t107;
t87 = t25 * qJ(4) - t38 * qJD(4);
t72 = t87 + t101;
t39 = t63 * t91;
t27 = -t43 * t62 - t39;
t82 = pkin(4) * t38 + t27;
t75 = t64 * t82;
t2 = t19 * t110 - qJD(5) * t75 - t64 * (t130 * t26 + t72) - t61 * (t16 - t127);
t128 = t2 * t61;
t126 = t26 * pkin(3);
t109 = qJD(5) * t64;
t100 = t62 * t129;
t28 = t100 * t63 - t124 * t43;
t23 = -t37 * pkin(4) + t28;
t113 = qJD(3) * t62;
t15 = -qJD(3) * t39 - t113 * t43 - t96 * t119 - t63 * t84;
t11 = -pkin(4) * t26 - t15;
t9 = t11 * t61;
t125 = t109 * t23 + t9;
t10 = t11 * t64;
t123 = t37 * t26;
t92 = pkin(2) * t95;
t46 = t92 + qJD(4);
t50 = pkin(2) * t62 + qJ(4);
t122 = t50 * t46;
t121 = t61 * t26;
t118 = t64 * t26;
t117 = t50 * t109 + t46 * t61;
t105 = qJ(4) * qJD(5);
t116 = qJD(4) * t61 + t64 * t105;
t7 = t64 * t19 + t61 * t82;
t114 = qJD(5) * t7;
t112 = qJD(5) * t23;
t108 = qJD(5) * t130;
t106 = t65 * qJD(2);
t104 = 0.2e1 * t123;
t22 = -0.2e1 * t38 * t25;
t103 = -0.2e1 * pkin(1) * qJD(2);
t102 = t61 * t118;
t55 = pkin(2) * t113;
t99 = t63 * t106;
t98 = t61 * t109;
t94 = t64 * t100;
t35 = t37 ^ 2;
t93 = t35 * t98;
t53 = -pkin(2) * t124 - pkin(3);
t6 = -t19 * t61 + t75;
t90 = t6 * t64 + t61 * t7;
t89 = -t6 * t61 + t64 * t7;
t88 = -t15 * t28 + t16 * t27;
t30 = (t59 + t60) * t55;
t86 = -qJ(4) * t26 - qJD(4) * t37;
t85 = t46 * qJ(4) + t50 * qJD(4);
t80 = t109 * t38 - t25 * t61;
t17 = -t110 * t38 - t25 * t64;
t79 = t109 * t37 + t121;
t78 = t110 * t37 - t118;
t49 = -pkin(8) + t53;
t77 = qJD(5) * (t37 * t50 - t38 * t49);
t76 = qJD(5) * (qJ(4) * t37 + t130 * t38);
t74 = 0.2e1 * t25 * t37 - 0.2e1 * t26 * t38;
t73 = t130 * t25 + t86;
t71 = -t26 * t50 - t37 * t46 + t38 * t55;
t3 = -t61 * (t26 * pkin(8) + t126 + t87) + t64 * (t81 - t127) - t114 + (qJD(3) * t94 + (-t61 * pkin(2) + t94) * qJD(2)) * t63;
t70 = qJD(5) * t89 + t3 * t64 - t128;
t69 = -t25 * t49 + t71;
t68 = 0.2e1 * t15 * t37 + 0.2e1 * t16 * t38 - 0.2e1 * t25 * t27 - 0.2e1 * t26 * t28;
t58 = qJ(4) * t67;
t57 = qJD(4) * t64;
t45 = -0.2e1 * t98;
t44 = 0.2e1 * t98;
t41 = t46 * t64;
t36 = 0.2e1 * t132;
t24 = pkin(3) * t37 + t83;
t14 = -t132 * t37 + t102;
t13 = t72 + t126;
t12 = -t115 * t26 - 0.4e1 * t37 * t98;
t1 = -t6 * t110 - t128 + (t3 + t114) * t64;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t99, 0.2e1 * (-t63 ^ 2 + t65 ^ 2) * qJD(2), 0, -0.2e1 * t99, 0, 0, t63 * t103, t65 * t103, 0, 0, t22, t74, 0, t104, 0, 0, 0.2e1 * t101 * t37 + 0.2e1 * t26 * t54, 0.2e1 * t101 * t38 - 0.2e1 * t25 * t54, t68, 0.2e1 * t101 * t54 + 0.2e1 * t88, 0, 0, 0, t22, t74, t104, t68, -0.2e1 * t13 * t37 - 0.2e1 * t24 * t26, -0.2e1 * t13 * t38 + 0.2e1 * t24 * t25, 0.2e1 * t13 * t24 + 0.2e1 * t88, 0.2e1 * t123 * t59 + 0.2e1 * t93, 0.4e1 * t102 * t37 - 0.2e1 * t132 * t35, 0.2e1 * t121 * t38 + 0.2e1 * t37 * t80, 0.2e1 * t123 * t60 - 0.2e1 * t93, 0.2e1 * t118 * t38 + 0.2e1 * t17 * t37, t22, -0.2e1 * t10 * t37 + 0.2e1 * t23 * t78 - 0.2e1 * t6 * t25 + 0.2e1 * t3 * t38, 0.2e1 * t2 * t38 + 0.2e1 * t23 * t79 + 0.2e1 * t7 * t25 + 0.2e1 * t37 * t9, 0.2e1 * t89 * t26 + 0.2e1 * (-qJD(5) * t90 - t2 * t64 - t3 * t61) * t37, 0.2e1 * t11 * t23 - 0.2e1 * t2 * t7 + 0.2e1 * t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, -t107, 0, -pkin(6) * t106, pkin(6) * t107, 0, 0, 0, 0, -t25, 0, -t26, 0, -t16, t15, (t124 * t25 - t26 * t62 + (-t124 * t37 + t38 * t62) * qJD(3)) * pkin(2), (-t124 * t16 - t15 * t62 + (t124 * t28 + t27 * t62) * qJD(3)) * pkin(2), 0, t25, t26, 0, 0, 0, -t25 * t53 + t71, t16, -t15, -t15 * t50 + t16 * t53 + t27 * t55 + t28 * t46, t14, t12, t17, -t14, -t80, 0, t61 * t77 + t64 * t69 + t125, t10 + t64 * t77 + (-t69 - t112) * t61, -t1, t11 * t50 + t23 * t46 + t49 * t70 + t55 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, -0.2e1 * t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55, 0.2e1 * t46, 0.2e1 * t53 * t55 + 0.2e1 * t122, t45, t36, 0, t44, 0, 0, 0.2e1 * t117, -0.2e1 * t110 * t50 + 0.2e1 * t41, -0.2e1 * t30, 0.2e1 * t30 * t49 + 0.2e1 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, -t26, 0, -t16, t15, 0, 0, 0, t25, t26, 0, 0, 0, pkin(3) * t25 + t86, t16, -t15, -pkin(3) * t16 - qJ(4) * t15 + qJD(4) * t28, t14, t12, t17, -t14, -t80, 0, t61 * t76 + t64 * t73 + t125, t10 + t64 * t76 + (-t73 - t112) * t61, -t1, t11 * qJ(4) + t23 * qJD(4) - t130 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t67 + t92, -pkin(3) * t55 + t85, t45, t36, 0, t44, 0, 0, t116 + t117, t41 + t57 + (-qJ(4) - t50) * t110, -t30, -t130 * t30 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t58, t45, t36, 0, t44, 0, 0, 0.2e1 * t116, -0.2e1 * t105 * t61 + 0.2e1 * t57, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, t16, 0, 0, 0, 0, 0, 0, t17, -t80, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t78, -t25, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, -t109, 0, -t110 * t49 + t55 * t64, -t109 * t49 - t55 * t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, -t109, 0, t61 * t108, t64 * t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
