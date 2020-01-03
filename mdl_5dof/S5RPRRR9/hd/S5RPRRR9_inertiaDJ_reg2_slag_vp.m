% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:00
% EndTime: 2019-12-31 19:08:05
% DurationCPUTime: 1.47s
% Computational Cost: add. (3219->148), mult. (7004->277), div. (0->0), fcn. (7186->8), ass. (0->105)
t68 = sin(qJ(5));
t64 = t68 ^ 2;
t70 = cos(qJ(5));
t65 = t70 ^ 2;
t116 = t64 - t65;
t134 = t116 * qJD(5);
t126 = cos(qJ(3));
t102 = qJD(3) * t126;
t119 = pkin(6) + qJ(2);
t67 = cos(pkin(9));
t51 = t119 * t67;
t66 = sin(pkin(9));
t124 = sin(qJ(3));
t91 = t124 * t119;
t98 = t124 * qJD(2);
t99 = t126 * qJD(2);
t26 = -t51 * t102 - (-qJD(3) * t91 + t99) * t66 - t67 * t98;
t101 = qJD(3) * t124;
t117 = t67 * t101 + t66 * t102;
t139 = 0.2e1 * t117;
t138 = t64 + t65;
t125 = cos(qJ(4));
t100 = qJD(4) * t125;
t69 = sin(qJ(4));
t114 = qJD(4) * t69;
t103 = t126 * t66;
t45 = t119 * t103;
t25 = qJD(3) * t45 + t51 * t101 + t66 * t98 - t67 * t99;
t23 = -t117 * pkin(7) - t25;
t34 = -t124 * t51 - t45;
t47 = t124 * t67 + t103;
t28 = -t47 * pkin(7) + t34;
t35 = t126 * t51 - t66 * t91;
t83 = t124 * t66 - t126 * t67;
t29 = -t83 * pkin(7) + t35;
t79 = t83 * qJD(3);
t72 = -pkin(7) * t79 - t26;
t71 = -t28 * t100 + t29 * t114 - t125 * t23 + t69 * t72;
t76 = t125 * t83;
t32 = t69 * t47 + t76;
t80 = t69 * t83;
t33 = t125 * t47 - t80;
t58 = -t67 * pkin(2) - pkin(1);
t36 = t83 * pkin(3) + t58;
t74 = t32 * pkin(4) - t33 * pkin(8) + t36;
t137 = -qJD(5) * t74 + t71;
t18 = t125 * t29 + t69 * t28;
t8 = -t68 * t18 + t70 * t74;
t9 = t70 * t18 + t68 * t74;
t135 = -t68 * t8 + t70 * t9;
t112 = qJD(5) * t68;
t104 = t117 * pkin(3);
t20 = t47 * t114 + t69 * t117 - (-qJD(3) - qJD(4)) * t76;
t21 = -qJD(4) * t80 + t47 * t100 + t125 * t117 - t69 * t79;
t81 = t21 * pkin(4) + t20 * pkin(8) + t104;
t2 = t18 * t112 + t137 * t70 - t68 * t81;
t61 = qJD(5) * t70;
t3 = t137 * t68 - t18 * t61 + t70 * t81;
t133 = -qJD(5) * t135 + t2 * t68 - t3 * t70;
t17 = -t125 * t28 + t69 * t29;
t6 = t18 * qJD(4) + t125 * t72 + t69 * t23;
t132 = t17 * t6;
t131 = t6 * t33;
t5 = t6 * t68;
t129 = t69 * pkin(3);
t127 = t17 * t61 + t5;
t123 = t33 * t20;
t122 = t68 * t21;
t121 = t70 * t20;
t120 = t70 * t21;
t107 = pkin(3) * t114;
t106 = t125 * pkin(3);
t60 = -t106 - pkin(4);
t118 = t68 * t107 + t60 * t61;
t115 = pkin(3) * qJD(4);
t111 = 0.2e1 * t32 * t21;
t110 = t68 * t121;
t109 = pkin(4) * t112;
t108 = pkin(4) * t61;
t105 = t68 * t61;
t31 = t33 ^ 2;
t97 = t31 * t105;
t96 = pkin(3) * t100;
t95 = 0.2e1 * (t66 ^ 2 + t67 ^ 2) * qJD(2);
t93 = t68 * t9 + t70 * t8;
t59 = pkin(8) + t129;
t90 = t32 * t59 - t33 * t60;
t89 = -t70 * t107 + t60 * t112;
t87 = -t68 * t20 + t33 * t61;
t86 = t33 * t112 + t121;
t13 = t32 * t61 + t122;
t85 = t32 * t112 - t120;
t84 = t138 * t125;
t82 = (-t125 * t32 + t33 * t69) * qJD(4);
t78 = -0.2e1 * t79;
t1 = -t93 * qJD(5) - t2 * t70 - t3 * t68;
t75 = pkin(3) * t82 - t20 * t60 - t21 * t59;
t54 = -0.2e1 * t105;
t53 = 0.2e1 * t105;
t48 = -0.2e1 * t134;
t41 = t84 * t115;
t14 = t17 * t112;
t11 = t33 * t134 + t110;
t7 = -0.4e1 * t33 * t105 + t116 * t20;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJ(2) * t95, t47 * t78, 0.2e1 * t83 ^ 2 * qJD(3) - 0.2e1 * t47 * t117, 0, t83 * t139, 0, 0, t58 * t139, t58 * t78, -0.2e1 * t35 * t117 - 0.2e1 * t26 * t47 + 0.2e1 * (qJD(3) * t34 + t25) * t83, -0.2e1 * t35 * t25 + 0.2e1 * t34 * t26, -0.2e1 * t123, 0.2e1 * t20 * t32 - 0.2e1 * t33 * t21, 0, t111, 0, 0, 0.2e1 * t32 * t104 + 0.2e1 * t36 * t21, 0.2e1 * t33 * t104 - 0.2e1 * t36 * t20, -0.2e1 * t17 * t20 - 0.2e1 * t18 * t21 + 0.2e1 * t32 * t71 + 0.2e1 * t131, 0.2e1 * t36 * t104 - 0.2e1 * t18 * t71 + 0.2e1 * t132, -0.2e1 * t65 * t123 - 0.2e1 * t97, 0.4e1 * t33 * t110 + 0.2e1 * t31 * t134, 0.2e1 * t33 * t120 - 0.2e1 * t86 * t32, -0.2e1 * t64 * t123 + 0.2e1 * t97, -0.2e1 * t33 * t122 - 0.2e1 * t87 * t32, t111, 0.2e1 * t87 * t17 + 0.2e1 * t8 * t21 + 0.2e1 * t3 * t32 + 0.2e1 * t33 * t5, 0.2e1 * t70 * t131 - 0.2e1 * t86 * t17 + 0.2e1 * t2 * t32 - 0.2e1 * t9 * t21, 0.2e1 * t133 * t33 + 0.2e1 * t93 * t20, -0.2e1 * t9 * t2 + 0.2e1 * t8 * t3 + 0.2e1 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, -t79, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t104, 0, 0, 0, 0, 0, 0, -t85, -t13, t138 * t20, -t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, -t117, 0, t26, t25, 0, 0, 0, 0, -t20, 0, -t21, 0, -t6, t71, (t125 * t20 - t21 * t69 + t82) * pkin(3), -t6 * t106 + t17 * t107 - t71 * t129 + t18 * t96, -t11, t7, t13, t11, -t85, 0, t14 + (-t90 * qJD(5) - t6) * t70 + t75 * t68, t90 * t112 + t75 * t70 + t127, t1, t6 * t60 + (t135 * t125 + t17 * t69) * t115 + t1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t107, -0.2e1 * t96, 0, 0, t53, t48, 0, t54, 0, 0, 0.2e1 * t89, 0.2e1 * t118, 0.2e1 * t41, 0.2e1 * (t84 * t59 + t60 * t69) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, 0, -t6, t71, 0, 0, -t11, t7, t13, t11, -t85, 0, t14 + (pkin(4) * t20 - pkin(8) * t21) * t68 + (-t6 + (-pkin(4) * t33 - pkin(8) * t32) * qJD(5)) * t70, t86 * pkin(4) + t85 * pkin(8) + t127, t1, -t6 * pkin(4) + t1 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t96, 0, 0, t53, t48, 0, t54, 0, 0, t89 - t109, -t108 + t118, t41, (-pkin(4) * t69 + t84 * pkin(8)) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t48, 0, t54, 0, 0, -0.2e1 * t109, -0.2e1 * t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, -t87, t21, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t112, 0, -t59 * t61 - t68 * t96, t59 * t112 - t70 * t96, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t112, 0, -pkin(8) * t61, pkin(8) * t112, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
