% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:17
% EndTime: 2019-12-31 21:14:23
% DurationCPUTime: 1.72s
% Computational Cost: add. (3560->169), mult. (7963->324), div. (0->0), fcn. (7571->8), ass. (0->116)
t140 = -pkin(7) - pkin(6);
t76 = sin(qJ(2));
t111 = t140 * t76;
t79 = cos(qJ(2));
t143 = t140 * t79;
t75 = sin(qJ(3));
t132 = t75 * t143;
t78 = cos(qJ(3));
t37 = t78 * t111 + t132;
t74 = sin(qJ(5));
t71 = t74 ^ 2;
t77 = cos(qJ(5));
t72 = t77 ^ 2;
t127 = t71 - t72;
t142 = t127 * qJD(5);
t145 = -t75 * t111 + t143 * t78;
t100 = qJD(3) * t111;
t122 = qJD(3) * t78;
t124 = cos(pkin(9));
t115 = t37 * qJD(2) + t78 * t100;
t55 = t75 * t79 + t78 * t76;
t125 = t55 * qJ(4);
t131 = t75 * t76;
t54 = -t78 * t79 + t131;
t88 = t55 * qJD(2);
t16 = -t54 * qJD(4) - qJ(4) * t88 + (-t125 + t132) * qJD(3) + t115;
t73 = sin(pkin(9));
t83 = t145 * qJD(2);
t118 = t79 * qJD(2);
t36 = -t78 * t118 - t79 * t122 + (qJD(2) + qJD(3)) * t131;
t92 = t36 * qJ(4) - t55 * qJD(4);
t7 = t124 * t16 + (-t75 * t100 + t122 * t143 + t83 + t92) * t73;
t34 = t124 * t54 + t73 * t55;
t35 = t124 * t55 - t73 * t54;
t69 = -t79 * pkin(2) - pkin(1);
t42 = t54 * pkin(3) + t69;
t84 = t34 * pkin(4) - t35 * pkin(8) + t42;
t144 = -qJD(5) * t84 - t7;
t120 = qJD(5) * t74;
t29 = -t54 * qJ(4) - t145;
t86 = t37 - t125;
t21 = t124 * t29 + t73 * t86;
t87 = qJD(3) * t55;
t81 = -t87 - t88;
t23 = -t124 * t81 - t73 * t36;
t24 = -t124 * t36 + t73 * t81;
t138 = t76 * pkin(2);
t30 = pkin(3) * t87 + (t55 * pkin(3) + t138) * qJD(2);
t80 = t23 * pkin(4) - t24 * pkin(8) + t30;
t2 = t21 * t120 + t144 * t77 - t74 * t80;
t70 = qJD(5) * t77;
t3 = t144 * t74 - t21 * t70 + t77 * t80;
t8 = -t74 * t21 + t77 * t84;
t9 = t77 * t21 + t74 * t84;
t97 = t74 * t8 - t77 * t9;
t141 = qJD(5) * t97 + t2 * t74 - t3 * t77;
t20 = -t124 * t86 + t73 * t29;
t26 = qJD(3) * t145 + t83;
t6 = t73 * t16 - t124 * (t26 + t92);
t139 = t20 * t6;
t5 = t6 * t74;
t137 = t20 * t70 + t5;
t103 = t124 * t75;
t45 = (t73 * t78 + t103) * qJD(3) * pkin(2);
t136 = t20 * t45;
t135 = t35 * t24;
t134 = t35 * t77;
t133 = t74 * t23;
t130 = t77 * t23;
t129 = t77 * t24;
t68 = t78 * pkin(2) + pkin(3);
t47 = -t73 * t75 * pkin(2) + t124 * t68;
t43 = -pkin(4) - t47;
t128 = t43 * t70 + t45 * t74;
t48 = pkin(2) * t103 + t73 * t68;
t126 = t71 + t72;
t123 = qJD(3) * t75;
t119 = t76 * qJD(2);
t117 = 0.2e1 * t34 * t23;
t116 = -0.2e1 * pkin(1) * qJD(2);
t114 = pkin(2) * t119;
t113 = pkin(2) * t123;
t112 = pkin(2) * t122;
t67 = -t124 * pkin(3) - pkin(4);
t109 = t67 * t120;
t108 = t67 * t70;
t107 = t74 * t70;
t106 = t76 * t118;
t46 = t124 * t112 - t73 * t113;
t28 = t126 * t46;
t105 = -0.4e1 * t74 * t134;
t104 = t43 * t120 - t45 * t77;
t32 = t35 ^ 2;
t101 = t32 * t107;
t98 = -t74 * t9 - t77 * t8;
t66 = t73 * pkin(3) + pkin(8);
t96 = -t23 * t66 + t24 * t67;
t44 = pkin(8) + t48;
t95 = t34 * t44 - t35 * t43;
t94 = -t46 * t34 + t45 * t35;
t93 = t34 * t66 - t35 * t67;
t14 = t34 * t70 + t133;
t91 = t74 * t24 + t35 * t70;
t90 = -t35 * t120 + t129;
t89 = t69 * t55;
t85 = -t23 * t44 + t24 * t43 + t94;
t1 = t98 * qJD(5) - t2 * t77 - t3 * t74;
t62 = -0.2e1 * t107;
t61 = 0.2e1 * t107;
t53 = -0.2e1 * t142;
t25 = -t123 * t143 - t115;
t18 = t20 * t120;
t13 = -t34 * t120 + t130;
t12 = -t74 * t129 + t35 * t142;
t10 = qJD(5) * t105 - t127 * t24;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t106, 0.2e1 * (-t76 ^ 2 + t79 ^ 2) * qJD(2), 0, -0.2e1 * t106, 0, 0, t76 * t116, t79 * t116, 0, 0, -0.2e1 * t55 * t36, 0.2e1 * t36 * t54 + 0.2e1 * t55 * t81, 0, -0.2e1 * t54 * t81, 0, 0, 0.2e1 * qJD(3) * t89 + 0.2e1 * (t54 * t138 + t89) * qJD(2), 0.2e1 * t55 * t114 - 0.2e1 * t69 * t36, -0.2e1 * t145 * t81 + 0.2e1 * t25 * t54 - 0.2e1 * t26 * t55 + 0.2e1 * t37 * t36, 0.2e1 * t69 * t114 + 0.2e1 * t145 * t25 + 0.2e1 * t37 * t26, 0.2e1 * t135, -0.2e1 * t35 * t23 - 0.2e1 * t24 * t34, 0, t117, 0, 0, 0.2e1 * t42 * t23 + 0.2e1 * t30 * t34, 0.2e1 * t42 * t24 + 0.2e1 * t30 * t35, 0.2e1 * t20 * t24 - 0.2e1 * t21 * t23 - 0.2e1 * t7 * t34 + 0.2e1 * t6 * t35, 0.2e1 * t21 * t7 + 0.2e1 * t42 * t30 + 0.2e1 * t139, 0.2e1 * t72 * t135 - 0.2e1 * t101, t24 * t105 + 0.2e1 * t32 * t142, 0.2e1 * t35 * t130 + 0.2e1 * t90 * t34, 0.2e1 * t71 * t135 + 0.2e1 * t101, -0.2e1 * t35 * t133 - 0.2e1 * t91 * t34, t117, 0.2e1 * t91 * t20 + 0.2e1 * t8 * t23 + 0.2e1 * t3 * t34 + 0.2e1 * t35 * t5, 0.2e1 * t6 * t134 + 0.2e1 * t2 * t34 + 0.2e1 * t90 * t20 - 0.2e1 * t9 * t23, 0.2e1 * t141 * t35 + 0.2e1 * t98 * t24, -0.2e1 * t9 * t2 + 0.2e1 * t8 * t3 + 0.2e1 * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, -t119, 0, -pkin(6) * t118, pkin(6) * t119, 0, 0, 0, 0, -t36, 0, t81, 0, t26, t25, (-t54 * t122 + t78 * t36 - t75 * t88) * pkin(2), (-t25 * t75 + t26 * t78 + (-t145 * t78 - t37 * t75) * qJD(3)) * pkin(2), 0, 0, t24, 0, -t23, 0, -t6, -t7, -t48 * t23 - t47 * t24 + t94, t21 * t46 - t6 * t47 + t7 * t48 + t136, -t12, t10, t14, t12, t13, 0, t18 + (-t95 * qJD(5) - t6) * t77 + t85 * t74, t95 * t120 + t85 * t77 + t137, t1, t1 * t44 + t6 * t43 - t46 * t97 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113, -0.2e1 * t112, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t45, -0.2e1 * t46, 0, -0.2e1 * t47 * t45 + 0.2e1 * t48 * t46, t61, t53, 0, t62, 0, 0, 0.2e1 * t104, 0.2e1 * t128, 0.2e1 * t28, 0.2e1 * t28 * t44 + 0.2e1 * t43 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, t81, 0, t26, t25, 0, 0, 0, 0, t24, 0, -t23, 0, -t6, -t7, (-t124 * t24 - t23 * t73) * pkin(3), (-t124 * t6 + t7 * t73) * pkin(3), -t12, t10, t14, t12, t13, 0, t18 + t96 * t74 + (-t93 * qJD(5) - t6) * t77, t93 * t120 + t77 * t96 + t137, t1, t1 * t66 + t6 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, (-t124 * t45 + t46 * t73) * pkin(3), t61, t53, 0, t62, 0, 0, t104 + t109, t108 + t128, t28, t28 * t66 + t45 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t53, 0, t62, 0, 0, 0.2e1 * t109, 0.2e1 * t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t24, 0, t30, 0, 0, 0, 0, 0, 0, t13, -t14, -t126 * t24, -t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, -t91, t23, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t120, 0, -t44 * t70 - t74 * t46, t120 * t44 - t77 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t120, 0, -t66 * t70, t66 * t120, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t70, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
