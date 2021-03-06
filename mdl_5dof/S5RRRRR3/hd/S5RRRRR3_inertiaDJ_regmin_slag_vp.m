% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:29
% EndTime: 2019-12-05 18:56:33
% DurationCPUTime: 1.20s
% Computational Cost: add. (1368->136), mult. (3541->249), div. (0->0), fcn. (3509->8), ass. (0->109)
t137 = qJD(2) + qJD(3);
t70 = sin(qJ(3));
t71 = sin(qJ(2));
t74 = cos(qJ(3));
t75 = cos(qJ(2));
t51 = t70 * t71 - t74 * t75;
t40 = t137 * t51;
t69 = sin(qJ(4));
t124 = t69 * t40;
t53 = t70 * t75 + t71 * t74;
t73 = cos(qJ(4));
t65 = qJD(4) * t73;
t84 = t53 * t65 - t124;
t136 = qJD(4) + qJD(5);
t68 = sin(qJ(5));
t126 = t68 * t69;
t72 = cos(qJ(5));
t50 = -t72 * t73 + t126;
t139 = t136 * t50;
t138 = t84 * pkin(3);
t122 = t72 * t69;
t52 = t68 * t73 + t122;
t29 = t52 * t53;
t111 = qJD(4) * t69;
t121 = t73 * t40;
t83 = t111 * t53 + t121;
t67 = t73 ^ 2;
t114 = t69 ^ 2 - t67;
t91 = t114 * qJD(4);
t135 = 2 * pkin(1);
t41 = t137 * t53;
t134 = t41 * pkin(3);
t133 = t51 * pkin(3);
t132 = t74 * pkin(1);
t131 = t50 * t41;
t130 = t52 * t41;
t129 = t53 * t40;
t128 = t53 * t69;
t127 = t53 * t73;
t109 = t71 * qJD(2);
t17 = pkin(1) * t109 + pkin(2) * t41 + pkin(5) * t40;
t125 = t69 * t17;
t123 = t69 * t41;
t15 = t73 * t17;
t120 = t73 * t41;
t39 = t136 * t52;
t113 = pkin(1) * qJD(3);
t101 = t70 * t113;
t99 = pkin(3) * t111;
t54 = t99 + t101;
t64 = -pkin(3) * t73 - pkin(2);
t57 = t64 - t132;
t119 = t39 * t57 + t50 * t54;
t118 = -t139 * t57 + t52 * t54;
t117 = t39 * t64 + t50 * t99;
t116 = -t139 * t64 + t52 * t99;
t63 = -pkin(2) - t132;
t115 = t101 * t69 + t63 * t65;
t112 = pkin(3) * qJD(5);
t110 = qJD(5) * t69;
t108 = t75 * qJD(2);
t107 = pkin(3) * t128;
t105 = t39 * t107 + t138 * t50;
t104 = -t107 * t139 + t138 * t52;
t103 = pkin(2) * t111;
t102 = pkin(2) * t65;
t100 = t74 * t113;
t98 = t68 * t112;
t97 = t72 * t112;
t35 = -pkin(1) * t75 + pkin(2) * t51 - pkin(5) * t53;
t94 = t35 * t110;
t93 = t69 * t65;
t92 = -0.4e1 * t69 * t127;
t89 = -t111 * t35 + t15;
t62 = pkin(1) * t70 + pkin(5);
t88 = t51 * t62 - t53 * t63;
t87 = -t101 * t73 + t111 * t63;
t86 = -t94 + t134;
t85 = -t35 * t65 - t125;
t19 = t51 * t65 + t123;
t82 = t111 * t51 - t120;
t20 = t35 * t73 + t133;
t79 = -qJD(5) * t20 + t85;
t78 = (-t20 - t133) * qJD(5) + t85;
t77 = -t40 * t63 - t41 * t62 + (-t51 * t74 + t53 * t70) * t113;
t59 = 0.2e1 * t93;
t49 = -0.2e1 * t91;
t48 = t53 ^ 2;
t34 = t139 * pkin(5);
t33 = t39 * pkin(5);
t30 = t50 * t53;
t26 = -0.2e1 * t52 * t139;
t25 = 0.2e1 * t51 * t41;
t16 = -t121 * t69 - t53 * t91;
t14 = -t100 * t52 + t139 * t62;
t13 = t100 * t50 + t39 * t62;
t12 = -t39 * t51 - t131;
t11 = -t139 * t51 + t130;
t10 = qJD(4) * t92 + t114 * t40;
t9 = 0.2e1 * t139 * t50 - 0.2e1 * t39 * t52;
t8 = t89 + t134;
t7 = t72 * t8;
t6 = (t127 * t136 - t124) * t72 + (-t110 * t53 - t83) * t68;
t5 = -t136 * t29 + t50 * t40;
t4 = t139 * t30 + t5 * t52;
t3 = t68 * t79 - t72 * t94 + t7;
t2 = (-t8 + t94) * t68 + t79 * t72;
t1 = t139 * t29 + t30 * t39 - t5 * t50 - t52 * t6;
t18 = [0, 0, 0, 0.2e1 * t71 * t108, 0.2e1 * (-t71 ^ 2 + t75 ^ 2) * qJD(2), 0, 0, 0, 0, 0, -0.2e1 * t129, 0.2e1 * t40 * t51 - 0.2e1 * t41 * t53, 0, 0, 0, (t109 * t51 - t41 * t75) * t135, (t109 * t53 + t40 * t75) * t135, -0.2e1 * t129 * t67 - 0.2e1 * t48 * t93, -t40 * t92 + 0.2e1 * t48 * t91, 0.2e1 * t120 * t53 - 0.2e1 * t51 * t83, -0.2e1 * t123 * t53 - 0.2e1 * t51 * t84, t25, 0.2e1 * t15 * t51 - 0.2e1 * t35 * t82, -0.2e1 * t125 * t51 - 0.2e1 * t19 * t35, -0.2e1 * t30 * t5, -0.2e1 * t29 * t5 + 0.2e1 * t30 * t6, -0.2e1 * t30 * t41 + 0.2e1 * t5 * t51, -0.2e1 * t29 * t41 - 0.2e1 * t51 * t6, t25, 0.2e1 * t3 * t51 + 0.2e1 * (-t126 * t35 + t20 * t72) * t41 + 0.2e1 * (t128 * t6 + t29 * t84) * pkin(3), 0.2e1 * t2 * t51 - 0.2e1 * (t122 * t35 + t20 * t68) * t41 + 0.2e1 * (t128 * t5 - t30 * t84) * pkin(3); 0, 0, 0, 0, 0, t108, -t109, 0, 0, 0, 0, 0, -t40, -t41, 0, 0, 0, t16, t10, t19, -t82, 0, -t65 * t88 + t69 * t77, t111 * t88 + t73 * t77, t4, t1, t11, t12, 0, -t130 * t62 + t14 * t51 + t54 * t29 + t57 * t6 + t105, t13 * t51 + t131 * t62 - t54 * t30 + t57 * t5 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t101, -0.2e1 * t100, t59, t49, 0, 0, 0, 0.2e1 * t87, 0.2e1 * t115, t26, t9, 0, 0, 0, 0.2e1 * t119, 0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t41, 0, 0, 0, t16, t10, t19, -t82, 0, -pkin(2) * t84 - pkin(5) * t19, pkin(2) * t83 + pkin(5) * t82, t4, t1, t11, t12, 0, -pkin(5) * t130 + t29 * t99 + t34 * t51 + t64 * t6 + t105, pkin(5) * t131 - t30 * t99 + t33 * t51 + t64 * t5 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100, t59, t49, 0, 0, 0, t87 - t103, -t102 + t115, t26, t9, 0, 0, 0, t117 + t119, t116 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t49, 0, 0, 0, -0.2e1 * t103, -0.2e1 * t102, t26, t9, 0, 0, 0, 0.2e1 * t117, 0.2e1 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t84, t41, t89, t85, 0, 0, t5, -t6, t41, t68 * t78 + t72 * t86 + t7, (-t8 - t86) * t68 + t78 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t111, 0, -t100 * t69 - t62 * t65, -t100 * t73 + t111 * t62, 0, 0, -t139, -t39, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t111, 0, -pkin(5) * t65, pkin(5) * t111, 0, 0, -t139, -t39, 0, t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t98, -0.2e1 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, t41, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t39, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t39, 0, t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
