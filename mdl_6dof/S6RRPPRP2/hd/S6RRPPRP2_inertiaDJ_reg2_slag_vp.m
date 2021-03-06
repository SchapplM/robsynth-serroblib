% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:57
% EndTime: 2019-03-09 08:32:03
% DurationCPUTime: 1.94s
% Computational Cost: add. (2523->197), mult. (5416->329), div. (0->0), fcn. (5018->6), ass. (0->116)
t136 = pkin(3) + pkin(8);
t124 = sin(pkin(9));
t125 = cos(pkin(9));
t75 = sin(qJ(2));
t77 = cos(qJ(2));
t58 = t124 * t77 + t125 * t75;
t52 = t58 * qJD(2);
t129 = -qJ(3) - pkin(7);
t42 = t58 * t129;
t86 = t58 * pkin(4) - t42;
t120 = t75 * qJD(2);
t112 = pkin(2) * t120;
t57 = t124 * t75 - t125 * t77;
t53 = t57 * qJD(2);
t87 = t53 * qJ(4) - t58 * qJD(4) + t112;
t139 = -qJD(5) * t86 - t136 * t52 - t87;
t74 = sin(qJ(5));
t72 = t74 ^ 2;
t76 = cos(qJ(5));
t73 = t76 ^ 2;
t128 = t72 - t73;
t122 = qJD(5) * t76;
t131 = t74 * t52;
t37 = t57 * t122 + t131;
t138 = t128 * qJD(5);
t137 = 0.2e1 * qJD(4);
t135 = t53 * pkin(5);
t134 = t57 * t74;
t133 = t57 * t76;
t132 = t72 * t52;
t48 = t73 * t52;
t130 = t76 * t52;
t70 = -t77 * pkin(2) - pkin(1);
t94 = -t58 * qJ(4) + t70;
t26 = t136 * t57 + t94;
t11 = t76 * t26 + t74 * t86;
t127 = t76 * qJ(6);
t69 = -t125 * pkin(2) - pkin(3);
t66 = -pkin(8) + t69;
t126 = qJ(6) - t66;
t123 = qJD(5) * t57;
t71 = qJD(5) * t74;
t121 = t74 * qJD(6);
t119 = t76 * qJD(6);
t118 = t77 * qJD(2);
t117 = 0.2e1 * t57 * t52;
t41 = -0.2e1 * t58 * t53;
t116 = -0.2e1 * pkin(1) * qJD(2);
t115 = t74 * t130;
t67 = t124 * pkin(2) + qJ(4);
t114 = t67 * t137;
t101 = qJD(2) * t129;
t89 = t77 * qJD(3) + t75 * t101;
t90 = -t75 * qJD(3) + t77 * t101;
t30 = t124 * t89 - t125 * t90;
t79 = -t53 * pkin(4) + t30;
t113 = t139 * t76 - t74 * t79;
t111 = pkin(5) * t71;
t110 = pkin(5) * t122;
t108 = t75 * t118;
t107 = t74 * t122;
t55 = t126 * t76;
t102 = qJ(6) * t57 + t26;
t56 = t57 ^ 2;
t100 = t56 * t107;
t31 = t124 * t90 + t125 * t89;
t19 = -t52 * pkin(4) + t31;
t34 = t57 * t71 - t130;
t12 = t34 * pkin(5) + t19;
t61 = t74 * pkin(5) + t67;
t99 = t61 * t123 + t12;
t29 = t76 * t86;
t6 = t58 * pkin(5) - t102 * t74 + t29;
t7 = t57 * t127 + t11;
t98 = -t6 * t74 + t7 * t76;
t10 = -t74 * t26 + t29;
t97 = -t10 * t74 + t11 * t76;
t43 = t57 * t129;
t96 = -t42 * t30 + t43 * t31;
t95 = -qJD(4) * t57 - t67 * t52;
t93 = t58 * t122 - t74 * t53;
t36 = t76 * t53 + t58 * t71;
t92 = -0.2e1 * t58 * t52 + 0.2e1 * t53 * t57;
t20 = (-t76 * pkin(5) - pkin(4)) * t57 + t43;
t64 = qJD(4) + t110;
t91 = qJD(5) * t20 - t52 * t61 - t57 * t64;
t88 = t19 + (t57 * t67 - t58 * t66) * qJD(5);
t32 = -t57 * pkin(4) + t43;
t85 = -qJD(5) * t32 + t53 * t66 - t95;
t5 = -t26 * t122 + t139 * t74 + t76 * t79;
t78 = -t37 * qJ(6) - t57 * t121 + t5;
t2 = t78 - t135;
t3 = t102 * t71 - t57 * t119 - t52 * t127 + t113;
t83 = -t2 * t74 - t3 * t76 + (-t6 * t76 - t7 * t74) * qJD(5);
t4 = t26 * t71 + t113;
t1 = t97 * qJD(5) - t4 * t74 + t5 * t76;
t82 = -t4 * t76 - t5 * t74 + (-t10 * t76 - t11 * t74) * qJD(5);
t81 = 0.2e1 * t30 * t58 - 0.2e1 * t31 * t57 + 0.2e1 * t42 * t53 - 0.2e1 * t43 * t52;
t63 = -0.2e1 * t107;
t62 = 0.2e1 * t107;
t60 = 0.2e1 * t138;
t54 = t126 * t74;
t45 = -qJD(5) * t55 - t121;
t44 = t126 * t71 - t119;
t40 = t57 * pkin(3) + t94;
t33 = t48 + t132;
t24 = 0.2e1 * t57 * t48 - 0.2e1 * t100;
t23 = 0.2e1 * t57 * t132 + 0.2e1 * t100;
t22 = t52 * pkin(3) + t87;
t21 = -t128 * t123 + t115;
t18 = -0.4e1 * t57 * t107 - t132 + t48;
t15 = 0.4e1 * t57 * t115 - 0.2e1 * t56 * t138;
t13 = t44 * t76 + t45 * t74 + (-t54 * t76 + t55 * t74) * qJD(5);
t9 = 0.2e1 * t58 * t131 + 0.2e1 * t93 * t57;
t8 = 0.2e1 * t58 * t130 - 0.2e1 * t36 * t57;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t108, 0.2e1 * (-t75 ^ 2 + t77 ^ 2) * qJD(2), 0, -0.2e1 * t108, 0, 0, t75 * t116, t77 * t116, 0, 0, t41, t92, 0, t117, 0, 0, 0.2e1 * t57 * t112 + 0.2e1 * t70 * t52, 0.2e1 * t58 * t112 - 0.2e1 * t70 * t53, t81, 0.2e1 * t70 * t112 + 0.2e1 * t96, 0, 0, 0, t41, t92, t117, t81, -0.2e1 * t22 * t57 - 0.2e1 * t40 * t52, -0.2e1 * t22 * t58 + 0.2e1 * t40 * t53, 0.2e1 * t40 * t22 + 0.2e1 * t96, t23, t15, t9, t24, t8, t41, -0.2e1 * t10 * t53 - 0.2e1 * t19 * t133 + 0.2e1 * t34 * t32 + 0.2e1 * t5 * t58, 0.2e1 * t11 * t53 + 0.2e1 * t19 * t134 + 0.2e1 * t37 * t32 + 0.2e1 * t4 * t58, 0.2e1 * t97 * t52 + 0.2e1 * t82 * t57, 0.2e1 * t10 * t5 - 0.2e1 * t11 * t4 + 0.2e1 * t32 * t19, t23, t15, t9, t24, t8, t41, -0.2e1 * t12 * t133 + 0.2e1 * t2 * t58 + 0.2e1 * t34 * t20 - 0.2e1 * t6 * t53, 0.2e1 * t12 * t134 + 0.2e1 * t37 * t20 + 0.2e1 * t3 * t58 + 0.2e1 * t7 * t53, 0.2e1 * t98 * t52 + 0.2e1 * t83 * t57, 0.2e1 * t20 * t12 + 0.2e1 * t6 * t2 - 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, -t120, 0, -pkin(7) * t118, pkin(7) * t120, 0, 0, 0, 0, -t53, 0, -t52, 0, -t30, -t31 (-t124 * t52 + t125 * t53) * pkin(2) (t124 * t31 - t125 * t30) * pkin(2), 0, t53, t52, 0, 0, 0, -t69 * t53 + t95, t30, t31, t43 * qJD(4) + t30 * t69 + t31 * t67, t21, t18, -t36, -t21, -t93, 0, t74 * t88 - t76 * t85, t74 * t85 + t76 * t88, -t1, t32 * qJD(4) + t1 * t66 + t19 * t67, t21, t18, -t36, -t21, -t93, 0, t44 * t58 + t55 * t53 + t99 * t74 + t91 * t76, -t45 * t58 - t54 * t53 - t91 * t74 + t99 * t76 (t45 * t57 - t52 * t54 - t2 + (t55 * t57 - t7) * qJD(5)) * t76 + (-t44 * t57 + t52 * t55 + t3 + (t54 * t57 + t6) * qJD(5)) * t74, t12 * t61 - t2 * t55 + t20 * t64 + t3 * t54 + t6 * t44 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t114, t63, t60, 0, t62, 0, 0, 0.2e1 * qJD(4) * t74 + 0.2e1 * t67 * t122, 0.2e1 * qJD(4) * t76 - 0.2e1 * t67 * t71, 0, t114, t63, t60, 0, t62, 0, 0, 0.2e1 * t61 * t122 + 0.2e1 * t64 * t74, -0.2e1 * t61 * t71 + 0.2e1 * t64 * t76, -0.2e1 * t13, -0.2e1 * t55 * t44 - 0.2e1 * t54 * t45 + 0.2e1 * t61 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t53, 0, t112, 0, 0, 0, 0, 0, 0, 0, -t52, t53, t22, 0, 0, 0, 0, 0, 0, -t93, t36, t33, t82, 0, 0, 0, 0, 0, 0, -t93, t36, t33, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t74 + t45 * t76 + (t54 * t74 + t55 * t76) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, t30, 0, 0, 0, 0, 0, 0, -t36, -t93, 0, t1, 0, 0, 0, 0, 0, 0, -t36, -t93, 0, t98 * qJD(5) + t2 * t76 - t3 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t34, -t53, t5, t4, 0, 0, 0, 0, t37, 0, -t34, -t53, t78 - 0.2e1 * t135, t3, -t37 * pkin(5), t2 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, -t122, 0, -t66 * t71, -t66 * t122, 0, 0, 0, 0, -t71, 0, -t122, 0, t44, -t45, t111, t44 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t71, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t71, 0, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t122, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t122, 0, -t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t37, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t71, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t14;
