% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:54
% EndTime: 2019-03-09 07:05:58
% DurationCPUTime: 1.46s
% Computational Cost: add. (4990->168), mult. (10879->277), div. (0->0), fcn. (11919->10), ass. (0->113)
t79 = cos(qJ(5));
t118 = qJD(5) * t79;
t76 = sin(qJ(4));
t128 = t76 * t79;
t75 = sin(qJ(5));
t80 = cos(qJ(4));
t137 = ((t75 * t80 + t128) * qJD(4) + t76 * t118) * pkin(3);
t78 = cos(qJ(6));
t71 = t78 ^ 2;
t74 = sin(qJ(6));
t122 = t74 ^ 2 - t71;
t105 = t122 * qJD(6);
t72 = sin(pkin(11));
t73 = cos(pkin(11));
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t93 = t77 * t72 - t81 * t73;
t125 = pkin(7) + qJ(2);
t58 = t125 * t72;
t59 = t125 * t73;
t95 = -t81 * t58 - t77 * t59;
t135 = t93 * qJD(2) - t95 * qJD(3);
t53 = t81 * t72 + t77 * t73;
t51 = t53 * qJD(3);
t32 = -t51 * pkin(8) - t135;
t50 = t93 * qJD(3);
t94 = t77 * t58 - t81 * t59;
t83 = -t53 * qJD(2) + t94 * qJD(3);
t33 = t50 * pkin(8) + t83;
t35 = -t53 * pkin(8) + t95;
t36 = -pkin(8) * t93 - t94;
t98 = t80 * t35 - t76 * t36;
t20 = -t98 * qJD(4) - t80 * t32 - t76 * t33;
t63 = -t73 * pkin(2) - pkin(1);
t134 = 0.2e1 * t63;
t41 = t80 * t53 - t76 * t93;
t24 = -t41 * pkin(9) + t98;
t96 = -t76 * t53 - t80 * t93;
t97 = -t76 * t35 - t80 * t36;
t25 = t96 * pkin(9) - t97;
t17 = t75 * t24 + t79 * t25;
t120 = qJD(4) * t80;
t121 = qJD(4) * t76;
t91 = t53 * t120 - t121 * t93 - t76 * t50 + t80 * t51;
t82 = -t91 * pkin(9) - t20;
t21 = t97 * qJD(4) - t76 * t32 + t80 * t33;
t29 = t96 * qJD(4) - t80 * t50 - t76 * t51;
t84 = t29 * pkin(9) - t21;
t5 = t17 * qJD(5) + t75 * t82 + t79 * t84;
t3 = t5 * t74;
t133 = t51 * pkin(3);
t16 = -t79 * t24 + t75 * t25;
t67 = qJD(6) * t78;
t132 = t16 * t67 + t3;
t30 = t75 * t41 - t79 * t96;
t18 = -t30 * qJD(5) + t79 * t29 - t75 * t91;
t31 = t79 * t41 + t75 * t96;
t131 = t31 * t18;
t130 = t31 * t78;
t19 = t31 * qJD(5) + t75 * t29 + t79 * t91;
t129 = t74 * t19;
t127 = t78 * t18;
t126 = t78 * t19;
t119 = qJD(5) * t75;
t66 = t80 * pkin(3) + pkin(4);
t40 = t66 * t119 + t137;
t116 = t75 * t76 * pkin(3);
t48 = -t79 * t66 - pkin(5) + t116;
t124 = t40 * t74 + t48 * t67;
t111 = pkin(4) * t119;
t65 = -t79 * pkin(4) - pkin(5);
t123 = t74 * t111 + t65 * t67;
t117 = qJD(6) * t74;
t115 = pkin(5) * t117;
t114 = pkin(5) * t67;
t113 = pkin(3) * t121;
t112 = pkin(3) * t120;
t110 = pkin(4) * t118;
t108 = t74 * t67;
t107 = -0.4e1 * t74 * t130;
t42 = t48 * t117;
t106 = -t40 * t78 + t42;
t104 = 0.2e1 * (t72 ^ 2 + t73 ^ 2) * qJD(2);
t44 = pkin(3) * t93 + t63;
t34 = -t96 * pkin(4) + t44;
t22 = t30 * pkin(5) - t31 * pkin(10) + t34;
t103 = t78 * t17 + t74 * t22;
t102 = t74 * t17 - t78 * t22;
t49 = pkin(3) * t128 + t75 * t66 + pkin(10);
t101 = t30 * t49 - t31 * t48;
t64 = t75 * pkin(4) + pkin(10);
t100 = t30 * t64 - t31 * t65;
t55 = t65 * t117;
t90 = -t78 * t111 + t55;
t89 = t74 * t18 + t31 * t67;
t88 = t31 * t117 - t127;
t10 = t30 * t67 + t129;
t87 = t30 * t117 - t126;
t39 = -t66 * t118 - t79 * t112 + (qJD(4) + qJD(5)) * t116;
t26 = t91 * pkin(4) + t133;
t86 = t18 * t48 - t19 * t49 + t30 * t39 + t31 * t40;
t85 = t18 * t65 - t19 * t64 + (-t30 * t79 + t31 * t75) * qJD(5) * pkin(4);
t62 = 0.2e1 * t108;
t54 = -0.2e1 * t105;
t28 = t31 ^ 2;
t14 = t16 * t117;
t8 = -t31 * t105 + t74 * t127;
t7 = qJD(6) * t107 - t122 * t18;
t6 = t19 * pkin(5) - t18 * pkin(10) + t26;
t4 = -t24 * t118 + t25 * t119 + t75 * t84 - t79 * t82;
t2 = -t103 * qJD(6) + t74 * t4 + t78 * t6;
t1 = t102 * qJD(6) + t78 * t4 - t74 * t6;
t9 = [0, 0, 0, 0, 0, t104, qJ(2) * t104, -0.2e1 * t53 * t50, 0.2e1 * t50 * t93 - 0.2e1 * t53 * t51, 0, 0, 0, t51 * t134, -t50 * t134, 0.2e1 * t41 * t29, 0.2e1 * t29 * t96 - 0.2e1 * t41 * t91, 0, 0, 0, -0.2e1 * t96 * t133 + 0.2e1 * t44 * t91, 0.2e1 * t41 * t133 + 0.2e1 * t44 * t29, 0.2e1 * t131, -0.2e1 * t18 * t30 - 0.2e1 * t31 * t19, 0, 0, 0, 0.2e1 * t34 * t19 + 0.2e1 * t26 * t30, 0.2e1 * t34 * t18 + 0.2e1 * t26 * t31, -0.2e1 * t28 * t108 + 0.2e1 * t71 * t131, 0.2e1 * t28 * t105 + t18 * t107, 0.2e1 * t31 * t126 - 0.2e1 * t88 * t30, -0.2e1 * t31 * t129 - 0.2e1 * t89 * t30, 0.2e1 * t30 * t19, -0.2e1 * t102 * t19 + 0.2e1 * t89 * t16 + 0.2e1 * t2 * t30 + 0.2e1 * t31 * t3, 0.2e1 * t1 * t30 - 0.2e1 * t103 * t19 + 0.2e1 * t5 * t130 - 0.2e1 * t88 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, 0, 0, 0, 0, 0, t91, t29, 0, 0, 0, 0, 0, t19, t18, 0, 0, 0, 0, 0, -t87, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t51, 0, t83, t135, 0, 0, t29, -t91, 0, t21, t20, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t87, 0, t14 + (-t101 * qJD(6) - t5) * t78 + t86 * t74, t101 * t117 + t86 * t78 + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113, -0.2e1 * t112, 0, 0, 0, 0, 0, -0.2e1 * t40, 0.2e1 * t39, t62, t54, 0, 0, 0, 0.2e1 * t106, 0.2e1 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t91, 0, t21, t20, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t87, 0, t14 + (-t100 * qJD(6) - t5) * t78 + t85 * t74, t100 * t117 + t85 * t78 + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, 0, 0, 0, 0 (-pkin(4) - t66) * t119 - t137, t39 - t110, t62, t54, 0, 0, 0, t42 + t55 + (-t40 - t111) * t78, t123 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t111, -0.2e1 * t110, t62, t54, 0, 0, 0, 0.2e1 * t90, 0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t87, 0, t14 + (-pkin(5) * t18 - pkin(10) * t19) * t74 + (-t5 + (-pkin(5) * t31 - pkin(10) * t30) * qJD(6)) * t78, t88 * pkin(5) + t87 * pkin(10) + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t39, t62, t54, 0, 0, 0, t106 - t115, -t114 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, t62, t54, 0, 0, 0, t90 - t115, -t114 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t54, 0, 0, 0, -0.2e1 * t115, -0.2e1 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t89, t19, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t117, 0, t74 * t39 - t49 * t67, t49 * t117 + t78 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t117, 0, -t74 * t110 - t64 * t67, -t78 * t110 + t64 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t117, 0, -pkin(10) * t67, pkin(10) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
