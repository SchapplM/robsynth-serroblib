% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:12
% EndTime: 2019-03-08 20:34:15
% DurationCPUTime: 0.96s
% Computational Cost: add. (1761->144), mult. (4366->266), div. (0->0), fcn. (4726->12), ass. (0->108)
t74 = cos(qJ(6));
t65 = t74 ^ 2;
t70 = sin(qJ(6));
t111 = t70 ^ 2 - t65;
t96 = t111 * qJD(6);
t66 = sin(pkin(12));
t68 = cos(pkin(12));
t72 = sin(qJ(4));
t76 = cos(qJ(4));
t49 = t72 * t66 - t76 * t68;
t114 = pkin(8) + qJ(3);
t54 = t114 * t66;
t55 = t114 * t68;
t90 = -t76 * t54 - t72 * t55;
t126 = t49 * qJD(3) - t90 * qJD(4);
t58 = -t68 * pkin(3) - pkin(2);
t125 = 0.2e1 * t58;
t50 = t76 * t66 + t72 * t68;
t48 = t50 * qJD(4);
t124 = t48 * pkin(4);
t35 = -t50 * pkin(9) + t90;
t89 = t72 * t54 - t76 * t55;
t36 = -t49 * pkin(9) - t89;
t71 = sin(qJ(5));
t75 = cos(qJ(5));
t20 = t71 * t35 + t75 * t36;
t78 = -t48 * pkin(9) - t126;
t47 = t49 * qJD(4);
t81 = -t50 * qJD(3) + t89 * qJD(4);
t79 = -t47 * pkin(9) - t81;
t11 = t20 * qJD(5) + t71 * t78 + t75 * t79;
t19 = -t75 * t35 + t71 * t36;
t61 = qJD(6) * t74;
t123 = t11 * t70 + t19 * t61;
t40 = t75 * t49 + t71 * t50;
t24 = -t40 * qJD(5) - t75 * t47 - t71 * t48;
t41 = -t71 * t49 + t75 * t50;
t122 = t41 * t24;
t121 = t41 * t70;
t120 = t41 * t74;
t67 = sin(pkin(6));
t73 = sin(qJ(2));
t119 = t67 * t73;
t77 = cos(qJ(2));
t118 = t67 * t77;
t25 = t41 * qJD(5) - t71 * t47 + t75 * t48;
t117 = t70 * t25;
t116 = t74 * t24;
t115 = t74 * t25;
t108 = qJD(5) * t71;
t103 = pkin(4) * t108;
t60 = -t75 * pkin(4) - pkin(5);
t113 = t70 * t103 + t60 * t61;
t112 = t66 ^ 2 + t68 ^ 2;
t110 = qJD(2) * t67;
t109 = qJD(2) * t73;
t107 = qJD(5) * t75;
t106 = qJD(6) * t70;
t105 = pkin(5) * t106;
t104 = pkin(5) * t61;
t102 = pkin(4) * t107;
t101 = t67 * t109;
t100 = t77 * t110;
t99 = t70 * t61;
t98 = t112 * t77;
t97 = -0.4e1 * t70 * t120;
t95 = 0.2e1 * t112 * qJD(3);
t42 = t49 * pkin(4) + t58;
t21 = t40 * pkin(5) - t41 * pkin(10) + t42;
t94 = t74 * t20 + t70 * t21;
t93 = t70 * t20 - t74 * t21;
t69 = cos(pkin(6));
t44 = -t66 * t119 + t69 * t68;
t45 = t68 * t119 + t69 * t66;
t37 = t76 * t44 - t72 * t45;
t38 = t72 * t44 + t76 * t45;
t23 = t71 * t37 + t75 * t38;
t59 = t71 * pkin(4) + pkin(10);
t92 = t40 * t59 - t41 * t60;
t91 = -t44 * t66 + t45 * t68;
t88 = -t74 * t103 + t60 * t106;
t87 = t74 * t118 + t70 * t23;
t86 = t70 * t118 - t74 * t23;
t85 = t70 * t24 + t41 * t61;
t84 = t41 * t106 - t116;
t16 = t40 * t61 + t117;
t83 = t40 * t106 - t115;
t82 = t24 * t60 - t25 * t59 + (-t40 * t75 + t41 * t71) * qJD(5) * pkin(4);
t80 = -t38 * qJD(4) - t50 * t100;
t57 = 0.2e1 * t99;
t51 = -0.2e1 * t96;
t39 = t41 ^ 2;
t29 = -t37 * qJD(4) + t49 * t100;
t22 = -t75 * t37 + t71 * t38;
t17 = t19 * t106;
t14 = t70 * t116 - t41 * t96;
t13 = t25 * pkin(5) - t24 * pkin(10) + t124;
t12 = qJD(6) * t97 - t111 * t24;
t10 = -t35 * t107 + t36 * t108 + t71 * t79 - t75 * t78;
t8 = t23 * qJD(5) - t71 * t29 - t75 * t80;
t7 = -t37 * t107 + t38 * t108 + t75 * t29 - t71 * t80;
t6 = t22 * t106 - t8 * t74;
t5 = t22 * t61 + t8 * t70;
t4 = t86 * qJD(6) + t74 * t101 + t70 * t7;
t3 = t87 * qJD(6) - t70 * t101 + t74 * t7;
t2 = -t94 * qJD(6) + t70 * t10 + t74 * t13;
t1 = t93 * qJD(6) + t74 * t10 - t70 * t13;
t9 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t91 - t119) * t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t101, -t100, -t68 * t101, t66 * t101, t98 * t110, t91 * qJD(3) + (-pkin(2) * t73 + qJ(3) * t98) * t110, 0, 0, 0, 0, 0 (t49 * t109 - t48 * t77) * t67 (t50 * t109 + t47 * t77) * t67, 0, 0, 0, 0, 0 (t40 * t109 - t25 * t77) * t67 (t41 * t109 - t24 * t77) * t67, 0, 0, 0, 0, 0, t8 * t121 + t85 * t22 - t87 * t25 + t4 * t40, t8 * t120 - t84 * t22 + t86 * t25 + t3 * t40; 0, 0, 0, 0, 0, 0, t95, qJ(3) * t95, -0.2e1 * t50 * t47, 0.2e1 * t47 * t49 - 0.2e1 * t50 * t48, 0, 0, 0, t48 * t125, -t47 * t125, 0.2e1 * t122, -0.2e1 * t24 * t40 - 0.2e1 * t41 * t25, 0, 0, 0, 0.2e1 * t40 * t124 + 0.2e1 * t42 * t25, 0.2e1 * t41 * t124 + 0.2e1 * t42 * t24, 0.2e1 * t65 * t122 - 0.2e1 * t39 * t99, t24 * t97 + 0.2e1 * t39 * t96, 0.2e1 * t41 * t115 - 0.2e1 * t84 * t40, -0.2e1 * t41 * t117 - 0.2e1 * t85 * t40, 0.2e1 * t40 * t25, 0.2e1 * t11 * t121 + 0.2e1 * t85 * t19 + 0.2e1 * t2 * t40 - 0.2e1 * t93 * t25, 0.2e1 * t1 * t40 + 0.2e1 * t11 * t120 - 0.2e1 * t84 * t19 - 0.2e1 * t94 * t25; 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, 0, 0, 0, 0, 0, t25, t24, 0, 0, 0, 0, 0, -t83, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t29, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, 0, t81, t126, 0, 0, t24, -t25, 0, -t11, t10, t14, t12, t16, -t83, 0, t17 + (-t92 * qJD(6) - t11) * t74 + t82 * t70, t92 * t106 + t82 * t74 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t103, -0.2e1 * t102, t57, t51, 0, 0, 0, 0.2e1 * t88, 0.2e1 * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25, 0, -t11, t10, t14, t12, t16, -t83, 0, t17 + (-pkin(5) * t24 - pkin(10) * t25) * t70 + (-t11 + (-pkin(5) * t41 - pkin(10) * t40) * qJD(6)) * t74, t84 * pkin(5) + t83 * pkin(10) + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t102, t57, t51, 0, 0, 0, t88 - t105, -t104 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t51, 0, 0, 0, -0.2e1 * t105, -0.2e1 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t85, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t106, 0, -t70 * t102 - t59 * t61, -t74 * t102 + t59 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t106, 0, -pkin(10) * t61, pkin(10) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
