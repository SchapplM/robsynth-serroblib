% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:46
% EndTime: 2019-03-09 03:22:52
% DurationCPUTime: 1.98s
% Computational Cost: add. (2401->200), mult. (4789->337), div. (0->0), fcn. (4361->6), ass. (0->123)
t124 = sin(pkin(9));
t125 = cos(pkin(9));
t143 = -pkin(1) - pkin(7);
t117 = qJ(4) - t143;
t73 = sin(qJ(3));
t121 = t73 * qJD(3);
t75 = cos(qJ(3));
t83 = -t75 * qJD(4) + t117 * t121;
t101 = t117 * t75;
t84 = -qJD(3) * t101 - t73 * qJD(4);
t78 = t124 * t83 + t125 * t84;
t54 = t124 * t75 + t125 * t73;
t55 = -t124 * t73 + t125 * t75;
t66 = t73 * pkin(3) + qJ(2);
t86 = t54 * pkin(4) - t55 * pkin(8) + t66;
t149 = -qJD(5) * t86 - t78;
t48 = t55 * qJD(3);
t49 = t54 * qJD(3);
t148 = (-t124 * t48 + t125 * t49) * pkin(3);
t72 = sin(qJ(5));
t70 = t72 ^ 2;
t74 = cos(qJ(5));
t71 = t74 ^ 2;
t129 = t70 - t71;
t138 = t49 * t72;
t68 = qJD(5) * t74;
t34 = t55 * t68 - t138;
t119 = t75 * qJD(3);
t61 = pkin(3) * t119 + qJD(2);
t82 = t48 * pkin(4) + t49 * pkin(8) + t61;
t114 = t149 * t74 - t72 * t82;
t122 = qJD(5) * t72;
t57 = t117 * t73;
t40 = -t124 * t101 - t125 * t57;
t4 = t40 * t122 + t114;
t5 = t149 * t72 - t40 * t68 + t74 * t82;
t14 = -t72 * t40 + t74 * t86;
t15 = t74 * t40 + t72 * t86;
t94 = t14 * t72 - t15 * t74;
t147 = qJD(5) * t94 + t4 * t72 - t5 * t74;
t142 = t48 * pkin(5);
t123 = qJD(5) * t55;
t107 = qJ(6) * t123;
t120 = t74 * qJD(6);
t127 = t74 * qJ(6);
t76 = t72 * t107 - t55 * t120 + t49 * t127 + t5;
t1 = t76 + t142;
t2 = t74 * t107 + (-t49 * qJ(6) + qJD(5) * t40 + t55 * qJD(6)) * t72 + t114;
t134 = t55 * t72;
t11 = -qJ(6) * t134 + t15;
t7 = t54 * pkin(5) - t55 * t127 + t14;
t97 = t11 * t74 - t7 * t72;
t146 = qJD(5) * t97 + t1 * t74 - t2 * t72;
t53 = t55 ^ 2;
t145 = 0.2e1 * qJD(2);
t144 = 0.2e1 * qJD(5);
t141 = t74 * pkin(5);
t24 = t124 * t84 - t125 * t83;
t39 = t125 * t101 - t124 * t57;
t140 = t39 * t24;
t65 = -t125 * pkin(3) - pkin(4);
t139 = t49 * t65;
t137 = t49 * t74;
t136 = t54 * t48;
t135 = t55 * t49;
t133 = t55 * t74;
t58 = t65 - t141;
t132 = t58 * t74;
t131 = t72 * t48;
t130 = t74 * t48;
t128 = t70 + t71;
t64 = t124 * pkin(3) + pkin(8);
t126 = qJ(6) + t64;
t118 = qJ(2) * qJD(3);
t38 = 0.2e1 * t136;
t116 = t72 * t137;
t115 = t65 * t144;
t113 = pkin(5) * t122;
t112 = t55 * t122;
t110 = t72 * t68;
t109 = t73 * t119;
t108 = t54 ^ 2 + t53;
t106 = -t58 + t141;
t105 = t143 * qJD(3);
t27 = t128 * t48;
t102 = t53 * t110;
t98 = pkin(5) * t70 + t132;
t96 = t11 * t72 + t7 * t74;
t95 = t14 * t74 + t15 * t72;
t93 = -t24 * t55 + t39 * t49;
t92 = t135 - t136;
t91 = -t48 * t64 - t139;
t50 = t126 * t72;
t51 = t126 * t74;
t90 = -t50 * t72 - t51 * t74;
t89 = t54 * t64 - t55 * t65;
t32 = t54 * t68 + t131;
t88 = t112 + t137;
t87 = 0.2e1 * t92;
t80 = -t95 * qJD(5) - t4 * t74 - t5 * t72;
t41 = t126 * t122 - t120;
t42 = -t72 * qJD(6) - t126 * t68;
t79 = -t41 * t74 - t42 * t72 + (t50 * t74 - t51 * t72) * qJD(5);
t77 = t40 * t48 + t78 * t54 + t93;
t69 = qJ(2) * t145;
t60 = -0.2e1 * t110;
t59 = 0.2e1 * t110;
t56 = -0.2e1 * t129 * qJD(5);
t30 = t54 * t122 - t130;
t26 = t128 * t49;
t21 = pkin(5) * t134 + t39;
t19 = -0.2e1 * t71 * t135 - 0.2e1 * t102;
t18 = -0.2e1 * t70 * t135 + 0.2e1 * t102;
t17 = t129 * t123 + t116;
t16 = -0.4e1 * t55 * t110 + t129 * t49;
t13 = t129 * t53 * t144 + 0.4e1 * t55 * t116;
t12 = t34 * pkin(5) + t24;
t10 = 0.2e1 * t54 * t27 - 0.2e1 * t135;
t9 = -0.2e1 * t55 * t131 - 0.2e1 * t34 * t54;
t8 = 0.2e1 * t55 * t130 - 0.2e1 * t88 * t54;
t6 = -t108 * t68 + t72 * t87;
t3 = t108 * t122 + t74 * t87;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t69, -0.2e1 * t109, 0.2e1 * (t73 ^ 2 - t75 ^ 2) * qJD(3), 0, 0.2e1 * t109, 0, 0, 0.2e1 * qJD(2) * t73 + 0.2e1 * t75 * t118, 0.2e1 * qJD(2) * t75 - 0.2e1 * t73 * t118, 0, t69, -0.2e1 * t135, -0.2e1 * t55 * t48 + 0.2e1 * t49 * t54, 0, t38, 0, 0, 0.2e1 * t66 * t48 + 0.2e1 * t61 * t54, -0.2e1 * t66 * t49 + 0.2e1 * t61 * t55, -0.2e1 * t77, 0.2e1 * t40 * t78 + 0.2e1 * t66 * t61 + 0.2e1 * t140, t19, t13, t8, t18, t9, t38, 0.2e1 * t134 * t24 + 0.2e1 * t14 * t48 + 0.2e1 * t34 * t39 + 0.2e1 * t5 * t54, 0.2e1 * t133 * t24 - 0.2e1 * t15 * t48 - 0.2e1 * t39 * t88 + 0.2e1 * t4 * t54, 0.2e1 * t147 * t55 + 0.2e1 * t95 * t49, 0.2e1 * t14 * t5 - 0.2e1 * t15 * t4 + 0.2e1 * t140, t19, t13, t8, t18, t9, t38, 0.2e1 * t1 * t54 + 0.2e1 * t12 * t134 + 0.2e1 * t21 * t34 + 0.2e1 * t7 * t48, -0.2e1 * t11 * t48 + 0.2e1 * t12 * t133 + 0.2e1 * t2 * t54 - 0.2e1 * t21 * t88, -0.2e1 * t146 * t55 + 0.2e1 * t96 * t49, 0.2e1 * t1 * t7 - 0.2e1 * t11 * t2 + 0.2e1 * t12 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t77, 0, 0, 0, 0, 0, 0, t6, t3, 0, -t48 * t94 + t54 * t80 + t93, 0, 0, 0, 0, 0, 0, t6, t3, 0, -t12 * t55 + t21 * t49 + t97 * t48 + (-qJD(5) * t96 - t1 * t72 - t2 * t74) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, 0, -t119, 0, -t73 * t105, -t75 * t105, 0, 0, 0, 0, -t49, 0, -t48, 0, -t24, -t78, t148 (t124 * t78 - t125 * t24) * pkin(3), -t17, t16, t32, t17, -t30, 0, -t24 * t74 + t91 * t72 + (t39 * t72 - t74 * t89) * qJD(5), t24 * t72 + t91 * t74 + (t39 * t74 + t72 * t89) * qJD(5), t80, t24 * t65 + t64 * t80, -t17, t16, t32, t17, -t30, 0, -t58 * t138 - t12 * t74 + t42 * t54 - t50 * t48 + (t21 * t72 + t55 * t98) * qJD(5), -t49 * t132 + t12 * t72 + t41 * t54 - t51 * t48 + (t106 * t134 + t21 * t74) * qJD(5) (-t42 * t55 - t49 * t50 - t2 + (-t51 * t55 - t7) * qJD(5)) * t74 + (t41 * t55 + t49 * t51 - t1 + (-t50 * t55 - t11) * qJD(5)) * t72, -t1 * t50 - t11 * t41 + t113 * t21 + t12 * t58 - t2 * t51 + t7 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t119, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, 0, -t148, 0, 0, 0, 0, 0, 0, -t88, -t34, t27, t27 * t64 + t139, 0, 0, 0, 0, 0, 0, -t88, -t34, t27, -pkin(5) * t112 - t48 * t90 + t49 * t58 + t54 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t56, 0, t60, 0, 0, t72 * t115, t74 * t115, 0, 0, t59, t56, 0, t60, 0, 0, -0.2e1 * t106 * t122, t98 * t144, 0.2e1 * t79, 0.2e1 * t113 * t58 - 0.2e1 * t51 * t41 - 0.2e1 * t50 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, 0, t61, 0, 0, 0, 0, 0, 0, -t30, -t32, t26, -t147, 0, 0, 0, 0, 0, 0, -t30, -t32, t26, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t90 - t41 * t72 + t42 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, 0, -t34, t48, t5, t4, 0, 0, 0, 0, -t88, 0, -t34, t48, t76 + 0.2e1 * t142, t2, t88 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, -t32 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, -t122, 0, -t64 * t68, t64 * t122, 0, 0, 0, 0, t68, 0, -t122, 0, t42, t41, -pkin(5) * t68, t42 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t68, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t68, 0, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t88, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t68, 0, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;