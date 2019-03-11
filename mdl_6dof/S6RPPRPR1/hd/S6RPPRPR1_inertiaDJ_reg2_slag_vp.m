% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:03
% EndTime: 2019-03-09 01:40:08
% DurationCPUTime: 1.80s
% Computational Cost: add. (2983->183), mult. (6530->327), div. (0->0), fcn. (6566->10), ass. (0->115)
t158 = cos(qJ(4));
t124 = qJD(4) * t158;
t87 = sin(qJ(4));
t133 = qJD(4) * t87;
t84 = sin(pkin(10));
t86 = cos(pkin(10));
t62 = -t86 * t124 + t84 * t133;
t69 = t158 * t84 + t87 * t86;
t63 = t69 * qJD(4);
t67 = -t158 * t86 + t87 * t84;
t163 = t67 * t62 - t69 * t63;
t165 = 0.2e1 * t163;
t75 = sin(pkin(9)) * pkin(1) + qJ(3);
t159 = pkin(7) + t75;
t113 = t159 * t158;
t120 = t158 * qJD(3);
t132 = t87 * qJD(3);
t64 = t159 * t86;
t162 = (-qJD(4) * t113 - t132) * t84 + t86 * t120 - t64 * t133;
t156 = sin(qJ(6));
t121 = qJD(6) * t156;
t157 = cos(qJ(6));
t122 = qJD(6) * t157;
t83 = sin(pkin(11));
t85 = cos(pkin(11));
t60 = t83 * t121 - t85 * t122;
t128 = t157 * t83;
t68 = t156 * t85 + t128;
t20 = -t60 * t69 - t68 * t62;
t39 = t68 * t69;
t138 = -t68 * t20 + t39 * t60;
t61 = t68 * qJD(6);
t66 = t156 * t83 - t157 * t85;
t19 = t69 * t61 - t66 * t62;
t40 = t66 * t69;
t143 = t66 * t19 + t40 * t61;
t164 = t143 - t138;
t161 = t85 * t162;
t48 = t67 * t63;
t45 = 0.2e1 * t48;
t160 = t63 * pkin(4);
t155 = t39 * t20;
t154 = t40 * t19;
t129 = t87 * t159;
t32 = t64 * t124 + t86 * t132 + (-qJD(4) * t129 + t120) * t84;
t41 = t84 * t113 + t87 * t64;
t153 = t41 * t32;
t152 = t66 * t61;
t151 = t68 * t60;
t150 = t69 * t62;
t79 = t83 ^ 2;
t149 = t79 * t62;
t42 = -t84 * t129 + t158 * t64;
t148 = t83 * t42;
t147 = t83 * t62;
t55 = t83 * t63;
t146 = t83 * t69;
t145 = t85 * t62;
t57 = t85 * t63;
t144 = pkin(8) + qJ(5);
t141 = -t19 * t68 + t40 * t60;
t140 = t20 * t66 + t39 * t61;
t70 = -cos(pkin(9)) * pkin(1) - pkin(2) - t86 * pkin(3);
t96 = t67 * pkin(4) + t70;
t94 = -t69 * qJ(5) + t96;
t22 = t85 * t42 + t83 * t94;
t137 = -t67 * t145 + t69 * t57;
t81 = t85 ^ 2;
t134 = t79 + t81;
t131 = -0.2e1 * t150;
t130 = t83 * t145;
t125 = t134 * t62;
t123 = qJD(5) * t156;
t119 = t157 * qJD(5);
t118 = t134 * qJD(5);
t117 = 0.2e1 * t118;
t116 = 0.2e1 * (t84 ^ 2 + t86 ^ 2) * qJD(3);
t112 = t144 * t156;
t103 = -t69 * qJD(5) + t160;
t100 = t62 * qJ(5) + t103;
t90 = t162 * t83;
t10 = t85 * t100 - t90;
t11 = t83 * t100 + t161;
t111 = t10 * t85 + t11 * t83;
t110 = -t10 * t83 + t11 * t85;
t21 = t85 * t94 - t148;
t109 = t21 * t83 - t22 * t85;
t108 = t32 * t69 - t41 * t62;
t107 = t32 * t67 + t41 * t63;
t106 = t60 * t67 - t68 * t63;
t105 = t67 * t61 + t63 * t66;
t102 = t144 * t128;
t101 = pkin(4) * t62 - qJ(5) * t63 - qJD(5) * t67;
t98 = t144 * t62 + t103;
t93 = t67 * pkin(5) - t148 + (-t144 * t69 + t96) * t85;
t92 = t157 * t93;
t91 = t156 * t93;
t89 = t98 * t83 + t161;
t88 = t63 * pkin(5) + t98 * t85 - t90;
t78 = -t85 * pkin(5) - pkin(4);
t71 = t144 * t85;
t56 = t81 * t62;
t52 = -t83 * t112 + t157 * t71;
t51 = -t156 * t71 - t102;
t38 = t56 + t149;
t35 = -t71 * t122 - t85 * t123 + (qJD(6) * t112 - t119) * t83;
t34 = qJD(6) * t102 - t85 * t119 + t71 * t121 + t83 * t123;
t33 = pkin(5) * t146 + t41;
t23 = -pkin(5) * t147 + t32;
t12 = -pkin(8) * t146 + t22;
t4 = t157 * t12 + t91;
t3 = -t156 * t12 + t92;
t2 = -qJD(6) * t91 - t12 * t122 - t156 * t89 + t157 * t88;
t1 = -qJD(6) * t92 + t12 * t121 - t156 * t88 - t157 * t89;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t75 * t116, t131, t165, 0, t45, 0, 0, 0.2e1 * t70 * t63, -0.2e1 * t70 * t62, -0.2e1 * t162 * t67 - 0.2e1 * t42 * t63 + 0.2e1 * t108, 0.2e1 * t162 * t42 + 0.2e1 * t153, t81 * t131, 0.4e1 * t69 * t130, 0.2e1 * t137, t79 * t131, t83 * t165, t45, 0.2e1 * t10 * t67 + 0.2e1 * t108 * t83 + 0.2e1 * t21 * t63, 0.2e1 * t108 * t85 - 0.2e1 * t11 * t67 - 0.2e1 * t22 * t63, -0.2e1 * t111 * t69 + 0.2e1 * (t21 * t85 + t22 * t83) * t62, 0.2e1 * t21 * t10 + 0.2e1 * t22 * t11 + 0.2e1 * t153, 0.2e1 * t154, 0.2e1 * t39 * t19 + 0.2e1 * t40 * t20, -0.2e1 * t19 * t67 - 0.2e1 * t40 * t63, 0.2e1 * t155, -0.2e1 * t67 * t20 - 0.2e1 * t63 * t39, t45, 0.2e1 * t2 * t67 + 0.2e1 * t33 * t20 + 0.2e1 * t23 * t39 + 0.2e1 * t3 * t63, 0.2e1 * t1 * t67 - 0.2e1 * t33 * t19 - 0.2e1 * t23 * t40 - 0.2e1 * t4 * t63, 0.2e1 * t1 * t39 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t40 - 0.2e1 * t4 * t20, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t33 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162 * t69 - t42 * t62 + t107, 0, 0, 0, 0, 0, 0, 0, t163 * t85 + t137, 0, t109 * t62 + t110 * t69 + t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t40 - t4 * t19 - t2 * t39 - t3 * t20 + t23 * t67 + t33 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t150 + 0.2e1 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t69 * t125 + 0.2e1 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t154 + 0.2e1 * t48 + 0.2e1 * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t62, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t55, t38, t111, 0, 0, 0, 0, 0, 0, -t105, t106, -t164, -t1 * t68 - t2 * t66 - t3 * t61 - t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t151 + 0.2e1 * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, -t63, 0, -t32, -t162, 0, 0, -t130, -t56 + t149, t55, t130, t57, 0, t101 * t83 - t32 * t85, t101 * t85 + t32 * t83, t110, -t32 * pkin(4) + t110 * qJ(5) - t109 * qJD(5), t141, t138 + t143, -t106, t140, -t105, 0, t78 * t20 + t23 * t66 + t33 * t61 + t35 * t67 + t51 * t63, -t78 * t19 + t23 * t68 - t33 * t60 + t34 * t67 - t52 * t63, t1 * t66 + t51 * t19 - t2 * t68 - t52 * t20 + t3 * t60 + t34 * t39 + t35 * t40 - t4 * t61, -t1 * t52 + t2 * t51 + t23 * t78 + t3 * t35 - t4 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t62, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t55, -t38, -qJ(5) * t125 + t118 * t69 - t160, 0, 0, 0, 0, 0, 0, t105, -t106, t164, -t19 * t52 - t20 * t51 + t40 * t34 - t39 * t35 + t63 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 * t34 - t66 * t35 - t61 * t51 - t60 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, qJ(5) * t117, -0.2e1 * t151, 0.2e1 * t66 * t60 - 0.2e1 * t68 * t61, 0, 0.2e1 * t152, 0, 0, 0.2e1 * t78 * t61, -0.2e1 * t78 * t60, 0.2e1 * t34 * t66 - 0.2e1 * t35 * t68 + 0.2e1 * t51 * t60 - 0.2e1 * t52 * t61, -0.2e1 * t52 * t34 + 0.2e1 * t51 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t145, 0, t32, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, -t20, t63, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, -t61, 0, t35, t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
