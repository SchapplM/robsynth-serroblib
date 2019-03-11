% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:19
% EndTime: 2019-03-09 03:16:23
% DurationCPUTime: 1.28s
% Computational Cost: add. (2962->185), mult. (6779->333), div. (0->0), fcn. (6789->8), ass. (0->93)
t89 = cos(pkin(10));
t81 = -t89 * pkin(4) - pkin(3);
t141 = 0.2e1 * t81;
t133 = cos(qJ(5));
t90 = cos(pkin(9));
t92 = sin(qJ(3));
t125 = t92 * t90;
t134 = cos(qJ(3));
t88 = sin(pkin(9));
t70 = t134 * t88 + t125;
t82 = -t90 * pkin(2) - pkin(1);
t95 = t134 * t90 - t92 * t88;
t45 = -pkin(3) * t95 - t70 * qJ(4) + t82;
t124 = pkin(7) + qJ(2);
t71 = t124 * t88;
t126 = t92 * t71;
t73 = t124 * t90;
t52 = t134 * t73 - t126;
t87 = sin(pkin(10));
t24 = t89 * t45 - t87 * t52;
t18 = -t89 * t70 * pkin(8) - pkin(4) * t95 + t24;
t131 = t70 * t87;
t25 = t87 * t45 + t89 * t52;
t21 = -pkin(8) * t131 + t25;
t91 = sin(qJ(5));
t140 = t133 * t21 + t91 * t18;
t109 = qJD(5) * t133;
t117 = qJD(5) * t91;
t139 = t89 * t109 - t87 * t117;
t138 = t133 * t89 - t91 * t87;
t60 = t95 * qJD(3);
t129 = t89 * t60;
t61 = t70 * qJD(3);
t33 = t61 * pkin(3) - t60 * qJ(4) - t70 * qJD(4);
t110 = qJD(3) * t134;
t111 = qJD(2) * t134;
t37 = t71 * t110 - t90 * t111 + (qJD(2) * t88 + qJD(3) * t73) * t92;
t16 = t89 * t33 + t87 * t37;
t10 = t61 * pkin(4) - pkin(8) * t129 + t16;
t130 = t87 * t60;
t17 = t87 * t33 - t89 * t37;
t13 = -pkin(8) * t130 + t17;
t4 = -qJD(5) * t140 + t133 * t10 - t91 * t13;
t137 = 0.2e1 * t82;
t136 = 2 * qJD(6);
t135 = t61 * pkin(5);
t127 = t91 * t89;
t69 = t133 * t87 + t127;
t132 = t69 * t139;
t123 = pkin(8) + qJ(4);
t23 = t139 * t70 + t69 * t60;
t40 = t69 * t70;
t121 = -t139 * t40 - t69 * t23;
t119 = t87 ^ 2 + t89 ^ 2;
t118 = t61 * qJ(6);
t116 = t95 * qJD(6);
t112 = t123 * t87;
t108 = t133 * qJD(4);
t38 = qJD(2) * t125 - qJD(3) * t126 + t73 * t110 + t88 * t111;
t107 = 0.2e1 * t119 * qJD(4);
t106 = 0.2e1 * (t88 ^ 2 + t90 ^ 2) * qJD(2);
t26 = pkin(4) * t130 + t38;
t50 = t134 * t71 + t92 * t73;
t104 = t16 * t89 + t17 * t87;
t103 = -t16 * t87 + t17 * t89;
t59 = t69 * qJD(5);
t22 = -t138 * t60 + t70 * t59;
t41 = t138 * t70;
t102 = t138 * t22 + t59 * t41;
t72 = t123 * t89;
t97 = t133 * t112;
t35 = qJD(5) * t97 - t89 * t108 + (qJD(4) * t87 + qJD(5) * t72) * t91;
t51 = -t91 * t112 + t133 * t72;
t101 = -t35 * t95 - t51 * t61;
t36 = t72 * t109 + qJD(4) * t127 + (-t123 * t117 + t108) * t87;
t49 = t91 * t72 + t97;
t100 = t36 * t95 - t49 * t61;
t99 = t38 * t70 + t50 * t60;
t98 = t139 * t95 - t69 * t61;
t28 = t138 * t61 + t59 * t95;
t39 = pkin(4) * t131 + t50;
t96 = t133 * t18 - t91 * t21;
t3 = -t91 * t10 - t18 * t109 + t21 * t117 - t133 * t13;
t93 = -pkin(3) * t60 - qJ(4) * t61 + qJD(4) * t95;
t32 = t59 * pkin(5) - qJ(6) * t139 - t69 * qJD(6);
t44 = -pkin(5) * t138 - t69 * qJ(6) + t81;
t12 = t40 * pkin(5) - t41 * qJ(6) + t39;
t7 = pkin(5) * t95 - t96;
t6 = -qJ(6) * t95 + t140;
t5 = t23 * pkin(5) + t22 * qJ(6) - t41 * qJD(6) + t26;
t2 = -t135 - t4;
t1 = -t116 - t3 + t118;
t8 = [0, 0, 0, 0, 0, t106, qJ(2) * t106, 0.2e1 * t70 * t60, 0.2e1 * t60 * t95 - 0.2e1 * t70 * t61, 0, 0, 0, t61 * t137, t60 * t137, -0.2e1 * t16 * t95 + 0.2e1 * t24 * t61 + 0.2e1 * t99 * t87, 0.2e1 * t17 * t95 - 0.2e1 * t25 * t61 + 0.2e1 * t99 * t89, -0.2e1 * t104 * t70 + 0.2e1 * (-t24 * t89 - t25 * t87) * t60, 0.2e1 * t24 * t16 + 0.2e1 * t25 * t17 + 0.2e1 * t50 * t38, -0.2e1 * t41 * t22, 0.2e1 * t22 * t40 - 0.2e1 * t41 * t23, 0.2e1 * t22 * t95 + 0.2e1 * t41 * t61, 0.2e1 * t23 * t95 - 0.2e1 * t40 * t61, -0.2e1 * t95 * t61, 0.2e1 * t39 * t23 + 0.2e1 * t26 * t40 - 0.2e1 * t4 * t95 + 0.2e1 * t96 * t61, -0.2e1 * t140 * t61 - 0.2e1 * t39 * t22 + 0.2e1 * t26 * t41 - 0.2e1 * t3 * t95, 0.2e1 * t12 * t23 + 0.2e1 * t2 * t95 + 0.2e1 * t5 * t40 - 0.2e1 * t7 * t61, -0.2e1 * t1 * t40 + 0.2e1 * t2 * t41 - 0.2e1 * t7 * t22 - 0.2e1 * t6 * t23, -0.2e1 * t1 * t95 + 0.2e1 * t12 * t22 - 0.2e1 * t5 * t41 + 0.2e1 * t6 * t61, 0.2e1 * t6 * t1 + 0.2e1 * t12 * t5 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, t89 * t61, -t87 * t61, -t119 * t60, t104, 0, 0, 0, 0, 0, t28, t98, t28, t102 + t121, -t98, t1 * t69 - t138 * t2 + t139 * t6 + t7 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t138 * t59 + 0.2e1 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t61, 0, -t38, t37, -t38 * t89 + t93 * t87, t38 * t87 + t93 * t89, t103, -t38 * pkin(3) + (-t24 * t87 + t25 * t89) * qJD(4) + t103 * qJ(4), t139 * t41 - t22 * t69, -t102 + t121, -t98, t28, 0, -t138 * t26 + t81 * t23 + t39 * t59 + t100, t139 * t39 - t81 * t22 + t26 * t69 + t101, t12 * t59 - t138 * t5 + t44 * t23 + t32 * t40 + t100, t1 * t138 + t139 * t7 + t2 * t69 - t49 * t22 - t51 * t23 + t35 * t40 + t36 * t41 - t6 * t59, -t12 * t139 + t44 * t22 - t32 * t41 - t5 * t69 - t101, t1 * t51 + t12 * t32 + t2 * t49 - t6 * t35 + t7 * t36 + t5 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138 * t36 + t139 * t51 - t69 * t35 + t59 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, qJ(4) * t107, 0.2e1 * t132, 0.2e1 * t138 * t139 - 0.2e1 * t59 * t69, 0, 0, 0, t59 * t141, t139 * t141, -0.2e1 * t138 * t32 + 0.2e1 * t44 * t59, -0.2e1 * t138 * t35 + 0.2e1 * t139 * t49 + 0.2e1 * t36 * t69 - 0.2e1 * t51 * t59, -0.2e1 * t139 * t44 - 0.2e1 * t32 * t69, 0.2e1 * t44 * t32 - 0.2e1 * t51 * t35 + 0.2e1 * t49 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t129, 0, t38, 0, 0, 0, 0, 0, t23, -t22, t23, 0, t22, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t139, t59, 0, -t139, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, t61, t4, t3, t4 + 0.2e1 * t135, pkin(5) * t22 - t23 * qJ(6) - t40 * qJD(6), -0.2e1 * t116 - t3 + 0.2e1 * t118, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t139, -t59, 0, t139, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t59, 0, -t36, t35, -t36, -pkin(5) * t139 - t59 * qJ(6) + qJD(6) * t138, -t35, -t36 * pkin(5) - t35 * qJ(6) + t51 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, qJ(6) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
