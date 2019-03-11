% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:18
% EndTime: 2019-03-08 20:12:22
% DurationCPUTime: 1.15s
% Computational Cost: add. (1692->163), mult. (4273->304), div. (0->0), fcn. (4349->10), ass. (0->107)
t130 = cos(qJ(4));
t62 = cos(pkin(11));
t100 = t130 * t62;
t60 = sin(pkin(11));
t64 = sin(qJ(4));
t78 = -t64 * t60 + t100;
t38 = t78 * qJD(4);
t44 = t130 * t60 + t64 * t62;
t65 = cos(qJ(5));
t55 = qJD(5) * t65;
t63 = sin(qJ(5));
t80 = t38 * t63 + t44 * t55;
t53 = -t62 * pkin(3) - pkin(2);
t33 = -pkin(4) * t78 - t44 * pkin(9) + t53;
t120 = pkin(8) + qJ(3);
t49 = t120 * t60;
t50 = t120 * t62;
t36 = t130 * t50 - t64 * t49;
t141 = t63 * t33 + t65 * t36;
t118 = t60 ^ 2 + t62 ^ 2;
t58 = t63 ^ 2;
t59 = t65 ^ 2;
t117 = t58 - t59;
t95 = t117 * qJD(5);
t125 = t44 * t65;
t66 = cos(qJ(2));
t114 = qJD(2) * t66;
t61 = sin(pkin(6));
t104 = t61 * t114;
t129 = sin(qJ(2));
t102 = t61 * t129;
t115 = cos(pkin(6));
t72 = t62 * t102 + t115 * t60;
t73 = -t60 * t102 + t115 * t62;
t24 = t130 * t72 + t64 * t73;
t14 = t24 * qJD(4) + t44 * t104;
t124 = t61 * t66;
t108 = t63 * t124;
t19 = t65 * t24 - t108;
t23 = -t130 * t73 + t64 * t72;
t39 = t44 * qJD(4);
t136 = t23 * qJD(4) - t78 * t104;
t18 = t65 * t124 + t63 * t24;
t96 = qJD(2) * t129;
t93 = t61 * t96;
t6 = t18 * qJD(5) + t136 * t65 - t63 * t93;
t113 = qJD(5) * t63;
t122 = t65 * t38;
t79 = t44 * t113 - t122;
t140 = -t14 * t125 + t19 * t39 + t79 * t23 + t6 * t78;
t101 = t130 * t49;
t21 = qJD(4) * t101 - qJD(3) * t100 + (qJD(3) * t60 + qJD(4) * t50) * t64;
t32 = t39 * pkin(4) - t38 * pkin(9);
t4 = -qJD(5) * t141 + t63 * t21 + t65 * t32;
t7 = -qJD(5) * t108 - t136 * t63 + t24 * t55 - t65 * t93;
t139 = -t6 * t63 - t7 * t65 + (t18 * t63 + t19 * t65) * qJD(5);
t112 = t78 * qJD(6);
t116 = qJ(6) * t39;
t3 = t36 * t113 + t65 * t21 - t63 * t32 - t33 * t55;
t1 = -t112 - t3 + t116;
t11 = -qJ(6) * t78 + t141;
t82 = t33 * t65 - t36 * t63;
t12 = pkin(5) * t78 - t82;
t131 = t39 * pkin(5);
t2 = -t131 - t4;
t138 = t1 * t63 - t2 * t65 + (t11 * t65 + t12 * t63) * qJD(5);
t88 = t65 * pkin(5) + t63 * qJ(6);
t137 = t88 * qJD(5) - t65 * qJD(6);
t135 = 0.2e1 * t53;
t134 = 0.2e1 * qJD(6);
t133 = pkin(9) * t39;
t132 = pkin(9) * t78;
t127 = t44 * t38;
t126 = t44 * t63;
t123 = t63 * t39;
t121 = t65 * t39;
t111 = t63 * qJD(6);
t109 = -0.2e1 * pkin(4) * qJD(5);
t107 = pkin(9) * t113;
t106 = pkin(9) * t55;
t103 = t63 * t55;
t99 = t118 * t66;
t98 = -0.4e1 * t63 * t125;
t94 = 0.2e1 * t118 * qJD(3);
t92 = -pkin(4) * t38 - t133;
t91 = pkin(4) * t44 - t132;
t35 = t64 * t50 + t101;
t87 = pkin(5) * t63 - qJ(6) * t65;
t85 = -t11 * t63 + t12 * t65;
t84 = t18 * t65 - t19 * t63;
t31 = -t55 * t78 + t123;
t29 = t113 * t78 + t121;
t77 = t14 * t126 - t18 * t39 + t80 * t23 + t78 * t7;
t76 = t118 * t129;
t45 = -pkin(4) - t88;
t22 = t44 * qJD(3) + t36 * qJD(4);
t5 = t137 * t44 + t87 * t38 + t22;
t75 = -t5 + (t44 * t45 + t132) * qJD(5);
t15 = t87 * t44 + t35;
t37 = -pkin(5) * t113 + qJ(6) * t55 + t111;
t71 = -qJD(5) * t15 + t37 * t44 - t38 * t45 + t133;
t69 = t85 * qJD(5) + t1 * t65 + t2 * t63;
t68 = t84 * qJD(5) - t6 * t65 + t7 * t63;
t41 = t44 ^ 2;
t9 = t23 * t113 - t14 * t65;
t8 = t14 * t63 + t23 * t55;
t10 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t129 + t76) * t61 ^ 2 * t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t23 + 0.2e1 * t18 * t7 - 0.2e1 * t19 * t6; 0, 0, -t93, -t104, -t62 * t93, t60 * t93, t61 * qJD(2) * t99 (t76 * qJD(3) + (-t129 * pkin(2) + qJ(3) * t99) * qJD(2)) * t61, 0, 0, 0, 0, 0 (-t39 * t66 - t78 * t96) * t61 (-t38 * t66 + t44 * t96) * t61, 0, 0, 0, 0, 0, t77, -t140, t77, -t139 * t44 + t84 * t38, t140, t1 * t19 - t11 * t6 + t12 * t7 + t14 * t15 + t18 * t2 + t23 * t5; 0, 0, 0, 0, 0, 0, t94, qJ(3) * t94, 0.2e1 * t127, 0.2e1 * t38 * t78 - 0.2e1 * t39 * t44, 0, 0, 0, t39 * t135, t38 * t135, -0.2e1 * t41 * t103 + 0.2e1 * t59 * t127, t38 * t98 + 0.2e1 * t41 * t95, 0.2e1 * t44 * t121 + 0.2e1 * t78 * t79, -0.2e1 * t44 * t123 + 0.2e1 * t78 * t80, -0.2e1 * t78 * t39, 0.2e1 * t22 * t126 + 0.2e1 * t80 * t35 + 0.2e1 * t82 * t39 - 0.2e1 * t4 * t78, 0.2e1 * t22 * t125 - 0.2e1 * t141 * t39 - 0.2e1 * t3 * t78 - 0.2e1 * t79 * t35, -0.2e1 * t12 * t39 + 0.2e1 * t5 * t126 + 0.2e1 * t80 * t15 + 0.2e1 * t2 * t78, -0.2e1 * t138 * t44 + 0.2e1 * t85 * t38, -0.2e1 * t1 * t78 + 0.2e1 * t11 * t39 - 0.2e1 * t5 * t125 + 0.2e1 * t79 * t15, 0.2e1 * t1 * t11 + 0.2e1 * t12 * t2 + 0.2e1 * t15 * t5; 0, 0, 0, 0, 0, 0, 0, t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t38, 0, 0, 0, 0, 0, t29, -t31, t29 (-t58 - t59) * t38, t31, t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t136, 0, 0, 0, 0, 0, t9, t8, t9, t68, -t8, t68 * pkin(9) + t14 * t45 - t23 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, -t22, t21, t63 * t122 - t44 * t95, qJD(5) * t98 - t117 * t38, t31, t29, 0, -t22 * t65 + t92 * t63 + (t35 * t63 - t91 * t65) * qJD(5), t22 * t63 + t92 * t65 + (t35 * t65 + t91 * t63) * qJD(5), -t71 * t63 + t75 * t65, t69, t75 * t63 + t71 * t65, t69 * pkin(9) - t15 * t37 + t5 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103, -0.2e1 * t95, 0, 0, 0, t63 * t109, t65 * t109, 0.2e1 * t45 * t113 + 0.2e1 * t37 * t65, 0, 0.2e1 * t37 * t63 - 0.2e1 * t45 * t55, -0.2e1 * t45 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t7, 0, -t6, -pkin(5) * t7 - qJ(6) * t6 + qJD(6) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t80, t39, t4, t3, t4 + 0.2e1 * t131, -t88 * t38 + (t87 * qJD(5) - t111) * t44, -0.2e1 * t112 - t3 + 0.2e1 * t116, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t55, -t113, 0, t55, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t113, 0, -t106, t107, -t106, -t137, -t107, -t137 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, qJ(6) * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t79, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
