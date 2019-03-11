% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:26
% EndTime: 2019-03-08 21:38:30
% DurationCPUTime: 1.56s
% Computational Cost: add. (1868->238), mult. (5011->422), div. (0->0), fcn. (4732->10), ass. (0->114)
t141 = cos(qJ(5));
t80 = cos(pkin(11));
t83 = sin(qJ(3));
t136 = t80 * t83;
t85 = cos(qJ(3));
t104 = -pkin(3) * t85 - qJ(4) * t83;
t61 = -pkin(2) + t104;
t57 = t80 * t61;
t78 = sin(pkin(11));
t34 = -pkin(9) * t136 + t57 + (-pkin(8) * t78 - pkin(4)) * t85;
t140 = t78 * t83;
t135 = t80 * t85;
t71 = pkin(8) * t135;
t44 = t78 * t61 + t71;
t40 = -pkin(9) * t140 + t44;
t82 = sin(qJ(5));
t147 = t141 * t40 + t82 * t34;
t146 = t78 ^ 2 + t80 ^ 2;
t109 = qJD(5) * t141;
t126 = qJD(5) * t82;
t145 = t80 * t109 - t78 * t126;
t86 = cos(qJ(2));
t131 = qJD(2) * t86;
t79 = sin(pkin(6));
t117 = t79 * t131;
t84 = sin(qJ(2));
t138 = t79 * t84;
t81 = cos(pkin(6));
t52 = t83 * t138 - t81 * t85;
t144 = t52 * qJD(3) - t85 * t117;
t130 = qJD(3) * t83;
t120 = pkin(8) * t130;
t49 = -qJD(4) * t83 + (pkin(3) * t83 - qJ(4) * t85) * qJD(3);
t36 = t78 * t120 + t80 * t49;
t25 = (pkin(4) * t83 - pkin(9) * t135) * qJD(3) + t36;
t139 = t78 * t85;
t45 = t78 * t49;
t29 = t45 + (-pkin(8) * t136 - pkin(9) * t139) * qJD(3);
t6 = -qJD(5) * t147 + t141 * t25 - t82 * t29;
t143 = 0.2e1 * t145;
t142 = 2 * qJD(6);
t75 = t83 * pkin(8);
t137 = t79 * t86;
t134 = pkin(9) + qJ(4);
t129 = qJD(3) * t85;
t118 = t78 * t129;
t74 = pkin(8) * t129;
t54 = pkin(4) * t118 + t74;
t60 = pkin(4) * t140 + t75;
t132 = qJD(2) * t84;
t128 = qJD(3) * t86;
t127 = qJD(4) * t78;
t125 = qJD(6) * t85;
t124 = t80 * qJD(4);
t123 = pkin(8) * t139;
t122 = -0.2e1 * pkin(2) * qJD(3);
t121 = pkin(5) * t130;
t116 = t80 * t129;
t115 = t81 * t129;
t114 = t83 * t129;
t73 = -pkin(4) * t80 - pkin(3);
t113 = t141 * t80;
t112 = qJ(6) * t130;
t111 = t134 * t78;
t108 = t141 * qJD(4);
t106 = 0.2e1 * t146 * qJD(4);
t37 = -t80 * t120 + t45;
t103 = -t36 * t78 + t37 * t80;
t53 = t85 * t138 + t81 * t83;
t39 = t53 * qJD(3) + t83 * t117;
t59 = t141 * t78 + t82 * t80;
t102 = t145 * t52 + t39 * t59;
t101 = t141 * t111;
t100 = t80 * t137 + t53 * t78;
t64 = t134 * t80;
t22 = qJD(5) * t101 - t80 * t108 + (qJD(5) * t64 + t127) * t82;
t42 = -t82 * t111 + t141 * t64;
t99 = -t42 * t130 - t22 * t85;
t23 = t64 * t109 + t82 * t124 + (-t134 * t126 + t108) * t78;
t41 = t82 * t64 + t101;
t98 = -t41 * t130 + t23 * t85;
t97 = t141 * t34 - t82 * t40;
t5 = -t34 * t109 + t40 * t126 - t141 * t29 - t82 * t25;
t95 = t82 * t100;
t94 = -t84 * t130 + t85 * t131;
t93 = t100 * qJD(3);
t38 = -t78 * t137 + t53 * t80;
t91 = t141 * t100;
t11 = t82 * t38 + t91;
t27 = t59 * t129 + t145 * t83;
t87 = t80 * t115 + (t78 * t132 + t94 * t80) * t79;
t88 = t78 * t115 - (t80 * t132 - t94 * t78) * t79;
t4 = -qJD(5) * t95 + t38 * t109 + t141 * t88 + t82 * t87;
t47 = t59 * t83;
t92 = -t11 * t130 + t52 * t27 + t39 * t47 + t4 * t85;
t51 = t59 * qJD(5);
t12 = t141 * t38 - t95;
t26 = -t113 * t129 + t82 * t118 + t83 * t51;
t3 = qJD(5) * t91 + t38 * t126 - t141 * t87 + t82 * t88;
t48 = t83 * t113 - t82 * t140;
t90 = t12 * t130 + t26 * t52 + t3 * t85 - t39 * t48;
t58 = t78 * t82 - t113;
t43 = t57 - t123;
t32 = pkin(5) * t58 - qJ(6) * t59 + t73;
t24 = t52 * t39;
t15 = pkin(5) * t51 - qJ(6) * t145 - qJD(6) * t59;
t14 = pkin(5) * t47 - qJ(6) * t48 + t60;
t10 = t85 * pkin(5) - t97;
t9 = -qJ(6) * t85 + t147;
t8 = t39 * t58 + t51 * t52;
t7 = pkin(5) * t27 + qJ(6) * t26 - qJD(6) * t48 + t54;
t2 = -t121 - t6;
t1 = t112 - t5 - t125;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100 * t88 + 0.2e1 * t38 * t87 + 0.2e1 * t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t11 * t4 - 0.2e1 * t12 * t3 + 0.2e1 * t24; 0, 0, -t79 * t132, -t117, 0, 0, 0, 0, 0 (-t83 * t128 - t85 * t132) * t79 (-t85 * t128 + t83 * t132) * t79, t52 * t118 + t39 * t140 - t83 * t93 + t88 * t85, t52 * t116 - t38 * t130 + t39 * t136 + t87 * t85, -t38 * t118 + t93 * t135 + t88 * t136 - t87 * t140, -t100 * t36 + t38 * t37 + t39 * t75 - t88 * t43 + t87 * t44 + t52 * t74, 0, 0, 0, 0, 0, t92, -t90, t92, -t11 * t26 - t12 * t27 + t3 * t47 + t4 * t48, t90, t1 * t12 + t10 * t4 + t11 * t2 + t14 * t39 - t3 * t9 + t52 * t7; 0, 0, 0, 0, 0.2e1 * t114, 0.2e1 * (-t83 ^ 2 + t85 ^ 2) * qJD(3), 0, 0, 0, t83 * t122, t85 * t122, -0.2e1 * t36 * t85 + 0.2e1 * (t43 + 0.2e1 * t123) * t130, 0.2e1 * t37 * t85 + 0.2e1 * (-t44 + 0.2e1 * t71) * t130, 0.2e1 * (-t36 * t80 - t37 * t78) * t83 + 0.2e1 * (-t43 * t80 - t44 * t78) * t129, 0.2e1 * pkin(8) ^ 2 * t114 + 0.2e1 * t36 * t43 + 0.2e1 * t37 * t44, -0.2e1 * t48 * t26, 0.2e1 * t26 * t47 - 0.2e1 * t27 * t48, 0.2e1 * t48 * t130 + 0.2e1 * t26 * t85, -0.2e1 * t47 * t130 + 0.2e1 * t27 * t85, -0.2e1 * t114, 0.2e1 * t97 * t130 + 0.2e1 * t60 * t27 + 0.2e1 * t54 * t47 - 0.2e1 * t6 * t85, -0.2e1 * t130 * t147 - 0.2e1 * t60 * t26 + 0.2e1 * t54 * t48 - 0.2e1 * t5 * t85, -0.2e1 * t10 * t130 + 0.2e1 * t14 * t27 + 0.2e1 * t2 * t85 + 0.2e1 * t47 * t7, -0.2e1 * t1 * t47 - 0.2e1 * t10 * t26 + 0.2e1 * t2 * t48 - 0.2e1 * t27 * t9, -0.2e1 * t1 * t85 + 0.2e1 * t9 * t130 + 0.2e1 * t14 * t26 - 0.2e1 * t48 * t7, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t14 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t144, -t39 * t80, t39 * t78, -t146 * t144, -t39 * pkin(3) + t100 * t127 + t38 * t124 + (t78 * t88 + t80 * t87) * qJ(4), 0, 0, 0, 0, 0, t8, t102, t8, t11 * t145 - t12 * t51 + t3 * t58 + t4 * t59, -t102, t11 * t23 - t12 * t22 + t15 * t52 - t3 * t42 + t32 * t39 + t4 * t41; 0, 0, 0, 0, 0, 0, t129, -t130, 0, -t74, t120, t85 * t127 + (t104 * t78 - t71) * qJD(3), t85 * t124 + (t104 * t80 + t123) * qJD(3), t103, -pkin(3) * t74 + (-t43 * t78 + t44 * t80) * qJD(4) + t103 * qJ(4), t145 * t48 - t26 * t59, -t145 * t47 + t26 * t58 - t27 * t59 - t48 * t51, t59 * t130 - t145 * t85, -t58 * t130 + t51 * t85, 0, t27 * t73 + t51 * t60 + t54 * t58 + t98, t145 * t60 - t26 * t73 + t54 * t59 + t99, t14 * t51 + t15 * t47 + t27 * t32 + t58 * t7 + t98, -t1 * t58 + t10 * t145 + t2 * t59 + t22 * t47 + t23 * t48 - t26 * t41 - t27 * t42 - t51 * t9, -t14 * t145 - t15 * t48 + t26 * t32 - t59 * t7 - t99, t1 * t42 + t10 * t23 + t14 * t15 + t2 * t41 - t22 * t9 + t32 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, qJ(4) * t106, t59 * t143, -0.2e1 * t145 * t58 - 0.2e1 * t51 * t59, 0, 0, 0, 0.2e1 * t73 * t51, t73 * t143, 0.2e1 * t15 * t58 + 0.2e1 * t32 * t51, 0.2e1 * t145 * t41 + 0.2e1 * t22 * t58 + 0.2e1 * t23 * t59 - 0.2e1 * t42 * t51, -0.2e1 * t145 * t32 - 0.2e1 * t15 * t59, 0.2e1 * t15 * t32 - 0.2e1 * t22 * t42 + 0.2e1 * t23 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t116, 0, t74, 0, 0, 0, 0, 0, t27, -t26, t27, 0, t26, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t145, t51, 0, -t145, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t4, 0, -t3, -pkin(5) * t4 - qJ(6) * t3 + qJD(6) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, t130, t6, t5, t6 + 0.2e1 * t121, pkin(5) * t26 - qJ(6) * t27 - qJD(6) * t47, 0.2e1 * t112 - t5 - 0.2e1 * t125, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, -t51, 0, -t23, t22, -t23, -pkin(5) * t145 - qJ(6) * t51 - qJD(6) * t58, -t22, -pkin(5) * t23 - qJ(6) * t22 + qJD(6) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, qJ(6) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
