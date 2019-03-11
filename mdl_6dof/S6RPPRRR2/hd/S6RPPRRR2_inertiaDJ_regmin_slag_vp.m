% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:27
% EndTime: 2019-03-09 02:21:29
% DurationCPUTime: 0.98s
% Computational Cost: add. (1558->156), mult. (3592->264), div. (0->0), fcn. (3603->10), ass. (0->106)
t131 = cos(qJ(4));
t70 = sin(pkin(10)) * pkin(1) + qJ(3);
t132 = pkin(7) + t70;
t78 = sin(pkin(11));
t54 = t132 * t78;
t79 = cos(pkin(11));
t55 = t132 * t79;
t83 = sin(qJ(4));
t39 = t131 * t55 - t83 * t54;
t85 = cos(qJ(5));
t28 = t85 * t39;
t103 = t131 * t79;
t123 = t83 * t78;
t57 = -t103 + t123;
t58 = t131 * t78 + t83 * t79;
t63 = -cos(pkin(10)) * pkin(1) - t79 * pkin(3) - pkin(2);
t29 = t57 * pkin(4) - t58 * pkin(8) + t63;
t82 = sin(qJ(5));
t119 = t82 * t29 + t28;
t116 = qJD(5) * t82;
t106 = t58 * t116;
t98 = qJD(4) * t131;
t52 = qJD(4) * t123 - t79 * t98;
t121 = t85 * t52;
t138 = -t121 - t106;
t77 = t85 ^ 2;
t117 = t82 ^ 2 - t77;
t97 = t117 * qJD(5);
t112 = qJD(5) + qJD(6);
t53 = t58 * qJD(4);
t137 = 0.2e1 * t53;
t136 = pkin(8) + pkin(9);
t135 = t53 * pkin(5);
t134 = t57 * pkin(5);
t115 = qJD(5) * t85;
t20 = t54 * t98 - qJD(3) * t103 + (qJD(3) * t78 + qJD(4) * t55) * t83;
t40 = t53 * pkin(4) + t52 * pkin(8);
t6 = -t29 * t115 + t39 * t116 + t85 * t20 - t82 * t40;
t125 = t82 * t52;
t86 = t58 * t115 - t125;
t5 = -t86 * pkin(9) - t6;
t84 = cos(qJ(6));
t133 = t84 * t5;
t130 = t58 * t52;
t129 = t58 * t82;
t128 = t58 * t85;
t13 = -pkin(9) * t129 + t119;
t81 = sin(qJ(6));
t127 = t81 * t13;
t126 = t81 * t82;
t124 = t82 * t53;
t122 = t84 * t13;
t120 = t85 * t53;
t118 = t58 * t120 - t57 * t121;
t114 = qJD(6) * t81;
t113 = qJD(6) * t84;
t111 = -0.2e1 * pkin(4) * qJD(5);
t110 = t82 * t121;
t109 = pkin(5) * t116;
t108 = pkin(5) * t114;
t107 = pkin(5) * t113;
t105 = t82 * t115;
t100 = t82 * t20 + t85 * t40;
t4 = pkin(9) * t121 + t135 + (-t28 + (pkin(9) * t58 - t29) * t82) * qJD(5) + t100;
t104 = t84 * t4 - t81 * t5;
t99 = t85 * t29 - t82 * t39;
t12 = -pkin(9) * t128 + t134 + t99;
t102 = -t12 - t134;
t101 = qJD(5) * t136;
t96 = 0.2e1 * (t78 ^ 2 + t79 ^ 2) * qJD(3);
t95 = pkin(4) * t52 - pkin(8) * t53;
t94 = pkin(4) * t58 + pkin(8) * t57;
t38 = t131 * t54 + t83 * t55;
t93 = t84 * t12 - t127;
t92 = t81 * t12 + t122;
t44 = t112 * t126 - t85 * t113 - t84 * t115;
t60 = t81 * t85 + t84 * t82;
t91 = t44 * t57 - t60 * t53;
t45 = t112 * t60;
t59 = -t84 * t85 + t126;
t90 = t57 * t45 + t53 * t59;
t89 = t57 * t52 - t53 * t58;
t64 = t136 * t82;
t65 = t136 * t85;
t88 = -t84 * t64 - t81 * t65;
t87 = -t81 * t64 + t84 * t65;
t36 = t57 * t115 + t124;
t21 = t58 * qJD(3) + t39 * qJD(4);
t73 = -t85 * pkin(5) - pkin(4);
t62 = t85 * t101;
t61 = t82 * t101;
t56 = t58 ^ 2;
t43 = t57 * t137;
t34 = t57 * t116 - t120;
t33 = t59 * t58;
t32 = t60 * t58;
t24 = pkin(5) * t129 + t38;
t19 = -t87 * qJD(6) + t81 * t61 - t84 * t62;
t18 = -t88 * qJD(6) + t84 * t61 + t81 * t62;
t14 = t86 * pkin(5) + t21;
t11 = -t114 * t129 + (t112 * t128 - t125) * t84 + t138 * t81;
t10 = t84 * t121 - t81 * t125 + t45 * t58;
t7 = -qJD(5) * t119 + t100;
t2 = -t92 * qJD(6) + t104;
t1 = -t93 * qJD(6) - t81 * t4 - t133;
t3 = [0, 0, 0, 0, 0, 0, t96, t70 * t96, -0.2e1 * t130, 0.2e1 * t89, 0, 0, 0, t63 * t137, -0.2e1 * t63 * t52, -0.2e1 * t56 * t105 - 0.2e1 * t77 * t130, 0.4e1 * t58 * t110 + 0.2e1 * t56 * t97, -0.2e1 * t57 * t106 + 0.2e1 * t118, -0.2e1 * t58 * t124 - 0.2e1 * t86 * t57, t43, 0.2e1 * t21 * t129 + 0.2e1 * t86 * t38 + 0.2e1 * t99 * t53 + 0.2e1 * t7 * t57, -0.2e1 * t119 * t53 + 0.2e1 * t21 * t128 + 0.2e1 * t138 * t38 + 0.2e1 * t6 * t57, 0.2e1 * t33 * t10, 0.2e1 * t10 * t32 + 0.2e1 * t33 * t11, -0.2e1 * t10 * t57 - 0.2e1 * t33 * t53, -0.2e1 * t57 * t11 - 0.2e1 * t53 * t32, t43, 0.2e1 * t24 * t11 + 0.2e1 * t14 * t32 + 0.2e1 * t2 * t57 + 0.2e1 * t93 * t53, 0.2e1 * t1 * t57 - 0.2e1 * t24 * t10 - 0.2e1 * t14 * t33 - 0.2e1 * t92 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t85 + t118, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, 0, 0, 0, 0, 0, -t34, -t36, 0, 0, 0, 0, 0, -t90, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, -t21, t20, -t58 * t97 - t110, -0.4e1 * t58 * t105 + t117 * t52, t36, -t34, 0, -t21 * t85 + t95 * t82 + (t38 * t82 - t94 * t85) * qJD(5), t21 * t82 + t95 * t85 + (t38 * t85 + t94 * t82) * qJD(5), -t10 * t60 + t33 * t44, t10 * t59 - t60 * t11 + t44 * t32 + t33 * t45, -t91, -t90, 0, t32 * t109 + t73 * t11 + t14 * t59 + t19 * t57 + t24 * t45 + t88 * t53, -t73 * t10 - t33 * t109 + t14 * t60 + t18 * t57 - t24 * t44 - t87 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t52, 0, 0, 0, 0, 0, t34, t36, 0, 0, 0, 0, 0, t90, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t105, -0.2e1 * t97, 0, 0, 0, t82 * t111, t85 * t111, -0.2e1 * t60 * t44, 0.2e1 * t44 * t59 - 0.2e1 * t60 * t45, 0, 0, 0, 0.2e1 * t59 * t109 + 0.2e1 * t73 * t45, 0.2e1 * t60 * t109 - 0.2e1 * t73 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, -t86, t53, t7, t6, 0, 0, -t10, -t11, t53, t84 * t135 + (t102 * t81 - t122) * qJD(6) + t104, -t133 + (-t4 - t135) * t81 + (t102 * t84 + t127) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t138, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t115, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t116, 0, -pkin(8) * t115, pkin(8) * t116, 0, 0, -t44, -t45, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108, -0.2e1 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, t53, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
