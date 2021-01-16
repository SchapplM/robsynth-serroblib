% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP1
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
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:49
% EndTime: 2021-01-16 02:38:57
% DurationCPUTime: 1.65s
% Computational Cost: add. (1894->207), mult. (4754->386), div. (0->0), fcn. (4607->10), ass. (0->116)
t118 = cos(pkin(11));
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t124 = qJ(4) + pkin(8);
t94 = qJD(3) * t124;
t43 = t71 * qJD(4) - t68 * t94;
t65 = sin(pkin(11));
t79 = -t68 * qJD(4) - t71 * t94;
t23 = t118 * t43 + t65 * t79;
t128 = t65 * t68;
t95 = t118 * t71;
t50 = -t95 + t128;
t96 = t118 * t68;
t51 = t65 * t71 + t96;
t60 = -pkin(3) * t71 - pkin(2);
t80 = -pkin(4) * t50 + pkin(9) * t51 - t60;
t139 = qJD(5) * t80 - t23;
t55 = t124 * t71;
t37 = t118 * t55 - t124 * t128;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t123 = t70 * t37 - t67 * t80;
t66 = sin(pkin(6));
t72 = cos(qJ(2));
t126 = t66 * t72;
t110 = t67 * t126;
t119 = cos(pkin(6));
t69 = sin(qJ(2));
t127 = t66 * t69;
t46 = t119 * t68 + t71 * t127;
t81 = t119 * t71 - t68 * t127;
t26 = t118 * t46 + t65 * t81;
t17 = t70 * t26 - t110;
t117 = qJD(2) * t69;
t103 = t66 * t117;
t116 = qJD(2) * t72;
t102 = t66 * t116;
t35 = t81 * qJD(3) + t71 * t102;
t75 = t46 * qJD(3) + t68 * t102;
t74 = t118 * t35 - t65 * t75;
t84 = t70 * t126 + t67 * t26;
t7 = t84 * qJD(5) - t67 * t103 - t70 * t74;
t62 = qJD(5) * t70;
t8 = qJD(5) * t110 + t70 * t103 - t26 * t62 - t67 * t74;
t138 = (-t17 * t70 - t67 * t84) * qJD(5) + t7 * t67 - t8 * t70;
t121 = qJ(6) * t51;
t97 = -t37 * t67 - t70 * t80;
t11 = pkin(5) * t50 - t70 * t121 + t97;
t12 = -t67 * t121 + t123;
t44 = t51 * qJD(3);
t135 = t44 * pkin(5);
t100 = qJD(5) * t121;
t115 = qJD(3) * t68;
t107 = pkin(3) * t115;
t45 = qJD(3) * t95 - t65 * t115;
t77 = pkin(4) * t44 - pkin(9) * t45 + t107;
t20 = t70 * t77;
t78 = qJ(6) * t45 + qJD(5) * t37 + qJD(6) * t51;
t73 = -t78 * t70 + t20 + (t100 + t139) * t67;
t3 = t73 + t135;
t108 = t139 * t70 - t67 * t77;
t4 = t70 * t100 + t78 * t67 + t108;
t137 = -t3 * t70 + t4 * t67 + (t11 * t67 - t12 * t70) * qJD(5);
t136 = 0.2e1 * qJD(5);
t134 = t70 * pkin(5);
t15 = t118 * t75 + t35 * t65;
t25 = -t118 * t81 + t46 * t65;
t133 = t25 * t15;
t132 = t45 * t67;
t131 = t51 * t67;
t130 = t51 * t70;
t59 = -t118 * pkin(3) - pkin(4);
t53 = t59 - t134;
t129 = t53 * t70;
t125 = t67 * t70;
t63 = t67 ^ 2;
t64 = t70 ^ 2;
t122 = t63 - t64;
t58 = pkin(3) * t65 + pkin(9);
t120 = qJ(6) + t58;
t114 = qJD(3) * t71;
t113 = qJD(3) * t72;
t112 = qJD(5) * t67;
t111 = -0.2e1 * pkin(2) * qJD(3);
t109 = t59 * t136;
t106 = pkin(5) * t112;
t105 = t68 * t113;
t104 = t25 * t112;
t101 = t67 * t62;
t99 = -t53 + t134;
t98 = -0.4e1 * t51 * t125;
t22 = -t118 * t79 + t65 * t43;
t36 = t124 * t96 + t55 * t65;
t93 = t122 * qJD(5);
t90 = pkin(5) * t63 + t129;
t88 = -t17 * t67 + t70 * t84;
t86 = -t44 * t58 + t45 * t59;
t85 = t50 * t58 - t51 * t59;
t83 = t44 * t67 + t50 * t62;
t32 = t51 * t62 + t132;
t82 = t51 * t112 - t45 * t70;
t49 = t51 ^ 2;
t48 = t120 * t70;
t47 = t120 * t67;
t39 = -t67 * qJD(6) - t120 * t62;
t38 = -t70 * qJD(6) + t120 * t112;
t29 = -t50 * t112 + t44 * t70;
t21 = pkin(5) * t131 + t36;
t13 = t32 * pkin(5) + t22;
t10 = -t15 * t70 + t104;
t9 = t15 * t67 + t25 * t62;
t6 = -t123 * qJD(5) - t67 * t23 + t20;
t5 = t37 * t112 + t108;
t2 = t15 * t131 + t32 * t25 - t44 * t84 + t8 * t50;
t1 = t15 * t130 - t17 * t44 - t82 * t25 + t7 * t50;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66 ^ 2 * t69 * t116 + 0.2e1 * t26 * t74 + 0.2e1 * t133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t7 - 0.2e1 * t8 * t84 + 0.2e1 * t133; 0, 0, -t103, -t102, 0, 0, 0, 0, 0, (-t71 * t117 - t105) * t66, (-t71 * t113 + t68 * t117) * t66, (t50 * t117 - t44 * t72) * t66, (t51 * t117 - t45 * t72) * t66, t15 * t51 + t25 * t45 - t26 * t44 - t74 * t50, -t66 * pkin(3) * t105 + t60 * t103 + t15 * t36 + t25 * t22 + t26 * t23 + t74 * t37, 0, 0, 0, 0, 0, t2, t1, t2, t1, t138 * t51 + t88 * t45, t11 * t8 - t12 * t7 + t13 * t25 + t15 * t21 - t17 * t4 - t3 * t84; 0, 0, 0, 0, 0.2e1 * t68 * t114, 0.2e1 * (-t68 ^ 2 + t71 ^ 2) * qJD(3), 0, 0, 0, t68 * t111, t71 * t111, 0.2e1 * t50 * t107 + 0.2e1 * t44 * t60, 0.2e1 * t51 * t107 + 0.2e1 * t45 * t60, 0.2e1 * t22 * t51 - 0.2e1 * t23 * t50 + 0.2e1 * t36 * t45 - 0.2e1 * t37 * t44, 0.2e1 * t60 * t107 + 0.2e1 * t22 * t36 + 0.2e1 * t23 * t37, 0.2e1 * t45 * t51 * t64 - 0.2e1 * t49 * t101, t122 * t49 * t136 + t45 * t98, 0.2e1 * t44 * t130 - 0.2e1 * t82 * t50, -0.2e1 * t44 * t131 - 0.2e1 * t32 * t50, 0.2e1 * t50 * t44, 0.2e1 * t22 * t131 + 0.2e1 * t32 * t36 + 0.2e1 * t97 * t44 + 0.2e1 * t6 * t50, -0.2e1 * t123 * t44 + 0.2e1 * t22 * t130 - 0.2e1 * t82 * t36 + 0.2e1 * t5 * t50, 0.2e1 * t11 * t44 + 0.2e1 * t13 * t131 + 0.2e1 * t32 * t21 + 0.2e1 * t3 * t50, -0.2e1 * t12 * t44 + 0.2e1 * t13 * t130 - 0.2e1 * t82 * t21 + 0.2e1 * t4 * t50, 0.2e1 * (-t11 * t70 - t12 * t67) * t45 + 0.2e1 * t137 * t51, 0.2e1 * t11 * t3 - 0.2e1 * t12 * t4 + 0.2e1 * t13 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t35, -t15, -t74, 0, (-t15 * t118 + t74 * t65) * pkin(3), 0, 0, 0, 0, 0, t10, t9, t10, t9, t88 * qJD(5) - t8 * t67 - t7 * t70, pkin(5) * t104 + t15 * t53 - t17 * t38 - t39 * t84 - t47 * t8 - t48 * t7; 0, 0, 0, 0, 0, 0, t114, -t115, 0, -pkin(8) * t114, pkin(8) * t115, -t22, -t23, (-t118 * t45 - t44 * t65) * pkin(3), (-t118 * t22 + t23 * t65) * pkin(3), t45 * t125 - t51 * t93, qJD(5) * t98 - t122 * t45, t83, t29, 0, -t22 * t70 + t86 * t67 + (t36 * t67 - t85 * t70) * qJD(5), t22 * t67 + t86 * t70 + (t36 * t70 + t85 * t67) * qJD(5), t53 * t132 - t13 * t70 + t39 * t50 - t44 * t47 + (t21 * t67 + t90 * t51) * qJD(5), t45 * t129 + t13 * t67 + t38 * t50 - t48 * t44 + (t99 * t131 + t21 * t70) * qJD(5), (-t39 * t51 + t45 * t47 - t4 + (-t48 * t51 - t11) * qJD(5)) * t70 + (t38 * t51 - t45 * t48 - t3 + (-t47 * t51 - t12) * qJD(5)) * t67, t21 * t106 + t11 * t39 - t12 * t38 + t13 * t53 - t3 * t47 - t4 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101, -0.2e1 * t93, 0, 0, 0, t67 * t109, t70 * t109, -0.2e1 * t99 * t112, t90 * t136, -0.2e1 * t38 * t70 - 0.2e1 * t39 * t67 + 0.2e1 * (t47 * t70 - t48 * t67) * qJD(5), 0.2e1 * t53 * t106 - 0.2e1 * t38 * t48 - 0.2e1 * t39 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, t107, 0, 0, 0, 0, 0, t29, -t83, t29, -t83, (-t63 - t64) * t45, -t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t67 + t39 * t70 + (t47 * t67 + t48 * t70) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, t8, t7, 0, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t32, t44, t6, t5, t73 + 0.2e1 * t135, t4, t82 * pkin(5), t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t112, 0, -t58 * t62, t58 * t112, t39, t38, -pkin(5) * t62, t39 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t62, -t112, -t62, 0, -t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t82, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t62, 0, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
