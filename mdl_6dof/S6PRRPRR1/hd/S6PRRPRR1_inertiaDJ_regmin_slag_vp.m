% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:22
% EndTime: 2021-01-16 03:28:28
% DurationCPUTime: 1.28s
% Computational Cost: add. (1957->169), mult. (4847->324), div. (0->0), fcn. (5031->12), ass. (0->118)
t85 = cos(pkin(12));
t113 = pkin(3) * t85 + pkin(4);
t140 = cos(qJ(5));
t83 = sin(pkin(12));
t142 = pkin(3) * t83;
t88 = sin(qJ(5));
t143 = -t140 * t113 + t88 * t142;
t91 = cos(qJ(6));
t82 = t91 ^ 2;
t87 = sin(qJ(6));
t129 = t87 ^ 2 - t82;
t108 = t129 * qJD(6);
t131 = qJ(4) + pkin(8);
t89 = sin(qJ(3));
t73 = t131 * t89;
t92 = cos(qJ(3));
t74 = t131 * t92;
t47 = -t85 * t73 - t74 * t83;
t67 = t83 * t92 + t85 * t89;
t38 = -pkin(9) * t67 + t47;
t48 = -t83 * t73 + t85 * t74;
t66 = t83 * t89 - t85 * t92;
t39 = -pkin(9) * t66 + t48;
t21 = t140 * t39 + t38 * t88;
t109 = qJD(3) * t131;
t59 = t92 * qJD(4) - t109 * t89;
t60 = -t89 * qJD(4) - t109 * t92;
t36 = -t59 * t83 + t60 * t85;
t125 = qJD(3) * t92;
t126 = qJD(3) * t89;
t63 = t125 * t85 - t126 * t83;
t96 = pkin(9) * t63 - t36;
t37 = t59 * t85 + t60 * t83;
t62 = t67 * qJD(3);
t97 = -pkin(9) * t62 + t37;
t11 = qJD(5) * t21 + t140 * t96 + t88 * t97;
t20 = -t140 * t38 + t39 * t88;
t80 = qJD(6) * t91;
t141 = t11 * t87 + t20 * t80;
t98 = -t140 * t66 - t67 * t88;
t24 = qJD(5) * t98 + t140 * t63 - t88 * t62;
t44 = t140 * t67 - t66 * t88;
t139 = t44 * t24;
t138 = t44 * t87;
t137 = t44 * t91;
t84 = sin(pkin(6));
t90 = sin(qJ(2));
t136 = t84 * t90;
t93 = cos(qJ(2));
t135 = t84 * t93;
t25 = qJD(5) * t44 + t140 * t62 + t88 * t63;
t134 = t87 * t25;
t133 = t91 * t24;
t132 = t91 * t25;
t95 = t113 * t88 + t140 * t142;
t55 = t95 * qJD(5);
t57 = -pkin(5) + t143;
t130 = t55 * t87 + t57 * t80;
t128 = qJD(2) * t90;
t127 = qJD(2) * t93;
t124 = qJD(3) * t93;
t123 = qJD(5) * t88;
t122 = qJD(6) * t87;
t120 = -0.2e1 * pkin(2) * qJD(3);
t119 = pkin(5) * t122;
t118 = pkin(5) * t80;
t79 = pkin(3) * t126;
t117 = t89 * t124;
t116 = t84 * t128;
t115 = t84 * t127;
t114 = t87 * t80;
t78 = -pkin(3) * t92 - pkin(2);
t49 = pkin(4) * t62 + t79;
t112 = -0.4e1 * t87 * t137;
t111 = t57 * t122 - t55 * t91;
t110 = qJD(5) * t140;
t53 = pkin(4) * t66 + t78;
t19 = -pkin(5) * t98 - pkin(10) * t44 + t53;
t107 = t19 * t91 - t21 * t87;
t106 = t19 * t87 + t21 * t91;
t58 = pkin(10) + t95;
t105 = -t44 * t57 - t58 * t98;
t86 = cos(pkin(6));
t64 = -t136 * t89 + t86 * t92;
t65 = t136 * t92 + t86 * t89;
t40 = t64 * t85 - t65 * t83;
t41 = t64 * t83 + t65 * t85;
t23 = t140 * t41 + t40 * t88;
t104 = t135 * t91 + t23 * t87;
t103 = t135 * t87 - t23 * t91;
t101 = t24 * t87 + t44 * t80;
t100 = t122 * t44 - t133;
t16 = -t80 * t98 + t134;
t99 = -t122 * t98 - t132;
t54 = t143 * qJD(5);
t94 = t24 * t57 - t25 * t58 + t44 * t55 - t54 * t98;
t75 = 0.2e1 * t114;
t71 = -0.2e1 * t108;
t46 = -qJD(3) * t65 - t115 * t89;
t45 = qJD(3) * t64 + t115 * t92;
t42 = t44 ^ 2;
t30 = t45 * t85 + t46 * t83;
t29 = -t45 * t83 + t46 * t85;
t22 = -t140 * t40 + t41 * t88;
t17 = t20 * t122;
t14 = -t108 * t44 + t133 * t87;
t13 = pkin(5) * t25 - pkin(10) * t24 + t49;
t12 = qJD(6) * t112 - t129 * t24;
t10 = -t110 * t38 + t123 * t39 - t140 * t97 + t88 * t96;
t8 = qJD(5) * t23 - t140 * t29 + t88 * t30;
t7 = -t110 * t40 + t123 * t41 - t140 * t30 - t29 * t88;
t6 = t122 * t22 - t8 * t91;
t5 = t22 * t80 + t8 * t87;
t4 = qJD(6) * t103 + t116 * t91 + t87 * t7;
t3 = qJD(6) * t104 - t116 * t87 + t91 * t7;
t2 = -qJD(6) * t106 + t87 * t10 + t91 * t13;
t1 = -qJD(6) * t107 + t91 * t10 - t87 * t13;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t127 * t84 ^ 2 * t90 + 0.2e1 * t29 * t40 + 0.2e1 * t30 * t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t116, -t115, 0, 0, 0, 0, 0, (-t128 * t92 - t117) * t84, (-t124 * t92 + t128 * t89) * t84, (t128 * t66 - t62 * t93) * t84, (t128 * t67 - t63 * t93) * t84, -t29 * t67 - t30 * t66 - t40 * t63 - t41 * t62, t29 * t47 + t30 * t48 + t40 * t36 + t41 * t37 + (-pkin(3) * t117 + t128 * t78) * t84, 0, 0, 0, 0, 0, (-t128 * t98 - t25 * t93) * t84, (t128 * t44 - t24 * t93) * t84, 0, 0, 0, 0, 0, t101 * t22 - t104 * t25 + t138 * t8 - t4 * t98, -t100 * t22 + t103 * t25 + t137 * t8 - t3 * t98; 0, 0, 0, 0, 0.2e1 * t89 * t125, 0.2e1 * (-t89 ^ 2 + t92 ^ 2) * qJD(3), 0, 0, 0, t89 * t120, t92 * t120, 0.2e1 * t62 * t78 + 0.2e1 * t66 * t79, 0.2e1 * t63 * t78 + 0.2e1 * t67 * t79, -0.2e1 * t36 * t67 - 0.2e1 * t37 * t66 - 0.2e1 * t47 * t63 - 0.2e1 * t48 * t62, 0.2e1 * t36 * t47 + 0.2e1 * t37 * t48 + 0.2e1 * t78 * t79, 0.2e1 * t139, 0.2e1 * t24 * t98 - 0.2e1 * t25 * t44, 0, 0, 0, 0.2e1 * t25 * t53 - 0.2e1 * t49 * t98, 0.2e1 * t24 * t53 + 0.2e1 * t44 * t49, -0.2e1 * t114 * t42 + 0.2e1 * t139 * t82, 0.2e1 * t108 * t42 + t112 * t24, 0.2e1 * t100 * t98 + 0.2e1 * t132 * t44, 0.2e1 * t101 * t98 - 0.2e1 * t134 * t44, -0.2e1 * t98 * t25, 0.2e1 * t101 * t20 + 0.2e1 * t107 * t25 + 0.2e1 * t11 * t138 - 0.2e1 * t2 * t98, -0.2e1 * t1 * t98 - 0.2e1 * t100 * t20 - 0.2e1 * t106 * t25 + 0.2e1 * t11 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, t29, -t30, 0, (t29 * t85 + t30 * t83) * pkin(3), 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, t125, -t126, 0, -pkin(8) * t125, pkin(8) * t126, t36, -t37, (-t62 * t83 - t63 * t85) * pkin(3), (t36 * t85 + t37 * t83) * pkin(3), 0, 0, t24, -t25, 0, -t11, t10, t14, t12, t16, -t99, 0, t17 + (-qJD(6) * t105 - t11) * t91 + t94 * t87, t105 * t122 + t91 * t94 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, 0.2e1 * t54, t75, t71, 0, 0, 0, 0.2e1 * t111, 0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t63, 0, t79, 0, 0, 0, 0, 0, t25, t24, 0, 0, 0, 0, 0, -t99, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25, 0, -t11, t10, t14, t12, t16, -t99, 0, t17 + (-pkin(5) * t24 - pkin(10) * t25) * t87 + (-t11 + (-pkin(5) * t44 + pkin(10) * t98) * qJD(6)) * t91, pkin(5) * t100 + pkin(10) * t99 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t54, t75, t71, 0, 0, 0, t111 - t119, -t118 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t71, 0, 0, 0, -0.2e1 * t119, -0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t101, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t122, 0, t54 * t87 - t58 * t80, t122 * t58 + t54 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t122, 0, -pkin(10) * t80, pkin(10) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
