% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:22
% DurationCPUTime: 1.19s
% Computational Cost: add. (2089->183), mult. (6499->320), div. (0->0), fcn. (6212->12), ass. (0->110)
t134 = cos(qJ(2));
t71 = sin(pkin(5));
t74 = sin(qJ(3));
t75 = sin(qJ(2));
t78 = cos(qJ(3));
t40 = (t134 * t74 + t75 * t78) * t71;
t136 = qJD(2) + qJD(3);
t120 = pkin(2) * qJD(3);
t111 = t74 * t120;
t73 = sin(qJ(4));
t112 = pkin(2) * t73 * t74;
t121 = qJD(4) * t112 + t73 * t111;
t77 = cos(qJ(4));
t110 = t78 * t120;
t67 = pkin(2) * t78 + pkin(3);
t86 = qJD(4) * t67 + t110;
t135 = -t86 * t77 + t121;
t104 = t134 * t78;
t125 = t74 * t75;
t39 = (t104 - t125) * t71;
t25 = t39 * t73 + t40 * t77;
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t118 = cos(pkin(6));
t24 = t39 * t77 - t40 * t73;
t119 = cos(pkin(5));
t70 = sin(pkin(6));
t96 = t119 * t70;
t82 = t118 * t24 + t96;
t15 = t76 * t25 + t82 * t72;
t117 = qJD(2) * t71;
t90 = t134 * t117;
t28 = -t78 * t90 + (-qJD(3) * t104 + t136 * t125) * t71;
t29 = t136 * t40;
t10 = -t25 * qJD(4) + t73 * t28 - t77 * t29;
t69 = t70 ^ 2;
t133 = t10 * t69;
t16 = t119 * t118 - t70 * t24;
t132 = t16 * t70;
t116 = qJD(4) * t73;
t115 = qJD(4) * t77;
t124 = t74 * t77;
t79 = (-t74 * t115 + (-t73 * t78 - t124) * qJD(3)) * pkin(2);
t31 = -t67 * t116 + t79;
t131 = t31 * t72;
t44 = t77 * t67 - t112;
t43 = pkin(4) + t44;
t130 = t43 * t69;
t129 = t69 * t72;
t128 = t69 * t76;
t127 = t70 * t72;
t126 = t70 * t76;
t45 = pkin(2) * t124 + t73 * t67;
t68 = t70 * pkin(10);
t37 = t68 + t45;
t99 = t72 * t118;
t85 = -t76 * t37 - t43 * t99;
t98 = t76 * t118;
t13 = t85 * qJD(5) + t135 * t72 + t31 * t98;
t122 = t13 * t118 + t31 * t128;
t114 = qJD(5) * t72;
t113 = qJD(5) * t76;
t109 = pkin(3) * t116;
t108 = pkin(3) * t115;
t107 = t69 * t114;
t106 = t69 * t113;
t105 = t70 * t114;
t64 = t70 * t113;
t103 = qJD(5) * (-pkin(4) - t43);
t66 = pkin(3) * t77 + pkin(4);
t102 = qJD(5) * (-pkin(4) - t66);
t94 = qJD(5) * t118;
t89 = t76 * t94;
t12 = t37 * t114 + t135 * t76 - t31 * t99 - t43 * t89;
t101 = t12 * t118;
t41 = -pkin(4) * t89 + pkin(10) * t105;
t100 = t41 * t118;
t97 = (-pkin(3) - t67) * qJD(4);
t95 = qJD(5) * (-t43 - t66);
t93 = t69 * t109;
t92 = t76 * t109;
t91 = t72 * t106;
t57 = pkin(3) * t73 + t68;
t22 = -t66 * t89 - t76 * t108 + (qJD(5) * t57 + t118 * t109) * t72;
t88 = t22 * t118 + t72 * t93;
t84 = -t76 * t57 - t66 * t99;
t83 = -pkin(4) * t99 - pkin(10) * t126;
t53 = -0.2e1 * t91;
t52 = 0.2e1 * t91;
t50 = 0.2e1 * t70 * t89;
t49 = -0.2e1 * t94 * t127;
t46 = pkin(4) * t98 - pkin(10) * t127;
t42 = t83 * qJD(5);
t38 = 0.2e1 * (-t72 ^ 2 + t76 ^ 2) * t69 * qJD(5);
t36 = t42 * t118;
t35 = t41 * t126;
t32 = -t72 * t57 + t66 * t98;
t23 = t84 * qJD(5) + (-t72 * t77 - t73 * t98) * qJD(4) * pkin(3);
t20 = -t72 * t37 + t43 * t98;
t19 = t23 * t118;
t17 = t22 * t126;
t14 = -t72 * t25 + t82 * t76;
t9 = -t24 * qJD(4) + t77 * t28 + t73 * t29;
t8 = t12 * t126;
t5 = -qJD(5) * t15 + t10 * t98 + t72 * t9;
t4 = -t10 * t99 - t96 * t113 + t25 * t114 - t24 * t89 + t76 * t9;
t3 = t10 * t128 + t16 * t105 + t5 * t118;
t2 = -t10 * t129 + t4 * t118 + t16 * t64;
t1 = (-t4 * t76 - t5 * t72 + (-t14 * t76 - t15 * t72) * qJD(5)) * t70;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t28 * t40 - 0.2e1 * t29 * t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t24 - 0.2e1 * t25 * t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10 * t132 + 0.2e1 * t14 * t5 - 0.2e1 * t15 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t117, -t90, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, 0, (-t28 * t74 - t29 * t78 + (-t39 * t74 + t40 * t78) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t10, t9, 0, t10 * t44 - t135 * t25 + t24 * t31 - t9 * t45, 0, 0, 0, 0, 0, 0, t3, t2, t1, t10 * t130 - t12 * t15 + t13 * t14 - t31 * t132 + t20 * t5 + t4 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t111, -0.2e1 * t110, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t31, 0.2e1 * t135, 0, -0.2e1 * t135 * t45 + 0.2e1 * t44 * t31, t52, t38, t50, t53, t49, 0, -0.2e1 * t43 * t107 + 0.2e1 * t122, 0.2e1 * t101 + 0.2e1 * (-t43 * t113 - t131) * t69, -0.2e1 * t8 + 0.2e1 * (-t13 * t72 + (-t20 * t76 + t72 * t85) * qJD(5)) * t70, 0.2e1 * t12 * t85 + 0.2e1 * t13 * t20 + 0.2e1 * t31 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, (t10 * t77 - t73 * t9 + (-t24 * t73 + t25 * t77) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, t3, t2, t1, t109 * t132 + t66 * t133 + t14 * t23 - t15 * t22 + t32 * t5 + t4 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t97 + t79, (t97 - t110) * t77 + t121, 0, ((-qJD(4) * t44 - t121) * t73 + (t45 * qJD(4) + t86 * t73 + t31) * t77) * pkin(3), t52, t38, t50, t53, t49, 0, t19 + (t72 * t95 - t92) * t69 + t122, t101 + (t76 * t95 - t131) * t69 + t88, -t17 - t8 + ((-t13 - t23) * t72 + ((-t20 - t32) * t76 + (t85 + t84) * t72) * qJD(5)) * t70, t12 * t84 + t13 * t32 + t20 * t23 + t85 * t22 + (-t43 * t109 + t31 * t66) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t109, -0.2e1 * t108, 0, 0, t52, t38, t50, t53, t49, 0, 0.2e1 * t19 + 0.2e1 * (-t66 * t114 - t92) * t69, -0.2e1 * t66 * t106 + 0.2e1 * t88, -0.2e1 * t17 + 0.2e1 * (-t23 * t72 + (-t32 * t76 + t72 * t84) * qJD(5)) * t70, 0.2e1 * t22 * t84 + 0.2e1 * t23 * t32 - 0.2e1 * t66 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, t1, pkin(4) * t133 + t14 * t42 - t15 * t41 + t4 * t83 + t46 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t135, 0, 0, t52, t38, t50, t53, t49, 0, t103 * t129 + t122 + t36, t101 + t100 + (t76 * t103 - t131) * t69, -t35 - t8 + ((-t13 - t42) * t72 + ((-t20 - t46) * t76 + (t85 + t83) * t72) * qJD(5)) * t70, pkin(4) * t31 * t69 + t12 * t83 + t13 * t46 + t20 * t42 + t41 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t108, 0, 0, t52, t38, t50, t53, t49, 0, t19 + t36 + (t72 * t102 - t92) * t69, t102 * t128 + t100 + t88, -t17 - t35 + ((-t23 - t42) * t72 + ((-t32 - t46) * t76 + (t84 + t83) * t72) * qJD(5)) * t70, -pkin(4) * t93 + t22 * t83 + t23 * t46 + t32 * t42 + t41 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t38, t50, t53, t49, 0, -0.2e1 * pkin(4) * t107 + 0.2e1 * t36, -0.2e1 * pkin(4) * t106 + 0.2e1 * t100, -0.2e1 * t35 + 0.2e1 * (-t42 * t72 + (-t46 * t76 + t72 * t83) * qJD(5)) * t70, 0.2e1 * t41 * t83 + 0.2e1 * t42 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t105, 0, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t105, 0, t23, t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t105, 0, t42, t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
