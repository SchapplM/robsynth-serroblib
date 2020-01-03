% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:39
% EndTime: 2019-12-31 21:11:43
% DurationCPUTime: 1.12s
% Computational Cost: add. (939->140), mult. (2162->220), div. (0->0), fcn. (1609->6), ass. (0->96)
t129 = pkin(1) * qJD(2);
t97 = cos(qJ(2));
t123 = t97 * t129;
t94 = sin(qJ(3));
t91 = t94 ^ 2;
t96 = cos(qJ(3));
t92 = t96 ^ 2;
t46 = (t91 + t92) * t123;
t146 = t96 * pkin(3) + t94 * qJ(4);
t93 = sin(qJ(5));
t138 = t94 * t93;
t140 = pkin(7) - pkin(8);
t145 = t140 * t138;
t139 = cos(qJ(5));
t117 = t94 * t139;
t87 = t96 * qJD(3);
t144 = qJD(3) * t117 - t93 * t87;
t143 = 0.2e1 * (-t91 + t92) * qJD(3);
t142 = 2 * qJD(4);
t141 = -pkin(3) - pkin(4);
t114 = qJD(5) * t139;
t128 = qJD(5) * t93;
t28 = t114 * t94 - t128 * t96 - t144;
t95 = sin(qJ(2));
t124 = t95 * t129;
t85 = t94 * qJD(3);
t48 = pkin(3) * t85 - qJ(4) * t87 - t94 * qJD(4);
t36 = -pkin(4) * t85 - t48;
t33 = t36 - t124;
t82 = -pkin(1) * t97 - pkin(2);
t53 = t82 - t146;
t89 = t96 * pkin(4);
t42 = -t53 + t89;
t61 = t139 * t96 + t138;
t137 = t28 * t42 + t33 * t61;
t29 = qJD(3) * t61 - t114 * t96 - t128 * t94;
t62 = -t93 * t96 + t117;
t136 = t29 * t42 + t33 * t62;
t127 = pkin(2) + t146;
t54 = t89 + t127;
t135 = t28 * t54 + t36 * t61;
t134 = t29 * t54 + t36 * t62;
t37 = t124 + t48;
t133 = -t37 - t48;
t118 = pkin(1) * t95 + pkin(7);
t132 = t118 * t46;
t131 = pkin(7) * t46;
t130 = t124 * t94 + t82 * t87;
t126 = pkin(2) * t85;
t125 = pkin(2) * t87;
t122 = pkin(7) * t85;
t121 = pkin(7) * t87;
t119 = t94 * t87;
t110 = -pkin(8) + t118;
t105 = t110 * t94;
t100 = t139 * t105;
t55 = t110 * t96;
t26 = -t93 * t55 + t100;
t102 = t93 * t105;
t27 = t139 * t55 + t102;
t111 = t96 * t123;
t98 = -t110 * t85 + t111;
t112 = t94 * t123;
t99 = qJD(3) * t55 + t112;
t5 = -qJD(5) * t100 + t128 * t55 - t139 * t98 - t93 * t99;
t6 = qJD(5) * t102 + t114 * t55 - t139 * t99 + t93 * t98;
t116 = -t26 * t29 - t27 * t28 + t5 * t61 + t6 * t62;
t109 = t140 * t139;
t104 = t94 * t109;
t69 = t140 * t96;
t12 = -qJD(5) * t104 + t69 * t128 + t140 * t144;
t13 = -t109 * t87 + t69 * t114 + (-qJD(3) + qJD(5)) * t145;
t34 = -t93 * t69 + t104;
t35 = t139 * t69 + t145;
t115 = t12 * t61 + t13 * t62 - t28 * t35 - t29 * t34;
t108 = t139 * t141;
t106 = t118 * qJD(3);
t103 = -t124 * t96 + t82 * t85;
t64 = qJ(4) * t139 + t141 * t93;
t47 = -qJD(3) * t146 + t96 * qJD(4);
t74 = -0.2e1 * t119;
t73 = 0.2e1 * t119;
t63 = -qJ(4) * t93 + t108;
t52 = t127 * t85;
t45 = t93 * qJD(4) + qJD(5) * t64;
t44 = qJ(4) * t128 - qJD(4) * t139 - qJD(5) * t108;
t41 = t53 * t85;
t40 = t106 * t96 + t112;
t39 = t106 * t94 - t111;
t38 = 0.2e1 * t46;
t21 = 0.2e1 * t62 * t29;
t20 = 0.2e1 * t61 * t28;
t7 = -0.2e1 * t28 * t62 - 0.2e1 * t29 * t61;
t4 = -t139 * t29 - t93 * t28 + (-t139 * t61 + t62 * t93) * qJD(5);
t1 = -t28 * t64 - t29 * t63 + t44 * t61 + t45 * t62;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t124, -0.2e1 * t123, 0, 0, t73, t143, 0, t74, 0, 0, 0.2e1 * t103, 0.2e1 * t130, t38, 0.2e1 * t124 * t82 + 0.2e1 * t132, t73, 0, -t143, 0, 0, t74, -0.2e1 * t37 * t96 + 0.2e1 * t41, t38, -0.2e1 * t37 * t94 - 0.2e1 * t53 * t87, 0.2e1 * t37 * t53 + 0.2e1 * t132, t21, t7, 0, t20, 0, 0, 0.2e1 * t137, 0.2e1 * t136, 0.2e1 * t116, -0.2e1 * t26 * t6 - 0.2e1 * t27 * t5 + 0.2e1 * t33 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t123, 0, 0, t73, t143, 0, t74, 0, 0, t103 - t126, -t125 + t130, t46, -pkin(2) * t124 + t131, t73, 0, -t143, 0, 0, t74, t133 * t96 + t41 - t52, t46, t133 * t94 + (-t53 + t127) * t87, -t127 * t37 + t48 * t53 + t131, t21, t7, 0, t20, 0, 0, t135 + t137, t134 + t136, t115 + t116, -t12 * t27 - t13 * t26 + t33 * t54 - t34 * t6 - t35 * t5 + t36 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t143, 0, t74, 0, 0, -0.2e1 * t126, -0.2e1 * t125, 0, 0, t73, 0, -t143, 0, 0, t74, -0.2e1 * t48 * t96 - 0.2e1 * t52, 0, 0.2e1 * t127 * t87 - 0.2e1 * t48 * t94, -0.2e1 * t127 * t48, t21, t7, 0, t20, 0, 0, 0.2e1 * t135, 0.2e1 * t134, 0.2e1 * t115, -0.2e1 * t12 * t35 - 0.2e1 * t13 * t34 + 0.2e1 * t36 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, -t85, 0, -t40, t39, 0, 0, 0, t87, 0, 0, t85, 0, -t40, t47, -t39, (-pkin(3) * t123 - qJ(4) * t106) * t94 + (-pkin(3) * t106 + qJ(4) * t123 + qJD(4) * t118) * t96, 0, 0, -t29, 0, t28, 0, t6, -t5, t1, -t26 * t45 - t27 * t44 - t5 * t64 - t6 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, -t85, 0, -t121, t122, 0, 0, 0, t87, 0, 0, t85, 0, -t121, t47, -t122, t47 * pkin(7), 0, 0, -t29, 0, t28, 0, t13, -t12, t1, -t12 * t64 - t13 * t63 - t34 * t45 - t35 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, qJ(4) * t142, 0, 0, 0, 0, 0, 0, 0.2e1 * t45, -0.2e1 * t44, 0, -0.2e1 * t44 * t64 - 0.2e1 * t45 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t6 * t139 - t5 * t93 + (t139 * t27 - t26 * t93) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, t121, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t13 * t139 - t12 * t93 + (t139 * t35 - t34 * t93) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t114, 0, -t45 * t139 - t44 * t93 + (t139 * t64 - t63 * t93) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t28, 0, -t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t28, 0, -t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t114, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
