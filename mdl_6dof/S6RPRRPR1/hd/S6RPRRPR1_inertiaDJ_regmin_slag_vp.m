% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:44
% EndTime: 2019-03-09 04:58:48
% DurationCPUTime: 1.12s
% Computational Cost: add. (2170->145), mult. (4684->261), div. (0->0), fcn. (4505->10), ass. (0->103)
t131 = qJD(3) + qJD(4);
t67 = cos(qJ(3));
t122 = sin(qJ(4));
t100 = t122 * t67;
t54 = sin(pkin(10)) * pkin(1) + pkin(7);
t125 = pkin(8) + t54;
t123 = cos(qJ(4));
t127 = t123 * t125;
t65 = sin(qJ(3));
t81 = t65 * t127;
t72 = t125 * t100 + t81;
t90 = t122 * t125;
t22 = -(-t67 * t90 - t81) * qJD(3) + t72 * qJD(4);
t23 = t131 * (-t127 * t67 + t65 * t90);
t102 = t123 * t67;
t96 = t122 * qJD(4);
t97 = qJD(4) * t123;
t34 = -qJD(3) * t102 - t67 * t97 + (t122 * qJD(3) + t96) * t65;
t46 = t123 * t65 + t100;
t129 = t34 * qJ(5) - t46 * qJD(5) + t23;
t66 = cos(qJ(6));
t61 = t66 ^ 2;
t64 = sin(qJ(6));
t113 = t64 ^ 2 - t61;
t95 = t113 * qJD(6);
t111 = cos(pkin(11));
t62 = sin(pkin(11));
t75 = t122 * t65 - t102;
t33 = t111 * t46 - t62 * t75;
t68 = t131 * t46;
t10 = -t68 * qJ(5) - t75 * qJD(5) - t22;
t4 = t10 * t62 - t111 * t129;
t126 = t33 * t4;
t28 = (t123 * qJ(5) + t127) * t67 + (-t122 * qJ(5) - t90) * t65;
t70 = -t46 * qJ(5) - t72;
t17 = -t111 * t70 + t28 * t62;
t59 = qJD(6) * t66;
t124 = t17 * t59 + t4 * t64;
t26 = -t111 * t34 - t62 * t68;
t121 = t26 * t33;
t120 = t26 * t66;
t25 = t111 * t68 - t34 * t62;
t32 = t111 * t75 + t46 * t62;
t119 = t32 * t25;
t118 = t33 * t64;
t117 = t64 * t66;
t116 = t66 * t25;
t115 = t33 * t116 + t32 * t120;
t101 = t122 * t62;
t57 = t123 * pkin(3) + pkin(4);
t42 = -pkin(3) * t101 + t111 * t57;
t38 = -pkin(5) - t42;
t112 = pkin(3) * qJD(4);
t89 = t111 * t122;
t40 = (t123 * t62 + t89) * t112;
t114 = t38 * t59 + t40 * t64;
t43 = pkin(3) * t89 + t62 * t57;
t110 = qJD(6) * t64;
t109 = t65 * qJD(3);
t108 = t67 * qJD(3);
t107 = 0.2e1 * t108;
t58 = pkin(3) * t109;
t106 = t33 * t110;
t55 = -t111 * pkin(4) - pkin(5);
t105 = t55 * t110;
t104 = t55 * t59;
t103 = t64 * t59;
t56 = -cos(pkin(10)) * pkin(1) - pkin(2);
t99 = -0.4e1 * t33 * t117;
t98 = t38 * t110 - t40 * t66;
t94 = pkin(3) * t97;
t93 = pkin(3) * t96;
t18 = t111 * t28 + t62 * t70;
t47 = -t67 * pkin(3) + t56;
t69 = t75 * pkin(4) + t47;
t21 = t32 * pkin(5) - t33 * pkin(9) + t69;
t88 = t18 * t66 + t21 * t64;
t87 = t18 * t64 - t21 * t66;
t53 = pkin(4) * t62 + pkin(9);
t86 = -t25 * t53 + t26 * t55;
t39 = pkin(9) + t43;
t85 = t32 * t39 - t33 * t38;
t41 = (t111 * t123 - t101) * t112;
t84 = -t32 * t41 + t33 * t40;
t83 = t32 * t53 - t33 * t55;
t15 = t25 * t64 + t32 * t59;
t76 = t26 * t64 + t33 * t59;
t14 = t106 - t120;
t74 = t47 * t46;
t73 = -t25 * t39 + t26 * t38 + t84;
t29 = t68 * pkin(4) + t58;
t49 = 0.2e1 * t103;
t45 = -0.2e1 * t95;
t30 = t33 ^ 2;
t13 = t32 * t110 - t116;
t11 = t17 * t110;
t8 = t26 * t117 - t33 * t95;
t7 = t25 * pkin(5) - t26 * pkin(9) + t29;
t6 = qJD(6) * t99 - t113 * t26;
t5 = t111 * t10 + t129 * t62;
t2 = -t88 * qJD(6) - t64 * t5 + t66 * t7;
t1 = t87 * qJD(6) - t66 * t5 - t64 * t7;
t3 = [0, 0, 0, 0, t65 * t107, 0.2e1 * (-t65 ^ 2 + t67 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t56 * t109, t56 * t107, -0.2e1 * t46 * t34, 0.2e1 * t34 * t75 - 0.2e1 * t46 * t68, 0, 0, 0, 0.2e1 * qJD(4) * t74 + 0.2e1 * (t65 * pkin(3) * t75 + t74) * qJD(3), -0.2e1 * t34 * t47 + 0.2e1 * t46 * t58, 0.2e1 * t17 * t26 - 0.2e1 * t18 * t25 - 0.2e1 * t32 * t5 + 0.2e1 * t126, 0.2e1 * t17 * t4 + 0.2e1 * t18 * t5 + 0.2e1 * t69 * t29, -0.2e1 * t30 * t103 + 0.2e1 * t61 * t121, t26 * t99 + 0.2e1 * t30 * t95, -0.2e1 * t32 * t106 + 0.2e1 * t115, -0.2e1 * t25 * t118 - 0.2e1 * t76 * t32, 0.2e1 * t119, 0.2e1 * t4 * t118 + 0.2e1 * t76 * t17 + 0.2e1 * t2 * t32 - 0.2e1 * t87 * t25, 0.2e1 * t1 * t32 + 0.2e1 * t66 * t126 - 0.2e1 * t14 * t17 - 0.2e1 * t88 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t25 + t18 * t26 + t32 * t4 + t33 * t5, 0, 0, 0, 0, 0, 0 (-t25 * t33 - t26 * t32) * t66 + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t119 + 0.2e1 * t121, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t108, -t109, 0, -t54 * t108, t54 * t109, 0, 0, -t34, -t68, 0, t23, t22, -t25 * t43 - t26 * t42 + t84, t17 * t40 + t18 * t41 - t4 * t42 + t43 * t5, t8, t6, t15, -t13, 0, t11 + (-t85 * qJD(6) - t4) * t66 + t73 * t64, t85 * t110 + t73 * t66 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t108, 0, 0, 0, 0, 0, -t68, t34, 0, -t25 * t42 + t26 * t43 + t32 * t40 + t33 * t41, 0, 0, 0, 0, 0, t13, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t93, -0.2e1 * t94, 0, -0.2e1 * t40 * t42 + 0.2e1 * t41 * t43, t49, t45, 0, 0, 0, 0.2e1 * t98, 0.2e1 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t68, 0, t23, t22 (-t111 * t26 - t25 * t62) * pkin(4) (-t111 * t4 + t5 * t62) * pkin(4), t8, t6, t15, -t13, 0, t11 + t86 * t64 + (-t83 * qJD(6) - t4) * t66, t83 * t110 + t86 * t66 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t34, 0 (-t111 * t25 + t26 * t62) * pkin(4), 0, 0, 0, 0, 0, t13, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t94, 0 (-t111 * t40 + t41 * t62) * pkin(4), t49, t45, 0, 0, 0, t98 + t105, t104 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t45, 0, 0, 0, 0.2e1 * t105, 0.2e1 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, -t13, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t76, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t110, 0, -t39 * t59 - t41 * t64, t39 * t110 - t41 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t110, 0, -t53 * t59, t53 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
