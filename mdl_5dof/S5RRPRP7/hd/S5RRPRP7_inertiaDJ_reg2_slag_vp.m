% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:32
% EndTime: 2019-12-31 20:01:36
% DurationCPUTime: 1.21s
% Computational Cost: add. (1610->145), mult. (3685->280), div. (0->0), fcn. (3284->6), ass. (0->99)
t64 = sin(qJ(4));
t61 = t64 ^ 2;
t66 = cos(qJ(4));
t62 = t66 ^ 2;
t107 = t61 - t62;
t48 = t107 * qJD(4);
t63 = sin(pkin(8));
t65 = sin(qJ(2));
t105 = cos(pkin(8));
t67 = cos(qJ(2));
t89 = t105 * t67;
t45 = t63 * t65 - t89;
t46 = t105 * t65 + t63 * t67;
t57 = -t67 * pkin(2) - pkin(1);
t30 = t45 * pkin(3) - t46 * pkin(7) + t57;
t108 = -qJ(3) - pkin(6);
t49 = t108 * t65;
t50 = t108 * t67;
t34 = -t105 * t50 + t63 * t49;
t123 = t64 * t30 + t66 * t34;
t88 = qJD(2) * t108;
t39 = t67 * qJD(3) + t65 * t88;
t74 = -t65 * qJD(3) + t67 * t88;
t22 = t105 * t39 + t63 * t74;
t41 = t46 * qJD(2);
t101 = t65 * qJD(2);
t42 = qJD(2) * t89 - t63 * t101;
t95 = pkin(2) * t101;
t72 = t41 * pkin(3) - t42 * pkin(7) + t95;
t4 = -qJD(4) * t123 - t64 * t22 + t66 * t72;
t11 = t66 * t30 - t64 * t34;
t59 = qJD(4) * t64;
t60 = qJD(4) * t66;
t3 = -t66 * t22 - t30 * t60 + t34 * t59 - t64 * t72;
t122 = t3 * t64 - t4 * t66 + (t11 * t64 - t123 * t66) * qJD(4);
t103 = t45 * qJD(5);
t106 = t41 * qJ(5);
t1 = t103 - t3 + t106;
t118 = t41 * pkin(4);
t2 = -t118 - t4;
t8 = t45 * qJ(5) + t123;
t9 = -t45 * pkin(4) - t11;
t121 = t1 * t64 - t2 * t66 + (t64 * t9 + t66 * t8) * qJD(4);
t82 = t66 * pkin(4) + t64 * qJ(5);
t120 = t82 * qJD(4) - t66 * qJD(5);
t119 = 0.2e1 * qJD(5);
t21 = -t105 * t74 + t63 * t39;
t33 = -t105 * t49 - t63 * t50;
t117 = t33 * t21;
t55 = t63 * pkin(2) + pkin(7);
t116 = t41 * t55;
t115 = t45 * t55;
t114 = t46 * t42;
t113 = t46 * t64;
t112 = t46 * t66;
t111 = t64 * t41;
t110 = t66 * t41;
t109 = t66 * t42;
t102 = t64 * qJD(5);
t99 = t67 * qJD(2);
t31 = 0.2e1 * t45 * t41;
t98 = -0.2e1 * pkin(1) * qJD(2);
t97 = t64 * t109;
t56 = -t105 * pkin(2) - pkin(3);
t96 = 0.2e1 * qJD(4) * t56;
t94 = t55 * t59;
t93 = t55 * t60;
t92 = t64 * t60;
t91 = t65 * t99;
t44 = t46 ^ 2;
t87 = t44 * t92;
t83 = -t64 * t8 + t66 * t9;
t81 = pkin(4) * t64 - qJ(5) * t66;
t80 = -t11 * t66 - t123 * t64;
t77 = t42 * t56 - t116;
t76 = -t46 * t56 + t115;
t28 = t45 * t60 + t111;
t26 = t45 * t59 - t110;
t29 = t64 * t42 + t46 * t60;
t75 = t46 * t59 - t109;
t43 = -t82 + t56;
t5 = t120 * t46 + t81 * t42 + t21;
t73 = -t5 + (t43 * t46 - t115) * qJD(4);
t13 = t81 * t46 + t33;
t40 = -pkin(4) * t59 + qJ(5) * t60 + t102;
t70 = qJD(4) * t13 - t40 * t46 + t42 * t43 - t116;
t69 = t83 * qJD(4) + t1 * t66 + t2 * t64;
t68 = t80 * qJD(4) - t3 * t66 - t4 * t64;
t52 = -0.2e1 * t92;
t51 = 0.2e1 * t92;
t24 = (-t61 - t62) * t42;
t17 = 0.2e1 * t62 * t114 - 0.2e1 * t87;
t16 = 0.2e1 * t61 * t114 + 0.2e1 * t87;
t15 = t46 * t48 - t97;
t14 = t107 * t42 + 0.4e1 * t46 * t92;
t10 = t44 * t48 - 0.2e1 * t46 * t97;
t7 = t46 * t111 + t29 * t45;
t6 = 0.2e1 * t46 * t110 - 0.2e1 * t75 * t45;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t91, 0.2e1 * (-t65 ^ 2 + t67 ^ 2) * qJD(2), 0, -0.2e1 * t91, 0, 0, t65 * t98, t67 * t98, 0, 0, 0.2e1 * t114, -0.2e1 * t46 * t41 - 0.2e1 * t42 * t45, 0, t31, 0, 0, 0.2e1 * t57 * t41 + 0.2e1 * t45 * t95, 0.2e1 * t57 * t42 + 0.2e1 * t46 * t95, 0.2e1 * t21 * t46 - 0.2e1 * t22 * t45 + 0.2e1 * t33 * t42 - 0.2e1 * t34 * t41, 0.2e1 * t34 * t22 + 0.2e1 * t57 * t95 + 0.2e1 * t117, t17, 0.2e1 * t10, t6, t16, -0.2e1 * t7, t31, 0.2e1 * t11 * t41 + 0.2e1 * t21 * t113 + 0.2e1 * t29 * t33 + 0.2e1 * t4 * t45, 0.2e1 * t21 * t112 - 0.2e1 * t123 * t41 + 0.2e1 * t3 * t45 - 0.2e1 * t75 * t33, 0.2e1 * t122 * t46 + 0.2e1 * t80 * t42, 0.2e1 * t11 * t4 - 0.2e1 * t123 * t3 + 0.2e1 * t117, t17, t6, -0.2e1 * t10, t31, 0.2e1 * t7, t16, 0.2e1 * t5 * t113 + 0.2e1 * t29 * t13 - 0.2e1 * t2 * t45 - 0.2e1 * t9 * t41, -0.2e1 * t121 * t46 + 0.2e1 * t83 * t42, 0.2e1 * t1 * t45 - 0.2e1 * t5 * t112 + 0.2e1 * t75 * t13 + 0.2e1 * t8 * t41, 0.2e1 * t8 * t1 + 0.2e1 * t13 * t5 + 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, -t101, 0, -pkin(6) * t99, pkin(6) * t101, 0, 0, 0, 0, t42, 0, -t41, 0, -t21, -t22, (-t105 * t42 - t41 * t63) * pkin(2), (-t105 * t21 + t22 * t63) * pkin(2), -t15, -t14, t28, t15, -t26, 0, -t21 * t66 + t77 * t64 + (t33 * t64 - t76 * t66) * qJD(4), t21 * t64 + t77 * t66 + (t33 * t66 + t76 * t64) * qJD(4), t68, t21 * t56 + t68 * t55, -t15, t28, t14, 0, t26, t15, t70 * t64 + t73 * t66, t69, t73 * t64 - t70 * t66, -t13 * t40 + t5 * t43 + t69 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -0.2e1 * t48, 0, t52, 0, 0, t64 * t96, t66 * t96, 0, 0, t51, 0, 0.2e1 * t48, 0, 0, t52, 0.2e1 * t40 * t66 + 0.2e1 * t43 * t59, 0, 0.2e1 * t40 * t64 - 0.2e1 * t43 * t60, -0.2e1 * t43 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, t95, 0, 0, 0, 0, 0, 0, -t26, -t28, t24, -t122, 0, 0, 0, 0, 0, 0, -t26, t24, t28, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, 0, -t29, t41, t4, t3, 0, 0, 0, -t75, 0, t41, t29, 0, t4 + 0.2e1 * t118, -t82 * t42 + (t81 * qJD(4) - t102) * t46, 0.2e1 * t103 - t3 + 0.2e1 * t106, -t2 * pkin(4) + t1 * qJ(5) + t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t59, 0, -t93, t94, 0, 0, 0, t60, 0, 0, t59, 0, -t93, -t120, -t94, -t120 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t60, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, t60, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, qJ(5) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t75, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
