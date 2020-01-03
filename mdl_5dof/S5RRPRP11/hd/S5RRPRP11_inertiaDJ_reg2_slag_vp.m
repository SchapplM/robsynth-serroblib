% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:04
% EndTime: 2019-12-31 20:14:09
% DurationCPUTime: 1.26s
% Computational Cost: add. (828->145), mult. (1858->261), div. (0->0), fcn. (1290->4), ass. (0->103)
t58 = sin(qJ(2));
t103 = t58 * qJ(3);
t114 = pkin(2) + pkin(7);
t60 = cos(qJ(2));
t107 = t114 * t60;
t118 = t103 + t107;
t104 = qJ(3) * t60;
t108 = t114 * t58;
t66 = -t104 + t108;
t113 = pkin(3) + pkin(6);
t31 = -pkin(1) - t118;
t41 = t113 * t58;
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t117 = t59 * t31 + t57 * t41;
t53 = t57 ^ 2;
t55 = t59 ^ 2;
t106 = t53 - t55;
t56 = t60 ^ 2;
t79 = qJD(2) * (t58 ^ 2 - t56);
t38 = t106 * qJD(4);
t48 = t60 * qJD(2);
t76 = t113 * t48;
t96 = t58 * qJD(3);
t62 = -qJD(4) * t117 + t57 * t96 + t59 * t76;
t116 = 0.2e1 * qJD(3);
t115 = 2 * qJD(5);
t112 = t60 * pkin(4);
t70 = -t57 * pkin(4) + t59 * qJ(5);
t23 = t70 * qJD(4) + t57 * qJD(5);
t71 = pkin(4) * t59 + qJ(5) * t57;
t97 = t58 * qJD(2);
t7 = t23 * t60 + (-t71 - t113) * t97;
t111 = t7 * t57;
t21 = t71 * qJD(4) - t59 * qJD(5) + qJD(3);
t110 = t21 * t60;
t35 = t113 * t97;
t109 = t35 * t57;
t42 = t113 * t60;
t102 = qJD(2) * t57;
t101 = qJD(2) * t59;
t100 = qJD(4) * t57;
t49 = qJD(4) * t59;
t99 = qJD(4) * t60;
t98 = qJD(4) * t114;
t95 = t58 * qJD(5);
t94 = t60 * qJD(3);
t93 = qJ(3) * qJD(4);
t92 = qJ(5) * qJD(2);
t91 = -0.2e1 * pkin(1) * qJD(2);
t90 = pkin(6) * t97;
t89 = pkin(6) * t48;
t88 = t57 * t99;
t87 = t57 * t98;
t86 = t59 * t99;
t85 = t59 * t98;
t84 = t58 * t48;
t83 = t59 * t97;
t82 = t59 * t48;
t81 = t57 * t49;
t80 = t60 * t92;
t78 = t57 * t83;
t77 = t56 * t81;
t8 = t58 * qJ(5) + t117;
t10 = -t57 * t31 + t59 * t41;
t9 = -t58 * pkin(4) - t10;
t74 = t57 * t9 + t59 * t8;
t73 = -t60 * pkin(2) - t103;
t72 = pkin(2) * t58 - t104;
t69 = -t10 * t57 + t117 * t59;
t36 = qJ(3) - t70;
t67 = t36 * t60 - t108;
t65 = t66 * qJD(2);
t5 = t31 * t100 - t59 * (t65 - t96) - t57 * t76 - t41 * t49;
t64 = t57 * (pkin(7) * t58 + t72);
t63 = t73 * qJD(2) + t94;
t3 = -t5 + t80 + t95;
t4 = (t64 - t112) * qJD(2) - t62;
t1 = t74 * qJD(4) + t3 * t57 - t4 * t59;
t6 = -t57 * t65 + t62;
t2 = t69 * qJD(4) - t5 * t57 + t6 * t59;
t52 = qJ(3) * t116;
t46 = -0.2e1 * t81;
t45 = 0.2e1 * t81;
t44 = -0.2e1 * t84;
t43 = 0.2e1 * t84;
t40 = t114 * t82;
t37 = -pkin(1) + t73;
t33 = -0.2e1 * t79;
t28 = t57 * t97 - t86;
t27 = t57 * t48 + t58 * t49;
t26 = t83 + t88;
t25 = -t58 * t100 + t82;
t24 = t72 * qJD(2) - t96;
t20 = t71 * t60 + t42;
t19 = -0.2e1 * t55 * t84 - 0.2e1 * t77;
t18 = -0.2e1 * t53 * t84 + 0.2e1 * t77;
t17 = -t106 * t99 - t78;
t16 = t58 * t88 + t59 * t79;
t15 = -t106 * t97 + 0.4e1 * t60 * t81;
t13 = 0.2e1 * t57 * t79 - 0.2e1 * t58 * t86;
t12 = -t56 * t38 - 0.2e1 * t60 * t78;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t33, 0, t44, 0, 0, t58 * t91, t60 * t91, 0, 0, 0, 0, 0, t43, t33, t44, 0, 0.2e1 * t24 * t60 - 0.2e1 * t37 * t97, -0.2e1 * t24 * t58 - 0.2e1 * t37 * t48, 0.2e1 * t37 * t24, t18, 0.2e1 * t12, t13, t19, 0.2e1 * t16, t43, 0.2e1 * (-t42 * t101 + t6) * t58 + 0.2e1 * (qJD(2) * t10 - t42 * t100 - t35 * t59) * t60, 0.2e1 * (t42 * t102 + t5) * t58 + 0.2e1 * (-qJD(2) * t117 - t42 * t49 + t109) * t60, 0.2e1 * t69 * t97 + 0.2e1 * (t5 * t59 + t57 * t6 + (t10 * t59 + t117 * t57) * qJD(4)) * t60, 0.2e1 * t10 * t6 - 0.2e1 * t117 * t5 - 0.2e1 * t42 * t35, t18, t13, -0.2e1 * t12, t43, -0.2e1 * t16, t19, 0.2e1 * (-t20 * t101 - t4) * t58 + 0.2e1 * (-qJD(2) * t9 - t20 * t100 + t7 * t59) * t60, 0.2e1 * t74 * t97 + 0.2e1 * (-t3 * t59 - t4 * t57 + (t57 * t8 - t59 * t9) * qJD(4)) * t60, 0.2e1 * (-t20 * t102 + t3) * t58 + 0.2e1 * (qJD(2) * t8 + t20 * t49 + t111) * t60, 0.2e1 * t20 * t7 + 0.2e1 * t8 * t3 + 0.2e1 * t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t97, 0, -t89, t90, 0, 0, 0, -t48, t97, 0, 0, 0, t63, t89, -t90, t63 * pkin(6), -t17, t15, t25, t17, -t27, 0, -t109 - t40 + (-qJ(3) * t97 + t94) * t59 + (t42 * t59 + t66 * t57) * qJD(4), (t66 * qJD(4) - t35) * t59 + (t118 * qJD(2) - qJD(4) * t42 - t94) * t57, -t2, -t35 * qJ(3) + t42 * qJD(3) - t114 * t2, -t17, t25, -t15, 0, t27, t17, t111 - t40 + (-t36 * t97 + t110) * t59 + (t20 * t59 - t67 * t57) * qJD(4), -t1, (t67 * qJD(4) - t7) * t59 + (qJD(4) * t20 + t110 + (-t36 * t58 - t107) * qJD(2)) * t57, -t1 * t114 + t20 * t21 + t7 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t52, t46, 0.2e1 * t38, 0, t45, 0, 0, 0.2e1 * qJD(3) * t57 + 0.2e1 * t59 * t93, 0.2e1 * qJD(3) * t59 - 0.2e1 * t57 * t93, 0, t52, t46, 0, -0.2e1 * t38, 0, 0, t45, 0.2e1 * t21 * t57 + 0.2e1 * t36 * t49, 0, 0.2e1 * t36 * t100 - 0.2e1 * t21 * t59, 0.2e1 * t36 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, t89, 0, 0, 0, 0, 0, 0, t25, -t27, 0, t2, 0, 0, 0, 0, 0, 0, t25, 0, t27, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t26, t48, t6, t5, 0, 0, 0, t28, 0, t48, -t26, 0, (-t64 + 0.2e1 * t112) * qJD(2) + t62, (-pkin(4) * t97 + qJ(5) * t99) * t57 + (t58 * t92 + (pkin(4) * qJD(4) - qJD(5)) * t60) * t59, -t5 + 0.2e1 * t80 + 0.2e1 * t95, -t4 * pkin(4) + t3 * qJ(5) + t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, 0, -t49, 0, t87, t85, 0, 0, 0, -t100, 0, 0, t49, 0, t87, -t23, -t85, -t23 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t100, 0, t49, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, qJ(5) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t28, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, 0, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
