% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:20
% DurationCPUTime: 1.06s
% Computational Cost: add. (777->150), mult. (1882->268), div. (0->0), fcn. (1369->6), ass. (0->109)
t57 = sin(qJ(3));
t121 = -0.4e1 * t57;
t58 = cos(qJ(4));
t51 = qJD(4) * t58;
t56 = sin(qJ(4));
t59 = cos(qJ(3));
t94 = t59 * qJD(3);
t27 = t57 * t51 + t56 * t94;
t120 = t27 * pkin(4);
t52 = t56 ^ 2;
t54 = t58 ^ 2;
t104 = t52 - t54;
t46 = sin(pkin(8)) * pkin(1) + pkin(6);
t119 = t57 * qJD(5) + (qJ(5) * qJD(3) + qJD(4) * t46) * t59;
t118 = pkin(3) * t57;
t117 = pkin(7) * t59;
t116 = t58 * pkin(4);
t78 = t46 * t94;
t11 = t78 + t120;
t115 = t11 * t56;
t114 = t11 * t58;
t107 = -qJ(5) - pkin(7);
t37 = t107 * t56;
t113 = t37 * t57;
t38 = t107 * t58;
t112 = t38 * t57;
t111 = t46 * t56;
t110 = t56 * t59;
t109 = t57 * t58;
t108 = t58 * t59;
t47 = -cos(pkin(8)) * pkin(1) - pkin(2);
t73 = -t59 * pkin(3) - t57 * pkin(7);
t64 = t47 + t73;
t63 = qJD(4) * t64;
t72 = -t117 + t118;
t65 = t72 * qJD(3);
t106 = -t56 * t65 - t58 * t63;
t49 = t57 * qJD(3);
t80 = t46 * t49;
t105 = t56 * t80 + t58 * t65;
t33 = t46 * t108;
t8 = t56 * t64 + t33;
t53 = t57 ^ 2;
t103 = -t59 ^ 2 + t53;
t102 = qJ(5) * t57;
t23 = (pkin(4) * t56 + t46) * t57;
t101 = qJD(3) * t23;
t100 = qJD(3) * t58;
t99 = qJD(4) * t53;
t98 = qJD(4) * t56;
t97 = qJD(4) * t57;
t96 = qJD(4) * t59;
t93 = t46 * t110;
t92 = 0.2e1 * qJD(3) * t47;
t91 = -0.2e1 * t98;
t90 = pkin(4) * t49;
t89 = pkin(4) * t98;
t88 = t46 * t99;
t87 = t56 * t97;
t86 = t56 * t96;
t85 = t58 * t96;
t84 = t23 * t98;
t83 = t52 * t94;
t82 = t56 * t51;
t81 = t57 * t94;
t79 = t58 * t94;
t48 = -pkin(3) - t116;
t77 = -t48 + t116;
t76 = t103 * qJD(3);
t75 = t56 * t79;
t74 = t53 * t82;
t20 = t58 * t64;
t5 = -t58 * t102 + t20 + (-pkin(4) - t111) * t59;
t6 = -t56 * t102 + t8;
t71 = -t5 * t58 - t56 * t6;
t70 = t5 * t56 - t58 * t6;
t7 = t20 - t93;
t69 = -t56 * t8 - t58 * t7;
t68 = t56 * t7 - t58 * t8;
t67 = pkin(4) * t52 + t48 * t58;
t25 = -t79 + t87;
t26 = t58 * t49 + t86;
t3 = t26 * t46 + t106;
t4 = -t8 * qJD(4) + t105;
t62 = t69 * qJD(4) - t3 * t58 - t4 * t56;
t21 = -t58 * qJD(5) - t107 * t98;
t22 = -t56 * qJD(5) + t107 * t51;
t61 = -t21 * t58 - t22 * t56 + (-t37 * t58 + t38 * t56) * qJD(4);
t60 = qJ(5) * t87 - t119 * t58 - t56 * t63 + t105;
t44 = t54 * t94;
t42 = -0.2e1 * t81;
t41 = -0.2e1 * t82;
t40 = 0.2e1 * t82;
t35 = t54 * t81;
t34 = t52 * t81;
t32 = -0.2e1 * t104 * qJD(4);
t28 = t56 * t49 - t85;
t24 = t44 + t83;
t17 = 0.2e1 * t35 - 0.2e1 * t74;
t16 = 0.2e1 * t34 + 0.2e1 * t74;
t15 = t104 * t97 - t75;
t14 = t82 * t121 + t44 - t83;
t13 = -0.2e1 * t56 * t76 + 0.2e1 * t57 * t85;
t12 = 0.2e1 * t103 * t100 + 0.2e1 * t57 * t86;
t10 = 0.2e1 * t104 * t99 + t75 * t121;
t9 = 0.2e1 * t34 + 0.2e1 * t35 - 0.2e1 * t81;
t2 = (qJ(5) * qJD(4) + qJD(3) * t46) * t109 + t119 * t56 + t106;
t1 = t60 + t90;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t81, -0.2e1 * t76, 0, t42, 0, 0, t57 * t92, t59 * t92, 0, 0, t17, t10, t12, t16, t13, t42, 0.2e1 * t58 * t88 - 0.2e1 * t4 * t59 + 0.2e1 * (t7 + 0.2e1 * t93) * t49, -0.2e1 * t56 * t88 - 0.2e1 * t3 * t59 + 0.2e1 * (-t8 + 0.2e1 * t33) * t49, 0.2e1 * t69 * t94 + 0.2e1 * (t68 * qJD(4) + t3 * t56 - t4 * t58) * t57, 0.2e1 * t46 ^ 2 * t81 - 0.2e1 * t8 * t3 + 0.2e1 * t7 * t4, t17, t10, t12, t16, t13, t42, 0.2e1 * (t56 * t101 - t1) * t59 + 0.2e1 * (qJD(3) * t5 + t23 * t51 + t115) * t57, 0.2e1 * (t23 * t100 - t2) * t59 + 0.2e1 * (-qJD(3) * t6 + t114 - t84) * t57, 0.2e1 * t71 * t94 + 0.2e1 * (t70 * qJD(4) - t1 * t58 + t2 * t56) * t57, 0.2e1 * t5 * t1 + 0.2e1 * t23 * t11 - 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t57 + (t103 * t46 - t68 * t59) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t70 * qJD(3) - t11) * t59 + (t71 * qJD(4) - t1 * t56 - t2 * t58 + t101) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, -t49, 0, -t78, t80, 0, 0, -t15, t14, t28, t15, t26, 0, (pkin(7) * t108 + (-pkin(3) * t58 + t111) * t57) * qJD(4) + (t73 * t56 - t33) * qJD(3), (t46 * t109 + t72 * t56) * qJD(4) + (t73 * t58 + t93) * qJD(3), t62, -pkin(3) * t78 + t62 * pkin(7), -t15, t14, t28, t15, t26, 0, -t114 - t22 * t59 + (t48 * t110 + t113) * qJD(3) + (t23 * t56 + t67 * t57) * qJD(4), t115 - t21 * t59 + (t48 * t108 + t112) * qJD(3) + (t77 * t57 * t56 + t23 * t58) * qJD(4), (-t37 * t94 - t22 * t57 - t2 + (-t5 + t112) * qJD(4)) * t58 + (t38 * t94 + t21 * t57 - t1 + (-t6 + t113) * qJD(4)) * t56, pkin(4) * t84 + t1 * t37 + t11 * t48 + t2 * t38 - t6 * t21 + t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t94, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t28, t24, (-t118 + (t52 + t54) * t117) * qJD(3), 0, 0, 0, 0, 0, 0, -t26, t28, t24, (-t89 + (-t37 * t56 - t38 * t58) * qJD(3)) * t59 + (qJD(3) * t48 + t61) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t32, 0, t41, 0, 0, pkin(3) * t91, -0.2e1 * pkin(3) * t51, 0, 0, t40, t32, 0, t41, 0, 0, t77 * t91, 0.2e1 * t67 * qJD(4), 0.2e1 * t61, 0.2e1 * t38 * t21 + 0.2e1 * t37 * t22 + 0.2e1 * t48 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, -t27, t49, t4, t3, 0, 0, 0, 0, -t25, 0, -t27, t49, t60 + 0.2e1 * t90, t2, t25 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t25, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t25, 0, -t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t98, 0, -pkin(7) * t51, pkin(7) * t98, 0, 0, 0, 0, t51, 0, -t98, 0, t22, t21, -pkin(4) * t51, t22 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t51, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
