% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP6
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:28
% EndTime: 2019-12-31 19:58:32
% DurationCPUTime: 1.37s
% Computational Cost: add. (1633->163), mult. (3765->307), div. (0->0), fcn. (3352->6), ass. (0->103)
t107 = -qJ(3) - pkin(6);
t64 = sin(pkin(8));
t69 = cos(qJ(2));
t112 = t64 * t69;
t65 = cos(pkin(8));
t67 = sin(qJ(2));
t46 = t64 * t67 - t65 * t69;
t51 = t107 * t67;
t48 = t65 * t51;
t71 = (t107 * t112 + t48) * qJD(2) - t46 * qJD(3);
t47 = t65 * t67 + t112;
t59 = -t69 * pkin(2) - pkin(1);
t76 = t46 * pkin(3) - t47 * pkin(7) + t59;
t124 = -qJD(4) * t76 - t71;
t66 = sin(qJ(4));
t62 = t66 ^ 2;
t68 = cos(qJ(4));
t63 = t68 ^ 2;
t106 = t62 - t63;
t52 = t107 * t69;
t34 = t64 * t51 - t65 * t52;
t11 = -t66 * t34 + t68 * t76;
t12 = t68 * t34 + t66 * t76;
t101 = qJD(4) * t66;
t41 = t47 * qJD(2);
t100 = t67 * qJD(2);
t98 = t69 * qJD(2);
t42 = -t64 * t100 + t65 * t98;
t94 = pkin(2) * t100;
t75 = t41 * pkin(3) - t42 * pkin(7) + t94;
t95 = t124 * t68 - t66 * t75;
t3 = t34 * t101 + t95;
t61 = qJD(4) * t68;
t4 = t124 * t66 - t34 * t61 + t68 * t75;
t123 = t3 * t66 - t4 * t68 + (t11 * t66 - t12 * t68) * qJD(4);
t120 = t41 * pkin(4);
t104 = t42 * qJ(5);
t102 = qJD(4) * t47;
t90 = qJ(5) * t102;
t99 = t68 * qJD(5);
t70 = -t68 * t104 - t47 * t99 + t66 * t90 + t4;
t1 = t70 + t120;
t2 = t68 * t90 + (qJD(4) * t34 + t47 * qJD(5) + t104) * t66 + t95;
t105 = qJ(5) * t47;
t5 = t46 * pkin(4) - t68 * t105 + t11;
t8 = -t66 * t105 + t12;
t122 = -t1 * t68 + t2 * t66 + (t5 * t66 - t68 * t8) * qJD(4);
t121 = 0.2e1 * qJD(4);
t119 = t64 * pkin(2);
t118 = t65 * pkin(2);
t117 = t68 * pkin(4);
t87 = qJD(2) * t107;
t21 = t64 * (t69 * qJD(3) + t67 * t87) - t65 * (-t67 * qJD(3) + t69 * t87);
t33 = -t64 * t52 - t48;
t116 = t33 * t21;
t115 = t47 * t42;
t114 = t47 * t66;
t113 = t47 * t68;
t111 = t66 * t41;
t110 = t66 * t42;
t109 = t68 * t41;
t108 = t68 * t42;
t57 = pkin(7) + t119;
t103 = qJ(5) + t57;
t31 = 0.2e1 * t46 * t41;
t97 = -0.2e1 * pkin(1) * qJD(2);
t58 = -pkin(3) - t118;
t96 = t58 * t121;
t93 = pkin(4) * t101;
t92 = t66 * t61;
t91 = t67 * t98;
t50 = t58 - t117;
t89 = -t50 + t117;
t88 = -0.4e1 * t66 * t113;
t45 = t47 ^ 2;
t86 = t45 * t92;
t82 = pkin(4) * t62 + t50 * t68;
t81 = -t11 * t68 - t12 * t66;
t79 = -t41 * t57 + t42 * t58;
t78 = t46 * t57 - t47 * t58;
t28 = t46 * t61 + t111;
t30 = t47 * t61 + t110;
t77 = t47 * t101 - t108;
t72 = t81 * qJD(4) - t3 * t68 - t4 * t66;
t54 = -0.2e1 * t92;
t53 = 0.2e1 * t92;
t49 = -0.2e1 * t106 * qJD(4);
t44 = t103 * t68;
t43 = t103 * t66;
t36 = -t66 * qJD(5) - t103 * t61;
t35 = t103 * t101 - t99;
t26 = -t46 * t101 + t109;
t23 = (-t62 - t63) * t42;
t20 = pkin(4) * t114 + t33;
t16 = 0.2e1 * t63 * t115 - 0.2e1 * t86;
t15 = 0.2e1 * t62 * t115 + 0.2e1 * t86;
t14 = t106 * t102 - t66 * t108;
t13 = qJD(4) * t88 - t106 * t42;
t10 = t106 * t45 * t121 + t42 * t88;
t9 = t30 * pkin(4) + t21;
t7 = -0.2e1 * t47 * t111 - 0.2e1 * t30 * t46;
t6 = 0.2e1 * t47 * t109 - 0.2e1 * t77 * t46;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t91, 0.2e1 * (-t67 ^ 2 + t69 ^ 2) * qJD(2), 0, -0.2e1 * t91, 0, 0, t67 * t97, t69 * t97, 0, 0, 0.2e1 * t115, -0.2e1 * t47 * t41 - 0.2e1 * t42 * t46, 0, t31, 0, 0, 0.2e1 * t41 * t59 + 0.2e1 * t46 * t94, 0.2e1 * t42 * t59 + 0.2e1 * t47 * t94, 0.2e1 * t21 * t47 + 0.2e1 * t33 * t42 - 0.2e1 * t34 * t41 - 0.2e1 * t71 * t46, 0.2e1 * t34 * t71 + 0.2e1 * t59 * t94 + 0.2e1 * t116, t16, t10, t6, t15, t7, t31, 0.2e1 * t11 * t41 + 0.2e1 * t21 * t114 + 0.2e1 * t30 * t33 + 0.2e1 * t4 * t46, 0.2e1 * t21 * t113 - 0.2e1 * t12 * t41 + 0.2e1 * t3 * t46 - 0.2e1 * t77 * t33, 0.2e1 * t123 * t47 + 0.2e1 * t81 * t42, 0.2e1 * t11 * t4 - 0.2e1 * t12 * t3 + 0.2e1 * t116, t16, t10, t6, t15, t7, t31, 0.2e1 * t1 * t46 + 0.2e1 * t9 * t114 + 0.2e1 * t30 * t20 + 0.2e1 * t5 * t41, 0.2e1 * t9 * t113 + 0.2e1 * t2 * t46 - 0.2e1 * t77 * t20 - 0.2e1 * t8 * t41, 0.2e1 * (-t5 * t68 - t66 * t8) * t42 + 0.2e1 * t122 * t47, 0.2e1 * t1 * t5 - 0.2e1 * t2 * t8 + 0.2e1 * t20 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, -t100, 0, -pkin(6) * t98, pkin(6) * t100, 0, 0, 0, 0, t42, 0, -t41, 0, -t21, -t71, (-t41 * t64 - t42 * t65) * pkin(2), -t21 * t118 + t71 * t119, -t14, t13, t28, t14, t26, 0, -t21 * t68 + t79 * t66 + (t33 * t66 - t78 * t68) * qJD(4), t21 * t66 + t79 * t68 + (t33 * t68 + t78 * t66) * qJD(4), t72, t21 * t58 + t72 * t57, -t14, t13, t28, t14, t26, 0, t50 * t110 + t36 * t46 - t43 * t41 - t9 * t68 + (t20 * t66 + t82 * t47) * qJD(4), t50 * t108 + t35 * t46 - t44 * t41 + t9 * t66 + (t89 * t114 + t20 * t68) * qJD(4), (-t36 * t47 + t42 * t43 - t2 + (-t44 * t47 - t5) * qJD(4)) * t68 + (t35 * t47 - t42 * t44 - t1 + (-t43 * t47 - t8) * qJD(4)) * t66, -t1 * t43 - t2 * t44 + t20 * t93 - t35 * t8 + t36 * t5 + t50 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t49, 0, t54, 0, 0, t66 * t96, t68 * t96, 0, 0, t53, t49, 0, t54, 0, 0, -0.2e1 * t89 * t101, t82 * t121, -0.2e1 * t35 * t68 - 0.2e1 * t36 * t66 + 0.2e1 * (t43 * t68 - t44 * t66) * qJD(4), -0.2e1 * t35 * t44 - 0.2e1 * t43 * t36 + 0.2e1 * t50 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, t94, 0, 0, 0, 0, 0, 0, t26, -t28, t23, -t123, 0, 0, 0, 0, 0, 0, t26, -t28, t23, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * t66 + t36 * t68 + (t43 * t66 + t44 * t68) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, -t30, t41, t4, t3, 0, 0, 0, 0, -t77, 0, -t30, t41, t70 + 0.2e1 * t120, t2, t77 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t101, 0, -t57 * t61, t57 * t101, 0, 0, 0, 0, t61, 0, -t101, 0, t36, t35, -pkin(4) * t61, t36 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t61, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t61, 0, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t77, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t61, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
