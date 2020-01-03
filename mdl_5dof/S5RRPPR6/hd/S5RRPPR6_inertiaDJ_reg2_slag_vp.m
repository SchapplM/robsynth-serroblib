% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:05
% EndTime: 2019-12-31 19:33:10
% DurationCPUTime: 1.20s
% Computational Cost: add. (2148->155), mult. (4906->312), div. (0->0), fcn. (4672->8), ass. (0->104)
t110 = cos(pkin(8));
t65 = sin(pkin(8));
t68 = sin(qJ(2));
t69 = cos(qJ(2));
t48 = t110 * t68 + t65 * t69;
t41 = t48 * qJD(2);
t108 = t68 * qJD(2);
t95 = t110 * t69;
t42 = qJD(2) * t95 - t65 * t108;
t46 = t65 * t68 - t95;
t87 = t41 * t48 + t42 * t46;
t128 = -0.2e1 * t87;
t67 = sin(qJ(5));
t109 = qJD(5) * t67;
t64 = sin(pkin(9));
t66 = cos(pkin(9));
t126 = cos(qJ(5));
t97 = qJD(5) * t126;
t43 = t64 * t109 - t66 * t97;
t47 = -t126 * t66 + t64 * t67;
t57 = pkin(2) * t65 + qJ(4);
t127 = pkin(7) + t57;
t114 = -qJ(3) - pkin(6);
t94 = qJD(2) * t114;
t40 = qJD(3) * t69 + t68 * t94;
t82 = -t68 * qJD(3) + t69 * t94;
t22 = -t110 * t82 + t40 * t65;
t52 = t114 * t68;
t53 = t114 * t69;
t32 = -t110 * t52 - t53 * t65;
t125 = t32 * t22;
t116 = t67 * t66;
t99 = t126 * t64;
t49 = t99 + t116;
t44 = t49 * qJD(5);
t124 = t47 * t44;
t123 = t49 * t43;
t33 = -t110 * t53 + t65 * t52;
t122 = t64 * t33;
t121 = t64 * t41;
t120 = t64 * t42;
t119 = t64 * t48;
t117 = t66 * t42;
t115 = -pkin(7) - qJ(4);
t13 = t49 * t42 - t43 * t48;
t23 = t49 * t48;
t113 = -t49 * t13 + t43 * t23;
t61 = -t69 * pkin(2) - pkin(1);
t85 = t46 * pkin(3) + t61;
t81 = -t48 * qJ(4) + t85;
t15 = t66 * t33 + t64 * t81;
t62 = t64 ^ 2;
t63 = t66 ^ 2;
t111 = t62 + t63;
t107 = t69 * qJD(2);
t28 = 0.2e1 * t46 * t41;
t106 = 0.2e1 * t48 * t42;
t105 = -0.2e1 * pkin(1) * qJD(2);
t104 = t64 * t117;
t103 = pkin(2) * t108;
t101 = t68 * t107;
t100 = t67 * t127;
t96 = t126 * qJD(4);
t93 = 0.2e1 * t111 * qJD(4);
t77 = t110 * t40 + t65 * t82;
t76 = t64 * t77;
t80 = t41 * pkin(3) - t48 * qJD(4) + t103;
t79 = -qJ(4) * t42 + t80;
t8 = t66 * t79 - t76;
t75 = t66 * t77;
t9 = t64 * t79 + t75;
t91 = t64 * t9 + t66 * t8;
t90 = -t64 * t8 + t66 * t9;
t60 = -t110 * pkin(2) - pkin(3);
t12 = t47 * t42 + t48 * t44;
t24 = t47 * t48;
t89 = -t12 * t47 - t24 * t44;
t88 = t22 * t48 + t32 * t42;
t86 = t41 * t49 - t43 * t46;
t84 = t127 * t99;
t83 = -qJD(4) * t46 - t41 * t57 + t42 * t60;
t78 = t115 * t42 + t80;
t74 = t46 * pkin(4) - t122 + (t115 * t48 + t85) * t66;
t73 = t67 * t74;
t72 = t126 * t74;
t71 = t78 * t64 + t75;
t70 = t41 * pkin(4) + t78 * t66 - t76;
t51 = -t66 * pkin(4) + t60;
t45 = t127 * t66;
t38 = t66 * t41;
t27 = -t64 * t100 + t126 * t45;
t26 = -t67 * t45 - t84;
t21 = pkin(4) * t119 + t32;
t20 = -t45 * t97 - qJD(4) * t116 + (qJD(5) * t100 - t96) * t64;
t19 = qJD(5) * t84 - t66 * t96 + (t64 * qJD(4) + qJD(5) * t45) * t67;
t17 = -t41 * t47 - t44 * t46;
t16 = pkin(4) * t120 + t22;
t14 = t66 * t81 - t122;
t11 = -pkin(7) * t119 + t15;
t4 = t126 * t11 + t73;
t3 = -t67 * t11 + t72;
t2 = -qJD(5) * t73 - t11 * t97 + t126 * t70 - t67 * t71;
t1 = -qJD(5) * t72 + t11 * t109 - t126 * t71 - t67 * t70;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101, 0.2e1 * (-t68 ^ 2 + t69 ^ 2) * qJD(2), 0, -0.2e1 * t101, 0, 0, t68 * t105, t69 * t105, 0, 0, t106, t128, 0, t28, 0, 0, 0.2e1 * t46 * t103 + 0.2e1 * t41 * t61, 0.2e1 * t48 * t103 + 0.2e1 * t42 * t61, -0.2e1 * t33 * t41 - 0.2e1 * t77 * t46 + 0.2e1 * t88, 0.2e1 * t61 * t103 + 0.2e1 * t33 * t77 + 0.2e1 * t125, t63 * t106, -0.4e1 * t48 * t104, 0.2e1 * t87 * t66, t62 * t106, t64 * t128, t28, 0.2e1 * t14 * t41 + 0.2e1 * t46 * t8 + 0.2e1 * t88 * t64, -0.2e1 * t15 * t41 - 0.2e1 * t46 * t9 + 0.2e1 * t88 * t66, -0.2e1 * t91 * t48 + 0.2e1 * (-t14 * t66 - t15 * t64) * t42, 0.2e1 * t14 * t8 + 0.2e1 * t15 * t9 + 0.2e1 * t125, 0.2e1 * t24 * t12, 0.2e1 * t12 * t23 + 0.2e1 * t13 * t24, -0.2e1 * t12 * t46 - 0.2e1 * t24 * t41, 0.2e1 * t23 * t13, -0.2e1 * t13 * t46 - 0.2e1 * t23 * t41, t28, 0.2e1 * t13 * t21 + 0.2e1 * t16 * t23 + 0.2e1 * t2 * t46 + 0.2e1 * t3 * t41, 0.2e1 * t1 * t46 - 0.2e1 * t12 * t21 - 0.2e1 * t16 * t24 - 0.2e1 * t4 * t41, 0.2e1 * t1 * t23 + 0.2e1 * t12 * t3 - 0.2e1 * t13 * t4 + 0.2e1 * t2 * t24, -0.2e1 * t1 * t4 + 0.2e1 * t16 * t21 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, -t108, 0, -pkin(6) * t107, pkin(6) * t108, 0, 0, 0, 0, t42, 0, -t41, 0, -t22, -t77, (-t110 * t42 - t41 * t65) * pkin(2), (-t22 * t110 + t77 * t65) * pkin(2), t104, (-t62 + t63) * t42, t121, -t104, t38, 0, -t22 * t66 + t83 * t64, t22 * t64 + t83 * t66, t90, t22 * t60 + t90 * t57 + (-t14 * t64 + t15 * t66) * qJD(4), -t12 * t49 + t24 * t43, -t89 + t113, t86, t13 * t47 + t23 * t44, t17, 0, t13 * t51 + t16 * t47 + t20 * t46 + t21 * t44 + t26 * t41, -t12 * t51 + t16 * t49 + t19 * t46 - t21 * t43 - t27 * t41, t1 * t47 + t12 * t26 - t13 * t27 + t19 * t23 - t2 * t49 + t20 * t24 + t3 * t43 - t4 * t44, -t1 * t27 + t16 * t51 - t19 * t4 + t2 * t26 + t20 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t57 * t93, -0.2e1 * t123, 0.2e1 * t47 * t43 - 0.2e1 * t49 * t44, 0, 0.2e1 * t124, 0, 0, 0.2e1 * t51 * t44, -0.2e1 * t51 * t43, 0.2e1 * t19 * t47 - 0.2e1 * t20 * t49 + 0.2e1 * t26 * t43 - 0.2e1 * t27 * t44, -0.2e1 * t19 * t27 + 0.2e1 * t20 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, t103, 0, 0, 0, 0, 0, 0, t38, -t121, -t111 * t42, t91, 0, 0, 0, 0, 0, 0, t17, -t86, t89 + t113, -t1 * t49 - t2 * t47 - t3 * t44 - t4 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t49 - t20 * t47 - t26 * t44 - t27 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t123 + 0.2e1 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t117, 0, t22, 0, 0, 0, 0, 0, 0, t13, -t12, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t13, t41, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, -t44, 0, t20, t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
