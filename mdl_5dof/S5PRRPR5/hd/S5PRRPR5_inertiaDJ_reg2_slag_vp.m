% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:28:01
% DurationCPUTime: 1.17s
% Computational Cost: add. (1256->148), mult. (3280->300), div. (0->0), fcn. (3187->10), ass. (0->97)
t48 = sin(qJ(5));
t43 = t48 ^ 2;
t51 = cos(qJ(5));
t44 = t51 ^ 2;
t100 = t43 - t44;
t79 = qJD(5) * t100;
t46 = sin(pkin(5));
t50 = sin(qJ(2));
t106 = t46 * t50;
t47 = cos(pkin(5));
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t29 = -t49 * t106 + t47 * t52;
t53 = cos(qJ(2));
t97 = qJD(2) * t53;
t86 = t46 * t97;
t20 = t29 * qJD(3) + t52 * t86;
t30 = t52 * t106 + t47 * t49;
t21 = -t30 * qJD(3) - t49 * t86;
t45 = sin(pkin(10));
t99 = cos(pkin(10));
t62 = t99 * t20 + t45 * t21;
t105 = t46 * t53;
t17 = t45 * t29 + t99 * t30;
t67 = t51 * t105 + t48 * t17;
t98 = qJD(2) * t50;
t87 = t46 * t98;
t3 = t67 * qJD(5) - t48 * t87 - t51 * t62;
t66 = t48 * t105 - t51 * t17;
t4 = t66 * qJD(5) - t48 * t62 + t51 * t87;
t114 = (-t48 * t67 + t51 * t66) * qJD(5) + t3 * t48 - t4 * t51;
t101 = -qJ(4) - pkin(7);
t107 = t45 * t49;
t36 = t101 * t52;
t23 = t101 * t107 - t99 * t36;
t80 = qJD(3) * t101;
t26 = t52 * qJD(4) + t49 * t80;
t60 = t49 * qJD(4) - t52 * t80;
t54 = t99 * t26 - t45 * t60;
t81 = t99 * t52;
t32 = -t81 + t107;
t82 = t99 * t49;
t33 = t45 * t52 + t82;
t41 = -t52 * pkin(3) - pkin(2);
t61 = -t32 * pkin(4) + t33 * pkin(8) - t41;
t58 = t51 * t61;
t27 = t33 * qJD(3);
t94 = t49 * qJD(3);
t28 = qJD(3) * t81 - t45 * t94;
t89 = pkin(3) * t94;
t59 = t27 * pkin(4) - t28 * pkin(8) + t89;
t96 = qJD(5) * t48;
t1 = qJD(5) * t58 + t23 * t96 - t48 * t59 - t51 * t54;
t6 = t51 * t23 - t48 * t61;
t2 = -t6 * qJD(5) - t48 * t54 + t51 * t59;
t5 = -t48 * t23 - t58;
t113 = t1 * t48 - t2 * t51 + (t48 * t5 - t51 * t6) * qJD(5);
t16 = -t99 * t29 + t45 * t30;
t8 = t45 * t20 - t99 * t21;
t112 = t16 * t8;
t14 = t45 * t26 + t99 * t60;
t22 = -t101 * t82 - t45 * t36;
t111 = t22 * t14;
t110 = t33 * t28;
t109 = t33 * t48;
t108 = t33 * t51;
t104 = t48 * t27;
t103 = t51 * t27;
t102 = t51 * t28;
t95 = qJD(5) * t51;
t93 = t52 * qJD(3);
t92 = 0.2e1 * t32 * t27;
t91 = -0.2e1 * pkin(2) * qJD(3);
t40 = -t99 * pkin(3) - pkin(4);
t90 = 0.2e1 * qJD(5) * t40;
t88 = t53 * t94;
t85 = t48 * t95;
t84 = t49 * t93;
t83 = -0.4e1 * t48 * t108;
t78 = t46 ^ 2 * t50 * t97;
t31 = t33 ^ 2;
t77 = t31 * t85;
t74 = -t48 * t6 - t5 * t51;
t72 = t16 * t14 + t8 * t22;
t71 = t48 * t66 + t51 * t67;
t39 = t45 * pkin(3) + pkin(8);
t69 = -t27 * t39 + t28 * t40;
t68 = t32 * t39 - t33 * t40;
t65 = t32 * t95 + t104;
t64 = t48 * t28 + t33 * t95;
t63 = -t33 * t96 + t102;
t57 = t74 * qJD(5) - t1 * t51 - t2 * t48;
t56 = t71 * qJD(5) - t3 * t51 - t4 * t48;
t55 = t20 * t52 - t21 * t49 + (-t29 * t52 - t30 * t49) * qJD(3);
t18 = -t32 * t96 + t103;
t9 = -t48 * t102 + t33 * t79;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t30 * t20 + 0.2e1 * t29 * t21 - 0.2e1 * t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t17 * t62 + 0.2e1 * t112 - 0.2e1 * t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t3 * t66 - 0.2e1 * t4 * t67 + 0.2e1 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86, 0, 0, 0, 0, 0, 0, 0, 0, (-t52 * t98 - t88) * t46, (t49 * t98 - t53 * t93) * t46, t55, -pkin(2) * t87 + t55 * pkin(7), 0, 0, 0, 0, 0, 0, (-t27 * t53 + t32 * t98) * t46, (-t28 * t53 + t33 * t98) * t46, t16 * t28 - t17 * t27 - t62 * t32 + t8 * t33, t62 * t23 + t17 * t54 + (-pkin(3) * t88 + t41 * t98) * t46 + t72, 0, 0, 0, 0, 0, 0, t8 * t109 + t64 * t16 - t27 * t67 + t4 * t32, t8 * t108 + t63 * t16 + t27 * t66 + t3 * t32, t114 * t33 + t71 * t28, t1 * t66 - t2 * t67 - t3 * t6 + t4 * t5 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t84, 0.2e1 * (-t49 ^ 2 + t52 ^ 2) * qJD(3), 0, -0.2e1 * t84, 0, 0, t49 * t91, t52 * t91, 0, 0, 0.2e1 * t110, -0.2e1 * t33 * t27 - 0.2e1 * t28 * t32, 0, t92, 0, 0, 0.2e1 * t41 * t27 + 0.2e1 * t32 * t89, 0.2e1 * t41 * t28 + 0.2e1 * t33 * t89, 0.2e1 * t14 * t33 + 0.2e1 * t22 * t28 - 0.2e1 * t23 * t27 - 0.2e1 * t54 * t32, 0.2e1 * t23 * t54 + 0.2e1 * t41 * t89 + 0.2e1 * t111, 0.2e1 * t44 * t110 - 0.2e1 * t77, t28 * t83 + 0.2e1 * t31 * t79, 0.2e1 * t33 * t103 + 0.2e1 * t63 * t32, 0.2e1 * t43 * t110 + 0.2e1 * t77, -0.2e1 * t33 * t104 - 0.2e1 * t64 * t32, t92, 0.2e1 * t14 * t109 + 0.2e1 * t2 * t32 + 0.2e1 * t64 * t22 + 0.2e1 * t5 * t27, 0.2e1 * t1 * t32 + 0.2e1 * t14 * t108 + 0.2e1 * t63 * t22 - 0.2e1 * t6 * t27, 0.2e1 * t113 * t33 + 0.2e1 * t74 * t28, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t62, 0, (t62 * t45 - t8 * t99) * pkin(3), 0, 0, 0, 0, 0, 0, t16 * t96 - t8 * t51, t16 * t95 + t8 * t48, t56, t56 * t39 + t8 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, -t94, 0, -pkin(7) * t93, pkin(7) * t94, 0, 0, 0, 0, t28, 0, -t27, 0, -t14, -t54, (-t27 * t45 - t99 * t28) * pkin(3), (-t14 * t99 + t54 * t45) * pkin(3), -t9, qJD(5) * t83 - t100 * t28, t65, t9, t18, 0, -t14 * t51 + t69 * t48 + (t22 * t48 - t68 * t51) * qJD(5), t14 * t48 + t69 * t51 + (t22 * t51 + t68 * t48) * qJD(5), t57, t14 * t40 + t57 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t85, -0.2e1 * t79, 0, -0.2e1 * t85, 0, 0, t48 * t90, t51 * t90, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28, 0, t89, 0, 0, 0, 0, 0, 0, t18, -t65, (-t43 - t44) * t28, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, -t64, t27, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, -t96, 0, -t39 * t95, t39 * t96, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
