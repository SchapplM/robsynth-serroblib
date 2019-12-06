% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaDJ_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:33
% EndTime: 2019-12-05 18:25:38
% DurationCPUTime: 1.08s
% Computational Cost: add. (1197->117), mult. (2780->239), div. (0->0), fcn. (2312->6), ass. (0->100)
t53 = sin(qJ(5));
t49 = t53 ^ 2;
t56 = cos(qJ(5));
t51 = t56 ^ 2;
t100 = t49 - t51;
t115 = t100 * qJD(5);
t117 = t49 + t51;
t55 = sin(qJ(2));
t50 = t55 ^ 2;
t57 = cos(qJ(2));
t52 = t57 ^ 2;
t116 = (t50 - t52) * qJD(2);
t108 = cos(qJ(4));
t101 = pkin(3) + qJ(3);
t36 = t101 * t55;
t37 = t101 * t57;
t54 = sin(qJ(4));
t21 = t108 * t37 - t54 * t36;
t34 = t108 * t55 + t54 * t57;
t58 = pkin(2) + pkin(1);
t38 = t58 * t57;
t72 = t34 * pkin(4) + t38;
t10 = t56 * t21 - t53 * t72;
t65 = t56 * t72;
t9 = -t53 * t21 - t65;
t74 = t10 * t56 - t53 * t9;
t114 = qJD(2) + qJD(4);
t79 = qJD(2) * t101;
t92 = t57 * qJD(3);
t27 = -t55 * t79 + t92;
t93 = t55 * qJD(3);
t63 = t57 * t79 + t93;
t98 = qJD(4) * t54;
t61 = t108 * t27 - t37 * t98 - t54 * t63;
t80 = t108 * qJD(4);
t59 = t36 * t80 - t61;
t104 = t54 * t55;
t82 = t108 * t57;
t17 = -qJD(2) * t82 + t114 * t104 - t57 * t80;
t94 = t55 * qJD(2);
t35 = t58 * t94;
t67 = t17 * pkin(4) + t35;
t95 = qJD(5) * t53;
t2 = qJD(5) * t65 + t21 * t95 - t53 * t67 + t56 * t59;
t3 = -t10 * qJD(5) + t53 * t59 + t56 * t67;
t113 = t74 * qJD(5) - t2 * t53 + t3 * t56;
t83 = t108 * t36;
t20 = t54 * t37 + t83;
t8 = t21 * qJD(4) + t108 * t63 + t54 * t27;
t112 = t20 * t8;
t6 = t8 * t53;
t110 = t8 * t56;
t48 = qJD(5) * t56;
t109 = t20 * t48 + t6;
t106 = t34 * t17;
t18 = t114 * t34;
t105 = t53 * t18;
t103 = t56 * t17;
t102 = t56 * t18;
t97 = qJD(4) * t58;
t47 = t57 * qJD(2);
t91 = qJ(3) * qJD(2);
t33 = -t82 + t104;
t90 = 0.2e1 * t33 * t18;
t89 = t53 * t103;
t88 = t8 * t108;
t87 = t54 * t97;
t86 = t53 * t48;
t85 = t55 * t47;
t84 = t58 * t108;
t81 = qJD(5) * t108;
t78 = (t50 + t52) * qJD(3);
t29 = t34 ^ 2;
t77 = t29 * t86;
t76 = t58 * t80;
t73 = t10 * t53 + t56 * t9;
t71 = -t53 * t17 + t34 * t48;
t70 = -t34 * t95 - t103;
t12 = t33 * t48 + t105;
t69 = t33 * t95 - t102;
t68 = -t57 * t91 - t93;
t43 = t54 * t58 + pkin(4);
t66 = t33 * t43 + t34 * t84;
t64 = t117 * t108;
t22 = t64 * t97;
t1 = -t73 * qJD(5) - t2 * t56 - t3 * t53;
t62 = t108 * t17 + (-t108 * t33 + t34 * t54) * qJD(4);
t60 = -t18 * t43 + t62 * t58;
t42 = -0.2e1 * t85;
t41 = 0.2e1 * t85;
t40 = -0.2e1 * t86;
t39 = 0.2e1 * t86;
t32 = -0.2e1 * t116;
t31 = -0.2e1 * t115;
t24 = (-t53 * t81 - t56 * t98) * t58;
t23 = (t53 * t98 - t56 * t81) * t58;
t15 = t20 * t95;
t7 = t34 * t115 + t89;
t4 = t100 * t17 - 0.4e1 * t34 * t86;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t32, 0, t42, 0, 0, 0, 0, 0, 0, t41, t32, 0, t42, 0, 0, -0.4e1 * pkin(1) * t85, 0.2e1 * pkin(1) * t116, 0.2e1 * t78, -0.2e1 * pkin(1) ^ 2 * t85 + 0.2e1 * qJ(3) * t78, -0.2e1 * t106, 0.2e1 * t17 * t33 - 0.2e1 * t34 * t18, 0, t90, 0, 0, -0.2e1 * t38 * t18 + 0.2e1 * t35 * t33, 0.2e1 * t38 * t17 + 0.2e1 * t35 * t34, -0.2e1 * t20 * t17 - 0.2e1 * t21 * t18 + 0.2e1 * t59 * t33 + 0.2e1 * t8 * t34, -0.2e1 * t21 * t59 - 0.2e1 * t38 * t35 + 0.2e1 * t112, -0.2e1 * t51 * t106 - 0.2e1 * t77, 0.2e1 * t29 * t115 + 0.4e1 * t34 * t89, 0.2e1 * t34 * t102 + 0.2e1 * t70 * t33, -0.2e1 * t49 * t106 + 0.2e1 * t77, -0.2e1 * t34 * t105 - 0.2e1 * t71 * t33, t90, 0.2e1 * t9 * t18 + 0.2e1 * t71 * t20 + 0.2e1 * t3 * t33 + 0.2e1 * t34 * t6, -0.2e1 * t10 * t18 + 0.2e1 * t34 * t110 + 0.2e1 * t2 * t33 + 0.2e1 * t70 * t20, -0.2e1 * t113 * t34 + 0.2e1 * t73 * t17, -0.2e1 * t10 * t2 + 0.2e1 * t9 * t3 + 0.2e1 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t94, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t94, 0, t68, t55 * t91 - t92, -pkin(1) * t47, t68 * pkin(1), 0, 0, -t17, 0, -t18, 0, -t8, t59, (-t18 * t54 + t62) * t58, (t61 * t54 - t88 + (t21 * t108 + (-t83 + t20) * t54) * qJD(4)) * t58, -t7, t4, t12, t7, -t69, 0, t15 + (-t66 * qJD(5) - t8) * t56 + t60 * t53, t60 * t56 + t66 * t95 + t109, t1, (-t88 + (t74 * t108 + t20 * t54) * qJD(4)) * t58 + t1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t87, -0.2e1 * t76, 0, 0, t39, t31, 0, t40, 0, 0, 0.2e1 * t24, 0.2e1 * t23, 0.2e1 * t22, 0.2e1 * (t64 * t43 - t54 * t84) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t47, 0, pkin(1) * t94, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t35, 0, 0, 0, 0, 0, 0, -t69, -t12, t117 * t17, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t18, 0, -t8, t59, 0, 0, -t7, t4, t12, t7, -t69, 0, -t12 * pkin(4) - t110 + t15, t69 * pkin(4) + t109, t1, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t76, 0, 0, t39, t31, 0, t40, 0, 0, t24, t23, t22, pkin(4) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t31, 0, t40, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t71, t18, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t95, 0, -t43 * t48 - t53 * t76, t43 * t95 - t56 * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t95, 0, -pkin(4) * t48, pkin(4) * t95, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
