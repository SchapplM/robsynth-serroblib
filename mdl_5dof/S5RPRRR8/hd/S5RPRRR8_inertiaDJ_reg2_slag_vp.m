% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:05
% EndTime: 2019-12-31 19:06:10
% DurationCPUTime: 1.23s
% Computational Cost: add. (1289->141), mult. (2579->253), div. (0->0), fcn. (2109->6), ass. (0->98)
t108 = cos(qJ(5));
t50 = sin(qJ(5));
t51 = sin(qJ(4));
t102 = t50 * t51;
t114 = qJD(4) + qJD(5);
t53 = cos(qJ(4));
t78 = qJD(5) * t108;
t81 = t108 * t53;
t12 = -qJD(4) * t81 + t102 * t114 - t53 * t78;
t32 = t108 * t51 + t50 * t53;
t13 = t114 * t32;
t67 = t81 - t102;
t115 = ((-t108 * t67 - t32 * t50) * qJD(5) - t108 * t12 + t13 * t50) * pkin(4);
t52 = sin(qJ(3));
t22 = t67 * t52;
t54 = cos(qJ(3));
t111 = pkin(1) + pkin(2);
t86 = t54 * t111;
t97 = qJD(3) * t52;
t19 = qJ(2) * t97 - qJD(2) * t54 + qJD(3) * t86;
t48 = t51 ^ 2;
t49 = t53 ^ 2;
t98 = t49 + t48;
t11 = t98 * t19;
t80 = t98 * t54;
t112 = 0.2e1 * qJD(2);
t110 = -pkin(8) - pkin(7);
t109 = t53 * pkin(4);
t107 = t12 * t32;
t35 = t54 * qJ(2) - t111 * t52;
t20 = t52 * qJD(2) + qJD(3) * t35;
t106 = t20 * t51;
t105 = t20 * t53;
t104 = t20 * t54;
t103 = t67 * t13;
t101 = t51 * t19;
t100 = t53 * t19;
t34 = -qJ(2) * t52 - t86;
t33 = pkin(3) - t34;
t23 = t33 + t109;
t45 = -pkin(3) - t109;
t99 = t23 - t45;
t96 = qJD(3) * t54;
t95 = qJD(5) * t50;
t94 = t51 * qJD(4);
t93 = t53 * qJD(4);
t92 = -0.2e1 * t107;
t91 = -0.2e1 * t103;
t90 = -0.2e1 * pkin(3) * qJD(4);
t89 = pkin(4) * t94;
t88 = pkin(4) * t95;
t87 = t50 * t110;
t85 = t51 * t96;
t84 = t51 * t93;
t83 = t52 * t96;
t82 = t54 * t94;
t79 = qJD(4) * (pkin(3) + t33);
t77 = t51 * t87;
t76 = pkin(4) * t78;
t75 = t110 * t108;
t17 = t20 - t89;
t74 = t17 - t89;
t73 = -t12 * t67 - t13 * t32;
t72 = t51 * t75;
t71 = t12 * t54 + t32 * t97;
t70 = -t13 * t54 - t67 * t97;
t66 = -pkin(7) + t35;
t65 = 0.2e1 * t73;
t64 = -pkin(8) + t66;
t38 = t110 * t53;
t15 = -t108 * t38 + t77;
t62 = t66 * qJD(4);
t61 = t64 * t51;
t21 = t32 * t52;
t3 = t13 * t52 + t50 * t85 - t81 * t96;
t4 = -t114 * t22 - t32 * t96;
t60 = -t12 * t21 - t13 * t22 - t3 * t67 - t32 * t4;
t59 = qJD(4) * t64;
t58 = t50 * t61;
t57 = t108 * t61;
t56 = -t53 * t59 + t101;
t55 = t51 * t59 + t100;
t40 = -0.2e1 * t84;
t39 = 0.2e1 * t84;
t36 = (-t48 + t49) * qJD(4);
t30 = 0.2e1 * t36;
t26 = t53 * t97 + t82;
t25 = t51 * t97 - t54 * t93;
t24 = qJD(3) * t80;
t18 = t64 * t53;
t14 = t50 * t38 + t72;
t10 = t108 * t18 - t58;
t9 = -t50 * t18 - t57;
t6 = -t15 * qJD(5) + (t53 * t75 - t77) * qJD(4);
t5 = -t114 * t72 - t38 * t95 - t87 * t93;
t2 = qJD(5) * t58 + t108 * t56 - t18 * t78 + t50 * t55;
t1 = qJD(5) * t57 + t108 * t55 + t18 * t95 - t50 * t56;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, qJ(2) * t112, 0, 0, 0, 0, 0, 0, 0.2e1 * t20, -0.2e1 * t19, 0, -0.2e1 * t19 * t35 - 0.2e1 * t20 * t34, t39, t30, 0, t40, 0, 0, -0.2e1 * t33 * t94 + 0.2e1 * t105, -0.2e1 * t33 * t93 - 0.2e1 * t106, 0.2e1 * t11, -0.2e1 * t11 * t66 + 0.2e1 * t33 * t20, t92, t65, 0, t91, 0, 0, -0.2e1 * t23 * t13 + 0.2e1 * t17 * t67, 0.2e1 * t23 * t12 - 0.2e1 * t17 * t32, 0.2e1 * t1 * t67 + 0.2e1 * t10 * t13 - 0.2e1 * t9 * t12 + 0.2e1 * t2 * t32, -0.2e1 * t1 * t10 + 0.2e1 * t17 * t23 + 0.2e1 * t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, 0, -t19 * t52 - t104 + (-t34 * t52 + t35 * t54) * qJD(3), 0, 0, 0, 0, 0, 0, t26, -t25, -t24, -t104 - t52 * t11 + (t33 * t52 + t66 * t80) * qJD(3), 0, 0, 0, 0, 0, 0, -t70, -t71, -t60, -t1 * t22 - t10 * t3 - t17 * t54 - t2 * t21 + t23 * t97 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t98) * t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21 * t4 - 0.2e1 * t22 * t3 - 0.2e1 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, 0, t40, -0.2e1 * t36, 0, t39, 0, 0, t51 * t79 - t105, t53 * t79 + t106, -t11, -t20 * pkin(3) - pkin(7) * t11, 0.2e1 * t107, -0.2e1 * t73, 0, 0.2e1 * t103, 0, 0, t13 * t99 - t67 * t74, -t12 * t99 + t32 * t74, (-t2 + t6) * t32 - (t1 - t5) * t67 + (-t10 + t15) * t13 + (-t14 + t9) * t12, -t1 * t15 - t10 * t5 + t14 * t2 + t17 * t45 + t23 * t89 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, t24, (-pkin(3) * t52 + pkin(7) * t80) * qJD(3), 0, 0, 0, 0, 0, 0, t70, t71, t60, -pkin(4) * t82 + t14 * t4 - t15 * t3 - t21 * t6 - t22 * t5 + t45 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t30, 0, t40, 0, 0, t51 * t90, t53 * t90, 0, 0, t92, t65, 0, t91, 0, 0, 0.2e1 * t45 * t13 - 0.2e1 * t67 * t89, -0.2e1 * t45 * t12 + 0.2e1 * t32 * t89, 0.2e1 * t14 * t12 - 0.2e1 * t15 * t13 - 0.2e1 * t6 * t32 - 0.2e1 * t5 * t67, 0.2e1 * t14 * t6 - 0.2e1 * t15 * t5 + 0.2e1 * t45 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, 0, t94, 0, -t53 * t62 + t101, t51 * t62 + t100, 0, 0, 0, 0, t12, 0, t13, 0, t2, t1, t115, (t108 * t2 - t1 * t50 + (t10 * t108 - t50 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t93 - t85, t52 * t94 - t53 * t96, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, (t108 * t4 - t3 * t50 + (t108 * t22 + t21 * t50) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, -t94, 0, -pkin(7) * t93, pkin(7) * t94, 0, 0, 0, 0, -t12, 0, -t13, 0, t6, t5, -t115, (t108 * t6 - t5 * t50 + (t108 * t15 - t14 * t50) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t88, -0.2e1 * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t13, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t13, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
