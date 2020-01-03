% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR12
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:09
% EndTime: 2019-12-31 19:13:13
% DurationCPUTime: 1.19s
% Computational Cost: add. (1763->146), mult. (3507->256), div. (0->0), fcn. (3076->6), ass. (0->98)
t115 = qJD(3) + qJD(4);
t54 = sin(qJ(4));
t55 = sin(qJ(3));
t57 = cos(qJ(4));
t58 = cos(qJ(3));
t34 = t54 * t58 + t57 * t55;
t23 = t115 * t34;
t101 = t57 * t58;
t92 = t55 * qJD(3);
t95 = qJD(4) * t54;
t24 = t115 * t101 - t54 * t92 - t55 * t95;
t35 = -t54 * t55 + t101;
t104 = t35 * t54;
t66 = (t34 * t57 - t104) * qJD(4);
t118 = (t23 * t57 - t24 * t54 - t66) * pkin(3);
t53 = sin(qJ(5));
t51 = t53 ^ 2;
t56 = cos(qJ(5));
t52 = t56 ^ 2;
t98 = t51 - t52;
t116 = t98 * qJD(5);
t32 = t35 ^ 2;
t114 = 2 * qJD(2);
t113 = t23 * pkin(4);
t59 = -pkin(1) - pkin(6);
t112 = pkin(7) - t59;
t36 = t112 * t55;
t80 = t112 * t58;
t26 = -t57 * t36 - t54 * t80;
t31 = qJD(3) * t80;
t77 = t112 * t92;
t13 = t26 * qJD(4) - t54 * t31 - t57 * t77;
t111 = t13 * t35;
t10 = t13 * t53;
t48 = -t57 * pkin(3) - pkin(4);
t110 = t23 * t48;
t109 = t23 * t56;
t25 = -t54 * t36 + t57 * t80;
t108 = t25 * t13;
t107 = t25 * t54;
t106 = t34 * t24;
t105 = t35 * t23;
t103 = t53 * t24;
t102 = t56 * t24;
t49 = qJD(5) * t56;
t100 = t25 * t49 + t10;
t85 = pkin(3) * t95;
t99 = t48 * t49 + t53 * t85;
t97 = t51 + t52;
t96 = pkin(3) * qJD(4);
t93 = qJD(5) * t53;
t91 = t58 * qJD(3);
t90 = qJ(2) * qJD(3);
t89 = 0.2e1 * t106;
t88 = t53 * t109;
t87 = pkin(4) * t93;
t86 = pkin(4) * t49;
t84 = t57 * t96;
t83 = t53 * t49;
t82 = t55 * t91;
t81 = t34 ^ 2 + t32;
t11 = t97 * t24;
t79 = t97 * t57;
t45 = t55 * pkin(3) + qJ(2);
t78 = t32 * t83;
t64 = t34 * pkin(4) - t35 * pkin(8) + t45;
t62 = t56 * t64;
t7 = -t53 * t26 + t62;
t8 = t56 * t26 + t53 * t64;
t76 = t53 * t8 + t56 * t7;
t75 = t53 * t7 - t56 * t8;
t42 = pkin(3) * t91 + qJD(2);
t74 = t25 * t23 - t111;
t73 = t105 - t106;
t47 = t54 * pkin(3) + pkin(8);
t71 = t34 * t47 - t35 * t48;
t69 = t48 * t93 - t56 * t85;
t68 = -t23 * t53 + t35 * t49;
t67 = t35 * t93 + t109;
t14 = t34 * t93 - t102;
t65 = 0.2e1 * t73;
t63 = t24 * pkin(4) + t23 * pkin(8) + t42;
t12 = t25 * qJD(4) + t57 * t31 - t54 * t77;
t61 = t12 * t34 - t26 * t24 - t74;
t2 = -qJD(5) * t62 + t56 * t12 + t26 * t93 - t53 * t63;
t3 = -t8 * qJD(5) + t53 * t12 + t56 * t63;
t1 = -t76 * qJD(5) - t2 * t56 - t3 * t53;
t60 = -pkin(3) * t66 - t24 * t47 - t110;
t50 = qJ(2) * t114;
t41 = -0.2e1 * t83;
t40 = 0.2e1 * t83;
t33 = -0.2e1 * t116;
t28 = t79 * t96;
t20 = t25 * t93;
t15 = t34 * t49 + t103;
t6 = t35 * t116 + t88;
t4 = t98 * t23 - 0.4e1 * t35 * t83;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t50, -0.2e1 * t82, 0.2e1 * (t55 ^ 2 - t58 ^ 2) * qJD(3), 0, 0.2e1 * t82, 0, 0, 0.2e1 * qJD(2) * t55 + 0.2e1 * t58 * t90, 0.2e1 * qJD(2) * t58 - 0.2e1 * t55 * t90, 0, t50, -0.2e1 * t105, 0.2e1 * t23 * t34 - 0.2e1 * t35 * t24, 0, t89, 0, 0, 0.2e1 * t45 * t24 + 0.2e1 * t42 * t34, -0.2e1 * t45 * t23 + 0.2e1 * t42 * t35, 0.2e1 * t61, -0.2e1 * t26 * t12 + 0.2e1 * t45 * t42 + 0.2e1 * t108, -0.2e1 * t52 * t105 - 0.2e1 * t78, 0.2e1 * t32 * t116 + 0.4e1 * t35 * t88, 0.2e1 * t35 * t102 - 0.2e1 * t67 * t34, -0.2e1 * t51 * t105 + 0.2e1 * t78, -0.2e1 * t35 * t103 - 0.2e1 * t68 * t34, t89, 0.2e1 * t35 * t10 + 0.2e1 * t7 * t24 + 0.2e1 * t68 * t25 + 0.2e1 * t3 * t34, 0.2e1 * t56 * t111 + 0.2e1 * t2 * t34 - 0.2e1 * t8 * t24 - 0.2e1 * t25 * t67, 0.2e1 * t76 * t23 + 0.2e1 * (qJD(5) * t75 + t2 * t53 - t3 * t56) * t35, -0.2e1 * t8 * t2 + 0.2e1 * t7 * t3 + 0.2e1 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t61, 0, 0, 0, 0, 0, 0, -t81 * t49 + t53 * t65, t56 * t65 + t81 * t93, 0, t1 * t34 - t24 * t75 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t11 * t34 - 0.2e1 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, 0, -t91, 0, -t59 * t92, -t59 * t91, 0, 0, 0, 0, -t23, 0, -t24, 0, -t13, t12, t118, (-t12 * t54 - t13 * t57 + (t26 * t57 + t107) * qJD(4)) * pkin(3), -t6, t4, t15, t6, -t14, 0, t20 + (-t71 * qJD(5) - t13) * t56 + t60 * t53, t56 * t60 + t71 * t93 + t100, t1, t13 * t48 + (-t57 * t75 + t107) * t96 + t1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, -t118, 0, 0, 0, 0, 0, 0, -t67, -t68, t11, t110 + t47 * t11 + (t34 * t79 - t104) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t85, -0.2e1 * t84, 0, 0, t40, t33, 0, t41, 0, 0, 0.2e1 * t69, 0.2e1 * t99, 0.2e1 * t28, 0.2e1 * (t47 * t79 + t48 * t54) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, 0, -t13, t12, 0, 0, -t6, t4, t15, t6, -t14, 0, t20 + (-pkin(8) * t24 + t113) * t53 + (-t13 + (-pkin(4) * t35 - pkin(8) * t34) * qJD(5)) * t56, pkin(4) * t67 + pkin(8) * t14 + t100, t1, -t13 * pkin(4) + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t68, t11, pkin(8) * t11 - t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t84, 0, 0, t40, t33, 0, t41, 0, 0, t69 - t87, -t86 + t99, t28, (-pkin(4) * t54 + pkin(8) * t79) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t33, 0, t41, 0, 0, -0.2e1 * t87, -0.2e1 * t86, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, -t68, t24, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, -t93, 0, -t47 * t49 - t53 * t84, t47 * t93 - t56 * t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, -t93, 0, -pkin(8) * t49, pkin(8) * t93, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
