% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:35
% EndTime: 2019-12-05 15:51:39
% DurationCPUTime: 0.94s
% Computational Cost: add. (752->122), mult. (2179->238), div. (0->0), fcn. (2136->10), ass. (0->90)
t36 = sin(pkin(10));
t37 = sin(pkin(5));
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t85 = cos(pkin(10));
t17 = (t36 * t43 + t85 * t40) * t37;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t86 = cos(pkin(5));
t49 = -t17 * t39 + t86 * t42;
t83 = t49 * qJD(4);
t96 = -0.4e1 * t39;
t38 = sin(qJ(5));
t32 = t38 ^ 2;
t41 = cos(qJ(5));
t34 = t41 ^ 2;
t62 = qJD(5) * (t32 - t34);
t95 = pkin(4) * t39;
t94 = pkin(8) * t42;
t12 = t17 * t42 + t86 * t39;
t63 = t85 * t43;
t84 = qJD(2) * t37;
t68 = t40 * t84;
t14 = t36 * t68 - t63 * t84;
t5 = t12 * qJD(4) - t14 * t39;
t93 = t49 * t5;
t92 = t5 * t39;
t15 = qJD(2) * t17;
t16 = (t36 * t40 - t63) * t37;
t91 = t16 * t15;
t30 = t36 * pkin(2) + pkin(7);
t90 = t38 * t30;
t89 = t41 * t42;
t33 = t39 ^ 2;
t87 = -t42 ^ 2 + t33;
t82 = qJD(4) * t41;
t81 = qJD(5) * t38;
t80 = qJD(5) * t41;
t79 = qJD(5) * t42;
t78 = t39 * qJD(4);
t77 = t42 * qJD(4);
t76 = -0.2e1 * pkin(4) * qJD(5);
t75 = t42 * t90;
t74 = t30 * t89;
t31 = -t85 * pkin(2) - pkin(3);
t73 = 0.2e1 * qJD(4) * t31;
t72 = qJD(5) * t30 * t33;
t71 = t38 * t79;
t70 = t41 * t79;
t69 = t32 * t77;
t67 = t38 * t80;
t66 = t39 * t77;
t65 = t41 * t77;
t64 = t30 * t77;
t61 = t87 * qJD(4);
t60 = t38 * t65;
t59 = t33 * t67;
t58 = -t42 * pkin(4) - t39 * pkin(8);
t57 = -t94 + t95;
t7 = -t12 * t38 + t16 * t41;
t8 = t12 * t41 + t16 * t38;
t56 = -t38 * t8 - t41 * t7;
t55 = t38 * t7 - t41 * t8;
t48 = -t31 - t58;
t10 = -t38 * t48 + t74;
t47 = t41 * t48;
t9 = -t47 - t75;
t54 = t10 * t41 - t38 * t9;
t53 = -t10 * t38 - t41 * t9;
t52 = t5 * t38 - t49 * t80;
t51 = -t5 * t41 - t49 * t81;
t50 = t57 * t38;
t20 = t41 * t78 + t71;
t6 = -t14 * t42 + t83;
t1 = -t8 * qJD(5) + t15 * t41 - t6 * t38;
t2 = t7 * qJD(5) + t15 * t38 + t6 * t41;
t46 = t56 * qJD(5) - t1 * t38 + t2 * t41;
t3 = -qJD(4) * t50 + qJD(5) * t47 + t20 * t30;
t4 = -t10 * qJD(5) + (t39 * t90 + t41 * t57) * qJD(4);
t45 = t53 * qJD(5) - t3 * t41 - t4 * t38;
t44 = t92 + t6 * t42 + (-t12 * t39 - t42 * t49) * qJD(4);
t28 = t34 * t77;
t27 = -0.2e1 * t66;
t26 = t34 * t66;
t25 = t32 * t66;
t22 = t38 * t78 - t70;
t21 = -t38 * t77 - t39 * t80;
t19 = t39 * t81 - t65;
t13 = t39 * t62 - t60;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t14 + 0.2e1 * t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 * t6 + 0.2e1 * t91 - 0.2e1 * t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t7 * t1 + 0.2e1 * t8 * t2 - 0.2e1 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t43 * t84, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t14, 0, (-t14 * t36 - t85 * t15) * pkin(2), 0, 0, 0, 0, 0, 0, -t15 * t42 + t16 * t78, t15 * t39 + t16 * t77, t44, t15 * t31 + t44 * t30, 0, 0, 0, 0, 0, 0, (-t38 * t83 - t1) * t42 + (qJD(4) * t7 + t52) * t39, (-t49 * t82 + t2) * t42 + (-qJD(4) * t8 - t51) * t39, t56 * t77 + (t55 * qJD(5) - t1 * t41 - t2 * t38) * t39, t1 * t9 + t2 * t10 - t8 * t3 + t7 * t4 + (-t49 * t77 + t92) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t66, -0.2e1 * t61, 0, t27, 0, 0, t39 * t73, t42 * t73, 0, 0, 0.2e1 * t26 - 0.2e1 * t59, 0.2e1 * t33 * t62 + t60 * t96, 0.2e1 * t39 * t71 + 0.2e1 * t87 * t82, 0.2e1 * t25 + 0.2e1 * t59, -0.2e1 * t38 * t61 + 0.2e1 * t39 * t70, t27, 0.2e1 * t41 * t72 - 0.2e1 * t4 * t42 + 0.2e1 * (t9 + 0.2e1 * t75) * t78, -0.2e1 * t38 * t72 - 0.2e1 * t3 * t42 + 0.2e1 * (-t10 + 0.2e1 * t74) * t78, 0.2e1 * t53 * t77 + 0.2e1 * (-t54 * qJD(5) + t3 * t38 - t4 * t41) * t39, 0.2e1 * t30 ^ 2 * t66 - 0.2e1 * t10 * t3 + 0.2e1 * t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t39 - t5 * t42 + (t12 * t42 - t39 * t49) * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t55 * qJD(4) - t5) * t42 + (t46 - t83) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t39 + (t87 * t30 + t54 * t42) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t25 + 0.2e1 * t26 - 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, 0, 0, 0, t51, t52, t46, -t5 * pkin(4) + t46 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, -t78, 0, -t64, t30 * t78, 0, 0, -t13, t67 * t96 + t28 - t69, t22, t13, t20, 0, (pkin(8) * t89 + (-t41 * pkin(4) + t90) * t39) * qJD(5) + (t58 * t38 - t74) * qJD(4), (t30 * t39 * t41 + t50) * qJD(5) + (t58 * t41 + t75) * qJD(4), t45, -pkin(4) * t64 + t45 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t22, t28 + t69, (-t95 + (t32 + t34) * t94) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67, -0.2e1 * t62, 0, -0.2e1 * t67, 0, 0, t38 * t76, t41 * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, t21, t78, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, -t81, 0, -pkin(8) * t80, pkin(8) * t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
