% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:37
% DurationCPUTime: 0.55s
% Computational Cost: add. (417->87), mult. (1008->157), div. (0->0), fcn. (717->6), ass. (0->72)
t87 = 2 * qJD(5);
t44 = cos(pkin(9));
t47 = cos(qJ(5));
t46 = sin(qJ(3));
t73 = pkin(2) * qJD(3);
t63 = t46 * t73;
t38 = t46 * pkin(2) + qJ(4);
t45 = sin(qJ(5));
t43 = sin(pkin(9));
t52 = t44 * pkin(4) + t43 * pkin(7) + pkin(3);
t48 = cos(qJ(3));
t84 = t48 * pkin(2);
t50 = t52 + t84;
t80 = t44 * t47;
t9 = t38 * t80 - t45 * t50;
t70 = qJD(5) * t9;
t62 = t48 * t73;
t31 = qJD(4) + t62;
t79 = t45 * t31;
t3 = -t44 * t79 + t47 * t63 - t70;
t65 = t45 * qJD(4);
t71 = qJ(4) * t44;
t12 = -t45 * t52 + t47 * t71;
t69 = qJD(5) * t12;
t7 = -t44 * t65 - t69;
t86 = -t3 - t7;
t68 = qJD(5) * t45;
t35 = t44 * t68;
t49 = t47 * t50;
t2 = qJD(5) * t49 - t31 * t80 + t38 * t35 - t45 * t63;
t85 = t2 * t45;
t51 = t47 * t52;
t57 = qJ(4) * t68;
t66 = t44 * qJD(4);
t6 = qJD(5) * t51 + t44 * t57 - t47 * t66;
t83 = t6 * t45;
t41 = t43 ^ 2;
t24 = t41 * t31;
t82 = -t2 * t44 + t47 * t24;
t39 = t41 * qJD(4);
t81 = t47 * t39 - t6 * t44;
t42 = t44 ^ 2;
t25 = t42 * t31;
t67 = qJD(5) * t47;
t60 = t41 * t67;
t78 = t38 * t60 + t41 * t79;
t72 = qJ(4) * t31;
t77 = t38 * t39 + t41 * t72;
t76 = t25 + t24;
t75 = qJ(4) * t60 + t41 * t65;
t74 = t42 * qJD(4) + t39;
t64 = qJ(4) * qJD(4);
t61 = t41 * t68;
t59 = t43 * t68;
t58 = t43 * t67;
t56 = t43 * t44 * t87;
t55 = t43 * t63;
t54 = t44 * t63;
t53 = t45 * t60;
t37 = t41 * t64;
t36 = t44 * t67;
t28 = -0.2e1 * t53;
t27 = 0.2e1 * t53;
t23 = t47 * t56;
t22 = t45 * t56;
t15 = t38 * t24;
t14 = (t45 ^ 2 - t47 ^ 2) * t41 * t87;
t11 = -t45 * t71 - t51;
t10 = t11 * t59;
t8 = -t45 * t44 * t38 - t49;
t4 = t8 * t59;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t2 * t47 - t3 * t45 - t44 * t31 + (-t45 * t9 - t47 * t8) * qJD(5)) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63, -0.2e1 * t62, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t54, 0.2e1 * t55, 0.2e1 * t76, 0.2e1 * t38 * t25 + 0.2e1 * t15 + 0.2e1 * (-pkin(3) - t84) * t63, t28, t14, t22, t27, t23, 0, -0.2e1 * t3 * t44 + 0.2e1 * t78, -0.2e1 * t38 * t61 + 0.2e1 * t82, 0.2e1 * t4 + 0.2e1 * (t85 + (-t3 - t70) * t47) * t43, -0.2e1 * t9 * t2 + 0.2e1 * t8 * t3 + 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t66 - t45 * t7 - t47 * t6 + (-t11 * t47 - t12 * t45) * qJD(5)) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t55, t74 + t76, -pkin(3) * t63 + (qJD(4) * t38 + t72) * t42 + t77, t28, t14, t22, t27, t23, 0, t86 * t44 + t75 + t78, (-qJ(4) - t38) * t61 + t81 + t82, t10 + t4 + ((t2 + t6) * t45 + ((-t12 - t9) * qJD(5) + t86) * t47) * t43, t3 * t11 - t2 * t12 - t9 * t6 + t8 * t7 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t74, 0.2e1 * t42 * t64 + 0.2e1 * t37, t28, t14, t22, t27, t23, 0, -0.2e1 * t7 * t44 + 0.2e1 * t75, -0.2e1 * t41 * t57 + 0.2e1 * t81, 0.2e1 * t10 + 0.2e1 * (t83 + (-t7 - t69) * t47) * t43, 0.2e1 * t11 * t7 - 0.2e1 * t12 * t6 + 0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, t35, t36, 0, -t85 + t3 * t47 + (-t45 * t8 + t47 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, 0, -t83 + t7 * t47 + (-t11 * t45 + t12 * t47) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, -t58, 0, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, -t58, 0, t7, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
