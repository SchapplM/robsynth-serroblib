% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:32
% EndTime: 2019-12-05 15:47:37
% DurationCPUTime: 0.67s
% Computational Cost: add. (628->87), mult. (1562->169), div. (0->0), fcn. (1515->8), ass. (0->58)
t32 = sin(pkin(9));
t33 = cos(pkin(9));
t36 = sin(qJ(2));
t38 = cos(qJ(2));
t22 = t32 * t38 + t33 * t36;
t37 = cos(qJ(4));
t65 = cos(qJ(5));
t51 = t65 * t37;
t34 = sin(qJ(5));
t35 = sin(qJ(4));
t62 = t34 * t35;
t43 = t51 - t62;
t8 = t43 * t22;
t67 = qJD(4) + qJD(5);
t21 = t32 * t36 - t33 * t38;
t19 = t21 * qJD(2);
t30 = t35 ^ 2;
t31 = t37 ^ 2;
t50 = (t31 + t30) * t19;
t66 = 2 * qJD(4);
t18 = t22 * qJD(2);
t13 = t21 * t18;
t24 = t34 * t37 + t65 * t35;
t12 = t67 * t24;
t64 = t43 * t12;
t49 = t65 * qJD(5);
t11 = -qJD(4) * t51 - t37 * t49 + t67 * t62;
t63 = t24 * t11;
t61 = t35 * t19;
t60 = qJD(5) * t34;
t59 = t35 * qJD(4);
t58 = t37 * qJD(4);
t28 = -t33 * pkin(2) - pkin(3);
t57 = t28 * t66;
t56 = pkin(4) * t59;
t55 = pkin(4) * t60;
t54 = t21 * t59;
t53 = t35 * t58;
t52 = t32 * pkin(2) + pkin(6);
t48 = pkin(7) + t52;
t47 = pkin(4) * t49;
t46 = t52 * qJD(4);
t45 = t34 * t48;
t44 = t35 * t45;
t42 = qJD(4) * t45;
t41 = t48 * t65;
t40 = t35 * t41;
t39 = qJD(4) * t41;
t25 = -t37 * pkin(4) + t28;
t20 = t48 * t37;
t10 = t65 * t20 - t44;
t9 = -t34 * t20 - t40;
t7 = t24 * t22;
t4 = qJD(5) * t44 - t20 * t49 + t35 * t42 - t37 * t39;
t3 = qJD(5) * t40 + t20 * t60 + t35 * t39 + t37 * t42;
t2 = t24 * t19 - t67 * t8;
t1 = t12 * t22 + t19 * t51 - t34 * t61;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t22 * t19 + 0.2e1 * t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t22 * t50 + 0.2e1 * t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8 * t1 - 0.2e1 * t7 * t2 + 0.2e1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * qJD(2), -t38 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t18, t19, 0, (-t18 * t33 - t19 * t32) * pkin(2), 0, 0, 0, 0, 0, 0, -t18 * t37 + t54, t18 * t35 + t21 * t58, -t50, t18 * t28 - t52 * t50, 0, 0, 0, 0, 0, 0, t21 * t12 - t18 * t43, -t21 * t11 + t18 * t24, -t1 * t43 - t7 * t11 - t8 * t12 - t2 * t24, pkin(4) * t54 - t1 * t10 + t18 * t25 + t2 * t9 - t8 * t3 - t7 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t53, (-t30 + t31) * t66, 0, -0.2e1 * t53, 0, 0, t35 * t57, t37 * t57, 0, 0, -0.2e1 * t63, -0.2e1 * t11 * t43 - 0.2e1 * t24 * t12, 0, -0.2e1 * t64, 0, 0, 0.2e1 * t25 * t12 - 0.2e1 * t43 * t56, -0.2e1 * t25 * t11 + 0.2e1 * t24 * t56, -0.2e1 * t10 * t12 + 0.2e1 * t9 * t11 - 0.2e1 * t4 * t24 - 0.2e1 * t3 * t43, -0.2e1 * t10 * t3 + 0.2e1 * t25 * t56 + 0.2e1 * t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t24 - t8 * t11 + t7 * t12 + t2 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t11 - t9 * t12 - t3 * t24 + t4 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63 - 0.2e1 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t58 + t61, t37 * t19 + t22 * t59, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, (t65 * t2 - t1 * t34 + (t34 * t7 + t65 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t59, 0, -t37 * t46, t35 * t46, 0, 0, 0, 0, -t11, 0, -t12, 0, t4, t3, (t65 * t11 - t12 * t34 + (t24 * t34 + t43 * t65) * qJD(5)) * pkin(4), (t65 * t4 - t3 * t34 + (t65 * t10 - t34 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, (-t65 * t12 - t11 * t34 + (t65 * t24 - t34 * t43) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, -0.2e1 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
