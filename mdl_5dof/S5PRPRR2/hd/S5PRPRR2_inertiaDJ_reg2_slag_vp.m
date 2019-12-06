% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR2
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:15
% DurationCPUTime: 0.45s
% Computational Cost: add. (391->51), mult. (1024->95), div. (0->0), fcn. (982->8), ass. (0->46)
t52 = cos(pkin(9));
t41 = t52 * pkin(2) + pkin(3);
t51 = sin(pkin(9));
t43 = t51 * pkin(2);
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t20 = t56 * t41 - t54 * t43;
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t23 = t51 * t57 + t52 * t55;
t38 = t51 * t55 - t52 * t57;
t36 = t56 * t38;
t7 = t54 * t23 + t36;
t35 = t54 * t38;
t8 = t56 * t23 - t35;
t37 = t23 * qJD(2);
t5 = -qJD(2) * t35 + t8 * qJD(4) + t56 * t37;
t60 = t7 * t5;
t21 = t54 * t41 + t56 * t43;
t13 = t21 * qJD(4);
t59 = t7 * t13;
t17 = -pkin(4) - t20;
t34 = cos(qJ(5));
t30 = t34 * qJD(5);
t33 = sin(qJ(5));
t58 = t13 * t33 + t17 * t30;
t31 = t33 ^ 2;
t32 = t34 ^ 2;
t53 = -t31 - t32;
t50 = t33 * qJD(5);
t49 = pkin(4) * t50;
t48 = pkin(4) * t30;
t47 = t33 * t30;
t4 = qJD(2) * t36 + t7 * qJD(4) + t54 * t37;
t1 = t53 * t4;
t12 = t20 * qJD(4);
t6 = t53 * t12;
t18 = pkin(7) + t21;
t44 = t53 * t18;
t42 = -t13 * t34 + t17 * t50;
t28 = -0.2e1 * t47;
t27 = 0.2e1 * t47;
t24 = 0.2e1 * (-t31 + t32) * qJD(5);
t3 = -t5 * t34 + t7 * t50;
t2 = t7 * t30 + t5 * t33;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8 * t4 + 0.2e1 * t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t8 * t1 + 0.2e1 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * qJD(2), -t57 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t37, t38 * qJD(2), 0, (-t23 * t52 - t38 * t51) * pkin(2) * qJD(2), 0, 0, 0, 0, 0, 0, -t5, t4, 0, t8 * t12 - t5 * t20 - t4 * t21 + t59, 0, 0, 0, 0, 0, 0, t3, t2, t1, t5 * t17 + t4 * t44 - t8 * t6 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13, -0.2e1 * t12, 0, 0.2e1 * t21 * t12 - 0.2e1 * t20 * t13, t27, t24, 0, t28, 0, 0, 0.2e1 * t42, 0.2e1 * t58, -0.2e1 * t6, -0.2e1 * t12 * t44 + 0.2e1 * t17 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, t1, -t5 * pkin(4) + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, 0, t27, t24, 0, t28, 0, 0, t42 - t49, -t48 + t58, -t6, -t13 * pkin(4) - pkin(7) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t24, 0, t28, 0, 0, -0.2e1 * t49, -0.2e1 * t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t30 + t33 * t4, t34 * t4 + t8 * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t50, 0, -t33 * t12 - t18 * t30, -t34 * t12 + t18 * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t50, 0, -pkin(7) * t30, pkin(7) * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
