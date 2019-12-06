% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:10
% EndTime: 2019-12-05 15:22:12
% DurationCPUTime: 0.41s
% Computational Cost: add. (377->61), mult. (992->136), div. (0->0), fcn. (937->6), ass. (0->43)
t32 = sin(pkin(9));
t53 = sin(qJ(5));
t42 = qJD(5) * t53;
t34 = cos(pkin(9));
t36 = cos(qJ(5));
t51 = t36 * t34;
t17 = -qJD(5) * t51 + t32 * t42;
t20 = -t53 * t32 + t51;
t35 = cos(pkin(8));
t56 = 0.2e1 * t35;
t33 = sin(pkin(8));
t55 = pkin(6) * t33;
t14 = t20 * t33;
t39 = -t36 * t32 - t53 * t34;
t13 = t39 * t33;
t9 = qJD(5) * t13;
t54 = t14 * t9;
t10 = t17 * t33;
t52 = t13 * t10;
t23 = -t35 * pkin(3) - t33 * qJ(4) - pkin(2);
t49 = qJ(3) * t35;
t50 = t32 * t23 + t34 * t49;
t48 = qJD(3) * t35;
t47 = qJD(4) * t33;
t30 = t33 ^ 2;
t46 = t30 * qJD(3);
t29 = t33 * qJD(3);
t45 = qJ(3) * qJD(3);
t15 = -t32 * t48 - t34 * t47;
t16 = -t32 * t47 + t34 * t48;
t40 = t15 * t34 + t16 * t32;
t38 = (-qJ(3) * t32 - pkin(4)) * t35 + (t23 - t55) * t34;
t37 = t36 * t38;
t8 = -t32 * t55 + t50;
t4 = t36 * t8 + t53 * t38;
t31 = t35 ^ 2;
t28 = t30 * t45;
t22 = (pkin(4) * t32 + qJ(3)) * t33;
t18 = t39 * qJD(5);
t3 = -t53 * t8 + t37;
t2 = -t4 * qJD(5) + t36 * t15 - t53 * t16;
t1 = -qJD(5) * t37 - t53 * t15 - t36 * t16 + t8 * t42;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52 + 0.2e1 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t15 * t32 + t16 * t34 - t48) * t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t1 + t10 * t3 + t13 * t2 - t35 * t29 + t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t30 + t31) * qJD(3), 0.2e1 * t31 * t45 + 0.2e1 * t28, 0, 0, 0, 0, 0, 0, -0.2e1 * t15 * t35 + 0.2e1 * t32 * t46, 0.2e1 * t16 * t35 + 0.2e1 * t34 * t46, -0.2e1 * t40 * t33, 0.2e1 * t50 * t16 + 0.2e1 * (t34 * t23 - t32 * t49) * t15 + 0.2e1 * t28, 0.2e1 * t54, 0.2e1 * t14 * t10 + 0.2e1 * t13 * t9, -t9 * t56, 0.2e1 * t52, -t10 * t56, 0, -0.2e1 * t22 * t10 - 0.2e1 * t13 * t29 - 0.2e1 * t2 * t35, -0.2e1 * t1 * t35 + 0.2e1 * t14 * t29 + 0.2e1 * t22 * t9, -0.2e1 * t1 * t13 + 0.2e1 * t4 * t10 - 0.2e1 * t2 * t14 - 0.2e1 * t3 * t9, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t22 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t20 + t13 * t18 - t14 * t17 - t39 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, -t18 * t35, -t17 * t35, -t10 * t39 - t17 * t13 - t18 * t14 - t20 * t9, t1 * t39 - t4 * t17 + t3 * t18 + t2 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t17 * t39 + 0.2e1 * t20 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, -t10, t9, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t10, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
