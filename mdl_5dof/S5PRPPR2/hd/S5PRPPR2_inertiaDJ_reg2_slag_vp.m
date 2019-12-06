% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPPR2
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:42
% EndTime: 2019-12-05 15:24:44
% DurationCPUTime: 0.40s
% Computational Cost: add. (335->57), mult. (871->114), div. (0->0), fcn. (869->8), ass. (0->44)
t31 = sin(pkin(8));
t33 = cos(pkin(8));
t35 = sin(qJ(2));
t36 = cos(qJ(2));
t21 = t31 * t36 + t33 * t35;
t32 = cos(pkin(9));
t30 = sin(pkin(9));
t34 = sin(qJ(5));
t50 = t34 * t30;
t53 = cos(qJ(5));
t37 = t53 * t32 - t50;
t6 = t37 * t21;
t26 = pkin(2) * t31 + qJ(4);
t54 = pkin(6) + t26;
t14 = t21 * qJD(2);
t19 = t31 * t35 - t33 * t36;
t10 = t19 * t14;
t45 = t53 * t30;
t49 = t34 * t32;
t22 = t45 + t49;
t17 = t22 * qJD(5);
t52 = t37 * t17;
t42 = qJD(5) * t53;
t16 = qJD(5) * t50 - t32 * t42;
t51 = t22 * t16;
t48 = t30 ^ 2 + t32 ^ 2;
t47 = -pkin(2) * t33 - pkin(3);
t46 = t34 * t54;
t15 = t19 * qJD(2);
t43 = t48 * t15;
t41 = t53 * qJD(4);
t40 = t48 * qJD(4);
t39 = 0.2e1 * t40;
t38 = t54 * t45;
t23 = -pkin(4) * t32 + t47;
t18 = t54 * t32;
t8 = t53 * t18 - t30 * t46;
t7 = -t34 * t18 - t38;
t5 = t22 * t21;
t4 = -t18 * t42 - qJD(4) * t49 + (qJD(5) * t46 - t41) * t30;
t3 = qJD(5) * t38 - t32 * t41 + (qJD(4) * t30 + qJD(5) * t18) * t34;
t2 = -qJD(5) * t6 + t22 * t15;
t1 = t37 * t15 + t21 * t17;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t15 * t21 + 0.2e1 * t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21 * t43 + 0.2e1 * t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t1 * t6 - 0.2e1 * t2 * t5 + 0.2e1 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * qJD(2), -t36 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t14, t15, 0, (-t14 * t33 - t15 * t31) * pkin(2), 0, 0, 0, 0, 0, 0, -t14 * t32, t14 * t30, -t43, t14 * t47 + t21 * t40 - t26 * t43, 0, 0, 0, 0, 0, 0, -t14 * t37 + t17 * t19, t14 * t22 - t16 * t19, -t1 * t37 - t16 * t5 - t17 * t6 - t2 * t22, -t1 * t8 + t14 * t23 + t2 * t7 - t3 * t6 - t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t26 * t39, -0.2e1 * t51, -0.2e1 * t16 * t37 - 0.2e1 * t22 * t17, 0, -0.2e1 * t52, 0, 0, 0.2e1 * t23 * t17, -0.2e1 * t23 * t16, 0.2e1 * t16 * t7 - 0.2e1 * t17 * t8 - 0.2e1 * t22 * t4 - 0.2e1 * t3 * t37, -0.2e1 * t3 * t8 + 0.2e1 * t4 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t22 - t16 * t6 + t17 * t5 + t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t8 - t17 * t7 - t22 * t3 + t37 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t51 - 0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, -t17, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
