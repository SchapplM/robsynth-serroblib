% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:47
% DurationCPUTime: 0.46s
% Computational Cost: add. (467->63), mult. (1082->138), div. (0->0), fcn. (1027->8), ass. (0->43)
t33 = sin(pkin(9));
t53 = sin(qJ(5));
t43 = qJD(5) * t53;
t35 = cos(pkin(9));
t37 = cos(qJ(5));
t50 = t37 * t35;
t19 = -qJD(5) * t50 + t33 * t43;
t22 = -t53 * t33 + t50;
t36 = cos(pkin(8));
t56 = 0.2e1 * t36;
t34 = sin(pkin(8));
t55 = pkin(6) * t34;
t14 = t22 * t34;
t40 = -t37 * t33 - t53 * t35;
t13 = t40 * t34;
t9 = qJD(5) * t13;
t54 = t14 * t9;
t10 = t19 * t34;
t52 = t13 * t10;
t29 = sin(pkin(7)) * pkin(1) + qJ(3);
t51 = t29 * t36;
t21 = -cos(pkin(7)) * pkin(1) - pkin(2) - t34 * qJ(4) - t36 * pkin(3);
t49 = t33 * t21 + t35 * t51;
t48 = qJD(3) * t36;
t47 = qJD(4) * t34;
t31 = t34 ^ 2;
t46 = t31 * qJD(3);
t30 = t34 * qJD(3);
t17 = -t33 * t48 - t35 * t47;
t18 = -t33 * t47 + t35 * t48;
t41 = t17 * t35 + t18 * t33;
t39 = (-t29 * t33 - pkin(4)) * t36 + (t21 - t55) * t35;
t38 = t37 * t39;
t8 = -t33 * t55 + t49;
t4 = t37 * t8 + t53 * t39;
t32 = t36 ^ 2;
t25 = t29 * t46;
t20 = t40 * qJD(5);
t16 = (pkin(4) * t33 + t29) * t34;
t3 = -t53 * t8 + t38;
t2 = -t4 * qJD(5) + t37 * t17 - t53 * t18;
t1 = -qJD(5) * t38 - t53 * t17 - t37 * t18 + t8 * t43;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t31 + t32) * qJD(3), 0.2e1 * t32 * t29 * qJD(3) + 0.2e1 * t25, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t36 + 0.2e1 * t33 * t46, 0.2e1 * t18 * t36 + 0.2e1 * t35 * t46, -0.2e1 * t41 * t34, 0.2e1 * t49 * t18 + 0.2e1 * (t35 * t21 - t33 * t51) * t17 + 0.2e1 * t25, 0.2e1 * t54, 0.2e1 * t14 * t10 + 0.2e1 * t13 * t9, -t9 * t56, 0.2e1 * t52, -t10 * t56, 0, -0.2e1 * t16 * t10 - 0.2e1 * t13 * t30 - 0.2e1 * t2 * t36, -0.2e1 * t1 * t36 + 0.2e1 * t14 * t30 + 0.2e1 * t16 * t9, -0.2e1 * t1 * t13 + 0.2e1 * t4 * t10 - 0.2e1 * t2 * t14 - 0.2e1 * t3 * t9, -0.2e1 * t4 * t1 + 0.2e1 * t16 * t30 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t17 * t33 + t18 * t35 - t48) * t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t14 + t3 * t10 + t2 * t13 - t36 * t30 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52 + 0.2e1 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, -t20 * t36, -t19 * t36, -t10 * t40 - t19 * t13 - t20 * t14 - t22 * t9, t1 * t40 - t4 * t19 + t2 * t22 + t3 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t22 + t13 * t20 - t14 * t19 - t40 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19 * t40 + 0.2e1 * t22 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, -t10, t9, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t10, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
