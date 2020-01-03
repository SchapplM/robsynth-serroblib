% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:49
% EndTime: 2020-01-03 11:25:52
% DurationCPUTime: 0.43s
% Computational Cost: add. (333->67), mult. (737->121), div. (0->0), fcn. (578->6), ass. (0->50)
t30 = sin(qJ(4));
t31 = cos(qJ(4));
t29 = cos(pkin(8));
t48 = qJD(4) * t30;
t23 = t29 * t48;
t25 = sin(pkin(7)) * pkin(1) + qJ(3);
t28 = sin(pkin(8));
t33 = -cos(pkin(7)) * pkin(1) - t29 * pkin(3) - pkin(2);
t32 = -t28 * pkin(6) + t33;
t12 = t31 * t32;
t45 = t29 * qJD(3);
t52 = -qJD(4) * t12 - t31 * t45;
t4 = t25 * t23 + t52;
t44 = t30 * qJD(3);
t40 = t29 * t44;
t15 = t31 * t29 * t25;
t8 = t30 * t32 + t15;
t5 = -qJD(4) * t8 - t40;
t51 = t25 * t30;
t43 = t29 * t51;
t7 = t12 - t43;
t56 = t4 * t30 - t5 * t31 + (t30 * t7 - t31 * t8) * qJD(4);
t50 = t28 * qJ(5);
t42 = t31 * t50;
t46 = t28 * qJD(5);
t1 = t30 * t46 + (t42 + t43) * qJD(4) + t52;
t2 = -t40 - t31 * t46 + (-t15 + ((pkin(6) + qJ(5)) * t28 - t33) * t30) * qJD(4);
t3 = -t42 + t12 + (-pkin(4) - t51) * t29;
t6 = -t30 * t50 + t8;
t55 = t1 * t30 - t2 * t31 + (t3 * t30 - t31 * t6) * qJD(4);
t54 = 0.2e1 * t28;
t53 = 0.2e1 * qJD(4);
t49 = qJD(3) * t25;
t47 = qJD(4) * t31;
t22 = t28 * t48;
t41 = t28 * t47;
t39 = t28 * t29 * t53;
t26 = t28 ^ 2;
t38 = t26 * t30 * t47;
t27 = t29 ^ 2;
t24 = t29 * t47;
t20 = -0.2e1 * t38;
t19 = 0.2e1 * t38;
t18 = t31 * t39;
t17 = t30 * t39;
t16 = t26 * t49;
t14 = (pkin(4) * t47 + qJD(3)) * t28;
t13 = (pkin(4) * t30 + t25) * t28;
t10 = (t30 ^ 2 - t31 ^ 2) * t26 * t53;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t26 + t27) * qJD(3), 0.2e1 * t27 * t49 + 0.2e1 * t16, t20, t10, t17, t19, t18, 0, -0.2e1 * t5 * t29 + 0.2e1 * (t25 * t47 + t44) * t26, -0.2e1 * t4 * t29 + 0.2e1 * (qJD(3) * t31 - t25 * t48) * t26, t56 * t54, -0.2e1 * t8 * t4 + 0.2e1 * t7 * t5 + 0.2e1 * t16, t20, t10, t17, t19, t18, 0, -0.2e1 * t2 * t29 + 0.2e1 * (t13 * t47 + t14 * t30) * t28, -0.2e1 * t1 * t29 + 0.2e1 * (-t13 * t48 + t14 * t31) * t28, t55 * t54, -0.2e1 * t6 * t1 + 0.2e1 * t13 * t14 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t45 - t30 * t5 - t31 * t4 + (-t30 * t8 - t31 * t7) * qJD(4)) * t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t29 + (-t1 * t31 - t2 * t30 + (-t3 * t31 - t30 * t6) * qJD(4)) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t24, 0, -t56, 0, 0, 0, 0, 0, 0, t23, t24, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, -t41, 0, t5, t4, 0, 0, 0, 0, -t22, 0, -t41, 0, t2, t1, pkin(4) * t22, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t22, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t22, 0, -pkin(4) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, -pkin(4) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t22, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
