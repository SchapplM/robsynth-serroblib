% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:38
% EndTime: 2021-01-15 10:56:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (325->69), mult. (842->159), div. (0->0), fcn. (710->6), ass. (0->56)
t60 = 2 * qJD(4);
t29 = sin(qJ(4));
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t53 = -qJ(3) - pkin(5);
t42 = qJD(2) * t53;
t11 = t32 * qJD(3) + t30 * t42;
t27 = sin(pkin(7));
t28 = cos(pkin(7));
t33 = -t30 * qJD(3) + t32 * t42;
t3 = t27 * t11 - t28 * t33;
t59 = t3 * t29;
t31 = cos(qJ(4));
t58 = t3 * t31;
t16 = t27 * t32 + t28 * t30;
t12 = t16 * qJD(2);
t57 = t29 * t12;
t48 = t32 * qJD(2);
t49 = t30 * qJD(2);
t13 = -t27 * t49 + t28 * t48;
t56 = t29 * t13;
t55 = t31 * t12;
t54 = t31 * t13;
t26 = t31 ^ 2;
t52 = t29 ^ 2 - t26;
t51 = qJD(4) * t29;
t50 = qJD(4) * t31;
t47 = -0.2e1 * pkin(1) * qJD(2);
t22 = -t28 * pkin(2) - pkin(3);
t46 = t22 * t60;
t24 = pkin(2) * t49;
t45 = t29 * t50;
t23 = -t32 * pkin(2) - pkin(1);
t44 = t53 * t30;
t43 = -0.4e1 * t16 * t29 * t31;
t41 = t52 * qJD(4);
t15 = t27 * t30 - t28 * t32;
t7 = t15 * pkin(3) - t16 * pkin(6) + t23;
t18 = t53 * t32;
t9 = -t28 * t18 + t27 * t44;
t40 = t29 * t9 - t31 * t7;
t39 = t29 * t7 + t31 * t9;
t21 = t27 * pkin(2) + pkin(6);
t38 = -t12 * t21 + t13 * t22;
t37 = t15 * t21 - t16 * t22;
t36 = t15 * t50 + t57;
t35 = -t16 * t50 - t56;
t34 = -t16 * t51 + t54;
t14 = t16 ^ 2;
t8 = -t27 * t18 - t28 * t44;
t6 = -t15 * t51 + t55;
t5 = t12 * pkin(3) - t13 * pkin(6) + t24;
t4 = t28 * t11 + t27 * t33;
t2 = -qJD(4) * t39 - t29 * t4 + t31 * t5;
t1 = qJD(4) * t40 - t29 * t5 - t31 * t4;
t10 = [0, 0, 0, 0.2e1 * t30 * t48, 0.2e1 * (-t30 ^ 2 + t32 ^ 2) * qJD(2), 0, 0, 0, t30 * t47, t32 * t47, 0.2e1 * t23 * t12 + 0.2e1 * t15 * t24, 0.2e1 * t23 * t13 + 0.2e1 * t16 * t24, -0.2e1 * t9 * t12 + 0.2e1 * t8 * t13 - 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16, 0.2e1 * t23 * t24 + 0.2e1 * t8 * t3 + 0.2e1 * t9 * t4, 0.2e1 * t26 * t16 * t13 - 0.2e1 * t14 * t45, t52 * t14 * t60 + t13 * t43, 0.2e1 * t15 * t34 + 0.2e1 * t16 * t55, 0.2e1 * t15 * t35 - 0.2e1 * t16 * t57, 0.2e1 * t15 * t12, 0.2e1 * t2 * t15 - 0.2e1 * t40 * t12 + 0.2e1 * t8 * t56 + 0.2e1 * (t8 * t50 + t59) * t16, 0.2e1 * t1 * t15 - 0.2e1 * t39 * t12 + 0.2e1 * t8 * t54 + 0.2e1 * (-t8 * t51 + t58) * t16; 0, 0, 0, 0, 0, t48, -t49, 0, -pkin(5) * t48, pkin(5) * t49, -t3, -t4, (-t12 * t27 - t13 * t28) * pkin(2), (t27 * t4 - t28 * t3) * pkin(2), -t16 * t41 + t29 * t54, qJD(4) * t43 - t52 * t13, t36, t6, 0, -t58 + t38 * t29 + (t29 * t8 - t31 * t37) * qJD(4), t59 + t38 * t31 + (t29 * t37 + t31 * t8) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t45, -0.2e1 * t41, 0, 0, 0, t29 * t46, t31 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, t24, 0, 0, 0, 0, 0, t6, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, t12, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t51, 0, -t21 * t50, t21 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
