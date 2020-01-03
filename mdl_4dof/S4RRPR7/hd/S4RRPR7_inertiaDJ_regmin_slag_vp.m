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
% MMD_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:45
% DurationCPUTime: 0.28s
% Computational Cost: add. (305->65), mult. (784->151), div. (0->0), fcn. (666->6), ass. (0->56)
t60 = 2 * qJD(4);
t28 = sin(qJ(4));
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t53 = -qJ(3) - pkin(5);
t41 = qJD(2) * t53;
t11 = t31 * qJD(3) + t29 * t41;
t26 = sin(pkin(7));
t27 = cos(pkin(7));
t32 = -t29 * qJD(3) + t31 * t41;
t3 = t26 * t11 - t27 * t32;
t59 = t3 * t28;
t30 = cos(qJ(4));
t58 = t3 * t30;
t16 = t26 * t31 + t27 * t29;
t12 = t16 * qJD(2);
t57 = t28 * t12;
t50 = qJD(2) * t31;
t51 = qJD(2) * t29;
t13 = -t26 * t51 + t27 * t50;
t56 = t28 * t13;
t55 = t30 * t12;
t54 = t30 * t13;
t25 = t30 ^ 2;
t52 = t28 ^ 2 - t25;
t49 = qJD(4) * t28;
t48 = qJD(4) * t30;
t47 = -0.2e1 * pkin(1) * qJD(2);
t22 = -t27 * pkin(2) - pkin(3);
t46 = t22 * t60;
t23 = pkin(2) * t51;
t45 = t28 * t48;
t44 = -t31 * pkin(2) - pkin(1);
t43 = t53 * t29;
t42 = -0.4e1 * t16 * t28 * t30;
t40 = t52 * qJD(4);
t15 = t26 * t29 - t27 * t31;
t7 = t15 * pkin(3) - t16 * pkin(6) + t44;
t18 = t53 * t31;
t9 = -t27 * t18 + t26 * t43;
t39 = t28 * t9 - t30 * t7;
t38 = t28 * t7 + t30 * t9;
t21 = t26 * pkin(2) + pkin(6);
t37 = -t12 * t21 + t13 * t22;
t36 = t15 * t21 - t16 * t22;
t35 = t15 * t48 + t57;
t34 = -t16 * t48 - t56;
t33 = -t16 * t49 + t54;
t14 = t16 ^ 2;
t8 = -t26 * t18 - t27 * t43;
t6 = -t15 * t49 + t55;
t5 = t12 * pkin(3) - t13 * pkin(6) + t23;
t4 = t27 * t11 + t26 * t32;
t2 = -t38 * qJD(4) - t28 * t4 + t30 * t5;
t1 = t39 * qJD(4) - t28 * t5 - t30 * t4;
t10 = [0, 0, 0, 0.2e1 * t29 * t50, 0.2e1 * (-t29 ^ 2 + t31 ^ 2) * qJD(2), 0, 0, 0, t29 * t47, t31 * t47, -0.2e1 * t9 * t12 + 0.2e1 * t8 * t13 - 0.2e1 * t4 * t15 + 0.2e1 * t3 * t16, 0.2e1 * t44 * t23 + 0.2e1 * t8 * t3 + 0.2e1 * t9 * t4, 0.2e1 * t25 * t16 * t13 - 0.2e1 * t14 * t45, t52 * t14 * t60 + t13 * t42, 0.2e1 * t33 * t15 + 0.2e1 * t16 * t55, 0.2e1 * t34 * t15 - 0.2e1 * t16 * t57, 0.2e1 * t15 * t12, 0.2e1 * t2 * t15 - 0.2e1 * t39 * t12 + 0.2e1 * t8 * t56 + 0.2e1 * (t8 * t48 + t59) * t16, 0.2e1 * t1 * t15 - 0.2e1 * t38 * t12 + 0.2e1 * t8 * t54 + 0.2e1 * (-t8 * t49 + t58) * t16; 0, 0, 0, 0, 0, t50, -t51, 0, -pkin(5) * t50, pkin(5) * t51, (-t12 * t26 - t13 * t27) * pkin(2), (t26 * t4 - t27 * t3) * pkin(2), -t16 * t40 + t28 * t54, qJD(4) * t42 - t52 * t13, t35, t6, 0, -t58 + t37 * t28 + (t28 * t8 - t36 * t30) * qJD(4), t59 + t37 * t30 + (t36 * t28 + t30 * t8) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t45, -0.2e1 * t40, 0, 0, 0, t28 * t46, t30 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t6, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t34, t12, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, 0, -t21 * t48, t21 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
