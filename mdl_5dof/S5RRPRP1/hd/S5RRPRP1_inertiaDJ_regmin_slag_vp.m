% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:12
% EndTime: 2020-01-03 11:59:12
% DurationCPUTime: 0.23s
% Computational Cost: add. (277->69), mult. (631->115), div. (0->0), fcn. (445->6), ass. (0->52)
t37 = cos(qJ(4));
t57 = t37 * pkin(4);
t38 = cos(qJ(2));
t28 = t38 * pkin(1) + pkin(2);
t34 = cos(pkin(8));
t33 = sin(pkin(8));
t36 = sin(qJ(2));
t55 = t33 * t36;
t41 = -pkin(1) * t55 + t34 * t28;
t14 = -pkin(3) - t41;
t51 = pkin(1) * qJD(2);
t54 = t34 * t36;
t16 = (t33 * t38 + t54) * t51;
t31 = qJD(4) * t37;
t35 = sin(qJ(4));
t56 = t14 * t31 + t16 * t35;
t17 = (t34 * t38 - t55) * t51;
t53 = t37 * t17;
t52 = pkin(1) * t54 + t33 * t28;
t15 = pkin(7) + t52;
t50 = -qJ(5) - t15;
t26 = t33 * pkin(2) + pkin(7);
t49 = -qJ(5) - t26;
t48 = qJD(4) * t35;
t47 = t36 * t51;
t46 = t38 * t51;
t29 = pkin(4) * t48;
t45 = pkin(4) * t31;
t27 = -t34 * pkin(2) - pkin(3);
t44 = t27 * t48;
t43 = t27 * t31;
t42 = t14 * t48 - t16 * t37;
t40 = qJD(4) * t50;
t39 = qJD(4) * t49;
t32 = t37 * qJ(5);
t30 = t37 * qJD(5);
t23 = 0.2e1 * t35 * t31;
t21 = t27 - t57;
t20 = 0.2e1 * (-t35 ^ 2 + t37 ^ 2) * qJD(4);
t19 = t37 * t26 + t32;
t18 = t49 * t35;
t12 = -t35 * qJD(5) + t37 * t39;
t11 = t35 * t39 + t30;
t10 = t14 - t57;
t7 = t29 + t16;
t6 = t11 * t37;
t5 = t37 * t15 + t32;
t4 = t50 * t35;
t3 = (-qJD(5) - t17) * t35 + t37 * t40;
t2 = t35 * t40 + t30 + t53;
t1 = t2 * t37;
t8 = [0, 0, 0, 0, -0.2e1 * t47, -0.2e1 * t46, -0.2e1 * t41 * t16 + 0.2e1 * t52 * t17, t23, t20, 0, 0, 0, 0.2e1 * t42, 0.2e1 * t56, -0.2e1 * t3 * t35 + 0.2e1 * t1 + 0.2e1 * (-t35 * t5 - t37 * t4) * qJD(4), 0.2e1 * t10 * t7 + 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, -t47, -t46, (-t16 * t34 + t17 * t33) * pkin(2), t23, t20, 0, 0, 0, t42 + t44, t43 + t56, t1 + t6 + (-t12 - t3) * t35 + ((-t18 - t4) * t37 + (-t19 - t5) * t35) * qJD(4), t10 * t29 + t5 * t11 + t4 * t12 + t3 * t18 + t2 * t19 + t7 * t21; 0, 0, 0, 0, 0, 0, 0, t23, t20, 0, 0, 0, 0.2e1 * t44, 0.2e1 * t43, -0.2e1 * t12 * t35 + 0.2e1 * t6 + 0.2e1 * (-t18 * t37 - t19 * t35) * qJD(4), 0.2e1 * t19 * t11 + 0.2e1 * t18 * t12 + 0.2e1 * t21 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t35 + t3 * t37 + (-t35 * t4 + t37 * t5) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t35 + t12 * t37 + (-t18 * t35 + t19 * t37) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t48, 0, -t15 * t31 - t35 * t17, t15 * t48 - t53, -t45, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t48, 0, -t26 * t31, t26 * t48, -t45, t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t31, 0, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
