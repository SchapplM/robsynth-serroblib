% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:35
% EndTime: 2020-01-03 12:00:36
% DurationCPUTime: 0.20s
% Computational Cost: add. (296->56), mult. (758->87), div. (0->0), fcn. (550->8), ass. (0->51)
t35 = cos(qJ(5));
t29 = qJD(5) * t35;
t37 = cos(qJ(2));
t28 = t37 * pkin(1) + pkin(2);
t31 = cos(pkin(9));
t30 = sin(pkin(9));
t34 = sin(qJ(2));
t56 = t30 * t34;
t40 = -pkin(1) * t56 + t31 * t28;
t17 = pkin(3) + t40;
t55 = t31 * t34;
t22 = pkin(1) * t55 + t30 * t28;
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t39 = t33 * t17 + t36 * t22;
t54 = pkin(1) * qJD(2);
t20 = (-t30 * t37 - t55) * t54;
t21 = (t31 * t37 - t56) * t54;
t41 = t36 * t20 - t33 * t21;
t3 = t39 * qJD(4) - t41;
t32 = sin(qJ(5));
t6 = -t36 * t17 + t33 * t22 - pkin(4);
t59 = t6 * t29 + t3 * t32;
t58 = pkin(2) * t30;
t27 = t31 * pkin(2) + pkin(3);
t38 = t33 * t27 + t36 * t58;
t16 = t38 * qJD(4);
t50 = t33 * t58;
t18 = -t36 * t27 - pkin(4) + t50;
t57 = t16 * t32 + t18 * t29;
t53 = qJD(4) * t33;
t52 = qJD(4) * t36;
t51 = qJD(5) * t32;
t49 = -t17 * t52 - t33 * t20 - t36 * t21;
t48 = pkin(4) * t51;
t47 = pkin(4) * t29;
t46 = t34 * t54;
t45 = t37 * t54;
t4 = t6 * t51;
t44 = -t3 * t35 + t4;
t43 = t22 + t58;
t9 = t18 * t51;
t42 = -t16 * t35 + t9;
t26 = 0.2e1 * t32 * t29;
t24 = t27 * t52;
t23 = 0.2e1 * (-t32 ^ 2 + t35 ^ 2) * qJD(5);
t19 = pkin(8) + t38;
t15 = qJD(4) * t50 - t24;
t7 = pkin(8) + t39;
t2 = t22 * t53 + t49;
t1 = [0, 0, 0, 0, -0.2e1 * t46, -0.2e1 * t45, 0.2e1 * t40 * t20 + 0.2e1 * t22 * t21, 0, -0.2e1 * t3, 0.2e1 * t2, t26, t23, 0, 0, 0, 0.2e1 * t44, 0.2e1 * t59; 0, 0, 0, 0, -t46, -t45, (t20 * t31 + t21 * t30) * pkin(2), 0, (-t43 * t36 + (-t17 - t27) * t33) * qJD(4) + t41, t43 * t53 - t24 + t49, t26, t23, 0, 0, 0, t4 + t9 + (-t16 - t3) * t35, t57 + t59; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t16, 0.2e1 * t15, t26, t23, 0, 0, 0, 0.2e1 * t42, 0.2e1 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, t26, t23, 0, 0, 0, t44 - t48, -t47 + t59; 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, t26, t23, 0, 0, 0, t42 - t48, -t47 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t23, 0, 0, 0, -0.2e1 * t48, -0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t51, 0, t32 * t2 - t7 * t29, t35 * t2 + t7 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t51, 0, t32 * t15 - t19 * t29, t35 * t15 + t19 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t51, 0, -pkin(8) * t29, pkin(8) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
