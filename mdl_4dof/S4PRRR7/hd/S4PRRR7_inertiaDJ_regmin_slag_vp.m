% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:49
% DurationCPUTime: 0.30s
% Computational Cost: add. (127->61), mult. (468->143), div. (0->0), fcn. (385->8), ass. (0->58)
t20 = cos(qJ(4));
t13 = t20 ^ 2;
t17 = sin(qJ(4));
t55 = t17 ^ 2 - t13;
t33 = t55 * qJD(4);
t59 = pkin(6) * t17;
t15 = sin(pkin(4));
t19 = sin(qJ(2));
t58 = t15 * t19;
t22 = cos(qJ(2));
t57 = t15 * t22;
t21 = cos(qJ(3));
t56 = t20 * t21;
t18 = sin(qJ(3));
t12 = t18 ^ 2;
t54 = -t21 ^ 2 + t12;
t16 = cos(pkin(4));
t7 = -t16 * t21 + t18 * t58;
t53 = qJD(3) * t7;
t52 = qJD(2) * t19;
t51 = qJD(3) * t18;
t50 = qJD(3) * t20;
t49 = qJD(3) * t21;
t48 = qJD(3) * t22;
t47 = qJD(4) * t17;
t46 = qJD(4) * t20;
t45 = qJD(4) * t21;
t44 = -0.2e1 * pkin(2) * qJD(3);
t43 = -0.2e1 * pkin(3) * qJD(4);
t42 = t17 * t45;
t41 = t20 * t45;
t40 = t15 * t52;
t39 = qJD(2) * t57;
t38 = t17 * t51;
t37 = t17 * t46;
t36 = t18 * t49;
t35 = t18 * t50;
t34 = t20 * t49;
t32 = t54 * qJD(3);
t31 = t18 * t34;
t30 = -t21 * pkin(3) - t18 * pkin(7);
t29 = pkin(3) * t18 - pkin(7) * t21;
t8 = t16 * t18 + t21 * t58;
t28 = -t8 * t17 - t20 * t57;
t27 = t17 * t57 - t8 * t20;
t6 = t8 * qJD(3) + t18 * t39;
t26 = t6 * t17 + t7 * t46;
t25 = -t6 * t20 + t7 * t47;
t24 = t35 + t42;
t23 = t38 - t41;
t10 = -pkin(2) + t30;
t9 = t29 * qJD(3);
t5 = t21 * t39 - t53;
t4 = t23 * pkin(6) - t10 * t47 + t20 * t9;
t3 = t24 * pkin(6) - t10 * t46 - t17 * t9;
t2 = t28 * qJD(4) + t17 * t40 + t5 * t20;
t1 = t27 * qJD(4) - t5 * t17 + t20 * t40;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t40, -t39, 0, 0, 0, 0, 0, (-t18 * t48 - t21 * t52) * t15, (t18 * t52 - t21 * t48) * t15, 0, 0, 0, 0, 0, (t17 * t53 - t1) * t21 + (t28 * qJD(3) + t26) * t18, (t7 * t50 + t2) * t21 + (t27 * qJD(3) - t25) * t18; 0, 0, 0, 0, 0.2e1 * t36, -0.2e1 * t32, 0, 0, 0, t18 * t44, t21 * t44, -0.2e1 * t12 * t37 + 0.2e1 * t13 * t36, 0.2e1 * t12 * t33 - 0.4e1 * t17 * t31, 0.2e1 * t18 * t42 + 0.2e1 * t54 * t50, -0.2e1 * t17 * t32 + 0.2e1 * t18 * t41, -0.2e1 * t36, 0.2e1 * t10 * t35 - 0.2e1 * t4 * t21 + 0.2e1 * (t12 * t46 + t17 * t36) * pkin(6), -0.2e1 * t10 * t38 - 0.2e1 * t3 * t21 + 0.2e1 * (-t12 * t47 + t31) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5, 0, 0, 0, 0, 0, t25, t26; 0, 0, 0, 0, 0, 0, t49, -t51, 0, -pkin(6) * t49, pkin(6) * t51, t17 * t34 - t18 * t33, -0.4e1 * t18 * t37 - t55 * t49, t23, t24, 0, (pkin(7) * t56 + (-pkin(3) * t20 + t59) * t18) * qJD(4) + (-pkin(6) * t56 + t30 * t17) * qJD(3), (pkin(6) * t18 * t20 + t29 * t17) * qJD(4) + (t30 * t20 + t21 * t59) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t37, -0.2e1 * t33, 0, 0, 0, t17 * t43, t20 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t47 + t34, -t17 * t49 - t18 * t46, t51, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t47, 0, -pkin(7) * t46, pkin(7) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
