% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:38
% DurationCPUTime: 0.30s
% Computational Cost: add. (461->60), mult. (1127->118), div. (0->0), fcn. (982->6), ass. (0->52)
t57 = qJD(2) + qJD(3);
t34 = sin(qJ(2));
t56 = pkin(5) + pkin(6);
t43 = qJD(2) * t56;
t23 = t34 * t43;
t37 = cos(qJ(2));
t24 = t37 * t43;
t33 = sin(qJ(3));
t36 = cos(qJ(3));
t25 = t56 * t34;
t26 = t56 * t37;
t40 = -t36 * t25 - t33 * t26;
t7 = -t40 * qJD(3) + t36 * t23 + t33 * t24;
t32 = sin(qJ(4));
t54 = pkin(2) * qJD(3);
t47 = t33 * t54;
t53 = qJD(4) * t32;
t55 = t33 * pkin(2) * t53 + t32 * t47;
t35 = cos(qJ(4));
t52 = qJD(4) * t35;
t51 = t34 * qJD(2);
t50 = t37 * qJD(2);
t49 = -0.2e1 * pkin(1) * qJD(2);
t48 = pkin(2) * t51;
t46 = t36 * t54;
t45 = pkin(3) * t53;
t44 = pkin(3) * t52;
t31 = -t37 * pkin(2) - pkin(1);
t30 = t36 * pkin(2) + pkin(3);
t42 = (-pkin(3) - t30) * qJD(4);
t21 = t33 * t34 - t36 * t37;
t22 = t33 * t37 + t36 * t34;
t15 = t35 * t21 + t32 * t22;
t16 = -t32 * t21 + t35 * t22;
t39 = t33 * t25 - t36 * t26;
t8 = t39 * qJD(3) + t33 * t23 - t36 * t24;
t38 = (-t33 * t52 + (-t32 * t36 - t33 * t35) * qJD(3)) * pkin(2);
t19 = t21 * pkin(3) + t31;
t18 = t57 * t22;
t17 = t57 * t21;
t13 = -t30 * t53 + t38;
t12 = (-qJD(4) * t30 - t46) * t35 + t55;
t11 = t18 * pkin(3) + t48;
t10 = -t21 * pkin(7) - t39;
t9 = -t22 * pkin(7) + t40;
t6 = t17 * pkin(7) + t8;
t5 = -t18 * pkin(7) - t7;
t4 = t16 * qJD(4) - t32 * t17 + t35 * t18;
t3 = -t15 * qJD(4) - t35 * t17 - t32 * t18;
t2 = -t32 * t5 + t35 * t6 + (-t10 * t35 - t32 * t9) * qJD(4);
t1 = -t32 * t6 - t35 * t5 + (t10 * t32 - t35 * t9) * qJD(4);
t14 = [0, 0, 0, 0.2e1 * t34 * t50, 0.2e1 * (-t34 ^ 2 + t37 ^ 2) * qJD(2), 0, 0, 0, t34 * t49, t37 * t49, -0.2e1 * t22 * t17, 0.2e1 * t17 * t21 - 0.2e1 * t22 * t18, 0, 0, 0, 0.2e1 * t31 * t18 + 0.2e1 * t21 * t48, -0.2e1 * t31 * t17 + 0.2e1 * t22 * t48, 0.2e1 * t16 * t3, -0.2e1 * t3 * t15 - 0.2e1 * t16 * t4, 0, 0, 0, 0.2e1 * t11 * t15 + 0.2e1 * t19 * t4, 0.2e1 * t11 * t16 + 0.2e1 * t19 * t3; 0, 0, 0, 0, 0, t50, -t51, 0, -pkin(5) * t50, pkin(5) * t51, 0, 0, -t17, -t18, 0, t8, t7, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t47, -0.2e1 * t46, 0, 0, 0, 0, 0, 0.2e1 * t13, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t8, t7, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, 0, 0, 0, 0, t32 * t42 + t38, (t42 - t46) * t35 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t45, -0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
