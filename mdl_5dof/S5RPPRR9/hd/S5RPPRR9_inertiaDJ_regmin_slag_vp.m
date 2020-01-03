% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:49
% DurationCPUTime: 0.35s
% Computational Cost: add. (197->64), mult. (481->147), div. (0->0), fcn. (367->6), ass. (0->56)
t24 = cos(qJ(5));
t18 = t24 ^ 2;
t22 = sin(qJ(5));
t59 = t22 ^ 2 - t18;
t37 = qJD(5) * t59;
t23 = sin(qJ(4));
t17 = t23 ^ 2;
t25 = cos(qJ(4));
t36 = (-t25 ^ 2 + t17) * qJD(4);
t62 = 2 * qJD(2);
t20 = sin(pkin(8));
t61 = t20 * t25;
t21 = cos(pkin(8));
t26 = -pkin(1) - pkin(2);
t60 = t21 * qJ(2) + t20 * t26;
t57 = qJD(5) * t22;
t56 = qJD(5) * t23;
t55 = qJD(5) * t24;
t54 = qJD(5) * t25;
t53 = t20 * qJD(2);
t52 = t21 * qJD(2);
t51 = t23 * qJD(4);
t50 = t25 * qJD(4);
t49 = -0.2e1 * pkin(4) * qJD(5);
t48 = pkin(7) * t55;
t47 = t22 * t54;
t46 = t24 * t54;
t45 = t24 * t52;
t44 = t22 * t51;
t43 = t22 * t55;
t42 = t23 * t50;
t41 = t20 * t51;
t40 = t24 * t51;
t39 = t24 * t50;
t38 = -t20 * qJ(2) + t21 * t26;
t35 = t23 * t39;
t11 = pkin(3) - t38;
t12 = -pkin(6) + t60;
t32 = -pkin(4) * t23 + pkin(7) * t25;
t34 = -qJD(4) * t32 + t12 * t54 - t53;
t33 = pkin(4) * t25 + pkin(7) * t23;
t31 = t12 * t51 - t25 * t52;
t30 = -t12 * t50 - t23 * t52;
t7 = t22 * t56 - t39;
t9 = t22 * t50 + t23 * t55;
t29 = t17 * t57 - t35;
t28 = -t17 * t55 - t22 * t42;
t5 = t11 + t33;
t27 = -qJD(5) * t5 + t31;
t10 = -t44 + t46;
t8 = -t40 - t47;
t4 = t22 * t41 + (t21 * t22 - t24 * t61) * qJD(5);
t3 = t20 * t40 + (t21 * t24 + t22 * t61) * qJD(5);
t2 = t22 * t27 - t24 * t34;
t1 = t22 * t34 + t24 * t27;
t6 = [0, 0, 0, 0, t62, qJ(2) * t62, 0.2e1 * t53, 0.2e1 * t52, (-t20 * t38 + t21 * t60) * t62, 0.2e1 * t42, -0.2e1 * t36, 0, 0, 0, -0.2e1 * t11 * t51 + 0.2e1 * t25 * t53, -0.2e1 * t11 * t50 - 0.2e1 * t23 * t53, -0.2e1 * t17 * t43 + 0.2e1 * t18 * t42, 0.2e1 * t17 * t37 - 0.4e1 * t22 * t35, 0.2e1 * t23 * t47 + 0.2e1 * t24 * t36, -0.2e1 * t22 * t36 + 0.2e1 * t23 * t46, -0.2e1 * t42, -0.2e1 * t17 * t22 * t52 + 0.2e1 * t12 * t28 + 0.2e1 * t2 * t25 - 0.2e1 * t40 * t5, 0.2e1 * t1 * t25 + 0.2e1 * t12 * t29 - 0.2e1 * t17 * t45 + 0.2e1 * t44 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t51, t21 * t50, 0, 0, 0, 0, 0, t20 * t28 + t21 * t40 + t4 * t25, t20 * t29 - t21 * t44 + t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t51, 0, t30, t31, -t22 * t39 + t23 * t37, 0.4e1 * t23 * t43 + t50 * t59, t10, t8, 0, (-t48 + (pkin(4) * t22 - t12 * t24) * qJD(4)) * t25 + (pkin(7) * qJD(4) * t22 - t45 + (pkin(4) * t24 + t12 * t22) * qJD(5)) * t23, (qJD(4) * t33 + t12 * t56) * t24 + (qJD(5) * t32 - t30) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t50, t41, 0, 0, 0, 0, 0, t7 * t20, t9 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, 0, 0, 0, 0, t8, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, -0.2e1 * t37, 0, 0, 0, t22 * t49, t24 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, -t51, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t57, 0, -t48, pkin(7) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
