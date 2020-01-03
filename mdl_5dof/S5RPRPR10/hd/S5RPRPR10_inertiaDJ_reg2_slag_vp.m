% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:07
% DurationCPUTime: 0.49s
% Computational Cost: add. (452->63), mult. (861->117), div. (0->0), fcn. (638->6), ass. (0->51)
t56 = 2 * qJD(2);
t55 = -pkin(1) - pkin(2);
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t18 = t35 * t38 - t36 * t40;
t31 = t40 * t55;
t49 = qJD(3) * t38;
t12 = qJ(2) * t49 - t40 * qJD(2) - qJD(3) * t31;
t23 = t40 * qJ(2) + t38 * t55;
t13 = t38 * qJD(2) + t23 * qJD(3);
t2 = -t35 * t12 + t36 * t13;
t54 = t2 * t18;
t37 = sin(qJ(5));
t53 = t2 * t37;
t39 = cos(qJ(5));
t52 = t2 * t39;
t19 = t35 * t40 + t36 * t38;
t14 = t19 * qJD(3);
t51 = t18 * t14;
t22 = -t38 * qJ(2) + t31;
t21 = -pkin(3) + t22;
t10 = t35 * t21 + t36 * t23;
t33 = t37 ^ 2;
t34 = t39 ^ 2;
t50 = t33 + t34;
t48 = qJD(3) * t40;
t47 = t37 * qJD(5);
t46 = t39 * qJD(5);
t30 = -t36 * pkin(3) - pkin(4);
t45 = 0.2e1 * qJD(5) * t30;
t44 = t37 * t46;
t3 = -t36 * t12 - t35 * t13;
t1 = t50 * t3;
t15 = t18 * qJD(3);
t4 = t50 * t15;
t43 = t50 * t19;
t29 = t35 * pkin(3) + pkin(7);
t42 = t50 * t29;
t9 = t36 * t21 - t35 * t23;
t7 = pkin(4) - t9;
t41 = qJD(5) * (-t30 + t7);
t26 = -0.2e1 * t44;
t25 = 0.2e1 * t44;
t24 = (-t33 + t34) * qJD(5);
t20 = 0.2e1 * t24;
t8 = -pkin(7) + t10;
t6 = t14 * t39 - t18 * t47;
t5 = t14 * t37 + t18 * t46;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, qJ(2) * t56, 0, 0, 0, 0, 0, 0, 0.2e1 * t13, -0.2e1 * t12, 0, -0.2e1 * t23 * t12 - 0.2e1 * t22 * t13, 0, 0, 0, 0, 0, 0, 0.2e1 * t2, 0.2e1 * t3, 0, 0.2e1 * t10 * t3 - 0.2e1 * t9 * t2, t25, t20, 0, t26, 0, 0, -0.2e1 * t7 * t47 + 0.2e1 * t52, -0.2e1 * t7 * t46 - 0.2e1 * t53, -0.2e1 * t1, 0.2e1 * t8 * t1 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, 0, -t12 * t38 - t13 * t40 + (-t22 * t38 + t23 * t40) * qJD(3), 0, 0, 0, 0, 0, 0, t14, -t15, 0, -t10 * t15 - t9 * t14 + t3 * t19 + t54, 0, 0, 0, 0, 0, 0, t6, -t5, t4, t7 * t14 + t3 * t43 - t8 * t4 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t19 * t15 + 0.2e1 * t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t15 * t43 + 0.2e1 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, 0, (-t2 * t36 + t3 * t35) * pkin(3), t26, -0.2e1 * t24, 0, t25, 0, 0, t37 * t41 - t52, t39 * t41 + t53, t1, t2 * t30 + t3 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t15, 0, (-t14 * t36 - t15 * t35) * pkin(3), 0, 0, 0, 0, 0, 0, -t6, t5, -t4, t14 * t30 - t15 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t20, 0, t26, 0, 0, t37 * t45, t39 * t45, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, t47, 0, -t37 * t3 - t8 * t46, -t39 * t3 + t8 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t15 - t19 * t46, t39 * t15 + t19 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t47, 0, -t29 * t46, t29 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
