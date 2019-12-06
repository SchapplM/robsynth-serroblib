% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:41
% EndTime: 2019-12-05 15:54:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (328->56), mult. (884->113), div. (0->0), fcn. (870->8), ass. (0->49)
t30 = sin(pkin(9));
t31 = cos(pkin(9));
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t21 = t33 * t30 - t36 * t31;
t19 = t21 * qJD(4);
t50 = pkin(6) + qJ(3);
t24 = t50 * t30;
t25 = t50 * t31;
t40 = -t36 * t24 - t33 * t25;
t55 = qJD(3) * t21 - qJD(4) * t40;
t27 = -t31 * pkin(3) - pkin(2);
t54 = 0.2e1 * t27;
t22 = t36 * t30 + t33 * t31;
t20 = t22 * qJD(4);
t53 = pkin(4) * t20;
t49 = t30 ^ 2 + t31 ^ 2;
t48 = pkin(4) * qJD(5);
t34 = sin(qJ(2));
t47 = t34 * qJD(2);
t37 = cos(qJ(2));
t46 = t37 * qJD(2);
t32 = sin(qJ(5));
t45 = t32 * t48;
t35 = cos(qJ(5));
t44 = t35 * t48;
t43 = t49 * t37;
t42 = t49 * qJD(3);
t41 = 0.2e1 * t42;
t13 = t35 * t21 + t32 * t22;
t14 = -t32 * t21 + t35 * t22;
t39 = t33 * t24 - t36 * t25;
t38 = -qJD(3) * t22 + qJD(4) * t39;
t18 = t21 * t34;
t17 = t22 * t34;
t15 = t21 * pkin(4) + t27;
t12 = t34 * t19 - t22 * t46;
t11 = t20 * t34 + t21 * t46;
t10 = -t21 * pkin(7) - t39;
t9 = -t22 * pkin(7) + t40;
t8 = t19 * pkin(7) + t38;
t7 = -t20 * pkin(7) - t55;
t6 = qJD(5) * t14 - t32 * t19 + t35 * t20;
t5 = -qJD(5) * t13 - t35 * t19 - t32 * t20;
t4 = t32 * t11 + t35 * t12 + (t17 * t32 + t18 * t35) * qJD(5);
t3 = t35 * t11 - t32 * t12 + (t17 * t35 - t18 * t32) * qJD(5);
t2 = -t32 * t7 + t35 * t8 + (-t10 * t35 - t32 * t9) * qJD(5);
t1 = -t32 * t8 - t35 * t7 + (t10 * t32 - t35 * t9) * qJD(5);
t16 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t49) * t34 * t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t47, -t46, -t31 * t47, t30 * t47, qJD(2) * t43, t34 * t42 + (-pkin(2) * t34 + qJ(3) * t43) * qJD(2), 0, 0, 0, 0, 0, -t37 * t20 + t21 * t47, t37 * t19 + t22 * t47, 0, 0, 0, 0, 0, t13 * t47 - t37 * t6, t14 * t47 - t37 * t5; 0, 0, 0, 0, 0, 0, t41, qJ(3) * t41, -0.2e1 * t22 * t19, 0.2e1 * t19 * t21 - 0.2e1 * t22 * t20, 0, 0, 0, t20 * t54, -t19 * t54, 0.2e1 * t14 * t5, -0.2e1 * t5 * t13 - 0.2e1 * t14 * t6, 0, 0, 0, 0.2e1 * t13 * t53 + 0.2e1 * t15 * t6, 0.2e1 * t14 * t53 + 0.2e1 * t15 * t5; 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, t38, t55, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t45, -0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
