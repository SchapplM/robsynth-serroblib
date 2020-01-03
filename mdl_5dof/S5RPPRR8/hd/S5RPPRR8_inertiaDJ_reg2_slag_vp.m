% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR8
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
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:14
% EndTime: 2019-12-31 18:01:15
% DurationCPUTime: 0.40s
% Computational Cost: add. (416->51), mult. (759->98), div. (0->0), fcn. (620->6), ass. (0->46)
t28 = sin(pkin(8));
t30 = sin(qJ(4));
t43 = cos(pkin(8));
t47 = cos(qJ(4));
t35 = t47 * t43;
t16 = t30 * t28 - t35;
t52 = 2 * qJD(2);
t51 = -pkin(1) - pkin(2);
t19 = t43 * qJ(2) + t28 * t51;
t34 = -t28 * qJ(2) + t43 * t51;
t33 = -pkin(3) + t34;
t10 = t47 * t19 + t30 * t33;
t17 = t47 * t28 + t30 * t43;
t3 = t17 * qJD(2) + t10 * qJD(4);
t50 = t3 * t16;
t29 = sin(qJ(5));
t49 = t3 * t29;
t31 = cos(qJ(5));
t48 = t3 * t31;
t13 = t17 * qJD(4);
t46 = t16 * t13;
t26 = t29 ^ 2;
t27 = t31 ^ 2;
t44 = -t26 - t27;
t42 = t28 * qJD(2);
t41 = t29 * qJD(5);
t40 = t31 * qJD(5);
t39 = -0.2e1 * pkin(4) * qJD(5);
t38 = t29 * t40;
t32 = t47 * t33;
t2 = -qJD(4) * t32 - qJD(2) * t35 + (qJD(4) * t19 + t42) * t30;
t1 = t44 * t2;
t9 = -t30 * t19 + t32;
t7 = pkin(4) - t9;
t37 = qJD(5) * (pkin(4) + t7);
t12 = t16 * qJD(4);
t4 = t44 * t12;
t36 = t44 * t17;
t22 = -0.2e1 * t38;
t21 = 0.2e1 * t38;
t20 = (-t26 + t27) * qJD(5);
t18 = 0.2e1 * t20;
t8 = -pkin(7) + t10;
t6 = t13 * t31 - t16 * t41;
t5 = t13 * t29 + t16 * t40;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, qJ(2) * t52, 0, 0, 0, 0, 0, 0, 0.2e1 * t42, t43 * t52, 0, (t19 * t43 - t34 * t28) * t52, 0, 0, 0, 0, 0, 0, 0.2e1 * t3, -0.2e1 * t2, 0, -0.2e1 * t10 * t2 - 0.2e1 * t9 * t3, t21, t18, 0, t22, 0, 0, -0.2e1 * t7 * t41 + 0.2e1 * t48, -0.2e1 * t7 * t40 - 0.2e1 * t49, -0.2e1 * t1, 0.2e1 * t8 * t1 + 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t10 * t12 - t9 * t13 - t2 * t17 + t50, 0, 0, 0, 0, 0, 0, t6, -t5, -t4, t7 * t13 + t2 * t36 + t8 * t4 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t12 + 0.2e1 * t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 * t36 + 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, 0, 0, t22, -0.2e1 * t20, 0, t21, 0, 0, t29 * t37 - t48, t31 * t37 + t49, t1, -t3 * pkin(4) + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, t4, -t13 * pkin(4) + pkin(7) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t18, 0, t22, 0, 0, t29 * t39, t31 * t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, t41, 0, t29 * t2 - t8 * t40, t31 * t2 + t8 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t12 - t17 * t40, t31 * t12 + t17 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t41, 0, -pkin(7) * t40, pkin(7) * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
