% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:06
% EndTime: 2019-12-05 15:36:08
% DurationCPUTime: 0.38s
% Computational Cost: add. (134->41), mult. (396->72), div. (0->0), fcn. (336->6), ass. (0->37)
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t28 = sin(qJ(2));
t30 = cos(qJ(2));
t13 = t25 * t28 - t26 * t30;
t11 = t13 * qJD(2);
t27 = sin(qJ(4));
t23 = t27 ^ 2;
t29 = cos(qJ(4));
t24 = t29 ^ 2;
t2 = (t23 + t24) * t11;
t40 = 2 * qJD(5);
t18 = t25 * pkin(2) + pkin(6);
t39 = t18 * t2;
t14 = t25 * t30 + t26 * t28;
t10 = t14 * qJD(2);
t37 = t13 * t10;
t21 = t27 * qJD(4);
t22 = t29 * qJD(4);
t19 = -t26 * pkin(2) - pkin(3);
t36 = 0.2e1 * qJD(4) * t19;
t35 = t27 * t22;
t34 = t18 * t21;
t33 = t18 * t22;
t32 = -t29 * pkin(4) - t27 * qJ(5);
t31 = t32 * qJD(4) + t29 * qJD(5);
t17 = -0.2e1 * t35;
t16 = 0.2e1 * t35;
t15 = (-t23 + t24) * qJD(4);
t12 = t19 + t32;
t9 = -pkin(4) * t21 + qJ(5) * t22 + t27 * qJD(5);
t6 = -t10 * t29 + t13 * t21;
t5 = t10 * t27 + t13 * t22;
t4 = -t27 * t11 + t14 * t22;
t3 = t29 * t11 + t14 * t21;
t1 = -0.2e1 * t14 * t2 + 0.2e1 * t37;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t14 * t11 + 0.2e1 * t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * qJD(2), -t30 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t10, t11, 0, (-t10 * t26 - t11 * t25) * pkin(2), 0, 0, 0, 0, 0, 0, t6, t5, -t2, t10 * t19 - t39, 0, 0, 0, 0, 0, 0, t6, -t2, -t5, t10 * t12 - t13 * t9 - t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0.2e1 * t15, 0, t17, 0, 0, t27 * t36, t29 * t36, 0, 0, t16, 0, -0.2e1 * t15, 0, 0, t17, 0.2e1 * t12 * t21 + 0.2e1 * t9 * t29, 0, -0.2e1 * t12 * t22 + 0.2e1 * t9 * t27, -0.2e1 * t12 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, -t3, -(-pkin(4) * t27 + qJ(5) * t29) * t11 + t31 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t21, 0, -t33, t34, 0, 0, 0, t22, 0, 0, t21, 0, -t33, t31, -t34, t31 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, t22, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, qJ(5) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
