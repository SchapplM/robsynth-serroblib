% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:43
% DurationCPUTime: 0.36s
% Computational Cost: add. (171->45), mult. (464->82), div. (0->0), fcn. (387->6), ass. (0->41)
t41 = 2 * qJD(4);
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t28 = sin(qJ(2));
t30 = cos(qJ(2));
t13 = t25 * t28 - t26 * t30;
t14 = t25 * t30 + t26 * t28;
t9 = t14 * qJD(2);
t40 = t13 * t9;
t29 = cos(qJ(4));
t39 = t29 * pkin(4);
t19 = t25 * pkin(2) + pkin(6);
t38 = qJ(5) + t19;
t27 = sin(qJ(4));
t37 = t27 * qJD(4);
t22 = t29 * qJD(4);
t36 = 0.2e1 * t37;
t35 = pkin(4) * t37;
t34 = t13 * t37;
t33 = t27 * t22;
t20 = -t26 * pkin(2) - pkin(3);
t10 = t13 * qJD(2);
t23 = t27 ^ 2;
t24 = t29 ^ 2;
t2 = (t23 + t24) * t10;
t11 = t38 * t27;
t12 = t38 * t29;
t32 = t11 * t27 + t12 * t29;
t4 = t27 * t10 - t14 * t22;
t7 = -t29 * qJD(5) + t38 * t37;
t8 = -t27 * qJD(5) - t38 * t22;
t31 = -t8 * t27 - t7 * t29 + (t11 * t29 - t12 * t27) * qJD(4);
t18 = -0.2e1 * t33;
t17 = 0.2e1 * t33;
t16 = t20 - t39;
t15 = (-t23 + t24) * t41;
t6 = -t9 * t29 + t34;
t5 = t13 * t22 + t9 * t27;
t3 = t29 * t10 + t14 * t37;
t1 = -0.2e1 * t14 * t2 + 0.2e1 * t40;
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t14 * t10 + 0.2e1 * t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * qJD(2), -t30 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t9, t10, 0, (-t10 * t25 - t26 * t9) * pkin(2), 0, 0, 0, 0, 0, 0, t6, t5, -t2, -t19 * t2 + t9 * t20, 0, 0, 0, 0, 0, 0, t6, t5, -t2, pkin(4) * t34 - t32 * t10 + t31 * t14 + t9 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t15, 0, t18, 0, 0, t20 * t36, 0.2e1 * t20 * t22, 0, 0, t17, t15, 0, t18, 0, 0, (t16 - t39) * t36, (pkin(4) * t23 + t16 * t29) * t41, 0.2e1 * t31, -0.2e1 * t11 * t8 - 0.2e1 * t12 * t7 + 0.2e1 * t16 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * qJD(4) - t7 * t27 + t8 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t37, 0, -t19 * t22, t19 * t37, 0, 0, 0, 0, t22, 0, -t37, 0, t8, t7, -pkin(4) * t22, t8 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t22, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t22, 0, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t22, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t21;
