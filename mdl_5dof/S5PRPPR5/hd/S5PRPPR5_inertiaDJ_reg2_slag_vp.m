% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:38
% EndTime: 2019-12-31 17:38:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (114->42), mult. (308->86), div. (0->0), fcn. (263->6), ass. (0->33)
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t3 = t17 * t14 + t19 * t15;
t2 = t3 * qJD(2);
t16 = sin(qJ(5));
t12 = t16 ^ 2;
t18 = cos(qJ(5));
t13 = t18 ^ 2;
t31 = t12 + t13;
t22 = t31 * t2;
t35 = 2 * qJD(3);
t25 = t19 * qJD(2);
t27 = t17 * qJD(2);
t1 = t14 * t25 - t15 * t27;
t34 = t3 * t1;
t33 = t1 * t15;
t32 = t14 * t3;
t20 = -pkin(2) - pkin(3);
t8 = t15 * qJ(3) + t14 * t20;
t30 = t14 * qJD(3);
t29 = t15 * qJD(3);
t28 = t16 * qJD(5);
t26 = t18 * qJD(5);
t24 = 0.2e1 * t29;
t23 = t16 * t26;
t21 = t15 * t31;
t7 = -t14 * qJ(3) + t15 * t20;
t6 = -pkin(6) + t8;
t5 = pkin(4) - t7;
t4 = -t19 * t14 + t17 * t15;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t4 * t2 + 0.2e1 * t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t4 * t22 + 0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, t25, t17 * qJD(3) + (-pkin(2) * t17 + qJ(3) * t19) * qJD(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1 * t7 + t2 * t8 + (t15 * t4 + t32) * qJD(3), 0, 0, 0, 0, 0, 0, t1 * t18 - t3 * t28, -t1 * t16 - t3 * t26, -t22, t1 * t5 + t6 * t22 + (t4 * t21 + t32) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, qJ(3) * t35, 0, 0, 0, 0, 0, 0, 0.2e1 * t30, t24, 0, (-t14 * t7 + t15 * t8) * t35, 0.2e1 * t23, 0.2e1 * (-t12 + t13) * qJD(5), 0, -0.2e1 * t23, 0, 0, 0.2e1 * t18 * t30 - 0.2e1 * t5 * t28, -0.2e1 * t16 * t30 - 0.2e1 * t5 * t26, -t31 * t24, (t14 * t5 + t6 * t21) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t14 - t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t22 - t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t28, t15 * t26, 0, (-0.1e1 + t31) * t14 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t2 - t4 * t26, -t18 * t2 + t4 * t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, t28, 0, -t16 * t29 - t6 * t26, -t18 * t29 + t6 * t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t26, t14 * t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
