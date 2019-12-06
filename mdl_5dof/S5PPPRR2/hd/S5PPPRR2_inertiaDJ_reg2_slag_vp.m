% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:38
% EndTime: 2019-12-05 14:59:39
% DurationCPUTime: 0.43s
% Computational Cost: add. (196->60), mult. (689->121), div. (0->0), fcn. (703->8), ass. (0->43)
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t40 = cos(pkin(8));
t17 = sin(pkin(8));
t18 = cos(pkin(9));
t42 = t17 * t18;
t10 = -t40 * t20 + t22 * t42;
t7 = qJD(4) * t10;
t9 = t20 * t42 + t22 * t40;
t46 = t9 * t7;
t45 = t20 * t7;
t16 = sin(pkin(9));
t44 = t16 * t17;
t43 = t16 * t22;
t19 = sin(qJ(5));
t14 = t19 ^ 2;
t21 = cos(qJ(5));
t15 = t21 ^ 2;
t41 = t14 + t15;
t8 = qJD(4) * t9;
t39 = t19 * qJD(5);
t38 = t20 * qJD(4);
t37 = t21 * qJD(5);
t36 = t22 * qJD(4);
t35 = -0.2e1 * pkin(4) * qJD(5);
t33 = t19 * t37;
t32 = t16 * t38;
t31 = t20 * t36;
t30 = t16 * t36;
t29 = t41 * t22;
t3 = -t10 * t19 + t21 * t44;
t4 = t10 * t21 + t19 * t44;
t12 = -t18 * t19 + t21 * t43;
t27 = t18 * t21 + t19 * t43;
t26 = t20 * t39 - t21 * t36;
t25 = t19 * t36 + t20 * t37;
t1 = -qJD(5) * t4 + t8 * t19;
t2 = qJD(5) * t3 - t8 * t21;
t24 = -t1 * t19 + t2 * t21 + (-t19 * t4 - t21 * t3) * qJD(5);
t5 = qJD(5) * t27 + t21 * t32;
t6 = -qJD(5) * t12 + t19 * t32;
t23 = -t6 * t19 - t5 * t21 + (-t12 * t19 + t21 * t27) * qJD(5);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10 * t8 + 0.2e1 * t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t1 * t3 + 0.2e1 * t2 * t4 + 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t45 - t22 * t8 + (-t10 * t20 + t22 * t9) * qJD(4)) * t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t27 + t2 * t12 + t3 * t6 - t4 * t5 + (t36 * t9 + t45) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16 ^ 2 * t31 - 0.2e1 * t12 * t5 - 0.2e1 * t27 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t20 - t7 * t22 + (t10 * t22 + t20 * t9) * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t7 + (-t19 * t3 + t21 * t4) * qJD(4)) * t22 + (t24 + t8) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t12 * t21 + t19 * t27 - t43) * t36 + (t23 + t32) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t41) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t7 + t39 * t9, t19 * t7 + t37 * t9, t24, -t7 * pkin(4) + pkin(6) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t32, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t16, t25 * t16, t23, -pkin(4) * t30 + pkin(6) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t36, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t38 - t22 * t39, t19 * t38 - t22 * t37, qJD(4) * t29, (-pkin(4) * t20 + pkin(6) * t29) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t33, 0.2e1 * (-t14 + t15) * qJD(5), 0, -0.2e1 * t33, 0, 0, t19 * t35, t21 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t39, 0, -pkin(6) * t37, pkin(6) * t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
