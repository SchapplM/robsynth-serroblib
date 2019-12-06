% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR3
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
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:32
% EndTime: 2019-12-05 15:47:34
% DurationCPUTime: 0.26s
% Computational Cost: add. (186->45), mult. (491->91), div. (0->0), fcn. (466->8), ass. (0->39)
t40 = qJD(4) + qJD(5);
t22 = sin(pkin(9));
t20 = t22 * pkin(2) + pkin(6);
t39 = pkin(7) + t20;
t24 = sin(qJ(5));
t25 = sin(qJ(4));
t38 = t24 * t25;
t27 = cos(qJ(5));
t37 = qJD(5) * t27;
t36 = t25 * qJD(4);
t28 = cos(qJ(4));
t35 = t28 * qJD(4);
t34 = 0.2e1 * t35;
t33 = pkin(4) * t36;
t32 = qJD(5) * t24 * pkin(4);
t31 = pkin(4) * t37;
t23 = cos(pkin(9));
t21 = -t23 * pkin(2) - pkin(3);
t30 = qJD(4) * t39;
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t14 = t22 * t29 + t23 * t26;
t13 = t22 * t26 - t23 * t29;
t16 = t24 * t28 + t27 * t25;
t15 = -t27 * t28 + t38;
t6 = t40 * t16;
t17 = -t28 * pkin(4) + t21;
t12 = t39 * t28;
t11 = t39 * t25;
t10 = t13 * qJD(2);
t9 = t14 * qJD(2);
t8 = t28 * t30;
t7 = t25 * t30;
t5 = -t27 * t35 - t28 * t37 + t40 * t38;
t4 = t24 * t7 - t27 * t8 + (t11 * t24 - t12 * t27) * qJD(5);
t3 = t24 * t8 + t27 * t7 + (t11 * t27 + t12 * t24) * qJD(5);
t2 = t40 * t14 * t15 + t16 * t10;
t1 = -t15 * t10 + t6 * t14;
t18 = [0, 0, 0, 0, -0.2e1 * t14 * t10 + 0.2e1 * t13 * t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t26 * qJD(2), -t29 * qJD(2), (-t10 * t22 - t23 * t9) * pkin(2), 0, 0, 0, 0, 0, t13 * t36 - t9 * t28, t13 * t35 + t9 * t25, 0, 0, 0, 0, 0, t13 * t6 + t9 * t15, -t13 * t5 + t9 * t16; 0, 0, 0, 0, 0, t25 * t34, 0.2e1 * (-t25 ^ 2 + t28 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t21 * t36, t21 * t34, -0.2e1 * t16 * t5, 0.2e1 * t5 * t15 - 0.2e1 * t16 * t6, 0, 0, 0, 0.2e1 * t15 * t33 + 0.2e1 * t17 * t6, 0.2e1 * t16 * t33 - 0.2e1 * t17 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t10 - t14 * t35, t28 * t10 + t14 * t36, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, t35, -t36, 0, -t20 * t35, t20 * t36, 0, 0, -t5, -t6, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t35, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t32, -0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
