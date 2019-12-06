% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:09
% EndTime: 2019-12-05 15:22:11
% DurationCPUTime: 0.23s
% Computational Cost: add. (115->39), mult. (375->94), div. (0->0), fcn. (334->6), ass. (0->37)
t22 = sin(pkin(9));
t24 = cos(pkin(9));
t26 = sin(qJ(5));
t27 = cos(qJ(5));
t43 = -t22 * t26 + t24 * t27;
t11 = t43 * qJD(5);
t25 = cos(pkin(8));
t42 = 0.2e1 * t25;
t23 = sin(pkin(8));
t41 = pkin(6) * t23;
t16 = -t25 * pkin(3) - t23 * qJ(4) - pkin(2);
t37 = qJ(3) * t25;
t38 = t22 * t16 + t24 * t37;
t36 = qJD(3) * t25;
t35 = qJD(4) * t23;
t20 = t23 ^ 2;
t33 = t20 * qJD(3);
t32 = t23 * qJD(3);
t31 = qJ(3) * qJD(3);
t10 = -t22 * t35 + t24 * t36;
t9 = -t22 * t36 - t24 * t35;
t30 = t10 * t22 + t9 * t24;
t29 = t22 * t27 + t24 * t26;
t7 = t29 * t23;
t21 = t25 ^ 2;
t19 = t20 * t31;
t15 = (pkin(4) * t22 + qJ(3)) * t23;
t14 = t24 * t16;
t12 = t29 * qJD(5);
t8 = t43 * t23;
t6 = t23 * t11;
t5 = qJD(5) * t7;
t4 = -t22 * t41 + t38;
t3 = -t24 * t41 + t14 + (-qJ(3) * t22 - pkin(4)) * t25;
t2 = -t26 * t10 + t27 * t9 + (-t26 * t3 - t27 * t4) * qJD(5);
t1 = -t27 * t10 - t26 * t9 + (t26 * t4 - t27 * t3) * qJD(5);
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t10 * t24 - t22 * t9 - t36) * t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0.2e1 * (t20 + t21) * qJD(3), 0.2e1 * t21 * t31 + 0.2e1 * t19, 0.2e1 * t22 * t33 - 0.2e1 * t9 * t25, 0.2e1 * t10 * t25 + 0.2e1 * t24 * t33, -0.2e1 * t30 * t23, 0.2e1 * t38 * t10 + 0.2e1 * (-t22 * t37 + t14) * t9 + 0.2e1 * t19, -0.2e1 * t8 * t5, 0.2e1 * t5 * t7 - 0.2e1 * t8 * t6, t5 * t42, t6 * t42, 0, 0.2e1 * t15 * t6 - 0.2e1 * t2 * t25 + 0.2e1 * t7 * t32, -0.2e1 * t1 * t25 - 0.2e1 * t15 * t5 + 0.2e1 * t8 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, t12 * t25, t11 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
