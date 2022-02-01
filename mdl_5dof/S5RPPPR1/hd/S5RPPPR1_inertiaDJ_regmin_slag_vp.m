% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:46
% DurationCPUTime: 0.24s
% Computational Cost: add. (137->41), mult. (386->94), div. (0->0), fcn. (347->8), ass. (0->36)
t23 = sin(pkin(9));
t25 = cos(pkin(9));
t27 = sin(qJ(5));
t28 = cos(qJ(5));
t42 = -t23 * t27 + t25 * t28;
t14 = t42 * qJD(5);
t26 = cos(pkin(8));
t41 = 0.2e1 * t26;
t24 = sin(pkin(8));
t40 = pkin(6) * t24;
t16 = -cos(pkin(7)) * pkin(1) - pkin(2) - t24 * qJ(4) - t26 * pkin(3);
t20 = sin(pkin(7)) * pkin(1) + qJ(3);
t38 = t20 * t26;
t39 = t23 * t16 + t25 * t38;
t35 = qJD(3) * t26;
t34 = qJD(4) * t24;
t21 = t24 ^ 2;
t32 = t21 * qJD(3);
t31 = t24 * qJD(3);
t30 = t23 * t28 + t25 * t27;
t7 = t30 * t24;
t22 = t26 ^ 2;
t18 = t20 * t32;
t15 = t30 * qJD(5);
t13 = -t23 * t34 + t25 * t35;
t12 = -t23 * t35 - t25 * t34;
t11 = (pkin(4) * t23 + t20) * t24;
t10 = t25 * t16;
t8 = t42 * t24;
t6 = t24 * t14;
t5 = qJD(5) * t7;
t4 = -t23 * t40 + t39;
t3 = -t25 * t40 + t10 + (-t20 * t23 - pkin(4)) * t26;
t2 = t28 * t12 - t27 * t13 + (-t27 * t3 - t28 * t4) * qJD(5);
t1 = -t27 * t12 - t28 * t13 + (t27 * t4 - t28 * t3) * qJD(5);
t9 = [0, 0, 0, 0, 0, 0.2e1 * (t21 + t22) * qJD(3), 0.2e1 * t22 * t20 * qJD(3) + 0.2e1 * t18, -0.2e1 * t12 * t26 + 0.2e1 * t23 * t32, 0.2e1 * t13 * t26 + 0.2e1 * t25 * t32, 0.2e1 * t39 * t13 + 0.2e1 * (-t23 * t38 + t10) * t12 + 0.2e1 * t18, -0.2e1 * t8 * t5, 0.2e1 * t5 * t7 - 0.2e1 * t8 * t6, t5 * t41, t6 * t41, 0, 0.2e1 * t11 * t6 - 0.2e1 * t2 * t26 + 0.2e1 * t7 * t31, -0.2e1 * t1 * t26 - 0.2e1 * t11 * t5 + 0.2e1 * t8 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t12 * t23 + t13 * t25 - t35) * t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t25 + t13 * t23, 0, 0, 0, 0, 0, t15 * t26, t14 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
