% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:50
% EndTime: 2019-12-05 16:15:52
% DurationCPUTime: 0.24s
% Computational Cost: add. (117->39), mult. (352->68), div. (0->0), fcn. (268->6), ass. (0->36)
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t40 = t27 ^ 2 + t28 ^ 2;
t51 = t40 * qJD(4);
t52 = 0.2e1 * t51;
t32 = cos(qJ(3));
t39 = pkin(2) * qJD(3);
t37 = t32 * t39;
t18 = qJD(4) + t37;
t50 = t40 * t18;
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t11 = t29 * t27 - t31 * t28;
t21 = -t28 * pkin(4) - pkin(3);
t45 = t32 * pkin(2);
t15 = t21 - t45;
t30 = sin(qJ(3));
t38 = t30 * t39;
t12 = t31 * t27 + t29 * t28;
t8 = t12 * qJD(5);
t49 = t11 * t38 + t15 * t8;
t7 = t11 * qJD(5);
t48 = t12 * t38 - t15 * t7;
t47 = t21 * t7;
t46 = t21 * t8;
t34 = t27 * t38;
t33 = t28 * t38;
t24 = t28 * pkin(7);
t20 = t30 * pkin(2) + qJ(4);
t17 = t28 * qJ(4) + t24;
t16 = (-pkin(7) - qJ(4)) * t27;
t10 = t28 * t20 + t24;
t9 = (-pkin(7) - t20) * t27;
t2 = -0.2e1 * t12 * t7;
t1 = 0.2e1 * t7 * t11 - 0.2e1 * t12 * t8;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t38, -0.2e1 * t37, -0.2e1 * t33, 0.2e1 * t34, 0.2e1 * t50, 0.2e1 * (-pkin(3) - t45) * t38 + 0.2e1 * t20 * t50, t2, t1, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t38, -t37, -t33, t34, t51 + t50, -pkin(3) * t38 + qJ(4) * t50 + t20 * t51, t2, t1, 0, 0, 0, t46 + t49, -t47 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, qJ(4) * t52, t2, t1, 0, 0, 0, 0.2e1 * t46, -0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, t8, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t12 * t18 + (-t10 * t31 - t29 * t9) * qJD(5), t11 * t18 + (t10 * t29 - t31 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, (-t16 * t29 - t17 * t31) * qJD(5) - t12 * qJD(4), (-t16 * t31 + t17 * t29) * qJD(5) + t11 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
