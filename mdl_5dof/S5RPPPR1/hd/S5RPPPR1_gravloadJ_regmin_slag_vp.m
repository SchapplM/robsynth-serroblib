% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (110->38), mult. (99->61), div. (0->0), fcn. (105->10), ass. (0->29)
t15 = qJ(1) + pkin(7);
t10 = sin(t15);
t33 = g(1) * t10;
t17 = sin(pkin(8));
t32 = g(3) * t17;
t19 = cos(pkin(8));
t31 = t10 * t19;
t12 = cos(t15);
t30 = t12 * t19;
t16 = sin(pkin(9));
t29 = t16 * t19;
t18 = cos(pkin(9));
t28 = t18 * t19;
t21 = cos(qJ(1));
t27 = t21 * pkin(1) + t12 * pkin(2) + t10 * qJ(3);
t20 = sin(qJ(1));
t26 = -t20 * pkin(1) + t12 * qJ(3);
t25 = -g(1) * t12 - g(2) * t10;
t24 = -g(2) * t12 + t33;
t23 = g(1) * t20 - g(2) * t21;
t22 = pkin(3) * t19 + qJ(4) * t17;
t14 = pkin(9) + qJ(5);
t11 = cos(t14);
t9 = sin(t14);
t4 = t10 * t9 + t11 * t30;
t3 = t10 * t11 - t9 * t30;
t2 = -t11 * t31 + t12 * t9;
t1 = t12 * t11 + t9 * t31;
t5 = [0, t23, g(1) * t21 + g(2) * t20, t23 * pkin(1), t24 * t19, t25, -g(1) * (-t10 * pkin(2) + t26) - g(2) * t27, -g(1) * (-t10 * t28 + t12 * t16) - g(2) * (t10 * t16 + t12 * t28), -g(1) * (t10 * t29 + t12 * t18) - g(2) * (t10 * t18 - t12 * t29), -g(1) * t26 - g(2) * (t22 * t12 + t27) - (-pkin(2) - t22) * t33, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t24, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t19 + t25 * t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t9 * t32, g(1) * t4 - g(2) * t2 + t11 * t32;];
taug_reg = t5;
