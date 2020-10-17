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
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:33
% EndTime: 2020-01-03 11:20:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (116->37), mult. (107->62), div. (0->0), fcn. (113->10), ass. (0->29)
t19 = sin(pkin(8));
t33 = g(1) * t19;
t17 = qJ(1) + pkin(7);
t11 = sin(t17);
t21 = cos(pkin(8));
t32 = t11 * t21;
t13 = cos(t17);
t31 = t13 * t21;
t18 = sin(pkin(9));
t30 = t18 * t21;
t20 = cos(pkin(9));
t29 = t20 * t21;
t23 = cos(qJ(1));
t28 = t23 * pkin(1) + t13 * pkin(2) + t11 * qJ(3);
t22 = sin(qJ(1));
t27 = t22 * pkin(1) + t11 * pkin(2) - t13 * qJ(3);
t6 = g(2) * t13 + g(3) * t11;
t26 = -g(2) * t11 + g(3) * t13;
t25 = -g(2) * t23 - g(3) * t22;
t24 = pkin(3) * t21 + qJ(4) * t19;
t16 = pkin(9) + qJ(5);
t12 = cos(t16);
t10 = sin(t16);
t5 = t6 * t19;
t4 = t11 * t10 + t12 * t31;
t3 = t10 * t31 - t11 * t12;
t2 = -t13 * t10 + t12 * t32;
t1 = -t10 * t32 - t13 * t12;
t7 = [0, t25, g(2) * t22 - g(3) * t23, t25 * pkin(1), -t6 * t21, t5, t26, -g(2) * t28 - g(3) * t27, -g(2) * (t11 * t18 + t13 * t29) - g(3) * (t11 * t29 - t13 * t18), -g(2) * (t11 * t20 - t13 * t30) - g(3) * (-t11 * t30 - t13 * t20), -t5, -g(2) * (t24 * t13 + t28) - g(3) * (t24 * t11 + t27), 0, 0, 0, 0, 0, -g(2) * t4 - g(3) * t2, g(2) * t3 - g(3) * t1; 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t21 + t26 * t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t1 - g(3) * t3 + t10 * t33, g(2) * t2 - g(3) * t4 + t12 * t33;];
taug_reg = t7;
