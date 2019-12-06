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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
t12 = pkin(9) + qJ(5);
t8 = sin(t12);
t13 = qJ(1) + pkin(7);
t9 = sin(t13);
t32 = t9 * t8;
t15 = sin(pkin(8));
t31 = g(1) * t15;
t11 = cos(t13);
t30 = g(2) * t11;
t19 = cos(qJ(1));
t29 = t19 * pkin(1);
t10 = cos(t12);
t28 = t9 * t10;
t17 = cos(pkin(8));
t27 = t11 * t17;
t14 = sin(pkin(9));
t26 = t14 * t17;
t16 = cos(pkin(9));
t25 = t16 * t17;
t18 = sin(qJ(1));
t24 = -t18 * pkin(1) + t11 * qJ(3);
t23 = g(2) * t9 - g(3) * t11;
t22 = g(3) * t9 + t30;
t21 = g(2) * t19 + g(3) * t18;
t20 = pkin(3) * t17 + qJ(4) * t15 + pkin(2);
t5 = t22 * t15;
t4 = -t10 * t27 - t32;
t3 = t8 * t27 - t28;
t2 = -t11 * t8 + t17 * t28;
t1 = t11 * t10 + t17 * t32;
t6 = [0, t21, -g(2) * t18 + g(3) * t19, t21 * pkin(1), t22 * t17, -t5, t23, -g(2) * (-t11 * pkin(2) - t9 * qJ(3) - t29) - g(3) * (-t9 * pkin(2) + t24), -g(2) * (-t11 * t25 - t9 * t14) - g(3) * (t11 * t14 - t9 * t25), -g(2) * (t11 * t26 - t9 * t16) - g(3) * (t11 * t16 + t9 * t26), t5, g(2) * t29 - g(3) * t24 + (g(2) * qJ(3) + g(3) * t20) * t9 + t20 * t30, 0, 0, 0, 0, 0, -g(2) * t4 + g(3) * t2, -g(2) * t3 - g(3) * t1; 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t17 + t23 * t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t1 + g(3) * t3 + t8 * t31, -g(2) * t2 - g(3) * t4 + t10 * t31;];
taug_reg = t6;
