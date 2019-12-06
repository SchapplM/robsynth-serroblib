% Calculate minimal parameter regressor of gravitation load for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% taug_reg [5x13]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = sin(pkin(9));
t13 = sin(pkin(8));
t28 = t12 * t13;
t19 = sin(qJ(4));
t27 = t13 * t19;
t21 = cos(qJ(4));
t26 = t13 * t21;
t14 = sin(pkin(7));
t16 = cos(pkin(8));
t25 = t14 * t16;
t17 = cos(pkin(7));
t24 = t17 * t12;
t15 = cos(pkin(9));
t23 = t17 * t15;
t6 = t15 * t25 - t24;
t8 = t14 * t12 + t16 * t23;
t22 = g(1) * (t17 * t26 - t8 * t19) + g(2) * (t14 * t26 - t6 * t19) + g(3) * (-t15 * t27 - t16 * t21);
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t11 = -g(1) * t14 + g(2) * t17;
t10 = t15 * t26 - t16 * t19;
t7 = -t14 * t15 + t16 * t24;
t5 = t12 * t25 + t23;
t4 = t17 * t27 + t8 * t21;
t2 = t14 * t27 + t6 * t21;
t1 = [-g(3), -g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t11, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(3) * t16 + (-g(1) * t17 - g(2) * t14) * t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t22, g(1) * t4 + g(2) * t2 + g(3) * t10, 0, 0, 0, 0, 0, -t22 * t20, t22 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t18 + t7 * t20) - g(2) * (-t2 * t18 + t5 * t20) - g(3) * (-t10 * t18 + t20 * t28), -g(1) * (-t7 * t18 - t4 * t20) - g(2) * (-t5 * t18 - t2 * t20) - g(3) * (-t10 * t20 - t18 * t28);];
taug_reg = t1;
