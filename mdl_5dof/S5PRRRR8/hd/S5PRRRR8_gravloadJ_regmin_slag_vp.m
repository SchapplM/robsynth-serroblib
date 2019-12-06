% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(pkin(5));
t42 = g(3) * t20;
t18 = qJ(3) + qJ(4);
t17 = cos(t18);
t23 = sin(qJ(5));
t41 = t17 * t23;
t26 = cos(qJ(5));
t40 = t17 * t26;
t19 = sin(pkin(10));
t39 = t19 * t20;
t21 = cos(pkin(10));
t38 = t20 * t21;
t24 = sin(qJ(3));
t37 = t20 * t24;
t25 = sin(qJ(2));
t36 = t20 * t25;
t27 = cos(qJ(3));
t35 = t20 * t27;
t22 = cos(pkin(5));
t34 = t22 * t25;
t28 = cos(qJ(2));
t33 = t22 * t28;
t32 = t23 * t28;
t31 = t26 * t28;
t12 = t19 * t28 + t21 * t34;
t14 = -t19 * t34 + t21 * t28;
t16 = sin(t18);
t30 = g(1) * (-t14 * t16 + t17 * t39) + g(2) * (-t12 * t16 - t17 * t38) + g(3) * (-t16 * t36 + t22 * t17);
t11 = t19 * t25 - t21 * t33;
t13 = t19 * t33 + t21 * t25;
t29 = -g(1) * t13 - g(2) * t11 + t28 * t42;
t10 = t22 * t16 + t17 * t36;
t8 = t14 * t17 + t16 * t39;
t6 = t12 * t17 - t16 * t38;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t10;
t2 = t30 * t26;
t1 = t30 * t23;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t29, g(1) * t14 + g(2) * t12 + g(3) * t36, 0, 0, 0, 0, 0, -t29 * t27, t29 * t24, 0, 0, 0, 0, 0, -t29 * t17, t29 * t16, 0, 0, 0, 0, 0, -g(1) * (-t13 * t40 + t14 * t23) - g(2) * (-t11 * t40 + t12 * t23) - (t17 * t31 + t23 * t25) * t42, -g(1) * (t13 * t41 + t14 * t26) - g(2) * (t11 * t41 + t12 * t26) - (-t17 * t32 + t25 * t26) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t24 + t19 * t35) - g(2) * (-t12 * t24 - t21 * t35) - g(3) * (t22 * t27 - t24 * t36), -g(1) * (-t14 * t27 - t19 * t37) - g(2) * (-t12 * t27 + t21 * t37) - g(3) * (-t22 * t24 - t25 * t35), 0, 0, 0, 0, 0, -t30, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t26 - t8 * t23) - g(2) * (t11 * t26 - t6 * t23) - g(3) * (-t10 * t23 - t20 * t31), -g(1) * (-t13 * t23 - t8 * t26) - g(2) * (-t11 * t23 - t6 * t26) - g(3) * (-t10 * t26 + t20 * t32);];
taug_reg = t3;
