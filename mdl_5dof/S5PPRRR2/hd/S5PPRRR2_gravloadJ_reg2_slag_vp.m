% Calculate inertial parameters regressor of gravitation load for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t20 = g(1) * t13 + g(2) * t12;
t10 = pkin(9) + qJ(3);
t6 = sin(t10);
t7 = cos(t10);
t19 = -g(3) * t7 + t20 * t6;
t34 = g(3) * t6;
t11 = qJ(4) + qJ(5);
t8 = sin(t11);
t30 = t12 * t8;
t9 = cos(t11);
t29 = t12 * t9;
t28 = t13 * t8;
t27 = t13 * t9;
t14 = sin(qJ(4));
t26 = t12 * t14;
t15 = cos(qJ(4));
t25 = t12 * t15;
t24 = t13 * t14;
t23 = t13 * t15;
t17 = -g(1) * (-t7 * t24 + t25) - g(2) * (-t7 * t26 - t23) + t14 * t34;
t16 = -pkin(7) - pkin(6);
t5 = t15 * pkin(4) + pkin(3);
t4 = -g(1) * t12 + g(2) * t13;
t3 = t20 * t7 + t34;
t2 = -g(1) * (-t7 * t27 - t30) - g(2) * (-t7 * t29 + t28) + t9 * t34;
t1 = -g(1) * (-t7 * t28 + t29) - g(2) * (-t7 * t30 - t27) + t8 * t34;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t3, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t15, -t19 * t14, -t3, -g(3) * (t7 * pkin(3) + t6 * pkin(6)) + t20 * (pkin(3) * t6 - pkin(6) * t7), 0, 0, 0, 0, 0, 0, t19 * t9, -t19 * t8, -t3, -g(3) * (-t6 * t16 + t7 * t5) + t20 * (t16 * t7 + t5 * t6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -g(1) * (-t7 * t23 - t26) - g(2) * (-t7 * t25 + t24) + t15 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t18;
