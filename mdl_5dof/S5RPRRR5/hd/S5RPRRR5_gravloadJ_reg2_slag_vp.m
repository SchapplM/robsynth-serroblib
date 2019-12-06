% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = qJ(1) + pkin(9);
t17 = qJ(3) + t20;
t12 = sin(t17);
t13 = cos(t17);
t34 = -t12 * pkin(3) + t13 * pkin(7);
t24 = cos(qJ(4));
t14 = t24 * pkin(4) + pkin(3);
t26 = -pkin(8) - pkin(7);
t33 = t12 * t26 - t13 * t14;
t15 = sin(t20);
t23 = sin(qJ(1));
t32 = -t23 * pkin(1) - pkin(2) * t15;
t16 = cos(t20);
t25 = cos(qJ(1));
t31 = -t25 * pkin(1) - pkin(2) * t16;
t30 = -t13 * pkin(3) - t12 * pkin(7);
t8 = g(2) * t13 + g(3) * t12;
t7 = g(2) * t12 - g(3) * t13;
t29 = g(2) * t25 + g(3) * t23;
t28 = -t12 * t14 - t13 * t26;
t22 = sin(qJ(4));
t27 = -g(1) * t24 - t7 * t22;
t21 = qJ(4) + qJ(5);
t19 = cos(t21);
t18 = sin(t21);
t6 = t8 * t24;
t5 = t8 * t22;
t4 = t8 * t19;
t3 = t8 * t18;
t2 = -g(1) * t19 - t7 * t18;
t1 = g(1) * t18 - t7 * t19;
t9 = [0, 0, 0, 0, 0, 0, t29, -g(2) * t23 + g(3) * t25, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t16 + g(3) * t15, -g(2) * t15 + g(3) * t16, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, 0, -g(2) * t31 - g(3) * t32, 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * (t30 + t31) - g(3) * (t32 + t34), 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * (t31 + t33) - g(3) * (t28 + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * t30 - g(3) * t34, 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t33 - g(3) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(1) * t22 - t7 * t24, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t9;
