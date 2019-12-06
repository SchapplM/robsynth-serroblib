% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = qJ(1) + qJ(2);
t18 = sin(t23);
t42 = pkin(2) * t18;
t20 = cos(t23);
t41 = pkin(2) * t20;
t25 = sin(qJ(1));
t40 = t25 * pkin(1);
t27 = cos(qJ(1));
t39 = t27 * pkin(1);
t21 = qJ(3) + t23;
t14 = sin(t21);
t15 = cos(t21);
t38 = -t14 * pkin(3) + t15 * pkin(8);
t26 = cos(qJ(4));
t16 = t26 * pkin(4) + pkin(3);
t28 = -pkin(9) - pkin(8);
t37 = t14 * t28 - t15 * t16;
t36 = -t15 * pkin(3) - t14 * pkin(8);
t8 = g(2) * t15 + g(3) * t14;
t7 = g(2) * t14 - g(3) * t15;
t10 = g(2) * t20 + g(3) * t18;
t35 = g(2) * t27 + g(3) * t25;
t34 = -t14 * t16 - t15 * t28;
t33 = t38 - t42;
t32 = t37 - t41;
t31 = t36 - t41;
t30 = t34 - t42;
t24 = sin(qJ(4));
t29 = -g(1) * t26 - t7 * t24;
t22 = qJ(4) + qJ(5);
t19 = cos(t22);
t17 = sin(t22);
t9 = -g(2) * t18 + g(3) * t20;
t6 = t8 * t26;
t5 = t8 * t24;
t4 = t8 * t19;
t3 = t8 * t17;
t2 = -g(1) * t19 - t7 * t17;
t1 = g(1) * t17 - t7 * t19;
t11 = [0, 0, 0, 0, 0, 0, t35, -g(2) * t25 + g(3) * t27, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, t35 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, 0, -g(2) * (-t39 - t41) - g(3) * (-t40 - t42), 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * (t31 - t39) - g(3) * (t33 - t40), 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * (t32 - t39) - g(3) * (t30 - t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t10 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * t31 - g(3) * t33, 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t32 - g(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * t36 - g(3) * t38, 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t37 - g(3) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(1) * t24 - t7 * t26, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t11;
