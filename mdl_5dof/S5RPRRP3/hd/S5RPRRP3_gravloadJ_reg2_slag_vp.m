% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = -pkin(7) - pkin(6);
t20 = sin(qJ(1));
t28 = t20 * pkin(1);
t22 = cos(qJ(1));
t27 = t22 * pkin(1);
t18 = qJ(3) + qJ(4);
t14 = cos(t18);
t21 = cos(qJ(3));
t15 = t21 * pkin(3);
t26 = pkin(4) * t14 + t15;
t17 = qJ(1) + pkin(8);
t11 = sin(t17);
t12 = cos(t17);
t6 = g(2) * t12 + g(3) * t11;
t5 = g(2) * t11 - g(3) * t12;
t25 = g(2) * t22 + g(3) * t20;
t13 = sin(t18);
t2 = -g(1) * t14 - t5 * t13;
t19 = sin(qJ(3));
t24 = -g(1) * t21 - t5 * t19;
t16 = -qJ(5) + t23;
t10 = t15 + pkin(2);
t7 = pkin(2) + t26;
t4 = t6 * t14;
t3 = t6 * t13;
t1 = g(1) * t13 - t5 * t14;
t8 = [0, 0, 0, 0, 0, 0, t25, -g(2) * t20 + g(3) * t22, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t25 * pkin(1), 0, 0, 0, 0, 0, 0, t6 * t21, -t6 * t19, t5, -g(2) * (-t12 * pkin(2) - t11 * pkin(6) - t27) - g(3) * (-t11 * pkin(2) + t12 * pkin(6) - t28), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (-t12 * t10 + t11 * t23 - t27) - g(3) * (-t11 * t10 - t12 * t23 - t28), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t11 * t16 - t12 * t7 - t27) - g(3) * (-t11 * t7 - t12 * t16 - t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, g(1) * t19 - t5 * t21, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t24 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t26 + t5 * (-t19 * pkin(3) - pkin(4) * t13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t8;
