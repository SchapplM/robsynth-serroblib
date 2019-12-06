% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t18 = pkin(8) + qJ(2);
t15 = sin(t18);
t27 = pkin(2) * t15;
t17 = qJ(3) + t18;
t12 = sin(t17);
t13 = cos(t17);
t26 = t13 * pkin(3) + t12 * pkin(7);
t25 = -t12 * pkin(3) + t13 * pkin(7);
t21 = cos(qJ(4));
t14 = t21 * pkin(4) + pkin(3);
t19 = -qJ(5) - pkin(7);
t24 = -t12 * t19 + t13 * t14;
t6 = g(1) * t13 + g(2) * t12;
t5 = g(1) * t12 - g(2) * t13;
t16 = cos(t18);
t23 = g(1) * t15 - g(2) * t16;
t22 = -t12 * t14 - t13 * t19;
t20 = sin(qJ(4));
t1 = -g(3) * t21 + t6 * t20;
t11 = pkin(2) * t16;
t4 = t5 * t21;
t3 = t5 * t20;
t2 = g(3) * t20 + t6 * t21;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(1) * t16 + g(2) * t15, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t23 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t25 - t27) - g(2) * (t11 + t26), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t22 - t27) - g(2) * (t11 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t25 - g(2) * t26, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t22 - g(2) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t7;
