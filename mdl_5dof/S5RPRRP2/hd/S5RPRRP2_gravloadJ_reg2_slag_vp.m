% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP2
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = qJ(1) + pkin(8);
t15 = qJ(3) + t16;
t10 = sin(t15);
t11 = cos(t15);
t28 = -t10 * pkin(3) + t11 * pkin(7);
t20 = cos(qJ(4));
t12 = t20 * pkin(4) + pkin(3);
t17 = -qJ(5) - pkin(7);
t27 = t10 * t17 - t11 * t12;
t13 = sin(t16);
t19 = sin(qJ(1));
t26 = -t19 * pkin(1) - pkin(2) * t13;
t14 = cos(t16);
t21 = cos(qJ(1));
t25 = -t21 * pkin(1) - pkin(2) * t14;
t24 = -t11 * pkin(3) - t10 * pkin(7);
t6 = g(2) * t11 + g(3) * t10;
t5 = g(2) * t10 - g(3) * t11;
t23 = g(2) * t21 + g(3) * t19;
t22 = -t10 * t12 - t11 * t17;
t18 = sin(qJ(4));
t2 = -g(1) * t20 - t5 * t18;
t4 = t6 * t20;
t3 = t6 * t18;
t1 = g(1) * t18 - t5 * t20;
t7 = [0, 0, 0, 0, 0, 0, t23, -g(2) * t19 + g(3) * t21, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t14 + g(3) * t13, -g(2) * t13 + g(3) * t14, 0, t23 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, 0, -g(2) * t25 - g(3) * t26, 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t24 + t25) - g(3) * (t26 + t28), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t25 + t27) - g(3) * (t22 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * t24 - g(3) * t28, 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * t27 - g(3) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t7;
