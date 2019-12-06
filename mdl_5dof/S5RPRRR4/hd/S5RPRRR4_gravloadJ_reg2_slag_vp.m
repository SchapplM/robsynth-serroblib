% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR4
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t17 = qJ(1) + pkin(9);
t16 = qJ(3) + t17;
t11 = sin(t16);
t30 = pkin(3) * t11;
t12 = cos(t16);
t29 = pkin(3) * t12;
t13 = qJ(4) + t16;
t10 = cos(t13);
t9 = sin(t13);
t28 = -t9 * pkin(4) + t10 * pkin(8);
t27 = -t10 * pkin(4) - t9 * pkin(8);
t3 = g(2) * t9 - g(3) * t10;
t4 = g(2) * t10 + g(3) * t9;
t14 = sin(t17);
t19 = sin(qJ(1));
t26 = -t19 * pkin(1) - pkin(2) * t14;
t15 = cos(t17);
t21 = cos(qJ(1));
t25 = -t21 * pkin(1) - pkin(2) * t15;
t6 = g(2) * t12 + g(3) * t11;
t24 = g(2) * t21 + g(3) * t19;
t23 = t26 - t30;
t22 = t25 - t29;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t5 = -g(2) * t11 + g(3) * t12;
t2 = t4 * t20;
t1 = t4 * t18;
t7 = [0, 0, 0, 0, 0, 0, t24, -g(2) * t19 + g(3) * t21, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t15 + g(3) * t14, -g(2) * t14 + g(3) * t15, 0, t24 * pkin(1), 0, 0, 0, 0, 0, 0, t6, t5, 0, -g(2) * t25 - g(3) * t26, 0, 0, 0, 0, 0, 0, t4, -t3, 0, -g(2) * t22 - g(3) * t23, 0, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * (t22 + t27) - g(3) * (t23 + t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t6 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * (t27 - t29) - g(3) * (t28 - t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * t27 - g(3) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 - t3 * t18, g(1) * t18 - t3 * t20, 0, 0;];
taug_reg = t7;
