% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = qJ(3) + pkin(9);
t15 = cos(t20);
t23 = cos(qJ(3));
t17 = t23 * pkin(3);
t25 = pkin(4) * t15 + t17;
t21 = -qJ(4) - pkin(6);
t19 = pkin(8) + qJ(2);
t12 = sin(t19);
t14 = cos(t19);
t4 = g(1) * t14 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t14;
t22 = sin(qJ(3));
t24 = -g(3) * t23 + t4 * t22;
t18 = -pkin(7) + t21;
t16 = qJ(5) + t20;
t13 = sin(t20);
t11 = t17 + pkin(2);
t10 = cos(t16);
t9 = sin(t16);
t5 = pkin(2) + t25;
t2 = g(3) * t9 + t4 * t10;
t1 = -g(3) * t10 + t4 * t9;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t23, -t3 * t22, -t4, -g(1) * (-t12 * pkin(2) + t14 * pkin(6)) - g(2) * (t14 * pkin(2) + t12 * pkin(6)), 0, 0, 0, 0, 0, 0, t3 * t15, -t3 * t13, -t4, -g(1) * (-t12 * t11 - t14 * t21) - g(2) * (t14 * t11 - t12 * t21), 0, 0, 0, 0, 0, 0, t3 * t10, -t3 * t9, -t4, -g(1) * (-t12 * t5 - t14 * t18) - g(2) * (-t12 * t18 + t14 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, g(3) * t22 + t4 * t23, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t15 + t4 * t13, g(3) * t13 + t4 * t15, 0, t24 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t25 - t4 * (-t22 * pkin(3) - pkin(4) * t13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
