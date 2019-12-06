% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = sin(qJ(1));
t26 = t21 * pkin(1);
t22 = cos(qJ(1));
t25 = t22 * pkin(1);
t19 = cos(pkin(9));
t8 = t19 * pkin(3) + pkin(2);
t20 = -pkin(6) - qJ(3);
t16 = pkin(9) + qJ(4);
t17 = qJ(1) + pkin(8);
t10 = sin(t17);
t12 = cos(t17);
t4 = g(2) * t12 + g(3) * t10;
t3 = g(2) * t10 - g(3) * t12;
t24 = g(2) * t22 + g(3) * t21;
t11 = cos(t16);
t9 = sin(t16);
t23 = -g(1) * t11 - t3 * t9;
t15 = -pkin(7) + t20;
t13 = qJ(5) + t16;
t7 = cos(t13);
t6 = sin(t13);
t5 = pkin(4) * t11 + t8;
t2 = -g(1) * t7 - t3 * t6;
t1 = g(1) * t6 - t3 * t7;
t14 = [0, 0, 0, 0, 0, 0, t24, -g(2) * t21 + g(3) * t22, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t24 * pkin(1), 0, 0, 0, 0, 0, 0, t4 * t19, -t4 * sin(pkin(9)), t3, -g(2) * (-t12 * pkin(2) - t10 * qJ(3) - t25) - g(3) * (-t10 * pkin(2) + t12 * qJ(3) - t26), 0, 0, 0, 0, 0, 0, t4 * t11, -t4 * t9, t3, -g(2) * (t10 * t20 - t12 * t8 - t25) - g(3) * (-t10 * t8 - t12 * t20 - t26), 0, 0, 0, 0, 0, 0, t4 * t7, -t4 * t6, t3, -g(2) * (t10 * t15 - t12 * t5 - t25) - g(3) * (-t10 * t5 - t12 * t15 - t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(1) * t9 - t3 * t11, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t14;
