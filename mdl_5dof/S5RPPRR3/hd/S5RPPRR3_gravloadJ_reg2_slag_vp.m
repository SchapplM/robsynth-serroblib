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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
t22 = cos(pkin(9));
t9 = t22 * pkin(3) + pkin(2);
t23 = -pkin(6) - qJ(3);
t19 = pkin(9) + qJ(4);
t20 = qJ(1) + pkin(8);
t11 = sin(t20);
t13 = cos(t20);
t4 = g(2) * t13 + g(3) * t11;
t3 = g(2) * t11 - g(3) * t13;
t24 = sin(qJ(1));
t25 = cos(qJ(1));
t27 = -g(2) * t25 - g(3) * t24;
t10 = sin(t19);
t12 = cos(t19);
t26 = -g(1) * t12 + t3 * t10;
t18 = -pkin(7) + t23;
t17 = t25 * pkin(1);
t16 = t24 * pkin(1);
t14 = qJ(5) + t19;
t8 = cos(t14);
t7 = sin(t14);
t5 = pkin(4) * t12 + t9;
t2 = -g(1) * t8 + t3 * t7;
t1 = g(1) * t7 + t3 * t8;
t6 = [0, 0, 0, 0, 0, 0, t27, g(2) * t24 - g(3) * t25, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, -t4 * t22, t4 * sin(pkin(9)), -t3, -g(2) * (t13 * pkin(2) + t11 * qJ(3) + t17) - g(3) * (t11 * pkin(2) - t13 * qJ(3) + t16), 0, 0, 0, 0, 0, 0, -t4 * t12, t4 * t10, -t3, -g(2) * (-t11 * t23 + t13 * t9 + t17) - g(3) * (t11 * t9 + t13 * t23 + t16), 0, 0, 0, 0, 0, 0, -t4 * t8, t4 * t7, -t3, -g(2) * (-t11 * t18 + t13 * t5 + t17) - g(3) * (t11 * t5 + t13 * t18 + t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, g(1) * t10 + t3 * t12, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t26 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t6;
