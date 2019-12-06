% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = pkin(9) + qJ(2);
t17 = sin(t22);
t33 = pkin(2) * t17;
t19 = qJ(3) + t22;
t14 = sin(t19);
t15 = cos(t19);
t32 = t15 * pkin(3) + t14 * pkin(7);
t31 = -t14 * pkin(3) + t15 * pkin(7);
t25 = cos(qJ(4));
t16 = t25 * pkin(4) + pkin(3);
t26 = -pkin(8) - pkin(7);
t30 = -t14 * t26 + t15 * t16;
t8 = g(1) * t15 + g(2) * t14;
t7 = g(1) * t14 - g(2) * t15;
t18 = cos(t22);
t29 = g(1) * t17 - g(2) * t18;
t28 = -t14 * t16 - t15 * t26;
t24 = sin(qJ(4));
t27 = -g(3) * t25 + t8 * t24;
t23 = qJ(4) + qJ(5);
t21 = cos(t23);
t20 = sin(t23);
t13 = pkin(2) * t18;
t6 = t7 * t25;
t5 = t7 * t24;
t4 = t7 * t21;
t3 = t7 * t20;
t2 = g(3) * t20 + t8 * t21;
t1 = -g(3) * t21 + t8 * t20;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(1) * t18 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t29 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t31 - t33) - g(2) * (t13 + t32), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t28 - t33) - g(2) * (t13 + t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t31 - g(2) * t32, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t28 - g(2) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t24 + t8 * t25, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
