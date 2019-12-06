% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = -pkin(7) - pkin(6);
t18 = qJ(3) + qJ(4);
t14 = cos(t18);
t20 = cos(qJ(3));
t15 = t20 * pkin(3);
t23 = pkin(4) * t14 + t15;
t17 = pkin(8) + qJ(2);
t11 = sin(t17);
t12 = cos(t17);
t6 = g(1) * t12 + g(2) * t11;
t5 = g(1) * t11 - g(2) * t12;
t13 = sin(t18);
t1 = -g(3) * t14 + t6 * t13;
t19 = sin(qJ(3));
t22 = -g(3) * t20 + t6 * t19;
t16 = -qJ(5) + t21;
t10 = t15 + pkin(2);
t7 = pkin(2) + t23;
t4 = t5 * t14;
t3 = t5 * t13;
t2 = g(3) * t13 + t6 * t14;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t20, -t5 * t19, -t6, -g(1) * (-t11 * pkin(2) + t12 * pkin(6)) - g(2) * (t12 * pkin(2) + t11 * pkin(6)), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t11 * t10 - t12 * t21) - g(2) * (t12 * t10 - t11 * t21), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t11 * t7 - t12 * t16) - g(2) * (-t11 * t16 + t12 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t19 + t6 * t20, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t22 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t23 - t6 * (-t19 * pkin(3) - pkin(4) * t13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t8;
