% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = qJ(1) + qJ(2);
t16 = sin(t19);
t36 = pkin(2) * t16;
t17 = cos(t19);
t35 = pkin(2) * t17;
t18 = qJ(3) + t19;
t14 = sin(t18);
t34 = pkin(3) * t14;
t15 = cos(t18);
t33 = pkin(3) * t15;
t21 = sin(qJ(1));
t32 = t21 * pkin(1);
t23 = cos(qJ(1));
t31 = t23 * pkin(1);
t30 = -t34 - t36;
t29 = -t33 - t35;
t13 = pkin(9) + t18;
t11 = sin(t13);
t12 = cos(t13);
t4 = g(2) * t12 + g(3) * t11;
t3 = g(2) * t11 - g(3) * t12;
t6 = g(2) * t15 + g(3) * t14;
t8 = g(2) * t17 + g(3) * t16;
t28 = g(2) * t23 + g(3) * t21;
t27 = -t11 * pkin(4) + t12 * pkin(8) - t34;
t26 = -t12 * pkin(4) - t11 * pkin(8) - t33;
t25 = t27 - t36;
t24 = t26 - t35;
t22 = cos(qJ(5));
t20 = sin(qJ(5));
t7 = -g(2) * t16 + g(3) * t17;
t5 = -g(2) * t14 + g(3) * t15;
t2 = t4 * t22;
t1 = t4 * t20;
t9 = [0, 0, 0, 0, 0, 0, t28, -g(2) * t21 + g(3) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t6, t5, 0, -g(2) * (-t31 - t35) - g(3) * (-t32 - t36), 0, 0, 0, 0, 0, 0, t4, -t3, 0, -g(2) * (t29 - t31) - g(3) * (t30 - t32), 0, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * (t24 - t31) - g(3) * (t25 - t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, 0, -g(2) * t29 - g(3) * t30, 0, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * t24 - g(3) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t6 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * t26 - g(3) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 - t3 * t20, g(1) * t20 - t3 * t22, 0, 0;];
taug_reg = t9;
