% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t35 = pkin(2) * t18;
t19 = cos(t21);
t34 = pkin(2) * t19;
t25 = sin(qJ(1));
t33 = t25 * pkin(1);
t26 = cos(qJ(1));
t32 = t26 * pkin(1);
t17 = pkin(8) + t21;
t12 = sin(t17);
t13 = cos(t17);
t6 = g(2) * t13 + g(3) * t12;
t5 = g(2) * t12 - g(3) * t13;
t8 = g(2) * t19 + g(3) * t18;
t31 = g(2) * t26 + g(3) * t25;
t30 = -t12 * pkin(3) + t13 * qJ(4) - t35;
t23 = cos(pkin(9));
t14 = t23 * pkin(4) + pkin(3);
t24 = -pkin(7) - qJ(4);
t29 = t12 * t24 - t13 * t14 - t34;
t28 = -t13 * pkin(3) - t12 * qJ(4) - t34;
t27 = -t12 * t14 - t13 * t24 - t35;
t20 = pkin(9) + qJ(5);
t16 = cos(t20);
t15 = sin(t20);
t7 = -g(2) * t18 + g(3) * t19;
t4 = t6 * t23;
t3 = t6 * sin(pkin(9));
t2 = t6 * t16;
t1 = t6 * t15;
t9 = [0, 0, 0, 0, 0, 0, t31, -g(2) * t25 + g(3) * t26, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t31 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, 0, -g(2) * (-t32 - t34) - g(3) * (-t33 - t35), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t28 - t32) - g(3) * (t30 - t33), 0, 0, 0, 0, 0, 0, t2, -t1, t5, -g(2) * (t29 - t32) - g(3) * (t27 - t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * t28 - g(3) * t30, 0, 0, 0, 0, 0, 0, t2, -t1, t5, -g(2) * t29 - g(3) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 - t5 * t15, g(1) * t15 - t5 * t16, 0, 0;];
taug_reg = t9;
