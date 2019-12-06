% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t18 = qJ(1) + qJ(2);
t16 = sin(t18);
t32 = pkin(2) * t16;
t17 = cos(t18);
t31 = pkin(2) * t17;
t21 = sin(qJ(1));
t30 = t21 * pkin(1);
t23 = cos(qJ(1));
t29 = t23 * pkin(1);
t15 = pkin(8) + t18;
t12 = sin(t15);
t13 = cos(t15);
t6 = g(2) * t13 + g(3) * t12;
t5 = g(2) * t12 - g(3) * t13;
t8 = g(2) * t17 + g(3) * t16;
t28 = g(2) * t23 + g(3) * t21;
t27 = -t12 * pkin(3) + t13 * pkin(7) - t32;
t22 = cos(qJ(4));
t14 = t22 * pkin(4) + pkin(3);
t19 = -qJ(5) - pkin(7);
t26 = t12 * t19 - t13 * t14 - t31;
t25 = -t13 * pkin(3) - t12 * pkin(7) - t31;
t24 = -t12 * t14 - t13 * t19 - t32;
t20 = sin(qJ(4));
t2 = -g(1) * t22 - t5 * t20;
t7 = -g(2) * t16 + g(3) * t17;
t4 = t6 * t22;
t3 = t6 * t20;
t1 = g(1) * t20 - t5 * t22;
t9 = [0, 0, 0, 0, 0, 0, t28, -g(2) * t21 + g(3) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, 0, -g(2) * (-t29 - t31) - g(3) * (-t30 - t32), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t25 - t29) - g(3) * (t27 - t30), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t26 - t29) - g(3) * (t24 - t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * t25 - g(3) * t27, 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * t26 - g(3) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6;];
taug_reg = t9;
