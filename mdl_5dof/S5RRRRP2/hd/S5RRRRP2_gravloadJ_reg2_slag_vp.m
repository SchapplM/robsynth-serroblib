% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t29 = -pkin(8) - pkin(7);
t26 = sin(qJ(1));
t40 = t26 * pkin(1);
t28 = cos(qJ(1));
t39 = t28 * pkin(1);
t23 = qJ(3) + qJ(4);
t19 = cos(t23);
t27 = cos(qJ(3));
t21 = t27 * pkin(3);
t38 = pkin(4) * t19 + t21;
t24 = qJ(1) + qJ(2);
t18 = sin(t24);
t20 = cos(t24);
t37 = -t18 * pkin(2) + t20 * pkin(7);
t22 = -qJ(5) + t29;
t9 = pkin(2) + t38;
t36 = t18 * t22 - t20 * t9;
t16 = t21 + pkin(2);
t35 = -t20 * t16 + t18 * t29;
t34 = -t20 * pkin(2) - t18 * pkin(7);
t8 = g(2) * t20 + g(3) * t18;
t7 = g(2) * t18 - g(3) * t20;
t33 = g(2) * t28 + g(3) * t26;
t32 = -t18 * t9 - t20 * t22;
t31 = -t18 * t16 - t20 * t29;
t17 = sin(t23);
t2 = -g(1) * t19 - t7 * t17;
t25 = sin(qJ(3));
t30 = -g(1) * t27 - t7 * t25;
t6 = t8 * t27;
t5 = t8 * t25;
t4 = t8 * t19;
t3 = t8 * t17;
t1 = g(1) * t17 - t7 * t19;
t10 = [0, 0, 0, 0, 0, 0, t33, -g(2) * t26 + g(3) * t28, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * (t34 - t39) - g(3) * (t37 - t40), 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * (t35 - t39) - g(3) * (t31 - t40), 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * (t36 - t39) - g(3) * (t32 - t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * t34 - g(3) * t37, 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t35 - g(3) * t31, 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t36 - g(3) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(1) * t25 - t7 * t27, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t30 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t38 + t7 * (-t25 * pkin(3) - pkin(4) * t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t10;
