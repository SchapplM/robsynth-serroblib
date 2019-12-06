% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t24 = pkin(8) + qJ(2);
t21 = sin(t24);
t40 = pkin(2) * t21;
t26 = cos(pkin(9));
t39 = pkin(4) * t26;
t23 = qJ(3) + t24;
t19 = sin(t23);
t38 = g(1) * t19;
t25 = sin(pkin(9));
t37 = g(3) * t25;
t20 = cos(t23);
t36 = t20 * t25;
t27 = sin(qJ(5));
t35 = t26 * t27;
t28 = cos(qJ(5));
t34 = t26 * t28;
t33 = t20 * pkin(3) + t19 * qJ(4);
t15 = t20 * qJ(4);
t32 = -t19 * pkin(3) + t15;
t31 = pkin(7) * t36 + t20 * t39 + t33;
t9 = -g(2) * t20 + t38;
t22 = cos(t24);
t30 = g(1) * t21 - g(2) * t22;
t29 = (-pkin(7) * t25 - pkin(3) - t39) * t38;
t18 = pkin(2) * t22;
t10 = g(1) * t20 + g(2) * t19;
t8 = t9 * t26;
t7 = -g(2) * t36 + t25 * t38;
t6 = t19 * t27 + t20 * t34;
t5 = t19 * t28 - t20 * t35;
t4 = -t19 * t34 + t20 * t27;
t3 = t19 * t35 + t20 * t28;
t2 = -g(1) * t4 - g(2) * t6;
t1 = -g(1) * t3 - g(2) * t5;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(1) * t22 + g(2) * t21, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t30 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t32 - t40) - g(2) * (t18 + t33), 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(1) * (t15 - t40) - g(2) * (t18 + t31) - t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * t32 - g(2) * t33, 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(1) * t15 - g(2) * t31 - t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t27 * t37, g(1) * t6 - g(2) * t4 + t28 * t37, 0, 0;];
taug_reg = t11;
