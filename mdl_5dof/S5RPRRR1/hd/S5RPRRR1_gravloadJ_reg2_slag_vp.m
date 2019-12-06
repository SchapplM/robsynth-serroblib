% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_gravloadJ_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = sin(qJ(1));
t23 = cos(qJ(1));
t44 = -g(1) * t23 - g(2) * t19;
t18 = sin(qJ(3));
t41 = g(3) * t18;
t16 = sin(qJ(5));
t40 = t18 * t16;
t20 = cos(qJ(5));
t39 = t18 * t20;
t21 = cos(qJ(4));
t38 = t18 * t21;
t37 = t18 * t23;
t22 = cos(qJ(3));
t36 = t19 * t22;
t35 = t22 * t16;
t34 = t22 * t20;
t17 = sin(qJ(4));
t33 = t23 * t17;
t32 = t23 * t21;
t5 = t17 * t36 + t32;
t7 = -t19 * t21 + t22 * t33;
t31 = g(1) * t5 - g(2) * t7;
t11 = g(1) * t19 - g(2) * t23;
t6 = t21 * t36 - t33;
t30 = -t6 * t16 + t19 * t39;
t29 = t19 * t40 + t6 * t20;
t28 = t20 * t38 - t35;
t27 = t16 * t38 + t34;
t26 = g(1) * t7 + g(2) * t5 + t17 * t41;
t8 = t19 * t17 + t22 * t32;
t25 = g(1) * t8 + g(2) * t6 + g(3) * t38;
t24 = -g(3) * t22 - t18 * t44;
t10 = t44 * qJ(2);
t9 = t11 * t18;
t4 = -t22 * t44 + t41;
t3 = t24 * t17;
t2 = t16 * t37 + t8 * t20;
t1 = -t8 * t16 + t20 * t37;
t12 = [0, 0, 0, 0, 0, 0, t11, -t44, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t44, t10, 0, 0, 0, 0, 0, 0, t11 * t22, -t9, t44, t10, 0, 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, -t31, t9, t10, 0, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t2, g(1) * t30 - g(2) * t1, t31, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t4, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t21, -t3, -t4, 0, 0, 0, 0, 0, 0, 0, -g(3) * (t21 * t34 + t40) - t44 * t28, -g(3) * (-t21 * t35 + t39) + t44 * t27, t3, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t25, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t20, -t26 * t16, -t25, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t30 + g(3) * t27, g(1) * t2 + g(2) * t29 + g(3) * t28, 0, 0;];
taug_reg = t12;
