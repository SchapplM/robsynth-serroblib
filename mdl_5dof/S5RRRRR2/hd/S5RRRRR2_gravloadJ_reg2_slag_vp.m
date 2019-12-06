% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_gravloadJ_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = qJ(1) + qJ(2);
t23 = sin(t27);
t25 = cos(t27);
t18 = g(1) * t23 - g(2) * t25;
t19 = g(1) * t25 + g(2) * t23;
t26 = qJ(3) + qJ(4);
t22 = sin(t26);
t24 = cos(t26);
t6 = -g(3) * t24 + t19 * t22;
t32 = cos(qJ(3));
t44 = pkin(2) * t32;
t41 = g(3) * t22;
t28 = sin(qJ(5));
t39 = t23 * t28;
t31 = cos(qJ(5));
t38 = t23 * t31;
t37 = t25 * t28;
t36 = t25 * t31;
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t35 = g(1) * t30 - g(2) * t33;
t17 = t18 * t32;
t29 = sin(qJ(3));
t34 = -g(3) * t32 + t19 * t29;
t20 = t35 * pkin(1);
t16 = t18 * t29;
t15 = pkin(2) * t17;
t14 = t24 * t36 + t39;
t13 = -t24 * t37 + t38;
t12 = -t24 * t38 + t37;
t11 = t24 * t39 + t36;
t10 = t18 * t24;
t9 = t18 * t22;
t8 = t34 * pkin(2);
t7 = t19 * t24 + t41;
t5 = -g(1) * (-t30 * pkin(1) - t23 * t44) - g(2) * (t33 * pkin(1) + t25 * t44);
t4 = t6 * t31;
t3 = t6 * t28;
t2 = -g(1) * t12 - g(2) * t14;
t1 = -g(1) * t11 - g(2) * t13;
t21 = [0, 0, 0, 0, 0, 0, t35, g(1) * t33 + g(2) * t30, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, t20, 0, 0, 0, 0, 0, 0, t17, -t16, -t19, t20, 0, 0, 0, 0, 0, 0, t10, -t9, -t19, t5, 0, 0, 0, 0, 0, 0, t2, t1, t9, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t19, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t19, t15, 0, 0, 0, 0, 0, 0, t2, t1, t9, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t29 + t19 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t8, 0, 0, 0, 0, 0, 0, t4, -t3, -t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t7, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t28 * t41, g(1) * t14 - g(2) * t12 + t31 * t41, 0, 0;];
taug_reg = t21;
