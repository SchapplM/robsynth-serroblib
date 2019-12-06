% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR3
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t32 = sin(qJ(1));
t45 = t32 * pkin(1);
t34 = cos(qJ(1));
t44 = t34 * pkin(1);
t30 = -qJ(4) - pkin(7);
t28 = qJ(3) + pkin(9);
t22 = cos(t28);
t33 = cos(qJ(3));
t26 = t33 * pkin(3);
t43 = pkin(4) * t22 + t26;
t29 = qJ(1) + qJ(2);
t24 = sin(t29);
t25 = cos(t29);
t42 = -t24 * pkin(2) + t25 * pkin(7);
t11 = pkin(2) + t43;
t27 = -pkin(8) + t30;
t41 = -t25 * t11 + t24 * t27;
t20 = t26 + pkin(2);
t40 = -t25 * t20 + t24 * t30;
t39 = -t25 * pkin(2) - t24 * pkin(7);
t10 = g(2) * t25 + g(3) * t24;
t9 = g(2) * t24 - g(3) * t25;
t38 = g(2) * t34 + g(3) * t32;
t37 = -t24 * t11 - t25 * t27;
t36 = -t24 * t20 - t25 * t30;
t31 = sin(qJ(3));
t35 = -g(1) * t33 - t31 * t9;
t23 = qJ(5) + t28;
t21 = sin(t28);
t17 = cos(t23);
t16 = sin(t23);
t8 = t10 * t33;
t7 = t10 * t31;
t6 = t10 * t22;
t5 = t10 * t21;
t4 = t10 * t17;
t3 = t10 * t16;
t2 = -g(1) * t17 - t16 * t9;
t1 = g(1) * t16 - t17 * t9;
t12 = [0, 0, 0, 0, 0, 0, t38, -g(2) * t32 + g(3) * t34, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * (t39 - t44) - g(3) * (t42 - t45), 0, 0, 0, 0, 0, 0, t6, -t5, t9, -g(2) * (t40 - t44) - g(3) * (t36 - t45), 0, 0, 0, 0, 0, 0, t4, -t3, t9, -g(2) * (t41 - t44) - g(3) * (t37 - t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * t39 - g(3) * t42, 0, 0, 0, 0, 0, 0, t6, -t5, t9, -g(2) * t40 - g(3) * t36, 0, 0, 0, 0, 0, 0, t4, -t3, t9, -g(2) * t41 - g(3) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(1) * t31 - t33 * t9, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 - t21 * t9, g(1) * t21 - t22 * t9, 0, t35 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t43 + t9 * (-t31 * pkin(3) - pkin(4) * t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t12;
