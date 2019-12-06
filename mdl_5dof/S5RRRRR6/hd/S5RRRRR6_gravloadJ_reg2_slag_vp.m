% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t36 = -pkin(8) - pkin(7);
t35 = cos(qJ(1));
t47 = pkin(1) * t35;
t33 = sin(qJ(1));
t46 = t33 * pkin(1);
t30 = qJ(3) + qJ(4);
t25 = cos(t30);
t34 = cos(qJ(3));
t28 = t34 * pkin(3);
t45 = pkin(4) * t25 + t28;
t31 = qJ(1) + qJ(2);
t24 = sin(t31);
t26 = cos(t31);
t44 = -pkin(2) * t24 + t26 * pkin(7);
t13 = pkin(2) + t45;
t29 = -pkin(9) + t36;
t43 = -t13 * t26 + t24 * t29;
t22 = t28 + pkin(2);
t42 = -t22 * t26 + t24 * t36;
t41 = -pkin(2) * t26 - pkin(7) * t24;
t12 = g(2) * t26 + g(3) * t24;
t11 = g(2) * t24 - g(3) * t26;
t40 = g(2) * t35 + g(3) * t33;
t39 = -t13 * t24 - t26 * t29;
t38 = -t24 * t22 - t26 * t36;
t23 = sin(t30);
t4 = -g(1) * t25 - t11 * t23;
t32 = sin(qJ(3));
t37 = -g(1) * t34 - t11 * t32;
t27 = qJ(5) + t30;
t21 = cos(t27);
t20 = sin(t27);
t10 = t12 * t34;
t9 = t12 * t32;
t8 = t12 * t25;
t7 = t12 * t23;
t6 = t12 * t21;
t5 = t12 * t20;
t3 = g(1) * t23 - t11 * t25;
t2 = -g(1) * t21 - t11 * t20;
t1 = g(1) * t20 - t11 * t21;
t14 = [0, 0, 0, 0, 0, 0, t40, -g(2) * t33 + g(3) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t40 * pkin(1), 0, 0, 0, 0, 0, 0, t10, -t9, t11, -g(2) * (t41 - t47) - g(3) * (t44 - t46), 0, 0, 0, 0, 0, 0, t8, -t7, t11, -g(2) * (t42 - t47) - g(3) * (t38 - t46), 0, 0, 0, 0, 0, 0, t6, -t5, t11, -g(2) * (t43 - t47) - g(3) * (t39 - t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t11, -g(2) * t41 - g(3) * t44, 0, 0, 0, 0, 0, 0, t8, -t7, t11, -g(2) * t42 - g(3) * t38, 0, 0, 0, 0, 0, 0, t6, -t5, t11, -g(2) * t43 - g(3) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(1) * t32 - t11 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t37 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t45 + t11 * (-pkin(3) * t32 - pkin(4) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t14;
