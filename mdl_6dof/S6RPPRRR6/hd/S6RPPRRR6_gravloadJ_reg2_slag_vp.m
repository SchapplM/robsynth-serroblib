% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t14 = g(1) * t30 + g(2) * t27;
t26 = sin(qJ(4));
t25 = sin(qJ(5));
t46 = t30 * t25;
t28 = cos(qJ(5));
t50 = t27 * t28;
t10 = -t26 * t46 - t50;
t29 = cos(qJ(4));
t56 = g(3) * t29;
t45 = t30 * t28;
t51 = t27 * t25;
t8 = t26 * t51 - t45;
t63 = -g(1) * t10 + g(2) * t8 + t25 * t56;
t33 = -g(3) * t26 + t14 * t29;
t60 = g(1) * t27;
t55 = t26 * pkin(4);
t54 = t30 * pkin(7);
t24 = qJ(5) + qJ(6);
t17 = sin(t24);
t53 = t27 * t17;
t18 = cos(t24);
t52 = t27 * t18;
t49 = t29 * t30;
t48 = t30 * t17;
t47 = t30 * t18;
t44 = -pkin(1) - qJ(3);
t43 = t30 * pkin(1) + t27 * qJ(2);
t41 = t30 * qJ(3) + t43;
t40 = pkin(5) * t25 + pkin(7);
t39 = g(2) * t41;
t37 = t29 * pkin(8) - t55;
t13 = -g(2) * t30 + t60;
t21 = t30 * qJ(2);
t36 = t44 * t27 + t21;
t16 = t28 * pkin(5) + pkin(4);
t31 = -pkin(9) - pkin(8);
t34 = t26 * t16 + t29 * t31;
t12 = -g(2) * t49 + t29 * t60;
t11 = t26 * t45 - t51;
t9 = -t26 * t50 - t46;
t7 = t14 * t26 + t56;
t6 = t26 * t47 - t53;
t5 = -t26 * t48 - t52;
t4 = -t26 * t52 - t48;
t3 = t26 * t53 - t47;
t2 = g(1) * t6 - g(2) * t4 + t18 * t56;
t1 = -g(1) * t5 + g(2) * t3 + t17 * t56;
t15 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -g(1) * (-t27 * pkin(1) + t21) - g(2) * t43, 0, 0, 0, 0, 0, 0, 0, -t14, t13, -g(1) * t36 - t39, 0, 0, 0, 0, 0, 0, t13 * t26, t12, t14, -g(1) * (t36 - t54) - g(2) * (-t27 * pkin(7) + t41) 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, -t12, -g(1) * (t21 - t54) - g(2) * (-pkin(8) * t49 + t30 * t55 + t41) + (-g(1) * (t37 + t44) + g(2) * pkin(7)) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, -t12, -g(1) * t21 - t39 + (g(1) * t40 - g(2) * t34) * t30 + (-g(1) * (-t34 + t44) + g(2) * t40) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t28, t33 * t25, -t7, -g(3) * t37 - t14 * (pkin(4) * t29 + pkin(8) * t26) 0, 0, 0, 0, 0, 0, -t33 * t18, t33 * t17, -t7, g(3) * t34 - t14 * (t16 * t29 - t26 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, g(1) * t11 - g(2) * t9 + t28 * t56, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t63 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;
