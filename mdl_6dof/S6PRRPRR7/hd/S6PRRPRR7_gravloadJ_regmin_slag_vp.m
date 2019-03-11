% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = sin(pkin(11));
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t41 = cos(pkin(11));
t42 = cos(pkin(6));
t37 = t42 * t41;
t10 = t24 * t28 - t31 * t37;
t40 = t24 * t42;
t12 = t41 * t28 + t31 * t40;
t55 = -g(1) * t12 - g(2) * t10;
t25 = sin(pkin(6));
t52 = g(3) * t25;
t23 = qJ(5) + qJ(6);
t21 = sin(t23);
t27 = sin(qJ(3));
t51 = t21 * t27;
t22 = cos(t23);
t50 = t22 * t27;
t49 = t25 * t28;
t30 = cos(qJ(3));
t48 = t25 * t30;
t47 = t25 * t31;
t26 = sin(qJ(5));
t46 = t26 * t27;
t29 = cos(qJ(5));
t45 = t27 * t29;
t44 = t27 * t31;
t43 = t29 * t31;
t39 = t25 * t41;
t11 = t24 * t31 + t28 * t37;
t13 = -t28 * t40 + t41 * t31;
t38 = -g(1) * t13 - g(2) * t11;
t14 = t27 * t49 - t42 * t30;
t6 = t11 * t27 + t30 * t39;
t8 = t13 * t27 - t24 * t48;
t35 = g(1) * t8 + g(2) * t6 + g(3) * t14;
t15 = t42 * t27 + t28 * t48;
t7 = t11 * t30 - t27 * t39;
t9 = t24 * t25 * t27 + t13 * t30;
t34 = g(1) * t9 + g(2) * t7 + g(3) * t15;
t33 = g(3) * t47 + t55;
t32 = g(3) * t49 - t38;
t5 = t33 * t30;
t4 = t33 * t27;
t2 = -g(1) * (-t12 * t22 - t8 * t21) - g(2) * (-t10 * t22 - t6 * t21) - g(3) * (-t14 * t21 + t22 * t47);
t1 = -g(1) * (-t12 * t21 + t8 * t22) - g(2) * (-t10 * t21 + t6 * t22) - g(3) * (t14 * t22 + t21 * t47);
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t33, t32, 0, 0, 0, 0, 0, -t5, t4, -t32, t5, -t4 (-t28 * t52 + t38) * pkin(8) + (-t31 * t52 - t55) * (pkin(3) * t30 + qJ(4) * t27 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t12 * t46 + t13 * t29) - g(2) * (-t10 * t46 + t11 * t29) - (t26 * t44 + t28 * t29) * t52, -g(1) * (-t12 * t45 - t13 * t26) - g(2) * (-t10 * t45 - t11 * t26) - (-t26 * t28 + t27 * t43) * t52, 0, 0, 0, 0, 0, -g(1) * (-t12 * t51 + t13 * t22) - g(2) * (-t10 * t51 + t11 * t22) - (t21 * t44 + t22 * t28) * t52, -g(1) * (-t12 * t50 - t13 * t21) - g(2) * (-t10 * t50 - t11 * t21) - (-t21 * t28 + t22 * t44) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, 0, -t35, -t34, -g(1) * (-t8 * pkin(3) + t9 * qJ(4)) - g(2) * (-t6 * pkin(3) + t7 * qJ(4)) - g(3) * (-t14 * pkin(3) + t15 * qJ(4)) 0, 0, 0, 0, 0, -t34 * t26, -t34 * t29, 0, 0, 0, 0, 0, -t34 * t21, -t34 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t26 + t8 * t29) - g(2) * (-t10 * t26 + t6 * t29) - g(3) * (t14 * t29 + t26 * t47) -g(1) * (-t12 * t29 - t8 * t26) - g(2) * (-t10 * t29 - t6 * t26) - g(3) * (-t14 * t26 + t25 * t43) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
