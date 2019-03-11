% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t59 = pkin(5) * t37 + qJ(6) * t35;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t19 = g(1) * t38 + g(2) * t36;
t32 = pkin(10) + qJ(3);
t29 = qJ(4) + t32;
t23 = sin(t29);
t41 = t19 * t23;
t24 = cos(t29);
t58 = t24 * pkin(4) + t23 * pkin(9);
t27 = sin(t32);
t57 = pkin(3) * t27;
t55 = pkin(9) * t24;
t52 = g(3) * t23;
t51 = g(3) * t35;
t50 = t36 * t35;
t49 = t36 * t37;
t48 = t38 * t35;
t47 = t38 * t37;
t45 = t59 * t24 + t58;
t10 = t24 * t48 - t49;
t8 = t24 * t50 + t47;
t44 = g(1) * t8 - g(2) * t10;
t18 = g(1) * t36 - g(2) * t38;
t28 = cos(t32);
t22 = pkin(3) * t28;
t34 = cos(pkin(10));
t43 = t34 * pkin(2) + pkin(1) + t22 + t58;
t1 = g(1) * t10 + g(2) * t8 + t23 * t51;
t11 = t24 * t47 + t50;
t9 = t24 * t49 - t48;
t40 = g(1) * t11 + g(2) * t9 + t37 * t52;
t5 = -g(3) * t24 + t41;
t39 = (pkin(4) + t59) * t41;
t31 = -pkin(8) - pkin(7) - qJ(2);
t16 = t38 * t55;
t14 = t36 * t55;
t7 = t18 * t23;
t6 = t19 * t24 + t52;
t4 = t5 * t37;
t3 = -t24 * t51 + t35 * t41;
t2 = g(1) * t9 - g(2) * t11;
t12 = [0, t18, t19, t18 * t34, -t18 * sin(pkin(10)) -t19, -g(1) * (-t36 * pkin(1) + t38 * qJ(2)) - g(2) * (t38 * pkin(1) + t36 * qJ(2)) 0, 0, 0, 0, 0, t18 * t28, -t18 * t27, 0, 0, 0, 0, 0, t18 * t24, -t7, 0, 0, 0, 0, 0, t2, -t44, t2, t7, t44, -g(1) * (-t9 * pkin(5) - t8 * qJ(6)) - g(2) * (t11 * pkin(5) + t10 * qJ(6)) + (g(1) * t31 - g(2) * t43) * t38 + (g(1) * t43 + g(2) * t31) * t36; 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t28 + t19 * t27, g(3) * t27 + t19 * t28, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-t38 * t57 + t16) - g(2) * (-t36 * t57 + t14) - g(3) * (t22 + t45) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t16 - g(2) * t14 - g(3) * t45 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t40, t1, 0, -t40, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - (-pkin(5) * t35 + qJ(6) * t37) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
