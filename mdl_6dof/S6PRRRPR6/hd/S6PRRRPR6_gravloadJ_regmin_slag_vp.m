% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = sin(pkin(11));
t42 = sin(qJ(2));
t46 = cos(qJ(2));
t66 = cos(pkin(11));
t67 = cos(pkin(6));
t54 = t67 * t66;
t26 = t37 * t42 - t46 * t54;
t64 = t37 * t67;
t28 = t66 * t42 + t46 * t64;
t81 = g(1) * t28 + g(2) * t26;
t27 = t37 * t46 + t42 * t54;
t29 = -t42 * t64 + t66 * t46;
t41 = sin(qJ(3));
t45 = cos(qJ(3));
t38 = sin(pkin(6));
t63 = t38 * t66;
t72 = t38 * t45;
t73 = t38 * t42;
t49 = g(3) * (-t41 * t73 + t67 * t45) + g(2) * (-t27 * t41 - t45 * t63) + g(1) * (-t29 * t41 + t37 * t72);
t31 = t67 * t41 + t42 * t72;
t40 = sin(qJ(4));
t44 = cos(qJ(4));
t68 = t44 * t46;
t18 = t31 * t40 + t38 * t68;
t71 = t38 * t46;
t65 = t40 * t71;
t19 = t31 * t44 - t65;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t15 = t27 * t45 - t41 * t63;
t5 = t15 * t40 - t26 * t44;
t6 = t15 * t44 + t26 * t40;
t17 = t37 * t38 * t41 + t29 * t45;
t7 = t17 * t40 - t28 * t44;
t8 = t17 * t44 + t28 * t40;
t80 = g(1) * (t7 * t39 + t8 * t43) + g(2) * (t5 * t39 + t6 * t43) + g(3) * (t18 * t39 + t19 * t43);
t79 = g(1) * (t8 * t39 - t7 * t43) + g(2) * (t6 * t39 - t5 * t43) - g(3) * (t18 * t43 - t19 * t39);
t70 = t40 * t45;
t69 = t44 * t45;
t53 = pkin(3) * t45 + pkin(9) * t41 + pkin(2);
t1 = g(1) * t7 + g(2) * t5 + g(3) * t18;
t51 = g(1) * t8 + g(2) * t6 + g(3) * t19;
t10 = -t26 * t70 - t27 * t44;
t12 = -t28 * t70 - t29 * t44;
t20 = -t44 * t73 + t45 * t65;
t50 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t48 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t47 = g(3) * t71 - t81;
t21 = (t40 * t42 + t45 * t68) * t38;
t13 = -t28 * t69 + t29 * t40;
t11 = -t26 * t69 + t27 * t40;
t9 = t47 * t41;
t4 = t49 * t44;
t3 = t49 * t40;
t2 = -g(1) * t13 - g(2) * t11 - g(3) * t21;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t47, g(1) * t29 + g(2) * t27 + g(3) * t73, 0, 0, 0, 0, 0, -t47 * t45, t9, 0, 0, 0, 0, 0, t2, t50, t2, -t9, -t50, -g(1) * (t13 * pkin(4) + t29 * pkin(8) + t12 * qJ(5)) - g(2) * (t11 * pkin(4) + t27 * pkin(8) + t10 * qJ(5)) + t81 * t53 + (-t21 * pkin(4) - t20 * qJ(5) - (pkin(8) * t42 + t53 * t46) * t38) * g(3), 0, 0, 0, 0, 0, -g(1) * (t12 * t39 + t13 * t43) - g(2) * (t10 * t39 + t11 * t43) - g(3) * (t20 * t39 + t21 * t43) -g(1) * (t12 * t43 - t13 * t39) - g(2) * (t10 * t43 - t11 * t39) - g(3) * (t20 * t43 - t21 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t48, 0, 0, 0, 0, 0, -t4, t3, -t4, -t48, -t3, -t48 * pkin(9) - t49 * (pkin(4) * t44 + qJ(5) * t40 + pkin(3)) 0, 0, 0, 0, 0, -t49 * (t39 * t40 + t43 * t44) t49 * (t39 * t44 - t40 * t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t51, t1, 0, -t51, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - g(3) * (-t18 * pkin(4) + t19 * qJ(5)) 0, 0, 0, 0, 0, -t79, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t80;];
taug_reg  = t14;
