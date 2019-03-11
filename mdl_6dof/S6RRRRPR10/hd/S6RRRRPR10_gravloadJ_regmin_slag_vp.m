% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(2));
t46 = cos(qJ(1));
t62 = cos(pkin(6));
t58 = t46 * t62;
t25 = t41 * t58 + t42 * t45;
t37 = qJ(3) + qJ(4);
t35 = sin(t37);
t36 = cos(t37);
t38 = sin(pkin(6));
t65 = t38 * t46;
t13 = t25 * t35 + t36 * t65;
t24 = t42 * t41 - t45 * t58;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t82 = t13 * t39 + t24 * t43;
t81 = t13 * t43 - t24 * t39;
t59 = t42 * t62;
t27 = -t41 * t59 + t46 * t45;
t40 = sin(qJ(3));
t44 = cos(qJ(3));
t66 = t38 * t44;
t19 = -t27 * t40 + t42 * t66;
t53 = t25 * t40 + t44 * t65;
t68 = t38 * t41;
t80 = g(2) * t53 - g(3) * (-t40 * t68 + t62 * t44) - g(1) * t19;
t79 = g(1) * t46 + g(2) * t42;
t26 = t46 * t41 + t45 * t59;
t73 = g(3) * t38;
t51 = -g(1) * t26 - g(2) * t24 + t45 * t73;
t70 = t35 * t39;
t69 = t35 * t43;
t67 = t38 * t42;
t64 = t39 * t45;
t63 = t43 * t45;
t14 = t25 * t36 - t35 * t65;
t60 = t25 * t44 - t40 * t65;
t17 = t27 * t35 - t36 * t67;
t57 = -g(1) * t13 + g(2) * t17;
t18 = t27 * t36 + t35 * t67;
t56 = -g(1) * t14 + g(2) * t18;
t55 = g(1) * t24 - g(2) * t26;
t54 = g(1) * t27 + g(2) * t25;
t22 = t35 * t68 - t36 * t62;
t4 = g(1) * t17 + g(2) * t13 + g(3) * t22;
t23 = t35 * t62 + t36 * t68;
t6 = g(1) * t18 + g(2) * t14 + g(3) * t23;
t50 = g(3) * t68 + t54;
t49 = -g(1) * (-t17 * pkin(4) + t18 * qJ(5)) - g(2) * (-t13 * pkin(4) + t14 * qJ(5)) - g(3) * (-t22 * pkin(4) + t23 * qJ(5));
t47 = -pkin(10) - pkin(9);
t34 = t44 * pkin(3) + pkin(2);
t20 = t27 * t44 + t40 * t67;
t10 = t17 * t39 + t26 * t43;
t9 = t17 * t43 - t26 * t39;
t8 = t51 * t36;
t7 = t51 * t35;
t2 = t6 * t43;
t1 = t6 * t39;
t3 = [0, g(1) * t42 - g(2) * t46, t79, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t27, -t55, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t20, -g(1) * t53 - g(2) * t19, 0, 0, 0, 0, 0, -t56, t57, t55, t56, -t57, -g(1) * (-t42 * pkin(1) - pkin(4) * t14 - qJ(5) * t13 + t24 * t47 - t25 * t34) - g(2) * (t46 * pkin(1) + t18 * pkin(4) + t17 * qJ(5) - t26 * t47 + t27 * t34) - t79 * t38 * (pkin(3) * t40 + pkin(8)) 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t10, g(1) * t81 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, 0, 0, 0, 0, -t51 * t44, t51 * t40, 0, 0, 0, 0, 0, -t8, t7, -t50, t8, -t7 (t41 * t73 + t54) * t47 - t51 * (pkin(4) * t36 + qJ(5) * t35 + t34) 0, 0, 0, 0, 0, -g(1) * (-t26 * t70 + t27 * t43) - g(2) * (-t24 * t70 + t25 * t43) - (t35 * t64 + t41 * t43) * t73, -g(1) * (-t26 * t69 - t27 * t39) - g(2) * (-t24 * t69 - t25 * t39) - (t35 * t63 - t39 * t41) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, g(1) * t20 + g(2) * t60 - g(3) * (-t40 * t62 - t41 * t66) 0, 0, 0, 0, 0, t4, t6, 0, -t4, -t6, t80 * pkin(3) + t49, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, -t4, -t6, t49, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t81 - g(3) * (t22 * t43 + t38 * t64) g(1) * t10 + g(2) * t82 - g(3) * (-t22 * t39 + t38 * t63);];
taug_reg  = t3;
