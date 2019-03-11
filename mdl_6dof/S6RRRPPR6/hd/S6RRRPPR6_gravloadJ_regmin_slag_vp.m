% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t49 = sin(qJ(2));
t50 = sin(qJ(1));
t53 = cos(qJ(2));
t54 = cos(qJ(1));
t69 = cos(pkin(6));
t62 = t54 * t69;
t25 = t49 * t50 - t53 * t62;
t47 = sin(qJ(6));
t51 = cos(qJ(6));
t26 = t49 * t62 + t50 * t53;
t44 = qJ(3) + pkin(11);
t41 = sin(t44);
t42 = cos(t44);
t45 = sin(pkin(6));
t75 = t45 * t54;
t7 = t26 * t41 + t42 * t75;
t87 = t25 * t51 + t47 * t7;
t86 = -t25 * t47 + t51 * t7;
t85 = g(3) * t45;
t63 = t50 * t69;
t28 = -t49 * t63 + t53 * t54;
t48 = sin(qJ(3));
t82 = t28 * t48;
t81 = t41 * t47;
t80 = t41 * t51;
t79 = t45 * t49;
t78 = t45 * t50;
t52 = cos(qJ(3));
t77 = t45 * t52;
t76 = t45 * t53;
t46 = -qJ(4) - pkin(9);
t74 = t46 * t49;
t73 = t47 * t53;
t72 = t51 * t53;
t40 = pkin(3) * t52 + pkin(2);
t71 = -t25 * t40 - t26 * t46;
t27 = t49 * t54 + t53 * t63;
t70 = -t27 * t40 - t28 * t46;
t68 = t48 * t79;
t67 = t48 * t78;
t66 = t50 * t77;
t35 = t48 * t75;
t65 = t52 * t75;
t64 = t26 * t52 - t35;
t61 = t69 * t52;
t60 = pkin(1) * t54 + pkin(3) * t67 + pkin(8) * t78 - t27 * t46 + t28 * t40;
t6 = g(1) * t25 - g(2) * t27;
t59 = pkin(4) * t42 + qJ(5) * t41;
t10 = -t26 * t42 + t41 * t75;
t58 = t26 * t48 + t65;
t57 = -pkin(1) * t50 + pkin(3) * t35 + pkin(8) * t75 + t25 * t46 - t26 * t40;
t12 = t28 * t42 + t41 * t78;
t20 = t41 * t69 + t42 * t79;
t56 = -g(1) * t12 + g(2) * t10 - g(3) * t20;
t4 = -g(1) * t27 - g(2) * t25 + g(3) * t76;
t55 = g(1) * t28 + g(2) * t26 + g(3) * t79;
t39 = pkin(3) * t61;
t32 = pkin(3) * t66;
t29 = t40 * t76;
t19 = t41 * t79 - t42 * t69;
t14 = t28 * t52 + t67;
t13 = t66 - t82;
t11 = t28 * t41 - t42 * t78;
t3 = t11 * t47 + t27 * t51;
t2 = t11 * t51 - t27 * t47;
t1 = -g(1) * t11 - g(2) * t7 - g(3) * t19;
t5 = [0, g(1) * t50 - g(2) * t54, g(1) * t54 + g(2) * t50, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t28, -t6, 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t14, -g(1) * t58 - g(2) * t13, t6, -g(1) * t57 - g(2) * t60, t6, g(1) * t10 + g(2) * t12, g(1) * t7 - g(2) * t11, -g(1) * (pkin(4) * t10 - qJ(5) * t7 + t57) - g(2) * (pkin(4) * t12 + qJ(5) * t11 + t60) 0, 0, 0, 0, 0, g(1) * t87 - g(2) * t3, g(1) * t86 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, -t4, t55, 0, 0, 0, 0, 0, -t4 * t52, t4 * t48, -t55, -g(1) * t70 - g(2) * t71 - g(3) * (-t45 * t74 + t29) -t55, t4 * t42, -t4 * t41, -g(1) * (-t27 * t59 + t70) - g(2) * (-t25 * t59 + t71) - g(3) * t29 - (t53 * t59 - t74) * t85, 0, 0, 0, 0, 0, -g(1) * (-t27 * t81 + t28 * t51) - g(2) * (-t25 * t81 + t26 * t51) - (t41 * t73 + t49 * t51) * t85, -g(1) * (-t27 * t80 - t28 * t47) - g(2) * (-t25 * t80 - t26 * t47) - (t41 * t72 - t47 * t49) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t58 - g(3) * (t61 - t68) g(1) * t14 + g(2) * t64 - g(3) * (-t48 * t69 - t49 * t77) 0, -g(1) * t32 - g(3) * t39 + (g(2) * t65 + t48 * t55) * pkin(3), 0, t1, t56, -g(1) * (-pkin(3) * t82 - pkin(4) * t11 + qJ(5) * t12 + t32) - g(2) * (-pkin(3) * t58 - t7 * pkin(4) - qJ(5) * t10) - g(3) * (-pkin(3) * t68 - pkin(4) * t19 + qJ(5) * t20 + t39) 0, 0, 0, 0, 0, t56 * t47, t56 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t86 - g(3) * (t19 * t51 + t45 * t73) g(1) * t3 + g(2) * t87 - g(3) * (-t19 * t47 + t45 * t72);];
taug_reg  = t5;
