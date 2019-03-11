% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t22 = g(1) * t46 + g(2) * t43;
t39 = qJ(2) + pkin(9);
t35 = sin(t39);
t84 = t22 * t35;
t36 = cos(t39);
t83 = -t36 * pkin(3) - t35 * pkin(8);
t82 = -pkin(4) - pkin(5);
t42 = sin(qJ(2));
t81 = pkin(2) * t42;
t80 = g(1) * t43;
t40 = -qJ(3) - pkin(7);
t78 = g(2) * t40;
t41 = sin(qJ(4));
t76 = t35 * t41;
t44 = cos(qJ(4));
t75 = t35 * t44;
t74 = t35 * t46;
t73 = t36 * t44;
t72 = t36 * t46;
t71 = t43 * t41;
t70 = t43 * t44;
t69 = t46 * t40;
t68 = t46 * t41;
t67 = t46 * t44;
t66 = qJ(5) * t41;
t65 = qJ(6) * t36;
t64 = t35 * qJ(6);
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t34 = t37 + pkin(1);
t29 = t46 * t34;
t63 = pkin(3) * t72 + pkin(8) * t74 + t29;
t62 = t37 - t83;
t61 = -pkin(3) - t66;
t15 = t36 * t71 + t67;
t16 = t36 * t70 - t68;
t60 = -t15 * pkin(4) + t16 * qJ(5);
t17 = t36 * t68 - t70;
t18 = t36 * t67 + t71;
t59 = -t17 * pkin(4) + t18 * qJ(5);
t23 = t43 * t36 * pkin(8);
t58 = -t43 * t81 + t23;
t26 = pkin(8) * t72;
t57 = -t46 * t81 + t26;
t56 = pkin(4) * t73 + t36 * t66 + t62;
t55 = -pkin(3) * t35 - t81;
t4 = g(1) * t15 - g(2) * t17;
t21 = -g(2) * t46 + t80;
t54 = -t34 + t83;
t52 = -t16 * pkin(4) - t15 * qJ(5) - t69;
t51 = t18 * pkin(4) + t17 * qJ(5) + t63;
t2 = g(1) * t17 + g(2) * t15 + g(3) * t76;
t49 = g(1) * t18 + g(2) * t16 + g(3) * t75;
t8 = -g(3) * t36 + t84;
t48 = -g(3) * t45 + t22 * t42;
t47 = (-g(1) * t54 + t78) * t43;
t19 = qJ(5) * t75;
t14 = -g(2) * t74 + t35 * t80;
t9 = g(3) * t35 + t22 * t36;
t7 = t8 * t44;
t6 = t8 * t41;
t5 = g(1) * t16 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t45, -t21 * t42, -t22, -g(1) * (-t43 * pkin(1) + t46 * pkin(7)) - g(2) * (t46 * pkin(1) + t43 * pkin(7)) 0, 0, 0, 0, 0, 0, t21 * t36, -t14, -t22, -g(1) * (-t43 * t34 - t69) - g(2) * (-t43 * t40 + t29) 0, 0, 0, 0, 0, 0, t5, -t4, t14, g(1) * t69 - g(2) * t63 + t47, 0, 0, 0, 0, 0, 0, t5, t14, t4, -g(1) * t52 - g(2) * t51 + t47, 0, 0, 0, 0, 0, 0, t5, t4, -t14, -g(1) * (-t16 * pkin(5) + t52) - g(2) * (t18 * pkin(5) - t46 * t64 + t51) + (-g(1) * (t54 + t64) + t78) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(3) * t42 + t22 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t48 * pkin(2), 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (t55 * t46 + t26) - g(2) * (t55 * t43 + t23) - g(3) * t62, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t57 - g(2) * t58 - g(3) * t56 + (pkin(4) * t44 - t61) * t84, 0, 0, 0, 0, 0, 0, t7, t6, t9, -g(1) * (-t46 * t65 + t57) - g(2) * (-t43 * t65 + t58) - g(3) * (pkin(5) * t73 + t56) + (g(3) * qJ(6) + t22 * (-t82 * t44 - t61)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t49, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t49, -g(1) * t59 - g(2) * t60 - g(3) * (-pkin(4) * t76 + t19) 0, 0, 0, 0, 0, 0, t2, -t49, 0, -g(1) * (-t17 * pkin(5) + t59) - g(2) * (-t15 * pkin(5) + t60) - g(3) * (t82 * t76 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t1;
