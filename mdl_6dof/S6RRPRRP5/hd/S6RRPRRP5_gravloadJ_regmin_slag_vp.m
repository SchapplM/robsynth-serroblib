% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:53:51
% EndTime: 2019-05-06 17:53:53
% DurationCPUTime: 0.57s
% Computational Cost: add. (446->98), mult. (1128->167), div. (0->0), fcn. (1437->12), ass. (0->68)
t40 = sin(pkin(11));
t46 = sin(qJ(2));
t50 = cos(qJ(2));
t69 = cos(pkin(11));
t30 = -t50 * t40 - t46 * t69;
t47 = sin(qJ(1));
t51 = cos(qJ(1));
t42 = cos(pkin(6));
t57 = -t46 * t40 + t50 * t69;
t53 = t57 * t42;
t11 = t47 * t30 + t51 * t53;
t44 = sin(qJ(5));
t48 = cos(qJ(5));
t23 = t30 * t42;
t12 = -t51 * t23 + t47 * t57;
t45 = sin(qJ(4));
t49 = cos(qJ(4));
t41 = sin(pkin(6));
t77 = t41 * t51;
t5 = t12 * t49 - t45 * t77;
t94 = t11 * t48 + t5 * t44;
t93 = -t11 * t44 + t5 * t48;
t22 = t30 * t41;
t17 = -t22 * t49 + t42 * t45;
t14 = t51 * t30 - t47 * t53;
t13 = -t47 * t23 - t51 * t57;
t79 = t41 * t47;
t9 = -t13 * t49 + t45 * t79;
t2 = -t14 * t48 - t9 * t44;
t21 = t57 * t41;
t92 = g(2) * t94 - g(3) * (-t17 * t44 - t21 * t48) - g(1) * t2;
t71 = t51 * t46;
t73 = t47 * t50;
t27 = -t42 * t73 - t71;
t78 = t41 * t50;
t91 = -g(1) * t27 - g(3) * t78;
t90 = -g(1) * t14 - g(2) * t11 - g(3) * t21;
t84 = t12 * t44;
t81 = t13 * t44;
t80 = t22 * t44;
t76 = t44 * t49;
t74 = t47 * t46;
t72 = t48 * t49;
t70 = t51 * t50;
t66 = t42 * t70;
t65 = pkin(5) * t44 + pkin(9);
t24 = t42 * t46 * pkin(2) + (-pkin(8) - qJ(3)) * t41;
t39 = t50 * pkin(2) + pkin(1);
t64 = -t47 * t24 + t51 * t39;
t4 = t12 * t45 + t49 * t77;
t8 = -t13 * t45 - t49 * t79;
t62 = -g(1) * t4 + g(2) * t8;
t61 = g(1) * t51 + g(2) * t47;
t60 = g(1) * t47 - g(2) * t51;
t59 = -t51 * t24 - t47 * t39;
t16 = -t22 * t45 - t42 * t49;
t56 = g(1) * t8 + g(2) * t4 + g(3) * t16;
t55 = g(1) * t9 + g(2) * t5 + g(3) * t17;
t43 = -qJ(6) - pkin(10);
t38 = t48 * pkin(5) + pkin(4);
t31 = pkin(2) * t66;
t28 = -t42 * t74 + t70;
t26 = -t42 * t71 - t73;
t25 = -t66 + t74;
t20 = -g(3) * t42 - t60 * t41;
t3 = -t14 * t44 + t9 * t48;
t1 = t90 * t45;
t6 = [0, t60, t61, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t28, -g(1) * t25 - g(2) * t27, -t61 * t41, -g(1) * t59 - g(2) * t64, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t62, 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t3, -g(1) * t94 - g(2) * t2, -t62, -g(1) * (-t12 * pkin(3) + t65 * t11 - t38 * t5 + t4 * t43 + t59) - g(2) * (-pkin(3) * t13 - t65 * t14 + t9 * t38 - t8 * t43 + t64); 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t25 + t91, g(3) * t41 * t46 + g(1) * t28 - g(2) * t26, 0, -g(2) * t31 + (g(2) * t74 + t91) * pkin(2), 0, 0, 0, 0, 0, t90 * t49, -t1, 0, 0, 0, 0, 0, -g(1) * (t14 * t72 - t81) - g(2) * (t11 * t72 + t84) - g(3) * (t21 * t72 - t80) -g(1) * (-t13 * t48 - t14 * t76) - g(2) * (-t11 * t76 + t12 * t48) - g(3) * (-t21 * t76 - t22 * t48) t1, -g(1) * (t27 * pkin(2) - pkin(5) * t81 - t13 * pkin(9)) - g(2) * (-pkin(2) * t74 + pkin(5) * t84 + pkin(9) * t12 + t31) - g(3) * (pkin(2) * t78 - pkin(5) * t80 - t22 * pkin(9)) + t90 * (t38 * t49 - t43 * t45 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, 0, 0, 0, 0, t56 * t48, -t56 * t44, -t55, -g(1) * (-t8 * t38 - t9 * t43) - g(2) * (-t4 * t38 - t5 * t43) - g(3) * (-t16 * t38 - t17 * t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, g(1) * t3 + g(2) * t93 - g(3) * (-t17 * t48 + t21 * t44) 0, t92 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56;];
taug_reg  = t6;
