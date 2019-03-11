% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t58 = sin(qJ(2));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t82 = cos(pkin(6));
t73 = t62 * t82;
t93 = sin(qJ(1));
t38 = t93 * t58 - t61 * t73;
t39 = t58 * t73 + t93 * t61;
t57 = sin(qJ(3));
t81 = cos(pkin(7));
t94 = cos(qJ(3));
t68 = t81 * t94;
t52 = sin(pkin(7));
t53 = sin(pkin(6));
t83 = t53 * t62;
t78 = t52 * t83;
t14 = t38 * t68 + t39 * t57 + t94 * t78;
t74 = t57 * t81;
t15 = -t38 * t74 + t39 * t94 - t57 * t78;
t75 = t53 * t81;
t28 = t38 * t52 - t62 * t75;
t51 = qJ(4) + pkin(13);
t49 = sin(t51);
t50 = cos(t51);
t4 = t15 * t50 + t28 * t49;
t55 = sin(qJ(6));
t59 = cos(qJ(6));
t106 = -t14 * t59 + t4 * t55;
t105 = t14 * t55 + t4 * t59;
t60 = cos(qJ(4));
t56 = sin(qJ(4));
t92 = t28 * t56;
t104 = t15 * t60 + t92;
t100 = t15 * t56 - t28 * t60;
t72 = t82 * t52;
t26 = t57 * t72 + (t94 * t58 + t61 * t74) * t53;
t84 = t53 * t61;
t37 = -t52 * t84 + t82 * t81;
t69 = t82 * t93;
t40 = -t62 * t58 - t61 * t69;
t41 = -t58 * t69 + t62 * t61;
t76 = t53 * t93;
t71 = t52 * t76;
t19 = t41 * t94 + (t81 * t40 + t71) * t57;
t29 = -t40 * t52 + t93 * t75;
t8 = -t19 * t56 + t29 * t60;
t99 = -g(1) * t8 + g(2) * t100 - g(3) * (-t26 * t56 + t37 * t60);
t85 = t53 * t58;
t98 = g(1) * t41 + g(2) * t39 + g(3) * t85;
t91 = t29 * t56;
t90 = t49 * t52;
t89 = t50 * t55;
t88 = t50 * t59;
t87 = t52 * t56;
t86 = t52 * t60;
t79 = t52 * t85;
t18 = -t40 * t68 + t41 * t57 - t94 * t71;
t70 = -g(1) * t14 + g(2) * t18;
t67 = g(1) * (-t19 * t49 + t29 * t50) + g(2) * (-t15 * t49 + t28 * t50) + g(3) * (-t26 * t49 + t37 * t50);
t25 = t57 * t85 - t68 * t84 - t94 * t72;
t66 = g(1) * t18 + g(2) * t14 + g(3) * t25;
t65 = g(1) * t19 + g(2) * t15 + g(3) * t26;
t21 = -t38 * t57 + t39 * t68;
t23 = t40 * t57 + t41 * t68;
t35 = (t57 * t61 + t58 * t68) * t53;
t64 = g(1) * t23 + g(2) * t21 + g(3) * t35;
t54 = -qJ(5) - pkin(11);
t48 = t60 * pkin(4) + pkin(3);
t36 = (-t58 * t74 + t94 * t61) * t53;
t24 = t40 * t94 - t41 * t74;
t22 = -t38 * t94 - t39 * t74;
t20 = t36 * t50 + t49 * t79;
t13 = t26 * t50 + t37 * t49;
t11 = t24 * t50 + t41 * t90;
t10 = t22 * t50 + t39 * t90;
t9 = t19 * t60 + t91;
t7 = t19 * t50 + t29 * t49;
t2 = t18 * t55 + t7 * t59;
t1 = t18 * t59 - t7 * t55;
t3 = [0, g(1) * t93 - g(2) * t62, g(1) * t62 + g(2) * t93, 0, 0, 0, 0, 0, g(1) * t39 - g(2) * t41, -g(1) * t38 - g(2) * t40, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t19, t70, 0, 0, 0, 0, 0, g(1) * t104 - g(2) * t9, -g(1) * t100 - g(2) * t8, -t70, -g(1) * (-t93 * pkin(1) - t39 * pkin(2) - pkin(4) * t92 + pkin(9) * t83 + t14 * t54 - t15 * t48) - g(2) * (t62 * pkin(1) + t41 * pkin(2) + pkin(4) * t91 + pkin(9) * t76 - t18 * t54 + t19 * t48) + (g(1) * t28 - g(2) * t29) * pkin(10), 0, 0, 0, 0, 0, g(1) * t105 - g(2) * t2, -g(1) * t106 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t40 + g(2) * t38 - g(3) * t84, t98, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t22 - g(3) * t36, t64, 0, 0, 0, 0, 0, -g(1) * (t24 * t60 + t41 * t87) - g(2) * (t22 * t60 + t39 * t87) - g(3) * (t36 * t60 + t56 * t79) -g(1) * (-t24 * t56 + t41 * t86) - g(2) * (-t22 * t56 + t39 * t86) - g(3) * (-t36 * t56 + t60 * t79) -t64, -g(1) * (t40 * pkin(2) - t23 * t54 + t24 * t48) - g(2) * (-t38 * pkin(2) - t21 * t54 + t22 * t48) - g(3) * (pkin(2) * t84 - t35 * t54 + t36 * t48) - t98 * t52 * (pkin(4) * t56 + pkin(10)) 0, 0, 0, 0, 0, -g(1) * (t11 * t59 + t23 * t55) - g(2) * (t10 * t59 + t21 * t55) - g(3) * (t20 * t59 + t35 * t55) -g(1) * (-t11 * t55 + t23 * t59) - g(2) * (-t10 * t55 + t21 * t59) - g(3) * (-t20 * t55 + t35 * t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, t66 * t60, -t66 * t56, -t65, -g(1) * (-t18 * t48 - t19 * t54) - g(2) * (-t14 * t48 - t15 * t54) - g(3) * (-t25 * t48 - t26 * t54) 0, 0, 0, 0, 0, -g(1) * (-t18 * t88 + t19 * t55) - g(2) * (-t14 * t88 + t15 * t55) - g(3) * (-t25 * t88 + t26 * t55) -g(1) * (t18 * t89 + t19 * t59) - g(2) * (t14 * t89 + t15 * t59) - g(3) * (t25 * t89 + t26 * t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, g(1) * t9 + g(2) * t104 - g(3) * (-t26 * t60 - t37 * t56) 0, t99 * pkin(4), 0, 0, 0, 0, 0, -t67 * t59, t67 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t106 - g(3) * (-t13 * t55 + t25 * t59) g(1) * t2 + g(2) * t105 - g(3) * (-t13 * t59 - t25 * t55);];
taug_reg  = t3;
