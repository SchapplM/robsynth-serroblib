% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR15_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t57 = cos(qJ(2));
t58 = cos(qJ(1));
t76 = cos(pkin(6));
t70 = t58 * t76;
t37 = t54 * t53 - t57 * t70;
t38 = t53 * t70 + t54 * t57;
t47 = sin(pkin(7));
t52 = sin(qJ(3));
t49 = cos(pkin(7));
t81 = t49 * t52;
t48 = sin(pkin(6));
t82 = t48 * t58;
t88 = cos(qJ(3));
t14 = -t47 * t52 * t82 - t37 * t81 + t38 * t88;
t50 = sin(qJ(6));
t55 = cos(qJ(6));
t68 = t47 * t48 * t88;
t74 = t49 * t88;
t13 = t37 * t74 + t38 * t52 + t58 * t68;
t26 = -t37 * t47 + t49 * t82;
t51 = sin(qJ(5));
t56 = cos(qJ(5));
t6 = t13 * t51 - t26 * t56;
t94 = -t14 * t55 + t6 * t50;
t93 = t14 * t50 + t6 * t55;
t90 = t13 * t56 + t26 * t51;
t89 = pkin(10) * t47;
t87 = t47 * t51;
t86 = t47 * t56;
t85 = t48 * t53;
t84 = t48 * t54;
t83 = t48 * t57;
t80 = t50 * t51;
t79 = t51 * t55;
t78 = t52 * t53;
t77 = t52 * t57;
t75 = t47 * t85;
t73 = t88 * t53;
t72 = t88 * t57;
t71 = t54 * t76;
t69 = t76 * t47;
t39 = -t58 * t53 - t57 * t71;
t40 = -t53 * t71 + t58 * t57;
t17 = -t39 * t74 + t40 * t52 - t54 * t68;
t67 = -g(1) * t13 + g(2) * t17;
t18 = t40 * t88 + (t39 * t49 + t47 * t84) * t52;
t66 = -g(1) * t14 + g(2) * t18;
t28 = -t39 * t47 + t49 * t84;
t65 = -g(1) * t26 - g(2) * t28;
t24 = -t88 * t69 + (-t49 * t72 + t78) * t48;
t36 = -t47 * t83 + t76 * t49;
t7 = t17 * t56 - t28 * t51;
t64 = g(1) * t7 + g(2) * t90 + g(3) * (t24 * t56 - t36 * t51);
t63 = g(1) * t17 + g(2) * t13 + g(3) * t24;
t25 = t52 * t69 + (t49 * t77 + t73) * t48;
t62 = g(1) * t18 + g(2) * t14 + g(3) * t25;
t19 = -t37 * t52 + t38 * t74;
t21 = t39 * t52 + t40 * t74;
t34 = (t49 * t73 + t77) * t48;
t61 = g(1) * t21 + g(2) * t19 + g(3) * t34;
t20 = -t37 * t88 - t38 * t81;
t22 = t39 * t88 - t40 * t81;
t35 = (-t49 * t78 + t72) * t48;
t60 = g(1) * t22 + g(2) * t20 + g(3) * t35;
t59 = g(1) * t40 + g(2) * t38 + g(3) * t85;
t23 = t34 * t51 + t56 * t75;
t12 = t24 * t51 + t36 * t56;
t10 = t21 * t51 + t40 * t86;
t9 = t19 * t51 + t38 * t86;
t8 = t17 * t51 + t28 * t56;
t2 = t18 * t50 + t8 * t55;
t1 = t18 * t55 - t8 * t50;
t3 = [0, g(1) * t54 - g(2) * t58, g(1) * t58 + g(2) * t54, 0, 0, 0, 0, 0, g(1) * t38 - g(2) * t40, -g(1) * t37 - g(2) * t39, 0, 0, 0, 0, 0, -t66, t67, t65, t66, -t67, -g(1) * (-t54 * pkin(1) - t38 * pkin(2) - pkin(3) * t14 + pkin(9) * t82 - qJ(4) * t13) - g(2) * (t58 * pkin(1) + t40 * pkin(2) + t18 * pkin(3) + pkin(9) * t84 + t17 * qJ(4)) + t65 * pkin(10), 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, g(1) * t90 - g(2) * t7, 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t2, -g(1) * t94 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t39 + g(2) * t37 - g(3) * t83, t59, 0, 0, 0, 0, 0, -t60, t61, -t59 * t47, t60, -t61, -g(1) * (t39 * pkin(2) + t22 * pkin(3) + t21 * qJ(4) + t40 * t89) - g(2) * (-t37 * pkin(2) + t20 * pkin(3) + t19 * qJ(4) + t38 * t89) - g(3) * (t35 * pkin(3) + t34 * qJ(4) + (pkin(2) * t57 + t53 * t89) * t48) 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t9 - g(3) * t23, -g(1) * (t21 * t56 - t40 * t87) - g(2) * (t19 * t56 - t38 * t87) - g(3) * (t34 * t56 - t51 * t75) 0, 0, 0, 0, 0, -g(1) * (t10 * t55 + t22 * t50) - g(2) * (t20 * t50 + t9 * t55) - g(3) * (t23 * t55 + t35 * t50) -g(1) * (-t10 * t50 + t22 * t55) - g(2) * (t20 * t55 - t9 * t50) - g(3) * (-t23 * t50 + t35 * t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t62, 0, -t63, -t62, -g(1) * (-t17 * pkin(3) + t18 * qJ(4)) - g(2) * (-t13 * pkin(3) + t14 * qJ(4)) - g(3) * (-t24 * pkin(3) + t25 * qJ(4)) 0, 0, 0, 0, 0, -t62 * t51, -t62 * t56, 0, 0, 0, 0, 0, -g(1) * (-t17 * t50 + t18 * t79) - g(2) * (-t13 * t50 + t14 * t79) - g(3) * (-t24 * t50 + t25 * t79) -g(1) * (-t17 * t55 - t18 * t80) - g(2) * (-t13 * t55 - t14 * t80) - g(3) * (-t24 * t55 - t25 * t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, g(1) * t8 + g(2) * t6 + g(3) * t12, 0, 0, 0, 0, 0, -t64 * t55, t64 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t94 - g(3) * (-t12 * t50 + t25 * t55) g(1) * t2 + g(2) * t93 - g(3) * (-t12 * t55 - t25 * t50);];
taug_reg  = t3;
