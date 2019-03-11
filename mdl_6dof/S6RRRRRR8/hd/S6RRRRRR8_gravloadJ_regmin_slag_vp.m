% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t59 = cos(qJ(1));
t74 = cos(pkin(6));
t68 = t59 * t74;
t84 = sin(qJ(1));
t39 = t84 * t55 - t58 * t68;
t40 = t55 * t68 + t84 * t58;
t54 = sin(qJ(3));
t73 = cos(pkin(7));
t85 = cos(qJ(3));
t64 = t73 * t85;
t50 = sin(pkin(7));
t51 = sin(pkin(6));
t79 = t50 * t51;
t71 = t59 * t79;
t18 = t39 * t64 + t40 * t54 + t85 * t71;
t52 = sin(qJ(6));
t56 = cos(qJ(6));
t69 = t54 * t73;
t19 = -t39 * t69 + t40 * t85 - t54 * t71;
t70 = t51 * t73;
t32 = t39 * t50 - t59 * t70;
t49 = qJ(4) + qJ(5);
t47 = sin(t49);
t48 = cos(t49);
t8 = t19 * t48 + t32 * t47;
t93 = -t18 * t56 + t8 * t52;
t92 = t18 * t52 + t8 * t56;
t89 = t19 * t47 - t32 * t48;
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t88 = t19 * t57 + t32 * t53;
t87 = t19 * t53 - t32 * t57;
t65 = t74 * t84;
t61 = t59 * t55 + t58 * t65;
t86 = t61 * t73 - t84 * t79;
t83 = t47 * t50;
t82 = t48 * t50;
t81 = t48 * t52;
t80 = t48 * t56;
t78 = t50 * t53;
t77 = t50 * t57;
t76 = t51 * t55;
t75 = t51 * t58;
t72 = t50 * t76;
t67 = t74 * t50;
t41 = -t55 * t65 + t59 * t58;
t23 = t41 * t85 - t86 * t54;
t33 = t61 * t50 + t84 * t70;
t10 = -t23 * t47 + t33 * t48;
t30 = t54 * t67 + (t85 * t55 + t58 * t69) * t51;
t38 = -t50 * t75 + t74 * t73;
t63 = g(1) * t10 - g(2) * t89 + g(3) * (-t30 * t47 + t38 * t48);
t22 = t41 * t54 + t86 * t85;
t29 = t54 * t76 - t64 * t75 - t85 * t67;
t62 = g(1) * t22 + g(2) * t18 + g(3) * t29;
t37 = (-t55 * t69 + t85 * t58) * t51;
t36 = (t54 * t58 + t55 * t64) * t51;
t28 = t37 * t48 + t47 * t72;
t27 = -t41 * t69 - t61 * t85;
t26 = t41 * t64 - t61 * t54;
t25 = -t39 * t85 - t40 * t69;
t24 = -t39 * t54 + t40 * t64;
t17 = t30 * t48 + t38 * t47;
t15 = t27 * t48 + t41 * t83;
t14 = t25 * t48 + t40 * t83;
t13 = t23 * t57 + t33 * t53;
t12 = -t23 * t53 + t33 * t57;
t11 = t23 * t48 + t33 * t47;
t6 = t11 * t56 + t22 * t52;
t5 = -t11 * t52 + t22 * t56;
t4 = g(1) * t11 + g(2) * t8 + g(3) * t17;
t2 = t63 * t56;
t1 = t63 * t52;
t3 = [0, g(1) * t84 - g(2) * t59, g(1) * t59 + g(2) * t84, 0, 0, 0, 0, 0, g(1) * t40 - g(2) * t41, -g(1) * t39 + g(2) * t61, 0, 0, 0, 0, 0, g(1) * t19 - g(2) * t23, -g(1) * t18 + g(2) * t22, 0, 0, 0, 0, 0, g(1) * t88 - g(2) * t13, -g(1) * t87 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t89 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t92 - g(2) * t6, -g(1) * t93 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t61 + g(2) * t39 - g(3) * t75, g(1) * t41 + g(2) * t40 + g(3) * t76, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t25 - g(3) * t37, g(1) * t26 + g(2) * t24 + g(3) * t36, 0, 0, 0, 0, 0, -g(1) * (t27 * t57 + t41 * t78) - g(2) * (t25 * t57 + t40 * t78) - g(3) * (t37 * t57 + t53 * t72) -g(1) * (-t27 * t53 + t41 * t77) - g(2) * (-t25 * t53 + t40 * t77) - g(3) * (-t37 * t53 + t57 * t72) 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t14 - g(3) * t28, -g(1) * (-t27 * t47 + t41 * t82) - g(2) * (-t25 * t47 + t40 * t82) - g(3) * (-t37 * t47 + t48 * t72) 0, 0, 0, 0, 0, -g(1) * (t15 * t56 + t26 * t52) - g(2) * (t14 * t56 + t24 * t52) - g(3) * (t28 * t56 + t36 * t52) -g(1) * (-t15 * t52 + t26 * t56) - g(2) * (-t14 * t52 + t24 * t56) - g(3) * (-t28 * t52 + t36 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, g(1) * t23 + g(2) * t19 + g(3) * t30, 0, 0, 0, 0, 0, t62 * t57, -t62 * t53, 0, 0, 0, 0, 0, t62 * t48, -t62 * t47, 0, 0, 0, 0, 0, -g(1) * (-t22 * t80 + t23 * t52) - g(2) * (-t18 * t80 + t19 * t52) - g(3) * (-t29 * t80 + t30 * t52) -g(1) * (t22 * t81 + t23 * t56) - g(2) * (t18 * t81 + t19 * t56) - g(3) * (t29 * t81 + t30 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t87 - g(3) * (-t30 * t53 + t38 * t57) g(1) * t13 + g(2) * t88 - g(3) * (-t30 * t57 - t38 * t53) 0, 0, 0, 0, 0, -t63, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t93 - g(3) * (-t17 * t52 + t29 * t56) g(1) * t6 + g(2) * t92 - g(3) * (-t17 * t56 - t29 * t52);];
taug_reg  = t3;
