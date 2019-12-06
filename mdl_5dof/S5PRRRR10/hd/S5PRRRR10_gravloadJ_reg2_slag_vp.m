% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t51 = sin(pkin(6));
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t53 = cos(pkin(11));
t92 = cos(pkin(5));
t83 = t53 * t92;
t90 = sin(pkin(11));
t63 = t90 * t57 - t60 * t83;
t91 = cos(pkin(6));
t52 = sin(pkin(5));
t98 = t52 * t53;
t104 = t51 * t98 + t63 * t91;
t74 = t92 * t90;
t64 = t53 * t57 + t60 * t74;
t84 = t52 * t90;
t103 = -t51 * t84 + t64 * t91;
t102 = pkin(8) * t51;
t101 = cos(qJ(3));
t55 = sin(qJ(4));
t100 = t51 * t55;
t59 = cos(qJ(4));
t99 = t51 * t59;
t97 = t52 * t57;
t96 = t52 * t60;
t54 = sin(qJ(5));
t95 = t54 * t59;
t58 = cos(qJ(5));
t94 = t58 * t59;
t88 = t51 * t97;
t93 = pkin(2) * t96 + pkin(8) * t88;
t41 = t57 * t83 + t90 * t60;
t56 = sin(qJ(3));
t12 = t104 * t101 + t41 * t56;
t13 = t41 * t101 - t104 * t56;
t87 = -t12 * pkin(3) + t13 * pkin(9);
t42 = t53 * t60 - t57 * t74;
t14 = t103 * t101 + t42 * t56;
t15 = t42 * t101 - t103 * t56;
t86 = -t14 * pkin(3) + t15 * pkin(9);
t76 = t91 * t101;
t81 = t92 * t51;
t27 = -t101 * t81 + t56 * t97 - t76 * t96;
t82 = t56 * t91;
t28 = t56 * t81 + (t101 * t57 + t60 * t82) * t52;
t85 = -t27 * pkin(3) + t28 * pkin(9);
t80 = -t63 * pkin(2) + t41 * t102;
t79 = -t64 * pkin(2) + t42 * t102;
t77 = -pkin(4) * t59 - pkin(10) * t55;
t36 = (t56 * t60 + t57 * t76) * t52;
t37 = (t101 * t60 - t57 * t82) * t52;
t75 = t37 * pkin(3) + t36 * pkin(9) + t93;
t40 = -t51 * t96 + t92 * t91;
t16 = -t28 * t55 + t40 * t59;
t29 = t63 * t51 - t91 * t98;
t2 = -t13 * t55 + t29 * t59;
t30 = t64 * t51 + t91 * t84;
t4 = -t15 * t55 + t30 * t59;
t73 = g(1) * t4 + g(2) * t2 + g(3) * t16;
t17 = t28 * t59 + t40 * t55;
t3 = t13 * t59 + t29 * t55;
t5 = t15 * t59 + t30 * t55;
t72 = g(1) * t5 + g(2) * t3 + g(3) * t17;
t24 = t37 * t55 - t59 * t88;
t21 = -t63 * t101 - t41 * t82;
t6 = t21 * t55 - t41 * t99;
t23 = -t64 * t101 - t42 * t82;
t8 = t23 * t55 - t42 * t99;
t71 = g(1) * t8 + g(2) * t6 + g(3) * t24;
t70 = g(1) * t14 + g(2) * t12 + g(3) * t27;
t69 = g(1) * t15 + g(2) * t13 + g(3) * t28;
t20 = t41 * t76 - t63 * t56;
t22 = t42 * t76 - t64 * t56;
t68 = g(1) * t22 + g(2) * t20 + g(3) * t36;
t67 = t21 * pkin(3) + t20 * pkin(9) + t80;
t66 = t23 * pkin(3) + t22 * pkin(9) + t79;
t65 = g(1) * t42 + g(2) * t41 + g(3) * t97;
t25 = t37 * t59 + t55 * t88;
t9 = t42 * t100 + t23 * t59;
t7 = t41 * t100 + t21 * t59;
t1 = t70 * t55;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t64 + g(2) * t63 - g(3) * t96, t65, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t21 - g(3) * t37, t68, -t65 * t51, -g(1) * t79 - g(2) * t80 - g(3) * t93, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t25, t71, -t68, -g(1) * t66 - g(2) * t67 - g(3) * t75, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t54 + t9 * t58) - g(2) * (t20 * t54 + t7 * t58) - g(3) * (t25 * t58 + t36 * t54), -g(1) * (t22 * t58 - t9 * t54) - g(2) * (t20 * t58 - t7 * t54) - g(3) * (-t25 * t54 + t36 * t58), -t71, -g(1) * (t9 * pkin(4) + t8 * pkin(10) + t66) - g(2) * (t7 * pkin(4) + t6 * pkin(10) + t67) - g(3) * (t25 * pkin(4) + t24 * pkin(10) + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t59, -t1, -t69, -g(1) * t86 - g(2) * t87 - g(3) * t85, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t94 + t15 * t54) - g(2) * (-t12 * t94 + t13 * t54) - g(3) * (-t27 * t94 + t28 * t54), -g(1) * (t14 * t95 + t15 * t58) - g(2) * (t12 * t95 + t13 * t58) - g(3) * (t27 * t95 + t28 * t58), t1, -g(1) * (t77 * t14 + t86) - g(2) * (t77 * t12 + t87) - g(3) * (t77 * t27 + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t72, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t58, t73 * t54, -t72, -g(1) * (t4 * pkin(4) + t5 * pkin(10)) - g(2) * (t2 * pkin(4) + t3 * pkin(10)) - g(3) * (t16 * pkin(4) + t17 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t58 - t5 * t54) - g(2) * (t12 * t58 - t3 * t54) - g(3) * (-t17 * t54 + t27 * t58), -g(1) * (-t14 * t54 - t5 * t58) - g(2) * (-t12 * t54 - t3 * t58) - g(3) * (-t17 * t58 - t27 * t54), 0, 0;];
taug_reg = t10;
