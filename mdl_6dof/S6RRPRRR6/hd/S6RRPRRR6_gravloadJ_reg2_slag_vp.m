% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t84 = qJ(4) + qJ(5);
t47 = sin(t84);
t61 = cos(qJ(1));
t56 = sin(qJ(2));
t78 = cos(t84);
t72 = t56 * t78;
t60 = cos(qJ(2));
t89 = t60 * t61;
t12 = t47 * t89 - t61 * t72;
t21 = t56 * t47 + t60 * t78;
t13 = t21 * t61;
t111 = -t12 * pkin(5) + t13 * pkin(10);
t22 = -t60 * t47 + t72;
t57 = sin(qJ(1));
t10 = t22 * t57;
t11 = t21 * t57;
t110 = t10 * pkin(5) + t11 * pkin(10);
t3 = g(1) * t12 - g(2) * t10 + g(3) * t21;
t54 = sin(qJ(6));
t1 = t3 * t54;
t58 = cos(qJ(6));
t2 = t3 * t58;
t55 = sin(qJ(4));
t82 = t55 * t89;
t34 = pkin(4) * t82;
t109 = -t34 + t111;
t105 = pkin(4) * t57;
t91 = t60 * t55;
t31 = t91 * t105;
t108 = -t31 + t110;
t35 = g(1) * t61 + g(2) * t57;
t107 = t35 * t56;
t79 = -t21 * pkin(5) + t22 * pkin(10);
t48 = t56 * qJ(3);
t87 = t60 * pkin(2) + t48;
t104 = g(1) * t57;
t101 = g(3) * t22;
t97 = t60 * pkin(3);
t59 = cos(qJ(4));
t44 = t59 * pkin(4) + pkin(3);
t96 = pkin(2) + t44;
t95 = t56 * t55;
t94 = t56 * t59;
t93 = t56 * t61;
t90 = t60 * t59;
t51 = t61 * pkin(7);
t62 = -pkin(9) - pkin(8);
t88 = t61 * t62 + t51;
t86 = t61 * pkin(1) + t57 * pkin(7);
t85 = qJ(3) * t60;
t42 = pkin(4) * t95;
t83 = t59 * t93;
t77 = t60 * t44 + t42 + t87;
t76 = pkin(2) * t89 + t61 * t48 + t86;
t75 = -pkin(4) * t90 - t42;
t74 = -g(1) * t10 - g(2) * t12;
t73 = -g(2) * t61 + t104;
t71 = t11 * t58 + t61 * t54;
t70 = t11 * t54 - t61 * t58;
t69 = t91 - t94;
t25 = t90 + t95;
t68 = -pkin(1) - t87;
t67 = t61 * t42 + t44 * t89 + t57 * t62 + t76;
t5 = g(1) * t13 + g(2) * t11 + t101;
t17 = t69 * t57;
t19 = t82 - t83;
t66 = g(1) * t19 + g(2) * t17 + g(3) * t25;
t18 = t25 * t57;
t20 = t25 * t61;
t65 = g(1) * t20 + g(2) * t18 - g(3) * t69;
t64 = t96 * t107;
t63 = (-pkin(1) - t96 * t60 + (-pkin(4) * t55 - qJ(3)) * t56) * t104;
t41 = t61 * t85;
t39 = t57 * t85;
t33 = pkin(4) * t83;
t30 = t94 * t105;
t24 = t73 * t60;
t23 = t73 * t56;
t15 = g(3) * t56 + t35 * t60;
t14 = -g(3) * t60 + t107;
t7 = t13 * t58 - t57 * t54;
t6 = -t13 * t54 - t57 * t58;
t4 = [0, 0, 0, 0, 0, 0, t73, t35, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, -t35, -g(1) * (-t57 * pkin(1) + t51) - g(2) * t86, 0, 0, 0, 0, 0, 0, t24, -t35, t23, -g(1) * t51 - g(2) * t76 - t68 * t104, 0, 0, 0, 0, 0, 0, g(1) * t18 - g(2) * t20, -g(1) * t17 + g(2) * t19, t35, -g(1) * (-t61 * pkin(8) + t51) - g(2) * (pkin(3) * t89 + t76) + (-g(1) * (t68 - t97) + g(2) * pkin(8)) * t57, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t13, -t74, t35, -g(1) * t88 - g(2) * t67 - t63, 0, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t7, -g(1) * t70 - g(2) * t6, t74, -g(1) * (-t11 * pkin(5) + t10 * pkin(10) + t88) - g(2) * (t13 * pkin(5) + t12 * pkin(10) + t67) - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t15, -g(1) * (-pkin(2) * t93 + t41) - g(2) * (-t57 * t56 * pkin(2) + t39) - g(3) * t87, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, -g(1) * t41 - g(2) * t39 - g(3) * (t87 + t97) + (pkin(2) + pkin(3)) * t107, 0, 0, 0, 0, 0, 0, -t3, -t5, 0, -g(1) * (t34 + t41) - g(2) * (t31 + t39) - g(3) * t77 + t64, 0, 0, 0, 0, 0, 0, -t2, t1, t5, -g(1) * (t41 - t109) - g(2) * (t39 - t108) - g(3) * (t77 - t79) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * (-t34 + t33) - g(2) * (-t31 + t30) - g(3) * t75, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t33 + t109) - g(2) * (t30 + t108) - g(3) * (t75 + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t111 - g(2) * t110 - g(3) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t70 + t54 * t101, g(1) * t7 + g(2) * t71 + t58 * t101, 0, 0;];
taug_reg  = t4;
