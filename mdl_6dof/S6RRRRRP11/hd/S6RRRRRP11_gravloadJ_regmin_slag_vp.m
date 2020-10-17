% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:07:50
% EndTime: 2019-05-08 07:07:54
% DurationCPUTime: 0.99s
% Computational Cost: add. (837->146), mult. (2316->251), div. (0->0), fcn. (2995->14), ass. (0->82)
t100 = cos(qJ(3));
t58 = cos(qJ(1));
t101 = cos(qJ(2));
t93 = cos(pkin(6));
t82 = t93 * t101;
t98 = sin(qJ(2));
t99 = sin(qJ(1));
t41 = -t58 * t82 + t99 * t98;
t81 = t93 * t98;
t42 = t99 * t101 + t58 * t81;
t55 = sin(qJ(3));
t90 = sin(pkin(7));
t91 = sin(pkin(6));
t76 = t91 * t90;
t73 = t58 * t76;
t92 = cos(pkin(7));
t85 = t55 * t92;
t17 = t42 * t100 - t41 * t85 - t55 * t73;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t77 = t92 * t91;
t64 = t41 * t90 - t58 * t77;
t5 = t17 * t57 + t64 * t54;
t56 = cos(qJ(5));
t80 = t92 * t100;
t16 = t100 * t73 + t41 * t80 + t42 * t55;
t53 = sin(qJ(5));
t97 = t16 * t53;
t115 = t5 * t56 + t97;
t114 = -t16 * t56 + t5 * t53;
t4 = t17 * t54 - t64 * t57;
t63 = t58 * t98 + t99 * t82;
t110 = t63 * t92 - t99 * t76;
t109 = t101 * t77 + t90 * t93;
t78 = t91 * t98;
t32 = t100 * t78 + t109 * t55;
t61 = -t101 * t76 + t93 * t92;
t15 = t32 * t57 + t61 * t54;
t43 = t58 * t101 - t99 * t81;
t20 = t110 * t100 + t43 * t55;
t21 = t43 * t100 - t110 * t55;
t59 = -t63 * t90 - t99 * t77;
t9 = t21 * t57 - t59 * t54;
t2 = t20 * t56 - t9 * t53;
t31 = -t109 * t100 + t55 * t78;
t108 = -g(1) * t2 + g(2) * t114 - g(3) * (-t15 * t53 + t31 * t56);
t107 = g(1) * t21 + g(2) * t17 + g(3) * t32;
t69 = g(1) * t20 + g(2) * t16 + g(3) * t31;
t96 = t20 * t53;
t95 = t53 * t57;
t94 = t56 * t57;
t89 = pkin(5) * t53 + pkin(11);
t88 = t91 * pkin(9);
t87 = t90 * pkin(10);
t86 = t54 * t90;
t84 = t57 * t90;
t8 = t21 * t54 + t59 * t57;
t83 = -g(1) * t4 + g(2) * t8;
t79 = t101 * t91;
t14 = t32 * t54 - t61 * t57;
t72 = g(1) * t8 + g(2) * t4 + g(3) * t14;
t71 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t23 = -t41 * t100 - t42 * t85;
t10 = t23 * t54 - t42 * t84;
t25 = -t63 * t100 - t43 * t85;
t12 = t25 * t54 - t43 * t84;
t68 = t98 * t77;
t39 = t100 * t79 - t55 * t68;
t66 = t98 * t76;
t26 = t39 * t54 - t57 * t66;
t70 = g(1) * t12 + g(2) * t10 + g(3) * t26;
t52 = -qJ(6) - pkin(12);
t51 = t56 * pkin(5) + pkin(4);
t38 = t100 * t68 + t55 * t79;
t27 = t39 * t57 + t54 * t66;
t24 = t43 * t80 - t63 * t55;
t22 = -t41 * t55 + t42 * t80;
t13 = t25 * t57 + t43 * t86;
t11 = t23 * t57 + t42 * t86;
t3 = t9 * t56 + t96;
t1 = t69 * t54;
t6 = [0, g(1) * t99 - g(2) * t58, g(1) * t58 + g(2) * t99, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t43, -g(1) * t41 + g(2) * t63, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t21, -g(1) * t16 + g(2) * t20, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t83, 0, 0, 0, 0, 0, g(1) * t115 - g(2) * t3, -g(1) * t114 - g(2) * t2, -t83, -g(1) * (-t99 * pkin(1) - t42 * pkin(2) - pkin(3) * t17 - pkin(5) * t97 - pkin(11) * t16 + t4 * t52 - t5 * t51 + t58 * t88) - g(2) * (t58 * pkin(1) + t43 * pkin(2) + t21 * pkin(3) + pkin(5) * t96 + t20 * pkin(11) + t9 * t51 - t8 * t52 + t99 * t88) + (g(1) * t64 + g(2) * t59) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t63 + g(2) * t41 - g(3) * t79, g(1) * t43 + g(2) * t42 + g(3) * t78, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t23 - g(3) * t39, g(1) * t24 + g(2) * t22 + g(3) * t38, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t27, t70, 0, 0, 0, 0, 0, -g(1) * (t13 * t56 + t24 * t53) - g(2) * (t11 * t56 + t22 * t53) - g(3) * (t27 * t56 + t38 * t53) -g(1) * (-t13 * t53 + t24 * t56) - g(2) * (-t11 * t53 + t22 * t56) - g(3) * (-t27 * t53 + t38 * t56) -t70, -g(1) * (-t63 * pkin(2) + t25 * pkin(3) - t12 * t52 + t13 * t51 + t89 * t24 + t43 * t87) - g(2) * (-t41 * pkin(2) + t23 * pkin(3) - t10 * t52 + t11 * t51 + t89 * t22 + t42 * t87) - g(3) * (pkin(2) * t79 + t39 * pkin(3) + pkin(10) * t66 - t26 * t52 + t27 * t51 + t89 * t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t107, 0, 0, 0, 0, 0, t69 * t57, -t1, 0, 0, 0, 0, 0, -g(1) * (-t20 * t94 + t21 * t53) - g(2) * (-t16 * t94 + t17 * t53) - g(3) * (-t31 * t94 + t32 * t53) -g(1) * (t20 * t95 + t21 * t56) - g(2) * (t16 * t95 + t17 * t56) - g(3) * (t31 * t95 + t32 * t56) t1, -t107 * t89 + t69 * (t51 * t57 - t52 * t54 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, 0, 0, 0, 0, 0, t72 * t56, -t72 * t53, -t71, -g(1) * (-t8 * t51 - t9 * t52) - g(2) * (-t4 * t51 - t5 * t52) - g(3) * (-t14 * t51 - t15 * t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, g(1) * t3 + g(2) * t115 - g(3) * (-t15 * t56 - t31 * t53) 0, t108 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72;];
taug_reg  = t6;
