% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t119 = pkin(4) * t71 + qJ(5) * t68;
t66 = sin(pkin(10));
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t100 = cos(pkin(10));
t101 = cos(pkin(6));
t84 = t101 * t100;
t52 = t66 * t73 + t70 * t84;
t69 = sin(qJ(3));
t72 = cos(qJ(3));
t67 = sin(pkin(6));
t90 = t67 * t100;
t29 = -t52 * t69 - t72 * t90;
t118 = t119 * t29;
t109 = t67 * t72;
t91 = t66 * t101;
t54 = t100 * t73 - t70 * t91;
t31 = t109 * t66 - t54 * t69;
t117 = t119 * t31;
t110 = t67 * t70;
t55 = t101 * t72 - t110 * t69;
t116 = t119 * t55;
t115 = pkin(5) + pkin(9);
t114 = pkin(3) * t72;
t51 = t66 * t70 - t73 * t84;
t112 = t51 * t69;
t53 = t100 * t70 + t73 * t91;
t111 = t53 * t69;
t108 = t67 * t73;
t107 = t68 * t72;
t106 = t71 * t72;
t105 = t72 * t73;
t104 = pkin(2) * t108 + pkin(8) * t110;
t102 = qJ(6) * t71;
t99 = t69 * t108;
t98 = t68 * t108;
t97 = -t51 * pkin(2) + t52 * pkin(8);
t96 = -t53 * pkin(2) + t54 * pkin(8);
t26 = t29 * pkin(3);
t30 = t52 * t72 - t69 * t90;
t95 = t30 * pkin(9) + t26;
t27 = t31 * pkin(3);
t32 = t66 * t67 * t69 + t54 * t72;
t94 = t32 * pkin(9) + t27;
t50 = t55 * pkin(3);
t56 = t101 * t69 + t109 * t70;
t93 = t56 * pkin(9) + t50;
t11 = t30 * t68 - t51 * t71;
t12 = t30 * t71 + t51 * t68;
t92 = -t11 * pkin(4) + t12 * qJ(5);
t13 = t32 * t68 - t53 * t71;
t14 = t32 * t71 + t53 * t68;
t89 = -t13 * pkin(4) + t14 * qJ(5);
t33 = t108 * t71 + t56 * t68;
t34 = t56 * t71 - t98;
t88 = -t33 * pkin(4) + t34 * qJ(5);
t87 = t67 * pkin(3) * t105 + pkin(9) * t99 + t104;
t86 = -pkin(9) * t112 - t51 * t114 + t97;
t85 = -pkin(9) * t111 - t53 * t114 + t96;
t2 = g(1) * t13 + g(2) * t11 + g(3) * t33;
t83 = g(1) * t14 + g(2) * t12 + g(3) * t34;
t18 = -t107 * t51 - t52 * t71;
t20 = -t107 * t53 - t54 * t71;
t36 = -t110 * t71 + t72 * t98;
t82 = g(1) * t20 + g(2) * t18 + g(3) * t36;
t19 = -t106 * t51 + t52 * t68;
t21 = -t106 * t53 + t54 * t68;
t37 = (t105 * t71 + t68 * t70) * t67;
t81 = g(1) * t21 + g(2) * t19 + g(3) * t37;
t80 = g(1) * t31 + g(2) * t29 + g(3) * t55;
t79 = g(1) * t32 + g(2) * t30 + g(3) * t56;
t78 = t37 * pkin(4) + t36 * qJ(5) + t87;
t77 = -g(1) * t53 - g(2) * t51 + g(3) * t108;
t76 = g(1) * t54 + g(2) * t52 + g(3) * t110;
t75 = t19 * pkin(4) + t18 * qJ(5) + t86;
t74 = t21 * pkin(4) + t20 * qJ(5) + t85;
t15 = t77 * t69;
t7 = t80 * t71;
t6 = t80 * t68;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t76, 0, 0, 0, 0, 0, 0, 0, 0, -t77 * t72, t15, -t76, -g(1) * t96 - g(2) * t97 - g(3) * t104, 0, 0, 0, 0, 0, 0, -t81, t82, -t15, -g(1) * t85 - g(2) * t86 - g(3) * t87, 0, 0, 0, 0, 0, 0, -t15, t81, -t82, -g(1) * t74 - g(2) * t75 - g(3) * t78, 0, 0, 0, 0, 0, 0, -t15, -t82, -t81, -g(1) * (-pkin(5) * t111 + t21 * qJ(6) + t74) - g(2) * (-pkin(5) * t112 + t19 * qJ(6) + t75) - g(3) * (pkin(5) * t99 + t37 * qJ(6) + t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t79, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t79, -g(1) * t94 - g(2) * t95 - g(3) * t93, 0, 0, 0, 0, 0, 0, -t79, t7, -t6, -g(1) * (t94 + t117) - g(2) * (t95 + t118) - g(3) * (t93 + t116) 0, 0, 0, 0, 0, 0, -t79, -t6, -t7, -g(1) * (t31 * t102 + t115 * t32 + t117 + t27) - g(2) * (t29 * t102 + t115 * t30 + t118 + t26) - g(3) * (t55 * t102 + t115 * t56 + t116 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t83, -g(1) * t89 - g(2) * t92 - g(3) * t88, 0, 0, 0, 0, 0, 0, 0, -t83, t2, -g(1) * (-t13 * qJ(6) + t89) - g(2) * (-t11 * qJ(6) + t92) - g(3) * (-t33 * qJ(6) + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83;];
taug_reg  = t1;
