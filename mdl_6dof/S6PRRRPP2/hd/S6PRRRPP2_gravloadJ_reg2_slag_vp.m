% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:46:16
% EndTime: 2019-05-05 06:46:18
% DurationCPUTime: 0.54s
% Computational Cost: add. (549->127), mult. (1447->168), div. (0->0), fcn. (1793->10), ass. (0->83)
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t120 = -pkin(3) * t74 - pkin(9) * t71;
t70 = sin(qJ(4));
t103 = qJ(5) * t70;
t68 = sin(pkin(10));
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t100 = cos(pkin(10));
t101 = cos(pkin(6));
t84 = t101 * t100;
t54 = t68 * t75 + t72 * t84;
t69 = sin(pkin(6));
t91 = t69 * t100;
t31 = -t54 * t71 - t74 * t91;
t73 = cos(qJ(4));
t114 = t31 * t73;
t119 = pkin(4) * t114 + t31 * t103;
t110 = t69 * t74;
t92 = t68 * t101;
t56 = t100 * t75 - t72 * t92;
t33 = t110 * t68 - t56 * t71;
t113 = t33 * t73;
t118 = pkin(4) * t113 + t33 * t103;
t111 = t69 * t72;
t57 = t101 * t74 - t111 * t71;
t112 = t57 * t73;
t117 = pkin(4) * t112 + t57 * t103;
t109 = t69 * t75;
t108 = t70 * t74;
t107 = t73 * t74;
t106 = t74 * t75;
t105 = pkin(9) - qJ(6);
t104 = pkin(2) * t109 + pkin(8) * t111;
t102 = qJ(6) * t71;
t99 = t71 * t109;
t98 = t70 * t109;
t53 = t68 * t72 - t75 * t84;
t97 = -t53 * pkin(2) + t54 * pkin(8);
t55 = t100 * t72 + t75 * t92;
t96 = -t55 * pkin(2) + t56 * pkin(8);
t28 = t31 * pkin(3);
t32 = t54 * t74 - t71 * t91;
t95 = t32 * pkin(9) + t28;
t29 = t33 * pkin(3);
t34 = t68 * t69 * t71 + t56 * t74;
t94 = t34 * pkin(9) + t29;
t52 = t57 * pkin(3);
t58 = t101 * t71 + t110 * t72;
t93 = t58 * pkin(9) + t52;
t13 = t32 * t70 - t53 * t73;
t14 = t32 * t73 + t53 * t70;
t90 = -t13 * pkin(4) + t14 * qJ(5);
t15 = t34 * t70 - t55 * t73;
t16 = t34 * t73 + t55 * t70;
t89 = -t15 * pkin(4) + t16 * qJ(5);
t35 = t109 * t73 + t58 * t70;
t36 = t58 * t73 - t98;
t88 = -t35 * pkin(4) + t36 * qJ(5);
t87 = t69 * pkin(3) * t106 + pkin(9) * t99 + t104;
t86 = t120 * t53 + t97;
t85 = t120 * t55 + t96;
t2 = g(1) * t15 + g(2) * t13 + g(3) * t35;
t83 = g(1) * t16 + g(2) * t14 + g(3) * t36;
t20 = -t108 * t53 - t54 * t73;
t22 = -t108 * t55 - t56 * t73;
t38 = -t111 * t73 + t74 * t98;
t82 = g(1) * t22 + g(2) * t20 + g(3) * t38;
t81 = g(1) * t33 + g(2) * t31 + g(3) * t57;
t10 = g(1) * t34 + g(2) * t32 + g(3) * t58;
t39 = (t106 * t73 + t70 * t72) * t69;
t80 = t39 * pkin(4) + t38 * qJ(5) + t87;
t79 = -g(1) * t55 - g(2) * t53 + g(3) * t109;
t78 = g(1) * t56 + g(2) * t54 + g(3) * t111;
t21 = -t107 * t53 + t54 * t70;
t77 = t21 * pkin(4) + t20 * qJ(5) + t86;
t23 = -t107 * t55 + t56 * t70;
t76 = t23 * pkin(4) + t22 * qJ(5) + t85;
t17 = t79 * t71;
t7 = t81 * t73;
t6 = t81 * t70;
t5 = -g(1) * t23 - g(2) * t21 - g(3) * t39;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t78, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t74, t17, -t78, -g(1) * t96 - g(2) * t97 - g(3) * t104, 0, 0, 0, 0, 0, 0, t5, t82, -t17, -g(1) * t85 - g(2) * t86 - g(3) * t87, 0, 0, 0, 0, 0, 0, t5, -t17, -t82, -g(1) * t76 - g(2) * t77 - g(3) * t80, 0, 0, 0, 0, 0, 0, t5, -t82, t17, -g(1) * (t23 * pkin(5) + t102 * t55 + t76) - g(2) * (t21 * pkin(5) + t102 * t53 + t77) - g(3) * (t39 * pkin(5) - qJ(6) * t99 + t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t10, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t10, -g(1) * t94 - g(2) * t95 - g(3) * t93, 0, 0, 0, 0, 0, 0, -t7, -t10, -t6, -g(1) * (t94 + t118) - g(2) * (t95 + t119) - g(3) * (t93 + t117) 0, 0, 0, 0, 0, 0, -t7, -t6, t10, -g(1) * (pkin(5) * t113 + t105 * t34 + t118 + t29) - g(2) * (pkin(5) * t114 + t105 * t32 + t119 + t28) - g(3) * (pkin(5) * t112 + t105 * t58 + t117 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t83, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t83, -g(1) * t89 - g(2) * t90 - g(3) * t88, 0, 0, 0, 0, 0, 0, t2, -t83, 0, -g(1) * (-t15 * pkin(5) + t89) - g(2) * (-t13 * pkin(5) + t90) - g(3) * (-t35 * pkin(5) + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81;];
taug_reg  = t1;
