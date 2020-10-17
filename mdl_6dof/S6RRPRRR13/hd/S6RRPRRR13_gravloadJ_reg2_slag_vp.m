% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 01:18:10
% EndTime: 2019-05-07 01:18:12
% DurationCPUTime: 0.78s
% Computational Cost: add. (568->149), mult. (1308->229), div. (0->0), fcn. (1588->12), ass. (0->86)
t51 = sin(qJ(2));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t87 = cos(pkin(6));
t80 = t56 * t87;
t29 = t51 * t80 + t52 * t55;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t28 = t51 * t52 - t55 * t80;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t48 = sin(pkin(6));
t95 = t48 * t56;
t69 = -t28 * t50 + t54 * t95;
t116 = t29 * t53 + t49 * t69;
t96 = t48 * t55;
t27 = -t50 * t96 + t54 * t87;
t81 = t52 * t87;
t30 = t51 * t56 + t55 * t81;
t97 = t48 * t52;
t14 = t30 * t50 + t54 * t97;
t31 = -t51 * t81 + t55 * t56;
t7 = -t14 * t49 + t31 * t53;
t90 = t51 * t53;
t118 = -g(2) * t116 - g(3) * (-t27 * t49 + t48 * t90) - g(1) * t7;
t117 = -g(1) * t31 - g(2) * t29;
t115 = -t29 * t49 + t53 * t69;
t47 = qJ(5) + qJ(6);
t44 = sin(t47);
t45 = cos(t47);
t114 = t29 * t45 + t44 * t69;
t113 = -t29 * t44 + t45 * t69;
t109 = g(3) * t48;
t108 = t28 * pkin(9);
t107 = t30 * pkin(9);
t106 = t28 * t49;
t101 = t30 * t49;
t100 = t44 * t50;
t99 = t45 * t50;
t98 = t48 * t51;
t94 = t49 * t50;
t93 = t49 * t55;
t92 = t50 * t51;
t91 = t50 * t53;
t89 = pkin(2) * t96 + qJ(3) * t98;
t88 = pkin(1) * t56 + pkin(8) * t97;
t86 = pkin(9) * t96 + t89;
t85 = pkin(5) * t49 + pkin(9);
t84 = -t52 * pkin(1) + pkin(8) * t95;
t22 = t28 * pkin(2);
t83 = -t22 - t108;
t24 = t30 * pkin(2);
t82 = -t24 - t107;
t79 = t29 * qJ(3) - t22;
t78 = t31 * qJ(3) - t24;
t77 = g(3) * t86;
t76 = pkin(4) * t50 - pkin(10) * t54;
t13 = -t30 * t54 + t50 * t97;
t68 = t28 * t54 + t50 * t95;
t75 = g(1) * t68 + g(2) * t13;
t74 = g(1) * t28 - g(2) * t30;
t12 = g(1) * t29 - g(2) * t31;
t73 = g(1) * t56 + g(2) * t52;
t43 = pkin(5) * t53 + pkin(4);
t57 = -pkin(11) - pkin(10);
t72 = t43 * t50 + t54 * t57;
t71 = pkin(2) * t31 + t30 * qJ(3) + t88;
t66 = pkin(3) * t97 + t71;
t65 = -pkin(2) * t29 - t28 * qJ(3) + t84;
t26 = -t50 * t87 - t54 * t96;
t64 = g(1) * t13 - g(2) * t68 - g(3) * t26;
t63 = g(1) * t14 - g(2) * t69 + g(3) * t27;
t62 = pkin(3) * t95 + t65;
t10 = -g(1) * t30 - g(2) * t28 + g(3) * t96;
t61 = g(3) * t98 - t117;
t60 = pkin(9) * t31 + t66;
t59 = -pkin(9) * t29 + t62;
t32 = t73 * t48;
t9 = t61 * t54;
t8 = t14 * t53 + t31 * t49;
t6 = t14 * t45 + t31 * t44;
t5 = -t14 * t44 + t31 * t45;
t2 = g(1) * t6 - g(2) * t113 - g(3) * (-t27 * t45 - t44 * t98);
t1 = -g(1) * t5 - g(2) * t114 - g(3) * (-t27 * t44 + t45 * t98);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t56, t73, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t74, -t32, -g(1) * t84 - g(2) * t88, 0, 0, 0, 0, 0, 0, -t32, -t12, t74, -g(1) * t65 - g(2) * t71, 0, 0, 0, 0, 0, 0, -g(1) * t69 - g(2) * t14, t75, t12, -g(1) * t59 - g(2) * t60, 0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t8, g(1) * t116 - g(2) * t7, -t75, -g(1) * (pkin(4) * t69 + pkin(10) * t68 + t59) - g(2) * (pkin(4) * t14 + pkin(10) * t13 + t60) 0, 0, 0, 0, 0, 0, -g(1) * t113 - g(2) * t6, g(1) * t114 - g(2) * t5, -t75, -g(1) * (-t29 * t85 + t43 * t69 - t57 * t68 + t62) - g(2) * (-t13 * t57 + t14 * t43 + t31 * t85 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t61, -g(1) * t78 - g(2) * t79 - g(3) * t89, 0, 0, 0, 0, 0, 0, -t61 * t50, -t9, -t10, -g(1) * (t78 - t107) - g(2) * (t79 - t108) - t77, 0, 0, 0, 0, 0, 0, -g(1) * (t31 * t91 - t101) - g(2) * (t29 * t91 - t106) - (t50 * t90 + t93) * t109, -g(1) * (-t30 * t53 - t31 * t94) - g(2) * (-t28 * t53 - t29 * t94) - (-t49 * t92 + t53 * t55) * t109, t9, -g(1) * t82 - g(2) * t83 - g(3) * (t76 * t98 + t86) + t117 * (qJ(3) + t76) 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t44 + t31 * t99) - g(2) * (-t28 * t44 + t29 * t99) - (t44 * t55 + t45 * t92) * t109, -g(1) * (-t100 * t31 - t30 * t45) - g(2) * (-t100 * t29 - t28 * t45) - (-t44 * t92 + t45 * t55) * t109, t9, -g(1) * (-pkin(5) * t101 + t82) - g(2) * (-pkin(5) * t106 + t83) - t77 - (pkin(5) * t93 + t51 * t72) * t109 + t117 * (qJ(3) + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t63, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t53, -t64 * t49, -t63, -g(1) * (-pkin(4) * t13 + pkin(10) * t14) - g(2) * (pkin(4) * t68 - pkin(10) * t69) - g(3) * (pkin(4) * t26 + pkin(10) * t27) 0, 0, 0, 0, 0, 0, t64 * t45, -t64 * t44, -t63, -g(1) * (-t13 * t43 - t14 * t57) - g(2) * (t43 * t68 + t57 * t69) - g(3) * (t26 * t43 - t27 * t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, g(1) * t8 - g(2) * t115 - g(3) * (-t27 * t53 - t49 * t98) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t118 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
