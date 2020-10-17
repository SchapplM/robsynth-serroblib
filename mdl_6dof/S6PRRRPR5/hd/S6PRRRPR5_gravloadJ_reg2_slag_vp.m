% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:05:11
% EndTime: 2019-05-05 08:05:16
% DurationCPUTime: 1.27s
% Computational Cost: add. (1020->186), mult. (2508->298), div. (0->0), fcn. (3195->16), ass. (0->100)
t115 = cos(pkin(6));
t68 = sin(pkin(7));
t103 = t115 * t68;
t114 = cos(pkin(7));
t73 = sin(qJ(3));
t104 = t73 * t114;
t132 = sin(qJ(2));
t133 = cos(qJ(3));
t134 = cos(qJ(2));
t69 = sin(pkin(6));
t38 = t73 * t103 + (t134 * t104 + t132 * t133) * t69;
t108 = t69 * t134;
t52 = -t68 * t108 + t115 * t114;
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t143 = -t38 * t72 + t52 * t75;
t112 = sin(pkin(12));
t105 = t69 * t112;
t113 = cos(pkin(12));
t88 = t115 * t112;
t79 = t113 * t132 + t134 * t88;
t139 = -t68 * t105 + t79 * t114;
t54 = t113 * t134 - t132 * t88;
t22 = t54 * t133 - t139 * t73;
t40 = t114 * t105 + t79 * t68;
t142 = -t22 * t72 + t40 * t75;
t106 = t69 * t113;
t89 = t115 * t113;
t78 = t112 * t132 - t134 * t89;
t140 = t68 * t106 + t78 * t114;
t53 = t112 * t134 + t132 * t89;
t20 = t53 * t133 - t140 * t73;
t39 = -t114 * t106 + t78 * t68;
t141 = -t20 * t72 + t39 * t75;
t138 = -g(1) * t54 - g(2) * t53;
t137 = pkin(9) * t68;
t67 = qJ(4) + pkin(13);
t65 = sin(t67);
t125 = t65 * t68;
t66 = cos(t67);
t124 = t66 * t68;
t71 = sin(qJ(6));
t123 = t66 * t71;
t74 = cos(qJ(6));
t122 = t66 * t74;
t121 = t68 * t72;
t120 = t68 * t75;
t19 = t140 * t133 + t53 * t73;
t64 = t75 * pkin(4) + pkin(3);
t70 = -qJ(5) - pkin(10);
t119 = -t19 * t64 - t20 * t70;
t21 = t139 * t133 + t54 * t73;
t118 = -t21 * t64 - t22 * t70;
t107 = t69 * t132;
t92 = t114 * t133;
t37 = -t133 * t103 + t73 * t107 - t92 * t108;
t117 = -t37 * t64 - t38 * t70;
t102 = t68 * t107;
t116 = pkin(2) * t108 + pkin(9) * t102;
t27 = t53 * t92 - t78 * t73;
t28 = -t53 * t104 - t78 * t133;
t50 = t78 * pkin(2);
t111 = -t27 * t70 + t28 * t64 - t50;
t29 = t54 * t92 - t79 * t73;
t30 = -t54 * t104 - t79 * t133;
t51 = t79 * pkin(2);
t110 = -t29 * t70 + t30 * t64 - t51;
t101 = t141 * pkin(4);
t100 = t142 * pkin(4);
t99 = t143 * pkin(4);
t98 = t53 * t137 - t50;
t97 = t54 * t137 - t51;
t91 = t114 * t132;
t47 = (t133 * t91 + t134 * t73) * t69;
t48 = (t134 * t133 - t73 * t91) * t69;
t90 = t72 * t102;
t94 = pkin(4) * t90 - t47 * t70 + t48 * t64 + t116;
t93 = -pkin(5) * t66 - pkin(11) * t65;
t15 = -t38 * t65 + t52 * t66;
t5 = -t20 * t65 + t39 * t66;
t7 = -t22 * t65 + t40 * t66;
t87 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t16 = t38 * t66 + t52 * t65;
t6 = t20 * t66 + t39 * t65;
t8 = t22 * t66 + t40 * t65;
t86 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t11 = -t54 * t124 + t30 * t65;
t31 = -t66 * t102 + t48 * t65;
t9 = -t53 * t124 + t28 * t65;
t85 = g(1) * t11 + g(2) * t9 + g(3) * t31;
t84 = g(1) * t21 + g(2) * t19 + g(3) * t37;
t83 = g(1) * t22 + g(2) * t20 + g(3) * t38;
t82 = g(1) * t29 + g(2) * t27 + g(3) * t47;
t81 = -g(3) * t107 + t138;
t80 = t138 * (pkin(4) * t72 + pkin(9)) * t68;
t32 = t65 * t102 + t48 * t66;
t12 = t54 * t125 + t30 * t66;
t10 = t53 * t125 + t28 * t66;
t1 = t84 * t65;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t79 + g(2) * t78 - g(3) * t108, -t81, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t28 - g(3) * t48, t82, t81 * t68, -g(1) * t97 - g(2) * t98 - g(3) * t116, 0, 0, 0, 0, 0, 0, -g(1) * (t54 * t121 + t30 * t75) - g(2) * (t53 * t121 + t28 * t75) - g(3) * (t48 * t75 + t90) -g(1) * (t54 * t120 - t30 * t72) - g(2) * (t53 * t120 - t28 * t72) - g(3) * (t75 * t102 - t48 * t72) -t82, -g(1) * (t30 * pkin(3) + t29 * pkin(10) + t97) - g(2) * (t28 * pkin(3) + t27 * pkin(10) + t98) - g(3) * (t48 * pkin(3) + t47 * pkin(10) + t116) 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10 - g(3) * t32, t85, -t82, -g(1) * t110 - g(2) * t111 - g(3) * t94 + t80, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t74 + t29 * t71) - g(2) * (t10 * t74 + t27 * t71) - g(3) * (t32 * t74 + t47 * t71) -g(1) * (-t12 * t71 + t29 * t74) - g(2) * (-t10 * t71 + t27 * t74) - g(3) * (-t32 * t71 + t47 * t74) -t85, -g(1) * (t12 * pkin(5) + t11 * pkin(11) + t110) - g(2) * (t10 * pkin(5) + t9 * pkin(11) + t111) - g(3) * (t32 * pkin(5) + t31 * pkin(11) + t94) + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t83, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t75, -t84 * t72, -t83, -g(1) * (-t21 * pkin(3) + t22 * pkin(10)) - g(2) * (-t19 * pkin(3) + t20 * pkin(10)) - g(3) * (-t37 * pkin(3) + t38 * pkin(10)) 0, 0, 0, 0, 0, 0, t84 * t66, -t1, -t83, -g(1) * t118 - g(2) * t119 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t122 + t22 * t71) - g(2) * (-t19 * t122 + t20 * t71) - g(3) * (-t37 * t122 + t38 * t71) -g(1) * (t21 * t123 + t22 * t74) - g(2) * (t19 * t123 + t20 * t74) - g(3) * (t37 * t123 + t38 * t74) t1, -g(1) * (t93 * t21 + t118) - g(2) * (t93 * t19 + t119) - g(3) * (t93 * t37 + t117); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t142 - g(2) * t141 - g(3) * t143, -g(1) * (-t22 * t75 - t40 * t72) - g(2) * (-t20 * t75 - t39 * t72) - g(3) * (-t38 * t75 - t52 * t72) 0, 0, 0, 0, 0, 0, 0, 0, -t87, t86, 0, -g(1) * t100 - g(2) * t101 - g(3) * t99, 0, 0, 0, 0, 0, 0, -t87 * t74, t87 * t71, -t86, -g(1) * (t7 * pkin(5) + t8 * pkin(11) + t100) - g(2) * (t5 * pkin(5) + t6 * pkin(11) + t101) - g(3) * (t15 * pkin(5) + t16 * pkin(11) + t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t21 * t74 - t8 * t71) - g(2) * (t19 * t74 - t6 * t71) - g(3) * (-t16 * t71 + t37 * t74) -g(1) * (-t21 * t71 - t8 * t74) - g(2) * (-t19 * t71 - t6 * t74) - g(3) * (-t16 * t74 - t37 * t71) 0, 0;];
taug_reg  = t2;
