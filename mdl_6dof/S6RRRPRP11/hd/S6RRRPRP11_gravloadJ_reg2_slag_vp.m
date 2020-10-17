% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:18:30
% EndTime: 2019-05-07 09:18:33
% DurationCPUTime: 0.77s
% Computational Cost: add. (597->151), mult. (1526->204), div. (0->0), fcn. (1860->10), ass. (0->92)
t62 = sin(qJ(2));
t63 = sin(qJ(1));
t66 = cos(qJ(2));
t112 = cos(qJ(1));
t95 = cos(pkin(6));
t77 = t95 * t112;
t41 = t62 * t77 + t63 * t66;
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t58 = sin(pkin(6));
t91 = t58 * t112;
t23 = t41 * t61 + t65 * t91;
t40 = t63 * t62 - t66 * t77;
t60 = sin(qJ(5));
t64 = cos(qJ(5));
t122 = t23 * t60 + t40 * t64;
t121 = t23 * t64 - t40 * t60;
t109 = t40 * t65;
t96 = qJ(4) * t61;
t120 = -pkin(3) * t109 - t40 * t96;
t86 = t63 * t95;
t42 = t112 * t62 + t66 * t86;
t108 = t42 * t65;
t119 = -pkin(3) * t108 - t42 * t96;
t102 = t60 * t66;
t105 = t58 * t65;
t43 = t112 * t66 - t62 * t86;
t27 = -t63 * t105 + t43 * t61;
t13 = t27 * t64 - t42 * t60;
t107 = t58 * t62;
t38 = t61 * t107 - t95 * t65;
t1 = -g(2) * t121 - g(3) * (t58 * t102 + t38 * t64) - g(1) * t13;
t118 = pkin(4) + pkin(9);
t116 = g(3) * t58;
t115 = t40 * pkin(9);
t114 = t42 * pkin(9);
t56 = t64 * pkin(5) + pkin(4);
t113 = pkin(9) + t56;
t106 = t58 * t63;
t104 = t58 * t66;
t103 = t60 * t61;
t101 = t61 * t64;
t100 = t64 * t66;
t99 = t65 * t66;
t98 = pkin(2) * t104 + pkin(9) * t107;
t97 = t112 * pkin(1) + pkin(8) * t106;
t34 = t40 * pkin(2);
t94 = -t34 + t120;
t36 = t42 * pkin(2);
t93 = -t36 + t119;
t92 = t43 * pkin(2) + t97;
t90 = -t63 * pkin(1) + pkin(8) * t91;
t89 = t41 * pkin(9) - t34;
t88 = t43 * pkin(9) - t36;
t87 = pkin(5) * t60 + qJ(4);
t24 = t41 * t65 - t61 * t91;
t19 = t23 * pkin(3);
t85 = t24 * qJ(4) - t19;
t21 = t27 * pkin(3);
t28 = t61 * t106 + t43 * t65;
t84 = t28 * qJ(4) - t21;
t33 = t38 * pkin(3);
t39 = t62 * t105 + t95 * t61;
t83 = t39 * qJ(4) - t33;
t82 = t28 * pkin(3) + t92;
t81 = t58 * pkin(3) * t99 + t96 * t104 + t98;
t80 = -t41 * pkin(2) + t90;
t79 = -g(1) * t23 + g(2) * t27;
t78 = -g(1) * t24 + g(2) * t28;
t18 = g(1) * t40 - g(2) * t42;
t76 = g(3) * t81;
t75 = -pkin(3) * t24 + t80;
t59 = -qJ(6) - pkin(10);
t74 = pkin(5) * t103 - t59 * t65;
t73 = g(1) * t112 + g(2) * t63;
t72 = t27 * qJ(4) + t82;
t10 = g(1) * t27 + g(2) * t23 + g(3) * t38;
t71 = g(1) * t28 + g(2) * t24 + g(3) * t39;
t70 = -qJ(4) * t23 + t75;
t69 = -g(1) * t42 - g(2) * t40 + g(3) * t104;
t68 = g(1) * t43 + g(2) * t41 + g(3) * t107;
t16 = t69 * t65;
t15 = t69 * t61;
t14 = t27 * t60 + t42 * t64;
t8 = t71 * t64;
t7 = t71 * t60;
t6 = g(1) * t122 - g(2) * t14;
t5 = g(1) * t121 - g(2) * t13;
t4 = -g(1) * (-t42 * t103 + t43 * t64) - g(2) * (-t40 * t103 + t41 * t64) - (t61 * t102 + t62 * t64) * t116;
t3 = -g(1) * (-t42 * t101 - t43 * t60) - g(2) * (-t40 * t101 - t41 * t60) - (t61 * t100 - t60 * t62) * t116;
t2 = g(1) * t14 + g(2) * t122 - g(3) * (t58 * t100 - t38 * t60);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t63 - g(2) * t112, t73, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t43, -t18, -t73 * t58, -g(1) * t90 - g(2) * t97, 0, 0, 0, 0, 0, 0, -t78, t79, t18, -g(1) * (t80 - t115) - g(2) * (t92 + t114) 0, 0, 0, 0, 0, 0, t18, t78, -t79, -g(1) * (t70 - t115) - g(2) * (t72 + t114) 0, 0, 0, 0, 0, 0, t6, t5, -t78, -g(1) * (-pkin(10) * t24 - t118 * t40 + t70) - g(2) * (t28 * pkin(10) + t118 * t42 + t72) 0, 0, 0, 0, 0, 0, t6, t5, -t78, -g(1) * (-t113 * t40 - t23 * t87 + t24 * t59 + t75) - g(2) * (t113 * t42 + t87 * t27 - t28 * t59 + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t68, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, -t68, -g(1) * t88 - g(2) * t89 - g(3) * t98, 0, 0, 0, 0, 0, 0, -t68, t16, -t15, -g(1) * (t88 + t119) - g(2) * (t89 + t120) - t76, 0, 0, 0, 0, 0, 0, t4, t3, -t16, -g(1) * (-pkin(10) * t108 + t118 * t43 + t93) - g(2) * (-pkin(10) * t109 + t118 * t41 + t94) - g(3) * ((pkin(4) * t62 + pkin(10) * t99) * t58 + t81) 0, 0, 0, 0, 0, 0, t4, t3, -t16, -g(1) * (t113 * t43 - t74 * t42 + t93) - g(2) * (t113 * t41 - t74 * t40 + t94) - t76 - (t56 * t62 + t74 * t66) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t71, -g(1) * t84 - g(2) * t85 - g(3) * t83, 0, 0, 0, 0, 0, 0, -t7, -t8, t10, -g(1) * (-t27 * pkin(10) + t84) - g(2) * (-t23 * pkin(10) + t85) - g(3) * (-t38 * pkin(10) + t83) 0, 0, 0, 0, 0, 0, -t7, -t8, t10, -g(1) * (t27 * t59 + t87 * t28 - t21) - g(2) * (t23 * t59 + t87 * t24 - t19) - g(3) * (t38 * t59 + t87 * t39 - t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71;];
taug_reg  = t9;
