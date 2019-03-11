% Calculate Gravitation load on the joints for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:09
% EndTime: 2019-03-09 13:38:13
% DurationCPUTime: 1.67s
% Computational Cost: add. (972->212), mult. (2237->320), div. (0->0), fcn. (2818->14), ass. (0->94)
t100 = sin(pkin(12));
t102 = cos(pkin(12));
t122 = cos(qJ(2));
t59 = sin(qJ(2));
t41 = -t100 * t122 - t102 * t59;
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t103 = cos(pkin(6));
t78 = t103 * t100;
t79 = t103 * t102;
t65 = t122 * t79 - t59 * t78;
t25 = t41 * t63 - t60 * t65;
t86 = t103 * t122;
t38 = -t63 * t59 - t60 * t86;
t68 = t38 * pkin(2);
t134 = t25 * pkin(3) + t68;
t133 = -t60 * t59 + t63 * t86;
t106 = -t122 * t78 - t59 * t79;
t40 = t100 * t59 - t102 * t122;
t21 = t106 * t63 + t60 * t40;
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t101 = sin(pkin(6));
t90 = t63 * t101;
t12 = -t21 * t62 - t58 * t90;
t22 = t60 * t41 + t63 * t65;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t132 = t12 * t57 + t22 * t61;
t131 = -t12 * t61 + t22 * t57;
t76 = t101 * t100;
t77 = t102 * t101;
t34 = t122 * t76 + t59 * t77;
t28 = t103 * t58 + t34 * t62;
t130 = g(3) * t28;
t26 = t106 * t60 - t63 * t40;
t91 = t60 * t101;
t15 = t26 * t58 - t62 * t91;
t27 = t103 * t62 - t34 * t58;
t93 = t21 * t58 - t62 * t90;
t129 = -g(1) * t15 + g(2) * t93 + g(3) * t27;
t33 = -t122 * t77 + t59 * t76;
t126 = g(3) * t33;
t125 = t62 * pkin(4);
t124 = rSges(5,3) + pkin(9);
t123 = pkin(10) + rSges(6,3);
t121 = t21 * t57;
t118 = t26 * t57;
t117 = t34 * t57;
t56 = qJ(5) + qJ(6);
t54 = sin(t56);
t116 = t54 * t62;
t55 = cos(t56);
t115 = t55 * t62;
t114 = t57 * t62;
t53 = pkin(2) * t122 + pkin(1);
t112 = t60 * t53;
t110 = t61 * t62;
t52 = pkin(5) * t61 + pkin(4);
t109 = t62 * t52;
t85 = t122 * t101;
t50 = pkin(2) * t85;
t105 = -t33 * pkin(3) + t50;
t104 = pkin(11) + pkin(10) + rSges(7,3);
t99 = g(1) * t123;
t98 = g(2) * t123;
t97 = -t57 * pkin(5) - pkin(9);
t96 = g(1) * t104;
t95 = g(2) * t104;
t94 = t101 * pkin(8);
t92 = t59 * t103;
t89 = t133 * pkin(2);
t88 = t34 * pkin(9) + t105;
t35 = pkin(2) * t92 - qJ(3) * t101 - t94;
t46 = t63 * t53;
t87 = t26 * pkin(3) - t60 * t35 + t46;
t84 = rSges(5,1) * t62 - rSges(5,2) * t58;
t16 = t26 * t62 + t58 * t91;
t7 = -t16 * t57 - t25 * t61;
t82 = -t28 * t57 + t33 * t61;
t81 = rSges(4,3) * t101 - t35;
t80 = t22 * pkin(3) + t89;
t72 = t55 * rSges(7,1) - t54 * rSges(7,2) + t52;
t71 = pkin(3) * t21 - t63 * t35 - t112;
t70 = -t21 * pkin(9) + t80;
t69 = rSges(3,3) * t101 + t94;
t5 = -t16 * t54 - t25 * t55;
t6 = t16 * t55 - t25 * t54;
t67 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * ((-t12 * t54 - t22 * t55) * rSges(7,1) + (-t12 * t55 + t22 * t54) * rSges(7,2)) + g(3) * ((-t28 * t54 + t33 * t55) * rSges(7,1) + (-t28 * t55 - t33 * t54) * rSges(7,2)));
t66 = pkin(9) * t26 + t134;
t39 = t122 * t63 - t60 * t92;
t37 = -t122 * t60 - t63 * t92;
t8 = t16 * t61 - t25 * t57;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t60 - rSges(2,2) * t63) + g(2) * (rSges(2,1) * t63 - rSges(2,2) * t60)) - m(3) * (g(1) * (t37 * rSges(3,1) - rSges(3,2) * t133 - t60 * pkin(1) + t63 * t69) + g(2) * (t39 * rSges(3,1) + t38 * rSges(3,2) + t63 * pkin(1) + t60 * t69)) - m(4) * (g(1) * (rSges(4,1) * t21 - t22 * rSges(4,2) + t63 * t81 - t112) + g(2) * (t26 * rSges(4,1) + t25 * rSges(4,2) + t60 * t81 + t46)) - m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t93 + t124 * t22 + t71) + g(2) * (rSges(5,1) * t16 - rSges(5,2) * t15 - t124 * t25 + t87)) - m(6) * (g(1) * (t131 * rSges(6,1) + t132 * rSges(6,2) - t12 * pkin(4) + t22 * pkin(9) + t123 * t93 + t71) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 + pkin(4) * t16 - t25 * pkin(9) + t123 * t15 + t87)) - m(7) * (g(1) * (-t72 * t12 + (t54 * rSges(7,1) + t55 * rSges(7,2) - t97) * t22 + t104 * t93 + t71) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t104 * t15 + t16 * t52 + t25 * t97 + t87)) -m(3) * (g(1) * (rSges(3,1) * t38 - rSges(3,2) * t39) + g(2) * (rSges(3,1) * t133 + rSges(3,2) * t37) + g(3) * (-rSges(3,2) * t101 * t59 + rSges(3,1) * t85)) - m(4) * (g(1) * (t25 * rSges(4,1) - rSges(4,2) * t26 + t68) + g(2) * (rSges(4,1) * t22 + rSges(4,2) * t21 + t89) + g(3) * (-rSges(4,1) * t33 - rSges(4,2) * t34 + t50)) - m(5) * (g(1) * (t124 * t26 + t25 * t84 + t134) + g(2) * (-t124 * t21 + t22 * t84 + t80) + g(3) * (t124 * t34 - t33 * t84 + t105)) - m(6) * (g(1) * (t25 * t125 + (t110 * t25 + t118) * rSges(6,1) + (-t114 * t25 + t26 * t61) * rSges(6,2) + t66) + g(2) * (t22 * t125 + (t110 * t22 - t121) * rSges(6,1) + (-t114 * t22 - t21 * t61) * rSges(6,2) + t70) + g(3) * (-t33 * t125 + (-t110 * t33 + t117) * rSges(6,1) + (t114 * t33 + t34 * t61) * rSges(6,2) + t88) + (-t123 * t126 + t22 * t98 + t25 * t99) * t58) - m(7) * (g(1) * (t25 * t109 + pkin(5) * t118 + (t115 * t25 + t26 * t54) * rSges(7,1) + (-t116 * t25 + t26 * t55) * rSges(7,2) + t66) + g(2) * (t22 * t109 - pkin(5) * t121 + (t115 * t22 - t21 * t54) * rSges(7,1) + (-t116 * t22 - t21 * t55) * rSges(7,2) + t70) + g(3) * (-t33 * t109 + pkin(5) * t117 + (-t115 * t33 + t34 * t54) * rSges(7,1) + (t116 * t33 + t34 * t55) * rSges(7,2) + t88) + (-t104 * t126 + t22 * t95 + t25 * t96) * t58) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t91 - g(2) * t90 + g(3) * t103) -m(5) * (g(1) * (-rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (rSges(5,1) * t93 - rSges(5,2) * t12) + g(3) * (rSges(5,1) * t27 - rSges(5,2) * t28)) - m(6) * (t123 * t130 + t12 * t98 + t16 * t99 + t129 * (t61 * rSges(6,1) - t57 * rSges(6,2) + pkin(4))) - m(7) * (t104 * t130 + t12 * t95 + t129 * t72 + t16 * t96) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t132 + t131 * rSges(6,2)) + g(3) * (t82 * rSges(6,1) + (-t28 * t61 - t33 * t57) * rSges(6,2))) - t67 - m(7) * (g(1) * t7 - g(2) * t132 + g(3) * t82) * pkin(5), -t67];
taug  = t1(:);
