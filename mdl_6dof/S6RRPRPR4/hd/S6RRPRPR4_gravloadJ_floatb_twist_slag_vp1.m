% Calculate Gravitation load on the joints for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:20
% EndTime: 2019-03-09 10:21:24
% DurationCPUTime: 1.45s
% Computational Cost: add. (829->189), mult. (1725->275), div. (0->0), fcn. (2132->14), ass. (0->82)
t101 = cos(pkin(11));
t62 = sin(pkin(11));
t68 = sin(qJ(2));
t72 = cos(qJ(2));
t41 = -t68 * t101 - t72 * t62;
t73 = cos(qJ(1));
t106 = t72 * t73;
t69 = sin(qJ(1));
t109 = t69 * t68;
t64 = cos(pkin(6));
t139 = t64 * t106 - t109;
t63 = sin(pkin(6));
t115 = t63 * t69;
t102 = t41 * t64;
t78 = t72 * t101 - t68 * t62;
t25 = t102 * t69 + t73 * t78;
t67 = sin(qJ(4));
t71 = cos(qJ(4));
t10 = t71 * t115 - t25 * t67;
t104 = t73 * t68;
t108 = t69 * t72;
t38 = -t64 * t108 - t104;
t134 = t38 * pkin(2);
t74 = t64 * t78;
t24 = t41 * t73 - t69 * t74;
t57 = pkin(4) * t71 + pkin(3);
t65 = -qJ(5) - pkin(9);
t138 = t24 * t57 - t25 * t65 + t134;
t34 = t41 * t63;
t137 = t34 * t67 + t64 * t71;
t21 = t69 * t41 + t73 * t74;
t114 = t63 * t73;
t20 = t102 * t73 - t69 * t78;
t61 = qJ(4) + pkin(12);
t59 = sin(t61);
t60 = cos(t61);
t5 = -t59 * t114 - t20 * t60;
t66 = sin(qJ(6));
t70 = cos(qJ(6));
t136 = t21 * t70 + t5 * t66;
t135 = t21 * t66 - t5 * t70;
t81 = t71 * t114 - t20 * t67;
t133 = t81 * pkin(4);
t124 = pkin(9) + rSges(5,3);
t123 = pkin(10) + rSges(7,3);
t33 = t78 * t63;
t132 = -g(1) * t24 - g(2) * t21 - g(3) * t33;
t131 = -m(6) - m(7);
t130 = pkin(2) * t72;
t126 = t60 * pkin(5);
t117 = t60 * t66;
t116 = t60 * t70;
t58 = pkin(1) + t130;
t110 = t69 * t58;
t100 = t67 * t115;
t53 = t67 * t114;
t54 = t63 * t130;
t95 = t33 * t57 + t34 * t65 + t54;
t35 = pkin(2) * t64 * t68 + (-pkin(8) - qJ(3)) * t63;
t94 = rSges(4,3) * t63 - t35;
t93 = -t60 * t114 + t20 * t59;
t92 = t20 * t71 + t53;
t52 = t73 * t58;
t91 = -t35 * t69 + t52;
t88 = t139 * pkin(2);
t87 = t10 * pkin(4);
t86 = t137 * pkin(4);
t85 = rSges(6,1) * t60 - rSges(6,2) * t59;
t84 = -t73 * t35 - t110;
t80 = t20 * t65 + t21 * t57 + t88;
t79 = pkin(4) * t100 + t24 * t65 + t25 * t57 + t91;
t75 = pkin(4) * t53 + t20 * t57 - t21 * t65 + t84;
t39 = -t109 * t64 + t106;
t37 = -t104 * t64 - t108;
t27 = -t34 * t60 + t59 * t64;
t26 = t34 * t59 + t60 * t64;
t11 = t25 * t71 + t100;
t9 = t115 * t59 + t25 * t60;
t8 = -t115 * t60 + t25 * t59;
t2 = -t24 * t66 + t70 * t9;
t1 = -t24 * t70 - t66 * t9;
t3 = [-m(2) * (g(1) * (-t69 * rSges(2,1) - rSges(2,2) * t73) + g(2) * (rSges(2,1) * t73 - t69 * rSges(2,2))) - m(3) * (g(1) * (t37 * rSges(3,1) - rSges(3,2) * t139 - t69 * pkin(1)) + g(2) * (t39 * rSges(3,1) + t38 * rSges(3,2) + pkin(1) * t73) + (g(1) * t73 + g(2) * t69) * t63 * (rSges(3,3) + pkin(8))) - m(4) * (g(1) * (rSges(4,1) * t20 - t21 * rSges(4,2) + t73 * t94 - t110) + g(2) * (rSges(4,1) * t25 + rSges(4,2) * t24 + t69 * t94 + t52)) - m(5) * (g(1) * (rSges(5,1) * t92 + rSges(5,2) * t81 + pkin(3) * t20 + t124 * t21 + t84) + g(2) * (rSges(5,1) * t11 + rSges(5,2) * t10 + pkin(3) * t25 - t124 * t24 + t91)) - m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t93 + t21 * rSges(6,3) + t75) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 - rSges(6,3) * t24 + t79)) - m(7) * (g(1) * (t135 * rSges(7,1) + t136 * rSges(7,2) - t5 * pkin(5) + t123 * t93 + t75) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t9 + t123 * t8 + t79)) -m(3) * (g(1) * (rSges(3,1) * t38 - rSges(3,2) * t39) + g(2) * (rSges(3,1) * t139 + rSges(3,2) * t37) + g(3) * (rSges(3,1) * t72 - rSges(3,2) * t68) * t63) - m(4) * (g(1) * (t24 * rSges(4,1) - rSges(4,2) * t25 + t134) + g(2) * (rSges(4,1) * t21 + rSges(4,2) * t20 + t88) + g(3) * (rSges(4,1) * t33 + rSges(4,2) * t34 + t54)) - m(5) * (g(1) * (t124 * t25 + t134) + g(2) * (-t124 * t20 + t88) + g(3) * (-t124 * t34 + t54) - t132 * (t71 * rSges(5,1) - t67 * rSges(5,2) + pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t25 + t24 * t85 + t138) + g(2) * (-rSges(6,3) * t20 + t21 * t85 + t80) + g(3) * (-rSges(6,3) * t34 + t33 * t85 + t95)) - m(7) * (g(1) * (t24 * t126 + (t116 * t24 + t25 * t66) * rSges(7,1) + (-t117 * t24 + t25 * t70) * rSges(7,2) + t138) + g(2) * (t21 * t126 + (t116 * t21 - t20 * t66) * rSges(7,1) + (-t117 * t21 - t20 * t70) * rSges(7,2) + t80) + g(3) * (t33 * t126 + (t116 * t33 - t34 * t66) * rSges(7,1) + (-t117 * t33 - t34 * t70) * rSges(7,2) + t95) - t132 * t59 * t123) (-m(4) - m(5) + t131) * (g(3) * t64 + (g(1) * t69 - g(2) * t73) * t63) -m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (-t81 * rSges(5,1) + t92 * rSges(5,2)) + g(3) * (t137 * rSges(5,1) + (t34 * t71 - t64 * t67) * rSges(5,2))) - m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9 + t87) + g(2) * (rSges(6,1) * t93 - t5 * rSges(6,2) - t133) + g(3) * (rSges(6,1) * t26 - rSges(6,2) * t27 + t86)) + (-g(1) * (t123 * t9 + t87) - g(2) * (t123 * t5 - t133) - g(3) * (t123 * t27 + t86) - (-g(1) * t8 + g(2) * t93 + g(3) * t26) * (t70 * rSges(7,1) - t66 * rSges(7,2) + pkin(5))) * m(7), t131 * t132, -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t136 * rSges(7,1) + t135 * rSges(7,2)) + g(3) * ((-t27 * t66 - t33 * t70) * rSges(7,1) + (-t27 * t70 + t33 * t66) * rSges(7,2)))];
taug  = t3(:);
