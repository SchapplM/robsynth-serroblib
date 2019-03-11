% Calculate Gravitation load on the joints for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:17
% EndTime: 2019-03-09 11:17:20
% DurationCPUTime: 1.03s
% Computational Cost: add. (578->178), mult. (1091->244), div. (0->0), fcn. (1260->12), ass. (0->72)
t57 = sin(qJ(2));
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t92 = cos(pkin(6));
t82 = t62 * t92;
t32 = t57 * t58 - t61 * t82;
t56 = sin(qJ(4));
t60 = cos(qJ(4));
t53 = sin(pkin(6));
t97 = t53 * t62;
t71 = t32 * t60 + t56 * t97;
t83 = t58 * t92;
t34 = t62 * t57 + t61 * t83;
t99 = t53 * t58;
t11 = t34 * t60 - t56 * t99;
t98 = t53 * t61;
t67 = -t92 * t56 - t60 * t98;
t120 = t67 * pkin(4);
t110 = rSges(7,3) + pkin(10);
t109 = -pkin(9) - rSges(5,3);
t55 = sin(qJ(6));
t59 = cos(qJ(6));
t119 = -rSges(7,1) * t55 - rSges(7,2) * t59;
t33 = t57 * t82 + t58 * t61;
t35 = -t57 * t83 + t61 * t62;
t118 = g(1) * t35 + g(2) * t33;
t52 = qJ(4) + pkin(11);
t49 = sin(t52);
t50 = cos(t52);
t73 = -t32 * t49 + t50 * t97;
t117 = t33 * t59 + t55 * t73;
t116 = -t33 * t55 + t59 * t73;
t115 = -m(6) - m(7);
t114 = pkin(4) * t56;
t111 = g(3) * t53;
t106 = t32 * t56;
t102 = t34 * t56;
t100 = t53 * t57;
t96 = t71 * pkin(4);
t95 = pkin(2) * t98 + qJ(3) * t100;
t94 = t62 * pkin(1) + pkin(8) * t99;
t93 = rSges(4,3) + qJ(3);
t28 = t32 * pkin(2);
t54 = -qJ(5) - pkin(9);
t88 = t33 * t114 + t32 * t54 - t28;
t30 = t34 * pkin(2);
t87 = t35 * t114 + t34 * t54 - t30;
t86 = t35 * pkin(2) + t94;
t85 = -t58 * pkin(1) + pkin(8) * t97;
t84 = t60 * t97 - t106;
t80 = g(3) * (t100 * t114 + t95);
t79 = -t33 * pkin(2) + t85;
t78 = rSges(5,1) * t56 + rSges(5,2) * t60;
t77 = rSges(6,1) * t49 + rSges(6,2) * t50;
t76 = t11 * pkin(4);
t75 = qJ(3) * t34 + t86;
t74 = rSges(7,1) * t59 - rSges(7,2) * t55 + pkin(5);
t72 = t32 * t50 + t49 * t97;
t68 = -t32 * qJ(3) + t79;
t48 = pkin(4) * t60 + pkin(3);
t66 = pkin(4) * t102 - t35 * t54 + t48 * t99 + t75;
t65 = -pkin(4) * t106 + t33 * t54 + t48 * t97 + t68;
t64 = -t110 * t50 + t74 * t49;
t17 = -t49 * t98 + t92 * t50;
t16 = -t92 * t49 - t50 * t98;
t12 = t60 * t99 + t102;
t6 = t34 * t49 + t50 * t99;
t5 = -t34 * t50 + t49 * t99;
t2 = t35 * t55 + t59 * t6;
t1 = t35 * t59 - t55 * t6;
t3 = [-m(2) * (g(1) * (-t58 * rSges(2,1) - rSges(2,2) * t62) + g(2) * (rSges(2,1) * t62 - t58 * rSges(2,2))) - m(3) * (g(1) * (-t33 * rSges(3,1) + t32 * rSges(3,2) + rSges(3,3) * t97 + t85) + g(2) * (rSges(3,1) * t35 - rSges(3,2) * t34 + rSges(3,3) * t99 + t94)) - m(4) * (g(1) * (rSges(4,1) * t97 + t33 * rSges(4,2) - t93 * t32 + t79) + g(2) * (rSges(4,1) * t99 - rSges(4,2) * t35 + t93 * t34 + t86)) - m(5) * (g(1) * (t84 * rSges(5,1) - t71 * rSges(5,2) + pkin(3) * t97 + t109 * t33 + t68) + g(2) * (rSges(5,1) * t12 + rSges(5,2) * t11 + pkin(3) * t99 - t109 * t35 + t75)) - m(6) * (g(1) * (rSges(6,1) * t73 - rSges(6,2) * t72 - rSges(6,3) * t33 + t65) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + rSges(6,3) * t35 + t66)) - m(7) * (g(1) * (t116 * rSges(7,1) - t117 * rSges(7,2) + t73 * pkin(5) + t110 * t72 + t65) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t6 + t110 * t5 + t66)) -m(3) * (g(1) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + g(2) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + (rSges(3,1) * t61 - rSges(3,2) * t57) * t111) - m(4) * (g(1) * (rSges(4,2) * t34 + t93 * t35 - t30) + g(2) * (rSges(4,2) * t32 + t93 * t33 - t28) + g(3) * ((-rSges(4,2) * t61 + rSges(4,3) * t57) * t53 + t95)) - m(5) * (g(1) * (t109 * t34 - t30) + g(2) * (t109 * t32 - t28) + g(3) * t95 + (-t109 * t61 + t78 * t57) * t111 + t118 * (qJ(3) + t78)) - m(6) * (g(1) * (-rSges(6,3) * t34 + t87) + g(2) * (-rSges(6,3) * t32 + t88) + t80 + ((rSges(6,3) - t54) * t61 + t77 * t57) * t111 + t118 * (qJ(3) + t77)) - m(7) * (g(1) * (t119 * t34 + t87) + g(2) * (t119 * t32 + t88) + t80 + ((-t54 - t119) * t61 + t64 * t57) * t111 + t118 * (qJ(3) + t64)) (-m(4) - m(5) + t115) * (g(1) * t34 + g(2) * t32 - g(3) * t98) -m(5) * (g(1) * (rSges(5,1) * t11 - rSges(5,2) * t12) + g(2) * (t71 * rSges(5,1) + t84 * rSges(5,2)) + g(3) * (t67 * rSges(5,1) + (t56 * t98 - t92 * t60) * rSges(5,2))) - m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6 + t76) + g(2) * (rSges(6,1) * t72 + rSges(6,2) * t73 + t96) + g(3) * (t16 * rSges(6,1) - t17 * rSges(6,2) + t120)) + (-g(1) * (t110 * t6 + t76) - g(2) * (-t110 * t73 + t96) - g(3) * (t110 * t17 + t120) - (-g(1) * t5 + g(2) * t72 + g(3) * t16) * t74) * m(7), t115 * (g(3) * t100 + t118) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (t117 * rSges(7,1) + t116 * rSges(7,2)) + g(3) * ((t59 * t100 - t17 * t55) * rSges(7,1) + (-t55 * t100 - t17 * t59) * rSges(7,2)))];
taug  = t3(:);
