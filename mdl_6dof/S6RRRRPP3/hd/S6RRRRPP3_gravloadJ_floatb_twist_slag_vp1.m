% Calculate Gravitation load on the joints for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:50
% EndTime: 2019-03-09 20:53:52
% DurationCPUTime: 0.96s
% Computational Cost: add. (563->170), mult. (715->202), div. (0->0), fcn. (696->8), ass. (0->78)
t117 = rSges(7,1) + pkin(5);
t110 = rSges(7,3) + qJ(6);
t48 = qJ(2) + qJ(3);
t45 = sin(t48);
t52 = cos(qJ(4));
t98 = t45 * t52;
t118 = g(3) * t98;
t51 = sin(qJ(1));
t49 = sin(qJ(4));
t99 = t45 * t49;
t82 = rSges(5,2) * t99;
t46 = cos(t48);
t95 = t46 * t51;
t116 = rSges(5,3) * t95 + t51 * t82;
t81 = rSges(6,2) * t98;
t115 = rSges(6,1) * t95 + t51 * t81;
t54 = cos(qJ(1));
t93 = t46 * t54;
t114 = rSges(5,3) * t93 + t54 * t82;
t113 = rSges(6,1) * t93 + t54 * t81;
t112 = t117 * t95;
t111 = t117 * t93;
t109 = t46 * rSges(4,1) - rSges(4,2) * t45;
t108 = g(1) * t54 + g(2) * t51;
t107 = t45 * t108;
t50 = sin(qJ(2));
t106 = pkin(2) * t50;
t55 = -pkin(8) - pkin(7);
t103 = g(2) * t55;
t102 = qJ(5) * t118;
t39 = t45 * pkin(9);
t41 = t46 * pkin(3);
t101 = rSges(3,3) + pkin(7);
t37 = t45 * rSges(6,1);
t35 = t45 * rSges(5,3);
t97 = t45 * t54;
t96 = t46 * t49;
t94 = t46 * t52;
t92 = t49 * t51;
t91 = t51 * t52;
t90 = t52 * t54;
t89 = t54 * t49;
t88 = t54 * t55;
t87 = rSges(4,3) - t55;
t86 = t39 + t41;
t85 = rSges(7,2) + qJ(5);
t84 = rSges(6,3) + qJ(5);
t80 = -pkin(4) - t110;
t53 = cos(qJ(2));
t47 = t53 * pkin(2);
t44 = t47 + pkin(1);
t17 = t54 * t44;
t79 = pkin(3) * t93 + pkin(9) * t97 + t17;
t77 = -t44 - t41;
t75 = pkin(4) * t94 + qJ(5) * t96 + t86;
t28 = pkin(9) * t95;
t74 = -t51 * t106 + t28;
t32 = pkin(9) * t93;
t73 = -t54 * t106 + t32;
t72 = rSges(3,1) * t53 - rSges(3,2) * t50;
t69 = -rSges(4,1) * t45 - rSges(4,2) * t46;
t68 = t77 - t39;
t67 = pkin(1) + t72;
t6 = t46 * t92 + t90;
t7 = t46 * t91 - t89;
t66 = -t7 * pkin(4) - t6 * qJ(5) - t88;
t8 = t46 * t89 - t91;
t9 = t46 * t90 + t92;
t65 = t9 * pkin(4) + t8 * qJ(5) + t79;
t64 = rSges(5,1) * t94 - rSges(5,2) * t96 + t35 + t86;
t61 = rSges(7,2) * t96 + t110 * t94 + t117 * t45 + t75;
t60 = -rSges(6,2) * t94 + rSges(6,3) * t96 + t37 + t75;
t58 = (-rSges(5,1) * t52 - pkin(3)) * t107;
t57 = (-pkin(4) * t52 - t84 * t49 - pkin(3)) * t107;
t56 = (-t85 * t49 + t80 * t52 - pkin(3)) * t107;
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t51 - rSges(2,2) * t54) + g(2) * (rSges(2,1) * t54 - rSges(2,2) * t51)) - m(3) * ((g(1) * t101 + g(2) * t67) * t54 + (-g(1) * t67 + g(2) * t101) * t51) - m(4) * (g(2) * t17 + (g(1) * t87 + g(2) * t109) * t54 + (g(1) * (-t44 - t109) + g(2) * t87) * t51) - m(5) * (g(1) * (-t7 * rSges(5,1) + t6 * rSges(5,2) - t88) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + rSges(5,3) * t97 + t79) + (g(1) * (t68 - t35) - t103) * t51) - m(6) * (g(1) * (t7 * rSges(6,2) - t6 * rSges(6,3) + t66) + g(2) * (rSges(6,1) * t97 - t9 * rSges(6,2) + t8 * rSges(6,3) + t65) + (g(1) * (t68 - t37) - t103) * t51) - m(7) * (g(1) * (-t6 * rSges(7,2) - t110 * t7 + t66) + g(2) * (t8 * rSges(7,2) + t110 * t9 + t117 * t97 + t65) + (-t103 + (t77 + (-pkin(9) - t117) * t45) * g(1)) * t51) -m(3) * (g(3) * t72 + t108 * (-rSges(3,1) * t50 - rSges(3,2) * t53)) - m(4) * (g(3) * (t47 + t109) + t108 * (t69 - t106)) - m(5) * (g(1) * (t73 + t114) + g(2) * (t74 + t116) + g(3) * (t47 + t64) + t58) - m(6) * (g(1) * (t73 + t113) + g(2) * (t74 + t115) + g(3) * (t47 + t60) + t57) - m(7) * (g(1) * (t73 + t111) + g(2) * (t74 + t112) + g(3) * (t47 + t61) + t56) -m(4) * (g(3) * t109 + t108 * t69) - m(5) * (g(1) * (t32 + t114) + g(2) * (t28 + t116) + g(3) * t64 + t58) - m(6) * (g(1) * (t32 + t113) + g(2) * (t28 + t115) + g(3) * t60 + t57) - m(7) * (g(1) * (t32 + t111) + g(2) * (t28 + t112) + g(3) * t61 + t56) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (rSges(6,2) * t8 + t84 * t9 - t4) + g(2) * (rSges(6,2) * t6 + t84 * t7 - t2) + t102) - m(7) * (g(1) * (-t110 * t8 + t85 * t9 - t4) + g(2) * (-t110 * t6 + t85 * t7 - t2) + t102) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t52 + (m(5) * rSges(5,1) - m(6) * (rSges(6,2) - pkin(4)) - m(7) * t80) * t49) * g(3) * t45 (-m(6) - m(7)) * (g(1) * t8 + g(2) * t6 + g(3) * t99) -m(7) * (g(1) * t9 + g(2) * t7 + t118)];
taug  = t1(:);
