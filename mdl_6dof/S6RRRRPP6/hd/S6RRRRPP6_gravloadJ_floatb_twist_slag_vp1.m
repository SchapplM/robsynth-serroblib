% Calculate Gravitation load on the joints for
% S6RRRRPP6
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:10
% EndTime: 2019-03-09 21:11:13
% DurationCPUTime: 1.08s
% Computational Cost: add. (575->175), mult. (771->225), div. (0->0), fcn. (782->8), ass. (0->75)
t83 = rSges(7,3) + qJ(6);
t47 = sin(qJ(3));
t49 = sin(qJ(1));
t50 = cos(qJ(3));
t51 = cos(qJ(2));
t52 = cos(qJ(1));
t92 = t51 * t52;
t24 = -t47 * t92 + t49 * t50;
t46 = qJ(3) + qJ(4);
t41 = sin(t46);
t42 = cos(t46);
t93 = t49 * t51;
t18 = t41 * t93 + t42 * t52;
t91 = t52 * t41;
t19 = t42 * t93 - t91;
t115 = t18 * rSges(6,2) + t19 * rSges(6,3);
t40 = pkin(3) * t50 + pkin(2);
t29 = t51 * t40;
t76 = -pkin(1) - t29;
t20 = -t49 * t42 + t51 * t91;
t21 = t41 * t49 + t42 * t92;
t114 = t20 * rSges(6,2) + t21 * rSges(6,3);
t113 = rSges(7,2) + qJ(5);
t104 = g(2) * t49;
t112 = g(1) * t52 + t104;
t111 = t19 * rSges(7,2) - t83 * t18;
t110 = t21 * rSges(7,2) - t83 * t20;
t100 = rSges(4,3) + pkin(8);
t48 = sin(qJ(2));
t109 = t51 * pkin(2) + t100 * t48;
t108 = -t18 * rSges(5,1) - t19 * rSges(5,2);
t107 = pkin(3) * t47;
t106 = g(1) * t49;
t103 = g(3) * t48;
t101 = -rSges(7,1) - pkin(5);
t99 = rSges(3,2) * t48;
t98 = t41 * t48;
t97 = t42 * t48;
t96 = t47 * t49;
t95 = t47 * t52;
t53 = -pkin(9) - pkin(8);
t90 = rSges(6,1) - t53;
t89 = rSges(5,3) - t53;
t88 = -t20 * rSges(5,1) - t21 * rSges(5,2);
t87 = t52 * pkin(1) + t49 * pkin(7);
t84 = rSges(6,3) + qJ(5);
t82 = -t53 - t101;
t80 = -pkin(4) - t83;
t27 = qJ(5) * t97;
t79 = rSges(6,2) * t98 + rSges(6,3) * t97 + t27;
t44 = t52 * pkin(7);
t78 = t49 * t48 * t53 + pkin(3) * t95 + t44;
t77 = g(3) * (rSges(7,2) * t97 + t27);
t75 = t90 * t52;
t74 = t89 * t52;
t73 = -t18 * pkin(4) + qJ(5) * t19;
t72 = -t20 * pkin(4) + qJ(5) * t21;
t71 = pkin(3) * t96 + t40 * t92 + t87;
t70 = g(3) * (t29 + (pkin(4) * t42 + qJ(5) * t41) * t51);
t69 = t82 * t52;
t68 = t80 * t41;
t67 = t21 * pkin(4) + t71;
t66 = rSges(3,1) * t51 - t99;
t64 = rSges(5,1) * t42 - rSges(5,2) * t41;
t63 = -rSges(5,1) * t41 - rSges(5,2) * t42;
t62 = t24 * pkin(3);
t61 = rSges(4,1) * t50 - rSges(4,2) * t47 + pkin(2);
t22 = t47 * t93 + t50 * t52;
t59 = -t19 * pkin(4) - qJ(5) * t18 + t78;
t58 = t22 * pkin(3);
t57 = t62 + t72;
t54 = -t58 + t73;
t25 = t50 * t92 + t96;
t23 = -t50 * t93 + t95;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t49 - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 - rSges(2,2) * t49)) - m(3) * (g(1) * (rSges(3,3) * t52 + t44) + g(2) * (rSges(3,1) * t92 - t52 * t99 + t87) + (g(1) * (-pkin(1) - t66) + g(2) * rSges(3,3)) * t49) - m(4) * (g(1) * (rSges(4,1) * t23 + rSges(4,2) * t22 + t44) + (-pkin(1) - t109) * t106 + (rSges(4,1) * t25 + rSges(4,2) * t24 + t109 * t52 + t87) * g(2)) - m(5) * (g(1) * (-rSges(5,1) * t19 + rSges(5,2) * t18 + t78) + g(2) * (t21 * rSges(5,1) - t20 * rSges(5,2) + t48 * t74 + t71) + (-t48 * rSges(5,3) + t76) * t106) - m(6) * (g(1) * (rSges(6,2) * t19 - rSges(6,3) * t18 + t59) + g(2) * (-t21 * rSges(6,2) + t84 * t20 + t48 * t75 + t67) + (-t48 * rSges(6,1) + t76) * t106) - m(7) * (g(1) * (-rSges(7,2) * t18 - t19 * t83 + t49 * t76 + t59) + g(2) * (t113 * t20 + t21 * t83 + t67) + (g(2) * t69 + t101 * t106) * t48) -m(3) * (g(3) * t66 + t112 * (-rSges(3,1) * t48 - rSges(3,2) * t51)) - m(4) * ((g(3) * t61 + t100 * t112) * t51 + (g(3) * t100 - t112 * t61) * t48) - m(5) * (g(3) * t29 + (g(1) * t74 + g(3) * t64 + t89 * t104) * t51 + (g(3) * t89 + t112 * (-t40 - t64)) * t48) - m(6) * (t70 + (g(3) * (-rSges(6,2) * t42 + rSges(6,3) * t41) + g(1) * t75 + t90 * t104) * t51 + (g(3) * t90 + t112 * (-t40 + (rSges(6,2) - pkin(4)) * t42 - t84 * t41)) * t48) - m(7) * (t70 + (g(3) * (rSges(7,2) * t41 + t83 * t42) + g(1) * t69 + t82 * t104) * t51 + (g(3) * t82 + t112 * (-t113 * t41 + t80 * t42 - t40)) * t48) -m(4) * (g(1) * (rSges(4,1) * t24 - rSges(4,2) * t25) + g(2) * (-rSges(4,1) * t22 + rSges(4,2) * t23) + (-rSges(4,1) * t47 - rSges(4,2) * t50) * t103) - m(5) * (g(1) * (t62 + t88) + g(2) * (-t58 + t108) + (t63 - t107) * t103) - m(6) * (g(1) * (t57 + t114) + g(2) * (t54 + t115) + g(3) * ((-pkin(4) * t41 - t107) * t48 + t79)) - m(7) * (g(1) * (t57 + t110) + g(2) * (t54 + t111) + t77 + (t68 - t107) * t103) -m(5) * (g(1) * t88 + g(2) * t108 + t63 * t103) - m(6) * (g(1) * (t72 + t114) + g(2) * (t73 + t115) + g(3) * (-pkin(4) * t98 + t79)) - m(7) * (g(1) * (t72 + t110) + g(2) * (t73 + t111) + t77 + t68 * t103) (-m(6) - m(7)) * (g(1) * t20 + g(2) * t18 + g(3) * t98) -m(7) * (g(1) * t21 + g(2) * t19 + g(3) * t97)];
taug  = t1(:);
