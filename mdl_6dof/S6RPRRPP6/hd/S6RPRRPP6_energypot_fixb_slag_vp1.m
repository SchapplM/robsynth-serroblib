% Calculate potential energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:16
% EndTime: 2019-03-09 04:46:17
% DurationCPUTime: 0.46s
% Computational Cost: add. (171->98), mult. (215->115), div. (0->0), fcn. (203->8), ass. (0->42)
t101 = rSges(7,1) + pkin(5);
t100 = rSges(5,3) + pkin(8);
t70 = -qJ(5) - pkin(8);
t99 = rSges(7,2) - t70;
t98 = rSges(7,3) + qJ(6);
t97 = pkin(2) + pkin(6);
t73 = sin(qJ(1));
t96 = g(1) * t73;
t76 = cos(qJ(1));
t95 = g(2) * t76;
t71 = sin(qJ(4));
t93 = t71 * t76;
t72 = sin(qJ(3));
t92 = t72 * t73;
t91 = t72 * t76;
t69 = qJ(4) + pkin(9);
t62 = sin(t69);
t90 = t73 * t62;
t63 = cos(t69);
t89 = t73 * t63;
t88 = t73 * t71;
t74 = cos(qJ(4));
t87 = t73 * t74;
t75 = cos(qJ(3));
t86 = t73 * t75;
t85 = t74 * t76;
t84 = t76 * t63;
t83 = rSges(6,3) - t70;
t82 = t76 * pkin(1) + t73 * qJ(2);
t66 = t73 * pkin(1);
t81 = t73 * pkin(7) + t66;
t61 = pkin(4) * t74 + pkin(3);
t80 = t75 * t61 + t97;
t79 = t76 * pkin(7) + t82;
t78 = -t76 * qJ(2) + t81;
t77 = pkin(4) * t93 + t61 * t92 + t70 * t86 + t79;
t59 = pkin(4) * t88;
t54 = -t72 * t84 + t90;
t53 = t62 * t91 + t89;
t52 = t62 * t76 + t72 * t89;
t51 = t72 * t90 - t84;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t76 - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) + rSges(2,2) * t76) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t76 + t73 * rSges(3,3) + t82) + g(2) * (-t73 * rSges(3,2) + t66 + (-rSges(3,3) - qJ(2)) * t76) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t92 + rSges(4,2) * t86 + t79) + g(2) * (t73 * rSges(4,3) + t81) + g(3) * (rSges(4,1) * t75 - rSges(4,2) * t72 + t97) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t72 - rSges(4,2) * t75 - qJ(2))) * t76) - m(5) * (g(1) * (pkin(3) * t92 + (t72 * t87 + t93) * rSges(5,1) + (-t72 * t88 + t85) * rSges(5,2) + t79) + g(2) * (-pkin(3) * t91 + (-t72 * t85 + t88) * rSges(5,1) + (t71 * t91 + t87) * rSges(5,2) + t78) + g(3) * (t100 * t72 + t97) + (g(3) * (rSges(5,1) * t74 - rSges(5,2) * t71 + pkin(3)) + (-t96 + t95) * t100) * t75) - m(6) * (g(1) * (rSges(6,1) * t52 - rSges(6,2) * t51 - rSges(6,3) * t86 + t77) + g(2) * (t54 * rSges(6,1) + t53 * rSges(6,2) + t59 + t81) + g(3) * ((rSges(6,1) * t63 - rSges(6,2) * t62) * t75 + t83 * t72 + t80) + (-t61 * t72 + t83 * t75 - qJ(2)) * t95) - m(7) * (g(1) * (t101 * t52 + t98 * t51 + t77) + g(2) * (t101 * t54 - t98 * t53 - t61 * t91 + t59 + t78) + g(3) * (t99 * t72 + t80) + (-rSges(7,2) * t96 + g(3) * (t101 * t63 + t98 * t62) + t99 * t95) * t75);
U  = t1;
