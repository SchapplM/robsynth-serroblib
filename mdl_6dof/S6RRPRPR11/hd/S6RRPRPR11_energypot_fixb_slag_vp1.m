% Calculate potential energy for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:35
% EndTime: 2019-03-09 11:11:36
% DurationCPUTime: 0.54s
% Computational Cost: add. (186->110), mult. (237->135), div. (0->0), fcn. (221->10), ass. (0->45)
t112 = rSges(5,3) + pkin(8);
t74 = -qJ(5) - pkin(8);
t111 = rSges(6,3) - t74;
t110 = rSges(7,3) + pkin(9) - t74;
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t109 = g(1) * t80 + g(2) * t77;
t75 = sin(qJ(4));
t106 = t75 * pkin(4);
t76 = sin(qJ(2));
t104 = t76 * pkin(2) + pkin(6);
t78 = cos(qJ(4));
t63 = t78 * pkin(4) + pkin(3);
t103 = t76 * t77;
t73 = qJ(4) + pkin(10);
t66 = qJ(6) + t73;
t61 = sin(t66);
t102 = t77 * t61;
t62 = cos(t66);
t101 = t77 * t62;
t64 = sin(t73);
t100 = t77 * t64;
t65 = cos(t73);
t99 = t77 * t65;
t98 = t77 * t75;
t97 = t77 * t78;
t79 = cos(qJ(2));
t96 = t77 * t79;
t95 = t80 * t61;
t94 = t80 * t62;
t93 = t80 * t64;
t92 = t80 * t65;
t91 = t80 * t75;
t90 = t80 * t78;
t87 = t80 * pkin(1) + t77 * pkin(7);
t86 = qJ(3) * t76;
t85 = t76 * t98;
t84 = t76 * t91;
t69 = t77 * pkin(1);
t83 = pkin(2) * t96 + t77 * t86 + t69;
t82 = t87 + (pkin(2) * t79 + t86) * t80;
t81 = -t80 * pkin(7) + t83;
t56 = pkin(5) * t64 + t106;
t55 = pkin(5) * t65 + t63;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t80 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + t80 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t77 * rSges(3,3) + t87) + g(2) * (rSges(3,1) * t96 - rSges(3,2) * t103 + t69) + g(3) * (t76 * rSges(3,1) + t79 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t79 - rSges(3,2) * t76) + g(2) * (-rSges(3,3) - pkin(7))) * t80) - m(4) * (g(1) * (t77 * rSges(4,1) + t82) + g(2) * (-rSges(4,2) * t96 + rSges(4,3) * t103 + t83) + g(3) * (-t76 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t79 + t104) + (g(1) * (-rSges(4,2) * t79 + rSges(4,3) * t76) + g(2) * (-rSges(4,1) - pkin(7))) * t80) - m(5) * (g(1) * (t77 * pkin(3) + (t84 + t97) * rSges(5,1) + (t76 * t90 - t98) * rSges(5,2) + t82) + g(2) * (-t80 * pkin(3) + (t85 - t90) * rSges(5,1) + (t76 * t97 + t91) * rSges(5,2) + t81) + g(3) * (t112 * t76 + t104) + (g(3) * (-rSges(5,1) * t75 - rSges(5,2) * t78 - qJ(3)) + t109 * t112) * t79) - m(6) * (g(1) * (t77 * t63 + pkin(4) * t84 + (t76 * t93 + t99) * rSges(6,1) + (t76 * t92 - t100) * rSges(6,2) + t82) + g(2) * (-t80 * t63 + pkin(4) * t85 + (t76 * t100 - t92) * rSges(6,1) + (t76 * t99 + t93) * rSges(6,2) + t81) + g(3) * (t111 * t76 + t104) + (g(3) * (-rSges(6,1) * t64 - rSges(6,2) * t65 - qJ(3) - t106) + t109 * t111) * t79) - m(7) * (g(1) * (t77 * t55 + t80 * t76 * t56 + (t76 * t95 + t101) * rSges(7,1) + (t76 * t94 - t102) * rSges(7,2) + t82) + g(2) * (-t80 * t55 + t56 * t103 + (t76 * t102 - t94) * rSges(7,1) + (t76 * t101 + t95) * rSges(7,2) + t81) + g(3) * (t110 * t76 + t104) + (g(3) * (-rSges(7,1) * t61 - rSges(7,2) * t62 - qJ(3) - t56) + t109 * t110) * t79);
U  = t1;
