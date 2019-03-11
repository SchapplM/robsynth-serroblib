% Calculate potential energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:10
% EndTime: 2019-03-09 17:14:10
% DurationCPUTime: 0.46s
% Computational Cost: add. (194->107), mult. (362->137), div. (0->0), fcn. (388->8), ass. (0->43)
t88 = sin(qJ(1));
t92 = cos(qJ(1));
t116 = g(1) * t92 + g(2) * t88;
t91 = cos(qJ(2));
t113 = g(3) * t91;
t112 = -rSges(6,3) - pkin(9);
t87 = sin(qJ(2));
t111 = t87 * pkin(2) + pkin(6);
t90 = cos(qJ(3));
t103 = t92 * t90;
t105 = t88 * t91;
t86 = sin(qJ(3));
t68 = t86 * t105 + t103;
t85 = sin(qJ(5));
t110 = t68 * t85;
t104 = t92 * t86;
t70 = t91 * t104 - t88 * t90;
t109 = t70 * t85;
t108 = t85 * t86;
t107 = t87 * t88;
t106 = t87 * t92;
t102 = -rSges(7,3) - qJ(6) - pkin(9);
t101 = t92 * pkin(1) + t88 * pkin(7);
t100 = rSges(5,3) + qJ(4);
t99 = t111 + (pkin(3) * t90 + qJ(4) * t86) * t87;
t98 = t92 * t91 * pkin(2) + pkin(8) * t106 + t101;
t71 = t91 * t103 + t88 * t86;
t97 = t71 * pkin(3) + t98;
t82 = t88 * pkin(1);
t96 = pkin(2) * t105 - t92 * pkin(7) + pkin(8) * t107 + t82;
t69 = t90 * t105 - t104;
t95 = t69 * pkin(3) + t96;
t94 = t70 * qJ(4) + t97;
t93 = t68 * qJ(4) + t95;
t89 = cos(qJ(5));
t79 = t89 * pkin(5) + pkin(4);
t65 = (t89 * t90 + t108) * t87;
t64 = (-t85 * t90 + t86 * t89) * t87;
t63 = t71 * t89 + t109;
t62 = t70 * t89 - t71 * t85;
t61 = t69 * t89 + t110;
t60 = t68 * t89 - t69 * t85;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t92 * rSges(2,1) - t88 * rSges(2,2)) + g(2) * (t88 * rSges(2,1) + t92 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t88 * rSges(3,3) + t101) + g(2) * (rSges(3,1) * t105 - rSges(3,2) * t107 + t82) + g(3) * (t87 * rSges(3,1) + t91 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t91 - rSges(3,2) * t87) + g(2) * (-rSges(3,3) - pkin(7))) * t92) - m(4) * (g(1) * (t71 * rSges(4,1) - t70 * rSges(4,2) + rSges(4,3) * t106 + t98) + g(2) * (t69 * rSges(4,1) - t68 * rSges(4,2) + rSges(4,3) * t107 + t96) + g(3) * ((-rSges(4,3) - pkin(8)) * t91 + (rSges(4,1) * t90 - rSges(4,2) * t86) * t87 + t111)) - m(5) * (g(1) * (t71 * rSges(5,1) + rSges(5,2) * t106 + t100 * t70 + t97) + g(2) * (t69 * rSges(5,1) + rSges(5,2) * t107 + t100 * t68 + t95) + g(3) * ((-rSges(5,2) - pkin(8)) * t91 + (rSges(5,1) * t90 + rSges(5,3) * t86) * t87 + t99)) - m(6) * (g(1) * (t63 * rSges(6,1) + t62 * rSges(6,2) + t71 * pkin(4) + t94) + g(2) * (t61 * rSges(6,1) + t60 * rSges(6,2) + t69 * pkin(4) + t93) + g(3) * (t65 * rSges(6,1) + t64 * rSges(6,2) + t99) + (-pkin(8) - t112) * t113 + (g(3) * pkin(4) * t90 + t116 * t112) * t87) - m(7) * (g(1) * (t63 * rSges(7,1) + t62 * rSges(7,2) + pkin(5) * t109 + t71 * t79 + t94) + g(2) * (t61 * rSges(7,1) + t60 * rSges(7,2) + pkin(5) * t110 + t69 * t79 + t93) + g(3) * (t65 * rSges(7,1) + t64 * rSges(7,2) + t99) + (-pkin(8) - t102) * t113 + (g(3) * (pkin(5) * t108 + t79 * t90) + t116 * t102) * t87);
U  = t1;
