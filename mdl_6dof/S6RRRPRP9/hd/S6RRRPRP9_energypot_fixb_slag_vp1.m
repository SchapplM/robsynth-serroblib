% Calculate potential energy for
% S6RRRPRP9
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:46
% EndTime: 2019-03-09 17:21:47
% DurationCPUTime: 0.38s
% Computational Cost: add. (203->103), mult. (392->127), div. (0->0), fcn. (430->8), ass. (0->42)
t129 = rSges(7,1) + pkin(5);
t128 = rSges(7,3) + qJ(6);
t103 = cos(qJ(1));
t98 = sin(qJ(2));
t99 = sin(qJ(1));
t127 = (g(1) * t103 + g(2) * t99) * t98;
t123 = t98 * pkin(2) + pkin(6);
t97 = sin(qJ(3));
t121 = t97 * t98;
t120 = t98 * t99;
t119 = t103 * pkin(1) + t99 * pkin(7);
t101 = cos(qJ(3));
t118 = t101 * t98;
t102 = cos(qJ(2));
t117 = t102 * t99;
t116 = t103 * t98;
t115 = t99 * t101;
t114 = rSges(5,3) + qJ(4);
t113 = t102 * t103;
t112 = t103 * t101;
t111 = pkin(3) * t118 + qJ(4) * t121 + t123;
t110 = pkin(2) * t113 + pkin(8) * t116 + t119;
t81 = t102 * t112 + t99 * t97;
t109 = t81 * pkin(3) + t110;
t93 = t99 * pkin(1);
t108 = pkin(2) * t117 - t103 * pkin(7) + pkin(8) * t120 + t93;
t107 = pkin(4) * t118 + t102 * pkin(9) + t111;
t79 = t102 * t115 - t103 * t97;
t106 = t79 * pkin(3) + t108;
t80 = t97 * t113 - t115;
t105 = t81 * pkin(4) + t80 * qJ(4) + t109;
t78 = t97 * t117 + t112;
t104 = t79 * pkin(4) + t78 * qJ(4) + t106;
t100 = cos(qJ(5));
t96 = sin(qJ(5));
t73 = (t100 * t101 + t96 * t97) * t98;
t72 = -t100 * t121 + t96 * t118;
t69 = t81 * t100 + t80 * t96;
t68 = -t80 * t100 + t81 * t96;
t67 = t79 * t100 + t78 * t96;
t66 = -t78 * t100 + t79 * t96;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t103 * rSges(2,1) - t99 * rSges(2,2)) + g(2) * (t99 * rSges(2,1) + t103 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t99 * rSges(3,3) + t119) + g(2) * (rSges(3,1) * t117 - rSges(3,2) * t120 + t93) + g(3) * (t98 * rSges(3,1) + t102 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t102 - rSges(3,2) * t98) + g(2) * (-rSges(3,3) - pkin(7))) * t103) - m(4) * (g(1) * (t81 * rSges(4,1) - t80 * rSges(4,2) + rSges(4,3) * t116 + t110) + g(2) * (t79 * rSges(4,1) - t78 * rSges(4,2) + rSges(4,3) * t120 + t108) + g(3) * ((rSges(4,1) * t101 - rSges(4,2) * t97) * t98 + (-rSges(4,3) - pkin(8)) * t102 + t123)) - m(5) * (g(1) * (t81 * rSges(5,1) + rSges(5,2) * t116 + t114 * t80 + t109) + g(2) * (t79 * rSges(5,1) + rSges(5,2) * t120 + t114 * t78 + t106) + g(3) * ((rSges(5,1) * t101 + rSges(5,3) * t97) * t98 + (-rSges(5,2) - pkin(8)) * t102 + t111)) - m(6) * (g(1) * (t69 * rSges(6,1) - t68 * rSges(6,2) + t105) + g(2) * (t67 * rSges(6,1) - t66 * rSges(6,2) + t104) + g(3) * (t73 * rSges(6,1) - t72 * rSges(6,2) + (rSges(6,3) - pkin(8)) * t102 + t107) + (-rSges(6,3) - pkin(9)) * t127) - m(7) * (g(1) * (t128 * t68 + t129 * t69 + t105) + g(2) * (t128 * t66 + t129 * t67 + t104) + g(3) * (t129 * t73 + t128 * t72 + (rSges(7,2) - pkin(8)) * t102 + t107) + (-rSges(7,2) - pkin(9)) * t127);
U  = t1;
