% Calculate potential energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:30
% EndTime: 2019-03-09 16:03:30
% DurationCPUTime: 0.60s
% Computational Cost: add. (278->91), mult. (599->95), div. (0->0), fcn. (697->10), ass. (0->49)
t100 = sin(qJ(6));
t104 = cos(qJ(6));
t138 = -t100 * mrSges(7,1) - t104 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t105 = cos(qJ(2));
t99 = sin(pkin(6));
t122 = t99 * t105;
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t129 = cos(qJ(3));
t119 = t99 * t129;
t121 = cos(pkin(6));
t84 = t121 * t101 + t102 * t119;
t137 = t84 * pkin(4) + qJ(5) * t122;
t136 = -m(7) * pkin(10) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t134 = -m(7) * (pkin(5) + qJ(4)) - t104 * mrSges(7,1) + t100 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t133 = mrSges(3,2) - t138 + (-m(6) - m(7)) * (pkin(9) - qJ(5));
t106 = cos(qJ(1));
t103 = sin(qJ(1));
t118 = t103 * t121;
t87 = t106 * t102 + t105 * t118;
t132 = pkin(9) * t87;
t117 = t106 * t121;
t85 = t102 * t103 - t105 * t117;
t131 = t85 * pkin(9);
t130 = t121 * pkin(8) + pkin(7);
t125 = t103 * t99;
t126 = t106 * pkin(1) + pkin(8) * t125;
t124 = t106 * t99;
t123 = t99 * t102;
t88 = -t102 * t118 + t106 * t105;
t120 = t88 * pkin(2) + t126;
t78 = t101 * t125 + t88 * t129;
t116 = t78 * pkin(3) + t120;
t115 = t103 * pkin(1) - pkin(8) * t124;
t86 = t102 * t117 + t103 * t105;
t113 = t86 * pkin(2) + t115;
t112 = pkin(2) * t123 - pkin(9) * t122 + t130;
t76 = -t101 * t124 + t86 * t129;
t111 = t76 * pkin(3) + t113;
t77 = t101 * t88 - t103 * t119;
t110 = qJ(4) * t77 + t116;
t109 = t84 * pkin(3) + t112;
t75 = t86 * t101 + t106 * t119;
t108 = t75 * qJ(4) + t111;
t83 = t101 * t123 - t121 * t129;
t107 = qJ(4) * t83 + t109;
t73 = t78 * pkin(4);
t71 = t76 * pkin(4);
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t130 - t121 * mrSges(3,3) - (t102 * mrSges(3,1) + t105 * mrSges(3,2)) * t99 - m(4) * t112 - m(5) * t107 - m(6) * (t107 + t137) - m(7) * (t109 + t137) + t134 * t83 + t138 * t122 + t136 * t84) * g(3) + (-mrSges(1,2) - t103 * mrSges(2,1) - t106 * mrSges(2,2) - m(3) * t115 - t86 * mrSges(3,1) + mrSges(3,3) * t124 - m(4) * (t113 + t131) - m(5) * (t108 + t131) - m(6) * (t108 + t71) - m(7) * (t111 + t71) + t134 * t75 + t133 * t85 + t136 * t76) * g(2) + (-mrSges(1,1) - t106 * mrSges(2,1) + t103 * mrSges(2,2) - m(3) * t126 - t88 * mrSges(3,1) - mrSges(3,3) * t125 - m(4) * (t120 + t132) - m(5) * (t110 + t132) - m(6) * (t110 + t73) - m(7) * (t116 + t73) + t134 * t77 + t133 * t87 + t136 * t78) * g(1);
U  = t1;
