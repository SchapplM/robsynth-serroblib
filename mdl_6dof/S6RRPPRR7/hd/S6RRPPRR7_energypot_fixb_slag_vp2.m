% Calculate potential energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:04
% EndTime: 2019-03-09 09:17:05
% DurationCPUTime: 0.67s
% Computational Cost: add. (232->85), mult. (478->97), div. (0->0), fcn. (530->10), ass. (0->41)
t133 = -m(6) - m(7);
t132 = -mrSges(4,3) - mrSges(5,1) + mrSges(3,2);
t131 = -mrSges(5,3) + mrSges(4,2) + mrSges(3,3);
t130 = -mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t129 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t102 = cos(qJ(6));
t98 = sin(qJ(6));
t128 = -m(7) * pkin(5) - t102 * mrSges(7,1) + t98 * mrSges(7,2) - mrSges(6,1);
t127 = -t98 * mrSges(7,1) - t102 * mrSges(7,2) + t130;
t97 = cos(pkin(6));
t126 = t97 * pkin(8) + pkin(7);
t105 = cos(qJ(1));
t101 = sin(qJ(1));
t96 = sin(pkin(6));
t123 = t101 * t96;
t125 = t105 * pkin(1) + pkin(8) * t123;
t100 = sin(qJ(2));
t124 = t100 * t96;
t104 = cos(qJ(2));
t122 = t104 * t96;
t121 = t105 * t96;
t120 = t100 * t101;
t119 = t100 * t105;
t118 = t101 * t104;
t117 = t104 * t105;
t116 = pkin(2) * t124 + t126;
t114 = t101 * pkin(1) - pkin(8) * t121;
t83 = t118 * t97 + t119;
t84 = -t120 * t97 + t117;
t113 = t84 * pkin(2) + qJ(3) * t83 + t125;
t112 = pkin(3) * t124 - qJ(4) * t97 + t116;
t111 = pkin(9) * t124 + t112;
t81 = -t117 * t97 + t120;
t82 = t119 * t97 + t118;
t110 = t82 * pkin(2) + t81 * qJ(3) + t114;
t109 = t82 * pkin(3) + qJ(4) * t121 + t110;
t108 = t84 * pkin(3) - qJ(4) * t123 + t113;
t103 = cos(qJ(5));
t99 = sin(qJ(5));
t80 = -t103 * t122 - t97 * t99;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t126 - m(4) * t116 - m(5) * t112 - m(6) * t111 - t80 * mrSges(6,1) - m(7) * (pkin(5) * t80 + t111) - (t102 * t80 + t124 * t98) * mrSges(7,1) - (t102 * t124 - t80 * t98) * mrSges(7,2) + t129 * (-t103 * t97 + t122 * t99) - t131 * t97 + (t130 * t100 + (t133 * (-pkin(4) - qJ(3)) + (m(4) + m(5)) * qJ(3) - t132) * t104) * t96) * g(3) + (-m(3) * t114 - m(4) * t110 - m(5) * t109 - t101 * mrSges(2,1) - mrSges(2,2) * t105 - mrSges(1,2) + t133 * (t81 * pkin(4) + t82 * pkin(9) + t109) - t129 * (-t103 * t121 + t81 * t99) + t132 * t81 + t128 * (t81 * t103 + t121 * t99) + t131 * t121 + t127 * t82) * g(2) + (-m(3) * t125 - m(4) * t113 - m(5) * t108 - mrSges(2,1) * t105 + t101 * mrSges(2,2) - mrSges(1,1) + t133 * (t83 * pkin(4) + pkin(9) * t84 + t108) - t129 * (t103 * t123 + t83 * t99) + t132 * t83 + t128 * (t103 * t83 - t123 * t99) - t131 * t123 + t127 * t84) * g(1);
U  = t1;
