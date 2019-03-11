% Calculate potential energy for
% S6RRPRPR13
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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:03
% EndTime: 2019-03-09 11:24:03
% DurationCPUTime: 0.72s
% Computational Cost: add. (272->101), mult. (545->110), div. (0->0), fcn. (622->12), ass. (0->47)
t138 = -m(5) - m(6);
t137 = mrSges(4,2) - mrSges(3,1);
t136 = mrSges(3,3) + mrSges(4,1);
t135 = -mrSges(4,3) + mrSges(3,2);
t134 = m(6) * qJ(5) - m(7) * (-pkin(10) - qJ(5)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t101 = cos(pkin(11));
t99 = sin(pkin(11));
t131 = pkin(5) * t99;
t98 = pkin(11) + qJ(6);
t92 = sin(t98);
t93 = cos(t98);
t133 = -m(7) * (pkin(9) + t131) - t99 * mrSges(6,1) - t92 * mrSges(7,1) - t101 * mrSges(6,2) - t93 * mrSges(7,2) - mrSges(5,3) + t137;
t91 = pkin(5) * t101 + pkin(4);
t132 = -m(6) * pkin(4) - m(7) * t91 - t101 * mrSges(6,1) - t93 * mrSges(7,1) + t99 * mrSges(6,2) + t92 * mrSges(7,2) - mrSges(5,1);
t102 = cos(pkin(6));
t130 = t102 * pkin(8) + pkin(7);
t109 = cos(qJ(1));
t100 = sin(pkin(6));
t106 = sin(qJ(1));
t127 = t100 * t106;
t129 = t109 * pkin(1) + pkin(8) * t127;
t105 = sin(qJ(2));
t128 = t100 * t105;
t108 = cos(qJ(2));
t126 = t100 * t108;
t125 = t100 * t109;
t124 = t105 * t106;
t123 = t105 * t109;
t122 = t106 * t108;
t121 = t108 * t109;
t120 = pkin(2) * t128 + t130;
t119 = pkin(8) * t125;
t117 = t102 * pkin(3) + pkin(9) * t128 + t120;
t81 = -t102 * t121 + t124;
t82 = t102 * t123 + t122;
t96 = t106 * pkin(1);
t116 = t82 * pkin(2) + t81 * qJ(3) + t96;
t83 = t102 * t122 + t123;
t84 = -t102 * t124 + t121;
t115 = t84 * pkin(2) + qJ(3) * t83 + t129;
t114 = pkin(3) * t127 + t115;
t113 = -qJ(3) * t126 + t117;
t111 = (-pkin(3) - pkin(8)) * t125 + t116;
t107 = cos(qJ(4));
t104 = sin(qJ(4));
t80 = t102 * t107 - t104 * t126;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t130 - m(4) * t120 - m(5) * t113 - t80 * mrSges(5,1) - mrSges(5,3) * t128 - m(6) * (pkin(4) * t80 + t113) - (t101 * t80 + t128 * t99) * mrSges(6,1) - (t101 * t128 - t80 * t99) * mrSges(6,2) - m(7) * (t80 * t91 + t117) - (t128 * t92 + t80 * t93) * mrSges(7,1) - (t128 * t93 - t80 * t92) * mrSges(7,2) - t136 * t102 + (((m(4) + m(7)) * qJ(3) - t135) * t108 + (-m(7) * t131 + t137) * t105) * t100 - t134 * (t102 * t104 + t107 * t126)) * g(3) + (-mrSges(1,2) - t106 * mrSges(2,1) - mrSges(2,2) * t109 - m(3) * (t96 - t119) - m(4) * (t116 - t119) - m(7) * t111 + t135 * t81 + t136 * t125 + t138 * (t82 * pkin(9) + t111) + t132 * (t81 * t104 - t107 * t125) + t133 * t82 + t134 * (t104 * t125 + t81 * t107)) * g(2) + (-m(3) * t129 - m(4) * t115 - m(7) * t114 - mrSges(2,1) * t109 + t106 * mrSges(2,2) - mrSges(1,1) + t135 * t83 - t136 * t127 + t138 * (pkin(9) * t84 + t114) + t132 * (t104 * t83 + t107 * t127) + t133 * t84 - t134 * (t104 * t127 - t83 * t107)) * g(1);
U  = t1;
