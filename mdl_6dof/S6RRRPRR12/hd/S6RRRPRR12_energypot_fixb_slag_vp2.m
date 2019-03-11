% Calculate potential energy for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:04
% EndTime: 2019-03-09 19:38:05
% DurationCPUTime: 0.73s
% Computational Cost: add. (335->105), mult. (613->121), div. (0->0), fcn. (717->14), ass. (0->46)
t134 = -m(4) - m(5);
t133 = -m(6) - m(7);
t100 = sin(pkin(12));
t117 = pkin(4) * t100 + pkin(9);
t103 = -pkin(10) - qJ(4);
t132 = -m(5) * qJ(4) + m(6) * t103 + m(7) * (-pkin(11) + t103) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t102 = cos(pkin(12));
t99 = pkin(12) + qJ(5);
t91 = sin(t99);
t129 = pkin(5) * t91 + t117;
t93 = qJ(6) + t99;
t88 = sin(t93);
t89 = cos(t93);
t92 = cos(t99);
t131 = -m(6) * t117 - m(7) * t129 - t100 * mrSges(5,1) - t91 * mrSges(6,1) - t88 * mrSges(7,1) - t102 * mrSges(5,2) - t92 * mrSges(6,2) - t89 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t90 = pkin(4) * t102 + pkin(3);
t81 = pkin(5) * t92 + t90;
t130 = -m(5) * pkin(3) - m(6) * t90 - m(7) * t81 - mrSges(5,1) * t102 - mrSges(6,1) * t92 - mrSges(7,1) * t89 + mrSges(5,2) * t100 + mrSges(6,2) * t91 + mrSges(7,2) * t88 - mrSges(4,1);
t124 = cos(pkin(6));
t128 = pkin(8) * t124 + pkin(7);
t127 = cos(qJ(3));
t108 = cos(qJ(1));
t101 = sin(pkin(6));
t106 = sin(qJ(1));
t122 = t101 * t106;
t125 = pkin(1) * t108 + pkin(8) * t122;
t105 = sin(qJ(2));
t123 = t101 * t105;
t107 = cos(qJ(2));
t121 = t101 * t107;
t120 = t108 * t101;
t119 = pkin(2) * t123 + t128;
t115 = t106 * t124;
t80 = -t105 * t115 + t107 * t108;
t118 = pkin(2) * t80 + t125;
t116 = t101 * t127;
t114 = t108 * t124;
t113 = pkin(1) * t106 - pkin(8) * t120;
t78 = t105 * t114 + t106 * t107;
t111 = pkin(2) * t78 + t113;
t110 = -pkin(9) * t121 + t119;
t104 = sin(qJ(3));
t79 = t105 * t108 + t107 * t115;
t77 = t105 * t106 - t107 * t114;
t76 = t104 * t124 + t105 * t116;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t128 - t124 * mrSges(3,3) - (t105 * mrSges(3,1) + t107 * mrSges(3,2)) * t101 - m(4) * t110 - t76 * mrSges(4,1) + mrSges(4,3) * t121 - m(5) * (pkin(3) * t76 + t110) - (-t100 * t121 + t102 * t76) * mrSges(5,1) - (-t100 * t76 - t102 * t121) * mrSges(5,2) - m(6) * (-t117 * t121 + t76 * t90 + t119) - (-t121 * t91 + t76 * t92) * mrSges(6,1) - (-t121 * t92 - t76 * t91) * mrSges(6,2) - m(7) * (-t121 * t129 + t76 * t81 + t119) - (-t121 * t88 + t76 * t89) * mrSges(7,1) - (-t121 * t89 - t76 * t88) * mrSges(7,2) + t132 * (t104 * t123 - t124 * t127)) * g(3) + (-m(3) * t113 - t106 * mrSges(2,1) - t78 * mrSges(3,1) - t108 * mrSges(2,2) + mrSges(3,3) * t120 - mrSges(1,2) + t133 * t111 + t134 * (pkin(9) * t77 + t111) + t130 * (-t104 * t120 + t127 * t78) + t131 * t77 + t132 * (t104 * t78 + t108 * t116)) * g(2) + (-m(3) * t125 - t108 * mrSges(2,1) - t80 * mrSges(3,1) + t106 * mrSges(2,2) - mrSges(3,3) * t122 - mrSges(1,1) + t133 * t118 + t134 * (pkin(9) * t79 + t118) + t130 * (t104 * t122 + t127 * t80) + t131 * t79 + t132 * (t104 * t80 - t106 * t116)) * g(1);
U  = t1;
