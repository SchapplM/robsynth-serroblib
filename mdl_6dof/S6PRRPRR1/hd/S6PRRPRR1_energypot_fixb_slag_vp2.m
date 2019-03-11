% Calculate potential energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:50:56
% EndTime: 2019-03-08 21:50:57
% DurationCPUTime: 0.87s
% Computational Cost: add. (357->119), mult. (480->138), div. (0->0), fcn. (534->14), ass. (0->47)
t139 = -m(6) - m(7);
t138 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t112 = -qJ(4) - pkin(8);
t137 = m(4) * pkin(8) - m(5) * t112 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t113 = sin(qJ(6));
t116 = cos(qJ(6));
t136 = -m(7) * pkin(5) - t116 * mrSges(7,1) + t113 * mrSges(7,2) - mrSges(6,1);
t135 = -t113 * mrSges(7,1) - t116 * mrSges(7,2) - mrSges(6,3) - t137;
t114 = sin(qJ(3));
t134 = pkin(3) * t114;
t117 = cos(qJ(3));
t98 = t117 * pkin(3) + pkin(2);
t110 = cos(pkin(11));
t108 = sin(pkin(11));
t109 = sin(pkin(6));
t132 = t108 * t109;
t133 = t110 * pkin(1) + pkin(7) * t132;
t131 = t109 * t110;
t130 = t109 * t114;
t115 = sin(qJ(2));
t129 = t109 * t115;
t128 = t109 * t117;
t118 = cos(qJ(2));
t127 = t109 * t118;
t111 = cos(pkin(6));
t126 = t111 * t115;
t125 = t111 * t118;
t124 = t111 * pkin(7) + qJ(1);
t107 = qJ(3) + pkin(12);
t123 = t108 * t130;
t102 = t108 * pkin(1);
t122 = -pkin(7) * t131 + t102;
t106 = -pkin(9) + t112;
t100 = cos(t107);
t90 = pkin(4) * t100 + t98;
t99 = sin(t107);
t91 = pkin(4) * t99 + t134;
t120 = t106 * t127 + t111 * t91 + t90 * t129 + t124;
t101 = qJ(5) + t107;
t97 = cos(t101);
t96 = sin(t101);
t87 = -t108 * t126 + t110 * t118;
t86 = t108 * t125 + t110 * t115;
t85 = t108 * t118 + t110 * t126;
t84 = t108 * t115 - t110 * t125;
t79 = t111 * t96 + t129 * t97;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(6) * t120 - t79 * mrSges(6,1) + mrSges(6,3) * t127 - m(7) * (pkin(5) * t79 + t120) - (-t113 * t127 + t79 * t116) * mrSges(7,1) - (-t79 * t113 - t116 * t127) * mrSges(7,2) + t138 * (-t111 * t97 + t129 * t96) + (-m(3) - m(4) - m(5)) * t124 + (-m(5) * t134 - t114 * mrSges(4,1) - t99 * mrSges(5,1) - t117 * mrSges(4,2) - t100 * mrSges(5,2) - mrSges(3,3)) * t111 + (t137 * t118 + (-m(4) * pkin(2) - m(5) * t98 - t117 * mrSges(4,1) - t100 * mrSges(5,1) + t114 * mrSges(4,2) + t99 * mrSges(5,2) - mrSges(3,1)) * t115) * t109) * g(3) + (-mrSges(1,2) - t108 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t122 - t85 * mrSges(3,1) + mrSges(3,3) * t131 - m(4) * (pkin(2) * t85 + t122) - (-t110 * t130 + t117 * t85) * mrSges(4,1) - (-t110 * t128 - t114 * t85) * mrSges(4,2) - m(5) * (t85 * t98 + t102 + (-pkin(7) - t134) * t131) - (t100 * t85 - t131 * t99) * mrSges(5,1) - (-t100 * t131 - t85 * t99) * mrSges(5,2) + t139 * (t102 + t85 * t90 - t84 * t106 + (-pkin(7) - t91) * t131) + t138 * (t131 * t97 + t85 * t96) + t136 * (-t131 * t96 + t85 * t97) + t135 * t84) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t133 - t87 * mrSges(3,1) - mrSges(3,3) * t132 - m(4) * (pkin(2) * t87 + t133) - (t117 * t87 + t123) * mrSges(4,1) - (t108 * t128 - t114 * t87) * mrSges(4,2) - m(5) * (pkin(3) * t123 + t87 * t98 + t133) - (t100 * t87 + t132 * t99) * mrSges(5,1) - (t100 * t132 - t87 * t99) * mrSges(5,2) + t139 * (-t86 * t106 + t91 * t132 + t87 * t90 + t133) + t138 * (-t132 * t97 + t87 * t96) + t136 * (t132 * t96 + t87 * t97) + t135 * t86) * g(1);
U  = t1;
