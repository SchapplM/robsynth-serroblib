% Calculate potential energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:25
% EndTime: 2019-03-09 22:26:26
% DurationCPUTime: 0.86s
% Computational Cost: add. (357->119), mult. (480->134), div. (0->0), fcn. (534->14), ass. (0->47)
t139 = -m(6) - m(7);
t138 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t118 = -pkin(10) - pkin(9);
t137 = m(4) * pkin(9) - m(5) * t118 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t136 = -m(7) * pkin(5) - t114 * mrSges(7,1) + t110 * mrSges(7,2) - mrSges(6,1);
t135 = -t110 * mrSges(7,1) - t114 * mrSges(7,2) - mrSges(6,3) - t137;
t111 = sin(qJ(3));
t134 = pkin(3) * t111;
t109 = cos(pkin(6));
t133 = t109 * pkin(8) + pkin(7);
t115 = cos(qJ(3));
t98 = t115 * pkin(3) + pkin(2);
t117 = cos(qJ(1));
t108 = sin(pkin(6));
t113 = sin(qJ(1));
t130 = t108 * t113;
t132 = t117 * pkin(1) + pkin(8) * t130;
t112 = sin(qJ(2));
t131 = t108 * t112;
t116 = cos(qJ(2));
t129 = t108 * t116;
t128 = t108 * t117;
t127 = t112 * t113;
t126 = t112 * t117;
t125 = t113 * t116;
t124 = t116 * t117;
t107 = qJ(3) + qJ(4);
t123 = t111 * t130;
t106 = -qJ(5) + t118;
t101 = cos(t107);
t90 = pkin(4) * t101 + t98;
t100 = sin(t107);
t91 = pkin(4) * t100 + t134;
t122 = t106 * t129 + t109 * t91 + t90 * t131 + t133;
t103 = t113 * pkin(1);
t121 = -pkin(8) * t128 + t103;
t99 = pkin(12) + t107;
t97 = cos(t99);
t96 = sin(t99);
t87 = -t109 * t127 + t124;
t86 = t109 * t125 + t126;
t85 = t109 * t126 + t125;
t84 = -t109 * t124 + t127;
t79 = t109 * t96 + t131 * t97;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(6) * t122 - t79 * mrSges(6,1) + mrSges(6,3) * t129 - m(7) * (pkin(5) * t79 + t122) - (-t110 * t129 + t114 * t79) * mrSges(7,1) - (-t110 * t79 - t114 * t129) * mrSges(7,2) + t138 * (-t109 * t97 + t131 * t96) + (-m(3) - m(4) - m(5)) * t133 + (-m(5) * t134 - t111 * mrSges(4,1) - t100 * mrSges(5,1) - t115 * mrSges(4,2) - t101 * mrSges(5,2) - mrSges(3,3)) * t109 + (t137 * t116 + (-m(4) * pkin(2) - m(5) * t98 - t115 * mrSges(4,1) - t101 * mrSges(5,1) + t111 * mrSges(4,2) + t100 * mrSges(5,2) - mrSges(3,1)) * t112) * t108) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t113 - mrSges(2,2) * t117 - m(3) * t121 - t85 * mrSges(3,1) + mrSges(3,3) * t128 - m(4) * (pkin(2) * t85 + t121) - (-t111 * t128 + t115 * t85) * mrSges(4,1) - (-t111 * t85 - t115 * t128) * mrSges(4,2) - m(5) * (t85 * t98 + t103 + (-pkin(8) - t134) * t128) - (-t100 * t128 + t101 * t85) * mrSges(5,1) - (-t100 * t85 - t101 * t128) * mrSges(5,2) + t139 * (t103 + t85 * t90 - t84 * t106 + (-pkin(8) - t91) * t128) + t138 * (t128 * t97 + t85 * t96) + t136 * (-t128 * t96 + t85 * t97) + t135 * t84) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t117 + mrSges(2,2) * t113 - m(3) * t132 - t87 * mrSges(3,1) - mrSges(3,3) * t130 - m(4) * (pkin(2) * t87 + t132) - (t115 * t87 + t123) * mrSges(4,1) - (-t111 * t87 + t115 * t130) * mrSges(4,2) - m(5) * (pkin(3) * t123 + t87 * t98 + t132) - (t100 * t130 + t101 * t87) * mrSges(5,1) - (-t100 * t87 + t101 * t130) * mrSges(5,2) + t139 * (-t86 * t106 + t91 * t130 + t87 * t90 + t132) + t138 * (-t130 * t97 + t87 * t96) + t136 * (t130 * t96 + t87 * t97) + t135 * t86) * g(1);
U  = t1;
