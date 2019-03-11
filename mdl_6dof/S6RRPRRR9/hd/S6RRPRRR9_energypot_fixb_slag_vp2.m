% Calculate potential energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:14
% EndTime: 2019-03-09 14:08:14
% DurationCPUTime: 0.87s
% Computational Cost: add. (357->119), mult. (480->134), div. (0->0), fcn. (534->14), ass. (0->47)
t139 = -m(6) - m(7);
t138 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t112 = -pkin(9) - qJ(3);
t137 = m(4) * qJ(3) - m(5) * t112 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t113 = sin(qJ(6));
t116 = cos(qJ(6));
t136 = -m(7) * pkin(5) - t116 * mrSges(7,1) + t113 * mrSges(7,2) - mrSges(6,1);
t135 = -t113 * mrSges(7,1) - t116 * mrSges(7,2) - mrSges(6,3) - t137;
t108 = sin(pkin(12));
t134 = pkin(3) * t108;
t110 = cos(pkin(12));
t98 = t110 * pkin(3) + pkin(2);
t111 = cos(pkin(6));
t133 = t111 * pkin(8) + pkin(7);
t118 = cos(qJ(1));
t109 = sin(pkin(6));
t115 = sin(qJ(1));
t130 = t109 * t115;
t132 = t118 * pkin(1) + pkin(8) * t130;
t114 = sin(qJ(2));
t131 = t109 * t114;
t117 = cos(qJ(2));
t129 = t109 * t117;
t128 = t109 * t118;
t127 = t114 * t115;
t126 = t114 * t118;
t125 = t115 * t117;
t124 = t117 * t118;
t107 = pkin(12) + qJ(4);
t123 = t108 * t130;
t106 = -pkin(10) + t112;
t100 = cos(t107);
t89 = pkin(4) * t100 + t98;
t99 = sin(t107);
t91 = pkin(4) * t99 + t134;
t122 = t106 * t129 + t111 * t91 + t89 * t131 + t133;
t104 = t115 * pkin(1);
t121 = -pkin(8) * t128 + t104;
t101 = qJ(5) + t107;
t97 = cos(t101);
t96 = sin(t101);
t87 = -t111 * t127 + t124;
t86 = t111 * t125 + t126;
t85 = t111 * t126 + t125;
t84 = -t111 * t124 + t127;
t79 = t111 * t96 + t131 * t97;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(6) * t122 - t79 * mrSges(6,1) + mrSges(6,3) * t129 - m(7) * (pkin(5) * t79 + t122) - (-t113 * t129 + t116 * t79) * mrSges(7,1) - (-t113 * t79 - t116 * t129) * mrSges(7,2) + t138 * (-t111 * t97 + t131 * t96) + (-m(3) - m(4) - m(5)) * t133 + (-m(5) * t134 - t108 * mrSges(4,1) - t99 * mrSges(5,1) - t110 * mrSges(4,2) - t100 * mrSges(5,2) - mrSges(3,3)) * t111 + (t137 * t117 + (-m(4) * pkin(2) - m(5) * t98 - t110 * mrSges(4,1) - t100 * mrSges(5,1) + t108 * mrSges(4,2) + t99 * mrSges(5,2) - mrSges(3,1)) * t114) * t109) * g(3) + (-mrSges(1,2) - t115 * mrSges(2,1) - t118 * mrSges(2,2) - m(3) * t121 - t85 * mrSges(3,1) + mrSges(3,3) * t128 - m(4) * (t85 * pkin(2) + t121) - (-t108 * t128 + t85 * t110) * mrSges(4,1) - (-t85 * t108 - t110 * t128) * mrSges(4,2) - m(5) * (t85 * t98 + t104 + (-pkin(8) - t134) * t128) - (t85 * t100 - t128 * t99) * mrSges(5,1) - (-t100 * t128 - t85 * t99) * mrSges(5,2) + t139 * (t104 + t85 * t89 - t84 * t106 + (-pkin(8) - t91) * t128) + t138 * (t128 * t97 + t85 * t96) + t136 * (-t128 * t96 + t85 * t97) + t135 * t84) * g(2) + (-mrSges(1,1) - t118 * mrSges(2,1) + t115 * mrSges(2,2) - m(3) * t132 - t87 * mrSges(3,1) - mrSges(3,3) * t130 - m(4) * (pkin(2) * t87 + t132) - (t110 * t87 + t123) * mrSges(4,1) - (-t108 * t87 + t110 * t130) * mrSges(4,2) - m(5) * (pkin(3) * t123 + t87 * t98 + t132) - (t100 * t87 + t130 * t99) * mrSges(5,1) - (t100 * t130 - t87 * t99) * mrSges(5,2) + t139 * (-t86 * t106 + t91 * t130 + t87 * t89 + t132) + t138 * (-t130 * t97 + t87 * t96) + t136 * (t130 * t96 + t87 * t97) + t135 * t86) * g(1);
U  = t1;
