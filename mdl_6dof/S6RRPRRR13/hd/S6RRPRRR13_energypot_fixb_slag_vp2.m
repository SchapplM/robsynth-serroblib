% Calculate potential energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:51
% EndTime: 2019-03-09 14:44:51
% DurationCPUTime: 0.72s
% Computational Cost: add. (272->91), mult. (545->96), div. (0->0), fcn. (622->12), ass. (0->44)
t140 = -m(5) - m(6);
t139 = -mrSges(3,1) + mrSges(4,2);
t138 = mrSges(3,3) + mrSges(4,1);
t137 = -mrSges(4,3) + mrSges(3,2);
t136 = m(6) * pkin(10) - m(7) * (-pkin(11) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t106 = cos(qJ(5));
t99 = qJ(5) + qJ(6);
t93 = sin(t99);
t94 = cos(t99);
t135 = -t93 * mrSges(7,1) - t106 * mrSges(6,2) - t94 * mrSges(7,2) - mrSges(5,3);
t102 = sin(qJ(5));
t134 = -m(7) * (pkin(5) * t102 + pkin(9)) - t102 * mrSges(6,1) + t135 + t139;
t133 = -m(6) * pkin(4) - m(7) * (pkin(5) * t106 + pkin(4)) - t106 * mrSges(6,1) - t94 * mrSges(7,1) + t102 * mrSges(6,2) + t93 * mrSges(7,2) - mrSges(5,1);
t101 = cos(pkin(6));
t132 = t101 * pkin(8) + pkin(7);
t109 = cos(qJ(1));
t100 = sin(pkin(6));
t105 = sin(qJ(1));
t129 = t100 * t105;
t131 = t109 * pkin(1) + pkin(8) * t129;
t104 = sin(qJ(2));
t130 = t100 * t104;
t108 = cos(qJ(2));
t128 = t100 * t108;
t127 = t100 * t109;
t125 = t104 * t105;
t124 = t104 * t109;
t123 = t105 * t108;
t122 = t108 * t109;
t121 = pkin(2) * t130 + t132;
t120 = pkin(8) * t127;
t118 = t101 * pkin(3) + pkin(9) * t130 + t121;
t82 = -t101 * t122 + t125;
t83 = t101 * t124 + t123;
t97 = t105 * pkin(1);
t117 = t83 * pkin(2) + t82 * qJ(3) + t97;
t84 = t101 * t123 + t124;
t85 = -t101 * t125 + t122;
t116 = t85 * pkin(2) + t84 * qJ(3) + t131;
t115 = pkin(3) * t129 + t116;
t112 = (-pkin(3) - pkin(8)) * t127 + t117;
t107 = cos(qJ(4));
t103 = sin(qJ(4));
t1 = (-m(2) * pkin(7) - m(3) * t132 - m(4) * t121 - m(7) * t118 - mrSges(1,3) - mrSges(2,3) + t140 * (-qJ(3) * t128 + t118) + t135 * t130 + t133 * (t101 * t107 - t103 * t128) - t138 * t101 + (((m(4) + m(7)) * qJ(3) - t137) * t108 + ((-m(7) * pkin(5) - mrSges(6,1)) * t102 + t139) * t104) * t100 - t136 * (t101 * t103 + t107 * t128)) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t105 - mrSges(2,2) * t109 - m(3) * (t97 - t120) - m(4) * (t117 - t120) - m(7) * t112 + t137 * t82 + t138 * t127 + t140 * (t83 * pkin(9) + t112) + t133 * (t103 * t82 - t107 * t127) + t134 * t83 + t136 * (t103 * t127 + t107 * t82)) * g(2) + (-m(3) * t131 - m(4) * t116 - m(7) * t115 - mrSges(2,1) * t109 + mrSges(2,2) * t105 - mrSges(1,1) + t137 * t84 - t138 * t129 + t140 * (t85 * pkin(9) + t115) + t133 * (t103 * t84 + t107 * t129) + t134 * t85 - t136 * (t103 * t129 - t84 * t107)) * g(1);
U  = t1;
