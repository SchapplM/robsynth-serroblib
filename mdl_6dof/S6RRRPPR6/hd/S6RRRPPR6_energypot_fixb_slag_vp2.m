% Calculate potential energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:40
% EndTime: 2019-03-09 15:47:41
% DurationCPUTime: 0.70s
% Computational Cost: add. (332->97), mult. (532->108), div. (0->0), fcn. (603->12), ass. (0->49)
t141 = -m(6) - m(7);
t140 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3);
t139 = -m(7) * pkin(10) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t109 = sin(qJ(6));
t113 = cos(qJ(6));
t138 = -t109 * mrSges(7,1) - t113 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t137 = m(7) * pkin(5) + t113 * mrSges(7,1) - t109 * mrSges(7,2) + mrSges(6,1) + mrSges(5,3);
t136 = -t137 - t140;
t110 = sin(qJ(3));
t135 = pkin(3) * t110;
t107 = cos(pkin(6));
t134 = t107 * pkin(8) + pkin(7);
t116 = cos(qJ(1));
t106 = sin(pkin(6));
t112 = sin(qJ(1));
t131 = t106 * t112;
t133 = t116 * pkin(1) + pkin(8) * t131;
t111 = sin(qJ(2));
t132 = t106 * t111;
t115 = cos(qJ(2));
t130 = t106 * t115;
t129 = t106 * t116;
t128 = t111 * t112;
t127 = t111 * t116;
t126 = t112 * t115;
t125 = t115 * t116;
t124 = t110 * t131;
t108 = -qJ(4) - pkin(9);
t114 = cos(qJ(3));
t99 = pkin(3) * t114 + pkin(2);
t123 = t107 * t135 + t108 * t130 + t99 * t132 + t134;
t103 = t112 * pkin(1);
t122 = -pkin(8) * t129 + t103;
t88 = t107 * t126 + t127;
t89 = -t107 * t128 + t125;
t121 = pkin(3) * t124 - t88 * t108 + t89 * t99 + t133;
t86 = -t107 * t125 + t128;
t87 = t107 * t127 + t126;
t118 = t103 + t87 * t99 - t86 * t108 + (-pkin(8) - t135) * t129;
t105 = qJ(3) + pkin(11);
t101 = cos(t105);
t100 = sin(t105);
t83 = t100 * t107 + t101 * t132;
t82 = t100 * t132 - t107 * t101;
t78 = t100 * t131 + t101 * t89;
t77 = t100 * t89 - t101 * t131;
t76 = -t100 * t129 + t101 * t87;
t75 = t87 * t100 + t101 * t129;
t1 = (-m(2) * pkin(7) - m(5) * t123 - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4)) * t134 + t141 * (t83 * pkin(4) + qJ(5) * t82 + t123) + t138 * t82 + t137 * t130 + (-t110 * mrSges(4,1) - t114 * mrSges(4,2) - mrSges(3,3)) * t107 + (t140 * t115 + (-m(4) * pkin(2) - t114 * mrSges(4,1) + t110 * mrSges(4,2) - mrSges(3,1)) * t111) * t106 + t139 * t83) * g(3) + (-mrSges(1,2) - t112 * mrSges(2,1) - t116 * mrSges(2,2) - m(3) * t122 - t87 * mrSges(3,1) + mrSges(3,3) * t129 - m(4) * (t87 * pkin(2) + t122) - (-t110 * t129 + t87 * t114) * mrSges(4,1) - (-t87 * t110 - t114 * t129) * mrSges(4,2) - m(5) * t118 + t141 * (t76 * pkin(4) + t75 * qJ(5) + t118) + t139 * t76 + t138 * t75 + t136 * t86) * g(2) + (-mrSges(1,1) - t116 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t133 - t89 * mrSges(3,1) - mrSges(3,3) * t131 - m(4) * (pkin(2) * t89 + t133) - (t114 * t89 + t124) * mrSges(4,1) - (-t110 * t89 + t114 * t131) * mrSges(4,2) - m(5) * t121 + t141 * (t78 * pkin(4) + qJ(5) * t77 + t121) + t139 * t78 + t138 * t77 + t136 * t88) * g(1);
U  = t1;
