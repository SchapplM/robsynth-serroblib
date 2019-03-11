% Calculate potential energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:28
% EndTime: 2019-03-09 11:32:28
% DurationCPUTime: 0.64s
% Computational Cost: add. (249->93), mult. (523->95), div. (0->0), fcn. (592->10), ass. (0->52)
t144 = mrSges(4,2) - mrSges(3,1);
t143 = mrSges(3,3) + mrSges(4,1);
t142 = -mrSges(4,3) + mrSges(3,2);
t104 = cos(pkin(6));
t106 = sin(qJ(4));
t110 = cos(qJ(4));
t103 = sin(pkin(6));
t111 = cos(qJ(2));
t128 = t103 * t111;
t86 = t104 * t106 + t110 * t128;
t87 = t104 * t110 - t106 * t128;
t141 = t87 * pkin(4) + qJ(5) * t86;
t140 = m(7) * pkin(10) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t105 = sin(qJ(6));
t109 = cos(qJ(6));
t139 = -t109 * mrSges(7,1) + t105 * mrSges(7,2) - mrSges(6,1) - mrSges(5,3);
t138 = t105 * mrSges(7,1) + t109 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t137 = -m(7) * (pkin(5) + pkin(9)) + t139 + t144;
t112 = cos(qJ(1));
t123 = t111 * t112;
t107 = sin(qJ(2));
t108 = sin(qJ(1));
t126 = t107 * t108;
t91 = -t104 * t126 + t123;
t135 = pkin(9) * t91;
t124 = t108 * t111;
t125 = t107 * t112;
t89 = t104 * t125 + t124;
t134 = t89 * pkin(9);
t133 = t104 * pkin(8) + pkin(7);
t129 = t103 * t108;
t131 = t112 * pkin(1) + pkin(8) * t129;
t130 = t103 * t107;
t127 = t103 * t112;
t122 = pkin(2) * t130 + t133;
t121 = pkin(8) * t127;
t120 = t104 * pkin(3) + pkin(9) * t130 + t122;
t101 = t108 * pkin(1);
t88 = -t104 * t123 + t126;
t119 = t89 * pkin(2) + t88 * qJ(3) + t101;
t90 = t104 * t124 + t125;
t118 = t91 * pkin(2) + qJ(3) * t90 + t131;
t117 = pkin(3) * t129 + t118;
t116 = -qJ(3) * t128 + t120;
t115 = (-pkin(3) - pkin(8)) * t127 + t119;
t77 = t106 * t129 - t90 * t110;
t78 = t106 * t90 + t110 * t129;
t114 = t78 * pkin(4) + t77 * qJ(5) + t117;
t79 = t106 * t127 + t88 * t110;
t80 = -t88 * t106 + t110 * t127;
t113 = -t80 * pkin(4) - t79 * qJ(5) + t115;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t133 - m(4) * t122 - m(5) * t116 - m(6) * (t116 + t141) - m(7) * (t120 + t141) - t138 * t86 + t139 * t130 - t143 * t104 + (((m(4) + m(7)) * qJ(3) - t142) * t111 + (-m(7) * pkin(5) + t144) * t107) * t103 - t140 * t87) * g(3) + (-mrSges(1,2) - t108 * mrSges(2,1) - t112 * mrSges(2,2) - m(3) * (t101 - t121) - m(4) * (t119 - t121) - m(5) * (t115 + t134) - m(6) * (t113 + t134) - m(7) * t113 + t142 * t88 + t143 * t127 + t140 * t80 + t138 * t79 + t137 * t89) * g(2) + (-mrSges(1,1) - t112 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t131 - m(4) * t118 - m(5) * (t117 + t135) - m(6) * (t114 + t135) - m(7) * t114 + t142 * t90 - t143 * t129 - t140 * t78 - t138 * t77 + t137 * t91) * g(1);
U  = t1;
