% Calculate potential energy for
% S6RRPPRR9
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:24
% EndTime: 2019-03-09 09:28:25
% DurationCPUTime: 0.83s
% Computational Cost: add. (228->95), mult. (468->101), div. (0->0), fcn. (516->10), ass. (0->48)
t143 = pkin(9) - qJ(3);
t142 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t141 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3);
t140 = -mrSges(3,3) - mrSges(4,1) - mrSges(5,1);
t139 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t138 = -m(7) * pkin(5) - t107 * mrSges(7,1) + t103 * mrSges(7,2) - mrSges(6,1);
t137 = m(6) * t143 + t103 * mrSges(7,1) + t107 * mrSges(7,2) + mrSges(6,3) - t141;
t136 = -pkin(3) - pkin(8);
t102 = cos(pkin(6));
t135 = t102 * pkin(8) + pkin(7);
t106 = sin(qJ(1));
t109 = cos(qJ(2));
t123 = t106 * t109;
t105 = sin(qJ(2));
t110 = cos(qJ(1));
t124 = t105 * t110;
t85 = t102 * t123 + t124;
t133 = qJ(3) * t85;
t122 = t109 * t110;
t125 = t105 * t106;
t83 = -t102 * t122 + t125;
t132 = t83 * qJ(3);
t101 = sin(pkin(6));
t129 = t101 * t106;
t131 = t110 * pkin(1) + pkin(8) * t129;
t130 = t101 * t105;
t108 = cos(qJ(5));
t128 = t101 * t108;
t127 = t101 * t109;
t126 = t101 * t110;
t121 = pkin(2) * t130 + t135;
t86 = -t102 * t125 + t122;
t120 = t86 * pkin(2) + t131;
t119 = t102 * pkin(3) + qJ(4) * t130 + t121;
t99 = t106 * pkin(1);
t118 = -pkin(8) * t126 + t99;
t84 = t102 * t124 + t123;
t79 = t84 * pkin(2);
t117 = t84 * qJ(4) + t79 + t99;
t115 = (-pkin(4) + t136) * t126;
t114 = pkin(3) * t129 + qJ(4) * t86 + t120;
t113 = t117 + t132;
t111 = t102 * pkin(4) + t143 * t127 + t119;
t104 = sin(qJ(5));
t82 = t102 * t108 + t104 * t130;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t135 - m(4) * t121 - m(5) * t119 - m(6) * t111 - t82 * mrSges(6,1) - mrSges(6,3) * t127 - m(7) * (pkin(5) * t82 + t111) - (t103 * t127 + t107 * t82) * mrSges(7,1) - (-t103 * t82 + t107 * t127) * mrSges(7,2) - t139 * (t102 * t104 - t105 * t128) + t140 * t102 + (((m(4) + m(5)) * qJ(3) + t141) * t109 + t142 * t105) * t101) * g(3) + (-mrSges(1,2) - t106 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t118 - m(4) * (t118 + t79 + t132) - m(5) * t113 - (t115 + t117) * m(6) - (t113 + t115) * m(7) + t139 * (t104 * t126 + t84 * t108) + t142 * t84 + t138 * (t84 * t104 - t108 * t126) + (-m(5) * t136 - t140) * t126 + (m(7) * pkin(9) + t137) * t83) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t106 * mrSges(2,2) - m(3) * t131 - m(4) * (t120 + t133) - m(5) * (t114 + t133) + (-m(6) - m(7)) * (pkin(4) * t129 + t114) - t139 * (t104 * t129 - t86 * t108) + t142 * t86 + t138 * (t104 * t86 + t106 * t128) + t140 * t129 + (m(7) * t143 + t137) * t85) * g(1);
U  = t1;
