% Calculate potential energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:05
% EndTime: 2019-03-08 19:51:06
% DurationCPUTime: 0.65s
% Computational Cost: add. (249->93), mult. (523->99), div. (0->0), fcn. (592->10), ass. (0->52)
t143 = mrSges(4,2) - mrSges(3,1);
t142 = mrSges(3,3) + mrSges(4,1);
t141 = -mrSges(4,3) + mrSges(3,2);
t105 = cos(pkin(6));
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t103 = sin(pkin(6));
t111 = cos(qJ(2));
t125 = t103 * t111;
t89 = t105 * t107 + t110 * t125;
t90 = t105 * t110 - t107 * t125;
t140 = t90 * pkin(4) + t89 * qJ(5);
t139 = m(7) * pkin(9) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t138 = -t109 * mrSges(7,1) + t106 * mrSges(7,2) - mrSges(6,1) - mrSges(5,3);
t137 = t106 * mrSges(7,1) + t109 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t136 = -m(7) * (pkin(5) + pkin(8)) + t138 + t143;
t102 = sin(pkin(10));
t104 = cos(pkin(10));
t108 = sin(qJ(2));
t124 = t105 * t108;
t86 = t102 * t111 + t104 * t124;
t134 = pkin(8) * t86;
t88 = -t102 * t124 + t104 * t111;
t133 = pkin(8) * t88;
t130 = t102 * t103;
t132 = t104 * pkin(1) + pkin(7) * t130;
t129 = t103 * t104;
t128 = t103 * t107;
t127 = t103 * t108;
t126 = t103 * t110;
t123 = t105 * t111;
t122 = t105 * pkin(7) + qJ(1);
t121 = pkin(7) * t129;
t120 = pkin(2) * t127 + t122;
t85 = t102 * t108 - t104 * t123;
t98 = t102 * pkin(1);
t119 = t86 * pkin(2) + qJ(3) * t85 + t98;
t118 = t105 * pkin(3) + pkin(8) * t127 + t120;
t87 = t102 * t123 + t104 * t108;
t117 = t88 * pkin(2) + qJ(3) * t87 + t132;
t116 = pkin(3) * t130 + t117;
t115 = -qJ(3) * t125 + t118;
t114 = (-pkin(3) - pkin(7)) * t129 + t119;
t76 = t102 * t128 - t87 * t110;
t77 = t102 * t126 + t107 * t87;
t113 = t77 * pkin(4) + t76 * qJ(5) + t116;
t78 = t104 * t128 + t110 * t85;
t79 = t104 * t126 - t85 * t107;
t112 = -t79 * pkin(4) - t78 * qJ(5) + t114;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t122 - m(4) * t120 - m(5) * t115 - m(6) * (t115 + t140) - m(7) * (t118 + t140) - t137 * t89 + t138 * t127 - t142 * t105 + (((m(4) + m(7)) * qJ(3) - t141) * t111 + (-m(7) * pkin(5) + t143) * t108) * t103 - t139 * t90) * g(3) + (-mrSges(1,2) - t102 * mrSges(2,1) - t104 * mrSges(2,2) - m(3) * (t98 - t121) - m(4) * (t119 - t121) - m(5) * (t114 + t134) - m(6) * (t112 + t134) - m(7) * t112 + t141 * t85 + t142 * t129 + t139 * t79 + t137 * t78 + t136 * t86) * g(2) + (-mrSges(1,1) - t104 * mrSges(2,1) + t102 * mrSges(2,2) - m(3) * t132 - m(4) * t117 - m(5) * (t116 + t133) - m(6) * (t113 + t133) - m(7) * t113 + t141 * t87 - t142 * t130 - t139 * t77 - t137 * t76 + t136 * t88) * g(1);
U  = t1;
