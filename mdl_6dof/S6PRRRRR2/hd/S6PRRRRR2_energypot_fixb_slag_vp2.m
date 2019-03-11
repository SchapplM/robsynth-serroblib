% Calculate potential energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:20
% EndTime: 2019-03-09 00:42:21
% DurationCPUTime: 0.92s
% Computational Cost: add. (355->108), mult. (554->127), div. (0->0), fcn. (633->14), ass. (0->44)
t106 = qJ(5) + qJ(6);
t101 = cos(t106);
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t119 = -pkin(9) - pkin(8);
t143 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t144 = -m(5) - m(6) - m(7);
t99 = sin(t106);
t145 = -t99 * mrSges(7,1) - t115 * mrSges(6,2) - t101 * mrSges(7,2) - t144 * t119 - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t112 - t143;
t141 = -m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t97 = pkin(5) * t115 + pkin(4);
t139 = -m(6) * pkin(4) - m(7) * t97 - t115 * mrSges(6,1) - t101 * mrSges(7,1) + t112 * mrSges(6,2) + t99 * mrSges(7,2) - mrSges(5,1);
t113 = sin(qJ(3));
t138 = pkin(3) * t113;
t110 = cos(pkin(12));
t108 = sin(pkin(12));
t109 = sin(pkin(6));
t134 = t108 * t109;
t135 = t110 * pkin(1) + pkin(7) * t134;
t133 = t109 * t110;
t132 = t109 * t113;
t114 = sin(qJ(2));
t131 = t109 * t114;
t116 = cos(qJ(3));
t130 = t109 * t116;
t117 = cos(qJ(2));
t129 = t109 * t117;
t111 = cos(pkin(6));
t128 = t111 * t114;
t127 = t111 * t117;
t126 = t111 * pkin(7) + qJ(1);
t125 = t108 * t132;
t124 = t112 * t129;
t103 = t108 * pkin(1);
t123 = -pkin(7) * t133 + t103;
t98 = pkin(3) * t116 + pkin(2);
t121 = t111 * t138 + t119 * t129 + t98 * t131 + t126;
t107 = qJ(3) + qJ(4);
t102 = cos(t107);
t100 = sin(t107);
t88 = -t108 * t128 + t110 * t117;
t86 = t108 * t117 + t110 * t128;
t82 = t100 * t111 + t102 * t131;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(5) * t121 - t82 * mrSges(5,1) + mrSges(5,3) * t129 - m(6) * (pkin(4) * t82 + t121) - (t115 * t82 - t124) * mrSges(6,1) - (-t112 * t82 - t115 * t129) * mrSges(6,2) - m(7) * (-pkin(5) * t124 + t82 * t97 + t121) - (t101 * t82 - t129 * t99) * mrSges(7,1) - (-t101 * t129 - t82 * t99) * mrSges(7,2) + (-m(3) - m(4)) * t126 + (-t113 * mrSges(4,1) - t116 * mrSges(4,2) - mrSges(3,3)) * t111 + (t143 * t117 + (-m(4) * pkin(2) - mrSges(4,1) * t116 + mrSges(4,2) * t113 - mrSges(3,1)) * t114) * t109 + t141 * (t100 * t131 - t111 * t102)) * g(3) + (-mrSges(1,2) - t108 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t123 - t86 * mrSges(3,1) + mrSges(3,3) * t133 - m(4) * (pkin(2) * t86 + t123) - (-t110 * t132 + t116 * t86) * mrSges(4,1) - (-t110 * t130 - t113 * t86) * mrSges(4,2) + t144 * (t103 + t86 * t98 + (-pkin(7) - t138) * t133) + t139 * (-t100 * t133 + t102 * t86) + t141 * (t100 * t86 + t102 * t133) + t145 * (t108 * t114 - t110 * t127)) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t108 * mrSges(2,2) - m(3) * t135 - t88 * mrSges(3,1) - mrSges(3,3) * t134 - m(4) * (pkin(2) * t88 + t135) - (t116 * t88 + t125) * mrSges(4,1) - (t108 * t130 - t113 * t88) * mrSges(4,2) + t144 * (pkin(3) * t125 + t88 * t98 + t135) + t139 * (t100 * t134 + t102 * t88) + t141 * (t100 * t88 - t102 * t134) + t145 * (t108 * t127 + t110 * t114)) * g(1);
U  = t1;
