% Calculate potential energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:09
% EndTime: 2019-03-08 22:30:09
% DurationCPUTime: 0.61s
% Computational Cost: add. (301->91), mult. (621->99), div. (0->0), fcn. (727->12), ass. (0->51)
t135 = pkin(4) + pkin(8);
t139 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t99 = qJ(5) + qJ(6);
t94 = sin(t99);
t95 = cos(t99);
t138 = -m(7) * (pkin(5) * t103 + qJ(4)) - t103 * mrSges(6,1) - t94 * mrSges(7,1) - t106 * mrSges(6,2) - t95 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t137 = m(6) * t135 + m(7) * (t106 * pkin(5) + t135) + t106 * mrSges(6,1) + t95 * mrSges(7,1) - t103 * mrSges(6,2) - t94 * mrSges(7,2) + mrSges(5,1) + mrSges(4,3);
t136 = mrSges(3,2) - t137;
t100 = sin(pkin(11));
t102 = cos(pkin(11));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t128 = cos(pkin(6));
t117 = t107 * t128;
t81 = t100 * t105 - t102 * t117;
t134 = t81 * pkin(8);
t83 = t100 * t117 + t102 * t105;
t133 = t83 * pkin(8);
t131 = cos(qJ(3));
t101 = sin(pkin(6));
t127 = t100 * t101;
t130 = t102 * pkin(1) + pkin(7) * t127;
t129 = t128 * pkin(7) + qJ(1);
t126 = t101 * t105;
t125 = t101 * t107;
t124 = t102 * t101;
t118 = t105 * t128;
t84 = -t100 * t118 + t102 * t107;
t123 = t84 * pkin(2) + t130;
t122 = pkin(8) * t125;
t121 = pkin(2) * t126 + t129;
t120 = t101 * t131;
t104 = sin(qJ(3));
t77 = t104 * t127 + t84 * t131;
t116 = t77 * pkin(3) + t123;
t86 = t128 * t104 + t105 * t120;
t115 = t86 * pkin(3) + t121;
t114 = t100 * pkin(1) - pkin(7) * t124;
t82 = t100 * t107 + t102 * t118;
t113 = t82 * pkin(2) + t114;
t75 = -t104 * t124 + t82 * t131;
t112 = t75 * pkin(3) + t113;
t76 = -t100 * t120 + t84 * t104;
t111 = t76 * qJ(4) + t116;
t85 = t104 * t126 - t128 * t131;
t110 = t85 * qJ(4) + t115;
t74 = t102 * t120 + t82 * t104;
t109 = t74 * qJ(4) + t112;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t129 - t128 * mrSges(3,3) - (t105 * mrSges(3,1) + t107 * mrSges(3,2)) * t101 - m(4) * (t121 - t122) - m(5) * (t110 - t122) - m(6) * t110 - m(7) * t115 + t138 * t85 + t137 * t125 + t139 * t86) * g(3) + (-mrSges(1,2) - t100 * mrSges(2,1) - t102 * mrSges(2,2) - m(3) * t114 - t82 * mrSges(3,1) + mrSges(3,3) * t124 - m(4) * (t113 + t134) - m(5) * (t109 + t134) - m(6) * t109 - m(7) * t112 + t138 * t74 + t136 * t81 + t139 * t75) * g(2) + (-mrSges(1,1) - t102 * mrSges(2,1) + t100 * mrSges(2,2) - m(3) * t130 - t84 * mrSges(3,1) - mrSges(3,3) * t127 - m(4) * (t123 + t133) - m(5) * (t111 + t133) - m(6) * t111 - m(7) * t116 + t138 * t76 + t136 * t83 + t139 * t77) * g(1);
U  = t1;
