% Calculate potential energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:34
% EndTime: 2019-03-08 19:46:34
% DurationCPUTime: 0.73s
% Computational Cost: add. (272->101), mult. (545->112), div. (0->0), fcn. (622->12), ass. (0->45)
t136 = -m(5) - m(6);
t135 = mrSges(4,2) - mrSges(3,1);
t134 = mrSges(3,3) + mrSges(4,1);
t133 = -mrSges(4,3) + mrSges(3,2);
t132 = m(6) * qJ(5) - m(7) * (-pkin(9) - qJ(5)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t102 = cos(pkin(11));
t99 = sin(pkin(11));
t129 = pkin(5) * t99;
t98 = pkin(11) + qJ(6);
t92 = sin(t98);
t93 = cos(t98);
t131 = -m(7) * (pkin(8) + t129) - t99 * mrSges(6,1) - t92 * mrSges(7,1) - t102 * mrSges(6,2) - t93 * mrSges(7,2) - mrSges(5,3) + t135;
t91 = pkin(5) * t102 + pkin(4);
t130 = -m(6) * pkin(4) - m(7) * t91 - t102 * mrSges(6,1) - t93 * mrSges(7,1) + t99 * mrSges(6,2) + t92 * mrSges(7,2) - mrSges(5,1);
t103 = cos(pkin(10));
t100 = sin(pkin(10));
t101 = sin(pkin(6));
t126 = t100 * t101;
t128 = t103 * pkin(1) + pkin(7) * t126;
t104 = cos(pkin(6));
t127 = t104 * pkin(7) + qJ(1);
t107 = sin(qJ(2));
t125 = t101 * t107;
t109 = cos(qJ(2));
t124 = t101 * t109;
t123 = t103 * t101;
t122 = t104 * t107;
t121 = t104 * t109;
t120 = pkin(7) * t123;
t119 = pkin(2) * t125 + t127;
t79 = t100 * t107 - t103 * t121;
t80 = t100 * t109 + t103 * t122;
t94 = t100 * pkin(1);
t117 = t80 * pkin(2) + qJ(3) * t79 + t94;
t116 = t104 * pkin(3) + pkin(8) * t125 + t119;
t81 = t100 * t121 + t103 * t107;
t82 = -t100 * t122 + t103 * t109;
t115 = t82 * pkin(2) + qJ(3) * t81 + t128;
t114 = pkin(3) * t126 + t115;
t112 = -qJ(3) * t124 + t116;
t111 = (-pkin(3) - pkin(7)) * t123 + t117;
t108 = cos(qJ(4));
t106 = sin(qJ(4));
t84 = t104 * t108 - t106 * t124;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t127 - m(4) * t119 - m(5) * t112 - t84 * mrSges(5,1) - mrSges(5,3) * t125 - m(6) * (t84 * pkin(4) + t112) - (t102 * t84 + t99 * t125) * mrSges(6,1) - (t102 * t125 - t84 * t99) * mrSges(6,2) - m(7) * (t84 * t91 + t116) - (t125 * t92 + t84 * t93) * mrSges(7,1) - (t125 * t93 - t84 * t92) * mrSges(7,2) - t134 * t104 + (((m(4) + m(7)) * qJ(3) - t133) * t109 + (-m(7) * t129 + t135) * t107) * t101 - t132 * (t104 * t106 + t108 * t124)) * g(3) + (-mrSges(1,2) - t100 * mrSges(2,1) - t103 * mrSges(2,2) - m(3) * (t94 - t120) - m(4) * (t117 - t120) - m(7) * t111 + t133 * t79 + t134 * t123 + t136 * (pkin(8) * t80 + t111) + t130 * (t106 * t79 - t108 * t123) + t131 * t80 + t132 * (t106 * t123 + t108 * t79)) * g(2) + (-m(3) * t128 - m(4) * t115 - m(7) * t114 - t103 * mrSges(2,1) + t100 * mrSges(2,2) - mrSges(1,1) + t133 * t81 - t134 * t126 + t136 * (pkin(8) * t82 + t114) + t130 * (t106 * t81 + t108 * t126) + t131 * t82 - t132 * (t106 * t126 - t81 * t108)) * g(1);
U  = t1;
