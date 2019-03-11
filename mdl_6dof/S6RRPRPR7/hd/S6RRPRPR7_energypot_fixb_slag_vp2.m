% Calculate potential energy for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:02
% EndTime: 2019-03-09 10:45:02
% DurationCPUTime: 0.59s
% Computational Cost: add. (208->80), mult. (299->77), div. (0->0), fcn. (293->10), ass. (0->39)
t140 = mrSges(3,2) - mrSges(4,3);
t139 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t102 = sin(qJ(2));
t106 = cos(qJ(2));
t105 = cos(qJ(4));
t101 = sin(qJ(4));
t125 = t101 * t102;
t110 = t105 * t106 + t125;
t124 = t101 * t106;
t111 = t102 * t105 - t124;
t100 = sin(qJ(6));
t104 = cos(qJ(6));
t131 = -m(7) * pkin(5) - mrSges(7,1) * t104 + mrSges(7,2) * t100 - mrSges(6,1);
t135 = -m(6) - m(7);
t98 = qJ(4) + pkin(10);
t92 = sin(t98);
t93 = cos(t98);
t77 = t102 * t92 + t106 * t93;
t138 = t135 * pkin(4) * t125 - t110 * mrSges(5,1) - t111 * mrSges(5,2) + t102 * t140 + t139 * t106 + t77 * t131 - mrSges(2,1);
t136 = -m(4) - m(5);
t103 = sin(qJ(1));
t121 = t103 * t106;
t123 = t102 * t103;
t134 = pkin(2) * t121 + qJ(3) * t123;
t133 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t129 = -t100 * mrSges(7,1) - t104 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3) - mrSges(6,3);
t128 = t102 * pkin(2) + pkin(6);
t107 = cos(qJ(1));
t126 = t107 * pkin(1) + t103 * pkin(7);
t122 = t102 * t107;
t120 = t106 * t107;
t96 = t103 * pkin(1);
t119 = t96 + t134;
t117 = -pkin(7) * t107 + t96;
t116 = pkin(2) * t120 + qJ(3) * t122 + t126;
t99 = -qJ(5) - pkin(8);
t91 = pkin(4) * t105 + pkin(3);
t85 = t102 * t91;
t1 = (-mrSges(1,3) - mrSges(2,3) - t111 * mrSges(5,1) + t110 * mrSges(5,2) - m(6) * (t85 + t128) - m(7) * (-pkin(4) * t124 + t85) + t131 * (t102 * t93 - t106 * t92) + t133 * t77 + (-m(7) + t136) * (-qJ(3) * t106 + t128) + (-m(6) * (-pkin(4) * t101 - qJ(3)) - t140) * t106 + t139 * t102 + (-m(2) - m(3)) * pkin(6)) * g(3) + (-mrSges(1,2) - m(3) * t117 - m(4) * (t117 + t134) - m(5) * t119 + t135 * (t91 * t121 + t119) + t133 * (t92 * t121 - t93 * t123) + (-m(5) * (-pkin(7) + pkin(8)) + t135 * (-pkin(7) - t99) + t129) * t107 + t138 * t103) * g(2) + (-m(3) * t126 - mrSges(1,1) + t135 * (t103 * t99 + t91 * t120 + t116) + t133 * (t92 * t120 - t93 * t122) + t136 * t116 + (m(5) * pkin(8) - t129) * t103 + t138 * t107) * g(1);
U  = t1;
