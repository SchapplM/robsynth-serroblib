% Calculate potential energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:22
% EndTime: 2019-03-09 16:21:23
% DurationCPUTime: 0.60s
% Computational Cost: add. (301->91), mult. (621->99), div. (0->0), fcn. (727->12), ass. (0->51)
t134 = pkin(4) + pkin(9);
t138 = -m(6) * qJ(5) + m(7) * (-pkin(10) - qJ(5)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t101 = cos(pkin(11));
t98 = pkin(11) + qJ(6);
t93 = sin(t98);
t94 = cos(t98);
t99 = sin(pkin(11));
t137 = -m(7) * (pkin(5) * t99 + qJ(4)) - t99 * mrSges(6,1) - t93 * mrSges(7,1) - t101 * mrSges(6,2) - t94 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t136 = m(6) * t134 + m(7) * (t101 * pkin(5) + t134) + t101 * mrSges(6,1) + t94 * mrSges(7,1) - t99 * mrSges(6,2) - t93 * mrSges(7,2) + mrSges(5,1) + mrSges(4,3);
t135 = mrSges(3,2) - t136;
t104 = sin(qJ(2));
t105 = sin(qJ(1));
t106 = cos(qJ(2));
t107 = cos(qJ(1));
t127 = cos(pkin(6));
t115 = t107 * t127;
t82 = t105 * t104 - t106 * t115;
t133 = t82 * pkin(9);
t116 = t105 * t127;
t84 = t107 * t104 + t106 * t116;
t132 = t84 * pkin(9);
t130 = t127 * pkin(8) + pkin(7);
t129 = cos(qJ(3));
t100 = sin(pkin(6));
t124 = t105 * t100;
t128 = t107 * pkin(1) + pkin(8) * t124;
t126 = t100 * t104;
t125 = t100 * t106;
t123 = t107 * t100;
t122 = pkin(2) * t126 + t130;
t85 = -t104 * t116 + t107 * t106;
t121 = t85 * pkin(2) + t128;
t120 = pkin(9) * t125;
t119 = t100 * t129;
t103 = sin(qJ(3));
t81 = t127 * t103 + t104 * t119;
t117 = t81 * pkin(3) + t122;
t76 = t103 * t124 + t85 * t129;
t114 = t76 * pkin(3) + t121;
t113 = t105 * pkin(1) - pkin(8) * t123;
t83 = t104 * t115 + t105 * t106;
t112 = t83 * pkin(2) + t113;
t80 = t103 * t126 - t127 * t129;
t111 = t80 * qJ(4) + t117;
t74 = -t103 * t123 + t83 * t129;
t110 = t74 * pkin(3) + t112;
t75 = t85 * t103 - t105 * t119;
t109 = t75 * qJ(4) + t114;
t73 = t83 * t103 + t107 * t119;
t108 = t73 * qJ(4) + t110;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t130 - t127 * mrSges(3,3) - (t104 * mrSges(3,1) + t106 * mrSges(3,2)) * t100 - m(4) * (-t120 + t122) - m(5) * (t111 - t120) - m(6) * t111 - m(7) * t117 + t137 * t80 + t136 * t125 + t138 * t81) * g(3) + (-mrSges(1,2) - t105 * mrSges(2,1) - t107 * mrSges(2,2) - m(3) * t113 - t83 * mrSges(3,1) + mrSges(3,3) * t123 - m(4) * (t112 + t133) - m(5) * (t108 + t133) - m(6) * t108 - m(7) * t110 + t137 * t73 + t135 * t82 + t138 * t74) * g(2) + (-mrSges(1,1) - t107 * mrSges(2,1) + t105 * mrSges(2,2) - m(3) * t128 - t85 * mrSges(3,1) - mrSges(3,3) * t124 - m(4) * (t121 + t132) - m(5) * (t109 + t132) - m(6) * t109 - m(7) * t114 + t137 * t75 + t135 * t84 + t138 * t76) * g(1);
U  = t1;
