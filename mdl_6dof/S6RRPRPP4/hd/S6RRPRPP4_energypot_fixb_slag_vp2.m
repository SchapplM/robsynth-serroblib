% Calculate potential energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:45
% EndTime: 2019-03-09 09:58:46
% DurationCPUTime: 0.61s
% Computational Cost: add. (183->82), mult. (261->80), div. (0->0), fcn. (240->8), ass. (0->39)
t121 = -mrSges(3,1) + mrSges(4,2);
t120 = mrSges(3,2) - mrSges(4,3);
t119 = m(4) + m(5);
t118 = -m(6) - m(7);
t85 = sin(qJ(2));
t101 = qJ(3) * t85;
t86 = sin(qJ(1));
t88 = cos(qJ(2));
t105 = t86 * t88;
t117 = pkin(2) * t105 + t86 * t101;
t116 = t120 * t85 + t121 * t88 - mrSges(2,1);
t115 = mrSges(4,1) + mrSges(3,3) - mrSges(2,2);
t114 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t113 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t112 = -m(5) * pkin(8) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t111 = t85 * pkin(2) + pkin(6);
t82 = qJ(4) + pkin(9);
t76 = sin(t82);
t110 = t76 * t86;
t89 = cos(qJ(1));
t109 = t85 * t89;
t77 = cos(t82);
t108 = t86 * t77;
t84 = sin(qJ(4));
t107 = t86 * t84;
t87 = cos(qJ(4));
t106 = t86 * t87;
t104 = t87 * t89;
t103 = t88 * t89;
t102 = t89 * pkin(1) + t86 * pkin(7);
t100 = t84 * t109;
t99 = t85 * t107;
t80 = t86 * pkin(1);
t98 = t80 + t117;
t97 = -pkin(7) * t89 + t80;
t95 = pkin(2) * t103 + t89 * t101 + t102;
t83 = -qJ(5) - pkin(8);
t75 = pkin(4) * t87 + pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + t118 * (-t83 * t85 + t111) - t119 * t111 + (-m(2) - m(3)) * pkin(6) + (t84 * mrSges(5,1) + t87 * mrSges(5,2) + t118 * (-pkin(4) * t84 - qJ(3)) - t113 * t77 + t114 * t76 + t119 * qJ(3) - t120) * t88 + (t112 + t121) * t85) * g(3) + (-mrSges(1,2) - m(3) * t97 - m(4) * (t97 + t117) - m(5) * t98 - (t99 - t104) * mrSges(5,1) - t85 * t106 * mrSges(5,2) + t118 * (-t83 * t105 + pkin(4) * t99 + (-pkin(7) - t75) * t89 + t98) - t114 * (t85 * t110 - t77 * t89) + t113 * (t85 * t108 + t76 * t89) + (-m(5) * (-pkin(3) - pkin(7)) - t84 * mrSges(5,2) + t115) * t89 + t116 * t86 + t112 * t105) * g(2) + (-mrSges(1,1) - m(3) * t102 - (t100 + t106) * mrSges(5,1) - (t85 * t104 - t107) * mrSges(5,2) - t119 * t95 + t118 * (pkin(4) * t100 - t83 * t103 + t86 * t75 + t95) - t114 * (t76 * t109 + t108) - t113 * (-t77 * t109 + t110) + t116 * t89 + (-m(5) * pkin(3) - t115) * t86 + t112 * t103) * g(1);
U  = t1;
