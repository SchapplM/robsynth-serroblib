% Calculate potential energy for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:16
% EndTime: 2019-03-09 12:50:17
% DurationCPUTime: 0.57s
% Computational Cost: add. (183->82), mult. (261->80), div. (0->0), fcn. (240->8), ass. (0->39)
t121 = -mrSges(3,1) + mrSges(4,2);
t120 = mrSges(3,2) - mrSges(4,3);
t119 = m(4) + m(5);
t118 = -m(6) - m(7);
t84 = sin(qJ(2));
t117 = t84 * mrSges(5,2);
t101 = qJ(3) * t84;
t85 = sin(qJ(1));
t87 = cos(qJ(2));
t106 = t85 * t87;
t116 = pkin(2) * t106 + t85 * t101;
t115 = t120 * t84 + t121 * t87 - mrSges(2,1);
t114 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t113 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t112 = -m(5) * pkin(8) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t83 = sin(qJ(4));
t111 = -t83 * mrSges(5,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t110 = t84 * pkin(2) + pkin(6);
t109 = t84 * t85;
t88 = cos(qJ(1));
t108 = t84 * t88;
t86 = cos(qJ(4));
t107 = t85 * t86;
t105 = t86 * t88;
t104 = t87 * t88;
t89 = -pkin(9) - pkin(8);
t103 = t87 * t89;
t102 = t88 * pkin(1) + t85 * pkin(7);
t100 = t83 * t109;
t99 = t83 * t108;
t80 = t85 * pkin(1);
t98 = t80 + t116;
t97 = -t88 * pkin(7) + t80;
t95 = pkin(2) * t104 + t88 * t101 + t102;
t82 = qJ(4) + qJ(5);
t77 = cos(t82);
t76 = sin(t82);
t75 = pkin(4) * t86 + pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + t118 * (-t84 * t89 + t110) - t119 * t110 + (-m(2) - m(3)) * pkin(6) + (t83 * mrSges(5,1) + t86 * mrSges(5,2) + t118 * (-pkin(4) * t83 - qJ(3)) - t113 * t77 + t114 * t76 + t119 * qJ(3) - t120) * t87 + (t112 + t121) * t84) * g(3) + (-mrSges(1,2) - m(3) * t97 - m(4) * (t97 + t116) - m(5) * t98 - (t100 - t105) * mrSges(5,1) - t107 * t117 + t118 * (-t85 * t103 + pkin(4) * t100 + (-pkin(7) - t75) * t88 + t98) - t114 * (t76 * t109 - t77 * t88) + t113 * (t77 * t109 + t76 * t88) + (-m(5) * (-pkin(3) - pkin(7)) + t111) * t88 + t115 * t85 + t112 * t106) * g(2) + (-mrSges(1,1) - m(3) * t102 - (t99 + t107) * mrSges(5,1) - t105 * t117 - t119 * t95 + t118 * (pkin(4) * t99 - t88 * t103 + t85 * t75 + t95) - t114 * (t76 * t108 + t77 * t85) - t113 * (-t77 * t108 + t76 * t85) + t115 * t88 + (-m(5) * pkin(3) - t111) * t85 + t112 * t104) * g(1);
U  = t1;
