% Calculate potential energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:20
% EndTime: 2019-03-09 08:41:20
% DurationCPUTime: 0.64s
% Computational Cost: add. (183->81), mult. (261->80), div. (0->0), fcn. (240->8), ass. (0->38)
t123 = mrSges(3,2) - mrSges(4,3);
t122 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2);
t120 = m(4) + m(5);
t119 = -m(6) - m(7);
t87 = sin(qJ(2));
t103 = qJ(3) * t87;
t88 = sin(qJ(1));
t89 = cos(qJ(2));
t106 = t88 * t89;
t118 = pkin(2) * t106 + t88 * t103;
t117 = mrSges(4,1) + mrSges(3,3) - mrSges(2,2);
t116 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,2);
t115 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t114 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t113 = t122 * t89 + t123 * t87 - mrSges(2,1);
t112 = t87 * pkin(2) + pkin(6);
t83 = pkin(9) + qJ(5);
t77 = sin(t83);
t111 = t77 * t88;
t90 = cos(qJ(1));
t110 = t87 * t90;
t78 = cos(t83);
t109 = t88 * t78;
t84 = sin(pkin(9));
t108 = t88 * t84;
t85 = cos(pkin(9));
t107 = t88 * t85;
t105 = t89 * t90;
t104 = t90 * pkin(1) + t88 * pkin(7);
t101 = t84 * t110;
t100 = t87 * t108;
t81 = t88 * pkin(1);
t99 = t81 + t118;
t98 = -pkin(7) * t90 + t81;
t96 = pkin(2) * t105 + t90 * t103 + t104;
t86 = -pkin(8) - qJ(4);
t76 = pkin(4) * t85 + pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + t119 * (-t86 * t87 + t112) - t120 * t112 + (-m(2) - m(3)) * pkin(6) + (t84 * mrSges(5,1) + t85 * mrSges(5,2) + t119 * (-pkin(4) * t84 - qJ(3)) - t114 * t78 + t115 * t77 + t120 * qJ(3) - t123) * t89 + (t116 + t122) * t87) * g(3) + (-mrSges(1,2) - m(3) * t98 - m(4) * (t98 + t118) - m(5) * t99 - t100 * mrSges(5,1) - t87 * t107 * mrSges(5,2) + t119 * (-t86 * t106 + pkin(4) * t100 + (-pkin(7) - t76) * t90 + t99) - t115 * (t87 * t111 - t78 * t90) + t114 * (t87 * t109 + t77 * t90) + (-m(5) * (-pkin(3) - pkin(7)) + t85 * mrSges(5,1) - t84 * mrSges(5,2) + t117) * t90 + t113 * t88 + t116 * t106) * g(2) + (-mrSges(1,1) - m(3) * t104 - (t101 + t107) * mrSges(5,1) - (t85 * t110 - t108) * mrSges(5,2) - t120 * t96 + t119 * (pkin(4) * t101 - t86 * t105 + t88 * t76 + t96) - t115 * (t77 * t110 + t109) - t114 * (-t78 * t110 + t111) + t113 * t90 + (-m(5) * pkin(3) - t117) * t88 + t116 * t105) * g(1);
U  = t1;
