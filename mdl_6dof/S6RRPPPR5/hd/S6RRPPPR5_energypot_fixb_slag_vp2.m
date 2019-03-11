% Calculate potential energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:15
% EndTime: 2019-03-09 08:21:16
% DurationCPUTime: 0.55s
% Computational Cost: add. (181->80), mult. (344->78), div. (0->0), fcn. (351->8), ass. (0->41)
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t122 = -t90 * mrSges(3,1) + t87 * mrSges(3,2) - mrSges(2,1);
t121 = mrSges(2,2) - mrSges(3,3);
t85 = cos(pkin(9));
t113 = t85 * t87;
t84 = sin(pkin(9));
t120 = t87 * t84 * qJ(4) + pkin(3) * t113;
t119 = mrSges(5,1) + mrSges(4,3) + mrSges(7,3) - mrSges(6,2);
t118 = -m(7) * pkin(8) - t119;
t86 = sin(qJ(6));
t89 = cos(qJ(6));
t117 = -t86 * mrSges(7,1) - t89 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t116 = -t89 * mrSges(7,1) + t86 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t115 = -m(7) * (pkin(5) + qJ(4)) + t116;
t114 = t87 * pkin(2) + pkin(6);
t88 = sin(qJ(1));
t112 = t87 * t88;
t91 = cos(qJ(1));
t111 = t87 * t91;
t110 = t88 * t90;
t109 = t90 * t91;
t108 = t91 * t84;
t107 = -pkin(4) - qJ(3);
t105 = t91 * pkin(1) + t88 * pkin(7);
t104 = qJ(3) * t87;
t65 = t110 * t84 + t85 * t91;
t103 = t65 * qJ(4);
t67 = t108 * t90 - t88 * t85;
t102 = t67 * qJ(4);
t101 = t88 * pkin(1) - pkin(7) * t91;
t100 = pkin(2) * t109 + t91 * t104 + t105;
t99 = -qJ(3) * t90 + t114;
t68 = t109 * t85 + t88 * t84;
t97 = t68 * pkin(3) + t100;
t95 = pkin(2) * t110 + t88 * t104 + t101;
t66 = t110 * t85 - t108;
t94 = t66 * pkin(3) + t95;
t93 = pkin(4) * t111 + t68 * qJ(5) + t97;
t92 = pkin(4) * t112 + t66 * qJ(5) + t94;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t99 - m(5) * (t99 + t120) + (-m(6) - m(7)) * (qJ(5) * t113 + t114 + t120) + (-m(2) - m(3)) * pkin(6) + (-mrSges(3,2) - m(6) * t107 - m(7) * (-pkin(8) + t107) + t119) * t90 + (-mrSges(3,1) + t117 * t85 + (-m(7) * pkin(5) + t116) * t84) * t87) * g(3) + (-mrSges(1,2) - m(3) * t101 - m(4) * t95 - m(5) * (t94 + t103) - m(6) * (t92 + t103) - m(7) * t92 - t121 * t91 + t122 * t88 + t117 * t66 + t115 * t65 + t118 * t112) * g(2) + (-mrSges(1,1) - m(3) * t105 - m(4) * t100 - m(5) * (t97 + t102) - m(6) * (t93 + t102) - m(7) * t93 + t122 * t91 + t121 * t88 + t117 * t68 + t115 * t67 + t118 * t111) * g(1);
U  = t1;
