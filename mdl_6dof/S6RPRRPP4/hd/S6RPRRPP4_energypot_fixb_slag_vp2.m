% Calculate potential energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:45
% EndTime: 2019-03-09 04:38:46
% DurationCPUTime: 0.43s
% Computational Cost: add. (239->73), mult. (236->70), div. (0->0), fcn. (215->10), ass. (0->37)
t89 = sin(qJ(4));
t91 = cos(qJ(4));
t121 = -m(5) * pkin(3) - t91 * mrSges(5,1) + t89 * mrSges(5,2) - mrSges(4,1);
t120 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(5,3);
t119 = -m(4) - m(5);
t118 = -m(6) - m(7);
t117 = -mrSges(6,3) - mrSges(7,2);
t83 = pkin(9) + qJ(3);
t78 = sin(t83);
t80 = cos(t83);
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t115 = -m(3) * pkin(1) - t86 * mrSges(3,1) + t85 * mrSges(3,2) + t120 * t78 + t121 * t80 - mrSges(2,1);
t114 = m(3) * qJ(2) + t89 * mrSges(5,1) + t91 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t113 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t112 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t111 = pkin(4) * t89;
t110 = t85 * pkin(2) + pkin(6);
t90 = sin(qJ(1));
t109 = t78 * t90;
t92 = cos(qJ(1));
t108 = t78 * t92;
t107 = t80 * t92;
t84 = qJ(4) + pkin(10);
t81 = cos(t84);
t106 = t81 * t92;
t79 = sin(t84);
t105 = t90 * t79;
t104 = t90 * t81;
t75 = pkin(2) * t86 + pkin(1);
t88 = -pkin(7) - qJ(2);
t103 = t90 * t75 + t92 * t88;
t71 = t92 * t75;
t101 = -t90 * t88 + t71;
t87 = -qJ(5) - pkin(8);
t77 = pkin(4) * t91 + pkin(3);
t1 = (-t85 * mrSges(3,1) - t86 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t119 * t110 + t118 * (t78 * t77 + t80 * t87 + t110) + (-m(2) - m(3)) * pkin(6) + (-t117 - t120) * t80 + (t112 * t79 + t113 * t81 + t121) * t78) * g(3) + (-mrSges(1,2) + t118 * (t90 * t80 * t77 - t87 * t109 - t92 * t111 + t103) + t113 * (t80 * t104 - t79 * t92) + t112 * (t80 * t105 + t106) + t117 * t109 + t119 * t103 + t114 * t92 + t115 * t90) * g(2) + (-m(4) * t101 - m(5) * t71 - mrSges(1,1) + t118 * (t77 * t107 - t87 * t108 + t90 * t111 + t101) + t113 * (t80 * t106 + t105) + t112 * (t79 * t107 - t104) + t117 * t108 + t115 * t92 + (m(5) * t88 - t114) * t90) * g(1);
U  = t1;
