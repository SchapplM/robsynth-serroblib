% Calculate potential energy for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:26
% EndTime: 2019-03-09 10:49:27
% DurationCPUTime: 0.64s
% Computational Cost: add. (244->84), mult. (315->78), div. (0->0), fcn. (314->10), ass. (0->41)
t86 = sin(pkin(10));
t87 = cos(pkin(10));
t124 = -m(4) * pkin(2) - mrSges(4,1) * t87 + mrSges(4,2) * t86 - mrSges(3,1);
t123 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + mrSges(7,3);
t122 = -m(6) - m(7);
t121 = -m(3) - m(4);
t120 = -mrSges(5,3) - mrSges(6,2);
t91 = sin(qJ(1));
t93 = cos(qJ(2));
t110 = t91 * t93;
t85 = pkin(10) + qJ(4);
t80 = sin(t85);
t81 = cos(t85);
t94 = cos(qJ(1));
t67 = t80 * t110 + t81 * t94;
t68 = t81 * t110 - t80 * t94;
t119 = t68 * pkin(4) + t67 * qJ(5);
t118 = -t86 * mrSges(4,1) - t87 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t88 = -pkin(8) - qJ(3);
t90 = sin(qJ(2));
t116 = -mrSges(2,1) + t124 * t93 + (-m(7) * (-pkin(9) - t88) + t123) * t90;
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t115 = -t89 * mrSges(7,1) - t92 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t114 = -m(7) * pkin(5) - t92 * mrSges(7,1) + t89 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t113 = pkin(3) * t86;
t112 = t90 * t91;
t111 = t90 * t94;
t109 = t93 * t94;
t108 = t94 * pkin(1) + t91 * pkin(7);
t105 = t88 * t111;
t83 = t91 * pkin(1);
t104 = -pkin(7) * t94 + t83;
t78 = pkin(3) * t87 + pkin(2);
t103 = t78 * t109 + t91 * t113 + t108;
t69 = t80 * t109 - t91 * t81;
t70 = t81 * t109 + t91 * t80;
t97 = t70 * pkin(4) + t69 * qJ(5) + t103;
t71 = t78 * t110;
t96 = -t88 * t112 + t71 + t83 + (-pkin(7) - t113) * t94;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t121) * pkin(6) + (-m(7) * pkin(9) - t120 - t123) * t93 + (-m(5) + t122) * (t90 * t78 + t93 * t88 + pkin(6)) + (t122 * (pkin(4) * t81 + qJ(5) * t80) + t114 * t81 + t115 * t80 + t124) * t90) * g(3) + (-mrSges(1,2) - m(3) * t104 - m(4) * t83 - m(5) * t96 - m(6) * (t96 + t119) - m(7) * (t104 + t71 + t119) + t114 * t68 + t115 * t67 + t120 * t112 + (m(4) * pkin(7) + m(7) * t113 - t118) * t94 + t116 * t91) * g(2) + (-mrSges(1,1) - m(5) * (t103 - t105) - m(6) * (t97 - t105) - m(7) * t97 + t114 * t70 + t115 * t69 + t120 * t111 + t121 * t108 + t118 * t91 + t116 * t94) * g(1);
U  = t1;
