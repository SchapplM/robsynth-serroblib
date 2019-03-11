% Calculate potential energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:03:52
% EndTime: 2019-03-08 21:03:52
% DurationCPUTime: 0.68s
% Computational Cost: add. (332->97), mult. (532->112), div. (0->0), fcn. (603->12), ass. (0->49)
t141 = -m(6) - m(7);
t140 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t139 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t111 = sin(qJ(6));
t114 = cos(qJ(6));
t138 = -t111 * mrSges(7,1) - t114 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t137 = m(7) * pkin(5) + mrSges(7,1) * t114 - mrSges(7,2) * t111 + mrSges(6,1) + mrSges(5,3);
t136 = -t137 - t140;
t112 = sin(qJ(3));
t135 = pkin(3) * t112;
t108 = cos(pkin(10));
t106 = sin(pkin(10));
t107 = sin(pkin(6));
t133 = t106 * t107;
t134 = t108 * pkin(1) + pkin(7) * t133;
t132 = t107 * t108;
t131 = t107 * t112;
t113 = sin(qJ(2));
t130 = t107 * t113;
t115 = cos(qJ(3));
t129 = t107 * t115;
t116 = cos(qJ(2));
t128 = t107 * t116;
t109 = cos(pkin(6));
t127 = t109 * t113;
t126 = t109 * t116;
t125 = t109 * pkin(7) + qJ(1);
t124 = t106 * t131;
t102 = t106 * pkin(1);
t123 = -pkin(7) * t132 + t102;
t110 = -qJ(4) - pkin(8);
t88 = t106 * t126 + t108 * t113;
t89 = -t106 * t127 + t108 * t116;
t99 = pkin(3) * t115 + pkin(2);
t122 = pkin(3) * t124 - t88 * t110 + t89 * t99 + t134;
t121 = t109 * t135 + t110 * t128 + t99 * t130 + t125;
t86 = t106 * t113 - t108 * t126;
t87 = t106 * t116 + t108 * t127;
t118 = t102 + t87 * t99 - t86 * t110 + (-pkin(7) - t135) * t132;
t105 = qJ(3) + pkin(11);
t101 = cos(t105);
t100 = sin(t105);
t83 = t100 * t109 + t101 * t130;
t82 = t100 * t130 - t109 * t101;
t78 = t100 * t133 + t101 * t89;
t77 = t100 * t89 - t101 * t133;
t76 = -t100 * t132 + t101 * t87;
t75 = t100 * t87 + t101 * t132;
t1 = (-m(2) * qJ(1) - m(5) * t121 - mrSges(1,3) - mrSges(2,3) + t141 * (t83 * pkin(4) + t82 * qJ(5) + t121) + t138 * t82 + t137 * t128 + (-m(3) - m(4)) * t125 + (-t112 * mrSges(4,1) - t115 * mrSges(4,2) - mrSges(3,3)) * t109 + (t140 * t116 + (-m(4) * pkin(2) - t115 * mrSges(4,1) + t112 * mrSges(4,2) - mrSges(3,1)) * t113) * t107 + t139 * t83) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t106 - mrSges(2,2) * t108 - m(3) * t123 - t87 * mrSges(3,1) + mrSges(3,3) * t132 - m(4) * (pkin(2) * t87 + t123) - (-t108 * t131 + t115 * t87) * mrSges(4,1) - (-t108 * t129 - t112 * t87) * mrSges(4,2) - m(5) * t118 + t141 * (t76 * pkin(4) + qJ(5) * t75 + t118) + t139 * t76 + t138 * t75 + t136 * t86) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t108 + mrSges(2,2) * t106 - m(3) * t134 - t89 * mrSges(3,1) - mrSges(3,3) * t133 - m(4) * (pkin(2) * t89 + t134) - (t115 * t89 + t124) * mrSges(4,1) - (t106 * t129 - t112 * t89) * mrSges(4,2) - m(5) * t122 + t141 * (t78 * pkin(4) + qJ(5) * t77 + t122) + t139 * t78 + t138 * t77 + t136 * t88) * g(1);
U  = t1;
