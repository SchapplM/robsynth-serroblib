% Calculate potential energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:03
% EndTime: 2019-03-08 21:08:04
% DurationCPUTime: 0.60s
% Computational Cost: add. (278->91), mult. (599->95), div. (0->0), fcn. (697->10), ass. (0->49)
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t139 = -t103 * mrSges(7,1) - t106 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t101 = sin(pkin(6));
t107 = cos(qJ(2));
t123 = t101 * t107;
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t131 = cos(qJ(3));
t120 = t101 * t131;
t126 = cos(pkin(6));
t89 = t126 * t104 + t105 * t120;
t138 = t89 * pkin(4) + qJ(5) * t123;
t137 = -m(7) * pkin(9) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t135 = -m(7) * (pkin(5) + qJ(4)) - t106 * mrSges(7,1) + t103 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t134 = mrSges(3,2) - t139 + (-m(6) - m(7)) * (pkin(8) - qJ(5));
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t118 = t107 * t126;
t84 = t100 * t105 - t102 * t118;
t133 = pkin(8) * t84;
t86 = t100 * t118 + t102 * t105;
t132 = pkin(8) * t86;
t125 = t100 * t101;
t128 = t102 * pkin(1) + pkin(7) * t125;
t127 = t126 * pkin(7) + qJ(1);
t124 = t101 * t105;
t122 = t102 * t101;
t119 = t105 * t126;
t87 = -t100 * t119 + t102 * t107;
t121 = t87 * pkin(2) + t128;
t79 = t104 * t125 + t87 * t131;
t117 = t79 * pkin(3) + t121;
t116 = t100 * pkin(1) - pkin(7) * t122;
t85 = t100 * t107 + t102 * t119;
t114 = t85 * pkin(2) + t116;
t77 = -t104 * t122 + t85 * t131;
t113 = t77 * pkin(3) + t114;
t78 = -t100 * t120 + t104 * t87;
t112 = qJ(4) * t78 + t117;
t111 = pkin(2) * t124 - pkin(8) * t123 + t127;
t110 = t89 * pkin(3) + t111;
t76 = t102 * t120 + t85 * t104;
t109 = t76 * qJ(4) + t113;
t88 = t104 * t124 - t126 * t131;
t108 = t88 * qJ(4) + t110;
t74 = t79 * pkin(4);
t72 = t77 * pkin(4);
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t127 - t126 * mrSges(3,3) - (t105 * mrSges(3,1) + t107 * mrSges(3,2)) * t101 - m(4) * t111 - m(5) * t108 - m(6) * (t108 + t138) - m(7) * (t110 + t138) + t135 * t88 + t139 * t123 + t137 * t89) * g(3) + (-mrSges(1,2) - t100 * mrSges(2,1) - t102 * mrSges(2,2) - m(3) * t116 - t85 * mrSges(3,1) + mrSges(3,3) * t122 - m(4) * (t114 + t133) - m(5) * (t109 + t133) - m(6) * (t109 + t72) - m(7) * (t113 + t72) + t135 * t76 + t134 * t84 + t137 * t77) * g(2) + (-mrSges(1,1) - t102 * mrSges(2,1) + t100 * mrSges(2,2) - m(3) * t128 - t87 * mrSges(3,1) - mrSges(3,3) * t125 - m(4) * (t121 + t132) - m(5) * (t112 + t132) - m(6) * (t112 + t74) - m(7) * (t117 + t74) + t135 * t78 + t134 * t86 + t137 * t79) * g(1);
U  = t1;
