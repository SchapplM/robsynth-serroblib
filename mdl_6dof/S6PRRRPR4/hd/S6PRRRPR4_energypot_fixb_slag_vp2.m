% Calculate potential energy for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:33
% EndTime: 2019-03-08 23:15:33
% DurationCPUTime: 0.74s
% Computational Cost: add. (335->105), mult. (613->121), div. (0->0), fcn. (717->14), ass. (0->46)
t135 = -m(4) - m(5);
t134 = -m(6) - m(7);
t105 = sin(qJ(4));
t118 = pkin(4) * t105 + pkin(8);
t104 = -qJ(5) - pkin(9);
t133 = -m(5) * pkin(9) + m(6) * t104 + m(7) * (-pkin(10) + t104) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t108 = cos(qJ(4));
t100 = qJ(4) + pkin(12);
t92 = sin(t100);
t130 = pkin(5) * t92 + t118;
t94 = qJ(6) + t100;
t89 = sin(t94);
t90 = cos(t94);
t93 = cos(t100);
t132 = -m(6) * t118 - m(7) * t130 - t105 * mrSges(5,1) - t92 * mrSges(6,1) - t89 * mrSges(7,1) - t108 * mrSges(5,2) - t93 * mrSges(6,2) - t90 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t91 = t108 * pkin(4) + pkin(3);
t82 = pkin(5) * t93 + t91;
t131 = -m(5) * pkin(3) - m(6) * t91 - m(7) * t82 - t108 * mrSges(5,1) - t93 * mrSges(6,1) - t90 * mrSges(7,1) + t105 * mrSges(5,2) + t92 * mrSges(6,2) + t89 * mrSges(7,2) - mrSges(4,1);
t129 = cos(qJ(3));
t103 = cos(pkin(11));
t101 = sin(pkin(11));
t102 = sin(pkin(6));
t124 = t101 * t102;
t127 = t103 * pkin(1) + pkin(7) * t124;
t125 = cos(pkin(6));
t126 = t125 * pkin(7) + qJ(1);
t107 = sin(qJ(2));
t123 = t102 * t107;
t109 = cos(qJ(2));
t122 = t102 * t109;
t121 = t103 * t102;
t116 = t107 * t125;
t79 = -t101 * t116 + t103 * t109;
t120 = t79 * pkin(2) + t127;
t119 = pkin(2) * t123 + t126;
t117 = t102 * t129;
t115 = t109 * t125;
t114 = t101 * pkin(1) - pkin(7) * t121;
t77 = t101 * t109 + t103 * t116;
t112 = t77 * pkin(2) + t114;
t111 = -pkin(8) * t122 + t119;
t106 = sin(qJ(3));
t81 = t106 * t125 + t107 * t117;
t78 = t101 * t115 + t103 * t107;
t76 = t101 * t107 - t103 * t115;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t126 - t125 * mrSges(3,3) - (t107 * mrSges(3,1) + t109 * mrSges(3,2)) * t102 - m(4) * t111 - t81 * mrSges(4,1) + mrSges(4,3) * t122 - m(5) * (t81 * pkin(3) + t111) - (-t105 * t122 + t81 * t108) * mrSges(5,1) - (-t81 * t105 - t108 * t122) * mrSges(5,2) - m(6) * (-t118 * t122 + t81 * t91 + t119) - (-t122 * t92 + t81 * t93) * mrSges(6,1) - (-t122 * t93 - t81 * t92) * mrSges(6,2) - m(7) * (-t122 * t130 + t81 * t82 + t119) - (-t122 * t89 + t81 * t90) * mrSges(7,1) - (-t122 * t90 - t81 * t89) * mrSges(7,2) + t133 * (t106 * t123 - t125 * t129)) * g(3) + (-m(3) * t114 - t101 * mrSges(2,1) - t77 * mrSges(3,1) - t103 * mrSges(2,2) + mrSges(3,3) * t121 - mrSges(1,2) + t134 * t112 + t135 * (pkin(8) * t76 + t112) + t131 * (-t106 * t121 + t129 * t77) + t132 * t76 + t133 * (t103 * t117 + t77 * t106)) * g(2) + (-m(3) * t127 - t103 * mrSges(2,1) - t79 * mrSges(3,1) + t101 * mrSges(2,2) - mrSges(3,3) * t124 - mrSges(1,1) + t134 * t120 + t135 * (pkin(8) * t78 + t120) + t131 * (t106 * t124 + t129 * t79) + t132 * t78 + t133 * (-t101 * t117 + t106 * t79)) * g(1);
U  = t1;
