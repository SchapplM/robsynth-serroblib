% Calculate potential energy for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:12
% EndTime: 2019-03-08 22:16:13
% DurationCPUTime: 0.73s
% Computational Cost: add. (335->105), mult. (613->121), div. (0->0), fcn. (717->14), ass. (0->46)
t135 = -m(4) - m(5);
t134 = -m(6) - m(7);
t101 = sin(pkin(12));
t118 = pkin(4) * t101 + pkin(8);
t106 = -pkin(9) - qJ(4);
t133 = -m(5) * qJ(4) + m(6) * t106 + m(7) * (-pkin(10) + t106) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t104 = cos(pkin(12));
t100 = pkin(12) + qJ(5);
t92 = sin(t100);
t130 = pkin(5) * t92 + t118;
t94 = qJ(6) + t100;
t89 = sin(t94);
t90 = cos(t94);
t93 = cos(t100);
t132 = -m(6) * t118 - m(7) * t130 - t101 * mrSges(5,1) - t92 * mrSges(6,1) - t89 * mrSges(7,1) - t104 * mrSges(5,2) - t93 * mrSges(6,2) - t90 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t91 = pkin(4) * t104 + pkin(3);
t82 = pkin(5) * t93 + t91;
t131 = -m(5) * pkin(3) - m(6) * t91 - m(7) * t82 - mrSges(5,1) * t104 - mrSges(6,1) * t93 - mrSges(7,1) * t90 + mrSges(5,2) * t101 + mrSges(6,2) * t92 + mrSges(7,2) * t89 - mrSges(4,1);
t129 = cos(qJ(3));
t105 = cos(pkin(11));
t102 = sin(pkin(11));
t103 = sin(pkin(6));
t124 = t102 * t103;
t127 = pkin(1) * t105 + pkin(7) * t124;
t125 = cos(pkin(6));
t126 = pkin(7) * t125 + qJ(1);
t108 = sin(qJ(2));
t123 = t103 * t108;
t109 = cos(qJ(2));
t122 = t103 * t109;
t121 = t105 * t103;
t116 = t108 * t125;
t79 = -t102 * t116 + t105 * t109;
t120 = pkin(2) * t79 + t127;
t119 = pkin(2) * t123 + t126;
t117 = t103 * t129;
t115 = t109 * t125;
t114 = pkin(1) * t102 - pkin(7) * t121;
t77 = t102 * t109 + t105 * t116;
t112 = pkin(2) * t77 + t114;
t111 = -pkin(8) * t122 + t119;
t107 = sin(qJ(3));
t81 = t107 * t125 + t108 * t117;
t78 = t102 * t115 + t105 * t108;
t76 = t102 * t108 - t105 * t115;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t126 - t125 * mrSges(3,3) - (t108 * mrSges(3,1) + t109 * mrSges(3,2)) * t103 - m(4) * t111 - t81 * mrSges(4,1) + mrSges(4,3) * t122 - m(5) * (pkin(3) * t81 + t111) - (-t101 * t122 + t104 * t81) * mrSges(5,1) - (-t101 * t81 - t104 * t122) * mrSges(5,2) - m(6) * (-t118 * t122 + t81 * t91 + t119) - (-t122 * t92 + t81 * t93) * mrSges(6,1) - (-t122 * t93 - t81 * t92) * mrSges(6,2) - m(7) * (-t122 * t130 + t81 * t82 + t119) - (-t122 * t89 + t81 * t90) * mrSges(7,1) - (-t122 * t90 - t81 * t89) * mrSges(7,2) + t133 * (t107 * t123 - t125 * t129)) * g(3) + (-m(3) * t114 - t102 * mrSges(2,1) - t77 * mrSges(3,1) - t105 * mrSges(2,2) + mrSges(3,3) * t121 - mrSges(1,2) + t134 * t112 + t135 * (pkin(8) * t76 + t112) + t131 * (-t107 * t121 + t129 * t77) + t132 * t76 + t133 * (t105 * t117 + t107 * t77)) * g(2) + (-m(3) * t127 - t105 * mrSges(2,1) - t79 * mrSges(3,1) + t102 * mrSges(2,2) - mrSges(3,3) * t124 - mrSges(1,1) + t134 * t120 + t135 * (pkin(8) * t78 + t120) + t131 * (t107 * t124 + t129 * t79) + t132 * t78 + t133 * (-t102 * t117 + t107 * t79)) * g(1);
U  = t1;
