% Calculate potential energy for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:05
% EndTime: 2019-03-08 20:45:06
% DurationCPUTime: 0.71s
% Computational Cost: add. (272->91), mult. (545->98), div. (0->0), fcn. (622->12), ass. (0->42)
t138 = -m(5) - m(6);
t137 = -mrSges(3,1) + mrSges(4,2);
t136 = -mrSges(3,3) - mrSges(4,1);
t135 = -mrSges(4,3) + mrSges(3,2);
t134 = m(6) * pkin(9) - m(7) * (-pkin(10) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t107 = cos(qJ(5));
t99 = qJ(5) + qJ(6);
t93 = sin(t99);
t94 = cos(t99);
t133 = -t93 * mrSges(7,1) - t107 * mrSges(6,2) - t94 * mrSges(7,2) - mrSges(5,3);
t104 = sin(qJ(5));
t132 = -m(7) * (pkin(5) * t104 + pkin(8)) - t104 * mrSges(6,1) + t133 + t137;
t131 = -m(6) * pkin(4) - m(7) * (pkin(5) * t107 + pkin(4)) - t107 * mrSges(6,1) - t94 * mrSges(7,1) + t104 * mrSges(6,2) + t93 * mrSges(7,2) - mrSges(5,1);
t102 = cos(pkin(11));
t100 = sin(pkin(11));
t101 = sin(pkin(6));
t128 = t100 * t101;
t130 = t102 * pkin(1) + pkin(7) * t128;
t103 = cos(pkin(6));
t129 = t103 * pkin(7) + qJ(1);
t106 = sin(qJ(2));
t127 = t101 * t106;
t109 = cos(qJ(2));
t126 = t101 * t109;
t125 = t102 * t101;
t124 = t103 * t106;
t123 = t103 * t109;
t121 = pkin(7) * t125;
t120 = pkin(2) * t127 + t129;
t80 = t100 * t106 - t102 * t123;
t81 = t100 * t109 + t102 * t124;
t95 = t100 * pkin(1);
t118 = t81 * pkin(2) + t80 * qJ(3) + t95;
t117 = t103 * pkin(3) + pkin(8) * t127 + t120;
t82 = t100 * t123 + t102 * t106;
t83 = -t100 * t124 + t102 * t109;
t116 = t83 * pkin(2) + t82 * qJ(3) + t130;
t115 = pkin(3) * t128 + t116;
t112 = (-pkin(3) - pkin(7)) * t125 + t118;
t108 = cos(qJ(4));
t105 = sin(qJ(4));
t1 = (-m(2) * qJ(1) - m(3) * t129 - m(4) * t120 - m(7) * t117 - mrSges(1,3) - mrSges(2,3) + t138 * (-qJ(3) * t126 + t117) + t133 * t127 + t131 * (t103 * t108 - t105 * t126) + t136 * t103 + (((m(4) + m(7)) * qJ(3) - t135) * t109 + ((-m(7) * pkin(5) - mrSges(6,1)) * t104 + t137) * t106) * t101 - t134 * (t103 * t105 + t108 * t126)) * g(3) + (-mrSges(1,2) - t100 * mrSges(2,1) - t102 * mrSges(2,2) - m(3) * (t95 - t121) - m(4) * (t118 - t121) - m(7) * t112 + t135 * t80 - t136 * t125 + t138 * (t81 * pkin(8) + t112) + t131 * (t105 * t80 - t108 * t125) + t132 * t81 + t134 * (t105 * t125 + t108 * t80)) * g(2) + (-m(3) * t130 - m(4) * t116 - m(7) * t115 - t102 * mrSges(2,1) + t100 * mrSges(2,2) - mrSges(1,1) + t135 * t82 + t136 * t128 + t138 * (t83 * pkin(8) + t115) + t131 * (t105 * t82 + t108 * t128) + t132 * t83 - t134 * (t105 * t128 - t82 * t108)) * g(1);
U  = t1;
