% Calculate potential energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:35:50
% EndTime: 2019-03-08 20:35:51
% DurationCPUTime: 0.93s
% Computational Cost: add. (355->108), mult. (554->125), div. (0->0), fcn. (633->14), ass. (0->42)
t105 = qJ(5) + qJ(6);
t100 = cos(t105);
t112 = -pkin(8) - qJ(3);
t113 = sin(qJ(5));
t115 = cos(qJ(5));
t139 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t140 = -m(5) - m(6) - m(7);
t99 = sin(t105);
t141 = -t99 * mrSges(7,1) - t115 * mrSges(6,2) - t100 * mrSges(7,2) - t140 * t112 - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t113 - t139;
t137 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t96 = pkin(5) * t115 + pkin(4);
t135 = -m(6) * pkin(4) - m(7) * t96 - t115 * mrSges(6,1) - t100 * mrSges(7,1) + t113 * mrSges(6,2) + t99 * mrSges(7,2) - mrSges(5,1);
t106 = sin(pkin(12));
t134 = pkin(3) * t106;
t110 = cos(pkin(11));
t107 = sin(pkin(11));
t108 = sin(pkin(6));
t130 = t107 * t108;
t131 = t110 * pkin(1) + pkin(7) * t130;
t129 = t108 * t110;
t114 = sin(qJ(2));
t128 = t108 * t114;
t116 = cos(qJ(2));
t127 = t108 * t116;
t111 = cos(pkin(6));
t126 = t111 * t114;
t125 = t111 * t116;
t124 = t111 * pkin(7) + qJ(1);
t123 = t106 * t130;
t122 = t113 * t127;
t101 = t107 * pkin(1);
t121 = -pkin(7) * t129 + t101;
t109 = cos(pkin(12));
t95 = pkin(3) * t109 + pkin(2);
t119 = t111 * t134 + t112 * t127 + t95 * t128 + t124;
t104 = pkin(12) + qJ(4);
t98 = cos(t104);
t97 = sin(t104);
t86 = -t107 * t126 + t110 * t116;
t84 = t107 * t116 + t110 * t126;
t80 = t111 * t97 + t128 * t98;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(5) * t119 - t80 * mrSges(5,1) + mrSges(5,3) * t127 - m(6) * (pkin(4) * t80 + t119) - (t115 * t80 - t122) * mrSges(6,1) - (-t113 * t80 - t115 * t127) * mrSges(6,2) - m(7) * (-pkin(5) * t122 + t80 * t96 + t119) - (t100 * t80 - t127 * t99) * mrSges(7,1) - (-t100 * t127 - t80 * t99) * mrSges(7,2) + (-m(3) - m(4)) * t124 + (-t106 * mrSges(4,1) - t109 * mrSges(4,2) - mrSges(3,3)) * t111 + (t139 * t116 + (-m(4) * pkin(2) - t109 * mrSges(4,1) + t106 * mrSges(4,2) - mrSges(3,1)) * t114) * t108 + t137 * (-t111 * t98 + t128 * t97)) * g(3) + (-mrSges(1,2) - t107 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t121 - t84 * mrSges(3,1) + mrSges(3,3) * t129 - m(4) * (pkin(2) * t84 + t121) - (-t106 * t129 + t109 * t84) * mrSges(4,1) - (-t106 * t84 - t109 * t129) * mrSges(4,2) + t140 * (t101 + t84 * t95 + (-pkin(7) - t134) * t129) + t135 * (-t129 * t97 + t84 * t98) + t137 * (t129 * t98 + t84 * t97) + t141 * (t107 * t114 - t110 * t125)) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t107 * mrSges(2,2) - m(3) * t131 - t86 * mrSges(3,1) - mrSges(3,3) * t130 - m(4) * (pkin(2) * t86 + t131) - (t109 * t86 + t123) * mrSges(4,1) - (-t106 * t86 + t109 * t130) * mrSges(4,2) + t140 * (pkin(3) * t123 + t86 * t95 + t131) + t135 * (t130 * t97 + t86 * t98) + t137 * (-t130 * t98 + t86 * t97) + t141 * (t107 * t125 + t110 * t114)) * g(1);
U  = t1;
