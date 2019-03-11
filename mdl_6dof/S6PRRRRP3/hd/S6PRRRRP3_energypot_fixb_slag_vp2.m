% Calculate potential energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:15
% EndTime: 2019-03-09 00:06:16
% DurationCPUTime: 0.68s
% Computational Cost: add. (323->102), mult. (613->120), div. (0->0), fcn. (717->12), ass. (0->48)
t146 = -m(4) - m(5);
t145 = -m(6) - m(7);
t144 = -mrSges(6,1) - mrSges(7,1);
t143 = -mrSges(6,2) - mrSges(7,2);
t112 = sin(qJ(4));
t126 = pkin(4) * t112 + pkin(8);
t117 = -pkin(10) - pkin(9);
t142 = -m(5) * pkin(9) + m(6) * t117 + m(7) * (-qJ(6) + t117) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t115 = cos(qJ(4));
t108 = qJ(4) + qJ(5);
t101 = sin(t108);
t139 = pkin(5) * t101 + t126;
t141 = -m(6) * t126 - m(7) * t139 - t112 * mrSges(5,1) - t115 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t100 = t115 * pkin(4) + pkin(3);
t102 = cos(t108);
t93 = pkin(5) * t102 + t100;
t140 = -m(5) * pkin(3) - m(6) * t100 - m(7) * t93 - t115 * mrSges(5,1) + t112 * mrSges(5,2) - mrSges(4,1);
t138 = cos(qJ(3));
t111 = cos(pkin(11));
t109 = sin(pkin(11));
t110 = sin(pkin(6));
t134 = t109 * t110;
t136 = t111 * pkin(1) + pkin(7) * t134;
t135 = cos(pkin(6));
t133 = t110 * t111;
t113 = sin(qJ(3));
t132 = t110 * t113;
t114 = sin(qJ(2));
t131 = t110 * t114;
t116 = cos(qJ(2));
t130 = t110 * t116;
t129 = t135 * pkin(7) + qJ(1);
t124 = t114 * t135;
t90 = -t109 * t124 + t111 * t116;
t128 = t90 * pkin(2) + t136;
t127 = pkin(2) * t131 + t129;
t125 = t110 * t138;
t123 = t116 * t135;
t122 = t109 * pkin(1) - pkin(7) * t133;
t88 = t109 * t116 + t111 * t124;
t120 = t88 * pkin(2) + t122;
t119 = -pkin(8) * t130 + t127;
t92 = t113 * t135 + t114 * t125;
t89 = t109 * t123 + t111 * t114;
t87 = t109 * t114 - t111 * t123;
t84 = t109 * t132 + t138 * t90;
t82 = -t111 * t132 + t138 * t88;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t129 - t135 * mrSges(3,3) - (t114 * mrSges(3,1) + t116 * mrSges(3,2)) * t110 - m(4) * t119 - t92 * mrSges(4,1) + mrSges(4,3) * t130 - m(5) * (pkin(3) * t92 + t119) - (-t112 * t130 + t115 * t92) * mrSges(5,1) - (-t112 * t92 - t115 * t130) * mrSges(5,2) - m(6) * (t92 * t100 - t126 * t130 + t127) - m(7) * (-t139 * t130 + t92 * t93 + t127) + t144 * (-t101 * t130 + t102 * t92) + t143 * (-t101 * t92 - t102 * t130) + t142 * (t113 * t131 - t135 * t138)) * g(3) + (-m(3) * t122 - mrSges(2,1) * t109 - t88 * mrSges(3,1) - mrSges(2,2) * t111 + mrSges(3,3) * t133 - mrSges(1,2) + t145 * t120 + t146 * (t87 * pkin(8) + t120) + t140 * t82 + t141 * t87 + t144 * (t101 * t87 + t102 * t82) + t143 * (-t101 * t82 + t102 * t87) + t142 * (t111 * t125 + t88 * t113)) * g(2) + (-m(3) * t136 - mrSges(2,1) * t111 - t90 * mrSges(3,1) + mrSges(2,2) * t109 - mrSges(3,3) * t134 - mrSges(1,1) + t145 * t128 + t146 * (t89 * pkin(8) + t128) + t140 * t84 + t141 * t89 + t144 * (t101 * t89 + t102 * t84) + t143 * (-t101 * t84 + t102 * t89) + t142 * (-t109 * t125 + t113 * t90)) * g(1);
U  = t1;
