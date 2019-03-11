% Calculate potential energy for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:21
% EndTime: 2019-03-10 02:00:22
% DurationCPUTime: 0.64s
% Computational Cost: add. (323->102), mult. (613->119), div. (0->0), fcn. (717->12), ass. (0->47)
t144 = -m(4) - m(5);
t143 = -m(6) - m(7);
t142 = -mrSges(6,1) - mrSges(7,1);
t141 = -mrSges(6,2) - mrSges(7,2);
t109 = sin(qJ(4));
t125 = pkin(4) * t109 + pkin(9);
t116 = -pkin(11) - pkin(10);
t140 = -m(5) * pkin(10) + m(6) * t116 + m(7) * (-qJ(6) + t116) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t113 = cos(qJ(4));
t107 = qJ(4) + qJ(5);
t100 = sin(t107);
t137 = pkin(5) * t100 + t125;
t139 = -m(6) * t125 - m(7) * t137 - t109 * mrSges(5,1) - t113 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t101 = cos(t107);
t99 = t113 * pkin(4) + pkin(3);
t92 = pkin(5) * t101 + t99;
t138 = -m(5) * pkin(3) - m(6) * t99 - m(7) * t92 - t113 * mrSges(5,1) + t109 * mrSges(5,2) - mrSges(4,1);
t136 = cos(qJ(3));
t132 = cos(pkin(6));
t134 = t132 * pkin(8) + pkin(7);
t115 = cos(qJ(1));
t108 = sin(pkin(6));
t112 = sin(qJ(1));
t130 = t108 * t112;
t133 = t115 * pkin(1) + pkin(8) * t130;
t111 = sin(qJ(2));
t131 = t108 * t111;
t114 = cos(qJ(2));
t129 = t108 * t114;
t128 = t108 * t115;
t127 = pkin(2) * t131 + t134;
t123 = t112 * t132;
t91 = -t111 * t123 + t115 * t114;
t126 = t91 * pkin(2) + t133;
t124 = t108 * t136;
t122 = t115 * t132;
t121 = t112 * pkin(1) - pkin(8) * t128;
t89 = t111 * t122 + t112 * t114;
t119 = t89 * pkin(2) + t121;
t118 = -pkin(9) * t129 + t127;
t110 = sin(qJ(3));
t90 = t115 * t111 + t114 * t123;
t88 = t111 * t112 - t114 * t122;
t87 = t110 * t132 + t111 * t124;
t83 = t110 * t130 + t136 * t91;
t81 = -t110 * t128 + t136 * t89;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t134 - t132 * mrSges(3,3) - (t111 * mrSges(3,1) + t114 * mrSges(3,2)) * t108 - m(4) * t118 - t87 * mrSges(4,1) + mrSges(4,3) * t129 - m(5) * (pkin(3) * t87 + t118) - (-t109 * t129 + t113 * t87) * mrSges(5,1) - (-t109 * t87 - t113 * t129) * mrSges(5,2) - m(6) * (-t125 * t129 + t87 * t99 + t127) - m(7) * (-t129 * t137 + t87 * t92 + t127) + t142 * (-t100 * t129 + t101 * t87) + t141 * (-t100 * t87 - t101 * t129) + t140 * (t110 * t131 - t132 * t136)) * g(3) + (-m(3) * t121 - mrSges(2,1) * t112 - t89 * mrSges(3,1) - mrSges(2,2) * t115 + mrSges(3,3) * t128 - mrSges(1,2) + t143 * t119 + t144 * (t88 * pkin(9) + t119) + t138 * t81 + t139 * t88 + t142 * (t100 * t88 + t101 * t81) + t141 * (-t100 * t81 + t101 * t88) + t140 * (t89 * t110 + t115 * t124)) * g(2) + (-m(3) * t133 - mrSges(2,1) * t115 - t91 * mrSges(3,1) + mrSges(2,2) * t112 - mrSges(3,3) * t130 - mrSges(1,1) + t143 * t126 + t144 * (t90 * pkin(9) + t126) + t138 * t83 + t139 * t90 + t142 * (t100 * t90 + t101 * t83) + t141 * (-t100 * t83 + t101 * t90) + t140 * (t110 * t91 - t112 * t124)) * g(1);
U  = t1;
