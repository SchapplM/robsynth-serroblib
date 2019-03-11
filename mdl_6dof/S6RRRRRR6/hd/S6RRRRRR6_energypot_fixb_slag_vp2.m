% Calculate potential energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:46
% EndTime: 2019-03-10 04:11:46
% DurationCPUTime: 0.87s
% Computational Cost: add. (355->108), mult. (554->123), div. (0->0), fcn. (633->14), ass. (0->44)
t106 = qJ(5) + qJ(6);
t101 = cos(t106);
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t119 = -pkin(10) - pkin(9);
t143 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3);
t144 = -m(5) - m(6) - m(7);
t99 = sin(t106);
t145 = -t99 * mrSges(7,1) - t114 * mrSges(6,2) - t101 * mrSges(7,2) - t144 * t119 - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t110 - t143;
t141 = -m(6) * pkin(11) + m(7) * (-pkin(12) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t97 = pkin(5) * t114 + pkin(4);
t139 = -m(6) * pkin(4) - m(7) * t97 - mrSges(6,1) * t114 - mrSges(7,1) * t101 + mrSges(6,2) * t110 + mrSges(7,2) * t99 - mrSges(5,1);
t111 = sin(qJ(3));
t138 = pkin(3) * t111;
t109 = cos(pkin(6));
t137 = pkin(8) * t109 + pkin(7);
t117 = cos(qJ(1));
t108 = sin(pkin(6));
t113 = sin(qJ(1));
t132 = t108 * t113;
t134 = pkin(1) * t117 + pkin(8) * t132;
t112 = sin(qJ(2));
t133 = t108 * t112;
t116 = cos(qJ(2));
t131 = t108 * t116;
t130 = t108 * t117;
t129 = t112 * t113;
t128 = t112 * t117;
t127 = t113 * t116;
t126 = t116 * t117;
t125 = t110 * t131;
t124 = t111 * t132;
t115 = cos(qJ(3));
t98 = pkin(3) * t115 + pkin(2);
t123 = t109 * t138 + t119 * t131 + t133 * t98 + t137;
t104 = t113 * pkin(1);
t122 = -pkin(8) * t130 + t104;
t107 = qJ(3) + qJ(4);
t102 = cos(t107);
t100 = sin(t107);
t88 = -t109 * t129 + t126;
t86 = t109 * t128 + t127;
t82 = t100 * t109 + t102 * t133;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(5) * t123 - t82 * mrSges(5,1) + mrSges(5,3) * t131 - m(6) * (pkin(4) * t82 + t123) - (t114 * t82 - t125) * mrSges(6,1) - (-t110 * t82 - t114 * t131) * mrSges(6,2) - m(7) * (-pkin(5) * t125 + t82 * t97 + t123) - (t101 * t82 - t131 * t99) * mrSges(7,1) - (-t101 * t131 - t82 * t99) * mrSges(7,2) + (-m(3) - m(4)) * t137 + (-t111 * mrSges(4,1) - t115 * mrSges(4,2) - mrSges(3,3)) * t109 + (t143 * t116 + (-m(4) * pkin(2) - mrSges(4,1) * t115 + mrSges(4,2) * t111 - mrSges(3,1)) * t112) * t108 + t141 * (t100 * t133 - t102 * t109)) * g(3) + (-mrSges(1,2) - t113 * mrSges(2,1) - t117 * mrSges(2,2) - m(3) * t122 - t86 * mrSges(3,1) + mrSges(3,3) * t130 - m(4) * (pkin(2) * t86 + t122) - (-t111 * t130 + t115 * t86) * mrSges(4,1) - (-t111 * t86 - t115 * t130) * mrSges(4,2) + t144 * (t104 + t86 * t98 + (-pkin(8) - t138) * t130) + t139 * (-t100 * t130 + t102 * t86) + t141 * (t100 * t86 + t102 * t130) + t145 * (-t109 * t126 + t129)) * g(2) + (-mrSges(1,1) - t117 * mrSges(2,1) + t113 * mrSges(2,2) - m(3) * t134 - t88 * mrSges(3,1) - mrSges(3,3) * t132 - m(4) * (pkin(2) * t88 + t134) - (t115 * t88 + t124) * mrSges(4,1) - (-t111 * t88 + t115 * t132) * mrSges(4,2) + t144 * (pkin(3) * t124 + t88 * t98 + t134) + t139 * (t100 * t132 + t102 * t88) + t141 * (t100 * t88 - t102 * t132) + t145 * (t109 * t127 + t128)) * g(1);
U  = t1;
