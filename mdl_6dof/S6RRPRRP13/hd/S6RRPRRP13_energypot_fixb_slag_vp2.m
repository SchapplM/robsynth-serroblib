% Calculate potential energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:28
% EndTime: 2019-03-09 12:55:29
% DurationCPUTime: 0.67s
% Computational Cost: add. (260->101), mult. (545->114), div. (0->0), fcn. (622->10), ass. (0->48)
t144 = -mrSges(3,1) + mrSges(4,2);
t143 = -mrSges(6,1) - mrSges(7,1);
t142 = -mrSges(6,2) - mrSges(7,2);
t141 = mrSges(3,3) + mrSges(4,1);
t140 = -mrSges(4,3) + mrSges(3,2);
t108 = sin(qJ(5));
t139 = -m(7) * (pkin(5) * t108 + pkin(9)) - mrSges(5,3) + t144;
t138 = m(6) * pkin(10) - m(7) * (-qJ(6) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t106 = cos(pkin(6));
t137 = t106 * pkin(8) + pkin(7);
t115 = cos(qJ(1));
t105 = sin(pkin(6));
t111 = sin(qJ(1));
t134 = t105 * t111;
t136 = t115 * pkin(1) + pkin(8) * t134;
t110 = sin(qJ(2));
t135 = t105 * t110;
t114 = cos(qJ(2));
t133 = t105 * t114;
t132 = t105 * t115;
t131 = t108 * t110;
t130 = t110 * t111;
t129 = t110 * t115;
t128 = t111 * t114;
t127 = t114 * t115;
t126 = pkin(2) * t135 + t137;
t125 = pkin(8) * t132;
t103 = t111 * pkin(1);
t90 = -t106 * t127 + t130;
t91 = t106 * t129 + t128;
t123 = t91 * pkin(2) + t90 * qJ(3) + t103;
t122 = t106 * pkin(3) + pkin(9) * t135 + t126;
t92 = t106 * t128 + t129;
t93 = -t106 * t130 + t127;
t121 = t93 * pkin(2) + qJ(3) * t92 + t136;
t120 = pkin(3) * t134 + t121;
t119 = -qJ(3) * t133 + t122;
t118 = pkin(9) * t93 + t120;
t117 = (-pkin(3) - pkin(8)) * t132 + t123;
t116 = t91 * pkin(9) + t117;
t113 = cos(qJ(4));
t112 = cos(qJ(5));
t109 = sin(qJ(4));
t100 = pkin(5) * t112 + pkin(4);
t89 = t106 * t113 - t109 * t133;
t84 = t90 * t109 - t113 * t132;
t82 = t109 * t92 + t113 * t134;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t137 - m(4) * t126 - m(5) * t119 - t89 * mrSges(5,1) - mrSges(5,3) * t135 - m(6) * (pkin(4) * t89 + t119) - m(7) * (t100 * t89 + t122) + t143 * (t105 * t131 + t112 * t89) + t142 * (-t108 * t89 + t112 * t135) - t141 * t106 + (-m(7) * pkin(5) * t131 + ((m(4) + m(7)) * qJ(3) - t140) * t114 + t144 * t110) * t105 - t138 * (t106 * t109 + t113 * t133)) * g(3) + (-mrSges(1,2) - t111 * mrSges(2,1) - mrSges(2,2) * t115 - m(3) * (t103 - t125) - m(4) * (t123 - t125) - m(5) * t116 - t84 * mrSges(5,1) - m(6) * (t84 * pkin(4) + t116) - m(7) * (t84 * t100 + t117) + t140 * t90 + t143 * (t108 * t91 + t112 * t84) + t142 * (-t108 * t84 + t112 * t91) + t141 * t132 + t139 * t91 + t138 * (t109 * t132 + t90 * t113)) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t115 + t111 * mrSges(2,2) - m(3) * t136 - m(4) * t121 - m(5) * t118 - t82 * mrSges(5,1) - m(6) * (pkin(4) * t82 + t118) - m(7) * (t100 * t82 + t120) + t140 * t92 + t143 * (t108 * t93 + t112 * t82) + t142 * (-t108 * t82 + t112 * t93) - t141 * t134 + t139 * t93 - t138 * (t109 * t134 - t92 * t113)) * g(1);
U  = t1;
