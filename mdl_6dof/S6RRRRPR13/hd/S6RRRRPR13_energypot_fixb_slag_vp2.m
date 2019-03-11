% Calculate potential energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:03
% EndTime: 2019-03-09 23:51:03
% DurationCPUTime: 0.58s
% Computational Cost: add. (361->98), mult. (816->114), div. (0->0), fcn. (997->12), ass. (0->55)
t148 = -mrSges(4,3) + mrSges(3,2);
t147 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(10) - pkin(11)) + mrSges(7,3);
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t146 = -t111 * mrSges(7,1) - t116 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t145 = -m(7) * pkin(5) - mrSges(7,1) * t116 + mrSges(7,2) * t111 - mrSges(5,1) - mrSges(6,1);
t114 = sin(qJ(2));
t118 = cos(qJ(2));
t119 = cos(qJ(1));
t115 = sin(qJ(1));
t139 = cos(pkin(6));
t132 = t115 * t139;
t101 = -t114 * t132 + t119 * t118;
t113 = sin(qJ(3));
t110 = sin(pkin(6));
t141 = cos(qJ(3));
t133 = t110 * t141;
t89 = t101 * t113 - t115 * t133;
t144 = pkin(10) * t89;
t138 = t110 * t114;
t96 = t113 * t138 - t139 * t141;
t143 = pkin(10) * t96;
t131 = t119 * t139;
t99 = t114 * t131 + t115 * t118;
t87 = t99 * t113 + t119 * t133;
t142 = t87 * pkin(10);
t140 = t139 * pkin(8) + pkin(7);
t137 = t110 * t115;
t136 = t110 * t118;
t135 = t110 * t119;
t134 = t119 * pkin(1) + pkin(8) * t137;
t129 = t115 * pkin(1) - pkin(8) * t135;
t100 = t119 * t114 + t118 * t132;
t128 = t101 * pkin(2) + pkin(9) * t100 + t134;
t127 = pkin(2) * t138 - pkin(9) * t136 + t140;
t90 = t101 * t141 + t113 * t137;
t126 = t90 * pkin(3) + t128;
t97 = t113 * t139 + t114 * t133;
t125 = t97 * pkin(3) + t127;
t98 = t114 * t115 - t118 * t131;
t124 = t99 * pkin(2) + t98 * pkin(9) + t129;
t88 = -t113 * t135 + t141 * t99;
t123 = t88 * pkin(3) + t124;
t112 = sin(qJ(4));
t117 = cos(qJ(4));
t80 = -t100 * t117 + t112 * t90;
t81 = t100 * t112 + t117 * t90;
t122 = t81 * pkin(4) + t80 * qJ(5) + t126;
t85 = t112 * t97 + t117 * t136;
t86 = -t112 * t136 + t117 * t97;
t121 = t86 * pkin(4) + qJ(5) * t85 + t125;
t78 = t112 * t88 - t98 * t117;
t79 = t112 * t98 + t117 * t88;
t120 = t79 * pkin(4) + t78 * qJ(5) + t123;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t140 - t139 * mrSges(3,3) - (t114 * mrSges(3,1) + t118 * mrSges(3,2)) * t110 - m(4) * t127 - t97 * mrSges(4,1) + mrSges(4,3) * t136 - m(5) * (t125 + t143) - m(6) * (t121 + t143) - m(7) * t121 + t145 * t86 + t146 * t85 + t147 * t96) * g(3) + (-mrSges(1,2) - t115 * mrSges(2,1) - t119 * mrSges(2,2) - m(3) * t129 - t99 * mrSges(3,1) + mrSges(3,3) * t135 - m(4) * t124 - t88 * mrSges(4,1) - m(5) * (t123 + t142) - m(6) * (t120 + t142) - m(7) * t120 + t148 * t98 + t145 * t79 + t146 * t78 + t147 * t87) * g(2) + (-mrSges(1,1) - t119 * mrSges(2,1) + t115 * mrSges(2,2) - m(3) * t134 - t101 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * t128 - t90 * mrSges(4,1) - m(5) * (t126 + t144) - m(6) * (t122 + t144) - m(7) * t122 + t145 * t81 + t146 * t80 + t148 * t100 + t147 * t89) * g(1);
U  = t1;
