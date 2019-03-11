% Calculate potential energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:02
% EndTime: 2019-03-09 10:38:02
% DurationCPUTime: 0.67s
% Computational Cost: add. (372->96), mult. (827->112), div. (0->0), fcn. (1011->12), ass. (0->57)
t154 = -mrSges(3,3) - mrSges(4,3);
t113 = sin(pkin(11));
t118 = sin(qJ(2));
t122 = cos(qJ(2));
t141 = cos(pkin(11));
t102 = -t118 * t113 + t122 * t141;
t153 = -m(3) * pkin(1) - mrSges(2,1);
t152 = m(3) * pkin(8) - t154;
t151 = -m(7) * pkin(10) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t150 = -t116 * mrSges(7,1) - t120 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t149 = m(7) * (pkin(5) + pkin(9)) + t120 * mrSges(7,1) - t116 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t101 = -t122 * t113 - t118 * t141;
t119 = sin(qJ(1));
t123 = cos(qJ(1));
t115 = cos(pkin(6));
t124 = t115 * t102;
t87 = t119 * t101 + t123 * t124;
t147 = pkin(9) * t87;
t89 = t101 * t123 - t119 * t124;
t146 = pkin(9) * t89;
t114 = sin(pkin(6));
t97 = t102 * t114;
t145 = pkin(9) * t97;
t144 = pkin(2) * t118;
t143 = t115 * pkin(8) + pkin(7);
t100 = t115 * t144 + (-pkin(8) - qJ(3)) * t114;
t110 = pkin(2) * t122 + pkin(1);
t142 = t123 * t100 + t119 * t110;
t140 = t114 * t119;
t139 = t114 * t123;
t137 = t118 * t123;
t136 = t119 * t118;
t135 = t119 * t122;
t134 = t122 * t123;
t99 = t101 * t115;
t88 = t119 * t102 - t123 * t99;
t133 = t88 * pkin(3) + t142;
t131 = -t100 * t119 + t123 * t110;
t130 = t115 * qJ(3) + t114 * t144 + t143;
t90 = t102 * t123 + t119 * t99;
t129 = t90 * pkin(3) + t131;
t98 = t101 * t114;
t128 = -t98 * pkin(3) + t130;
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t81 = t88 * t117 + t121 * t139;
t82 = -t117 * t139 + t121 * t88;
t127 = t82 * pkin(4) + t81 * qJ(5) + t133;
t83 = t117 * t90 - t121 * t140;
t84 = t117 * t140 + t121 * t90;
t126 = t84 * pkin(4) + t83 * qJ(5) + t129;
t92 = -t115 * t121 - t117 * t98;
t93 = t115 * t117 - t121 * t98;
t125 = t93 * pkin(4) + qJ(5) * t92 + t128;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t143 - (t118 * mrSges(3,1) + t122 * mrSges(3,2)) * t114 - m(4) * t130 + t98 * mrSges(4,1) - m(5) * (t128 - t145) - m(6) * (t125 - t145) - m(7) * t125 + t154 * t115 + t150 * t92 + t149 * t97 + t151 * t93) * g(3) + (-mrSges(1,2) - mrSges(2,2) * t123 - (t115 * t137 + t135) * mrSges(3,1) - (t115 * t134 - t136) * mrSges(3,2) - m(4) * t142 - t88 * mrSges(4,1) - m(5) * (t133 - t147) - m(6) * (t127 - t147) - m(7) * t127 + t153 * t119 + t152 * t139 + t150 * t81 + t149 * t87 + t151 * t82) * g(2) + (-mrSges(1,1) + t119 * mrSges(2,2) - (-t115 * t136 + t134) * mrSges(3,1) - (-t115 * t135 - t137) * mrSges(3,2) - m(4) * t131 - t90 * mrSges(4,1) - m(5) * (t129 - t146) - m(6) * (t126 - t146) - m(7) * t126 + t153 * t123 - t152 * t140 + t150 * t83 + t149 * t89 + t151 * t84) * g(1);
U  = t1;
