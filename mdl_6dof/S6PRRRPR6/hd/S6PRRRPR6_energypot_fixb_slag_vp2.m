% Calculate potential energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:03
% EndTime: 2019-03-08 23:31:03
% DurationCPUTime: 0.59s
% Computational Cost: add. (361->98), mult. (816->115), div. (0->0), fcn. (997->12), ass. (0->56)
t150 = -mrSges(4,3) + mrSges(3,2);
t149 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(9) - pkin(10)) + mrSges(7,3);
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t148 = -t114 * mrSges(7,1) - t118 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t147 = -m(7) * pkin(5) - mrSges(7,1) * t118 + mrSges(7,2) * t114 - mrSges(5,1) - mrSges(6,1);
t113 = cos(pkin(11));
t116 = sin(qJ(3));
t112 = sin(pkin(6));
t144 = cos(qJ(3));
t134 = t112 * t144;
t111 = sin(pkin(11));
t120 = cos(qJ(2));
t117 = sin(qJ(2));
t142 = cos(pkin(6));
t133 = t117 * t142;
t98 = t111 * t120 + t113 * t133;
t86 = t113 * t134 + t98 * t116;
t146 = pkin(9) * t86;
t100 = -t111 * t133 + t113 * t120;
t88 = t100 * t116 - t111 * t134;
t145 = pkin(9) * t88;
t138 = t112 * t117;
t101 = t116 * t138 - t142 * t144;
t143 = t101 * pkin(9);
t141 = t111 * t112;
t140 = t112 * t113;
t139 = t112 * t116;
t137 = t112 * t120;
t136 = t113 * pkin(1) + pkin(7) * t141;
t135 = t142 * pkin(7) + qJ(1);
t132 = t120 * t142;
t130 = t111 * pkin(1) - pkin(7) * t140;
t99 = t111 * t132 + t113 * t117;
t129 = t100 * pkin(2) + pkin(8) * t99 + t136;
t89 = t100 * t144 + t111 * t139;
t128 = t89 * pkin(3) + t129;
t127 = pkin(2) * t138 - pkin(8) * t137 + t135;
t102 = t116 * t142 + t117 * t134;
t126 = t102 * pkin(3) + t127;
t97 = t111 * t117 - t113 * t132;
t125 = t98 * pkin(2) + pkin(8) * t97 + t130;
t87 = -t113 * t139 + t144 * t98;
t124 = t87 * pkin(3) + t125;
t115 = sin(qJ(4));
t119 = cos(qJ(4));
t81 = t115 * t89 - t99 * t119;
t82 = t115 * t99 + t119 * t89;
t123 = t82 * pkin(4) + t81 * qJ(5) + t128;
t90 = t102 * t115 + t119 * t137;
t91 = t102 * t119 - t115 * t137;
t122 = t91 * pkin(4) + t90 * qJ(5) + t126;
t79 = t115 * t87 - t97 * t119;
t80 = t115 * t97 + t119 * t87;
t121 = t80 * pkin(4) + t79 * qJ(5) + t124;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t135 - t142 * mrSges(3,3) - (t117 * mrSges(3,1) + t120 * mrSges(3,2)) * t112 - m(4) * t127 - t102 * mrSges(4,1) + mrSges(4,3) * t137 - m(5) * (t126 + t143) - m(6) * (t122 + t143) - m(7) * t122 + t147 * t91 + t148 * t90 + t149 * t101) * g(3) + (-mrSges(1,2) - t111 * mrSges(2,1) - t113 * mrSges(2,2) - m(3) * t130 - t98 * mrSges(3,1) + mrSges(3,3) * t140 - m(4) * t125 - t87 * mrSges(4,1) - m(5) * (t124 + t146) - m(6) * (t121 + t146) - m(7) * t121 + t150 * t97 + t147 * t80 + t148 * t79 + t149 * t86) * g(2) + (-mrSges(1,1) - t113 * mrSges(2,1) + t111 * mrSges(2,2) - m(3) * t136 - t100 * mrSges(3,1) - mrSges(3,3) * t141 - m(4) * t129 - t89 * mrSges(4,1) - m(5) * (t128 + t145) - m(6) * (t123 + t145) - m(7) * t123 + t150 * t99 + t147 * t82 + t148 * t81 + t149 * t88) * g(1);
U  = t1;
