% Calculate potential energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:25
% EndTime: 2019-03-08 21:12:26
% DurationCPUTime: 0.60s
% Computational Cost: add. (361->98), mult. (816->115), div. (0->0), fcn. (997->12), ass. (0->56)
t150 = -mrSges(4,3) + mrSges(3,2);
t149 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (-pkin(9) + qJ(4)) + mrSges(7,3);
t116 = sin(qJ(6));
t119 = cos(qJ(6));
t148 = -mrSges(7,1) * t116 - mrSges(7,2) * t119 + mrSges(5,2) - mrSges(6,3);
t147 = -m(7) * pkin(5) - mrSges(7,1) * t119 + mrSges(7,2) * t116 - mrSges(5,1) - mrSges(6,1);
t146 = cos(qJ(3));
t115 = cos(pkin(10));
t117 = sin(qJ(3));
t113 = sin(pkin(6));
t134 = t113 * t146;
t112 = sin(pkin(10));
t120 = cos(qJ(2));
t118 = sin(qJ(2));
t143 = cos(pkin(6));
t133 = t118 * t143;
t98 = t112 * t120 + t115 * t133;
t88 = t115 * t134 + t98 * t117;
t145 = qJ(4) * t88;
t100 = -t112 * t133 + t115 * t120;
t90 = t100 * t117 - t112 * t134;
t144 = qJ(4) * t90;
t138 = t113 * t118;
t101 = t117 * t138 - t143 * t146;
t142 = t101 * qJ(4);
t141 = t112 * t113;
t140 = t113 * t115;
t139 = t113 * t117;
t137 = t113 * t120;
t136 = t115 * pkin(1) + pkin(7) * t141;
t135 = t143 * pkin(7) + qJ(1);
t132 = t120 * t143;
t131 = t112 * pkin(1) - pkin(7) * t140;
t99 = t112 * t132 + t115 * t118;
t129 = t100 * pkin(2) + pkin(8) * t99 + t136;
t91 = t100 * t146 + t112 * t139;
t128 = t91 * pkin(3) + t129;
t127 = pkin(2) * t138 - pkin(8) * t137 + t135;
t102 = t117 * t143 + t118 * t134;
t126 = t102 * pkin(3) + t127;
t97 = t112 * t118 - t115 * t132;
t125 = t98 * pkin(2) + pkin(8) * t97 + t131;
t89 = -t115 * t139 + t146 * t98;
t124 = t89 * pkin(3) + t125;
t111 = sin(pkin(11));
t114 = cos(pkin(11));
t81 = t111 * t91 - t99 * t114;
t82 = t111 * t99 + t114 * t91;
t123 = t82 * pkin(4) + t81 * qJ(5) + t128;
t86 = t102 * t111 + t114 * t137;
t87 = t102 * t114 - t111 * t137;
t122 = t87 * pkin(4) + t86 * qJ(5) + t126;
t79 = t111 * t89 - t97 * t114;
t80 = t111 * t97 + t114 * t89;
t121 = t80 * pkin(4) + t79 * qJ(5) + t124;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t135 - t143 * mrSges(3,3) - (t118 * mrSges(3,1) + t120 * mrSges(3,2)) * t113 - m(4) * t127 - t102 * mrSges(4,1) + mrSges(4,3) * t137 - m(5) * (t126 + t142) - m(6) * (t122 + t142) - m(7) * t122 + t147 * t87 + t148 * t86 + t149 * t101) * g(3) + (-mrSges(1,2) - t112 * mrSges(2,1) - t115 * mrSges(2,2) - m(3) * t131 - t98 * mrSges(3,1) + mrSges(3,3) * t140 - m(4) * t125 - t89 * mrSges(4,1) - m(5) * (t124 + t145) - m(6) * (t121 + t145) - m(7) * t121 + t150 * t97 + t147 * t80 + t148 * t79 + t149 * t88) * g(2) + (-mrSges(1,1) - t115 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t136 - t100 * mrSges(3,1) - mrSges(3,3) * t141 - m(4) * t129 - t91 * mrSges(4,1) - m(5) * (t128 + t144) - m(6) * (t123 + t144) - m(7) * t123 + t150 * t99 + t147 * t82 + t148 * t81 + t149 * t90) * g(1);
U  = t1;
