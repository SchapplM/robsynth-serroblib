% Calculate potential energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:23
% EndTime: 2019-03-09 15:14:24
% DurationCPUTime: 0.57s
% Computational Cost: add. (292->94), mult. (677->114), div. (0->0), fcn. (786->10), ass. (0->48)
t156 = -m(6) - m(7);
t155 = mrSges(2,2) - mrSges(3,3);
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t154 = -t123 * mrSges(3,1) + t120 * mrSges(3,2) - mrSges(2,1);
t153 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t152 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1) - mrSges(5,3);
t151 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t150 = t120 * pkin(2) + pkin(7);
t149 = cos(pkin(10));
t117 = sin(pkin(6));
t148 = qJ(4) * t117;
t147 = t117 * t123;
t119 = sin(qJ(3));
t146 = t120 * t119;
t124 = cos(qJ(1));
t145 = t120 * t124;
t121 = sin(qJ(1));
t144 = t121 * t120;
t143 = t121 * t123;
t142 = t123 * t124;
t141 = t124 * pkin(1) + t121 * pkin(8);
t140 = t117 * t146;
t118 = cos(pkin(6));
t139 = t118 * t144;
t138 = t118 * t145;
t137 = t121 * pkin(1) - pkin(8) * t124;
t136 = t118 * t149;
t135 = t120 * t149;
t134 = pkin(2) * t142 + pkin(9) * t145 + t141;
t133 = t117 * t135;
t131 = pkin(2) * t143 + pkin(9) * t144 + t137;
t122 = cos(qJ(3));
t130 = qJ(4) * t140 + t120 * t122 * pkin(3) + (-qJ(4) * t118 - pkin(9)) * t123 + t150;
t98 = -t119 * t142 + t121 * t122;
t99 = t121 * t119 + t122 * t142;
t129 = t99 * pkin(3) + qJ(4) * t138 - t98 * t148 + t134;
t96 = -t119 * t143 - t122 * t124;
t97 = -t119 * t124 + t122 * t143;
t128 = t97 * pkin(3) + qJ(4) * t139 - t96 * t148 + t131;
t116 = sin(pkin(10));
t88 = t122 * t135 + (-t118 * t146 - t147) * t116;
t87 = t149 * t147 + (t116 * t122 + t119 * t136) * t120;
t85 = t99 * t149 + (t117 * t145 + t118 * t98) * t116;
t84 = t116 * t99 - t124 * t133 - t98 * t136;
t83 = t97 * t149 + (t117 * t144 + t118 * t96) * t116;
t82 = t116 * t97 - t121 * t133 - t96 * t136;
t1 = (-m(4) * t150 - m(5) * t130 - mrSges(1,3) - mrSges(2,3) + t156 * (t88 * pkin(4) + qJ(5) * t87 + t130) + (m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3)) * t123 + (-t122 * mrSges(4,1) + t119 * mrSges(4,2) - mrSges(3,1)) * t120 + (-m(2) - m(3)) * pkin(7) + t152 * (-t118 * t123 + t140) + t151 * t88 + t153 * t87) * g(3) + (-m(3) * t137 - m(4) * t131 - m(5) * t128 - t97 * mrSges(4,1) - t96 * mrSges(4,2) - mrSges(4,3) * t144 - mrSges(1,2) + t156 * (t83 * pkin(4) + t82 * qJ(5) + t128) - t155 * t124 + t154 * t121 + t152 * (-t117 * t96 + t139) + t151 * t83 + t153 * t82) * g(2) + (-m(3) * t141 - m(4) * t134 - m(5) * t129 - t99 * mrSges(4,1) - t98 * mrSges(4,2) - mrSges(4,3) * t145 - mrSges(1,1) + t156 * (t85 * pkin(4) + qJ(5) * t84 + t129) + t154 * t124 + t155 * t121 + t152 * (-t98 * t117 + t138) + t151 * t85 + t153 * t84) * g(1);
U  = t1;
