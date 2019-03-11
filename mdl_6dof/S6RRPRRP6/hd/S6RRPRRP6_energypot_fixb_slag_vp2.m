% Calculate potential energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:39
% EndTime: 2019-03-09 12:07:40
% DurationCPUTime: 0.62s
% Computational Cost: add. (418->98), mult. (947->125), div. (0->0), fcn. (1177->12), ass. (0->54)
t166 = -m(6) - m(7);
t165 = mrSges(4,2) - mrSges(5,3);
t164 = -mrSges(3,3) - mrSges(4,3);
t163 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t162 = -m(3) * pkin(1) - mrSges(2,1);
t161 = m(3) * pkin(8) - t164;
t160 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t159 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t130 = cos(pkin(6));
t158 = t130 * pkin(8) + pkin(7);
t128 = sin(pkin(6));
t133 = sin(qJ(2));
t157 = t128 * t133;
t134 = sin(qJ(1));
t156 = t128 * t134;
t138 = cos(qJ(1));
t155 = t128 * t138;
t129 = cos(pkin(11));
t137 = cos(qJ(2));
t154 = t129 * t137;
t153 = t133 * t138;
t152 = t134 * t133;
t151 = t134 * t137;
t150 = t137 * t138;
t115 = pkin(2) * t130 * t133 + (-pkin(8) - qJ(3)) * t128;
t124 = pkin(2) * t137 + pkin(1);
t149 = t138 * t115 + t134 * t124;
t148 = -t115 * t134 + t138 * t124;
t147 = pkin(2) * t157 + t130 * qJ(3) + t158;
t127 = sin(pkin(11));
t146 = t127 * t137 + t129 * t133;
t117 = -t127 * t133 + t154;
t144 = t117 * t130;
t101 = -t134 * t146 + t138 * t144;
t114 = t146 * t130;
t102 = t114 * t138 + t134 * t117;
t145 = t102 * pkin(3) - pkin(9) * t101 + t149;
t103 = -t134 * t144 - t138 * t146;
t104 = -t134 * t114 + t117 * t138;
t143 = t104 * pkin(3) - pkin(9) * t103 + t148;
t112 = t127 * t157 - t128 * t154;
t113 = t146 * t128;
t142 = t113 * pkin(3) + pkin(9) * t112 + t147;
t136 = cos(qJ(4));
t135 = cos(qJ(5));
t132 = sin(qJ(4));
t131 = sin(qJ(5));
t107 = t113 * t136 + t130 * t132;
t106 = t113 * t132 - t130 * t136;
t96 = t104 * t136 + t132 * t156;
t95 = t104 * t132 - t136 * t156;
t94 = t102 * t136 - t132 * t155;
t93 = t102 * t132 + t136 * t155;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t158 - (mrSges(3,1) * t133 + mrSges(3,2) * t137) * t128 - m(4) * t147 - t113 * mrSges(4,1) - m(5) * t142 - t107 * mrSges(5,1) + t166 * (t107 * pkin(4) + pkin(10) * t106 + t142) + t160 * (t107 * t135 + t112 * t131) + t159 * (t107 * t131 - t112 * t135) + t164 * t130 + t165 * t112 + t163 * t106) * g(3) + (-mrSges(1,2) - t138 * mrSges(2,2) - (t130 * t153 + t151) * mrSges(3,1) - (t130 * t150 - t152) * mrSges(3,2) - m(4) * t149 - t102 * mrSges(4,1) - m(5) * t145 - t94 * mrSges(5,1) + t166 * (t94 * pkin(4) + pkin(10) * t93 + t145) + t160 * (-t101 * t131 + t135 * t94) + t159 * (t101 * t135 + t131 * t94) + t162 * t134 + t161 * t155 - t165 * t101 + t163 * t93) * g(2) + (-mrSges(1,1) + t134 * mrSges(2,2) - (-t130 * t152 + t150) * mrSges(3,1) - (-t130 * t151 - t153) * mrSges(3,2) - m(4) * t148 - t104 * mrSges(4,1) - m(5) * t143 - t96 * mrSges(5,1) + t166 * (t96 * pkin(4) + pkin(10) * t95 + t143) + t160 * (-t103 * t131 + t135 * t96) + t159 * (t103 * t135 + t131 * t96) + t162 * t138 - t161 * t156 - t165 * t103 + t163 * t95) * g(1);
U  = t1;
