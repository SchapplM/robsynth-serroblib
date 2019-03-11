% Calculate potential energy for
% S6RRPRRP5
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:41
% EndTime: 2019-03-09 11:57:42
% DurationCPUTime: 0.72s
% Computational Cost: add. (388->104), mult. (862->129), div. (0->0), fcn. (1059->12), ass. (0->53)
t158 = -mrSges(6,1) - mrSges(7,1);
t157 = -mrSges(6,2) - mrSges(7,2);
t156 = -mrSges(3,3) - mrSges(4,3);
t118 = sin(pkin(11));
t124 = sin(qJ(2));
t128 = cos(qJ(2));
t149 = cos(pkin(11));
t107 = -t124 * t118 + t128 * t149;
t155 = -m(3) * pkin(1) - mrSges(2,1);
t154 = m(3) * pkin(8) - t156;
t122 = sin(qJ(5));
t153 = m(7) * (pkin(5) * t122 + pkin(9)) - mrSges(4,2) + mrSges(5,3);
t152 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t151 = pkin(2) * t124;
t120 = cos(pkin(6));
t150 = t120 * pkin(8) + pkin(7);
t119 = sin(pkin(6));
t125 = sin(qJ(1));
t148 = t119 * t125;
t129 = cos(qJ(1));
t147 = t119 * t129;
t145 = t124 * t129;
t144 = t125 * t124;
t143 = t125 * t128;
t142 = t128 * t129;
t105 = t120 * t151 + (-pkin(8) - qJ(3)) * t119;
t115 = pkin(2) * t128 + pkin(1);
t141 = t129 * t105 + t125 * t115;
t106 = -t128 * t118 - t124 * t149;
t104 = t106 * t120;
t94 = -t104 * t129 + t125 * t107;
t140 = t94 * pkin(3) + t141;
t137 = -t105 * t125 + t129 * t115;
t136 = t120 * qJ(3) + t119 * t151 + t150;
t96 = t125 * t104 + t107 * t129;
t135 = t96 * pkin(3) + t137;
t103 = t106 * t119;
t134 = -t103 * pkin(3) + t136;
t130 = t120 * t107;
t93 = t125 * t106 + t129 * t130;
t133 = -pkin(9) * t93 + t140;
t95 = t106 * t129 - t125 * t130;
t132 = -pkin(9) * t95 + t135;
t102 = t107 * t119;
t131 = -pkin(9) * t102 + t134;
t127 = cos(qJ(4));
t126 = cos(qJ(5));
t123 = sin(qJ(4));
t114 = pkin(5) * t126 + pkin(4);
t98 = -t103 * t127 + t120 * t123;
t90 = t123 * t148 + t127 * t96;
t88 = -t123 * t147 + t94 * t127;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t150 - (t124 * mrSges(3,1) + t128 * mrSges(3,2)) * t119 - m(4) * t136 + t103 * mrSges(4,1) - m(5) * t131 - t98 * mrSges(5,1) - m(6) * (pkin(4) * t98 + t131) - m(7) * (t114 * t98 + t134) + t158 * (-t102 * t122 + t126 * t98) + t157 * (-t102 * t126 - t122 * t98) + t156 * t120 + t153 * t102 + t152 * (-t103 * t123 - t120 * t127)) * g(3) + (-mrSges(1,2) - t129 * mrSges(2,2) - (t120 * t145 + t143) * mrSges(3,1) - (t120 * t142 - t144) * mrSges(3,2) - m(4) * t141 - t94 * mrSges(4,1) - m(5) * t133 - t88 * mrSges(5,1) - m(6) * (pkin(4) * t88 + t133) - m(7) * (t114 * t88 + t140) + t153 * t93 + t158 * (-t122 * t93 + t126 * t88) + t157 * (-t122 * t88 - t126 * t93) + t155 * t125 + t154 * t147 + t152 * (t94 * t123 + t127 * t147)) * g(2) + (-mrSges(1,1) + t125 * mrSges(2,2) - (-t120 * t144 + t142) * mrSges(3,1) - (-t120 * t143 - t145) * mrSges(3,2) - m(4) * t137 - t96 * mrSges(4,1) - m(5) * t132 - t90 * mrSges(5,1) - m(6) * (pkin(4) * t90 + t132) - m(7) * (t114 * t90 + t135) + t153 * t95 + t158 * (-t122 * t95 + t126 * t90) + t157 * (-t122 * t90 - t126 * t95) + t155 * t129 - t154 * t148 + t152 * (t123 * t96 - t127 * t148)) * g(1);
U  = t1;
