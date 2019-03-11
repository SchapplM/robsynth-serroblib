% Calculate potential energy for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:04
% EndTime: 2019-03-09 00:01:04
% DurationCPUTime: 0.70s
% Computational Cost: add. (370->104), mult. (600->126), div. (0->0), fcn. (697->12), ass. (0->49)
t155 = -m(6) - m(7);
t154 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t153 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t152 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t151 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t150 = -mrSges(5,3) - t153;
t125 = sin(qJ(3));
t149 = pkin(3) * t125;
t120 = sin(pkin(11));
t121 = sin(pkin(6));
t148 = t120 * t121;
t122 = cos(pkin(11));
t147 = t121 * t122;
t146 = t121 * t125;
t126 = sin(qJ(2));
t145 = t121 * t126;
t128 = cos(qJ(3));
t144 = t121 * t128;
t129 = cos(qJ(2));
t143 = t121 * t129;
t123 = cos(pkin(6));
t142 = t123 * t126;
t141 = t123 * t129;
t140 = t122 * pkin(1) + pkin(7) * t148;
t139 = t123 * pkin(7) + qJ(1);
t138 = t120 * t146;
t116 = t120 * pkin(1);
t137 = -pkin(7) * t147 + t116;
t103 = t120 * t141 + t122 * t126;
t104 = -t120 * t142 + t122 * t129;
t113 = pkin(3) * t128 + pkin(2);
t130 = -pkin(9) - pkin(8);
t136 = pkin(3) * t138 - t103 * t130 + t104 * t113 + t140;
t135 = t113 * t145 + t123 * t149 + t130 * t143 + t139;
t101 = t120 * t126 - t122 * t141;
t102 = t120 * t129 + t122 * t142;
t132 = t116 + t102 * t113 - t101 * t130 + (-pkin(7) - t149) * t147;
t127 = cos(qJ(5));
t124 = sin(qJ(5));
t119 = qJ(3) + qJ(4);
t115 = cos(t119);
t114 = sin(t119);
t96 = t114 * t123 + t115 * t145;
t95 = t114 * t145 - t123 * t115;
t89 = t104 * t115 + t114 * t148;
t88 = t104 * t114 - t115 * t148;
t87 = t102 * t115 - t114 * t147;
t86 = t102 * t114 + t115 * t147;
t1 = (-m(2) * qJ(1) - m(5) * t135 - t96 * mrSges(5,1) + mrSges(5,3) * t143 - mrSges(1,3) - mrSges(2,3) + t155 * (t96 * pkin(4) + pkin(10) * t95 + t135) + t152 * (-t124 * t143 + t127 * t96) + t151 * (t124 * t96 + t127 * t143) + (-m(3) - m(4)) * t139 + (-t125 * mrSges(4,1) - t128 * mrSges(4,2) - mrSges(3,3)) * t123 + (t153 * t129 + (-m(4) * pkin(2) - t128 * mrSges(4,1) + t125 * mrSges(4,2) - mrSges(3,1)) * t126) * t121 + t154 * t95) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t120 - mrSges(2,2) * t122 - m(3) * t137 - t102 * mrSges(3,1) + mrSges(3,3) * t147 - m(4) * (pkin(2) * t102 + t137) - (t102 * t128 - t122 * t146) * mrSges(4,1) - (-t102 * t125 - t122 * t144) * mrSges(4,2) - m(5) * t132 - t87 * mrSges(5,1) + t155 * (t87 * pkin(4) + pkin(10) * t86 + t132) + t152 * (t101 * t124 + t127 * t87) + t151 * (-t101 * t127 + t124 * t87) + t154 * t86 + t150 * t101) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t122 + mrSges(2,2) * t120 - m(3) * t140 - t104 * mrSges(3,1) - mrSges(3,3) * t148 - m(4) * (pkin(2) * t104 + t140) - (t104 * t128 + t138) * mrSges(4,1) - (-t104 * t125 + t120 * t144) * mrSges(4,2) - m(5) * t136 - t89 * mrSges(5,1) + t155 * (t89 * pkin(4) + pkin(10) * t88 + t136) + t152 * (t103 * t124 + t127 * t89) + t151 * (-t103 * t127 + t124 * t89) + t154 * t88 + t150 * t103) * g(1);
U  = t1;
