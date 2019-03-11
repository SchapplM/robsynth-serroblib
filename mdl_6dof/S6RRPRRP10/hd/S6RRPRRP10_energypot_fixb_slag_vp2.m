% Calculate potential energy for
% S6RRPRRP10
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:26
% EndTime: 2019-03-09 12:36:27
% DurationCPUTime: 0.70s
% Computational Cost: add. (370->104), mult. (600->122), div. (0->0), fcn. (697->12), ass. (0->49)
t155 = -m(6) - m(7);
t154 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t153 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t152 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t151 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t150 = -mrSges(5,3) - t152;
t120 = sin(pkin(11));
t149 = pkin(3) * t120;
t123 = cos(pkin(6));
t148 = t123 * pkin(8) + pkin(7);
t121 = sin(pkin(6));
t126 = sin(qJ(2));
t147 = t121 * t126;
t127 = sin(qJ(1));
t146 = t121 * t127;
t129 = cos(qJ(2));
t145 = t121 * t129;
t130 = cos(qJ(1));
t144 = t121 * t130;
t143 = t126 * t127;
t142 = t126 * t130;
t141 = t127 * t129;
t140 = t129 * t130;
t139 = t130 * pkin(1) + pkin(8) * t146;
t138 = t120 * t146;
t117 = t127 * pkin(1);
t137 = -pkin(8) * t144 + t117;
t122 = cos(pkin(11));
t113 = pkin(3) * t122 + pkin(2);
t124 = -pkin(9) - qJ(3);
t136 = t113 * t147 + t123 * t149 + t124 * t145 + t148;
t103 = t123 * t141 + t142;
t104 = -t123 * t143 + t140;
t135 = pkin(3) * t138 - t103 * t124 + t104 * t113 + t139;
t101 = -t123 * t140 + t143;
t102 = t123 * t142 + t141;
t132 = t117 + t102 * t113 - t101 * t124 + (-pkin(8) - t149) * t144;
t128 = cos(qJ(5));
t125 = sin(qJ(5));
t119 = pkin(11) + qJ(4);
t115 = cos(t119);
t114 = sin(t119);
t96 = t114 * t123 + t115 * t147;
t95 = t114 * t147 - t123 * t115;
t91 = t104 * t115 + t114 * t146;
t90 = t104 * t114 - t115 * t146;
t89 = t102 * t115 - t114 * t144;
t88 = t102 * t114 + t115 * t144;
t1 = (-m(2) * pkin(7) - m(5) * t136 - t96 * mrSges(5,1) + mrSges(5,3) * t145 - mrSges(1,3) - mrSges(2,3) + t155 * (t96 * pkin(4) + t95 * pkin(10) + t136) + t153 * (-t125 * t145 + t128 * t96) + t151 * (t125 * t96 + t128 * t145) + (-m(3) - m(4)) * t148 + (-t120 * mrSges(4,1) - t122 * mrSges(4,2) - mrSges(3,3)) * t123 + (t152 * t129 + (-m(4) * pkin(2) - t122 * mrSges(4,1) + t120 * mrSges(4,2) - mrSges(3,1)) * t126) * t121 + t154 * t95) * g(3) + (-mrSges(1,2) - t127 * mrSges(2,1) - t130 * mrSges(2,2) - m(3) * t137 - t102 * mrSges(3,1) + mrSges(3,3) * t144 - m(4) * (t102 * pkin(2) + t137) - (t102 * t122 - t120 * t144) * mrSges(4,1) - (-t102 * t120 - t122 * t144) * mrSges(4,2) - m(5) * t132 - t89 * mrSges(5,1) + t155 * (t89 * pkin(4) + t88 * pkin(10) + t132) + t153 * (t101 * t125 + t128 * t89) + t151 * (-t101 * t128 + t125 * t89) + t154 * t88 + t150 * t101) * g(2) + (-mrSges(1,1) - t130 * mrSges(2,1) + t127 * mrSges(2,2) - m(3) * t139 - t104 * mrSges(3,1) - mrSges(3,3) * t146 - m(4) * (pkin(2) * t104 + t139) - (t104 * t122 + t138) * mrSges(4,1) - (-t104 * t120 + t122 * t146) * mrSges(4,2) - m(5) * t135 - t91 * mrSges(5,1) + t155 * (t91 * pkin(4) + pkin(10) * t90 + t135) + t153 * (t103 * t125 + t128 * t91) + t151 * (-t103 * t128 + t125 * t91) + t154 * t90 + t150 * t103) * g(1);
U  = t1;
