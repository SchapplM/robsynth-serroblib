% Calculate potential energy for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:10
% EndTime: 2019-03-10 01:47:10
% DurationCPUTime: 0.70s
% Computational Cost: add. (370->104), mult. (600->122), div. (0->0), fcn. (697->12), ass. (0->49)
t155 = -m(6) - m(7);
t154 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t153 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3);
t152 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t151 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t150 = -mrSges(5,3) - t153;
t123 = sin(qJ(3));
t149 = pkin(3) * t123;
t121 = cos(pkin(6));
t148 = t121 * pkin(8) + pkin(7);
t120 = sin(pkin(6));
t124 = sin(qJ(2));
t147 = t120 * t124;
t125 = sin(qJ(1));
t146 = t120 * t125;
t128 = cos(qJ(2));
t145 = t120 * t128;
t129 = cos(qJ(1));
t144 = t120 * t129;
t143 = t124 * t125;
t142 = t124 * t129;
t141 = t125 * t128;
t140 = t128 * t129;
t139 = t129 * pkin(1) + pkin(8) * t146;
t138 = t123 * t146;
t117 = t125 * pkin(1);
t137 = -pkin(8) * t144 + t117;
t127 = cos(qJ(3));
t113 = pkin(3) * t127 + pkin(2);
t130 = -pkin(10) - pkin(9);
t136 = t113 * t147 + t121 * t149 + t130 * t145 + t148;
t103 = t121 * t141 + t142;
t104 = -t121 * t143 + t140;
t135 = pkin(3) * t138 - t103 * t130 + t104 * t113 + t139;
t101 = -t121 * t140 + t143;
t102 = t121 * t142 + t141;
t132 = t117 + t102 * t113 - t101 * t130 + (-pkin(8) - t149) * t144;
t126 = cos(qJ(5));
t122 = sin(qJ(5));
t119 = qJ(3) + qJ(4);
t115 = cos(t119);
t114 = sin(t119);
t96 = t114 * t121 + t115 * t147;
t95 = t114 * t147 - t121 * t115;
t91 = t104 * t115 + t114 * t146;
t90 = t104 * t114 - t115 * t146;
t89 = t102 * t115 - t114 * t144;
t88 = t102 * t114 + t115 * t144;
t1 = (-m(2) * pkin(7) - m(5) * t136 - t96 * mrSges(5,1) + mrSges(5,3) * t145 - mrSges(1,3) - mrSges(2,3) + t155 * (t96 * pkin(4) + pkin(11) * t95 + t136) + t152 * (-t122 * t145 + t126 * t96) + t151 * (t122 * t96 + t126 * t145) + (-m(3) - m(4)) * t148 + (-t123 * mrSges(4,1) - t127 * mrSges(4,2) - mrSges(3,3)) * t121 + (t153 * t128 + (-m(4) * pkin(2) - t127 * mrSges(4,1) + t123 * mrSges(4,2) - mrSges(3,1)) * t124) * t120 + t154 * t95) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t125 - mrSges(2,2) * t129 - m(3) * t137 - t102 * mrSges(3,1) + mrSges(3,3) * t144 - m(4) * (pkin(2) * t102 + t137) - (t102 * t127 - t123 * t144) * mrSges(4,1) - (-t102 * t123 - t127 * t144) * mrSges(4,2) - m(5) * t132 - t89 * mrSges(5,1) + t155 * (t89 * pkin(4) + pkin(11) * t88 + t132) + t152 * (t101 * t122 + t126 * t89) + t151 * (-t101 * t126 + t122 * t89) + t154 * t88 + t150 * t101) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t129 + mrSges(2,2) * t125 - m(3) * t139 - t104 * mrSges(3,1) - mrSges(3,3) * t146 - m(4) * (pkin(2) * t104 + t139) - (t104 * t127 + t138) * mrSges(4,1) - (-t104 * t123 + t127 * t146) * mrSges(4,2) - m(5) * t135 - t91 * mrSges(5,1) + t155 * (t91 * pkin(4) + pkin(11) * t90 + t135) + t152 * (t103 * t122 + t126 * t91) + t151 * (-t103 * t126 + t122 * t91) + t154 * t90 + t150 * t103) * g(1);
U  = t1;
