% Calculate potential energy for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:22
% EndTime: 2019-03-09 17:29:23
% DurationCPUTime: 0.66s
% Computational Cost: add. (346->100), mult. (669->118), div. (0->0), fcn. (793->12), ass. (0->51)
t160 = -m(4) - m(5);
t159 = -m(6) - m(7);
t158 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t157 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t124 = cos(pkin(11));
t156 = -t124 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t155 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t122 = sin(pkin(11));
t154 = -m(5) * pkin(3) - t124 * mrSges(5,1) + t122 * mrSges(5,2) - mrSges(4,1);
t125 = cos(pkin(6));
t153 = t125 * pkin(8) + pkin(7);
t131 = cos(qJ(2));
t132 = cos(qJ(1));
t142 = t131 * t132;
t128 = sin(qJ(2));
t129 = sin(qJ(1));
t145 = t128 * t129;
t106 = -t125 * t142 + t145;
t152 = t106 * t122;
t143 = t129 * t131;
t144 = t128 * t132;
t108 = t125 * t143 + t144;
t151 = t108 * t122;
t123 = sin(pkin(6));
t150 = t123 * t128;
t149 = t123 * t129;
t130 = cos(qJ(3));
t148 = t123 * t130;
t147 = t123 * t131;
t146 = t123 * t132;
t141 = t132 * pkin(1) + pkin(8) * t149;
t140 = pkin(2) * t150 + t153;
t139 = t129 * pkin(1) - pkin(8) * t146;
t109 = -t125 * t145 + t142;
t138 = t109 * pkin(2) + pkin(9) * t108 + t141;
t137 = -pkin(9) * t147 + t140;
t107 = t125 * t144 + t143;
t136 = t107 * pkin(2) + t106 * pkin(9) + t139;
t127 = sin(qJ(3));
t126 = -pkin(10) - qJ(4);
t121 = pkin(11) + qJ(5);
t117 = cos(t121);
t116 = sin(t121);
t115 = pkin(4) * t124 + pkin(3);
t105 = t125 * t127 + t128 * t148;
t104 = -t125 * t130 + t127 * t150;
t95 = t109 * t130 + t127 * t149;
t94 = t109 * t127 - t129 * t148;
t93 = t107 * t130 - t127 * t146;
t92 = t107 * t127 + t130 * t146;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t153 - t125 * mrSges(3,3) - (t128 * mrSges(3,1) + t131 * mrSges(3,2)) * t123 - m(4) * t137 - t105 * mrSges(4,1) + mrSges(4,3) * t147 - m(5) * (pkin(3) * t105 + t137) - (t105 * t124 - t122 * t147) * mrSges(5,1) - (-t105 * t122 - t124 * t147) * mrSges(5,2) + t159 * (-t104 * t126 + t105 * t115 + (-pkin(4) * t122 - pkin(9)) * t147 + t140) + t158 * (t105 * t117 - t116 * t147) + t157 * (t105 * t116 + t117 * t147) + t155 * t104) * g(3) + (-m(3) * t139 - t129 * mrSges(2,1) - t107 * mrSges(3,1) - t152 * mrSges(5,1) - mrSges(2,2) * t132 + mrSges(3,3) * t146 - mrSges(1,2) + t160 * t136 + t154 * t93 + t159 * (pkin(4) * t152 + t93 * t115 - t92 * t126 + t136) + t158 * (t106 * t116 + t117 * t93) + t157 * (-t106 * t117 + t116 * t93) + t156 * t106 + t155 * t92) * g(2) + (-m(3) * t141 - mrSges(2,1) * t132 - t109 * mrSges(3,1) - t151 * mrSges(5,1) + t129 * mrSges(2,2) - mrSges(3,3) * t149 - mrSges(1,1) + t160 * t138 + t154 * t95 + t159 * (pkin(4) * t151 + t95 * t115 - t94 * t126 + t138) + t158 * (t108 * t116 + t117 * t95) + t157 * (-t108 * t117 + t116 * t95) + t156 * t108 + t155 * t94) * g(1);
U  = t1;
