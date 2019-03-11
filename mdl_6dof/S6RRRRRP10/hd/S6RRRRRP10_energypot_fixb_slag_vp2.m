% Calculate potential energy for
% S6RRRRRP10
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:02
% EndTime: 2019-03-10 02:16:03
% DurationCPUTime: 0.65s
% Computational Cost: add. (346->100), mult. (669->118), div. (0->0), fcn. (793->12), ass. (0->51)
t160 = -m(4) - m(5);
t159 = -m(6) - m(7);
t158 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t157 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t128 = cos(qJ(4));
t156 = -t128 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t155 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t124 = sin(qJ(4));
t154 = -m(5) * pkin(3) - t128 * mrSges(5,1) + t124 * mrSges(5,2) - mrSges(4,1);
t123 = cos(pkin(6));
t153 = t123 * pkin(8) + pkin(7);
t130 = cos(qJ(2));
t131 = cos(qJ(1));
t142 = t130 * t131;
t126 = sin(qJ(2));
t127 = sin(qJ(1));
t145 = t126 * t127;
t106 = -t123 * t142 + t145;
t152 = t106 * t124;
t143 = t127 * t130;
t144 = t126 * t131;
t108 = t123 * t143 + t144;
t151 = t108 * t124;
t122 = sin(pkin(6));
t150 = t122 * t126;
t149 = t122 * t127;
t129 = cos(qJ(3));
t148 = t122 * t129;
t147 = t122 * t130;
t146 = t122 * t131;
t141 = t131 * pkin(1) + pkin(8) * t149;
t140 = pkin(2) * t150 + t153;
t139 = t127 * pkin(1) - pkin(8) * t146;
t109 = -t123 * t145 + t142;
t138 = t109 * pkin(2) + pkin(9) * t108 + t141;
t137 = -pkin(9) * t147 + t140;
t107 = t123 * t144 + t143;
t136 = t107 * pkin(2) + pkin(9) * t106 + t139;
t132 = -pkin(11) - pkin(10);
t125 = sin(qJ(3));
t121 = qJ(4) + qJ(5);
t117 = cos(t121);
t116 = sin(t121);
t115 = pkin(4) * t128 + pkin(3);
t105 = t123 * t125 + t126 * t148;
t104 = -t123 * t129 + t125 * t150;
t95 = t109 * t129 + t125 * t149;
t94 = t109 * t125 - t127 * t148;
t93 = t107 * t129 - t125 * t146;
t92 = t107 * t125 + t129 * t146;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t153 - t123 * mrSges(3,3) - (t126 * mrSges(3,1) + t130 * mrSges(3,2)) * t122 - m(4) * t137 - t105 * mrSges(4,1) + mrSges(4,3) * t147 - m(5) * (pkin(3) * t105 + t137) - (t105 * t128 - t124 * t147) * mrSges(5,1) - (-t105 * t124 - t128 * t147) * mrSges(5,2) + t159 * (-t104 * t132 + t105 * t115 + (-pkin(4) * t124 - pkin(9)) * t147 + t140) + t158 * (t105 * t117 - t116 * t147) + t157 * (t105 * t116 + t117 * t147) + t155 * t104) * g(3) + (-m(3) * t139 - mrSges(2,1) * t127 - t107 * mrSges(3,1) - t152 * mrSges(5,1) - mrSges(2,2) * t131 + mrSges(3,3) * t146 - mrSges(1,2) + t160 * t136 + t154 * t93 + t159 * (pkin(4) * t152 + t93 * t115 - t92 * t132 + t136) + t158 * (t106 * t116 + t117 * t93) + t157 * (-t106 * t117 + t116 * t93) + t156 * t106 + t155 * t92) * g(2) + (-m(3) * t141 - mrSges(2,1) * t131 - t109 * mrSges(3,1) - t151 * mrSges(5,1) + mrSges(2,2) * t127 - mrSges(3,3) * t149 - mrSges(1,1) + t160 * t138 + t154 * t95 + t159 * (pkin(4) * t151 + t95 * t115 - t94 * t132 + t138) + t158 * (t108 * t116 + t117 * t95) + t157 * (-t108 * t117 + t116 * t95) + t156 * t108 + t155 * t94) * g(1);
U  = t1;
