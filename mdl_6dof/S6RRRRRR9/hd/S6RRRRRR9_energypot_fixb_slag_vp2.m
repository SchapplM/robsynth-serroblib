% Calculate potential energy for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:37
% EndTime: 2019-03-10 05:16:38
% DurationCPUTime: 0.78s
% Computational Cost: add. (563->105), mult. (1397->133), div. (0->0), fcn. (1753->16), ass. (0->55)
t136 = cos(pkin(6));
t141 = sin(qJ(1));
t143 = cos(qJ(2));
t167 = t141 * t143;
t140 = sin(qJ(2));
t144 = cos(qJ(1));
t168 = t140 * t144;
t117 = -t136 * t167 - t168;
t133 = sin(pkin(7));
t135 = cos(pkin(7));
t134 = sin(pkin(6));
t173 = t134 * t141;
t156 = -t117 * t133 + t135 * t173;
t172 = t134 * t143;
t155 = -t133 * t172 + t135 * t136;
t182 = -m(5) - m(6);
t181 = -m(6) * pkin(12) + m(7) * (-pkin(13) - pkin(12)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t132 = qJ(5) + qJ(6);
t127 = sin(t132);
t128 = cos(t132);
t137 = sin(qJ(5));
t142 = cos(qJ(5));
t180 = -m(7) * (pkin(5) * t137 + pkin(11)) - t137 * mrSges(6,1) - t127 * mrSges(7,1) - t142 * mrSges(6,2) - t128 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t179 = -m(6) * pkin(4) - m(7) * (pkin(5) * t142 + pkin(4)) - mrSges(6,1) * t142 - mrSges(7,1) * t128 + mrSges(6,2) * t137 + mrSges(7,2) * t127 - mrSges(5,1);
t178 = cos(qJ(3));
t177 = cos(qJ(4));
t176 = t136 * pkin(9) + pkin(8);
t174 = t134 * t140;
t171 = t134 * t144;
t169 = t140 * t141;
t166 = t143 * t144;
t165 = t144 * pkin(1) + pkin(9) * t173;
t161 = t133 * t178;
t160 = t135 * t178;
t159 = t134 * t161;
t158 = t141 * pkin(1) - pkin(9) * t171;
t115 = t136 * t166 - t169;
t157 = -t115 * t133 - t135 * t171;
t118 = -t136 * t169 + t166;
t154 = t118 * pkin(2) + t156 * pkin(10) + t165;
t153 = pkin(2) * t174 + t155 * pkin(10) + t176;
t139 = sin(qJ(3));
t104 = t118 * t178 + (t117 * t135 + t133 * t173) * t139;
t152 = t104 * pkin(3) + t154;
t109 = t136 * t133 * t139 + (t135 * t139 * t143 + t140 * t178) * t134;
t151 = t109 * pkin(3) + t153;
t116 = t136 * t168 + t167;
t148 = t116 * pkin(2) + pkin(10) * t157 + t158;
t102 = t116 * t178 + (t115 * t135 - t133 * t171) * t139;
t147 = t102 * pkin(3) + t148;
t138 = sin(qJ(4));
t108 = -t136 * t161 + t139 * t174 - t160 * t172;
t103 = -t117 * t160 + t118 * t139 - t141 * t159;
t101 = -t115 * t160 + t116 * t139 + t144 * t159;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t176 - t136 * mrSges(3,3) - (t140 * mrSges(3,1) + t143 * mrSges(3,2)) * t134 - m(4) * t153 - t109 * mrSges(4,1) - t155 * mrSges(4,3) - m(7) * t151 + t182 * (t108 * pkin(11) + t151) + t179 * (t109 * t177 + t138 * t155) + t180 * t108 + t181 * (t109 * t138 - t155 * t177)) * g(3) + (-m(3) * t158 - m(4) * t148 - m(7) * t147 - t141 * mrSges(2,1) - t116 * mrSges(3,1) - t102 * mrSges(4,1) - t144 * mrSges(2,2) - t115 * mrSges(3,2) + mrSges(3,3) * t171 - t157 * mrSges(4,3) - mrSges(1,2) + t182 * (t101 * pkin(11) + t147) + t179 * (t102 * t177 + t138 * t157) + t180 * t101 + t181 * (t102 * t138 - t157 * t177)) * g(2) + (-m(3) * t165 - m(4) * t154 - m(7) * t152 - t144 * mrSges(2,1) - t118 * mrSges(3,1) - t104 * mrSges(4,1) + t141 * mrSges(2,2) - t117 * mrSges(3,2) - mrSges(3,3) * t173 - t156 * mrSges(4,3) - mrSges(1,1) + t182 * (t103 * pkin(11) + t152) + t179 * (t104 * t177 + t138 * t156) + t180 * t103 + t181 * (t104 * t138 - t156 * t177)) * g(1);
U  = t1;
