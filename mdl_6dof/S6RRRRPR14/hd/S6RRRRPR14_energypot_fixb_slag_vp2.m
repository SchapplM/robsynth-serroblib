% Calculate potential energy for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:35
% EndTime: 2019-03-10 00:08:36
% DurationCPUTime: 0.77s
% Computational Cost: add. (563->105), mult. (1397->133), div. (0->0), fcn. (1753->16), ass. (0->55)
t138 = cos(pkin(6));
t143 = sin(qJ(1));
t144 = cos(qJ(2));
t167 = t143 * t144;
t142 = sin(qJ(2));
t145 = cos(qJ(1));
t169 = t142 * t145;
t117 = -t138 * t167 - t169;
t134 = sin(pkin(7));
t137 = cos(pkin(7));
t135 = sin(pkin(6));
t173 = t135 * t143;
t156 = -t117 * t134 + t137 * t173;
t172 = t135 * t144;
t155 = -t134 * t172 + t137 * t138;
t182 = -m(5) - m(6);
t181 = -m(6) * qJ(5) + m(7) * (-pkin(12) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t132 = pkin(13) + qJ(6);
t127 = sin(t132);
t128 = cos(t132);
t133 = sin(pkin(13));
t136 = cos(pkin(13));
t180 = -m(7) * (pkin(5) * t133 + pkin(11)) - t133 * mrSges(6,1) - t127 * mrSges(7,1) - t136 * mrSges(6,2) - t128 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t179 = -m(6) * pkin(4) - m(7) * (pkin(5) * t136 + pkin(4)) - t136 * mrSges(6,1) - t128 * mrSges(7,1) + t133 * mrSges(6,2) + t127 * mrSges(7,2) - mrSges(5,1);
t178 = cos(qJ(3));
t177 = cos(qJ(4));
t176 = t138 * pkin(9) + pkin(8);
t174 = t135 * t142;
t171 = t135 * t145;
t168 = t143 * t142;
t166 = t144 * t145;
t165 = t145 * pkin(1) + pkin(9) * t173;
t161 = t134 * t178;
t160 = t137 * t178;
t159 = t135 * t161;
t158 = t143 * pkin(1) - pkin(9) * t171;
t115 = t138 * t166 - t168;
t157 = -t115 * t134 - t137 * t171;
t118 = -t138 * t168 + t166;
t154 = t118 * pkin(2) + t156 * pkin(10) + t165;
t153 = pkin(2) * t174 + t155 * pkin(10) + t176;
t141 = sin(qJ(3));
t104 = t118 * t178 + (t117 * t137 + t134 * t173) * t141;
t152 = t104 * pkin(3) + t154;
t109 = t138 * t134 * t141 + (t137 * t141 * t144 + t142 * t178) * t135;
t151 = t109 * pkin(3) + t153;
t116 = t138 * t169 + t167;
t148 = t116 * pkin(2) + pkin(10) * t157 + t158;
t102 = t116 * t178 + (t115 * t137 - t134 * t171) * t141;
t147 = t102 * pkin(3) + t148;
t140 = sin(qJ(4));
t108 = -t138 * t161 + t141 * t174 - t160 * t172;
t103 = -t117 * t160 + t118 * t141 - t143 * t159;
t101 = -t115 * t160 + t116 * t141 + t145 * t159;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t176 - t138 * mrSges(3,3) - (t142 * mrSges(3,1) + t144 * mrSges(3,2)) * t135 - m(4) * t153 - t109 * mrSges(4,1) - t155 * mrSges(4,3) - m(7) * t151 + t182 * (pkin(11) * t108 + t151) + t179 * (t109 * t177 + t140 * t155) + t180 * t108 + t181 * (t109 * t140 - t155 * t177)) * g(3) + (-m(3) * t158 - m(4) * t148 - m(7) * t147 - t143 * mrSges(2,1) - t116 * mrSges(3,1) - t102 * mrSges(4,1) - t145 * mrSges(2,2) - t115 * mrSges(3,2) + mrSges(3,3) * t171 - t157 * mrSges(4,3) - mrSges(1,2) + t182 * (t101 * pkin(11) + t147) + t179 * (t102 * t177 + t140 * t157) + t180 * t101 + t181 * (t102 * t140 - t157 * t177)) * g(2) + (-m(3) * t165 - m(4) * t154 - m(7) * t152 - t145 * mrSges(2,1) - t118 * mrSges(3,1) - t104 * mrSges(4,1) + t143 * mrSges(2,2) - t117 * mrSges(3,2) - mrSges(3,3) * t173 - t156 * mrSges(4,3) - mrSges(1,1) + t182 * (pkin(11) * t103 + t152) + t179 * (t104 * t177 + t140 * t156) + t180 * t103 + t181 * (t104 * t140 - t156 * t177)) * g(1);
U  = t1;
