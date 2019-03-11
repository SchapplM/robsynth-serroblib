% Calculate potential energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:36
% EndTime: 2019-03-10 05:49:37
% DurationCPUTime: 0.87s
% Computational Cost: add. (951->121), mult. (2574->164), div. (0->0), fcn. (3324->18), ass. (0->66)
t174 = cos(pkin(6));
t180 = sin(qJ(1));
t184 = cos(qJ(2));
t204 = t180 * t184;
t179 = sin(qJ(2));
t185 = cos(qJ(1));
t206 = t179 * t185;
t160 = -t174 * t204 - t206;
t170 = sin(pkin(7));
t173 = cos(pkin(7));
t171 = sin(pkin(6));
t210 = t171 * t180;
t152 = -t160 * t170 + t173 * t210;
t157 = -t170 * t171 * t184 + t173 * t174;
t178 = sin(qJ(3));
t183 = cos(qJ(3));
t207 = t173 * t184;
t211 = t170 * t174;
t149 = t183 * t211 + (-t178 * t179 + t183 * t207) * t171;
t169 = sin(pkin(8));
t172 = cos(pkin(8));
t136 = -t149 * t169 + t157 * t172;
t203 = t184 * t185;
t205 = t180 * t179;
t161 = -t174 * t205 + t203;
t195 = t160 * t173 + t170 * t210;
t141 = -t161 * t178 + t183 * t195;
t132 = -t141 * t169 + t152 * t172;
t159 = t174 * t206 + t204;
t158 = t174 * t203 - t205;
t209 = t171 * t185;
t196 = t158 * t173 - t170 * t209;
t139 = -t159 * t178 + t183 * t196;
t151 = -t158 * t170 - t173 * t209;
t131 = -t139 * t169 + t151 * t172;
t224 = -m(6) - m(7);
t223 = -m(7) * pkin(14) + mrSges(6,2) - mrSges(7,3);
t175 = sin(qJ(6));
t181 = cos(qJ(6));
t222 = -t175 * mrSges(7,1) - t181 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t221 = -m(7) * pkin(5) - t181 * mrSges(7,1) + t175 * mrSges(7,2) - mrSges(6,1);
t220 = cos(qJ(4));
t219 = t174 * pkin(10) + pkin(9);
t202 = t185 * pkin(1) + pkin(10) * t210;
t199 = t169 * t220;
t198 = t172 * t220;
t197 = t180 * pkin(1) - pkin(10) * t209;
t194 = t161 * pkin(2) + t152 * pkin(11) + t202;
t193 = t171 * t179 * pkin(2) + t157 * pkin(11) + t219;
t142 = t161 * t183 + t178 * t195;
t192 = t142 * pkin(3) + t132 * pkin(12) + t194;
t191 = t159 * pkin(2) + pkin(11) * t151 + t197;
t150 = t178 * t211 + (t178 * t207 + t179 * t183) * t171;
t190 = t150 * pkin(3) + t136 * pkin(12) + t193;
t140 = t159 * t183 + t178 * t196;
t187 = t140 * pkin(3) + t131 * pkin(12) + t191;
t182 = cos(qJ(5));
t177 = sin(qJ(4));
t176 = sin(qJ(5));
t128 = t150 * t220 + (t149 * t172 + t157 * t169) * t177;
t127 = -t149 * t198 + t150 * t177 - t157 * t199;
t125 = t142 * t220 + (t141 * t172 + t152 * t169) * t177;
t124 = -t141 * t198 + t142 * t177 - t152 * t199;
t123 = t140 * t220 + (t139 * t172 + t151 * t169) * t177;
t122 = -t139 * t198 + t140 * t177 - t151 * t199;
t1 = (-mrSges(1,3) - m(2) * pkin(9) - mrSges(2,3) - m(3) * t219 - t174 * mrSges(3,3) - (t179 * mrSges(3,1) + t184 * mrSges(3,2)) * t171 - m(4) * t193 - t150 * mrSges(4,1) - t149 * mrSges(4,2) - t157 * mrSges(4,3) - m(5) * t190 - t128 * mrSges(5,1) - t136 * mrSges(5,3) + t224 * (t128 * pkin(4) + pkin(13) * t127 + t190) + t221 * (t128 * t182 + t136 * t176) + t222 * t127 + t223 * (t128 * t176 - t136 * t182)) * g(3) + (-m(3) * t197 - m(4) * t191 - m(5) * t187 - t180 * mrSges(2,1) - t159 * mrSges(3,1) - t140 * mrSges(4,1) - t123 * mrSges(5,1) - t185 * mrSges(2,2) - t158 * mrSges(3,2) - t139 * mrSges(4,2) + mrSges(3,3) * t209 - t151 * mrSges(4,3) - t131 * mrSges(5,3) - mrSges(1,2) + t224 * (t123 * pkin(4) + t122 * pkin(13) + t187) + t221 * (t123 * t182 + t131 * t176) + t222 * t122 + t223 * (t123 * t176 - t131 * t182)) * g(2) + (-m(3) * t202 - m(4) * t194 - m(5) * t192 - t185 * mrSges(2,1) - t161 * mrSges(3,1) - t142 * mrSges(4,1) - t125 * mrSges(5,1) + t180 * mrSges(2,2) - t160 * mrSges(3,2) - t141 * mrSges(4,2) - mrSges(3,3) * t210 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,1) + t224 * (t125 * pkin(4) + pkin(13) * t124 + t192) + t221 * (t125 * t182 + t132 * t176) + t222 * t124 + t223 * (t125 * t176 - t132 * t182)) * g(1);
U  = t1;
