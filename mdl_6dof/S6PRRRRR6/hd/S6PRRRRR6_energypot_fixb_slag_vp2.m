% Calculate potential energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:13:58
% EndTime: 2019-03-09 01:13:59
% DurationCPUTime: 0.88s
% Computational Cost: add. (951->121), mult. (2574->167), div. (0->0), fcn. (3324->18), ass. (0->66)
t175 = cos(pkin(7));
t176 = cos(pkin(6));
t185 = cos(qJ(2));
t171 = sin(pkin(7));
t172 = sin(pkin(6));
t211 = t171 * t172;
t157 = t176 * t175 - t185 * t211;
t169 = sin(pkin(14));
t173 = cos(pkin(14));
t181 = sin(qJ(2));
t204 = t176 * t185;
t160 = -t169 * t204 - t173 * t181;
t208 = t172 * t175;
t152 = -t160 * t171 + t169 * t208;
t180 = sin(qJ(3));
t184 = cos(qJ(3));
t207 = t175 * t185;
t210 = t171 * t176;
t149 = t184 * t210 + (-t180 * t181 + t184 * t207) * t172;
t170 = sin(pkin(8));
t174 = cos(pkin(8));
t138 = -t149 * t170 + t157 * t174;
t205 = t176 * t181;
t161 = -t169 * t205 + t173 * t185;
t195 = t160 * t175 + t169 * t211;
t141 = -t161 * t180 + t184 * t195;
t132 = -t141 * t170 + t152 * t174;
t159 = t169 * t185 + t173 * t205;
t158 = -t169 * t181 + t173 * t204;
t209 = t172 * t173;
t196 = t158 * t175 - t171 * t209;
t139 = -t159 * t180 + t184 * t196;
t151 = -t158 * t171 - t173 * t208;
t131 = -t139 * t170 + t151 * t174;
t224 = -m(6) - m(7);
t223 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t177 = sin(qJ(6));
t182 = cos(qJ(6));
t222 = -t177 * mrSges(7,1) - t182 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t221 = -m(7) * pkin(5) - t182 * mrSges(7,1) + t177 * mrSges(7,2) - mrSges(6,1);
t220 = cos(qJ(4));
t212 = t169 * t172;
t203 = t173 * pkin(1) + pkin(9) * t212;
t202 = t176 * pkin(9) + qJ(1);
t199 = t170 * t220;
t198 = t174 * t220;
t197 = t169 * pkin(1) - pkin(9) * t209;
t194 = t161 * pkin(2) + t152 * pkin(10) + t203;
t193 = t172 * t181 * pkin(2) + t157 * pkin(10) + t202;
t142 = t161 * t184 + t180 * t195;
t192 = t142 * pkin(3) + t132 * pkin(11) + t194;
t191 = t159 * pkin(2) + pkin(10) * t151 + t197;
t150 = t180 * t210 + (t180 * t207 + t181 * t184) * t172;
t190 = t150 * pkin(3) + t138 * pkin(11) + t193;
t140 = t159 * t184 + t180 * t196;
t187 = t140 * pkin(3) + t131 * pkin(11) + t191;
t183 = cos(qJ(5));
t179 = sin(qJ(4));
t178 = sin(qJ(5));
t130 = t150 * t220 + (t149 * t174 + t157 * t170) * t179;
t129 = -t149 * t198 + t150 * t179 - t157 * t199;
t125 = t142 * t220 + (t141 * t174 + t152 * t170) * t179;
t124 = -t141 * t198 + t142 * t179 - t152 * t199;
t123 = t140 * t220 + (t139 * t174 + t151 * t170) * t179;
t122 = -t139 * t198 + t140 * t179 - t151 * t199;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t202 - t176 * mrSges(3,3) - (t181 * mrSges(3,1) + t185 * mrSges(3,2)) * t172 - m(4) * t193 - t150 * mrSges(4,1) - t149 * mrSges(4,2) - t157 * mrSges(4,3) - m(5) * t190 - t130 * mrSges(5,1) - t138 * mrSges(5,3) + t224 * (t130 * pkin(4) + t129 * pkin(12) + t190) + t221 * (t130 * t183 + t138 * t178) + t222 * t129 + t223 * (t130 * t178 - t138 * t183)) * g(3) + (-m(3) * t197 - m(4) * t191 - m(5) * t187 - t169 * mrSges(2,1) - t159 * mrSges(3,1) - t140 * mrSges(4,1) - t123 * mrSges(5,1) - t173 * mrSges(2,2) - t158 * mrSges(3,2) - t139 * mrSges(4,2) + mrSges(3,3) * t209 - t151 * mrSges(4,3) - t131 * mrSges(5,3) - mrSges(1,2) + t224 * (t123 * pkin(4) + pkin(12) * t122 + t187) + t221 * (t123 * t183 + t131 * t178) + t222 * t122 + t223 * (t123 * t178 - t131 * t183)) * g(2) + (-m(3) * t203 - m(4) * t194 - m(5) * t192 - t173 * mrSges(2,1) - t161 * mrSges(3,1) - t142 * mrSges(4,1) - t125 * mrSges(5,1) + t169 * mrSges(2,2) - t160 * mrSges(3,2) - t141 * mrSges(4,2) - mrSges(3,3) * t212 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,1) + t224 * (t125 * pkin(4) + pkin(12) * t124 + t192) + t221 * (t125 * t183 + t132 * t178) + t222 * t124 + t223 * (t125 * t178 - t132 * t183)) * g(1);
U  = t1;
