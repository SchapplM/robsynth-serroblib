% Calculate potential energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:18
% EndTime: 2019-03-09 07:50:19
% DurationCPUTime: 0.87s
% Computational Cost: add. (951->121), mult. (2574->165), div. (0->0), fcn. (3324->18), ass. (0->65)
t169 = sin(pkin(14));
t176 = cos(pkin(6));
t185 = cos(qJ(1));
t173 = cos(pkin(14));
t181 = sin(qJ(1));
t203 = t181 * t173;
t160 = -t169 * t185 - t176 * t203;
t171 = sin(pkin(7));
t175 = cos(pkin(7));
t172 = sin(pkin(6));
t209 = t172 * t181;
t152 = -t160 * t171 + t175 * t209;
t157 = -t171 * t172 * t173 + t175 * t176;
t180 = sin(qJ(3));
t184 = cos(qJ(3));
t207 = t173 * t175;
t210 = t171 * t176;
t149 = t184 * t210 + (-t169 * t180 + t184 * t207) * t172;
t170 = sin(pkin(8));
t174 = cos(pkin(8));
t136 = -t149 * t170 + t157 * t174;
t204 = t181 * t169;
t161 = t173 * t185 - t176 * t204;
t195 = t160 * t175 + t171 * t209;
t141 = -t161 * t180 + t195 * t184;
t132 = -t141 * t170 + t152 * t174;
t205 = t176 * t185;
t159 = t169 * t205 + t203;
t158 = t173 * t205 - t204;
t208 = t172 * t185;
t196 = t158 * t175 - t171 * t208;
t139 = -t159 * t180 + t196 * t184;
t151 = -t158 * t171 - t175 * t208;
t131 = -t139 * t170 + t151 * t174;
t223 = -m(6) - m(7);
t222 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t177 = sin(qJ(6));
t182 = cos(qJ(6));
t221 = -t177 * mrSges(7,1) - t182 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t220 = -m(7) * pkin(5) - t182 * mrSges(7,1) + t177 * mrSges(7,2) - mrSges(6,1);
t219 = cos(qJ(4));
t218 = t176 * qJ(2) + pkin(9);
t202 = t185 * pkin(1) + qJ(2) * t209;
t199 = t170 * t219;
t198 = t174 * t219;
t197 = t181 * pkin(1) - qJ(2) * t208;
t194 = t161 * pkin(2) + t152 * pkin(10) + t202;
t193 = t172 * t169 * pkin(2) + t157 * pkin(10) + t218;
t142 = t161 * t184 + t195 * t180;
t192 = t142 * pkin(3) + t132 * pkin(11) + t194;
t150 = t180 * t210 + (t169 * t184 + t180 * t207) * t172;
t191 = t150 * pkin(3) + t136 * pkin(11) + t193;
t190 = t159 * pkin(2) + t151 * pkin(10) + t197;
t140 = t159 * t184 + t196 * t180;
t187 = t140 * pkin(3) + t131 * pkin(11) + t190;
t183 = cos(qJ(5));
t179 = sin(qJ(4));
t178 = sin(qJ(5));
t128 = t150 * t219 + (t149 * t174 + t157 * t170) * t179;
t127 = -t149 * t198 + t150 * t179 - t157 * t199;
t125 = t142 * t219 + (t141 * t174 + t152 * t170) * t179;
t124 = -t141 * t198 + t142 * t179 - t152 * t199;
t123 = t140 * t219 + (t139 * t174 + t151 * t170) * t179;
t122 = -t139 * t198 + t140 * t179 - t151 * t199;
t1 = (-mrSges(1,3) - m(2) * pkin(9) - mrSges(2,3) - m(3) * t218 - t176 * mrSges(3,3) - (t169 * mrSges(3,1) + t173 * mrSges(3,2)) * t172 - m(4) * t193 - t150 * mrSges(4,1) - t149 * mrSges(4,2) - t157 * mrSges(4,3) - m(5) * t191 - t128 * mrSges(5,1) - t136 * mrSges(5,3) + t223 * (t128 * pkin(4) + pkin(12) * t127 + t191) + t220 * (t128 * t183 + t136 * t178) + t221 * t127 + t222 * (t128 * t178 - t136 * t183)) * g(3) + (-m(3) * t197 - m(4) * t190 - m(5) * t187 - t181 * mrSges(2,1) - t159 * mrSges(3,1) - t140 * mrSges(4,1) - t123 * mrSges(5,1) - t185 * mrSges(2,2) - t158 * mrSges(3,2) - t139 * mrSges(4,2) + mrSges(3,3) * t208 - t151 * mrSges(4,3) - t131 * mrSges(5,3) - mrSges(1,2) + t223 * (t123 * pkin(4) + t122 * pkin(12) + t187) + t220 * (t123 * t183 + t131 * t178) + t221 * t122 + t222 * (t123 * t178 - t131 * t183)) * g(2) + (-m(3) * t202 - m(4) * t194 - m(5) * t192 - t185 * mrSges(2,1) - t161 * mrSges(3,1) - t142 * mrSges(4,1) - t125 * mrSges(5,1) + t181 * mrSges(2,2) - t160 * mrSges(3,2) - t141 * mrSges(4,2) - mrSges(3,3) * t209 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,1) + t223 * (t125 * pkin(4) + pkin(12) * t124 + t192) + t220 * (t125 * t183 + t132 * t178) + t221 * t124 + t222 * (t125 * t178 - t132 * t183)) * g(1);
U  = t1;
