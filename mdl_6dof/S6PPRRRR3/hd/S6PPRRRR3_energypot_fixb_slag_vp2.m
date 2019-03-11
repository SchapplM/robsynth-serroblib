% Calculate potential energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:12
% EndTime: 2019-03-08 19:07:13
% DurationCPUTime: 0.90s
% Computational Cost: add. (951->121), mult. (2574->167), div. (0->0), fcn. (3324->18), ass. (0->66)
t174 = cos(pkin(14));
t177 = cos(pkin(7));
t178 = cos(pkin(6));
t172 = sin(pkin(7));
t173 = sin(pkin(6));
t210 = t172 * t173;
t157 = -t174 * t210 + t177 * t178;
t169 = sin(pkin(14));
t175 = cos(pkin(13));
t170 = sin(pkin(13));
t211 = t170 * t178;
t160 = -t169 * t175 - t174 * t211;
t207 = t173 * t177;
t152 = -t160 * t172 + t170 * t207;
t182 = sin(qJ(3));
t185 = cos(qJ(3));
t206 = t174 * t177;
t209 = t172 * t178;
t149 = t185 * t209 + (-t169 * t182 + t185 * t206) * t173;
t171 = sin(pkin(8));
t176 = cos(pkin(8));
t142 = -t149 * t171 + t157 * t176;
t161 = -t169 * t211 + t174 * t175;
t195 = t160 * t177 + t170 * t210;
t140 = -t161 * t182 + t195 * t185;
t132 = -t140 * t171 + t152 * t176;
t205 = t175 * t178;
t159 = t169 * t205 + t170 * t174;
t158 = -t169 * t170 + t174 * t205;
t208 = t173 * t175;
t196 = t158 * t177 - t172 * t208;
t138 = -t159 * t182 + t196 * t185;
t151 = -t158 * t172 - t175 * t207;
t131 = -t138 * t171 + t151 * t176;
t224 = -m(6) - m(7);
t223 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t179 = sin(qJ(6));
t183 = cos(qJ(6));
t222 = -mrSges(7,1) * t179 - mrSges(7,2) * t183 + mrSges(5,2) - mrSges(6,3);
t221 = -m(7) * pkin(5) - t183 * mrSges(7,1) + t179 * mrSges(7,2) - mrSges(6,1);
t220 = cos(qJ(4));
t212 = t170 * t173;
t203 = t178 * qJ(2) + qJ(1);
t202 = t175 * pkin(1) + qJ(2) * t212;
t199 = t171 * t220;
t198 = t176 * t220;
t197 = t170 * pkin(1) - qJ(2) * t208;
t194 = t161 * pkin(2) + t152 * pkin(9) + t202;
t193 = t173 * t169 * pkin(2) + t157 * pkin(9) + t203;
t141 = t161 * t185 + t195 * t182;
t192 = t141 * pkin(3) + t132 * pkin(10) + t194;
t191 = t159 * pkin(2) + t151 * pkin(9) + t197;
t150 = t182 * t209 + (t169 * t185 + t182 * t206) * t173;
t190 = t150 * pkin(3) + t142 * pkin(10) + t193;
t139 = t159 * t185 + t196 * t182;
t187 = t139 * pkin(3) + t131 * pkin(10) + t191;
t184 = cos(qJ(5));
t181 = sin(qJ(4));
t180 = sin(qJ(5));
t130 = t150 * t220 + (t149 * t176 + t157 * t171) * t181;
t129 = -t149 * t198 + t150 * t181 - t157 * t199;
t123 = t141 * t220 + (t140 * t176 + t152 * t171) * t181;
t122 = -t140 * t198 + t141 * t181 - t152 * t199;
t121 = t139 * t220 + (t138 * t176 + t151 * t171) * t181;
t120 = -t138 * t198 + t139 * t181 - t151 * t199;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t203 - t178 * mrSges(3,3) - (mrSges(3,1) * t169 + mrSges(3,2) * t174) * t173 - m(4) * t193 - t150 * mrSges(4,1) - t149 * mrSges(4,2) - t157 * mrSges(4,3) - m(5) * t190 - t130 * mrSges(5,1) - t142 * mrSges(5,3) + t224 * (t130 * pkin(4) + pkin(11) * t129 + t190) + t221 * (t130 * t184 + t142 * t180) + t222 * t129 + t223 * (t130 * t180 - t142 * t184)) * g(3) + (-m(3) * t197 - m(4) * t191 - m(5) * t187 - t170 * mrSges(2,1) - t159 * mrSges(3,1) - t139 * mrSges(4,1) - t121 * mrSges(5,1) - t175 * mrSges(2,2) - t158 * mrSges(3,2) - t138 * mrSges(4,2) + mrSges(3,3) * t208 - t151 * mrSges(4,3) - t131 * mrSges(5,3) - mrSges(1,2) + t224 * (t121 * pkin(4) + pkin(11) * t120 + t187) + t221 * (t121 * t184 + t131 * t180) + t222 * t120 + t223 * (t121 * t180 - t131 * t184)) * g(2) + (-m(3) * t202 - m(4) * t194 - m(5) * t192 - t175 * mrSges(2,1) - t161 * mrSges(3,1) - t141 * mrSges(4,1) - t123 * mrSges(5,1) + t170 * mrSges(2,2) - t160 * mrSges(3,2) - t140 * mrSges(4,2) - mrSges(3,3) * t212 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,1) + t224 * (t123 * pkin(4) + pkin(11) * t122 + t192) + t221 * (t123 * t184 + t132 * t180) + t222 * t122 + t223 * (t123 * t180 - t132 * t184)) * g(1);
U  = t1;
