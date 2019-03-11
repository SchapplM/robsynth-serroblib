% Calculate potential energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPPRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:03
% EndTime: 2019-03-08 18:39:04
% DurationCPUTime: 0.87s
% Computational Cost: add. (951->121), mult. (2574->166), div. (0->0), fcn. (3324->18), ass. (0->67)
t177 = cos(pkin(13));
t180 = cos(pkin(7));
t181 = cos(pkin(6));
t174 = sin(pkin(7));
t175 = sin(pkin(6));
t211 = t174 * t175;
t158 = -t177 * t211 + t180 * t181;
t171 = sin(pkin(13));
t178 = cos(pkin(12));
t172 = sin(pkin(12));
t212 = t172 * t181;
t161 = -t171 * t178 - t177 * t212;
t208 = t175 * t180;
t153 = -t161 * t174 + t172 * t208;
t170 = sin(pkin(14));
t176 = cos(pkin(14));
t207 = t177 * t180;
t210 = t174 * t181;
t150 = t176 * t210 + (-t170 * t171 + t176 * t207) * t175;
t173 = sin(pkin(8));
t179 = cos(pkin(8));
t143 = -t150 * t173 + t158 * t179;
t162 = -t171 * t212 + t177 * t178;
t196 = t161 * t180 + t172 * t211;
t141 = -t162 * t170 + t196 * t176;
t133 = -t141 * t173 + t153 * t179;
t206 = t178 * t181;
t160 = t171 * t206 + t172 * t177;
t159 = -t171 * t172 + t177 * t206;
t209 = t175 * t178;
t197 = t159 * t180 - t174 * t209;
t139 = -t160 * t170 + t197 * t176;
t152 = -t159 * t174 - t178 * t208;
t132 = -t139 * t173 + t152 * t179;
t226 = -m(6) - m(7);
t225 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t182 = sin(qJ(6));
t185 = cos(qJ(6));
t224 = -t182 * mrSges(7,1) - t185 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t223 = -m(7) * pkin(5) - t185 * mrSges(7,1) + t182 * mrSges(7,2) - mrSges(6,1);
t222 = cos(qJ(4));
t214 = t171 * t175;
t213 = t172 * t175;
t204 = t181 * qJ(2) + qJ(1);
t203 = t178 * pkin(1) + qJ(2) * t213;
t200 = t173 * t222;
t199 = t179 * t222;
t198 = t172 * pkin(1) - qJ(2) * t209;
t195 = t162 * pkin(2) + t153 * qJ(3) + t203;
t194 = pkin(2) * t214 + t158 * qJ(3) + t204;
t142 = t162 * t176 + t196 * t170;
t193 = t142 * pkin(3) + t133 * pkin(9) + t195;
t192 = t160 * pkin(2) + t152 * qJ(3) + t198;
t151 = t176 * t214 + (t175 * t207 + t210) * t170;
t191 = t151 * pkin(3) + t143 * pkin(9) + t194;
t140 = t160 * t176 + t197 * t170;
t188 = t140 * pkin(3) + t132 * pkin(9) + t192;
t186 = cos(qJ(5));
t184 = sin(qJ(4));
t183 = sin(qJ(5));
t131 = t151 * t222 + (t150 * t179 + t158 * t173) * t184;
t130 = -t150 * t199 + t151 * t184 - t158 * t200;
t126 = t142 * t222 + (t141 * t179 + t153 * t173) * t184;
t125 = -t141 * t199 + t142 * t184 - t153 * t200;
t124 = t140 * t222 + (t139 * t179 + t152 * t173) * t184;
t123 = -t139 * t199 + t140 * t184 - t152 * t200;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t204 - t181 * mrSges(3,3) - (t171 * mrSges(3,1) + t177 * mrSges(3,2)) * t175 - m(4) * t194 - t151 * mrSges(4,1) - t150 * mrSges(4,2) - t158 * mrSges(4,3) - m(5) * t191 - t131 * mrSges(5,1) - t143 * mrSges(5,3) + t226 * (t131 * pkin(4) + pkin(10) * t130 + t191) + t223 * (t131 * t186 + t143 * t183) + t224 * t130 + t225 * (t131 * t183 - t143 * t186)) * g(3) + (-m(3) * t198 - m(4) * t192 - m(5) * t188 - t172 * mrSges(2,1) - t160 * mrSges(3,1) - t140 * mrSges(4,1) - t124 * mrSges(5,1) - t178 * mrSges(2,2) - t159 * mrSges(3,2) - t139 * mrSges(4,2) + mrSges(3,3) * t209 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,2) + t226 * (t124 * pkin(4) + pkin(10) * t123 + t188) + t223 * (t124 * t186 + t132 * t183) + t224 * t123 + t225 * (t124 * t183 - t132 * t186)) * g(2) + (-m(3) * t203 - m(4) * t195 - m(5) * t193 - t178 * mrSges(2,1) - t162 * mrSges(3,1) - t142 * mrSges(4,1) - t126 * mrSges(5,1) + t172 * mrSges(2,2) - t161 * mrSges(3,2) - t141 * mrSges(4,2) - mrSges(3,3) * t213 - t153 * mrSges(4,3) - t133 * mrSges(5,3) - mrSges(1,1) + t226 * (t126 * pkin(4) + pkin(10) * t125 + t193) + t223 * (t126 * t186 + t133 * t183) + t224 * t125 + t225 * (t126 * t183 - t133 * t186)) * g(1);
U  = t1;
