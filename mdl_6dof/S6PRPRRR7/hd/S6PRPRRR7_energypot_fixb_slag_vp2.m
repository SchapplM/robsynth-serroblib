% Calculate potential energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:29
% EndTime: 2019-03-08 20:50:30
% DurationCPUTime: 0.88s
% Computational Cost: add. (951->121), mult. (2574->166), div. (0->0), fcn. (3324->18), ass. (0->67)
t178 = cos(pkin(7));
t179 = cos(pkin(6));
t186 = cos(qJ(2));
t173 = sin(pkin(7));
t174 = sin(pkin(6));
t213 = t173 * t174;
t158 = t179 * t178 - t186 * t213;
t171 = sin(pkin(13));
t176 = cos(pkin(13));
t183 = sin(qJ(2));
t205 = t179 * t186;
t161 = -t171 * t205 - t176 * t183;
t210 = t174 * t178;
t153 = -t161 * t173 + t171 * t210;
t170 = sin(pkin(14));
t175 = cos(pkin(14));
t208 = t178 * t186;
t212 = t173 * t179;
t150 = t175 * t212 + (-t170 * t183 + t175 * t208) * t174;
t172 = sin(pkin(8));
t177 = cos(pkin(8));
t139 = -t150 * t172 + t158 * t177;
t206 = t179 * t183;
t162 = -t171 * t206 + t176 * t186;
t196 = t161 * t178 + t171 * t213;
t142 = -t162 * t170 + t196 * t175;
t133 = -t142 * t172 + t153 * t177;
t160 = t171 * t186 + t176 * t206;
t159 = -t171 * t183 + t176 * t205;
t211 = t174 * t176;
t197 = t159 * t178 - t173 * t211;
t140 = -t160 * t170 + t197 * t175;
t152 = -t159 * t173 - t176 * t210;
t132 = -t140 * t172 + t152 * t177;
t226 = -m(6) - m(7);
t225 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t180 = sin(qJ(6));
t184 = cos(qJ(6));
t224 = -t180 * mrSges(7,1) - t184 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t223 = -m(7) * pkin(5) - t184 * mrSges(7,1) + t180 * mrSges(7,2) - mrSges(6,1);
t222 = cos(qJ(4));
t214 = t171 * t174;
t209 = t174 * t183;
t204 = t176 * pkin(1) + pkin(9) * t214;
t203 = t179 * pkin(9) + qJ(1);
t200 = t172 * t222;
t199 = t177 * t222;
t198 = t171 * pkin(1) - pkin(9) * t211;
t195 = t162 * pkin(2) + t153 * qJ(3) + t204;
t194 = pkin(2) * t209 + t158 * qJ(3) + t203;
t143 = t162 * t175 + t196 * t170;
t193 = t143 * pkin(3) + t133 * pkin(10) + t195;
t192 = t160 * pkin(2) + t152 * qJ(3) + t198;
t151 = t175 * t209 + (t174 * t208 + t212) * t170;
t191 = t151 * pkin(3) + t139 * pkin(10) + t194;
t141 = t160 * t175 + t197 * t170;
t188 = t141 * pkin(3) + t132 * pkin(10) + t192;
t185 = cos(qJ(5));
t182 = sin(qJ(4));
t181 = sin(qJ(5));
t131 = t151 * t222 + (t150 * t177 + t158 * t172) * t182;
t130 = -t150 * t199 + t151 * t182 - t158 * t200;
t126 = t143 * t222 + (t142 * t177 + t153 * t172) * t182;
t125 = -t142 * t199 + t143 * t182 - t153 * t200;
t124 = t141 * t222 + (t140 * t177 + t152 * t172) * t182;
t123 = -t140 * t199 + t141 * t182 - t152 * t200;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t203 - t179 * mrSges(3,3) - (t183 * mrSges(3,1) + t186 * mrSges(3,2)) * t174 - m(4) * t194 - t151 * mrSges(4,1) - t150 * mrSges(4,2) - t158 * mrSges(4,3) - m(5) * t191 - t131 * mrSges(5,1) - t139 * mrSges(5,3) + t226 * (t131 * pkin(4) + t130 * pkin(11) + t191) + t223 * (t131 * t185 + t139 * t181) + t224 * t130 + t225 * (t131 * t181 - t139 * t185)) * g(3) + (-m(3) * t198 - m(4) * t192 - m(5) * t188 - t171 * mrSges(2,1) - t160 * mrSges(3,1) - t141 * mrSges(4,1) - t124 * mrSges(5,1) - t176 * mrSges(2,2) - t159 * mrSges(3,2) - t140 * mrSges(4,2) + mrSges(3,3) * t211 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,2) + t226 * (t124 * pkin(4) + pkin(11) * t123 + t188) + t223 * (t124 * t185 + t132 * t181) + t224 * t123 + t225 * (t124 * t181 - t132 * t185)) * g(2) + (-m(3) * t204 - m(4) * t195 - m(5) * t193 - t176 * mrSges(2,1) - t162 * mrSges(3,1) - t143 * mrSges(4,1) - t126 * mrSges(5,1) + t171 * mrSges(2,2) - t161 * mrSges(3,2) - t142 * mrSges(4,2) - mrSges(3,3) * t214 - t153 * mrSges(4,3) - t133 * mrSges(5,3) - mrSges(1,1) + t226 * (t126 * pkin(4) + pkin(11) * t125 + t193) + t223 * (t126 * t185 + t133 * t181) + t224 * t125 + t225 * (t126 * t181 - t133 * t185)) * g(1);
U  = t1;
