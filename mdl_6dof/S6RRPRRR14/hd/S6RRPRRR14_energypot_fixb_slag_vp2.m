% Calculate potential energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:07:43
% EndTime: 2019-01-03 10:07:44
% DurationCPUTime: 0.89s
% Computational Cost: add. (951->121), mult. (2574->163), div. (0->0), fcn. (3324->18), ass. (0->67)
t177 = cos(pkin(6));
t182 = sin(qJ(1));
t185 = cos(qJ(2));
t205 = t182 * t185;
t181 = sin(qJ(2));
t186 = cos(qJ(1));
t207 = t181 * t186;
t161 = -t177 * t205 - t207;
t172 = sin(pkin(7));
t176 = cos(pkin(7));
t173 = sin(pkin(6));
t211 = t173 * t182;
t153 = -t161 * t172 + t176 * t211;
t158 = -t172 * t173 * t185 + t176 * t177;
t170 = sin(pkin(14));
t174 = cos(pkin(14));
t208 = t176 * t185;
t213 = t172 * t177;
t150 = t174 * t213 + (-t170 * t181 + t174 * t208) * t173;
t171 = sin(pkin(8));
t175 = cos(pkin(8));
t137 = -t150 * t171 + t158 * t175;
t204 = t185 * t186;
t206 = t182 * t181;
t162 = -t177 * t206 + t204;
t196 = t161 * t176 + t172 * t211;
t142 = -t162 * t170 + t196 * t174;
t133 = -t142 * t171 + t153 * t175;
t160 = t177 * t207 + t205;
t159 = t177 * t204 - t206;
t210 = t173 * t186;
t197 = t159 * t176 - t172 * t210;
t140 = -t160 * t170 + t197 * t174;
t152 = -t159 * t172 - t176 * t210;
t132 = -t140 * t171 + t152 * t175;
t226 = -m(6) - m(7);
t225 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t178 = sin(qJ(6));
t183 = cos(qJ(6));
t224 = -t178 * mrSges(7,1) - t183 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t223 = -m(7) * pkin(5) - t183 * mrSges(7,1) + t178 * mrSges(7,2) - mrSges(6,1);
t222 = cos(qJ(4));
t221 = t177 * pkin(10) + pkin(9);
t212 = t173 * t181;
t203 = t186 * pkin(1) + pkin(10) * t211;
t200 = t171 * t222;
t199 = t175 * t222;
t198 = t182 * pkin(1) - pkin(10) * t210;
t195 = t162 * pkin(2) + t153 * qJ(3) + t203;
t194 = pkin(2) * t212 + t158 * qJ(3) + t221;
t143 = t162 * t174 + t196 * t170;
t193 = t143 * pkin(3) + t133 * pkin(11) + t195;
t192 = t160 * pkin(2) + t152 * qJ(3) + t198;
t151 = t174 * t212 + (t173 * t208 + t213) * t170;
t191 = t151 * pkin(3) + t137 * pkin(11) + t194;
t141 = t160 * t174 + t197 * t170;
t188 = t141 * pkin(3) + t132 * pkin(11) + t192;
t184 = cos(qJ(5));
t180 = sin(qJ(4));
t179 = sin(qJ(5));
t129 = t151 * t222 + (t150 * t175 + t158 * t171) * t180;
t128 = -t150 * t199 + t151 * t180 - t158 * t200;
t126 = t143 * t222 + (t142 * t175 + t153 * t171) * t180;
t125 = -t142 * t199 + t143 * t180 - t153 * t200;
t124 = t141 * t222 + (t140 * t175 + t152 * t171) * t180;
t123 = -t140 * t199 + t141 * t180 - t152 * t200;
t1 = (-mrSges(1,3) - m(2) * pkin(9) - mrSges(2,3) - m(3) * t221 - t177 * mrSges(3,3) - (t181 * mrSges(3,1) + t185 * mrSges(3,2)) * t173 - m(4) * t194 - t151 * mrSges(4,1) - t150 * mrSges(4,2) - t158 * mrSges(4,3) - m(5) * t191 - t129 * mrSges(5,1) - t137 * mrSges(5,3) + t226 * (t129 * pkin(4) + pkin(12) * t128 + t191) + t223 * (t129 * t184 + t137 * t179) + t224 * t128 + t225 * (t129 * t179 - t137 * t184)) * g(3) + (-m(3) * t198 - m(4) * t192 - m(5) * t188 - t182 * mrSges(2,1) - t160 * mrSges(3,1) - t141 * mrSges(4,1) - t124 * mrSges(5,1) - t186 * mrSges(2,2) - t159 * mrSges(3,2) - t140 * mrSges(4,2) + mrSges(3,3) * t210 - t152 * mrSges(4,3) - t132 * mrSges(5,3) - mrSges(1,2) + t226 * (t124 * pkin(4) + t123 * pkin(12) + t188) + t223 * (t124 * t184 + t132 * t179) + t224 * t123 + t225 * (t124 * t179 - t132 * t184)) * g(2) + (-m(3) * t203 - m(4) * t195 - m(5) * t193 - t186 * mrSges(2,1) - t162 * mrSges(3,1) - t143 * mrSges(4,1) - t126 * mrSges(5,1) + t182 * mrSges(2,2) - t161 * mrSges(3,2) - t142 * mrSges(4,2) - mrSges(3,3) * t211 - t153 * mrSges(4,3) - t133 * mrSges(5,3) - mrSges(1,1) + t226 * (t126 * pkin(4) + pkin(12) * t125 + t193) + t223 * (t126 * t184 + t133 * t179) + t224 * t125 + t225 * (t126 * t179 - t133 * t184)) * g(1);
U  = t1;
