% Calculate potential energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:08
% EndTime: 2019-03-09 04:01:09
% DurationCPUTime: 0.96s
% Computational Cost: add. (598->126), mult. (1473->172), div. (0->0), fcn. (1852->16), ass. (0->64)
t186 = -m(6) - m(7);
t185 = -mrSges(4,3) - mrSges(5,3);
t184 = -m(4) * pkin(9) + t185;
t183 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t150 = sin(qJ(6));
t154 = cos(qJ(6));
t182 = t150 * mrSges(7,1) + t154 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t181 = -m(7) * pkin(5) - t154 * mrSges(7,1) + t150 * mrSges(7,2) - mrSges(6,1);
t180 = pkin(9) + qJ(4);
t149 = cos(pkin(6));
t179 = t149 * qJ(2) + pkin(8);
t145 = sin(pkin(6));
t147 = cos(pkin(12));
t178 = t145 * t147;
t153 = sin(qJ(1));
t177 = t145 * t153;
t157 = cos(qJ(1));
t176 = t145 * t157;
t148 = cos(pkin(7));
t175 = t148 * t149;
t152 = sin(qJ(3));
t174 = t148 * t152;
t173 = t149 * t157;
t143 = sin(pkin(12));
t172 = t153 * t143;
t171 = t153 * t147;
t170 = t157 * pkin(1) + qJ(2) * t177;
t140 = t153 * pkin(1);
t169 = -qJ(2) * t176 + t140;
t144 = sin(pkin(7));
t127 = pkin(3) * t144 * t152 + t148 * t180;
t128 = pkin(3) * t174 - t144 * t180;
t156 = cos(qJ(3));
t138 = pkin(3) * t156 + pkin(2);
t168 = t145 * t143 * t138 + t149 * t127 + t128 * t178 + t179;
t132 = -t143 * t157 - t149 * t171;
t133 = t147 * t157 - t149 * t172;
t167 = t127 * t177 + t132 * t128 + t133 * t138 + t170;
t142 = sin(pkin(13));
t146 = cos(pkin(13));
t166 = t142 * t156 + t146 * t152;
t135 = -t142 * t152 + t146 * t156;
t130 = t147 * t173 - t172;
t115 = -t130 * t144 - t148 * t176;
t165 = t130 * t148 - t144 * t176;
t116 = -t132 * t144 + t148 * t177;
t164 = t132 * t148 + t144 * t177;
t163 = t135 * t144;
t162 = t145 * t163;
t131 = t143 * t173 + t171;
t160 = t130 * t128 + t131 * t138 + t140 + (-qJ(2) - t127) * t176;
t155 = cos(qJ(5));
t151 = sin(qJ(5));
t129 = -t144 * t178 + t175;
t126 = t166 * t148;
t125 = t135 * t148;
t124 = t166 * t144;
t110 = t124 * t149 + (t126 * t147 + t135 * t143) * t145;
t109 = (t125 * t147 - t143 * t166) * t145 + t149 * t163;
t107 = t124 * t177 + t126 * t132 + t133 * t135;
t106 = t132 * t125 - t133 * t166 + t153 * t162;
t105 = -t124 * t176 + t130 * t126 + t131 * t135;
t104 = t125 * t130 - t131 * t166 - t157 * t162;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t179 - m(4) * (pkin(9) * t175 + t179) - m(5) * t168 - t110 * mrSges(5,1) + t186 * (t110 * pkin(4) - pkin(10) * t109 + t168) + (-mrSges(3,3) - (t152 * mrSges(4,1) + t156 * mrSges(4,2)) * t144) * t149 + (-t143 * mrSges(3,1) - t147 * mrSges(3,2) - m(4) * (-pkin(9) * t144 * t147 + pkin(2) * t143) - (t143 * t156 + t147 * t174) * mrSges(4,1) - (t147 * t148 * t156 - t143 * t152) * mrSges(4,2)) * t145 + t185 * t129 + t181 * (t110 * t155 + t129 * t151) + t182 * t109 + t183 * (t110 * t151 - t129 * t155)) * g(3) + (-mrSges(1,2) - t153 * mrSges(2,1) - t157 * mrSges(2,2) - m(3) * t169 - t131 * mrSges(3,1) - t130 * mrSges(3,2) + mrSges(3,3) * t176 - m(4) * (t131 * pkin(2) + t169) - (t131 * t156 + t152 * t165) * mrSges(4,1) - (-t131 * t152 + t156 * t165) * mrSges(4,2) - m(5) * t160 - t105 * mrSges(5,1) + t186 * (t105 * pkin(4) - t104 * pkin(10) + t160) + t183 * (t105 * t151 - t115 * t155) + t184 * t115 + t181 * (t105 * t155 + t115 * t151) + t182 * t104) * g(2) + (-mrSges(1,1) - t157 * mrSges(2,1) + t153 * mrSges(2,2) - m(3) * t170 - t133 * mrSges(3,1) - t132 * mrSges(3,2) - mrSges(3,3) * t177 - m(4) * (pkin(2) * t133 + t170) - (t133 * t156 + t152 * t164) * mrSges(4,1) - (-t133 * t152 + t156 * t164) * mrSges(4,2) - m(5) * t167 - t107 * mrSges(5,1) + t186 * (t107 * pkin(4) - pkin(10) * t106 + t167) + t183 * (t107 * t151 - t116 * t155) + t184 * t116 + t181 * (t107 * t155 + t116 * t151) + t182 * t106) * g(1);
U  = t1;
