% Calculate potential energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:41:55
% EndTime: 2019-03-08 18:41:56
% DurationCPUTime: 0.99s
% Computational Cost: add. (598->126), mult. (1473->174), div. (0->0), fcn. (1852->16), ass. (0->64)
t186 = -m(6) - m(7);
t185 = -mrSges(4,3) - mrSges(5,3);
t184 = -m(4) * pkin(8) + t185;
t183 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t152 = sin(qJ(6));
t155 = cos(qJ(6));
t182 = t152 * mrSges(7,1) + t155 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t181 = -m(7) * pkin(5) - t155 * mrSges(7,1) + t152 * mrSges(7,2) - mrSges(6,1);
t180 = pkin(8) + qJ(4);
t144 = sin(pkin(11));
t146 = sin(pkin(6));
t179 = t144 * t146;
t151 = cos(pkin(6));
t178 = t144 * t151;
t148 = cos(pkin(12));
t177 = t146 * t148;
t149 = cos(pkin(11));
t176 = t146 * t149;
t150 = cos(pkin(7));
t175 = t146 * t150;
t174 = t149 * t151;
t173 = t150 * t151;
t154 = sin(qJ(3));
t172 = t150 * t154;
t171 = t151 * qJ(2) + qJ(1);
t170 = t149 * pkin(1) + qJ(2) * t179;
t140 = t144 * pkin(1);
t169 = -qJ(2) * t176 + t140;
t145 = sin(pkin(7));
t127 = pkin(3) * t145 * t154 + t150 * t180;
t128 = pkin(3) * t172 - t145 * t180;
t143 = sin(pkin(12));
t132 = -t143 * t149 - t148 * t178;
t133 = -t143 * t178 + t148 * t149;
t157 = cos(qJ(3));
t138 = pkin(3) * t157 + pkin(2);
t168 = t127 * t179 + t132 * t128 + t133 * t138 + t170;
t167 = t146 * t143 * t138 + t151 * t127 + t128 * t177 + t171;
t142 = sin(pkin(13));
t147 = cos(pkin(13));
t166 = t142 * t157 + t154 * t147;
t135 = -t154 * t142 + t147 * t157;
t130 = -t143 * t144 + t148 * t174;
t115 = -t130 * t145 - t149 * t175;
t165 = t130 * t150 - t145 * t176;
t116 = -t132 * t145 + t144 * t175;
t164 = t132 * t150 + t145 * t179;
t163 = t135 * t145;
t162 = t146 * t163;
t131 = t143 * t174 + t144 * t148;
t161 = t130 * t128 + t131 * t138 + t140 + (-qJ(2) - t127) * t176;
t156 = cos(qJ(5));
t153 = sin(qJ(5));
t129 = -t145 * t177 + t173;
t126 = t166 * t150;
t125 = t135 * t150;
t124 = t166 * t145;
t110 = t124 * t151 + (t126 * t148 + t135 * t143) * t146;
t109 = (t125 * t148 - t143 * t166) * t146 + t151 * t163;
t105 = t124 * t179 + t126 * t132 + t133 * t135;
t104 = t132 * t125 - t133 * t166 + t144 * t162;
t103 = -t124 * t176 + t126 * t130 + t131 * t135;
t102 = t125 * t130 - t131 * t166 - t149 * t162;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t171 - m(4) * (pkin(8) * t173 + t171) - m(5) * t167 - t110 * mrSges(5,1) + t186 * (t110 * pkin(4) - pkin(9) * t109 + t167) + (-mrSges(3,3) - (mrSges(4,1) * t154 + mrSges(4,2) * t157) * t145) * t151 + (-t143 * mrSges(3,1) - t148 * mrSges(3,2) - m(4) * (-pkin(8) * t145 * t148 + pkin(2) * t143) - (t143 * t157 + t148 * t172) * mrSges(4,1) - (t148 * t150 * t157 - t143 * t154) * mrSges(4,2)) * t146 + t185 * t129 + t181 * (t110 * t156 + t129 * t153) + t182 * t109 + t183 * (t110 * t153 - t129 * t156)) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t144 - mrSges(2,2) * t149 - m(3) * t169 - t131 * mrSges(3,1) - t130 * mrSges(3,2) + mrSges(3,3) * t176 - m(4) * (pkin(2) * t131 + t169) - (t131 * t157 + t154 * t165) * mrSges(4,1) - (-t131 * t154 + t157 * t165) * mrSges(4,2) - m(5) * t161 - t103 * mrSges(5,1) + t186 * (t103 * pkin(4) - pkin(9) * t102 + t161) + t183 * (t103 * t153 - t115 * t156) + t184 * t115 + t181 * (t103 * t156 + t115 * t153) + t182 * t102) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t149 + mrSges(2,2) * t144 - m(3) * t170 - t133 * mrSges(3,1) - t132 * mrSges(3,2) - mrSges(3,3) * t179 - m(4) * (pkin(2) * t133 + t170) - (t133 * t157 + t154 * t164) * mrSges(4,1) - (-t133 * t154 + t157 * t164) * mrSges(4,2) - m(5) * t168 - t105 * mrSges(5,1) + t186 * (t105 * pkin(4) - pkin(9) * t104 + t168) + t183 * (t105 * t153 - t116 * t156) + t184 * t116 + t181 * (t105 * t156 + t116 * t153) + t182 * t104) * g(1);
U  = t1;
