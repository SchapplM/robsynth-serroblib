% Calculate potential energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:26
% EndTime: 2019-03-09 18:56:27
% DurationCPUTime: 1.00s
% Computational Cost: add. (598->126), mult. (1473->170), div. (0->0), fcn. (1852->16), ass. (0->66)
t189 = -m(6) - m(7);
t188 = -mrSges(4,3) - mrSges(5,3);
t187 = -m(4) * pkin(10) + t188;
t186 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t149 = sin(qJ(6));
t154 = cos(qJ(6));
t185 = t149 * mrSges(7,1) + t154 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t184 = -m(7) * pkin(5) - t154 * mrSges(7,1) + t149 * mrSges(7,2) - mrSges(6,1);
t151 = sin(qJ(3));
t183 = pkin(3) * t151;
t182 = pkin(10) + qJ(4);
t148 = cos(pkin(6));
t181 = t148 * pkin(9) + pkin(8);
t145 = sin(pkin(6));
t153 = sin(qJ(1));
t180 = t145 * t153;
t157 = cos(qJ(2));
t179 = t145 * t157;
t158 = cos(qJ(1));
t178 = t145 * t158;
t147 = cos(pkin(7));
t177 = t147 * t148;
t176 = t147 * t157;
t152 = sin(qJ(2));
t175 = t152 * t158;
t174 = t153 * t152;
t173 = t153 * t157;
t172 = t157 * t158;
t171 = t158 * pkin(1) + pkin(9) * t180;
t141 = t153 * pkin(1);
t170 = -pkin(9) * t178 + t141;
t144 = sin(pkin(7));
t128 = t144 * t183 + t147 * t182;
t129 = -t144 * t182 + t147 * t183;
t156 = cos(qJ(3));
t139 = pkin(3) * t156 + pkin(2);
t169 = t145 * t152 * t139 + t148 * t128 + t129 * t179 + t181;
t133 = -t148 * t173 - t175;
t134 = -t148 * t174 + t172;
t168 = t128 * t180 + t133 * t129 + t134 * t139 + t171;
t143 = sin(pkin(13));
t146 = cos(pkin(13));
t167 = t143 * t156 + t146 * t151;
t136 = -t143 * t151 + t146 * t156;
t131 = t148 * t172 - t174;
t116 = -t131 * t144 - t147 * t178;
t166 = t131 * t147 - t144 * t178;
t117 = -t133 * t144 + t147 * t180;
t165 = t133 * t147 + t144 * t180;
t164 = t136 * t144;
t163 = t145 * t164;
t132 = t148 * t175 + t173;
t162 = t131 * t129 + t132 * t139 + t141 + (-pkin(9) - t128) * t178;
t155 = cos(qJ(5));
t150 = sin(qJ(5));
t130 = -t144 * t179 + t177;
t127 = t167 * t147;
t126 = t136 * t147;
t125 = t167 * t144;
t111 = t125 * t148 + (t127 * t157 + t136 * t152) * t145;
t110 = (t126 * t157 - t152 * t167) * t145 + t148 * t164;
t108 = t125 * t180 + t127 * t133 + t134 * t136;
t107 = t133 * t126 - t134 * t167 + t153 * t163;
t106 = -t125 * t178 + t131 * t127 + t132 * t136;
t105 = t126 * t131 - t132 * t167 - t158 * t163;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t181 - m(4) * (pkin(10) * t177 + t181) - m(5) * t169 - t111 * mrSges(5,1) + t189 * (t111 * pkin(4) - pkin(11) * t110 + t169) + (-mrSges(3,3) - (mrSges(4,1) * t151 + mrSges(4,2) * t156) * t144) * t148 + (-t152 * mrSges(3,1) - t157 * mrSges(3,2) - m(4) * (-pkin(10) * t144 * t157 + pkin(2) * t152) - (t151 * t176 + t152 * t156) * mrSges(4,1) - (-t151 * t152 + t156 * t176) * mrSges(4,2)) * t145 + t188 * t130 + t184 * (t111 * t155 + t130 * t150) + t185 * t110 + t186 * (t111 * t150 - t130 * t155)) * g(3) + (-mrSges(1,2) - t153 * mrSges(2,1) - t158 * mrSges(2,2) - m(3) * t170 - t132 * mrSges(3,1) - t131 * mrSges(3,2) + mrSges(3,3) * t178 - m(4) * (t132 * pkin(2) + t170) - (t132 * t156 + t151 * t166) * mrSges(4,1) - (-t132 * t151 + t156 * t166) * mrSges(4,2) - m(5) * t162 - t106 * mrSges(5,1) + t189 * (t106 * pkin(4) - t105 * pkin(11) + t162) + t186 * (t106 * t150 - t116 * t155) + t187 * t116 + t184 * (t106 * t155 + t116 * t150) + t185 * t105) * g(2) + (-mrSges(1,1) - t158 * mrSges(2,1) + t153 * mrSges(2,2) - m(3) * t171 - t134 * mrSges(3,1) - t133 * mrSges(3,2) - mrSges(3,3) * t180 - m(4) * (pkin(2) * t134 + t171) - (t134 * t156 + t151 * t165) * mrSges(4,1) - (-t134 * t151 + t156 * t165) * mrSges(4,2) - m(5) * t168 - t108 * mrSges(5,1) + t189 * (t108 * pkin(4) - pkin(11) * t107 + t168) + t186 * (t108 * t150 - t117 * t155) + t187 * t117 + t184 * (t108 * t155 + t117 * t150) + t185 * t107) * g(1);
U  = t1;
