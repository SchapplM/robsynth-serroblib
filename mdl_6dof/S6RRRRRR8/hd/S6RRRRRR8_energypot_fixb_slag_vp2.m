% Calculate potential energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:12
% EndTime: 2019-03-10 04:48:13
% DurationCPUTime: 0.88s
% Computational Cost: add. (526->116), mult. (1200->150), div. (0->0), fcn. (1486->16), ass. (0->57)
t150 = cos(pkin(6));
t155 = sin(qJ(1));
t158 = cos(qJ(2));
t175 = t155 * t158;
t154 = sin(qJ(2));
t159 = cos(qJ(1));
t176 = t154 * t159;
t131 = -t150 * t175 - t176;
t147 = sin(pkin(7));
t149 = cos(pkin(7));
t148 = sin(pkin(6));
t181 = t148 * t155;
t121 = -t131 * t147 + t149 * t181;
t180 = t148 * t158;
t128 = -t147 * t180 + t149 * t150;
t192 = -m(6) - m(7);
t191 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t151 = sin(qJ(6));
t156 = cos(qJ(6));
t190 = -m(7) * pkin(5) - t156 * mrSges(7,1) + t151 * mrSges(7,2) - mrSges(6,1);
t189 = -m(5) * pkin(11) - t151 * mrSges(7,1) - t156 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t188 = cos(qJ(3));
t187 = t150 * pkin(9) + pkin(8);
t174 = t158 * t159;
t177 = t154 * t155;
t129 = t150 * t174 - t177;
t179 = t148 * t159;
t120 = -t129 * t147 - t149 * t179;
t152 = sin(qJ(4));
t186 = t120 * t152;
t185 = t121 * t152;
t184 = t128 * t152;
t182 = t148 * t154;
t173 = t159 * pkin(1) + pkin(9) * t181;
t170 = t147 * t188;
t169 = t149 * t188;
t168 = t148 * t170;
t167 = t155 * pkin(1) - pkin(9) * t179;
t132 = -t150 * t177 + t174;
t166 = t132 * pkin(2) + t121 * pkin(10) + t173;
t165 = pkin(2) * t182 + t128 * pkin(10) + t187;
t130 = t150 * t176 + t175;
t162 = t130 * pkin(2) + pkin(10) * t120 + t167;
t160 = -pkin(12) - pkin(11);
t157 = cos(qJ(4));
t153 = sin(qJ(3));
t146 = qJ(4) + qJ(5);
t142 = cos(t146);
t141 = sin(t146);
t140 = pkin(4) * t157 + pkin(3);
t119 = t150 * t147 * t153 + (t149 * t153 * t158 + t154 * t188) * t148;
t118 = -t150 * t170 + t153 * t182 - t169 * t180;
t111 = t132 * t188 + (t131 * t149 + t147 * t181) * t153;
t110 = -t131 * t169 + t132 * t153 - t155 * t168;
t109 = t130 * t188 + (t129 * t149 - t147 * t179) * t153;
t108 = -t129 * t169 + t130 * t153 + t159 * t168;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t187 - t150 * mrSges(3,3) - (t154 * mrSges(3,1) + t158 * mrSges(3,2)) * t148 - m(4) * t165 - t119 * mrSges(4,1) - t128 * mrSges(4,3) - m(5) * (pkin(3) * t119 + t165) - (t119 * t157 + t184) * mrSges(5,1) - (-t119 * t152 + t128 * t157) * mrSges(5,2) + t192 * (pkin(4) * t184 - t118 * t160 + t119 * t140 + t165) + t191 * (t119 * t141 - t128 * t142) + t190 * (t119 * t142 + t128 * t141) + t189 * t118) * g(3) + (-mrSges(1,2) - t155 * mrSges(2,1) - t159 * mrSges(2,2) - m(3) * t167 - t130 * mrSges(3,1) - t129 * mrSges(3,2) + mrSges(3,3) * t179 - m(4) * t162 - t109 * mrSges(4,1) - t120 * mrSges(4,3) - m(5) * (pkin(3) * t109 + t162) - (t109 * t157 + t186) * mrSges(5,1) - (-t109 * t152 + t120 * t157) * mrSges(5,2) + t192 * (pkin(4) * t186 - t108 * t160 + t109 * t140 + t162) + t191 * (t109 * t141 - t120 * t142) + t190 * (t109 * t142 + t120 * t141) + t189 * t108) * g(2) + (-mrSges(1,1) - t159 * mrSges(2,1) + t155 * mrSges(2,2) - m(3) * t173 - t132 * mrSges(3,1) - t131 * mrSges(3,2) - mrSges(3,3) * t181 - m(4) * t166 - t111 * mrSges(4,1) - t121 * mrSges(4,3) - m(5) * (pkin(3) * t111 + t166) - (t111 * t157 + t185) * mrSges(5,1) - (-t111 * t152 + t121 * t157) * mrSges(5,2) + t192 * (pkin(4) * t185 - t110 * t160 + t111 * t140 + t166) + t191 * (t111 * t141 - t121 * t142) + t190 * (t111 * t142 + t121 * t141) + t189 * t110) * g(1);
U  = t1;
