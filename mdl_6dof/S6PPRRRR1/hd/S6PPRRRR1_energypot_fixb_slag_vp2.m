% Calculate potential energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:38
% EndTime: 2019-03-08 18:59:39
% DurationCPUTime: 0.90s
% Computational Cost: add. (526->116), mult. (1200->153), div. (0->0), fcn. (1486->16), ass. (0->56)
t149 = sin(pkin(7));
t153 = cos(pkin(7));
t154 = cos(pkin(6));
t150 = sin(pkin(6));
t151 = cos(pkin(13));
t179 = t150 * t151;
t128 = -t149 * t179 + t153 * t154;
t147 = sin(pkin(13));
t152 = cos(pkin(12));
t148 = sin(pkin(12));
t180 = t148 * t154;
t131 = -t147 * t152 - t151 * t180;
t177 = t150 * t153;
t121 = -t131 * t149 + t148 * t177;
t191 = -m(6) - m(7);
t190 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t155 = sin(qJ(6));
t158 = cos(qJ(6));
t189 = -m(7) * pkin(5) - t158 * mrSges(7,1) + t155 * mrSges(7,2) - mrSges(6,1);
t188 = -m(5) * pkin(9) - t155 * mrSges(7,1) - t158 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t187 = cos(qJ(3));
t176 = t152 * t154;
t129 = -t147 * t148 + t151 * t176;
t120 = -t129 * t149 - t152 * t177;
t156 = sin(qJ(4));
t186 = t120 * t156;
t185 = t121 * t156;
t184 = t128 * t156;
t182 = t147 * t150;
t181 = t148 * t150;
t178 = t150 * t152;
t174 = t154 * qJ(2) + qJ(1);
t173 = t152 * pkin(1) + qJ(2) * t181;
t170 = t149 * t187;
t169 = t153 * t187;
t168 = t150 * t170;
t167 = t148 * pkin(1) - qJ(2) * t178;
t132 = -t147 * t180 + t151 * t152;
t166 = t132 * pkin(2) + t121 * pkin(8) + t173;
t165 = pkin(2) * t182 + t128 * pkin(8) + t174;
t130 = t147 * t176 + t148 * t151;
t162 = t130 * pkin(2) + pkin(8) * t120 + t167;
t160 = -pkin(10) - pkin(9);
t159 = cos(qJ(4));
t157 = sin(qJ(3));
t146 = qJ(4) + qJ(5);
t143 = cos(t146);
t142 = sin(t146);
t140 = pkin(4) * t159 + pkin(3);
t119 = t154 * t149 * t157 + (t151 * t153 * t157 + t147 * t187) * t150;
t118 = -t154 * t170 + t157 * t182 - t169 * t179;
t111 = t132 * t187 + (t131 * t153 + t149 * t181) * t157;
t110 = -t131 * t169 + t132 * t157 - t148 * t168;
t109 = t130 * t187 + (t129 * t153 - t149 * t178) * t157;
t108 = -t129 * t169 + t130 * t157 + t152 * t168;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t174 - t154 * mrSges(3,3) - (t147 * mrSges(3,1) + t151 * mrSges(3,2)) * t150 - m(4) * t165 - t119 * mrSges(4,1) - t128 * mrSges(4,3) - m(5) * (pkin(3) * t119 + t165) - (t119 * t159 + t184) * mrSges(5,1) - (-t119 * t156 + t128 * t159) * mrSges(5,2) + t191 * (pkin(4) * t184 - t118 * t160 + t119 * t140 + t165) + t190 * (t119 * t142 - t128 * t143) + t189 * (t119 * t143 + t128 * t142) + t188 * t118) * g(3) + (-mrSges(1,2) - t148 * mrSges(2,1) - t152 * mrSges(2,2) - m(3) * t167 - t130 * mrSges(3,1) - t129 * mrSges(3,2) + mrSges(3,3) * t178 - m(4) * t162 - t109 * mrSges(4,1) - t120 * mrSges(4,3) - m(5) * (pkin(3) * t109 + t162) - (t109 * t159 + t186) * mrSges(5,1) - (-t109 * t156 + t120 * t159) * mrSges(5,2) + t191 * (pkin(4) * t186 - t108 * t160 + t109 * t140 + t162) + t190 * (t109 * t142 - t120 * t143) + t189 * (t109 * t143 + t120 * t142) + t188 * t108) * g(2) + (-mrSges(1,1) - t152 * mrSges(2,1) + t148 * mrSges(2,2) - m(3) * t173 - t132 * mrSges(3,1) - t131 * mrSges(3,2) - mrSges(3,3) * t181 - m(4) * t166 - t111 * mrSges(4,1) - t121 * mrSges(4,3) - m(5) * (pkin(3) * t111 + t166) - (t111 * t159 + t185) * mrSges(5,1) - (-t111 * t156 + t121 * t159) * mrSges(5,2) + t191 * (pkin(4) * t185 - t110 * t160 + t111 * t140 + t166) + t190 * (t111 * t142 - t121 * t143) + t189 * (t111 * t143 + t121 * t142) + t188 * t110) * g(1);
U  = t1;
