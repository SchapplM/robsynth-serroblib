% Calculate potential energy for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:44
% EndTime: 2019-03-09 00:54:45
% DurationCPUTime: 0.87s
% Computational Cost: add. (526->116), mult. (1200->153), div. (0->0), fcn. (1486->16), ass. (0->56)
t147 = sin(pkin(7));
t150 = cos(pkin(7));
t151 = cos(pkin(6));
t148 = sin(pkin(6));
t158 = cos(qJ(2));
t177 = t148 * t158;
t127 = -t147 * t177 + t150 * t151;
t146 = sin(pkin(13));
t149 = cos(pkin(13));
t155 = sin(qJ(2));
t174 = t151 * t158;
t130 = -t146 * t174 - t149 * t155;
t179 = t148 * t150;
t120 = -t130 * t147 + t146 * t179;
t190 = -m(6) - m(7);
t189 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t152 = sin(qJ(6));
t156 = cos(qJ(6));
t188 = -m(7) * pkin(5) - t156 * mrSges(7,1) + t152 * mrSges(7,2) - mrSges(6,1);
t187 = -m(5) * pkin(10) - t152 * mrSges(7,1) - t156 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t186 = cos(qJ(3));
t128 = -t146 * t155 + t149 * t174;
t119 = -t128 * t147 - t149 * t179;
t153 = sin(qJ(4));
t185 = t119 * t153;
t184 = t120 * t153;
t183 = t127 * t153;
t181 = t146 * t148;
t180 = t148 * t149;
t178 = t148 * t155;
t175 = t151 * t155;
t173 = t149 * pkin(1) + pkin(8) * t181;
t172 = t151 * pkin(8) + qJ(1);
t169 = t147 * t186;
t168 = t150 * t186;
t167 = t148 * t169;
t166 = t146 * pkin(1) - pkin(8) * t180;
t131 = -t146 * t175 + t149 * t158;
t165 = t131 * pkin(2) + t120 * pkin(9) + t173;
t164 = pkin(2) * t178 + t127 * pkin(9) + t172;
t129 = t146 * t158 + t149 * t175;
t161 = t129 * pkin(2) + pkin(9) * t119 + t166;
t159 = -pkin(11) - pkin(10);
t157 = cos(qJ(4));
t154 = sin(qJ(3));
t145 = qJ(4) + qJ(5);
t141 = cos(t145);
t140 = sin(t145);
t139 = pkin(4) * t157 + pkin(3);
t118 = t151 * t147 * t154 + (t150 * t154 * t158 + t155 * t186) * t148;
t117 = -t151 * t169 + t154 * t178 - t168 * t177;
t110 = t131 * t186 + (t130 * t150 + t147 * t181) * t154;
t109 = -t130 * t168 + t131 * t154 - t146 * t167;
t108 = t129 * t186 + (t128 * t150 - t147 * t180) * t154;
t107 = -t128 * t168 + t129 * t154 + t149 * t167;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t172 - t151 * mrSges(3,3) - (t155 * mrSges(3,1) + t158 * mrSges(3,2)) * t148 - m(4) * t164 - t118 * mrSges(4,1) - t127 * mrSges(4,3) - m(5) * (pkin(3) * t118 + t164) - (t118 * t157 + t183) * mrSges(5,1) - (-t118 * t153 + t127 * t157) * mrSges(5,2) + t190 * (pkin(4) * t183 - t117 * t159 + t118 * t139 + t164) + t189 * (t118 * t140 - t127 * t141) + t188 * (t118 * t141 + t127 * t140) + t187 * t117) * g(3) + (-mrSges(1,2) - t146 * mrSges(2,1) - t149 * mrSges(2,2) - m(3) * t166 - t129 * mrSges(3,1) - t128 * mrSges(3,2) + mrSges(3,3) * t180 - m(4) * t161 - t108 * mrSges(4,1) - t119 * mrSges(4,3) - m(5) * (pkin(3) * t108 + t161) - (t108 * t157 + t185) * mrSges(5,1) - (-t108 * t153 + t119 * t157) * mrSges(5,2) + t190 * (pkin(4) * t185 - t107 * t159 + t108 * t139 + t161) + t189 * (t108 * t140 - t119 * t141) + t188 * (t108 * t141 + t119 * t140) + t187 * t107) * g(2) + (-mrSges(1,1) - t149 * mrSges(2,1) + t146 * mrSges(2,2) - m(3) * t173 - t131 * mrSges(3,1) - t130 * mrSges(3,2) - mrSges(3,3) * t181 - m(4) * t165 - t110 * mrSges(4,1) - t120 * mrSges(4,3) - m(5) * (pkin(3) * t110 + t165) - (t110 * t157 + t184) * mrSges(5,1) - (-t110 * t153 + t120 * t157) * mrSges(5,2) + t190 * (pkin(4) * t184 - t109 * t159 + t110 * t139 + t165) + t189 * (t110 * t140 - t120 * t141) + t188 * (t110 * t141 + t120 * t140) + t187 * t109) * g(1);
U  = t1;
