% Calculate potential energy for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:55:56
% EndTime: 2019-03-08 18:55:57
% DurationCPUTime: 0.75s
% Computational Cost: add. (600->109), mult. (1541->147), div. (0->0), fcn. (1949->14), ass. (0->58)
t150 = sin(pkin(7));
t154 = cos(pkin(7));
t155 = cos(pkin(6));
t151 = sin(pkin(6));
t152 = cos(pkin(12));
t184 = t151 * t152;
t169 = -t150 * t184 + t154 * t155;
t148 = sin(pkin(12));
t153 = cos(pkin(11));
t149 = sin(pkin(11));
t185 = t149 * t155;
t136 = -t148 * t153 - t152 * t185;
t182 = t151 * t154;
t170 = -t136 * t150 + t149 * t182;
t195 = -m(6) - m(7);
t194 = mrSges(4,2) - mrSges(5,3);
t193 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t192 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t191 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t190 = cos(qJ(3));
t189 = cos(qJ(4));
t187 = t148 * t151;
t186 = t149 * t151;
t183 = t151 * t153;
t181 = t153 * t155;
t179 = t155 * qJ(2) + qJ(1);
t178 = t153 * pkin(1) + qJ(2) * t186;
t175 = t150 * t190;
t174 = t154 * t190;
t173 = t151 * t175;
t172 = t149 * pkin(1) - qJ(2) * t183;
t134 = -t148 * t149 + t152 * t181;
t171 = -t134 * t150 - t153 * t182;
t137 = -t148 * t185 + t152 * t153;
t168 = t137 * pkin(2) + t170 * pkin(8) + t178;
t167 = pkin(2) * t187 + t169 * pkin(8) + t179;
t158 = sin(qJ(3));
t119 = -t136 * t174 + t137 * t158 - t149 * t173;
t120 = t137 * t190 + (t136 * t154 + t150 * t186) * t158;
t166 = t120 * pkin(3) + pkin(9) * t119 + t168;
t127 = -t155 * t175 + t158 * t187 - t174 * t184;
t128 = t155 * t150 * t158 + (t152 * t154 * t158 + t148 * t190) * t151;
t165 = t128 * pkin(3) + pkin(9) * t127 + t167;
t135 = t148 * t181 + t149 * t152;
t164 = t135 * pkin(2) + pkin(8) * t171 + t172;
t117 = -t134 * t174 + t135 * t158 + t153 * t173;
t118 = t135 * t190 + (t134 * t154 - t150 * t183) * t158;
t161 = t118 * pkin(3) + pkin(9) * t117 + t164;
t159 = cos(qJ(5));
t157 = sin(qJ(4));
t156 = sin(qJ(5));
t122 = t128 * t189 + t157 * t169;
t121 = t128 * t157 - t169 * t189;
t109 = t120 * t189 + t157 * t170;
t108 = t120 * t157 - t170 * t189;
t107 = t118 * t189 + t157 * t171;
t106 = t118 * t157 - t171 * t189;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t179 - t155 * mrSges(3,3) - (t148 * mrSges(3,1) + t152 * mrSges(3,2)) * t151 - m(4) * t167 - t128 * mrSges(4,1) - t169 * mrSges(4,3) - m(5) * t165 - t122 * mrSges(5,1) + t195 * (t122 * pkin(4) + pkin(10) * t121 + t165) + t194 * t127 + t192 * (t122 * t159 + t127 * t156) + t191 * (t122 * t156 - t127 * t159) + t193 * t121) * g(3) + (-m(3) * t172 - m(4) * t164 - m(5) * t161 - t149 * mrSges(2,1) - t135 * mrSges(3,1) - t118 * mrSges(4,1) - t107 * mrSges(5,1) - t153 * mrSges(2,2) - t134 * mrSges(3,2) + mrSges(3,3) * t183 - t171 * mrSges(4,3) - mrSges(1,2) + t195 * (t107 * pkin(4) + pkin(10) * t106 + t161) + t194 * t117 + t192 * (t107 * t159 + t117 * t156) + t191 * (t107 * t156 - t117 * t159) + t193 * t106) * g(2) + (-m(3) * t178 - m(4) * t168 - m(5) * t166 - t153 * mrSges(2,1) - t137 * mrSges(3,1) - t120 * mrSges(4,1) - t109 * mrSges(5,1) + t149 * mrSges(2,2) - t136 * mrSges(3,2) - mrSges(3,3) * t186 - t170 * mrSges(4,3) - mrSges(1,1) + t195 * (t109 * pkin(4) + pkin(10) * t108 + t166) + t194 * t119 + t192 * (t109 * t159 + t119 * t156) + t191 * (t109 * t156 - t119 * t159) + t193 * t108) * g(1);
U  = t1;
