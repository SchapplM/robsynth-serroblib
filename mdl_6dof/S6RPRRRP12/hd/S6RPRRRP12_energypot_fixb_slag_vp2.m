% Calculate potential energy for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:03
% EndTime: 2019-03-09 06:44:04
% DurationCPUTime: 0.71s
% Computational Cost: add. (600->109), mult. (1541->145), div. (0->0), fcn. (1949->14), ass. (0->58)
t148 = sin(pkin(12));
t153 = cos(pkin(6));
t159 = cos(qJ(1));
t151 = cos(pkin(12));
t157 = sin(qJ(1));
t179 = t157 * t151;
t136 = -t148 * t159 - t153 * t179;
t149 = sin(pkin(7));
t152 = cos(pkin(7));
t150 = sin(pkin(6));
t184 = t150 * t157;
t170 = -t136 * t149 + t152 * t184;
t185 = t150 * t151;
t169 = -t149 * t185 + t152 * t153;
t195 = -m(6) - m(7);
t194 = mrSges(4,2) - mrSges(5,3);
t193 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t192 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t191 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t190 = cos(qJ(3));
t189 = cos(qJ(4));
t188 = t153 * qJ(2) + pkin(8);
t186 = t148 * t150;
t183 = t150 * t159;
t181 = t153 * t159;
t180 = t157 * t148;
t178 = t159 * pkin(1) + qJ(2) * t184;
t175 = t149 * t190;
t174 = t152 * t190;
t173 = t150 * t175;
t172 = t157 * pkin(1) - qJ(2) * t183;
t134 = t151 * t181 - t180;
t171 = -t134 * t149 - t152 * t183;
t137 = t151 * t159 - t153 * t180;
t168 = t137 * pkin(2) + t170 * pkin(9) + t178;
t167 = pkin(2) * t186 + t169 * pkin(9) + t188;
t156 = sin(qJ(3));
t121 = -t136 * t174 + t137 * t156 - t157 * t173;
t122 = t137 * t190 + (t136 * t152 + t149 * t184) * t156;
t166 = t122 * pkin(3) + pkin(10) * t121 + t168;
t127 = -t153 * t175 + t156 * t186 - t174 * t185;
t128 = t153 * t149 * t156 + (t151 * t152 * t156 + t148 * t190) * t150;
t165 = t128 * pkin(3) + pkin(10) * t127 + t167;
t135 = t148 * t181 + t179;
t164 = t135 * pkin(2) + pkin(9) * t171 + t172;
t119 = -t134 * t174 + t135 * t156 + t159 * t173;
t120 = t135 * t190 + (t134 * t152 - t149 * t183) * t156;
t161 = t120 * pkin(3) + t119 * pkin(10) + t164;
t158 = cos(qJ(5));
t155 = sin(qJ(4));
t154 = sin(qJ(5));
t118 = t128 * t189 + t155 * t169;
t117 = t128 * t155 - t169 * t189;
t111 = t122 * t189 + t155 * t170;
t110 = t122 * t155 - t170 * t189;
t109 = t120 * t189 + t155 * t171;
t108 = t120 * t155 - t171 * t189;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t188 - t153 * mrSges(3,3) - (t148 * mrSges(3,1) + t151 * mrSges(3,2)) * t150 - m(4) * t167 - t128 * mrSges(4,1) - t169 * mrSges(4,3) - m(5) * t165 - t118 * mrSges(5,1) + t195 * (t118 * pkin(4) + pkin(11) * t117 + t165) + t194 * t127 + t192 * (t118 * t158 + t127 * t154) + t191 * (t118 * t154 - t127 * t158) + t193 * t117) * g(3) + (-m(3) * t172 - m(4) * t164 - m(5) * t161 - t157 * mrSges(2,1) - t135 * mrSges(3,1) - t120 * mrSges(4,1) - t109 * mrSges(5,1) - t159 * mrSges(2,2) - t134 * mrSges(3,2) + mrSges(3,3) * t183 - t171 * mrSges(4,3) - mrSges(1,2) + t195 * (t109 * pkin(4) + t108 * pkin(11) + t161) + t194 * t119 + t192 * (t109 * t158 + t119 * t154) + t191 * (t109 * t154 - t119 * t158) + t193 * t108) * g(2) + (-m(3) * t178 - m(4) * t168 - m(5) * t166 - t159 * mrSges(2,1) - t137 * mrSges(3,1) - t122 * mrSges(4,1) - t111 * mrSges(5,1) + t157 * mrSges(2,2) - t136 * mrSges(3,2) - mrSges(3,3) * t184 - t170 * mrSges(4,3) - mrSges(1,1) + t195 * (t111 * pkin(4) + pkin(11) * t110 + t166) + t194 * t121 + t192 * (t111 * t158 + t121 * t154) + t191 * (t111 * t154 - t121 * t158) + t193 * t110) * g(1);
U  = t1;
