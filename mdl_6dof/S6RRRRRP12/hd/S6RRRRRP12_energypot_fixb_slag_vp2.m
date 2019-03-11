% Calculate potential energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:32
% EndTime: 2019-03-10 02:59:33
% DurationCPUTime: 0.74s
% Computational Cost: add. (600->109), mult. (1541->144), div. (0->0), fcn. (1949->14), ass. (0->59)
t151 = cos(pkin(6));
t156 = sin(qJ(1));
t158 = cos(qJ(2));
t180 = t156 * t158;
t155 = sin(qJ(2));
t159 = cos(qJ(1));
t182 = t155 * t159;
t136 = -t151 * t180 - t182;
t148 = sin(pkin(7));
t150 = cos(pkin(7));
t149 = sin(pkin(6));
t186 = t149 * t156;
t170 = -t136 * t148 + t150 * t186;
t185 = t149 * t158;
t169 = -t148 * t185 + t150 * t151;
t196 = -m(6) - m(7);
t195 = mrSges(4,2) - mrSges(5,3);
t194 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t193 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t192 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t191 = cos(qJ(3));
t190 = cos(qJ(4));
t189 = t151 * pkin(9) + pkin(8);
t187 = t149 * t155;
t184 = t149 * t159;
t181 = t156 * t155;
t179 = t158 * t159;
t178 = t159 * pkin(1) + pkin(9) * t186;
t175 = t148 * t191;
t174 = t150 * t191;
t173 = t149 * t175;
t172 = t156 * pkin(1) - pkin(9) * t184;
t134 = t151 * t179 - t181;
t171 = -t134 * t148 - t150 * t184;
t137 = -t151 * t181 + t179;
t168 = t137 * pkin(2) + t170 * pkin(10) + t178;
t167 = pkin(2) * t187 + t169 * pkin(10) + t189;
t154 = sin(qJ(3));
t121 = -t136 * t174 + t137 * t154 - t156 * t173;
t122 = t137 * t191 + (t136 * t150 + t148 * t186) * t154;
t166 = t122 * pkin(3) + pkin(11) * t121 + t168;
t127 = -t151 * t175 + t154 * t187 - t174 * t185;
t128 = t151 * t148 * t154 + (t150 * t154 * t158 + t155 * t191) * t149;
t165 = t128 * pkin(3) + pkin(11) * t127 + t167;
t135 = t151 * t182 + t180;
t164 = t135 * pkin(2) + pkin(10) * t171 + t172;
t119 = -t134 * t174 + t135 * t154 + t159 * t173;
t120 = t135 * t191 + (t134 * t150 - t148 * t184) * t154;
t161 = t120 * pkin(3) + t119 * pkin(11) + t164;
t157 = cos(qJ(5));
t153 = sin(qJ(4));
t152 = sin(qJ(5));
t118 = t128 * t190 + t153 * t169;
t117 = t128 * t153 - t169 * t190;
t111 = t122 * t190 + t153 * t170;
t110 = t122 * t153 - t170 * t190;
t109 = t120 * t190 + t153 * t171;
t108 = t120 * t153 - t171 * t190;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t189 - t151 * mrSges(3,3) - (t155 * mrSges(3,1) + t158 * mrSges(3,2)) * t149 - m(4) * t167 - t128 * mrSges(4,1) - t169 * mrSges(4,3) - m(5) * t165 - t118 * mrSges(5,1) + t196 * (t118 * pkin(4) + pkin(12) * t117 + t165) + t195 * t127 + t193 * (t118 * t157 + t127 * t152) + t192 * (t118 * t152 - t127 * t157) + t194 * t117) * g(3) + (-m(3) * t172 - m(4) * t164 - m(5) * t161 - t156 * mrSges(2,1) - t135 * mrSges(3,1) - t120 * mrSges(4,1) - t109 * mrSges(5,1) - t159 * mrSges(2,2) - t134 * mrSges(3,2) + mrSges(3,3) * t184 - t171 * mrSges(4,3) - mrSges(1,2) + t196 * (t109 * pkin(4) + t108 * pkin(12) + t161) + t195 * t119 + t193 * (t109 * t157 + t119 * t152) + t192 * (t109 * t152 - t119 * t157) + t194 * t108) * g(2) + (-m(3) * t178 - m(4) * t168 - m(5) * t166 - t159 * mrSges(2,1) - t137 * mrSges(3,1) - t122 * mrSges(4,1) - t111 * mrSges(5,1) + t156 * mrSges(2,2) - t136 * mrSges(3,2) - mrSges(3,3) * t186 - t170 * mrSges(4,3) - mrSges(1,1) + t196 * (t111 * pkin(4) + pkin(12) * t110 + t166) + t195 * t121 + t193 * (t111 * t157 + t121 * t152) + t192 * (t111 * t152 - t121 * t157) + t194 * t110) * g(1);
U  = t1;
