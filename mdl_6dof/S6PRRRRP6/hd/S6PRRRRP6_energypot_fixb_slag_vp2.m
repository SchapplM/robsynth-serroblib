% Calculate potential energy for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:38
% EndTime: 2019-03-09 00:28:38
% DurationCPUTime: 0.71s
% Computational Cost: add. (600->109), mult. (1541->147), div. (0->0), fcn. (1949->14), ass. (0->58)
t148 = sin(pkin(7));
t151 = cos(pkin(7));
t152 = cos(pkin(6));
t149 = sin(pkin(6));
t158 = cos(qJ(2));
t182 = t149 * t158;
t168 = -t148 * t182 + t152 * t151;
t147 = sin(pkin(12));
t150 = cos(pkin(12));
t156 = sin(qJ(2));
t179 = t152 * t158;
t135 = -t147 * t179 - t150 * t156;
t184 = t149 * t151;
t169 = -t135 * t148 + t147 * t184;
t194 = -m(6) - m(7);
t193 = mrSges(4,2) - mrSges(5,3);
t192 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t191 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t190 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t189 = cos(qJ(3));
t188 = cos(qJ(4));
t186 = t147 * t149;
t185 = t149 * t150;
t183 = t149 * t156;
t180 = t152 * t156;
t178 = t150 * pkin(1) + pkin(8) * t186;
t177 = t152 * pkin(8) + qJ(1);
t174 = t148 * t189;
t173 = t151 * t189;
t172 = t149 * t174;
t171 = t147 * pkin(1) - pkin(8) * t185;
t133 = -t147 * t156 + t150 * t179;
t170 = -t133 * t148 - t150 * t184;
t136 = -t147 * t180 + t150 * t158;
t167 = t136 * pkin(2) + t169 * pkin(9) + t178;
t166 = pkin(2) * t183 + t168 * pkin(9) + t177;
t155 = sin(qJ(3));
t118 = -t135 * t173 + t136 * t155 - t147 * t172;
t119 = t136 * t189 + (t135 * t151 + t148 * t186) * t155;
t165 = t119 * pkin(3) + pkin(10) * t118 + t167;
t126 = -t152 * t174 + t155 * t183 - t173 * t182;
t127 = t152 * t148 * t155 + (t151 * t155 * t158 + t156 * t189) * t149;
t164 = t127 * pkin(3) + t126 * pkin(10) + t166;
t134 = t147 * t158 + t150 * t180;
t163 = t134 * pkin(2) + pkin(9) * t170 + t171;
t116 = -t133 * t173 + t134 * t155 + t150 * t172;
t117 = t134 * t189 + (t133 * t151 - t148 * t185) * t155;
t160 = t117 * pkin(3) + pkin(10) * t116 + t163;
t157 = cos(qJ(5));
t154 = sin(qJ(4));
t153 = sin(qJ(5));
t121 = t127 * t188 + t154 * t168;
t120 = t127 * t154 - t168 * t188;
t110 = t119 * t188 + t154 * t169;
t109 = t119 * t154 - t169 * t188;
t108 = t117 * t188 + t154 * t170;
t107 = t117 * t154 - t170 * t188;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t177 - t152 * mrSges(3,3) - (t156 * mrSges(3,1) + t158 * mrSges(3,2)) * t149 - m(4) * t166 - t127 * mrSges(4,1) - t168 * mrSges(4,3) - m(5) * t164 - t121 * mrSges(5,1) + t194 * (t121 * pkin(4) + t120 * pkin(11) + t164) + t193 * t126 + t191 * (t121 * t157 + t126 * t153) + t190 * (t121 * t153 - t126 * t157) + t192 * t120) * g(3) + (-m(3) * t171 - m(4) * t163 - m(5) * t160 - t147 * mrSges(2,1) - t134 * mrSges(3,1) - t117 * mrSges(4,1) - t108 * mrSges(5,1) - t150 * mrSges(2,2) - t133 * mrSges(3,2) + mrSges(3,3) * t185 - t170 * mrSges(4,3) - mrSges(1,2) + t194 * (t108 * pkin(4) + pkin(11) * t107 + t160) + t190 * (t108 * t153 - t116 * t157) + t193 * t116 + t191 * (t108 * t157 + t116 * t153) + t192 * t107) * g(2) + (-m(3) * t178 - m(4) * t167 - m(5) * t165 - t150 * mrSges(2,1) - t136 * mrSges(3,1) - t119 * mrSges(4,1) - t110 * mrSges(5,1) + t147 * mrSges(2,2) - t135 * mrSges(3,2) - mrSges(3,3) * t186 - t169 * mrSges(4,3) - mrSges(1,1) + t194 * (t110 * pkin(4) + pkin(11) * t109 + t165) + t193 * t118 + t191 * (t110 * t157 + t118 * t153) + t190 * (t110 * t153 - t118 * t157) + t192 * t109) * g(1);
U  = t1;
