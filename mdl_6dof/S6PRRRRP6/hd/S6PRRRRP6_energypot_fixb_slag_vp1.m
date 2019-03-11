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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:27
% EndTime: 2019-03-09 00:28:28
% DurationCPUTime: 0.49s
% Computational Cost: add. (600->123), mult. (1532->172), div. (0->0), fcn. (1949->14), ass. (0->67)
t147 = sin(pkin(7));
t150 = cos(pkin(7));
t151 = cos(pkin(6));
t148 = sin(pkin(6));
t157 = cos(qJ(2));
t180 = t148 * t157;
t167 = -t147 * t180 + t151 * t150;
t146 = sin(pkin(12));
t149 = cos(pkin(12));
t155 = sin(qJ(2));
t177 = t151 * t157;
t134 = -t146 * t177 - t149 * t155;
t182 = t148 * t150;
t168 = -t134 * t147 + t146 * t182;
t192 = rSges(7,1) + pkin(5);
t191 = rSges(7,2) + pkin(11);
t190 = rSges(5,3) + pkin(10);
t189 = rSges(6,3) + pkin(11);
t188 = cos(qJ(3));
t187 = cos(qJ(4));
t186 = rSges(7,3) + qJ(6);
t184 = t146 * t148;
t183 = t148 * t149;
t181 = t148 * t155;
t178 = t151 * t155;
t176 = t151 * pkin(8) + qJ(1);
t175 = t149 * pkin(1) + pkin(8) * t184;
t172 = t147 * t188;
t171 = t150 * t188;
t170 = t148 * t172;
t132 = -t146 * t155 + t149 * t177;
t169 = -t132 * t147 - t149 * t182;
t135 = -t146 * t178 + t149 * t157;
t166 = t135 * pkin(2) + t168 * pkin(9) + t175;
t165 = pkin(2) * t181 + t167 * pkin(9) + t176;
t154 = sin(qJ(3));
t118 = t135 * t188 + (t134 * t150 + t147 * t184) * t154;
t164 = t118 * pkin(3) + t166;
t126 = t151 * t147 * t154 + (t150 * t154 * t157 + t155 * t188) * t148;
t163 = t126 * pkin(3) + t165;
t153 = sin(qJ(4));
t109 = t118 * t187 + t153 * t168;
t117 = -t134 * t171 + t135 * t154 - t146 * t170;
t162 = t109 * pkin(4) + pkin(10) * t117 + t164;
t120 = t126 * t187 + t153 * t167;
t125 = -t151 * t172 + t154 * t181 - t171 * t180;
t161 = t120 * pkin(4) + t125 * pkin(10) + t163;
t133 = t146 * t157 + t149 * t178;
t143 = t146 * pkin(1);
t160 = t133 * pkin(2) - pkin(8) * t183 + pkin(9) * t169 + t143;
t116 = t133 * t188 + (t132 * t150 - t147 * t183) * t154;
t159 = t116 * pkin(3) + t160;
t107 = t116 * t187 + t153 * t169;
t115 = -t132 * t171 + t133 * t154 + t149 * t170;
t158 = t107 * pkin(4) + pkin(10) * t115 + t159;
t156 = cos(qJ(5));
t152 = sin(qJ(5));
t119 = t126 * t153 - t167 * t187;
t108 = t118 * t153 - t168 * t187;
t106 = t116 * t153 - t169 * t187;
t105 = t120 * t156 + t125 * t152;
t104 = t120 * t152 - t125 * t156;
t101 = t109 * t156 + t117 * t152;
t100 = t109 * t152 - t117 * t156;
t99 = t107 * t156 + t115 * t152;
t98 = t107 * t152 - t115 * t156;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t149 - rSges(2,2) * t146) + g(2) * (rSges(2,1) * t146 + rSges(2,2) * t149) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t135 + rSges(3,2) * t134 + t175) + g(2) * (rSges(3,1) * t133 + rSges(3,2) * t132 + t143) + g(3) * (t151 * rSges(3,3) + t176) + (g(1) * rSges(3,3) * t146 + g(3) * (rSges(3,1) * t155 + rSges(3,2) * t157) + g(2) * (-rSges(3,3) - pkin(8)) * t149) * t148) - m(4) * (g(1) * (t118 * rSges(4,1) - t117 * rSges(4,2) + rSges(4,3) * t168 + t166) + g(2) * (t116 * rSges(4,1) - t115 * rSges(4,2) + rSges(4,3) * t169 + t160) + g(3) * (t126 * rSges(4,1) - t125 * rSges(4,2) + rSges(4,3) * t167 + t165)) - m(5) * (g(1) * (rSges(5,1) * t109 - rSges(5,2) * t108 + t117 * t190 + t164) + g(2) * (rSges(5,1) * t107 - rSges(5,2) * t106 + t115 * t190 + t159) + g(3) * (t120 * rSges(5,1) - t119 * rSges(5,2) + t125 * t190 + t163)) - m(6) * (g(1) * (rSges(6,1) * t101 - rSges(6,2) * t100 + t108 * t189 + t162) + g(2) * (rSges(6,1) * t99 - rSges(6,2) * t98 + t106 * t189 + t158) + g(3) * (t105 * rSges(6,1) - t104 * rSges(6,2) + t119 * t189 + t161)) - m(7) * (g(1) * (t186 * t100 + t101 * t192 + t191 * t108 + t162) + g(2) * (t191 * t106 + t186 * t98 + t192 * t99 + t158) + g(3) * (t186 * t104 + t105 * t192 + t191 * t119 + t161));
U  = t1;
