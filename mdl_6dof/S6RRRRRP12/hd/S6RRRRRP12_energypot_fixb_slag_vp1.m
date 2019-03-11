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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:32
% EndTime: 2019-03-10 02:59:32
% DurationCPUTime: 0.50s
% Computational Cost: add. (600->123), mult. (1532->169), div. (0->0), fcn. (1949->14), ass. (0->68)
t147 = sin(pkin(7));
t149 = cos(pkin(7));
t150 = cos(pkin(6));
t148 = sin(pkin(6));
t157 = cos(qJ(2));
t183 = t148 * t157;
t168 = -t147 * t183 + t149 * t150;
t155 = sin(qJ(1));
t178 = t155 * t157;
t154 = sin(qJ(2));
t158 = cos(qJ(1));
t180 = t154 * t158;
t135 = -t150 * t178 - t180;
t184 = t148 * t155;
t169 = -t135 * t147 + t149 * t184;
t194 = rSges(7,1) + pkin(5);
t193 = rSges(7,2) + pkin(12);
t192 = rSges(5,3) + pkin(11);
t191 = rSges(6,3) + pkin(12);
t190 = cos(qJ(3));
t189 = cos(qJ(4));
t188 = t150 * pkin(9) + pkin(8);
t187 = rSges(7,3) + qJ(6);
t185 = t148 * t154;
t182 = t148 * t158;
t179 = t155 * t154;
t177 = t157 * t158;
t176 = t158 * pkin(1) + pkin(9) * t184;
t173 = t147 * t190;
t172 = t149 * t190;
t171 = t148 * t173;
t133 = t150 * t177 - t179;
t170 = -t133 * t147 - t149 * t182;
t136 = -t150 * t179 + t177;
t167 = t136 * pkin(2) + t169 * pkin(10) + t176;
t166 = pkin(2) * t185 + t168 * pkin(10) + t188;
t153 = sin(qJ(3));
t121 = t136 * t190 + (t135 * t149 + t147 * t184) * t153;
t165 = t121 * pkin(3) + t167;
t127 = t150 * t147 * t153 + (t149 * t153 * t157 + t154 * t190) * t148;
t164 = t127 * pkin(3) + t166;
t152 = sin(qJ(4));
t110 = t121 * t189 + t152 * t169;
t120 = -t135 * t172 + t136 * t153 - t155 * t171;
t163 = t110 * pkin(4) + pkin(11) * t120 + t165;
t117 = t127 * t189 + t152 * t168;
t126 = -t150 * t173 + t153 * t185 - t172 * t183;
t162 = t117 * pkin(4) + pkin(11) * t126 + t164;
t134 = t150 * t180 + t178;
t145 = t155 * pkin(1);
t161 = t134 * pkin(2) - pkin(9) * t182 + pkin(10) * t170 + t145;
t119 = t134 * t190 + (t133 * t149 - t147 * t182) * t153;
t160 = t119 * pkin(3) + t161;
t108 = t119 * t189 + t152 * t170;
t118 = -t133 * t172 + t134 * t153 + t158 * t171;
t159 = t108 * pkin(4) + t118 * pkin(11) + t160;
t156 = cos(qJ(5));
t151 = sin(qJ(5));
t116 = t127 * t152 - t168 * t189;
t109 = t121 * t152 - t169 * t189;
t107 = t119 * t152 - t170 * t189;
t104 = t117 * t156 + t126 * t151;
t103 = t117 * t151 - t126 * t156;
t102 = t110 * t156 + t120 * t151;
t101 = t110 * t151 - t120 * t156;
t100 = t108 * t156 + t118 * t151;
t99 = t108 * t151 - t118 * t156;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t158 - t155 * rSges(2,2)) + g(2) * (t155 * rSges(2,1) + rSges(2,2) * t158) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t136 + rSges(3,2) * t135 + t176) + g(2) * (t134 * rSges(3,1) + t133 * rSges(3,2) + t145) + g(3) * (rSges(3,3) * t150 + t188) + (g(1) * rSges(3,3) * t155 + g(3) * (rSges(3,1) * t154 + rSges(3,2) * t157) + g(2) * (-rSges(3,3) - pkin(9)) * t158) * t148) - m(4) * (g(1) * (t121 * rSges(4,1) - t120 * rSges(4,2) + rSges(4,3) * t169 + t167) + g(2) * (t119 * rSges(4,1) - t118 * rSges(4,2) + rSges(4,3) * t170 + t161) + g(3) * (t127 * rSges(4,1) - t126 * rSges(4,2) + rSges(4,3) * t168 + t166)) - m(5) * (g(1) * (rSges(5,1) * t110 - rSges(5,2) * t109 + t120 * t192 + t165) + g(2) * (t108 * rSges(5,1) - t107 * rSges(5,2) + t118 * t192 + t160) + g(3) * (rSges(5,1) * t117 - rSges(5,2) * t116 + t126 * t192 + t164)) - m(6) * (g(1) * (rSges(6,1) * t102 - rSges(6,2) * t101 + t109 * t191 + t163) + g(2) * (t100 * rSges(6,1) - t99 * rSges(6,2) + t107 * t191 + t159) + g(3) * (rSges(6,1) * t104 - rSges(6,2) * t103 + t116 * t191 + t162)) - m(7) * (g(1) * (t187 * t101 + t102 * t194 + t193 * t109 + t163) + g(2) * (t100 * t194 + t193 * t107 + t187 * t99 + t159) + g(3) * (t187 * t103 + t104 * t194 + t193 * t116 + t162));
U  = t1;
