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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:43:53
% EndTime: 2019-03-09 06:43:53
% DurationCPUTime: 0.47s
% Computational Cost: add. (600->123), mult. (1532->170), div. (0->0), fcn. (1949->14), ass. (0->67)
t148 = sin(pkin(12));
t153 = cos(pkin(6));
t159 = cos(qJ(1));
t151 = cos(pkin(12));
t157 = sin(qJ(1));
t178 = t157 * t151;
t136 = -t148 * t159 - t153 * t178;
t149 = sin(pkin(7));
t152 = cos(pkin(7));
t150 = sin(pkin(6));
t183 = t150 * t157;
t170 = -t136 * t149 + t152 * t183;
t184 = t150 * t151;
t169 = -t149 * t184 + t152 * t153;
t194 = rSges(7,1) + pkin(5);
t193 = rSges(7,2) + pkin(11);
t192 = rSges(5,3) + pkin(10);
t191 = rSges(6,3) + pkin(11);
t190 = cos(qJ(3));
t189 = cos(qJ(4));
t188 = t153 * qJ(2) + pkin(8);
t187 = rSges(7,3) + qJ(6);
t185 = t148 * t150;
t182 = t150 * t159;
t180 = t153 * t159;
t179 = t157 * t148;
t177 = t159 * pkin(1) + qJ(2) * t183;
t174 = t149 * t190;
t173 = t152 * t190;
t172 = t150 * t174;
t134 = t151 * t180 - t179;
t171 = -t134 * t149 - t152 * t182;
t137 = t151 * t159 - t153 * t179;
t168 = t137 * pkin(2) + t170 * pkin(9) + t177;
t167 = pkin(2) * t185 + t169 * pkin(9) + t188;
t156 = sin(qJ(3));
t122 = t137 * t190 + (t136 * t152 + t149 * t183) * t156;
t166 = t122 * pkin(3) + t168;
t128 = t153 * t149 * t156 + (t151 * t152 * t156 + t190 * t148) * t150;
t165 = t128 * pkin(3) + t167;
t155 = sin(qJ(4));
t111 = t122 * t189 + t170 * t155;
t121 = -t136 * t173 + t137 * t156 - t157 * t172;
t164 = t111 * pkin(4) + pkin(10) * t121 + t166;
t118 = t128 * t189 + t169 * t155;
t127 = -t153 * t174 + t156 * t185 - t173 * t184;
t163 = t118 * pkin(4) + pkin(10) * t127 + t165;
t135 = t148 * t180 + t178;
t146 = t157 * pkin(1);
t162 = t135 * pkin(2) + t171 * pkin(9) - qJ(2) * t182 + t146;
t120 = t135 * t190 + (t134 * t152 - t149 * t182) * t156;
t161 = t120 * pkin(3) + t162;
t109 = t120 * t189 + t171 * t155;
t119 = -t134 * t173 + t135 * t156 + t159 * t172;
t160 = t109 * pkin(4) + t119 * pkin(10) + t161;
t158 = cos(qJ(5));
t154 = sin(qJ(5));
t117 = t128 * t155 - t169 * t189;
t110 = t122 * t155 - t170 * t189;
t108 = t120 * t155 - t171 * t189;
t105 = t118 * t158 + t127 * t154;
t104 = t118 * t154 - t127 * t158;
t103 = t111 * t158 + t121 * t154;
t102 = t111 * t154 - t121 * t158;
t101 = t109 * t158 + t119 * t154;
t100 = t109 * t154 - t119 * t158;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t159 - t157 * rSges(2,2)) + g(2) * (t157 * rSges(2,1) + rSges(2,2) * t159) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t137 + rSges(3,2) * t136 + t177) + g(2) * (t135 * rSges(3,1) + t134 * rSges(3,2) + t146) + g(3) * (rSges(3,3) * t153 + t188) + (g(1) * rSges(3,3) * t157 + g(3) * (rSges(3,1) * t148 + rSges(3,2) * t151) + g(2) * (-rSges(3,3) - qJ(2)) * t159) * t150) - m(4) * (g(1) * (t122 * rSges(4,1) - t121 * rSges(4,2) + t170 * rSges(4,3) + t168) + g(2) * (t120 * rSges(4,1) - t119 * rSges(4,2) + t171 * rSges(4,3) + t162) + g(3) * (t128 * rSges(4,1) - t127 * rSges(4,2) + t169 * rSges(4,3) + t167)) - m(5) * (g(1) * (rSges(5,1) * t111 - rSges(5,2) * t110 + t192 * t121 + t166) + g(2) * (t109 * rSges(5,1) - t108 * rSges(5,2) + t192 * t119 + t161) + g(3) * (rSges(5,1) * t118 - rSges(5,2) * t117 + t192 * t127 + t165)) - m(6) * (g(1) * (rSges(6,1) * t103 - rSges(6,2) * t102 + t191 * t110 + t164) + g(2) * (t101 * rSges(6,1) - t100 * rSges(6,2) + t191 * t108 + t160) + g(3) * (rSges(6,1) * t105 - rSges(6,2) * t104 + t191 * t117 + t163)) - m(7) * (g(1) * (t187 * t102 + t194 * t103 + t193 * t110 + t164) + g(2) * (t187 * t100 + t194 * t101 + t193 * t108 + t160) + g(3) * (t187 * t104 + t194 * t105 + t193 * t117 + t163));
U  = t1;
