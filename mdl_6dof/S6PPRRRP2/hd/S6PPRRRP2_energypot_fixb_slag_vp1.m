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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:55:56
% EndTime: 2019-03-08 18:55:56
% DurationCPUTime: 0.51s
% Computational Cost: add. (600->123), mult. (1532->172), div. (0->0), fcn. (1949->14), ass. (0->67)
t148 = sin(pkin(7));
t152 = cos(pkin(7));
t153 = cos(pkin(6));
t149 = sin(pkin(6));
t150 = cos(pkin(12));
t181 = t149 * t150;
t167 = -t148 * t181 + t152 * t153;
t146 = sin(pkin(12));
t151 = cos(pkin(11));
t147 = sin(pkin(11));
t182 = t147 * t153;
t134 = -t146 * t151 - t150 * t182;
t179 = t149 * t152;
t168 = -t134 * t148 + t147 * t179;
t192 = rSges(7,1) + pkin(5);
t191 = rSges(7,2) + pkin(10);
t190 = rSges(5,3) + pkin(9);
t189 = rSges(6,3) + pkin(10);
t188 = cos(qJ(3));
t187 = cos(qJ(4));
t186 = rSges(7,3) + qJ(6);
t184 = t146 * t149;
t183 = t147 * t149;
t180 = t149 * t151;
t178 = t151 * t153;
t176 = t153 * qJ(2) + qJ(1);
t175 = t151 * pkin(1) + qJ(2) * t183;
t172 = t148 * t188;
t171 = t152 * t188;
t170 = t149 * t172;
t132 = -t146 * t147 + t150 * t178;
t169 = -t132 * t148 - t151 * t179;
t135 = -t146 * t182 + t150 * t151;
t166 = t135 * pkin(2) + t168 * pkin(8) + t175;
t165 = pkin(2) * t184 + t167 * pkin(8) + t176;
t156 = sin(qJ(3));
t118 = t135 * t188 + (t134 * t152 + t148 * t183) * t156;
t164 = t118 * pkin(3) + t166;
t126 = t153 * t148 * t156 + (t150 * t152 * t156 + t188 * t146) * t149;
t163 = t126 * pkin(3) + t165;
t155 = sin(qJ(4));
t107 = t118 * t187 + t168 * t155;
t117 = -t134 * t171 + t135 * t156 - t147 * t170;
t162 = t107 * pkin(4) + pkin(9) * t117 + t164;
t120 = t126 * t187 + t167 * t155;
t125 = -t153 * t172 + t156 * t184 - t171 * t181;
t161 = t120 * pkin(4) + pkin(9) * t125 + t163;
t133 = t146 * t178 + t147 * t150;
t144 = t147 * pkin(1);
t160 = t133 * pkin(2) + t169 * pkin(8) - qJ(2) * t180 + t144;
t116 = t133 * t188 + (t132 * t152 - t148 * t180) * t156;
t159 = t116 * pkin(3) + t160;
t105 = t116 * t187 + t169 * t155;
t115 = -t132 * t171 + t133 * t156 + t151 * t170;
t158 = t105 * pkin(4) + pkin(9) * t115 + t159;
t157 = cos(qJ(5));
t154 = sin(qJ(5));
t119 = t126 * t155 - t167 * t187;
t109 = t120 * t157 + t125 * t154;
t108 = t120 * t154 - t125 * t157;
t106 = t118 * t155 - t168 * t187;
t104 = t116 * t155 - t169 * t187;
t101 = t107 * t157 + t117 * t154;
t100 = t107 * t154 - t117 * t157;
t99 = t105 * t157 + t115 * t154;
t98 = t105 * t154 - t115 * t157;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t151 - rSges(2,2) * t147) + g(2) * (rSges(2,1) * t147 + rSges(2,2) * t151) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t135 + rSges(3,2) * t134 + t175) + g(2) * (rSges(3,1) * t133 + rSges(3,2) * t132 + t144) + g(3) * (rSges(3,3) * t153 + t176) + (g(1) * rSges(3,3) * t147 + g(3) * (rSges(3,1) * t146 + rSges(3,2) * t150) + g(2) * (-rSges(3,3) - qJ(2)) * t151) * t149) - m(4) * (g(1) * (t118 * rSges(4,1) - t117 * rSges(4,2) + t168 * rSges(4,3) + t166) + g(2) * (t116 * rSges(4,1) - t115 * rSges(4,2) + t169 * rSges(4,3) + t160) + g(3) * (t126 * rSges(4,1) - t125 * rSges(4,2) + t167 * rSges(4,3) + t165)) - m(5) * (g(1) * (rSges(5,1) * t107 - rSges(5,2) * t106 + t190 * t117 + t164) + g(2) * (rSges(5,1) * t105 - rSges(5,2) * t104 + t190 * t115 + t159) + g(3) * (rSges(5,1) * t120 - rSges(5,2) * t119 + t190 * t125 + t163)) - m(6) * (g(1) * (rSges(6,1) * t101 - rSges(6,2) * t100 + t189 * t106 + t162) + g(2) * (rSges(6,1) * t99 - rSges(6,2) * t98 + t189 * t104 + t158) + g(3) * (rSges(6,1) * t109 - rSges(6,2) * t108 + t189 * t119 + t161)) - m(7) * (g(1) * (t186 * t100 + t192 * t101 + t191 * t106 + t162) + g(2) * (t191 * t104 + t186 * t98 + t192 * t99 + t158) + g(3) * (t186 * t108 + t192 * t109 + t191 * t119 + t161));
U  = t1;
