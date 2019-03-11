% Calculate potential energy for
% S6RPRRRP11
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:21
% EndTime: 2019-03-09 06:34:22
% DurationCPUTime: 0.50s
% Computational Cost: add. (551->129), mult. (1388->178), div. (0->0), fcn. (1753->14), ass. (0->64)
t135 = sin(pkin(12));
t140 = cos(pkin(6));
t147 = cos(qJ(1));
t138 = cos(pkin(12));
t145 = sin(qJ(1));
t164 = t145 * t138;
t122 = -t135 * t147 - t140 * t164;
t136 = sin(pkin(7));
t139 = cos(pkin(7));
t137 = sin(pkin(6));
t169 = t137 * t145;
t155 = -t122 * t136 + t139 * t169;
t170 = t137 * t138;
t154 = -t136 * t170 + t139 * t140;
t178 = rSges(5,3) + pkin(10);
t177 = rSges(6,3) + pkin(11);
t176 = cos(qJ(3));
t175 = cos(qJ(4));
t174 = t140 * qJ(2) + pkin(8);
t173 = rSges(7,3) + qJ(6) + pkin(11);
t171 = t135 * t137;
t168 = t137 * t147;
t166 = t140 * t147;
t165 = t145 * t135;
t163 = t147 * pkin(1) + qJ(2) * t169;
t142 = sin(qJ(5));
t160 = pkin(5) * t142 + pkin(10);
t159 = t136 * t176;
t158 = t139 * t176;
t157 = t137 * t159;
t120 = t138 * t166 - t165;
t156 = -t120 * t136 - t139 * t168;
t123 = t138 * t147 - t140 * t165;
t153 = t123 * pkin(2) + t155 * pkin(9) + t163;
t152 = pkin(2) * t171 + t154 * pkin(9) + t174;
t144 = sin(qJ(3));
t109 = t123 * t176 + (t122 * t139 + t136 * t169) * t144;
t151 = t109 * pkin(3) + t153;
t114 = t140 * t136 * t144 + (t138 * t139 * t144 + t176 * t135) * t137;
t150 = t114 * pkin(3) + t152;
t121 = t135 * t166 + t164;
t133 = t145 * pkin(1);
t149 = t121 * pkin(2) + t156 * pkin(9) - qJ(2) * t168 + t133;
t107 = t121 * t176 + (t120 * t139 - t136 * t168) * t144;
t148 = t107 * pkin(3) + t149;
t146 = cos(qJ(5));
t143 = sin(qJ(4));
t131 = pkin(5) * t146 + pkin(4);
t113 = -t140 * t159 + t144 * t171 - t158 * t170;
t108 = -t122 * t158 + t123 * t144 - t145 * t157;
t106 = -t120 * t158 + t121 * t144 + t147 * t157;
t105 = t114 * t175 + t154 * t143;
t104 = t114 * t143 - t154 * t175;
t101 = t109 * t175 + t155 * t143;
t100 = t109 * t143 - t155 * t175;
t99 = t107 * t175 + t156 * t143;
t98 = t107 * t143 - t156 * t175;
t97 = t105 * t146 + t113 * t142;
t96 = -t105 * t142 + t113 * t146;
t95 = t101 * t146 + t108 * t142;
t94 = -t101 * t142 + t108 * t146;
t93 = t106 * t142 + t146 * t99;
t92 = t106 * t146 - t142 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t147 - t145 * rSges(2,2)) + g(2) * (t145 * rSges(2,1) + rSges(2,2) * t147) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t123 + rSges(3,2) * t122 + t163) + g(2) * (t121 * rSges(3,1) + t120 * rSges(3,2) + t133) + g(3) * (rSges(3,3) * t140 + t174) + (g(1) * rSges(3,3) * t145 + g(3) * (rSges(3,1) * t135 + rSges(3,2) * t138) + g(2) * (-rSges(3,3) - qJ(2)) * t147) * t137) - m(4) * (g(1) * (t109 * rSges(4,1) - t108 * rSges(4,2) + t155 * rSges(4,3) + t153) + g(2) * (t107 * rSges(4,1) - t106 * rSges(4,2) + t156 * rSges(4,3) + t149) + g(3) * (t114 * rSges(4,1) - t113 * rSges(4,2) + t154 * rSges(4,3) + t152)) - m(5) * (g(1) * (rSges(5,1) * t101 - rSges(5,2) * t100 + t178 * t108 + t151) + g(2) * (t99 * rSges(5,1) - t98 * rSges(5,2) + t178 * t106 + t148) + g(3) * (rSges(5,1) * t105 - rSges(5,2) * t104 + t178 * t113 + t150)) - m(6) * (g(1) * (rSges(6,1) * t95 + rSges(6,2) * t94 + pkin(4) * t101 + pkin(10) * t108 + t177 * t100 + t151) + g(2) * (t93 * rSges(6,1) + t92 * rSges(6,2) + t99 * pkin(4) + t106 * pkin(10) + t177 * t98 + t148) + g(3) * (rSges(6,1) * t97 + rSges(6,2) * t96 + pkin(4) * t105 + pkin(10) * t113 + t177 * t104 + t150)) - m(7) * (g(1) * (rSges(7,1) * t95 + rSges(7,2) * t94 + t173 * t100 + t101 * t131 + t160 * t108 + t151) + g(2) * (t93 * rSges(7,1) + t92 * rSges(7,2) + t160 * t106 + t99 * t131 + t173 * t98 + t148) + g(3) * (rSges(7,1) * t97 + rSges(7,2) * t96 + t173 * t104 + t105 * t131 + t160 * t113 + t150));
U  = t1;
