% Calculate potential energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:54
% EndTime: 2019-03-08 22:35:55
% DurationCPUTime: 0.50s
% Computational Cost: add. (477->127), mult. (1189->177), div. (0->0), fcn. (1483->14), ass. (0->61)
t134 = sin(pkin(7));
t137 = cos(pkin(7));
t138 = cos(pkin(6));
t135 = sin(pkin(6));
t145 = cos(qJ(2));
t166 = t135 * t145;
t117 = -t134 * t166 + t138 * t137;
t133 = sin(pkin(12));
t136 = cos(pkin(12));
t142 = sin(qJ(2));
t162 = t138 * t145;
t120 = -t133 * t162 - t136 * t142;
t168 = t135 * t137;
t111 = -t120 * t134 + t133 * t168;
t176 = rSges(6,3) + pkin(10);
t175 = pkin(11) + rSges(7,3);
t174 = cos(qJ(3));
t173 = rSges(5,3) + qJ(4);
t171 = t133 * t135;
t141 = sin(qJ(3));
t170 = t134 * t141;
t169 = t135 * t136;
t167 = t135 * t142;
t165 = t137 * t141;
t163 = t138 * t142;
t161 = t138 * pkin(8) + qJ(1);
t160 = t136 * pkin(1) + pkin(8) * t171;
t157 = t134 * t174;
t156 = t137 * t174;
t155 = t135 * t157;
t118 = -t133 * t142 + t136 * t162;
t110 = -t118 * t134 - t136 * t168;
t121 = -t133 * t163 + t136 * t145;
t154 = t121 * pkin(2) + t111 * pkin(9) + t160;
t101 = t121 * t174 + (t120 * t137 + t134 * t171) * t141;
t153 = t101 * pkin(3) + t154;
t152 = pkin(2) * t167 + t117 * pkin(9) + t161;
t109 = t138 * t170 + (t142 * t174 + t145 * t165) * t135;
t151 = t109 * pkin(3) + t152;
t100 = -t120 * t156 + t121 * t141 - t133 * t155;
t150 = t111 * pkin(4) + t100 * qJ(4) + t153;
t108 = -t138 * t157 + t141 * t167 - t156 * t166;
t149 = t117 * pkin(4) + t108 * qJ(4) + t151;
t119 = t133 * t145 + t136 * t163;
t130 = t133 * pkin(1);
t148 = t119 * pkin(2) - pkin(8) * t169 + pkin(9) * t110 + t130;
t99 = t118 * t165 + t119 * t174 - t169 * t170;
t147 = t99 * pkin(3) + t148;
t98 = -t118 * t156 + t119 * t141 + t136 * t155;
t146 = t110 * pkin(4) + t98 * qJ(4) + t147;
t144 = cos(qJ(5));
t143 = cos(qJ(6));
t140 = sin(qJ(5));
t139 = sin(qJ(6));
t103 = t108 * t140 + t117 * t144;
t102 = -t108 * t144 + t117 * t140;
t93 = t100 * t140 + t111 * t144;
t92 = -t100 * t144 + t111 * t140;
t91 = t110 * t144 + t140 * t98;
t90 = t110 * t140 - t98 * t144;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t136 - rSges(2,2) * t133) + g(2) * (rSges(2,1) * t133 + rSges(2,2) * t136) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t121 + rSges(3,2) * t120 + t160) + g(2) * (rSges(3,1) * t119 + rSges(3,2) * t118 + t130) + g(3) * (t138 * rSges(3,3) + t161) + (g(1) * rSges(3,3) * t133 + g(3) * (rSges(3,1) * t142 + rSges(3,2) * t145) + g(2) * (-rSges(3,3) - pkin(8)) * t136) * t135) - m(4) * (g(1) * (rSges(4,1) * t101 - rSges(4,2) * t100 + rSges(4,3) * t111 + t154) + g(2) * (rSges(4,1) * t99 - rSges(4,2) * t98 + rSges(4,3) * t110 + t148) + g(3) * (t109 * rSges(4,1) - t108 * rSges(4,2) + t117 * rSges(4,3) + t152)) - m(5) * (g(1) * (rSges(5,1) * t111 - rSges(5,2) * t101 + t100 * t173 + t153) + g(2) * (rSges(5,1) * t110 - rSges(5,2) * t99 + t173 * t98 + t147) + g(3) * (t117 * rSges(5,1) - t109 * rSges(5,2) + t108 * t173 + t151)) - m(6) * (g(1) * (rSges(6,1) * t93 - rSges(6,2) * t92 + t101 * t176 + t150) + g(2) * (rSges(6,1) * t91 - rSges(6,2) * t90 + t176 * t99 + t146) + g(3) * (t103 * rSges(6,1) - t102 * rSges(6,2) + t109 * t176 + t149)) - m(7) * (g(1) * (t93 * pkin(5) + t101 * pkin(10) + (t101 * t139 + t143 * t93) * rSges(7,1) + (t101 * t143 - t139 * t93) * rSges(7,2) + t175 * t92 + t150) + g(2) * (t91 * pkin(5) + t99 * pkin(10) + (t139 * t99 + t143 * t91) * rSges(7,1) + (-t139 * t91 + t143 * t99) * rSges(7,2) + t175 * t90 + t146) + g(3) * (t103 * pkin(5) + t109 * pkin(10) + (t103 * t143 + t109 * t139) * rSges(7,1) + (-t103 * t139 + t109 * t143) * rSges(7,2) + t175 * t102 + t149));
U  = t1;
