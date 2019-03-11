% Calculate potential energy for
% S6PRRRRP5
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:05
% EndTime: 2019-03-09 00:20:05
% DurationCPUTime: 0.51s
% Computational Cost: add. (551->129), mult. (1388->180), div. (0->0), fcn. (1753->14), ass. (0->64)
t134 = sin(pkin(7));
t137 = cos(pkin(7));
t138 = cos(pkin(6));
t135 = sin(pkin(6));
t145 = cos(qJ(2));
t166 = t135 * t145;
t152 = -t134 * t166 + t138 * t137;
t133 = sin(pkin(12));
t136 = cos(pkin(12));
t143 = sin(qJ(2));
t163 = t138 * t145;
t120 = -t133 * t163 - t136 * t143;
t168 = t135 * t137;
t153 = -t120 * t134 + t133 * t168;
t176 = rSges(5,3) + pkin(10);
t175 = rSges(6,3) + pkin(11);
t174 = cos(qJ(3));
t173 = cos(qJ(4));
t172 = rSges(7,3) + qJ(6) + pkin(11);
t170 = t133 * t135;
t169 = t135 * t136;
t167 = t135 * t143;
t164 = t138 * t143;
t162 = t138 * pkin(8) + qJ(1);
t161 = t136 * pkin(1) + pkin(8) * t170;
t140 = sin(qJ(5));
t158 = pkin(5) * t140 + pkin(10);
t157 = t134 * t174;
t156 = t137 * t174;
t155 = t135 * t157;
t118 = -t133 * t143 + t136 * t163;
t154 = -t118 * t134 - t136 * t168;
t121 = -t133 * t164 + t136 * t145;
t151 = t121 * pkin(2) + t153 * pkin(9) + t161;
t150 = pkin(2) * t167 + t152 * pkin(9) + t162;
t142 = sin(qJ(3));
t105 = t121 * t174 + (t120 * t137 + t134 * t170) * t142;
t149 = t105 * pkin(3) + t151;
t112 = t138 * t134 * t142 + (t137 * t142 * t145 + t174 * t143) * t135;
t148 = t112 * pkin(3) + t150;
t119 = t133 * t145 + t136 * t164;
t130 = t133 * pkin(1);
t147 = t119 * pkin(2) - pkin(8) * t169 + t154 * pkin(9) + t130;
t103 = t119 * t174 + (t118 * t137 - t134 * t169) * t142;
t146 = t103 * pkin(3) + t147;
t144 = cos(qJ(5));
t141 = sin(qJ(4));
t129 = pkin(5) * t144 + pkin(4);
t111 = -t138 * t157 + t142 * t167 - t156 * t166;
t107 = t112 * t173 + t152 * t141;
t106 = t112 * t141 - t152 * t173;
t104 = -t120 * t156 + t121 * t142 - t133 * t155;
t102 = -t118 * t156 + t119 * t142 + t136 * t155;
t99 = t105 * t173 + t153 * t141;
t98 = t105 * t141 - t153 * t173;
t97 = t103 * t173 + t141 * t154;
t96 = t103 * t141 - t154 * t173;
t95 = t107 * t144 + t111 * t140;
t94 = -t107 * t140 + t111 * t144;
t93 = t104 * t140 + t144 * t99;
t92 = t104 * t144 - t140 * t99;
t91 = t102 * t140 + t144 * t97;
t90 = t102 * t144 - t140 * t97;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t136 - rSges(2,2) * t133) + g(2) * (rSges(2,1) * t133 + rSges(2,2) * t136) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t121 + rSges(3,2) * t120 + t161) + g(2) * (rSges(3,1) * t119 + rSges(3,2) * t118 + t130) + g(3) * (t138 * rSges(3,3) + t162) + (g(1) * rSges(3,3) * t133 + g(3) * (rSges(3,1) * t143 + rSges(3,2) * t145) + g(2) * (-rSges(3,3) - pkin(8)) * t136) * t135) - m(4) * (g(1) * (t105 * rSges(4,1) - t104 * rSges(4,2) + rSges(4,3) * t153 + t151) + g(2) * (t103 * rSges(4,1) - t102 * rSges(4,2) + rSges(4,3) * t154 + t147) + g(3) * (t112 * rSges(4,1) - t111 * rSges(4,2) + rSges(4,3) * t152 + t150)) - m(5) * (g(1) * (rSges(5,1) * t99 - rSges(5,2) * t98 + t104 * t176 + t149) + g(2) * (rSges(5,1) * t97 - rSges(5,2) * t96 + t102 * t176 + t146) + g(3) * (t107 * rSges(5,1) - t106 * rSges(5,2) + t111 * t176 + t148)) - m(6) * (g(1) * (rSges(6,1) * t93 + rSges(6,2) * t92 + pkin(4) * t99 + pkin(10) * t104 + t175 * t98 + t149) + g(2) * (rSges(6,1) * t91 + rSges(6,2) * t90 + pkin(4) * t97 + pkin(10) * t102 + t175 * t96 + t146) + g(3) * (t95 * rSges(6,1) + t94 * rSges(6,2) + t107 * pkin(4) + t111 * pkin(10) + t106 * t175 + t148)) - m(7) * (g(1) * (rSges(7,1) * t93 + rSges(7,2) * t92 + t104 * t158 + t129 * t99 + t172 * t98 + t149) + g(2) * (rSges(7,1) * t91 + rSges(7,2) * t90 + t102 * t158 + t129 * t97 + t172 * t96 + t146) + g(3) * (t95 * rSges(7,1) + t94 * rSges(7,2) + t106 * t172 + t107 * t129 + t111 * t158 + t148));
U  = t1;
