% Calculate potential energy for
% S6PPRRRP1
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:14
% EndTime: 2019-03-08 18:52:14
% DurationCPUTime: 0.51s
% Computational Cost: add. (551->129), mult. (1388->180), div. (0->0), fcn. (1753->14), ass. (0->64)
t135 = sin(pkin(7));
t139 = cos(pkin(7));
t140 = cos(pkin(6));
t136 = sin(pkin(6));
t137 = cos(pkin(12));
t167 = t136 * t137;
t152 = -t135 * t167 + t140 * t139;
t133 = sin(pkin(12));
t138 = cos(pkin(11));
t134 = sin(pkin(11));
t168 = t134 * t140;
t120 = -t133 * t138 - t137 * t168;
t165 = t136 * t139;
t153 = -t120 * t135 + t134 * t165;
t176 = rSges(5,3) + pkin(9);
t175 = rSges(6,3) + pkin(10);
t174 = cos(qJ(3));
t173 = cos(qJ(4));
t172 = rSges(7,3) + qJ(6) + pkin(10);
t170 = t133 * t136;
t169 = t134 * t136;
t166 = t136 * t138;
t164 = t138 * t140;
t162 = t140 * qJ(2) + qJ(1);
t161 = t138 * pkin(1) + qJ(2) * t169;
t142 = sin(qJ(5));
t158 = pkin(5) * t142 + pkin(9);
t157 = t135 * t174;
t156 = t139 * t174;
t155 = t136 * t157;
t118 = -t133 * t134 + t137 * t164;
t154 = -t118 * t135 - t138 * t165;
t121 = -t133 * t168 + t137 * t138;
t151 = t121 * pkin(2) + t153 * pkin(8) + t161;
t150 = pkin(2) * t170 + t152 * pkin(8) + t162;
t144 = sin(qJ(3));
t105 = t121 * t174 + (t120 * t139 + t135 * t169) * t144;
t149 = t105 * pkin(3) + t151;
t112 = t140 * t135 * t144 + (t137 * t139 * t144 + t133 * t174) * t136;
t148 = t112 * pkin(3) + t150;
t119 = t133 * t164 + t134 * t137;
t131 = t134 * pkin(1);
t147 = t119 * pkin(2) + t154 * pkin(8) - qJ(2) * t166 + t131;
t103 = t119 * t174 + (t118 * t139 - t135 * t166) * t144;
t146 = t103 * pkin(3) + t147;
t145 = cos(qJ(5));
t143 = sin(qJ(4));
t129 = pkin(5) * t145 + pkin(4);
t111 = -t140 * t157 + t144 * t170 - t156 * t167;
t107 = t112 * t173 + t143 * t152;
t106 = t112 * t143 - t152 * t173;
t104 = -t120 * t156 + t121 * t144 - t134 * t155;
t102 = -t118 * t156 + t119 * t144 + t138 * t155;
t99 = t107 * t145 + t111 * t142;
t98 = -t107 * t142 + t111 * t145;
t97 = t105 * t173 + t143 * t153;
t96 = t105 * t143 - t153 * t173;
t95 = t103 * t173 + t143 * t154;
t94 = t103 * t143 - t154 * t173;
t93 = t104 * t142 + t145 * t97;
t92 = t104 * t145 - t97 * t142;
t91 = t102 * t142 + t145 * t95;
t90 = t102 * t145 - t95 * t142;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t138 - rSges(2,2) * t134) + g(2) * (rSges(2,1) * t134 + rSges(2,2) * t138) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t121 + rSges(3,2) * t120 + t161) + g(2) * (rSges(3,1) * t119 + rSges(3,2) * t118 + t131) + g(3) * (rSges(3,3) * t140 + t162) + (g(1) * rSges(3,3) * t134 + g(3) * (rSges(3,1) * t133 + rSges(3,2) * t137) + g(2) * (-rSges(3,3) - qJ(2)) * t138) * t136) - m(4) * (g(1) * (t105 * rSges(4,1) - t104 * rSges(4,2) + rSges(4,3) * t153 + t151) + g(2) * (t103 * rSges(4,1) - t102 * rSges(4,2) + rSges(4,3) * t154 + t147) + g(3) * (t112 * rSges(4,1) - t111 * rSges(4,2) + rSges(4,3) * t152 + t150)) - m(5) * (g(1) * (rSges(5,1) * t97 - rSges(5,2) * t96 + t104 * t176 + t149) + g(2) * (rSges(5,1) * t95 - rSges(5,2) * t94 + t102 * t176 + t146) + g(3) * (rSges(5,1) * t107 - rSges(5,2) * t106 + t111 * t176 + t148)) - m(6) * (g(1) * (rSges(6,1) * t93 + rSges(6,2) * t92 + pkin(4) * t97 + pkin(9) * t104 + t175 * t96 + t149) + g(2) * (rSges(6,1) * t91 + rSges(6,2) * t90 + pkin(4) * t95 + t102 * pkin(9) + t175 * t94 + t146) + g(3) * (rSges(6,1) * t99 + rSges(6,2) * t98 + pkin(4) * t107 + pkin(9) * t111 + t106 * t175 + t148)) - m(7) * (g(1) * (rSges(7,1) * t93 + rSges(7,2) * t92 + t104 * t158 + t129 * t97 + t172 * t96 + t149) + g(2) * (rSges(7,1) * t91 + rSges(7,2) * t90 + t102 * t158 + t129 * t95 + t172 * t94 + t146) + g(3) * (rSges(7,1) * t99 + rSges(7,2) * t98 + t106 * t172 + t107 * t129 + t111 * t158 + t148));
U  = t1;
