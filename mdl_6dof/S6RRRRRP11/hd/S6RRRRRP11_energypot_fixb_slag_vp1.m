% Calculate potential energy for
% S6RRRRRP11
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:29
% EndTime: 2019-03-10 02:32:29
% DurationCPUTime: 0.50s
% Computational Cost: add. (551->129), mult. (1388->177), div. (0->0), fcn. (1753->14), ass. (0->65)
t137 = cos(pkin(6));
t143 = sin(qJ(1));
t145 = cos(qJ(2));
t164 = t143 * t145;
t142 = sin(qJ(2));
t146 = cos(qJ(1));
t166 = t142 * t146;
t121 = -t137 * t164 - t166;
t134 = sin(pkin(7));
t136 = cos(pkin(7));
t135 = sin(pkin(6));
t170 = t135 * t143;
t154 = -t121 * t134 + t136 * t170;
t169 = t135 * t145;
t153 = -t134 * t169 + t136 * t137;
t178 = rSges(5,3) + pkin(11);
t177 = rSges(6,3) + pkin(12);
t176 = cos(qJ(3));
t175 = cos(qJ(4));
t174 = t137 * pkin(9) + pkin(8);
t173 = rSges(7,3) + qJ(6) + pkin(12);
t171 = t135 * t142;
t168 = t135 * t146;
t165 = t143 * t142;
t163 = t145 * t146;
t162 = t146 * pkin(1) + pkin(9) * t170;
t139 = sin(qJ(5));
t159 = pkin(5) * t139 + pkin(11);
t158 = t134 * t176;
t157 = t136 * t176;
t156 = t135 * t158;
t119 = t137 * t163 - t165;
t155 = -t119 * t134 - t136 * t168;
t122 = -t137 * t165 + t163;
t152 = t122 * pkin(2) + t154 * pkin(10) + t162;
t151 = pkin(2) * t171 + t153 * pkin(10) + t174;
t141 = sin(qJ(3));
t108 = t122 * t176 + (t121 * t136 + t134 * t170) * t141;
t150 = t108 * pkin(3) + t152;
t113 = t137 * t134 * t141 + (t136 * t141 * t145 + t176 * t142) * t135;
t149 = t113 * pkin(3) + t151;
t120 = t137 * t166 + t164;
t132 = t143 * pkin(1);
t148 = t120 * pkin(2) - pkin(9) * t168 + t155 * pkin(10) + t132;
t106 = t120 * t176 + (t119 * t136 - t134 * t168) * t141;
t147 = t106 * pkin(3) + t148;
t144 = cos(qJ(5));
t140 = sin(qJ(4));
t130 = pkin(5) * t144 + pkin(4);
t112 = -t137 * t158 + t141 * t171 - t157 * t169;
t107 = -t121 * t157 + t122 * t141 - t143 * t156;
t105 = -t119 * t157 + t120 * t141 + t146 * t156;
t104 = t113 * t175 + t153 * t140;
t103 = t113 * t140 - t153 * t175;
t100 = t108 * t175 + t154 * t140;
t99 = t108 * t140 - t154 * t175;
t98 = t106 * t175 + t155 * t140;
t97 = t106 * t140 - t155 * t175;
t96 = t104 * t144 + t112 * t139;
t95 = -t104 * t139 + t112 * t144;
t94 = t100 * t144 + t107 * t139;
t93 = -t100 * t139 + t107 * t144;
t92 = t105 * t139 + t144 * t98;
t91 = t105 * t144 - t98 * t139;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t146 - t143 * rSges(2,2)) + g(2) * (t143 * rSges(2,1) + rSges(2,2) * t146) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t122 + rSges(3,2) * t121 + t162) + g(2) * (t120 * rSges(3,1) + t119 * rSges(3,2) + t132) + g(3) * (rSges(3,3) * t137 + t174) + (g(1) * rSges(3,3) * t143 + g(3) * (rSges(3,1) * t142 + rSges(3,2) * t145) + g(2) * (-rSges(3,3) - pkin(9)) * t146) * t135) - m(4) * (g(1) * (t108 * rSges(4,1) - t107 * rSges(4,2) + t154 * rSges(4,3) + t152) + g(2) * (t106 * rSges(4,1) - t105 * rSges(4,2) + t155 * rSges(4,3) + t148) + g(3) * (t113 * rSges(4,1) - t112 * rSges(4,2) + t153 * rSges(4,3) + t151)) - m(5) * (g(1) * (rSges(5,1) * t100 - rSges(5,2) * t99 + t178 * t107 + t150) + g(2) * (t98 * rSges(5,1) - t97 * rSges(5,2) + t178 * t105 + t147) + g(3) * (rSges(5,1) * t104 - rSges(5,2) * t103 + t178 * t112 + t149)) - m(6) * (g(1) * (rSges(6,1) * t94 + rSges(6,2) * t93 + pkin(4) * t100 + pkin(11) * t107 + t177 * t99 + t150) + g(2) * (t92 * rSges(6,1) + t91 * rSges(6,2) + t98 * pkin(4) + t105 * pkin(11) + t177 * t97 + t147) + g(3) * (rSges(6,1) * t96 + rSges(6,2) * t95 + pkin(4) * t104 + pkin(11) * t112 + t177 * t103 + t149)) - m(7) * (g(1) * (rSges(7,1) * t94 + rSges(7,2) * t93 + t100 * t130 + t159 * t107 + t173 * t99 + t150) + g(2) * (t92 * rSges(7,1) + t91 * rSges(7,2) + t159 * t105 + t98 * t130 + t173 * t97 + t147) + g(3) * (rSges(7,1) * t96 + rSges(7,2) * t95 + t173 * t103 + t104 * t130 + t159 * t112 + t149));
U  = t1;
