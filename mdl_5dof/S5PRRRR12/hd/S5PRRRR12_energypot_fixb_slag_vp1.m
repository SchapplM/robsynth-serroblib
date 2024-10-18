% Calculate potential energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:09
% EndTime: 2024-09-28 18:07:10
% DurationCPUTime: 0.53s
% Computational Cost: add. (316->142), mult. (398->208), div. (0->0), fcn. (399->26), ass. (0->71)
t147 = pkin(7) + pkin(8);
t131 = pkin(9) + t147;
t179 = t131 + rSges(5,3);
t178 = t147 + rSges(4,3);
t133 = sin(pkin(11));
t136 = cos(pkin(11));
t177 = g(1) * t133 - g(2) * t136;
t176 = rSges(3,3) + pkin(7);
t137 = cos(pkin(6));
t175 = t131 + (rSges(6,3) + pkin(10)) * t137;
t134 = sin(pkin(6));
t172 = pkin(10) * t134;
t146 = cos(qJ(2));
t122 = t146 * pkin(2) + pkin(1);
t139 = sin(qJ(5));
t168 = t137 * t139;
t143 = cos(qJ(5));
t167 = t137 * t143;
t138 = cos(pkin(5));
t166 = t138 * t133;
t165 = t138 * t136;
t164 = t138 * t139;
t142 = sin(qJ(2));
t163 = t138 * t142;
t162 = t138 * t143;
t161 = t138 * t146;
t132 = qJ(2) + qJ(3);
t160 = t133 * t172;
t159 = t136 * t172;
t158 = t137 * t164;
t157 = t137 * t162;
t124 = pkin(5) - t132;
t123 = pkin(5) + t132;
t155 = rSges(6,1) * t139 + rSges(6,2) * t143;
t141 = sin(qJ(3));
t145 = cos(qJ(3));
t154 = pkin(3) * t141 * t146 + (pkin(3) * t145 + pkin(2)) * t142;
t126 = cos(t132);
t153 = rSges(4,1) * t126 - rSges(4,2) * sin(t132) + t122;
t129 = qJ(4) + t132;
t119 = sin(t129);
t120 = cos(t129);
t152 = rSges(5,1) * t120 - rSges(5,2) * t119 + pkin(3) * t126 + t122;
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t97 = pkin(4) * t136 + t138 * t160;
t99 = pkin(4) * t166 - t159;
t151 = pkin(3) * t136 - t140 * t99 + t144 * t97;
t100 = pkin(4) * t165 + t160;
t98 = -pkin(4) * t133 + t138 * t159;
t150 = pkin(3) * t133 + t100 * t140 - t144 * t98;
t117 = qJ(4) + t123;
t107 = sin(t117) / 0.2e1;
t118 = -qJ(4) + t124;
t108 = cos(t118) / 0.2e1;
t113 = sin(t118);
t114 = cos(t117);
t135 = sin(pkin(5));
t149 = rSges(5,1) * (t107 - t113 / 0.2e1) + rSges(5,2) * (t108 + t114 / 0.2e1) + t138 * t154 - t179 * t135;
t110 = sin(t123) / 0.2e1;
t111 = cos(t124) / 0.2e1;
t115 = sin(t124);
t116 = cos(t123);
t148 = rSges(4,1) * (t110 - t115 / 0.2e1) + rSges(4,2) * (t111 + t116 / 0.2e1) + pkin(2) * t163 - t178 * t135;
t128 = t136 * pkin(1);
t127 = t133 * pkin(1);
t105 = -pkin(4) * t140 + t144 * t172;
t104 = pkin(4) * t144 + t140 * t172 + pkin(3);
t93 = pkin(3) * t165 + t100 * t144 + t140 * t98;
t92 = -pkin(3) * t166 - t140 * t97 - t144 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + rSges(1,2) * g(2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t136 - rSges(2,2) * t133) + g(2) * (rSges(2,1) * t133 + rSges(2,2) * t136) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t128 + (-t133 * t163 + t136 * t146) * rSges(3,1) + (-t133 * t161 - t136 * t142) * rSges(3,2)) + g(2) * (t127 + (t133 * t146 + t136 * t163) * rSges(3,1) + (-t133 * t142 + t136 * t161) * rSges(3,2)) + g(3) * (t176 * t138 + qJ(1)) + (g(3) * (rSges(3,1) * t142 + rSges(3,2) * t146) + t177 * t176) * t135) - m(4) * (g(1) * (-t133 * t148 + t136 * t153) + g(2) * (t133 * t153 + t136 * t148) + g(3) * (t135 * t142 * pkin(2) + qJ(1) + (t111 - t116 / 0.2e1) * rSges(4,1) + (t110 + t115 / 0.2e1) * rSges(4,2) + t178 * t138)) - m(5) * (g(1) * (-t133 * t149 + t136 * t152) + g(2) * (t133 * t152 + t136 * t149) + g(3) * (qJ(1) + (t108 - t114 / 0.2e1) * rSges(5,1) + (t107 + t113 / 0.2e1) * rSges(5,2) + t179 * t138 + t154 * t135)) - m(6) * (g(1) * ((t136 * pkin(2) + t92 * t141 + t145 * t151) * t146 + (-pkin(2) * t166 - t141 * t151 + t92 * t145) * t142 + t128 + ((-t133 * t158 + t136 * t143) * t120 + (-t133 * t162 - t136 * t168) * t119) * rSges(6,1) + ((-t133 * t157 - t136 * t139) * t120 + (t133 * t164 - t136 * t167) * t119) * rSges(6,2)) + g(2) * ((t133 * pkin(2) + t93 * t141 + t145 * t150) * t146 + (pkin(2) * t165 - t141 * t150 + t93 * t145) * t142 + t127 + ((t133 * t143 + t136 * t158) * t120 + (-t133 * t168 + t136 * t162) * t119) * rSges(6,1) + ((-t133 * t139 + t136 * t157) * t120 + (-t133 * t167 - t136 * t164) * t119) * rSges(6,2)) + g(3) * (t175 * t138 + qJ(1)) + (g(3) * t155 * t138 + (g(1) * (t119 * t136 + t120 * t166) + g(2) * (t119 * t133 - t120 * t165)) * rSges(6,3)) * t134 + (g(3) * ((t104 * t145 + t105 * t141 + pkin(2)) * t142 - (-t104 * t141 + t105 * t145) * t146 + (t119 * t143 + t120 * t168) * rSges(6,1) + (-t119 * t139 + t120 * t167) * rSges(6,2) - t120 * t134 * rSges(6,3)) + t177 * (t134 * t155 + t175)) * t135);
U = t1;
