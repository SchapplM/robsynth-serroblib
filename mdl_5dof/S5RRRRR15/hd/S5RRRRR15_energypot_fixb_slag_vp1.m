% Calculate potential energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR15_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR15_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:52
% EndTime: 2024-09-27 22:23:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (316->106), mult. (360->137), div. (0->0), fcn. (365->26), ass. (0->62)
t137 = pkin(8) + pkin(9);
t121 = pkin(10) + t137;
t123 = sin(pkin(6));
t125 = cos(pkin(6));
t127 = sin(qJ(5));
t132 = cos(qJ(5));
t171 = t121 + (rSges(6,1) * t127 + rSges(6,2) * t132) * t123 + (pkin(11) + rSges(6,3)) * t125;
t122 = qJ(2) + qJ(3);
t119 = qJ(4) + t122;
t112 = cos(t119);
t130 = sin(qJ(2));
t135 = cos(qJ(2));
t152 = t132 * t112;
t156 = t127 * t112;
t111 = sin(t119);
t157 = t111 * t132;
t158 = t111 * t127;
t160 = t123 * rSges(6,3);
t129 = sin(qJ(3));
t134 = cos(qJ(3));
t128 = sin(qJ(4));
t133 = cos(qJ(4));
t161 = pkin(11) * t123;
t95 = pkin(4) * t133 + t128 * t161 + pkin(3);
t98 = -pkin(4) * t128 + t133 * t161;
t89 = t129 * t98 + t134 * t95 + pkin(2);
t90 = -t129 * t95 + t134 * t98;
t170 = (t125 * t156 + t157) * rSges(6,1) + (t125 * t152 - t158) * rSges(6,2) - t112 * t160 + t89 * t130 - t90 * t135;
t169 = t121 + rSges(5,3);
t168 = t137 + rSges(4,3);
t124 = sin(pkin(5));
t126 = cos(pkin(5));
t167 = t171 * t124 - t170 * t126;
t166 = rSges(3,3) + pkin(8);
t162 = pkin(2) * t130;
t114 = t135 * pkin(2) + pkin(1);
t131 = sin(qJ(1));
t155 = t130 * t131;
t136 = cos(qJ(1));
t154 = t130 * t136;
t153 = t131 * t135;
t151 = t135 * t136;
t116 = pkin(5) - t122;
t115 = pkin(5) + t122;
t118 = cos(t122);
t147 = t112 * rSges(5,1) - t111 * rSges(5,2) + pkin(3) * t118 + t114;
t146 = t135 * t129 * pkin(3) + t130 * (pkin(3) * t134 + pkin(2));
t145 = t118 * rSges(4,1) - sin(t122) * rSges(4,2) + t114;
t109 = qJ(4) + t115;
t100 = sin(t109) / 0.2e1;
t110 = -qJ(4) + t116;
t101 = cos(t110) / 0.2e1;
t103 = sin(t110);
t104 = cos(t109);
t142 = (t100 - t103 / 0.2e1) * rSges(5,1) + (t101 + t104 / 0.2e1) * rSges(5,2) + t146 * t126 - t169 * t124;
t105 = sin(t115);
t106 = sin(t116);
t107 = cos(t115);
t108 = cos(t116);
t141 = t126 * t162 + (t105 - t106) * rSges(4,1) / 0.2e1 + (t108 + t107) * rSges(4,2) / 0.2e1 - t168 * t124;
t138 = t130 * t90 + t135 * t89 + pkin(1) + (-t125 * t158 + t152) * rSges(6,1) + (-t125 * t157 - t156) * rSges(6,2) + t111 * t160;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t136 - rSges(2,2) * t131) + g(2) * (rSges(2,1) * t131 + rSges(2,2) * t136) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t136 * pkin(1) + (-t126 * t155 + t151) * rSges(3,1) + (-t126 * t153 - t154) * rSges(3,2)) + g(2) * (t131 * pkin(1) + (t126 * t154 + t153) * rSges(3,1) + (t126 * t151 - t155) * rSges(3,2)) + g(3) * (t166 * t126 + pkin(7)) + (g(3) * (rSges(3,1) * t130 + rSges(3,2) * t135) + (g(1) * t131 - g(2) * t136) * t166) * t124) - m(4) * (g(1) * (-t141 * t131 + t145 * t136) + g(2) * (t145 * t131 + t141 * t136) + g(3) * (t124 * t162 + pkin(7) + (t108 / 0.2e1 - t107 / 0.2e1) * rSges(4,1) + (t105 / 0.2e1 + t106 / 0.2e1) * rSges(4,2) + t168 * t126)) - m(5) * (g(1) * (-t142 * t131 + t147 * t136) + g(2) * (t147 * t131 + t142 * t136) + g(3) * (pkin(7) + (t101 - t104 / 0.2e1) * rSges(5,1) + (t100 + t103 / 0.2e1) * rSges(5,2) + t169 * t126 + t146 * t124)) - m(6) * (g(1) * (t131 * t167 + t138 * t136) + g(2) * (t138 * t131 - t136 * t167) + g(3) * (t170 * t124 + t171 * t126 + pkin(7)));
U = t1;
