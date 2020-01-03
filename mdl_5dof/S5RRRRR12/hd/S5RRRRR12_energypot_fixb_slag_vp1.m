% Calculate potential energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR12_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:29
% EndTime: 2019-12-31 22:46:30
% DurationCPUTime: 0.42s
% Computational Cost: add. (356->108), mult. (889->156), div. (0->0), fcn. (1108->14), ass. (0->56)
t123 = cos(pkin(5));
t128 = sin(qJ(1));
t131 = cos(qJ(2));
t146 = t128 * t131;
t127 = sin(qJ(2));
t132 = cos(qJ(1));
t148 = t127 * t132;
t108 = -t123 * t146 - t148;
t120 = sin(pkin(6));
t122 = cos(pkin(6));
t121 = sin(pkin(5));
t152 = t121 * t128;
t99 = -t108 * t120 + t122 * t152;
t151 = t121 * t131;
t105 = -t120 * t151 + t122 * t123;
t158 = rSges(5,3) + pkin(10);
t157 = pkin(11) + rSges(6,3);
t156 = cos(qJ(3));
t155 = t123 * pkin(8) + pkin(7);
t153 = t121 * t127;
t150 = t121 * t132;
t147 = t128 * t127;
t145 = t131 * t132;
t144 = t132 * pkin(1) + pkin(8) * t152;
t141 = t120 * t156;
t140 = t122 * t156;
t139 = t121 * t141;
t106 = t123 * t145 - t147;
t98 = -t106 * t120 - t122 * t150;
t109 = -t123 * t147 + t145;
t138 = t109 * pkin(2) + t99 * pkin(9) + t144;
t137 = pkin(2) * t153 + t105 * pkin(9) + t155;
t126 = sin(qJ(3));
t92 = t109 * t156 + (t108 * t122 + t120 * t152) * t126;
t136 = t92 * pkin(3) + t138;
t97 = t123 * t120 * t126 + (t122 * t126 * t131 + t156 * t127) * t121;
t135 = t97 * pkin(3) + t137;
t107 = t123 * t148 + t146;
t118 = t128 * pkin(1);
t134 = t107 * pkin(2) - pkin(8) * t150 + t98 * pkin(9) + t118;
t90 = t107 * t156 + (t106 * t122 - t120 * t150) * t126;
t133 = t90 * pkin(3) + t134;
t130 = cos(qJ(4));
t129 = cos(qJ(5));
t125 = sin(qJ(4));
t124 = sin(qJ(5));
t96 = -t123 * t141 + t126 * t153 - t140 * t151;
t91 = -t108 * t140 + t109 * t126 - t128 * t139;
t89 = -t106 * t140 + t107 * t126 + t132 * t139;
t88 = t105 * t125 + t130 * t97;
t87 = -t105 * t130 + t125 * t97;
t84 = t125 * t99 + t130 * t92;
t83 = t125 * t92 - t99 * t130;
t82 = t125 * t98 + t130 * t90;
t81 = t125 * t90 - t98 * t130;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t132 - t128 * rSges(2,2)) + g(2) * (t128 * rSges(2,1) + rSges(2,2) * t132) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t109 + rSges(3,2) * t108 + t144) + g(2) * (t107 * rSges(3,1) + t106 * rSges(3,2) + t118) + g(3) * (rSges(3,3) * t123 + t155) + (g(1) * rSges(3,3) * t128 + g(3) * (rSges(3,1) * t127 + rSges(3,2) * t131) + g(2) * (-rSges(3,3) - pkin(8)) * t132) * t121) - m(4) * (g(1) * (rSges(4,1) * t92 - rSges(4,2) * t91 + rSges(4,3) * t99 + t138) + g(2) * (t90 * rSges(4,1) - t89 * rSges(4,2) + t98 * rSges(4,3) + t134) + g(3) * (rSges(4,1) * t97 - rSges(4,2) * t96 + rSges(4,3) * t105 + t137)) - m(5) * (g(1) * (rSges(5,1) * t84 - rSges(5,2) * t83 + t158 * t91 + t136) + g(2) * (t82 * rSges(5,1) - t81 * rSges(5,2) + t158 * t89 + t133) + g(3) * (t88 * rSges(5,1) - rSges(5,2) * t87 + t158 * t96 + t135)) - m(6) * (g(1) * (t84 * pkin(4) + t91 * pkin(10) + (t124 * t91 + t129 * t84) * rSges(6,1) + (-t124 * t84 + t129 * t91) * rSges(6,2) + t157 * t83 + t136) + g(2) * (t82 * pkin(4) + t89 * pkin(10) + (t124 * t89 + t129 * t82) * rSges(6,1) + (-t124 * t82 + t129 * t89) * rSges(6,2) + t157 * t81 + t133) + g(3) * (t88 * pkin(4) + t96 * pkin(10) + (t124 * t96 + t129 * t88) * rSges(6,1) + (-t124 * t88 + t129 * t96) * rSges(6,2) + t157 * t87 + t135));
U = t1;
