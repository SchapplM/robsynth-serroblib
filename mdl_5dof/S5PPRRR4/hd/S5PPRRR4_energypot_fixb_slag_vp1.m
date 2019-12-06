% Calculate potential energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:31
% EndTime: 2019-12-05 15:18:32
% DurationCPUTime: 0.44s
% Computational Cost: add. (356->108), mult. (889->159), div. (0->0), fcn. (1108->14), ass. (0->55)
t121 = sin(pkin(6));
t125 = cos(pkin(6));
t126 = cos(pkin(5));
t122 = sin(pkin(5));
t123 = cos(pkin(11));
t149 = t122 * t123;
t104 = -t121 * t149 + t125 * t126;
t119 = sin(pkin(11));
t124 = cos(pkin(10));
t120 = sin(pkin(10));
t150 = t120 * t126;
t107 = -t119 * t124 - t123 * t150;
t147 = t122 * t125;
t98 = -t107 * t121 + t120 * t147;
t156 = rSges(5,3) + pkin(8);
t155 = pkin(9) + rSges(6,3);
t154 = cos(qJ(3));
t152 = t119 * t122;
t151 = t120 * t122;
t148 = t122 * t124;
t146 = t124 * t126;
t144 = t126 * qJ(2) + qJ(1);
t143 = t124 * pkin(1) + qJ(2) * t151;
t140 = t121 * t154;
t139 = t125 * t154;
t138 = t122 * t140;
t105 = -t119 * t120 + t123 * t146;
t97 = -t105 * t121 - t124 * t147;
t108 = -t119 * t150 + t123 * t124;
t137 = t108 * pkin(2) + t98 * pkin(7) + t143;
t129 = sin(qJ(3));
t89 = t108 * t154 + (t107 * t125 + t121 * t151) * t129;
t136 = t89 * pkin(3) + t137;
t135 = pkin(2) * t152 + t104 * pkin(7) + t144;
t96 = t126 * t121 * t129 + (t123 * t125 * t129 + t119 * t154) * t122;
t134 = t96 * pkin(3) + t135;
t106 = t119 * t146 + t120 * t123;
t117 = t120 * pkin(1);
t133 = t106 * pkin(2) + pkin(7) * t97 - qJ(2) * t148 + t117;
t87 = t106 * t154 + (t105 * t125 - t121 * t148) * t129;
t132 = t87 * pkin(3) + t133;
t131 = cos(qJ(4));
t130 = cos(qJ(5));
t128 = sin(qJ(4));
t127 = sin(qJ(5));
t95 = -t126 * t140 + t129 * t152 - t139 * t149;
t91 = t104 * t128 + t131 * t96;
t90 = -t104 * t131 + t128 * t96;
t88 = -t107 * t139 + t108 * t129 - t120 * t138;
t86 = -t105 * t139 + t106 * t129 + t124 * t138;
t83 = t98 * t128 + t131 * t89;
t82 = t128 * t89 - t98 * t131;
t81 = t97 * t128 + t131 * t87;
t80 = t128 * t87 - t97 * t131;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t124 - rSges(2,2) * t120) + g(2) * (rSges(2,1) * t120 + rSges(2,2) * t124) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t108 + rSges(3,2) * t107 + t143) + g(2) * (rSges(3,1) * t106 + rSges(3,2) * t105 + t117) + g(3) * (rSges(3,3) * t126 + t144) + (g(1) * rSges(3,3) * t120 + g(3) * (rSges(3,1) * t119 + rSges(3,2) * t123) + g(2) * (-rSges(3,3) - qJ(2)) * t124) * t122) - m(4) * (g(1) * (rSges(4,1) * t89 - rSges(4,2) * t88 + rSges(4,3) * t98 + t137) + g(2) * (rSges(4,1) * t87 - rSges(4,2) * t86 + rSges(4,3) * t97 + t133) + g(3) * (rSges(4,1) * t96 - rSges(4,2) * t95 + rSges(4,3) * t104 + t135)) - m(5) * (g(1) * (rSges(5,1) * t83 - rSges(5,2) * t82 + t156 * t88 + t136) + g(2) * (rSges(5,1) * t81 - rSges(5,2) * t80 + t156 * t86 + t132) + g(3) * (rSges(5,1) * t91 - rSges(5,2) * t90 + t156 * t95 + t134)) - m(6) * (g(1) * (t83 * pkin(4) + t88 * pkin(8) + (t127 * t88 + t130 * t83) * rSges(6,1) + (-t83 * t127 + t130 * t88) * rSges(6,2) + t155 * t82 + t136) + g(2) * (t81 * pkin(4) + t86 * pkin(8) + (t127 * t86 + t81 * t130) * rSges(6,1) + (-t127 * t81 + t130 * t86) * rSges(6,2) + t155 * t80 + t132) + g(3) * (t91 * pkin(4) + t95 * pkin(8) + (t127 * t95 + t130 * t91) * rSges(6,1) + (-t127 * t91 + t130 * t95) * rSges(6,2) + t155 * t90 + t134));
U = t1;
