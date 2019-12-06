% Calculate potential energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:10
% EndTime: 2019-12-05 17:23:10
% DurationCPUTime: 0.44s
% Computational Cost: add. (356->108), mult. (889->159), div. (0->0), fcn. (1108->14), ass. (0->55)
t120 = sin(pkin(6));
t123 = cos(pkin(6));
t124 = cos(pkin(5));
t121 = sin(pkin(5));
t131 = cos(qJ(2));
t148 = t121 * t131;
t104 = -t120 * t148 + t124 * t123;
t119 = sin(pkin(11));
t122 = cos(pkin(11));
t128 = sin(qJ(2));
t145 = t124 * t131;
t107 = -t119 * t145 - t122 * t128;
t150 = t121 * t123;
t98 = -t107 * t120 + t119 * t150;
t156 = rSges(5,3) + pkin(9);
t155 = pkin(10) + rSges(6,3);
t154 = cos(qJ(3));
t152 = t119 * t121;
t151 = t121 * t122;
t149 = t121 * t128;
t146 = t124 * t128;
t144 = t124 * pkin(7) + qJ(1);
t143 = t122 * pkin(1) + pkin(7) * t152;
t140 = t120 * t154;
t139 = t123 * t154;
t138 = t121 * t140;
t105 = -t119 * t128 + t122 * t145;
t97 = -t105 * t120 - t122 * t150;
t108 = -t119 * t146 + t122 * t131;
t137 = t108 * pkin(2) + t98 * pkin(8) + t143;
t127 = sin(qJ(3));
t89 = t108 * t154 + (t107 * t123 + t120 * t152) * t127;
t136 = t89 * pkin(3) + t137;
t135 = pkin(2) * t149 + t104 * pkin(8) + t144;
t96 = t124 * t120 * t127 + (t123 * t127 * t131 + t154 * t128) * t121;
t134 = t96 * pkin(3) + t135;
t106 = t119 * t131 + t122 * t146;
t116 = t119 * pkin(1);
t133 = t106 * pkin(2) - pkin(7) * t151 + t97 * pkin(8) + t116;
t87 = t106 * t154 + (t105 * t123 - t120 * t151) * t127;
t132 = t87 * pkin(3) + t133;
t130 = cos(qJ(4));
t129 = cos(qJ(5));
t126 = sin(qJ(4));
t125 = sin(qJ(5));
t95 = -t124 * t140 + t127 * t149 - t139 * t148;
t91 = t104 * t126 + t130 * t96;
t90 = -t104 * t130 + t126 * t96;
t88 = -t107 * t139 + t108 * t127 - t119 * t138;
t86 = -t105 * t139 + t106 * t127 + t122 * t138;
t83 = t126 * t98 + t130 * t89;
t82 = t126 * t89 - t98 * t130;
t81 = t126 * t97 + t130 * t87;
t80 = t126 * t87 - t97 * t130;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t122 - rSges(2,2) * t119) + g(2) * (rSges(2,1) * t119 + rSges(2,2) * t122) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t108 + rSges(3,2) * t107 + t143) + g(2) * (rSges(3,1) * t106 + rSges(3,2) * t105 + t116) + g(3) * (t124 * rSges(3,3) + t144) + (g(1) * rSges(3,3) * t119 + g(3) * (rSges(3,1) * t128 + rSges(3,2) * t131) + g(2) * (-rSges(3,3) - pkin(7)) * t122) * t121) - m(4) * (g(1) * (rSges(4,1) * t89 - rSges(4,2) * t88 + rSges(4,3) * t98 + t137) + g(2) * (rSges(4,1) * t87 - rSges(4,2) * t86 + rSges(4,3) * t97 + t133) + g(3) * (t96 * rSges(4,1) - t95 * rSges(4,2) + t104 * rSges(4,3) + t135)) - m(5) * (g(1) * (rSges(5,1) * t83 - rSges(5,2) * t82 + t156 * t88 + t136) + g(2) * (rSges(5,1) * t81 - rSges(5,2) * t80 + t156 * t86 + t132) + g(3) * (t91 * rSges(5,1) - t90 * rSges(5,2) + t156 * t95 + t134)) - m(6) * (g(1) * (t83 * pkin(4) + t88 * pkin(9) + (t125 * t88 + t129 * t83) * rSges(6,1) + (-t125 * t83 + t129 * t88) * rSges(6,2) + t155 * t82 + t136) + g(2) * (t81 * pkin(4) + t86 * pkin(9) + (t125 * t86 + t129 * t81) * rSges(6,1) + (-t125 * t81 + t129 * t86) * rSges(6,2) + t155 * t80 + t132) + g(3) * (t91 * pkin(4) + t95 * pkin(9) + (t125 * t95 + t129 * t91) * rSges(6,1) + (-t125 * t91 + t129 * t95) * rSges(6,2) + t155 * t90 + t134));
U = t1;
