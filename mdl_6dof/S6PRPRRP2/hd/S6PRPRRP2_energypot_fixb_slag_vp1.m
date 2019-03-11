% Calculate potential energy for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:50
% EndTime: 2019-03-08 19:59:51
% DurationCPUTime: 0.49s
% Computational Cost: add. (418->114), mult. (938->155), div. (0->0), fcn. (1177->12), ass. (0->56)
t122 = sin(pkin(11));
t125 = cos(pkin(11));
t130 = sin(qJ(2));
t133 = cos(qJ(2));
t157 = t130 * t122 - t125 * t133;
t156 = rSges(7,1) + pkin(5);
t155 = rSges(7,2) + pkin(9);
t153 = rSges(5,3) + pkin(8);
t152 = rSges(6,3) + pkin(9);
t151 = rSges(7,3) + qJ(6);
t123 = sin(pkin(10));
t124 = sin(pkin(6));
t150 = t123 * t124;
t126 = cos(pkin(10));
t148 = t126 * t124;
t127 = cos(pkin(6));
t147 = t127 * t130;
t146 = t127 * t133;
t144 = t127 * pkin(7) + qJ(1);
t110 = pkin(2) * t147 + (-pkin(7) - qJ(3)) * t124;
t119 = pkin(2) * t133 + pkin(1);
t143 = t126 * t110 + t123 * t119;
t138 = t122 * t133 + t130 * t125;
t109 = t138 * t127;
t97 = t109 * t126 - t123 * t157;
t142 = t97 * pkin(3) + t143;
t141 = t124 * t130 * pkin(2) + t127 * qJ(3) + t144;
t114 = t126 * t119;
t99 = -t109 * t123 - t126 * t157;
t140 = t99 * pkin(3) - t110 * t123 + t114;
t108 = t138 * t124;
t139 = t108 * pkin(3) + t141;
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t89 = -t129 * t148 + t132 * t97;
t136 = t157 * t127;
t96 = -t123 * t138 - t126 * t136;
t137 = t89 * pkin(4) - pkin(8) * t96 + t142;
t91 = t129 * t150 + t132 * t99;
t98 = t123 * t136 - t126 * t138;
t135 = t91 * pkin(4) - pkin(8) * t98 + t140;
t102 = t108 * t132 + t127 * t129;
t107 = t157 * t124;
t134 = t102 * pkin(4) + pkin(8) * t107 + t139;
t131 = cos(qJ(5));
t128 = sin(qJ(5));
t101 = t108 * t129 - t127 * t132;
t90 = t129 * t99 - t132 * t150;
t88 = t129 * t97 + t132 * t148;
t85 = t102 * t131 + t107 * t128;
t84 = t102 * t128 - t107 * t131;
t83 = -t128 * t98 + t131 * t91;
t82 = t128 * t91 + t98 * t131;
t81 = -t128 * t96 + t131 * t89;
t80 = t128 * t89 + t96 * t131;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t126 - rSges(2,2) * t123) + g(2) * (rSges(2,1) * t123 + rSges(2,2) * t126) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t126 * pkin(1) + (-t123 * t147 + t126 * t133) * rSges(3,1) + (-t123 * t146 - t126 * t130) * rSges(3,2)) + g(2) * (t123 * pkin(1) + (t123 * t133 + t126 * t147) * rSges(3,1) + (-t123 * t130 + t126 * t146) * rSges(3,2)) + g(3) * (t127 * rSges(3,3) + t144) + (g(3) * (rSges(3,1) * t130 + rSges(3,2) * t133) + (g(1) * t123 - g(2) * t126) * (rSges(3,3) + pkin(7))) * t124) - m(4) * (g(1) * (rSges(4,1) * t99 + rSges(4,2) * t98 + t114 + (rSges(4,3) * t124 - t110) * t123) + g(2) * (rSges(4,1) * t97 + rSges(4,2) * t96 - rSges(4,3) * t148 + t143) + g(3) * (rSges(4,1) * t108 - rSges(4,2) * t107 + rSges(4,3) * t127 + t141)) - m(5) * (g(1) * (rSges(5,1) * t91 - rSges(5,2) * t90 - t153 * t98 + t140) + g(2) * (rSges(5,1) * t89 - rSges(5,2) * t88 - t153 * t96 + t142) + g(3) * (rSges(5,1) * t102 - rSges(5,2) * t101 + t107 * t153 + t139)) - m(6) * (g(1) * (rSges(6,1) * t83 - rSges(6,2) * t82 + t152 * t90 + t135) + g(2) * (rSges(6,1) * t81 - rSges(6,2) * t80 + t152 * t88 + t137) + g(3) * (rSges(6,1) * t85 - rSges(6,2) * t84 + t101 * t152 + t134)) - m(7) * (g(1) * (t151 * t82 + t155 * t90 + t156 * t83 + t135) + g(2) * (t151 * t80 + t155 * t88 + t156 * t81 + t137) + g(3) * (t155 * t101 + t151 * t84 + t156 * t85 + t134));
U  = t1;
