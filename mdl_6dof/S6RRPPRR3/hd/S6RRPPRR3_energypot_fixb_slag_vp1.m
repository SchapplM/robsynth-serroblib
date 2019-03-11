% Calculate potential energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:24
% EndTime: 2019-03-09 08:54:25
% DurationCPUTime: 0.60s
% Computational Cost: add. (383->126), mult. (724->172), div. (0->0), fcn. (880->14), ass. (0->56)
t119 = sin(pkin(11));
t122 = cos(pkin(11));
t126 = sin(qJ(2));
t129 = cos(qJ(2));
t102 = -t126 * t119 + t129 * t122;
t153 = pkin(10) + rSges(7,3);
t152 = pkin(2) * t126;
t123 = cos(pkin(6));
t151 = t123 * pkin(8) + pkin(7);
t150 = qJ(4) + rSges(5,3);
t120 = sin(pkin(6));
t100 = t123 * t152 + (-pkin(8) - qJ(3)) * t120;
t112 = t129 * pkin(2) + pkin(1);
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t149 = t130 * t100 + t127 * t112;
t118 = sin(pkin(12));
t148 = t123 * t118;
t146 = t127 * t120;
t145 = t127 * t126;
t144 = t127 * t129;
t142 = t130 * t120;
t141 = t130 * t126;
t140 = t130 * t129;
t139 = t118 * t146;
t138 = t118 * t142;
t107 = t130 * t112;
t137 = -t127 * t100 + t107;
t136 = t123 * qJ(3) + t120 * t152 + t151;
t135 = t129 * t119 + t126 * t122;
t121 = cos(pkin(12));
t111 = t121 * pkin(4) + pkin(3);
t124 = -pkin(9) - qJ(4);
t132 = t102 * t123;
t89 = -t127 * t132 - t130 * t135;
t99 = t135 * t123;
t90 = t130 * t102 - t127 * t99;
t134 = pkin(4) * t139 + t90 * t111 + t89 * t124 + t137;
t97 = t102 * t120;
t98 = t135 * t120;
t133 = pkin(4) * t148 + t98 * t111 + t97 * t124 + t136;
t87 = -t127 * t135 + t130 * t132;
t88 = t127 * t102 + t130 * t99;
t131 = -pkin(4) * t138 + t88 * t111 + t87 * t124 + t149;
t128 = cos(qJ(6));
t125 = sin(qJ(6));
t117 = pkin(12) + qJ(5);
t114 = cos(t117);
t113 = sin(t117);
t92 = t123 * t113 + t98 * t114;
t91 = t98 * t113 - t123 * t114;
t82 = t113 * t146 + t90 * t114;
t81 = t90 * t113 - t114 * t146;
t80 = -t113 * t142 + t88 * t114;
t79 = t88 * t113 + t114 * t142;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t130 * rSges(2,1) - t127 * rSges(2,2)) + g(2) * (t127 * rSges(2,1) + t130 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t130 * pkin(1) + (-t123 * t145 + t140) * rSges(3,1) + (-t123 * t144 - t141) * rSges(3,2)) + g(2) * (t127 * pkin(1) + (t123 * t141 + t144) * rSges(3,1) + (t123 * t140 - t145) * rSges(3,2)) + g(3) * (t123 * rSges(3,3) + t151) + (g(3) * (rSges(3,1) * t126 + rSges(3,2) * t129) + (g(1) * t127 - g(2) * t130) * (rSges(3,3) + pkin(8))) * t120) - m(4) * (g(1) * (t90 * rSges(4,1) + t89 * rSges(4,2) + t107 + (rSges(4,3) * t120 - t100) * t127) + g(2) * (t88 * rSges(4,1) + t87 * rSges(4,2) - rSges(4,3) * t142 + t149) + g(3) * (t98 * rSges(4,1) + t97 * rSges(4,2) + t123 * rSges(4,3) + t136)) - m(5) * (g(1) * (t90 * pkin(3) + (t90 * t121 + t139) * rSges(5,1) + (-t90 * t118 + t121 * t146) * rSges(5,2) - t150 * t89 + t137) + g(2) * (t88 * pkin(3) + (t88 * t121 - t138) * rSges(5,1) + (-t88 * t118 - t121 * t142) * rSges(5,2) - t150 * t87 + t149) + g(3) * (t98 * pkin(3) + (t98 * t121 + t148) * rSges(5,1) + (-t98 * t118 + t123 * t121) * rSges(5,2) - t150 * t97 + t136)) - m(6) * (g(1) * (t82 * rSges(6,1) - t81 * rSges(6,2) - t89 * rSges(6,3) + t134) + g(2) * (t80 * rSges(6,1) - t79 * rSges(6,2) - t87 * rSges(6,3) + t131) + g(3) * (t92 * rSges(6,1) - t91 * rSges(6,2) - t97 * rSges(6,3) + t133)) - m(7) * (g(1) * (t82 * pkin(5) + (-t89 * t125 + t82 * t128) * rSges(7,1) + (-t82 * t125 - t89 * t128) * rSges(7,2) + t153 * t81 + t134) + g(2) * (t80 * pkin(5) + (-t87 * t125 + t80 * t128) * rSges(7,1) + (-t80 * t125 - t87 * t128) * rSges(7,2) + t153 * t79 + t131) + g(3) * (t92 * pkin(5) + (-t97 * t125 + t92 * t128) * rSges(7,1) + (-t92 * t125 - t97 * t128) * rSges(7,2) + t153 * t91 + t133));
U  = t1;
