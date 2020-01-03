% Calculate potential energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR14_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:48
% EndTime: 2019-12-31 19:16:48
% DurationCPUTime: 0.43s
% Computational Cost: add. (356->108), mult. (889->157), div. (0->0), fcn. (1108->14), ass. (0->55)
t121 = sin(pkin(11));
t126 = cos(pkin(5));
t133 = cos(qJ(1));
t124 = cos(pkin(11));
t130 = sin(qJ(1));
t146 = t130 * t124;
t109 = -t121 * t133 - t126 * t146;
t122 = sin(pkin(6));
t125 = cos(pkin(6));
t123 = sin(pkin(5));
t151 = t123 * t130;
t100 = -t109 * t122 + t125 * t151;
t152 = t123 * t124;
t106 = -t122 * t152 + t125 * t126;
t158 = rSges(5,3) + pkin(9);
t157 = pkin(10) + rSges(6,3);
t156 = cos(qJ(3));
t155 = t126 * qJ(2) + pkin(7);
t153 = t121 * t123;
t150 = t123 * t133;
t148 = t126 * t133;
t147 = t130 * t121;
t145 = t133 * pkin(1) + qJ(2) * t151;
t142 = t122 * t156;
t141 = t125 * t156;
t140 = t123 * t142;
t107 = t124 * t148 - t147;
t99 = -t107 * t122 - t125 * t150;
t110 = t124 * t133 - t126 * t147;
t139 = t110 * pkin(2) + t100 * pkin(8) + t145;
t138 = pkin(2) * t153 + t106 * pkin(8) + t155;
t129 = sin(qJ(3));
t93 = t110 * t156 + (t109 * t125 + t122 * t151) * t129;
t137 = t93 * pkin(3) + t139;
t98 = t126 * t122 * t129 + (t124 * t125 * t129 + t156 * t121) * t123;
t136 = t98 * pkin(3) + t138;
t108 = t121 * t148 + t146;
t119 = t130 * pkin(1);
t135 = t108 * pkin(2) + t99 * pkin(8) - qJ(2) * t150 + t119;
t91 = t108 * t156 + (t107 * t125 - t122 * t150) * t129;
t134 = t91 * pkin(3) + t135;
t132 = cos(qJ(4));
t131 = cos(qJ(5));
t128 = sin(qJ(4));
t127 = sin(qJ(5));
t97 = -t126 * t142 + t129 * t153 - t141 * t152;
t92 = -t109 * t141 + t110 * t129 - t130 * t140;
t90 = -t107 * t141 + t108 * t129 + t133 * t140;
t89 = t106 * t128 + t132 * t98;
t88 = -t106 * t132 + t128 * t98;
t85 = t100 * t128 + t132 * t93;
t84 = -t100 * t132 + t128 * t93;
t83 = t128 * t99 + t132 * t91;
t82 = t128 * t91 - t99 * t132;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t133 - t130 * rSges(2,2)) + g(2) * (t130 * rSges(2,1) + rSges(2,2) * t133) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t110 + rSges(3,2) * t109 + t145) + g(2) * (t108 * rSges(3,1) + t107 * rSges(3,2) + t119) + g(3) * (rSges(3,3) * t126 + t155) + (g(1) * rSges(3,3) * t130 + g(3) * (rSges(3,1) * t121 + rSges(3,2) * t124) + g(2) * (-rSges(3,3) - qJ(2)) * t133) * t123) - m(4) * (g(1) * (rSges(4,1) * t93 - rSges(4,2) * t92 + rSges(4,3) * t100 + t139) + g(2) * (t91 * rSges(4,1) - t90 * rSges(4,2) + t99 * rSges(4,3) + t135) + g(3) * (rSges(4,1) * t98 - rSges(4,2) * t97 + rSges(4,3) * t106 + t138)) - m(5) * (g(1) * (rSges(5,1) * t85 - rSges(5,2) * t84 + t158 * t92 + t137) + g(2) * (t83 * rSges(5,1) - t82 * rSges(5,2) + t158 * t90 + t134) + g(3) * (rSges(5,1) * t89 - rSges(5,2) * t88 + t158 * t97 + t136)) - m(6) * (g(1) * (t85 * pkin(4) + t92 * pkin(9) + (t127 * t92 + t131 * t85) * rSges(6,1) + (-t127 * t85 + t131 * t92) * rSges(6,2) + t157 * t84 + t137) + g(2) * (t83 * pkin(4) + t90 * pkin(9) + (t127 * t90 + t131 * t83) * rSges(6,1) + (-t127 * t83 + t131 * t90) * rSges(6,2) + t157 * t82 + t134) + g(3) * (t89 * pkin(4) + t97 * pkin(9) + (t127 * t97 + t131 * t89) * rSges(6,1) + (-t127 * t89 + t131 * t97) * rSges(6,2) + t157 * t88 + t136));
U = t1;
