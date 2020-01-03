% Calculate joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:22
% EndTime: 2019-12-31 19:23:26
% DurationCPUTime: 1.64s
% Computational Cost: add. (1938->213), mult. (5169->315), div. (0->0), fcn. (5979->8), ass. (0->102)
t107 = sin(qJ(2));
t109 = cos(qJ(2));
t185 = Icges(3,5) * t107 + Icges(3,6) * t109;
t110 = cos(qJ(1));
t146 = t107 * t110;
t145 = t109 * t110;
t182 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t181 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t180 = -Icges(4,4) - Icges(5,6) + Icges(6,6);
t179 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t178 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t177 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t176 = rSges(6,1) + pkin(4);
t154 = rSges(6,3) + qJ(5);
t106 = sin(pkin(5));
t108 = sin(qJ(1));
t153 = cos(pkin(5));
t83 = t106 * t146 + t108 * t153;
t143 = m(5) / 0.2e1 + m(6) / 0.2e1;
t175 = 0.2e1 * t143;
t152 = cos(pkin(8));
t127 = t153 * t152;
t151 = sin(pkin(8));
t113 = t107 * t127 + t109 * t151;
t132 = t110 * t152;
t61 = t106 * t132 + t108 * t113;
t126 = t153 * t151;
t118 = t107 * t126;
t134 = t108 * t152;
t136 = t106 * t151;
t62 = -t108 * t118 + t109 * t134 - t110 * t136;
t133 = t110 * t153;
t148 = t106 * t107;
t82 = t108 * t148 - t133;
t174 = t177 * t61 + t178 * t82 + t180 * t62;
t63 = -t106 * t134 + t110 * t113;
t64 = t108 * t136 + t109 * t132 - t110 * t118;
t173 = t177 * t63 + t178 * t83 + t180 * t64;
t172 = -t179 * t82 + t180 * t61 + t182 * t62;
t171 = -t179 * t83 + t180 * t63 + t182 * t64;
t170 = t178 * t61 - t179 * t62 + t181 * t82;
t169 = t178 * t63 - t179 * t64 + t181 * t83;
t147 = t106 * t109;
t80 = t107 * t151 - t109 * t127;
t81 = t107 * t152 + t109 * t126;
t168 = t179 * t147 + t180 * t80 + t182 * t81;
t167 = -t178 * t147 + t177 * t80 + t180 * t81;
t166 = -t181 * t147 + t178 * t80 - t179 * t81;
t165 = -t61 * rSges(6,2) - t176 * t82;
t164 = t108 ^ 2;
t163 = t110 ^ 2;
t84 = pkin(2) * t107 - qJ(3) * t147;
t158 = -rSges(4,1) * t81 + rSges(4,2) * t80 + rSges(4,3) * t147 - t84;
t157 = -pkin(3) * t81 - qJ(4) * t80 - t84;
t117 = pkin(2) * t109 + qJ(3) * t148;
t139 = pkin(2) * t145 + t83 * qJ(3);
t99 = qJ(3) * t133;
t156 = t108 * (t108 * t117 - t99) + t110 * t139;
t155 = t110 * rSges(3,3);
t144 = t110 * pkin(1) + t108 * pkin(7);
t142 = rSges(5,1) * t147 + rSges(5,2) * t81 - rSges(5,3) * t80 + t157;
t141 = t64 * rSges(4,1) - t63 * rSges(4,2) + t83 * rSges(4,3);
t140 = t83 * rSges(5,1) - t64 * rSges(5,2) + t63 * rSges(5,3);
t137 = t64 * pkin(3) + qJ(4) * t63;
t53 = t61 * qJ(4);
t131 = t108 * (pkin(3) * t62 + t53) + t110 * t137 + t156;
t130 = -rSges(6,2) * t80 + t176 * t147 - t154 * t81 + t157;
t129 = -rSges(5,1) * t82 - rSges(5,3) * t61;
t128 = t139 + t144;
t125 = rSges(3,1) * t109 - rSges(3,2) * t107;
t120 = Icges(3,5) * t109 - Icges(3,6) * t107;
t119 = rSges(3,1) * t145 - rSges(3,2) * t146 + t108 * rSges(3,3);
t116 = t63 * rSges(6,2) + t154 * t64 + t176 * t83;
t115 = rSges(4,1) * t62 - rSges(4,2) * t61 + rSges(4,3) * t82;
t114 = t128 + t137;
t104 = t110 * pkin(7);
t112 = t104 + t99 + (-pkin(1) - t117) * t108;
t111 = -t53 + t112;
t94 = rSges(2,1) * t110 - t108 * rSges(2,2);
t93 = -t108 * rSges(2,1) - rSges(2,2) * t110;
t92 = rSges(3,1) * t107 + rSges(3,2) * t109;
t71 = Icges(3,3) * t108 + t110 * t120;
t70 = -Icges(3,3) * t110 + t108 * t120;
t68 = t119 + t144;
t67 = t155 + t104 + (-pkin(1) - t125) * t108;
t38 = t110 * t119 + (t108 * t125 - t155) * t108;
t37 = t158 * t110;
t36 = t158 * t108;
t14 = t142 * t110;
t13 = t142 * t108;
t12 = t128 + t141;
t11 = t112 - t115;
t9 = t130 * t110;
t8 = t130 * t108;
t7 = t114 + t140;
t6 = (rSges(5,2) - pkin(3)) * t62 + t111 + t129;
t5 = t108 * t115 + t110 * t141 + t156;
t4 = t114 + t116;
t3 = (-pkin(3) - t154) * t62 + t111 + t165;
t2 = t108 * (-rSges(5,2) * t62 - t129) + t110 * t140 + t131;
t1 = t116 * t110 + (t154 * t62 - t165) * t108 + t131;
t10 = [Icges(2,3) + t168 * t81 + t167 * t80 + (Icges(3,2) * t109 - t166 * t106) * t109 + m(6) * (t3 ^ 2 + t4 ^ 2) + m(4) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + m(3) * (t67 ^ 2 + t68 ^ 2) + m(2) * (t93 ^ 2 + t94 ^ 2) + (Icges(3,1) * t107 + 0.2e1 * Icges(3,4) * t109) * t107; m(6) * (t3 * t9 + t4 * t8) + m(4) * (t11 * t37 + t12 * t36) + m(5) * (t13 * t7 + t14 * t6) + m(3) * (-t108 * t68 - t110 * t67) * t92 + (t164 / 0.2e1 + t163 / 0.2e1) * t185 + (t185 * t108 - t169 * t147 + t166 * t83 + t167 * t63 + t168 * t64 + t171 * t81 + t173 * t80) * t108 / 0.2e1 - (-Icges(3,5) * t146 - Icges(3,6) * t145 - t170 * t147 + t166 * t82 + t167 * t61 + t168 * t62 + t172 * t81 + t174 * t80) * t110 / 0.2e1; m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2 + t2 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2 + t5 ^ 2) + m(3) * (t38 ^ 2 + (t163 + t164) * t92 ^ 2) + (-t163 * t70 + (t170 * t82 + t172 * t62 + t174 * t61) * t110) * t110 + (t164 * t71 + (t110 * t71 - t169 * t82 - t170 * t83 - t171 * t62 - t172 * t64 - t173 * t61 - t174 * t63) * t110 + (-t110 * t70 + t169 * t83 + t171 * t64 + t173 * t63) * t108) * t108; m(6) * (t3 * t83 + t4 * t82) + m(4) * (t11 * t83 + t12 * t82) + m(5) * (t6 * t83 + t7 * t82); m(6) * (-t1 * t147 + t8 * t82 + t83 * t9) + m(5) * (t13 * t82 + t14 * t83 - t147 * t2) + m(4) * (-t147 * t5 + t36 * t82 + t37 * t83); 0.2e1 * (m(4) / 0.2e1 + t143) * (t106 ^ 2 * t109 ^ 2 + t82 ^ 2 + t83 ^ 2); m(6) * (t3 * t63 + t4 * t61) + m(5) * (t6 * t63 + t61 * t7); m(6) * (t1 * t80 + t61 * t8 + t63 * t9) + m(5) * (t13 * t61 + t14 * t63 + t2 * t80); (-t147 * t80 + t61 * t82 + t63 * t83) * t175; (t61 ^ 2 + t63 ^ 2 + t80 ^ 2) * t175; m(6) * (t3 * t64 + t4 * t62); m(6) * (t1 * t81 + t62 * t8 + t64 * t9); m(6) * (-t147 * t81 + t62 * t82 + t64 * t83); m(6) * (t61 * t62 + t63 * t64 + t80 * t81); m(6) * (t62 ^ 2 + t64 ^ 2 + t81 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
