% Calculate joint inertia matrix for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP13_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP13_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:36
% DurationCPUTime: 1.63s
% Computational Cost: add. (1848->267), mult. (4486->402), div. (0->0), fcn. (4832->6), ass. (0->131)
t124 = cos(qJ(3));
t182 = Icges(4,5) * t124;
t121 = sin(qJ(3));
t181 = Icges(4,6) * t121;
t178 = rSges(6,1) + pkin(4);
t177 = rSges(6,3) + qJ(5);
t120 = sin(qJ(4));
t125 = cos(qJ(1));
t146 = t125 * t120;
t122 = sin(qJ(1));
t123 = cos(qJ(4));
t150 = t122 * t123;
t93 = t121 * t146 + t150;
t145 = t125 * t123;
t151 = t122 * t120;
t95 = -t121 * t145 + t151;
t180 = -t177 * t93 + t178 * t95;
t179 = t182 / 0.2e1 - t181 / 0.2e1;
t153 = t120 * t124;
t69 = Icges(6,6) * t121 + (Icges(6,5) * t123 + Icges(6,3) * t120) * t124;
t70 = Icges(5,3) * t121 + (Icges(5,5) * t123 - Icges(5,6) * t120) * t124;
t73 = Icges(6,2) * t121 + (Icges(6,4) * t123 + Icges(6,6) * t120) * t124;
t77 = Icges(6,4) * t121 + (Icges(6,1) * t123 + Icges(6,5) * t120) * t124;
t78 = Icges(5,5) * t121 + (Icges(5,1) * t123 - Icges(5,4) * t120) * t124;
t176 = t69 * t153 + (t77 + t78) * t123 * t124 + (t70 + t73) * t121;
t149 = t122 * t124;
t91 = t121 * t151 - t145;
t92 = t121 * t150 + t146;
t42 = Icges(6,5) * t92 - Icges(6,6) * t149 + Icges(6,3) * t91;
t46 = Icges(6,4) * t92 - Icges(6,2) * t149 + Icges(6,6) * t91;
t50 = Icges(6,1) * t92 - Icges(6,4) * t149 + Icges(6,5) * t91;
t11 = -t46 * t149 + t42 * t91 + t50 * t92;
t147 = t124 * t125;
t43 = Icges(6,5) * t95 + Icges(6,6) * t147 - Icges(6,3) * t93;
t47 = Icges(6,4) * t95 + Icges(6,2) * t147 - Icges(6,6) * t93;
t51 = Icges(6,1) * t95 + Icges(6,4) * t147 - Icges(6,5) * t93;
t12 = -t47 * t149 + t43 * t91 + t51 * t92;
t44 = Icges(5,5) * t92 - Icges(5,6) * t91 - Icges(5,3) * t149;
t48 = Icges(5,4) * t92 - Icges(5,2) * t91 - Icges(5,6) * t149;
t52 = Icges(5,1) * t92 - Icges(5,4) * t91 - Icges(5,5) * t149;
t13 = -t44 * t149 - t48 * t91 + t52 * t92;
t45 = Icges(5,5) * t95 + Icges(5,6) * t93 + Icges(5,3) * t147;
t49 = Icges(5,4) * t95 + Icges(5,2) * t93 + Icges(5,6) * t147;
t53 = Icges(5,1) * t95 + Icges(5,4) * t93 + Icges(5,5) * t147;
t14 = -t45 * t149 - t49 * t91 + t53 * t92;
t28 = -t73 * t149 + t91 * t69 + t92 * t77;
t74 = Icges(5,6) * t121 + (Icges(5,4) * t123 - Icges(5,2) * t120) * t124;
t29 = -t70 * t149 - t91 * t74 + t92 * t78;
t175 = ((t12 + t14) * t125 + (-t11 - t13) * t122) * t124 + (t28 + t29) * t121;
t15 = t46 * t147 - t93 * t42 + t95 * t50;
t16 = t47 * t147 - t93 * t43 + t95 * t51;
t17 = t44 * t147 + t93 * t48 + t95 * t52;
t18 = t45 * t147 + t93 * t49 + t95 * t53;
t30 = t73 * t147 - t93 * t69 + t95 * t77;
t31 = t70 * t147 + t93 * t74 + t95 * t78;
t174 = ((t16 + t18) * t125 + (-t15 - t17) * t122) * t124 + (t30 + t31) * t121;
t19 = t121 * t46 + (t120 * t42 + t123 * t50) * t124;
t21 = t121 * t44 + (-t120 * t48 + t123 * t52) * t124;
t173 = t19 + t21;
t20 = t121 * t47 + (t120 * t43 + t123 * t51) * t124;
t22 = t121 * t45 + (-t120 * t49 + t123 * t53) * t124;
t172 = t20 + t22;
t171 = t177 * t91 + t178 * t92;
t170 = (rSges(4,1) * t121 + rSges(4,2) * t124) * t125;
t118 = t122 ^ 2;
t119 = t125 ^ 2;
t169 = -pkin(1) - pkin(6);
t166 = t122 / 0.2e1;
t164 = t125 / 0.2e1;
t104 = t124 * rSges(4,1) - t121 * rSges(4,2);
t163 = m(4) * t104;
t162 = (-t74 * t153 + t176) * t121;
t161 = -rSges(6,2) * t149 + t171;
t160 = rSges(6,2) * t147 + t180;
t158 = t121 * rSges(6,2) + (t177 * t120 + t178 * t123) * t124;
t156 = t92 * rSges(5,1) - t91 * rSges(5,2);
t152 = t121 * t122;
t144 = t125 * pkin(1) + t122 * qJ(2);
t143 = t118 + t119;
t141 = rSges(4,1) * t152 + rSges(4,2) * t149 + t125 * rSges(4,3);
t140 = t125 * pkin(6) + t144;
t139 = (-rSges(6,2) - pkin(7)) * t124;
t138 = (-rSges(5,3) - pkin(7)) * t124;
t137 = t158 * t124;
t110 = pkin(3) * t152;
t136 = t110 + t140;
t135 = -t95 * rSges(5,1) - t93 * rSges(5,2);
t129 = Icges(4,5) * t121 + Icges(4,6) * t124;
t111 = t125 * t121 * pkin(3);
t114 = t125 * qJ(2);
t128 = t169 * t122 + t111 + t114;
t127 = t20 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1 + t22 / 0.2e1;
t126 = -t29 / 0.2e1 - t28 / 0.2e1 - t21 / 0.2e1 - t19 / 0.2e1;
t106 = t124 * pkin(3) + t121 * pkin(7);
t105 = t125 * rSges(2,1) - t122 * rSges(2,2);
t103 = -t122 * rSges(2,1) - t125 * rSges(2,2);
t100 = -t181 + t182;
t98 = t122 * t106;
t97 = -pkin(7) * t149 + t110;
t85 = t125 * (pkin(7) * t147 - t111);
t84 = -t125 * rSges(3,2) + t122 * rSges(3,3) + t144;
t83 = t125 * rSges(3,3) + t114 + (rSges(3,2) - pkin(1)) * t122;
t82 = t121 * rSges(5,3) + (rSges(5,1) * t123 - rSges(5,2) * t120) * t124;
t72 = Icges(4,3) * t122 - t129 * t125;
t71 = Icges(4,3) * t125 + t129 * t122;
t63 = t140 + t141;
t62 = t114 + t170 + (-rSges(4,3) + t169) * t122;
t59 = (-t106 - t82) * t125;
t58 = t122 * t82 + t98;
t57 = rSges(5,3) * t147 - t135;
t55 = -rSges(5,3) * t149 + t156;
t41 = -t122 * t141 + (t122 * rSges(4,3) - t170) * t125;
t40 = (-t106 - t158) * t125;
t39 = t158 * t122 + t98;
t38 = t122 * t138 + t136 + t156;
t37 = t125 * t138 + t128 + t135;
t36 = -t121 * t57 + t82 * t147;
t35 = t121 * t55 + t82 * t149;
t32 = (-t122 * t57 - t125 * t55) * t124;
t27 = t122 * t139 + t136 + t171;
t26 = t125 * t139 + t128 - t180;
t25 = t125 * t57 + t85 + (-t55 - t97) * t122;
t24 = -t160 * t121 + t125 * t137;
t23 = t161 * t121 + t122 * t137;
t10 = (-t160 * t122 - t161 * t125) * t124;
t9 = t85 + t160 * t125 + (-t97 - t161) * t122;
t8 = t18 * t122 + t125 * t17;
t7 = t16 * t122 + t125 * t15;
t6 = t14 * t122 + t125 * t13;
t5 = t11 * t125 + t12 * t122;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t124 - t120 * t74) * t124 + m(6) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t83 ^ 2 + t84 ^ 2) + m(2) * (t103 ^ 2 + t105 ^ 2) + t176 + (-0.2e1 * Icges(4,4) * t124 + Icges(4,2) * t121) * t121; m(6) * (t122 * t26 - t125 * t27) + m(5) * (t122 * t37 - t125 * t38) + m(4) * (t122 * t62 - t125 * t63) + m(3) * (t122 * t83 - t125 * t84); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t143; m(6) * (t26 * t39 + t27 * t40) + m(5) * (t37 * t58 + t38 * t59) + (t100 * t164 + t179 * t125 - t63 * t163 - t126) * t125 + (t100 * t166 + t179 * t122 + t62 * t163 + t127) * t122; m(5) * (t58 * t122 - t125 * t59) + m(6) * (t39 * t122 - t125 * t40) + t143 * t163; m(6) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t143 * t104 ^ 2 + t41 ^ 2) + (t118 * t72 + t7 + t8) * t122 + (t119 * t71 + t5 + t6 + (t122 * t71 + t125 * t72) * t122) * t125; m(6) * (t23 * t27 + t24 * t26) + m(5) * (t35 * t38 + t36 * t37) + (t126 * t122 + t127 * t125) * t124 + t162; m(5) * (t36 * t122 - t125 * t35) + m(6) * (t24 * t122 - t125 * t23); m(6) * (t10 * t9 + t23 * t40 + t24 * t39) + m(5) * (t25 * t32 + t35 * t59 + t36 * t58) + ((t7 / 0.2e1 + t8 / 0.2e1) * t125 + (-t5 / 0.2e1 - t6 / 0.2e1) * t122) * t124 + (t172 * t122 + t173 * t125) * t121 / 0.2e1 + t174 * t166 + t175 * t164; t162 * t121 + m(6) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t32 ^ 2 + t35 ^ 2 + t36 ^ 2) + (t174 * t125 - t175 * t122 + (-t173 * t122 + t172 * t125) * t121) * t124; m(6) * (t26 * t91 - t27 * t93); m(6) * (t91 * t122 + t125 * t93); m(6) * (t9 * t153 + t39 * t91 - t40 * t93); m(6) * (t10 * t153 - t23 * t93 + t24 * t91); m(6) * (t120 ^ 2 * t124 ^ 2 + t91 ^ 2 + t93 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
