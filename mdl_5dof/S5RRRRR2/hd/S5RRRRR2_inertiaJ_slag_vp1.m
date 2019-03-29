% Calculate joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:52
% EndTime: 2019-03-29 15:25:55
% DurationCPUTime: 1.23s
% Computational Cost: add. (5357->250), mult. (4728->386), div. (0->0), fcn. (4920->10), ass. (0->142)
t140 = qJ(1) + qJ(2);
t135 = sin(t140);
t132 = t135 ^ 2;
t137 = cos(t140);
t133 = t137 ^ 2;
t171 = t132 + t133;
t139 = qJ(3) + qJ(4);
t134 = sin(t139);
t136 = cos(t139);
t181 = Icges(5,4) * t136;
t157 = -Icges(5,2) * t134 + t181;
t72 = Icges(5,6) * t135 + t137 * t157;
t182 = Icges(5,4) * t134;
t159 = Icges(5,1) * t136 - t182;
t74 = Icges(5,5) * t135 + t137 * t159;
t163 = -t134 * t72 + t136 * t74;
t71 = -Icges(5,6) * t137 + t135 * t157;
t73 = -Icges(5,5) * t137 + t135 * t159;
t164 = t134 * t71 - t136 * t73;
t155 = Icges(5,5) * t136 - Icges(5,6) * t134;
t69 = -Icges(5,3) * t137 + t135 * t155;
t70 = Icges(5,3) * t135 + t137 * t155;
t180 = t134 * t135;
t144 = cos(qJ(5));
t175 = t137 * t144;
t141 = sin(qJ(5));
t178 = t135 * t141;
t95 = -t136 * t178 - t175;
t176 = t137 * t141;
t177 = t135 * t144;
t96 = t136 * t177 - t176;
t52 = Icges(6,5) * t96 + Icges(6,6) * t95 + Icges(6,3) * t180;
t54 = Icges(6,4) * t96 + Icges(6,2) * t95 + Icges(6,6) * t180;
t56 = Icges(6,1) * t96 + Icges(6,4) * t95 + Icges(6,5) * t180;
t13 = t180 * t52 + t54 * t95 + t56 * t96;
t179 = t134 * t137;
t97 = -t136 * t176 + t177;
t98 = t136 * t175 + t178;
t53 = Icges(6,5) * t98 + Icges(6,6) * t97 + Icges(6,3) * t179;
t55 = Icges(6,4) * t98 + Icges(6,2) * t97 + Icges(6,6) * t179;
t57 = Icges(6,1) * t98 + Icges(6,4) * t97 + Icges(6,5) * t179;
t14 = t180 * t53 + t55 * t95 + t57 * t96;
t8 = -t13 * t137 + t135 * t14;
t198 = -t133 * t69 - (t163 * t135 + (t164 - t70) * t137) * t135 - t8;
t89 = -rSges(6,3) * t136 + (rSges(6,1) * t144 - rSges(6,2) * t141) * t134;
t197 = m(6) * t89;
t196 = t135 / 0.2e1;
t195 = -t137 / 0.2e1;
t15 = t179 * t52 + t54 * t97 + t56 * t98;
t16 = t179 * t53 + t55 * t97 + t57 * t98;
t9 = t135 * t16 - t137 * t15;
t194 = (t132 * t70 + t9 + (t164 * t137 + (t163 - t69) * t135) * t137) * t135;
t142 = sin(qJ(3));
t145 = cos(qJ(3));
t119 = rSges(4,1) * t142 + rSges(4,2) * t145;
t193 = m(4) * t119;
t107 = rSges(5,1) * t134 + rSges(5,2) * t136;
t192 = m(5) * t107;
t143 = sin(qJ(1));
t191 = pkin(1) * t143;
t190 = pkin(2) * t142;
t189 = pkin(2) * t145;
t165 = -rSges(6,1) * t96 - rSges(6,2) * t95;
t58 = rSges(6,3) * t180 - t165;
t59 = rSges(6,1) * t98 + rSges(6,2) * t97 + rSges(6,3) * t179;
t31 = t135 * t58 + t137 * t59;
t188 = rSges(5,1) * t136;
t151 = -rSges(5,2) * t179 + rSges(5,3) * t135 + t137 * t188;
t173 = rSges(5,2) * t180 + rSges(5,3) * t137;
t42 = t135 * (t135 * t188 - t173) + t137 * t151;
t187 = rSges(4,2) * t142;
t21 = -t136 * t52 + (-t141 * t54 + t144 * t56) * t134;
t186 = t21 * t137;
t22 = -t136 * t53 + (-t141 * t55 + t144 * t57) * t134;
t185 = t22 * t135;
t184 = Icges(4,4) * t142;
t183 = Icges(4,4) * t145;
t174 = t137 * t145;
t172 = t171 * t189;
t170 = t180 / 0.2e1;
t169 = t179 / 0.2e1;
t168 = -t89 - t190;
t167 = -t107 - t190;
t109 = rSges(3,1) * t137 - rSges(3,2) * t135;
t126 = pkin(2) * t174;
t44 = t126 + t59;
t80 = -Icges(6,3) * t136 + (Icges(6,5) * t144 - Icges(6,6) * t141) * t134;
t83 = -Icges(6,6) * t136 + (Icges(6,4) * t144 - Icges(6,2) * t141) * t134;
t86 = -Icges(6,5) * t136 + (Icges(6,1) * t144 - Icges(6,4) * t141) * t134;
t29 = t180 * t80 + t83 * t95 + t86 * t96;
t3 = -t136 * t29 + (t13 * t135 + t137 * t14) * t134;
t30 = t179 * t80 + t83 * t97 + t86 * t98;
t4 = -t136 * t30 + (t135 * t15 + t137 * t16) * t134;
t166 = t9 * t169 + t8 * t170 + t3 * t195 + t4 * t196 - t136 * (t185 - t186) / 0.2e1;
t108 = -rSges(3,1) * t135 - rSges(3,2) * t137;
t160 = Icges(4,1) * t145 - t184;
t158 = -Icges(4,2) * t142 + t183;
t156 = Icges(4,5) * t145 - Icges(4,6) * t142;
t105 = Icges(5,2) * t136 + t182;
t106 = Icges(5,1) * t134 + t181;
t154 = -t105 * t134 + t106 * t136;
t117 = Icges(4,2) * t145 + t184;
t118 = Icges(4,1) * t142 + t183;
t153 = -t117 * t142 + t118 * t145;
t152 = t137 * t198 + t194;
t90 = -t137 * rSges(4,3) + (rSges(4,1) * t145 - t187) * t135;
t92 = rSges(4,1) * t174 + rSges(4,3) * t135 - t137 * t187;
t32 = -t136 * t80 + (-t141 * t83 + t144 * t86) * t134;
t150 = -t32 * t136 + (t21 + t29) * t170 + (t22 + t30) * t169;
t65 = t126 + t151;
t104 = Icges(5,5) * t134 + Icges(5,6) * t136;
t149 = -t186 / 0.2e1 + t185 / 0.2e1 + (t104 * t135 + t134 * t74 + t136 * t72 + t137 * t154 + t30) * t196 + (-t104 * t137 + t134 * t73 + t135 * t154 + t136 * t71 + t29) * t195;
t64 = (-t188 - t189) * t135 + t173;
t116 = Icges(4,5) * t142 + Icges(4,6) * t145;
t148 = t149 + (t116 * t135 + t153 * t137 + t142 * (Icges(4,5) * t135 + t137 * t160) + t145 * (Icges(4,6) * t135 + t137 * t158)) * t196 + (-t116 * t137 + t153 * t135 + t142 * (-Icges(4,5) * t137 + t135 * t160) + t145 * (-Icges(4,6) * t137 + t135 * t158)) * t195;
t43 = (-rSges(6,3) * t134 - t189) * t135 + t165;
t147 = t105 * t136 + t106 * t134 + t117 * t145 + t118 * t142 + Icges(3,3) + t32;
t146 = cos(qJ(1));
t138 = t146 * pkin(1);
t121 = rSges(2,1) * t146 - rSges(2,2) * t143;
t120 = -rSges(2,1) * t143 - rSges(2,2) * t146;
t102 = t109 + t138;
t101 = t108 - t191;
t82 = Icges(4,3) * t135 + t137 * t156;
t81 = -Icges(4,3) * t137 + t135 * t156;
t78 = t167 * t137;
t77 = t167 * t135;
t76 = t138 + t92;
t75 = -t90 - t191;
t63 = t168 * t137;
t62 = t168 * t135;
t61 = t138 + t65;
t60 = t64 - t191;
t49 = t135 * t90 + t137 * t92;
t41 = t138 + t44;
t40 = t43 - t191;
t35 = t172 + t42;
t34 = -t136 * t59 - t179 * t89;
t33 = t136 * t58 + t180 * t89;
t28 = (-t135 * t59 + t137 * t58) * t134;
t23 = t172 + t31;
t1 = [Icges(2,3) + m(6) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t75 ^ 2 + t76 ^ 2) + m(3) * (t101 ^ 2 + t102 ^ 2) + m(2) * (t120 ^ 2 + t121 ^ 2) + t147; m(6) * (t40 * t43 + t41 * t44) + m(5) * (t60 * t64 + t61 * t65) + m(4) * (-t75 * t90 + t76 * t92) + m(3) * (t101 * t108 + t102 * t109) + t147; m(6) * (t43 ^ 2 + t44 ^ 2) + m(5) * (t64 ^ 2 + t65 ^ 2) + m(4) * (t90 ^ 2 + t92 ^ 2) + m(3) * (t108 ^ 2 + t109 ^ 2) + t147; (-t135 * t76 - t137 * t75) * t193 + t148 + m(6) * (t40 * t63 + t41 * t62) + m(5) * (t60 * t78 + t61 * t77); (-t135 * t92 + t137 * t90) * t193 + t148 + m(6) * (t43 * t63 + t44 * t62) + m(5) * (t64 * t78 + t65 * t77); m(6) * (t23 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t35 ^ 2 + t77 ^ 2 + t78 ^ 2) + t135 * t132 * t82 + m(4) * (t119 ^ 2 * t171 + t49 ^ 2) + t194 + (-t133 * t81 + (-t135 * t81 + t137 * t82) * t135 + t198) * t137; (-t135 * t41 - t137 * t40) * t197 + (-t135 * t61 - t137 * t60) * t192 + t149; (-t135 * t44 - t137 * t43) * t197 + (-t135 * t65 - t137 * t64) * t192 + t149; m(6) * (t23 * t31 + (-t135 * t62 - t137 * t63) * t89) + m(5) * (t35 * t42 + (-t135 * t77 - t137 * t78) * t107) + t152; m(5) * (t107 ^ 2 * t171 + t42 ^ 2) + m(6) * (t171 * t89 ^ 2 + t31 ^ 2) + t152; m(6) * (t33 * t40 + t34 * t41) + t150; m(6) * (t33 * t43 + t34 * t44) + t150; m(6) * (t23 * t28 + t33 * t63 + t34 * t62) + t166; m(6) * (t28 * t31 + (-t135 * t34 - t137 * t33) * t89) + t166; m(6) * (t28 ^ 2 + t33 ^ 2 + t34 ^ 2) + t136 ^ 2 * t32 + (t137 * t4 + t135 * t3 - t136 * (t135 * t21 + t137 * t22)) * t134;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11); t1(2) t1(3) t1(5) t1(8) t1(12); t1(4) t1(5) t1(6) t1(9) t1(13); t1(7) t1(8) t1(9) t1(10) t1(14); t1(11) t1(12) t1(13) t1(14) t1(15);];
Mq  = res;
