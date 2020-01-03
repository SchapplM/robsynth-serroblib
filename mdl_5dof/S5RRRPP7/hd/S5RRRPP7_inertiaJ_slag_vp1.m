% Calculate joint inertia matrix for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:58
% EndTime: 2019-12-31 21:04:05
% DurationCPUTime: 2.61s
% Computational Cost: add. (2836->364), mult. (6940->524), div. (0->0), fcn. (7606->6), ass. (0->171)
t163 = sin(qJ(3));
t165 = sin(qJ(1));
t166 = cos(qJ(3));
t167 = cos(qJ(2));
t168 = cos(qJ(1));
t197 = t167 * t168;
t138 = t163 * t197 - t165 * t166;
t139 = t165 * t163 + t166 * t197;
t164 = sin(qJ(2));
t199 = t164 * t168;
t228 = rSges(6,3) + qJ(5);
t229 = rSges(6,1) + pkin(4);
t206 = t138 * rSges(6,2) + t229 * t139 - t199 * t228;
t231 = Icges(3,5) * t164;
t230 = t231 / 0.2e1;
t109 = -Icges(4,3) * t167 + (Icges(4,5) * t166 - Icges(4,6) * t163) * t164;
t113 = -Icges(5,2) * t167 + (Icges(5,4) * t166 + Icges(5,6) * t163) * t164;
t227 = t109 + t113;
t198 = t165 * t167;
t136 = t163 * t198 + t166 * t168;
t137 = -t163 * t168 + t166 * t198;
t201 = t164 * t165;
t63 = Icges(6,5) * t137 + Icges(6,6) * t136 - Icges(6,3) * t201;
t69 = Icges(6,4) * t137 + Icges(6,2) * t136 - Icges(6,6) * t201;
t75 = Icges(6,1) * t137 + Icges(6,4) * t136 - Icges(6,5) * t201;
t19 = t136 * t69 + t137 * t75 - t63 * t201;
t64 = Icges(6,5) * t139 + Icges(6,6) * t138 - Icges(6,3) * t199;
t70 = Icges(6,4) * t139 + Icges(6,2) * t138 - Icges(6,6) * t199;
t76 = Icges(6,1) * t139 + Icges(6,4) * t138 - Icges(6,5) * t199;
t20 = t136 * t70 + t137 * t76 - t64 * t201;
t65 = Icges(5,5) * t137 + Icges(5,6) * t201 + Icges(5,3) * t136;
t71 = Icges(5,4) * t137 + Icges(5,2) * t201 + Icges(5,6) * t136;
t77 = Icges(5,1) * t137 + Icges(5,4) * t201 + Icges(5,5) * t136;
t21 = t136 * t65 + t137 * t77 + t71 * t201;
t66 = Icges(5,5) * t139 + Icges(5,6) * t199 + Icges(5,3) * t138;
t72 = Icges(5,4) * t139 + Icges(5,2) * t199 + Icges(5,6) * t138;
t78 = Icges(5,1) * t139 + Icges(5,4) * t199 + Icges(5,5) * t138;
t22 = t136 * t66 + t137 * t78 + t72 * t201;
t67 = Icges(4,5) * t137 - Icges(4,6) * t136 + Icges(4,3) * t201;
t73 = Icges(4,4) * t137 - Icges(4,2) * t136 + Icges(4,6) * t201;
t79 = Icges(4,1) * t137 - Icges(4,4) * t136 + Icges(4,5) * t201;
t23 = -t136 * t73 + t137 * t79 + t67 * t201;
t68 = Icges(4,5) * t139 - Icges(4,6) * t138 + Icges(4,3) * t199;
t74 = Icges(4,4) * t139 - Icges(4,2) * t138 + Icges(4,6) * t199;
t80 = Icges(4,1) * t139 - Icges(4,4) * t138 + Icges(4,5) * t199;
t24 = -t136 * t74 + t137 * t80 + t68 * t201;
t107 = Icges(6,3) * t167 + (Icges(6,5) * t166 + Icges(6,6) * t163) * t164;
t112 = Icges(6,6) * t167 + (Icges(6,4) * t166 + Icges(6,2) * t163) * t164;
t117 = Icges(6,5) * t167 + (Icges(6,1) * t166 + Icges(6,4) * t163) * t164;
t44 = -t107 * t201 + t112 * t136 + t117 * t137;
t108 = -Icges(5,6) * t167 + (Icges(5,5) * t166 + Icges(5,3) * t163) * t164;
t118 = -Icges(5,4) * t167 + (Icges(5,1) * t166 + Icges(5,5) * t163) * t164;
t45 = t108 * t136 + t113 * t201 + t118 * t137;
t114 = -Icges(4,6) * t167 + (Icges(4,4) * t166 - Icges(4,2) * t163) * t164;
t119 = -Icges(4,5) * t167 + (Icges(4,1) * t166 - Icges(4,4) * t163) * t164;
t46 = t109 * t201 - t114 * t136 + t119 * t137;
t225 = (-t44 - t45 - t46) * t167 + ((t20 + t22 + t24) * t168 + (t19 + t21 + t23) * t165) * t164;
t25 = t138 * t69 + t139 * t75 - t63 * t199;
t26 = t138 * t70 + t139 * t76 - t64 * t199;
t27 = t138 * t65 + t139 * t77 + t71 * t199;
t28 = t138 * t66 + t139 * t78 + t72 * t199;
t29 = -t138 * t73 + t139 * t79 + t67 * t199;
t30 = -t138 * t74 + t139 * t80 + t68 * t199;
t47 = -t107 * t199 + t138 * t112 + t139 * t117;
t48 = t138 * t108 + t113 * t199 + t139 * t118;
t49 = t109 * t199 - t138 * t114 + t139 * t119;
t224 = (-t47 - t48 - t49) * t167 + ((t26 + t28 + t30) * t168 + (t25 + t27 + t29) * t165) * t164;
t31 = t167 * t63 + (t163 * t69 + t166 * t75) * t164;
t33 = -t167 * t71 + (t163 * t65 + t166 * t77) * t164;
t35 = -t167 * t67 + (-t163 * t73 + t166 * t79) * t164;
t223 = -t31 - t33 - t35;
t32 = t167 * t64 + (t163 * t70 + t166 * t76) * t164;
t34 = -t167 * t72 + (t163 * t66 + t166 * t78) * t164;
t36 = -t167 * t68 + (-t163 * t74 + t166 * t80) * t164;
t222 = t32 + t34 + t36;
t200 = t164 * t166;
t202 = t163 * t164;
t221 = t167 * t107 + (t108 + t112) * t202 + (t117 + t118 + t119) * t200;
t161 = t165 ^ 2;
t220 = t167 ^ 2;
t162 = t168 ^ 2;
t219 = t165 / 0.2e1;
t218 = -t167 / 0.2e1;
t146 = rSges(3,1) * t164 + rSges(3,2) * t167;
t216 = m(3) * t146;
t215 = m(6) * t164;
t214 = pkin(2) * t167;
t85 = t139 * rSges(5,1) + rSges(5,2) * t199 + t138 * rSges(5,3);
t95 = t139 * pkin(3) + t138 * qJ(4);
t213 = -t85 - t95;
t212 = rSges(6,2) * t136;
t211 = rSges(5,3) * t136;
t210 = t168 * rSges(3,3);
t207 = t229 * t137 - t228 * t201 + t212;
t140 = (pkin(3) * t166 + qJ(4) * t163) * t164;
t127 = t136 * qJ(4);
t94 = pkin(3) * t137 + t127;
t205 = t140 * t201 + t167 * t94;
t203 = Icges(3,4) * t167;
t196 = (rSges(6,1) * t166 + rSges(6,2) * t163) * t164 + pkin(4) * t200 + t228 * t167;
t123 = -rSges(5,2) * t167 + (rSges(5,1) * t166 + rSges(5,3) * t163) * t164;
t195 = -t123 - t140;
t124 = -rSges(4,3) * t167 + (rSges(4,1) * t166 - rSges(4,2) * t163) * t164;
t149 = pkin(2) * t164 - pkin(7) * t167;
t194 = -t124 - t149;
t191 = pkin(2) * t197 + pkin(7) * t199;
t193 = t161 * (pkin(7) * t164 + t214) + t168 * t191;
t190 = t168 * pkin(1) + t165 * pkin(6);
t158 = t168 * pkin(6);
t189 = t158 - t127;
t188 = t161 + t162;
t187 = t114 * t202 + t167 * t227 - t221;
t186 = -t95 - t206;
t185 = -t140 - t196;
t184 = -t149 + t195;
t86 = t139 * rSges(4,1) - t138 * rSges(4,2) + rSges(4,3) * t199;
t183 = -pkin(1) - t214;
t182 = t165 * t94 + t168 * t95 + t193;
t181 = -t149 + t185;
t180 = t190 + t191;
t179 = rSges(3,1) * t167 - rSges(3,2) * t164;
t178 = -rSges(4,1) * t137 + rSges(4,2) * t136;
t176 = -Icges(3,2) * t164 + t203;
t175 = Icges(3,5) * t167 - Icges(3,6) * t164;
t172 = rSges(3,1) * t197 - rSges(3,2) * t199 + t165 * rSges(3,3);
t171 = t180 + t95;
t170 = t36 / 0.2e1 + t34 / 0.2e1 + t32 / 0.2e1 + t49 / 0.2e1 + t48 / 0.2e1 + t47 / 0.2e1;
t169 = t46 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + t35 / 0.2e1 + t33 / 0.2e1 + t31 / 0.2e1;
t160 = t164 ^ 2;
t148 = rSges(2,1) * t168 - t165 * rSges(2,2);
t147 = -t165 * rSges(2,1) - rSges(2,2) * t168;
t143 = Icges(3,6) * t167 + t231;
t111 = Icges(3,3) * t165 + t175 * t168;
t110 = -Icges(3,3) * t168 + t175 * t165;
t97 = t172 + t190;
t96 = t210 + t158 + (-pkin(1) - t179) * t165;
t93 = t194 * t168;
t92 = t194 * t165;
t87 = t94 * t199;
t83 = rSges(4,3) * t201 - t178;
t82 = rSges(5,1) * t137 + rSges(5,2) * t201 + t211;
t62 = t168 * t172 + (t179 * t165 - t210) * t165;
t61 = t184 * t168;
t60 = t184 * t165;
t59 = t180 + t86;
t58 = t158 + ((-rSges(4,3) - pkin(7)) * t164 + t183) * t165 + t178;
t57 = t181 * t168;
t56 = t181 * t165;
t55 = -t124 * t199 - t167 * t86;
t54 = t124 * t201 + t167 * t83;
t50 = (-t165 * t86 + t168 * t83) * t164;
t43 = t171 + t85;
t42 = -t211 + (-rSges(5,1) - pkin(3)) * t137 + ((-rSges(5,2) - pkin(7)) * t164 + t183) * t165 + t189;
t41 = t171 + t206;
t40 = -t212 + (-pkin(3) - t229) * t137 + ((-pkin(7) + t228) * t164 + t183) * t165 + t189;
t39 = t165 * t83 + t168 * t86 + t193;
t38 = t213 * t167 + t195 * t199;
t37 = t123 * t201 + t167 * t82 + t205;
t18 = t186 * t167 + t185 * t199;
t17 = t207 * t167 + t196 * t201 + t205;
t16 = t87 + (t213 * t165 + t168 * t82) * t164;
t15 = t165 * t82 + t168 * t85 + t182;
t14 = t87 + (t186 * t165 + t207 * t168) * t164;
t13 = t207 * t165 + t206 * t168 + t182;
t12 = t30 * t165 - t168 * t29;
t11 = t28 * t165 - t168 * t27;
t10 = t26 * t165 - t168 * t25;
t9 = t24 * t165 - t168 * t23;
t8 = t22 * t165 - t168 * t21;
t7 = t20 * t165 - t168 * t19;
t1 = [Icges(2,3) + (Icges(3,1) * t164 - t163 * t114 + t203) * t164 + (Icges(3,4) * t164 + Icges(3,2) * t167 - t227) * t167 + m(6) * (t40 ^ 2 + t41 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(3) * (t96 ^ 2 + t97 ^ 2) + m(2) * (t147 ^ 2 + t148 ^ 2) + t221; m(6) * (t40 * t57 + t41 * t56) + m(4) * (t58 * t93 + t59 * t92) + m(5) * (t42 * t61 + t43 * t60) + (t176 * t165 * t218 - t96 * t216 - t169 + (-Icges(3,6) * t218 + t230 + t143 / 0.2e1) * t168) * t168 + ((Icges(3,6) * t165 + t176 * t168) * t167 / 0.2e1 + t165 * t230 - t97 * t216 + t143 * t219 + t170) * t165; m(6) * (t13 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t39 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(3) * (t188 * t146 ^ 2 + t62 ^ 2) + (-t162 * t110 - t7 - t8 - t9) * t168 + (t161 * t111 + t10 + t11 + t12 + (-t165 * t110 + t168 * t111) * t168) * t165; t187 * t167 + m(6) * (t17 * t40 + t18 * t41) + m(4) * (t54 * t58 + t55 * t59) + m(5) * (t37 * t42 + t38 * t43) + (t169 * t165 + t170 * t168) * t164; m(6) * (t13 * t14 + t17 * t57 + t18 * t56) + m(5) * (t15 * t16 + t37 * t61 + t38 * t60) + m(4) * (t39 * t50 + t54 * t93 + t55 * t92) + ((t10 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t168 + (t7 / 0.2e1 + t8 / 0.2e1 + t9 / 0.2e1) * t165) * t164 + t224 * t219 + (t222 * t165 + t223 * t168) * t218 - t225 * t168 / 0.2e1; m(5) * (t16 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(6) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(4) * (t50 ^ 2 + t54 ^ 2 + t55 ^ 2) - t187 * t220 + ((-t222 * t167 + t224) * t168 + (t223 * t167 + t225) * t165) * t164; m(6) * (t136 * t41 + t138 * t40) + m(5) * (t136 * t43 + t138 * t42); m(6) * (t13 * t202 + t136 * t56 + t138 * t57) + m(5) * (t136 * t60 + t138 * t61 + t15 * t202); m(5) * (t136 * t38 + t138 * t37 + t16 * t202) + m(6) * (t136 * t18 + t138 * t17 + t14 * t202); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t160 * t163 ^ 2 + t136 ^ 2 + t138 ^ 2); (-t165 * t41 - t168 * t40) * t215; m(6) * (t167 * t13 + (-t165 * t56 - t168 * t57) * t164); m(6) * (t167 * t14 + (-t165 * t18 - t168 * t17) * t164); (-t136 * t165 - t138 * t168 + t163 * t167) * t215; m(6) * (t188 * t160 + t220);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
