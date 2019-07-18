% Calculate joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:45
% EndTime: 2019-07-18 17:20:54
% DurationCPUTime: 2.35s
% Computational Cost: add. (3080->267), mult. (4307->402), div. (0->0), fcn. (4477->8), ass. (0->145)
t224 = Icges(3,4) + Icges(4,4);
t223 = Icges(3,1) + Icges(4,1);
t222 = Icges(3,5) + Icges(4,5);
t221 = Icges(3,2) + Icges(4,2);
t220 = Icges(3,6) + Icges(4,6);
t144 = cos(qJ(2));
t219 = t224 * t144;
t141 = sin(qJ(2));
t218 = t224 * t141;
t213 = t222 * t141 + t220 * t144;
t215 = Icges(3,3) + Icges(4,3);
t214 = -t220 * t141 + t222 * t144;
t212 = t219 * t144 + (-t218 + (-t221 + t223) * t144) * t141;
t142 = sin(qJ(1));
t136 = t142 ^ 2;
t145 = cos(qJ(1));
t137 = t145 ^ 2;
t177 = t136 + t137;
t138 = qJ(2) + qJ(4);
t129 = sin(t138);
t130 = cos(t138);
t190 = Icges(5,4) * t130;
t155 = -Icges(5,2) * t129 + t190;
t73 = Icges(5,6) * t142 + t145 * t155;
t191 = Icges(5,4) * t129;
t158 = Icges(5,1) * t130 - t191;
t75 = Icges(5,5) * t142 + t145 * t158;
t165 = -t129 * t73 + t130 * t75;
t72 = -Icges(5,6) * t145 + t142 * t155;
t74 = -Icges(5,5) * t145 + t142 * t158;
t166 = t129 * t72 - t130 * t74;
t152 = Icges(5,5) * t130 - Icges(5,6) * t129;
t70 = -Icges(5,3) * t145 + t142 * t152;
t71 = Icges(5,3) * t142 + t145 * t152;
t189 = t129 * t142;
t143 = cos(qJ(5));
t182 = t143 * t145;
t140 = sin(qJ(5));
t186 = t140 * t142;
t95 = -t130 * t186 - t182;
t184 = t142 * t143;
t185 = t140 * t145;
t96 = t130 * t184 - t185;
t44 = Icges(6,5) * t96 + Icges(6,6) * t95 + Icges(6,3) * t189;
t46 = Icges(6,4) * t96 + Icges(6,2) * t95 + Icges(6,6) * t189;
t48 = Icges(6,1) * t96 + Icges(6,4) * t95 + Icges(6,5) * t189;
t13 = t189 * t44 + t46 * t95 + t48 * t96;
t188 = t129 * t145;
t97 = -t130 * t185 + t184;
t98 = t130 * t182 + t186;
t45 = Icges(6,5) * t98 + Icges(6,6) * t97 + Icges(6,3) * t188;
t47 = Icges(6,4) * t98 + Icges(6,2) * t97 + Icges(6,6) * t188;
t49 = Icges(6,1) * t98 + Icges(6,4) * t97 + Icges(6,5) * t188;
t14 = t189 * t45 + t47 * t95 + t49 * t96;
t8 = -t13 * t145 + t14 * t142;
t211 = -t137 * t70 - (t165 * t142 + (t166 - t71) * t145) * t142 - t8;
t210 = t142 * t215 + t214 * t145;
t209 = -t142 * t214 + t215 * t145;
t208 = t142 / 0.2e1;
t207 = -t145 / 0.2e1;
t15 = t188 * t44 + t46 * t97 + t48 * t98;
t16 = t188 * t45 + t47 * t97 + t49 * t98;
t9 = t142 * t16 - t145 * t15;
t206 = (t136 * t71 + t9 + (t166 * t145 + (t165 - t70) * t142) * t145) * t142;
t204 = t141 * pkin(1);
t146 = pkin(2) + pkin(1);
t203 = -pkin(1) + t146;
t187 = t130 * t145;
t148 = rSges(5,1) * t187 - rSges(5,2) * t188 + t142 * rSges(5,3);
t167 = rSges(5,1) * t130 - rSges(5,2) * t129;
t41 = t142 * (-t145 * rSges(5,3) + t142 * t167) + t145 * t148;
t202 = rSges(3,2) * t141;
t201 = rSges(4,2) * t141;
t67 = -Icges(6,6) * t130 + (Icges(6,4) * t143 - Icges(6,2) * t140) * t129;
t200 = t140 * t67;
t199 = t145 * rSges(4,3);
t20 = -t130 * t44 + (-t140 * t46 + t143 * t48) * t129;
t198 = t20 * t145;
t21 = -t130 * t45 + (-t140 * t47 + t143 * t49) * t129;
t197 = t21 * t142;
t132 = t145 * qJ(3);
t181 = t144 * t145;
t178 = pkin(1) * t181 + t142 * qJ(3);
t183 = t142 * t144;
t196 = t145 * t178 + t142 * (pkin(1) * t183 - t132);
t180 = t144 * t146;
t139 = -pkin(3) - qJ(3);
t179 = t145 * t139;
t51 = t98 * rSges(6,1) + t97 * rSges(6,2) + rSges(6,3) * t188;
t176 = -t203 * t141 - t204;
t175 = -rSges(4,1) * t141 - rSges(4,2) * t144 - t204;
t173 = -t139 * t142 + t145 * t180;
t174 = t142 * (t203 * t183 + t132 + t179) + t145 * (t173 - t178) + t196;
t168 = -t96 * rSges(6,1) - t95 * rSges(6,2);
t50 = rSges(6,3) * t189 - t168;
t27 = t177 * pkin(4) * t129 + t142 * t50 + t145 * t51;
t66 = -Icges(6,3) * t130 + (Icges(6,5) * t143 - Icges(6,6) * t140) * t129;
t68 = -Icges(6,5) * t130 + (Icges(6,1) * t143 - Icges(6,4) * t140) * t129;
t25 = t189 * t66 + t67 * t95 + t68 * t96;
t3 = -t130 * t25 + (t13 * t142 + t14 * t145) * t129;
t26 = t188 * t66 + t67 * t97 + t68 * t98;
t4 = -t130 * t26 + (t142 * t15 + t145 * t16) * t129;
t172 = t3 * t207 + t4 * t208 - t130 * (t197 - t198) / 0.2e1 + t8 * t189 / 0.2e1 + t9 * t188 / 0.2e1;
t69 = -rSges(6,3) * t130 + (rSges(6,1) * t143 - rSges(6,2) * t140) * t129;
t171 = t176 - t69;
t105 = rSges(5,1) * t129 + rSges(5,2) * t130;
t169 = -t105 + t176;
t103 = Icges(5,2) * t130 + t191;
t104 = Icges(5,1) * t129 + t190;
t151 = -t103 * t129 + t104 * t130;
t150 = t211 * t145 + t206;
t149 = rSges(4,1) * t181 + t142 * rSges(4,3) - t145 * t201;
t102 = Icges(5,5) * t129 + Icges(5,6) * t130;
t147 = -t198 / 0.2e1 + t197 / 0.2e1 + (t102 * t142 + t129 * t75 + t130 * t73 + t145 * t151 + t26) * t208 + (-t102 * t145 + t129 * t74 + t130 * t72 + t142 * t151 + t25) * t207;
t123 = pkin(4) * t187;
t121 = t142 * t130 * pkin(4);
t117 = rSges(2,1) * t145 - rSges(2,2) * t142;
t116 = -rSges(2,1) * t142 - rSges(2,2) * t145;
t115 = rSges(3,1) * t141 + rSges(3,2) * t144;
t92 = rSges(3,3) * t142 + (rSges(3,1) * t144 - t202) * t145;
t90 = rSges(3,1) * t183 - t145 * rSges(3,3) - t142 * t202;
t77 = t175 * t145;
t76 = t175 * t142;
t61 = t149 + t178;
t60 = t199 + t132 + (t201 + (-rSges(4,1) - pkin(1)) * t144) * t142;
t59 = t129 * t143 * t68;
t58 = t148 + t173;
t57 = (rSges(5,3) - t139) * t145 + (-t167 - t180) * t142;
t56 = -t145 * t69 + t123;
t55 = -t142 * t69 + t121;
t54 = t169 * t145;
t53 = t169 * t142;
t52 = t142 * t90 + t145 * t92;
t40 = t145 * t171 + t123;
t39 = t142 * t171 + t121;
t34 = pkin(4) * t188 + t173 + t51;
t33 = -t179 + (-t180 + (-rSges(6,3) - pkin(4)) * t129) * t142 + t168;
t32 = t145 * t149 + (-t199 + (rSges(4,1) * t144 - t201) * t142) * t142 + t196;
t31 = -t130 * t51 - t188 * t69;
t30 = t130 * t50 + t189 * t69;
t29 = -t129 * t200 - t130 * t66 + t59;
t28 = (-t142 * t51 + t145 * t50) * t129;
t22 = t174 + t41;
t11 = t27 + t174;
t1 = [Icges(2,3) + t59 + (t103 - t66) * t130 + (t104 - t200) * t129 + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t90 ^ 2 + t92 ^ 2) + m(2) * (t116 ^ 2 + t117 ^ 2) + (t221 * t144 + t218) * t144 + (t223 * t141 + t219) * t141; m(6) * (t33 * t40 + t34 * t39) + m(5) * (t53 * t58 + t54 * t57) + m(4) * (t60 * t77 + t61 * t76) + t147 + m(3) * (-t142 * t92 + t145 * t90) * t115 + (t142 * t213 + t212 * t145) * t208 + (t212 * t142 - t213 * t145) * t207 + t213 * (t137 / 0.2e1 + t136 / 0.2e1); m(6) * (t11 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t22 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(4) * (t32 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (t115 ^ 2 * t177 + t52 ^ 2) + t206 + t210 * t142 * t136 + (t209 * t137 + (t209 * t142 + t210 * t145) * t142 + t211) * t145; m(6) * (t142 * t33 - t145 * t34) + m(5) * (t142 * t57 - t145 * t58) + m(4) * (t142 * t60 - t145 * t61); m(6) * (t142 * t40 - t145 * t39) + m(5) * (t142 * t54 - t145 * t53) + m(4) * (t142 * t77 - t145 * t76); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t177; m(6) * (t33 * t56 + t34 * t55) + m(5) * (-t142 * t58 - t145 * t57) * t105 + t147; m(6) * (t11 * t27 + t39 * t55 + t40 * t56) + m(5) * (t22 * t41 + (-t142 * t53 - t145 * t54) * t105) + t150; m(6) * (t142 * t56 - t145 * t55); m(5) * (t105 ^ 2 * t177 + t41 ^ 2) + m(6) * (t27 ^ 2 + t55 ^ 2 + t56 ^ 2) + t150; m(6) * (t30 * t33 + t31 * t34) - t29 * t130 + ((t21 / 0.2e1 + t26 / 0.2e1) * t145 + (t20 / 0.2e1 + t25 / 0.2e1) * t142) * t129; m(6) * (t11 * t28 + t30 * t40 + t31 * t39) + t172; m(6) * (t142 * t30 - t145 * t31); m(6) * (t27 * t28 + t30 * t56 + t31 * t55) + t172; m(6) * (t28 ^ 2 + t30 ^ 2 + t31 ^ 2) + t130 ^ 2 * t29 + (t145 * t4 + t142 * t3 - t130 * (t142 * t20 + t145 * t21)) * t129;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
