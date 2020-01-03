% Calculate joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:04
% EndTime: 2019-12-31 20:29:11
% DurationCPUTime: 2.46s
% Computational Cost: add. (3475->306), mult. (8767->460), div. (0->0), fcn. (10359->8), ass. (0->162)
t236 = Icges(3,1) + Icges(4,1);
t234 = Icges(3,5) + Icges(4,4);
t135 = sin(qJ(2));
t235 = (Icges(3,4) - Icges(4,5)) * t135;
t233 = Icges(4,2) + Icges(3,3);
t138 = cos(qJ(2));
t232 = t234 * t138 + (-Icges(3,6) + Icges(4,6)) * t135;
t231 = t236 * t138 - t235;
t136 = sin(qJ(1));
t230 = -t136 / 0.2e1;
t229 = t136 / 0.2e1;
t139 = cos(qJ(1));
t228 = -t139 / 0.2e1;
t224 = t139 / 0.2e1;
t134 = sin(qJ(4));
t214 = cos(qJ(4));
t174 = t135 * t214;
t104 = -t138 * t134 + t174;
t96 = t104 * t136;
t227 = t96 / 0.2e1;
t186 = Icges(5,3) * t139;
t196 = Icges(5,6) * t96;
t198 = Icges(5,2) * t96;
t103 = t135 * t134 + t138 * t214;
t97 = t103 * t136;
t201 = Icges(5,5) * t97;
t203 = Icges(5,4) * t97;
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t99 = t103 * t139;
t72 = -t99 * t133 - t136 * t137;
t73 = -t136 * t133 + t99 * t137;
t182 = t138 * t139;
t98 = t134 * t182 - t139 * t174;
t37 = Icges(6,5) * t73 + Icges(6,6) * t72 + Icges(6,3) * t98;
t38 = Icges(6,4) * t73 + Icges(6,2) * t72 + Icges(6,6) * t98;
t39 = Icges(6,1) * t73 + Icges(6,4) * t72 + Icges(6,5) * t98;
t70 = -t97 * t133 + t139 * t137;
t71 = t139 * t133 + t97 * t137;
t11 = -t96 * t37 + t70 * t38 + t71 * t39;
t195 = Icges(6,6) * t96;
t197 = Icges(6,2) * t70;
t202 = Icges(6,4) * t71;
t143 = (-0.2e1 * t195 + t197 + 0.2e1 * t202) * t70;
t199 = Icges(6,5) * t96;
t204 = Icges(6,1) * t71;
t3 = t11 * t136 - (Icges(6,3) * t96 ^ 2 + (-0.2e1 * t199 + t204) * t71 + t143) * t139;
t53 = Icges(5,5) * t99 - Icges(5,6) * t98 - Icges(5,3) * t136;
t54 = Icges(5,4) * t99 - Icges(5,2) * t98 - Icges(5,6) * t136;
t55 = Icges(5,1) * t99 - Icges(5,4) * t98 - Icges(5,5) * t136;
t218 = -(t139 * t53 + t96 * t54 + t97 * t55) * t136 + (Icges(5,1) * t97 ^ 2 - (-t198 - 0.2e1 * t203) * t96 + (t186 + 0.2e1 * t196 + 0.2e1 * t201) * t139) * t139 - t3;
t144 = Icges(5,6) * t139 + t198 + t203;
t145 = Icges(5,1) * t97 + Icges(5,4) * t96 + Icges(5,5) * t139;
t12 = t98 * t37 + t72 * t38 + t73 * t39;
t194 = Icges(6,3) * t96;
t200 = Icges(6,5) * t71;
t146 = Icges(6,6) * t70 - t194 + t200;
t147 = -t195 + t197 + t202;
t148 = Icges(6,4) * t70 - t199 + t204;
t140 = t98 * t146 + t72 * t147 + t73 * t148;
t4 = t12 * t136 - t140 * t139;
t217 = (-t136 * t53 - t98 * t54 + t99 * t55) * t136 - (t99 * t145 - t98 * t144 - t136 * (t186 + t196 + t201)) * t139 + t4;
t226 = -t103 / 0.2e1;
t223 = -t232 * t136 + t233 * t139;
t222 = t233 * t136 + t232 * t139;
t41 = t73 * rSges(6,1) + t72 * rSges(6,2) + t98 * rSges(6,3);
t209 = t99 * pkin(4) + t98 * pkin(8) + t41;
t216 = t97 * pkin(4);
t167 = -t71 * rSges(6,1) - t70 * rSges(6,2);
t40 = -t96 * rSges(6,3) - t167;
t15 = -(-t96 * pkin(8) + t216 + t40) * t136 - t209 * t139;
t221 = -t217 * t136 - t218 * t139;
t131 = t136 ^ 2;
t132 = t139 ^ 2;
t220 = m(4) / 0.2e1;
t219 = m(6) / 0.2e1;
t215 = -rSges(5,3) - pkin(7);
t114 = t135 * rSges(3,1) + t138 * rSges(3,2);
t213 = m(3) * t114;
t212 = t136 * pkin(7);
t211 = t139 * pkin(7);
t48 = Icges(6,3) * t103 + (t137 * Icges(6,5) - t133 * Icges(6,6)) * t104;
t50 = Icges(6,5) * t103 + (t137 * Icges(6,1) - t133 * Icges(6,4)) * t104;
t208 = t104 * t137 * t50 + t103 * t48;
t51 = t103 * rSges(6,3) + (rSges(6,1) * t137 - rSges(6,2) * t133) * t104;
t207 = t104 * pkin(4) + t103 * pkin(8) + t51;
t206 = t99 * rSges(5,1) - t98 * rSges(5,2);
t184 = t135 * t139;
t180 = pkin(2) * t182 + qJ(3) * t184;
t185 = qJ(3) * t135;
t205 = t131 * (pkin(2) * t138 + t185) + t139 * t180;
t49 = Icges(6,6) * t103 + (t137 * Icges(6,4) - t133 * Icges(6,2)) * t104;
t193 = t133 * t49;
t192 = t139 * rSges(4,2);
t191 = t139 * rSges(3,3);
t189 = Icges(3,4) * t138;
t187 = Icges(4,5) * t138;
t112 = t135 * pkin(2) - t138 * qJ(3);
t181 = -t135 * rSges(4,1) + t138 * rSges(4,3) - t112;
t179 = t139 * pkin(1) + t136 * pkin(6);
t178 = t131 + t132;
t14 = t103 * t37 + (-t133 * t38 + t137 * t39) * t104;
t17 = t98 * t48 + t72 * t49 + t73 * t50;
t177 = t14 / 0.2e1 + t17 / 0.2e1;
t13 = t103 * t146 + (-t133 * t147 + t137 * t148) * t104;
t16 = -t96 * t48 + t70 * t49 + t71 * t50;
t176 = t16 / 0.2e1 + t13 / 0.2e1;
t175 = rSges(4,1) * t182 + t136 * rSges(4,2) + rSges(4,3) * t184;
t173 = t234 * t135 / 0.2e1 + (Icges(3,6) / 0.2e1 - Icges(4,6) / 0.2e1) * t138;
t172 = -pkin(3) * t135 - t112;
t123 = pkin(3) * t182;
t171 = t136 * (t136 * t138 * pkin(3) + t211) + t139 * (t123 - t212) + t205;
t170 = t179 + t180;
t68 = t104 * rSges(5,1) - t103 * rSges(5,2);
t169 = t172 - t68;
t168 = -t97 * rSges(5,1) - t96 * rSges(5,2);
t166 = rSges(3,1) * t138 - rSges(3,2) * t135;
t46 = t169 * t136;
t47 = t169 * t139;
t161 = t136 * t46 + t139 * t47;
t33 = t136 * t168 - t139 * t206;
t160 = t123 + t170;
t159 = t172 - t207;
t156 = -Icges(3,2) * t135 + t189;
t153 = Icges(4,3) * t135 + t187;
t152 = rSges(3,1) * t182 - rSges(3,2) * t184 + t136 * rSges(3,3);
t63 = Icges(5,5) * t104 - Icges(5,6) * t103;
t64 = Icges(5,4) * t104 - Icges(5,2) * t103;
t65 = Icges(5,1) * t104 - Icges(5,4) * t103;
t151 = t63 * t229 + t98 * t64 / 0.2e1 - t99 * t65 / 0.2e1 + t103 * t54 / 0.2e1 - t104 * t55 / 0.2e1 - t177;
t150 = t63 * t224 + t64 * t227 + t97 * t65 / 0.2e1 + t144 * t226 + t104 * t145 / 0.2e1 + t176;
t128 = t139 * pkin(6);
t142 = t128 + (-t185 - pkin(1) + (-pkin(2) - pkin(3)) * t138) * t136;
t42 = t215 * t139 + t142 + t168;
t43 = t215 * t136 + t160 + t206;
t149 = m(5) * (t136 * t43 + t139 * t42);
t1 = t16 * t103 + t11 * t98 - (Icges(6,1) * t71 ^ 2 - (-t194 + 0.2e1 * t200) * t96 + t143) * t96;
t2 = t17 * t103 + t12 * t98 - t140 * t96;
t141 = (-t13 * t139 + t14 * t136) * t226 + t2 * t230 + t1 * t224 + t3 * t227 - t98 * t4 / 0.2e1;
t116 = t139 * rSges(2,1) - t136 * rSges(2,2);
t115 = -t136 * rSges(2,1) - t139 * rSges(2,2);
t77 = t181 * t139;
t76 = t181 * t136;
t75 = t152 + t179;
t74 = t191 + t128 + (-pkin(1) - t166) * t136;
t59 = t170 + t175;
t58 = t192 + t128 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t138 + (-rSges(4,3) - qJ(3)) * t135) * t136;
t52 = t139 * t152 + (t166 * t136 - t191) * t136;
t36 = t139 * t175 + (-t192 + (rSges(4,1) * t138 + rSges(4,3) * t135) * t136) * t136 + t205;
t35 = t207 * t139;
t34 = t207 * t136;
t32 = t159 * t139;
t31 = t159 * t136;
t26 = t160 + t209 - t212;
t25 = -t216 - t211 - (-rSges(6,3) - pkin(8)) * t96 + t142 + t167;
t24 = t103 * t41 - t98 * t51;
t23 = -t103 * t40 - t96 * t51;
t22 = -t33 + t171;
t19 = t98 * t40 + t96 * t41;
t18 = (-t104 * t193 + t208) * t103;
t10 = -t15 + t171;
t5 = [-t103 * t64 + Icges(2,3) + (t65 - t193) * t104 + m(6) * (t25 ^ 2 + t26 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2) + m(3) * (t74 ^ 2 + t75 ^ 2) + m(2) * (t115 ^ 2 + t116 ^ 2) + t208 + ((Icges(3,2) + Icges(4,3)) * t138 + t235) * t138 + (t236 * t135 - t187 + t189) * t135; m(6) * (t32 * t25 + t31 * t26) + m(5) * (t47 * t42 + t46 * t43) + m(4) * (t77 * t58 + t76 * t59) + (-t74 * t213 + t173 * t139 + (Icges(3,6) * t224 + Icges(4,6) * t228 + t153 * t229 + t156 * t230) * t138 + (t234 * t224 + t231 * t230) * t135 - t150) * t139 + (-t75 * t213 + t173 * t136 + (Icges(3,6) * t229 + Icges(4,6) * t230 + t153 * t228 + t156 * t224) * t138 + (t231 * t224 + t234 * t229) * t135 - t151) * t136; m(6) * (t10 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t22 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(4) * (t36 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (t178 * t114 ^ 2 + t52 ^ 2) + (t223 * t132 + t218) * t139 + (t222 * t131 + (t223 * t136 + t222 * t139) * t139 + t217) * t136; 0.2e1 * ((t136 * t26 + t139 * t25) * t219 + t149 / 0.2e1 + (t136 * t59 + t139 * t58) * t220) * t135; m(6) * (-t138 * t10 + (t136 * t31 + t139 * t32) * t135) + m(5) * (t161 * t135 - t138 * t22) + m(4) * (-t138 * t36 + (t136 * t76 + t139 * t77) * t135); 0.2e1 * (t220 + m(5) / 0.2e1 + t219) * (t178 * t135 ^ 2 + t138 ^ 2); t150 * t139 + t151 * t136 + m(6) * (t35 * t25 + t34 * t26) + t68 * t149; m(6) * (t15 * t10 + t34 * t31 + t35 * t32) + m(5) * (t161 * t68 + t33 * t22) + t221; m(5) * (t178 * t68 * t135 - t33 * t138) + m(6) * (-t15 * t138 + (t136 * t34 + t139 * t35) * t135); m(5) * (t178 * t68 ^ 2 + t33 ^ 2) + m(6) * (t15 ^ 2 + t34 ^ 2 + t35 ^ 2) - t221; m(6) * (t23 * t25 + t24 * t26) + t18 + t177 * t98 - t176 * t96; m(6) * (t19 * t10 + t23 * t32 + t24 * t31) - t141; m(6) * (-t19 * t138 + (t136 * t24 + t139 * t23) * t135); m(6) * (t19 * t15 + t23 * t35 + t24 * t34) + t141; m(6) * (t19 ^ 2 + t23 ^ 2 + t24 ^ 2) + t98 * t2 - t96 * t1 + t103 * (-t13 * t96 + t14 * t98 + t18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
