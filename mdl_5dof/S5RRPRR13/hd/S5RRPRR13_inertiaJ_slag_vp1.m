% Calculate joint inertia matrix for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR13_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR13_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:09
% EndTime: 2019-12-31 20:32:16
% DurationCPUTime: 2.56s
% Computational Cost: add. (6493->394), mult. (7195->575), div. (0->0), fcn. (7735->10), ass. (0->194)
t180 = sin(qJ(1));
t243 = t180 / 0.2e1;
t182 = cos(qJ(1));
t240 = t182 / 0.2e1;
t177 = cos(pkin(9));
t162 = t177 * pkin(3) + pkin(2);
t172 = pkin(9) + qJ(4);
t164 = cos(t172);
t140 = pkin(4) * t164 + t162;
t163 = sin(t172);
t176 = sin(pkin(9));
t142 = t176 * pkin(3) + pkin(4) * t163;
t181 = cos(qJ(2));
t223 = t181 * t182;
t165 = qJ(5) + t172;
t161 = cos(t165);
t160 = sin(t165);
t221 = t182 * t160;
t111 = t180 * t161 - t181 * t221;
t220 = t182 * t161;
t112 = t180 * t160 + t181 * t220;
t179 = sin(qJ(2));
t226 = t179 * t182;
t72 = t112 * rSges(6,1) + t111 * rSges(6,2) + rSges(6,3) * t226;
t248 = t140 * t223 + t180 * t142 + t72;
t174 = t180 ^ 2;
t175 = t182 ^ 2;
t247 = m(4) / 0.2e1;
t246 = m(5) / 0.2e1;
t245 = m(6) / 0.2e1;
t227 = t179 * t180;
t224 = t180 * t181;
t109 = -t160 * t224 - t220;
t110 = t161 * t224 - t221;
t65 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t227;
t67 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t227;
t69 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t227;
t20 = t109 * t67 + t110 * t69 + t65 * t227;
t66 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t226;
t68 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t226;
t70 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t226;
t21 = t109 * t68 + t110 * t70 + t66 * t227;
t100 = -Icges(6,5) * t181 + (Icges(6,1) * t161 - Icges(6,4) * t160) * t179;
t98 = -Icges(6,3) * t181 + (Icges(6,5) * t161 - Icges(6,6) * t160) * t179;
t99 = -Icges(6,6) * t181 + (Icges(6,4) * t161 - Icges(6,2) * t160) * t179;
t38 = t110 * t100 + t109 * t99 + t98 * t227;
t5 = -t38 * t181 + (t180 * t20 + t182 * t21) * t179;
t22 = t111 * t67 + t112 * t69 + t65 * t226;
t23 = t111 * t68 + t112 * t70 + t66 * t226;
t39 = t112 * t100 + t111 * t99 + t98 * t226;
t6 = -t39 * t181 + (t180 * t22 + t182 * t23) * t179;
t244 = t6 * t226 + t5 * t227;
t242 = -t181 / 0.2e1;
t241 = -t182 / 0.2e1;
t148 = t179 * rSges(3,1) + t181 * rSges(3,2);
t239 = m(3) * t148;
t238 = pkin(2) * t181;
t237 = -pkin(2) + t162;
t178 = -pkin(7) - qJ(3);
t171 = -pkin(8) + t178;
t206 = t171 - t178;
t225 = t180 * t176;
t210 = -pkin(3) * t225 - t162 * t223;
t236 = -t206 * t226 + t210 + t248;
t101 = -t181 * rSges(6,3) + (rSges(6,1) * t161 - rSges(6,2) * t160) * t179;
t190 = -t110 * rSges(6,1) - t109 * rSges(6,2);
t71 = rSges(6,3) * t227 - t190;
t51 = t101 * t227 + t181 * t71;
t235 = t160 * t99;
t234 = t182 * rSges(3,3);
t93 = t179 * t161 * t100;
t49 = -t179 * t235 - t181 * t98 + t93;
t233 = t49 * t181;
t211 = t140 - t162;
t90 = t211 * t179 + t206 * t181;
t232 = -t101 - t90;
t231 = Icges(3,4) * t179;
t230 = Icges(3,4) * t181;
t229 = qJ(3) * t179;
t104 = -Icges(5,6) * t181 + (Icges(5,4) * t164 - Icges(5,2) * t163) * t179;
t228 = t163 * t104;
t222 = t182 * t142;
t219 = t182 * t163;
t218 = t182 * t164;
t217 = t182 * t176;
t216 = t182 * t177;
t147 = t179 * pkin(2) - t181 * qJ(3);
t215 = -(qJ(3) + t178) * t181 - t237 * t179 - t147;
t214 = t181 * rSges(4,3) - (rSges(4,1) * t177 - rSges(4,2) * t176) * t179 - t147;
t208 = pkin(2) * t223 + qJ(3) * t226;
t212 = t174 * (t229 + t238) + t182 * t208;
t209 = pkin(3) * t217 + t178 * t227;
t207 = t182 * pkin(1) + t180 * pkin(6);
t205 = t174 + t175;
t119 = -t163 * t224 - t218;
t120 = t164 * t224 - t219;
t73 = Icges(5,5) * t120 + Icges(5,6) * t119 + Icges(5,3) * t227;
t75 = Icges(5,4) * t120 + Icges(5,2) * t119 + Icges(5,6) * t227;
t77 = Icges(5,1) * t120 + Icges(5,4) * t119 + Icges(5,5) * t227;
t34 = -t181 * t73 + (-t163 * t75 + t164 * t77) * t179;
t103 = -Icges(5,3) * t181 + (Icges(5,5) * t164 - Icges(5,6) * t163) * t179;
t105 = -Icges(5,5) * t181 + (Icges(5,1) * t164 - Icges(5,4) * t163) * t179;
t42 = t103 * t227 + t119 * t104 + t120 * t105;
t204 = t34 / 0.2e1 + t42 / 0.2e1;
t121 = t180 * t164 - t181 * t219;
t122 = t180 * t163 + t181 * t218;
t74 = Icges(5,5) * t122 + Icges(5,6) * t121 + Icges(5,3) * t226;
t76 = Icges(5,4) * t122 + Icges(5,2) * t121 + Icges(5,6) * t226;
t78 = Icges(5,1) * t122 + Icges(5,4) * t121 + Icges(5,5) * t226;
t35 = -t181 * t74 + (-t163 * t76 + t164 * t78) * t179;
t43 = t103 * t226 + t121 * t104 + t122 * t105;
t203 = t35 / 0.2e1 + t43 / 0.2e1;
t108 = -t181 * rSges(5,3) + (rSges(5,1) * t164 - rSges(5,2) * t163) * t179;
t202 = -t108 + t215;
t80 = t122 * rSges(5,1) + t121 * rSges(5,2) + rSges(5,3) * t226;
t137 = t180 * t177 - t181 * t217;
t138 = t181 * t216 + t225;
t201 = t138 * rSges(4,1) + t137 * rSges(4,2) + rSges(4,3) * t226;
t200 = t227 / 0.2e1;
t199 = t226 / 0.2e1;
t32 = -t181 * t65 + (-t160 * t67 + t161 * t69) * t179;
t33 = -t181 * t66 + (-t160 * t68 + t161 * t70) * t179;
t198 = (t32 + t38) * t200 + (t33 + t39) * t199;
t183 = -t178 * t226 - t210;
t197 = t180 * ((t237 * t181 - t229) * t180 - t209) + t182 * (t183 - t208) + t212;
t9 = -t233 + (t180 * t32 + t182 * t33) * t179;
t196 = -t181 * t9 + t244;
t195 = t215 + t232;
t12 = t21 * t180 - t20 * t182;
t13 = t23 * t180 - t22 * t182;
t194 = t12 * t200 + t13 * t199 + t5 * t241 + t6 * t243 + (t33 * t180 - t32 * t182) * t242;
t193 = rSges(3,1) * t181 - rSges(3,2) * t179;
t135 = -t176 * t224 - t216;
t136 = t177 * t224 - t217;
t192 = -t136 * rSges(4,1) - t135 * rSges(4,2);
t191 = -t120 * rSges(5,1) - t119 * rSges(5,2);
t189 = Icges(3,1) * t181 - t231;
t188 = -Icges(3,2) * t179 + t230;
t187 = Icges(3,5) * t181 - Icges(3,6) * t179;
t184 = rSges(3,1) * t223 - rSges(3,2) * t226 + t180 * rSges(3,3);
t169 = t182 * pkin(6);
t150 = t182 * rSges(2,1) - t180 * rSges(2,2);
t149 = -t180 * rSges(2,1) - t182 * rSges(2,2);
t144 = Icges(3,5) * t179 + Icges(3,6) * t181;
t124 = Icges(3,3) * t180 + t187 * t182;
t123 = -Icges(3,3) * t182 + t187 * t180;
t117 = -Icges(4,5) * t181 + (Icges(4,1) * t177 - Icges(4,4) * t176) * t179;
t116 = -Icges(4,6) * t181 + (Icges(4,4) * t177 - Icges(4,2) * t176) * t179;
t97 = t184 + t207;
t96 = t234 + t169 + (-pkin(1) - t193) * t180;
t94 = t179 * t164 * t105;
t92 = t214 * t182;
t91 = t214 * t180;
t89 = Icges(4,1) * t138 + Icges(4,4) * t137 + Icges(4,5) * t226;
t88 = Icges(4,1) * t136 + Icges(4,4) * t135 + Icges(4,5) * t227;
t87 = Icges(4,4) * t138 + Icges(4,2) * t137 + Icges(4,6) * t226;
t86 = Icges(4,4) * t136 + Icges(4,2) * t135 + Icges(4,6) * t227;
t85 = Icges(4,5) * t138 + Icges(4,6) * t137 + Icges(4,3) * t226;
t84 = Icges(4,5) * t136 + Icges(4,6) * t135 + Icges(4,3) * t227;
t83 = t182 * t184 + (t193 * t180 - t234) * t180;
t79 = rSges(5,3) * t227 - t191;
t63 = t71 * t226;
t62 = t201 + t207 + t208;
t61 = t169 + (-t238 - pkin(1) + (-rSges(4,3) - qJ(3)) * t179) * t180 + t192;
t60 = t202 * t182;
t59 = t202 * t180;
t57 = -t222 + (-t171 * t179 + t211 * t181) * t180 + t209;
t56 = -t108 * t226 - t181 * t80;
t55 = t108 * t227 + t181 * t79;
t54 = t183 + t80 + t207;
t53 = t169 + (-rSges(5,3) * t179 - t162 * t181 - pkin(1)) * t180 + t191 + t209;
t52 = -t101 * t226 - t181 * t72;
t50 = -t181 * t103 - t179 * t228 + t94;
t48 = -t171 * t226 + t207 + t248;
t47 = t222 + t169 + (-t140 * t181 - pkin(1) + (-rSges(6,3) + t171) * t179) * t180 + t190;
t46 = t195 * t182;
t45 = t195 * t180;
t44 = (-t180 * t80 + t182 * t79) * t179;
t41 = -t72 * t227 + t63;
t40 = t180 * (rSges(4,3) * t227 - t192) + t182 * t201 + t212;
t29 = t121 * t76 + t122 * t78 + t74 * t226;
t28 = t121 * t75 + t122 * t77 + t73 * t226;
t27 = t119 * t76 + t120 * t78 + t74 * t227;
t26 = t119 * t75 + t120 * t77 + t73 * t227;
t25 = -t236 * t181 + t232 * t226;
t24 = t181 * t57 + t90 * t227 + t51;
t19 = t180 * t79 + t182 * t80 + t197;
t18 = t63 + (-t236 * t180 + t182 * t57) * t179;
t16 = t236 * t182 + (t57 + t71) * t180 + t197;
t15 = t29 * t180 - t28 * t182;
t14 = t27 * t180 - t26 * t182;
t8 = -t43 * t181 + (t180 * t28 + t182 * t29) * t179;
t7 = -t42 * t181 + (t180 * t26 + t182 * t27) * t179;
t1 = [Icges(2,3) + t93 + t94 + (-t98 - t103 + t231 - (Icges(4,5) * t177 - Icges(4,6) * t176) * t179 + (Icges(3,2) + Icges(4,3)) * t181) * t181 + (Icges(3,1) * t179 - t176 * t116 + t177 * t117 - t228 + t230 - t235) * t179 + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t53 ^ 2 + t54 ^ 2) + m(4) * (t61 ^ 2 + t62 ^ 2) + m(3) * (t96 ^ 2 + t97 ^ 2) + m(2) * (t149 ^ 2 + t150 ^ 2); m(6) * (t45 * t48 + t46 * t47) + m(5) * (t60 * t53 + t59 * t54) + m(4) * (t92 * t61 + t91 * t62) + (-t38 / 0.2e1 - t32 / 0.2e1 - t135 * t116 / 0.2e1 - t136 * t117 / 0.2e1 - t96 * t239 + t144 * t240 + (t84 / 0.2e1 + Icges(3,6) * t240 - t188 * t180 / 0.2e1) * t181 - t204) * t182 + (t33 / 0.2e1 + t39 / 0.2e1 + t137 * t116 / 0.2e1 + t138 * t117 / 0.2e1 - t97 * t239 + t144 * t243 + (Icges(3,6) * t243 + t188 * t240 - t85 / 0.2e1) * t181 + t203) * t180 + ((Icges(3,5) * t180 - t176 * t87 + t177 * t89 + t189 * t182) * t243 + (-Icges(3,5) * t182 - t176 * t86 + t177 * t88 + t189 * t180) * t241) * t179; m(6) * (t16 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t19 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t40 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(3) * (t205 * t148 ^ 2 + t83 ^ 2) + (-t175 * t123 - t12 - t14 + (t135 * t86 + t136 * t88 + t84 * t227) * t182) * t182 + (t13 + t15 + t174 * t124 + (t137 * t87 + t138 * t89 + t85 * t226) * t180 + (-t180 * t123 + t182 * t124 - t135 * t87 - t136 * t89 - t137 * t86 - t138 * t88 - t84 * t226 - t85 * t227) * t182) * t180; 0.2e1 * ((t180 * t48 + t182 * t47) * t245 + (t180 * t54 + t182 * t53) * t246 + (t180 * t62 + t182 * t61) * t247) * t179; m(6) * (-t181 * t16 + (t180 * t45 + t182 * t46) * t179) + m(5) * (-t181 * t19 + (t180 * t59 + t182 * t60) * t179) + m(4) * (-t181 * t40 + (t180 * t91 + t182 * t92) * t179); 0.2e1 * (t247 + t246 + t245) * (t205 * t179 ^ 2 + t181 ^ 2); (-t49 - t50) * t181 + m(6) * (t24 * t47 + t25 * t48) + m(5) * (t55 * t53 + t56 * t54) + (t204 * t180 + t203 * t182) * t179 + t198; t8 * t243 + (t35 * t180 - t34 * t182) * t242 + t7 * t241 + (t14 * t243 + t15 * t240) * t179 + m(6) * (t18 * t16 + t24 * t46 + t25 * t45) + m(5) * (t44 * t19 + t55 * t60 + t56 * t59) + t194; m(5) * (-t44 * t181 + (t180 * t56 + t182 * t55) * t179) + m(6) * (-t18 * t181 + (t180 * t25 + t182 * t24) * t179); (t50 * t181 - t9) * t181 + m(6) * (t18 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t44 ^ 2 + t55 ^ 2 + t56 ^ 2) + (t182 * t8 + t180 * t7 - t181 * (t180 * t34 + t182 * t35)) * t179 + t244; -t233 + m(6) * (t51 * t47 + t52 * t48) + t198; m(6) * (t41 * t16 + t52 * t45 + t51 * t46) + t194; m(6) * (-t41 * t181 + (t180 * t52 + t182 * t51) * t179); m(6) * (t41 * t18 + t51 * t24 + t52 * t25) + t196; m(6) * (t41 ^ 2 + t51 ^ 2 + t52 ^ 2) + t196;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
