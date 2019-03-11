% Calculate joint inertia matrix for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:59
% EndTime: 2019-03-09 03:37:02
% DurationCPUTime: 2.54s
% Computational Cost: add. (9119->374), mult. (6633->530), div. (0->0), fcn. (7072->12), ass. (0->191)
t251 = Icges(4,3) + Icges(5,3);
t164 = qJ(3) + pkin(11);
t157 = sin(t164);
t159 = cos(t164);
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t250 = Icges(4,5) * t172 + Icges(5,5) * t159 - Icges(4,6) * t169 - Icges(5,6) * t157;
t249 = t157 / 0.2e1;
t248 = t159 / 0.2e1;
t247 = t169 / 0.2e1;
t246 = t172 / 0.2e1;
t165 = qJ(1) + pkin(10);
t158 = sin(t165);
t160 = cos(t165);
t245 = -t250 * t158 + t251 * t160;
t244 = t251 * t158 + t250 * t160;
t171 = cos(qJ(5));
t153 = t171 * pkin(5) + pkin(4);
t174 = -pkin(9) - pkin(8);
t214 = t159 * t160;
t168 = sin(qJ(5));
t216 = t158 * t168;
t219 = t157 * t160;
t166 = qJ(5) + qJ(6);
t161 = sin(t166);
t213 = t160 * t161;
t162 = cos(t166);
t217 = t158 * t162;
t118 = -t159 * t213 + t217;
t212 = t160 * t162;
t218 = t158 * t161;
t119 = t159 * t212 + t218;
t66 = t119 * rSges(7,1) + t118 * rSges(7,2) + rSges(7,3) * t219;
t243 = pkin(5) * t216 + t153 * t214 - t174 * t219 + t66;
t155 = t158 ^ 2;
t156 = t160 ^ 2;
t220 = t157 * t158;
t116 = -t159 * t218 - t212;
t117 = t159 * t217 - t213;
t59 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t220;
t61 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t220;
t63 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t220;
t19 = t116 * t61 + t117 * t63 + t220 * t59;
t60 = Icges(7,5) * t119 + Icges(7,6) * t118 + Icges(7,3) * t219;
t62 = Icges(7,4) * t119 + Icges(7,2) * t118 + Icges(7,6) * t219;
t64 = Icges(7,1) * t119 + Icges(7,4) * t118 + Icges(7,5) * t219;
t20 = t116 * t62 + t117 * t64 + t220 * t60;
t96 = -Icges(7,3) * t159 + (Icges(7,5) * t162 - Icges(7,6) * t161) * t157;
t97 = -Icges(7,6) * t159 + (Icges(7,4) * t162 - Icges(7,2) * t161) * t157;
t98 = -Icges(7,5) * t159 + (Icges(7,1) * t162 - Icges(7,4) * t161) * t157;
t39 = t116 * t97 + t117 * t98 + t220 * t96;
t5 = -t39 * t159 + (t158 * t19 + t160 * t20) * t157;
t21 = t118 * t61 + t119 * t63 + t219 * t59;
t22 = t118 * t62 + t119 * t64 + t219 * t60;
t40 = t118 * t97 + t119 * t98 + t219 * t96;
t6 = -t40 * t159 + (t158 * t21 + t160 * t22) * t157;
t242 = t6 * t219 + t5 * t220;
t241 = t158 / 0.2e1;
t240 = -t159 / 0.2e1;
t239 = -t160 / 0.2e1;
t145 = t169 * rSges(4,1) + t172 * rSges(4,2);
t238 = m(4) * t145;
t237 = pkin(3) * t169;
t236 = pkin(4) * t159;
t170 = sin(qJ(1));
t235 = t170 * pkin(1);
t234 = -pkin(4) + t153;
t233 = pkin(8) + t174;
t207 = pkin(4) * t214 + pkin(8) * t219;
t232 = -t207 + t243;
t188 = -t117 * rSges(7,1) - t116 * rSges(7,2);
t65 = rSges(7,3) * t220 - t188;
t99 = -t159 * rSges(7,3) + (rSges(7,1) * t162 - rSges(7,2) * t161) * t157;
t48 = t159 * t65 + t99 * t220;
t154 = t172 * pkin(3) + pkin(2);
t140 = t160 * t154;
t152 = t160 * pkin(7);
t167 = -qJ(4) - pkin(7);
t211 = t160 * t167;
t231 = t158 * (t211 + t152 + (-pkin(2) + t154) * t158) + t160 * (-t160 * pkin(2) + t140 + (-pkin(7) - t167) * t158);
t89 = t157 * t234 + t159 * t233;
t230 = -t89 - t99;
t229 = rSges(4,1) * t172;
t228 = rSges(4,2) * t169;
t227 = t160 * rSges(4,3);
t226 = t161 * t97;
t84 = t157 * t162 * t98;
t47 = -t157 * t226 - t159 * t96 + t84;
t225 = t47 * t159;
t224 = Icges(4,4) * t169;
t223 = Icges(4,4) * t172;
t222 = Icges(5,4) * t157;
t221 = Icges(5,4) * t159;
t215 = t158 * t171;
t210 = t160 * t168;
t209 = t160 * t171;
t109 = -Icges(6,6) * t159 + (Icges(6,4) * t171 - Icges(6,2) * t168) * t157;
t208 = t168 * t109;
t206 = t158 * rSges(4,3) + t160 * t229;
t205 = t155 + t156;
t122 = -t159 * t216 - t209;
t123 = t159 * t215 - t210;
t72 = Icges(6,5) * t123 + Icges(6,6) * t122 + Icges(6,3) * t220;
t74 = Icges(6,4) * t123 + Icges(6,2) * t122 + Icges(6,6) * t220;
t76 = Icges(6,1) * t123 + Icges(6,4) * t122 + Icges(6,5) * t220;
t34 = -t159 * t72 + (-t168 * t74 + t171 * t76) * t157;
t106 = -Icges(6,3) * t159 + (Icges(6,5) * t171 - Icges(6,6) * t168) * t157;
t112 = -Icges(6,5) * t159 + (Icges(6,1) * t171 - Icges(6,4) * t168) * t157;
t43 = t106 * t220 + t122 * t109 + t123 * t112;
t204 = t43 / 0.2e1 + t34 / 0.2e1;
t124 = -t159 * t210 + t215;
t125 = t159 * t209 + t216;
t73 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t219;
t75 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t219;
t77 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t219;
t35 = -t159 * t73 + (-t168 * t75 + t171 * t77) * t157;
t44 = t106 * t219 + t124 * t109 + t125 * t112;
t203 = t44 / 0.2e1 + t35 / 0.2e1;
t79 = t125 * rSges(6,1) + t124 * rSges(6,2) + rSges(6,3) * t219;
t202 = t220 / 0.2e1;
t201 = t219 / 0.2e1;
t200 = Icges(4,5) * t247 + Icges(5,5) * t249 + Icges(4,6) * t246 + Icges(5,6) * t248;
t199 = -t157 * rSges(5,1) - t159 * rSges(5,2) - t237;
t198 = -t157 * pkin(4) + t159 * pkin(8) - t237;
t30 = -t159 * t59 + (-t161 * t61 + t162 * t63) * t157;
t31 = -t159 * t60 + (-t161 * t62 + t162 * t64) * t157;
t197 = (t30 + t39) * t202 + (t31 + t40) * t201;
t196 = t155 * (pkin(8) * t157 + t236) + t160 * t207 + t231;
t9 = -t225 + (t158 * t30 + t160 * t31) * t157;
t195 = -t159 * t9 + t242;
t12 = t20 * t158 - t19 * t160;
t13 = t22 * t158 - t21 * t160;
t194 = t12 * t202 + t13 * t201 + t5 * t239 + t6 * t241 + (t31 * t158 - t30 * t160) * t240;
t115 = -t159 * rSges(6,3) + (rSges(6,1) * t171 - rSges(6,2) * t168) * t157;
t193 = -t115 + t198;
t173 = cos(qJ(1));
t163 = t173 * pkin(1);
t192 = -t158 * t167 + t140 + t163;
t191 = -t228 + t229;
t190 = rSges(5,1) * t159 - rSges(5,2) * t157;
t189 = -t123 * rSges(6,1) - t122 * rSges(6,2);
t185 = t198 + t230;
t184 = Icges(4,1) * t172 - t224;
t183 = Icges(5,1) * t159 - t222;
t182 = -Icges(4,2) * t169 + t223;
t181 = -Icges(5,2) * t157 + t221;
t176 = rSges(5,1) * t214 - rSges(5,2) * t219 + t158 * rSges(5,3);
t147 = t173 * rSges(2,1) - t170 * rSges(2,2);
t146 = -t170 * rSges(2,1) - t173 * rSges(2,2);
t127 = t160 * rSges(3,1) - t158 * rSges(3,2) + t163;
t126 = -t158 * rSges(3,1) - t160 * rSges(3,2) - t235;
t101 = t199 * t160;
t100 = t199 * t158;
t88 = t157 * t171 * t112;
t83 = t158 * pkin(7) + t163 + (pkin(2) - t228) * t160 + t206;
t82 = t227 - t235 + t152 + (-pkin(2) - t191) * t158;
t81 = t176 + t192;
t80 = -t235 + (rSges(5,3) - t167) * t160 + (-t154 - t190) * t158;
t78 = rSges(6,3) * t220 - t189;
t71 = t193 * t160;
t70 = t193 * t158;
t68 = -pkin(5) * t210 + (-t157 * t233 + t159 * t234) * t158;
t67 = t160 * (-t160 * t228 + t206) + (t158 * t191 - t227) * t158;
t57 = t65 * t219;
t56 = -t159 * t106 - t157 * t208 + t88;
t55 = -t115 * t219 - t159 * t79;
t54 = t115 * t220 + t159 * t78;
t53 = t185 * t160;
t52 = t185 * t158;
t51 = t192 + t79 + t207;
t50 = -t235 - t211 + (-t236 - t154 + (-rSges(6,3) - pkin(8)) * t157) * t158 + t189;
t49 = -t159 * t66 - t219 * t99;
t46 = t192 + t243;
t45 = -t235 + (pkin(5) * t168 - t167) * t160 + (-t153 * t159 - t154 + (-rSges(7,3) + t174) * t157) * t158 + t188;
t42 = (-t158 * t79 + t160 * t78) * t157;
t41 = t160 * t176 + (-t160 * rSges(5,3) + t158 * t190) * t158 + t231;
t38 = -t220 * t66 + t57;
t33 = -t159 * t232 + t219 * t230;
t32 = t159 * t68 + t220 * t89 + t48;
t29 = t124 * t75 + t125 * t77 + t219 * t73;
t28 = t124 * t74 + t125 * t76 + t219 * t72;
t27 = t122 * t75 + t123 * t77 + t220 * t73;
t26 = t122 * t74 + t123 * t76 + t220 * t72;
t23 = t158 * t78 + t160 * t79 + t196;
t18 = t57 + (-t158 * t232 + t160 * t68) * t157;
t17 = t232 * t160 + (t65 + t68) * t158 + t196;
t15 = t29 * t158 - t28 * t160;
t14 = t27 * t158 - t26 * t160;
t8 = -t44 * t159 + (t158 * t28 + t160 * t29) * t157;
t7 = -t43 * t159 + (t158 * t26 + t160 * t27) * t157;
t1 = [t172 * (Icges(4,2) * t172 + t224) + t169 * (Icges(4,1) * t169 + t223) + Icges(2,3) + Icges(3,3) + t84 + t88 + (Icges(5,2) * t159 - t106 + t222 - t96) * t159 + (Icges(5,1) * t157 - t208 + t221 - t226) * t157 + m(7) * (t45 ^ 2 + t46 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(4) * (t82 ^ 2 + t83 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2) + m(2) * (t146 ^ 2 + t147 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t53 * t45 + t52 * t46) + m(6) * (t71 * t50 + t70 * t51) + m(5) * (t100 * t81 + t101 * t80) + (-t30 / 0.2e1 - t39 / 0.2e1 - t82 * t238 - t172 * (-Icges(4,6) * t160 + t158 * t182) / 0.2e1 - t169 * (-Icges(4,5) * t160 + t158 * t184) / 0.2e1 - t157 * (-Icges(5,5) * t160 + t158 * t183) / 0.2e1 + (-Icges(5,6) * t160 + t158 * t181) * t240 + t200 * t160 - t204) * t160 + (t31 / 0.2e1 + t40 / 0.2e1 - t83 * t238 + (Icges(4,6) * t158 + t160 * t182) * t246 + (Icges(4,5) * t158 + t160 * t184) * t247 + (Icges(5,5) * t158 + t160 * t183) * t249 + (Icges(5,6) * t158 + t160 * t181) * t248 + t200 * t158 + t203) * t158; m(4) * t67 + m(5) * t41 + m(6) * t23 + m(7) * t17; m(7) * (t17 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t23 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t41 ^ 2) + m(4) * (t145 ^ 2 * t205 + t67 ^ 2) + (t156 * t245 - t12 - t14) * t160 + (t13 + t15 + t244 * t155 + (t158 * t245 + t244 * t160) * t160) * t158; m(7) * (t158 * t45 - t160 * t46) + m(6) * (t158 * t50 - t160 * t51) + m(5) * (t158 * t80 - t160 * t81); 0; m(7) * (t158 * t53 - t160 * t52) + m(6) * (t158 * t71 - t160 * t70) + m(5) * (-t160 * t100 + t158 * t101); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t205; (-t47 - t56) * t159 + m(7) * (t32 * t45 + t33 * t46) + m(6) * (t54 * t50 + t55 * t51) + (t158 * t204 + t160 * t203) * t157 + t197; m(6) * t42 + m(7) * t18; t8 * t241 + t7 * t239 + (t35 * t158 - t34 * t160) * t240 + (t160 * t15 / 0.2e1 + t14 * t241) * t157 + m(7) * (t17 * t18 + t32 * t53 + t33 * t52) + m(6) * (t23 * t42 + t54 * t71 + t55 * t70) + t194; m(6) * (t158 * t54 - t160 * t55) + m(7) * (t158 * t32 - t160 * t33); (t56 * t159 - t9) * t159 + m(7) * (t18 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t42 ^ 2 + t54 ^ 2 + t55 ^ 2) + (t160 * t8 + t158 * t7 - t159 * (t158 * t34 + t160 * t35)) * t157 + t242; -t225 + m(7) * (t45 * t48 + t46 * t49) + t197; m(7) * t38; m(7) * (t17 * t38 + t48 * t53 + t49 * t52) + t194; m(7) * (t158 * t48 - t160 * t49); m(7) * (t18 * t38 + t32 * t48 + t33 * t49) + t195; m(7) * (t38 ^ 2 + t48 ^ 2 + t49 ^ 2) + t195;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
