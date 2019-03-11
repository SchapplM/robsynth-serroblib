% Calculate joint inertia matrix for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:44
% EndTime: 2019-03-08 19:17:53
% DurationCPUTime: 5.00s
% Computational Cost: add. (13743->433), mult. (35218->654), div. (0->0), fcn. (45884->12), ass. (0->192)
t242 = Icges(4,1) + Icges(5,2);
t241 = Icges(4,4) + Icges(5,6);
t238 = Icges(5,4) - Icges(4,5);
t237 = Icges(5,5) - Icges(4,6);
t240 = Icges(4,2) + Icges(5,3);
t239 = Icges(5,1) + Icges(4,3);
t176 = sin(pkin(10));
t178 = cos(pkin(10));
t182 = sin(qJ(2));
t211 = sin(pkin(11));
t212 = cos(pkin(11));
t218 = cos(qJ(2));
t171 = -t182 * t211 + t218 * t212;
t179 = cos(pkin(6));
t185 = t179 * t171;
t186 = t182 * t212 + t211 * t218;
t146 = -t176 * t186 + t178 * t185;
t184 = t179 * t186;
t148 = t176 * t171 + t178 * t184;
t177 = sin(pkin(6));
t209 = t177 * t178;
t234 = -t240 * t146 - t241 * t148 - t237 * t209;
t149 = -t176 * t185 - t178 * t186;
t151 = t171 * t178 - t176 * t184;
t210 = t176 * t177;
t236 = t240 * t149 + t241 * t151 - t237 * t210;
t230 = -t241 * t146 - t242 * t148 - t238 * t209;
t229 = t241 * t149 + t242 * t151 - t238 * t210;
t162 = t171 * t177;
t163 = t186 * t177;
t228 = t240 * t162 + t241 * t163 - t237 * t179;
t227 = -t241 * t162 - t242 * t163 + t238 * t179;
t193 = m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1;
t235 = 0.2e1 * t193;
t233 = t237 * t146 + t238 * t148 + t239 * t209;
t232 = -t237 * t149 - t238 * t151 + t239 * t210;
t197 = t179 * t218;
t166 = -t176 * t182 + t178 * t197;
t207 = t179 * t182;
t167 = t176 * t218 + t178 * t207;
t134 = Icges(3,5) * t167 + Icges(3,6) * t166 - Icges(3,3) * t209;
t226 = -t134 + t233;
t168 = -t176 * t197 - t178 * t182;
t169 = -t176 * t207 + t178 * t218;
t135 = Icges(3,5) * t169 + Icges(3,6) * t168 + Icges(3,3) * t210;
t225 = t135 + t232;
t224 = (Icges(3,5) * t182 + Icges(3,6) * t218) * t177 - t238 * t163 - t237 * t162 + (Icges(3,3) + t239) * t179;
t181 = sin(qJ(5));
t208 = t177 * t181;
t217 = cos(qJ(5));
t117 = t149 * t217 + t176 * t208;
t119 = -t146 * t217 + t178 * t208;
t152 = t162 * t217 + t179 * t181;
t198 = t177 * t217;
t118 = -t149 * t181 + t176 * t198;
t180 = sin(qJ(6));
t183 = cos(qJ(6));
t80 = -t118 * t180 + t151 * t183;
t81 = t118 * t183 + t151 * t180;
t53 = Icges(7,5) * t81 + Icges(7,6) * t80 + Icges(7,3) * t117;
t55 = Icges(7,4) * t81 + Icges(7,2) * t80 + Icges(7,6) * t117;
t57 = Icges(7,1) * t81 + Icges(7,4) * t80 + Icges(7,5) * t117;
t17 = t117 * t53 + t55 * t80 + t57 * t81;
t120 = -t146 * t181 - t178 * t198;
t82 = -t120 * t180 + t148 * t183;
t83 = t120 * t183 + t148 * t180;
t54 = Icges(7,5) * t83 + Icges(7,6) * t82 - Icges(7,3) * t119;
t56 = Icges(7,4) * t83 + Icges(7,2) * t82 - Icges(7,6) * t119;
t58 = Icges(7,1) * t83 + Icges(7,4) * t82 - Icges(7,5) * t119;
t18 = t117 * t54 + t56 * t80 + t58 * t81;
t153 = -t162 * t181 + t179 * t217;
t114 = -t153 * t180 + t163 * t183;
t115 = t153 * t183 + t163 * t180;
t71 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t152;
t72 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t152;
t73 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t152;
t29 = t117 * t71 + t72 * t80 + t73 * t81;
t1 = t117 * t17 - t119 * t18 + t152 * t29;
t223 = t1 / 0.2e1;
t21 = t114 * t55 + t115 * t57 + t152 * t53;
t22 = t114 * t56 + t115 * t58 + t152 * t54;
t36 = t114 * t72 + t115 * t73 + t152 * t71;
t7 = t117 * t21 - t119 * t22 + t152 * t36;
t222 = t7 / 0.2e1;
t221 = t117 / 0.2e1;
t220 = -t119 / 0.2e1;
t219 = t152 / 0.2e1;
t216 = pkin(2) * t218;
t59 = rSges(7,1) * t81 + rSges(7,2) * t80 + rSges(7,3) * t117;
t215 = pkin(5) * t118 + pkin(9) * t117 + t59;
t60 = rSges(7,1) * t83 + rSges(7,2) * t82 - rSges(7,3) * t119;
t214 = pkin(5) * t120 - pkin(9) * t119 + t60;
t74 = rSges(7,1) * t115 + rSges(7,2) * t114 + rSges(7,3) * t152;
t213 = pkin(5) * t153 + pkin(9) * t152 + t74;
t104 = pkin(3) * t151 - qJ(4) * t149;
t196 = pkin(2) * t207 - qJ(3) * t177;
t144 = -t176 * t196 + t178 * t216;
t133 = t179 * t144;
t206 = t179 * t104 + t133;
t103 = pkin(3) * t148 - qJ(4) * t146;
t143 = t176 * t216 + t178 * t196;
t205 = -t103 - t143;
t204 = t143 * t210 + t144 * t209;
t172 = pkin(2) * t177 * t182 + qJ(3) * t179;
t203 = -pkin(3) * t163 + qJ(4) * t162 - t172;
t202 = m(4) + m(5) + m(6) + m(7);
t121 = pkin(4) * t210 + pkin(8) * t151;
t201 = t179 * t121 + t206;
t122 = -pkin(4) * t209 + pkin(8) * t148;
t200 = -t122 + t205;
t199 = -pkin(4) * t179 - pkin(8) * t163 + t203;
t194 = (-rSges(4,1) * t163 - rSges(4,2) * t162 - rSges(4,3) * t179 - t172) * t177;
t192 = t103 * t210 + t104 * t209 + t204;
t191 = (-rSges(5,1) * t179 + rSges(5,2) * t163 + rSges(5,3) * t162 + t203) * t177;
t108 = rSges(6,1) * t153 - rSges(6,2) * t152 + rSges(6,3) * t163;
t189 = (-t108 + t199) * t177;
t188 = t121 * t209 + t122 * t210 + t192;
t187 = (t199 - t213) * t177;
t161 = t179 * rSges(3,3) + (rSges(3,1) * t182 + rSges(3,2) * t218) * t177;
t160 = Icges(3,5) * t179 + (Icges(3,1) * t182 + Icges(3,4) * t218) * t177;
t159 = Icges(3,6) * t179 + (Icges(3,4) * t182 + Icges(3,2) * t218) * t177;
t141 = rSges(3,1) * t169 + rSges(3,2) * t168 + rSges(3,3) * t210;
t140 = rSges(3,1) * t167 + rSges(3,2) * t166 - rSges(3,3) * t209;
t139 = Icges(3,1) * t169 + Icges(3,4) * t168 + Icges(3,5) * t210;
t138 = Icges(3,1) * t167 + Icges(3,4) * t166 - Icges(3,5) * t209;
t137 = Icges(3,4) * t169 + Icges(3,2) * t168 + Icges(3,6) * t210;
t136 = Icges(3,4) * t167 + Icges(3,2) * t166 - Icges(3,6) * t209;
t110 = -t140 * t179 - t161 * t209;
t109 = t141 * t179 - t161 * t210;
t107 = Icges(6,1) * t153 - Icges(6,4) * t152 + Icges(6,5) * t163;
t106 = Icges(6,4) * t153 - Icges(6,2) * t152 + Icges(6,6) * t163;
t105 = Icges(6,5) * t153 - Icges(6,6) * t152 + Icges(6,3) * t163;
t99 = rSges(4,1) * t151 + rSges(4,2) * t149 + rSges(4,3) * t210;
t98 = rSges(4,1) * t148 + rSges(4,2) * t146 - rSges(4,3) * t209;
t97 = -rSges(5,1) * t209 - rSges(5,2) * t148 - rSges(5,3) * t146;
t96 = rSges(5,1) * t210 - rSges(5,2) * t151 - rSges(5,3) * t149;
t79 = (t140 * t176 + t141 * t178) * t177;
t70 = rSges(6,1) * t120 + rSges(6,2) * t119 + rSges(6,3) * t148;
t69 = rSges(6,1) * t118 - rSges(6,2) * t117 + rSges(6,3) * t151;
t68 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t148;
t67 = Icges(6,1) * t118 - Icges(6,4) * t117 + Icges(6,5) * t151;
t66 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t148;
t65 = Icges(6,4) * t118 - Icges(6,2) * t117 + Icges(6,6) * t151;
t64 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t148;
t63 = Icges(6,5) * t118 - Icges(6,6) * t117 + Icges(6,3) * t151;
t62 = (-t143 - t98) * t179 + t178 * t194;
t61 = t176 * t194 + t179 * t99 + t133;
t52 = (t176 * t98 + t178 * t99) * t177 + t204;
t51 = -t108 * t151 + t163 * t69;
t50 = t108 * t148 - t163 * t70;
t49 = t105 * t163 - t106 * t152 + t107 * t153;
t48 = (-t97 + t205) * t179 + t178 * t191;
t47 = t176 * t191 + t179 * t96 + t206;
t46 = -t148 * t69 + t151 * t70;
t45 = t105 * t148 + t106 * t119 + t107 * t120;
t44 = t105 * t151 - t106 * t117 + t107 * t118;
t43 = (t176 * t97 + t178 * t96) * t177 + t192;
t42 = -t119 * t74 - t152 * t60;
t41 = -t117 * t74 + t152 * t59;
t40 = (-t70 + t200) * t179 + t178 * t189;
t39 = t176 * t189 + t179 * t69 + t201;
t38 = -t152 * t66 + t153 * t68 + t163 * t64;
t37 = -t152 * t65 + t153 * t67 + t163 * t63;
t35 = t119 * t66 + t120 * t68 + t148 * t64;
t34 = t119 * t65 + t120 * t67 + t148 * t63;
t33 = -t117 * t66 + t118 * t68 + t151 * t64;
t32 = -t117 * t65 + t118 * t67 + t151 * t63;
t31 = t117 * t60 + t119 * t59;
t30 = -t119 * t71 + t72 * t82 + t73 * t83;
t28 = -t151 * t213 + t163 * t215;
t27 = t148 * t213 - t163 * t214;
t26 = (t176 * t70 + t178 * t69) * t177 + t188;
t25 = (t200 - t214) * t179 + t178 * t187;
t24 = t176 * t187 + t179 * t215 + t201;
t23 = -t148 * t215 + t151 * t214;
t20 = -t119 * t54 + t56 * t82 + t58 * t83;
t19 = -t119 * t53 + t55 * t82 + t57 * t83;
t16 = (t176 * t214 + t178 * t215) * t177 + t188;
t15 = t179 * t49 + (t176 * t37 - t178 * t38) * t177;
t14 = t148 * t38 + t151 * t37 + t163 * t49;
t13 = t179 * t45 + (t176 * t34 - t178 * t35) * t177;
t12 = t179 * t44 + (t176 * t32 - t178 * t33) * t177;
t11 = t148 * t35 + t151 * t34 + t163 * t45;
t10 = t148 * t33 + t151 * t32 + t163 * t44;
t9 = t179 * t36 + (t176 * t21 - t178 * t22) * t177;
t8 = t148 * t22 + t151 * t21 + t163 * t36;
t6 = t179 * t30 + (t176 * t19 - t178 * t20) * t177;
t5 = t179 * t29 + (t17 * t176 - t178 * t18) * t177;
t4 = t148 * t20 + t151 * t19 + t163 * t30;
t3 = t148 * t18 + t151 * t17 + t163 * t29;
t2 = t117 * t19 - t119 * t20 + t152 * t30;
t75 = [m(2) + m(3) + t202; m(3) * t79 + m(4) * t52 + m(5) * t43 + m(6) * t26 + m(7) * t16; m(7) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t26 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t43 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t52 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(3) * (t109 ^ 2 + t110 ^ 2 + t79 ^ 2) + (t9 + t15 + (t228 * t162 - t227 * t163 + t224 * t179) * t179 + ((t159 * t218 + t160 * t182) * t179 + (t162 * t234 + t163 * t230 + t179 * t233) * t178 + (t162 * t236 + t163 * t229 + t179 * t232) * t176) * t177) * t179 + (t5 + t12 + (t137 * t168 + t139 * t169 + t236 * t149 + t229 * t151 + t225 * t210) * t210 + (t159 * t168 + t160 * t169 + (t137 * t218 + t139 * t182) * t177 + t224 * t210 + t179 * t135 - t227 * t151 + t228 * t149) * t179) * t210 + (-t6 - t13 + (t136 * t166 + t138 * t167 - t234 * t146 - t230 * t148 + t226 * t209) * t209 + (-t159 * t166 - t160 * t167 - (t136 * t218 + t138 * t182) * t177 + t224 * t209 - t179 * t134 + t227 * t148 - t228 * t146) * t179 + (-t136 * t168 - t137 * t166 - t138 * t169 - t139 * t167 - t146 * t236 - t148 * t229 + t149 * t234 + t151 * t230 + t209 * t225 + t210 * t226) * t210) * t209; t202 * t179; m(7) * (t16 * t179 + (t176 * t25 - t178 * t24) * t177) + m(6) * (t179 * t26 + (t176 * t40 - t178 * t39) * t177) + m(5) * (t179 * t43 + (t176 * t48 - t178 * t47) * t177) + m(4) * (t179 * t52 + (t176 * t62 - t178 * t61) * t177); 0.2e1 * (m(4) / 0.2e1 + t193) * (t179 ^ 2 + (t176 ^ 2 + t178 ^ 2) * t177 ^ 2); -t162 * t235; m(7) * (-t146 * t24 - t149 * t25 - t16 * t162) + m(6) * (-t146 * t39 - t149 * t40 - t162 * t26) + m(5) * (-t146 * t47 - t149 * t48 - t162 * t43); (-t179 * t162 + (t178 * t146 - t149 * t176) * t177) * t235; (t146 ^ 2 + t149 ^ 2 + t162 ^ 2) * t235; m(6) * t46 + m(7) * t23; (t8 / 0.2e1 + t14 / 0.2e1) * t179 + (t9 / 0.2e1 + t15 / 0.2e1) * t163 + (t5 / 0.2e1 + t12 / 0.2e1) * t151 + (t6 / 0.2e1 + t13 / 0.2e1) * t148 + m(7) * (t16 * t23 + t24 * t28 + t25 * t27) + m(6) * (t26 * t46 + t39 * t51 + t40 * t50) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t178 + (t3 / 0.2e1 + t10 / 0.2e1) * t176) * t177; m(6) * (t179 * t46 + (t176 * t50 - t178 * t51) * t177) + m(7) * (t179 * t23 + (t176 * t27 - t178 * t28) * t177); m(6) * (-t146 * t51 - t149 * t50 - t162 * t46) + m(7) * (-t146 * t28 - t149 * t27 - t162 * t23); (t8 + t14) * t163 + (t3 + t10) * t151 + (t4 + t11) * t148 + m(7) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t46 ^ 2 + t50 ^ 2 + t51 ^ 2); m(7) * t31; m(7) * (t16 * t31 + t24 * t41 + t25 * t42) + t6 * t220 + t5 * t221 + t9 * t219 + t179 * t222 + (t176 * t223 - t178 * t2 / 0.2e1) * t177; m(7) * (t179 * t31 + (t176 * t42 - t178 * t41) * t177); m(7) * (-t146 * t41 - t149 * t42 - t162 * t31); t148 * t2 / 0.2e1 + t4 * t220 + t3 * t221 + t163 * t222 + t151 * t223 + m(7) * (t23 * t31 + t27 * t42 + t28 * t41) + t8 * t219; m(7) * (t31 ^ 2 + t41 ^ 2 + t42 ^ 2) + t117 * t1 - t119 * t2 + t152 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t75(1) t75(2) t75(4) t75(7) t75(11) t75(16); t75(2) t75(3) t75(5) t75(8) t75(12) t75(17); t75(4) t75(5) t75(6) t75(9) t75(13) t75(18); t75(7) t75(8) t75(9) t75(10) t75(14) t75(19); t75(11) t75(12) t75(13) t75(14) t75(15) t75(20); t75(16) t75(17) t75(18) t75(19) t75(20) t75(21);];
Mq  = res;
