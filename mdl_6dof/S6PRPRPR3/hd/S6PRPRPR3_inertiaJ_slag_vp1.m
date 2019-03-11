% Calculate joint inertia matrix for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:32
% EndTime: 2019-03-08 19:34:41
% DurationCPUTime: 4.84s
% Computational Cost: add. (16383->508), mult. (42028->754), div. (0->0), fcn. (54905->12), ass. (0->218)
t221 = m(6) / 0.2e1 + m(7) / 0.2e1;
t246 = 0.2e1 * t221;
t196 = sin(pkin(10));
t198 = cos(pkin(10));
t229 = sin(pkin(11));
t230 = cos(pkin(11));
t235 = sin(qJ(2));
t237 = cos(qJ(2));
t188 = -t235 * t229 + t237 * t230;
t199 = cos(pkin(6));
t202 = t199 * t188;
t203 = t237 * t229 + t235 * t230;
t167 = -t196 * t203 + t198 * t202;
t180 = t203 * t199;
t168 = t180 * t198 + t188 * t196;
t197 = sin(pkin(6));
t227 = t197 * t198;
t118 = Icges(4,5) * t168 + Icges(4,6) * t167 - Icges(4,3) * t227;
t214 = t199 * t237;
t183 = -t196 * t235 + t198 * t214;
t213 = t199 * t235;
t184 = t196 * t237 + t198 * t213;
t156 = Icges(3,5) * t184 + Icges(3,6) * t183 - Icges(3,3) * t227;
t245 = -t156 - t118;
t169 = -t196 * t202 - t198 * t203;
t170 = -t180 * t196 + t188 * t198;
t228 = t196 * t197;
t119 = Icges(4,5) * t170 + Icges(4,6) * t169 + Icges(4,3) * t228;
t185 = -t196 * t214 - t198 * t235;
t186 = -t196 * t213 + t198 * t237;
t157 = Icges(3,5) * t186 + Icges(3,6) * t185 + Icges(3,3) * t228;
t244 = t157 + t119;
t178 = t188 * t197;
t179 = t203 * t197;
t243 = Icges(4,5) * t179 + Icges(4,6) * t178 + (t235 * Icges(3,5) + t237 * Icges(3,6)) * t197 + (Icges(4,3) + Icges(3,3)) * t199;
t234 = sin(qJ(4));
t215 = t197 * t234;
t236 = cos(qJ(4));
t145 = t168 * t236 - t198 * t215;
t147 = t170 * t236 + t196 * t215;
t172 = t179 * t236 + t199 * t234;
t216 = t197 * t236;
t144 = t168 * t234 + t198 * t216;
t200 = sin(qJ(6));
t201 = cos(qJ(6));
t112 = t144 * t201 + t167 * t200;
t113 = t144 * t200 - t167 * t201;
t71 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t145;
t73 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t145;
t75 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t145;
t24 = t112 * t73 + t113 * t75 + t145 * t71;
t146 = t170 * t234 - t196 * t216;
t114 = t146 * t201 + t169 * t200;
t115 = t146 * t200 - t169 * t201;
t72 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t147;
t74 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t147;
t76 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t147;
t25 = t112 * t74 + t113 * t76 + t145 * t72;
t171 = t179 * t234 - t199 * t236;
t142 = t171 * t201 + t178 * t200;
t143 = t171 * t200 - t178 * t201;
t100 = Icges(7,4) * t143 + Icges(7,2) * t142 + Icges(7,6) * t172;
t101 = Icges(7,1) * t143 + Icges(7,4) * t142 + Icges(7,5) * t172;
t99 = Icges(7,5) * t143 + Icges(7,6) * t142 + Icges(7,3) * t172;
t35 = t100 * t112 + t101 * t113 + t145 * t99;
t1 = t145 * t24 + t147 * t25 + t172 * t35;
t242 = -t1 / 0.2e1;
t30 = t142 * t73 + t143 * t75 + t172 * t71;
t31 = t142 * t74 + t143 * t76 + t172 * t72;
t50 = t100 * t142 + t101 * t143 + t172 * t99;
t7 = t145 * t30 + t147 * t31 + t172 * t50;
t241 = t7 / 0.2e1;
t240 = t145 / 0.2e1;
t239 = t147 / 0.2e1;
t238 = t172 / 0.2e1;
t233 = t237 * pkin(2);
t77 = rSges(7,1) * t113 + rSges(7,2) * t112 + rSges(7,3) * t145;
t232 = -pkin(5) * t167 + pkin(9) * t145 + t77;
t78 = rSges(7,1) * t115 + rSges(7,2) * t114 + rSges(7,3) * t147;
t231 = -pkin(5) * t169 + pkin(9) * t147 + t78;
t102 = rSges(7,1) * t143 + rSges(7,2) * t142 + rSges(7,3) * t172;
t226 = pkin(5) * t178 - pkin(9) * t172 - t102;
t130 = pkin(3) * t170 - pkin(8) * t169;
t212 = pkin(2) * t213 - qJ(3) * t197;
t165 = -t212 * t196 + t233 * t198;
t155 = t199 * t165;
t225 = t199 * t130 + t155;
t129 = pkin(3) * t168 - pkin(8) * t167;
t164 = t233 * t196 + t212 * t198;
t224 = -t129 - t164;
t223 = t164 * t228 + t165 * t227;
t189 = t197 * t235 * pkin(2) + t199 * qJ(3);
t222 = -pkin(3) * t179 + pkin(8) * t178 - t189;
t220 = m(4) + m(5) + m(6) + m(7);
t109 = pkin(4) * t147 + qJ(5) * t146;
t219 = t199 * t109 + t225;
t108 = pkin(4) * t145 + qJ(5) * t144;
t218 = -t108 + t224;
t141 = pkin(4) * t172 + qJ(5) * t171;
t217 = -t141 + t222;
t211 = (-rSges(4,1) * t179 - rSges(4,2) * t178 - rSges(4,3) * t199 - t189) * t197;
t210 = t129 * t228 + t130 * t227 + t223;
t138 = rSges(5,1) * t172 - rSges(5,2) * t171 - rSges(5,3) * t178;
t209 = (-t138 + t222) * t197;
t137 = -rSges(6,1) * t178 - rSges(6,2) * t172 + rSges(6,3) * t171;
t206 = (-t137 + t217) * t197;
t205 = t108 * t228 + t109 * t227 + t210;
t204 = (t217 + t226) * t197;
t177 = t199 * rSges(3,3) + (t235 * rSges(3,1) + t237 * rSges(3,2)) * t197;
t176 = Icges(3,5) * t199 + (t235 * Icges(3,1) + t237 * Icges(3,4)) * t197;
t175 = Icges(3,6) * t199 + (t235 * Icges(3,4) + t237 * Icges(3,2)) * t197;
t163 = rSges(3,1) * t186 + rSges(3,2) * t185 + rSges(3,3) * t228;
t162 = rSges(3,1) * t184 + rSges(3,2) * t183 - rSges(3,3) * t227;
t161 = Icges(3,1) * t186 + Icges(3,4) * t185 + Icges(3,5) * t228;
t160 = Icges(3,1) * t184 + Icges(3,4) * t183 - Icges(3,5) * t227;
t159 = Icges(3,4) * t186 + Icges(3,2) * t185 + Icges(3,6) * t228;
t158 = Icges(3,4) * t184 + Icges(3,2) * t183 - Icges(3,6) * t227;
t153 = Icges(4,1) * t179 + Icges(4,4) * t178 + Icges(4,5) * t199;
t152 = Icges(4,4) * t179 + Icges(4,2) * t178 + Icges(4,6) * t199;
t140 = -t162 * t199 - t177 * t227;
t139 = t163 * t199 - t177 * t228;
t136 = Icges(5,1) * t172 - Icges(5,4) * t171 - Icges(5,5) * t178;
t135 = -Icges(6,1) * t178 - Icges(6,4) * t172 + Icges(6,5) * t171;
t134 = Icges(5,4) * t172 - Icges(5,2) * t171 - Icges(5,6) * t178;
t133 = -Icges(6,4) * t178 - Icges(6,2) * t172 + Icges(6,6) * t171;
t132 = Icges(5,5) * t172 - Icges(5,6) * t171 - Icges(5,3) * t178;
t131 = -Icges(6,5) * t178 - Icges(6,6) * t172 + Icges(6,3) * t171;
t125 = rSges(4,1) * t170 + rSges(4,2) * t169 + rSges(4,3) * t228;
t124 = rSges(4,1) * t168 + rSges(4,2) * t167 - rSges(4,3) * t227;
t123 = Icges(4,1) * t170 + Icges(4,4) * t169 + Icges(4,5) * t228;
t122 = Icges(4,1) * t168 + Icges(4,4) * t167 - Icges(4,5) * t227;
t121 = Icges(4,4) * t170 + Icges(4,2) * t169 + Icges(4,6) * t228;
t120 = Icges(4,4) * t168 + Icges(4,2) * t167 - Icges(4,6) * t227;
t111 = (t162 * t196 + t163 * t198) * t197;
t110 = t167 * t141;
t104 = t178 * t109;
t98 = t169 * t108;
t97 = rSges(5,1) * t147 - rSges(5,2) * t146 - rSges(5,3) * t169;
t96 = rSges(5,1) * t145 - rSges(5,2) * t144 - rSges(5,3) * t167;
t95 = -rSges(6,1) * t169 - rSges(6,2) * t147 + rSges(6,3) * t146;
t94 = -rSges(6,1) * t167 - rSges(6,2) * t145 + rSges(6,3) * t144;
t93 = Icges(5,1) * t147 - Icges(5,4) * t146 - Icges(5,5) * t169;
t92 = Icges(5,1) * t145 - Icges(5,4) * t144 - Icges(5,5) * t167;
t91 = -Icges(6,1) * t169 - Icges(6,4) * t147 + Icges(6,5) * t146;
t90 = -Icges(6,1) * t167 - Icges(6,4) * t145 + Icges(6,5) * t144;
t89 = Icges(5,4) * t147 - Icges(5,2) * t146 - Icges(5,6) * t169;
t88 = Icges(5,4) * t145 - Icges(5,2) * t144 - Icges(5,6) * t167;
t87 = -Icges(6,4) * t169 - Icges(6,2) * t147 + Icges(6,6) * t146;
t86 = -Icges(6,4) * t167 - Icges(6,2) * t145 + Icges(6,6) * t144;
t85 = Icges(5,5) * t147 - Icges(5,6) * t146 - Icges(5,3) * t169;
t84 = Icges(5,5) * t145 - Icges(5,6) * t144 - Icges(5,3) * t167;
t83 = -Icges(6,5) * t169 - Icges(6,6) * t147 + Icges(6,3) * t146;
t82 = -Icges(6,5) * t167 - Icges(6,6) * t145 + Icges(6,3) * t144;
t80 = (-t124 - t164) * t199 + t198 * t211;
t79 = t125 * t199 + t196 * t211 + t155;
t70 = (t124 * t196 + t125 * t198) * t197 + t223;
t69 = t138 * t169 - t178 * t97;
t68 = -t138 * t167 + t178 * t96;
t67 = -t132 * t178 - t134 * t171 + t136 * t172;
t66 = t131 * t171 - t133 * t172 - t135 * t178;
t65 = t167 * t97 - t169 * t96;
t64 = -t132 * t169 - t134 * t146 + t136 * t147;
t63 = -t132 * t167 - t134 * t144 + t136 * t145;
t62 = t131 * t146 - t133 * t147 - t135 * t169;
t61 = t131 * t144 - t133 * t145 - t135 * t167;
t60 = (-t96 + t224) * t199 + t198 * t209;
t59 = t196 * t209 + t199 * t97 + t225;
t58 = -t102 * t147 + t172 * t78;
t57 = t102 * t145 - t172 * t77;
t56 = -t178 * t95 - t104 + (t137 + t141) * t169;
t55 = -t137 * t167 - t110 - (-t108 - t94) * t178;
t54 = -t171 * t89 + t172 * t93 - t178 * t85;
t53 = -t171 * t88 + t172 * t92 - t178 * t84;
t52 = t171 * t83 - t172 * t87 - t178 * t91;
t51 = t171 * t82 - t172 * t86 - t178 * t90;
t49 = (t196 * t96 + t198 * t97) * t197 + t210;
t48 = -t146 * t89 + t147 * t93 - t169 * t85;
t47 = -t146 * t88 + t147 * t92 - t169 * t84;
t46 = -t144 * t89 + t145 * t93 - t167 * t85;
t45 = -t144 * t88 + t145 * t92 - t167 * t84;
t44 = t146 * t83 - t147 * t87 - t169 * t91;
t43 = t146 * t82 - t147 * t86 - t169 * t90;
t42 = t144 * t83 - t145 * t87 - t167 * t91;
t41 = t144 * t82 - t145 * t86 - t167 * t90;
t40 = -t145 * t78 + t147 * t77;
t39 = (-t94 + t218) * t199 + t198 * t206;
t38 = t196 * t206 + t199 * t95 + t219;
t37 = -t169 * t94 - t98 + (t109 + t95) * t167;
t36 = t100 * t114 + t101 * t115 + t147 * t99;
t34 = (t196 * t94 + t198 * t95) * t197 + t205;
t33 = -t104 - t231 * t178 + (t141 - t226) * t169;
t32 = -t110 + t226 * t167 - (-t108 - t232) * t178;
t29 = (t218 - t232) * t199 + t198 * t204;
t28 = t196 * t204 + t199 * t231 + t219;
t27 = t114 * t74 + t115 * t76 + t147 * t72;
t26 = t114 * t73 + t115 * t75 + t147 * t71;
t23 = -t98 - t232 * t169 + (t109 + t231) * t167;
t22 = (t196 * t232 + t198 * t231) * t197 + t205;
t21 = t199 * t67 + (t196 * t54 - t198 * t53) * t197;
t20 = t199 * t66 + (t196 * t52 - t198 * t51) * t197;
t19 = -t167 * t53 - t169 * t54 - t178 * t67;
t18 = -t167 * t51 - t169 * t52 - t178 * t66;
t17 = t199 * t64 + (t196 * t48 - t198 * t47) * t197;
t16 = t199 * t63 + (t196 * t46 - t198 * t45) * t197;
t15 = t199 * t62 + (t196 * t44 - t198 * t43) * t197;
t14 = t199 * t61 + (t196 * t42 - t198 * t41) * t197;
t13 = -t167 * t47 - t169 * t48 - t178 * t64;
t12 = -t167 * t45 - t169 * t46 - t178 * t63;
t11 = -t167 * t43 - t169 * t44 - t178 * t62;
t10 = -t167 * t41 - t169 * t42 - t178 * t61;
t9 = t199 * t50 + (t196 * t31 - t198 * t30) * t197;
t8 = -t167 * t30 - t169 * t31 - t178 * t50;
t6 = t199 * t36 + (t196 * t27 - t198 * t26) * t197;
t5 = t199 * t35 + (t196 * t25 - t198 * t24) * t197;
t4 = -t167 * t26 - t169 * t27 - t178 * t36;
t3 = -t167 * t24 - t169 * t25 - t178 * t35;
t2 = t145 * t26 + t147 * t27 + t172 * t36;
t81 = [m(2) + m(3) + t220; m(3) * t111 + m(4) * t70 + m(5) * t49 + m(6) * t34 + m(7) * t22; m(7) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t34 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t49 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t70 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(3) * (t111 ^ 2 + t139 ^ 2 + t140 ^ 2) + (t9 + t21 + t20 + (t152 * t178 + t153 * t179 + t199 * t243) * t199 + ((t121 * t178 + t123 * t179) * t196 - (t120 * t178 + t122 * t179) * t198 + (-t118 * t198 + t119 * t196 + t175 * t237 + t176 * t235) * t199) * t197) * t199 + (t6 + t17 + t15 + (t121 * t169 + t123 * t170 + t159 * t185 + t161 * t186 + t228 * t244) * t228 + (t175 * t185 + t176 * t186 + t152 * t169 + t153 * t170 + (t159 * t237 + t161 * t235) * t197 + t243 * t228 + t199 * t157) * t199) * t228 + (-t5 - t16 - t14 + (t120 * t167 + t122 * t168 + t158 * t183 + t160 * t184 + t227 * t245) * t227 + (-t175 * t183 - t176 * t184 - t152 * t167 - t153 * t168 - (t158 * t237 + t160 * t235) * t197 + t243 * t227 - t199 * t156) * t199 + (-t120 * t169 - t121 * t167 - t122 * t170 - t123 * t168 - t158 * t185 - t159 * t183 - t160 * t186 - t161 * t184 + t244 * t227 + t228 * t245) * t228) * t227; t220 * t199; m(7) * (t199 * t22 + (t196 * t29 - t198 * t28) * t197) + m(6) * (t199 * t34 + (t196 * t39 - t198 * t38) * t197) + m(5) * (t199 * t49 + (t196 * t60 - t198 * t59) * t197) + m(4) * (t199 * t70 + (t196 * t80 - t198 * t79) * t197); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t221) * (t199 ^ 2 + (t196 ^ 2 + t198 ^ 2) * t197 ^ 2); m(5) * t65 + m(6) * t37 + m(7) * t23; (t8 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t199 - (t9 / 0.2e1 + t20 / 0.2e1 + t21 / 0.2e1) * t178 + (-t6 / 0.2e1 - t17 / 0.2e1 - t15 / 0.2e1) * t169 + (-t5 / 0.2e1 - t16 / 0.2e1 - t14 / 0.2e1) * t167 + m(7) * (t22 * t23 + t28 * t33 + t29 * t32) + m(6) * (t34 * t37 + t38 * t56 + t39 * t55) + m(5) * (t49 * t65 + t59 * t69 + t60 * t68) + ((-t3 / 0.2e1 - t12 / 0.2e1 - t10 / 0.2e1) * t198 + (t4 / 0.2e1 + t11 / 0.2e1 + t13 / 0.2e1) * t196) * t197; m(5) * (t199 * t65 + (t196 * t68 - t198 * t69) * t197) + m(6) * (t199 * t37 + (t196 * t55 - t198 * t56) * t197) + m(7) * (t199 * t23 + (t196 * t32 - t198 * t33) * t197); -(t8 + t18 + t19) * t178 + (-t4 - t13 - t11) * t169 + (-t3 - t10 - t12) * t167 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t37 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t65 ^ 2 + t68 ^ 2 + t69 ^ 2); t171 * t246; m(7) * (t144 * t28 + t146 * t29 + t171 * t22) + m(6) * (t144 * t38 + t146 * t39 + t171 * t34); (t171 * t199 + (-t144 * t198 + t146 * t196) * t197) * t246; m(7) * (t144 * t33 + t146 * t32 + t171 * t23) + m(6) * (t144 * t56 + t146 * t55 + t171 * t37); (t144 ^ 2 + t146 ^ 2 + t171 ^ 2) * t246; m(7) * t40; m(7) * (t22 * t40 + t28 * t58 + t29 * t57) + t5 * t240 + t6 * t239 + t199 * t241 + t9 * t238 + (t196 * t2 / 0.2e1 + t198 * t242) * t197; m(7) * (t199 * t40 + (t196 * t57 - t198 * t58) * t197); m(7) * (t23 * t40 + t32 * t57 + t33 * t58) + t8 * t238 + t3 * t240 - t169 * t2 / 0.2e1 + t167 * t242 + t4 * t239 - t178 * t241; m(7) * (t144 * t58 + t146 * t57 + t171 * t40); t147 * t2 + t145 * t1 + t172 * t7 + m(7) * (t40 ^ 2 + t57 ^ 2 + t58 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t81(1) t81(2) t81(4) t81(7) t81(11) t81(16); t81(2) t81(3) t81(5) t81(8) t81(12) t81(17); t81(4) t81(5) t81(6) t81(9) t81(13) t81(18); t81(7) t81(8) t81(9) t81(10) t81(14) t81(19); t81(11) t81(12) t81(13) t81(14) t81(15) t81(20); t81(16) t81(17) t81(18) t81(19) t81(20) t81(21);];
Mq  = res;
