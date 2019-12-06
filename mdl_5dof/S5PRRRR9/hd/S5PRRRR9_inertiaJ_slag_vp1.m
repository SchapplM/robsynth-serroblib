% Calculate joint inertia matrix for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:04
% EndTime: 2019-12-05 17:19:15
% DurationCPUTime: 3.04s
% Computational Cost: add. (17185->452), mult. (37274->661), div. (0->0), fcn. (47850->12), ass. (0->215)
t195 = sin(pkin(10));
t197 = cos(pkin(10));
t203 = cos(qJ(2));
t198 = cos(pkin(5));
t201 = sin(qJ(2));
t225 = t198 * t201;
t182 = t195 * t203 + t197 * t225;
t200 = sin(qJ(3));
t196 = sin(pkin(5));
t236 = cos(qJ(3));
t212 = t196 * t236;
t171 = t182 * t200 + t197 * t212;
t242 = t171 / 0.2e1;
t184 = -t195 * t225 + t197 * t203;
t173 = t184 * t200 - t195 * t212;
t241 = t173 / 0.2e1;
t224 = t203 * t198;
t181 = t195 * t201 - t197 * t224;
t240 = t181 / 0.2e1;
t183 = t195 * t224 + t197 * t201;
t239 = t183 / 0.2e1;
t228 = t196 * t200;
t185 = -t198 * t236 + t201 * t228;
t238 = t185 / 0.2e1;
t237 = t198 / 0.2e1;
t202 = cos(qJ(4));
t235 = t202 * pkin(4);
t226 = t197 * t196;
t172 = t182 * t236 - t200 * t226;
t194 = qJ(4) + qJ(5);
t192 = sin(t194);
t193 = cos(t194);
t142 = -t172 * t192 + t181 * t193;
t143 = t172 * t193 + t181 * t192;
t102 = t143 * rSges(6,1) + t142 * rSges(6,2) + t171 * rSges(6,3);
t199 = sin(qJ(4));
t231 = t181 * t199;
t94 = pkin(4) * t231 + pkin(9) * t171 + t172 * t235;
t233 = t102 + t94;
t174 = t184 * t236 + t195 * t228;
t144 = -t174 * t192 + t183 * t193;
t145 = t174 * t193 + t183 * t192;
t103 = t145 * rSges(6,1) + t144 * rSges(6,2) + t173 * rSges(6,3);
t230 = t183 * t199;
t95 = pkin(4) * t230 + pkin(9) * t173 + t174 * t235;
t232 = t103 + t95;
t229 = t195 * t196;
t227 = t196 * t203;
t148 = -t174 * t199 + t183 * t202;
t149 = t174 * t202 + t230;
t113 = t149 * rSges(5,1) + t148 * rSges(5,2) + t173 * rSges(5,3);
t141 = t174 * pkin(3) + t173 * pkin(8);
t223 = -t113 - t141;
t186 = t198 * t200 + t201 * t212;
t169 = -t186 * t192 - t193 * t227;
t170 = t186 * t193 - t192 * t227;
t119 = t170 * rSges(6,1) + t169 * rSges(6,2) + t185 * rSges(6,3);
t214 = t199 * t227;
t120 = -pkin(4) * t214 + pkin(9) * t185 + t186 * t235;
t222 = t119 + t120;
t175 = -t186 * t199 - t202 * t227;
t176 = t186 * t202 - t214;
t133 = t176 * rSges(5,1) + t175 * rSges(5,2) + t185 * rSges(5,3);
t168 = t186 * pkin(3) + t185 * pkin(8);
t221 = -t133 - t168;
t140 = t172 * pkin(3) + t171 * pkin(8);
t220 = t140 * t227 + t181 * t168;
t167 = t184 * pkin(2) + t183 * pkin(7);
t165 = t198 * t167;
t219 = t198 * t141 + t165;
t166 = t182 * pkin(2) + t181 * pkin(7);
t218 = -t140 - t166;
t217 = t166 * t229 + t167 * t226;
t100 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t171;
t96 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t171;
t98 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t171;
t57 = t170 * t100 + t169 * t98 + t185 * t96;
t101 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t173;
t97 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t173;
t99 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t173;
t58 = t170 * t101 + t169 * t99 + t185 * t97;
t116 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t185;
t117 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t185;
t118 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t185;
t72 = t185 * t116 + t169 * t117 + t170 * t118;
t26 = t57 * t171 + t58 * t173 + t72 * t185;
t47 = t143 * t100 + t142 * t98 + t171 * t96;
t48 = t143 * t101 + t142 * t99 + t171 * t97;
t63 = t171 * t116 + t142 * t117 + t143 * t118;
t7 = t47 * t171 + t48 * t173 + t63 * t185;
t49 = t145 * t100 + t144 * t98 + t173 * t96;
t50 = t145 * t101 + t144 * t99 + t173 * t97;
t64 = t173 * t116 + t144 * t117 + t145 * t118;
t8 = t49 * t171 + t50 * t173 + t64 * t185;
t216 = t171 * t7 + t173 * t8 + t185 * t26;
t215 = -t141 - t232;
t213 = -t168 - t222;
t211 = -t227 / 0.2e1;
t162 = t186 * rSges(4,1) - t185 * rSges(4,2) - rSges(4,3) * t227;
t187 = (pkin(2) * t201 - pkin(7) * t203) * t196;
t210 = (-t162 - t187) * t196;
t209 = t140 * t229 + t141 * t226 + t217;
t208 = (-t187 + t221) * t196;
t13 = t47 * t181 + t48 * t183 - t227 * t63;
t14 = t49 * t181 + t50 * t183 - t227 * t64;
t28 = t57 * t181 + t58 * t183 - t227 * t72;
t207 = t13 * t242 + t14 * t241 + t26 * t211 + t28 * t238 + t8 * t239 + t7 * t240;
t15 = t63 * t198 + (t195 * t48 - t197 * t47) * t196;
t16 = t64 * t198 + (t195 * t50 - t197 * t49) * t196;
t30 = t72 * t198 + (t195 * t58 - t197 * t57) * t196;
t206 = t15 * t242 + t16 * t241 + t26 * t237 + t30 * t238 + t8 * t229 / 0.2e1 - t7 * t226 / 0.2e1;
t205 = (-t187 + t213) * t196;
t180 = t198 * rSges(3,3) + (rSges(3,1) * t201 + rSges(3,2) * t203) * t196;
t179 = Icges(3,5) * t198 + (Icges(3,1) * t201 + Icges(3,4) * t203) * t196;
t178 = Icges(3,6) * t198 + (Icges(3,4) * t201 + Icges(3,2) * t203) * t196;
t177 = Icges(3,3) * t198 + (Icges(3,5) * t201 + Icges(3,6) * t203) * t196;
t161 = Icges(4,1) * t186 - Icges(4,4) * t185 - Icges(4,5) * t227;
t160 = Icges(4,4) * t186 - Icges(4,2) * t185 - Icges(4,6) * t227;
t159 = Icges(4,5) * t186 - Icges(4,6) * t185 - Icges(4,3) * t227;
t158 = t184 * rSges(3,1) - t183 * rSges(3,2) + rSges(3,3) * t229;
t157 = t182 * rSges(3,1) - t181 * rSges(3,2) - rSges(3,3) * t226;
t156 = Icges(3,1) * t184 - Icges(3,4) * t183 + Icges(3,5) * t229;
t155 = Icges(3,1) * t182 - Icges(3,4) * t181 - Icges(3,5) * t226;
t154 = Icges(3,4) * t184 - Icges(3,2) * t183 + Icges(3,6) * t229;
t153 = Icges(3,4) * t182 - Icges(3,2) * t181 - Icges(3,6) * t226;
t152 = Icges(3,5) * t184 - Icges(3,6) * t183 + Icges(3,3) * t229;
t151 = Icges(3,5) * t182 - Icges(3,6) * t181 - Icges(3,3) * t226;
t147 = t172 * t202 + t231;
t146 = -t172 * t199 + t181 * t202;
t135 = -t198 * t157 - t180 * t226;
t134 = t198 * t158 - t180 * t229;
t132 = t183 * t140;
t131 = Icges(5,1) * t176 + Icges(5,4) * t175 + Icges(5,5) * t185;
t130 = Icges(5,4) * t176 + Icges(5,2) * t175 + Icges(5,6) * t185;
t129 = Icges(5,5) * t176 + Icges(5,6) * t175 + Icges(5,3) * t185;
t128 = t174 * rSges(4,1) - t173 * rSges(4,2) + t183 * rSges(4,3);
t127 = t172 * rSges(4,1) - t171 * rSges(4,2) + t181 * rSges(4,3);
t126 = Icges(4,1) * t174 - Icges(4,4) * t173 + Icges(4,5) * t183;
t125 = Icges(4,1) * t172 - Icges(4,4) * t171 + Icges(4,5) * t181;
t124 = Icges(4,4) * t174 - Icges(4,2) * t173 + Icges(4,6) * t183;
t123 = Icges(4,4) * t172 - Icges(4,2) * t171 + Icges(4,6) * t181;
t122 = Icges(4,5) * t174 - Icges(4,6) * t173 + Icges(4,3) * t183;
t121 = Icges(4,5) * t172 - Icges(4,6) * t171 + Icges(4,3) * t181;
t115 = (t157 * t195 + t158 * t197) * t196;
t114 = t171 * t119;
t112 = t147 * rSges(5,1) + t146 * rSges(5,2) + t171 * rSges(5,3);
t111 = Icges(5,1) * t149 + Icges(5,4) * t148 + Icges(5,5) * t173;
t110 = Icges(5,1) * t147 + Icges(5,4) * t146 + Icges(5,5) * t171;
t109 = Icges(5,4) * t149 + Icges(5,2) * t148 + Icges(5,6) * t173;
t108 = Icges(5,4) * t147 + Icges(5,2) * t146 + Icges(5,6) * t171;
t107 = Icges(5,5) * t149 + Icges(5,6) * t148 + Icges(5,3) * t173;
t106 = Icges(5,5) * t147 + Icges(5,6) * t146 + Icges(5,3) * t171;
t105 = -t128 * t227 - t183 * t162;
t104 = t127 * t227 + t181 * t162;
t93 = t185 * t103;
t92 = -t159 * t227 - t185 * t160 + t186 * t161;
t91 = t173 * t102;
t90 = t183 * t127 - t181 * t128;
t89 = (-t127 - t166) * t198 + t197 * t210;
t88 = t198 * t128 + t195 * t210 + t165;
t87 = t183 * t159 - t173 * t160 + t174 * t161;
t86 = t181 * t159 - t171 * t160 + t172 * t161;
t85 = (t127 * t195 + t128 * t197) * t196 + t217;
t84 = t185 * t113 - t173 * t133;
t83 = -t185 * t112 + t171 * t133;
t82 = -t173 * t119 + t93;
t81 = -t185 * t102 + t114;
t80 = -t122 * t227 - t185 * t124 + t186 * t126;
t79 = -t121 * t227 - t185 * t123 + t186 * t125;
t78 = t185 * t129 + t175 * t130 + t176 * t131;
t77 = t183 * t122 - t173 * t124 + t174 * t126;
t76 = t183 * t121 - t173 * t123 + t174 * t125;
t75 = t181 * t122 - t171 * t124 + t172 * t126;
t74 = t181 * t121 - t171 * t123 + t172 * t125;
t73 = t173 * t112 - t171 * t113;
t71 = -t171 * t103 + t91;
t70 = t183 * t221 + t223 * t227;
t69 = t112 * t227 + t181 * t133 + t220;
t68 = t173 * t129 + t148 * t130 + t149 * t131;
t67 = t171 * t129 + t146 * t130 + t147 * t131;
t66 = (-t112 + t218) * t198 + t197 * t208;
t65 = t198 * t113 + t195 * t208 + t219;
t62 = t183 * t112 + t181 * t223 + t132;
t61 = t185 * t107 + t175 * t109 + t176 * t111;
t60 = t185 * t106 + t175 * t108 + t176 * t110;
t59 = (t112 * t195 + t113 * t197) * t196 + t209;
t56 = t173 * t107 + t148 * t109 + t149 * t111;
t55 = t173 * t106 + t148 * t108 + t149 * t110;
t54 = t171 * t107 + t146 * t109 + t147 * t111;
t53 = t171 * t106 + t146 * t108 + t147 * t110;
t52 = -t173 * t222 + t185 * t95 + t93;
t51 = t171 * t120 - t185 * t233 + t114;
t46 = t183 * t213 + t215 * t227;
t45 = t181 * t222 + t227 * t233 + t220;
t44 = (t218 - t233) * t198 + t197 * t205;
t43 = t195 * t205 + t198 * t232 + t219;
t42 = -t171 * t232 + t173 * t94 + t91;
t41 = t92 * t198 + (t195 * t80 - t197 * t79) * t196;
t40 = t79 * t181 + t80 * t183 - t227 * t92;
t39 = t181 * t215 + t183 * t233 + t132;
t38 = (t195 * t233 + t197 * t232) * t196 + t209;
t37 = t87 * t198 + (t195 * t77 - t197 * t76) * t196;
t36 = t86 * t198 + (t195 * t75 - t197 * t74) * t196;
t35 = t76 * t181 + t77 * t183 - t227 * t87;
t34 = t74 * t181 + t75 * t183 - t227 * t86;
t33 = t78 * t198 + (t195 * t61 - t197 * t60) * t196;
t32 = t60 * t181 + t61 * t183 - t227 * t78;
t31 = t60 * t171 + t61 * t173 + t78 * t185;
t22 = t68 * t198 + (t195 * t56 - t197 * t55) * t196;
t21 = t67 * t198 + (t195 * t54 - t197 * t53) * t196;
t20 = t55 * t181 + t56 * t183 - t227 * t68;
t19 = t53 * t181 + t54 * t183 - t227 * t67;
t18 = t55 * t171 + t56 * t173 + t68 * t185;
t17 = t53 * t171 + t54 * t173 + t67 * t185;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t115 + m(4) * t85 + m(5) * t59 + m(6) * t38; m(6) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t59 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t85 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(3) * (t115 ^ 2 + t134 ^ 2 + t135 ^ 2) + (t16 + t22 + t37 + (t152 * t229 - t183 * t154 + t184 * t156) * t229) * t229 + (-t15 - t21 - t36 + (-t151 * t226 - t181 * t153 + t182 * t155) * t226 + (-t151 * t229 + t152 * t226 + t183 * t153 + t181 * t154 - t184 * t155 - t182 * t156) * t229) * t226 + (-(-t177 * t226 - t181 * t178 + t182 * t179) * t226 + (t177 * t229 - t183 * t178 + t184 * t179) * t229 + t30 + t33 + t41 + ((t154 * t203 + t156 * t201) * t195 - (t153 * t203 + t155 * t201) * t197) * t196 ^ 2 + ((-t151 * t197 + t152 * t195 + t178 * t203 + t179 * t201) * t196 + t198 * t177) * t198) * t198; m(4) * t90 + m(5) * t62 + m(6) * t39; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1) * t198 + (t16 / 0.2e1 + t22 / 0.2e1 + t37 / 0.2e1) * t183 + (t15 / 0.2e1 + t21 / 0.2e1 + t36 / 0.2e1) * t181 + m(6) * (t38 * t39 + t43 * t46 + t44 * t45) + m(5) * (t59 * t62 + t65 * t70 + t66 * t69) + m(4) * (t104 * t89 + t105 * t88 + t90 * t85) + ((-t30 / 0.2e1 - t33 / 0.2e1 - t41 / 0.2e1) * t203 + (-t13 / 0.2e1 - t19 / 0.2e1 - t34 / 0.2e1) * t197 + (t14 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1) * t195) * t196; (-t28 - t32 - t40) * t227 + (t14 + t20 + t35) * t183 + (t13 + t19 + t34) * t181 + m(6) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t62 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2 + t90 ^ 2); m(5) * t73 + m(6) * t42; t33 * t238 + t22 * t241 + t21 * t242 + t31 * t237 + (t195 * t18 / 0.2e1 - t197 * t17 / 0.2e1) * t196 + m(6) * (t38 * t42 + t43 * t52 + t44 * t51) + m(5) * (t59 * t73 + t65 * t84 + t66 * t83) + t206; t18 * t239 + t17 * t240 + t19 * t242 + t20 * t241 + t32 * t238 + t31 * t211 + m(6) * (t39 * t42 + t45 * t51 + t46 * t52) + m(5) * (t62 * t73 + t69 * t83 + t70 * t84) + t207; t171 * t17 + t173 * t18 + t185 * t31 + m(6) * (t42 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t73 ^ 2 + t83 ^ 2 + t84 ^ 2) + t216; m(6) * t71; m(6) * (t38 * t71 + t43 * t82 + t44 * t81) + t206; m(6) * (t39 * t71 + t45 * t81 + t46 * t82) + t207; m(6) * (t42 * t71 + t51 * t81 + t52 * t82) + t216; m(6) * (t71 ^ 2 + t81 ^ 2 + t82 ^ 2) + t216;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
