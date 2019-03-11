% Calculate joint inertia matrix for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:12
% EndTime: 2019-03-08 21:08:22
% DurationCPUTime: 5.11s
% Computational Cost: add. (10688->530), mult. (27602->769), div. (0->0), fcn. (34651->10), ass. (0->230)
t258 = m(6) + m(7);
t256 = m(5) + t258;
t207 = sin(pkin(6));
t257 = t207 ^ 2;
t206 = sin(pkin(10));
t208 = cos(pkin(10));
t213 = cos(qJ(2));
t209 = cos(pkin(6));
t211 = sin(qJ(2));
t242 = t209 * t211;
t195 = t206 * t213 + t208 * t242;
t248 = sin(qJ(3));
t224 = t207 * t248;
t249 = cos(qJ(3));
t183 = t195 * t249 - t208 * t224;
t198 = -t206 * t242 + t208 * t213;
t185 = t198 * t249 + t206 * t224;
t225 = t207 * t249;
t200 = t209 * t248 + t211 * t225;
t184 = t198 * t248 - t206 * t225;
t241 = t209 * t213;
t196 = t206 * t241 + t208 * t211;
t210 = sin(qJ(6));
t212 = cos(qJ(6));
t148 = -t184 * t210 - t196 * t212;
t149 = t184 * t212 - t196 * t210;
t182 = t195 * t248 + t208 * t225;
t193 = t206 * t211 - t208 * t241;
t146 = -t182 * t210 - t193 * t212;
t147 = t182 * t212 - t193 * t210;
t91 = Icges(7,5) * t147 + Icges(7,6) * t146 + Icges(7,3) * t183;
t93 = Icges(7,4) * t147 + Icges(7,2) * t146 + Icges(7,6) * t183;
t95 = Icges(7,1) * t147 + Icges(7,4) * t146 + Icges(7,5) * t183;
t36 = t148 * t93 + t149 * t95 + t185 * t91;
t92 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t185;
t94 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t185;
t96 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t185;
t37 = t148 * t94 + t149 * t96 + t185 * t92;
t199 = -t209 * t249 + t211 * t224;
t243 = t207 * t213;
t186 = -t199 * t210 + t212 * t243;
t187 = t199 * t212 + t210 * t243;
t126 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t200;
t127 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t200;
t128 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t200;
t46 = t126 * t185 + t127 * t148 + t128 * t149;
t2 = t183 * t36 + t185 * t37 + t200 * t46;
t254 = t2 / 0.2e1;
t252 = t183 / 0.2e1;
t251 = t185 / 0.2e1;
t250 = t200 / 0.2e1;
t97 = rSges(7,1) * t147 + rSges(7,2) * t146 + rSges(7,3) * t183;
t247 = pkin(5) * t182 + pkin(9) * t183 + t97;
t98 = rSges(7,1) * t149 + rSges(7,2) * t148 + rSges(7,3) * t185;
t246 = pkin(5) * t184 + pkin(9) * t185 + t98;
t245 = t206 * t207;
t244 = t207 * t208;
t123 = rSges(5,1) * t185 + rSges(5,2) * t196 + rSges(5,3) * t184;
t139 = pkin(3) * t185 + qJ(4) * t184;
t240 = -t123 - t139;
t138 = pkin(3) * t183 + qJ(4) * t182;
t125 = t196 * t138;
t151 = pkin(4) * t183 - qJ(5) * t193;
t239 = t196 * t151 + t125;
t129 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t200;
t238 = pkin(5) * t199 + pkin(9) * t200 + t129;
t180 = pkin(3) * t200 + qJ(4) * t199;
t237 = t138 * t243 + t193 * t180;
t179 = pkin(2) * t198 + pkin(8) * t196;
t177 = t209 * t179;
t236 = t209 * t139 + t177;
t178 = pkin(2) * t195 + pkin(8) * t193;
t235 = -t138 - t178;
t152 = pkin(4) * t185 - qJ(5) * t196;
t234 = -t139 - t152;
t171 = t200 * rSges(5,1) - rSges(5,2) * t243 + t199 * rSges(5,3);
t233 = -t171 - t180;
t232 = t178 * t245 + t179 * t244;
t188 = t200 * pkin(4) + qJ(5) * t243;
t231 = -t180 - t188;
t122 = rSges(6,1) * t184 - rSges(6,2) * t185 - rSges(6,3) * t196;
t229 = -t122 + t234;
t228 = t209 * t152 + t236;
t227 = -t151 + t235;
t170 = t199 * rSges(6,1) - t200 * rSges(6,2) + rSges(6,3) * t243;
t226 = -t170 + t231;
t172 = t200 * rSges(4,1) - t199 * rSges(4,2) - rSges(4,3) * t243;
t201 = (pkin(2) * t211 - pkin(8) * t213) * t207;
t223 = (-t172 - t201) * t207;
t221 = t234 - t246;
t220 = t231 - t238;
t219 = t138 * t245 + t139 * t244 + t232;
t218 = t151 * t243 + t193 * t188 + t237;
t217 = (-t201 + t233) * t207;
t216 = (-t201 + t226) * t207;
t215 = t151 * t245 + t152 * t244 + t219;
t214 = (-t201 + t220) * t207;
t192 = t209 * rSges(3,3) + (rSges(3,1) * t211 + rSges(3,2) * t213) * t207;
t191 = Icges(3,5) * t209 + (Icges(3,1) * t211 + Icges(3,4) * t213) * t207;
t190 = Icges(3,6) * t209 + (Icges(3,4) * t211 + Icges(3,2) * t213) * t207;
t189 = Icges(3,3) * t209 + (Icges(3,5) * t211 + Icges(3,6) * t213) * t207;
t169 = Icges(4,1) * t200 - Icges(4,4) * t199 - Icges(4,5) * t243;
t168 = Icges(5,1) * t200 - Icges(5,4) * t243 + Icges(5,5) * t199;
t167 = Icges(6,1) * t199 - Icges(6,4) * t200 + Icges(6,5) * t243;
t166 = Icges(4,4) * t200 - Icges(4,2) * t199 - Icges(4,6) * t243;
t165 = Icges(5,4) * t200 - Icges(5,2) * t243 + Icges(5,6) * t199;
t164 = Icges(6,4) * t199 - Icges(6,2) * t200 + Icges(6,6) * t243;
t163 = Icges(4,5) * t200 - Icges(4,6) * t199 - Icges(4,3) * t243;
t162 = Icges(5,5) * t200 - Icges(5,6) * t243 + Icges(5,3) * t199;
t161 = Icges(6,5) * t199 - Icges(6,6) * t200 + Icges(6,3) * t243;
t160 = rSges(3,1) * t198 - rSges(3,2) * t196 + rSges(3,3) * t245;
t159 = rSges(3,1) * t195 - rSges(3,2) * t193 - rSges(3,3) * t244;
t158 = Icges(3,1) * t198 - Icges(3,4) * t196 + Icges(3,5) * t245;
t157 = Icges(3,1) * t195 - Icges(3,4) * t193 - Icges(3,5) * t244;
t156 = Icges(3,4) * t198 - Icges(3,2) * t196 + Icges(3,6) * t245;
t155 = Icges(3,4) * t195 - Icges(3,2) * t193 - Icges(3,6) * t244;
t154 = Icges(3,5) * t198 - Icges(3,6) * t196 + Icges(3,3) * t245;
t153 = Icges(3,5) * t195 - Icges(3,6) * t193 - Icges(3,3) * t244;
t132 = -t159 * t209 - t192 * t244;
t131 = t160 * t209 - t192 * t245;
t124 = rSges(4,1) * t185 - rSges(4,2) * t184 + rSges(4,3) * t196;
t121 = rSges(4,1) * t183 - rSges(4,2) * t182 + rSges(4,3) * t193;
t120 = rSges(5,1) * t183 + rSges(5,2) * t193 + rSges(5,3) * t182;
t119 = rSges(6,1) * t182 - rSges(6,2) * t183 - rSges(6,3) * t193;
t118 = Icges(4,1) * t185 - Icges(4,4) * t184 + Icges(4,5) * t196;
t117 = Icges(4,1) * t183 - Icges(4,4) * t182 + Icges(4,5) * t193;
t116 = Icges(5,1) * t185 + Icges(5,4) * t196 + Icges(5,5) * t184;
t115 = Icges(5,1) * t183 + Icges(5,4) * t193 + Icges(5,5) * t182;
t114 = Icges(6,1) * t184 - Icges(6,4) * t185 - Icges(6,5) * t196;
t113 = Icges(6,1) * t182 - Icges(6,4) * t183 - Icges(6,5) * t193;
t112 = Icges(4,4) * t185 - Icges(4,2) * t184 + Icges(4,6) * t196;
t111 = Icges(4,4) * t183 - Icges(4,2) * t182 + Icges(4,6) * t193;
t110 = Icges(5,4) * t185 + Icges(5,2) * t196 + Icges(5,6) * t184;
t109 = Icges(5,4) * t183 + Icges(5,2) * t193 + Icges(5,6) * t182;
t108 = Icges(6,4) * t184 - Icges(6,2) * t185 - Icges(6,6) * t196;
t107 = Icges(6,4) * t182 - Icges(6,2) * t183 - Icges(6,6) * t193;
t106 = Icges(4,5) * t185 - Icges(4,6) * t184 + Icges(4,3) * t196;
t105 = Icges(4,5) * t183 - Icges(4,6) * t182 + Icges(4,3) * t193;
t104 = Icges(5,5) * t185 + Icges(5,6) * t196 + Icges(5,3) * t184;
t103 = Icges(5,5) * t183 + Icges(5,6) * t193 + Icges(5,3) * t182;
t102 = Icges(6,5) * t184 - Icges(6,6) * t185 - Icges(6,3) * t196;
t101 = Icges(6,5) * t182 - Icges(6,6) * t183 - Icges(6,3) * t193;
t100 = (t159 * t206 + t160 * t208) * t207;
t90 = -t124 * t243 - t196 * t172;
t89 = t121 * t243 + t193 * t172;
t88 = -t163 * t243 - t199 * t166 + t200 * t169;
t87 = t199 * t162 - t165 * t243 + t200 * t168;
t86 = t161 * t243 - t200 * t164 + t199 * t167;
t85 = t121 * t196 - t124 * t193;
t84 = (-t121 - t178) * t209 + t208 * t223;
t83 = t124 * t209 + t206 * t223 + t177;
t82 = t163 * t196 - t166 * t184 + t169 * t185;
t81 = t162 * t184 + t165 * t196 + t168 * t185;
t80 = -t161 * t196 - t164 * t185 + t167 * t184;
t79 = t163 * t193 - t166 * t182 + t169 * t183;
t78 = t162 * t182 + t165 * t193 + t168 * t183;
t77 = -t161 * t193 - t164 * t183 + t167 * t182;
t76 = (t121 * t206 + t124 * t208) * t207 + t232;
t75 = -t129 * t185 + t200 * t98;
t74 = t129 * t183 - t200 * t97;
t73 = t233 * t196 + t240 * t243;
t72 = t120 * t243 + t193 * t171 + t237;
t71 = -t106 * t243 - t199 * t112 + t200 * t118;
t70 = -t105 * t243 - t199 * t111 + t200 * t117;
t69 = t199 * t104 - t110 * t243 + t200 * t116;
t68 = t199 * t103 - t109 * t243 + t200 * t115;
t67 = t102 * t243 - t200 * t108 + t199 * t114;
t66 = t101 * t243 - t200 * t107 + t199 * t113;
t65 = (-t120 + t235) * t209 + t208 * t217;
t64 = t123 * t209 + t206 * t217 + t236;
t63 = t126 * t200 + t127 * t186 + t128 * t187;
t62 = t106 * t196 - t112 * t184 + t118 * t185;
t61 = t105 * t196 - t111 * t184 + t117 * t185;
t60 = t104 * t184 + t110 * t196 + t116 * t185;
t59 = t103 * t184 + t109 * t196 + t115 * t185;
t58 = -t102 * t196 - t108 * t185 + t114 * t184;
t57 = -t101 * t196 - t107 * t185 + t113 * t184;
t56 = t106 * t193 - t112 * t182 + t118 * t183;
t55 = t105 * t193 - t111 * t182 + t117 * t183;
t54 = t104 * t182 + t110 * t193 + t116 * t183;
t53 = t103 * t182 + t109 * t193 + t115 * t183;
t52 = -t102 * t193 - t108 * t183 + t114 * t182;
t51 = -t101 * t193 - t107 * t183 + t113 * t182;
t50 = -t183 * t98 + t185 * t97;
t49 = t120 * t196 + t240 * t193 + t125;
t48 = t226 * t196 + t229 * t243;
t47 = t119 * t243 + t193 * t170 + t218;
t45 = t126 * t183 + t127 * t146 + t128 * t147;
t44 = (-t119 + t227) * t209 + t208 * t216;
t43 = t122 * t209 + t206 * t216 + t228;
t42 = (t120 * t206 + t123 * t208) * t207 + t219;
t41 = t119 * t196 + t229 * t193 + t239;
t40 = t186 * t94 + t187 * t96 + t200 * t92;
t39 = t186 * t93 + t187 * t95 + t200 * t91;
t38 = (t119 * t206 + t122 * t208) * t207 + t215;
t35 = t146 * t94 + t147 * t96 + t183 * t92;
t34 = t146 * t93 + t147 * t95 + t183 * t91;
t33 = t220 * t196 + t221 * t243;
t32 = t238 * t193 + t247 * t243 + t218;
t31 = (t227 - t247) * t209 + t208 * t214;
t30 = t206 * t214 + t246 * t209 + t228;
t29 = t209 * t88 + (t206 * t71 - t208 * t70) * t207;
t28 = t209 * t87 + (t206 * t69 - t208 * t68) * t207;
t27 = t209 * t86 + (t206 * t67 - t208 * t66) * t207;
t26 = t221 * t193 + t247 * t196 + t239;
t25 = (t247 * t206 + t246 * t208) * t207 + t215;
t24 = t70 * t193 + t71 * t196 - t88 * t243;
t23 = t68 * t193 + t69 * t196 - t87 * t243;
t22 = t66 * t193 + t67 * t196 - t86 * t243;
t21 = t209 * t82 + (t206 * t62 - t208 * t61) * t207;
t20 = t209 * t81 + (t206 * t60 - t208 * t59) * t207;
t19 = t209 * t80 + (t206 * t58 - t208 * t57) * t207;
t18 = t209 * t79 + (t206 * t56 - t208 * t55) * t207;
t17 = t209 * t78 + (t206 * t54 - t208 * t53) * t207;
t16 = t209 * t77 + (t206 * t52 - t208 * t51) * t207;
t15 = t61 * t193 + t62 * t196 - t82 * t243;
t14 = t59 * t193 + t60 * t196 - t81 * t243;
t13 = t57 * t193 + t58 * t196 - t80 * t243;
t12 = t55 * t193 + t56 * t196 - t79 * t243;
t11 = t53 * t193 + t54 * t196 - t78 * t243;
t10 = t51 * t193 + t52 * t196 - t77 * t243;
t9 = t209 * t63 + (t206 * t40 - t208 * t39) * t207;
t8 = t39 * t193 + t40 * t196 - t63 * t243;
t7 = t183 * t39 + t185 * t40 + t200 * t63;
t6 = t209 * t46 + (t206 * t37 - t208 * t36) * t207;
t5 = t209 * t45 + (t206 * t35 - t208 * t34) * t207;
t4 = t36 * t193 + t37 * t196 - t46 * t243;
t3 = t34 * t193 + t35 * t196 - t45 * t243;
t1 = t183 * t34 + t185 * t35 + t200 * t45;
t99 = [m(2) + m(3) + m(4) + t256; m(3) * t100 + m(4) * t76 + m(5) * t42 + m(6) * t38 + m(7) * t25; m(7) * (t25 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t42 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t76 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(3) * (t100 ^ 2 + t131 ^ 2 + t132 ^ 2) + (t6 + t21 + t20 + t19 + (t154 * t245 - t156 * t196 + t158 * t198) * t245) * t245 + (-t5 - t16 - t18 - t17 + (-t153 * t244 - t155 * t193 + t157 * t195) * t244 + (-t153 * t245 + t154 * t244 + t155 * t196 + t156 * t193 - t157 * t198 - t158 * t195) * t245) * t244 + ((t189 * t245 - t190 * t196 + t191 * t198) * t245 - (-t189 * t244 - t190 * t193 + t191 * t195) * t244 + t9 + t29 + t28 + t27 + ((t156 * t213 + t158 * t211) * t206 - (t155 * t213 + t157 * t211) * t208) * t257 + ((-t153 * t208 + t154 * t206 + t190 * t213 + t191 * t211) * t207 + t209 * t189) * t209) * t209; m(4) * t85 + m(5) * t49 + m(6) * t41 + m(7) * t26; (t8 / 0.2e1 + t22 / 0.2e1 + t23 / 0.2e1 + t24 / 0.2e1) * t209 + (t6 / 0.2e1 + t19 / 0.2e1 + t20 / 0.2e1 + t21 / 0.2e1) * t196 + (t5 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1 + t18 / 0.2e1) * t193 + m(7) * (t25 * t26 + t30 * t33 + t31 * t32) + m(6) * (t38 * t41 + t43 * t48 + t44 * t47) + m(5) * (t42 * t49 + t64 * t73 + t65 * t72) + m(4) * (t76 * t85 + t83 * t90 + t84 * t89) + ((-t9 / 0.2e1 - t27 / 0.2e1 - t28 / 0.2e1 - t29 / 0.2e1) * t213 + (-t3 / 0.2e1 - t10 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1) * t208 + (t4 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1) * t206) * t207; (-t22 - t23 - t24 - t8) * t243 + (t4 + t14 + t13 + t15) * t196 + (t3 + t11 + t12 + t10) * t193 + m(6) * (t41 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(5) * (t49 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(4) * (t85 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(7) * (t26 ^ 2 + t32 ^ 2 + t33 ^ 2); t199 * t256; m(7) * (t182 * t30 + t184 * t31 + t199 * t25) + m(6) * (t182 * t43 + t184 * t44 + t199 * t38) + m(5) * (t182 * t64 + t184 * t65 + t199 * t42); m(7) * (t182 * t33 + t184 * t32 + t199 * t26) + m(6) * (t182 * t48 + t184 * t47 + t199 * t41) + m(5) * (t182 * t73 + t184 * t72 + t199 * t49); (t182 ^ 2 + t184 ^ 2 + t199 ^ 2) * t256; t258 * t243; m(7) * (-t193 * t30 - t196 * t31 + t25 * t243) + m(6) * (-t193 * t43 - t196 * t44 + t38 * t243); m(7) * (-t193 * t33 - t196 * t32 + t26 * t243) + m(6) * (-t193 * t48 - t196 * t47 + t41 * t243); (-t193 * t182 - t196 * t184 + t199 * t243) * t258; (t257 * t213 ^ 2 + t193 ^ 2 + t196 ^ 2) * t258; m(7) * t50; t6 * t251 + t5 * t252 + m(7) * (t25 * t50 + t30 * t75 + t31 * t74) + t9 * t250 + t209 * t7 / 0.2e1 + (t206 * t254 - t208 * t1 / 0.2e1) * t207; m(7) * (t26 * t50 + t32 * t74 + t33 * t75) + t8 * t250 + t196 * t254 + t4 * t251 + t193 * t1 / 0.2e1 + t3 * t252 - t7 * t243 / 0.2e1; m(7) * (t182 * t75 + t184 * t74 + t199 * t50); m(7) * (-t75 * t193 - t74 * t196 + t243 * t50); t185 * t2 + t183 * t1 + t200 * t7 + m(7) * (t50 ^ 2 + t74 ^ 2 + t75 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t99(1) t99(2) t99(4) t99(7) t99(11) t99(16); t99(2) t99(3) t99(5) t99(8) t99(12) t99(17); t99(4) t99(5) t99(6) t99(9) t99(13) t99(18); t99(7) t99(8) t99(9) t99(10) t99(14) t99(19); t99(11) t99(12) t99(13) t99(14) t99(15) t99(20); t99(16) t99(17) t99(18) t99(19) t99(20) t99(21);];
Mq  = res;
