% Calculate joint inertia matrix for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:16
% EndTime: 2019-03-08 20:18:26
% DurationCPUTime: 4.95s
% Computational Cost: add. (12895->484), mult. (33352->697), div. (0->0), fcn. (42504->10), ass. (0->211)
t252 = Icges(3,1) + Icges(4,2);
t250 = Icges(3,4) + Icges(4,6);
t249 = Icges(3,5) - Icges(4,4);
t251 = Icges(3,2) + Icges(4,3);
t248 = Icges(3,6) - Icges(4,5);
t247 = Icges(3,3) + Icges(4,1);
t246 = rSges(7,1) + pkin(5);
t245 = rSges(7,3) + qJ(6);
t192 = cos(pkin(6));
t189 = sin(pkin(10));
t195 = sin(qJ(2));
t223 = t195 * t189;
t191 = cos(pkin(10));
t196 = cos(qJ(2));
t224 = t191 * t196;
t178 = -t192 * t224 + t223;
t222 = t196 * t189;
t225 = t191 * t195;
t179 = t192 * t225 + t222;
t190 = sin(pkin(6));
t229 = t190 * t191;
t244 = t251 * t178 - t250 * t179 + t248 * t229;
t243 = t250 * t178 - t252 * t179 + t249 * t229;
t180 = t192 * t222 + t225;
t181 = -t192 * t223 + t224;
t230 = t189 * t190;
t242 = t251 * t180 - t250 * t181 - t248 * t230;
t241 = -t250 * t180 + t252 * t181 + t249 * t230;
t240 = t248 * t178 - t249 * t179 + t247 * t229;
t239 = -t248 * t180 + t249 * t181 + t247 * t230;
t238 = t247 * t192 + (t249 * t195 + t248 * t196) * t190;
t237 = t248 * t192 + (t250 * t195 + t251 * t196) * t190;
t236 = t249 * t192 + (t252 * t195 + t250 * t196) * t190;
t234 = cos(qJ(4));
t233 = cos(qJ(5));
t194 = sin(qJ(4));
t203 = t190 * t234;
t160 = t180 * t194 + t189 * t203;
t193 = sin(qJ(5));
t124 = t160 * t193 - t181 * t233;
t125 = t160 * t233 + t181 * t193;
t228 = t190 * t194;
t159 = -t180 * t234 + t189 * t228;
t232 = rSges(7,2) * t159 + t245 * t124 + t246 * t125;
t162 = t178 * t194 - t191 * t203;
t126 = t162 * t193 - t179 * t233;
t127 = t162 * t233 + t179 * t193;
t161 = t178 * t234 + t191 * t228;
t231 = -rSges(7,2) * t161 + t245 * t126 + t246 * t127;
t227 = t190 * t195;
t226 = t190 * t196;
t183 = t192 * t234 - t194 * t226;
t163 = t183 * t193 - t227 * t233;
t164 = t183 * t233 + t193 * t227;
t182 = t192 * t194 + t196 * t203;
t221 = rSges(7,2) * t182 + t245 * t163 + t246 * t164;
t153 = pkin(2) * t179 + qJ(3) * t178;
t154 = pkin(2) * t181 + qJ(3) * t180;
t220 = t153 * t230 + t154 * t229;
t152 = t192 * t154;
t165 = pkin(3) * t230 + pkin(8) * t181;
t219 = t192 * t165 + t152;
t166 = -pkin(3) * t229 + pkin(8) * t179;
t218 = -t153 - t166;
t184 = (pkin(2) * t195 - qJ(3) * t196) * t190;
t217 = -pkin(3) * t192 - pkin(8) * t227 - t184;
t77 = Icges(7,5) * t125 + Icges(7,6) * t159 + Icges(7,3) * t124;
t81 = Icges(7,4) * t125 + Icges(7,2) * t159 + Icges(7,6) * t124;
t85 = Icges(7,1) * t125 + Icges(7,4) * t159 + Icges(7,5) * t124;
t31 = t124 * t77 + t125 * t85 + t159 * t81;
t78 = Icges(7,5) * t127 - Icges(7,6) * t161 + Icges(7,3) * t126;
t82 = Icges(7,4) * t127 - Icges(7,2) * t161 + Icges(7,6) * t126;
t86 = Icges(7,1) * t127 - Icges(7,4) * t161 + Icges(7,5) * t126;
t32 = t124 * t78 + t125 * t86 + t159 * t82;
t106 = Icges(7,5) * t164 + Icges(7,6) * t182 + Icges(7,3) * t163;
t108 = Icges(7,4) * t164 + Icges(7,2) * t182 + Icges(7,6) * t163;
t110 = Icges(7,1) * t164 + Icges(7,4) * t182 + Icges(7,5) * t163;
t50 = t106 * t124 + t108 * t159 + t110 * t125;
t1 = t159 * t31 - t161 * t32 + t182 * t50;
t79 = Icges(6,5) * t125 - Icges(6,6) * t124 + Icges(6,3) * t159;
t83 = Icges(6,4) * t125 - Icges(6,2) * t124 + Icges(6,6) * t159;
t87 = Icges(6,1) * t125 - Icges(6,4) * t124 + Icges(6,5) * t159;
t33 = -t124 * t83 + t125 * t87 + t159 * t79;
t80 = Icges(6,5) * t127 - Icges(6,6) * t126 - Icges(6,3) * t161;
t84 = Icges(6,4) * t127 - Icges(6,2) * t126 - Icges(6,6) * t161;
t88 = Icges(6,1) * t127 - Icges(6,4) * t126 - Icges(6,5) * t161;
t34 = -t124 * t84 + t125 * t88 + t159 * t80;
t107 = Icges(6,5) * t164 - Icges(6,6) * t163 + Icges(6,3) * t182;
t109 = Icges(6,4) * t164 - Icges(6,2) * t163 + Icges(6,6) * t182;
t111 = Icges(6,1) * t164 - Icges(6,4) * t163 + Icges(6,5) * t182;
t51 = t107 * t159 - t109 * t124 + t111 * t125;
t2 = t159 * t33 - t161 * t34 + t182 * t51;
t216 = t2 / 0.2e1 + t1 / 0.2e1;
t35 = t126 * t77 + t127 * t85 - t161 * t81;
t36 = t126 * t78 + t127 * t86 - t161 * t82;
t52 = t106 * t126 - t108 * t161 + t110 * t127;
t3 = t159 * t35 - t161 * t36 + t182 * t52;
t37 = -t126 * t83 + t127 * t87 - t161 * t79;
t38 = -t126 * t84 + t127 * t88 - t161 * t80;
t53 = -t107 * t161 - t109 * t126 + t111 * t127;
t4 = t159 * t37 - t161 * t38 + t182 * t53;
t215 = -t4 / 0.2e1 - t3 / 0.2e1;
t5 = t179 * t32 + t181 * t31 + t227 * t50;
t6 = t179 * t34 + t181 * t33 + t227 * t51;
t214 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t179 * t36 + t181 * t35 + t227 * t52;
t8 = t179 * t38 + t181 * t37 + t227 * t53;
t213 = -t8 / 0.2e1 - t7 / 0.2e1;
t10 = t192 * t51 + (t189 * t33 - t191 * t34) * t190;
t9 = t192 * t50 + (t189 * t31 - t191 * t32) * t190;
t212 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t192 * t52 + (t189 * t35 - t191 * t36) * t190;
t12 = t192 * t53 + (t189 * t37 - t191 * t38) * t190;
t211 = -t12 / 0.2e1 - t11 / 0.2e1;
t41 = t163 * t77 + t164 * t85 + t182 * t81;
t42 = t163 * t78 + t164 * t86 + t182 * t82;
t62 = t106 * t163 + t108 * t182 + t110 * t164;
t13 = t159 * t41 - t161 * t42 + t182 * t62;
t43 = -t163 * t83 + t164 * t87 + t182 * t79;
t44 = -t163 * t84 + t164 * t88 + t182 * t80;
t63 = t107 * t182 - t109 * t163 + t111 * t164;
t14 = t159 * t43 - t161 * t44 + t182 * t63;
t210 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t179 * t42 + t181 * t41 + t227 * t62;
t16 = t179 * t44 + t181 * t43 + t227 * t63;
t209 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t192 * t62 + (t189 * t41 - t191 * t42) * t190;
t18 = t192 * t63 + (t189 * t43 - t191 * t44) * t190;
t208 = t18 / 0.2e1 + t17 / 0.2e1;
t207 = -m(4) - m(5) - m(6) - m(7);
t121 = pkin(4) * t160 + pkin(9) * t159;
t206 = t192 * t121 + t219;
t122 = pkin(4) * t162 - pkin(9) * t161;
t205 = -t122 + t218;
t155 = pkin(4) * t183 + pkin(9) * t182;
t204 = -t155 + t217;
t202 = (-t192 * rSges(4,1) - (-rSges(4,2) * t195 - rSges(4,3) * t196) * t190 - t184) * t190;
t201 = t165 * t229 + t166 * t230 + t220;
t148 = rSges(5,1) * t183 - rSges(5,2) * t182 + rSges(5,3) * t227;
t200 = (-t148 + t217) * t190;
t114 = rSges(6,1) * t164 - rSges(6,2) * t163 + rSges(6,3) * t182;
t199 = (-t114 + t204) * t190;
t198 = t121 * t229 + t122 * t230 + t201;
t197 = (t204 - t221) * t190;
t173 = t192 * rSges(3,3) + (rSges(3,1) * t195 + rSges(3,2) * t196) * t190;
t147 = Icges(5,1) * t183 - Icges(5,4) * t182 + Icges(5,5) * t227;
t146 = Icges(5,4) * t183 - Icges(5,2) * t182 + Icges(5,6) * t227;
t145 = Icges(5,5) * t183 - Icges(5,6) * t182 + Icges(5,3) * t227;
t144 = rSges(3,1) * t181 - rSges(3,2) * t180 + rSges(3,3) * t230;
t143 = rSges(3,1) * t179 - rSges(3,2) * t178 - rSges(3,3) * t229;
t142 = -rSges(4,1) * t229 - rSges(4,2) * t179 + rSges(4,3) * t178;
t141 = rSges(4,1) * t230 - rSges(4,2) * t181 + rSges(4,3) * t180;
t128 = t179 * t155;
t119 = t121 * t227;
t116 = -t143 * t192 - t173 * t229;
t115 = t144 * t192 - t173 * t230;
t112 = t181 * t122;
t105 = rSges(5,1) * t162 + rSges(5,2) * t161 + rSges(5,3) * t179;
t104 = rSges(5,1) * t160 - rSges(5,2) * t159 + rSges(5,3) * t181;
t103 = Icges(5,1) * t162 + Icges(5,4) * t161 + Icges(5,5) * t179;
t102 = Icges(5,1) * t160 - Icges(5,4) * t159 + Icges(5,5) * t181;
t101 = Icges(5,4) * t162 + Icges(5,2) * t161 + Icges(5,6) * t179;
t100 = Icges(5,4) * t160 - Icges(5,2) * t159 + Icges(5,6) * t181;
t99 = Icges(5,5) * t162 + Icges(5,6) * t161 + Icges(5,3) * t179;
t98 = Icges(5,5) * t160 - Icges(5,6) * t159 + Icges(5,3) * t181;
t97 = (t143 * t189 + t144 * t191) * t190;
t94 = (-t142 - t153) * t192 + t191 * t202;
t93 = t141 * t192 + t189 * t202 + t152;
t92 = rSges(6,1) * t127 - rSges(6,2) * t126 - rSges(6,3) * t161;
t90 = rSges(6,1) * t125 - rSges(6,2) * t124 + rSges(6,3) * t159;
t76 = t104 * t227 - t148 * t181;
t75 = -t105 * t227 + t148 * t179;
t74 = t145 * t227 - t146 * t182 + t147 * t183;
t73 = (t141 * t191 + t142 * t189) * t190 + t220;
t72 = -t104 * t179 + t105 * t181;
t71 = t145 * t179 + t146 * t161 + t147 * t162;
t70 = t145 * t181 - t146 * t159 + t147 * t160;
t69 = (-t105 + t218) * t192 + t191 * t200;
t68 = t104 * t192 + t189 * t200 + t219;
t67 = -t114 * t161 - t182 * t92;
t66 = -t114 * t159 + t182 * t90;
t65 = -t101 * t182 + t103 * t183 + t227 * t99;
t64 = -t100 * t182 + t102 * t183 + t227 * t98;
t61 = t101 * t161 + t103 * t162 + t179 * t99;
t60 = t100 * t161 + t102 * t162 + t179 * t98;
t59 = -t101 * t159 + t103 * t160 + t181 * t99;
t58 = -t100 * t159 + t102 * t160 + t181 * t98;
t57 = (t104 * t191 + t105 * t189) * t190 + t201;
t56 = t159 * t92 + t161 * t90;
t55 = t90 * t227 + t119 + (-t114 - t155) * t181;
t54 = t114 * t179 + t128 + (-t122 - t92) * t227;
t49 = (-t92 + t205) * t192 + t191 * t199;
t48 = t189 * t199 + t192 * t90 + t206;
t47 = t181 * t92 + t112 + (-t121 - t90) * t179;
t46 = -t221 * t161 - t182 * t231;
t45 = -t221 * t159 + t182 * t232;
t40 = t119 + t232 * t227 + (-t155 - t221) * t181;
t39 = t128 + t221 * t179 + (-t122 - t231) * t227;
t30 = (t189 * t92 + t191 * t90) * t190 + t198;
t29 = t159 * t231 + t161 * t232;
t28 = (t205 - t231) * t192 + t191 * t197;
t27 = t189 * t197 + t192 * t232 + t206;
t26 = t112 + t231 * t181 + (-t121 - t232) * t179;
t25 = t192 * t74 + (t189 * t64 - t191 * t65) * t190;
t24 = t179 * t65 + t181 * t64 + t227 * t74;
t23 = (t189 * t231 + t191 * t232) * t190 + t198;
t22 = t192 * t71 + (t189 * t60 - t191 * t61) * t190;
t21 = t192 * t70 + (t189 * t58 - t191 * t59) * t190;
t20 = t179 * t61 + t181 * t60 + t227 * t71;
t19 = t179 * t59 + t181 * t58 + t227 * t70;
t89 = [m(2) + m(3) - t207; m(3) * t97 + m(4) * t73 + m(5) * t57 + m(6) * t30 + m(7) * t23; m(7) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t30 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t57 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t73 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(3) * (t115 ^ 2 + t116 ^ 2 + t97 ^ 2) + (t17 + t18 + t25 + t238 * t192 ^ 2 + ((t195 * t236 + t196 * t237) * t192 + (t240 * t192 + (t243 * t195 + t244 * t196) * t190) * t191 + (t239 * t192 + (t241 * t195 - t242 * t196) * t190) * t189) * t190) * t192 + (t9 + t10 + t21 + (t242 * t180 + t241 * t181 + t239 * t230) * t230 + (-t180 * t237 + t181 * t236 + t230 * t238) * t192) * t230 + (-t12 - t11 - t22 + (t244 * t178 - t243 * t179 + t240 * t229) * t229 + (t178 * t237 - t179 * t236 + t229 * t238) * t192 + (-t242 * t178 - t241 * t179 - t244 * t180 + t243 * t181 + t239 * t229 + t240 * t230) * t230) * t229; t207 * t226; m(7) * (t178 * t27 + t180 * t28 - t226 * t23) + m(6) * (t178 * t48 + t180 * t49 - t226 * t30) + m(5) * (t178 * t68 + t180 * t69 - t226 * t57) + m(4) * (t178 * t93 + t180 * t94 - t226 * t73); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t190 ^ 2 * t196 ^ 2 + t178 ^ 2 + t180 ^ 2); m(5) * t72 + m(6) * t47 + m(7) * t26; (t24 / 0.2e1 + t209) * t192 + (t21 / 0.2e1 + t212) * t181 + (t22 / 0.2e1 - t211) * t179 + m(7) * (t23 * t26 + t27 * t40 + t28 * t39) + m(6) * (t30 * t47 + t48 * t55 + t49 * t54) + m(5) * (t57 * t72 + t68 * t76 + t69 * t75) + ((t25 / 0.2e1 + t208) * t195 + (-t20 / 0.2e1 + t213) * t191 + (t19 / 0.2e1 + t214) * t189) * t190; m(5) * (t76 * t178 + t75 * t180 - t226 * t72) + m(6) * (t55 * t178 + t54 * t180 - t226 * t47) + m(7) * (t40 * t178 + t39 * t180 - t226 * t26); (t15 + t16 + t24) * t227 + (t6 + t5 + t19) * t181 + (t7 + t8 + t20) * t179 + m(7) * (t26 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2); m(6) * t56 + m(7) * t29; t210 * t192 + t208 * t182 + t211 * t161 + t212 * t159 + m(7) * (t23 * t29 + t27 * t45 + t28 * t46) + m(6) * (t30 * t56 + t48 * t66 + t49 * t67) + (t189 * t216 + t191 * t215) * t190; m(6) * (t66 * t178 + t67 * t180 - t226 * t56) + m(7) * (t45 * t178 + t46 * t180 - t226 * t29); t210 * t227 + t209 * t182 + t216 * t181 - t215 * t179 + t213 * t161 + t214 * t159 + m(7) * (t26 * t29 + t39 * t46 + t40 * t45) + m(6) * (t47 * t56 + t54 * t67 + t55 * t66); (t14 + t13) * t182 + (-t4 - t3) * t161 + (t2 + t1) * t159 + m(7) * (t29 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * t163; m(7) * (t124 * t28 + t126 * t27 + t163 * t23); m(7) * (t124 * t180 + t126 * t178 - t163 * t226); m(7) * (t124 * t39 + t126 * t40 + t163 * t26); m(7) * (t124 * t46 + t126 * t45 + t163 * t29); m(7) * (t124 ^ 2 + t126 ^ 2 + t163 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t89(1) t89(2) t89(4) t89(7) t89(11) t89(16); t89(2) t89(3) t89(5) t89(8) t89(12) t89(17); t89(4) t89(5) t89(6) t89(9) t89(13) t89(18); t89(7) t89(8) t89(9) t89(10) t89(14) t89(19); t89(11) t89(12) t89(13) t89(14) t89(15) t89(20); t89(16) t89(17) t89(18) t89(19) t89(20) t89(21);];
Mq  = res;
