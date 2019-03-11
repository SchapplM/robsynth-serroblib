% Calculate joint inertia matrix for
% S6PRPRRP5
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:05
% EndTime: 2019-03-08 20:14:14
% DurationCPUTime: 4.87s
% Computational Cost: add. (13415->489), mult. (34124->704), div. (0->0), fcn. (43316->10), ass. (0->214)
t256 = Icges(3,1) + Icges(4,2);
t253 = Icges(3,4) + Icges(4,6);
t252 = Icges(4,4) - Icges(3,5);
t251 = Icges(4,5) - Icges(3,6);
t255 = Icges(3,2) + Icges(4,3);
t254 = Icges(4,1) + Icges(3,3);
t250 = rSges(7,3) + qJ(6);
t191 = cos(pkin(6));
t190 = cos(pkin(10));
t197 = cos(qJ(2));
t223 = t197 * t190;
t188 = sin(pkin(10));
t195 = sin(qJ(2));
t225 = t195 * t188;
t177 = -t191 * t223 + t225;
t224 = t197 * t188;
t226 = t190 * t195;
t178 = t191 * t226 + t224;
t189 = sin(pkin(6));
t230 = t189 * t190;
t249 = t255 * t177 - t253 * t178 - t251 * t230;
t248 = t253 * t177 - t256 * t178 - t252 * t230;
t179 = t191 * t224 + t226;
t180 = -t191 * t225 + t223;
t231 = t188 * t189;
t247 = t255 * t179 - t253 * t180 + t251 * t231;
t246 = -t253 * t179 + t256 * t180 - t252 * t231;
t245 = -t251 * t177 + t252 * t178 + t254 * t230;
t244 = t251 * t179 - t252 * t180 + t254 * t231;
t243 = -t251 * t191 + (t253 * t195 + t255 * t197) * t189;
t242 = -t252 * t191 + (t256 * t195 + t253 * t197) * t189;
t241 = t254 * t191 + (-t252 * t195 - t251 * t197) * t189;
t239 = cos(qJ(4));
t196 = cos(qJ(5));
t238 = pkin(5) * t196;
t194 = sin(qJ(4));
t204 = t189 * t239;
t160 = t179 * t194 + t188 * t204;
t193 = sin(qJ(5));
t124 = -t160 * t193 + t180 * t196;
t232 = t180 * t193;
t125 = t160 * t196 + t232;
t229 = t189 * t194;
t159 = -t179 * t239 + t188 * t229;
t236 = rSges(7,1) * t125 + rSges(7,2) * t124 + pkin(5) * t232 + t250 * t159 + t160 * t238;
t163 = t177 * t194 - t190 * t204;
t126 = -t163 * t193 + t178 * t196;
t233 = t178 * t193;
t127 = t163 * t196 + t233;
t161 = t177 * t239 + t190 * t229;
t235 = rSges(7,1) * t127 + rSges(7,2) * t126 + pkin(5) * t233 - t250 * t161 + t163 * t238;
t227 = t189 * t197;
t182 = t191 * t239 - t194 * t227;
t228 = t189 * t195;
t164 = -t182 * t193 + t196 * t228;
t208 = t193 * t228;
t165 = t182 * t196 + t208;
t181 = t191 * t194 + t197 * t204;
t234 = rSges(7,1) * t165 + rSges(7,2) * t164 + pkin(5) * t208 + t250 * t181 + t182 * t238;
t153 = pkin(2) * t178 + qJ(3) * t177;
t154 = pkin(2) * t180 + qJ(3) * t179;
t222 = t153 * t231 + t154 * t230;
t152 = t191 * t154;
t166 = pkin(3) * t231 + pkin(8) * t180;
t221 = t191 * t166 + t152;
t167 = -pkin(3) * t230 + pkin(8) * t178;
t220 = -t153 - t167;
t183 = (pkin(2) * t195 - qJ(3) * t197) * t189;
t219 = -pkin(3) * t191 - pkin(8) * t228 - t183;
t79 = Icges(7,5) * t125 + Icges(7,6) * t124 + Icges(7,3) * t159;
t83 = Icges(7,4) * t125 + Icges(7,2) * t124 + Icges(7,6) * t159;
t87 = Icges(7,1) * t125 + Icges(7,4) * t124 + Icges(7,5) * t159;
t33 = t124 * t83 + t125 * t87 + t159 * t79;
t80 = Icges(7,5) * t127 + Icges(7,6) * t126 - Icges(7,3) * t161;
t84 = Icges(7,4) * t127 + Icges(7,2) * t126 - Icges(7,6) * t161;
t88 = Icges(7,1) * t127 + Icges(7,4) * t126 - Icges(7,5) * t161;
t34 = t124 * t84 + t125 * t88 + t159 * t80;
t107 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t181;
t109 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t181;
t111 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t181;
t50 = t107 * t159 + t109 * t124 + t111 * t125;
t1 = t159 * t33 - t161 * t34 + t181 * t50;
t81 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t159;
t85 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t159;
t89 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t159;
t35 = t124 * t85 + t125 * t89 + t159 * t81;
t82 = Icges(6,5) * t127 + Icges(6,6) * t126 - Icges(6,3) * t161;
t86 = Icges(6,4) * t127 + Icges(6,2) * t126 - Icges(6,6) * t161;
t90 = Icges(6,1) * t127 + Icges(6,4) * t126 - Icges(6,5) * t161;
t36 = t124 * t86 + t125 * t90 + t159 * t82;
t108 = Icges(6,5) * t165 + Icges(6,6) * t164 + Icges(6,3) * t181;
t110 = Icges(6,4) * t165 + Icges(6,2) * t164 + Icges(6,6) * t181;
t112 = Icges(6,1) * t165 + Icges(6,4) * t164 + Icges(6,5) * t181;
t51 = t108 * t159 + t110 * t124 + t112 * t125;
t2 = t159 * t35 - t161 * t36 + t181 * t51;
t218 = t1 / 0.2e1 + t2 / 0.2e1;
t37 = t126 * t83 + t127 * t87 - t161 * t79;
t38 = t126 * t84 + t127 * t88 - t161 * t80;
t52 = -t107 * t161 + t109 * t126 + t111 * t127;
t3 = t159 * t37 - t161 * t38 + t181 * t52;
t39 = t126 * t85 + t127 * t89 - t161 * t81;
t40 = t126 * t86 + t127 * t90 - t161 * t82;
t53 = -t108 * t161 + t110 * t126 + t112 * t127;
t4 = t159 * t39 - t161 * t40 + t181 * t53;
t217 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t178 * t34 + t180 * t33 + t228 * t50;
t6 = t178 * t36 + t180 * t35 + t228 * t51;
t216 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t178 * t38 + t180 * t37 + t228 * t52;
t8 = t178 * t40 + t180 * t39 + t228 * t53;
t215 = -t8 / 0.2e1 - t7 / 0.2e1;
t10 = t191 * t51 + (t188 * t35 - t190 * t36) * t189;
t9 = t191 * t50 + (t188 * t33 - t190 * t34) * t189;
t214 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t191 * t52 + (t188 * t37 - t190 * t38) * t189;
t12 = t191 * t53 + (t188 * t39 - t190 * t40) * t189;
t213 = -t11 / 0.2e1 - t12 / 0.2e1;
t43 = t164 * t83 + t165 * t87 + t181 * t79;
t44 = t164 * t84 + t165 * t88 + t181 * t80;
t62 = t107 * t181 + t109 * t164 + t111 * t165;
t13 = t159 * t43 - t161 * t44 + t181 * t62;
t45 = t164 * t85 + t165 * t89 + t181 * t81;
t46 = t164 * t86 + t165 * t90 + t181 * t82;
t63 = t108 * t181 + t110 * t164 + t112 * t165;
t14 = t159 * t45 - t161 * t46 + t181 * t63;
t212 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t178 * t44 + t180 * t43 + t228 * t62;
t16 = t178 * t46 + t180 * t45 + t228 * t63;
t211 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t191 * t62 + (t188 * t43 - t190 * t44) * t189;
t18 = t191 * t63 + (t188 * t45 - t190 * t46) * t189;
t210 = t17 / 0.2e1 + t18 / 0.2e1;
t209 = -m(4) - m(5) - m(6) - m(7);
t122 = pkin(4) * t160 + pkin(9) * t159;
t207 = t191 * t122 + t221;
t123 = pkin(4) * t163 - pkin(9) * t161;
t206 = -t123 + t220;
t155 = pkin(4) * t182 + pkin(9) * t181;
t205 = -t155 + t219;
t203 = (-t191 * rSges(4,1) - (-rSges(4,2) * t195 - rSges(4,3) * t197) * t189 - t183) * t189;
t202 = t166 * t230 + t167 * t231 + t222;
t148 = rSges(5,1) * t182 - rSges(5,2) * t181 + rSges(5,3) * t228;
t201 = (-t148 + t219) * t189;
t115 = rSges(6,1) * t165 + rSges(6,2) * t164 + rSges(6,3) * t181;
t200 = (-t115 + t205) * t189;
t199 = t122 * t230 + t123 * t231 + t202;
t198 = (t205 - t234) * t189;
t174 = t191 * rSges(3,3) + (rSges(3,1) * t195 + rSges(3,2) * t197) * t189;
t147 = Icges(5,1) * t182 - Icges(5,4) * t181 + Icges(5,5) * t228;
t146 = Icges(5,4) * t182 - Icges(5,2) * t181 + Icges(5,6) * t228;
t145 = Icges(5,5) * t182 - Icges(5,6) * t181 + Icges(5,3) * t228;
t144 = rSges(3,1) * t180 - rSges(3,2) * t179 + rSges(3,3) * t231;
t143 = rSges(3,1) * t178 - rSges(3,2) * t177 - rSges(3,3) * t230;
t142 = -rSges(4,1) * t230 - rSges(4,2) * t178 + rSges(4,3) * t177;
t141 = rSges(4,1) * t231 - rSges(4,2) * t180 + rSges(4,3) * t179;
t128 = t178 * t155;
t120 = t122 * t228;
t117 = -t143 * t191 - t174 * t230;
t116 = t144 * t191 - t174 * t231;
t113 = t180 * t123;
t106 = rSges(5,1) * t163 + rSges(5,2) * t161 + rSges(5,3) * t178;
t105 = rSges(5,1) * t160 - rSges(5,2) * t159 + rSges(5,3) * t180;
t104 = Icges(5,1) * t163 + Icges(5,4) * t161 + Icges(5,5) * t178;
t103 = Icges(5,1) * t160 - Icges(5,4) * t159 + Icges(5,5) * t180;
t102 = Icges(5,4) * t163 + Icges(5,2) * t161 + Icges(5,6) * t178;
t101 = Icges(5,4) * t160 - Icges(5,2) * t159 + Icges(5,6) * t180;
t100 = Icges(5,5) * t163 + Icges(5,6) * t161 + Icges(5,3) * t178;
t99 = Icges(5,5) * t160 - Icges(5,6) * t159 + Icges(5,3) * t180;
t97 = (t143 * t188 + t144 * t190) * t189;
t96 = (-t142 - t153) * t191 + t190 * t203;
t95 = t141 * t191 + t188 * t203 + t152;
t94 = rSges(6,1) * t127 + rSges(6,2) * t126 - rSges(6,3) * t161;
t92 = rSges(6,1) * t125 + rSges(6,2) * t124 + rSges(6,3) * t159;
t78 = t105 * t228 - t148 * t180;
t77 = -t106 * t228 + t148 * t178;
t74 = t145 * t228 - t146 * t181 + t147 * t182;
t73 = (t141 * t190 + t142 * t188) * t189 + t222;
t72 = -t105 * t178 + t106 * t180;
t71 = t145 * t178 + t146 * t161 + t147 * t163;
t70 = t145 * t180 - t146 * t159 + t147 * t160;
t69 = (-t106 + t220) * t191 + t190 * t201;
t68 = t105 * t191 + t188 * t201 + t221;
t67 = -t115 * t161 - t181 * t94;
t66 = -t115 * t159 + t181 * t92;
t65 = t100 * t228 - t102 * t181 + t104 * t182;
t64 = -t101 * t181 + t103 * t182 + t228 * t99;
t61 = t100 * t178 + t102 * t161 + t104 * t163;
t60 = t101 * t161 + t103 * t163 + t178 * t99;
t59 = t100 * t180 - t102 * t159 + t104 * t160;
t58 = -t101 * t159 + t103 * t160 + t180 * t99;
t57 = (t105 * t190 + t106 * t188) * t189 + t202;
t56 = t159 * t94 + t161 * t92;
t55 = t92 * t228 + t120 + (-t115 - t155) * t180;
t54 = t115 * t178 + t128 + (-t123 - t94) * t228;
t49 = (-t94 + t206) * t191 + t190 * t200;
t48 = t188 * t200 + t191 * t92 + t207;
t47 = t180 * t94 + t113 + (-t122 - t92) * t178;
t42 = -t161 * t234 - t181 * t235;
t41 = -t159 * t234 + t181 * t236;
t32 = (t188 * t94 + t190 * t92) * t189 + t199;
t31 = t120 + t236 * t228 + (-t155 - t234) * t180;
t30 = t128 + t234 * t178 + (-t123 - t235) * t228;
t29 = (t206 - t235) * t191 + t190 * t198;
t28 = t188 * t198 + t191 * t236 + t207;
t27 = t159 * t235 + t161 * t236;
t26 = t191 * t74 + (t188 * t64 - t190 * t65) * t189;
t25 = t178 * t65 + t180 * t64 + t228 * t74;
t24 = t113 + t235 * t180 + (-t122 - t236) * t178;
t23 = t191 * t71 + (t188 * t60 - t190 * t61) * t189;
t22 = t191 * t70 + (t188 * t58 - t190 * t59) * t189;
t21 = t178 * t61 + t180 * t60 + t228 * t71;
t20 = t178 * t59 + t180 * t58 + t228 * t70;
t19 = (t188 * t235 + t190 * t236) * t189 + t199;
t75 = [m(2) + m(3) - t209; m(3) * t97 + m(4) * t73 + m(5) * t57 + m(6) * t32 + m(7) * t19; m(7) * (t19 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t57 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t73 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(3) * (t116 ^ 2 + t117 ^ 2 + t97 ^ 2) + (t17 + t18 + t26 + t241 * t191 ^ 2 + ((t195 * t242 + t197 * t243) * t191 + (t245 * t191 + (t248 * t195 + t249 * t197) * t189) * t190 + (t244 * t191 + (t246 * t195 - t247 * t197) * t189) * t188) * t189) * t191 + (t10 + t9 + t22 + (t247 * t179 + t246 * t180 + t244 * t231) * t231 + (-t179 * t243 + t180 * t242 + t231 * t241) * t191) * t231 + (-t11 - t12 - t23 + (t249 * t177 - t248 * t178 + t245 * t230) * t230 + (t177 * t243 - t178 * t242 + t230 * t241) * t191 + (-t247 * t177 - t246 * t178 - t249 * t179 + t248 * t180 + t244 * t230 + t245 * t231) * t231) * t230; t209 * t227; m(7) * (t177 * t28 + t179 * t29 - t19 * t227) + m(6) * (t177 * t48 + t179 * t49 - t227 * t32) + m(5) * (t177 * t68 + t179 * t69 - t227 * t57) + m(4) * (t177 * t95 + t179 * t96 - t227 * t73); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t189 ^ 2 * t197 ^ 2 + t177 ^ 2 + t179 ^ 2); m(5) * t72 + m(6) * t47 + m(7) * t24; (t25 / 0.2e1 + t211) * t191 + (t22 / 0.2e1 + t214) * t180 + (t23 / 0.2e1 - t213) * t178 + m(7) * (t19 * t24 + t28 * t31 + t29 * t30) + m(6) * (t32 * t47 + t48 * t55 + t49 * t54) + m(5) * (t57 * t72 + t68 * t78 + t69 * t77) + ((t26 / 0.2e1 + t210) * t195 + (-t21 / 0.2e1 + t215) * t190 + (t20 / 0.2e1 + t216) * t188) * t189; m(5) * (t78 * t177 + t77 * t179 - t227 * t72) + m(6) * (t55 * t177 + t54 * t179 - t227 * t47) + m(7) * (t31 * t177 + t30 * t179 - t227 * t24); (t15 + t16 + t25) * t228 + (t5 + t6 + t20) * t180 + (t7 + t8 + t21) * t178 + m(7) * (t24 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t72 ^ 2 + t77 ^ 2 + t78 ^ 2); m(6) * t56 + m(7) * t27; t212 * t191 + t210 * t181 + t213 * t161 + t214 * t159 + m(7) * (t19 * t27 + t28 * t41 + t29 * t42) + m(6) * (t32 * t56 + t48 * t66 + t49 * t67) + (t188 * t218 - t190 * t217) * t189; m(6) * (t66 * t177 + t67 * t179 - t227 * t56) + m(7) * (t41 * t177 + t42 * t179 - t227 * t27); t212 * t228 + t211 * t181 + t218 * t180 + t217 * t178 + t215 * t161 + t216 * t159 + m(7) * (t24 * t27 + t30 * t42 + t31 * t41) + m(6) * (t47 * t56 + t54 * t67 + t55 * t66); (t13 + t14) * t181 + (-t4 - t3) * t161 + (t1 + t2) * t159 + m(7) * (t27 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * t181; m(7) * (t159 * t29 - t161 * t28 + t181 * t19); m(7) * (t159 * t179 - t161 * t177 - t181 * t227); m(7) * (t159 * t30 - t161 * t31 + t181 * t24); m(7) * (t159 * t42 - t161 * t41 + t181 * t27); m(7) * (t159 ^ 2 + t161 ^ 2 + t181 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t75(1) t75(2) t75(4) t75(7) t75(11) t75(16); t75(2) t75(3) t75(5) t75(8) t75(12) t75(17); t75(4) t75(5) t75(6) t75(9) t75(13) t75(18); t75(7) t75(8) t75(9) t75(10) t75(14) t75(19); t75(11) t75(12) t75(13) t75(14) t75(15) t75(20); t75(16) t75(17) t75(18) t75(19) t75(20) t75(21);];
Mq  = res;
