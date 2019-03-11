% Calculate joint inertia matrix for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:13
% EndTime: 2019-03-09 05:48:25
% DurationCPUTime: 5.41s
% Computational Cost: add. (30327->556), mult. (83985->773), div. (0->0), fcn. (110435->14), ass. (0->252)
t240 = sin(qJ(1));
t236 = cos(pkin(6));
t290 = cos(pkin(12));
t265 = t240 * t290;
t234 = sin(pkin(12));
t242 = cos(qJ(1));
t284 = t242 * t234;
t255 = t236 * t265 + t284;
t235 = sin(pkin(6));
t291 = cos(pkin(7));
t267 = t235 * t291;
t289 = sin(pkin(7));
t213 = t240 * t267 + t255 * t289;
t264 = t242 * t290;
t285 = t240 * t234;
t256 = -t236 * t264 + t285;
t212 = -t242 * t267 + t256 * t289;
t303 = m(7) / 0.2e1;
t304 = m(6) / 0.2e1;
t275 = t304 + t303;
t308 = 0.2e1 * t275;
t307 = m(3) / 0.2e1;
t306 = m(4) / 0.2e1;
t305 = m(5) / 0.2e1;
t221 = t236 * t284 + t265;
t239 = sin(qJ(3));
t251 = t256 * t291;
t266 = t235 * t289;
t296 = cos(qJ(3));
t204 = t221 * t296 + (-t242 * t266 - t251) * t239;
t238 = sin(qJ(4));
t295 = cos(qJ(4));
t183 = t204 * t295 + t212 * t238;
t222 = -t236 * t285 + t264;
t249 = t255 * t291;
t206 = t222 * t296 + (t240 * t266 - t249) * t239;
t185 = t206 * t295 + t213 * t238;
t260 = t291 * t290;
t211 = t236 * t289 * t239 + (t234 * t296 + t239 * t260) * t235;
t220 = t236 * t291 - t266 * t290;
t202 = t211 * t295 + t220 * t238;
t182 = t204 * t238 - t212 * t295;
t261 = t296 * t289;
t257 = t235 * t261;
t203 = t221 * t239 + t242 * t257 + t251 * t296;
t237 = sin(qJ(6));
t241 = cos(qJ(6));
t141 = t182 * t241 - t203 * t237;
t142 = t182 * t237 + t203 * t241;
t85 = Icges(7,5) * t142 + Icges(7,6) * t141 + Icges(7,3) * t183;
t87 = Icges(7,4) * t142 + Icges(7,2) * t141 + Icges(7,6) * t183;
t89 = Icges(7,1) * t142 + Icges(7,4) * t141 + Icges(7,5) * t183;
t24 = t141 * t87 + t142 * t89 + t183 * t85;
t184 = t206 * t238 - t213 * t295;
t205 = t222 * t239 - t240 * t257 + t249 * t296;
t143 = t184 * t241 - t205 * t237;
t144 = t184 * t237 + t205 * t241;
t86 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t185;
t88 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t185;
t90 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t185;
t25 = t141 * t88 + t142 * t90 + t183 * t86;
t201 = t211 * t238 - t220 * t295;
t288 = t234 * t235;
t210 = -t235 * t260 * t296 - t236 * t261 + t239 * t288;
t173 = t201 * t241 - t210 * t237;
t174 = t201 * t237 + t210 * t241;
t103 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t202;
t104 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t202;
t105 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t202;
t35 = t103 * t183 + t104 * t141 + t105 * t142;
t1 = t183 * t24 + t185 * t25 + t202 * t35;
t302 = t1 / 0.2e1;
t26 = t143 * t87 + t144 * t89 + t185 * t85;
t27 = t143 * t88 + t144 * t90 + t185 * t86;
t36 = t103 * t185 + t104 * t143 + t105 * t144;
t2 = t183 * t26 + t185 * t27 + t202 * t36;
t301 = t2 / 0.2e1;
t30 = t173 * t87 + t174 * t89 + t202 * t85;
t31 = t173 * t88 + t174 * t90 + t202 * t86;
t43 = t202 * t103 + t173 * t104 + t174 * t105;
t38 = t43 * t202;
t7 = t30 * t183 + t31 * t185 + t38;
t300 = t7 / 0.2e1;
t299 = t183 / 0.2e1;
t298 = t185 / 0.2e1;
t297 = t202 / 0.2e1;
t294 = pkin(5) * t203;
t258 = -rSges(7,1) * t142 - rSges(7,2) * t141;
t91 = rSges(7,3) * t183 - t258;
t293 = pkin(11) * t183 + t294 + t91;
t92 = t144 * rSges(7,1) + t143 * rSges(7,2) + t185 * rSges(7,3);
t292 = t205 * pkin(5) + t185 * pkin(11) + t92;
t287 = t235 * t240;
t286 = t235 * t242;
t106 = rSges(7,1) * t174 + rSges(7,2) * t173 + rSges(7,3) * t202;
t283 = pkin(5) * t210 + pkin(11) * t202 + t106;
t259 = -rSges(6,1) * t203 - rSges(6,3) * t182;
t119 = -rSges(6,2) * t183 - t259;
t175 = t182 * qJ(5);
t134 = pkin(4) * t183 + t175;
t282 = -t119 - t134;
t120 = t205 * rSges(6,1) - t185 * rSges(6,2) + t184 * rSges(6,3);
t135 = t185 * pkin(4) + t184 * qJ(5);
t281 = -t120 - t135;
t169 = pkin(3) * t204 + t203 * pkin(10);
t164 = t213 * t169;
t280 = t213 * t134 + t164;
t170 = t206 * pkin(3) + t205 * pkin(10);
t165 = t220 * t170;
t279 = t220 * t135 + t165;
t153 = rSges(6,1) * t210 - rSges(6,2) * t202 + rSges(6,3) * t201;
t168 = pkin(4) * t202 + qJ(5) * t201;
t278 = -t153 - t168;
t191 = pkin(3) * t211 + pkin(10) * t210;
t172 = t212 * t191;
t277 = t212 * t168 + t172;
t276 = t242 * pkin(1) + qJ(2) * t287;
t274 = t30 / 0.2e1 + t35 / 0.2e1;
t273 = t36 / 0.2e1 + t31 / 0.2e1;
t272 = -t134 - t293;
t271 = -t135 - t292;
t270 = -t168 - t283;
t148 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t210;
t150 = Icges(5,4) * t202 - Icges(5,2) * t201 + Icges(5,6) * t210;
t152 = Icges(5,1) * t202 - Icges(5,4) * t201 + Icges(5,5) * t210;
t78 = t210 * t148 - t201 * t150 + t202 * t152;
t147 = Icges(6,5) * t210 - Icges(6,6) * t202 + Icges(6,3) * t201;
t149 = Icges(6,4) * t210 - Icges(6,2) * t202 + Icges(6,6) * t201;
t151 = Icges(6,1) * t210 - Icges(6,4) * t202 + Icges(6,5) * t201;
t77 = t201 * t147 - t202 * t149 + t210 * t151;
t187 = Icges(4,5) * t211 - Icges(4,6) * t210 + Icges(4,3) * t220;
t188 = Icges(4,4) * t211 - Icges(4,2) * t210 + Icges(4,6) * t220;
t189 = Icges(4,1) * t211 - Icges(4,4) * t210 + Icges(4,5) * t220;
t269 = t220 * t187 - t210 * t188 + t211 * t189;
t122 = t185 * rSges(5,1) - t184 * rSges(5,2) + t205 * rSges(5,3);
t163 = t206 * rSges(4,1) - t205 * rSges(4,2) + t213 * rSges(4,3);
t268 = -t240 * pkin(1) + qJ(2) * t286;
t162 = rSges(4,1) * t204 - rSges(4,2) * t203 + rSges(4,3) * t212;
t121 = rSges(5,1) * t183 - rSges(5,2) * t182 + rSges(5,3) * t203;
t254 = -pkin(2) * t221 - pkin(9) * t212 + t268;
t107 = Icges(6,5) * t203 - Icges(6,6) * t183 + Icges(6,3) * t182;
t111 = Icges(6,4) * t203 - Icges(6,2) * t183 + Icges(6,6) * t182;
t115 = Icges(6,1) * t203 - Icges(6,4) * t183 + Icges(6,5) * t182;
t55 = t107 * t201 - t111 * t202 + t115 * t210;
t109 = Icges(5,5) * t183 - Icges(5,6) * t182 + Icges(5,3) * t203;
t113 = Icges(5,4) * t183 - Icges(5,2) * t182 + Icges(5,6) * t203;
t117 = Icges(5,1) * t183 - Icges(5,4) * t182 + Icges(5,5) * t203;
t57 = t109 * t210 - t113 * t201 + t117 * t202;
t64 = t147 * t182 - t149 * t183 + t151 * t203;
t66 = t148 * t203 - t150 * t182 + t152 * t183;
t253 = t57 / 0.2e1 + t55 / 0.2e1 + t66 / 0.2e1 + t64 / 0.2e1 + t274;
t108 = Icges(6,5) * t205 - Icges(6,6) * t185 + Icges(6,3) * t184;
t112 = Icges(6,4) * t205 - Icges(6,2) * t185 + Icges(6,6) * t184;
t116 = Icges(6,1) * t205 - Icges(6,4) * t185 + Icges(6,5) * t184;
t56 = t108 * t201 - t112 * t202 + t116 * t210;
t110 = Icges(5,5) * t185 - Icges(5,6) * t184 + Icges(5,3) * t205;
t114 = Icges(5,4) * t185 - Icges(5,2) * t184 + Icges(5,6) * t205;
t118 = Icges(5,1) * t185 - Icges(5,4) * t184 + Icges(5,5) * t205;
t58 = t110 * t210 - t114 * t201 + t118 * t202;
t65 = t147 * t184 - t149 * t185 + t151 * t205;
t67 = t148 * t205 - t150 * t184 + t152 * t185;
t252 = t67 / 0.2e1 + t65 / 0.2e1 + t58 / 0.2e1 + t56 / 0.2e1 + t273;
t247 = -t169 + t254;
t246 = -t175 + t247;
t245 = t222 * pkin(2) + pkin(9) * t213 + t276;
t244 = t170 + t245;
t243 = t135 + t244;
t228 = rSges(2,1) * t242 - t240 * rSges(2,2);
t227 = -t240 * rSges(2,1) - rSges(2,2) * t242;
t200 = t222 * rSges(3,1) - rSges(3,2) * t255 + rSges(3,3) * t287 + t276;
t199 = -t221 * rSges(3,1) + rSges(3,2) * t256 + rSges(3,3) * t286 + t268;
t190 = rSges(4,1) * t211 - rSges(4,2) * t210 + rSges(4,3) * t220;
t161 = Icges(4,1) * t206 - Icges(4,4) * t205 + Icges(4,5) * t213;
t160 = Icges(4,1) * t204 - Icges(4,4) * t203 + Icges(4,5) * t212;
t159 = Icges(4,4) * t206 - Icges(4,2) * t205 + Icges(4,6) * t213;
t158 = Icges(4,4) * t204 - Icges(4,2) * t203 + Icges(4,6) * t212;
t157 = Icges(4,5) * t206 - Icges(4,6) * t205 + Icges(4,3) * t213;
t156 = Icges(4,5) * t204 - Icges(4,6) * t203 + Icges(4,3) * t212;
t154 = rSges(5,1) * t202 - rSges(5,2) * t201 + rSges(5,3) * t210;
t138 = t203 * t168;
t132 = t245 + t163;
t131 = -t162 + t254;
t125 = t210 * t135;
t123 = t205 * t134;
t101 = t163 * t220 - t190 * t213;
t100 = -t162 * t220 + t190 * t212;
t96 = t162 * t213 - t163 * t212;
t95 = t269 * t220;
t94 = t187 * t213 - t188 * t205 + t189 * t206;
t93 = t187 * t212 - t188 * t203 + t189 * t204;
t84 = t244 + t122;
t83 = -t121 + t247;
t82 = t122 * t210 - t154 * t205;
t81 = -t121 * t210 + t154 * t203;
t80 = t157 * t220 - t159 * t210 + t161 * t211;
t79 = t156 * t220 - t158 * t210 + t160 * t211;
t76 = t121 * t205 - t122 * t203;
t75 = t243 + t120;
t74 = (rSges(6,2) - pkin(4)) * t183 + t246 + t259;
t73 = t78 * t220;
t72 = t77 * t220;
t71 = t78 * t210;
t70 = t77 * t210;
t69 = t122 * t220 + t165 + (-t154 - t191) * t213;
t68 = t154 * t212 + t172 + (-t121 - t169) * t220;
t63 = -t106 * t185 + t202 * t92;
t62 = t106 * t183 - t202 * t91;
t61 = t121 * t213 + t164 + (-t122 - t170) * t212;
t60 = t120 * t210 + t205 * t278 + t125;
t59 = t153 * t203 + t210 * t282 + t138;
t54 = t243 + t292;
t53 = -t294 + (-rSges(7,3) - pkin(4) - pkin(11)) * t183 + t246 + t258;
t52 = t110 * t205 - t114 * t184 + t118 * t185;
t51 = t109 * t205 - t113 * t184 + t117 * t185;
t50 = t110 * t203 - t114 * t182 + t118 * t183;
t49 = t109 * t203 - t113 * t182 + t117 * t183;
t48 = t108 * t184 - t112 * t185 + t116 * t205;
t47 = t107 * t184 - t111 * t185 + t115 * t205;
t46 = t108 * t182 - t112 * t183 + t116 * t203;
t45 = t107 * t182 - t111 * t183 + t115 * t203;
t44 = -t183 * t92 + t185 * t91;
t42 = t43 * t220;
t41 = t120 * t220 + (-t191 + t278) * t213 + t279;
t40 = t153 * t212 + (-t169 + t282) * t220 + t277;
t39 = t43 * t210;
t37 = t119 * t205 + t203 * t281 + t123;
t34 = t119 * t213 + (-t170 + t281) * t212 + t280;
t33 = t205 * t270 + t210 * t292 + t125;
t32 = t203 * t283 + t210 * t272 + t138;
t29 = t292 * t220 + (-t191 + t270) * t213 + t279;
t28 = t283 * t212 + (-t169 + t272) * t220 + t277;
t23 = t203 * t271 + t205 * t293 + t123;
t22 = t293 * t213 + (-t170 + t271) * t212 + t280;
t21 = t57 * t212 + t58 * t213 + t73;
t20 = t55 * t212 + t56 * t213 + t72;
t19 = t57 * t203 + t58 * t205 + t71;
t18 = t55 * t203 + t56 * t205 + t70;
t17 = t212 * t51 + t213 * t52 + t220 * t67;
t16 = t212 * t49 + t213 * t50 + t220 * t66;
t15 = t212 * t47 + t213 * t48 + t220 * t65;
t14 = t212 * t45 + t213 * t46 + t220 * t64;
t13 = t203 * t51 + t205 * t52 + t210 * t67;
t12 = t203 * t49 + t205 * t50 + t210 * t66;
t11 = t203 * t47 + t205 * t48 + t210 * t65;
t10 = t203 * t45 + t205 * t46 + t210 * t64;
t9 = t30 * t212 + t31 * t213 + t42;
t8 = t30 * t203 + t31 * t205 + t39;
t6 = t212 * t26 + t213 * t27 + t220 * t36;
t5 = t212 * t24 + t213 * t25 + t220 * t35;
t4 = t203 * t26 + t205 * t27 + t210 * t36;
t3 = t203 * t24 + t205 * t25 + t210 * t35;
t97 = [t269 + t236 * (Icges(3,3) * t236 + (Icges(3,5) * t234 + Icges(3,6) * t290) * t235) + t235 * t290 * (Icges(3,6) * t236 + (Icges(3,4) * t234 + Icges(3,2) * t290) * t235) + m(7) * (t53 ^ 2 + t54 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t83 ^ 2 + t84 ^ 2) + m(4) * (t131 ^ 2 + t132 ^ 2) + m(3) * (t199 ^ 2 + t200 ^ 2) + m(2) * (t227 ^ 2 + t228 ^ 2) + Icges(2,3) + (Icges(3,5) * t236 + (Icges(3,1) * t234 + Icges(3,4) * t290) * t235) * t288 + t77 + t78 + t43; 0.2e1 * ((t240 * t53 - t242 * t54) * t303 + (t240 * t74 - t242 * t75) * t304 + (t240 * t83 - t242 * t84) * t305 + (t131 * t240 - t132 * t242) * t306 + (t199 * t240 - t200 * t242) * t307) * t235; 0.2e1 * (t306 + t305 + t307 + t275) * (t236 ^ 2 + (t240 ^ 2 + t242 ^ 2) * t235 ^ 2); t42 + t72 + t73 + t95 + m(7) * (t28 * t53 + t29 * t54) + m(6) * (t40 * t74 + t41 * t75) + m(5) * (t68 * t83 + t69 * t84) + m(4) * (t100 * t131 + t101 * t132) + (t80 / 0.2e1 + t94 / 0.2e1 + t252) * t213 + (t93 / 0.2e1 + t79 / 0.2e1 + t253) * t212; m(4) * (t96 * t236 + (t100 * t240 - t101 * t242) * t235) + m(5) * (t61 * t236 + (t240 * t68 - t242 * t69) * t235) + m(6) * (t34 * t236 + (t240 * t40 - t242 * t41) * t235) + m(7) * (t22 * t236 + (t240 * t28 - t242 * t29) * t235); (t20 + t21 + t9 + t95) * t220 + (t15 + t17 + t6 + (t213 * t157 - t205 * t159 + t206 * t161) * t213 + (t80 + t94) * t220) * t213 + (t5 + t16 + t14 + (t212 * t156 - t203 * t158 + t204 * t160) * t212 + (t79 + t93) * t220 + (t213 * t156 + t212 * t157 - t205 * t158 - t203 * t159 + t206 * t160 + t204 * t161) * t213) * t212 + m(7) * (t22 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t34 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t61 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t100 ^ 2 + t101 ^ 2 + t96 ^ 2); t39 + t70 + t71 + m(7) * (t32 * t53 + t33 * t54) + m(6) * (t59 * t74 + t60 * t75) + m(5) * (t81 * t83 + t82 * t84) + t252 * t205 + t253 * t203; m(5) * (t76 * t236 + (t240 * t81 - t242 * t82) * t235) + m(6) * (t37 * t236 + (t240 * t59 - t242 * t60) * t235) + m(7) * (t23 * t236 + (t240 * t32 - t242 * t33) * t235); (t8 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t220 + (t4 / 0.2e1 + t13 / 0.2e1 + t11 / 0.2e1) * t213 + (t3 / 0.2e1 + t12 / 0.2e1 + t10 / 0.2e1) * t212 + (t9 / 0.2e1 + t21 / 0.2e1 + t20 / 0.2e1) * t210 + (t6 / 0.2e1 + t17 / 0.2e1 + t15 / 0.2e1) * t205 + (t5 / 0.2e1 + t16 / 0.2e1 + t14 / 0.2e1) * t203 + m(7) * (t22 * t23 + t28 * t32 + t29 * t33) + m(6) * (t34 * t37 + t40 * t59 + t41 * t60) + m(5) * (t61 * t76 + t68 * t81 + t69 * t82); (t8 + t18 + t19) * t210 + (t4 + t13 + t11) * t205 + (t3 + t12 + t10) * t203 + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t37 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t76 ^ 2 + t81 ^ 2 + t82 ^ 2); m(7) * (t182 * t54 + t184 * t53) + m(6) * (t182 * t75 + t184 * t74); (t201 * t236 + (-t182 * t242 + t184 * t240) * t235) * t308; m(7) * (t182 * t29 + t184 * t28 + t201 * t22) + m(6) * (t182 * t41 + t184 * t40 + t201 * t34); m(7) * (t182 * t33 + t184 * t32 + t201 * t23) + m(6) * (t182 * t60 + t184 * t59 + t201 * t37); (t182 ^ 2 + t184 ^ 2 + t201 ^ 2) * t308; m(7) * (t53 * t62 + t54 * t63) + t38 + t273 * t185 + t274 * t183; m(7) * (t44 * t236 + (t240 * t62 - t242 * t63) * t235); m(7) * (t22 * t44 + t28 * t62 + t29 * t63) + t6 * t298 + t212 * t302 + t9 * t297 + t220 * t300 + t5 * t299 + t213 * t301; m(7) * (t23 * t44 + t32 * t62 + t33 * t63) + t8 * t297 + t4 * t298 + t205 * t301 + t203 * t302 + t3 * t299 + t210 * t300; m(7) * (t182 * t63 + t184 * t62 + t201 * t44); t202 * t7 + t185 * t2 + t183 * t1 + m(7) * (t44 ^ 2 + t62 ^ 2 + t63 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t97(1) t97(2) t97(4) t97(7) t97(11) t97(16); t97(2) t97(3) t97(5) t97(8) t97(12) t97(17); t97(4) t97(5) t97(6) t97(9) t97(13) t97(18); t97(7) t97(8) t97(9) t97(10) t97(14) t97(19); t97(11) t97(12) t97(13) t97(14) t97(15) t97(20); t97(16) t97(17) t97(18) t97(19) t97(20) t97(21);];
Mq  = res;
