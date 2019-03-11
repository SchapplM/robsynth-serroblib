% Calculate joint inertia matrix for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:01:57
% EndTime: 2019-03-09 22:02:11
% DurationCPUTime: 5.72s
% Computational Cost: add. (9898->408), mult. (8943->597), div. (0->0), fcn. (9163->10), ass. (0->204)
t316 = Icges(5,4) + Icges(6,6);
t315 = Icges(5,1) + Icges(6,2);
t314 = -Icges(5,2) - Icges(6,3);
t195 = qJ(2) + qJ(3);
t184 = qJ(4) + t195;
t180 = cos(t184);
t313 = t316 * t180;
t179 = sin(t184);
t312 = t316 * t179;
t311 = Icges(6,4) - Icges(5,5);
t310 = Icges(6,5) - Icges(5,6);
t309 = t314 * t179 + t313;
t308 = t315 * t180 - t312;
t307 = Icges(6,1) + Icges(5,3);
t198 = sin(qJ(1));
t201 = cos(qJ(1));
t306 = t198 * t309 + t201 * t310;
t305 = -t198 * t310 + t201 * t309;
t304 = t308 * t198 + t201 * t311;
t303 = -t198 * t311 + t308 * t201;
t302 = t310 * t179 - t180 * t311;
t301 = t314 * t180 - t312;
t300 = t315 * t179 + t313;
t299 = -t302 * t198 + t307 * t201;
t298 = t307 * t198 + t302 * t201;
t297 = t306 * t179 - t304 * t180;
t296 = t305 * t179 - t303 * t180;
t193 = t198 ^ 2;
t295 = t198 * pkin(7);
t265 = t180 * t201;
t199 = cos(qJ(6));
t261 = t199 * t201;
t196 = sin(qJ(6));
t264 = t196 * t198;
t139 = t179 * t261 - t264;
t262 = t198 * t199;
t263 = t196 * t201;
t140 = t179 * t263 + t262;
t76 = t140 * rSges(7,1) + t139 * rSges(7,2) + rSges(7,3) * t265;
t294 = t198 * pkin(5) + pkin(10) * t265 + t76;
t293 = -t179 * t311 - t310 * t180;
t292 = t301 * t179 + t300 * t180;
t182 = sin(t195);
t183 = cos(t195);
t241 = rSges(4,1) * t183 - rSges(4,2) * t182;
t229 = Icges(4,5) * t183 - Icges(4,6) * t182;
t125 = -Icges(4,3) * t201 + t229 * t198;
t126 = Icges(4,3) * t198 + t229 * t201;
t194 = t201 ^ 2;
t273 = Icges(4,4) * t183;
t233 = -Icges(4,2) * t182 + t273;
t128 = Icges(4,6) * t198 + t233 * t201;
t274 = Icges(4,4) * t182;
t236 = Icges(4,1) * t183 - t274;
t130 = Icges(4,5) * t198 + t236 * t201;
t220 = -t128 * t182 + t130 * t183;
t127 = -Icges(4,6) * t201 + t233 * t198;
t129 = -Icges(4,5) * t201 + t236 * t198;
t221 = t127 * t182 - t129 * t183;
t141 = t179 * t262 + t263;
t142 = t179 * t264 - t261;
t266 = t180 * t198;
t70 = Icges(7,5) * t140 + Icges(7,6) * t139 + Icges(7,3) * t265;
t72 = Icges(7,4) * t140 + Icges(7,2) * t139 + Icges(7,6) * t265;
t74 = Icges(7,1) * t140 + Icges(7,4) * t139 + Icges(7,5) * t265;
t20 = t141 * t72 + t142 * t74 + t70 * t266;
t71 = Icges(7,5) * t142 + Icges(7,6) * t141 + Icges(7,3) * t266;
t73 = Icges(7,4) * t142 + Icges(7,2) * t141 + Icges(7,6) * t266;
t75 = Icges(7,1) * t142 + Icges(7,4) * t141 + Icges(7,5) * t266;
t21 = t141 * t73 + t142 * t75 + t71 * t266;
t9 = t198 * t20 - t201 * t21;
t251 = -t9 + t299 * t194 + (t296 * t198 + (-t297 + t298) * t201) * t198;
t291 = -t194 * t125 - (t220 * t198 + (-t126 + t221) * t201) * t198 + t251;
t290 = m(6) / 0.2e1;
t289 = m(7) / 0.2e1;
t202 = -pkin(8) - pkin(7);
t288 = t198 / 0.2e1;
t287 = -t201 / 0.2e1;
t197 = sin(qJ(2));
t286 = pkin(2) * t197;
t285 = pkin(3) * t182;
t200 = cos(qJ(2));
t181 = t200 * pkin(2) + pkin(1);
t161 = pkin(3) * t183 + t181;
t155 = t201 * t161;
t175 = t201 * t181;
t284 = t201 * (t155 - t175) + (t161 - t181) * t193;
t267 = t179 * t201;
t213 = rSges(5,1) * t265 - rSges(5,2) * t267 + t198 * rSges(5,3);
t240 = rSges(5,1) * t180 - rSges(5,2) * t179;
t64 = t198 * (-t201 * rSges(5,3) + t240 * t198) + t201 * t213;
t283 = rSges(3,1) * t200;
t281 = rSges(3,2) * t197;
t279 = t201 * rSges(3,3);
t29 = t179 * t70 + (-t196 * t74 - t199 * t72) * t180;
t278 = t29 * t198;
t30 = t179 * t71 + (-t196 * t75 - t199 * t73) * t180;
t277 = t30 * t201;
t276 = Icges(3,4) * t197;
t275 = Icges(3,4) * t200;
t268 = qJ(5) * t179;
t191 = t201 * pkin(7);
t260 = t198 * (t191 + (-pkin(1) + t181) * t198) + t201 * (-t201 * pkin(1) + t175 - t295);
t214 = t198 * rSges(4,3) + t241 * t201;
t69 = t198 * (-t201 * rSges(4,3) + t241 * t198) + t201 * t214;
t257 = pkin(4) * t265 + qJ(5) * t267;
t259 = t193 * (pkin(4) * t180 + t268) + t201 * t257;
t152 = pkin(4) * t179 - qJ(5) * t180;
t258 = rSges(6,2) * t179 + rSges(6,3) * t180 - t152;
t256 = t198 * rSges(3,3) + t201 * t283;
t253 = t193 + t194;
t18 = t139 * t72 + t140 * t74 + t70 * t265;
t19 = t139 * t73 + t140 * t75 + t71 * t265;
t8 = t18 * t198 - t19 * t201;
t252 = (t8 + t298 * t193 + ((-t296 + t299) * t198 + t297 * t201) * t201) * t198;
t250 = t198 * (t193 * t126 + (t221 * t201 + (-t125 + t220) * t198) * t201) + t252;
t160 = rSges(4,1) * t182 + rSges(4,2) * t183;
t249 = -t160 - t286;
t154 = rSges(5,1) * t179 + rSges(5,2) * t180;
t248 = -t154 - t285;
t33 = t64 + t284;
t192 = -pkin(9) + t202;
t247 = -t198 * t192 + t155;
t212 = t198 * rSges(6,1) - rSges(6,2) * t265 + rSges(6,3) * t267;
t41 = t198 * (-t201 * rSges(6,1) + (-rSges(6,2) * t180 + rSges(6,3) * t179) * t198) + t201 * t212 + t259;
t92 = Icges(7,3) * t179 + (-Icges(7,5) * t196 - Icges(7,6) * t199) * t180;
t93 = Icges(7,6) * t179 + (-Icges(7,4) * t196 - Icges(7,2) * t199) * t180;
t94 = Icges(7,5) * t179 + (-Icges(7,1) * t196 - Icges(7,4) * t199) * t180;
t34 = t139 * t93 + t140 * t94 + t92 * t265;
t3 = t34 * t179 + (t18 * t201 + t19 * t198) * t180;
t35 = t141 * t93 + t142 * t94 + t92 * t266;
t4 = t35 * t179 + (t198 * t21 + t20 * t201) * t180;
t246 = t4 * t287 + t3 * t288 + t179 * (-t277 + t278) / 0.2e1 + t9 * t266 / 0.2e1 + t8 * t265 / 0.2e1;
t99 = t179 * rSges(7,3) + (-rSges(7,1) * t196 - rSges(7,2) * t199) * t180;
t245 = -pkin(10) * t179 - t152 - t99;
t244 = t258 - t285;
t243 = -t285 - t286;
t242 = -t281 + t283;
t239 = -rSges(7,1) * t142 - rSges(7,2) * t141;
t238 = -t196 * t94 - t199 * t93;
t237 = Icges(3,1) * t200 - t276;
t234 = -Icges(3,2) * t197 + t275;
t230 = Icges(3,5) * t200 - Icges(3,6) * t197;
t158 = Icges(4,2) * t183 + t274;
t159 = Icges(4,1) * t182 + t273;
t215 = -t158 * t182 + t159 * t183;
t25 = t41 + t284;
t211 = -t154 + t243;
t210 = t247 + t257;
t77 = rSges(7,3) * t266 - t239;
t22 = t259 + t294 * t201 + (-pkin(5) * t201 + pkin(10) * t266 + t77) * t198;
t209 = t245 - t285;
t208 = t243 + t258;
t207 = t251 * t201 + t252;
t12 = t22 + t284;
t206 = t201 * t291 + t250;
t205 = t209 - t286;
t204 = t278 / 0.2e1 - t277 / 0.2e1 + (t303 * t179 + t305 * t180 + t293 * t198 + t292 * t201 + t34) * t288 + (t304 * t179 + t306 * t180 + t292 * t198 - t293 * t201 + t35) * t287;
t157 = Icges(4,5) * t182 + Icges(4,6) * t183;
t203 = t204 + (t128 * t183 + t130 * t182 + t198 * t157 + t215 * t201) * t288 + (t127 * t183 + t129 * t182 - t201 * t157 + t215 * t198) * t287;
t174 = rSges(2,1) * t201 - rSges(2,2) * t198;
t173 = -rSges(2,1) * t198 - rSges(2,2) * t201;
t172 = rSges(3,1) * t197 + rSges(3,2) * t200;
t134 = Icges(3,3) * t198 + t230 * t201;
t133 = -Icges(3,3) * t201 + t230 * t198;
t124 = t249 * t201;
t123 = t249 * t198;
t103 = t295 + (pkin(1) - t281) * t201 + t256;
t102 = t279 + t191 + (-pkin(1) - t242) * t198;
t101 = t248 * t201;
t100 = t248 * t198;
t91 = t211 * t201;
t90 = t211 * t198;
t89 = t179 * t92;
t88 = -t198 * t202 + t175 + t214;
t87 = (rSges(4,3) - t202) * t201 + (-t181 - t241) * t198;
t86 = t258 * t201;
t85 = t258 * t198;
t82 = t201 * (-t201 * t281 + t256) + (t242 * t198 - t279) * t198;
t81 = t213 + t247;
t80 = (rSges(5,3) - t192) * t201 + (-t161 - t240) * t198;
t79 = t244 * t201;
t78 = t244 * t198;
t66 = t208 * t201;
t65 = t208 * t198;
t63 = t245 * t201;
t62 = t245 * t198;
t57 = t210 + t212;
t56 = (rSges(6,1) - t192) * t201 + (-t161 + (rSges(6,2) - pkin(4)) * t180 + (-rSges(6,3) - qJ(5)) * t179) * t198;
t51 = t209 * t201;
t50 = t209 * t198;
t45 = t205 * t201;
t44 = t205 * t198;
t43 = t179 * t76 - t99 * t265;
t42 = -t179 * t77 + t99 * t266;
t40 = t210 + t294;
t39 = (pkin(5) - t192) * t201 + (-t268 - t161 + (-rSges(7,3) - pkin(4) - pkin(10)) * t180) * t198 + t239;
t38 = t69 + t260;
t37 = (t238 * t180 + t89) * t179;
t36 = (-t198 * t76 + t201 * t77) * t180;
t23 = t33 + t260;
t13 = t25 + t260;
t11 = t12 + t260;
t1 = [t183 * t158 + t182 * t159 + t200 * (Icges(3,2) * t200 + t276) + t197 * (Icges(3,1) * t197 + t275) + Icges(2,3) + t89 + t300 * t179 + (t238 - t301) * t180 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(6) * (t56 ^ 2 + t57 ^ 2) + m(4) * (t87 ^ 2 + t88 ^ 2) + m(3) * (t102 ^ 2 + t103 ^ 2) + m(2) * (t173 ^ 2 + t174 ^ 2); m(7) * (t39 * t45 + t40 * t44) + m(6) * (t56 * t66 + t57 * t65) + m(5) * (t80 * t91 + t81 * t90) + m(4) * (t123 * t88 + t124 * t87) + m(3) * (-t102 * t201 - t103 * t198) * t172 + t203 + (t194 / 0.2e1 + t193 / 0.2e1) * (Icges(3,5) * t197 + Icges(3,6) * t200) + (t200 * (Icges(3,6) * t198 + t234 * t201) + t197 * (Icges(3,5) * t198 + t237 * t201)) * t288 + (t200 * (-Icges(3,6) * t201 + t234 * t198) + t197 * (-Icges(3,5) * t201 + t237 * t198)) * t287; m(7) * (t11 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t13 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t23 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(4) * (t123 ^ 2 + t124 ^ 2 + t38 ^ 2) + t198 * t193 * t134 + m(3) * (t253 * t172 ^ 2 + t82 ^ 2) + t250 + (-t194 * t133 + (-t198 * t133 + t201 * t134) * t198 + t291) * t201; m(4) * (-t198 * t88 - t201 * t87) * t160 + t203 + m(7) * (t39 * t51 + t40 * t50) + m(5) * (t100 * t81 + t101 * t80) + m(6) * (t56 * t79 + t57 * t78); m(7) * (t11 * t12 + t44 * t50 + t45 * t51) + m(6) * (t13 * t25 + t65 * t78 + t66 * t79) + m(5) * (t100 * t90 + t101 * t91 + t23 * t33) + m(4) * (t69 * t38 + (-t123 * t198 - t124 * t201) * t160) + t206; m(7) * (t12 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t25 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t33 ^ 2) + m(4) * (t253 * t160 ^ 2 + t69 ^ 2) + t206; m(7) * (t39 * t63 + t40 * t62) + m(6) * (t56 * t86 + t57 * t85) + m(5) * (-t198 * t81 - t201 * t80) * t154 + t204; m(7) * (t11 * t22 + t44 * t62 + t45 * t63) + m(6) * (t13 * t41 + t65 * t85 + t66 * t86) + m(5) * (t64 * t23 + (-t198 * t90 - t201 * t91) * t154) + t207; m(7) * (t12 * t22 + t50 * t62 + t51 * t63) + m(6) * (t25 * t41 + t78 * t85 + t79 * t86) + m(5) * (t64 * t33 + (-t100 * t198 - t101 * t201) * t154) + t207; m(7) * (t22 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t41 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(5) * (t253 * t154 ^ 2 + t64 ^ 2) + t207; 0.2e1 * ((t198 * t40 + t201 * t39) * t289 + (t198 * t57 + t201 * t56) * t290) * t179; m(7) * (-t180 * t11 + (t198 * t44 + t201 * t45) * t179) + m(6) * (-t180 * t13 + (t198 * t65 + t201 * t66) * t179); m(7) * (-t180 * t12 + (t198 * t50 + t201 * t51) * t179) + m(6) * (-t180 * t25 + (t198 * t78 + t201 * t79) * t179); m(7) * (-t180 * t22 + (t198 * t62 + t201 * t63) * t179) + m(6) * (-t180 * t41 + (t198 * t85 + t201 * t86) * t179); 0.2e1 * (t290 + t289) * (t253 * t179 ^ 2 + t180 ^ 2); m(7) * (t39 * t42 + t40 * t43) + t37 + ((t34 / 0.2e1 + t29 / 0.2e1) * t201 + (t35 / 0.2e1 + t30 / 0.2e1) * t198) * t180; m(7) * (t11 * t36 + t42 * t45 + t43 * t44) + t246; m(7) * (t12 * t36 + t42 * t51 + t43 * t50) + t246; m(7) * (t22 * t36 + t42 * t63 + t43 * t62) + t246; m(7) * (-t36 * t180 + (t198 * t43 + t201 * t42) * t179); t179 * t37 + m(7) * (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) + (t201 * t3 + t198 * t4 + t179 * (t198 * t30 + t201 * t29)) * t180;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
