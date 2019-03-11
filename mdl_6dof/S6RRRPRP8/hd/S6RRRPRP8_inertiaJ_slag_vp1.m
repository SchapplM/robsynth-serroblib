% Calculate joint inertia matrix for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:17
% EndTime: 2019-03-09 17:14:30
% DurationCPUTime: 5.42s
% Computational Cost: add. (9235->552), mult. (22571->780), div. (0->0), fcn. (26955->8), ass. (0->260)
t339 = rSges(7,3) + qJ(6) + pkin(9);
t235 = sin(qJ(5));
t236 = sin(qJ(3));
t237 = sin(qJ(2));
t239 = cos(qJ(5));
t240 = cos(qJ(3));
t192 = (-t235 * t240 + t236 * t239) * t237;
t308 = t235 * t236;
t193 = (t239 * t240 + t308) * t237;
t241 = cos(qJ(2));
t122 = Icges(7,5) * t193 + Icges(7,6) * t192 + Icges(7,3) * t241;
t123 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t241;
t124 = Icges(7,4) * t193 + Icges(7,2) * t192 + Icges(7,6) * t241;
t125 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t241;
t126 = Icges(7,1) * t193 + Icges(7,4) * t192 + Icges(7,5) * t241;
t127 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t241;
t342 = (t122 + t123) * t241 + (t126 + t127) * t193 + (t124 + t125) * t192;
t242 = cos(qJ(1));
t300 = t242 * t240;
t238 = sin(qJ(1));
t303 = t238 * t241;
t204 = t236 * t303 + t300;
t301 = t242 * t236;
t205 = t240 * t303 - t301;
t156 = t204 * t239 - t205 * t235;
t310 = t204 * t235;
t157 = t205 * t239 + t310;
t306 = t237 * t238;
t89 = Icges(7,5) * t157 + Icges(7,6) * t156 - Icges(7,3) * t306;
t93 = Icges(7,4) * t157 + Icges(7,2) * t156 - Icges(7,6) * t306;
t97 = Icges(7,1) * t157 + Icges(7,4) * t156 - Icges(7,5) * t306;
t36 = t192 * t93 + t193 * t97 + t241 * t89;
t91 = Icges(6,5) * t157 + Icges(6,6) * t156 - Icges(6,3) * t306;
t95 = Icges(6,4) * t157 + Icges(6,2) * t156 - Icges(6,6) * t306;
t99 = Icges(6,1) * t157 + Icges(6,4) * t156 - Icges(6,5) * t306;
t38 = t192 * t95 + t193 * t99 + t241 * t91;
t341 = t36 + t38;
t206 = -t238 * t240 + t241 * t301;
t207 = t238 * t236 + t241 * t300;
t158 = t206 * t239 - t207 * t235;
t309 = t206 * t235;
t159 = t207 * t239 + t309;
t304 = t237 * t242;
t90 = Icges(7,5) * t159 + Icges(7,6) * t158 - Icges(7,3) * t304;
t94 = Icges(7,4) * t159 + Icges(7,2) * t158 - Icges(7,6) * t304;
t98 = Icges(7,1) * t159 + Icges(7,4) * t158 - Icges(7,5) * t304;
t37 = t192 * t94 + t193 * t98 + t241 * t90;
t100 = Icges(6,1) * t159 + Icges(6,4) * t158 - Icges(6,5) * t304;
t92 = Icges(6,5) * t159 + Icges(6,6) * t158 - Icges(6,3) * t304;
t96 = Icges(6,4) * t159 + Icges(6,2) * t158 - Icges(6,6) * t304;
t39 = t193 * t100 + t192 * t96 + t241 * t92;
t340 = t37 + t39;
t178 = -Icges(4,3) * t241 + (Icges(4,5) * t240 - Icges(4,6) * t236) * t237;
t181 = -Icges(5,2) * t241 + (Icges(5,4) * t240 + Icges(5,6) * t236) * t237;
t338 = t178 + t181;
t316 = t342 * t241;
t321 = -t316 + (t341 * t238 + t340 * t242) * t237;
t27 = t158 * t93 + t159 * t97 - t89 * t304;
t28 = t158 * t94 + t159 * t98 - t90 * t304;
t50 = -t122 * t304 + t158 * t124 + t159 * t126;
t3 = -t50 * t241 + (t238 * t27 + t242 * t28) * t237;
t29 = t158 * t95 + t159 * t99 - t91 * t304;
t30 = t159 * t100 + t158 * t96 - t92 * t304;
t51 = -t123 * t304 + t158 * t125 + t159 * t127;
t4 = -t51 * t241 + (t238 * t29 + t242 * t30) * t237;
t322 = t3 + t4;
t23 = t156 * t93 + t157 * t97 - t89 * t306;
t24 = t156 * t94 + t157 * t98 - t90 * t306;
t48 = -t122 * t306 + t156 * t124 + t157 * t126;
t1 = -t48 * t241 + (t23 * t238 + t24 * t242) * t237;
t25 = t156 * t95 + t157 * t99 - t91 * t306;
t26 = t157 * t100 + t156 * t96 - t92 * t306;
t49 = -t123 * t306 + t156 * t125 + t157 * t127;
t2 = -t49 * t241 + (t238 * t25 + t242 * t26) * t237;
t323 = t1 + t2;
t337 = (t323 * t238 + t322 * t242) * t237 - t321 * t241;
t336 = t237 / 0.2e1;
t334 = t241 / 0.2e1;
t332 = (t340 * t238 / 0.2e1 - t341 * t242 / 0.2e1) * t241;
t182 = -Icges(4,6) * t241 + (Icges(4,4) * t240 - Icges(4,2) * t236) * t237;
t307 = t236 * t237;
t177 = -Icges(5,6) * t241 + (Icges(5,5) * t240 + Icges(5,3) * t236) * t237;
t185 = -Icges(5,4) * t241 + (Icges(5,1) * t240 + Icges(5,5) * t236) * t237;
t186 = -Icges(4,5) * t241 + (Icges(4,1) * t240 - Icges(4,4) * t236) * t237;
t305 = t237 * t240;
t326 = t177 * t307 + (t185 + t186) * t305;
t331 = (t182 * t307 + t338 * t241 - t326) * t241;
t222 = pkin(9) * t306;
t226 = t239 * pkin(5) + pkin(4);
t317 = -pkin(4) + t226;
t325 = -t157 * rSges(7,1) - t156 * rSges(7,2) - pkin(5) * t310;
t299 = t317 * t205 - t339 * t306 + t222 - t325;
t330 = t299 * t241;
t329 = t299 * t242;
t327 = t159 * rSges(7,1) + t158 * rSges(7,2) + pkin(5) * t309 + t207 * t226 - t339 * t304;
t232 = t238 ^ 2;
t233 = t242 ^ 2;
t320 = Icges(3,5) * t336 + Icges(3,6) * t334;
t319 = m(7) * t237;
t318 = pkin(2) * t241;
t315 = t204 * rSges(5,3);
t314 = t242 * rSges(3,3);
t313 = Icges(3,4) * t237;
t312 = Icges(3,4) * t241;
t257 = -t157 * rSges(6,1) - t156 * rSges(6,2);
t102 = -rSges(6,3) * t306 - t257;
t311 = t102 * t242;
t302 = t241 * t242;
t202 = t207 * pkin(4);
t174 = -pkin(9) * t304 + t202;
t298 = -t174 + t327;
t296 = t193 * rSges(7,1) + t192 * rSges(7,2) + (pkin(5) * t308 + t317 * t240) * t237 + (-pkin(9) + t339) * t241;
t144 = t207 * rSges(5,1) + rSges(5,2) * t304 + t206 * rSges(5,3);
t163 = t207 * pkin(3) + t206 * qJ(4);
t295 = -t144 - t163;
t197 = t204 * qJ(4);
t162 = t205 * pkin(3) + t197;
t146 = t162 * t304;
t173 = t205 * pkin(4) - t222;
t294 = t173 * t304 + t146;
t293 = t159 * rSges(6,1) + t158 * rSges(6,2);
t208 = (pkin(3) * t240 + qJ(4) * t236) * t237;
t292 = t241 * t162 + t208 * t306;
t291 = -t163 - t174;
t190 = -t241 * rSges(5,2) + (rSges(5,1) * t240 + rSges(5,3) * t236) * t237;
t289 = -t190 - t208;
t191 = -t241 * rSges(4,3) + (rSges(4,1) * t240 - rSges(4,2) * t236) * t237;
t217 = t237 * pkin(2) - t241 * pkin(8);
t288 = -t191 - t217;
t285 = pkin(2) * t302 + pkin(8) * t304;
t287 = t232 * (pkin(8) * t237 + t318) + t242 * t285;
t209 = pkin(4) * t305 + t241 * pkin(9);
t286 = -t208 - t209;
t284 = t242 * pkin(1) + t238 * pkin(7);
t229 = t242 * pkin(7);
t283 = t229 - t197;
t282 = t232 + t233;
t281 = t2 / 0.2e1 + t1 / 0.2e1;
t280 = -t3 / 0.2e1 - t4 / 0.2e1;
t7 = -t23 * t242 + t24 * t238;
t8 = t26 * t238 - t25 * t242;
t279 = -t7 / 0.2e1 - t8 / 0.2e1;
t10 = t30 * t238 - t29 * t242;
t9 = t28 * t238 - t27 * t242;
t277 = -t10 / 0.2e1 - t9 / 0.2e1;
t130 = Icges(5,5) * t205 + Icges(5,6) * t306 + Icges(5,3) * t204;
t134 = Icges(5,4) * t205 + Icges(5,2) * t306 + Icges(5,6) * t204;
t138 = Icges(5,1) * t205 + Icges(5,4) * t306 + Icges(5,5) * t204;
t69 = -t241 * t134 + (t130 * t236 + t138 * t240) * t237;
t132 = Icges(4,5) * t205 - Icges(4,6) * t204 + Icges(4,3) * t306;
t136 = Icges(4,4) * t205 - Icges(4,2) * t204 + Icges(4,6) * t306;
t140 = Icges(4,1) * t205 - Icges(4,4) * t204 + Icges(4,5) * t306;
t71 = -t241 * t132 + (-t136 * t236 + t140 * t240) * t237;
t275 = t69 / 0.2e1 + t71 / 0.2e1;
t131 = Icges(5,5) * t207 + Icges(5,6) * t304 + Icges(5,3) * t206;
t135 = Icges(5,4) * t207 + Icges(5,2) * t304 + Icges(5,6) * t206;
t139 = Icges(5,1) * t207 + Icges(5,4) * t304 + Icges(5,5) * t206;
t70 = -t241 * t135 + (t131 * t236 + t139 * t240) * t237;
t133 = Icges(4,5) * t207 - Icges(4,6) * t206 + Icges(4,3) * t304;
t137 = Icges(4,4) * t207 - Icges(4,2) * t206 + Icges(4,6) * t304;
t141 = Icges(4,1) * t207 - Icges(4,4) * t206 + Icges(4,5) * t304;
t72 = -t241 * t133 + (-t137 * t236 + t141 * t240) * t237;
t274 = -t70 / 0.2e1 - t72 / 0.2e1;
t104 = -rSges(6,3) * t304 + t293;
t273 = -t104 + t291;
t129 = t193 * rSges(6,1) + t192 * rSges(6,2) + t241 * rSges(6,3);
t270 = -t129 + t286;
t268 = -t217 + t289;
t145 = t207 * rSges(4,1) - t206 * rSges(4,2) + rSges(4,3) * t304;
t267 = -pkin(1) - t318;
t266 = t237 * t296;
t265 = t291 - t298;
t264 = t286 - t296;
t263 = -t217 + t270;
t262 = t238 * t162 + t242 * t163 + t287;
t261 = t241 * t173 + t209 * t306 + t292;
t260 = t284 + t285;
t259 = rSges(3,1) * t241 - rSges(3,2) * t237;
t258 = -t205 * rSges(4,1) + t204 * rSges(4,2);
t255 = -t217 + t264;
t254 = Icges(3,1) * t241 - t313;
t253 = -Icges(3,2) * t237 + t312;
t252 = Icges(3,5) * t241 - Icges(3,6) * t237;
t249 = rSges(3,1) * t302 - rSges(3,2) * t304 + t238 * rSges(3,3);
t248 = -t38 / 0.2e1 - t36 / 0.2e1 - t49 / 0.2e1 - t48 / 0.2e1;
t247 = -t39 / 0.2e1 - t37 / 0.2e1 - t51 / 0.2e1 - t50 / 0.2e1;
t246 = t238 * t173 + t242 * t174 + t262;
t73 = -t241 * t102 - t129 * t306;
t245 = t163 + t260;
t82 = t204 * t177 + t181 * t306 + t205 * t185;
t83 = t178 * t306 - t204 * t182 + t205 * t186;
t244 = t83 / 0.2e1 + t82 / 0.2e1 - t248 + t275;
t84 = t206 * t177 + t181 * t304 + t207 * t185;
t85 = t178 * t304 - t206 * t182 + t207 * t186;
t243 = t85 / 0.2e1 + t84 / 0.2e1 - t247 - t274;
t231 = t237 ^ 2;
t216 = t242 * rSges(2,1) - t238 * rSges(2,2);
t215 = -t238 * rSges(2,1) - t242 * rSges(2,2);
t214 = t237 * rSges(3,1) + t241 * rSges(3,2);
t180 = Icges(3,3) * t238 + t252 * t242;
t179 = -Icges(3,3) * t242 + t252 * t238;
t169 = t249 + t284;
t168 = t314 + t229 + (-pkin(1) - t259) * t238;
t161 = t288 * t242;
t160 = t288 * t238;
t143 = rSges(4,3) * t306 - t258;
t142 = t205 * rSges(5,1) + rSges(5,2) * t306 + t315;
t121 = t242 * t249 + (t259 * t238 - t314) * t238;
t118 = t268 * t242;
t117 = t268 * t238;
t112 = t260 + t145;
t111 = t229 + ((-rSges(4,3) - pkin(8)) * t237 + t267) * t238 + t258;
t110 = -t241 * t145 - t191 * t304;
t109 = t241 * t143 + t191 * t306;
t88 = t263 * t242;
t87 = t263 * t238;
t86 = (t143 * t242 - t145 * t238) * t237;
t81 = t245 + t144;
t80 = -t315 + (-rSges(5,1) - pkin(3)) * t205 + ((-rSges(5,2) - pkin(8)) * t237 + t267) * t238 + t283;
t79 = t238 * t143 + t242 * t145 + t287;
t78 = t295 * t241 + t289 * t304;
t77 = t241 * t142 + t190 * t306 + t292;
t76 = t255 * t242;
t75 = t255 * t238;
t74 = t241 * t104 + t129 * t304;
t68 = t133 * t304 - t206 * t137 + t207 * t141;
t67 = t132 * t304 - t206 * t136 + t207 * t140;
t66 = t206 * t131 + t135 * t304 + t207 * t139;
t65 = t206 * t130 + t134 * t304 + t207 * t138;
t64 = t133 * t306 - t204 * t137 + t205 * t141;
t63 = t132 * t306 - t204 * t136 + t205 * t140;
t62 = t204 * t131 + t135 * t306 + t205 * t139;
t61 = t204 * t130 + t134 * t306 + t205 * t138;
t60 = t202 + (-rSges(6,3) - pkin(9)) * t304 + t245 + t293;
t59 = t222 + (-pkin(3) - pkin(4)) * t205 + ((rSges(6,3) - pkin(8)) * t237 + t267) * t238 + t257 + t283;
t58 = t146 + (t142 * t242 + t295 * t238) * t237;
t55 = (t104 * t238 - t311) * t237;
t54 = t245 + t327;
t53 = (-pkin(3) - t226) * t205 + ((-pkin(8) + t339) * t237 + t267) * t238 + t283 + t325;
t52 = t238 * t142 + t242 * t144 + t262;
t43 = t273 * t241 + t270 * t304;
t42 = -t73 + t261;
t41 = t298 * t241 + t242 * t266;
t40 = -t296 * t306 - t330;
t35 = (t273 * t238 + t311) * t237 + t294;
t34 = t68 * t238 - t67 * t242;
t33 = t66 * t238 - t65 * t242;
t32 = t64 * t238 - t63 * t242;
t31 = t62 * t238 - t61 * t242;
t22 = (t298 * t238 - t329) * t237;
t21 = t238 * t102 + t242 * t104 + t246;
t20 = t265 * t241 + t264 * t304;
t19 = t238 * t266 + t261 + t330;
t18 = -t85 * t241 + (t238 * t67 + t242 * t68) * t237;
t17 = -t84 * t241 + (t238 * t65 + t242 * t66) * t237;
t16 = -t83 * t241 + (t238 * t63 + t242 * t64) * t237;
t15 = -t82 * t241 + (t238 * t61 + t242 * t62) * t237;
t14 = (t265 * t238 + t329) * t237 + t294;
t13 = t299 * t238 + t298 * t242 + t246;
t5 = [Icges(2,3) + (Icges(3,1) * t237 - t236 * t182 + t312) * t237 + (Icges(3,2) * t241 + t313 - t338) * t241 + m(7) * (t53 ^ 2 + t54 ^ 2) + m(6) * (t59 ^ 2 + t60 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2) + m(3) * (t168 ^ 2 + t169 ^ 2) + m(2) * (t215 ^ 2 + t216 ^ 2) + t326 + t342; (-t241 * (-Icges(3,6) * t242 + t253 * t238) / 0.2e1 - t237 * (-Icges(3,5) * t242 + t254 * t238) / 0.2e1 + t242 * t320 - t244) * t242 + ((Icges(3,6) * t238 + t253 * t242) * t334 + (Icges(3,5) * t238 + t254 * t242) * t336 + t238 * t320 + t243) * t238 + m(7) * (t53 * t76 + t54 * t75) + m(6) * (t59 * t88 + t60 * t87) + m(5) * (t117 * t81 + t118 * t80) + m(4) * (t111 * t161 + t112 * t160) + m(3) * (-t168 * t242 - t169 * t238) * t214; m(4) * (t160 ^ 2 + t161 ^ 2 + t79 ^ 2) + m(3) * (t282 * t214 ^ 2 + t121 ^ 2) + m(7) * (t13 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t21 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t117 ^ 2 + t118 ^ 2 + t52 ^ 2) + (-t233 * t179 - t31 - t32 - t7 - t8) * t242 + (t232 * t180 + t10 + t33 + t34 + t9 + (-t238 * t179 + t242 * t180) * t242) * t238; t331 + m(7) * (t19 * t53 + t20 * t54) + m(6) * (t42 * t59 + t43 * t60) + m(5) * (t77 * t80 + t78 * t81) + m(4) * (t109 * t111 + t110 * t112) + (t238 * t244 + t242 * t243) * t237 - t316; -t332 + (-t15 / 0.2e1 - t16 / 0.2e1 + t275 * t241 - t281) * t242 + (t18 / 0.2e1 + t17 / 0.2e1 + t274 * t241 - t280) * t238 + m(7) * (t14 * t13 + t19 * t76 + t20 * t75) + m(6) * (t35 * t21 + t42 * t88 + t43 * t87) + m(5) * (t117 * t78 + t118 * t77 + t58 * t52) + m(4) * (t109 * t161 + t110 * t160 + t79 * t86) + ((t33 / 0.2e1 + t34 / 0.2e1 - t277) * t242 + (t32 / 0.2e1 + t31 / 0.2e1 - t279) * t238) * t237; (-t321 - t331) * t241 + m(6) * (t35 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t58 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t109 ^ 2 + t110 ^ 2 + t86 ^ 2) + m(7) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + ((t17 + t18 + (-t70 - t72) * t241 + t322) * t242 + (t15 + t16 + (-t69 - t71) * t241 + t323) * t238) * t237; m(7) * (t204 * t54 + t206 * t53) + m(6) * (t204 * t60 + t206 * t59) + m(5) * (t204 * t81 + t206 * t80); m(7) * (t13 * t307 + t204 * t75 + t206 * t76) + m(6) * (t204 * t87 + t206 * t88 + t21 * t307) + m(5) * (t204 * t117 + t206 * t118 + t52 * t307); m(7) * (t14 * t307 + t206 * t19 + t204 * t20) + m(6) * (t204 * t43 + t206 * t42 + t35 * t307) + m(5) * (t204 * t78 + t206 * t77 + t58 * t307); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t231 * t236 ^ 2 + t204 ^ 2 + t206 ^ 2); m(7) * (t40 * t53 + t41 * t54) + m(6) * (t59 * t73 + t60 * t74) + (t248 * t238 + t247 * t242) * t237 + t316; t281 * t242 + t332 + t280 * t238 + m(7) * (t22 * t13 + t40 * t76 + t41 * t75) + m(6) * (t55 * t21 + t73 * t88 + t74 * t87) + (t279 * t238 + t277 * t242) * t237; m(7) * (t14 * t22 + t19 * t40 + t20 * t41) + m(6) * (t55 * t35 + t42 * t73 + t43 * t74) - t337; m(6) * (t74 * t204 + t73 * t206 + t307 * t55) + m(7) * (t41 * t204 + t40 * t206 + t22 * t307); m(7) * (t22 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t55 ^ 2 + t73 ^ 2 + t74 ^ 2) + t337; (-t238 * t54 - t242 * t53) * t319; m(7) * (t241 * t13 + (-t238 * t75 - t242 * t76) * t237); m(7) * (t241 * t14 + (-t19 * t242 - t20 * t238) * t237); (-t204 * t238 - t206 * t242 + t236 * t241) * t319; m(7) * (t241 * t22 + (-t238 * t41 - t242 * t40) * t237); m(7) * (t231 * t282 + t241 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
