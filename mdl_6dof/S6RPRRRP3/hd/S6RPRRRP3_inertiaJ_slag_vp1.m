% Calculate joint inertia matrix for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:03:01
% EndTime: 2019-03-09 06:03:09
% DurationCPUTime: 4.03s
% Computational Cost: add. (11750->405), mult. (10804->569), div. (0->0), fcn. (11796->10), ass. (0->215)
t211 = qJ(1) + pkin(10);
t206 = sin(t211);
t207 = cos(t211);
t212 = qJ(4) + qJ(5);
t209 = cos(t212);
t208 = sin(t212);
t217 = cos(qJ(3));
t264 = t208 * t217;
t157 = t206 * t264 + t207 * t209;
t262 = t209 * t217;
t158 = t206 * t262 - t207 * t208;
t214 = sin(qJ(3));
t269 = t206 * t214;
t101 = Icges(7,4) * t158 + Icges(7,2) * t269 + Icges(7,6) * t157;
t99 = Icges(6,5) * t158 - Icges(6,6) * t157 + Icges(6,3) * t269;
t316 = t101 + t99;
t103 = Icges(6,4) * t158 - Icges(6,2) * t157 + Icges(6,6) * t269;
t97 = Icges(7,5) * t158 + Icges(7,6) * t269 + Icges(7,3) * t157;
t315 = -t103 + t97;
t159 = -t206 * t209 + t207 * t264;
t160 = t206 * t208 + t207 * t262;
t267 = t207 * t214;
t104 = Icges(6,4) * t160 - Icges(6,2) * t159 + Icges(6,6) * t267;
t98 = Icges(7,5) * t160 + Icges(7,6) * t267 + Icges(7,3) * t159;
t314 = -t104 + t98;
t100 = Icges(6,5) * t160 - Icges(6,6) * t159 + Icges(6,3) * t267;
t102 = Icges(7,4) * t160 + Icges(7,2) * t267 + Icges(7,6) * t159;
t313 = t100 + t102;
t105 = Icges(7,1) * t158 + Icges(7,4) * t269 + Icges(7,5) * t157;
t107 = Icges(6,1) * t158 - Icges(6,4) * t157 + Icges(6,5) * t269;
t312 = t105 + t107;
t106 = Icges(7,1) * t160 + Icges(7,4) * t267 + Icges(7,5) * t159;
t108 = Icges(6,1) * t160 - Icges(6,4) * t159 + Icges(6,5) * t267;
t311 = t106 + t108;
t290 = rSges(7,3) + qJ(6);
t298 = rSges(7,1) + pkin(5);
t310 = -t290 * t157 - t298 * t158;
t309 = t315 * t157 + t312 * t158 + t316 * t269;
t308 = t314 * t157 + t311 * t158 + t313 * t269;
t307 = t315 * t159 + t312 * t160 + t316 * t267;
t306 = t314 * t159 + t311 * t160 + t313 * t267;
t161 = -Icges(7,6) * t217 + (Icges(7,5) * t209 + Icges(7,3) * t208) * t214;
t163 = -Icges(7,2) * t217 + (Icges(7,4) * t209 + Icges(7,6) * t208) * t214;
t165 = -Icges(7,4) * t217 + (Icges(7,1) * t209 + Icges(7,5) * t208) * t214;
t70 = t157 * t161 + t158 * t165 + t163 * t269;
t162 = -Icges(6,3) * t217 + (Icges(6,5) * t209 - Icges(6,6) * t208) * t214;
t164 = -Icges(6,6) * t217 + (Icges(6,4) * t209 - Icges(6,2) * t208) * t214;
t166 = -Icges(6,5) * t217 + (Icges(6,1) * t209 - Icges(6,4) * t208) * t214;
t71 = -t157 * t164 + t158 * t166 + t162 * t269;
t305 = -t71 - t70;
t72 = t159 * t161 + t160 * t165 + t163 * t267;
t73 = -t159 * t164 + t160 * t166 + t162 * t267;
t304 = -t72 - t73;
t265 = t208 * t214;
t300 = t161 * t265 + (t165 + t166) * t209 * t214;
t302 = -t162 - t163;
t291 = -t164 * t265 + t217 * t302 + t300;
t303 = t291 * t217;
t301 = Icges(4,5) * t214;
t299 = t301 / 0.2e1;
t297 = t305 * t217 + (t206 * t309 + t207 * t308) * t214;
t296 = t304 * t217 + (t206 * t307 + t207 * t306) * t214;
t295 = t206 * t308 - t207 * t309;
t294 = t206 * t306 - t207 * t307;
t52 = -t101 * t217 + (t105 * t209 + t208 * t97) * t214;
t54 = -t217 * t99 + (-t103 * t208 + t107 * t209) * t214;
t293 = -t52 - t54;
t53 = -t102 * t217 + (t106 * t209 + t208 * t98) * t214;
t55 = -t100 * t217 + (-t104 * t208 + t108 * t209) * t214;
t292 = t53 + t55;
t289 = rSges(7,2) * t269 - t310;
t288 = -rSges(7,2) * t217 + (t290 * t208 + t209 * t298) * t214;
t287 = t206 ^ 2;
t286 = t207 ^ 2;
t285 = t206 / 0.2e1;
t284 = -t207 / 0.2e1;
t283 = t207 / 0.2e1;
t282 = -t217 / 0.2e1;
t189 = rSges(4,1) * t214 + rSges(4,2) * t217;
t281 = m(4) * t189;
t215 = sin(qJ(1));
t280 = pkin(1) * t215;
t279 = pkin(3) * t217;
t278 = pkin(8) * t214;
t216 = cos(qJ(4));
t205 = pkin(4) * t216 + pkin(3);
t277 = -pkin(3) + t205;
t276 = t303 + (t293 * t206 - t292 * t207) * t214;
t274 = t207 * rSges(4,3);
t273 = t289 * t267;
t232 = -rSges(6,1) * t158 + rSges(6,2) * t157;
t111 = rSges(6,3) * t269 - t232;
t168 = -rSges(6,3) * t217 + (rSges(6,1) * t209 - rSges(6,2) * t208) * t214;
t81 = t217 * t111 + t168 * t269;
t271 = Icges(4,4) * t217;
t213 = sin(qJ(4));
t270 = t206 * t213;
t268 = t207 * t213;
t266 = t207 * t217;
t176 = -Icges(5,6) * t217 + (Icges(5,4) * t216 - Icges(5,2) * t213) * t214;
t261 = t213 * t176;
t260 = t213 * t217;
t219 = -pkin(9) - pkin(8);
t259 = t214 * t219;
t258 = t216 * t217;
t257 = rSges(7,2) * t267 + t290 * t159 + t298 * t160;
t113 = t160 * rSges(6,1) - t159 * rSges(6,2) + rSges(6,3) * t267;
t224 = pkin(4) * t270 + t205 * t266 - t207 * t259;
t248 = pkin(3) * t266 + pkin(8) * t267;
t128 = t224 - t248;
t256 = -t113 - t128;
t249 = -pkin(4) * t268 - t206 * t259;
t127 = (t217 * t277 - t278) * t206 + t249;
t156 = (pkin(8) + t219) * t217 + t277 * t214;
t255 = t217 * t127 + t156 * t269;
t253 = t287 * (t278 + t279) + t207 * t248;
t252 = -t156 - t168;
t178 = -rSges(5,3) * t217 + (rSges(5,1) * t216 - rSges(5,2) * t213) * t214;
t196 = t214 * pkin(3) - t217 * pkin(8);
t250 = -t178 - t196;
t171 = -t206 * t260 - t207 * t216;
t172 = t206 * t258 - t268;
t119 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t269;
t121 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t269;
t123 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t269;
t60 = -t119 * t217 + (-t121 * t213 + t123 * t216) * t214;
t175 = -Icges(5,3) * t217 + (Icges(5,5) * t216 - Icges(5,6) * t213) * t214;
t177 = -Icges(5,5) * t217 + (Icges(5,1) * t216 - Icges(5,4) * t213) * t214;
t77 = t171 * t176 + t172 * t177 + t175 * t269;
t247 = t60 / 0.2e1 + t77 / 0.2e1;
t173 = t206 * t216 - t207 * t260;
t174 = t207 * t258 + t270;
t120 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t267;
t122 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t267;
t124 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t267;
t61 = -t120 * t217 + (-t122 * t213 + t124 * t216) * t214;
t78 = t173 * t176 + t174 * t177 + t175 * t267;
t246 = t78 / 0.2e1 + t61 / 0.2e1;
t245 = t296 * t267 + t297 * t269;
t244 = -t128 - t257;
t243 = -t156 - t288;
t242 = -t196 + t252;
t126 = t174 * rSges(5,1) + t173 * rSges(5,2) + rSges(5,3) * t267;
t218 = cos(qJ(1));
t210 = t218 * pkin(1);
t241 = t207 * pkin(2) + t206 * pkin(7) + t210;
t240 = t269 / 0.2e1;
t239 = t267 / 0.2e1;
t238 = t207 * pkin(7) - t280;
t237 = -t205 * t217 - pkin(2);
t58 = t289 * t217 + t288 * t269;
t236 = t206 * t127 + t207 * t128 + t253;
t235 = -t196 + t243;
t234 = rSges(4,1) * t217 - rSges(4,2) * t214;
t233 = -rSges(5,1) * t172 - rSges(5,2) * t171;
t230 = -Icges(4,2) * t214 + t271;
t229 = Icges(4,5) * t217 - Icges(4,6) * t214;
t226 = t238 - t249;
t225 = rSges(4,1) * t266 - rSges(4,2) * t267 + t206 * rSges(4,3);
t223 = t276 * t217 + t245;
t222 = (-t293 - t305) * t240 + (t292 - t304) * t239;
t221 = t296 * t285 + t297 * t284 + (t292 * t206 + t293 * t207) * t282 + t295 * t240 + t294 * t239;
t220 = t224 + t241;
t191 = rSges(2,1) * t218 - rSges(2,2) * t215;
t190 = -rSges(2,1) * t215 - rSges(2,2) * t218;
t186 = Icges(4,6) * t217 + t301;
t181 = rSges(3,1) * t207 - rSges(3,2) * t206 + t210;
t180 = -rSges(3,1) * t206 - rSges(3,2) * t207 - t280;
t148 = t214 * t216 * t177;
t143 = Icges(4,3) * t206 + t207 * t229;
t142 = -Icges(4,3) * t207 + t206 * t229;
t134 = t250 * t207;
t133 = t250 * t206;
t132 = t225 + t241;
t131 = t274 + (-pkin(2) - t234) * t206 + t238;
t125 = rSges(5,3) * t269 - t233;
t109 = t127 * t267;
t96 = t207 * t225 + (t206 * t234 - t274) * t206;
t93 = t111 * t267;
t91 = t242 * t207;
t90 = t242 * t206;
t89 = -t217 * t175 - t214 * t261 + t148;
t88 = -t126 * t217 - t178 * t267;
t87 = t125 * t217 + t178 * t269;
t86 = t241 + t126 + t248;
t85 = (-t279 - pkin(2) + (-rSges(5,3) - pkin(8)) * t214) * t206 + t233 + t238;
t82 = -t113 * t217 - t168 * t267;
t80 = t235 * t207;
t79 = t235 * t206;
t76 = t220 + t113;
t75 = (-rSges(6,3) * t214 + t237) * t206 + t226 + t232;
t74 = (t125 * t207 - t126 * t206) * t214;
t65 = -t113 * t269 + t93;
t64 = t220 + t257;
t63 = (-rSges(7,2) * t214 + t237) * t206 + t226 + t310;
t62 = t125 * t206 + t126 * t207 + t253;
t59 = -t217 * t257 - t267 * t288;
t57 = t217 * t256 + t252 * t267;
t56 = t255 + t81;
t51 = t120 * t267 + t122 * t173 + t124 * t174;
t50 = t119 * t267 + t121 * t173 + t123 * t174;
t49 = t120 * t269 + t122 * t171 + t124 * t172;
t48 = t119 * t269 + t121 * t171 + t123 * t172;
t35 = -t257 * t269 + t273;
t34 = t256 * t269 + t109 + t93;
t33 = t217 * t244 + t243 * t267;
t32 = t58 + t255;
t31 = t111 * t206 + t113 * t207 + t236;
t30 = t244 * t269 + t109 + t273;
t29 = t289 * t206 + t257 * t207 + t236;
t26 = t206 * t51 - t207 * t50;
t25 = t206 * t49 - t207 * t48;
t14 = -t217 * t78 + (t206 * t50 + t207 * t51) * t214;
t13 = -t217 * t77 + (t206 * t48 + t207 * t49) * t214;
t1 = [Icges(2,3) + Icges(3,3) + t148 + (Icges(4,1) * t214 - t208 * t164 - t261 + t271) * t214 + (Icges(4,4) * t214 + Icges(4,2) * t217 - t175 + t302) * t217 + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t75 ^ 2 + t76 ^ 2) + m(5) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t131 ^ 2 + t132 ^ 2) + m(2) * (t190 ^ 2 + t191 ^ 2) + m(3) * (t180 ^ 2 + t181 ^ 2) + t300; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t63 * t80 + t64 * t79) + m(6) * (t75 * t91 + t76 * t90) + m(5) * (t133 * t86 + t134 * t85) + (-t54 / 0.2e1 - t52 / 0.2e1 - t71 / 0.2e1 - t70 / 0.2e1 + (-Icges(4,6) * t207 + t206 * t230) * t282 + t207 * t299 - t131 * t281 + t186 * t283 - t247) * t207 + (t73 / 0.2e1 + t72 / 0.2e1 + t55 / 0.2e1 + t53 / 0.2e1 + (Icges(4,6) * t206 + t207 * t230) * t217 / 0.2e1 + t206 * t299 - t132 * t281 + t186 * t285 + t246) * t206; m(4) * t96 + m(5) * t62 + m(6) * t31 + m(7) * t29; m(7) * (t29 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t31 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(5) * (t133 ^ 2 + t134 ^ 2 + t62 ^ 2) + m(4) * (t96 ^ 2 + (t286 + t287) * t189 ^ 2) + (-t286 * t142 - t25 - t295) * t207 + (t287 * t143 + t26 + (-t206 * t142 + t207 * t143) * t207 + t294) * t206; (-t89 - t291) * t217 + m(7) * (t32 * t63 + t33 * t64) + m(6) * (t56 * t75 + t57 * t76) + m(5) * (t85 * t87 + t86 * t88) + (t206 * t247 + t207 * t246) * t214 + t222; m(5) * t74 + m(6) * t34 + m(7) * t30; t221 + (t61 * t206 - t60 * t207) * t282 + m(7) * (t30 * t29 + t32 * t80 + t33 * t79) + m(6) * (t34 * t31 + t56 * t91 + t57 * t90) + m(5) * (t133 * t88 + t134 * t87 + t74 * t62) + (t25 * t285 + t26 * t283) * t214 + t14 * t285 + t13 * t284; (t89 * t217 + t276) * t217 + (t207 * t14 + t206 * t13 - t217 * (t60 * t206 + t61 * t207)) * t214 + m(7) * (t30 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t34 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t74 ^ 2 + t87 ^ 2 + t88 ^ 2) + t245; -t303 + m(7) * (t58 * t63 + t59 * t64) + m(6) * (t75 * t81 + t76 * t82) + t222; m(6) * t65 + m(7) * t35; m(7) * (t35 * t29 + t58 * t80 + t59 * t79) + m(6) * (t65 * t31 + t81 * t91 + t82 * t90) + t221; m(7) * (t30 * t35 + t32 * t58 + t33 * t59) + m(6) * (t65 * t34 + t56 * t81 + t57 * t82) + t223; m(7) * (t35 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(6) * (t65 ^ 2 + t81 ^ 2 + t82 ^ 2) + t223; m(7) * (t157 * t64 + t159 * t63); m(7) * t265; m(7) * (t157 * t79 + t159 * t80 + t265 * t29); m(7) * (t157 * t33 + t159 * t32 + t265 * t30); m(7) * (t157 * t59 + t159 * t58 + t265 * t35); m(7) * (t208 ^ 2 * t214 ^ 2 + t157 ^ 2 + t159 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
