% Calculate joint inertia matrix for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR15_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR15_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR15_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:58
% EndTime: 2024-09-27 22:24:01
% DurationCPUTime: 1.55s
% Computational Cost: add. (29924->499), mult. (34363->694), div. (0->0), fcn. (39578->26), ass. (0->276)
t272 = sin(qJ(4));
t277 = cos(qJ(4));
t267 = sin(pkin(6));
t343 = pkin(11) * t267;
t231 = pkin(4) * t277 + t272 * t343 + pkin(3);
t234 = -pkin(4) * t272 + t277 * t343;
t273 = sin(qJ(3));
t278 = cos(qJ(3));
t192 = t231 * t278 + t234 * t273 + pkin(2);
t193 = -t231 * t273 + t234 * t278;
t281 = pkin(8) + pkin(9);
t265 = pkin(10) + t281;
t269 = cos(pkin(6));
t241 = pkin(11) * t269 + t265;
t268 = sin(pkin(5));
t270 = cos(pkin(5));
t279 = cos(qJ(2));
t274 = sin(qJ(2));
t333 = t270 * t274;
t115 = t193 * t270 * t279 - t192 * t333 + t241 * t268;
t128 = t192 * t279 + t193 * t274 + pkin(1);
t275 = sin(qJ(1));
t280 = cos(qJ(1));
t266 = qJ(2) + qJ(3);
t263 = qJ(4) + t266;
t255 = sin(t263);
t256 = cos(t263);
t271 = sin(qJ(5));
t276 = cos(qJ(5));
t327 = t275 * t276;
t306 = t270 * t327;
t332 = t271 * t275;
t308 = t270 * t332;
t328 = t275 * t268;
t310 = t267 * t328;
t324 = t276 * t280;
t331 = t271 * t280;
t150 = (-t269 * t308 + t324) * t256 + (-t269 * t331 - t306) * t255 + t271 * t310;
t151 = (-t269 * t306 - t331) * t256 + (-t269 * t324 + t308) * t255 + t276 * t310;
t337 = t256 * t275;
t339 = t255 * t280;
t188 = t269 * t328 + (t270 * t337 + t339) * t267;
t90 = t150 * rSges(6,1) + t151 * rSges(6,2) + t188 * rSges(6,3);
t47 = t115 * t275 + t128 * t280 + t90;
t257 = pkin(3) * t278 + pkin(2);
t344 = pkin(3) * t273;
t313 = t279 * t344;
t200 = -t265 * t268 + (t257 * t274 + t313) * t270;
t258 = t279 * pkin(2) + pkin(1);
t262 = cos(t266);
t235 = pkin(3) * t262 + t258;
t305 = t270 * t324;
t307 = t270 * t331;
t322 = t280 * t268;
t309 = t267 * t322;
t152 = (t269 * t307 + t327) * t256 + (-t269 * t332 + t305) * t255 - t271 * t309;
t153 = (t269 * t305 - t332) * t256 + (-t269 * t327 - t307) * t255 - t276 * t309;
t336 = t256 * t280;
t340 = t255 * t275;
t189 = -t269 * t322 + (-t270 * t336 + t340) * t267;
t91 = rSges(6,1) * t152 + rSges(6,2) * t153 + rSges(6,3) * t189;
t350 = (-t115 - t200) * t280 + (t128 - t235) * t275 + t91;
t294 = t200 * t275 - t280 * t235;
t349 = t294 + t47;
t260 = pkin(5) - t266;
t254 = -qJ(4) + t260;
t240 = cos(t254) / 0.2e1;
t259 = pkin(5) + t266;
t253 = qJ(4) + t259;
t243 = cos(t253);
t221 = t240 + t243 / 0.2e1;
t196 = t221 * t280 - t340;
t239 = sin(t253) / 0.2e1;
t242 = sin(t254);
t220 = t239 - t242 / 0.2e1;
t197 = t220 * t280 + t337;
t131 = Icges(5,5) * t197 + Icges(5,6) * t196 - Icges(5,3) * t322;
t348 = t131 * t322;
t347 = -t275 / 0.2e1;
t346 = t280 / 0.2e1;
t198 = -t221 * t275 - t339;
t199 = -t220 * t275 + t336;
t132 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t328;
t133 = Icges(5,4) * t197 + Icges(5,2) * t196 - Icges(5,6) * t322;
t134 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t328;
t135 = Icges(5,1) * t197 + Icges(5,4) * t196 - Icges(5,5) * t322;
t136 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t328;
t219 = t239 + t242 / 0.2e1;
t222 = t240 - t243 / 0.2e1;
t171 = Icges(5,5) * t222 + Icges(5,6) * t219 + Icges(5,3) * t270;
t172 = Icges(5,4) * t222 + Icges(5,2) * t219 + Icges(5,6) * t270;
t173 = Icges(5,1) * t222 + Icges(5,4) * t219 + Icges(5,5) * t270;
t72 = -t171 * t322 + t172 * t196 + t173 * t197;
t21 = (-t132 * t322 + t134 * t196 + t136 * t197) * t328 - (t133 * t196 + t135 * t197 - t348) * t322 + t72 * t270;
t84 = Icges(6,5) * t150 + Icges(6,6) * t151 + Icges(6,3) * t188;
t86 = Icges(6,4) * t150 + Icges(6,2) * t151 + Icges(6,6) * t188;
t88 = Icges(6,1) * t150 + Icges(6,4) * t151 + Icges(6,5) * t188;
t29 = t152 * t88 + t153 * t86 + t189 * t84;
t85 = Icges(6,5) * t152 + Icges(6,6) * t153 + Icges(6,3) * t189;
t87 = Icges(6,4) * t152 + Icges(6,2) * t153 + Icges(6,6) * t189;
t89 = Icges(6,1) * t152 + Icges(6,4) * t153 + Icges(6,5) * t189;
t30 = t152 * t89 + t153 * t87 + t189 * t85;
t335 = t267 * t270;
t338 = t256 * t269;
t190 = t276 * t335 + (-t255 * t271 + t276 * t338) * t268;
t191 = t271 * t335 + (t255 * t276 + t271 * t338) * t268;
t214 = -t256 * t267 * t268 + t269 * t270;
t109 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t214;
t110 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t214;
t111 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t214;
t41 = t109 * t189 + t110 * t153 + t111 * t152;
t9 = t270 * t41 + (t275 * t29 - t280 * t30) * t268;
t345 = -t21 - t9;
t342 = t349 * t270;
t334 = t268 * t274;
t330 = t274 * t275;
t329 = t274 * t280;
t326 = t275 * t279;
t325 = t275 * t280;
t323 = t279 * t280;
t112 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t214;
t321 = -(t241 - t265) * t270 - ((-t193 - t344) * t279 + (t192 - t257) * t274) * t268 - t112;
t230 = pkin(2) * t333 - t268 * t281;
t125 = (t200 - t230) * t280 + (t235 - t258) * t275;
t314 = -t230 * t275 + t280 * t258;
t126 = -t294 - t314;
t319 = t125 * t328 + t126 * t322;
t121 = t270 * t126;
t138 = t199 * rSges(5,1) + t198 * rSges(5,2) + rSges(5,3) * t328;
t127 = t270 * t138;
t318 = t121 + t127;
t289 = -rSges(5,1) * t197 - rSges(5,2) * t196;
t137 = -rSges(5,3) * t322 - t289;
t92 = t137 * t328 + t138 * t322;
t317 = -t125 - t137;
t250 = cos(t259);
t251 = cos(t260);
t233 = t251 + t250;
t261 = sin(t266);
t205 = t233 * t346 - t275 * t261;
t248 = sin(t259);
t249 = sin(t260);
t232 = t248 - t249;
t206 = t232 * t346 + t275 * t262;
t290 = -rSges(4,1) * t206 - rSges(4,2) * t205;
t154 = -rSges(4,3) * t322 - t290;
t207 = t233 * t347 - t280 * t261;
t208 = t232 * t347 + t280 * t262;
t155 = t208 * rSges(4,1) + t207 * rSges(4,2) + rSges(4,3) * t328;
t96 = t154 * t328 + t155 * t322;
t247 = pkin(8) * t322;
t182 = t230 * t280 + t247 + (-pkin(1) + t258) * t275;
t286 = -pkin(1) * t280 - pkin(8) * t328;
t183 = t286 + t314;
t316 = t182 * t328 + t183 * t322;
t174 = (t265 - t281) * t270 + (t313 + (-pkin(2) + t257) * t274) * t268;
t177 = rSges(5,1) * t222 + rSges(5,2) * t219 + rSges(5,3) * t270;
t315 = -t174 - t177;
t312 = t121 + t342;
t311 = -t125 - t350;
t46 = t214 * t109 + t190 * t110 + t191 * t111;
t304 = -t174 + t321;
t303 = t270 * t171 + t219 * t172 + t222 * t173;
t228 = t248 / 0.2e1 + t249 / 0.2e1;
t229 = t251 / 0.2e1 - t250 / 0.2e1;
t184 = Icges(4,5) * t229 + Icges(4,6) * t228 + Icges(4,3) * t270;
t185 = Icges(4,4) * t229 + Icges(4,2) * t228 + Icges(4,6) * t270;
t186 = Icges(4,1) * t229 + Icges(4,4) * t228 + Icges(4,5) * t270;
t302 = t270 * t184 + t228 * t185 + t229 * t186;
t210 = Icges(3,3) * t270 + (Icges(3,5) * t274 + Icges(3,6) * t279) * t268;
t211 = Icges(3,6) * t270 + (Icges(3,4) * t274 + Icges(3,2) * t279) * t268;
t212 = Icges(3,5) * t270 + (Icges(3,1) * t274 + Icges(3,4) * t279) * t268;
t301 = t268 * t279 * t211 + t270 * t210 + t212 * t334;
t225 = -t270 * t326 - t329;
t226 = -t270 * t330 + t323;
t176 = t226 * rSges(3,1) + t225 * rSges(3,2) + rSges(3,3) * t328;
t300 = t328 / 0.2e1;
t299 = -t322 / 0.2e1;
t36 = t190 * t86 + t191 * t88 + t214 * t84;
t37 = t190 * t87 + t191 * t89 + t214 * t85;
t45 = t46 * t270;
t14 = t45 + (t275 * t36 - t280 * t37) * t268;
t56 = t131 * t270 + t133 * t219 + t135 * t222;
t57 = t132 * t270 + t134 * t219 + t136 * t222;
t73 = t171 * t328 + t172 * t198 + t173 * t199;
t27 = t150 * t88 + t151 * t86 + t188 * t84;
t28 = t150 * t89 + t151 * t87 + t188 * t85;
t40 = t109 * t188 + t110 * t151 + t111 * t150;
t8 = t270 * t40 + (t27 * t275 - t28 * t280) * t268;
t82 = t303 * t270;
t298 = (t8 - (t133 * t198 + t135 * t199) * t322 + (t132 * t328 + t134 * t198 + t136 * t199 - t348) * t328) * t328 + (t14 + t73 * t328 + t82 + (t275 * t57 - t280 * t56) * t268) * t270;
t19 = t349 * t322 + t350 * t328;
t297 = t268 * t321;
t296 = t268 * t315;
t187 = rSges(4,1) * t229 + rSges(4,2) * t228 + rSges(4,3) * t270;
t215 = pkin(2) * t334 + (-pkin(8) + t281) * t270;
t295 = t268 * (-t187 - t215);
t53 = t92 + t319;
t293 = t268 * t304;
t292 = t268 * (-t215 + t315);
t42 = t46 * t214;
t11 = t36 * t188 + t37 * t189 + t42;
t3 = t188 * t27 + t189 * t28 + t214 * t40;
t4 = t188 * t29 + t189 * t30 + t214 * t41;
t291 = t188 * t8 / 0.2e1 + t4 * t299 + t3 * t300 + t270 * t11 / 0.2e1 + t214 * t14 / 0.2e1 + t189 * t9 / 0.2e1;
t144 = Icges(4,5) * t206 + Icges(4,6) * t205 - Icges(4,3) * t322;
t145 = Icges(4,5) * t208 + Icges(4,6) * t207 + Icges(4,3) * t328;
t146 = Icges(4,4) * t206 + Icges(4,2) * t205 - Icges(4,6) * t322;
t147 = Icges(4,4) * t208 + Icges(4,2) * t207 + Icges(4,6) * t328;
t148 = Icges(4,1) * t206 + Icges(4,4) * t205 - Icges(4,5) * t322;
t149 = Icges(4,1) * t208 + Icges(4,4) * t207 + Icges(4,5) * t328;
t65 = t144 * t270 + t146 * t228 + t148 * t229;
t66 = t145 * t270 + t147 * t228 + t149 * t229;
t81 = t184 * t328 + t185 * t207 + t186 * t208;
t93 = t302 * t270;
t288 = ((t145 * t328 + t147 * t207 + t149 * t208) * t328 - (t144 * t328 + t146 * t207 + t148 * t208) * t322 + t81 * t270) * t328 + t270 * (t93 + (t275 * t66 - t280 * t65) * t268) + t298;
t16 = t19 + t319;
t287 = t268 * (-t215 + t304);
t285 = t345 * t322 + t298;
t223 = t270 * t323 - t330;
t224 = t270 * t329 + t326;
t175 = rSges(3,1) * t224 + rSges(3,2) * t223 - rSges(3,3) * t322;
t284 = t45 + t82 + (t36 + t40 + t57 + t73) * t300 + (t37 + t41 + t56 + t72) * t299;
t80 = -t184 * t322 + t185 * t205 + t186 * t206;
t26 = (-t145 * t322 + t147 * t205 + t149 * t206) * t328 - (-t144 * t322 + t146 * t205 + t148 * t206) * t322 + t80 * t270;
t283 = (-t26 + t345) * t322 + t288;
t282 = t284 + t93 + (t66 + t81) * t300 + (t65 + t80) * t299;
t237 = rSges(2,1) * t280 - rSges(2,2) * t275;
t236 = -rSges(2,1) * t275 - rSges(2,2) * t280;
t213 = rSges(3,3) * t270 + (rSges(3,1) * t274 + rSges(3,2) * t279) * t268;
t178 = t270 * t183;
t170 = Icges(3,1) * t226 + Icges(3,4) * t225 + Icges(3,5) * t328;
t169 = Icges(3,1) * t224 + Icges(3,4) * t223 - Icges(3,5) * t322;
t168 = Icges(3,4) * t226 + Icges(3,2) * t225 + Icges(3,6) * t328;
t167 = Icges(3,4) * t224 + Icges(3,2) * t223 - Icges(3,6) * t322;
t166 = Icges(3,5) * t226 + Icges(3,6) * t225 + Icges(3,3) * t328;
t165 = Icges(3,5) * t224 + Icges(3,6) * t223 - Icges(3,3) * t322;
t159 = -t286 + t176;
t158 = -pkin(1) * t275 - t175 + t247;
t141 = t270 * t155;
t130 = -t175 * t270 - t213 * t322;
t129 = t176 * t270 - t213 * t328;
t118 = t301 * t270;
t117 = t155 + t314;
t116 = -t258 * t275 + (rSges(4,3) * t268 - t230) * t280 + t290;
t113 = (t175 * t275 + t176 * t280) * t268;
t108 = t210 * t328 + t211 * t225 + t212 * t226;
t107 = -t210 * t322 + t211 * t223 + t212 * t224;
t106 = -t154 * t270 - t187 * t322;
t105 = -t187 * t328 + t141;
t103 = -t294 + t138;
t102 = -t235 * t275 + (rSges(5,3) * t268 - t200) * t280 + t289;
t100 = -t137 * t270 - t177 * t322;
t99 = -t177 * t328 + t127;
t95 = t166 * t270 + (t168 * t279 + t170 * t274) * t268;
t94 = t165 * t270 + (t167 * t279 + t169 * t274) * t268;
t77 = (-t154 - t182) * t270 + t280 * t295;
t76 = t275 * t295 + t141 + t178;
t62 = t316 + t96;
t59 = t317 * t270 + t280 * t296;
t58 = t275 * t296 + t318;
t52 = (-t182 + t317) * t270 + t280 * t292;
t51 = t275 * t292 + t178 + t318;
t50 = t112 * t189 - t214 * t91;
t49 = -t112 * t188 + t214 * t90;
t48 = t115 * t280 - t128 * t275 - t91;
t44 = t53 + t316;
t43 = t188 * t91 - t189 * t90;
t32 = -t270 * t350 + t280 * t297;
t31 = t275 * t297 + t342;
t23 = t311 * t270 + t280 * t293;
t22 = t275 * t293 + t312;
t18 = (-t182 + t311) * t270 + t280 * t287;
t17 = t275 * t287 + t178 + t312;
t15 = t16 + t316;
t1 = [Icges(2,3) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2) + m(4) * (t116 ^ 2 + t117 ^ 2) + m(3) * (t158 ^ 2 + t159 ^ 2) + m(2) * (t236 ^ 2 + t237 ^ 2) + t301 + t302 + t303 + t46; t282 + t118 + ((-t107 / 0.2e1 - t94 / 0.2e1) * t280 + (t108 / 0.2e1 + t95 / 0.2e1) * t275) * t268 + m(6) * (t17 * t47 + t18 * t48) + m(5) * (t102 * t52 + t103 * t51) + m(4) * (t116 * t77 + t117 * t76) + m(3) * (t129 * t159 + t130 * t158); t270 * t118 + (-t280 * t9 - t280 * t21 - t280 * t26 + (t275 * ((t168 * t225 + t170 * t226) * t275 - (t167 * t225 + t169 * t226) * t280) - t280 * ((t168 * t223 + t170 * t224) * t275 - (t167 * t223 + t169 * t224) * t280) + (t275 * (t166 * t275 ^ 2 - t165 * t325) - t280 * (t165 * t280 ^ 2 - t166 * t325)) * t268) * t268 + ((-t107 - t94) * t280 + (t108 + t95) * t275) * t270) * t268 + m(6) * (t15 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(5) * (t44 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(4) * (t62 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (t113 ^ 2 + t129 ^ 2 + t130 ^ 2) + t288; t282 + m(6) * (t22 * t47 + t23 * t48) + m(5) * (t102 * t59 + t103 * t58) + m(4) * (t105 * t117 + t106 * t116); m(6) * (t15 * t16 + t17 * t22 + t18 * t23) + m(5) * (t44 * t53 + t51 * t58 + t52 * t59) + m(4) * (t105 * t76 + t106 * t77 + t62 * t96) + t283; m(6) * (t16 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(5) * (t53 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t105 ^ 2 + t106 ^ 2 + t96 ^ 2) + t283; m(6) * (t31 * t47 + t32 * t48) + m(5) * (t100 * t102 + t103 * t99) + t284; m(6) * (t15 * t19 + t17 * t31 + t18 * t32) + m(5) * (t100 * t52 + t44 * t92 + t51 * t99) + t285; m(6) * (t16 * t19 + t22 * t31 + t23 * t32) + m(5) * (t100 * t59 + t53 * t92 + t58 * t99) + t285; m(5) * (t100 ^ 2 + t92 ^ 2 + t99 ^ 2) + m(6) * (t19 ^ 2 + t31 ^ 2 + t32 ^ 2) + t285; t42 + m(6) * (t47 * t49 + t48 * t50) + (t41 / 0.2e1 + t37 / 0.2e1) * t189 + (t40 / 0.2e1 + t36 / 0.2e1) * t188; m(6) * (t15 * t43 + t17 * t49 + t18 * t50) + t291; m(6) * (t16 * t43 + t22 * t49 + t23 * t50) + t291; m(6) * (t19 * t43 + t31 * t49 + t32 * t50) + t291; m(6) * (t43 ^ 2 + t49 ^ 2 + t50 ^ 2) + t188 * t3 + t189 * t4 + t214 * t11;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
