% Calculate joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:25
% EndTime: 2019-03-08 23:05:45
% DurationCPUTime: 8.49s
% Computational Cost: add. (31303->547), mult. (48411->782), div. (0->0), fcn. (61283->14), ass. (0->262)
t367 = Icges(5,2) + Icges(6,3);
t271 = sin(pkin(11));
t274 = cos(pkin(11));
t280 = cos(qJ(2));
t275 = cos(pkin(6));
t278 = sin(qJ(2));
t326 = t275 * t278;
t255 = t271 * t280 + t274 * t326;
t324 = qJ(3) + qJ(4);
t268 = sin(t324);
t295 = cos(t324);
t272 = sin(pkin(6));
t332 = t272 * t274;
t241 = t255 * t295 - t268 * t332;
t325 = t275 * t280;
t254 = t271 * t278 - t274 * t325;
t270 = sin(pkin(12));
t273 = cos(pkin(12));
t209 = -t241 * t270 + t254 * t273;
t335 = t254 * t270;
t210 = t241 * t273 + t335;
t288 = t272 * t295;
t240 = t255 * t268 + t274 * t288;
t366 = -Icges(5,4) * t241 + Icges(6,5) * t210 - Icges(5,6) * t254 + Icges(6,6) * t209 + t240 * t367;
t257 = -t271 * t326 + t274 * t280;
t333 = t271 * t272;
t243 = t257 * t295 + t268 * t333;
t256 = t271 * t325 + t274 * t278;
t211 = -t243 * t270 + t256 * t273;
t334 = t256 * t270;
t212 = t243 * t273 + t334;
t242 = t257 * t268 - t271 * t288;
t365 = -Icges(5,4) * t243 + Icges(6,5) * t212 - Icges(5,6) * t256 + Icges(6,6) * t211 + t242 * t367;
t253 = t275 * t268 + t278 * t288;
t328 = t272 * t280;
t238 = -t253 * t270 - t273 * t328;
t308 = t270 * t328;
t239 = t253 * t273 - t308;
t330 = t272 * t278;
t252 = t268 * t330 - t275 * t295;
t364 = -Icges(5,4) * t253 + Icges(6,5) * t239 + Icges(5,6) * t328 + Icges(6,6) * t238 + t252 * t367;
t146 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t240;
t148 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t240;
t174 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t254;
t178 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t254;
t363 = t146 * t209 + t148 * t210 + t174 * t254 + t178 * t241 + t240 * t366;
t147 = Icges(6,4) * t212 + Icges(6,2) * t211 + Icges(6,6) * t242;
t149 = Icges(6,1) * t212 + Icges(6,4) * t211 + Icges(6,5) * t242;
t175 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t256;
t179 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t256;
t362 = t147 * t209 + t149 * t210 + t175 * t254 + t179 * t241 + t240 * t365;
t361 = t146 * t211 + t148 * t212 + t174 * t256 + t178 * t243 + t242 * t366;
t360 = t147 * t211 + t149 * t212 + t175 * t256 + t179 * t243 + t242 * t365;
t359 = t146 * t238 + t148 * t239 - t174 * t328 + t178 * t253 + t252 * t366;
t358 = t147 * t238 + t149 * t239 - t175 * t328 + t179 * t253 + t252 * t365;
t170 = Icges(6,4) * t239 + Icges(6,2) * t238 + Icges(6,6) * t252;
t171 = Icges(6,1) * t239 + Icges(6,4) * t238 + Icges(6,5) * t252;
t213 = Icges(5,5) * t253 - Icges(5,6) * t252 - Icges(5,3) * t328;
t215 = Icges(5,1) * t253 - Icges(5,4) * t252 - Icges(5,5) * t328;
t357 = t170 * t209 + t171 * t210 + t213 * t254 + t215 * t241 + t240 * t364;
t356 = t170 * t211 + t171 * t212 + t213 * t256 + t215 * t243 + t242 * t364;
t355 = t170 * t238 + t171 * t239 - t213 * t328 + t215 * t253 + t252 * t364;
t347 = m(6) + m(7);
t269 = pkin(12) + qJ(6);
t266 = sin(t269);
t267 = cos(t269);
t203 = -t241 * t266 + t254 * t267;
t204 = t241 * t267 + t254 * t266;
t141 = rSges(7,1) * t204 + rSges(7,2) * t203 + rSges(7,3) * t240;
t337 = pkin(5) * t273;
t321 = pkin(5) * t335 + pkin(10) * t240 + t241 * t337 + t141;
t234 = -t253 * t266 - t267 * t328;
t235 = t253 * t267 - t266 * t328;
t164 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t252;
t354 = -pkin(5) * t308 + pkin(10) * t252 + t253 * t337 + t164;
t135 = Icges(7,5) * t204 + Icges(7,6) * t203 + Icges(7,3) * t240;
t137 = Icges(7,4) * t204 + Icges(7,2) * t203 + Icges(7,6) * t240;
t139 = Icges(7,1) * t204 + Icges(7,4) * t203 + Icges(7,5) * t240;
t69 = t135 * t240 + t137 * t203 + t139 * t204;
t205 = -t243 * t266 + t256 * t267;
t206 = t243 * t267 + t256 * t266;
t136 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t242;
t138 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t242;
t140 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t242;
t70 = t136 * t240 + t138 * t203 + t140 * t204;
t161 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t252;
t162 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t252;
t163 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t252;
t87 = t161 * t240 + t162 * t203 + t163 * t204;
t11 = t254 * t69 + t256 * t70 - t328 * t87;
t353 = t254 * t363 + t362 * t256 - t357 * t328 + t11;
t71 = t135 * t242 + t137 * t205 + t139 * t206;
t72 = t136 * t242 + t138 * t205 + t140 * t206;
t88 = t161 * t242 + t162 * t205 + t163 * t206;
t12 = t254 * t71 + t256 * t72 - t328 * t88;
t352 = t254 * t361 + t256 * t360 - t328 * t356 + t12;
t15 = t87 * t275 + (t271 * t70 - t274 * t69) * t272;
t351 = t15 + t357 * t275 + (t362 * t271 - t274 * t363) * t272;
t16 = t88 * t275 + (t271 * t72 - t274 * t71) * t272;
t350 = t16 + t356 * t275 + (t271 * t360 - t274 * t361) * t272;
t79 = t135 * t252 + t137 * t234 + t139 * t235;
t80 = t136 * t252 + t138 * t234 + t140 * t235;
t94 = t161 * t252 + t162 * t234 + t163 * t235;
t31 = t254 * t79 + t256 * t80 - t328 * t94;
t349 = t254 * t359 + t256 * t358 - t328 * t355 + t31;
t33 = t94 * t275 + (t271 * t80 - t274 * t79) * t272;
t348 = t33 + t355 * t275 + (t271 * t358 - t274 * t359) * t272;
t346 = t240 / 0.2e1;
t345 = t242 / 0.2e1;
t344 = t252 / 0.2e1;
t343 = t254 / 0.2e1;
t342 = t256 / 0.2e1;
t341 = t271 / 0.2e1;
t340 = -t274 / 0.2e1;
t339 = t275 / 0.2e1;
t279 = cos(qJ(3));
t338 = pkin(3) * t279;
t277 = sin(qJ(3));
t331 = t272 * t277;
t329 = t272 * t279;
t327 = t275 * t277;
t150 = rSges(6,1) * t210 + rSges(6,2) * t209 + rSges(6,3) * t240;
t200 = t241 * pkin(4) + t240 * qJ(5);
t185 = t256 * t200;
t322 = t256 * t150 + t185;
t142 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t242;
t320 = pkin(5) * t334 + pkin(10) * t242 + t243 * t337 + t142;
t151 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t242;
t201 = t243 * pkin(4) + t242 * qJ(5);
t319 = -t151 - t201;
t306 = t274 * t331;
t183 = -pkin(3) * t306 + pkin(9) * t254 + t255 * t338;
t230 = pkin(3) * t327 + (-pkin(9) * t280 + t278 * t338) * t272;
t318 = t183 * t328 + t254 * t230;
t172 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t252;
t228 = t253 * pkin(4) + t252 * qJ(5);
t317 = -t172 - t228;
t307 = t271 * t331;
t184 = pkin(3) * t307 + pkin(9) * t256 + t257 * t338;
t237 = t257 * pkin(2) + t256 * pkin(8);
t233 = t275 * t237;
t316 = t275 * t184 + t233;
t182 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t256;
t315 = -t182 - t184;
t236 = t255 * pkin(2) + t254 * pkin(8);
t314 = -t183 - t236;
t313 = t200 * t328 + t254 * t228;
t181 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t254;
t216 = rSges(5,1) * t253 - rSges(5,2) * t252 - rSges(5,3) * t328;
t133 = t181 * t328 + t254 * t216;
t312 = -t216 - t230;
t311 = t236 * t333 + t237 * t332;
t305 = t256 * t321 + t185;
t304 = -t201 - t320;
t303 = -t184 + t319;
t302 = -t228 - t354;
t301 = -t230 + t317;
t300 = t275 * t201 + t316;
t299 = -t200 + t314;
t296 = -t328 / 0.2e1;
t258 = t275 * t279 - t277 * t330;
t259 = t278 * t329 + t327;
t229 = rSges(4,1) * t259 + rSges(4,2) * t258 - rSges(4,3) * t328;
t260 = (pkin(2) * t278 - pkin(8) * t280) * t272;
t294 = (-t229 - t260) * t272;
t293 = -t184 + t304;
t292 = -t230 + t302;
t291 = t183 * t333 + t184 * t332 + t311;
t92 = t150 * t328 + t254 * t172 + t313;
t290 = (-t260 + t312) * t272;
t28 = t240 * t79 + t242 * t80 + t252 * t94;
t3 = t240 * t69 + t242 * t70 + t252 * t87;
t4 = t240 * t71 + t242 * t72 + t252 * t88;
t289 = t11 * t346 + t12 * t345 + t28 * t296 + t3 * t343 + t31 * t344 + t4 * t342;
t287 = t254 * t353 + t256 * t352;
t286 = (-t260 + t301) * t272;
t285 = t200 * t333 + t201 * t332 + t291;
t67 = t254 * t354 + t321 * t328 + t313;
t284 = (-t260 + t292) * t272;
t283 = -t328 * t349 + t287;
t282 = t351 * t343 + t350 * t342 + t349 * t339 + t352 * t333 / 0.2e1 - t353 * t332 / 0.2e1 + t348 * t296;
t251 = t275 * rSges(3,3) + (rSges(3,1) * t278 + rSges(3,2) * t280) * t272;
t250 = Icges(3,5) * t275 + (Icges(3,1) * t278 + Icges(3,4) * t280) * t272;
t249 = Icges(3,6) * t275 + (Icges(3,4) * t278 + Icges(3,2) * t280) * t272;
t248 = Icges(3,3) * t275 + (Icges(3,5) * t278 + Icges(3,6) * t280) * t272;
t247 = t257 * t279 + t307;
t246 = -t257 * t277 + t271 * t329;
t245 = t255 * t279 - t306;
t244 = -t255 * t277 - t274 * t329;
t227 = Icges(4,1) * t259 + Icges(4,4) * t258 - Icges(4,5) * t328;
t226 = Icges(4,4) * t259 + Icges(4,2) * t258 - Icges(4,6) * t328;
t225 = Icges(4,5) * t259 + Icges(4,6) * t258 - Icges(4,3) * t328;
t224 = rSges(3,1) * t257 - rSges(3,2) * t256 + rSges(3,3) * t333;
t223 = rSges(3,1) * t255 - rSges(3,2) * t254 - rSges(3,3) * t332;
t222 = Icges(3,1) * t257 - Icges(3,4) * t256 + Icges(3,5) * t333;
t221 = Icges(3,1) * t255 - Icges(3,4) * t254 - Icges(3,5) * t332;
t220 = Icges(3,4) * t257 - Icges(3,2) * t256 + Icges(3,6) * t333;
t219 = Icges(3,4) * t255 - Icges(3,2) * t254 - Icges(3,6) * t332;
t218 = Icges(3,5) * t257 - Icges(3,6) * t256 + Icges(3,3) * t333;
t217 = Icges(3,5) * t255 - Icges(3,6) * t254 - Icges(3,3) * t332;
t198 = -t223 * t275 - t251 * t332;
t197 = t224 * t275 - t251 * t333;
t193 = rSges(4,1) * t247 + rSges(4,2) * t246 + rSges(4,3) * t256;
t192 = rSges(4,1) * t245 + rSges(4,2) * t244 + rSges(4,3) * t254;
t191 = Icges(4,1) * t247 + Icges(4,4) * t246 + Icges(4,5) * t256;
t190 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t254;
t189 = Icges(4,4) * t247 + Icges(4,2) * t246 + Icges(4,6) * t256;
t188 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t254;
t187 = Icges(4,5) * t247 + Icges(4,6) * t246 + Icges(4,3) * t256;
t186 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t254;
t160 = (t223 * t271 + t224 * t274) * t272;
t158 = t256 * t183;
t157 = t256 * t181;
t153 = -t193 * t328 - t229 * t256;
t152 = t192 * t328 + t229 * t254;
t134 = -t182 * t328 - t256 * t216;
t128 = -t225 * t328 + t226 * t258 + t227 * t259;
t124 = t192 * t256 - t193 * t254;
t123 = (-t192 - t236) * t275 + t274 * t294;
t122 = t275 * t193 + t271 * t294 + t233;
t120 = t225 * t256 + t226 * t246 + t227 * t247;
t119 = t225 * t254 + t226 * t244 + t227 * t245;
t118 = -t182 * t254 + t157;
t115 = (t192 * t271 + t193 * t274) * t272 + t311;
t114 = -t187 * t328 + t189 * t258 + t191 * t259;
t113 = -t186 * t328 + t188 * t258 + t190 * t259;
t110 = t142 * t252 - t164 * t242;
t109 = -t141 * t252 + t164 * t240;
t108 = t256 * t312 + t315 * t328;
t107 = t133 + t318;
t106 = t187 * t256 + t189 * t246 + t191 * t247;
t105 = t186 * t256 + t188 * t246 + t190 * t247;
t104 = t187 * t254 + t189 * t244 + t191 * t245;
t103 = t186 * t254 + t188 * t244 + t190 * t245;
t98 = (-t181 + t314) * t275 + t274 * t290;
t97 = t275 * t182 + t271 * t290 + t316;
t95 = t141 * t242 - t142 * t240;
t93 = t256 * t317 + t319 * t328;
t91 = t254 * t315 + t157 + t158;
t86 = (t181 * t271 + t182 * t274) * t272 + t291;
t85 = t254 * t319 + t322;
t82 = t256 * t301 + t303 * t328;
t81 = t92 + t318;
t78 = (-t150 + t299) * t275 + t274 * t286;
t77 = t275 * t151 + t271 * t286 + t300;
t68 = t256 * t302 + t304 * t328;
t66 = t254 * t303 + t158 + t322;
t65 = (t150 * t271 + t151 * t274) * t272 + t285;
t64 = t128 * t275 + (-t113 * t274 + t114 * t271) * t272;
t63 = t113 * t254 + t114 * t256 - t128 * t328;
t60 = t120 * t275 + (-t105 * t274 + t106 * t271) * t272;
t59 = t119 * t275 + (-t103 * t274 + t104 * t271) * t272;
t58 = t256 * t292 + t293 * t328;
t57 = t67 + t318;
t54 = (t299 - t321) * t275 + t274 * t284;
t53 = t271 * t284 + t275 * t320 + t300;
t52 = t105 * t254 + t106 * t256 - t120 * t328;
t51 = t103 * t254 + t104 * t256 - t119 * t328;
t50 = t254 * t304 + t305;
t39 = t254 * t293 + t158 + t305;
t38 = (t271 * t321 + t274 * t320) * t272 + t285;
t1 = [m(2) + m(3) + m(4) + m(5) + t347; m(3) * t160 + m(4) * t115 + m(5) * t86 + m(6) * t65 + m(7) * t38; m(7) * (t38 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t65 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t86 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (t115 ^ 2 + t122 ^ 2 + t123 ^ 2) + m(3) * (t160 ^ 2 + t197 ^ 2 + t198 ^ 2) + (t60 + (t218 * t333 - t220 * t256 + t222 * t257) * t333 + t350) * t333 + (-t59 + (-t217 * t332 - t219 * t254 + t221 * t255) * t332 + (-t217 * t333 + t218 * t332 + t219 * t256 + t220 * t254 - t221 * t257 - t222 * t255) * t333 - t351) * t332 + (-(-t248 * t332 - t254 * t249 + t255 * t250) * t332 + (t248 * t333 - t256 * t249 + t257 * t250) * t333 + t64 + ((t220 * t280 + t222 * t278) * t271 - (t219 * t280 + t221 * t278) * t274) * t272 ^ 2 + ((-t217 * t274 + t218 * t271 + t249 * t280 + t250 * t278) * t272 + t275 * t248) * t275 + t348) * t275; m(4) * t124 + m(5) * t91 + m(6) * t66 + m(7) * t39; m(7) * (t38 * t39 + t53 * t58 + t54 * t57) + m(6) * (t65 * t66 + t77 * t82 + t78 * t81) + m(5) * (t107 * t98 + t108 * t97 + t86 * t91) + m(4) * (t115 * t124 + t122 * t153 + t123 * t152) + (t51 * t340 + t52 * t341 - t280 * t64 / 0.2e1) * t272 + t282 + t59 * t343 + t60 * t342 + t63 * t339; t254 * t51 + t256 * t52 + (-t63 - t349) * t328 + m(7) * (t39 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t66 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t107 ^ 2 + t108 ^ 2 + t91 ^ 2) + m(4) * (t124 ^ 2 + t152 ^ 2 + t153 ^ 2) + t287; m(5) * t118 + m(6) * t85 + m(7) * t50; m(7) * (t38 * t50 + t53 * t68 + t54 * t67) + m(6) * (t65 * t85 + t77 * t93 + t78 * t92) + m(5) * (t118 * t86 + t133 * t98 + t134 * t97) + t282; m(7) * (t39 * t50 + t57 * t67 + t58 * t68) + m(6) * (t66 * t85 + t81 * t92 + t82 * t93) + m(5) * (t107 * t133 + t108 * t134 + t118 * t91) + t283; m(7) * (t50 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t85 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(5) * (t118 ^ 2 + t133 ^ 2 + t134 ^ 2) + t283; t252 * t347; m(7) * (t240 * t53 + t242 * t54 + t252 * t38) + m(6) * (t240 * t77 + t242 * t78 + t252 * t65); m(7) * (t240 * t58 + t242 * t57 + t252 * t39) + m(6) * (t240 * t82 + t242 * t81 + t252 * t66); m(7) * (t240 * t68 + t242 * t67 + t252 * t50) + m(6) * (t240 * t93 + t242 * t92 + t252 * t85); (t240 ^ 2 + t242 ^ 2 + t252 ^ 2) * t347; m(7) * t95; t16 * t345 + t28 * t339 + m(7) * (t109 * t54 + t110 * t53 + t38 * t95) + t33 * t344 + t15 * t346 + (t3 * t340 + t341 * t4) * t272; m(7) * (t109 * t57 + t110 * t58 + t39 * t95) + t289; m(7) * (t109 * t67 + t110 * t68 + t50 * t95) + t289; m(7) * (t109 * t242 + t110 * t240 + t252 * t95); t242 * t4 + t240 * t3 + t252 * t28 + m(7) * (t109 ^ 2 + t110 ^ 2 + t95 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
