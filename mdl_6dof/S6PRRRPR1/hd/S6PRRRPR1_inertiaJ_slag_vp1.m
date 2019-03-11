% Calculate joint inertia matrix for
% S6PRRRPR1
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:15
% EndTime: 2019-03-08 23:00:34
% DurationCPUTime: 8.96s
% Computational Cost: add. (32442->548), mult. (44365->788), div. (0->0), fcn. (55124->14), ass. (0->261)
t368 = Icges(5,3) + Icges(6,3);
t273 = sin(pkin(11));
t275 = cos(pkin(11));
t282 = cos(qJ(2));
t276 = cos(pkin(6));
t279 = sin(qJ(2));
t328 = t276 * t279;
t255 = t273 * t282 + t275 * t328;
t272 = qJ(3) + qJ(4);
t299 = pkin(12) + t272;
t266 = sin(t299);
t274 = sin(pkin(6));
t293 = cos(t299);
t288 = t274 * t293;
t230 = t255 * t266 + t275 * t288;
t334 = t274 * t275;
t231 = t255 * t293 - t266 * t334;
t268 = sin(t272);
t269 = cos(t272);
t238 = -t255 * t268 - t269 * t334;
t239 = t255 * t269 - t268 * t334;
t327 = t276 * t282;
t254 = t273 * t279 - t275 * t327;
t367 = Icges(5,5) * t239 + Icges(6,5) * t231 + Icges(5,6) * t238 - Icges(6,6) * t230 + t254 * t368;
t257 = -t273 * t328 + t275 * t282;
t232 = t257 * t266 - t273 * t288;
t335 = t273 * t274;
t233 = t257 * t293 + t266 * t335;
t240 = -t257 * t268 + t269 * t335;
t241 = t257 * t269 + t268 * t335;
t256 = t273 * t327 + t275 * t279;
t366 = Icges(5,5) * t241 + Icges(6,5) * t233 + Icges(5,6) * t240 - Icges(6,6) * t232 + t256 * t368;
t332 = t274 * t279;
t246 = t266 * t332 - t276 * t293;
t247 = t276 * t266 + t279 * t288;
t252 = -t268 * t332 + t269 * t276;
t253 = t268 * t276 + t269 * t332;
t330 = t274 * t282;
t365 = Icges(5,5) * t253 + Icges(6,5) * t247 + Icges(5,6) * t252 - Icges(6,6) * t246 - t330 * t368;
t159 = Icges(6,4) * t231 - Icges(6,2) * t230 + Icges(6,6) * t254;
t161 = Icges(6,1) * t231 - Icges(6,4) * t230 + Icges(6,5) * t254;
t173 = Icges(5,4) * t239 + Icges(5,2) * t238 + Icges(5,6) * t254;
t175 = Icges(5,1) * t239 + Icges(5,4) * t238 + Icges(5,5) * t254;
t364 = -t159 * t230 + t161 * t231 + t173 * t238 + t175 * t239 + t254 * t367;
t160 = Icges(6,4) * t233 - Icges(6,2) * t232 + Icges(6,6) * t256;
t162 = Icges(6,1) * t233 - Icges(6,4) * t232 + Icges(6,5) * t256;
t174 = Icges(5,4) * t241 + Icges(5,2) * t240 + Icges(5,6) * t256;
t176 = Icges(5,1) * t241 + Icges(5,4) * t240 + Icges(5,5) * t256;
t363 = -t160 * t230 + t162 * t231 + t174 * t238 + t176 * t239 + t254 * t366;
t362 = -t159 * t232 + t161 * t233 + t173 * t240 + t175 * t241 + t256 * t367;
t361 = -t160 * t232 + t162 * t233 + t174 * t240 + t176 * t241 + t256 * t366;
t360 = -t159 * t246 + t161 * t247 + t173 * t252 + t175 * t253 - t330 * t367;
t359 = -t160 * t246 + t162 * t247 + t174 * t252 + t176 * t253 - t330 * t366;
t205 = Icges(6,4) * t247 - Icges(6,2) * t246 - Icges(6,6) * t330;
t206 = Icges(6,1) * t247 - Icges(6,4) * t246 - Icges(6,5) * t330;
t209 = Icges(5,4) * t253 + Icges(5,2) * t252 - Icges(5,6) * t330;
t210 = Icges(5,1) * t253 + Icges(5,4) * t252 - Icges(5,5) * t330;
t358 = -t205 * t230 + t206 * t231 + t209 * t238 + t210 * t239 + t254 * t365;
t357 = -t205 * t232 + t206 * t233 + t209 * t240 + t210 * t241 + t256 * t365;
t356 = -t205 * t246 + t206 * t247 + t209 * t252 + t210 * t253 - t330 * t365;
t277 = sin(qJ(6));
t280 = cos(qJ(6));
t199 = -t231 * t277 + t254 * t280;
t200 = t231 * t280 + t254 * t277;
t136 = rSges(7,1) * t200 + rSges(7,2) * t199 + rSges(7,3) * t230;
t326 = pkin(5) * t231 + pkin(10) * t230 + t136;
t236 = -t247 * t277 - t280 * t330;
t237 = t247 * t280 - t277 * t330;
t156 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t246;
t355 = pkin(5) * t247 + pkin(10) * t246 + t156;
t130 = Icges(7,5) * t200 + Icges(7,6) * t199 + Icges(7,3) * t230;
t132 = Icges(7,4) * t200 + Icges(7,2) * t199 + Icges(7,6) * t230;
t134 = Icges(7,1) * t200 + Icges(7,4) * t199 + Icges(7,5) * t230;
t69 = t130 * t230 + t132 * t199 + t134 * t200;
t201 = -t233 * t277 + t256 * t280;
t202 = t233 * t280 + t256 * t277;
t131 = Icges(7,5) * t202 + Icges(7,6) * t201 + Icges(7,3) * t232;
t133 = Icges(7,4) * t202 + Icges(7,2) * t201 + Icges(7,6) * t232;
t135 = Icges(7,1) * t202 + Icges(7,4) * t201 + Icges(7,5) * t232;
t70 = t131 * t230 + t133 * t199 + t135 * t200;
t153 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t246;
t154 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t246;
t155 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t246;
t81 = t153 * t230 + t154 * t199 + t155 * t200;
t11 = t254 * t69 + t256 * t70 - t330 * t81;
t354 = t364 * t254 + t363 * t256 - t358 * t330 + t11;
t71 = t130 * t232 + t132 * t201 + t134 * t202;
t72 = t131 * t232 + t133 * t201 + t135 * t202;
t82 = t153 * t232 + t154 * t201 + t155 * t202;
t12 = t254 * t71 + t256 * t72 - t330 * t82;
t353 = t362 * t254 + t361 * t256 - t357 * t330 + t12;
t15 = t276 * t81 + (t273 * t70 - t275 * t69) * t274;
t352 = t15 + t358 * t276 + (t363 * t273 - t364 * t275) * t274;
t16 = t276 * t82 + (t273 * t72 - t275 * t71) * t274;
t351 = t16 + t357 * t276 + (t361 * t273 - t362 * t275) * t274;
t73 = t130 * t246 + t132 * t236 + t134 * t237;
t74 = t131 * t246 + t133 * t236 + t135 * t237;
t85 = t153 * t246 + t154 * t236 + t155 * t237;
t21 = t254 * t73 + t256 * t74 - t330 * t85;
t350 = t360 * t254 + t359 * t256 - t356 * t330 + t21;
t23 = t276 * t85 + (t273 * t74 - t275 * t73) * t274;
t349 = t23 + t356 * t276 + (t359 * t273 - t360 * t275) * t274;
t348 = t274 ^ 2;
t347 = -m(6) - m(7);
t346 = t230 / 0.2e1;
t345 = t232 / 0.2e1;
t344 = t246 / 0.2e1;
t343 = t254 / 0.2e1;
t342 = t256 / 0.2e1;
t341 = t273 / 0.2e1;
t340 = -t275 / 0.2e1;
t339 = t276 / 0.2e1;
t281 = cos(qJ(3));
t337 = t281 * pkin(3);
t278 = sin(qJ(3));
t333 = t274 * t278;
t331 = t274 * t281;
t329 = t276 * t278;
t137 = rSges(7,1) * t202 + rSges(7,2) * t201 + rSges(7,3) * t232;
t325 = pkin(5) * t233 + pkin(10) * t232 + t137;
t298 = pkin(4) * t268;
t314 = pkin(4) * t269;
t145 = qJ(5) * t254 + t255 * t314 - t298 * t334;
t140 = t256 * t145;
t163 = rSges(6,1) * t231 - rSges(6,2) * t230 + rSges(6,3) * t254;
t324 = t256 * t163 + t140;
t196 = t298 * t276 + (-qJ(5) * t282 + t279 * t314) * t274;
t323 = t145 * t330 + t254 * t196;
t146 = qJ(5) * t256 + t257 * t314 + t298 * t335;
t164 = rSges(6,1) * t233 - rSges(6,2) * t232 + rSges(6,3) * t256;
t322 = -t146 - t164;
t310 = t275 * t333;
t179 = -pkin(3) * t310 + pkin(9) * t254 + t255 * t337;
t225 = pkin(3) * t329 + (-pkin(9) * t282 + t279 * t337) * t274;
t321 = t179 * t330 + t254 * t225;
t311 = t273 * t333;
t180 = pkin(3) * t311 + pkin(9) * t256 + t257 * t337;
t235 = t257 * pkin(2) + t256 * pkin(8);
t229 = t276 * t235;
t320 = t276 * t180 + t229;
t178 = rSges(5,1) * t241 + rSges(5,2) * t240 + rSges(5,3) * t256;
t319 = -t178 - t180;
t234 = t255 * pkin(2) + t254 * pkin(8);
t318 = -t179 - t234;
t207 = rSges(6,1) * t247 - rSges(6,2) * t246 - rSges(6,3) * t330;
t317 = -t196 - t207;
t177 = rSges(5,1) * t239 + rSges(5,2) * t238 + rSges(5,3) * t254;
t211 = rSges(5,1) * t253 + rSges(5,2) * t252 - rSges(5,3) * t330;
t128 = t177 * t330 + t254 * t211;
t316 = -t211 - t225;
t315 = t234 * t335 + t235 * t334;
t309 = t256 * t326 + t140;
t308 = -t146 - t325;
t307 = t276 * t146 + t320;
t306 = -t145 + t318;
t305 = -t180 + t322;
t304 = -t196 - t355;
t303 = -t225 + t317;
t300 = -t330 / 0.2e1;
t258 = t276 * t281 - t278 * t332;
t259 = t279 * t331 + t329;
t224 = rSges(4,1) * t259 + rSges(4,2) * t258 - rSges(4,3) * t330;
t260 = (pkin(2) * t279 - pkin(8) * t282) * t274;
t297 = (-t224 - t260) * t274;
t296 = -t180 + t308;
t295 = -t225 + t304;
t294 = t179 * t335 + t180 * t334 + t315;
t86 = t163 * t330 + t254 * t207 + t323;
t292 = (-t260 + t316) * t274;
t18 = t230 * t73 + t232 * t74 + t246 * t85;
t3 = t230 * t69 + t232 * t70 + t246 * t81;
t4 = t230 * t71 + t232 * t72 + t246 * t82;
t291 = t11 * t346 + t12 * t345 + t18 * t300 + t21 * t344 + t3 * t343 + t4 * t342;
t290 = t254 * t354 + t256 * t353;
t289 = (-t260 + t303) * t274;
t287 = t145 * t335 + t146 * t334 + t294;
t67 = t355 * t254 + t326 * t330 + t323;
t286 = (-t260 + t295) * t274;
t285 = -t330 * t350 + t290;
t284 = t352 * t343 + t351 * t342 + t350 * t339 + t353 * t335 / 0.2e1 - t354 * t334 / 0.2e1 + t349 * t300;
t251 = rSges(3,3) * t276 + (rSges(3,1) * t279 + rSges(3,2) * t282) * t274;
t250 = Icges(3,5) * t276 + (Icges(3,1) * t279 + Icges(3,4) * t282) * t274;
t249 = Icges(3,6) * t276 + (Icges(3,4) * t279 + Icges(3,2) * t282) * t274;
t248 = Icges(3,3) * t276 + (Icges(3,5) * t279 + Icges(3,6) * t282) * t274;
t245 = t257 * t281 + t311;
t244 = -t257 * t278 + t273 * t331;
t243 = t255 * t281 - t310;
t242 = -t255 * t278 - t275 * t331;
t223 = Icges(4,1) * t259 + Icges(4,4) * t258 - Icges(4,5) * t330;
t222 = Icges(4,4) * t259 + Icges(4,2) * t258 - Icges(4,6) * t330;
t221 = Icges(4,5) * t259 + Icges(4,6) * t258 - Icges(4,3) * t330;
t220 = rSges(3,1) * t257 - rSges(3,2) * t256 + rSges(3,3) * t335;
t219 = rSges(3,1) * t255 - rSges(3,2) * t254 - rSges(3,3) * t334;
t218 = Icges(3,1) * t257 - Icges(3,4) * t256 + Icges(3,5) * t335;
t217 = Icges(3,1) * t255 - Icges(3,4) * t254 - Icges(3,5) * t334;
t216 = Icges(3,4) * t257 - Icges(3,2) * t256 + Icges(3,6) * t335;
t215 = Icges(3,4) * t255 - Icges(3,2) * t254 - Icges(3,6) * t334;
t214 = Icges(3,5) * t257 - Icges(3,6) * t256 + Icges(3,3) * t335;
t213 = Icges(3,5) * t255 - Icges(3,6) * t254 - Icges(3,3) * t334;
t194 = -t219 * t276 - t251 * t334;
t193 = t220 * t276 - t251 * t335;
t190 = rSges(4,1) * t245 + rSges(4,2) * t244 + rSges(4,3) * t256;
t189 = rSges(4,1) * t243 + rSges(4,2) * t242 + rSges(4,3) * t254;
t188 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t256;
t187 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t254;
t186 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t256;
t185 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t254;
t184 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t256;
t183 = Icges(4,5) * t243 + Icges(4,6) * t242 + Icges(4,3) * t254;
t151 = (t219 * t273 + t220 * t275) * t274;
t150 = t256 * t179;
t149 = t256 * t177;
t139 = -t190 * t330 - t224 * t256;
t138 = t189 * t330 + t224 * t254;
t129 = -t178 * t330 - t211 * t256;
t126 = -t221 * t330 + t222 * t258 + t223 * t259;
t124 = t189 * t256 - t190 * t254;
t123 = (-t189 - t234) * t276 + t275 * t297;
t122 = t190 * t276 + t273 * t297 + t229;
t120 = t221 * t256 + t222 * t244 + t223 * t245;
t119 = t221 * t254 + t222 * t242 + t223 * t243;
t118 = -t178 * t254 + t149;
t114 = (t189 * t273 + t190 * t275) * t274 + t315;
t111 = -t184 * t330 + t186 * t258 + t188 * t259;
t110 = -t183 * t330 + t185 * t258 + t187 * t259;
t107 = t256 * t316 + t319 * t330;
t106 = t128 + t321;
t105 = t184 * t256 + t186 * t244 + t188 * t245;
t104 = t183 * t256 + t185 * t244 + t187 * t245;
t103 = t184 * t254 + t186 * t242 + t188 * t243;
t102 = t183 * t254 + t185 * t242 + t187 * t243;
t101 = t137 * t246 - t156 * t232;
t100 = -t136 * t246 + t156 * t230;
t93 = (-t177 + t318) * t276 + t275 * t292;
t92 = t178 * t276 + t273 * t292 + t320;
t87 = t256 * t317 + t322 * t330;
t84 = t136 * t232 - t137 * t230;
t83 = t254 * t319 + t149 + t150;
t80 = (t177 * t273 + t178 * t275) * t274 + t294;
t79 = t254 * t322 + t324;
t78 = t256 * t303 + t305 * t330;
t77 = t86 + t321;
t76 = (-t163 + t306) * t276 + t275 * t289;
t75 = t164 * t276 + t273 * t289 + t307;
t68 = t256 * t304 + t308 * t330;
t66 = t254 * t305 + t150 + t324;
t65 = (t163 * t273 + t164 * t275) * t274 + t287;
t64 = t126 * t276 + (-t110 * t275 + t111 * t273) * t274;
t63 = t110 * t254 + t111 * t256 - t126 * t330;
t62 = t256 * t295 + t296 * t330;
t61 = t67 + t321;
t60 = (t306 - t326) * t276 + t275 * t286;
t59 = t273 * t286 + t276 * t325 + t307;
t58 = t254 * t308 + t309;
t55 = t120 * t276 + (-t104 * t275 + t105 * t273) * t274;
t54 = t119 * t276 + (-t102 * t275 + t103 * t273) * t274;
t51 = t104 * t254 + t105 * t256 - t120 * t330;
t50 = t102 * t254 + t103 * t256 - t119 * t330;
t35 = t254 * t296 + t150 + t309;
t34 = (t273 * t326 + t275 * t325) * t274 + t287;
t1 = [m(2) + m(3) + m(4) + m(5) - t347; m(3) * t151 + m(4) * t114 + m(5) * t80 + m(6) * t65 + m(7) * t34; m(7) * (t34 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(6) * (t65 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t80 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(4) * (t114 ^ 2 + t122 ^ 2 + t123 ^ 2) + m(3) * (t151 ^ 2 + t193 ^ 2 + t194 ^ 2) + (t55 + (t214 * t335 - t256 * t216 + t257 * t218) * t335 + t351) * t335 + (-t54 + (-t213 * t334 - t215 * t254 + t217 * t255) * t334 + (-t213 * t335 + t214 * t334 + t215 * t256 + t216 * t254 - t217 * t257 - t218 * t255) * t335 - t352) * t334 + (-(-t248 * t334 - t249 * t254 + t250 * t255) * t334 + (t248 * t335 - t249 * t256 + t250 * t257) * t335 + t64 + ((t216 * t282 + t218 * t279) * t273 - (t215 * t282 + t217 * t279) * t275) * t348 + ((-t213 * t275 + t214 * t273 + t249 * t282 + t250 * t279) * t274 + t248 * t276) * t276 + t349) * t276; m(4) * t124 + m(5) * t83 + m(6) * t66 + m(7) * t35; m(7) * (t34 * t35 + t59 * t62 + t60 * t61) + m(6) * (t65 * t66 + t75 * t78 + t76 * t77) + m(5) * (t106 * t93 + t107 * t92 + t80 * t83) + m(4) * (t114 * t124 + t122 * t139 + t123 * t138) + (t51 * t341 + t50 * t340 - t282 * t64 / 0.2e1) * t274 + t63 * t339 + t284 + t54 * t343 + t55 * t342; t254 * t50 + t256 * t51 + (-t63 - t350) * t330 + m(7) * (t35 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t66 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t106 ^ 2 + t107 ^ 2 + t83 ^ 2) + m(4) * (t124 ^ 2 + t138 ^ 2 + t139 ^ 2) + t290; m(5) * t118 + m(6) * t79 + m(7) * t58; m(7) * (t34 * t58 + t59 * t68 + t60 * t67) + m(6) * (t65 * t79 + t75 * t87 + t76 * t86) + m(5) * (t118 * t80 + t128 * t93 + t129 * t92) + t284; m(7) * (t35 * t58 + t61 * t67 + t62 * t68) + m(6) * (t66 * t79 + t77 * t86 + t78 * t87) + m(5) * (t106 * t128 + t107 * t129 + t118 * t83) + t285; m(7) * (t58 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t79 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(5) * (t118 ^ 2 + t128 ^ 2 + t129 ^ 2) + t285; t347 * t330; m(7) * (t254 * t59 + t256 * t60 - t330 * t34) + m(6) * (t254 * t75 + t256 * t76 - t330 * t65); m(7) * (t254 * t62 + t256 * t61 - t330 * t35) + m(6) * (t254 * t78 + t256 * t77 - t330 * t66); m(7) * (t254 * t68 + t256 * t67 - t330 * t58) + m(6) * (t254 * t87 + t256 * t86 - t330 * t79); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t282 ^ 2 * t348 + t254 ^ 2 + t256 ^ 2); m(7) * t84; t23 * t344 + t18 * t339 + t15 * t346 + t16 * t345 + m(7) * (t100 * t60 + t101 * t59 + t34 * t84) + (t3 * t340 + t341 * t4) * t274; m(7) * (t100 * t61 + t101 * t62 + t35 * t84) + t291; m(7) * (t100 * t67 + t101 * t68 + t58 * t84) + t291; m(7) * (t100 * t256 + t101 * t254 - t330 * t84); t232 * t4 + t230 * t3 + t246 * t18 + m(7) * (t100 ^ 2 + t101 ^ 2 + t84 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
