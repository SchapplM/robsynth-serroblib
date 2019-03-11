% Calculate joint inertia matrix for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:21
% EndTime: 2019-03-09 18:44:34
% DurationCPUTime: 5.90s
% Computational Cost: add. (35151->698), mult. (53836->953), div. (0->0), fcn. (68118->14), ass. (0->330)
t362 = qJ(3) + pkin(12);
t310 = sin(t362);
t316 = cos(pkin(6));
t320 = sin(qJ(2));
t315 = sin(pkin(6));
t346 = cos(t362);
t334 = t315 * t346;
t281 = t316 * t310 + t320 * t334;
t314 = qJ(5) + qJ(6);
t311 = sin(t314);
t312 = cos(t314);
t324 = cos(qJ(2));
t385 = t315 * t324;
t250 = -t281 * t311 - t312 * t385;
t251 = t281 * t312 - t311 * t385;
t388 = t315 * t320;
t280 = t310 * t388 - t316 * t346;
t166 = Icges(7,5) * t251 + Icges(7,6) * t250 + Icges(7,3) * t280;
t167 = Icges(7,4) * t251 + Icges(7,2) * t250 + Icges(7,6) * t280;
t168 = Icges(7,1) * t251 + Icges(7,4) * t250 + Icges(7,5) * t280;
t85 = t280 * t166 + t250 * t167 + t251 * t168;
t318 = sin(qJ(5));
t322 = cos(qJ(5));
t261 = -t281 * t318 - t322 * t385;
t358 = t318 * t385;
t262 = t281 * t322 - t358;
t171 = Icges(6,5) * t262 + Icges(6,6) * t261 + Icges(6,3) * t280;
t172 = Icges(6,4) * t262 + Icges(6,2) * t261 + Icges(6,6) * t280;
t173 = Icges(6,1) * t262 + Icges(6,4) * t261 + Icges(6,5) * t280;
t89 = t280 * t171 + t261 * t172 + t262 * t173;
t404 = -t85 - t89;
t325 = cos(qJ(1));
t384 = t315 * t325;
t403 = t315 ^ 2;
t321 = sin(qJ(1));
t380 = t321 * t324;
t381 = t320 * t325;
t294 = t316 * t381 + t380;
t263 = t294 * t310 + t325 * t334;
t402 = t263 / 0.2e1;
t379 = t324 * t325;
t382 = t320 * t321;
t296 = -t316 * t382 + t379;
t265 = t296 * t310 - t321 * t334;
t401 = t265 / 0.2e1;
t400 = t280 / 0.2e1;
t293 = -t316 * t379 + t382;
t399 = t293 / 0.2e1;
t295 = t316 * t380 + t381;
t398 = t295 / 0.2e1;
t397 = t316 / 0.2e1;
t323 = cos(qJ(3));
t309 = pkin(3) * t323 + pkin(2);
t396 = -pkin(2) + t309;
t308 = pkin(5) * t322 + pkin(4);
t395 = -pkin(4) + t308;
t264 = t294 * t346 - t310 * t384;
t214 = -t264 * t311 + t293 * t312;
t215 = t264 * t312 + t293 * t311;
t131 = Icges(7,5) * t215 + Icges(7,6) * t214 + Icges(7,3) * t263;
t133 = Icges(7,4) * t215 + Icges(7,2) * t214 + Icges(7,6) * t263;
t135 = Icges(7,1) * t215 + Icges(7,4) * t214 + Icges(7,5) * t263;
t69 = t131 * t280 + t133 * t250 + t135 * t251;
t394 = t69 * t263;
t387 = t315 * t321;
t266 = t296 * t346 + t310 * t387;
t216 = -t266 * t311 + t295 * t312;
t217 = t266 * t312 + t295 * t311;
t132 = Icges(7,5) * t217 + Icges(7,6) * t216 + Icges(7,3) * t265;
t134 = Icges(7,4) * t217 + Icges(7,2) * t216 + Icges(7,6) * t265;
t136 = Icges(7,1) * t217 + Icges(7,4) * t216 + Icges(7,5) * t265;
t70 = t132 * t280 + t134 * t250 + t136 * t251;
t393 = t70 * t265;
t237 = Icges(3,5) * t294 - Icges(3,6) * t293 - Icges(3,3) * t384;
t392 = t237 * t325;
t317 = -qJ(4) - pkin(9);
t391 = t293 * t317;
t390 = t293 * t318;
t389 = t295 * t318;
t386 = t315 * t323;
t319 = sin(qJ(3));
t383 = t316 * t319;
t259 = t263 * pkin(10);
t326 = -pkin(11) - pkin(10);
t129 = pkin(5) * t390 - t263 * t326 + t264 * t395 - t259;
t335 = -t215 * rSges(7,1) - t214 * rSges(7,2);
t137 = t263 * rSges(7,3) - t335;
t378 = t129 + t137;
t208 = t266 * pkin(4) + pkin(10) * t265;
t351 = pkin(5) * t389 - t265 * t326 + t266 * t308;
t130 = -t208 + t351;
t138 = t217 * rSges(7,1) + t216 * rSges(7,2) + t265 * rSges(7,3);
t377 = t130 + t138;
t162 = -pkin(5) * t358 + t395 * t281 + (-pkin(10) - t326) * t280;
t169 = rSges(7,1) * t251 + rSges(7,2) * t250 + rSges(7,3) * t280;
t376 = t162 + t169;
t289 = t293 * pkin(9);
t356 = t319 * t384;
t300 = pkin(3) * t356;
t188 = t294 * t396 - t289 - t300 - t391;
t163 = t295 * t188;
t207 = t264 * pkin(4) + t259;
t375 = t295 * t207 + t163;
t235 = pkin(3) * t383 + ((pkin(9) + t317) * t324 + t396 * t320) * t315;
t374 = t188 * t385 + t293 * t235;
t254 = t296 * pkin(2) + pkin(9) * t295;
t357 = t319 * t387;
t349 = pkin(3) * t357 - t295 * t317 + t296 * t309;
t189 = -t254 + t349;
t249 = t316 * t254;
t373 = t316 * t189 + t249;
t186 = t266 * rSges(5,1) - t265 * rSges(5,2) + t295 * rSges(5,3);
t372 = -t186 - t189;
t253 = pkin(2) * t294 + t289;
t371 = -t188 - t253;
t370 = -t189 - t208;
t269 = -t294 * t319 - t323 * t384;
t270 = t294 * t323 - t356;
t197 = rSges(4,1) * t270 + rSges(4,2) * t269 + rSges(4,3) * t293;
t369 = -t197 - t253;
t228 = Icges(5,4) * t281 - Icges(5,2) * t280 - Icges(5,6) * t385;
t229 = Icges(5,1) * t281 - Icges(5,4) * t280 - Icges(5,5) * t385;
t368 = -t280 * t228 + t281 * t229;
t291 = t316 * t323 - t319 * t388;
t292 = t320 * t386 + t383;
t233 = Icges(4,4) * t292 + Icges(4,2) * t291 - Icges(4,6) * t385;
t234 = Icges(4,1) * t292 + Icges(4,4) * t291 - Icges(4,5) * t385;
t367 = t291 * t233 + t292 * t234;
t230 = rSges(5,1) * t281 - rSges(5,2) * t280 - rSges(5,3) * t385;
t366 = -t230 - t235;
t231 = t281 * pkin(4) + t280 * pkin(10);
t365 = -t231 - t235;
t364 = t253 * t387 + t254 * t384;
t363 = t325 * pkin(1) + pkin(8) * t387;
t82 = t85 * t280;
t26 = t393 + t82 + t394;
t57 = t131 * t263 + t133 * t214 + t135 * t215;
t58 = t132 * t263 + t134 * t214 + t136 * t215;
t77 = t166 * t263 + t167 * t214 + t168 * t215;
t7 = t263 * t57 + t265 * t58 + t280 * t77;
t59 = t131 * t265 + t133 * t216 + t135 * t217;
t60 = t132 * t265 + t134 * t216 + t136 * t217;
t78 = t166 * t265 + t167 * t216 + t168 * t217;
t8 = t263 * t59 + t265 * t60 + t280 * t78;
t361 = t280 * t26 + t263 * t7 + t265 * t8;
t223 = -t266 * t318 + t295 * t322;
t224 = t266 * t322 + t389;
t140 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t265;
t142 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t265;
t144 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t265;
t74 = t140 * t280 + t142 * t261 + t144 * t262;
t81 = t171 * t265 + t172 * t223 + t173 * t224;
t360 = t74 / 0.2e1 + t81 / 0.2e1;
t221 = -t264 * t318 + t293 * t322;
t222 = t264 * t322 + t390;
t139 = Icges(6,5) * t222 + Icges(6,6) * t221 + Icges(6,3) * t263;
t141 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t263;
t143 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t263;
t73 = t139 * t280 + t141 * t261 + t143 * t262;
t80 = t171 * t263 + t172 * t221 + t173 * t222;
t359 = t80 / 0.2e1 + t73 / 0.2e1;
t146 = t224 * rSges(6,1) + t223 * rSges(6,2) + t265 * rSges(6,3);
t355 = -t146 + t370;
t174 = rSges(6,1) * t262 + rSges(6,2) * t261 + rSges(6,3) * t280;
t354 = -t174 + t365;
t353 = t316 * t208 + t373;
t352 = -t207 + t371;
t271 = -t296 * t319 + t321 * t386;
t272 = t296 * t323 + t357;
t198 = t272 * rSges(4,1) + t271 * rSges(4,2) + t295 * rSges(4,3);
t277 = Icges(3,3) * t316 + (Icges(3,5) * t320 + Icges(3,6) * t324) * t315;
t278 = Icges(3,6) * t316 + (Icges(3,4) * t320 + Icges(3,2) * t324) * t315;
t279 = Icges(3,5) * t316 + (Icges(3,1) * t320 + Icges(3,4) * t324) * t315;
t350 = t316 * t277 + t278 * t385 + t279 * t388;
t244 = t296 * rSges(3,1) - t295 * rSges(3,2) + rSges(3,3) * t387;
t348 = -t385 / 0.2e1;
t347 = -t321 * pkin(1) + pkin(8) * t384;
t236 = rSges(4,1) * t292 + rSges(4,2) * t291 - rSges(4,3) * t385;
t297 = (pkin(2) * t320 - pkin(9) * t324) * t315;
t345 = t315 * (-t236 - t297);
t344 = t370 - t377;
t343 = t365 - t376;
t342 = t188 * t387 + t189 * t384 + t364;
t341 = t207 * t385 + t293 * t231 + t374;
t340 = t394 / 0.2e1 + t393 / 0.2e1 + t77 * t402 + t78 * t401 + t82;
t339 = t315 * (-t297 + t366);
t13 = t293 * t57 + t295 * t58 - t385 * t77;
t14 = t293 * t59 + t295 * t60 - t385 * t78;
t28 = t69 * t293 + t70 * t295 - t385 * t85;
t338 = t13 * t402 + t14 * t401 + t26 * t348 + t28 * t400 + t8 * t398 + t7 * t399;
t15 = t77 * t316 + (t321 * t58 - t325 * t57) * t315;
t16 = t78 * t316 + (t321 * t60 - t325 * t59) * t315;
t84 = t85 * t316;
t30 = t84 + (t70 * t321 - t69 * t325) * t315;
t337 = t15 * t402 + t16 * t401 + t26 * t397 + t30 * t400 + t8 * t387 / 0.2e1 - t7 * t384 / 0.2e1;
t336 = -t264 * rSges(5,1) + t263 * rSges(5,2);
t333 = t349 + t363;
t332 = t315 * (-t297 + t354);
t331 = t207 * t387 + t208 * t384 + t342;
t330 = t315 * (-t297 + t343);
t329 = -t294 * t309 + t300 + t347;
t145 = rSges(6,1) * t222 + rSges(6,2) * t221 + rSges(6,3) * t263;
t243 = rSges(3,1) * t294 - rSges(3,2) * t293 - rSges(3,3) * t384;
t180 = Icges(5,5) * t266 - Icges(5,6) * t265 + Icges(5,3) * t295;
t182 = Icges(5,4) * t266 - Icges(5,2) * t265 + Icges(5,6) * t295;
t184 = Icges(5,1) * t266 - Icges(5,4) * t265 + Icges(5,5) * t295;
t104 = -t180 * t385 - t182 * t280 + t184 * t281;
t192 = Icges(4,5) * t272 + Icges(4,6) * t271 + Icges(4,3) * t295;
t194 = Icges(4,4) * t272 + Icges(4,2) * t271 + Icges(4,6) * t295;
t196 = Icges(4,1) * t272 + Icges(4,4) * t271 + Icges(4,5) * t295;
t114 = -t192 * t385 + t194 * t291 + t196 * t292;
t227 = Icges(5,5) * t281 - Icges(5,6) * t280 - Icges(5,3) * t385;
t117 = t227 * t295 - t228 * t265 + t229 * t266;
t232 = Icges(4,5) * t292 + Icges(4,6) * t291 - Icges(4,3) * t385;
t121 = t232 * t295 + t233 * t271 + t234 * t272;
t328 = t70 / 0.2e1 + t78 / 0.2e1 + t121 / 0.2e1 + t117 / 0.2e1 + t114 / 0.2e1 + t104 / 0.2e1 + t360;
t179 = Icges(5,5) * t264 - Icges(5,6) * t263 + Icges(5,3) * t293;
t181 = Icges(5,4) * t264 - Icges(5,2) * t263 + Icges(5,6) * t293;
t183 = Icges(5,1) * t264 - Icges(5,4) * t263 + Icges(5,5) * t293;
t103 = -t179 * t385 - t181 * t280 + t183 * t281;
t191 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t293;
t193 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t293;
t195 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t293;
t113 = -t191 * t385 + t193 * t291 + t195 * t292;
t116 = t227 * t293 - t228 * t263 + t229 * t264;
t120 = t232 * t293 + t233 * t269 + t234 * t270;
t327 = t77 / 0.2e1 + t69 / 0.2e1 + t113 / 0.2e1 + t103 / 0.2e1 + t116 / 0.2e1 + t120 / 0.2e1 + t359;
t302 = rSges(2,1) * t325 - rSges(2,2) * t321;
t301 = -rSges(2,1) * t321 - rSges(2,2) * t325;
t282 = t316 * rSges(3,3) + (rSges(3,1) * t320 + rSges(3,2) * t324) * t315;
t242 = Icges(3,1) * t296 - Icges(3,4) * t295 + Icges(3,5) * t387;
t241 = Icges(3,1) * t294 - Icges(3,4) * t293 - Icges(3,5) * t384;
t240 = Icges(3,4) * t296 - Icges(3,2) * t295 + Icges(3,6) * t387;
t239 = Icges(3,4) * t294 - Icges(3,2) * t293 - Icges(3,6) * t384;
t238 = Icges(3,5) * t296 - Icges(3,6) * t295 + Icges(3,3) * t387;
t226 = t244 + t363;
t225 = -t243 + t347;
t206 = -t243 * t316 - t282 * t384;
t205 = t244 * t316 - t282 * t387;
t190 = t350 * t316;
t185 = rSges(5,3) * t293 - t336;
t170 = (t243 * t321 + t244 * t325) * t315;
t165 = t277 * t387 - t278 * t295 + t279 * t296;
t164 = -t277 * t384 - t278 * t293 + t279 * t294;
t159 = t254 + t198 + t363;
t158 = t347 + t369;
t155 = t263 * t169;
t152 = -t198 * t385 - t236 * t295;
t151 = t197 * t385 + t236 * t293;
t150 = t316 * t238 + (t240 * t324 + t242 * t320) * t315;
t149 = t316 * t237 + (t239 * t324 + t241 * t320) * t315;
t148 = t333 + t186;
t147 = (-rSges(5,3) + t317) * t293 + t329 + t336;
t128 = -t232 * t385 + t367;
t127 = t128 * t316;
t126 = t280 * t138;
t125 = t197 * t295 - t198 * t293;
t124 = t265 * t137;
t123 = t316 * t369 + t325 * t345;
t122 = t316 * t198 + t321 * t345 + t249;
t119 = -t227 * t385 + t368;
t118 = t119 * t316;
t115 = (t197 * t321 + t198 * t325) * t315 + t364;
t112 = t208 + t333 + t146;
t111 = -t145 - t207 + t329 + t391;
t110 = t146 * t280 - t174 * t265;
t109 = -t145 * t280 + t174 * t263;
t108 = t192 * t295 + t194 * t271 + t196 * t272;
t107 = t191 * t295 + t193 * t271 + t195 * t272;
t106 = t192 * t293 + t194 * t269 + t196 * t270;
t105 = t191 * t293 + t193 * t269 + t195 * t270;
t102 = -t169 * t265 + t126;
t101 = -t137 * t280 + t155;
t100 = t295 * t366 + t372 * t385;
t99 = t185 * t385 + t230 * t293 + t374;
t98 = t333 + t351 + t138;
t97 = -t264 * t308 + (-pkin(5) * t318 + t317) * t293 + (-rSges(7,3) + t326) * t263 + t329 + t335;
t96 = t180 * t295 - t182 * t265 + t184 * t266;
t95 = t179 * t295 - t181 * t265 + t183 * t266;
t94 = t180 * t293 - t182 * t263 + t184 * t264;
t93 = t179 * t293 - t181 * t263 + t183 * t264;
t92 = (-t185 + t371) * t316 + t325 * t339;
t91 = t316 * t186 + t321 * t339 + t373;
t90 = t145 * t265 - t146 * t263;
t88 = t89 * t316;
t87 = -t138 * t263 + t124;
t86 = t89 * t280;
t83 = t295 * t185 + t293 * t372 + t163;
t79 = (t185 * t321 + t186 * t325) * t315 + t342;
t72 = t295 * t354 + t355 * t385;
t71 = t145 * t385 + t174 * t293 + t341;
t68 = (-t145 + t352) * t316 + t325 * t332;
t67 = t316 * t146 + t321 * t332 + t353;
t66 = t140 * t265 + t142 * t223 + t144 * t224;
t65 = t139 * t265 + t141 * t223 + t143 * t224;
t64 = t140 * t263 + t142 * t221 + t144 * t222;
t63 = t139 * t263 + t141 * t221 + t143 * t222;
t56 = t280 * t130 - t265 * t376 + t126;
t55 = t263 * t162 - t280 * t378 + t155;
t54 = t295 * t145 + t293 * t355 + t375;
t53 = (t145 * t321 + t146 * t325) * t315 + t331;
t52 = t127 + (-t113 * t325 + t114 * t321) * t315;
t51 = t113 * t293 + t114 * t295 - t128 * t385;
t50 = t129 * t265 - t263 * t377 + t124;
t49 = t121 * t316 + (-t107 * t325 + t108 * t321) * t315;
t48 = t120 * t316 + (-t105 * t325 + t106 * t321) * t315;
t47 = t118 + (-t103 * t325 + t104 * t321) * t315;
t46 = t295 * t343 + t344 * t385;
t45 = t293 * t376 + t378 * t385 + t341;
t44 = t107 * t293 + t108 * t295 - t121 * t385;
t43 = t105 * t293 + t106 * t295 - t120 * t385;
t42 = (t352 - t378) * t316 + t325 * t330;
t41 = t316 * t377 + t321 * t330 + t353;
t40 = t103 * t293 + t104 * t295 - t119 * t385;
t39 = t117 * t316 + (t321 * t96 - t325 * t95) * t315;
t38 = t116 * t316 + (t321 * t94 - t325 * t93) * t315;
t37 = -t117 * t385 + t293 * t95 + t295 * t96;
t36 = -t116 * t385 + t293 * t93 + t295 * t94;
t35 = t293 * t344 + t295 * t378 + t375;
t34 = (t321 * t378 + t325 * t377) * t315 + t331;
t33 = t88 + (t74 * t321 - t73 * t325) * t315;
t32 = t73 * t293 + t74 * t295 - t385 * t89;
t31 = t73 * t263 + t74 * t265 + t86;
t23 = t81 * t316 + (t321 * t66 - t325 * t65) * t315;
t22 = t80 * t316 + (t321 * t64 - t325 * t63) * t315;
t20 = t293 * t65 + t295 * t66 - t385 * t81;
t19 = t293 * t63 + t295 * t64 - t385 * t80;
t18 = t263 * t65 + t265 * t66 + t280 * t81;
t17 = t263 * t63 + t265 * t64 + t280 * t80;
t1 = [(-t227 - t232) * t385 + m(7) * (t97 ^ 2 + t98 ^ 2) + m(6) * (t111 ^ 2 + t112 ^ 2) + m(5) * (t147 ^ 2 + t148 ^ 2) + m(4) * (t158 ^ 2 + t159 ^ 2) + m(3) * (t225 ^ 2 + t226 ^ 2) + m(2) * (t301 ^ 2 + t302 ^ 2) + t350 + Icges(2,3) + t367 + t368 - t404; t84 + t190 + t118 + t127 + t88 + m(7) * (t41 * t98 + t42 * t97) + m(6) * (t111 * t68 + t112 * t67) + m(5) * (t147 * t92 + t148 * t91) + m(4) * (t122 * t159 + t123 * t158) + m(3) * (t205 * t226 + t206 * t225) + ((-t149 / 0.2e1 - t164 / 0.2e1 - t327) * t325 + (t165 / 0.2e1 + t150 / 0.2e1 + t328) * t321) * t315; (t30 + t33 + t52 + t47 + t190) * t316 + m(7) * (t34 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t53 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(5) * (t79 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t115 ^ 2 + t122 ^ 2 + t123 ^ 2) + m(3) * (t170 ^ 2 + t205 ^ 2 + t206 ^ 2) + ((-t15 - t22 - t48 - t38 + ((-t239 * t293 + t241 * t294) * t315 - t403 * t392) * t325) * t325 + (t16 + t23 + t49 + t39 + ((-t240 * t295 + t242 * t296 + (t238 * t321 - t392) * t315) * t321 + (t238 * t384 + t239 * t295 + t240 * t293 - t241 * t296 - t242 * t294) * t325) * t315) * t321 + ((-t149 - t164) * t325 + (t150 + t165) * t321) * t316) * t315; (-t119 - t128 + t404) * t385 + m(7) * (t45 * t97 + t46 * t98) + m(6) * (t111 * t71 + t112 * t72) + m(5) * (t100 * t148 + t147 * t99) + m(4) * (t151 * t158 + t152 * t159) + t328 * t295 + t327 * t293; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1 + t51 / 0.2e1) * t316 + (t16 / 0.2e1 + t23 / 0.2e1 + t39 / 0.2e1 + t49 / 0.2e1) * t295 + (t15 / 0.2e1 + t22 / 0.2e1 + t38 / 0.2e1 + t48 / 0.2e1) * t293 + m(5) * (t100 * t91 + t79 * t83 + t92 * t99) + m(7) * (t34 * t35 + t41 * t46 + t42 * t45) + m(6) * (t53 * t54 + t67 * t72 + t68 * t71) + m(4) * (t115 * t125 + t122 * t152 + t123 * t151) + ((-t36 / 0.2e1 - t43 / 0.2e1 - t13 / 0.2e1 - t19 / 0.2e1) * t325 + (-t47 / 0.2e1 - t52 / 0.2e1 - t30 / 0.2e1 - t33 / 0.2e1) * t324 + (t37 / 0.2e1 + t44 / 0.2e1 + t14 / 0.2e1 + t20 / 0.2e1) * t321) * t315; (-t28 - t32 - t40 - t51) * t385 + (t14 + t20 + t44 + t37) * t295 + (t13 + t19 + t43 + t36) * t293 + m(7) * (t35 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t54 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t100 ^ 2 + t83 ^ 2 + t99 ^ 2) + m(4) * (t125 ^ 2 + t151 ^ 2 + t152 ^ 2); m(7) * (t293 * t98 + t295 * t97) + m(6) * (t111 * t295 + t112 * t293) + m(5) * (t147 * t295 + t148 * t293); m(7) * (t293 * t41 + t295 * t42 - t34 * t385) + m(6) * (t293 * t67 + t295 * t68 - t385 * t53) + m(5) * (t293 * t91 + t295 * t92 - t385 * t79); m(7) * (t293 * t46 + t295 * t45 - t35 * t385) + m(6) * (t293 * t72 + t295 * t71 - t385 * t54) + m(5) * (t100 * t293 + t295 * t99 - t385 * t83); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t324 ^ 2 * t403 + t293 ^ 2 + t295 ^ 2); t86 + t360 * t265 + t359 * t263 + m(7) * (t55 * t97 + t56 * t98) + m(6) * (t109 * t111 + t110 * t112) + t340; t31 * t397 + t33 * t400 + t22 * t402 + t23 * t401 + (t321 * t18 / 0.2e1 - t325 * t17 / 0.2e1) * t315 + m(7) * (t34 * t50 + t41 * t56 + t42 * t55) + m(6) * (t109 * t68 + t110 * t67 + t53 * t90) + t337; t31 * t348 + t20 * t401 + t17 * t399 + t18 * t398 + t19 * t402 + t32 * t400 + m(7) * (t35 * t50 + t45 * t55 + t46 * t56) + m(6) * (t109 * t71 + t110 * t72 + t54 * t90) + t338; m(6) * (t109 * t295 + t110 * t293 - t385 * t90) + m(7) * (t293 * t56 + t295 * t55 - t385 * t50); t263 * t17 + t265 * t18 + t280 * t31 + m(7) * (t50 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2 + t90 ^ 2) + t361; m(7) * (t101 * t97 + t102 * t98) + t340; m(7) * (t101 * t42 + t102 * t41 + t34 * t87) + t337; m(7) * (t101 * t45 + t102 * t46 + t35 * t87) + t338; m(7) * (t101 * t295 + t102 * t293 - t385 * t87); m(7) * (t101 * t55 + t102 * t56 + t50 * t87) + t361; m(7) * (t101 ^ 2 + t102 ^ 2 + t87 ^ 2) + t361;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
