% Calculate joint inertia matrix for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:38
% EndTime: 2019-03-09 00:42:53
% DurationCPUTime: 7.00s
% Computational Cost: add. (43894->603), mult. (66614->856), div. (0->0), fcn. (84964->14), ass. (0->298)
t296 = sin(pkin(12));
t298 = cos(pkin(12));
t305 = cos(qJ(2));
t299 = cos(pkin(6));
t302 = sin(qJ(2));
t354 = t299 * t302;
t281 = t296 * t305 + t298 * t354;
t352 = qJ(3) + qJ(4);
t293 = sin(t352);
t297 = sin(pkin(6));
t323 = cos(t352);
t317 = t297 * t323;
t264 = t281 * t293 + t298 * t317;
t360 = t297 * t298;
t265 = t281 * t323 - t293 * t360;
t353 = t299 * t305;
t280 = t296 * t302 - t298 * t353;
t300 = sin(qJ(5));
t363 = t280 * t300;
t303 = cos(qJ(5));
t366 = pkin(5) * t303;
t156 = pkin(5) * t363 + pkin(11) * t264 + t265 * t366;
t295 = qJ(5) + qJ(6);
t292 = sin(t295);
t294 = cos(t295);
t229 = -t265 * t292 + t280 * t294;
t230 = t265 * t294 + t280 * t292;
t167 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t264;
t350 = t156 + t167;
t358 = t297 * t302;
t278 = t293 * t358 - t299 * t323;
t279 = t299 * t293 + t302 * t317;
t356 = t297 * t305;
t334 = t300 * t356;
t187 = -pkin(5) * t334 + pkin(11) * t278 + t279 * t366;
t260 = -t279 * t292 - t294 * t356;
t261 = t279 * t294 - t292 * t356;
t191 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t278;
t382 = t187 + t191;
t201 = Icges(5,5) * t265 - Icges(5,6) * t264 + Icges(5,3) * t280;
t203 = Icges(5,4) * t265 - Icges(5,2) * t264 + Icges(5,6) * t280;
t205 = Icges(5,1) * t265 - Icges(5,4) * t264 + Icges(5,5) * t280;
t121 = t201 * t280 - t203 * t264 + t205 * t265;
t283 = -t296 * t354 + t298 * t305;
t266 = t283 * t293 - t296 * t317;
t361 = t296 * t297;
t267 = t283 * t323 + t293 * t361;
t282 = t296 * t353 + t298 * t302;
t202 = Icges(5,5) * t267 - Icges(5,6) * t266 + Icges(5,3) * t282;
t204 = Icges(5,4) * t267 - Icges(5,2) * t266 + Icges(5,6) * t282;
t206 = Icges(5,1) * t267 - Icges(5,4) * t266 + Icges(5,5) * t282;
t122 = t202 * t280 - t204 * t264 + t206 * t265;
t239 = Icges(5,5) * t279 - Icges(5,6) * t278 - Icges(5,3) * t356;
t240 = Icges(5,4) * t279 - Icges(5,2) * t278 - Icges(5,6) * t356;
t241 = Icges(5,1) * t279 - Icges(5,4) * t278 - Icges(5,5) * t356;
t140 = t239 * t280 - t240 * t264 + t241 * t265;
t188 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t278;
t189 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t278;
t190 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t278;
t108 = t188 * t264 + t189 * t229 + t190 * t230;
t161 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t264;
t163 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t264;
t165 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t264;
t88 = t161 * t264 + t163 * t229 + t165 * t230;
t231 = -t267 * t292 + t282 * t294;
t232 = t267 * t294 + t282 * t292;
t162 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t266;
t164 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t266;
t166 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t266;
t89 = t162 * t264 + t164 * t229 + t166 * t230;
t17 = -t108 * t356 + t280 * t88 + t282 * t89;
t268 = -t279 * t300 - t303 * t356;
t269 = t279 * t303 - t334;
t197 = Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t278;
t198 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t278;
t199 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t278;
t235 = -t265 * t300 + t280 * t303;
t236 = t265 * t303 + t363;
t110 = t197 * t264 + t198 * t235 + t199 * t236;
t170 = Icges(6,5) * t236 + Icges(6,6) * t235 + Icges(6,3) * t264;
t172 = Icges(6,4) * t236 + Icges(6,2) * t235 + Icges(6,6) * t264;
t174 = Icges(6,1) * t236 + Icges(6,4) * t235 + Icges(6,5) * t264;
t94 = t170 * t264 + t172 * t235 + t174 * t236;
t237 = -t267 * t300 + t282 * t303;
t362 = t282 * t300;
t238 = t267 * t303 + t362;
t171 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t266;
t173 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t266;
t175 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t266;
t95 = t171 * t264 + t173 * t235 + t175 * t236;
t33 = -t110 * t356 + t280 * t94 + t282 * t95;
t381 = t121 * t280 + t122 * t282 - t140 * t356 + t17 + t33;
t123 = t201 * t282 - t203 * t266 + t205 * t267;
t124 = t202 * t282 - t204 * t266 + t206 * t267;
t141 = t239 * t282 - t240 * t266 + t241 * t267;
t109 = t188 * t266 + t189 * t231 + t190 * t232;
t90 = t161 * t266 + t163 * t231 + t165 * t232;
t91 = t162 * t266 + t164 * t231 + t166 * t232;
t18 = -t109 * t356 + t280 * t90 + t282 * t91;
t111 = t197 * t266 + t198 * t237 + t199 * t238;
t96 = t170 * t266 + t172 * t237 + t174 * t238;
t97 = t171 * t266 + t173 * t237 + t175 * t238;
t34 = -t111 * t356 + t280 * t96 + t282 * t97;
t380 = t123 * t280 + t124 * t282 - t141 * t356 + t18 + t34;
t21 = t108 * t299 + (t296 * t89 - t298 * t88) * t297;
t37 = t110 * t299 + (t296 * t95 - t298 * t94) * t297;
t379 = t21 + t37 + t140 * t299 + (-t121 * t298 + t122 * t296) * t297;
t22 = t109 * t299 + (t296 * t91 - t298 * t90) * t297;
t38 = t111 * t299 + (t296 * t97 - t298 * t96) * t297;
t378 = t22 + t38 + t141 * t299 + (-t123 * t298 + t124 * t296) * t297;
t131 = -t201 * t356 - t203 * t278 + t205 * t279;
t132 = -t202 * t356 - t204 * t278 + t206 * t279;
t145 = -t239 * t356 - t240 * t278 + t241 * t279;
t100 = t161 * t278 + t163 * t260 + t165 * t261;
t101 = t162 * t278 + t164 * t260 + t166 * t261;
t115 = t188 * t278 + t189 * t260 + t190 * t261;
t45 = t100 * t280 + t101 * t282 - t115 * t356;
t104 = t170 * t278 + t172 * t268 + t174 * t269;
t105 = t171 * t278 + t173 * t268 + t175 * t269;
t118 = t197 * t278 + t198 * t268 + t199 * t269;
t53 = t104 * t280 + t105 * t282 - t118 * t356;
t377 = t131 * t280 + t132 * t282 - t145 * t356 + t45 + t53;
t48 = t115 * t299 + (-t100 * t298 + t101 * t296) * t297;
t55 = t118 * t299 + (-t104 * t298 + t105 * t296) * t297;
t376 = t48 + t55 + t145 * t299 + (-t131 * t298 + t132 * t296) * t297;
t375 = t264 / 0.2e1;
t374 = t266 / 0.2e1;
t373 = t278 / 0.2e1;
t372 = t280 / 0.2e1;
t371 = t282 / 0.2e1;
t370 = t296 / 0.2e1;
t369 = -t298 / 0.2e1;
t368 = t299 / 0.2e1;
t304 = cos(qJ(3));
t367 = pkin(3) * t304;
t301 = sin(qJ(3));
t359 = t297 * t301;
t357 = t297 * t304;
t355 = t299 * t301;
t176 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t264;
t226 = t265 * pkin(4) + t264 * pkin(10);
t211 = t282 * t226;
t351 = t282 * t176 + t211;
t157 = pkin(5) * t362 + pkin(11) * t266 + t267 * t366;
t168 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t266;
t349 = t157 + t168;
t177 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t266;
t227 = t267 * pkin(4) + t266 * pkin(10);
t348 = -t177 - t227;
t335 = t298 * t359;
t209 = -pkin(3) * t335 + pkin(9) * t280 + t281 * t367;
t256 = pkin(3) * t355 + (-pkin(9) * t305 + t302 * t367) * t297;
t346 = t209 * t356 + t280 * t256;
t336 = t296 * t359;
t210 = pkin(3) * t336 + pkin(9) * t282 + t283 * t367;
t263 = t283 * pkin(2) + t282 * pkin(8);
t259 = t299 * t263;
t345 = t299 * t210 + t259;
t200 = rSges(6,1) * t269 + rSges(6,2) * t268 + rSges(6,3) * t278;
t255 = t279 * pkin(4) + t278 * pkin(10);
t344 = -t200 - t255;
t208 = rSges(5,1) * t267 - rSges(5,2) * t266 + rSges(5,3) * t282;
t343 = -t208 - t210;
t262 = t281 * pkin(2) + t280 * pkin(8);
t342 = -t209 - t262;
t341 = t226 * t356 + t280 * t255;
t207 = rSges(5,1) * t265 - rSges(5,2) * t264 + rSges(5,3) * t280;
t242 = rSges(5,1) * t279 - rSges(5,2) * t278 - rSges(5,3) * t356;
t159 = t207 * t356 + t280 * t242;
t340 = -t242 - t256;
t339 = t262 * t361 + t263 * t360;
t42 = t100 * t264 + t101 * t266 + t115 * t278;
t7 = t108 * t278 + t264 * t88 + t266 * t89;
t8 = t109 * t278 + t264 * t90 + t266 * t91;
t338 = t264 * t7 + t266 * t8 + t278 * t42;
t333 = t282 * t350 + t211;
t332 = -t227 - t349;
t331 = -t210 + t348;
t330 = -t255 - t382;
t329 = t299 * t227 + t345;
t328 = -t256 + t344;
t327 = -t226 + t342;
t326 = t361 / 0.2e1;
t325 = -t360 / 0.2e1;
t324 = -t356 / 0.2e1;
t284 = t299 * t304 - t301 * t358;
t285 = t302 * t357 + t355;
t254 = rSges(4,1) * t285 + rSges(4,2) * t284 - rSges(4,3) * t356;
t286 = (pkin(2) * t302 - pkin(8) * t305) * t297;
t322 = (-t254 - t286) * t297;
t321 = -t210 + t332;
t320 = -t256 + t330;
t319 = t209 * t361 + t210 * t360 + t339;
t113 = t176 * t356 + t280 * t200 + t341;
t318 = (-t286 + t340) * t297;
t316 = t17 * t375 + t18 * t374 + t42 * t324 + t8 * t371 + t7 * t372 + t45 * t373;
t315 = t21 * t375 + t22 * t374 + t7 * t325 + t8 * t326 + t42 * t368 + t48 * t373;
t314 = t280 * t381 + t282 * t380;
t313 = (-t286 + t328) * t297;
t312 = t226 * t361 + t227 * t360 + t319;
t86 = t382 * t280 + t350 * t356 + t341;
t311 = (-t286 + t320) * t297;
t25 = t110 * t278 + t264 * t94 + t266 * t95;
t26 = t111 * t278 + t264 * t96 + t266 * t97;
t50 = t104 * t264 + t105 * t266 + t118 * t278;
t310 = t25 * t372 + t26 * t371 + t50 * t324 + t33 * t375 + t34 * t374 + t53 * t373 + t316;
t309 = -t356 * t377 + t314;
t308 = t324 * t376 + t325 * t381 + t326 * t380 + t368 * t377 + t371 * t378 + t372 * t379;
t277 = t299 * rSges(3,3) + (rSges(3,1) * t302 + rSges(3,2) * t305) * t297;
t276 = Icges(3,5) * t299 + (Icges(3,1) * t302 + Icges(3,4) * t305) * t297;
t275 = Icges(3,6) * t299 + (Icges(3,4) * t302 + Icges(3,2) * t305) * t297;
t274 = Icges(3,3) * t299 + (Icges(3,5) * t302 + Icges(3,6) * t305) * t297;
t273 = t283 * t304 + t336;
t272 = -t283 * t301 + t296 * t357;
t271 = t281 * t304 - t335;
t270 = -t281 * t301 - t298 * t357;
t253 = Icges(4,1) * t285 + Icges(4,4) * t284 - Icges(4,5) * t356;
t252 = Icges(4,4) * t285 + Icges(4,2) * t284 - Icges(4,6) * t356;
t251 = Icges(4,5) * t285 + Icges(4,6) * t284 - Icges(4,3) * t356;
t250 = rSges(3,1) * t283 - rSges(3,2) * t282 + rSges(3,3) * t361;
t249 = rSges(3,1) * t281 - rSges(3,2) * t280 - rSges(3,3) * t360;
t248 = Icges(3,1) * t283 - Icges(3,4) * t282 + Icges(3,5) * t361;
t247 = Icges(3,1) * t281 - Icges(3,4) * t280 - Icges(3,5) * t360;
t246 = Icges(3,4) * t283 - Icges(3,2) * t282 + Icges(3,6) * t361;
t245 = Icges(3,4) * t281 - Icges(3,2) * t280 - Icges(3,6) * t360;
t244 = Icges(3,5) * t283 - Icges(3,6) * t282 + Icges(3,3) * t361;
t243 = Icges(3,5) * t281 - Icges(3,6) * t280 - Icges(3,3) * t360;
t224 = -t249 * t299 - t277 * t360;
t223 = t250 * t299 - t277 * t361;
t219 = rSges(4,1) * t273 + rSges(4,2) * t272 + rSges(4,3) * t282;
t218 = rSges(4,1) * t271 + rSges(4,2) * t270 + rSges(4,3) * t280;
t217 = Icges(4,1) * t273 + Icges(4,4) * t272 + Icges(4,5) * t282;
t216 = Icges(4,1) * t271 + Icges(4,4) * t270 + Icges(4,5) * t280;
t215 = Icges(4,4) * t273 + Icges(4,2) * t272 + Icges(4,6) * t282;
t214 = Icges(4,4) * t271 + Icges(4,2) * t270 + Icges(4,6) * t280;
t213 = Icges(4,5) * t273 + Icges(4,6) * t272 + Icges(4,3) * t282;
t212 = Icges(4,5) * t271 + Icges(4,6) * t270 + Icges(4,3) * t280;
t186 = (t249 * t296 + t250 * t298) * t297;
t185 = t282 * t209;
t184 = t282 * t207;
t180 = t264 * t191;
t179 = -t219 * t356 - t254 * t282;
t178 = t218 * t356 + t254 * t280;
t160 = -t208 * t356 - t282 * t242;
t154 = -t251 * t356 + t252 * t284 + t253 * t285;
t151 = t278 * t168;
t149 = t266 * t167;
t148 = t218 * t282 - t219 * t280;
t147 = (-t218 - t262) * t299 + t298 * t322;
t146 = t299 * t219 + t296 * t322 + t259;
t144 = t251 * t282 + t252 * t272 + t253 * t273;
t143 = t251 * t280 + t252 * t270 + t253 * t271;
t142 = -t208 * t280 + t184;
t139 = (t218 * t296 + t219 * t298) * t297 + t339;
t138 = -t213 * t356 + t215 * t284 + t217 * t285;
t137 = -t212 * t356 + t214 * t284 + t216 * t285;
t136 = t177 * t278 - t200 * t266;
t135 = -t176 * t278 + t200 * t264;
t134 = -t191 * t266 + t151;
t133 = -t167 * t278 + t180;
t130 = t282 * t340 + t343 * t356;
t129 = t159 + t346;
t128 = t213 * t282 + t215 * t272 + t217 * t273;
t127 = t212 * t282 + t214 * t272 + t216 * t273;
t126 = t213 * t280 + t215 * t270 + t217 * t271;
t125 = t212 * t280 + t214 * t270 + t216 * t271;
t120 = (-t207 + t342) * t299 + t298 * t318;
t119 = t299 * t208 + t296 * t318 + t345;
t117 = t176 * t266 - t177 * t264;
t116 = -t168 * t264 + t149;
t114 = t282 * t344 + t348 * t356;
t112 = t280 * t343 + t184 + t185;
t107 = (t207 * t296 + t208 * t298) * t297 + t319;
t106 = t280 * t348 + t351;
t103 = t282 * t328 + t331 * t356;
t102 = t113 + t346;
t99 = (-t176 + t327) * t299 + t298 * t313;
t98 = t299 * t177 + t296 * t313 + t329;
t93 = t278 * t157 - t266 * t382 + t151;
t92 = t264 * t187 - t278 * t350 + t180;
t87 = t282 * t330 + t332 * t356;
t85 = t280 * t331 + t185 + t351;
t84 = (t176 * t296 + t177 * t298) * t297 + t312;
t83 = t154 * t299 + (-t137 * t298 + t138 * t296) * t297;
t82 = t137 * t280 + t138 * t282 - t154 * t356;
t81 = t266 * t156 - t264 * t349 + t149;
t78 = t282 * t320 + t321 * t356;
t77 = t86 + t346;
t76 = t144 * t299 + (-t127 * t298 + t128 * t296) * t297;
t75 = t143 * t299 + (-t125 * t298 + t126 * t296) * t297;
t73 = (t327 - t350) * t299 + t298 * t311;
t72 = t296 * t311 + t299 * t349 + t329;
t70 = t127 * t280 + t128 * t282 - t144 * t356;
t69 = t125 * t280 + t126 * t282 - t143 * t356;
t68 = t280 * t332 + t333;
t57 = t280 * t321 + t185 + t333;
t56 = (t296 * t350 + t298 * t349) * t297 + t312;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(3) * t186 + m(4) * t139 + m(5) * t107 + m(6) * t84 + m(7) * t56; m(7) * (t56 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(6) * (t84 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(5) * (t107 ^ 2 + t119 ^ 2 + t120 ^ 2) + m(4) * (t139 ^ 2 + t146 ^ 2 + t147 ^ 2) + m(3) * (t186 ^ 2 + t223 ^ 2 + t224 ^ 2) + (t76 + (t244 * t361 - t246 * t282 + t248 * t283) * t361 + t378) * t361 + (-t75 + (-t243 * t360 - t245 * t280 + t247 * t281) * t360 + (-t243 * t361 + t244 * t360 + t245 * t282 + t246 * t280 - t247 * t283 - t248 * t281) * t361 - t379) * t360 + ((t274 * t361 - t282 * t275 + t283 * t276) * t361 - (-t274 * t360 - t280 * t275 + t281 * t276) * t360 + t83 + ((t246 * t305 + t248 * t302) * t296 - (t245 * t305 + t247 * t302) * t298) * t297 ^ 2 + ((-t243 * t298 + t244 * t296 + t275 * t305 + t276 * t302) * t297 + t299 * t274) * t299 + t376) * t299; m(4) * t148 + m(5) * t112 + m(6) * t85 + m(7) * t57; t82 * t368 + m(7) * (t56 * t57 + t72 * t78 + t73 * t77) + m(6) * (t102 * t99 + t103 * t98 + t84 * t85) + m(5) * (t107 * t112 + t119 * t130 + t120 * t129) + m(4) * (t139 * t148 + t146 * t179 + t147 * t178) + (t70 * t370 + t69 * t369 - t305 * t83 / 0.2e1) * t297 + t76 * t371 + t75 * t372 + t308; t280 * t69 + t282 * t70 + (-t82 - t377) * t356 + m(7) * (t57 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(6) * (t102 ^ 2 + t103 ^ 2 + t85 ^ 2) + m(5) * (t112 ^ 2 + t129 ^ 2 + t130 ^ 2) + m(4) * (t148 ^ 2 + t178 ^ 2 + t179 ^ 2) + t314; m(5) * t142 + m(6) * t106 + m(7) * t68; m(7) * (t56 * t68 + t72 * t87 + t73 * t86) + m(6) * (t106 * t84 + t113 * t99 + t114 * t98) + m(5) * (t107 * t142 + t119 * t160 + t120 * t159) + t308; m(7) * (t57 * t68 + t77 * t86 + t78 * t87) + m(6) * (t102 * t113 + t103 * t114 + t106 * t85) + m(5) * (t112 * t142 + t129 * t159 + t130 * t160) + t309; m(7) * (t68 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(6) * (t106 ^ 2 + t113 ^ 2 + t114 ^ 2) + m(5) * (t142 ^ 2 + t159 ^ 2 + t160 ^ 2) + t309; m(6) * t117 + m(7) * t81; t37 * t375 + t38 * t374 + t55 * t373 + t50 * t368 + (t25 * t369 + t26 * t370) * t297 + m(7) * (t56 * t81 + t72 * t93 + t73 * t92) + m(6) * (t117 * t84 + t135 * t99 + t136 * t98) + t315; m(7) * (t57 * t81 + t77 * t92 + t78 * t93) + m(6) * (t102 * t135 + t103 * t136 + t117 * t85) + t310; m(7) * (t68 * t81 + t86 * t92 + t87 * t93) + m(6) * (t106 * t117 + t113 * t135 + t114 * t136) + t310; t264 * t25 + t266 * t26 + t278 * t50 + m(7) * (t81 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(6) * (t117 ^ 2 + t135 ^ 2 + t136 ^ 2) + t338; m(7) * t116; m(7) * (t116 * t56 + t133 * t73 + t134 * t72) + t315; m(7) * (t116 * t57 + t133 * t77 + t134 * t78) + t316; m(7) * (t116 * t68 + t133 * t86 + t134 * t87) + t316; m(7) * (t116 * t81 + t133 * t92 + t134 * t93) + t338; m(7) * (t116 ^ 2 + t133 ^ 2 + t134 ^ 2) + t338;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
