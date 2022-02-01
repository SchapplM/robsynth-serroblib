% Calculate time derivative of joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:57
% EndTime: 2022-01-23 09:25:17
% DurationCPUTime: 11.12s
% Computational Cost: add. (14182->683), mult. (17234->954), div. (0->0), fcn. (16161->10), ass. (0->298)
t268 = qJ(3) + pkin(9);
t257 = sin(t268);
t270 = cos(pkin(8));
t275 = cos(qJ(1));
t258 = cos(t268);
t273 = sin(qJ(1));
t339 = t273 * t258;
t216 = -t257 * t275 + t270 * t339;
t340 = t273 * t257;
t284 = t258 * t275 + t270 * t340;
t269 = sin(pkin(8));
t348 = t269 * t273;
t138 = Icges(5,5) * t216 - Icges(5,6) * t284 + Icges(5,3) * t348;
t274 = cos(qJ(3));
t337 = t273 * t274;
t272 = sin(qJ(3));
t343 = t272 * t275;
t232 = t270 * t337 - t343;
t336 = t274 * t275;
t338 = t273 * t272;
t283 = t270 * t338 + t336;
t154 = Icges(4,5) * t232 - Icges(4,6) * t283 + Icges(4,3) * t348;
t380 = t138 + t154;
t345 = t270 * t275;
t217 = -t257 * t345 + t339;
t218 = t258 * t345 + t340;
t346 = t269 * t275;
t139 = Icges(5,5) * t218 + Icges(5,6) * t217 + Icges(5,3) * t346;
t233 = -t270 * t343 + t337;
t234 = t270 * t336 + t338;
t155 = Icges(4,5) * t234 + Icges(4,6) * t233 + Icges(4,3) * t346;
t373 = t139 + t155;
t323 = qJD(3) * t269;
t357 = Icges(5,4) * t257;
t199 = (-Icges(5,2) * t258 - t357) * t323;
t359 = Icges(4,4) * t272;
t222 = (-Icges(4,2) * t274 - t359) * t323;
t382 = -t257 * t199 - t272 * t222;
t170 = qJD(1) * t233 - qJD(3) * t232;
t171 = qJD(1) * t234 - qJD(3) * t283;
t324 = qJD(1) * t275;
t312 = t269 * t324;
t105 = Icges(4,5) * t171 + Icges(4,6) * t170 + Icges(4,3) * t312;
t152 = qJD(1) * t217 - qJD(3) * t216;
t153 = qJD(1) * t218 - qJD(3) * t284;
t89 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t312;
t381 = t105 + t89;
t198 = (-Icges(5,5) * t257 - Icges(5,6) * t258) * t323;
t356 = Icges(5,4) * t258;
t200 = (-Icges(5,1) * t257 - t356) * t323;
t221 = (-Icges(4,5) * t272 - Icges(4,6) * t274) * t323;
t358 = Icges(4,4) * t274;
t223 = (-Icges(4,1) * t272 - t358) * t323;
t347 = t269 * t274;
t379 = t269 * t258 * t200 + t223 * t347 + (-t198 - t221) * t270;
t140 = Icges(5,4) * t216 - Icges(5,2) * t284 + Icges(5,6) * t348;
t142 = Icges(5,1) * t216 - Icges(5,4) * t284 + Icges(5,5) * t348;
t156 = Icges(4,4) * t232 - Icges(4,2) * t283 + Icges(4,6) * t348;
t158 = Icges(4,1) * t232 - Icges(4,4) * t283 + Icges(4,5) * t348;
t378 = -t140 * t284 + t142 * t216 - t156 * t283 + t158 * t232 + t348 * t380;
t141 = Icges(5,4) * t218 + Icges(5,2) * t217 + Icges(5,6) * t346;
t143 = Icges(5,1) * t218 + Icges(5,4) * t217 + Icges(5,5) * t346;
t157 = Icges(4,4) * t234 + Icges(4,2) * t233 + Icges(4,6) * t346;
t159 = Icges(4,1) * t234 + Icges(4,4) * t233 + Icges(4,5) * t346;
t377 = -t217 * t141 - t218 * t143 - t233 * t157 - t234 * t159 - t346 * t373;
t376 = t141 * t284 - t143 * t216 + t157 * t283 - t159 * t232 - t348 * t373;
t375 = t217 * t140 + t218 * t142 + t233 * t156 + t234 * t158 + t346 * t380;
t210 = -Icges(4,6) * t270 + (-Icges(4,2) * t272 + t358) * t269;
t211 = -Icges(4,5) * t270 + (Icges(4,1) * t274 - t359) * t269;
t191 = -Icges(5,6) * t270 + (-Icges(5,2) * t257 + t356) * t269;
t192 = -Icges(5,5) * t270 + (Icges(5,1) * t258 - t357) * t269;
t370 = -t191 * t258 - t192 * t257;
t374 = t379 + ((-t210 * t274 - t211 * t272 + t370) * qJD(3) + t382) * t269;
t168 = qJD(1) * t283 - qJD(3) * t234;
t169 = -qJD(1) * t232 + qJD(3) * t233;
t325 = qJD(1) * t273;
t313 = t269 * t325;
t104 = Icges(4,5) * t169 + Icges(4,6) * t168 - Icges(4,3) * t313;
t150 = qJD(1) * t284 - qJD(3) * t218;
t151 = -qJD(1) * t216 + qJD(3) * t217;
t88 = Icges(5,5) * t151 + Icges(5,6) * t150 - Icges(5,3) * t313;
t372 = (t104 + t88) * t275;
t364 = pkin(3) * t272;
t241 = pkin(4) * t257 + t364;
t371 = t273 * (qJ(2) + t241);
t259 = qJ(5) + t268;
t253 = sin(t259);
t254 = cos(t259);
t341 = t273 * t254;
t205 = -t253 * t345 + t341;
t342 = t273 * t253;
t206 = t254 * t345 + t342;
t130 = rSges(6,1) * t206 + rSges(6,2) * t205 + rSges(6,3) * t346;
t363 = pkin(3) * t274;
t297 = pkin(4) * t258 + t363;
t240 = pkin(2) + t297;
t271 = qJ(4) + pkin(6);
t266 = -pkin(7) - t271;
t262 = t273 * qJ(2);
t326 = pkin(1) * t275 + t262;
t87 = t240 * t345 + t241 * t273 - t266 * t346 + t130 + t326;
t369 = 2 * m(4);
t368 = 2 * m(5);
t367 = 2 * m(6);
t366 = m(5) / 0.2e1;
t365 = m(6) / 0.2e1;
t362 = pkin(2) * t270 + pkin(1);
t267 = qJD(3) + qJD(5);
t349 = t267 * t269;
t355 = Icges(6,4) * t253;
t178 = (-Icges(6,2) * t254 - t355) * t349;
t185 = -Icges(6,5) * t270 + (Icges(6,1) * t254 - t355) * t269;
t177 = (-Icges(6,5) * t253 - Icges(6,6) * t254) * t349;
t354 = Icges(6,4) * t254;
t179 = (-Icges(6,1) * t253 - t354) * t349;
t305 = t179 * t254 * t269 - t270 * t177;
t184 = -Icges(6,6) * t270 + (-Icges(6,2) * t253 + t354) * t269;
t317 = t267 * t254 * t184;
t48 = (-t317 + (-t185 * t267 - t178) * t253) * t269 + t305;
t361 = t48 * t270;
t322 = qJD(3) * t272;
t318 = pkin(3) * t322;
t320 = qJD(4) * t269;
t236 = -t270 * t318 + t320;
t321 = qJD(3) * t274;
t248 = pkin(3) * t321 + qJD(2);
t255 = qJ(2) + t364;
t314 = t236 * t275 + t248 * t273 + t255 * t324;
t327 = qJ(2) * t324 + qJD(2) * t273;
t219 = t269 * t271 + t270 * t363 + t362;
t246 = pkin(6) * t269 + t362;
t328 = t219 - t246;
t102 = -t325 * t328 + t314 - t327;
t331 = rSges(5,1) * t151 + rSges(5,2) * t150;
t360 = rSges(5,3) * t313 - t102 - t331;
t351 = t255 * t275;
t285 = t254 * t275 + t270 * t342;
t134 = qJD(1) * t285 - t206 * t267;
t204 = -t253 * t275 + t270 * t341;
t135 = -qJD(1) * t204 + t205 * t267;
t334 = rSges(6,1) * t135 + rSges(6,2) * t134;
t263 = t275 * qJ(2);
t144 = t273 * t328 + t263 - t351;
t220 = pkin(3) * t347 + (pkin(6) - t271) * t270;
t333 = t144 * t270 + t220 * t348;
t301 = -t246 * t275 - t262;
t329 = t219 * t275 + t255 * t273;
t145 = t301 + t329;
t147 = rSges(5,1) * t218 + rSges(5,2) * t217 + rSges(5,3) * t346;
t332 = -t145 - t147;
t330 = rSges(4,1) * t169 + rSges(4,2) * t168;
t290 = -rSges(6,1) * t204 + rSges(6,2) * t285;
t129 = rSges(6,3) * t348 - t290;
t186 = -rSges(6,3) * t270 + (rSges(6,1) * t254 - rSges(6,2) * t253) * t269;
t99 = t129 * t270 + t186 * t348;
t238 = t241 * qJD(3);
t239 = t297 * qJD(3);
t280 = -t238 * t345 + t239 * t273 + t241 * t324 + t266 * t313 + t275 * t320 + t327;
t308 = -t240 * t270 - pkin(1);
t298 = t219 + t308;
t82 = -rSges(6,3) * t313 + t334;
t319 = -t298 * t325 - t102 - t280 + t314 - t82;
t180 = (-rSges(6,1) * t253 - rSges(6,2) * t254) * t349;
t136 = qJD(1) * t205 - t204 * t267;
t137 = qJD(1) * t206 - t267 * t285;
t291 = -rSges(6,1) * t137 - rSges(6,2) * t136;
t83 = rSges(6,3) * t312 - t291;
t45 = t180 * t348 + t186 * t312 + t270 * t83;
t261 = qJD(2) * t275;
t245 = t255 * t325;
t302 = -t236 * t273 - t245;
t103 = -t248 * t275 + t261 + (t275 * t328 - t262) * qJD(1) - t302;
t311 = t269 * t322;
t237 = -pkin(3) * t311 - qJD(4) * t270;
t316 = t103 * t270 + t220 * t312 + t237 * t348;
t315 = t329 - t145 - t87;
t161 = rSges(4,1) * t234 + rSges(4,2) * t233 + rSges(4,3) * t346;
t307 = -rSges(4,3) * t269 - t246;
t306 = -rSges(5,3) * t269 - t219;
t296 = rSges(3,1) * t270 - rSges(3,2) * t269;
t295 = -rSges(4,1) * t171 - rSges(4,2) * t170;
t294 = -rSges(4,1) * t232 + rSges(4,2) * t283;
t293 = -rSges(5,1) * t153 - rSges(5,2) * t152;
t292 = -rSges(5,1) * t216 + rSges(5,2) * t284;
t289 = t307 * t273;
t288 = t306 * qJD(1);
t287 = t238 * t270 - t320;
t286 = -pkin(1) - t296;
t125 = Icges(6,4) * t204 - Icges(6,2) * t285 + Icges(6,6) * t348;
t127 = Icges(6,1) * t204 - Icges(6,4) * t285 + Icges(6,5) * t348;
t77 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t312;
t79 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t312;
t81 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t312;
t12 = -t270 * t77 + ((-t125 * t267 + t81) * t254 + (-t127 * t267 - t79) * t253) * t269;
t126 = Icges(6,4) * t206 + Icges(6,2) * t205 + Icges(6,6) * t346;
t128 = Icges(6,1) * t206 + Icges(6,4) * t205 + Icges(6,5) * t346;
t76 = Icges(6,5) * t135 + Icges(6,6) * t134 - Icges(6,3) * t313;
t78 = Icges(6,4) * t135 + Icges(6,2) * t134 - Icges(6,6) * t313;
t80 = Icges(6,1) * t135 + Icges(6,4) * t134 - Icges(6,5) * t313;
t13 = -t270 * t76 + ((-t126 * t267 + t80) * t254 + (-t128 * t267 - t78) * t253) * t269;
t183 = -Icges(6,3) * t270 + (Icges(6,5) * t254 - Icges(6,6) * t253) * t269;
t22 = t134 * t184 + t135 * t185 + t205 * t178 + t206 * t179 + (t177 * t275 - t183 * t325) * t269;
t23 = t136 * t184 + t137 * t185 - t285 * t178 + t204 * t179 + (t177 * t273 + t183 * t324) * t269;
t123 = Icges(6,5) * t204 - Icges(6,6) * t285 + Icges(6,3) * t348;
t46 = -t123 * t270 + (-t125 * t253 + t127 * t254) * t269;
t124 = Icges(6,5) * t206 + Icges(6,6) * t205 + Icges(6,3) * t346;
t47 = -t124 * t270 + (-t126 * t253 + t128 * t254) * t269;
t68 = t183 * t348 - t184 * t285 + t185 * t204;
t69 = t183 * t346 + t184 * t205 + t185 * t206;
t282 = (t12 + t23) * t348 / 0.2e1 - (t47 + t69) * t313 / 0.2e1 + (t13 + t22 + (t46 + t68) * qJD(1)) * t346 / 0.2e1;
t281 = -t266 * t269 - t298;
t279 = (-rSges(6,3) + t266) * t269 + t308;
t34 = t123 * t348 - t125 * t285 + t127 * t204;
t35 = t124 * t348 - t126 * t285 + t128 * t204;
t36 = t123 * t346 + t125 * t205 + t127 * t206;
t37 = t124 * t346 + t126 * t205 + t128 * t206;
t277 = -t270 * (-t361 + (t12 * t273 + t13 * t275 + (-t273 * t47 + t275 * t46) * qJD(1)) * t269) - t313 * (-t69 * t270 + (t273 * t36 + t275 * t37) * t269) + (-t23 * t270 + ((t136 * t126 + t137 * t128 + t204 * t80 - t285 * t78) * t275 - t35 * t325 + (t136 * t125 + t137 * t127 + t204 * t81 - t285 * t79) * t273 + t34 * t324 + ((t124 * t324 + t273 * t76) * t275 + (t123 * t324 + t273 * t77) * t273) * t269) * t269) * t348 + (-t22 * t270 + ((t134 * t126 + t135 * t128 + t205 * t78 + t206 * t80) * t275 - t37 * t325 + (t134 * t125 + t135 * t127 + t205 * t79 + t206 * t81) * t273 + t36 * t324 + ((-t124 * t325 + t275 * t76) * t275 + (-t123 * t325 + t275 * t77) * t273) * t269) * t269) * t346 + (-t68 * t270 + (t273 * t34 + t275 * t35) * t269) * t312;
t276 = rSges(3,3) * t275 + t273 * t286;
t225 = (-rSges(4,1) * t272 - rSges(4,2) * t274) * t323;
t214 = -rSges(4,3) * t270 + (rSges(4,1) * t274 - rSges(4,2) * t272) * t269;
t209 = -Icges(4,3) * t270 + (Icges(4,5) * t274 - Icges(4,6) * t272) * t269;
t202 = (-rSges(5,1) * t257 - rSges(5,2) * t258) * t323;
t195 = -rSges(5,3) * t270 + (rSges(5,1) * t258 - rSges(5,2) * t257) * t269;
t190 = -Icges(5,3) * t270 + (Icges(5,5) * t258 - Icges(5,6) * t257) * t269;
t188 = t220 * t313;
t187 = (-t238 + t318) * t269;
t182 = rSges(3,3) * t273 + t275 * t296 + t326;
t181 = t263 + t276;
t173 = t186 * t313;
t164 = t261 + ((-rSges(3,3) - qJ(2)) * t273 + t286 * t275) * qJD(1);
t163 = qJD(1) * t276 + t327;
t162 = (t266 + t271) * t270 + (-pkin(2) + t240 - t363) * t269;
t160 = rSges(4,3) * t348 - t294;
t146 = rSges(5,3) * t348 - t292;
t122 = t144 * t346;
t120 = t263 + t289 + t294;
t119 = -t301 + t161;
t118 = t129 * t346;
t117 = -t161 * t270 - t214 * t346;
t116 = t160 * t270 + t214 * t348;
t115 = t147 + t329;
t114 = t273 * t306 + t292 + t351;
t113 = rSges(4,3) * t312 - t295;
t112 = -rSges(4,3) * t313 + t330;
t110 = -t263 + (-t241 + t255) * t275 + t281 * t273;
t109 = Icges(4,1) * t171 + Icges(4,4) * t170 + Icges(4,5) * t312;
t108 = Icges(4,1) * t169 + Icges(4,4) * t168 - Icges(4,5) * t313;
t107 = Icges(4,4) * t171 + Icges(4,2) * t170 + Icges(4,6) * t312;
t106 = Icges(4,4) * t169 + Icges(4,2) * t168 - Icges(4,6) * t313;
t100 = -t130 * t270 - t186 * t346;
t98 = t103 * t346;
t97 = t209 * t346 + t210 * t233 + t211 * t234;
t96 = t209 * t348 - t210 * t283 + t211 * t232;
t95 = rSges(5,3) * t312 - t293;
t93 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t312;
t92 = Icges(5,1) * t151 + Icges(5,4) * t150 - Icges(5,5) * t313;
t91 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t312;
t90 = Icges(5,4) * t151 + Icges(5,2) * t150 - Icges(5,6) * t313;
t86 = t241 * t275 + t273 * t279 + t263 + t290;
t85 = qJD(1) * t289 + t327 + t330;
t84 = t261 + (t275 * t307 - t262) * qJD(1) + t295;
t75 = t190 * t346 + t191 * t217 + t192 * t218;
t74 = t190 * t348 - t191 * t284 + t192 * t216;
t72 = -t130 * t348 + t118;
t71 = t83 * t346;
t67 = t270 * t113 + (t214 * t324 + t225 * t273) * t269;
t66 = -t270 * t112 + (t214 * t325 - t225 * t275) * t269;
t63 = -t155 * t270 + (-t157 * t272 + t159 * t274) * t269;
t62 = -t154 * t270 + (-t156 * t272 + t158 * t274) * t269;
t60 = t273 * t288 + t314 + t331;
t59 = (t248 + t288) * t275 + t293 + t302;
t58 = t332 * t270 + (-t195 - t220) * t346;
t57 = t146 * t270 + t195 * t348 + t333;
t51 = -t245 - t261 + (-t239 + t248) * t275 + (-t236 - t287) * t273 + (t275 * t281 + t371) * qJD(1);
t50 = -t139 * t270 + (-t141 * t257 + t143 * t258) * t269;
t49 = -t138 * t270 + (-t140 * t257 + t142 * t258) * t269;
t44 = -t180 * t346 - t270 * t82 + t173;
t33 = t239 * t275 + t261 + t287 * t273 + (t275 * t279 - t371) * qJD(1) + t291;
t32 = (-rSges(6,3) * t269 + t308) * t325 + t280 + t334;
t31 = t170 * t210 + t171 * t211 - t283 * t222 + t232 * t223 + (t209 * t324 + t221 * t273) * t269;
t30 = t168 * t210 + t169 * t211 + t233 * t222 + t234 * t223 + (-t209 * t325 + t221 * t275) * t269;
t29 = t152 * t191 + t153 * t192 - t284 * t199 + t216 * t200 + (t190 * t324 + t198 * t273) * t269;
t28 = t150 * t191 + t151 * t192 + t217 * t199 + t218 * t200 + (-t190 * t325 + t198 * t275) * t269;
t27 = t270 * t95 + (t195 * t324 + t202 * t273) * t269 + t316;
t26 = t188 + t360 * t270 + (t195 * t325 + (-t202 - t237) * t275) * t269;
t25 = t315 * t270 + (-t162 - t186 - t220) * t346;
t24 = t110 * t270 + t162 * t348 + t333 + t99;
t19 = -t104 * t270 + (-t106 * t272 + t108 * t274 + (-t157 * t274 - t159 * t272) * qJD(3)) * t269;
t18 = -t105 * t270 + (-t107 * t272 + t109 * t274 + (-t156 * t274 - t158 * t272) * qJD(3)) * t269;
t17 = t71 + (-t273 * t82 + (-t129 * t273 - t130 * t275) * qJD(1)) * t269;
t16 = t118 + t122 + (t110 * t275 + t273 * t315) * t269;
t15 = -t270 * t88 + (-t257 * t90 + t258 * t92 + (-t141 * t258 - t143 * t257) * qJD(3)) * t269;
t14 = -t270 * t89 + (-t257 * t91 + t258 * t93 + (-t140 * t258 - t142 * t257) * qJD(3)) * t269;
t9 = t270 * t51 + (t162 * t324 + t187 * t273) * t269 + t316 + t45;
t8 = t173 + t188 + t319 * t270 + (t162 * t325 + (-t180 - t187 - t237) * t275) * t269;
t5 = t98 + (t275 * t95 + t360 * t273 + (t332 * t275 + (-t144 - t146) * t273) * qJD(1)) * t269;
t4 = t71 + t98 + (t275 * t51 + t319 * t273 + (t315 * t275 + (-t110 - t129 - t144) * t273) * qJD(1)) * t269;
t1 = [(t32 * t87 + t33 * t86) * t367 - t253 * t185 * t349 - t211 * t311 + (t114 * t59 + t115 * t60) * t368 + (t119 * t85 + t120 * t84) * t369 + 0.2e1 * m(3) * (t163 * t182 + t164 * t181) + t305 + t370 * t323 + (-t178 * t253 - t210 * t321 - t317 + t382) * t269 + t379; m(6) * (t273 * t33 - t275 * t32 + (t273 * t87 + t275 * t86) * qJD(1)) + m(5) * (t273 * t59 - t275 * t60 + (t114 * t275 + t115 * t273) * qJD(1)) + m(4) * (t273 * t84 - t275 * t85 + (t119 * t273 + t120 * t275) * qJD(1)) + m(3) * (-t163 * t275 + t273 * t164 + (t181 * t275 + t182 * t273) * qJD(1)); 0; (-t48 - t374) * t270 + m(5) * (t114 * t27 + t115 * t26 + t57 * t59 + t58 * t60) + m(4) * (t116 * t84 + t117 * t85 + t119 * t66 + t120 * t67) + m(6) * (t24 * t33 + t25 * t32 + t8 * t87 + t86 * t9) + ((t19 / 0.2e1 + t15 / 0.2e1 + t28 / 0.2e1 + t30 / 0.2e1) * t275 + (t18 / 0.2e1 + t14 / 0.2e1 + t29 / 0.2e1 + t31 / 0.2e1) * t273 + ((t62 / 0.2e1 + t49 / 0.2e1 + t74 / 0.2e1 + t96 / 0.2e1) * t275 + (-t75 / 0.2e1 - t97 / 0.2e1 - t63 / 0.2e1 - t50 / 0.2e1) * t273) * qJD(1)) * t269 + t282; m(4) * (t67 * t273 - t275 * t66 + (t116 * t275 + t117 * t273) * qJD(1)) + m(5) * (-t26 * t275 + t27 * t273 + (t273 * t58 + t275 * t57) * qJD(1)) + m(6) * (t9 * t273 - t275 * t8 + (t24 * t275 + t25 * t273) * qJD(1)); (t16 * t4 + t24 * t9 + t25 * t8) * t367 + (t122 * t5 + t58 * t26 + t57 * t27) * t368 + (t116 * t67 + t117 * t66) * t369 + t277 + ((t146 * t275 + t273 * t332) * t5 * t368 + (t160 * t275 - t161 * t273) * (-t112 * t273 + t113 * t275 + (-t160 * t273 - t161 * t275) * qJD(1)) * t369 * t269 + (-t273 * t375 + t275 * t377) * t313 + (t273 * t378 - t275 * t376) * t312 + (t376 * t325 + t378 * t324 + (-t106 * t283 + t232 * t108 + t152 * t141 + t153 * t143 + t170 * t157 + t171 * t159 + t216 * t92 - t284 * t90 + t312 * t373) * t275 + (-t283 * t107 + t232 * t109 + t152 * t140 + t153 * t142 + t170 * t156 + t171 * t158 - t284 * t91 + t216 * t93 + (t273 * t381 + t324 * t380 + t372) * t269) * t273) * t348 + (t377 * t325 + t375 * t324 + (t233 * t106 + t234 * t108 + t150 * t141 + t151 * t143 + t168 * t157 + t169 * t159 + t217 * t90 + t218 * t92 + (-t325 * t373 + t372) * t269) * t275 + (t233 * t107 + t234 * t109 + t168 * t156 + t169 * t158 + t150 * t140 + t151 * t142 + t217 * t91 + t218 * t93 + (t275 * t381 - t325 * t380) * t269) * t273) * t346) * t269 + (t374 * t270 + (-t31 - t29) * t348 + (-t30 - t28) * t346 + (t97 + t75) * t313 + (-t96 - t74) * t312 + ((-t15 - t19) * t275 + (-t14 - t18) * t273 + ((-t49 - t62) * t275 + (t50 + t63) * t273) * qJD(1)) * t269) * t270; 0.2e1 * ((t273 * t32 + t275 * t33 + t324 * t87 - t325 * t86) * t365 + (-t114 * t325 + t115 * t324 + t273 * t60 + t275 * t59) * t366) * t269; 0; 0.2e1 * (-m(6) * t4 / 0.2e1 - m(5) * t5 / 0.2e1) * t270 + 0.2e1 * ((-t24 * t325 + t25 * t324 + t273 * t8 + t275 * t9) * t365 + (t26 * t273 + t27 * t275 + t324 * t58 - t325 * t57) * t366) * t269; 0; m(6) * (t100 * t32 + t33 * t99 + t44 * t87 + t45 * t86) - t361 + t282; m(6) * (t45 * t273 - t275 * t44 + (t100 * t273 + t275 * t99) * qJD(1)); m(6) * (t100 * t8 + t16 * t17 + t24 * t45 + t25 * t44 + t4 * t72 + t9 * t99) + t277; m(6) * (-t17 * t270 + (t273 * t44 + t275 * t45 + (t100 * t275 - t273 * t99) * qJD(1)) * t269); (t100 * t44 + t17 * t72 + t45 * t99) * t367 + t277;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
