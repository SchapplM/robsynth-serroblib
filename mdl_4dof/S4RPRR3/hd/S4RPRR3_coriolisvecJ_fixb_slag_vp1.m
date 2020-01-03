% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:15
% DurationCPUTime: 8.02s
% Computational Cost: add. (9272->484), mult. (8640->664), div. (0->0), fcn. (6762->8), ass. (0->291)
t223 = qJD(1) ^ 2;
t216 = qJ(3) + qJ(4);
t211 = cos(t216);
t206 = Icges(5,4) * t211;
t210 = sin(t216);
t153 = Icges(5,1) * t210 + t206;
t259 = -Icges(5,2) * t210 + t206;
t411 = t153 + t259;
t215 = qJ(1) + pkin(7);
t208 = sin(t215);
t346 = t208 * t210;
t165 = Icges(5,4) * t346;
t345 = t208 * t211;
t209 = cos(t215);
t362 = Icges(5,5) * t209;
t103 = Icges(5,1) * t345 - t165 - t362;
t358 = Icges(5,6) * t209;
t101 = Icges(5,4) * t345 - Icges(5,2) * t346 - t358;
t355 = t101 * t210;
t258 = -t103 * t211 + t355;
t245 = t258 * t208;
t150 = Icges(5,5) * t211 - Icges(5,6) * t210;
t100 = Icges(5,3) * t208 + t150 * t209;
t364 = Icges(5,4) * t210;
t154 = Icges(5,1) * t211 - t364;
t104 = Icges(5,5) * t208 + t154 * t209;
t340 = t209 * t211;
t380 = t208 * t100 + t104 * t340;
t410 = -t245 - t380;
t409 = 2 * qJD(3);
t168 = rSges(5,2) * t346;
t105 = rSges(5,1) * t345 - t209 * rSges(5,3) - t168;
t218 = sin(qJ(1));
t386 = pkin(1) * t218;
t408 = -t105 - t386;
t199 = t208 * rSges(4,3);
t219 = cos(qJ(3));
t338 = t209 * t219;
t217 = sin(qJ(3));
t339 = t209 * t217;
t115 = rSges(4,1) * t338 - rSges(4,2) * t339 + t199;
t148 = t209 * pkin(2) + t208 * pkin(5);
t220 = cos(qJ(1));
t213 = t220 * pkin(1);
t290 = t148 + t213;
t407 = t115 + t290;
t287 = t209 * rSges(3,1) - rSges(3,2) * t208;
t212 = Icges(4,4) * t219;
t260 = -Icges(4,2) * t217 + t212;
t177 = Icges(4,1) * t217 + t212;
t406 = t213 + t287;
t214 = qJD(3) + qJD(4);
t373 = rSges(5,2) * t211;
t305 = t214 * t373;
t337 = t210 * t214;
t405 = rSges(5,1) * t337 + t305;
t174 = Icges(4,5) * t219 - Icges(4,6) * t217;
t173 = Icges(4,5) * t217 + Icges(4,6) * t219;
t242 = qJD(3) * t173;
t365 = Icges(4,4) * t217;
t178 = Icges(4,1) * t219 - t365;
t113 = Icges(4,5) * t208 + t178 * t209;
t111 = Icges(4,6) * t208 + t209 * t260;
t352 = t111 * t217;
t255 = -t113 * t219 + t352;
t357 = Icges(4,3) * t209;
t403 = -t209 * t242 + (-t174 * t208 + t255 + t357) * qJD(1);
t344 = t208 * t217;
t185 = Icges(4,4) * t344;
t343 = t208 * t219;
t363 = Icges(4,5) * t209;
t112 = Icges(4,1) * t343 - t185 - t363;
t359 = Icges(4,6) * t209;
t110 = Icges(4,4) * t343 - Icges(4,2) * t344 - t359;
t353 = t110 * t217;
t256 = -t112 * t219 + t353;
t109 = Icges(4,3) * t208 + t174 * t209;
t319 = qJD(1) * t109;
t402 = qJD(1) * t256 - t208 * t242 + t319;
t149 = Icges(5,5) * t210 + Icges(5,6) * t211;
t246 = t149 * t214;
t102 = Icges(5,6) * t208 + t209 * t259;
t354 = t102 * t210;
t356 = Icges(5,3) * t209;
t401 = -t209 * t246 + (-t104 * t211 - t150 * t208 + t354 + t356) * qJD(1);
t320 = qJD(1) * t100;
t400 = qJD(1) * t258 - t208 * t246 + t320;
t108 = Icges(4,5) * t343 - Icges(4,6) * t344 - t357;
t46 = -t209 * t108 - t208 * t256;
t151 = Icges(5,2) * t211 + t364;
t253 = t151 * t210 - t153 * t211;
t399 = qJD(1) * t253 + t150 * t214;
t175 = Icges(4,2) * t219 + t365;
t252 = t217 * t175 - t219 * t177;
t398 = t252 * qJD(1) + t174 * qJD(3);
t198 = t208 * rSges(5,3);
t325 = rSges(5,1) * t340 + t198;
t341 = t209 * t210;
t106 = -rSges(5,2) * t341 + t325;
t155 = rSges(5,1) * t210 + t373;
t123 = t155 * t208;
t124 = t155 * t209;
t144 = t208 * t214;
t145 = t209 * t214;
t374 = rSges(5,2) * t210;
t376 = rSges(5,1) * t211;
t156 = -t374 + t376;
t204 = t209 * pkin(5);
t147 = pkin(2) * t208 - t204;
t221 = -pkin(6) - pkin(5);
t190 = t209 * t221;
t385 = pkin(3) * t219;
t207 = pkin(2) + t385;
t326 = -t208 * t207 - t190;
t95 = t147 + t326;
t277 = t209 * t207 - t208 * t221;
t96 = t277 - t148;
t35 = t105 * t144 + t106 * t145 + qJD(2) + (-t208 * t95 + t209 * t96) * qJD(3);
t312 = qJD(3) * t217;
t299 = t209 * t312;
t275 = pkin(3) * t299;
t241 = -t145 * t155 - t275;
t273 = t95 + t408;
t38 = (-t147 + t273) * qJD(1) + t241;
t304 = pkin(3) * t312;
t171 = t208 * t304;
t368 = -t106 - t96;
t39 = -t144 * t155 - t171 + (t290 - t368) * qJD(1);
t274 = qJD(1) * t214;
t140 = t208 * t274;
t141 = t209 * t274;
t317 = qJD(1) * t208;
t316 = qJD(1) * t209;
t327 = rSges(5,3) * t316 + qJD(1) * t168;
t66 = -t209 * t305 + (-t209 * t337 - t211 * t317) * rSges(5,1) + t327;
t67 = -t214 * t123 + (t156 * t209 + t198) * qJD(1);
t191 = pkin(5) * t316;
t382 = pkin(2) - t207;
t81 = -t275 - t191 + (t208 * t382 - t190) * qJD(1);
t82 = -t171 + (-t382 * t209 + (-pkin(5) - t221) * t208) * qJD(1);
t8 = t105 * t141 - t106 * t140 + t144 * t67 + t145 * t66 + (t208 * t82 + t209 * t81 + (-t208 * t96 - t209 * t95) * qJD(1)) * qJD(3);
t397 = -t38 * (qJD(1) * t123 - t145 * t156) - t35 * (-t144 * t123 - t124 * t145) - t39 * (-qJD(1) * t124 - t144 * t156) + t8 * (t208 * t105 + t209 * t106);
t330 = -Icges(4,2) * t343 + t112 - t185;
t332 = t177 * t208 + t110;
t396 = -t217 * t330 - t219 * t332;
t395 = qJD(1) * t411 + t144 * (-t151 * t209 + t104) - t145 * (-Icges(5,2) * t345 + t103 - t165);
t394 = t140 / 0.2e1;
t393 = t141 / 0.2e1;
t392 = -t144 / 0.2e1;
t391 = t144 / 0.2e1;
t390 = -t145 / 0.2e1;
t389 = t145 / 0.2e1;
t388 = t208 / 0.2e1;
t387 = -t209 / 0.2e1;
t384 = -qJD(1) / 0.2e1;
t383 = qJD(1) / 0.2e1;
t99 = Icges(5,5) * t345 - Icges(5,6) * t346 - t356;
t381 = -t103 * t340 - t208 * t99;
t379 = -t208 * t108 - t112 * t338;
t378 = t208 * t109 + t113 * t338;
t377 = rSges(4,1) * t219;
t179 = rSges(4,1) * t217 + rSges(4,2) * t219;
t138 = t179 * t209;
t314 = qJD(3) * t208;
t301 = t179 * t314;
t59 = qJD(1) * t407 - t301;
t372 = t138 * t59;
t321 = rSges(4,2) * t344 + t209 * rSges(4,3);
t114 = rSges(4,1) * t343 - t321;
t291 = -t114 - t386;
t313 = qJD(3) * t209;
t300 = t179 * t313;
t58 = -t300 + (-t147 + t291) * qJD(1);
t371 = t208 * t58;
t370 = t209 * t39;
t369 = t209 * t58;
t351 = t149 * t208;
t350 = t149 * t209;
t349 = t151 * t214;
t348 = t173 * t208;
t347 = t173 * t209;
t56 = -t208 * t253 - t350;
t336 = t56 * qJD(1);
t77 = -t208 * t252 - t347;
t335 = t77 * qJD(1);
t331 = -t177 * t209 - t111;
t329 = -t175 * t209 + t113;
t315 = qJD(1) * t217;
t324 = rSges(4,2) * t208 * t315 + rSges(4,3) * t316;
t323 = -t175 + t178;
t322 = t177 + t260;
t318 = qJD(1) * t174;
t311 = qJD(3) * t219;
t310 = t105 * t316 + t208 * t67 + t209 * t66;
t309 = (qJD(3) ^ 2) * t385;
t308 = t223 * t386;
t307 = t223 * t213;
t298 = -pkin(2) - t377;
t297 = t317 / 0.2e1;
t296 = t316 / 0.2e1;
t295 = -t314 / 0.2e1;
t292 = t313 / 0.2e1;
t289 = -pkin(3) * t217 - t155;
t248 = t153 * t214;
t286 = qJD(1) * t104 - t101 * t214 - t208 * t248;
t285 = -t99 + t354;
t284 = -t102 * t214 - t209 * t248 + (-t154 * t208 + t362) * qJD(1);
t283 = qJD(1) * t102 + t103 * t214 - t208 * t349;
t282 = t104 * t214 - t209 * t349 + (-t208 * t259 + t358) * qJD(1);
t92 = t113 * t343;
t281 = t209 * t109 - t92;
t280 = -t108 + t352;
t279 = t411 * t214;
t278 = t154 * t214 - t349;
t272 = qJD(1) * (-pkin(2) * t317 + t191) - t308;
t130 = t156 * t214;
t269 = -pkin(3) * t311 - t130;
t83 = t104 * t345;
t268 = t102 * t346 - t83;
t146 = rSges(3,1) * t208 + rSges(3,2) * t209;
t265 = -rSges(4,2) * t217 + t377;
t264 = -t208 * t39 - t209 * t38;
t263 = -t208 * t59 - t369;
t50 = t101 * t211 + t103 * t210;
t68 = t219 * t110 + t112 * t217;
t69 = t111 * t219 + t113 * t217;
t254 = t114 * t208 + t115 * t209;
t137 = t179 * t208;
t47 = -t111 * t344 - t281;
t250 = (t208 * t47 - t209 * t46) * qJD(3);
t48 = -t110 * t339 - t379;
t49 = -t111 * t339 + t378;
t249 = (t208 * t49 - t209 * t48) * qJD(3);
t244 = qJD(3) * t177;
t243 = qJD(3) * t175;
t239 = qJD(1) * t150 - t144 * t350 + t145 * t351;
t238 = -t217 * t329 + t219 * t331;
t229 = -t210 * t282 + t211 * t284 + t320;
t10 = t208 * t401 + t229 * t209;
t230 = qJD(1) * t99 - t210 * t283 + t211 * t286;
t11 = t230 * t208 - t209 * t400;
t12 = t229 * t208 - t209 * t401;
t40 = -t209 * t99 - t245;
t41 = -t209 * t100 - t268;
t20 = t144 * t41 - t145 * t40 + t336;
t42 = -t101 * t341 - t381;
t43 = -t102 * t341 + t380;
t57 = -t209 * t253 + t351;
t52 = t57 * qJD(1);
t21 = t144 * t43 - t145 * t42 + t52;
t233 = (-t153 * t209 - t102) * t144 - (-t153 * t208 - t101) * t145 + (-t151 + t154) * qJD(1);
t224 = -t210 * t395 + t233 * t211;
t28 = t210 * t286 + t211 * t283;
t29 = t210 * t284 + t211 * t282;
t228 = qJD(1) * t149 - t210 * t279 + t211 * t278;
t30 = t208 * t399 + t228 * t209;
t31 = t228 * t208 - t209 * t399;
t51 = t102 * t211 + t104 * t210;
t9 = t208 * t400 + t230 * t209;
t237 = (qJD(1) * t30 + t10 * t144 + t140 * t42 + t141 * t43 - t145 * t9) * t388 + (t233 * t210 + t211 * t395) * t384 + t20 * t297 + t21 * t296 + (qJD(1) * t31 - t11 * t145 + t12 * t144 + t140 * t40 + t141 * t41) * t387 + (t208 * t41 - t209 * t40) * t394 + (t208 * t43 - t209 * t42) * t393 + (t10 * t208 - t209 * t9 + (t208 * t42 + t209 * t43) * qJD(1)) * t391 + (-t11 * t209 + t12 * t208 + (t208 * t40 + t209 * t41) * qJD(1)) * t390 + (t208 * t29 - t209 * t28 + (t208 * t50 + t209 * t51) * qJD(1)) * t383 + (t208 * t239 + t209 * t224) * t392 + (t208 * t224 - t209 * t239) * t389;
t236 = (-t217 * t322 + t219 * t323) * qJD(1);
t79 = -rSges(4,2) * t209 * t311 + (-t219 * t317 - t299) * rSges(4,1) + t324;
t80 = -qJD(3) * t137 + (t209 * t265 + t199) * qJD(1);
t234 = t208 * t80 + t209 * t79 + (t114 * t209 - t115 * t208) * qJD(1);
t74 = qJD(1) * t111 - t208 * t243;
t76 = qJD(1) * t113 - t208 * t244;
t227 = qJD(1) * t108 - qJD(3) * t68 - t217 * t74 + t219 * t76;
t73 = -t209 * t243 + (-t208 * t260 + t359) * qJD(1);
t75 = -t209 * t244 + (-t178 * t208 + t363) * qJD(1);
t226 = -qJD(3) * t69 - t217 * t73 + t219 * t75 + t319;
t158 = t260 * qJD(3);
t159 = t178 * qJD(3);
t225 = qJD(1) * t173 - t158 * t217 + t159 * t219 + (-t175 * t219 - t177 * t217) * qJD(3);
t161 = t265 * qJD(3);
t143 = qJD(1) * t147;
t142 = t148 * qJD(1);
t78 = -t209 * t252 + t348;
t70 = t78 * qJD(1);
t55 = qJD(3) * t254 + qJD(2);
t45 = -t307 - t161 * t313 + (-t142 - t80 + t301) * qJD(1);
t44 = -t161 * t314 + (t79 - t300) * qJD(1) + t272;
t37 = t225 * t208 - t209 * t398;
t36 = t208 * t398 + t225 * t209;
t34 = -qJD(3) * t255 + t217 * t75 + t219 * t73;
t33 = -t256 * qJD(3) + t217 * t76 + t219 * t74;
t32 = t234 * qJD(3);
t27 = -t209 * t309 - t307 - t130 * t145 + t140 * t155 + (-t142 - t67 - t82 + t171) * qJD(1);
t26 = -t208 * t309 - t130 * t144 - t141 * t155 + (t66 + t81 - t275) * qJD(1) + t272;
t23 = t70 + t249;
t22 = t250 + t335;
t1 = [m(3) * ((-t146 * t223 - t308) * t406 + (-t307 + (-0.2e1 * t287 - t213 + t406) * t223) * (-t146 - t386)) + (t52 + (t41 + (t100 + t355) * t209 + t268 + t381) * t145 + (-t209 * t285 + t40 - t410) * t144) * t389 + (t70 + ((t47 - t92 + (t109 + t353) * t209 + t379) * t209 + t378 * t208) * qJD(3)) * t292 + (t50 + t56) * t394 + (t51 + t57) * t393 + (t20 - t336 + (t43 + t410) * t145 + (t285 * t208 + t42 - t83) * t144 + ((t100 + t258) * t144 + t285 * t145) * t209) * t392 + (t29 + t30) * t391 + (-t335 + ((t209 * t280 - t378 + t49) * t209 + (t208 * t280 + t281 + t48) * t208) * qJD(3) + t22) * t295 + (t34 + t36) * t314 / 0.2e1 + (-qJD(3) * t252 + t158 * t219 + t159 * t217 + t210 * t278 + t211 * t279) * qJD(1) + (-(qJD(1) * t273 - t143 + t241 - t38) * t39 + t27 * (t326 + t408) + t38 * (t405 * t208 + t171) + t26 * (t213 + t277 + t325) + t39 * t327 + (-t26 * t374 + t39 * (-t304 - t405)) * t209 + ((-t218 * t39 - t220 * t38) * pkin(1) + (t38 * (-t156 - t207) - t39 * t221) * t209 + (t38 * (-rSges(5,3) + t221) + t39 * (-t207 - t376)) * t208) * qJD(1)) * m(5) + (-(qJD(1) * t291 - t143 - t300 - t58) * t59 + t45 * (t208 * t298 + t204 + t321 - t386) + t44 * t407 + t59 * (t191 + t324) + (t179 * t371 - t372) * qJD(3) + ((-t218 * t59 - t220 * t58) * pkin(1) + (-pkin(2) - t265) * t369 + (t58 * (-rSges(4,3) - pkin(5)) + t59 * t298) * t208) * qJD(1)) * m(4) + (t28 + t31 + t21) * t390 - (t33 + t37 + t23) * t313 / 0.2e1 + ((t68 + t77) * t208 + (t69 + t78) * t209) * qJD(3) * t383; m(4) * t32 + m(5) * t8; ((-t314 * t347 + t318) * t208 + (t236 + (-t396 * t209 + (t348 + t238) * t208) * qJD(3)) * t209) * t295 + ((-t313 * t348 - t318) * t209 + (t236 + (t238 * t208 + (t347 - t396) * t209) * qJD(3)) * t208) * t292 + ((t217 * t323 + t219 * t322) * qJD(1) + ((t208 * t329 - t209 * t330) * t219 + (t208 * t331 + t209 * t332) * t217) * qJD(3)) * t384 + (t34 * t208 - t33 * t209 + (t68 * t208 + t209 * t69) * qJD(1)) * t383 + t237 + (qJD(1) * t36 + (-(t208 * t402 + t227 * t209) * t209 + (t208 * t403 + t226 * t209) * t208 + (t48 * t208 + t49 * t209) * qJD(1)) * t409) * t388 + (qJD(1) * t37 + (-(t227 * t208 - t209 * t402) * t209 + (t226 * t208 - t209 * t403) * t208 + (t46 * t208 + t47 * t209) * qJD(1)) * t409) * t387 + (t250 + t22) * t297 + (t249 + t23) * t296 + (t35 * t310 + (t27 * t289 + t38 * t269 + t8 * t96 + t35 * t81 + (t289 * t39 - t35 * t95) * qJD(1)) * t209 + (t26 * t289 + t39 * t269 - t8 * t95 + t35 * t82 + (t38 * t155 + t35 * t368) * qJD(1)) * t208 - (-t315 * t370 + (t264 * t219 + t35 * (-t208 ^ 2 - t209 ^ 2) * t217) * qJD(3)) * pkin(3) + t397) * m(5) + (-(t137 * t58 - t372) * qJD(1) - (t55 * (-t137 * t208 - t138 * t209) + t263 * t265) * qJD(3) + t32 * t254 + t55 * t234 + t263 * t161 + (-t44 * t208 - t45 * t209 + (-t209 * t59 + t371) * qJD(1)) * t179) * m(4); t237 + (t35 * (-t106 * t317 + t310) + t264 * t130 + (-t26 * t208 - t27 * t209 + (t208 * t38 - t370) * qJD(1)) * t155 + t397) * m(5);];
tauc = t1(:);
