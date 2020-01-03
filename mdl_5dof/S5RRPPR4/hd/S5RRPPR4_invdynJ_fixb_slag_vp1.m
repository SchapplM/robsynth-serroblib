% Calculate vector of inverse dynamics joint torques for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:47
% DurationCPUTime: 5.20s
% Computational Cost: add. (9958->401), mult. (9962->480), div. (0->0), fcn. (10012->8), ass. (0->226)
t312 = qJ(1) + qJ(2);
t286 = cos(t312);
t181 = qJD(3) * t286;
t197 = qJD(1) + qJD(2);
t285 = sin(t312);
t302 = t286 * pkin(2) + t285 * qJ(3);
t111 = t197 * t302 - t181;
t196 = qJDD(1) + qJDD(2);
t188 = t286 * pkin(3);
t195 = t197 ^ 2;
t199 = sin(qJ(1));
t201 = cos(qJ(1));
t202 = qJD(1) ^ 2;
t232 = (-qJDD(1) * t199 - t201 * t202) * pkin(1);
t367 = qJDD(3) * t285 + t197 * t181 + t232;
t212 = -t195 * t188 + t367;
t183 = t286 * qJ(3);
t273 = t285 * pkin(2);
t146 = t273 - t183;
t272 = t285 * pkin(3);
t241 = -t146 - t272;
t330 = sin(pkin(8));
t331 = cos(pkin(8));
t144 = -t285 * t330 - t286 * t331;
t145 = -t285 * t331 + t286 * t330;
t310 = t145 * rSges(5,1) - t144 * rSges(5,2);
t230 = t310 + t241;
t119 = t144 * t197;
t120 = t145 * t197;
t263 = t119 * rSges(5,1) + t120 * rSges(5,2);
t28 = (-t111 + t263) * t197 + t230 * t196 + t212;
t375 = -g(1) + t28;
t198 = sin(qJ(5));
t200 = cos(qJ(5));
t328 = Icges(6,4) * t200;
t250 = -Icges(6,2) * t198 + t328;
t71 = Icges(6,6) * t144 + t145 * t250;
t329 = Icges(6,4) * t198;
t252 = Icges(6,1) * t200 - t329;
t74 = Icges(6,5) * t144 + t145 * t252;
t369 = t198 * t71 - t200 * t74;
t248 = Icges(6,5) * t200 - Icges(6,6) * t198;
t68 = Icges(6,3) * t144 + t145 * t248;
t24 = t144 * t68 - t369 * t145;
t136 = t144 * pkin(7);
t317 = t145 * t198;
t277 = -rSges(6,2) * t317 + t144 * rSges(6,3);
t316 = t145 * t200;
t77 = rSges(6,1) * t316 + t277;
t371 = -t145 * pkin(4) - t136 - t77;
t70 = Icges(6,3) * t145 - t144 * t248;
t372 = t145 * t70;
t225 = -t273 - t272;
t221 = t183 + t225;
t370 = t221 + t310;
t256 = -t198 * t74 - t200 * t71;
t297 = qJD(5) * t198;
t234 = t120 * t200 + t144 * t297;
t296 = qJD(5) * t200;
t235 = -t120 * t198 + t144 * t296;
t48 = rSges(6,1) * t234 + rSges(6,2) * t235 + t119 * rSges(6,3);
t368 = t120 * pkin(4) + t119 * pkin(7) + t48;
t247 = Icges(6,5) * t198 + Icges(6,6) * t200;
t91 = t247 * t144;
t90 = t247 * t145;
t101 = -t144 * rSges(5,1) - t145 * rSges(5,2);
t350 = t188 + t302;
t81 = t101 + t350;
t237 = t119 * t198 + t145 * t296;
t236 = -t119 * t200 + t145 * t297;
t334 = pkin(1) * qJD(1);
t294 = t199 * t334;
t148 = t285 * rSges(3,1) + t286 * rSges(3,2);
t315 = t148 * t197;
t121 = -t294 - t315;
t318 = t144 * t200;
t319 = t144 * t198;
t79 = -rSges(6,1) * t318 + rSges(6,2) * t319 + t145 * rSges(6,3);
t332 = -t144 * pkin(4) + pkin(7) * t145 + t79;
t52 = t350 + t332;
t180 = qJD(3) * t285;
t265 = t197 * t285;
t194 = t201 * pkin(1);
t342 = t199 * pkin(1);
t274 = qJDD(1) * t194 - t202 * t342;
t163 = t197 * t183;
t307 = t163 + t180;
t216 = -qJDD(3) * t286 + t196 * t302 + t274 + (-pkin(2) * t265 + t180 + t307) * t197;
t208 = t196 * t188 - t195 * t272 + t216;
t311 = t120 * rSges(5,1) - t119 * rSges(5,2);
t29 = t196 * t101 + t197 * t311 + t208;
t366 = -g(2) + t29;
t262 = rSges(6,1) * t198 + rSges(6,2) * t200;
t335 = rSges(6,1) * t200;
t174 = rSges(6,2) * t198 - t335;
t156 = t174 * qJD(5);
t299 = qJD(5) * t156;
t88 = qJD(5) * t119 + qJDD(5) * t145;
t14 = -t145 * t299 + t332 * t196 + t197 * t368 + t88 * t262 + t208;
t365 = t14 - g(2);
t303 = t286 * rSges(4,1) + t285 * rSges(4,3);
t185 = t286 * rSges(4,3);
t269 = t285 * rSges(4,1);
t147 = t269 - t185;
t351 = -t146 - t147;
t50 = -t197 * t111 - t195 * t303 + t196 * t351 + t367;
t364 = t50 - g(1);
t266 = t197 * t286;
t171 = rSges(4,3) * t266;
t51 = t196 * t303 + t197 * (-rSges(4,1) * t265 + t171) + t216;
t363 = t51 - g(2);
t129 = rSges(3,1) * t266 - rSges(3,2) * t265;
t362 = -t129 * t197 - t148 * t196 - g(1) + t232;
t151 = t286 * rSges(3,1) - t285 * rSges(3,2);
t361 = t151 * t196 - t197 * t315 - g(2) + t274;
t249 = Icges(6,2) * t200 + t329;
t251 = Icges(6,1) * t198 + t328;
t245 = -t198 * t249 + t200 * t251;
t56 = t245 * t145 + t91;
t359 = t197 * t56;
t227 = t241 - t371;
t298 = qJD(5) * t262;
t291 = t144 * t298;
t264 = t180 + t291;
t34 = t197 * t227 + t264 - t294;
t293 = t201 * t334;
t270 = -t181 + t293;
t346 = t145 * t298 + t197 * t52;
t35 = t270 + t346;
t358 = (t35 * t225 - t34 * t350) * t197;
t47 = rSges(6,1) * t236 + rSges(6,2) * t237 + t120 * rSges(6,3);
t357 = t119 * pkin(4) - t120 * pkin(7) - t47;
t356 = -t197 * t310 + t311;
t355 = t197 * t81;
t354 = t197 * t272;
t349 = t302 + t303;
t353 = t197 * t349;
t352 = t197 * t147 + t171;
t142 = t197 * t146;
t348 = t307 - t180 + t142;
t347 = t371 * t197 + t307 + t368;
t337 = t249 * t145 - t74;
t339 = -t251 * t145 - t71;
t345 = t198 * t337 + t200 * t339;
t344 = t88 / 0.2e1;
t89 = qJD(5) * t120 - qJDD(5) * t144;
t343 = t89 / 0.2e1;
t76 = Icges(6,5) * t145 - t144 * t252;
t341 = -t144 * t70 - t76 * t316;
t340 = -t76 * t318 + t372;
t73 = Icges(6,6) * t145 - t144 * t250;
t338 = -t251 * t144 + t73;
t336 = t249 * t144 + t76;
t333 = t198 * t73;
t314 = t248 * t197;
t306 = -t249 + t252;
t305 = -t250 - t251;
t301 = qJD(5) * t144;
t300 = qJD(5) * t145;
t284 = -t301 / 0.2e1;
t283 = t301 / 0.2e1;
t282 = -t300 / 0.2e1;
t281 = t300 / 0.2e1;
t276 = -t286 / 0.2e1;
t275 = t285 / 0.2e1;
t271 = t180 - t294;
t254 = t198 * t76 + t200 * t73;
t44 = Icges(6,4) * t234 + Icges(6,2) * t235 + Icges(6,6) * t119;
t46 = Icges(6,1) * t234 + Icges(6,4) * t235 + Icges(6,5) * t119;
t218 = t254 * qJD(5) + t198 * t44 - t200 * t46;
t43 = Icges(6,4) * t236 + Icges(6,2) * t237 + Icges(6,6) * t120;
t45 = Icges(6,1) * t236 + Icges(6,4) * t237 + Icges(6,5) * t120;
t219 = t256 * qJD(5) + t198 * t43 - t200 * t45;
t253 = -t200 * t76 + t333;
t41 = Icges(6,5) * t236 + Icges(6,6) * t237 + Icges(6,3) * t120;
t42 = Icges(6,5) * t234 + Icges(6,6) * t235 + Icges(6,3) * t119;
t268 = -t144 * (-t119 * t369 - t120 * t68 - t144 * t41 + t219 * t145) + t145 * (t253 * t119 + t120 * t70 - t144 * t42 + t218 * t145);
t267 = -t144 * (-t119 * t68 + t120 * t369 + t219 * t144 + t145 * t41) + t145 * (t119 * t70 - t253 * t120 + t218 * t144 + t145 * t42);
t175 = rSges(2,1) * t201 - rSges(2,2) * t199;
t173 = rSges(2,1) * t199 + rSges(2,2) * t201;
t25 = t317 * t73 + t341;
t261 = -t144 * t24 + t145 * t25;
t26 = -t145 * t68 + t318 * t74 - t319 * t71;
t27 = t319 * t73 + t340;
t260 = -t144 * t26 + t145 * t27;
t259 = -t144 * t34 - t145 * t35;
t258 = -t144 * t48 - t145 * t47;
t257 = t144 * t79 - t145 * t77;
t246 = -t198 * t251 - t200 * t249;
t244 = -t142 + t271;
t243 = t163 + t271;
t242 = t181 + t263;
t229 = t336 * t198 + t338 * t200;
t226 = (t305 * t198 + t306 * t200) * t197;
t224 = -t269 - t273;
t220 = -pkin(3) * t265 + t244;
t112 = t183 + t185 + t224;
t154 = t250 * qJD(5);
t155 = t252 * qJD(5);
t217 = t246 * qJD(5) - t154 * t198 + t155 * t200;
t214 = t181 + t357;
t57 = t144 * t245 - t90;
t53 = t57 * t197;
t10 = qJD(5) * t260 + t53;
t17 = -qJD(5) * t369 - t198 * t45 - t200 * t43;
t18 = qJD(5) * t253 - t198 * t46 - t200 * t44;
t153 = t248 * qJD(5);
t22 = t119 * t245 - t120 * t247 + t144 * t153 + t145 * t217;
t23 = -t119 * t247 - t120 * t245 + t144 * t217 - t145 * t153;
t9 = qJD(5) * t261 + t359;
t211 = (t53 + ((t25 - t26 - t341) * t144 + t340 * t145) * qJD(5)) * t283 + (t245 * qJD(5) + t154 * t200 + t155 * t198) * t197 + (-t254 + t57) * t344 + (-t256 + t56) * t343 + (t18 + t23) * t281 + (t17 + t22 + t10) * t284 + (-t359 + ((t26 + (-t253 + t68) * t145) * t145 + (t27 + (t68 - t333) * t144 + t372 - t340) * t144) * qJD(5) + t9) * t282 + (-t246 + Icges(5,3) + Icges(3,3) + Icges(4,2)) * t196;
t209 = t221 - t371;
t206 = t136 + (pkin(4) + t335) * t145 + t221 + t277;
t61 = t197 * t230 + t271;
t62 = t270 + t355;
t204 = (t62 * t225 - t350 * t61) * t197;
t86 = t197 * t351 + t271;
t87 = t270 + t353;
t203 = (t87 * t224 - t349 * t86) * t197;
t122 = t151 * t197 + t293;
t97 = t262 * t144;
t96 = t262 * t145;
t36 = qJD(5) * t257 - qJD(4);
t13 = -t144 * t299 - t89 * t262 + (-t111 + t357) * t197 + t227 * t196 + t212;
t12 = qJD(5) * t258 + t77 * t88 + t79 * t89 + qJDD(4);
t1 = [Icges(2,3) * qJDD(1) + t211 + (t361 * (t194 + t151) + t362 * (-t148 - t342) + (-t129 - t293 + t122) * t121) * m(3) + ((t173 ^ 2 + t175 ^ 2) * qJDD(1) + g(1) * t173 - g(2) * t175) * m(2) + (t13 * (t206 - t342) + t34 * (t214 - t293) - g(1) * (t209 - t342) + t365 * (t194 + t52) + (-t220 + t34 - t291 - t294 + t347) * t35 + t358) * m(6) + (t61 * (t242 - t293) + t204 + t366 * (t194 + t81) + (t243 + t61 - t220 + t356) * t62 + t375 * (t370 - t342)) * m(5) + (-t86 * t270 + t203 + (t243 - t244 + t86 + t352) * t87 + t363 * (t194 + t349) + t364 * (t112 - t342)) * m(4); t211 + (-g(1) * t209 + t13 * t206 + t365 * t52 + (t142 - t264 + t347 + t354) * t35 + (-t181 + t214 + t346) * t34 + t358) * m(6) + (t204 + t366 * t81 + (t348 + t354 + t356) * t62 + (-t181 + t242 + t355) * t61 + t375 * t370) * m(5) + (t86 * t353 + t203 + (t348 + t352) * t87 + t363 * t349 + t364 * t112) * m(4) + (-t121 * t129 - t122 * t315 + (t121 * t197 + t361) * t151 + (t122 * t197 - t362) * t148) * m(3); (-m(4) - m(5) - m(6)) * (g(1) * t285 - g(2) * t286) + 0.2e1 * (t13 * t275 + t14 * t276) * m(6) + 0.2e1 * (t275 * t28 + t276 * t29) * m(5) + 0.2e1 * (t275 * t50 + t276 * t51) * m(4); (t12 + g(3)) * m(6) + (qJDD(4) + g(3)) * m(5); t119 * t10 / 0.2e1 + t145 * (t267 * qJD(5) + t196 * t57 + t197 * t23 + t26 * t89 + t27 * t88) / 0.2e1 + t260 * t344 + (t119 * t27 + t120 * t26 + t267) * t281 + t120 * t9 / 0.2e1 - t144 * (t268 * qJD(5) + t196 * t56 + t197 * t22 + t24 * t89 + t25 * t88) / 0.2e1 + t261 * t343 + (t119 * t25 + t120 * t24 + t268) * t284 + t196 * (t144 * t256 - t145 * t254) / 0.2e1 + t197 * (-t119 * t254 - t120 * t256 - t144 * t17 + t145 * t18) / 0.2e1 + ((t91 * t300 - t314) * t145 + (t226 + (-t345 * t144 + (-t90 + t229) * t145) * qJD(5)) * t144) * t282 + ((t90 * t301 + t314) * t144 + (t226 + (t229 * t145 + (-t345 - t91) * t144) * qJD(5)) * t145) * t283 - t197 * ((t306 * t198 - t305 * t200) * t197 + ((t144 * t337 - t145 * t336) * t200 + (-t144 * t339 + t145 * t338) * t198) * qJD(5)) / 0.2e1 + (-t12 * t257 + t36 * (-t119 * t77 - t120 * t79 - t258) + t259 * t156 - (-t119 * t35 + t34 * t120 - t13 * t144 - t14 * t145) * t262 - (-t34 * t96 + t35 * t97) * t197 - (t36 * (t144 * t97 + t145 * t96) + t259 * t174) * qJD(5) - g(1) * t97 - g(2) * t96 - g(3) * t174) * m(6);];
tau = t1;
