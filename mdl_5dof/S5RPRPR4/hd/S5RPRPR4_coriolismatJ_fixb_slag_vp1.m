% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR4
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:10
% EndTime: 2019-12-05 17:53:24
% DurationCPUTime: 5.70s
% Computational Cost: add. (22019->409), mult. (14381->522), div. (0->0), fcn. (13192->10), ass. (0->244)
t287 = qJ(1) + pkin(8);
t281 = cos(t287);
t286 = qJ(3) + pkin(9);
t282 = qJ(5) + t286;
t273 = sin(t282);
t274 = cos(t282);
t222 = rSges(6,1) * t273 + rSges(6,2) * t274;
t278 = sin(t286);
t399 = pkin(4) * t278;
t289 = sin(qJ(3));
t400 = pkin(3) * t289;
t248 = -t399 - t400;
t435 = (t222 - t248) * t281;
t120 = t281 * t435;
t279 = sin(t287);
t353 = t279 * t289;
t265 = pkin(3) * t353;
t125 = t265 + (t222 + t399) * t279;
t280 = cos(t286);
t229 = rSges(5,1) * t278 + rSges(5,2) * t280;
t173 = t229 * t279 + t265;
t425 = m(6) / 0.2e1;
t426 = m(5) / 0.2e1;
t310 = (t229 + t400) * t281;
t450 = t281 * t310;
t384 = (-t173 * t279 - t450) * t426 + (-t125 * t279 - t120) * t425;
t359 = t274 * t279;
t363 = t273 * t279;
t189 = rSges(6,1) * t363 + rSges(6,2) * t359;
t129 = -t248 * t279 + t189;
t354 = t279 * t280;
t357 = t278 * t279;
t324 = rSges(5,1) * t357 + rSges(5,2) * t354;
t158 = t265 + t324;
t395 = (t129 * t279 + t120) * t425 + (t158 * t279 + t450) * t426;
t23 = t395 - t384;
t461 = t23 * qJD(1);
t190 = t222 * t281;
t165 = t281 * t190;
t446 = t279 * t189 + t165;
t295 = m(6) * t446;
t358 = t274 * t281;
t362 = t273 * t281;
t153 = Icges(6,5) * t358 - Icges(6,6) * t362 + Icges(6,3) * t279;
t267 = Icges(6,4) * t274;
t440 = Icges(6,2) * t273 - t267;
t154 = Icges(6,6) * t281 + t440 * t279;
t217 = Icges(6,5) * t274 - Icges(6,6) * t273;
t371 = t217 * t279;
t346 = t281 * (Icges(6,3) * t281 - t371) + t154 * t363;
t155 = Icges(6,4) * t358 - Icges(6,2) * t362 + Icges(6,6) * t279;
t234 = Icges(6,4) * t362;
t157 = Icges(6,1) * t358 + Icges(6,5) * t279 - t234;
t447 = (t155 * t273 - t157 * t274) * t281;
t460 = t153 * t279 + t346 - t447;
t459 = t154 * t362;
t276 = t279 ^ 2;
t277 = t281 ^ 2;
t322 = t276 + t277;
t458 = t222 * t322;
t291 = cos(qJ(3));
t457 = Icges(4,5) * t289 + Icges(5,5) * t278 + Icges(4,6) * t291 + Icges(5,6) * t280;
t456 = -t279 / 0.2e1;
t409 = t279 / 0.2e1;
t455 = -t280 / 0.2e1;
t406 = t281 / 0.2e1;
t454 = -t291 / 0.2e1;
t337 = -Icges(6,2) * t358 + t157 - t234;
t233 = Icges(6,4) * t363;
t156 = -Icges(6,1) * t359 + Icges(6,5) * t281 + t233;
t338 = Icges(6,2) * t359 + t156 + t233;
t443 = Icges(6,1) * t273 + t267;
t339 = t281 * t443 + t155;
t340 = -t279 * t443 + t154;
t453 = -(t337 * t279 + t338 * t281) * t273 - (t339 * t279 + t340 * t281) * t274;
t452 = -(Icges(5,1) * t280 / 0.2e1 - Icges(5,4) * t278 + Icges(5,2) * t455) * t278 - (Icges(4,1) * t291 / 0.2e1 - Icges(4,4) * t289 + Icges(4,2) * t454) * t289;
t390 = rSges(6,1) * t274;
t223 = -rSges(6,2) * t273 + t390;
t369 = t223 * t279;
t302 = Icges(6,5) * t273 + Icges(6,6) * t274;
t177 = t279 * t302;
t178 = t302 * t281;
t397 = (-t276 * t178 + (t279 * t177 + t453) * t281) * t409 + (t277 * t177 + (-t281 * t178 - t453) * t279) * t406;
t284 = t291 * pkin(3);
t275 = t284 + pkin(2);
t249 = t281 * t275;
t288 = -qJ(4) - pkin(6);
t137 = t281 * (t281 * pkin(2) - t249 + (pkin(6) + t288) * t279);
t312 = -rSges(6,1) * t358 + rSges(6,2) * t362;
t142 = t281 * (t279 * rSges(6,3) - t312);
t272 = t281 * pkin(6);
t350 = t281 * t288;
t150 = -t350 - t272 + (pkin(2) - t275) * t279;
t326 = -rSges(6,2) * t363 - t281 * rSges(6,3);
t159 = -rSges(6,1) * t359 - t326;
t398 = pkin(4) * t280;
t239 = t275 + t398;
t216 = t281 * t239;
t47 = -t281 * (-t216 + t249) - t137 + t142 + (-t150 - t159 + (t239 - t275) * t279) * t279;
t6 = t397 + m(6) * (t223 * t120 + t125 * t369 - t446 * t47);
t451 = t6 * qJD(5);
t348 = t281 * t291;
t349 = t281 * t289;
t186 = Icges(4,4) * t348 - Icges(4,2) * t349 + Icges(4,6) * t279;
t262 = Icges(4,4) * t349;
t188 = Icges(4,1) * t348 + Icges(4,5) * t279 - t262;
t449 = (t186 * t289 - t188 * t291) * t281;
t351 = t280 * t281;
t356 = t278 * t281;
t169 = Icges(5,4) * t351 - Icges(5,2) * t356 + Icges(5,6) * t279;
t243 = Icges(5,4) * t356;
t171 = Icges(5,1) * t351 + Icges(5,5) * t279 - t243;
t448 = (t169 * t278 - t171 * t280) * t281;
t445 = t457 * t279;
t444 = t457 * t281;
t268 = Icges(5,4) * t280;
t442 = Icges(5,1) * t278 + t268;
t283 = Icges(4,4) * t291;
t441 = Icges(4,1) * t289 + t283;
t439 = Icges(5,2) * t278 - t268;
t438 = Icges(4,2) * t289 - t283;
t378 = Icges(6,4) * t273;
t218 = Icges(6,2) * t274 + t378;
t221 = Icges(6,1) * t274 - t378;
t434 = (-t440 + t443) * t273 + (t218 - t221) * t274;
t315 = (-t440 / 0.2e1 + t443 / 0.2e1) * t274 + (-t218 / 0.2e1 + t221 / 0.2e1) * t273;
t318 = -t156 * t274 - t153;
t374 = t154 * t273;
t65 = -t156 * t359 + t346;
t66 = t281 * t153 + t155 * t363 - t157 * t359;
t5 = ((t447 + t460) * t281 + ((t318 - t374) * t281 + t66 + t459) * t279) * t409 + (t279 * t66 + t281 * t65) * t456 + ((t66 + (-t153 + t374) * t281 - t459) * t281 + (t318 * t279 + t460 - t65) * t279) * t406;
t329 = -Icges(4,2) * t348 + t188 - t262;
t331 = t281 * t441 + t186;
t333 = -Icges(5,2) * t351 + t171 - t243;
t335 = t281 * t442 + t169;
t431 = t333 * t278 + t335 * t280 + t329 * t289 + t331 * t291;
t261 = Icges(4,4) * t353;
t352 = t279 * t291;
t187 = -Icges(4,1) * t352 + Icges(4,5) * t281 + t261;
t330 = Icges(4,2) * t352 + t187 + t261;
t185 = Icges(4,6) * t281 + t438 * t279;
t332 = -t279 * t441 + t185;
t242 = Icges(5,4) * t357;
t170 = -Icges(5,1) * t354 + Icges(5,5) * t281 + t242;
t334 = Icges(5,2) * t354 + t170 + t242;
t168 = Icges(5,6) * t281 + t439 * t279;
t336 = -t279 * t442 + t168;
t430 = -t334 * t278 - t336 * t280 - t330 * t289 - t332 * t291;
t429 = 0.4e1 * qJD(1);
t428 = 2 * qJD(3);
t427 = m(4) / 0.2e1;
t392 = rSges(4,1) * t291;
t320 = pkin(2) + t392;
t323 = -rSges(4,2) * t353 - t281 * rSges(4,3);
t402 = sin(qJ(1)) * pkin(1);
t131 = t320 * t279 - t272 + t323 + t402;
t264 = rSges(4,2) * t349;
t401 = cos(qJ(1)) * pkin(1);
t132 = -t401 + t264 - t320 * t281 + (-rSges(4,3) - pkin(6)) * t279;
t255 = rSges(4,1) * t289 + rSges(4,2) * t291;
t212 = t255 * t279;
t213 = t255 * t281;
t424 = m(4) * (-t131 * t212 + t132 * t213);
t325 = -rSges(5,2) * t357 - t281 * rSges(5,3);
t391 = rSges(5,1) * t280;
t115 = t402 + t350 + (t275 + t391) * t279 + t325;
t313 = -rSges(5,1) * t351 + rSges(5,2) * t356;
t116 = -t401 - t249 + (-rSges(5,3) + t288) * t279 + t313;
t423 = m(5) * (-t115 * t158 + t116 * t310);
t422 = m(5) * (-t281 * t115 - t116 * t279);
t285 = -pkin(7) + t288;
t110 = t402 + t281 * t285 + (t239 + t390) * t279 + t326;
t105 = t281 * t110;
t111 = -t401 - t216 + (-rSges(6,3) + t285) * t279 + t312;
t311 = t223 * t105 + t111 * t369;
t418 = m(6) * (t125 * t190 - t189 * t435 + t311);
t417 = m(6) * ((-t129 * t281 + t279 * t435) * t222 + t311);
t416 = m(6) * (-t110 * t129 + t111 * t435);
t415 = m(6) * (-t110 * t189 + t111 * t190);
t414 = m(6) * (-t111 * t279 - t105);
t407 = -t281 / 0.2e1;
t393 = m(6) * qJD(5);
t373 = t168 * t278;
t372 = t185 * t289;
t297 = -m(6) * t458 / 0.2e1;
t75 = -t295 / 0.2e1 + t297;
t347 = t75 * qJD(1);
t304 = Icges(5,5) * t280 - Icges(5,6) * t278;
t166 = Icges(5,3) * t281 - t304 * t279;
t344 = t281 * t166 + t168 * t357;
t343 = t279 * t166 + t170 * t351;
t306 = Icges(4,5) * t291 - Icges(4,6) * t289;
t183 = Icges(4,3) * t281 - t306 * t279;
t342 = t281 * t183 + t185 * t353;
t341 = t279 * t183 + t187 * t348;
t319 = t223 + t398;
t167 = Icges(5,5) * t351 - Icges(5,6) * t356 + Icges(5,3) * t279;
t317 = -t170 * t280 - t167;
t184 = Icges(4,5) * t348 - Icges(4,6) * t349 + Icges(4,3) * t279;
t316 = -t187 * t291 - t184;
t72 = t281 * t167 + t169 * t357 - t171 * t354;
t81 = t281 * t184 + t186 * t353 - t188 * t352;
t294 = -t5 + (-t339 * t273 + t337 * t274 - t434 * t281 + t371) * t409 + (t217 * t281 - t340 * t273 + t338 * t274 + t434 * t279) * t406;
t293 = -t315 + (t406 + t407) * (t155 * t274 + t157 * t273);
t266 = pkin(3) * t352;
t257 = -rSges(4,2) * t289 + t392;
t230 = -rSges(5,2) * t278 + t391;
t176 = (-t230 - t284) * t281;
t174 = t230 * t279 + t266;
t128 = (-t319 - t284) * t281;
t126 = t319 * t279 + t266;
t118 = -t212 * t279 - t213 * t281;
t103 = t446 * t393;
t98 = -t279 * t159 + t142;
t95 = -t229 * t277 - t279 * t324 - t322 * t400;
t83 = t184 * t279 - t449;
t82 = -t185 * t349 + t341;
t80 = -t187 * t352 + t342;
t76 = t295 / 0.2e1 + t297;
t74 = t167 * t279 - t448;
t73 = -t168 * t356 + t343;
t71 = -t170 * t354 + t344;
t60 = t248 * t277 - t165 + (t265 - t129) * t279 + (-t322 + t277) * t400;
t53 = t279 * t83 + t281 * t82;
t52 = t279 * t81 + t281 * t80;
t49 = t279 * t74 + t281 * t73;
t48 = t279 * t72 + t281 * t71;
t44 = t414 + t422;
t37 = t315 + t415;
t35 = t417 / 0.2e1;
t33 = t418 / 0.2e1;
t25 = t384 + t395;
t22 = (t81 + (-t184 + t372) * t281 - t341) * t281 + (t316 * t279 + t342 - t80) * t279;
t21 = (t342 + t83 + t449) * t281 + (-t82 + (t316 - t372) * t281 + t81 + t341) * t279;
t19 = (t72 + (-t167 + t373) * t281 - t343) * t281 + (t317 * t279 + t344 - t71) * t279;
t18 = (t344 + t74 + t448) * t281 + (-t73 + (t317 - t373) * t281 + t72 + t343) * t279;
t17 = (t441 / 0.2e1 - t438 / 0.2e1) * t291 + (t442 / 0.2e1 - t439 / 0.2e1) * t280 + t424 + t423 + t416 + t315 - t452;
t8 = m(6) * (t223 * t458 - t446 * t98) + t397;
t7 = t8 * qJD(5);
t4 = t33 - t417 / 0.2e1 + t5;
t3 = t35 - t418 / 0.2e1 + t5;
t2 = t33 + t35 + t294;
t1 = (t19 / 0.2e1 + t49 / 0.2e1 + t53 / 0.2e1 + t22 / 0.2e1) * t281 + (-t48 / 0.2e1 + t18 / 0.2e1 + t21 / 0.2e1 - t52 / 0.2e1) * t279 + t5;
t9 = [qJD(3) * t17 + qJD(4) * t44 + qJD(5) * t37, 0, t17 * qJD(1) + t25 * qJD(4) + t2 * qJD(5) + ((-t115 * t176 + t116 * t174 + (-t158 + t173) * t310) * t426 + (-t110 * t128 + t111 * t126 + (t125 - t129) * t435) * t425 + ((t131 * t281 + t132 * t279) * t257 + (-t212 * t281 + t213 * t279) * t255) * t427) * t428 + (t294 + (t18 + t21) * t456 + (-t336 * t278 + t334 * t280 - t332 * t289 + t330 * t291) * t406 + (t306 + t304) * (t277 / 0.2e1 + t276 / 0.2e1) + (-t335 * t278 + t333 * t280 - t331 * t289 + t329 * t291 + t48 + t52) * t409 + (t19 + t49 + t53 + t22) * t407) * qJD(3), qJD(1) * t44 + qJD(3) * t25 + qJD(5) * t76, t37 * qJD(1) + t2 * qJD(3) + t76 * qJD(4) + t294 * qJD(5) + ((-t189 * t281 + t190 * t279) * t222 + t311) * t393; 0, 0, -t103 + (t118 * t427 + t60 * t425 + t95 * t426) * t428, 0, -qJD(3) * t295 - t103; (t293 + (-t439 + t442) * t455 + (-t438 + t441) * t454 + t452) * qJD(1) + t1 * qJD(3) - t23 * qJD(4) + t4 * qJD(5) + (-t424 / 0.4e1 - t416 / 0.4e1 - t423 / 0.4e1) * t429, 0, t1 * qJD(1) + (m(6) * (t125 * t126 - t128 * t435 + t47 * t60) + m(5) * (t173 * t174 - t310 * t176 + (t281 * (t279 * rSges(5,3) - t313) - t137 + (rSges(5,1) * t354 - t150 + t325) * t279) * t95) + m(4) * (t322 * t257 * t255 + (t281 * (rSges(4,1) * t348 + t279 * rSges(4,3) - t264) - t279 * (-rSges(4,1) * t352 - t323)) * t118) + t397 + ((t430 * t281 + (-t431 + t445) * t279) * t281 - t444 * t276) * t409 + ((t431 * t279 + (-t430 - t444) * t281) * t279 + t445 * t277) * t406) * qJD(3) + t451, -t461, t4 * qJD(1) + t6 * qJD(3) + t451; t23 * qJD(3) - t75 * qJD(5) + (-t422 / 0.4e1 - t414 / 0.4e1) * t429, 0, t461 + ((t126 * t281 + t128 * t279) * t425 + (t174 * t281 + t176 * t279) * t426) * t428, 0, -t347; (t293 - t415) * qJD(1) + t3 * qJD(3) + t75 * qJD(4) + t5 * qJD(5), 0, t3 * qJD(1) + ((t60 * t98 + (t126 * t279 - t128 * t281) * t222) * m(6) + t397) * qJD(3) + t7, t347, qJD(1) * t5 + qJD(3) * t8 + t7;];
Cq = t9;
