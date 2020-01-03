% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:30
% DurationCPUTime: 6.71s
% Computational Cost: add. (7005->342), mult. (9893->461), div. (0->0), fcn. (8765->6), ass. (0->206)
t269 = sin(qJ(1));
t271 = cos(qJ(1));
t265 = qJ(3) + pkin(7);
t249 = sin(t265);
t250 = cos(t265);
t367 = rSges(6,1) + pkin(4);
t432 = rSges(6,3) + qJ(5);
t321 = t432 * t249 + t367 * t250;
t270 = cos(qJ(3));
t363 = pkin(3) * t270;
t442 = (t321 + t363) * t271;
t447 = t269 * t442;
t99 = t271 * t442;
t446 = Icges(5,3) + Icges(4,3);
t286 = Icges(5,5) * t249 + Icges(5,6) * t250;
t268 = sin(qJ(3));
t288 = Icges(4,5) * t268 + Icges(4,6) * t270;
t445 = t286 + t288;
t340 = t250 * t271;
t342 = t249 * t271;
t144 = Icges(6,5) * t342 - Icges(6,6) * t269 - Icges(6,3) * t340;
t349 = Icges(5,4) * t249;
t292 = Icges(5,2) * t250 + t349;
t150 = -Icges(5,6) * t269 + t292 * t271;
t444 = t144 - t150;
t148 = Icges(6,4) * t342 - Icges(6,2) * t269 - Icges(6,6) * t340;
t430 = -t446 * t269 + t445 * t271;
t439 = -t148 - t430;
t233 = Icges(6,5) * t340;
t152 = Icges(6,1) * t342 - Icges(6,4) * t269 - t233;
t348 = Icges(5,4) * t250;
t294 = Icges(5,1) * t249 + t348;
t154 = -Icges(5,5) * t269 + t294 * t271;
t443 = -t152 - t154;
t338 = t269 * t270;
t245 = pkin(3) * t338;
t341 = t250 * t269;
t343 = t249 * t269;
t300 = rSges(5,1) * t341 - rSges(5,2) * t343;
t133 = t245 + t300;
t225 = rSges(4,1) * t270 - rSges(4,2) * t268;
t196 = t225 * t269;
t197 = t225 * t271;
t303 = -t367 * t341 - t432 * t343;
t103 = t245 - t303;
t285 = t103 * t271 - t447;
t397 = m(5) / 0.2e1;
t213 = rSges(5,1) * t250 - rSges(5,2) * t249;
t298 = (t213 + t363) * t271;
t413 = t269 * t298;
t434 = -m(6) / 0.2e1;
t317 = m(4) * (-t196 * t271 + t269 * t197) / 0.2e1 + t285 * t434 + (-t133 * t271 + t413) * t397;
t140 = t213 * t269 + t245;
t109 = t321 * t269 + t245;
t415 = -t109 * t271 + t447;
t360 = t415 * t434 + (t140 * t271 - t413) * t397;
t12 = t360 - t317;
t441 = t12 * qJD(1);
t396 = m(6) / 0.2e1;
t412 = t271 * t298;
t358 = (-t109 * t269 - t99) * t396 + (-t140 * t269 - t412) * t397;
t361 = (t269 * t103 + t99) * t396 + (t269 * t133 + t412) * t397;
t13 = t361 - t358;
t440 = t13 * qJD(1);
t149 = Icges(5,6) * t271 + t292 * t269;
t234 = Icges(5,4) * t341;
t153 = Icges(5,1) * t343 + Icges(5,5) * t271 + t234;
t351 = Icges(4,4) * t268;
t293 = Icges(4,2) * t270 + t351;
t166 = Icges(4,6) * t271 + t293 * t269;
t244 = Icges(4,4) * t338;
t339 = t268 * t269;
t168 = Icges(4,1) * t339 + Icges(4,5) * t271 + t244;
t429 = -t149 * t250 - t153 * t249 - t166 * t270 - t168 * t268;
t421 = t269 / 0.2e1;
t438 = -t271 / 0.2e1;
t419 = t271 / 0.2e1;
t266 = t269 ^ 2;
t267 = t271 ^ 2;
t243 = t266 + t267;
t113 = (0.1e1 - t243) * t250 * t249;
t355 = m(6) * qJD(5);
t437 = t113 * t355;
t315 = t367 * t249;
t436 = Icges(4,5) * t270 - Icges(4,6) * t268 + (Icges(6,4) + Icges(5,5)) * t250 + (-Icges(5,6) + Icges(6,6)) * t249;
t167 = -Icges(4,6) * t269 + t293 * t271;
t350 = Icges(4,4) * t270;
t295 = Icges(4,1) * t268 + t350;
t169 = -Icges(4,5) * t269 + t295 * t271;
t435 = (-t167 * t270 - t169 * t268 + t443 * t249 + t444 * t250) * t271;
t433 = t149 * t341 + t153 * t343 + t166 * t338 + t168 * t339 + (t445 * t269 + t446 * t271) * t271;
t296 = rSges(5,1) * t249 + rSges(5,2) * t250;
t431 = t296 * t271;
t428 = -t167 * t338 - t169 * t339 + t439 * t271 + t444 * t341 + t443 * t343;
t247 = Icges(6,5) * t250;
t352 = Icges(6,1) * t249;
t151 = Icges(6,4) * t271 + (-t247 + t352) * t269;
t290 = Icges(6,4) * t249 - Icges(6,6) * t250;
t336 = t271 * (Icges(6,2) * t271 + t290 * t269) + t151 * t343;
t426 = t439 * t269 + t336 - t435;
t203 = -Icges(5,2) * t249 + t348;
t347 = Icges(6,5) * t249;
t205 = Icges(6,1) * t250 + t347;
t207 = Icges(5,1) * t250 - t349;
t220 = -Icges(4,2) * t268 + t350;
t346 = Icges(6,3) * t249;
t425 = (t295 / 0.2e1 + t220 / 0.2e1) * t270 - (-t294 / 0.2e1 - t203 / 0.2e1 + t247 - t352 / 0.2e1 + t346 / 0.2e1) * t250 + (t207 / 0.2e1 - t292 / 0.2e1 + t205 / 0.2e1 - Icges(6,3) * t250 / 0.2e1 + t347 / 0.2e1) * t249;
t424 = -t151 * t342 + t429 * t271;
t364 = pkin(3) * t268;
t246 = t271 * t364;
t251 = t271 * qJ(2);
t362 = -qJ(4) - pkin(6);
t276 = t246 + t251 + (-pkin(1) + t362) * t269;
t316 = -t269 * rSges(6,2) - t432 * t340;
t77 = t271 * t315 + t276 + t316;
t248 = t271 * t362;
t314 = qJ(2) + t364;
t320 = t432 * t341;
t78 = -t248 + (rSges(6,2) + pkin(1)) * t271 + (t315 + t314) * t269 - t320;
t416 = t269 * t78 + t271 * t77;
t414 = t243 * t250;
t406 = t436 * t271;
t405 = t436 * t269;
t400 = 0.2e1 * t243;
t399 = 0.4e1 * qJD(1);
t398 = 2 * qJD(3);
t395 = pkin(1) + pkin(6);
t297 = rSges(4,1) * t268 + rSges(4,2) * t270;
t274 = -t269 * rSges(4,3) + t297 * t271;
t123 = -t395 * t269 + t251 + t274;
t124 = (rSges(4,3) + t395) * t271 + (qJ(2) + t297) * t269;
t394 = m(4) * (t123 * t197 + t124 * t196);
t393 = m(4) * (t123 * t271 + t124 * t269);
t275 = -t269 * rSges(5,3) + t431;
t100 = t275 + t276;
t101 = -t248 + (rSges(5,3) + pkin(1)) * t271 + (t296 + t314) * t269;
t392 = m(5) * (t100 * t298 + t101 * t133);
t391 = m(5) * (-t100 * t269 + t271 * t101);
t390 = m(5) * (t100 * t271 + t269 * t101);
t301 = -t78 * t342 + t77 * t343;
t383 = m(6) * (t250 * t285 + t301);
t381 = m(6) * (t415 * t250 + t301);
t380 = m(6) * (t103 * t78 + t442 * t77);
t379 = m(6) * (-t269 * t77 + t271 * t78);
t378 = m(6) * t416;
t366 = m(3) * ((rSges(3,3) * t271 + t251) * t271 + (rSges(3,3) + qJ(2)) * t266);
t356 = m(6) * qJD(3);
t344 = t151 * t249;
t95 = (t434 - m(5) / 0.2e1) * t400;
t337 = t95 * qJD(1);
t232 = Icges(6,5) * t343;
t143 = Icges(6,6) * t271 - Icges(6,3) * t341 + t232;
t334 = Icges(6,1) * t341 + t143 + t232;
t333 = t205 * t271 + t144;
t332 = -t207 * t269 + t149;
t331 = -t207 * t271 + t150;
t330 = t151 - (t247 + t346) * t269;
t329 = -Icges(6,3) * t342 + t152 - t233;
t328 = -Icges(5,2) * t343 + t153 + t234;
t327 = t203 * t271 + t154;
t322 = -t432 * t250 + t315;
t107 = m(6) * t414;
t319 = t107 * qJD(1);
t36 = t416 * t250;
t318 = m(6) * t36 * qJD(1);
t313 = t334 * t271;
t312 = t333 * t269;
t311 = t332 * t271;
t310 = t331 * t269;
t309 = t330 * t271;
t308 = t329 * t269;
t307 = t328 * t271;
t306 = t327 * t269;
t305 = t243 * t297;
t304 = -t143 * t250 + t148;
t302 = t243 * t363;
t299 = -t288 / 0.2e1 - t290 / 0.2e1 - t286 / 0.2e1;
t222 = Icges(4,1) * t270 - t351;
t108 = (-t322 - t364) * t269;
t110 = t322 * t271 + t246;
t284 = t108 * t269 - t110 * t271;
t55 = -t143 * t341 + t336;
t273 = t430 * t266 / 0.2e1 + (t304 * t269 + t426 - t55) * t421 + ((t344 - t429 - t439) * t271 + t424 + t428) * t419;
t272 = -t428 * t269 / 0.2e1 + ((t304 - t344) * t271 - t424 + t428) * t421 + (t55 + t433) * t438 + ((t429 + t430) * t269 + t426 + t433 + t435) * t419;
t188 = t243 * t249;
t187 = pkin(3) * t339 - pkin(6) * t271 - t248;
t163 = t271 * (-qJ(4) * t269 + t246);
t141 = t246 + t431;
t139 = (-t296 - t364) * t269;
t94 = (-m(6) / 0.4e1 - m(5) / 0.4e1) * t400 + (m(5) + m(6)) * t243 / 0.2e1;
t91 = t109 * t343;
t42 = -t267 * t321 + t303 * t269 - t302;
t39 = -t163 + (-t367 * t342 - t316) * t271 + (-rSges(6,2) * t271 - t367 * t343 - t187 + t320) * t269;
t28 = t379 + t391;
t25 = t381 / 0.2e1;
t22 = t383 / 0.2e1;
t21 = t249 * t99 + t39 * t414 + t91;
t19 = t366 + t378 + t390 + t393;
t14 = t358 + t361;
t10 = t317 + t360;
t5 = (-t222 / 0.2e1 + t293 / 0.2e1) * t268 + t394 + t392 + t380 - t425;
t4 = t25 - t383 / 0.2e1;
t3 = t25 + t22;
t2 = t22 - t381 / 0.2e1;
t1 = t272 * t269 + t273 * t271;
t6 = [t19 * qJD(2) + t5 * qJD(3) + t28 * qJD(4) - t36 * t355, qJD(1) * t19 + qJD(3) * t10 + qJD(4) * t94, t5 * qJD(1) + t10 * qJD(2) + t14 * qJD(4) + t3 * qJD(5) + ((t108 * t77 + t110 * t78 + (-t103 + t109) * t442) * t396 + (t100 * t139 + t101 * t141 + (-t133 + t140) * t298) * t397) * t398 + ((m(4) * (t124 * t297 - t196 * t225) + t299 * t271 + (-t166 / 0.2e1 + t222 * t421) * t270 + (-t168 / 0.2e1 + Icges(4,2) * t339 / 0.2e1 - t244 / 0.2e1) * t268 - t273) * t271 + (m(4) * (-t123 * t297 + t197 * t225) + t299 * t269 + (t167 / 0.2e1 + t222 * t438) * t270 + (t169 / 0.2e1 + t220 * t419) * t268 - t272) * t269 + (-t309 / 0.2e1 - t307 / 0.2e1 + t308 / 0.2e1 + t306 / 0.2e1) * t249 + (t313 / 0.2e1 - t311 / 0.2e1 - t312 / 0.2e1 + t310 / 0.2e1) * t250) * qJD(3), qJD(1) * t28 + qJD(2) * t94 + qJD(3) * t14, t3 * qJD(3) - t318; -t12 * qJD(3) + t95 * qJD(4) + (-t378 / 0.4e1 - t390 / 0.4e1 - t393 / 0.4e1 - t366 / 0.4e1) * t399, 0, -t441 + (-m(4) * t305 / 0.2e1 + (t139 * t269 - t141 * t271) * t397 + t284 * t396) * t398 + t188 * t355, t337, t188 * t356; t12 * qJD(2) + t1 * qJD(3) - t13 * qJD(4) + t4 * qJD(5) + (-t380 / 0.4e1 - t392 / 0.4e1 - t394 / 0.4e1) * t399 + ((-t293 + t222) * t268 / 0.2e1 + t425) * qJD(1), t441, t1 * qJD(1) + (m(5) * (t140 * t139 - t298 * t141 + (-t271 * t275 - t163 + (-rSges(5,3) * t271 - t296 * t269 - t187) * t269) * (-t213 * t267 - t269 * t300 - t302)) + m(4) * ((-t271 * t274 + (-t271 * rSges(4,3) - t297 * t269) * t269) * (-t269 * t196 - t197 * t271) - t225 * t305) + m(6) * (t108 * t109 - t110 * t442 + t39 * t42)) * qJD(3) + t21 * t355 + ((((t332 - t334) * t342 + (t306 - t307 + t308 - t309) * t250 + ((-t331 + t333) * t249 + t405) * t269) * t271 - t406 * t266) * t269 + (((-t327 - t329) * t341 + (t313 - t312 - t311 + t310) * t249 + ((t328 + t330) * t250 - t406) * t271) * t269 + t405 * t267) * t271) * qJD(3) / 0.2e1, -t440, t4 * qJD(1) + t21 * t356 - t437; -t95 * qJD(2) + t13 * qJD(3) - t107 * qJD(5) + (-t379 / 0.4e1 - t391 / 0.4e1) * t399, -t337, t440 + ((t108 * t271 + t269 * t110) * t396 + (t139 * t271 + t269 * t141) * t397) * t398, 0, -t319; t2 * qJD(3) + t107 * qJD(4) + t318, 0, t2 * qJD(1) + (t91 + (t42 + t99) * t249 + (-t284 + t39) * t250 - t21) * t356 + t437, t319, t113 * t356;];
Cq = t6;
