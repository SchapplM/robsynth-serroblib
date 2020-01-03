% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:09
% EndTime: 2019-12-31 18:16:18
% DurationCPUTime: 6.00s
% Computational Cost: add. (4754->376), mult. (10885->493), div. (0->0), fcn. (9730->4), ass. (0->224)
t273 = sin(qJ(1));
t274 = cos(qJ(3));
t275 = cos(qJ(1));
t348 = t274 * t275;
t272 = sin(qJ(3));
t350 = t272 * t275;
t151 = Icges(5,5) * t350 - Icges(5,6) * t273 - Icges(5,3) * t348;
t363 = Icges(4,4) * t272;
t295 = Icges(4,2) * t274 + t363;
t159 = -Icges(4,6) * t273 + t295 * t275;
t454 = t151 - t159;
t264 = Icges(6,4) * t274;
t364 = Icges(6,1) * t272;
t160 = -Icges(6,5) * t275 + (-t264 + t364) * t273;
t262 = Icges(5,5) * t274;
t365 = Icges(5,1) * t272;
t162 = Icges(5,4) * t275 + (-t262 + t365) * t273;
t453 = t160 + t162;
t239 = Icges(5,5) * t348;
t163 = Icges(5,1) * t350 - Icges(5,4) * t273 - t239;
t362 = Icges(4,4) * t274;
t296 = Icges(4,1) * t272 + t362;
t165 = -Icges(4,5) * t273 + t296 * t275;
t452 = -t163 - t165;
t351 = t272 * t273;
t233 = qJ(4) * t348;
t257 = t275 * qJ(2);
t406 = pkin(1) + pkin(6);
t302 = -t406 * t273 + t257;
t322 = (rSges(5,1) + pkin(3)) * t272;
t331 = -t273 * rSges(5,2) - rSges(5,3) * t348;
t95 = t275 * t322 - t233 + t302 + t331;
t349 = t273 * t274;
t231 = qJ(4) * t349;
t245 = rSges(5,3) * t349;
t96 = -t231 - t245 + (rSges(5,2) + t406) * t275 + (qJ(2) + t322) * t273;
t303 = -t96 * t350 + t95 * t351;
t250 = rSges(6,2) * t348;
t376 = rSges(6,1) + pkin(4);
t305 = (pkin(3) + t376) * t272;
t437 = rSges(6,3) + qJ(5);
t90 = -t233 - t250 + t257 + t275 * t305 + (-t406 + t437) * t273;
t449 = t437 * t275;
t323 = rSges(6,2) * t349 + t449;
t91 = -t231 + t406 * t275 + (qJ(2) + t305) * t273 - t323;
t304 = -t91 * t350 + t90 * t351;
t407 = m(6) / 0.2e1;
t408 = m(5) / 0.2e1;
t225 = pkin(3) * t274 + qJ(4) * t272;
t199 = t273 * t225;
t227 = rSges(5,1) * t274 + rSges(5,3) * t272;
t125 = t227 * t273 + t199;
t127 = (t225 + t227) * t275;
t416 = -t125 * t275 + t127 * t273;
t368 = rSges(6,2) * t272;
t226 = rSges(6,1) * t274 + t368;
t251 = pkin(4) * t349;
t115 = t226 * t273 + t199 + t251;
t117 = (pkin(4) * t274 + t225 + t226) * t275;
t417 = -t115 * t275 + t117 * t273;
t372 = (t274 * t417 + t304) * t407 + (t274 * t416 + t303) * t408;
t194 = pkin(3) * t349 + qJ(4) * t351;
t122 = rSges(5,1) * t349 + rSges(5,3) * t351 + t194;
t330 = -pkin(3) * t348 - qJ(4) * t350;
t123 = t227 * t275 - t330;
t286 = t122 * t275 - t273 * t123;
t415 = -rSges(6,1) * t349 - rSges(6,2) * t351 - t194;
t107 = t251 - t415;
t108 = (t376 * t274 + t368) * t275 - t330;
t288 = t107 * t275 - t273 * t108;
t373 = (t288 * t274 + t304) * t407 + (t286 * t274 + t303) * t408;
t2 = t373 - t372;
t451 = t2 * qJD(1);
t228 = rSges(4,1) * t274 - rSges(4,2) * t272;
t195 = t228 * t273;
t197 = t228 * t275;
t440 = -m(6) / 0.2e1;
t441 = -m(5) / 0.2e1;
t442 = m(4) / 0.2e1;
t324 = (-t195 * t275 + t273 * t197) * t442 + t288 * t440 + t286 * t441;
t370 = t416 * t441 + t417 * t440;
t6 = t370 - t324;
t450 = t6 * qJD(1);
t291 = Icges(4,5) * t272 + Icges(4,6) * t274;
t153 = -Icges(4,3) * t273 + t291 * t275;
t157 = Icges(5,4) * t350 - Icges(5,2) * t273 - Icges(5,6) * t348;
t448 = -t153 - t157;
t149 = Icges(6,5) * t350 - Icges(6,6) * t348 + Icges(6,3) * t273;
t447 = t149 + t448;
t293 = Icges(5,4) * t272 - Icges(5,6) * t274;
t446 = t275 * (Icges(5,2) * t275 + t293 * t273) + t453 * t351;
t445 = (t452 * t272 + t454 * t274) * t275;
t241 = Icges(6,4) * t348;
t161 = Icges(6,1) * t350 + Icges(6,5) * t273 - t241;
t444 = t454 * t349 + t448 * t275 + (-t161 + t452) * t351;
t443 = (-Icges(5,4) - Icges(4,5) + Icges(6,5)) * t274 + (Icges(4,6) - Icges(5,6) + Icges(6,6)) * t272;
t270 = t273 ^ 2;
t396 = m(5) * (t96 * t273 + t275 * t95);
t384 = m(6) * (t273 * t91 + t275 * t90);
t289 = Icges(6,5) * t272 - Icges(6,6) * t274;
t148 = -Icges(6,3) * t275 + t289 * t273;
t238 = Icges(5,5) * t351;
t150 = Icges(5,6) * t275 - Icges(5,3) * t349 + t238;
t240 = Icges(6,4) * t351;
t154 = -Icges(6,2) * t349 - Icges(6,6) * t275 + t240;
t439 = -t148 * t275 + t446 + (-t150 - t154) * t349;
t158 = Icges(4,6) * t275 + t295 * t273;
t242 = Icges(4,4) * t349;
t164 = Icges(4,1) * t351 + Icges(4,5) * t275 + t242;
t282 = -t158 * t274 - t164 * t272;
t436 = -t160 * t272 + t282;
t155 = Icges(6,4) * t350 - Icges(6,2) * t348 + Icges(6,6) * t273;
t435 = t149 * t275 + t155 * t349 + t444;
t354 = t155 * t274;
t283 = t161 * t272 - t354;
t434 = t447 * t273 + t283 * t275 - t445 + t446;
t433 = -t273 * t148 + t282 * t275 - t453 * t350;
t432 = t270 / 0.2e1;
t430 = t273 / 0.2e1;
t427 = t275 / 0.2e1;
t271 = t275 ^ 2;
t329 = -t271 - t270;
t419 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (0.1e1 + t329) * t274 * t272;
t418 = t329 * t274;
t414 = t443 * t273;
t413 = t443 * t275;
t213 = -Icges(4,2) * t272 + t362;
t361 = Icges(6,4) * t272;
t215 = Icges(6,1) * t274 + t361;
t360 = Icges(5,5) * t272;
t217 = Icges(5,1) * t274 + t360;
t219 = Icges(4,1) * t274 - t363;
t358 = Icges(5,3) * t272;
t359 = Icges(6,2) * t272;
t412 = (t219 / 0.2e1 - t295 / 0.2e1 + t217 / 0.2e1 + t360 / 0.2e1 + t215 / 0.2e1 + t361 / 0.2e1 - (Icges(5,3) + Icges(6,2)) * t274 / 0.2e1) * t272 - (-t296 / 0.2e1 - t213 / 0.2e1 + t262 - t365 / 0.2e1 + t358 / 0.2e1 + t264 - t364 / 0.2e1 + t359 / 0.2e1) * t274;
t411 = 0.4e1 * qJD(1);
t410 = 2 * qJD(3);
t297 = rSges(4,1) * t272 + rSges(4,2) * t274;
t280 = -t273 * rSges(4,3) + t297 * t275;
t119 = t280 + t302;
t120 = (rSges(4,3) + t406) * t275 + (qJ(2) + t297) * t273;
t405 = m(4) * (t119 * t197 + t120 * t195);
t404 = m(4) * (t119 * t275 + t120 * t273);
t300 = t125 * t351 + t127 * t350;
t171 = t275 * (pkin(3) * t350 - t233);
t193 = pkin(3) * t351 - t231;
t57 = -t275 * (rSges(5,1) * t350 + t331) - t171 + (-rSges(5,1) * t351 - rSges(5,2) * t275 - t193 + t245) * t273;
t400 = m(5) * (-t418 * t57 + t300);
t399 = m(5) * (t122 * t96 + t123 * t95);
t397 = t274 * t396;
t301 = t115 * t351 + t117 * t350;
t49 = -t171 + (-t376 * t350 + t250) * t275 + (-t376 * t351 - t193 + t323 - t449) * t273;
t389 = m(6) * (-t418 * t49 + t301);
t387 = m(6) * (t107 * t91 + t108 * t90);
t386 = t274 * t384;
t385 = m(6) * (t273 * t90 - t275 * t91);
t382 = m(6) * (-t273 * t107 - t108 * t275);
t380 = m(6) * (t115 * t273 + t117 * t275);
t375 = m(3) * ((rSges(3,3) * t275 + t257) * t275 + (rSges(3,3) + qJ(2)) * t270);
t355 = t154 * t274;
t352 = t162 * t272;
t345 = Icges(5,1) * t349 + t150 + t238;
t344 = t217 * t275 + t151;
t343 = Icges(6,1) * t349 + t154 + t240;
t342 = t215 * t275 + t155;
t341 = -t219 * t273 + t158;
t340 = -t219 * t275 + t159;
t339 = t160 - (t264 + t359) * t273;
t338 = -Icges(6,2) * t350 + t161 - t241;
t337 = t162 - (t262 + t358) * t273;
t336 = -Icges(5,3) * t350 + t163 - t239;
t335 = -Icges(4,2) * t351 + t164 + t242;
t334 = t213 * t275 + t165;
t111 = m(6) * t418;
t328 = t111 * qJD(1);
t140 = 0.2e1 * (t432 + t271 / 0.2e1) * m(6);
t327 = t140 * qJD(1);
t325 = qJD(1) * t385;
t68 = t275 * (Icges(4,3) * t275 + t291 * t273) + t158 * t349 + t164 * t351;
t321 = rSges(6,2) * t274 - t376 * t272;
t320 = t345 * t275;
t319 = t344 * t273;
t318 = t343 * t275;
t317 = t342 * t273;
t316 = t341 * t275;
t315 = t340 * t273;
t314 = t339 * t275;
t313 = t338 * t273;
t312 = t337 * t275;
t311 = t336 * t273;
t310 = t335 * t275;
t309 = t334 * t273;
t308 = t329 * t297;
t306 = -t150 * t274 + t157;
t298 = t289 / 0.2e1 - t293 / 0.2e1 - t291 / 0.2e1;
t220 = -pkin(3) * t272 + qJ(4) * t274;
t198 = t273 * t220;
t114 = t321 * t273 + t198;
t116 = (-t220 - t321) * t275;
t287 = t114 * t273 - t116 * t275;
t222 = -rSges(5,1) * t272 + rSges(5,3) * t274;
t124 = t222 * t273 + t198;
t126 = (-t220 - t222) * t275;
t285 = t124 * t273 - t126 * t275;
t279 = -t435 * t273 / 0.2e1 + ((-t148 + t354) * t273 + (t306 - t352) * t275 - t433 + t444) * t430 - (t68 + t439) * t275 / 0.2e1 + ((-t148 - t283) * t275 + t68 + (-t355 + t153 + t436) * t273 + t434 + t445) * t427;
t278 = t153 * t432 + ((-t149 - t355 + t306) * t273 + t434 - t439) * t430 + ((t352 - t436 - t447) * t275 + t433 + t435) * t427;
t172 = t275 * t330;
t112 = 0.2e1 * t329 * t272 * (t408 + t407);
t92 = 0.4e1 * t419;
t80 = -t122 * t273 - t227 * t271 + t172;
t62 = t380 / 0.2e1;
t55 = t382 / 0.2e1;
t53 = pkin(4) * t418 - t226 * t271 + t273 * t415 + t172;
t22 = t62 - t382 / 0.2e1;
t21 = t62 + t55;
t20 = t55 - t380 / 0.2e1;
t19 = -t386 - t397;
t17 = t375 + t384 + t396 + t404;
t15 = t389 + t400;
t7 = t324 + t370;
t5 = t387 + t399 + t405 - t412;
t4 = t372 + t373;
t1 = t279 * t273 + t278 * t275;
t3 = [t17 * qJD(2) + t5 * qJD(3) + t19 * qJD(4) + qJD(5) * t385, qJD(1) * t17 + qJD(3) * t7, t5 * qJD(1) + t7 * qJD(2) + t4 * qJD(4) + t21 * qJD(5) + ((-t107 * t117 + t108 * t115 + t114 * t90 + t116 * t91) * t407 + (-t122 * t127 + t123 * t125 + t124 * t95 + t126 * t96) * t408) * t410 + ((m(4) * (t120 * t297 - t195 * t228) + t298 * t275 - t278) * t275 + (m(4) * (-t119 * t297 + t197 * t228) + t298 * t273 - t279) * t273 + (-t314 / 0.2e1 - t312 / 0.2e1 - t310 / 0.2e1 + t309 / 0.2e1 + t313 / 0.2e1 + t311 / 0.2e1) * t272 + (t318 / 0.2e1 + t320 / 0.2e1 - t316 / 0.2e1 + t315 / 0.2e1 - t317 / 0.2e1 - t319 / 0.2e1) * t274) * qJD(3), qJD(1) * t19 + qJD(3) * t4, t21 * qJD(3) + t325; -t6 * qJD(3) + t140 * qJD(5) + (-t375 / 0.4e1 - t404 / 0.4e1 - t396 / 0.4e1 - t384 / 0.4e1) * t411, 0, -t450 - t112 * qJD(4) + (t285 * t408 + t287 * t407 + t308 * t442) * t410, -t112 * qJD(3), t327; t6 * qJD(2) + t1 * qJD(3) - t2 * qJD(4) + t22 * qJD(5) + (-t387 / 0.4e1 - t399 / 0.4e1 - t405 / 0.4e1) * t411 + t412 * qJD(1), t450, t1 * qJD(1) + (m(6) * (t114 * t115 - t116 * t117 + t49 * t53) + m(5) * (t124 * t125 - t126 * t127 + t57 * t80) + m(4) * ((-t275 * t280 + (-t275 * rSges(4,3) - t297 * t273) * t273) * (-t273 * t195 - t197 * t275) + t228 * t308)) * qJD(3) + t15 * qJD(4) + ((((t313 - t314 + t311 - t312 + t309 - t310) * t274 - t414 * t273 + ((t341 - t343 - t345) * t275 + (t342 + t344 - t340) * t273) * t272) * t275 + t413 * t270) * t273 + ((t413 * t275 + (t318 - t317 + t320 - t319 - t316 + t315) * t272 + ((-t334 - t336 - t338) * t273 + (t339 + t337 + t335) * t275) * t274) * t273 - t414 * t271) * t275) * qJD(3) / 0.2e1, t15 * qJD(3) - 0.4e1 * t419 * qJD(4) - t451, t22 * qJD(1); (t397 / 0.4e1 + t386 / 0.4e1) * t411 + t2 * qJD(3) + t111 * qJD(5), 0, t451 + t92 * qJD(4) + 0.4e1 * (-t389 / 0.4e1 - t400 / 0.4e1) * qJD(3) + ((t272 * t53 + t301) * t407 + (t272 * t80 + t300) * t408 + ((-t287 + t49) * t407 + (-t285 + t57) * t408) * t274) * t410, t92 * qJD(3), t328; -t140 * qJD(2) + t20 * qJD(3) - t111 * qJD(4) - t325, -t327, t20 * qJD(1) + m(6) * (-t114 * t275 - t273 * t116) * qJD(3), -t328, 0;];
Cq = t3;
