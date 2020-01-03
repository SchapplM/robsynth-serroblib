% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:11
% EndTime: 2019-12-31 17:00:20
% DurationCPUTime: 6.16s
% Computational Cost: add. (4484->354), mult. (10504->464), div. (0->0), fcn. (9499->4), ass. (0->214)
t419 = Icges(4,1) + Icges(3,3);
t418 = Icges(4,4) - Icges(3,5);
t417 = -Icges(4,5) + Icges(3,6);
t253 = sin(qJ(2));
t255 = cos(qJ(2));
t196 = Icges(3,5) * t255 - Icges(3,6) * t253;
t198 = -Icges(4,4) * t255 + Icges(4,5) * t253;
t254 = sin(qJ(1));
t256 = cos(qJ(1));
t416 = (t196 + t198) * t256 + t419 * t254;
t344 = Icges(3,4) * t253;
t202 = Icges(3,1) * t255 - t344;
t142 = Icges(3,5) * t254 + t202 * t256;
t243 = Icges(5,6) * t253;
t338 = Icges(5,3) * t255;
t143 = Icges(5,5) * t254 + (t243 + t338) * t256;
t415 = t143 + t142;
t328 = t254 * t255;
t330 = t253 * t254;
t412 = t419 * t256 + t418 * t328 + t417 * t330;
t203 = pkin(2) * t253 - qJ(3) * t255;
t274 = rSges(5,2) * t255 - rSges(5,3) * t253;
t336 = qJ(4) * t253;
t278 = t203 - t274 + t336;
t100 = t278 * t256;
t350 = rSges(4,2) * t253;
t275 = rSges(4,3) * t255 + t350;
t308 = t203 - t275;
t114 = t308 * t254;
t116 = t308 * t256;
t329 = t253 * t256;
t327 = t255 * t256;
t248 = t256 * pkin(5);
t296 = rSges(4,1) * t256 - rSges(4,3) * t330;
t337 = qJ(3) * t253;
t88 = t248 + (-t337 - pkin(1) + (rSges(4,2) - pkin(2)) * t255) * t254 + t296;
t237 = rSges(4,2) * t327;
t347 = rSges(4,3) + qJ(3);
t360 = pkin(2) * t255;
t89 = -t237 + (rSges(4,1) + pkin(5)) * t254 + (t253 * t347 + pkin(1) + t360) * t256;
t355 = t88 * t327 + t89 * t328;
t208 = t337 + t360;
t213 = qJ(4) * t328;
t298 = -rSges(5,2) * t330 - rSges(5,3) * t328 - t213;
t364 = rSges(5,1) + pkin(3);
t78 = t248 + t364 * t256 + (-pkin(1) - t208) * t254 + t298;
t346 = rSges(5,3) + qJ(4);
t300 = pkin(2) + t346;
t348 = rSges(5,2) + qJ(3);
t79 = (pkin(5) + t364) * t254 + (t253 * t348 + t255 * t300 + pkin(1)) * t256;
t356 = t78 * t327 + t79 * t328;
t385 = m(5) / 0.2e1;
t386 = m(4) / 0.2e1;
t98 = t278 * t254;
t358 = (t100 * t330 - t98 * t329 + t356) * t385 + (-t114 * t329 + t116 * t330 + t355) * t386;
t238 = pkin(2) * t330;
t108 = t238 + (-t255 * t347 - t350) * t254;
t214 = qJ(3) * t327;
t279 = -pkin(2) * t329 + t214;
t304 = rSges(4,2) * t329 + rSges(4,3) * t327;
t109 = t279 + t304;
t95 = t238 + (t253 * t346 - t255 * t348) * t254;
t236 = rSges(5,2) * t327;
t96 = -t300 * t329 + t214 + t236;
t273 = t254 * t96 + t256 * t95;
t359 = (t253 * t273 + t356) * t385 + ((t108 * t256 + t109 * t254) * t253 + t355) * t386;
t2 = t359 - t358;
t414 = t2 * qJD(1);
t340 = Icges(4,6) * t255;
t191 = Icges(4,3) * t253 - t340;
t145 = Icges(4,5) * t254 + t191 * t256;
t413 = -t142 * t328 - t145 * t330;
t250 = t254 ^ 2;
t252 = t256 ^ 2;
t302 = t250 + t252;
t411 = t416 * t256 + t413;
t199 = Icges(3,2) * t255 + t344;
t342 = Icges(5,2) * t255;
t410 = (t243 - t342 - t199) * t256 + t415;
t228 = Icges(3,4) * t330;
t141 = Icges(3,1) * t328 - Icges(3,5) * t256 - t228;
t146 = Icges(4,5) * t256 + Icges(4,6) * t328 - Icges(4,3) * t330;
t409 = -t141 * t327 + t146 * t329 + t412 * t254;
t139 = Icges(3,4) * t328 - Icges(3,2) * t330 - Icges(3,6) * t256;
t220 = Icges(4,6) * t330;
t150 = Icges(4,4) * t256 + Icges(4,2) * t328 - t220;
t408 = t139 * t253 - t150 * t255;
t219 = Icges(5,6) * t327;
t147 = Icges(5,4) * t254 + Icges(5,2) * t329 + t219;
t197 = Icges(5,4) * t253 + Icges(5,5) * t255;
t151 = Icges(5,1) * t254 + t197 * t256;
t407 = (t145 + t147) * t329 + t415 * t327 + (t151 + t416) * t254;
t406 = (Icges(5,4) - t417) * t255 + (-Icges(5,5) + t418) * t253;
t186 = t302 * t255;
t130 = (t186 - t255) * t253;
t352 = m(5) * qJD(4);
t405 = t130 * t352;
t246 = Icges(3,4) * t255;
t343 = Icges(3,2) * t253;
t140 = Icges(3,6) * t254 + (t246 - t343) * t256;
t221 = Icges(4,6) * t329;
t149 = Icges(4,4) * t254 - Icges(4,2) * t327 + t221;
t404 = t140 * t253 + t149 * t255 + t412;
t403 = -t140 * t329 - t149 * t327 + t407;
t332 = (Icges(5,1) * t256 - Icges(5,4) * t330 - Icges(5,5) * t328) * t256;
t402 = -t332 + t407;
t401 = -t139 * t329 - t140 * t330 - t149 * t328 + t150 * t327 - t409 - t411;
t400 = -t254 / 0.2e1;
t366 = t254 / 0.2e1;
t399 = -t256 / 0.2e1;
t217 = Icges(5,6) * t330;
t144 = Icges(5,5) * t256 - Icges(5,3) * t328 - t217;
t218 = Icges(5,6) * t328;
t148 = Icges(5,4) * t256 - Icges(5,2) * t330 - t218;
t397 = (t144 * t255 + t148 * t253) * t254;
t396 = (m(4) / 0.4e1 + m(5) / 0.4e1) * t130;
t285 = (Icges(4,3) * t327 + t149 + t221) * t254;
t395 = qJD(2) * t256;
t394 = t406 * t254;
t393 = t406 * t256;
t392 = -t410 * t254 + t285;
t287 = (-Icges(5,3) * t329 + t147 + t219) * t254;
t267 = Icges(4,2) * t253 + t340;
t289 = (-t267 * t256 + t145) * t254;
t345 = Icges(3,1) * t253;
t271 = -t246 - t345;
t294 = (t271 * t256 - t140) * t254;
t391 = t287 + t289 + t294;
t390 = -t253 * (t202 / 0.2e1 - t199 / 0.2e1 - Icges(4,6) * t253 + t243 - t342 / 0.2e1 + t338 / 0.2e1 + (Icges(4,2) / 0.2e1 - Icges(4,3) / 0.2e1) * t255) - t255 * (t246 + t345 / 0.2e1 - t343 / 0.2e1 + t267 / 0.2e1 - t191 / 0.2e1 - Icges(5,6) * t255 + (-Icges(5,2) / 0.2e1 + Icges(5,3) / 0.2e1) * t253);
t389 = 0.4e1 * qJD(1);
t388 = 0.2e1 * qJD(2);
t351 = rSges(3,1) * t255;
t297 = pkin(1) + t351;
t305 = rSges(3,2) * t330 + t256 * rSges(3,3);
t112 = -t254 * t297 + t248 + t305;
t234 = rSges(3,2) * t329;
t113 = -t234 + t297 * t256 + (rSges(3,3) + pkin(5)) * t254;
t206 = rSges(3,1) * t253 + rSges(3,2) * t255;
t182 = t206 * t254;
t184 = t206 * t256;
t384 = m(3) * (t112 * t182 - t113 * t184);
t185 = t302 * t253;
t325 = -t114 * t328 - t116 * t327;
t310 = t302 * t208;
t58 = -t254 * (rSges(4,2) * t328 + t296) + t256 * (t254 * rSges(4,1) + rSges(4,3) * t329 - t237) + t310;
t380 = m(4) * (t185 * t58 + t325);
t378 = m(4) * (t108 * t88 + t109 * t89);
t377 = m(4) * (t89 * t329 - t330 * t88);
t354 = -t100 * t327 - t98 * t328;
t349 = rSges(5,2) * t253;
t50 = -t298 * t254 + t310 + (t255 * t346 + t349) * t252;
t375 = m(5) * (t185 * t50 + t354);
t261 = (-t254 * t79 - t256 * t78) * t253;
t374 = m(5) * (t255 * t273 + t261);
t372 = m(5) * (t100 * t328 - t98 * t327 + t261);
t368 = m(5) * (t78 * t95 + t79 * t96);
t367 = m(5) * (t79 * t329 - t330 * t78);
t249 = t253 ^ 2;
t251 = t255 ^ 2;
t303 = t302 * t251;
t363 = m(5) * (-t251 + (0.1e1 - t302) * t249 + t303);
t362 = m(5) * (-t186 * t255 - t249 * t302);
t361 = m(5) * (t185 * t253 + t303);
t353 = m(5) * qJD(2);
t320 = -t271 * t254 + t139;
t319 = -Icges(3,2) * t328 + t141 - t228;
t316 = Icges(5,2) * t328 + t144 - t217;
t315 = t267 * t254 + t146;
t314 = Icges(5,3) * t330 + t148 - t218;
t312 = -Icges(4,3) * t328 + t150 - t220;
t311 = t254 * (qJ(3) * t328 - t238) + t256 * t279;
t307 = -rSges(5,3) * t255 - t208 - t349;
t306 = rSges(4,2) * t255 - rSges(4,3) * t253 - t208;
t49 = t79 * t327 - t328 * t78;
t299 = m(5) * t49 * qJD(1);
t295 = t320 * t256;
t293 = t319 * t256;
t290 = t316 * t256;
t288 = t315 * t256;
t286 = t314 * t256;
t284 = t312 * t256;
t277 = -t143 * t328 - t147 * t330 + t151 * t256;
t276 = t197 / 0.2e1 + t198 / 0.2e1 + t196 / 0.2e1;
t272 = t100 * t256 + t254 * t98;
t101 = (-qJ(4) * t255 + t307) * t256;
t99 = t254 * t307 - t213;
t262 = t101 * t256 + t254 * t99 + t50;
t69 = t332 - t397;
t260 = t277 * t400 + (t69 + t412 * t256 + (t141 * t255 - t146 * t253 - t408) * t254) * t399 + (t404 * t256 - t402 + t403) * t256 / 0.2e1 + (-t144 * t327 - t148 * t329 + t404 * t254 + t277 + t401 + t411) * t366;
t259 = (t69 + t397 + t402) * t400 + t403 * t366 + ((t408 + t416) * t256 + t401 + t409 + t413) * t399;
t211 = -rSges(3,2) * t253 + t351;
t117 = t306 * t256;
t115 = t306 * t254;
t103 = t361 / 0.2e1;
t102 = t362 / 0.2e1;
t94 = t363 / 0.2e1;
t77 = 0.4e1 * t396;
t72 = t250 * t275 + t256 * t304 + t311;
t56 = t256 * (-rSges(5,3) * t329 + t236) + t274 * t250 - t302 * t336 + t311;
t54 = t103 + t94 - t362 / 0.2e1;
t53 = t102 + t103 - t363 / 0.2e1;
t52 = t102 + t94 - t361 / 0.2e1;
t25 = t372 / 0.2e1;
t21 = t50 * t186 + t253 * t272;
t19 = t374 / 0.2e1;
t18 = t367 + t377;
t15 = t375 + t380;
t8 = t368 + t378 + t384 - t390;
t7 = t25 - t374 / 0.2e1;
t6 = t25 + t19;
t5 = t19 - t372 / 0.2e1;
t3 = t358 + t359;
t1 = t254 * t260 + t256 * t259;
t4 = [t8 * qJD(2) + t18 * qJD(3) + t352 * t49, t8 * qJD(1) + t3 * qJD(3) + t6 * qJD(4) + ((-t100 * t95 + t101 * t78 + t79 * t99 - t96 * t98) * t385 + (-t108 * t116 - t109 * t114 + t115 * t89 + t117 * t88) * t386) * t388 + (m(3) * (-t112 * t211 - t182 * t206) + t276 * t256 - t259) * t395 + ((m(3) * (-t113 * t211 + t184 * t206) + t276 * t254 - t260) * t254 + (t287 / 0.2e1 + t289 / 0.2e1 + t294 / 0.2e1 + t286 / 0.2e1 + t288 / 0.2e1 + t295 / 0.2e1) * t253 + (-t285 / 0.2e1 + t290 / 0.2e1 - t284 / 0.2e1 - t293 / 0.2e1 + t410 * t366) * t255) * qJD(2), qJD(1) * t18 + qJD(2) * t3, t6 * qJD(2) + t299; t1 * qJD(2) - t2 * qJD(3) + t7 * qJD(4) + (-t368 / 0.4e1 - t378 / 0.4e1 - t384 / 0.4e1) * t389 + t390 * qJD(1), t1 * qJD(1) + t15 * qJD(3) + t21 * t352 - ((-t393 * t256 + (t286 + t288 + t295 + t391) * t255 + ((t312 - t316 + t319) * t256 + t392) * t253) * t254 + t394 * t252) * t395 / 0.2e1 + (m(4) * (-t114 * t115 - t116 * t117 + t58 * t72) + m(5) * (-t100 * t101 + t50 * t56 - t98 * t99) + m(3) * ((t254 * (rSges(3,1) * t328 - t305) + t256 * (rSges(3,1) * t327 + t254 * rSges(3,3) - t234)) * (-t254 * t182 - t184 * t256) + t302 * t211 * t206) + ((-t394 * t254 + (-t290 + t293 + t284 + t392) * t253 + ((t314 + t315 + t320) * t256 + t391) * t255) * t256 + t393 * t250) * t366) * qJD(2), t15 * qJD(2) - 0.4e1 * qJD(3) * t396 + t53 * qJD(4) - t414, t7 * qJD(1) + t53 * qJD(3) + t21 * t353 + t405; t2 * qJD(2) + (-t377 / 0.4e1 - t367 / 0.4e1) * t389, t414 + t77 * qJD(3) + t52 * qJD(4) + 0.4e1 * (-t380 / 0.4e1 - t375 / 0.4e1) * qJD(2) + ((-t255 * t72 + t325) * t386 + (-t255 * t56 + t354) * t385 + ((t115 * t254 + t117 * t256 + t58) * t386 + t262 * t385) * t253) * t388, t77 * qJD(2), t52 * qJD(2); t5 * qJD(2) - t299, t5 * qJD(1) + (t262 * t255 + (t272 + t56) * t253 - t21) * t353 + t54 * qJD(3) - t405, t54 * qJD(2), -t130 * t353;];
Cq = t4;
