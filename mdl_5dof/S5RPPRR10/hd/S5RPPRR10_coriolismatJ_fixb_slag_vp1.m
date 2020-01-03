% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:00
% EndTime: 2019-12-31 18:04:08
% DurationCPUTime: 5.02s
% Computational Cost: add. (17333->348), mult. (26444->491), div. (0->0), fcn. (29299->8), ass. (0->216)
t276 = sin(pkin(8));
t277 = cos(pkin(8));
t356 = qJ(4) + qJ(5);
t323 = sin(t356);
t324 = cos(t356);
t253 = t276 * t324 - t277 * t323;
t280 = cos(qJ(1));
t233 = t253 * t280;
t284 = t276 * t323 + t277 * t324;
t234 = t284 * t280;
t279 = sin(qJ(1));
t172 = t234 * rSges(6,1) + t233 * rSges(6,2) - t279 * rSges(6,3);
t397 = cos(qJ(4));
t269 = t397 * pkin(4) + pkin(3);
t278 = sin(qJ(4));
t281 = -pkin(7) - pkin(6);
t110 = (qJ(2) + t281) * t279 + (pkin(1) + (pkin(2) + t269) * t277 + (pkin(4) * t278 + qJ(3)) * t276) * t280 + t172;
t231 = t253 * t279;
t232 = t284 * t279;
t181 = t231 * rSges(6,1) - t232 * rSges(6,2);
t182 = t233 * rSges(6,1) - t234 * rSges(6,2);
t270 = t280 * qJ(2);
t320 = qJ(3) * t276 + pkin(1);
t363 = t276 * t278;
t328 = pkin(4) * t363;
t359 = t277 * t279;
t333 = -t269 * t359 - t279 * t328;
t433 = -t232 * rSges(6,1) - t231 * rSges(6,2);
t444 = t270 + (-rSges(6,3) + t281) * t280 + (-pkin(2) * t277 - t320) * t279 + t333 + t433;
t458 = m(6) * (t110 * t182 - t181 * t444);
t326 = t276 * t397;
t360 = t277 * t278;
t257 = t326 - t360;
t292 = t277 * t397 + t363;
t221 = rSges(5,1) * t257 - rSges(5,2) * t292;
t274 = t279 ^ 2;
t332 = t280 ^ 2 + t274;
t138 = t332 * t221;
t212 = rSges(6,1) * t253 - rSges(6,2) * t284;
t329 = pkin(4) * t360;
t440 = pkin(3) - t269;
t338 = -t276 * t440 + t212 - t329;
t141 = t338 * t279;
t142 = t338 * t280;
t421 = m(6) / 0.2e1;
t422 = m(5) / 0.2e1;
t375 = t138 * t422 + (t141 * t279 + t142 * t280) * t421;
t441 = pkin(4) * t257;
t139 = -t279 * t441 - t181;
t295 = (pkin(4) * t326 - t329) * t280;
t140 = t295 + t182;
t242 = t278 * t359 - t279 * t326;
t243 = t292 * t279;
t202 = -rSges(5,1) * t242 - rSges(5,2) * t243;
t244 = t257 * t280;
t245 = t292 * t280;
t203 = rSges(5,1) * t244 - rSges(5,2) * t245;
t296 = -t202 * t279 - t203 * t280;
t376 = t296 * t422 + (t139 * t279 - t140 * t280) * t421;
t26 = t376 - t375;
t457 = t26 * qJD(1);
t361 = t276 * t280;
t372 = Icges(6,4) * t234;
t166 = Icges(6,2) * t233 - Icges(6,6) * t279 + t372;
t226 = Icges(6,4) * t233;
t169 = Icges(6,1) * t234 - Icges(6,5) * t279 + t226;
t312 = -t233 * t166 - t234 * t169 + t279 * (Icges(6,5) * t234 + Icges(6,6) * t233 - Icges(6,3) * t279);
t456 = t212 * t332;
t241 = Icges(6,4) * t284;
t210 = Icges(6,1) * t253 - t241;
t341 = -Icges(6,2) * t253 + t210 - t241;
t273 = t279 * pkin(6);
t293 = rSges(5,1) * t245 + rSges(5,2) * t244 - t279 * rSges(5,3);
t430 = (pkin(2) + pkin(3)) * t277 + t320;
t126 = t279 * qJ(2) + t280 * t430 - t273 + t293;
t432 = -t243 * rSges(5,1) + t242 * rSges(5,2);
t443 = t270 + (-rSges(5,3) - pkin(6)) * t280 - t430 * t279 + t432;
t78 = t126 * t279 + t280 * t443;
t225 = Icges(6,4) * t232;
t164 = Icges(6,2) * t231 + Icges(6,6) * t280 + t225;
t224 = Icges(6,4) * t231;
t168 = -Icges(6,1) * t232 - Icges(6,5) * t280 - t224;
t353 = t233 * t164 - t168 * t234;
t237 = Icges(5,4) * t243;
t187 = -Icges(5,2) * t242 + Icges(5,6) * t280 + t237;
t236 = Icges(5,4) * t242;
t191 = -Icges(5,1) * t243 - Icges(5,5) * t280 + t236;
t453 = -t187 * t242 - t243 * t191;
t351 = t244 * t187 - t191 * t245;
t401 = -t279 / 0.2e1;
t451 = -t280 / 0.2e1;
t399 = t280 / 0.2e1;
t186 = Icges(5,5) * t245 + Icges(5,6) * t244 - Icges(5,3) * t279;
t374 = Icges(5,4) * t245;
t189 = Icges(5,2) * t244 - Icges(5,6) * t279 + t374;
t238 = Icges(5,4) * t244;
t192 = Icges(5,1) * t245 - Icges(5,5) * t279 + t238;
t350 = t244 * t189 + t245 * t192;
t311 = t279 * t186 - t350;
t184 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t280;
t71 = t184 * t280 + t453;
t449 = t311 - t71;
t362 = t276 * t279;
t373 = Icges(5,4) * t257;
t217 = -Icges(5,2) * t292 + t373;
t218 = -Icges(5,1) * t292 - t373;
t445 = (t218 / 0.2e1 - t217 / 0.2e1) * t257;
t170 = t280 * rSges(6,3) - t433;
t106 = -t279 * t181 - t280 * t182;
t211 = -rSges(6,1) * t284 - rSges(6,2) * t253;
t368 = t211 * t280;
t369 = t211 * t279;
t348 = Icges(6,1) * t233 - t166 - t372;
t349 = Icges(6,1) * t231 - t164 - t225;
t286 = -t348 * t279 + t349 * t280;
t300 = (Icges(6,5) * t231 - Icges(6,6) * t232) * t280 - (Icges(6,5) * t233 - Icges(6,6) * t234) * t279;
t346 = -Icges(6,2) * t234 + t169 + t226;
t347 = -Icges(6,2) * t232 - t168 + t224;
t428 = t346 * t279 - t347 * t280;
t388 = (-t233 * t428 + t286 * t234 - t300 * t279) * t401 + (-t231 * t428 + t286 * t232 + t300 * t280) * t399;
t66 = (pkin(3) * t359 + t280 * pkin(6) - t170 + t333) * t279 + (-t172 - t273 + (t277 * t440 - t328) * t280) * t280;
t6 = t388 + m(6) * (t66 * t106 + t141 * t369 + t142 * t368);
t438 = t6 * qJD(5);
t436 = t332 * t276;
t434 = -t242 * t189 + t243 * t192;
t254 = Icges(5,4) * t292;
t216 = -Icges(5,2) * t257 - t254;
t219 = Icges(5,1) * t257 - t254;
t337 = t216 + t219;
t297 = -t202 * t280 + t203 * t279;
t302 = t139 * t280 + t140 * t279;
t377 = (m(5) * t297 + m(6) * t302) * t276 / 0.2e1;
t342 = -Icges(5,2) * t245 + t192 + t238;
t343 = -Icges(5,2) * t243 - t191 - t236;
t429 = t279 * t342 - t280 * t343;
t427 = (rSges(4,3) + qJ(3)) * t276 + (rSges(4,1) + pkin(2)) * t277 + pkin(1);
t371 = Icges(6,4) * t253;
t208 = -Icges(6,2) * t284 + t371;
t209 = -Icges(6,1) * t284 - t371;
t315 = -t341 * t284 / 0.2e1 + (t209 / 0.2e1 - t208 / 0.2e1) * t253;
t161 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t280;
t319 = (t353 * t280 + (-t161 * t280 + t312) * t279) * t399 + (t279 * t312 + t280 * (-t161 * t279 + t353)) * t451;
t426 = -0.2e1 * t436;
t424 = 0.4e1 * qJD(1);
t423 = 2 * qJD(4);
t419 = m(5) * (t126 * t203 - t202 * t443);
t418 = m(5) * (t126 * t361 - t362 * t443);
t417 = m(5) * t78;
t387 = t110 * t369 + t368 * t444;
t415 = m(6) * (t141 * t182 - t142 * t181 + t387);
t413 = m(6) * (t302 * t212 + t387);
t411 = m(6) * (t110 * t140 + t139 * t444);
t410 = m(6) * (t110 * t361 - t362 * t444);
t409 = m(6) * (t110 * t279 + t280 * t444);
t407 = m(6) * (t141 * t361 - t142 * t362);
t400 = t279 / 0.2e1;
t396 = m(3) * ((rSges(3,2) * t362 + rSges(3,3) * t280 + t270) * t280 + (-rSges(3,2) * t361 + (rSges(3,3) + qJ(2)) * t279) * t279);
t194 = rSges(4,2) * t280 - t279 * t427 + t270;
t195 = (rSges(4,2) + qJ(2)) * t279 + t427 * t280;
t395 = m(4) * (-t194 * t362 + t195 * t361);
t394 = m(4) * (t194 * t280 + t195 * t279);
t386 = m(6) * qJD(4);
t299 = -t181 * t280 + t182 * t279;
t291 = t299 * t276;
t56 = -t291 * m(6) / 0.2e1;
t358 = t56 * qJD(3);
t288 = t106 * t421;
t314 = m(6) * t456;
t70 = t288 - t314 / 0.2e1;
t357 = t70 * qJD(1);
t345 = -Icges(5,1) * t242 - t187 - t237;
t344 = Icges(5,1) * t244 - t189 - t374;
t340 = -t208 + t209;
t339 = -t292 * pkin(4) + t211;
t336 = -t217 + t218;
t105 = (t421 + t422 + m(4) / 0.2e1) * t426;
t331 = t105 * qJD(1);
t150 = t339 * t279;
t151 = t339 * t280;
t301 = t150 * t279 + t151 * t280;
t298 = t280 * (-Icges(5,5) * t242 - Icges(5,6) * t243) - t279 * (Icges(5,5) * t244 - Icges(5,6) * t245);
t205 = -Icges(6,5) * t284 - Icges(6,6) * t253;
t290 = -t319 + (-t205 * t279 + t341 * t233 + t340 * t234 + t348 * t253 - t284 * t346) * t401 + (t205 * t280 + t341 * t231 + t340 * t232 + t349 * t253 - t284 * t347) * t399;
t289 = -t315 + (t400 + t401) * ((Icges(6,5) * t253 - Icges(6,6) * t284) * t280 + t208 * t231 + t210 * t232);
t285 = -t344 * t279 + t345 * t280;
t220 = -rSges(5,1) * t292 - rSges(5,2) * t257;
t214 = -Icges(5,5) * t292 - Icges(5,6) * t257;
t104 = (m(6) / 0.4e1 + m(5) / 0.4e1 + m(4) / 0.4e1) * t426 + (m(6) + m(5) + m(4)) * t436 / 0.2e1;
t96 = -t279 * t170 - t280 * t172;
t81 = -t274 * t441 - t280 * t295 + t106;
t80 = t407 / 0.2e1;
t69 = t288 + t314 / 0.2e1;
t68 = -t106 * t277 + t211 * t436;
t65 = m(6) * t68 * qJD(5);
t57 = t291 * t421;
t45 = t279 * t311 + t280 * (-t184 * t279 + t351);
t44 = -t279 * (t280 * t186 + t434) + t280 * t71;
t36 = t395 + t410 + t418;
t34 = t413 / 0.2e1;
t31 = t415 / 0.2e1;
t30 = t315 + t458;
t29 = t394 + t396 + t409 + t417;
t27 = t375 + t376;
t20 = t80 - t377;
t19 = -t407 / 0.2e1 + t377;
t18 = t80 + t377;
t17 = t445 - (t219 / 0.2e1 + t216 / 0.2e1) * t292 + t419 + t411 + t315;
t16 = t351 * t280 + (t449 + t453) * t279;
t15 = (t350 + t449) * t280 + t279 * t434;
t8 = m(6) * (t96 * t106 + t211 * t456) + t388;
t7 = t8 * qJD(5);
t4 = t31 - t413 / 0.2e1 + t319;
t3 = t34 - t415 / 0.2e1 + t319;
t2 = t31 + t34 + t290;
t1 = (-t45 / 0.2e1 + t16 / 0.2e1) * t280 + (-t15 / 0.2e1 - t44 / 0.2e1) * t279 + t319;
t5 = [t29 * qJD(2) + t36 * qJD(3) + t17 * qJD(4) + t30 * qJD(5), qJD(1) * t29 + qJD(3) * t104 + qJD(4) * t27 + qJD(5) * t69, qJD(1) * t36 + qJD(2) * t104 + qJD(4) * t18 - qJD(5) * t56, t17 * qJD(1) + t27 * qJD(2) + t18 * qJD(3) + t2 * qJD(5) + ((t110 * t150 + t139 * t142 + t140 * t141 + t151 * t444) * t421 + (t220 * t78 + t297 * t221) * t422) * t423 + (t16 * t451 + t290 + (t214 * t280 - t242 * t337 + t243 * t336 + t257 * t345 - t292 * t343 + t45) * t399 + (-t214 * t279 + t244 * t337 + t245 * t336 + t257 * t344 - t292 * t342) * t401 + (t15 + t44) * t400) * qJD(4), t30 * qJD(1) + t69 * qJD(2) - t358 + t2 * qJD(4) + ((t212 * t299 + t387) * m(6) + t290) * qJD(5); t105 * qJD(3) + t26 * qJD(4) + t70 * qJD(5) + (-t409 / 0.4e1 - t417 / 0.4e1 - t394 / 0.4e1 - t396 / 0.4e1) * t424, 0, t331, t457 + (-t150 * t280 + t151 * t279) * t386, t357; -t105 * qJD(2) + t19 * qJD(4) + t57 * qJD(5) + (-t410 / 0.4e1 - t418 / 0.4e1 - t395 / 0.4e1) * t424, -t331, 0, t19 * qJD(1) + ((t220 * t436 - t277 * t296) * t422 + (t276 * t301 - t277 * t81) * t421) * t423 + t65, t57 * qJD(1) + t386 * t68 + t65; -t26 * qJD(2) + t20 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t411 / 0.4e1 - t419 / 0.4e1) * t424 + (t289 + t337 * t292 / 0.2e1 - t445) * qJD(1), -t457, t20 * qJD(1), t1 * qJD(1) + (m(5) * ((-t279 * (rSges(5,3) * t280 - t432) - t280 * t293) * t296 + t138 * t220) + (-t244 * t429 + t285 * t245 - t298 * t279) * t401 + (t242 * t429 + t285 * t243 + t298 * t280) * t399 + m(6) * (t141 * t150 + t142 * t151 + t66 * t81) + t388) * qJD(4) + t438, t4 * qJD(1) + t6 * qJD(4) + t438; (t289 - t458) * qJD(1) - t70 * qJD(2) + t358 + t3 * qJD(4) + t319 * qJD(5), -t357, t56 * qJD(1), t3 * qJD(1) + ((t212 * t301 + t81 * t96) * m(6) + t388) * qJD(4) + t7, qJD(1) * t319 + qJD(4) * t8 + t7;];
Cq = t5;
