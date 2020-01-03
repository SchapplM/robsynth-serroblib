% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:33
% DurationCPUTime: 5.50s
% Computational Cost: add. (13683->368), mult. (13477->499), div. (0->0), fcn. (12476->8), ass. (0->236)
t280 = cos(qJ(1));
t273 = qJ(2) + pkin(7);
t257 = qJ(4) + t273;
t251 = sin(t257);
t252 = cos(t257);
t202 = rSges(5,1) * t251 + rSges(5,2) * t252;
t255 = sin(t273);
t277 = sin(qJ(2));
t402 = pkin(2) * t277;
t222 = -pkin(3) * t255 - t402;
t289 = t202 - t222;
t113 = t289 * t280;
t278 = sin(qJ(1));
t440 = t289 * t278;
t448 = t278 * t440;
t439 = t113 * t280 + t448;
t256 = cos(t273);
t283 = rSges(4,1) * t255 + rSges(4,2) * t256 + t402;
t441 = t283 * t280;
t442 = t283 * t278;
t446 = t278 * t442 + t280 * t441;
t395 = -m(5) * t439 / 0.2e1 - m(4) * t446 / 0.2e1;
t183 = t202 * t280;
t213 = t280 * t222;
t116 = t213 - t183;
t427 = m(5) / 0.2e1;
t428 = m(4) / 0.2e1;
t396 = (-t116 * t280 + t448) * t427 + t446 * t428;
t23 = t396 - t395;
t454 = t23 * qJD(1);
t274 = t278 ^ 2;
t275 = t280 ^ 2;
t324 = t274 + t275;
t453 = t202 * t324;
t242 = Icges(5,4) * t252;
t199 = -Icges(5,2) * t251 + t242;
t200 = Icges(5,1) * t251 + t242;
t452 = t199 + t200;
t279 = cos(qJ(2));
t451 = -Icges(3,5) * t277 - Icges(4,5) * t255 - Icges(3,6) * t279 - Icges(4,6) * t256;
t406 = t278 / 0.2e1;
t404 = -t280 / 0.2e1;
t449 = t280 / 0.2e1;
t359 = t255 * t278;
t357 = t256 * t278;
t354 = t277 * t278;
t352 = t278 * t279;
t382 = Icges(4,4) * t255;
t215 = Icges(4,2) * t256 + t382;
t218 = Icges(4,1) * t256 - t382;
t383 = Icges(3,4) * t277;
t229 = Icges(3,2) * t279 + t383;
t232 = Icges(3,1) * t279 - t383;
t447 = -(t218 / 0.2e1 - t215 / 0.2e1) * t255 - (t232 / 0.2e1 - t229 / 0.2e1) * t277;
t392 = rSges(5,1) * t252;
t203 = -rSges(5,2) * t251 + t392;
t294 = Icges(5,5) * t251 + Icges(5,6) * t252;
t176 = t294 * t278;
t177 = t280 * t294;
t381 = Icges(5,4) * t251;
t201 = Icges(5,1) * t252 - t381;
t145 = Icges(5,5) * t278 + t201 * t280;
t198 = Icges(5,2) * t252 + t381;
t338 = -t198 * t280 + t145;
t365 = t251 * t278;
t225 = Icges(5,4) * t365;
t361 = t252 * t278;
t144 = Icges(5,1) * t361 - Icges(5,5) * t280 - t225;
t339 = -Icges(5,2) * t361 + t144 - t225;
t143 = Icges(5,6) * t278 + t199 * t280;
t340 = -t200 * t280 - t143;
t142 = Icges(5,4) * t361 - Icges(5,2) * t365 - Icges(5,6) * t280;
t341 = t200 * t278 + t142;
t431 = (-t338 * t278 + t339 * t280) * t251 + (t340 * t278 + t341 * t280) * t252;
t398 = (-t274 * t177 + (t278 * t176 + t431) * t280) * t406 + (-t275 * t176 + (t280 * t177 + t431) * t278) * t404;
t271 = t279 * pkin(2);
t254 = t271 + pkin(1);
t401 = pkin(3) * t256;
t221 = t254 + t401;
t276 = -qJ(3) - pkin(5);
t323 = -pkin(6) + t276;
t248 = t278 * t323;
t253 = t278 * t276;
t327 = -t278 * t221 - t280 * t323;
t272 = t280 * pkin(5);
t355 = t276 * t280;
t399 = pkin(1) - t254;
t342 = -t278 * (t399 * t278 - t272 - t355) + t280 * (-t278 * pkin(5) - t399 * t280 - t253);
t147 = rSges(5,1) * t361 - rSges(5,2) * t365 - t280 * rSges(5,3);
t364 = t251 * t280;
t317 = -rSges(5,2) * t364 + t278 * rSges(5,3);
t360 = t252 * t280;
t93 = t278 * t147 + t280 * (rSges(5,1) * t360 + t317);
t43 = -t278 * (t278 * t254 + t327 + t355) + (-t248 + t253 + (t221 - t254) * t280) * t280 + t342 + t93;
t182 = t202 * t278;
t438 = t278 * t182 + t280 * t183;
t6 = t398 + m(5) * (t439 * t203 - t43 * t438);
t444 = t6 * qJD(4);
t100 = -t147 + t327;
t101 = -t248 + (t221 + t392) * t280 + t317;
t443 = t100 * t280 + t101 * t278;
t437 = t451 * t278;
t436 = t451 * t280;
t250 = Icges(4,4) * t256;
t216 = -Icges(4,2) * t255 + t250;
t217 = Icges(4,1) * t255 + t250;
t267 = Icges(3,4) * t279;
t230 = -Icges(3,2) * t277 + t267;
t231 = Icges(3,1) * t277 + t267;
t304 = t452 * t252 / 0.2e1 + (-t198 / 0.2e1 + t201 / 0.2e1) * t251;
t119 = t145 * t361;
t140 = Icges(5,5) * t361 - Icges(5,6) * t365 - Icges(5,3) * t280;
t197 = Icges(5,5) * t252 - Icges(5,6) * t251;
t369 = t197 * t280;
t141 = Icges(5,3) * t278 + t369;
t308 = t143 * t251 - t140;
t310 = t141 * t280 - t119;
t346 = t278 * t141 + t145 * t360;
t347 = -t278 * t140 - t144 * t360;
t374 = t142 * t251;
t62 = -t143 * t365 - t310;
t63 = -t142 * t364 - t347;
t64 = -t143 * t364 + t346;
t320 = ((t62 - t119 + (t141 + t374) * t280 + t347) * t280 + t346 * t278) * t404 + (t64 * t278 - t280 * t63) * t449 + ((t308 * t278 + t310 + t62 + t63) * t278 + (-t346 + t64 + (-t144 * t252 + t374) * t278 + (t308 + t140) * t280) * t280) * t406;
t189 = Icges(3,5) * t278 + t232 * t280;
t328 = -t229 * t280 + t189;
t187 = Icges(3,6) * t278 + t230 * t280;
t330 = -t231 * t280 - t187;
t172 = Icges(4,5) * t278 + t218 * t280;
t332 = -t215 * t280 + t172;
t170 = Icges(4,6) * t278 + t216 * t280;
t334 = -t217 * t280 - t170;
t433 = -t328 * t354 + t330 * t352 - t332 * t359 + t334 * t357;
t245 = Icges(3,4) * t354;
t188 = Icges(3,1) * t352 - Icges(3,5) * t280 - t245;
t329 = -Icges(3,2) * t352 + t188 - t245;
t186 = Icges(3,4) * t352 - Icges(3,2) * t354 - Icges(3,6) * t280;
t331 = t231 * t278 + t186;
t239 = Icges(4,4) * t359;
t171 = Icges(4,1) * t357 - Icges(4,5) * t280 - t239;
t333 = -Icges(4,2) * t357 + t171 - t239;
t169 = Icges(4,4) * t357 - Icges(4,2) * t359 - Icges(4,6) * t280;
t335 = t217 * t278 + t169;
t432 = t333 * t255 + t335 * t256 + t329 * t277 + t331 * t279;
t430 = 0.4e1 * qJD(1);
t429 = 2 * qJD(2);
t394 = rSges(3,1) * t279;
t322 = pkin(1) + t394;
t325 = rSges(3,2) * t354 + t280 * rSges(3,3);
t136 = -t322 * t278 + t272 + t325;
t353 = t277 * t280;
t247 = rSges(3,2) * t353;
t137 = -t247 + t322 * t280 + (rSges(3,3) + pkin(5)) * t278;
t233 = rSges(3,1) * t277 + rSges(3,2) * t279;
t211 = t233 * t278;
t212 = t233 * t280;
t426 = m(3) * (t136 * t211 - t137 * t212);
t393 = rSges(4,1) * t256;
t319 = t254 + t393;
t326 = rSges(4,2) * t359 + t280 * rSges(4,3);
t108 = -t319 * t278 + t326 - t355;
t358 = t255 * t280;
t318 = -rSges(4,2) * t358 + t278 * rSges(4,3);
t109 = t319 * t280 - t253 + t318;
t424 = m(4) * (t108 * t442 - t109 * t441);
t423 = m(4) * (t108 * t280 + t109 * t278);
t285 = t443 * t203;
t416 = m(5) * (-t113 * t182 + t183 * t440 - t285);
t415 = m(5) * (-t285 + (-t116 * t278 - t280 * t440) * t202);
t414 = m(5) * (t100 * t440 + t101 * t116);
t413 = m(5) * (t100 * t182 - t101 * t183);
t412 = m(5) * t443;
t407 = -t278 / 0.2e1;
t372 = t169 * t255;
t370 = t186 * t277;
t356 = t256 * t280;
t351 = t279 * t280;
t286 = t438 * t427;
t303 = m(5) * t453;
t66 = t286 + t303 / 0.2e1;
t348 = t66 * qJD(1);
t167 = Icges(4,5) * t357 - Icges(4,6) * t359 - Icges(4,3) * t280;
t345 = -t278 * t167 - t171 * t356;
t296 = Icges(4,5) * t256 - Icges(4,6) * t255;
t168 = Icges(4,3) * t278 + t296 * t280;
t344 = t278 * t168 + t172 * t356;
t184 = Icges(3,5) * t352 - Icges(3,6) * t354 - Icges(3,3) * t280;
t337 = -t278 * t184 - t188 * t351;
t298 = Icges(3,5) * t279 - Icges(3,6) * t277;
t185 = Icges(3,3) * t278 + t298 * t280;
t336 = t278 * t185 + t189 * t351;
t321 = rSges(4,2) * t255 - t271 - t393;
t122 = t172 * t357;
t309 = t168 * t280 - t122;
t149 = t189 * t352;
t307 = t185 * t280 - t149;
t306 = t170 * t255 - t167;
t305 = t187 * t277 - t184;
t288 = -t203 - t271 - t401;
t281 = (-t198 + t201) * t252 - t452 * t251;
t284 = -t320 + (t278 * t197 + t340 * t251 + t338 * t252 + t280 * t281) * t406 + (-t341 * t251 + t339 * t252 + t281 * t278 - t369) * t404;
t282 = -t304 + (t406 + t407) * (t142 * t252 + t144 * t251);
t235 = -rSges(3,2) * t277 + t394;
t165 = t321 * t280;
t163 = t321 * t278;
t114 = t288 * t280;
t112 = t288 * t278;
t84 = -t187 * t353 + t336;
t83 = -t186 * t353 - t337;
t82 = -t187 * t354 - t307;
t70 = -t170 * t358 + t344;
t69 = -t169 * t358 - t345;
t68 = -t170 * t359 - t309;
t65 = t286 - t303 / 0.2e1;
t57 = t280 * t213 + t222 * t274 - t438;
t50 = t84 * t278 - t280 * t83;
t49 = t82 * t278 - t280 * (-(-t188 * t279 + t370) * t278 - t184 * t280);
t46 = t70 * t278 - t280 * t69;
t45 = t68 * t278 - t280 * (-(-t171 * t256 + t372) * t278 - t167 * t280);
t42 = t412 + t423;
t35 = t304 + t413;
t32 = t415 / 0.2e1;
t30 = t416 / 0.2e1;
t25 = t395 + t396;
t22 = (t82 - t149 + (t185 + t370) * t280 + t337) * t280 + t336 * t278;
t21 = (t305 * t280 - t336 + t84) * t280 + (t305 * t278 + t307 + t83) * t278;
t19 = (t68 - t122 + (t168 + t372) * t280 + t345) * t280 + t344 * t278;
t18 = (t306 * t280 - t344 + t70) * t280 + (t306 * t278 + t309 + t69) * t278;
t9 = (t231 / 0.2e1 + t230 / 0.2e1) * t279 + (t217 / 0.2e1 + t216 / 0.2e1) * t256 + t426 + t424 + t414 + t304 - t447;
t8 = m(5) * (t203 * t453 - t438 * t93) + t398;
t7 = t8 * qJD(4);
t4 = t30 - t415 / 0.2e1 + t320;
t3 = t32 - t416 / 0.2e1 + t320;
t2 = t30 + t32 + t284;
t1 = (t50 / 0.2e1 - t19 / 0.2e1 - t22 / 0.2e1 + t46 / 0.2e1) * t280 + (t21 / 0.2e1 + t45 / 0.2e1 + t49 / 0.2e1 + t18 / 0.2e1) * t278 + t320;
t5 = [qJD(2) * t9 + qJD(3) * t42 + qJD(4) * t35, t9 * qJD(1) + t25 * qJD(3) + t2 * qJD(4) + (m(3) * ((-t136 * t280 - t137 * t278) * t235 + (-t211 * t280 + t212 * t278) * t233) / 0.2e1 + (t108 * t165 + t109 * t163) * t428 + (t100 * t114 + t101 * t112 + (-t113 - t116) * t440) * t427) * t429 + ((t298 + t296) * (t274 / 0.2e1 + t275 / 0.2e1) + t284 + (t334 * t255 + t332 * t256 + t330 * t277 + t328 * t279) * t406 + (t19 + t22) * t449 + (t21 + t45 + t49 + t18) * t407 + (-t335 * t255 + t333 * t256 - t331 * t277 + t329 * t279 + t46 + t50) * t404) * qJD(2), qJD(1) * t42 + qJD(2) * t25 + qJD(4) * t65, t35 * qJD(1) + t2 * qJD(2) + t65 * qJD(3) + ((-t285 + (-t182 * t280 + t183 * t278) * t202) * m(5) + t284) * qJD(4); (t282 - (t217 + t216) * t256 / 0.2e1 - (t231 + t230) * t279 / 0.2e1 + t447) * qJD(1) + t1 * qJD(2) - t23 * qJD(3) + t4 * qJD(4) + (-t414 / 0.4e1 - t424 / 0.4e1 - t426 / 0.4e1) * t430, t1 * qJD(1) + (m(5) * (-t112 * t440 - t113 * t114 + t43 * t57) + m(4) * (-t442 * t163 - t441 * t165 - (t278 * (rSges(4,1) * t357 - t326) + t280 * (rSges(4,1) * t356 + t318) + t342) * t283 * t324) + m(3) * ((t278 * (rSges(3,1) * t352 - t325) + t280 * (rSges(3,1) * t351 + t278 * rSges(3,3) - t247)) * (-t278 * t211 - t212 * t280) + t324 * t235 * t233) + t398 + ((-t437 * t278 + t432 * t280 + t433) * t280 + t436 * t274) * t406 + (((t432 - t436) * t280 + t433) * t278 + t437 * t275) * t404) * qJD(2) + t444, -t454, t4 * qJD(1) + t6 * qJD(2) + t444; t23 * qJD(2) + t66 * qJD(4) + (-t423 / 0.4e1 - t412 / 0.4e1) * t430, t454 + ((-t112 * t280 + t278 * t114) * t427 + (-t163 * t280 + t278 * t165) * t428) * t429, 0, t348; (t282 - t413) * qJD(1) + t3 * qJD(2) - t66 * qJD(3) + t320 * qJD(4), t3 * qJD(1) + ((t93 * t57 + (-t112 * t278 - t114 * t280) * t202) * m(5) + t398) * qJD(2) + t7, -t348, qJD(1) * t320 + qJD(2) * t8 + t7;];
Cq = t5;
