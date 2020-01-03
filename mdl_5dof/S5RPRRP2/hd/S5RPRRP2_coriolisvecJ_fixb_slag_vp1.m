% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:24
% DurationCPUTime: 12.05s
% Computational Cost: add. (12066->423), mult. (8983->520), div. (0->0), fcn. (6890->8), ass. (0->253)
t240 = qJ(1) + pkin(8);
t233 = qJ(3) + t240;
t226 = sin(t233);
t227 = cos(t233);
t244 = cos(qJ(4));
t386 = t227 * t244;
t242 = sin(qJ(4));
t387 = t227 * t242;
t120 = Icges(6,4) * t386 - Icges(6,2) * t387 + Icges(6,6) * t226;
t122 = Icges(5,4) * t386 - Icges(5,2) * t387 + Icges(5,6) * t226;
t511 = t120 + t122;
t409 = Icges(6,4) * t242;
t200 = Icges(6,2) * t244 + t409;
t410 = Icges(5,4) * t242;
t202 = Icges(5,2) * t244 + t410;
t509 = t200 + t202;
t520 = Icges(5,5) + Icges(6,5);
t519 = -Icges(5,6) - Icges(6,6);
t518 = Icges(5,3) + Icges(6,3);
t190 = Icges(6,4) * t387;
t124 = Icges(6,1) * t386 + Icges(6,5) * t226 - t190;
t191 = Icges(5,4) * t387;
t126 = Icges(5,1) * t386 + Icges(5,5) * t226 - t191;
t517 = t124 + t126;
t234 = Icges(6,4) * t244;
t304 = -Icges(6,2) * t242 + t234;
t119 = -Icges(6,6) * t227 + t226 * t304;
t235 = Icges(5,4) * t244;
t305 = -Icges(5,2) * t242 + t235;
t121 = -Icges(5,6) * t227 + t226 * t305;
t512 = t119 + t121;
t205 = Icges(6,1) * t244 - t409;
t282 = t205 * t226;
t123 = -Icges(6,5) * t227 + t282;
t207 = Icges(5,1) * t244 - t410;
t283 = t207 * t226;
t125 = -Icges(5,5) * t227 + t283;
t510 = t123 + t125;
t516 = t205 + t207;
t451 = Icges(5,1) * t242 + t235;
t515 = -t509 * t242 + t244 * t451;
t514 = t304 + t305;
t513 = t511 * t242;
t239 = qJD(1) + qJD(3);
t383 = t239 * t244;
t199 = Icges(5,5) * t244 - Icges(5,6) * t242;
t197 = Icges(6,5) * t244 - Icges(6,6) * t242;
t278 = t197 * t226;
t466 = t199 * t226 - t518 * t227 + t278;
t465 = t518 * t226 + t520 * t386 + t519 * t387;
t452 = Icges(6,1) * t242 + t234;
t502 = t244 * t452 + t515;
t508 = -t517 * t244 + t513;
t507 = t451 + t452;
t389 = t226 * t244;
t468 = t510 * t389;
t506 = t517 * t389;
t467 = t512 * t387;
t505 = t514 * qJD(4);
t504 = t516 * qJD(4);
t196 = Icges(6,5) * t242 + Icges(6,6) * t244;
t198 = Icges(5,5) * t242 + Icges(5,6) * t244;
t503 = t196 + t198;
t501 = (t507 * qJD(4) - t520 * t239) * t242 + (t509 * qJD(4) + t519 * t239) * t244;
t390 = t226 * t242;
t475 = -t466 * t227 - t512 * t390 + t468;
t474 = -t465 * t227 - t511 * t390 + t506;
t473 = -t466 * t226 - t510 * t386 + t467;
t472 = -t465 * t226 + t508 * t227;
t138 = t196 * t227;
t140 = t198 * t227;
t500 = t502 * t226 - t138 - t140;
t392 = t198 * t226;
t394 = t196 * t226;
t499 = t515 * t227 + t386 * t452 + t392 + t394;
t241 = -qJ(5) - pkin(7);
t498 = rSges(6,3) - t241;
t384 = t239 * t242;
t195 = rSges(5,2) * t387;
t130 = rSges(5,1) * t386 + rSges(5,3) * t226 - t195;
t114 = t239 * t130;
t426 = pkin(3) * t227;
t159 = pkin(7) * t226 + t426;
t153 = t239 * t159;
t388 = t227 * t239;
t391 = t226 * t239;
t362 = pkin(3) * t388 + pkin(7) * t391;
t343 = t227 * t383;
t365 = rSges(5,1) * t343 + rSges(5,3) * t391;
t497 = t362 + t365 - t114 - t153;
t425 = pkin(4) * t244;
t228 = pkin(3) + t425;
t340 = rSges(6,1) * t343 + rSges(6,3) * t391 + t228 * t388;
t363 = -rSges(6,2) * t387 + t227 * t228;
t378 = -rSges(6,1) * t386 - t363 + t426 + (pkin(7) - t498) * t226;
t486 = t378 * t239;
t496 = t340 - t153 + t486;
t460 = t510 * t244;
t459 = t512 * t242;
t493 = t514 * t383;
t492 = t502 * t239 + (-t197 - t199) * qJD(4);
t491 = t459 - t460;
t489 = -t504 * t244 + t505 * t242 - t503 * t239 + (t507 * t242 + t509 * t244) * qJD(4);
t488 = t503 * qJD(4) - t518 * t239;
t487 = t500 * t239;
t279 = t199 * t239;
t485 = (t197 * t388 - t488 * t226 + t227 * t279 + t491 * t239) * t227;
t484 = (t472 * t226 - t473 * t227) * qJD(4);
t483 = (t474 * t226 - t475 * t227) * qJD(4);
t454 = t226 * rSges(4,1) + t227 * rSges(4,2);
t132 = t454 * t239;
t231 = sin(t240);
t243 = sin(qJ(1));
t237 = t243 * pkin(1);
t453 = pkin(2) * t231 + t237;
t289 = t453 * qJD(1);
t111 = t289 + t132;
t482 = t499 * t239;
t481 = t483 + t487;
t480 = -t482 + t484;
t479 = -t516 * t242 * t388 + t491 * qJD(4) + t501 * t226 - t227 * t493;
t478 = -t508 * qJD(4) + (-t282 - t283) * t384 - t226 * t493 - t501 * t227;
t477 = t492 * t226 + t489 * t227;
t476 = -t489 * t226 + t492 * t227;
t471 = t510 * t242 + t512 * t244;
t359 = t452 + t304;
t360 = t200 - t205;
t470 = (t242 * t359 + t244 * t360) * t239;
t357 = t451 + t305;
t358 = t202 - t207;
t469 = (t242 * t357 + t244 * t358) * t239;
t463 = t226 * t279 + (t278 - t508) * t239 + t488 * t227;
t462 = t517 * t242 + t511 * t244;
t456 = 0.2e1 * qJD(4);
t208 = t227 * t241;
t337 = rSges(6,1) * t389 + t226 * t228 + t208;
t215 = qJD(5) * t226;
t344 = t226 * t384;
t339 = rSges(6,2) * t344 + rSges(6,3) * t388 + t215;
t221 = t226 * pkin(3);
t158 = -pkin(7) * t227 + t221;
t379 = -rSges(6,2) * t390 - rSges(6,3) * t227 - t158 + t337;
t445 = -t239 * (t158 + t379) + t215;
t193 = rSges(5,1) * t389;
t128 = -rSges(5,2) * t390 - rSges(5,3) * t227 + t193;
t210 = rSges(5,1) * t242 + rSges(5,2) * t244;
t348 = qJD(4) * t227;
t335 = t210 * t348;
t444 = -t239 * (t128 + t158) - t335;
t247 = qJD(1) ^ 2;
t428 = t239 / 0.2e1;
t427 = rSges(5,3) + pkin(7);
t232 = cos(t240);
t225 = pkin(2) * t232;
t245 = cos(qJ(1));
t238 = t245 * pkin(1);
t184 = pkin(7) * t388;
t347 = qJD(4) * t242;
t333 = t227 * t347;
t276 = t226 * t383 + t333;
t346 = qJD(4) * t244;
t332 = t227 * t346;
t424 = -pkin(4) * t333 - t184 - (t208 + (-pkin(3) + t228) * t226) * t239 - rSges(6,1) * t276 - rSges(6,2) * t332 + t339;
t275 = -t226 * t346 - t227 * t384;
t334 = t226 * t347;
t345 = qJD(5) * t227;
t385 = t239 * t241;
t423 = -t345 + (-pkin(4) * t347 - t385) * t226 - t362 - rSges(6,1) * t334 + rSges(6,2) * t275 + t340;
t422 = rSges(5,1) * t244;
t421 = rSges(6,1) * t244;
t420 = rSges(5,2) * t242;
t419 = rSges(6,2) * t242;
t418 = pkin(1) * qJD(1);
t349 = qJD(4) * t226;
t336 = t210 * t349;
t350 = qJD(1) * t232;
t355 = pkin(2) * t350 + t245 * t418;
t292 = -t336 + t355;
t60 = (t130 + t159) * t239 + t292;
t417 = t239 * t60;
t416 = t244 * rSges(6,2);
t393 = t197 * t239;
t377 = t226 * t452 + t119;
t376 = t227 * t452 + t120;
t375 = t226 * t451 + t121;
t374 = t227 * t451 + t122;
t373 = -t200 * t226 + t123;
t372 = -Icges(6,2) * t386 + t124 - t190;
t371 = -t202 * t226 + t125;
t370 = -Icges(5,2) * t386 + t126 - t191;
t367 = rSges(5,2) * t344 + rSges(5,3) * t388;
t361 = t193 + t221;
t229 = t247 * t238;
t356 = t247 * t225 + t229;
t352 = t225 + t238;
t351 = qJD(1) * t231;
t341 = t239 * t362 + t356;
t338 = t184 + t367;
t329 = pkin(3) + t422;
t327 = t349 / 0.2e1;
t326 = -t348 / 0.2e1;
t325 = t348 / 0.2e1;
t209 = rSges(6,1) * t242 + t416;
t323 = pkin(4) * t242 + t209;
t318 = -pkin(4) * t390 - t209 * t226;
t317 = -pkin(4) * t387 - t209 * t227;
t133 = rSges(4,1) * t388 - rSges(4,2) * t391;
t157 = rSges(4,1) * t227 - rSges(4,2) * t226;
t112 = t157 * t239 + t355;
t49 = t323 * t348 + t289 - t445;
t314 = t49 * t323;
t308 = qJD(4) * t323;
t257 = -t226 * t308 - t345;
t255 = t257 + t355;
t50 = (t159 - t378) * t239 + t255;
t313 = t50 * t323;
t161 = rSges(3,1) * t232 - rSges(3,2) * t231;
t310 = -t420 + t422;
t212 = -t419 + t421;
t59 = t289 - t444;
t309 = -t226 * t60 + t227 * t59;
t297 = t128 * t226 + t130 * t227;
t174 = t212 * qJD(4);
t291 = (t425 * qJD(4) + t174) * qJD(4);
t290 = t453 * t247;
t288 = qJD(4) * t210;
t277 = -pkin(2) * t351 - t243 * t418;
t268 = qJD(4) * (-t416 + (-rSges(6,1) - pkin(4)) * t242);
t261 = t242 * t373 + t244 * t377;
t260 = t242 * t372 + t244 * t376;
t259 = t242 * t371 + t244 * t375;
t258 = t242 * t370 + t244 * t374;
t91 = rSges(5,1) * t276 + rSges(5,2) * t332 - t367;
t93 = -rSges(5,1) * t334 + rSges(5,2) * t275 + t365;
t256 = (t128 * t239 - t91) * t227 + (t93 - t114) * t226;
t250 = ((((t466 - t508) * t227 - t468 + t472) * t227 + ((t466 - t513) * t226 + (t459 + t460) * t227 - t467 + t473 + t506) * t226) * qJD(4) + t487) * t327 + (t476 - t479) * t326 + ((((-t460 + t465) * t227 + t467 + t474) * t227 + ((t459 + t465) * t226 - t468 + t475) * t226) * qJD(4) + t480 + t482) * t325 - (t477 - t478 + t481) * t349 / 0.2e1 + (t499 * t227 + (t471 + t500) * t226) * qJD(4) * t428 + (t502 * qJD(4) + t504 * t242 + t505 * t244 - t462 * t326) * t239;
t134 = pkin(3) * t391 - t184;
t175 = t310 * qJD(4);
t47 = -t175 * t349 - t290 + (-t134 - t91 - t335) * t239;
t48 = t239 * t93 + (t175 * t227 - t210 * t391) * qJD(4) + t341;
t249 = (-t288 * t59 - t329 * t417 - t420 * t48 + t427 * t47) * t226 + (-rSges(5,2) * t384 * t59 - t288 * t60 + t329 * t47 - t427 * t48) * t227;
t23 = -t290 - t291 * t226 + (-t227 * t308 - t134 + t215 + t424) * t239;
t24 = t291 * t227 + (t257 + t423) * t239 + t341;
t248 = (t23 * t498 - t24 * t419 + (t50 * (-t228 - t421) - t49 * t241) * t239 + t49 * t268) * t226 + (t23 * t421 - t24 * rSges(6,3) + t49 * (-rSges(6,2) * t384 - qJD(5)) + (-t385 + t268) * t50) * t227;
t152 = t210 * t227;
t150 = t210 * t226;
t102 = t133 * t239 + t356;
t101 = -t132 * t239 - t290;
t65 = qJD(4) * t297 + qJD(2);
t42 = qJD(2) + (t379 * t226 - t378 * t227) * qJD(4);
t25 = t256 * qJD(4);
t5 = ((t379 * t239 + t424) * t227 + (t423 + t486) * t226) * qJD(4);
t1 = [t250 + m(4) * (t101 * (t157 + t352) + t102 * (t453 + t454) + (-t112 + t133 + t355) * t111) + m(3) * (-t229 + t247 * (t161 + t238) + (-0.2e1 * rSges(3,1) * t350 + 0.2e1 * rSges(3,2) * t351 + qJD(1) * t161) * qJD(1)) * (-t231 * rSges(3,1) - t232 * rSges(3,2) - t237) + (t23 * (t352 + t363) + t50 * (t277 + t339) + t24 * (t337 + t453) + t248 + (t355 + t50 - t255 + t496) * t49) * m(6) + (t47 * (-t195 + t352) + t60 * (t277 + t338) + t48 * (t453 + t361) + t249 + (t355 - t292 + t60 + t497) * t59) * m(5); m(5) * t25 + m(6) * t5; t250 + (t23 * t363 + t24 * t337 + t248 - (-t226 * t314 - t227 * t313) * qJD(4) + (t339 - t445) * t50 + (t345 + t496) * t49) * m(6) + (-t47 * t195 + t48 * t361 + t249 + (t338 - t444) * t60 + (t336 + t497) * t59) * m(5) + (t101 * t157 + t102 * t454 + t111 * t133 - t112 * t132 - (t111 * t157 - t112 * t454) * t239) * m(4); -(((t357 + t359) * t244 + (-t358 - t360) * t242) * t239 + (((-t371 - t373) * t227 + (t370 + t372) * t226) * t244 + ((t375 + t377) * t227 + (-t374 - t376) * t226) * t242) * qJD(4)) * t239 / 0.2e1 + ((t239 * t462 + t479) * t227 + (t239 * t471 + t478) * t226) * t428 + ((t140 * t349 - t279) * t226 + (t469 + (-t259 * t227 + (-t392 + t258) * t226) * qJD(4)) * t227 + (t138 * t349 - t393) * t226 + (t470 + (-t261 * t227 + (-t394 + t260) * t226) * qJD(4)) * t227) * t327 + ((-t348 * t392 - t279) * t227 + (-t469 + (-t258 * t226 + (t140 + t259) * t227) * qJD(4)) * t226 + (-t348 * t394 - t393) * t227 + (-t470 + (-t260 * t226 + (t138 + t261) * t227) * qJD(4)) * t226) * t325 + (t25 * t297 + t65 * t256 + t309 * t175 + ((t48 - t417) * t227 + (-t239 * t59 - t47) * t226) * t210 - (-t150 * t59 - t152 * t60) * t239 - (t65 * (-t150 * t226 - t152 * t227) + t309 * t310) * qJD(4)) * m(5) - (t477 * t239 + (t472 * t388 + (t463 * t226 + t473 * t239 + t485) * t226) * t456) * t226 / 0.2e1 - (t476 * t239 + ((t474 * t239 + t485) * t227 + (t227 * t463 + t475 * t239) * t226) * t456) * t227 / 0.2e1 + ((-t5 * t378 + t42 * t424 + t24 * t323 + t49 * t174 + (t379 * t42 - t313) * t239) * t227 + (t5 * t379 + t42 * t423 - t23 * t323 + t50 * (-pkin(4) * t346 - t174) + (t378 * t42 - t314) * t239) * t226 - (t317 * t50 + t318 * t49) * t239 - ((t49 * t212 + t317 * t42) * t227 + (t50 * (-t212 - t425) + t318 * t42) * t226) * qJD(4)) * m(6) + (t481 + t483) * t391 / 0.2e1 - (t480 + t484) * t388 / 0.2e1; m(6) * (-t226 * t24 - t227 * t23);];
tauc = t1(:);
