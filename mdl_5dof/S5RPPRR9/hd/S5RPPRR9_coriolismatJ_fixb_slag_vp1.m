% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR9_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:20
% EndTime: 2019-12-31 18:02:38
% DurationCPUTime: 14.32s
% Computational Cost: add. (31051->561), mult. (69639->833), div. (0->0), fcn. (88709->8), ass. (0->338)
t315 = sin(qJ(4));
t456 = sin(pkin(8));
t457 = cos(pkin(8));
t484 = sin(qJ(1));
t485 = cos(qJ(1));
t293 = -t484 * t456 - t485 * t457;
t294 = t485 * t456 - t484 * t457;
t314 = sin(qJ(5));
t316 = cos(qJ(5));
t317 = cos(qJ(4));
t423 = t316 * t317;
t252 = t293 * t314 + t294 * t423;
t433 = t314 * t317;
t253 = -t293 * t316 + t294 * t433;
t438 = t294 * t315;
t178 = Icges(6,5) * t252 - Icges(6,6) * t253 + Icges(6,3) * t438;
t444 = t178 * t317;
t453 = Icges(6,4) * t252;
t181 = -Icges(6,2) * t253 + Icges(6,6) * t438 + t453;
t249 = Icges(6,4) * t253;
t184 = Icges(6,1) * t252 + Icges(6,5) * t438 - t249;
t539 = t181 * t314 - t184 * t316;
t104 = t315 * t539 + t444;
t370 = t252 * rSges(6,1) - t253 * rSges(6,2);
t395 = rSges(6,3) * t438;
t187 = -t370 - t395;
t442 = t187 * t317;
t464 = rSges(6,1) * t316;
t368 = -rSges(6,2) * t314 + t464;
t517 = -rSges(6,3) * t317 + t368 * t315;
t543 = t517 * t438 - t442;
t372 = t315 * pkin(4) - pkin(7) * t317;
t401 = -t517 - t372;
t227 = t401 * t294;
t371 = t315 * rSges(5,1) + rSges(5,2) * t317;
t375 = t401 * t293;
t508 = m(6) / 0.2e1;
t510 = m(5) / 0.2e1;
t414 = (-t227 * t484 - t375 * t485) * t508 + (t293 * t485 + t294 * t484) * t371 * t510;
t486 = rSges(6,3) + pkin(7);
t189 = (t486 * t317 + (-pkin(4) - t368) * t315) * t294;
t439 = t293 * t317;
t284 = pkin(7) * t439;
t434 = t314 * t315;
t404 = -t293 * rSges(6,2) * t434 - rSges(6,3) * t439;
t440 = t293 * t315;
t190 = -t284 + (pkin(4) + t464) * t440 + t404;
t263 = t371 * t294;
t264 = t371 * t293;
t417 = (t484 * t189 - t485 * t190) * t508 + (-t484 * t263 - t485 * t264) * t510;
t529 = t414 - t417;
t542 = t529 * qJD(1);
t541 = -t181 * t253 + t184 * t252;
t255 = t293 * t433 + t294 * t316;
t256 = -t293 * t423 + t294 * t314;
t540 = t181 * t255 + t184 * t256;
t309 = Icges(6,4) * t434;
t426 = t315 * t316;
t449 = Icges(6,5) * t317;
t273 = -Icges(6,1) * t426 + t309 + t449;
t451 = Icges(6,4) * t316;
t360 = -Icges(6,2) * t314 + t451;
t331 = -Icges(6,6) * t317 + t360 * t315;
t357 = Icges(6,5) * t316 - Icges(6,6) * t314;
t330 = -Icges(6,3) * t317 + t357 * t315;
t427 = t315 * t330;
t147 = t252 * t273 + t253 * t331 - t294 * t427;
t538 = t147 * t317;
t454 = Icges(5,4) * t317;
t362 = -Icges(5,2) * t315 + t454;
t236 = Icges(5,6) * t293 + t362 * t294;
t429 = t236 * t315;
t432 = t315 * t178;
t525 = t294 * t432 + t541;
t534 = t293 * t525;
t359 = Icges(5,5) * t317 - Icges(5,6) * t315;
t235 = Icges(5,3) * t294 - t359 * t293;
t533 = t294 * t235;
t532 = t294 * t525;
t509 = -m(6) / 0.2e1;
t94 = t293 * t432 - t540;
t531 = t293 * t94;
t530 = t294 * t94;
t452 = Icges(6,4) * t256;
t183 = Icges(6,2) * t255 - Icges(6,6) * t440 + t452;
t250 = Icges(6,4) * t255;
t186 = Icges(6,1) * t256 - Icges(6,5) * t440 + t250;
t526 = -t253 * t183 + t186 * t252;
t493 = -t293 / 0.2e1;
t491 = -t294 / 0.2e1;
t490 = t294 / 0.2e1;
t288 = t293 * rSges(5,3);
t465 = rSges(5,1) * t317;
t305 = t315 * rSges(5,2) - t465;
t374 = -t484 * pkin(1) + t485 * qJ(2);
t340 = -t484 * pkin(2) + t374;
t334 = t293 * pkin(6) + t340;
t191 = t288 + (pkin(3) - t305) * t294 + t334;
t320 = t294 * pkin(3) + t334;
t437 = t294 * t317;
t351 = rSges(5,1) * t437 - rSges(5,2) * t438 + t288;
t193 = t320 + t351;
t524 = m(5) * (t191 - t193);
t329 = pkin(4) * t437 + pkin(7) * t438 - t187;
t138 = t320 + t329;
t471 = pkin(4) * t317;
t516 = t486 * t315 + pkin(3) + t471;
t139 = t516 * t294 + t334 + t370;
t416 = t138 - t139;
t523 = m(6) * t416;
t242 = t305 * t294 - t288;
t522 = t293 * (-t242 - t351);
t308 = -t315 * pkin(7) - t471;
t405 = t308 * t294 + t187;
t102 = (t329 + t405) * t293;
t369 = t256 * rSges(6,1) + t255 * rSges(6,2);
t460 = t315 * rSges(6,3);
t103 = t405 * t294 + ((t308 - t460) * t293 + t369) * t293;
t376 = rSges(5,2) * t440 + t294 * rSges(5,3);
t163 = t294 * t242 + t293 * (-rSges(5,1) * t439 + t376);
t519 = m(5) * t163 * t522 - m(6) * t102 * t103;
t430 = t315 * t187;
t124 = -t187 * t438 + t294 * t430;
t125 = (t294 * t369 + (-t187 - t395) * t293) * t315;
t518 = (t102 * t125 + t103 * t124) * t509;
t233 = Icges(5,3) * t293 + t359 * t294;
t455 = Icges(5,4) * t315;
t365 = Icges(5,1) * t317 - t455;
t239 = Icges(5,5) * t293 + t365 * t294;
t286 = Icges(6,2) * t426 + t309;
t287 = (Icges(6,1) * t314 + t451) * t315;
t515 = (t273 / 0.2e1 + t286 / 0.2e1) * t314 - (t287 / 0.2e1 + t331 / 0.2e1) * t316;
t514 = t293 ^ 2;
t513 = t294 ^ 2;
t512 = 0.4e1 * qJD(1);
t511 = 2 * qJD(4);
t149 = -t255 * t331 + t256 * t273 + t293 * t427;
t145 = t149 * t317;
t180 = Icges(6,5) * t256 + Icges(6,6) * t255 - Icges(6,3) * t440;
t431 = t315 * t180;
t95 = t255 * t183 + t256 * t186 - t293 * t431;
t366 = -t293 * t95 - t530;
t43 = t366 * t315 + t145;
t507 = t43 / 0.2e1;
t218 = t331 * t294;
t363 = Icges(6,1) * t316 - Icges(6,4) * t314;
t332 = t363 * t315 - t449;
t220 = t332 * t294;
t348 = t330 * t294 - t539;
t80 = t348 * t317 + (t218 * t314 - t220 * t316 + t178) * t315;
t506 = -t80 / 0.2e1;
t219 = t331 * t293;
t221 = t332 * t293;
t353 = t183 * t314 - t186 * t316;
t347 = t330 * t293 + t353;
t81 = t347 * t317 + (t219 * t314 - t221 * t316 - t180) * t315;
t505 = t81 / 0.2e1;
t197 = Icges(6,5) * t253 + Icges(6,6) * t252;
t411 = Icges(6,2) * t252 - t184 + t249;
t413 = -Icges(6,1) * t253 - t181 - t453;
t87 = t197 * t317 + (t411 * t314 + t413 * t316) * t315;
t504 = -t87 / 0.2e1;
t222 = t517 * t294;
t276 = -t368 * t317 - t460;
t422 = t317 * t517;
t121 = t430 - t222 * t317 + (-t276 * t315 + t422) * t294;
t223 = rSges(6,1) * t293 * t426 + t404;
t122 = -t315 * t369 + t317 * t223 + (-t422 + (t276 + t460) * t315) * t293;
t333 = t317 * (-rSges(6,3) * t440 + t369);
t156 = t440 * t517 - t333;
t96 = (-t222 * t315 - t442) * t293 + (t315 * t223 + t333) * t294;
t501 = m(6) * (t121 * t543 - t122 * t156 + t125 * t96);
t342 = t485 * pkin(1) + t484 * qJ(2);
t328 = t485 * pkin(2) + t342;
t319 = t294 * pkin(6) + t328;
t140 = -t293 * t516 + t319 + t369;
t500 = m(6) * (t121 * t139 + t122 * t140 - t156 * t190 + t189 * t543);
t204 = rSges(6,1) * t253 + rSges(6,2) * t252;
t205 = rSges(6,1) * t255 - rSges(6,2) * t256;
t292 = (rSges(6,1) * t314 + rSges(6,2) * t316) * t315;
t498 = m(6) * (t204 * t375 - t205 * t227 + (-t139 * t293 - t140 * t294) * t292);
t496 = m(6) * (t139 * t189 + t140 * t190);
t136 = t139 * t485;
t495 = m(6) * (t140 * t484 + t136);
t494 = t236 / 0.2e1;
t492 = -t293 / 0.4e1;
t489 = t294 / 0.4e1;
t488 = -t359 / 0.2e1;
t487 = t317 / 0.2e1;
t483 = m(3) * ((t485 * rSges(3,3) + t374) * t485 + (t484 * rSges(3,3) + t342) * t484);
t482 = m(4) * ((-t293 * rSges(4,1) - t294 * rSges(4,2) + t328) * t484 + (t294 * rSges(4,1) - t293 * rSges(4,2) + t340) * t485);
t192 = (-pkin(3) - t465) * t293 + t319 + t376;
t481 = m(5) * (-t191 * t263 + t192 * t264);
t176 = t191 * t485;
t480 = m(5) * (t192 * t484 + t176);
t476 = m(6) * (-t156 * t484 + t543 * t485);
t473 = m(6) * (-t484 * t204 - t485 * t205);
t469 = (-t178 * t293 + t180 * t294) * t315 + t94 + t526 + t540;
t467 = -t178 * t438 + t525 - t541 + t95;
t466 = m(6) * qJD(5);
t462 = t293 * t43;
t93 = -t294 * t431 - t526;
t367 = -t293 * t93 - t532;
t42 = t367 * t315 - t538;
t461 = t294 * t42;
t443 = t180 * t317;
t441 = t330 * t317;
t436 = t314 * t331;
t272 = -Icges(6,6) * t315 - t360 * t317;
t435 = t314 * t272;
t238 = Icges(5,6) * t294 - t362 * t293;
t428 = t315 * t238;
t425 = t316 * t273;
t274 = -Icges(6,5) * t315 - t363 * t317;
t424 = t316 * t274;
t285 = (Icges(6,5) * t314 + Icges(6,6) * t316) * t315;
t421 = t317 * t285;
t86 = t102 * t509 + t510 * t522;
t420 = t86 * qJD(3);
t412 = -Icges(6,1) * t255 + t183 + t452;
t410 = -Icges(6,2) * t256 + t186 + t250;
t241 = Icges(5,5) * t294 - t365 * t293;
t408 = -t293 * t235 - t241 * t437;
t407 = -t241 * t439 + t533;
t403 = -t331 - t287;
t402 = t273 + t286;
t400 = -t276 - t308;
t399 = qJD(1) * t315;
t398 = qJD(1) * t317;
t397 = qJD(5) * t315;
t396 = t124 * qJD(3);
t394 = t124 * t509;
t10 = t538 + (-t469 * t293 + t532) * t315;
t393 = -t10 / 0.2e1 - t42 / 0.2e1;
t270 = -Icges(6,3) * t315 - t357 * t317;
t352 = -t425 - t436;
t346 = t270 + t352;
t323 = -t346 * t315 + t441;
t113 = -t252 * t274 + t253 * t272 + t294 * t323;
t325 = -t348 * t315 + t444;
t62 = t253 * t218 - t220 * t252 + t294 * t325;
t324 = -t347 * t315 - t443;
t63 = t253 * t219 - t221 * t252 + t294 * t324;
t12 = (t113 + t367) * t317 + (-t293 * t63 - t294 * t62 + t147) * t315;
t66 = -t197 * t438 + t252 * t413 + t411 * t253;
t198 = Icges(6,5) * t255 - Icges(6,6) * t256;
t67 = -t198 * t438 + t252 * t412 + t410 * t253;
t33 = -t293 * t66 + t294 * t67;
t392 = -t33 / 0.2e1 + t12 / 0.2e1;
t114 = t255 * t272 + t256 * t274 + t293 * t323;
t64 = t255 * t218 + t256 * t220 + t293 * t325;
t65 = t255 * t219 + t256 * t221 + t293 * t324;
t13 = (t114 + t366) * t317 + (-t293 * t65 - t294 * t64 - t149) * t315;
t68 = -t197 * t440 + t411 * t255 - t413 * t256;
t69 = -t198 * t440 + t410 * t255 - t412 * t256;
t34 = -t293 * t68 + t294 * t69;
t391 = -t34 / 0.2e1 + t13 / 0.2e1;
t11 = t145 + (-t467 * t293 - t530) * t315;
t390 = t507 - t11 / 0.2e1;
t382 = -t440 / 0.4e1;
t380 = -t438 / 0.4e1;
t116 = t252 * t403 + t402 * t253 - t285 * t438;
t117 = t402 * t255 - t403 * t256 - t285 * t440;
t88 = t198 * t317 + (t410 * t314 + t412 * t316) * t315;
t373 = t498 / 0.2e1 + (t116 + t87) * t492 + (t117 + t88) * t489;
t364 = Icges(5,1) * t315 + t454;
t361 = Icges(5,2) * t317 + t455;
t358 = Icges(5,5) * t315 + Icges(5,6) * t317;
t135 = t204 * t294 + t205 * t293;
t70 = 0.2e1 * (-t96 / 0.4e1 + t135 / 0.4e1) * m(6);
t341 = -t484 * t293 + t485 * t294;
t326 = t341 * t292 * t508;
t335 = m(6) * (t121 * t484 - t122 * t485);
t75 = -t335 / 0.2e1 + t326;
t356 = -t75 * qJD(2) + t70 * qJD(3);
t106 = t353 * t315 + t443;
t355 = t104 * t294 - t106 * t293;
t127 = t294 * t428 + t408;
t129 = t293 * t428 + t407;
t16 = t469 * t294 + t534;
t47 = t294 * t93 - t534;
t350 = (t129 + (t233 - t428) * t293 + t533 - t407) * t293 / 0.2e1 + (t293 * t233 - (-t239 * t317 + t429) * t294) * t493 + t16 / 0.2e1 + t47 / 0.2e1 + (t239 * t439 - t293 * t429 + t127 + (t241 * t317 - t428) * t294) * t490;
t17 = t467 * t294 - t531;
t48 = t294 * t95 - t531;
t349 = t129 * t490 + (t127 - t408) * t493 + t407 * t491 + t48 / 0.2e1 - t17 / 0.2e1;
t343 = t86 * qJD(4) + qJD(5) * t394;
t259 = t361 * t294;
t261 = t364 * t294;
t339 = (t236 + t261) * t317 + (t239 - t259) * t315;
t260 = t361 * t293;
t262 = t364 * t293;
t338 = (t238 - t262) * t317 + (t241 + t260) * t315;
t337 = t10 * t489 + t11 * t492 + t17 * t380 + t462 / 0.4e1 + t461 / 0.4e1 + t48 * t438 / 0.4e1 - t518 + (t16 + t47) * t382;
t150 = t346 * t317 + (t330 - t424 + t435) * t315;
t206 = t352 * t315 - t441;
t327 = t150 * t487 - t206 * t315 / 0.2e1 + t500 / 0.2e1 + (t114 + t81) * t382 - (t106 + t149) * t439 / 0.4e1 + (t113 + t80) * t380 - (-t104 - t147) * t437 / 0.4e1;
t322 = -t425 / 0.2e1 - t436 / 0.2e1 + t270 / 0.2e1 + t364 / 0.2e1 + t362 / 0.2e1;
t321 = -t424 / 0.2e1 + t435 / 0.2e1 + t330 / 0.2e1 + t365 / 0.2e1 - t361 / 0.2e1;
t258 = t358 * t293;
t257 = t358 * t294;
t230 = t400 * t293;
t228 = t400 * t294;
t174 = t263 * t294 + t264 * t293;
t168 = t205 * t317 + t292 * t440;
t167 = -t204 * t317 - t292 * t438;
t166 = (t421 + (t402 * t314 + t403 * t316) * t315) * t317;
t143 = t473 / 0.2e1;
t133 = (-t204 * t293 + t205 * t294) * t315;
t123 = qJD(1) * t394;
t119 = (pkin(4) * t440 + t223 - t284) * t293 + (t372 * t294 + t222) * t294;
t109 = t476 / 0.2e1;
t85 = t86 * qJD(1);
t79 = -t139 * t204 + t140 * t205;
t74 = t335 / 0.2e1 + t326;
t71 = (t135 + t96) * t509;
t61 = m(6) * t79 + t421 / 0.2e1 + t515 * t315;
t59 = t103 * t135 + (t227 * t294 + t293 * t375) * t292;
t57 = t480 + t482 + t483 + t495;
t56 = t206 * t317 + t355 * t315;
t49 = t414 + t417;
t39 = t124 * t125;
t38 = t109 - t473 / 0.2e1;
t37 = t143 + t109;
t36 = t143 - t476 / 0.2e1;
t35 = t315 * t321 + t317 * t322 + t481 + t496;
t32 = -t293 * t64 + t294 * t65;
t31 = -t293 * t62 + t294 * t63;
t22 = t117 * t317 + (-t293 * t69 - t294 * t68) * t315;
t21 = t116 * t317 + (-t293 * t67 - t294 * t66) * t315;
t18 = (t150 + t355) * t317 + (-t81 * t293 - t80 * t294 - t206) * t315;
t7 = m(6) * t59 + t33 * t493 + t34 * t490;
t6 = m(6) * t39 + (t393 * t293 + t390 * t294) * t315;
t5 = t501 + (-t462 / 0.2e1 - t461 / 0.2e1 + t18 / 0.2e1) * t317 + (t13 * t493 + t12 * t491 - t56 / 0.2e1) * t315;
t4 = t349 * t293 + t350 * t294 - t519;
t3 = t327 + (-t42 / 0.4e1 - t10 / 0.4e1 + (t17 / 0.4e1 - t48 / 0.4e1) * t315) * t294 + (t11 / 0.4e1 - t43 / 0.4e1 + (t16 / 0.4e1 + t47 / 0.4e1) * t315) * t293 + t373 + t518;
t2 = t337 + (-t117 / 0.4e1 - t88 / 0.4e1) * t294 + (t116 / 0.4e1 + t87 / 0.4e1) * t293 + t327 - t498 / 0.2e1;
t1 = t337 + (-t150 / 0.2e1 + (-t147 / 0.4e1 - t104 / 0.4e1) * t294 + (t149 / 0.4e1 + t106 / 0.4e1) * t293) * t317 + (t206 / 0.2e1 + (t113 / 0.4e1 + t80 / 0.4e1) * t294 + (t114 / 0.4e1 + t81 / 0.4e1) * t293) * t315 - t500 / 0.2e1 + t373;
t8 = [(-t192 * t524 / 0.4e1 + t140 * t523 / 0.4e1) * t512 + t57 * qJD(2) + t35 * qJD(4) + t61 * qJD(5), qJD(1) * t57 + qJD(4) * t49 + qJD(5) * t37, -t343, t35 * qJD(1) + t49 * qJD(2) + t3 * qJD(5) - t420 + (m(6) * (t139 * t230 + t140 * t228 - t189 * t375 - t190 * t227) + (t505 + t114 / 0.2e1 + m(5) * (-t192 * t305 + t264 * t371) + t294 * t488 + (-t241 / 0.2e1 - t260 / 0.2e1) * t317 + (t238 / 0.2e1 - t262 / 0.2e1) * t315 - t350) * t294 + (m(5) * (-t191 * t305 - t263 * t371) + t506 - t113 / 0.2e1 + t293 * t488 + (-t239 / 0.2e1 + t259 / 0.2e1) * t317 + (t494 + t261 / 0.2e1) * t315 - t349) * t293 + t519) * qJD(4), t61 * qJD(1) + t37 * qJD(2) + t3 * qJD(4) + t166 * qJD(5) + (t396 / 0.2e1 + (t139 * t167 + t140 * t168 - t156 * t205 - t204 * t543 - t39) * qJD(5)) * m(6) + ((t504 - t116 / 0.2e1 - t390) * t294 + (-t88 / 0.2e1 - t117 / 0.2e1 - t393) * t293) * t397; -t529 * qJD(4) + t36 * qJD(5) + (-t495 / 0.4e1 - t480 / 0.4e1 - t482 / 0.4e1 - t483 / 0.4e1) * t512 + 0.2e1 * ((-t485 * t138 + t136) * t508 + (-t485 * t193 + t176) * t510) * qJD(1), 0, 0, -t542 + ((-t228 * t485 + t230 * t484) * t508 + t341 * t305 * t510) * t511 + t74 * qJD(5), t36 * qJD(1) + t74 * qJD(4) + (t167 * t484 - t168 * t485) * t466; t343, 0, 0, t85 + (-m(5) * t174 / 0.2e1 + t119 * t509) * t511 + t71 * qJD(5), t71 * qJD(4) - t133 * t466 + t123; t529 * qJD(2) + t420 + t4 * qJD(4) + t1 * qJD(5) + (-t481 / 0.4e1 - t496 / 0.4e1) * t512 - t322 * t398 - t321 * t399 + (-t227 * t523 + ((t494 - t236 / 0.2e1) * t317 - t371 * t524) * t294) * qJD(1), qJD(5) * t75 + t542, -qJD(5) * t70 + t85, t4 * qJD(1) + (m(5) * (t163 * t174 - (t513 + t514) * t305 * t371) + m(6) * (t103 * t119 - t227 * t228 - t230 * t375) + (t514 * t257 + (t338 * t294 + (-t258 + t339) * t293) * t294 + t31) * t493 + (t513 * t258 + (t339 * t293 + (-t257 + t338) * t294) * t293 + t32) * t490) * qJD(4) + t7 * qJD(5), t1 * qJD(1) + t7 * qJD(4) + (t56 / 0.2e1 + t392 * t294 + t391 * t293) * t397 - t356 + (m(6) * (t103 * t133 + t125 * t135 - t167 * t375 - t168 * t227 + (t156 * t294 - t293 * t543) * t292) + t22 * t490 + t21 * t493 - t501 + (-t18 / 0.2e1 + (t88 / 0.2e1 + t42 / 0.2e1) * t294 + (t504 + t507) * t293) * t317) * qJD(5); -t285 * t398 / 0.2e1 + t38 * qJD(2) + t2 * qJD(4) + t6 * qJD(5) + ((-t416 * t156 - t79) * qJD(1) - t396 / 0.2e1) * m(6) - t515 * t399, qJD(1) * t38 - qJD(4) * t75, qJD(4) * t70 + t123, t2 * qJD(1) + (((-t47 / 0.2e1 + t505) * t317 + (-t31 / 0.2e1 - t106 / 0.2e1) * t315 + t391) * t294 + ((-t48 / 0.2e1 + t506) * t317 + (-t32 / 0.2e1 - t104 / 0.2e1) * t315 - t392) * t293 + (t103 * t96 + t119 * t125 - t121 * t375 - t122 * t227 - t156 * t228 + t230 * t543 - t59) * m(6)) * qJD(4) + t5 * qJD(5) + t356, t6 * qJD(1) + t5 * qJD(4) + (m(6) * (t125 * t133 - t156 * t168 + t167 * t543) + t166 * t487 + (t22 * t493 + t21 * t491 + (-t293 * t88 - t294 * t87) * t487) * t315) * qJD(5);];
Cq = t8;
