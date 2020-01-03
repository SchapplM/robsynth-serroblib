% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:20
% EndTime: 2019-12-31 21:22:44
% DurationCPUTime: 11.22s
% Computational Cost: add. (18218->606), mult. (39136->857), div. (0->0), fcn. (41448->8), ass. (0->317)
t361 = sin(pkin(9));
t362 = cos(pkin(9));
t364 = sin(qJ(3));
t367 = cos(qJ(3));
t320 = -t361 * t364 + t362 * t367;
t321 = -t361 * t367 - t362 * t364;
t363 = sin(qJ(5));
t366 = cos(qJ(5));
t245 = t320 * t363 - t321 * t366;
t421 = t366 * t320 + t321 * t363;
t439 = Ifges(6,5) * t421 - Ifges(6,6) * t245;
t510 = -qJ(4) - pkin(7);
t337 = t510 * t364;
t338 = t510 * t367;
t265 = t361 * t337 - t362 * t338;
t210 = -pkin(8) * t320 - t265;
t584 = t362 * t337 + t361 * t338;
t598 = pkin(8) * t321 + t584;
t117 = t210 * t366 - t363 * t598;
t608 = t210 * t363 + t366 * t598;
t624 = t117 * mrSges(6,1) - t608 * mrSges(6,2);
t16 = t439 + t624;
t625 = t16 * qJD(5);
t365 = sin(qJ(2));
t298 = t321 * t365;
t359 = t365 * pkin(6);
t446 = t364 * t365;
t334 = pkin(3) * t446 + t359;
t248 = -pkin(4) * t298 + t334;
t525 = pkin(3) * t367;
t354 = -pkin(2) - t525;
t284 = -pkin(4) * t320 + t354;
t444 = t365 * t367;
t297 = t361 * t446 - t362 * t444;
t422 = t297 * t363 + t366 * t298;
t368 = cos(qJ(2));
t202 = t297 * t366 - t298 * t363;
t481 = t202 * mrSges(6,3);
t171 = -mrSges(6,1) * t368 + t481;
t565 = -t171 / 0.2e1;
t196 = Ifges(6,4) * t422;
t106 = -Ifges(6,1) * t202 - t368 * Ifges(6,5) + t196;
t585 = Ifges(6,2) * t202 + t106 + t196;
t600 = t202 * mrSges(6,1);
t609 = t422 * mrSges(6,2) - t600;
t599 = t245 * mrSges(6,1);
t610 = t421 * mrSges(6,2) + t599;
t235 = Ifges(6,4) * t421;
t140 = t245 * Ifges(6,1) + t235;
t597 = -Ifges(6,2) * t245 + t235;
t611 = t597 + t140;
t622 = -t117 * t565 + t585 * t421 / 0.4e1 + t611 * t422 / 0.4e1 + t248 * t610 / 0.2e1 + t284 * t609 / 0.2e1;
t620 = t248 * t609;
t619 = t284 * t610;
t617 = t265 / 0.2e1;
t616 = t365 / 0.2e1;
t555 = -t599 / 0.2e1;
t564 = t600 / 0.2e1;
t615 = t265 * mrSges(5,3);
t517 = t365 * pkin(7);
t336 = -pkin(2) * t368 - pkin(1) - t517;
t322 = t367 * t336;
t417 = -qJ(4) * t444 + t322;
t247 = (-pkin(6) * t364 - pkin(3)) * t368 + t417;
t442 = t367 * t368;
t283 = pkin(6) * t442 + t364 * t336;
t259 = -qJ(4) * t446 + t283;
t449 = t362 * t259;
t148 = t361 * t247 + t449;
t445 = t364 * t368;
t436 = pkin(6) * t445;
t258 = t417 - t436;
t153 = -t258 * t361 - t449;
t614 = t153 + t148;
t474 = t245 * Ifges(6,4);
t613 = Ifges(6,1) * t421 - t474;
t501 = Ifges(6,4) * t202;
t612 = Ifges(6,1) * t422 + t501;
t589 = Ifges(6,5) * t422;
t602 = Ifges(6,6) * t202;
t395 = t602 + t589;
t431 = -t602 / 0.2e1 - t589 / 0.2e1;
t606 = 0.2e1 * mrSges(6,2);
t605 = -t202 / 0.4e1;
t559 = t202 / 0.2e1;
t604 = t245 / 0.4e1;
t603 = -t584 / 0.2e1;
t299 = t368 * t321;
t300 = t368 * t320;
t204 = t299 * t366 - t300 * t363;
t207 = t299 * t363 + t300 * t366;
t386 = -Ifges(6,5) * t207 / 0.2e1 - Ifges(6,6) * t204 / 0.2e1;
t518 = t365 * pkin(2);
t341 = -pkin(7) * t368 + t518;
t289 = pkin(6) * t446 + t367 * t341;
t254 = t365 * pkin(3) - qJ(4) * t442 + t289;
t290 = -pkin(6) * t444 + t364 * t341;
t261 = -qJ(4) * t445 + t290;
t149 = t362 * t254 - t261 * t361;
t121 = pkin(4) * t365 - pkin(8) * t300 + t149;
t150 = t361 * t254 + t362 * t261;
t125 = pkin(8) * t299 + t150;
t48 = t121 * t366 - t125 * t363;
t49 = t121 * t363 + t125 * t366;
t416 = Ifges(6,3) * t616 - t49 * mrSges(6,2) / 0.2e1 + t48 * mrSges(6,1) / 0.2e1 - t386;
t596 = t368 / 0.2e1;
t595 = -t421 / 0.2e1;
t563 = -t422 / 0.2e1;
t562 = t422 / 0.2e1;
t496 = Ifges(4,6) * t364;
t500 = Ifges(4,5) * t367;
t388 = t500 / 0.2e1 - t496 / 0.2e1;
t586 = Ifges(3,4) - t388;
t583 = Ifges(5,5) * t320 + Ifges(5,6) * t321 + t439;
t582 = Ifges(5,5) * t298 + Ifges(5,6) * t297 + t395;
t581 = -t289 * t364 + t290 * t367;
t580 = -mrSges(4,1) * t367 + mrSges(4,2) * t364;
t250 = t361 * t259;
t147 = t362 * t247 - t250;
t154 = t362 * t258 - t250;
t579 = t154 / 0.2e1 - t147 / 0.2e1;
t577 = m(5) / 0.2e1;
t576 = m(6) / 0.2e1;
t575 = pkin(7) / 0.2e1;
t574 = m(5) * pkin(3);
t573 = -mrSges(4,2) / 0.2e1;
t570 = pkin(2) * mrSges(4,1);
t569 = pkin(2) * mrSges(4,2);
t113 = -mrSges(6,1) * t422 - mrSges(6,2) * t202;
t568 = t113 / 0.2e1;
t138 = Ifges(6,2) * t421 + t474;
t567 = -t138 / 0.2e1;
t485 = t422 * mrSges(6,3);
t169 = mrSges(6,2) * t368 + t485;
t566 = t169 / 0.2e1;
t560 = t204 / 0.2e1;
t558 = t202 / 0.4e1;
t557 = -t202 / 0.2e1;
t556 = t207 / 0.2e1;
t553 = t421 / 0.2e1;
t550 = -t245 / 0.2e1;
t549 = -t245 / 0.4e1;
t548 = -t297 / 0.2e1;
t547 = t298 / 0.2e1;
t545 = t299 / 0.2e1;
t544 = t300 / 0.2e1;
t353 = pkin(3) * t362 + pkin(4);
t527 = pkin(3) * t361;
t304 = t353 * t366 - t363 * t527;
t543 = t304 / 0.2e1;
t305 = t353 * t363 + t366 * t527;
t542 = -t305 / 0.2e1;
t541 = -t320 / 0.4e1;
t540 = t320 / 0.2e1;
t539 = -t321 / 0.4e1;
t538 = t321 / 0.2e1;
t537 = -t354 / 0.2e1;
t536 = t354 / 0.2e1;
t535 = -t364 / 0.2e1;
t534 = t364 / 0.2e1;
t532 = -t367 / 0.2e1;
t531 = t367 / 0.2e1;
t529 = -t368 / 0.4e1;
t526 = pkin(3) * t364;
t524 = pkin(4) * t321;
t523 = pkin(7) * t364;
t522 = pkin(7) * t367;
t521 = pkin(8) * t298;
t519 = t297 * pkin(8);
t360 = t368 * pkin(6);
t115 = -pkin(4) * t368 + t147 + t519;
t122 = t148 + t521;
t44 = t115 * t366 - t122 * t363;
t516 = t44 * mrSges(6,2);
t45 = t115 * t363 + t122 * t366;
t515 = t45 * mrSges(6,1);
t128 = t153 - t521;
t129 = t154 + t519;
t53 = t128 * t366 - t129 * t363;
t512 = t53 * mrSges(6,1);
t54 = t128 * t363 + t129 * t366;
t511 = t54 * mrSges(6,2);
t507 = mrSges(6,3) * t304;
t506 = Ifges(4,4) * t364;
t505 = Ifges(4,4) * t367;
t504 = Ifges(5,4) * t297;
t503 = Ifges(5,4) * t320;
t502 = Ifges(5,4) * t321;
t491 = pkin(3) * qJD(3);
t484 = t204 * mrSges(6,1);
t480 = t207 * mrSges(6,2);
t479 = t245 * mrSges(6,3);
t473 = t584 * mrSges(5,3);
t472 = t297 * mrSges(5,3);
t471 = t298 * mrSges(5,3);
t470 = t299 * mrSges(5,1);
t104 = Ifges(6,2) * t422 - t368 * Ifges(6,6) - t501;
t105 = Ifges(6,4) * t207 + Ifges(6,2) * t204 + Ifges(6,6) * t365;
t107 = Ifges(6,1) * t207 + Ifges(6,4) * t204 + Ifges(6,5) * t365;
t114 = t480 - t484;
t170 = -mrSges(6,2) * t365 + mrSges(6,3) * t204;
t172 = mrSges(6,1) * t365 - mrSges(6,3) * t207;
t197 = Ifges(5,2) * t298 - t368 * Ifges(5,6) - t504;
t198 = Ifges(5,4) * t300 + Ifges(5,2) * t299 + Ifges(5,6) * t365;
t288 = Ifges(5,4) * t298;
t199 = -Ifges(5,1) * t297 - t368 * Ifges(5,5) + t288;
t200 = Ifges(5,1) * t300 + Ifges(5,4) * t299 + Ifges(5,5) * t365;
t214 = -mrSges(5,1) * t298 - mrSges(5,2) * t297;
t468 = t300 * mrSges(5,2);
t215 = t468 - t470;
t335 = pkin(3) * t445 + t360;
t249 = -pkin(4) * t299 + t335;
t268 = mrSges(5,2) * t368 + t471;
t269 = -mrSges(5,2) * t365 + mrSges(5,3) * t299;
t270 = -mrSges(5,1) * t368 + t472;
t271 = mrSges(5,1) * t365 - mrSges(5,3) * t300;
t282 = t322 - t436;
t403 = -Ifges(4,2) * t364 + t505;
t294 = Ifges(4,6) * t365 + t368 * t403;
t410 = Ifges(4,1) * t367 - t506;
t296 = Ifges(4,5) * t365 + t368 * t410;
t414 = mrSges(4,1) * t364 + mrSges(4,2) * t367;
t314 = t368 * t414;
t330 = mrSges(4,2) * t368 - mrSges(4,3) * t446;
t331 = -t365 * mrSges(4,2) - mrSges(4,3) * t445;
t332 = -mrSges(4,1) * t368 - mrSges(4,3) * t444;
t333 = t365 * mrSges(4,1) - mrSges(4,3) * t442;
t387 = -Ifges(5,5) * t300 / 0.2e1 - Ifges(5,6) * t299 / 0.2e1;
t382 = t410 * t365;
t465 = t368 * Ifges(4,5);
t295 = t382 - t465;
t443 = t367 * t295;
t381 = t403 * t365;
t464 = t368 * Ifges(4,6);
t293 = t381 - t464;
t447 = t364 * t293;
t3 = (-pkin(1) * mrSges(3,2) - t447 / 0.2e1 + t443 / 0.2e1 + t386 + t387 + t586 * t368) * t368 + (Ifges(5,5) * t548 + Ifges(5,6) * t547 + Ifges(6,5) * t557 + Ifges(6,6) * t562 - pkin(1) * mrSges(3,1) + pkin(6) * t314 + t294 * t535 + t296 * t531 - t586 * t365 + (-Ifges(3,2) + Ifges(3,1) - Ifges(5,3) - Ifges(6,3) - Ifges(4,3) + (m(4) * pkin(6) + t414) * pkin(6)) * t368) * t365 + m(4) * (t282 * t289 + t283 * t290) + m(6) * (t248 * t249 + t44 * t48 + t45 * t49) + m(5) * (t147 * t149 + t148 * t150 + t334 * t335) + t107 * t557 + t104 * t560 + t105 * t562 + t199 * t544 + t197 * t545 + t198 * t547 + t200 * t548 + t106 * t556 + t290 * t330 + t283 * t331 + t289 * t332 + t282 * t333 + t334 * t215 + t335 * t214 + t150 * t268 + t148 * t269 + t149 * t270 + t147 * t271 + t248 * t114 + t249 * t113 + t48 * t171 + t44 * t172 + t49 * t169 + t45 * t170;
t469 = t3 * qJD(1);
t467 = t320 * mrSges(5,3);
t466 = t321 * mrSges(5,3);
t437 = pkin(3) * t444;
t260 = -pkin(4) * t297 + t437;
t397 = Ifges(4,5) * t364 + Ifges(4,6) * t367;
t402 = Ifges(4,2) * t367 + t506;
t408 = Ifges(5,1) * t298 + t504;
t409 = Ifges(4,1) * t364 + t505;
t424 = Ifges(5,2) * t297 + t288;
t285 = t297 * mrSges(5,1);
t428 = t298 * mrSges(5,2) - t285;
t4 = t282 * t330 - t283 * t332 + t334 * t428 + t408 * t548 + t297 * t197 / 0.2e1 + t154 * t268 + t153 * t270 + t260 * t113 + t620 + t612 * t557 + t104 * t559 + t53 * t171 + t54 * t169 + (t202 * t45 - t422 * t44) * mrSges(6,3) + (-t147 * t298 + t148 * t297) * mrSges(5,3) + m(5) * (t147 * t153 + t148 * t154) + m(6) * (t248 * t260 + t44 * t53 + t45 * t54) + (t397 * t596 + t295 * t535 + t293 * t532 + (-pkin(6) * t580 + t402 * t534 + t409 * t532) * t365 + (m(5) * t334 + t214) * t525 + (t282 * t364 - t283 * t367) * mrSges(4,3)) * t365 + (t424 + t199) * t547 - t582 * t368 / 0.2e1 + t585 * t562;
t463 = t4 * qJD(1);
t462 = t44 * t421;
t461 = t45 * t245;
t460 = t53 * t245;
t459 = t54 * t421;
t7 = t395 * t596 - t620 + t612 * t559 + t104 * t557 + t585 * t563 + (-t481 + t171) * t45 + (-t169 + t485) * t44;
t458 = t7 * qJD(1);
t18 = t422 * t169 + t202 * t171 + t298 * t268 + t297 * t270 + m(6) * (t202 * t44 + t422 * t45) + m(5) * (t147 * t297 + t148 * t298);
t457 = qJD(1) * t18;
t456 = t608 * t422;
t455 = t117 * t202;
t454 = t202 * t304;
t453 = t422 * t305;
t450 = t334 * t364;
t94 = -mrSges(6,1) * t305 - mrSges(6,2) * t304;
t441 = t94 * qJD(5);
t435 = mrSges(4,3) * t575;
t434 = t526 / 0.2e1;
t433 = Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t427 = -t321 * mrSges(5,1) + t320 * mrSges(5,2);
t420 = -mrSges(4,3) * t517 / 0.2e1;
t419 = t437 * t577;
t418 = -t524 + t526;
t413 = -t320 * mrSges(5,1) - t321 * mrSges(5,2);
t411 = -mrSges(6,1) * t421 + t245 * mrSges(6,2);
t407 = -Ifges(5,1) * t321 + t503;
t401 = Ifges(5,2) * t320 - t502;
t398 = -t496 + t500;
t369 = -m(6) * (t248 * t418 + t284 * t260 + (t45 + t53) * t608 + (t44 - t54) * t117) / 0.2e1 - t608 * t169 / 0.2e1 + t138 * t605 + t613 * t558 + t104 * t604 + t612 * t549 - t260 * t411 / 0.2e1 + t270 * t617 + t268 * t603 + t199 * t541 + t197 * t539 - t334 * t427 / 0.2e1 + t583 * t368 / 0.4e1 - t622;
t372 = (t304 * t48 + t305 * t49) * t576 + t149 * mrSges(5,1) / 0.2e1 - t150 * mrSges(5,2) / 0.2e1 + t289 * mrSges(4,1) / 0.2e1 + t290 * t573 + t172 * t543 + t305 * t170 / 0.2e1 - t387 + t416;
t377 = t614 * t584 + (-t147 + t154) * t265;
t379 = (t149 * t362 + t150 * t361) * t577;
t383 = t361 * t269 / 0.2e1 + t362 * t271 / 0.2e1;
t384 = pkin(3) * t413;
t1 = t369 + t288 * t541 + t285 * t536 + (t455 / 0.2e1 + t456 / 0.2e1 + t461 / 0.2e1 - t459 / 0.2e1 + t462 / 0.2e1 + t460 / 0.2e1) * mrSges(6,3) + pkin(3) * t379 + (-t295 / 0.4e1 + 0.3e1 / 0.4e1 * t465 + t332 * t575 + (t570 / 0.2e1 + 0.3e1 / 0.2e1 * t506 + pkin(6) * t573 + t537 * t574 - t384 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.4e1 + t435) * t367) * t365) * t367 + ((-t569 / 0.2e1 - pkin(6) * mrSges(4,1) / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.4e1 + t435) * t364) * t364 + t433) * t365 + (0.3e1 / 0.4e1 * t502 - t615 / 0.2e1 + (-Ifges(5,2) / 0.2e1 + Ifges(5,1) / 0.4e1) * t320) * t297 + (mrSges(5,2) * t537 - t503 / 0.2e1 + t473 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.4e1) * t321) * t298 + t524 * t568 + ((-t148 / 0.2e1 - t153 / 0.2e1) * t321 - t579 * t320) * mrSges(5,3) + (t293 / 0.4e1 + t330 * t575 - 0.3e1 / 0.4e1 * t464 + (-t113 / 0.2e1 - t214 / 0.2e1) * pkin(3)) * t364 - m(5) * (pkin(3) * t450 + t377) / 0.2e1 + t383 * pkin(3) + t372;
t8 = -t619 - t265 * t466 - t117 * t479 + t613 * t550 + (t502 + t615) * t321 + (t569 - t505) * t367 + (t473 - t503 + (Ifges(5,1) - Ifges(5,2)) * t321) * t320 + (-t384 + t506 + t570 + (-Ifges(4,1) + Ifges(4,2)) * t367) * t364 + t611 * t595 + (-m(5) * t526 - t427) * t354 - t584 * t467 + (-m(6) * t284 - t411) * t418 + (t117 * mrSges(6,3) - t567) * t245;
t393 = -t1 * qJD(1) - t8 * qJD(2);
t14 = t619 + (t613 / 0.2e1 + t567) * t245 + (t140 / 0.2e1 + t597 / 0.2e1) * t421;
t370 = (-t117 * t559 + t563 * t608) * mrSges(6,3) + t608 * t566 + t138 * t558 + t613 * t605 + t104 * t549 + t612 * t604 + t439 * t529 + t622;
t6 = t370 - t416;
t392 = t6 * qJD(1) + t14 * qJD(2);
t371 = (t202 * t550 + t422 * t553) * mrSges(6,3) + (t297 * t538 + t298 * t540) * mrSges(5,3) + (t147 * t321 + t148 * t320 + t265 * t298 + t297 * t584) * t577 + (-t117 * t422 + t202 * t608 - t245 * t44 + t421 * t45) * t576 + t245 * t565 + t421 * t566 + t268 * t540 + t270 * t538;
t374 = t335 * t577 + t249 * t576 - t484 / 0.2e1 + t480 / 0.2e1 - t470 / 0.2e1 + t468 / 0.2e1;
t12 = t371 - t374;
t22 = (t245 ^ 2 + t421 ^ 2) * mrSges(6,3) + (t320 ^ 2 + t321 ^ 2) * mrSges(5,3) + m(6) * (-t117 * t421 - t245 * t608) + m(5) * (t265 * t320 + t321 * t584);
t391 = -qJD(1) * t12 - qJD(2) * t22;
t380 = (t297 * t362 + t298 * t361) * t574;
t27 = t419 + 0.2e1 * (t260 / 0.4e1 - t454 / 0.4e1 - t453 / 0.4e1) * m(6) - t380 / 0.2e1 + t609 + t428;
t375 = (-t245 * t304 + t305 * t421) * t576 + (t320 * t361 + t321 * t362) * t574 / 0.2e1;
t376 = m(5) * t434 + t418 * t576;
t28 = -t375 + t376 + t610 + t427;
t390 = qJD(1) * t27 + qJD(2) * t28;
t42 = t563 * t606 + 0.2e1 * t564;
t56 = t595 * t606 + 0.2e1 * t555;
t389 = qJD(1) * t42 + qJD(2) * t56;
t373 = (-t202 * t542 + t304 * t563) * mrSges(6,3) + t169 * t543 + t171 * t542 - t431;
t10 = (-t44 / 0.2e1 + t54 / 0.2e1) * mrSges(6,2) + (-t45 / 0.2e1 - t53 / 0.2e1) * mrSges(6,1) + t373 + t431;
t378 = t10 * qJD(1) + t94 * qJD(3);
t57 = t599 / 0.2e1 + t555;
t50 = t375 + t376;
t43 = -t600 / 0.2e1 + t564;
t37 = t380 / 0.2e1 + t419 + (t453 + t454 + t260) * t576;
t11 = t371 + t374;
t9 = -t516 / 0.2e1 - t515 / 0.2e1 - t511 / 0.2e1 + t512 / 0.2e1 + t373 - t431;
t5 = t370 + t416;
t2 = t443 / 0.4e1 - t369 + t388 * t368 + t297 * t401 / 0.4e1 + t472 * t617 + (t379 + t383) * pkin(3) + t579 * t467 + (t382 / 0.4e1 + t420 * t367) * t367 + (-t381 / 0.4e1 + t420 * t364) * t364 + (-t461 + t459 - t455) * mrSges(6,3) / 0.2e1 + (Ifges(5,2) * t321 + t407 + t503) * t298 / 0.4e1 - (t462 + t460 + t456) * mrSges(6,3) / 0.2e1 + t214 * t434 + ((t354 * t444 + t450) * pkin(3) + t377) * t577 + t418 * t568 + t408 * t539 - t330 * t523 / 0.2e1 + t414 * t359 / 0.2e1 - t332 * t522 / 0.2e1 - t297 * (Ifges(5,1) * t320 + t502) / 0.4e1 + t580 * t518 / 0.2e1 - t409 * t446 / 0.2e1 - t402 * t444 / 0.2e1 + t384 * t444 / 0.2e1 + t433 * t365 - t447 / 0.4e1 + t614 * t466 / 0.2e1 + t372 + t471 * t603 + t398 * t529 + t428 * t536 + t320 * t424 / 0.4e1;
t13 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t18 - qJD(5) * t7, t2 * qJD(3) + t11 * qJD(4) + t5 * qJD(5) + t469 + (t249 * t411 + t335 * t413 + t49 * t421 * mrSges(6,3) + (-Ifges(5,5) * t321 + Ifges(6,5) * t245 + Ifges(5,6) * t320 + Ifges(6,6) * t421 + t397) * t616 + 0.2e1 * (-t117 * t49 + t249 * t284 + t48 * t608) * t576 - t117 * t170 + t581 * mrSges(4,3) + m(4) * (-pkin(2) * t360 + t581 * pkin(7)) + (t107 / 0.2e1 - t48 * mrSges(6,3)) * t245 + 0.2e1 * (t149 * t584 + t150 * t265 + t335 * t354) * t577 + t584 * t271 + t608 * t172 + t138 * t560 + t198 * t540 + t407 * t544 + t401 * t545 + t105 * t553 + t140 * t556 - Ifges(3,6) * t365 + t354 * t215 - t321 * t200 / 0.2e1 - pkin(2) * t314 + t284 * t114 + t265 * t269 - t333 * t523 + (Ifges(3,5) + t409 * t531 + t402 * t535 + (-mrSges(3,1) + t580) * pkin(6)) * t368 + t149 * t466 + t150 * t467 + mrSges(3,2) * t359 + t331 * t522 + t294 * t531 + t296 * t534) * qJD(2), t463 + t2 * qJD(2) + (-t282 * mrSges(4,2) - t283 * mrSges(4,1) - Ifges(4,5) * t446 - Ifges(4,6) * t444 + t305 * t481 - t422 * t507 + m(6) * (t304 * t53 + t305 * t54) - t511 + t512 + t153 * mrSges(5,1) - t154 * mrSges(5,2) + t582) * qJD(3) + t37 * qJD(4) + t9 * qJD(5) + (m(5) * (t153 * t362 + t154 * t361) + (t297 * t361 - t298 * t362) * mrSges(5,3)) * t491, qJD(2) * t11 + qJD(3) * t37 + qJD(5) * t43 + t457, -t458 + t5 * qJD(2) + t9 * qJD(3) + t43 * qJD(4) + (t395 - t515 - t516) * qJD(5); -qJD(3) * t1 + qJD(4) * t12 + qJD(5) * t6 - t469, -qJD(3) * t8 + qJD(4) * t22 + qJD(5) * t14, (-t305 * t479 - t421 * t507 + m(6) * (t117 * t304 + t305 * t608) - t584 * mrSges(5,2) - t265 * mrSges(5,1) + t398 + t580 * pkin(7) + t583 + t624) * qJD(3) + t50 * qJD(4) + t625 + (m(5) * (-t265 * t362 + t361 * t584) + (-t320 * t362 + t321 * t361) * mrSges(5,3)) * t491 + t393, qJD(3) * t50 + qJD(5) * t57 - t391, t16 * qJD(3) + t57 * qJD(4) + t392 + t625; qJD(2) * t1 - qJD(4) * t27 + qJD(5) * t10 - t463, -qJD(4) * t28 - t393, t441, -t390, t378 + t441; -qJD(2) * t12 + qJD(3) * t27 - qJD(5) * t42 - t457, qJD(3) * t28 - qJD(5) * t56 + t391, t390, 0, -t389; -qJD(2) * t6 - qJD(3) * t10 + qJD(4) * t42 + t458, qJD(4) * t56 - t392, -t378, t389, 0;];
Cq = t13;
