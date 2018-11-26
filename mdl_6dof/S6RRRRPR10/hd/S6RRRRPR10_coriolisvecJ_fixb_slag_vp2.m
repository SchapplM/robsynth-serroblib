% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:19:18
% EndTime: 2018-11-23 18:19:53
% DurationCPUTime: 35.99s
% Computational Cost: add. (20223->859), mult. (52101->1159), div. (0->0), fcn. (40696->10), ass. (0->417)
t316 = cos(qJ(2));
t314 = sin(qJ(2));
t310 = sin(pkin(6));
t452 = qJD(1) * t310;
t425 = t314 * t452;
t311 = cos(pkin(6));
t451 = qJD(1) * t311;
t440 = pkin(1) * t451;
t263 = -pkin(8) * t425 + t316 * t440;
t372 = (pkin(2) * t314 - pkin(9) * t316) * t310;
t264 = qJD(1) * t372;
t487 = sin(qJ(3));
t489 = cos(qJ(3));
t202 = -t263 * t487 + t489 * t264;
t443 = t489 * pkin(9);
t290 = pkin(10) * t489 + t443;
t429 = t316 * t489;
t609 = -(pkin(3) * t314 - pkin(10) * t429) * t452 - t202 - qJD(3) * t290;
t203 = t489 * t263 + t487 * t264;
t441 = t487 * pkin(9);
t289 = -pkin(10) * t487 - t441;
t430 = t310 * t487;
t401 = qJD(1) * t430;
t378 = t316 * t401;
t608 = -pkin(10) * t378 - t289 * qJD(3) + t203;
t313 = sin(qJ(4));
t488 = cos(qJ(4));
t403 = t488 * t489;
t417 = t487 * qJD(3);
t550 = qJD(3) + qJD(4);
t222 = (qJD(4) * t487 + t417) * t313 - t550 * t403;
t431 = t310 * t489;
t402 = qJD(1) * t431;
t379 = t316 * t402;
t231 = -t313 * t378 + t379 * t488;
t607 = t222 + t231;
t278 = t313 * t489 + t487 * t488;
t223 = t550 * t278;
t424 = t316 * t452;
t230 = t278 * t424;
t601 = t223 - t230;
t294 = t314 * t440;
t266 = pkin(8) * t424 + t294;
t556 = -t266 + (-t378 + t417) * pkin(3);
t411 = qJD(2) + t451;
t363 = t487 * t411;
t246 = t314 * t402 + t363;
t364 = t489 * t411;
t335 = t314 * t401 - t364;
t331 = t313 * t335;
t324 = t246 * t488 - t331;
t500 = -t324 / 0.2e1;
t606 = Ifges(6,2) * t500;
t605 = qJ(5) * t607 - qJD(5) * t278 + t556;
t234 = t313 * t289 + t290 * t488;
t562 = -qJD(4) * t234 + t608 * t313 + t488 * t609;
t235 = t488 * t335;
t428 = t316 * t487;
t399 = qJD(2) * t428;
t376 = t310 * t399;
t322 = qJD(1) * t376 + qJD(3) * t246;
t330 = t335 * qJD(3);
t400 = qJD(2) * t429;
t377 = t310 * t400;
t323 = qJD(1) * t377 - t330;
t449 = qJD(4) * t313;
t109 = qJD(4) * t235 + t246 * t449 + t313 * t322 - t323 * t488;
t521 = -t109 / 0.2e1;
t418 = qJD(4) * t488;
t110 = -qJD(4) * t331 + t246 * t418 + t313 * t323 + t322 * t488;
t519 = -t110 / 0.2e1;
t191 = t246 * t313 + t235;
t503 = -t191 / 0.2e1;
t288 = qJD(3) - t424;
t281 = -qJD(4) - t288;
t493 = -t281 / 0.2e1;
t492 = t281 / 0.2e1;
t499 = t324 / 0.2e1;
t583 = Ifges(6,1) + Ifges(5,3);
t582 = Ifges(6,4) - Ifges(5,5);
t581 = -Ifges(5,6) + Ifges(6,5);
t458 = t310 * t314;
t522 = pkin(4) + pkin(11);
t409 = t522 * t458;
t604 = -pkin(5) * t607 + qJD(1) * t409 - t562;
t229 = pkin(9) * t411 + t266;
t259 = (-pkin(2) * t316 - pkin(9) * t314 - pkin(1)) * t310;
t241 = qJD(1) * t259;
t180 = -t229 * t487 + t241 * t489;
t153 = -t246 * pkin(10) + t180;
t352 = -t229 * t489 - t241 * t487;
t154 = -pkin(10) * t335 - t352;
t455 = t313 * t154;
t81 = t153 * t488 - t455;
t603 = pkin(3) * t418 - t81;
t561 = t289 * t418 - t290 * t449 + t313 * t609 - t608 * t488;
t602 = t522 * t601 + t605;
t518 = t110 / 0.2e1;
t600 = t518 - t519;
t484 = t191 * pkin(5);
t188 = Ifges(5,4) * t191;
t126 = Ifges(5,1) * t324 - t281 * Ifges(5,5) - t188;
t189 = qJD(6) + t324;
t466 = t189 * Ifges(7,3);
t312 = sin(qJ(6));
t315 = cos(qJ(6));
t162 = t191 * t312 - t281 * t315;
t467 = t162 * Ifges(7,5);
t161 = t191 * t315 + t281 * t312;
t468 = t161 * Ifges(7,6);
t68 = t466 + t467 + t468;
t573 = t126 + t68;
t571 = -qJD(5) - t603;
t462 = qJ(5) * t191;
t564 = qJ(5) * t425 - t561;
t382 = qJD(2) * t294;
t445 = qJD(1) * qJD(2);
t416 = t316 * t445;
t398 = t310 * t416;
t35 = pkin(3) * t322 + t110 * pkin(4) + pkin(8) * t398 + t109 * qJ(5) - qJD(5) * t324 + t382;
t16 = t110 * pkin(11) + t35;
t585 = pkin(5) * t324;
t337 = t288 * pkin(3) + t153;
t140 = t488 * t337;
t74 = -t140 + t455;
t373 = t74 + t585;
t595 = qJD(5) + t373;
t51 = t281 * t522 + t595;
t228 = -pkin(2) * t411 - t263;
t194 = pkin(3) * t335 + t228;
t318 = -qJ(5) * t324 + t194;
t64 = t191 * t522 + t318;
t23 = -t312 * t64 + t315 * t51;
t334 = t313 * t337;
t265 = qJD(2) * t372;
t253 = qJD(1) * t265;
t300 = pkin(8) * t458;
t483 = t316 * pkin(1);
t273 = t311 * t483 - t300;
t267 = t273 * qJD(2);
t254 = qJD(1) * t267;
t390 = t489 * t253 - t254 * t487;
t450 = qJD(2) * t314;
t423 = t310 * t450;
t397 = qJD(1) * t423;
t419 = qJD(3) * t489;
t86 = pkin(3) * t397 - pkin(10) * t323 - t229 * t419 - t241 * t417 + t390;
t129 = -t229 * t417 + t241 * t419 + t487 * t253 + t489 * t254;
t92 = -pkin(10) * t322 + t129;
t366 = qJD(4) * t334 + t154 * t418 + t313 * t92 - t488 * t86;
t381 = qJD(2) * t409;
t7 = -t109 * pkin(5) - qJD(1) * t381 + t366;
t1 = qJD(6) * t23 + t16 * t315 + t312 * t7;
t24 = t312 * t51 + t315 * t64;
t461 = qJD(6) * t24;
t2 = -t16 * t312 + t315 * t7 - t461;
t383 = t23 * t312 - t24 * t315;
t336 = m(7) * (-qJD(6) * t383 + t1 * t312 + t2 * t315);
t62 = qJD(6) * t161 + t110 * t312 + t315 * t397;
t38 = -mrSges(7,1) * t109 - mrSges(7,3) * t62;
t63 = -qJD(6) * t162 + t110 * t315 - t312 * t397;
t39 = mrSges(7,2) * t109 + mrSges(7,3) * t63;
t599 = t312 * t39 + t315 * t38 + t336;
t388 = mrSges(7,1) * t315 - mrSges(7,2) * t312;
t471 = Ifges(7,4) * t162;
t69 = t161 * Ifges(7,2) + t189 * Ifges(7,6) + t471;
t525 = -t69 / 0.2e1;
t150 = t488 * t154;
t75 = t150 + t334;
t72 = t281 * qJ(5) - t75;
t53 = -t72 - t484;
t598 = t315 * t525 + t53 * t388;
t565 = -qJD(5) - t74;
t71 = pkin(4) * t281 - t565;
t597 = t71 * mrSges(6,1) + t74 * mrSges(5,3);
t596 = t72 * mrSges(6,1) - t75 * mrSges(5,3);
t187 = Ifges(6,6) * t191;
t87 = t191 * pkin(4) + t318;
t540 = Ifges(6,4) * t493 + t606 + t187 / 0.2e1 - t23 * mrSges(7,1) - t194 * mrSges(5,2) + t24 * mrSges(7,2) + t87 * mrSges(6,3);
t502 = t191 / 0.2e1;
t544 = -t194 * mrSges(5,1) + t87 * mrSges(6,2) + Ifges(6,5) * t492 + Ifges(5,6) * t493 + (Ifges(5,2) + Ifges(6,3)) * t503 + (Ifges(5,4) + Ifges(6,6)) * t499;
t594 = -Ifges(5,4) * t499 - Ifges(5,2) * t503 + Ifges(6,6) * t500 + Ifges(6,3) * t502 + t493 * t581 - t544 + t596;
t593 = t573 / 0.2e1 + t597;
t10 = Ifges(7,5) * t62 + Ifges(7,6) * t63 - Ifges(7,3) * t109;
t345 = -t314 * t419 - t399;
t351 = qJD(3) * t363;
t446 = t316 * qJD(2);
t179 = t382 + pkin(3) * t351 + (-pkin(3) * t345 + pkin(8) * t446) * t452;
t18 = -pkin(4) * t397 + t366;
t584 = qJD(2) / 0.2e1;
t413 = t310 * t584;
t392 = qJD(1) * t413;
t375 = t314 * t392;
t520 = t109 / 0.2e1;
t526 = t63 / 0.2e1;
t527 = t62 / 0.2e1;
t549 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t587 = -t397 / 0.2e1;
t592 = 0.2e1 * Ifges(5,4) * t519 + Ifges(5,5) * t397 / 0.2e1 + t10 / 0.2e1 - t582 * t375 + mrSges(6,1) * t18 + mrSges(5,2) * t179 + mrSges(5,3) * t366 - mrSges(6,3) * t35 + Ifges(7,5) * t527 + Ifges(7,6) * t526 + Ifges(6,4) * t587 + t549 - t600 * Ifges(6,6) - Ifges(6,2) * t520 + (0.2e1 * Ifges(5,1) + Ifges(7,3) + Ifges(6,2)) * t521;
t385 = Ifges(7,5) * t312 + Ifges(7,6) * t315;
t470 = Ifges(7,4) * t312;
t386 = Ifges(7,2) * t315 + t470;
t469 = Ifges(7,4) * t315;
t387 = Ifges(7,1) * t312 + t469;
t490 = -t312 / 0.2e1;
t505 = -t189 / 0.2e1;
t508 = -t162 / 0.2e1;
t510 = -t161 / 0.2e1;
t160 = Ifges(7,4) * t161;
t70 = t162 * Ifges(7,1) + t189 * Ifges(7,5) + t160;
t591 = t385 * t505 + t386 * t510 + t387 * t508 + t70 * t490 - t596 + t598;
t590 = Ifges(6,5) / 0.2e1;
t589 = -t322 / 0.2e1;
t588 = t323 / 0.2e1;
t586 = pkin(4) * t324;
t277 = t313 * t487 - t403;
t308 = -pkin(3) * t489 - pkin(2);
t362 = -t278 * qJ(5) + t308;
t195 = t277 * t522 + t362;
t233 = -t488 * t289 + t290 * t313;
t206 = pkin(5) * t278 + t233;
t138 = -t195 * t312 + t206 * t315;
t580 = qJD(6) * t138 + t312 * t604 + t315 * t602;
t139 = t195 * t315 + t206 * t312;
t579 = -qJD(6) * t139 - t312 * t602 + t315 * t604;
t25 = -mrSges(7,1) * t63 + mrSges(7,2) * t62;
t93 = mrSges(6,1) * t110 - mrSges(6,3) * t397;
t577 = -t93 + t25;
t576 = Ifges(4,5) * t323;
t575 = Ifges(4,6) * t322;
t574 = t254 * mrSges(3,2);
t572 = t585 - t571;
t570 = -pkin(5) * t601 - t564;
t569 = t324 * t522;
t563 = pkin(4) * t425 - t562;
t560 = pkin(4) * t601 + t605;
t169 = mrSges(6,1) * t324 - mrSges(6,2) * t281;
t171 = -mrSges(5,1) * t281 - mrSges(5,3) * t324;
t559 = -t169 + t171;
t204 = t230 * t315 - t312 * t425;
t448 = qJD(6) * t312;
t454 = t315 * t223;
t368 = t277 * t448 - t454;
t558 = t204 + t368;
t205 = t230 * t312 + t315 * t425;
t447 = qJD(6) * t315;
t456 = t312 * t223;
t369 = t277 * t447 + t456;
t557 = t205 - t369;
t555 = t246 * mrSges(4,1) - t335 * mrSges(4,2);
t130 = qJD(3) * t352 + t390;
t554 = t489 * t129 - t487 * t130;
t553 = t109 * t582 + t110 * t581 + t397 * t583;
t303 = t311 * t314 * pkin(1);
t457 = t310 * t316;
t274 = pkin(8) * t457 + t303;
t258 = pkin(9) * t311 + t274;
t404 = t314 * t430;
t350 = -t311 * t489 + t404;
t329 = -qJD(3) * t350 + t377;
t389 = t489 * t265 - t267 * t487;
t111 = pkin(3) * t423 - pkin(10) * t329 - t258 * t419 - t259 * t417 + t389;
t141 = -t258 * t417 + t259 * t419 + t487 * t265 + t489 * t267;
t269 = t311 * t487 + t314 * t431;
t328 = qJD(3) * t269 + t376;
t119 = -pkin(10) * t328 + t141;
t200 = -t258 * t487 + t489 * t259;
t159 = -pkin(3) * t457 - t269 * pkin(10) + t200;
t201 = t489 * t258 + t487 * t259;
t172 = -pkin(10) * t350 + t201;
t359 = -t313 * t111 - t488 * t119 - t159 * t418 + t172 * t449;
t26 = -t310 * (qJ(5) * t450 - qJD(5) * t316) + t359;
t408 = mrSges(3,3) * t425;
t548 = -m(4) * t228 + mrSges(3,1) * t411 - mrSges(4,1) * t335 - t246 * mrSges(4,2) - t408;
t547 = t314 * (Ifges(4,5) * t246 - Ifges(4,6) * t335 + Ifges(4,3) * t288 + t191 * t581 - t281 * t583 - t324 * t582);
t545 = -t130 * mrSges(4,1) + t129 * mrSges(4,2);
t20 = qJD(4) * t140 - t154 * t449 + t313 * t86 + t488 * t92;
t15 = -qJ(5) * t397 + qJD(5) * t281 - t20;
t543 = mrSges(5,1) * t366 + t20 * mrSges(5,2) - t18 * mrSges(6,2) + t15 * mrSges(6,3);
t542 = -mrSges(5,1) * t74 - mrSges(5,2) * t75 + mrSges(6,2) * t71 - mrSges(6,3) * t72;
t539 = mrSges(5,1) * t179 + mrSges(6,1) * t15 - mrSges(6,2) * t35 - mrSges(5,3) * t20 + Ifges(5,6) * t587 + 0.2e1 * Ifges(6,3) * t518 + t397 * t590 + (Ifges(5,4) + 0.2e1 * Ifges(6,6)) * t520 + t600 * Ifges(5,2);
t538 = -Ifges(5,4) * t500 - Ifges(5,2) * t502 + Ifges(6,6) * t499 + Ifges(6,3) * t503 + t492 * t581 + t544;
t535 = Ifges(5,1) * t500 + Ifges(5,4) * t502 + Ifges(7,5) * t508 - Ifges(6,2) * t499 - Ifges(6,6) * t503 + Ifges(7,6) * t510 + Ifges(7,3) * t505 - t492 * t582 + t540;
t504 = t189 / 0.2e1;
t507 = t162 / 0.2e1;
t509 = t161 / 0.2e1;
t534 = -Ifges(5,1) * t499 - Ifges(5,4) * t503 - Ifges(7,5) * t507 + Ifges(6,6) * t502 - Ifges(7,6) * t509 - Ifges(7,3) * t504 + t493 * t582 + t540 + t606;
t533 = t310 ^ 2;
t531 = Ifges(7,1) * t527 + Ifges(7,4) * t526 + Ifges(7,5) * t521;
t524 = t69 / 0.2e1;
t523 = t70 / 0.2e1;
t344 = -t314 * t417 + t400;
t320 = (Ifges(4,1) * t344 + Ifges(4,4) * t345 + Ifges(4,5) * t450) * t310;
t511 = (Ifges(4,1) * t364 - Ifges(4,4) * t363) * qJD(3) / 0.2e1 + qJD(1) * t320 / 0.2e1;
t242 = Ifges(4,4) * t335;
t178 = Ifges(4,1) * t246 + Ifges(4,5) * t288 - t242;
t506 = t178 / 0.2e1;
t495 = t246 / 0.2e1;
t491 = t288 / 0.2e1;
t486 = pkin(3) * t246;
t485 = pkin(3) * t313;
t477 = mrSges(4,3) * t246;
t476 = mrSges(7,3) * t312;
t475 = mrSges(7,3) * t315;
t474 = Ifges(3,4) * t314;
t473 = Ifges(3,4) * t316;
t472 = Ifges(4,4) * t246;
t460 = t277 * t312;
t459 = t277 * t315;
t100 = -mrSges(7,1) * t161 + mrSges(7,2) * t162;
t168 = mrSges(6,1) * t191 + mrSges(6,3) * t281;
t453 = t100 - t168;
t99 = t313 * t159 + t488 * t172;
t422 = t310 * t446;
t268 = pkin(8) * t422 + qJD(2) * t303;
t442 = t488 * pkin(3);
t438 = Ifges(4,4) * t489;
t437 = Ifges(4,4) * t487;
t433 = Ifges(4,3) * t397 - t575 + t576;
t415 = -t452 / 0.2e1;
t412 = -t448 / 0.2e1;
t94 = -t109 * mrSges(6,1) + mrSges(6,2) * t397;
t80 = t153 * t313 + t150;
t407 = mrSges(3,3) * t424;
t406 = mrSges(4,3) * t419;
t405 = mrSges(4,3) * t417;
t307 = -t442 - pkin(4);
t395 = t316 * t415;
t98 = t159 * t488 - t313 * t172;
t384 = t486 + t462;
t343 = t313 * t350;
t212 = t269 * t488 - t343;
t89 = pkin(4) * t457 - t98;
t67 = t212 * pkin(5) + pkin(11) * t457 + t89;
t339 = t488 * t350;
t211 = t269 * t313 + t339;
t215 = pkin(3) * t404 + t300 + (t308 - t483) * t311;
t325 = -t212 * qJ(5) + t215;
t83 = t211 * t522 + t325;
t36 = -t312 * t83 + t315 * t67;
t37 = t312 * t67 + t315 * t83;
t115 = -mrSges(7,2) * t189 + mrSges(7,3) * t161;
t116 = mrSges(7,1) * t189 - mrSges(7,3) * t162;
t380 = t115 * t315 - t116 * t312;
t88 = qJ(5) * t457 - t99;
t371 = -t211 * t312 + t315 * t457;
t185 = t211 * t315 + t312 * t457;
t361 = -t111 * t488 + t313 * t119 + t159 * t449 + t172 * t418;
t360 = Ifges(3,5) * t398 - Ifges(3,6) * t397;
t358 = pkin(1) * t533 * (mrSges(3,1) * t314 + mrSges(3,2) * t316);
t357 = mrSges(4,1) * t487 + mrSges(4,2) * t489;
t356 = Ifges(4,1) * t489 - t437;
t355 = -Ifges(4,2) * t487 + t438;
t354 = Ifges(4,5) * t489 - Ifges(4,6) * t487;
t353 = t314 * t533 * (Ifges(3,1) * t316 - t474);
t348 = qJD(3) * t356;
t347 = qJD(3) * t354;
t346 = t357 * t457;
t342 = t310 * t411 * (Ifges(3,5) * t316 - Ifges(3,6) * t314);
t341 = qJD(2) * t346;
t340 = (Ifges(3,6) * t311 + (Ifges(3,2) * t316 + t474) * t310) * qJD(1);
t332 = t335 * mrSges(4,3);
t327 = (-mrSges(4,2) * t450 + mrSges(4,3) * t345) * t310;
t326 = (mrSges(4,1) * t450 - mrSges(4,3) * t344) * t310;
t198 = pkin(3) * t328 + t268;
t11 = Ifges(7,4) * t62 + Ifges(7,2) * t63 - Ifges(7,6) * t109;
t6 = -pkin(5) * t110 - t15;
t321 = t23 * mrSges(7,3) * t448 + t6 * (mrSges(7,1) * t312 + mrSges(7,2) * t315) + (Ifges(7,1) * t315 - t470) * t527 + (-Ifges(7,2) * t312 + t469) * t526 + t70 * t412 + t11 * t490 + t315 * t531 + (Ifges(7,5) * t315 - Ifges(7,6) * t312) * t521 - t543 + t553 + t598 * qJD(6) - (t161 * t386 + t162 * t387 + t189 * t385) * qJD(6) / 0.2e1;
t131 = qJD(4) * t339 + t269 * t449 + t313 * t328 - t329 * t488;
t132 = -qJD(4) * t343 + t269 * t418 + t313 * t329 + t328 * t488;
t42 = t132 * pkin(4) + t131 * qJ(5) - t212 * qJD(5) + t198;
t319 = (Ifges(4,4) * t344 + Ifges(4,2) * t345 + Ifges(4,6) * t450) * t310;
t305 = qJ(5) + t485;
t291 = Ifges(3,4) * t424;
t262 = -mrSges(3,2) * t411 + t407;
t257 = t300 + (-pkin(2) - t483) * t311;
t255 = t274 * t445;
t227 = Ifges(3,1) * t425 + Ifges(3,5) * t411 + t291;
t226 = Ifges(3,6) * qJD(2) + t340;
t218 = t277 * pkin(4) + t362;
t217 = mrSges(4,1) * t288 - t477;
t216 = -t288 * mrSges(4,2) - t332;
t207 = -t277 * pkin(5) + t234;
t197 = -mrSges(4,3) * t351 + qJD(1) * t327;
t196 = -mrSges(4,3) * qJD(3) * t364 + qJD(1) * t326;
t177 = -Ifges(4,2) * t335 + Ifges(4,6) * t288 + t472;
t170 = mrSges(5,2) * t281 - mrSges(5,3) * t191;
t155 = qJD(1) * t341 + qJD(3) * t555;
t145 = (Ifges(4,4) * t364 - Ifges(4,2) * t363) * qJD(3) + qJD(1) * t319;
t142 = -qJD(3) * t201 + t389;
t135 = -mrSges(6,2) * t191 - mrSges(6,3) * t324;
t134 = mrSges(5,1) * t191 + mrSges(5,2) * t324;
t133 = t462 + t586;
t120 = t211 * pkin(4) + t325;
t112 = t384 + t586;
t96 = -mrSges(5,2) * t397 - mrSges(5,3) * t110;
t95 = mrSges(5,1) * t397 + mrSges(5,3) * t109;
t85 = t462 + t569;
t78 = qJD(6) * t371 + t132 * t315 - t312 * t423;
t77 = qJD(6) * t185 + t132 * t312 + t315 * t423;
t76 = t384 + t569;
t73 = -pkin(5) * t211 - t88;
t58 = t80 - t484;
t55 = t75 - t484;
t50 = mrSges(5,1) * t110 - mrSges(5,2) * t109;
t49 = -mrSges(6,2) * t110 + mrSges(6,3) * t109;
t34 = t312 * t55 + t315 * t85;
t33 = -t312 * t85 + t315 * t55;
t30 = t312 * t58 + t315 * t76;
t29 = -t312 * t76 + t315 * t58;
t28 = -pkin(4) * t423 + t361;
t27 = t132 * pkin(11) + t42;
t14 = -pkin(5) * t132 - t26;
t13 = -t131 * pkin(5) + t361 - t381;
t4 = -qJD(6) * t37 + t13 * t315 - t27 * t312;
t3 = qJD(6) * t36 + t13 * t312 + t27 * t315;
t5 = [(Ifges(7,4) * t77 + Ifges(7,2) * t78) * t509 + m(3) * (t254 * t274 + t266 * t267) + (-m(3) * t273 + m(4) * t257 - mrSges(3,1) * t311 + mrSges(4,1) * t350 + t269 * mrSges(4,2)) * t255 + (Ifges(7,1) * t77 + Ifges(7,4) * t78) * t507 - (t340 + t226) * t423 / 0.2e1 + (-Ifges(5,4) * t521 + t539) * t211 + (-t129 * t350 - t130 * t269) * mrSges(4,3) + (t254 * t457 + t255 * t458 - t263 * t422 - t266 * t423 - t273 * t398 - t274 * t397) * mrSges(3,3) + t592 * t212 + t77 * t523 + t78 * t524 + t533 * (-Ifges(3,2) * t314 + t473) * t416 + t311 * t360 / 0.2e1 + (Ifges(7,5) * t77 + Ifges(7,6) * t78) * t504 + (-0.2e1 * t358 + t353) * t445 - t335 * (qJD(3) * t311 * t355 + t319) / 0.2e1 - t350 * t145 / 0.2e1 + t228 * (t341 + (mrSges(4,1) * t269 - mrSges(4,2) * t350) * qJD(3)) - t359 * t170 + m(5) * (t179 * t215 + t194 * t198 + t20 * t99 - t359 * t75 + t361 * t74 - t366 * t98) - t361 * t171 + m(4) * (t129 * t201 + t130 * t200 - t141 * t352 + t142 * t180) - t352 * (-t311 * t405 + t327) - t371 * t531 + (-Ifges(7,4) * t371 + Ifges(7,2) * t185) * t526 + (-Ifges(7,1) * t371 + Ifges(7,4) * t185) * t527 + (t1 * t185 + t2 * t371 - t23 * t77 + t24 * t78) * mrSges(7,3) + (-Ifges(7,5) * t371 + Ifges(7,6) * t185) * t521 + t6 * (-mrSges(7,1) * t185 - mrSges(7,2) * t371) + t269 * t511 + t329 * t506 + (t311 * t348 + t320) * t495 + (t311 * t347 + (Ifges(4,5) * t344 + Ifges(4,6) * t345 + Ifges(4,3) * t450) * t310) * t491 + t413 * t547 + t180 * (-t311 * t406 + t326) + t342 * t584 + (Ifges(4,1) * t269 - Ifges(4,4) * t350) * t588 + t36 * t38 + t37 * t39 + (t534 - t593) * t131 + t594 * t132 + (-m(3) * t263 - t548) * t268 + t73 * t25 + t53 * (-mrSges(7,1) * t78 + mrSges(7,2) * t77) + (Ifges(4,4) * t269 - Ifges(4,2) * t350) * t589 + t88 * t93 + t89 * t94 + t98 * t95 + t99 * t96 + t14 * t100 - (t433 + t553) * t457 / 0.2e1 + t3 * t115 + t4 * t116 + t120 * t49 + (Ifges(4,5) * t269 - Ifges(4,6) * t350 + t581 * t211 + (-Ifges(4,3) - t583) * t457) * t375 + t42 * t135 - t311 * t574 + (-t576 / 0.2e1 + t575 / 0.2e1 - Ifges(6,5) * t518 - Ifges(5,6) * t519 - Ifges(6,4) * t520 - Ifges(5,5) * t521 + t543 + t545) * t457 + t26 * t168 + t28 * t169 + t185 * t11 / 0.2e1 + (Ifges(6,4) * t500 + Ifges(5,5) * t499 + Ifges(6,5) * t502 + Ifges(5,6) * t503 + t493 * t583 + t542) * t423 + t198 * t134 + t200 * t196 + t201 * t197 + ((Ifges(3,5) * t311 + (Ifges(3,1) * t314 + t473) * t310) * t392 + t227 * t413) * t316 + t215 * t50 + t141 * t216 + t142 * t217 + t257 * t155 + m(6) * (t120 * t35 + t15 * t88 + t18 * t89 + t26 * t72 + t28 * t71 + t42 * t87) + m(7) * (t1 * t37 + t14 * t53 + t2 * t36 + t23 * t4 + t24 * t3 + t6 * t73) + t267 * t262 - t328 * t177 / 0.2e1; (t407 - t262) * t263 + (t335 * (Ifges(4,6) * t314 + t316 * t355) + t314 * t226) * t452 / 0.2e1 + (qJD(6) * t70 + t11) * t459 / 0.2e1 + (t291 + t227) * t395 + (Ifges(7,4) * t369 - Ifges(7,2) * t368) * t509 + (Ifges(7,4) * t205 + Ifges(7,2) * t204) * t510 + t592 * t278 + (Ifges(7,5) * t369 - Ifges(7,6) * t368) * t504 + (Ifges(7,5) * t205 + Ifges(7,6) * t204) * t505 + (-t228 * t346 - t342 / 0.2e1 + (t358 - t353 / 0.2e1) * qJD(1)) * qJD(1) + (-(t314 * mrSges(4,1) - mrSges(4,3) * t429) * t452 - t406) * t180 - t574 + t456 * t523 + t454 * t524 + t204 * t525 + t460 * t531 + (t94 - t95) * t233 + (t96 - t93) * t234 - t178 * t379 / 0.2e1 + t360 + (Ifges(7,1) * t369 - Ifges(7,4) * t368) * t507 + (Ifges(7,1) * t205 + Ifges(7,4) * t204) * t508 - t355 * t330 / 0.2e1 + t228 * t357 * qJD(3) + (Ifges(4,5) * t487 + Ifges(4,6) * t489) * t375 + t489 * t145 / 0.2e1 + (t288 * (Ifges(4,3) * t314 + t316 * t354) + t246 * (Ifges(4,5) * t314 + t316 * t356) + t547) * t415 + (t179 * t308 + t194 * t556 + t20 * t234 + t233 * t366 + t561 * t75 - t562 * t74) * m(5) - (-t405 - (-t314 * mrSges(4,2) - mrSges(4,3) * t428) * t452) * t352 + (-t180 * t202 + t352 * t203 - t255 * pkin(2) + ((-t180 * t489 + t352 * t487) * qJD(3) + t554) * pkin(9)) * m(4) - t196 * t441 + t419 * t506 + t487 * t511 + t348 * t495 + t347 * t491 + (-t417 / 0.2e1 + t378 / 0.2e1) * t177 + (Ifges(4,1) * t487 + t438) * t588 + t594 * t223 + (t538 - t596) * t230 + (-t216 * t417 - t217 * t419) * pkin(9) + (t535 - t597) * t231 + (t534 - t597) * t222 + (t408 + t548) * t266 + (Ifges(4,2) * t489 + t437) * t589 + t573 * (-t231 / 0.2e1 - t222 / 0.2e1) + t554 * mrSges(4,3) + t556 * t134 + (mrSges(7,1) * t558 - mrSges(7,2) * t557) * t53 + (t1 * t459 - t2 * t460 + t23 * t557 - t24 * t558) * mrSges(7,3) + t560 * t135 + t561 * t170 + t562 * t171 + t563 * t169 + t564 * t168 + (-t15 * t234 + t18 * t233 + t218 * t35 + t560 * t87 + t563 * t71 + t564 * t72) * m(6) + t570 * t100 + t138 * t38 + t139 * t39 - pkin(2) * t155 + t579 * t116 + t580 * t115 + (t1 * t139 + t138 * t2 + t207 * t6 + t23 * t579 + t24 * t580 + t53 * t570) * m(7) + (t539 + (-Ifges(5,4) + t385) * t521 + t386 * t526 + t387 * t527 - t388 * t6 + t412 * t69 + t581 * t375) * t277 + (Ifges(6,4) * t499 + Ifges(5,5) * t500 + Ifges(6,5) * t503 - Ifges(3,2) * t395 + Ifges(5,6) * t502 + t492 * t583 - t542) * t425 - t205 * t70 / 0.2e1 + t207 * t25 - t203 * t216 - t202 * t217 + t218 * t49 + t197 * t443 + t308 * t50 + (-mrSges(4,1) * t489 + mrSges(4,2) * t487 - mrSges(3,1)) * t255; (t23 * t476 - t24 * t475 + t538 + t591) * t324 + t433 - t246 * (-Ifges(4,1) * t335 - t472) / 0.2e1 + t335 * (-Ifges(4,2) * t246 - t242) / 0.2e1 + (-t216 - t332) * t180 + t603 * t170 - t2 * t475 - t1 * t476 - t24 * mrSges(7,3) * t447 + ((-t488 * t366 + t20 * t313 + (t313 * t74 + t488 * t75) * qJD(4)) * pkin(3) - t194 * t486 - t74 * t80 - t75 * t81) * m(5) - (t477 + t217) * t352 + t335 * t506 + t177 * t495 + t96 * t485 - t288 * (-Ifges(4,5) * t335 - Ifges(4,6) * t246) / 0.2e1 - t545 + (m(7) * (t23 * t315 + t24 * t312) + m(6) * t71 + t315 * t116 + t312 * t115 - t559) * pkin(3) * t449 + (-t535 + t593) * t191 + (t115 * t447 - t116 * t448 + t599) * (-pkin(11) + t307) + t321 - t228 * t555 + t559 * t80 - t30 * t115 - t29 * t116 - t112 * t135 + (-t112 * t87 - t15 * t305 + t18 * t307 + t571 * t72 - t71 * t80) * m(6) + t571 * t168 + (-t23 * t29 - t24 * t30 + t305 * t6 + t53 * t572) * m(7) + t572 * t100 + t577 * t305 + t95 * t442 + t307 * t94 - t134 * t486; (-(qJD(6) * t115 + t38) * t522 + (-t2 - t461) * mrSges(7,3)) * t315 + (-t1 * mrSges(7,3) - (-qJD(6) * t116 + t39) * t522) * t312 + (-t168 + t170) * t74 + (-t187 / 0.2e1 - t188 / 0.2e1 + t466 / 0.2e1 + t468 / 0.2e1 + t467 / 0.2e1 + t68 / 0.2e1 + t126 / 0.2e1 + (Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t281 - t540 + t597) * t191 - t522 * t336 + ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t324 + (t590 - Ifges(5,6) / 0.2e1) * t281 + t383 * mrSges(7,3) + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t191 + t544 + t591) * t324 + t453 * qJD(5) + t559 * t75 + t577 * qJ(5) - pkin(4) * t94 + t373 * t100 + t321 - t34 * t115 - t33 * t116 - t133 * t135 + (t6 * qJ(5) - t23 * t33 - t24 * t34 + t53 * t595) * m(7) + (-pkin(4) * t18 - qJ(5) * t15 - t133 * t87 + t565 * t72 - t71 * t75) * m(6); t453 * t281 + t380 * qJD(6) + (t135 + t380) * t324 - m(7) * (-t281 * t53 + t324 * t383) + t94 + (-t281 * t72 + t324 * t87 + t18) * m(6) + t599; -t53 * (mrSges(7,1) * t162 + mrSges(7,2) * t161) + (Ifges(7,1) * t161 - t471) * t508 + t69 * t507 + (Ifges(7,5) * t161 - Ifges(7,6) * t162) * t505 - t23 * t115 + t24 * t116 + (t161 * t23 + t162 * t24) * mrSges(7,3) + t10 + (-Ifges(7,2) * t162 + t160 + t70) * t510 + t549;];
tauc  = t5(:);
