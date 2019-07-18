% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:15
% EndTime: 2019-07-18 17:17:35
% DurationCPUTime: 8.36s
% Computational Cost: add. (14946->485), mult. (33604->668), div. (0->0), fcn. (37325->8), ass. (0->296)
t313 = sin(qJ(5));
t315 = cos(qJ(4));
t314 = sin(qJ(4));
t534 = cos(qJ(5));
t406 = t534 * t314;
t360 = t313 * t315 + t406;
t446 = t313 * t314;
t561 = t534 * t315 - t446;
t226 = -mrSges(6,1) * t561 + mrSges(6,2) * t360;
t517 = Ifges(5,4) * t314;
t297 = Ifges(5,2) * t315 + t517;
t309 = Ifges(5,4) * t315;
t298 = Ifges(5,1) * t314 + t309;
t528 = pkin(3) * t314;
t587 = -t315 / 0.2e1;
t588 = t314 / 0.2e1;
t599 = -t226 * t528 + t297 * t588 + t298 * t587;
t511 = Ifges(5,2) * t314;
t376 = t511 - t309;
t377 = Ifges(5,1) * t315 - t517;
t551 = t376 * t587 + t377 * t588 - t599;
t535 = cos(qJ(3));
t427 = t535 * pkin(1);
t307 = -t427 - pkin(2);
t527 = t315 * pkin(3);
t294 = t307 - t527;
t308 = -pkin(2) - t527;
t424 = t528 / 0.2e1;
t519 = mrSges(5,2) * t315;
t520 = mrSges(5,1) * t314;
t296 = t519 + t520;
t452 = t307 * t296;
t530 = pkin(2) * t296;
t598 = -t530 / 0.2e1 + t452 / 0.2e1 + m(6) * (t294 + t308) * t424 + t551;
t532 = sin(qJ(3));
t533 = sin(qJ(2));
t536 = cos(qJ(2));
t289 = -t532 * t536 - t533 * t535;
t538 = -t289 / 0.2e1;
t288 = t532 * t533 - t535 * t536;
t204 = t561 * t288;
t540 = -t204 / 0.2e1;
t202 = t360 * t288;
t542 = t202 / 0.2e1;
t565 = Ifges(6,5) * t540 + Ifges(6,6) * t542;
t408 = Ifges(6,3) * t538 + t565;
t235 = -t289 * pkin(2) + t288 * pkin(5);
t284 = t289 * pkin(3);
t157 = t315 * t235 - t284;
t90 = t157 * t534 - t235 * t446;
t593 = t90 / 0.2e1;
t595 = -mrSges(6,2) / 0.2e1;
t91 = t313 * t157 + t235 * t406;
t567 = mrSges(6,1) * t593 + t595 * t91;
t597 = -t408 - t567;
t426 = t533 * pkin(1);
t217 = t426 + t235;
t148 = t315 * t217 - t284;
t73 = t313 * t148 + t217 * t406;
t594 = -t73 / 0.2e1;
t71 = t148 * t534 - t217 * t446;
t568 = t71 * mrSges(6,1) / 0.2e1 + mrSges(6,2) * t594;
t596 = -t408 - t568;
t203 = t360 * t289;
t503 = t203 * mrSges(6,3);
t140 = -mrSges(6,2) * t288 + t503;
t592 = -t140 / 0.2e1;
t425 = t532 * pkin(1);
t383 = t425 + pkin(5);
t354 = t383 * t534;
t364 = t315 * t383;
t247 = -t313 * t364 - t314 * t354;
t591 = t247 / 0.2e1;
t272 = t561 * pkin(5);
t590 = -t272 / 0.2e1;
t589 = t294 / 0.2e1;
t311 = t314 ^ 2;
t312 = t315 ^ 2;
t562 = t311 + t312;
t586 = t427 * t562;
t278 = Ifges(6,4) * t561;
t231 = Ifges(6,1) * t360 + t278;
t585 = -Ifges(6,2) * t360 + t231 + t278;
t493 = t315 * mrSges(5,1);
t494 = t314 * mrSges(5,2);
t584 = t493 / 0.2e1 - t494 / 0.2e1;
t583 = -pkin(3) / 0.2e1;
t381 = t535 * t534;
t407 = t313 * t535;
t254 = (-t314 * t381 - t315 * t407) * pkin(1);
t581 = t254 / 0.2e1;
t580 = t561 / 0.2e1;
t578 = t360 / 0.2e1;
t496 = t360 * mrSges(6,3);
t411 = -t496 / 0.2e1;
t412 = t496 / 0.2e1;
t365 = t314 * t383;
t246 = t313 * t365 - t315 * t354;
t577 = -t272 + t246;
t515 = Ifges(6,4) * t360;
t229 = Ifges(6,2) * t561 + t515;
t230 = Ifges(6,1) * t561 - t515;
t384 = -t360 * t229 / 0.2e1 + t230 * t578 + t585 * t580;
t225 = mrSges(6,1) * t360 + mrSges(6,2) * t561;
t449 = t308 * t225;
t454 = t294 * t225;
t576 = t384 + t246 * t411 + t272 * t412 + t449 / 0.2e1 + t454 / 0.2e1;
t227 = Ifges(6,5) * t561 - Ifges(6,6) * t360;
t273 = t360 * pkin(5);
t563 = -t272 * mrSges(6,1) + t273 * mrSges(6,2);
t366 = t227 + t563;
t574 = qJD(5) * t366;
t564 = t246 * mrSges(6,1) - t247 * mrSges(6,2);
t367 = t227 + t564;
t573 = qJD(5) * t367;
t255 = (-t314 * t407 + t315 * t381) * pkin(1);
t571 = (-t254 * t360 + t255 * t561) * mrSges(6,3);
t201 = t561 * t289;
t101 = -mrSges(6,1) * t201 + mrSges(6,2) * t203;
t504 = t201 * mrSges(6,3);
t142 = mrSges(6,1) * t288 + t504;
t570 = t308 * t101 / 0.2e1 + t273 * t592 + t142 * t590;
t569 = t101 * t589 + t140 * t591 + t246 * t142 / 0.2e1;
t513 = Ifges(5,5) * t315;
t375 = -Ifges(5,6) * t314 + t513;
t358 = t288 * t375;
t438 = mrSges(6,1) * t581 + t255 * t595;
t559 = mrSges(5,3) * t538 * t562;
t434 = t294 * t528;
t244 = t288 * t434;
t546 = m(6) / 0.2e1;
t428 = t536 * pkin(1);
t218 = t288 * pkin(2) + t289 * pkin(5) - t428;
t480 = t218 * t315;
t149 = pkin(3) * t288 + t480;
t72 = t149 * t534 - t218 * t446;
t74 = t313 * t149 + t218 * t406;
t558 = (-t246 * t91 + t247 * t90 + t254 * t72 + t255 * t74 - t244) * t546 + t142 * t581 + t255 * t140 / 0.2e1;
t295 = -t493 + t494;
t556 = -t535 * mrSges(4,2) + (-mrSges(4,1) + t226 + t295) * t532;
t416 = -t503 / 0.2e1;
t418 = -t504 / 0.2e1;
t555 = t246 * t418 + t247 * t416 + t569;
t544 = m(6) * pkin(3);
t436 = t544 / 0.2e1;
t554 = (t313 * t91 + t534 * t90) * t436 + t567 + t584 * t235;
t553 = (t313 * t73 + t534 * t71) * t436 + t568 + t584 * t217;
t552 = -t272 * t418 - t273 * t416 + t570;
t457 = t289 * t315;
t458 = t289 * t314;
t510 = Ifges(5,6) * t288;
t514 = Ifges(5,5) * t288;
t549 = t315 * (-t289 * t377 + t514) / 0.4e1 - t314 * (t289 * t376 + t510) / 0.4e1 + t358 / 0.4e1 + (t298 / 0.2e1 - t376 / 0.4e1) * t458 + (t297 / 0.2e1 - t377 / 0.4e1 + t226 * t583) * t457;
t195 = Ifges(6,4) * t203;
t104 = Ifges(6,2) * t201 + t195;
t516 = Ifges(6,4) * t201;
t105 = Ifges(6,1) * t203 + t516;
t435 = pkin(3) * t458;
t81 = Ifges(6,2) * t203 + t288 * Ifges(6,6) - t516;
t83 = -Ifges(6,1) * t201 + t288 * Ifges(6,5) + t195;
t548 = -t225 * t435 / 0.2e1 + t288 * t227 / 0.4e1 + t585 * t203 / 0.4e1 + (t104 + t83) * t561 / 0.4e1 + (t105 / 0.4e1 - t81 / 0.4e1) * t360 + (-t230 / 0.4e1 + t229 / 0.4e1) * t201;
t547 = -m(6) / 0.2e1;
t545 = pkin(2) / 0.2e1;
t543 = -t201 / 0.2e1;
t541 = t203 / 0.2e1;
t220 = -mrSges(5,2) * t288 + mrSges(5,3) * t458;
t539 = t220 / 0.2e1;
t537 = t315 / 0.2e1;
t214 = t296 * t288;
t531 = pkin(2) * t214;
t529 = pkin(3) * t313;
t525 = t72 * mrSges(6,2);
t523 = t74 * mrSges(6,1);
t509 = Ifges(5,6) * t289;
t106 = t360 * t218;
t506 = t106 * mrSges(6,1);
t107 = t561 * t218;
t505 = t107 * mrSges(6,2);
t498 = t561 * mrSges(6,3);
t316 = pkin(3) ^ 2;
t448 = t311 * t316;
t232 = t289 * t288 * t448;
t103 = -mrSges(6,1) * t203 - mrSges(6,2) * t201;
t432 = t103 * t528;
t321 = -t432 + Ifges(4,4) * t288 - t358 + (-Ifges(5,3) - Ifges(6,3) - Ifges(4,2) + Ifges(4,1) + t312 * Ifges(5,1) / 0.2e1 + (-t309 + t511 / 0.2e1) * t314) * t289 + t565;
t102 = -mrSges(6,1) * t202 - mrSges(6,2) * t204;
t150 = t288 * t376 - t509;
t151 = -Ifges(5,5) * t289 - t288 * t377;
t322 = Ifges(6,5) * t201 / 0.2e1 - Ifges(6,6) * t203 / 0.2e1 + t151 * t587 + (t513 / 0.2e1 - Ifges(4,4)) * t289 + (-pkin(3) * t102 + t150 / 0.2e1 - t509 / 0.2e1) * t314;
t139 = mrSges(6,2) * t289 + mrSges(6,3) * t202;
t141 = -mrSges(6,1) * t289 + mrSges(6,3) * t204;
t80 = -Ifges(6,4) * t204 + Ifges(6,2) * t202 - t289 * Ifges(6,6);
t82 = -Ifges(6,1) * t204 + Ifges(6,4) * t202 - t289 * Ifges(6,5);
t328 = t74 * t139 + t72 * t141 + t83 * t540 + t80 * t541 + t81 * t542 + t82 * t543 - (-mrSges(4,1) * t289 - mrSges(4,2) * t288) * t428;
t222 = mrSges(5,1) * t288 + mrSges(5,3) * t457;
t445 = t314 * t220;
t368 = t315 * t222 + t445;
t460 = t288 * t314;
t219 = mrSges(5,2) * t289 + mrSges(5,3) * t460;
t459 = t288 * t315;
t221 = -mrSges(5,1) * t289 + mrSges(5,3) * t459;
t369 = t314 * t219 + t315 * t221;
t400 = m(5) * t562;
t5 = (-mrSges(4,2) * t426 + t322) * t289 + (mrSges(4,1) * t426 + t321) * t288 + m(6) * (t71 * t72 + t73 * t74 + t232) + t328 + t73 * t140 + t71 * t142 + t369 * t218 + (t218 * t400 + t368) * t217 + (-t533 ^ 2 + t536 ^ 2) * Ifges(3,4) + (-m(4) * pkin(1) ^ 2 + Ifges(3,1) - Ifges(3,2)) * t536 * t533;
t492 = t5 * qJD(1);
t7 = t91 * t140 + t90 * t142 + m(6) * (t72 * t90 + t74 * t91 + t232) + t368 * t235 + (t235 * t400 + t369) * t218 + t322 * t289 + t321 * t288 + t328;
t491 = t7 * qJD(1);
t490 = t90 * t360;
t489 = t91 * t561;
t488 = -t106 + t74;
t487 = -t107 + t72;
t440 = Ifges(6,5) * t203 + Ifges(6,6) * t201;
t327 = (t83 / 0.2e1 + t104 / 0.2e1) * t203 + (-t105 / 0.2e1 + t81 / 0.2e1 + t74 * mrSges(6,3)) * t201 - t72 * t503 + t288 * t440 / 0.2e1;
t11 = -t106 * t142 + t107 * t140 + m(6) * (-t106 * t72 + t107 * t74) + (t220 * t315 - t222 * t314) * t218 + ((-Ifges(5,4) * t457 - pkin(3) * t103 + t510) * t315 + (Ifges(5,4) * t458 + t514 - pkin(3) * t101 + (m(6) * t316 - Ifges(5,1) + Ifges(5,2)) * t457) * t314) * t289 + t327;
t486 = t11 * qJD(1);
t14 = -t101 * t435 + t72 * t140 - t74 * t142 + t327;
t485 = t14 * qJD(1);
t479 = t246 * t139;
t476 = t247 * t141;
t472 = t272 * t139;
t469 = t273 * t141;
t455 = t294 * t102;
t453 = t307 * t214;
t450 = t308 * t102;
t447 = t313 * t142;
t444 = t314 * t221;
t437 = mrSges(6,3) * t529;
t433 = t308 * t528;
t431 = -t532 / 0.2e1;
t430 = t532 / 0.2e1;
t429 = pkin(3) * t534;
t423 = pkin(5) * t537;
t417 = t504 / 0.2e1;
t415 = t503 / 0.2e1;
t414 = -t498 / 0.2e1;
t413 = t498 / 0.2e1;
t398 = t201 * t437;
t396 = mrSges(6,3) * t429;
t395 = t288 * t433;
t391 = t532 * t528;
t52 = t74 * t412;
t390 = t429 / 0.2e1;
t389 = -t427 / 0.2e1;
t380 = t535 * t539;
t379 = t103 * t430;
t378 = t203 * t396;
t373 = t314 * t389;
t324 = Ifges(5,3) * t538 - Ifges(5,5) * t459 / 0.2e1 + Ifges(5,6) * t460 / 0.2e1 + t139 * t529 / 0.2e1 + t141 * t390 + t408;
t317 = -t548 - t549 - t432 / 0.2e1 + t52 + t324 - t106 * t412 + t72 * t413 + t107 * t414;
t353 = t364 / 0.2e1;
t356 = t295 * t289;
t318 = ((-t294 * t527 - t448) * t289 + t488 * t247 + t487 * t246) * t547 - t307 * t356 / 0.2e1 + t365 * t539 + t222 * t353 + t383 * t559;
t1 = t246 * t417 + t247 * t415 + t317 + t318 + t553 - t569;
t333 = t384 + t454;
t21 = m(6) * t434 + t333 + t452 + t551;
t372 = -t1 * qJD(1) + t21 * qJD(2);
t349 = t562 * t535;
t30 = t571 + (mrSges(5,3) * t349 + t556) * pkin(1) + m(6) * (-t246 * t255 + t247 * t254 + t294 * t425) + m(5) * (t307 * t425 + t383 * t586);
t338 = m(6) * (t272 * t73 - t273 * t71 - t395);
t343 = t296 * t431;
t8 = -(t545 + t307 / 0.2e1) * t214 + (t273 / 0.2e1 + t591) * t141 + (t590 - t246 / 0.2e1) * t139 + (-t308 / 0.2e1 + t589) * t102 - t338 / 0.2e1 + (-(t593 - t71 / 0.2e1) * t360 + (t91 / 0.2e1 + t594) * t561) * mrSges(6,3) + (t379 + (t219 * t430 + t380) * t315 + (-t535 * t222 / 0.2e1 + t221 * t431) * t314 + (t391 * t547 + t343) * t289) * pkin(1) + t558;
t371 = t8 * qJD(1) + t30 * qJD(2);
t339 = t74 * t411 + t548;
t329 = t52 + t339;
t326 = t329 + t555;
t10 = t326 + t596;
t370 = t10 * qJD(1) + qJD(2) * t333;
t355 = t221 * t365;
t323 = (t254 * t534 + t255 * t313) * t436 + mrSges(5,1) * t373 + t389 * t519 + t438;
t16 = t577 * t411 + t323 - t576 - t598;
t332 = t384 + t449;
t27 = m(6) * t433 + t332 - t530 + t551;
t320 = ((-t308 * t527 - t448) * t289 - t488 * t273 - t487 * t272) * t547 + t356 * t545 + t222 * t423 + (t445 / 0.2e1 + t559) * pkin(5);
t3 = -t272 * t417 - t273 * t415 + t317 + t320 + t554 - t570;
t347 = -t3 * qJD(1) - t16 * qJD(2) + t27 * qJD(3);
t325 = t329 + t552;
t13 = t325 + t597;
t340 = t577 * t412 + t576;
t19 = t340 - t438;
t345 = t13 * qJD(1) + t19 * qJD(2) + qJD(3) * t332;
t344 = -t360 * t437 - t396 * t561 + t227;
t18 = (-t107 / 0.2e1 + t72 / 0.2e1) * mrSges(6,2) + (-t106 / 0.2e1 + t74 / 0.2e1) * mrSges(6,1) + (t534 * t592 + t447 / 0.2e1 + (t313 * t543 + t534 * t541) * mrSges(6,3)) * pkin(3);
t293 = (mrSges(6,1) * t313 + mrSges(6,2) * t534) * pkin(3);
t337 = qJD(1) * t18 + qJD(4) * t293;
t331 = t229 * t542 + t231 * t540 + t151 * t588 + t150 * t537 + Ifges(4,6) * t289 + t80 * t580 + t82 * t578 + (Ifges(5,5) * t314 + Ifges(6,5) * t360 + Ifges(5,6) * t315 + Ifges(6,6) * t561) * t538 + (-Ifges(4,5) + t599) * t288;
t319 = t103 * t424 - t106 * t411 + t107 * t413 + t72 * t414 + t324 + t339 + t549;
t285 = t293 * qJD(5);
t20 = t340 + t438;
t17 = t323 + t340 + t598;
t15 = -t523 / 0.2e1 - t525 / 0.2e1 + t140 * t390 + t398 / 0.2e1 + t447 * t583 - t378 / 0.2e1 - t506 / 0.2e1 - t505 / 0.2e1 + t440;
t12 = t325 - t597;
t9 = t326 - t596;
t6 = t331 - t355 / 0.2e1 + t222 * t373 - t453 / 0.2e1 + t450 / 0.2e1 + t455 / 0.2e1 - t469 / 0.2e1 + t472 / 0.2e1 + t476 / 0.2e1 - t479 / 0.2e1 + t531 / 0.2e1 + (t489 / 0.2e1 - t490 / 0.2e1) * mrSges(6,3) + t338 / 0.2e1 + t73 * t413 + t71 * t411 - pkin(5) * t444 / 0.2e1 + (t353 + t423) * t219 + (t315 * t380 + t379 + (-t391 * t546 + t343) * t289) * pkin(1) + t558;
t4 = t319 - t320 + t552 + t554;
t2 = -t318 + t319 + t553 + t555;
t22 = [qJD(2) * t5 + qJD(3) * t7 + qJD(4) * t11 + qJD(5) * t14, t492 + (t331 + Ifges(3,5) * t536 - Ifges(3,6) * t533 + (-t360 * t71 + t561 * t73) * mrSges(6,3) + (t288 * t535 + t289 * t532) * pkin(1) * mrSges(4,3) - t453 + t455 + t476 - t479 + m(6) * (-t246 * t73 + t247 * t71 - t244) - t355 + t219 * t364) * qJD(2) + t6 * qJD(3) + t2 * qJD(4) + t9 * qJD(5), t491 + t6 * qJD(2) + (t331 + t450 - t469 + t472 + t531 + (t315 * t219 - t444) * pkin(5) + (t489 - t490) * mrSges(6,3) + m(6) * (t272 * t91 - t273 * t90 - t395)) * qJD(3) + t4 * qJD(4) + t12 * qJD(5), t486 + t2 * qJD(2) + t4 * qJD(3) + (-mrSges(5,2) * t480 - t218 * t520 + Ifges(5,5) * t458 + Ifges(5,6) * t457 - t506 - t505 + (-t106 * t534 + t107 * t313) * t544 + t398 - t378 + t440) * qJD(4) + t15 * qJD(5), t485 + t9 * qJD(2) + t12 * qJD(3) + t15 * qJD(4) + (t440 - t523 - t525) * qJD(5); qJD(3) * t8 - qJD(4) * t1 + qJD(5) * t10 - t492, qJD(3) * t30 + qJD(4) * t21 + qJD(5) * t333, t17 * qJD(4) + t20 * qJD(5) + t371 + (m(6) * (-t273 * t254 + t272 * t255 + t308 * t425) + mrSges(5,3) * t586 + t571 + (m(5) * (-pkin(2) * t532 + pkin(5) * t349) + t556) * pkin(1)) * qJD(3), t17 * qJD(3) + (m(6) * (t246 * t429 + t247 * t529) + (-mrSges(5,1) * t383 + Ifges(5,5)) * t315 + (mrSges(5,2) * t383 - Ifges(5,6)) * t314 + t344 + t564) * qJD(4) + t573 + t372, t20 * qJD(3) + qJD(4) * t367 + t370 + t573; -qJD(2) * t8 - qJD(4) * t3 + qJD(5) * t13 - t491, -qJD(4) * t16 + qJD(5) * t19 - t371, qJD(4) * t27 + qJD(5) * t332, (m(6) * (-t272 * t429 - t273 * t529) + t295 * pkin(5) + t344 + t375 + t563) * qJD(4) + t574 + t347, qJD(4) * t366 + t345 + t574; qJD(2) * t1 + qJD(3) * t3 - qJD(5) * t18 - t486, qJD(3) * t16 - t372, -t347, -t285, -t285 - t337; -qJD(2) * t10 - qJD(3) * t13 + qJD(4) * t18 - t485, -qJD(3) * t19 - t370, -t345, t337, 0;];
Cq  = t22;
