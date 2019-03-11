% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:11
% EndTime: 2019-03-09 05:40:33
% DurationCPUTime: 49.29s
% Computational Cost: add. (37955->1019), mult. (121121->1393), div. (0->0), fcn. (105526->18), ass. (0->458)
t385 = sin(pkin(7));
t386 = sin(pkin(6));
t384 = sin(pkin(12));
t391 = sin(qJ(3));
t522 = cos(pkin(12));
t523 = cos(pkin(7));
t453 = t523 * t522;
t546 = cos(qJ(3));
t421 = -t384 * t391 + t546 * t453;
t524 = cos(pkin(6));
t460 = t524 * t546;
t399 = t385 * t460 + t386 * t421;
t282 = t399 * qJD(1);
t277 = qJD(4) - t282;
t475 = t391 * t523;
t413 = t386 * (-t384 * t475 + t522 * t546);
t313 = qJD(1) * t413;
t483 = qJD(3) * t546;
t463 = t385 * t483;
t683 = t463 - t313;
t544 = pkin(9) * t384;
t317 = (-pkin(2) * t522 - t385 * t544 - pkin(1)) * t386;
t306 = qJD(1) * t317 + qJD(2);
t454 = t524 * t522;
t374 = pkin(1) * t454;
t367 = qJD(1) * t374;
t491 = t524 * pkin(2);
t515 = t384 * t386;
t414 = t491 + (-pkin(9) * t523 - qJ(2)) * t515;
t278 = qJD(1) * t414 + t367;
t473 = t523 * t278;
t478 = t384 * t524;
t373 = pkin(1) * t478;
t476 = t386 * t522;
t458 = qJD(1) * t476;
t328 = qJ(2) * t458 + qJD(1) * t373;
t474 = t524 * t385;
t411 = (t386 * t453 + t474) * pkin(9);
t271 = qJD(1) * t411 + t328;
t490 = t546 * t271;
t189 = t490 + (t306 * t385 + t473) * t391;
t390 = sin(qJ(4));
t393 = cos(qJ(4));
t687 = -t390 * qJD(5) - t189 + t277 * (pkin(4) * t390 - qJ(5) * t393);
t513 = t385 * t391;
t343 = t390 * t513 - t393 * t523;
t509 = qJD(1) * t386;
t489 = t384 * t509;
t465 = t385 * t489;
t612 = -qJD(4) * t343 - t390 * t465 + t393 * t683;
t459 = t523 * t546;
t412 = t386 * (t384 * t459 + t391 * t522);
t312 = qJD(1) * t412;
t508 = qJD(3) * t391;
t487 = t385 * t508;
t686 = -t487 + t312;
t438 = t278 * t459;
t494 = t385 * t546;
t188 = -t391 * t271 + t306 * t494 + t438;
t420 = t384 * t546 + t391 * t453;
t295 = t386 * t420 + t391 * t474;
t285 = t295 * qJD(1);
t227 = pkin(3) * t285 - pkin(10) * t282;
t132 = t393 * t188 + t390 * t227;
t118 = qJ(5) * t285 + t132;
t383 = sin(pkin(13));
t387 = cos(pkin(13));
t506 = qJD(4) * t390;
t502 = pkin(10) * t506;
t631 = t687 * t387 + (t118 + t502) * t383;
t685 = -t387 * t118 + t383 * t687;
t378 = pkin(5) * t387 + pkin(4);
t449 = -t387 * mrSges(6,1) + t383 * mrSges(6,2);
t684 = -m(6) * pkin(4) - m(7) * t378 + t449;
t545 = sin(qJ(1));
t547 = cos(qJ(1));
t426 = t545 * t384 - t454 * t547;
t493 = t386 * t547;
t682 = t385 * t493 + t426 * t523;
t335 = t478 * t547 + t522 * t545;
t249 = -t335 * t546 + t391 * t682;
t477 = t386 * t523;
t403 = t426 * t385 - t547 * t477;
t681 = t249 * t393 - t390 * t403;
t680 = t249 * t390 + t393 * t403;
t511 = t387 * t393;
t216 = t282 * t511 + t383 * t285;
t519 = t282 * t390;
t679 = -pkin(5) * t519 + pkin(11) * t216 + (pkin(5) * t390 - pkin(11) * t511) * qJD(4) + t631;
t516 = t383 * t393;
t215 = -t282 * t516 + t387 * t285;
t512 = t387 * t390;
t678 = -pkin(11) * t215 + (-pkin(10) * t512 - pkin(11) * t516) * qJD(4) + t685;
t614 = -t383 * t612 - t686 * t387;
t613 = -t686 * t383 + t387 * t612;
t382 = pkin(13) + qJ(6);
t379 = sin(t382);
t380 = cos(t382);
t656 = -t380 * mrSges(7,1) + t379 * mrSges(7,2);
t608 = -t656 - t684;
t541 = pkin(11) + qJ(5);
t607 = -m(6) * qJ(5) - m(7) * t541 - mrSges(6,3) - mrSges(7,3);
t419 = t385 * t476 - t523 * t524;
t408 = -qJD(1) * t419 + qJD(3);
t163 = -pkin(3) * t408 - t188;
t235 = t390 * t285 - t393 * t408;
t236 = t393 * t285 + t390 * t408;
t105 = t235 * pkin(4) - t236 * qJ(5) + t163;
t228 = -t278 * t385 + t523 * t306;
t161 = -pkin(3) * t282 - pkin(10) * t285 + t228;
t462 = t391 * t473;
t164 = pkin(10) * t408 + t306 * t513 + t462 + t490;
t90 = t390 * t161 + t164 * t393;
t86 = t277 * qJ(5) + t90;
t49 = t387 * t105 - t383 * t86;
t50 = t383 * t105 + t387 * t86;
t677 = t49 * mrSges(6,1) - t50 * mrSges(6,2);
t193 = t236 * t387 + t277 * t383;
t33 = pkin(5) * t235 - pkin(11) * t193 + t49;
t389 = sin(qJ(6));
t471 = -t236 * t383 + t387 * t277;
t39 = pkin(11) * t471 + t50;
t392 = cos(qJ(6));
t12 = t33 * t392 - t389 * t39;
t13 = t33 * t389 + t39 * t392;
t676 = t163 * mrSges(5,1) + t12 * mrSges(7,1) - t13 * mrSges(7,2);
t530 = t236 * Ifges(5,4);
t664 = t277 * Ifges(5,6);
t143 = -t235 * Ifges(5,2) + t530 + t664;
t127 = t193 * t392 + t389 * t471;
t234 = qJD(6) + t235;
t655 = -t193 * t389 + t392 * t471;
t666 = t193 * Ifges(6,5);
t669 = Ifges(6,6) * t471;
t643 = t127 * Ifges(7,5) + Ifges(7,6) * t655 + t235 * Ifges(6,3) + t234 * Ifges(7,3) + t666 + t669;
t675 = t143 / 0.2e1 - t643 / 0.2e1;
t276 = Ifges(4,4) * t282;
t670 = Ifges(4,2) * t282;
t667 = t163 * mrSges(5,2);
t665 = t277 * Ifges(5,5);
t663 = t277 * Ifges(5,3);
t662 = t408 * Ifges(4,5);
t661 = t408 * Ifges(4,6);
t528 = t285 * mrSges(4,3);
t660 = mrSges(4,1) * t408 - mrSges(5,1) * t235 - mrSges(5,2) * t236 - t528;
t344 = t390 * t523 + t393 * t513;
t611 = qJD(4) * t344 + t390 * t683 + t393 * t465;
t417 = t384 * t547 + t454 * t545;
t396 = t417 * t385 + t545 * t477;
t448 = t383 * mrSges(6,1) + t387 * mrSges(6,2);
t496 = pkin(5) * t383 + pkin(10);
t659 = -m(7) * t496 - t448;
t450 = -mrSges(5,1) * t393 + t390 * mrSges(5,2);
t599 = -t390 * t607 + t393 * t608 - t450;
t649 = m(6) + m(7);
t658 = mrSges(4,1) + t599 + pkin(3) * (m(5) + t649);
t657 = -t379 * mrSges(7,1) - t380 * mrSges(7,2);
t504 = qJD(1) * qJD(3);
t457 = t524 * t504;
t226 = t386 * (qJDD(1) * t421 - t420 * t504) - (-qJDD(1) * t460 + t391 * t457) * t385;
t246 = t335 * t391 + t546 * t682;
t221 = qJDD(4) - t226;
t348 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t386;
t310 = qJDD(1) * t373 + t522 * t348;
t255 = qJDD(1) * t411 + t310;
t309 = qJDD(1) * t374 - t384 * t348;
t258 = (-t477 * t544 + t491) * qJDD(1) + t309;
t299 = qJDD(1) * t317 + qJDD(2);
t106 = qJD(3) * t438 + t546 * t255 + t258 * t475 - t271 * t508 + t299 * t513 + t306 * t463;
t323 = -qJDD(1) * t419 + qJDD(3);
t101 = pkin(10) * t323 + t106;
t214 = -t258 * t385 + t523 * t299;
t470 = qJDD(1) * t524;
t225 = (t391 * t470 + t457 * t546) * t385 + (qJDD(1) * t420 + t421 * t504) * t386;
t122 = -pkin(3) * t226 - pkin(10) * t225 + t214;
t505 = qJD(4) * t393;
t31 = t393 * t101 + t390 * t122 + t161 * t505 - t164 * t506;
t27 = qJ(5) * t221 + qJD(5) * t277 + t31;
t107 = -qJD(3) * t462 - t391 * t255 + t258 * t459 - t271 * t483 + t299 * t494 - t306 * t487;
t102 = -t323 * pkin(3) - t107;
t615 = qJD(4) * t408 + t225;
t147 = -t285 * t506 + t390 * t323 + t393 * t615;
t148 = t285 * t505 - t393 * t323 + t390 * t615;
t42 = t148 * pkin(4) - t147 * qJ(5) - t236 * qJD(5) + t102;
t10 = -t27 * t383 + t387 * t42;
t112 = t147 * t387 + t221 * t383;
t5 = pkin(5) * t148 - pkin(11) * t112 + t10;
t11 = t387 * t27 + t383 * t42;
t111 = -t147 * t383 + t221 * t387;
t6 = pkin(11) * t111 + t11;
t1 = qJD(6) * t12 + t389 * t5 + t392 * t6;
t2 = -qJD(6) * t13 - t389 * t6 + t392 * t5;
t654 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t564 = t148 / 0.2e1;
t573 = t112 / 0.2e1;
t653 = Ifges(6,1) * t573 + Ifges(6,5) * t564;
t560 = t221 / 0.2e1;
t565 = -t148 / 0.2e1;
t566 = t147 / 0.2e1;
t577 = Ifges(5,1) * t566 + Ifges(5,4) * t565 + Ifges(5,5) * t560;
t558 = t234 / 0.2e1;
t569 = t127 / 0.2e1;
t571 = t655 / 0.2e1;
t652 = Ifges(7,5) * t569 + Ifges(7,6) * t571 + Ifges(7,3) * t558 - t675;
t557 = -t235 / 0.2e1;
t559 = -t234 / 0.2e1;
t562 = -t193 / 0.2e1;
t563 = -t471 / 0.2e1;
t570 = -t127 / 0.2e1;
t572 = -t655 / 0.2e1;
t651 = Ifges(6,5) * t562 + Ifges(7,5) * t570 + Ifges(6,6) * t563 + Ifges(7,6) * t572 + Ifges(6,3) * t557 + Ifges(7,3) * t559 - t677;
t37 = qJD(6) * t655 + t111 * t389 + t112 * t392;
t585 = t37 / 0.2e1;
t38 = -qJD(6) * t127 + t111 * t392 - t112 * t389;
t584 = t38 / 0.2e1;
t650 = -m(6) - m(5);
t574 = t111 / 0.2e1;
t141 = qJDD(6) + t148;
t568 = t141 / 0.2e1;
t7 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t141;
t648 = t112 * Ifges(6,5) + t111 * Ifges(6,6) + t148 * Ifges(6,3) + t7;
t89 = t161 * t393 - t390 * t164;
t645 = t89 * mrSges(5,1);
t644 = t90 * mrSges(5,2);
t642 = Ifges(4,5) * t225;
t641 = Ifges(4,6) * t226;
t640 = Ifges(4,3) * t323;
t356 = -pkin(4) * t393 - t390 * qJ(5) - pkin(3);
t347 = t387 * t356;
t287 = -pkin(11) * t512 + t347 + (-pkin(10) * t383 - pkin(5)) * t393;
t322 = pkin(10) * t511 + t383 * t356;
t517 = t383 * t390;
t305 = -pkin(11) * t517 + t322;
t229 = t287 * t392 - t305 * t389;
t635 = qJD(6) * t229 + t389 * t679 + t678 * t392;
t230 = t287 * t389 + t305 * t392;
t634 = -qJD(6) * t230 - t678 * t389 + t392 * t679;
t357 = t541 * t383;
t358 = t541 * t387;
t307 = -t357 * t392 - t358 * t389;
t440 = t383 * t389 - t387 * t392;
t520 = t235 * t387;
t177 = pkin(4) * t236 + qJ(5) * t235;
t72 = t387 * t177 - t383 * t89;
t52 = pkin(5) * t236 + pkin(11) * t520 + t72;
t521 = t235 * t383;
t73 = t383 * t177 + t387 * t89;
t63 = pkin(11) * t521 + t73;
t633 = -qJD(5) * t440 + qJD(6) * t307 - t389 * t52 - t392 * t63;
t308 = -t357 * t389 + t358 * t392;
t351 = t383 * t392 + t387 * t389;
t632 = -qJD(5) * t351 - qJD(6) * t308 + t389 * t63 - t392 * t52;
t630 = -t387 * t502 + t685;
t131 = -t390 * t188 + t227 * t393;
t119 = -t285 * pkin(4) - t131;
t629 = t215 * pkin(5) + t496 * t505 - t119;
t113 = mrSges(5,1) * t221 - mrSges(5,3) * t147;
t62 = -t111 * mrSges(6,1) + t112 * mrSges(6,2);
t628 = t62 - t113;
t145 = t215 * t392 - t216 * t389;
t339 = t440 * qJD(6);
t269 = t339 * t390 - t351 * t505;
t622 = t145 - t269;
t146 = t215 * t389 + t216 * t392;
t340 = t351 * qJD(6);
t268 = -t340 * t390 - t440 * t505;
t621 = t146 - t268;
t300 = -t383 * t344 - t387 * t494;
t301 = t387 * t344 - t383 * t494;
t232 = t300 * t389 + t301 * t392;
t620 = -qJD(6) * t232 - t389 * t613 + t392 * t614;
t231 = t300 * t392 - t301 * t389;
t619 = qJD(6) * t231 + t389 * t614 + t392 * t613;
t165 = t351 * t235;
t618 = t165 + t340;
t166 = t440 * t235;
t617 = t166 + t339;
t128 = -mrSges(6,1) * t471 + mrSges(6,2) * t193;
t539 = mrSges(5,3) * t236;
t195 = mrSges(5,1) * t277 - t539;
t616 = t195 - t128;
t338 = qJ(2) * t476 + t373;
t286 = t411 + t338;
t296 = t374 + t414;
t196 = -t391 * t286 + t296 * t459 + t317 * t494;
t492 = t386 * t545;
t610 = -t385 * t492 + t417 * t523;
t609 = -m(6) * pkin(10) + t657 + t659;
t32 = -t390 * t101 + t122 * t393 - t161 * t506 - t164 * t505;
t30 = -t221 * pkin(4) + qJDD(5) - t32;
t85 = -t277 * pkin(4) + qJD(5) - t89;
t605 = t30 * t390 + t85 * t505;
t604 = t506 - t519;
t603 = t31 * t393 - t32 * t390;
t602 = -t10 * t383 + t11 * t387;
t601 = mrSges(5,1) + t608;
t590 = mrSges(5,2) + t607;
t598 = -g(1) * t492 + g(2) * t493 - g(3) * t524;
t595 = -m(5) * t163 + t660;
t594 = t32 * mrSges(5,1) - t31 * mrSges(5,2);
t593 = -t107 * mrSges(4,1) + t106 * mrSges(4,2);
t592 = -m(5) * pkin(10) + mrSges(4,2) + t609;
t589 = mrSges(4,2) - mrSges(5,3) + t659;
t588 = -mrSges(5,1) + t684;
t587 = Ifges(7,4) * t585 + Ifges(7,2) * t584 + Ifges(7,6) * t568;
t586 = Ifges(7,1) * t585 + Ifges(7,4) * t584 + Ifges(7,5) * t568;
t45 = Ifges(6,4) * t112 + Ifges(6,2) * t111 + Ifges(6,6) * t148;
t583 = t45 / 0.2e1;
t582 = Ifges(6,4) * t574 + t653;
t533 = Ifges(7,4) * t127;
t60 = Ifges(7,2) * t655 + Ifges(7,6) * t234 + t533;
t581 = -t60 / 0.2e1;
t580 = t60 / 0.2e1;
t123 = Ifges(7,4) * t655;
t61 = Ifges(7,1) * t127 + Ifges(7,5) * t234 + t123;
t579 = -t61 / 0.2e1;
t578 = t61 / 0.2e1;
t99 = t193 * Ifges(6,4) + Ifges(6,2) * t471 + Ifges(6,6) * t235;
t576 = -t99 / 0.2e1;
t100 = t193 * Ifges(6,1) + Ifges(6,4) * t471 + Ifges(6,5) * t235;
t575 = t100 / 0.2e1;
t556 = t235 / 0.2e1;
t555 = -t236 / 0.2e1;
t554 = t236 / 0.2e1;
t552 = -t277 / 0.2e1;
t549 = t285 / 0.2e1;
t284 = t295 * qJD(3);
t173 = qJD(2) * t413 + qJD(3) * t196;
t237 = -t296 * t385 + t523 * t317;
t292 = t399 * pkin(3);
t479 = t295 * pkin(10) + t292;
t180 = t237 - t479;
t275 = t546 * t286;
t472 = t523 * t296;
t197 = t317 * t513 + t391 * t472 + t275;
t187 = -pkin(10) * t419 + t197;
t283 = t399 * qJD(3);
t488 = qJD(2) * t515;
t464 = t385 * t488;
t213 = pkin(3) * t284 - pkin(10) * t283 + t464;
t66 = t393 * t173 + t180 * t505 - t187 * t506 + t390 * t213;
t54 = qJ(5) * t284 - qJD(5) * t399 + t66;
t174 = qJD(2) * t412 + (t275 + (t317 * t385 + t472) * t391) * qJD(3);
t245 = t295 * t393 - t390 * t419;
t202 = qJD(4) * t245 + t283 * t390;
t244 = t295 * t390 + t393 * t419;
t203 = -qJD(4) * t244 + t283 * t393;
t80 = t202 * pkin(4) - t203 * qJ(5) - t245 * qJD(5) + t174;
t25 = t383 * t80 + t387 * t54;
t540 = mrSges(5,3) * t235;
t538 = Ifges(4,4) * t285;
t537 = Ifges(5,4) * t390;
t536 = Ifges(5,4) * t393;
t535 = Ifges(6,4) * t383;
t534 = Ifges(6,4) * t387;
t529 = t282 * mrSges(4,3);
t186 = pkin(3) * t419 - t196;
t121 = t244 * pkin(4) - t245 * qJ(5) + t186;
t109 = t390 * t180 + t393 * t187;
t94 = -qJ(5) * t399 + t109;
t56 = t383 * t121 + t387 * t94;
t518 = t282 * t393;
t510 = t547 * pkin(1) + qJ(2) * t492;
t507 = qJD(4) * t235;
t503 = qJDD(1) * t386;
t501 = pkin(10) * t505;
t499 = t383 * t576;
t498 = Ifges(5,5) * t147 - Ifges(5,6) * t148 + Ifges(5,3) * t221;
t497 = t640 + t641 + t642;
t14 = -t38 * mrSges(7,1) + t37 * mrSges(7,2);
t482 = t384 * t503;
t480 = t505 / 0.2e1;
t24 = -t383 * t54 + t387 * t80;
t55 = t387 * t121 - t383 * t94;
t108 = t180 * t393 - t390 * t187;
t461 = -pkin(1) * t545 + qJ(2) * t493;
t456 = qJDD(1) * t476;
t452 = -mrSges(4,1) * t399 + mrSges(4,2) * t295;
t451 = mrSges(5,1) * t244 + mrSges(5,2) * t245;
t447 = Ifges(5,1) * t393 - t537;
t446 = Ifges(6,1) * t387 - t535;
t445 = -Ifges(5,2) * t390 + t536;
t444 = -Ifges(6,2) * t383 + t534;
t443 = Ifges(5,5) * t393 - Ifges(5,6) * t390;
t442 = Ifges(6,5) * t387 - Ifges(6,6) * t383;
t205 = t245 * t387 - t383 * t399;
t43 = pkin(5) * t244 - pkin(11) * t205 + t55;
t204 = -t245 * t383 - t387 * t399;
t48 = pkin(11) * t204 + t56;
t15 = -t389 * t48 + t392 * t43;
t16 = t389 * t43 + t392 * t48;
t134 = t204 * t392 - t205 * t389;
t135 = t204 * t389 + t205 * t392;
t67 = -t390 * t173 - t180 * t506 - t187 * t505 + t213 * t393;
t95 = pkin(4) * t399 - t108;
t435 = t163 * (mrSges(5,1) * t390 + mrSges(5,2) * t393);
t434 = -(-qJ(2) * t489 + t367) * t384 + t328 * t522;
t336 = -t478 * t545 + t522 * t547;
t251 = t336 * t546 - t391 * t610;
t210 = t251 * t390 - t393 * t396;
t431 = -g(1) * t210 + g(2) * t680 - g(3) * t244;
t429 = -mrSges(3,1) * t456 + mrSges(3,2) * t482;
t427 = mrSges(3,1) * t524 - mrSges(3,3) * t515;
t58 = -t284 * pkin(4) - t67;
t425 = -mrSges(3,2) * t524 + mrSges(3,3) * t476;
t401 = -t335 * pkin(2) - pkin(9) * t403 + t461;
t400 = t336 * pkin(2) + pkin(9) * t396 + t510;
t398 = t249 * pkin(3) + t401;
t397 = t251 * pkin(3) + t400;
t369 = -pkin(1) * t503 + qJDD(2);
t353 = t496 * t390;
t342 = t425 * qJD(1);
t341 = t427 * qJD(1);
t337 = -qJ(2) * t515 + t374;
t330 = t440 * t390;
t329 = t351 * t390;
t321 = -pkin(10) * t516 + t347;
t250 = t336 * t391 + t546 * t610;
t238 = -mrSges(4,2) * t408 + t529;
t233 = Ifges(5,4) * t235;
t224 = -mrSges(4,1) * t282 + mrSges(4,2) * t285;
t211 = t251 * t393 + t390 * t396;
t201 = Ifges(4,1) * t285 + t276 + t662;
t200 = t538 + t661 + t670;
t199 = -mrSges(4,2) * t323 + mrSges(4,3) * t226;
t198 = mrSges(4,1) * t323 - mrSges(4,3) * t225;
t194 = -mrSges(5,2) * t277 - t540;
t168 = t203 * t387 + t284 * t383;
t167 = -t203 * t383 + t284 * t387;
t157 = -mrSges(4,1) * t226 + mrSges(4,2) * t225;
t156 = t211 * t380 + t250 * t379;
t155 = -t211 * t379 + t250 * t380;
t144 = t236 * Ifges(5,1) - t233 + t665;
t142 = t236 * Ifges(5,5) - t235 * Ifges(5,6) + t663;
t137 = mrSges(6,1) * t235 - mrSges(6,3) * t193;
t136 = -mrSges(6,2) * t235 + mrSges(6,3) * t471;
t114 = -mrSges(5,2) * t221 - mrSges(5,3) * t148;
t93 = mrSges(7,1) * t234 - mrSges(7,3) * t127;
t92 = -mrSges(7,2) * t234 + mrSges(7,3) * t655;
t83 = mrSges(5,1) * t148 + mrSges(5,2) * t147;
t82 = -pkin(5) * t521 + t90;
t81 = -t204 * pkin(5) + t95;
t76 = t147 * Ifges(5,4) - t148 * Ifges(5,2) + t221 * Ifges(5,6);
t75 = mrSges(6,1) * t148 - mrSges(6,3) * t112;
t74 = -mrSges(6,2) * t148 + mrSges(6,3) * t111;
t71 = -pkin(5) * t471 + t85;
t68 = -mrSges(7,1) * t655 + mrSges(7,2) * t127;
t65 = -qJD(6) * t135 + t167 * t392 - t168 * t389;
t64 = qJD(6) * t134 + t167 * t389 + t168 * t392;
t47 = -t167 * pkin(5) + t58;
t29 = -mrSges(7,2) * t141 + mrSges(7,3) * t38;
t28 = mrSges(7,1) * t141 - mrSges(7,3) * t37;
t21 = pkin(11) * t167 + t25;
t18 = -t111 * pkin(5) + t30;
t17 = pkin(5) * t202 - pkin(11) * t168 + t24;
t4 = -qJD(6) * t16 + t17 * t392 - t21 * t389;
t3 = qJD(6) * t15 + t17 * t389 + t21 * t392;
t8 = [(-Ifges(4,4) * t549 + Ifges(5,5) * t554 + Ifges(5,6) * t557 - t189 * mrSges(4,3) + t228 * mrSges(4,1) - t200 / 0.2e1 + t663 / 0.2e1 + t142 / 0.2e1 - t670 / 0.2e1 - t661 / 0.2e1 + t645 - t644) * t284 + t471 * (Ifges(6,4) * t168 + Ifges(6,2) * t167) / 0.2e1 + (Ifges(6,4) * t205 + Ifges(6,2) * t204) * t574 + (Ifges(7,5) * t64 + Ifges(7,6) * t65) * t558 + (Ifges(7,5) * t135 + Ifges(7,6) * t134) * t568 + m(3) * (t309 * t337 + t310 * t338) + ((Ifges(3,4) * t384 + Ifges(3,2) * t522) * t456 + (Ifges(3,5) * t384 + Ifges(3,6) * t522) * t470 + (Ifges(3,1) * t384 + Ifges(3,4) * t522) * t482 - pkin(1) * t429 + t369 * (-mrSges(3,1) * t522 + mrSges(3,2) * t384) + m(3) * (-pkin(1) * t369 + qJD(2) * t434)) * t386 + (Ifges(5,1) * t554 - t89 * mrSges(5,3) + Ifges(5,4) * t557 + t665 / 0.2e1 + t667 + t144 / 0.2e1) * t203 + (Ifges(6,5) * t168 + Ifges(6,6) * t167) * t556 + (Ifges(6,5) * t205 + Ifges(6,6) * t204) * t564 + (-t31 * mrSges(5,3) - Ifges(5,4) * t566 + Ifges(6,3) * t564 + Ifges(7,3) * t568 - Ifges(5,6) * t560 - t76 / 0.2e1 - t11 * mrSges(6,2) + t10 * mrSges(6,1) - Ifges(5,2) * t565 + t648 / 0.2e1 + Ifges(6,5) * t573 + Ifges(6,6) * t574 + Ifges(7,6) * t584 + Ifges(7,5) * t585 + t654) * t244 + (-m(4) * t188 - t595) * t174 + (Ifges(7,1) * t64 + Ifges(7,4) * t65) * t569 + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t585 + (-m(3) * t461 - m(4) * t401 - m(7) * t398 + mrSges(2,1) * t545 + t335 * mrSges(3,1) - t249 * mrSges(4,1) + mrSges(2,2) * t547 - mrSges(3,2) * t426 - mrSges(3,3) * t493 + mrSges(4,3) * t403 + t650 * (-pkin(10) * t246 + t398) + (t588 + t656) * t681 - (t589 + t657) * t246 + t590 * t680) * g(1) + (Ifges(4,1) * t549 - t188 * mrSges(4,3) + t228 * mrSges(4,2) + t201 / 0.2e1 + t276 / 0.2e1 + t662 / 0.2e1) * t283 + (-mrSges(5,3) * t32 + 0.2e1 * t577) * t245 + (-mrSges(4,3) * t107 + Ifges(4,1) * t225 + Ifges(4,4) * t226 + Ifges(4,5) * t323) * t295 + (Ifges(3,5) * t482 + Ifges(3,6) * t456 + Ifges(3,3) * t470) * t524 + (t337 * t427 + t338 * t425 + Ifges(2,3)) * qJDD(1) + m(6) * (t10 * t55 + t11 * t56 + t24 * t49 + t25 * t50 + t30 * t95 + t58 * t85) + m(7) * (t1 * t16 + t12 * t4 + t13 * t3 + t15 * t2 + t18 * t81 + t47 * t71) + (-Ifges(5,4) * t554 - t90 * mrSges(5,3) + t669 / 0.2e1 + Ifges(6,3) * t556 - Ifges(5,2) * t557 - t664 / 0.2e1 + t666 / 0.2e1 + t652 + t676 + t677) * t202 - t341 * t488 + t193 * (Ifges(6,1) * t168 + Ifges(6,4) * t167) / 0.2e1 + (Ifges(6,1) * t205 + Ifges(6,4) * t204) * t573 + m(4) * (t106 * t197 + t107 * t196 + t173 * t189 + t214 * t237 + t228 * t464) + (t1 * t134 - t12 * t64 + t13 * t65 - t135 * t2) * mrSges(7,3) + t309 * t427 + t102 * t451 + t214 * t452 + t224 * t464 + (-t10 * t205 + t11 * t204 + t167 * t50 - t168 * t49) * mrSges(6,3) - (Ifges(5,3) * t560 + Ifges(5,6) * t565 + Ifges(5,5) * t566 + t498 / 0.2e1 - Ifges(4,4) * t225 - Ifges(4,2) * t226 - Ifges(4,6) * t323 - t106 * mrSges(4,3) + t594) * t399 + m(5) * (t102 * t186 + t108 * t32 + t109 * t31 + t66 * t90 + t67 * t89) + (Ifges(7,4) * t64 + Ifges(7,2) * t65) * t571 + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t584 + t310 * t425 + t237 * t157 + t173 * t238 + t30 * (-mrSges(6,1) * t204 + mrSges(6,2) * t205) + t66 * t194 + t67 * t195 + t196 * t198 + t197 * t199 + t186 * t83 + t85 * (-mrSges(6,1) * t167 + mrSges(6,2) * t168) + t167 * t99 / 0.2e1 + t25 * t136 + t24 * t137 + t18 * (-mrSges(7,1) * t134 + mrSges(7,2) * t135) + t58 * t128 + t109 * t114 + t108 * t113 + t95 * t62 + t3 * t92 + t4 * t93 + t81 * t14 + t56 * t74 + t55 * t75 + t71 * (-mrSges(7,1) * t65 + mrSges(7,2) * t64) + t47 * t68 + (-m(3) * t510 - m(4) * t400 - m(7) * t397 - mrSges(2,1) * t547 - t336 * mrSges(3,1) - t251 * mrSges(4,1) - t156 * mrSges(7,1) + mrSges(2,2) * t545 + mrSges(3,2) * t417 - t155 * mrSges(7,2) - mrSges(3,3) * t492 - mrSges(4,3) * t396 + t650 * (t250 * pkin(10) + t397) + t588 * t211 + t589 * t250 + t590 * t210) * g(2) + t168 * t575 + qJD(2) * t342 * t476 + t64 * t578 + t65 * t580 + t205 * t582 + t204 * t583 + t135 * t586 + t134 * t587 + t15 * t28 + t16 * t29 + (-t497 / 0.2e1 - t640 / 0.2e1 - t642 / 0.2e1 - t641 / 0.2e1 + t593) * t419; t611 * (t68 - t616) + (t198 - t83) * t494 + t683 * t238 + t523 * t157 + t429 + t199 * t513 - t342 * t458 - t224 * t465 + t344 * t114 + t300 * t75 + t301 * t74 + t231 * t28 + t232 * t29 + (t214 * t523 + (t546 * t107 + t106 * t391 + (-t188 * t391 + t189 * t546) * qJD(3)) * t385 + t188 * t312 - t189 * t313 - t228 * t465 + t598) * m(4) + (-t434 * t509 + t369 + t598) * m(3) + t341 * t489 + t612 * t194 + (t31 * t344 - t32 * t343 + (-t102 * t546 + t163 * t508) * t385 - t163 * t312 + t612 * t90 - t611 * t89 + t598) * m(5) + t613 * t136 + t614 * t137 + (t10 * t300 + t11 * t301 + t30 * t343 + t49 * t614 + t50 * t613 + t611 * t85 + t598) * m(6) + t619 * t92 + t620 * t93 + (t1 * t232 + t12 * t620 + t13 * t619 + t18 * t343 + t2 * t231 + t611 * t71 + t598) * m(7) + t660 * t686 + (t14 + t628) * t343; -(Ifges(4,1) * t282 + t142 - t538) * t285 / 0.2e1 - (-Ifges(4,2) * t285 + t201 + t276) * t282 / 0.2e1 + (t49 * (mrSges(6,1) * t390 - mrSges(6,3) * t511) + t50 * (-mrSges(6,2) * t390 - mrSges(6,3) * t516) + t435) * qJD(4) + t200 * t549 + (Ifges(5,3) * t285 + t282 * t443) * t552 + (-t102 * pkin(3) - t131 * t89 - t132 * t90) * m(5) + (t10 * t321 + t11 * t322 - t119 * t85 + t49 * t631 + t50 * t630) * m(6) + (t393 * t114 + t628 * t390 + ((-t390 * t90 - t393 * t89) * qJD(4) + t603) * m(5) + t605 * m(6)) * pkin(10) + (t85 * mrSges(6,1) - mrSges(6,3) * t50 + Ifges(6,4) * t562 + Ifges(6,2) * t563 + Ifges(6,6) * t557 + t576) * t215 + (-t502 - t132) * t194 + (t193 * (Ifges(6,5) * t390 + t393 * t446) + t471 * (Ifges(6,6) * t390 + t393 * t444) + t277 * t443 + t236 * t447) * qJD(4) / 0.2e1 + t1 * (mrSges(7,2) * t393 - t329 * mrSges(7,3)) + (-Ifges(7,5) * t330 - Ifges(7,6) * t329 - Ifges(7,3) * t393) * t568 + t18 * (mrSges(7,1) * t329 - mrSges(7,2) * t330) + (-Ifges(7,4) * t330 - Ifges(7,2) * t329 - Ifges(7,6) * t393) * t584 + (-Ifges(7,1) * t330 - Ifges(7,4) * t329 - Ifges(7,5) * t393) * t585 + t2 * (-mrSges(7,1) * t393 + t330 * mrSges(7,3)) + (t528 + t595) * t189 + (t651 + t675) * t519 + t652 * t506 + (t246 * t658 - t249 * t592) * g(2) + (-g(1) * t251 + g(2) * t249 - g(3) * t295 - t604 * t90 + (-t505 + t518) * t89 + t603) * mrSges(5,3) + (t529 - t238) * t188 + t11 * (mrSges(6,2) * t393 - mrSges(6,3) * t517) - t45 * t517 / 0.2e1 + (-t501 - t131) * t195 + (t250 * t658 + t251 * t592) * g(1) + (Ifges(5,5) * t285 + t282 * t447) * t555 + (Ifges(5,6) * t285 + t282 * t445) * t556 + (Ifges(5,5) * t390 + Ifges(5,6) * t393) * t560 + (-Ifges(6,3) * t393 + t390 * t442) * t564 + (Ifges(5,2) * t393 + t537) * t565 + (Ifges(5,1) * t390 + t536) * t566 + (-t119 + t501) * t128 + t10 * (-mrSges(6,1) * t393 - mrSges(6,3) * t512) + (Ifges(6,5) * t557 + Ifges(6,1) * t562 + Ifges(6,4) * t563 + t49 * mrSges(6,3) - t85 * mrSges(6,2) - t100 / 0.2e1) * t216 + t499 * t505 - t282 * t435 + (Ifges(6,3) * t390 + t393 * t442) * t507 / 0.2e1 - t445 * t507 / 0.2e1 + t102 * t450 + (Ifges(7,5) * t268 + Ifges(7,6) * t269) * t558 + (Ifges(7,5) * t146 + Ifges(7,6) * t145) * t559 + (-m(5) * t479 - t292 * t649 + t295 * t609 - t399 * t599 + t452) * g(3) + (-t518 / 0.2e1 + t480) * t144 + t497 - t593 + (Ifges(7,4) * t268 + Ifges(7,2) * t269) * t571 + (Ifges(7,4) * t146 + Ifges(7,2) * t145) * t572 + (Ifges(7,1) * t146 + Ifges(7,4) * t145) * t570 + (Ifges(7,1) * t268 + Ifges(7,4) * t269) * t569 + t353 * t14 + t321 * t75 + t322 * t74 - t228 * (mrSges(4,1) * t285 + mrSges(4,2) * t282) + t229 * t28 + t230 * t29 - t408 * (Ifges(4,5) * t282 - Ifges(4,6) * t285) / 0.2e1 - pkin(3) * t83 + t285 * t644 - t648 * t393 / 0.2e1 - t285 * t645 + (-Ifges(6,5) * t393 + t390 * t446) * t573 + (-Ifges(6,6) * t393 + t390 * t444) * t574 + t390 * t577 + t387 * t100 * t480 + t605 * t448 + (mrSges(7,1) * t604 + mrSges(7,3) * t621) * t12 + (mrSges(7,1) * t622 - mrSges(7,2) * t621) * t71 + (-mrSges(7,2) * t604 - mrSges(7,3) * t622) * t13 + t393 * t76 / 0.2e1 + t268 * t578 + t146 * t579 + t269 * t580 + t145 * t581 + t512 * t582 - t330 * t586 - t329 * t587 + t629 * t68 + t630 * t136 + t631 * t137 + t634 * t93 + t635 * t92 + (t1 * t230 + t12 * t634 + t13 * t635 + t18 * t353 + t2 * t229 + t629 * t71) * m(7); (Ifges(6,2) * t574 + Ifges(6,6) * t564 + qJ(5) * t74 + qJD(5) * t136 + t583) * t387 + (Ifges(7,4) * t166 + Ifges(7,2) * t165) * t572 + t143 * t554 + (Ifges(7,5) * t351 - Ifges(7,6) * t440) * t568 + t18 * (mrSges(7,1) * t440 + mrSges(7,2) * t351) + (Ifges(7,4) * t351 - Ifges(7,2) * t440) * t584 + (Ifges(7,1) * t351 - Ifges(7,4) * t440) * t585 + (-t1 * t440 + t12 * t617 - t13 * t618 - t2 * t351) * mrSges(7,3) - t440 * t587 + t535 * t574 + (-Ifges(5,5) * t552 - t442 * t557 - t444 * t563 - t446 * t562 + t448 * t85 + t499 + t667) * t235 + (-Ifges(5,2) * t556 - Ifges(5,6) * t552 + t651 - t676) * t236 + (Ifges(7,5) * t166 + Ifges(7,6) * t165) * t559 + (Ifges(7,1) * t166 + Ifges(7,4) * t165) * t570 + (-qJ(5) * t75 - qJD(5) * t137 + t582 + t653) * t383 + t534 * t573 + (-t590 * t681 - t601 * t680) * g(2) + (-t540 - t194) * t89 + (-Ifges(7,5) * t339 - Ifges(7,6) * t340) * t558 + (-Ifges(7,1) * t339 - Ifges(7,4) * t340) * t569 + (-Ifges(7,4) * t339 - Ifges(7,2) * t340) * t571 + (-t233 + t144) * t556 + t594 - t378 * t14 + t30 * t449 + t498 + t307 * t28 + t308 * t29 - t73 * t136 - t72 * t137 - t82 * t68 - pkin(4) * t62 + (-Ifges(5,1) * t235 - t530 + t643) * t555 + t520 * t575 + (t210 * t601 + t211 * t590) * g(1) + (-t49 * t520 - t50 * t521 + t602) * mrSges(6,3) + (-pkin(4) * t30 + (-t383 * t49 + t387 * t50) * qJD(5) + t602 * qJ(5) - t49 * t72 - t50 * t73) * m(6) + (t244 * t608 + t245 * t607 + t451) * g(3) + (-m(6) * t85 + t539 + t616) * t90 + (mrSges(7,1) * t618 - mrSges(7,2) * t617) * t71 - t339 * t578 + t166 * t579 - t340 * t580 + t165 * t581 + t351 * t586 + t632 * t93 + t633 * t92 + (t1 * t308 + t12 * t632 + t13 * t633 - t18 * t378 + t2 * t307 - t71 * t82) * m(7); t127 * t93 - t655 * t92 - t471 * t136 + t193 * t137 + t14 + t62 + (t12 * t127 - t13 * t655 + t18 + t431) * m(7) + (t193 * t49 - t471 * t50 + t30 + t431) * m(6); -t71 * (mrSges(7,1) * t127 + mrSges(7,2) * t655) + (Ifges(7,1) * t655 - t533) * t570 + t60 * t569 + (Ifges(7,5) * t655 - Ifges(7,6) * t127) * t559 - t12 * t92 + t13 * t93 - g(1) * (mrSges(7,1) * t155 - mrSges(7,2) * t156) - g(2) * ((t246 * t380 + t379 * t681) * mrSges(7,1) + (-t246 * t379 + t380 * t681) * mrSges(7,2)) - g(3) * ((-t245 * t379 - t380 * t399) * mrSges(7,1) + (-t245 * t380 + t379 * t399) * mrSges(7,2)) + (t12 * t655 + t127 * t13) * mrSges(7,3) + t7 + (-Ifges(7,2) * t127 + t123 + t61) * t572 + t654;];
tau  = t8;
