% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:57
% EndTime: 2019-03-09 01:05:26
% DurationCPUTime: 51.50s
% Computational Cost: add. (22860->1035), mult. (56246->1465), div. (0->0), fcn. (47571->18), ass. (0->465)
t384 = sin(pkin(6));
t385 = cos(pkin(7));
t395 = cos(qJ(3));
t575 = cos(qJ(2));
t491 = t575 * t395;
t390 = sin(qJ(3));
t391 = sin(qJ(2));
t524 = t390 * t391;
t415 = -t385 * t524 + t491;
t288 = t415 * t384;
t266 = qJD(1) * t288;
t383 = sin(pkin(7));
t533 = t383 * t390;
t372 = pkin(9) * t533;
t528 = t385 * t395;
t330 = pkin(2) * t528 - t372;
t313 = t330 * qJD(3);
t639 = t313 - t266;
t529 = t385 * t390;
t531 = t383 * t395;
t331 = pkin(2) * t529 + pkin(9) * t531;
t300 = pkin(10) * t385 + t331;
t450 = -pkin(3) * t395 - pkin(10) * t390;
t301 = (-pkin(2) + t450) * t383;
t428 = t383 * (pkin(3) * t390 - pkin(10) * t395);
t312 = qJD(3) * t428;
t389 = sin(qJ(4));
t394 = cos(qJ(4));
t530 = t384 * t391;
t489 = qJD(1) * t530;
t464 = t383 * t489;
t510 = qJD(4) * t394;
t511 = qJD(4) * t389;
t644 = -t300 * t511 + t301 * t510 + t639 * t394 + (t312 - t464) * t389;
t492 = t575 * t390;
t523 = t391 * t395;
t414 = t385 * t523 + t492;
t287 = t414 * t384;
t694 = -qJD(1) * t287 + t331 * qJD(3);
t513 = qJD(3) * t390;
t484 = t383 * t513;
t693 = pkin(11) * t484 + t644;
t326 = -t394 * t385 + t389 * t533;
t512 = qJD(3) * t395;
t483 = t383 * t512;
t252 = -qJD(4) * t326 + t394 * t483;
t532 = t383 * t394;
t327 = t385 * t389 + t390 * t532;
t253 = qJD(4) * t327 + t389 * t483;
t692 = pkin(4) * t253 - pkin(11) * t252 + t694;
t493 = t384 * t575;
t460 = qJD(1) * t493;
t352 = qJD(2) * pkin(2) + t460;
t449 = pkin(4) * t389 - pkin(11) * t394;
t514 = qJD(2) * t395;
t386 = cos(pkin(6));
t516 = qJD(1) * t386;
t515 = qJD(2) * t383;
t333 = pkin(9) * t515 + t489;
t536 = t333 * t395;
t691 = -t352 * t529 - t536 - (t390 * t516 + t449 * t514) * t383 + t449 * qJD(4);
t388 = sin(qJ(5));
t393 = cos(qJ(5));
t520 = t394 * t395;
t270 = (-t388 * t520 + t390 * t393) * t515;
t507 = qJD(5) * t393;
t642 = t388 * t510 + t389 * t507 + t270;
t299 = t372 + (-pkin(2) * t395 - pkin(3)) * t385;
t206 = pkin(4) * t326 - pkin(11) * t327 + t299;
t220 = t394 * t300 + t389 * t301;
t208 = -pkin(11) * t531 + t220;
t103 = t388 * t206 + t393 * t208;
t655 = -qJD(5) * t103 - t388 * t693 + t692 * t393;
t509 = qJD(5) * t388;
t654 = t206 * t507 - t208 * t509 + t692 * t388 + t393 * t693;
t426 = t352 * t385 + t383 * t516;
t217 = -t390 * t333 + t395 * t426;
t311 = qJD(2) * t428;
t155 = t394 * t217 + t389 * t311;
t488 = t390 * t515;
t142 = pkin(11) * t488 + t155;
t356 = -pkin(4) * t394 - pkin(11) * t389 - pkin(3);
t651 = -t393 * t142 + t356 * t507 + (-t393 * t511 - t394 * t509) * pkin(10) + t691 * t388;
t502 = pkin(10) * t511;
t690 = t691 * t393 + (t142 + t502) * t388;
t371 = qJD(2) * t385 + qJD(3);
t285 = t371 * t394 - t389 * t488;
t461 = t394 * t488;
t286 = t371 * t389 + t461;
t549 = t286 * Ifges(5,4);
t487 = t383 * t514;
t357 = qJD(4) - t487;
t658 = t357 * Ifges(5,6);
t179 = t285 * Ifges(5,2) + t549 + t658;
t276 = qJD(5) - t285;
t264 = qJD(6) + t276;
t585 = t264 / 0.2e1;
t235 = -t286 * t388 + t357 * t393;
t236 = t286 * t393 + t357 * t388;
t387 = sin(qJ(6));
t392 = cos(qJ(6));
t135 = t235 * t387 + t236 * t392;
t597 = t135 / 0.2e1;
t469 = t392 * t235 - t236 * t387;
t599 = t469 / 0.2e1;
t218 = t390 * t426 + t536;
t197 = pkin(10) * t371 + t218;
t367 = t385 * t516;
t232 = t367 + (qJD(2) * t450 - t352) * t383;
t116 = t197 * t394 + t232 * t389;
t105 = pkin(11) * t357 + t116;
t196 = -pkin(3) * t371 - t217;
t118 = -pkin(4) * t285 - pkin(11) * t286 + t196;
t54 = -t105 * t388 + t393 * t118;
t55 = t105 * t393 + t118 * t388;
t623 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t689 = Ifges(7,5) * t597 + Ifges(7,6) * t599 + Ifges(7,3) * t585 - t179 / 0.2e1 - t623;
t670 = -m(6) - m(5);
t254 = -t327 * t388 - t393 * t531;
t147 = qJD(5) * t254 + t252 * t393 + t388 * t484;
t688 = pkin(5) * t253 - pkin(12) * t147 + t655;
t427 = -t327 * t393 + t388 * t531;
t148 = qJD(5) * t427 - t252 * t388 + t393 * t484;
t687 = -pkin(12) * t148 - t654;
t271 = (t388 * t390 + t393 * t520) * t515;
t521 = t393 * t394;
t377 = pkin(10) * t521;
t462 = t389 * t487;
t686 = -pkin(5) * t462 + pkin(12) * t271 + (pkin(5) * t389 - pkin(12) * t521) * qJD(4) + (-t377 + (pkin(12) * t389 - t356) * t388) * qJD(5) + t690;
t685 = -pkin(12) * t642 + t651;
t396 = -pkin(12) - pkin(11);
t494 = qJD(5) * t396;
t540 = t285 * t388;
t115 = -t389 * t197 + t232 * t394;
t216 = pkin(4) * t286 - pkin(11) * t285;
t79 = t393 * t115 + t388 * t216;
t684 = pkin(12) * t540 + t388 * t494 - t79;
t539 = t285 * t393;
t78 = -t115 * t388 + t393 * t216;
t683 = -pkin(5) * t286 + pkin(12) * t539 + t393 * t494 - t78;
t379 = pkin(5) * t393 + pkin(4);
t382 = qJ(5) + qJ(6);
t380 = sin(t382);
t381 = cos(t382);
t447 = -mrSges(6,1) * t393 + mrSges(6,2) * t388;
t682 = m(6) * pkin(4) + m(7) * t379 + mrSges(7,1) * t381 - mrSges(7,2) * t380 - t447;
t681 = m(6) * pkin(11) - m(7) * t396 + mrSges(6,3) + mrSges(7,3);
t680 = t509 - t540;
t506 = qJD(2) * qJD(3);
t319 = (-qJDD(2) * t395 + t390 * t506) * t383;
t304 = qJDD(4) + t319;
t486 = qJD(2) * t530;
t459 = qJD(1) * t486;
t321 = qJDD(1) * t493 - t459;
t295 = qJDD(2) * pkin(2) + t321;
t505 = qJDD(1) * t386;
t241 = -t295 * t383 + t385 * t505;
t320 = (qJDD(2) * t390 + t395 * t506) * t383;
t153 = pkin(3) * t319 - pkin(10) * t320 + t241;
t479 = t383 * t505;
t363 = qJD(2) * t460;
t322 = qJDD(1) * t530 + t363;
t675 = pkin(9) * qJDD(2) * t383 + qJD(3) * t426 + t322;
t100 = t295 * t529 - t333 * t513 + t390 * t479 + t395 * t675;
t370 = qJDD(2) * t385 + qJDD(3);
t95 = pkin(10) * t370 + t100;
t41 = t389 * t153 - t197 * t511 + t232 * t510 + t394 * t95;
t34 = pkin(11) * t304 + t41;
t192 = qJD(4) * t285 + t320 * t394 + t370 * t389;
t193 = -qJD(4) * t461 - t389 * t320 + t370 * t394 - t371 * t511;
t101 = t395 * (t295 * t385 + t479) - t333 * t512 - t675 * t390;
t96 = -pkin(3) * t370 - t101;
t53 = -pkin(4) * t193 - pkin(11) * t192 + t96;
t9 = -qJD(5) * t55 - t34 * t388 + t393 * t53;
t668 = t9 * mrSges(6,1);
t8 = -t105 * t509 + t118 * t507 + t393 * t34 + t388 * t53;
t669 = t8 * mrSges(6,2);
t679 = t668 - t669;
t48 = -pkin(12) * t236 + t54;
t47 = pkin(5) * t276 + t48;
t49 = pkin(12) * t235 + t55;
t547 = t387 * t49;
t16 = t392 * t47 - t547;
t545 = t392 * t49;
t17 = t387 * t47 + t545;
t551 = t116 * mrSges(5,3);
t678 = t551 - t196 * mrSges(5,1) - t16 * mrSges(7,1) + t17 * mrSges(7,2) + t658 / 0.2e1;
t188 = qJDD(5) - t193;
t174 = qJDD(6) + t188;
t596 = t174 / 0.2e1;
t97 = qJD(5) * t235 + t192 * t393 + t304 * t388;
t98 = -qJD(5) * t236 - t192 * t388 + t304 * t393;
t32 = -qJD(6) * t135 - t387 * t97 + t392 * t98;
t611 = t32 / 0.2e1;
t31 = qJD(6) * t469 + t387 * t98 + t392 * t97;
t612 = t31 / 0.2e1;
t613 = Ifges(7,1) * t612 + Ifges(7,4) * t611 + Ifges(7,5) * t596;
t614 = Ifges(7,4) * t612 + Ifges(7,2) * t611 + Ifges(7,6) * t596;
t656 = t236 * Ifges(6,5) + t135 * Ifges(7,5) + t235 * Ifges(6,6) + Ifges(7,6) * t469 + t276 * Ifges(6,3) + t264 * Ifges(7,3);
t154 = -t389 * t217 + t311 * t394;
t141 = -pkin(4) * t488 - t154;
t501 = pkin(10) * t510;
t677 = -t141 + t501;
t448 = -mrSges(5,1) * t394 + mrSges(5,2) * t389;
t632 = m(7) - t670;
t676 = pkin(3) * t632 + t389 * t681 + t394 * t682 + mrSges(4,1) - t448;
t446 = mrSges(6,1) * t388 + mrSges(6,2) * t393;
t495 = pkin(5) * t388 + pkin(10);
t619 = -m(7) * t495 - t380 * mrSges(7,1) - t381 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - t446;
t468 = mrSges(4,3) * t488;
t624 = -m(5) * t196 + mrSges(4,1) * t371 + mrSges(5,1) * t285 - mrSges(5,2) * t286 - t468;
t674 = -m(4) * t217 - t624;
t6 = pkin(5) * t188 - pkin(12) * t97 + t9;
t7 = pkin(12) * t98 + t8;
t2 = qJD(6) * t16 + t387 * t6 + t392 * t7;
t3 = -qJD(6) * t17 - t387 * t7 + t392 * t6;
t673 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t603 = t97 / 0.2e1;
t602 = t98 / 0.2e1;
t594 = t188 / 0.2e1;
t593 = t192 / 0.2e1;
t592 = t193 / 0.2e1;
t578 = t304 / 0.2e1;
t102 = t393 * t206 - t208 * t388;
t73 = pkin(5) * t326 + pkin(12) * t427 + t102;
t85 = pkin(12) * t254 + t103;
t37 = t387 * t73 + t392 * t85;
t667 = -qJD(6) * t37 + t387 * t687 + t392 * t688;
t36 = -t387 * t85 + t392 * t73;
t666 = qJD(6) * t36 + t387 * t688 - t392 * t687;
t337 = t393 * t356;
t525 = t389 * t393;
t243 = -pkin(12) * t525 + t337 + (-pkin(10) * t388 - pkin(5)) * t394;
t297 = t388 * t356 + t377;
t526 = t388 * t389;
t256 = -pkin(12) * t526 + t297;
t163 = t243 * t387 + t256 * t392;
t664 = -qJD(6) * t163 - t387 * t685 + t392 * t686;
t162 = t243 * t392 - t256 * t387;
t663 = qJD(6) * t162 + t387 * t686 + t392 * t685;
t12 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t174;
t38 = Ifges(6,5) * t97 + Ifges(6,6) * t98 + Ifges(6,3) * t188;
t662 = t38 + t12;
t661 = t100 * mrSges(4,2);
t660 = t101 * mrSges(4,1);
t659 = t357 * Ifges(5,5);
t615 = m(7) * pkin(5);
t657 = -mrSges(6,1) - t615;
t361 = t396 * t388;
t362 = t396 * t393;
t262 = t361 * t392 + t362 * t387;
t653 = qJD(6) * t262 + t387 * t683 + t392 * t684;
t263 = t361 * t387 - t362 * t392;
t652 = -qJD(6) * t263 - t387 * t684 + t392 * t683;
t650 = -qJD(5) * t297 + t690;
t143 = mrSges(5,1) * t304 - mrSges(5,3) * t192;
t50 = -mrSges(6,1) * t98 + mrSges(6,2) * t97;
t649 = t50 - t143;
t648 = pkin(5) * t680 - t116;
t434 = t387 * t388 - t392 * t393;
t316 = t434 * t389;
t339 = t387 * t393 + t388 * t392;
t631 = qJD(5) + qJD(6);
t249 = t631 * t339;
t170 = -t249 * t389 - t434 * t510;
t199 = t270 * t387 + t271 * t392;
t647 = t170 - t199;
t171 = t316 * t631 - t339 * t510;
t198 = t270 * t392 - t271 * t387;
t646 = t171 - t198;
t233 = Ifges(6,4) * t235;
t112 = Ifges(6,1) * t236 + Ifges(6,5) * t276 + t233;
t274 = Ifges(5,4) * t285;
t180 = t286 * Ifges(5,1) + t274 + t659;
t645 = t393 * t112 + t180;
t140 = -mrSges(6,1) * t235 + mrSges(6,2) * t236;
t240 = mrSges(5,1) * t357 - mrSges(5,3) * t286;
t643 = t240 - t140;
t508 = qJD(5) * t389;
t641 = t388 * t508 - t393 * t510 + t271;
t640 = pkin(5) * t642 + t677;
t637 = t462 - t511;
t636 = t385 * t491 - t524;
t42 = t153 * t394 - t197 * t510 - t232 * t511 - t389 * t95;
t634 = -t389 * t42 + t394 * t41;
t633 = -t388 * t9 + t393 * t8;
t273 = -t352 * t383 + t367;
t630 = (-t371 * (Ifges(4,5) * t395 - Ifges(4,6) * t390) / 0.2e1 - t273 * (mrSges(4,1) * t390 + mrSges(4,2) * t395)) * t383;
t629 = mrSges(5,1) + t682;
t621 = mrSges(5,2) - t681;
t104 = -pkin(4) * t357 - t115;
t82 = -pkin(5) * t235 + t104;
t628 = -mrSges(7,1) * t82 + mrSges(7,3) * t17;
t627 = mrSges(7,2) * t82 - mrSges(7,3) * t16;
t622 = pkin(10) * t670 + t619;
t72 = -mrSges(7,1) * t469 + mrSges(7,2) * t135;
t620 = m(5) * t115 - m(6) * t104 - m(7) * t82 + t643 - t72;
t584 = -t276 / 0.2e1;
t586 = -t264 / 0.2e1;
t589 = -t236 / 0.2e1;
t591 = -t235 / 0.2e1;
t598 = -t135 / 0.2e1;
t600 = -t469 / 0.2e1;
t617 = Ifges(6,5) * t589 + Ifges(7,5) * t598 + Ifges(6,6) * t591 + Ifges(7,6) * t600 + Ifges(6,3) * t584 + Ifges(7,3) * t586 + t623;
t616 = t383 ^ 2;
t397 = qJD(2) ^ 2;
t39 = t97 * Ifges(6,4) + t98 * Ifges(6,2) + t188 * Ifges(6,6);
t610 = t39 / 0.2e1;
t609 = Ifges(6,1) * t603 + Ifges(6,4) * t602 + Ifges(6,5) * t594;
t554 = Ifges(7,4) * t135;
t62 = Ifges(7,2) * t469 + Ifges(7,6) * t264 + t554;
t608 = -t62 / 0.2e1;
t607 = t62 / 0.2e1;
t131 = Ifges(7,4) * t469;
t63 = Ifges(7,1) * t135 + Ifges(7,5) * t264 + t131;
t606 = -t63 / 0.2e1;
t605 = t63 / 0.2e1;
t604 = Ifges(5,1) * t593 + Ifges(5,4) * t592 + Ifges(5,5) * t578;
t601 = t112 / 0.2e1;
t590 = t235 / 0.2e1;
t588 = t236 / 0.2e1;
t583 = t276 / 0.2e1;
t581 = t285 / 0.2e1;
t579 = t286 / 0.2e1;
t572 = pkin(2) * t383;
t571 = pkin(5) * t236;
t569 = pkin(10) * t394;
t542 = sin(pkin(13));
t451 = t542 * t575;
t543 = cos(pkin(13));
t471 = t543 * t391;
t324 = t386 * t471 + t451;
t452 = t543 * t575;
t470 = t542 * t391;
t403 = -t386 * t452 + t470;
t401 = t403 * t390;
t473 = t384 * t543;
t454 = t383 * t473;
t201 = t324 * t395 - t385 * t401 - t390 * t454;
t246 = t383 * t403 - t385 * t473;
t126 = t201 * t394 + t246 * t389;
t400 = t403 * t395;
t200 = t324 * t390 + t385 * t400 + t395 * t454;
t564 = (-t126 * t380 + t200 * t381) * mrSges(7,1) + (-t126 * t381 - t200 * t380) * mrSges(7,2);
t325 = -t386 * t470 + t452;
t404 = t386 * t451 + t471;
t472 = t384 * t542;
t453 = t383 * t472;
t203 = t325 * t395 + (-t385 * t404 + t453) * t390;
t247 = t383 * t404 + t385 * t472;
t128 = t203 * t394 + t247 * t389;
t402 = t404 * t395;
t202 = t325 * t390 + t385 * t402 - t395 * t453;
t563 = (-t128 * t380 + t202 * t381) * mrSges(7,1) + (-t128 * t381 - t202 * t380) * mrSges(7,2);
t562 = mrSges(6,3) * t235;
t561 = mrSges(6,3) * t236;
t560 = Ifges(4,4) * t390;
t559 = Ifges(4,4) * t395;
t558 = Ifges(5,4) * t389;
t557 = Ifges(5,4) * t394;
t556 = Ifges(6,4) * t388;
t555 = Ifges(6,4) * t393;
t553 = pkin(5) * qJD(6);
t552 = t115 * mrSges(5,3);
t550 = t236 * Ifges(6,4);
t35 = -pkin(4) * t304 - t42;
t548 = t35 * t389;
t538 = t324 * t383;
t537 = t325 * t383;
t534 = t383 * t389;
t111 = t235 * Ifges(6,2) + t276 * Ifges(6,6) + t550;
t527 = t388 * t111;
t413 = t385 * t492 + t523;
t245 = t384 * t413 + t386 * t533;
t416 = -t383 * t493 + t386 * t385;
t205 = t245 * t394 + t389 * t416;
t244 = -t384 * t636 - t386 * t531;
t519 = (-t205 * t380 + t244 * t381) * mrSges(7,1) + (-t205 * t381 - t244 * t380) * mrSges(7,2);
t499 = t383 * t530;
t517 = pkin(2) * t493 + pkin(9) * t499;
t498 = Ifges(5,5) * t192 + Ifges(5,6) * t193 + Ifges(5,3) * t304;
t497 = t288 * pkin(3) + t517;
t496 = Ifges(4,5) * t320 - Ifges(4,6) * t319 + Ifges(4,3) * t370;
t482 = t533 / 0.2e1;
t481 = -t527 / 0.2e1;
t478 = -t515 / 0.2e1;
t474 = -t508 / 0.2e1;
t219 = -t389 * t300 + t301 * t394;
t467 = mrSges(4,3) * t487;
t463 = t383 * t486;
t457 = t395 * t478;
t456 = -t403 * pkin(2) + pkin(9) * t538;
t455 = -t404 * pkin(2) + pkin(9) * t537;
t207 = pkin(4) * t531 - t219;
t445 = Ifges(5,1) * t394 - t558;
t444 = Ifges(6,1) * t393 - t556;
t443 = Ifges(6,1) * t388 + t555;
t442 = -Ifges(5,2) * t389 + t557;
t441 = -Ifges(6,2) * t388 + t555;
t440 = Ifges(6,2) * t393 + t556;
t439 = Ifges(5,5) * t394 - Ifges(5,6) * t389;
t438 = Ifges(6,5) * t393 - Ifges(6,6) * t388;
t437 = Ifges(6,5) * t388 + Ifges(6,6) * t393;
t123 = -t205 * t388 + t244 * t393;
t124 = t205 * t393 + t244 * t388;
t65 = t123 * t392 - t124 * t387;
t66 = t123 * t387 + t124 * t392;
t166 = t254 * t392 + t387 * t427;
t167 = t254 * t387 - t392 * t427;
t226 = -t324 * t529 - t400;
t432 = t226 * pkin(3) + t456;
t228 = -t325 * t529 - t402;
t431 = t228 * pkin(3) + t455;
t130 = -t300 * t510 - t301 * t511 + t312 * t394 - t389 * t313;
t429 = t12 + t673;
t425 = t104 * t446;
t423 = (-mrSges(4,1) * t395 + mrSges(4,2) * t390) * t383;
t422 = (t395 * Ifges(4,2) + t560) * t383;
t410 = t390 * t616 * (Ifges(4,1) * t395 - t560);
t120 = -pkin(4) * t484 - t130;
t204 = t245 * t389 - t394 * t416;
t364 = Ifges(4,4) * t487;
t350 = t495 * t389;
t315 = t339 * t389;
t310 = qJD(2) * t423;
t309 = -mrSges(4,2) * t371 + t467;
t296 = -t388 * t569 + t337;
t258 = Ifges(4,1) * t488 + t371 * Ifges(4,5) + t364;
t257 = t371 * Ifges(4,6) + qJD(2) * t422;
t251 = mrSges(4,1) * t370 - mrSges(4,3) * t320;
t250 = -mrSges(4,2) * t370 - mrSges(4,3) * t319;
t239 = -mrSges(5,2) * t357 + mrSges(5,3) * t285;
t234 = mrSges(4,1) * t319 + mrSges(4,2) * t320;
t227 = t325 * t528 - t390 * t404;
t225 = t324 * t528 - t401;
t195 = t434 * t285;
t194 = t339 * t285;
t183 = t386 * t483 + (t415 * qJD(2) + qJD(3) * t636) * t384;
t182 = t386 * t484 + (qJD(2) * t414 + qJD(3) * t413) * t384;
t178 = t286 * Ifges(5,5) + t285 * Ifges(5,6) + t357 * Ifges(5,3);
t165 = mrSges(6,1) * t276 - t561;
t164 = -mrSges(6,2) * t276 + t562;
t144 = -mrSges(5,2) * t304 + mrSges(5,3) * t193;
t136 = -pkin(5) * t254 + t207;
t107 = mrSges(7,1) * t264 - mrSges(7,3) * t135;
t106 = -mrSges(7,2) * t264 + mrSges(7,3) * t469;
t99 = -mrSges(5,1) * t193 + mrSges(5,2) * t192;
t89 = -qJD(4) * t204 + t183 * t394 + t389 * t463;
t86 = t192 * Ifges(5,4) + t193 * Ifges(5,2) + t304 * Ifges(5,6);
t71 = -pkin(5) * t148 + t120;
t68 = -mrSges(6,2) * t188 + mrSges(6,3) * t98;
t67 = mrSges(6,1) * t188 - mrSges(6,3) * t97;
t57 = -qJD(6) * t167 - t147 * t387 + t148 * t392;
t56 = qJD(6) * t166 + t147 * t392 + t148 * t387;
t46 = qJD(5) * t123 + t182 * t388 + t393 * t89;
t45 = -qJD(5) * t124 + t182 * t393 - t388 * t89;
t23 = -mrSges(7,2) * t174 + mrSges(7,3) * t32;
t22 = mrSges(7,1) * t174 - mrSges(7,3) * t31;
t20 = -pkin(5) * t98 + t35;
t19 = t392 * t48 - t547;
t18 = -t387 * t48 - t545;
t15 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t11 = -qJD(6) * t66 - t387 * t46 + t392 * t45;
t10 = qJD(6) * t65 + t387 * t45 + t392 * t46;
t1 = [m(6) * (t123 * t9 + t124 * t8 + t45 * t54 + t46 * t55) + m(7) * (t10 * t17 + t11 * t16 + t2 * t66 + t3 * t65) + t674 * t182 + (-m(4) * t101 + m(5) * t96 - t251 + t99) * t244 + (-qJDD(2) * t530 - t397 * t493) * mrSges(3,2) + (qJDD(2) * t493 - t397 * t530) * mrSges(3,1) + (-m(2) - m(3) - m(4) - t632) * g(3) - t620 * (qJD(4) * t205 + t183 * t389 - t394 * t463) + (-m(5) * t42 + m(6) * t35 + m(7) * t20 + t15 + t649) * t204 + t310 * t463 + t416 * t234 + m(4) * (t100 * t245 + t218 * t183 + t241 * t416 + t273 * t463) + t183 * t309 + t245 * t250 + t89 * t239 + t205 * t144 + t46 * t164 + t45 * t165 + t123 * t67 + t124 * t68 + t65 * t22 + t66 * t23 + m(3) * (t386 ^ 2 * qJDD(1) + (t321 * t575 + t322 * t391) * t384) + m(5) * (t116 * t89 + t205 * t41) + m(2) * qJDD(1) + t10 * t106 + t11 * t107; (-t116 * t484 + t196 * t252 + t327 * t96 + t41 * t531) * mrSges(5,2) + (t100 * t531 - t101 * t533 - t217 * t483 - t218 * t484) * mrSges(4,3) + (Ifges(5,5) * t327 - Ifges(5,3) * t531) * t578 + t357 * (Ifges(5,5) * t252 + Ifges(5,3) * t484) / 0.2e1 + (t100 * t331 + t101 * t330 + t218 * t639 - t241 * t572 - t273 * t464) * m(4) + t35 * (-mrSges(6,1) * t254 - mrSges(6,2) * t427) + (-t147 * t54 + t148 * t55 + t254 * t8 + t427 * t9) * mrSges(6,3) + (-Ifges(6,1) * t427 + Ifges(6,4) * t254) * t603 + (-Ifges(6,4) * t427 + Ifges(6,2) * t254) * t602 + (-Ifges(6,5) * t427 + Ifges(6,6) * t254) * t594 - t427 * t609 + (t178 * t482 - t630) * qJD(3) + (Ifges(6,1) * t147 + Ifges(6,4) * t148) * t588 + (Ifges(6,4) * t147 + Ifges(6,2) * t148) * t590 + t639 * t309 + t644 * t239 + (t616 * qJD(2) * (-Ifges(4,2) * t390 + t559) + t383 * t258) * t512 / 0.2e1 + t620 * (t266 * t389 - t394 * t464) + (Ifges(6,5) * t147 + Ifges(6,6) * t148) * t583 - t385 * t661 + (Ifges(7,1) * t167 + Ifges(7,4) * t166) * t612 + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t597 + t654 * t164 + (t102 * t9 + t103 * t8 + t104 * t120 + t207 * t35 + t54 * t655 + t55 * t654) * m(6) + t655 * t165 + (-t16 * t56 + t166 * t2 - t167 * t3 + t17 * t57) * mrSges(7,3) + t320 * (Ifges(4,5) * t385 + (t390 * Ifges(4,1) + t559) * t383) / 0.2e1 + (Ifges(4,4) * t320 - Ifges(4,2) * t319 + Ifges(4,6) * t370) * t531 / 0.2e1 + (Ifges(4,1) * t320 - Ifges(4,4) * t319 + Ifges(4,5) * t370) * t482 + (-m(4) * t455 - m(7) * t431 + mrSges(3,1) * t404 - t228 * mrSges(4,1) + t325 * mrSges(3,2) - mrSges(4,3) * t537 + t670 * (pkin(10) * t227 + t431) - t629 * (t228 * t394 + t325 * t534) + t619 * t227 + t621 * (t228 * t389 - t325 * t532)) * g(1) + (-m(4) * t456 - m(7) * t432 + mrSges(3,1) * t403 - t226 * mrSges(4,1) + t324 * mrSges(3,2) - mrSges(4,3) * t538 + t670 * (pkin(10) * t225 + t432) - t629 * (t226 * t394 + t324 * t534) + t619 * t225 + t621 * (t226 * t389 - t324 * t532)) * g(2) + (-(mrSges(3,1) * t575 - mrSges(3,2) * t391) * t384 - m(7) * t497 - m(4) * t517 - t288 * mrSges(4,1) - mrSges(4,3) * t499 + t670 * (pkin(10) * t287 + t497) - t629 * (t288 * t394 + t389 * t499) + t619 * t287 + t621 * (t288 * t389 - t394 * t499)) * g(3) + t385 * t660 + (Ifges(5,4) * t252 + Ifges(5,6) * t484) * t581 + (Ifges(5,4) * t327 - Ifges(5,6) * t531) * t592 - t257 * t484 / 0.2e1 + t115 * (mrSges(5,1) * t484 - mrSges(5,3) * t252) + (t115 * t130 + t116 * t644 + t219 * t42 + t220 * t41 + t299 * t96) * m(5) - t310 * t464 - t234 * t572 + t385 * t496 / 0.2e1 - t498 * t531 / 0.2e1 + t42 * (-mrSges(5,1) * t531 - mrSges(5,3) * t327) + (t363 - t322) * mrSges(3,2) + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t611 + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t599 + t666 * t106 + (t136 * t20 + t16 * t667 + t17 * t666 + t2 * t37 + t3 * t36 + t71 * t82) * m(7) + t667 * t107 + (Ifges(5,1) * t252 + Ifges(5,5) * t484) * t579 + (Ifges(5,1) * t327 - Ifges(5,5) * t531) * t593 + (t459 + t321) * mrSges(3,1) + t410 * t506 / 0.2e1 - t319 * (Ifges(4,6) * t385 + t422) / 0.2e1 + Ifges(3,3) * qJDD(2) + t241 * t423 + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t585 + (Ifges(7,5) * t167 + Ifges(7,6) * t166) * t596 + t674 * t694 + t330 * t251 + t331 * t250 + t299 * t99 + t252 * t180 / 0.2e1 + t130 * t240 + t219 * t143 + t220 * t144 + t207 * t50 + t20 * (-mrSges(7,1) * t166 + mrSges(7,2) * t167) + t104 * (-mrSges(6,1) * t148 + mrSges(6,2) * t147) + t148 * t111 / 0.2e1 + t120 * t140 + t136 * t15 + t57 * t607 + t254 * t610 + t167 * t613 + t166 * t614 + t370 * (Ifges(4,3) * t385 + (Ifges(4,5) * t390 + Ifges(4,6) * t395) * t383) / 0.2e1 + t36 * t22 + t37 * t23 + (Ifges(6,3) * t594 + Ifges(6,6) * t602 + Ifges(6,5) * t603 - Ifges(5,6) * t578 - t41 * mrSges(5,3) + t662 / 0.2e1 - t86 / 0.2e1 + t96 * mrSges(5,1) + Ifges(7,6) * t611 + Ifges(7,5) * t612 - Ifges(5,2) * t592 - Ifges(5,4) * t593 + Ifges(7,3) * t596 + t673 + t679) * t326 + (t656 / 0.2e1 - Ifges(5,4) * t579 - Ifges(5,2) * t581 + Ifges(6,3) * t583 + Ifges(6,5) * t588 + Ifges(6,6) * t590 - t678 + t689) * t253 + t147 * t601 + t327 * t604 + t56 * t605 + t71 * t72 + t82 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t102 * t67 + t103 * t68; (-mrSges(7,1) * t646 + mrSges(7,2) * t647) * t82 + (-mrSges(7,1) * t637 - mrSges(7,3) * t647) * t16 + (Ifges(7,5) * t170 + Ifges(7,6) * t171) * t585 + (Ifges(7,5) * t199 + Ifges(7,6) * t198) * t586 + (t467 - t309) * t217 + (t200 * t676 + t201 * t622) * g(2) + (t202 * t676 + t203 * t622) * g(1) + (t244 * t676 + t245 * t622) * g(3) + t677 * t140 + t394 * t669 + (Ifges(6,5) * t271 + Ifges(6,6) * t270) * t584 + (Ifges(7,4) * t170 + Ifges(7,2) * t171) * t599 + (Ifges(7,4) * t199 + Ifges(7,2) * t198) * t600 + (t257 * t482 + t630) * qJD(2) + (Ifges(6,4) * t271 + Ifges(6,2) * t270) * t591 + t649 * pkin(10) * t389 + t640 * t72 + (mrSges(6,1) * t642 - mrSges(6,2) * t641) * t104 + (-t525 * t9 - t526 * t8 + t54 * t641 - t55 * t642) * mrSges(6,3) + t645 * t510 / 0.2e1 + (mrSges(7,2) * t637 + mrSges(7,3) * t646) * t17 + (-t501 - t154) * t240 + (t179 / 0.2e1 + t617) * t462 + t357 * t196 * (mrSges(5,1) * t389 + mrSges(5,2) * t394) + (t285 * t442 + t286 * t445 + t357 * t439) * qJD(4) / 0.2e1 + (-t552 + t481) * t510 + (Ifges(7,1) * t170 + Ifges(7,4) * t171) * t597 + (Ifges(7,1) * t199 + Ifges(7,4) * t198) * t598 + (Ifges(6,1) * t271 + Ifges(6,4) * t270) * t589 + t650 * t165 + t651 * t164 + (-t104 * t141 + t296 * t9 + t297 * t8 + (t104 * t510 + t548) * pkin(10) + t651 * t55 + t650 * t54) * m(6) + t496 + t660 - t661 + t2 * (mrSges(7,2) * t394 - mrSges(7,3) * t315) + t20 * (mrSges(7,1) * t315 - mrSges(7,2) * t316) + (-Ifges(7,4) * t316 - Ifges(7,2) * t315 - Ifges(7,6) * t394) * t611 + (-Ifges(7,1) * t316 - Ifges(7,4) * t315 - Ifges(7,5) * t394) * t612 + (-Ifges(7,5) * t316 - Ifges(7,6) * t315 - Ifges(7,3) * t394) * t596 + t3 * (-mrSges(7,1) * t394 + mrSges(7,3) * t316) + (-t270 / 0.2e1 + t393 * t474) * t111 + (-Ifges(4,2) * t488 + t394 * t180 + t258 + t364) * t457 + t144 * t569 + (Ifges(5,5) * t389 + Ifges(5,6) * t394) * t578 + (-t437 * t508 + (Ifges(6,3) * t389 + t394 * t438) * qJD(4)) * t583 + (t468 + t624) * t218 + (-t271 / 0.2e1 + t388 * t474) * t112 - t397 * t410 / 0.2e1 + t394 * t86 / 0.2e1 + (-t116 * (-mrSges(5,3) * t389 * t395 - mrSges(5,2) * t390) - t115 * (mrSges(5,1) * t390 - mrSges(5,3) * t520)) * t515 - t662 * t394 / 0.2e1 + t663 * t106 + (t16 * t664 + t162 * t3 + t163 * t2 + t17 * t663 + t20 * t350 + t640 * t82) * m(7) + t664 * t107 + (t357 * (Ifges(5,3) * t390 + t395 * t439) + t286 * (Ifges(5,5) * t390 + t395 * t445) + t285 * (Ifges(5,6) * t390 + t395 * t442) + t390 * t178) * t478 + t96 * t448 - t394 * t668 - t39 * t526 / 0.2e1 + t446 * t548 + (-t115 * t154 - t116 * t155 - pkin(3) * t96 + ((-t115 * t394 - t116 * t389) * qJD(4) + t634) * pkin(10)) * m(5) + t634 * mrSges(5,3) + t656 * (t389 * t457 + t511 / 0.2e1) + (-t502 - t155) * t239 + t350 * t15 + t296 * t67 + t297 * t68 + t162 * t22 + t163 * t23 + t171 * t607 + t198 * t608 + t525 * t609 - t316 * t613 - t315 * t614 + (-t443 * t508 + (Ifges(6,5) * t389 + t394 * t444) * qJD(4)) * t588 + (-t440 * t508 + (Ifges(6,6) * t389 + t394 * t441) * qJD(4)) * t590 + (Ifges(5,2) * t394 + t558) * t592 + (Ifges(5,1) * t389 + t557) * t593 + (-Ifges(6,3) * t394 + t389 * t438) * t594 + (-t551 + t689) * t511 + (-Ifges(6,6) * t394 + t389 * t441) * t602 + (-Ifges(6,5) * t394 + t389 * t444) * t603 + t389 * t604 + t170 * t605 + t199 * t606 - pkin(3) * t99; (-pkin(4) * t35 - t104 * t116 - t54 * t78 - t55 * t79) * m(6) + (-Ifges(7,5) * t195 - Ifges(7,6) * t194) * t586 + (-Ifges(7,1) * t195 - Ifges(7,4) * t194) * t598 + (-Ifges(7,4) * t195 - Ifges(7,2) * t194) * t600 - t82 * (mrSges(7,1) * t194 - mrSges(7,2) * t195) + (t204 * t629 + t205 * t621) * g(3) + (t621 * t126 - t629 * (-t201 * t389 + t246 * t394)) * g(2) + (t621 * t128 - t629 * (-t203 * t389 + t247 * t394)) * g(1) + t643 * t116 - (-Ifges(5,2) * t286 + t274 + t645) * t285 / 0.2e1 + (t235 * t441 + t236 * t444 + t276 * t438) * qJD(5) / 0.2e1 + t648 * t72 + t652 * t107 + (t16 * t652 + t17 * t653 + t2 * t263 - t20 * t379 + t262 * t3 + t648 * t82) * m(7) + t653 * t106 - (Ifges(5,1) * t285 - t549 + t656) * t286 / 0.2e1 + (-t425 + t438 * t584 + t552 - t659 / 0.2e1 - t196 * mrSges(5,2) + t444 * t589 + t441 * t591) * t285 + t179 * t579 + t527 * t581 + (-t2 * mrSges(7,3) + t20 * mrSges(7,1) - 0.2e1 * t614 - (Ifges(7,1) * t597 + Ifges(7,4) * t599 + Ifges(7,5) * t585 + t605 + t627) * t631) * t434 + (mrSges(7,2) * t20 - mrSges(7,3) * t3 + 0.2e1 * t613) * t339 + (-t16 * t195 + t17 * t194) * mrSges(7,3) + t35 * t447 + (t481 + t425) * qJD(5) + (m(6) * ((-t388 * t55 - t393 * t54) * qJD(5) + t633) - t164 * t509 - t165 * t507 - t388 * t67 + t393 * t68) * pkin(11) - t379 * t15 + t262 * t22 + t263 * t23 - t115 * t239 - t79 * t164 - t78 * t165 - t194 * t608 + t388 * t609 + t393 * t610 - (Ifges(7,4) * t597 + Ifges(7,2) * t599 + Ifges(7,6) * t585 + t607 + t628) * t249 + (t617 + t678) * t286 - t41 * mrSges(5,2) + t42 * mrSges(5,1) + (-t680 * t55 + (-t507 + t539) * t54 + t633) * mrSges(6,3) + t437 * t594 + t507 * t601 + t440 * t602 + t443 * t603 - t195 * t606 - pkin(4) * t50 + t498; (mrSges(6,2) * t124 + t123 * t657 - t519) * g(3) + (-t564 - (-t126 * t393 - t200 * t388) * mrSges(6,2) + t657 * (-t126 * t388 + t200 * t393)) * g(2) + (-t563 - (-t128 * t393 - t202 * t388) * mrSges(6,2) + t657 * (-t128 * t388 + t202 * t393)) * g(1) + (t392 * t553 - t19) * t106 + t679 + (Ifges(6,5) * t235 - Ifges(6,6) * t236) * t584 - t72 * t571 - m(7) * (t16 * t18 + t17 * t19 + t571 * t82) + (-t387 * t553 - t18) * t107 + (t22 * t392 + t23 * t387) * pkin(5) + t38 - t104 * (mrSges(6,1) * t236 + mrSges(6,2) * t235) + (t2 * t387 + t3 * t392 + (-t16 * t387 + t17 * t392) * qJD(6)) * t615 + (Ifges(7,1) * t598 + Ifges(7,4) * t600 + Ifges(7,5) * t586 + t606 - t627) * t469 - (Ifges(7,4) * t598 + Ifges(7,2) * t600 + Ifges(7,6) * t586 + t608 - t628) * t135 + t111 * t588 + (Ifges(6,1) * t235 - t550) * t589 + (t562 - t164) * t54 + (t561 + t165) * t55 + (-Ifges(6,2) * t236 + t112 + t233) * t591 + t429; -t82 * (mrSges(7,1) * t135 + mrSges(7,2) * t469) + (Ifges(7,1) * t469 - t554) * t598 + t62 * t597 + (Ifges(7,5) * t469 - Ifges(7,6) * t135) * t586 - t16 * t106 + t17 * t107 - g(1) * t563 - g(2) * t564 - g(3) * t519 + (t135 * t17 + t16 * t469) * mrSges(7,3) + t429 + (-Ifges(7,2) * t135 + t131 + t63) * t600;];
tau  = t1;
