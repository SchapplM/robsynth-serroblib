% Calculate vector of inverse dynamics joint torques for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR15_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:24:03
% EndTime: 2024-09-27 22:24:23
% DurationCPUTime: 11.70s
% Computational Cost: add. (22997->808), mult. (69258->1153), div. (0->0), fcn. (56140->26), ass. (0->433)
t396 = cos(pkin(5));
t400 = sin(qJ(2));
t501 = t396 * t400;
t368 = pkin(1) * t501;
t405 = cos(qJ(2));
t394 = sin(pkin(5));
t572 = pkin(8) + pkin(9);
t596 = t394 * t572;
t609 = t405 * t596 + t368;
t271 = t609 * qJD(1);
t404 = cos(qJ(3));
t255 = t404 * t271;
t558 = pkin(1) * t396;
t369 = t405 * t558;
t362 = qJD(1) * t369;
t441 = t400 * t596;
t270 = -qJD(1) * t441 + t362;
t399 = sin(qJ(3));
t194 = -t270 * t399 - t255;
t498 = t399 * t400;
t420 = t404 * t405 - t498;
t295 = t420 * t394;
t285 = qJD(1) * t295;
t551 = pkin(10) * t285;
t161 = t194 - t551;
t252 = t399 * t271;
t195 = t404 * t270 - t252;
t496 = t399 * t405;
t421 = t400 * t404 + t496;
t296 = t421 * t394;
t286 = qJD(1) * t296;
t282 = t286 * pkin(10);
t162 = -t282 + t195;
t557 = pkin(2) * t404;
t381 = pkin(3) + t557;
t398 = sin(qJ(4));
t403 = cos(qJ(4));
t472 = qJD(4) * t403;
t473 = qJD(4) * t398;
t499 = t398 * t399;
t592 = -t398 * t161 - t403 * t162 + t381 * t472 + (-t399 * t473 + (t403 * t404 - t499) * qJD(3)) * pkin(2);
t213 = t285 * t398 + t286 * t403;
t395 = cos(pkin(6));
t550 = pkin(11) * t395;
t203 = t213 * t550;
t617 = t203 + t592;
t581 = m(4) * pkin(2);
t616 = -t581 - mrSges(3,1);
t397 = sin(qJ(5));
t402 = cos(qJ(5));
t366 = t396 * qJD(1) + qJD(2);
t359 = qJD(3) + t366;
t344 = qJD(4) + t359;
t393 = sin(pkin(6));
t443 = t403 * t285 - t286 * t398;
t519 = t443 * t395;
t427 = t344 * t393 + t519;
t127 = -t213 * t397 + t402 * t427;
t128 = t213 * t402 + t397 * t427;
t520 = t443 * t393;
t181 = t344 * t395 + qJD(5) - t520;
t382 = t405 * pkin(2) + pkin(1);
t336 = t382 * t394;
t326 = qJD(1) * t336;
t229 = -pkin(3) * t285 - t326;
t613 = t213 * t393;
t118 = -pkin(4) * t443 - pkin(11) * t613 + t229;
t233 = pkin(2) * t366 + t270;
t182 = t404 * t233 - t252;
t151 = t182 - t282;
t145 = pkin(3) * t359 + t151;
t183 = t233 * t399 + t255;
t152 = t183 + t551;
t146 = t398 * t152;
t91 = t403 * t145 - t146;
t73 = t91 - t203;
t71 = pkin(4) * t344 + t73;
t433 = t118 * t393 + t395 * t71;
t148 = t403 * t152;
t92 = t145 * t398 + t148;
t69 = pkin(11) * t427 + t92;
t24 = -t397 * t69 + t402 * t433;
t25 = t397 * t433 + t402 * t69;
t469 = qJD(1) * qJD(2);
t311 = (qJDD(1) * t405 - t400 * t469) * t394;
t312 = (qJDD(1) * t400 + t405 * t469) * t394;
t187 = qJD(3) * t285 + t311 * t399 + t312 * t404;
t383 = t396 * qJDD(1);
t364 = t383 + qJDD(2);
t349 = qJDD(3) + t364;
t464 = qJD(2) * t558;
t438 = qJD(1) * t464;
t460 = pkin(1) * t383;
t228 = -pkin(8) * t312 - t400 * t438 + t405 * t460;
t186 = pkin(2) * t364 - pkin(9) * t312 + t228;
t504 = t394 * t405;
t367 = pkin(8) * t504;
t505 = t394 * t400;
t447 = qJD(2) * t505;
t442 = pkin(8) * t447;
t227 = -qJD(1) * t442 + qJDD(1) * t367 + t400 * t460 + t405 * t438;
t197 = pkin(9) * t311 + t227;
t90 = -qJD(3) * t183 + t404 * t186 - t197 * t399;
t67 = pkin(3) * t349 - pkin(10) * t187 + t90;
t412 = t421 * qJD(3);
t477 = qJD(1) * t394;
t188 = t311 * t404 - t312 * t399 - t412 * t477;
t474 = qJD(3) * t404;
t475 = qJD(3) * t399;
t89 = t399 * t186 + t404 * t197 + t233 * t474 - t271 * t475;
t70 = pkin(10) * t188 + t89;
t20 = t145 * t472 - t152 * t473 + t398 * t67 + t403 * t70;
t343 = qJDD(4) + t349;
t98 = -qJD(4) * t213 - t187 * t398 + t188 * t403;
t428 = t343 * t393 + t395 * t98;
t17 = pkin(11) * t428 + t20;
t21 = -qJD(4) * t92 - t398 * t70 + t403 * t67;
t97 = qJD(4) * t443 + t187 * t403 + t188 * t398;
t18 = pkin(4) * t343 - t550 * t97 + t21;
t559 = pkin(1) * t394;
t274 = -pkin(2) * t311 - qJDD(1) * t559;
t149 = -pkin(3) * t188 + t274;
t389 = t393 * pkin(11);
t44 = -pkin(4) * t98 - t389 * t97 + t149;
t435 = t18 * t395 + t393 * t44;
t597 = qJD(5) * t24;
t3 = t17 * t402 + t397 * t435 + t597;
t4 = -qJD(5) * t25 - t17 * t397 + t435 * t402;
t506 = t393 * t402;
t465 = mrSges(6,3) * t506;
t471 = qJD(5) * t393;
t48 = t118 * t395 - t393 * t71;
t507 = t393 * t397;
t534 = Ifges(6,4) * t402;
t535 = Ifges(6,4) * t397;
t8 = -t18 * t393 + t395 * t44;
t546 = t393 * t8;
t568 = -t181 / 0.2e1;
t570 = -t128 / 0.2e1;
t571 = -t127 / 0.2e1;
t85 = t343 * t395 - t393 * t98 + qJDD(5);
t573 = t85 / 0.2e1;
t536 = Ifges(6,4) * t128;
t63 = Ifges(6,2) * t127 + Ifges(6,6) * t181 + t536;
t575 = -t63 / 0.2e1;
t41 = -qJD(5) * t128 - t397 * t97 + t428 * t402;
t576 = t41 / 0.2e1;
t40 = qJD(5) * t127 + t397 * t428 + t402 * t97;
t577 = t40 / 0.2e1;
t578 = Ifges(6,1) * t577 + Ifges(6,4) * t576 + Ifges(6,5) * t573;
t579 = Ifges(6,4) * t577 + Ifges(6,2) * t576 + Ifges(6,6) * t573;
t12 = Ifges(6,5) * t40 + Ifges(6,6) * t41 + Ifges(6,3) * t85;
t580 = t12 / 0.2e1;
t584 = t21 * mrSges(5,1) - t20 * mrSges(5,2) + Ifges(5,5) * t97 + Ifges(5,6) * t98 + Ifges(5,3) * t343;
t503 = t395 * t397;
t610 = -t213 * t503 + t402 * t443;
t502 = t395 * t402;
t611 = -t213 * t502 - t397 * t443;
t125 = Ifges(6,4) * t127;
t64 = Ifges(6,1) * t128 + Ifges(6,5) * t181 + t125;
t615 = (-mrSges(6,1) * t402 + mrSges(6,2) * t397) * t546 + (Ifges(6,3) * t395 + (Ifges(6,5) * t397 + Ifges(6,6) * t402) * t393) * t573 + (Ifges(6,6) * t395 + (Ifges(6,2) * t402 + t535) * t393) * t576 + (Ifges(6,5) * t395 + (Ifges(6,1) * t397 + t534) * t393) * t577 + t507 * t578 + t506 * t579 + t395 * t580 + t3 * (-mrSges(6,2) * t395 + t465) + t584 + t4 * (mrSges(6,1) * t395 - mrSges(6,3) * t507) - t465 * t597 + (-mrSges(6,3) * t25 + t575) * t397 * t471 + (t127 * (-Ifges(6,2) * t397 + t534) + t128 * (Ifges(6,1) * t402 - t535) + t181 * (Ifges(6,5) * t402 - Ifges(6,6) * t397) + t402 * t64) * t471 / 0.2e1 - t610 * t64 / 0.2e1 + t611 * t575 + (Ifges(6,1) * t610 + Ifges(6,4) * t611 + Ifges(6,5) * t613) * t570 + (Ifges(6,4) * t610 + Ifges(6,2) * t611 + Ifges(6,6) * t613) * t571 + (Ifges(6,5) * t610 + Ifges(6,6) * t611 + Ifges(6,3) * t613) * t568 - t24 * (mrSges(6,1) * t613 - mrSges(6,3) * t610) - t25 * (-mrSges(6,2) * t613 + mrSges(6,3) * t611) + ((mrSges(6,1) * t397 + mrSges(6,2) * t402) * t471 + mrSges(6,1) * t611 - mrSges(6,2) * t610) * t48;
t614 = pkin(4) * t213;
t209 = Ifges(5,4) * t443;
t133 = Ifges(5,1) * t213 + Ifges(5,5) * t344 + t209;
t612 = t133 + t209;
t608 = Ifges(6,5) * t570 + Ifges(6,6) * t571 + Ifges(6,3) * t568;
t537 = Ifges(5,4) * t213;
t132 = Ifges(5,2) * t443 + Ifges(5,6) * t344 + t537;
t607 = t132 / 0.2e1;
t567 = -t443 / 0.2e1;
t606 = mrSges(5,2) * t443;
t543 = mrSges(5,3) * t443;
t605 = Ifges(5,1) * t443;
t604 = Ifges(5,5) * t443;
t553 = pkin(3) * t398;
t346 = t389 + t553;
t552 = pkin(3) * t403;
t379 = pkin(4) + t552;
t277 = t346 * t402 + t379 * t503;
t556 = pkin(3) * t286;
t124 = -t389 * t443 + t556 + t614;
t468 = t443 * t550;
t99 = -t151 * t398 - t148;
t76 = t99 - t468;
t431 = t124 * t393 + t395 * t76;
t533 = pkin(3) * qJD(4);
t100 = t403 * t151 - t146;
t77 = -t203 + t100;
t603 = t397 * t77 - t402 * t431 - t277 * qJD(5) + (-t397 * t403 - t398 * t502) * t533;
t275 = -t346 * t397 + t379 * t502;
t602 = -t397 * t431 - t402 * t77 + t275 * qJD(5) + (-t398 * t503 + t402 * t403) * t533;
t329 = pkin(4) * t503 + pkin(11) * t506;
t140 = -pkin(11) * t520 + t614;
t74 = -pkin(11) * t519 - t92;
t429 = t140 * t393 + t395 * t74;
t601 = -t329 * qJD(5) + t397 * t73 - t402 * t429;
t327 = pkin(4) * t502 - pkin(11) * t507;
t600 = t327 * qJD(5) - t397 * t429 - t402 * t73;
t497 = t399 * t403;
t317 = pkin(2) * t497 + t398 * t381;
t287 = t389 + t317;
t316 = -pkin(2) * t499 + t403 * t381;
t313 = pkin(4) + t316;
t215 = t287 * t402 + t313 * t503;
t231 = -t381 * t473 + (-t399 * t472 + (-t398 * t404 - t497) * qJD(3)) * pkin(2);
t449 = t400 * t477;
t360 = pkin(2) * t449;
t119 = t124 + t360;
t101 = t403 * t161 - t162 * t398;
t79 = t101 - t468;
t432 = t119 * t393 + t395 * t79;
t599 = -qJD(5) * t215 + t231 * t502 - t397 * t617 - t402 * t432;
t214 = -t287 * t397 + t313 * t502;
t598 = qJD(5) * t214 + t231 * t503 - t397 * t432 + t402 * t617;
t593 = -t101 + t231;
t269 = pkin(2) * t396 + t369 - t441;
t478 = t367 + t368;
t284 = pkin(9) * t504 + t478;
t200 = t404 * t269 - t284 * t399;
t554 = pkin(3) * t396;
t163 = -pkin(10) * t296 + t200 + t554;
t201 = t399 * t269 + t404 * t284;
t171 = pkin(10) * t295 + t201;
t111 = t398 * t163 + t403 * t171;
t339 = -pkin(4) * t398 + t389 * t403;
t325 = t339 * t404;
t338 = -pkin(4) * t403 - t389 * t398;
t333 = pkin(3) - t338;
t239 = -t333 * t399 + t325;
t232 = t239 * t405;
t511 = t339 * t399;
t238 = -t333 * t404 - t511;
t236 = pkin(2) - t238;
t174 = -t236 * t400 + t232;
t392 = qJ(2) + qJ(3);
t385 = pkin(5) + t392;
t375 = qJ(4) + t385;
t354 = cos(t375) / 0.2e1;
t386 = pkin(5) - t392;
t436 = -qJ(4) + t386;
t365 = cos(t436);
t390 = qJ(4) + t392;
t377 = sin(t390);
t378 = cos(t390);
t419 = sin(t436);
t539 = mrSges(6,3) * t393;
t466 = t377 * t539;
t467 = sin(t375) / 0.2e1;
t591 = -(t467 + t419 / 0.2e1) * mrSges(5,1) - (t354 - t365 / 0.2e1) * mrSges(5,2) - ((-t377 * t503 + t378 * t402) * mrSges(6,1) + (-t377 * t502 - t378 * t397) * mrSges(6,2) + t466) * t394;
t315 = t365 / 0.2e1 + t354;
t401 = sin(qJ(1));
t406 = cos(qJ(1));
t510 = t377 * t406;
t251 = -t401 * t315 - t510;
t491 = t401 * t397;
t455 = t396 * t491;
t485 = t406 * t402;
t291 = t395 * t455 - t485;
t488 = t402 * t401;
t454 = t396 * t488;
t500 = t397 * t406;
t292 = -t395 * t454 - t500;
t299 = t395 * t500 + t454;
t300 = -t395 * t485 + t455;
t314 = t467 - t419 / 0.2e1;
t494 = t401 * t377;
t509 = t378 * t406;
t590 = -t251 * mrSges(5,1) - (t401 * t314 - t509) * mrSges(5,2) - (t291 * t377 - t299 * t378) * mrSges(6,1) - (-t292 * t377 + t300 * t378) * mrSges(6,2) - (-t396 * t494 + t509) * t539;
t250 = -t315 * t406 + t494;
t457 = t396 * t500;
t293 = t395 * t457 + t488;
t456 = t396 * t485;
t294 = t395 * t456 - t491;
t297 = t395 * t491 - t456;
t298 = t395 * t488 + t457;
t493 = t401 * t378;
t589 = t250 * mrSges(5,1) - (-t314 * t406 - t493) * mrSges(5,2) - (-t293 * t377 - t297 * t378) * mrSges(6,1) - (-t294 * t377 - t298 * t378) * mrSges(6,2) - (t396 * t510 + t493) * t539;
t538 = Ifges(3,4) * t400;
t588 = pkin(1) * (mrSges(3,1) * t400 + mrSges(3,2) * t405) - t400 * (Ifges(3,1) * t405 - t538) / 0.2e1;
t371 = sin(t385);
t372 = sin(t386);
t373 = cos(t385);
t374 = cos(t386);
t587 = -(t371 / 0.2e1 + t372 / 0.2e1) * mrSges(4,1) - (t373 / 0.2e1 - t374 / 0.2e1) * mrSges(4,2) + t591;
t335 = t374 + t373;
t387 = sin(t392);
t492 = t401 * t387;
t560 = -t406 / 0.2e1;
t276 = t335 * t560 + t492;
t334 = t371 - t372;
t388 = cos(t392);
t586 = t276 * mrSges(4,1) - (t334 * t560 - t401 * t388) * mrSges(4,2) + t589;
t486 = t406 * t387;
t278 = -t486 - t401 * t335 / 0.2e1;
t562 = t334 / 0.2e1;
t585 = -t278 * mrSges(4,1) - (-t406 * t388 + t401 * t562) * mrSges(4,2) + t590;
t221 = t295 * t398 + t296 * t403;
t220 = t295 * t403 - t296 * t398;
t426 = t220 * t395 + t393 * t396;
t139 = t221 * t402 + t397 * t426;
t555 = pkin(3) * t295;
t237 = -t336 - t555;
t126 = -pkin(4) * t220 - t221 * t389 + t237;
t110 = t403 * t163 - t171 * t398;
t83 = pkin(4) * t396 - t221 * t550 + t110;
t430 = t126 * t393 + t395 * t83;
t82 = pkin(11) * t426 + t111;
t35 = t397 * t430 + t402 * t82;
t583 = t90 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,5) * t187 + Ifges(4,6) * t188 + Ifges(4,3) * t349;
t582 = t228 * mrSges(3,1) - t227 * mrSges(3,2) + Ifges(3,5) * t312 + Ifges(3,6) * t311 + Ifges(3,3) * t364;
t569 = t128 / 0.2e1;
t566 = -t213 / 0.2e1;
t565 = t213 / 0.2e1;
t563 = t286 / 0.2e1;
t561 = -t344 / 0.2e1;
t549 = g(1) * t406;
t548 = g(2) * t401;
t544 = mrSges(4,3) * t285;
t542 = mrSges(5,3) * t213;
t541 = mrSges(6,3) * t127;
t540 = mrSges(6,3) * t128;
t15 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t532 = t15 * t393;
t529 = t286 * mrSges(4,3);
t528 = t286 * Ifges(4,4);
t527 = t366 * Ifges(3,5);
t526 = t366 * Ifges(3,6);
t72 = -mrSges(6,1) * t127 + mrSges(6,2) * t128;
t525 = t393 * t72;
t222 = (qJD(2) + qJD(3)) * t295;
t223 = (-qJD(2) * t421 - t412) * t394;
t121 = -qJD(4) * t221 - t222 * t398 + t223 * t403;
t522 = t121 * t393;
t516 = t221 * t397;
t513 = t239 * t400;
t508 = t393 * t394;
t495 = t400 * t406;
t490 = t401 * t400;
t489 = t401 * t405;
t487 = t405 * t406;
t470 = pkin(10) + t572;
t463 = pkin(3) * t473;
t459 = t401 * t508;
t458 = t406 * t508;
t453 = t393 * t608;
t448 = t405 * t477;
t204 = pkin(2) * t447 - pkin(3) * t223;
t440 = mrSges(3,3) * t449;
t439 = mrSges(3,3) * t448;
t120 = qJD(4) * t220 + t222 * t403 + t223 * t398;
t363 = t405 * t464;
t272 = -qJD(2) * t441 + t363;
t273 = t609 * qJD(2);
t129 = t269 * t474 + t404 * t272 - t399 * t273 - t284 * t475;
t112 = pkin(10) * t223 + t129;
t130 = -qJD(3) * t201 - t272 * t399 - t404 * t273;
t113 = -pkin(10) * t222 + t130;
t43 = -qJD(4) * t111 - t112 * t398 + t403 * t113;
t27 = -t120 * t550 + t43;
t61 = -pkin(4) * t121 - t120 * t389 + t204;
t434 = t27 * t395 + t393 * t61;
t425 = t236 * t405 + t513;
t245 = t338 * t404 - t511;
t246 = t338 * t399 + t325;
t424 = t245 * t405 - t246 * t400;
t423 = -t291 * t378 - t299 * t377;
t422 = -t294 * t378 + t298 * t377;
t380 = pkin(3) * t404 + pkin(2);
t418 = -pkin(3) * t498 + t380 * t405;
t416 = t396 * t487 - t490;
t322 = -t396 * t489 - t495;
t42 = t403 * t112 + t398 * t113 + t163 * t472 - t171 * t473;
t411 = -m(3) * pkin(8) - mrSges(6,3) * t395 - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t34 = -t397 * t82 + t402 * t430;
t410 = m(4) * (pkin(2) * t501 - t596) + m(5) * (-t394 * t470 + (pkin(3) * t496 + t380 * t400) * t396) - m(6) * (t394 * (t470 + t550) + t174 * t396) + mrSges(4,1) * t562 + t314 * mrSges(5,1) + mrSges(2,2);
t409 = m(3) * pkin(1) + m(4) * t382 + m(5) * (pkin(3) * t388 + t382) + m(6) * (pkin(1) + t425) + mrSges(4,1) * t388 + mrSges(5,1) * t378 + mrSges(2,1) + t466;
t198 = t285 * Ifges(4,2) + t359 * Ifges(4,6) + t528;
t281 = Ifges(4,4) * t285;
t199 = t286 * Ifges(4,1) + t359 * Ifges(4,5) + t281;
t407 = t91 * t543 + t182 * t544 + t198 * t563 + t183 * t529 + t605 * t566 + t604 * t561 - t229 * t606 + t612 * t567 + t583 - t286 * (Ifges(4,1) * t285 - t528) / 0.2e1 - t359 * (Ifges(4,5) * t285 - Ifges(4,6) * t286) / 0.2e1 + t326 * (mrSges(4,1) * t286 + mrSges(4,2) * t285) + t613 * t608 + t615 - (-Ifges(4,2) * t286 + t199 + t281) * t285 / 0.2e1 + (-t229 * mrSges(5,1) + mrSges(5,3) * t92 - Ifges(5,4) * t566 - Ifges(5,2) * t567 - Ifges(5,6) * t561 + t607) * t213;
t358 = Ifges(3,4) * t448;
t347 = t396 * t506;
t341 = -pkin(2) * t400 - pkin(3) * t387;
t328 = -pkin(8) * t505 + t369;
t324 = (-mrSges(3,1) * t405 + mrSges(3,2) * t400) * t394;
t323 = -t396 * t490 + t487;
t321 = -t396 * t495 - t489;
t310 = t478 * qJD(2);
t309 = t363 - t442;
t306 = t478 * qJD(1);
t305 = -pkin(8) * t449 + t362;
t304 = -mrSges(3,2) * t366 + t439;
t303 = mrSges(3,1) * t366 - t440;
t288 = t420 * t554;
t283 = t418 * t396;
t268 = Ifges(3,1) * t449 + t358 + t527;
t267 = t526 + (Ifges(3,2) * t405 + t538) * t477;
t240 = t360 + t556;
t235 = mrSges(4,1) * t359 - t529;
t234 = -mrSges(4,2) * t359 + t544;
t225 = t239 * t505;
t219 = -mrSges(4,1) * t285 + mrSges(4,2) * t286;
t202 = -t220 * t393 + t395 * t396;
t192 = -t293 * t378 + t297 * t377 + t397 * t458;
t191 = t292 * t378 + t300 * t377 + t402 * t459;
t190 = mrSges(5,1) * t344 - t542;
t189 = -mrSges(5,2) * t344 + t543;
t176 = t245 * t400 + t246 * t405;
t175 = t238 * t400 + t232;
t172 = t424 * t396;
t170 = (t238 * t405 - t513) * t396;
t169 = t425 * t396;
t168 = -mrSges(4,2) * t349 + mrSges(4,3) * t188;
t167 = mrSges(4,1) * t349 - mrSges(4,3) * t187;
t144 = -mrSges(5,1) * t443 + mrSges(5,2) * t213;
t138 = t220 * t502 + t347 - t516;
t94 = mrSges(6,1) * t181 - t540;
t93 = -mrSges(6,2) * t181 + t541;
t87 = -mrSges(5,2) * t343 + mrSges(5,3) * t98;
t86 = mrSges(5,1) * t343 - mrSges(5,3) * t97;
t57 = t126 * t395 - t393 * t83;
t56 = t140 * t395 - t393 * t74;
t55 = t119 * t395 - t393 * t79;
t54 = t124 * t395 - t393 * t76;
t52 = -qJD(5) * t139 - t120 * t397 + t121 * t502;
t51 = t121 * t503 + t120 * t402 + (t402 * t426 - t516) * qJD(5);
t26 = t121 * t550 + t42;
t23 = -mrSges(6,2) * t85 + mrSges(6,3) * t41;
t22 = mrSges(6,1) * t85 - mrSges(6,3) * t40;
t16 = -t27 * t393 + t395 * t61;
t6 = -qJD(5) * t35 - t26 * t397 + t434 * t402;
t5 = qJD(5) * t34 + t26 * t402 + t397 * t434;
t1 = [t121 * t453 + t34 * t22 + t35 * t23 + t48 * (-mrSges(6,1) * t52 + mrSges(6,2) * t51) + (Ifges(4,1) * t222 + Ifges(4,4) * t223) * t563 + (Ifges(5,1) * t120 + Ifges(5,4) * t121) * t565 + (Ifges(6,1) * t51 + Ifges(6,4) * t52 - Ifges(6,5) * t522) * t569 + (Ifges(6,5) * t139 + Ifges(6,6) * t138 + Ifges(6,3) * t202) * t573 + (Ifges(6,4) * t139 + Ifges(6,2) * t138 + Ifges(6,6) * t202) * t576 + (Ifges(6,1) * t139 + Ifges(6,4) * t138 + Ifges(6,5) * t202) * t577 + t139 * t578 + t138 * t579 + t202 * t580 + t57 * t15 + t52 * t63 / 0.2e1 + t51 * t64 / 0.2e1 + t16 * t72 + t5 * t93 + t6 * t94 + t110 * t86 + t111 * t87 + t120 * t133 / 0.2e1 + t8 * (-mrSges(6,1) * t138 + mrSges(6,2) * t139) + t42 * t189 + t43 * t190 + t200 * t167 + t201 * t168 + t3 * (-mrSges(6,2) * t202 + mrSges(6,3) * t138) + t4 * (mrSges(6,1) * t202 - mrSges(6,3) * t139) + t204 * t144 + t222 * t199 / 0.2e1 + t223 * t198 / 0.2e1 + t229 * (-mrSges(5,1) * t121 + mrSges(5,2) * t120) + t129 * t234 + t130 * t235 + t237 * (-mrSges(5,1) * t98 + mrSges(5,2) * t97) + t285 * (Ifges(4,4) * t222 + Ifges(4,2) * t223) / 0.2e1 + (-t120 * t91 + t121 * t92) * mrSges(5,3) + (-t182 * t222 + t183 * t223) * mrSges(4,3) + m(6) * (t16 * t48 + t24 * t6 + t25 * t5 + t3 * t35 + t34 * t4 + t57 * t8) + m(5) * (t110 * t21 + t111 * t20 + t149 * t237 + t204 * t229 + t42 * t92 + t43 * t91) + (mrSges(4,2) * t274 - mrSges(4,3) * t90 + Ifges(4,1) * t187 + Ifges(4,4) * t188 + Ifges(4,5) * t349) * t296 + (mrSges(5,2) * t149 - mrSges(5,3) * t21 + Ifges(5,1) * t97 + Ifges(5,4) * t98 + Ifges(5,5) * t343) * t221 + (-mrSges(4,1) * t274 + mrSges(4,3) * t89 + Ifges(4,4) * t187 + Ifges(4,2) * t188 + Ifges(4,6) * t349) * t295 + (-t323 * mrSges(3,1) - mrSges(6,1) * t423 - t322 * mrSges(3,2) - t278 * mrSges(4,2) - t251 * mrSges(5,2) - t191 * mrSges(6,2) + t401 * t410 - t406 * t409) * g(2) + (-mrSges(5,1) * t149 + mrSges(5,3) * t20 + Ifges(5,4) * t97 + Ifges(5,2) * t98 + Ifges(5,6) * t343) * t220 + Ifges(2,3) * qJDD(1) + t328 * (mrSges(3,1) * t364 - mrSges(3,3) * t312) + t359 * (Ifges(4,5) * t222 + Ifges(4,6) * t223) / 0.2e1 + t25 * (mrSges(6,2) * t522 + mrSges(6,3) * t52) + t344 * (Ifges(5,5) * t120 + Ifges(5,6) * t121) / 0.2e1 + t181 * (Ifges(6,5) * t51 + Ifges(6,6) * t52 - Ifges(6,3) * t522) / 0.2e1 + t127 * (Ifges(6,4) * t51 + Ifges(6,2) * t52 - Ifges(6,6) * t522) / 0.2e1 + t24 * (-mrSges(6,1) * t522 - mrSges(6,3) * t51) - t336 * (-mrSges(4,1) * t188 + mrSges(4,2) * t187) - t326 * (-mrSges(4,1) * t223 + mrSges(4,2) * t222) + t309 * t304 - t310 * t303 + t478 * (-mrSges(3,2) * t364 + mrSges(3,3) * t311) + m(3) * (t227 * t478 + t228 * t328 - t305 * t310 + t306 * t309) + t443 * (Ifges(5,4) * t120 + Ifges(5,2) * t121) / 0.2e1 + ((mrSges(3,1) * t311 - mrSges(3,2) * t312 + (m(3) * t559 - t324) * qJDD(1)) * pkin(1) + (-mrSges(6,2) * t506 + t411) * t549 + (-mrSges(6,1) * t507 + t411) * t548 + (t227 * mrSges(3,3) + Ifges(3,4) * t312 + Ifges(3,2) * t311 + Ifges(3,6) * t364) * t405 + (-t228 * mrSges(3,3) + Ifges(3,1) * t312 + Ifges(3,4) * t311 + Ifges(3,5) * t364) * t400 + ((-t305 * mrSges(3,3) + t527 / 0.2e1 + t268 / 0.2e1) * t405 + (-t306 * mrSges(3,3) - t526 / 0.2e1 - t267 / 0.2e1 + (-m(4) * t326 + t219) * pkin(2)) * t400 + (t405 * (Ifges(3,4) * t405 - Ifges(3,2) * t400) / 0.2e1 - t588) * t477) * qJD(2)) * t394 + ((-t548 - t549) * t378 * t539 + t582 + t583 + t584) * t396 + m(4) * (t129 * t183 + t130 * t182 + t200 * t90 + t201 * t89 - t274 * t336) + (-t321 * mrSges(3,1) - t192 * mrSges(6,1) + mrSges(3,2) * t416 - t276 * mrSges(4,2) - t250 * mrSges(5,2) - mrSges(6,2) * t422 + t401 * t409 + t406 * t410) * g(1) + t121 * t607; -t219 * t360 + t167 * t557 + (t399 * t89 + t404 * t90 + (-t182 * t399 + t183 * t404) * qJD(3)) * t581 + t267 * t449 / 0.2e1 + t588 * qJD(1) ^ 2 * t394 ^ 2 - t55 * t72 + t214 * t22 + t215 * t23 - t195 * t234 - t194 * t235 - t240 * t144 + t407 + t582 + (t399 * t168 + t234 * t474 - t235 * t475) * pkin(2) + t592 * t189 + t593 * t190 + (t20 * t317 + t21 * t316 - t229 * t240 + t592 * t92 + t593 * t91) * m(5) + t598 * t93 + (t214 * t4 + t215 * t3 + (-t231 * t48 - t313 * t8) * t393 - t48 * t55 + t598 * t25 + t599 * t24) * m(6) + t599 * t94 - t313 * t532 - t231 * t525 + t316 * t86 + t317 * t87 - m(4) * (t182 * t194 + t183 * t195 - t326 * t360) + (-m(5) * t394 * t418 - m(6) * (t236 * t504 + t225) - t504 * t581 + t324 + t587) * g(3) - (t366 * (Ifges(3,5) * t405 - Ifges(3,6) * t400) + (-Ifges(3,2) * t449 + t268 + t358) * t405) * t477 / 0.2e1 + (t440 + t303) * t306 + (t439 - t304) * t305 + (-m(6) * (t169 * t406 + t174 * t401) - m(5) * (t283 * t406 + t401 * t341) - mrSges(3,2) * t321 + t586 + t616 * t416) * g(2) + (-m(6) * (-t169 * t401 + t174 * t406) - m(5) * (-t401 * t283 + t341 * t406) + mrSges(3,2) * t323 + t616 * t322 + t585) * g(1); -t144 * t556 - t182 * t234 + t183 * t235 + t275 * t22 + t277 * t23 - t379 * t532 + t463 * t525 - t54 * t72 + t86 * t552 + t87 * t553 + t407 + t603 * t94 + t602 * t93 + (-t463 - t99) * t190 + (pkin(3) * t472 - t100) * t189 + (t275 * t4 + t277 * t3 + (-t379 * t8 + t463 * t48) * t393 - t48 * t54 + t602 * t25 + t603 * t24) * m(6) + ((t20 * t398 + t21 * t403 + (-t398 * t91 + t403 * t92) * qJD(4)) * pkin(3) - t100 * t92 - t229 * t556 - t91 * t99) * m(5) + (-m(6) * (-t238 * t504 + t225) - m(5) * t555 + t587) * g(3) + (-m(6) * (-t170 * t406 + t175 * t401) - m(5) * (-pkin(3) * t492 + t288 * t406) + t586) * g(2) + (-m(6) * (t170 * t401 + t175 * t406) - m(5) * (-pkin(3) * t486 - t401 * t288) + t585) * g(1); t213 * t453 + t132 * t565 - t56 * t72 + (t542 + t190) * t92 + (-Ifges(5,6) * t213 + t604) * t561 + (t543 - t189) * t91 + t589 * g(2) + t590 * g(1) + t591 * g(3) + (-Ifges(5,2) * t213 + t612) * t567 - t229 * (mrSges(5,1) * t213 + t606) + (-t537 + t605) * t566 + t600 * t93 + (g(3) * t394 * t424 - g(2) * (-t172 * t406 + t176 * t401) - g(1) * (t172 * t401 + t176 * t406) - t48 * t56 - pkin(4) * t546 + t3 * t329 + t327 * t4 + t600 * t25 + t601 * t24) * m(6) + t601 * t94 - pkin(4) * t532 + t327 * t22 + t329 * t23 + t615; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t48 * (mrSges(6,1) * t128 + mrSges(6,2) * t127) + (Ifges(6,1) * t127 - t536) * t570 + t63 * t569 + (Ifges(6,5) * t127 - Ifges(6,6) * t128) * t568 - g(1) * (t191 * mrSges(6,1) + (-t397 * t459 - t423) * mrSges(6,2)) - g(2) * ((-t402 * t458 - t422) * mrSges(6,1) + t192 * mrSges(6,2)) - g(3) * (-t396 * mrSges(6,2) * t507 + t347 * mrSges(6,1) + ((-t377 * t397 + t378 * t502) * mrSges(6,1) + (-t377 * t402 - t378 * t503) * mrSges(6,2)) * t394) + t12 + (-Ifges(6,2) * t128 + t125 + t64) * t571 + (t540 + t94) * t25 + (-t93 + t541) * t24;];
tau = t1;
