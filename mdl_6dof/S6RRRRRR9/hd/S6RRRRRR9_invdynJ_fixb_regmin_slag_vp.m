% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x38]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:30
% EndTime: 2019-03-10 05:40:26
% DurationCPUTime: 28.51s
% Computational Cost: add. (37998->934), mult. (106683->1329), div. (0->0), fcn. (92002->18), ass. (0->417)
t433 = cos(qJ(2));
t623 = cos(pkin(6));
t562 = pkin(1) * t623;
t410 = t433 * t562;
t400 = qJD(1) * t410;
t429 = sin(qJ(2));
t423 = sin(pkin(6));
t424 = cos(pkin(7));
t638 = pkin(10) * t424;
t507 = t423 * (-pkin(9) - t638);
t489 = t429 * t507;
t303 = qJD(1) * t489 + t400;
t409 = t429 * t562;
t444 = t433 * t507 - t409;
t304 = t444 * qJD(1);
t422 = sin(pkin(7));
t639 = pkin(10) * t422;
t487 = pkin(2) * t429 - t433 * t639;
t582 = qJD(1) * t423;
t342 = t487 * t582;
t428 = sin(qJ(3));
t608 = t424 * t428;
t612 = t422 * t428;
t643 = cos(qJ(3));
t404 = pkin(10) * t612;
t558 = t424 * t643;
t651 = pkin(2) * t558 - t404;
t690 = t651 * qJD(3) - t643 * t303 - t304 * t608 - t342 * t612;
t554 = t643 * t433;
t604 = t428 * t429;
t468 = -t424 * t604 + t554;
t338 = t468 * t423;
t326 = qJD(1) * t338;
t545 = qJD(3) * t643;
t514 = t422 * t545;
t689 = t514 - t326;
t231 = -t304 * t422 + t424 * t342;
t555 = t643 * t429;
t603 = t428 * t433;
t465 = t424 * t555 + t603;
t337 = t465 * t423;
t325 = qJD(1) * t337;
t688 = pkin(3) * t325 - pkin(11) * t326 + t231 - (pkin(3) * t428 - pkin(11) * t643) * t422 * qJD(3);
t553 = t429 * t582;
t519 = t422 * t553;
t687 = pkin(11) * t519 - t690;
t427 = sin(qJ(4));
t432 = cos(qJ(4));
t354 = -t432 * t424 + t427 * t612;
t593 = qJD(4) * t354 + t427 * t519 - t432 * t689;
t355 = t424 * t427 + t432 * t612;
t590 = qJD(4) * t355 + t427 * t689 + t432 * t519;
t580 = qJD(3) * t428;
t550 = t422 * t580;
t674 = t325 - t550;
t560 = t422 * t643;
t585 = pkin(2) * t608 + pkin(10) * t560;
t686 = t585 * qJD(3) - t428 * t303 + t304 * t558;
t534 = t623 * qJD(1);
t492 = t534 + qJD(2);
t474 = t492 * t422;
t552 = t433 * t582;
t493 = t643 * t552;
t518 = t428 * t553;
t661 = t424 * t493 + t643 * t474 - t518;
t685 = t661 - qJD(4);
t557 = t643 * t342;
t588 = -(-pkin(3) * t553 - t557) * t422 + t686;
t346 = pkin(11) * t424 + t585;
t347 = (-pkin(3) * t643 - pkin(11) * t428 - pkin(2)) * t422;
t577 = qJD(4) * t432;
t579 = qJD(4) * t427;
t684 = t346 * t579 - t347 * t577 + t427 * t688 + t687 * t432;
t434 = cos(qJ(1));
t536 = t434 * t623;
t642 = sin(qJ(1));
t356 = t429 * t642 - t433 * t536;
t357 = t429 * t536 + t433 * t642;
t609 = t423 * t434;
t245 = -t356 * t608 + t357 * t643 - t609 * t612;
t307 = -t356 * t422 + t424 * t609;
t187 = t245 * t432 - t307 * t427;
t244 = t356 * t558 + t357 * t428 + t560 * t609;
t421 = qJ(5) + qJ(6);
t416 = sin(t421);
t417 = cos(t421);
t683 = t187 * t416 - t244 * t417;
t682 = t187 * t417 + t244 * t416;
t426 = sin(qJ(5));
t431 = cos(qJ(5));
t681 = t187 * t426 - t244 * t431;
t680 = t187 * t431 + t244 * t426;
t679 = -t346 * t577 - t347 * t579 + t687 * t427 - t432 * t688;
t418 = t423 ^ 2;
t678 = 0.2e1 * t418;
t670 = -pkin(12) * t674 - t684;
t610 = t423 * t433;
t457 = pkin(9) * t610 + t409;
t265 = t457 * qJD(1) + (t424 * t552 + t474) * pkin(10);
t449 = pkin(2) * t623 + t489;
t270 = qJD(2) * pkin(2) + qJD(1) * t449 + t400;
t482 = -pkin(2) * t433 - t429 * t639 - pkin(1);
t334 = t482 * t423;
t321 = qJD(1) * t334;
t157 = t643 * t265 + t270 * t608 + t321 * t612;
t677 = t157 + t685 * (pkin(4) * t427 - pkin(12) * t432);
t676 = -t590 * pkin(4) - pkin(12) * t593 - t588;
t467 = t424 * t603 + t555;
t451 = t467 * t423;
t461 = t428 * t474;
t273 = qJD(1) * t451 + t461;
t192 = t426 * t432 * t661 - t431 * t273;
t663 = -t426 * t577 + t192;
t602 = t431 * t432;
t193 = t273 * t426 + t602 * t661;
t498 = t431 * t577 - t193;
t575 = qJD(5) * t426;
t673 = -t427 * t575 + t498;
t672 = pkin(1) * t678;
t450 = t422 * t552 - t424 * t492 - qJD(3);
t218 = t432 * t273 - t427 * t450;
t163 = t218 * t426 + t431 * t685;
t165 = t218 * t431 - t426 * t685;
t425 = sin(qJ(6));
t430 = cos(qJ(6));
t494 = t163 * t425 - t430 * t165;
t88 = t430 * t163 + t165 * t425;
t671 = t494 * t88;
t574 = qJD(5) * t431;
t664 = t427 * t574 - t663;
t466 = -t426 * t355 - t431 * t560;
t598 = qJD(5) * t466 - t674 * t426 - t431 * t593;
t522 = t426 * t560;
t597 = -qJD(5) * t522 + t355 * t574 - t426 * t593 + t674 * t431;
t320 = t432 * t450;
t216 = t273 * t427 + t320;
t647 = qJD(5) + qJD(6);
t662 = t216 + t647;
t660 = t494 ^ 2 - t88 ^ 2;
t524 = t623 * qJDD(1);
t403 = t524 + qJDD(2);
t569 = qJDD(1) * t433;
t541 = t423 * t569;
t570 = qJDD(1) * t429;
t542 = t423 * t570;
t176 = t403 * t612 + t541 * t608 + t643 * t542 + (-t424 * t518 + t493) * qJD(2) + t661 * qJD(3);
t571 = qJD(1) * qJD(2);
t544 = t429 * t571;
t513 = t423 * t544;
t278 = t403 * t424 + qJDD(3) + (t513 - t541) * t422;
t101 = -qJD(4) * t320 + t432 * t176 - t273 * t579 + t427 * t278;
t436 = qJD(2) * t465 + qJD(3) * t467;
t521 = t424 * t554;
t495 = t423 * t521;
t177 = qJD(3) * t461 - qJDD(1) * t495 - t403 * t560 + t423 * (qJD(1) * t436 + t428 * t570);
t173 = qJDD(4) + t177;
t50 = t431 * t101 + t426 * t173 - t218 * t575 - t574 * t685;
t51 = qJD(5) * t165 + t101 * t426 - t431 * t173;
t572 = qJD(6) * t430;
t573 = qJD(6) * t425;
t13 = -t163 * t572 - t165 * t573 - t425 * t51 + t430 * t50;
t214 = qJD(5) + t216;
t213 = qJD(6) + t214;
t659 = t213 * t88 + t13;
t506 = t623 * t642;
t358 = -t429 * t506 + t433 * t434;
t455 = t434 * t429 + t433 * t506;
t559 = t423 * t642;
t520 = t422 * t559;
t249 = t358 * t643 + (-t424 * t455 + t520) * t428;
t309 = t422 * t455 + t424 * t559;
t190 = t249 * t432 + t309 * t427;
t441 = t455 * t643;
t248 = t358 * t428 + t424 * t441 - t520 * t643;
t127 = t190 * t417 + t248 * t416;
t535 = t623 * t422;
t510 = t428 * t535;
t296 = t510 + t451;
t353 = t422 * t610 - t424 * t623;
t239 = t296 * t432 - t353 * t427;
t211 = -t270 * t422 + t424 * t321;
t131 = -pkin(3) * t661 - pkin(11) * t273 + t211;
t139 = -pkin(11) * t450 + t157;
t71 = t427 * t131 + t432 * t139;
t63 = -pkin(12) * t685 + t71;
t156 = -t428 * t265 + t270 * t558 + t321 * t560;
t138 = pkin(3) * t450 - t156;
t77 = t216 * pkin(4) - t218 * pkin(12) + t138;
t35 = t426 * t77 + t431 * t63;
t29 = -pkin(13) * t163 + t35;
t27 = t29 * t573;
t488 = t643 * t535;
t295 = t423 * t604 - t488 - t495;
t70 = t131 * t432 - t427 * t139;
t62 = pkin(4) * t685 - t70;
t46 = pkin(5) * t163 + t62;
t658 = t46 * t88 + g(1) * t127 + g(2) * t682 - g(3) * (-t239 * t417 - t295 * t416) + t27;
t126 = -t190 * t416 + t248 * t417;
t34 = -t426 * t63 + t431 * t77;
t28 = -pkin(13) * t165 + t34;
t21 = pkin(5) * t214 + t28;
t633 = t29 * t430;
t12 = t21 * t425 + t633;
t102 = qJD(4) * t218 + t427 * t176 - t432 * t278;
t100 = qJDD(5) + t102;
t516 = qJD(2) * t562;
t490 = qJD(1) * t516;
t511 = pkin(1) * t524;
t563 = pkin(9) * t541 + t429 * t511 + t433 * t490;
t456 = -pkin(9) * t513 + t563;
t649 = -t544 + t569;
t206 = (t423 * t424 * t649 + t403 * t422) * pkin(10) + t456;
t458 = -t429 * t490 + t433 * t511;
t543 = t433 * t571;
t479 = -t543 - t570;
t459 = t479 * pkin(9);
t212 = t403 * pkin(2) + (t479 * t638 + t459) * t423 + t458;
t469 = t487 * qJD(2);
t250 = (qJD(1) * t469 + qJDD(1) * t482) * t423;
t515 = t424 * t545;
t453 = -t643 * t206 - t212 * t608 - t250 * t612 + t265 * t580 - t270 * t515 - t321 * t514;
t58 = pkin(11) * t278 - t453;
t146 = -t212 * t422 + t424 * t250;
t64 = pkin(3) * t177 - pkin(11) * t176 + t146;
t18 = t131 * t577 - t139 * t579 + t427 * t64 + t432 * t58;
t16 = pkin(12) * t173 + t18;
t549 = t424 * t580;
t66 = -t428 * t206 + t212 * t558 + t250 * t560 - t265 * t545 - t270 * t549 - t321 * t550;
t59 = -pkin(3) * t278 - t66;
t24 = pkin(4) * t102 - pkin(12) * t101 + t59;
t7 = -qJD(5) * t35 - t426 * t16 + t431 * t24;
t4 = pkin(5) * t100 - pkin(13) * t50 + t7;
t484 = -t431 * t16 - t426 * t24 - t77 * t574 + t575 * t63;
t5 = -pkin(13) * t51 - t484;
t2 = -qJD(6) * t12 + t430 * t4 - t425 * t5;
t657 = t46 * t494 - g(1) * t126 + g(2) * t683 - g(3) * (-t239 * t416 + t295 * t417) + t2;
t445 = qJD(6) * t494 - t425 * t50 - t430 * t51;
t656 = -t213 * t494 + t445;
t655 = t245 * t427 + t307 * t432;
t630 = pkin(4) * t674 - t679;
t654 = t685 * t427;
t375 = t425 * t431 + t426 * t430;
t351 = t375 * t427;
t637 = pkin(11) * t426;
t653 = -t431 * t677 + t579 * t637;
t652 = t676 * t431;
t345 = t404 + (-pkin(2) * t643 - pkin(3)) * t424;
t240 = t354 * pkin(4) - t355 * pkin(12) + t345;
t589 = t432 * t346 + t427 * t347;
t243 = -pkin(12) * t560 + t589;
t594 = t426 * t240 + t431 * t243;
t650 = -t240 * t574 + t243 * t575 + t426 * t676 - t431 * t670;
t389 = -pkin(4) * t432 - pkin(12) * t427 - pkin(3);
t198 = pkin(3) * t273 - pkin(11) * t661;
t601 = t432 * t156 + t427 * t198;
t86 = pkin(12) * t273 + t601;
t648 = -t389 * t574 + t426 * t677 + t431 * t86;
t539 = qJD(6) * t21 + t5;
t646 = t425 * t4 + t430 * t539;
t435 = qJD(1) ^ 2;
t644 = pkin(12) + pkin(13);
t640 = pkin(5) * t427;
t302 = t410 + t449;
t227 = -t302 * t422 + t424 * t334;
t152 = pkin(3) * t295 - pkin(11) * t296 + t227;
t288 = (t424 * t610 + t535) * pkin(10) + t457;
t566 = t643 * t288 + t302 * t608 + t334 * t612;
t160 = -pkin(11) * t353 + t566;
t600 = t427 * t152 + t432 * t160;
t80 = pkin(12) * t295 + t600;
t447 = -t428 * t288 + t302 * t558 + t334 * t560;
t159 = t353 * pkin(3) - t447;
t238 = t296 * t427 + t353 * t432;
t95 = t238 * pkin(4) - t239 * pkin(12) + t159;
t635 = t426 * t95 + t431 * t80;
t634 = pkin(11) * qJD(4);
t632 = t426 * t50;
t631 = pkin(5) * t597 + t630;
t150 = t427 * t156;
t85 = -pkin(4) * t273 - t198 * t432 + t150;
t629 = pkin(5) * t664 + pkin(11) * t577 - t85;
t143 = pkin(4) * t218 + pkin(12) * t216;
t628 = t426 * t143 + t431 * t70;
t313 = t431 * t355 - t522;
t221 = t313 * t425 - t430 * t466;
t625 = -qJD(6) * t221 - t425 * t597 + t430 * t598;
t222 = t313 * t430 + t425 * t466;
t624 = qJD(6) * t222 + t425 * t598 + t430 * t597;
t622 = t100 * t431;
t621 = t163 * t214;
t620 = t165 * t214;
t619 = t216 * t685;
t618 = t216 * t426;
t617 = t218 * t685;
t616 = t416 * t432;
t615 = t417 * t432;
t614 = t418 * t435;
t613 = t422 * t427;
t611 = t423 * t429;
t607 = t426 * t100;
t606 = t426 * t427;
t605 = t427 * t431;
t374 = t425 * t426 - t430 * t431;
t596 = t192 * t425 - t193 * t430 - t351 * t647 - t374 * t577;
t595 = -t573 * t606 + (t605 * t647 - t663) * t430 + t673 * t425;
t592 = t662 * t374;
t591 = t662 * t375;
t411 = pkin(11) * t602;
t584 = t426 * t389 + t411;
t419 = t429 ^ 2;
t583 = -t433 ^ 2 + t419;
t581 = qJD(2) * t423;
t578 = qJD(4) * t431;
t576 = qJD(5) * t214;
t568 = t433 * t614;
t567 = t422 * t611;
t561 = qJD(5) * t644;
t343 = t423 * t469;
t556 = t643 * t343;
t551 = t429 * t581;
t537 = -t426 * t80 + t431 * t95;
t532 = t152 * t432 - t427 * t160;
t530 = t431 * t240 - t243 * t426;
t306 = t444 * qJD(2);
t232 = -t306 * t422 + t424 * t343;
t529 = -t427 * t346 + t347 * t432;
t526 = t685 * t432;
t525 = t214 * t431;
t19 = -t131 * t579 - t139 * t577 - t427 * t58 + t432 * t64;
t517 = t422 * t551;
t512 = -t71 + (t575 + t618) * pkin(5);
t120 = pkin(5) * t354 - pkin(13) * t313 + t530;
t509 = pkin(13) * t597 - qJD(6) * t120 + t650;
t140 = pkin(13) * t466 + t594;
t508 = -pkin(5) * t590 + pkin(13) * t598 + qJD(5) * t594 + qJD(6) * t140 + t426 * t670 + t652;
t505 = t423 * t435 * t623;
t316 = -pkin(13) * t606 + t584;
t503 = -pkin(13) * t193 + qJD(6) * t316 + t661 * t640 - t426 * t86 - (-pkin(13) * t602 + t640) * qJD(4) - (-t411 + (pkin(13) * t427 - t389) * t426) * qJD(5) - t653;
t373 = t431 * t389;
t294 = -pkin(13) * t605 + t373 + (-pkin(5) - t637) * t432;
t502 = -qJD(6) * t294 - (-t427 * t578 - t432 * t575) * pkin(11) + t648 + t664 * pkin(13);
t142 = t431 * t143;
t395 = t644 * t431;
t501 = pkin(5) * t218 + qJD(6) * t395 - t426 * t70 + t142 + (pkin(13) * t216 + t561) * t431;
t394 = t644 * t426;
t500 = pkin(13) * t618 + qJD(6) * t394 + t426 * t561 + t628;
t242 = pkin(4) * t560 - t529;
t185 = t239 * t431 + t295 * t426;
t33 = pkin(5) * t238 - pkin(13) * t185 + t537;
t184 = t239 * t426 - t295 * t431;
t36 = -pkin(13) * t184 + t635;
t497 = t33 * t430 - t36 * t425;
t496 = t33 * t425 + t36 * t430;
t113 = t430 * t184 + t185 * t425;
t114 = -t184 * t425 + t185 * t430;
t491 = 0.2e1 * t534 + qJD(2);
t401 = t433 * t516;
t305 = qJD(2) * t489 + t401;
t452 = -t288 * t580 + t302 * t515 + t643 * t305 + t306 * t608 + t334 * t514 + t343 * t612;
t105 = pkin(11) * t517 + t452;
t229 = qJD(3) * t510 + t423 * t436;
t230 = qJD(3) * t488 + ((t521 - t604) * qJD(3) + t468 * qJD(2)) * t423;
t116 = pkin(3) * t229 - pkin(11) * t230 + t232;
t486 = -t427 * t105 + t116 * t432 - t152 * t579 - t160 * t577;
t79 = -pkin(4) * t295 - t532;
t17 = -pkin(4) * t173 - t19;
t478 = t432 * t105 + t427 * t116 + t152 * t577 - t160 * t579;
t31 = pkin(12) * t229 + t478;
t454 = -t288 * t545 - t302 * t549 - t428 * t305 + t306 * t558 - t334 * t550;
t106 = (-pkin(3) * t551 - t556) * t422 - t454;
t144 = qJD(4) * t239 + t230 * t427 - t432 * t517;
t145 = -qJD(4) * t238 + t230 * t432 + t427 * t517;
t45 = t144 * pkin(4) - t145 * pkin(12) + t106;
t483 = t431 * t31 + t426 * t45 + t95 * t574 - t575 * t80;
t481 = -pkin(12) * t100 + t214 * t62;
t480 = -pkin(11) * t173 - t138 * t685;
t475 = t423 * (t524 + t403);
t189 = -t249 * t427 + t309 * t432;
t473 = g(1) * t189 - g(2) * t655 - g(3) * t238;
t472 = g(1) * t248 + g(2) * t244 + g(3) * t295;
t471 = g(1) * t249 + g(2) * t245 + g(3) * t296;
t262 = -t356 * t643 - t357 * t608;
t264 = -t358 * t608 - t441;
t470 = g(1) * t264 + g(2) * t262 + g(3) * t338;
t462 = t472 - t59;
t460 = g(1) * t358 + g(2) * t357 + g(3) * t611;
t32 = -pkin(4) * t229 - t486;
t448 = pkin(11) * t576 - t472;
t446 = -qJD(5) * t635 - t31 * t426 + t431 * t45;
t443 = pkin(12) * t576 + t17 + t473;
t442 = qJD(3) * t450;
t439 = t450 * t567;
t438 = t492 * t457;
t415 = -pkin(5) * t431 - pkin(4);
t382 = (pkin(5) * t426 + pkin(11)) * t427;
t352 = t374 * t427;
t274 = t338 * t432 + t427 * t567;
t263 = t358 * t558 - t428 * t455;
t261 = -t356 * t428 + t357 * t558;
t220 = t264 * t432 + t358 * t613;
t219 = t262 * t432 + t357 * t613;
t194 = -pkin(5) * t466 + t242;
t135 = t190 * t431 + t248 * t426;
t134 = -t190 * t426 + t248 * t431;
t97 = qJDD(6) + t100;
t74 = -qJD(5) * t184 + t145 * t431 + t229 * t426;
t73 = qJD(5) * t185 + t145 * t426 - t229 * t431;
t57 = pkin(5) * t184 + t79;
t26 = qJD(6) * t114 + t425 * t74 + t430 * t73;
t25 = -qJD(6) * t113 - t425 * t73 + t430 * t74;
t20 = pkin(5) * t73 + t32;
t11 = t21 * t430 - t29 * t425;
t10 = pkin(5) * t51 + t17;
t9 = -pkin(13) * t73 + t483;
t8 = pkin(5) * t144 - pkin(13) * t74 + t446;
t1 = -t27 + t646;
t3 = [qJDD(1), g(1) * t642 - g(2) * t434, g(1) * t434 + g(2) * t642 (qJDD(1) * t419 + 0.2e1 * t429 * t543) * t418 (t429 * t569 - t571 * t583) * t678, t433 * t491 * t581 + t429 * t475, t433 * t475 - t491 * t551, t403 * t623, -qJD(2) * t438 + (-pkin(9) * t611 + t410) * t403 + (t423 * t459 + t458) * t623 + g(1) * t357 - g(2) * t358 + t649 * t672 -(-pkin(9) * t551 + t401) * t492 - t457 * t403 - t456 * t623 - g(1) * t356 + g(2) * t455 + t479 * t672, t176 * t296 + t230 * t273, -t176 * t295 - t177 * t296 - t229 * t273 + t230 * t661, -t176 * t353 - t230 * t450 + t273 * t517 + t296 * t278, t177 * t353 + t229 * t450 - t295 * t278 + t517 * t661, -qJD(2) * t439 - t278 * t353 -(t422 * t556 + t454) * t450 + t447 * t278 - t66 * t353 + t156 * t517 - t232 * t661 + t227 * t177 + t146 * t295 + t211 * t229 + g(1) * t245 - g(2) * t249, -g(1) * t244 + g(2) * t248 + t146 * t296 - t157 * t517 + t227 * t176 + t211 * t230 + t232 * t273 - t278 * t566 - t353 * t453 + t450 * t452, t101 * t239 + t145 * t218, -t101 * t238 - t102 * t239 - t144 * t218 - t145 * t216, t101 * t295 - t145 * t685 + t173 * t239 + t218 * t229, -t102 * t295 + t144 * t685 - t173 * t238 - t216 * t229, t173 * t295 - t229 * t685, g(1) * t187 - g(2) * t190 + t159 * t102 + t106 * t216 + t138 * t144 + t173 * t532 + t19 * t295 + t70 * t229 + t59 * t238 - t486 * t685, -g(1) * t655 - g(2) * t189 + t159 * t101 + t106 * t218 + t138 * t145 - t600 * t173 - t18 * t295 - t71 * t229 + t59 * t239 + t478 * t685, t165 * t74 + t185 * t50, -t163 * t74 - t165 * t73 - t184 * t50 - t185 * t51, t100 * t185 + t144 * t165 + t214 * t74 + t238 * t50, -t100 * t184 - t144 * t163 - t214 * t73 - t238 * t51, t100 * t238 + t144 * t214, g(1) * t680 - g(2) * t135 + t537 * t100 + t34 * t144 + t32 * t163 + t17 * t184 + t446 * t214 + t7 * t238 + t79 * t51 + t62 * t73, -g(1) * t681 - g(2) * t134 - t635 * t100 - t35 * t144 + t32 * t165 + t17 * t185 - t483 * t214 + t484 * t238 + t79 * t50 + t62 * t74, t114 * t13 - t25 * t494, -t113 * t13 + t114 * t445 - t25 * t88 + t26 * t494, t114 * t97 + t13 * t238 - t144 * t494 + t213 * t25, -t113 * t97 - t144 * t88 - t213 * t26 + t238 * t445, t144 * t213 + t238 * t97 (-qJD(6) * t496 - t425 * t9 + t430 * t8) * t213 + t497 * t97 + t2 * t238 + t11 * t144 + t20 * t88 - t57 * t445 + t10 * t113 + t46 * t26 + g(1) * t682 - g(2) * t127 -(qJD(6) * t497 + t425 * t8 + t430 * t9) * t213 - t496 * t97 - t1 * t238 - t12 * t144 - t20 * t494 + t57 * t13 + t10 * t114 + t46 * t25 - g(1) * t683 - g(2) * t126; 0, 0, 0, -t429 * t568, t583 * t614, -t433 * t505 + t542, t429 * t505 + t541, t403, pkin(1) * t429 * t614 + g(1) * t455 + g(2) * t356 - g(3) * t610 + qJD(1) * t438 + t458 + (-t423 * t543 - t542) * pkin(9) (-pkin(9) * t553 + t400) * t534 + pkin(1) * t568 + t400 * qJD(2) + t460 - t563, t176 * t612 + t273 * t689, -t326 * t661 + t273 * t325 + (t643 * t176 - t177 * t428 + (-t273 * t428 + t643 * t661) * qJD(3)) * t422, t176 * t424 + t326 * t450 + (-t273 * t553 + t428 * t278 - t442 * t643) * t422, -t177 * t424 - t325 * t450 + (t278 * t643 + t428 * t442 - t553 * t661) * t422, qJD(1) * t439 + t278 * t424, t651 * t278 + t66 * t424 - t422 * pkin(2) * t177 - t146 * t560 - t156 * t519 + t231 * t661 - t674 * t211 - t470 + (t422 * t557 + t686) * t450, -t585 * t278 + t453 * t424 - t231 * t273 - t211 * t326 + g(1) * t263 + g(2) * t261 + g(3) * t337 + (-pkin(2) * t176 + t146 * t428 + t157 * t553 + t211 * t545) * t422 + t690 * t450, t101 * t355 - t218 * t593, -t101 * t354 - t102 * t355 + t216 * t593 - t218 * t590, -t101 * t560 + t173 * t355 - t218 * t674 + t593 * t685, t102 * t560 - t354 * t173 + t216 * t674 + t590 * t685, -t173 * t560 + t674 * t685, t529 * t173 + t345 * t102 + t59 * t354 - t70 * t325 - g(1) * t220 - g(2) * t219 - g(3) * t274 + (-t19 * t643 + t580 * t70) * t422 - t679 * t685 + t588 * t216 + t590 * t138, -t589 * t173 + t345 * t101 + t59 * t355 + t71 * t325 + t470 * t427 + (t18 * t643 - t432 * t460 - t580 * t71) * t422 - t684 * t685 + t588 * t218 - t593 * t138, t165 * t598 + t313 * t50, -t163 * t598 - t165 * t597 - t313 * t51 + t466 * t50, t100 * t313 + t165 * t590 + t214 * t598 + t354 * t50, t100 * t466 - t163 * t590 - t214 * t597 - t354 * t51, t100 * t354 + t214 * t590, t530 * t100 + t7 * t354 + t242 * t51 - t17 * t466 - g(1) * (t220 * t431 + t263 * t426) - g(2) * (t219 * t431 + t261 * t426) - g(3) * (t274 * t431 + t337 * t426) + t597 * t62 + t590 * t34 + (-t243 * t574 + (-qJD(5) * t240 - t670) * t426 - t652) * t214 + t630 * t163, -t594 * t100 + t484 * t354 + t242 * t50 + t17 * t313 - g(1) * (-t220 * t426 + t263 * t431) - g(2) * (-t219 * t426 + t261 * t431) - g(3) * (-t274 * t426 + t337 * t431) + t598 * t62 - t590 * t35 + t650 * t214 + t630 * t165, t13 * t222 - t494 * t625, -t13 * t221 + t222 * t445 + t494 * t624 - t625 * t88, t13 * t354 + t213 * t625 + t222 * t97 - t494 * t590, -t213 * t624 - t221 * t97 + t354 * t445 - t590 * t88, t213 * t590 + t354 * t97 (t120 * t430 - t140 * t425) * t97 + t2 * t354 - t194 * t445 + t10 * t221 - g(1) * (t220 * t417 + t263 * t416) - g(2) * (t219 * t417 + t261 * t416) - g(3) * (t274 * t417 + t337 * t416) + t631 * t88 + t624 * t46 + (t425 * t509 - t430 * t508) * t213 + t590 * t11 -(t120 * t425 + t140 * t430) * t97 - t1 * t354 + t194 * t13 + t10 * t222 - g(1) * (-t220 * t416 + t263 * t417) - g(2) * (-t219 * t416 + t261 * t417) - g(3) * (-t274 * t416 + t337 * t417) - t631 * t494 + t625 * t46 + (t425 * t508 + t430 * t509) * t213 - t590 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273 * t661, t273 ^ 2 - t661 ^ 2, t450 * t661 + t176, -t273 * t450 - t177, t278, -t157 * t450 - t211 * t273 + t472 + t66, -t156 * t450 - t211 * t661 + t453 + t471, t101 * t427 - t218 * t526 (t101 + t619) * t432 + (-t102 + t617) * t427, t173 * t427 - t218 * t273 + t526 * t685, t173 * t432 + t216 * t273 - t654 * t685, t685 * t273, -pkin(3) * t102 - t150 * t685 - t157 * t216 - t70 * t273 + t480 * t427 + (-(-t198 - t634) * t685 + t462) * t432, -pkin(3) * t101 - t601 * t685 + t71 * t273 - t157 * t218 + t480 * t432 + (-t634 * t685 - t462) * t427, t165 * t673 + t50 * t605, t163 * t193 + t165 * t192 + (-t163 * t431 - t165 * t426) * t577 + (-t632 - t431 * t51 + (t163 * t426 - t165 * t431) * qJD(5)) * t427, -t432 * t50 + t498 * t214 + (-t165 * t685 - t214 * t575 + t622) * t427, t432 * t51 + t663 * t214 + (t163 * t685 - t214 * t574 - t607) * t427, -t100 * t432 - t214 * t654, t373 * t100 - t85 * t163 - t62 * t192 + t653 * t214 + ((-qJD(5) * t389 + t86) * t214 - t471) * t426 + (t62 * t426 * qJD(4) - t7 + (qJD(4) * t163 - t607) * pkin(11) - t448 * t431) * t432 + (pkin(11) * t51 + t17 * t426 - t34 * t685 + t574 * t62) * t427, -t584 * t100 - t85 * t165 - t62 * t193 + t648 * t214 - t471 * t431 + (-t484 + (pkin(11) * t165 + t431 * t62) * qJD(4) + t448 * t426) * t432 + (-t62 * t575 + t17 * t431 + t685 * t35 + (t214 * t578 + t50) * pkin(11)) * t427, -t13 * t352 - t494 * t596, -t13 * t351 - t352 * t445 + t494 * t595 - t596 * t88, -t13 * t432 + t213 * t596 - t352 * t97 + t494 * t654, -t213 * t595 - t351 * t97 - t432 * t445 + t654 * t88, -t213 * t654 - t432 * t97 (t294 * t430 - t316 * t425) * t97 - t2 * t432 - t382 * t445 + t10 * t351 - g(1) * (-t248 * t615 + t249 * t416) - g(2) * (-t244 * t615 + t245 * t416) - g(3) * (-t295 * t615 + t296 * t416) + t629 * t88 + t595 * t46 + (t425 * t502 - t430 * t503) * t213 - t11 * t654 -(t294 * t425 + t316 * t430) * t97 + t1 * t432 + t382 * t13 - t10 * t352 - g(1) * (t248 * t616 + t249 * t417) - g(2) * (t244 * t616 + t245 * t417) - g(3) * (t295 * t616 + t296 * t417) - t629 * t494 + t596 * t46 + (t425 * t503 + t430 * t502) * t213 + t12 * t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216 * t218, -t216 ^ 2 + t218 ^ 2, t101 - t619, -t102 - t617, t173, -t138 * t218 - t685 * t71 + t19 - t473, g(1) * t190 + g(2) * t187 + g(3) * t239 + t138 * t216 - t685 * t70 - t18, t165 * t525 + t632 (t50 - t621) * t431 + (-t51 - t620) * t426, -t165 * t218 + t214 * t525 + t607, -t214 ^ 2 * t426 + t163 * t218 + t622, -t214 * t218, -pkin(4) * t51 - t142 * t214 - t71 * t163 - t34 * t218 + (t70 * t214 + t481) * t426 - t443 * t431, -pkin(4) * t50 - t71 * t165 + t214 * t628 + t35 * t218 + t426 * t443 + t431 * t481, t13 * t375 + t494 * t592, -t13 * t374 + t375 * t445 + t494 * t591 + t592 * t88, -t213 * t592 + t218 * t494 + t375 * t97, -t213 * t591 + t218 * t88 - t374 * t97, -t213 * t218 (-t394 * t430 - t395 * t425) * t97 - t415 * t445 + t10 * t374 - t11 * t218 + t512 * t88 + t591 * t46 + (t425 * t500 - t430 * t501) * t213 - t473 * t417 -(-t394 * t425 + t395 * t430) * t97 + t415 * t13 + t10 * t375 + t12 * t218 - t512 * t494 - t592 * t46 + (t425 * t501 + t430 * t500) * t213 + t473 * t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165 * t163, -t163 ^ 2 + t165 ^ 2, t50 + t621, -t51 + t620, t100, -g(1) * t134 + g(2) * t681 + g(3) * t184 - t62 * t165 + t35 * t214 + t7, g(1) * t135 + g(2) * t680 + g(3) * t185 + t62 * t163 + t34 * t214 + t484, -t671, t660, t659, t656, t97 -(-t28 * t425 - t633) * t213 + (-t165 * t88 - t213 * t573 + t430 * t97) * pkin(5) + t657 (-t213 * t29 - t4) * t425 + (t213 * t28 - t539) * t430 + (t165 * t494 - t213 * t572 - t425 * t97) * pkin(5) + t658; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t671, t660, t659, t656, t97, t12 * t213 + t657, t11 * t213 - t646 + t658;];
tau_reg  = t3;
