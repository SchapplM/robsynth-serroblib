% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:38
% EndTime: 2019-03-09 13:53:49
% DurationCPUTime: 7.45s
% Computational Cost: add. (5640->417), mult. (13144->564), div. (0->0), fcn. (9496->8), ass. (0->197)
t513 = sin(qJ(4));
t517 = cos(qJ(4));
t518 = cos(qJ(2));
t585 = qJD(1) * t518;
t514 = sin(qJ(2));
t586 = qJD(1) * t514;
t446 = -t513 * t586 - t517 * t585;
t448 = -t513 * t585 + t517 * t586;
t512 = sin(qJ(5));
t516 = cos(qJ(5));
t539 = t446 * t512 + t516 * t448;
t511 = sin(qJ(6));
t578 = qJD(6) * t511;
t515 = cos(qJ(6));
t461 = t513 * t514 + t517 * t518;
t529 = t461 * qJD(4);
t574 = qJD(1) * qJD(2);
t566 = t518 * t574;
t567 = t514 * t574;
t592 = t513 * t567 + t517 * t566;
t413 = -qJD(1) * t529 + t592;
t581 = qJD(4) * t517;
t582 = qJD(4) * t513;
t583 = qJD(2) * t518;
t662 = t513 * t583 + t514 * t581 - t518 * t582;
t414 = qJD(1) * t662 - t517 * t567;
t579 = qJD(5) * t516;
t580 = qJD(5) * t512;
t530 = -t516 * t413 + t512 * t414 - t446 * t579 + t448 * t580;
t577 = qJD(6) * t515;
t505 = qJD(2) - qJD(4);
t659 = -qJD(5) + t505;
t598 = -t515 * t530 - t577 * t659;
t346 = -t539 * t578 + t598;
t344 = t346 * t511;
t345 = t346 * t515;
t395 = -t511 * t659 + t515 * t539;
t616 = t530 * t511;
t347 = qJD(6) * t395 - t616;
t366 = qJD(5) * t539 + t413 * t512 + t516 * t414;
t362 = t515 * t366;
t609 = t539 * t511;
t393 = t515 * t659 + t609;
t360 = t511 * t366;
t632 = -t516 * t446 + t448 * t512;
t642 = qJD(6) + t632;
t599 = t577 * t642 + t360;
t649 = t632 * t515;
t669 = t642 * t511;
t674 = t577 + t649;
t677 = -(t511 * t347 + t674 * t393 + t395 * t669 - t345) * MDP(30) - (-t393 * t539 + t642 * t669 - t362) * MDP(32) - (-t539 ^ 2 + t632 ^ 2) * MDP(23) + (t674 * t395 + t344) * MDP(29) + (-t395 * t539 + t642 * t649 + t599) * MDP(31) - (t632 * t659 + t530) * MDP(24) - t366 * MDP(25) + (t632 * MDP(22) - t659 * MDP(25) - MDP(33) * t642) * t539;
t676 = -(t446 ^ 2 - t448 ^ 2) * MDP(16) - t448 * t446 * MDP(15) - (t448 * t505 + t414) * MDP(18) + t677;
t496 = pkin(7) * t586;
t671 = -pkin(8) * t586 + qJD(3) + t496;
t497 = pkin(7) * t585;
t468 = -pkin(8) * t585 + t497;
t519 = -pkin(2) - pkin(3);
t545 = -qJ(3) * t513 + t517 * t519;
t664 = qJD(4) * t545 - t513 * t468 + t671 * t517;
t471 = qJ(3) * t517 + t513 * t519;
t663 = qJD(4) * t471 + t517 * t468 + t671 * t513;
t379 = pkin(5) * t539 + pkin(10) * t632;
t571 = t519 * qJD(2);
t434 = t571 + t671;
t507 = qJD(2) * qJ(3);
t449 = t468 + t507;
t561 = t517 * t434 - t449 * t513;
t619 = pkin(9) * t448;
t388 = t561 - t619;
t385 = -pkin(4) * t505 + t388;
t540 = t434 * t513 + t449 * t517;
t620 = pkin(9) * t446;
t389 = t540 + t620;
t614 = t389 * t512;
t352 = t385 * t516 - t614;
t350 = pkin(5) * t659 - t352;
t651 = t350 * t632;
t647 = t620 + t663;
t646 = t619 - t664;
t450 = -qJD(1) * pkin(1) - pkin(2) * t585 - qJ(3) * t586;
t432 = pkin(3) * t585 - t450;
t584 = qJD(2) * t514;
t622 = pkin(7) - pkin(8);
t467 = t622 * t584;
t506 = qJD(2) * qJD(3);
t439 = -qJD(1) * t467 + t506;
t489 = pkin(7) * t566;
t457 = -pkin(8) * t566 + t489;
t525 = -t540 * qJD(4) - t513 * t439 + t517 * t457;
t645 = t432 * t448 - t525;
t532 = -t434 * t581 - t517 * t439 + t449 * t582 - t513 * t457;
t363 = -pkin(9) * t414 - t532;
t364 = -pkin(9) * t413 + t525;
t335 = t512 * t363 - t516 * t364 + t385 * t580 + t389 * t579;
t412 = -pkin(4) * t446 + t432;
t630 = -t539 * t412 - t335;
t613 = t389 * t516;
t353 = t385 * t512 + t613;
t351 = -pkin(10) * t659 + t353;
t367 = pkin(5) * t632 - pkin(10) * t539 + t412;
t342 = t351 * t515 + t367 * t511;
t628 = t335 * t511 + t342 * t539 + t350 * t577;
t549 = -t516 * t363 - t512 * t364 - t385 * t579 + t389 * t580;
t640 = t632 * t412 + t549;
t542 = t351 * t511 - t367 * t515;
t629 = -t335 * t515 + t350 * t578 + t542 * t539;
t636 = -0.2e1 * t574;
t508 = t514 ^ 2;
t634 = MDP(5) * (-t518 ^ 2 + t508);
t423 = -t517 * t584 + t662;
t477 = t622 * t518;
t469 = qJD(2) * t477;
t476 = t622 * t514;
t531 = -t517 * t467 + t513 * t469 + t476 * t581 - t477 * t582;
t375 = -pkin(9) * t423 + t531;
t424 = qJD(2) * t461 - t529;
t536 = -t476 * t513 - t477 * t517;
t524 = qJD(4) * t536 + t467 * t513 + t517 * t469;
t376 = -pkin(9) * t424 + t524;
t462 = -t513 * t518 + t514 * t517;
t399 = -pkin(9) * t462 + t476 * t517 - t477 * t513;
t400 = -pkin(9) * t461 - t536;
t541 = t399 * t516 - t400 * t512;
t338 = qJD(5) * t541 + t375 * t516 + t376 * t512;
t417 = t516 * t461 + t462 * t512;
t371 = -qJD(5) * t417 - t423 * t512 + t424 * t516;
t418 = -t461 * t512 + t462 * t516;
t473 = -t518 * pkin(2) - t514 * qJ(3) - pkin(1);
t458 = t518 * pkin(3) - t473;
t425 = pkin(4) * t461 + t458;
t374 = pkin(5) * t417 - pkin(10) * t418 + t425;
t378 = t399 * t512 + t400 * t516;
t624 = t335 * t418 + t350 * t371 - t378 * t366 - (qJD(6) * t374 + t338) * t642 - (qJD(6) * t367 - t549) * t417;
t621 = pkin(4) * t448;
t618 = qJD(2) * pkin(2);
t617 = t350 * t418;
t615 = t374 * t366;
t607 = t446 * t505;
t520 = qJD(2) ^ 2;
t604 = t514 * t520;
t603 = t518 * t520;
t521 = qJD(1) ^ 2;
t602 = t518 * t521;
t465 = -pkin(4) + t545;
t538 = t465 * t516 - t471 * t512;
t597 = -qJD(5) * t538 + t512 * t647 + t516 * t646;
t537 = t465 * t512 + t471 * t516;
t596 = qJD(5) * t537 - t512 * t646 + t516 * t647;
t460 = t512 * t513 - t516 * t517;
t595 = t659 * t460;
t463 = t512 * t517 + t513 * t516;
t594 = t659 * t463;
t501 = t514 * qJD(3);
t590 = qJ(3) * t566 + qJD(1) * t501;
t589 = qJ(3) * t583 + t501;
t573 = t514 * t602;
t565 = pkin(1) * t636;
t564 = qJD(3) - t618;
t555 = qJD(1) * t473 + t450;
t492 = qJ(3) * t585;
t442 = t519 * t586 + t492;
t415 = t442 - t621;
t420 = -pkin(10) + t537;
t552 = qJD(6) * t420 - t379 + t415;
t494 = pkin(4) * t512 + pkin(10);
t551 = qJD(6) * t494 + t379 + t621;
t550 = t505 ^ 2;
t548 = t514 * t571;
t354 = t388 * t512 + t613;
t547 = pkin(4) * t580 - t354;
t355 = t388 * t516 - t614;
t546 = -pkin(4) * t579 + t355;
t544 = -t366 * t420 - t651;
t543 = -t366 * t494 + t651;
t535 = qJD(6) * t463 + t586;
t431 = pkin(2) * t567 - t590;
t443 = pkin(2) * t584 - t589;
t534 = -pkin(7) * t520 - qJD(1) * t443 - t431;
t533 = t371 * t515 - t418 * t578;
t428 = t548 + t589;
t426 = qJD(1) * t548 + t590;
t392 = pkin(4) * t423 + t428;
t526 = -t432 * t446 + t532;
t384 = pkin(4) * t414 + t426;
t470 = -pkin(7) * t567 + t506;
t472 = t496 + t564;
t474 = t497 + t507;
t522 = t470 * t518 + (t472 * t518 + (-t474 + t497) * t514) * qJD(2);
t495 = -pkin(4) * t516 - pkin(5);
t464 = pkin(2) * t586 - t492;
t419 = pkin(5) - t538;
t372 = qJD(5) * t418 + t516 * t423 + t424 * t512;
t340 = pkin(5) * t372 - pkin(10) * t371 + t392;
t339 = qJD(5) * t378 + t375 * t512 - t376 * t516;
t337 = pkin(5) * t366 + pkin(10) * t530 + t384;
t336 = t515 * t337;
t1 = [(-t424 * MDP(17) + t423 * MDP(18) - MDP(20) * t524 + MDP(21) * t531) * t505 + (t339 * t395 - t342 * t372 - t541 * t346 + (-(-qJD(6) * t378 + t340) * t642 - t615 - (-qJD(6) * t351 + t337) * t417 - qJD(6) * t617) * t511 + t624 * t515) * MDP(35) + (t336 * t417 + t339 * t393 - t542 * t372 - t541 * t347 + (t340 * t642 + t615 + (-t351 * t417 - t378 * t642 + t617) * qJD(6)) * t515 + t624 * t511) * MDP(34) + (-t418 * t360 - t347 * t417 - t372 * t393 + (-t371 * t511 - t418 * t577) * t642) * MDP(32) + (t346 * t417 + t362 * t418 + t372 * t395 + t533 * t642) * MDP(31) + (t366 * t417 + t372 * t642) * MDP(33) + (-MDP(24) * t371 + MDP(25) * t372 + MDP(27) * t339 + MDP(28) * t338) * t659 + (t366 * t425 + t372 * t412 + t384 * t417 + t392 * t632) * MDP(27) + (-t366 * t418 - t371 * t632 - t372 * t539 + t417 * t530) * MDP(23) + (t371 * t412 + t384 * t418 + t392 * t539 - t425 * t530) * MDP(28) + (t371 * t539 - t418 * t530) * MDP(22) + (t345 * t418 + t395 * t533) * MDP(29) + ((-t393 * t515 - t395 * t511) * t371 + (-t344 - t347 * t515 + (t393 * t511 - t395 * t515) * qJD(6)) * t418) * MDP(30) - MDP(7) * t604 + (pkin(7) * t604 + t518 * t565) * MDP(10) + (-pkin(7) * t603 + t514 * t565) * MDP(9) + (t514 * t534 - t555 * t583) * MDP(13) + (t518 * t534 + t555 * t584) * MDP(11) + (-t413 * t461 - t414 * t462 - t423 * t448 + t424 * t446) * MDP(16) + (t413 * t462 + t424 * t448) * MDP(15) + MDP(6) * t603 + (t458 * t414 + t432 * t423 + t426 * t461 - t428 * t446) * MDP(20) + (t458 * t413 + t432 * t424 + t426 * t462 + t428 * t448) * MDP(21) + t522 * MDP(12) + (pkin(7) * t522 + t431 * t473 + t443 * t450) * MDP(14) + 0.2e1 * t514 * MDP(4) * t566 + t634 * t636; (t442 * t446 + t505 * t663 + t645) * MDP(20) + (-t442 * t448 + t505 * t664 - t526) * MDP(21) + (t419 * t346 + t544 * t515 + t596 * t395 + (t511 * t552 + t515 * t597) * t642 - t628) * MDP(35) + (t419 * t347 + t544 * t511 + t596 * t393 + (t511 * t597 - t515 * t552) * t642 - t629) * MDP(34) + (-t415 * t632 + t596 * t659 - t630) * MDP(27) + (-t415 * t539 - t597 * t659 - t640) * MDP(28) + ((t474 * t514 + (-t472 - t618) * t518) * pkin(7) * MDP(14) + MDP(17) * t529 + (t450 * t518 + t464 * t514) * MDP(13) + (-t450 * t514 + t464 * t518) * MDP(11) + ((t474 - t507) * t514 + (-t472 + t564) * t518) * MDP(12)) * qJD(1) + (qJ(3) * t470 + qJD(3) * t474 - t450 * t464) * MDP(14) + (MDP(9) * t514 * t521 + MDP(10) * t602) * pkin(1) + 0.2e1 * t506 * MDP(13) + (-t592 - t607) * MDP(17) - MDP(4) * t573 + t521 * t634 - t676; -MDP(11) * t573 + (-t508 * t521 - t520) * MDP(13) + (-qJD(2) * t474 + t450 * t586 + t489) * MDP(14) + (t446 * t586 - t513 * t550) * MDP(20) + (-t448 * t586 - t517 * t550) * MDP(21) + (-t586 * t632 - t594 * t659) * MDP(27) + (-t539 * t586 + t595 * t659) * MDP(28) + (-t463 * t360 + t460 * t347 - t594 * t393 + (-t511 * t595 - t515 * t535) * t642) * MDP(34) + (-t463 * t362 + t460 * t346 - t594 * t395 + (t511 * t535 - t515 * t595) * t642) * MDP(35); (t413 + t607) * MDP(17) + (-t505 * t540 - t645) * MDP(20) + (-t505 * t561 + t526) * MDP(21) + (-t354 * t659 + (-t448 * t632 + t580 * t659) * pkin(4) + t630) * MDP(27) + (-t355 * t659 + (-t448 * t539 + t579 * t659) * pkin(4) + t640) * MDP(28) + (t495 * t347 + t543 * t511 + t547 * t393 + (t511 * t546 - t515 * t551) * t642 + t629) * MDP(34) + (t495 * t346 + t543 * t515 + t547 * t395 + (t511 * t551 + t515 * t546) * t642 + t628) * MDP(35) + t676; (-t353 * t659 + t630) * MDP(27) + (-t352 * t659 + t640) * MDP(28) + (-pkin(5) * t347 - (-t352 * t511 + t379 * t515) * t642 - t353 * t393 + t511 * t651 - t599 * pkin(10) + t629) * MDP(34) + (-pkin(5) * t346 + (t352 * t515 + t379 * t511) * t642 - t353 * t395 + t350 * t649 + (t578 * t642 - t362) * pkin(10) + t628) * MDP(35) + t677; t395 * t393 * MDP(29) + (-t393 ^ 2 + t395 ^ 2) * MDP(30) + (t393 * t642 + t598) * MDP(31) + (t395 * t642 + t616) * MDP(32) + t366 * MDP(33) + (t342 * t642 - t350 * t395 + t511 * t549 + t336) * MDP(34) + (-t337 * t511 + t350 * t393 + t515 * t549 - t542 * t642) * MDP(35) + (-MDP(31) * t609 - MDP(32) * t395 - MDP(34) * t342 + MDP(35) * t542) * qJD(6);];
tauc  = t1;
