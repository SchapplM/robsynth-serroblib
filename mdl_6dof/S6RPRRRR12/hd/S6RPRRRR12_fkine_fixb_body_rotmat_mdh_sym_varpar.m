% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:58
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:12
	% EndTime: 2020-11-04 21:58:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:12
	% EndTime: 2020-11-04 21:58:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t306 = cos(qJ(1));
	t305 = sin(qJ(1));
	t1 = [t306, -t305, 0, 0; t305, t306, 0, 0; 0, 0, 1, pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:12
	% EndTime: 2020-11-04 21:58:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t307 = sin(pkin(14));
	t311 = sin(qJ(1));
	t318 = t311 * t307;
	t308 = sin(pkin(6));
	t317 = t311 * t308;
	t309 = cos(pkin(14));
	t316 = t311 * t309;
	t312 = cos(qJ(1));
	t315 = t312 * t307;
	t314 = t312 * t308;
	t313 = t312 * t309;
	t310 = cos(pkin(6));
	t1 = [-t310 * t318 + t313, -t310 * t316 - t315, t317, t312 * pkin(1) + qJ(2) * t317 + 0; t310 * t315 + t316, t310 * t313 - t318, -t314, t311 * pkin(1) - qJ(2) * t314 + 0; t308 * t307, t308 * t309, t310, t310 * qJ(2) + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:12
	% EndTime: 2020-11-04 21:58:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t327 = sin(pkin(7));
	t350 = pkin(10) * t327;
	t326 = sin(pkin(14));
	t332 = sin(qJ(3));
	t349 = t326 * t332;
	t334 = cos(qJ(3));
	t348 = t326 * t334;
	t328 = sin(pkin(6));
	t347 = t328 * t327;
	t330 = cos(pkin(7));
	t346 = t328 * t330;
	t331 = cos(pkin(6));
	t345 = t331 * t330;
	t344 = t331 * t332;
	t343 = t331 * t334;
	t329 = cos(pkin(14));
	t342 = t332 * t329;
	t333 = sin(qJ(1));
	t341 = t333 * t326;
	t340 = t334 * t329;
	t335 = cos(qJ(1));
	t339 = t335 * t326;
	t338 = t335 * t329;
	t322 = -t326 * pkin(2) + t329 * t350;
	t324 = t330 * pkin(10) + qJ(2);
	t337 = t322 * t331 + t328 * t324;
	t319 = t329 * t345 - t347;
	t336 = t319 * t332 + t326 * t343;
	t321 = t329 * pkin(2) + t326 * t350 + pkin(1);
	t320 = t330 * t349 - t340;
	t1 = [-t335 * t320 - t336 * t333, (-t319 * t333 - t330 * t339) * t334 + t332 * (t331 * t341 - t338), (t333 * t331 * t329 + t339) * t327 + t333 * t346, t321 * t335 + t337 * t333 + 0; -t333 * t320 + t336 * t335, (t319 * t334 - t326 * t344) * t335 - t333 * (t330 * t348 + t342), -(t331 * t338 - t341) * t327 - t335 * t346, t321 * t333 - t337 * t335 + 0; t327 * t344 + (t330 * t342 + t348) * t328, t327 * t343 + (t330 * t340 - t349) * t328, -t329 * t347 + t345, -t322 * t328 + t324 * t331 + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:12
	% EndTime: 2020-11-04 21:58:12
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (133->58), mult. (337->106), div. (0->0), fcn. (426->14), ass. (0->53)
	t377 = cos(pkin(14));
	t379 = cos(pkin(7));
	t380 = cos(pkin(6));
	t397 = t380 * t379;
	t375 = sin(pkin(7));
	t376 = sin(pkin(6));
	t401 = t376 * t375;
	t364 = t377 * t397 - t401;
	t385 = cos(qJ(3));
	t373 = sin(pkin(14));
	t382 = sin(qJ(3));
	t406 = t373 * t382;
	t411 = t364 * t385 - t380 * t406;
	t378 = cos(pkin(8));
	t370 = t378 * pkin(11) + pkin(10);
	t403 = t375 * t370;
	t361 = -t373 * pkin(2) + t377 * t403;
	t374 = sin(pkin(8));
	t409 = pkin(11) * t374;
	t394 = t377 * t409;
	t408 = t373 * pkin(3);
	t366 = t379 * t394 - t408;
	t395 = t373 * t409;
	t400 = t377 * t379;
	t367 = pkin(3) * t400 + t395;
	t369 = t370 * t379 + qJ(2);
	t404 = t374 * t375;
	t396 = pkin(11) * t404;
	t410 = (pkin(3) * t401 - t367 * t380) * t382 - (-t366 * t380 + t376 * t396) * t385 + t361 * t380 + t369 * t376;
	t405 = t373 * t385;
	t402 = t375 * t380;
	t399 = t377 * t382;
	t398 = t379 * t385;
	t388 = t376 * t379 + t377 * t402;
	t354 = -t374 * t388 + t411 * t378;
	t359 = t364 * t382 + t380 * t405;
	t381 = sin(qJ(4));
	t384 = cos(qJ(4));
	t392 = t354 * t381 + t359 * t384;
	t363 = t376 * t400 + t402;
	t389 = -t363 * t385 + t376 * t406;
	t386 = cos(qJ(1));
	t383 = sin(qJ(1));
	t365 = -t385 * t377 + t379 * t406;
	t362 = t377 * t401 - t397;
	t358 = t382 * t363 + t376 * t405;
	t357 = t378 * t399 + (t378 * t398 - t404) * t373;
	t356 = t374 * t399 + (t374 * t398 + t375 * t378) * t373;
	t355 = t357 * t381 + t384 * t365;
	t353 = -t374 * t362 - t389 * t378;
	t352 = t411 * t374 + t378 * t388;
	t351 = (t377 * pkin(3) + t379 * t395) * t385 + (-t379 * t408 + t394) * t382 + t377 * pkin(2) + t373 * t403 + pkin(1);
	t1 = [-t355 * t386 - t392 * t383, (-t354 * t383 - t386 * t357) * t384 + (t359 * t383 + t386 * t365) * t381, t352 * t383 + t386 * t356, t351 * t386 + t410 * t383 + 0; -t383 * t355 + t392 * t386, (t354 * t384 - t381 * t359) * t386 - t383 * (t357 * t384 - t381 * t365), -t352 * t386 + t383 * t356, t351 * t383 - t410 * t386 + 0; t353 * t381 + t384 * t358, t353 * t384 - t358 * t381, -t378 * t362 + t389 * t374, (-t366 * t376 - t380 * t396) * t385 + (pkin(3) * t402 + t367 * t376) * t382 - t361 * t376 + t369 * t380 + 0 + pkin(9); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:12
	% EndTime: 2020-11-04 21:58:13
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (256->103), mult. (695->183), div. (0->0), fcn. (845->16), ass. (0->74)
	t450 = cos(pkin(14));
	t452 = cos(pkin(7));
	t453 = cos(pkin(6));
	t486 = t453 * t452;
	t448 = sin(pkin(7));
	t449 = sin(pkin(6));
	t491 = t449 * t448;
	t431 = t450 * t486 - t491;
	t451 = cos(pkin(8));
	t455 = sin(qJ(4));
	t460 = cos(qJ(3));
	t446 = sin(pkin(14));
	t456 = sin(qJ(3));
	t498 = t446 * t456;
	t481 = t453 * t498;
	t459 = cos(qJ(4));
	t484 = t459 * t456;
	t488 = t451 * t455;
	t492 = t448 * t453;
	t429 = t449 * t452 + t450 * t492;
	t447 = sin(pkin(8));
	t494 = t447 * t429;
	t497 = t446 * t459;
	t414 = -(t451 * t481 + t494) * t455 + (t431 * t488 + t453 * t497) * t460 + t431 * t484;
	t467 = t431 * t460 - t481;
	t418 = t451 * t429 + t467 * t447;
	t454 = sin(qJ(5));
	t458 = cos(qJ(5));
	t506 = t414 * t454 + t458 * t418;
	t442 = t451 * pkin(11) + pkin(10);
	t427 = t448 * t442 * t450 - t446 * pkin(2);
	t490 = t450 * t452;
	t500 = t447 * pkin(11);
	t433 = -t446 * pkin(3) + t490 * t500;
	t477 = t451 * t490;
	t434 = -t446 * pkin(4) + pkin(12) * t477;
	t435 = pkin(4) * t477 + t446 * pkin(12);
	t499 = t446 * t451;
	t436 = -pkin(4) * t499 + pkin(12) * t490;
	t437 = pkin(4) * t490 + pkin(12) * t499;
	t438 = pkin(3) * t490 + t446 * t500;
	t439 = t442 * t452 + qJ(2);
	t479 = t451 * t491;
	t493 = t447 * t448;
	t483 = pkin(11) * t493;
	t501 = pkin(12) * t459;
	t504 = (pkin(4) * t455 - t501) * t447;
	t505 = t427 * t453 + t429 * t504 + t439 * t449 + ((pkin(12) * t491 - t436 * t453) * t455 + (pkin(4) * t491 - t437 * t453) * t459 + pkin(3) * t491 - t438 * t453) * t456 + ((pkin(4) * t479 - t435 * t453) * t455 - (pkin(12) * t479 - t434 * t453) * t459 + t433 * t453 - t449 * t483) * t460;
	t496 = t446 * t460;
	t428 = t450 * t491 - t486;
	t495 = t447 * t428;
	t489 = t450 * t456;
	t487 = t452 * t460;
	t482 = t449 * t498;
	t480 = t452 * t498;
	t478 = t451 * t492;
	t474 = pkin(4) * t459 + pkin(12) * t455 + pkin(3);
	t430 = t449 * t490 + t492;
	t468 = -t430 * t460 + t482;
	t463 = -(t451 * t482 + t495) * t455 + (t430 * t488 + t449 * t497) * t460 + t430 * t484;
	t462 = -pkin(4) * t488 + t451 * t501 + t500;
	t461 = cos(qJ(1));
	t457 = sin(qJ(1));
	t440 = t451 * t489;
	t432 = -t460 * t450 + t480;
	t426 = t431 * t456 + t453 * t496;
	t425 = t440 + (t451 * t487 - t493) * t446;
	t424 = t447 * t489 + (t447 * t487 + t448 * t451) * t446;
	t419 = t467 * t451 - t494;
	t417 = -t428 * t451 + t468 * t447;
	t416 = (t452 * t446 * t488 - t450 * t459) * t460 + (-t446 * t493 + t440) * t455 + t459 * t480;
	t413 = t416 * t454 + t458 * t424;
	t412 = pkin(1) + (t462 * t456 + t474 * t460 + pkin(2)) * t450 + ((t442 + t504) * t448 + (-t474 * t456 + t462 * t460) * t452) * t446;
	t1 = [(-t414 * t457 - t461 * t416) * t458 + (t418 * t457 + t461 * t424) * t454, t413 * t461 + t506 * t457, (t419 * t457 + t461 * t425) * t459 - (t426 * t457 + t461 * t432) * t455, t412 * t461 + t505 * t457 + 0; (t414 * t458 - t454 * t418) * t461 - t457 * (t416 * t458 - t454 * t424), t413 * t457 - t506 * t461, (-t419 * t459 + t455 * t426) * t461 + t457 * (t425 * t459 - t455 * t432), t412 * t457 - t505 * t461 + 0; t417 * t454 + t463 * t458, t417 * t458 - t463 * t454, (t468 * t451 + t495) * t459 + (t456 * t430 + t449 * t496) * t455, ((-pkin(12) * t478 - t434 * t449) * t459 + (pkin(4) * t478 + t435 * t449) * t455 - t433 * t449 - t453 * t483) * t460 + ((pkin(4) * t492 + t437 * t449) * t459 + (pkin(12) * t492 + t436 * t449) * t455 + t438 * t449 + pkin(3) * t492) * t456 - t427 * t449 + t439 * t453 + 0 + pkin(9) - t428 * t504; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:13
	% EndTime: 2020-11-04 21:58:14
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (443->133), mult. (1262->232), div. (0->0), fcn. (1519->18), ass. (0->96)
	t550 = sin(pkin(14));
	t554 = cos(pkin(14));
	t555 = cos(pkin(8));
	t546 = t555 * pkin(11) + pkin(10);
	t552 = sin(pkin(7));
	t616 = t552 * t546;
	t529 = -t550 * pkin(2) + t554 * t616;
	t553 = sin(pkin(6));
	t556 = cos(pkin(7));
	t557 = cos(pkin(6));
	t615 = t552 * t557;
	t531 = t553 * t556 + t554 * t615;
	t604 = t557 * t556;
	t614 = t553 * t552;
	t533 = t554 * t604 - t614;
	t609 = t555 * t556;
	t591 = t554 * t609;
	t536 = -t550 * pkin(4) + pkin(12) * t591;
	t612 = t554 * t556;
	t623 = t550 * t555;
	t538 = -pkin(4) * t623 + pkin(12) * t612;
	t542 = t546 * t556 + qJ(2);
	t551 = sin(pkin(8));
	t560 = sin(qJ(4));
	t561 = sin(qJ(3));
	t565 = cos(qJ(4));
	t566 = cos(qJ(3));
	t559 = sin(qJ(5));
	t564 = cos(qJ(5));
	t584 = pkin(5) * t564 + pkin(13) * t559;
	t569 = pkin(4) * t614 - t584 * t533;
	t629 = pkin(13) * t564;
	t583 = pkin(5) * t559 - t629;
	t624 = t550 * t551;
	t586 = (pkin(4) * t612 + pkin(12) * t623) * t565 + pkin(3) * t612 + pkin(11) * t624;
	t628 = t551 * pkin(11);
	t587 = t550 * pkin(3) + (pkin(4) * t591 + t550 * pkin(12)) * t560 - t612 * t628;
	t597 = pkin(12) * t614;
	t606 = t555 * t564;
	t608 = t555 * t559;
	t618 = t551 * t559;
	t543 = pkin(4) + t584;
	t631 = pkin(12) * t565;
	t642 = t543 * t560 - t631;
	t653 = ((((pkin(5) * t606 + pkin(13) * t608) * t560 - pkin(5) * t618 + t551 * t629) * t561 - t584 * t565 * t566) * t550 - (t538 * t560 + t586) * t561 + (t536 * t565 - t587) * t566 + t529) * t557 + ((-pkin(11) * t614 + t583 * t533) * t566 + t642 * t531) * t551 + ((t569 * t560 - t565 * t597) * t566 + t583 * t531) * t555 + t542 * t553 + (pkin(3) * t614 + t560 * t597 + t569 * t565) * t561;
	t525 = t555 * t531;
	t534 = t560 * t606 - t618;
	t603 = t557 * t561;
	t594 = t550 * t603;
	t521 = t551 * t531 + t555 * t594;
	t600 = t565 * t561;
	t578 = -t521 * t560 + t533 * t600;
	t620 = t550 * t565;
	t589 = t557 * t620;
	t509 = t578 * t564 + (t534 * t533 + t564 * t589) * t566 - t559 * (-t551 * t594 + t525);
	t602 = t561 * t560;
	t605 = t555 * t565;
	t621 = t550 * t560;
	t514 = (t533 * t605 - t557 * t621) * t566 - t521 * t565 - t533 * t602;
	t558 = sin(qJ(6));
	t563 = cos(qJ(6));
	t652 = t509 * t558 + t514 * t563;
	t640 = pkin(12) * t560 + t543 * t565 + pkin(3);
	t607 = t555 * t560;
	t639 = ((t533 * t607 + t589) * t566 + t578) * t559 + t564 * (t525 + (t533 * t566 - t594) * t551);
	t544 = t552 * t623;
	t611 = t554 * t561;
	t527 = t551 * t611 + t544;
	t593 = t552 * t624;
	t528 = t555 * t611 - t593;
	t622 = t550 * t556;
	t590 = t561 * t622;
	t572 = t528 * t560 + t565 * t590;
	t610 = t554 * t565;
	t511 = t572 * t564 - (-t534 * t622 + t564 * t610) * t566 - t559 * t527;
	t530 = -t554 * t614 + t604;
	t613 = t553 * t561;
	t596 = t550 * t613;
	t520 = t530 * t555 + t551 * t596;
	t532 = t553 * t612 + t615;
	t619 = t551 * t530;
	t522 = -t555 * t596 + t619;
	t577 = t522 * t560 + t532 * t600;
	t592 = t553 * t620;
	t636 = t577 * t564 + (t534 * t532 + t564 * t592) * t566 + t520 * t559;
	t617 = t551 * t566;
	t595 = t550 * t609;
	t568 = (pkin(11) + t583) * t551 - t642 * t555;
	t567 = cos(qJ(1));
	t562 = sin(qJ(1));
	t515 = (t554 * t560 + t565 * t595) * t566 + t528 * t565 - t560 * t590;
	t513 = (t532 * t605 - t553 * t621) * t566 + t522 * t565 - t532 * t602;
	t512 = ((t560 * t595 - t610) * t566 + t572) * t559 + t564 * (t617 * t622 + t527);
	t508 = t511 * t558 + t563 * t515;
	t507 = (t640 * t554 + t568 * t622) * t566 + (t568 * t554 - t640 * t622) * t561 + t554 * pkin(2) + t550 * t616 + pkin(1) + t642 * t593 + t583 * t544;
	t1 = [(-t509 * t562 - t511 * t567) * t563 + (t514 * t562 + t515 * t567) * t558, t508 * t567 + t652 * t562, -t512 * t567 - t639 * t562, t507 * t567 + t653 * t562 + 0; (t509 * t563 - t514 * t558) * t567 + (-t511 * t563 + t515 * t558) * t562, t562 * t508 - t652 * t567, -t512 * t562 + t639 * t567, t507 * t562 - t653 * t567 + 0; -t513 * t558 + t636 * t563, -t563 * t513 - t636 * t558, ((t532 * t607 + t592) * t566 + t577) * t559 - (-t532 * t617 + t520) * t564, ((-t628 + (pkin(4) * t560 - t631) * t555) * t615 + ((t584 * t550 - t536) * t565 + t587) * t553 + (-t583 * t551 + t584 * t607) * t532) * t566 + (t552 * pkin(12) * t603 + t543 * t619 + (-t584 * t623 + t538) * t613) * t560 + ((pkin(4) * t615 + t584 * t532) * t565 + pkin(3) * t615 + (t583 * t624 + t586) * t553) * t561 - t619 * t631 - t529 * t553 + t542 * t557 + 0 + pkin(9) + (pkin(5) * t608 - pkin(13) * t606) * t530; 0, 0, 0, 1;];
	Tc_mdh = t1;
end