% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRR15_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:25
	% EndTime: 2019-10-10 12:18:25
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (178->54), mult. (587->106), div. (0->0), fcn. (611->10), ass. (0->49)
	t364 = cos(pkin(6));
	t366 = sin(qJ(2));
	t367 = sin(qJ(1));
	t392 = t367 * t366;
	t385 = t364 * t392;
	t369 = cos(qJ(2));
	t370 = cos(qJ(1));
	t388 = t370 * t369;
	t351 = -qJD(1) * t385 - qJD(2) * t392 + (qJD(2) * t364 + qJD(1)) * t388;
	t389 = t370 * t366;
	t391 = t367 * t369;
	t353 = t364 * t389 + t391;
	t361 = sin(pkin(7));
	t365 = sin(qJ(3));
	t368 = cos(qJ(3));
	t375 = t364 * t391 + t389;
	t350 = t375 * qJD(1) + t353 * qJD(2);
	t363 = cos(pkin(7));
	t362 = sin(pkin(6));
	t387 = qJD(1) * t362;
	t384 = t367 * t387;
	t372 = -t350 * t363 + t361 * t384;
	t398 = t362 * t370;
	t352 = -t364 * t388 + t392;
	t401 = t352 * t363;
	t405 = ((t361 * t398 + t401) * t365 - t353 * t368) * qJD(3) - t351 * t365 + t372 * t368;
	t399 = t361 * t362;
	t397 = t363 * t365;
	t396 = t363 * t368;
	t395 = t365 * t366;
	t394 = t365 * t369;
	t393 = t366 * t368;
	t390 = t368 * t369;
	t386 = qJD(3) * t361;
	t383 = t370 * t387;
	t382 = t368 * t386;
	t380 = -t363 * t375 + t367 * t399;
	t379 = -t363 * t390 + t395;
	t378 = -t363 * t393 - t394;
	t377 = -t363 * t394 - t393;
	t376 = t363 * t395 - t390;
	t374 = t385 - t388;
	t348 = t352 * qJD(1) + t374 * qJD(2);
	t373 = t348 * t363 + t361 * t383;
	t371 = t382 * t398 + (qJD(3) * t401 - t351) * t368 + (qJD(3) * t353 - t372) * t365;
	t349 = t353 * qJD(1) + t375 * qJD(2);
	t347 = -t349 * t368 + t373 * t365 + (t365 * t374 + t380 * t368) * qJD(3);
	t346 = t349 * t365 + t373 * t368 + (-t380 * t365 + t368 * t374) * qJD(3);
	t1 = [t371, t349 * t397 + t348 * t368 + (t365 * t375 + t374 * t396) * qJD(3), t346, 0, 0, 0; t347, -t351 * t397 - t350 * t368 + (t352 * t365 - t353 * t396) * qJD(3), t405, 0, 0, 0; 0, (t377 * qJD(2) + t378 * qJD(3)) * t362, -t364 * t365 * t386 + (t378 * qJD(2) + t377 * qJD(3)) * t362, 0, 0, 0; -t405, t349 * t396 - t348 * t365 + (t368 * t375 - t374 * t397) * qJD(3), -t347, 0, 0, 0; t346, -t351 * t396 + t350 * t365 + (t352 * t368 + t353 * t397) * qJD(3), t371, 0, 0, 0; 0, (t379 * qJD(2) + t376 * qJD(3)) * t362, -t364 * t382 + (t376 * qJD(2) + t379 * qJD(3)) * t362, 0, 0, 0; -t350 * t361 - t363 * t384, -t349 * t361, 0, 0, 0, 0; -t348 * t361 + t363 * t383, t351 * t361, 0, 0, 0, 0; 0, qJD(2) * t369 * t399, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:26
	% EndTime: 2019-10-10 12:18:27
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (178->54), mult. (587->104), div. (0->0), fcn. (611->10), ass. (0->47)
	t473 = cos(pkin(6));
	t476 = sin(qJ(1));
	t478 = cos(qJ(2));
	t499 = t476 * t478;
	t475 = sin(qJ(2));
	t479 = cos(qJ(1));
	t500 = t475 * t479;
	t462 = t473 * t500 + t499;
	t484 = t473 * t499 + t500;
	t459 = t484 * qJD(1) + t462 * qJD(2);
	t502 = t475 * t476;
	t495 = t473 * t502;
	t497 = t478 * t479;
	t460 = -qJD(1) * t495 - qJD(2) * t502 + (qJD(2) * t473 + qJD(1)) * t497;
	t472 = cos(pkin(7));
	t474 = sin(qJ(3));
	t477 = cos(qJ(3));
	t461 = -t473 * t497 + t502;
	t470 = sin(pkin(7));
	t471 = sin(pkin(6));
	t507 = t470 * t471;
	t490 = t461 * t472 + t479 * t507;
	t496 = qJD(1) * t471;
	t494 = t476 * t496;
	t491 = t470 * t494;
	t512 = (-t462 * t477 + t490 * t474) * qJD(3) + (-t459 * t472 + t491) * t477 - t460 * t474;
	t506 = t472 * t474;
	t505 = t472 * t477;
	t504 = t474 * t475;
	t503 = t474 * t478;
	t501 = t475 * t477;
	t498 = t477 * t478;
	t493 = t479 * t496;
	t492 = qJD(3) * t470 * t473;
	t489 = -t472 * t484 + t476 * t507;
	t488 = t472 * t498 - t504;
	t487 = t472 * t501 + t503;
	t486 = t472 * t503 + t501;
	t485 = -t472 * t504 + t498;
	t483 = t495 - t497;
	t457 = t461 * qJD(1) + t483 * qJD(2);
	t482 = t457 * t472 + t470 * t493;
	t480 = -t459 * t506 + t460 * t477 + t474 * t491 + (-t462 * t474 - t490 * t477) * qJD(3);
	t458 = t462 * qJD(1) + t484 * qJD(2);
	t456 = -t458 * t477 + t482 * t474 + (t474 * t483 + t489 * t477) * qJD(3);
	t455 = -t458 * t474 - t482 * t477 + (t489 * t474 - t477 * t483) * qJD(3);
	t1 = [-t459 * t470 - t472 * t494, -t458 * t470, 0, 0, 0, 0; -t457 * t470 + t472 * t493, t460 * t470, 0, 0, 0, 0; 0, qJD(2) * t478 * t507, 0, 0, 0, 0; t480, -t458 * t506 - t457 * t477 + (-t474 * t484 - t483 * t505) * qJD(3), t455, 0, 0, 0; -t456, t460 * t506 + t459 * t477 + (-t461 * t474 + t462 * t505) * qJD(3), -t512, 0, 0, 0; 0, (t486 * qJD(2) + t487 * qJD(3)) * t471, t474 * t492 + (t487 * qJD(2) + t486 * qJD(3)) * t471, 0, 0, 0; t512, -t458 * t505 + t457 * t474 + (-t477 * t484 + t483 * t506) * qJD(3), t456, 0, 0, 0; t455, t460 * t505 - t459 * t474 + (-t461 * t477 - t462 * t506) * qJD(3), t480, 0, 0, 0; 0, (t488 * qJD(2) + t485 * qJD(3)) * t471, t477 * t492 + (t485 * qJD(2) + t488 * qJD(3)) * t471, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:28
	% EndTime: 2019-10-10 12:18:29
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (493->104), mult. (1577->198), div. (0->0), fcn. (1713->12), ass. (0->85)
	t598 = sin(qJ(2));
	t599 = sin(qJ(1));
	t602 = cos(qJ(2));
	t603 = cos(qJ(1));
	t649 = cos(pkin(6));
	t622 = t603 * t649;
	t584 = t598 * t622 + t599 * t602;
	t623 = t599 * t649;
	t605 = t603 * t598 + t602 * t623;
	t574 = t605 * qJD(1) + t584 * qJD(2);
	t617 = t598 * t623;
	t632 = t603 * t602;
	t634 = t599 * t598;
	t575 = -qJD(1) * t617 - qJD(2) * t634 + (qJD(2) * t649 + qJD(1)) * t632;
	t583 = -t602 * t622 + t634;
	t597 = sin(qJ(3));
	t601 = cos(qJ(3));
	t593 = sin(pkin(7));
	t594 = sin(pkin(6));
	t631 = qJD(1) * t594;
	t626 = t599 * t631;
	t620 = t593 * t626;
	t640 = t594 * t603;
	t627 = t593 * t640;
	t595 = cos(pkin(7));
	t638 = t595 * t601;
	t639 = t595 * t597;
	t646 = t584 * t601;
	t556 = (t583 * t639 - t646) * qJD(3) - t574 * t638 + t601 * t620 + (qJD(3) * t627 - t575) * t597;
	t614 = t583 * t595 + t627;
	t561 = t584 * t597 + t614 * t601;
	t568 = t574 * t593 + t595 * t626;
	t578 = -t583 * t593 + t595 * t640;
	t596 = sin(qJ(5));
	t600 = cos(qJ(5));
	t660 = (t561 * t596 - t578 * t600) * qJD(5) + t556 * t600 + t568 * t596;
	t657 = t561 * qJD(3) - (-t574 * t595 + t620) * t597 - t575 * t601;
	t656 = t556 * t596 - t568 * t600 + (-t561 * t600 - t578 * t596) * qJD(5);
	t606 = t617 - t632;
	t641 = t594 * t599;
	t613 = t593 * t641 - t595 * t605;
	t650 = t597 * t606 + t613 * t601;
	t644 = t593 * t594;
	t643 = t593 * t596;
	t642 = t593 * t600;
	t637 = t597 * t598;
	t636 = t597 * t602;
	t635 = t598 * t601;
	t633 = t601 * t602;
	t630 = qJD(5) * t596;
	t629 = qJD(5) * t600;
	t628 = t598 * t644;
	t625 = t603 * t631;
	t624 = qJD(2) * t644;
	t621 = t649 * t593;
	t619 = t598 * t624;
	t618 = t602 * t624;
	t616 = qJD(3) * t621;
	t570 = -t583 * t597 + t584 * t638;
	t571 = -t597 * t605 - t606 * t638;
	t612 = t595 * t633 - t637;
	t611 = t595 * t635 + t636;
	t610 = t595 * t636 + t635;
	t609 = -t595 * t637 + t633;
	t572 = t583 * qJD(1) + t606 * qJD(2);
	t608 = t572 * t595 + t593 * t625;
	t565 = t613 * t597 - t601 * t606;
	t582 = t649 * t595 - t602 * t644;
	t581 = t611 * t594;
	t580 = t593 * t605 + t595 * t641;
	t577 = t610 * t594 + t597 * t621;
	t576 = -t612 * t594 - t601 * t621;
	t573 = t584 * qJD(1) + t605 * qJD(2);
	t569 = (t612 * qJD(2) + t609 * qJD(3)) * t594;
	t566 = -t572 * t593 + t595 * t625;
	t562 = -t614 * t597 + t646;
	t560 = t601 * t616 + (t609 * qJD(2) + t612 * qJD(3)) * t594;
	t559 = t597 * t616 + (t611 * qJD(2) + t610 * qJD(3)) * t594;
	t558 = t575 * t638 - t574 * t597 + (-t583 * t601 - t584 * t639) * qJD(3);
	t557 = -t573 * t638 + t572 * t597 + (-t601 * t605 + t606 * t639) * qJD(3);
	t553 = t650 * qJD(3) - t573 * t601 + t608 * t597;
	t552 = t565 * qJD(3) - t573 * t597 - t608 * t601;
	t551 = t552 * t596 + t566 * t600 + (-t580 * t596 - t600 * t650) * qJD(5);
	t550 = t552 * t600 - t566 * t596 + (-t580 * t600 + t596 * t650) * qJD(5);
	t1 = [t656, -t573 * t642 + t557 * t596 + (t571 * t600 + t606 * t643) * qJD(5), t553 * t596 + t565 * t629, 0, t550, 0; t551, t575 * t642 + t558 * t596 + (t570 * t600 - t584 * t643) * qJD(5), t562 * t629 - t596 * t657, 0, -t660, 0; 0, t600 * t618 + t569 * t596 + (t581 * t600 - t596 * t628) * qJD(5), t560 * t596 + t577 * t629, 0, -t596 * t619 + t559 * t600 + (-t576 * t596 - t582 * t600) * qJD(5), 0; t660, t573 * t643 + t557 * t600 + (-t571 * t596 + t606 * t642) * qJD(5), t553 * t600 - t565 * t630, 0, -t551, 0; t550, -t575 * t643 + t558 * t600 + (-t570 * t596 - t584 * t642) * qJD(5), -t562 * t630 - t600 * t657, 0, t656, 0; 0, -t596 * t618 + t569 * t600 + (-t581 * t596 - t600 * t628) * qJD(5), t560 * t600 - t577 * t630, 0, -t600 * t619 - t559 * t596 + (-t576 * t600 + t582 * t596) * qJD(5), 0; t657, -t571 * qJD(3) + t572 * t601 + t573 * t639, -t552, 0, 0, 0; t553, -t570 * qJD(3) - t574 * t601 - t575 * t639, t556, 0, 0, 0; 0, (-t610 * qJD(2) - t611 * qJD(3)) * t594, -t559, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:35
	% EndTime: 2019-10-10 12:18:37
	% DurationCPUTime: 2.26s
	% Computational Cost: add. (1250->170), mult. (3901->307), div. (0->0), fcn. (4371->14), ass. (0->132)
	t854 = sin(qJ(2));
	t855 = sin(qJ(1));
	t859 = cos(qJ(2));
	t860 = cos(qJ(1));
	t917 = cos(pkin(6));
	t884 = t860 * t917;
	t835 = t854 * t884 + t855 * t859;
	t885 = t855 * t917;
	t865 = t860 * t854 + t859 * t885;
	t822 = t865 * qJD(1) + t835 * qJD(2);
	t879 = t854 * t885;
	t900 = t860 * t859;
	t902 = t855 * t854;
	t823 = -qJD(1) * t879 - qJD(2) * t902 + (qJD(2) * t917 + qJD(1)) * t900;
	t850 = cos(pkin(7));
	t853 = sin(qJ(3));
	t858 = cos(qJ(3));
	t848 = sin(pkin(7));
	t849 = sin(pkin(6));
	t899 = qJD(1) * t849;
	t889 = t855 * t899;
	t883 = t848 * t889;
	t834 = -t859 * t884 + t902;
	t908 = t849 * t860;
	t891 = t848 * t908;
	t873 = t834 * t850 + t891;
	t915 = t835 * t853;
	t918 = t873 * t858 + t915;
	t783 = t918 * qJD(3) - (-t822 * t850 + t883) * t853 - t823 * t858;
	t851 = sin(qJ(6));
	t929 = t783 * t851;
	t856 = cos(qJ(6));
	t928 = t783 * t856;
	t914 = t835 * t858;
	t805 = t873 * t853 - t914;
	t927 = t805 * t851;
	t926 = t805 * t856;
	t906 = t850 * t858;
	t907 = t850 * t853;
	t782 = (t834 * t907 - t914) * qJD(3) - t822 * t906 + t858 * t883 - (-qJD(3) * t891 + t823) * t853;
	t810 = t822 * t848 + t850 * t889;
	t852 = sin(qJ(5));
	t857 = cos(qJ(5));
	t925 = t782 * t852 - t810 * t857;
	t924 = -t782 * t857 - t810 * t852;
	t826 = -t834 * t848 + t850 * t908;
	t921 = t826 * t852;
	t920 = t826 * t857;
	t901 = t858 * t859;
	t905 = t853 * t854;
	t871 = t850 * t901 - t905;
	t866 = t879 - t900;
	t913 = t866 * t853;
	t912 = t848 * t849;
	t911 = t848 * t852;
	t910 = t848 * t857;
	t909 = t849 * t855;
	t904 = t853 * t859;
	t903 = t854 * t858;
	t898 = qJD(5) * t852;
	t897 = qJD(5) * t857;
	t896 = qJD(6) * t851;
	t895 = qJD(6) * t852;
	t894 = qJD(6) * t856;
	t893 = t854 * t912;
	t892 = t848 * t909;
	t888 = t860 * t899;
	t887 = qJD(2) * t912;
	t886 = t848 * t917;
	t882 = t848 * t888;
	t881 = t854 * t887;
	t880 = t859 * t887;
	t878 = qJD(3) * t886;
	t872 = -t850 * t865 + t892;
	t807 = t872 * t853 - t858 * t866;
	t820 = t834 * qJD(1) + t866 * qJD(2);
	t821 = t835 * qJD(1) + t865 * qJD(2);
	t778 = t807 * qJD(3) - t820 * t906 - t821 * t853 - t858 * t882;
	t877 = -t807 * t895 - t778;
	t876 = t805 * t895 + t782;
	t869 = t850 * t904 + t903;
	t870 = t850 * t903 + t904;
	t798 = t853 * t878 + (t870 * qJD(2) + t869 * qJD(3)) * t849;
	t825 = t869 * t849 + t853 * t886;
	t875 = -t825 * t895 - t798;
	t802 = t834 * t906 + t858 * t891 + t915;
	t792 = t802 * t857 + t921;
	t793 = t802 * t852 - t920;
	t791 = -t852 * t918 + t920;
	t806 = -t858 * t892 + t865 * t906 - t913;
	t828 = t848 * t865 + t850 * t909;
	t794 = t806 * t857 - t828 * t852;
	t795 = t806 * t852 + t828 * t857;
	t824 = -t871 * t849 - t858 * t886;
	t833 = t917 * t850 - t859 * t912;
	t800 = t824 * t857 - t833 * t852;
	t801 = t824 * t852 + t833 * t857;
	t815 = -t834 * t853 + t835 * t906;
	t796 = t815 * t852 + t835 * t910;
	t817 = -t853 * t865 - t866 * t906;
	t797 = t817 * t852 - t866 * t910;
	t816 = -t834 * t858 - t835 * t907;
	t818 = -t858 * t865 + t866 * t907;
	t868 = -t850 * t905 + t901;
	t831 = t870 * t849;
	t819 = t831 * t852 + t857 * t893;
	t779 = -t821 * t858 + (t820 * t850 + t882) * t853 + (t872 * t858 + t913) * qJD(3);
	t863 = -qJD(6) * t806 + t779 * t852 + t807 * t897;
	t862 = -qJD(6) * t802 - t783 * t852 - t805 * t897;
	t799 = t858 * t878 + (t868 * qJD(2) + t871 * qJD(3)) * t849;
	t861 = -qJD(6) * t824 + t799 * t852 + t825 * t897;
	t832 = t868 * t849;
	t812 = (-t869 * qJD(2) - t870 * qJD(3)) * t849;
	t811 = (t871 * qJD(2) + t868 * qJD(3)) * t849;
	t808 = -t820 * t848 + t850 * t888;
	t790 = t857 * t880 + t811 * t852 + (t831 * t857 - t852 * t893) * qJD(5);
	t789 = -t815 * qJD(3) - t822 * t858 - t823 * t907;
	t788 = t816 * qJD(3) - t822 * t853 + t823 * t906;
	t787 = -t817 * qJD(3) + t820 * t858 + t821 * t907;
	t786 = t818 * qJD(3) + t820 * t853 - t821 * t906;
	t785 = t800 * qJD(5) + t798 * t852 + t857 * t881;
	t784 = -t801 * qJD(5) + t798 * t857 - t852 * t881;
	t776 = t823 * t910 + t788 * t852 + (t815 * t857 - t835 * t911) * qJD(5);
	t775 = -t821 * t910 + t786 * t852 + (t817 * t857 + t866 * t911) * qJD(5);
	t774 = t792 * qJD(5) - t925;
	t773 = -t793 * qJD(5) + t924;
	t772 = (-t857 * t918 - t921) * qJD(5) + t925;
	t771 = t794 * qJD(5) + t778 * t852 + t808 * t857;
	t770 = t795 * qJD(5) - t778 * t857 + t808 * t852;
	t769 = t771 * t856 + t779 * t851 + (-t795 * t851 + t807 * t856) * qJD(6);
	t768 = -t771 * t851 + t779 * t856 + (-t795 * t856 - t807 * t851) * qJD(6);
	t1 = [t772 * t856 + t929 + (-t791 * t851 + t926) * qJD(6), t775 * t856 + t787 * t851 + (-t797 * t851 + t818 * t856) * qJD(6), t877 * t851 + t863 * t856, 0, -t770 * t856 - t794 * t896, t768; t769, t776 * t856 + t789 * t851 + (-t796 * t851 + t816 * t856) * qJD(6), t876 * t851 + t862 * t856, 0, t773 * t856 - t792 * t896, -t774 * t851 - t928 + (-t793 * t856 + t927) * qJD(6); 0, t790 * t856 + t812 * t851 + (-t819 * t851 + t832 * t856) * qJD(6), t875 * t851 + t861 * t856, 0, t784 * t856 - t800 * t896, -t785 * t851 + t799 * t856 + (-t801 * t856 - t825 * t851) * qJD(6); -t772 * t851 + t928 + (-t791 * t856 - t927) * qJD(6), -t775 * t851 + t787 * t856 + (-t797 * t856 - t818 * t851) * qJD(6), -t863 * t851 + t877 * t856, 0, t770 * t851 - t794 * t894, -t769; t768, -t776 * t851 + t789 * t856 + (-t796 * t856 - t816 * t851) * qJD(6), -t862 * t851 + t876 * t856, 0, -t773 * t851 - t792 * t894, -t774 * t856 + t929 + (t793 * t851 + t926) * qJD(6); 0, -t790 * t851 + t812 * t856 + (-t819 * t856 - t832 * t851) * qJD(6), -t861 * t851 + t875 * t856, 0, -t784 * t851 - t800 * t894, -t785 * t856 - t799 * t851 + (t801 * t851 - t825 * t856) * qJD(6); t791 * qJD(5) + t924, t797 * qJD(5) - t786 * t857 - t821 * t911, -t779 * t857 + t807 * t898, 0, t771, 0; t770, t796 * qJD(5) - t788 * t857 + t823 * t911, t783 * t857 - t805 * t898, 0, t774, 0; 0, t819 * qJD(5) - t811 * t857 + t852 * t880, -t799 * t857 + t825 * t898, 0, t785, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end