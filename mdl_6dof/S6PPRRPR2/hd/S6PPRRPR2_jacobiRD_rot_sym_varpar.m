% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PPRRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->14), mult. (80->39), div. (0->0), fcn. (88->10), ass. (0->20)
	t132 = sin(pkin(11));
	t138 = cos(pkin(6));
	t147 = t132 * t138;
	t133 = sin(pkin(7));
	t134 = sin(pkin(6));
	t146 = t133 * t134;
	t145 = t133 * t138;
	t135 = cos(pkin(12));
	t137 = cos(pkin(7));
	t144 = t135 * t137;
	t136 = cos(pkin(11));
	t143 = t136 * t138;
	t131 = sin(pkin(12));
	t142 = -(-t132 * t131 + t135 * t143) * t137 + t136 * t146;
	t141 = -(-t136 * t131 - t135 * t147) * t137 - t132 * t146;
	t140 = cos(qJ(3));
	t139 = sin(qJ(3));
	t130 = -t131 * t147 + t136 * t135;
	t128 = t131 * t143 + t132 * t135;
	t1 = [0, 0, (-t130 * t140 + t141 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t128 * t140 + t142 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t139 * t145 + (-t131 * t140 - t139 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, (t130 * t139 + t141 * t140) * qJD(3), 0, 0, 0; 0, 0, (t128 * t139 + t142 * t140) * qJD(3), 0, 0, 0; 0, 0, (-t140 * t145 + (t131 * t139 - t140 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:35
	% EndTime: 2019-10-09 21:12:35
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (118->35), mult. (402->84), div. (0->0), fcn. (472->12), ass. (0->42)
	t336 = sin(pkin(12));
	t339 = sin(pkin(6));
	t345 = sin(qJ(3));
	t347 = cos(qJ(3));
	t340 = cos(pkin(12));
	t342 = cos(pkin(7));
	t355 = t340 * t342;
	t338 = sin(pkin(7));
	t343 = cos(pkin(6));
	t357 = t338 * t343;
	t328 = (t336 * t347 + t345 * t355) * t339 + t345 * t357;
	t341 = cos(pkin(11));
	t337 = sin(pkin(11));
	t359 = t337 * t343;
	t335 = -t336 * t359 + t341 * t340;
	t334 = -t341 * t336 - t340 * t359;
	t358 = t338 * t339;
	t349 = t334 * t342 + t337 * t358;
	t324 = t335 * t347 + t349 * t345;
	t354 = t341 * t343;
	t333 = t336 * t354 + t337 * t340;
	t332 = -t337 * t336 + t340 * t354;
	t350 = -t332 * t342 + t341 * t358;
	t362 = -t333 * t347 + t350 * t345;
	t356 = t339 * t342;
	t344 = sin(qJ(4));
	t353 = qJD(4) * t344;
	t346 = cos(qJ(4));
	t352 = qJD(4) * t346;
	t321 = -t333 * t345 - t350 * t347;
	t323 = -t335 * t345 + t349 * t347;
	t327 = t347 * t357 + (-t336 * t345 + t347 * t355) * t339;
	t331 = -t340 * t358 + t343 * t342;
	t330 = -t334 * t338 + t337 * t356;
	t329 = -t332 * t338 - t341 * t356;
	t326 = t328 * qJD(3);
	t325 = t327 * qJD(3);
	t320 = t324 * qJD(3);
	t319 = t323 * qJD(3);
	t318 = t362 * qJD(3);
	t317 = t321 * qJD(3);
	t1 = [0, 0, -t320 * t346 - t323 * t353, -t319 * t344 + (-t324 * t346 - t330 * t344) * qJD(4), 0, 0; 0, 0, t318 * t346 - t321 * t353, -t317 * t344 + (-t329 * t344 + t346 * t362) * qJD(4), 0, 0; 0, 0, -t326 * t346 - t327 * t353, -t325 * t344 + (-t328 * t346 - t331 * t344) * qJD(4), 0, 0; 0, 0, t320 * t344 - t323 * t352, -t319 * t346 + (t324 * t344 - t330 * t346) * qJD(4), 0, 0; 0, 0, -t318 * t344 - t321 * t352, -t317 * t346 + (-t329 * t346 - t344 * t362) * qJD(4), 0, 0; 0, 0, t326 * t344 - t327 * t352, -t325 * t346 + (t328 * t344 - t331 * t346) * qJD(4), 0, 0; 0, 0, t319, 0, 0, 0; 0, 0, t317, 0, 0, 0; 0, 0, t325, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:36
	% EndTime: 2019-10-09 21:12:36
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (118->35), mult. (402->84), div. (0->0), fcn. (472->12), ass. (0->42)
	t384 = sin(pkin(12));
	t387 = sin(pkin(6));
	t393 = sin(qJ(3));
	t395 = cos(qJ(3));
	t388 = cos(pkin(12));
	t390 = cos(pkin(7));
	t403 = t388 * t390;
	t386 = sin(pkin(7));
	t391 = cos(pkin(6));
	t405 = t386 * t391;
	t376 = t387 * (t384 * t395 + t393 * t403) + t393 * t405;
	t389 = cos(pkin(11));
	t385 = sin(pkin(11));
	t407 = t385 * t391;
	t383 = -t384 * t407 + t389 * t388;
	t382 = -t389 * t384 - t388 * t407;
	t406 = t386 * t387;
	t397 = t382 * t390 + t385 * t406;
	t372 = t383 * t395 + t393 * t397;
	t402 = t389 * t391;
	t381 = t384 * t402 + t385 * t388;
	t380 = -t385 * t384 + t388 * t402;
	t398 = -t380 * t390 + t389 * t406;
	t410 = -t381 * t395 + t393 * t398;
	t404 = t387 * t390;
	t392 = sin(qJ(4));
	t401 = qJD(4) * t392;
	t394 = cos(qJ(4));
	t400 = qJD(4) * t394;
	t369 = -t381 * t393 - t395 * t398;
	t371 = -t383 * t393 + t395 * t397;
	t375 = t395 * t405 + (-t384 * t393 + t395 * t403) * t387;
	t379 = -t388 * t406 + t391 * t390;
	t378 = -t382 * t386 + t385 * t404;
	t377 = -t380 * t386 - t389 * t404;
	t374 = t376 * qJD(3);
	t373 = t375 * qJD(3);
	t368 = t372 * qJD(3);
	t367 = t371 * qJD(3);
	t366 = t410 * qJD(3);
	t365 = t369 * qJD(3);
	t1 = [0, 0, t367, 0, 0, 0; 0, 0, t365, 0, 0, 0; 0, 0, t373, 0, 0, 0; 0, 0, t368 * t394 + t371 * t401, t367 * t392 + (t372 * t394 + t378 * t392) * qJD(4), 0, 0; 0, 0, -t366 * t394 + t369 * t401, t365 * t392 + (t377 * t392 - t394 * t410) * qJD(4), 0, 0; 0, 0, t374 * t394 + t375 * t401, t373 * t392 + (t376 * t394 + t379 * t392) * qJD(4), 0, 0; 0, 0, -t368 * t392 + t371 * t400, t367 * t394 + (-t372 * t392 + t378 * t394) * qJD(4), 0, 0; 0, 0, t366 * t392 + t369 * t400, t365 * t394 + (t377 * t394 + t392 * t410) * qJD(4), 0, 0; 0, 0, -t374 * t392 + t375 * t400, t373 * t394 + (-t376 * t392 + t379 * t394) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:38
	% EndTime: 2019-10-09 21:12:38
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (438->76), mult. (1401->151), div. (0->0), fcn. (1686->14), ass. (0->66)
	t518 = sin(pkin(12));
	t521 = sin(pkin(6));
	t525 = cos(pkin(6));
	t528 = sin(qJ(3));
	t520 = sin(pkin(7));
	t558 = cos(qJ(3));
	t544 = t520 * t558;
	t522 = cos(pkin(12));
	t524 = cos(pkin(7));
	t553 = t522 * t524;
	t560 = (-t518 * t528 + t558 * t553) * t521 + t525 * t544;
	t519 = sin(pkin(11));
	t523 = cos(pkin(11));
	t552 = t523 * t525;
	t512 = t518 * t552 + t519 * t522;
	t537 = -t519 * t518 + t522 * t552;
	t535 = t537 * t524;
	t555 = t521 * t520;
	t498 = t512 * t558 + (-t523 * t555 + t535) * t528;
	t556 = t519 * t525;
	t554 = t521 * t524;
	t527 = sin(qJ(4));
	t551 = qJD(4) * t527;
	t530 = cos(qJ(4));
	t550 = qJD(4) * t530;
	t526 = sin(qJ(6));
	t549 = qJD(6) * t526;
	t548 = qJD(6) * t527;
	t529 = cos(qJ(6));
	t547 = qJD(6) * t529;
	t543 = t521 * t544;
	t497 = t512 * t528 + t523 * t543 - t558 * t535;
	t493 = t497 * qJD(3);
	t540 = t497 * t548 + t493;
	t513 = -t518 * t556 + t523 * t522;
	t536 = t523 * t518 + t522 * t556;
	t534 = t536 * t524;
	t499 = t513 * t528 - t519 * t543 + t558 * t534;
	t495 = t499 * qJD(3);
	t539 = t499 * t548 + t495;
	t503 = t560 * qJD(3);
	t538 = -t548 * t560 - t503;
	t507 = -t537 * t520 - t523 * t554;
	t490 = t498 * t530 + t507 * t527;
	t489 = t498 * t527 - t507 * t530;
	t500 = t513 * t558 + (t519 * t555 - t534) * t528;
	t508 = t519 * t554 + t536 * t520;
	t492 = t500 * t530 + t508 * t527;
	t491 = t500 * t527 - t508 * t530;
	t506 = t525 * t520 * t528 + (t558 * t518 + t528 * t553) * t521;
	t511 = -t522 * t555 + t525 * t524;
	t502 = t506 * t530 + t511 * t527;
	t501 = t506 * t527 - t511 * t530;
	t494 = t498 * qJD(3);
	t533 = -qJD(6) * t498 - t494 * t527 - t497 * t550;
	t496 = t500 * qJD(3);
	t532 = -qJD(6) * t500 - t496 * t527 - t499 * t550;
	t504 = t506 * qJD(3);
	t531 = -qJD(6) * t506 - t504 * t527 + t550 * t560;
	t488 = -t501 * qJD(4) + t503 * t530;
	t487 = t502 * qJD(4) + t503 * t527;
	t486 = -t491 * qJD(4) - t495 * t530;
	t485 = t492 * qJD(4) - t495 * t527;
	t484 = -t489 * qJD(4) - t493 * t530;
	t483 = t490 * qJD(4) - t493 * t527;
	t1 = [0, 0, t532 * t526 - t539 * t529, t486 * t526 + t492 * t547, 0, t485 * t529 - t496 * t526 + (-t491 * t526 - t499 * t529) * qJD(6); 0, 0, t533 * t526 - t540 * t529, t484 * t526 + t490 * t547, 0, t483 * t529 - t494 * t526 + (-t489 * t526 - t497 * t529) * qJD(6); 0, 0, t531 * t526 - t538 * t529, t488 * t526 + t502 * t547, 0, t487 * t529 - t504 * t526 + (-t501 * t526 + t529 * t560) * qJD(6); 0, 0, t539 * t526 + t532 * t529, t486 * t529 - t492 * t549, 0, -t485 * t526 - t496 * t529 + (-t491 * t529 + t499 * t526) * qJD(6); 0, 0, t540 * t526 + t533 * t529, t484 * t529 - t490 * t549, 0, -t483 * t526 - t494 * t529 + (-t489 * t529 + t497 * t526) * qJD(6); 0, 0, t538 * t526 + t531 * t529, t488 * t529 - t502 * t549, 0, -t487 * t526 - t504 * t529 + (-t501 * t529 - t526 * t560) * qJD(6); 0, 0, -t496 * t530 + t499 * t551, -t485, 0, 0; 0, 0, -t494 * t530 + t497 * t551, -t483, 0, 0; 0, 0, -t504 * t530 - t551 * t560, -t487, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end