% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
	t139 = sin(pkin(6));
	t152 = t139 * cos(pkin(7));
	t142 = sin(qJ(2));
	t143 = sin(qJ(1));
	t151 = t142 * t143;
	t145 = cos(qJ(1));
	t150 = t142 * t145;
	t144 = cos(qJ(2));
	t149 = t143 * t144;
	t148 = t144 * t145;
	t147 = qJD(1) * t139;
	t138 = sin(pkin(7));
	t146 = qJD(2) * t138;
	t141 = cos(pkin(6));
	t1 = [0, t145 * t147, -(t141 * t151 - t148) * t146 + (-(-t141 * t148 + t151) * t138 + t145 * t152) * qJD(1), 0, 0, 0; 0, t143 * t147, -(-t141 * t150 - t149) * t146 + (-(-t141 * t149 - t150) * t138 + t143 * t152) * qJD(1), 0, 0, 0; 0, 0, t139 * t142 * t146, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:54
	% EndTime: 2019-10-10 13:34:55
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (49->27), mult. (173->62), div. (0->0), fcn. (177->12), ass. (0->34)
	t270 = sin(pkin(7));
	t271 = sin(pkin(6));
	t297 = t271 * t270;
	t273 = cos(pkin(7));
	t278 = cos(qJ(3));
	t296 = t273 * t278;
	t275 = sin(qJ(3));
	t279 = cos(qJ(2));
	t295 = t275 * t279;
	t276 = sin(qJ(2));
	t294 = t276 * t278;
	t277 = sin(qJ(1));
	t293 = t277 * t276;
	t292 = t277 * t279;
	t280 = cos(qJ(1));
	t291 = t280 * t276;
	t290 = t280 * t279;
	t289 = qJD(1) * t271;
	t269 = sin(pkin(8));
	t288 = qJD(3) * t269;
	t287 = t277 * t289;
	t286 = t280 * t289;
	t285 = t270 * t278 * t289;
	t274 = cos(pkin(6));
	t284 = t274 * t290 - t293;
	t283 = -t274 * t292 - t291;
	t282 = t274 * t291 + t292;
	t281 = t274 * t293 - t290;
	t272 = cos(pkin(8));
	t268 = t283 * qJD(1) - t282 * qJD(2);
	t267 = -t284 * qJD(1) + t281 * qJD(2);
	t266 = -t268 * t270 + t273 * t287;
	t265 = -t267 * t270 + t273 * t286;
	t1 = [0, t286, t265, -(-(-t282 * qJD(1) + t283 * qJD(2)) * t275 + t267 * t296 + t280 * t285) * t269 + t265 * t272 - (t281 * t278 + (-t283 * t273 - t277 * t297) * t275) * t288, 0, 0; 0, t287, t266, -(-(-t281 * qJD(1) + t284 * qJD(2)) * t275 + t268 * t296 + t277 * t285) * t269 + t266 * t272 - (-t282 * t278 + (-t284 * t273 + t280 * t297) * t275) * t288, 0, 0; 0, 0, qJD(2) * t276 * t297, t274 * t270 * t275 * t288 + (-(-t273 * t295 - t294) * t288 + (-(-t273 * t294 - t295) * t269 + t276 * t270 * t272) * qJD(2)) * t271, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:56
	% EndTime: 2019-10-10 13:34:57
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (157->57), mult. (522->120), div. (0->0), fcn. (559->14), ass. (0->56)
	t387 = sin(qJ(3));
	t391 = cos(qJ(3));
	t385 = cos(pkin(6));
	t392 = cos(qJ(2));
	t393 = cos(qJ(1));
	t409 = t393 * t392;
	t388 = sin(qJ(2));
	t389 = sin(qJ(1));
	t413 = t389 * t388;
	t398 = t385 * t413 - t409;
	t410 = t393 * t388;
	t412 = t389 * t392;
	t378 = -t385 * t412 - t410;
	t381 = sin(pkin(7));
	t384 = cos(pkin(7));
	t382 = sin(pkin(6));
	t418 = t382 * t389;
	t401 = t378 * t384 + t381 * t418;
	t424 = t401 * t387 - t391 * t398;
	t377 = t385 * t410 + t412;
	t376 = t385 * t409 - t413;
	t417 = t382 * t393;
	t402 = -t376 * t384 + t381 * t417;
	t423 = -t377 * t391 + t402 * t387;
	t420 = t381 * t385;
	t419 = t381 * t388;
	t416 = t387 * t388;
	t415 = t387 * t392;
	t414 = t388 * t391;
	t411 = t391 * t392;
	t408 = qJD(1) * t382;
	t386 = sin(qJ(4));
	t407 = qJD(3) * t386;
	t406 = t389 * t408;
	t405 = t393 * t408;
	t404 = qJD(3) * t420;
	t403 = t382 * qJD(2) * t419;
	t400 = t384 * t411 - t416;
	t399 = t384 * t415 + t414;
	t372 = -t376 * qJD(1) + t398 * qJD(2);
	t397 = t372 * t384 + t381 * t405;
	t374 = t378 * qJD(1) - t377 * qJD(2);
	t396 = t374 * t384 + t381 * t406;
	t395 = -t377 * t387 - t402 * t391;
	t394 = t387 * t398 + t401 * t391;
	t390 = cos(qJ(4));
	t383 = cos(pkin(8));
	t380 = sin(pkin(8));
	t375 = -t398 * qJD(1) + t376 * qJD(2);
	t373 = -t377 * qJD(1) + t378 * qJD(2);
	t371 = -t374 * t381 + t384 * t406;
	t370 = -t372 * t381 + t384 * t405;
	t369 = -t387 * t404 + (-t399 * qJD(3) + (-t384 * t414 - t415) * qJD(2)) * t382;
	t368 = t423 * qJD(3) - t375 * t387 + t396 * t391;
	t367 = -t424 * qJD(3) - t373 * t387 + t397 * t391;
	t1 = [0, t405, t370, -t367 * t380 + t370 * t383, (t373 * t391 + t387 * t397) * t386 + (-t367 * t383 - t370 * t380) * t390 + t394 * t407 + (t424 * t390 + (t394 * t383 + (-t378 * t381 + t384 * t418) * t380) * t386) * qJD(4), 0; 0, t406, t371, -t368 * t380 + t371 * t383, (t375 * t391 + t387 * t396) * t386 + (-t368 * t383 - t371 * t380) * t390 + t395 * t407 + (-t423 * t390 + (t395 * t383 + (-t376 * t381 - t384 * t417) * t380) * t386) * qJD(4), 0; 0, 0, t403, -t369 * t380 + t383 * t403, t391 * t386 * t404 - t369 * t383 * t390 + (t400 * t407 + ((-t384 * t416 + t411) * t386 - t380 * t390 * t419) * qJD(2)) * t382 + ((t382 * t399 + t387 * t420) * t390 + ((t382 * t400 + t391 * t420) * t383 + (-t382 * t392 * t381 + t385 * t384) * t380) * t386) * qJD(4), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:59
	% EndTime: 2019-10-10 13:35:00
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (364->81), mult. (1171->158), div. (0->0), fcn. (1294->16), ass. (0->78)
	t504 = sin(qJ(3));
	t509 = cos(qJ(3));
	t501 = cos(pkin(6));
	t510 = cos(qJ(2));
	t511 = cos(qJ(1));
	t532 = t511 * t510;
	t505 = sin(qJ(2));
	t506 = sin(qJ(1));
	t536 = t506 * t505;
	t517 = t501 * t536 - t532;
	t533 = t511 * t505;
	t535 = t506 * t510;
	t494 = -t501 * t535 - t533;
	t497 = sin(pkin(7));
	t500 = cos(pkin(7));
	t498 = sin(pkin(6));
	t542 = t498 * t506;
	t520 = t494 * t500 + t497 * t542;
	t480 = t504 * t520 - t509 * t517;
	t493 = t501 * t533 + t535;
	t492 = t501 * t532 - t536;
	t541 = t498 * t511;
	t521 = -t492 * t500 + t497 * t541;
	t548 = -t493 * t509 + t504 * t521;
	t496 = sin(pkin(8));
	t503 = sin(qJ(4));
	t545 = t496 * t503;
	t544 = t497 * t498;
	t543 = t497 * t501;
	t499 = cos(pkin(8));
	t540 = t499 * t503;
	t539 = t504 * t505;
	t538 = t504 * t510;
	t537 = t505 * t509;
	t534 = t509 * t510;
	t531 = qJD(1) * t498;
	t502 = sin(qJ(5));
	t530 = qJD(4) * t502;
	t529 = t506 * t531;
	t528 = t511 * t531;
	t527 = qJD(3) * t543;
	t526 = qJD(2) * t505 * t544;
	t525 = t496 * t526;
	t477 = -t493 * t504 - t509 * t521;
	t489 = -t492 * t497 - t500 * t541;
	t524 = t477 * t499 + t489 * t496;
	t479 = t504 * t517 + t509 * t520;
	t490 = -t494 * t497 + t500 * t542;
	t523 = t479 * t499 + t490 * t496;
	t519 = t500 * t534 - t539;
	t487 = t498 * t519 + t509 * t543;
	t491 = t501 * t500 - t510 * t544;
	t522 = t487 * t499 + t491 * t496;
	t518 = t500 * t538 + t537;
	t483 = -qJD(1) * t492 + qJD(2) * t517;
	t516 = t483 * t500 + t497 * t528;
	t485 = qJD(1) * t494 - qJD(2) * t493;
	t515 = t485 * t500 + t497 * t529;
	t508 = cos(qJ(4));
	t514 = t503 * t524 - t508 * t548;
	t513 = t480 * t508 + t503 * t523;
	t488 = t498 * t518 + t504 * t543;
	t512 = t488 * t508 + t503 * t522;
	t507 = cos(qJ(5));
	t486 = -qJD(1) * t517 + qJD(2) * t492;
	t484 = -qJD(1) * t493 + qJD(2) * t494;
	t482 = -t485 * t497 + t500 * t529;
	t481 = -t483 * t497 + t500 * t528;
	t476 = t509 * t527 + (t519 * qJD(3) + (-t500 * t539 + t534) * qJD(2)) * t498;
	t475 = -t504 * t527 + (-t518 * qJD(3) + (-t500 * t537 - t538) * qJD(2)) * t498;
	t474 = -t475 * t496 + t499 * t526;
	t473 = qJD(3) * t477 + t486 * t509 + t504 * t515;
	t472 = t548 * qJD(3) - t486 * t504 + t515 * t509;
	t471 = qJD(3) * t479 + t484 * t509 + t504 * t516;
	t470 = -t480 * qJD(3) - t484 * t504 + t516 * t509;
	t469 = -t472 * t496 + t482 * t499;
	t468 = -t470 * t496 + t481 * t499;
	t1 = [0, t528, t481, t468, t471 * t503 + (-t470 * t499 - t481 * t496) * t508 + t513 * qJD(4), (t470 * t540 + t471 * t508 + t481 * t545) * t502 - t468 * t507 + (t513 * t507 + (-t479 * t496 + t490 * t499) * t502) * qJD(5) + (-t480 * t503 + t508 * t523) * t530; 0, t529, t482, t469, t473 * t503 + (-t472 * t499 - t482 * t496) * t508 + t514 * qJD(4), (t472 * t540 + t473 * t508 + t482 * t545) * t502 - t469 * t507 + (t514 * t507 + (-t477 * t496 + t489 * t499) * t502) * qJD(5) + (t503 * t548 + t508 * t524) * t530; 0, 0, t526, t474, t476 * t503 + (-t475 * t499 - t525) * t508 + t512 * qJD(4), (t475 * t540 + t476 * t508 + t503 * t525) * t502 - t474 * t507 + (t512 * t507 + (-t487 * t496 + t491 * t499) * t502) * qJD(5) + (-t488 * t503 + t508 * t522) * t530;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end