% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:45
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:45
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:45
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(10));
	t58 = sin(pkin(10));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:45
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->7), mult. (60->18), div. (0->0), fcn. (60->8), ass. (0->13)
	t105 = sin(pkin(11));
	t108 = cos(pkin(11));
	t111 = sin(qJ(2));
	t112 = cos(qJ(2));
	t115 = (t105 * t112 + t108 * t111) * qJD(2);
	t103 = (t105 * t111 - t108 * t112) * qJD(2);
	t110 = cos(pkin(6));
	t109 = cos(pkin(10));
	t107 = sin(pkin(6));
	t106 = sin(pkin(10));
	t102 = t110 * t115;
	t101 = t110 * t103;
	t1 = [0, t106 * t102 + t109 * t103, 0, 0, 0, 0; 0, -t109 * t102 + t106 * t103, 0, 0, 0, 0; 0, -t107 * t115, 0, 0, 0, 0; 0, -t106 * t101 + t109 * t115, 0, 0, 0, 0; 0, t109 * t101 + t106 * t115, 0, 0, 0, 0; 0, t107 * t103, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:46
	% EndTime: 2019-10-09 21:42:46
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (81->31), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->38)
	t294 = sin(pkin(6));
	t298 = sin(qJ(4));
	t310 = t294 * t298;
	t300 = cos(qJ(4));
	t309 = t294 * t300;
	t299 = sin(qJ(2));
	t308 = qJD(2) * t299;
	t301 = cos(qJ(2));
	t307 = qJD(2) * t301;
	t306 = qJD(4) * t298;
	t305 = qJD(4) * t300;
	t292 = sin(pkin(11));
	t295 = cos(pkin(11));
	t297 = cos(pkin(6));
	t279 = (t292 * t308 - t295 * t307) * t297;
	t286 = -t292 * t307 - t295 * t308;
	t293 = sin(pkin(10));
	t296 = cos(pkin(10));
	t270 = -t296 * t279 + t293 * t286;
	t272 = t293 * t279 + t296 * t286;
	t304 = t301 * t292 + t299 * t295;
	t303 = t299 * t292 - t301 * t295;
	t281 = t303 * t294;
	t302 = qJD(2) * t304;
	t285 = t303 * qJD(2);
	t284 = t304 * t297;
	t283 = t303 * t297;
	t282 = t304 * t294;
	t280 = t297 * t302;
	t278 = qJD(2) * t281;
	t277 = t294 * t302;
	t276 = -t293 * t284 - t296 * t303;
	t275 = t293 * t283 - t296 * t304;
	t274 = t296 * t284 - t293 * t303;
	t273 = -t296 * t283 - t293 * t304;
	t271 = t293 * t280 + t296 * t285;
	t269 = -t296 * t280 + t293 * t285;
	t1 = [0, t271 * t300 - t275 * t306, 0, -t272 * t298 + (-t276 * t300 - t293 * t310) * qJD(4), 0, 0; 0, t269 * t300 - t273 * t306, 0, -t270 * t298 + (-t274 * t300 + t296 * t310) * qJD(4), 0, 0; 0, -t277 * t300 + t281 * t306, 0, t278 * t298 + (-t282 * t300 - t297 * t298) * qJD(4), 0, 0; 0, -t271 * t298 - t275 * t305, 0, -t272 * t300 + (t276 * t298 - t293 * t309) * qJD(4), 0, 0; 0, -t269 * t298 - t273 * t305, 0, -t270 * t300 + (t274 * t298 + t296 * t309) * qJD(4), 0, 0; 0, t277 * t298 + t281 * t305, 0, t278 * t300 + (t282 * t298 - t297 * t300) * qJD(4), 0, 0; 0, t272, 0, 0, 0, 0; 0, t270, 0, 0, 0, 0; 0, -t278, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:48
	% EndTime: 2019-10-09 21:42:48
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (303->66), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->60)
	t463 = cos(pkin(6));
	t458 = sin(pkin(11));
	t461 = cos(pkin(11));
	t466 = sin(qJ(2));
	t469 = cos(qJ(2));
	t476 = t469 * t458 + t466 * t461;
	t448 = t476 * t463;
	t451 = t466 * t458 - t469 * t461;
	t449 = t451 * qJD(2);
	t460 = sin(pkin(6));
	t465 = sin(qJ(4));
	t494 = t460 * t465;
	t468 = cos(qJ(4));
	t493 = t460 * t468;
	t490 = qJD(2) * t466;
	t489 = qJD(2) * t469;
	t488 = qJD(4) * t465;
	t487 = qJD(4) * t468;
	t464 = sin(qJ(5));
	t486 = qJD(5) * t464;
	t467 = cos(qJ(5));
	t485 = qJD(5) * t467;
	t484 = qJD(5) * t468;
	t459 = sin(pkin(10));
	t462 = cos(pkin(10));
	t474 = t451 * t463;
	t434 = -t459 * t476 - t462 * t474;
	t445 = (t458 * t490 - t461 * t489) * t463;
	t450 = -t458 * t489 - t461 * t490;
	t480 = -t462 * t445 + t459 * t450;
	t483 = -t434 * t484 + t480;
	t437 = t459 * t474 - t462 * t476;
	t479 = t459 * t445 + t462 * t450;
	t482 = -t437 * t484 + t479;
	t444 = t460 * t449;
	t446 = t451 * t460;
	t481 = t446 * t484 - t444;
	t447 = t476 * t460;
	t440 = t447 * t468 + t463 * t465;
	t439 = -t447 * t465 + t463 * t468;
	t478 = t462 * t448 - t459 * t451;
	t477 = -t459 * t448 - t462 * t451;
	t423 = -t462 * t493 - t465 * t478;
	t475 = t462 * t494 - t468 * t478;
	t425 = t459 * t493 - t465 * t477;
	t426 = t459 * t494 + t468 * t477;
	t473 = qJD(2) * t448;
	t427 = t459 * t449 - t462 * t473;
	t472 = -qJD(5) * t478 - t427 * t468 + t434 * t488;
	t430 = t462 * t449 + t459 * t473;
	t471 = -qJD(5) * t477 - t430 * t468 + t437 * t488;
	t443 = qJD(2) * t447;
	t470 = qJD(5) * t447 - t443 * t468 + t446 * t488;
	t422 = t439 * qJD(4) - t444 * t468;
	t421 = -t440 * qJD(4) + t444 * t465;
	t420 = t425 * qJD(4) + t468 * t479;
	t419 = -t426 * qJD(4) - t465 * t479;
	t418 = t423 * qJD(4) + t468 * t480;
	t417 = t475 * qJD(4) - t465 * t480;
	t1 = [0, t482 * t464 - t471 * t467, 0, t419 * t467 - t425 * t486, -t420 * t464 - t430 * t467 + (-t426 * t467 + t437 * t464) * qJD(5), 0; 0, t483 * t464 - t472 * t467, 0, t417 * t467 - t423 * t486, -t418 * t464 - t427 * t467 + (t434 * t464 + t467 * t475) * qJD(5), 0; 0, t481 * t464 + t470 * t467, 0, t421 * t467 - t439 * t486, -t422 * t464 + t443 * t467 + (-t440 * t467 - t446 * t464) * qJD(5), 0; 0, t471 * t464 + t482 * t467, 0, -t419 * t464 - t425 * t485, -t420 * t467 + t430 * t464 + (t426 * t464 + t437 * t467) * qJD(5), 0; 0, t472 * t464 + t483 * t467, 0, -t417 * t464 - t423 * t485, -t418 * t467 + t427 * t464 + (t434 * t467 - t464 * t475) * qJD(5), 0; 0, -t470 * t464 + t481 * t467, 0, -t421 * t464 - t439 * t485, -t422 * t467 - t443 * t464 + (t440 * t464 - t446 * t467) * qJD(5), 0; 0, t430 * t465 + t437 * t487, 0, t420, 0, 0; 0, t427 * t465 + t434 * t487, 0, t418, 0, 0; 0, -t443 * t465 - t446 * t487, 0, t422, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:48
	% EndTime: 2019-10-09 21:42:48
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (303->66), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->60)
	t491 = cos(pkin(6));
	t486 = sin(pkin(11));
	t489 = cos(pkin(11));
	t494 = sin(qJ(2));
	t497 = cos(qJ(2));
	t504 = t497 * t486 + t494 * t489;
	t476 = t504 * t491;
	t479 = t494 * t486 - t497 * t489;
	t477 = t479 * qJD(2);
	t488 = sin(pkin(6));
	t493 = sin(qJ(4));
	t522 = t488 * t493;
	t496 = cos(qJ(4));
	t521 = t488 * t496;
	t518 = qJD(2) * t494;
	t517 = qJD(2) * t497;
	t516 = qJD(4) * t493;
	t515 = qJD(4) * t496;
	t492 = sin(qJ(5));
	t514 = qJD(5) * t492;
	t495 = cos(qJ(5));
	t513 = qJD(5) * t495;
	t512 = qJD(5) * t496;
	t487 = sin(pkin(10));
	t490 = cos(pkin(10));
	t502 = t479 * t491;
	t462 = -t487 * t504 - t490 * t502;
	t473 = (t486 * t518 - t489 * t517) * t491;
	t478 = -t486 * t517 - t489 * t518;
	t508 = -t490 * t473 + t487 * t478;
	t511 = -t462 * t512 + t508;
	t465 = t487 * t502 - t490 * t504;
	t507 = t487 * t473 + t490 * t478;
	t510 = -t465 * t512 + t507;
	t472 = t488 * t477;
	t474 = t479 * t488;
	t509 = t474 * t512 - t472;
	t475 = t504 * t488;
	t468 = t475 * t496 + t491 * t493;
	t467 = -t475 * t493 + t491 * t496;
	t506 = t490 * t476 - t487 * t479;
	t505 = -t487 * t476 - t490 * t479;
	t451 = -t490 * t521 - t493 * t506;
	t503 = t490 * t522 - t496 * t506;
	t453 = t487 * t521 - t493 * t505;
	t454 = t487 * t522 + t496 * t505;
	t501 = qJD(2) * t476;
	t455 = t487 * t477 - t490 * t501;
	t500 = -qJD(5) * t506 - t455 * t496 + t462 * t516;
	t458 = t490 * t477 + t487 * t501;
	t499 = -qJD(5) * t505 - t458 * t496 + t465 * t516;
	t471 = qJD(2) * t475;
	t498 = qJD(5) * t475 - t471 * t496 + t474 * t516;
	t450 = t467 * qJD(4) - t472 * t496;
	t449 = -t468 * qJD(4) + t472 * t493;
	t448 = t453 * qJD(4) + t496 * t507;
	t447 = -t454 * qJD(4) - t493 * t507;
	t446 = t451 * qJD(4) + t496 * t508;
	t445 = t503 * qJD(4) - t493 * t508;
	t1 = [0, t510 * t492 - t499 * t495, 0, t447 * t495 - t453 * t514, -t448 * t492 - t458 * t495 + (-t454 * t495 + t465 * t492) * qJD(5), 0; 0, t511 * t492 - t500 * t495, 0, t445 * t495 - t451 * t514, -t446 * t492 - t455 * t495 + (t462 * t492 + t495 * t503) * qJD(5), 0; 0, t509 * t492 + t498 * t495, 0, t449 * t495 - t467 * t514, -t450 * t492 + t471 * t495 + (-t468 * t495 - t474 * t492) * qJD(5), 0; 0, t499 * t492 + t510 * t495, 0, -t447 * t492 - t453 * t513, -t448 * t495 + t458 * t492 + (t454 * t492 + t465 * t495) * qJD(5), 0; 0, t500 * t492 + t511 * t495, 0, -t445 * t492 - t451 * t513, -t446 * t495 + t455 * t492 + (t462 * t495 - t492 * t503) * qJD(5), 0; 0, -t498 * t492 + t509 * t495, 0, -t449 * t492 - t467 * t513, -t450 * t495 - t471 * t492 + (t468 * t492 - t474 * t495) * qJD(5), 0; 0, t458 * t493 + t465 * t515, 0, t448, 0, 0; 0, t455 * t493 + t462 * t515, 0, t446, 0, 0; 0, -t471 * t493 - t474 * t515, 0, t450, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end