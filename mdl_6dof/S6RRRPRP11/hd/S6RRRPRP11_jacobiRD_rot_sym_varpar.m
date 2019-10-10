% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:51
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRP11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:32
	% EndTime: 2019-10-10 11:51:32
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
	% StartTime: 2019-10-10 11:51:33
	% EndTime: 2019-10-10 11:51:33
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t284 = cos(pkin(6));
	t287 = sin(qJ(1));
	t286 = sin(qJ(2));
	t307 = t287 * t286;
	t298 = t284 * t307;
	t302 = qJD(2) * t286;
	t289 = cos(qJ(2));
	t290 = cos(qJ(1));
	t304 = t290 * t289;
	t275 = -qJD(1) * t298 - t287 * t302 + (qJD(2) * t284 + qJD(1)) * t304;
	t305 = t290 * t286;
	t306 = t287 * t289;
	t277 = t284 * t305 + t306;
	t285 = sin(qJ(3));
	t288 = cos(qJ(3));
	t283 = sin(pkin(6));
	t303 = qJD(1) * t283;
	t297 = t287 * t303;
	t308 = t283 * t290;
	t311 = (-t277 * t288 + t285 * t308) * qJD(3) - t275 * t285 + t288 * t297;
	t310 = t283 * t285;
	t309 = t283 * t288;
	t301 = qJD(3) * t285;
	t300 = qJD(3) * t288;
	t299 = qJD(3) * t289;
	t296 = t290 * t303;
	t295 = t283 * qJD(2) * t289;
	t276 = t284 * t304 - t307;
	t278 = -t284 * t306 - t305;
	t293 = t298 - t304;
	t291 = -t275 * t288 + t300 * t308 + (qJD(3) * t277 - t297) * t285;
	t274 = t278 * qJD(1) - t277 * qJD(2);
	t273 = -t277 * qJD(1) + t278 * qJD(2);
	t272 = -t276 * qJD(1) + t293 * qJD(2);
	t271 = t285 * t296 + t273 * t288 + (t285 * t293 + t287 * t309) * qJD(3);
	t270 = t288 * t296 - t273 * t285 + (-t287 * t310 + t288 * t293) * qJD(3);
	t1 = [t291, t272 * t288 - t278 * t301, t270, 0, 0, 0; t271, t274 * t288 - t276 * t301, t311, 0, 0, 0; 0, (-t285 * t299 - t288 * t302) * t283, -t285 * t295 + (-t284 * t285 - t286 * t309) * qJD(3), 0, 0, 0; -t311, -t272 * t285 - t278 * t300, -t271, 0, 0, 0; t270, -t274 * t285 - t276 * t300, t291, 0, 0, 0; 0, (t285 * t302 - t288 * t299) * t283, -t288 * t295 + (-t284 * t288 + t286 * t310) * qJD(3), 0, 0, 0; t274, t273, 0, 0, 0, 0; -t272, t275, 0, 0, 0, 0; 0, t295, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:34
	% EndTime: 2019-10-10 11:51:34
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (94->35), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->37)
	t337 = cos(pkin(6));
	t340 = sin(qJ(1));
	t339 = sin(qJ(2));
	t360 = t340 * t339;
	t351 = t337 * t360;
	t355 = qJD(2) * t339;
	t342 = cos(qJ(2));
	t343 = cos(qJ(1));
	t357 = t343 * t342;
	t328 = -qJD(1) * t351 - t340 * t355 + (qJD(2) * t337 + qJD(1)) * t357;
	t358 = t343 * t339;
	t359 = t340 * t342;
	t330 = t337 * t358 + t359;
	t338 = sin(qJ(3));
	t341 = cos(qJ(3));
	t336 = sin(pkin(6));
	t356 = qJD(1) * t336;
	t350 = t340 * t356;
	t361 = t336 * t343;
	t364 = (-t330 * t341 + t338 * t361) * qJD(3) - t328 * t338 + t341 * t350;
	t363 = t336 * t338;
	t362 = t336 * t341;
	t354 = qJD(3) * t338;
	t353 = qJD(3) * t341;
	t352 = qJD(3) * t342;
	t349 = t343 * t356;
	t348 = t336 * qJD(2) * t342;
	t329 = t337 * t357 - t360;
	t331 = -t337 * t359 - t358;
	t346 = t351 - t357;
	t344 = t328 * t341 + t338 * t350 + (-t330 * t338 - t341 * t361) * qJD(3);
	t327 = t331 * qJD(1) - t330 * qJD(2);
	t326 = -t330 * qJD(1) + t331 * qJD(2);
	t325 = -t329 * qJD(1) + t346 * qJD(2);
	t324 = t338 * t349 + t326 * t341 + (t338 * t346 + t340 * t362) * qJD(3);
	t323 = -t341 * t349 + t326 * t338 + (t340 * t363 - t341 * t346) * qJD(3);
	t1 = [t327, t326, 0, 0, 0, 0; -t325, t328, 0, 0, 0, 0; 0, t348, 0, 0, 0, 0; t344, -t325 * t341 + t331 * t354, t323, 0, 0, 0; -t324, -t327 * t341 + t329 * t354, -t364, 0, 0, 0; 0, (t338 * t352 + t341 * t355) * t336, t338 * t348 + (t337 * t338 + t339 * t362) * qJD(3), 0, 0, 0; t364, t325 * t338 + t331 * t353, t324, 0, 0, 0; t323, t327 * t338 + t329 * t353, t344, 0, 0, 0; 0, (-t338 * t355 + t341 * t352) * t336, t341 * t348 + (t337 * t341 - t339 * t363) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:35
	% EndTime: 2019-10-10 11:51:35
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (289->75), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t460 = sin(qJ(1));
	t456 = cos(pkin(6));
	t473 = qJD(2) * t456 + qJD(1);
	t459 = sin(qJ(2));
	t494 = t460 * t459;
	t480 = t456 * t494;
	t488 = qJD(2) * t459;
	t463 = cos(qJ(2));
	t464 = cos(qJ(1));
	t490 = t464 * t463;
	t434 = -qJD(1) * t480 - t460 * t488 + t473 * t490;
	t491 = t464 * t459;
	t493 = t460 * t463;
	t445 = t456 * t491 + t493;
	t458 = sin(qJ(3));
	t462 = cos(qJ(3));
	t455 = sin(pkin(6));
	t489 = qJD(1) * t455;
	t477 = t460 * t489;
	t496 = t455 * t464;
	t479 = t458 * t496;
	t485 = qJD(3) * t462;
	t428 = -qJD(3) * t479 + t434 * t458 + t445 * t485 - t462 * t477;
	t446 = t456 * t493 + t491;
	t433 = t446 * qJD(1) + t445 * qJD(2);
	t437 = t445 * t458 + t462 * t496;
	t478 = t456 * t490;
	t444 = -t478 + t494;
	t457 = sin(qJ(5));
	t461 = cos(qJ(5));
	t507 = (t437 * t457 + t444 * t461) * qJD(5) - t428 * t461 + t433 * t457;
	t504 = -t428 * t457 - t433 * t461 + (-t437 * t461 + t444 * t457) * qJD(5);
	t503 = t437 * qJD(3) - t434 * t462 - t458 * t477;
	t498 = t455 * t458;
	t497 = t455 * t462;
	t495 = t457 * t463;
	t492 = t461 * t463;
	t487 = qJD(2) * t463;
	t486 = qJD(3) * t458;
	t484 = qJD(3) * t463;
	t483 = qJD(5) * t457;
	t482 = qJD(5) * t458;
	t481 = qJD(5) * t461;
	t476 = t464 * t489;
	t475 = t455 * t488;
	t474 = t455 * t487;
	t472 = qJD(2) + t482;
	t432 = -t445 * qJD(1) - t446 * qJD(2);
	t471 = t446 * t482 - t432;
	t470 = t444 * t482 - t434;
	t447 = -t480 + t490;
	t469 = -t447 * t458 + t460 * t497;
	t441 = t447 * t462 + t460 * t498;
	t443 = t456 * t458 + t459 * t497;
	t442 = -t456 * t462 + t459 * t498;
	t431 = -qJD(1) * t478 - t464 * t487 + t473 * t494;
	t467 = -qJD(5) * t447 + t431 * t458 - t446 * t485;
	t466 = -qJD(5) * t445 - t433 * t458 - t444 * t485;
	t465 = t462 * t484 + (-qJD(2) * t458 - qJD(5)) * t459;
	t438 = t445 * t462 - t479;
	t436 = -t442 * qJD(3) + t462 * t474;
	t435 = t443 * qJD(3) + t458 * t474;
	t427 = t469 * qJD(3) + t432 * t462 + t458 * t476;
	t426 = t441 * qJD(3) + t432 * t458 - t462 * t476;
	t425 = t426 * t457 - t431 * t461 + (-t446 * t457 - t461 * t469) * qJD(5);
	t424 = t426 * t461 + t431 * t457 + (-t446 * t461 + t457 * t469) * qJD(5);
	t1 = [t504, t457 * t467 - t461 * t471, t427 * t457 + t441 * t481, 0, t424, 0; t425, t457 * t466 - t461 * t470, t438 * t481 - t457 * t503, 0, -t507, 0; 0, (t457 * t465 + t472 * t492) * t455, t436 * t457 + t443 * t481, 0, -t457 * t475 + t435 * t461 + (-t442 * t457 + t455 * t492) * qJD(5), 0; t507, t457 * t471 + t461 * t467, t427 * t461 - t441 * t483, 0, -t425, 0; t424, t457 * t470 + t461 * t466, -t438 * t483 - t461 * t503, 0, t504, 0; 0, (t461 * t465 - t472 * t495) * t455, t436 * t461 - t443 * t483, 0, -t461 * t475 - t435 * t457 + (-t442 * t461 - t455 * t495) * qJD(5), 0; t503, t431 * t462 + t446 * t486, -t426, 0, 0, 0; t427, -t433 * t462 + t444 * t486, -t428, 0, 0, 0; 0, (-t458 * t484 - t462 * t488) * t455, -t435, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:35
	% EndTime: 2019-10-10 11:51:36
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (289->75), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t479 = sin(qJ(1));
	t475 = cos(pkin(6));
	t492 = qJD(2) * t475 + qJD(1);
	t478 = sin(qJ(2));
	t513 = t478 * t479;
	t498 = t475 * t513;
	t507 = qJD(2) * t478;
	t482 = cos(qJ(2));
	t483 = cos(qJ(1));
	t509 = t482 * t483;
	t453 = -qJD(1) * t498 - t479 * t507 + t492 * t509;
	t511 = t479 * t482;
	t512 = t478 * t483;
	t464 = t475 * t512 + t511;
	t477 = sin(qJ(3));
	t481 = cos(qJ(3));
	t474 = sin(pkin(6));
	t508 = qJD(1) * t474;
	t496 = t479 * t508;
	t515 = t474 * t483;
	t499 = t477 * t515;
	t504 = qJD(3) * t481;
	t447 = -qJD(3) * t499 + t453 * t477 + t464 * t504 - t481 * t496;
	t465 = t475 * t511 + t512;
	t452 = t465 * qJD(1) + t464 * qJD(2);
	t456 = t464 * t477 + t481 * t515;
	t497 = t475 * t509;
	t463 = -t497 + t513;
	t476 = sin(qJ(5));
	t480 = cos(qJ(5));
	t526 = (t456 * t476 + t463 * t480) * qJD(5) - t447 * t480 + t452 * t476;
	t523 = -t447 * t476 - t452 * t480 + (-t456 * t480 + t463 * t476) * qJD(5);
	t522 = t456 * qJD(3) - t453 * t481 - t477 * t496;
	t517 = t474 * t477;
	t516 = t474 * t481;
	t514 = t476 * t482;
	t510 = t480 * t482;
	t506 = qJD(2) * t482;
	t505 = qJD(3) * t477;
	t503 = qJD(3) * t482;
	t502 = qJD(5) * t476;
	t501 = qJD(5) * t477;
	t500 = qJD(5) * t480;
	t495 = t483 * t508;
	t494 = t474 * t507;
	t493 = t474 * t506;
	t491 = qJD(2) + t501;
	t451 = -t464 * qJD(1) - t465 * qJD(2);
	t490 = t465 * t501 - t451;
	t489 = t463 * t501 - t453;
	t466 = -t498 + t509;
	t488 = -t466 * t477 + t479 * t516;
	t460 = t466 * t481 + t479 * t517;
	t462 = t475 * t477 + t478 * t516;
	t461 = -t475 * t481 + t478 * t517;
	t450 = -qJD(1) * t497 - t483 * t506 + t492 * t513;
	t486 = -qJD(5) * t466 + t450 * t477 - t465 * t504;
	t485 = -qJD(5) * t464 - t452 * t477 - t463 * t504;
	t484 = t481 * t503 + (-qJD(2) * t477 - qJD(5)) * t478;
	t457 = t464 * t481 - t499;
	t455 = -t461 * qJD(3) + t481 * t493;
	t454 = t462 * qJD(3) + t477 * t493;
	t446 = t488 * qJD(3) + t451 * t481 + t477 * t495;
	t445 = t460 * qJD(3) + t451 * t477 - t481 * t495;
	t444 = t445 * t476 - t450 * t480 + (-t465 * t476 - t480 * t488) * qJD(5);
	t443 = t445 * t480 + t450 * t476 + (-t465 * t480 + t476 * t488) * qJD(5);
	t1 = [t523, t486 * t476 - t490 * t480, t446 * t476 + t460 * t500, 0, t443, 0; t444, t485 * t476 - t489 * t480, t457 * t500 - t476 * t522, 0, -t526, 0; 0, (t484 * t476 + t491 * t510) * t474, t455 * t476 + t462 * t500, 0, -t476 * t494 + t454 * t480 + (-t461 * t476 + t474 * t510) * qJD(5), 0; t526, t490 * t476 + t486 * t480, t446 * t480 - t460 * t502, 0, -t444, 0; t443, t489 * t476 + t485 * t480, -t457 * t502 - t480 * t522, 0, t523, 0; 0, (t484 * t480 - t491 * t514) * t474, t455 * t480 - t462 * t502, 0, -t480 * t494 - t454 * t476 + (-t461 * t480 - t474 * t514) * qJD(5), 0; t522, t450 * t481 + t465 * t505, -t445, 0, 0, 0; t446, -t452 * t481 + t463 * t505, -t447, 0, 0, 0; 0, (-t477 * t503 - t481 * t507) * t474, -t454, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end