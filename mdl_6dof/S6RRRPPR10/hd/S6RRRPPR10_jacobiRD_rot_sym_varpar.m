% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
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
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
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
	% StartTime: 2019-10-10 11:33:21
	% EndTime: 2019-10-10 11:33:21
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
	% StartTime: 2019-10-10 11:33:22
	% EndTime: 2019-10-10 11:33:22
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
	% StartTime: 2019-10-10 11:33:22
	% EndTime: 2019-10-10 11:33:23
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (157->47), mult. (518->103), div. (0->0), fcn. (536->10), ass. (0->46)
	t391 = sin(qJ(1));
	t388 = cos(pkin(6));
	t399 = qJD(2) * t388 + qJD(1);
	t390 = sin(qJ(2));
	t415 = t391 * t390;
	t405 = t388 * t415;
	t410 = qJD(2) * t390;
	t393 = cos(qJ(2));
	t394 = cos(qJ(1));
	t412 = t394 * t393;
	t373 = -qJD(1) * t405 - t391 * t410 + t399 * t412;
	t413 = t394 * t390;
	t414 = t391 * t393;
	t376 = t388 * t413 + t414;
	t389 = sin(qJ(3));
	t392 = cos(qJ(3));
	t386 = sin(pkin(6));
	t411 = qJD(1) * t386;
	t403 = t391 * t411;
	t417 = t386 * t394;
	t419 = (t376 * t389 + t392 * t417) * qJD(3) - t373 * t392 - t389 * t403;
	t418 = t386 * t392;
	t416 = t389 * t390;
	t409 = qJD(2) * t393;
	t408 = qJD(3) * t389;
	t407 = qJD(3) * t392;
	t406 = qJD(3) * t393;
	t404 = t388 * t412;
	t402 = t394 * t411;
	t401 = t386 * t409;
	t400 = t392 * t406;
	t377 = -t388 * t414 - t413;
	t370 = -qJD(1) * t404 - t394 * t409 + t399 * t415;
	t397 = t370 * t389 + t377 * t407;
	t372 = t377 * qJD(1) - t376 * qJD(2);
	t375 = t404 - t415;
	t396 = t372 * t389 + t375 * t407;
	t369 = -t373 * t389 - t376 * t407 + t392 * t403 + t408 * t417;
	t387 = cos(pkin(11));
	t385 = sin(pkin(11));
	t378 = -t405 + t412;
	t374 = t392 * t401 + (-t386 * t416 + t388 * t392) * qJD(3);
	t371 = -t376 * qJD(1) + t377 * qJD(2);
	t367 = t389 * t402 + t371 * t392 + (-t378 * t389 + t391 * t418) * qJD(3);
	t366 = -t392 * t402 + t371 * t389 + (t386 * t389 * t391 + t378 * t392) * qJD(3);
	t1 = [t369 * t385 + t372 * t387, t371 * t387 + t397 * t385, t367 * t385, 0, 0, 0; t366 * t385 - t370 * t387, t373 * t387 + t396 * t385, -t419 * t385, 0, 0, 0; 0, (t385 * t400 + (-t385 * t416 + t387 * t393) * qJD(2)) * t386, t374 * t385, 0, 0, 0; t369 * t387 - t372 * t385, -t371 * t385 + t397 * t387, t367 * t387, 0, 0, 0; t366 * t387 + t370 * t385, -t373 * t385 + t396 * t387, -t419 * t387, 0, 0, 0; 0, (t387 * t400 + (-t385 * t393 - t387 * t416) * qJD(2)) * t386, t374 * t387, 0, 0, 0; t419, t370 * t392 - t377 * t408, -t366, 0, 0, 0; t367, t372 * t392 - t375 * t408, t369, 0, 0, 0; 0, (-t389 * t406 - t392 * t410) * t386, -t389 * t401 + (-t388 * t389 - t390 * t418) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:23
	% EndTime: 2019-10-10 11:33:23
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (371->76), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t474 = sin(qJ(1));
	t471 = cos(pkin(6));
	t486 = qJD(2) * t471 + qJD(1);
	t473 = sin(qJ(2));
	t506 = t474 * t473;
	t493 = t471 * t506;
	t501 = qJD(2) * t473;
	t476 = cos(qJ(2));
	t477 = cos(qJ(1));
	t503 = t477 * t476;
	t446 = -qJD(1) * t493 - t474 * t501 + t486 * t503;
	t504 = t477 * t473;
	t505 = t474 * t476;
	t457 = t471 * t504 + t505;
	t472 = sin(qJ(3));
	t475 = cos(qJ(3));
	t470 = sin(pkin(6));
	t502 = qJD(1) * t470;
	t490 = t474 * t502;
	t507 = t470 * t477;
	t492 = t472 * t507;
	t498 = qJD(3) * t475;
	t440 = -qJD(3) * t492 + t446 * t472 + t457 * t498 - t475 * t490;
	t458 = t471 * t505 + t504;
	t445 = t458 * qJD(1) + t457 * qJD(2);
	t449 = t457 * t472 + t475 * t507;
	t491 = t471 * t503;
	t456 = -t491 + t506;
	t469 = pkin(11) + qJ(6);
	t467 = sin(t469);
	t468 = cos(t469);
	t520 = (t449 * t467 + t456 * t468) * qJD(6) - t440 * t468 + t445 * t467;
	t517 = -t440 * t467 - t445 * t468 + (-t449 * t468 + t456 * t467) * qJD(6);
	t494 = qJD(6) * t472;
	t516 = t476 * (qJD(2) + t494);
	t515 = t449 * qJD(3) - t446 * t475 - t472 * t490;
	t510 = t470 * t472;
	t509 = t470 * t475;
	t508 = t470 * t476;
	t500 = qJD(2) * t476;
	t499 = qJD(3) * t472;
	t497 = qJD(3) * t476;
	t496 = qJD(6) * t467;
	t495 = qJD(6) * t468;
	t489 = t477 * t502;
	t488 = t470 * t501;
	t487 = t470 * t500;
	t444 = -t457 * qJD(1) - t458 * qJD(2);
	t484 = t458 * t494 - t444;
	t483 = t456 * t494 - t446;
	t459 = -t493 + t503;
	t482 = -t459 * t472 + t474 * t509;
	t453 = t459 * t475 + t474 * t510;
	t455 = t471 * t472 + t473 * t509;
	t454 = -t471 * t475 + t473 * t510;
	t443 = -qJD(1) * t491 - t477 * t500 + t486 * t506;
	t480 = -qJD(6) * t459 + t443 * t472 - t458 * t498;
	t479 = -qJD(6) * t457 - t445 * t472 - t456 * t498;
	t478 = t475 * t497 + (-qJD(2) * t472 - qJD(6)) * t473;
	t450 = t457 * t475 - t492;
	t448 = -t454 * qJD(3) + t475 * t487;
	t447 = t455 * qJD(3) + t472 * t487;
	t439 = t482 * qJD(3) + t444 * t475 + t472 * t489;
	t438 = t453 * qJD(3) + t444 * t472 - t475 * t489;
	t437 = t438 * t467 - t443 * t468 + (-t458 * t467 - t468 * t482) * qJD(6);
	t436 = t438 * t468 + t443 * t467 + (-t458 * t468 + t467 * t482) * qJD(6);
	t1 = [t517, t467 * t480 - t468 * t484, t439 * t467 + t453 * t495, 0, 0, t436; t437, t467 * t479 - t468 * t483, t450 * t495 - t467 * t515, 0, 0, -t520; 0, (t467 * t478 + t468 * t516) * t470, t448 * t467 + t455 * t495, 0, 0, -t467 * t488 + t447 * t468 + (-t454 * t467 + t468 * t508) * qJD(6); t520, t467 * t484 + t468 * t480, t439 * t468 - t453 * t496, 0, 0, -t437; t436, t467 * t483 + t468 * t479, -t450 * t496 - t468 * t515, 0, 0, t517; 0, (-t467 * t516 + t468 * t478) * t470, t448 * t468 - t455 * t496, 0, 0, -t468 * t488 - t447 * t467 + (-t454 * t468 - t467 * t508) * qJD(6); t515, t443 * t475 + t458 * t499, -t438, 0, 0, 0; t439, -t445 * t475 + t456 * t499, -t440, 0, 0, 0; 0, (-t472 * t497 - t475 * t501) * t470, -t447, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end