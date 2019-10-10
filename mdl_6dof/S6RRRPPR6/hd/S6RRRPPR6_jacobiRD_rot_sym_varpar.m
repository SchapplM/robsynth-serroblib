% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:57
	% EndTime: 2019-10-10 11:25:57
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
	% StartTime: 2019-10-10 11:25:57
	% EndTime: 2019-10-10 11:25:58
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (144->36), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->38)
	t310 = cos(pkin(6));
	t312 = sin(qJ(1));
	t311 = sin(qJ(2));
	t331 = t312 * t311;
	t322 = t310 * t331;
	t326 = qJD(2) * t311;
	t313 = cos(qJ(2));
	t314 = cos(qJ(1));
	t328 = t314 * t313;
	t298 = -qJD(1) * t322 - t312 * t326 + (qJD(2) * t310 + qJD(1)) * t328;
	t329 = t314 * t311;
	t330 = t312 * t313;
	t300 = t310 * t329 + t330;
	t308 = qJ(3) + pkin(11);
	t306 = sin(t308);
	t307 = cos(t308);
	t309 = sin(pkin(6));
	t327 = qJD(1) * t309;
	t321 = t312 * t327;
	t332 = t309 * t314;
	t335 = (-t300 * t307 + t306 * t332) * qJD(3) - t298 * t306 + t307 * t321;
	t334 = t309 * t311;
	t333 = t309 * t312;
	t325 = qJD(3) * t306;
	t324 = qJD(3) * t307;
	t323 = qJD(3) * t313;
	t320 = t314 * t327;
	t319 = t309 * qJD(2) * t313;
	t299 = t310 * t328 - t331;
	t301 = -t310 * t330 - t329;
	t317 = t322 - t328;
	t315 = -t298 * t307 + t324 * t332 + (qJD(3) * t300 - t321) * t306;
	t297 = t301 * qJD(1) - t300 * qJD(2);
	t296 = -t300 * qJD(1) + t301 * qJD(2);
	t295 = -t299 * qJD(1) + t317 * qJD(2);
	t294 = t306 * t320 + t296 * t307 + (t306 * t317 + t307 * t333) * qJD(3);
	t293 = t307 * t320 - t296 * t306 + (-t306 * t333 + t307 * t317) * qJD(3);
	t1 = [t315, t295 * t307 - t301 * t325, t293, 0, 0, 0; t294, t297 * t307 - t299 * t325, t335, 0, 0, 0; 0, (-t306 * t323 - t307 * t326) * t309, -t306 * t319 + (-t306 * t310 - t307 * t334) * qJD(3), 0, 0, 0; -t335, -t295 * t306 - t301 * t324, -t294, 0, 0, 0; t293, -t297 * t306 - t299 * t324, t315, 0, 0, 0; 0, (t306 * t326 - t307 * t323) * t309, -t307 * t319 + (t306 * t334 - t307 * t310) * qJD(3), 0, 0, 0; t297, t296, 0, 0, 0, 0; -t295, t298, 0, 0, 0, 0; 0, t319, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:58
	% EndTime: 2019-10-10 11:25:58
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (144->36), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->38)
	t365 = cos(pkin(6));
	t367 = sin(qJ(1));
	t366 = sin(qJ(2));
	t386 = t367 * t366;
	t377 = t365 * t386;
	t381 = qJD(2) * t366;
	t368 = cos(qJ(2));
	t369 = cos(qJ(1));
	t383 = t369 * t368;
	t353 = -qJD(1) * t377 - t367 * t381 + (qJD(2) * t365 + qJD(1)) * t383;
	t384 = t369 * t366;
	t385 = t367 * t368;
	t355 = t365 * t384 + t385;
	t363 = qJ(3) + pkin(11);
	t361 = sin(t363);
	t362 = cos(t363);
	t364 = sin(pkin(6));
	t382 = qJD(1) * t364;
	t376 = t367 * t382;
	t387 = t364 * t369;
	t390 = (-t355 * t362 + t361 * t387) * qJD(3) - t353 * t361 + t362 * t376;
	t389 = t364 * t366;
	t388 = t364 * t367;
	t380 = qJD(3) * t361;
	t379 = qJD(3) * t362;
	t378 = qJD(3) * t368;
	t375 = t369 * t382;
	t374 = t364 * qJD(2) * t368;
	t354 = t365 * t383 - t386;
	t356 = -t365 * t385 - t384;
	t372 = t377 - t383;
	t370 = t353 * t362 + t361 * t376 + (-t355 * t361 - t362 * t387) * qJD(3);
	t352 = t356 * qJD(1) - t355 * qJD(2);
	t351 = -t355 * qJD(1) + t356 * qJD(2);
	t350 = -t354 * qJD(1) + t372 * qJD(2);
	t349 = t361 * t375 + t351 * t362 + (t361 * t372 + t362 * t388) * qJD(3);
	t348 = -t362 * t375 + t351 * t361 + (t361 * t388 - t362 * t372) * qJD(3);
	t1 = [t352, t351, 0, 0, 0, 0; -t350, t353, 0, 0, 0, 0; 0, t374, 0, 0, 0, 0; t370, -t350 * t362 + t356 * t380, t348, 0, 0, 0; -t349, -t352 * t362 + t354 * t380, -t390, 0, 0, 0; 0, (t361 * t378 + t362 * t381) * t364, t361 * t374 + (t361 * t365 + t362 * t389) * qJD(3), 0, 0, 0; t390, t350 * t361 + t356 * t379, t349, 0, 0, 0; t348, t352 * t361 + t354 * t379, t370, 0, 0, 0; 0, (-t361 * t381 + t362 * t378) * t364, t362 * t374 + (-t361 * t389 + t362 * t365) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:59
	% EndTime: 2019-10-10 11:26:00
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (424->76), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->68)
	t482 = sin(qJ(1));
	t479 = cos(pkin(6));
	t494 = qJD(2) * t479 + qJD(1);
	t481 = sin(qJ(2));
	t515 = t482 * t481;
	t501 = t479 * t515;
	t509 = qJD(2) * t481;
	t484 = cos(qJ(2));
	t485 = cos(qJ(1));
	t511 = t485 * t484;
	t454 = -qJD(1) * t501 - t482 * t509 + t494 * t511;
	t512 = t485 * t481;
	t514 = t482 * t484;
	t465 = t479 * t512 + t514;
	t477 = qJ(3) + pkin(11);
	t475 = sin(t477);
	t476 = cos(t477);
	t478 = sin(pkin(6));
	t510 = qJD(1) * t478;
	t498 = t482 * t510;
	t517 = t478 * t485;
	t500 = t475 * t517;
	t506 = qJD(3) * t476;
	t448 = -qJD(3) * t500 + t454 * t475 + t465 * t506 - t476 * t498;
	t466 = t479 * t514 + t512;
	t453 = t466 * qJD(1) + t465 * qJD(2);
	t457 = t465 * t475 + t476 * t517;
	t499 = t479 * t511;
	t464 = -t499 + t515;
	t480 = sin(qJ(6));
	t483 = cos(qJ(6));
	t528 = (t457 * t480 + t464 * t483) * qJD(6) - t448 * t483 + t453 * t480;
	t525 = -t448 * t480 - t453 * t483 + (-t457 * t483 + t464 * t480) * qJD(6);
	t524 = t457 * qJD(3) - t454 * t476 - t475 * t498;
	t519 = t478 * t481;
	t518 = t478 * t482;
	t516 = t480 * t484;
	t513 = t483 * t484;
	t508 = qJD(2) * t484;
	t507 = qJD(3) * t475;
	t505 = qJD(3) * t484;
	t504 = qJD(6) * t475;
	t503 = qJD(6) * t480;
	t502 = qJD(6) * t483;
	t497 = t485 * t510;
	t496 = t478 * t509;
	t495 = t478 * t508;
	t493 = qJD(2) + t504;
	t452 = -t465 * qJD(1) - t466 * qJD(2);
	t492 = t466 * t504 - t452;
	t491 = t464 * t504 - t454;
	t467 = -t501 + t511;
	t490 = -t467 * t475 + t476 * t518;
	t461 = t467 * t476 + t475 * t518;
	t463 = t479 * t475 + t476 * t519;
	t462 = t475 * t519 - t479 * t476;
	t451 = -qJD(1) * t499 - t485 * t508 + t494 * t515;
	t488 = -qJD(6) * t467 + t451 * t475 - t466 * t506;
	t487 = -qJD(6) * t465 - t453 * t475 - t464 * t506;
	t486 = t476 * t505 + (-qJD(2) * t475 - qJD(6)) * t481;
	t458 = t465 * t476 - t500;
	t456 = -t462 * qJD(3) + t476 * t495;
	t455 = t463 * qJD(3) + t475 * t495;
	t447 = t490 * qJD(3) + t452 * t476 + t475 * t497;
	t446 = t461 * qJD(3) + t452 * t475 - t476 * t497;
	t445 = t446 * t480 - t451 * t483 + (-t466 * t480 - t483 * t490) * qJD(6);
	t444 = t446 * t483 + t451 * t480 + (-t466 * t483 + t480 * t490) * qJD(6);
	t1 = [t525, t488 * t480 - t492 * t483, t447 * t480 + t461 * t502, 0, 0, t444; t445, t487 * t480 - t491 * t483, t458 * t502 - t480 * t524, 0, 0, -t528; 0, (t486 * t480 + t493 * t513) * t478, t456 * t480 + t463 * t502, 0, 0, -t480 * t496 + t455 * t483 + (-t462 * t480 + t478 * t513) * qJD(6); t528, t492 * t480 + t488 * t483, t447 * t483 - t461 * t503, 0, 0, -t445; t444, t491 * t480 + t487 * t483, -t458 * t503 - t483 * t524, 0, 0, t525; 0, (t486 * t483 - t493 * t516) * t478, t456 * t483 - t463 * t503, 0, 0, -t483 * t496 - t455 * t480 + (-t462 * t483 - t478 * t516) * qJD(6); t524, t451 * t476 + t466 * t507, -t446, 0, 0, 0; t447, -t453 * t476 + t464 * t507, -t448, 0, 0, 0; 0, (-t475 * t505 - t476 * t509) * t478, -t455, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end