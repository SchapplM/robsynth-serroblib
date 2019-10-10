% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
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
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
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
	% StartTime: 2019-10-09 21:33:41
	% EndTime: 2019-10-09 21:33:42
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
	% StartTime: 2019-10-09 21:33:42
	% EndTime: 2019-10-09 21:33:42
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (81->31), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->38)
	t348 = sin(pkin(6));
	t352 = sin(qJ(4));
	t364 = t348 * t352;
	t354 = cos(qJ(4));
	t363 = t348 * t354;
	t353 = sin(qJ(2));
	t362 = qJD(2) * t353;
	t355 = cos(qJ(2));
	t361 = qJD(2) * t355;
	t360 = qJD(4) * t352;
	t359 = qJD(4) * t354;
	t346 = sin(pkin(11));
	t349 = cos(pkin(11));
	t351 = cos(pkin(6));
	t333 = (t346 * t362 - t349 * t361) * t351;
	t340 = -t346 * t361 - t349 * t362;
	t347 = sin(pkin(10));
	t350 = cos(pkin(10));
	t324 = -t350 * t333 + t347 * t340;
	t326 = t347 * t333 + t350 * t340;
	t358 = t355 * t346 + t353 * t349;
	t357 = t353 * t346 - t355 * t349;
	t335 = t357 * t348;
	t356 = qJD(2) * t358;
	t339 = t357 * qJD(2);
	t338 = t358 * t351;
	t337 = t357 * t351;
	t336 = t358 * t348;
	t334 = t351 * t356;
	t332 = qJD(2) * t335;
	t331 = t348 * t356;
	t330 = -t347 * t338 - t350 * t357;
	t329 = t347 * t337 - t350 * t358;
	t328 = t350 * t338 - t347 * t357;
	t327 = -t350 * t337 - t347 * t358;
	t325 = t347 * t334 + t350 * t339;
	t323 = -t350 * t334 + t347 * t339;
	t1 = [0, t326, 0, 0, 0, 0; 0, t324, 0, 0, 0, 0; 0, -t332, 0, 0, 0, 0; 0, -t325 * t354 + t329 * t360, 0, t326 * t352 + (t330 * t354 + t347 * t364) * qJD(4), 0, 0; 0, -t323 * t354 + t327 * t360, 0, t324 * t352 + (t328 * t354 - t350 * t364) * qJD(4), 0, 0; 0, t331 * t354 - t335 * t360, 0, -t332 * t352 + (t336 * t354 + t351 * t352) * qJD(4), 0, 0; 0, t325 * t352 + t329 * t359, 0, t326 * t354 + (-t330 * t352 + t347 * t363) * qJD(4), 0, 0; 0, t323 * t352 + t327 * t359, 0, t324 * t354 + (-t328 * t352 - t350 * t363) * qJD(4), 0, 0; 0, -t331 * t352 - t335 * t359, 0, -t332 * t354 + (-t336 * t352 + t351 * t354) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:43
	% EndTime: 2019-10-09 21:33:44
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (303->69), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->60)
	t477 = cos(pkin(6));
	t472 = sin(pkin(11));
	t475 = cos(pkin(11));
	t480 = sin(qJ(2));
	t483 = cos(qJ(2));
	t490 = t483 * t472 + t480 * t475;
	t462 = t490 * t477;
	t465 = t480 * t472 - t483 * t475;
	t463 = t465 * qJD(2);
	t474 = sin(pkin(6));
	t479 = sin(qJ(4));
	t508 = t474 * t479;
	t482 = cos(qJ(4));
	t507 = t474 * t482;
	t504 = qJD(2) * t480;
	t503 = qJD(2) * t483;
	t502 = qJD(4) * t479;
	t501 = qJD(4) * t482;
	t478 = sin(qJ(6));
	t500 = qJD(6) * t478;
	t499 = qJD(6) * t479;
	t481 = cos(qJ(6));
	t498 = qJD(6) * t481;
	t473 = sin(pkin(10));
	t476 = cos(pkin(10));
	t488 = t465 * t477;
	t448 = -t473 * t490 - t476 * t488;
	t459 = (t472 * t504 - t475 * t503) * t477;
	t464 = -t472 * t503 - t475 * t504;
	t494 = -t476 * t459 + t473 * t464;
	t497 = t448 * t499 + t494;
	t451 = t473 * t488 - t476 * t490;
	t493 = t473 * t459 + t476 * t464;
	t496 = t451 * t499 + t493;
	t458 = t474 * t463;
	t460 = t465 * t474;
	t495 = t460 * t499 + t458;
	t461 = t490 * t474;
	t454 = t461 * t482 + t477 * t479;
	t453 = t461 * t479 - t477 * t482;
	t492 = t476 * t462 - t473 * t465;
	t491 = -t473 * t462 - t476 * t465;
	t437 = t476 * t507 + t479 * t492;
	t438 = -t476 * t508 + t482 * t492;
	t489 = t473 * t507 - t479 * t491;
	t440 = t473 * t508 + t482 * t491;
	t487 = qJD(2) * t462;
	t441 = t473 * t463 - t476 * t487;
	t486 = -qJD(6) * t492 + t441 * t479 + t448 * t501;
	t444 = t476 * t463 + t473 * t487;
	t485 = -qJD(6) * t491 + t444 * t479 + t451 * t501;
	t457 = qJD(2) * t461;
	t484 = -qJD(6) * t461 - t457 * t479 - t460 * t501;
	t436 = -t453 * qJD(4) - t458 * t482;
	t435 = t454 * qJD(4) - t458 * t479;
	t434 = t489 * qJD(4) + t482 * t493;
	t433 = t440 * qJD(4) + t479 * t493;
	t432 = -t437 * qJD(4) + t482 * t494;
	t431 = t438 * qJD(4) + t479 * t494;
	t1 = [0, t485 * t478 + t496 * t481, 0, t434 * t478 + t440 * t498, 0, t433 * t481 + t444 * t478 + (t451 * t481 + t478 * t489) * qJD(6); 0, t486 * t478 + t497 * t481, 0, t432 * t478 + t438 * t498, 0, t431 * t481 + t441 * t478 + (-t437 * t478 + t448 * t481) * qJD(6); 0, t484 * t478 - t495 * t481, 0, t436 * t478 + t454 * t498, 0, t435 * t481 - t457 * t478 + (-t453 * t478 - t460 * t481) * qJD(6); 0, -t496 * t478 + t485 * t481, 0, t434 * t481 - t440 * t500, 0, -t433 * t478 + t444 * t481 + (-t451 * t478 + t481 * t489) * qJD(6); 0, -t497 * t478 + t486 * t481, 0, t432 * t481 - t438 * t500, 0, -t431 * t478 + t441 * t481 + (-t437 * t481 - t448 * t478) * qJD(6); 0, t495 * t478 + t484 * t481, 0, t436 * t481 - t454 * t500, 0, -t435 * t478 - t457 * t481 + (-t453 * t481 + t460 * t478) * qJD(6); 0, t444 * t482 - t451 * t502, 0, -t433, 0, 0; 0, t441 * t482 - t448 * t502, 0, -t431, 0, 0; 0, -t457 * t482 + t460 * t502, 0, -t435, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end