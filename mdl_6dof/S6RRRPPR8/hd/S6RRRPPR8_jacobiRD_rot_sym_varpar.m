% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:34
	% EndTime: 2019-10-10 11:29:34
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
	% StartTime: 2019-10-10 11:29:35
	% EndTime: 2019-10-10 11:29:35
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t350 = cos(pkin(6));
	t353 = sin(qJ(1));
	t352 = sin(qJ(2));
	t373 = t353 * t352;
	t364 = t350 * t373;
	t368 = qJD(2) * t352;
	t355 = cos(qJ(2));
	t356 = cos(qJ(1));
	t370 = t356 * t355;
	t340 = -qJD(1) * t364 - t353 * t368 + (qJD(2) * t350 + qJD(1)) * t370;
	t371 = t356 * t352;
	t372 = t353 * t355;
	t342 = t350 * t371 + t372;
	t351 = sin(qJ(3));
	t354 = cos(qJ(3));
	t349 = sin(pkin(6));
	t369 = qJD(1) * t349;
	t363 = t353 * t369;
	t374 = t349 * t356;
	t377 = (t342 * t351 + t354 * t374) * qJD(3) - t340 * t354 - t351 * t363;
	t376 = t349 * t351;
	t375 = t349 * t354;
	t367 = qJD(3) * t351;
	t366 = qJD(3) * t354;
	t365 = qJD(3) * t355;
	t362 = t356 * t369;
	t361 = t349 * qJD(2) * t355;
	t341 = t350 * t370 - t373;
	t343 = -t350 * t372 - t371;
	t359 = t364 - t370;
	t357 = -t340 * t351 - t342 * t366 + t354 * t363 + t367 * t374;
	t339 = t343 * qJD(1) - t342 * qJD(2);
	t338 = -t342 * qJD(1) + t343 * qJD(2);
	t337 = -t341 * qJD(1) + t359 * qJD(2);
	t336 = t351 * t362 + t338 * t354 + (t351 * t359 + t353 * t375) * qJD(3);
	t335 = -t354 * t362 + t338 * t351 + (t353 * t376 - t354 * t359) * qJD(3);
	t1 = [t377, t337 * t354 - t343 * t367, -t335, 0, 0, 0; t336, t339 * t354 - t341 * t367, t357, 0, 0, 0; 0, (-t351 * t365 - t354 * t368) * t349, -t351 * t361 + (-t350 * t351 - t352 * t375) * qJD(3), 0, 0, 0; t339, t338, 0, 0, 0, 0; -t337, t340, 0, 0, 0, 0; 0, t361, 0, 0, 0, 0; t357, t337 * t351 + t343 * t366, t336, 0, 0, 0; t335, t339 * t351 + t341 * t366, -t377, 0, 0, 0; 0, (-t351 * t368 + t354 * t365) * t349, t354 * t361 + (t350 * t354 - t352 * t376) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:34
	% EndTime: 2019-10-10 11:29:35
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (95->36), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->37)
	t316 = cos(pkin(6));
	t319 = sin(qJ(1));
	t318 = sin(qJ(2));
	t340 = t319 * t318;
	t331 = t316 * t340;
	t335 = qJD(2) * t318;
	t321 = cos(qJ(2));
	t322 = cos(qJ(1));
	t337 = t322 * t321;
	t307 = -qJD(1) * t331 - t319 * t335 + (qJD(2) * t316 + qJD(1)) * t337;
	t338 = t322 * t318;
	t339 = t319 * t321;
	t309 = t316 * t338 + t339;
	t317 = sin(qJ(3));
	t320 = cos(qJ(3));
	t315 = sin(pkin(6));
	t336 = qJD(1) * t315;
	t330 = t319 * t336;
	t341 = t315 * t322;
	t344 = (-t309 * t320 + t317 * t341) * qJD(3) - t307 * t317 + t320 * t330;
	t343 = t315 * t317;
	t342 = t315 * t320;
	t334 = qJD(3) * t317;
	t333 = qJD(3) * t320;
	t332 = qJD(3) * t321;
	t329 = t322 * t336;
	t328 = t315 * qJD(2) * t321;
	t308 = t316 * t337 - t340;
	t326 = t316 * t339 + t338;
	t325 = t331 - t337;
	t323 = t307 * t320 + t317 * t330 + (-t309 * t317 - t320 * t341) * qJD(3);
	t306 = t326 * qJD(1) + t309 * qJD(2);
	t305 = t309 * qJD(1) + t326 * qJD(2);
	t304 = -t308 * qJD(1) + t325 * qJD(2);
	t303 = t317 * t329 - t305 * t320 + (t317 * t325 + t319 * t342) * qJD(3);
	t302 = -t320 * t329 - t305 * t317 + (t319 * t343 - t320 * t325) * qJD(3);
	t1 = [t344, t304 * t317 - t326 * t333, t303, 0, 0, 0; t302, -t306 * t317 + t308 * t333, t323, 0, 0, 0; 0, (-t317 * t335 + t320 * t332) * t315, t320 * t328 + (t316 * t320 - t318 * t343) * qJD(3), 0, 0, 0; t323, -t304 * t320 - t326 * t334, t302, 0, 0, 0; -t303, t306 * t320 + t308 * t334, -t344, 0, 0, 0; 0, (t317 * t332 + t320 * t335) * t315, t317 * t328 + (t316 * t317 + t318 * t342) * qJD(3), 0, 0, 0; t306, t305, 0, 0, 0, 0; t304, -t307, 0, 0, 0, 0; 0, -t328, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:36
	% EndTime: 2019-10-10 11:29:37
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (289->74), mult. (898->138), div. (0->0), fcn. (964->10), ass. (0->63)
	t453 = cos(pkin(6));
	t457 = sin(qJ(1));
	t456 = sin(qJ(2));
	t490 = t457 * t456;
	t478 = t453 * t490;
	t485 = qJD(2) * t456;
	t460 = cos(qJ(2));
	t461 = cos(qJ(1));
	t487 = t461 * t460;
	t432 = -qJD(1) * t478 - t457 * t485 + (qJD(2) * t453 + qJD(1)) * t487;
	t488 = t461 * t456;
	t489 = t457 * t460;
	t443 = t453 * t488 + t489;
	t455 = sin(qJ(3));
	t459 = cos(qJ(3));
	t452 = sin(pkin(6));
	t486 = qJD(1) * t452;
	t476 = t457 * t486;
	t491 = t452 * t461;
	t477 = t455 * t491;
	t483 = qJD(3) * t459;
	t426 = -qJD(3) * t477 + t432 * t455 + t443 * t483 - t459 * t476;
	t466 = t453 * t489 + t488;
	t431 = t466 * qJD(1) + t443 * qJD(2);
	t435 = t443 * t455 + t459 * t491;
	t442 = -t453 * t487 + t490;
	t454 = sin(qJ(6));
	t458 = cos(qJ(6));
	t504 = (t435 * t458 - t442 * t454) * qJD(6) + t426 * t454 + t431 * t458;
	t501 = -t426 * t458 + t431 * t454 + (t435 * t454 + t442 * t458) * qJD(6);
	t482 = qJD(3) * t460;
	t500 = (qJD(2) * t455 + qJD(6)) * t456 - t459 * t482;
	t499 = t435 * qJD(3) - t432 * t459 - t455 * t476;
	t494 = t452 * t455;
	t493 = t452 * t459;
	t492 = t452 * t460;
	t484 = qJD(3) * t455;
	t481 = qJD(6) * t454;
	t480 = qJD(6) * t455;
	t479 = qJD(6) * t458;
	t475 = t461 * t486;
	t474 = t452 * t485;
	t473 = qJD(2) * t492;
	t430 = t443 * qJD(1) + t466 * qJD(2);
	t470 = t466 * t480 + t430;
	t469 = t442 * t480 - t432;
	t468 = (-qJD(2) - t480) * t460;
	t465 = t478 - t487;
	t467 = t455 * t465 + t457 * t493;
	t439 = t457 * t494 - t459 * t465;
	t441 = t453 * t455 + t456 * t493;
	t440 = -t453 * t459 + t456 * t494;
	t429 = t442 * qJD(1) + t465 * qJD(2);
	t463 = qJD(6) * t465 + t429 * t455 - t466 * t483;
	t462 = qJD(6) * t443 + t431 * t455 + t442 * t483;
	t436 = t443 * t459 - t477;
	t434 = -t440 * qJD(3) + t459 * t473;
	t433 = t441 * qJD(3) + t455 * t473;
	t425 = t467 * qJD(3) - t430 * t459 + t455 * t475;
	t424 = t439 * qJD(3) - t430 * t455 - t459 * t475;
	t423 = t424 * t458 + t429 * t454 + (t454 * t467 - t458 * t466) * qJD(6);
	t422 = -t424 * t454 + t429 * t458 + (t454 * t466 + t458 * t467) * qJD(6);
	t1 = [t501, t470 * t454 + t463 * t458, t425 * t458 - t439 * t481, 0, 0, t422; t423, t469 * t454 - t462 * t458, -t436 * t481 - t458 * t499, 0, 0, -t504; 0, (t454 * t468 - t500 * t458) * t452, t434 * t458 - t441 * t481, 0, 0, -t458 * t474 - t433 * t454 + (-t440 * t458 - t454 * t492) * qJD(6); t504, -t463 * t454 + t470 * t458, -t425 * t454 - t439 * t479, 0, 0, -t423; t422, t462 * t454 + t469 * t458, -t436 * t479 + t454 * t499, 0, 0, t501; 0, (t500 * t454 + t458 * t468) * t452, -t434 * t454 - t441 * t479, 0, 0, t454 * t474 - t433 * t458 + (t440 * t454 - t458 * t492) * qJD(6); t499, t429 * t459 + t466 * t484, -t424, 0, 0, 0; t425, -t431 * t459 + t442 * t484, -t426, 0, 0, 0; 0, (-t455 * t482 - t459 * t485) * t452, -t433, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end