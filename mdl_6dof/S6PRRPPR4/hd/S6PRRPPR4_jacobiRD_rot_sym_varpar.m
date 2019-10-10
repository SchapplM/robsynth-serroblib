% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
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
	% StartTime: 2019-10-09 22:12:45
	% EndTime: 2019-10-09 22:12:45
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(6));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(10));
	t212 = cos(pkin(10));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0, 0; 0, t204, 0, 0, 0, 0; 0, t202, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:46
	% EndTime: 2019-10-09 22:12:46
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (64->32), mult. (240->82), div. (0->0), fcn. (252->10), ass. (0->35)
	t299 = sin(pkin(6));
	t303 = sin(qJ(3));
	t320 = t299 * t303;
	t305 = cos(qJ(3));
	t319 = t299 * t305;
	t302 = cos(pkin(6));
	t304 = sin(qJ(2));
	t318 = t302 * t304;
	t306 = cos(qJ(2));
	t317 = t302 * t306;
	t316 = t303 * t304;
	t315 = t304 * t305;
	t314 = qJD(3) * t303;
	t313 = qJD(3) * t305;
	t312 = qJD(3) * t306;
	t311 = qJD(2) * t299 * t306;
	t310 = t303 * t312;
	t298 = sin(pkin(10));
	t301 = cos(pkin(10));
	t293 = -t298 * t304 + t301 * t317;
	t294 = t298 * t306 + t301 * t318;
	t295 = -t298 * t317 - t301 * t304;
	t309 = t298 * t318 - t301 * t306;
	t290 = t294 * qJD(2);
	t308 = t290 * t305 + t293 * t314;
	t292 = t309 * qJD(2);
	t307 = -t292 * t305 + t295 * t314;
	t300 = cos(pkin(11));
	t297 = sin(pkin(11));
	t291 = t295 * qJD(2);
	t289 = t293 * qJD(2);
	t288 = -t303 * t311 + (-t299 * t315 - t302 * t303) * qJD(3);
	t287 = -t291 * t303 + (-t298 * t320 + t305 * t309) * qJD(3);
	t286 = -t289 * t303 + (-t294 * t305 + t301 * t320) * qJD(3);
	t1 = [0, t291 * t297 - t307 * t300, t287 * t300, 0, 0, 0; 0, t289 * t297 - t308 * t300, t286 * t300, 0, 0, 0; 0, (-t300 * t310 + (t297 * t306 - t300 * t315) * qJD(2)) * t299, t288 * t300, 0, 0, 0; 0, t291 * t300 + t307 * t297, -t287 * t297, 0, 0, 0; 0, t289 * t300 + t308 * t297, -t286 * t297, 0, 0, 0; 0, (t297 * t310 + (t297 * t315 + t300 * t306) * qJD(2)) * t299, -t288 * t297, 0, 0, 0; 0, t292 * t303 + t295 * t313, t291 * t305 + (t298 * t319 + t303 * t309) * qJD(3), 0, 0, 0; 0, -t290 * t303 + t293 * t313, t289 * t305 + (-t294 * t303 - t301 * t319) * qJD(3), 0, 0, 0; 0, (-qJD(2) * t316 + t305 * t312) * t299, t305 * t311 + (-t299 * t316 + t302 * t305) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:46
	% EndTime: 2019-10-09 22:12:46
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (61->29), mult. (240->82), div. (0->0), fcn. (252->10), ass. (0->35)
	t339 = sin(pkin(6));
	t343 = sin(qJ(3));
	t360 = t339 * t343;
	t345 = cos(qJ(3));
	t359 = t339 * t345;
	t342 = cos(pkin(6));
	t344 = sin(qJ(2));
	t358 = t342 * t344;
	t346 = cos(qJ(2));
	t357 = t342 * t346;
	t356 = t343 * t344;
	t355 = t344 * t345;
	t354 = qJD(3) * t343;
	t353 = qJD(3) * t345;
	t352 = qJD(3) * t346;
	t351 = qJD(2) * t339 * t346;
	t350 = t343 * t352;
	t338 = sin(pkin(10));
	t341 = cos(pkin(10));
	t333 = -t338 * t344 + t341 * t357;
	t334 = t338 * t346 + t341 * t358;
	t335 = -t338 * t357 - t341 * t344;
	t349 = t338 * t358 - t341 * t346;
	t330 = t334 * qJD(2);
	t348 = -t330 * t345 - t333 * t354;
	t332 = t349 * qJD(2);
	t347 = t332 * t345 - t335 * t354;
	t340 = cos(pkin(11));
	t337 = sin(pkin(11));
	t331 = t335 * qJD(2);
	t329 = t333 * qJD(2);
	t328 = -t343 * t351 + (-t339 * t355 - t342 * t343) * qJD(3);
	t327 = -t331 * t343 + (-t338 * t360 + t345 * t349) * qJD(3);
	t326 = -t329 * t343 + (-t334 * t345 + t341 * t360) * qJD(3);
	t1 = [0, t331 * t337 + t340 * t347, t327 * t340, 0, 0, 0; 0, t329 * t337 + t340 * t348, t326 * t340, 0, 0, 0; 0, (-t340 * t350 + (t337 * t346 - t340 * t355) * qJD(2)) * t339, t328 * t340, 0, 0, 0; 0, t332 * t343 + t335 * t353, t331 * t345 + (t338 * t359 + t343 * t349) * qJD(3), 0, 0, 0; 0, -t330 * t343 + t333 * t353, t329 * t345 + (-t334 * t343 - t341 * t359) * qJD(3), 0, 0, 0; 0, (-qJD(2) * t356 + t345 * t352) * t339, t345 * t351 + (-t339 * t356 + t342 * t345) * qJD(3), 0, 0, 0; 0, -t331 * t340 + t337 * t347, t327 * t337, 0, 0, 0; 0, -t329 * t340 + t337 * t348, t326 * t337, 0, 0, 0; 0, (-t337 * t350 + (-t337 * t355 - t340 * t346) * qJD(2)) * t339, t328 * t337, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:47
	% EndTime: 2019-10-09 22:12:47
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (307->95), mult. (1016->198), div. (0->0), fcn. (1132->12), ass. (0->82)
	t462 = sin(pkin(11));
	t472 = cos(qJ(3));
	t500 = t462 * t472;
	t473 = cos(qJ(2));
	t499 = t462 * t473;
	t464 = sin(pkin(6));
	t469 = sin(qJ(3));
	t498 = t464 * t469;
	t497 = t464 * t472;
	t465 = cos(pkin(11));
	t496 = t465 * t472;
	t495 = t465 * t473;
	t467 = cos(pkin(6));
	t470 = sin(qJ(2));
	t494 = t467 * t470;
	t493 = t467 * t473;
	t492 = t470 * t472;
	t491 = t472 * t473;
	t490 = qJD(2) * t470;
	t489 = qJD(2) * t473;
	t488 = qJD(3) * t469;
	t487 = qJD(3) * t472;
	t486 = qJD(3) * t473;
	t463 = sin(pkin(10));
	t485 = t463 * t494;
	t484 = t464 * t490;
	t483 = t464 * t489;
	t482 = t469 * t486;
	t468 = sin(qJ(6));
	t471 = cos(qJ(6));
	t481 = t462 * t471 - t465 * t468;
	t480 = t462 * t468 + t465 * t471;
	t466 = cos(pkin(10));
	t455 = t463 * t473 + t466 * t494;
	t442 = -t455 * t469 - t466 * t497;
	t479 = -t455 * t472 + t466 * t498;
	t457 = t466 * t473 - t485;
	t444 = -t457 * t469 + t463 * t497;
	t445 = t457 * t472 + t463 * t498;
	t478 = -t463 * t470 + t466 * t493;
	t456 = t463 * t493 + t466 * t470;
	t459 = t464 * t492 + t467 * t469;
	t458 = t467 * t472 - t470 * t498;
	t451 = t455 * qJD(2);
	t477 = -t451 * t472 - t478 * t488;
	t453 = -qJD(2) * t485 + t466 * t489;
	t476 = -t453 * t472 + t456 * t488;
	t475 = qJD(6) * t481;
	t474 = qJD(6) * t480;
	t452 = t456 * qJD(2);
	t450 = t478 * qJD(2);
	t449 = (t462 * t470 + t465 * t491) * t464;
	t448 = (t462 * t491 - t465 * t470) * t464;
	t447 = t458 * qJD(3) + t472 * t483;
	t446 = -t459 * qJD(3) - t469 * t483;
	t441 = t459 * t465 - t464 * t499;
	t440 = t459 * t462 + t464 * t495;
	t439 = (-t465 * t482 + (-t465 * t492 + t499) * qJD(2)) * t464;
	t438 = (-t462 * t482 + (-t462 * t492 - t495) * qJD(2)) * t464;
	t437 = -t456 * t496 + t457 * t462;
	t436 = -t456 * t500 - t457 * t465;
	t435 = t455 * t462 + t478 * t496;
	t434 = -t455 * t465 + t478 * t500;
	t433 = t447 * t465 + t462 * t484;
	t432 = t447 * t462 - t465 * t484;
	t431 = t445 * t465 + t456 * t462;
	t430 = t445 * t462 - t456 * t465;
	t429 = -t462 * t478 - t465 * t479;
	t428 = -t462 * t479 + t465 * t478;
	t427 = t444 * qJD(3) - t452 * t472;
	t426 = -t445 * qJD(3) + t452 * t469;
	t425 = t442 * qJD(3) + t450 * t472;
	t424 = t479 * qJD(3) - t450 * t469;
	t423 = -t452 * t462 + t476 * t465;
	t422 = t452 * t465 + t476 * t462;
	t421 = t450 * t462 + t477 * t465;
	t420 = -t450 * t465 + t477 * t462;
	t419 = t427 * t465 + t453 * t462;
	t418 = t427 * t462 - t453 * t465;
	t417 = t425 * t465 + t451 * t462;
	t416 = t425 * t462 - t451 * t465;
	t1 = [0, t422 * t468 + t423 * t471 + (t436 * t471 - t437 * t468) * qJD(6), t480 * t426 + t444 * t475, 0, 0, t418 * t471 - t419 * t468 + (-t430 * t468 - t431 * t471) * qJD(6); 0, t420 * t468 + t421 * t471 + (t434 * t471 - t435 * t468) * qJD(6), t480 * t424 + t442 * t475, 0, 0, t416 * t471 - t417 * t468 + (-t428 * t468 - t429 * t471) * qJD(6); 0, t438 * t468 + t439 * t471 + (t448 * t471 - t449 * t468) * qJD(6), t480 * t446 + t458 * t475, 0, 0, t432 * t471 - t433 * t468 + (-t440 * t468 - t441 * t471) * qJD(6); 0, t422 * t471 - t423 * t468 + (-t436 * t468 - t437 * t471) * qJD(6), t481 * t426 - t444 * t474, 0, 0, -t418 * t468 - t419 * t471 + (-t430 * t471 + t431 * t468) * qJD(6); 0, t420 * t471 - t421 * t468 + (-t434 * t468 - t435 * t471) * qJD(6), t481 * t424 - t442 * t474, 0, 0, -t416 * t468 - t417 * t471 + (-t428 * t471 + t429 * t468) * qJD(6); 0, t438 * t471 - t439 * t468 + (-t448 * t468 - t449 * t471) * qJD(6), t481 * t446 - t458 * t474, 0, 0, -t432 * t468 - t433 * t471 + (-t440 * t471 + t441 * t468) * qJD(6); 0, t453 * t469 + t456 * t487, -t427, 0, 0, 0; 0, t451 * t469 - t478 * t487, -t425, 0, 0, 0; 0, (t469 * t490 - t472 * t486) * t464, -t447, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end