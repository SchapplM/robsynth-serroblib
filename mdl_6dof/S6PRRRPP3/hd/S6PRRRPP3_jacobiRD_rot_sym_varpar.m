% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
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
	% StartTime: 2019-10-09 22:44:51
	% EndTime: 2019-10-09 22:44:51
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
	% StartTime: 2019-10-09 22:44:52
	% EndTime: 2019-10-09 22:44:52
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->53)
	t354 = sin(qJ(3));
	t355 = sin(qJ(2));
	t357 = cos(qJ(3));
	t358 = cos(qJ(2));
	t374 = qJD(3) * t358;
	t384 = (qJD(2) * t357 - qJD(4)) * t355 + t354 * t374;
	t350 = sin(pkin(6));
	t383 = t350 * t354;
	t382 = t350 * t357;
	t381 = t350 * t358;
	t352 = cos(pkin(6));
	t380 = t352 * t355;
	t379 = t352 * t358;
	t378 = qJD(2) * t355;
	t377 = qJD(2) * t358;
	t376 = qJD(3) * t354;
	t375 = qJD(3) * t357;
	t353 = sin(qJ(4));
	t373 = qJD(4) * t353;
	t356 = cos(qJ(4));
	t372 = qJD(4) * t356;
	t371 = qJD(4) * t357;
	t349 = sin(pkin(10));
	t370 = t349 * t380;
	t369 = t350 * t378;
	t368 = t350 * t377;
	t351 = cos(pkin(10));
	t361 = -t349 * t355 + t351 * t379;
	t337 = t361 * qJD(2);
	t365 = -t361 * t371 + t337;
	t343 = t349 * t379 + t351 * t355;
	t339 = t343 * qJD(2);
	t364 = t343 * t371 - t339;
	t363 = (qJD(2) - t371) * t358;
	t342 = t349 * t358 + t351 * t380;
	t331 = -t342 * t354 - t351 * t382;
	t362 = -t342 * t357 + t351 * t383;
	t344 = t351 * t358 - t370;
	t333 = -t344 * t354 + t349 * t382;
	t334 = t344 * t357 + t349 * t383;
	t346 = t352 * t354 + t355 * t382;
	t345 = t352 * t357 - t355 * t383;
	t338 = t342 * qJD(2);
	t360 = qJD(4) * t342 - t338 * t357 - t361 * t376;
	t340 = -qJD(2) * t370 + t351 * t377;
	t359 = qJD(4) * t344 - t340 * t357 + t343 * t376;
	t336 = qJD(3) * t345 + t357 * t368;
	t335 = -qJD(3) * t346 - t354 * t368;
	t330 = qJD(3) * t333 - t339 * t357;
	t329 = -qJD(3) * t334 + t339 * t354;
	t328 = qJD(3) * t331 + t337 * t357;
	t327 = qJD(3) * t362 - t337 * t354;
	t1 = [0, t353 * t364 + t356 * t359, t329 * t356 - t333 * t373, -t330 * t353 + t340 * t356 + (-t334 * t356 - t343 * t353) * qJD(4), 0, 0; 0, t353 * t365 + t356 * t360, t327 * t356 - t331 * t373, -t328 * t353 + t338 * t356 + (t353 * t361 + t356 * t362) * qJD(4), 0, 0; 0, (t353 * t363 - t356 * t384) * t350, t335 * t356 - t345 * t373, t356 * t369 - t336 * t353 + (-t346 * t356 + t353 * t381) * qJD(4), 0, 0; 0, -t353 * t359 + t356 * t364, -t329 * t353 - t333 * t372, -t330 * t356 - t340 * t353 + (t334 * t353 - t343 * t356) * qJD(4), 0, 0; 0, -t353 * t360 + t356 * t365, -t327 * t353 - t331 * t372, -t328 * t356 - t338 * t353 + (-t353 * t362 + t356 * t361) * qJD(4), 0, 0; 0, (t353 * t384 + t356 * t363) * t350, -t335 * t353 - t345 * t372, -t353 * t369 - t336 * t356 + (t346 * t353 + t356 * t381) * qJD(4), 0, 0; 0, -t340 * t354 - t343 * t375, t330, 0, 0, 0; 0, -t338 * t354 + t361 * t375, t328, 0, 0, 0; 0, (-t354 * t378 + t357 * t374) * t350, t336, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:53
	% EndTime: 2019-10-09 22:44:54
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->53)
	t424 = sin(qJ(3));
	t425 = sin(qJ(2));
	t427 = cos(qJ(3));
	t428 = cos(qJ(2));
	t444 = qJD(3) * t428;
	t454 = (qJD(2) * t427 - qJD(4)) * t425 + t424 * t444;
	t420 = sin(pkin(6));
	t453 = t420 * t424;
	t452 = t420 * t427;
	t451 = t420 * t428;
	t422 = cos(pkin(6));
	t450 = t422 * t425;
	t449 = t422 * t428;
	t448 = qJD(2) * t425;
	t447 = qJD(2) * t428;
	t446 = qJD(3) * t424;
	t445 = qJD(3) * t427;
	t423 = sin(qJ(4));
	t443 = qJD(4) * t423;
	t426 = cos(qJ(4));
	t442 = qJD(4) * t426;
	t441 = qJD(4) * t427;
	t419 = sin(pkin(10));
	t440 = t419 * t450;
	t439 = t420 * t448;
	t438 = t420 * t447;
	t421 = cos(pkin(10));
	t431 = -t419 * t425 + t421 * t449;
	t407 = t431 * qJD(2);
	t435 = t431 * t441 - t407;
	t413 = t419 * t449 + t421 * t425;
	t409 = t413 * qJD(2);
	t434 = -t413 * t441 + t409;
	t433 = (-qJD(2) + t441) * t428;
	t412 = t419 * t428 + t421 * t450;
	t401 = -t412 * t424 - t421 * t452;
	t432 = -t412 * t427 + t421 * t453;
	t414 = t421 * t428 - t440;
	t403 = -t414 * t424 + t419 * t452;
	t404 = t414 * t427 + t419 * t453;
	t416 = t422 * t424 + t425 * t452;
	t415 = t422 * t427 - t425 * t453;
	t408 = t412 * qJD(2);
	t430 = qJD(4) * t412 - t408 * t427 - t431 * t446;
	t410 = -qJD(2) * t440 + t421 * t447;
	t429 = qJD(4) * t414 - t410 * t427 + t413 * t446;
	t406 = t415 * qJD(3) + t427 * t438;
	t405 = -t416 * qJD(3) - t424 * t438;
	t400 = t403 * qJD(3) - t409 * t427;
	t399 = -t404 * qJD(3) + t409 * t424;
	t398 = t401 * qJD(3) + t407 * t427;
	t397 = t432 * qJD(3) - t407 * t424;
	t1 = [0, -t410 * t424 - t413 * t445, t400, 0, 0, 0; 0, -t408 * t424 + t431 * t445, t398, 0, 0, 0; 0, (-t424 * t448 + t427 * t444) * t420, t406, 0, 0, 0; 0, t434 * t423 - t429 * t426, -t399 * t426 + t403 * t443, t400 * t423 - t410 * t426 + (t404 * t426 + t413 * t423) * qJD(4), 0, 0; 0, t435 * t423 - t430 * t426, -t397 * t426 + t401 * t443, t398 * t423 - t408 * t426 + (-t423 * t431 - t426 * t432) * qJD(4), 0, 0; 0, (t423 * t433 + t454 * t426) * t420, -t405 * t426 + t415 * t443, -t426 * t439 + t406 * t423 + (t416 * t426 - t423 * t451) * qJD(4), 0, 0; 0, t429 * t423 + t434 * t426, t399 * t423 + t403 * t442, t400 * t426 + t410 * t423 + (-t404 * t423 + t413 * t426) * qJD(4), 0, 0; 0, t430 * t423 + t435 * t426, t397 * t423 + t401 * t442, t398 * t426 + t408 * t423 + (t423 * t432 - t426 * t431) * qJD(4), 0, 0; 0, (-t454 * t423 + t426 * t433) * t420, t405 * t423 + t415 * t442, t423 * t439 + t406 * t426 + (-t416 * t423 - t426 * t451) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:54
	% EndTime: 2019-10-09 22:44:54
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->54)
	t426 = sin(pkin(6));
	t430 = sin(qJ(3));
	t459 = t426 * t430;
	t433 = cos(qJ(3));
	t458 = t426 * t433;
	t428 = cos(pkin(6));
	t431 = sin(qJ(2));
	t457 = t428 * t431;
	t434 = cos(qJ(2));
	t456 = t428 * t434;
	t429 = sin(qJ(4));
	t455 = t429 * t434;
	t432 = cos(qJ(4));
	t454 = t432 * t434;
	t453 = qJD(2) * t431;
	t452 = qJD(2) * t434;
	t451 = qJD(3) * t430;
	t450 = qJD(3) * t433;
	t449 = qJD(3) * t434;
	t448 = qJD(4) * t429;
	t447 = qJD(4) * t432;
	t446 = qJD(4) * t433;
	t425 = sin(pkin(10));
	t445 = t425 * t457;
	t444 = t426 * t453;
	t443 = t426 * t452;
	t442 = -qJD(2) + t446;
	t427 = cos(pkin(10));
	t438 = -t425 * t431 + t427 * t456;
	t413 = t438 * qJD(2);
	t441 = -t438 * t446 + t413;
	t419 = t425 * t456 + t427 * t431;
	t415 = t419 * qJD(2);
	t440 = t419 * t446 - t415;
	t418 = t425 * t434 + t427 * t457;
	t407 = -t418 * t430 - t427 * t458;
	t439 = -t418 * t433 + t427 * t459;
	t420 = t427 * t434 - t445;
	t409 = -t420 * t430 + t425 * t458;
	t410 = t420 * t433 + t425 * t459;
	t422 = t428 * t430 + t431 * t458;
	t421 = t428 * t433 - t431 * t459;
	t414 = t418 * qJD(2);
	t437 = qJD(4) * t418 - t414 * t433 - t438 * t451;
	t416 = -qJD(2) * t445 + t427 * t452;
	t436 = qJD(4) * t420 - t416 * t433 + t419 * t451;
	t435 = -t430 * t449 + (-qJD(2) * t433 + qJD(4)) * t431;
	t412 = t421 * qJD(3) + t433 * t443;
	t411 = -t422 * qJD(3) - t430 * t443;
	t406 = t409 * qJD(3) - t415 * t433;
	t405 = -t410 * qJD(3) + t415 * t430;
	t404 = t407 * qJD(3) + t413 * t433;
	t403 = t439 * qJD(3) - t413 * t430;
	t1 = [0, -t416 * t430 - t419 * t450, t406, 0, 0, 0; 0, -t414 * t430 + t438 * t450, t404, 0, 0, 0; 0, (-t430 * t453 + t433 * t449) * t426, t412, 0, 0, 0; 0, t436 * t429 - t440 * t432, t405 * t429 + t409 * t447, t406 * t432 + t416 * t429 + (-t410 * t429 + t419 * t432) * qJD(4), 0, 0; 0, t437 * t429 - t441 * t432, t403 * t429 + t407 * t447, t404 * t432 + t414 * t429 + (t429 * t439 - t432 * t438) * qJD(4), 0, 0; 0, (t435 * t429 + t442 * t454) * t426, t411 * t429 + t421 * t447, t429 * t444 + t412 * t432 + (-t422 * t429 - t426 * t454) * qJD(4), 0, 0; 0, t440 * t429 + t436 * t432, t405 * t432 - t409 * t448, -t406 * t429 + t416 * t432 + (-t410 * t432 - t419 * t429) * qJD(4), 0, 0; 0, t441 * t429 + t437 * t432, t403 * t432 - t407 * t448, -t404 * t429 + t414 * t432 + (t429 * t438 + t432 * t439) * qJD(4), 0, 0; 0, (t435 * t432 - t442 * t455) * t426, t411 * t432 - t421 * t448, t432 * t444 - t412 * t429 + (-t422 * t432 + t426 * t455) * qJD(4), 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end