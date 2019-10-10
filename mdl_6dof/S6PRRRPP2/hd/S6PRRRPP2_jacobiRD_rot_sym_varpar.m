% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP2
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
% Datum: 2019-10-09 22:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
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
	% StartTime: 2019-10-09 22:43:01
	% EndTime: 2019-10-09 22:43:01
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
	% StartTime: 2019-10-09 22:43:03
	% EndTime: 2019-10-09 22:43:03
	% DurationCPUTime: 0.37s
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
	% StartTime: 2019-10-09 22:43:04
	% EndTime: 2019-10-09 22:43:04
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->54)
	t431 = sin(pkin(6));
	t435 = sin(qJ(3));
	t464 = t431 * t435;
	t438 = cos(qJ(3));
	t463 = t431 * t438;
	t433 = cos(pkin(6));
	t436 = sin(qJ(2));
	t462 = t433 * t436;
	t439 = cos(qJ(2));
	t461 = t433 * t439;
	t434 = sin(qJ(4));
	t460 = t434 * t439;
	t437 = cos(qJ(4));
	t459 = t437 * t439;
	t458 = qJD(2) * t436;
	t457 = qJD(2) * t439;
	t456 = qJD(3) * t435;
	t455 = qJD(3) * t438;
	t454 = qJD(3) * t439;
	t453 = qJD(4) * t434;
	t452 = qJD(4) * t437;
	t451 = qJD(4) * t438;
	t430 = sin(pkin(10));
	t450 = t430 * t462;
	t449 = t431 * t458;
	t448 = t431 * t457;
	t447 = -qJD(2) + t451;
	t432 = cos(pkin(10));
	t443 = -t430 * t436 + t432 * t461;
	t418 = t443 * qJD(2);
	t446 = -t443 * t451 + t418;
	t424 = t430 * t461 + t432 * t436;
	t420 = t424 * qJD(2);
	t445 = t424 * t451 - t420;
	t423 = t430 * t439 + t432 * t462;
	t412 = -t423 * t435 - t432 * t463;
	t444 = -t423 * t438 + t432 * t464;
	t425 = t432 * t439 - t450;
	t414 = -t425 * t435 + t430 * t463;
	t415 = t425 * t438 + t430 * t464;
	t427 = t433 * t435 + t436 * t463;
	t426 = t433 * t438 - t436 * t464;
	t419 = t423 * qJD(2);
	t442 = qJD(4) * t423 - t419 * t438 - t443 * t456;
	t421 = -qJD(2) * t450 + t432 * t457;
	t441 = qJD(4) * t425 - t421 * t438 + t424 * t456;
	t440 = -t435 * t454 + (-qJD(2) * t438 + qJD(4)) * t436;
	t417 = t426 * qJD(3) + t438 * t448;
	t416 = -t427 * qJD(3) - t435 * t448;
	t411 = t414 * qJD(3) - t420 * t438;
	t410 = -t415 * qJD(3) + t420 * t435;
	t409 = t412 * qJD(3) + t418 * t438;
	t408 = t444 * qJD(3) - t418 * t435;
	t1 = [0, t445 * t434 + t441 * t437, t410 * t437 - t414 * t453, -t411 * t434 + t421 * t437 + (-t415 * t437 - t424 * t434) * qJD(4), 0, 0; 0, t446 * t434 + t442 * t437, t408 * t437 - t412 * t453, -t409 * t434 + t419 * t437 + (t434 * t443 + t437 * t444) * qJD(4), 0, 0; 0, (t440 * t437 - t447 * t460) * t431, t416 * t437 - t426 * t453, t437 * t449 - t417 * t434 + (-t427 * t437 + t431 * t460) * qJD(4), 0, 0; 0, -t421 * t435 - t424 * t455, t411, 0, 0, 0; 0, -t419 * t435 + t443 * t455, t409, 0, 0, 0; 0, (-t435 * t458 + t438 * t454) * t431, t417, 0, 0, 0; 0, t441 * t434 - t445 * t437, t410 * t434 + t414 * t452, t411 * t437 + t421 * t434 + (-t415 * t434 + t424 * t437) * qJD(4), 0, 0; 0, t442 * t434 - t446 * t437, t408 * t434 + t412 * t452, t409 * t437 + t419 * t434 + (t434 * t444 - t437 * t443) * qJD(4), 0, 0; 0, (t440 * t434 + t447 * t459) * t431, t416 * t434 + t426 * t452, t434 * t449 + t417 * t437 + (-t427 * t434 - t431 * t459) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:03
	% EndTime: 2019-10-09 22:43:03
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (153->62), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->54)
	t401 = sin(pkin(6));
	t405 = sin(qJ(3));
	t434 = t401 * t405;
	t408 = cos(qJ(3));
	t433 = t401 * t408;
	t403 = cos(pkin(6));
	t406 = sin(qJ(2));
	t432 = t403 * t406;
	t409 = cos(qJ(2));
	t431 = t403 * t409;
	t404 = sin(qJ(4));
	t430 = t404 * t409;
	t407 = cos(qJ(4));
	t429 = t407 * t409;
	t428 = qJD(2) * t406;
	t427 = qJD(2) * t409;
	t426 = qJD(3) * t405;
	t425 = qJD(3) * t408;
	t424 = qJD(3) * t409;
	t423 = qJD(4) * t404;
	t422 = qJD(4) * t407;
	t421 = qJD(4) * t408;
	t400 = sin(pkin(10));
	t420 = t400 * t432;
	t419 = t401 * t428;
	t418 = t401 * t427;
	t417 = -qJD(2) + t421;
	t402 = cos(pkin(10));
	t413 = -t400 * t406 + t402 * t431;
	t388 = t413 * qJD(2);
	t416 = -t413 * t421 + t388;
	t394 = t400 * t431 + t402 * t406;
	t390 = t394 * qJD(2);
	t415 = t394 * t421 - t390;
	t393 = t400 * t409 + t402 * t432;
	t382 = -t393 * t405 - t402 * t433;
	t414 = -t393 * t408 + t402 * t434;
	t395 = t402 * t409 - t420;
	t384 = -t395 * t405 + t400 * t433;
	t385 = t395 * t408 + t400 * t434;
	t397 = t403 * t405 + t406 * t433;
	t396 = t403 * t408 - t406 * t434;
	t389 = t393 * qJD(2);
	t412 = qJD(4) * t393 - t389 * t408 - t413 * t426;
	t391 = -qJD(2) * t420 + t402 * t427;
	t411 = qJD(4) * t395 - t391 * t408 + t394 * t426;
	t410 = -t405 * t424 + (-qJD(2) * t408 + qJD(4)) * t406;
	t387 = t396 * qJD(3) + t408 * t418;
	t386 = -t397 * qJD(3) - t405 * t418;
	t381 = t384 * qJD(3) - t390 * t408;
	t380 = -t385 * qJD(3) + t390 * t405;
	t379 = t382 * qJD(3) + t388 * t408;
	t378 = t414 * qJD(3) - t388 * t405;
	t1 = [0, t415 * t404 + t411 * t407, t380 * t407 - t384 * t423, -t381 * t404 + t391 * t407 + (-t385 * t407 - t394 * t404) * qJD(4), 0, 0; 0, t416 * t404 + t412 * t407, t378 * t407 - t382 * t423, -t379 * t404 + t389 * t407 + (t404 * t413 + t407 * t414) * qJD(4), 0, 0; 0, (t410 * t407 - t417 * t430) * t401, t386 * t407 - t396 * t423, t407 * t419 - t387 * t404 + (-t397 * t407 + t401 * t430) * qJD(4), 0, 0; 0, t411 * t404 - t415 * t407, t380 * t404 + t384 * t422, t381 * t407 + t391 * t404 + (-t385 * t404 + t394 * t407) * qJD(4), 0, 0; 0, t412 * t404 - t416 * t407, t378 * t404 + t382 * t422, t379 * t407 + t389 * t404 + (t404 * t414 - t407 * t413) * qJD(4), 0, 0; 0, (t410 * t404 + t417 * t429) * t401, t386 * t404 + t396 * t422, t404 * t419 + t387 * t407 + (-t397 * t404 - t401 * t429) * qJD(4), 0, 0; 0, t391 * t405 + t394 * t425, -t381, 0, 0, 0; 0, t389 * t405 - t413 * t425, -t379, 0, 0, 0; 0, (t405 * t428 - t408 * t424) * t401, -t387, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end