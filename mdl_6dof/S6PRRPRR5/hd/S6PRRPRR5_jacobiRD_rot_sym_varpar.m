% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:26
	% EndTime: 2019-10-09 22:33:26
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
	t210 = sin(pkin(11));
	t212 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:33:27
	% EndTime: 2019-10-09 22:33:27
	% DurationCPUTime: 0.23s
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
	t298 = sin(pkin(11));
	t301 = cos(pkin(11));
	t293 = -t298 * t304 + t301 * t317;
	t294 = t298 * t306 + t301 * t318;
	t295 = -t298 * t317 - t301 * t304;
	t309 = t298 * t318 - t301 * t306;
	t290 = t294 * qJD(2);
	t308 = t290 * t305 + t293 * t314;
	t292 = t309 * qJD(2);
	t307 = -t292 * t305 + t295 * t314;
	t300 = cos(pkin(12));
	t297 = sin(pkin(12));
	t291 = t295 * qJD(2);
	t289 = t293 * qJD(2);
	t288 = -t303 * t311 + (-t299 * t315 - t302 * t303) * qJD(3);
	t287 = -t291 * t303 + (-t298 * t320 + t305 * t309) * qJD(3);
	t286 = -t289 * t303 + (-t294 * t305 + t301 * t320) * qJD(3);
	t1 = [0, t291 * t297 - t300 * t307, t287 * t300, 0, 0, 0; 0, t289 * t297 - t300 * t308, t286 * t300, 0, 0, 0; 0, (-t300 * t310 + (t297 * t306 - t300 * t315) * qJD(2)) * t299, t288 * t300, 0, 0, 0; 0, t291 * t300 + t297 * t307, -t287 * t297, 0, 0, 0; 0, t289 * t300 + t297 * t308, -t286 * t297, 0, 0, 0; 0, (t297 * t310 + (t297 * t315 + t300 * t306) * qJD(2)) * t299, -t288 * t297, 0, 0, 0; 0, t292 * t303 + t295 * t313, t291 * t305 + (t298 * t319 + t303 * t309) * qJD(3), 0, 0, 0; 0, -t290 * t303 + t293 * t313, t289 * t305 + (-t294 * t303 - t301 * t319) * qJD(3), 0, 0, 0; 0, (-qJD(2) * t316 + t305 * t312) * t299, t305 * t311 + (-t299 * t316 + t302 * t305) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:27
	% EndTime: 2019-10-09 22:33:28
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (219->60), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->54)
	t369 = sin(qJ(3));
	t370 = sin(qJ(2));
	t371 = cos(qJ(3));
	t372 = cos(qJ(2));
	t388 = qJD(3) * t372;
	t398 = (qJD(2) * t371 - qJD(5)) * t370 + t369 * t388;
	t366 = sin(pkin(6));
	t397 = t366 * t369;
	t396 = t366 * t371;
	t395 = t366 * t372;
	t368 = cos(pkin(6));
	t394 = t368 * t370;
	t393 = t368 * t372;
	t392 = qJD(2) * t370;
	t391 = qJD(2) * t372;
	t390 = qJD(3) * t369;
	t389 = qJD(3) * t371;
	t364 = pkin(12) + qJ(5);
	t362 = sin(t364);
	t387 = qJD(5) * t362;
	t363 = cos(t364);
	t386 = qJD(5) * t363;
	t385 = qJD(5) * t371;
	t365 = sin(pkin(11));
	t384 = t365 * t394;
	t383 = t366 * t392;
	t382 = t366 * t391;
	t367 = cos(pkin(11));
	t375 = -t365 * t370 + t367 * t393;
	t350 = t375 * qJD(2);
	t379 = -t375 * t385 + t350;
	t356 = t365 * t393 + t367 * t370;
	t352 = t356 * qJD(2);
	t378 = t356 * t385 - t352;
	t377 = (qJD(2) - t385) * t372;
	t355 = t365 * t372 + t367 * t394;
	t344 = -t355 * t369 - t367 * t396;
	t376 = -t355 * t371 + t367 * t397;
	t357 = t367 * t372 - t384;
	t346 = -t357 * t369 + t365 * t396;
	t347 = t357 * t371 + t365 * t397;
	t359 = t368 * t369 + t370 * t396;
	t358 = t368 * t371 - t370 * t397;
	t351 = t355 * qJD(2);
	t374 = qJD(5) * t355 - t351 * t371 - t375 * t390;
	t353 = -qJD(2) * t384 + t367 * t391;
	t373 = qJD(5) * t357 - t353 * t371 + t356 * t390;
	t349 = t358 * qJD(3) + t371 * t382;
	t348 = -t359 * qJD(3) - t369 * t382;
	t343 = t346 * qJD(3) - t352 * t371;
	t342 = -t347 * qJD(3) + t352 * t369;
	t341 = t344 * qJD(3) + t350 * t371;
	t340 = t376 * qJD(3) - t350 * t369;
	t1 = [0, t378 * t362 + t373 * t363, t342 * t363 - t346 * t387, 0, -t343 * t362 + t353 * t363 + (-t347 * t363 - t356 * t362) * qJD(5), 0; 0, t379 * t362 + t374 * t363, t340 * t363 - t344 * t387, 0, -t341 * t362 + t351 * t363 + (t362 * t375 + t363 * t376) * qJD(5), 0; 0, (t362 * t377 - t398 * t363) * t366, t348 * t363 - t358 * t387, 0, t363 * t383 - t349 * t362 + (-t359 * t363 + t362 * t395) * qJD(5), 0; 0, -t373 * t362 + t378 * t363, -t342 * t362 - t346 * t386, 0, -t343 * t363 - t353 * t362 + (t347 * t362 - t356 * t363) * qJD(5), 0; 0, -t374 * t362 + t379 * t363, -t340 * t362 - t344 * t386, 0, -t341 * t363 - t351 * t362 + (-t362 * t376 + t363 * t375) * qJD(5), 0; 0, (t398 * t362 + t363 * t377) * t366, -t348 * t362 - t358 * t386, 0, -t362 * t383 - t349 * t363 + (t359 * t362 + t363 * t395) * qJD(5), 0; 0, -t353 * t369 - t356 * t389, t343, 0, 0, 0; 0, -t351 * t369 + t375 * t389, t341, 0, 0, 0; 0, (-t369 * t392 + t371 * t388) * t366, t349, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:28
	% EndTime: 2019-10-09 22:33:28
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (431->55), mult. (692->115), div. (0->0), fcn. (758->10), ass. (0->65)
	t421 = qJD(5) + qJD(6);
	t426 = sin(qJ(3));
	t427 = sin(qJ(2));
	t428 = cos(qJ(3));
	t429 = cos(qJ(2));
	t447 = qJD(3) * t429;
	t459 = (qJD(2) * t428 - t421) * t427 + t426 * t447;
	t420 = pkin(12) + qJ(5) + qJ(6);
	t418 = sin(t420);
	t458 = t418 * t421;
	t419 = cos(t420);
	t457 = t419 * t421;
	t456 = t421 * t428;
	t423 = sin(pkin(6));
	t455 = t423 * t426;
	t454 = t423 * t428;
	t425 = cos(pkin(6));
	t453 = t425 * t427;
	t452 = t425 * t429;
	t451 = qJD(2) * t427;
	t450 = qJD(2) * t429;
	t449 = qJD(3) * t426;
	t448 = qJD(3) * t428;
	t422 = sin(pkin(11));
	t446 = t422 * t453;
	t445 = t423 * t450;
	t424 = cos(pkin(11));
	t411 = t422 * t429 + t424 * t453;
	t400 = -t411 * t426 - t424 * t454;
	t433 = -t422 * t427 + t424 * t452;
	t406 = t433 * qJD(2);
	t397 = t400 * qJD(3) + t406 * t428;
	t443 = t421 * t433 - t397;
	t413 = t424 * t429 - t446;
	t402 = -t413 * t426 + t422 * t454;
	t412 = t422 * t452 + t424 * t427;
	t408 = t412 * qJD(2);
	t399 = t402 * qJD(3) - t408 * t428;
	t442 = -t412 * t421 - t399;
	t407 = t411 * qJD(2);
	t434 = -t411 * t428 + t424 * t455;
	t441 = -t421 * t434 - t407;
	t403 = t413 * t428 + t422 * t455;
	t409 = -qJD(2) * t446 + t424 * t450;
	t440 = t403 * t421 - t409;
	t414 = t425 * t428 - t427 * t455;
	t405 = t414 * qJD(3) + t428 * t445;
	t438 = t421 * t423 * t429 - t405;
	t437 = -t433 * t456 + t406;
	t436 = t412 * t456 - t408;
	t435 = (qJD(2) - t456) * t429;
	t415 = t425 * t426 + t427 * t454;
	t432 = -t415 * t421 + t423 * t451;
	t431 = -t407 * t428 + t411 * t421 - t433 * t449;
	t430 = -t409 * t428 + t412 * t449 + t413 * t421;
	t404 = -t415 * qJD(3) - t426 * t445;
	t398 = -t403 * qJD(3) + t408 * t426;
	t396 = t434 * qJD(3) - t406 * t426;
	t395 = -t432 * t418 + t438 * t419;
	t394 = t438 * t418 + t432 * t419;
	t393 = t440 * t418 + t442 * t419;
	t392 = t442 * t418 - t440 * t419;
	t391 = t441 * t418 + t443 * t419;
	t390 = t443 * t418 - t441 * t419;
	t1 = [0, t436 * t418 + t430 * t419, t398 * t419 - t402 * t458, 0, t392, t392; 0, t437 * t418 + t431 * t419, t396 * t419 - t400 * t458, 0, t390, t390; 0, (t418 * t435 - t459 * t419) * t423, t404 * t419 - t414 * t458, 0, t394, t394; 0, -t430 * t418 + t436 * t419, -t398 * t418 - t402 * t457, 0, t393, t393; 0, -t431 * t418 + t437 * t419, -t396 * t418 - t400 * t457, 0, t391, t391; 0, (t459 * t418 + t419 * t435) * t423, -t404 * t418 - t414 * t457, 0, t395, t395; 0, -t409 * t426 - t412 * t448, t399, 0, 0, 0; 0, -t407 * t426 + t433 * t448, t397, 0, 0, 0; 0, (-t426 * t451 + t428 * t447) * t423, t405, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end