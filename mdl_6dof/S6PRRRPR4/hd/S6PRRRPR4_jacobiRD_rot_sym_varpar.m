% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 22:52:16
	% EndTime: 2019-10-09 22:52:16
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
	% StartTime: 2019-10-09 22:52:17
	% EndTime: 2019-10-09 22:52:17
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
	t349 = sin(pkin(11));
	t370 = t349 * t380;
	t369 = t350 * t378;
	t368 = t350 * t377;
	t351 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:52:17
	% EndTime: 2019-10-09 22:52:18
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (219->60), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->54)
	t378 = sin(qJ(3));
	t379 = sin(qJ(2));
	t380 = cos(qJ(3));
	t381 = cos(qJ(2));
	t397 = qJD(3) * t381;
	t407 = (qJD(2) * t380 - qJD(4)) * t379 + t378 * t397;
	t375 = sin(pkin(6));
	t406 = t375 * t378;
	t405 = t375 * t380;
	t404 = t375 * t381;
	t377 = cos(pkin(6));
	t403 = t377 * t379;
	t402 = t377 * t381;
	t401 = qJD(2) * t379;
	t400 = qJD(2) * t381;
	t399 = qJD(3) * t378;
	t398 = qJD(3) * t380;
	t373 = qJ(4) + pkin(12);
	t371 = sin(t373);
	t396 = qJD(4) * t371;
	t372 = cos(t373);
	t395 = qJD(4) * t372;
	t394 = qJD(4) * t380;
	t374 = sin(pkin(11));
	t393 = t374 * t403;
	t392 = t375 * t401;
	t391 = t375 * t400;
	t376 = cos(pkin(11));
	t384 = -t374 * t379 + t376 * t402;
	t359 = t384 * qJD(2);
	t388 = -t384 * t394 + t359;
	t365 = t374 * t402 + t376 * t379;
	t361 = t365 * qJD(2);
	t387 = t365 * t394 - t361;
	t386 = (qJD(2) - t394) * t381;
	t364 = t374 * t381 + t376 * t403;
	t353 = -t364 * t378 - t376 * t405;
	t385 = -t364 * t380 + t376 * t406;
	t366 = t376 * t381 - t393;
	t355 = -t366 * t378 + t374 * t405;
	t356 = t366 * t380 + t374 * t406;
	t368 = t377 * t378 + t379 * t405;
	t367 = t377 * t380 - t379 * t406;
	t360 = t364 * qJD(2);
	t383 = qJD(4) * t364 - t360 * t380 - t384 * t399;
	t362 = -qJD(2) * t393 + t376 * t400;
	t382 = qJD(4) * t366 - t362 * t380 + t365 * t399;
	t358 = t367 * qJD(3) + t380 * t391;
	t357 = -t368 * qJD(3) - t378 * t391;
	t352 = t355 * qJD(3) - t361 * t380;
	t351 = -t356 * qJD(3) + t361 * t378;
	t350 = t353 * qJD(3) + t359 * t380;
	t349 = t385 * qJD(3) - t359 * t378;
	t1 = [0, t387 * t371 + t382 * t372, t351 * t372 - t355 * t396, -t352 * t371 + t362 * t372 + (-t356 * t372 - t365 * t371) * qJD(4), 0, 0; 0, t388 * t371 + t383 * t372, t349 * t372 - t353 * t396, -t350 * t371 + t360 * t372 + (t371 * t384 + t372 * t385) * qJD(4), 0, 0; 0, (t371 * t386 - t407 * t372) * t375, t357 * t372 - t367 * t396, t372 * t392 - t358 * t371 + (-t368 * t372 + t371 * t404) * qJD(4), 0, 0; 0, -t382 * t371 + t387 * t372, -t351 * t371 - t355 * t395, -t352 * t372 - t362 * t371 + (t356 * t371 - t365 * t372) * qJD(4), 0, 0; 0, -t383 * t371 + t388 * t372, -t349 * t371 - t353 * t395, -t350 * t372 - t360 * t371 + (-t371 * t385 + t372 * t384) * qJD(4), 0, 0; 0, (t407 * t371 + t372 * t386) * t375, -t357 * t371 - t367 * t395, -t371 * t392 - t358 * t372 + (t368 * t371 + t372 * t404) * qJD(4), 0, 0; 0, -t362 * t378 - t365 * t398, t352, 0, 0, 0; 0, -t360 * t378 + t384 * t398, t350, 0, 0, 0; 0, (-t378 * t401 + t380 * t397) * t375, t358, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:17
	% EndTime: 2019-10-09 22:52:18
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (431->55), mult. (692->115), div. (0->0), fcn. (758->10), ass. (0->65)
	t426 = qJD(4) + qJD(6);
	t431 = sin(qJ(3));
	t432 = sin(qJ(2));
	t433 = cos(qJ(3));
	t434 = cos(qJ(2));
	t452 = qJD(3) * t434;
	t464 = (qJD(2) * t433 - t426) * t432 + t431 * t452;
	t425 = qJ(4) + pkin(12) + qJ(6);
	t423 = sin(t425);
	t463 = t423 * t426;
	t424 = cos(t425);
	t462 = t424 * t426;
	t461 = t426 * t433;
	t428 = sin(pkin(6));
	t460 = t428 * t431;
	t459 = t428 * t433;
	t430 = cos(pkin(6));
	t458 = t430 * t432;
	t457 = t430 * t434;
	t456 = qJD(2) * t432;
	t455 = qJD(2) * t434;
	t454 = qJD(3) * t431;
	t453 = qJD(3) * t433;
	t427 = sin(pkin(11));
	t451 = t427 * t458;
	t450 = t428 * t455;
	t429 = cos(pkin(11));
	t416 = t427 * t434 + t429 * t458;
	t405 = -t416 * t431 - t429 * t459;
	t438 = -t427 * t432 + t429 * t457;
	t411 = t438 * qJD(2);
	t402 = t405 * qJD(3) + t411 * t433;
	t448 = t426 * t438 - t402;
	t418 = t429 * t434 - t451;
	t407 = -t418 * t431 + t427 * t459;
	t417 = t427 * t457 + t429 * t432;
	t413 = t417 * qJD(2);
	t404 = t407 * qJD(3) - t413 * t433;
	t447 = -t417 * t426 - t404;
	t412 = t416 * qJD(2);
	t439 = -t416 * t433 + t429 * t460;
	t446 = -t426 * t439 - t412;
	t408 = t418 * t433 + t427 * t460;
	t414 = -qJD(2) * t451 + t429 * t455;
	t445 = t408 * t426 - t414;
	t419 = t430 * t433 - t432 * t460;
	t410 = t419 * qJD(3) + t433 * t450;
	t443 = t426 * t428 * t434 - t410;
	t442 = -t438 * t461 + t411;
	t441 = t417 * t461 - t413;
	t440 = (qJD(2) - t461) * t434;
	t420 = t430 * t431 + t432 * t459;
	t437 = -t420 * t426 + t428 * t456;
	t436 = -t412 * t433 + t416 * t426 - t438 * t454;
	t435 = -t414 * t433 + t417 * t454 + t418 * t426;
	t409 = -t420 * qJD(3) - t431 * t450;
	t403 = -t408 * qJD(3) + t413 * t431;
	t401 = t439 * qJD(3) - t411 * t431;
	t400 = -t437 * t423 + t443 * t424;
	t399 = t443 * t423 + t437 * t424;
	t398 = t445 * t423 + t447 * t424;
	t397 = t447 * t423 - t445 * t424;
	t396 = t446 * t423 + t448 * t424;
	t395 = t448 * t423 - t446 * t424;
	t1 = [0, t441 * t423 + t435 * t424, t403 * t424 - t407 * t463, t397, 0, t397; 0, t442 * t423 + t436 * t424, t401 * t424 - t405 * t463, t395, 0, t395; 0, (t423 * t440 - t464 * t424) * t428, t409 * t424 - t419 * t463, t399, 0, t399; 0, -t435 * t423 + t441 * t424, -t403 * t423 - t407 * t462, t398, 0, t398; 0, -t436 * t423 + t442 * t424, -t401 * t423 - t405 * t462, t396, 0, t396; 0, (t464 * t423 + t424 * t440) * t428, -t409 * t423 - t419 * t462, t400, 0, t400; 0, -t414 * t431 - t417 * t453, t404, 0, 0, 0; 0, -t412 * t431 + t438 * t453, t402, 0, 0, 0; 0, (-t431 * t456 + t433 * t452) * t428, t410, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end