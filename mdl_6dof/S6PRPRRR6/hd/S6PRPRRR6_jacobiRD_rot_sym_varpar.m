% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
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
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t139 = cos(pkin(6));
	t140 = sin(qJ(2));
	t144 = t139 * t140;
	t141 = cos(qJ(2));
	t143 = t139 * t141;
	t142 = qJD(2) * sin(pkin(6));
	t138 = cos(pkin(11));
	t136 = sin(pkin(11));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, (-t136 * t144 + t138 * t141) * qJD(2), 0, 0, 0, 0; 0, (t136 * t141 + t138 * t144) * qJD(2), 0, 0, 0, 0; 0, t140 * t142, 0, 0, 0, 0; 0, (-t136 * t143 - t138 * t140) * qJD(2), 0, 0, 0, 0; 0, (-t136 * t140 + t138 * t143) * qJD(2), 0, 0, 0, 0; 0, t141 * t142, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (37->25), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t209 = sin(pkin(6));
	t212 = sin(qJ(4));
	t225 = t209 * t212;
	t214 = cos(qJ(4));
	t224 = t209 * t214;
	t215 = cos(qJ(2));
	t223 = t209 * t215;
	t211 = cos(pkin(6));
	t213 = sin(qJ(2));
	t222 = t211 * t213;
	t221 = t211 * t215;
	t220 = qJD(2) * t215;
	t219 = qJD(4) * t212;
	t218 = qJD(4) * t214;
	t217 = t209 * qJD(2) * t213;
	t208 = sin(pkin(11));
	t210 = cos(pkin(11));
	t216 = -t208 * t213 + t210 * t221;
	t205 = t208 * t215 + t210 * t222;
	t206 = t208 * t221 + t210 * t213;
	t207 = -t208 * t222 + t210 * t215;
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t201 = t205 * qJD(2);
	t200 = t216 * qJD(2);
	t1 = [0, -t202 * t212 + t207 * t218, 0, t203 * t214 + (-t206 * t212 - t208 * t224) * qJD(4), 0, 0; 0, t200 * t212 + t205 * t218, 0, t201 * t214 + (t210 * t224 + t212 * t216) * qJD(4), 0, 0; 0, (t212 * t220 + t213 * t218) * t209, 0, t214 * t217 + (-t211 * t214 + t212 * t223) * qJD(4), 0, 0; 0, -t202 * t214 - t207 * t219, 0, -t203 * t212 + (-t206 * t214 + t208 * t225) * qJD(4), 0, 0; 0, t200 * t214 - t205 * t219, 0, -t201 * t212 + (-t210 * t225 + t214 * t216) * qJD(4), 0, 0; 0, (-t213 * t219 + t214 * t220) * t209, 0, -t212 * t217 + (t211 * t212 + t214 * t223) * qJD(4), 0, 0; 0, -t203, 0, 0, 0, 0; 0, -t201, 0, 0, 0, 0; 0, -t217, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:09
	% EndTime: 2019-10-09 22:03:09
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (153->60), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->53)
	t359 = sin(qJ(4));
	t363 = cos(qJ(2));
	t389 = (qJD(2) * t359 + qJD(5)) * t363;
	t355 = sin(pkin(6));
	t388 = t355 * t359;
	t360 = sin(qJ(2));
	t387 = t355 * t360;
	t362 = cos(qJ(4));
	t386 = t355 * t362;
	t385 = t355 * t363;
	t357 = cos(pkin(6));
	t384 = t357 * t360;
	t383 = t357 * t363;
	t361 = cos(qJ(5));
	t382 = t360 * t361;
	t381 = qJD(2) * t363;
	t380 = qJD(4) * t359;
	t379 = qJD(4) * t362;
	t358 = sin(qJ(5));
	t378 = qJD(5) * t358;
	t377 = qJD(5) * t359;
	t376 = qJD(5) * t361;
	t354 = sin(pkin(11));
	t375 = t354 * t384;
	t374 = qJD(2) * t387;
	t373 = t355 * t381;
	t372 = -qJD(2) - t377;
	t356 = cos(pkin(11));
	t347 = t354 * t363 + t356 * t384;
	t343 = t347 * qJD(2);
	t370 = -t347 * t377 - t343;
	t345 = -qJD(2) * t375 + t356 * t381;
	t349 = t356 * t363 - t375;
	t369 = -t349 * t377 - t345;
	t367 = -t354 * t360 + t356 * t383;
	t368 = t356 * t386 + t359 * t367;
	t338 = t356 * t388 - t362 * t367;
	t348 = t354 * t383 + t356 * t360;
	t337 = t348 * t359 + t354 * t386;
	t336 = t348 * t362 - t354 * t388;
	t350 = -t357 * t359 - t362 * t385;
	t366 = -t357 * t362 + t359 * t385;
	t342 = t367 * qJD(2);
	t365 = qJD(5) * t367 + t342 * t359 + t347 * t379;
	t344 = t348 * qJD(2);
	t364 = -qJD(5) * t348 - t344 * t359 + t349 * t379;
	t341 = t366 * qJD(4) + t362 * t374;
	t340 = t350 * qJD(4) + t359 * t374;
	t335 = t368 * qJD(4) + t343 * t362;
	t334 = t338 * qJD(4) + t343 * t359;
	t333 = -t337 * qJD(4) + t345 * t362;
	t332 = t336 * qJD(4) + t345 * t359;
	t1 = [0, t369 * t358 + t364 * t361, 0, t333 * t361 - t336 * t378, -t332 * t358 - t344 * t361 + (-t337 * t361 - t349 * t358) * qJD(5), 0; 0, t370 * t358 + t365 * t361, 0, t335 * t361 - t338 * t378, -t334 * t358 + t342 * t361 + (-t347 * t358 + t361 * t368) * qJD(5), 0; 0, (t361 * t389 + (t372 * t358 + t361 * t379) * t360) * t355, 0, t341 * t361 - t350 * t378, t361 * t373 - t340 * t358 + (-t358 * t387 + t361 * t366) * qJD(5), 0; 0, -t364 * t358 + t369 * t361, 0, -t333 * t358 - t336 * t376, -t332 * t361 + t344 * t358 + (t337 * t358 - t349 * t361) * qJD(5), 0; 0, -t365 * t358 + t370 * t361, 0, -t335 * t358 - t338 * t376, -t334 * t361 - t342 * t358 + (-t347 * t361 - t358 * t368) * qJD(5), 0; 0, (t372 * t382 + (-t360 * t379 - t389) * t358) * t355, 0, -t341 * t358 - t350 * t376, -t358 * t373 - t340 * t361 + (-t355 * t382 - t358 * t366) * qJD(5), 0; 0, t344 * t362 + t349 * t380, 0, t332, 0, 0; 0, -t342 * t362 + t347 * t380, 0, t334, 0, 0; 0, (t360 * t380 - t362 * t381) * t355, 0, t340, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:09
	% EndTime: 2019-10-09 22:03:09
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (341->55), mult. (692->114), div. (0->0), fcn. (758->10), ass. (0->65)
	t416 = qJD(5) + qJD(6);
	t422 = sin(qJ(4));
	t423 = sin(qJ(2));
	t425 = cos(qJ(2));
	t424 = cos(qJ(4));
	t444 = qJD(4) * t424;
	t456 = (qJD(2) * t422 + t416) * t425 + t423 * t444;
	t417 = qJ(5) + qJ(6);
	t414 = sin(t417);
	t455 = t414 * t416;
	t415 = cos(t417);
	t454 = t415 * t416;
	t453 = t416 * t422;
	t419 = sin(pkin(6));
	t452 = t419 * t422;
	t451 = t419 * t423;
	t450 = t419 * t424;
	t449 = t419 * t425;
	t421 = cos(pkin(6));
	t448 = t421 * t423;
	t447 = t421 * t425;
	t446 = qJD(2) * t425;
	t445 = qJD(4) * t422;
	t418 = sin(pkin(11));
	t443 = t418 * t448;
	t442 = qJD(2) * t451;
	t420 = cos(pkin(11));
	t408 = t418 * t447 + t420 * t423;
	t396 = t408 * t424 - t418 * t452;
	t405 = -qJD(2) * t443 + t420 * t446;
	t392 = t396 * qJD(4) + t405 * t422;
	t409 = t420 * t425 - t443;
	t440 = -t409 * t416 - t392;
	t430 = -t418 * t423 + t420 * t447;
	t398 = t420 * t452 - t424 * t430;
	t407 = t418 * t425 + t420 * t448;
	t403 = t407 * qJD(2);
	t394 = t398 * qJD(4) + t403 * t422;
	t439 = -t407 * t416 - t394;
	t397 = t408 * t422 + t418 * t450;
	t404 = t408 * qJD(2);
	t438 = t397 * t416 + t404;
	t402 = t430 * qJD(2);
	t431 = t420 * t450 + t422 * t430;
	t437 = -t416 * t431 - t402;
	t410 = -t421 * t422 - t424 * t449;
	t400 = t410 * qJD(4) + t422 * t442;
	t435 = -t416 * t451 - t400;
	t434 = -t407 * t453 - t403;
	t433 = -t409 * t453 - t405;
	t432 = (-qJD(2) - t453) * t423;
	t429 = -t421 * t424 + t422 * t449;
	t428 = t416 * t429 + t419 * t446;
	t427 = t402 * t422 + t407 * t444 + t416 * t430;
	t426 = -t404 * t422 - t408 * t416 + t409 * t444;
	t401 = t429 * qJD(4) + t424 * t442;
	t395 = t431 * qJD(4) + t403 * t424;
	t393 = -t397 * qJD(4) + t405 * t424;
	t391 = -t428 * t414 + t435 * t415;
	t390 = t435 * t414 + t428 * t415;
	t389 = t437 * t414 + t439 * t415;
	t388 = t439 * t414 - t437 * t415;
	t387 = t438 * t414 + t440 * t415;
	t386 = t440 * t414 - t438 * t415;
	t1 = [0, t433 * t414 + t426 * t415, 0, t393 * t415 - t396 * t455, t386, t386; 0, t434 * t414 + t427 * t415, 0, t395 * t415 - t398 * t455, t388, t388; 0, (t414 * t432 + t456 * t415) * t419, 0, t401 * t415 - t410 * t455, t390, t390; 0, -t426 * t414 + t433 * t415, 0, -t393 * t414 - t396 * t454, t387, t387; 0, -t427 * t414 + t434 * t415, 0, -t395 * t414 - t398 * t454, t389, t389; 0, (-t456 * t414 + t415 * t432) * t419, 0, -t401 * t414 - t410 * t454, t391, t391; 0, t404 * t424 + t409 * t445, 0, t392, 0, 0; 0, -t402 * t424 + t407 * t445, 0, t394, 0, 0; 0, (t423 * t445 - t424 * t446) * t419, 0, t400, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end