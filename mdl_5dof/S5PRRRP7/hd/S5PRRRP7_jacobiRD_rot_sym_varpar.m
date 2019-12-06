% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRP7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:07
	% EndTime: 2019-12-05 16:57:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(5));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(5));
	t60 = cos(pkin(9));
	t58 = sin(pkin(9));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0; 0, -t62 * t64, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0; 0, -t63 * t64, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:07
	% EndTime: 2019-12-05 16:57:07
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(5));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(5));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(9));
	t212 = cos(pkin(9));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0; 0, t204, 0, 0, 0; 0, t202, 0, 0, 0; 0, t219, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:08
	% EndTime: 2019-12-05 16:57:08
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->53)
	t354 = sin(qJ(3));
	t355 = sin(qJ(2));
	t357 = cos(qJ(3));
	t358 = cos(qJ(2));
	t374 = qJD(3) * t358;
	t384 = (qJD(2) * t357 - qJD(4)) * t355 + t354 * t374;
	t350 = sin(pkin(5));
	t383 = t350 * t354;
	t382 = t350 * t357;
	t381 = t350 * t358;
	t352 = cos(pkin(5));
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
	t349 = sin(pkin(9));
	t370 = t349 * t380;
	t369 = t350 * t378;
	t368 = t350 * t377;
	t351 = cos(pkin(9));
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
	t336 = t345 * qJD(3) + t357 * t368;
	t335 = -t346 * qJD(3) - t354 * t368;
	t330 = t333 * qJD(3) - t339 * t357;
	t329 = -t334 * qJD(3) + t339 * t354;
	t328 = t331 * qJD(3) + t337 * t357;
	t327 = t362 * qJD(3) - t337 * t354;
	t1 = [0, t364 * t353 + t359 * t356, t329 * t356 - t333 * t373, -t330 * t353 + t340 * t356 + (-t334 * t356 - t343 * t353) * qJD(4), 0; 0, t365 * t353 + t360 * t356, t327 * t356 - t331 * t373, -t328 * t353 + t338 * t356 + (t353 * t361 + t356 * t362) * qJD(4), 0; 0, (t353 * t363 - t384 * t356) * t350, t335 * t356 - t345 * t373, t356 * t369 - t336 * t353 + (-t346 * t356 + t353 * t381) * qJD(4), 0; 0, -t359 * t353 + t364 * t356, -t329 * t353 - t333 * t372, -t330 * t356 - t340 * t353 + (t334 * t353 - t343 * t356) * qJD(4), 0; 0, -t360 * t353 + t365 * t356, -t327 * t353 - t331 * t372, -t328 * t356 - t338 * t353 + (-t353 * t362 + t356 * t361) * qJD(4), 0; 0, (t384 * t353 + t356 * t363) * t350, -t335 * t353 - t345 * t372, -t353 * t369 - t336 * t356 + (t346 * t353 + t356 * t381) * qJD(4), 0; 0, -t340 * t354 - t343 * t375, t330, 0, 0; 0, -t338 * t354 + t361 * t375, t328, 0, 0; 0, (-t354 * t378 + t357 * t374) * t350, t336, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:08
	% EndTime: 2019-12-05 16:57:08
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->53)
	t376 = sin(qJ(3));
	t377 = sin(qJ(2));
	t379 = cos(qJ(3));
	t380 = cos(qJ(2));
	t396 = qJD(3) * t380;
	t406 = (qJD(2) * t379 - qJD(4)) * t377 + t376 * t396;
	t372 = sin(pkin(5));
	t405 = t372 * t376;
	t404 = t372 * t379;
	t403 = t372 * t380;
	t374 = cos(pkin(5));
	t402 = t374 * t377;
	t401 = t374 * t380;
	t400 = qJD(2) * t377;
	t399 = qJD(2) * t380;
	t398 = qJD(3) * t376;
	t397 = qJD(3) * t379;
	t375 = sin(qJ(4));
	t395 = qJD(4) * t375;
	t378 = cos(qJ(4));
	t394 = qJD(4) * t378;
	t393 = qJD(4) * t379;
	t371 = sin(pkin(9));
	t392 = t371 * t402;
	t391 = t372 * t400;
	t390 = t372 * t399;
	t373 = cos(pkin(9));
	t383 = -t371 * t377 + t373 * t401;
	t359 = t383 * qJD(2);
	t387 = -t383 * t393 + t359;
	t365 = t371 * t401 + t373 * t377;
	t361 = t365 * qJD(2);
	t386 = t365 * t393 - t361;
	t385 = (qJD(2) - t393) * t380;
	t364 = t371 * t380 + t373 * t402;
	t353 = -t364 * t376 - t373 * t404;
	t384 = -t364 * t379 + t373 * t405;
	t366 = t373 * t380 - t392;
	t355 = -t366 * t376 + t371 * t404;
	t356 = t366 * t379 + t371 * t405;
	t368 = t374 * t376 + t377 * t404;
	t367 = t374 * t379 - t377 * t405;
	t360 = t364 * qJD(2);
	t382 = qJD(4) * t364 - t360 * t379 - t383 * t398;
	t362 = -qJD(2) * t392 + t373 * t399;
	t381 = qJD(4) * t366 - t362 * t379 + t365 * t398;
	t358 = t367 * qJD(3) + t379 * t390;
	t357 = -t368 * qJD(3) - t376 * t390;
	t352 = t355 * qJD(3) - t361 * t379;
	t351 = -t356 * qJD(3) + t361 * t376;
	t350 = t353 * qJD(3) + t359 * t379;
	t349 = t384 * qJD(3) - t359 * t376;
	t1 = [0, t386 * t375 + t381 * t378, t351 * t378 - t355 * t395, -t352 * t375 + t362 * t378 + (-t356 * t378 - t365 * t375) * qJD(4), 0; 0, t387 * t375 + t382 * t378, t349 * t378 - t353 * t395, -t350 * t375 + t360 * t378 + (t375 * t383 + t378 * t384) * qJD(4), 0; 0, (t375 * t385 - t406 * t378) * t372, t357 * t378 - t367 * t395, t378 * t391 - t358 * t375 + (-t368 * t378 + t375 * t403) * qJD(4), 0; 0, -t381 * t375 + t386 * t378, -t351 * t375 - t355 * t394, -t352 * t378 - t362 * t375 + (t356 * t375 - t365 * t378) * qJD(4), 0; 0, -t382 * t375 + t387 * t378, -t349 * t375 - t353 * t394, -t350 * t378 - t360 * t375 + (-t375 * t384 + t378 * t383) * qJD(4), 0; 0, (t406 * t375 + t378 * t385) * t372, -t357 * t375 - t367 * t394, -t375 * t391 - t358 * t378 + (t368 * t375 + t378 * t403) * qJD(4), 0; 0, -t362 * t376 - t365 * t397, t352, 0, 0; 0, -t360 * t376 + t383 * t397, t350, 0, 0; 0, (-t376 * t400 + t379 * t396) * t372, t358, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end