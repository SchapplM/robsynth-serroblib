% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(2) + qJ(3);
	t68 = sin(t71);
	t70 = qJD(2) + qJD(3);
	t79 = t70 * t68;
	t69 = cos(t71);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t67 = t68 * t77 - t69 * t74;
	t66 = t68 * t74 + t69 * t77;
	t65 = t68 * t76 + t69 * t75;
	t64 = t68 * t75 - t69 * t76;
	t1 = [t67, t64, t64, 0, 0, 0; -t65, -t66, -t66, 0, 0, 0; 0, -t79, -t79, 0, 0, 0; t66, t65, t65, 0, 0, 0; t64, t67, t67, 0, 0, 0; 0, -t78, -t78, 0, 0, 0; -t75, 0, 0, 0, 0, 0; t74, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:30
	% EndTime: 2019-10-10 12:20:30
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t329 = qJD(2) + qJD(3);
	t331 = sin(qJ(4));
	t353 = t329 * t331;
	t332 = sin(qJ(1));
	t352 = t329 * t332;
	t333 = cos(qJ(4));
	t351 = t329 * t333;
	t334 = cos(qJ(1));
	t350 = t329 * t334;
	t349 = t333 * t334;
	t348 = qJD(1) * t332;
	t347 = qJD(1) * t334;
	t346 = qJD(4) * t331;
	t345 = qJD(4) * t333;
	t344 = qJD(4) * t334;
	t330 = qJ(2) + qJ(3);
	t327 = sin(t330);
	t343 = t327 * t351;
	t342 = t327 * t352;
	t328 = cos(t330);
	t341 = t328 * t352;
	t340 = t327 * t350;
	t339 = t328 * t350;
	t338 = qJD(4) * t328 - qJD(1);
	t337 = qJD(1) * t328 - qJD(4);
	t336 = t338 * t331;
	t335 = t337 * t332 + t340;
	t326 = t329 * t328;
	t325 = t328 * t347 - t342;
	t324 = -t328 * t348 - t340;
	t323 = -t328 * t346 - t343;
	t322 = t327 * t353 - t328 * t345;
	t321 = -t333 * t341 + (t332 * t346 - t333 * t347) * t327;
	t320 = t331 * t341 + (t331 * t347 + t332 * t345) * t327;
	t319 = -t333 * t339 + (t331 * t344 + t333 * t348) * t327;
	t318 = t331 * t339 + (-t331 * t348 + t333 * t344) * t327;
	t317 = -t337 * t349 + (t336 + t343) * t332;
	t316 = t338 * t333 * t332 + (t337 * t334 - t342) * t331;
	t315 = t335 * t333 + t334 * t336;
	t314 = t335 * t331 - t338 * t349;
	t1 = [t317, t319, t319, t314, 0, 0; -t315, t321, t321, -t316, 0, 0; 0, t323, t323, -t327 * t345 - t328 * t353, 0, 0; t316, t318, t318, t315, 0, 0; t314, t320, t320, t317, 0, 0; 0, t322, t322, t327 * t346 - t328 * t351, 0, 0; -t327 * t347 - t341, t324, t324, 0, 0, 0; -t327 * t348 + t339, t325, t325, 0, 0, 0; 0, t326, t326, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:30
	% EndTime: 2019-10-10 12:20:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (240->27), mult. (226->45), div. (0->0), fcn. (226->6), ass. (0->42)
	t356 = cos(qJ(1));
	t354 = qJ(2) + qJ(3);
	t351 = cos(t354);
	t367 = qJD(4) * t351;
	t362 = -qJD(1) + t367;
	t375 = t356 * t362;
	t361 = qJD(1) * t351 - qJD(4);
	t350 = sin(t354);
	t352 = qJD(2) + qJD(3);
	t355 = sin(qJ(1));
	t372 = t352 * t355;
	t366 = t350 * t372;
	t374 = t361 * t356 - t366;
	t373 = t350 * t352;
	t347 = t352 * t351;
	t371 = t352 * t356;
	t370 = qJD(1) * t355;
	t369 = qJD(1) * t356;
	t368 = qJD(4) * t350;
	t365 = t350 * t371;
	t353 = qJ(4) + pkin(10);
	t348 = sin(t353);
	t364 = t348 * t368;
	t349 = cos(t353);
	t363 = t349 * t368;
	t360 = t362 * t355;
	t359 = t350 * t369 + t351 * t372;
	t358 = t350 * t370 - t351 * t371;
	t357 = t361 * t355 + t365;
	t346 = t351 * t369 - t366;
	t345 = -t351 * t370 - t365;
	t344 = -t348 * t367 - t349 * t373;
	t343 = t348 * t373 - t349 * t367;
	t342 = -t359 * t349 + t355 * t364;
	t341 = t359 * t348 + t355 * t363;
	t340 = t358 * t349 + t356 * t364;
	t339 = -t358 * t348 + t356 * t363;
	t338 = t348 * t360 - t374 * t349;
	t337 = t374 * t348 + t349 * t360;
	t336 = t348 * t375 + t357 * t349;
	t335 = t357 * t348 - t349 * t375;
	t1 = [t338, t340, t340, t335, 0, 0; -t336, t342, t342, -t337, 0, 0; 0, t344, t344, -t348 * t347 - t363, 0, 0; t337, t339, t339, t336, 0, 0; t335, t341, t341, t338, 0, 0; 0, t343, t343, -t349 * t347 + t364, 0, 0; -t359, t345, t345, 0, 0, 0; -t358, t346, t346, 0, 0, 0; 0, t347, t347, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:31
	% EndTime: 2019-10-10 12:20:31
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (240->26), mult. (226->45), div. (0->0), fcn. (226->6), ass. (0->42)
	t417 = sin(qJ(1));
	t416 = qJ(2) + qJ(3);
	t413 = cos(t416);
	t423 = qJD(1) * t413 - qJD(4);
	t412 = sin(t416);
	t414 = qJD(2) + qJD(3);
	t418 = cos(qJ(1));
	t433 = t414 * t418;
	t427 = t412 * t433;
	t437 = t423 * t417 + t427;
	t434 = t414 * t417;
	t428 = t412 * t434;
	t436 = t423 * t418 - t428;
	t435 = t412 * t414;
	t409 = t414 * t413;
	t432 = qJD(1) * t417;
	t431 = qJD(1) * t418;
	t430 = qJD(4) * t412;
	t429 = qJD(4) * t413;
	t415 = qJ(4) + pkin(10);
	t410 = sin(t415);
	t426 = t410 * t430;
	t411 = cos(t415);
	t425 = t411 * t430;
	t424 = qJD(1) - t429;
	t422 = t424 * t417;
	t421 = t424 * t418;
	t420 = -t412 * t431 - t413 * t434;
	t419 = t412 * t432 - t413 * t433;
	t408 = t413 * t431 - t428;
	t407 = -t413 * t432 - t427;
	t406 = -t410 * t429 - t411 * t435;
	t405 = -t410 * t435 + t411 * t429;
	t404 = t420 * t411 + t417 * t426;
	t403 = t420 * t410 - t417 * t425;
	t402 = t419 * t411 + t418 * t426;
	t401 = t419 * t410 - t418 * t425;
	t400 = t410 * t422 + t436 * t411;
	t399 = -t436 * t410 + t411 * t422;
	t398 = t410 * t421 - t437 * t411;
	t397 = t437 * t410 + t411 * t421;
	t1 = [-t400, t402, t402, t397, 0, 0; t398, t404, t404, t399, 0, 0; 0, t406, t406, -t410 * t409 - t425, 0, 0; t420, t407, t407, 0, 0, 0; -t419, t408, t408, 0, 0, 0; 0, t409, t409, 0, 0, 0; t399, t401, t401, t398, 0, 0; -t397, t403, t403, t400, 0, 0; 0, t405, t405, t411 * t409 - t426, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end