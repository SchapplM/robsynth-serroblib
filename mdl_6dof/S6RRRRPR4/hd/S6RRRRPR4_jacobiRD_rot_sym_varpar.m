% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:36
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
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
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 12:36:52
	% EndTime: 2019-10-10 12:36:52
	% DurationCPUTime: 0.16s
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
	% StartTime: 2019-10-10 12:36:53
	% EndTime: 2019-10-10 12:36:53
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
	t353 = qJ(4) + pkin(11);
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
	% StartTime: 2019-10-10 12:36:53
	% EndTime: 2019-10-10 12:36:53
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (444->32), mult. (286->49), div. (0->0), fcn. (286->6), ass. (0->45)
	t395 = qJ(2) + qJ(3);
	t391 = sin(t395);
	t393 = qJD(4) + qJD(6);
	t416 = t391 * t393;
	t394 = qJD(2) + qJD(3);
	t415 = t391 * t394;
	t392 = cos(t395);
	t414 = t392 * t393;
	t387 = t394 * t392;
	t396 = sin(qJ(1));
	t413 = t394 * t396;
	t397 = cos(qJ(1));
	t412 = t394 * t397;
	t411 = qJD(1) * t396;
	t410 = qJD(1) * t397;
	t390 = qJ(4) + pkin(11) + qJ(6);
	t388 = sin(t390);
	t409 = t388 * t416;
	t389 = cos(t390);
	t408 = t389 * t416;
	t407 = t396 * t393 * t389;
	t406 = t391 * t413;
	t405 = t397 * t393 * t388;
	t404 = t391 * t412;
	t403 = -qJD(1) + t414;
	t402 = qJD(1) * t392 - t393;
	t401 = t388 * t403;
	t400 = t391 * t410 + t392 * t413;
	t399 = t391 * t411 - t392 * t412;
	t398 = t402 * t396 + t404;
	t383 = t392 * t410 - t406;
	t382 = -t392 * t411 - t404;
	t381 = -t389 * t387 + t409;
	t380 = -t388 * t387 - t408;
	t379 = -t388 * t414 - t389 * t415;
	t378 = t388 * t415 - t389 * t414;
	t377 = -t400 * t389 + t396 * t409;
	t376 = t400 * t388 + t391 * t407;
	t375 = t399 * t389 + t391 * t405;
	t374 = -t399 * t388 + t397 * t408;
	t373 = t396 * t401 + (-t402 * t397 + t406) * t389;
	t372 = -t388 * t406 - t405 - t389 * t411 + (t388 * t410 + t407) * t392;
	t371 = t398 * t389 + t397 * t401;
	t370 = -t403 * t397 * t389 + t398 * t388;
	t1 = [t373, t375, t375, t370, 0, t370; -t371, t377, t377, -t372, 0, -t372; 0, t379, t379, t380, 0, t380; t372, t374, t374, t371, 0, t371; t370, t376, t376, t373, 0, t373; 0, t378, t378, t381, 0, t381; -t400, t382, t382, 0, 0, 0; -t399, t383, t383, 0, 0, 0; 0, t387, t387, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end