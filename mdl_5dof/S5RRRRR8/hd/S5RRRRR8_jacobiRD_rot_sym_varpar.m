% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.06s
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
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.06s
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
	t1 = [t67, t64, t64, 0, 0; -t65, -t66, -t66, 0, 0; 0, -t79, -t79, 0, 0; t66, t65, t65, 0, 0; t64, t67, t67, 0, 0; 0, -t78, -t78, 0, 0; -t75, 0, 0, 0, 0; t74, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:56
	% EndTime: 2019-12-29 20:55:56
	% DurationCPUTime: 0.20s
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
	t1 = [t317, t319, t319, t314, 0; -t315, t321, t321, -t316, 0; 0, t323, t323, -t327 * t345 - t328 * t353, 0; t316, t318, t318, t315, 0; t314, t320, t320, t317, 0; 0, t322, t322, t327 * t346 - t328 * t351, 0; -t327 * t347 - t341, t324, t324, 0, 0; -t327 * t348 + t339, t325, t325, 0, 0; 0, t326, t326, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:56
	% EndTime: 2019-12-29 20:55:56
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (344->32), mult. (286->53), div. (0->0), fcn. (286->6), ass. (0->47)
	t390 = qJ(2) + qJ(3);
	t384 = sin(t390);
	t387 = qJD(4) + qJD(5);
	t413 = t384 * t387;
	t388 = qJD(2) + qJD(3);
	t412 = t384 * t388;
	t386 = cos(t390);
	t411 = t386 * t387;
	t391 = sin(qJ(1));
	t410 = t387 * t391;
	t382 = t388 * t386;
	t409 = t388 * t391;
	t392 = cos(qJ(1));
	t408 = t388 * t392;
	t407 = qJD(1) * t391;
	t406 = qJD(1) * t392;
	t389 = qJ(4) + qJ(5);
	t385 = cos(t389);
	t405 = t385 * t413;
	t404 = t385 * t410;
	t403 = t384 * t409;
	t402 = t386 * t409;
	t401 = t384 * t408;
	t400 = t386 * t408;
	t399 = -qJD(1) + t411;
	t398 = qJD(1) * t386 - t387;
	t383 = sin(t389);
	t397 = t383 * t399;
	t396 = t392 * t387 * t383 + t385 * t407;
	t395 = t384 * t406 + t402;
	t394 = -t384 * t407 + t400;
	t393 = t398 * t391 + t401;
	t378 = t386 * t406 - t403;
	t377 = -t386 * t407 - t401;
	t376 = -t385 * t382 + t383 * t413;
	t375 = -t383 * t382 - t405;
	t374 = -t383 * t411 - t385 * t412;
	t373 = t383 * t412 - t385 * t411;
	t372 = -t385 * t402 + (t383 * t410 - t385 * t406) * t384;
	t371 = t395 * t383 + t384 * t404;
	t370 = t396 * t384 - t385 * t400;
	t369 = t394 * t383 + t392 * t405;
	t368 = t391 * t397 + (-t398 * t392 + t403) * t385;
	t367 = -t383 * t403 + (t383 * t406 + t404) * t386 - t396;
	t366 = t393 * t385 + t392 * t397;
	t365 = -t399 * t392 * t385 + t393 * t383;
	t1 = [t368, t370, t370, t365, t365; -t366, t372, t372, -t367, -t367; 0, t374, t374, t375, t375; t367, t369, t369, t366, t366; t365, t371, t371, t368, t368; 0, t373, t373, t376, t376; -t395, t377, t377, 0, 0; t394, t378, t378, 0, 0; 0, t382, t382, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end