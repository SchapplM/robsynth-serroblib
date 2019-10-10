% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:02
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:02
	% EndTime: 2019-10-10 10:54:02
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
	% StartTime: 2019-10-10 10:54:02
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(1));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t54 = qJD(1) * t51;
	t53 = qJD(2) * t50;
	t52 = qJD(2) * t51;
	t49 = qJ(2) + pkin(11);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = t47 * t53 - t48 * t54;
	t45 = t47 * t54 + t48 * t53;
	t44 = t47 * t52 + t48 * t55;
	t43 = t47 * t55 - t48 * t52;
	t1 = [t46, t43, 0, 0, 0, 0; -t44, -t45, 0, 0, 0, 0; 0, -qJD(2) * t47, 0, 0, 0, 0; t45, t44, 0, 0, 0, 0; t43, t46, 0, 0, 0, 0; 0, -qJD(2) * t48, 0, 0, 0, 0; -t55, 0, 0, 0, 0, 0; t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:02
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (89->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t77 = qJ(2) + pkin(11) + qJ(4);
	t75 = sin(t77);
	t78 = qJD(2) + qJD(4);
	t86 = t78 * t75;
	t76 = cos(t77);
	t85 = t78 * t76;
	t79 = sin(qJ(1));
	t84 = t78 * t79;
	t80 = cos(qJ(1));
	t83 = t78 * t80;
	t82 = qJD(1) * t79;
	t81 = qJD(1) * t80;
	t74 = t75 * t84 - t76 * t81;
	t73 = t75 * t81 + t76 * t84;
	t72 = t75 * t83 + t76 * t82;
	t71 = t75 * t82 - t76 * t83;
	t1 = [t74, t71, 0, t71, 0, 0; -t72, -t73, 0, -t73, 0, 0; 0, -t86, 0, -t86, 0, 0; t73, t72, 0, t72, 0, 0; t71, t74, 0, t74, 0, 0; 0, -t85, 0, -t85, 0, 0; -t82, 0, 0, 0, 0, 0; t81, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:03
	% EndTime: 2019-10-10 10:54:03
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (238->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t337 = qJD(2) + qJD(4);
	t338 = sin(qJ(5));
	t360 = t337 * t338;
	t339 = sin(qJ(1));
	t359 = t337 * t339;
	t340 = cos(qJ(5));
	t358 = t337 * t340;
	t341 = cos(qJ(1));
	t357 = t337 * t341;
	t356 = t340 * t341;
	t355 = qJD(1) * t339;
	t354 = qJD(1) * t341;
	t353 = qJD(5) * t338;
	t352 = qJD(5) * t340;
	t351 = qJD(5) * t341;
	t336 = qJ(2) + pkin(11) + qJ(4);
	t334 = sin(t336);
	t350 = t334 * t358;
	t349 = t334 * t359;
	t335 = cos(t336);
	t348 = t335 * t359;
	t347 = t334 * t357;
	t346 = t335 * t357;
	t345 = qJD(5) * t335 - qJD(1);
	t344 = qJD(1) * t335 - qJD(5);
	t343 = t345 * t338;
	t342 = t344 * t339 + t347;
	t333 = t337 * t335;
	t332 = t335 * t354 - t349;
	t331 = -t335 * t355 - t347;
	t330 = -t335 * t353 - t350;
	t329 = t334 * t360 - t335 * t352;
	t328 = -t340 * t348 + (t339 * t353 - t340 * t354) * t334;
	t327 = t338 * t348 + (t338 * t354 + t339 * t352) * t334;
	t326 = -t340 * t346 + (t338 * t351 + t340 * t355) * t334;
	t325 = t338 * t346 + (-t338 * t355 + t340 * t351) * t334;
	t324 = -t344 * t356 + (t343 + t350) * t339;
	t323 = t345 * t340 * t339 + (t344 * t341 - t349) * t338;
	t322 = t342 * t340 + t341 * t343;
	t321 = t342 * t338 - t345 * t356;
	t1 = [t324, t326, 0, t326, t321, 0; -t322, t328, 0, t328, -t323, 0; 0, t330, 0, t330, -t334 * t352 - t335 * t360, 0; t323, t325, 0, t325, t322, 0; t321, t327, 0, t327, t324, 0; 0, t329, 0, t329, t334 * t353 - t335 * t358, 0; -t334 * t354 - t348, t331, 0, t331, 0, 0; -t334 * t355 + t346, t332, 0, t332, 0, 0; 0, t333, 0, t333, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:04
	% EndTime: 2019-10-10 10:54:04
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (434->32), mult. (286->56), div. (0->0), fcn. (286->6), ass. (0->45)
	t393 = qJ(2) + pkin(11) + qJ(4);
	t391 = sin(t393);
	t397 = qJD(2) + qJD(4);
	t419 = t391 * t397;
	t398 = qJ(5) + qJ(6);
	t394 = sin(t398);
	t396 = qJD(5) + qJD(6);
	t418 = t394 * t396;
	t395 = cos(t398);
	t417 = t395 * t396;
	t399 = sin(qJ(1));
	t416 = t396 * t399;
	t400 = cos(qJ(1));
	t415 = t396 * t400;
	t392 = cos(t393);
	t389 = t397 * t392;
	t414 = t397 * t399;
	t413 = t397 * t400;
	t412 = qJD(1) * t399;
	t411 = qJD(1) * t400;
	t410 = t391 * t414;
	t409 = t392 * t414;
	t408 = t391 * t413;
	t407 = t392 * t413;
	t406 = t392 * t396 - qJD(1);
	t405 = qJD(1) * t392 - t396;
	t404 = t394 * t406;
	t403 = t394 * t415 + t395 * t412;
	t402 = t394 * t411 + t395 * t416;
	t401 = t405 * t399 + t408;
	t386 = t392 * t411 - t410;
	t385 = -t392 * t412 - t408;
	t384 = -t395 * t389 + t391 * t418;
	t383 = -t394 * t389 - t391 * t417;
	t382 = -t392 * t418 - t395 * t419;
	t381 = -t392 * t417 + t394 * t419;
	t380 = -t395 * t409 + (t394 * t416 - t395 * t411) * t391;
	t379 = t402 * t391 + t394 * t409;
	t378 = t403 * t391 - t395 * t407;
	t377 = t394 * t407 + (-t394 * t412 + t395 * t415) * t391;
	t376 = t399 * t404 + (-t405 * t400 + t410) * t395;
	t375 = t402 * t392 - t394 * t410 - t403;
	t374 = t401 * t395 + t400 * t404;
	t373 = -t406 * t400 * t395 + t401 * t394;
	t1 = [t376, t378, 0, t378, t373, t373; -t374, t380, 0, t380, -t375, -t375; 0, t382, 0, t382, t383, t383; t375, t377, 0, t377, t374, t374; t373, t379, 0, t379, t376, t376; 0, t381, 0, t381, t384, t384; -t391 * t411 - t409, t385, 0, t385, 0, 0; -t391 * t412 + t407, t386, 0, t386, 0, 0; 0, t389, 0, t389, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end