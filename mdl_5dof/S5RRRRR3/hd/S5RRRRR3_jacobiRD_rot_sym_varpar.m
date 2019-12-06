% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:57:18
	% EndTime: 2019-12-05 18:57:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:57:18
	% EndTime: 2019-12-05 18:57:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:57:18
	% EndTime: 2019-12-05 18:57:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t31 = sin(qJ(1));
	t38 = qJD(1) * t31;
	t33 = cos(qJ(1));
	t37 = qJD(1) * t33;
	t30 = sin(qJ(2));
	t36 = qJD(2) * t30;
	t32 = cos(qJ(2));
	t35 = qJD(2) * t32;
	t34 = qJD(2) * t33;
	t29 = t31 * t36 - t32 * t37;
	t28 = t30 * t37 + t31 * t35;
	t27 = t30 * t34 + t32 * t38;
	t26 = t30 * t38 - t32 * t34;
	t1 = [t29, t26, 0, 0, 0; -t27, -t28, 0, 0, 0; 0, -t36, 0, 0, 0; t28, t27, 0, 0, 0; t26, t29, 0, 0, 0; 0, -t35, 0, 0, 0; -t38, 0, 0, 0, 0; t37, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:57:18
	% EndTime: 2019-12-05 18:57:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t65 = qJ(2) + qJ(3);
	t62 = sin(t65);
	t64 = qJD(2) + qJD(3);
	t73 = t64 * t62;
	t63 = cos(t65);
	t72 = t64 * t63;
	t66 = sin(qJ(1));
	t71 = t64 * t66;
	t67 = cos(qJ(1));
	t70 = t64 * t67;
	t69 = qJD(1) * t66;
	t68 = qJD(1) * t67;
	t61 = t62 * t71 - t63 * t68;
	t60 = t62 * t68 + t63 * t71;
	t59 = t62 * t70 + t63 * t69;
	t58 = t62 * t69 - t63 * t70;
	t1 = [t61, t58, t58, 0, 0; -t59, -t60, -t60, 0, 0; 0, -t73, -t73, 0, 0; t60, t59, t59, 0, 0; t58, t61, t61, 0, 0; 0, -t72, -t72, 0, 0; -t69, 0, 0, 0, 0; t68, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:57:19
	% EndTime: 2019-12-05 18:57:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t326 = qJD(2) + qJD(3);
	t328 = sin(qJ(4));
	t350 = t326 * t328;
	t329 = sin(qJ(1));
	t349 = t326 * t329;
	t330 = cos(qJ(4));
	t348 = t326 * t330;
	t331 = cos(qJ(1));
	t347 = t326 * t331;
	t346 = t330 * t331;
	t345 = qJD(1) * t329;
	t344 = qJD(1) * t331;
	t343 = qJD(4) * t328;
	t342 = qJD(4) * t330;
	t341 = qJD(4) * t331;
	t327 = qJ(2) + qJ(3);
	t324 = sin(t327);
	t340 = t324 * t348;
	t339 = t324 * t349;
	t325 = cos(t327);
	t338 = t325 * t349;
	t337 = t324 * t347;
	t336 = t325 * t347;
	t335 = qJD(4) * t325 - qJD(1);
	t334 = qJD(1) * t325 - qJD(4);
	t333 = t335 * t328;
	t332 = t334 * t329 + t337;
	t323 = t326 * t325;
	t322 = t325 * t344 - t339;
	t321 = -t325 * t345 - t337;
	t320 = -t325 * t343 - t340;
	t319 = t324 * t350 - t325 * t342;
	t318 = -t330 * t338 + (t329 * t343 - t330 * t344) * t324;
	t317 = t328 * t338 + (t328 * t344 + t329 * t342) * t324;
	t316 = -t330 * t336 + (t328 * t341 + t330 * t345) * t324;
	t315 = t328 * t336 + (-t328 * t345 + t330 * t341) * t324;
	t314 = -t334 * t346 + (t333 + t340) * t329;
	t313 = t335 * t330 * t329 + (t334 * t331 - t339) * t328;
	t312 = t332 * t330 + t331 * t333;
	t311 = t332 * t328 - t335 * t346;
	t1 = [t314, t316, t316, t311, 0; -t312, t318, t318, -t313, 0; 0, t320, t320, -t324 * t342 - t325 * t350, 0; t313, t315, t315, t312, 0; t311, t317, t317, t314, 0; 0, t319, t319, t324 * t343 - t325 * t348, 0; -t324 * t344 - t338, t321, t321, 0, 0; -t324 * t345 + t336, t322, t322, 0, 0; 0, t323, t323, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:57:19
	% EndTime: 2019-12-05 18:57:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (344->32), mult. (286->53), div. (0->0), fcn. (286->6), ass. (0->47)
	t384 = qJ(2) + qJ(3);
	t378 = sin(t384);
	t381 = qJD(4) + qJD(5);
	t407 = t378 * t381;
	t382 = qJD(2) + qJD(3);
	t406 = t378 * t382;
	t380 = cos(t384);
	t405 = t380 * t381;
	t385 = sin(qJ(1));
	t404 = t381 * t385;
	t376 = t382 * t380;
	t403 = t382 * t385;
	t386 = cos(qJ(1));
	t402 = t382 * t386;
	t401 = qJD(1) * t385;
	t400 = qJD(1) * t386;
	t383 = qJ(4) + qJ(5);
	t379 = cos(t383);
	t399 = t379 * t407;
	t398 = t379 * t404;
	t397 = t378 * t403;
	t396 = t380 * t403;
	t395 = t378 * t402;
	t394 = t380 * t402;
	t393 = -qJD(1) + t405;
	t392 = qJD(1) * t380 - t381;
	t377 = sin(t383);
	t391 = t377 * t393;
	t390 = t386 * t381 * t377 + t379 * t401;
	t389 = t378 * t400 + t396;
	t388 = -t378 * t401 + t394;
	t387 = t392 * t385 + t395;
	t372 = t380 * t400 - t397;
	t371 = -t380 * t401 - t395;
	t370 = -t379 * t376 + t377 * t407;
	t369 = -t377 * t376 - t399;
	t368 = -t377 * t405 - t379 * t406;
	t367 = t377 * t406 - t379 * t405;
	t366 = -t379 * t396 + (t377 * t404 - t379 * t400) * t378;
	t365 = t389 * t377 + t378 * t398;
	t364 = t390 * t378 - t379 * t394;
	t363 = t388 * t377 + t386 * t399;
	t362 = t385 * t391 + (-t392 * t386 + t397) * t379;
	t361 = -t377 * t397 + (t377 * t400 + t398) * t380 - t390;
	t360 = t387 * t379 + t386 * t391;
	t359 = -t393 * t386 * t379 + t387 * t377;
	t1 = [t362, t364, t364, t359, t359; -t360, t366, t366, -t361, -t361; 0, t368, t368, t369, t369; t361, t363, t363, t360, t360; t359, t365, t365, t362, t362; 0, t367, t367, t370, t370; -t389, t371, t371, 0, 0; t388, t372, t372, 0, 0; 0, t376, t376, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end