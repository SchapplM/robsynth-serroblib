% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:02
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-10 12:24:02
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t396 = sin(qJ(1));
	t394 = qJ(2) + qJ(3);
	t392 = cos(t394);
	t401 = qJD(1) * t392 - qJD(4);
	t391 = sin(t394);
	t393 = qJD(2) + qJD(3);
	t398 = cos(qJ(1));
	t413 = t393 * t398;
	t404 = t391 * t413;
	t417 = t401 * t396 + t404;
	t395 = sin(qJ(4));
	t416 = t393 * t395;
	t415 = t393 * t396;
	t397 = cos(qJ(4));
	t414 = t393 * t397;
	t412 = qJD(1) * t396;
	t411 = qJD(1) * t398;
	t410 = qJD(4) * t395;
	t409 = qJD(4) * t397;
	t408 = qJD(4) * t398;
	t407 = t391 * t414;
	t406 = t391 * t415;
	t405 = t392 * t415;
	t403 = t392 * t413;
	t402 = qJD(4) * t392 - qJD(1);
	t400 = t402 * t398;
	t399 = t401 * t398;
	t390 = t393 * t392;
	t389 = t392 * t411 - t406;
	t388 = -t392 * t412 - t404;
	t387 = t392 * t410 + t407;
	t386 = -t391 * t416 + t392 * t409;
	t385 = t397 * t405 + (-t396 * t410 + t397 * t411) * t391;
	t384 = -t395 * t405 + (-t395 * t411 - t396 * t409) * t391;
	t383 = t397 * t403 + (-t395 * t408 - t397 * t412) * t391;
	t382 = -t395 * t403 + (t395 * t412 - t397 * t408) * t391;
	t381 = t397 * t399 + (-t402 * t395 - t407) * t396;
	t380 = t402 * t397 * t396 + (t399 - t406) * t395;
	t379 = t395 * t400 + t417 * t397;
	t378 = -t417 * t395 + t397 * t400;
	t1 = [-t391 * t411 - t405, t388, t388, 0, 0, 0; -t391 * t412 + t403, t389, t389, 0, 0, 0; 0, t390, t390, 0, 0, 0; t381, t383, t383, t378, 0, 0; t379, t385, t385, t380, 0, 0; 0, t387, t387, t391 * t409 + t392 * t416, 0, 0; -t380, t382, t382, -t379, 0, 0; t378, t384, t384, t381, 0, 0; 0, t386, t386, -t391 * t410 + t392 * t414, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:02
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t401 = cos(qJ(1));
	t397 = qJ(2) + qJ(3);
	t395 = cos(t397);
	t403 = qJD(1) * t395 - qJD(4);
	t420 = t403 * t401;
	t399 = sin(qJ(1));
	t394 = sin(t397);
	t396 = qJD(2) + qJD(3);
	t415 = t396 * t401;
	t406 = t394 * t415;
	t419 = t403 * t399 + t406;
	t398 = sin(qJ(4));
	t418 = t396 * t398;
	t417 = t396 * t399;
	t400 = cos(qJ(4));
	t416 = t396 * t400;
	t414 = qJD(1) * t399;
	t413 = qJD(1) * t401;
	t412 = qJD(4) * t398;
	t411 = qJD(4) * t400;
	t410 = qJD(4) * t401;
	t409 = t394 * t416;
	t408 = t394 * t417;
	t407 = t395 * t417;
	t405 = t395 * t415;
	t404 = -qJD(4) * t395 + qJD(1);
	t402 = t404 * t401;
	t393 = t396 * t395;
	t392 = t395 * t413 - t408;
	t391 = -t395 * t414 - t406;
	t390 = -t395 * t412 - t409;
	t389 = -t394 * t418 + t395 * t411;
	t388 = -t400 * t407 + (t399 * t412 - t400 * t413) * t394;
	t387 = -t398 * t407 + (-t398 * t413 - t399 * t411) * t394;
	t386 = -t400 * t405 + (t398 * t410 + t400 * t414) * t394;
	t385 = -t398 * t405 + (t398 * t414 - t400 * t410) * t394;
	t384 = t400 * t420 + (t404 * t398 - t409) * t399;
	t383 = t404 * t400 * t399 + (t408 - t420) * t398;
	t382 = t398 * t402 - t419 * t400;
	t381 = t419 * t398 + t400 * t402;
	t1 = [-t394 * t413 - t407, t391, t391, 0, 0, 0; -t394 * t414 + t405, t392, t392, 0, 0, 0; 0, t393, t393, 0, 0, 0; t383, t385, t385, t382, 0, 0; -t381, t387, t387, t384, 0, 0; 0, t389, t389, -t394 * t412 + t395 * t416, 0, 0; -t384, t386, t386, t381, 0, 0; t382, t388, t388, t383, 0, 0; 0, t390, t390, -t394 * t411 - t395 * t418, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end