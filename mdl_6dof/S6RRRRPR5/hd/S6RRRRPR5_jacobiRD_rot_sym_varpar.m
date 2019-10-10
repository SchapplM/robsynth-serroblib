% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:39
	% EndTime: 2019-10-10 12:38:39
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
	% StartTime: 2019-10-10 12:38:40
	% EndTime: 2019-10-10 12:38:40
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t395 = cos(qJ(1));
	t391 = qJ(2) + qJ(3);
	t389 = cos(t391);
	t397 = qJD(1) * t389 - qJD(4);
	t414 = t397 * t395;
	t393 = sin(qJ(1));
	t388 = sin(t391);
	t390 = qJD(2) + qJD(3);
	t409 = t390 * t395;
	t400 = t388 * t409;
	t413 = t397 * t393 + t400;
	t392 = sin(qJ(4));
	t412 = t390 * t392;
	t411 = t390 * t393;
	t394 = cos(qJ(4));
	t410 = t390 * t394;
	t408 = qJD(1) * t393;
	t407 = qJD(1) * t395;
	t406 = qJD(4) * t392;
	t405 = qJD(4) * t394;
	t404 = qJD(4) * t395;
	t403 = t388 * t410;
	t402 = t388 * t411;
	t401 = t389 * t411;
	t399 = t389 * t409;
	t398 = -qJD(4) * t389 + qJD(1);
	t396 = t398 * t395;
	t387 = t390 * t389;
	t386 = t389 * t407 - t402;
	t385 = -t389 * t408 - t400;
	t384 = -t389 * t406 - t403;
	t383 = -t388 * t412 + t389 * t405;
	t382 = -t394 * t401 + (t393 * t406 - t394 * t407) * t388;
	t381 = -t392 * t401 + (-t392 * t407 - t393 * t405) * t388;
	t380 = -t394 * t399 + (t392 * t404 + t394 * t408) * t388;
	t379 = -t392 * t399 + (t392 * t408 - t394 * t404) * t388;
	t378 = t394 * t414 + (t398 * t392 - t403) * t393;
	t377 = t398 * t394 * t393 + (t402 - t414) * t392;
	t376 = t392 * t396 - t413 * t394;
	t375 = t413 * t392 + t394 * t396;
	t1 = [-t378, t380, t380, t375, 0, 0; t376, t382, t382, t377, 0, 0; 0, t384, t384, -t388 * t405 - t389 * t412, 0, 0; -t388 * t407 - t401, t385, t385, 0, 0, 0; -t388 * t408 + t399, t386, t386, 0, 0, 0; 0, t387, t387, 0, 0, 0; t377, t379, t379, t376, 0, 0; -t375, t381, t381, t378, 0, 0; 0, t383, t383, -t388 * t406 + t389 * t410, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:40
	% EndTime: 2019-10-10 12:38:41
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (546->60), mult. (882->92), div. (0->0), fcn. (906->8), ass. (0->57)
	t484 = sin(qJ(6));
	t485 = sin(qJ(4));
	t487 = cos(qJ(6));
	t488 = cos(qJ(4));
	t499 = t484 * t485 + t487 * t488;
	t524 = qJD(4) - qJD(6);
	t490 = t524 * t499;
	t500 = t484 * t488 - t485 * t487;
	t525 = t524 * t500;
	t483 = qJ(2) + qJ(3);
	t481 = cos(t483);
	t480 = sin(t483);
	t482 = qJD(2) + qJD(3);
	t486 = sin(qJ(1));
	t519 = t482 * t486;
	t509 = t480 * t519;
	t489 = cos(qJ(1));
	t512 = qJD(1) * t489;
	t465 = -t481 * t512 + t509;
	t514 = t489 * t488;
	t517 = t486 * t485;
	t466 = t481 * t517 + t514;
	t510 = qJD(4) * t489;
	t506 = t488 * t510;
	t511 = qJD(4) * t486;
	t507 = t485 * t511;
	t518 = t482 * t489;
	t508 = t480 * t518;
	t460 = t466 * qJD(1) - t481 * t506 + t485 * t508 - t507;
	t515 = t489 * t485;
	t461 = (-qJD(4) * t481 + qJD(1)) * t515 + (-t508 + (-qJD(1) * t481 + qJD(4)) * t486) * t488;
	t516 = t486 * t488;
	t468 = t481 * t515 - t516;
	t469 = t481 * t514 + t517;
	t523 = t460 * t487 + t461 * t484 + (t468 * t484 + t469 * t487) * qJD(6);
	t513 = qJD(1) * t486;
	t462 = (t511 * t481 - t513) * t488 + (-t465 - t510) * t485;
	t463 = t469 * qJD(1) - t481 * t507 - t488 * t509 - t506;
	t467 = t481 * t516 - t515;
	t522 = t462 * t484 + t463 * t487 + (t466 * t487 - t467 * t484) * qJD(6);
	t494 = (t466 * t484 + t467 * t487) * qJD(6) - t462 * t487 + t463 * t484;
	t520 = t482 * t481;
	t498 = t500 * t482;
	t497 = t499 * t482;
	t496 = t481 * t497;
	t495 = t481 * t498;
	t450 = -t460 * t484 + t461 * t487 + (t468 * t487 - t469 * t484) * qJD(6);
	t464 = t481 * t513 + t508;
	t458 = t480 * t525 + t499 * t520;
	t457 = -t490 * t480 + t495;
	t456 = -t480 * t497 + t481 * t525;
	t455 = t480 * t498 + t481 * t490;
	t454 = -t486 * t496 + (-t486 * t525 - t499 * t512) * t480;
	t453 = t486 * t495 + (-t490 * t486 + t500 * t512) * t480;
	t452 = -t489 * t496 + (-t489 * t525 + t499 * t513) * t480;
	t451 = t489 * t495 + (-t490 * t489 - t500 * t513) * t480;
	t1 = [-t522, t452, t452, t523, 0, -t523; t450, t454, t454, t494, 0, -t494; 0, t456, t456, t457, 0, -t457; t494, t451, t451, t450, 0, -t450; -t523, t453, t453, t522, 0, -t522; 0, t455, t455, t458, 0, -t458; t480 * t512 + t481 * t519, t464, t464, 0, 0, 0; t480 * t513 - t481 * t518, t465, t465, 0, 0, 0; 0, -t520, -t520, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end