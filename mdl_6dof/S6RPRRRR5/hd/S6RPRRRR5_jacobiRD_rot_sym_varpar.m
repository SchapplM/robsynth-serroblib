% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
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
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(11));
	t17 = sin(pkin(11));
	t1 = [-t18 * t21, 0, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(11) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0, 0; -t37, 0, -t38, 0, 0, 0; 0, 0, -qJD(3) * t40, 0, 0, 0; t38, 0, t37, 0, 0, 0; t36, 0, t39, 0, 0, 0; 0, 0, -qJD(3) * t41, 0, 0, 0; -t48, 0, 0, 0, 0, 0; t47, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (89->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = pkin(11) + qJ(3) + qJ(4);
	t70 = sin(t72);
	t73 = qJD(3) + qJD(4);
	t81 = t73 * t70;
	t71 = cos(t72);
	t80 = t73 * t71;
	t74 = sin(qJ(1));
	t79 = t73 * t74;
	t75 = cos(qJ(1));
	t78 = t73 * t75;
	t77 = qJD(1) * t74;
	t76 = qJD(1) * t75;
	t69 = t70 * t79 - t71 * t76;
	t68 = t70 * t76 + t71 * t79;
	t67 = t70 * t78 + t71 * t77;
	t66 = t70 * t77 - t71 * t78;
	t1 = [t69, 0, t66, t66, 0, 0; -t67, 0, -t68, -t68, 0, 0; 0, 0, -t81, -t81, 0, 0; t68, 0, t67, t67, 0, 0; t66, 0, t69, t69, 0, 0; 0, 0, -t80, -t80, 0, 0; -t77, 0, 0, 0, 0, 0; t76, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:09
	% EndTime: 2019-10-10 09:06:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (238->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t332 = qJD(3) + qJD(4);
	t333 = sin(qJ(5));
	t355 = t332 * t333;
	t334 = sin(qJ(1));
	t354 = t332 * t334;
	t335 = cos(qJ(5));
	t353 = t332 * t335;
	t336 = cos(qJ(1));
	t352 = t332 * t336;
	t351 = t335 * t336;
	t350 = qJD(1) * t334;
	t349 = qJD(1) * t336;
	t348 = qJD(5) * t333;
	t347 = qJD(5) * t335;
	t346 = qJD(5) * t336;
	t331 = pkin(11) + qJ(3) + qJ(4);
	t329 = sin(t331);
	t345 = t329 * t353;
	t344 = t329 * t354;
	t330 = cos(t331);
	t343 = t330 * t354;
	t342 = t329 * t352;
	t341 = t330 * t352;
	t340 = qJD(5) * t330 - qJD(1);
	t339 = qJD(1) * t330 - qJD(5);
	t338 = t340 * t333;
	t337 = t339 * t334 + t342;
	t328 = t332 * t330;
	t327 = t330 * t349 - t344;
	t326 = -t330 * t350 - t342;
	t325 = -t330 * t348 - t345;
	t324 = t329 * t355 - t330 * t347;
	t323 = -t335 * t343 + (t334 * t348 - t335 * t349) * t329;
	t322 = t333 * t343 + (t333 * t349 + t334 * t347) * t329;
	t321 = -t335 * t341 + (t333 * t346 + t335 * t350) * t329;
	t320 = t333 * t341 + (-t333 * t350 + t335 * t346) * t329;
	t319 = -t339 * t351 + (t338 + t345) * t334;
	t318 = t340 * t335 * t334 + (t339 * t336 - t344) * t333;
	t317 = t337 * t335 + t336 * t338;
	t316 = t337 * t333 - t340 * t351;
	t1 = [t319, 0, t321, t321, t316, 0; -t317, 0, t323, t323, -t318, 0; 0, 0, t325, t325, -t329 * t347 - t330 * t355, 0; t318, 0, t320, t320, t317, 0; t316, 0, t322, t322, t319, 0; 0, 0, t324, t324, t329 * t348 - t330 * t353, 0; -t329 * t349 - t343, 0, t326, t326, 0, 0; -t329 * t350 + t341, 0, t327, t327, 0, 0; 0, 0, t328, t328, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:09
	% EndTime: 2019-10-10 09:06:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (434->32), mult. (286->56), div. (0->0), fcn. (286->6), ass. (0->45)
	t388 = pkin(11) + qJ(3) + qJ(4);
	t386 = sin(t388);
	t392 = qJD(3) + qJD(4);
	t414 = t386 * t392;
	t393 = qJ(5) + qJ(6);
	t389 = sin(t393);
	t391 = qJD(5) + qJD(6);
	t413 = t389 * t391;
	t390 = cos(t393);
	t412 = t390 * t391;
	t394 = sin(qJ(1));
	t411 = t391 * t394;
	t395 = cos(qJ(1));
	t410 = t391 * t395;
	t387 = cos(t388);
	t384 = t392 * t387;
	t409 = t392 * t394;
	t408 = t392 * t395;
	t407 = qJD(1) * t394;
	t406 = qJD(1) * t395;
	t405 = t386 * t409;
	t404 = t387 * t409;
	t403 = t386 * t408;
	t402 = t387 * t408;
	t401 = t387 * t391 - qJD(1);
	t400 = qJD(1) * t387 - t391;
	t399 = t389 * t401;
	t398 = t389 * t410 + t390 * t407;
	t397 = t389 * t406 + t390 * t411;
	t396 = t400 * t394 + t403;
	t381 = t387 * t406 - t405;
	t380 = -t387 * t407 - t403;
	t379 = -t390 * t384 + t386 * t413;
	t378 = -t389 * t384 - t386 * t412;
	t377 = -t387 * t413 - t390 * t414;
	t376 = -t387 * t412 + t389 * t414;
	t375 = -t390 * t404 + (t389 * t411 - t390 * t406) * t386;
	t374 = t397 * t386 + t389 * t404;
	t373 = t398 * t386 - t390 * t402;
	t372 = t389 * t402 + (-t389 * t407 + t390 * t410) * t386;
	t371 = t394 * t399 + (-t400 * t395 + t405) * t390;
	t370 = t397 * t387 - t389 * t405 - t398;
	t369 = t396 * t390 + t395 * t399;
	t368 = -t401 * t395 * t390 + t396 * t389;
	t1 = [t371, 0, t373, t373, t368, t368; -t369, 0, t375, t375, -t370, -t370; 0, 0, t377, t377, t378, t378; t370, 0, t372, t372, t369, t369; t368, 0, t374, t374, t371, t371; 0, 0, t376, t376, t379, t379; -t386 * t406 - t404, 0, t380, t380, 0, 0; -t386 * t407 + t402, 0, t381, t381, 0, 0; 0, 0, t384, t384, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end