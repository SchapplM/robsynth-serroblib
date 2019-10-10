% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
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
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 11:40:27
	% EndTime: 2019-10-10 11:40:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t235 = qJD(2) + qJD(3);
	t237 = sin(qJ(1));
	t242 = t235 * t237;
	t238 = cos(qJ(1));
	t241 = t235 * t238;
	t240 = qJD(1) * t237;
	t239 = qJD(1) * t238;
	t236 = qJ(2) + qJ(3);
	t234 = cos(t236);
	t233 = sin(t236);
	t232 = t235 * t234;
	t231 = t235 * t233;
	t230 = -t233 * t242 + t234 * t239;
	t229 = t233 * t239 + t234 * t242;
	t228 = t233 * t241 + t234 * t240;
	t227 = -t233 * t240 + t234 * t241;
	t1 = [-t240, 0, 0, 0, 0, 0; t239, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t230, t227, t227, 0, 0, 0; t228, t229, t229, 0, 0, 0; 0, t231, t231, 0, 0, 0; -t229, -t228, -t228, 0, 0, 0; t227, t230, t230, 0, 0, 0; 0, t232, t232, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:27
	% EndTime: 2019-10-10 11:40:28
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (166->31), mult. (226->56), div. (0->0), fcn. (226->6), ass. (0->42)
	t333 = cos(qJ(5));
	t330 = qJ(2) + qJ(3);
	t327 = sin(t330);
	t338 = qJD(5) * t327 + qJD(1);
	t355 = t333 * t338;
	t331 = sin(qJ(5));
	t354 = t338 * t331;
	t329 = qJD(2) + qJD(3);
	t353 = t329 * t327;
	t352 = t329 * t331;
	t332 = sin(qJ(1));
	t351 = t329 * t332;
	t350 = t329 * t333;
	t334 = cos(qJ(1));
	t349 = t329 * t334;
	t348 = qJD(1) * t332;
	t347 = qJD(1) * t334;
	t346 = qJD(5) * t331;
	t345 = qJD(5) * t333;
	t344 = qJD(5) * t334;
	t328 = cos(t330);
	t343 = t328 * t350;
	t342 = t327 * t351;
	t341 = t328 * t351;
	t340 = t327 * t349;
	t339 = t328 * t349;
	t337 = -qJD(1) * t327 - qJD(5);
	t336 = t337 * t334;
	t335 = t332 * t337 + t339;
	t326 = -t327 * t347 - t341;
	t325 = t327 * t348 - t339;
	t324 = -t327 * t346 + t343;
	t323 = t327 * t345 + t328 * t352;
	t322 = -t333 * t342 + (-t332 * t346 + t333 * t347) * t328;
	t321 = -t331 * t342 + (t331 * t347 + t332 * t345) * t328;
	t320 = -t333 * t340 + (-t331 * t344 - t333 * t348) * t328;
	t319 = -t331 * t340 + (-t331 * t348 + t333 * t344) * t328;
	t318 = t331 * t335 + t334 * t355;
	t317 = t333 * t335 - t334 * t354;
	t316 = -t332 * t355 + (t336 - t341) * t331;
	t315 = t333 * t336 + (-t343 + t354) * t332;
	t1 = [t316, t319, t319, 0, t317, 0; t318, t321, t321, 0, -t315, 0; 0, t323, t323, 0, t327 * t350 + t328 * t346, 0; t315, t320, t320, 0, -t318, 0; t317, t322, t322, 0, t316, 0; 0, t324, t324, 0, -t327 * t352 + t328 * t345, 0; -t328 * t347 + t342, t325, t325, 0, 0, 0; -t328 * t348 - t340, t326, t326, 0, 0, 0; 0, -t353, -t353, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:28
	% EndTime: 2019-10-10 11:40:28
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (166->31), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t408 = sin(qJ(1));
	t406 = qJ(2) + qJ(3);
	t403 = sin(t406);
	t413 = qJD(1) * t403 + qJD(5);
	t404 = cos(t406);
	t405 = qJD(2) + qJD(3);
	t410 = cos(qJ(1));
	t425 = t405 * t410;
	t415 = t404 * t425;
	t430 = t413 * t408 - t415;
	t429 = t405 * t403;
	t407 = sin(qJ(5));
	t428 = t405 * t407;
	t427 = t405 * t408;
	t409 = cos(qJ(5));
	t426 = t405 * t409;
	t424 = qJD(1) * t408;
	t423 = qJD(1) * t410;
	t422 = qJD(5) * t407;
	t421 = qJD(5) * t409;
	t420 = qJD(5) * t410;
	t419 = t404 * t426;
	t418 = t403 * t427;
	t417 = t404 * t427;
	t416 = t403 * t425;
	t414 = qJD(5) * t403 + qJD(1);
	t412 = t414 * t410;
	t411 = t413 * t410;
	t402 = -t403 * t423 - t417;
	t401 = t403 * t424 - t415;
	t400 = t403 * t422 - t419;
	t399 = t403 * t421 + t404 * t428;
	t398 = t409 * t418 + (t408 * t422 - t409 * t423) * t404;
	t397 = -t407 * t418 + (t407 * t423 + t408 * t421) * t404;
	t396 = t409 * t416 + (t407 * t420 + t409 * t424) * t404;
	t395 = -t407 * t416 + (-t407 * t424 + t409 * t420) * t404;
	t394 = -t430 * t407 + t409 * t412;
	t393 = t407 * t412 + t430 * t409;
	t392 = t414 * t409 * t408 + (t411 + t417) * t407;
	t391 = t409 * t411 + (-t414 * t407 + t419) * t408;
	t1 = [-t392, t395, t395, 0, -t393, 0; t394, t397, t397, 0, t391, 0; 0, t399, t399, 0, t403 * t426 + t404 * t422, 0; -t404 * t423 + t418, t401, t401, 0, 0, 0; -t404 * t424 - t416, t402, t402, 0, 0, 0; 0, -t429, -t429, 0, 0, 0; t391, t396, t396, 0, t394, 0; t393, t398, t398, 0, t392, 0; 0, t400, t400, 0, t403 * t428 - t404 * t421, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end