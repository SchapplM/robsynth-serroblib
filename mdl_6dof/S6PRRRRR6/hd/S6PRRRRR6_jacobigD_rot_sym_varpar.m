% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t106 = sin(qJ(2));
	t109 = cos(pkin(6)) * t106;
	t108 = qJD(2) * sin(pkin(7));
	t107 = cos(qJ(2));
	t104 = cos(pkin(14));
	t102 = sin(pkin(14));
	t1 = [0, 0, -(t102 * t109 - t104 * t107) * t108, 0, 0, 0; 0, 0, -(-t102 * t107 - t104 * t109) * t108, 0, 0, 0; 0, 0, sin(pkin(6)) * t106 * t108, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (29->20), mult. (109->51), div. (0->0), fcn. (113->12), ass. (0->27)
	t222 = sin(pkin(7));
	t244 = t222 * cos(pkin(8));
	t223 = sin(pkin(6));
	t243 = t223 * t222;
	t226 = cos(pkin(7));
	t230 = cos(qJ(3));
	t242 = t226 * t230;
	t227 = cos(pkin(6));
	t229 = sin(qJ(2));
	t241 = t227 * t229;
	t231 = cos(qJ(2));
	t240 = t227 * t231;
	t228 = sin(qJ(3));
	t239 = t228 * t231;
	t238 = t229 * t230;
	t237 = qJD(2) * t228;
	t221 = sin(pkin(8));
	t236 = qJD(3) * t221;
	t220 = sin(pkin(14));
	t224 = cos(pkin(14));
	t235 = -t220 * t229 + t224 * t240;
	t234 = t220 * t231 + t224 * t241;
	t233 = t220 * t240 + t224 * t229;
	t232 = t220 * t241 - t224 * t231;
	t219 = t232 * qJD(2);
	t218 = t234 * qJD(2);
	t1 = [0, 0, -t219 * t222, -(t219 * t242 + t233 * t237) * t221 - t219 * t244 - (t232 * t230 + (-t220 * t243 + t233 * t226) * t228) * t236, 0, 0; 0, 0, t218 * t222, -(-t218 * t242 - t235 * t237) * t221 + t218 * t244 - (-t234 * t230 + (t224 * t243 - t235 * t226) * t228) * t236, 0, 0; 0, 0, qJD(2) * t229 * t243, t227 * t222 * t228 * t236 + (-(-t226 * t239 - t238) * t236 + (-(-t226 * t238 - t239) * t221 + t229 * t244) * qJD(2)) * t223, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:37
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (111->50), mult. (386->112), div. (0->0), fcn. (423->14), ass. (0->51)
	t325 = sin(qJ(3));
	t328 = cos(qJ(3));
	t316 = sin(pkin(14));
	t320 = cos(pkin(14));
	t329 = cos(qJ(2));
	t323 = cos(pkin(6));
	t326 = sin(qJ(2));
	t345 = t323 * t326;
	t334 = t316 * t345 - t320 * t329;
	t344 = t323 * t329;
	t314 = -t316 * t344 - t320 * t326;
	t322 = cos(pkin(7));
	t318 = sin(pkin(7));
	t319 = sin(pkin(6));
	t352 = t318 * t319;
	t335 = t314 * t322 + t316 * t352;
	t357 = t335 * t325 - t328 * t334;
	t313 = t316 * t329 + t320 * t345;
	t312 = -t316 * t326 + t320 * t344;
	t336 = -t312 * t322 + t320 * t352;
	t356 = -t313 * t328 + t336 * t325;
	t317 = sin(pkin(8));
	t353 = t318 * t317;
	t321 = cos(pkin(8));
	t351 = t318 * t321;
	t350 = t318 * t323;
	t349 = t318 * t326;
	t348 = t319 * t322;
	t347 = t322 * t325;
	t346 = t322 * t328;
	t343 = t325 * t326;
	t342 = t325 * t329;
	t341 = t326 * t328;
	t340 = t328 * t329;
	t324 = sin(qJ(4));
	t339 = qJD(3) * t324;
	t338 = qJD(3) * t350;
	t337 = t319 * qJD(2) * t349;
	t333 = t322 * t340 - t343;
	t332 = t322 * t342 + t341;
	t331 = -t313 * t325 - t336 * t328;
	t330 = t325 * t334 + t335 * t328;
	t327 = cos(qJ(4));
	t311 = t334 * qJD(2);
	t310 = t314 * qJD(2);
	t309 = t313 * qJD(2);
	t308 = t312 * qJD(2);
	t307 = -t325 * t338 + (-t332 * qJD(3) + (-t322 * t341 - t342) * qJD(2)) * t319;
	t306 = -t357 * qJD(3) - t310 * t325 + t311 * t346;
	t305 = t356 * qJD(3) - t308 * t325 - t309 * t346;
	t1 = [0, 0, -t311 * t318, -t306 * t317 - t311 * t351, (t310 * t328 + t311 * t347) * t324 + (-t306 * t321 + t311 * t353) * t327 + t330 * t339 + (t357 * t327 + (t330 * t321 + (-t314 * t318 + t316 * t348) * t317) * t324) * qJD(4), 0; 0, 0, t309 * t318, -t305 * t317 + t309 * t351, (t308 * t328 - t309 * t347) * t324 + (-t305 * t321 - t309 * t353) * t327 + t331 * t339 + (-t356 * t327 + (t331 * t321 + (-t312 * t318 - t320 * t348) * t317) * t324) * qJD(4), 0; 0, 0, t337, -t307 * t317 + t321 * t337, t328 * t324 * t338 - t307 * t321 * t327 + (t333 * t339 + ((-t322 * t343 + t340) * t324 - t317 * t327 * t349) * qJD(2)) * t319 + ((t332 * t319 + t325 * t350) * t327 + ((t333 * t319 + t328 * t350) * t321 + (t323 * t322 - t329 * t352) * t317) * t324) * qJD(4), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:38
	% EndTime: 2019-10-09 23:23:39
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (276->74), mult. (919->150), div. (0->0), fcn. (1042->16), ass. (0->72)
	t428 = sin(qJ(3));
	t432 = cos(qJ(3));
	t418 = sin(pkin(14));
	t422 = cos(pkin(14));
	t433 = cos(qJ(2));
	t425 = cos(pkin(6));
	t429 = sin(qJ(2));
	t455 = t425 * t429;
	t439 = t418 * t455 - t422 * t433;
	t454 = t425 * t433;
	t416 = -t418 * t454 - t422 * t429;
	t424 = cos(pkin(7));
	t420 = sin(pkin(7));
	t421 = sin(pkin(6));
	t462 = t420 * t421;
	t440 = t416 * t424 + t418 * t462;
	t404 = t440 * t428 - t432 * t439;
	t415 = t418 * t433 + t422 * t455;
	t414 = -t418 * t429 + t422 * t454;
	t441 = -t414 * t424 + t422 * t462;
	t466 = -t415 * t432 + t441 * t428;
	t419 = sin(pkin(8));
	t463 = t419 * t420;
	t423 = cos(pkin(8));
	t461 = t420 * t423;
	t460 = t420 * t425;
	t459 = t421 * t424;
	t427 = sin(qJ(4));
	t458 = t423 * t427;
	t457 = t424 * t428;
	t456 = t424 * t432;
	t453 = t428 * t429;
	t452 = t428 * t433;
	t451 = t429 * t432;
	t450 = t432 * t433;
	t426 = sin(qJ(5));
	t449 = qJD(4) * t426;
	t448 = t427 * t463;
	t447 = qJD(3) * t460;
	t446 = qJD(2) * t429 * t462;
	t445 = t419 * t446;
	t401 = -t415 * t428 - t441 * t432;
	t407 = -t414 * t420 - t422 * t459;
	t444 = t401 * t423 + t407 * t419;
	t403 = t428 * t439 + t440 * t432;
	t408 = -t416 * t420 + t418 * t459;
	t443 = t403 * t423 + t408 * t419;
	t438 = t424 * t450 - t453;
	t405 = t438 * t421 + t432 * t460;
	t413 = t425 * t424 - t433 * t462;
	t442 = t405 * t423 + t413 * t419;
	t437 = t424 * t452 + t451;
	t431 = cos(qJ(4));
	t436 = t444 * t427 - t431 * t466;
	t435 = t404 * t431 + t443 * t427;
	t406 = t437 * t421 + t428 * t460;
	t434 = t406 * t431 + t442 * t427;
	t430 = cos(qJ(5));
	t412 = t439 * qJD(2);
	t411 = t416 * qJD(2);
	t410 = t415 * qJD(2);
	t409 = t414 * qJD(2);
	t400 = t432 * t447 + (t438 * qJD(3) + (-t424 * t453 + t450) * qJD(2)) * t421;
	t399 = -t428 * t447 + (-t437 * qJD(3) + (-t424 * t451 - t452) * qJD(2)) * t421;
	t398 = -t399 * t419 + t423 * t446;
	t397 = t403 * qJD(3) + t411 * t432 + t412 * t457;
	t396 = -t404 * qJD(3) - t411 * t428 + t412 * t456;
	t395 = t401 * qJD(3) + t409 * t432 - t410 * t457;
	t394 = t466 * qJD(3) - t409 * t428 - t410 * t456;
	t393 = -t396 * t419 - t412 * t461;
	t392 = -t394 * t419 + t410 * t461;
	t1 = [0, 0, -t412 * t420, t393, t397 * t427 + (-t396 * t423 + t412 * t463) * t431 + t435 * qJD(4), (t396 * t458 + t397 * t431 - t412 * t448) * t426 - t393 * t430 + (t435 * t430 + (-t403 * t419 + t408 * t423) * t426) * qJD(5) + (-t404 * t427 + t443 * t431) * t449; 0, 0, t410 * t420, t392, t395 * t427 + (-t394 * t423 - t410 * t463) * t431 + t436 * qJD(4), (t394 * t458 + t395 * t431 + t410 * t448) * t426 - t392 * t430 + (t436 * t430 + (-t401 * t419 + t407 * t423) * t426) * qJD(5) + (t427 * t466 + t444 * t431) * t449; 0, 0, t446, t398, t400 * t427 + (-t399 * t423 - t445) * t431 + t434 * qJD(4), (t399 * t458 + t400 * t431 + t427 * t445) * t426 - t398 * t430 + (t434 * t430 + (-t405 * t419 + t413 * t423) * t426) * qJD(5) + (-t406 * t427 + t442 * t431) * t449;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end