% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR12_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t127 = sin(pkin(6)) * cos(pkin(7));
	t126 = cos(pkin(6)) * cos(pkin(14));
	t125 = cos(qJ(1));
	t124 = sin(qJ(1));
	t119 = sin(pkin(7));
	t118 = sin(pkin(14));
	t1 = [0, 0, (-(t118 * t124 - t125 * t126) * t119 + t125 * t127) * qJD(1), 0, 0, 0; 0, 0, (-(-t118 * t125 - t124 * t126) * t119 + t124 * t127) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (31->21), mult. (115->46), div. (0->0), fcn. (119->12), ass. (0->32)
	t220 = sin(pkin(7));
	t221 = sin(pkin(6));
	t244 = t221 * t220;
	t218 = sin(pkin(14));
	t227 = sin(qJ(1));
	t243 = t227 * t218;
	t222 = cos(pkin(14));
	t242 = t227 * t222;
	t229 = cos(qJ(1));
	t241 = t229 * t218;
	t240 = t229 * t222;
	t219 = sin(pkin(8));
	t239 = qJD(1) * t219;
	t238 = qJD(3) * t219;
	t224 = cos(pkin(7));
	t228 = cos(qJ(3));
	t237 = t219 * t224 * t228;
	t236 = t227 * t244;
	t235 = t229 * t244;
	t234 = qJD(1) * t221 * t224;
	t225 = cos(pkin(6));
	t233 = t225 * t240 - t243;
	t232 = -t225 * t242 - t241;
	t231 = t225 * t241 + t242;
	t230 = t225 * t243 - t240;
	t226 = sin(qJ(3));
	t223 = cos(pkin(8));
	t217 = t232 * qJD(1);
	t216 = t233 * qJD(1);
	t215 = -t217 * t220 + t227 * t234;
	t214 = t216 * t220 + t229 * t234;
	t1 = [0, 0, t214, t216 * t237 + t214 * t223 - (t230 * t228 + (-t232 * t224 - t236) * t226) * t238 - (t231 * t226 + t228 * t235) * t239, 0, 0; 0, 0, t215, -t217 * t237 + t215 * t223 - (-t231 * t228 + (-t233 * t224 + t235) * t226) * t238 - (t230 * t226 + t228 * t236) * t239, 0, 0; 0, 0, 0, -(-t220 * t225 * t226 + (-t222 * t224 * t226 - t218 * t228) * t221) * t238, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:48
	% EndTime: 2019-10-10 09:15:48
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (114->45), mult. (395->99), div. (0->0), fcn. (432->14), ass. (0->50)
	t311 = sin(pkin(14));
	t314 = sin(pkin(6));
	t320 = sin(qJ(3));
	t323 = cos(qJ(3));
	t315 = cos(pkin(14));
	t317 = cos(pkin(7));
	t342 = t315 * t317;
	t313 = sin(pkin(7));
	t318 = cos(pkin(6));
	t345 = t313 * t318;
	t350 = (t311 * t323 + t320 * t342) * t314 + t320 * t345;
	t324 = cos(qJ(1));
	t338 = t324 * t315;
	t321 = sin(qJ(1));
	t341 = t321 * t311;
	t310 = -t318 * t341 + t338;
	t339 = t324 * t311;
	t340 = t321 * t315;
	t309 = -t318 * t340 - t339;
	t344 = t314 * t321;
	t331 = t309 * t317 + t313 * t344;
	t349 = t310 * t323 + t320 * t331;
	t308 = t318 * t339 + t340;
	t307 = t318 * t338 - t341;
	t343 = t314 * t324;
	t332 = -t307 * t317 + t313 * t343;
	t348 = -t308 * t323 + t320 * t332;
	t337 = qJD(1) * t314;
	t319 = sin(qJ(4));
	t336 = qJD(3) * t319;
	t334 = t321 * t337;
	t333 = t324 * t337;
	t303 = t307 * qJD(1);
	t329 = -t303 * t317 + t313 * t333;
	t305 = t309 * qJD(1);
	t328 = t305 * t317 + t313 * t334;
	t327 = -t308 * t320 - t323 * t332;
	t326 = -t310 * t320 + t323 * t331;
	t325 = t323 * t345 + (-t311 * t320 + t323 * t342) * t314;
	t322 = cos(qJ(4));
	t316 = cos(pkin(8));
	t312 = sin(pkin(8));
	t306 = t310 * qJD(1);
	t304 = t308 * qJD(1);
	t302 = -t305 * t313 + t317 * t334;
	t301 = t303 * t313 + t317 * t333;
	t300 = t350 * qJD(3);
	t299 = qJD(3) * t348 - t306 * t320 + t328 * t323;
	t298 = -qJD(3) * t349 + t304 * t320 + t329 * t323;
	t1 = [0, 0, t301, -t298 * t312 + t301 * t316, (-t304 * t323 + t329 * t320) * t319 + (-t298 * t316 - t301 * t312) * t322 + t326 * t336 + (t349 * t322 + (t326 * t316 + (-t309 * t313 + t317 * t344) * t312) * t319) * qJD(4), 0; 0, 0, t302, -t299 * t312 + t302 * t316, (t306 * t323 + t328 * t320) * t319 + (-t299 * t316 - t302 * t312) * t322 + t327 * t336 + (-t348 * t322 + (t327 * t316 + (-t307 * t313 - t317 * t343) * t312) * t319) * qJD(4), 0; 0, 0, 0, t300 * t312, t300 * t316 * t322 + (t350 * t322 + (t325 * t316 + (-t313 * t314 * t315 + t317 * t318) * t312) * t319) * qJD(4) + t325 * t336, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:50
	% EndTime: 2019-10-10 09:15:51
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (281->68), mult. (934->140), div. (0->0), fcn. (1057->16), ass. (0->69)
	t416 = sin(pkin(14));
	t419 = sin(pkin(6));
	t426 = sin(qJ(3));
	t430 = cos(qJ(3));
	t420 = cos(pkin(14));
	t422 = cos(pkin(7));
	t453 = t420 * t422;
	t418 = sin(pkin(7));
	t423 = cos(pkin(6));
	t456 = t418 * t423;
	t404 = (t416 * t430 + t426 * t453) * t419 + t426 * t456;
	t431 = cos(qJ(1));
	t448 = t431 * t420;
	t427 = sin(qJ(1));
	t451 = t427 * t416;
	t415 = -t423 * t451 + t448;
	t449 = t431 * t416;
	t450 = t427 * t420;
	t414 = -t423 * t450 - t449;
	t455 = t419 * t427;
	t438 = t414 * t422 + t418 * t455;
	t398 = t415 * t430 + t438 * t426;
	t413 = t423 * t449 + t450;
	t412 = t423 * t448 - t451;
	t454 = t419 * t431;
	t439 = -t412 * t422 + t418 * t454;
	t461 = -t413 * t430 + t439 * t426;
	t400 = t404 * qJD(3);
	t417 = sin(pkin(8));
	t460 = t400 * t417;
	t425 = sin(qJ(4));
	t457 = t417 * t425;
	t421 = cos(pkin(8));
	t452 = t421 * t425;
	t447 = qJD(1) * t419;
	t424 = sin(qJ(5));
	t446 = qJD(4) * t424;
	t444 = t427 * t447;
	t443 = t431 * t447;
	t395 = -t413 * t426 - t439 * t430;
	t405 = -t412 * t418 - t422 * t454;
	t442 = t395 * t421 + t405 * t417;
	t397 = -t415 * t426 + t438 * t430;
	t406 = -t414 * t418 + t422 * t455;
	t441 = t397 * t421 + t406 * t417;
	t403 = t430 * t456 + (-t416 * t426 + t430 * t453) * t419;
	t411 = -t419 * t420 * t418 + t423 * t422;
	t440 = t403 * t421 + t411 * t417;
	t407 = t412 * qJD(1);
	t436 = -t407 * t422 + t418 * t443;
	t409 = t414 * qJD(1);
	t435 = t409 * t422 + t418 * t444;
	t429 = cos(qJ(4));
	t434 = t442 * t425 - t429 * t461;
	t433 = t398 * t429 + t441 * t425;
	t432 = t404 * t429 + t440 * t425;
	t428 = cos(qJ(5));
	t410 = t415 * qJD(1);
	t408 = t413 * qJD(1);
	t402 = -t409 * t418 + t422 * t444;
	t401 = t407 * t418 + t422 * t443;
	t399 = t403 * qJD(3);
	t394 = t395 * qJD(3) + t410 * t430 + t435 * t426;
	t393 = t461 * qJD(3) - t410 * t426 + t435 * t430;
	t392 = t397 * qJD(3) - t408 * t430 + t436 * t426;
	t391 = -t398 * qJD(3) + t408 * t426 + t436 * t430;
	t390 = -t393 * t417 + t402 * t421;
	t389 = -t391 * t417 + t401 * t421;
	t1 = [0, 0, t401, t389, t392 * t425 + (-t391 * t421 - t401 * t417) * t429 + t433 * qJD(4), (t391 * t452 + t392 * t429 + t401 * t457) * t424 - t389 * t428 + (t433 * t428 + (-t397 * t417 + t406 * t421) * t424) * qJD(5) + (-t398 * t425 + t441 * t429) * t446; 0, 0, t402, t390, t394 * t425 + (-t393 * t421 - t402 * t417) * t429 + t434 * qJD(4), (t393 * t452 + t394 * t429 + t402 * t457) * t424 - t390 * t428 + (t434 * t428 + (-t395 * t417 + t405 * t421) * t424) * qJD(5) + (t425 * t461 + t442 * t429) * t446; 0, 0, 0, t460, t400 * t421 * t429 + t432 * qJD(4) + t399 * t425, (t399 * t429 - t400 * t452) * t424 - t428 * t460 + (t432 * t428 + (-t403 * t417 + t411 * t421) * t424) * qJD(5) + (-t404 * t425 + t440 * t429) * t446;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end