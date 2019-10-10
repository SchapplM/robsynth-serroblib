% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% JgD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S7RRRRRRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobigD_rot_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t71 = sin(qJ(1));
	t76 = qJD(1) * t71;
	t73 = cos(qJ(1));
	t75 = qJD(1) * t73;
	t74 = qJD(2) * cos(qJ(2));
	t70 = sin(qJ(2));
	t1 = [0, t75, t70 * t76 - t73 * t74, 0, 0, 0, 0; 0, t76, -t70 * t75 - t71 * t74, 0, 0, 0, 0; 0, 0, -qJD(2) * t70, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (41->23), div. (0->0), fcn. (41->6), ass. (0->14)
	t115 = sin(qJ(1));
	t125 = qJD(1) * t115;
	t118 = cos(qJ(1));
	t124 = qJD(1) * t118;
	t114 = sin(qJ(2));
	t123 = qJD(2) * t114;
	t117 = cos(qJ(2));
	t122 = qJD(2) * t117;
	t121 = qJD(2) * t118;
	t120 = qJD(1) * t117 + qJD(3);
	t116 = cos(qJ(3));
	t119 = (-qJD(3) * t117 - qJD(1)) * t116;
	t113 = sin(qJ(3));
	t1 = [0, t124, t114 * t125 - t117 * t121, t118 * t119 + (t114 * t121 + t120 * t115) * t113, 0, 0, 0; 0, t125, -t114 * t124 - t115 * t122, t115 * t119 + (t115 * t123 - t120 * t118) * t113, 0, 0, 0; 0, 0, -t123, -t114 * qJD(3) * t116 - t113 * t122, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (34->26), mult. (110->52), div. (0->0), fcn. (112->8), ass. (0->24)
	t188 = sin(qJ(4));
	t193 = cos(qJ(3));
	t194 = cos(qJ(2));
	t199 = qJD(1) * t194 + qJD(3);
	t212 = t199 * t188 * t193;
	t189 = sin(qJ(3));
	t200 = -qJD(3) * t194 - qJD(1);
	t211 = t200 * t189;
	t209 = t193 * t194;
	t191 = sin(qJ(1));
	t208 = qJD(1) * t191;
	t195 = cos(qJ(1));
	t207 = qJD(1) * t195;
	t190 = sin(qJ(2));
	t206 = qJD(2) * t190;
	t205 = qJD(2) * t193;
	t204 = qJD(2) * t194;
	t203 = qJD(2) * t195;
	t202 = qJD(3) * t190;
	t198 = -qJD(4) + t205;
	t197 = t200 * t193;
	t196 = t190 * t208 - t194 * t203;
	t192 = cos(qJ(4));
	t1 = [0, t207, t196, t195 * t197 + (t190 * t203 + t199 * t191) * t189, -t191 * t212 + (-t198 * t190 + t211) * t188 * t195 + ((-t191 * t189 + t195 * t209) * qJD(4) + t196) * t192, 0, 0; 0, t208, -t190 * t207 - t191 * t204, t191 * t197 + (t191 * t206 - t199 * t195) * t189, (t212 + (-qJD(1) * t190 + t189 * qJD(4)) * t192) * t195 + ((-t190 * t205 + t211) * t188 - t192 * t204 + (t190 * t188 + t192 * t209) * qJD(4)) * t191, 0, 0; 0, 0, -t206, -t189 * t204 - t193 * t202, (qJD(4) * t193 - qJD(2)) * t192 * t190 + (-t189 * t202 + t198 * t194) * t188, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:06
	% EndTime: 2019-10-10 17:10:06
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (81->44), mult. (250->87), div. (0->0), fcn. (261->10), ass. (0->43)
	t279 = cos(qJ(3));
	t311 = -qJD(4) * t279 + qJD(2);
	t275 = sin(qJ(2));
	t276 = sin(qJ(1));
	t280 = cos(qJ(2));
	t290 = qJD(1) * t280 + qJD(3);
	t281 = cos(qJ(1));
	t297 = qJD(2) * t281;
	t310 = t275 * t297 + t290 * t276;
	t272 = sin(qJ(5));
	t273 = sin(qJ(4));
	t309 = t272 * t273;
	t278 = cos(qJ(4));
	t308 = t272 * t278;
	t307 = t273 * t275;
	t306 = t276 * t279;
	t305 = t276 * t280;
	t274 = sin(qJ(3));
	t304 = t281 * t274;
	t303 = t281 * t279;
	t302 = qJD(1) * t276;
	t301 = qJD(1) * t281;
	t300 = qJD(2) * t275;
	t299 = qJD(2) * t279;
	t298 = qJD(2) * t280;
	t296 = qJD(3) * t274;
	t295 = qJD(3) * t279;
	t294 = qJD(4) * t275;
	t291 = -qJD(3) * t280 - qJD(1);
	t289 = -qJD(4) + t299;
	t286 = t291 * t281;
	t288 = t274 * t286 - t310 * t279 + t281 * t294;
	t287 = t290 * t303 + (t291 * t274 - t275 * t299 + t294) * t276;
	t285 = -t275 * t301 - t276 * t298;
	t284 = t275 * t302 - t280 * t297;
	t270 = t279 * t305 + t304;
	t283 = -qJD(4) * t270 - t285;
	t271 = -t276 * t274 + t280 * t303;
	t282 = qJD(4) * t271 + t284;
	t277 = cos(qJ(5));
	t268 = t291 * t306 + (t276 * t300 - t290 * t281) * t274;
	t266 = t310 * t274 + t279 * t286;
	t1 = [0, t301, t284, t266, t288 * t273 + t282 * t278, -t266 * t277 + t288 * t308 - t282 * t309 + ((t271 * t278 + t281 * t307) * t277 + (-t280 * t304 - t306) * t272) * qJD(5), 0; 0, t302, t285, t268, t287 * t273 - t283 * t278, -t268 * t277 + t287 * t308 + t283 * t309 + ((t270 * t278 + t276 * t307) * t277 + (-t274 * t305 + t303) * t272) * qJD(5), 0; 0, 0, -t300, -t274 * t298 - t275 * t295, -t311 * t278 * t275 + (-t275 * t296 + t289 * t280) * t273, (t289 * t308 + (qJD(2) * t274 - t273 * qJD(5)) * t277) * t280 + ((t311 * t273 - t278 * t296) * t272 + t277 * t295 + (t278 * t279 * t277 - t274 * t272) * qJD(5)) * t275, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:08
	% EndTime: 2019-10-10 17:10:09
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (174->63), mult. (520->121), div. (0->0), fcn. (557->12), ass. (0->60)
	t400 = sin(qJ(2));
	t401 = sin(qJ(1));
	t406 = cos(qJ(2));
	t422 = qJD(1) * t406 + qJD(3);
	t407 = cos(qJ(1));
	t427 = qJD(2) * t407;
	t441 = t400 * t427 + t422 * t401;
	t399 = sin(qJ(3));
	t440 = t399 * t400;
	t404 = cos(qJ(4));
	t439 = t400 * t404;
	t405 = cos(qJ(3));
	t438 = t400 * t405;
	t437 = t400 * t407;
	t436 = t401 * t405;
	t435 = t401 * t406;
	t434 = t405 * t407;
	t433 = t406 * t407;
	t432 = qJD(1) * t401;
	t431 = qJD(1) * t407;
	t430 = qJD(2) * t400;
	t429 = qJD(2) * t405;
	t428 = qJD(2) * t406;
	t426 = qJD(3) * t400;
	t425 = qJD(4) * t400;
	t423 = -qJD(3) * t406 - qJD(1);
	t421 = qJD(4) * t405 - qJD(2);
	t415 = t423 * t407;
	t420 = t399 * t415 - t441 * t405 + t407 * t425;
	t419 = t422 * t434 + (t423 * t399 - t400 * t429 + t425) * t401;
	t398 = sin(qJ(4));
	t414 = (-qJD(4) + t429) * t406;
	t418 = -qJD(5) * t440 + t404 * t414 + (-qJD(3) * t399 * t404 - t421 * t398) * t400;
	t393 = t399 * t407 + t405 * t435;
	t389 = t398 * t400 * t401 + t393 * t404;
	t392 = -t399 * t435 + t434;
	t397 = sin(qJ(5));
	t403 = cos(qJ(5));
	t417 = t389 * t403 + t392 * t397;
	t395 = -t401 * t399 + t405 * t433;
	t390 = t395 * t404 + t398 * t437;
	t394 = -t399 * t433 - t436;
	t416 = t390 * t403 + t394 * t397;
	t413 = -t400 * t431 - t401 * t428;
	t412 = t400 * t432 - t406 * t427;
	t411 = -t399 * t428 - t405 * t426;
	t410 = -qJD(4) * t393 - t413;
	t409 = qJD(4) * t395 + t412;
	t391 = -t398 * t406 + t404 * t438;
	t408 = qJD(5) * t391 - t411;
	t402 = cos(qJ(6));
	t396 = sin(qJ(6));
	t387 = t423 * t436 + (t401 * t430 - t422 * t407) * t399;
	t385 = t441 * t399 + t405 * t415;
	t383 = t421 * t439 + (-t399 * t426 + t414) * t398;
	t382 = t410 * t398 + t419 * t404;
	t381 = t419 * t398 - t410 * t404;
	t380 = -t409 * t398 + t420 * t404;
	t379 = t420 * t398 + t409 * t404;
	t1 = [0, t431, t412, t385, t379, t416 * qJD(5) + t380 * t397 - t385 * t403, -(t380 * t403 + t385 * t397 + (-t390 * t397 + t394 * t403) * qJD(5)) * t396 + t379 * t402 + (-t416 * t402 - (t395 * t398 - t404 * t437) * t396) * qJD(6); 0, t432, t413, t387, t381, t417 * qJD(5) + t382 * t397 - t387 * t403, -(t382 * t403 + t387 * t397 + (-t389 * t397 + t392 * t403) * qJD(5)) * t396 + t381 * t402 + (-t417 * t402 - (t393 * t398 - t401 * t439) * t396) * qJD(6); 0, 0, -t430, t411, t383, t418 * t397 + t408 * t403, (t383 - (t391 * t403 - t397 * t440) * qJD(6)) * t402 + (-t418 * t403 + t408 * t397 + (-t398 * t438 - t404 * t406) * qJD(6)) * t396;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,7);
end