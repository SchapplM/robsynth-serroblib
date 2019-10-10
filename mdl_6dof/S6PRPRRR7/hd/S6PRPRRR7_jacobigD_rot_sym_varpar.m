% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR7_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:03
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->10), mult. (54->29), div. (0->0), fcn. (54->12), ass. (0->15)
	t205 = sin(pkin(7)) * cos(pkin(8));
	t204 = cos(pkin(7)) * cos(pkin(14));
	t198 = cos(pkin(6));
	t199 = sin(qJ(2));
	t203 = t198 * t199;
	t200 = cos(qJ(2));
	t202 = t198 * t200;
	t189 = sin(pkin(14));
	t201 = qJD(2) * t189;
	t195 = cos(pkin(13));
	t191 = sin(pkin(8));
	t190 = sin(pkin(13));
	t188 = (t190 * t203 - t195 * t200) * qJD(2);
	t187 = (-t190 * t200 - t195 * t203) * qJD(2);
	t1 = [0, 0, 0, -(t188 * t204 - (-t190 * t202 - t195 * t199) * t201) * t191 - t188 * t205, 0, 0; 0, 0, 0, -(t187 * t204 - (-t190 * t199 + t195 * t202) * t201) * t191 - t187 * t205, 0, 0; 0, 0, 0, (-(-t189 * t200 - t199 * t204) * t191 + t199 * t205) * sin(pkin(6)) * qJD(2), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:03
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (68->41), mult. (251->100), div. (0->0), fcn. (276->14), ass. (0->41)
	t285 = sin(pkin(14));
	t293 = cos(pkin(7));
	t314 = t285 * t293;
	t287 = sin(pkin(8));
	t288 = sin(pkin(7));
	t313 = t288 * t287;
	t289 = sin(pkin(6));
	t312 = t288 * t289;
	t292 = cos(pkin(8));
	t311 = t288 * t292;
	t294 = cos(pkin(6));
	t310 = t288 * t294;
	t296 = sin(qJ(2));
	t309 = t288 * t296;
	t308 = t289 * t293;
	t290 = cos(pkin(14));
	t307 = t290 * t293;
	t306 = t293 * t296;
	t298 = cos(qJ(2));
	t305 = t293 * t298;
	t304 = t294 * t296;
	t303 = t294 * t298;
	t302 = qJD(2) * t289;
	t286 = sin(pkin(13));
	t291 = cos(pkin(13));
	t281 = -t286 * t296 + t291 * t303;
	t301 = t281 * t293 - t291 * t312;
	t283 = -t286 * t303 - t291 * t296;
	t300 = t283 * t293 + t286 * t312;
	t282 = t286 * t298 + t291 * t304;
	t299 = t286 * t304 - t291 * t298;
	t297 = cos(qJ(4));
	t295 = sin(qJ(4));
	t280 = t299 * qJD(2);
	t279 = t283 * qJD(2);
	t278 = t282 * qJD(2);
	t277 = t281 * qJD(2);
	t276 = (-t285 * t298 - t290 * t306) * t302;
	t275 = -t279 * t285 + t280 * t307;
	t274 = -t277 * t285 - t278 * t307;
	t1 = [0, 0, 0, -t275 * t287 - t280 * t311, (t279 * t290 + t280 * t314) * t295 + (-t275 * t292 + t280 * t313) * t297 + ((t300 * t285 - t290 * t299) * t297 + ((t285 * t299 + t300 * t290) * t292 + (-t283 * t288 + t286 * t308) * t287) * t295) * qJD(4), 0; 0, 0, 0, -t274 * t287 + t278 * t311, (t277 * t290 - t278 * t314) * t295 + (-t274 * t292 - t278 * t313) * t297 + ((t282 * t290 + t301 * t285) * t297 + ((-t282 * t285 + t301 * t290) * t292 + (-t281 * t288 - t291 * t308) * t287) * t295) * qJD(4), 0; 0, 0, 0, t292 * t302 * t309 - t276 * t287, -t276 * t292 * t297 + ((t289 * t296 * t290 + (t289 * t305 + t310) * t285) * t297 + ((t290 * t310 + (-t285 * t296 + t290 * t305) * t289) * t292 + (t294 * t293 - t298 * t312) * t287) * t295) * qJD(4) + ((-t285 * t306 + t290 * t298) * t295 - t287 * t297 * t309) * t302, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:04
	% EndTime: 2019-10-09 22:05:05
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (194->65), mult. (664->140), div. (0->0), fcn. (763->16), ass. (0->68)
	t384 = sin(pkin(14));
	t392 = cos(pkin(7));
	t425 = t384 * t392;
	t386 = sin(pkin(8));
	t387 = sin(pkin(7));
	t424 = t386 * t387;
	t388 = sin(pkin(6));
	t423 = t387 * t388;
	t391 = cos(pkin(8));
	t422 = t387 * t391;
	t393 = cos(pkin(6));
	t421 = t387 * t393;
	t420 = t388 * t392;
	t389 = cos(pkin(14));
	t419 = t389 * t392;
	t395 = sin(qJ(4));
	t418 = t391 * t395;
	t396 = sin(qJ(2));
	t417 = t392 * t396;
	t399 = cos(qJ(2));
	t416 = t392 * t399;
	t415 = t393 * t396;
	t414 = t393 * t399;
	t413 = qJD(2) * t388;
	t394 = sin(qJ(5));
	t412 = qJD(4) * t394;
	t411 = t395 * t424;
	t410 = t387 * t396 * t413;
	t409 = t386 * t410;
	t385 = sin(pkin(13));
	t390 = cos(pkin(13));
	t381 = t385 * t399 + t390 * t415;
	t380 = -t385 * t396 + t390 * t414;
	t405 = t380 * t392 - t390 * t423;
	t360 = -t381 * t384 + t405 * t389;
	t371 = -t380 * t387 - t390 * t420;
	t408 = t360 * t391 + t371 * t386;
	t403 = t385 * t415 - t390 * t399;
	t382 = -t385 * t414 - t390 * t396;
	t404 = t382 * t392 + t385 * t423;
	t362 = t384 * t403 + t404 * t389;
	t372 = -t382 * t387 + t385 * t420;
	t407 = t362 * t391 + t372 * t386;
	t369 = t389 * t421 + (-t384 * t396 + t389 * t416) * t388;
	t379 = t393 * t392 - t399 * t423;
	t406 = t369 * t391 + t379 * t386;
	t361 = t381 * t389 + t405 * t384;
	t398 = cos(qJ(4));
	t402 = t361 * t398 + t408 * t395;
	t363 = t404 * t384 - t389 * t403;
	t401 = t363 * t398 + t407 * t395;
	t370 = t388 * t396 * t389 + (t388 * t416 + t421) * t384;
	t400 = t370 * t398 + t406 * t395;
	t397 = cos(qJ(5));
	t378 = t403 * qJD(2);
	t377 = t382 * qJD(2);
	t376 = t381 * qJD(2);
	t375 = t380 * qJD(2);
	t374 = (-t384 * t417 + t389 * t399) * t413;
	t373 = (-t384 * t399 - t389 * t417) * t413;
	t368 = -t373 * t386 + t391 * t410;
	t367 = t377 * t389 + t378 * t425;
	t366 = -t377 * t384 + t378 * t419;
	t365 = t375 * t389 - t376 * t425;
	t364 = -t375 * t384 - t376 * t419;
	t359 = -t366 * t386 - t378 * t422;
	t358 = -t364 * t386 + t376 * t422;
	t1 = [0, 0, 0, t359, t367 * t395 + (-t366 * t391 + t378 * t424) * t398 + t401 * qJD(4), (t366 * t418 + t367 * t398 - t378 * t411) * t394 - t359 * t397 + (t401 * t397 + (-t362 * t386 + t372 * t391) * t394) * qJD(5) + (-t363 * t395 + t407 * t398) * t412; 0, 0, 0, t358, t365 * t395 + (-t364 * t391 - t376 * t424) * t398 + t402 * qJD(4), (t364 * t418 + t365 * t398 + t376 * t411) * t394 - t358 * t397 + (t402 * t397 + (-t360 * t386 + t371 * t391) * t394) * qJD(5) + (-t361 * t395 + t408 * t398) * t412; 0, 0, 0, t368, t374 * t395 + (-t373 * t391 - t409) * t398 + t400 * qJD(4), (t373 * t418 + t374 * t398 + t395 * t409) * t394 - t368 * t397 + (t400 * t397 + (-t369 * t386 + t379 * t391) * t394) * qJD(5) + (-t370 * t395 + t406 * t398) * t412;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end