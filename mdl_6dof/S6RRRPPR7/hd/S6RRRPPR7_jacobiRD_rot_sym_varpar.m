% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:46
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
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:46
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
	% StartTime: 2019-10-10 11:27:47
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(3));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(3));
	t249 = qJD(3) * t229;
	t248 = qJD(3) * t230;
	t247 = qJD(3) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(3) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0, 0, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0, 0, 0; t227, t236 * t229 + t234 * t245, t226, 0, 0, 0; t225, t237 * t229 + t231 * t245, t228, 0, 0, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0, 0; 0, t251, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:47
	% EndTime: 2019-10-10 11:27:48
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t297 = cos(qJ(1));
	t296 = cos(qJ(2));
	t316 = qJD(1) * t296;
	t301 = -qJD(3) + t316;
	t319 = t301 * t297;
	t294 = sin(qJ(1));
	t293 = sin(qJ(2));
	t312 = qJD(2) * t297;
	t304 = t293 * t312;
	t318 = t301 * t294 + t304;
	t317 = qJD(1) * t294;
	t315 = qJD(1) * t297;
	t314 = qJD(2) * t293;
	t313 = qJD(2) * t296;
	t292 = sin(qJ(3));
	t311 = qJD(3) * t292;
	t310 = qJD(3) * t293;
	t309 = qJD(3) * t296;
	t295 = cos(qJ(3));
	t308 = t295 * t314;
	t307 = t295 * t310;
	t306 = t294 * t314;
	t305 = t294 * t313;
	t303 = t296 * t312;
	t302 = qJD(1) - t309;
	t300 = t302 * t297;
	t299 = -t293 * t315 - t305;
	t298 = t293 * t317 - t303;
	t291 = t295 * t319 + (t302 * t292 - t308) * t294;
	t290 = t302 * t295 * t294 + (t306 - t319) * t292;
	t289 = t292 * t300 - t318 * t295;
	t288 = t318 * t292 + t295 * t300;
	t1 = [-t291, -t295 * t303 + (t295 * t317 + t297 * t311) * t293, t288, 0, 0, 0; t289, -t295 * t305 + (t294 * t311 - t295 * t315) * t293, t290, 0, 0, 0; 0, -t292 * t309 - t308, -t292 * t313 - t307, 0, 0, 0; t299, -t294 * t316 - t304, 0, 0, 0, 0; -t298, t296 * t315 - t306, 0, 0, 0, 0; 0, t313, 0, 0, 0, 0; t290, t298 * t292 - t297 * t307, t289, 0, 0, 0; -t288, t299 * t292 - t294 * t307, t291, 0, 0, 0; 0, -t292 * t314 + t295 * t309, -t292 * t310 + t295 * t313, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:48
	% EndTime: 2019-10-10 11:27:48
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (109->39), mult. (385->76), div. (0->0), fcn. (385->8), ass. (0->39)
	t319 = sin(pkin(10));
	t320 = cos(pkin(10));
	t321 = sin(qJ(3));
	t324 = cos(qJ(3));
	t336 = t319 * t324 - t320 * t321;
	t354 = qJD(3) * t336;
	t323 = sin(qJ(1));
	t353 = t321 * t323;
	t326 = cos(qJ(1));
	t352 = t324 * t326;
	t351 = qJD(1) * t323;
	t325 = cos(qJ(2));
	t350 = qJD(1) * t325;
	t349 = qJD(1) * t326;
	t348 = qJD(2) * t325;
	t347 = qJD(2) * t326;
	t322 = sin(qJ(2));
	t346 = t323 * qJD(2) * t322;
	t345 = qJD(3) * t353;
	t344 = t322 * t347;
	t343 = qJD(3) * t352;
	t310 = t321 * t344 - t325 * t343 - t345 + (t325 * t353 + t352) * qJD(1);
	t339 = qJD(3) - t350;
	t340 = -qJD(3) * t325 + qJD(1);
	t311 = t340 * t326 * t321 + (t339 * t323 - t344) * t324;
	t342 = -t310 * t319 + t311 * t320;
	t312 = t340 * t324 * t323 + (t339 * t326 + t346) * t321;
	t313 = -t325 * t345 - t324 * t346 - t343 + (t325 * t352 + t353) * qJD(1);
	t341 = t312 * t320 + t313 * t319;
	t338 = t310 * t320 + t311 * t319;
	t337 = t312 * t319 - t313 * t320;
	t335 = t319 * t321 + t320 * t324;
	t334 = t335 * t325;
	t333 = qJD(2) * t336;
	t332 = qJD(2) * t335;
	t329 = qJD(3) * t335;
	t328 = t325 * t332;
	t327 = t325 * t333;
	t1 = [t337, -t326 * t328 + (-t326 * t354 + t335 * t351) * t322, t338, 0, 0, 0; t342, -t323 * t328 + (-t323 * t354 - t335 * t349) * t322, t341, 0, 0, 0; 0, -t322 * t332 + t325 * t354, -t322 * t329 + t327, 0, 0, 0; t341, t326 * t327 + (-t326 * t329 - t336 * t351) * t322, t342, 0, 0, 0; -t338, t323 * t327 + (-t323 * t329 + t336 * t349) * t322, -t337, 0, 0, 0; 0, qJD(3) * t334 + t322 * t333, qJD(2) * t334 + t322 * t354, 0, 0, 0; t322 * t349 + t323 * t348, t323 * t350 + t344, 0, 0, 0, 0; t322 * t351 - t325 * t347, -t325 * t349 + t346, 0, 0, 0, 0; 0, -t348, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:48
	% EndTime: 2019-10-10 11:27:49
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (335->58), mult. (709->92), div. (0->0), fcn. (733->8), ass. (0->48)
	t383 = pkin(10) + qJ(6);
	t381 = sin(t383);
	t382 = cos(t383);
	t384 = sin(qJ(3));
	t387 = cos(qJ(3));
	t399 = t381 * t384 + t382 * t387;
	t424 = qJD(3) - qJD(6);
	t390 = t424 * t399;
	t400 = t381 * t387 - t382 * t384;
	t425 = t424 * t400;
	t388 = cos(qJ(2));
	t385 = sin(qJ(2));
	t386 = sin(qJ(1));
	t409 = t386 * qJD(2) * t385;
	t389 = cos(qJ(1));
	t414 = qJD(1) * t389;
	t423 = -t388 * t414 + t409;
	t417 = t389 * t387;
	t419 = t386 * t388;
	t367 = t384 * t419 + t417;
	t410 = qJD(3) * t389;
	t406 = t387 * t410;
	t412 = qJD(2) * t389;
	t407 = t385 * t412;
	t411 = qJD(3) * t386;
	t408 = t384 * t411;
	t363 = t367 * qJD(1) + t384 * t407 - t388 * t406 - t408;
	t415 = qJD(1) * t388;
	t418 = t389 * t384;
	t364 = (-qJD(3) * t388 + qJD(1)) * t418 + (-t407 + (qJD(3) - t415) * t386) * t387;
	t369 = -t386 * t387 + t388 * t418;
	t370 = t386 * t384 + t388 * t417;
	t422 = t363 * t382 + t364 * t381 + (t369 * t381 + t370 * t382) * qJD(6);
	t416 = qJD(1) * t386;
	t365 = (t411 * t388 - t416) * t387 + (-t410 - t423) * t384;
	t366 = t370 * qJD(1) - t387 * t409 - t388 * t408 - t406;
	t368 = t387 * t419 - t418;
	t421 = t365 * t381 + t366 * t382 + (t367 * t382 - t368 * t381) * qJD(6);
	t394 = (t367 * t381 + t368 * t382) * qJD(6) - t365 * t382 + t366 * t381;
	t413 = qJD(2) * t388;
	t398 = qJD(2) * t400;
	t397 = qJD(2) * t399;
	t396 = t388 * t397;
	t395 = t388 * t398;
	t359 = -t363 * t381 + t364 * t382 + (t369 * t382 - t370 * t381) * qJD(6);
	t361 = t385 * t425 + t399 * t413;
	t360 = -t390 * t385 + t395;
	t1 = [-t421, -t389 * t396 + (-t389 * t425 + t399 * t416) * t385, t422, 0, 0, -t422; t359, -t386 * t396 + (-t386 * t425 - t399 * t414) * t385, t394, 0, 0, -t394; 0, -t385 * t397 + t388 * t425, t360, 0, 0, -t360; t394, t389 * t395 + (-t390 * t389 - t400 * t416) * t385, t359, 0, 0, -t359; -t422, t386 * t395 + (-t390 * t386 + t400 * t414) * t385, t421, 0, 0, -t421; 0, t385 * t398 + t388 * t390, t361, 0, 0, -t361; t385 * t414 + t386 * t413, t386 * t415 + t407, 0, 0, 0, 0; t385 * t416 - t388 * t412, t423, 0, 0, 0, 0; 0, -t413, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end