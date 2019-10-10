% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(12));
	t50 = sin(pkin(12));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(12));
	t182 = cos(pkin(12));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:30
	% EndTime: 2019-10-09 23:17:31
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (224->70), mult. (725->134), div. (0->0), fcn. (724->10), ass. (0->50)
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t290 = sin(qJ(4));
	t293 = cos(qJ(4));
	t305 = r_i_i_C(1) * t293 - r_i_i_C(2) * t290;
	t303 = pkin(3) + t305;
	t324 = pkin(9) + r_i_i_C(3);
	t326 = (t303 * t291 - t324 * t294) * qJD(3);
	t304 = t290 * r_i_i_C(1) + t293 * r_i_i_C(2);
	t321 = cos(pkin(6));
	t288 = sin(pkin(6));
	t320 = t288 * t291;
	t319 = t288 * t294;
	t295 = cos(qJ(2));
	t318 = t288 * t295;
	t292 = sin(qJ(2));
	t317 = qJD(2) * t292;
	t316 = qJD(2) * t295;
	t315 = qJD(4) * t290;
	t314 = qJD(4) * t293;
	t313 = qJD(4) * t294;
	t312 = t288 * t317;
	t311 = t288 * t316;
	t310 = t292 * t321;
	t309 = t295 * t321;
	t287 = sin(pkin(12));
	t307 = t287 * t310;
	t289 = cos(pkin(12));
	t306 = t289 * t309;
	t279 = t287 * t295 + t289 * t310;
	t302 = -t279 * t291 - t289 * t319;
	t301 = -t279 * t294 + t289 * t320;
	t281 = t289 * t295 - t307;
	t300 = -t281 * t291 + t287 * t319;
	t271 = t281 * t294 + t287 * t320;
	t299 = qJD(4) * t304;
	t280 = t287 * t309 + t289 * t292;
	t298 = -t292 * t320 + t321 * t294;
	t283 = t321 * t291 + t292 * t319;
	t297 = -t324 * t291 - t303 * t294 - pkin(2);
	t296 = t304 * t313 + t326;
	t278 = t287 * t292 - t306;
	t277 = -qJD(2) * t307 + t289 * t316;
	t276 = t280 * qJD(2);
	t275 = t279 * qJD(2);
	t274 = -qJD(2) * t306 + t287 * t317;
	t273 = t298 * qJD(3) + t294 * t311;
	t267 = t300 * qJD(3) - t276 * t294;
	t265 = t302 * qJD(3) - t274 * t294;
	t1 = [0, (-t276 * t290 + t281 * t314) * r_i_i_C(1) + (-t276 * t293 - t281 * t315) * r_i_i_C(2) - t276 * pkin(8) + t297 * t277 + t296 * t280, t324 * t267 - t300 * t299 + t303 * (-t271 * qJD(3) + t276 * t291), (-t267 * t290 + t277 * t293) * r_i_i_C(1) + (-t267 * t293 - t277 * t290) * r_i_i_C(2) + ((-t271 * t293 - t280 * t290) * r_i_i_C(1) + (t271 * t290 - t280 * t293) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (-t274 * t290 + t279 * t314) * r_i_i_C(1) + (-t274 * t293 - t279 * t315) * r_i_i_C(2) - t274 * pkin(8) + t297 * t275 + t296 * t278, t324 * t265 - t302 * t299 + t303 * (t301 * qJD(3) + t274 * t291), (-t265 * t290 + t275 * t293) * r_i_i_C(1) + (-t265 * t293 - t275 * t290) * r_i_i_C(2) + ((-t278 * t290 + t293 * t301) * r_i_i_C(1) + (-t278 * t293 - t290 * t301) * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((t297 * qJD(2) + t305 * qJD(4)) * t292 + (qJD(2) * pkin(8) - t326 + t304 * (qJD(2) - t313)) * t295) * t288, t324 * t273 - t298 * t299 + t303 * (-t283 * qJD(3) - t291 * t311), (-t273 * t290 + t293 * t312) * r_i_i_C(1) + (-t273 * t293 - t290 * t312) * r_i_i_C(2) + ((-t283 * t293 + t290 * t318) * r_i_i_C(1) + (t283 * t290 + t293 * t318) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:30
	% EndTime: 2019-10-09 23:17:31
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (487->86), mult. (1093->149), div. (0->0), fcn. (1098->12), ass. (0->63)
	t334 = sin(qJ(3));
	t337 = cos(qJ(3));
	t336 = cos(qJ(4));
	t330 = qJ(4) + qJ(5);
	t327 = sin(t330);
	t328 = cos(t330);
	t353 = t328 * r_i_i_C(1) - t327 * r_i_i_C(2);
	t349 = t336 * pkin(4) + pkin(3) + t353;
	t378 = r_i_i_C(3) + pkin(10) + pkin(9);
	t384 = (t349 * t334 - t378 * t337) * qJD(3);
	t335 = sin(qJ(2));
	t338 = cos(qJ(2));
	t331 = sin(pkin(12));
	t377 = cos(pkin(6));
	t362 = t331 * t377;
	t376 = cos(pkin(12));
	t320 = -t335 * t362 + t376 * t338;
	t352 = t327 * r_i_i_C(1) + t328 * r_i_i_C(2);
	t329 = qJD(4) + qJD(5);
	t333 = sin(qJ(4));
	t379 = t333 * pkin(4);
	t383 = qJD(4) * t379 + t352 * t329;
	t375 = t327 * t329;
	t374 = t328 * t329;
	t332 = sin(pkin(6));
	t373 = t332 * t334;
	t372 = t332 * t337;
	t371 = t332 * t338;
	t350 = t377 * t376;
	t318 = t331 * t338 + t335 * t350;
	t314 = t318 * qJD(2);
	t361 = t332 * t376;
	t343 = -t318 * t337 + t334 * t361;
	t357 = -t329 * t343 - t314;
	t348 = t338 * t350;
	t367 = qJD(2) * t335;
	t313 = -qJD(2) * t348 + t331 * t367;
	t345 = -t318 * t334 - t337 * t361;
	t304 = t345 * qJD(3) - t313 * t337;
	t317 = t331 * t335 - t348;
	t359 = -t317 * t329 - t304;
	t370 = (t359 * t327 - t357 * t328) * r_i_i_C(1) + (t357 * t327 + t359 * t328) * r_i_i_C(2);
	t310 = t320 * t337 + t331 * t373;
	t316 = t320 * qJD(2);
	t356 = t310 * t329 - t316;
	t319 = t376 * t335 + t338 * t362;
	t315 = t319 * qJD(2);
	t347 = -t320 * t334 + t331 * t372;
	t306 = t347 * qJD(3) - t315 * t337;
	t358 = -t319 * t329 - t306;
	t369 = (t358 * t327 - t356 * t328) * r_i_i_C(1) + (t356 * t327 + t358 * t328) * r_i_i_C(2);
	t322 = t377 * t334 + t335 * t372;
	t364 = t332 * t367;
	t346 = -t322 * t329 + t364;
	t344 = -t335 * t373 + t377 * t337;
	t363 = qJD(2) * t371;
	t312 = t344 * qJD(3) + t337 * t363;
	t351 = t329 * t371 - t312;
	t368 = (t351 * t327 + t346 * t328) * r_i_i_C(1) + (-t346 * t327 + t351 * t328) * r_i_i_C(2);
	t366 = qJD(4) * t336;
	t341 = -t378 * t334 - t349 * t337 - pkin(2);
	t340 = t383 * t337 + t384;
	t1 = [0, (-t315 * t327 + t320 * t374) * r_i_i_C(1) + (-t315 * t328 - t320 * t375) * r_i_i_C(2) - t315 * pkin(8) + (-t315 * t333 + t320 * t366) * pkin(4) + t341 * t316 + t340 * t319, t378 * t306 - t383 * t347 + t349 * (-t310 * qJD(3) + t315 * t334), (-t306 * t333 + t316 * t336 + (-t310 * t336 - t319 * t333) * qJD(4)) * pkin(4) + t369, t369, 0; 0, (-t313 * t327 + t318 * t374) * r_i_i_C(1) + (-t313 * t328 - t318 * t375) * r_i_i_C(2) - t313 * pkin(8) + (-t313 * t333 + t318 * t366) * pkin(4) + t341 * t314 + t340 * t317, t378 * t304 - t383 * t345 + t349 * (t343 * qJD(3) + t313 * t334), (-t304 * t333 + t314 * t336 + (-t317 * t333 + t336 * t343) * qJD(4)) * pkin(4) + t370, t370, 0; 0, ((pkin(4) * t366 + t341 * qJD(2) + t353 * t329) * t335 + (qJD(2) * pkin(8) + (-qJD(4) * t337 + qJD(2)) * t379 - t384 + t352 * (-t329 * t337 + qJD(2))) * t338) * t332, t378 * t312 - t383 * t344 + t349 * (-t322 * qJD(3) - t334 * t363), (t336 * t364 - t312 * t333 + (-t322 * t336 + t333 * t371) * qJD(4)) * pkin(4) + t368, t368, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:30
	% EndTime: 2019-10-09 23:17:31
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (889->90), mult. (1423->142), div. (0->0), fcn. (1422->14), ass. (0->70)
	t350 = sin(qJ(3));
	t353 = cos(qJ(3));
	t346 = qJ(4) + qJ(5);
	t341 = cos(t346);
	t352 = cos(qJ(4));
	t332 = t352 * pkin(4) + pkin(5) * t341;
	t342 = qJ(6) + t346;
	t337 = sin(t342);
	t338 = cos(t342);
	t370 = t338 * r_i_i_C(1) - t337 * r_i_i_C(2);
	t366 = pkin(3) + t332 + t370;
	t391 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9);
	t340 = sin(t346);
	t349 = sin(qJ(4));
	t390 = pkin(4) * qJD(4);
	t344 = qJD(4) + qJD(5);
	t392 = pkin(5) * t344;
	t328 = -t340 * t392 - t349 * t390;
	t339 = qJD(6) + t344;
	t369 = t337 * r_i_i_C(1) + t338 * r_i_i_C(2);
	t394 = t369 * t339 - t328;
	t355 = t394 * t353 + (t366 * t350 - t391 * t353) * qJD(3);
	t351 = sin(qJ(2));
	t354 = cos(qJ(2));
	t347 = sin(pkin(12));
	t389 = cos(pkin(6));
	t378 = t347 * t389;
	t388 = cos(pkin(12));
	t325 = -t351 * t378 + t388 * t354;
	t348 = sin(pkin(6));
	t387 = t348 * t350;
	t386 = t348 * t353;
	t385 = t348 * t354;
	t367 = t389 * t388;
	t323 = t347 * t354 + t351 * t367;
	t319 = t323 * qJD(2);
	t377 = t348 * t388;
	t359 = -t323 * t353 + t350 * t377;
	t373 = -t339 * t359 - t319;
	t365 = t354 * t367;
	t381 = qJD(2) * t351;
	t318 = -qJD(2) * t365 + t347 * t381;
	t361 = -t323 * t350 - t353 * t377;
	t309 = t361 * qJD(3) - t318 * t353;
	t322 = t347 * t351 - t365;
	t375 = -t322 * t339 - t309;
	t384 = (t375 * t337 - t373 * t338) * r_i_i_C(1) + (t373 * t337 + t375 * t338) * r_i_i_C(2);
	t315 = t325 * t353 + t347 * t387;
	t321 = t325 * qJD(2);
	t372 = t315 * t339 - t321;
	t324 = t388 * t351 + t354 * t378;
	t320 = t324 * qJD(2);
	t364 = -t325 * t350 + t347 * t386;
	t311 = t364 * qJD(3) - t320 * t353;
	t374 = -t324 * t339 - t311;
	t383 = (t374 * t337 - t372 * t338) * r_i_i_C(1) + (t372 * t337 + t374 * t338) * r_i_i_C(2);
	t327 = t389 * t350 + t351 * t386;
	t380 = t348 * t381;
	t363 = -t327 * t339 + t380;
	t360 = -t351 * t387 + t389 * t353;
	t379 = qJD(2) * t385;
	t317 = t360 * qJD(3) + t353 * t379;
	t368 = t339 * t385 - t317;
	t382 = (t368 * t337 + t363 * t338) * r_i_i_C(1) + (-t363 * t337 + t368 * t338) * r_i_i_C(2);
	t331 = t349 * pkin(4) + pkin(5) * t340;
	t362 = pkin(8) + t331 + t369;
	t329 = t341 * t392 + t352 * t390;
	t357 = t370 * t339 + t329;
	t356 = -t391 * t350 - t366 * t353 - pkin(2);
	t1 = [0, -t362 * t320 + t356 * t321 + t355 * t324 + t357 * t325, t391 * t311 - t394 * t364 + t366 * (-t315 * qJD(3) + t320 * t350), -t311 * t331 - t315 * t329 + t321 * t332 + t324 * t328 + t383, ((-t315 * t344 + t321) * t341 + (-t324 * t344 - t311) * t340) * pkin(5) + t383, t383; 0, -t362 * t318 + t356 * t319 + t355 * t322 + t357 * t323, t391 * t309 - t394 * t361 + t366 * (t359 * qJD(3) + t318 * t350), -t309 * t331 + t319 * t332 + t322 * t328 + t329 * t359 + t384, ((t344 * t359 + t319) * t341 + (-t322 * t344 - t309) * t340) * pkin(5) + t384, t384; 0, ((t356 * qJD(2) + t357) * t351 + (t362 * qJD(2) - t355) * t354) * t348, t391 * t317 - t394 * t360 + t366 * (-t327 * qJD(3) - t350 * t379), -t317 * t331 - t327 * t329 + (-t354 * t328 + t332 * t381) * t348 + t382, ((-t327 * t344 + t380) * t341 + (t344 * t385 - t317) * t340) * pkin(5) + t382, t382;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end