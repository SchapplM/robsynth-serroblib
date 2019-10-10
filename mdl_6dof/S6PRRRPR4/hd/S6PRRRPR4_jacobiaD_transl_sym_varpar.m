% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:16
	% EndTime: 2019-10-09 22:52:16
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:52:17
	% EndTime: 2019-10-09 22:52:17
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (224->70), mult. (725->134), div. (0->0), fcn. (724->10), ass. (0->50)
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t290 = sin(qJ(4));
	t293 = cos(qJ(4));
	t305 = t293 * r_i_i_C(1) - t290 * r_i_i_C(2);
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
	t287 = sin(pkin(11));
	t307 = t287 * t310;
	t289 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:52:17
	% EndTime: 2019-10-09 22:52:17
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (385->92), mult. (978->164), div. (0->0), fcn. (980->12), ass. (0->56)
	t305 = sin(qJ(3));
	t308 = cos(qJ(3));
	t300 = qJ(4) + pkin(12);
	t298 = sin(t300);
	t299 = cos(t300);
	t307 = cos(qJ(4));
	t347 = pkin(4) * t307 + r_i_i_C(1) * t299 - r_i_i_C(2) * t298;
	t318 = pkin(3) + t347;
	t341 = r_i_i_C(3) + qJ(5) + pkin(9);
	t348 = -(t318 * t305 - t341 * t308) * qJD(3) + t305 * qJD(5);
	t301 = sin(pkin(11));
	t306 = sin(qJ(2));
	t309 = cos(qJ(2));
	t339 = cos(pkin(11));
	t340 = cos(pkin(6));
	t319 = t340 * t339;
	t288 = t301 * t309 + t306 * t319;
	t302 = sin(pkin(6));
	t325 = t302 * t339;
	t278 = t288 * t308 - t305 * t325;
	t326 = t301 * t340;
	t290 = -t306 * t326 + t339 * t309;
	t304 = sin(qJ(4));
	t315 = pkin(4) * t304 + r_i_i_C(1) * t298 + r_i_i_C(2) * t299;
	t337 = t302 * t305;
	t336 = t302 * t308;
	t335 = t302 * t309;
	t334 = qJD(2) * t306;
	t333 = qJD(4) * t298;
	t332 = qJD(4) * t299;
	t331 = qJD(4) * t307;
	t330 = qJD(4) * t308;
	t328 = t302 * t334;
	t327 = qJD(2) * t335;
	t317 = t309 * t319;
	t316 = -t290 * t305 + t301 * t336;
	t280 = t290 * t308 + t301 * t337;
	t314 = -t288 * t305 - t308 * t325;
	t313 = -t306 * t337 + t340 * t308;
	t292 = t340 * t305 + t306 * t336;
	t312 = qJD(4) * t315;
	t289 = t339 * t306 + t309 * t326;
	t311 = -t341 * t305 - t318 * t308 - pkin(2);
	t310 = t315 * t330 - t348;
	t287 = t301 * t306 - t317;
	t286 = t290 * qJD(2);
	t285 = t289 * qJD(2);
	t284 = t288 * qJD(2);
	t283 = -qJD(2) * t317 + t301 * t334;
	t282 = t313 * qJD(3) + t308 * t327;
	t281 = t292 * qJD(3) + t305 * t327;
	t276 = t316 * qJD(3) - t285 * t308;
	t275 = t280 * qJD(3) - t285 * t305;
	t274 = t314 * qJD(3) - t283 * t308;
	t273 = t278 * qJD(3) - t283 * t305;
	t1 = [0, (-t285 * t298 + t290 * t332) * r_i_i_C(1) + (-t285 * t299 - t290 * t333) * r_i_i_C(2) - t285 * pkin(8) + (-t285 * t304 + t290 * t331) * pkin(4) + t311 * t286 + t310 * t289, qJD(5) * t280 - t318 * t275 + t341 * t276 - t316 * t312, (-t276 * t298 + t286 * t299) * r_i_i_C(1) + (-t276 * t299 - t286 * t298) * r_i_i_C(2) + ((-t280 * t299 - t289 * t298) * r_i_i_C(1) + (t280 * t298 - t289 * t299) * r_i_i_C(2)) * qJD(4) + (-t276 * t304 + t286 * t307 + (-t280 * t307 - t289 * t304) * qJD(4)) * pkin(4), t275, 0; 0, (-t283 * t298 + t288 * t332) * r_i_i_C(1) + (-t283 * t299 - t288 * t333) * r_i_i_C(2) - t283 * pkin(8) + (-t283 * t304 + t288 * t331) * pkin(4) + t311 * t284 + t310 * t287, qJD(5) * t278 - t318 * t273 + t341 * t274 - t314 * t312, (-t274 * t298 + t284 * t299) * r_i_i_C(1) + (-t274 * t299 - t284 * t298) * r_i_i_C(2) + ((-t278 * t299 - t287 * t298) * r_i_i_C(1) + (t278 * t298 - t287 * t299) * r_i_i_C(2)) * qJD(4) + (-t274 * t304 + t284 * t307 + (-t278 * t307 - t287 * t304) * qJD(4)) * pkin(4), t273, 0; 0, ((t311 * qJD(2) + t347 * qJD(4)) * t306 + (qJD(2) * pkin(8) + t315 * (qJD(2) - t330) + t348) * t309) * t302, qJD(5) * t292 - t318 * t281 + t341 * t282 - t313 * t312, (-t282 * t298 + t299 * t328) * r_i_i_C(1) + (-t282 * t299 - t298 * t328) * r_i_i_C(2) + ((-t292 * t299 + t298 * t335) * r_i_i_C(1) + (t292 * t298 + t299 * t335) * r_i_i_C(2)) * qJD(4) + (t307 * t328 - t282 * t304 + (-t292 * t307 + t304 * t335) * qJD(4)) * pkin(4), t281, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:17
	% EndTime: 2019-10-09 22:52:17
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (675->79), mult. (1208->127), div. (0->0), fcn. (1207->14), ass. (0->65)
	t350 = sin(qJ(3));
	t353 = cos(qJ(3));
	t346 = qJ(4) + pkin(12);
	t332 = pkin(5) * cos(t346) + cos(qJ(4)) * pkin(4);
	t342 = qJ(6) + t346;
	t338 = sin(t342);
	t339 = cos(t342);
	t369 = t339 * r_i_i_C(1) - t338 * r_i_i_C(2);
	t365 = pkin(3) + t332 + t369;
	t391 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9);
	t331 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t346);
	t328 = t331 * qJD(4);
	t345 = qJD(4) + qJD(6);
	t368 = t338 * r_i_i_C(1) + t339 * r_i_i_C(2);
	t393 = t368 * t345 + t328;
	t355 = t393 * t353 + (t365 * t350 - t391 * t353) * qJD(3) - t350 * qJD(5);
	t347 = sin(pkin(11));
	t351 = sin(qJ(2));
	t354 = cos(qJ(2));
	t389 = cos(pkin(11));
	t390 = cos(pkin(6));
	t366 = t390 * t389;
	t323 = t347 * t354 + t351 * t366;
	t348 = sin(pkin(6));
	t377 = t348 * t389;
	t313 = t323 * t353 - t350 * t377;
	t378 = t347 * t390;
	t325 = -t351 * t378 + t389 * t354;
	t387 = t348 * t350;
	t386 = t348 * t353;
	t385 = t348 * t354;
	t319 = t323 * qJD(2);
	t373 = t313 * t345 - t319;
	t364 = t354 * t366;
	t381 = qJD(2) * t351;
	t318 = -qJD(2) * t364 + t347 * t381;
	t360 = -t323 * t350 - t353 * t377;
	t309 = t360 * qJD(3) - t318 * t353;
	t322 = t347 * t351 - t364;
	t375 = -t322 * t345 - t309;
	t384 = (t375 * t338 - t373 * t339) * r_i_i_C(1) + (t373 * t338 + t375 * t339) * r_i_i_C(2);
	t315 = t325 * t353 + t347 * t387;
	t321 = t325 * qJD(2);
	t372 = t315 * t345 - t321;
	t324 = t389 * t351 + t354 * t378;
	t320 = t324 * qJD(2);
	t363 = -t325 * t350 + t347 * t386;
	t311 = t363 * qJD(3) - t320 * t353;
	t374 = -t324 * t345 - t311;
	t383 = (t374 * t338 - t372 * t339) * r_i_i_C(1) + (t372 * t338 + t374 * t339) * r_i_i_C(2);
	t327 = t390 * t350 + t351 * t386;
	t362 = -t327 * t345 + t348 * t381;
	t359 = -t351 * t387 + t390 * t353;
	t379 = qJD(2) * t385;
	t317 = t359 * qJD(3) + t353 * t379;
	t367 = t345 * t385 - t317;
	t382 = (t367 * t338 + t362 * t339) * r_i_i_C(1) + (-t362 * t338 + t367 * t339) * r_i_i_C(2);
	t361 = pkin(8) + t331 + t368;
	t329 = t332 * qJD(4);
	t357 = t369 * t345 + t329;
	t356 = -t391 * t350 - t365 * t353 - pkin(2);
	t316 = t327 * qJD(3) + t350 * t379;
	t310 = t315 * qJD(3) - t320 * t350;
	t308 = t313 * qJD(3) - t318 * t350;
	t1 = [0, -t361 * t320 + t356 * t321 + t355 * t324 + t357 * t325, t315 * qJD(5) - t365 * t310 + t391 * t311 - t363 * t393, -t311 * t331 - t315 * t329 + t321 * t332 - t324 * t328 + t383, t310, t383; 0, -t361 * t318 + t356 * t319 + t355 * t322 + t357 * t323, t313 * qJD(5) - t365 * t308 + t391 * t309 - t360 * t393, -t309 * t331 - t313 * t329 + t319 * t332 - t322 * t328 + t384, t308, t384; 0, ((t356 * qJD(2) + t357) * t351 + (t361 * qJD(2) - t355) * t354) * t348, t327 * qJD(5) - t365 * t316 + t391 * t317 - t359 * t393, -t317 * t331 - t327 * t329 + (t354 * t328 + t332 * t381) * t348 + t382, t316, t382;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end