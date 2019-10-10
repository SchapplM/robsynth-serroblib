% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:09
	% EndTime: 2019-10-09 22:41:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(10));
	t182 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:41:10
	% EndTime: 2019-10-09 22:41:10
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
	t304 = r_i_i_C(1) * t290 + r_i_i_C(2) * t293;
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
	t287 = sin(pkin(10));
	t307 = t287 * t310;
	t289 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:41:10
	% EndTime: 2019-10-09 22:41:10
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (385->92), mult. (978->164), div. (0->0), fcn. (980->12), ass. (0->56)
	t305 = sin(qJ(3));
	t308 = cos(qJ(3));
	t300 = qJ(4) + pkin(11);
	t298 = sin(t300);
	t299 = cos(t300);
	t307 = cos(qJ(4));
	t347 = t307 * pkin(4) + t299 * r_i_i_C(1) - t298 * r_i_i_C(2);
	t318 = pkin(3) + t347;
	t341 = r_i_i_C(3) + qJ(5) + pkin(9);
	t348 = -(t318 * t305 - t341 * t308) * qJD(3) + t305 * qJD(5);
	t301 = sin(pkin(10));
	t306 = sin(qJ(2));
	t309 = cos(qJ(2));
	t339 = cos(pkin(10));
	t340 = cos(pkin(6));
	t319 = t340 * t339;
	t288 = t301 * t309 + t306 * t319;
	t302 = sin(pkin(6));
	t325 = t302 * t339;
	t278 = t288 * t308 - t305 * t325;
	t326 = t301 * t340;
	t290 = -t306 * t326 + t339 * t309;
	t304 = sin(qJ(4));
	t315 = t304 * pkin(4) + t298 * r_i_i_C(1) + t299 * r_i_i_C(2);
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
	t1 = [0, (-t285 * t298 + t290 * t332) * r_i_i_C(1) + (-t285 * t299 - t290 * t333) * r_i_i_C(2) - t285 * pkin(8) + (-t285 * t304 + t290 * t331) * pkin(4) + t311 * t286 + t310 * t289, t280 * qJD(5) - t318 * t275 + t341 * t276 - t316 * t312, (-t276 * t298 + t286 * t299) * r_i_i_C(1) + (-t276 * t299 - t286 * t298) * r_i_i_C(2) + ((-t280 * t299 - t289 * t298) * r_i_i_C(1) + (t280 * t298 - t289 * t299) * r_i_i_C(2)) * qJD(4) + (-t276 * t304 + t286 * t307 + (-t280 * t307 - t289 * t304) * qJD(4)) * pkin(4), t275, 0; 0, (-t283 * t298 + t288 * t332) * r_i_i_C(1) + (-t283 * t299 - t288 * t333) * r_i_i_C(2) - t283 * pkin(8) + (-t283 * t304 + t288 * t331) * pkin(4) + t311 * t284 + t310 * t287, t278 * qJD(5) - t318 * t273 + t341 * t274 - t314 * t312, (-t274 * t298 + t284 * t299) * r_i_i_C(1) + (-t274 * t299 - t284 * t298) * r_i_i_C(2) + ((-t278 * t299 - t287 * t298) * r_i_i_C(1) + (t278 * t298 - t287 * t299) * r_i_i_C(2)) * qJD(4) + (-t274 * t304 + t284 * t307 + (-t278 * t307 - t287 * t304) * qJD(4)) * pkin(4), t273, 0; 0, ((t311 * qJD(2) + t347 * qJD(4)) * t306 + (qJD(2) * pkin(8) + t315 * (qJD(2) - t330) + t348) * t309) * t302, t292 * qJD(5) - t318 * t281 + t341 * t282 - t313 * t312, (-t282 * t298 + t299 * t328) * r_i_i_C(1) + (-t282 * t299 - t298 * t328) * r_i_i_C(2) + ((-t292 * t299 + t298 * t335) * r_i_i_C(1) + (t292 * t298 + t299 * t335) * r_i_i_C(2)) * qJD(4) + (t307 * t328 - t282 * t304 + (-t292 * t307 + t304 * t335) * qJD(4)) * pkin(4), t281, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:11
	% EndTime: 2019-10-09 22:41:12
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (691->114), mult. (1620->184), div. (0->0), fcn. (1673->12), ass. (0->70)
	t379 = cos(qJ(4));
	t369 = t379 * pkin(4) + pkin(3);
	t377 = sin(qJ(3));
	t380 = cos(qJ(3));
	t372 = qJ(4) + pkin(11);
	t370 = sin(t372);
	t411 = t370 * qJD(6);
	t424 = r_i_i_C(2) + qJ(5) + pkin(9);
	t376 = sin(qJ(4));
	t425 = t376 * pkin(4);
	t382 = -(t424 * qJD(3) - qJD(4) * t425 + t411) * t380 + (qJD(3) * t369 - qJD(5)) * t377;
	t373 = sin(pkin(10));
	t378 = sin(qJ(2));
	t381 = cos(qJ(2));
	t421 = cos(pkin(10));
	t422 = cos(pkin(6));
	t396 = t422 * t421;
	t358 = t373 * t381 + t378 * t396;
	t374 = sin(pkin(6));
	t404 = t374 * t421;
	t346 = t358 * t380 - t377 * t404;
	t405 = t373 * t422;
	t360 = -t378 * t405 + t421 * t381;
	t417 = t374 * t380;
	t362 = t422 * t377 + t378 * t417;
	t371 = cos(t372);
	t416 = t374 * t381;
	t429 = -t362 * t371 + t370 * t416;
	t414 = qJD(3) * t377;
	t428 = (qJD(2) * t380 - qJD(4)) * t378 + t381 * t414;
	t426 = pkin(5) + r_i_i_C(1);
	t423 = r_i_i_C(3) + qJ(6);
	t418 = t374 * t377;
	t415 = qJD(2) * t378;
	t413 = qJD(4) * t379;
	t412 = qJD(4) * t380;
	t410 = t371 * qJD(6);
	t408 = t374 * t415;
	t407 = qJD(2) * t416;
	t392 = t381 * t396;
	t353 = -qJD(2) * t392 + t373 * t415;
	t357 = t373 * t378 - t392;
	t398 = t357 * t412 - t353;
	t359 = t421 * t378 + t381 * t405;
	t355 = t359 * qJD(2);
	t397 = t359 * t412 - t355;
	t395 = t346 * t371 + t357 * t370;
	t348 = t360 * t380 + t373 * t418;
	t394 = t348 * t371 + t359 * t370;
	t393 = (qJD(2) - t412) * t381;
	t390 = -t360 * t377 + t373 * t417;
	t389 = -t380 * t369 - t424 * t377 - pkin(2);
	t388 = -t358 * t377 - t380 * t404;
	t387 = -t378 * t418 + t422 * t380;
	t386 = -t423 * t370 - t426 * t371 - t369;
	t354 = t358 * qJD(2);
	t385 = qJD(4) * t358 - t354 * t380 + t357 * t414;
	t356 = t360 * qJD(2);
	t384 = qJD(4) * t360 - t356 * t380 + t359 * t414;
	t383 = t411 + (-t426 * t370 + t423 * t371 - t425) * qJD(4);
	t350 = t387 * qJD(3) + t380 * t407;
	t349 = t362 * qJD(3) + t377 * t407;
	t344 = t390 * qJD(3) - t355 * t380;
	t343 = t348 * qJD(3) - t355 * t377;
	t342 = t388 * qJD(3) - t353 * t380;
	t341 = t346 * qJD(3) - t353 * t377;
	t337 = -t429 * qJD(4) + t350 * t370 - t371 * t408;
	t331 = t394 * qJD(4) + t344 * t370 - t356 * t371;
	t329 = t395 * qJD(4) + t342 * t370 - t354 * t371;
	t1 = [0, -t360 * t410 - t355 * pkin(8) + t426 * (t397 * t370 + t384 * t371) + t423 * (t384 * t370 - t397 * t371) + (-t355 * t376 + t360 * t413) * pkin(4) + t389 * t356 + t382 * t359, t348 * qJD(5) + t386 * t343 + t424 * t344 + t383 * t390, t394 * qJD(6) + t423 * (t344 * t371 + t356 * t370 + (-t348 * t370 + t359 * t371) * qJD(4)) - t426 * t331 + (-t344 * t376 + t356 * t379 + (-t348 * t379 - t359 * t376) * qJD(4)) * pkin(4), t343, t331; 0, -t358 * t410 - t353 * pkin(8) + t426 * (t398 * t370 + t385 * t371) + t423 * (t385 * t370 - t398 * t371) + (-t353 * t376 + t358 * t413) * pkin(4) + t389 * t354 + t382 * t357, t346 * qJD(5) + t386 * t341 + t424 * t342 + t383 * t388, t395 * qJD(6) + t423 * (t342 * t371 + t354 * t370 + (-t346 * t370 + t357 * t371) * qJD(4)) - t426 * t329 + (-t342 * t376 + t354 * t379 + (-t346 * t379 - t357 * t376) * qJD(4)) * pkin(4), t341, t329; 0, (t426 * (t370 * t393 - t428 * t371) - t423 * (t428 * t370 + t371 * t393) + (pkin(4) * t413 + t389 * qJD(2) - t410) * t378 + ((pkin(8) + t425) * qJD(2) - t382) * t381) * t374, t362 * qJD(5) + t386 * t349 + t424 * t350 + t383 * t387, -t429 * qJD(6) + t423 * (t370 * t408 + t350 * t371 + (-t362 * t370 - t371 * t416) * qJD(4)) - t426 * t337 + (t379 * t408 - t350 * t376 + (-t362 * t379 + t376 * t416) * qJD(4)) * pkin(4), t349, t337;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end