% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:21
	% EndTime: 2019-10-10 11:33:21
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(6));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(9);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(8) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(8) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:21
	% EndTime: 2019-10-10 11:33:22
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (278->59), mult. (828->98), div. (0->0), fcn. (802->8), ass. (0->45)
	t285 = cos(pkin(6));
	t287 = sin(qJ(2));
	t288 = sin(qJ(1));
	t311 = t288 * t287;
	t305 = t285 * t311;
	t290 = cos(qJ(2));
	t291 = cos(qJ(1));
	t307 = t291 * t290;
	t273 = -qJD(1) * t305 - qJD(2) * t311 + (qJD(2) * t285 + qJD(1)) * t307;
	t289 = cos(qJ(3));
	t308 = t291 * t287;
	t309 = t288 * t290;
	t278 = t285 * t308 + t309;
	t286 = sin(qJ(3));
	t284 = sin(pkin(6));
	t313 = t284 * t291;
	t302 = t278 * t286 + t289 * t313;
	t315 = t284 * t288;
	t306 = t286 * t315;
	t322 = -qJD(1) * t306 + t302 * qJD(3) - t273 * t289;
	t312 = t286 * t291;
	t301 = -t278 * t289 + t284 * t312;
	t314 = t284 * t289;
	t304 = qJD(1) * t314;
	t321 = t301 * qJD(3) - t273 * t286 + t288 * t304;
	t317 = r_i_i_C(3) + qJ(4);
	t319 = pkin(3) - r_i_i_C(2);
	t320 = t317 * t286 + t319 * t289 + pkin(2);
	t318 = pkin(9) + r_i_i_C(1);
	t296 = t305 - t307;
	t316 = t296 * t286;
	t310 = t288 * t289;
	t303 = qJD(2) * t284 * t290;
	t300 = -t289 * t296 + t306;
	t299 = t285 * t286 + t287 * t314;
	t298 = t285 * t307 - t311;
	t297 = t285 * t309 + t308;
	t292 = qJD(4) * t286 + (-t319 * t286 + t317 * t289) * qJD(3);
	t274 = t299 * qJD(3) + t286 * t303;
	t272 = t297 * qJD(1) + t278 * qJD(2);
	t271 = t278 * qJD(1) + t297 * qJD(2);
	t270 = -t298 * qJD(1) + t296 * qJD(2);
	t265 = -t271 * t289 + qJD(3) * t316 + (qJD(1) * t312 + qJD(3) * t310) * t284;
	t264 = t300 * qJD(3) - t271 * t286 - t291 * t304;
	t1 = [-t302 * qJD(4) - t273 * pkin(2) - t318 * t272 + t319 * t322 + t317 * t321 + (-t291 * pkin(1) - pkin(8) * t315) * qJD(1), t320 * t270 - t318 * t271 - t292 * t297, t300 * qJD(4) - t319 * t264 + t317 * t265, t264, 0, 0; -(t284 * t310 + t316) * qJD(4) - t271 * pkin(2) - t318 * t270 + t319 * t265 + t317 * t264 + (-t288 * pkin(1) + pkin(8) * t313) * qJD(1), -t272 * t320 + t318 * t273 + t292 * t298, -t301 * qJD(4) - t317 * t322 + t319 * t321, -t321, 0, 0; 0, (t292 * t290 + (-t287 * t320 + t318 * t290) * qJD(2)) * t284, t299 * qJD(4) + t317 * (t289 * t303 + (-t284 * t286 * t287 + t285 * t289) * qJD(3)) - t319 * t274, t274, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:22
	% EndTime: 2019-10-10 11:33:22
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (445->70), mult. (1325->111), div. (0->0), fcn. (1308->10), ass. (0->54)
	t323 = cos(pkin(6));
	t325 = sin(qJ(2));
	t329 = cos(qJ(1));
	t354 = t329 * t325;
	t326 = sin(qJ(1));
	t328 = cos(qJ(2));
	t355 = t326 * t328;
	t309 = t323 * t354 + t355;
	t324 = sin(qJ(3));
	t327 = cos(qJ(3));
	t321 = sin(pkin(6));
	t358 = t321 * t329;
	t365 = t309 * t324 + t327 * t358;
	t341 = qJD(2) * t323 + qJD(1);
	t356 = t326 * t325;
	t347 = t323 * t356;
	t352 = qJD(2) * t325;
	t353 = t329 * t328;
	t302 = -qJD(1) * t347 - t326 * t352 + t341 * t353;
	t360 = t321 * t326;
	t348 = t324 * t360;
	t364 = -qJD(1) * t348 + qJD(3) * t365 - t302 * t327;
	t320 = sin(pkin(11));
	t322 = cos(pkin(11));
	t340 = t320 * r_i_i_C(1) + t322 * r_i_i_C(2) + qJ(4);
	t349 = r_i_i_C(3) + qJ(5) + pkin(3);
	t363 = t340 * t324 + t349 * t327 + pkin(2);
	t311 = -t347 + t353;
	t361 = t311 * t324;
	t359 = t321 * t327;
	t357 = t324 * t329;
	t351 = qJD(2) * t328;
	t350 = qJD(3) * t327;
	t346 = t321 * t357;
	t345 = t323 * t353;
	t344 = qJD(1) * t359;
	t343 = t321 * t351;
	t339 = t322 * r_i_i_C(1) - t320 * r_i_i_C(2) + pkin(4) + pkin(9);
	t337 = -t309 * t327 + t346;
	t336 = t311 * t327 + t348;
	t335 = t323 * t324 + t325 * t359;
	t334 = t321 * t325 * t324 - t323 * t327;
	t333 = t323 * t355 + t354;
	t295 = -qJD(3) * t346 + t302 * t324 + t309 * t350 - t326 * t344;
	t330 = t324 * qJD(4) + t327 * qJD(5) + (-t349 * t324 + t340 * t327) * qJD(3);
	t305 = -t326 * t359 + t361;
	t304 = -qJD(3) * t334 + t327 * t343;
	t303 = qJD(3) * t335 + t324 * t343;
	t301 = qJD(1) * t333 + qJD(2) * t309;
	t300 = qJD(1) * t309 + qJD(2) * t333;
	t299 = -qJD(1) * t345 - t329 * t351 + t341 * t356;
	t294 = -t300 * t327 - qJD(3) * t361 + (qJD(1) * t357 + t326 * t350) * t321;
	t293 = qJD(3) * t336 - t300 * t324 - t329 * t344;
	t1 = [t337 * qJD(5) - t365 * qJD(4) - t302 * pkin(2) - t340 * t295 - t339 * t301 + (-t329 * pkin(1) - pkin(8) * t360) * qJD(1) + t349 * t364, t363 * t299 - t339 * t300 - t330 * t333, qJD(4) * t336 - t305 * qJD(5) - t349 * t293 + t340 * t294, t293, t294, 0; -t300 * pkin(2) + t305 * qJD(4) + t336 * qJD(5) + t340 * t293 - t339 * t299 + (-t326 * pkin(1) + pkin(8) * t358) * qJD(1) + t349 * t294, t339 * t302 - t363 * t301 + t330 * (t345 - t356), -t337 * qJD(4) - qJD(5) * t365 - t349 * t295 - t340 * t364, t295, -t364, 0; 0, (-t363 * t352 + (qJD(2) * t339 + t330) * t328) * t321, t335 * qJD(4) - t334 * qJD(5) - t349 * t303 + t340 * t304, t303, t304, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:22
	% EndTime: 2019-10-10 11:33:23
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (727->102), mult. (1878->165), div. (0->0), fcn. (1892->12), ass. (0->68)
	t386 = cos(qJ(2));
	t387 = cos(qJ(1));
	t427 = cos(pkin(6));
	t403 = t387 * t427;
	t401 = t386 * t403;
	t383 = sin(qJ(2));
	t384 = sin(qJ(1));
	t417 = t383 * t384;
	t362 = -t401 + t417;
	t378 = pkin(11) + qJ(6);
	t376 = sin(t378);
	t377 = cos(t378);
	t363 = t383 * t403 + t384 * t386;
	t382 = sin(qJ(3));
	t385 = cos(qJ(3));
	t380 = sin(pkin(6));
	t419 = t380 * t387;
	t432 = t363 * t382 + t385 * t419;
	t435 = t362 * t377 + t376 * t432;
	t434 = t362 * t376 - t377 * t432;
	t400 = t377 * r_i_i_C(1) - t376 * r_i_i_C(2);
	t390 = t400 * qJD(6) + qJD(4);
	t399 = -t376 * r_i_i_C(1) - t377 * r_i_i_C(2);
	t405 = pkin(5) * sin(pkin(11)) + qJ(4);
	t391 = -t399 + t405;
	t411 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(3);
	t388 = (t411 * t382 - t391 * t385) * qJD(3) - t385 * qJD(5) - t390 * t382;
	t398 = qJD(2) * t427 + qJD(1);
	t404 = t384 * t427;
	t402 = t383 * t404;
	t415 = qJD(2) * t383;
	t416 = t387 * t386;
	t351 = -qJD(1) * t402 - t384 * t415 + t398 * t416;
	t422 = t380 * t384;
	t410 = t382 * t422;
	t431 = -qJD(1) * t410 + qJD(3) * t432 - t351 * t385;
	t429 = t391 * t382 + t411 * t385 + pkin(2);
	t428 = -pkin(9) - cos(pkin(11)) * pkin(5) - pkin(4);
	t365 = -t402 + t416;
	t423 = t365 * t382;
	t421 = t380 * t385;
	t420 = t380 * t386;
	t418 = t382 * t387;
	t414 = qJD(2) * t386;
	t413 = qJD(3) * t385;
	t409 = t380 * t418;
	t408 = qJD(1) * t421;
	t407 = t380 * t415;
	t406 = t380 * t414;
	t396 = t363 * t385 - t409;
	t358 = t365 * t385 + t410;
	t395 = t399 * qJD(6);
	t394 = t400 - t428;
	t364 = t387 * t383 + t386 * t404;
	t360 = t380 * t383 * t382 - t427 * t385;
	t392 = t427 * t382 + t383 * t421;
	t344 = -qJD(3) * t409 + t351 * t382 + t363 * t413 - t384 * t408;
	t357 = -t384 * t421 + t423;
	t353 = -t360 * qJD(3) + t385 * t406;
	t352 = t392 * qJD(3) + t382 * t406;
	t350 = t364 * qJD(1) + t363 * qJD(2);
	t349 = t363 * qJD(1) + t364 * qJD(2);
	t348 = -qJD(1) * t401 - t387 * t414 + t398 * t417;
	t343 = -t349 * t385 - qJD(3) * t423 + (qJD(1) * t418 + t384 * t413) * t380;
	t342 = t358 * qJD(3) - t349 * t382 - t387 * t408;
	t341 = t342 * t376 - t348 * t377 + (t357 * t377 - t364 * t376) * qJD(6);
	t340 = t342 * t377 + t348 * t376 + (-t357 * t376 - t364 * t377) * qJD(6);
	t1 = [-t396 * qJD(5) - t432 * qJD(4) - t351 * pkin(2) - t394 * t350 - t391 * t344 + (t434 * r_i_i_C(1) + t435 * r_i_i_C(2)) * qJD(6) + (-pkin(1) * t387 - pkin(8) * t422) * qJD(1) + t411 * t431, t429 * t348 - t394 * t349 + t388 * t364 + t365 * t395, -qJD(5) * t357 - t411 * t342 + t391 * t343 + t390 * t358, t342, t343, r_i_i_C(1) * t340 - r_i_i_C(2) * t341; -t349 * pkin(2) + t341 * r_i_i_C(1) + t340 * r_i_i_C(2) + t357 * qJD(4) + t358 * qJD(5) + t428 * t348 + t405 * t342 + (-pkin(1) * t384 + pkin(8) * t419) * qJD(1) + t411 * t343, -t350 * t429 + t394 * t351 + t388 * t362 + t363 * t395, -qJD(5) * t432 - t411 * t344 + t390 * t396 - t391 * t431, t344, -t431, (t344 * t377 - t350 * t376) * r_i_i_C(1) + (-t344 * t376 - t350 * t377) * r_i_i_C(2) + (-t435 * r_i_i_C(1) + t434 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t429 + t395) * t383 + (t394 * qJD(2) - t388) * t386) * t380, -qJD(5) * t360 - t411 * t352 + t391 * t353 + t390 * t392, t352, t353, (t352 * t377 - t376 * t407) * r_i_i_C(1) + (-t352 * t376 - t377 * t407) * r_i_i_C(2) + ((-t360 * t376 + t377 * t420) * r_i_i_C(1) + (-t360 * t377 - t376 * t420) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end