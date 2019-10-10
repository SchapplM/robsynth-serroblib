% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:51
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
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
	% StartTime: 2019-10-10 11:51:32
	% EndTime: 2019-10-10 11:51:32
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 11:51:33
	% EndTime: 2019-10-10 11:51:33
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-10 11:51:33
	% EndTime: 2019-10-10 11:51:34
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
	t318 = r_i_i_C(2) - pkin(3);
	t320 = t317 * t286 - t318 * t289 + pkin(2);
	t319 = -r_i_i_C(1) - pkin(9);
	t296 = t305 - t307;
	t316 = t296 * t286;
	t310 = t288 * t289;
	t303 = qJD(2) * t284 * t290;
	t300 = -t289 * t296 + t306;
	t299 = t285 * t286 + t287 * t314;
	t298 = t285 * t307 - t311;
	t297 = t285 * t309 + t308;
	t292 = qJD(4) * t286 + (t318 * t286 + t317 * t289) * qJD(3);
	t274 = t299 * qJD(3) + t286 * t303;
	t272 = t297 * qJD(1) + t278 * qJD(2);
	t271 = t278 * qJD(1) + t297 * qJD(2);
	t270 = -t298 * qJD(1) + t296 * qJD(2);
	t265 = -t271 * t289 + qJD(3) * t316 + (qJD(1) * t312 + qJD(3) * t310) * t284;
	t264 = t300 * qJD(3) - t271 * t286 - t291 * t304;
	t1 = [-t302 * qJD(4) - t273 * pkin(2) + t319 * t272 - t318 * t322 + t317 * t321 + (-t291 * pkin(1) - pkin(8) * t315) * qJD(1), t320 * t270 + t319 * t271 - t292 * t297, t300 * qJD(4) + t318 * t264 + t317 * t265, t264, 0, 0; -(t284 * t310 + t316) * qJD(4) - t271 * pkin(2) + t319 * t270 - t318 * t265 + t317 * t264 + (-t288 * pkin(1) + pkin(8) * t313) * qJD(1), -t272 * t320 - t319 * t273 + t292 * t298, -t301 * qJD(4) - t317 * t322 - t318 * t321, -t321, 0, 0; 0, (t292 * t290 + (-t287 * t320 - t319 * t290) * qJD(2)) * t284, t299 * qJD(4) + t317 * (t289 * t303 + (-t284 * t286 * t287 + t285 * t289) * qJD(3)) + t318 * t274, t274, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:34
	% EndTime: 2019-10-10 11:51:35
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (545->92), mult. (1620->157), div. (0->0), fcn. (1626->10), ass. (0->64)
	t376 = cos(qJ(2));
	t377 = cos(qJ(1));
	t415 = cos(pkin(6));
	t393 = t377 * t415;
	t391 = t376 * t393;
	t372 = sin(qJ(2));
	t373 = sin(qJ(1));
	t405 = t373 * t372;
	t356 = -t391 + t405;
	t370 = sin(qJ(5));
	t374 = cos(qJ(5));
	t357 = t372 * t393 + t373 * t376;
	t371 = sin(qJ(3));
	t375 = cos(qJ(3));
	t369 = sin(pkin(6));
	t407 = t369 * t377;
	t420 = t357 * t371 + t375 * t407;
	t423 = t356 * t374 + t370 * t420;
	t422 = t356 * t370 - t374 * t420;
	t390 = t374 * r_i_i_C(1) - t370 * r_i_i_C(2);
	t380 = t390 * qJD(5) + qJD(4);
	t389 = -t370 * r_i_i_C(1) - t374 * r_i_i_C(2);
	t387 = qJ(4) - t389;
	t400 = r_i_i_C(3) + pkin(10) + pkin(3);
	t378 = (t400 * t371 - t387 * t375) * qJD(3) - t380 * t371;
	t388 = qJD(2) * t415 + qJD(1);
	t394 = t373 * t415;
	t392 = t372 * t394;
	t403 = qJD(2) * t372;
	t404 = t377 * t376;
	t345 = -qJD(1) * t392 - t373 * t403 + t388 * t404;
	t410 = t369 * t373;
	t399 = t371 * t410;
	t419 = -qJD(1) * t399 + qJD(3) * t420 - t345 * t375;
	t417 = t387 * t371 + t400 * t375 + pkin(2);
	t416 = -pkin(4) - pkin(9);
	t359 = -t392 + t404;
	t411 = t359 * t371;
	t409 = t369 * t375;
	t408 = t369 * t376;
	t406 = t371 * t377;
	t402 = qJD(2) * t376;
	t401 = qJD(3) * t375;
	t398 = t369 * t406;
	t397 = qJD(1) * t409;
	t396 = t369 * t403;
	t395 = t369 * t402;
	t386 = t390 - t416;
	t384 = t359 * t375 + t399;
	t383 = t389 * qJD(5);
	t358 = t377 * t372 + t376 * t394;
	t354 = t369 * t372 * t371 - t415 * t375;
	t381 = t415 * t371 + t372 * t409;
	t338 = -qJD(3) * t398 + t345 * t371 + t357 * t401 - t373 * t397;
	t351 = -t373 * t409 + t411;
	t346 = t381 * qJD(3) + t371 * t395;
	t344 = t358 * qJD(1) + t357 * qJD(2);
	t343 = t357 * qJD(1) + t358 * qJD(2);
	t342 = -qJD(1) * t391 - t377 * t402 + t388 * t405;
	t337 = -t343 * t375 - qJD(3) * t411 + (qJD(1) * t406 + t373 * t401) * t369;
	t336 = t384 * qJD(3) - t343 * t371 - t377 * t397;
	t335 = t336 * t370 - t342 * t374 + (t351 * t374 - t358 * t370) * qJD(5);
	t334 = t336 * t374 + t342 * t370 + (-t351 * t370 - t358 * t374) * qJD(5);
	t1 = [-t345 * pkin(2) - t420 * qJD(4) - t387 * t338 - t386 * t344 + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(5) + (-t377 * pkin(1) - pkin(8) * t410) * qJD(1) + t400 * t419, t417 * t342 - t386 * t343 + t378 * t358 + t359 * t383, -t400 * t336 + t387 * t337 + t380 * t384, t336, t334 * r_i_i_C(1) - t335 * r_i_i_C(2), 0; -t343 * pkin(2) + t335 * r_i_i_C(1) + t334 * r_i_i_C(2) + t336 * qJ(4) + t351 * qJD(4) + t416 * t342 + (-pkin(1) * t373 + pkin(8) * t407) * qJD(1) + t400 * t337, -t344 * t417 + t386 * t345 + t378 * t356 + t357 * t383, t380 * (t357 * t375 - t398) - t387 * t419 - t400 * t338, t338, (t338 * t374 - t344 * t370) * r_i_i_C(1) + (-t338 * t370 - t344 * t374) * r_i_i_C(2) + (-t423 * r_i_i_C(1) + t422 * r_i_i_C(2)) * qJD(5), 0; 0, ((-qJD(2) * t417 + t383) * t372 + (t386 * qJD(2) - t378) * t376) * t369, t380 * t381 + t387 * (-t354 * qJD(3) + t375 * t395) - t400 * t346, t346, (t346 * t374 - t370 * t396) * r_i_i_C(1) + (-t346 * t370 - t374 * t396) * r_i_i_C(2) + ((-t354 * t370 + t374 * t408) * r_i_i_C(1) + (-t354 * t374 - t370 * t408) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:34
	% EndTime: 2019-10-10 11:51:35
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (714->113), mult. (2093->176), div. (0->0), fcn. (2106->10), ass. (0->72)
	t380 = sin(qJ(3));
	t384 = cos(qJ(3));
	t379 = sin(qJ(5));
	t383 = cos(qJ(5));
	t434 = pkin(5) + r_i_i_C(1);
	t400 = -t383 * r_i_i_C(2) - t434 * t379;
	t397 = -qJ(4) + t400;
	t417 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(10);
	t442 = -(t417 * t380 + t397 * t384) * qJD(3) + t384 * qJD(6);
	t387 = t397 * t380 - t417 * t384 - pkin(2);
	t385 = cos(qJ(2));
	t429 = cos(pkin(6));
	t433 = cos(qJ(1));
	t405 = t429 * t433;
	t401 = t385 * t405;
	t381 = sin(qJ(2));
	t382 = sin(qJ(1));
	t422 = t381 * t382;
	t363 = -t401 + t422;
	t364 = t381 * t405 + t382 * t385;
	t377 = sin(pkin(6));
	t416 = t377 * t433;
	t372 = t384 * t416;
	t409 = t364 * t380 + t372;
	t440 = t363 * t383 + t379 * t409;
	t439 = t363 * t379 - t383 * t409;
	t432 = t379 * r_i_i_C(2);
	t389 = qJD(4) + (t434 * t383 - t432) * qJD(5);
	t410 = t382 * t429;
	t406 = t381 * t410;
	t411 = t433 * qJD(1);
	t421 = qJD(2) * t381;
	t352 = -qJD(1) * t406 - t382 * t421 + (qJD(2) * t405 + t411) * t385;
	t425 = t377 * t382;
	t414 = qJD(1) * t425;
	t435 = qJD(3) * t409 - t352 * t384 - t380 * t414;
	t430 = -pkin(5) * t383 - pkin(4) - pkin(9);
	t415 = t433 * t385;
	t366 = t415 - t406;
	t426 = t366 * t380;
	t424 = t377 * t384;
	t423 = t377 * t385;
	t420 = qJD(3) * t384;
	t419 = qJD(5) * t380;
	t413 = t377 * t421;
	t412 = qJD(2) * t423;
	t408 = -t430 - t432;
	t407 = t380 * t416;
	t345 = -qJD(3) * t407 + t352 * t380 + t364 * t420 - t384 * t414;
	t365 = t433 * t381 + t385 * t410;
	t351 = t365 * qJD(1) + t364 * qJD(2);
	t404 = t345 * t383 - t351 * t379;
	t361 = t377 * t381 * t380 - t429 * t384;
	t399 = -t361 * t379 + t383 * t423;
	t359 = t366 * t384 + t380 * t425;
	t398 = t383 * r_i_i_C(1) + t408;
	t395 = t364 * t384 - t407;
	t390 = t429 * t380 + t381 * t424;
	t353 = t390 * qJD(3) + t380 * t412;
	t392 = t353 * t383 - t379 * t413;
	t391 = t400 * qJD(5);
	t350 = t364 * qJD(1) + t365 * qJD(2);
	t343 = -qJD(1) * t372 + t359 * qJD(3) - t350 * t380;
	t358 = -t382 * t424 + t426;
	t388 = t343 * t379 + (t358 * t383 - t365 * t379) * qJD(5);
	t349 = -qJD(1) * t401 - qJD(2) * t415 + (qJD(2) * t429 + qJD(1)) * t422;
	t341 = t343 * t383 + t349 * t379 + (-t358 * t379 - t365 * t383) * qJD(5);
	t386 = -t389 * t380 - t442;
	t354 = -t361 * qJD(3) + t384 * t412;
	t344 = -t350 * t384 - qJD(3) * t426 + (t380 * t411 + t382 * t420) * t377;
	t342 = -t349 * t383 + t388;
	t1 = [-t395 * qJD(6) - t409 * qJD(4) - t352 * pkin(2) - t398 * t351 + t397 * t345 + (-t433 * pkin(1) - pkin(8) * t425) * qJD(1) + t417 * t435 + (t440 * r_i_i_C(2) + t434 * t439) * qJD(5), -t387 * t349 - t398 * t350 + t386 * t365 + t366 * t391, -qJD(6) * t358 - t417 * t343 - t344 * t397 + t389 * t359, t343, -r_i_i_C(2) * t342 + t434 * t341, t344; -t350 * pkin(2) + t342 * r_i_i_C(1) + t341 * r_i_i_C(2) + t343 * qJ(4) + t358 * qJD(4) + t359 * qJD(6) + t430 * t349 + (-pkin(1) * t382 + pkin(8) * t416) * qJD(1) + t417 * t344 + t388 * pkin(5), t387 * t351 + t398 * t352 + t386 * t363 + t364 * t391, -qJD(6) * t409 - t417 * t345 + t389 * t395 + t397 * t435, t345, t404 * r_i_i_C(1) + (-t345 * t379 - t351 * t383) * r_i_i_C(2) + (-r_i_i_C(1) * t440 + t439 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t440 + t404) * pkin(5), -t435; 0, ((t387 * qJD(2) + t391) * t381 + ((-qJD(5) * t432 + qJD(4)) * t380 + t408 * qJD(2) + ((qJD(2) + t419) * r_i_i_C(1) + pkin(5) * t419) * t383 + t442) * t385) * t377, -qJD(6) * t361 - t417 * t353 - t354 * t397 + t389 * t390, t353, t392 * r_i_i_C(1) + (-t353 * t379 - t383 * t413) * r_i_i_C(2) + (t399 * r_i_i_C(1) + (-t361 * t383 - t379 * t423) * r_i_i_C(2)) * qJD(5) + (t399 * qJD(5) + t392) * pkin(5), t354;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end