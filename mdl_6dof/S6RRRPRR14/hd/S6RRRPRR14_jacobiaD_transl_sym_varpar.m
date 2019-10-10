% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:16:27
	% EndTime: 2019-10-10 12:16:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:16:27
	% EndTime: 2019-10-10 12:16:27
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
	% StartTime: 2019-10-10 12:16:28
	% EndTime: 2019-10-10 12:16:28
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
	% StartTime: 2019-10-10 12:16:29
	% EndTime: 2019-10-10 12:16:29
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
	% StartTime: 2019-10-10 12:16:29
	% EndTime: 2019-10-10 12:16:30
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
	t318 = -r_i_i_C(2) + pkin(3);
	t320 = t317 * t286 + t318 * t289 + pkin(2);
	t319 = pkin(9) + r_i_i_C(1);
	t296 = t305 - t307;
	t316 = t296 * t286;
	t310 = t288 * t289;
	t303 = qJD(2) * t284 * t290;
	t300 = -t289 * t296 + t306;
	t299 = t285 * t286 + t287 * t314;
	t298 = t285 * t307 - t311;
	t297 = t285 * t309 + t308;
	t292 = qJD(4) * t286 + (-t318 * t286 + t317 * t289) * qJD(3);
	t274 = t299 * qJD(3) + t286 * t303;
	t272 = t297 * qJD(1) + t278 * qJD(2);
	t271 = t278 * qJD(1) + t297 * qJD(2);
	t270 = -t298 * qJD(1) + t296 * qJD(2);
	t265 = -t271 * t289 + qJD(3) * t316 + (qJD(1) * t312 + qJD(3) * t310) * t284;
	t264 = t300 * qJD(3) - t271 * t286 - t291 * t304;
	t1 = [-t302 * qJD(4) - t273 * pkin(2) - t319 * t272 + t318 * t322 + t317 * t321 + (-t291 * pkin(1) - pkin(8) * t315) * qJD(1), t320 * t270 - t319 * t271 - t292 * t297, t300 * qJD(4) - t318 * t264 + t317 * t265, t264, 0, 0; -(t284 * t310 + t316) * qJD(4) - t271 * pkin(2) - t319 * t270 + t318 * t265 + t317 * t264 + (-t288 * pkin(1) + pkin(8) * t313) * qJD(1), -t272 * t320 + t319 * t273 + t292 * t298, -t301 * qJD(4) - t317 * t322 + t318 * t321, -t321, 0, 0; 0, (t292 * t290 + (-t287 * t320 + t319 * t290) * qJD(2)) * t284, t299 * qJD(4) + t317 * (t289 * t303 + (-t284 * t286 * t287 + t285 * t289) * qJD(3)) - t318 * t274, t274, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:16:30
	% EndTime: 2019-10-10 12:16:31
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
	% StartTime: 2019-10-10 12:16:30
	% EndTime: 2019-10-10 12:16:31
	% DurationCPUTime: 0.87s
	% Computational Cost: add. (914->113), mult. (2230->180), div. (0->0), fcn. (2242->12), ass. (0->75)
	t411 = sin(qJ(3));
	t415 = cos(qJ(3));
	t407 = qJD(5) + qJD(6);
	t414 = cos(qJ(5));
	t408 = qJ(5) + qJ(6);
	t405 = sin(t408);
	t406 = cos(t408);
	t433 = t406 * r_i_i_C(1) - t405 * r_i_i_C(2);
	t462 = qJD(5) * pkin(5);
	t420 = t433 * t407 + t414 * t462 + qJD(4);
	t410 = sin(qJ(5));
	t432 = -t405 * r_i_i_C(1) - t406 * r_i_i_C(2);
	t429 = qJ(4) - t432;
	t422 = pkin(5) * t410 + t429;
	t450 = pkin(3) + r_i_i_C(3) + pkin(11) + pkin(10);
	t418 = (t450 * t411 - t422 * t415) * qJD(3) - t420 * t411;
	t413 = sin(qJ(1));
	t416 = cos(qJ(2));
	t461 = cos(pkin(6));
	t464 = cos(qJ(1));
	t435 = t461 * t464;
	t412 = sin(qJ(2));
	t443 = t413 * t461;
	t436 = t412 * t443;
	t444 = t464 * qJD(1);
	t452 = qJD(2) * t412;
	t380 = -qJD(1) * t436 - t413 * t452 + (qJD(2) * t435 + t444) * t416;
	t392 = t412 * t435 + t413 * t416;
	t409 = sin(pkin(6));
	t449 = t409 * t464;
	t400 = t415 * t449;
	t438 = t392 * t411 + t400;
	t459 = t409 * t413;
	t447 = qJD(1) * t459;
	t467 = qJD(3) * t438 - t380 * t415 - t411 * t447;
	t465 = t422 * t411 + t450 * t415 + pkin(2);
	t463 = -t414 * pkin(5) - pkin(4) - pkin(9);
	t448 = t464 * t416;
	t394 = t448 - t436;
	t460 = t394 * t411;
	t458 = t409 * t415;
	t457 = t409 * t416;
	t456 = t413 * t412;
	t430 = t416 * t435;
	t377 = -qJD(1) * t430 - qJD(2) * t448 + (qJD(2) * t461 + qJD(1)) * t456;
	t386 = -t413 * t458 + t460;
	t440 = t386 * t407 - t377;
	t393 = t464 * t412 + t416 * t443;
	t378 = t392 * qJD(1) + t393 * qJD(2);
	t428 = t394 * t415 + t411 * t459;
	t371 = -qJD(1) * t400 + t428 * qJD(3) - t378 * t411;
	t442 = -t393 * t407 + t371;
	t367 = -t440 * t405 + t442 * t406;
	t368 = t442 * t405 + t440 * t406;
	t455 = t367 * r_i_i_C(1) - t368 * r_i_i_C(2);
	t379 = t393 * qJD(1) + t392 * qJD(2);
	t439 = -t407 * t438 - t379;
	t437 = t411 * t449;
	t451 = qJD(3) * t415;
	t373 = -qJD(3) * t437 + t380 * t411 + t392 * t451 - t415 * t447;
	t391 = -t430 + t456;
	t441 = t391 * t407 - t373;
	t454 = (t439 * t405 - t441 * t406) * r_i_i_C(1) + (t441 * t405 + t439 * t406) * r_i_i_C(2);
	t389 = t409 * t412 * t411 - t461 * t415;
	t446 = t409 * t452;
	t427 = -t389 * t407 - t446;
	t423 = t461 * t411 + t412 * t458;
	t445 = qJD(2) * t457;
	t381 = t423 * qJD(3) + t411 * t445;
	t431 = t407 * t457 + t381;
	t453 = (t427 * t405 + t431 * t406) * r_i_i_C(1) + (-t431 * t405 + t427 * t406) * r_i_i_C(2);
	t426 = t433 - t463;
	t421 = t432 * t407 - t410 * t462;
	t372 = -t378 * t415 - qJD(3) * t460 + (t411 * t444 + t413 * t451) * t409;
	t1 = [-t380 * pkin(2) - t438 * qJD(4) - t429 * t373 + ((t391 * t405 - t406 * t438) * r_i_i_C(1) + (t391 * t406 + t405 * t438) * r_i_i_C(2)) * t407 - t426 * t379 + (-t464 * pkin(1) - pkin(8) * t459) * qJD(1) + t450 * t467 + (-t373 * t410 + (t391 * t410 - t414 * t438) * qJD(5)) * pkin(5), t465 * t377 - t426 * t378 + t418 * t393 + t421 * t394, -t450 * t371 + t422 * t372 + t420 * t428, t371, (t371 * t414 + t377 * t410 + (-t386 * t410 - t393 * t414) * qJD(5)) * pkin(5) + t455, t455; -t378 * pkin(2) + t368 * r_i_i_C(1) + t367 * r_i_i_C(2) + t371 * qJ(4) + t386 * qJD(4) + t463 * t377 + (-pkin(1) * t413 + pkin(8) * t449) * qJD(1) + t450 * t372 + (t371 * t410 + (t386 * t414 - t393 * t410) * qJD(5)) * pkin(5), -t379 * t465 + t426 * t380 + t418 * t391 + t421 * t392, -t450 * t373 + t420 * (t392 * t415 - t437) - t422 * t467, t373, (t373 * t414 - t379 * t410 + (-t391 * t414 - t410 * t438) * qJD(5)) * pkin(5) + t454, t454; 0, ((-qJD(2) * t465 + t421) * t412 + (t426 * qJD(2) - t418) * t416) * t409, -t450 * t381 + t420 * t423 + t422 * (-t389 * qJD(3) + t415 * t445), t381, (-t410 * t446 + t381 * t414 + (-t389 * t410 + t414 * t457) * qJD(5)) * pkin(5) + t453, t453;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end