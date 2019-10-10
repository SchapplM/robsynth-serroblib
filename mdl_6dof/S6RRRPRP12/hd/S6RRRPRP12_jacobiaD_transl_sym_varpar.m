% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP12
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
% Datum: 2019-10-10 11:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
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
	% StartTime: 2019-10-10 11:53:25
	% EndTime: 2019-10-10 11:53:25
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
	% StartTime: 2019-10-10 11:53:26
	% EndTime: 2019-10-10 11:53:26
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
	% StartTime: 2019-10-10 11:53:26
	% EndTime: 2019-10-10 11:53:27
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
	t318 = -r_i_i_C(1) - pkin(9);
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
	t1 = [-t302 * qJD(4) - t273 * pkin(2) + t318 * t272 + t319 * t322 + t317 * t321 + (-t291 * pkin(1) - pkin(8) * t315) * qJD(1), t320 * t270 + t318 * t271 - t292 * t297, t300 * qJD(4) - t319 * t264 + t317 * t265, t264, 0, 0; -(t284 * t310 + t316) * qJD(4) - t271 * pkin(2) + t318 * t270 + t319 * t265 + t317 * t264 + (-t288 * pkin(1) + pkin(8) * t313) * qJD(1), -t272 * t320 - t318 * t273 + t292 * t298, -t301 * qJD(4) - t317 * t322 + t319 * t321, -t321, 0, 0; 0, (t292 * t290 + (-t287 * t320 - t318 * t290) * qJD(2)) * t284, t299 * qJD(4) + t317 * (t289 * t303 + (-t284 * t286 * t287 + t285 * t289) * qJD(3)) - t319 * t274, t274, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:27
	% EndTime: 2019-10-10 11:53:28
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
	% StartTime: 2019-10-10 11:53:29
	% EndTime: 2019-10-10 11:53:30
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (894->122), mult. (2628->192), div. (0->0), fcn. (2705->10), ass. (0->77)
	t466 = sin(qJ(1));
	t523 = cos(pkin(6));
	t482 = qJD(2) * t523 + qJD(1);
	t465 = sin(qJ(2));
	t494 = t466 * t523;
	t490 = t465 * t494;
	t509 = qJD(2) * t465;
	t469 = cos(qJ(2));
	t470 = cos(qJ(1));
	t510 = t470 * t469;
	t436 = -qJD(1) * t490 - t466 * t509 + t482 * t510;
	t493 = t470 * t523;
	t449 = t465 * t493 + t466 * t469;
	t464 = sin(qJ(3));
	t462 = sin(pkin(6));
	t468 = cos(qJ(3));
	t515 = t462 * t468;
	t500 = qJD(1) * t515;
	t512 = t464 * t470;
	t502 = t462 * t512;
	t507 = qJD(3) * t468;
	t425 = -qJD(3) * t502 + t436 * t464 + t449 * t507 - t466 * t500;
	t450 = t470 * t465 + t469 * t494;
	t435 = t450 * qJD(1) + t449 * qJD(2);
	t463 = sin(qJ(5));
	t467 = cos(qJ(5));
	t489 = t469 * t493;
	t511 = t465 * t466;
	t448 = -t489 + t511;
	t513 = t462 * t470;
	t531 = t449 * t464 + t468 * t513;
	t532 = -t448 * t463 + t467 * t531;
	t537 = t532 * qJD(5) + t425 * t463 + t435 * t467;
	t486 = t448 * t467 + t463 * t531;
	t536 = -t486 * qJD(5) + t425 * t467 - t435 * t463;
	t491 = -t467 * qJD(6) + qJD(4);
	t504 = r_i_i_C(2) + pkin(10) + pkin(3);
	t471 = -qJ(4) * t507 + (t504 * qJD(3) - t491) * t464;
	t530 = (qJD(2) * t464 + qJD(5)) * t465 - t469 * t507;
	t529 = t464 * qJ(4) + t504 * t468 + pkin(2);
	t516 = t462 * t466;
	t503 = t464 * t516;
	t528 = -qJD(1) * t503 + qJD(3) * t531 - t436 * t468;
	t526 = pkin(4) + pkin(9);
	t525 = r_i_i_C(1) + pkin(5);
	t524 = r_i_i_C(3) + qJ(6);
	t451 = -t490 + t510;
	t517 = t451 * t464;
	t514 = t462 * t469;
	t508 = qJD(2) * t469;
	t506 = qJD(5) * t464;
	t505 = t463 * qJD(6);
	t501 = t467 * t514;
	t499 = t462 * t509;
	t498 = t462 * t508;
	t434 = t449 * qJD(1) + t450 * qJD(2);
	t488 = -t450 * t506 - t434;
	t487 = -t448 * t506 + t436;
	t443 = -t466 * t515 + t517;
	t484 = t443 * t467 - t450 * t463;
	t483 = t443 * t463 + t450 * t467;
	t481 = (qJD(2) + t506) * t469;
	t479 = t451 * t468 + t503;
	t446 = t462 * t465 * t464 - t523 * t468;
	t477 = t523 * t464 + t465 * t515;
	t475 = t525 * t463 - t524 * t467 + qJ(4);
	t433 = -qJD(1) * t489 - t470 * t508 + t482 * t511;
	t474 = qJD(5) * t451 - t433 * t464 + t450 * t507;
	t473 = qJD(5) * t449 + t435 * t464 + t448 * t507;
	t472 = (t524 * t463 + t525 * t467) * qJD(5) + t491;
	t438 = t477 * qJD(3) + t464 * t498;
	t429 = -t438 * t467 - qJD(5) * t501 + (qJD(5) * t446 + t499) * t463;
	t424 = -t434 * t468 - qJD(3) * t517 + (qJD(1) * t512 + t466 * t507) * t462;
	t423 = t479 * qJD(3) - t434 * t464 - t470 * t500;
	t412 = t484 * qJD(5) + t423 * t463 - t433 * t467;
	t411 = t483 * qJD(5) - t423 * t467 - t433 * t463;
	t1 = [t532 * qJD(6) - t425 * qJ(4) - t531 * qJD(4) - t436 * pkin(2) - t526 * t435 - t525 * t537 + t524 * t536 + (-pkin(1) * t470 - pkin(8) * t516) * qJD(1) + t504 * t528, t451 * t505 - t526 * t434 + t525 * (-t474 * t463 + t488 * t467) + t524 * (t488 * t463 + t474 * t467) + t471 * t450 + t529 * t433, -t504 * t423 + t475 * t424 + t472 * t479, t423, t483 * qJD(6) - t525 * t411 + t524 * t412, t411; -t484 * qJD(6) + t423 * qJ(4) + t443 * qJD(4) - t434 * pkin(2) - t526 * t433 + t525 * t412 + t524 * t411 + (-pkin(1) * t466 + pkin(8) * t513) * qJD(1) + t504 * t424, t449 * t505 + t526 * t436 + t525 * (-t473 * t463 + t487 * t467) + t524 * (t487 * t463 + t473 * t467) + t471 * t448 - t529 * t435, -t504 * t425 - t475 * t528 + t472 * (t449 * t468 - t502), t425, t486 * qJD(6) + t524 * t537 + t525 * t536, -t536; 0, (t525 * (-t530 * t463 + t467 * t481) + t524 * (t463 * t481 + t530 * t467) + t465 * t505 - t471 * t469 + (-t465 * t529 + t526 * t469) * qJD(2)) * t462, -t504 * t438 + t475 * (-t446 * qJD(3) + t468 * t498) + t472 * t477, t438, -(-t446 * t463 + t501) * qJD(6) + t524 * (t467 * t499 + t438 * t463 + (t446 * t467 + t463 * t514) * qJD(5)) - t525 * t429, t429;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end