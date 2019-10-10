% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:07
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
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:07
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (332->69), mult. (663->114), div. (0->0), fcn. (618->10), ass. (0->55)
	t280 = cos(pkin(6));
	t282 = sin(qJ(2));
	t283 = sin(qJ(1));
	t314 = t283 * t282;
	t302 = t280 * t314;
	t285 = cos(qJ(2));
	t286 = cos(qJ(1));
	t311 = t286 * t285;
	t266 = -qJD(1) * t302 - qJD(2) * t314 + (qJD(2) * t280 + qJD(1)) * t311;
	t277 = qJD(3) + qJD(4);
	t279 = sin(pkin(6));
	t315 = t279 * t286;
	t322 = t277 * t315 - t266;
	t321 = r_i_i_C(3) + pkin(10) + pkin(9);
	t320 = pkin(3) * qJD(3);
	t312 = t286 * t282;
	t313 = t283 * t285;
	t268 = t280 * t312 + t313;
	t319 = t268 * t277;
	t278 = qJ(3) + qJ(4);
	t275 = sin(t278);
	t318 = t275 * t277;
	t317 = t279 * t283;
	t284 = cos(qJ(3));
	t316 = t279 * t284;
	t276 = cos(t278);
	t292 = t302 - t311;
	t306 = qJD(1) * t286;
	t300 = t279 * t306;
	t290 = t277 * t292 + t300;
	t293 = t280 * t313 + t312;
	t264 = t268 * qJD(1) + t293 * qJD(2);
	t296 = t277 * t317 - t264;
	t259 = -t296 * t275 + t290 * t276;
	t260 = t290 * t275 + t296 * t276;
	t310 = t259 * r_i_i_C(1) - t260 * r_i_i_C(2);
	t307 = qJD(1) * t283;
	t301 = t279 * t307;
	t291 = t301 - t319;
	t298 = t322 * t276;
	t309 = (t322 * t275 + t291 * t276) * r_i_i_C(1) + (-t291 * t275 + t298) * r_i_i_C(2);
	t299 = qJD(2) * t279 * t285;
	t289 = -t277 * t280 - t299;
	t304 = t277 * t279 * t282;
	t308 = (t289 * t275 - t276 * t304) * r_i_i_C(1) + (t275 * t304 + t289 * t276) * r_i_i_C(2);
	t281 = sin(qJ(3));
	t305 = t281 * t320;
	t297 = -r_i_i_C(1) * t275 - r_i_i_C(2) * t276;
	t274 = t284 * pkin(3) + pkin(2);
	t295 = t276 * r_i_i_C(1) - t275 * r_i_i_C(2) + t274;
	t294 = t280 * t311 - t314;
	t288 = t297 * t277 - t305;
	t265 = t293 * qJD(1) + t268 * qJD(2);
	t263 = -t294 * qJD(1) + t292 * qJD(2);
	t1 = [(t268 * t318 + t298) * r_i_i_C(1) + (t266 * t275 + t276 * t319) * r_i_i_C(2) - t266 * t274 + t268 * t305 - pkin(1) * t306 - t321 * t265 + ((-r_i_i_C(2) * t318 + t284 * t320) * t286 + (-pkin(3) * t281 - pkin(8) + t297) * t307) * t279, t295 * t263 - t321 * t264 - t288 * t293, (t284 * t300 + t264 * t281 + (-t281 * t317 + t284 * t292) * qJD(3)) * pkin(3) + t310, t310, 0, 0; t260 * r_i_i_C(1) + t259 * r_i_i_C(2) - t264 * t274 - t321 * t263 + (-pkin(1) * t283 + pkin(8) * t315) * qJD(1) + (t281 * t300 + (t281 * t292 + t283 * t316) * qJD(3)) * pkin(3), -t295 * t265 + t321 * t266 + t288 * t294, (t284 * t301 - t266 * t281 + (-t268 * t284 + t281 * t315) * qJD(3)) * pkin(3) + t309, t309, 0, 0; 0, (t288 * t285 + (-t295 * t282 + t321 * t285) * qJD(2)) * t279, (-t281 * t299 + (-t280 * t281 - t282 * t316) * qJD(3)) * pkin(3) + t308, t308, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:07
	% EndTime: 2019-10-10 12:46:08
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (657->83), mult. (1178->130), div. (0->0), fcn. (1139->10), ass. (0->56)
	t383 = pkin(4) - r_i_i_C(2);
	t381 = r_i_i_C(3) + qJ(5);
	t342 = cos(pkin(6));
	t345 = sin(qJ(1));
	t344 = sin(qJ(2));
	t374 = t345 * t344;
	t367 = t342 * t374;
	t369 = qJD(2) * t344;
	t347 = cos(qJ(2));
	t348 = cos(qJ(1));
	t371 = t348 * t347;
	t321 = -qJD(1) * t367 - t345 * t369 + (qJD(2) * t342 + qJD(1)) * t371;
	t340 = qJ(3) + qJ(4);
	t337 = sin(t340);
	t338 = cos(t340);
	t339 = qJD(3) + qJD(4);
	t372 = t348 * t344;
	t373 = t345 * t347;
	t328 = t342 * t372 + t373;
	t341 = sin(pkin(6));
	t370 = qJD(1) * t341;
	t365 = t345 * t370;
	t356 = -t328 * t339 + t365;
	t375 = t341 * t348;
	t385 = t337 * (t339 * t375 - t321) + t338 * t356;
	t346 = cos(qJ(3));
	t336 = t346 * pkin(3) + pkin(2);
	t384 = t381 * t337 + t338 * t383 + t336;
	t382 = r_i_i_C(1) + pkin(10) + pkin(9);
	t357 = t367 - t371;
	t380 = t357 * t338;
	t379 = t337 * t339;
	t378 = t341 * t344;
	t377 = t341 * t345;
	t376 = t341 * t346;
	t368 = t338 * t378;
	t366 = t338 * t375;
	t364 = t348 * t370;
	t363 = qJD(2) * t341 * t347;
	t362 = -t321 * t338 + t339 * t366;
	t358 = t342 * t373 + t372;
	t319 = qJD(1) * t328 + qJD(2) * t358;
	t361 = t339 * t377 - t319;
	t359 = t342 * t371 - t374;
	t355 = t339 * t342 + t363;
	t307 = t337 * t361 - t338 * t364 - t339 * t380;
	t308 = t337 * t364 + t338 * t361 + t357 * t379;
	t354 = -(-t337 * t377 + t380) * qJD(5) + t381 * t308 - t383 * t307;
	t353 = -(-t328 * t338 + t337 * t375) * qJD(5) + t381 * (-t328 * t379 + t337 * t365 - t362) + t383 * t385;
	t316 = t337 * t355 + t339 * t368;
	t352 = -(-t342 * t337 - t368) * qJD(5) + t381 * (t338 * t355 - t378 * t379) - t383 * t316;
	t343 = sin(qJ(3));
	t350 = -qJD(3) * t343 * pkin(3) + t337 * qJD(5) + (-t337 * t383 + t381 * t338) * t339;
	t320 = qJD(1) * t358 + qJD(2) * t328;
	t318 = -qJD(1) * t359 + qJD(2) * t357;
	t1 = [-(t328 * t337 + t366) * qJD(5) - t321 * t336 - t382 * t320 + t383 * (-t337 * t356 + t362) + t381 * t385 + (-t348 * pkin(1) - pkin(8) * t377) * qJD(1) + (-t343 * t365 + (t328 * t343 + t346 * t375) * qJD(3)) * pkin(3), t384 * t318 - t382 * t319 - t350 * t358, (t346 * t364 + t319 * t343 + (-t343 * t377 + t346 * t357) * qJD(3)) * pkin(3) + t354, t354, t307, 0; -(t337 * t357 + t338 * t377) * qJD(5) - t319 * t336 - t382 * t318 + t383 * t308 + t381 * t307 + (-t345 * pkin(1) + pkin(8) * t375) * qJD(1) + (t343 * t364 + (t343 * t357 + t345 * t376) * qJD(3)) * pkin(3), -t320 * t384 + t321 * t382 + t350 * t359, (t346 * t365 - t321 * t343 + (-t328 * t346 + t343 * t375) * qJD(3)) * pkin(3) + t353, t353, -t385, 0; 0, (-t384 * t369 + (qJD(2) * t382 + t350) * t347) * t341, (-t343 * t363 + (-t342 * t343 - t344 * t376) * qJD(3)) * pkin(3) + t352, t352, t316, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:08
	% EndTime: 2019-10-10 12:46:09
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (1172->117), mult. (2134->191), div. (0->0), fcn. (2129->12), ass. (0->76)
	t442 = sin(qJ(6));
	t446 = cos(qJ(6));
	t464 = t446 * r_i_i_C(1) - t442 * r_i_i_C(2);
	t506 = t464 * qJD(6) + qJD(5);
	t463 = -r_i_i_C(1) * t442 - r_i_i_C(2) * t446;
	t507 = qJ(5) - t463;
	t479 = pkin(4) + pkin(11) + r_i_i_C(3);
	t444 = sin(qJ(2));
	t445 = sin(qJ(1));
	t448 = cos(qJ(2));
	t449 = cos(qJ(1));
	t494 = cos(pkin(6));
	t468 = t449 * t494;
	t421 = t444 * t468 + t445 * t448;
	t440 = qJ(3) + qJ(4);
	t438 = cos(t440);
	t441 = sin(pkin(6));
	t485 = t441 * t449;
	t430 = t438 * t485;
	t437 = sin(t440);
	t410 = t421 * t437 + t430;
	t465 = t448 * t468;
	t484 = t444 * t445;
	t420 = -t465 + t484;
	t504 = t410 * t442 + t420 * t446;
	t503 = -t410 * t446 + t420 * t442;
	t439 = qJD(3) + qJD(4);
	t443 = sin(qJ(3));
	t490 = t438 * t439;
	t451 = -(-t479 * t439 + t506) * t437 + qJD(3) * t443 * pkin(3) - t507 * t490;
	t447 = cos(qJ(3));
	t436 = pkin(3) * t447 + pkin(2);
	t501 = t437 * t507 + t479 * t438 + t436;
	t495 = -pkin(5) - pkin(10) - pkin(9);
	t491 = t437 * t439;
	t489 = t441 * t444;
	t488 = t441 * t445;
	t487 = t441 * t447;
	t486 = t441 * t448;
	t483 = t449 * t448;
	t482 = qJD(1) * t441;
	t481 = qJD(2) * t444;
	t480 = qJD(2) * t448;
	t476 = t437 * t489;
	t475 = t438 * t489;
	t474 = t437 * t485;
	t473 = t445 * t482;
	t472 = t449 * t482;
	t471 = t441 * t481;
	t470 = t441 * t480;
	t469 = t445 * t494;
	t460 = qJD(2) * t494 + qJD(1);
	t466 = t444 * t469;
	t409 = -qJD(1) * t466 - t445 * t481 + t460 * t483;
	t467 = -t409 * t438 + t439 * t430;
	t422 = t449 * t444 + t448 * t469;
	t407 = t421 * qJD(1) + t422 * qJD(2);
	t461 = t439 * t488 - t407;
	t458 = t463 * qJD(6);
	t457 = t464 - t495;
	t456 = t494 * t439 + t470;
	t394 = t409 * t437 + t421 * t490 - t438 * t473 - t439 * t474;
	t423 = -t466 + t483;
	t392 = t423 * t490 + t461 * t437 - t438 * t472;
	t393 = -t423 * t491 + t437 * t472 + t461 * t438;
	t454 = t506 * (t423 * t438 + t437 * t488) + t507 * t393 - t479 * t392;
	t453 = t506 * (t421 * t438 - t474) + t507 * (-t421 * t491 + t437 * t473 - t467) - t479 * t394;
	t402 = t456 * t437 + t439 * t475;
	t452 = t506 * (t494 * t437 + t475) + t507 * (t456 * t438 - t439 * t476) - t479 * t402;
	t418 = -t494 * t438 + t476;
	t413 = t423 * t437 - t438 * t488;
	t408 = t422 * qJD(1) + t421 * qJD(2);
	t406 = -qJD(1) * t465 - t449 * t480 + t460 * t484;
	t383 = t392 * t442 - t406 * t446 + (t413 * t446 - t422 * t442) * qJD(6);
	t382 = t392 * t446 + t406 * t442 + (-t413 * t442 - t422 * t446) * qJD(6);
	t1 = [-t410 * qJD(5) - t409 * t436 - t507 * t394 - t457 * t408 + (t503 * r_i_i_C(1) + t504 * r_i_i_C(2)) * qJD(6) + (-pkin(1) * t449 - pkin(8) * t488) * qJD(1) + t479 * ((t421 * t439 - t473) * t437 + t467) + (-t443 * t473 + (t421 * t443 + t447 * t485) * qJD(3)) * pkin(3), t501 * t406 - t457 * t407 + t451 * t422 + t423 * t458, (t447 * t472 + t407 * t443 + (-t423 * t447 - t443 * t488) * qJD(3)) * pkin(3) + t454, t454, t392, r_i_i_C(1) * t382 - r_i_i_C(2) * t383; t383 * r_i_i_C(1) + t382 * r_i_i_C(2) + t392 * qJ(5) + t413 * qJD(5) - t407 * t436 + t495 * t406 + (-pkin(1) * t445 + pkin(8) * t485) * qJD(1) + t479 * t393 + (t443 * t472 + (-t423 * t443 + t445 * t487) * qJD(3)) * pkin(3), -t408 * t501 + t457 * t409 + t451 * t420 + t421 * t458, (t447 * t473 - t409 * t443 + (-t421 * t447 + t443 * t485) * qJD(3)) * pkin(3) + t453, t453, t394, (t394 * t446 - t408 * t442) * r_i_i_C(1) + (-t394 * t442 - t408 * t446) * r_i_i_C(2) + (-t504 * r_i_i_C(1) + t503 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t501 + t458) * t444 + (t457 * qJD(2) - t451) * t448) * t441, (-t443 * t470 + (-t494 * t443 - t444 * t487) * qJD(3)) * pkin(3) + t452, t452, t402, (t402 * t446 - t442 * t471) * r_i_i_C(1) + (-t402 * t442 - t446 * t471) * r_i_i_C(2) + ((-t418 * t442 + t446 * t486) * r_i_i_C(1) + (-t418 * t446 - t442 * t486) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end