% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(6));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(6));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(12));
	t113 = sin(pkin(12));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(7));
	t278 = r_i_i_C(3) + pkin(9);
	t279 = t278 * t248 + qJ(2);
	t250 = sin(qJ(3));
	t277 = t250 * r_i_i_C(2);
	t246 = sin(pkin(6));
	t251 = sin(qJ(1));
	t276 = t246 * t251;
	t253 = cos(qJ(1));
	t275 = t246 * t253;
	t274 = t248 * t250;
	t252 = cos(qJ(3));
	t273 = t248 * t252;
	t244 = sin(pkin(12));
	t272 = t251 * t244;
	t247 = cos(pkin(12));
	t271 = t251 * t247;
	t270 = t253 * t244;
	t269 = t253 * t247;
	t268 = qJD(1) * t251;
	t267 = qJD(1) * t253;
	t245 = sin(pkin(7));
	t266 = qJD(3) * t245;
	t265 = t278 * t245;
	t249 = cos(pkin(6));
	t264 = t249 * t272;
	t263 = t246 * t268;
	t262 = t246 * t267;
	t261 = t266 * t275;
	t260 = r_i_i_C(1) * t250 + r_i_i_C(2) * t252;
	t236 = -qJD(1) * t264 + t247 * t267;
	t237 = -t249 * t269 + t272;
	t259 = qJD(3) * t237 * t248 - t236;
	t258 = t260 * t245;
	t256 = t249 * t271 + t270;
	t257 = t245 * t276 - t248 * t256;
	t238 = t249 * t270 + t271;
	t233 = t237 * qJD(1);
	t255 = t233 * t248 + t245 * t262;
	t235 = t256 * qJD(1);
	t254 = -qJD(3) * t238 - t235 * t248 + t245 * t263;
	t241 = t252 * t261;
	t240 = -t264 + t269;
	t234 = t238 * qJD(1);
	t232 = -t234 * t252 + t255 * t250 + (-t240 * t250 + t257 * t252) * qJD(3);
	t231 = t234 * t250 + t255 * t252 + (-t240 * t252 - t257 * t250) * qJD(3);
	t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0, 0; 0, 0, (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (208->86), mult. (693->152), div. (0->0), fcn. (698->12), ass. (0->61)
	t302 = sin(pkin(13));
	t306 = cos(pkin(13));
	t312 = cos(qJ(3));
	t329 = qJD(3) * t312;
	t310 = sin(qJ(3));
	t330 = qJD(3) * t310;
	t340 = t302 * t330 - t306 * t329;
	t339 = pkin(9) + qJ(4);
	t304 = sin(pkin(7));
	t305 = sin(pkin(6));
	t338 = t304 * t305;
	t308 = cos(pkin(7));
	t337 = t308 * t310;
	t303 = sin(pkin(12));
	t311 = sin(qJ(1));
	t336 = t311 * t303;
	t307 = cos(pkin(12));
	t335 = t311 * t307;
	t313 = cos(qJ(1));
	t334 = t313 * t303;
	t333 = t313 * t307;
	t332 = qJD(1) * t311;
	t331 = qJD(1) * t313;
	t328 = pkin(3) * t330;
	t327 = pkin(3) * t329;
	t309 = cos(pkin(6));
	t326 = t309 * t336;
	t325 = t305 * t332;
	t324 = t305 * t331;
	t319 = t312 * t302 + t310 * t306;
	t275 = t319 * t304;
	t321 = t304 * t310 * pkin(3) + r_i_i_C(1) * t275 + t339 * t308 + qJ(2);
	t270 = t340 * t304;
	t320 = -t270 * r_i_i_C(1) + t308 * qJD(4) + t304 * t327 + qJD(2);
	t292 = t310 * t302 - t312 * t306;
	t284 = -t309 * t333 + t336;
	t318 = t309 * t335 + t334;
	t285 = t309 * t334 + t335;
	t317 = qJD(3) * t319;
	t272 = t340 * t308;
	t281 = -qJD(1) * t326 + t307 * t331;
	t289 = -t302 * t329 - t306 * t330;
	t316 = -t284 * t272 + t281 * t292 - t285 * t289;
	t273 = t308 * t317;
	t288 = t292 * qJD(3);
	t315 = -t284 * t273 + t281 * t319 - t285 * t288;
	t277 = t319 * t308;
	t278 = t284 * qJD(1);
	t279 = t285 * qJD(1);
	t287 = -t326 + t333;
	t314 = -t272 * t318 - t278 * t277 - t279 * t292 - t287 * t289;
	t301 = t312 * pkin(3) + pkin(2);
	t291 = -t304 * qJD(4) + t308 * t327;
	t283 = pkin(3) * t337 - t339 * t304;
	t280 = t318 * qJD(1);
	t276 = t292 * t308;
	t274 = t292 * t304;
	t271 = t304 * t317;
	t269 = -t278 * t304 + t308 * t324;
	t268 = t318 * t273 - t278 * t276 + t279 * t319 + t287 * t288 + (-t271 * t311 - t274 * t331) * t305;
	t1 = [t316 * r_i_i_C(1) + t315 * r_i_i_C(2) - t281 * t301 + t285 * t328 + t284 * t291 - pkin(1) * t331 + (t277 * r_i_i_C(1) - t276 * r_i_i_C(2) - t304 * r_i_i_C(3) + t283) * t280 + ((-t271 * r_i_i_C(2) + t320) * t313 + (r_i_i_C(2) * t274 - r_i_i_C(3) * t308 - t321) * t332) * t305, t324, t268 * r_i_i_C(1) + ((t270 * t311 - t275 * t331) * t305 + t314) * r_i_i_C(2) + (t279 * t310 + (t278 * t308 + t304 * t324) * t312 + (-t287 * t312 + (t308 * t318 - t311 * t338) * t310) * qJD(3)) * pkin(3), t269, 0, 0; -t314 * r_i_i_C(1) + t268 * r_i_i_C(2) + t269 * r_i_i_C(3) - t279 * t301 - t287 * t328 + t278 * t283 - t318 * t291 - pkin(1) * t332 + (t320 * t311 + t321 * t331) * t305, t325, (t280 * t276 - t315) * r_i_i_C(1) + (t280 * t277 + t316) * r_i_i_C(2) + ((t271 * t313 - t274 * t332) * r_i_i_C(1) + (-t270 * t313 - t275 * t332) * r_i_i_C(2)) * t305 + (-t281 * t310 + (-t280 * t308 + t304 * t325) * t312 + (-t285 * t312 + (t284 * t308 + t313 * t338) * t310) * qJD(3)) * pkin(3), t280 * t304 + t308 * t325, 0, 0; 0, 0, (-t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t304 * t328) * t309 + ((-t273 * t307 + t288 * t303) * r_i_i_C(1) + (t272 * t307 - t289 * t303) * r_i_i_C(2) + (-t303 * t312 - t307 * t337) * qJD(3) * pkin(3)) * t305, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:13
	% EndTime: 2019-10-10 01:00:14
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (657->123), mult. (2123->219), div. (0->0), fcn. (2305->14), ass. (0->84)
	t457 = sin(pkin(7));
	t455 = sin(pkin(13));
	t466 = cos(qJ(3));
	t503 = cos(pkin(13));
	t481 = qJD(3) * t503;
	t463 = sin(qJ(3));
	t490 = qJD(3) * t463;
	t507 = t455 * t490 - t466 * t481;
	t422 = t507 * t457;
	t460 = cos(pkin(7));
	t424 = t507 * t460;
	t446 = -t466 * t455 - t463 * t503;
	t427 = t446 * t457;
	t429 = t446 * t460;
	t461 = cos(pkin(6));
	t456 = sin(pkin(12));
	t467 = cos(qJ(1));
	t495 = t467 * t456;
	t459 = cos(pkin(12));
	t464 = sin(qJ(1));
	t496 = t464 * t459;
	t476 = t461 * t496 + t495;
	t432 = t476 * qJD(1);
	t497 = t464 * t456;
	t485 = t461 * t497;
	t491 = qJD(1) * t467;
	t433 = -qJD(1) * t485 + t459 * t491;
	t494 = t467 * t459;
	t437 = -t461 * t494 + t497;
	t438 = t461 * t495 + t496;
	t489 = qJD(3) * t466;
	t442 = -t455 * t489 - t463 * t481;
	t458 = sin(pkin(6));
	t472 = -t463 * t455 + t466 * t503;
	t492 = qJD(1) * t464;
	t402 = t437 * t424 + t432 * t429 + t433 * t472 + t438 * t442 + (t422 * t467 - t427 * t492) * t458;
	t484 = t458 * t492;
	t417 = t432 * t457 + t460 * t484;
	t462 = sin(qJ(5));
	t465 = cos(qJ(5));
	t513 = t402 * t462 - t417 * t465;
	t512 = -t402 * t465 - t417 * t462;
	t500 = t458 * t467;
	t407 = -t427 * t500 - t437 * t429 - t438 * t472;
	t418 = -t437 * t457 + t460 * t500;
	t509 = t407 * t465 + t418 * t462;
	t508 = -t407 * t462 + t418 * t465;
	t430 = t437 * qJD(1);
	t431 = t438 * qJD(1);
	t440 = -t485 + t494;
	t399 = t476 * t424 - t430 * t429 - t431 * t472 + t440 * t442 + (-t422 * t464 - t427 * t491) * t458;
	t506 = t461 * t422 + (t424 * t459 - t442 * t456) * t458;
	t505 = -r_i_i_C(3) - pkin(10);
	t504 = pkin(9) + qJ(4);
	t502 = t457 * t463;
	t501 = t458 * t464;
	t499 = t460 * t463;
	t493 = pkin(3) * t502 + t504 * t460 + qJ(2);
	t486 = pkin(3) * t489;
	t488 = t460 * qJD(4) + t457 * t486 + qJD(2);
	t487 = pkin(3) * t490;
	t483 = t458 * t491;
	t479 = -t462 * r_i_i_C(1) - t465 * r_i_i_C(2);
	t477 = t465 * r_i_i_C(1) - t462 * r_i_i_C(2) + pkin(4);
	t475 = qJD(5) * t479;
	t471 = qJD(3) * t446;
	t423 = t457 * t471;
	t425 = t460 * t471;
	t426 = t472 * t457;
	t428 = t472 * t460;
	t441 = t472 * qJD(3);
	t468 = t423 * t500 + t437 * t425 - t426 * t484 + t432 * t428 - t433 * t446 + t438 * t441;
	t454 = t466 * pkin(3) + pkin(2);
	t444 = -t457 * qJD(4) + t460 * t486;
	t436 = -t458 * t459 * t457 + t461 * t460;
	t435 = pkin(3) * t499 - t504 * t457;
	t420 = t457 * t476 + t460 * t501;
	t415 = -t430 * t457 + t460 * t483;
	t414 = -t461 * t427 + (-t429 * t459 + t456 * t472) * t458;
	t409 = -t427 * t501 + t429 * t476 + t440 * t472;
	t398 = -t476 * t425 + t430 * t428 - t431 * t446 - t440 * t441 + (t423 * t464 + t426 * t491) * t458;
	t396 = t399 * t465 + t415 * t462 + (-t409 * t462 + t420 * t465) * qJD(5);
	t395 = -t399 * t462 + t415 * t465 + (-t409 * t465 - t420 * t462) * qJD(5);
	t1 = [t512 * r_i_i_C(1) + t513 * r_i_i_C(2) - t402 * pkin(4) - t433 * t454 + t438 * t487 + t432 * t435 + t437 * t444 - pkin(1) * t491 + t505 * t468 + (t508 * r_i_i_C(1) - t509 * r_i_i_C(2)) * qJD(5) + (t488 * t467 - t493 * t492) * t458, t483, -t505 * t399 + (t426 * t501 - t428 * t476 + t440 * t446) * t475 + t477 * t398 + (t431 * t463 + (t430 * t460 + t457 * t483) * t466 + (-t440 * t466 + (-t457 * t501 + t460 * t476) * t463) * qJD(3)) * pkin(3), t415, t395 * r_i_i_C(1) - t396 * r_i_i_C(2), 0; -t440 * t487 - pkin(1) * t492 + t399 * pkin(4) + t396 * r_i_i_C(1) + t395 * r_i_i_C(2) + t430 * t435 - t431 * t454 - t476 * t444 + t505 * t398 + (t488 * t464 + t493 * t491) * t458, t484, -t505 * t402 + (-t426 * t500 - t437 * t428 + t438 * t446) * t475 - t477 * t468 + (-t433 * t463 + (-t432 * t460 + t457 * t484) * t466 + (-t438 * t466 + (t437 * t460 + t457 * t500) * t463) * qJD(3)) * pkin(3), t417, -t513 * r_i_i_C(1) + t512 * r_i_i_C(2) + (t509 * r_i_i_C(1) + t508 * r_i_i_C(2)) * qJD(5), 0; 0, 0, t505 * t506 + (t461 * t426 + (t428 * t459 + t446 * t456) * t458) * t475 + t477 * (t461 * t423 + (t425 * t459 - t441 * t456) * t458) + (-t461 * t502 + (-t456 * t466 - t459 * t499) * t458) * qJD(3) * pkin(3), 0, -t479 * t506 + ((-t414 * t465 - t436 * t462) * r_i_i_C(1) + (t414 * t462 - t436 * t465) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:17
	% EndTime: 2019-10-10 01:00:19
	% DurationCPUTime: 2.14s
	% Computational Cost: add. (1858->170), mult. (5833->294), div. (0->0), fcn. (6652->16), ass. (0->106)
	t631 = sin(pkin(7));
	t629 = sin(pkin(13));
	t642 = cos(qJ(3));
	t690 = cos(pkin(13));
	t666 = qJD(3) * t690;
	t638 = sin(qJ(3));
	t677 = qJD(3) * t638;
	t693 = t629 * t677 - t642 * t666;
	t596 = t693 * t631;
	t620 = -t642 * t629 - t638 * t690;
	t601 = t620 * t631;
	t632 = sin(pkin(6));
	t643 = cos(qJ(1));
	t634 = cos(pkin(7));
	t598 = t693 * t634;
	t603 = t620 * t634;
	t635 = cos(pkin(6));
	t630 = sin(pkin(12));
	t682 = t643 * t630;
	t633 = cos(pkin(12));
	t639 = sin(qJ(1));
	t683 = t639 * t633;
	t659 = t635 * t683 + t682;
	t606 = t659 * qJD(1);
	t684 = t639 * t630;
	t670 = t635 * t684;
	t678 = qJD(1) * t643;
	t607 = -qJD(1) * t670 + t633 * t678;
	t681 = t643 * t633;
	t611 = -t635 * t681 + t684;
	t612 = t635 * t682 + t683;
	t676 = qJD(3) * t642;
	t616 = -t629 * t676 - t638 * t666;
	t653 = -t638 * t629 + t642 * t690;
	t647 = t611 * t598 + t606 * t603 + t607 * t653 + t612 * t616;
	t679 = qJD(1) * t639;
	t562 = t632 * (t596 * t643 - t601 * t679) + t647;
	t669 = t632 * t679;
	t590 = t606 * t631 + t634 * t669;
	t637 = sin(qJ(5));
	t641 = cos(qJ(5));
	t687 = t632 * t643;
	t578 = -t601 * t687 - t611 * t603 - t612 * t653;
	t591 = -t611 * t631 + t634 * t687;
	t664 = t578 * t637 - t591 * t641;
	t554 = qJD(5) * t664 + t562 * t641 + t590 * t637;
	t651 = qJD(3) * t620;
	t597 = t631 * t651;
	t599 = t634 * t651;
	t600 = t653 * t631;
	t602 = t653 * t634;
	t615 = t653 * qJD(3);
	t561 = -t611 * t599 - t606 * t602 + t607 * t620 - t612 * t615 + t632 * (-t597 * t643 + t600 * t679);
	t636 = sin(qJ(6));
	t640 = cos(qJ(6));
	t706 = t554 * t636 + t561 * t640;
	t705 = -t554 * t640 + t561 * t636;
	t569 = t578 * t641 + t591 * t637;
	t658 = -t600 * t687 - t611 * t602 + t612 * t620;
	t704 = -t569 * t636 + t640 * t658;
	t703 = t569 * t640 + t636 * t658;
	t702 = qJD(5) * t569 - t562 * t637 + t590 * t641;
	t604 = t611 * qJD(1);
	t605 = t612 * qJD(1);
	t614 = -t670 + t681;
	t559 = t659 * t598 - t604 * t603 - t605 * t653 + t614 * t616 + t632 * (-t596 * t639 - t601 * t678);
	t584 = t635 * t596 + (t598 * t633 - t616 * t630) * t632;
	t692 = pkin(11) + r_i_i_C(3);
	t691 = pkin(9) + qJ(4);
	t689 = t631 * t638;
	t688 = t632 * t639;
	t686 = t634 * t638;
	t680 = pkin(3) * t689 + t634 * t691 + qJ(2);
	t675 = qJD(6) * t636;
	t674 = qJD(6) * t640;
	t671 = pkin(3) * t676;
	t673 = t634 * qJD(4) + t631 * t671 + qJD(2);
	t672 = pkin(3) * t677;
	t668 = t632 * t678;
	t593 = t631 * t659 + t634 * t688;
	t650 = -t601 * t688 + t603 * t659 + t614 * t653;
	t571 = t593 * t637 + t641 * t650;
	t663 = t593 * t641 - t637 * t650;
	t610 = -t632 * t633 * t631 + t635 * t634;
	t649 = -t635 * t601 + (-t603 * t633 + t630 * t653) * t632;
	t573 = t610 * t637 + t641 * t649;
	t662 = t610 * t641 - t637 * t649;
	t660 = t640 * r_i_i_C(1) - t636 * r_i_i_C(2) + pkin(5);
	t657 = qJD(6) * (-t636 * r_i_i_C(1) - t640 * r_i_i_C(2));
	t652 = -t604 * t631 + t634 * t668;
	t645 = t692 * t637 + t660 * t641 + pkin(4);
	t644 = t641 * t657 + (-t660 * t637 + t641 * t692) * qJD(5);
	t628 = t642 * pkin(3) + pkin(2);
	t618 = -t631 * qJD(4) + t634 * t671;
	t609 = pkin(3) * t686 - t631 * t691;
	t586 = t635 * t600 + (t602 * t633 + t620 * t630) * t632;
	t582 = t635 * t597 + (t599 * t633 - t615 * t630) * t632;
	t580 = t600 * t688 - t602 * t659 + t614 * t620;
	t566 = qJD(5) * t662 - t584 * t641;
	t560 = -t596 * t687 + t601 * t669 - t647;
	t558 = -t659 * t599 + t604 * t602 - t605 * t620 - t614 * t615 + (t597 * t639 + t600 * t678) * t632;
	t552 = qJD(5) * t663 + t559 * t641 + t637 * t652;
	t551 = qJD(5) * t571 + t559 * t637 - t641 * t652;
	t550 = t552 * t640 - t558 * t636 + (-t571 * t636 - t580 * t640) * qJD(6);
	t549 = -t552 * t636 - t558 * t640 + (-t571 * t640 + t580 * t636) * qJD(6);
	t1 = [t705 * r_i_i_C(1) + t706 * r_i_i_C(2) - t554 * pkin(5) - t562 * pkin(4) + t561 * pkin(10) - t607 * t628 + t612 * t672 + t606 * t609 + t611 * t618 - pkin(1) * t678 + t692 * t702 + (t704 * r_i_i_C(1) - t703 * r_i_i_C(2)) * qJD(6) + (t643 * t673 - t679 * t680) * t632, t668, (t559 * t636 + t650 * t674) * r_i_i_C(1) + (t559 * t640 - t650 * t675) * r_i_i_C(2) + t559 * pkin(10) + t645 * t558 + t644 * t580 + (t605 * t638 + (t604 * t634 + t631 * t668) * t642 + (-t614 * t642 + (-t631 * t688 + t634 * t659) * t638) * qJD(3)) * pkin(3), t652, -t660 * t551 + t552 * t692 + t663 * t657, t549 * r_i_i_C(1) - t550 * r_i_i_C(2); -t614 * t672 - pkin(1) * t679 + t559 * pkin(4) + t552 * pkin(5) - t558 * pkin(10) + t550 * r_i_i_C(1) + t549 * r_i_i_C(2) + t604 * t609 - t605 * t628 - t659 * t618 + t692 * t551 + (t639 * t673 + t678 * t680) * t632, t669, (-t560 * t636 - t578 * t674) * r_i_i_C(1) + (-t560 * t640 + t578 * t675) * r_i_i_C(2) - t560 * pkin(10) + t645 * t561 + t644 * t658 + (-t607 * t638 + (-t606 * t634 + t631 * t669) * t642 + (-t612 * t642 + (t611 * t634 + t631 * t687) * t638) * qJD(3)) * pkin(3), t590, t692 * t554 + t664 * t657 + t660 * t702, -t706 * r_i_i_C(1) + t705 * r_i_i_C(2) + (t703 * r_i_i_C(1) + t704 * r_i_i_C(2)) * qJD(6); 0, 0, (-t584 * t636 + t649 * t674) * r_i_i_C(1) + (-t584 * t640 - t649 * t675) * r_i_i_C(2) - t584 * pkin(10) + (-t635 * t689 + (-t630 * t642 - t633 * t686) * t632) * qJD(3) * pkin(3) + t645 * t582 + t644 * t586, 0, t692 * t566 + t662 * t657 + t660 * (-qJD(5) * t573 + t584 * t637), (-t566 * t636 - t582 * t640) * r_i_i_C(1) + (-t566 * t640 + t582 * t636) * r_i_i_C(2) + ((-t573 * t640 + t586 * t636) * r_i_i_C(1) + (t573 * t636 + t586 * t640) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end