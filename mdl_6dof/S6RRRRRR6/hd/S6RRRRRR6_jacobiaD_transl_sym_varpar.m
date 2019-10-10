% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
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
	% StartTime: 2019-10-10 13:24:47
	% EndTime: 2019-10-10 13:24:47
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
	% StartTime: 2019-10-10 13:24:48
	% EndTime: 2019-10-10 13:24:48
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
	% StartTime: 2019-10-10 13:24:48
	% EndTime: 2019-10-10 13:24:48
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
	% StartTime: 2019-10-10 13:24:50
	% EndTime: 2019-10-10 13:24:51
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (953->124), mult. (1772->205), div. (0->0), fcn. (1760->12), ass. (0->81)
	t491 = pkin(11) + r_i_i_C(3);
	t435 = sin(qJ(5));
	t439 = cos(qJ(5));
	t453 = t439 * r_i_i_C(1) - t435 * r_i_i_C(2);
	t451 = pkin(4) + t453;
	t438 = sin(qJ(1));
	t434 = cos(pkin(6));
	t455 = qJD(2) * t434 + qJD(1);
	t437 = sin(qJ(2));
	t476 = t437 * t438;
	t462 = t434 * t476;
	t471 = qJD(2) * t437;
	t441 = cos(qJ(2));
	t442 = cos(qJ(1));
	t473 = t441 * t442;
	t406 = -qJD(1) * t462 - t438 * t471 + t455 * t473;
	t431 = qJD(3) + qJD(4);
	t433 = sin(pkin(6));
	t477 = t433 * t442;
	t500 = t431 * t477 - t406;
	t467 = qJD(5) * t439;
	t468 = qJD(5) * t435;
	t499 = -r_i_i_C(1) * t468 - t467 * r_i_i_C(2);
	t474 = t438 * t441;
	t475 = t437 * t442;
	t417 = t434 * t475 + t474;
	t432 = qJ(3) + qJ(4);
	t429 = sin(t432);
	t430 = cos(t432);
	t409 = -t417 * t430 + t429 * t477;
	t461 = t434 * t473;
	t416 = -t461 + t476;
	t498 = -t409 * t435 - t416 * t439;
	t497 = t409 * t439 - t416 * t435;
	t436 = sin(qJ(3));
	t496 = -qJD(3) * t436 * pkin(3) - (t451 * t429 - t491 * t430) * t431;
	t494 = t435 * r_i_i_C(1) + t439 * r_i_i_C(2);
	t440 = cos(qJ(3));
	t428 = pkin(3) * t440 + pkin(2);
	t492 = t491 * t429 + t451 * t430 + t428;
	t418 = t434 * t474 + t475;
	t405 = t418 * qJD(1) + t417 * qJD(2);
	t487 = t405 * t435;
	t486 = t405 * t439;
	t419 = -t462 + t473;
	t483 = t419 * t430;
	t482 = t429 * t431;
	t481 = t433 * t437;
	t480 = t433 * t438;
	t479 = t433 * t440;
	t478 = t433 * t441;
	t472 = qJD(1) * t433;
	t470 = qJD(2) * t441;
	t469 = qJD(5) * t430;
	t464 = t431 * t481;
	t460 = t438 * t472;
	t459 = t442 * t472;
	t458 = t433 * t471;
	t457 = t433 * t470;
	t456 = t500 * t430;
	t404 = t417 * qJD(1) + t418 * qJD(2);
	t452 = t431 * t480 - t404;
	t450 = -t417 * t431 + t460;
	t449 = t431 * t434 + t457;
	t394 = t500 * t429 + t450 * t430;
	t395 = -t417 * t482 + t429 * t460 - t456;
	t447 = t499 * (-t417 * t429 - t430 * t477) + t491 * t395 + t451 * t394;
	t392 = t452 * t429 - t430 * t459 + t431 * t483;
	t393 = -t419 * t482 + t429 * t459 + t452 * t430;
	t446 = t499 * (-t419 * t429 + t430 * t480) + t491 * t393 - t451 * t392;
	t402 = -t429 * t464 + t449 * t430;
	t445 = t499 * (-t429 * t481 + t430 * t434) + t491 * t402 + t451 * (-t449 * t429 - t430 * t464);
	t444 = t494 * t469 - t496;
	t443 = -pkin(10) - pkin(9);
	t415 = t429 * t434 + t430 * t481;
	t411 = t429 * t480 + t483;
	t403 = -qJD(1) * t461 - t442 * t470 + t455 * t476;
	t397 = -t450 * t429 + t456;
	t385 = t393 * t439 - t403 * t435 + (-t411 * t435 + t418 * t439) * qJD(5);
	t384 = -t393 * t435 - t403 * t439 + (-t411 * t439 - t418 * t435) * qJD(5);
	t1 = [(t397 * t439 - t487) * r_i_i_C(1) + (-t397 * t435 - t486) * r_i_i_C(2) + t397 * pkin(4) - t406 * t428 + t405 * t443 + t491 * t394 + (t498 * r_i_i_C(1) - t497 * r_i_i_C(2)) * qJD(5) + (-t442 * pkin(1) - pkin(8) * t480) * qJD(1) + (-t436 * t460 + (t417 * t436 + t440 * t477) * qJD(3)) * pkin(3), (-t404 * t435 + t419 * t467) * r_i_i_C(1) + (-t404 * t439 - t419 * t468) * r_i_i_C(2) + t404 * t443 + t492 * t403 + t444 * t418, (t440 * t459 + t404 * t436 + (-t419 * t440 - t436 * t480) * qJD(3)) * pkin(3) + t446, t446, r_i_i_C(1) * t384 - r_i_i_C(2) * t385, 0; t393 * pkin(4) + t385 * r_i_i_C(1) + t384 * r_i_i_C(2) + t403 * t443 - t404 * t428 + t491 * t392 + (-pkin(1) * t438 + pkin(8) * t477) * qJD(1) + (t436 * t459 + (-t419 * t436 + t438 * t479) * qJD(3)) * pkin(3), (t406 * t435 + t417 * t467) * r_i_i_C(1) + (t406 * t439 - t417 * t468) * r_i_i_C(2) - t406 * t443 - t492 * t405 + t444 * t416, (t440 * t460 - t406 * t436 + (-t417 * t440 + t436 * t477) * qJD(3)) * pkin(3) + t447, t447, (-t395 * t435 + t486) * r_i_i_C(1) + (-t395 * t439 - t487) * r_i_i_C(2) + (t497 * r_i_i_C(1) + t498 * r_i_i_C(2)) * qJD(5), 0; 0, ((-qJD(2) * t492 + t453 * qJD(5)) * t437 + (-qJD(2) * t443 + t494 * (qJD(2) - t469) + t496) * t441) * t433, (-t436 * t457 + (-t434 * t436 - t437 * t479) * qJD(3)) * pkin(3) + t445, t445, (-t402 * t435 + t439 * t458) * r_i_i_C(1) + (-t402 * t439 - t435 * t458) * r_i_i_C(2) + ((-t415 * t439 + t435 * t478) * r_i_i_C(1) + (t415 * t435 + t439 * t478) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:50
	% EndTime: 2019-10-10 13:24:51
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (1416->140), mult. (2311->224), div. (0->0), fcn. (2307->14), ass. (0->89)
	t477 = sin(qJ(1));
	t538 = cos(pkin(6));
	t494 = qJD(2) * t538 + qJD(1);
	t476 = sin(qJ(2));
	t509 = t477 * t538;
	t500 = t476 * t509;
	t521 = qJD(2) * t476;
	t480 = cos(qJ(2));
	t481 = cos(qJ(1));
	t526 = t481 * t480;
	t444 = -qJD(1) * t500 - t477 * t521 + t494 * t526;
	t472 = qJ(3) + qJ(4);
	t466 = sin(t472);
	t468 = cos(t472);
	t470 = qJD(3) + qJD(4);
	t508 = t481 * t538;
	t453 = t476 * t508 + t477 * t480;
	t473 = sin(pkin(6));
	t522 = qJD(1) * t473;
	t513 = t477 * t522;
	t492 = -t453 * t470 + t513;
	t528 = t473 * t481;
	t514 = t468 * t528;
	t432 = t444 * t468 + t492 * t466 - t470 * t514;
	t499 = t480 * t508;
	t527 = t477 * t476;
	t452 = -t499 + t527;
	t469 = qJD(5) + qJD(6);
	t505 = t452 * t469 + t432;
	t475 = sin(qJ(3));
	t471 = qJ(5) + qJ(6);
	t465 = sin(t471);
	t467 = cos(t471);
	t497 = r_i_i_C(1) * t465 + r_i_i_C(2) * t467;
	t474 = sin(qJ(5));
	t540 = t474 * pkin(5);
	t489 = -qJD(5) * t540 - t497 * t469;
	t539 = r_i_i_C(3) + pkin(12) + pkin(11);
	t478 = cos(qJ(5));
	t463 = t478 * pkin(5) + pkin(4);
	t498 = t467 * r_i_i_C(1) - t465 * r_i_i_C(2);
	t550 = t463 + t498;
	t484 = qJD(3) * t475 * pkin(3) + t466 * t470 * t550 - (t539 * t470 + t489) * t468;
	t454 = t481 * t476 + t480 * t509;
	t443 = t454 * qJD(1) + t453 * qJD(2);
	t447 = -t453 * t468 + t466 * t528;
	t547 = -t447 * t469 - t443;
	t455 = -t500 + t526;
	t512 = t481 * t522;
	t545 = -t455 * t470 + t512;
	t479 = cos(qJ(3));
	t464 = t479 * pkin(3) + pkin(2);
	t543 = t539 * t466 + t468 * t550 + t464;
	t535 = t465 * t469;
	t533 = t467 * t469;
	t532 = t473 * t476;
	t531 = t473 * t477;
	t530 = t473 * t479;
	t529 = t473 * t480;
	t520 = qJD(2) * t480;
	t441 = -qJD(1) * t499 - t481 * t520 + t494 * t527;
	t449 = t455 * t468 + t466 * t531;
	t504 = -t449 * t469 - t441;
	t442 = t453 * qJD(1) + t454 * qJD(2);
	t495 = t470 * t531 - t442;
	t430 = t545 * t466 + t495 * t468;
	t507 = t454 * t469 + t430;
	t417 = -t507 * t465 + t504 * t467;
	t418 = t504 * t465 + t507 * t467;
	t525 = t417 * r_i_i_C(1) - t418 * r_i_i_C(2);
	t524 = (-t505 * t465 - t467 * t547) * r_i_i_C(1) + (t465 * t547 - t505 * t467) * r_i_i_C(2);
	t515 = t468 * t532;
	t451 = t538 * t466 + t515;
	t511 = t473 * t521;
	t491 = -t451 * t469 + t511;
	t510 = t473 * t520;
	t490 = t538 * t470 + t510;
	t516 = t466 * t532;
	t440 = t490 * t468 - t470 * t516;
	t496 = t469 * t529 - t440;
	t523 = (t496 * t465 + t491 * t467) * r_i_i_C(1) + (-t491 * t465 + t496 * t467) * r_i_i_C(2);
	t519 = qJD(5) * t478;
	t431 = t492 * t468 + (t470 * t528 - t444) * t466;
	t429 = t495 * t466 - t545 * t468;
	t487 = t489 * (-t455 * t466 + t468 * t531) + t539 * t430 - t550 * t429;
	t486 = t489 * (-t453 * t466 - t514) + t539 * t432 + t550 * t431;
	t485 = t489 * (t538 * t468 - t516) + t539 * t440 + t550 * (-t490 * t466 - t470 * t515);
	t483 = -pkin(10) - pkin(9);
	t1 = [-t432 * t463 + t443 * t483 - t444 * t464 + (-t505 * r_i_i_C(1) + r_i_i_C(2) * t547) * t467 + (r_i_i_C(1) * t547 + t505 * r_i_i_C(2)) * t465 + t539 * t431 + (-t481 * pkin(1) - pkin(8) * t531) * qJD(1) + (-t443 * t474 + (-t447 * t474 - t452 * t478) * qJD(5)) * pkin(5) + (-t475 * t513 + (t453 * t475 + t479 * t528) * qJD(3)) * pkin(3), (-t442 * t465 + t455 * t533) * r_i_i_C(1) + (-t442 * t467 - t455 * t535) * r_i_i_C(2) + t442 * t483 + (-t442 * t474 + t455 * t519) * pkin(5) + t543 * t441 + t484 * t454, (t479 * t512 + t442 * t475 + (-t455 * t479 - t475 * t531) * qJD(3)) * pkin(3) + t487, t487, (-t430 * t474 - t441 * t478 + (-t449 * t478 - t454 * t474) * qJD(5)) * pkin(5) + t525, t525; t418 * r_i_i_C(1) + t417 * r_i_i_C(2) + t430 * t463 + t441 * t483 - t442 * t464 + t539 * t429 + (-pkin(1) * t477 + pkin(8) * t528) * qJD(1) + (-t441 * t474 + (-t449 * t474 + t454 * t478) * qJD(5)) * pkin(5) + (t475 * t512 + (-t455 * t475 + t477 * t530) * qJD(3)) * pkin(3), (t444 * t465 + t453 * t533) * r_i_i_C(1) + (t444 * t467 - t453 * t535) * r_i_i_C(2) - t444 * t483 + (t444 * t474 + t453 * t519) * pkin(5) - t543 * t443 + t484 * t452, (t479 * t513 - t444 * t475 + (-t453 * t479 + t475 * t528) * qJD(3)) * pkin(3) + t486, t486, (-t432 * t474 + t443 * t478 + (t447 * t478 - t452 * t474) * qJD(5)) * pkin(5) + t524, t524; 0, ((pkin(5) * t519 - qJD(2) * t543 + t498 * t469) * t476 + ((-t483 + t497 + t540) * qJD(2) - t484) * t480) * t473, (-t475 * t510 + (-t538 * t475 - t476 * t530) * qJD(3)) * pkin(3) + t485, t485, (t478 * t511 - t440 * t474 + (-t451 * t478 + t474 * t529) * qJD(5)) * pkin(5) + t523, t523;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end