% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (103->31), mult. (307->49), div. (0->0), fcn. (282->8), ass. (0->26)
	t199 = sin(pkin(12));
	t200 = sin(pkin(6));
	t201 = cos(pkin(12));
	t220 = t200 * (r_i_i_C(1) * t199 + r_i_i_C(2) * t201 + pkin(8));
	t219 = r_i_i_C(3) + qJ(3);
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t218 = t204 * t203;
	t205 = cos(qJ(2));
	t217 = t204 * t205;
	t206 = cos(qJ(1));
	t216 = t206 * t203;
	t215 = t206 * t205;
	t214 = qJD(2) * t203;
	t202 = cos(pkin(6));
	t213 = t202 * t218;
	t212 = t202 * t215;
	t211 = qJD(2) * t202 + qJD(1);
	t210 = t201 * r_i_i_C(1) - t199 * r_i_i_C(2) + pkin(2);
	t208 = t202 * t217 + t216;
	t207 = t202 * t216 + t217;
	t194 = -qJD(1) * t213 - t204 * t214 + t211 * t215;
	t193 = t208 * qJD(1) + t207 * qJD(2);
	t192 = t207 * qJD(1) + t208 * qJD(2);
	t191 = -qJD(1) * t212 - qJD(2) * t215 + t211 * t218;
	t1 = [-(-t212 + t218) * qJD(3) - t219 * t193 - t210 * t194 + (-t206 * pkin(1) - t204 * t220) * qJD(1), -(t213 - t215) * qJD(3) - t219 * t192 + t210 * t191, -t191, 0, 0, 0; t208 * qJD(3) - t219 * t191 - t210 * t192 + (-t204 * pkin(1) + t206 * t220) * qJD(1), t207 * qJD(3) - t210 * t193 + t219 * t194, t193, 0, 0, 0; 0, (t203 * qJD(3) + (-t210 * t203 + t219 * t205) * qJD(2)) * t200, t200 * t214, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (222->61), mult. (494->100), div. (0->0), fcn. (467->10), ass. (0->45)
	t259 = sin(qJ(1));
	t256 = cos(pkin(6));
	t267 = qJD(2) * t256 + qJD(1);
	t258 = sin(qJ(2));
	t280 = t259 * t258;
	t272 = t256 * t280;
	t275 = qJD(2) * t258;
	t260 = cos(qJ(2));
	t261 = cos(qJ(1));
	t277 = t261 * t260;
	t239 = -qJD(1) * t272 - t259 * t275 + t267 * t277;
	t255 = sin(pkin(6));
	t281 = t255 * t261;
	t268 = qJD(4) * t281;
	t285 = t239 - t268;
	t253 = pkin(12) + qJ(4);
	t251 = sin(t253);
	t252 = cos(t253);
	t266 = r_i_i_C(1) * t251 + r_i_i_C(2) * t252;
	t263 = qJD(4) * t266;
	t284 = r_i_i_C(3) + pkin(9) + qJ(3);
	t283 = t255 * t258;
	t282 = t255 * t259;
	t279 = t259 * t260;
	t278 = t261 * t258;
	t276 = qJD(1) * t255;
	t274 = qJD(2) * t260;
	t241 = t256 * t278 + t279;
	t273 = qJD(4) * t241;
	t271 = t256 * t277;
	t270 = sin(pkin(12)) * pkin(3) + pkin(8);
	t269 = t261 * t276;
	t250 = cos(pkin(12)) * pkin(3) + pkin(2);
	t265 = t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + t250;
	t264 = t256 * t279 + t278;
	t262 = t259 * t276 - t273;
	t244 = t252 * t268;
	t243 = -t272 + t277;
	t240 = -t271 + t280;
	t238 = t264 * qJD(1) + t241 * qJD(2);
	t237 = t241 * qJD(1) + t264 * qJD(2);
	t236 = -qJD(1) * t271 - t261 * t274 + t267 * t280;
	t235 = t251 * t269 - t237 * t252 + (-t243 * t251 + t252 * t282) * qJD(4);
	t234 = t252 * t269 + t237 * t251 + (-t243 * t252 - t251 * t282) * qJD(4);
	t1 = [(-t239 * t252 + t251 * t273 + t244) * r_i_i_C(1) + (t285 * t251 + t252 * t273) * r_i_i_C(2) - t239 * t250 - t240 * qJD(3) - t284 * t238 + (-t261 * pkin(1) + (-t266 - t270) * t282) * qJD(1), t243 * qJD(3) + t265 * t236 - t284 * t237 + t263 * t264, -t236, t234 * r_i_i_C(1) - t235 * r_i_i_C(2), 0, 0; t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t264 * qJD(3) - t237 * t250 - t284 * t236 + (-pkin(1) * t259 + t270 * t281) * qJD(1), t241 * qJD(3) - t265 * t238 + t284 * t239 + t240 * t263, t238, t244 * r_i_i_C(2) + (t262 * r_i_i_C(1) - t239 * r_i_i_C(2)) * t252 + (-t285 * r_i_i_C(1) - t262 * r_i_i_C(2)) * t251, 0, 0; 0, (qJD(3) * t258 - t260 * t263 + (-t265 * t258 + t284 * t260) * qJD(2)) * t255, t255 * t275, -t266 * t255 * t274 + ((-t251 * t256 - t252 * t283) * r_i_i_C(1) + (t251 * t283 - t252 * t256) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:15
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (456->78), mult. (710->123), div. (0->0), fcn. (669->12), ass. (0->61)
	t292 = sin(qJ(1));
	t290 = cos(pkin(6));
	t303 = qJD(2) * t290 + qJD(1);
	t291 = sin(qJ(2));
	t323 = t292 * t291;
	t309 = t290 * t323;
	t314 = qJD(2) * t291;
	t293 = cos(qJ(2));
	t294 = cos(qJ(1));
	t320 = t294 * t293;
	t268 = -qJD(1) * t309 - t292 * t314 + t303 * t320;
	t288 = qJD(4) + qJD(5);
	t289 = sin(pkin(6));
	t324 = t289 * t294;
	t333 = t288 * t324 - t268;
	t287 = pkin(12) + qJ(4);
	t285 = qJ(5) + t287;
	t281 = sin(t285);
	t282 = cos(t285);
	t302 = r_i_i_C(1) * t281 + r_i_i_C(2) * t282;
	t283 = sin(t287);
	t329 = pkin(4) * qJD(4);
	t312 = t283 * t329;
	t332 = t302 * t288 + t312;
	t331 = pkin(8) + pkin(4) * t283 + sin(pkin(12)) * pkin(3);
	t330 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(3);
	t321 = t294 * t291;
	t322 = t292 * t293;
	t271 = t290 * t321 + t322;
	t328 = t271 * t288;
	t327 = t281 * t288;
	t326 = t289 * t291;
	t325 = t289 * t292;
	t273 = -t309 + t320;
	t315 = qJD(1) * t294;
	t306 = t289 * t315;
	t297 = -t273 * t288 + t306;
	t299 = t290 * t322 + t321;
	t266 = t271 * qJD(1) + t299 * qJD(2);
	t301 = t288 * t325 - t266;
	t261 = -t301 * t281 + t297 * t282;
	t262 = t297 * t281 + t301 * t282;
	t319 = t261 * r_i_i_C(1) - t262 * r_i_i_C(2);
	t316 = qJD(1) * t292;
	t307 = t289 * t316;
	t298 = t307 - t328;
	t304 = t333 * t282;
	t318 = (t333 * t281 + t298 * t282) * r_i_i_C(1) + (-t298 * t281 + t304) * r_i_i_C(2);
	t313 = qJD(2) * t293;
	t305 = t289 * t313;
	t296 = -t288 * t290 - t305;
	t311 = t288 * t326;
	t317 = (t296 * t281 - t282 * t311) * r_i_i_C(1) + (t281 * t311 + t296 * t282) * r_i_i_C(2);
	t308 = t290 * t320;
	t284 = cos(t287);
	t274 = pkin(4) * t284 + cos(pkin(12)) * pkin(3) + pkin(2);
	t300 = t282 * r_i_i_C(1) - t281 * r_i_i_C(2) + t274;
	t270 = -t308 + t323;
	t267 = t299 * qJD(1) + t271 * qJD(2);
	t265 = -qJD(1) * t308 - t294 * t313 + t303 * t323;
	t1 = [(t271 * t327 + t304) * r_i_i_C(1) + (t268 * t281 + t282 * t328) * r_i_i_C(2) - t268 * t274 + t271 * t312 - t270 * qJD(3) - pkin(1) * t315 - t330 * t267 + ((-r_i_i_C(2) * t327 + t284 * t329) * t294 + (-t302 - t331) * t316) * t289, t273 * qJD(3) + t300 * t265 - t330 * t266 + t299 * t332, -t265, (t284 * t306 + t266 * t283 + (-t273 * t284 - t283 * t325) * qJD(4)) * pkin(4) + t319, t319, 0; t262 * r_i_i_C(1) + t261 * r_i_i_C(2) + t299 * qJD(3) - t266 * t274 - t330 * t265 + (-t273 * t283 + t284 * t325) * t329 + (-pkin(1) * t292 + t331 * t324) * qJD(1), t271 * qJD(3) - t300 * t267 + t330 * t268 + t332 * t270, t267, (t284 * t307 - t268 * t283 + (-t271 * t284 + t283 * t324) * qJD(4)) * pkin(4) + t318, t318, 0; 0, (qJD(3) * t291 - t332 * t293 + (-t300 * t291 + t330 * t293) * qJD(2)) * t289, t289 * t314, (-t283 * t305 + (-t283 * t290 - t284 * t326) * qJD(4)) * pkin(4) + t317, t317, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:16
	% EndTime: 2019-10-10 11:03:17
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1257->132), mult. (1819->207), div. (0->0), fcn. (1811->14), ass. (0->83)
	t496 = pkin(11) + r_i_i_C(3);
	t442 = sin(qJ(6));
	t445 = cos(qJ(6));
	t457 = t445 * r_i_i_C(1) - t442 * r_i_i_C(2);
	t455 = pkin(5) + t457;
	t444 = sin(qJ(1));
	t441 = cos(pkin(6));
	t459 = qJD(2) * t441 + qJD(1);
	t443 = sin(qJ(2));
	t480 = t443 * t444;
	t466 = t441 * t480;
	t475 = qJD(2) * t443;
	t446 = cos(qJ(2));
	t447 = cos(qJ(1));
	t477 = t446 * t447;
	t408 = -qJD(1) * t466 - t444 * t475 + t459 * t477;
	t439 = qJD(4) + qJD(5);
	t440 = sin(pkin(6));
	t481 = t440 * t447;
	t505 = t439 * t481 - t408;
	t471 = qJD(6) * t445;
	t472 = qJD(6) * t442;
	t504 = -r_i_i_C(1) * t472 - t471 * r_i_i_C(2);
	t478 = t444 * t446;
	t479 = t443 * t447;
	t421 = t441 * t479 + t478;
	t438 = pkin(12) + qJ(4);
	t436 = qJ(5) + t438;
	t432 = sin(t436);
	t433 = cos(t436);
	t411 = -t421 * t433 + t432 * t481;
	t465 = t441 * t477;
	t420 = -t465 + t480;
	t503 = -t411 * t442 - t420 * t445;
	t502 = t411 * t445 - t420 * t442;
	t434 = sin(t438);
	t491 = pkin(4) * qJD(4);
	t501 = -t434 * t491 - (t455 * t432 - t496 * t433) * t439;
	t499 = t442 * r_i_i_C(1) + t445 * r_i_i_C(2);
	t435 = cos(t438);
	t425 = pkin(4) * t435 + cos(pkin(12)) * pkin(3) + pkin(2);
	t497 = t496 * t432 + t455 * t433 + t425;
	t492 = pkin(8) + pkin(4) * t434 + sin(pkin(12)) * pkin(3);
	t422 = t441 * t478 + t479;
	t407 = t422 * qJD(1) + t421 * qJD(2);
	t490 = t407 * t442;
	t489 = t407 * t445;
	t423 = -t466 + t477;
	t486 = t423 * t433;
	t485 = t432 * t439;
	t484 = t440 * t443;
	t483 = t440 * t444;
	t482 = t440 * t446;
	t476 = qJD(1) * t440;
	t474 = qJD(2) * t446;
	t473 = qJD(6) * t433;
	t468 = t439 * t484;
	t464 = t444 * t476;
	t463 = t447 * t476;
	t462 = t440 * t474;
	t461 = t440 * t475;
	t460 = t505 * t433;
	t406 = t421 * qJD(1) + t422 * qJD(2);
	t456 = t439 * t483 - t406;
	t454 = -t421 * t439 + t464;
	t453 = t439 * t441 + t462;
	t396 = t505 * t432 + t454 * t433;
	t404 = -t432 * t468 + t453 * t433;
	t451 = t504 * (-t432 * t484 + t433 * t441) + t496 * t404 + t455 * (-t453 * t432 - t433 * t468);
	t397 = -t421 * t485 + t432 * t464 - t460;
	t450 = t504 * (-t421 * t432 - t433 * t481) + t496 * t397 + t455 * t396;
	t394 = t456 * t432 - t433 * t463 + t439 * t486;
	t395 = -t423 * t485 + t432 * t463 + t456 * t433;
	t449 = t504 * (-t423 * t432 + t433 * t483) + t496 * t395 - t455 * t394;
	t448 = t499 * t473 - t501;
	t437 = -pkin(10) - pkin(9) - qJ(3);
	t417 = t432 * t441 + t433 * t484;
	t413 = t432 * t483 + t486;
	t405 = -qJD(1) * t465 - t447 * t474 + t459 * t480;
	t399 = -t454 * t432 + t460;
	t387 = t395 * t445 - t405 * t442 + (-t413 * t442 + t422 * t445) * qJD(6);
	t386 = -t395 * t442 - t405 * t445 + (-t413 * t445 - t422 * t442) * qJD(6);
	t1 = [(t399 * t445 - t490) * r_i_i_C(1) + (-t399 * t442 - t489) * r_i_i_C(2) + t399 * pkin(5) - t408 * t425 + t407 * t437 - t420 * qJD(3) + t496 * t396 + (t503 * r_i_i_C(1) - t502 * r_i_i_C(2)) * qJD(6) + (t421 * t434 + t435 * t481) * t491 + (-t447 * pkin(1) - t492 * t483) * qJD(1), (-t406 * t442 + t423 * t471) * r_i_i_C(1) + (-t406 * t445 - t423 * t472) * r_i_i_C(2) + t406 * t437 + t423 * qJD(3) + t497 * t405 + t448 * t422, -t405, (t435 * t463 + t406 * t434 + (-t423 * t435 - t434 * t483) * qJD(4)) * pkin(4) + t449, t449, r_i_i_C(1) * t386 - r_i_i_C(2) * t387; t395 * pkin(5) + t387 * r_i_i_C(1) + t386 * r_i_i_C(2) + t422 * qJD(3) + t405 * t437 - t406 * t425 + t496 * t394 + (-t423 * t434 + t435 * t483) * t491 + (-pkin(1) * t444 + t492 * t481) * qJD(1), (t408 * t442 + t421 * t471) * r_i_i_C(1) + (t408 * t445 - t421 * t472) * r_i_i_C(2) - t408 * t437 + t421 * qJD(3) - t497 * t407 + t448 * t420, t407, (t435 * t464 - t408 * t434 + (-t421 * t435 + t434 * t481) * qJD(4)) * pkin(4) + t450, t450, (-t397 * t442 + t489) * r_i_i_C(1) + (-t397 * t445 - t490) * r_i_i_C(2) + (t502 * r_i_i_C(1) + t503 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t497 + t457 * qJD(6) + qJD(3)) * t443 + (-qJD(2) * t437 + t499 * (qJD(2) - t473) + t501) * t446) * t440, t461, (-t434 * t462 + (-t434 * t441 - t435 * t484) * qJD(4)) * pkin(4) + t451, t451, (-t404 * t442 + t445 * t461) * r_i_i_C(1) + (-t404 * t445 - t442 * t461) * r_i_i_C(2) + ((-t417 * t445 + t442 * t482) * r_i_i_C(1) + (t417 * t442 + t445 * t482) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end