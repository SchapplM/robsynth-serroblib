% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:06
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:06:07
	% EndTime: 2019-10-10 12:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:06:07
	% EndTime: 2019-10-10 12:06:07
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
	% StartTime: 2019-10-10 12:06:08
	% EndTime: 2019-10-10 12:06:08
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
	% StartTime: 2019-10-10 12:06:09
	% EndTime: 2019-10-10 12:06:09
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
	% StartTime: 2019-10-10 12:06:09
	% EndTime: 2019-10-10 12:06:09
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (252->77), mult. (593->125), div. (0->0), fcn. (554->10), ass. (0->52)
	t254 = qJ(3) + pkin(12);
	t252 = sin(t254);
	t253 = cos(t254);
	t258 = sin(qJ(3));
	t292 = pkin(3) * t258;
	t265 = r_i_i_C(1) * t252 + r_i_i_C(2) * t253 + t292;
	t264 = qJD(3) * t265;
	t291 = t252 * r_i_i_C(2);
	t261 = cos(qJ(3));
	t290 = t261 * pkin(3);
	t289 = r_i_i_C(3) + qJ(4) + pkin(9);
	t255 = sin(pkin(6));
	t259 = sin(qJ(2));
	t288 = t255 * t259;
	t260 = sin(qJ(1));
	t287 = t255 * t260;
	t286 = t255 * t261;
	t263 = cos(qJ(1));
	t285 = t255 * t263;
	t284 = t260 * t259;
	t262 = cos(qJ(2));
	t283 = t260 * t262;
	t282 = t263 * t259;
	t281 = t263 * t262;
	t280 = qJD(1) * t260;
	t279 = qJD(1) * t263;
	t278 = qJD(2) * t259;
	t277 = qJD(2) * t262;
	t256 = cos(pkin(6));
	t242 = t256 * t282 + t283;
	t276 = qJD(3) * t242;
	t275 = qJD(3) * t263;
	t274 = t256 * t284;
	t273 = t256 * t281;
	t272 = t255 * t280;
	t271 = t255 * t279;
	t270 = t255 * t275;
	t269 = qJD(2) * t256 + qJD(1);
	t251 = pkin(2) + t290;
	t268 = t253 * r_i_i_C(1) + t251 - t291;
	t267 = t256 * t283 + t282;
	t266 = t272 - t276;
	t245 = t253 * t270;
	t244 = -t274 + t281;
	t241 = -t273 + t284;
	t240 = -qJD(1) * t274 - t260 * t278 + t269 * t281;
	t239 = t267 * qJD(1) + t242 * qJD(2);
	t238 = t242 * qJD(1) + t267 * qJD(2);
	t237 = -qJD(1) * t273 - t263 * t277 + t269 * t284;
	t236 = t252 * t271 - t238 * t253 + (-t244 * t252 + t253 * t287) * qJD(3);
	t235 = t253 * t271 + t238 * t252 + (-t244 * t253 - t252 * t287) * qJD(3);
	t1 = [(-t240 * t253 + t252 * t276 + t245) * r_i_i_C(1) + (t240 * t252 + t253 * t276) * r_i_i_C(2) - t240 * t251 + t276 * t292 - t241 * qJD(4) - pkin(1) * t279 - t289 * t239 + ((t290 - t291) * t275 + (-pkin(8) - t265) * t280) * t255, t244 * qJD(4) + t268 * t237 - t289 * t238 + t264 * t267, t235 * r_i_i_C(1) - t236 * r_i_i_C(2) + (t261 * t271 + t238 * t258 + (-t244 * t261 - t258 * t287) * qJD(3)) * pkin(3), -t237, 0, 0; t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t267 * qJD(4) - t238 * t251 - t289 * t237 + (-pkin(1) * t260 + pkin(8) * t285) * qJD(1) + (t258 * t271 + (-t244 * t258 + t260 * t286) * qJD(3)) * pkin(3), t242 * qJD(4) - t268 * t239 + t289 * t240 + t241 * t264, t245 * r_i_i_C(2) + (t266 * r_i_i_C(1) - t240 * r_i_i_C(2)) * t253 + ((-t240 + t270) * r_i_i_C(1) - t266 * r_i_i_C(2)) * t252 + (t261 * t272 - t240 * t258 + (-t242 * t261 + t258 * t285) * qJD(3)) * pkin(3), t239, 0, 0; 0, (qJD(4) * t259 - t262 * t264 + (-t268 * t259 + t289 * t262) * qJD(2)) * t255, -t265 * t255 * t277 + ((-t252 * t256 - t253 * t288) * r_i_i_C(1) + (t252 * t288 - t253 * t256) * r_i_i_C(2) + (-t256 * t258 - t259 * t286) * pkin(3)) * qJD(3), t255 * t278, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:06:10
	% EndTime: 2019-10-10 12:06:11
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (686->123), mult. (1483->200), div. (0->0), fcn. (1476->12), ass. (0->71)
	t393 = sin(qJ(1));
	t388 = cos(pkin(6));
	t409 = qJD(2) * t388 + qJD(1);
	t392 = sin(qJ(2));
	t427 = t393 * t392;
	t417 = t388 * t427;
	t422 = qJD(2) * t392;
	t396 = cos(qJ(2));
	t397 = cos(qJ(1));
	t424 = t397 * t396;
	t362 = -qJD(1) * t417 - t393 * t422 + t409 * t424;
	t425 = t397 * t392;
	t426 = t393 * t396;
	t373 = t388 * t425 + t426;
	t386 = qJ(3) + pkin(12);
	t384 = sin(t386);
	t385 = cos(t386);
	t387 = sin(pkin(6));
	t423 = qJD(1) * t387;
	t414 = t393 * t423;
	t428 = t387 * t397;
	t416 = t385 * t428;
	t356 = (-qJD(3) * t373 + t414) * t384 - qJD(3) * t416 + t362 * t385;
	t374 = t388 * t426 + t425;
	t361 = t374 * qJD(1) + t373 * qJD(2);
	t390 = sin(qJ(5));
	t394 = cos(qJ(5));
	t447 = t356 * t390 - t361 * t394;
	t446 = -t356 * t394 - t361 * t390;
	t391 = sin(qJ(3));
	t407 = t394 * r_i_i_C(1) - t390 * r_i_i_C(2);
	t405 = pkin(4) + t407;
	t440 = pkin(10) + r_i_i_C(3);
	t445 = (t391 * pkin(3) + t405 * t384 - t440 * t385) * qJD(3);
	t367 = -t373 * t385 + t384 * t428;
	t415 = t388 * t424;
	t372 = -t415 + t427;
	t444 = -t367 * t390 - t372 * t394;
	t443 = t367 * t394 - t372 * t390;
	t406 = t390 * r_i_i_C(1) + t394 * r_i_i_C(2);
	t395 = cos(qJ(3));
	t383 = t395 * pkin(3) + pkin(2);
	t442 = t440 * t384 + t405 * t385 + t383;
	t432 = t387 * t392;
	t431 = t387 * t393;
	t430 = t387 * t395;
	t429 = t387 * t396;
	t421 = qJD(2) * t396;
	t420 = qJD(5) * t385;
	t419 = qJD(5) * t390;
	t418 = qJD(5) * t394;
	t413 = t397 * t423;
	t412 = t387 * t421;
	t411 = t387 * t422;
	t375 = -t417 + t424;
	t404 = -t375 * t384 + t385 * t431;
	t369 = t375 * t385 + t384 * t431;
	t371 = t388 * t384 + t385 * t432;
	t403 = -t384 * t432 + t388 * t385;
	t402 = qJD(5) * t406;
	t399 = t367 * qJD(3) - t362 * t384 + t385 * t414;
	t398 = t406 * t420 + t445;
	t389 = -qJ(4) - pkin(9);
	t364 = t403 * qJD(3) + t385 * t412;
	t360 = t373 * qJD(1) + t374 * qJD(2);
	t359 = -qJD(1) * t415 - t397 * t421 + t409 * t427;
	t354 = t404 * qJD(3) - t360 * t385 + t384 * t413;
	t353 = t369 * qJD(3) - t360 * t384 - t385 * t413;
	t352 = t354 * t394 - t359 * t390 + (-t369 * t390 + t374 * t394) * qJD(5);
	t351 = -t354 * t390 - t359 * t394 + (-t369 * t394 - t374 * t390) * qJD(5);
	t1 = [t446 * r_i_i_C(1) + t447 * r_i_i_C(2) - t356 * pkin(4) - t362 * t383 + t361 * t389 - t372 * qJD(4) + t440 * t399 + (t444 * r_i_i_C(1) - t443 * r_i_i_C(2)) * qJD(5) + (-t397 * pkin(1) - pkin(8) * t431) * qJD(1) + (-t391 * t414 + (t373 * t391 + t395 * t428) * qJD(3)) * pkin(3), (-t360 * t390 + t375 * t418) * r_i_i_C(1) + (-t360 * t394 - t375 * t419) * r_i_i_C(2) + t360 * t389 + t375 * qJD(4) + t442 * t359 + t398 * t374, t440 * t354 - t404 * t402 - t405 * t353 + (t395 * t413 + t360 * t391 + (-t375 * t395 - t391 * t431) * qJD(3)) * pkin(3), -t359, t351 * r_i_i_C(1) - t352 * r_i_i_C(2), 0; t354 * pkin(4) + t352 * r_i_i_C(1) + t351 * r_i_i_C(2) + t374 * qJD(4) + t359 * t389 - t360 * t383 + t440 * t353 + (-pkin(1) * t393 + pkin(8) * t428) * qJD(1) + (t391 * t413 + (-t375 * t391 + t393 * t430) * qJD(3)) * pkin(3), (t362 * t390 + t373 * t418) * r_i_i_C(1) + (t362 * t394 - t373 * t419) * r_i_i_C(2) - t362 * t389 + t373 * qJD(4) - t442 * t361 + t398 * t372, t440 * t356 - (-t373 * t384 - t416) * t402 + t405 * t399 + (t395 * t414 - t362 * t391 + (-t373 * t395 + t391 * t428) * qJD(3)) * pkin(3), t361, -t447 * r_i_i_C(1) + t446 * r_i_i_C(2) + (t443 * r_i_i_C(1) + t444 * r_i_i_C(2)) * qJD(5), 0; 0, ((-qJD(2) * t442 + t407 * qJD(5) + qJD(4)) * t392 + (-qJD(2) * t389 - t445 + t406 * (qJD(2) - t420)) * t396) * t387, t440 * t364 - t403 * t402 + t405 * (-t371 * qJD(3) - t384 * t412) + (-t391 * t412 + (-t388 * t391 - t392 * t430) * qJD(3)) * pkin(3), t411, (-t364 * t390 + t394 * t411) * r_i_i_C(1) + (-t364 * t394 - t390 * t411) * r_i_i_C(2) + ((-t371 * t394 + t390 * t429) * r_i_i_C(1) + (t371 * t390 + t394 * t429) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:06:11
	% EndTime: 2019-10-10 12:06:12
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (1093->146), mult. (1995->230), div. (0->0), fcn. (1994->14), ass. (0->82)
	t432 = sin(qJ(1));
	t486 = cos(pkin(6));
	t447 = qJD(2) * t486 + qJD(1);
	t431 = sin(qJ(2));
	t462 = t432 * t486;
	t452 = t431 * t462;
	t471 = qJD(2) * t431;
	t435 = cos(qJ(2));
	t436 = cos(qJ(1));
	t476 = t436 * t435;
	t397 = -qJD(1) * t452 - t432 * t471 + t447 * t476;
	t461 = t436 * t486;
	t408 = t431 * t461 + t432 * t435;
	t425 = qJ(3) + pkin(12);
	t420 = sin(t425);
	t421 = cos(t425);
	t427 = sin(pkin(6));
	t472 = qJD(1) * t427;
	t466 = t432 * t472;
	t478 = t427 * t436;
	t467 = t421 * t478;
	t391 = (-qJD(3) * t408 + t466) * t420 - qJD(3) * t467 + t397 * t421;
	t451 = t435 * t461;
	t477 = t432 * t431;
	t407 = -t451 + t477;
	t424 = qJD(5) + qJD(6);
	t458 = t407 * t424 + t391;
	t430 = sin(qJ(3));
	t433 = cos(qJ(5));
	t418 = t433 * pkin(5) + pkin(4);
	t426 = qJ(5) + qJ(6);
	t422 = sin(t426);
	t423 = cos(t426);
	t450 = r_i_i_C(1) * t423 - r_i_i_C(2) * t422;
	t446 = t418 + t450;
	t487 = r_i_i_C(3) + pkin(11) + pkin(10);
	t496 = (t430 * pkin(3) + t446 * t420 - t487 * t421) * qJD(3);
	t409 = t436 * t431 + t435 * t462;
	t396 = t409 * qJD(1) + t408 * qJD(2);
	t402 = -t408 * t421 + t420 * t478;
	t495 = -t402 * t424 - t396;
	t449 = r_i_i_C(1) * t422 + r_i_i_C(2) * t423;
	t429 = sin(qJ(5));
	t489 = t429 * pkin(5);
	t494 = qJD(5) * t489 + t449 * t424;
	t434 = cos(qJ(3));
	t419 = t434 * pkin(3) + pkin(2);
	t493 = t487 * t420 + t446 * t421 + t419;
	t484 = t422 * t424;
	t483 = t423 * t424;
	t482 = t427 * t431;
	t481 = t427 * t432;
	t480 = t427 * t434;
	t479 = t427 * t435;
	t470 = qJD(2) * t435;
	t394 = -qJD(1) * t451 - t436 * t470 + t447 * t477;
	t410 = -t452 + t476;
	t404 = t410 * t421 + t420 * t481;
	t457 = -t404 * t424 - t394;
	t395 = t408 * qJD(1) + t409 * qJD(2);
	t445 = -t410 * t420 + t421 * t481;
	t465 = t436 * t472;
	t389 = t445 * qJD(3) - t395 * t421 + t420 * t465;
	t460 = t409 * t424 + t389;
	t384 = -t460 * t422 + t457 * t423;
	t385 = t457 * t422 + t460 * t423;
	t475 = t384 * r_i_i_C(1) - t385 * r_i_i_C(2);
	t474 = (-t458 * t422 - t423 * t495) * r_i_i_C(1) + (t422 * t495 - t458 * t423) * r_i_i_C(2);
	t406 = t486 * t420 + t421 * t482;
	t463 = t427 * t471;
	t444 = -t406 * t424 + t463;
	t442 = -t420 * t482 + t486 * t421;
	t464 = t427 * t470;
	t399 = t442 * qJD(3) + t421 * t464;
	t448 = t424 * t479 - t399;
	t473 = (t448 * t422 + t444 * t423) * r_i_i_C(1) + (-t444 * t422 + t448 * t423) * r_i_i_C(2);
	t469 = qJD(5) * t433;
	t439 = t402 * qJD(3) - t397 * t420 + t421 * t466;
	t438 = t494 * t421 + t496;
	t428 = -qJ(4) - pkin(9);
	t388 = t404 * qJD(3) - t395 * t420 - t421 * t465;
	t1 = [-t407 * qJD(4) - t391 * t418 + t396 * t428 - t397 * t419 + (-t458 * r_i_i_C(1) + r_i_i_C(2) * t495) * t423 + (r_i_i_C(1) * t495 + t458 * r_i_i_C(2)) * t422 + t487 * t439 + (-t436 * pkin(1) - pkin(8) * t481) * qJD(1) + (-t396 * t429 + (-t402 * t429 - t407 * t433) * qJD(5)) * pkin(5) + (-t430 * t466 + (t408 * t430 + t434 * t478) * qJD(3)) * pkin(3), (-t395 * t422 + t410 * t483) * r_i_i_C(1) + (-t395 * t423 - t410 * t484) * r_i_i_C(2) + t395 * t428 + t410 * qJD(4) + (-t395 * t429 + t410 * t469) * pkin(5) + t493 * t394 + t438 * t409, t487 * t389 - t494 * t445 - t446 * t388 + (t434 * t465 + t395 * t430 + (-t410 * t434 - t430 * t481) * qJD(3)) * pkin(3), -t394, (-t389 * t429 - t394 * t433 + (-t404 * t433 - t409 * t429) * qJD(5)) * pkin(5) + t475, t475; t385 * r_i_i_C(1) + t384 * r_i_i_C(2) + t409 * qJD(4) + t389 * t418 + t394 * t428 - t395 * t419 + t487 * t388 + (-t432 * pkin(1) + pkin(8) * t478) * qJD(1) + (-t394 * t429 + (-t404 * t429 + t409 * t433) * qJD(5)) * pkin(5) + (t430 * t465 + (-t410 * t430 + t432 * t480) * qJD(3)) * pkin(3), (t397 * t422 + t408 * t483) * r_i_i_C(1) + (t397 * t423 - t408 * t484) * r_i_i_C(2) - t397 * t428 + t408 * qJD(4) + (t397 * t429 + t408 * t469) * pkin(5) - t493 * t396 + t438 * t407, t487 * t391 - t494 * (-t408 * t420 - t467) + t446 * t439 + (t434 * t466 - t397 * t430 + (-t408 * t434 + t430 * t478) * qJD(3)) * pkin(3), t396, (-t391 * t429 + t396 * t433 + (t402 * t433 - t407 * t429) * qJD(5)) * pkin(5) + t474, t474; 0, ((pkin(5) * t469 - qJD(2) * t493 + t450 * t424 + qJD(4)) * t431 + (-qJD(2) * t428 + (-qJD(5) * t421 + qJD(2)) * t489 - t496 + t449 * (-t421 * t424 + qJD(2))) * t435) * t427, t487 * t399 - t494 * t442 + t446 * (-t406 * qJD(3) - t420 * t464) + (-t430 * t464 + (-t486 * t430 - t431 * t480) * qJD(3)) * pkin(3), t463, (t433 * t463 - t399 * t429 + (-t406 * t433 + t429 * t479) * qJD(5)) * pkin(5) + t473, t473;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end