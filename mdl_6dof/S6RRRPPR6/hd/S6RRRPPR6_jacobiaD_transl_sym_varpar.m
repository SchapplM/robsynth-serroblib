% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:57
	% EndTime: 2019-10-10 11:25:57
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
	% StartTime: 2019-10-10 11:25:57
	% EndTime: 2019-10-10 11:25:57
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (252->77), mult. (593->125), div. (0->0), fcn. (554->10), ass. (0->52)
	t254 = qJ(3) + pkin(11);
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
	% StartTime: 2019-10-10 11:25:58
	% EndTime: 2019-10-10 11:25:58
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (459->86), mult. (980->134), div. (0->0), fcn. (944->10), ass. (0->53)
	t301 = qJ(3) + pkin(11);
	t299 = sin(t301);
	t300 = cos(t301);
	t305 = sin(qJ(3));
	t339 = r_i_i_C(3) + qJ(5);
	t342 = pkin(4) - r_i_i_C(2);
	t347 = (pkin(3) * t305 + t342 * t299 - t339 * t300) * qJD(3) - t299 * qJD(5);
	t307 = sin(qJ(1));
	t303 = cos(pkin(6));
	t320 = qJD(2) * t303 + qJD(1);
	t306 = sin(qJ(2));
	t333 = t307 * t306;
	t325 = t303 * t333;
	t328 = qJD(2) * t306;
	t309 = cos(qJ(2));
	t310 = cos(qJ(1));
	t330 = t310 * t309;
	t284 = -qJD(1) * t325 - t307 * t328 + t320 * t330;
	t331 = t310 * t306;
	t332 = t307 * t309;
	t289 = t303 * t331 + t332;
	t302 = sin(pkin(6));
	t334 = t302 * t310;
	t319 = t289 * t299 + t300 * t334;
	t329 = qJD(1) * t302;
	t323 = t307 * t329;
	t346 = t319 * qJD(3) - t284 * t300 - t299 * t323;
	t318 = -t289 * t300 + t299 * t334;
	t345 = t318 * qJD(3) - t284 * t299 + t300 * t323;
	t308 = cos(qJ(3));
	t298 = t308 * pkin(3) + pkin(2);
	t344 = t339 * t299 + t342 * t300 + t298;
	t340 = r_i_i_C(1) + qJ(4) + pkin(9);
	t291 = -t325 + t330;
	t338 = t291 * t299;
	t337 = t302 * t306;
	t336 = t302 * t307;
	t335 = t302 * t308;
	t327 = qJD(2) * t309;
	t324 = t303 * t330;
	t322 = t310 * t329;
	t321 = t302 * t327;
	t317 = t291 * t300 + t299 * t336;
	t316 = t303 * t299 + t300 * t337;
	t315 = t303 * t332 + t331;
	t288 = -t324 + t333;
	t285 = t316 * qJD(3) + t299 * t321;
	t283 = t315 * qJD(1) + t289 * qJD(2);
	t282 = t289 * qJD(1) + t315 * qJD(2);
	t281 = -qJD(1) * t324 - t310 * t327 + t320 * t333;
	t276 = t299 * t322 - qJD(3) * t338 + (qJD(3) * t336 - t282) * t300;
	t275 = t317 * qJD(3) - t282 * t299 - t300 * t322;
	t1 = [-t319 * qJD(5) - t284 * t298 - t288 * qJD(4) - t340 * t283 + t342 * t346 + t339 * t345 + (-t310 * pkin(1) - pkin(8) * t336) * qJD(1) + (-t305 * t323 + (t289 * t305 + t308 * t334) * qJD(3)) * pkin(3), t291 * qJD(4) + t344 * t281 - t340 * t282 + t315 * t347, t317 * qJD(5) + t339 * t276 - t342 * t275 + (t308 * t322 + t282 * t305 + (-t291 * t308 - t305 * t336) * qJD(3)) * pkin(3), -t281, t275, 0; -(t300 * t336 - t338) * qJD(5) - t282 * t298 + t315 * qJD(4) - t340 * t281 + t342 * t276 + t339 * t275 + (-t307 * pkin(1) + pkin(8) * t334) * qJD(1) + (t305 * t322 + (-t291 * t305 + t307 * t335) * qJD(3)) * pkin(3), t289 * qJD(4) - t283 * t344 + t340 * t284 + t347 * t288, -t318 * qJD(5) - t339 * t346 + t342 * t345 + (t308 * t323 - t284 * t305 + (-t289 * t308 + t305 * t334) * qJD(3)) * pkin(3), t283, -t345, 0; 0, ((-qJD(2) * t344 + qJD(4)) * t306 + (t340 * qJD(2) - t347) * t309) * t302, t316 * qJD(5) + t339 * (t300 * t321 + (-t299 * t337 + t300 * t303) * qJD(3)) - t342 * t285 + (-t305 * t321 + (-t303 * t305 - t306 * t335) * qJD(3)) * pkin(3), t302 * t328, t285, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:58
	% EndTime: 2019-10-10 11:25:59
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (836->115), mult. (1772->186), div. (0->0), fcn. (1768->12), ass. (0->68)
	t389 = sin(qJ(2));
	t390 = sin(qJ(1));
	t393 = cos(qJ(2));
	t394 = cos(qJ(1));
	t433 = cos(pkin(6));
	t410 = t394 * t433;
	t369 = t389 * t410 + t390 * t393;
	t384 = qJ(3) + pkin(11);
	t382 = sin(t384);
	t383 = cos(t384);
	t385 = sin(pkin(6));
	t423 = t385 * t394;
	t443 = t369 * t383 - t382 * t423;
	t408 = t393 * t410;
	t422 = t390 * t389;
	t368 = -t408 + t422;
	t387 = sin(qJ(6));
	t391 = cos(qJ(6));
	t439 = t369 * t382 + t383 * t423;
	t442 = t368 * t391 + t387 * t439;
	t441 = t368 * t387 - t391 * t439;
	t388 = sin(qJ(3));
	t407 = t391 * r_i_i_C(1) - t387 * r_i_i_C(2);
	t397 = t407 * qJD(6) + qJD(5);
	t406 = -t387 * r_i_i_C(1) - t391 * r_i_i_C(2);
	t404 = qJ(5) - t406;
	t417 = r_i_i_C(3) + pkin(10) + pkin(4);
	t395 = (t388 * pkin(3) + t417 * t382 - t404 * t383) * qJD(3) - t397 * t382;
	t405 = qJD(2) * t433 + qJD(1);
	t411 = t390 * t433;
	t409 = t389 * t411;
	t419 = qJD(2) * t389;
	t421 = t394 * t393;
	t357 = -qJD(1) * t409 - t390 * t419 + t405 * t421;
	t420 = qJD(1) * t385;
	t415 = t390 * t420;
	t438 = qJD(3) * t439 - t357 * t383 - t382 * t415;
	t392 = cos(qJ(3));
	t381 = t392 * pkin(3) + pkin(2);
	t437 = t404 * t382 + t417 * t383 + t381;
	t434 = -pkin(5) - qJ(4) - pkin(9);
	t371 = -t409 + t421;
	t428 = t371 * t382;
	t427 = t385 * t389;
	t426 = t385 * t390;
	t425 = t385 * t392;
	t424 = t385 * t393;
	t418 = qJD(2) * t393;
	t414 = t394 * t420;
	t413 = t385 * t418;
	t412 = t385 * t419;
	t402 = t371 * t383 + t382 * t426;
	t401 = t407 - t434;
	t370 = t394 * t389 + t393 * t411;
	t366 = t382 * t427 - t433 * t383;
	t399 = t433 * t382 + t383 * t427;
	t398 = t406 * qJD(6) + qJD(4);
	t350 = t443 * qJD(3) + t357 * t382 - t383 * t415;
	t363 = -t383 * t426 + t428;
	t358 = t399 * qJD(3) + t382 * t413;
	t356 = t370 * qJD(1) + t369 * qJD(2);
	t355 = t369 * qJD(1) + t370 * qJD(2);
	t354 = -qJD(1) * t408 - t394 * t418 + t405 * t422;
	t349 = t382 * t414 - qJD(3) * t428 + (qJD(3) * t426 - t355) * t383;
	t348 = t402 * qJD(3) - t355 * t382 - t383 * t414;
	t347 = t348 * t387 - t354 * t391 + (t363 * t391 - t370 * t387) * qJD(6);
	t346 = t348 * t391 + t354 * t387 + (-t363 * t387 - t370 * t391) * qJD(6);
	t1 = [-t368 * qJD(4) - t439 * qJD(5) - t357 * t381 - t404 * t350 - t401 * t356 + (t441 * r_i_i_C(1) + t442 * r_i_i_C(2)) * qJD(6) + (-t394 * pkin(1) - pkin(8) * t426) * qJD(1) + t417 * t438 + (-t388 * t415 + (t369 * t388 + t392 * t423) * qJD(3)) * pkin(3), t437 * t354 - t401 * t355 + t395 * t370 + t398 * t371, t397 * t402 + t404 * t349 - t417 * t348 + (t392 * t414 + t355 * t388 + (-t371 * t392 - t388 * t426) * qJD(3)) * pkin(3), -t354, t348, t346 * r_i_i_C(1) - t347 * r_i_i_C(2); t347 * r_i_i_C(1) + t346 * r_i_i_C(2) + t348 * qJ(5) + t370 * qJD(4) + t363 * qJD(5) - t355 * t381 + t434 * t354 + (-pkin(1) * t390 + pkin(8) * t423) * qJD(1) + t417 * t349 + (t388 * t414 + (-t371 * t388 + t390 * t425) * qJD(3)) * pkin(3), -t356 * t437 + t401 * t357 + t395 * t368 + t398 * t369, t397 * t443 - t404 * t438 - t417 * t350 + (t392 * t415 - t357 * t388 + (-t369 * t392 + t388 * t423) * qJD(3)) * pkin(3), t356, t350, (t350 * t391 - t356 * t387) * r_i_i_C(1) + (-t350 * t387 - t356 * t391) * r_i_i_C(2) + (-t442 * r_i_i_C(1) + t441 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t437 + t398) * t389 + (t401 * qJD(2) - t395) * t393) * t385, t397 * t399 + t404 * (-t366 * qJD(3) + t383 * t413) - t417 * t358 + (-t388 * t413 + (-t433 * t388 - t389 * t425) * qJD(3)) * pkin(3), t412, t358, (t358 * t391 - t387 * t412) * r_i_i_C(1) + (-t358 * t387 - t391 * t412) * r_i_i_C(2) + ((-t366 * t387 + t391 * t424) * r_i_i_C(1) + (-t366 * t391 - t387 * t424) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end