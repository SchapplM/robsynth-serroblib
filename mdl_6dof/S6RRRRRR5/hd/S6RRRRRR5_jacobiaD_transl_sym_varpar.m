% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR5
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
% Datum: 2019-10-10 13:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:27
	% EndTime: 2019-10-10 13:22:27
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
	% StartTime: 2019-10-10 13:22:27
	% EndTime: 2019-10-10 13:22:28
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
	% StartTime: 2019-10-10 13:22:27
	% EndTime: 2019-10-10 13:22:28
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (649->85), mult. (874->128), div. (0->0), fcn. (805->12), ass. (0->66)
	t296 = cos(pkin(6));
	t298 = sin(qJ(2));
	t299 = sin(qJ(1));
	t329 = t299 * t298;
	t317 = t296 * t329;
	t301 = cos(qJ(2));
	t302 = cos(qJ(1));
	t326 = t302 * t301;
	t271 = -qJD(1) * t317 - qJD(2) * t329 + (qJD(2) * t296 + qJD(1)) * t326;
	t292 = qJD(3) + qJD(4);
	t287 = qJD(5) + t292;
	t295 = sin(pkin(6));
	t330 = t295 * t302;
	t340 = t287 * t330 - t271;
	t294 = qJ(3) + qJ(4);
	t288 = sin(t294);
	t339 = pkin(4) * t288;
	t297 = sin(qJ(3));
	t280 = pkin(3) * t297 + t339;
	t338 = pkin(8) + t280;
	t337 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9);
	t336 = pkin(3) * qJD(3);
	t327 = t302 * t298;
	t328 = t299 * t301;
	t276 = t296 * t327 + t328;
	t335 = t276 * t287;
	t290 = qJ(5) + t294;
	t285 = sin(t290);
	t334 = t285 * t287;
	t289 = cos(t294);
	t333 = t289 * t292;
	t332 = t295 * t298;
	t331 = t295 * t299;
	t286 = cos(t290);
	t307 = t317 - t326;
	t321 = qJD(1) * t302;
	t315 = t295 * t321;
	t305 = t287 * t307 + t315;
	t308 = t296 * t328 + t327;
	t269 = qJD(1) * t276 + qJD(2) * t308;
	t311 = t287 * t331 - t269;
	t264 = -t285 * t311 + t286 * t305;
	t265 = t285 * t305 + t286 * t311;
	t325 = t264 * r_i_i_C(1) - t265 * r_i_i_C(2);
	t322 = qJD(1) * t299;
	t316 = t295 * t322;
	t306 = t316 - t335;
	t313 = t340 * t286;
	t324 = (t285 * t340 + t306 * t286) * r_i_i_C(1) + (-t285 * t306 + t313) * r_i_i_C(2);
	t320 = qJD(2) * t301;
	t314 = t295 * t320;
	t304 = -t287 * t296 - t314;
	t319 = t287 * t332;
	t323 = (t285 * t304 - t286 * t319) * r_i_i_C(1) + (t285 * t319 + t286 * t304) * r_i_i_C(2);
	t300 = cos(qJ(3));
	t281 = t300 * pkin(3) + pkin(4) * t289;
	t312 = -r_i_i_C(1) * t285 - r_i_i_C(2) * t286;
	t279 = pkin(2) + t281;
	t310 = r_i_i_C(1) * t286 - r_i_i_C(2) * t285 + t279;
	t309 = t296 * t326 - t329;
	t273 = -t292 * t339 - t297 * t336;
	t303 = t287 * t312 + t273;
	t274 = pkin(4) * t333 + t300 * t336;
	t270 = qJD(1) * t308 + qJD(2) * t276;
	t268 = -qJD(1) * t309 + qJD(2) * t307;
	t1 = [(t276 * t334 + t313) * r_i_i_C(1) + (t271 * t285 + t286 * t335) * r_i_i_C(2) - t271 * t279 - t276 * t273 - pkin(1) * t321 - t337 * t270 + ((-r_i_i_C(2) * t334 + t274) * t302 + (t312 - t338) * t322) * t295, t268 * t310 - t269 * t337 - t303 * t308, t269 * t280 + t307 * t274 + (t273 * t299 + t281 * t321) * t295 + t325, ((t292 * t307 + t315) * t289 + (-t292 * t331 + t269) * t288) * pkin(4) + t325, t325, 0; t274 * t331 + t265 * r_i_i_C(1) + t264 * r_i_i_C(2) - t269 * t279 - t307 * t273 - t337 * t268 + (-pkin(1) * t299 + t330 * t338) * qJD(1), -t270 * t310 + t271 * t337 + t303 * t309, -t271 * t280 - t276 * t274 + (-t273 * t302 + t281 * t322) * t295 + t324, ((-t276 * t292 + t316) * t289 + (t292 * t330 - t271) * t288) * pkin(4) + t324, t324, 0; 0, (t303 * t301 + (-t298 * t310 + t301 * t337) * qJD(2)) * t295, t296 * t273 + (-t274 * t298 - t280 * t320) * t295 + t323, (-t332 * t333 + (-t292 * t296 - t314) * t288) * pkin(4) + t323, t323, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:30
	% EndTime: 2019-10-10 13:22:31
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (1718->141), mult. (2202->217), div. (0->0), fcn. (2167->14), ass. (0->94)
	t510 = pkin(12) + r_i_i_C(3);
	t451 = sin(qJ(6));
	t455 = cos(qJ(6));
	t468 = t455 * r_i_i_C(1) - t451 * r_i_i_C(2);
	t466 = pkin(5) + t468;
	t482 = qJD(6) * t455;
	t483 = qJD(6) * t451;
	t518 = -r_i_i_C(1) * t483 - t482 * r_i_i_C(2);
	t450 = cos(pkin(6));
	t453 = sin(qJ(2));
	t458 = cos(qJ(1));
	t490 = t458 * t453;
	t454 = sin(qJ(1));
	t457 = cos(qJ(2));
	t491 = t454 * t457;
	t426 = t450 * t490 + t491;
	t448 = qJ(3) + qJ(4);
	t444 = qJ(5) + t448;
	t439 = sin(t444);
	t440 = cos(t444);
	t449 = sin(pkin(6));
	t493 = t449 * t458;
	t414 = -t426 * t440 + t439 * t493;
	t489 = t458 * t457;
	t476 = t450 * t489;
	t492 = t454 * t453;
	t425 = -t476 + t492;
	t517 = -t414 * t451 - t425 * t455;
	t516 = t414 * t455 - t425 * t451;
	t446 = qJD(3) + qJD(4);
	t452 = sin(qJ(3));
	t504 = pkin(3) * qJD(3);
	t442 = sin(t448);
	t509 = pkin(4) * t442;
	t423 = -t446 * t509 - t452 * t504;
	t441 = qJD(5) + t446;
	t515 = -(t466 * t439 - t510 * t440) * t441 + t423;
	t513 = t451 * r_i_i_C(1) + t455 * r_i_i_C(2);
	t443 = cos(t448);
	t456 = cos(qJ(3));
	t432 = t456 * pkin(3) + pkin(4) * t443;
	t430 = pkin(2) + t432;
	t511 = t510 * t439 + t466 * t440 + t430;
	t431 = t452 * pkin(3) + t509;
	t505 = pkin(8) + t431;
	t427 = t450 * t491 + t490;
	t410 = t427 * qJD(1) + t426 * qJD(2);
	t503 = t410 * t451;
	t502 = t410 * t455;
	t478 = t450 * t492;
	t428 = -t478 + t489;
	t499 = t428 * t440;
	t498 = t439 * t441;
	t497 = t443 * t446;
	t496 = t449 * t453;
	t495 = t449 * t454;
	t494 = t449 * t457;
	t488 = qJD(1) * t454;
	t487 = qJD(1) * t458;
	t486 = qJD(2) * t453;
	t485 = qJD(2) * t457;
	t484 = qJD(6) * t440;
	t480 = t439 * t496;
	t479 = t440 * t496;
	t477 = t440 * t493;
	t475 = t449 * t488;
	t474 = t449 * t487;
	t473 = t449 * t486;
	t472 = t449 * t485;
	t470 = qJD(2) * t450 + qJD(1);
	t411 = -qJD(1) * t478 - t454 * t486 + t470 * t489;
	t471 = -t411 * t440 + t441 * t477;
	t409 = t426 * qJD(1) + t427 * qJD(2);
	t467 = t441 * t495 - t409;
	t465 = -t426 * t441 + t475;
	t464 = t441 * t450 + t472;
	t399 = t465 * t440 + (t441 * t493 - t411) * t439;
	t400 = -t426 * t498 + t439 * t475 - t471;
	t462 = t518 * (-t426 * t439 - t477) + t510 * t400 + t466 * t399;
	t397 = t467 * t439 - t440 * t474 + t441 * t499;
	t398 = -t428 * t498 + t439 * t474 + t467 * t440;
	t461 = t518 * (-t428 * t439 + t440 * t495) + t510 * t398 - t466 * t397;
	t407 = t464 * t440 - t441 * t480;
	t460 = t518 * (t450 * t440 - t480) + t510 * t407 + t466 * (-t464 * t439 - t441 * t479);
	t459 = t513 * t484 - t515;
	t447 = -pkin(11) - pkin(10) - pkin(9);
	t424 = pkin(4) * t497 + t456 * t504;
	t420 = t450 * t439 + t479;
	t416 = t439 * t495 + t499;
	t408 = -qJD(1) * t476 - t458 * t485 + t470 * t492;
	t402 = -t465 * t439 + t471;
	t390 = t398 * t455 - t408 * t451 + (-t416 * t451 + t427 * t455) * qJD(6);
	t389 = -t398 * t451 - t408 * t455 + (-t416 * t455 - t427 * t451) * qJD(6);
	t1 = [(t402 * t455 - t503) * r_i_i_C(1) + (-t402 * t451 - t502) * r_i_i_C(2) + t402 * pkin(5) - t411 * t430 - t426 * t423 + t410 * t447 + t424 * t493 + t510 * t399 + (t517 * r_i_i_C(1) - t516 * r_i_i_C(2)) * qJD(6) + (-t458 * pkin(1) - t505 * t495) * qJD(1), (-t409 * t451 + t428 * t482) * r_i_i_C(1) + (-t409 * t455 - t428 * t483) * r_i_i_C(2) + t409 * t447 + t511 * t408 + t459 * t427, t409 * t431 - t428 * t424 + (t423 * t454 + t432 * t487) * t449 + t461, ((-t428 * t446 + t474) * t443 + (-t446 * t495 + t409) * t442) * pkin(4) + t461, t461, t389 * r_i_i_C(1) - t390 * r_i_i_C(2); t424 * t495 + t398 * pkin(5) + t390 * r_i_i_C(1) + t389 * r_i_i_C(2) + t408 * t447 - t409 * t430 + t428 * t423 + t510 * t397 + (-pkin(1) * t454 + t505 * t493) * qJD(1), (t411 * t451 + t426 * t482) * r_i_i_C(1) + (t411 * t455 - t426 * t483) * r_i_i_C(2) - t411 * t447 - t511 * t410 + t459 * t425, -t411 * t431 - t426 * t424 + (-t423 * t458 + t432 * t488) * t449 + t462, ((-t426 * t446 + t475) * t443 + (t446 * t493 - t411) * t442) * pkin(4) + t462, t462, (-t400 * t451 + t502) * r_i_i_C(1) + (-t400 * t455 - t503) * r_i_i_C(2) + (t516 * r_i_i_C(1) + t517 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t511 + t468 * qJD(6)) * t453 + (-qJD(2) * t447 + t513 * (qJD(2) - t484) + t515) * t457) * t449, t450 * t423 + (-t424 * t453 - t431 * t485) * t449 + t460, (-t496 * t497 + (-t446 * t450 - t472) * t442) * pkin(4) + t460, t460, (-t407 * t451 + t455 * t473) * r_i_i_C(1) + (-t407 * t455 - t451 * t473) * r_i_i_C(2) + ((-t420 * t455 + t451 * t494) * r_i_i_C(1) + (t420 * t451 + t455 * t494) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end