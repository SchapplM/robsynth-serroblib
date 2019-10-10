% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
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
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:15
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
	% StartTime: 2019-10-10 12:42:16
	% EndTime: 2019-10-10 12:42:16
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
	% StartTime: 2019-10-10 12:42:16
	% EndTime: 2019-10-10 12:42:16
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
	% StartTime: 2019-10-10 12:42:16
	% EndTime: 2019-10-10 12:42:16
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (523->86), mult. (804->130), div. (0->0), fcn. (741->12), ass. (0->63)
	t301 = sin(qJ(1));
	t298 = cos(pkin(6));
	t314 = qJD(2) * t298 + qJD(1);
	t300 = sin(qJ(2));
	t330 = t301 * t300;
	t317 = t298 * t330;
	t321 = qJD(2) * t300;
	t303 = cos(qJ(2));
	t304 = cos(qJ(1));
	t327 = t304 * t303;
	t271 = -qJD(1) * t317 - t301 * t321 + t314 * t327;
	t295 = qJD(3) + qJD(4);
	t297 = sin(pkin(6));
	t331 = t297 * t304;
	t311 = t295 * t331 - t271;
	t296 = qJ(3) + qJ(4);
	t291 = sin(t296);
	t299 = sin(qJ(3));
	t335 = pkin(3) * qJD(3);
	t338 = pkin(4) * t295;
	t272 = -t291 * t338 - t299 * t335;
	t290 = pkin(12) + t296;
	t287 = sin(t290);
	t288 = cos(t290);
	t313 = r_i_i_C(1) * t287 + r_i_i_C(2) * t288;
	t339 = t313 * t295 - t272;
	t280 = t299 * pkin(3) + pkin(4) * t291;
	t337 = pkin(8) + t280;
	t336 = r_i_i_C(3) + qJ(5) + pkin(10) + pkin(9);
	t328 = t304 * t300;
	t329 = t301 * t303;
	t276 = t298 * t328 + t329;
	t334 = t276 * t295;
	t333 = t287 * t295;
	t332 = t297 * t301;
	t278 = -t317 + t327;
	t322 = qJD(1) * t304;
	t307 = -t278 * t295 + t297 * t322;
	t309 = t298 * t329 + t328;
	t269 = t276 * qJD(1) + t309 * qJD(2);
	t312 = t295 * t332 - t269;
	t264 = -t312 * t287 + t307 * t288;
	t265 = t307 * t287 + t312 * t288;
	t326 = t264 * r_i_i_C(1) - t265 * r_i_i_C(2);
	t323 = qJD(1) * t301;
	t308 = t297 * t323 - t334;
	t315 = t311 * t288;
	t325 = (t311 * t287 + t308 * t288) * r_i_i_C(1) + (-t308 * t287 + t315) * r_i_i_C(2);
	t320 = qJD(2) * t303;
	t306 = -t295 * t298 - t297 * t320;
	t319 = t295 * t297 * t300;
	t324 = (t306 * t287 - t288 * t319) * r_i_i_C(1) + (t287 * t319 + t306 * t288) * r_i_i_C(2);
	t292 = cos(t296);
	t302 = cos(qJ(3));
	t281 = t302 * pkin(3) + pkin(4) * t292;
	t316 = t298 * t327;
	t279 = pkin(2) + t281;
	t310 = t288 * r_i_i_C(1) - t287 * r_i_i_C(2) + t279;
	t275 = -t316 + t330;
	t273 = t292 * t338 + t302 * t335;
	t270 = t309 * qJD(1) + t276 * qJD(2);
	t268 = -qJD(1) * t316 - t304 * t320 + t314 * t330;
	t1 = [(t276 * t333 + t315) * r_i_i_C(1) + (t271 * t287 + t288 * t334) * r_i_i_C(2) - t271 * t279 - t276 * t272 - t275 * qJD(5) - pkin(1) * t322 - t336 * t270 + ((-r_i_i_C(2) * t333 + t273) * t304 + (-t313 - t337) * t323) * t297, t278 * qJD(5) + t310 * t268 - t336 * t269 + t309 * t339, t269 * t280 - t278 * t273 + (t272 * t301 + t281 * t322) * t297 + t326, (-t312 * t291 + t307 * t292) * pkin(4) + t326, -t268, 0; t273 * t332 + t265 * r_i_i_C(1) + t264 * r_i_i_C(2) + t309 * qJD(5) - t269 * t279 + t278 * t272 - t336 * t268 + (-pkin(1) * t301 + t337 * t331) * qJD(1), t276 * qJD(5) - t310 * t270 + t336 * t271 + t339 * t275, -t271 * t280 - t276 * t273 + (-t272 * t304 + t281 * t323) * t297 + t325, (t311 * t291 + t308 * t292) * pkin(4) + t325, t270, 0; 0, (qJD(5) * t300 - t339 * t303 + (-t310 * t300 + t336 * t303) * qJD(2)) * t297, t298 * t272 + (-t273 * t300 - t280 * t320) * t297 + t324, (t306 * t291 - t292 * t319) * pkin(4) + t324, t297 * t321, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:18
	% EndTime: 2019-10-10 12:42:19
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (1324->140), mult. (1913->216), div. (0->0), fcn. (1883->14), ass. (0->90)
	t507 = pkin(11) + r_i_i_C(3);
	t450 = sin(qJ(6));
	t454 = cos(qJ(6));
	t468 = r_i_i_C(1) * t454 - r_i_i_C(2) * t450;
	t465 = pkin(5) + t468;
	t453 = sin(qJ(1));
	t449 = cos(pkin(6));
	t470 = qJD(2) * t449 + qJD(1);
	t452 = sin(qJ(2));
	t490 = t452 * t453;
	t476 = t449 * t490;
	t484 = qJD(2) * t452;
	t456 = cos(qJ(2));
	t457 = cos(qJ(1));
	t487 = t456 * t457;
	t411 = -qJD(1) * t476 - t453 * t484 + t470 * t487;
	t446 = qJD(3) + qJD(4);
	t448 = sin(pkin(6));
	t491 = t448 * t457;
	t466 = t446 * t491 - t411;
	t480 = qJD(6) * t454;
	t481 = qJD(6) * t450;
	t516 = -r_i_i_C(1) * t481 - t480 * r_i_i_C(2);
	t488 = t453 * t456;
	t489 = t452 * t457;
	t426 = t449 * t489 + t488;
	t447 = qJ(3) + qJ(4);
	t441 = pkin(12) + t447;
	t438 = sin(t441);
	t439 = cos(t441);
	t414 = -t426 * t439 + t438 * t491;
	t475 = t449 * t487;
	t425 = -t475 + t490;
	t515 = -t414 * t450 - t425 * t454;
	t514 = t414 * t454 - t425 * t450;
	t442 = sin(t447);
	t451 = sin(qJ(3));
	t501 = pkin(3) * qJD(3);
	t506 = pkin(4) * t446;
	t421 = -t442 * t506 - t451 * t501;
	t513 = -(t465 * t438 - t507 * t439) * t446 + t421;
	t428 = -t476 + t487;
	t485 = qJD(1) * t457;
	t473 = t448 * t485;
	t511 = -t428 * t446 + t473;
	t510 = r_i_i_C(1) * t450 + r_i_i_C(2) * t454;
	t443 = cos(t447);
	t455 = cos(qJ(3));
	t432 = t455 * pkin(3) + pkin(4) * t443;
	t430 = pkin(2) + t432;
	t508 = t507 * t438 + t465 * t439 + t430;
	t431 = pkin(3) * t451 + pkin(4) * t442;
	t502 = pkin(8) + t431;
	t427 = t449 * t488 + t489;
	t410 = t427 * qJD(1) + t426 * qJD(2);
	t500 = t410 * t450;
	t499 = t410 * t454;
	t495 = t438 * t446;
	t494 = t448 * t452;
	t493 = t448 * t453;
	t492 = t448 * t456;
	t486 = qJD(1) * t453;
	t483 = qJD(2) * t456;
	t482 = qJD(6) * t439;
	t478 = t446 * t494;
	t474 = t448 * t486;
	t472 = t448 * t484;
	t471 = t466 * t439;
	t409 = t426 * qJD(1) + t427 * qJD(2);
	t467 = t446 * t493 - t409;
	t464 = -t426 * t446 + t474;
	t463 = t446 * t449 + t448 * t483;
	t399 = t466 * t438 + t464 * t439;
	t397 = t467 * t438 - t511 * t439;
	t398 = -t428 * t495 + t438 * t473 + t467 * t439;
	t461 = t516 * (-t428 * t438 + t439 * t493) + t507 * t398 - t465 * t397;
	t407 = -t438 * t478 + t463 * t439;
	t460 = t516 * (-t438 * t494 + t439 * t449) + t507 * t407 + t465 * (-t463 * t438 - t439 * t478);
	t400 = -t426 * t495 + t438 * t474 - t471;
	t459 = t516 * (-t426 * t438 - t439 * t491) + t507 * t400 + t465 * t399;
	t458 = t510 * t482 - t513;
	t445 = -qJ(5) - pkin(10) - pkin(9);
	t422 = t443 * t506 + t455 * t501;
	t420 = t438 * t449 + t439 * t494;
	t416 = t428 * t439 + t438 * t493;
	t408 = -qJD(1) * t475 - t457 * t483 + t470 * t490;
	t402 = -t464 * t438 + t471;
	t390 = t398 * t454 - t408 * t450 + (-t416 * t450 + t427 * t454) * qJD(6);
	t389 = -t398 * t450 - t408 * t454 + (-t416 * t454 - t427 * t450) * qJD(6);
	t1 = [(t402 * t454 - t500) * r_i_i_C(1) + (-t402 * t450 - t499) * r_i_i_C(2) + t402 * pkin(5) - t411 * t430 - t426 * t421 + t410 * t445 - t425 * qJD(5) + t422 * t491 + t507 * t399 + (t515 * r_i_i_C(1) - t514 * r_i_i_C(2)) * qJD(6) + (-t457 * pkin(1) - t502 * t493) * qJD(1), (-t409 * t450 + t428 * t480) * r_i_i_C(1) + (-t409 * t454 - t428 * t481) * r_i_i_C(2) + t409 * t445 + t428 * qJD(5) + t508 * t408 + t458 * t427, t409 * t431 - t428 * t422 + (t421 * t453 + t432 * t485) * t448 + t461, (-t467 * t442 + t511 * t443) * pkin(4) + t461, -t408, r_i_i_C(1) * t389 - r_i_i_C(2) * t390; t422 * t493 + t398 * pkin(5) + t390 * r_i_i_C(1) + t389 * r_i_i_C(2) + t427 * qJD(5) + t408 * t445 - t409 * t430 + t428 * t421 + t507 * t397 + (-pkin(1) * t453 + t502 * t491) * qJD(1), (t411 * t450 + t426 * t480) * r_i_i_C(1) + (t411 * t454 - t426 * t481) * r_i_i_C(2) - t411 * t445 + t426 * qJD(5) - t508 * t410 + t458 * t425, -t411 * t431 - t426 * t422 + (-t421 * t457 + t432 * t486) * t448 + t459, (t466 * t442 + t464 * t443) * pkin(4) + t459, t410, (-t400 * t450 + t499) * r_i_i_C(1) + (-t400 * t454 - t500) * r_i_i_C(2) + (t514 * r_i_i_C(1) + t515 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t508 + t468 * qJD(6) + qJD(5)) * t452 + (-qJD(2) * t445 + t510 * (qJD(2) - t482) + t513) * t456) * t448, t421 * t449 + (-t422 * t452 - t431 * t483) * t448 + t460, (-t463 * t442 - t443 * t478) * pkin(4) + t460, t472, (-t407 * t450 + t454 * t472) * r_i_i_C(1) + (-t407 * t454 - t450 * t472) * r_i_i_C(2) + ((-t420 * t454 + t450 * t492) * r_i_i_C(1) + (t420 * t450 + t454 * t492) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end