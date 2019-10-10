% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR7
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
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:13
	% EndTime: 2019-10-10 12:04:13
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
	% StartTime: 2019-10-10 12:04:13
	% EndTime: 2019-10-10 12:04:14
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
	% StartTime: 2019-10-10 12:04:13
	% EndTime: 2019-10-10 12:04:14
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (474->78), mult. (741->117), div. (0->0), fcn. (687->12), ass. (0->57)
	t299 = sin(qJ(1));
	t296 = cos(pkin(6));
	t311 = qJD(2) * t296 + qJD(1);
	t298 = sin(qJ(2));
	t327 = t299 * t298;
	t314 = t296 * t327;
	t318 = qJD(2) * t298;
	t301 = cos(qJ(2));
	t302 = cos(qJ(1));
	t324 = t302 * t301;
	t269 = -qJD(1) * t314 - t299 * t318 + t311 * t324;
	t293 = qJD(3) + qJD(5);
	t295 = sin(pkin(6));
	t328 = t295 * t302;
	t335 = t293 * t328 - t269;
	t294 = qJ(3) + pkin(12);
	t279 = pkin(4) * cos(t294) + cos(qJ(3)) * pkin(3);
	t278 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t294);
	t275 = t278 * qJD(3);
	t290 = qJ(5) + t294;
	t286 = sin(t290);
	t287 = cos(t290);
	t310 = r_i_i_C(1) * t286 + r_i_i_C(2) * t287;
	t334 = t310 * t293 + t275;
	t333 = pkin(8) + t278;
	t332 = r_i_i_C(3) + pkin(10) + qJ(4) + pkin(9);
	t325 = t302 * t298;
	t326 = t299 * t301;
	t272 = t296 * t325 + t326;
	t331 = t272 * t293;
	t330 = t286 * t293;
	t329 = t295 * t299;
	t274 = -t314 + t324;
	t319 = qJD(1) * t302;
	t305 = -t274 * t293 + t295 * t319;
	t307 = t296 * t326 + t325;
	t267 = t272 * qJD(1) + t307 * qJD(2);
	t309 = t293 * t329 - t267;
	t262 = -t309 * t286 + t305 * t287;
	t263 = t305 * t286 + t309 * t287;
	t323 = t262 * r_i_i_C(1) - t263 * r_i_i_C(2);
	t320 = qJD(1) * t299;
	t306 = t295 * t320 - t331;
	t312 = t335 * t287;
	t322 = (t335 * t286 + t306 * t287) * r_i_i_C(1) + (-t306 * t286 + t312) * r_i_i_C(2);
	t317 = qJD(2) * t301;
	t304 = -t293 * t296 - t295 * t317;
	t316 = t293 * t295 * t298;
	t321 = (t304 * t286 - t287 * t316) * r_i_i_C(1) + (t286 * t316 + t304 * t287) * r_i_i_C(2);
	t313 = t296 * t324;
	t277 = pkin(2) + t279;
	t308 = t287 * r_i_i_C(1) - t286 * r_i_i_C(2) + t277;
	t276 = t279 * qJD(3);
	t271 = -t313 + t327;
	t268 = t307 * qJD(1) + t272 * qJD(2);
	t266 = -qJD(1) * t313 - t302 * t317 + t311 * t327;
	t1 = [(t272 * t330 + t312) * r_i_i_C(1) + (t269 * t286 + t287 * t331) * r_i_i_C(2) - t269 * t277 + t272 * t275 - t271 * qJD(4) - pkin(1) * t319 - t332 * t268 + ((-r_i_i_C(2) * t330 + t276) * t302 + (-t310 - t333) * t320) * t295, t274 * qJD(4) + t308 * t266 - t332 * t267 + t307 * t334, t267 * t278 - t274 * t276 + (-t275 * t299 + t279 * t319) * t295 + t323, -t266, t323, 0; t276 * t329 + t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t307 * qJD(4) - t267 * t277 - t274 * t275 - t332 * t266 + (-pkin(1) * t299 + t333 * t328) * qJD(1), t272 * qJD(4) - t308 * t268 + t332 * t269 + t334 * t271, -t269 * t278 - t272 * t276 + (t275 * t302 + t279 * t320) * t295 + t322, t268, t322, 0; 0, (qJD(4) * t298 - t334 * t301 + (-t308 * t298 + t332 * t301) * qJD(2)) * t295, -t296 * t275 + (-t276 * t298 - t278 * t317) * t295 + t321, t295 * t318, t321, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:15
	% EndTime: 2019-10-10 12:04:16
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (1275->132), mult. (1850->203), div. (0->0), fcn. (1829->14), ass. (0->84)
	t502 = pkin(11) + r_i_i_C(3);
	t448 = sin(qJ(6));
	t452 = cos(qJ(6));
	t465 = t452 * r_i_i_C(1) - t448 * r_i_i_C(2);
	t463 = pkin(5) + t465;
	t451 = sin(qJ(1));
	t447 = cos(pkin(6));
	t467 = qJD(2) * t447 + qJD(1);
	t450 = sin(qJ(2));
	t487 = t450 * t451;
	t473 = t447 * t487;
	t481 = qJD(2) * t450;
	t454 = cos(qJ(2));
	t455 = cos(qJ(1));
	t484 = t454 * t455;
	t409 = -qJD(1) * t473 - t451 * t481 + t467 * t484;
	t444 = qJD(3) + qJD(5);
	t446 = sin(pkin(6));
	t488 = t446 * t455;
	t511 = t444 * t488 - t409;
	t477 = qJD(6) * t452;
	t478 = qJD(6) * t448;
	t510 = -r_i_i_C(1) * t478 - t477 * r_i_i_C(2);
	t485 = t451 * t454;
	t486 = t450 * t455;
	t422 = t447 * t486 + t485;
	t445 = qJ(3) + pkin(12);
	t441 = qJ(5) + t445;
	t437 = sin(t441);
	t438 = cos(t441);
	t412 = -t422 * t438 + t437 * t488;
	t472 = t447 * t484;
	t421 = -t472 + t487;
	t509 = -t412 * t448 - t421 * t452;
	t508 = t412 * t452 - t421 * t448;
	t429 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t445);
	t425 = t429 * qJD(3);
	t507 = -(t463 * t437 - t502 * t438) * t444 - t425;
	t430 = pkin(4) * cos(t445) + cos(qJ(3)) * pkin(3);
	t505 = t448 * r_i_i_C(1) + t452 * r_i_i_C(2);
	t428 = pkin(2) + t430;
	t503 = t502 * t437 + t463 * t438 + t428;
	t498 = pkin(8) + t429;
	t423 = t447 * t485 + t486;
	t408 = t423 * qJD(1) + t422 * qJD(2);
	t497 = t408 * t448;
	t496 = t408 * t452;
	t424 = -t473 + t484;
	t493 = t424 * t438;
	t492 = t437 * t444;
	t491 = t446 * t450;
	t490 = t446 * t451;
	t489 = t446 * t454;
	t483 = qJD(1) * t451;
	t482 = qJD(1) * t455;
	t480 = qJD(2) * t454;
	t479 = qJD(6) * t438;
	t475 = t444 * t491;
	t471 = t446 * t483;
	t470 = t446 * t482;
	t469 = t446 * t481;
	t468 = t511 * t438;
	t407 = t422 * qJD(1) + t423 * qJD(2);
	t464 = t444 * t490 - t407;
	t462 = -t422 * t444 + t471;
	t461 = t444 * t447 + t446 * t480;
	t397 = t511 * t437 + t462 * t438;
	t405 = -t437 * t475 + t461 * t438;
	t459 = t510 * (-t437 * t491 + t438 * t447) + t502 * t405 + t463 * (-t461 * t437 - t438 * t475);
	t398 = -t422 * t492 + t437 * t471 - t468;
	t458 = t510 * (-t422 * t437 - t438 * t488) + t502 * t398 + t463 * t397;
	t395 = t464 * t437 - t438 * t470 + t444 * t493;
	t396 = -t424 * t492 + t437 * t470 + t464 * t438;
	t457 = t510 * (-t424 * t437 + t438 * t490) + t502 * t396 - t463 * t395;
	t456 = t505 * t479 - t507;
	t443 = -pkin(10) - qJ(4) - pkin(9);
	t426 = t430 * qJD(3);
	t418 = t437 * t447 + t438 * t491;
	t414 = t437 * t490 + t493;
	t406 = -qJD(1) * t472 - t455 * t480 + t467 * t487;
	t400 = -t462 * t437 + t468;
	t388 = t396 * t452 - t406 * t448 + (-t414 * t448 + t423 * t452) * qJD(6);
	t387 = -t396 * t448 - t406 * t452 + (-t414 * t452 - t423 * t448) * qJD(6);
	t1 = [(t400 * t452 - t497) * r_i_i_C(1) + (-t400 * t448 - t496) * r_i_i_C(2) + t400 * pkin(5) - t409 * t428 + t422 * t425 + t408 * t443 - t421 * qJD(4) + t426 * t488 + t502 * t397 + (t509 * r_i_i_C(1) - t508 * r_i_i_C(2)) * qJD(6) + (-t455 * pkin(1) - t498 * t490) * qJD(1), (-t407 * t448 + t424 * t477) * r_i_i_C(1) + (-t407 * t452 - t424 * t478) * r_i_i_C(2) + t407 * t443 + t424 * qJD(4) + t503 * t406 + t456 * t423, t407 * t429 - t424 * t426 + (-t425 * t451 + t430 * t482) * t446 + t457, -t406, t457, r_i_i_C(1) * t387 - r_i_i_C(2) * t388; t426 * t490 + t396 * pkin(5) + t388 * r_i_i_C(1) + t387 * r_i_i_C(2) + t423 * qJD(4) + t406 * t443 - t407 * t428 - t424 * t425 + t502 * t395 + (-pkin(1) * t451 + t498 * t488) * qJD(1), (t409 * t448 + t422 * t477) * r_i_i_C(1) + (t409 * t452 - t422 * t478) * r_i_i_C(2) - t409 * t443 + t422 * qJD(4) - t503 * t408 + t456 * t421, -t409 * t429 - t422 * t426 + (t425 * t455 + t430 * t483) * t446 + t458, t408, t458, (-t398 * t448 + t496) * r_i_i_C(1) + (-t398 * t452 - t497) * r_i_i_C(2) + (t508 * r_i_i_C(1) + t509 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t503 + t465 * qJD(6) + qJD(4)) * t450 + (-qJD(2) * t443 + t505 * (qJD(2) - t479) + t507) * t454) * t446, -t425 * t447 + (-t426 * t450 - t429 * t480) * t446 + t459, t469, t459, (-t405 * t448 + t452 * t469) * r_i_i_C(1) + (-t405 * t452 - t448 * t469) * r_i_i_C(2) + ((-t418 * t452 + t448 * t489) * r_i_i_C(1) + (t418 * t448 + t452 * t489) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end