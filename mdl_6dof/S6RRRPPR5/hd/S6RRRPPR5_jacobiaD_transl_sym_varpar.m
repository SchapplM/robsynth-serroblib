% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:04
	% EndTime: 2019-10-10 11:24:04
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
	% StartTime: 2019-10-10 11:24:04
	% EndTime: 2019-10-10 11:24:05
	% DurationCPUTime: 0.37s
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
	t239 = qJD(1) * t267 + qJD(2) * t242;
	t238 = qJD(1) * t242 + qJD(2) * t267;
	t237 = -qJD(1) * t273 - t263 * t277 + t269 * t284;
	t236 = t252 * t271 - t238 * t253 + (-t244 * t252 + t253 * t287) * qJD(3);
	t235 = t253 * t271 + t238 * t252 + (-t244 * t253 - t252 * t287) * qJD(3);
	t1 = [(-t240 * t253 + t252 * t276 + t245) * r_i_i_C(1) + (t240 * t252 + t253 * t276) * r_i_i_C(2) - t240 * t251 + t276 * t292 - t241 * qJD(4) - pkin(1) * t279 - t289 * t239 + ((t290 - t291) * t275 + (-pkin(8) - t265) * t280) * t255, t244 * qJD(4) + t237 * t268 - t238 * t289 + t264 * t267, t235 * r_i_i_C(1) - t236 * r_i_i_C(2) + (t261 * t271 + t238 * t258 + (-t244 * t261 - t258 * t287) * qJD(3)) * pkin(3), -t237, 0, 0; t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t267 * qJD(4) - t238 * t251 - t289 * t237 + (-pkin(1) * t260 + pkin(8) * t285) * qJD(1) + (t258 * t271 + (-t244 * t258 + t260 * t286) * qJD(3)) * pkin(3), t242 * qJD(4) - t239 * t268 + t240 * t289 + t241 * t264, t245 * r_i_i_C(2) + (r_i_i_C(1) * t266 - t240 * r_i_i_C(2)) * t253 + ((-t240 + t270) * r_i_i_C(1) - t266 * r_i_i_C(2)) * t252 + (t261 * t272 - t240 * t258 + (-t242 * t261 + t258 * t285) * qJD(3)) * pkin(3), t239, 0, 0; 0, (qJD(4) * t259 - t262 * t264 + (-t259 * t268 + t262 * t289) * qJD(2)) * t255, -t265 * t255 * t277 + ((-t252 * t256 - t253 * t288) * r_i_i_C(1) + (t252 * t288 - t253 * t256) * r_i_i_C(2) + (-t256 * t258 - t259 * t286) * pkin(3)) * qJD(3), t255 * t278, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:05
	% EndTime: 2019-10-10 11:24:06
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (550->95), mult. (1188->149), div. (0->0), fcn. (1158->12), ass. (0->56)
	t345 = sin(qJ(1));
	t341 = cos(pkin(6));
	t358 = qJD(2) * t341 + qJD(1);
	t344 = sin(qJ(2));
	t374 = t345 * t344;
	t365 = t341 * t374;
	t369 = qJD(2) * t344;
	t347 = cos(qJ(2));
	t348 = cos(qJ(1));
	t371 = t348 * t347;
	t318 = -qJD(1) * t365 - t345 * t369 + t358 * t371;
	t339 = sin(pkin(6));
	t375 = t339 * t348;
	t386 = -qJD(3) * t375 + t318;
	t372 = t348 * t344;
	t373 = t345 * t347;
	t322 = t341 * t372 + t373;
	t370 = qJD(1) * t339;
	t363 = t345 * t370;
	t385 = -qJD(3) * t322 + t363;
	t337 = qJ(3) + pkin(11);
	t335 = sin(t337);
	t336 = cos(t337);
	t343 = sin(qJ(3));
	t338 = sin(pkin(12));
	t340 = cos(pkin(12));
	t357 = t340 * r_i_i_C(1) - t338 * r_i_i_C(2) + pkin(4);
	t379 = r_i_i_C(3) + qJ(5);
	t384 = (pkin(3) * t343 + t357 * t335 - t379 * t336) * qJD(3) - t335 * qJD(5);
	t383 = t335 * t385 + t386 * t336;
	t346 = cos(qJ(3));
	t334 = t346 * pkin(3) + pkin(2);
	t382 = t379 * t335 + t357 * t336 + t334;
	t378 = t339 * t344;
	t377 = t339 * t345;
	t376 = t339 * t346;
	t368 = qJD(2) * t347;
	t364 = t341 * t371;
	t362 = t348 * t370;
	t361 = t339 * t368;
	t342 = -qJ(4) - pkin(9);
	t356 = t338 * r_i_i_C(1) + t340 * r_i_i_C(2) - t342;
	t324 = -t365 + t371;
	t355 = -t324 * t335 + t336 * t377;
	t354 = t324 * t336 + t335 * t377;
	t353 = t341 * t335 + t336 * t378;
	t352 = t341 * t373 + t372;
	t311 = t386 * t335 - t385 * t336;
	t321 = -t364 + t374;
	t319 = t353 * qJD(3) + t335 * t361;
	t317 = t352 * qJD(1) + t322 * qJD(2);
	t316 = t322 * qJD(1) + t352 * qJD(2);
	t315 = -qJD(1) * t364 - t348 * t368 + t358 * t374;
	t310 = t355 * qJD(3) - t316 * t336 + t335 * t362;
	t309 = t354 * qJD(3) - t316 * t335 - t336 * t362;
	t1 = [(-t317 * t338 - t340 * t383) * r_i_i_C(1) + (-t317 * t340 + t338 * t383) * r_i_i_C(2) - t383 * pkin(4) - (t322 * t335 + t336 * t375) * qJD(5) - t318 * t334 + t317 * t342 - t321 * qJD(4) - t379 * t311 + (-t348 * pkin(1) - pkin(8) * t377) * qJD(1) + (-t343 * t363 + (t322 * t343 + t346 * t375) * qJD(3)) * pkin(3), t324 * qJD(4) + t382 * t315 - t356 * t316 + t352 * t384, t354 * qJD(5) + t379 * t310 - t357 * t309 + (t346 * t362 + t316 * t343 + (-t324 * t346 - t343 * t377) * qJD(3)) * pkin(3), -t315, t309, 0; (t310 * t340 - t315 * t338) * r_i_i_C(1) + (-t310 * t338 - t315 * t340) * r_i_i_C(2) + t310 * pkin(4) - t355 * qJD(5) - t316 * t334 + t315 * t342 + t352 * qJD(4) + t379 * t309 + (-t345 * pkin(1) + pkin(8) * t375) * qJD(1) + (t343 * t362 + (-t324 * t343 + t345 * t376) * qJD(3)) * pkin(3), t322 * qJD(4) - t317 * t382 + t356 * t318 + t384 * t321, -(-t322 * t336 + t335 * t375) * qJD(5) + t379 * t383 - t357 * t311 + (t346 * t363 - t318 * t343 + (-t322 * t346 + t343 * t375) * qJD(3)) * pkin(3), t317, t311, 0; 0, ((-qJD(2) * t382 + qJD(4)) * t344 + (t356 * qJD(2) - t384) * t347) * t339, t353 * qJD(5) + t379 * (t336 * t361 + (-t335 * t378 + t336 * t341) * qJD(3)) - t357 * t319 + (-t343 * t361 + (-t341 * t343 - t344 * t376) * qJD(3)) * pkin(3), t339 * t369, t319, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:06
	% EndTime: 2019-10-10 11:24:06
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (866->120), mult. (1643->191), div. (0->0), fcn. (1644->14), ass. (0->72)
	t402 = sin(qJ(1));
	t443 = cos(pkin(6));
	t415 = qJD(2) * t443 + qJD(1);
	t401 = sin(qJ(2));
	t422 = t402 * t443;
	t419 = t401 * t422;
	t432 = qJD(2) * t401;
	t404 = cos(qJ(2));
	t405 = cos(qJ(1));
	t434 = t405 * t404;
	t364 = -qJD(1) * t419 - t402 * t432 + t415 * t434;
	t397 = sin(pkin(6));
	t436 = t397 * t405;
	t452 = -qJD(3) * t436 + t364;
	t421 = t405 * t443;
	t375 = t401 * t421 + t402 * t404;
	t433 = qJD(1) * t397;
	t428 = t402 * t433;
	t451 = -qJD(3) * t375 + t428;
	t395 = qJ(3) + pkin(11);
	t391 = sin(t395);
	t393 = cos(t395);
	t368 = t375 * t393 - t391 * t436;
	t418 = t404 * t421;
	t435 = t402 * t401;
	t374 = -t418 + t435;
	t394 = pkin(12) + qJ(6);
	t390 = sin(t394);
	t392 = cos(t394);
	t450 = t368 * t390 - t374 * t392;
	t449 = t368 * t392 + t374 * t390;
	t400 = sin(qJ(3));
	t416 = t390 * r_i_i_C(1) + t392 * r_i_i_C(2);
	t412 = qJD(6) * t416;
	t388 = cos(pkin(12)) * pkin(5) + pkin(4);
	t417 = t392 * r_i_i_C(1) - t390 * r_i_i_C(2);
	t414 = t388 + t417;
	t444 = r_i_i_C(3) + pkin(10) + qJ(5);
	t406 = (t400 * pkin(3) + t414 * t391 - t444 * t393) * qJD(3) - t391 * qJD(5) + t393 * t412;
	t358 = t451 * t391 + t452 * t393;
	t403 = cos(qJ(3));
	t389 = t403 * pkin(3) + pkin(2);
	t447 = t444 * t391 + t414 * t393 + t389;
	t440 = t397 * t401;
	t439 = t397 * t402;
	t438 = t397 * t403;
	t437 = t397 * t404;
	t431 = qJD(2) * t404;
	t427 = t405 * t433;
	t426 = t397 * t431;
	t424 = t397 * t432;
	t423 = -sin(pkin(12)) * pkin(5) - qJ(4) - pkin(9);
	t413 = t375 * t391 + t393 * t436;
	t377 = -t419 + t434;
	t370 = -t377 * t391 + t393 * t439;
	t371 = t377 * t393 + t391 * t439;
	t376 = t405 * t401 + t404 * t422;
	t410 = -t391 * t440 + t443 * t393;
	t373 = t443 * t391 + t393 * t440;
	t409 = t416 - t423;
	t408 = t417 * qJD(6) + qJD(4);
	t357 = t452 * t391 - t451 * t393;
	t366 = t410 * qJD(3) + t393 * t426;
	t365 = t373 * qJD(3) + t391 * t426;
	t363 = t376 * qJD(1) + t375 * qJD(2);
	t362 = t375 * qJD(1) + t376 * qJD(2);
	t361 = -qJD(1) * t418 - t405 * t431 + t415 * t435;
	t356 = t370 * qJD(3) - t362 * t393 + t391 * t427;
	t355 = t371 * qJD(3) - t362 * t391 - t393 * t427;
	t354 = t356 * t392 - t361 * t390 + (-t371 * t390 + t376 * t392) * qJD(6);
	t353 = -t356 * t390 - t361 * t392 + (-t371 * t392 - t376 * t390) * qJD(6);
	t1 = [-t413 * qJD(5) - t364 * t389 - t374 * qJD(4) - t414 * t358 - t409 * t363 - t444 * t357 + (t450 * r_i_i_C(1) + t449 * r_i_i_C(2)) * qJD(6) + (-t405 * pkin(1) - pkin(8) * t439) * qJD(1) + (-t400 * t428 + (t375 * t400 + t403 * t436) * qJD(3)) * pkin(3), t447 * t361 - t409 * t362 + t406 * t376 + t408 * t377, t371 * qJD(5) + t444 * t356 - t370 * t412 - t414 * t355 + (t403 * t427 + t362 * t400 + (-t377 * t403 - t400 * t439) * qJD(3)) * pkin(3), -t361, t355, t353 * r_i_i_C(1) - t354 * r_i_i_C(2); t354 * r_i_i_C(1) + t353 * r_i_i_C(2) + t376 * qJD(4) - t370 * qJD(5) + t356 * t388 - t362 * t389 + t423 * t361 + t444 * t355 + (-pkin(1) * t402 + pkin(8) * t436) * qJD(1) + (t400 * t427 + (-t377 * t400 + t402 * t438) * qJD(3)) * pkin(3), -t363 * t447 + t409 * t364 + t406 * t374 + t408 * t375, t368 * qJD(5) + t444 * t358 + t413 * t412 - t414 * t357 + (t403 * t428 - t364 * t400 + (-t375 * t403 + t400 * t436) * qJD(3)) * pkin(3), t363, t357, (-t358 * t390 + t363 * t392) * r_i_i_C(1) + (-t358 * t392 - t363 * t390) * r_i_i_C(2) + (-t449 * r_i_i_C(1) + t450 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t447 + t408) * t401 + (t409 * qJD(2) - t406) * t404) * t397, t373 * qJD(5) + t444 * t366 - t410 * t412 - t414 * t365 + (-t400 * t426 + (-t443 * t400 - t401 * t438) * qJD(3)) * pkin(3), t424, t365, (-t366 * t390 + t392 * t424) * r_i_i_C(1) + (-t366 * t392 - t390 * t424) * r_i_i_C(2) + ((-t373 * t392 + t390 * t437) * r_i_i_C(1) + (t373 * t390 + t392 * t437) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end