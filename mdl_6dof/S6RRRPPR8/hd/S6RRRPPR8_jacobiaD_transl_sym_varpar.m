% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:34
	% EndTime: 2019-10-10 11:29:34
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
	% StartTime: 2019-10-10 11:29:35
	% EndTime: 2019-10-10 11:29:35
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (278->57), mult. (828->96), div. (0->0), fcn. (802->8), ass. (0->43)
	t297 = cos(pkin(6));
	t299 = sin(qJ(2));
	t300 = sin(qJ(1));
	t324 = t300 * t299;
	t318 = t297 * t324;
	t302 = cos(qJ(2));
	t303 = cos(qJ(1));
	t321 = t303 * t302;
	t283 = -qJD(1) * t318 - qJD(2) * t324 + (qJD(2) * t297 + qJD(1)) * t321;
	t296 = sin(pkin(6));
	t325 = t296 * t303;
	t334 = -qJD(3) * t325 + t283;
	t322 = t303 * t299;
	t323 = t300 * t302;
	t287 = t297 * t322 + t323;
	t320 = qJD(1) * t296;
	t333 = -qJD(3) * t287 + t300 * t320;
	t298 = sin(qJ(3));
	t301 = cos(qJ(3));
	t332 = t333 * t298 + t334 * t301;
	t328 = r_i_i_C(3) + qJ(4);
	t330 = pkin(3) + r_i_i_C(1);
	t331 = t328 * t298 + t330 * t301 + pkin(2);
	t329 = pkin(9) + r_i_i_C(2);
	t327 = t296 * t300;
	t326 = t296 * t301;
	t316 = t303 * t320;
	t315 = qJD(2) * t296 * t302;
	t307 = t318 - t321;
	t312 = t298 * t307 + t300 * t326;
	t311 = t298 * t327 - t301 * t307;
	t310 = t297 * t298 + t299 * t326;
	t309 = t297 * t321 - t324;
	t308 = t297 * t323 + t322;
	t276 = t334 * t298 - t333 * t301;
	t304 = qJD(4) * t298 + (-t330 * t298 + t328 * t301) * qJD(3);
	t284 = t310 * qJD(3) + t298 * t315;
	t282 = t308 * qJD(1) + t287 * qJD(2);
	t281 = t287 * qJD(1) + t308 * qJD(2);
	t280 = -t309 * qJD(1) + t307 * qJD(2);
	t275 = t312 * qJD(3) - t281 * t301 + t298 * t316;
	t274 = t311 * qJD(3) - t281 * t298 - t301 * t316;
	t1 = [-(t287 * t298 + t301 * t325) * qJD(4) - t283 * pkin(2) - t329 * t282 - t330 * t332 - t328 * t276 + (-t303 * pkin(1) - pkin(8) * t327) * qJD(1), t331 * t280 - t329 * t281 - t304 * t308, t311 * qJD(4) - t330 * t274 + t328 * t275, t274, 0, 0; -t312 * qJD(4) - t281 * pkin(2) - t329 * t280 + t330 * t275 + t328 * t274 + (-t300 * pkin(1) + pkin(8) * t325) * qJD(1), -t282 * t331 + t329 * t283 + t304 * t309, -(-t287 * t301 + t298 * t325) * qJD(4) + t328 * t332 - t330 * t276, t276, 0, 0; 0, (t304 * t302 + (-t299 * t331 + t329 * t302) * qJD(2)) * t296, t310 * qJD(4) + t328 * (t301 * t315 + (-t296 * t298 * t299 + t297 * t301) * qJD(3)) - t330 * t284, t284, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:34
	% EndTime: 2019-10-10 11:29:34
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (360->66), mult. (1060->106), div. (0->0), fcn. (1028->8), ass. (0->46)
	t266 = sin(qJ(3));
	t269 = cos(qJ(3));
	t285 = r_i_i_C(2) - pkin(4) - pkin(3);
	t298 = r_i_i_C(1) + qJ(4);
	t272 = (t285 * t266 + t298 * t269) * qJD(3) + qJD(4) * t266;
	t265 = cos(pkin(6));
	t268 = sin(qJ(1));
	t267 = sin(qJ(2));
	t292 = t268 * t267;
	t282 = t265 * t292;
	t288 = qJD(2) * t267;
	t270 = cos(qJ(2));
	t271 = cos(qJ(1));
	t289 = t271 * t270;
	t251 = -qJD(1) * t282 - t268 * t288 + (qJD(2) * t265 + qJD(1)) * t289;
	t290 = t271 * t267;
	t291 = t268 * t270;
	t256 = t265 * t290 + t291;
	t264 = sin(pkin(6));
	t294 = t264 * t271;
	t278 = t256 * t266 + t269 * t294;
	t296 = t264 * t268;
	t283 = t266 * t296;
	t301 = -qJD(1) * t283 + t278 * qJD(3) - t251 * t269;
	t299 = t298 * t266 - t285 * t269 + pkin(2);
	t275 = t282 - t289;
	t297 = t275 * t266;
	t295 = t264 * t269;
	t293 = t266 * t271;
	t287 = qJD(3) * t269;
	t284 = pkin(9) - r_i_i_C(3) - qJ(5);
	t281 = t264 * t293;
	t280 = qJD(1) * t295;
	t279 = qJD(2) * t264 * t270;
	t277 = -t269 * t275 + t283;
	t276 = t265 * t266 + t267 * t295;
	t255 = t265 * t289 - t292;
	t257 = t265 * t291 + t290;
	t244 = -qJD(3) * t281 + t251 * t266 + t256 * t287 - t268 * t280;
	t252 = t276 * qJD(3) + t266 * t279;
	t250 = t257 * qJD(1) + t256 * qJD(2);
	t249 = t256 * qJD(1) + t257 * qJD(2);
	t248 = -t255 * qJD(1) + t275 * qJD(2);
	t243 = -t249 * t269 + qJD(3) * t297 + (qJD(1) * t293 + t268 * t287) * t264;
	t242 = t277 * qJD(3) - t249 * t266 - t271 * t280;
	t1 = [-t255 * qJD(5) - t278 * qJD(4) - t251 * pkin(2) - t298 * t244 + (-t271 * pkin(1) - pkin(8) * t296) * qJD(1) - t284 * t250 - t285 * t301, t275 * qJD(5) + t299 * t248 - t284 * t249 - t272 * t257, t277 * qJD(4) + t285 * t242 + t298 * t243, t242, t248, 0; -t257 * qJD(5) - (t268 * t295 + t297) * qJD(4) - t249 * pkin(2) + t298 * t242 + (-t268 * pkin(1) + pkin(8) * t294) * qJD(1) - t284 * t248 - t285 * t243, -t256 * qJD(5) - t250 * t299 + t284 * t251 + t272 * t255, -(-t256 * t269 + t281) * qJD(4) - t298 * t301 + t285 * t244, t244, -t250, 0; 0, (-qJD(5) * t267 + t272 * t270 + (-t267 * t299 + t284 * t270) * qJD(2)) * t264, t276 * qJD(4) + t298 * (t269 * t279 + (-t264 * t266 * t267 + t265 * t269) * qJD(3)) + t285 * t252, t252, -t264 * t288, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:35
	% EndTime: 2019-10-10 11:29:36
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (659->97), mult. (1950->157), div. (0->0), fcn. (1950->10), ass. (0->61)
	t370 = sin(qJ(3));
	t374 = cos(qJ(3));
	t369 = sin(qJ(6));
	t373 = cos(qJ(6));
	t389 = -t373 * r_i_i_C(1) + t369 * r_i_i_C(2);
	t413 = pkin(5) + qJ(4);
	t385 = -t389 + t413;
	t398 = pkin(3) + pkin(4) + pkin(10) + r_i_i_C(3);
	t388 = t369 * r_i_i_C(1) + t373 * r_i_i_C(2);
	t416 = t388 * qJD(6) - qJD(4);
	t377 = t416 * t370 + (t398 * t370 - t385 * t374) * qJD(3);
	t371 = sin(qJ(2));
	t372 = sin(qJ(1));
	t375 = cos(qJ(2));
	t376 = cos(qJ(1));
	t411 = cos(pkin(6));
	t391 = t376 * t411;
	t356 = t372 * t371 - t375 * t391;
	t357 = t371 * t391 + t372 * t375;
	t368 = sin(pkin(6));
	t403 = t368 * t376;
	t418 = t357 * t370 + t374 * t403;
	t420 = t356 * t373 + t369 * t418;
	t419 = -t356 * t369 + t373 * t418;
	t392 = t372 * t411;
	t390 = t371 * t392;
	t400 = qJD(2) * t371;
	t401 = t376 * t375;
	t345 = -qJD(1) * t390 - t372 * t400 + (qJD(2) * t411 + qJD(1)) * t401;
	t406 = t368 * t372;
	t397 = t370 * t406;
	t417 = -qJD(1) * t397 + qJD(3) * t418 - t345 * t374;
	t414 = t385 * t370 + t398 * t374 + pkin(2);
	t412 = -pkin(9) + qJ(5);
	t381 = t390 - t401;
	t407 = t381 * t370;
	t405 = t368 * t374;
	t404 = t368 * t375;
	t402 = t370 * t376;
	t399 = qJD(3) * t374;
	t396 = t368 * t402;
	t395 = qJD(1) * t405;
	t394 = qJD(2) * t404;
	t393 = t368 * t400;
	t386 = -t374 * t381 + t397;
	t384 = t388 + t412;
	t358 = t376 * t371 + t375 * t392;
	t354 = t368 * t371 * t370 - t411 * t374;
	t382 = t411 * t370 + t371 * t405;
	t379 = t389 * qJD(6) - qJD(5);
	t338 = -qJD(3) * t396 + t345 * t370 + t357 * t399 - t372 * t395;
	t351 = -t372 * t405 - t407;
	t346 = t382 * qJD(3) + t370 * t394;
	t344 = t358 * qJD(1) + t357 * qJD(2);
	t343 = t357 * qJD(1) + t358 * qJD(2);
	t342 = t356 * qJD(1) + t381 * qJD(2);
	t337 = -t343 * t374 + qJD(3) * t407 + (qJD(1) * t402 + t372 * t399) * t368;
	t336 = t386 * qJD(3) - t343 * t370 - t376 * t395;
	t335 = t336 * t373 + t342 * t369 + (-t351 * t369 - t358 * t373) * qJD(6);
	t334 = -t336 * t369 + t342 * t373 + (-t351 * t373 + t358 * t369) * qJD(6);
	t1 = [-t345 * pkin(2) - t418 * qJD(4) + t356 * qJD(5) + t384 * t344 - t385 * t338 + (t420 * r_i_i_C(1) + t419 * r_i_i_C(2)) * qJD(6) + (-pkin(1) * t376 - pkin(8) * t406) * qJD(1) + t398 * t417, t414 * t342 + t384 * t343 + t377 * t358 - t379 * t381, -t398 * t336 + t385 * t337 - t386 * t416, t336, t342, r_i_i_C(1) * t334 - r_i_i_C(2) * t335; -t343 * pkin(2) + t335 * r_i_i_C(1) + t334 * r_i_i_C(2) + t351 * qJD(4) - t358 * qJD(5) + t412 * t342 + t413 * t336 + (-pkin(1) * t372 + pkin(8) * t403) * qJD(1) + t398 * t337, -t344 * t414 - t384 * t345 + t377 * t356 + t379 * t357, -t416 * (t357 * t374 - t396) - t385 * t417 - t398 * t338, t338, -t344, (-t338 * t369 - t344 * t373) * r_i_i_C(1) + (-t338 * t373 + t344 * t369) * r_i_i_C(2) + (-t419 * r_i_i_C(1) + t420 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t414 + t379) * t371 + (-t384 * qJD(2) - t377) * t375) * t368, -t416 * t382 + t385 * (-t354 * qJD(3) + t374 * t394) - t398 * t346, t346, -t393, (-t346 * t369 - t373 * t393) * r_i_i_C(1) + (-t346 * t373 + t369 * t393) * r_i_i_C(2) + ((-t354 * t373 - t369 * t404) * r_i_i_C(1) + (t354 * t369 - t373 * t404) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end