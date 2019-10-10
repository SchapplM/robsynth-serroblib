% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t175 = sin(pkin(6));
	t194 = t175 * (pkin(8) + r_i_i_C(1));
	t193 = pkin(2) - r_i_i_C(2);
	t191 = r_i_i_C(3) + qJ(3);
	t177 = sin(qJ(2));
	t178 = sin(qJ(1));
	t190 = t178 * t177;
	t179 = cos(qJ(2));
	t189 = t178 * t179;
	t180 = cos(qJ(1));
	t188 = t180 * t177;
	t187 = t180 * t179;
	t186 = qJD(2) * t177;
	t176 = cos(pkin(6));
	t185 = t176 * t190;
	t184 = t176 * t187;
	t183 = qJD(2) * t176 + qJD(1);
	t182 = t176 * t189 + t188;
	t181 = t176 * t188 + t189;
	t170 = -qJD(1) * t185 - t178 * t186 + t183 * t187;
	t169 = t182 * qJD(1) + t181 * qJD(2);
	t168 = t181 * qJD(1) + t182 * qJD(2);
	t167 = -qJD(1) * t184 - qJD(2) * t187 + t183 * t190;
	t1 = [-(-t184 + t190) * qJD(3) - t193 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1), -(t185 - t187) * qJD(3) - t191 * t168 + t193 * t167, -t167, 0, 0, 0; t182 * qJD(3) - t193 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), t181 * qJD(3) - t193 * t169 + t191 * t170, t169, 0, 0, 0; 0, (t177 * qJD(3) + (-t193 * t177 + t191 * t179) * qJD(2)) * t175, t175 * t186, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:17
	% EndTime: 2019-10-10 10:28:17
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
	t275 = pkin(3) + pkin(8);
	t243 = cos(pkin(6));
	t245 = sin(qJ(2));
	t249 = cos(qJ(1));
	t267 = t249 * t245;
	t246 = sin(qJ(1));
	t248 = cos(qJ(2));
	t268 = t246 * t248;
	t235 = t243 * t268 + t267;
	t251 = t243 * t267 + t268;
	t231 = t235 * qJD(1) + t251 * qJD(2);
	t244 = sin(qJ(4));
	t274 = t231 * t244;
	t247 = cos(qJ(4));
	t273 = t231 * t247;
	t242 = sin(pkin(6));
	t272 = t242 * t246;
	t271 = t242 * t248;
	t270 = t242 * t249;
	t269 = t246 * t245;
	t266 = t249 * t248;
	t265 = qJD(1) * t246;
	t264 = qJD(1) * t249;
	t263 = qJD(2) * t245;
	t259 = t243 * t266;
	t233 = -t259 + t269;
	t262 = qJD(4) * t233;
	t261 = -r_i_i_C(3) - pkin(9) - pkin(2);
	t260 = t243 * t269;
	t258 = t242 * t265;
	t257 = t242 * t264;
	t256 = t242 * t263;
	t255 = qJD(2) * t243 + qJD(1);
	t254 = r_i_i_C(1) * t247 - r_i_i_C(2) * t244;
	t253 = -t244 * r_i_i_C(1) - t247 * r_i_i_C(2);
	t252 = qJ(3) - t253;
	t250 = t254 * qJD(4) + qJD(3);
	t232 = -qJD(1) * t260 - t246 * t263 + t255 * t266;
	t230 = t251 * qJD(1) + t235 * qJD(2);
	t229 = -qJD(1) * t259 - qJD(2) * t266 + t255 * t269;
	t228 = t247 * t257 - t229 * t244 + (t235 * t247 - t244 * t272) * qJD(4);
	t227 = -t244 * t257 - t229 * t247 + (-t235 * t244 - t247 * t272) * qJD(4);
	t1 = [(-t247 * t262 - t274) * r_i_i_C(1) + (t244 * t262 - t273) * r_i_i_C(2) - t231 * qJ(3) - t233 * qJD(3) - pkin(1) * t264 + t261 * t232 + (t253 * t249 * qJD(4) + (-t254 - t275) * t265) * t242, t250 * (-t260 + t266) - t252 * t230 - t261 * t229, -t229, t227 * r_i_i_C(1) - t228 * r_i_i_C(2), 0, 0; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t229 * qJ(3) + t235 * qJD(3) + t261 * t230 + (-pkin(1) * t246 + t275 * t270) * qJD(1), t261 * t231 + t252 * t232 + t250 * t251, t231, (-t244 * t258 + t273) * r_i_i_C(1) + (-t247 * t258 - t274) * r_i_i_C(2) + ((-t233 * t244 + t247 * t270) * r_i_i_C(1) + (-t233 * t247 - t244 * t270) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t250 * t245 + (t261 * t245 + t252 * t248) * qJD(2)) * t242, t256, t254 * t256 + ((-t243 * t247 + t244 * t271) * r_i_i_C(1) + (t243 * t244 + t247 * t271) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:18
	% EndTime: 2019-10-10 10:28:18
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (312->65), mult. (921->101), div. (0->0), fcn. (892->8), ass. (0->48)
	t298 = sin(qJ(4));
	t301 = cos(qJ(4));
	t330 = r_i_i_C(3) + qJ(5);
	t331 = r_i_i_C(2) - pkin(4);
	t333 = t331 * t298 + t330 * t301 - qJ(3);
	t332 = pkin(3) + pkin(8);
	t296 = sin(pkin(6));
	t300 = sin(qJ(1));
	t329 = t296 * t300;
	t302 = cos(qJ(2));
	t328 = t296 * t302;
	t303 = cos(qJ(1));
	t327 = t296 * t303;
	t299 = sin(qJ(2));
	t326 = t300 * t299;
	t325 = t300 * t302;
	t324 = t303 * t299;
	t323 = t303 * t302;
	t322 = qJD(1) * t296;
	t321 = qJD(2) * t299;
	t320 = qJD(2) * t302;
	t319 = pkin(2) + pkin(9) + r_i_i_C(1);
	t318 = t298 * t329;
	t297 = cos(pkin(6));
	t317 = t297 * t326;
	t316 = t301 * t327;
	t315 = t297 * t323;
	t314 = t300 * t322;
	t313 = t303 * t322;
	t312 = t296 * t321;
	t311 = qJD(2) * t297 + qJD(1);
	t285 = -t315 + t326;
	t310 = t285 * t301 + t298 * t327;
	t287 = t297 * t325 + t324;
	t309 = t287 * t298 + t301 * t329;
	t308 = -t297 * t301 + t298 * t328;
	t307 = t297 * t324 + t325;
	t281 = t287 * qJD(1) + t307 * qJD(2);
	t272 = -t281 * t301 - qJD(4) * t316 + (qJD(4) * t285 + t314) * t298;
	t305 = -t301 * qJD(5) + qJD(3) + (t330 * t298 - t331 * t301) * qJD(4);
	t304 = t310 * qJD(4) + t281 * t298 + t301 * t314;
	t284 = -t308 * qJD(4) - t301 * t312;
	t282 = -qJD(1) * t317 - t300 * t321 + t311 * t323;
	t280 = t307 * qJD(1) + t287 * qJD(2);
	t279 = -qJD(1) * t315 - t303 * t320 + t311 * t326;
	t277 = -t279 * t298 - qJD(4) * t318 + (qJD(4) * t287 + t313) * t301;
	t276 = t309 * qJD(4) + t279 * t301 + t298 * t313;
	t1 = [t310 * qJD(5) - t281 * qJ(3) - t285 * qJD(3) + t331 * t304 - t330 * t272 - t319 * t282 + (-t303 * pkin(1) - t332 * t329) * qJD(1), t319 * t279 + t333 * t280 + t305 * (-t317 + t323), -t279, t309 * qJD(5) + t331 * t276 + t330 * t277, t276, 0; -(t287 * t301 - t318) * qJD(5) - t279 * qJ(3) + t287 * qJD(3) - t331 * t277 + t330 * t276 - t319 * t280 + (-t300 * pkin(1) + t332 * t327) * qJD(1), -t319 * t281 - t282 * t333 + t305 * t307, t281, -(-t285 * t298 + t316) * qJD(5) + t330 * t304 + t331 * t272, t272, 0; 0, (-t333 * t320 + (-t319 * qJD(2) + t305) * t299) * t296, t312, -t308 * qJD(5) + t331 * t284 - t330 * (-t298 * t312 + (t297 * t298 + t301 * t328) * qJD(4)), t284, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:18
	% EndTime: 2019-10-10 10:28:19
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (579->103), mult. (1713->165), div. (0->0), fcn. (1716->10), ass. (0->67)
	t380 = cos(pkin(6));
	t387 = cos(qJ(2));
	t388 = cos(qJ(1));
	t417 = t388 * t387;
	t407 = t380 * t417;
	t383 = sin(qJ(2));
	t384 = sin(qJ(1));
	t420 = t384 * t383;
	t365 = -t407 + t420;
	t382 = sin(qJ(4));
	t386 = cos(qJ(4));
	t379 = sin(pkin(6));
	t421 = t379 * t388;
	t357 = t365 * t386 + t382 * t421;
	t418 = t388 * t383;
	t419 = t384 * t387;
	t366 = t380 * t418 + t419;
	t381 = sin(qJ(6));
	t385 = cos(qJ(6));
	t431 = t357 * t381 - t366 * t385;
	t430 = t357 * t385 + t366 * t381;
	t397 = t381 * r_i_i_C(1) + t385 * r_i_i_C(2);
	t394 = qJD(6) * t397;
	t398 = t385 * r_i_i_C(1) - t381 * r_i_i_C(2);
	t391 = qJD(6) * t398 + qJD(5);
	t396 = qJ(5) + t397;
	t410 = -r_i_i_C(3) - pkin(10) - pkin(4);
	t429 = t382 * t410 + t386 * t396 - qJ(3);
	t428 = pkin(3) + pkin(8);
	t367 = t380 * t419 + t418;
	t351 = qJD(1) * t367 + qJD(2) * t366;
	t427 = t351 * t382;
	t424 = t379 * t383;
	t423 = t379 * t384;
	t422 = t379 * t387;
	t416 = qJD(1) * t384;
	t415 = qJD(2) * t383;
	t414 = qJD(2) * t387;
	t413 = qJD(4) * t382;
	t412 = qJD(4) * t386;
	t411 = -pkin(2) - pkin(5) - pkin(9);
	t409 = t380 * t420;
	t408 = t386 * t421;
	t406 = t379 * t416;
	t405 = qJD(1) * t421;
	t404 = t386 * t416;
	t403 = t379 * t414;
	t402 = t379 * t413;
	t401 = t379 * t415;
	t400 = -qJD(4) * t408 - t351 * t386;
	t399 = qJD(2) * t380 + qJD(1);
	t395 = t367 * t382 + t386 * t423;
	t363 = t380 * t382 + t386 * t422;
	t393 = t398 - t411;
	t389 = qJD(3) - t391 * t386 + (t382 * t396 - t386 * t410) * qJD(4);
	t368 = -t409 + t417;
	t355 = -t367 * t386 + t382 * t423;
	t354 = t380 * t412 - t386 * t401 - t387 * t402;
	t352 = -qJD(1) * t409 - t384 * t415 + t399 * t417;
	t350 = qJD(1) * t366 + qJD(2) * t367;
	t349 = -qJD(1) * t407 - t388 * t414 + t399 * t420;
	t345 = -t349 * t382 - t384 * t402 + (qJD(4) * t367 + t405) * t386;
	t344 = qJD(4) * t395 + t349 * t386 + t382 * t405;
	t340 = (qJD(4) * t365 + t406) * t382 + t400;
	t339 = t344 * t381 - t350 * t385 + (t355 * t385 - t368 * t381) * qJD(6);
	t338 = t344 * t385 + t350 * t381 + (-t355 * t381 - t368 * t385) * qJD(6);
	t1 = [-t351 * qJ(3) - t365 * qJD(3) + t357 * qJD(5) - t396 * (t365 * t413 + t382 * t406 + t400) + (t430 * r_i_i_C(1) - t431 * r_i_i_C(2)) * qJD(6) - t393 * t352 + t410 * (t365 * t412 + t427 + (t388 * t413 + t404) * t379) + (-t388 * pkin(1) - t423 * t428) * qJD(1), t393 * t349 + t429 * t350 + t367 * t394 + t389 * t368, -t349, t344 * t410 + t345 * t396 + t391 * t395, t344, t338 * r_i_i_C(1) - t339 * r_i_i_C(2); t339 * r_i_i_C(1) + t338 * r_i_i_C(2) - t349 * qJ(3) + t344 * qJ(5) + t367 * qJD(3) + t355 * qJD(5) + t411 * t350 - t410 * t345 + (-pkin(1) * t384 + t421 * t428) * qJD(1), -t351 * t393 - t352 * t429 + t365 * t394 + t366 * t389, t351, -t391 * (-t365 * t382 + t408) + t396 * (qJD(4) * t357 + t404 * t379 + t427) + t410 * t340, t340, (t340 * t385 - t352 * t381) * r_i_i_C(1) + (-t340 * t381 - t352 * t385) * r_i_i_C(2) + (t431 * r_i_i_C(1) + t430 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t429 - t394) * t387 + (-qJD(2) * t393 + t389) * t383) * t379, t401, t391 * (t380 * t386 - t382 * t422) + t410 * t354 - t396 * (qJD(4) * t363 - t382 * t401), t354, (t354 * t385 - t381 * t403) * r_i_i_C(1) + (-t354 * t381 - t385 * t403) * r_i_i_C(2) + ((-t363 * t381 - t385 * t424) * r_i_i_C(1) + (-t363 * t385 + t381 * t424) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end