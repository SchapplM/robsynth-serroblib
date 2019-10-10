% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->14), mult. (83->33), div. (0->0), fcn. (72->8), ass. (0->16)
	t89 = sin(pkin(11));
	t92 = cos(pkin(11));
	t95 = sin(qJ(2));
	t96 = cos(qJ(2));
	t98 = t89 * t96 + t92 * t95;
	t88 = t98 * qJD(2);
	t94 = cos(pkin(6));
	t100 = t94 * t95;
	t99 = pkin(2) * qJD(2);
	t97 = t89 * t95 - t92 * t96;
	t87 = t97 * qJD(2);
	t93 = cos(pkin(10));
	t90 = sin(pkin(10));
	t86 = t94 * t88;
	t85 = t94 * t87;
	t1 = [0, (t90 * t86 + t93 * t87) * r_i_i_C(1) + (-t90 * t85 + t93 * t88) * r_i_i_C(2) + (t90 * t100 - t93 * t96) * t99, 0, 0, 0, 0; 0, (-t93 * t86 + t90 * t87) * r_i_i_C(1) + (t93 * t85 + t90 * t88) * r_i_i_C(2) + (-t93 * t100 - t90 * t96) * t99, 0, 0, 0, 0; 0, (-t95 * pkin(2) - t98 * r_i_i_C(1) + t97 * r_i_i_C(2)) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:41
	% EndTime: 2019-10-09 21:33:41
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (115->39), mult. (386->84), div. (0->0), fcn. (382->10), ass. (0->37)
	t241 = sin(pkin(11));
	t244 = cos(pkin(11));
	t250 = cos(qJ(2));
	t261 = qJD(2) * t250;
	t248 = sin(qJ(2));
	t262 = qJD(2) * t248;
	t268 = t241 * t262 - t244 * t261;
	t267 = -pkin(8) - r_i_i_C(3);
	t266 = pkin(2) * qJD(2);
	t243 = sin(pkin(6));
	t247 = sin(qJ(4));
	t265 = t243 * t247;
	t249 = cos(qJ(4));
	t264 = t243 * t249;
	t246 = cos(pkin(6));
	t263 = t246 * t248;
	t258 = r_i_i_C(1) * t247 + r_i_i_C(2) * t249;
	t227 = t268 * t246;
	t234 = -t241 * t261 - t244 * t262;
	t242 = sin(pkin(10));
	t245 = cos(pkin(10));
	t257 = t227 * t245 - t234 * t242;
	t256 = t227 * t242 + t234 * t245;
	t255 = t241 * t250 + t248 * t244;
	t254 = t248 * t241 - t244 * t250;
	t253 = r_i_i_C(1) * t249 - r_i_i_C(2) * t247 + pkin(3);
	t252 = qJD(4) * t258;
	t251 = qJD(2) * t255;
	t233 = t254 * qJD(2);
	t232 = t255 * t246;
	t231 = t254 * t246;
	t230 = t255 * t243;
	t228 = t246 * t251;
	t225 = t268 * t243;
	t224 = -t232 * t242 - t245 * t254;
	t222 = t232 * t245 - t242 * t254;
	t1 = [0, -t267 * t256 - (t231 * t242 - t245 * t255) * t252 + (t242 * t263 - t245 * t250) * t266 + t253 * (t228 * t242 + t233 * t245), 0, -t258 * t256 + ((-t224 * t249 - t242 * t265) * r_i_i_C(1) + (t224 * t247 - t242 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t257 - (-t231 * t245 - t242 * t255) * t252 + (-t242 * t250 - t245 * t263) * t266 + t253 * (-t228 * t245 + t233 * t242), 0, t258 * t257 + ((-t222 * t249 + t245 * t265) * r_i_i_C(1) + (t222 * t247 + t245 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t225 + (-pkin(2) * t262 - t251 * t253 + t252 * t254) * t243, 0, t258 * t225 + ((-t230 * t249 - t246 * t247) * r_i_i_C(1) + (t230 * t247 - t246 * t249) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:42
	% EndTime: 2019-10-09 21:33:42
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (244->48), mult. (780->94), div. (0->0), fcn. (810->10), ass. (0->44)
	t292 = sin(pkin(11));
	t295 = cos(pkin(11));
	t301 = cos(qJ(2));
	t313 = qJD(2) * t301;
	t299 = sin(qJ(2));
	t314 = qJD(2) * t299;
	t322 = t292 * t314 - t295 * t313;
	t321 = pkin(4) - r_i_i_C(2);
	t320 = -pkin(8) - r_i_i_C(1);
	t319 = r_i_i_C(3) + qJ(5);
	t318 = pkin(2) * qJD(2);
	t294 = sin(pkin(6));
	t298 = sin(qJ(4));
	t317 = t294 * t298;
	t300 = cos(qJ(4));
	t316 = t294 * t300;
	t297 = cos(pkin(6));
	t315 = t297 * t299;
	t278 = t322 * t297;
	t285 = -t292 * t313 - t295 * t314;
	t293 = sin(pkin(10));
	t296 = cos(pkin(10));
	t310 = t296 * t278 - t293 * t285;
	t270 = t293 * t278 + t296 * t285;
	t308 = t301 * t292 + t299 * t295;
	t281 = t308 * t294;
	t309 = t281 * t300 + t297 * t298;
	t307 = t299 * t292 - t301 * t295;
	t283 = t308 * t297;
	t273 = t296 * t283 - t293 * t307;
	t306 = -t273 * t300 + t296 * t317;
	t275 = -t293 * t283 - t296 * t307;
	t305 = t275 * t300 + t293 * t317;
	t304 = qJD(2) * t308;
	t303 = t319 * t298 + t321 * t300 + pkin(3);
	t302 = qJD(5) * t298 + (-t321 * t298 + t319 * t300) * qJD(4);
	t284 = t307 * qJD(2);
	t282 = t307 * t297;
	t279 = t297 * t304;
	t276 = t322 * t294;
	t264 = t309 * qJD(4) - t276 * t298;
	t262 = t305 * qJD(4) + t270 * t298;
	t260 = -t306 * qJD(4) - t298 * t310;
	t1 = [0, -t320 * t270 + (t293 * t315 - t296 * t301) * t318 + t302 * (t293 * t282 - t296 * t308) + t303 * (t293 * t279 + t296 * t284), 0, t305 * qJD(5) + t319 * (t270 * t300 + (-t275 * t298 + t293 * t316) * qJD(4)) - t321 * t262, t262, 0; 0, t320 * t310 + (-t293 * t301 - t296 * t315) * t318 + t302 * (-t296 * t282 - t293 * t308) + t303 * (-t296 * t279 + t293 * t284), 0, -t306 * qJD(5) + t319 * (-t310 * t300 + (-t273 * t298 - t296 * t316) * qJD(4)) - t321 * t260, t260, 0; 0, t320 * t276 + (-pkin(2) * t314 - t302 * t307 - t303 * t304) * t294, 0, t309 * qJD(5) + t319 * (-t276 * t300 + (-t281 * t298 + t297 * t300) * qJD(4)) - t321 * t264, t264, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:43
	% EndTime: 2019-10-09 21:33:43
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (525->74), mult. (1644->141), div. (0->0), fcn. (1774->12), ass. (0->58)
	t390 = sin(qJ(4));
	t393 = cos(qJ(4));
	t389 = sin(qJ(6));
	t392 = cos(qJ(6));
	t410 = t392 * r_i_i_C(1) - t389 * r_i_i_C(2);
	t398 = t410 * qJD(6) + qJD(5);
	t409 = -t389 * r_i_i_C(1) - t392 * r_i_i_C(2);
	t406 = qJ(5) - t409;
	t416 = pkin(4) + pkin(9) + r_i_i_C(3);
	t430 = (t416 * t390 - t406 * t393) * qJD(4) - t398 * t390;
	t388 = cos(pkin(6));
	t384 = sin(pkin(11));
	t391 = sin(qJ(2));
	t423 = cos(pkin(11));
	t425 = cos(qJ(2));
	t400 = t425 * t384 + t391 * t423;
	t372 = t400 * t388;
	t411 = t425 * t423;
	t417 = qJD(2) * t391;
	t428 = -qJD(2) * t411 + t384 * t417;
	t399 = -t391 * t384 + t411;
	t385 = sin(pkin(10));
	t387 = cos(pkin(10));
	t358 = t387 * t372 + t385 * t399;
	t386 = sin(pkin(6));
	t421 = t386 * t390;
	t427 = t358 * t393 - t387 * t421;
	t395 = t406 * t390 + t416 * t393 + pkin(3);
	t424 = pkin(2) * qJD(2);
	t420 = t386 * t393;
	t419 = t388 * t391;
	t369 = t428 * t388;
	t374 = t400 * qJD(2);
	t408 = t387 * t369 + t385 * t374;
	t354 = t385 * t369 - t387 * t374;
	t371 = t400 * t386;
	t407 = t371 * t393 + t388 * t390;
	t362 = t371 * t390 - t388 * t393;
	t361 = -t385 * t372 + t387 * t399;
	t404 = -pkin(5) - pkin(8) - t410;
	t346 = t358 * t390 + t387 * t420;
	t403 = -t361 * t390 + t385 * t420;
	t402 = t361 * t393 + t385 * t421;
	t401 = qJD(6) * t409;
	t397 = t399 * t388;
	t396 = qJD(2) * t372;
	t373 = t399 * qJD(2);
	t370 = t399 * t386;
	t368 = qJD(2) * t371;
	t367 = t428 * t386;
	t360 = -t385 * t397 - t387 * t400;
	t357 = -t385 * t400 + t387 * t397;
	t353 = -t387 * t373 + t385 * t396;
	t350 = -t385 * t373 - t387 * t396;
	t344 = t407 * qJD(4) - t367 * t390;
	t342 = t402 * qJD(4) + t354 * t390;
	t340 = t427 * qJD(4) - t390 * t408;
	t1 = [0, t361 * t401 - t404 * t354 + (t385 * t419 - t425 * t387) * t424 + t395 * t353 - t430 * t360, 0, t398 * t402 + t406 * (t403 * qJD(4) + t354 * t393) - t416 * t342, t342, (t342 * t392 + t353 * t389) * r_i_i_C(1) + (-t342 * t389 + t353 * t392) * r_i_i_C(2) + ((t360 * t392 + t389 * t403) * r_i_i_C(1) + (-t360 * t389 + t392 * t403) * r_i_i_C(2)) * qJD(6); 0, t358 * t401 + t404 * t408 + (-t425 * t385 - t387 * t419) * t424 + t395 * t350 - t430 * t357, 0, t398 * t427 + t406 * (-t346 * qJD(4) - t393 * t408) - t416 * t340, t340, (t340 * t392 + t350 * t389) * r_i_i_C(1) + (-t340 * t389 + t350 * t392) * r_i_i_C(2) + ((-t346 * t389 + t357 * t392) * r_i_i_C(1) + (-t346 * t392 - t357 * t389) * r_i_i_C(2)) * qJD(6); 0, -t386 * pkin(2) * t417 + t404 * t367 - t395 * t368 - t430 * t370 + t371 * t401, 0, t398 * t407 + t406 * (-t362 * qJD(4) - t367 * t393) - t416 * t344, t344, (t344 * t392 - t368 * t389) * r_i_i_C(1) + (-t344 * t389 - t368 * t392) * r_i_i_C(2) + ((-t362 * t389 + t370 * t392) * r_i_i_C(1) + (-t362 * t392 - t370 * t389) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end