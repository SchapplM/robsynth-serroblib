% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
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
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:45
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(10));
	t182 = cos(pkin(10));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:45
	% EndTime: 2019-10-09 22:12:46
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (152->37), mult. (504->70), div. (0->0), fcn. (488->10), ass. (0->34)
	t256 = sin(pkin(10));
	t259 = cos(pkin(10));
	t262 = sin(qJ(2));
	t260 = cos(pkin(6));
	t264 = cos(qJ(2));
	t276 = t260 * t264;
	t285 = -t256 * t262 + t259 * t276;
	t277 = t260 * t262;
	t250 = t256 * t264 + t259 * t277;
	t263 = cos(qJ(3));
	t257 = sin(pkin(6));
	t261 = sin(qJ(3));
	t279 = t257 * t261;
	t284 = -t250 * t263 + t259 * t279;
	t255 = sin(pkin(11));
	t258 = cos(pkin(11));
	t272 = t258 * r_i_i_C(1) - t255 * r_i_i_C(2) + pkin(3);
	t282 = r_i_i_C(3) + qJ(4);
	t283 = t282 * t261 + t272 * t263 + pkin(2);
	t278 = t257 * t263;
	t273 = qJD(2) * t257 * t264;
	t271 = t255 * r_i_i_C(1) + t258 * r_i_i_C(2) + pkin(8);
	t268 = t256 * t277 - t259 * t264;
	t270 = t256 * t279 - t263 * t268;
	t269 = t256 * t276 + t259 * t262;
	t267 = t260 * t261 + t262 * t278;
	t266 = qJD(2) * t283;
	t265 = t261 * qJD(4) + (-t272 * t261 + t282 * t263) * qJD(3);
	t247 = t269 * qJD(2);
	t245 = t285 * qJD(2);
	t243 = t267 * qJD(3) + t261 * t273;
	t241 = t270 * qJD(3) - t247 * t261;
	t239 = -t284 * qJD(3) + t245 * t261;
	t1 = [0, -t271 * t247 - t265 * t269 + t268 * t266, t270 * qJD(4) + t282 * (-t247 * t263 + (t256 * t278 + t261 * t268) * qJD(3)) - t272 * t241, t241, 0, 0; 0, t271 * t245 - t250 * t266 + t265 * t285, -t284 * qJD(4) + t282 * (t245 * t263 + (-t250 * t261 - t259 * t278) * qJD(3)) - t272 * t239, t239, 0, 0; 0, (t265 * t264 + (-t283 * t262 + t271 * t264) * qJD(2)) * t257, t267 * qJD(4) + t282 * (t263 * t273 + (t260 * t263 - t262 * t279) * qJD(3)) - t272 * t243, t243, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:46
	% EndTime: 2019-10-09 22:12:46
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (241->65), mult. (793->118), div. (0->0), fcn. (787->10), ass. (0->53)
	t286 = sin(pkin(10));
	t289 = cos(pkin(10));
	t294 = cos(qJ(2));
	t290 = cos(pkin(6));
	t292 = sin(qJ(2));
	t318 = t290 * t292;
	t279 = t286 * t294 + t289 * t318;
	t293 = cos(qJ(3));
	t287 = sin(pkin(6));
	t291 = sin(qJ(3));
	t320 = t287 * t291;
	t325 = -t279 * t293 + t289 * t320;
	t324 = pkin(4) + r_i_i_C(1);
	t323 = r_i_i_C(2) + qJ(4);
	t322 = r_i_i_C(3) + qJ(5);
	t319 = t287 * t293;
	t317 = t290 * t294;
	t316 = t292 * t293;
	t315 = qJD(2) * t292;
	t314 = qJD(2) * t294;
	t313 = qJD(3) * t291;
	t285 = sin(pkin(11));
	t312 = qJD(5) * t285;
	t288 = cos(pkin(11));
	t311 = t288 * qJD(5);
	t309 = t286 * t315;
	t308 = t287 * t314;
	t307 = t289 * t314;
	t306 = t294 * t313;
	t305 = -t279 * t291 - t289 * t319;
	t281 = -t286 * t318 + t289 * t294;
	t304 = -t281 * t291 + t286 * t319;
	t303 = t281 * t293 + t286 * t320;
	t302 = t286 * t317 + t289 * t292;
	t301 = t287 * t316 + t290 * t291;
	t300 = t290 * t293 - t292 * t320;
	t275 = t279 * qJD(2);
	t278 = -t286 * t292 + t289 * t317;
	t299 = -t275 * t293 - t278 * t313;
	t277 = -t290 * t309 + t307;
	t298 = -t277 * t293 + t302 * t313;
	t297 = -pkin(3) * t293 - t323 * t291 - pkin(2);
	t296 = -t322 * t285 - t324 * t288 - pkin(3);
	t295 = t293 * t312 + t291 * qJD(4) + (-pkin(3) * t291 + t323 * t293) * qJD(3);
	t276 = t302 * qJD(2);
	t274 = -t290 * t307 + t309;
	t273 = t300 * qJD(3) + t293 * t308;
	t272 = t301 * qJD(3) + t291 * t308;
	t269 = t304 * qJD(3) - t276 * t293;
	t268 = t303 * qJD(3) - t276 * t291;
	t267 = t305 * qJD(3) - t274 * t293;
	t266 = -t325 * qJD(3) - t274 * t291;
	t1 = [0, -t281 * t311 - t276 * pkin(8) + t324 * (-t276 * t285 + t298 * t288) + t322 * (t276 * t288 + t298 * t285) - t295 * t302 + t297 * t277, t303 * qJD(4) + t296 * t268 + t323 * t269 + t304 * t312, t268, t269 * t285 - t277 * t288, 0; 0, -t279 * t311 - t274 * pkin(8) + t324 * (-t274 * t285 + t299 * t288) + t322 * (t274 * t288 + t299 * t285) + t295 * t278 + t297 * t275, -t325 * qJD(4) + t296 * t266 + t323 * t267 + t305 * t312, t266, t267 * t285 - t275 * t288, 0; 0, (t324 * (-t288 * t306 + (t285 * t294 - t288 * t316) * qJD(2)) - t322 * (t285 * t306 + (t285 * t316 + t288 * t294) * qJD(2)) - t292 * t311 + t295 * t294 + (t294 * pkin(8) + t297 * t292) * qJD(2)) * t287, t301 * qJD(4) + t296 * t272 + t323 * t273 + t300 * t312, t272, -t287 * t288 * t315 + t273 * t285, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:46
	% EndTime: 2019-10-09 22:12:47
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (541->121), mult. (1742->225), div. (0->0), fcn. (1827->12), ass. (0->77)
	t373 = sin(qJ(3));
	t376 = cos(qJ(3));
	t397 = -r_i_i_C(3) - pkin(9) + qJ(4);
	t415 = t373 * qJD(4) + (-pkin(3) * t373 + t397 * t376) * qJD(3);
	t369 = sin(pkin(10));
	t377 = cos(qJ(2));
	t410 = cos(pkin(10));
	t411 = cos(pkin(6));
	t387 = t411 * t410;
	t385 = t377 * t387;
	t374 = sin(qJ(2));
	t400 = qJD(2) * t374;
	t352 = -qJD(2) * t385 + t369 * t400;
	t370 = sin(pkin(6));
	t391 = t370 * t410;
	t414 = -qJD(3) * t391 - t352;
	t392 = t369 * t411;
	t359 = -t374 * t392 + t410 * t377;
	t357 = t369 * t377 + t374 * t387;
	t409 = t357 * t376;
	t368 = sin(pkin(11));
	t408 = t368 * t376;
	t407 = t369 * t370;
	t406 = t370 * t373;
	t405 = t370 * t377;
	t371 = cos(pkin(11));
	t404 = t371 * t376;
	t403 = t371 * t377;
	t402 = t374 * t376;
	t401 = t376 * t377;
	t399 = qJD(3) * t373;
	t396 = t374 * t406;
	t395 = t370 * t400;
	t394 = qJD(2) * t405;
	t393 = t377 * t399;
	t372 = sin(qJ(6));
	t375 = cos(qJ(6));
	t386 = t372 * r_i_i_C(1) + t375 * r_i_i_C(2) + qJ(5);
	t384 = t375 * r_i_i_C(1) - t372 * r_i_i_C(2) + pkin(4) + pkin(5);
	t345 = t359 * t376 + t369 * t406;
	t353 = t357 * qJD(2);
	t356 = t369 * t374 - t385;
	t383 = -t353 * t376 + t356 * t399;
	t355 = t359 * qJD(2);
	t358 = t410 * t374 + t377 * t392;
	t382 = -t355 * t376 + t358 * t399;
	t361 = t370 * t402 + t411 * t373;
	t381 = -t376 * pkin(3) - t397 * t373 - pkin(2);
	t379 = -t386 * t368 - t384 * t371 - pkin(3);
	t378 = t368 * qJD(5) + ((t368 * t375 - t371 * t372) * r_i_i_C(1) + (-t368 * t372 - t371 * t375) * r_i_i_C(2)) * qJD(6);
	t354 = t358 * qJD(2);
	t349 = (t368 * t374 + t371 * t401) * t370;
	t348 = (t368 * t401 - t371 * t374) * t370;
	t347 = -qJD(3) * t396 + (t411 * qJD(3) + t394) * t376;
	t346 = t361 * qJD(3) + t373 * t394;
	t343 = -t373 * t391 + t409;
	t341 = t361 * t371 - t368 * t405;
	t340 = t361 * t368 + t370 * t403;
	t337 = -t358 * t404 + t359 * t368;
	t336 = -t358 * t408 - t359 * t371;
	t335 = -t356 * t404 + t357 * t368;
	t334 = -t356 * t408 - t357 * t371;
	t333 = t347 * t371 + t368 * t395;
	t332 = t347 * t368 - t371 * t395;
	t331 = t345 * t371 + t358 * t368;
	t330 = t345 * t368 - t358 * t371;
	t329 = t343 * t371 + t356 * t368;
	t328 = t343 * t368 - t356 * t371;
	t327 = -t359 * t399 + (qJD(3) * t407 - t354) * t376;
	t326 = t345 * qJD(3) - t354 * t373;
	t325 = -t357 * t399 + t414 * t376;
	t324 = qJD(3) * t409 + t414 * t373;
	t319 = t327 * t371 + t355 * t368;
	t318 = t327 * t368 - t355 * t371;
	t317 = t325 * t371 + t353 * t368;
	t316 = t325 * t368 - t353 * t371;
	t1 = [0, -t354 * pkin(8) + t336 * qJD(5) + t386 * (t354 * t371 + t382 * t368) + t384 * (-t354 * t368 + t382 * t371) + ((t336 * t375 - t337 * t372) * r_i_i_C(1) + (-t336 * t372 - t337 * t375) * r_i_i_C(2)) * qJD(6) - t415 * t358 + t381 * t355, t345 * qJD(4) + t397 * t327 + t378 * (-t359 * t373 + t376 * t407) + t379 * t326, t326, t318, (t318 * t375 - t319 * t372) * r_i_i_C(1) + (-t318 * t372 - t319 * t375) * r_i_i_C(2) + ((-t330 * t372 - t331 * t375) * r_i_i_C(1) + (-t330 * t375 + t331 * t372) * r_i_i_C(2)) * qJD(6); 0, -t352 * pkin(8) + t334 * qJD(5) + t386 * (t352 * t371 + t383 * t368) + t384 * (-t352 * t368 + t383 * t371) + ((t334 * t375 - t335 * t372) * r_i_i_C(1) + (-t334 * t372 - t335 * t375) * r_i_i_C(2)) * qJD(6) - t415 * t356 + t381 * t353, t343 * qJD(4) + t397 * t325 + t378 * (-t357 * t373 - t376 * t391) + t379 * t324, t324, t316, (t316 * t375 - t317 * t372) * r_i_i_C(1) + (-t316 * t372 - t317 * t375) * r_i_i_C(2) + ((-t328 * t372 - t329 * t375) * r_i_i_C(1) + (-t328 * t375 + t329 * t372) * r_i_i_C(2)) * qJD(6); 0, t348 * qJD(5) + ((t348 * t375 - t349 * t372) * r_i_i_C(1) + (-t348 * t372 - t349 * t375) * r_i_i_C(2)) * qJD(6) + (-t386 * (t368 * t393 + (t368 * t402 + t403) * qJD(2)) + t384 * (-t371 * t393 + (t368 * t377 - t371 * t402) * qJD(2)) + t415 * t377 + (pkin(8) * t377 + t381 * t374) * qJD(2)) * t370, t361 * qJD(4) + t397 * t347 + t378 * (t411 * t376 - t396) + t379 * t346, t346, t332, (t332 * t375 - t333 * t372) * r_i_i_C(1) + (-t332 * t372 - t333 * t375) * r_i_i_C(2) + ((-t340 * t372 - t341 * t375) * r_i_i_C(1) + (-t340 * t375 + t341 * t372) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end