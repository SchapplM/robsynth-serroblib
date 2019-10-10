% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
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
	% StartTime: 2019-10-09 22:20:06
	% EndTime: 2019-10-09 22:20:06
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-09 22:20:06
	% EndTime: 2019-10-09 22:20:07
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
	% StartTime: 2019-10-09 22:20:07
	% EndTime: 2019-10-09 22:20:07
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (336->69), mult. (819->125), div. (0->0), fcn. (824->12), ass. (0->51)
	t305 = sin(qJ(3));
	t307 = cos(qJ(3));
	t300 = pkin(11) + qJ(5);
	t298 = sin(t300);
	t299 = cos(t300);
	t320 = t298 * r_i_i_C(1) + t299 * r_i_i_C(2);
	t315 = qJD(5) * t320;
	t321 = r_i_i_C(1) * t299 - r_i_i_C(2) * t298;
	t318 = cos(pkin(11)) * pkin(4) + pkin(3) + t321;
	t337 = r_i_i_C(3) + pkin(9) + qJ(4);
	t309 = (t318 * t305 - t337 * t307) * qJD(3) - t305 * qJD(4) + t307 * t315;
	t302 = sin(pkin(10));
	t306 = sin(qJ(2));
	t308 = cos(qJ(2));
	t335 = cos(pkin(10));
	t336 = cos(pkin(6));
	t319 = t336 * t335;
	t288 = t302 * t308 + t306 * t319;
	t303 = sin(pkin(6));
	t325 = t303 * t335;
	t278 = t288 * t307 - t305 * t325;
	t326 = t302 * t336;
	t290 = -t306 * t326 + t335 * t308;
	t333 = t303 * t305;
	t332 = t303 * t307;
	t331 = t303 * t308;
	t330 = qJD(2) * t306;
	t328 = t303 * t330;
	t327 = qJD(2) * t331;
	t317 = t308 * t319;
	t316 = -t290 * t305 + t302 * t332;
	t280 = t290 * t307 + t302 * t333;
	t314 = t321 * qJD(5);
	t313 = -t288 * t305 - t307 * t325;
	t312 = -t306 * t333 + t336 * t307;
	t292 = t336 * t305 + t306 * t332;
	t311 = sin(pkin(11)) * pkin(4) + pkin(8) + t320;
	t289 = t335 * t306 + t308 * t326;
	t310 = -t337 * t305 - t318 * t307 - pkin(2);
	t287 = t302 * t306 - t317;
	t286 = t290 * qJD(2);
	t285 = t289 * qJD(2);
	t284 = t288 * qJD(2);
	t283 = -qJD(2) * t317 + t302 * t330;
	t282 = t312 * qJD(3) + t307 * t327;
	t281 = t292 * qJD(3) + t305 * t327;
	t276 = t316 * qJD(3) - t285 * t307;
	t275 = t280 * qJD(3) - t285 * t305;
	t274 = t313 * qJD(3) - t283 * t307;
	t273 = t278 * qJD(3) - t283 * t305;
	t1 = [0, -t311 * t285 + t310 * t286 + t309 * t289 + t290 * t314, qJD(4) * t280 - t318 * t275 + t337 * t276 - t316 * t315, t275, (-t276 * t298 + t286 * t299) * r_i_i_C(1) + (-t276 * t299 - t286 * t298) * r_i_i_C(2) + ((-t280 * t299 - t289 * t298) * r_i_i_C(1) + (t280 * t298 - t289 * t299) * r_i_i_C(2)) * qJD(5), 0; 0, -t311 * t283 + t310 * t284 + t309 * t287 + t288 * t314, qJD(4) * t278 - t318 * t273 + t337 * t274 - t313 * t315, t273, (-t274 * t298 + t284 * t299) * r_i_i_C(1) + (-t274 * t299 - t284 * t298) * r_i_i_C(2) + ((-t278 * t299 - t287 * t298) * r_i_i_C(1) + (t278 * t298 - t287 * t299) * r_i_i_C(2)) * qJD(5), 0; 0, ((t310 * qJD(2) + t314) * t306 + (t311 * qJD(2) - t309) * t308) * t303, qJD(4) * t292 - t318 * t281 + t337 * t282 - t312 * t315, t281, (-t282 * t298 + t299 * t328) * r_i_i_C(1) + (-t282 * t299 - t298 * t328) * r_i_i_C(2) + ((-t292 * t299 + t298 * t331) * r_i_i_C(1) + (t292 * t298 + t299 * t331) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:08
	% EndTime: 2019-10-09 22:20:09
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (642->95), mult. (1461->159), div. (0->0), fcn. (1517->12), ass. (0->67)
	t368 = cos(pkin(11)) * pkin(4) + pkin(3);
	t376 = sin(qJ(3));
	t378 = cos(qJ(3));
	t371 = pkin(11) + qJ(5);
	t369 = sin(t371);
	t409 = qJD(6) * t369;
	t422 = r_i_i_C(2) + pkin(9) + qJ(4);
	t427 = (-t368 * t376 + t422 * t378) * qJD(3) + t376 * qJD(4) + t378 * t409;
	t373 = sin(pkin(10));
	t377 = sin(qJ(2));
	t379 = cos(qJ(2));
	t419 = cos(pkin(10));
	t420 = cos(pkin(6));
	t394 = t420 * t419;
	t357 = t373 * t379 + t377 * t394;
	t374 = sin(pkin(6));
	t401 = t374 * t419;
	t345 = t357 * t378 - t376 * t401;
	t402 = t373 * t420;
	t359 = -t377 * t402 + t419 * t379;
	t414 = t374 * t378;
	t361 = t420 * t376 + t377 * t414;
	t370 = cos(t371);
	t413 = t374 * t379;
	t426 = -t361 * t370 + t369 * t413;
	t411 = qJD(3) * t376;
	t425 = (qJD(2) * t378 - qJD(5)) * t377 + t379 * t411;
	t423 = r_i_i_C(1) + pkin(5);
	t421 = r_i_i_C(3) + qJ(6);
	t415 = t374 * t376;
	t412 = qJD(2) * t377;
	t410 = qJD(5) * t378;
	t408 = t370 * qJD(6);
	t406 = sin(pkin(11)) * pkin(4) + pkin(8);
	t405 = t374 * t412;
	t404 = qJD(2) * t413;
	t390 = t379 * t394;
	t352 = -qJD(2) * t390 + t373 * t412;
	t356 = t373 * t377 - t390;
	t396 = t356 * t410 - t352;
	t358 = t419 * t377 + t379 * t402;
	t354 = t358 * qJD(2);
	t395 = t358 * t410 - t354;
	t393 = t345 * t370 + t356 * t369;
	t347 = t359 * t378 + t373 * t415;
	t392 = t347 * t370 + t358 * t369;
	t391 = (qJD(2) - t410) * t379;
	t389 = -t359 * t376 + t373 * t414;
	t387 = -t378 * t368 - t422 * t376 - pkin(2);
	t386 = -t357 * t376 - t378 * t401;
	t385 = -t377 * t415 + t420 * t378;
	t384 = -t421 * t369 - t423 * t370 - t368;
	t353 = t357 * qJD(2);
	t383 = qJD(5) * t357 - t353 * t378 + t356 * t411;
	t355 = t359 * qJD(2);
	t382 = qJD(5) * t359 - t355 * t378 + t358 * t411;
	t381 = t409 + (-t423 * t369 + t421 * t370) * qJD(5);
	t349 = t385 * qJD(3) + t378 * t404;
	t348 = t361 * qJD(3) + t376 * t404;
	t343 = t389 * qJD(3) - t354 * t378;
	t342 = t347 * qJD(3) - t354 * t376;
	t341 = t386 * qJD(3) - t352 * t378;
	t340 = t345 * qJD(3) - t352 * t376;
	t336 = -t426 * qJD(5) + t349 * t369 - t370 * t405;
	t330 = t392 * qJD(5) + t343 * t369 - t355 * t370;
	t328 = t393 * qJD(5) + t341 * t369 - t353 * t370;
	t1 = [0, -t359 * t408 - t406 * t354 + t423 * (t395 * t369 + t382 * t370) + t421 * (t382 * t369 - t395 * t370) - t427 * t358 + t387 * t355, t347 * qJD(4) + t384 * t342 + t422 * t343 + t381 * t389, t342, t392 * qJD(6) + t421 * (t343 * t370 + t355 * t369 + (-t347 * t369 + t358 * t370) * qJD(5)) - t423 * t330, t330; 0, -t357 * t408 - t406 * t352 + t423 * (t396 * t369 + t383 * t370) + t421 * (t383 * t369 - t396 * t370) - t427 * t356 + t387 * t353, t345 * qJD(4) + t384 * t340 + t422 * t341 + t381 * t386, t340, t393 * qJD(6) + t421 * (t341 * t370 + t353 * t369 + (-t345 * t369 + t356 * t370) * qJD(5)) - t423 * t328, t328; 0, (t423 * (t369 * t391 - t425 * t370) - t421 * (t425 * t369 + t370 * t391) - t377 * t408 + t427 * t379 + (t387 * t377 + t406 * t379) * qJD(2)) * t374, t361 * qJD(4) + t384 * t348 + t422 * t349 + t381 * t385, t348, -t426 * qJD(6) + t421 * (t369 * t405 + t349 * t370 + (-t361 * t369 - t370 * t413) * qJD(5)) - t423 * t336, t336;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end