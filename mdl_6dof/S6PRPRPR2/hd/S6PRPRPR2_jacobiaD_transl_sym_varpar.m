% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
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
	% StartTime: 2019-10-09 21:31:52
	% EndTime: 2019-10-09 21:31:52
	% DurationCPUTime: 0.25s
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
	% StartTime: 2019-10-09 21:31:53
	% EndTime: 2019-10-09 21:31:53
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (297->50), mult. (960->98), div. (0->0), fcn. (1004->12), ass. (0->46)
	t335 = sin(pkin(11));
	t339 = cos(pkin(11));
	t345 = cos(qJ(2));
	t359 = qJD(2) * t345;
	t343 = sin(qJ(2));
	t360 = qJD(2) * t343;
	t368 = t335 * t360 - t339 * t359;
	t341 = cos(pkin(6));
	t353 = t335 * t345 + t343 * t339;
	t324 = t353 * t341;
	t336 = sin(pkin(10));
	t340 = cos(pkin(10));
	t352 = t343 * t335 - t339 * t345;
	t313 = t324 * t340 - t336 * t352;
	t344 = cos(qJ(4));
	t337 = sin(pkin(6));
	t342 = sin(qJ(4));
	t363 = t337 * t342;
	t367 = -t313 * t344 + t340 * t363;
	t366 = r_i_i_C(3) + qJ(5);
	t365 = pkin(2) * qJD(2);
	t362 = t337 * t344;
	t361 = t341 * t343;
	t319 = t368 * t341;
	t326 = -t335 * t359 - t339 * t360;
	t355 = t319 * t340 - t326 * t336;
	t310 = t336 * t319 + t326 * t340;
	t322 = t353 * t337;
	t354 = t322 * t344 + t341 * t342;
	t334 = sin(pkin(12));
	t338 = cos(pkin(12));
	t351 = r_i_i_C(1) * t338 - r_i_i_C(2) * t334 + pkin(4);
	t350 = -r_i_i_C(1) * t334 - r_i_i_C(2) * t338 - pkin(8);
	t315 = -t324 * t336 - t340 * t352;
	t349 = t315 * t344 + t336 * t363;
	t348 = qJD(2) * t353;
	t347 = t366 * t342 + t351 * t344 + pkin(3);
	t346 = t342 * qJD(5) + (-t351 * t342 + t366 * t344) * qJD(4);
	t325 = t352 * qJD(2);
	t323 = t352 * t341;
	t320 = t341 * t348;
	t317 = t368 * t337;
	t304 = t354 * qJD(4) - t317 * t342;
	t302 = t349 * qJD(4) + t310 * t342;
	t300 = -t367 * qJD(4) - t342 * t355;
	t1 = [0, (t336 * t361 - t340 * t345) * t365 - t350 * t310 + t346 * (t323 * t336 - t340 * t353) + t347 * (t320 * t336 + t325 * t340), 0, t349 * qJD(5) + t366 * (t310 * t344 + (-t315 * t342 + t336 * t362) * qJD(4)) - t351 * t302, t302, 0; 0, (-t336 * t345 - t340 * t361) * t365 + t350 * t355 + t346 * (-t323 * t340 - t336 * t353) + t347 * (-t320 * t340 + t325 * t336), 0, -t367 * qJD(5) + t366 * (-t355 * t344 + (-t313 * t342 - t340 * t362) * qJD(4)) - t351 * t300, t300, 0; 0, t350 * t317 + (-pkin(2) * t360 - t346 * t352 - t347 * t348) * t337, 0, t354 * qJD(5) + t366 * (-t317 * t344 + (-t322 * t342 + t341 * t344) * qJD(4)) - t351 * t304, t304, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:53
	% EndTime: 2019-10-09 21:31:54
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (566->80), mult. (1521->147), div. (0->0), fcn. (1646->14), ass. (0->62)
	t395 = sin(qJ(4));
	t397 = cos(qJ(4));
	t387 = pkin(12) + qJ(6);
	t385 = sin(t387);
	t386 = cos(t387);
	t413 = t385 * r_i_i_C(1) + t386 * r_i_i_C(2);
	t405 = qJD(6) * t413;
	t414 = t386 * r_i_i_C(1) - t385 * r_i_i_C(2);
	t410 = cos(pkin(12)) * pkin(5) + pkin(4) + t414;
	t429 = r_i_i_C(3) + pkin(9) + qJ(5);
	t434 = (t410 * t395 - t429 * t397) * qJD(4) - t395 * qJD(5) + t397 * t405;
	t393 = cos(pkin(6));
	t389 = sin(pkin(11));
	t396 = sin(qJ(2));
	t427 = cos(pkin(11));
	t430 = cos(qJ(2));
	t403 = t430 * t389 + t396 * t427;
	t372 = t403 * t393;
	t415 = t430 * t427;
	t421 = qJD(2) * t396;
	t432 = -qJD(2) * t415 + t389 * t421;
	t402 = -t396 * t389 + t415;
	t390 = sin(pkin(10));
	t392 = cos(pkin(10));
	t358 = t392 * t372 + t390 * t402;
	t391 = sin(pkin(6));
	t425 = t391 * t395;
	t347 = t358 * t397 - t392 * t425;
	t399 = t429 * t395 + t410 * t397 + pkin(3);
	t428 = pkin(2) * qJD(2);
	t424 = t391 * t397;
	t423 = t393 * t396;
	t369 = t432 * t393;
	t374 = t403 * qJD(2);
	t412 = t392 * t369 + t390 * t374;
	t354 = t390 * t369 - t392 * t374;
	t371 = t403 * t391;
	t363 = t371 * t397 + t393 * t395;
	t411 = -t371 * t395 + t393 * t397;
	t361 = -t390 * t372 + t392 * t402;
	t408 = -t358 * t395 - t392 * t424;
	t407 = -t361 * t395 + t390 * t424;
	t349 = t361 * t397 + t390 * t425;
	t406 = qJD(6) * t414;
	t404 = -sin(pkin(12)) * pkin(5) - pkin(8) - t413;
	t401 = t402 * t393;
	t400 = qJD(2) * t372;
	t373 = t402 * qJD(2);
	t370 = t402 * t391;
	t368 = qJD(2) * t371;
	t367 = t432 * t391;
	t360 = -t390 * t401 - t392 * t403;
	t357 = -t390 * t403 + t392 * t401;
	t353 = -t392 * t373 + t390 * t400;
	t350 = -t390 * t373 - t392 * t400;
	t345 = t411 * qJD(4) - t367 * t397;
	t344 = t363 * qJD(4) - t367 * t395;
	t343 = t407 * qJD(4) + t354 * t397;
	t342 = t349 * qJD(4) + t354 * t395;
	t341 = t408 * qJD(4) - t397 * t412;
	t340 = t347 * qJD(4) - t395 * t412;
	t1 = [0, t361 * t406 - t404 * t354 + (t390 * t423 - t430 * t392) * t428 + t399 * t353 - t434 * t360, 0, t349 * qJD(5) - t410 * t342 + t429 * t343 - t407 * t405, t342, (-t343 * t385 - t353 * t386) * r_i_i_C(1) + (-t343 * t386 + t353 * t385) * r_i_i_C(2) + ((-t349 * t386 + t360 * t385) * r_i_i_C(1) + (t349 * t385 + t360 * t386) * r_i_i_C(2)) * qJD(6); 0, t358 * t406 + t404 * t412 + (-t430 * t390 - t392 * t423) * t428 + t399 * t350 - t434 * t357, 0, t347 * qJD(5) - t410 * t340 + t429 * t341 - t408 * t405, t340, (-t341 * t385 - t350 * t386) * r_i_i_C(1) + (-t341 * t386 + t350 * t385) * r_i_i_C(2) + ((-t347 * t386 + t357 * t385) * r_i_i_C(1) + (t347 * t385 + t357 * t386) * r_i_i_C(2)) * qJD(6); 0, -t391 * pkin(2) * t421 + t404 * t367 - t399 * t368 - t434 * t370 + t371 * t406, 0, t363 * qJD(5) - t410 * t344 + t429 * t345 - t411 * t405, t344, (-t345 * t385 + t368 * t386) * r_i_i_C(1) + (-t345 * t386 - t368 * t385) * r_i_i_C(2) + ((-t363 * t386 + t370 * t385) * r_i_i_C(1) + (t363 * t385 + t370 * t386) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end