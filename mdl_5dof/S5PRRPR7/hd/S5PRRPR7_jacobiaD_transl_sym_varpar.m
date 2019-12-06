% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(5));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(9));
	t50 = sin(pkin(9));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(5)) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:33
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(9));
	t182 = cos(pkin(9));
	t185 = sin(qJ(2));
	t183 = cos(pkin(5));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(7) + r_i_i_C(3);
	t181 = sin(pkin(5));
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
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:33
	% EndTime: 2019-12-05 16:38:33
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (152->37), mult. (504->70), div. (0->0), fcn. (488->10), ass. (0->34)
	t256 = sin(pkin(9));
	t259 = cos(pkin(9));
	t262 = sin(qJ(2));
	t260 = cos(pkin(5));
	t264 = cos(qJ(2));
	t276 = t260 * t264;
	t285 = -t256 * t262 + t259 * t276;
	t277 = t260 * t262;
	t250 = t256 * t264 + t259 * t277;
	t263 = cos(qJ(3));
	t257 = sin(pkin(5));
	t261 = sin(qJ(3));
	t279 = t257 * t261;
	t284 = -t250 * t263 + t259 * t279;
	t255 = sin(pkin(10));
	t258 = cos(pkin(10));
	t272 = t258 * r_i_i_C(1) - t255 * r_i_i_C(2) + pkin(3);
	t282 = r_i_i_C(3) + qJ(4);
	t283 = t282 * t261 + t272 * t263 + pkin(2);
	t278 = t257 * t263;
	t273 = qJD(2) * t257 * t264;
	t271 = t255 * r_i_i_C(1) + t258 * r_i_i_C(2) + pkin(7);
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
	t1 = [0, -t271 * t247 - t265 * t269 + t268 * t266, t270 * qJD(4) + t282 * (-t247 * t263 + (t256 * t278 + t261 * t268) * qJD(3)) - t272 * t241, t241, 0; 0, t271 * t245 - t250 * t266 + t265 * t285, -t284 * qJD(4) + t282 * (t245 * t263 + (-t250 * t261 - t259 * t278) * qJD(3)) - t272 * t239, t239, 0; 0, (t265 * t264 + (-t283 * t262 + t271 * t264) * qJD(2)) * t257, t267 * qJD(4) + t282 * (t263 * t273 + (t260 * t263 - t262 * t279) * qJD(3)) - t272 * t243, t243, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:33
	% EndTime: 2019-12-05 16:38:34
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (422->111), mult. (1359->210), div. (0->0), fcn. (1420->12), ass. (0->70)
	t386 = sin(qJ(3));
	t389 = cos(qJ(3));
	t385 = sin(qJ(5));
	t388 = cos(qJ(5));
	t399 = t385 * r_i_i_C(1) + t388 * r_i_i_C(2) + qJ(4);
	t401 = t388 * r_i_i_C(1) - t385 * r_i_i_C(2);
	t391 = -(-t386 * pkin(3) + t399 * t389) * qJD(3) - (t401 * qJD(5) + qJD(4)) * t386;
	t382 = sin(pkin(9));
	t387 = sin(qJ(2));
	t390 = cos(qJ(2));
	t421 = cos(pkin(9));
	t422 = cos(pkin(5));
	t400 = t422 * t421;
	t370 = t382 * t390 + t387 * t400;
	t383 = sin(pkin(5));
	t405 = t383 * t421;
	t359 = t370 * t389 - t386 * t405;
	t406 = t382 * t422;
	t372 = -t387 * t406 + t421 * t390;
	t424 = r_i_i_C(3) + pkin(8);
	t419 = t383 * t386;
	t418 = t383 * t390;
	t384 = cos(pkin(10));
	t417 = t384 * t385;
	t416 = t384 * t388;
	t415 = t384 * t389;
	t414 = t384 * t390;
	t413 = t387 * t389;
	t412 = qJD(2) * t387;
	t411 = qJD(3) * t386;
	t410 = qJD(5) * t385;
	t409 = qJD(5) * t388;
	t408 = qJD(2) * t418;
	t407 = t390 * t411;
	t398 = t390 * t400;
	t397 = t382 * t383 * t389 - t372 * t386;
	t361 = t372 * t389 + t382 * t419;
	t366 = t370 * qJD(2);
	t369 = t382 * t387 - t398;
	t396 = -t366 * t389 + t369 * t411;
	t368 = t372 * qJD(2);
	t371 = t421 * t387 + t390 * t406;
	t395 = -t368 * t389 + t371 * t411;
	t373 = t387 * t419 - t422 * t389;
	t374 = t383 * t413 + t422 * t386;
	t394 = -t370 * t386 - t389 * t405;
	t393 = -t389 * pkin(3) - t399 * t386 - pkin(2);
	t381 = sin(pkin(10));
	t392 = -pkin(3) - t424 * t381 + (-pkin(4) - t401) * t384;
	t367 = t371 * qJD(2);
	t365 = -qJD(2) * t398 + t382 * t412;
	t364 = (t381 * t387 + t389 * t414) * t383;
	t363 = -t373 * qJD(3) + t389 * t408;
	t362 = t374 * qJD(3) + t386 * t408;
	t357 = t374 * t384 - t381 * t418;
	t356 = (-t384 * t407 + (t381 * t390 - t384 * t413) * qJD(2)) * t383;
	t354 = -t371 * t415 + t372 * t381;
	t353 = -t369 * t415 + t370 * t381;
	t352 = t381 * t383 * t412 + t363 * t384;
	t351 = t361 * t384 + t371 * t381;
	t350 = t359 * t384 + t369 * t381;
	t349 = t397 * qJD(3) - t367 * t389;
	t348 = t361 * qJD(3) - t367 * t386;
	t347 = t394 * qJD(3) - t365 * t389;
	t346 = t359 * qJD(3) - t365 * t386;
	t345 = -t367 * t381 + t395 * t384;
	t343 = -t365 * t381 + t396 * t384;
	t341 = t349 * t384 + t368 * t381;
	t340 = t347 * t384 + t366 * t381;
	t1 = [0, (t345 * t388 - t354 * t410) * r_i_i_C(1) + (-t345 * t385 - t354 * t409) * r_i_i_C(2) + t345 * pkin(4) - t367 * pkin(7) + t424 * (t367 * t384 + t395 * t381) + t393 * t368 + t391 * t371, t361 * qJD(4) + t399 * t349 + ((t361 * t388 - t397 * t417) * r_i_i_C(1) + (-t361 * t385 - t397 * t416) * r_i_i_C(2)) * qJD(5) + t392 * t348, t348, (-t341 * t385 + t348 * t388) * r_i_i_C(1) + (-t341 * t388 - t348 * t385) * r_i_i_C(2) + ((-t351 * t388 + t385 * t397) * r_i_i_C(1) + (t351 * t385 + t388 * t397) * r_i_i_C(2)) * qJD(5); 0, (t343 * t388 - t353 * t410) * r_i_i_C(1) + (-t343 * t385 - t353 * t409) * r_i_i_C(2) + t343 * pkin(4) - t365 * pkin(7) + t424 * (t365 * t384 + t396 * t381) + t393 * t366 + t391 * t369, t359 * qJD(4) + t399 * t347 + ((t359 * t388 - t394 * t417) * r_i_i_C(1) + (-t359 * t385 - t394 * t416) * r_i_i_C(2)) * qJD(5) + t392 * t346, t346, (-t340 * t385 + t346 * t388) * r_i_i_C(1) + (-t340 * t388 - t346 * t385) * r_i_i_C(2) + ((-t350 * t388 + t385 * t394) * r_i_i_C(1) + (t350 * t385 + t388 * t394) * r_i_i_C(2)) * qJD(5); 0, (t356 * t388 - t364 * t410) * r_i_i_C(1) + (-t356 * t385 - t364 * t409) * r_i_i_C(2) + t356 * pkin(4) + (-t424 * (t381 * t407 + (t381 * t413 + t414) * qJD(2)) + t393 * t412 + (qJD(2) * pkin(7) - t391) * t390) * t383, t374 * qJD(4) + t399 * t363 + ((t373 * t417 + t374 * t388) * r_i_i_C(1) + (t373 * t416 - t374 * t385) * r_i_i_C(2)) * qJD(5) + t392 * t362, t362, (-t352 * t385 + t362 * t388) * r_i_i_C(1) + (-t352 * t388 - t362 * t385) * r_i_i_C(2) + ((-t357 * t388 - t373 * t385) * r_i_i_C(1) + (t357 * t385 - t373 * t388) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end