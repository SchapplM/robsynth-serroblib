% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:30
	% EndTime: 2019-10-09 22:27:30
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:30
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (116->40), mult. (277->78), div. (0->0), fcn. (252->10), ass. (0->33)
	t218 = r_i_i_C(3) + qJ(4) + pkin(8);
	t195 = sin(pkin(11));
	t196 = sin(pkin(6));
	t217 = t195 * t196;
	t197 = cos(pkin(11));
	t216 = t196 * t197;
	t200 = sin(qJ(3));
	t215 = t196 * t200;
	t201 = sin(qJ(2));
	t214 = t196 * t201;
	t198 = cos(pkin(6));
	t213 = t198 * t201;
	t203 = cos(qJ(2));
	t212 = t198 * t203;
	t211 = qJD(2) * t201;
	t210 = qJD(2) * t203;
	t209 = t195 * t211;
	t208 = t197 * t210;
	t194 = qJ(3) + pkin(12);
	t192 = sin(t194);
	t193 = cos(t194);
	t202 = cos(qJ(3));
	t207 = -pkin(3) * t202 - r_i_i_C(1) * t193 + r_i_i_C(2) * t192 - pkin(2);
	t186 = t195 * t203 + t197 * t213;
	t206 = t195 * t212 + t197 * t201;
	t205 = pkin(3) * t200 + r_i_i_C(1) * t192 + r_i_i_C(2) * t193;
	t204 = qJD(3) * t205;
	t188 = -t195 * t213 + t197 * t203;
	t184 = -t198 * t209 + t208;
	t183 = t206 * qJD(2);
	t182 = t186 * qJD(2);
	t181 = -t198 * t208 + t209;
	t1 = [0, qJD(4) * t188 - t218 * t183 + t207 * t184 + t206 * t204, t205 * t183 + ((-t188 * t193 - t192 * t217) * r_i_i_C(1) + (t188 * t192 - t193 * t217) * r_i_i_C(2) + (-t188 * t202 - t195 * t215) * pkin(3)) * qJD(3), t184, 0, 0; 0, qJD(4) * t186 - t218 * t181 + t207 * t182 - (-t195 * t201 + t197 * t212) * t204, t205 * t181 + ((-t186 * t193 + t192 * t216) * r_i_i_C(1) + (t186 * t192 + t193 * t216) * r_i_i_C(2) + (-t186 * t202 + t197 * t215) * pkin(3)) * qJD(3), t182, 0, 0; 0, (qJD(4) * t201 - t203 * t204 + (t207 * t201 + t218 * t203) * qJD(2)) * t196, -t205 * t196 * t210 + ((-t192 * t198 - t193 * t214) * r_i_i_C(1) + (t192 * t214 - t193 * t198) * r_i_i_C(2) + (-t198 * t200 - t202 * t214) * pkin(3)) * qJD(3), t196 * t211, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:31
	% EndTime: 2019-10-09 22:27:31
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (378->86), mult. (811->155), div. (0->0), fcn. (804->12), ass. (0->56)
	t307 = qJ(3) + pkin(12);
	t305 = sin(t307);
	t306 = cos(t307);
	t313 = sin(qJ(3));
	t312 = sin(qJ(5));
	t315 = cos(qJ(5));
	t327 = t315 * r_i_i_C(1) - t312 * r_i_i_C(2);
	t325 = pkin(4) + t327;
	t349 = pkin(9) + r_i_i_C(3);
	t351 = (t313 * pkin(3) + t325 * t305 - t349 * t306) * qJD(3);
	t326 = t312 * r_i_i_C(1) + t315 * r_i_i_C(2);
	t345 = cos(pkin(6));
	t308 = sin(pkin(11));
	t309 = sin(pkin(6));
	t344 = t308 * t309;
	t310 = cos(pkin(11));
	t343 = t309 * t310;
	t342 = t309 * t313;
	t314 = sin(qJ(2));
	t341 = t309 * t314;
	t317 = cos(qJ(2));
	t340 = t309 * t317;
	t339 = qJD(2) * t314;
	t338 = qJD(2) * t317;
	t337 = qJD(5) * t306;
	t336 = qJD(5) * t312;
	t335 = qJD(5) * t315;
	t334 = t309 * t338;
	t333 = t309 * t339;
	t332 = t314 * t345;
	t331 = t317 * t345;
	t329 = t308 * t332;
	t328 = t310 * t331;
	t298 = t308 * t317 + t310 * t332;
	t324 = -t298 * t305 - t306 * t343;
	t323 = -t298 * t306 + t305 * t343;
	t300 = t310 * t317 - t329;
	t322 = -t300 * t305 + t306 * t344;
	t290 = t300 * t306 + t305 * t344;
	t321 = qJD(5) * t326;
	t299 = t308 * t331 + t310 * t314;
	t320 = -t305 * t341 + t345 * t306;
	t292 = t345 * t305 + t306 * t341;
	t316 = cos(qJ(3));
	t319 = -t316 * pkin(3) - t349 * t305 - t325 * t306 - pkin(2);
	t318 = t326 * t337 + t351;
	t311 = -qJ(4) - pkin(8);
	t297 = t308 * t314 - t328;
	t296 = -qJD(2) * t329 + t310 * t338;
	t295 = t299 * qJD(2);
	t294 = t298 * qJD(2);
	t293 = -qJD(2) * t328 + t308 * t339;
	t286 = t320 * qJD(3) + t306 * t334;
	t284 = t322 * qJD(3) - t295 * t306;
	t282 = t324 * qJD(3) - t293 * t306;
	t1 = [0, (-t295 * t312 + t300 * t335) * r_i_i_C(1) + (-t295 * t315 - t300 * t336) * r_i_i_C(2) + t295 * t311 + t300 * qJD(4) + t319 * t296 + t318 * t299, t349 * t284 - t322 * t321 + t325 * (-t290 * qJD(3) + t295 * t305) + (t295 * t313 + (-t300 * t316 - t308 * t342) * qJD(3)) * pkin(3), t296, (-t284 * t312 + t296 * t315) * r_i_i_C(1) + (-t284 * t315 - t296 * t312) * r_i_i_C(2) + ((-t290 * t315 - t299 * t312) * r_i_i_C(1) + (t290 * t312 - t299 * t315) * r_i_i_C(2)) * qJD(5), 0; 0, (-t293 * t312 + t298 * t335) * r_i_i_C(1) + (-t293 * t315 - t298 * t336) * r_i_i_C(2) + t293 * t311 + t298 * qJD(4) + t319 * t294 + t318 * t297, t349 * t282 - t324 * t321 + t325 * (t323 * qJD(3) + t293 * t305) + (t293 * t313 + (-t298 * t316 + t310 * t342) * qJD(3)) * pkin(3), t294, (-t282 * t312 + t294 * t315) * r_i_i_C(1) + (-t282 * t315 - t294 * t312) * r_i_i_C(2) + ((-t297 * t312 + t315 * t323) * r_i_i_C(1) + (-t297 * t315 - t312 * t323) * r_i_i_C(2)) * qJD(5), 0; 0, ((t319 * qJD(2) + t327 * qJD(5) + qJD(4)) * t314 + (-qJD(2) * t311 - t351 + t326 * (qJD(2) - t337)) * t317) * t309, t349 * t286 - t320 * t321 + t325 * (-t292 * qJD(3) - t305 * t334) + (-t313 * t334 + (-t345 * t313 - t316 * t341) * qJD(3)) * pkin(3), t333, (-t286 * t312 + t315 * t333) * r_i_i_C(1) + (-t286 * t315 - t312 * t333) * r_i_i_C(2) + ((-t292 * t315 + t312 * t340) * r_i_i_C(1) + (t292 * t312 + t315 * t340) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:31
	% EndTime: 2019-10-09 22:27:31
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (695->102), mult. (1179->168), div. (0->0), fcn. (1178->14), ass. (0->67)
	t350 = qJ(3) + pkin(12);
	t345 = sin(t350);
	t346 = cos(t350);
	t356 = sin(qJ(3));
	t358 = cos(qJ(5));
	t351 = qJ(5) + qJ(6);
	t347 = sin(t351);
	t348 = cos(t351);
	t375 = r_i_i_C(1) * t348 - r_i_i_C(2) * t347;
	t371 = t358 * pkin(5) + pkin(4) + t375;
	t400 = r_i_i_C(3) + pkin(10) + pkin(9);
	t407 = (t356 * pkin(3) + t371 * t345 - t400 * t346) * qJD(3);
	t357 = sin(qJ(2));
	t360 = cos(qJ(2));
	t352 = sin(pkin(11));
	t399 = cos(pkin(6));
	t384 = t352 * t399;
	t398 = cos(pkin(11));
	t339 = -t357 * t384 + t398 * t360;
	t374 = r_i_i_C(1) * t347 + r_i_i_C(2) * t348;
	t349 = qJD(5) + qJD(6);
	t355 = sin(qJ(5));
	t402 = t355 * pkin(5);
	t406 = qJD(5) * t402 + t374 * t349;
	t397 = t347 * t349;
	t396 = t348 * t349;
	t353 = sin(pkin(6));
	t395 = t352 * t353;
	t394 = t353 * t357;
	t393 = t353 * t360;
	t372 = t399 * t398;
	t337 = t352 * t360 + t357 * t372;
	t333 = t337 * qJD(2);
	t383 = t353 * t398;
	t365 = -t337 * t346 + t345 * t383;
	t379 = -t349 * t365 - t333;
	t370 = t360 * t372;
	t389 = qJD(2) * t357;
	t332 = -qJD(2) * t370 + t352 * t389;
	t367 = -t337 * t345 - t346 * t383;
	t321 = t367 * qJD(3) - t332 * t346;
	t336 = t352 * t357 - t370;
	t381 = -t336 * t349 - t321;
	t392 = (t381 * t347 - t379 * t348) * r_i_i_C(1) + (t379 * t347 + t381 * t348) * r_i_i_C(2);
	t329 = t339 * t346 + t345 * t395;
	t335 = t339 * qJD(2);
	t378 = t329 * t349 - t335;
	t338 = t398 * t357 + t360 * t384;
	t334 = t338 * qJD(2);
	t369 = -t339 * t345 + t346 * t395;
	t323 = t369 * qJD(3) - t334 * t346;
	t380 = -t338 * t349 - t323;
	t391 = (t380 * t347 - t378 * t348) * r_i_i_C(1) + (t378 * t347 + t380 * t348) * r_i_i_C(2);
	t331 = t399 * t345 + t346 * t394;
	t385 = t353 * t389;
	t368 = -t331 * t349 + t385;
	t366 = -t345 * t394 + t399 * t346;
	t386 = qJD(2) * t393;
	t325 = t366 * qJD(3) + t346 * t386;
	t373 = t349 * t393 - t325;
	t390 = (t373 * t347 + t368 * t348) * r_i_i_C(1) + (-t368 * t347 + t373 * t348) * r_i_i_C(2);
	t388 = qJD(5) * t358;
	t359 = cos(qJ(3));
	t363 = -t359 * pkin(3) - t400 * t345 - t371 * t346 - pkin(2);
	t362 = t406 * t346 + t407;
	t354 = -qJ(4) - pkin(8);
	t1 = [0, (-t334 * t347 + t339 * t396) * r_i_i_C(1) + (-t334 * t348 - t339 * t397) * r_i_i_C(2) + t334 * t354 + t339 * qJD(4) + (-t334 * t355 + t339 * t388) * pkin(5) + t363 * t335 + t362 * t338, t400 * t323 - t406 * t369 + t371 * (-t329 * qJD(3) + t334 * t345) + (t334 * t356 + (-t339 * t359 - t356 * t395) * qJD(3)) * pkin(3), t335, (-t323 * t355 + t335 * t358 + (-t329 * t358 - t338 * t355) * qJD(5)) * pkin(5) + t391, t391; 0, (-t332 * t347 + t337 * t396) * r_i_i_C(1) + (-t332 * t348 - t337 * t397) * r_i_i_C(2) + t332 * t354 + t337 * qJD(4) + (-t332 * t355 + t337 * t388) * pkin(5) + t363 * t333 + t362 * t336, t400 * t321 - t406 * t367 + t371 * (t365 * qJD(3) + t332 * t345) + (t332 * t356 + (-t337 * t359 + t356 * t383) * qJD(3)) * pkin(3), t333, (-t321 * t355 + t333 * t358 + (-t336 * t355 + t358 * t365) * qJD(5)) * pkin(5) + t392, t392; 0, ((pkin(5) * t388 + t363 * qJD(2) + t375 * t349 + qJD(4)) * t357 + (-qJD(2) * t354 + (-qJD(5) * t346 + qJD(2)) * t402 - t407 + t374 * (-t346 * t349 + qJD(2))) * t360) * t353, t400 * t325 - t406 * t366 + t371 * (-t331 * qJD(3) - t345 * t386) + (-t356 * t386 + (-t399 * t356 - t359 * t394) * qJD(3)) * pkin(3), t385, (t358 * t385 - t325 * t355 + (-t331 * t358 + t355 * t393) * qJD(5)) * pkin(5) + t390, t390;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end