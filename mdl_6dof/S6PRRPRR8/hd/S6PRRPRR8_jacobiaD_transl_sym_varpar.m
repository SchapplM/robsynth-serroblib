% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(12));
	t50 = sin(pkin(12));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:11
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (95->48), mult. (333->96), div. (0->0), fcn. (328->10), ass. (0->36)
	t218 = sin(pkin(7));
	t248 = (pkin(9) + r_i_i_C(3)) * t218;
	t217 = sin(pkin(12));
	t220 = cos(pkin(12));
	t224 = sin(qJ(2));
	t222 = cos(pkin(6));
	t226 = cos(qJ(2));
	t241 = t222 * t226;
	t211 = -t217 * t224 + t220 * t241;
	t219 = sin(pkin(6));
	t245 = t218 * t219;
	t221 = cos(pkin(7));
	t223 = sin(qJ(3));
	t244 = t221 * t223;
	t225 = cos(qJ(3));
	t243 = t221 * t225;
	t242 = t222 * t224;
	t240 = t223 * t224;
	t239 = t223 * t226;
	t238 = t224 * t225;
	t237 = t225 * t226;
	t236 = t223 * t245;
	t235 = t225 * t245;
	t233 = r_i_i_C(1) * t223 + r_i_i_C(2) * t225;
	t232 = r_i_i_C(1) * t225 - r_i_i_C(2) * t223 + pkin(2);
	t212 = t217 * t226 + t220 * t242;
	t231 = t217 * t241 + t220 * t224;
	t230 = t217 * t242 - t220 * t226;
	t229 = t233 * t221 - t248;
	t228 = (-t221 * t238 - t239) * r_i_i_C(1) + (t221 * t240 - t237) * r_i_i_C(2);
	t227 = (-t221 * t239 - t238) * r_i_i_C(1) + (-t221 * t237 + t240) * r_i_i_C(2);
	t210 = t230 * qJD(2);
	t209 = t231 * qJD(2);
	t208 = t212 * qJD(2);
	t207 = t211 * qJD(2);
	t1 = [0, t232 * t210 + t229 * t209 + ((t223 * t231 + t230 * t243) * r_i_i_C(1) + (t225 * t231 - t230 * t244) * r_i_i_C(2)) * qJD(3), (t209 * t223 + t210 * t243) * r_i_i_C(1) + (t209 * t225 - t210 * t244) * r_i_i_C(2) + ((-t217 * t236 + t225 * t230 + t231 * t244) * r_i_i_C(1) + (-t217 * t235 - t223 * t230 + t231 * t243) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, -t232 * t208 - t229 * t207 + ((-t211 * t223 - t212 * t243) * r_i_i_C(1) + (-t211 * t225 + t212 * t244) * r_i_i_C(2)) * qJD(3), (-t207 * t223 - t208 * t243) * r_i_i_C(1) + (-t207 * t225 + t208 * t244) * r_i_i_C(2) + ((-t211 * t244 - t212 * t225 + t220 * t236) * r_i_i_C(1) + (-t211 * t243 + t212 * t223 + t220 * t235) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (t228 * qJD(3) + (-t224 * pkin(2) + t226 * t248 + t227) * qJD(2)) * t219, -t233 * t222 * t218 * qJD(3) + (t228 * qJD(2) + t227 * qJD(3)) * t219, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:11
	% EndTime: 2019-10-09 22:39:12
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (224->62), mult. (748->112), div. (0->0), fcn. (768->10), ass. (0->46)
	t309 = sin(pkin(7));
	t330 = (pkin(9) + r_i_i_C(1)) * t309;
	t308 = sin(pkin(12));
	t311 = cos(pkin(12));
	t315 = sin(qJ(2));
	t313 = cos(pkin(6));
	t317 = cos(qJ(2));
	t335 = t313 * t317;
	t302 = -t308 * t315 + t311 * t335;
	t314 = sin(qJ(3));
	t316 = cos(qJ(3));
	t336 = t313 * t315;
	t322 = t308 * t336 - t311 * t317;
	t312 = cos(pkin(7));
	t323 = t308 * t335 + t311 * t315;
	t310 = sin(pkin(6));
	t340 = t309 * t310;
	t324 = t308 * t340 - t312 * t323;
	t348 = t324 * t314 - t316 * t322;
	t303 = t308 * t317 + t311 * t336;
	t326 = -t302 * t312 + t311 * t340;
	t347 = -t303 * t316 + t326 * t314;
	t346 = pkin(3) - r_i_i_C(2);
	t344 = r_i_i_C(3) + qJ(4);
	t339 = t309 * t313;
	t338 = t312 * t314;
	t337 = t312 * t316;
	t334 = t314 * t315;
	t333 = t314 * t317;
	t332 = t315 * t316;
	t331 = t316 * t317;
	t328 = qJD(3) * t339;
	t327 = -t302 * t314 - t303 * t337;
	t325 = t314 * t323 + t322 * t337;
	t321 = t312 * t331 - t334;
	t320 = t312 * t332 + t333;
	t319 = t312 * t333 + t332;
	t318 = t312 * t334 - t331;
	t301 = t322 * qJD(2);
	t300 = t323 * qJD(2);
	t299 = t303 * qJD(2);
	t298 = t302 * qJD(2);
	t294 = t314 * t328 + (t320 * qJD(2) + t319 * qJD(3)) * t310;
	t288 = t348 * qJD(3) - t300 * t314 - t301 * t337;
	t286 = -t347 * qJD(3) + t298 * t314 + t299 * t337;
	t1 = [0, -t325 * qJD(4) + t301 * pkin(2) - t300 * t330 + t346 * (t325 * qJD(3) + t300 * t338 + t301 * t316) + t344 * (-t300 * t337 + t301 * t314 + (-t316 * t323 + t322 * t338) * qJD(3)), t348 * qJD(4) + t344 * (t301 * t338 - t300 * t316 + (t314 * t322 + t324 * t316) * qJD(3)) - t346 * t288, t288, 0, 0; 0, -t327 * qJD(4) - t299 * pkin(2) + t298 * t330 + t346 * (t327 * qJD(3) - t298 * t338 - t299 * t316) + t344 * (t298 * t337 - t299 * t314 + (t302 * t316 - t303 * t338) * qJD(3)), -t347 * qJD(4) + t344 * (-t299 * t338 + t298 * t316 + (-t303 * t314 - t326 * t316) * qJD(3)) - t346 * t286, t286, 0, 0; 0, (-t346 * (t319 * qJD(2) + t320 * qJD(3)) - t344 * (-t321 * qJD(2) + t318 * qJD(3)) + t320 * qJD(4) + (-t315 * pkin(2) + t317 * t330) * qJD(2)) * t310, -(-t319 * t310 - t314 * t339) * qJD(4) + t344 * (t316 * t328 + (-t318 * qJD(2) + t321 * qJD(3)) * t310) - t346 * t294, t294, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:12
	% EndTime: 2019-10-09 22:39:13
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (471->110), mult. (1558->203), div. (0->0), fcn. (1646->12), ass. (0->72)
	t403 = sin(qJ(3));
	t406 = cos(qJ(3));
	t396 = sin(pkin(12));
	t399 = cos(pkin(12));
	t407 = cos(qJ(2));
	t401 = cos(pkin(6));
	t404 = sin(qJ(2));
	t436 = t401 * t404;
	t414 = t396 * t436 - t399 * t407;
	t400 = cos(pkin(7));
	t435 = t401 * t407;
	t415 = t396 * t435 + t399 * t404;
	t397 = sin(pkin(7));
	t398 = sin(pkin(6));
	t444 = t397 * t398;
	t416 = t396 * t444 - t400 * t415;
	t448 = t403 * t414 + t416 * t406;
	t388 = t396 * t407 + t399 * t436;
	t426 = t399 * t435;
	t387 = -t396 * t404 + t426;
	t440 = t398 * t399;
	t417 = -t387 * t400 + t397 * t440;
	t367 = t388 * t403 + t417 * t406;
	t446 = t388 * t406;
	t402 = sin(qJ(5));
	t443 = t397 * t402;
	t442 = t397 * t403;
	t405 = cos(qJ(5));
	t441 = t397 * t405;
	t439 = t398 * t400;
	t438 = t400 * t403;
	t437 = t400 * t406;
	t434 = t403 * t404;
	t433 = t403 * t407;
	t432 = t404 * t406;
	t431 = t406 * t407;
	t430 = qJD(2) * t404;
	t429 = qJD(5) * t402;
	t428 = qJD(5) * t405;
	t427 = r_i_i_C(3) + pkin(10) + pkin(3);
	t425 = t400 * t431;
	t424 = t401 * t397 * t406;
	t423 = qJD(3) * t442;
	t422 = t430 * t444;
	t421 = t405 * r_i_i_C(1) - t402 * r_i_i_C(2);
	t420 = -t402 * r_i_i_C(1) - t405 * r_i_i_C(2);
	t419 = qJ(4) - t420;
	t418 = pkin(4) + pkin(9) + t421;
	t373 = t387 * t403 + t388 * t437;
	t374 = -t403 * t415 - t414 * t437;
	t413 = t425 - t434;
	t412 = t400 * t432 + t433;
	t411 = t400 * t433 + t432;
	t410 = qJD(5) * t420;
	t409 = t421 * qJD(5) + qJD(4);
	t408 = t416 * t403 - t406 * t414;
	t386 = t401 * t400 - t407 * t444;
	t385 = t414 * qJD(2);
	t384 = t415 * qJD(2);
	t383 = t388 * qJD(2);
	t382 = -qJD(2) * t426 + t396 * t430;
	t381 = t412 * t398;
	t378 = t396 * t439 + t397 * t415;
	t377 = -t387 * t397 - t399 * t439;
	t375 = -t413 * t398 - t424;
	t371 = (-qJD(2) * t425 - qJD(3) * t431 + (qJD(3) * t400 + qJD(2)) * t434) * t398;
	t365 = t401 * t423 + (t412 * qJD(2) + t411 * qJD(3)) * t398;
	t363 = -t384 * t437 + t385 * t403 + (-t406 * t415 + t414 * t438) * qJD(3);
	t361 = -t382 * t437 - t383 * t403 + (t387 * t406 - t388 * t438) * qJD(3);
	t359 = t408 * qJD(3) - t384 * t403 - t385 * t437;
	t357 = -t382 * t403 + t383 * t437 - t423 * t440 + (t387 * t438 + t446) * qJD(3);
	t1 = [0, (t363 * t402 + t374 * t428) * r_i_i_C(1) + (t363 * t405 - t374 * t429) * r_i_i_C(2) + t363 * qJ(4) + t374 * qJD(4) + t385 * pkin(2) + t427 * (-t374 * qJD(3) + t384 * t438 + t385 * t406) + (-t418 * t384 - t410 * t414) * t397, t409 * t408 + t419 * (t448 * qJD(3) - t384 * t406 + t385 * t438) - t427 * t359, t359, (t359 * t405 + t385 * t443) * r_i_i_C(1) + (-t359 * t402 + t385 * t441) * r_i_i_C(2) + ((-t378 * t405 + t402 * t448) * r_i_i_C(1) + (t378 * t402 + t405 * t448) * r_i_i_C(2)) * qJD(5), 0; 0, (t361 * t402 + t373 * t428) * r_i_i_C(1) + (t361 * t405 - t373 * t429) * r_i_i_C(2) + t361 * qJ(4) + t373 * qJD(4) - t383 * pkin(2) + t427 * (-t373 * qJD(3) + t382 * t438 - t383 * t406) + (-t418 * t382 + t388 * t410) * t397, t409 * (-t417 * t403 + t446) + t419 * (-t367 * qJD(3) - t382 * t406 - t383 * t438) - t427 * t357, t357, (t357 * t405 - t383 * t443) * r_i_i_C(1) + (-t357 * t402 - t383 * t441) * r_i_i_C(2) + ((-t367 * t402 - t377 * t405) * r_i_i_C(1) + (-t367 * t405 + t377 * t402) * r_i_i_C(2)) * qJD(5), 0; 0, (-t371 * t402 + t381 * t428) * r_i_i_C(1) + (-t371 * t405 - t381 * t429) * r_i_i_C(2) - t371 * qJ(4) + t381 * qJD(4) + (-t427 * (t411 * qJD(2) + t412 * qJD(3)) - pkin(2) * t430 + (t418 * t407 * qJD(2) + t404 * t410) * t397) * t398, t409 * (t411 * t398 + t401 * t442) + t419 * (qJD(3) * t424 + (t413 * qJD(3) + (-t400 * t434 + t431) * qJD(2)) * t398) - t427 * t365, t365, (t365 * t405 - t402 * t422) * r_i_i_C(1) + (-t365 * t402 - t405 * t422) * r_i_i_C(2) + ((-t375 * t402 - t386 * t405) * r_i_i_C(1) + (-t375 * t405 + t386 * t402) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:15
	% EndTime: 2019-10-09 22:39:16
	% DurationCPUTime: 1.36s
	% Computational Cost: add. (1242->168), mult. (4014->298), div. (0->0), fcn. (4428->14), ass. (0->105)
	t568 = sin(pkin(12));
	t571 = cos(pkin(12));
	t577 = sin(qJ(2));
	t573 = cos(pkin(6));
	t581 = cos(qJ(2));
	t617 = t573 * t581;
	t556 = -t568 * t577 + t571 * t617;
	t634 = pkin(4) + pkin(9);
	t633 = r_i_i_C(3) + pkin(11);
	t618 = t573 * t577;
	t557 = t568 * t581 + t571 * t618;
	t576 = sin(qJ(3));
	t632 = t557 * t576;
	t580 = cos(qJ(3));
	t631 = t557 * t580;
	t589 = t568 * t618 - t571 * t581;
	t630 = t589 * t576;
	t569 = sin(pkin(7));
	t570 = sin(pkin(6));
	t628 = t569 * t570;
	t575 = sin(qJ(5));
	t627 = t569 * t575;
	t626 = t569 * t576;
	t579 = cos(qJ(5));
	t625 = t569 * t579;
	t624 = t569 * t580;
	t623 = t569 * t581;
	t622 = t570 * t571;
	t572 = cos(pkin(7));
	t621 = t570 * t572;
	t620 = t572 * t576;
	t619 = t572 * t580;
	t616 = t576 * t577;
	t615 = t576 * t581;
	t614 = t577 * t580;
	t613 = t580 * t581;
	t612 = qJD(2) * t570;
	t611 = qJD(5) * t575;
	t610 = qJD(5) * t579;
	t609 = t634 * t569;
	t608 = t577 * t628;
	t607 = t570 * t624;
	t606 = t570 * t616;
	t604 = t572 * t613;
	t603 = t573 * t624;
	t602 = t569 * t612;
	t601 = qJD(3) * t626;
	t600 = t570 * t604;
	t599 = t577 * t602;
	t598 = t581 * t602;
	t574 = sin(qJ(6));
	t578 = cos(qJ(6));
	t597 = -t578 * r_i_i_C(1) + t574 * r_i_i_C(2);
	t596 = -t574 * r_i_i_C(1) - t578 * r_i_i_C(2);
	t528 = -t556 * t619 + t571 * t607 + t632;
	t543 = -t556 * t569 - t571 * t621;
	t519 = t528 * t575 + t543 * t579;
	t590 = t568 * t617 + t571 * t577;
	t530 = -t568 * t607 + t590 * t619 - t630;
	t544 = t568 * t621 + t569 * t590;
	t521 = t530 * t575 + t544 * t579;
	t541 = -t600 - t603 + t606;
	t555 = -t570 * t623 + t573 * t572;
	t595 = t541 * t579 - t555 * t575;
	t533 = t541 * t575 + t555 * t579;
	t594 = pkin(5) - t597;
	t593 = pkin(3) + pkin(10) - t596;
	t536 = t556 * t576 + t557 * t619;
	t522 = t536 * t575 + t557 * t625;
	t538 = -t576 * t590 - t589 * t619;
	t523 = t538 * t575 - t589 * t625;
	t537 = t556 * t580 - t557 * t620;
	t592 = t556 * t572 - t569 * t622;
	t539 = -t580 * t590 + t589 * t620;
	t591 = t568 * t628 - t572 * t590;
	t588 = t572 * t614 + t615;
	t587 = t572 * t615 + t614;
	t586 = -t572 * t616 + t613;
	t585 = qJD(6) * t597;
	t584 = qJD(6) * t596;
	t549 = t588 * t570;
	t540 = t549 * t575 + t579 * t608;
	t531 = t591 * t576 - t580 * t589;
	t583 = t594 * t575 - t633 * t579 + qJ(4);
	t582 = qJD(4) + t575 * t584 + (t633 * t575 + t594 * t579) * qJD(5);
	t554 = t589 * qJD(2);
	t553 = t590 * qJD(2);
	t552 = t557 * qJD(2);
	t551 = t556 * qJD(2);
	t550 = t586 * t570;
	t542 = t587 * t570 + t573 * t626;
	t534 = -qJD(2) * t600 - t570 * qJD(3) * t613 + (qJD(3) * t572 + qJD(2)) * t606;
	t529 = t592 * t576 + t631;
	t527 = qJD(3) * t603 + ((t604 - t616) * qJD(3) + t586 * qJD(2)) * t570;
	t526 = t573 * t601 + (t588 * qJD(2) + t587 * qJD(3)) * t570;
	t516 = t539 * qJD(3) - t553 * t619 + t554 * t576;
	t514 = t537 * qJD(3) + t551 * t619 - t552 * t576;
	t511 = t554 * t620 - t553 * t580 + (t591 * t580 + t630) * qJD(3);
	t510 = t531 * qJD(3) - t553 * t576 - t554 * t619;
	t509 = -t552 * t620 + t551 * t580 + (t592 * t580 - t632) * qJD(3);
	t508 = t551 * t576 + t552 * t619 - t601 * t622 + (t556 * t620 + t631) * qJD(3);
	t505 = t595 * qJD(5) + t526 * t575 + t579 * t599;
	t499 = -t510 * t575 - t530 * t610 + t544 * t611 + t554 * t625;
	t497 = -t508 * t575 - t528 * t610 + t543 * t611 - t552 * t625;
	t1 = [0, t554 * pkin(2) + t516 * qJ(4) + t538 * qJD(4) - t553 * t609 + t594 * (-t553 * t625 + t516 * t575 + (t538 * t579 + t589 * t627) * qJD(5)) + t593 * (-t538 * qJD(3) + t553 * t620 + t554 * t580) - t633 * (-t523 * qJD(5) + t516 * t579 + t553 * t627) + ((-t523 * t574 + t539 * t578) * r_i_i_C(1) + (-t523 * t578 - t539 * t574) * r_i_i_C(2)) * qJD(6), -t593 * t510 + t583 * t511 + t530 * t585 + t582 * t531, t510, -t633 * t499 + (t530 * t579 - t544 * t575) * t584 + t594 * (-t521 * qJD(5) + t510 * t579 + t554 * t627), (t499 * t574 + t511 * t578) * r_i_i_C(1) + (t499 * t578 - t511 * t574) * r_i_i_C(2) + ((-t521 * t578 - t531 * t574) * r_i_i_C(1) + (t521 * t574 - t531 * t578) * r_i_i_C(2)) * qJD(6); 0, -t552 * pkin(2) + t514 * qJ(4) + t536 * qJD(4) + t551 * t609 + t594 * (t551 * t625 + t514 * t575 + (t536 * t579 - t557 * t627) * qJD(5)) + t593 * (-t536 * qJD(3) - t551 * t620 - t552 * t580) - t633 * (-t522 * qJD(5) + t514 * t579 - t551 * t627) + ((-t522 * t574 + t537 * t578) * r_i_i_C(1) + (-t522 * t578 - t537 * t574) * r_i_i_C(2)) * qJD(6), -t593 * t508 + t583 * t509 + t528 * t585 + t582 * t529, t508, -t633 * t497 + (t528 * t579 - t543 * t575) * t584 + t594 * (-t519 * qJD(5) + t508 * t579 - t552 * t627), (t497 * t574 + t509 * t578) * r_i_i_C(1) + (t497 * t578 - t509 * t574) * r_i_i_C(2) + ((-t519 * t578 - t529 * t574) * r_i_i_C(1) + (t519 * t574 - t529 * t578) * r_i_i_C(2)) * qJD(6); 0, -t534 * qJ(4) + t549 * qJD(4) + t594 * (t579 * t598 - t534 * t575 + (t549 * t579 - t575 * t608) * qJD(5)) - t593 * (t587 * qJD(2) + t588 * qJD(3)) * t570 + t633 * (t540 * qJD(5) + t534 * t579 + t575 * t598) + ((-t540 * t574 + t550 * t578) * r_i_i_C(1) + (-t540 * t578 - t550 * t574) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t577 + t634 * t623) * t612, -t593 * t526 + t583 * t527 + t541 * t585 + t582 * t542, t526, t633 * t505 + t595 * t584 + t594 * (-t533 * qJD(5) + t526 * t579 - t575 * t599), (-t505 * t574 + t527 * t578) * r_i_i_C(1) + (-t505 * t578 - t527 * t574) * r_i_i_C(2) + ((-t533 * t578 - t542 * t574) * r_i_i_C(1) + (t533 * t574 - t542 * t578) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end