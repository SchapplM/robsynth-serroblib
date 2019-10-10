% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
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
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:21
	% DurationCPUTime: 0.36s
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
	% StartTime: 2019-10-09 22:35:21
	% EndTime: 2019-10-09 22:35:22
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (272->74), mult. (919->136), div. (0->0), fcn. (944->12), ass. (0->52)
	t354 = sin(pkin(12));
	t358 = cos(pkin(12));
	t362 = sin(qJ(2));
	t360 = cos(pkin(6));
	t364 = cos(qJ(2));
	t383 = t360 * t364;
	t346 = -t354 * t362 + t358 * t383;
	t361 = sin(qJ(3));
	t363 = cos(qJ(3));
	t384 = t360 * t362;
	t369 = t354 * t384 - t358 * t364;
	t359 = cos(pkin(7));
	t370 = t354 * t383 + t358 * t362;
	t355 = sin(pkin(7));
	t356 = sin(pkin(6));
	t388 = t355 * t356;
	t371 = t354 * t388 - t359 * t370;
	t395 = t371 * t361 - t363 * t369;
	t394 = pkin(9) * t355;
	t393 = r_i_i_C(3) + qJ(4);
	t347 = t354 * t364 + t358 * t384;
	t392 = t347 * t363;
	t353 = sin(pkin(13));
	t390 = t353 * t355;
	t357 = cos(pkin(13));
	t387 = t355 * t357;
	t386 = t359 * t361;
	t385 = t359 * t363;
	t382 = t361 * t362;
	t381 = t361 * t364;
	t380 = t362 * t363;
	t379 = t363 * t364;
	t378 = qJD(3) * t355;
	t376 = t361 * t378;
	t375 = r_i_i_C(1) * t357 - r_i_i_C(2) * t353 + pkin(3);
	t374 = -t346 * t361 - t347 * t385;
	t373 = -t346 * t359 + t358 * t388;
	t372 = t361 * t370 + t369 * t385;
	t368 = t359 * t379 - t382;
	t367 = t359 * t380 + t381;
	t366 = t359 * t381 + t380;
	t365 = t359 * t382 - t379;
	t345 = t369 * qJD(2);
	t344 = t370 * qJD(2);
	t343 = t347 * qJD(2);
	t342 = t346 * qJD(2);
	t336 = t360 * t376 + (t367 * qJD(2) + t366 * qJD(3)) * t356;
	t335 = t372 * qJD(3) + t344 * t386 + t345 * t363;
	t333 = t374 * qJD(3) - t342 * t386 - t343 * t363;
	t330 = t395 * qJD(3) - t344 * t361 - t345 * t385;
	t328 = t342 * t361 + t343 * t385 - t358 * t356 * t376 + (t346 * t386 + t392) * qJD(3);
	t1 = [0, (t335 * t357 - t344 * t390) * r_i_i_C(1) + (-t335 * t353 - t344 * t387) * r_i_i_C(2) + t335 * pkin(3) - t372 * qJD(4) + t345 * pkin(2) - t344 * t394 + t393 * (-t344 * t385 + t345 * t361 + (-t363 * t370 + t369 * t386) * qJD(3)), t395 * qJD(4) + t393 * (t345 * t386 - t344 * t363 + (t361 * t369 + t371 * t363) * qJD(3)) - t375 * t330, t330, 0, 0; 0, (t333 * t357 + t342 * t390) * r_i_i_C(1) + (-t333 * t353 + t342 * t387) * r_i_i_C(2) + t333 * pkin(3) - t374 * qJD(4) - t343 * pkin(2) + t342 * t394 + t393 * (t342 * t385 - t343 * t361 + (t346 * t363 - t347 * t386) * qJD(3)), -(t373 * t361 - t392) * qJD(4) + t393 * (-t343 * t386 + t342 * t363 + (-t347 * t361 - t373 * t363) * qJD(3)) - t375 * t328, t328, 0, 0; 0, (-t393 * (-t368 * qJD(2) + t365 * qJD(3)) + t375 * (-t366 * qJD(2) - t367 * qJD(3)) + t367 * qJD(4) + (-t362 * pkin(2) + (r_i_i_C(1) * t353 + r_i_i_C(2) * t357 + pkin(9)) * t364 * t355) * qJD(2)) * t356, -(-t360 * t355 * t361 - t366 * t356) * qJD(4) + t393 * (t360 * t363 * t378 + (-t365 * qJD(2) + t368 * qJD(3)) * t356) - t375 * t336, t336, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:22
	% EndTime: 2019-10-09 22:35:22
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (500->113), mult. (1426->206), div. (0->0), fcn. (1512->14), ass. (0->80)
	t461 = r_i_i_C(3) + pkin(10) + qJ(4);
	t410 = sin(pkin(12));
	t413 = cos(pkin(12));
	t420 = cos(qJ(2));
	t415 = cos(pkin(6));
	t418 = sin(qJ(2));
	t451 = t415 * t418;
	t399 = t410 * t420 + t413 * t451;
	t419 = cos(qJ(3));
	t460 = t399 * t419;
	t408 = pkin(13) + qJ(5);
	t406 = sin(t408);
	t411 = sin(pkin(7));
	t459 = t406 * t411;
	t407 = cos(t408);
	t458 = t407 * t411;
	t412 = sin(pkin(6));
	t457 = t411 * t412;
	t417 = sin(qJ(3));
	t456 = t411 * t417;
	t455 = t412 * t413;
	t414 = cos(pkin(7));
	t454 = t412 * t414;
	t453 = t414 * t417;
	t452 = t414 * t419;
	t450 = t415 * t420;
	t449 = t417 * t418;
	t448 = t417 * t420;
	t447 = t418 * t419;
	t446 = t419 * t420;
	t445 = qJD(2) * t418;
	t444 = qJD(5) * t406;
	t443 = qJD(5) * t407;
	t442 = t413 * t450;
	t441 = t415 * t411 * t419;
	t440 = qJD(3) * t456;
	t439 = t445 * t457;
	t438 = t407 * r_i_i_C(1) - t406 * r_i_i_C(2);
	t437 = -r_i_i_C(1) * t406 - r_i_i_C(2) * t407;
	t405 = cos(pkin(13)) * pkin(4) + pkin(3);
	t436 = -t405 - t438;
	t398 = -t410 * t418 + t442;
	t435 = -t398 * t417 - t399 * t452;
	t384 = t398 * t419 - t399 * t453;
	t434 = t398 * t414 - t411 * t455;
	t430 = t410 * t451 - t413 * t420;
	t431 = t410 * t450 + t413 * t418;
	t433 = t417 * t431 + t430 * t452;
	t385 = -t419 * t431 + t430 * t453;
	t432 = t410 * t457 - t414 * t431;
	t429 = t414 * t446 - t449;
	t428 = t414 * t447 + t448;
	t427 = t414 * t448 + t447;
	t426 = t414 * t449 - t446;
	t425 = qJD(5) * t438;
	t424 = qJD(5) * t437;
	t423 = sin(pkin(13)) * pkin(4) + pkin(9) - t437;
	t422 = -t399 * t417 + t434 * t419;
	t421 = t417 * t430 + t432 * t419;
	t381 = t432 * t417 - t419 * t430;
	t397 = t415 * t414 - t420 * t457;
	t396 = t430 * qJD(2);
	t395 = t431 * qJD(2);
	t394 = t399 * qJD(2);
	t393 = -qJD(2) * t442 + t410 * t445;
	t392 = t426 * t412;
	t389 = t410 * t454 + t411 * t431;
	t388 = -t398 * t411 - t413 * t454;
	t387 = t427 * t412 + t415 * t456;
	t383 = (-t427 * qJD(2) - t428 * qJD(3)) * t412;
	t379 = t434 * t417 + t460;
	t377 = qJD(3) * t441 + (-t426 * qJD(2) + t429 * qJD(3)) * t412;
	t376 = t415 * t440 + (t428 * qJD(2) + t427 * qJD(3)) * t412;
	t375 = t433 * qJD(3) + t395 * t453 + t396 * t419;
	t373 = t435 * qJD(3) + t393 * t453 - t394 * t419;
	t371 = t421 * qJD(3) - t395 * t419 + t396 * t453;
	t370 = t381 * qJD(3) - t395 * t417 - t396 * t452;
	t369 = t422 * qJD(3) - t393 * t419 - t394 * t453;
	t368 = -t393 * t417 + t394 * t452 - t440 * t455 + (t398 * t453 + t460) * qJD(3);
	t1 = [0, (t375 * t407 - t385 * t444) * r_i_i_C(1) + (-t375 * t406 - t385 * t443) * r_i_i_C(2) + t375 * t405 - t433 * qJD(4) + t396 * pkin(2) + t461 * (t385 * qJD(3) - t395 * t452 + t396 * t417) + (-t423 * t395 - t425 * t430) * t411, qJD(4) * t381 + t436 * t370 + t461 * t371 + t421 * t424, t370, (-t371 * t406 - t396 * t458) * r_i_i_C(1) + (-t371 * t407 + t396 * t459) * r_i_i_C(2) + ((-t381 * t407 - t389 * t406) * r_i_i_C(1) + (t381 * t406 - t389 * t407) * r_i_i_C(2)) * qJD(5), 0; 0, (t373 * t407 - t384 * t444) * r_i_i_C(1) + (-t373 * t406 - t384 * t443) * r_i_i_C(2) + t373 * t405 - t435 * qJD(4) - t394 * pkin(2) + t461 * (t384 * qJD(3) - t393 * t452 - t394 * t417) + (-t423 * t393 + t399 * t425) * t411, qJD(4) * t379 + t436 * t368 + t461 * t369 + t422 * t424, t368, (-t369 * t406 + t394 * t458) * r_i_i_C(1) + (-t369 * t407 - t394 * t459) * r_i_i_C(2) + ((-t379 * t407 - t388 * t406) * r_i_i_C(1) + (t379 * t406 - t388 * t407) * r_i_i_C(2)) * qJD(5), 0; 0, (t383 * t407 + t392 * t444) * r_i_i_C(1) + (-t383 * t406 + t392 * t443) * r_i_i_C(2) + t383 * t405 + (-t461 * (-t429 * qJD(2) + t426 * qJD(3)) + t428 * qJD(4) - pkin(2) * t445 + (t423 * t420 * qJD(2) + t418 * t425) * t411) * t412, qJD(4) * t387 + t461 * t377 + (t429 * t412 + t441) * t424 + t436 * t376, t376, (-t377 * t406 + t407 * t439) * r_i_i_C(1) + (-t377 * t407 - t406 * t439) * r_i_i_C(2) + ((-t387 * t407 - t397 * t406) * r_i_i_C(1) + (t387 * t406 - t397 * t407) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:25
	% EndTime: 2019-10-09 22:35:26
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (1427->189), mult. (3882->327), div. (0->0), fcn. (4294->16), ass. (0->107)
	t569 = sin(qJ(2));
	t571 = cos(qJ(2));
	t564 = cos(pkin(12));
	t624 = cos(pkin(6));
	t601 = t564 * t624;
	t623 = sin(pkin(12));
	t545 = -t623 * t569 + t571 * t601;
	t567 = sin(qJ(6));
	t570 = cos(qJ(6));
	t584 = (r_i_i_C(1) * t567 + r_i_i_C(2) * t570) * qJD(6);
	t626 = r_i_i_C(3) + pkin(11);
	t625 = cos(qJ(3));
	t560 = pkin(13) + qJ(5);
	t558 = sin(t560);
	t562 = sin(pkin(7));
	t622 = t558 * t562;
	t559 = cos(t560);
	t621 = t559 * t562;
	t620 = t562 * t571;
	t563 = sin(pkin(6));
	t619 = t563 * t564;
	t565 = cos(pkin(7));
	t568 = sin(qJ(3));
	t618 = t565 * t568;
	t617 = t568 * t569;
	t616 = t568 * t571;
	t615 = qJD(2) * t563;
	t614 = qJD(6) * t567;
	t613 = qJD(6) * t570;
	t612 = t562 * t619;
	t611 = t562 * t563 * t569;
	t610 = t563 * t617;
	t609 = pkin(4) * sin(pkin(13)) + pkin(9);
	t578 = -t569 * t601 - t623 * t571;
	t608 = t578 * t625;
	t607 = t565 * t625;
	t606 = t625 * t569;
	t605 = t625 * t571;
	t604 = t562 * t615;
	t603 = t562 * t624;
	t602 = t563 * t623;
	t599 = t565 * t605;
	t598 = t569 * t604;
	t597 = t571 * t604;
	t595 = t568 * t603;
	t594 = t562 * t602;
	t593 = t562 * t609;
	t591 = t624 * t623;
	t590 = t563 * t599;
	t520 = -t608 + (t545 * t565 - t612) * t568;
	t532 = -t545 * t562 - t565 * t619;
	t506 = t520 * t559 + t532 * t558;
	t589 = -t520 * t558 + t532 * t559;
	t576 = t564 * t569 + t571 * t591;
	t577 = -t564 * t571 + t569 * t591;
	t522 = -t577 * t625 + (-t565 * t576 + t594) * t568;
	t533 = t562 * t576 + t565 * t602;
	t508 = t522 * t559 + t533 * t558;
	t588 = -t522 * t558 + t533 * t559;
	t580 = t565 * t616 + t606;
	t531 = t580 * t563 + t595;
	t544 = -t563 * t620 + t624 * t565;
	t516 = t531 * t559 + t544 * t558;
	t587 = -t531 * t558 + t544 * t559;
	t586 = r_i_i_C(1) * t570 - r_i_i_C(2) * t567 + pkin(5);
	t585 = t625 * t603;
	t526 = t545 * t625 + t578 * t618;
	t513 = t526 * t559 - t578 * t622;
	t528 = -t576 * t625 + t577 * t618;
	t514 = t528 * t559 - t577 * t622;
	t581 = -t565 * t617 + t605;
	t539 = t581 * t563;
	t529 = t539 * t559 + t558 * t611;
	t583 = -t545 * t568 + t578 * t607;
	t582 = t568 * t576 + t577 * t607;
	t579 = t565 * t606 + t616;
	t557 = cos(pkin(13)) * pkin(4) + pkin(3);
	t575 = -t626 * t558 - t586 * t559 - t557;
	t574 = t545 * t607 + t568 * t578 - t625 * t612;
	t573 = t568 * t577 - t576 * t607 + t625 * t594;
	t572 = t559 * t584 + (t586 * t558 - t626 * t559) * qJD(5);
	t566 = -pkin(10) - qJ(4);
	t543 = t577 * qJD(2);
	t542 = t576 * qJD(2);
	t541 = t578 * qJD(2);
	t540 = t545 * qJD(2);
	t538 = t579 * t563;
	t530 = -t585 - t590 + t610;
	t524 = (-t580 * qJD(2) - t579 * qJD(3)) * t563;
	t523 = -qJD(2) * t590 - t563 * qJD(3) * t605 + (qJD(3) * t565 + qJD(2)) * t610;
	t518 = qJD(3) * t585 + ((t599 - t617) * qJD(3) + t581 * qJD(2)) * t563;
	t517 = qJD(3) * t595 + (t579 * qJD(2) + t580 * qJD(3)) * t563;
	t512 = t582 * qJD(3) + t542 * t618 + t543 * t625;
	t511 = t528 * qJD(3) - t542 * t607 + t543 * t568;
	t510 = t583 * qJD(3) - t540 * t618 + t541 * t625;
	t509 = t526 * qJD(3) + t540 * t607 + t541 * t568;
	t504 = t573 * qJD(3) - t542 * t625 + t543 * t618;
	t503 = t522 * qJD(3) - t542 * t568 - t543 * t607;
	t502 = t574 * qJD(3) + t540 * t625 + t541 * t618;
	t501 = t540 * t568 - t541 * t607 + (t545 * t618 - t568 * t612 - t608) * qJD(3);
	t500 = t558 * t597 + t524 * t559 + (-t539 * t558 + t559 * t611) * qJD(5);
	t498 = t587 * qJD(5) + t518 * t559 + t558 * t598;
	t496 = -t542 * t622 + t512 * t559 + (-t528 * t558 - t577 * t621) * qJD(5);
	t494 = t540 * t622 + t510 * t559 + (-t526 * t558 - t578 * t621) * qJD(5);
	t492 = t588 * qJD(5) + t504 * t559 - t543 * t622;
	t490 = t589 * qJD(5) + t502 * t559 - t541 * t622;
	t1 = [0, (t496 * t570 + t511 * t567) * r_i_i_C(1) + (-t496 * t567 + t511 * t570) * r_i_i_C(2) + t496 * pkin(5) + t512 * t557 - t511 * t566 - t582 * qJD(4) + t543 * pkin(2) - t542 * t593 + t626 * (t514 * qJD(5) + t512 * t558 + t542 * t621) + ((-t514 * t567 - t570 * t582) * r_i_i_C(1) + (-t514 * t570 + t567 * t582) * r_i_i_C(2)) * qJD(6), (t504 * t567 + t522 * t613) * r_i_i_C(1) + (t504 * t570 - t522 * t614) * r_i_i_C(2) - t504 * t566 + t522 * qJD(4) + t575 * t503 - t572 * t573, t503, t626 * t492 - t588 * t584 + t586 * (-t508 * qJD(5) - t504 * t558 - t543 * t621), (-t492 * t567 + t503 * t570) * r_i_i_C(1) + (-t492 * t570 - t503 * t567) * r_i_i_C(2) + ((-t508 * t570 + t567 * t573) * r_i_i_C(1) + (t508 * t567 + t570 * t573) * r_i_i_C(2)) * qJD(6); 0, (t494 * t570 + t509 * t567) * r_i_i_C(1) + (-t494 * t567 + t509 * t570) * r_i_i_C(2) + t494 * pkin(5) + t510 * t557 - t509 * t566 - t583 * qJD(4) + t541 * pkin(2) + t540 * t593 + t626 * (t513 * qJD(5) + t510 * t558 - t540 * t621) + ((-t513 * t567 - t570 * t583) * r_i_i_C(1) + (-t513 * t570 + t567 * t583) * r_i_i_C(2)) * qJD(6), (t502 * t567 + t520 * t613) * r_i_i_C(1) + (t502 * t570 - t520 * t614) * r_i_i_C(2) - t502 * t566 + t520 * qJD(4) + t575 * t501 - t572 * t574, t501, t626 * t490 - t589 * t584 + t586 * (-t506 * qJD(5) - t502 * t558 - t541 * t621), (-t490 * t567 + t501 * t570) * r_i_i_C(1) + (-t490 * t570 - t501 * t567) * r_i_i_C(2) + ((-t506 * t570 + t567 * t574) * r_i_i_C(1) + (t506 * t567 + t570 * t574) * r_i_i_C(2)) * qJD(6); 0, (t500 * t570 - t523 * t567) * r_i_i_C(1) + (-t500 * t567 - t523 * t570) * r_i_i_C(2) + t500 * pkin(5) + t524 * t557 + t523 * t566 + t538 * qJD(4) + t626 * (t529 * qJD(5) + t524 * t558 - t559 * t597) + ((-t529 * t567 + t538 * t570) * r_i_i_C(1) + (-t529 * t570 - t538 * t567) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t569 + t609 * t620) * t615, (t518 * t567 + t531 * t613) * r_i_i_C(1) + (t518 * t570 - t531 * t614) * r_i_i_C(2) - t518 * t566 + t531 * qJD(4) + t575 * t517 + t572 * t530, t517, t626 * t498 - t587 * t584 + t586 * (-t516 * qJD(5) - t518 * t558 + t559 * t598), (-t498 * t567 + t517 * t570) * r_i_i_C(1) + (-t498 * t570 - t517 * t567) * r_i_i_C(2) + ((-t516 * t570 - t530 * t567) * r_i_i_C(1) + (t516 * t567 - t530 * t570) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end