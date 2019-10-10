% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
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
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:09
	% DurationCPUTime: 0.35s
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
	% StartTime: 2019-10-09 22:54:10
	% EndTime: 2019-10-09 22:54:11
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (380->110), mult. (1266->209), div. (0->0), fcn. (1334->12), ass. (0->67)
	t400 = sin(qJ(3));
	t403 = cos(qJ(3));
	t393 = sin(pkin(12));
	t396 = cos(pkin(12));
	t404 = cos(qJ(2));
	t398 = cos(pkin(6));
	t401 = sin(qJ(2));
	t429 = t398 * t401;
	t412 = t393 * t429 - t396 * t404;
	t397 = cos(pkin(7));
	t428 = t398 * t404;
	t413 = t393 * t428 + t396 * t401;
	t394 = sin(pkin(7));
	t395 = sin(pkin(6));
	t436 = t394 * t395;
	t414 = t393 * t436 - t397 * t413;
	t372 = t414 * t400 - t403 * t412;
	t388 = t393 * t404 + t396 * t429;
	t419 = t396 * t428;
	t387 = -t393 * t401 + t419;
	t415 = -t387 * t397 + t396 * t436;
	t440 = -t388 * t403 + t415 * t400;
	t439 = r_i_i_C(3) + pkin(10);
	t435 = t394 * t398;
	t399 = sin(qJ(4));
	t434 = t394 * t399;
	t402 = cos(qJ(4));
	t433 = t394 * t402;
	t432 = t395 * t397;
	t431 = t397 * t400;
	t430 = t397 * t403;
	t427 = t400 * t401;
	t426 = t400 * t404;
	t425 = t401 * t403;
	t424 = t403 * t404;
	t423 = qJD(2) * t401;
	t422 = qJD(2) * t404;
	t421 = qJD(4) * t399;
	t420 = qJD(4) * t402;
	t418 = qJD(3) * t435;
	t417 = t423 * t436;
	t416 = t402 * r_i_i_C(1) - t399 * r_i_i_C(2) + pkin(3);
	t375 = t387 * t403 - t388 * t431;
	t376 = -t403 * t413 + t412 * t431;
	t411 = t397 * t424 - t427;
	t410 = -t397 * t425 - t426;
	t409 = t397 * t426 + t425;
	t408 = t397 * t427 - t424;
	t407 = qJD(4) * (-t399 * r_i_i_C(1) - t402 * r_i_i_C(2));
	t406 = -t388 * t400 - t415 * t403;
	t405 = t400 * t412 + t414 * t403;
	t386 = t398 * t397 - t404 * t436;
	t385 = t412 * qJD(2);
	t384 = t413 * qJD(2);
	t383 = t388 * qJD(2);
	t382 = -qJD(2) * t419 + t393 * t423;
	t381 = t408 * t395;
	t380 = t393 * t432 + t394 * t413;
	t379 = -t387 * t394 - t396 * t432;
	t378 = t409 * t395 + t400 * t435;
	t374 = (-t409 * qJD(2) + t410 * qJD(3)) * t395;
	t368 = t403 * t418 + (-t408 * qJD(2) + t411 * qJD(3)) * t395;
	t366 = t384 * t431 + t385 * t403 + (t400 * t413 + t412 * t430) * qJD(3);
	t364 = t382 * t431 - t383 * t403 + (-t387 * t400 - t388 * t430) * qJD(3);
	t362 = t405 * qJD(3) - t384 * t403 + t385 * t431;
	t360 = t406 * qJD(3) - t382 * t403 - t383 * t431;
	t1 = [0, (t366 * t402 - t376 * t421) * r_i_i_C(1) + (-t366 * t399 - t376 * t420) * r_i_i_C(2) + t366 * pkin(3) + t385 * pkin(2) + t439 * (t376 * qJD(3) - t384 * t430 + t385 * t400) + ((-t384 * t399 - t412 * t420) * r_i_i_C(1) + (-t384 * t402 + t412 * t421) * r_i_i_C(2) - t384 * pkin(9)) * t394, t439 * t362 + t405 * t407 + t416 * (-t372 * qJD(3) + t384 * t400 + t385 * t430), (-t362 * t399 - t385 * t433) * r_i_i_C(1) + (-t362 * t402 + t385 * t434) * r_i_i_C(2) + ((-t372 * t402 - t380 * t399) * r_i_i_C(1) + (t372 * t399 - t380 * t402) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t364 * t402 - t375 * t421) * r_i_i_C(1) + (-t364 * t399 - t375 * t420) * r_i_i_C(2) + t364 * pkin(3) - t383 * pkin(2) + t439 * (t375 * qJD(3) - t382 * t430 - t383 * t400) + ((-t382 * t399 + t388 * t420) * r_i_i_C(1) + (-t382 * t402 - t388 * t421) * r_i_i_C(2) - t382 * pkin(9)) * t394, t439 * t360 + t406 * t407 + t416 * (t440 * qJD(3) + t382 * t400 - t383 * t430), (-t360 * t399 + t383 * t433) * r_i_i_C(1) + (-t360 * t402 - t383 * t434) * r_i_i_C(2) + ((-t379 * t399 + t402 * t440) * r_i_i_C(1) + (-t379 * t402 - t399 * t440) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t374 * t402 + t381 * t421) * r_i_i_C(1) + (-t374 * t399 + t381 * t420) * r_i_i_C(2) + t374 * pkin(3) + (-t439 * (-t411 * qJD(2) + t408 * qJD(3)) - pkin(2) * t423 + ((t399 * t422 + t401 * t420) * r_i_i_C(1) + (-t401 * t421 + t402 * t422) * r_i_i_C(2) + pkin(9) * t422) * t394) * t395, t439 * t368 + (t411 * t395 + t403 * t435) * t407 + t416 * (-t400 * t418 + (t410 * qJD(2) - t409 * qJD(3)) * t395), (-t368 * t399 + t402 * t417) * r_i_i_C(1) + (-t368 * t402 - t399 * t417) * r_i_i_C(2) + ((-t378 * t402 - t386 * t399) * r_i_i_C(1) + (t378 * t399 - t386 * t402) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:10
	% EndTime: 2019-10-09 22:54:11
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (577->142), mult. (1681->254), div. (0->0), fcn. (1784->14), ass. (0->84)
	t419 = sin(qJ(4));
	t467 = t419 * pkin(4);
	t466 = r_i_i_C(3) + qJ(5) + pkin(10);
	t412 = sin(pkin(12));
	t415 = cos(pkin(12));
	t424 = cos(qJ(2));
	t417 = cos(pkin(6));
	t421 = sin(qJ(2));
	t455 = t417 * t421;
	t402 = t412 * t424 + t415 * t455;
	t423 = cos(qJ(3));
	t465 = t402 * t423;
	t411 = qJ(4) + pkin(13);
	t409 = sin(t411);
	t413 = sin(pkin(7));
	t464 = t409 * t413;
	t410 = cos(t411);
	t463 = t410 * t413;
	t414 = sin(pkin(6));
	t462 = t413 * t414;
	t420 = sin(qJ(3));
	t461 = t413 * t420;
	t422 = cos(qJ(4));
	t460 = t413 * t422;
	t459 = t414 * t415;
	t416 = cos(pkin(7));
	t458 = t414 * t416;
	t457 = t416 * t420;
	t456 = t416 * t423;
	t454 = t417 * t424;
	t453 = t420 * t421;
	t452 = t420 * t424;
	t451 = t421 * t423;
	t450 = t423 * t424;
	t449 = qJD(2) * t421;
	t448 = qJD(2) * t424;
	t447 = qJD(4) * t409;
	t446 = qJD(4) * t410;
	t445 = qJD(4) * t421;
	t444 = qJD(4) * t422;
	t443 = qJD(4) * t467;
	t442 = t415 * t454;
	t441 = t417 * t413 * t423;
	t440 = qJD(3) * t461;
	t439 = t449 * t462;
	t408 = t422 * pkin(4) + pkin(3);
	t438 = -t410 * r_i_i_C(1) + t409 * r_i_i_C(2) - t408;
	t401 = -t412 * t421 + t442;
	t437 = -t401 * t420 - t402 * t456;
	t387 = t401 * t423 - t402 * t457;
	t436 = t401 * t416 - t413 * t459;
	t432 = t412 * t455 - t415 * t424;
	t433 = t412 * t454 + t415 * t421;
	t435 = t420 * t433 + t432 * t456;
	t388 = -t423 * t433 + t432 * t457;
	t434 = t412 * t462 - t416 * t433;
	t431 = t416 * t450 - t453;
	t430 = t416 * t451 + t452;
	t429 = t416 * t452 + t451;
	t428 = t416 * t453 - t450;
	t427 = qJD(4) * (-t409 * r_i_i_C(1) - t410 * r_i_i_C(2) - t467);
	t426 = -t402 * t420 + t436 * t423;
	t425 = t420 * t432 + t434 * t423;
	t384 = t434 * t420 - t423 * t432;
	t400 = t417 * t416 - t424 * t462;
	t399 = t432 * qJD(2);
	t398 = t433 * qJD(2);
	t397 = t402 * qJD(2);
	t396 = -qJD(2) * t442 + t412 * t449;
	t395 = t428 * t414;
	t392 = t412 * t458 + t413 * t433;
	t391 = -t401 * t413 - t415 * t458;
	t390 = t429 * t414 + t417 * t461;
	t386 = (-t429 * qJD(2) - t430 * qJD(3)) * t414;
	t382 = t436 * t420 + t465;
	t380 = qJD(3) * t441 + (-t428 * qJD(2) + t431 * qJD(3)) * t414;
	t379 = t417 * t440 + (t430 * qJD(2) + t429 * qJD(3)) * t414;
	t378 = t435 * qJD(3) + t398 * t457 + t399 * t423;
	t376 = t437 * qJD(3) + t396 * t457 - t397 * t423;
	t374 = t425 * qJD(3) - t398 * t423 + t399 * t457;
	t373 = t384 * qJD(3) - t398 * t420 - t399 * t456;
	t372 = t426 * qJD(3) - t396 * t423 - t397 * t457;
	t371 = -t396 * t420 + t397 * t456 - t440 * t459 + (t401 * t457 + t465) * qJD(3);
	t1 = [0, (t378 * t410 - t388 * t447) * r_i_i_C(1) + (-t378 * t409 - t388 * t446) * r_i_i_C(2) + t378 * t408 - t388 * t443 - t435 * qJD(5) + t399 * pkin(2) + t466 * (t388 * qJD(3) - t398 * t456 + t399 * t420) + ((-t398 * t409 - t432 * t446) * r_i_i_C(1) + (-t398 * t410 + t432 * t447) * r_i_i_C(2) - t398 * pkin(9) + (-t398 * t419 - t432 * t444) * pkin(4)) * t413, t384 * qJD(5) + t438 * t373 + t466 * t374 + t425 * t427, (-t374 * t409 - t399 * t463) * r_i_i_C(1) + (-t374 * t410 + t399 * t464) * r_i_i_C(2) + ((-t384 * t410 - t392 * t409) * r_i_i_C(1) + (t384 * t409 - t392 * t410) * r_i_i_C(2)) * qJD(4) + (-t399 * t460 - t374 * t419 + (-t384 * t422 - t392 * t419) * qJD(4)) * pkin(4), t373, 0; 0, (t376 * t410 - t387 * t447) * r_i_i_C(1) + (-t376 * t409 - t387 * t446) * r_i_i_C(2) + t376 * t408 - t387 * t443 - t437 * qJD(5) - t397 * pkin(2) + t466 * (t387 * qJD(3) - t396 * t456 - t397 * t420) + ((-t396 * t409 + t402 * t446) * r_i_i_C(1) + (-t396 * t410 - t402 * t447) * r_i_i_C(2) - t396 * pkin(9) + (-t396 * t419 + t402 * t444) * pkin(4)) * t413, t382 * qJD(5) + t438 * t371 + t466 * t372 + t426 * t427, (-t372 * t409 + t397 * t463) * r_i_i_C(1) + (-t372 * t410 - t397 * t464) * r_i_i_C(2) + ((-t382 * t410 - t391 * t409) * r_i_i_C(1) + (t382 * t409 - t391 * t410) * r_i_i_C(2)) * qJD(4) + (t397 * t460 - t372 * t419 + (-t382 * t422 - t391 * t419) * qJD(4)) * pkin(4), t371, 0; 0, (t386 * t410 + t395 * t447) * r_i_i_C(1) + (-t386 * t409 + t395 * t446) * r_i_i_C(2) + t386 * t408 + t395 * t443 + (-t466 * (-t431 * qJD(2) + t428 * qJD(3)) + t430 * qJD(5) - pkin(2) * t449 + ((t409 * t448 + t410 * t445) * r_i_i_C(1) + (-t409 * t445 + t410 * t448) * r_i_i_C(2) + pkin(9) * t448 + (t419 * t448 + t421 * t444) * pkin(4)) * t413) * t414, t390 * qJD(5) + t466 * t380 + t438 * t379 + (t431 * t414 + t441) * t427, (-t380 * t409 + t410 * t439) * r_i_i_C(1) + (-t380 * t410 - t409 * t439) * r_i_i_C(2) + ((-t390 * t410 - t400 * t409) * r_i_i_C(1) + (t390 * t409 - t400 * t410) * r_i_i_C(2)) * qJD(4) + (t422 * t439 - t380 * t419 + (-t390 * t422 - t400 * t419) * qJD(4)) * pkin(4), t379, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:13
	% EndTime: 2019-10-09 22:54:15
	% DurationCPUTime: 1.69s
	% Computational Cost: add. (1504->210), mult. (4137->363), div. (0->0), fcn. (4566->16), ass. (0->112)
	t570 = sin(qJ(2));
	t573 = cos(qJ(2));
	t564 = cos(pkin(12));
	t626 = cos(pkin(6));
	t602 = t564 * t626;
	t625 = sin(pkin(12));
	t546 = -t625 * t570 + t573 * t602;
	t567 = sin(qJ(6));
	t571 = cos(qJ(6));
	t586 = (r_i_i_C(1) * t567 + r_i_i_C(2) * t571) * qJD(6);
	t631 = r_i_i_C(3) + pkin(11);
	t630 = cos(qJ(3));
	t562 = sin(pkin(7));
	t629 = pkin(9) * t562;
	t568 = sin(qJ(4));
	t628 = t568 * pkin(4);
	t627 = pkin(4) * qJD(4);
	t561 = qJ(4) + pkin(13);
	t559 = sin(t561);
	t624 = t559 * t562;
	t560 = cos(t561);
	t623 = t560 * t562;
	t563 = sin(pkin(6));
	t622 = t562 * t563;
	t621 = t562 * t568;
	t572 = cos(qJ(4));
	t620 = t562 * t572;
	t619 = t562 * t573;
	t618 = t563 * t564;
	t565 = cos(pkin(7));
	t569 = sin(qJ(3));
	t617 = t565 * t569;
	t616 = t569 * t570;
	t615 = t569 * t573;
	t614 = qJD(6) * t567;
	t613 = qJD(6) * t571;
	t612 = t562 * t618;
	t611 = t570 * t622;
	t610 = t563 * t616;
	t580 = -t570 * t602 - t625 * t573;
	t609 = t580 * t630;
	t608 = t565 * t630;
	t607 = t630 * t570;
	t606 = t630 * t573;
	t605 = qJD(2) * t622;
	t604 = t562 * t626;
	t603 = t563 * t625;
	t600 = t565 * t606;
	t599 = t570 * t605;
	t598 = t573 * t605;
	t596 = t569 * t604;
	t595 = t562 * t603;
	t593 = t626 * t625;
	t592 = t563 * t600;
	t521 = -t609 + (t546 * t565 - t612) * t569;
	t533 = -t546 * t562 - t565 * t618;
	t507 = t521 * t560 + t533 * t559;
	t591 = -t521 * t559 + t533 * t560;
	t578 = t564 * t570 + t573 * t593;
	t579 = -t564 * t573 + t570 * t593;
	t523 = -t579 * t630 + (-t565 * t578 + t595) * t569;
	t534 = t562 * t578 + t565 * t603;
	t509 = t523 * t560 + t534 * t559;
	t590 = -t523 * t559 + t534 * t560;
	t582 = t565 * t615 + t607;
	t532 = t582 * t563 + t596;
	t545 = -t563 * t619 + t626 * t565;
	t517 = t532 * t560 + t545 * t559;
	t589 = -t532 * t559 + t545 * t560;
	t588 = r_i_i_C(1) * t571 - r_i_i_C(2) * t567 + pkin(5);
	t587 = t630 * t604;
	t527 = t546 * t630 + t580 * t617;
	t514 = t527 * t560 - t580 * t624;
	t529 = -t578 * t630 + t579 * t617;
	t515 = t529 * t560 - t579 * t624;
	t583 = -t565 * t616 + t606;
	t540 = t583 * t563;
	t530 = t540 * t560 + t559 * t611;
	t585 = -t546 * t569 + t580 * t608;
	t584 = t569 * t578 + t579 * t608;
	t581 = t565 * t607 + t615;
	t558 = pkin(4) * t572 + pkin(3);
	t577 = -t631 * t559 - t588 * t560 - t558;
	t576 = t546 * t608 + t569 * t580 - t630 * t612;
	t575 = t569 * t579 - t578 * t608 + t630 * t595;
	t574 = t560 * t586 + (t588 * t559 - t631 * t560 + t628) * qJD(4);
	t566 = -qJ(5) - pkin(10);
	t544 = t579 * qJD(2);
	t543 = t578 * qJD(2);
	t542 = t580 * qJD(2);
	t541 = t546 * qJD(2);
	t539 = t581 * t563;
	t531 = -t587 - t592 + t610;
	t525 = (-t582 * qJD(2) - t581 * qJD(3)) * t563;
	t524 = -qJD(2) * t592 - t563 * qJD(3) * t606 + (qJD(3) * t565 + qJD(2)) * t610;
	t519 = qJD(3) * t587 + ((t600 - t616) * qJD(3) + t583 * qJD(2)) * t563;
	t518 = qJD(3) * t596 + (t581 * qJD(2) + t582 * qJD(3)) * t563;
	t513 = t584 * qJD(3) + t543 * t617 + t544 * t630;
	t512 = t529 * qJD(3) - t543 * t608 + t544 * t569;
	t511 = t585 * qJD(3) - t541 * t617 + t542 * t630;
	t510 = t527 * qJD(3) + t541 * t608 + t542 * t569;
	t505 = t575 * qJD(3) - t543 * t630 + t544 * t617;
	t504 = t523 * qJD(3) - t543 * t569 - t544 * t608;
	t503 = t576 * qJD(3) + t541 * t630 + t542 * t617;
	t502 = t541 * t569 - t542 * t608 + (t546 * t617 - t569 * t612 - t609) * qJD(3);
	t501 = t559 * t598 + t525 * t560 + (-t540 * t559 + t560 * t611) * qJD(4);
	t499 = t589 * qJD(4) + t519 * t560 + t559 * t599;
	t497 = -t543 * t624 + t513 * t560 + (-t529 * t559 - t579 * t623) * qJD(4);
	t495 = t541 * t624 + t511 * t560 + (-t527 * t559 - t580 * t623) * qJD(4);
	t493 = t590 * qJD(4) + t505 * t560 - t544 * t624;
	t491 = t591 * qJD(4) + t503 * t560 - t542 * t624;
	t1 = [0, (t497 * t571 + t512 * t567) * r_i_i_C(1) + (-t497 * t567 + t512 * t571) * r_i_i_C(2) + t497 * pkin(5) + t513 * t558 - t512 * t566 - t584 * qJD(5) + t544 * pkin(2) - t543 * t629 + t631 * (t515 * qJD(4) + t513 * t559 + t543 * t623) + ((-t515 * t567 - t571 * t584) * r_i_i_C(1) + (-t515 * t571 + t567 * t584) * r_i_i_C(2)) * qJD(6) + (-t543 * t621 + (-t529 * t568 - t579 * t620) * qJD(4)) * pkin(4), (t505 * t567 + t523 * t613) * r_i_i_C(1) + (t505 * t571 - t523 * t614) * r_i_i_C(2) - t505 * t566 + t523 * qJD(5) + t577 * t504 - t574 * t575, t631 * t493 - t590 * t586 + t588 * (-t509 * qJD(4) - t505 * t559 - t544 * t623) + (-t544 * t620 - t505 * t568 + (-t523 * t572 - t534 * t568) * qJD(4)) * pkin(4), t504, (-t493 * t567 + t504 * t571) * r_i_i_C(1) + (-t493 * t571 - t504 * t567) * r_i_i_C(2) + ((-t509 * t571 + t567 * t575) * r_i_i_C(1) + (t509 * t567 + t571 * t575) * r_i_i_C(2)) * qJD(6); 0, (t495 * t571 + t510 * t567) * r_i_i_C(1) + (-t495 * t567 + t510 * t571) * r_i_i_C(2) + t495 * pkin(5) + t511 * t558 - t510 * t566 - t585 * qJD(5) + t542 * pkin(2) + t541 * t629 + t631 * (t514 * qJD(4) + t511 * t559 - t541 * t623) + ((-t514 * t567 - t571 * t585) * r_i_i_C(1) + (-t514 * t571 + t567 * t585) * r_i_i_C(2)) * qJD(6) + (t541 * t621 + (-t527 * t568 - t580 * t620) * qJD(4)) * pkin(4), (t503 * t567 + t521 * t613) * r_i_i_C(1) + (t503 * t571 - t521 * t614) * r_i_i_C(2) - t503 * t566 + t521 * qJD(5) + t577 * t502 - t574 * t576, t631 * t491 - t591 * t586 + t588 * (-t507 * qJD(4) - t503 * t559 - t542 * t623) + (-t542 * t620 - t503 * t568 + (-t521 * t572 - t533 * t568) * qJD(4)) * pkin(4), t502, (-t491 * t567 + t502 * t571) * r_i_i_C(1) + (-t491 * t571 - t502 * t567) * r_i_i_C(2) + ((-t507 * t571 + t567 * t576) * r_i_i_C(1) + (t507 * t567 + t571 * t576) * r_i_i_C(2)) * qJD(6); 0, (t501 * t571 - t524 * t567) * r_i_i_C(1) + (-t501 * t567 - t524 * t571) * r_i_i_C(2) + t501 * pkin(5) + t525 * t558 - t540 * t568 * t627 + t524 * t566 + t539 * qJD(5) + t631 * (t530 * qJD(4) + t525 * t559 - t560 * t598) + ((-t530 * t567 + t539 * t571) * r_i_i_C(1) + (-t530 * t571 - t539 * t567) * r_i_i_C(2)) * qJD(6) + (t570 * t620 * t627 + (-pkin(2) * t570 + (pkin(9) + t628) * t619) * qJD(2)) * t563, (t519 * t567 + t532 * t613) * r_i_i_C(1) + (t519 * t571 - t532 * t614) * r_i_i_C(2) - t519 * t566 + t532 * qJD(5) + t577 * t518 + t574 * t531, t631 * t499 - t589 * t586 + t588 * (-t517 * qJD(4) - t519 * t559 + t560 * t599) + (t572 * t599 - t519 * t568 + (-t532 * t572 - t545 * t568) * qJD(4)) * pkin(4), t518, (-t499 * t567 + t518 * t571) * r_i_i_C(1) + (-t499 * t571 - t518 * t567) * r_i_i_C(2) + ((-t517 * t571 - t531 * t567) * r_i_i_C(1) + (t517 * t567 - t531 * t571) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end