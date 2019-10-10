% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:36
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
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
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
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:37
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
	% StartTime: 2019-10-09 21:44:38
	% EndTime: 2019-10-09 21:44:38
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (429->82), mult. (1352->159), div. (0->0), fcn. (1456->12), ass. (0->57)
	t379 = sin(qJ(4));
	t382 = cos(qJ(4));
	t378 = sin(qJ(5));
	t381 = cos(qJ(5));
	t389 = (t378 * r_i_i_C(1) + t381 * r_i_i_C(2)) * qJD(5);
	t393 = t381 * r_i_i_C(1) - t378 * r_i_i_C(2) + pkin(4);
	t412 = pkin(9) + r_i_i_C(3);
	t416 = (t393 * t379 - t412 * t382) * qJD(4) + t382 * t389;
	t380 = sin(qJ(2));
	t373 = sin(pkin(11));
	t383 = cos(qJ(2));
	t405 = t383 * t373;
	t410 = cos(pkin(11));
	t366 = -t380 * t410 - t405;
	t377 = cos(pkin(6));
	t362 = t366 * t377;
	t399 = qJD(2) * t410;
	t404 = qJD(2) * t380;
	t414 = t373 * t404 - t383 * t399;
	t388 = -t380 * t373 + t383 * t410;
	t385 = t412 * t379 + t393 * t382 + pkin(3);
	t411 = pkin(2) * qJD(2);
	t375 = sin(pkin(6));
	t409 = t375 * t379;
	t408 = t375 * t382;
	t407 = t377 * t380;
	t403 = qJD(5) * t378;
	t402 = qJD(5) * t381;
	t359 = t414 * t377;
	t364 = -qJD(2) * t405 - t380 * t399;
	t374 = sin(pkin(10));
	t376 = cos(pkin(10));
	t342 = t376 * t359 - t374 * t364;
	t344 = t374 * t359 + t376 * t364;
	t361 = t366 * t375;
	t353 = -t361 * t382 + t377 * t379;
	t396 = t361 * t379 + t377 * t382;
	t395 = -t376 * t362 + t374 * t388;
	t394 = t374 * t362 + t376 * t388;
	t392 = -t376 * t408 - t379 * t395;
	t391 = t376 * t409 - t382 * t395;
	t390 = t374 * t408 - t379 * t394;
	t339 = t374 * t409 + t382 * t394;
	t387 = t388 * t377;
	t386 = qJD(2) * t362;
	t363 = t388 * qJD(2);
	t360 = t388 * t375;
	t358 = qJD(2) * t361;
	t357 = t414 * t375;
	t350 = t376 * t366 - t374 * t387;
	t347 = t374 * t366 + t376 * t387;
	t343 = -t376 * t363 - t374 * t386;
	t340 = -t374 * t363 + t376 * t386;
	t335 = qJD(4) * t396 - t357 * t382;
	t333 = qJD(4) * t390 + t344 * t382;
	t331 = qJD(4) * t392 - t342 * t382;
	t1 = [0, (t344 * t378 + t394 * t402) * r_i_i_C(1) + (t344 * t381 - t394 * t403) * r_i_i_C(2) + t344 * pkin(8) + (t374 * t407 - t376 * t383) * t411 + t385 * t343 - t416 * t350, 0, t412 * t333 - t390 * t389 + t393 * (-qJD(4) * t339 - t344 * t379), (-t333 * t378 - t343 * t381) * r_i_i_C(1) + (-t333 * t381 + t343 * t378) * r_i_i_C(2) + ((-t339 * t381 + t350 * t378) * r_i_i_C(1) + (t339 * t378 + t350 * t381) * r_i_i_C(2)) * qJD(5), 0; 0, (-t342 * t378 + t395 * t402) * r_i_i_C(1) + (-t342 * t381 - t395 * t403) * r_i_i_C(2) - t342 * pkin(8) + (-t374 * t383 - t376 * t407) * t411 + t385 * t340 - t416 * t347, 0, t412 * t331 - t392 * t389 + t393 * (qJD(4) * t391 + t342 * t379), (-t331 * t378 - t340 * t381) * r_i_i_C(1) + (-t331 * t381 + t340 * t378) * r_i_i_C(2) + ((t347 * t378 + t381 * t391) * r_i_i_C(1) + (t347 * t381 - t378 * t391) * r_i_i_C(2)) * qJD(5), 0; 0, (-t357 * t378 - t361 * t402) * r_i_i_C(1) + (-t357 * t381 + t361 * t403) * r_i_i_C(2) - t357 * pkin(8) - t375 * pkin(2) * t404 + t385 * t358 - t416 * t360, 0, t412 * t335 - t396 * t389 + t393 * (-qJD(4) * t353 + t357 * t379), (-t335 * t378 - t358 * t381) * r_i_i_C(1) + (-t335 * t381 + t358 * t378) * r_i_i_C(2) + ((-t353 * t381 + t360 * t378) * r_i_i_C(1) + (t353 * t378 + t360 * t381) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:40
	% EndTime: 2019-10-09 21:44:40
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (822->101), mult. (2522->179), div. (0->0), fcn. (2797->12), ass. (0->71)
	t451 = sin(qJ(4));
	t454 = cos(qJ(4));
	t497 = pkin(9) + r_i_i_C(2);
	t463 = qJD(4) * (-pkin(4) * t451 + t497 * t454);
	t449 = cos(pkin(6));
	t445 = sin(pkin(11));
	t452 = sin(qJ(2));
	t491 = cos(pkin(11));
	t495 = cos(qJ(2));
	t462 = t495 * t445 + t452 * t491;
	t434 = t462 * t449;
	t478 = t495 * t491;
	t485 = qJD(2) * t452;
	t499 = -qJD(2) * t478 + t445 * t485;
	t461 = -t452 * t445 + t478;
	t464 = t454 * pkin(4) + t497 * t451 + pkin(3);
	t496 = r_i_i_C(1) + pkin(5);
	t493 = r_i_i_C(3) + qJ(6);
	t492 = pkin(2) * qJD(2);
	t447 = sin(pkin(6));
	t490 = t447 * t451;
	t489 = t447 * t454;
	t488 = t449 * t452;
	t450 = sin(qJ(5));
	t487 = t450 * t454;
	t484 = qJD(4) * t451;
	t483 = qJD(5) * t454;
	t430 = t499 * t449;
	t436 = t462 * qJD(2);
	t446 = sin(pkin(10));
	t448 = cos(pkin(10));
	t413 = t448 * t430 + t446 * t436;
	t456 = t461 * t449;
	t418 = -t446 * t462 + t448 * t456;
	t477 = t418 * t483 + t413;
	t415 = t446 * t430 - t448 * t436;
	t421 = -t446 * t456 - t448 * t462;
	t476 = t421 * t483 - t415;
	t428 = t499 * t447;
	t432 = t461 * t447;
	t475 = -t432 * t483 - t428;
	t453 = cos(qJ(5));
	t470 = t448 * t434 + t446 * t461;
	t466 = t448 * t490 - t454 * t470;
	t474 = -t418 * t450 - t453 * t466;
	t469 = -t446 * t434 + t448 * t461;
	t408 = t446 * t490 + t454 * t469;
	t473 = t408 * t453 - t421 * t450;
	t433 = t462 * t447;
	t424 = t433 * t454 + t449 * t451;
	t472 = t424 * t453 - t432 * t450;
	t471 = -t433 * t451 + t449 * t454;
	t467 = -t448 * t489 - t451 * t470;
	t465 = t446 * t489 - t451 * t469;
	t460 = t493 * t450 + t496 * t453 + pkin(4);
	t431 = qJD(2) * t434;
	t435 = t461 * qJD(2);
	t411 = -t448 * t431 - t446 * t435;
	t459 = qJD(5) * t470 + t411 * t454 - t418 * t484;
	t414 = t446 * t431 - t448 * t435;
	t458 = qJD(5) * t469 + t414 * t454 - t421 * t484;
	t429 = qJD(2) * t433;
	t457 = qJD(5) * t433 - t429 * t454 - t432 * t484;
	t455 = qJD(6) * t450 + (-t496 * t450 + t493 * t453) * qJD(5);
	t404 = t471 * qJD(4) - t428 * t454;
	t402 = t465 * qJD(4) + t415 * t454;
	t400 = t467 * qJD(4) - t413 * t454;
	t395 = t472 * qJD(5) + t404 * t450 - t429 * t453;
	t389 = t473 * qJD(5) + t402 * t450 + t414 * t453;
	t387 = t474 * qJD(5) + t400 * t450 + t411 * t453;
	t1 = [0, -(-t421 * t487 + t453 * t469) * qJD(6) + t415 * pkin(8) + t496 * (-t476 * t450 + t458 * t453) + t493 * (t458 * t450 + t476 * t453) + (t446 * t488 - t495 * t448) * t492 + t421 * t463 + t464 * t414, 0, t497 * t402 + t455 * t465 + t460 * (-t408 * qJD(4) - t415 * t451), t473 * qJD(6) + t493 * (t402 * t453 - t414 * t450 + (-t408 * t450 - t421 * t453) * qJD(5)) - t496 * t389, t389; 0, -(-t418 * t487 + t453 * t470) * qJD(6) - t413 * pkin(8) + t496 * (-t477 * t450 + t459 * t453) + t493 * (t459 * t450 + t477 * t453) + (-t495 * t446 - t448 * t488) * t492 + t418 * t463 + t464 * t411, 0, t497 * t400 + t455 * t467 + t460 * (t466 * qJD(4) + t413 * t451), t474 * qJD(6) + t493 * (t400 * t453 - t411 * t450 + (-t418 * t453 + t450 * t466) * qJD(5)) - t496 * t387, t387; 0, -(-t432 * t487 + t433 * t453) * qJD(6) - t428 * pkin(8) - t447 * pkin(2) * t485 + t496 * (t475 * t450 + t457 * t453) + t493 * (t457 * t450 - t475 * t453) + t432 * t463 - t464 * t429, 0, t497 * t404 + t455 * t471 + t460 * (-t424 * qJD(4) + t428 * t451), t472 * qJD(6) + t493 * (t404 * t453 + t429 * t450 + (-t424 * t450 - t432 * t453) * qJD(5)) - t496 * t395, t395;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end