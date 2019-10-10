% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (180->39), mult. (333->75), div. (0->0), fcn. (304->10), ass. (0->38)
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t229 = sin(qJ(2));
	t227 = cos(pkin(6));
	t231 = cos(qJ(2));
	t247 = t227 * t231;
	t256 = -t224 * t229 + t226 * t247;
	t255 = r_i_i_C(3) + pkin(9) + pkin(8);
	t223 = qJ(3) + qJ(4);
	t220 = sin(t223);
	t222 = qJD(3) + qJD(4);
	t254 = t220 * t222;
	t221 = cos(t223);
	t253 = t221 * t222;
	t225 = sin(pkin(6));
	t252 = t222 * t225;
	t228 = sin(qJ(3));
	t250 = t225 * t228;
	t249 = t225 * t229;
	t248 = t227 * t229;
	t215 = t224 * t231 + t226 * t248;
	t210 = t256 * qJD(2);
	t240 = t226 * t252 - t210;
	t246 = (-t215 * t253 + t240 * t220) * r_i_i_C(1) + (t215 * t254 + t240 * t221) * r_i_i_C(2);
	t236 = t224 * t248 - t226 * t231;
	t237 = t224 * t247 + t226 * t229;
	t212 = t237 * qJD(2);
	t239 = -t224 * t252 + t212;
	t245 = (t239 * t220 + t236 * t253) * r_i_i_C(1) + (t239 * t221 - t236 * t254) * r_i_i_C(2);
	t241 = qJD(2) * t225 * t231;
	t235 = -t222 * t227 - t241;
	t243 = t222 * t249;
	t244 = (t235 * t220 - t221 * t243) * r_i_i_C(1) + (t220 * t243 + t235 * t221) * r_i_i_C(2);
	t230 = cos(qJ(3));
	t238 = t230 * pkin(3) + r_i_i_C(1) * t221 - r_i_i_C(2) * t220 + pkin(2);
	t234 = qJD(2) * t238;
	t233 = -pkin(3) * qJD(3) * t228 + (-r_i_i_C(1) * t220 - r_i_i_C(2) * t221) * t222;
	t1 = [0, -t255 * t212 - t233 * t237 + t236 * t234, (t212 * t228 + (-t224 * t250 + t230 * t236) * qJD(3)) * pkin(3) + t245, t245, 0, 0; 0, t255 * t210 - t215 * t234 + t233 * t256, (-t210 * t228 + (-t215 * t230 + t226 * t250) * qJD(3)) * pkin(3) + t246, t246, 0, 0; 0, (t233 * t231 + (-t238 * t229 + t255 * t231) * qJD(2)) * t225, (-t228 * t241 + (-t227 * t228 - t230 * t249) * qJD(3)) * pkin(3) + t244, t244, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:58
	% EndTime: 2019-10-09 23:03:59
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (593->88), mult. (1038->158), div. (0->0), fcn. (1028->12), ass. (0->61)
	t406 = pkin(10) + r_i_i_C(3);
	t362 = sin(qJ(5));
	t365 = cos(qJ(5));
	t379 = t365 * r_i_i_C(1) - t362 * r_i_i_C(2);
	t412 = pkin(4) + t379;
	t392 = qJD(5) * t365;
	t393 = qJD(5) * t362;
	t411 = -r_i_i_C(1) * t393 - t392 * r_i_i_C(2);
	t359 = qJ(3) + qJ(4);
	t356 = sin(t359);
	t357 = cos(t359);
	t358 = qJD(3) + qJD(4);
	t363 = sin(qJ(3));
	t410 = -qJD(3) * t363 * pkin(3) - (t356 * t412 - t406 * t357) * t358;
	t364 = sin(qJ(2));
	t367 = cos(qJ(2));
	t360 = sin(pkin(11));
	t402 = cos(pkin(6));
	t385 = t360 * t402;
	t401 = cos(pkin(11));
	t349 = -t364 * t385 + t401 * t367;
	t408 = t362 * r_i_i_C(1) + t365 * r_i_i_C(2);
	t400 = t356 * t358;
	t399 = t357 * t358;
	t361 = sin(pkin(6));
	t398 = t360 * t361;
	t397 = t361 * t364;
	t396 = t361 * t367;
	t395 = qJD(2) * t364;
	t394 = qJD(5) * t357;
	t389 = t356 * t397;
	t388 = t357 * t397;
	t387 = t361 * t395;
	t386 = qJD(2) * t396;
	t384 = t361 * t401;
	t380 = t357 * t384;
	t348 = t401 * t364 + t367 * t385;
	t344 = t348 * qJD(2);
	t378 = t358 * t398 - t344;
	t377 = t402 * t401;
	t375 = t367 * t377;
	t374 = t402 * t358 + t386;
	t347 = t360 * t367 + t364 * t377;
	t366 = cos(qJ(3));
	t373 = -t366 * pkin(3) - t406 * t356 - t357 * t412 - pkin(2);
	t342 = -qJD(2) * t375 + t360 * t395;
	t326 = -t342 * t357 - t347 * t400 - t358 * t380;
	t372 = t411 * (-t347 * t356 - t380) + t406 * t326 + t412 * (-t347 * t399 + (t358 * t384 + t342) * t356);
	t328 = -t349 * t400 + t378 * t357;
	t371 = t411 * (-t349 * t356 + t357 * t398) + t406 * t328 + t412 * (-t349 * t399 - t378 * t356);
	t333 = t374 * t357 - t358 * t389;
	t370 = t411 * (t402 * t357 - t389) + t406 * t333 + t412 * (-t374 * t356 - t358 * t388);
	t369 = t408 * t394 - t410;
	t368 = -pkin(9) - pkin(8);
	t346 = t360 * t364 - t375;
	t345 = t349 * qJD(2);
	t343 = t347 * qJD(2);
	t341 = t402 * t356 + t388;
	t337 = t349 * t357 + t356 * t398;
	t335 = t347 * t357 - t356 * t384;
	t1 = [0, (-t344 * t362 + t349 * t392) * r_i_i_C(1) + (-t344 * t365 - t349 * t393) * r_i_i_C(2) + t344 * t368 + t373 * t345 + t369 * t348, (t344 * t363 + (-t349 * t366 - t363 * t398) * qJD(3)) * pkin(3) + t371, t371, (-t328 * t362 + t345 * t365) * r_i_i_C(1) + (-t328 * t365 - t345 * t362) * r_i_i_C(2) + ((-t337 * t365 - t348 * t362) * r_i_i_C(1) + (t337 * t362 - t348 * t365) * r_i_i_C(2)) * qJD(5), 0; 0, (-t342 * t362 + t347 * t392) * r_i_i_C(1) + (-t342 * t365 - t347 * t393) * r_i_i_C(2) + t342 * t368 + t373 * t343 + t369 * t346, (t342 * t363 + (-t347 * t366 + t363 * t384) * qJD(3)) * pkin(3) + t372, t372, (-t326 * t362 + t343 * t365) * r_i_i_C(1) + (-t326 * t365 - t343 * t362) * r_i_i_C(2) + ((-t335 * t365 - t346 * t362) * r_i_i_C(1) + (t335 * t362 - t346 * t365) * r_i_i_C(2)) * qJD(5), 0; 0, ((t373 * qJD(2) + t379 * qJD(5)) * t364 + (-qJD(2) * t368 + t408 * (qJD(2) - t394) + t410) * t367) * t361, (-t363 * t386 + (-t402 * t363 - t366 * t397) * qJD(3)) * pkin(3) + t370, t370, (-t333 * t362 + t365 * t387) * r_i_i_C(1) + (-t333 * t365 - t362 * t387) * r_i_i_C(2) + ((-t341 * t365 + t362 * t396) * r_i_i_C(1) + (t341 * t362 + t365 * t396) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:04:00
	% EndTime: 2019-10-09 23:04:01
	% DurationCPUTime: 0.87s
	% Computational Cost: add. (1048->111), mult. (1833->184), div. (0->0), fcn. (1881->12), ass. (0->79)
	t462 = cos(qJ(5));
	t518 = -r_i_i_C(1) - pkin(5);
	t528 = -t518 * t462 + pkin(4);
	t517 = r_i_i_C(3) + qJ(6);
	t519 = pkin(10) + r_i_i_C(2);
	t459 = sin(qJ(5));
	t500 = qJD(6) * t459;
	t527 = qJD(5) * t459 * t518 + t500;
	t456 = qJ(3) + qJ(4);
	t454 = cos(t456);
	t455 = qJD(3) + qJD(4);
	t460 = sin(qJ(3));
	t453 = sin(t456);
	t508 = t453 * t455;
	t522 = -qJD(3) * t460 * pkin(3) - pkin(4) * t508 + (t519 * t455 + t500) * t454;
	t461 = sin(qJ(2));
	t464 = cos(qJ(2));
	t457 = sin(pkin(11));
	t516 = cos(pkin(6));
	t489 = t457 * t516;
	t515 = cos(pkin(11));
	t444 = -t461 * t489 + t515 * t464;
	t520 = (qJD(2) * t454 - qJD(5)) * t461 + t464 * t508;
	t507 = t454 * t455;
	t458 = sin(pkin(6));
	t506 = t457 * t458;
	t505 = t458 * t461;
	t504 = t458 * t464;
	t503 = qJD(2) * t461;
	t502 = qJD(5) * t454;
	t501 = qJD(5) * t462;
	t499 = t462 * qJD(6);
	t498 = t453 * t505;
	t497 = t454 * t505;
	t496 = t459 * t504;
	t494 = t458 * t503;
	t493 = qJD(2) * t504;
	t488 = t458 * t515;
	t483 = t454 * t488;
	t443 = t515 * t461 + t464 * t489;
	t439 = t443 * qJD(2);
	t482 = t455 * t506 - t439;
	t479 = t516 * t515;
	t475 = t464 * t479;
	t437 = -qJD(2) * t475 + t457 * t503;
	t441 = t457 * t461 - t475;
	t481 = t441 * t502 - t437;
	t480 = t443 * t502 - t439;
	t442 = t457 * t464 + t461 * t479;
	t425 = t442 * t454 - t453 * t488;
	t478 = t425 * t462 + t441 * t459;
	t427 = t444 * t454 + t453 * t506;
	t477 = t427 * t462 + t443 * t459;
	t476 = (qJD(2) - t502) * t464;
	t463 = cos(qJ(3));
	t473 = -t463 * pkin(3) - t454 * pkin(4) - t519 * t453 - pkin(2);
	t472 = t516 * t455 + t493;
	t438 = t442 * qJD(2);
	t471 = qJD(5) * t442 - t438 * t454 + t441 * t508;
	t440 = t444 * qJD(2);
	t470 = qJD(5) * t444 - t440 * t454 + t443 * t508;
	t406 = -t442 * t507 + (t455 * t488 + t437) * t453;
	t407 = -t437 * t454 - t442 * t508 - t455 * t483;
	t424 = -t442 * t453 - t483;
	t468 = t519 * t407 + t527 * t424 + t517 * (t406 * t459 + t424 * t501) + t528 * t406;
	t408 = -t444 * t507 - t482 * t453;
	t409 = -t444 * t508 + t482 * t454;
	t426 = -t444 * t453 + t454 * t506;
	t467 = t519 * t409 + t527 * t426 + t517 * (t408 * t459 + t426 * t501) + t528 * t408;
	t418 = -t472 * t453 - t455 * t497;
	t419 = t472 * t454 - t455 * t498;
	t435 = t516 * t454 - t498;
	t466 = t519 * t419 + t527 * t435 + t517 * (t418 * t459 + t435 * t501) + t528 * t418;
	t465 = -pkin(9) - pkin(8);
	t436 = t516 * t453 + t497;
	t390 = -qJD(5) * t496 + t419 * t459 + t436 * t501 - t462 * t494;
	t384 = t477 * qJD(5) + t409 * t459 - t440 * t462;
	t382 = t478 * qJD(5) + t407 * t459 - t438 * t462;
	t1 = [0, -t444 * t499 + t439 * t465 - t518 * (t480 * t459 + t470 * t462) + t517 * (t470 * t459 - t480 * t462) - t522 * t443 + t473 * t440, (t439 * t460 + (-t444 * t463 - t460 * t506) * qJD(3)) * pkin(3) + t467, t467, t477 * qJD(6) + t517 * (t409 * t462 + t440 * t459 + (-t427 * t459 + t443 * t462) * qJD(5)) + t518 * t384, t384; 0, -t442 * t499 + t437 * t465 - t518 * (t481 * t459 + t471 * t462) + t517 * (t471 * t459 - t481 * t462) - t522 * t441 + t473 * t438, (t437 * t460 + (-t442 * t463 + t460 * t488) * qJD(3)) * pkin(3) + t468, t468, t478 * qJD(6) + t517 * (t407 * t462 + t438 * t459 + (-t425 * t459 + t441 * t462) * qJD(5)) + t518 * t382, t382; 0, (-t518 * (t459 * t476 - t520 * t462) - t517 * (t520 * t459 + t462 * t476) - t461 * t499 + t522 * t464 + (t473 * t461 - t464 * t465) * qJD(2)) * t458, (-t460 * t493 + (-t516 * t460 - t463 * t505) * qJD(3)) * pkin(3) + t466, t466, -(-t436 * t462 + t496) * qJD(6) + t517 * (t459 * t494 + t419 * t462 + (-t436 * t459 - t462 * t504) * qJD(5)) + t518 * t390, t390;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end