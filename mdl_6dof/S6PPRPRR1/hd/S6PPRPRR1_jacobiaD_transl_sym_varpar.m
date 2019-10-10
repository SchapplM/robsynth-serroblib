% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (23->19), mult. (86->45), div. (0->0), fcn. (88->10), ass. (0->22)
	t111 = sin(pkin(11));
	t117 = cos(pkin(6));
	t126 = t111 * t117;
	t112 = sin(pkin(7));
	t113 = sin(pkin(6));
	t125 = t112 * t113;
	t115 = cos(pkin(11));
	t124 = t115 * t117;
	t116 = cos(pkin(7));
	t118 = sin(qJ(3));
	t123 = t116 * t118;
	t119 = cos(qJ(3));
	t122 = t116 * t119;
	t121 = t118 * t125;
	t120 = t119 * t125;
	t114 = cos(pkin(12));
	t110 = sin(pkin(12));
	t109 = -t110 * t126 + t115 * t114;
	t108 = -t115 * t110 - t114 * t126;
	t107 = t110 * t124 + t111 * t114;
	t106 = -t111 * t110 + t114 * t124;
	t1 = [0, 0, ((-t108 * t123 - t109 * t119 - t111 * t121) * r_i_i_C(1) + (-t108 * t122 + t109 * t118 - t111 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-t106 * t123 - t107 * t119 + t115 * t121) * r_i_i_C(1) + (-t106 * t122 + t107 * t118 + t115 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-r_i_i_C(1) * t118 - r_i_i_C(2) * t119) * t117 * t112 + ((-t110 * t119 - t114 * t123) * r_i_i_C(1) + (t110 * t118 - t114 * t122) * r_i_i_C(2)) * t113) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (54->30), mult. (201->65), div. (0->0), fcn. (198->12), ass. (0->28)
	t150 = sin(pkin(13));
	t155 = cos(pkin(13));
	t160 = sin(qJ(3));
	t161 = cos(qJ(3));
	t169 = (t150 * t161 + t155 * t160) * qJD(3);
	t168 = pkin(3) * qJD(3);
	t152 = sin(pkin(11));
	t154 = sin(pkin(6));
	t167 = t152 * t154;
	t159 = cos(pkin(6));
	t166 = t152 * t159;
	t157 = cos(pkin(11));
	t165 = t154 * t157;
	t164 = t157 * t159;
	t148 = (t150 * t160 - t155 * t161) * qJD(3);
	t158 = cos(pkin(7));
	t156 = cos(pkin(12));
	t153 = sin(pkin(7));
	t151 = sin(pkin(12));
	t147 = -t151 * t166 + t157 * t156;
	t146 = -t157 * t151 - t156 * t166;
	t145 = t151 * t164 + t152 * t156;
	t144 = -t152 * t151 + t156 * t164;
	t143 = t158 * t169;
	t142 = t158 * t148;
	t141 = t153 * t169;
	t140 = t153 * t148;
	t1 = [0, 0, (-t141 * t167 - t146 * t143 + t147 * t148) * r_i_i_C(1) + (t140 * t167 + t146 * t142 + t147 * t169) * r_i_i_C(2) + (-t147 * t161 + (-t146 * t158 - t153 * t167) * t160) * t168, 0, 0, 0; 0, 0, (t141 * t165 - t144 * t143 + t145 * t148) * r_i_i_C(1) + (-t140 * t165 + t144 * t142 + t145 * t169) * r_i_i_C(2) + (-t145 * t161 + (-t144 * t158 + t153 * t165) * t160) * t168, 0, 0, 0; 0, 0, (-t153 * t160 * t168 - t141 * r_i_i_C(1) + t140 * r_i_i_C(2)) * t159 + ((-t143 * t156 + t148 * t151) * r_i_i_C(1) + (t142 * t156 + t151 * t169) * r_i_i_C(2) + (-t156 * t158 * t160 - t151 * t161) * t168) * t154, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:45
	% EndTime: 2019-10-09 21:08:45
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (265->66), mult. (885->136), div. (0->0), fcn. (967->14), ass. (0->55)
	t337 = sin(pkin(13));
	t342 = cos(pkin(13));
	t350 = cos(qJ(3));
	t362 = qJD(3) * t350;
	t348 = sin(qJ(3));
	t363 = qJD(3) * t348;
	t374 = t337 * t363 - t342 * t362;
	t340 = sin(pkin(7));
	t316 = t374 * t340;
	t345 = cos(pkin(7));
	t318 = t374 * t345;
	t330 = -t337 * t362 - t342 * t363;
	t338 = sin(pkin(12));
	t341 = sin(pkin(6));
	t343 = cos(pkin(12));
	t346 = cos(pkin(6));
	t373 = t346 * t316 + (t318 * t343 - t330 * t338) * t341;
	t372 = -pkin(9) - r_i_i_C(3);
	t371 = pkin(3) * qJD(3);
	t339 = sin(pkin(11));
	t370 = t339 * t341;
	t369 = t339 * t346;
	t368 = t340 * t341;
	t344 = cos(pkin(11));
	t367 = t341 * t344;
	t366 = t341 * t345;
	t365 = t344 * t346;
	t347 = sin(qJ(5));
	t349 = cos(qJ(5));
	t359 = -t347 * r_i_i_C(1) - t349 * r_i_i_C(2);
	t357 = t350 * t337 + t348 * t342;
	t356 = t348 * t337 - t350 * t342;
	t355 = t349 * r_i_i_C(1) - t347 * r_i_i_C(2) + pkin(4);
	t354 = qJD(5) * t359;
	t353 = qJD(3) * t357;
	t325 = -t339 * t338 + t343 * t365;
	t326 = t338 * t365 + t339 * t343;
	t352 = t316 * t367 - t325 * t318 + t326 * t330;
	t327 = -t344 * t338 - t343 * t369;
	t328 = -t338 * t369 + t344 * t343;
	t351 = t316 * t370 + t327 * t318 - t328 * t330;
	t329 = t356 * qJD(3);
	t324 = -t343 * t368 + t346 * t345;
	t323 = t357 * t345;
	t322 = t356 * t345;
	t321 = t357 * t340;
	t320 = t356 * t340;
	t319 = t345 * t353;
	t317 = t340 * t353;
	t315 = -t327 * t340 + t339 * t366;
	t314 = -t325 * t340 - t344 * t366;
	t313 = t346 * t321 + (t323 * t343 - t338 * t356) * t341;
	t308 = t321 * t370 + t327 * t323 - t328 * t356;
	t306 = -t321 * t367 + t325 * t323 - t326 * t356;
	t1 = [0, 0, t372 * t351 + (-t320 * t370 - t327 * t322 - t328 * t357) * t354 + t355 * (-t317 * t370 - t327 * t319 + t328 * t329) + (-t328 * t350 + (-t327 * t345 - t339 * t368) * t348) * t371, 0, -t359 * t351 + ((-t308 * t349 - t315 * t347) * r_i_i_C(1) + (t308 * t347 - t315 * t349) * r_i_i_C(2)) * qJD(5), 0; 0, 0, -t372 * t352 + (t320 * t367 - t325 * t322 - t326 * t357) * t354 + t355 * (t317 * t367 - t325 * t319 + t326 * t329) + (-t326 * t350 + (-t325 * t345 + t340 * t367) * t348) * t371, 0, t359 * t352 + ((-t306 * t349 - t314 * t347) * r_i_i_C(1) + (t306 * t347 - t314 * t349) * r_i_i_C(2)) * qJD(5), 0; 0, 0, t372 * t373 + (-t346 * t320 + (-t322 * t343 - t338 * t357) * t341) * t354 + t355 * (-t346 * t317 + (-t319 * t343 + t329 * t338) * t341) + (-t340 * t346 * t348 + (-t343 * t345 * t348 - t338 * t350) * t341) * t371, 0, -t359 * t373 + ((-t313 * t349 - t324 * t347) * r_i_i_C(1) + (t313 * t347 - t324 * t349) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:47
	% EndTime: 2019-10-09 21:08:48
	% DurationCPUTime: 0.85s
	% Computational Cost: add. (972->109), mult. (3071->210), div. (0->0), fcn. (3566->16), ass. (0->75)
	t483 = sin(pkin(13));
	t488 = cos(pkin(13));
	t498 = cos(qJ(3));
	t516 = qJD(3) * t498;
	t495 = sin(qJ(3));
	t517 = qJD(3) * t495;
	t527 = t483 * t517 - t488 * t516;
	t486 = sin(pkin(7));
	t461 = t527 * t486;
	t491 = cos(pkin(7));
	t463 = t527 * t491;
	t475 = -t483 * t516 - t488 * t517;
	t484 = sin(pkin(12));
	t487 = sin(pkin(6));
	t489 = cos(pkin(12));
	t492 = cos(pkin(6));
	t452 = t492 * t461 + (t463 * t489 - t475 * t484) * t487;
	t526 = pkin(10) + r_i_i_C(3);
	t525 = pkin(3) * qJD(3);
	t485 = sin(pkin(11));
	t524 = t485 * t487;
	t523 = t485 * t492;
	t522 = t486 * t487;
	t490 = cos(pkin(11));
	t521 = t487 * t490;
	t520 = t487 * t491;
	t519 = t490 * t492;
	t493 = sin(qJ(6));
	t515 = qJD(6) * t493;
	t496 = cos(qJ(6));
	t514 = qJD(6) * t496;
	t470 = -t485 * t484 + t489 * t519;
	t456 = -t470 * t486 - t490 * t520;
	t494 = sin(qJ(5));
	t497 = cos(qJ(5));
	t507 = t498 * t483 + t495 * t488;
	t466 = t507 * t486;
	t468 = t507 * t491;
	t471 = t484 * t519 + t485 * t489;
	t476 = t495 * t483 - t498 * t488;
	t503 = -t466 * t521 + t470 * t468 - t471 * t476;
	t433 = t456 * t494 + t497 * t503;
	t511 = t456 * t497 - t494 * t503;
	t472 = -t490 * t484 - t489 * t523;
	t457 = -t472 * t486 + t485 * t520;
	t473 = -t484 * t523 + t490 * t489;
	t502 = t466 * t524 + t472 * t468 - t473 * t476;
	t435 = t457 * t494 + t497 * t502;
	t510 = t457 * t497 - t494 * t502;
	t469 = -t489 * t522 + t492 * t491;
	t501 = t492 * t466 + (t468 * t489 - t476 * t484) * t487;
	t449 = t469 * t494 + t497 * t501;
	t509 = t469 * t497 - t494 * t501;
	t506 = t496 * r_i_i_C(1) - t493 * r_i_i_C(2) + pkin(5);
	t505 = qJD(6) * (-t493 * r_i_i_C(1) - t496 * r_i_i_C(2));
	t504 = qJD(3) * t507;
	t437 = t461 * t521 - t470 * t463 + t471 * t475;
	t441 = t461 * t524 + t472 * t463 - t473 * t475;
	t500 = t526 * t494 + t506 * t497 + pkin(4);
	t499 = t497 * t505 + (-t506 * t494 + t526 * t497) * qJD(5);
	t474 = t476 * qJD(3);
	t467 = t476 * t491;
	t465 = t476 * t486;
	t464 = t491 * t504;
	t462 = t486 * t504;
	t454 = -t492 * t465 + (-t467 * t489 - t484 * t507) * t487;
	t450 = -t492 * t462 + (-t464 * t489 + t474 * t484) * t487;
	t446 = -t465 * t524 - t472 * t467 - t473 * t507;
	t443 = t465 * t521 - t470 * t467 - t471 * t507;
	t439 = -t462 * t524 - t472 * t464 + t473 * t474;
	t436 = t462 * t521 - t470 * t464 + t471 * t474;
	t431 = t509 * qJD(5) - t452 * t497;
	t429 = t510 * qJD(5) - t441 * t497;
	t427 = t511 * qJD(5) + t437 * t497;
	t1 = [0, 0, (-t441 * t493 + t502 * t514) * r_i_i_C(1) + (-t441 * t496 - t502 * t515) * r_i_i_C(2) - t441 * pkin(9) + (-t473 * t498 + (-t472 * t491 - t485 * t522) * t495) * t525 + t500 * t439 + t499 * t446, 0, t526 * t429 + t510 * t505 + t506 * (-t435 * qJD(5) + t441 * t494), (-t429 * t493 - t439 * t496) * r_i_i_C(1) + (-t429 * t496 + t439 * t493) * r_i_i_C(2) + ((-t435 * t496 + t446 * t493) * r_i_i_C(1) + (t435 * t493 + t446 * t496) * r_i_i_C(2)) * qJD(6); 0, 0, (t437 * t493 + t503 * t514) * r_i_i_C(1) + (t437 * t496 - t503 * t515) * r_i_i_C(2) + t437 * pkin(9) + (-t471 * t498 + (-t470 * t491 + t486 * t521) * t495) * t525 + t500 * t436 + t499 * t443, 0, t526 * t427 + t511 * t505 + t506 * (-t433 * qJD(5) - t437 * t494), (-t427 * t493 - t436 * t496) * r_i_i_C(1) + (-t427 * t496 + t436 * t493) * r_i_i_C(2) + ((-t433 * t496 + t443 * t493) * r_i_i_C(1) + (t433 * t493 + t443 * t496) * r_i_i_C(2)) * qJD(6); 0, 0, (-t452 * t493 + t501 * t514) * r_i_i_C(1) + (-t452 * t496 - t501 * t515) * r_i_i_C(2) - t452 * pkin(9) + (-t486 * t492 * t495 + (-t489 * t491 * t495 - t484 * t498) * t487) * t525 + t500 * t450 + t499 * t454, 0, t526 * t431 + t509 * t505 + t506 * (-t449 * qJD(5) + t452 * t494), (-t431 * t493 - t450 * t496) * r_i_i_C(1) + (-t431 * t496 + t450 * t493) * r_i_i_C(2) + ((-t449 * t496 + t454 * t493) * r_i_i_C(1) + (t449 * t493 + t454 * t496) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end