% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR13_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR13_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR13_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:10
	% EndTime: 2019-10-10 12:52:10
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:11
	% EndTime: 2019-10-10 12:52:11
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(6));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(9);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(8) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(8) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:12
	% EndTime: 2019-10-10 12:52:13
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (444->96), mult. (1331->167), div. (0->0), fcn. (1334->10), ass. (0->65)
	t374 = sin(qJ(1));
	t370 = cos(pkin(6));
	t390 = qJD(2) * t370 + qJD(1);
	t373 = sin(qJ(2));
	t408 = t374 * t373;
	t398 = t370 * t408;
	t403 = qJD(2) * t373;
	t377 = cos(qJ(2));
	t378 = cos(qJ(1));
	t405 = t378 * t377;
	t348 = -qJD(1) * t398 - t374 * t403 + t390 * t405;
	t406 = t378 * t373;
	t407 = t374 * t377;
	t359 = t370 * t406 + t407;
	t372 = sin(qJ(3));
	t376 = cos(qJ(3));
	t369 = sin(pkin(6));
	t404 = qJD(1) * t369;
	t395 = t374 * t404;
	t409 = t369 * t378;
	t397 = t376 * t409;
	t342 = (-qJD(3) * t359 + t395) * t372 - qJD(3) * t397 + t348 * t376;
	t360 = t370 * t407 + t406;
	t347 = t360 * qJD(1) + t359 * qJD(2);
	t371 = sin(qJ(4));
	t375 = cos(qJ(4));
	t426 = t342 * t371 - t347 * t375;
	t425 = -t342 * t375 - t347 * t371;
	t388 = t375 * r_i_i_C(1) - t371 * r_i_i_C(2);
	t386 = pkin(3) + t388;
	t419 = pkin(10) + r_i_i_C(3);
	t424 = (t386 * t372 - t419 * t376) * qJD(3);
	t353 = -t359 * t376 + t372 * t409;
	t396 = t370 * t405;
	t358 = -t396 + t408;
	t423 = -t353 * t371 - t358 * t375;
	t422 = t353 * t375 - t358 * t371;
	t387 = t371 * r_i_i_C(1) + t375 * r_i_i_C(2);
	t420 = t419 * t372 + t386 * t376 + pkin(2);
	t412 = t369 * t374;
	t411 = t369 * t376;
	t410 = t369 * t377;
	t402 = qJD(2) * t377;
	t401 = qJD(4) * t371;
	t400 = qJD(4) * t375;
	t399 = qJD(4) * t376;
	t394 = t378 * t404;
	t393 = t369 * t403;
	t392 = t369 * t402;
	t361 = -t398 + t405;
	t385 = -t361 * t372 + t374 * t411;
	t355 = t361 * t376 + t372 * t412;
	t357 = t370 * t372 + t373 * t411;
	t384 = -t369 * t373 * t372 + t370 * t376;
	t383 = qJD(4) * t387;
	t380 = t353 * qJD(3) - t348 * t372 + t376 * t395;
	t379 = t387 * t399 + t424;
	t350 = t384 * qJD(3) + t376 * t392;
	t346 = t359 * qJD(1) + t360 * qJD(2);
	t345 = -qJD(1) * t396 - t378 * t402 + t390 * t408;
	t340 = t385 * qJD(3) - t346 * t376 + t372 * t394;
	t339 = t355 * qJD(3) - t346 * t372 - t376 * t394;
	t338 = t340 * t375 - t345 * t371 + (-t355 * t371 + t360 * t375) * qJD(4);
	t337 = -t340 * t371 - t345 * t375 + (-t355 * t375 - t360 * t371) * qJD(4);
	t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(9) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-t378 * pkin(1) - pkin(8) * t412) * qJD(1), (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(9) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, t337 * r_i_i_C(1) - t338 * r_i_i_C(2), 0, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(9) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(8) * t409) * qJD(1), (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(9) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t397) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(9) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-t357 * qJD(3) - t372 * t392), (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:15
	% EndTime: 2019-10-10 12:52:16
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (793->117), mult. (2339->192), div. (0->0), fcn. (2413->10), ass. (0->74)
	t455 = cos(pkin(6));
	t459 = sin(qJ(1));
	t458 = sin(qJ(2));
	t501 = t459 * t458;
	t491 = t455 * t501;
	t496 = qJD(2) * t458;
	t462 = cos(qJ(2));
	t463 = cos(qJ(1));
	t498 = t463 * t462;
	t433 = -qJD(1) * t491 - t459 * t496 + (qJD(2) * t455 + qJD(1)) * t498;
	t499 = t463 * t458;
	t500 = t459 * t462;
	t444 = t455 * t499 + t500;
	t457 = sin(qJ(3));
	t461 = cos(qJ(3));
	t454 = sin(pkin(6));
	t497 = qJD(1) * t454;
	t488 = t459 * t497;
	t504 = t454 * t463;
	t490 = t461 * t504;
	t421 = (-qJD(3) * t444 + t488) * t457 - qJD(3) * t490 + t433 * t461;
	t445 = t455 * t500 + t499;
	t432 = t445 * qJD(1) + t444 * qJD(2);
	t456 = sin(qJ(4));
	t460 = cos(qJ(4));
	t438 = -t444 * t461 + t457 * t504;
	t443 = -t455 * t498 + t501;
	t521 = t438 * t460 - t443 * t456;
	t526 = t521 * qJD(4) - t421 * t456 + t432 * t460;
	t522 = t438 * t456 + t443 * t460;
	t525 = t522 * qJD(4) + t421 * t460 + t432 * t456;
	t506 = t454 * t461;
	t442 = t455 * t457 + t458 * t506;
	t502 = t456 * t462;
	t520 = -t442 * t460 + t454 * t502;
	t515 = r_i_i_C(2) + pkin(10);
	t519 = t461 * pkin(3) + t515 * t457 + pkin(2);
	t494 = qJD(3) * t462;
	t518 = (qJD(2) * t461 - qJD(4)) * t458 + t457 * t494;
	t517 = -pkin(3) * t457 + t515 * t461;
	t513 = r_i_i_C(3) + qJ(5);
	t516 = r_i_i_C(1) + pkin(4);
	t468 = t513 * t456 + t516 * t460 + pkin(3);
	t507 = t454 * t459;
	t505 = t454 * t462;
	t503 = t456 * t461;
	t495 = qJD(3) * t457;
	t493 = qJD(4) * t461;
	t487 = t463 * t497;
	t486 = t454 * t496;
	t485 = qJD(2) * t505;
	t431 = t444 * qJD(1) + t445 * qJD(2);
	t481 = t445 * t493 - t431;
	t480 = t443 * t493 + t433;
	t472 = t491 - t498;
	t440 = t457 * t507 - t461 * t472;
	t477 = t440 * t460 + t445 * t456;
	t476 = -t440 * t456 + t445 * t460;
	t475 = (qJD(2) - t493) * t462;
	t474 = t457 * t472 + t459 * t506;
	t473 = -t454 * t458 * t457 + t455 * t461;
	t469 = qJD(3) * t517;
	t430 = t443 * qJD(1) + t472 * qJD(2);
	t467 = -qJD(4) * t472 + t430 * t461 + t445 * t495;
	t466 = qJD(4) * t444 - t432 * t461 + t443 * t495;
	t465 = qJD(5) * t456 + (-t516 * t456 + t513 * t460) * qJD(4);
	t464 = t438 * qJD(3) - t433 * t457 + t461 * t488;
	t435 = t473 * qJD(3) + t461 * t485;
	t424 = -t520 * qJD(4) + t435 * t456 - t460 * t486;
	t419 = t474 * qJD(3) - t431 * t461 + t457 * t487;
	t418 = t440 * qJD(3) - t431 * t457 - t461 * t487;
	t409 = t476 * qJD(4) + t419 * t460 - t430 * t456;
	t408 = t477 * qJD(4) + t419 * t456 + t430 * t460;
	t1 = [t522 * qJD(5) - t421 * pkin(3) - t433 * pkin(2) - t432 * pkin(9) + t515 * t464 - t516 * t525 + t513 * t526 + (-t463 * pkin(1) - pkin(8) * t507) * qJD(1), -(t445 * t503 - t460 * t472) * qJD(5) - t431 * pkin(9) + t516 * (t481 * t456 + t467 * t460) + t513 * (t467 * t456 - t481 * t460) - t445 * t469 + t519 * t430, -t468 * t418 + t515 * t419 + t465 * t474, t477 * qJD(5) - t516 * t408 + t513 * t409, t408, 0; -t476 * qJD(5) + t419 * pkin(3) - t431 * pkin(2) - t430 * pkin(9) + t515 * t418 + t516 * t409 + t513 * t408 + (-t459 * pkin(1) + pkin(8) * t504) * qJD(1), -(t443 * t503 + t444 * t460) * qJD(5) + t433 * pkin(9) + t516 * (t480 * t456 + t466 * t460) + t513 * (t466 * t456 - t480 * t460) - t443 * t469 - t519 * t432, t515 * t421 + t465 * (-t444 * t457 - t490) + t468 * t464, -t521 * qJD(5) + t513 * t525 + t516 * t526, -t526, 0; 0, (t516 * (t456 * t475 - t518 * t460) - t513 * (t518 * t456 + t460 * t475) - (t458 * t460 - t461 * t502) * qJD(5) + t517 * t494 + (t462 * pkin(9) - t458 * t519) * qJD(2)) * t454, t515 * t435 + t465 * t473 + t468 * (-t442 * qJD(3) - t457 * t485), -t520 * qJD(5) + t513 * (t456 * t486 + t435 * t460 + (-t442 * t456 - t460 * t505) * qJD(4)) - t516 * t424, t424, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:52:14
	% EndTime: 2019-10-10 12:52:16
	% DurationCPUTime: 1.93s
	% Computational Cost: add. (1624->179), mult. (4756->297), div. (0->0), fcn. (5097->12), ass. (0->100)
	t503 = sin(qJ(2));
	t504 = sin(qJ(1));
	t508 = cos(qJ(2));
	t561 = cos(pkin(6));
	t563 = cos(qJ(1));
	t526 = t561 * t563;
	t484 = t503 * t526 + t504 * t508;
	t502 = sin(qJ(3));
	t507 = cos(qJ(3));
	t499 = sin(pkin(6));
	t537 = t499 * t563;
	t471 = t484 * t507 - t502 * t537;
	t483 = t504 * t503 - t508 * t526;
	t501 = sin(qJ(4));
	t506 = cos(qJ(4));
	t450 = t471 * t501 - t483 * t506;
	t530 = t504 * t561;
	t485 = t563 * t503 + t508 * t530;
	t464 = t485 * qJD(1) + t484 * qJD(2);
	t574 = -t450 * qJD(4) + t464 * t501;
	t451 = t471 * t506 + t483 * t501;
	t573 = t451 * qJD(4) - t464 * t506;
	t500 = sin(qJ(6));
	t505 = cos(qJ(6));
	t572 = t450 * t505 - t451 * t500;
	t571 = t450 * t500 + t451 * t505;
	t570 = (t571 * r_i_i_C(1) + t572 * r_i_i_C(2)) * qJD(6);
	t551 = t503 * t507;
	t482 = t499 * t551 + t561 * t502;
	t548 = t506 * t508;
	t468 = t482 * t501 + t499 * t548;
	t552 = t501 * t508;
	t558 = t482 * t506;
	t469 = -t499 * t552 + t558;
	t569 = ((t468 * t500 + t469 * t505) * r_i_i_C(1) + (t468 * t505 - t469 * t500) * r_i_i_C(2)) * qJD(6);
	t540 = -r_i_i_C(3) - pkin(11) + pkin(10);
	t568 = t507 * pkin(3) + t540 * t502 + pkin(2);
	t567 = -pkin(3) * t502 + t540 * t507;
	t564 = pkin(4) + pkin(5);
	t516 = r_i_i_C(1) * t505 - r_i_i_C(2) * t500 + t564;
	t517 = r_i_i_C(1) * t500 + r_i_i_C(2) * t505 + qJ(5);
	t511 = t517 * t501 + t516 * t506 + pkin(3);
	t555 = t499 * t504;
	t554 = t501 * t503;
	t553 = t501 * t507;
	t550 = t504 * t507;
	t549 = t506 * t507;
	t547 = t507 * t508;
	t546 = qJD(1) * t504;
	t545 = qJD(2) * t503;
	t544 = qJD(3) * t502;
	t543 = qJD(3) * t508;
	t542 = qJD(4) * t499;
	t541 = qJD(4) * t507;
	t539 = t499 * t503 * t502;
	t536 = t499 * t546;
	t535 = t499 * t545;
	t534 = qJD(2) * t499 * t508;
	t533 = t502 * t543;
	t532 = t508 * t542;
	t531 = t563 * qJD(1);
	t527 = t503 * t530;
	t465 = -qJD(1) * t527 - t504 * t545 + (qJD(2) * t526 + t531) * t508;
	t528 = t507 * t537;
	t529 = qJD(3) * t528 - t465 * t507;
	t463 = t484 * qJD(1) + t485 * qJD(2);
	t525 = t485 * t541 - t463;
	t524 = t483 * t541 + t465;
	t486 = t563 * t508 - t527;
	t474 = t486 * t507 + t502 * t555;
	t454 = t474 * t501 - t485 * t506;
	t455 = t474 * t506 + t485 * t501;
	t521 = t454 * t505 - t455 * t500;
	t520 = t454 * t500 + t455 * t505;
	t514 = qJD(3) * t567;
	t462 = qJD(2) * t527 + t503 * t546 + (-qJD(1) * t526 - t563 * qJD(2)) * t508;
	t513 = qJD(4) * t486 + t462 * t507 + t485 * t544;
	t512 = qJD(4) * t484 - t464 * t507 + t483 * t544;
	t510 = -t471 * qJD(3) - t465 * t502 + t507 * t536;
	t509 = t501 * qJD(5) + ((-t500 * t506 + t501 * t505) * r_i_i_C(1) + (-t500 * t501 - t505 * t506) * r_i_i_C(2)) * qJD(6) + (-t516 * t501 + t517 * t506) * qJD(4);
	t476 = (t506 * t547 + t554) * t499;
	t475 = (t501 * t547 - t503 * t506) * t499;
	t467 = -qJD(3) * t539 + (t561 * qJD(3) + t534) * t507;
	t459 = -t485 * t549 + t486 * t501;
	t458 = -t485 * t553 - t486 * t506;
	t457 = -t483 * t549 + t484 * t501;
	t456 = -t483 * t553 - t484 * t506;
	t447 = -t468 * qJD(4) + t467 * t506 + t501 * t535;
	t446 = qJD(4) * t558 - t506 * t535 + (t467 - t532) * t501;
	t445 = (qJD(3) * t484 - t536) * t502 + t529;
	t443 = -t484 * t544 + t502 * t536 - t529;
	t441 = -t463 * t507 - t486 * t544 + (qJD(3) * t550 + t502 * t531) * t499;
	t440 = -qJD(1) * t528 + t474 * qJD(3) - t463 * t502;
	t433 = t443 * t506 + t574;
	t432 = t443 * t501 + t573;
	t431 = -t454 * qJD(4) + t441 * t506 - t462 * t501;
	t430 = t455 * qJD(4) + t441 * t501 + t462 * t506;
	t429 = t521 * qJD(6) + t430 * t500 + t431 * t505;
	t428 = -t520 * qJD(6) + t430 * t505 - t431 * t500;
	t1 = [-t465 * pkin(2) + t445 * pkin(3) - t464 * pkin(9) - t450 * qJD(5) + t517 * (t445 * t501 - t573) + t516 * (t445 * t506 - t574) + (-t572 * r_i_i_C(1) + t571 * r_i_i_C(2)) * qJD(6) + (-t563 * pkin(1) - pkin(8) * t555) * qJD(1) + t540 * t510, -t463 * pkin(9) + t458 * qJD(5) + t517 * (t513 * t501 - t525 * t506) + t516 * (t525 * t501 + t513 * t506) + ((t458 * t505 - t459 * t500) * r_i_i_C(1) + (-t458 * t500 - t459 * t505) * r_i_i_C(2)) * qJD(6) - t485 * t514 + t568 * t462, t540 * t441 - t511 * t440 + t509 * (-t486 * t502 + t499 * t550), t455 * qJD(5) + t517 * t431 - t516 * t430 + (t520 * r_i_i_C(1) + t521 * r_i_i_C(2)) * qJD(6), t430, r_i_i_C(1) * t428 - r_i_i_C(2) * t429; -t463 * pkin(2) + t441 * pkin(3) - t462 * pkin(9) + t429 * r_i_i_C(1) + t428 * r_i_i_C(2) + t430 * qJ(5) + t454 * qJD(5) + t564 * t431 + (-pkin(1) * t504 + pkin(8) * t537) * qJD(1) + t540 * t440, t465 * pkin(9) + t456 * qJD(5) + t517 * (t512 * t501 - t524 * t506) + t516 * (t524 * t501 + t512 * t506) + ((t456 * t505 - t457 * t500) * r_i_i_C(1) + (-t456 * t500 - t457 * t505) * r_i_i_C(2)) * qJD(6) - t483 * t514 - t568 * t464, t540 * t443 + t511 * t510 + t509 * (-t484 * t502 - t528), t451 * qJD(5) - t516 * t432 + t517 * t433 + t570, t432, (t432 * t505 - t433 * t500) * r_i_i_C(1) + (-t432 * t500 - t433 * t505) * r_i_i_C(2) - t570; 0, t475 * qJD(5) - t517 * (-t532 * t549 - t542 * t554) + ((t475 * t505 - t476 * t500) * r_i_i_C(1) + (-t475 * t500 - t476 * t505) * r_i_i_C(2)) * qJD(6) + (-t517 * (t501 * t533 + (t501 * t551 + t548) * qJD(2)) + t516 * ((qJD(2) - t541) * t552 + (-t533 + (-qJD(2) * t507 + qJD(4)) * t503) * t506) + t567 * t543 + (pkin(9) * t508 - t503 * t568) * qJD(2)) * t499, t540 * t467 + t511 * (-t482 * qJD(3) - t502 * t534) + t509 * (t561 * t507 - t539), t469 * qJD(5) - t516 * t446 + t517 * t447 + t569, t446, (t446 * t505 - t447 * t500) * r_i_i_C(1) + (-t446 * t500 - t447 * t505) * r_i_i_C(2) - t569;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end