% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
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
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 12:31:17
	% EndTime: 2019-10-10 12:31:17
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
	% StartTime: 2019-10-10 12:31:18
	% EndTime: 2019-10-10 12:31:19
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (444->96), mult. (1331->167), div. (0->0), fcn. (1334->10), ass. (0->65)
	t374 = sin(qJ(1));
	t370 = cos(pkin(6));
	t390 = qJD(2) * t370 + qJD(1);
	t373 = sin(qJ(2));
	t408 = t373 * t374;
	t397 = t370 * t408;
	t403 = qJD(2) * t373;
	t377 = cos(qJ(2));
	t378 = cos(qJ(1));
	t405 = t377 * t378;
	t348 = -qJD(1) * t397 - t374 * t403 + t390 * t405;
	t406 = t374 * t377;
	t407 = t373 * t378;
	t359 = t370 * t407 + t406;
	t372 = sin(qJ(3));
	t376 = cos(qJ(3));
	t369 = sin(pkin(6));
	t404 = qJD(1) * t369;
	t395 = t374 * t404;
	t409 = t369 * t378;
	t398 = t376 * t409;
	t342 = (-qJD(3) * t359 + t395) * t372 - qJD(3) * t398 + t348 * t376;
	t360 = t370 * t406 + t407;
	t347 = t360 * qJD(1) + t359 * qJD(2);
	t371 = sin(qJ(4));
	t375 = cos(qJ(4));
	t426 = t342 * t371 - t347 * t375;
	t425 = -t342 * t375 - t347 * t371;
	t388 = r_i_i_C(1) * t375 - r_i_i_C(2) * t371;
	t386 = pkin(3) + t388;
	t419 = pkin(10) + r_i_i_C(3);
	t424 = (t386 * t372 - t419 * t376) * qJD(3);
	t353 = -t359 * t376 + t372 * t409;
	t396 = t370 * t405;
	t358 = -t396 + t408;
	t423 = -t353 * t371 - t358 * t375;
	t422 = t353 * t375 - t358 * t371;
	t387 = r_i_i_C(1) * t371 + r_i_i_C(2) * t375;
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
	t361 = -t397 + t405;
	t385 = -t361 * t372 + t374 * t411;
	t355 = t361 * t376 + t372 * t412;
	t357 = t370 * t372 + t373 * t411;
	t384 = -t369 * t372 * t373 + t370 * t376;
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
	t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(9) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-pkin(1) * t378 - pkin(8) * t412) * qJD(1), (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(9) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, r_i_i_C(1) * t337 - r_i_i_C(2) * t338, 0, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(9) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(8) * t409) * qJD(1), (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(9) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t398) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(9) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-t357 * qJD(3) - t372 * t392), (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:21
	% EndTime: 2019-10-10 12:31:22
	% DurationCPUTime: 1.04s
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
	% StartTime: 2019-10-10 12:31:18
	% EndTime: 2019-10-10 12:31:19
	% DurationCPUTime: 1.14s
	% Computational Cost: add. (1012->129), mult. (2977->203), div. (0->0), fcn. (3075->10), ass. (0->78)
	t416 = sin(qJ(2));
	t417 = sin(qJ(1));
	t420 = cos(qJ(2));
	t421 = cos(qJ(1));
	t471 = cos(pkin(6));
	t440 = t421 * t471;
	t399 = t416 * t440 + t417 * t420;
	t441 = t417 * t471;
	t400 = t421 * t416 + t420 * t441;
	t385 = t400 * qJD(1) + t399 * qJD(2);
	t418 = cos(qJ(4));
	t413 = sin(pkin(6));
	t419 = cos(qJ(3));
	t415 = sin(qJ(3));
	t461 = t415 * t421;
	t391 = -t399 * t419 + t413 * t461;
	t398 = t416 * t417 - t420 * t440;
	t414 = sin(qJ(4));
	t481 = t391 * t418 - t398 * t414;
	t484 = -t481 * qJD(4) - t385 * t418;
	t482 = t391 * t414 + t398 * t418;
	t483 = t482 * qJD(4) + t385 * t414;
	t451 = pkin(10) - r_i_i_C(3) - qJ(6);
	t454 = qJD(5) * t414;
	t480 = (-pkin(3) * t415 + t451 * t419) * qJD(3) - t415 * qJD(6) + t419 * t454;
	t452 = pkin(4) + pkin(5) + r_i_i_C(1);
	t472 = r_i_i_C(2) + qJ(5);
	t479 = -t454 + (t452 * t414 - t472 * t418) * qJD(4);
	t464 = t413 * t419;
	t397 = t471 * t415 + t416 * t464;
	t463 = t413 * t420;
	t478 = -t397 * t418 + t414 * t463;
	t477 = t419 * pkin(3) + t451 * t415 + pkin(2);
	t456 = qJD(3) * t415;
	t476 = (qJD(2) * t419 - qJD(4)) * t416 + t420 * t456;
	t426 = t472 * t414 + t452 * t418 + pkin(3);
	t465 = t413 * t417;
	t462 = t413 * t421;
	t460 = t417 * t419;
	t459 = t421 * t420;
	t458 = qJD(1) * t417;
	t457 = qJD(2) * t416;
	t455 = qJD(4) * t419;
	t453 = qJD(5) * t418;
	t449 = t413 * t416 * t415;
	t448 = t419 * t462;
	t446 = t413 * t458;
	t445 = qJD(1) * t464;
	t444 = t413 * t457;
	t443 = qJD(2) * t463;
	t437 = t416 * t441;
	t386 = -qJD(1) * t437 - t417 * t457 + (qJD(2) * t471 + qJD(1)) * t459;
	t439 = qJD(3) * t448 - t386 * t419;
	t384 = t399 * qJD(1) + t400 * qJD(2);
	t436 = t400 * t455 - t384;
	t435 = t398 * t455 + t386;
	t401 = -t437 + t459;
	t393 = t401 * t419 + t415 * t465;
	t432 = t393 * t418 + t400 * t414;
	t431 = -t393 * t414 + t400 * t418;
	t430 = (qJD(2) - t455) * t420;
	t429 = -t399 * t415 - t448;
	t383 = qJD(2) * t437 + t416 * t458 + (-t471 * qJD(1) - qJD(2)) * t459;
	t425 = qJD(4) * t401 + t383 * t419 + t400 * t456;
	t424 = qJD(4) * t399 - t385 * t419 + t398 * t456;
	t373 = t391 * qJD(3) - t386 * t415 + t417 * t445;
	t392 = t401 * t415 - t413 * t460;
	t388 = -qJD(3) * t449 + (t471 * qJD(3) + t443) * t419;
	t387 = -t397 * qJD(3) - t415 * t443;
	t377 = -t478 * qJD(4) + t388 * t414 - t418 * t444;
	t376 = (qJD(3) * t399 - t446) * t415 + t439;
	t374 = -t399 * t456 + t415 * t446 - t439;
	t372 = -t384 * t419 - t401 * t456 + (qJD(1) * t461 + qJD(3) * t460) * t413;
	t371 = t393 * qJD(3) - t384 * t415 - t421 * t445;
	t363 = t374 * t414 + t484;
	t362 = t431 * qJD(4) + t372 * t418 - t383 * t414;
	t361 = t432 * qJD(4) + t372 * t414 + t383 * t418;
	t1 = [-t429 * qJD(6) + t482 * qJD(5) + t376 * pkin(3) - t386 * pkin(2) - t385 * pkin(9) + t472 * (t376 * t414 - t484) + (-pkin(1) * t421 - pkin(8) * t465) * qJD(1) + t451 * t373 + t452 * (t376 * t418 - t483), -t401 * t453 - t384 * pkin(9) + t472 * (t425 * t414 - t436 * t418) + t452 * (t436 * t414 + t425 * t418) - t480 * t400 + t477 * t383, -t393 * qJD(6) - t426 * t371 + t451 * t372 + t479 * t392, t432 * qJD(5) - t452 * t361 + t472 * t362, t361, -t371; -t392 * qJD(6) - t431 * qJD(5) + t372 * pkin(3) - t384 * pkin(2) - t383 * pkin(9) + t472 * t361 + (-pkin(1) * t417 + pkin(8) * t462) * qJD(1) + t451 * t371 + t452 * t362, -t399 * t453 + t386 * pkin(9) + t472 * (t424 * t414 - t435 * t418) + t452 * (t435 * t414 + t424 * t418) - t480 * t398 - t477 * t385, qJD(6) * t391 + t426 * t373 + t451 * t374 - t429 * t479, -t481 * qJD(5) + t472 * (t374 * t418 + t483) - t452 * t363, t363, t373; 0, (-t472 * (t476 * t414 + t418 * t430) + t452 * (t414 * t430 - t476 * t418) - t416 * t453 + t480 * t420 + (t420 * pkin(9) - t416 * t477) * qJD(2)) * t413, -qJD(6) * t397 + t451 * t388 - t479 * (t471 * t419 - t449) + t426 * t387, -t478 * qJD(5) + t472 * (t414 * t444 + t388 * t418 + (-t397 * t414 - t418 * t463) * qJD(4)) - t452 * t377, t377, t387;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end