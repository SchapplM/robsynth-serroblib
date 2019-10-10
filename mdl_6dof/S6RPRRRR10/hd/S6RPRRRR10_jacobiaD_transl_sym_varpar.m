% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
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
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(6));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(6));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(13));
	t113 = sin(pkin(13));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:31
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(7));
	t278 = r_i_i_C(3) + pkin(9);
	t279 = t278 * t248 + qJ(2);
	t250 = sin(qJ(3));
	t277 = t250 * r_i_i_C(2);
	t246 = sin(pkin(6));
	t251 = sin(qJ(1));
	t276 = t246 * t251;
	t253 = cos(qJ(1));
	t275 = t246 * t253;
	t274 = t248 * t250;
	t252 = cos(qJ(3));
	t273 = t248 * t252;
	t244 = sin(pkin(13));
	t272 = t251 * t244;
	t247 = cos(pkin(13));
	t271 = t251 * t247;
	t270 = t253 * t244;
	t269 = t253 * t247;
	t268 = qJD(1) * t251;
	t267 = qJD(1) * t253;
	t245 = sin(pkin(7));
	t266 = qJD(3) * t245;
	t265 = t278 * t245;
	t249 = cos(pkin(6));
	t264 = t249 * t272;
	t263 = t246 * t268;
	t262 = t246 * t267;
	t261 = t266 * t275;
	t260 = r_i_i_C(1) * t250 + r_i_i_C(2) * t252;
	t236 = -qJD(1) * t264 + t247 * t267;
	t237 = -t249 * t269 + t272;
	t259 = qJD(3) * t237 * t248 - t236;
	t258 = t260 * t245;
	t256 = t249 * t271 + t270;
	t257 = t245 * t276 - t248 * t256;
	t238 = t249 * t270 + t271;
	t233 = t237 * qJD(1);
	t255 = t233 * t248 + t245 * t262;
	t235 = t256 * qJD(1);
	t254 = -qJD(3) * t238 - t235 * t248 + t245 * t263;
	t241 = t252 * t261;
	t240 = -t264 + t269;
	t234 = t238 * qJD(1);
	t232 = -t234 * t252 + t255 * t250 + (-t240 * t250 + t257 * t252) * qJD(3);
	t231 = t234 * t250 + t255 * t252 + (-t240 * t252 - t257 * t250) * qJD(3);
	t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0, 0; 0, 0, (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:32
	% EndTime: 2019-10-10 09:11:33
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (393->78), mult. (1301->142), div. (0->0), fcn. (1376->12), ass. (0->64)
	t408 = sin(qJ(3));
	t402 = sin(pkin(7));
	t403 = sin(pkin(6));
	t412 = cos(qJ(1));
	t442 = t403 * t412;
	t431 = t402 * t442;
	t405 = cos(pkin(7));
	t404 = cos(pkin(13));
	t406 = cos(pkin(6));
	t440 = t406 * t412;
	t401 = sin(pkin(13));
	t409 = sin(qJ(1));
	t445 = t401 * t409;
	t421 = t404 * t440 - t445;
	t454 = t421 * t405;
	t424 = -t454 + t431;
	t463 = t424 * t408;
	t432 = t406 * t445;
	t438 = qJD(1) * t412;
	t389 = -qJD(1) * t432 + t404 * t438;
	t439 = t409 * t404;
	t392 = t401 * t440 + t439;
	t411 = cos(qJ(3));
	t420 = t401 * t412 + t406 * t439;
	t388 = t420 * qJD(1);
	t443 = t403 * t409;
	t429 = qJD(1) * t443;
	t418 = -t388 * t405 + t402 * t429;
	t428 = t411 * t431;
	t368 = -qJD(3) * t428 + (-qJD(3) * t392 + t418) * t408 - (-qJD(3) * t454 - t389) * t411;
	t449 = t388 * t402;
	t380 = t405 * t429 + t449;
	t407 = sin(qJ(4));
	t410 = cos(qJ(4));
	t462 = t368 * t407 - t380 * t410;
	t461 = -t368 * t410 - t380 * t407;
	t373 = -t392 * t411 + t463;
	t453 = t421 * t402 + t405 * t442;
	t460 = -t373 * t407 + t453 * t410;
	t459 = t373 * t410 + t453 * t407;
	t441 = t404 * t405;
	t444 = t402 * t406;
	t452 = (-t401 * t408 + t411 * t441) * t403 + t411 * t444;
	t382 = (t401 * t411 + t408 * t441) * t403 + t408 * t444;
	t450 = r_i_i_C(3) + pkin(10);
	t436 = t403 * qJD(2);
	t427 = t407 * r_i_i_C(1) + t410 * r_i_i_C(2);
	t425 = t410 * r_i_i_C(1) - t407 * r_i_i_C(2) + pkin(3);
	t423 = t402 * t443 - t405 * t420;
	t419 = qJD(4) * t427;
	t394 = t404 * t412 - t432;
	t415 = -t394 * t408 + t423 * t411;
	t375 = t394 * t411 + t423 * t408;
	t413 = t373 * qJD(3) - t389 * t408 + t418 * t411;
	t390 = -t402 * t403 * t404 + t405 * t406;
	t387 = t392 * qJD(1);
	t385 = t402 * t420 + t405 * t443;
	t378 = t453 * qJD(1);
	t376 = t452 * qJD(3);
	t366 = qJD(1) * t463 + t415 * qJD(3) - t387 * t411;
	t365 = t375 * qJD(3) - t387 * t408 + (t411 * t454 - t428) * qJD(1);
	t364 = t366 * t410 + t378 * t407 + (-t375 * t407 + t385 * t410) * qJD(4);
	t363 = -t366 * t407 + t378 * t410 + (-t375 * t410 - t385 * t407) * qJD(4);
	t1 = [t461 * r_i_i_C(1) + t462 * r_i_i_C(2) - t368 * pkin(3) - t389 * pkin(2) - pkin(9) * t449 + t412 * t436 + t450 * t413 + (t460 * r_i_i_C(1) - t459 * r_i_i_C(2)) * qJD(4) + (-t412 * pkin(1) + (-pkin(9) * t405 - qJ(2)) * t443) * qJD(1), t403 * t438, -t425 * t365 + t450 * t366 - t415 * t419, r_i_i_C(1) * t363 - r_i_i_C(2) * t364, 0, 0; t409 * t436 - t387 * pkin(2) + t366 * pkin(3) + t364 * r_i_i_C(1) + t363 * r_i_i_C(2) + t450 * t365 + (-t409 * pkin(1) + pkin(9) * t453 + qJ(2) * t442) * qJD(1), t429, t450 * t368 - (-t392 * t408 - t424 * t411) * t419 + t425 * t413, -t462 * r_i_i_C(1) + t461 * r_i_i_C(2) + (t459 * r_i_i_C(1) + t460 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t425 * t382 * qJD(3) + t450 * t376 - t452 * t419, -t427 * t376 + ((-t382 * t410 - t390 * t407) * r_i_i_C(1) + (t382 * t407 - t390 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:32
	% EndTime: 2019-10-10 09:11:33
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (695->95), mult. (1898->163), div. (0->0), fcn. (2026->14), ass. (0->74)
	t450 = sin(qJ(3));
	t444 = sin(pkin(7));
	t445 = sin(pkin(6));
	t454 = cos(qJ(1));
	t496 = t445 * t454;
	t480 = t444 * t496;
	t447 = cos(pkin(7));
	t448 = cos(pkin(6));
	t446 = cos(pkin(13));
	t491 = t454 * t446;
	t443 = sin(pkin(13));
	t451 = sin(qJ(1));
	t494 = t451 * t443;
	t464 = t448 * t491 - t494;
	t508 = t464 * t447;
	t467 = -t508 + t480;
	t512 = t467 * t450;
	t481 = t448 * t494;
	t487 = qJD(1) * t454;
	t426 = -qJD(1) * t481 + t446 * t487;
	t492 = t454 * t443;
	t493 = t451 * t446;
	t429 = t448 * t492 + t493;
	t453 = cos(qJ(3));
	t463 = t448 * t493 + t492;
	t425 = t463 * qJD(1);
	t497 = t445 * t451;
	t478 = qJD(1) * t497;
	t462 = -t425 * t447 + t444 * t478;
	t470 = t453 * t480;
	t403 = -qJD(3) * t470 + (-qJD(3) * t429 + t462) * t450 - (-qJD(3) * t508 - t426) * t453;
	t441 = qJD(4) + qJD(5);
	t507 = t464 * t444 + t447 * t496;
	t511 = t507 * t441 - t403;
	t410 = -t429 * t453 + t512;
	t502 = t425 * t444;
	t417 = t447 * t478 + t502;
	t510 = -t410 * t441 - t417;
	t495 = t446 * t447;
	t498 = t444 * t448;
	t506 = (-t443 * t450 + t453 * t495) * t445 + t453 * t498;
	t419 = (t443 * t453 + t450 * t495) * t445 + t450 * t498;
	t504 = r_i_i_C(3) + pkin(11) + pkin(10);
	t503 = t419 * t441;
	t442 = qJ(4) + qJ(5);
	t439 = sin(t442);
	t440 = cos(t442);
	t431 = -t481 + t491;
	t466 = t444 * t497 - t447 * t463;
	t412 = t431 * t453 + t466 * t450;
	t415 = t507 * qJD(1);
	t472 = -t412 * t441 + t415;
	t424 = t429 * qJD(1);
	t458 = -t431 * t450 + t466 * t453;
	t401 = qJD(1) * t512 + t458 * qJD(3) - t424 * t453;
	t422 = t444 * t463 + t447 * t497;
	t477 = t422 * t441 + t401;
	t398 = -t477 * t439 + t472 * t440;
	t399 = t472 * t439 + t477 * t440;
	t490 = t398 * r_i_i_C(1) - t399 * r_i_i_C(2);
	t489 = (t439 * t511 - t440 * t510) * r_i_i_C(1) + (t439 * t510 + t440 * t511) * r_i_i_C(2);
	t413 = t506 * qJD(3);
	t427 = -t445 * t446 * t444 + t448 * t447;
	t471 = -t427 * t441 - t413;
	t488 = (t471 * t439 - t440 * t503) * r_i_i_C(1) + (t439 * t503 + t471 * t440) * r_i_i_C(2);
	t485 = t445 * qJD(2);
	t452 = cos(qJ(4));
	t438 = t452 * pkin(4) + pkin(3);
	t468 = t440 * r_i_i_C(1) - t439 * r_i_i_C(2) + t438;
	t449 = sin(qJ(4));
	t459 = -qJD(4) * t449 * pkin(4) + (-t439 * r_i_i_C(1) - t440 * r_i_i_C(2)) * t441;
	t456 = t410 * qJD(3) - t426 * t450 + t462 * t453;
	t400 = t412 * qJD(3) - t424 * t450 + (t453 * t508 - t470) * qJD(1);
	t1 = [-pkin(9) * t502 + t454 * t485 - t426 * pkin(2) - t403 * t438 + (r_i_i_C(1) * t511 + r_i_i_C(2) * t510) * t440 + (r_i_i_C(1) * t510 - r_i_i_C(2) * t511) * t439 + t504 * t456 + (-t454 * pkin(1) + (-pkin(9) * t447 - qJ(2)) * t497) * qJD(1) + (-t417 * t449 + (-t410 * t449 + t452 * t507) * qJD(4)) * pkin(4), t445 * t487, -t468 * t400 + t504 * t401 + t459 * t458, (-t401 * t449 + t415 * t452 + (-t412 * t452 - t422 * t449) * qJD(4)) * pkin(4) + t490, t490, 0; t451 * t485 - t424 * pkin(2) + t399 * r_i_i_C(1) + t398 * r_i_i_C(2) + t401 * t438 + t504 * t400 + (t415 * t449 + (-t412 * t449 + t422 * t452) * qJD(4)) * pkin(4) + (-t451 * pkin(1) + pkin(9) * t507 + qJ(2) * t496) * qJD(1), t478, t504 * t403 + t459 * (-t429 * t450 - t467 * t453) + t468 * t456, (-t403 * t449 + t417 * t452 + (t410 * t452 + t449 * t507) * qJD(4)) * pkin(4) + t489, t489, 0; 0, 0, -t468 * t419 * qJD(3) + t504 * t413 + t459 * t506, (-t413 * t449 + (-t419 * t452 - t427 * t449) * qJD(4)) * pkin(4) + t488, t488, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:38
	% EndTime: 2019-10-10 09:11:39
	% DurationCPUTime: 1.96s
	% Computational Cost: add. (1920->158), mult. (4975->262), div. (0->0), fcn. (5524->16), ass. (0->105)
	t639 = qJ(4) + qJ(5);
	t637 = cos(t639);
	t646 = cos(qJ(1));
	t713 = cos(pkin(13));
	t715 = cos(pkin(6));
	t681 = t715 * t713;
	t711 = sin(pkin(13));
	t719 = sin(qJ(1));
	t662 = t646 * t711 + t719 * t681;
	t616 = t662 * qJD(1);
	t679 = t715 * t711;
	t621 = t646 * t713 - t719 * t679;
	t617 = t621 * qJD(1);
	t643 = sin(qJ(3));
	t660 = -t646 * t679 - t719 * t713;
	t640 = sin(pkin(6));
	t712 = sin(pkin(7));
	t694 = t640 * t712;
	t673 = t719 * t694;
	t714 = cos(pkin(7));
	t720 = cos(qJ(3));
	t661 = t646 * t681 - t719 * t711;
	t684 = t646 * t694;
	t672 = t720 * t684;
	t683 = t714 * t720;
	t727 = t661 * t683 - t672;
	t577 = (qJD(1) * t673 + qJD(3) * t660 - t714 * t616) * t643 + t617 * t720 + t727 * qJD(3);
	t638 = qJD(4) + qJD(5);
	t706 = t640 * t646;
	t725 = t661 * t712 + t714 * t706;
	t690 = -t638 * t725 + t577;
	t739 = t690 * t637;
	t691 = t714 * t661;
	t699 = t660 * t720;
	t597 = (t684 - t691) * t643 + t699;
	t636 = sin(t639);
	t587 = t597 * t637 + t636 * t725;
	t594 = -t643 * t660 - t727;
	t641 = sin(qJ(6));
	t644 = cos(qJ(6));
	t738 = -t587 * t641 - t594 * t644;
	t737 = t587 * t644 - t594 * t641;
	t682 = t714 * t719;
	t674 = t640 * t682;
	t695 = t616 * t712;
	t605 = qJD(1) * t674 + t695;
	t688 = t597 * t638 + t605;
	t736 = -t690 * t636 + t688 * t637;
	t668 = t720 * t673;
	t676 = t643 * t684;
	t578 = qJD(1) * t668 - t616 * t683 - t617 * t643 + (-t643 * t691 + t676 + t699) * qJD(3);
	t735 = t578 * t641;
	t734 = t578 * t644;
	t721 = r_i_i_C(3) + pkin(12);
	t732 = -r_i_i_C(1) * t644 - pkin(5);
	t704 = qJD(6) * t641;
	t702 = r_i_i_C(1) * t704;
	t703 = qJD(6) * t644;
	t726 = -t703 * r_i_i_C(2) - t702;
	t677 = -r_i_i_C(2) * t641 - t732;
	t656 = t662 * t714;
	t599 = t621 * t720 + (-t656 + t673) * t643;
	t650 = qJD(1) * t725;
	t717 = pkin(4) * qJD(4);
	t724 = pkin(4) * t650 - t599 * t717;
	t723 = -t621 * t643 - t720 * t656 + t668;
	t678 = t714 * t713;
	t680 = t715 * t712;
	t606 = -t720 * t680 + (t643 * t711 - t678 * t720) * t640;
	t615 = t660 * qJD(1);
	t654 = qJD(1) * t691;
	t575 = qJD(1) * t676 + t723 * qJD(3) + t615 * t720 - t643 * t654;
	t610 = t662 * t712 + t674;
	t707 = t637 * t638;
	t568 = t599 * t707 - t637 * t650 + (t610 * t638 + t575) * t636;
	t708 = t636 * t638;
	t569 = t575 * t637 - t599 * t708 + t610 * t707 + t636 * t650;
	t588 = -t599 * t636 + t610 * t637;
	t722 = (t568 * t641 - t588 * t703) * r_i_i_C(2) - t588 * t702 + t721 * t569 + t732 * t568;
	t705 = qJD(2) * t640;
	t642 = sin(qJ(4));
	t701 = t642 * t717;
	t698 = qJD(1) * t706;
	t696 = qJD(1) * t719;
	t602 = t606 * qJD(3);
	t618 = -t713 * t694 + t715 * t714;
	t687 = t618 * t638 - t602;
	t645 = cos(qJ(4));
	t635 = pkin(4) * t645 + pkin(3);
	t659 = -t721 * t636 - t677 * t637 - t635;
	t571 = t597 * t708 + t605 * t636 + t739;
	t652 = t726 * (t597 * t636 - t637 * t725) + t721 * t571 + t677 * t736;
	t607 = t643 * t680 + (t643 * t678 + t720 * t711) * t640;
	t584 = -t607 * t708 + t687 * t637;
	t651 = t726 * (-t607 * t636 + t618 * t637) + t721 * t584 + t677 * (-t607 * t707 - t687 * t636);
	t648 = t701 + (r_i_i_C(1) * t641 + r_i_i_C(2) * t644) * t637 * qJD(6) + (t677 * t636 - t721 * t637) * t638;
	t647 = -pkin(11) - pkin(10);
	t603 = t607 * qJD(3);
	t593 = t607 * t637 + t618 * t636;
	t589 = t599 * t637 + t610 * t636;
	t574 = -qJD(1) * t672 + t599 * qJD(3) + t615 * t643 + t720 * t654;
	t573 = -t688 * t636 - t739;
	t561 = t569 * t644 + t574 * t641 + (-t589 * t641 - t644 * t723) * qJD(6);
	t560 = -t569 * t641 + t574 * t644 + (-t589 * t644 + t641 * t723) * qJD(6);
	t1 = [(t573 * t644 + t735) * r_i_i_C(1) + (-t573 * t641 + t734) * r_i_i_C(2) + t573 * pkin(5) - t577 * t635 - t578 * t647 - t617 * pkin(2) - pkin(9) * t695 + t646 * t705 + t721 * t736 + (t738 * r_i_i_C(1) - t737 * r_i_i_C(2)) * qJD(6) + (-t646 * pkin(1) + (-pkin(9) * t682 - t719 * qJ(2)) * t640) * qJD(1) + (-t605 * t642 + (-t597 * t642 + t645 * t725) * qJD(4)) * pkin(4), t698, (t575 * t641 + t599 * t703) * r_i_i_C(1) + (t575 * t644 - t599 * t704) * r_i_i_C(2) - t575 * t647 + t659 * t574 - t648 * t723, -t575 * t642 * pkin(4) - t610 * t701 + t724 * t645 + t722, t722, r_i_i_C(1) * t560 - r_i_i_C(2) * t561; t610 * t645 * t717 - pkin(1) * t696 + t615 * pkin(2) + t569 * pkin(5) + pkin(9) * t650 + t561 * r_i_i_C(1) + t560 * r_i_i_C(2) + qJ(2) * t698 + t721 * t568 - t574 * t647 + t575 * t635 + t724 * t642 + t719 * t705, t640 * t696, (t577 * t641 - t597 * t703) * r_i_i_C(1) + (t577 * t644 + t597 * t704) * r_i_i_C(2) - t577 * t647 - t659 * t578 + t648 * t594, (-t577 * t642 + t605 * t645 + (t597 * t645 + t642 * t725) * qJD(4)) * pkin(4) + t652, t652, (-t571 * t641 - t734) * r_i_i_C(1) + (-t571 * t644 + t735) * r_i_i_C(2) + (t737 * r_i_i_C(1) + t738 * r_i_i_C(2)) * qJD(6); 0, 0, (-t602 * t641 + t607 * t703) * r_i_i_C(1) + (-t602 * t644 - t607 * t704) * r_i_i_C(2) + t602 * t647 + t659 * t603 + t648 * t606, (t602 * t642 + (-t607 * t645 - t618 * t642) * qJD(4)) * pkin(4) + t651, t651, (-t584 * t641 + t603 * t644) * r_i_i_C(1) + (-t584 * t644 - t603 * t641) * r_i_i_C(2) + ((-t593 * t644 - t606 * t641) * r_i_i_C(1) + (t593 * t641 - t606 * t644) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end