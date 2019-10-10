% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
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
	t115 = cos(pkin(12));
	t113 = sin(pkin(12));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:33
	% DurationCPUTime: 0.23s
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
	t244 = sin(pkin(12));
	t272 = t251 * t244;
	t247 = cos(pkin(12));
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
	% StartTime: 2019-10-10 01:37:34
	% EndTime: 2019-10-10 01:37:35
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (393->78), mult. (1301->142), div. (0->0), fcn. (1376->12), ass. (0->64)
	t408 = sin(qJ(3));
	t402 = sin(pkin(7));
	t403 = sin(pkin(6));
	t412 = cos(qJ(1));
	t442 = t403 * t412;
	t431 = t402 * t442;
	t405 = cos(pkin(7));
	t404 = cos(pkin(12));
	t406 = cos(pkin(6));
	t440 = t406 * t412;
	t401 = sin(pkin(12));
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
	% StartTime: 2019-10-10 01:37:34
	% EndTime: 2019-10-10 01:37:35
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (581->105), mult. (1716->177), div. (0->0), fcn. (1834->14), ass. (0->75)
	t424 = cos(pkin(12));
	t426 = cos(pkin(6));
	t421 = sin(pkin(12));
	t430 = sin(qJ(1));
	t462 = t430 * t421;
	t452 = t426 * t462;
	t433 = cos(qJ(1));
	t457 = qJD(1) * t433;
	t402 = -qJD(1) * t452 + t424 * t457;
	t422 = sin(pkin(7));
	t423 = sin(pkin(6));
	t465 = t423 * t433;
	t451 = t422 * t465;
	t488 = -qJD(3) * t451 + t402;
	t425 = cos(pkin(7));
	t459 = t433 * t424;
	t441 = t426 * t459 - t462;
	t482 = t441 * t425;
	t487 = -t482 + t451;
	t460 = t433 * t421;
	t461 = t430 * t424;
	t440 = t426 * t461 + t460;
	t401 = t440 * qJD(1);
	t405 = t426 * t460 + t461;
	t429 = sin(qJ(3));
	t432 = cos(qJ(3));
	t466 = t423 * t430;
	t449 = qJD(1) * t466;
	t380 = (-qJD(3) * t405 - t401 * t425 + t422 * t449) * t429 + (qJD(3) * t482 + t488) * t432;
	t472 = t401 * t422;
	t392 = t425 * t449 + t472;
	t420 = qJ(4) + pkin(13);
	t418 = sin(t420);
	t419 = cos(t420);
	t486 = t380 * t418 - t392 * t419;
	t485 = -t380 * t419 - t392 * t418;
	t481 = t441 * t422 + t425 * t465;
	t463 = t425 * t432;
	t468 = t422 * t426;
	t480 = (-t421 * t429 + t424 * t463) * t423 + t432 * t468;
	t464 = t425 * t429;
	t444 = -t405 * t432 - t441 * t464;
	t385 = t429 * t451 + t444;
	t479 = t385 * t419 + t418 * t481;
	t478 = -t385 * t418 + t419 * t481;
	t476 = t405 * t429 + t487 * t432;
	t458 = qJD(1) * t432;
	t467 = t423 * t422;
	t448 = t458 * t467;
	t475 = t444 * qJD(3) - t401 * t463 - t488 * t429 + t430 * t448;
	t428 = sin(qJ(4));
	t474 = t428 * pkin(4);
	t473 = r_i_i_C(3) + qJ(5) + pkin(10);
	t455 = t423 * qJD(2);
	t431 = cos(qJ(4));
	t417 = t431 * pkin(4) + pkin(3);
	t445 = -t419 * r_i_i_C(1) + t418 * r_i_i_C(2) - t417;
	t442 = t422 * t466 - t425 * t440;
	t439 = t418 * r_i_i_C(1) + t419 * r_i_i_C(2) + t474;
	t436 = qJD(4) * t439;
	t407 = -t452 + t459;
	t386 = -t407 * t429 + t442 * t432;
	t387 = t407 * t432 + t442 * t429;
	t394 = t429 * t468 + (t421 * t432 + t424 * t464) * t423;
	t403 = -t424 * t467 + t426 * t425;
	t400 = t405 * qJD(1);
	t397 = t422 * t440 + t425 * t466;
	t390 = t481 * qJD(1);
	t389 = t394 * qJD(3);
	t388 = t480 * qJD(3);
	t378 = t487 * t429 * qJD(1) + t386 * qJD(3) - t400 * t432;
	t377 = t387 * qJD(3) - t400 * t429 - t433 * t448 + t458 * t482;
	t376 = t378 * t419 + t390 * t418 + (-t387 * t418 + t397 * t419) * qJD(4);
	t375 = -t378 * t418 + t390 * t419 + (-t387 * t419 - t397 * t418) * qJD(4);
	t1 = [t485 * r_i_i_C(1) + t486 * r_i_i_C(2) - t380 * t417 - t476 * qJD(5) - t392 * t474 - t402 * pkin(2) - pkin(9) * t472 + t433 * t455 + t473 * t475 + (-t433 * pkin(1) + (-pkin(9) * t425 - qJ(2)) * t466) * qJD(1) + (t478 * r_i_i_C(1) - t479 * r_i_i_C(2) + (-t385 * t428 + t431 * t481) * pkin(4)) * qJD(4), t423 * t457, t387 * qJD(5) + t445 * t377 + t473 * t378 - t386 * t436, t375 * r_i_i_C(1) - t376 * r_i_i_C(2) + (-t378 * t428 + t390 * t431 + (-t387 * t431 - t397 * t428) * qJD(4)) * pkin(4), t377, 0; t430 * t455 - t400 * pkin(2) + t376 * r_i_i_C(1) + t375 * r_i_i_C(2) - t386 * qJD(5) + t378 * t417 + t473 * t377 + (t390 * t428 + (-t387 * t428 + t397 * t431) * qJD(4)) * pkin(4) + (-t430 * pkin(1) + pkin(9) * t481 + qJ(2) * t465) * qJD(1), t449, -qJD(5) * t385 + t473 * t380 + t476 * t436 - t445 * t475, -t486 * r_i_i_C(1) + t485 * r_i_i_C(2) + (t479 * r_i_i_C(1) + t478 * r_i_i_C(2)) * qJD(4) + (-t380 * t428 + t392 * t431 + (t385 * t431 + t428 * t481) * qJD(4)) * pkin(4), -t475, 0; 0, 0, t394 * qJD(5) + t473 * t388 + t445 * t389 - t480 * t436, -t439 * t388 + ((-t394 * t419 - t403 * t418) * r_i_i_C(1) + (t394 * t418 - t403 * t419) * r_i_i_C(2) + (-t394 * t431 - t403 * t428) * pkin(4)) * qJD(4), t389, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:38
	% EndTime: 2019-10-10 01:37:40
	% DurationCPUTime: 1.94s
	% Computational Cost: add. (1500->152), mult. (4191->252), div. (0->0), fcn. (4652->16), ass. (0->94)
	t601 = cos(qJ(1));
	t657 = cos(pkin(12));
	t659 = cos(pkin(6));
	t633 = t659 * t657;
	t655 = sin(pkin(12));
	t664 = sin(qJ(1));
	t614 = t601 * t655 + t664 * t633;
	t570 = t614 * qJD(1);
	t632 = t659 * t655;
	t575 = t601 * t657 - t664 * t632;
	t571 = t575 * qJD(1);
	t598 = sin(qJ(3));
	t612 = -t601 * t632 - t664 * t657;
	t594 = sin(pkin(6));
	t656 = sin(pkin(7));
	t642 = t594 * t656;
	t624 = t664 * t642;
	t658 = cos(pkin(7));
	t665 = cos(qJ(3));
	t613 = t601 * t633 - t664 * t655;
	t639 = t601 * t642;
	t622 = t665 * t639;
	t636 = t658 * t665;
	t671 = t613 * t636 - t622;
	t537 = (qJD(1) * t624 + qJD(3) * t612 - t658 * t570) * t598 + t571 * t665 + t671 * qJD(3);
	t635 = t658 * t664;
	t625 = t594 * t635;
	t643 = t570 * t656;
	t559 = qJD(1) * t625 + t643;
	t593 = qJ(4) + pkin(13);
	t591 = sin(t593);
	t592 = cos(t593);
	t653 = t594 * t601;
	t670 = t613 * t656 + t658 * t653;
	t583 = t598 * t639;
	t672 = t613 * t658;
	t680 = -t598 * t672 + t612 * t665 + t583;
	t669 = -t591 * t680 + t592 * t670;
	t533 = t669 * qJD(4) - t537 * t592 - t559 * t591;
	t618 = t665 * t624;
	t538 = qJD(1) * t618 + t680 * qJD(3) - t570 * t636 - t571 * t598;
	t596 = sin(qJ(6));
	t599 = cos(qJ(6));
	t687 = t533 * t596 - t538 * t599;
	t686 = t533 * t599 + t538 * t596;
	t543 = -t591 * t670 - t592 * t680;
	t549 = -t598 * t612 - t671;
	t685 = t543 * t596 - t549 * t599;
	t684 = t543 * t599 + t549 * t596;
	t683 = -t543 * qJD(4) - t537 * t591 + t559 * t592;
	t608 = t614 * t658;
	t554 = t575 * t665 + (-t608 + t624) * t598;
	t604 = t670 * qJD(1);
	t679 = -t554 * qJD(4) + t604;
	t666 = r_i_i_C(3) + pkin(11);
	t668 = t679 * pkin(4);
	t623 = qJD(6) * (t596 * r_i_i_C(1) + t599 * r_i_i_C(2));
	t627 = t599 * r_i_i_C(1) - t596 * r_i_i_C(2) + pkin(5);
	t667 = -t575 * t598 - t665 * t608 + t618;
	t630 = t656 * t659;
	t631 = t658 * t657;
	t560 = -t665 * t630 + (t598 * t655 - t631 * t665) * t594;
	t597 = sin(qJ(4));
	t662 = t597 * pkin(4);
	t652 = qJD(6) * t596;
	t651 = qJD(6) * t599;
	t564 = t614 * t656 + t625;
	t649 = t564 * qJD(4);
	t648 = t594 * qJD(2);
	t647 = pkin(4) * t649;
	t645 = qJD(1) * t653;
	t644 = qJD(1) * t664;
	t546 = t554 * t592 + t564 * t591;
	t561 = t598 * t630 + (t598 * t631 + t665 * t655) * t594;
	t572 = -t657 * t642 + t659 * t658;
	t548 = t561 * t592 + t572 * t591;
	t628 = -t561 * t591 + t572 * t592;
	t600 = cos(qJ(4));
	t590 = t600 * pkin(4) + pkin(3);
	t611 = -t666 * t591 - t627 * t592 - t590;
	t606 = qJD(1) * t672;
	t602 = t592 * t623 + (t627 * t591 - t666 * t592 + t662) * qJD(4);
	t595 = -qJ(5) - pkin(10);
	t569 = t612 * qJD(1);
	t557 = t561 * qJD(3);
	t556 = t560 * qJD(3);
	t541 = t628 * qJD(4) - t556 * t592;
	t535 = qJD(1) * t583 + t667 * qJD(3) + t569 * t665 - t598 * t606;
	t534 = -qJD(1) * t622 + t554 * qJD(3) + t569 * t598 + t665 * t606;
	t529 = (t535 + t649) * t592 + t679 * t591;
	t528 = t546 * qJD(4) + t535 * t591 - t592 * t604;
	t527 = t529 * t599 + t534 * t596 + (-t546 * t596 - t599 * t667) * qJD(6);
	t526 = -t529 * t596 + t534 * t599 + (-t546 * t599 + t596 * t667) * qJD(6);
	t1 = [t686 * r_i_i_C(1) - t687 * r_i_i_C(2) + t533 * pkin(5) - t537 * t590 - t538 * t595 - t549 * qJD(5) - t571 * pkin(2) - pkin(9) * t643 + t601 * t648 + t666 * t683 + (t685 * r_i_i_C(1) + t684 * r_i_i_C(2)) * qJD(6) + (-t601 * pkin(1) + (-pkin(9) * t635 - t664 * qJ(2)) * t594) * qJD(1) + (-t559 * t597 + (-t597 * t680 + t600 * t670) * qJD(4)) * pkin(4), t645, (t535 * t596 + t554 * t651) * r_i_i_C(1) + (t535 * t599 - t554 * t652) * r_i_i_C(2) - t535 * t595 + t554 * qJD(5) + t611 * t534 - t602 * t667, -t535 * t662 - t597 * t647 + t668 * t600 + (-r_i_i_C(1) * t652 - r_i_i_C(2) * t651) * (-t554 * t591 + t564 * t592) + t666 * t529 - t627 * t528, t534, t526 * r_i_i_C(1) - t527 * r_i_i_C(2); -pkin(1) * t644 + t569 * pkin(2) + t529 * pkin(5) + pkin(9) * t604 + t527 * r_i_i_C(1) + t526 * r_i_i_C(2) + qJ(2) * t645 - qJD(5) * t667 + t666 * t528 - t534 * t595 + t535 * t590 + t668 * t597 + t600 * t647 + t664 * t648, t594 * t644, (t537 * t596 - t651 * t680) * r_i_i_C(1) + (t537 * t599 + t652 * t680) * r_i_i_C(2) - t537 * t595 - t680 * qJD(5) - t611 * t538 + t602 * t549, -t666 * t533 + t669 * t623 + t627 * t683 + (-t537 * t597 + t559 * t600 + (t597 * t670 + t600 * t680) * qJD(4)) * pkin(4), -t538, t687 * r_i_i_C(1) + t686 * r_i_i_C(2) + (-t684 * r_i_i_C(1) + t685 * r_i_i_C(2)) * qJD(6); 0, 0, (-t556 * t596 + t561 * t651) * r_i_i_C(1) + (-t556 * t599 - t561 * t652) * r_i_i_C(2) + t556 * t595 + t561 * qJD(5) + t611 * t557 + t602 * t560, t666 * t541 - t628 * t623 + t627 * (-t548 * qJD(4) + t556 * t591) + (t556 * t597 + (-t561 * t600 - t572 * t597) * qJD(4)) * pkin(4), t557, (-t541 * t596 + t557 * t599) * r_i_i_C(1) + (-t541 * t599 - t557 * t596) * r_i_i_C(2) + ((-t548 * t599 - t560 * t596) * r_i_i_C(1) + (t548 * t596 - t560 * t599) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end