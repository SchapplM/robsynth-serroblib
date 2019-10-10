% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
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
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 08:58:25
	% EndTime: 2019-10-10 08:58:25
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(7));
	t278 = pkin(9) + r_i_i_C(3);
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
	% StartTime: 2019-10-10 08:58:26
	% EndTime: 2019-10-10 08:58:27
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
	t427 = r_i_i_C(1) * t407 + r_i_i_C(2) * t410;
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
	% StartTime: 2019-10-10 08:58:30
	% EndTime: 2019-10-10 08:58:32
	% DurationCPUTime: 1.62s
	% Computational Cost: add. (1177->127), mult. (3776->220), div. (0->0), fcn. (4194->14), ass. (0->86)
	t584 = cos(qJ(1));
	t641 = cos(pkin(12));
	t643 = cos(pkin(6));
	t617 = t643 * t641;
	t639 = sin(pkin(12));
	t645 = sin(qJ(1));
	t596 = t584 * t639 + t645 * t617;
	t559 = t596 * qJD(1);
	t615 = t643 * t639;
	t564 = t584 * t641 - t645 * t615;
	t560 = t564 * qJD(1);
	t581 = sin(qJ(3));
	t594 = -t584 * t615 - t645 * t641;
	t578 = sin(pkin(6));
	t640 = sin(pkin(7));
	t624 = t640 * t578;
	t607 = t645 * t624;
	t642 = cos(pkin(7));
	t646 = cos(qJ(3));
	t595 = t584 * t617 - t645 * t639;
	t621 = t584 * t624;
	t605 = t646 * t621;
	t620 = t642 * t646;
	t651 = t595 * t620 - t605;
	t526 = (qJD(1) * t607 + qJD(3) * t594 - t642 * t559) * t581 + t560 * t646 + t651 * qJD(3);
	t619 = t642 * t645;
	t608 = t578 * t619;
	t628 = t559 * t640;
	t548 = qJD(1) * t608 + t628;
	t580 = sin(qJ(4));
	t583 = cos(qJ(4));
	t625 = t642 * t595;
	t631 = t594 * t646;
	t541 = (t621 - t625) * t581 + t631;
	t635 = t578 * t584;
	t650 = t595 * t640 + t642 * t635;
	t613 = t541 * t580 - t583 * t650;
	t522 = -t613 * qJD(4) - t526 * t583 - t548 * t580;
	t601 = t646 * t607;
	t610 = t581 * t621;
	t527 = qJD(1) * t601 - t559 * t620 - t560 * t581 + (-t581 * t625 + t610 + t631) * qJD(3);
	t579 = sin(qJ(5));
	t582 = cos(qJ(5));
	t667 = t522 * t579 - t527 * t582;
	t666 = t522 * t582 + t527 * t579;
	t533 = t541 * t583 + t580 * t650;
	t538 = -t581 * t594 - t651;
	t665 = -t533 * t579 - t538 * t582;
	t664 = t533 * t582 - t538 * t579;
	t663 = t533 * qJD(4) - t526 * t580 + t548 * t583;
	t647 = r_i_i_C(3) + pkin(11);
	t590 = t596 * t642;
	t543 = t564 * t646 + (-t590 + t607) * t581;
	t553 = t596 * t640 + t608;
	t649 = -t543 * t580 + t553 * t583;
	t606 = qJD(5) * (t579 * r_i_i_C(1) + t582 * r_i_i_C(2));
	t648 = -t564 * t581 - t646 * t590 + t601;
	t614 = t642 * t641;
	t616 = t643 * t640;
	t549 = -t646 * t616 + (t581 * t639 - t614 * t646) * t578;
	t585 = t650 * qJD(1);
	t634 = qJD(5) * t579;
	t633 = qJD(5) * t582;
	t632 = t578 * qJD(2);
	t630 = qJD(1) * t635;
	t629 = qJD(1) * t645;
	t535 = t543 * t583 + t553 * t580;
	t550 = t581 * t616 + (t581 * t614 + t646 * t639) * t578;
	t561 = -t641 * t624 + t643 * t642;
	t537 = t550 * t583 + t561 * t580;
	t612 = -t550 * t580 + t561 * t583;
	t611 = t582 * r_i_i_C(1) - t579 * r_i_i_C(2) + pkin(4);
	t593 = -t647 * t580 - t611 * t583 - pkin(3);
	t588 = qJD(1) * t625;
	t586 = t583 * t606 + (t611 * t580 - t647 * t583) * qJD(4);
	t558 = t594 * qJD(1);
	t546 = t550 * qJD(3);
	t545 = t549 * qJD(3);
	t530 = t612 * qJD(4) - t545 * t583;
	t524 = qJD(1) * t610 + t648 * qJD(3) + t558 * t646 - t581 * t588;
	t523 = -qJD(1) * t605 + t543 * qJD(3) + t558 * t581 + t646 * t588;
	t518 = t649 * qJD(4) + t524 * t583 + t580 * t585;
	t517 = t535 * qJD(4) + t524 * t580 - t583 * t585;
	t516 = t518 * t582 + t523 * t579 + (-t535 * t579 - t582 * t648) * qJD(5);
	t515 = -t518 * t579 + t523 * t582 + (-t535 * t582 + t579 * t648) * qJD(5);
	t1 = [t666 * r_i_i_C(1) - t667 * r_i_i_C(2) + t522 * pkin(4) - t526 * pkin(3) + t527 * pkin(10) - t560 * pkin(2) - pkin(9) * t628 + t584 * t632 + t647 * t663 + (t665 * r_i_i_C(1) - t664 * r_i_i_C(2)) * qJD(5) + (-t584 * pkin(1) + (-pkin(9) * t619 - t645 * qJ(2)) * t578) * qJD(1), t630, (t524 * t579 + t543 * t633) * r_i_i_C(1) + (t524 * t582 - t543 * t634) * r_i_i_C(2) + t524 * pkin(10) + t593 * t523 - t586 * t648, -t611 * t517 + t647 * t518 - t649 * t606, r_i_i_C(1) * t515 - r_i_i_C(2) * t516, 0; -pkin(1) * t629 + t558 * pkin(2) + t524 * pkin(3) + t518 * pkin(4) + pkin(9) * t585 + t523 * pkin(10) + t516 * r_i_i_C(1) + t515 * r_i_i_C(2) + qJ(2) * t630 + t647 * t517 + t645 * t632, t578 * t629, (t526 * t579 - t541 * t633) * r_i_i_C(1) + (t526 * t582 + t541 * t634) * r_i_i_C(2) + t526 * pkin(10) - t593 * t527 + t586 * t538, -t522 * t647 - t613 * t606 + t611 * t663, t667 * r_i_i_C(1) + t666 * r_i_i_C(2) + (t664 * r_i_i_C(1) + t665 * r_i_i_C(2)) * qJD(5), 0; 0, 0, (-t545 * t579 + t550 * t633) * r_i_i_C(1) + (-t545 * t582 - t550 * t634) * r_i_i_C(2) - t545 * pkin(10) + t593 * t546 + t586 * t549, t647 * t530 - t612 * t606 + t611 * (-t537 * qJD(4) + t545 * t580), (-t530 * t579 + t546 * t582) * r_i_i_C(1) + (-t530 * t582 - t546 * t579) * r_i_i_C(2) + ((-t537 * t582 - t549 * t579) * r_i_i_C(1) + (t537 * t579 - t549 * t582) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:34
	% EndTime: 2019-10-10 08:58:36
	% DurationCPUTime: 2.21s
	% Computational Cost: add. (2090->149), mult. (6608->243), div. (0->0), fcn. (7481->14), ass. (0->98)
	t719 = cos(pkin(12));
	t721 = cos(pkin(6));
	t696 = t721 * t719;
	t717 = sin(pkin(12));
	t724 = sin(qJ(1));
	t726 = cos(qJ(1));
	t662 = t724 * t696 + t726 * t717;
	t631 = t662 * qJD(1);
	t694 = t721 * t717;
	t636 = -t724 * t694 + t726 * t719;
	t632 = t636 * qJD(1);
	t635 = t694 * t726 + t719 * t724;
	t653 = sin(qJ(3));
	t650 = sin(pkin(6));
	t718 = sin(pkin(7));
	t707 = t650 * t718;
	t681 = t724 * t707;
	t720 = cos(pkin(7));
	t725 = cos(qJ(3));
	t634 = -t696 * t726 + t717 * t724;
	t683 = t726 * t707;
	t674 = t725 * t683;
	t703 = t720 * t725;
	t730 = -t634 * t703 - t674;
	t597 = t653 * (qJD(1) * t681 - qJD(3) * t635 - t631 * t720) + t632 * t725 + t730 * qJD(3);
	t702 = t720 * t724;
	t682 = t650 * t702;
	t708 = t631 * t718;
	t620 = qJD(1) * t682 + t708;
	t652 = sin(qJ(4));
	t655 = cos(qJ(4));
	t706 = t720 * t634;
	t711 = t635 * t725;
	t612 = t653 * (t683 + t706) - t711;
	t710 = t650 * t726;
	t624 = t634 * t718 - t720 * t710;
	t687 = t612 * t652 + t624 * t655;
	t587 = -qJD(4) * t687 - t597 * t655 - t620 * t652;
	t673 = t725 * t681;
	t680 = t653 * t683;
	t598 = qJD(1) * t673 - t631 * t703 - t632 * t653 + (t653 * t706 + t680 - t711) * qJD(3);
	t651 = sin(qJ(5));
	t654 = cos(qJ(5));
	t604 = t612 * t655 - t624 * t652;
	t609 = t635 * t653 - t730;
	t743 = t604 * t654 - t609 * t651;
	t748 = t743 * qJD(5) + t587 * t651 - t598 * t654;
	t744 = t604 * t651 + t609 * t654;
	t747 = t744 * qJD(5) - t587 * t654 - t598 * t651;
	t742 = qJD(4) * t604 - t597 * t652 + t620 * t655;
	t728 = r_i_i_C(1) + pkin(5);
	t727 = r_i_i_C(2) + pkin(11);
	t723 = r_i_i_C(3) + qJ(6);
	t660 = t662 * t720;
	t729 = -t636 * t653 - t725 * t660 + t673;
	t693 = t720 * t719;
	t695 = t721 * t718;
	t621 = -t725 * t695 + (t653 * t717 - t693 * t725) * t650;
	t656 = t624 * qJD(1);
	t671 = t723 * t651 + t654 * t728 + pkin(4);
	t625 = t662 * t718 + t682;
	t716 = t625 * t655;
	t714 = t651 * t655;
	t713 = qJD(4) * t652;
	t712 = qJD(5) * t655;
	t709 = qJD(1) * t724;
	t705 = qJD(1) * t710;
	t630 = t635 * qJD(1);
	t658 = qJD(1) * t706;
	t595 = qJD(1) * t680 + t729 * qJD(3) - t630 * t725 + t653 * t658;
	t699 = -t712 * t729 + t595;
	t698 = t609 * t712 + t597;
	t617 = t621 * qJD(3);
	t697 = t621 * t712 - t617;
	t614 = t636 * t725 + (-t660 + t681) * t653;
	t606 = t614 * t655 + t625 * t652;
	t690 = t606 * t654 - t651 * t729;
	t689 = -t606 * t651 - t654 * t729;
	t622 = t653 * t695 + (t653 * t693 + t717 * t725) * t650;
	t633 = -t707 * t719 + t720 * t721;
	t608 = t622 * t655 + t633 * t652;
	t688 = t608 * t654 + t621 * t651;
	t686 = -t622 * t652 + t633 * t655;
	t679 = -t655 * pkin(4) - t727 * t652 - pkin(3);
	t675 = qJD(4) * (pkin(4) * t652 - t655 * t727);
	t594 = -qJD(1) * t674 + qJD(3) * t614 - t630 * t653 - t658 * t725;
	t669 = qJD(5) * t614 - t594 * t655 - t713 * t729;
	t668 = -qJD(5) * t612 + t598 * t655 + t609 * t713;
	t618 = t622 * qJD(3);
	t667 = qJD(5) * t622 - t618 * t655 + t621 * t713;
	t663 = qJD(6) * t651 + (-t728 * t651 + t723 * t654) * qJD(5);
	t601 = qJD(4) * t686 - t617 * t655;
	t588 = qJD(5) * t688 + t601 * t651 - t618 * t654;
	t583 = qJD(4) * t716 + t595 * t655 - t614 * t713 - t652 * t656;
	t582 = qJD(4) * t606 + t595 * t652 + t655 * t656;
	t573 = qJD(5) * t689 + t583 * t654 + t594 * t651;
	t572 = qJD(5) * t690 + t583 * t651 - t594 * t654;
	t1 = [t744 * qJD(6) + t587 * pkin(4) - t597 * pkin(3) + t598 * pkin(10) - t632 * pkin(2) - pkin(9) * t708 + qJD(2) * t710 + t727 * t742 - t728 * t747 + t723 * t748 + (-t726 * pkin(1) + (-pkin(9) * t702 - qJ(2) * t724) * t650) * qJD(1), t705, -(t614 * t654 - t714 * t729) * qJD(6) + t595 * pkin(10) + t728 * (t651 * t699 + t654 * t669) + t723 * (t651 * t669 - t654 * t699) - t729 * t675 + t679 * t594, t727 * t583 + t663 * (-t614 * t652 + t716) - t671 * t582, t690 * qJD(6) - t572 * t728 + t723 * t573, t572; qJD(2) * t650 * t724 - pkin(1) * t709 - t630 * pkin(2) + t595 * pkin(3) + t583 * pkin(4) - pkin(9) * t656 + t594 * pkin(10) + qJ(2) * t705 - qJD(6) * t689 + t723 * t572 + t728 * t573 + t727 * t582, t650 * t709, -(t609 * t714 - t612 * t654) * qJD(6) + t597 * pkin(10) + t728 * (t651 * t698 + t654 * t668) + t723 * (t651 * t668 - t654 * t698) + t609 * t675 - t679 * t598, -t587 * t727 + t663 * t687 + t671 * t742, -t743 * qJD(6) + t723 * t747 + t728 * t748, -t748; 0, 0, -(t621 * t714 + t622 * t654) * qJD(6) - t617 * pkin(10) + t728 * (t651 * t697 + t654 * t667) + t723 * (t651 * t667 - t654 * t697) + t621 * t675 + t679 * t618, t727 * t601 + t663 * t686 + t671 * (-qJD(4) * t608 + t617 * t652), t688 * qJD(6) + t723 * (t601 * t654 + t618 * t651 + (-t608 * t651 + t621 * t654) * qJD(5)) - t728 * t588, t588;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end