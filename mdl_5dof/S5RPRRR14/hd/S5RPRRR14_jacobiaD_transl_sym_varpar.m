% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR14
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:38
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(5));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(5));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(11));
	t113 = sin(pkin(11));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:39
	% EndTime: 2019-12-29 18:05:39
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(6));
	t278 = r_i_i_C(3) + pkin(8);
	t279 = t278 * t248 + qJ(2);
	t250 = sin(qJ(3));
	t277 = t250 * r_i_i_C(2);
	t246 = sin(pkin(5));
	t251 = sin(qJ(1));
	t276 = t246 * t251;
	t253 = cos(qJ(1));
	t275 = t246 * t253;
	t274 = t248 * t250;
	t252 = cos(qJ(3));
	t273 = t248 * t252;
	t244 = sin(pkin(11));
	t272 = t251 * t244;
	t247 = cos(pkin(11));
	t271 = t251 * t247;
	t270 = t253 * t244;
	t269 = t253 * t247;
	t268 = qJD(1) * t251;
	t267 = qJD(1) * t253;
	t245 = sin(pkin(6));
	t266 = qJD(3) * t245;
	t265 = t278 * t245;
	t249 = cos(pkin(5));
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
	t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0; 0, 0, (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:41
	% EndTime: 2019-12-29 18:05:42
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (393->78), mult. (1301->142), div. (0->0), fcn. (1376->12), ass. (0->64)
	t408 = sin(qJ(3));
	t402 = sin(pkin(6));
	t403 = sin(pkin(5));
	t412 = cos(qJ(1));
	t442 = t403 * t412;
	t431 = t402 * t442;
	t405 = cos(pkin(6));
	t404 = cos(pkin(11));
	t406 = cos(pkin(5));
	t440 = t406 * t412;
	t401 = sin(pkin(11));
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
	t450 = r_i_i_C(3) + pkin(9);
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
	t1 = [t461 * r_i_i_C(1) + t462 * r_i_i_C(2) - t368 * pkin(3) - t389 * pkin(2) - pkin(8) * t449 + t412 * t436 + t450 * t413 + (t460 * r_i_i_C(1) - t459 * r_i_i_C(2)) * qJD(4) + (-t412 * pkin(1) + (-pkin(8) * t405 - qJ(2)) * t443) * qJD(1), t403 * t438, -t425 * t365 + t450 * t366 - t415 * t419, r_i_i_C(1) * t363 - r_i_i_C(2) * t364, 0; t409 * t436 - t387 * pkin(2) + t366 * pkin(3) + t364 * r_i_i_C(1) + t363 * r_i_i_C(2) + t450 * t365 + (-t409 * pkin(1) + pkin(8) * t453 + qJ(2) * t442) * qJD(1), t429, t450 * t368 - (-t392 * t408 - t424 * t411) * t419 + t425 * t413, -t462 * r_i_i_C(1) + t461 * r_i_i_C(2) + (t459 * r_i_i_C(1) + t460 * r_i_i_C(2)) * qJD(4), 0; 0, 0, -t425 * t382 * qJD(3) + t450 * t376 - t452 * t419, -t427 * t376 + ((-t382 * t410 - t390 * t407) * r_i_i_C(1) + (t382 * t407 - t390 * t410) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:47
	% EndTime: 2019-12-29 18:05:49
	% DurationCPUTime: 2.40s
	% Computational Cost: add. (1177->127), mult. (3776->220), div. (0->0), fcn. (4194->14), ass. (0->86)
	t584 = cos(qJ(1));
	t641 = cos(pkin(11));
	t643 = cos(pkin(5));
	t617 = t643 * t641;
	t639 = sin(pkin(11));
	t645 = sin(qJ(1));
	t596 = t584 * t639 + t645 * t617;
	t559 = t596 * qJD(1);
	t615 = t643 * t639;
	t564 = t584 * t641 - t645 * t615;
	t560 = t564 * qJD(1);
	t581 = sin(qJ(3));
	t594 = -t584 * t615 - t645 * t641;
	t578 = sin(pkin(5));
	t640 = sin(pkin(6));
	t627 = t578 * t640;
	t607 = t645 * t627;
	t642 = cos(pkin(6));
	t646 = cos(qJ(3));
	t595 = t584 * t617 - t645 * t639;
	t621 = t584 * t627;
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
	t624 = t642 * t595;
	t631 = t594 * t646;
	t541 = (t621 - t624) * t581 + t631;
	t635 = t578 * t584;
	t650 = t595 * t640 + t642 * t635;
	t613 = t541 * t580 - t583 * t650;
	t522 = -t613 * qJD(4) - t526 * t583 - t548 * t580;
	t601 = t646 * t607;
	t610 = t581 * t621;
	t527 = qJD(1) * t601 - t559 * t620 - t560 * t581 + (-t581 * t624 + t610 + t631) * qJD(3);
	t579 = sin(qJ(5));
	t582 = cos(qJ(5));
	t667 = t522 * t579 - t527 * t582;
	t666 = t522 * t582 + t527 * t579;
	t533 = t541 * t583 + t580 * t650;
	t538 = -t581 * t594 - t651;
	t665 = -t533 * t579 - t538 * t582;
	t664 = t533 * t582 - t538 * t579;
	t663 = t533 * qJD(4) - t526 * t580 + t548 * t583;
	t647 = r_i_i_C(3) + pkin(10);
	t590 = t596 * t642;
	t543 = t564 * t646 + (-t590 + t607) * t581;
	t553 = t596 * t640 + t608;
	t649 = -t543 * t580 + t553 * t583;
	t606 = qJD(5) * (t579 * r_i_i_C(1) + t582 * r_i_i_C(2));
	t648 = -t564 * t581 - t646 * t590 + t601;
	t614 = t642 * t641;
	t616 = t643 * t640;
	t549 = -t646 * t616 + (t581 * t639 - t614 * t646) * t578;
	t585 = qJD(1) * t650;
	t634 = qJD(2) * t578;
	t633 = qJD(5) * t579;
	t632 = qJD(5) * t582;
	t630 = qJD(1) * t635;
	t629 = qJD(1) * t645;
	t535 = t543 * t583 + t553 * t580;
	t550 = t581 * t616 + (t581 * t614 + t646 * t639) * t578;
	t561 = -t641 * t627 + t643 * t642;
	t537 = t550 * t583 + t561 * t580;
	t612 = -t550 * t580 + t561 * t583;
	t611 = t582 * r_i_i_C(1) - t579 * r_i_i_C(2) + pkin(4);
	t593 = -t647 * t580 - t611 * t583 - pkin(3);
	t588 = qJD(1) * t624;
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
	t1 = [t666 * r_i_i_C(1) - t667 * r_i_i_C(2) + t522 * pkin(4) - t526 * pkin(3) + t527 * pkin(9) - t560 * pkin(2) - pkin(8) * t628 + t584 * t634 + t647 * t663 + (t665 * r_i_i_C(1) - t664 * r_i_i_C(2)) * qJD(5) + (-t584 * pkin(1) + (-pkin(8) * t619 - t645 * qJ(2)) * t578) * qJD(1), t630, (t524 * t579 + t543 * t632) * r_i_i_C(1) + (t524 * t582 - t543 * t633) * r_i_i_C(2) + t524 * pkin(9) + t593 * t523 - t586 * t648, -t611 * t517 + t647 * t518 - t649 * t606, r_i_i_C(1) * t515 - r_i_i_C(2) * t516; -pkin(1) * t629 + t558 * pkin(2) + t524 * pkin(3) + t518 * pkin(4) + pkin(8) * t585 + t523 * pkin(9) + t516 * r_i_i_C(1) + t515 * r_i_i_C(2) + qJ(2) * t630 + t647 * t517 + t645 * t634, t578 * t629, (t526 * t579 - t541 * t632) * r_i_i_C(1) + (t526 * t582 + t541 * t633) * r_i_i_C(2) + t526 * pkin(9) - t593 * t527 + t586 * t538, -t522 * t647 - t613 * t606 + t611 * t663, t667 * r_i_i_C(1) + t666 * r_i_i_C(2) + (t664 * r_i_i_C(1) + t665 * r_i_i_C(2)) * qJD(5); 0, 0, (-t545 * t579 + t550 * t632) * r_i_i_C(1) + (-t545 * t582 - t550 * t633) * r_i_i_C(2) - t545 * pkin(9) + t593 * t546 + t586 * t549, t647 * t530 - t612 * t606 + t611 * (-t537 * qJD(4) + t545 * t580), (-t530 * t579 + t546 * t582) * r_i_i_C(1) + (-t530 * t582 - t546 * t579) * r_i_i_C(2) + ((-t537 * t582 - t549 * t579) * r_i_i_C(1) + (t537 * t579 - t549 * t582) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end