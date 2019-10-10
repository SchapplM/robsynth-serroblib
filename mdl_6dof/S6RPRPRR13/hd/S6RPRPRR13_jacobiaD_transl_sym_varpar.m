% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR13_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:45
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
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
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
	% StartTime: 2019-10-10 01:07:46
	% EndTime: 2019-10-10 01:07:46
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (232->51), mult. (778->91), div. (0->0), fcn. (800->10), ass. (0->51)
	t324 = cos(pkin(12));
	t326 = cos(pkin(6));
	t321 = sin(pkin(12));
	t328 = sin(qJ(1));
	t349 = t328 * t321;
	t341 = t326 * t349;
	t330 = cos(qJ(1));
	t345 = qJD(1) * t330;
	t313 = -qJD(1) * t341 + t324 * t345;
	t327 = sin(qJ(3));
	t329 = cos(qJ(3));
	t347 = t330 * t321;
	t348 = t328 * t324;
	t315 = t326 * t347 + t348;
	t346 = t330 * t324;
	t314 = -t326 * t346 + t349;
	t322 = sin(pkin(7));
	t325 = cos(pkin(7));
	t323 = sin(pkin(6));
	t351 = t323 * t330;
	t337 = t314 * t325 + t322 * t351;
	t331 = -t315 * t329 + t337 * t327;
	t334 = t326 * t348 + t347;
	t312 = t334 * qJD(1);
	t352 = t323 * t328;
	t340 = qJD(1) * t352;
	t333 = -t312 * t325 + t322 * t340;
	t302 = -t331 * qJD(3) + t313 * t327 - t333 * t329;
	t332 = t315 * t327 + t337 * t329;
	t370 = -t332 * qJD(3) + t313 * t329 + t333 * t327;
	t362 = r_i_i_C(1) + pkin(9);
	t367 = t362 * t325 + qJ(2);
	t353 = t322 * t326;
	t366 = (t324 * t325 * t327 + t321 * t329) * t323 + t327 * t353;
	t317 = -t341 + t346;
	t336 = t322 * t352 - t325 * t334;
	t364 = t317 * t329 + t336 * t327;
	t361 = r_i_i_C(2) - pkin(3);
	t360 = r_i_i_C(3) + qJ(4);
	t355 = t317 * t327;
	t350 = t325 * t329;
	t344 = t323 * qJD(2);
	t343 = t362 * t322;
	t339 = t323 * t345;
	t338 = t322 * t339;
	t311 = t315 * qJD(1);
	t310 = t314 * qJD(1);
	t307 = t366 * qJD(3);
	t301 = -qJD(3) * t355 + (t310 * t325 + t338) * t327 + (t336 * qJD(3) - t311) * t329;
	t300 = t364 * qJD(3) - t310 * t350 - t311 * t327 - t329 * t338;
	t1 = [-t332 * qJD(4) - t313 * pkin(2) + t330 * t344 - t312 * t343 + t361 * t370 - t360 * t302 + (-t330 * pkin(1) - t367 * t352) * qJD(1), t339, t364 * qJD(4) + t361 * t300 + t360 * t301, t300, 0, 0; -(t336 * t329 - t355) * qJD(4) - t311 * pkin(2) + t328 * t344 - t310 * t343 - t361 * t301 + t360 * t300 + (-t328 * pkin(1) + t367 * t351) * qJD(1), t340, -t331 * qJD(4) + t361 * t302 + t360 * t370, t302, 0, 0; 0, 0, t366 * qJD(4) + t361 * t307 - t360 * (-t329 * t353 + (t321 * t327 - t324 * t350) * t323) * qJD(3), t307, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:47
	% EndTime: 2019-10-10 01:07:47
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (485->92), mult. (1598->158), div. (0->0), fcn. (1696->12), ass. (0->74)
	t408 = cos(pkin(6));
	t403 = sin(pkin(12));
	t414 = cos(qJ(1));
	t441 = t414 * t403;
	t406 = cos(pkin(12));
	t411 = sin(qJ(1));
	t442 = t411 * t406;
	t420 = t408 * t442 + t441;
	t388 = t420 * qJD(1);
	t443 = t411 * t403;
	t432 = t408 * t443;
	t439 = qJD(1) * t414;
	t389 = -qJD(1) * t432 + t406 * t439;
	t440 = t414 * t406;
	t391 = -t408 * t440 + t443;
	t404 = sin(pkin(7));
	t410 = sin(qJ(3));
	t413 = cos(qJ(3));
	t405 = sin(pkin(6));
	t447 = t405 * t411;
	t431 = qJD(1) * t447;
	t427 = t404 * t431;
	t438 = qJD(3) * t410;
	t429 = t405 * t438;
	t407 = cos(pkin(7));
	t444 = t407 * t413;
	t445 = t407 * t410;
	t392 = t408 * t441 + t442;
	t449 = t392 * t413;
	t366 = t414 * t404 * t429 + (t391 * t445 - t449) * qJD(3) - t388 * t444 - t389 * t410 + t413 * t427;
	t452 = t388 * t404;
	t377 = t407 * t431 + t452;
	t409 = sin(qJ(5));
	t412 = cos(qJ(5));
	t462 = t366 * t409 - t377 * t412;
	t461 = t366 * t412 + t377 * t409;
	t446 = t405 * t414;
	t435 = t404 * t446;
	t422 = t391 * t407 + t435;
	t450 = t392 * t410;
	t458 = (t422 * t413 + t450) * qJD(3) - (-t388 * t407 + t427) * t410 - t389 * t413;
	t368 = t391 * t444 + t413 * t435 + t450;
	t380 = -t391 * t404 + t407 * t446;
	t457 = -t368 * t412 - t380 * t409;
	t456 = t368 * t409 - t380 * t412;
	t386 = t391 * qJD(1);
	t453 = t386 * t404;
	t448 = t404 * t408;
	t437 = t405 * qJD(2);
	t436 = r_i_i_C(3) + pkin(10) + pkin(3);
	t434 = t406 * t444;
	t433 = t413 * t448;
	t430 = t405 * t439;
	t428 = pkin(9) * t407 + qJ(2);
	t426 = t404 * t430;
	t425 = t412 * r_i_i_C(1) - t409 * r_i_i_C(2);
	t424 = t409 * r_i_i_C(1) + t412 * r_i_i_C(2) + qJ(4);
	t421 = t404 * t447 - t407 * t420;
	t418 = t425 * qJD(5) + qJD(4);
	t394 = -t432 + t440;
	t416 = t394 * t413 + t421 * t410;
	t415 = t410 * t448 + (t403 * t413 + t406 * t445) * t405;
	t390 = -t405 * t406 * t404 + t408 * t407;
	t387 = t392 * qJD(1);
	t382 = t404 * t420 + t407 * t447;
	t378 = -t433 + (t403 * t410 - t434) * t405;
	t375 = t407 * t430 - t453;
	t374 = t415 * qJD(3);
	t371 = t394 * t410 - t421 * t413;
	t363 = -t394 * t438 + (t386 * t407 + t426) * t410 + (t421 * qJD(3) - t387) * t413;
	t362 = t416 * qJD(3) - t386 * t444 - t387 * t410 - t413 * t426;
	t361 = t362 * t409 + t375 * t412 + (t371 * t412 - t382 * t409) * qJD(5);
	t360 = t362 * t412 - t375 * t409 + (-t371 * t409 - t382 * t412) * qJD(5);
	t1 = [t462 * r_i_i_C(1) + t461 * r_i_i_C(2) - t377 * pkin(4) + t366 * qJ(4) - t368 * qJD(4) - t389 * pkin(2) - pkin(9) * t452 + t414 * t437 + (t457 * r_i_i_C(1) + t456 * r_i_i_C(2)) * qJD(5) + t436 * t458 + (-t414 * pkin(1) - t428 * t447) * qJD(1), t430, -t436 * t362 + t424 * t363 + t418 * t416, t362, t360 * r_i_i_C(1) - t361 * r_i_i_C(2), 0; -pkin(9) * t453 + t411 * t437 - t387 * pkin(2) + t375 * pkin(4) + t361 * r_i_i_C(1) + t360 * r_i_i_C(2) + t362 * qJ(4) + t371 * qJD(4) + t436 * t363 + (-pkin(1) * t411 + t428 * t446) * qJD(1), t431, t418 * (-t422 * t410 + t449) - t424 * t458 + t436 * t366, -t366, -t461 * r_i_i_C(1) + t462 * r_i_i_C(2) + (-t456 * r_i_i_C(1) + t457 * r_i_i_C(2)) * qJD(5), 0; 0, 0, t418 * t415 - t436 * t374 - t424 * (t403 * t429 + (-t405 * t434 - t433) * qJD(3)), t374, t425 * t374 + ((-t378 * t409 - t390 * t412) * r_i_i_C(1) + (-t378 * t412 + t390 * t409) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:51
	% EndTime: 2019-10-10 01:07:52
	% DurationCPUTime: 1.31s
	% Computational Cost: add. (1269->129), mult. (4073->213), div. (0->0), fcn. (4514->14), ass. (0->91)
	t581 = sin(qJ(3));
	t585 = cos(qJ(3));
	t574 = sin(pkin(12));
	t577 = cos(pkin(12));
	t582 = sin(qJ(1));
	t586 = cos(qJ(1));
	t636 = cos(pkin(6));
	t612 = t586 * t636;
	t592 = -t582 * t574 + t577 * t612;
	t593 = -t574 * t612 - t582 * t577;
	t575 = sin(pkin(7));
	t576 = sin(pkin(6));
	t627 = t576 * t586;
	t618 = t575 * t627;
	t578 = cos(pkin(7));
	t626 = t578 * t581;
	t535 = t581 * t618 + t585 * t593 - t592 * t626;
	t613 = t582 * t636;
	t594 = t586 * t574 + t577 * t613;
	t556 = t594 * qJD(1);
	t607 = t574 * t613;
	t623 = qJD(1) * t586;
	t557 = -qJD(1) * t607 + t577 * t623;
	t624 = qJD(1) * t582;
	t616 = t576 * t624;
	t609 = t575 * t616;
	t625 = t578 * t585;
	t521 = t535 * qJD(3) - t556 * t625 - t557 * t581 + t585 * t609;
	t634 = t556 * t575;
	t544 = t578 * t616 + t634;
	t580 = sin(qJ(5));
	t584 = cos(qJ(5));
	t654 = t521 * t580 - t544 * t584;
	t653 = -t521 * t584 - t544 * t580;
	t631 = t593 * t581;
	t639 = (-t578 * t592 + t618) * t585 - t631;
	t652 = t639 * qJD(3) - (-t556 * t578 + t609) * t581 - t557 * t585;
	t579 = sin(qJ(6));
	t651 = t535 * t579;
	t583 = cos(qJ(6));
	t650 = t535 * t583;
	t547 = t592 * t575 + t578 * t627;
	t647 = t547 * t580;
	t646 = t547 * t584;
	t645 = pkin(10) + pkin(3);
	t638 = r_i_i_C(3) + pkin(11);
	t644 = pkin(9) * t578 + qJ(2);
	t590 = qJD(1) * t592;
	t615 = t576 * t623;
	t643 = t575 * t615 - t578 * t590;
	t614 = t575 * t636;
	t629 = t576 * t577;
	t545 = t576 * t574 * t581 - t585 * t614 - t625 * t629;
	t558 = -t575 * t629 + t636 * t578;
	t642 = -t545 * t584 + t558 * t580;
	t562 = t586 * t577 - t607;
	t628 = t576 * t582;
	t620 = t575 * t628;
	t641 = -t562 * t581 + t585 * t620 - t594 * t625;
	t605 = -t583 * r_i_i_C(1) + t579 * r_i_i_C(2);
	t601 = pkin(5) - t605;
	t588 = t601 * t580 - t638 * t584 + qJ(4);
	t622 = qJD(2) * t576;
	t604 = -t579 * r_i_i_C(1) - t583 * r_i_i_C(2);
	t532 = t585 * t618 - t592 * t625 - t631;
	t603 = t532 * t584 + t647;
	t527 = t532 * t580 - t646;
	t525 = -t580 * t639 + t646;
	t549 = t575 * t594 + t578 * t628;
	t602 = -t549 * t580 - t584 * t641;
	t529 = t549 * t584 - t580 * t641;
	t531 = t545 * t580 + t558 * t584;
	t600 = -t604 + t645;
	t597 = qJD(6) * t605;
	t596 = qJD(6) * t604;
	t537 = t562 * t585 + (-t578 * t594 + t620) * t581;
	t546 = t581 * t614 + (t574 * t585 + t577 * t626) * t576;
	t587 = qJD(4) + t580 * t596 + (t638 * t580 + t601 * t584) * qJD(5);
	t555 = t593 * qJD(1);
	t542 = t547 * qJD(1);
	t541 = t546 * qJD(3);
	t540 = t545 * qJD(3);
	t523 = t642 * qJD(5) - t541 * t580;
	t518 = t641 * qJD(3) + t555 * t585 + t643 * t581;
	t517 = t537 * qJD(3) + t555 * t581 - t643 * t585;
	t515 = t603 * qJD(5) - t654;
	t511 = t602 * qJD(5) + t517 * t580 + t542 * t584;
	t510 = t529 * qJD(5) - t517 * t584 + t542 * t580;
	t509 = t511 * t583 + t518 * t579 + (-t529 * t579 + t537 * t583) * qJD(6);
	t508 = -t511 * t579 + t518 * t583 + (-t529 * t583 - t537 * t579) * qJD(6);
	t1 = [-pkin(9) * t634 + t586 * t622 - t557 * pkin(2) - t544 * pkin(4) + t521 * qJ(4) - t639 * qJD(4) + t601 * ((-t584 * t639 - t647) * qJD(5) + t654) + t600 * t652 + t638 * (t525 * qJD(5) + t653) + ((-t525 * t579 + t650) * r_i_i_C(1) + (-t525 * t583 - t651) * r_i_i_C(2)) * qJD(6) + (-t586 * pkin(1) - t644 * t628) * qJD(1), t615, -t600 * t517 + t588 * t518 + t587 * t537 - t597 * t641, t517, -t601 * t510 + t638 * t511 + t602 * t596, t508 * r_i_i_C(1) - t509 * r_i_i_C(2); t575 * pkin(9) * t590 - pkin(1) * t624 + t555 * pkin(2) + t542 * pkin(4) + t511 * pkin(5) + t509 * r_i_i_C(1) + t508 * r_i_i_C(2) + t517 * qJ(4) - qJD(4) * t641 + t638 * t510 + t645 * t518 + t582 * t622 + t644 * t615, t616, t521 * t600 + t532 * t597 - t535 * t587 - t588 * t652, -t521, t638 * t515 + t603 * t596 + t601 * (-t527 * qJD(5) + t653), (-t515 * t579 - t583 * t652) * r_i_i_C(1) + (-t515 * t583 + t579 * t652) * r_i_i_C(2) + ((-t527 * t583 + t651) * r_i_i_C(1) + (t527 * t579 + t650) * r_i_i_C(2)) * qJD(6); 0, 0, -t588 * t540 - t600 * t541 + t545 * t597 + t587 * t546, t541, -t638 * t523 - t642 * t596 + t601 * (-t531 * qJD(5) + t541 * t584), (t523 * t579 - t540 * t583) * r_i_i_C(1) + (t523 * t583 + t540 * t579) * r_i_i_C(2) + ((-t531 * t583 - t546 * t579) * r_i_i_C(1) + (t531 * t579 - t546 * t583) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end