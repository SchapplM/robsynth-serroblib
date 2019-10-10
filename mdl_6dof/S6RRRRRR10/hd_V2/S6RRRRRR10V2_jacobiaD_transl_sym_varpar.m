% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10V2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
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
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(1);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0, 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0, 0, 0; 0, -t20, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (77->25), mult. (110->37), div. (0->0), fcn. (71->6), ass. (0->26)
	t36 = qJD(2) + qJD(3);
	t37 = qJ(2) + qJ(3);
	t35 = cos(t37);
	t55 = r_i_i_C(2) * t35;
	t34 = sin(t37);
	t57 = r_i_i_C(1) * t34;
	t46 = t55 + t57;
	t44 = t46 * t36;
	t38 = sin(qJ(2));
	t58 = pkin(2) * t38;
	t59 = qJD(2) * t58 + t44;
	t56 = r_i_i_C(2) * t34;
	t54 = t35 * t36;
	t39 = sin(qJ(1));
	t53 = qJD(1) * t39;
	t41 = cos(qJ(1));
	t52 = qJD(1) * t41;
	t40 = cos(qJ(2));
	t51 = qJD(2) * t40;
	t50 = r_i_i_C(1) * t54;
	t49 = t36 * t56;
	t48 = qJD(1) * t55;
	t45 = -t40 * pkin(2) - r_i_i_C(1) * t35 - pkin(1) + t56;
	t43 = t39 * t48 + t53 * t57 + (t49 - t50) * t41;
	t29 = t39 * t49;
	t1 = [t59 * t39 + (-r_i_i_C(3) * t39 + t45 * t41) * qJD(1), (t38 * t53 - t41 * t51) * pkin(2) + t43, t43, 0, 0, 0; -t59 * t41 + (r_i_i_C(3) * t41 + t45 * t39) * qJD(1), t29 + (-pkin(2) * t51 - t50) * t39 + (-t46 - t58) * t52, -t41 * t48 + t29 + (-t34 * t52 - t39 * t54) * r_i_i_C(1), 0, 0, 0; 0, -t59, -t44, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:03
	% EndTime: 2019-10-10 13:38:03
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (266->55), mult. (390->92), div. (0->0), fcn. (295->8), ass. (0->53)
	t257 = qJ(2) + qJ(3);
	t254 = sin(t257);
	t261 = cos(qJ(4));
	t308 = r_i_i_C(1) * t261 + pkin(3);
	t311 = t254 * t308;
	t258 = sin(qJ(4));
	t291 = qJD(4) * t261;
	t255 = cos(t257);
	t256 = qJD(2) + qJD(3);
	t299 = t255 * t256;
	t310 = t254 * t291 + t258 * t299;
	t305 = pkin(5) + r_i_i_C(3);
	t286 = t305 * t255;
	t259 = sin(qJ(2));
	t300 = pkin(2) * qJD(2);
	t288 = t259 * t300;
	t303 = pkin(3) * t254;
	t309 = -t288 + (t286 - t303) * t256;
	t292 = qJD(4) * t258;
	t279 = t254 * t292;
	t306 = r_i_i_C(1) * t279 + t310 * r_i_i_C(2);
	t304 = pkin(2) * t259;
	t301 = r_i_i_C(2) * t258;
	t260 = sin(qJ(1));
	t298 = t256 * t260;
	t297 = t256 * t261;
	t263 = cos(qJ(1));
	t296 = t256 * t263;
	t295 = t261 * t263;
	t294 = qJD(1) * t260;
	t293 = qJD(1) * t263;
	t290 = t254 * t301;
	t289 = qJD(1) * t301;
	t287 = t305 * t254;
	t285 = t305 * t260;
	t284 = t254 * t297;
	t274 = qJD(4) * t255 - qJD(1);
	t273 = qJD(1) * t255 - qJD(4);
	t272 = t308 * t255;
	t271 = t308 * t263;
	t270 = t306 * t263 + t294 * t311;
	t269 = t274 * t258;
	t268 = t263 * t254 * t289 + t306 * t260 + t293 * t286;
	t262 = cos(qJ(2));
	t267 = qJD(1) * (-t262 * pkin(2) - pkin(3) * t255 - pkin(1) - t287);
	t266 = t254 * t296 + t273 * t260;
	t265 = -t262 * t300 + (-t272 - t287) * t256;
	t264 = -t255 * r_i_i_C(2) * t291 + (-t255 * t292 - t284) * r_i_i_C(1) + t305 * t299 + (-t303 + t290) * t256;
	t238 = -t273 * t295 + (t269 + t284) * t260;
	t237 = t274 * t261 * t260 + (-t254 * t298 + t273 * t263) * t258;
	t236 = t266 * t261 + t263 * t269;
	t235 = t266 * t258 - t274 * t295;
	t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) - t309 * t260 + t263 * t267, (-t286 - t290 + t304) * t294 + t265 * t263 + t270, (-t260 * t289 - t305 * t296) * t254 + (-qJD(1) * t285 - t256 * t271) * t255 + t270, t235 * r_i_i_C(1) + t236 * r_i_i_C(2), 0, 0; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t260 * t267 + t309 * t263, (-t304 - t311) * t293 + t265 * t260 + t268, -t272 * t298 + (-qJD(1) * t271 - t256 * t285) * t254 + t268, -t237 * r_i_i_C(1) + t238 * r_i_i_C(2), 0, 0; 0, t264 - t288, t264, (-t255 * t297 + t279) * r_i_i_C(2) - t310 * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:04
	% EndTime: 2019-10-10 13:38:05
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (562->105), mult. (860->177), div. (0->0), fcn. (777->10), ass. (0->81)
	t366 = sin(qJ(5));
	t370 = cos(qJ(5));
	t367 = sin(qJ(4));
	t414 = qJD(4) * t367;
	t364 = qJD(2) + qJD(3);
	t371 = cos(qJ(4));
	t436 = -qJD(5) * t371 + t364;
	t435 = -t436 * t366 + t370 * t414;
	t438 = t366 * t414 + t436 * t370;
	t430 = r_i_i_C(3) * t367;
	t437 = pkin(3) + t430;
	t365 = qJ(2) + qJ(3);
	t362 = sin(t365);
	t369 = sin(qJ(1));
	t363 = cos(t365);
	t417 = qJD(1) * t363;
	t393 = -qJD(4) + t417;
	t373 = cos(qJ(1));
	t422 = t364 * t373;
	t434 = t362 * t422 + t393 * t369;
	t368 = sin(qJ(2));
	t432 = pkin(2) * t368;
	t431 = pkin(5) * t362;
	t429 = t362 * pkin(3);
	t428 = t363 * pkin(5);
	t427 = t366 * r_i_i_C(1);
	t426 = pkin(2) * qJD(2);
	t425 = t362 * t371;
	t424 = t363 * t364;
	t423 = t364 * t369;
	t421 = t369 * t367;
	t420 = t369 * t371;
	t419 = t373 * t367;
	t418 = t373 * t371;
	t416 = qJD(1) * t369;
	t415 = qJD(1) * t373;
	t413 = qJD(4) * t371;
	t412 = qJD(5) * t362;
	t411 = qJD(5) * t366;
	t410 = qJD(5) * t370;
	t408 = pkin(5) * t417;
	t407 = r_i_i_C(3) * t413;
	t406 = t368 * t426;
	t404 = t363 * t421;
	t403 = t362 * t423;
	t396 = t364 * t371 - qJD(5);
	t385 = t396 * t366;
	t374 = -t362 * t438 + t363 * t385;
	t386 = t396 * t370;
	t375 = t435 * t362 - t363 * t386;
	t381 = t363 * t370 + t366 * t425;
	t382 = -t363 * t366 + t370 * t425;
	t402 = (t374 * t369 + t381 * t415) * r_i_i_C(2) + (t375 * t369 - t382 * t415) * r_i_i_C(1) + t373 * t408;
	t400 = t362 * t416;
	t372 = cos(qJ(2));
	t397 = -t372 * pkin(2) - pkin(3) * t363 - pkin(1);
	t394 = -qJD(4) * t363 + qJD(1);
	t392 = (t374 * t373 - t381 * t416) * r_i_i_C(2) + (t375 * t373 + t382 * t416) * r_i_i_C(1) + t437 * t400;
	t391 = t437 * t364;
	t390 = -t370 * r_i_i_C(1) + t366 * r_i_i_C(2);
	t389 = t370 * r_i_i_C(2) + t427;
	t383 = t394 * t373;
	t345 = t367 * t383 - t434 * t371;
	t388 = t373 * t412 + t345;
	t351 = t363 * t418 + t421;
	t347 = t351 * qJD(1) - qJD(4) * t404 - t371 * t403 - t373 * t413;
	t387 = -t369 * t412 - t347;
	t384 = -pkin(5) - t389;
	t349 = t363 * t420 - t419;
	t380 = -qJD(5) * t349 + t362 * t415 + t363 * t423;
	t379 = -qJD(5) * t351 + t363 * t422 - t400;
	t378 = -t362 * t391 + (t362 * t385 + t363 * t438) * r_i_i_C(2) + (-t362 * t386 - t435 * t363) * r_i_i_C(1) + t363 * t407 + pkin(5) * t424;
	t377 = -t362 * t407 + (-t363 * t437 - t431) * t364;
	t376 = -t372 * t426 + t377;
	t350 = -t363 * t419 + t420;
	t348 = -t404 - t418;
	t346 = t394 * t420 + (-t393 * t373 + t403) * t367;
	t344 = t434 * t367 + t371 * t383;
	t337 = t379 * t366 + t388 * t370;
	t336 = -t388 * t366 + t379 * t370;
	t1 = [(-t347 * t370 + t349 * t411) * r_i_i_C(1) + (t347 * t366 + t349 * t410) * r_i_i_C(2) + t346 * r_i_i_C(3) + (t384 * t362 + t397) * t415 + (t406 + t390 * t412 + (t384 * t363 + t429) * t364) * t369, (-t428 + t432) * t416 + t376 * t373 + t392, -t369 * t408 + t377 * t373 + t392, t345 * r_i_i_C(3) + (-t344 * t366 - t350 * t410) * r_i_i_C(2) + (t344 * t370 - t350 * t411) * r_i_i_C(1), t336 * r_i_i_C(1) - t337 * r_i_i_C(2), 0; t337 * r_i_i_C(1) + t336 * r_i_i_C(2) - t344 * r_i_i_C(3) + (-t406 + (t428 - t429) * t364) * t373 + (t397 - t431) * t416, (-t362 * t437 - t432) * t415 + t376 * t369 + t402, -t369 * t363 * t391 + (-pkin(3) * t415 - pkin(5) * t423 + (-t367 * t415 - t369 * t413) * r_i_i_C(3)) * t362 + t402, t347 * r_i_i_C(3) + (-t346 * t366 - t348 * t410) * r_i_i_C(2) + (t346 * t370 - t348 * t411) * r_i_i_C(1), (t380 * r_i_i_C(1) + t387 * r_i_i_C(2)) * t370 + (t387 * r_i_i_C(1) - t380 * r_i_i_C(2)) * t366, 0; 0, t378 - t406, t378, (r_i_i_C(3) * t371 + t390 * t367) * t424 + (t389 * t367 * qJD(5) + (t390 * t371 - t430) * qJD(4)) * t362, (-r_i_i_C(2) * t386 - t396 * t427) * t363 + (t438 * r_i_i_C(1) + t435 * r_i_i_C(2)) * t362, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:06
	% EndTime: 2019-10-10 13:38:08
	% DurationCPUTime: 1.64s
	% Computational Cost: add. (1406->167), mult. (2206->278), div. (0->0), fcn. (2177->12), ass. (0->122)
	t524 = qJ(2) + qJ(3);
	t521 = sin(t524);
	t526 = sin(qJ(5));
	t531 = cos(qJ(5));
	t523 = qJD(2) + qJD(3);
	t532 = cos(qJ(4));
	t566 = qJD(5) * t532 - t523;
	t527 = sin(qJ(4));
	t592 = qJD(4) * t527;
	t621 = -t526 * t592 + t566 * t531;
	t634 = t621 * t521;
	t622 = t566 * t526 + t531 * t592;
	t633 = t622 * t521;
	t522 = cos(t524);
	t534 = cos(qJ(1));
	t596 = t534 * t527;
	t529 = sin(qJ(1));
	t599 = t529 * t532;
	t501 = t522 * t599 - t596;
	t603 = t523 * t529;
	t578 = t522 * t603;
	t593 = qJD(1) * t534;
	t550 = -t521 * t593 - t578;
	t544 = -qJD(5) * t501 - t550;
	t595 = t534 * t532;
	t576 = t522 * t595;
	t600 = t529 * t527;
	t503 = t576 + t600;
	t590 = qJD(4) * t534;
	t568 = t532 * t590;
	t569 = t529 * t592;
	t579 = t521 * t603;
	t489 = t503 * qJD(1) - t522 * t569 - t532 * t579 - t568;
	t608 = t521 * t529;
	t561 = qJD(5) * t608 + t489;
	t472 = t544 * t526 + t561 * t531;
	t591 = qJD(4) * t532;
	t548 = t527 * t593 + t529 * t591;
	t594 = qJD(1) * t529;
	t549 = t527 * t590 + t532 * t594;
	t488 = t548 * t522 - t527 * t579 - t549;
	t525 = sin(qJ(6));
	t530 = cos(qJ(6));
	t632 = t472 * t525 - t488 * t530;
	t631 = -t472 * t530 - t488 * t525;
	t620 = pkin(6) + r_i_i_C(3);
	t567 = t523 * t532 - qJD(5);
	t558 = t567 * t531;
	t587 = qJD(6) * t527;
	t570 = t521 * t587;
	t630 = t522 * t558 + t570 - t633;
	t562 = t530 * r_i_i_C(1) - t525 * r_i_i_C(2);
	t583 = t620 * t526;
	t629 = t562 * t531 + t583;
	t491 = t501 * t531 + t526 * t608;
	t500 = t522 * t600 + t595;
	t628 = t491 * t525 - t500 * t530;
	t627 = t491 * t530 + t500 * t525;
	t625 = t525 * r_i_i_C(1) + t530 * r_i_i_C(2);
	t598 = t531 * t532;
	t606 = t522 * t526;
	t498 = t521 * t598 - t606;
	t624 = qJD(6) * t498;
	t528 = sin(qJ(2));
	t614 = pkin(2) * qJD(2);
	t584 = t528 * t614;
	t617 = pkin(5) * t522;
	t618 = pkin(3) * t521;
	t623 = (-t617 + t618) * t523 + t584;
	t619 = pkin(2) * t528;
	t609 = t521 * t526;
	t607 = t522 * t523;
	t605 = t522 * t531;
	t604 = t523 * t527;
	t602 = t523 * t534;
	t601 = t527 * t531;
	t597 = t531 * t534;
	t589 = qJD(5) * t534;
	t588 = qJD(6) * t525;
	t586 = qJD(6) * t530;
	t585 = qJD(6) * t531;
	t582 = t521 * t602;
	t581 = t522 * t602;
	t580 = t523 * t597;
	t577 = t522 * t596;
	t575 = t521 * t594;
	t574 = t531 * t594;
	t571 = t521 * t589;
	t564 = -pkin(3) * t522 - pkin(5) * t521;
	t560 = -t521 * t558 + (-t622 + t587) * t522;
	t485 = t567 * t605 - t633;
	t559 = -t485 - t570;
	t557 = t567 * t526;
	t556 = -t498 * t594 + t630 * t534;
	t496 = t498 * t534;
	t555 = qJD(1) * t496 + t630 * t529;
	t494 = t503 * t531 + t534 * t609;
	t554 = t532 * t609 + t605;
	t537 = -t523 * t577 + qJD(6) * t496 + (t527 * t594 - t568) * t521;
	t553 = (t556 * t525 + t537 * t530) * r_i_i_C(2) + (t537 * t525 - t556 * t530) * r_i_i_C(1) + pkin(3) * t575 + t620 * (-t523 * t526 * t576 - t571 * t598 + (t526 * t589 + t574) * t522 + (t549 * t526 + t580) * t521);
	t538 = -t548 * t521 - t527 * t578 + t529 * t624;
	t552 = (t555 * t525 + t538 * t530) * r_i_i_C(2) + (t538 * t525 - t555 * t530) * r_i_i_C(1) + t593 * t617 - t620 * (t554 * t593 + (t522 * t557 + t634) * t529);
	t533 = cos(qJ(2));
	t551 = qJD(1) * (-t533 * pkin(2) - pkin(1) + t564);
	t542 = -qJD(6) * (t522 * t598 + t609) - t521 * t604 + t522 * t591;
	t545 = -t523 * t618 + (-t560 * t525 + t542 * t530) * r_i_i_C(2) + (t542 * t525 + t560 * t530) * r_i_i_C(1) + pkin(5) * t607 - t620 * (t521 * t557 - t522 * t621);
	t543 = t521 * t591 + t522 * t604 - t624;
	t541 = t564 * t523 - t533 * t614;
	t539 = (t562 * t526 - t620 * t531) * qJD(5);
	t471 = -t561 * t526 + t544 * t531;
	t535 = t625 * t585 + t539;
	t502 = t577 - t599;
	t493 = -t503 * t526 + t521 * t597;
	t490 = -t501 * t526 + t531 * t608;
	t487 = (-qJD(4) * t522 + qJD(1)) * t596 + (-t582 + (-qJD(1) * t522 + qJD(4)) * t529) * t532;
	t486 = t500 * qJD(1) - t522 * t568 + t527 * t582 - t569;
	t484 = -t567 * t606 - t634;
	t470 = (t487 + t571) * t531 + (-qJD(5) * t503 - t575 + t581) * t526;
	t469 = t494 * qJD(5) + t487 * t526 + t521 * t574 - t522 * t580;
	t462 = t470 * t530 - t486 * t525 + (-t494 * t525 + t502 * t530) * qJD(6);
	t461 = -t470 * t525 - t486 * t530 + (-t494 * t530 - t502 * t525) * qJD(6);
	t1 = [t631 * r_i_i_C(1) + t632 * r_i_i_C(2) + t620 * t471 + (t628 * r_i_i_C(1) + t627 * r_i_i_C(2)) * qJD(6) + t623 * t529 + t534 * t551, (-t617 + t619) * t594 + t541 * t534 + t553, -pkin(3) * t581 + (-t522 * t594 - t582) * pkin(5) + t553, (t487 * t525 + t503 * t586) * r_i_i_C(1) + (t487 * t530 - t503 * t588) * r_i_i_C(2) + t629 * t486 + t535 * t502, t620 * t470 + (t469 * t525 - t493 * t586) * r_i_i_C(2) + (-t469 * t530 - t493 * t588) * r_i_i_C(1), t461 * r_i_i_C(1) - t462 * r_i_i_C(2); t462 * r_i_i_C(1) + t461 * r_i_i_C(2) + t620 * t469 + t529 * t551 - t623 * t534, (-t618 - t619) * t593 + t541 * t529 + t552, t550 * pkin(3) - pkin(5) * t579 + t552, (t489 * t525 + t501 * t586) * r_i_i_C(1) + (t489 * t530 - t501 * t588) * r_i_i_C(2) - t629 * t488 + t535 * t500, t620 * t472 + (-t471 * t525 - t490 * t586) * r_i_i_C(2) + (t471 * t530 - t490 * t588) * r_i_i_C(1), -t632 * r_i_i_C(1) + t631 * r_i_i_C(2) + (-t627 * r_i_i_C(1) + t628 * r_i_i_C(2)) * qJD(6); 0, t545 - t584, t545, ((t525 * t532 - t530 * t601) * r_i_i_C(1) + (t525 * t601 + t530 * t532) * r_i_i_C(2) - t527 * t583) * t607 + ((-qJD(4) * t629 + t562 * qJD(6)) * t532 + (t539 + t625 * (-qJD(4) + t585)) * t527) * t521, t620 * t485 + (-t484 * t525 + t554 * t586) * r_i_i_C(2) + (t484 * t530 + t554 * t588) * r_i_i_C(1), (t543 * r_i_i_C(1) + t559 * r_i_i_C(2)) * t530 + (t559 * r_i_i_C(1) - t543 * r_i_i_C(2)) * t525;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end