% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
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
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
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
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t175 = sin(pkin(6));
	t194 = t175 * (pkin(8) + r_i_i_C(1));
	t193 = pkin(2) - r_i_i_C(2);
	t191 = r_i_i_C(3) + qJ(3);
	t177 = sin(qJ(2));
	t178 = sin(qJ(1));
	t190 = t178 * t177;
	t179 = cos(qJ(2));
	t189 = t178 * t179;
	t180 = cos(qJ(1));
	t188 = t180 * t177;
	t187 = t180 * t179;
	t186 = qJD(2) * t177;
	t176 = cos(pkin(6));
	t185 = t176 * t190;
	t184 = t176 * t187;
	t183 = qJD(2) * t176 + qJD(1);
	t182 = t176 * t189 + t188;
	t181 = t176 * t188 + t189;
	t170 = -qJD(1) * t185 - t178 * t186 + t183 * t187;
	t169 = t182 * qJD(1) + t181 * qJD(2);
	t168 = t181 * qJD(1) + t182 * qJD(2);
	t167 = -qJD(1) * t184 - qJD(2) * t187 + t183 * t190;
	t1 = [-(-t184 + t190) * qJD(3) - t193 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1), -(t185 - t187) * qJD(3) - t191 * t168 + t193 * t167, -t167, 0, 0, 0; t182 * qJD(3) - t193 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), t181 * qJD(3) - t193 * t169 + t191 * t170, t169, 0, 0, 0; 0, (t177 * qJD(3) + (-t193 * t177 + t191 * t179) * qJD(2)) * t175, t175 * t186, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:27
	% EndTime: 2019-10-10 10:50:27
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
	t275 = pkin(3) + pkin(8);
	t243 = cos(pkin(6));
	t245 = sin(qJ(2));
	t249 = cos(qJ(1));
	t267 = t249 * t245;
	t246 = sin(qJ(1));
	t248 = cos(qJ(2));
	t268 = t246 * t248;
	t235 = t243 * t268 + t267;
	t251 = t243 * t267 + t268;
	t231 = t235 * qJD(1) + t251 * qJD(2);
	t244 = sin(qJ(4));
	t274 = t231 * t244;
	t247 = cos(qJ(4));
	t273 = t231 * t247;
	t242 = sin(pkin(6));
	t272 = t242 * t246;
	t271 = t242 * t248;
	t270 = t242 * t249;
	t269 = t246 * t245;
	t266 = t249 * t248;
	t265 = qJD(1) * t246;
	t264 = qJD(1) * t249;
	t263 = qJD(2) * t245;
	t259 = t243 * t266;
	t233 = -t259 + t269;
	t262 = qJD(4) * t233;
	t261 = -r_i_i_C(3) - pkin(9) - pkin(2);
	t260 = t243 * t269;
	t258 = t242 * t265;
	t257 = t242 * t264;
	t256 = t242 * t263;
	t255 = qJD(2) * t243 + qJD(1);
	t254 = r_i_i_C(1) * t247 - r_i_i_C(2) * t244;
	t253 = -t244 * r_i_i_C(1) - t247 * r_i_i_C(2);
	t252 = qJ(3) - t253;
	t250 = t254 * qJD(4) + qJD(3);
	t232 = -qJD(1) * t260 - t246 * t263 + t255 * t266;
	t230 = t251 * qJD(1) + t235 * qJD(2);
	t229 = -qJD(1) * t259 - qJD(2) * t266 + t255 * t269;
	t228 = t247 * t257 - t229 * t244 + (t235 * t247 - t244 * t272) * qJD(4);
	t227 = -t244 * t257 - t229 * t247 + (-t235 * t244 - t247 * t272) * qJD(4);
	t1 = [(-t247 * t262 - t274) * r_i_i_C(1) + (t244 * t262 - t273) * r_i_i_C(2) - t231 * qJ(3) - t233 * qJD(3) - pkin(1) * t264 + t261 * t232 + (t253 * t249 * qJD(4) + (-t254 - t275) * t265) * t242, t250 * (-t260 + t266) - t252 * t230 - t261 * t229, -t229, t227 * r_i_i_C(1) - t228 * r_i_i_C(2), 0, 0; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t229 * qJ(3) + t235 * qJD(3) + t261 * t230 + (-pkin(1) * t246 + t275 * t270) * qJD(1), t261 * t231 + t252 * t232 + t250 * t251, t231, (-t244 * t258 + t273) * r_i_i_C(1) + (-t247 * t258 - t274) * r_i_i_C(2) + ((-t233 * t244 + t247 * t270) * r_i_i_C(1) + (-t233 * t247 - t244 * t270) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t250 * t245 + (t261 * t245 + t252 * t248) * qJD(2)) * t242, t256, t254 * t256 + ((-t243 * t247 + t244 * t271) * r_i_i_C(1) + (t243 * t244 + t247 * t271) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:29
	% EndTime: 2019-10-10 10:50:29
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (478->92), mult. (1424->150), div. (0->0), fcn. (1424->10), ass. (0->61)
	t380 = cos(pkin(6));
	t387 = cos(qJ(2));
	t388 = cos(qJ(1));
	t413 = t388 * t387;
	t383 = sin(qJ(2));
	t384 = sin(qJ(1));
	t416 = t384 * t383;
	t370 = -t380 * t413 + t416;
	t382 = sin(qJ(4));
	t386 = cos(qJ(4));
	t379 = sin(pkin(6));
	t417 = t379 * t388;
	t409 = t386 * t417;
	t366 = -t370 * t382 + t409;
	t414 = t388 * t383;
	t415 = t384 * t387;
	t371 = t380 * t414 + t415;
	t381 = sin(qJ(5));
	t385 = cos(qJ(5));
	t428 = -t366 * t381 - t371 * t385;
	t427 = t366 * t385 - t371 * t381;
	t403 = t385 * r_i_i_C(1) - t381 * r_i_i_C(2);
	t394 = t403 * qJD(5);
	t372 = t380 * t415 + t414;
	t359 = t372 * qJD(1) + t371 * qJD(2);
	t412 = qJD(1) * t379;
	t408 = t384 * t412;
	t426 = (qJD(4) * t370 + t408) * t382 - qJD(4) * t409 - t359 * t386;
	t399 = t370 * t386 + t382 * t417;
	t352 = t399 * qJD(4) + t359 * t382 + t386 * t408;
	t401 = pkin(4) + t403;
	t423 = -r_i_i_C(3) - pkin(10);
	t390 = t401 * t382 + t423 * t386 + qJ(3);
	t425 = -pkin(2) - pkin(9);
	t424 = pkin(3) + pkin(8);
	t420 = t379 * t383;
	t419 = t379 * t384;
	t418 = t379 * t387;
	t411 = qJD(2) * t383;
	t410 = t380 * t416;
	t407 = t388 * t412;
	t406 = qJD(2) * t418;
	t405 = t379 * t411;
	t402 = -t381 * r_i_i_C(1) - t385 * r_i_i_C(2);
	t400 = -t402 - t425;
	t364 = t372 * t382 + t386 * t419;
	t398 = t372 * t386 - t382 * t419;
	t397 = t380 * t382 + t386 * t418;
	t396 = -t380 * t386 + t382 * t418;
	t395 = t410 - t413;
	t393 = qJD(5) * t402;
	t389 = qJD(3) + t382 * t393 + (-t423 * t382 + t401 * t386) * qJD(4);
	t361 = t397 * qJD(4) - t382 * t405;
	t360 = -qJD(1) * t410 - t384 * t411 + (qJD(2) * t380 + qJD(1)) * t413;
	t358 = t371 * qJD(1) + t372 * qJD(2);
	t357 = t370 * qJD(1) + t395 * qJD(2);
	t355 = t398 * qJD(4) - t357 * t382 + t386 * t407;
	t354 = t364 * qJD(4) + t357 * t386 + t382 * t407;
	t349 = t355 * t385 - t358 * t381 + (-t364 * t381 - t385 * t395) * qJD(5);
	t348 = -t355 * t381 - t358 * t385 + (-t364 * t385 + t381 * t395) * qJD(5);
	t1 = [-t359 * qJ(3) - t370 * qJD(3) - t401 * t352 - t400 * t360 + t423 * t426 + (t428 * r_i_i_C(1) - t427 * r_i_i_C(2)) * qJD(5) + (-t388 * pkin(1) - t424 * t419) * qJD(1), t400 * t357 - t390 * t358 - t372 * t394 - t389 * t395, -t357, -t401 * t354 - t423 * t355 + t398 * t393, t348 * r_i_i_C(1) - t349 * r_i_i_C(2), 0; t355 * pkin(4) + t349 * r_i_i_C(1) + t348 * r_i_i_C(2) - t357 * qJ(3) + t372 * qJD(3) + t425 * t358 - t423 * t354 + (-pkin(1) * t384 + t424 * t417) * qJD(1), -t400 * t359 + t390 * t360 - t370 * t394 + t389 * t371, t359, -t423 * t352 + t399 * t393 - t401 * t426, (-t352 * t381 + t360 * t385) * r_i_i_C(1) + (-t352 * t385 - t360 * t381) * r_i_i_C(2) + (t427 * r_i_i_C(1) + t428 * r_i_i_C(2)) * qJD(5), 0; 0, ((t390 * qJD(2) + t394) * t387 + (-t400 * qJD(2) + t389) * t383) * t379, t405, t423 * t361 - t397 * t393 + t401 * (t396 * qJD(4) + t386 * t405), (t361 * t381 + t385 * t406) * r_i_i_C(1) + (t361 * t385 - t381 * t406) * r_i_i_C(2) + ((-t381 * t420 + t385 * t396) * r_i_i_C(1) + (-t381 * t396 - t385 * t420) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:31
	% EndTime: 2019-10-10 10:50:32
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (827->124), mult. (2432->190), div. (0->0), fcn. (2503->10), ass. (0->76)
	t468 = cos(pkin(6));
	t471 = sin(qJ(2));
	t476 = cos(qJ(1));
	t514 = t476 * t471;
	t472 = sin(qJ(1));
	t475 = cos(qJ(2));
	t515 = t472 * t475;
	t456 = t468 * t514 + t515;
	t457 = t468 * t515 + t514;
	t444 = t457 * qJD(1) + t456 * qJD(2);
	t470 = sin(qJ(4));
	t474 = cos(qJ(4));
	t513 = t476 * t475;
	t516 = t472 * t471;
	t455 = -t468 * t513 + t516;
	t467 = sin(pkin(6));
	t518 = t467 * t476;
	t488 = t455 * t474 + t470 * t518;
	t512 = qJD(1) * t472;
	t505 = t467 * t512;
	t431 = t488 * qJD(4) + t444 * t470 + t474 * t505;
	t511 = qJD(2) * t471;
	t500 = t472 * t511;
	t503 = t471 * t512;
	t445 = -t468 * t503 - t500 + (qJD(2) * t468 + qJD(1)) * t513;
	t469 = sin(qJ(5));
	t473 = cos(qJ(5));
	t506 = t474 * t518;
	t451 = -t455 * t470 + t506;
	t532 = t451 * t473 - t456 * t469;
	t537 = t532 * qJD(5) - t431 * t469 + t445 * t473;
	t533 = t451 * t469 + t456 * t473;
	t536 = t533 * qJD(5) + t431 * t473 + t445 * t469;
	t531 = (qJD(4) * t455 + t505) * t470 - qJD(4) * t506 - t444 * t474;
	t526 = r_i_i_C(2) + pkin(10);
	t530 = -t470 * pkin(4) + t526 * t474 - qJ(3);
	t525 = r_i_i_C(3) + qJ(6);
	t527 = r_i_i_C(1) + pkin(5);
	t481 = t525 * t469 + t527 * t473 + pkin(4);
	t529 = pkin(3) + pkin(8);
	t528 = pkin(9) + pkin(2);
	t520 = t467 * t472;
	t519 = t467 * t475;
	t517 = t471 * t473;
	t510 = qJD(4) * t474;
	t509 = qJD(5) * t470;
	t508 = qJD(6) * t469;
	t507 = t473 * qJD(6);
	t504 = qJD(1) * t518;
	t502 = qJD(2) * t519;
	t501 = t467 * t511;
	t497 = qJD(2) + t509;
	t442 = t468 * t500 + t503 + (-qJD(1) * t468 - qJD(2)) * t513;
	t458 = -t468 * t516 + t513;
	t496 = t458 * t509 - t442;
	t495 = t456 * t509 + t444;
	t449 = t457 * t470 + t474 * t520;
	t494 = t449 * t473 + t458 * t469;
	t493 = -t449 * t469 + t458 * t473;
	t490 = (qJD(2) * t470 + qJD(5)) * t475;
	t485 = -t468 * t474 + t470 * t519;
	t489 = t467 * t471 * t469 - t473 * t485;
	t487 = t457 * t474 - t470 * t520;
	t486 = t468 * t470 + t474 * t519;
	t443 = t456 * qJD(1) + t457 * qJD(2);
	t480 = -qJD(5) * t457 - t443 * t470 + t458 * t510;
	t479 = -qJD(5) * t455 + t445 * t470 + t456 * t510;
	t478 = t508 + (-t527 * t469 + t525 * t473) * qJD(5);
	t477 = t470 * t508 + qJD(3) + (pkin(4) * t474 + t526 * t470) * qJD(4);
	t446 = t486 * qJD(4) - t470 * t501;
	t435 = t489 * qJD(5) - t446 * t469 - t473 * t502;
	t434 = t487 * qJD(4) - t442 * t470 + t474 * t504;
	t433 = t449 * qJD(4) + t442 * t474 + t470 * t504;
	t424 = t493 * qJD(5) + t434 * t473 - t443 * t469;
	t423 = t494 * qJD(5) + t434 * t469 + t443 * t473;
	t1 = [t533 * qJD(6) - t431 * pkin(4) - t444 * qJ(3) - t455 * qJD(3) - t528 * t445 - t526 * t531 - t527 * t536 + t525 * t537 + (-t476 * pkin(1) - t529 * t520) * qJD(1), t457 * t507 + t528 * t442 + t527 * (-t496 * t469 + t480 * t473) + t525 * (t480 * t469 + t496 * t473) + t477 * t458 + t530 * t443, -t442, -t481 * t433 + t526 * t434 + t478 * t487, t494 * qJD(6) - t527 * t423 + t525 * t424, t423; -t493 * qJD(6) + t434 * pkin(4) - t442 * qJ(3) + t457 * qJD(3) - t528 * t443 + t526 * t433 + t527 * t424 + t525 * t423 + (-t472 * pkin(1) + t529 * t518) * qJD(1), t455 * t507 - t528 * t444 + t527 * (-t495 * t469 + t479 * t473) + t525 * (t479 * t469 + t495 * t473) + t477 * t456 - t530 * t445, t444, t526 * t431 + t478 * t488 - t481 * t531, -t532 * qJD(6) + t525 * t536 + t527 * t537, -t537; 0, (t527 * (t473 * t490 + (-t497 * t469 + t473 * t510) * t471) + t525 * (t497 * t517 + (t471 * t510 + t490) * t469) - t475 * t507 + t477 * t471 + (-t528 * t471 - t475 * t530) * qJD(2)) * t467, t501, -t526 * t446 - t478 * t486 + t481 * (t485 * qJD(4) + t474 * t501), t489 * qJD(6) + t525 * (t469 * t502 - t446 * t473 + (t467 * t517 + t469 * t485) * qJD(5)) - t527 * t435, t435;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end