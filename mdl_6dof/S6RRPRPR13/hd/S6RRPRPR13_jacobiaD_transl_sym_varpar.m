% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR13_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR13_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR13_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
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
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 10:26:25
	% EndTime: 2019-10-10 10:26:25
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
	% StartTime: 2019-10-10 10:26:25
	% EndTime: 2019-10-10 10:26:25
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
	% StartTime: 2019-10-10 10:26:27
	% EndTime: 2019-10-10 10:26:27
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (378->66), mult. (1129->103), div. (0->0), fcn. (1106->10), ass. (0->48)
	t336 = cos(pkin(6));
	t337 = sin(qJ(4));
	t340 = cos(qJ(4));
	t334 = sin(pkin(6));
	t341 = cos(qJ(2));
	t367 = t334 * t341;
	t372 = -t336 * t340 + t337 * t367;
	t338 = sin(qJ(2));
	t342 = cos(qJ(1));
	t362 = t342 * t338;
	t339 = sin(qJ(1));
	t363 = t339 * t341;
	t325 = t336 * t363 + t362;
	t347 = t336 * t362 + t363;
	t319 = t325 * qJD(1) + t347 * qJD(2);
	t361 = t342 * t341;
	t364 = t339 * t338;
	t323 = -t336 * t361 + t364;
	t366 = t334 * t342;
	t350 = t323 * t340 + t337 * t366;
	t360 = qJD(1) * t334;
	t355 = t339 * t360;
	t371 = t350 * qJD(4) + t319 * t337 + t340 * t355;
	t333 = sin(pkin(11));
	t335 = cos(pkin(11));
	t352 = t335 * r_i_i_C(1) - t333 * r_i_i_C(2) + pkin(4);
	t369 = r_i_i_C(3) + qJ(5);
	t344 = t352 * t337 - t369 * t340 + qJ(3);
	t370 = pkin(3) + pkin(8);
	t368 = t334 * t339;
	t359 = qJD(2) * t338;
	t357 = t336 * t364;
	t356 = t340 * t366;
	t354 = t342 * t360;
	t353 = t334 * t359;
	t351 = t333 * r_i_i_C(1) + t335 * r_i_i_C(2) + pkin(2) + pkin(9);
	t349 = t325 * t337 + t340 * t368;
	t348 = t325 * t340 - t337 * t368;
	t346 = t357 - t361;
	t309 = -t319 * t340 - qJD(4) * t356 + (qJD(4) * t323 + t355) * t337;
	t343 = -t340 * qJD(5) + qJD(3) + (t369 * t337 + t352 * t340) * qJD(4);
	t322 = -qJD(4) * t372 - t340 * t353;
	t320 = -qJD(1) * t357 - t339 * t359 + (qJD(2) * t336 + qJD(1)) * t361;
	t318 = t347 * qJD(1) + t325 * qJD(2);
	t317 = t323 * qJD(1) + t346 * qJD(2);
	t314 = t348 * qJD(4) - t317 * t337 + t340 * t354;
	t313 = t349 * qJD(4) + t317 * t340 + t337 * t354;
	t1 = [t350 * qJD(5) - t319 * qJ(3) - t323 * qJD(3) - t352 * t371 - t351 * t320 - t369 * t309 + (-t342 * pkin(1) - t370 * t368) * qJD(1), t351 * t317 - t318 * t344 - t343 * t346, -t317, t349 * qJD(5) - t352 * t313 + t369 * t314, t313, 0; -t348 * qJD(5) - t317 * qJ(3) + t325 * qJD(3) + t352 * t314 - t351 * t318 + t369 * t313 + (-t339 * pkin(1) + t370 * t366) * qJD(1), -t351 * t319 + t344 * t320 + t343 * t347, t319, -(-t323 * t337 + t356) * qJD(5) + t369 * t371 - t352 * t309, t309, 0; 0, (t344 * t341 * qJD(2) + (-t351 * qJD(2) + t343) * t338) * t334, t353, -t372 * qJD(5) - t369 * (-t337 * t353 + (t336 * t337 + t340 * t367) * qJD(4)) - t352 * t322, t322, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:27
	% EndTime: 2019-10-10 10:26:28
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (634->101), mult. (1584->158), div. (0->0), fcn. (1592->12), ass. (0->64)
	t398 = cos(pkin(6));
	t400 = sin(qJ(4));
	t403 = cos(qJ(4));
	t397 = sin(pkin(6));
	t404 = cos(qJ(2));
	t435 = t397 * t404;
	t380 = t398 * t403 - t400 * t435;
	t395 = pkin(11) + qJ(6);
	t393 = sin(t395);
	t394 = cos(t395);
	t418 = t394 * r_i_i_C(1) - t393 * r_i_i_C(2);
	t411 = t418 * qJD(6);
	t401 = sin(qJ(2));
	t405 = cos(qJ(1));
	t430 = t405 * t401;
	t402 = sin(qJ(1));
	t431 = t402 * t404;
	t382 = t398 * t430 + t431;
	t429 = t405 * t404;
	t432 = t402 * t401;
	t381 = -t398 * t429 + t432;
	t434 = t397 * t405;
	t424 = t403 * t434;
	t415 = -t381 * t400 + t424;
	t443 = -t382 * t394 - t393 * t415;
	t442 = -t382 * t393 + t394 * t415;
	t383 = t398 * t431 + t430;
	t369 = t383 * qJD(1) + t382 * qJD(2);
	t414 = t381 * t403 + t400 * t434;
	t428 = qJD(1) * t397;
	t423 = t402 * t428;
	t361 = t414 * qJD(4) + t369 * t400 + t403 * t423;
	t392 = cos(pkin(11)) * pkin(5) + pkin(4);
	t416 = t392 + t418;
	t440 = r_i_i_C(3) + pkin(10) + qJ(5);
	t407 = t416 * t400 - t440 * t403 + qJ(3);
	t441 = pkin(3) + pkin(8);
	t437 = t397 * t401;
	t436 = t397 * t402;
	t427 = qJD(2) * t401;
	t425 = t398 * t432;
	t422 = t405 * t428;
	t421 = qJD(2) * t435;
	t420 = t397 * t427;
	t419 = -sin(pkin(11)) * pkin(5) - pkin(2) - pkin(9);
	t417 = -t393 * r_i_i_C(1) - t394 * r_i_i_C(2);
	t374 = t383 * t400 + t403 * t436;
	t373 = t383 * t403 - t400 * t436;
	t413 = t398 * t400 + t403 * t435;
	t412 = t425 - t429;
	t410 = qJD(6) * t417;
	t408 = -t417 - t419;
	t359 = -t369 * t403 - qJD(4) * t424 + (qJD(4) * t381 + t423) * t400;
	t406 = -t403 * qJD(5) + qJD(3) + t400 * t410 + (t440 * t400 + t416 * t403) * qJD(4);
	t372 = t380 * qJD(4) - t403 * t420;
	t371 = t413 * qJD(4) - t400 * t420;
	t370 = -qJD(1) * t425 - t402 * t427 + (qJD(2) * t398 + qJD(1)) * t429;
	t368 = t382 * qJD(1) + t383 * qJD(2);
	t367 = t381 * qJD(1) + t412 * qJD(2);
	t364 = t373 * qJD(4) - t367 * t400 + t403 * t422;
	t363 = t374 * qJD(4) + t367 * t403 + t400 * t422;
	t358 = t364 * t394 - t368 * t393 + (-t374 * t393 - t394 * t412) * qJD(6);
	t357 = -t364 * t393 - t368 * t394 + (-t374 * t394 + t393 * t412) * qJD(6);
	t1 = [t414 * qJD(5) - t369 * qJ(3) - t381 * qJD(3) - t440 * t359 - t416 * t361 + (t443 * r_i_i_C(1) - t442 * r_i_i_C(2)) * qJD(6) - t408 * t370 + (-t405 * pkin(1) - t441 * t436) * qJD(1), t408 * t367 - t407 * t368 - t383 * t411 - t406 * t412, -t367, t374 * qJD(5) - t416 * t363 + t440 * t364 + t373 * t410, t363, t357 * r_i_i_C(1) - t358 * r_i_i_C(2); t358 * r_i_i_C(1) + t357 * r_i_i_C(2) - t367 * qJ(3) + t383 * qJD(3) - t373 * qJD(5) + t364 * t392 + t440 * t363 + t419 * t368 + (-pkin(1) * t402 + t441 * t434) * qJD(1), -t408 * t369 + t407 * t370 - t381 * t411 + t406 * t382, t369, -qJD(5) * t415 - t416 * t359 + t440 * t361 + t414 * t410, t359, (-t361 * t393 + t370 * t394) * r_i_i_C(1) + (-t361 * t394 - t370 * t393) * r_i_i_C(2) + (t442 * r_i_i_C(1) + t443 * r_i_i_C(2)) * qJD(6); 0, ((t407 * qJD(2) + t411) * t404 + (-t408 * qJD(2) + t406) * t401) * t397, t420, t380 * qJD(5) - t440 * t371 - t416 * t372 - t413 * t410, t372, (t371 * t393 + t394 * t421) * r_i_i_C(1) + (t371 * t394 - t393 * t421) * r_i_i_C(2) + ((-t380 * t394 - t393 * t437) * r_i_i_C(1) + (t380 * t393 - t394 * t437) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end