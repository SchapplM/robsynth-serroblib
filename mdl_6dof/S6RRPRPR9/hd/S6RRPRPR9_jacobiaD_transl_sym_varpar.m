% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
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
	% StartTime: 2019-10-10 10:19:00
	% EndTime: 2019-10-10 10:19:00
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (103->31), mult. (307->49), div. (0->0), fcn. (282->8), ass. (0->26)
	t199 = sin(pkin(11));
	t200 = sin(pkin(6));
	t201 = cos(pkin(11));
	t220 = t200 * (r_i_i_C(1) * t199 + r_i_i_C(2) * t201 + pkin(8));
	t219 = r_i_i_C(3) + qJ(3);
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t218 = t204 * t203;
	t205 = cos(qJ(2));
	t217 = t204 * t205;
	t206 = cos(qJ(1));
	t216 = t206 * t203;
	t215 = t206 * t205;
	t214 = qJD(2) * t203;
	t202 = cos(pkin(6));
	t213 = t202 * t218;
	t212 = t202 * t215;
	t211 = qJD(2) * t202 + qJD(1);
	t210 = t201 * r_i_i_C(1) - t199 * r_i_i_C(2) + pkin(2);
	t208 = t202 * t217 + t216;
	t207 = t202 * t216 + t217;
	t194 = -qJD(1) * t213 - t204 * t214 + t211 * t215;
	t193 = t208 * qJD(1) + t207 * qJD(2);
	t192 = t207 * qJD(1) + t208 * qJD(2);
	t191 = -qJD(1) * t212 - qJD(2) * t215 + t211 * t218;
	t1 = [-(-t212 + t218) * qJD(3) - t219 * t193 - t210 * t194 + (-t206 * pkin(1) - t204 * t220) * qJD(1), -(t213 - t215) * qJD(3) - t219 * t192 + t210 * t191, -t191, 0, 0, 0; t208 * qJD(3) - t219 * t191 - t210 * t192 + (-t204 * pkin(1) + t206 * t220) * qJD(1), t207 * qJD(3) - t210 * t193 + t219 * t194, t193, 0, 0, 0; 0, (t203 * qJD(3) + (-t210 * t203 + t219 * t205) * qJD(2)) * t200, t200 * t214, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:00
	% EndTime: 2019-10-10 10:19:00
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (222->61), mult. (494->100), div. (0->0), fcn. (467->10), ass. (0->45)
	t259 = sin(qJ(1));
	t256 = cos(pkin(6));
	t267 = qJD(2) * t256 + qJD(1);
	t258 = sin(qJ(2));
	t280 = t259 * t258;
	t272 = t256 * t280;
	t275 = qJD(2) * t258;
	t260 = cos(qJ(2));
	t261 = cos(qJ(1));
	t277 = t261 * t260;
	t239 = -qJD(1) * t272 - t259 * t275 + t267 * t277;
	t255 = sin(pkin(6));
	t281 = t255 * t261;
	t268 = qJD(4) * t281;
	t285 = t239 - t268;
	t253 = pkin(11) + qJ(4);
	t251 = sin(t253);
	t252 = cos(t253);
	t266 = r_i_i_C(1) * t251 + r_i_i_C(2) * t252;
	t263 = qJD(4) * t266;
	t284 = r_i_i_C(3) + pkin(9) + qJ(3);
	t283 = t255 * t258;
	t282 = t255 * t259;
	t279 = t259 * t260;
	t278 = t261 * t258;
	t276 = qJD(1) * t255;
	t274 = qJD(2) * t260;
	t241 = t256 * t278 + t279;
	t273 = qJD(4) * t241;
	t271 = t256 * t277;
	t270 = pkin(3) * sin(pkin(11)) + pkin(8);
	t269 = t261 * t276;
	t250 = cos(pkin(11)) * pkin(3) + pkin(2);
	t265 = t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + t250;
	t264 = t256 * t279 + t278;
	t262 = t259 * t276 - t273;
	t244 = t252 * t268;
	t243 = -t272 + t277;
	t240 = -t271 + t280;
	t238 = t264 * qJD(1) + t241 * qJD(2);
	t237 = t241 * qJD(1) + t264 * qJD(2);
	t236 = -qJD(1) * t271 - t261 * t274 + t267 * t280;
	t235 = t251 * t269 - t237 * t252 + (-t243 * t251 + t252 * t282) * qJD(4);
	t234 = t252 * t269 + t237 * t251 + (-t243 * t252 - t251 * t282) * qJD(4);
	t1 = [(-t239 * t252 + t251 * t273 + t244) * r_i_i_C(1) + (t285 * t251 + t252 * t273) * r_i_i_C(2) - t239 * t250 - t240 * qJD(3) - t284 * t238 + (-t261 * pkin(1) + (-t266 - t270) * t282) * qJD(1), t243 * qJD(3) + t265 * t236 - t284 * t237 + t263 * t264, -t236, t234 * r_i_i_C(1) - t235 * r_i_i_C(2), 0, 0; t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t264 * qJD(3) - t237 * t250 - t284 * t236 + (-pkin(1) * t259 + t270 * t281) * qJD(1), t241 * qJD(3) - t265 * t238 + t284 * t239 + t240 * t263, t238, t244 * r_i_i_C(2) + (t262 * r_i_i_C(1) - t239 * r_i_i_C(2)) * t252 + (-t285 * r_i_i_C(1) - t262 * r_i_i_C(2)) * t251, 0, 0; 0, (qJD(3) * t258 - t260 * t263 + (-t265 * t258 + t284 * t260) * qJD(2)) * t255, t255 * t275, -t266 * t255 * t274 + ((-t251 * t256 - t252 * t283) * r_i_i_C(1) + (t251 * t283 - t252 * t256) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:01
	% EndTime: 2019-10-10 10:19:02
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (520->78), mult. (1089->122), div. (0->0), fcn. (1071->12), ass. (0->53)
	t344 = sin(qJ(1));
	t341 = cos(pkin(6));
	t356 = qJD(2) * t341 + qJD(1);
	t343 = sin(qJ(2));
	t373 = t344 * t343;
	t364 = t341 * t373;
	t368 = qJD(2) * t343;
	t345 = cos(qJ(2));
	t346 = cos(qJ(1));
	t370 = t346 * t345;
	t317 = -qJD(1) * t364 - t344 * t368 + t356 * t370;
	t339 = sin(pkin(6));
	t374 = t339 * t346;
	t383 = -qJD(4) * t374 + t317;
	t371 = t346 * t343;
	t372 = t344 * t345;
	t321 = t341 * t371 + t372;
	t369 = qJD(1) * t339;
	t382 = -qJD(4) * t321 + t344 * t369;
	t336 = pkin(11) + qJ(4);
	t334 = sin(t336);
	t335 = cos(t336);
	t337 = sin(pkin(12));
	t340 = cos(pkin(12));
	t355 = t340 * r_i_i_C(1) - t337 * r_i_i_C(2) + pkin(4);
	t377 = r_i_i_C(3) + qJ(5);
	t381 = (t355 * t334 - t377 * t335) * qJD(4) - t334 * qJD(5);
	t380 = t382 * t334 + t383 * t335;
	t333 = cos(pkin(11)) * pkin(3) + pkin(2);
	t378 = t377 * t334 + t355 * t335 + t333;
	t376 = t339 * t343;
	t375 = t339 * t344;
	t367 = qJD(2) * t345;
	t363 = t341 * t370;
	t362 = pkin(3) * sin(pkin(11)) + pkin(8);
	t360 = t346 * t369;
	t359 = t339 * t367;
	t342 = -pkin(9) - qJ(3);
	t354 = t337 * r_i_i_C(1) + t340 * r_i_i_C(2) - t342;
	t323 = -t364 + t370;
	t353 = -t323 * t334 + t335 * t375;
	t352 = t323 * t335 + t334 * t375;
	t351 = t341 * t334 + t335 * t376;
	t350 = t341 * t372 + t371;
	t310 = t383 * t334 - t382 * t335;
	t320 = -t363 + t373;
	t318 = t351 * qJD(4) + t334 * t359;
	t316 = t350 * qJD(1) + t321 * qJD(2);
	t315 = t321 * qJD(1) + t350 * qJD(2);
	t314 = -qJD(1) * t363 - t346 * t367 + t356 * t373;
	t309 = t353 * qJD(4) - t315 * t335 + t334 * t360;
	t308 = t352 * qJD(4) - t315 * t334 - t335 * t360;
	t1 = [(-t316 * t337 - t340 * t380) * r_i_i_C(1) + (-t316 * t340 + t337 * t380) * r_i_i_C(2) - t380 * pkin(4) - (t321 * t334 + t335 * t374) * qJD(5) - t317 * t333 + t316 * t342 - t320 * qJD(3) - t377 * t310 + (-t346 * pkin(1) - t362 * t375) * qJD(1), t323 * qJD(3) + t378 * t314 - t354 * t315 + t350 * t381, -t314, t352 * qJD(5) - t355 * t308 + t377 * t309, t308, 0; (t309 * t340 - t314 * t337) * r_i_i_C(1) + (-t309 * t337 - t314 * t340) * r_i_i_C(2) + t309 * pkin(4) - t353 * qJD(5) - t315 * t333 + t314 * t342 + t350 * qJD(3) + t377 * t308 + (-t344 * pkin(1) + t362 * t374) * qJD(1), t321 * qJD(3) - t316 * t378 + t354 * t317 + t381 * t320, t316, -(-t321 * t335 + t334 * t374) * qJD(5) + t377 * t380 - t355 * t310, t310, 0; 0, (t343 * qJD(3) - t381 * t345 + (-t343 * t378 + t354 * t345) * qJD(2)) * t339, t339 * t368, t351 * qJD(5) + t377 * (t335 * t359 + (-t334 * t376 + t335 * t341) * qJD(4)) - t355 * t318, t318, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:19:02
	% EndTime: 2019-10-10 10:19:02
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (836->103), mult. (1544->163), div. (0->0), fcn. (1557->14), ass. (0->69)
	t401 = sin(qJ(1));
	t441 = cos(pkin(6));
	t413 = qJD(2) * t441 + qJD(1);
	t400 = sin(qJ(2));
	t420 = t401 * t441;
	t417 = t400 * t420;
	t431 = qJD(2) * t400;
	t402 = cos(qJ(2));
	t403 = cos(qJ(1));
	t433 = t403 * t402;
	t363 = -qJD(1) * t417 - t401 * t431 + t413 * t433;
	t397 = sin(pkin(6));
	t435 = t397 * t403;
	t449 = -qJD(4) * t435 + t363;
	t419 = t403 * t441;
	t374 = t400 * t419 + t401 * t402;
	t432 = qJD(1) * t397;
	t448 = -qJD(4) * t374 + t401 * t432;
	t394 = pkin(11) + qJ(4);
	t390 = sin(t394);
	t392 = cos(t394);
	t367 = t374 * t392 - t390 * t435;
	t416 = t402 * t419;
	t434 = t400 * t401;
	t373 = -t416 + t434;
	t393 = pkin(12) + qJ(6);
	t389 = sin(t393);
	t391 = cos(t393);
	t447 = t367 * t389 - t373 * t391;
	t446 = t367 * t391 + t373 * t389;
	t414 = t389 * r_i_i_C(1) + t391 * r_i_i_C(2);
	t410 = qJD(6) * t414;
	t387 = cos(pkin(12)) * pkin(5) + pkin(4);
	t415 = t391 * r_i_i_C(1) - t389 * r_i_i_C(2);
	t412 = t387 + t415;
	t442 = r_i_i_C(3) + pkin(10) + qJ(5);
	t404 = (t412 * t390 - t442 * t392) * qJD(4) - t390 * qJD(5) + t392 * t410;
	t357 = t448 * t390 + t449 * t392;
	t388 = cos(pkin(11)) * pkin(3) + pkin(2);
	t443 = t442 * t390 + t412 * t392 + t388;
	t438 = t397 * t400;
	t437 = t397 * t401;
	t436 = t397 * t402;
	t430 = qJD(2) * t402;
	t427 = pkin(3) * sin(pkin(11)) + pkin(8);
	t425 = t403 * t432;
	t424 = t397 * t430;
	t422 = t397 * t431;
	t421 = -sin(pkin(12)) * pkin(5) - pkin(9) - qJ(3);
	t411 = t374 * t390 + t392 * t435;
	t376 = -t417 + t433;
	t369 = -t376 * t390 + t392 * t437;
	t370 = t376 * t392 + t390 * t437;
	t375 = t403 * t400 + t402 * t420;
	t408 = -t390 * t438 + t441 * t392;
	t372 = t441 * t390 + t392 * t438;
	t407 = t414 - t421;
	t406 = t415 * qJD(6) + qJD(3);
	t356 = t449 * t390 - t448 * t392;
	t365 = t408 * qJD(4) + t392 * t424;
	t364 = t372 * qJD(4) + t390 * t424;
	t362 = t375 * qJD(1) + t374 * qJD(2);
	t361 = t374 * qJD(1) + t375 * qJD(2);
	t360 = -qJD(1) * t416 - t403 * t430 + t413 * t434;
	t355 = t369 * qJD(4) - t361 * t392 + t390 * t425;
	t354 = t370 * qJD(4) - t361 * t390 - t392 * t425;
	t353 = t355 * t391 - t360 * t389 + (-t370 * t389 + t375 * t391) * qJD(6);
	t352 = -t355 * t389 - t360 * t391 + (-t370 * t391 - t375 * t389) * qJD(6);
	t1 = [-t411 * qJD(5) - t363 * t388 - t373 * qJD(3) - t412 * t357 - t407 * t362 - t442 * t356 + (t447 * r_i_i_C(1) + t446 * r_i_i_C(2)) * qJD(6) + (-pkin(1) * t403 - t427 * t437) * qJD(1), t443 * t360 - t407 * t361 + t404 * t375 + t406 * t376, -t360, qJD(5) * t370 - t412 * t354 + t442 * t355 - t369 * t410, t354, r_i_i_C(1) * t352 - r_i_i_C(2) * t353; t353 * r_i_i_C(1) + t352 * r_i_i_C(2) + t375 * qJD(3) - t369 * qJD(5) + t355 * t387 - t361 * t388 + t421 * t360 + t442 * t354 + (-pkin(1) * t401 + t427 * t435) * qJD(1), -t362 * t443 + t407 * t363 + t404 * t373 + t406 * t374, t362, qJD(5) * t367 - t412 * t356 + t442 * t357 + t411 * t410, t356, (-t357 * t389 + t362 * t391) * r_i_i_C(1) + (-t357 * t391 - t362 * t389) * r_i_i_C(2) + (-t446 * r_i_i_C(1) + t447 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t443 + t406) * t400 + (t407 * qJD(2) - t404) * t402) * t397, t422, qJD(5) * t372 - t412 * t364 + t442 * t365 - t408 * t410, t364, (-t365 * t389 + t391 * t422) * r_i_i_C(1) + (-t365 * t391 - t389 * t422) * r_i_i_C(2) + ((-t372 * t391 + t389 * t436) * r_i_i_C(1) + (t372 * t389 + t391 * t436) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end