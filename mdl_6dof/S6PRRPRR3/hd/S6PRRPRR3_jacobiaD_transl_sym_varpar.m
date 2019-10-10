% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(12));
	t50 = sin(pkin(12));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (95->48), mult. (333->96), div. (0->0), fcn. (328->10), ass. (0->36)
	t218 = sin(pkin(7));
	t248 = (pkin(9) + r_i_i_C(3)) * t218;
	t217 = sin(pkin(12));
	t220 = cos(pkin(12));
	t224 = sin(qJ(2));
	t222 = cos(pkin(6));
	t226 = cos(qJ(2));
	t241 = t222 * t226;
	t211 = -t217 * t224 + t220 * t241;
	t219 = sin(pkin(6));
	t245 = t218 * t219;
	t221 = cos(pkin(7));
	t223 = sin(qJ(3));
	t244 = t221 * t223;
	t225 = cos(qJ(3));
	t243 = t221 * t225;
	t242 = t222 * t224;
	t240 = t223 * t224;
	t239 = t223 * t226;
	t238 = t224 * t225;
	t237 = t225 * t226;
	t236 = t223 * t245;
	t235 = t225 * t245;
	t233 = r_i_i_C(1) * t223 + r_i_i_C(2) * t225;
	t232 = r_i_i_C(1) * t225 - r_i_i_C(2) * t223 + pkin(2);
	t212 = t217 * t226 + t220 * t242;
	t231 = t217 * t241 + t220 * t224;
	t230 = t217 * t242 - t220 * t226;
	t229 = t233 * t221 - t248;
	t228 = (-t221 * t238 - t239) * r_i_i_C(1) + (t221 * t240 - t237) * r_i_i_C(2);
	t227 = (-t221 * t239 - t238) * r_i_i_C(1) + (-t221 * t237 + t240) * r_i_i_C(2);
	t210 = t230 * qJD(2);
	t209 = t231 * qJD(2);
	t208 = t212 * qJD(2);
	t207 = t211 * qJD(2);
	t1 = [0, t232 * t210 + t229 * t209 + ((t223 * t231 + t230 * t243) * r_i_i_C(1) + (t225 * t231 - t230 * t244) * r_i_i_C(2)) * qJD(3), (t209 * t223 + t210 * t243) * r_i_i_C(1) + (t209 * t225 - t210 * t244) * r_i_i_C(2) + ((-t217 * t236 + t225 * t230 + t231 * t244) * r_i_i_C(1) + (-t217 * t235 - t223 * t230 + t231 * t243) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, -t232 * t208 - t229 * t207 + ((-t211 * t223 - t212 * t243) * r_i_i_C(1) + (-t211 * t225 + t212 * t244) * r_i_i_C(2)) * qJD(3), (-t207 * t223 - t208 * t243) * r_i_i_C(1) + (-t207 * t225 + t208 * t244) * r_i_i_C(2) + ((-t211 * t244 - t212 * t225 + t220 * t236) * r_i_i_C(1) + (-t211 * t243 + t212 * t223 + t220 * t235) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (t228 * qJD(3) + (-t224 * pkin(2) + t226 * t248 + t227) * qJD(2)) * t219, -t233 * t222 * t218 * qJD(3) + (t228 * qJD(2) + t227 * qJD(3)) * t219, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (198->86), mult. (654->154), div. (0->0), fcn. (660->12), ass. (0->49)
	t275 = sin(pkin(12));
	t277 = sin(pkin(6));
	t304 = t275 * t277;
	t276 = sin(pkin(7));
	t303 = t276 * t277;
	t279 = cos(pkin(12));
	t302 = t277 * t279;
	t280 = cos(pkin(7));
	t284 = cos(qJ(3));
	t301 = t280 * t284;
	t281 = cos(pkin(6));
	t283 = sin(qJ(2));
	t300 = t281 * t283;
	t285 = cos(qJ(2));
	t299 = t281 * t285;
	t282 = sin(qJ(3));
	t298 = t282 * t285;
	t297 = t283 * t284;
	t296 = qJD(2) * t283;
	t295 = qJD(3) * t282;
	t294 = qJD(3) * t284;
	t293 = pkin(3) * t295;
	t292 = t279 * t299;
	t291 = t280 * t294;
	t274 = sin(pkin(13));
	t278 = cos(pkin(13));
	t290 = t284 * t274 + t282 * t278;
	t265 = t282 * t274 - t284 * t278;
	t259 = t275 * t285 + t279 * t300;
	t289 = t275 * t299 + t279 * t283;
	t288 = t275 * t300 - t279 * t285;
	t287 = qJD(3) * t290;
	t262 = t265 * qJD(3);
	t251 = t265 * t280;
	t252 = t290 * t280;
	t286 = t280 * t282 * pkin(3) + t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + (-r_i_i_C(3) - pkin(9) - qJ(4)) * t276;
	t273 = t284 * pkin(3) + pkin(2);
	t264 = pkin(3) * t291 - t276 * qJD(4);
	t263 = -t274 * t294 - t278 * t295;
	t258 = -t275 * t283 + t292;
	t256 = t288 * qJD(2);
	t255 = t289 * qJD(2);
	t254 = t259 * qJD(2);
	t253 = -qJD(2) * t292 + t275 * t296;
	t250 = t280 * t287;
	t249 = t280 * t274 * t295 - t278 * t291;
	t248 = t276 * t287;
	t247 = t276 * t262;
	t1 = [0, (-t249 * t288 - t256 * t265 - t263 * t289) * r_i_i_C(1) + (-t250 * t288 - t256 * t290 - t262 * t289) * r_i_i_C(2) + t256 * t273 + t289 * t293 + t288 * t264 + t286 * t255, (-t248 * t304 + t250 * t289 - t256 * t251 + t255 * t290 - t262 * t288) * r_i_i_C(1) + (t247 * t304 - t249 * t289 - t256 * t252 - t255 * t265 + t263 * t288) * r_i_i_C(2) + (t256 * t301 + t255 * t282 + (t288 * t284 + (-t275 * t303 + t280 * t289) * t282) * qJD(3)) * pkin(3), -t256 * t276, 0, 0; 0, (t259 * t249 + t254 * t265 + t258 * t263) * r_i_i_C(1) + (t259 * t250 + t254 * t290 + t258 * t262) * r_i_i_C(2) - t254 * t273 - t258 * t293 - t259 * t264 + t286 * t253, (t248 * t302 - t258 * t250 + t254 * t251 + t253 * t290 + t259 * t262) * r_i_i_C(1) + (-t247 * t302 + t258 * t249 + t254 * t252 - t253 * t265 - t259 * t263) * r_i_i_C(2) + (-t254 * t301 + t253 * t282 + (-t259 * t284 + (-t258 * t280 + t276 * t302) * t282) * qJD(3)) * pkin(3), t254 * t276, 0, 0; 0, ((t249 * t283 + t263 * t285) * r_i_i_C(1) + (t250 * t283 + t262 * t285) * r_i_i_C(2) - t285 * t293 - t283 * t264 + ((t265 * r_i_i_C(1) + r_i_i_C(2) * t290 - t273) * t283 - t286 * t285) * qJD(2)) * t277, (-t248 * r_i_i_C(1) + t247 * r_i_i_C(2) - t276 * t293) * t281 + ((-t250 * t285 + t262 * t283) * r_i_i_C(1) + (t249 * t285 - t263 * t283) * r_i_i_C(2) + ((t251 * t283 - t285 * t290) * r_i_i_C(1) + (t252 * t283 + t265 * t285) * r_i_i_C(2)) * qJD(2) + ((-t280 * t298 - t297) * qJD(3) + (-t280 * t297 - t298) * qJD(2)) * pkin(3)) * t277, t296 * t303, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:29
	% EndTime: 2019-10-09 22:29:31
	% DurationCPUTime: 1.49s
	% Computational Cost: add. (642->156), mult. (2061->282), div. (0->0), fcn. (2236->14), ass. (0->84)
	t449 = sin(pkin(7));
	t447 = sin(pkin(13));
	t495 = cos(pkin(13));
	t497 = cos(qJ(3));
	t471 = t497 * t495;
	t455 = sin(qJ(3));
	t482 = qJD(3) * t455;
	t500 = -qJD(3) * t471 + t447 * t482;
	t416 = t500 * t449;
	t452 = cos(pkin(7));
	t418 = t500 * t452;
	t473 = t455 * t495;
	t474 = qJD(3) * t497;
	t435 = -qJD(3) * t473 - t447 * t474;
	t450 = sin(pkin(6));
	t453 = cos(pkin(6));
	t456 = sin(qJ(2));
	t458 = cos(qJ(2));
	t438 = -t497 * t447 - t473;
	t423 = t438 * t452;
	t461 = -t455 * t447 + t471;
	t467 = -t423 * t456 - t458 * t461;
	t501 = t453 * t416 + (t467 * qJD(2) + t418 * t458 - t435 * t456) * t450;
	t498 = r_i_i_C(3) + pkin(10);
	t496 = pkin(3) * t452;
	t448 = sin(pkin(12));
	t494 = t448 * t450;
	t493 = t449 * t450;
	t454 = sin(qJ(5));
	t492 = t449 * t454;
	t457 = cos(qJ(5));
	t491 = t449 * t457;
	t451 = cos(pkin(12));
	t490 = t450 * t451;
	t489 = t450 * t452;
	t487 = t453 * t456;
	t486 = t453 * t458;
	t485 = t455 * t458;
	t484 = qJD(2) * t456;
	t483 = qJD(2) * t458;
	t481 = qJD(5) * t454;
	t480 = qJD(5) * t457;
	t479 = pkin(3) * t482;
	t478 = t451 * t486;
	t477 = t452 * t497;
	t476 = t497 * t456;
	t472 = t484 * t493;
	t422 = t461 * t452;
	t469 = t422 * t458 + t438 * t456;
	t468 = -t423 * t458 + t456 * t461;
	t466 = t457 * r_i_i_C(1) - t454 * r_i_i_C(2) + pkin(4);
	t431 = t448 * t458 + t451 * t487;
	t464 = t448 * t486 + t451 * t456;
	t463 = t448 * t487 - t451 * t458;
	t462 = qJD(5) * (-t454 * r_i_i_C(1) - t457 * r_i_i_C(2));
	t460 = qJD(3) * t438;
	t424 = -qJD(2) * t478 + t448 * t484;
	t425 = t431 * qJD(2);
	t430 = -t448 * t456 + t478;
	t392 = t416 * t490 - t430 * t418 + t425 * t423 - t424 * t461 + t431 * t435;
	t426 = t464 * qJD(2);
	t427 = t463 * qJD(2);
	t459 = t416 * t494 - t418 * t464 + t427 * t423 + t426 * t461 + t435 * t463;
	t446 = t497 * pkin(3) + pkin(2);
	t436 = -t449 * qJD(4) + t474 * t496;
	t434 = t461 * qJD(3);
	t429 = t453 * t452 - t458 * t493;
	t428 = t455 * t496 + (-pkin(9) - qJ(4)) * t449;
	t421 = t438 * t449;
	t420 = t461 * t449;
	t419 = t452 * t460;
	t417 = t449 * t460;
	t415 = t448 * t489 + t449 * t464;
	t414 = -t430 * t449 - t451 * t489;
	t413 = t467 * t450;
	t412 = -t423 * t463 - t461 * t464;
	t411 = t431 * t423 + t430 * t461;
	t410 = -t453 * t421 + t468 * t450;
	t408 = -t421 * t494 + t423 * t464 - t461 * t463;
	t406 = t421 * t490 - t430 * t423 + t431 * t461;
	t404 = (-t468 * qJD(2) + t418 * t456 + t435 * t458) * t450;
	t399 = -t418 * t463 - t426 * t423 + t427 * t461 - t435 * t464;
	t397 = t431 * t418 - t424 * t423 - t425 * t461 + t430 * t435;
	t1 = [0, (t399 * t457 - t426 * t492) * r_i_i_C(1) + (-t399 * t454 - t426 * t491) * r_i_i_C(2) + t399 * pkin(4) + t427 * t446 + t464 * t479 + t426 * t428 + t463 * t436 - t498 * (t419 * t463 + t426 * t422 + t427 * t438 + t434 * t464) + ((-t412 * t454 - t463 * t491) * r_i_i_C(1) + (-t412 * t457 + t463 * t492) * r_i_i_C(2)) * qJD(5), -t498 * t459 + (t420 * t494 - t422 * t464 - t438 * t463) * t462 + t466 * (t417 * t494 - t419 * t464 + t427 * t422 - t426 * t438 + t434 * t463) + (t427 * t477 + t426 * t455 + (t497 * t463 + (-t448 * t493 + t452 * t464) * t455) * qJD(3)) * pkin(3), -t427 * t449, (-t427 * t491 + t454 * t459) * r_i_i_C(1) + (t427 * t492 + t457 * t459) * r_i_i_C(2) + ((-t408 * t457 - t415 * t454) * r_i_i_C(1) + (t408 * t454 - t415 * t457) * r_i_i_C(2)) * qJD(5), 0; 0, (t397 * t457 - t424 * t492) * r_i_i_C(1) + (-t397 * t454 - t424 * t491) * r_i_i_C(2) + t397 * pkin(4) - t425 * t446 - t430 * t479 + t424 * t428 - t431 * t436 - t498 * (-t431 * t419 + t424 * t422 - t425 * t438 - t430 * t434) + ((-t411 * t454 + t431 * t491) * r_i_i_C(1) + (-t411 * t457 - t431 * t492) * r_i_i_C(2)) * qJD(5), t498 * t392 + (-t420 * t490 + t430 * t422 + t431 * t438) * t462 + t466 * (-t417 * t490 + t430 * t419 - t425 * t422 - t424 * t438 - t431 * t434) + (-t425 * t477 + t424 * t455 + (-t497 * t431 + (-t430 * t452 + t449 * t490) * t455) * qJD(3)) * pkin(3), t425 * t449, (-t392 * t454 + t425 * t491) * r_i_i_C(1) + (-t392 * t457 - t425 * t492) * r_i_i_C(2) + ((-t406 * t457 - t414 * t454) * r_i_i_C(1) + (t406 * t454 - t414 * t457) * r_i_i_C(2)) * qJD(5), 0; 0, (t404 * t457 + t413 * t481) * r_i_i_C(1) + (-t404 * t454 + t413 * t480) * r_i_i_C(2) + t404 * pkin(4) + (t498 * (t469 * qJD(2) + t419 * t456 + t434 * t458) - t458 * t479 - t456 * t436 + (-t458 * t428 - t456 * t446) * qJD(2) + ((t454 * t483 + t456 * t480) * r_i_i_C(1) + (-t456 * t481 + t457 * t483) * r_i_i_C(2)) * t449) * t450, -t498 * t501 + (t453 * t420 + t469 * t450) * t462 + t466 * (t453 * t417 + (t419 * t458 - t434 * t456 + (-t422 * t456 + t438 * t458) * qJD(2)) * t450) + (-t449 * t453 * t482 + ((-t452 * t485 - t476) * qJD(3) + (-t452 * t476 - t485) * qJD(2)) * t450) * pkin(3), t472, (t454 * t501 + t457 * t472) * r_i_i_C(1) + (-t454 * t472 + t457 * t501) * r_i_i_C(2) + ((-t410 * t457 - t429 * t454) * r_i_i_C(1) + (t410 * t454 - t429 * t457) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:33
	% EndTime: 2019-10-09 22:29:35
	% DurationCPUTime: 1.80s
	% Computational Cost: add. (1830->227), mult. (5726->394), div. (0->0), fcn. (6521->16), ass. (0->121)
	t601 = sin(pkin(13));
	t661 = cos(pkin(13));
	t663 = cos(qJ(3));
	t633 = t663 * t661;
	t610 = sin(qJ(3));
	t648 = qJD(3) * t610;
	t665 = -qJD(3) * t633 + t601 * t648;
	t664 = r_i_i_C(3) + pkin(11);
	t606 = cos(pkin(7));
	t662 = pkin(3) * t606;
	t636 = t610 * t661;
	t637 = t663 * qJD(3);
	t589 = -qJD(3) * t636 - t601 * t637;
	t611 = sin(qJ(2));
	t660 = t589 * t611;
	t602 = sin(pkin(12));
	t604 = sin(pkin(6));
	t659 = t602 * t604;
	t603 = sin(pkin(7));
	t658 = t603 * t604;
	t609 = sin(qJ(5));
	t657 = t603 * t609;
	t613 = cos(qJ(5));
	t656 = t603 * t613;
	t605 = cos(pkin(12));
	t655 = t604 * t605;
	t654 = t604 * t606;
	t570 = t665 * t603;
	t607 = cos(pkin(6));
	t653 = t607 * t570;
	t652 = t607 * t611;
	t614 = cos(qJ(2));
	t651 = t607 * t614;
	t650 = t610 * t614;
	t649 = qJD(2) * t611;
	t608 = sin(qJ(6));
	t647 = qJD(6) * t608;
	t612 = cos(qJ(6));
	t646 = qJD(6) * t612;
	t645 = pkin(3) * t648;
	t644 = t611 * t658;
	t643 = t614 * t658;
	t642 = t605 * t651;
	t641 = t606 * t663;
	t640 = t663 * t611;
	t639 = t604 * t649;
	t635 = qJD(2) * t643;
	t634 = t603 * t639;
	t584 = -t602 * t611 + t642;
	t563 = -t584 * t603 - t605 * t654;
	t592 = -t663 * t601 - t636;
	t575 = t592 * t603;
	t577 = t592 * t606;
	t585 = t602 * t614 + t605 * t652;
	t621 = -t610 * t601 + t633;
	t620 = t575 * t655 - t584 * t577 + t585 * t621;
	t535 = t563 * t609 + t613 * t620;
	t632 = t563 * t613 - t609 * t620;
	t624 = t602 * t651 + t605 * t611;
	t564 = t602 * t654 + t603 * t624;
	t623 = t602 * t652 - t605 * t614;
	t619 = -t575 * t659 + t577 * t624 - t621 * t623;
	t537 = t564 * t609 + t613 * t619;
	t631 = t564 * t613 - t609 * t619;
	t583 = t607 * t606 - t643;
	t628 = -t577 * t614 + t611 * t621;
	t617 = -t607 * t575 + t628 * t604;
	t552 = t583 * t609 + t613 * t617;
	t630 = t583 * t613 - t609 * t617;
	t576 = t621 * t606;
	t629 = t576 * t614 + t592 * t611;
	t627 = t577 * t611 + t614 * t621;
	t626 = t612 * r_i_i_C(1) - t608 * r_i_i_C(2) + pkin(5);
	t557 = t585 * t577 + t584 * t621;
	t543 = t557 * t613 + t585 * t657;
	t559 = -t577 * t623 - t621 * t624;
	t544 = t559 * t613 - t623 * t657;
	t622 = qJD(6) * (-t608 * r_i_i_C(1) - t612 * r_i_i_C(2));
	t562 = t627 * t604;
	t560 = t562 * t613 + t609 * t644;
	t618 = qJD(3) * t592;
	t616 = t664 * t609 + t626 * t613 + pkin(4);
	t572 = t665 * t606;
	t578 = -qJD(2) * t642 + t602 * t649;
	t579 = t585 * qJD(2);
	t526 = t570 * t655 - t584 * t572 + t579 * t577 - t578 * t621 + t585 * t589;
	t580 = t624 * qJD(2);
	t581 = t623 * qJD(2);
	t527 = t570 * t659 - t572 * t624 + t581 * t577 + t580 * t621 + t589 * t623;
	t615 = t613 * t622 + (-t626 * t609 + t664 * t613) * qJD(5);
	t600 = t663 * pkin(3) + pkin(2);
	t590 = -t603 * qJD(4) + t637 * t662;
	t588 = t621 * qJD(3);
	t582 = t610 * t662 + (-pkin(9) - qJ(4)) * t603;
	t574 = t621 * t603;
	t573 = t606 * t618;
	t571 = t603 * t618;
	t561 = (t576 * t611 - t592 * t614) * t604;
	t558 = -t576 * t623 + t592 * t624;
	t556 = t585 * t576 - t584 * t592;
	t554 = t607 * t574 + t629 * t604;
	t549 = t574 * t659 - t576 * t624 - t592 * t623;
	t546 = -t574 * t655 + t584 * t576 + t585 * t592;
	t542 = (-t628 * qJD(2) + t572 * t611 + t589 * t614) * t604;
	t541 = (t629 * qJD(2) + t573 * t611 + t588 * t614) * t604;
	t540 = -t653 + (t627 * qJD(2) - t572 * t614 + t660) * t604;
	t539 = t607 * t571 - t576 * t639 + (-t588 * t611 + (qJD(2) * t592 + t573) * t614) * t604;
	t538 = t653 - t577 * t639 + (-t660 + (-qJD(2) * t621 + t572) * t614) * t604;
	t533 = -t572 * t623 - t580 * t577 + t581 * t621 - t589 * t624;
	t532 = t573 * t623 + t580 * t576 + t581 * t592 + t588 * t624;
	t531 = t585 * t572 - t578 * t577 - t579 * t621 + t584 * t589;
	t530 = -t585 * t573 + t578 * t576 - t579 * t592 - t584 * t588;
	t528 = t571 * t659 - t573 * t624 + t581 * t576 - t580 * t592 + t588 * t623;
	t525 = -t571 * t655 + t584 * t573 - t579 * t576 - t578 * t592 - t585 * t588;
	t523 = t609 * t635 + t542 * t613 + (-t562 * t609 + t613 * t644) * qJD(5);
	t521 = t630 * qJD(5) + t540 * t613 + t609 * t634;
	t519 = -t580 * t657 + t533 * t613 + (-t559 * t609 - t623 * t656) * qJD(5);
	t517 = -t578 * t657 + t531 * t613 + (-t557 * t609 + t585 * t656) * qJD(5);
	t515 = t631 * qJD(5) - t527 * t613 - t581 * t657;
	t513 = t632 * qJD(5) + t526 * t613 + t579 * t657;
	t1 = [0, (t519 * t612 - t532 * t608) * r_i_i_C(1) + (-t519 * t608 - t532 * t612) * r_i_i_C(2) + t519 * pkin(5) + t533 * pkin(4) - t532 * pkin(10) + t581 * t600 + t624 * t645 + t580 * t582 + t623 * t590 + t664 * (t544 * qJD(5) + t533 * t609 + t580 * t656) + ((-t544 * t608 + t558 * t612) * r_i_i_C(1) + (-t544 * t612 - t558 * t608) * r_i_i_C(2)) * qJD(6), (-t527 * t608 + t619 * t646) * r_i_i_C(1) + (-t527 * t612 - t619 * t647) * r_i_i_C(2) - t527 * pkin(10) + t616 * t528 + t615 * t549 + (t581 * t641 + t580 * t610 + (t663 * t623 + (-t602 * t658 + t606 * t624) * t610) * qJD(3)) * pkin(3), -t581 * t603, t664 * t515 + t631 * t622 + t626 * (-t537 * qJD(5) + t527 * t609 - t581 * t656), (-t515 * t608 - t528 * t612) * r_i_i_C(1) + (-t515 * t612 + t528 * t608) * r_i_i_C(2) + ((-t537 * t612 + t549 * t608) * r_i_i_C(1) + (t537 * t608 + t549 * t612) * r_i_i_C(2)) * qJD(6); 0, (t517 * t612 - t530 * t608) * r_i_i_C(1) + (-t517 * t608 - t530 * t612) * r_i_i_C(2) + t517 * pkin(5) + t531 * pkin(4) - t530 * pkin(10) - t579 * t600 - t584 * t645 + t578 * t582 - t585 * t590 + t664 * (t543 * qJD(5) + t531 * t609 + t578 * t656) + ((-t543 * t608 + t556 * t612) * r_i_i_C(1) + (-t543 * t612 - t556 * t608) * r_i_i_C(2)) * qJD(6), (t526 * t608 + t620 * t646) * r_i_i_C(1) + (t526 * t612 - t620 * t647) * r_i_i_C(2) + t526 * pkin(10) + t616 * t525 + t615 * t546 + (-t579 * t641 + t578 * t610 + (-t663 * t585 + (-t584 * t606 + t603 * t655) * t610) * qJD(3)) * pkin(3), t579 * t603, t664 * t513 + t632 * t622 + t626 * (-t535 * qJD(5) - t526 * t609 + t579 * t656), (-t513 * t608 - t525 * t612) * r_i_i_C(1) + (-t513 * t612 + t525 * t608) * r_i_i_C(2) + ((-t535 * t612 + t546 * t608) * r_i_i_C(1) + (t535 * t608 + t546 * t612) * r_i_i_C(2)) * qJD(6); 0, (t523 * t612 + t541 * t608) * r_i_i_C(1) + (-t523 * t608 + t541 * t612) * r_i_i_C(2) + t523 * pkin(5) + t542 * pkin(4) + t541 * pkin(10) + t664 * (t560 * qJD(5) + t542 * t609 - t613 * t635) + ((-t560 * t608 + t561 * t612) * r_i_i_C(1) + (-t560 * t612 - t561 * t608) * r_i_i_C(2)) * qJD(6) + (-t614 * t645 - t611 * t590 + (-t582 * t614 - t600 * t611) * qJD(2)) * t604, (-t538 * t608 + t617 * t646) * r_i_i_C(1) + (-t538 * t612 - t617 * t647) * r_i_i_C(2) - t538 * pkin(10) + t616 * t539 + t615 * t554 + (-t607 * t603 * t648 + ((-t606 * t650 - t640) * qJD(3) + (-t606 * t640 - t650) * qJD(2)) * t604) * pkin(3), t634, t664 * t521 + t630 * t622 + t626 * (-t552 * qJD(5) - t540 * t609 + t613 * t634), (-t521 * t608 - t539 * t612) * r_i_i_C(1) + (-t521 * t612 + t539 * t608) * r_i_i_C(2) + ((-t552 * t612 + t554 * t608) * r_i_i_C(1) + (t552 * t608 + t554 * t612) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end