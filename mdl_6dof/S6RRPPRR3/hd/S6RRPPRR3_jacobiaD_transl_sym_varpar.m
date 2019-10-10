% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
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
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
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
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (90->45), mult. (273->84), div. (0->0), fcn. (260->8), ass. (0->34)
	t190 = cos(pkin(6));
	t187 = sin(pkin(11));
	t189 = cos(pkin(11));
	t191 = sin(qJ(2));
	t193 = cos(qJ(2));
	t196 = t187 * t193 + t189 * t191;
	t175 = t196 * t190;
	t209 = t191 * pkin(2);
	t192 = sin(qJ(1));
	t208 = t191 * t192;
	t194 = cos(qJ(1));
	t207 = t191 * t194;
	t206 = t192 * t193;
	t205 = t193 * t194;
	t188 = sin(pkin(6));
	t204 = qJD(1) * t188;
	t203 = qJD(2) * t191;
	t202 = qJD(2) * t193;
	t201 = pkin(2) * t203;
	t200 = t190 * t202;
	t180 = t187 * t191 - t193 * t189;
	t199 = -pkin(2) * t193 + r_i_i_C(1) * t180 - pkin(1);
	t172 = t187 * t190 * t203 - t189 * t200;
	t178 = -t187 * t202 - t189 * t203;
	t198 = t194 * t172 - t192 * t178;
	t197 = t192 * t172 + t194 * t178;
	t195 = t175 * r_i_i_C(1) + t190 * t209 + (-r_i_i_C(3) - pkin(8) - qJ(3)) * t188;
	t179 = pkin(2) * t200 - qJD(3) * t188;
	t177 = t180 * qJD(2);
	t174 = t180 * t190;
	t173 = qJD(2) * t175;
	t171 = -t194 * t173 + t192 * t177 + (t174 * t192 - t194 * t196) * qJD(1);
	t170 = t192 * t173 + t194 * t177 + (t174 * t194 + t192 * t196) * qJD(1);
	t1 = [t198 * r_i_i_C(1) - t171 * r_i_i_C(2) + t192 * t201 - t194 * t179 + (t192 * t195 + t194 * t199) * qJD(1), t170 * r_i_i_C(1) + ((t175 * t194 - t180 * t192) * qJD(1) - t197) * r_i_i_C(2) + ((t190 * t208 - t205) * qJD(2) + (-t190 * t205 + t208) * qJD(1)) * pkin(2), t194 * t204, 0, 0, 0; t197 * r_i_i_C(1) + t170 * r_i_i_C(2) - t194 * t201 - t192 * t179 + (t192 * t199 - t194 * t195) * qJD(1), t171 * r_i_i_C(1) + ((t175 * t192 + t180 * t194) * qJD(1) + t198) * r_i_i_C(2) + ((-t190 * t207 - t206) * qJD(2) + (-t190 * t206 - t207) * qJD(1)) * pkin(2), t192 * t204, 0, 0, 0; 0, (-r_i_i_C(1) * t196 + r_i_i_C(2) * t180 - t209) * t188 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:20
	% EndTime: 2019-10-10 09:39:20
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (227->54), mult. (677->94), div. (0->0), fcn. (692->10), ass. (0->42)
	t269 = sin(pkin(11));
	t272 = cos(pkin(11));
	t274 = sin(qJ(2));
	t276 = cos(qJ(2));
	t282 = t276 * t269 + t274 * t272;
	t298 = qJD(2) * t282;
	t268 = sin(pkin(12));
	t270 = sin(pkin(6));
	t271 = cos(pkin(12));
	t273 = cos(pkin(6));
	t297 = -t273 * t274 * pkin(2) + (r_i_i_C(1) * t268 + r_i_i_C(2) * t271 + pkin(8) + qJ(3)) * t270;
	t296 = r_i_i_C(3) + qJ(4);
	t275 = sin(qJ(1));
	t295 = t274 * t275;
	t277 = cos(qJ(1));
	t294 = t274 * t277;
	t293 = t275 * t276;
	t292 = t276 * t277;
	t291 = qJD(1) * t275;
	t290 = qJD(2) * t274;
	t289 = qJD(2) * t276;
	t288 = pkin(2) * t290;
	t287 = t273 * t289;
	t261 = t274 * t269 - t276 * t272;
	t256 = t261 * t273;
	t285 = -t277 * t256 - t275 * t282;
	t257 = t282 * t273;
	t284 = -t277 * t257 + t275 * t261;
	t283 = t275 * t257 + t277 * t261;
	t281 = t271 * r_i_i_C(1) - t268 * r_i_i_C(2) + pkin(3);
	t280 = t261 * qJD(2);
	t253 = t273 * t269 * t290 - t272 * t287;
	t259 = -t269 * t289 - t272 * t290;
	t279 = t284 * qJD(1) + t275 * t253 + t277 * t259;
	t278 = t283 * qJD(1) + t277 * t253 - t275 * t259;
	t267 = t276 * pkin(2) + pkin(1);
	t260 = pkin(2) * t287 - t270 * qJD(3);
	t254 = t273 * t298;
	t252 = t270 * t298;
	t248 = t256 * t291 + (-qJD(1) * t282 - t254) * t277 + t275 * t280;
	t245 = t285 * qJD(1) - t275 * t254 - t277 * t280;
	t1 = [t285 * qJD(4) + t275 * t288 - t277 * t260 + t296 * t248 + t281 * t278 + (-t277 * t267 - t297 * t275) * qJD(1), -t283 * qJD(4) + t296 * t279 - t281 * t245 + ((t273 * t295 - t292) * qJD(2) + (-t273 * t292 + t295) * qJD(1)) * pkin(2), qJD(1) * t277 * t270, t245, 0, 0; -(t275 * t256 - t277 * t282) * qJD(4) - t277 * t288 - t275 * t260 + t296 * t245 + t281 * t279 + (-t275 * t267 + t297 * t277) * qJD(1), -t284 * qJD(4) - t296 * t278 + t281 * t248 + ((-t273 * t294 - t293) * qJD(2) + (-t273 * t293 - t294) * qJD(1)) * pkin(2), t270 * t291, -t248, 0, 0; 0, -t281 * t252 + (t282 * qJD(4) - t296 * t280 - t288) * t270, 0, t252, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:20
	% EndTime: 2019-10-10 09:39:21
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (392->89), mult. (992->146), div. (0->0), fcn. (1041->12), ass. (0->62)
	t332 = cos(pkin(6));
	t329 = sin(pkin(11));
	t331 = cos(pkin(11));
	t334 = sin(qJ(2));
	t336 = cos(qJ(2));
	t342 = t336 * t329 + t334 * t331;
	t311 = t342 * t332;
	t353 = qJD(2) * t334;
	t346 = t329 * t353;
	t352 = qJD(2) * t336;
	t347 = t332 * t352;
	t305 = -t331 * t347 + t332 * t346;
	t315 = t334 * t329 - t336 * t331;
	t337 = cos(qJ(1));
	t335 = sin(qJ(1));
	t354 = qJD(1) * t335;
	t313 = -t329 * t352 - t331 * t353;
	t357 = t335 * t313;
	t296 = t357 - t311 * t354 + (-qJD(1) * t315 - t305) * t337;
	t330 = sin(pkin(6));
	t360 = t330 * t337;
	t345 = qJD(5) * t360;
	t364 = t296 - t345;
	t363 = pkin(4) * sin(pkin(12));
	t362 = r_i_i_C(3) + pkin(9) + qJ(4);
	t361 = t330 * t335;
	t359 = t334 * t335;
	t358 = t334 * t337;
	t356 = t335 * t336;
	t355 = t336 * t337;
	t298 = t337 * t311 - t335 * t315;
	t351 = qJD(5) * t298;
	t350 = pkin(2) * t353;
	t349 = t330 * t354;
	t348 = qJD(1) * t360;
	t327 = pkin(12) + qJ(5);
	t325 = sin(t327);
	t326 = cos(t327);
	t344 = t325 * r_i_i_C(1) + t326 * r_i_i_C(2);
	t310 = t315 * t332;
	t297 = -t337 * t310 - t335 * t342;
	t343 = t335 * t311 + t337 * t315;
	t323 = cos(pkin(12)) * pkin(4) + pkin(3);
	t341 = t326 * r_i_i_C(1) - t325 * r_i_i_C(2) + t323;
	t340 = qJD(5) * t344;
	t339 = t349 - t351;
	t309 = t342 * t330;
	t338 = t315 * qJD(2);
	t293 = -t298 * qJD(1) + t335 * t305 + t337 * t313;
	t324 = t336 * pkin(2) + pkin(1);
	t317 = t326 * t345;
	t314 = pkin(2) * t347 - t330 * qJD(3);
	t312 = t332 * t334 * pkin(2) + (-pkin(8) - qJ(3)) * t330;
	t306 = qJD(2) * t311;
	t304 = qJD(2) * t309;
	t303 = (-t331 * t352 + t346) * t330;
	t299 = t335 * t310 - t337 * t342;
	t295 = t310 * t354 + (-qJD(1) * t342 - t306) * t337 + t335 * t338;
	t292 = t297 * qJD(1) - t335 * t306 - t337 * t338;
	t290 = t325 * t348 + t293 * t326 + (t325 * t343 + t326 * t361) * qJD(5);
	t289 = t326 * t348 - t293 * t325 + (-t325 * t361 + t326 * t343) * qJD(5);
	t1 = [(-t296 * t326 + t325 * t351 + t317) * r_i_i_C(1) + (t364 * t325 + t326 * t351) * r_i_i_C(2) - t296 * t323 + t297 * qJD(4) + t335 * t350 - t337 * t314 + t362 * t295 + (-t337 * t324 + (t312 + (-t344 - t363) * t330) * t335) * qJD(1), -t343 * qJD(4) + t362 * t293 - t299 * t340 - t341 * t292 + ((t332 * t359 - t355) * qJD(2) + (-t332 * t355 + t359) * qJD(1)) * pkin(2), t348, t292, t289 * r_i_i_C(1) - t290 * r_i_i_C(2), 0; -t337 * t350 + t290 * r_i_i_C(1) + t289 * r_i_i_C(2) - t299 * qJD(4) + t293 * t323 - t335 * t314 + t362 * t292 + (-t324 * t335 + (t330 * t363 - t312) * t337) * qJD(1), t298 * qJD(4) - t362 * (t343 * qJD(1) + t337 * t305 - t357) - t297 * t340 + t341 * t295 + ((-t332 * t358 - t356) * qJD(2) + (-t332 * t356 - t358) * qJD(1)) * pkin(2), t349, -t295, t317 * r_i_i_C(2) + (t339 * r_i_i_C(1) - t296 * r_i_i_C(2)) * t326 + (-t364 * r_i_i_C(1) - t339 * r_i_i_C(2)) * t325, 0; 0, t309 * qJD(4) - t362 * t303 - t341 * t304 + (t315 * t340 - t350) * t330, 0, t304, t344 * t303 + ((-t309 * t326 - t325 * t332) * r_i_i_C(1) + (t309 * t325 - t326 * t332) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:22
	% EndTime: 2019-10-10 09:39:23
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (1075->131), mult. (2590->217), div. (0->0), fcn. (2849->14), ass. (0->80)
	t469 = cos(pkin(6));
	t467 = sin(pkin(11));
	t520 = cos(pkin(11));
	t522 = cos(qJ(2));
	t493 = t522 * t520;
	t472 = sin(qJ(2));
	t508 = qJD(2) * t472;
	t526 = -qJD(2) * t493 + t467 * t508;
	t442 = t526 * t469;
	t497 = t472 * t520;
	t484 = t522 * t467 + t497;
	t447 = t484 * t469;
	t498 = qJD(2) * t522;
	t449 = -qJD(2) * t497 - t467 * t498;
	t473 = sin(qJ(1));
	t475 = cos(qJ(1));
	t483 = -t472 * t467 + t493;
	t510 = qJD(1) * t473;
	t417 = t447 * t510 - t473 * t449 + (-qJD(1) * t483 + t442) * t475;
	t429 = t447 * t475 + t473 * t483;
	t465 = pkin(12) + qJ(5);
	t463 = sin(t465);
	t464 = cos(t465);
	t468 = sin(pkin(6));
	t501 = t468 * t510;
	t514 = t468 * t475;
	t504 = t464 * t514;
	t411 = (-qJD(5) * t429 + t501) * t463 - qJD(5) * t504 - t417 * t464;
	t443 = qJD(2) * t447;
	t509 = qJD(1) * t475;
	t482 = t469 * t483;
	t527 = qJD(1) * t482 + qJD(2) * t483;
	t418 = -t475 * t443 - t527 * t473 - t484 * t509;
	t471 = sin(qJ(6));
	t474 = cos(qJ(6));
	t533 = t411 * t471 + t418 * t474;
	t532 = -t411 * t474 + t418 * t471;
	t424 = -t429 * t464 + t463 * t514;
	t428 = -t473 * t484 + t475 * t482;
	t531 = -t424 * t471 + t428 * t474;
	t530 = t424 * t474 + t428 * t471;
	t487 = qJD(6) * (r_i_i_C(1) * t471 + r_i_i_C(2) * t474);
	t490 = r_i_i_C(1) * t474 - r_i_i_C(2) * t471 + pkin(5);
	t523 = r_i_i_C(3) + pkin(10);
	t529 = (t490 * t463 - t523 * t464) * qJD(5) + t464 * t487;
	t461 = cos(pkin(12)) * pkin(4) + pkin(3);
	t524 = t523 * t463 + t490 * t464 + t461;
	t521 = pkin(2) * t469;
	t515 = t468 * t473;
	t512 = t472 * t473;
	t511 = t472 * t475;
	t507 = qJD(6) * t471;
	t506 = qJD(6) * t474;
	t505 = pkin(2) * t508;
	t503 = t522 * t473;
	t502 = t522 * t475;
	t500 = t468 * t509;
	t494 = -t472 * t521 + (pkin(4) * sin(pkin(12)) + pkin(8) + qJ(3)) * t468;
	t446 = t484 * t468;
	t434 = t446 * t464 + t463 * t469;
	t491 = -t446 * t463 + t464 * t469;
	t430 = t473 * t447 - t475 * t483;
	t488 = t430 * t463 + t464 * t515;
	t426 = -t430 * t464 + t463 * t515;
	t478 = -t429 * qJD(1) + t473 * t442 + t475 * t449;
	t477 = t424 * qJD(5) + t417 * t463 + t464 * t501;
	t470 = -pkin(9) - qJ(4);
	t462 = t522 * pkin(2) + pkin(1);
	t450 = -t468 * qJD(3) + t498 * t521;
	t445 = t483 * t468;
	t441 = qJD(2) * t446;
	t440 = t526 * t468;
	t431 = -t473 * t482 - t475 * t484;
	t421 = t491 * qJD(5) - t440 * t464;
	t415 = -t473 * t443 + t527 * t475 - t484 * t510;
	t409 = t488 * qJD(5) + t463 * t500 + t464 * t478;
	t408 = t426 * qJD(5) + t463 * t478 - t464 * t500;
	t407 = t409 * t474 + t415 * t471 + (-t426 * t471 - t431 * t474) * qJD(6);
	t406 = -t409 * t471 + t415 * t474 + (-t426 * t474 + t431 * t471) * qJD(6);
	t1 = [t532 * r_i_i_C(1) + t533 * r_i_i_C(2) - t411 * pkin(5) + t417 * t461 - t418 * t470 + t428 * qJD(4) + t473 * t505 - t475 * t450 + t523 * t477 + (t531 * r_i_i_C(1) - t530 * r_i_i_C(2)) * qJD(6) + (-t475 * t462 - t494 * t473) * qJD(1), (-t430 * t506 + t471 * t478) * r_i_i_C(1) + (t430 * t507 + t474 * t478) * r_i_i_C(2) - t478 * t470 - t430 * qJD(4) - t524 * t415 + ((t469 * t512 - t502) * qJD(2) + (-t469 * t502 + t512) * qJD(1)) * pkin(2) - t529 * t431, t500, t415, -t490 * t408 + t523 * t409 - t488 * t487, r_i_i_C(1) * t406 - r_i_i_C(2) * t407; -t475 * t505 + t409 * pkin(5) + t407 * r_i_i_C(1) + t406 * r_i_i_C(2) - t431 * qJD(4) - t415 * t470 + t478 * t461 - t473 * t450 + t523 * t408 + (-t462 * t473 + t494 * t475) * qJD(1), (-t417 * t471 + t429 * t506) * r_i_i_C(1) + (-t417 * t474 - t429 * t507) * r_i_i_C(2) + t417 * t470 + t429 * qJD(4) + t524 * t418 + ((-t469 * t511 - t503) * qJD(2) + (-t469 * t503 - t511) * qJD(1)) * pkin(2) - t529 * t428, t501, -t418, t523 * t411 - (-t429 * t463 - t504) * t487 + t490 * t477, -t533 * r_i_i_C(1) + t532 * r_i_i_C(2) + (t530 * r_i_i_C(1) + t531 * r_i_i_C(2)) * qJD(6); 0, (-t440 * t471 + t446 * t506) * r_i_i_C(1) + (-t440 * t474 - t446 * t507) * r_i_i_C(2) + t440 * t470 + t446 * qJD(4) - t468 * t505 - t524 * t441 - t529 * t445, 0, t441, t523 * t421 - t491 * t487 + t490 * (-t434 * qJD(5) + t440 * t463), (-t421 * t471 + t441 * t474) * r_i_i_C(1) + (-t421 * t474 - t441 * t471) * r_i_i_C(2) + ((-t434 * t474 + t445 * t471) * r_i_i_C(1) + (t434 * t471 + t445 * t474) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end