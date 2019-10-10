% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR4
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
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (90->45), mult. (273->84), div. (0->0), fcn. (260->8), ass. (0->34)
	t190 = cos(pkin(6));
	t187 = sin(pkin(11));
	t189 = cos(pkin(11));
	t191 = sin(qJ(2));
	t193 = cos(qJ(2));
	t196 = t193 * t187 + t191 * t189;
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
	t180 = t191 * t187 - t193 * t189;
	t199 = -t193 * pkin(2) + t180 * r_i_i_C(1) - pkin(1);
	t172 = t190 * t187 * t203 - t189 * t200;
	t178 = -t187 * t202 - t189 * t203;
	t198 = t194 * t172 - t192 * t178;
	t197 = t192 * t172 + t194 * t178;
	t195 = t175 * r_i_i_C(1) + t190 * t209 + (-r_i_i_C(3) - pkin(8) - qJ(3)) * t188;
	t179 = pkin(2) * t200 - t188 * qJD(3);
	t177 = t180 * qJD(2);
	t174 = t180 * t190;
	t173 = qJD(2) * t175;
	t171 = -t194 * t173 + t192 * t177 + (t174 * t192 - t194 * t196) * qJD(1);
	t170 = t192 * t173 + t194 * t177 + (t174 * t194 + t192 * t196) * qJD(1);
	t1 = [t198 * r_i_i_C(1) - t171 * r_i_i_C(2) + t192 * t201 - t194 * t179 + (t195 * t192 + t199 * t194) * qJD(1), t170 * r_i_i_C(1) + ((t175 * t194 - t180 * t192) * qJD(1) - t197) * r_i_i_C(2) + ((t190 * t208 - t205) * qJD(2) + (-t190 * t205 + t208) * qJD(1)) * pkin(2), t194 * t204, 0, 0, 0; t197 * r_i_i_C(1) + t170 * r_i_i_C(2) - t194 * t201 - t192 * t179 + (t199 * t192 - t195 * t194) * qJD(1), t171 * r_i_i_C(1) + ((t175 * t192 + t180 * t194) * qJD(1) + t198) * r_i_i_C(2) + ((-t190 * t207 - t206) * qJD(2) + (-t190 * t206 - t207) * qJD(1)) * pkin(2), t192 * t204, 0, 0, 0; 0, (-t196 * r_i_i_C(1) + t180 * r_i_i_C(2) - t209) * t188 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:43
	% EndTime: 2019-10-10 10:09:44
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (297->77), mult. (888->135), div. (0->0), fcn. (922->10), ass. (0->59)
	t316 = cos(pkin(6));
	t313 = sin(pkin(11));
	t315 = cos(pkin(11));
	t321 = cos(qJ(2));
	t340 = qJD(2) * t321;
	t318 = sin(qJ(2));
	t341 = qJD(2) * t318;
	t352 = t313 * t341 - t315 * t340;
	t294 = t352 * t316;
	t328 = t321 * t313 + t318 * t315;
	t300 = t328 * t316;
	t304 = t318 * t313 - t321 * t315;
	t322 = cos(qJ(1));
	t319 = sin(qJ(1));
	t342 = qJD(1) * t319;
	t302 = -t313 * t340 - t315 * t341;
	t345 = t319 * t302;
	t286 = t345 - t300 * t342 + (-qJD(1) * t304 - t294) * t322;
	t314 = sin(pkin(6));
	t348 = t314 * t322;
	t333 = qJD(4) * t348;
	t353 = t286 - t333;
	t351 = r_i_i_C(3) + pkin(9);
	t350 = pkin(2) * t316;
	t349 = t314 * t319;
	t347 = t318 * t319;
	t346 = t318 * t322;
	t344 = t319 * t321;
	t343 = t321 * t322;
	t288 = t322 * t300 - t319 * t304;
	t339 = qJD(4) * t288;
	t338 = pkin(2) * t341;
	t337 = t314 * t342;
	t336 = qJD(1) * t348;
	t317 = sin(qJ(4));
	t320 = cos(qJ(4));
	t332 = t317 * r_i_i_C(1) + t320 * r_i_i_C(2);
	t299 = t304 * t316;
	t331 = -t322 * t299 - t319 * t328;
	t330 = t319 * t299 - t322 * t328;
	t329 = t319 * t300 + t322 * t304;
	t327 = t320 * r_i_i_C(1) - t317 * r_i_i_C(2) + pkin(3);
	t326 = qJD(4) * t332;
	t325 = t337 - t339;
	t324 = qJD(2) * t328;
	t323 = t304 * qJD(2);
	t283 = -t288 * qJD(1) + t319 * t294 + t322 * t302;
	t312 = t321 * pkin(2) + pkin(1);
	t308 = t320 * t333;
	t303 = -t314 * qJD(3) + t340 * t350;
	t301 = t318 * t350 + (-pkin(8) - qJ(3)) * t314;
	t298 = t328 * t314;
	t295 = t316 * t324;
	t292 = t352 * t314;
	t285 = t330 * qJD(1) - t322 * t295 + t319 * t323;
	t282 = t331 * qJD(1) - t319 * t295 - t322 * t323;
	t280 = t317 * t336 + t283 * t320 + (t317 * t329 + t320 * t349) * qJD(4);
	t279 = t320 * t336 - t283 * t317 + (-t317 * t349 + t320 * t329) * qJD(4);
	t1 = [(-t286 * t320 + t317 * t339 + t308) * r_i_i_C(1) + (t353 * t317 + t320 * t339) * r_i_i_C(2) - t286 * pkin(3) + t319 * t338 - t322 * t303 + t351 * t285 + (-t322 * t312 + (-t332 * t314 + t301) * t319) * qJD(1), t351 * t283 - t330 * t326 - t327 * t282 + ((t316 * t347 - t343) * qJD(2) + (-t316 * t343 + t347) * qJD(1)) * pkin(2), t336, t279 * r_i_i_C(1) - t280 * r_i_i_C(2), 0, 0; -t322 * t338 + t283 * pkin(3) + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) - t319 * t303 + t351 * t282 + (-t301 * t322 - t312 * t319) * qJD(1), -t351 * (qJD(1) * t329 + t322 * t294 - t345) - t331 * t326 + t327 * t285 + ((-t316 * t346 - t344) * qJD(2) + (-t316 * t344 - t346) * qJD(1)) * pkin(2), t337, t308 * r_i_i_C(2) + (r_i_i_C(1) * t325 - t286 * r_i_i_C(2)) * t320 + (-t353 * r_i_i_C(1) - t325 * r_i_i_C(2)) * t317, 0, 0; 0, -t351 * t292 + (t304 * t326 - t324 * t327 - t338) * t314, 0, t332 * t292 + ((-t298 * t320 - t316 * t317) * r_i_i_C(1) + (t298 * t317 - t316 * t320) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:43
	% EndTime: 2019-10-10 10:09:44
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (445->104), mult. (1155->173), div. (0->0), fcn. (1210->12), ass. (0->64)
	t331 = cos(pkin(6));
	t328 = sin(pkin(11));
	t330 = cos(pkin(11));
	t334 = sin(qJ(2));
	t337 = cos(qJ(2));
	t344 = t328 * t337 + t330 * t334;
	t311 = t344 * t331;
	t354 = qJD(2) * t337;
	t348 = t331 * t354;
	t355 = qJD(2) * t334;
	t305 = t328 * t331 * t355 - t330 * t348;
	t315 = t328 * t334 - t337 * t330;
	t338 = cos(qJ(1));
	t335 = sin(qJ(1));
	t356 = qJD(1) * t335;
	t313 = -t328 * t354 - t330 * t355;
	t359 = t335 * t313;
	t296 = t359 - t311 * t356 + (-qJD(1) * t315 - t305) * t338;
	t329 = sin(pkin(6));
	t362 = t329 * t338;
	t347 = qJD(4) * t362;
	t365 = t296 - t347;
	t364 = r_i_i_C(3) + qJ(5) + pkin(9);
	t363 = t329 * t335;
	t361 = t334 * t335;
	t360 = t334 * t338;
	t358 = t335 * t337;
	t357 = t337 * t338;
	t298 = t311 * t338 - t335 * t315;
	t353 = qJD(4) * t298;
	t352 = pkin(2) * t355;
	t351 = t329 * t356;
	t350 = qJD(1) * t362;
	t349 = t329 * t355;
	t327 = qJ(4) + pkin(12);
	t325 = sin(t327);
	t326 = cos(t327);
	t346 = -r_i_i_C(1) * t325 - r_i_i_C(2) * t326;
	t310 = t315 * t331;
	t297 = -t310 * t338 - t335 * t344;
	t345 = t335 * t311 + t315 * t338;
	t336 = cos(qJ(4));
	t323 = pkin(4) * t336 + pkin(3);
	t343 = r_i_i_C(1) * t326 - r_i_i_C(2) * t325 + t323;
	t342 = t351 - t353;
	t309 = t344 * t329;
	t341 = t315 * qJD(2);
	t333 = sin(qJ(4));
	t340 = pkin(4) * t333 - t346;
	t339 = qJD(4) * t340;
	t293 = -t298 * qJD(1) + t335 * t305 + t313 * t338;
	t324 = pkin(2) * t337 + pkin(1);
	t317 = t326 * t347;
	t314 = pkin(2) * t348 - qJD(3) * t329;
	t312 = pkin(2) * t331 * t334 + (-pkin(8) - qJ(3)) * t329;
	t306 = qJD(2) * t311;
	t304 = qJD(2) * t309;
	t303 = -t329 * t330 * t354 + t328 * t349;
	t299 = t335 * t310 - t338 * t344;
	t295 = t310 * t356 + (-qJD(1) * t344 - t306) * t338 + t335 * t341;
	t292 = t297 * qJD(1) - t335 * t306 - t338 * t341;
	t290 = t325 * t350 + t293 * t326 + (t325 * t345 + t326 * t363) * qJD(4);
	t289 = t326 * t350 - t293 * t325 + (-t325 * t363 + t326 * t345) * qJD(4);
	t1 = [(-t296 * t326 + t325 * t353 + t317) * r_i_i_C(1) + (t365 * t325 + t326 * t353) * r_i_i_C(2) - t296 * t323 + t297 * qJD(5) + t335 * t352 - t338 * t314 + t364 * t295 + (-t338 * t324 + (t346 * t329 + t312) * t335) * qJD(1) + (-t333 * t351 + (t298 * t333 + t336 * t362) * qJD(4)) * pkin(4), -t345 * qJD(5) + t364 * t293 - t343 * t292 - t299 * t339 + ((t331 * t361 - t357) * qJD(2) + (-t331 * t357 + t361) * qJD(1)) * pkin(2), t350, t289 * r_i_i_C(1) - t290 * r_i_i_C(2) + (t336 * t350 - t293 * t333 + (-t333 * t363 + t336 * t345) * qJD(4)) * pkin(4), t292, 0; -t338 * t352 + t290 * r_i_i_C(1) + t289 * r_i_i_C(2) - t299 * qJD(5) + t293 * t323 - t335 * t314 + t364 * t292 + (-t312 * t338 - t324 * t335) * qJD(1) + (t333 * t350 + (t333 * t345 + t336 * t363) * qJD(4)) * pkin(4), t298 * qJD(5) - t364 * (t345 * qJD(1) + t305 * t338 - t359) + t343 * t295 - t297 * t339 + ((-t331 * t360 - t358) * qJD(2) + (-t331 * t358 - t360) * qJD(1)) * pkin(2), t351, t317 * r_i_i_C(2) + (t342 * r_i_i_C(1) - t296 * r_i_i_C(2)) * t326 + (-t365 * r_i_i_C(1) - t342 * r_i_i_C(2)) * t325 + (t336 * t351 - t296 * t333 + (-t298 * t336 + t333 * t362) * qJD(4)) * pkin(4), -t295, 0; 0, t315 * t329 * t339 - pkin(2) * t349 + qJD(5) * t309 - t364 * t303 - t343 * t304, 0, t340 * t303 + ((-t309 * t326 - t325 * t331) * r_i_i_C(1) + (t309 * t325 - t326 * t331) * r_i_i_C(2) + (-t309 * t336 - t331 * t333) * pkin(4)) * qJD(4), t304, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:45
	% EndTime: 2019-10-10 10:09:46
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (1128->148), mult. (2753->244), div. (0->0), fcn. (3018->14), ass. (0->82)
	t469 = cos(pkin(6));
	t467 = sin(pkin(11));
	t521 = cos(pkin(11));
	t524 = cos(qJ(2));
	t495 = t524 * t521;
	t473 = sin(qJ(2));
	t509 = qJD(2) * t473;
	t528 = -qJD(2) * t495 + t467 * t509;
	t443 = t528 * t469;
	t498 = t473 * t521;
	t486 = t524 * t467 + t498;
	t448 = t486 * t469;
	t499 = t524 * qJD(2);
	t450 = -qJD(2) * t498 - t467 * t499;
	t474 = sin(qJ(1));
	t477 = cos(qJ(1));
	t485 = -t473 * t467 + t495;
	t511 = qJD(1) * t474;
	t418 = t448 * t511 - t474 * t450 + (-qJD(1) * t485 + t443) * t477;
	t430 = t477 * t448 + t474 * t485;
	t466 = qJ(4) + pkin(12);
	t464 = sin(t466);
	t465 = cos(t466);
	t468 = sin(pkin(6));
	t502 = t468 * t511;
	t515 = t468 * t477;
	t505 = t465 * t515;
	t412 = (-qJD(4) * t430 + t502) * t464 - qJD(4) * t505 - t418 * t465;
	t444 = qJD(2) * t448;
	t510 = qJD(1) * t477;
	t484 = t485 * t469;
	t529 = qJD(1) * t484 + t485 * qJD(2);
	t419 = -t477 * t444 - t474 * t529 - t486 * t510;
	t471 = sin(qJ(6));
	t475 = cos(qJ(6));
	t535 = t412 * t471 + t419 * t475;
	t534 = -t412 * t475 + t419 * t471;
	t425 = -t430 * t465 + t464 * t515;
	t429 = -t474 * t486 + t477 * t484;
	t533 = -t425 * t471 + t429 * t475;
	t532 = t425 * t475 + t429 * t471;
	t472 = sin(qJ(4));
	t489 = qJD(6) * (t471 * r_i_i_C(1) + t475 * r_i_i_C(2));
	t492 = t475 * r_i_i_C(1) - t471 * r_i_i_C(2) + pkin(5);
	t525 = pkin(10) + r_i_i_C(3);
	t531 = qJD(4) * (t472 * pkin(4) + t492 * t464 - t525 * t465) + t465 * t489;
	t476 = cos(qJ(4));
	t462 = t476 * pkin(4) + pkin(3);
	t527 = t525 * t464 + t492 * t465 + t462;
	t523 = pkin(2) * t469;
	t516 = t468 * t474;
	t513 = t473 * t474;
	t512 = t473 * t477;
	t508 = qJD(6) * t471;
	t507 = qJD(6) * t475;
	t506 = pkin(2) * t509;
	t504 = t524 * t474;
	t503 = t524 * t477;
	t501 = t468 * t510;
	t447 = t486 * t468;
	t435 = t447 * t465 + t469 * t464;
	t493 = -t447 * t464 + t469 * t465;
	t431 = t474 * t448 - t477 * t485;
	t490 = t431 * t464 + t465 * t516;
	t427 = -t431 * t465 + t464 * t516;
	t480 = -t430 * qJD(1) + t474 * t443 + t477 * t450;
	t479 = t425 * qJD(4) + t418 * t464 + t465 * t502;
	t470 = -qJ(5) - pkin(9);
	t463 = t524 * pkin(2) + pkin(1);
	t451 = -t468 * qJD(3) + t499 * t523;
	t449 = t473 * t523 + (-pkin(8) - qJ(3)) * t468;
	t446 = t485 * t468;
	t442 = qJD(2) * t447;
	t441 = t528 * t468;
	t432 = -t474 * t484 - t477 * t486;
	t422 = t493 * qJD(4) - t441 * t465;
	t416 = -t474 * t444 + t477 * t529 - t486 * t511;
	t410 = t490 * qJD(4) + t464 * t501 + t465 * t480;
	t409 = t427 * qJD(4) + t464 * t480 - t465 * t501;
	t408 = t410 * t475 + t416 * t471 + (-t427 * t471 - t432 * t475) * qJD(6);
	t407 = -t410 * t471 + t416 * t475 + (-t427 * t475 + t432 * t471) * qJD(6);
	t1 = [t534 * r_i_i_C(1) + t535 * r_i_i_C(2) - t412 * pkin(5) + t418 * t462 - t419 * t470 + t429 * qJD(5) + t474 * t506 - t477 * t451 + t525 * t479 + (t533 * r_i_i_C(1) - t532 * r_i_i_C(2)) * qJD(6) + (t474 * t449 - t477 * t463) * qJD(1) + (-t472 * t502 + (t430 * t472 + t476 * t515) * qJD(4)) * pkin(4), (-t431 * t507 + t471 * t480) * r_i_i_C(1) + (t431 * t508 + t475 * t480) * r_i_i_C(2) - t480 * t470 - t431 * qJD(5) - t527 * t416 + ((t469 * t513 - t503) * qJD(2) + (-t469 * t503 + t513) * qJD(1)) * pkin(2) - t531 * t432, t501, t525 * t410 - t490 * t489 - t492 * t409 + (t476 * t501 - t480 * t472 + (t431 * t476 - t472 * t516) * qJD(4)) * pkin(4), t416, t407 * r_i_i_C(1) - t408 * r_i_i_C(2); -t477 * t506 + t410 * pkin(5) + t408 * r_i_i_C(1) + t407 * r_i_i_C(2) - t432 * qJD(5) - t416 * t470 + t480 * t462 - t474 * t451 + t525 * t409 + (-t449 * t477 - t463 * t474) * qJD(1) + (t472 * t501 + (t431 * t472 + t476 * t516) * qJD(4)) * pkin(4), (-t418 * t471 + t430 * t507) * r_i_i_C(1) + (-t418 * t475 - t430 * t508) * r_i_i_C(2) + t418 * t470 + t430 * qJD(5) + t527 * t419 + ((-t469 * t512 - t504) * qJD(2) + (-t469 * t504 - t512) * qJD(1)) * pkin(2) - t531 * t429, t502, t525 * t412 - (-t430 * t464 - t505) * t489 + t492 * t479 + (t476 * t502 + t418 * t472 + (-t430 * t476 + t472 * t515) * qJD(4)) * pkin(4), -t419, -t535 * r_i_i_C(1) + t534 * r_i_i_C(2) + (t532 * r_i_i_C(1) + t533 * r_i_i_C(2)) * qJD(6); 0, (-t441 * t471 + t447 * t507) * r_i_i_C(1) + (-t441 * t475 - t447 * t508) * r_i_i_C(2) + t441 * t470 + t447 * qJD(5) - t468 * t506 - t527 * t442 - t531 * t446, 0, t525 * t422 - t493 * t489 + t492 * (-t435 * qJD(4) + t441 * t464) + (t441 * t472 + (-t447 * t476 - t469 * t472) * qJD(4)) * pkin(4), t442, (-t422 * t471 + t442 * t475) * r_i_i_C(1) + (-t422 * t475 - t442 * t471) * r_i_i_C(2) + ((-t435 * t475 + t446 * t471) * r_i_i_C(1) + (t435 * t471 + t446 * t475) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end