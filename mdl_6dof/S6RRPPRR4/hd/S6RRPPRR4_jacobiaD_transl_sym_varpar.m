% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:14
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
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:14
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
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (191->52), mult. (562->91), div. (0->0), fcn. (572->8), ass. (0->40)
	t238 = cos(pkin(6));
	t235 = sin(pkin(11));
	t237 = cos(pkin(11));
	t239 = sin(qJ(2));
	t241 = cos(qJ(2));
	t244 = t241 * t235 + t239 * t237;
	t224 = t244 * t238;
	t253 = qJD(2) * t241;
	t251 = t238 * t253;
	t254 = qJD(2) * t239;
	t220 = t238 * t235 * t254 - t237 * t251;
	t226 = -t235 * t253 - t237 * t254;
	t240 = sin(qJ(1));
	t242 = cos(qJ(1));
	t228 = t239 * t235 - t241 * t237;
	t245 = t240 * t224 + t242 * t228;
	t264 = t245 * qJD(1) + t242 * t220 - t240 * t226;
	t246 = t242 * t224 - t240 * t228;
	t263 = t246 * qJD(1) - t240 * t220 - t242 * t226;
	t262 = pkin(3) - r_i_i_C(2);
	t261 = r_i_i_C(3) + qJ(4);
	t260 = t239 * t240;
	t259 = t239 * t242;
	t258 = t240 * t241;
	t257 = t241 * t242;
	t256 = qJD(1) * t240;
	t236 = sin(pkin(6));
	t255 = qJD(2) * t236;
	t252 = pkin(2) * t254;
	t250 = -t238 * t239 * pkin(2) + (r_i_i_C(1) + pkin(8) + qJ(3)) * t236;
	t223 = t228 * t238;
	t247 = -t242 * t223 - t240 * t244;
	t243 = t228 * qJD(2);
	t234 = t241 * pkin(2) + pkin(1);
	t227 = pkin(2) * t251 - t236 * qJD(3);
	t221 = qJD(2) * t224;
	t218 = t244 * t255;
	t215 = t223 * t256 + (-qJD(1) * t244 - t221) * t242 + t240 * t243;
	t212 = t247 * qJD(1) - t240 * t221 - t242 * t243;
	t1 = [t247 * qJD(4) + t240 * t252 - t242 * t227 + t262 * t264 + t261 * t215 + (-t242 * t234 - t250 * t240) * qJD(1), -t245 * qJD(4) - t262 * t212 - t261 * t263 + ((t238 * t260 - t257) * qJD(2) + (-t238 * t257 + t260) * qJD(1)) * pkin(2), qJD(1) * t242 * t236, t212, 0, 0; -(t240 * t223 - t242 * t244) * qJD(4) - t242 * t252 - t240 * t227 - t262 * t263 + t261 * t212 + (-t240 * t234 + t250 * t242) * qJD(1), t246 * qJD(4) + t262 * t215 - t261 * t264 + ((-t238 * t259 - t258) * qJD(2) + (-t238 * t258 - t259) * qJD(1)) * pkin(2), t236 * t256, -t215, 0, 0; 0, -t261 * t228 * t255 - t262 * t218 + (t244 * qJD(4) - t252) * t236, 0, t218, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:15
	% EndTime: 2019-10-10 09:41:15
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (366->80), mult. (1084->138), div. (0->0), fcn. (1136->10), ass. (0->54)
	t317 = cos(pkin(6));
	t314 = sin(pkin(11));
	t316 = cos(pkin(11));
	t319 = sin(qJ(2));
	t322 = cos(qJ(2));
	t333 = t322 * t314 + t319 * t316;
	t303 = t333 * t317;
	t307 = t319 * t314 - t322 * t316;
	t326 = t307 * qJD(2);
	t343 = qJD(2) * t322;
	t338 = t317 * t343;
	t344 = qJD(2) * t319;
	t297 = t317 * t314 * t344 - t316 * t338;
	t305 = -t314 * t343 - t316 * t344;
	t320 = sin(qJ(1));
	t323 = cos(qJ(1));
	t345 = qJD(1) * t320;
	t354 = t303 * t345 - t320 * t305 + (qJD(1) * t307 + t297) * t323;
	t334 = t323 * t303 - t320 * t307;
	t353 = t334 * qJD(1) - t320 * t297 - t323 * t305;
	t318 = sin(qJ(5));
	t321 = cos(qJ(5));
	t336 = t321 * r_i_i_C(1) - t318 * r_i_i_C(2);
	t352 = t336 * qJD(5) + qJD(4);
	t315 = sin(pkin(6));
	t351 = t315 * t320;
	t350 = t315 * t323;
	t349 = t319 * t320;
	t348 = t319 * t323;
	t347 = t320 * t322;
	t346 = t322 * t323;
	t342 = -r_i_i_C(3) - pkin(9) - pkin(3);
	t341 = pkin(2) * t344;
	t340 = t315 * t345;
	t339 = qJD(1) * t350;
	t332 = t318 * r_i_i_C(1) + t321 * r_i_i_C(2) + qJ(4);
	t327 = t307 * t317;
	t287 = -t320 * t333 - t323 * t327;
	t331 = t287 * t318 + t321 * t350;
	t330 = t287 * t321 - t318 * t350;
	t328 = t333 * t315;
	t325 = t320 * t327;
	t313 = t322 * pkin(2) + pkin(1);
	t306 = pkin(2) * t338 - t315 * qJD(3);
	t304 = t317 * t319 * pkin(2) + (-pkin(8) - qJ(3)) * t315;
	t301 = t307 * t315;
	t298 = qJD(2) * t303;
	t295 = qJD(2) * t328;
	t289 = -t323 * t333 + t325;
	t284 = qJD(1) * t325 + (-qJD(1) * t333 - t298) * t323 + t320 * t326;
	t281 = t287 * qJD(1) - t320 * t298 - t323 * t326;
	t279 = t321 * t339 + t281 * t318 + (-t289 * t321 - t318 * t351) * qJD(5);
	t278 = -t318 * t339 + t281 * t321 + (t289 * t318 - t321 * t351) * qJD(5);
	t1 = [t320 * t341 + t287 * qJD(4) - t323 * t306 + t332 * t284 + (t330 * r_i_i_C(1) - t331 * r_i_i_C(2)) * qJD(5) - t342 * t354 + (-t323 * t313 + (t304 + (-pkin(4) - t336) * t315) * t320) * qJD(1), -t352 * (t320 * t303 + t323 * t307) + t342 * t281 - t332 * t353 + ((t317 * t349 - t346) * qJD(2) + (-t317 * t346 + t349) * qJD(1)) * pkin(2), t339, t281, t278 * r_i_i_C(1) - t279 * r_i_i_C(2), 0; -t323 * t341 + t279 * r_i_i_C(1) + t278 * r_i_i_C(2) + t281 * qJ(4) - t289 * qJD(4) - t320 * t306 + t342 * t353 + (-t313 * t320 + (pkin(4) * t315 - t304) * t323) * qJD(1), t352 * t334 - t342 * t284 - t332 * t354 + ((-t317 * t348 - t347) * qJD(2) + (-t317 * t347 - t348) * qJD(1)) * pkin(2), t340, -t284, (-t284 * t321 - t318 * t340) * r_i_i_C(1) + (t284 * t318 - t321 * t340) * r_i_i_C(2) + (t331 * r_i_i_C(1) + t330 * r_i_i_C(2)) * qJD(5), 0; 0, t352 * t328 + t342 * t295 + (-t332 * t326 - t341) * t315, 0, t295, t336 * t295 + ((-t301 * t318 - t317 * t321) * r_i_i_C(1) + (-t301 * t321 + t317 * t318) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:17
	% EndTime: 2019-10-10 09:41:18
	% DurationCPUTime: 1.18s
	% Computational Cost: add. (914->117), mult. (2682->194), div. (0->0), fcn. (2944->12), ass. (0->75)
	t474 = sin(qJ(2));
	t527 = sin(pkin(11));
	t529 = cos(pkin(6));
	t497 = t529 * t527;
	t528 = cos(pkin(11));
	t498 = t529 * t528;
	t531 = cos(qJ(2));
	t481 = -t474 * t497 + t531 * t498;
	t451 = t481 * qJD(2);
	t484 = t474 * t528 + t531 * t527;
	t458 = t484 * qJD(2);
	t475 = sin(qJ(1));
	t478 = cos(qJ(1));
	t460 = t474 * t527 - t531 * t528;
	t518 = -t474 * t498 - t531 * t497;
	t487 = t475 * t460 + t478 * t518;
	t543 = t487 * qJD(1) - t475 * t451 - t478 * t458;
	t538 = t475 * t518;
	t439 = -t478 * t460 + t538;
	t435 = -t475 * t484 + t478 * t481;
	t473 = sin(qJ(5));
	t477 = cos(qJ(5));
	t471 = sin(pkin(6));
	t522 = t471 * t478;
	t515 = t477 * t522;
	t432 = t435 * t473 + t515;
	t472 = sin(qJ(6));
	t476 = cos(qJ(6));
	t542 = -t432 * t472 + t476 * t487;
	t541 = t432 * t476 + t472 * t487;
	t426 = qJD(1) * t538 - (qJD(1) * t460 - t451) * t478 - t475 * t458;
	t499 = t472 * r_i_i_C(1) + t476 * r_i_i_C(2);
	t489 = qJD(6) * t499;
	t500 = t476 * r_i_i_C(1) - t472 * r_i_i_C(2);
	t495 = pkin(5) + t500;
	t532 = pkin(10) + r_i_i_C(3);
	t540 = qJD(4) - t473 * t489 + (t532 * t473 + t495 * t477) * qJD(5);
	t455 = t460 * t471;
	t537 = -t455 * t477 + t529 * t473;
	t490 = qJD(6) * t500;
	t483 = t460 * qJD(2);
	t452 = t518 * qJD(2);
	t480 = t475 * t481;
	t517 = qJD(1) * t478;
	t425 = -qJD(1) * t480 + t478 * t452 + t475 * t483 - t484 * t517;
	t523 = t471 * t475;
	t514 = qJD(1) * t523;
	t536 = (-qJD(5) * t435 + t514) * t473 - qJD(5) * t515 + t425 * t477;
	t492 = -t435 * t477 + t473 * t522;
	t416 = t492 * qJD(5) - t425 * t473 + t477 * t514;
	t482 = t495 * t473 - t532 * t477 + qJ(4);
	t533 = -pkin(3) - pkin(9);
	t530 = pkin(2) * qJD(2);
	t516 = t474 * t530;
	t513 = t471 * t517;
	t511 = t474 * t529;
	t512 = -pkin(2) * t511 + (pkin(4) + pkin(8) + qJ(3)) * t471;
	t503 = t529 * t531;
	t493 = t499 - t533;
	t438 = -t478 * t484 - t480;
	t430 = -t438 * t473 + t477 * t523;
	t491 = -t438 * t477 - t473 * t523;
	t441 = t455 * t473 + t529 * t477;
	t456 = t484 * t471;
	t470 = t531 * pkin(2) + pkin(1);
	t459 = -t471 * qJD(3) + t503 * t530;
	t450 = t471 * t483;
	t449 = qJD(2) * t456;
	t427 = t537 * qJD(5) - t449 * t473;
	t422 = t435 * qJD(1) + t475 * t452 - t478 * t483;
	t419 = t491 * qJD(5) + t422 * t473 + t477 * t513;
	t418 = t430 * qJD(5) - t422 * t477 + t473 * t513;
	t413 = t419 * t476 + t543 * t472 + (-t430 * t472 + t439 * t476) * qJD(6);
	t412 = -t419 * t472 + t543 * t476 + (-t430 * t476 - t439 * t472) * qJD(6);
	t1 = [t475 * t516 + t425 * qJ(4) + t435 * qJD(4) - t478 * t459 - t495 * t416 - t493 * t426 - t532 * t536 + (t542 * r_i_i_C(1) - t541 * r_i_i_C(2)) * qJD(6) + (-t478 * t470 - t512 * t475) * qJD(1), t438 * t490 - t493 * t422 + t482 * t543 + ((t475 * t511 - t531 * t478) * qJD(2) + (t474 * t475 - t478 * t503) * qJD(1)) * pkin(2) + t540 * t439, t513, t422, -t495 * t418 + t532 * t419 - t491 * t489, t412 * r_i_i_C(1) - t413 * r_i_i_C(2); -t478 * t516 + t419 * pkin(5) + t413 * r_i_i_C(1) + t412 * r_i_i_C(2) + t422 * qJ(4) - t438 * qJD(4) - t475 * t459 - t533 * t543 + t532 * t418 + (-t470 * t475 + t512 * t478) * qJD(1), t435 * t490 + t493 * t425 + t482 * t426 + ((-t531 * t475 - t478 * t511) * qJD(2) + (-t474 * t478 - t475 * t503) * qJD(1)) * pkin(2) - t540 * t487, t514, -t425, t532 * t416 - t492 * t489 - t495 * t536, (-t416 * t472 + t426 * t476) * r_i_i_C(1) + (-t416 * t476 - t426 * t472) * r_i_i_C(2) + (t541 * r_i_i_C(1) + t542 * r_i_i_C(2)) * qJD(6); 0, -t493 * t449 - t482 * t450 - t455 * t490 + t540 * t456 - t471 * t516, 0, t449, -t532 * t427 + t537 * t489 + t495 * (-t441 * qJD(5) + t449 * t477), (t427 * t472 - t450 * t476) * r_i_i_C(1) + (t427 * t476 + t450 * t472) * r_i_i_C(2) + ((-t441 * t476 - t456 * t472) * r_i_i_C(1) + (t441 * t472 - t456 * t476) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end