% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(5));
	t151 = t136 * (pkin(7) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(5));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (90->45), mult. (273->84), div. (0->0), fcn. (260->8), ass. (0->34)
	t190 = cos(pkin(5));
	t187 = sin(pkin(10));
	t189 = cos(pkin(10));
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
	t188 = sin(pkin(5));
	t204 = qJD(1) * t188;
	t203 = qJD(2) * t191;
	t202 = qJD(2) * t193;
	t201 = pkin(2) * t203;
	t200 = t190 * t202;
	t180 = t187 * t191 - t193 * t189;
	t199 = -pkin(2) * t193 + t180 * r_i_i_C(1) - pkin(1);
	t172 = t187 * t190 * t203 - t189 * t200;
	t178 = -t187 * t202 - t189 * t203;
	t198 = t172 * t194 - t192 * t178;
	t197 = t192 * t172 + t178 * t194;
	t195 = t175 * r_i_i_C(1) + t190 * t209 + (-r_i_i_C(3) - pkin(7) - qJ(3)) * t188;
	t179 = pkin(2) * t200 - qJD(3) * t188;
	t177 = t180 * qJD(2);
	t174 = t180 * t190;
	t173 = qJD(2) * t175;
	t171 = -t173 * t194 + t192 * t177 + (t174 * t192 - t194 * t196) * qJD(1);
	t170 = t192 * t173 + t177 * t194 + (t174 * t194 + t192 * t196) * qJD(1);
	t1 = [t198 * r_i_i_C(1) - t171 * r_i_i_C(2) + t192 * t201 - t194 * t179 + (t195 * t192 + t199 * t194) * qJD(1), t170 * r_i_i_C(1) + ((t175 * t194 - t180 * t192) * qJD(1) - t197) * r_i_i_C(2) + ((t190 * t208 - t205) * qJD(2) + (-t190 * t205 + t208) * qJD(1)) * pkin(2), t194 * t204, 0, 0; t197 * r_i_i_C(1) + t170 * r_i_i_C(2) - t194 * t201 - t192 * t179 + (t199 * t192 - t195 * t194) * qJD(1), t171 * r_i_i_C(1) + ((t175 * t192 + t180 * t194) * qJD(1) + t198) * r_i_i_C(2) + ((-t190 * t207 - t206) * qJD(2) + (-t190 * t206 - t207) * qJD(1)) * pkin(2), t192 * t204, 0, 0; 0, (-t196 * r_i_i_C(1) + t180 * r_i_i_C(2) - t209) * t188 * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:06
	% EndTime: 2019-12-31 20:28:06
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (297->77), mult. (888->135), div. (0->0), fcn. (922->10), ass. (0->59)
	t316 = cos(pkin(5));
	t313 = sin(pkin(10));
	t315 = cos(pkin(10));
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
	t314 = sin(pkin(5));
	t348 = t314 * t322;
	t333 = qJD(4) * t348;
	t353 = t286 - t333;
	t351 = r_i_i_C(3) + pkin(8);
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
	t283 = -qJD(1) * t288 + t319 * t294 + t322 * t302;
	t312 = t321 * pkin(2) + pkin(1);
	t308 = t320 * t333;
	t303 = -t314 * qJD(3) + t340 * t350;
	t301 = t318 * t350 + (-pkin(7) - qJ(3)) * t314;
	t298 = t328 * t314;
	t295 = t316 * t324;
	t292 = t352 * t314;
	t285 = qJD(1) * t330 - t322 * t295 + t319 * t323;
	t282 = qJD(1) * t331 - t319 * t295 - t322 * t323;
	t280 = t317 * t336 + t283 * t320 + (t317 * t329 + t320 * t349) * qJD(4);
	t279 = t320 * t336 - t283 * t317 + (-t317 * t349 + t320 * t329) * qJD(4);
	t1 = [(-t286 * t320 + t317 * t339 + t308) * r_i_i_C(1) + (t353 * t317 + t320 * t339) * r_i_i_C(2) - t286 * pkin(3) + t319 * t338 - t322 * t303 + t351 * t285 + (-t322 * t312 + (-t314 * t332 + t301) * t319) * qJD(1), t351 * t283 - t330 * t326 - t327 * t282 + ((t316 * t347 - t343) * qJD(2) + (-t316 * t343 + t347) * qJD(1)) * pkin(2), t336, t279 * r_i_i_C(1) - t280 * r_i_i_C(2), 0; -t322 * t338 + t283 * pkin(3) + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) - t319 * t303 + t351 * t282 + (-t301 * t322 - t312 * t319) * qJD(1), -t351 * (qJD(1) * t329 + t322 * t294 - t345) - t331 * t326 + t327 * t285 + ((-t316 * t346 - t344) * qJD(2) + (-t316 * t344 - t346) * qJD(1)) * pkin(2), t337, t308 * r_i_i_C(2) + (r_i_i_C(1) * t325 - t286 * r_i_i_C(2)) * t320 + (-t353 * r_i_i_C(1) - t325 * r_i_i_C(2)) * t317, 0; 0, -t292 * t351 + (t304 * t326 - t324 * t327 - t338) * t314, 0, t332 * t292 + ((-t298 * t320 - t316 * t317) * r_i_i_C(1) + (t298 * t317 - t316 * t320) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:07
	% EndTime: 2019-12-31 20:28:08
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (845->121), mult. (2486->210), div. (0->0), fcn. (2730->12), ass. (0->77)
	t453 = cos(pkin(5));
	t451 = sin(pkin(10));
	t506 = cos(pkin(10));
	t508 = cos(qJ(2));
	t480 = t508 * t506;
	t456 = sin(qJ(2));
	t494 = qJD(2) * t456;
	t512 = -qJD(2) * t480 + t451 * t494;
	t431 = t512 * t453;
	t483 = t456 * t506;
	t469 = t508 * t451 + t483;
	t436 = t469 * t453;
	t484 = t508 * qJD(2);
	t438 = -qJD(2) * t483 - t451 * t484;
	t457 = sin(qJ(1));
	t460 = cos(qJ(1));
	t468 = -t456 * t451 + t480;
	t496 = qJD(1) * t457;
	t407 = t436 * t496 - t457 * t438 + (-qJD(1) * t468 + t431) * t460;
	t455 = sin(qJ(4));
	t459 = cos(qJ(4));
	t477 = t460 * t436 + t457 * t468;
	t452 = sin(pkin(5));
	t487 = t452 * t496;
	t500 = t452 * t460;
	t490 = t459 * t500;
	t401 = (-qJD(4) * t477 + t487) * t455 - qJD(4) * t490 - t407 * t459;
	t432 = qJD(2) * t436;
	t495 = qJD(1) * t460;
	t467 = t468 * t453;
	t513 = qJD(1) * t467 + t468 * qJD(2);
	t408 = -t460 * t432 - t513 * t457 - t469 * t495;
	t454 = sin(qJ(5));
	t458 = cos(qJ(5));
	t519 = t401 * t454 + t408 * t458;
	t518 = -t401 * t458 + t408 * t454;
	t414 = t455 * t500 - t459 * t477;
	t418 = -t457 * t469 + t460 * t467;
	t517 = -t414 * t454 + t418 * t458;
	t516 = t414 * t458 + t418 * t454;
	t472 = qJD(5) * (t454 * r_i_i_C(1) + t458 * r_i_i_C(2));
	t475 = t458 * r_i_i_C(1) - t454 * r_i_i_C(2) + pkin(4);
	t509 = r_i_i_C(3) + pkin(9);
	t515 = (t475 * t455 - t509 * t459) * qJD(4) + t459 * t472;
	t510 = t509 * t455 + t475 * t459 + pkin(3);
	t507 = pkin(2) * t453;
	t501 = t452 * t457;
	t498 = t456 * t457;
	t497 = t456 * t460;
	t493 = qJD(5) * t454;
	t492 = qJD(5) * t458;
	t491 = pkin(2) * t494;
	t489 = t508 * t457;
	t488 = t508 * t460;
	t486 = t452 * t495;
	t435 = t469 * t452;
	t424 = t435 * t459 + t453 * t455;
	t478 = -t435 * t455 + t453 * t459;
	t476 = -t457 * t436 + t460 * t468;
	t473 = -t455 * t476 + t459 * t501;
	t416 = t455 * t501 + t459 * t476;
	t463 = -t477 * qJD(1) + t457 * t431 + t460 * t438;
	t462 = t414 * qJD(4) + t407 * t455 + t459 * t487;
	t450 = t508 * pkin(2) + pkin(1);
	t439 = -t452 * qJD(3) + t484 * t507;
	t437 = t456 * t507 + (-pkin(7) - qJ(3)) * t452;
	t434 = t468 * t452;
	t430 = qJD(2) * t435;
	t429 = t512 * t452;
	t421 = -t457 * t467 - t460 * t469;
	t411 = t478 * qJD(4) - t429 * t459;
	t405 = -t457 * t432 + t513 * t460 - t469 * t496;
	t399 = t473 * qJD(4) + t455 * t486 + t459 * t463;
	t398 = t416 * qJD(4) + t455 * t463 - t459 * t486;
	t397 = t399 * t458 + t405 * t454 + (-t416 * t454 - t421 * t458) * qJD(5);
	t396 = -t399 * t454 + t405 * t458 + (-t416 * t458 + t421 * t454) * qJD(5);
	t1 = [t518 * r_i_i_C(1) + t519 * r_i_i_C(2) - t401 * pkin(4) + t407 * pkin(3) + t408 * pkin(8) + t457 * t491 - t460 * t439 + t509 * t462 + (t517 * r_i_i_C(1) - t516 * r_i_i_C(2)) * qJD(5) + (t457 * t437 - t460 * t450) * qJD(1), (t454 * t463 + t476 * t492) * r_i_i_C(1) + (t458 * t463 - t476 * t493) * r_i_i_C(2) + t463 * pkin(8) - t510 * t405 + ((t453 * t498 - t488) * qJD(2) + (-t453 * t488 + t498) * qJD(1)) * pkin(2) - t515 * t421, t486, -t475 * t398 + t509 * t399 - t473 * t472, t396 * r_i_i_C(1) - t397 * r_i_i_C(2); -t460 * t491 + t463 * pkin(3) + t399 * pkin(4) + t405 * pkin(8) + t397 * r_i_i_C(1) + t396 * r_i_i_C(2) - t457 * t439 + t509 * t398 + (-t437 * t460 - t450 * t457) * qJD(1), (-t407 * t454 + t477 * t492) * r_i_i_C(1) + (-t407 * t458 - t477 * t493) * r_i_i_C(2) - t407 * pkin(8) + t510 * t408 + ((-t453 * t497 - t489) * qJD(2) + (-t453 * t489 - t497) * qJD(1)) * pkin(2) - t515 * t418, t487, t509 * t401 - (-t455 * t477 - t490) * t472 + t475 * t462, -t519 * r_i_i_C(1) + t518 * r_i_i_C(2) + (t516 * r_i_i_C(1) + t517 * r_i_i_C(2)) * qJD(5); 0, (-t429 * t454 + t435 * t492) * r_i_i_C(1) + (-t429 * t458 - t435 * t493) * r_i_i_C(2) - t429 * pkin(8) - t452 * t491 - t510 * t430 - t515 * t434, 0, t509 * t411 - t478 * t472 + t475 * (-t424 * qJD(4) + t429 * t455), (-t411 * t454 + t430 * t458) * r_i_i_C(1) + (-t411 * t458 - t430 * t454) * r_i_i_C(2) + ((-t424 * t458 + t434 * t454) * r_i_i_C(1) + (t424 * t454 + t434 * t458) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end