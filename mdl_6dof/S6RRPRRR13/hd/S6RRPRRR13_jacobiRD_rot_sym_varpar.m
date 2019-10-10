% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:56
	% EndTime: 2019-10-10 11:08:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:57
	% EndTime: 2019-10-10 11:08:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:58
	% EndTime: 2019-10-10 11:08:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->11), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t201 = sin(qJ(2));
	t202 = sin(qJ(1));
	t214 = t201 * t202;
	t204 = cos(qJ(1));
	t213 = t201 * t204;
	t203 = cos(qJ(2));
	t212 = t202 * t203;
	t211 = t203 * t204;
	t199 = sin(pkin(6));
	t210 = qJD(1) * t199;
	t209 = qJD(2) * t199;
	t200 = cos(pkin(6));
	t208 = t200 * t211 - t214;
	t207 = t200 * t212 + t213;
	t206 = t200 * t213 + t212;
	t205 = -t200 * t214 + t211;
	t198 = t205 * qJD(1) + t208 * qJD(2);
	t197 = t207 * qJD(1) + t206 * qJD(2);
	t196 = t206 * qJD(1) + t207 * qJD(2);
	t195 = t208 * qJD(1) + t205 * qJD(2);
	t1 = [-t202 * t210, 0, 0, 0, 0, 0; t204 * t210, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t198, t195, 0, 0, 0, 0; t196, t197, 0, 0, 0, 0; 0, t201 * t209, 0, 0, 0, 0; -t197, -t196, 0, 0, 0, 0; t195, t198, 0, 0, 0, 0; 0, t203 * t209, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:08:58
	% EndTime: 2019-10-10 11:08:58
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (95->37), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->36)
	t290 = cos(pkin(6));
	t292 = sin(qJ(2));
	t296 = cos(qJ(1));
	t310 = t296 * t292;
	t293 = sin(qJ(1));
	t295 = cos(qJ(2));
	t311 = t293 * t295;
	t283 = t290 * t310 + t311;
	t284 = t290 * t311 + t310;
	t280 = t284 * qJD(1) + t283 * qJD(2);
	t291 = sin(qJ(4));
	t294 = cos(qJ(4));
	t309 = t296 * t295;
	t304 = t290 * t309;
	t312 = t293 * t292;
	t299 = t304 - t312;
	t289 = sin(pkin(6));
	t308 = qJD(1) * t289;
	t303 = t293 * t308;
	t313 = t289 * t296;
	t316 = (t291 * t299 + t294 * t313) * qJD(4) + t280 * t294 - t291 * t303;
	t315 = t289 * t293;
	t314 = t289 * t295;
	t307 = qJD(2) * t295;
	t306 = qJD(4) * t291;
	t305 = qJD(4) * t294;
	t302 = t296 * t308;
	t301 = t289 * qJD(2) * t292;
	t285 = -t290 * t312 + t309;
	t297 = -t294 * t303 - t280 * t291 + (-t291 * t313 + t294 * t299) * qJD(4);
	t281 = t285 * qJD(1) + t299 * qJD(2);
	t279 = -t283 * qJD(1) - t284 * qJD(2);
	t278 = -qJD(1) * t304 - t296 * t307 + (qJD(2) * t290 + qJD(1)) * t312;
	t277 = t294 * t302 - t278 * t291 + (t284 * t294 - t291 * t315) * qJD(4);
	t276 = -t291 * t302 - t278 * t294 + (-t284 * t291 - t294 * t315) * qJD(4);
	t1 = [t297, t279 * t291 + t285 * t305, 0, t276, 0, 0; t277, t281 * t291 + t283 * t305, 0, t316, 0, 0; 0, (t291 * t307 + t292 * t305) * t289, 0, t294 * t301 + (-t290 * t294 + t291 * t314) * qJD(4), 0, 0; -t316, t279 * t294 - t285 * t306, 0, -t277, 0, 0; t276, t281 * t294 - t283 * t306, 0, t297, 0, 0; 0, (-t292 * t306 + t294 * t307) * t289, 0, -t291 * t301 + (t290 * t291 + t294 * t314) * qJD(4), 0, 0; -t281, t278, 0, 0, 0, 0; t279, -t280, 0, 0, 0, 0; 0, -t301, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:09:00
	% EndTime: 2019-10-10 11:09:01
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (289->72), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->63)
	t459 = cos(pkin(6));
	t462 = sin(qJ(2));
	t467 = cos(qJ(1));
	t492 = t467 * t462;
	t463 = sin(qJ(1));
	t466 = cos(qJ(2));
	t493 = t463 * t466;
	t451 = t459 * t492 + t493;
	t452 = t459 * t493 + t492;
	t439 = t452 * qJD(1) + t451 * qJD(2);
	t491 = t467 * t466;
	t494 = t463 * t462;
	t450 = -t459 * t491 + t494;
	t461 = sin(qJ(4));
	t465 = cos(qJ(4));
	t458 = sin(pkin(6));
	t496 = t458 * t467;
	t445 = t450 * t465 + t461 * t496;
	t490 = qJD(1) * t458;
	t480 = t463 * t490;
	t432 = t445 * qJD(4) + t439 * t461 + t465 * t480;
	t482 = t459 * t494;
	t489 = qJD(2) * t462;
	t440 = -qJD(1) * t482 - t463 * t489 + (qJD(2) * t459 + qJD(1)) * t491;
	t481 = t465 * t496;
	t446 = -t450 * t461 + t481;
	t460 = sin(qJ(5));
	t464 = cos(qJ(5));
	t507 = -t432 * t464 + (-t446 * t460 - t451 * t464) * qJD(5) - t440 * t460;
	t506 = (t446 * t464 - t451 * t460) * qJD(5) - t432 * t460 + t440 * t464;
	t503 = (qJD(2) * t461 + qJD(5)) * t466;
	t498 = t458 * t463;
	t497 = t458 * t466;
	t495 = t462 * t464;
	t488 = qJD(2) * t466;
	t487 = qJD(4) * t461;
	t486 = qJD(4) * t465;
	t485 = qJD(5) * t460;
	t484 = qJD(5) * t461;
	t483 = qJD(5) * t464;
	t479 = t467 * t490;
	t478 = t458 * t489;
	t477 = t458 * t488;
	t476 = -qJD(2) - t484;
	t471 = t482 - t491;
	t437 = t450 * qJD(1) + t471 * qJD(2);
	t474 = t471 * t484 + t437;
	t473 = -t451 * t484 - t439;
	t444 = t452 * t461 + t465 * t498;
	t443 = t452 * t465 - t461 * t498;
	t448 = -t459 * t461 - t465 * t497;
	t472 = -t459 * t465 + t461 * t497;
	t438 = -t451 * qJD(1) - t452 * qJD(2);
	t469 = -qJD(5) * t452 + t438 * t461 - t471 * t486;
	t468 = -qJD(5) * t450 + t440 * t461 + t451 * t486;
	t431 = t439 * t465 + qJD(4) * t481 + (-qJD(4) * t450 - t480) * t461;
	t442 = t472 * qJD(4) + t465 * t478;
	t441 = t448 * qJD(4) + t461 * t478;
	t435 = t443 * qJD(4) - t437 * t461 + t465 * t479;
	t434 = t444 * qJD(4) + t437 * t465 + t461 * t479;
	t430 = t435 * t464 + t438 * t460 + (-t444 * t460 - t464 * t471) * qJD(5);
	t429 = -t435 * t460 + t438 * t464 + (-t444 * t464 + t460 * t471) * qJD(5);
	t1 = [t507, t474 * t460 + t469 * t464, 0, -t434 * t464 - t443 * t485, t429, 0; t430, t473 * t460 + t468 * t464, 0, t431 * t464 - t445 * t485, t506, 0; 0, (t464 * t503 + (t476 * t460 + t464 * t486) * t462) * t458, 0, t442 * t464 - t448 * t485, t464 * t477 - t441 * t460 + (-t458 * t460 * t462 + t464 * t472) * qJD(5), 0; -t506, -t469 * t460 + t474 * t464, 0, t434 * t460 - t443 * t483, -t430, 0; t429, -t468 * t460 + t473 * t464, 0, -t431 * t460 - t445 * t483, t507, 0; 0, (t476 * t495 + (-t462 * t486 - t503) * t460) * t458, 0, -t442 * t460 - t448 * t483, -t460 * t477 - t441 * t464 + (-t458 * t495 - t460 * t472) * qJD(5), 0; t431, -t438 * t465 - t471 * t487, 0, t435, 0, 0; t434, -t440 * t465 + t451 * t487, 0, t432, 0, 0; 0, (t462 * t487 - t465 * t488) * t458, 0, t441, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:09:00
	% EndTime: 2019-10-10 11:09:01
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (521->68), mult. (1126->126), div. (0->0), fcn. (1210->10), ass. (0->70)
	t516 = qJ(5) + qJ(6);
	t513 = sin(t516);
	t514 = cos(t516);
	t518 = cos(pkin(6));
	t520 = sin(qJ(2));
	t524 = cos(qJ(1));
	t555 = t524 * t520;
	t521 = sin(qJ(1));
	t523 = cos(qJ(2));
	t556 = t521 * t523;
	t506 = t518 * t555 + t556;
	t507 = t518 * t556 + t555;
	t495 = t507 * qJD(1) + t506 * qJD(2);
	t554 = t524 * t523;
	t557 = t521 * t520;
	t505 = -t518 * t554 + t557;
	t519 = sin(qJ(4));
	t522 = cos(qJ(4));
	t517 = sin(pkin(6));
	t558 = t517 * t524;
	t500 = t505 * t522 + t519 * t558;
	t553 = qJD(1) * t517;
	t546 = t521 * t553;
	t487 = t500 * qJD(4) + t495 * t519 + t522 * t546;
	t515 = qJD(5) + qJD(6);
	t542 = -t506 * t515 - t487;
	t548 = t518 * t557;
	t552 = qJD(2) * t520;
	t496 = -qJD(1) * t548 - t521 * t552 + (qJD(2) * t518 + qJD(1)) * t554;
	t547 = t522 * t558;
	t566 = -(-t505 * t519 + t547) * t515 - t496;
	t480 = t542 * t513 - t514 * t566;
	t481 = t513 * t566 + t542 * t514;
	t494 = -t506 * qJD(1) - t507 * qJD(2);
	t560 = t517 * t521;
	t532 = t507 * t519 + t522 * t560;
	t567 = -t532 * t515 + t494;
	t549 = qJD(4) * t522;
	t565 = (qJD(2) * t519 + t515) * t523 + t520 * t549;
	t563 = t513 * t515;
	t562 = t514 * t515;
	t561 = t515 * t519;
	t559 = t517 * t523;
	t551 = qJD(2) * t523;
	t550 = qJD(4) * t519;
	t545 = t524 * t553;
	t544 = t517 * t552;
	t530 = t548 - t554;
	t493 = t505 * qJD(1) + t530 * qJD(2);
	t499 = t507 * t522 - t519 * t560;
	t490 = t499 * qJD(4) - t493 * t519 + t522 * t545;
	t540 = -t515 * t530 + t490;
	t536 = t530 * t561 + t493;
	t535 = -t506 * t561 - t495;
	t503 = -t518 * t519 - t522 * t559;
	t497 = t503 * qJD(4) + t519 * t544;
	t534 = -t515 * t517 * t520 - t497;
	t533 = (-qJD(2) - t561) * t520;
	t531 = -t518 * t522 + t519 * t559;
	t529 = t515 * t531 + t517 * t551;
	t526 = t494 * t519 - t507 * t515 - t530 * t549;
	t525 = t496 * t519 - t505 * t515 + t506 * t549;
	t486 = t495 * t522 + qJD(4) * t547 + (-qJD(4) * t505 - t546) * t519;
	t498 = t531 * qJD(4) + t522 * t544;
	t489 = t532 * qJD(4) + t493 * t522 + t519 * t545;
	t485 = -t529 * t513 + t534 * t514;
	t484 = t534 * t513 + t529 * t514;
	t483 = t567 * t513 + t540 * t514;
	t482 = -t540 * t513 + t567 * t514;
	t1 = [t481, t536 * t513 + t526 * t514, 0, -t489 * t514 - t499 * t563, t482, t482; t483, t535 * t513 + t525 * t514, 0, t486 * t514 - t500 * t563, t480, t480; 0, (t513 * t533 + t565 * t514) * t517, 0, t498 * t514 - t503 * t563, t484, t484; -t480, -t526 * t513 + t536 * t514, 0, t489 * t513 - t499 * t562, -t483, -t483; t482, -t525 * t513 + t514 * t535, 0, -t486 * t513 - t500 * t562, t481, t481; 0, (-t565 * t513 + t514 * t533) * t517, 0, -t498 * t513 - t503 * t562, t485, t485; t486, -t494 * t522 - t530 * t550, 0, t490, 0, 0; t489, -t496 * t522 + t506 * t550, 0, t487, 0, 0; 0, (t520 * t550 - t522 * t551) * t517, 0, t497, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end