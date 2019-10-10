% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:04
	% EndTime: 2019-10-10 11:24:04
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t284 = cos(pkin(6));
	t287 = sin(qJ(1));
	t286 = sin(qJ(2));
	t307 = t287 * t286;
	t298 = t284 * t307;
	t302 = qJD(2) * t286;
	t289 = cos(qJ(2));
	t290 = cos(qJ(1));
	t304 = t290 * t289;
	t275 = -qJD(1) * t298 - t287 * t302 + (qJD(2) * t284 + qJD(1)) * t304;
	t305 = t290 * t286;
	t306 = t287 * t289;
	t277 = t284 * t305 + t306;
	t285 = sin(qJ(3));
	t288 = cos(qJ(3));
	t283 = sin(pkin(6));
	t303 = qJD(1) * t283;
	t297 = t287 * t303;
	t308 = t283 * t290;
	t311 = (-t277 * t288 + t285 * t308) * qJD(3) - t275 * t285 + t288 * t297;
	t310 = t283 * t285;
	t309 = t283 * t288;
	t301 = qJD(3) * t285;
	t300 = qJD(3) * t288;
	t299 = qJD(3) * t289;
	t296 = t290 * t303;
	t295 = t283 * qJD(2) * t289;
	t276 = t284 * t304 - t307;
	t278 = -t284 * t306 - t305;
	t293 = t298 - t304;
	t291 = -t275 * t288 + t300 * t308 + (qJD(3) * t277 - t297) * t285;
	t274 = t278 * qJD(1) - t277 * qJD(2);
	t273 = -t277 * qJD(1) + t278 * qJD(2);
	t272 = -t276 * qJD(1) + t293 * qJD(2);
	t271 = t285 * t296 + t273 * t288 + (t285 * t293 + t287 * t309) * qJD(3);
	t270 = t288 * t296 - t273 * t285 + (-t287 * t310 + t288 * t293) * qJD(3);
	t1 = [t291, t272 * t288 - t278 * t301, t270, 0, 0, 0; t271, t274 * t288 - t276 * t301, t311, 0, 0, 0; 0, (-t285 * t299 - t288 * t302) * t283, -t285 * t295 + (-t284 * t285 - t286 * t309) * qJD(3), 0, 0, 0; -t311, -t272 * t285 - t278 * t300, -t271, 0, 0, 0; t270, -t274 * t285 - t276 * t300, t291, 0, 0, 0; 0, (t285 * t302 - t288 * t299) * t283, -t288 * t295 + (-t284 * t288 + t286 * t310) * qJD(3), 0, 0, 0; t274, t273, 0, 0, 0, 0; -t272, t275, 0, 0, 0, 0; 0, t295, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:05
	% EndTime: 2019-10-10 11:24:05
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (144->36), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->38)
	t310 = cos(pkin(6));
	t312 = sin(qJ(1));
	t311 = sin(qJ(2));
	t331 = t312 * t311;
	t322 = t310 * t331;
	t326 = qJD(2) * t311;
	t313 = cos(qJ(2));
	t314 = cos(qJ(1));
	t328 = t314 * t313;
	t298 = -qJD(1) * t322 - t312 * t326 + (qJD(2) * t310 + qJD(1)) * t328;
	t329 = t314 * t311;
	t330 = t312 * t313;
	t300 = t310 * t329 + t330;
	t308 = qJ(3) + pkin(11);
	t306 = sin(t308);
	t307 = cos(t308);
	t309 = sin(pkin(6));
	t327 = qJD(1) * t309;
	t321 = t312 * t327;
	t332 = t309 * t314;
	t335 = (-t300 * t307 + t306 * t332) * qJD(3) - t298 * t306 + t307 * t321;
	t334 = t309 * t311;
	t333 = t309 * t312;
	t325 = qJD(3) * t306;
	t324 = qJD(3) * t307;
	t323 = qJD(3) * t313;
	t320 = t314 * t327;
	t319 = t309 * qJD(2) * t313;
	t299 = t310 * t328 - t331;
	t301 = -t310 * t330 - t329;
	t317 = t322 - t328;
	t315 = -t298 * t307 + t324 * t332 + (qJD(3) * t300 - t321) * t306;
	t297 = t301 * qJD(1) - t300 * qJD(2);
	t296 = -t300 * qJD(1) + t301 * qJD(2);
	t295 = -t299 * qJD(1) + t317 * qJD(2);
	t294 = t306 * t320 + t296 * t307 + (t306 * t317 + t307 * t333) * qJD(3);
	t293 = t307 * t320 - t296 * t306 + (-t306 * t333 + t307 * t317) * qJD(3);
	t1 = [t315, t295 * t307 - t301 * t325, t293, 0, 0, 0; t294, t297 * t307 - t299 * t325, t335, 0, 0, 0; 0, (-t306 * t323 - t307 * t326) * t309, -t306 * t319 + (-t306 * t310 - t307 * t334) * qJD(3), 0, 0, 0; -t335, -t295 * t306 - t301 * t324, -t294, 0, 0, 0; t293, -t297 * t306 - t299 * t324, t315, 0, 0, 0; 0, (t306 * t326 - t307 * t323) * t309, -t307 * t319 + (t306 * t334 - t307 * t310) * qJD(3), 0, 0, 0; t297, t296, 0, 0, 0, 0; -t295, t298, 0, 0, 0, 0; 0, t319, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:06
	% EndTime: 2019-10-10 11:24:06
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (235->49), mult. (518->103), div. (0->0), fcn. (536->10), ass. (0->48)
	t415 = sin(qJ(1));
	t413 = cos(pkin(6));
	t422 = qJD(2) * t413 + qJD(1);
	t414 = sin(qJ(2));
	t438 = t415 * t414;
	t428 = t413 * t438;
	t433 = qJD(2) * t414;
	t416 = cos(qJ(2));
	t417 = cos(qJ(1));
	t435 = t417 * t416;
	t394 = -qJD(1) * t428 - t415 * t433 + t422 * t435;
	t436 = t417 * t414;
	t437 = t415 * t416;
	t397 = t413 * t436 + t437;
	t409 = qJ(3) + pkin(11);
	t407 = sin(t409);
	t408 = cos(t409);
	t411 = sin(pkin(6));
	t434 = qJD(1) * t411;
	t426 = t415 * t434;
	t439 = t411 * t417;
	t390 = (t397 * t407 + t408 * t439) * qJD(3) - t394 * t408 - t407 * t426;
	t442 = t408 * t414;
	t441 = t411 * t414;
	t440 = t411 * t415;
	t432 = qJD(2) * t416;
	t431 = qJD(3) * t407;
	t430 = qJD(3) * t408;
	t429 = qJD(3) * t416;
	t427 = t413 * t435;
	t425 = t417 * t434;
	t424 = t411 * t432;
	t423 = t407 * t429;
	t398 = -t413 * t437 - t436;
	t391 = -qJD(1) * t427 - t417 * t432 + t422 * t438;
	t420 = -t391 * t408 + t398 * t431;
	t393 = t398 * qJD(1) - t397 * qJD(2);
	t396 = t427 - t438;
	t419 = -t393 * t408 + t396 * t431;
	t389 = -t394 * t407 - t397 * t430 + t408 * t426 + t431 * t439;
	t412 = cos(pkin(12));
	t410 = sin(pkin(12));
	t399 = -t428 + t435;
	t395 = -t407 * t424 + (-t407 * t413 - t408 * t441) * qJD(3);
	t392 = -t397 * qJD(1) + t398 * qJD(2);
	t388 = t407 * t425 + t392 * t408 + (-t399 * t407 + t408 * t440) * qJD(3);
	t387 = t392 * t407 - t408 * t425 + (t399 * t408 + t407 * t440) * qJD(3);
	t1 = [t390 * t412 + t393 * t410, t392 * t410 - t420 * t412, -t387 * t412, 0, 0, 0; t388 * t412 - t391 * t410, t394 * t410 - t419 * t412, t389 * t412, 0, 0, 0; 0, (-t412 * t423 + (t410 * t416 - t412 * t442) * qJD(2)) * t411, t395 * t412, 0, 0, 0; -t390 * t410 + t393 * t412, t392 * t412 + t420 * t410, t387 * t410, 0, 0, 0; -t388 * t410 - t391 * t412, t394 * t412 + t419 * t410, -t389 * t410, 0, 0, 0; 0, (t410 * t423 + (t410 * t442 + t412 * t416) * qJD(2)) * t411, -t395 * t410, 0, 0, 0; t389, t391 * t407 + t398 * t430, t388, 0, 0, 0; t387, t393 * t407 + t396 * t430, -t390, 0, 0, 0; 0, (-t407 * t433 + t408 * t429) * t411, t408 * t424 + (-t407 * t441 + t408 * t413) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:06
	% EndTime: 2019-10-10 11:24:07
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (506->75), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->68)
	t494 = sin(qJ(1));
	t492 = cos(pkin(6));
	t505 = qJD(2) * t492 + qJD(1);
	t493 = sin(qJ(2));
	t526 = t493 * t494;
	t512 = t492 * t526;
	t521 = qJD(2) * t493;
	t495 = cos(qJ(2));
	t496 = cos(qJ(1));
	t523 = t495 * t496;
	t463 = -qJD(1) * t512 - t494 * t521 + t505 * t523;
	t490 = qJ(3) + pkin(11);
	t486 = sin(t490);
	t488 = cos(t490);
	t524 = t494 * t495;
	t525 = t493 * t496;
	t474 = t492 * t525 + t524;
	t491 = sin(pkin(6));
	t527 = t491 * t496;
	t500 = t474 * t486 + t488 * t527;
	t522 = qJD(1) * t491;
	t510 = t494 * t522;
	t459 = t500 * qJD(3) - t463 * t488 - t486 * t510;
	t475 = t492 * t524 + t525;
	t462 = t475 * qJD(1) + t474 * qJD(2);
	t513 = t486 * t527;
	t468 = -t474 * t488 + t513;
	t511 = t492 * t523;
	t473 = -t511 + t526;
	t489 = pkin(12) + qJ(6);
	t485 = sin(t489);
	t487 = cos(t489);
	t539 = t459 * t487 - t462 * t485 + (-t468 * t485 - t473 * t487) * qJD(6);
	t538 = (t468 * t487 - t473 * t485) * qJD(6) + t459 * t485 + t462 * t487;
	t517 = qJD(3) * t495;
	t535 = (qJD(2) * t488 - qJD(6)) * t493 + t486 * t517;
	t530 = t491 * t493;
	t529 = t491 * t494;
	t528 = t491 * t495;
	t520 = qJD(2) * t495;
	t519 = qJD(3) * t486;
	t518 = qJD(3) * t488;
	t516 = qJD(6) * t485;
	t515 = qJD(6) * t487;
	t514 = qJD(6) * t488;
	t509 = t496 * t522;
	t508 = t491 * t521;
	t507 = t491 * t520;
	t461 = -t474 * qJD(1) - t475 * qJD(2);
	t503 = t475 * t514 + t461;
	t502 = t473 * t514 + t463;
	t501 = (qJD(2) - t514) * t495;
	t476 = -t512 + t523;
	t469 = -t476 * t486 + t488 * t529;
	t470 = t476 * t488 + t486 * t529;
	t472 = t486 * t492 + t488 * t530;
	t471 = -t486 * t530 + t488 * t492;
	t457 = qJD(3) * t513 - t463 * t486 - t474 * t518 + t488 * t510;
	t460 = -qJD(1) * t511 - t496 * t520 + t505 * t526;
	t498 = qJD(6) * t476 + t460 * t488 + t475 * t519;
	t497 = qJD(6) * t474 - t462 * t488 + t473 * t519;
	t465 = t471 * qJD(3) + t488 * t507;
	t464 = -t472 * qJD(3) - t486 * t507;
	t456 = t469 * qJD(3) + t461 * t488 + t486 * t509;
	t455 = t470 * qJD(3) + t461 * t486 - t488 * t509;
	t454 = t456 * t487 - t460 * t485 + (-t470 * t485 + t475 * t487) * qJD(6);
	t453 = -t456 * t485 - t460 * t487 + (-t470 * t487 - t475 * t485) * qJD(6);
	t1 = [t539, t503 * t485 + t498 * t487, -t455 * t487 - t469 * t516, 0, 0, t453; t454, t502 * t485 + t497 * t487, t457 * t487 + t500 * t516, 0, 0, t538; 0, (t485 * t501 - t535 * t487) * t491, t464 * t487 - t471 * t516, 0, 0, t487 * t508 - t465 * t485 + (-t472 * t487 + t485 * t528) * qJD(6); -t538, -t498 * t485 + t503 * t487, t455 * t485 - t469 * t515, 0, 0, -t454; t453, -t497 * t485 + t502 * t487, -t457 * t485 + t500 * t515, 0, 0, t539; 0, (t535 * t485 + t487 * t501) * t491, -t464 * t485 - t471 * t515, 0, 0, -t485 * t508 - t465 * t487 + (t472 * t485 + t487 * t528) * qJD(6); t457, t460 * t486 - t475 * t518, t456, 0, 0, 0; t455, -t462 * t486 - t473 * t518, -t459, 0, 0, 0; 0, (-t486 * t521 + t488 * t517) * t491, t465, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end