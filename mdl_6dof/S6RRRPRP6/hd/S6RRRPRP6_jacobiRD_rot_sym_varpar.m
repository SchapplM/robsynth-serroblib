% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
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
	% StartTime: 2019-10-10 11:44:00
	% EndTime: 2019-10-10 11:44:00
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
	% StartTime: 2019-10-10 11:44:01
	% EndTime: 2019-10-10 11:44:01
	% DurationCPUTime: 0.19s
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
	% StartTime: 2019-10-10 11:44:01
	% EndTime: 2019-10-10 11:44:01
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
	% StartTime: 2019-10-10 11:44:03
	% EndTime: 2019-10-10 11:44:04
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (424->74), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t487 = sin(qJ(1));
	t484 = cos(pkin(6));
	t499 = qJD(2) * t484 + qJD(1);
	t486 = sin(qJ(2));
	t520 = t487 * t486;
	t507 = t484 * t520;
	t515 = qJD(2) * t486;
	t489 = cos(qJ(2));
	t490 = cos(qJ(1));
	t517 = t490 * t489;
	t458 = -qJD(1) * t507 - t487 * t515 + t499 * t517;
	t482 = qJ(3) + pkin(11);
	t480 = sin(t482);
	t481 = cos(t482);
	t518 = t490 * t486;
	t519 = t487 * t489;
	t469 = t484 * t518 + t519;
	t483 = sin(pkin(6));
	t521 = t483 * t490;
	t494 = t469 * t480 + t481 * t521;
	t516 = qJD(1) * t483;
	t504 = t487 * t516;
	t454 = t494 * qJD(3) - t458 * t481 - t480 * t504;
	t470 = t484 * t519 + t518;
	t457 = t470 * qJD(1) + t469 * qJD(2);
	t506 = t480 * t521;
	t463 = -t469 * t481 + t506;
	t505 = t484 * t517;
	t468 = -t505 + t520;
	t485 = sin(qJ(5));
	t488 = cos(qJ(5));
	t533 = t454 * t488 - t457 * t485 + (-t463 * t485 - t468 * t488) * qJD(5);
	t532 = (t463 * t488 - t468 * t485) * qJD(5) + t454 * t485 + t457 * t488;
	t511 = qJD(3) * t489;
	t529 = (qJD(2) * t481 - qJD(5)) * t486 + t480 * t511;
	t524 = t483 * t486;
	t523 = t483 * t487;
	t522 = t483 * t489;
	t514 = qJD(2) * t489;
	t513 = qJD(3) * t480;
	t512 = qJD(3) * t481;
	t510 = qJD(5) * t481;
	t509 = qJD(5) * t485;
	t508 = qJD(5) * t488;
	t503 = t490 * t516;
	t502 = t483 * t515;
	t501 = t483 * t514;
	t456 = -t469 * qJD(1) - t470 * qJD(2);
	t497 = t470 * t510 + t456;
	t496 = t468 * t510 + t458;
	t495 = (qJD(2) - t510) * t489;
	t471 = -t507 + t517;
	t464 = -t471 * t480 + t481 * t523;
	t465 = t471 * t481 + t480 * t523;
	t467 = t484 * t480 + t481 * t524;
	t466 = -t480 * t524 + t484 * t481;
	t452 = qJD(3) * t506 - t458 * t480 - t469 * t512 + t481 * t504;
	t455 = -qJD(1) * t505 - t490 * t514 + t499 * t520;
	t492 = qJD(5) * t471 + t455 * t481 + t470 * t513;
	t491 = qJD(5) * t469 - t457 * t481 + t468 * t513;
	t460 = t466 * qJD(3) + t481 * t501;
	t459 = -t467 * qJD(3) - t480 * t501;
	t451 = t464 * qJD(3) + t456 * t481 + t480 * t503;
	t450 = t465 * qJD(3) + t456 * t480 - t481 * t503;
	t449 = t451 * t488 - t455 * t485 + (-t465 * t485 + t470 * t488) * qJD(5);
	t448 = -t451 * t485 - t455 * t488 + (-t465 * t488 - t470 * t485) * qJD(5);
	t1 = [t533, t497 * t485 + t488 * t492, -t450 * t488 - t464 * t509, 0, t448, 0; t449, t485 * t496 + t488 * t491, t452 * t488 + t494 * t509, 0, t532, 0; 0, (t485 * t495 - t529 * t488) * t483, t459 * t488 - t466 * t509, 0, t488 * t502 - t460 * t485 + (-t467 * t488 + t485 * t522) * qJD(5), 0; -t532, -t485 * t492 + t488 * t497, t450 * t485 - t464 * t508, 0, -t449, 0; t448, -t485 * t491 + t488 * t496, -t452 * t485 + t494 * t508, 0, t533, 0; 0, (t529 * t485 + t488 * t495) * t483, -t459 * t485 - t466 * t508, 0, -t485 * t502 - t460 * t488 + (t467 * t485 + t488 * t522) * qJD(5), 0; t452, t455 * t480 - t470 * t512, t451, 0, 0, 0; t450, -t457 * t480 - t468 * t512, -t454, 0, 0, 0; 0, (-t480 * t515 + t481 * t511) * t483, t460, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:44:03
	% EndTime: 2019-10-10 11:44:04
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (424->74), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t500 = sin(qJ(1));
	t497 = cos(pkin(6));
	t512 = qJD(2) * t497 + qJD(1);
	t499 = sin(qJ(2));
	t533 = t500 * t499;
	t520 = t497 * t533;
	t528 = qJD(2) * t499;
	t502 = cos(qJ(2));
	t503 = cos(qJ(1));
	t530 = t503 * t502;
	t471 = -qJD(1) * t520 - t500 * t528 + t512 * t530;
	t495 = qJ(3) + pkin(11);
	t493 = sin(t495);
	t494 = cos(t495);
	t531 = t503 * t499;
	t532 = t500 * t502;
	t482 = t497 * t531 + t532;
	t496 = sin(pkin(6));
	t534 = t496 * t503;
	t507 = t482 * t493 + t494 * t534;
	t529 = qJD(1) * t496;
	t517 = t500 * t529;
	t467 = t507 * qJD(3) - t471 * t494 - t493 * t517;
	t483 = t497 * t532 + t531;
	t470 = t483 * qJD(1) + t482 * qJD(2);
	t519 = t493 * t534;
	t476 = -t482 * t494 + t519;
	t518 = t497 * t530;
	t481 = -t518 + t533;
	t498 = sin(qJ(5));
	t501 = cos(qJ(5));
	t546 = t467 * t501 - t470 * t498 + (-t476 * t498 - t481 * t501) * qJD(5);
	t545 = (t476 * t501 - t481 * t498) * qJD(5) + t467 * t498 + t470 * t501;
	t524 = qJD(3) * t502;
	t542 = (qJD(2) * t494 - qJD(5)) * t499 + t493 * t524;
	t537 = t496 * t499;
	t536 = t496 * t500;
	t535 = t496 * t502;
	t527 = qJD(2) * t502;
	t526 = qJD(3) * t493;
	t525 = qJD(3) * t494;
	t523 = qJD(5) * t494;
	t522 = qJD(5) * t498;
	t521 = qJD(5) * t501;
	t516 = t503 * t529;
	t515 = t496 * t528;
	t514 = t496 * t527;
	t469 = -t482 * qJD(1) - t483 * qJD(2);
	t510 = t483 * t523 + t469;
	t509 = t481 * t523 + t471;
	t508 = (qJD(2) - t523) * t502;
	t484 = -t520 + t530;
	t477 = -t484 * t493 + t494 * t536;
	t478 = t484 * t494 + t493 * t536;
	t480 = t497 * t493 + t494 * t537;
	t479 = -t493 * t537 + t497 * t494;
	t465 = qJD(3) * t519 - t471 * t493 - t482 * t525 + t494 * t517;
	t468 = -qJD(1) * t518 - t503 * t527 + t512 * t533;
	t505 = qJD(5) * t484 + t468 * t494 + t483 * t526;
	t504 = qJD(5) * t482 - t470 * t494 + t481 * t526;
	t473 = t479 * qJD(3) + t494 * t514;
	t472 = -t480 * qJD(3) - t493 * t514;
	t464 = t477 * qJD(3) + t469 * t494 + t493 * t516;
	t463 = t478 * qJD(3) + t469 * t493 - t494 * t516;
	t462 = t464 * t501 - t468 * t498 + (-t478 * t498 + t483 * t501) * qJD(5);
	t461 = -t464 * t498 - t468 * t501 + (-t478 * t501 - t483 * t498) * qJD(5);
	t1 = [t546, t510 * t498 + t505 * t501, -t463 * t501 - t477 * t522, 0, t461, 0; t462, t509 * t498 + t504 * t501, t465 * t501 + t507 * t522, 0, t545, 0; 0, (t498 * t508 - t542 * t501) * t496, t472 * t501 - t479 * t522, 0, t501 * t515 - t473 * t498 + (-t480 * t501 + t498 * t535) * qJD(5), 0; -t545, -t505 * t498 + t510 * t501, t463 * t498 - t477 * t521, 0, -t462, 0; t461, -t504 * t498 + t509 * t501, -t465 * t498 + t507 * t521, 0, t546, 0; 0, (t542 * t498 + t501 * t508) * t496, -t472 * t498 - t479 * t521, 0, -t498 * t515 - t473 * t501 + (t480 * t498 + t501 * t535) * qJD(5), 0; t465, t468 * t493 - t483 * t525, t464, 0, 0, 0; t463, -t470 * t493 - t481 * t525, -t467, 0, 0, 0; 0, (-t493 * t528 + t494 * t524) * t496, t473, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end