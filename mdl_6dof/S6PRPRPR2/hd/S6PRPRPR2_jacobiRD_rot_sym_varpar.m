% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(10));
	t58 = sin(pkin(10));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->7), mult. (60->18), div. (0->0), fcn. (60->8), ass. (0->13)
	t105 = sin(pkin(11));
	t108 = cos(pkin(11));
	t111 = sin(qJ(2));
	t112 = cos(qJ(2));
	t115 = (t105 * t112 + t108 * t111) * qJD(2);
	t103 = (t105 * t111 - t108 * t112) * qJD(2);
	t110 = cos(pkin(6));
	t109 = cos(pkin(10));
	t107 = sin(pkin(6));
	t106 = sin(pkin(10));
	t102 = t110 * t115;
	t101 = t110 * t103;
	t1 = [0, t106 * t102 + t109 * t103, 0, 0, 0, 0; 0, -t109 * t102 + t106 * t103, 0, 0, 0, 0; 0, -t107 * t115, 0, 0, 0, 0; 0, -t106 * t101 + t109 * t115, 0, 0, 0, 0; 0, t109 * t101 + t106 * t115, 0, 0, 0, 0; 0, t107 * t103, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:52
	% EndTime: 2019-10-09 21:31:52
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (81->31), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->38)
	t294 = sin(pkin(6));
	t298 = sin(qJ(4));
	t310 = t294 * t298;
	t300 = cos(qJ(4));
	t309 = t294 * t300;
	t299 = sin(qJ(2));
	t308 = qJD(2) * t299;
	t301 = cos(qJ(2));
	t307 = qJD(2) * t301;
	t306 = qJD(4) * t298;
	t305 = qJD(4) * t300;
	t292 = sin(pkin(11));
	t295 = cos(pkin(11));
	t297 = cos(pkin(6));
	t279 = (t292 * t308 - t295 * t307) * t297;
	t286 = -t292 * t307 - t295 * t308;
	t293 = sin(pkin(10));
	t296 = cos(pkin(10));
	t270 = -t296 * t279 + t293 * t286;
	t272 = t293 * t279 + t296 * t286;
	t304 = t301 * t292 + t299 * t295;
	t303 = t299 * t292 - t301 * t295;
	t281 = t303 * t294;
	t302 = qJD(2) * t304;
	t285 = t303 * qJD(2);
	t284 = t304 * t297;
	t283 = t303 * t297;
	t282 = t304 * t294;
	t280 = t297 * t302;
	t278 = qJD(2) * t281;
	t277 = t294 * t302;
	t276 = -t293 * t284 - t296 * t303;
	t275 = t293 * t283 - t296 * t304;
	t274 = t296 * t284 - t293 * t303;
	t273 = -t296 * t283 - t293 * t304;
	t271 = t293 * t280 + t296 * t285;
	t269 = -t296 * t280 + t293 * t285;
	t1 = [0, t271 * t300 - t275 * t306, 0, -t272 * t298 + (-t276 * t300 - t293 * t310) * qJD(4), 0, 0; 0, t269 * t300 - t273 * t306, 0, -t270 * t298 + (-t274 * t300 + t296 * t310) * qJD(4), 0, 0; 0, -t277 * t300 + t281 * t306, 0, t278 * t298 + (-t282 * t300 - t297 * t298) * qJD(4), 0, 0; 0, -t271 * t298 - t275 * t305, 0, -t272 * t300 + (t276 * t298 - t293 * t309) * qJD(4), 0, 0; 0, -t269 * t298 - t273 * t305, 0, -t270 * t300 + (t274 * t298 + t296 * t309) * qJD(4), 0, 0; 0, t277 * t298 + t281 * t305, 0, t278 * t300 + (t282 * t298 - t297 * t300) * qJD(4), 0, 0; 0, t272, 0, 0, 0, 0; 0, t270, 0, 0, 0, 0; 0, -t278, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:53
	% EndTime: 2019-10-09 21:31:54
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (134->39), mult. (462->94), div. (0->0), fcn. (504->12), ass. (0->46)
	t400 = sin(pkin(6));
	t405 = sin(qJ(4));
	t422 = t400 * t405;
	t407 = cos(qJ(4));
	t421 = t400 * t407;
	t406 = sin(qJ(2));
	t420 = qJD(2) * t406;
	t408 = cos(qJ(2));
	t419 = qJD(2) * t408;
	t418 = qJD(4) * t405;
	t417 = qJD(4) * t407;
	t398 = sin(pkin(11));
	t402 = cos(pkin(11));
	t404 = cos(pkin(6));
	t384 = (t398 * t420 - t402 * t419) * t404;
	t391 = -t398 * t419 - t402 * t420;
	t399 = sin(pkin(10));
	t403 = cos(pkin(10));
	t416 = -t403 * t384 + t399 * t391;
	t415 = t399 * t384 + t403 * t391;
	t414 = t408 * t398 + t406 * t402;
	t413 = t406 * t398 - t408 * t402;
	t409 = qJD(2) * t414;
	t385 = t404 * t409;
	t390 = t413 * qJD(2);
	t372 = -t403 * t385 + t399 * t390;
	t388 = t413 * t404;
	t378 = -t403 * t388 - t399 * t414;
	t412 = -t372 * t407 + t378 * t418;
	t375 = t399 * t385 + t403 * t390;
	t380 = t399 * t388 - t403 * t414;
	t411 = -t375 * t407 + t380 * t418;
	t382 = t400 * t409;
	t386 = t413 * t400;
	t410 = t382 * t407 - t386 * t418;
	t401 = cos(pkin(12));
	t397 = sin(pkin(12));
	t389 = t414 * t404;
	t387 = t414 * t400;
	t383 = qJD(2) * t386;
	t381 = -t399 * t389 - t403 * t413;
	t379 = t403 * t389 - t399 * t413;
	t371 = t383 * t405 + (-t387 * t407 - t404 * t405) * qJD(4);
	t370 = -t415 * t405 + (-t381 * t407 - t399 * t422) * qJD(4);
	t369 = -t416 * t405 + (-t379 * t407 + t403 * t422) * qJD(4);
	t1 = [0, t397 * t415 - t411 * t401, 0, t370 * t401, 0, 0; 0, t397 * t416 - t412 * t401, 0, t369 * t401, 0, 0; 0, -t383 * t397 - t410 * t401, 0, t371 * t401, 0, 0; 0, t411 * t397 + t401 * t415, 0, -t370 * t397, 0, 0; 0, t412 * t397 + t401 * t416, 0, -t369 * t397, 0, 0; 0, -t383 * t401 + t410 * t397, 0, -t371 * t397, 0, 0; 0, t375 * t405 + t380 * t417, 0, t415 * t407 + (-t381 * t405 + t399 * t421) * qJD(4), 0, 0; 0, t372 * t405 + t378 * t417, 0, t416 * t407 + (-t379 * t405 - t403 * t421) * qJD(4), 0, 0; 0, -t382 * t405 - t386 * t417, 0, -t383 * t407 + (-t387 * t405 + t404 * t407) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:54
	% EndTime: 2019-10-09 21:31:54
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (369->67), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->61)
	t484 = cos(pkin(6));
	t479 = sin(pkin(11));
	t482 = cos(pkin(11));
	t486 = sin(qJ(2));
	t488 = cos(qJ(2));
	t495 = t488 * t479 + t486 * t482;
	t466 = t495 * t484;
	t469 = t486 * t479 - t488 * t482;
	t467 = t469 * qJD(2);
	t481 = sin(pkin(6));
	t485 = sin(qJ(4));
	t513 = t481 * t485;
	t487 = cos(qJ(4));
	t512 = t481 * t487;
	t509 = qJD(2) * t486;
	t508 = qJD(2) * t488;
	t507 = qJD(4) * t485;
	t506 = qJD(4) * t487;
	t478 = pkin(12) + qJ(6);
	t476 = sin(t478);
	t505 = qJD(6) * t476;
	t477 = cos(t478);
	t504 = qJD(6) * t477;
	t503 = qJD(6) * t487;
	t480 = sin(pkin(10));
	t483 = cos(pkin(10));
	t493 = t469 * t484;
	t452 = -t480 * t495 - t483 * t493;
	t463 = (t479 * t509 - t482 * t508) * t484;
	t468 = -t479 * t508 - t482 * t509;
	t499 = -t483 * t463 + t480 * t468;
	t502 = -t452 * t503 + t499;
	t455 = t480 * t493 - t483 * t495;
	t498 = t480 * t463 + t483 * t468;
	t501 = -t455 * t503 + t498;
	t462 = t481 * t467;
	t464 = t469 * t481;
	t500 = t464 * t503 - t462;
	t465 = t495 * t481;
	t458 = t465 * t487 + t484 * t485;
	t457 = -t465 * t485 + t484 * t487;
	t497 = t483 * t466 - t480 * t469;
	t496 = -t480 * t466 - t483 * t469;
	t441 = -t483 * t512 - t485 * t497;
	t494 = t483 * t513 - t487 * t497;
	t443 = t480 * t512 - t485 * t496;
	t444 = t480 * t513 + t487 * t496;
	t492 = qJD(2) * t466;
	t445 = t480 * t467 - t483 * t492;
	t491 = -qJD(6) * t497 - t445 * t487 + t452 * t507;
	t448 = t483 * t467 + t480 * t492;
	t490 = -qJD(6) * t496 - t448 * t487 + t455 * t507;
	t461 = qJD(2) * t465;
	t489 = qJD(6) * t465 - t461 * t487 + t464 * t507;
	t440 = t457 * qJD(4) - t462 * t487;
	t439 = -t458 * qJD(4) + t462 * t485;
	t438 = t443 * qJD(4) + t487 * t498;
	t437 = -t444 * qJD(4) - t485 * t498;
	t436 = t441 * qJD(4) + t487 * t499;
	t435 = t494 * qJD(4) - t485 * t499;
	t1 = [0, t501 * t476 - t490 * t477, 0, t437 * t477 - t443 * t505, 0, -t438 * t476 - t448 * t477 + (-t444 * t477 + t455 * t476) * qJD(6); 0, t502 * t476 - t491 * t477, 0, t435 * t477 - t441 * t505, 0, -t436 * t476 - t445 * t477 + (t452 * t476 + t477 * t494) * qJD(6); 0, t500 * t476 + t489 * t477, 0, t439 * t477 - t457 * t505, 0, -t440 * t476 + t461 * t477 + (-t458 * t477 - t464 * t476) * qJD(6); 0, t490 * t476 + t501 * t477, 0, -t437 * t476 - t443 * t504, 0, -t438 * t477 + t448 * t476 + (t444 * t476 + t455 * t477) * qJD(6); 0, t491 * t476 + t502 * t477, 0, -t435 * t476 - t441 * t504, 0, -t436 * t477 + t445 * t476 + (t452 * t477 - t476 * t494) * qJD(6); 0, -t489 * t476 + t500 * t477, 0, -t439 * t476 - t457 * t504, 0, -t440 * t477 - t461 * t476 + (t458 * t476 - t464 * t477) * qJD(6); 0, t448 * t485 + t455 * t506, 0, t438, 0, 0; 0, t445 * t485 + t452 * t506, 0, t436, 0, 0; 0, -t461 * t485 - t464 * t506, 0, t440, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end