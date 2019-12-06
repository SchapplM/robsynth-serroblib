% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(5));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(5));
	t60 = cos(pkin(9));
	t58 = sin(pkin(9));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0; 0, -t62 * t64, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0; 0, -t63 * t64, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:16
	% EndTime: 2019-12-05 15:52:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->7), mult. (60->18), div. (0->0), fcn. (60->8), ass. (0->13)
	t105 = sin(pkin(10));
	t108 = cos(pkin(10));
	t111 = sin(qJ(2));
	t112 = cos(qJ(2));
	t115 = (t105 * t112 + t108 * t111) * qJD(2);
	t103 = (t105 * t111 - t108 * t112) * qJD(2);
	t110 = cos(pkin(5));
	t109 = cos(pkin(9));
	t107 = sin(pkin(5));
	t106 = sin(pkin(9));
	t102 = t110 * t115;
	t101 = t110 * t103;
	t1 = [0, t106 * t102 + t109 * t103, 0, 0, 0; 0, -t109 * t102 + t106 * t103, 0, 0, 0; 0, -t107 * t115, 0, 0, 0; 0, -t106 * t101 + t109 * t115, 0, 0, 0; 0, t109 * t101 + t106 * t115, 0, 0, 0; 0, t107 * t103, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:17
	% EndTime: 2019-12-05 15:52:17
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (81->31), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->38)
	t294 = sin(pkin(5));
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
	t292 = sin(pkin(10));
	t295 = cos(pkin(10));
	t297 = cos(pkin(5));
	t279 = (t292 * t308 - t295 * t307) * t297;
	t286 = -t292 * t307 - t295 * t308;
	t293 = sin(pkin(9));
	t296 = cos(pkin(9));
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
	t1 = [0, t271 * t300 - t275 * t306, 0, -t272 * t298 + (-t276 * t300 - t293 * t310) * qJD(4), 0; 0, t269 * t300 - t273 * t306, 0, -t270 * t298 + (-t274 * t300 + t296 * t310) * qJD(4), 0; 0, -t277 * t300 + t281 * t306, 0, t278 * t298 + (-t282 * t300 - t297 * t298) * qJD(4), 0; 0, -t271 * t298 - t275 * t305, 0, -t272 * t300 + (t276 * t298 - t293 * t309) * qJD(4), 0; 0, -t269 * t298 - t273 * t305, 0, -t270 * t300 + (t274 * t298 + t296 * t309) * qJD(4), 0; 0, t277 * t298 + t281 * t305, 0, t278 * t300 + (t282 * t298 - t297 * t300) * qJD(4), 0; 0, t272, 0, 0, 0; 0, t270, 0, 0, 0; 0, -t278, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:52:19
	% EndTime: 2019-12-05 15:52:19
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (303->66), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->60)
	t463 = cos(pkin(5));
	t458 = sin(pkin(10));
	t461 = cos(pkin(10));
	t466 = sin(qJ(2));
	t469 = cos(qJ(2));
	t476 = t469 * t458 + t466 * t461;
	t448 = t476 * t463;
	t451 = t466 * t458 - t469 * t461;
	t449 = t451 * qJD(2);
	t460 = sin(pkin(5));
	t465 = sin(qJ(4));
	t494 = t460 * t465;
	t468 = cos(qJ(4));
	t493 = t460 * t468;
	t490 = qJD(2) * t466;
	t489 = qJD(2) * t469;
	t488 = qJD(4) * t465;
	t487 = qJD(4) * t468;
	t464 = sin(qJ(5));
	t486 = qJD(5) * t464;
	t467 = cos(qJ(5));
	t485 = qJD(5) * t467;
	t484 = qJD(5) * t468;
	t459 = sin(pkin(9));
	t462 = cos(pkin(9));
	t474 = t451 * t463;
	t434 = -t459 * t476 - t462 * t474;
	t445 = (t458 * t490 - t461 * t489) * t463;
	t450 = -t458 * t489 - t461 * t490;
	t480 = -t462 * t445 + t459 * t450;
	t483 = -t434 * t484 + t480;
	t437 = t459 * t474 - t462 * t476;
	t479 = t459 * t445 + t462 * t450;
	t482 = -t437 * t484 + t479;
	t444 = t460 * t449;
	t446 = t451 * t460;
	t481 = t446 * t484 - t444;
	t447 = t476 * t460;
	t440 = t447 * t468 + t463 * t465;
	t439 = -t447 * t465 + t463 * t468;
	t478 = t462 * t448 - t459 * t451;
	t477 = -t459 * t448 - t462 * t451;
	t423 = -t462 * t493 - t465 * t478;
	t475 = t462 * t494 - t468 * t478;
	t425 = t459 * t493 - t465 * t477;
	t426 = t459 * t494 + t468 * t477;
	t473 = qJD(2) * t448;
	t427 = t459 * t449 - t462 * t473;
	t472 = -qJD(5) * t478 - t427 * t468 + t434 * t488;
	t430 = t462 * t449 + t459 * t473;
	t471 = -qJD(5) * t477 - t430 * t468 + t437 * t488;
	t443 = qJD(2) * t447;
	t470 = qJD(5) * t447 - t443 * t468 + t446 * t488;
	t422 = t439 * qJD(4) - t444 * t468;
	t421 = -t440 * qJD(4) + t444 * t465;
	t420 = t425 * qJD(4) + t468 * t479;
	t419 = -t426 * qJD(4) - t465 * t479;
	t418 = t423 * qJD(4) + t468 * t480;
	t417 = t475 * qJD(4) - t465 * t480;
	t1 = [0, t482 * t464 - t471 * t467, 0, t419 * t467 - t425 * t486, -t420 * t464 - t430 * t467 + (-t426 * t467 + t437 * t464) * qJD(5); 0, t483 * t464 - t472 * t467, 0, t417 * t467 - t423 * t486, -t418 * t464 - t427 * t467 + (t434 * t464 + t467 * t475) * qJD(5); 0, t481 * t464 + t470 * t467, 0, t421 * t467 - t439 * t486, -t422 * t464 + t443 * t467 + (-t440 * t467 - t446 * t464) * qJD(5); 0, t471 * t464 + t482 * t467, 0, -t419 * t464 - t425 * t485, -t420 * t467 + t430 * t464 + (t426 * t464 + t437 * t467) * qJD(5); 0, t472 * t464 + t483 * t467, 0, -t417 * t464 - t423 * t485, -t418 * t467 + t427 * t464 + (t434 * t467 - t464 * t475) * qJD(5); 0, -t470 * t464 + t481 * t467, 0, -t421 * t464 - t439 * t485, -t422 * t467 - t443 * t464 + (t440 * t464 - t446 * t467) * qJD(5); 0, t430 * t465 + t437 * t487, 0, t420, 0; 0, t427 * t465 + t434 * t487, 0, t418, 0; 0, -t443 * t465 - t446 * t487, 0, t422, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end