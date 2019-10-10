% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR1
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
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
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
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
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
	% StartTime: 2019-10-09 21:30:00
	% EndTime: 2019-10-09 21:30:00
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-09 21:30:00
	% EndTime: 2019-10-09 21:30:00
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (111->32), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->39)
	t312 = sin(pkin(10));
	t313 = sin(pkin(6));
	t327 = t312 * t313;
	t315 = cos(pkin(10));
	t326 = t313 * t315;
	t317 = sin(qJ(2));
	t325 = qJD(2) * t317;
	t318 = cos(qJ(2));
	t324 = qJD(2) * t318;
	t310 = qJ(4) + pkin(12);
	t308 = sin(t310);
	t323 = qJD(4) * t308;
	t309 = cos(t310);
	t322 = qJD(4) * t309;
	t311 = sin(pkin(11));
	t314 = cos(pkin(11));
	t316 = cos(pkin(6));
	t295 = (t311 * t325 - t314 * t324) * t316;
	t302 = -t311 * t324 - t314 * t325;
	t286 = -t295 * t315 + t302 * t312;
	t288 = t295 * t312 + t302 * t315;
	t321 = t311 * t318 + t314 * t317;
	t320 = t311 * t317 - t314 * t318;
	t297 = t320 * t313;
	t319 = qJD(2) * t321;
	t301 = t320 * qJD(2);
	t300 = t321 * t316;
	t299 = t320 * t316;
	t298 = t321 * t313;
	t296 = t316 * t319;
	t294 = qJD(2) * t297;
	t293 = t313 * t319;
	t292 = -t300 * t312 - t315 * t320;
	t291 = t299 * t312 - t315 * t321;
	t290 = t300 * t315 - t312 * t320;
	t289 = -t299 * t315 - t312 * t321;
	t287 = t296 * t312 + t301 * t315;
	t285 = -t296 * t315 + t301 * t312;
	t1 = [0, t287 * t309 - t291 * t323, 0, -t288 * t308 + (-t292 * t309 - t308 * t327) * qJD(4), 0, 0; 0, t285 * t309 - t289 * t323, 0, -t286 * t308 + (-t290 * t309 + t308 * t326) * qJD(4), 0, 0; 0, -t293 * t309 + t297 * t323, 0, t294 * t308 + (-t298 * t309 - t308 * t316) * qJD(4), 0, 0; 0, -t287 * t308 - t291 * t322, 0, -t288 * t309 + (t292 * t308 - t309 * t327) * qJD(4), 0, 0; 0, -t285 * t308 - t289 * t322, 0, -t286 * t309 + (t290 * t308 + t309 * t326) * qJD(4), 0, 0; 0, t293 * t308 + t297 * t322, 0, t294 * t309 + (t298 * t308 - t309 * t316) * qJD(4), 0, 0; 0, t288, 0, 0, 0, 0; 0, t286, 0, 0, 0, 0; 0, -t294, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:30:02
	% EndTime: 2019-10-09 21:30:03
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (396->67), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->61)
	t486 = cos(pkin(6));
	t481 = sin(pkin(11));
	t484 = cos(pkin(11));
	t488 = sin(qJ(2));
	t490 = cos(qJ(2));
	t497 = t490 * t481 + t488 * t484;
	t468 = t497 * t486;
	t471 = t488 * t481 - t490 * t484;
	t469 = t471 * qJD(2);
	t482 = sin(pkin(10));
	t483 = sin(pkin(6));
	t515 = t482 * t483;
	t485 = cos(pkin(10));
	t514 = t483 * t485;
	t511 = qJD(2) * t488;
	t510 = qJD(2) * t490;
	t480 = qJ(4) + pkin(12);
	t478 = sin(t480);
	t509 = qJD(4) * t478;
	t479 = cos(t480);
	t508 = qJD(4) * t479;
	t507 = qJD(6) * t479;
	t487 = sin(qJ(6));
	t506 = qJD(6) * t487;
	t489 = cos(qJ(6));
	t505 = qJD(6) * t489;
	t495 = t471 * t486;
	t454 = -t482 * t497 - t485 * t495;
	t465 = (t481 * t511 - t484 * t510) * t486;
	t470 = -t481 * t510 - t484 * t511;
	t501 = -t485 * t465 + t482 * t470;
	t504 = -t454 * t507 + t501;
	t457 = t482 * t495 - t485 * t497;
	t500 = t482 * t465 + t485 * t470;
	t503 = -t457 * t507 + t500;
	t464 = t483 * t469;
	t466 = t471 * t483;
	t502 = t466 * t507 - t464;
	t467 = t497 * t483;
	t460 = t467 * t479 + t486 * t478;
	t459 = -t467 * t478 + t486 * t479;
	t499 = t485 * t468 - t482 * t471;
	t498 = -t482 * t468 - t485 * t471;
	t443 = -t478 * t499 - t479 * t514;
	t496 = t478 * t514 - t479 * t499;
	t445 = -t478 * t498 + t479 * t515;
	t446 = t478 * t515 + t479 * t498;
	t494 = qJD(2) * t468;
	t447 = t482 * t469 - t485 * t494;
	t493 = -qJD(6) * t499 - t447 * t479 + t454 * t509;
	t450 = t485 * t469 + t482 * t494;
	t492 = -qJD(6) * t498 - t450 * t479 + t457 * t509;
	t463 = qJD(2) * t467;
	t491 = qJD(6) * t467 - t463 * t479 + t466 * t509;
	t442 = t459 * qJD(4) - t464 * t479;
	t441 = -t460 * qJD(4) + t464 * t478;
	t440 = t445 * qJD(4) + t479 * t500;
	t439 = -t446 * qJD(4) - t478 * t500;
	t438 = t443 * qJD(4) + t479 * t501;
	t437 = t496 * qJD(4) - t478 * t501;
	t1 = [0, t503 * t487 - t492 * t489, 0, t439 * t489 - t445 * t506, 0, -t440 * t487 - t450 * t489 + (-t446 * t489 + t457 * t487) * qJD(6); 0, t504 * t487 - t493 * t489, 0, t437 * t489 - t443 * t506, 0, -t438 * t487 - t447 * t489 + (t454 * t487 + t489 * t496) * qJD(6); 0, t502 * t487 + t491 * t489, 0, t441 * t489 - t459 * t506, 0, -t442 * t487 + t463 * t489 + (-t460 * t489 - t466 * t487) * qJD(6); 0, t492 * t487 + t503 * t489, 0, -t439 * t487 - t445 * t505, 0, -t440 * t489 + t450 * t487 + (t446 * t487 + t457 * t489) * qJD(6); 0, t493 * t487 + t504 * t489, 0, -t437 * t487 - t443 * t505, 0, -t438 * t489 + t447 * t487 + (t454 * t489 - t487 * t496) * qJD(6); 0, -t491 * t487 + t502 * t489, 0, -t441 * t487 - t459 * t505, 0, -t442 * t489 - t463 * t487 + (t460 * t487 - t466 * t489) * qJD(6); 0, t450 * t478 + t457 * t508, 0, t440, 0, 0; 0, t447 * t478 + t454 * t508, 0, t438, 0, 0; 0, -t463 * t478 - t466 * t508, 0, t442, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end