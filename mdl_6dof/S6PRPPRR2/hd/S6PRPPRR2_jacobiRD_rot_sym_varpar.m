% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPPRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->7), mult. (60->18), div. (0->0), fcn. (60->8), ass. (0->13)
	t180 = sin(pkin(11));
	t183 = cos(pkin(11));
	t186 = sin(qJ(2));
	t187 = cos(qJ(2));
	t178 = (t180 * t186 - t183 * t187) * qJD(2);
	t179 = (t180 * t187 + t183 * t186) * qJD(2);
	t185 = cos(pkin(6));
	t184 = cos(pkin(10));
	t182 = sin(pkin(6));
	t181 = sin(pkin(10));
	t177 = t185 * t179;
	t176 = t185 * t178;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -t181 * t177 - t184 * t178, 0, 0, 0, 0; 0, t184 * t177 - t181 * t178, 0, 0, 0, 0; 0, t182 * t179, 0, 0, 0, 0; 0, t181 * t176 - t184 * t179, 0, 0, 0, 0; 0, -t184 * t176 - t181 * t179, 0, 0, 0, 0; 0, -t182 * t178, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:20
	% EndTime: 2019-10-09 21:26:20
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (81->29), mult. (282->70), div. (0->0), fcn. (310->10), ass. (0->35)
	t281 = sin(pkin(11));
	t284 = cos(pkin(11));
	t288 = sin(qJ(2));
	t290 = cos(qJ(2));
	t279 = t288 * t281 - t290 * t284;
	t277 = t279 * qJD(2);
	t283 = sin(pkin(6));
	t287 = sin(qJ(5));
	t296 = t283 * t287;
	t289 = cos(qJ(5));
	t295 = t283 * t289;
	t294 = qJD(5) * t287;
	t293 = qJD(5) * t289;
	t292 = t290 * t281 + t288 * t284;
	t286 = cos(pkin(6));
	t276 = t292 * t286;
	t278 = t292 * qJD(2);
	t291 = qJD(2) * t276;
	t285 = cos(pkin(10));
	t282 = sin(pkin(10));
	t275 = t279 * t286;
	t274 = t292 * t283;
	t273 = t279 * t283;
	t272 = t286 * t277;
	t271 = t283 * t277;
	t270 = t283 * t278;
	t268 = -t282 * t275 + t285 * t292;
	t267 = -t282 * t276 - t285 * t279;
	t266 = t285 * t275 + t282 * t292;
	t265 = t285 * t276 - t282 * t279;
	t264 = t282 * t272 - t285 * t278;
	t263 = t285 * t277 + t282 * t291;
	t262 = -t285 * t272 - t282 * t278;
	t261 = t282 * t277 - t285 * t291;
	t1 = [0, t264 * t287 + t267 * t293, 0, 0, -t263 * t289 + (-t268 * t287 - t282 * t295) * qJD(5), 0; 0, t262 * t287 + t265 * t293, 0, 0, -t261 * t289 + (-t266 * t287 + t285 * t295) * qJD(5), 0; 0, -t271 * t287 + t274 * t293, 0, 0, t270 * t289 + (-t273 * t287 - t286 * t289) * qJD(5), 0; 0, t264 * t289 - t267 * t294, 0, 0, t263 * t287 + (-t268 * t289 + t282 * t296) * qJD(5), 0; 0, t262 * t289 - t265 * t294, 0, 0, t261 * t287 + (-t266 * t289 - t285 * t296) * qJD(5), 0; 0, -t271 * t289 - t274 * t294, 0, 0, -t270 * t287 + (-t273 * t289 + t286 * t287) * qJD(5), 0; 0, t263, 0, 0, 0, 0; 0, t261, 0, 0, 0, 0; 0, -t270, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:22
	% EndTime: 2019-10-09 21:26:22
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (303->66), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->60)
	t474 = cos(pkin(6));
	t469 = sin(pkin(11));
	t472 = cos(pkin(11));
	t477 = sin(qJ(2));
	t480 = cos(qJ(2));
	t487 = t480 * t469 + t477 * t472;
	t459 = t487 * t474;
	t462 = t477 * t469 - t480 * t472;
	t460 = t462 * qJD(2);
	t471 = sin(pkin(6));
	t476 = sin(qJ(5));
	t505 = t471 * t476;
	t479 = cos(qJ(5));
	t504 = t471 * t479;
	t501 = qJD(2) * t477;
	t500 = qJD(2) * t480;
	t499 = qJD(5) * t476;
	t498 = qJD(5) * t479;
	t475 = sin(qJ(6));
	t497 = qJD(6) * t475;
	t496 = qJD(6) * t476;
	t478 = cos(qJ(6));
	t495 = qJD(6) * t478;
	t470 = sin(pkin(10));
	t473 = cos(pkin(10));
	t484 = qJD(2) * t459;
	t438 = t470 * t460 - t473 * t484;
	t489 = t473 * t459 - t470 * t462;
	t494 = -t489 * t496 + t438;
	t441 = t473 * t460 + t470 * t484;
	t488 = -t470 * t459 - t473 * t462;
	t493 = -t488 * t496 + t441;
	t458 = t487 * t471;
	t454 = qJD(2) * t458;
	t492 = -t458 * t496 - t454;
	t456 = (t469 * t501 - t472 * t500) * t474;
	t461 = -t469 * t500 - t472 * t501;
	t491 = -t473 * t456 + t470 * t461;
	t490 = t470 * t456 + t473 * t461;
	t457 = t462 * t471;
	t450 = t457 * t479 - t474 * t476;
	t451 = t457 * t476 + t474 * t479;
	t485 = t462 * t474;
	t445 = -t470 * t487 - t473 * t485;
	t486 = t445 * t476 + t473 * t504;
	t436 = -t445 * t479 + t473 * t505;
	t448 = t470 * t485 - t473 * t487;
	t435 = -t448 * t476 + t470 * t504;
	t434 = -t448 * t479 - t470 * t505;
	t483 = qJD(6) * t445 + t476 * t491 + t489 * t498;
	t482 = qJD(6) * t448 + t476 * t490 + t488 * t498;
	t455 = t471 * t460;
	t481 = -qJD(6) * t457 - t455 * t476 + t458 * t498;
	t433 = -t451 * qJD(5) + t454 * t479;
	t432 = t450 * qJD(5) + t454 * t476;
	t431 = t486 * qJD(5) - t438 * t479;
	t430 = t436 * qJD(5) - t438 * t476;
	t429 = -t435 * qJD(5) - t441 * t479;
	t428 = t434 * qJD(5) - t441 * t476;
	t1 = [0, t493 * t475 + t482 * t478, 0, 0, t429 * t478 - t434 * t497, -t428 * t475 + t490 * t478 + (-t435 * t478 - t475 * t488) * qJD(6); 0, t494 * t475 + t483 * t478, 0, 0, t431 * t478 - t436 * t497, -t430 * t475 + t491 * t478 + (-t475 * t489 + t478 * t486) * qJD(6); 0, t492 * t475 + t481 * t478, 0, 0, t433 * t478 - t450 * t497, -t432 * t475 - t455 * t478 + (-t451 * t478 - t458 * t475) * qJD(6); 0, -t482 * t475 + t493 * t478, 0, 0, -t429 * t475 - t434 * t495, -t428 * t478 - t490 * t475 + (t435 * t475 - t478 * t488) * qJD(6); 0, -t483 * t475 + t494 * t478, 0, 0, -t431 * t475 - t436 * t495, -t430 * t478 - t491 * t475 + (-t475 * t486 - t478 * t489) * qJD(6); 0, -t481 * t475 + t492 * t478, 0, 0, -t433 * t475 - t450 * t495, -t432 * t478 + t455 * t475 + (t451 * t475 - t458 * t478) * qJD(6); 0, -t479 * t490 + t488 * t499, 0, 0, t428, 0; 0, -t479 * t491 + t489 * t499, 0, 0, t430, 0; 0, t455 * t479 + t458 * t499, 0, 0, t432, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end