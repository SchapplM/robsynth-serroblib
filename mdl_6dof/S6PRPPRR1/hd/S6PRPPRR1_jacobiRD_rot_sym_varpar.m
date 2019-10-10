% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
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
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
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
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->10), mult. (96->24), div. (0->0), fcn. (96->10), ass. (0->18)
	t212 = sin(pkin(11));
	t216 = cos(pkin(11));
	t219 = sin(qJ(2));
	t220 = cos(qJ(2));
	t209 = (t212 * t219 - t216 * t220) * qJD(2);
	t224 = (t212 * t220 + t216 * t219) * qJD(2);
	t218 = cos(pkin(6));
	t217 = cos(pkin(10));
	t215 = cos(pkin(12));
	t214 = sin(pkin(6));
	t213 = sin(pkin(10));
	t211 = sin(pkin(12));
	t208 = t218 * t224;
	t207 = t218 * t209;
	t206 = t214 * t224;
	t205 = t213 * t208 + t217 * t209;
	t204 = -t217 * t208 + t213 * t209;
	t1 = [0, t205 * t215, 0, 0, 0, 0; 0, t204 * t215, 0, 0, 0, 0; 0, -t206 * t215, 0, 0, 0, 0; 0, -t205 * t211, 0, 0, 0, 0; 0, -t204 * t211, 0, 0, 0, 0; 0, t206 * t211, 0, 0, 0, 0; 0, t213 * t207 - t217 * t224, 0, 0, 0, 0; 0, -t217 * t207 - t213 * t224, 0, 0, 0, 0; 0, -t214 * t209, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:29
	% EndTime: 2019-10-09 21:24:29
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (111->32), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->39)
	t307 = sin(pkin(10));
	t308 = sin(pkin(6));
	t322 = t307 * t308;
	t310 = cos(pkin(10));
	t321 = t308 * t310;
	t312 = sin(qJ(2));
	t320 = qJD(2) * t312;
	t313 = cos(qJ(2));
	t319 = qJD(2) * t313;
	t305 = pkin(12) + qJ(5);
	t303 = sin(t305);
	t318 = qJD(5) * t303;
	t304 = cos(t305);
	t317 = qJD(5) * t304;
	t306 = sin(pkin(11));
	t309 = cos(pkin(11));
	t311 = cos(pkin(6));
	t290 = (t306 * t320 - t309 * t319) * t311;
	t297 = -t306 * t319 - t309 * t320;
	t281 = -t310 * t290 + t307 * t297;
	t283 = t307 * t290 + t310 * t297;
	t316 = t313 * t306 + t312 * t309;
	t315 = t312 * t306 - t313 * t309;
	t292 = t315 * t308;
	t314 = qJD(2) * t316;
	t296 = t315 * qJD(2);
	t295 = t316 * t311;
	t294 = t315 * t311;
	t293 = t316 * t308;
	t291 = t311 * t314;
	t289 = qJD(2) * t292;
	t288 = t308 * t314;
	t287 = -t307 * t295 - t310 * t315;
	t286 = t307 * t294 - t310 * t316;
	t285 = t310 * t295 - t307 * t315;
	t284 = -t310 * t294 - t307 * t316;
	t282 = t307 * t291 + t310 * t296;
	t280 = -t310 * t291 + t307 * t296;
	t1 = [0, t282 * t304 - t286 * t318, 0, 0, -t283 * t303 + (-t287 * t304 - t303 * t322) * qJD(5), 0; 0, t280 * t304 - t284 * t318, 0, 0, -t281 * t303 + (-t285 * t304 + t303 * t321) * qJD(5), 0; 0, -t288 * t304 + t292 * t318, 0, 0, t289 * t303 + (-t293 * t304 - t303 * t311) * qJD(5), 0; 0, -t282 * t303 - t286 * t317, 0, 0, -t283 * t304 + (t287 * t303 - t304 * t322) * qJD(5), 0; 0, -t280 * t303 - t284 * t317, 0, 0, -t281 * t304 + (t285 * t303 + t304 * t321) * qJD(5), 0; 0, t288 * t303 + t292 * t317, 0, 0, t289 * t304 + (t293 * t303 - t304 * t311) * qJD(5), 0; 0, t283, 0, 0, 0, 0; 0, t281, 0, 0, 0, 0; 0, -t289, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:31
	% EndTime: 2019-10-09 21:24:31
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (396->67), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->61)
	t479 = cos(pkin(6));
	t474 = sin(pkin(11));
	t477 = cos(pkin(11));
	t481 = sin(qJ(2));
	t483 = cos(qJ(2));
	t490 = t483 * t474 + t481 * t477;
	t461 = t490 * t479;
	t464 = t481 * t474 - t483 * t477;
	t462 = t464 * qJD(2);
	t475 = sin(pkin(10));
	t476 = sin(pkin(6));
	t508 = t475 * t476;
	t478 = cos(pkin(10));
	t507 = t476 * t478;
	t504 = qJD(2) * t481;
	t503 = qJD(2) * t483;
	t473 = pkin(12) + qJ(5);
	t471 = sin(t473);
	t502 = qJD(5) * t471;
	t472 = cos(t473);
	t501 = qJD(5) * t472;
	t500 = qJD(6) * t472;
	t480 = sin(qJ(6));
	t499 = qJD(6) * t480;
	t482 = cos(qJ(6));
	t498 = qJD(6) * t482;
	t488 = t464 * t479;
	t447 = -t475 * t490 - t478 * t488;
	t458 = (t474 * t504 - t477 * t503) * t479;
	t463 = -t474 * t503 - t477 * t504;
	t494 = -t478 * t458 + t475 * t463;
	t497 = -t447 * t500 + t494;
	t450 = t475 * t488 - t478 * t490;
	t493 = t475 * t458 + t478 * t463;
	t496 = -t450 * t500 + t493;
	t457 = t476 * t462;
	t459 = t464 * t476;
	t495 = t459 * t500 - t457;
	t460 = t490 * t476;
	t453 = t460 * t472 + t479 * t471;
	t452 = -t460 * t471 + t479 * t472;
	t492 = t478 * t461 - t475 * t464;
	t491 = -t475 * t461 - t478 * t464;
	t436 = -t471 * t492 - t472 * t507;
	t489 = t471 * t507 - t472 * t492;
	t438 = -t471 * t491 + t472 * t508;
	t439 = t471 * t508 + t472 * t491;
	t487 = qJD(2) * t461;
	t440 = t475 * t462 - t478 * t487;
	t486 = -qJD(6) * t492 - t440 * t472 + t447 * t502;
	t443 = t478 * t462 + t475 * t487;
	t485 = -qJD(6) * t491 - t443 * t472 + t450 * t502;
	t456 = qJD(2) * t460;
	t484 = qJD(6) * t460 - t456 * t472 + t459 * t502;
	t435 = t452 * qJD(5) - t457 * t472;
	t434 = -t453 * qJD(5) + t457 * t471;
	t433 = t438 * qJD(5) + t472 * t493;
	t432 = -t439 * qJD(5) - t471 * t493;
	t431 = t436 * qJD(5) + t472 * t494;
	t430 = t489 * qJD(5) - t471 * t494;
	t1 = [0, t496 * t480 - t485 * t482, 0, 0, t432 * t482 - t438 * t499, -t433 * t480 - t443 * t482 + (-t439 * t482 + t450 * t480) * qJD(6); 0, t497 * t480 - t486 * t482, 0, 0, t430 * t482 - t436 * t499, -t431 * t480 - t440 * t482 + (t447 * t480 + t482 * t489) * qJD(6); 0, t495 * t480 + t484 * t482, 0, 0, t434 * t482 - t452 * t499, -t435 * t480 + t456 * t482 + (-t453 * t482 - t459 * t480) * qJD(6); 0, t485 * t480 + t496 * t482, 0, 0, -t432 * t480 - t438 * t498, -t433 * t482 + t443 * t480 + (t439 * t480 + t450 * t482) * qJD(6); 0, t486 * t480 + t497 * t482, 0, 0, -t430 * t480 - t436 * t498, -t431 * t482 + t440 * t480 + (t447 * t482 - t480 * t489) * qJD(6); 0, -t484 * t480 + t495 * t482, 0, 0, -t434 * t480 - t452 * t498, -t435 * t482 - t456 * t480 + (t453 * t480 - t459 * t482) * qJD(6); 0, t443 * t471 + t450 * t501, 0, 0, t433, 0; 0, t440 * t471 + t447 * t501, 0, 0, t431, 0; 0, -t456 * t471 - t459 * t501, 0, 0, t435, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end