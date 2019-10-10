% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR14_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR14_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
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
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
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
	% StartTime: 2019-10-10 10:28:17
	% EndTime: 2019-10-10 10:28:17
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
	t280 = qJD(1) * t284 + qJD(2) * t283;
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
	t316 = qJD(4) * (t291 * t299 + t294 * t313) + t280 * t294 - t291 * t303;
	t315 = t289 * t293;
	t314 = t289 * t295;
	t307 = qJD(2) * t295;
	t306 = qJD(4) * t291;
	t305 = qJD(4) * t294;
	t302 = t296 * t308;
	t301 = t289 * qJD(2) * t292;
	t285 = -t290 * t312 + t309;
	t297 = -t294 * t303 - t280 * t291 + (-t291 * t313 + t294 * t299) * qJD(4);
	t281 = qJD(1) * t285 + qJD(2) * t299;
	t279 = -qJD(1) * t283 - qJD(2) * t284;
	t278 = -qJD(1) * t304 - t296 * t307 + (qJD(2) * t290 + qJD(1)) * t312;
	t277 = t294 * t302 - t278 * t291 + (t284 * t294 - t291 * t315) * qJD(4);
	t276 = -t291 * t302 - t278 * t294 + (-t284 * t291 - t294 * t315) * qJD(4);
	t1 = [t297, t279 * t291 + t285 * t305, 0, t276, 0, 0; t277, t281 * t291 + t283 * t305, 0, t316, 0, 0; 0, (t291 * t307 + t292 * t305) * t289, 0, t294 * t301 + (-t290 * t294 + t291 * t314) * qJD(4), 0, 0; -t316, t279 * t294 - t285 * t306, 0, -t277, 0, 0; t276, t281 * t294 - t283 * t306, 0, t297, 0, 0; 0, (-t292 * t306 + t294 * t307) * t289, 0, -t291 * t301 + (t290 * t291 + t294 * t314) * qJD(4), 0, 0; -t281, t278, 0, 0, 0, 0; t279, -t280, 0, 0, 0, 0; 0, -t301, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:18
	% EndTime: 2019-10-10 10:28:18
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (95->37), mult. (310->69), div. (0->0), fcn. (322->8), ass. (0->36)
	t350 = cos(pkin(6));
	t352 = sin(qJ(2));
	t356 = cos(qJ(1));
	t370 = t356 * t352;
	t353 = sin(qJ(1));
	t355 = cos(qJ(2));
	t371 = t353 * t355;
	t341 = t350 * t370 + t371;
	t342 = t350 * t371 + t370;
	t338 = qJD(1) * t342 + qJD(2) * t341;
	t351 = sin(qJ(4));
	t354 = cos(qJ(4));
	t369 = t356 * t355;
	t364 = t350 * t369;
	t372 = t353 * t352;
	t359 = t364 - t372;
	t349 = sin(pkin(6));
	t368 = qJD(1) * t349;
	t363 = t353 * t368;
	t373 = t349 * t356;
	t376 = qJD(4) * (t351 * t359 + t354 * t373) + t338 * t354 - t351 * t363;
	t375 = t349 * t353;
	t374 = t349 * t355;
	t367 = qJD(2) * t355;
	t366 = qJD(4) * t351;
	t365 = qJD(4) * t354;
	t362 = t356 * t368;
	t361 = t349 * qJD(2) * t352;
	t343 = -t350 * t372 + t369;
	t357 = t338 * t351 + t354 * t363 - t359 * t365 + t366 * t373;
	t339 = qJD(1) * t343 + qJD(2) * t359;
	t337 = -qJD(1) * t341 - qJD(2) * t342;
	t336 = -qJD(1) * t364 - t356 * t367 + (qJD(2) * t350 + qJD(1)) * t372;
	t335 = t354 * t362 - t336 * t351 + (t342 * t354 - t351 * t375) * qJD(4);
	t334 = t351 * t362 + t336 * t354 + (t342 * t351 + t354 * t375) * qJD(4);
	t1 = [-t339, t336, 0, 0, 0, 0; t337, -t338, 0, 0, 0, 0; 0, -t361, 0, 0, 0, 0; t357, -t337 * t351 - t343 * t365, 0, t334, 0, 0; -t335, -t339 * t351 - t341 * t365, 0, -t376, 0, 0; 0, (-t351 * t367 - t352 * t365) * t349, 0, -t354 * t361 + (t350 * t354 - t351 * t374) * qJD(4), 0, 0; t376, -t337 * t354 + t343 * t366, 0, t335, 0, 0; t334, -t339 * t354 + t341 * t366, 0, t357, 0, 0; 0, (t352 * t366 - t354 * t367) * t349, 0, t351 * t361 + (-t350 * t351 - t354 * t374) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:19
	% EndTime: 2019-10-10 10:28:19
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (289->77), mult. (898->141), div. (0->0), fcn. (964->10), ass. (0->65)
	t463 = cos(pkin(6));
	t466 = sin(qJ(2));
	t471 = cos(qJ(1));
	t498 = t471 * t466;
	t467 = sin(qJ(1));
	t470 = cos(qJ(2));
	t499 = t467 * t470;
	t453 = t463 * t498 + t499;
	t454 = t463 * t499 + t498;
	t441 = t454 * qJD(1) + t453 * qJD(2);
	t497 = t471 * t470;
	t486 = t463 * t497;
	t500 = t467 * t466;
	t452 = -t486 + t500;
	t465 = sin(qJ(4));
	t469 = cos(qJ(4));
	t462 = sin(pkin(6));
	t496 = qJD(1) * t462;
	t485 = t467 * t496;
	t502 = t462 * t471;
	t487 = t469 * t502;
	t434 = (qJD(4) * t452 + t485) * t465 - qJD(4) * t487 - t441 * t469;
	t480 = qJD(2) * t463 + qJD(1);
	t488 = t463 * t500;
	t495 = qJD(2) * t466;
	t442 = -qJD(1) * t488 - t467 * t495 + t480 * t497;
	t447 = t452 * t469 + t465 * t502;
	t464 = sin(qJ(6));
	t468 = cos(qJ(6));
	t512 = -t434 * t464 + (t447 * t468 + t453 * t464) * qJD(6) - t442 * t468;
	t511 = (t447 * t464 - t453 * t468) * qJD(6) + t434 * t468 - t442 * t464;
	t436 = t447 * qJD(4) + t441 * t465 + t469 * t485;
	t504 = t462 * t467;
	t503 = t462 * t470;
	t501 = t466 * t468;
	t494 = qJD(2) * t470;
	t493 = qJD(4) * t465;
	t492 = qJD(4) * t469;
	t491 = qJD(6) * t464;
	t490 = qJD(6) * t468;
	t489 = qJD(6) * t469;
	t484 = t471 * t496;
	t483 = t462 * t495;
	t482 = t462 * t494;
	t479 = qJD(2) + t489;
	t439 = -qJD(1) * t486 - t471 * t494 + t480 * t500;
	t455 = -t488 + t497;
	t478 = t455 * t489 - t439;
	t477 = t453 * t489 + t441;
	t476 = (-qJD(2) * t469 - qJD(6)) * t470;
	t446 = t454 * t465 + t469 * t504;
	t445 = -t454 * t469 + t465 * t504;
	t450 = t463 * t465 + t469 * t503;
	t451 = t463 * t469 - t465 * t503;
	t440 = -t453 * qJD(1) - t454 * qJD(2);
	t473 = qJD(6) * t454 - t440 * t469 + t455 * t493;
	t472 = qJD(6) * t452 - t442 * t469 + t453 * t493;
	t449 = t452 * t465 - t487;
	t444 = t451 * qJD(4) - t469 * t483;
	t443 = -t450 * qJD(4) + t465 * t483;
	t438 = -t445 * qJD(4) - t439 * t465 + t469 * t484;
	t437 = t446 * qJD(4) + t439 * t469 + t465 * t484;
	t433 = t437 * t464 + t440 * t468 + (t445 * t468 - t455 * t464) * qJD(6);
	t432 = t437 * t468 - t440 * t464 + (-t445 * t464 - t455 * t468) * qJD(6);
	t1 = [t512, t473 * t464 - t478 * t468, 0, t438 * t464 + t446 * t490, 0, t432; t433, t472 * t464 - t477 * t468, 0, t436 * t464 + t449 * t490, 0, t511; 0, (-t479 * t501 + (t466 * t493 + t476) * t464) * t462, 0, t443 * t464 + t451 * t490, 0, -t464 * t482 + t444 * t468 + (-t450 * t464 - t462 * t501) * qJD(6); -t511, t478 * t464 + t473 * t468, 0, t438 * t468 - t446 * t491, 0, -t433; t432, t477 * t464 + t472 * t468, 0, t436 * t468 - t449 * t491, 0, t512; 0, (t468 * t476 + (t479 * t464 + t468 * t493) * t466) * t462, 0, t443 * t468 - t451 * t491, 0, -t468 * t482 - t444 * t464 + (t462 * t464 * t466 - t450 * t468) * qJD(6); -t436, t440 * t465 + t455 * t492, 0, -t437, 0, 0; t438, t442 * t465 + t453 * t492, 0, -t434, 0, 0; 0, (t465 * t494 + t466 * t492) * t462, 0, -t444, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end