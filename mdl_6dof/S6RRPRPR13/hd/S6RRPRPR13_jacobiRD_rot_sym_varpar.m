% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
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
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:25
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
	% StartTime: 2019-10-10 10:26:25
	% EndTime: 2019-10-10 10:26:25
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
	% StartTime: 2019-10-10 10:26:25
	% EndTime: 2019-10-10 10:26:26
	% DurationCPUTime: 0.23s
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
	% StartTime: 2019-10-10 10:26:27
	% EndTime: 2019-10-10 10:26:27
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (160->46), mult. (518->100), div. (0->0), fcn. (536->10), ass. (0->43)
	t394 = cos(pkin(6));
	t396 = sin(qJ(2));
	t400 = cos(qJ(1));
	t416 = t400 * t396;
	t397 = sin(qJ(1));
	t399 = cos(qJ(2));
	t418 = t397 * t399;
	t384 = t394 * t416 + t418;
	t385 = t394 * t418 + t416;
	t380 = qJD(1) * t385 + qJD(2) * t384;
	t415 = t400 * t399;
	t419 = t397 * t396;
	t383 = -t394 * t415 + t419;
	t395 = sin(qJ(4));
	t398 = cos(qJ(4));
	t392 = sin(pkin(6));
	t414 = qJD(1) * t392;
	t409 = t397 * t414;
	t421 = t392 * t400;
	t423 = (t383 * t398 + t395 * t421) * qJD(4) + t380 * t395 + t398 * t409;
	t422 = t392 * t397;
	t420 = t395 * t399;
	t417 = t398 * t399;
	t413 = qJD(2) * t396;
	t412 = qJD(4) * t395;
	t411 = qJD(4) * t398;
	t410 = t394 * t419;
	t408 = t400 * t414;
	t407 = t392 * t413;
	t406 = t396 * t411;
	t404 = t410 - t415;
	t379 = -qJD(1) * t384 - qJD(2) * t385;
	t403 = t379 * t395 - t404 * t411;
	t381 = -qJD(1) * t410 - t397 * t413 + (qJD(2) * t394 + qJD(1)) * t415;
	t402 = t381 * t395 + t384 * t411;
	t373 = t380 * t398 + t411 * t421 + (-qJD(4) * t383 - t409) * t395;
	t393 = cos(pkin(11));
	t391 = sin(pkin(11));
	t382 = t398 * t407 + (t392 * t420 - t394 * t398) * qJD(4);
	t378 = qJD(1) * t383 + qJD(2) * t404;
	t376 = t398 * t408 - t378 * t395 + (t385 * t398 - t395 * t422) * qJD(4);
	t375 = t395 * t408 + t378 * t398 + (t385 * t395 + t398 * t422) * qJD(4);
	t1 = [-t381 * t391 - t393 * t423, t378 * t391 + t393 * t403, 0, -t375 * t393, 0, 0; t376 * t393 + t379 * t391, -t380 * t391 + t393 * t402, 0, t373 * t393, 0, 0; 0, (t393 * t406 + (-t391 * t396 + t393 * t420) * qJD(2)) * t392, 0, t382 * t393, 0, 0; -t381 * t393 + t391 * t423, t378 * t393 - t391 * t403, 0, t375 * t391, 0, 0; -t376 * t391 + t379 * t393, -t380 * t393 - t391 * t402, 0, -t373 * t391, 0, 0; 0, (-t391 * t406 + (-t391 * t420 - t393 * t396) * qJD(2)) * t392, 0, -t382 * t391, 0, 0; t373, -t379 * t398 - t404 * t412, 0, t376, 0, 0; t375, -t381 * t398 + t384 * t412, 0, t423, 0, 0; 0, (-qJD(2) * t417 + t396 * t412) * t392, 0, t395 * t407 + (-t392 * t417 - t394 * t395) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:28
	% EndTime: 2019-10-10 10:26:28
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (371->72), mult. (898->138), div. (0->0), fcn. (964->10), ass. (0->64)
	t477 = cos(pkin(6));
	t479 = sin(qJ(2));
	t483 = cos(qJ(1));
	t509 = t483 * t479;
	t480 = sin(qJ(1));
	t482 = cos(qJ(2));
	t510 = t480 * t482;
	t466 = t477 * t509 + t510;
	t467 = t477 * t510 + t509;
	t454 = t467 * qJD(1) + t466 * qJD(2);
	t508 = t483 * t482;
	t511 = t480 * t479;
	t465 = -t477 * t508 + t511;
	t478 = sin(qJ(4));
	t481 = cos(qJ(4));
	t476 = sin(pkin(6));
	t512 = t476 * t483;
	t460 = t465 * t481 + t478 * t512;
	t507 = qJD(1) * t476;
	t497 = t480 * t507;
	t447 = t460 * qJD(4) + t454 * t478 + t481 * t497;
	t499 = t477 * t511;
	t506 = qJD(2) * t479;
	t455 = -qJD(1) * t499 - t480 * t506 + (qJD(2) * t477 + qJD(1)) * t508;
	t498 = t481 * t512;
	t461 = -t465 * t478 + t498;
	t475 = pkin(11) + qJ(6);
	t473 = sin(t475);
	t474 = cos(t475);
	t524 = -t447 * t474 + (-t461 * t473 - t466 * t474) * qJD(6) - t455 * t473;
	t523 = (t461 * t474 - t466 * t473) * qJD(6) - t447 * t473 + t455 * t474;
	t503 = qJD(4) * t481;
	t520 = (qJD(2) * t478 + qJD(6)) * t482 + t479 * t503;
	t515 = t476 * t479;
	t514 = t476 * t480;
	t513 = t476 * t482;
	t505 = qJD(2) * t482;
	t504 = qJD(4) * t478;
	t502 = qJD(6) * t473;
	t501 = qJD(6) * t474;
	t500 = qJD(6) * t478;
	t496 = t483 * t507;
	t495 = t476 * t506;
	t494 = t476 * t505;
	t487 = t499 - t508;
	t452 = t465 * qJD(1) + t487 * qJD(2);
	t491 = t487 * t500 + t452;
	t490 = -t466 * t500 - t454;
	t489 = (-qJD(2) - t500) * t479;
	t459 = t467 * t478 + t481 * t514;
	t458 = t467 * t481 - t478 * t514;
	t463 = -t477 * t478 - t481 * t513;
	t488 = -t477 * t481 + t478 * t513;
	t453 = -t466 * qJD(1) - t467 * qJD(2);
	t485 = -qJD(6) * t467 + t453 * t478 - t487 * t503;
	t484 = -qJD(6) * t465 + t455 * t478 + t466 * t503;
	t446 = t454 * t481 + qJD(4) * t498 + (-qJD(4) * t465 - t497) * t478;
	t457 = t488 * qJD(4) + t481 * t495;
	t456 = t463 * qJD(4) + t478 * t495;
	t450 = t458 * qJD(4) - t452 * t478 + t481 * t496;
	t449 = t459 * qJD(4) + t452 * t481 + t478 * t496;
	t445 = t450 * t474 + t453 * t473 + (-t459 * t473 - t474 * t487) * qJD(6);
	t444 = -t450 * t473 + t453 * t474 + (-t459 * t474 + t473 * t487) * qJD(6);
	t1 = [t524, t491 * t473 + t485 * t474, 0, -t449 * t474 - t458 * t502, 0, t444; t445, t490 * t473 + t484 * t474, 0, t446 * t474 - t460 * t502, 0, t523; 0, (t473 * t489 + t520 * t474) * t476, 0, t457 * t474 - t463 * t502, 0, t474 * t494 - t456 * t473 + (-t473 * t515 + t474 * t488) * qJD(6); -t523, -t485 * t473 + t491 * t474, 0, t449 * t473 - t458 * t501, 0, -t445; t444, -t484 * t473 + t490 * t474, 0, -t446 * t473 - t460 * t501, 0, t524; 0, (-t520 * t473 + t474 * t489) * t476, 0, -t457 * t473 - t463 * t501, 0, -t473 * t494 - t456 * t474 + (-t473 * t488 - t474 * t515) * qJD(6); t446, -t453 * t481 - t487 * t504, 0, t450, 0, 0; t449, -t455 * t481 + t466 * t504, 0, t447, 0, 0; 0, (t479 * t504 - t481 * t505) * t476, 0, t456, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end