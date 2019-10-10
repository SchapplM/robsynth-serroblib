% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:15
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:37
	% EndTime: 2019-10-09 23:15:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:37
	% EndTime: 2019-10-09 23:15:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:37
	% EndTime: 2019-10-09 23:15:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(12));
	t58 = sin(pkin(12));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:38
	% EndTime: 2019-10-09 23:15:38
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(6));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(12));
	t212 = cos(pkin(12));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0, 0; 0, t204, 0, 0, 0, 0; 0, t202, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:38
	% EndTime: 2019-10-09 23:15:38
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (134->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t274 = qJ(3) + qJ(4);
	t271 = sin(t274);
	t273 = qJD(3) + qJD(4);
	t293 = t271 * t273;
	t272 = cos(t274);
	t292 = t272 * t273;
	t276 = sin(pkin(6));
	t291 = t273 * t276;
	t280 = cos(qJ(2));
	t290 = t273 * t280;
	t278 = cos(pkin(6));
	t279 = sin(qJ(2));
	t289 = t278 * t279;
	t288 = t278 * t280;
	t287 = qJD(2) * t279;
	t286 = t279 * t291;
	t285 = t276 * qJD(2) * t280;
	t275 = sin(pkin(12));
	t277 = cos(pkin(12));
	t267 = -t275 * t279 + t277 * t288;
	t263 = t267 * qJD(2);
	t284 = t277 * t291 - t263;
	t269 = -t275 * t288 - t277 * t279;
	t265 = t269 * qJD(2);
	t283 = -t275 * t291 - t265;
	t268 = t275 * t280 + t277 * t289;
	t282 = t275 * t289 - t277 * t280;
	t281 = -t273 * t278 - t285;
	t266 = t282 * qJD(2);
	t264 = t268 * qJD(2);
	t262 = t271 * t286 + t281 * t272;
	t261 = t281 * t271 - t272 * t286;
	t260 = t283 * t272 - t282 * t293;
	t259 = t283 * t271 + t282 * t292;
	t258 = t268 * t293 + t284 * t272;
	t257 = -t268 * t292 + t284 * t271;
	t1 = [0, t266 * t272 - t269 * t293, t259, t259, 0, 0; 0, -t264 * t272 - t267 * t293, t257, t257, 0, 0; 0, (-t271 * t290 - t272 * t287) * t276, t261, t261, 0, 0; 0, -t266 * t271 - t269 * t292, t260, t260, 0, 0; 0, t264 * t271 - t267 * t292, t258, t258, 0, 0; 0, (t271 * t287 - t272 * t290) * t276, t262, t262, 0, 0; 0, t265, 0, 0, 0, 0; 0, t263, 0, 0, 0, 0; 0, t285, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:40
	% EndTime: 2019-10-09 23:15:40
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (388->64), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->66)
	t444 = qJ(3) + qJ(4);
	t441 = sin(t444);
	t442 = cos(t444);
	t450 = sin(qJ(2));
	t443 = qJD(3) + qJD(4);
	t452 = cos(qJ(2));
	t479 = t443 * t452;
	t482 = (qJD(2) * t442 - qJD(5)) * t450 + t441 * t479;
	t481 = t441 * t443;
	t480 = t442 * t443;
	t445 = sin(pkin(12));
	t446 = sin(pkin(6));
	t478 = t445 * t446;
	t447 = cos(pkin(12));
	t477 = t446 * t447;
	t476 = t446 * t450;
	t475 = t446 * t452;
	t448 = cos(pkin(6));
	t474 = t448 * t450;
	t473 = t448 * t452;
	t472 = qJD(2) * t450;
	t471 = qJD(2) * t452;
	t470 = qJD(5) * t442;
	t449 = sin(qJ(5));
	t469 = qJD(5) * t449;
	t451 = cos(qJ(5));
	t468 = qJD(5) * t451;
	t466 = t445 * t474;
	t465 = t441 * t476;
	t464 = t442 * t476;
	t463 = t446 * t472;
	t456 = -t445 * t450 + t447 * t473;
	t431 = t456 * qJD(2);
	t461 = t443 * t477 - t431;
	t437 = t445 * t473 + t447 * t450;
	t433 = t437 * qJD(2);
	t460 = t443 * t478 - t433;
	t459 = -t456 * t470 + t431;
	t458 = t437 * t470 - t433;
	t457 = (qJD(2) - t470) * t452;
	t436 = t445 * t452 + t447 * t474;
	t455 = t443 * t448 + t446 * t471;
	t432 = t436 * qJD(2);
	t454 = qJD(5) * t436 - t432 * t442 - t456 * t481;
	t434 = -qJD(2) * t466 + t447 * t471;
	t438 = t447 * t452 - t466;
	t453 = qJD(5) * t438 - t434 * t442 + t437 * t481;
	t430 = t448 * t441 + t464;
	t429 = t448 * t442 - t465;
	t428 = t438 * t442 + t441 * t478;
	t427 = -t438 * t441 + t442 * t478;
	t426 = t436 * t442 - t441 * t477;
	t425 = -t436 * t441 - t442 * t477;
	t424 = t455 * t442 - t443 * t465;
	t423 = -t455 * t441 - t443 * t464;
	t422 = -t438 * t481 + t460 * t442;
	t421 = -t438 * t480 - t460 * t441;
	t420 = -t436 * t481 - t461 * t442;
	t419 = -t436 * t480 + t461 * t441;
	t418 = t423 * t451 - t429 * t469;
	t417 = -t423 * t449 - t429 * t468;
	t416 = t421 * t451 - t427 * t469;
	t415 = -t421 * t449 - t427 * t468;
	t414 = t419 * t451 - t425 * t469;
	t413 = -t419 * t449 - t425 * t468;
	t1 = [0, t458 * t449 + t453 * t451, t416, t416, -t422 * t449 + t434 * t451 + (-t428 * t451 - t437 * t449) * qJD(5), 0; 0, t459 * t449 + t454 * t451, t414, t414, -t420 * t449 + t432 * t451 + (-t426 * t451 + t449 * t456) * qJD(5), 0; 0, (t449 * t457 - t482 * t451) * t446, t418, t418, t451 * t463 - t424 * t449 + (-t430 * t451 + t449 * t475) * qJD(5), 0; 0, -t453 * t449 + t458 * t451, t415, t415, -t422 * t451 - t434 * t449 + (t428 * t449 - t437 * t451) * qJD(5), 0; 0, -t454 * t449 + t459 * t451, t413, t413, -t420 * t451 - t432 * t449 + (t426 * t449 + t451 * t456) * qJD(5), 0; 0, (t482 * t449 + t451 * t457) * t446, t417, t417, -t449 * t463 - t424 * t451 + (t430 * t449 + t451 * t475) * qJD(5), 0; 0, -t434 * t441 - t437 * t480, t422, t422, 0, 0; 0, -t432 * t441 + t456 * t480, t420, t420, 0, 0; 0, (-t441 * t472 + t442 * t479) * t446, t424, t424, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:15:40
	% EndTime: 2019-10-09 23:15:40
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (636->60), mult. (848->119), div. (0->0), fcn. (930->10), ass. (0->75)
	t507 = qJ(3) + qJ(4);
	t501 = sin(t507);
	t503 = cos(t507);
	t504 = qJD(5) + qJD(6);
	t512 = sin(qJ(2));
	t505 = qJD(3) + qJD(4);
	t513 = cos(qJ(2));
	t541 = t505 * t513;
	t547 = (qJD(2) * t503 - t504) * t512 + t501 * t541;
	t506 = qJ(5) + qJ(6);
	t500 = sin(t506);
	t546 = t500 * t504;
	t545 = t501 * t505;
	t502 = cos(t506);
	t544 = t502 * t504;
	t543 = t503 * t504;
	t542 = t503 * t505;
	t508 = sin(pkin(12));
	t509 = sin(pkin(6));
	t540 = t508 * t509;
	t510 = cos(pkin(12));
	t539 = t509 * t510;
	t538 = t509 * t512;
	t511 = cos(pkin(6));
	t537 = t511 * t512;
	t536 = t511 * t513;
	t535 = qJD(2) * t512;
	t534 = qJD(2) * t513;
	t532 = t508 * t537;
	t531 = t501 * t538;
	t530 = t503 * t538;
	t495 = t508 * t513 + t510 * t537;
	t518 = -t508 * t512 + t510 * t536;
	t490 = t518 * qJD(2);
	t522 = t505 * t539 - t490;
	t479 = -t495 * t545 - t522 * t503;
	t529 = t504 * t518 - t479;
	t497 = t510 * t513 - t532;
	t496 = t508 * t536 + t510 * t512;
	t492 = t496 * qJD(2);
	t520 = t505 * t540 - t492;
	t481 = -t497 * t545 + t520 * t503;
	t528 = -t496 * t504 - t481;
	t491 = t495 * qJD(2);
	t527 = (t495 * t503 - t501 * t539) * t504 - t491;
	t493 = -qJD(2) * t532 + t510 * t534;
	t526 = (t497 * t503 + t501 * t540) * t504 - t493;
	t516 = t505 * t511 + t509 * t534;
	t483 = t516 * t503 - t505 * t531;
	t524 = t504 * t509 * t513 - t483;
	t523 = -t518 * t543 + t490;
	t521 = t496 * t543 - t492;
	t519 = (qJD(2) - t543) * t513;
	t517 = -(t511 * t501 + t530) * t504 + t509 * t535;
	t515 = -t491 * t503 + t495 * t504 - t518 * t545;
	t514 = -t493 * t503 + t496 * t545 + t497 * t504;
	t488 = t511 * t503 - t531;
	t486 = -t497 * t501 + t503 * t540;
	t484 = -t495 * t501 - t503 * t539;
	t482 = -t516 * t501 - t505 * t530;
	t480 = -t497 * t542 - t520 * t501;
	t478 = -t495 * t542 + t522 * t501;
	t477 = t482 * t502 - t488 * t546;
	t476 = -t482 * t500 - t488 * t544;
	t475 = -t517 * t500 + t524 * t502;
	t474 = t524 * t500 + t517 * t502;
	t473 = t480 * t502 - t486 * t546;
	t472 = -t480 * t500 - t486 * t544;
	t471 = t478 * t502 - t484 * t546;
	t470 = -t478 * t500 - t484 * t544;
	t469 = t526 * t500 + t528 * t502;
	t468 = t528 * t500 - t526 * t502;
	t467 = t527 * t500 + t529 * t502;
	t466 = t529 * t500 - t527 * t502;
	t1 = [0, t521 * t500 + t514 * t502, t473, t473, t468, t468; 0, t523 * t500 + t515 * t502, t471, t471, t466, t466; 0, (t500 * t519 - t547 * t502) * t509, t477, t477, t474, t474; 0, -t514 * t500 + t521 * t502, t472, t472, t469, t469; 0, -t515 * t500 + t523 * t502, t470, t470, t467, t467; 0, (t547 * t500 + t502 * t519) * t509, t476, t476, t475, t475; 0, -t493 * t501 - t496 * t542, t481, t481, 0, 0; 0, -t491 * t501 + t518 * t542, t479, t479, 0, 0; 0, (-t501 * t535 + t503 * t541) * t509, t483, t483, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end