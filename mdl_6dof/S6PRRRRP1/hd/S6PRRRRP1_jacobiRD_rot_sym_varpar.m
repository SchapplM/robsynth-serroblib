% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRRP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:06
	% EndTime: 2019-10-09 23:02:06
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
	t210 = sin(pkin(11));
	t212 = cos(pkin(11));
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
	% StartTime: 2019-10-09 23:02:06
	% EndTime: 2019-10-09 23:02:06
	% DurationCPUTime: 0.12s
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
	t275 = sin(pkin(11));
	t277 = cos(pkin(11));
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
	% StartTime: 2019-10-09 23:02:08
	% EndTime: 2019-10-09 23:02:08
	% DurationCPUTime: 0.37s
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
	t445 = sin(pkin(11));
	t446 = sin(pkin(6));
	t478 = t445 * t446;
	t447 = cos(pkin(11));
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
	% StartTime: 2019-10-09 23:02:08
	% EndTime: 2019-10-09 23:02:09
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (388->64), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->66)
	t462 = qJ(3) + qJ(4);
	t459 = sin(t462);
	t460 = cos(t462);
	t468 = sin(qJ(2));
	t461 = qJD(3) + qJD(4);
	t470 = cos(qJ(2));
	t497 = t461 * t470;
	t500 = (qJD(2) * t460 - qJD(5)) * t468 + t459 * t497;
	t499 = t459 * t461;
	t498 = t460 * t461;
	t463 = sin(pkin(11));
	t464 = sin(pkin(6));
	t496 = t463 * t464;
	t465 = cos(pkin(11));
	t495 = t464 * t465;
	t494 = t464 * t468;
	t493 = t464 * t470;
	t466 = cos(pkin(6));
	t492 = t466 * t468;
	t491 = t466 * t470;
	t490 = qJD(2) * t468;
	t489 = qJD(2) * t470;
	t488 = qJD(5) * t460;
	t467 = sin(qJ(5));
	t487 = qJD(5) * t467;
	t469 = cos(qJ(5));
	t486 = qJD(5) * t469;
	t484 = t463 * t492;
	t483 = t459 * t494;
	t482 = t460 * t494;
	t481 = t464 * t490;
	t474 = -t463 * t468 + t465 * t491;
	t449 = t474 * qJD(2);
	t479 = t461 * t495 - t449;
	t455 = t463 * t491 + t465 * t468;
	t451 = t455 * qJD(2);
	t478 = t461 * t496 - t451;
	t477 = -t474 * t488 + t449;
	t476 = t455 * t488 - t451;
	t475 = (qJD(2) - t488) * t470;
	t454 = t463 * t470 + t465 * t492;
	t473 = t461 * t466 + t464 * t489;
	t450 = t454 * qJD(2);
	t472 = qJD(5) * t454 - t450 * t460 - t474 * t499;
	t452 = -qJD(2) * t484 + t465 * t489;
	t456 = t465 * t470 - t484;
	t471 = qJD(5) * t456 - t452 * t460 + t455 * t499;
	t448 = t466 * t459 + t482;
	t447 = t466 * t460 - t483;
	t446 = t456 * t460 + t459 * t496;
	t445 = -t456 * t459 + t460 * t496;
	t444 = t454 * t460 - t459 * t495;
	t443 = -t454 * t459 - t460 * t495;
	t442 = t460 * t473 - t461 * t483;
	t441 = -t459 * t473 - t461 * t482;
	t440 = -t456 * t499 + t478 * t460;
	t439 = -t456 * t498 - t478 * t459;
	t438 = -t454 * t499 - t479 * t460;
	t437 = -t454 * t498 + t479 * t459;
	t436 = t441 * t469 - t447 * t487;
	t435 = -t441 * t467 - t447 * t486;
	t434 = t439 * t469 - t445 * t487;
	t433 = -t439 * t467 - t445 * t486;
	t432 = t437 * t469 - t443 * t487;
	t431 = -t437 * t467 - t443 * t486;
	t1 = [0, t476 * t467 + t471 * t469, t434, t434, -t440 * t467 + t452 * t469 + (-t446 * t469 - t455 * t467) * qJD(5), 0; 0, t477 * t467 + t472 * t469, t432, t432, -t438 * t467 + t450 * t469 + (-t444 * t469 + t467 * t474) * qJD(5), 0; 0, (t467 * t475 - t500 * t469) * t464, t436, t436, t469 * t481 - t442 * t467 + (-t448 * t469 + t467 * t493) * qJD(5), 0; 0, -t471 * t467 + t476 * t469, t433, t433, -t440 * t469 - t452 * t467 + (t446 * t467 - t455 * t469) * qJD(5), 0; 0, -t472 * t467 + t477 * t469, t431, t431, -t438 * t469 - t450 * t467 + (t444 * t467 + t469 * t474) * qJD(5), 0; 0, (t500 * t467 + t469 * t475) * t464, t435, t435, -t467 * t481 - t442 * t469 + (t448 * t467 + t469 * t493) * qJD(5), 0; 0, -t452 * t459 - t455 * t498, t440, t440, 0, 0; 0, -t450 * t459 + t474 * t498, t438, t438, 0, 0; 0, (-t459 * t490 + t460 * t497) * t464, t442, t442, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end