% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP2
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
% Datum: 2019-10-09 23:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:03:57
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
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:03:57
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
	% StartTime: 2019-10-09 23:03:59
	% EndTime: 2019-10-09 23:03:59
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
	% StartTime: 2019-10-09 23:04:01
	% EndTime: 2019-10-09 23:04:01
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (388->64), mult. (672->131), div. (0->0), fcn. (736->10), ass. (0->67)
	t560 = qJ(3) + qJ(4);
	t557 = sin(t560);
	t559 = qJD(3) + qJD(4);
	t597 = t557 * t559;
	t558 = cos(t560);
	t596 = t558 * t559;
	t568 = cos(qJ(2));
	t595 = t559 * t568;
	t561 = sin(pkin(11));
	t562 = sin(pkin(6));
	t594 = t561 * t562;
	t563 = cos(pkin(11));
	t593 = t562 * t563;
	t566 = sin(qJ(2));
	t592 = t562 * t566;
	t564 = cos(pkin(6));
	t591 = t564 * t566;
	t590 = t564 * t568;
	t565 = sin(qJ(5));
	t589 = t565 * t568;
	t567 = cos(qJ(5));
	t588 = t567 * t568;
	t587 = qJD(2) * t566;
	t586 = qJD(2) * t568;
	t585 = qJD(5) * t558;
	t584 = qJD(5) * t565;
	t583 = qJD(5) * t567;
	t582 = t561 * t591;
	t581 = t557 * t592;
	t580 = t558 * t592;
	t579 = t562 * t587;
	t578 = -qJD(2) + t585;
	t573 = -t561 * t566 + t563 * t590;
	t547 = t573 * qJD(2);
	t577 = t559 * t593 - t547;
	t553 = t561 * t590 + t563 * t566;
	t549 = t553 * qJD(2);
	t576 = t559 * t594 - t549;
	t575 = -t573 * t585 + t547;
	t574 = t553 * t585 - t549;
	t552 = t561 * t568 + t563 * t591;
	t572 = t559 * t564 + t562 * t586;
	t548 = t552 * qJD(2);
	t571 = qJD(5) * t552 - t548 * t558 - t573 * t597;
	t550 = -qJD(2) * t582 + t563 * t586;
	t554 = t563 * t568 - t582;
	t570 = qJD(5) * t554 - t550 * t558 + t553 * t597;
	t569 = -t557 * t595 + (-qJD(2) * t558 + qJD(5)) * t566;
	t546 = t564 * t557 + t580;
	t545 = t564 * t558 - t581;
	t544 = t554 * t558 + t557 * t594;
	t543 = -t554 * t557 + t558 * t594;
	t542 = t552 * t558 - t557 * t593;
	t541 = -t552 * t557 - t558 * t593;
	t540 = t572 * t558 - t559 * t581;
	t539 = -t572 * t557 - t559 * t580;
	t538 = -t554 * t597 + t576 * t558;
	t537 = -t554 * t596 - t576 * t557;
	t536 = -t552 * t597 - t577 * t558;
	t535 = -t552 * t596 + t577 * t557;
	t534 = t539 * t567 - t545 * t584;
	t533 = t539 * t565 + t545 * t583;
	t532 = t537 * t567 - t543 * t584;
	t531 = t537 * t565 + t543 * t583;
	t530 = t535 * t567 - t541 * t584;
	t529 = t535 * t565 + t541 * t583;
	t1 = [0, t574 * t565 + t570 * t567, t532, t532, -t538 * t565 + t550 * t567 + (-t544 * t567 - t553 * t565) * qJD(5), 0; 0, t575 * t565 + t571 * t567, t530, t530, -t536 * t565 + t548 * t567 + (-t542 * t567 + t565 * t573) * qJD(5), 0; 0, (t569 * t567 - t578 * t589) * t562, t534, t534, t567 * t579 - t540 * t565 + (-t546 * t567 + t562 * t589) * qJD(5), 0; 0, -t550 * t557 - t553 * t596, t538, t538, 0, 0; 0, -t548 * t557 + t573 * t596, t536, t536, 0, 0; 0, (-t557 * t587 + t558 * t595) * t562, t540, t540, 0, 0; 0, t570 * t565 - t574 * t567, t531, t531, t538 * t567 + t550 * t565 + (-t544 * t565 + t553 * t567) * qJD(5), 0; 0, t571 * t565 - t575 * t567, t529, t529, t536 * t567 + t548 * t565 + (-t542 * t565 - t567 * t573) * qJD(5), 0; 0, (t569 * t565 + t578 * t588) * t562, t533, t533, t565 * t579 + t540 * t567 + (-t546 * t565 - t562 * t588) * qJD(5), 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end