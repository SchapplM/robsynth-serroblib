% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
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
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
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
	% StartTime: 2019-10-10 11:07:03
	% EndTime: 2019-10-10 11:07:03
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
	% StartTime: 2019-10-10 11:07:03
	% EndTime: 2019-10-10 11:07:03
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
	% StartTime: 2019-10-10 11:07:03
	% EndTime: 2019-10-10 11:07:04
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (237->36), mult. (414->59), div. (0->0), fcn. (430->8), ass. (0->43)
	t350 = cos(pkin(6));
	t354 = cos(qJ(1));
	t353 = cos(qJ(2));
	t366 = t354 * t353;
	t361 = t350 * t366;
	t364 = qJD(2) * t353;
	t351 = sin(qJ(2));
	t352 = sin(qJ(1));
	t369 = t352 * t351;
	t333 = -qJD(1) * t361 - t354 * t364 + (qJD(2) * t350 + qJD(1)) * t369;
	t347 = qJD(4) + qJD(5);
	t349 = sin(pkin(6));
	t371 = t347 * t349;
	t375 = -t352 * t371 - t333;
	t348 = qJ(4) + qJ(5);
	t345 = sin(t348);
	t346 = cos(t348);
	t358 = t361 - t369;
	t365 = qJD(1) * t349;
	t357 = -t347 * t358 + t352 * t365;
	t367 = t354 * t351;
	t368 = t352 * t353;
	t338 = t350 * t367 + t368;
	t339 = t350 * t368 + t367;
	t335 = t339 * qJD(1) + t338 * qJD(2);
	t359 = t354 * t371 + t335;
	t374 = t357 * t345 - t359 * t346;
	t373 = t345 * t347;
	t372 = t346 * t347;
	t370 = t347 * t351;
	t362 = t353 * t371;
	t360 = t349 * qJD(2) * t351;
	t340 = -t350 * t369 + t366;
	t356 = t339 * t347 + t354 * t365;
	t355 = -t347 * t350 + t360;
	t328 = -t359 * t345 - t357 * t346;
	t336 = t340 * qJD(1) + t358 * qJD(2);
	t334 = -t338 * qJD(1) - t339 * qJD(2);
	t332 = t345 * t362 + t355 * t346;
	t331 = -t355 * t345 + t346 * t362;
	t330 = t375 * t345 + t356 * t346;
	t329 = -t356 * t345 + t375 * t346;
	t1 = [t328, t334 * t345 + t340 * t372, 0, t329, t329, 0; t330, t336 * t345 + t338 * t372, 0, -t374, -t374, 0; 0, (t345 * t364 + t346 * t370) * t349, 0, t332, t332, 0; t374, t334 * t346 - t340 * t373, 0, -t330, -t330, 0; t329, t336 * t346 - t338 * t373, 0, t328, t328, 0; 0, (-t345 * t370 + t346 * t364) * t349, 0, t331, t331, 0; -t336, t333, 0, 0, 0, 0; t334, -t335, 0, 0, 0, 0; 0, -t360, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:07
	% EndTime: 2019-10-10 11:07:07
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (602->76), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->75)
	t549 = cos(pkin(6));
	t551 = sin(qJ(2));
	t555 = cos(qJ(1));
	t580 = t555 * t551;
	t552 = sin(qJ(1));
	t554 = cos(qJ(2));
	t581 = t552 * t554;
	t537 = t549 * t580 + t581;
	t538 = t549 * t581 + t580;
	t527 = t538 * qJD(1) + t537 * qJD(2);
	t547 = qJ(4) + qJ(5);
	t544 = sin(t547);
	t545 = cos(t547);
	t546 = qJD(4) + qJD(5);
	t579 = t555 * t554;
	t582 = t552 * t551;
	t536 = -t549 * t579 + t582;
	t548 = sin(pkin(6));
	t578 = qJD(1) * t548;
	t560 = t536 * t546 + t552 * t578;
	t584 = t548 * t555;
	t516 = (t546 * t584 + t527) * t544 + t560 * t545;
	t569 = t549 * t582;
	t577 = qJD(2) * t551;
	t528 = -qJD(1) * t569 - t552 * t577 + (qJD(2) * t549 + qJD(1)) * t579;
	t568 = t545 * t584;
	t532 = -t536 * t544 + t568;
	t550 = sin(qJ(6));
	t553 = cos(qJ(6));
	t598 = -t516 * t553 + (-t532 * t550 - t537 * t553) * qJD(6) - t528 * t550;
	t597 = (t532 * t553 - t537 * t550) * qJD(6) - t516 * t550 + t528 * t553;
	t594 = (qJD(2) * t544 + qJD(6)) * t554;
	t589 = t544 * t546;
	t588 = t545 * t546;
	t587 = t546 * t551;
	t586 = t548 * t552;
	t585 = t548 * t554;
	t583 = t551 * t553;
	t576 = qJD(2) * t554;
	t575 = qJD(6) * t544;
	t574 = qJD(6) * t550;
	t573 = qJD(6) * t553;
	t572 = t544 * t585;
	t571 = t545 * t585;
	t570 = t545 * t586;
	t567 = t548 * t576;
	t566 = -qJD(2) - t575;
	t561 = t569 - t579;
	t525 = t536 * qJD(1) + t561 * qJD(2);
	t563 = t561 * t575 + t525;
	t562 = -t537 * t575 - t527;
	t559 = t538 * t546 + t555 * t578;
	t558 = -t546 * t549 + t548 * t577;
	t526 = -t537 * qJD(1) - t538 * qJD(2);
	t557 = -qJD(6) * t538 + t526 * t544 - t561 * t588;
	t556 = -qJD(6) * t536 + t528 * t544 + t537 * t588;
	t515 = t527 * t545 - t560 * t544 + t546 * t568;
	t535 = t549 * t545 - t572;
	t534 = -t549 * t544 - t571;
	t531 = t536 * t545 + t544 * t584;
	t530 = t538 * t544 + t570;
	t529 = t538 * t545 - t544 * t586;
	t524 = t558 * t545 + t546 * t572;
	t523 = t558 * t544 - t546 * t571;
	t521 = t524 * t553 - t534 * t574;
	t520 = -t524 * t550 - t534 * t573;
	t519 = t559 * t545 + (-t546 * t586 - t525) * t544;
	t518 = t525 * t545 + t559 * t544 + t546 * t570;
	t514 = -t518 * t553 - t529 * t574;
	t513 = t518 * t550 - t529 * t573;
	t512 = t515 * t553 - t531 * t574;
	t511 = -t515 * t550 - t531 * t573;
	t510 = t519 * t553 + t526 * t550 + (-t530 * t550 - t553 * t561) * qJD(6);
	t509 = -t519 * t550 + t526 * t553 + (-t530 * t553 + t550 * t561) * qJD(6);
	t1 = [t598, t563 * t550 + t557 * t553, 0, t514, t514, t509; t510, t562 * t550 + t556 * t553, 0, t512, t512, t597; 0, (t553 * t594 + (t566 * t550 + t553 * t588) * t551) * t548, 0, t521, t521, t553 * t567 - t523 * t550 + (-t548 * t550 * t551 - t535 * t553) * qJD(6); -t597, -t557 * t550 + t563 * t553, 0, t513, t513, -t510; t509, -t556 * t550 + t562 * t553, 0, t511, t511, t598; 0, (t566 * t583 + (-t545 * t587 - t594) * t550) * t548, 0, t520, t520, -t550 * t567 - t523 * t553 + (t535 * t550 - t548 * t583) * qJD(6); t515, -t526 * t545 - t561 * t589, 0, t519, t519, 0; t518, -t528 * t545 + t537 * t589, 0, t516, t516, 0; 0, (t544 * t587 - t545 * t576) * t548, 0, t523, t523, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end