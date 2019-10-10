% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
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
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->20), mult. (192->41), div. (0->0), fcn. (208->8), ass. (0->25)
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t227 = cos(pkin(6));
	t230 = cos(qJ(2));
	t237 = qJD(2) * t230;
	t228 = sin(qJ(2));
	t238 = qJD(2) * t228;
	t212 = (t224 * t238 - t226 * t237) * t227;
	t234 = t230 * t224 + t228 * t226;
	t215 = t234 * t227;
	t217 = -t224 * t237 - t226 * t238;
	t218 = t228 * t224 - t230 * t226;
	t229 = sin(qJ(1));
	t231 = cos(qJ(1));
	t240 = -t229 * t212 - t231 * t217 + (t215 * t231 - t218 * t229) * qJD(1);
	t225 = sin(pkin(6));
	t239 = qJD(1) * t225;
	t233 = qJD(2) * t234;
	t216 = t218 * qJD(2);
	t232 = t231 * t212 - t229 * t217 + (t215 * t229 + t218 * t231) * qJD(1);
	t214 = t218 * t227;
	t213 = t227 * t233;
	t211 = -t231 * t213 + t229 * t216 + (t214 * t229 - t231 * t234) * qJD(1);
	t210 = t229 * t213 + t231 * t216 + (t214 * t231 + t229 * t234) * qJD(1);
	t1 = [t232, t210, 0, 0, 0, 0; -t240, t211, 0, 0, 0, 0; 0, -t225 * t233, 0, 0, 0, 0; -t211, t240, 0, 0, 0, 0; t210, t232, 0, 0, 0, 0; 0, t225 * t216, 0, 0, 0, 0; -t229 * t239, 0, 0, 0, 0, 0; t231 * t239, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:34
	% EndTime: 2019-10-10 10:13:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (191->44), mult. (594->89), div. (0->0), fcn. (666->10), ass. (0->46)
	t382 = sin(pkin(11));
	t384 = cos(pkin(11));
	t385 = cos(pkin(6));
	t390 = cos(qJ(2));
	t402 = qJD(2) * t390;
	t387 = sin(qJ(2));
	t403 = qJD(2) * t387;
	t367 = (t382 * t403 - t384 * t402) * t385;
	t397 = t390 * t382 + t387 * t384;
	t373 = t397 * t385;
	t375 = t387 * t382 - t390 * t384;
	t391 = cos(qJ(1));
	t388 = sin(qJ(1));
	t404 = qJD(1) * t388;
	t374 = -t382 * t402 - t384 * t403;
	t405 = t388 * t374;
	t359 = -t373 * t404 + t405 + (-qJD(1) * t375 - t367) * t391;
	t361 = t391 * t373 - t388 * t375;
	t386 = sin(qJ(4));
	t389 = cos(qJ(4));
	t383 = sin(pkin(6));
	t399 = t383 * t404;
	t406 = t383 * t391;
	t408 = (-t361 * t389 + t386 * t406) * qJD(4) - t359 * t386 + t389 * t399;
	t407 = t383 * t388;
	t401 = qJD(4) * t386;
	t400 = qJD(4) * t389;
	t398 = qJD(1) * t406;
	t372 = t375 * t385;
	t360 = -t391 * t372 - t388 * t397;
	t362 = t388 * t372 - t391 * t397;
	t363 = -t388 * t373 - t391 * t375;
	t370 = t375 * t383;
	t395 = qJD(2) * t397;
	t394 = t375 * qJD(2);
	t392 = -t359 * t389 + t400 * t406 + (qJD(4) * t361 - t399) * t386;
	t357 = -t361 * qJD(1) + t388 * t367 + t391 * t374;
	t371 = t397 * t383;
	t368 = t385 * t395;
	t366 = qJD(2) * t370;
	t365 = t383 * t395;
	t358 = t362 * qJD(1) - t391 * t368 + t388 * t394;
	t356 = t360 * qJD(1) - t388 * t368 - t391 * t394;
	t355 = t386 * t398 + t357 * t389 + (-t363 * t386 + t389 * t407) * qJD(4);
	t354 = t389 * t398 - t357 * t386 + (-t363 * t389 - t386 * t407) * qJD(4);
	t1 = [t392, -t356 * t389 - t362 * t401, 0, t354, 0, 0; t355, t358 * t389 - t360 * t401, 0, t408, 0, 0; 0, -t365 * t389 + t370 * t401, 0, t366 * t386 + (-t371 * t389 - t385 * t386) * qJD(4), 0, 0; -t408, t356 * t386 - t362 * t400, 0, -t355, 0, 0; t354, -t358 * t386 - t360 * t400, 0, t392, 0, 0; 0, t365 * t386 + t370 * t400, 0, t366 * t389 + (t371 * t386 - t385 * t389) * qJD(4), 0, 0; t358, t357, 0, 0, 0, 0; t356, t363 * qJD(1) - t391 * t367 + t405, 0, 0, 0, 0; 0, -t366, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:35
	% EndTime: 2019-10-10 10:13:36
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (191->44), mult. (594->90), div. (0->0), fcn. (666->10), ass. (0->46)
	t429 = sin(pkin(11));
	t431 = cos(pkin(11));
	t432 = cos(pkin(6));
	t437 = cos(qJ(2));
	t449 = qJD(2) * t437;
	t434 = sin(qJ(2));
	t450 = qJD(2) * t434;
	t414 = (t429 * t450 - t431 * t449) * t432;
	t444 = t437 * t429 + t434 * t431;
	t420 = t444 * t432;
	t422 = t434 * t429 - t437 * t431;
	t438 = cos(qJ(1));
	t435 = sin(qJ(1));
	t451 = qJD(1) * t435;
	t421 = -t429 * t449 - t431 * t450;
	t452 = t435 * t421;
	t406 = -t420 * t451 + t452 + (-qJD(1) * t422 - t414) * t438;
	t408 = t438 * t420 - t435 * t422;
	t433 = sin(qJ(4));
	t436 = cos(qJ(4));
	t430 = sin(pkin(6));
	t446 = t430 * t451;
	t453 = t430 * t438;
	t455 = (-t408 * t436 + t433 * t453) * qJD(4) - t406 * t433 + t436 * t446;
	t454 = t430 * t435;
	t448 = qJD(4) * t433;
	t447 = qJD(4) * t436;
	t445 = qJD(1) * t453;
	t419 = t422 * t432;
	t407 = -t438 * t419 - t435 * t444;
	t409 = t435 * t419 - t438 * t444;
	t410 = -t435 * t420 - t438 * t422;
	t417 = t422 * t430;
	t442 = qJD(2) * t444;
	t441 = t422 * qJD(2);
	t439 = t406 * t436 + t433 * t446 + (-t408 * t433 - t436 * t453) * qJD(4);
	t404 = -t408 * qJD(1) + t435 * t414 + t438 * t421;
	t418 = t444 * t430;
	t415 = t432 * t442;
	t413 = qJD(2) * t417;
	t412 = t430 * t442;
	t405 = t409 * qJD(1) - t438 * t415 + t435 * t441;
	t403 = t407 * qJD(1) - t435 * t415 - t438 * t441;
	t402 = t433 * t445 + t404 * t436 + (-t410 * t433 + t436 * t454) * qJD(4);
	t401 = -t436 * t445 + t404 * t433 + (t410 * t436 + t433 * t454) * qJD(4);
	t1 = [t405, t404, 0, 0, 0, 0; t403, t410 * qJD(1) - t438 * t414 + t452, 0, 0, 0, 0; 0, -t413, 0, 0, 0, 0; t439, t403 * t436 + t409 * t448, 0, t401, 0, 0; -t402, -t405 * t436 + t407 * t448, 0, -t455, 0, 0; 0, t412 * t436 - t417 * t448, 0, -t413 * t433 + (t418 * t436 + t432 * t433) * qJD(4), 0, 0; t455, -t403 * t433 + t409 * t447, 0, t402, 0, 0; t401, t405 * t433 + t407 * t447, 0, t439, 0, 0; 0, -t412 * t433 - t417 * t447, 0, -t413 * t436 + (-t418 * t433 + t432 * t436) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:37
	% EndTime: 2019-10-10 10:13:37
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (555->84), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->69)
	t566 = sin(pkin(11));
	t568 = cos(pkin(11));
	t569 = cos(pkin(6));
	t576 = cos(qJ(2));
	t599 = qJD(2) * t576;
	t572 = sin(qJ(2));
	t600 = qJD(2) * t572;
	t549 = (t566 * t600 - t568 * t599) * t569;
	t585 = t576 * t566 + t572 * t568;
	t555 = t585 * t569;
	t557 = t572 * t566 - t576 * t568;
	t577 = cos(qJ(1));
	t573 = sin(qJ(1));
	t601 = qJD(1) * t573;
	t556 = -t566 * t599 - t568 * t600;
	t603 = t573 * t556;
	t528 = -t555 * t601 + t603 + (-qJD(1) * t557 - t549) * t577;
	t571 = sin(qJ(4));
	t575 = cos(qJ(4));
	t587 = t577 * t555 - t573 * t557;
	t567 = sin(pkin(6));
	t592 = t567 * t601;
	t605 = t567 * t577;
	t593 = t571 * t605;
	t597 = qJD(4) * t575;
	t520 = -qJD(4) * t593 + t528 * t571 - t575 * t592 + t587 * t597;
	t550 = qJD(2) * t555;
	t554 = t557 * t569;
	t583 = t557 * qJD(2);
	t527 = t554 * t601 + (-qJD(1) * t585 - t550) * t577 + t573 * t583;
	t531 = t571 * t587 + t575 * t605;
	t537 = -t577 * t554 - t573 * t585;
	t570 = sin(qJ(6));
	t574 = cos(qJ(6));
	t615 = (t531 * t570 - t537 * t574) * qJD(6) - t520 * t574 - t527 * t570;
	t612 = -t520 * t570 + t527 * t574 + (-t531 * t574 - t537 * t570) * qJD(6);
	t611 = t531 * qJD(4) - t528 * t575 - t571 * t592;
	t606 = t567 * t573;
	t598 = qJD(4) * t571;
	t596 = qJD(6) * t570;
	t595 = qJD(6) * t571;
	t594 = qJD(6) * t574;
	t591 = qJD(1) * t605;
	t540 = t573 * t554 - t577 * t585;
	t578 = -t587 * qJD(1) + t573 * t549 + t577 * t556;
	t590 = t540 * t595 + t578;
	t586 = -t573 * t555 - t577 * t557;
	t589 = t586 * qJD(1) + t537 * t595 - t577 * t549 + t603;
	t548 = t567 * t583;
	t552 = t557 * t567;
	t588 = t552 * t595 + t548;
	t553 = t585 * t567;
	t543 = t553 * t575 + t569 * t571;
	t542 = t553 * t571 - t569 * t575;
	t584 = -t571 * t586 + t575 * t606;
	t535 = t571 * t606 + t575 * t586;
	t524 = t537 * qJD(1) - t573 * t550 - t577 * t583;
	t581 = -qJD(6) * t586 - t524 * t571 + t540 * t597;
	t580 = -qJD(6) * t587 + t527 * t571 + t537 * t597;
	t547 = qJD(2) * t553;
	t579 = -qJD(6) * t553 - t547 * t571 - t552 * t597;
	t532 = t575 * t587 - t593;
	t530 = -t542 * qJD(4) - t548 * t575;
	t529 = t543 * qJD(4) - t548 * t571;
	t519 = t584 * qJD(4) + t571 * t591 + t575 * t578;
	t518 = t535 * qJD(4) + t571 * t578 - t575 * t591;
	t517 = t518 * t570 + t524 * t574 + (t540 * t570 - t574 * t584) * qJD(6);
	t516 = t518 * t574 - t524 * t570 + (t540 * t574 + t570 * t584) * qJD(6);
	t1 = [t612, t581 * t570 + t590 * t574, 0, t519 * t570 + t535 * t594, 0, t516; t517, t580 * t570 + t589 * t574, 0, t532 * t594 - t570 * t611, 0, -t615; 0, t579 * t570 - t588 * t574, 0, t530 * t570 + t543 * t594, 0, t529 * t574 - t547 * t570 + (-t542 * t570 - t552 * t574) * qJD(6); t615, -t590 * t570 + t581 * t574, 0, t519 * t574 - t535 * t596, 0, -t517; t516, -t589 * t570 + t580 * t574, 0, -t532 * t596 - t574 * t611, 0, t612; 0, t588 * t570 + t579 * t574, 0, t530 * t574 - t543 * t596, 0, -t529 * t570 - t547 * t574 + (-t542 * t574 + t552 * t570) * qJD(6); t611, -t524 * t575 - t540 * t598, 0, -t518, 0, 0; t519, t527 * t575 - t537 * t598, 0, -t520, 0, 0; 0, -t547 * t575 + t552 * t598, 0, -t529, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end