% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(5));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(5));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0; -t153, -t154, 0, 0, 0; 0, -t158 * t166, 0, 0, 0; t154, t153, 0, 0, 0; t152, t155, 0, 0, 0; 0, -t160 * t166, 0, 0, 0; -t159 * t167, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:05
	% EndTime: 2019-12-31 20:28:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (59->20), mult. (192->41), div. (0->0), fcn. (208->8), ass. (0->25)
	t224 = sin(pkin(10));
	t226 = cos(pkin(10));
	t227 = cos(pkin(5));
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
	t225 = sin(pkin(5));
	t239 = qJD(1) * t225;
	t233 = qJD(2) * t234;
	t216 = t218 * qJD(2);
	t232 = t231 * t212 - t229 * t217 + (t215 * t229 + t218 * t231) * qJD(1);
	t214 = t218 * t227;
	t213 = t227 * t233;
	t211 = -t231 * t213 + t229 * t216 + (t214 * t229 - t231 * t234) * qJD(1);
	t210 = t229 * t213 + t231 * t216 + (t214 * t231 + t229 * t234) * qJD(1);
	t1 = [t232, t210, 0, 0, 0; -t240, t211, 0, 0, 0; 0, -t225 * t233, 0, 0, 0; -t211, t240, 0, 0, 0; t210, t232, 0, 0, 0; 0, t225 * t216, 0, 0, 0; -t229 * t239, 0, 0, 0, 0; t231 * t239, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:06
	% EndTime: 2019-12-31 20:28:06
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (191->44), mult. (594->89), div. (0->0), fcn. (666->10), ass. (0->46)
	t382 = sin(pkin(10));
	t384 = cos(pkin(10));
	t385 = cos(pkin(5));
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
	t383 = sin(pkin(5));
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
	t1 = [t392, -t356 * t389 - t362 * t401, 0, t354, 0; t355, t358 * t389 - t360 * t401, 0, t408, 0; 0, -t365 * t389 + t370 * t401, 0, t366 * t386 + (-t371 * t389 - t385 * t386) * qJD(4), 0; -t408, t356 * t386 - t362 * t400, 0, -t355, 0; t354, -t358 * t386 - t360 * t400, 0, t392, 0; 0, t365 * t386 + t370 * t400, 0, t366 * t389 + (t371 * t386 - t385 * t389) * qJD(4), 0; t358, t357, 0, 0, 0; t356, t363 * qJD(1) - t391 * t367 + t405, 0, 0, 0; 0, -t366, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:28:08
	% EndTime: 2019-12-31 20:28:08
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (555->82), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->69)
	t571 = sin(pkin(10));
	t573 = cos(pkin(10));
	t574 = cos(pkin(5));
	t581 = cos(qJ(2));
	t604 = qJD(2) * t581;
	t577 = sin(qJ(2));
	t605 = qJD(2) * t577;
	t553 = (t571 * t605 - t573 * t604) * t574;
	t590 = t581 * t571 + t577 * t573;
	t559 = t590 * t574;
	t561 = t577 * t571 - t581 * t573;
	t582 = cos(qJ(1));
	t578 = sin(qJ(1));
	t606 = qJD(1) * t578;
	t560 = -t571 * t604 - t573 * t605;
	t608 = t578 * t560;
	t532 = -t559 * t606 + t608 + (-qJD(1) * t561 - t553) * t582;
	t576 = sin(qJ(4));
	t580 = cos(qJ(4));
	t592 = t582 * t559 - t578 * t561;
	t572 = sin(pkin(5));
	t610 = t572 * t582;
	t589 = t576 * t592 + t580 * t610;
	t597 = t572 * t606;
	t526 = t589 * qJD(4) - t532 * t580 - t576 * t597;
	t554 = qJD(2) * t559;
	t558 = t561 * t574;
	t588 = t561 * qJD(2);
	t531 = t558 * t606 + (-qJD(1) * t590 - t554) * t582 + t578 * t588;
	t598 = t576 * t610;
	t537 = -t580 * t592 + t598;
	t541 = -t582 * t558 - t578 * t590;
	t575 = sin(qJ(5));
	t579 = cos(qJ(5));
	t619 = t526 * t579 + t531 * t575 + (-t537 * t575 + t541 * t579) * qJD(5);
	t618 = (t537 * t579 + t541 * t575) * qJD(5) + t526 * t575 - t531 * t579;
	t611 = t572 * t578;
	t603 = qJD(4) * t576;
	t602 = qJD(4) * t580;
	t601 = qJD(5) * t575;
	t600 = qJD(5) * t579;
	t599 = qJD(5) * t580;
	t596 = qJD(1) * t610;
	t544 = t578 * t558 - t582 * t590;
	t583 = -t592 * qJD(1) + t578 * t553 + t582 * t560;
	t595 = -t544 * t599 + t583;
	t591 = -t578 * t559 - t582 * t561;
	t594 = t591 * qJD(1) - t541 * t599 - t582 * t553 + t608;
	t552 = t572 * t588;
	t556 = t561 * t572;
	t593 = t556 * t599 - t552;
	t557 = t590 * t572;
	t547 = t557 * t580 + t574 * t576;
	t546 = -t557 * t576 + t574 * t580;
	t538 = -t576 * t591 + t580 * t611;
	t539 = t576 * t611 + t580 * t591;
	t524 = qJD(4) * t598 - t532 * t576 + t580 * t597 - t592 * t602;
	t528 = t541 * qJD(1) - t578 * t554 - t582 * t588;
	t586 = -qJD(5) * t591 + t528 * t580 + t544 * t603;
	t585 = -qJD(5) * t592 - t531 * t580 + t541 * t603;
	t551 = qJD(2) * t557;
	t584 = qJD(5) * t557 - t551 * t580 + t556 * t603;
	t534 = t546 * qJD(4) - t552 * t580;
	t533 = -t547 * qJD(4) + t552 * t576;
	t523 = t538 * qJD(4) + t576 * t596 + t580 * t583;
	t522 = t539 * qJD(4) + t576 * t583 - t580 * t596;
	t521 = t523 * t579 + t528 * t575 + (-t539 * t575 - t544 * t579) * qJD(5);
	t520 = -t523 * t575 + t528 * t579 + (-t539 * t579 + t544 * t575) * qJD(5);
	t1 = [t619, t595 * t575 - t586 * t579, 0, -t522 * t579 - t538 * t601, t520; t521, t594 * t575 - t585 * t579, 0, t524 * t579 + t589 * t601, t618; 0, t593 * t575 + t584 * t579, 0, t533 * t579 - t546 * t601, -t534 * t575 + t551 * t579 + (-t547 * t579 - t556 * t575) * qJD(5); -t618, t586 * t575 + t595 * t579, 0, t522 * t575 - t538 * t600, -t521; t520, t585 * t575 + t594 * t579, 0, -t524 * t575 + t589 * t600, t619; 0, -t584 * t575 + t593 * t579, 0, -t533 * t575 - t546 * t600, -t534 * t579 - t551 * t575 + (t547 * t575 - t556 * t579) * qJD(5); t524, -t528 * t576 + t544 * t602, 0, t523, 0; t522, t531 * t576 + t541 * t602, 0, -t526, 0; 0, -t551 * t576 - t556 * t602, 0, t534, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end