% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
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
	% StartTime: 2019-10-10 10:57:46
	% EndTime: 2019-10-10 10:57:46
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->20), mult. (192->41), div. (0->0), fcn. (208->8), ass. (0->25)
	t224 = sin(pkin(12));
	t226 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:57:47
	% EndTime: 2019-10-10 10:57:47
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (191->44), mult. (594->89), div. (0->0), fcn. (666->10), ass. (0->46)
	t382 = sin(pkin(12));
	t384 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:57:49
	% EndTime: 2019-10-10 10:57:50
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (555->82), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->69)
	t571 = sin(pkin(12));
	t573 = cos(pkin(12));
	t574 = cos(pkin(6));
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
	t572 = sin(pkin(6));
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
	t1 = [t619, t595 * t575 - t586 * t579, 0, -t522 * t579 - t538 * t601, t520, 0; t521, t594 * t575 - t585 * t579, 0, t524 * t579 + t589 * t601, t618, 0; 0, t593 * t575 + t584 * t579, 0, t533 * t579 - t546 * t601, -t534 * t575 + t551 * t579 + (-t547 * t579 - t556 * t575) * qJD(5), 0; -t618, t586 * t575 + t595 * t579, 0, t522 * t575 - t538 * t600, -t521, 0; t520, t585 * t575 + t594 * t579, 0, -t524 * t575 + t589 * t600, t619, 0; 0, -t584 * t575 + t593 * t579, 0, -t533 * t575 - t546 * t600, -t534 * t579 - t551 * t575 + (t547 * t575 - t556 * t579) * qJD(5), 0; t524, -t528 * t576 + t544 * t602, 0, t523, 0, 0; t522, t531 * t576 + t541 * t602, 0, -t526, 0, 0; 0, -t551 * t576 - t556 * t602, 0, t534, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:50
	% EndTime: 2019-10-10 10:57:50
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (853->78), mult. (2076->142), div. (0->0), fcn. (2392->12), ass. (0->78)
	t626 = qJ(5) + qJ(6);
	t623 = sin(t626);
	t624 = cos(t626);
	t630 = cos(pkin(6));
	t627 = sin(pkin(12));
	t629 = cos(pkin(12));
	t632 = sin(qJ(2));
	t635 = cos(qJ(2));
	t646 = t635 * t627 + t632 * t629;
	t611 = t646 * t630;
	t606 = qJD(2) * t611;
	t613 = t632 * t627 - t635 * t629;
	t610 = t613 * t630;
	t633 = sin(qJ(1));
	t636 = cos(qJ(1));
	t643 = t613 * qJD(2);
	t666 = qJD(1) * t633;
	t583 = t610 * t666 + (-qJD(1) * t646 - t606) * t636 + t633 * t643;
	t625 = qJD(5) + qJD(6);
	t634 = cos(qJ(4));
	t648 = t636 * t611 - t633 * t613;
	t631 = sin(qJ(4));
	t628 = sin(pkin(6));
	t670 = t628 * t636;
	t661 = t631 * t670;
	t676 = -(-t634 * t648 + t661) * t625 + t583;
	t664 = qJD(2) * t635;
	t665 = qJD(2) * t632;
	t605 = (t627 * t665 - t629 * t664) * t630;
	t612 = -t627 * t664 - t629 * t665;
	t668 = t633 * t612;
	t584 = -t611 * t666 + t668 + (-qJD(1) * t613 - t605) * t636;
	t645 = t631 * t648 + t634 * t670;
	t660 = t628 * t666;
	t577 = -t645 * qJD(4) + t584 * t634 + t631 * t660;
	t593 = -t636 * t610 - t633 * t646;
	t678 = t593 * t625 - t577;
	t570 = t623 * t678 - t624 * t676;
	t571 = t623 * t676 + t624 * t678;
	t580 = t593 * qJD(1) - t633 * t606 - t636 * t643;
	t647 = -t633 * t611 - t636 * t613;
	t671 = t628 * t633;
	t644 = t631 * t671 + t634 * t647;
	t677 = -t644 * t625 + t580;
	t674 = t623 * t625;
	t673 = t624 * t625;
	t672 = t625 * t634;
	t663 = qJD(4) * t631;
	t662 = qJD(4) * t634;
	t659 = qJD(1) * t670;
	t591 = -t631 * t647 + t634 * t671;
	t637 = -t648 * qJD(1) + t633 * t605 + t636 * t612;
	t575 = t591 * qJD(4) + t631 * t659 + t634 * t637;
	t596 = t633 * t610 - t636 * t646;
	t658 = t596 * t625 - t575;
	t609 = t646 * t628;
	t598 = -t609 * t631 + t630 * t634;
	t604 = t628 * t643;
	t587 = t598 * qJD(4) - t604 * t634;
	t608 = t613 * t628;
	t653 = -t608 * t625 - t587;
	t599 = t609 * t634 + t630 * t631;
	t603 = qJD(2) * t609;
	t652 = t599 * t625 - t603;
	t651 = -t596 * t672 + t637;
	t650 = t647 * qJD(1) - t593 * t672 - t636 * t605 + t668;
	t649 = t608 * t672 - t604;
	t576 = qJD(4) * t661 - t584 * t631 + t634 * t660 - t648 * t662;
	t640 = t580 * t634 + t596 * t663 - t625 * t647;
	t639 = -t583 * t634 + t593 * t663 - t625 * t648;
	t638 = -t603 * t634 + t608 * t663 + t609 * t625;
	t586 = -t599 * qJD(4) + t604 * t631;
	t574 = t644 * qJD(4) + t631 * t637 - t634 * t659;
	t573 = t652 * t623 + t653 * t624;
	t572 = t653 * t623 - t652 * t624;
	t569 = t677 * t623 - t658 * t624;
	t568 = t658 * t623 + t677 * t624;
	t1 = [t571, t651 * t623 - t640 * t624, 0, -t574 * t624 - t591 * t674, t568, t568; t569, t650 * t623 - t639 * t624, 0, t576 * t624 + t645 * t674, t570, t570; 0, t649 * t623 + t638 * t624, 0, t586 * t624 - t598 * t674, t572, t572; -t570, t640 * t623 + t651 * t624, 0, t574 * t623 - t591 * t673, -t569, -t569; t568, t639 * t623 + t650 * t624, 0, -t576 * t623 + t645 * t673, t571, t571; 0, -t638 * t623 + t649 * t624, 0, -t586 * t623 - t598 * t673, t573, t573; t576, -t580 * t631 + t596 * t662, 0, t575, 0, 0; t574, t583 * t631 + t593 * t662, 0, t577, 0, 0; 0, -t603 * t631 - t608 * t662, 0, t587, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end