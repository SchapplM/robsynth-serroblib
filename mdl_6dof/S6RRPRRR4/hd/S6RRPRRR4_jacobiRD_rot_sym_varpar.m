% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR4
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
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
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
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
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
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
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
	% StartTime: 2019-10-10 10:55:50
	% EndTime: 2019-10-10 10:55:50
	% DurationCPUTime: 0.23s
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
	% StartTime: 2019-10-10 10:55:50
	% EndTime: 2019-10-10 10:55:51
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (361->42), mult. (780->79), div. (0->0), fcn. (874->10), ass. (0->51)
	t438 = sin(pkin(12));
	t440 = cos(pkin(12));
	t441 = cos(pkin(6));
	t444 = cos(qJ(2));
	t457 = qJD(2) * t444;
	t442 = sin(qJ(2));
	t458 = qJD(2) * t442;
	t419 = (t438 * t458 - t440 * t457) * t441;
	t450 = t438 * t444 + t440 * t442;
	t425 = t450 * t441;
	t428 = t438 * t442 - t440 * t444;
	t436 = qJD(4) + qJD(5);
	t445 = cos(qJ(1));
	t443 = sin(qJ(1));
	t459 = qJD(1) * t443;
	t426 = -t438 * t457 - t440 * t458;
	t460 = t443 * t426;
	t439 = sin(pkin(6));
	t461 = t439 * t445;
	t465 = t425 * t459 - t460 - (-qJD(1) * t428 - t419) * t445 + t436 * t461;
	t451 = -t425 * t443 - t428 * t445;
	t464 = qJD(1) * t461 - t451 * t436;
	t437 = qJ(4) + qJ(5);
	t434 = sin(t437);
	t435 = cos(t437);
	t414 = t425 * t445 - t428 * t443;
	t449 = -t414 * t436 + t439 * t459;
	t404 = t434 * t465 + t449 * t435;
	t463 = t434 * t436;
	t462 = t435 * t436;
	t422 = t428 * t439;
	t418 = qJD(2) * t422;
	t454 = -t436 * t441 + t418;
	t407 = -qJD(1) * t414 + t443 * t419 + t445 * t426;
	t453 = t436 * t439 * t443 + t407;
	t424 = t428 * t441;
	t413 = -t424 * t445 - t443 * t450;
	t415 = t424 * t443 - t445 * t450;
	t447 = qJD(2) * t450;
	t446 = t428 * qJD(2);
	t405 = -t449 * t434 + t435 * t465;
	t423 = t450 * t439;
	t420 = t441 * t447;
	t417 = t439 * t447;
	t411 = t423 * t463 + t435 * t454;
	t410 = -t423 * t462 + t434 * t454;
	t408 = qJD(1) * t415 - t445 * t420 + t443 * t446;
	t406 = qJD(1) * t413 - t443 * t420 - t445 * t446;
	t403 = t434 * t464 + t453 * t435;
	t402 = -t453 * t434 + t435 * t464;
	t1 = [t405, -t406 * t435 - t415 * t463, 0, t402, t402, 0; t403, t408 * t435 - t413 * t463, 0, t404, t404, 0; 0, -t417 * t435 + t422 * t463, 0, t410, t410, 0; -t404, t406 * t434 - t415 * t462, 0, -t403, -t403, 0; t402, -t408 * t434 - t413 * t462, 0, t405, t405, 0; 0, t417 * t434 + t422 * t462, 0, t411, t411, 0; t408, t407, 0, 0, 0, 0; t406, qJD(1) * t451 - t445 * t419 + t460, 0, 0, 0, 0; 0, -t418, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:54
	% EndTime: 2019-10-10 10:55:55
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (920->86), mult. (2011->157), div. (0->0), fcn. (2316->12), ass. (0->79)
	t651 = sin(pkin(12));
	t653 = cos(pkin(12));
	t654 = cos(pkin(6));
	t659 = cos(qJ(2));
	t682 = qJD(2) * t659;
	t656 = sin(qJ(2));
	t683 = qJD(2) * t656;
	t629 = (t651 * t683 - t653 * t682) * t654;
	t667 = t659 * t651 + t656 * t653;
	t635 = t667 * t654;
	t638 = t656 * t651 - t659 * t653;
	t660 = cos(qJ(1));
	t657 = sin(qJ(1));
	t684 = qJD(1) * t657;
	t636 = -t651 * t682 - t653 * t683;
	t686 = t657 * t636;
	t608 = -t635 * t684 + t686 + (-qJD(1) * t638 - t629) * t660;
	t650 = qJ(4) + qJ(5);
	t647 = sin(t650);
	t648 = cos(t650);
	t649 = qJD(4) + qJD(5);
	t669 = t660 * t635 - t657 * t638;
	t652 = sin(pkin(6));
	t677 = t652 * t684;
	t688 = t652 * t660;
	t599 = (-t649 * t669 + t677) * t647 - (t649 * t688 - t608) * t648;
	t630 = qJD(2) * t635;
	t634 = t638 * t654;
	t665 = t638 * qJD(2);
	t607 = t634 * t684 + (-qJD(1) * t667 - t630) * t660 + t657 * t665;
	t678 = t647 * t688;
	t613 = -t648 * t669 + t678;
	t617 = -t660 * t634 - t657 * t667;
	t655 = sin(qJ(6));
	t658 = cos(qJ(6));
	t699 = -t599 * t658 + t607 * t655 + (-t613 * t655 + t617 * t658) * qJD(6);
	t698 = (t613 * t658 + t617 * t655) * qJD(6) - t599 * t655 - t607 * t658;
	t691 = t647 * t649;
	t690 = t648 * t649;
	t689 = t652 * t657;
	t681 = qJD(6) * t648;
	t680 = qJD(6) * t655;
	t679 = qJD(6) * t658;
	t676 = qJD(1) * t688;
	t628 = t652 * t665;
	t675 = t649 * t654 - t628;
	t661 = -t669 * qJD(1) + t657 * t629 + t660 * t636;
	t674 = t649 * t689 + t661;
	t620 = t657 * t634 - t660 * t667;
	t672 = -t620 * t681 + t661;
	t668 = -t657 * t635 - t660 * t638;
	t671 = t668 * qJD(1) - t617 * t681 - t660 * t629 + t686;
	t632 = t638 * t652;
	t670 = t632 * t681 - t628;
	t633 = t667 * t652;
	t598 = -t608 * t647 + t648 * t677 + t649 * t678 - t669 * t690;
	t604 = t617 * qJD(1) - t657 * t630 - t660 * t665;
	t664 = -qJD(6) * t668 + t604 * t648 + t620 * t691;
	t663 = -qJD(6) * t669 - t607 * t648 + t617 * t691;
	t627 = qJD(2) * t633;
	t662 = qJD(6) * t633 - t627 * t648 + t632 * t691;
	t623 = t633 * t648 + t654 * t647;
	t622 = -t633 * t647 + t654 * t648;
	t615 = t647 * t689 + t648 * t668;
	t614 = -t647 * t668 + t648 * t689;
	t611 = -t647 * t669 - t648 * t688;
	t610 = -t633 * t691 + t675 * t648;
	t609 = -t633 * t690 - t675 * t647;
	t602 = t609 * t658 - t622 * t680;
	t601 = -t609 * t655 - t622 * t679;
	t597 = t674 * t648 + (-t649 * t668 + t676) * t647;
	t596 = t674 * t647 - t648 * t676 + t668 * t690;
	t595 = t598 * t658 - t611 * t680;
	t594 = -t598 * t655 - t611 * t679;
	t593 = -t596 * t658 - t614 * t680;
	t592 = t596 * t655 - t614 * t679;
	t591 = t597 * t658 + t604 * t655 + (-t615 * t655 - t620 * t658) * qJD(6);
	t590 = -t597 * t655 + t604 * t658 + (-t615 * t658 + t620 * t655) * qJD(6);
	t1 = [t699, t672 * t655 - t664 * t658, 0, t593, t593, t590; t591, t671 * t655 - t663 * t658, 0, t595, t595, t698; 0, t670 * t655 + t662 * t658, 0, t602, t602, -t610 * t655 + t627 * t658 + (-t623 * t658 - t632 * t655) * qJD(6); -t698, t664 * t655 + t672 * t658, 0, t592, t592, -t591; t590, t663 * t655 + t671 * t658, 0, t594, t594, t699; 0, -t662 * t655 + t670 * t658, 0, t601, t601, -t610 * t658 - t627 * t655 + (t623 * t655 - t632 * t658) * qJD(6); t598, -t604 * t647 + t620 * t690, 0, t597, t597, 0; t596, t607 * t647 + t617 * t690, 0, t599, t599, 0; 0, -t627 * t647 - t632 * t690, 0, t610, t610, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end