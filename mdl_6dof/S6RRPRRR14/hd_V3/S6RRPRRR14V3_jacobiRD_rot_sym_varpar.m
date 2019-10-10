% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR14V3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiRD_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
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
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t31 = sin(qJ(1));
	t38 = qJD(1) * t31;
	t33 = cos(qJ(1));
	t37 = qJD(1) * t33;
	t30 = sin(qJ(2));
	t36 = qJD(2) * t30;
	t32 = cos(qJ(2));
	t35 = qJD(2) * t32;
	t34 = qJD(2) * t33;
	t29 = t31 * t36 - t32 * t37;
	t28 = t30 * t37 + t31 * t35;
	t27 = t30 * t34 + t32 * t38;
	t26 = t30 * t38 - t32 * t34;
	t1 = [t29, t26, 0, 0, 0, 0; -t27, -t28, 0, 0, 0, 0; 0, -t36, 0, 0, 0, 0; t28, t27, 0, 0, 0, 0; t26, t29, 0, 0, 0, 0; 0, -t35, 0, 0, 0, 0; -t38, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:24
	% EndTime: 2019-10-10 11:13:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t154 = sin(qJ(1));
	t161 = qJD(1) * t154;
	t156 = cos(qJ(1));
	t160 = qJD(1) * t156;
	t153 = sin(qJ(2));
	t159 = qJD(2) * t153;
	t155 = cos(qJ(2));
	t158 = qJD(2) * t155;
	t157 = qJD(2) * t156;
	t152 = -t154 * t159 + t155 * t160;
	t151 = -t153 * t160 - t154 * t158;
	t150 = -t153 * t157 - t155 * t161;
	t149 = t153 * t161 - t155 * t157;
	t1 = [-t152, t149, 0, 0, 0, 0; t150, t151, 0, 0, 0, 0; 0, -t159, 0, 0, 0, 0; -t161, 0, 0, 0, 0, 0; t160, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t151, t150, 0, 0, 0, 0; -t149, t152, 0, 0, 0, 0; 0, t158, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:24
	% EndTime: 2019-10-10 11:13:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(4));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(4));
	t249 = qJD(4) * t229;
	t248 = qJD(4) * t230;
	t247 = qJD(4) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(4) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, 0, t225, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, 0, -t227, 0, 0; 0, -t229 * t247 - t246, 0, -t229 * t251 - t245, 0, 0; t227, t236 * t229 + t234 * t245, 0, t226, 0, 0; t225, t237 * t229 + t231 * t245, 0, t228, 0, 0; 0, t229 * t252 - t232 * t247, 0, t229 * t248 - t232 * t251, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0, 0; 0, t251, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:26
	% EndTime: 2019-10-10 11:13:26
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (162->55), mult. (548->107), div. (0->0), fcn. (560->8), ass. (0->60)
	t376 = sin(qJ(5));
	t378 = sin(qJ(2));
	t380 = cos(qJ(5));
	t382 = cos(qJ(2));
	t381 = cos(qJ(4));
	t398 = qJD(2) * t381 - qJD(5);
	t393 = t398 * t382;
	t399 = qJD(5) * t381 - qJD(2);
	t377 = sin(qJ(4));
	t413 = qJD(4) * t377;
	t403 = t378 * t413;
	t386 = t399 * t380 * t378 + (t393 - t403) * t376;
	t411 = qJD(4) * t382;
	t405 = t377 * t411;
	t427 = t398 * t378 + t405;
	t379 = sin(qJ(1));
	t400 = qJD(1) * t382 - qJD(4);
	t383 = cos(qJ(1));
	t414 = qJD(2) * t383;
	t426 = t378 * t414 + t400 * t379;
	t415 = qJD(2) * t382;
	t417 = qJD(1) * t383;
	t390 = -t378 * t417 - t379 * t415;
	t419 = t383 * t377;
	t421 = t379 * t382;
	t388 = -qJD(5) * (t381 * t421 - t419) - t390;
	t420 = t381 * t383;
	t371 = t379 * t377 + t382 * t420;
	t412 = qJD(4) * t381;
	t402 = t383 * t412;
	t416 = qJD(2) * t378;
	t407 = t379 * t416;
	t367 = t371 * qJD(1) - t379 * t405 - t381 * t407 - t402;
	t409 = qJD(5) * t378;
	t396 = t379 * t409 + t367;
	t425 = t396 * t376 - t388 * t380;
	t423 = t378 * t381;
	t422 = t379 * t381;
	t418 = qJD(1) * t379;
	t410 = qJD(5) * t376;
	t408 = qJD(5) * t380;
	t404 = t378 * t412;
	t401 = qJD(1) - t411;
	t395 = t401 * t383;
	t365 = t377 * t395 - t426 * t381;
	t397 = t383 * t409 + t365;
	t394 = t399 * t382;
	t392 = -t376 * t382 + t380 * t423;
	t391 = t376 * t423 + t380 * t382;
	t389 = t378 * t418 - t382 * t414;
	t387 = -qJD(5) * t371 - t389;
	t385 = -t380 * t393 + (t399 * t376 + t380 * t413) * t378;
	t384 = -t388 * t376 - t396 * t380;
	t370 = -t382 * t419 + t422;
	t368 = -t377 * t421 - t420;
	t366 = t401 * t422 + (-t400 * t383 + t407) * t377;
	t364 = t426 * t377 + t381 * t395;
	t363 = t387 * t376 + t397 * t380;
	t362 = -t397 * t376 + t387 * t380;
	t1 = [t384, t385 * t383 + t392 * t418, 0, t364 * t380 - t370 * t410, t362, 0; t363, t385 * t379 - t392 * t417, 0, t366 * t380 - t368 * t410, -t425, 0; 0, -t376 * t394 - t427 * t380, 0, -t380 * t404 + (t376 * t409 - t380 * t415) * t377, -t386, 0; t425, t386 * t383 - t391 * t418, 0, -t364 * t376 - t370 * t408, -t363, 0; t362, t386 * t379 + t391 * t417, 0, -t366 * t376 - t368 * t408, t384, 0; 0, t427 * t376 - t380 * t394, 0, t378 * t377 * t408 + (t377 * t415 + t404) * t376, t385, 0; t366, t389 * t377 - t378 * t402, 0, t365, 0, 0; -t364, t390 * t377 - t379 * t404, 0, t367, 0, 0; 0, -t377 * t416 + t381 * t411, 0, t381 * t415 - t403, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:28
	% EndTime: 2019-10-10 11:13:30
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (438->113), mult. (1411->208), div. (0->0), fcn. (1483->10), ass. (0->94)
	t592 = sin(qJ(4));
	t593 = sin(qJ(2));
	t591 = sin(qJ(5));
	t597 = cos(qJ(4));
	t625 = qJD(5) * t597 - qJD(2);
	t616 = t625 * t591;
	t596 = cos(qJ(5));
	t647 = qJD(4) * t596;
	t676 = t593 * (t592 * t647 + t616);
	t594 = sin(qJ(1));
	t598 = cos(qJ(2));
	t599 = cos(qJ(1));
	t654 = t599 * t597;
	t576 = t594 * t592 + t598 * t654;
	t644 = qJD(4) * t599;
	t626 = t597 * t644;
	t648 = qJD(4) * t592;
	t628 = t594 * t648;
	t651 = qJD(2) * t593;
	t634 = t594 * t651;
	t562 = t576 * qJD(1) - t597 * t634 - t598 * t628 - t626;
	t655 = t599 * t592;
	t658 = t594 * t598;
	t573 = t597 * t658 - t655;
	t642 = qJD(5) * t593;
	t650 = qJD(2) * t598;
	t652 = qJD(1) * t599;
	t671 = t593 * t652 + t594 * t650;
	t552 = (-qJD(5) * t573 + t671) * t591 + (t594 * t642 + t562) * t596;
	t646 = qJD(4) * t597;
	t627 = t594 * t646;
	t653 = qJD(1) * t594;
	t561 = -t597 * t653 + t627 * t598 + (t652 * t598 - t634 - t644) * t592;
	t662 = t593 * t591;
	t564 = t573 * t596 + t594 * t662;
	t572 = t592 * t658 + t654;
	t590 = sin(qJ(6));
	t595 = cos(qJ(6));
	t675 = -t552 * t595 - t561 * t590 + (t564 * t590 - t572 * t595) * qJD(6);
	t674 = qJD(6) * (t564 * t595 + t572 * t590) + t552 * t590 - t561 * t595;
	t624 = qJD(2) * t597 - qJD(5);
	t615 = t624 * t598;
	t639 = qJD(6) * t592;
	t629 = t593 * t639;
	t673 = t596 * t615 + t629 - t676;
	t649 = qJD(2) * t599;
	t610 = t593 * t653 - t598 * t649;
	t660 = t593 * t597;
	t571 = -t598 * t591 + t596 * t660;
	t668 = qJD(6) * t571;
	t663 = t592 * t596;
	t661 = t593 * t596;
	t659 = t593 * t599;
	t657 = t595 * t597;
	t656 = t598 * t596;
	t645 = qJD(4) * t598;
	t643 = qJD(5) * t591;
	t641 = qJD(5) * t596;
	t640 = qJD(6) * t590;
	t638 = qJD(6) * t595;
	t637 = qJD(6) * t596;
	t633 = t593 * t649;
	t623 = -qJD(4) + t637;
	t622 = -qJD(6) + t647;
	t606 = -t592 * t645 - t624 * t593;
	t620 = t606 * t596 + (-t616 + t639) * t598;
	t558 = t624 * t656 - t676;
	t619 = -t558 - t629;
	t560 = (qJD(1) - t645) * t655 + (-t633 + (-qJD(1) * t598 + qJD(4)) * t594) * t597;
	t575 = -t594 * t597 + t598 * t655;
	t618 = t575 * t637 + t560;
	t617 = t572 * t637 + t562;
	t614 = -t571 * t653 + t599 * t673;
	t569 = t571 * t599;
	t613 = qJD(1) * t569 + t594 * t673;
	t567 = t576 * t596 + t591 * t659;
	t612 = t591 * t660 + t656;
	t609 = -t592 * t650 - t593 * t646;
	t559 = t572 * qJD(1) + t592 * t633 - t598 * t626 - t628;
	t608 = qJD(6) * t576 + t559 * t596 + t575 * t643;
	t607 = qJD(6) * t573 - t561 * t596 + t572 * t643;
	t604 = -t609 - t668;
	t603 = -qJD(6) * (t597 * t656 + t662) - t592 * t651 + t597 * t645;
	t551 = -t564 * qJD(5) - t562 * t591 + t596 * t671;
	t602 = -t592 * t671 - t593 * t627 + t594 * t668;
	t601 = qJD(6) * t569 + t610 * t592 - t593 * t626;
	t557 = -t625 * t661 + (t593 * t648 - t615) * t591;
	t566 = -t576 * t591 + t596 * t659;
	t563 = -t573 * t591 + t594 * t661;
	t550 = (t599 * t642 + t560) * t596 + (-qJD(5) * t576 - t610) * t591;
	t549 = t567 * qJD(5) + t560 * t591 + t610 * t596;
	t548 = t550 * t595 - t559 * t590 + (-t567 * t590 + t575 * t595) * qJD(6);
	t547 = -t550 * t590 - t559 * t595 + (-t567 * t595 - t575 * t590) * qJD(6);
	t1 = [t675, t601 * t590 - t614 * t595, 0, t618 * t590 + t608 * t595, -t549 * t595 - t566 * t640, t547; t548, t602 * t590 - t613 * t595, 0, t617 * t590 + t607 * t595, t551 * t595 - t563 * t640, -t674; 0, t603 * t590 + t620 * t595, 0, (t590 * t597 - t595 * t663) * t650 + (-t622 * t657 + (t623 * t590 + t595 * t643) * t592) * t593, t557 * t595 + t612 * t640, t619 * t590 + t604 * t595; t674, t614 * t590 + t601 * t595, 0, -t608 * t590 + t618 * t595, t549 * t590 - t566 * t638, -t548; t547, t613 * t590 + t602 * t595, 0, -t607 * t590 + t617 * t595, -t551 * t590 - t563 * t638, t675; 0, -t620 * t590 + t603 * t595, 0, (t590 * t663 + t657) * t650 + (t623 * t595 * t592 + (-t592 * t643 + t622 * t597) * t590) * t593, -t557 * t590 + t612 * t638, -t604 * t590 + t619 * t595; t551, t557 * t599 + t612 * t653, 0, t559 * t591 - t575 * t641, t550, 0; t549, t557 * t594 - t612 * t652, 0, -t561 * t591 - t572 * t641, t552, 0; 0, t606 * t591 + t625 * t656, 0, -t593 * t592 * t641 + t609 * t591, t558, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end