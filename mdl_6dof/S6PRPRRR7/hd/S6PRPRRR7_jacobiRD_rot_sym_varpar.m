% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiRD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(13));
	t58 = sin(pkin(13));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:03
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->12), mult. (81->36), div. (0->0), fcn. (81->10), ass. (0->20)
	t198 = sin(pkin(14));
	t204 = cos(pkin(7));
	t213 = t198 * t204;
	t202 = cos(pkin(14));
	t212 = t202 * t204;
	t207 = cos(qJ(2));
	t211 = t204 * t207;
	t205 = cos(pkin(6));
	t206 = sin(qJ(2));
	t210 = t205 * t206;
	t209 = t205 * t207;
	t208 = qJD(2) * sin(pkin(6));
	t203 = cos(pkin(13));
	t200 = sin(pkin(7));
	t199 = sin(pkin(13));
	t197 = (t199 * t210 - t203 * t207) * qJD(2);
	t196 = (t199 * t209 + t203 * t206) * qJD(2);
	t195 = (-t199 * t207 - t203 * t210) * qJD(2);
	t194 = (t199 * t206 - t203 * t209) * qJD(2);
	t1 = [0, t196 * t213 + t197 * t202, 0, 0, 0, 0; 0, t194 * t213 + t195 * t202, 0, 0, 0, 0; 0, (-t198 * t211 - t202 * t206) * t208, 0, 0, 0, 0; 0, t196 * t212 - t197 * t198, 0, 0, 0, 0; 0, t194 * t212 - t195 * t198, 0, 0, 0, 0; 0, (t198 * t206 - t202 * t211) * t208, 0, 0, 0, 0; 0, -t196 * t200, 0, 0, 0, 0; 0, -t194 * t200, 0, 0, 0, 0; 0, t207 * t200 * t208, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:05
	% EndTime: 2019-10-09 22:05:06
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (218->82), mult. (800->175), div. (0->0), fcn. (882->14), ass. (0->68)
	t451 = sin(pkin(14));
	t453 = sin(pkin(8));
	t455 = sin(pkin(6));
	t456 = cos(pkin(14));
	t458 = cos(pkin(8));
	t462 = sin(qJ(2));
	t464 = cos(qJ(2));
	t459 = cos(pkin(7));
	t490 = t459 * t462;
	t454 = sin(pkin(7));
	t496 = t454 * t455;
	t467 = t453 * t462 * t496 + (-t451 * t464 - t456 * t490) * t455 * t458;
	t452 = sin(pkin(13));
	t457 = cos(pkin(13));
	t460 = cos(pkin(6));
	t487 = t460 * t464;
	t445 = -t452 * t462 + t457 * t487;
	t499 = t451 * t459;
	t497 = t453 * t454;
	t495 = t454 * t458;
	t494 = t454 * t460;
	t493 = t455 * t459;
	t492 = t456 * t459;
	t491 = t456 * t462;
	t489 = t459 * t464;
	t488 = t460 * t462;
	t486 = qJD(2) * t455;
	t485 = t464 * t496;
	t482 = qJD(2) * t485;
	t446 = t452 * t464 + t457 * t488;
	t472 = t445 * t459 - t457 * t496;
	t481 = -(-t446 * t451 + t472 * t456) * t458 - (-t445 * t454 - t457 * t493) * t453;
	t468 = t452 * t488 - t457 * t464;
	t469 = t452 * t487 + t457 * t462;
	t471 = t452 * t496 - t459 * t469;
	t480 = -(t451 * t468 + t471 * t456) * t458 - (t452 * t493 + t454 * t469) * t453;
	t470 = -t451 * t462 + t456 * t489;
	t479 = -(t470 * t455 + t456 * t494) * t458 - (t460 * t459 - t485) * t453;
	t440 = t445 * qJD(2);
	t441 = t446 * qJD(2);
	t478 = -(-t440 * t451 - t441 * t492) * t458 - t441 * t497;
	t420 = -t440 * t492 + t441 * t451;
	t477 = t420 * t458 + t440 * t497;
	t442 = t469 * qJD(2);
	t443 = t468 * qJD(2);
	t476 = -(t442 * t451 + t443 * t492) * t458 + t443 * t497;
	t424 = t442 * t492 - t443 * t451;
	t475 = t424 * t458 - t442 * t497;
	t474 = (-t445 * t451 - t446 * t492) * t458 + t446 * t497;
	t473 = (t451 * t469 + t468 * t492) * t458 - t468 * t497;
	t439 = (-t451 * t490 + t456 * t464) * t455;
	t434 = t470 * t486;
	t466 = -t434 * t458 + t453 * t482;
	t465 = t467 * qJD(2);
	t463 = cos(qJ(4));
	t461 = sin(qJ(4));
	t437 = qJD(2) * t439;
	t435 = (-t451 * t489 - t491) * t486;
	t431 = t455 * t491 + (t455 * t489 + t494) * t451;
	t429 = -t456 * t469 + t468 * t499;
	t427 = t445 * t456 - t446 * t499;
	t425 = t442 * t499 + t443 * t456;
	t423 = -t442 * t456 + t443 * t499;
	t421 = -t440 * t499 - t441 * t456;
	t419 = t440 * t456 - t441 * t499;
	t417 = t471 * t451 - t456 * t468;
	t415 = t446 * t456 + t472 * t451;
	t1 = [0, t425 * t463 + t475 * t461 + (-t429 * t461 + t473 * t463) * qJD(4), 0, -t423 * t461 - t476 * t463 + (-t417 * t463 + t480 * t461) * qJD(4), 0, 0; 0, t421 * t463 + t477 * t461 + (-t427 * t461 + t474 * t463) * qJD(4), 0, -t419 * t461 - t478 * t463 + (-t415 * t463 + t481 * t461) * qJD(4), 0, 0; 0, t435 * t463 + t466 * t461 + (-t439 * t461 + t467 * t463) * qJD(4), 0, -t437 * t461 + t465 * t463 + (-t431 * t463 + t479 * t461) * qJD(4), 0, 0; 0, -t425 * t461 + t475 * t463 + (-t429 * t463 - t473 * t461) * qJD(4), 0, -t423 * t463 + t476 * t461 + (t417 * t461 + t480 * t463) * qJD(4), 0, 0; 0, -t421 * t461 + t477 * t463 + (-t427 * t463 - t474 * t461) * qJD(4), 0, -t419 * t463 + t478 * t461 + (t415 * t461 + t481 * t463) * qJD(4), 0, 0; 0, -t435 * t461 + t466 * t463 + (-t439 * t463 - t467 * t461) * qJD(4), 0, -t437 * t463 - t465 * t461 + (t431 * t461 + t479 * t463) * qJD(4), 0, 0; 0, -t424 * t453 - t442 * t495, 0, 0, 0, 0; 0, -t420 * t453 + t440 * t495, 0, 0, 0, 0; 0, t434 * t453 + t458 * t482, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:09
	% EndTime: 2019-10-09 22:05:10
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (739->134), mult. (2511->271), div. (0->0), fcn. (2912->16), ass. (0->116)
	t717 = sin(pkin(13));
	t722 = cos(pkin(13));
	t728 = sin(qJ(2));
	t725 = cos(pkin(6));
	t731 = cos(qJ(2));
	t757 = t725 * t731;
	t710 = -t717 * t728 + t722 * t757;
	t716 = sin(pkin(14));
	t720 = sin(pkin(6));
	t724 = cos(pkin(7));
	t759 = t724 * t731;
	t721 = cos(pkin(14));
	t761 = t721 * t728;
	t719 = sin(pkin(7));
	t764 = t719 * t725;
	t696 = t720 * t761 + (t720 * t759 + t764) * t716;
	t727 = sin(qJ(4));
	t730 = cos(qJ(4));
	t737 = -t716 * t728 + t721 * t759;
	t695 = t737 * t720 + t721 * t764;
	t766 = t719 * t720;
	t709 = t725 * t724 - t731 * t766;
	t718 = sin(pkin(8));
	t723 = cos(pkin(8));
	t746 = t695 * t723 + t709 * t718;
	t665 = t696 * t730 + t746 * t727;
	t758 = t725 * t728;
	t735 = t717 * t758 - t722 * t731;
	t736 = t717 * t757 + t722 * t728;
	t738 = t717 * t766 - t724 * t736;
	t679 = t738 * t716 - t721 * t735;
	t678 = t716 * t735 + t738 * t721;
	t763 = t720 * t724;
	t698 = t717 * t763 + t719 * t736;
	t747 = t678 * t723 + t698 * t718;
	t660 = t679 * t730 + t747 * t727;
	t711 = t717 * t731 + t722 * t758;
	t739 = t710 * t724 - t722 * t766;
	t677 = t711 * t721 + t739 * t716;
	t676 = -t711 * t716 + t739 * t721;
	t697 = -t710 * t719 - t722 * t763;
	t748 = t676 * t723 + t697 * t718;
	t658 = t677 * t730 + t748 * t727;
	t769 = t716 * t724;
	t767 = t718 * t719;
	t765 = t719 * t723;
	t762 = t721 * t724;
	t760 = t724 * t728;
	t756 = qJD(2) * t720;
	t726 = sin(qJ(5));
	t755 = qJD(5) * t726;
	t729 = cos(qJ(5));
	t754 = qJD(5) * t729;
	t753 = t728 * t766;
	t751 = t719 * t756;
	t750 = t728 * t751;
	t749 = t731 * t751;
	t705 = t710 * qJD(2);
	t706 = t711 * qJD(2);
	t680 = -t705 * t716 - t706 * t762;
	t745 = t680 * t723 + t706 * t767;
	t682 = -t705 * t762 + t706 * t716;
	t744 = -t682 * t723 - t705 * t767;
	t707 = t736 * qJD(2);
	t708 = t735 * qJD(2);
	t684 = t707 * t716 + t708 * t762;
	t743 = t684 * t723 - t708 * t767;
	t686 = t707 * t762 - t708 * t716;
	t742 = -t686 * t723 + t707 * t767;
	t688 = -t710 * t716 - t711 * t762;
	t741 = t688 * t723 + t711 * t767;
	t690 = t716 * t736 + t735 * t762;
	t740 = t690 * t723 - t735 * t767;
	t703 = (-t716 * t731 - t721 * t760) * t720;
	t734 = t703 * t723 + t718 * t753;
	t704 = (-t716 * t760 + t721 * t731) * t720;
	t699 = t737 * t756;
	t733 = -t699 * t723 + t718 * t749;
	t701 = qJD(2) * t703;
	t732 = t701 * t723 + t718 * t750;
	t657 = -t677 * t727 + t748 * t730;
	t659 = -t679 * t727 + t747 * t730;
	t664 = -t696 * t727 + t746 * t730;
	t689 = t710 * t721 - t711 * t769;
	t662 = t689 * t730 + t741 * t727;
	t691 = -t721 * t736 + t735 * t769;
	t663 = t691 * t730 + t740 * t727;
	t674 = t704 * t730 + t734 * t727;
	t702 = qJD(2) * t704;
	t700 = (-t716 * t759 - t761) * t756;
	t694 = -t703 * t718 + t723 * t753;
	t693 = -t701 * t718 + t723 * t750;
	t692 = t699 * t718 + t723 * t749;
	t687 = t707 * t769 + t708 * t721;
	t685 = -t707 * t721 + t708 * t769;
	t683 = -t705 * t769 - t706 * t721;
	t681 = t705 * t721 - t706 * t769;
	t675 = -t695 * t718 + t709 * t723;
	t673 = -t690 * t718 - t735 * t765;
	t672 = -t688 * t718 + t711 * t765;
	t671 = -t686 * t718 - t707 * t765;
	t670 = -t684 * t718 - t708 * t765;
	t669 = -t682 * t718 + t705 * t765;
	t668 = -t680 * t718 + t706 * t765;
	t667 = -t678 * t718 + t698 * t723;
	t666 = -t676 * t718 + t697 * t723;
	t661 = t700 * t730 + t733 * t727 + (-t704 * t727 + t734 * t730) * qJD(4);
	t656 = t664 * qJD(4) + t702 * t730 + t732 * t727;
	t655 = -t665 * qJD(4) - t702 * t727 + t732 * t730;
	t654 = t687 * t730 - t742 * t727 + (-t691 * t727 + t740 * t730) * qJD(4);
	t653 = t683 * t730 - t744 * t727 + (-t689 * t727 + t741 * t730) * qJD(4);
	t652 = t659 * qJD(4) + t685 * t730 + t743 * t727;
	t651 = -t660 * qJD(4) - t685 * t727 + t743 * t730;
	t650 = t657 * qJD(4) + t681 * t730 + t745 * t727;
	t649 = -t658 * qJD(4) - t681 * t727 + t745 * t730;
	t1 = [0, t654 * t729 + t671 * t726 + (-t663 * t726 + t673 * t729) * qJD(5), 0, t651 * t729 - t659 * t755, -t652 * t726 + t670 * t729 + (-t660 * t729 - t667 * t726) * qJD(5), 0; 0, t653 * t729 + t669 * t726 + (-t662 * t726 + t672 * t729) * qJD(5), 0, t649 * t729 - t657 * t755, -t650 * t726 + t668 * t729 + (-t658 * t729 - t666 * t726) * qJD(5), 0; 0, t661 * t729 + t692 * t726 + (-t674 * t726 + t694 * t729) * qJD(5), 0, t655 * t729 - t664 * t755, -t656 * t726 + t693 * t729 + (-t665 * t729 - t675 * t726) * qJD(5), 0; 0, -t654 * t726 + t671 * t729 + (-t663 * t729 - t673 * t726) * qJD(5), 0, -t651 * t726 - t659 * t754, -t652 * t729 - t670 * t726 + (t660 * t726 - t667 * t729) * qJD(5), 0; 0, -t653 * t726 + t669 * t729 + (-t662 * t729 - t672 * t726) * qJD(5), 0, -t649 * t726 - t657 * t754, -t650 * t729 - t668 * t726 + (t658 * t726 - t666 * t729) * qJD(5), 0; 0, -t661 * t726 + t692 * t729 + (-t674 * t729 - t694 * t726) * qJD(5), 0, -t655 * t726 - t664 * t754, -t656 * t729 - t693 * t726 + (t665 * t726 - t675 * t729) * qJD(5), 0; 0, t663 * qJD(4) + t687 * t727 + t742 * t730, 0, t652, 0, 0; 0, t662 * qJD(4) + t683 * t727 + t744 * t730, 0, t650, 0, 0; 0, t674 * qJD(4) + t700 * t727 - t733 * t730, 0, t656, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:17
	% EndTime: 2019-10-09 22:05:19
	% DurationCPUTime: 2.44s
	% Computational Cost: add. (2096->191), mult. (6837->368), div. (0->0), fcn. (8169->18), ass. (0->154)
	t973 = cos(pkin(6));
	t981 = cos(qJ(2));
	t1033 = t973 * t981;
	t965 = sin(pkin(13));
	t970 = cos(pkin(13));
	t977 = sin(qJ(2));
	t1004 = t965 * t1033 + t970 * t977;
	t968 = sin(pkin(6));
	t972 = cos(pkin(7));
	t1039 = t968 * t972;
	t967 = sin(pkin(7));
	t1008 = t1004 * t967 + t965 * t1039;
	t966 = sin(pkin(8));
	t971 = cos(pkin(8));
	t1034 = t973 * t977;
	t1003 = t965 * t1034 - t970 * t981;
	t1043 = t967 * t968;
	t1007 = -t1004 * t972 + t965 * t1043;
	t964 = sin(pkin(14));
	t969 = cos(pkin(14));
	t989 = -t1003 * t964 - t1007 * t969;
	t1057 = -t1008 * t966 + t989 * t971;
	t926 = -t1003 * t969 + t1007 * t964;
	t976 = sin(qJ(4));
	t980 = cos(qJ(4));
	t901 = t1057 * t980 + t926 * t976;
	t957 = t970 * t1033 - t965 * t977;
	t1010 = -t970 * t1039 - t957 * t967;
	t1009 = -t970 * t1043 + t957 * t972;
	t958 = t970 * t1034 + t965 * t981;
	t990 = -t1009 * t969 + t958 * t964;
	t1056 = -t1010 * t966 + t990 * t971;
	t925 = t1009 * t964 + t958 * t969;
	t899 = t1056 * t980 + t925 * t976;
	t1035 = t972 * t981;
	t1037 = t969 * t977;
	t1041 = t967 * t973;
	t945 = t968 * t1037 + (t968 * t1035 + t1041) * t964;
	t1025 = t981 * t1043;
	t1002 = t973 * t972 - t1025;
	t1006 = t969 * t1035 - t964 * t977;
	t988 = t1006 * t968 + t969 * t1041;
	t982 = t1002 * t966 + t988 * t971;
	t910 = t945 * t976 - t980 * t982;
	t1040 = t967 * t977;
	t1023 = t966 * t1040;
	t1021 = t968 * t1023;
	t1036 = t972 * t977;
	t1005 = t969 * t1036 + t964 * t981;
	t951 = t1005 * t968;
	t1001 = -t951 * t971 + t1021;
	t952 = (-t964 * t1036 + t969 * t981) * t968;
	t1055 = t1001 * t980 - t952 * t976;
	t1044 = t966 * t967;
	t1038 = t969 * t972;
	t937 = t1003 * t1038 + t1004 * t964;
	t1013 = -t1003 * t1044 + t937 * t971;
	t1046 = t964 * t972;
	t938 = t1003 * t1046 - t1004 * t969;
	t1054 = t1013 * t980 - t938 * t976;
	t935 = -t958 * t1038 - t957 * t964;
	t1014 = t958 * t1044 + t935 * t971;
	t936 = -t958 * t1046 + t957 * t969;
	t1053 = t1014 * t980 - t936 * t976;
	t1042 = t967 * t971;
	t1032 = qJD(2) * t968;
	t975 = sin(qJ(5));
	t1031 = qJD(5) * t975;
	t979 = cos(qJ(5));
	t1030 = qJD(5) * t979;
	t974 = sin(qJ(6));
	t1029 = qJD(6) * t974;
	t978 = cos(qJ(6));
	t1028 = qJD(6) * t978;
	t1027 = qJD(6) * t979;
	t1026 = t980 * t1044;
	t1022 = t971 * t1040;
	t1020 = qJD(2) * t1025;
	t953 = t957 * qJD(2);
	t954 = t958 * qJD(2);
	t1012 = -t954 * t1038 - t953 * t964;
	t1000 = t1012 * t971;
	t929 = -t954 * t1046 + t953 * t969;
	t882 = t929 * t980 + (t954 * t1044 + t1000) * t976 - t899 * qJD(4);
	t1019 = t899 * t1027 + t882;
	t955 = t1004 * qJD(2);
	t956 = t1003 * qJD(2);
	t932 = t956 * t1046 - t955 * t969;
	t1011 = t956 * t1038 + t955 * t964;
	t999 = t1011 * t971;
	t884 = t932 * t980 + (-t956 * t1044 + t999) * t976 - t901 * qJD(4);
	t1018 = t901 * t1027 + t884;
	t949 = qJD(2) * t952;
	t996 = t1005 * t971;
	t896 = t949 * t980 - t910 * qJD(4) + (-t996 + t1023) * t976 * t1032;
	t1017 = t910 * t1027 + t896;
	t900 = -t1056 * t976 + t925 * t980;
	t912 = t1010 * t971 + t966 * t990;
	t890 = t900 * t979 + t912 * t975;
	t889 = -t900 * t975 + t912 * t979;
	t902 = -t1057 * t976 + t926 * t980;
	t913 = t1008 * t971 + t966 * t989;
	t892 = t902 * t979 + t913 * t975;
	t891 = -t902 * t975 + t913 * t979;
	t906 = t1014 * t976 + t936 * t980;
	t918 = t958 * t1042 - t935 * t966;
	t893 = t906 * t979 + t918 * t975;
	t908 = t1013 * t976 + t938 * t980;
	t919 = -t1003 * t1042 - t937 * t966;
	t894 = t908 * t979 + t919 * t975;
	t911 = t945 * t980 + t976 * t982;
	t924 = t1002 * t971 - t966 * t988;
	t898 = t911 * t979 + t924 * t975;
	t897 = -t911 * t975 + t924 * t979;
	t921 = t1001 * t976 + t952 * t980;
	t941 = t968 * t1022 + t951 * t966;
	t909 = t921 * t979 + t941 * t975;
	t930 = -t953 * t1038 + t954 * t964;
	t1016 = -t953 * t1044 - t930 * t971;
	t933 = t955 * t1038 - t956 * t964;
	t1015 = t955 * t1044 - t933 * t971;
	t947 = t1006 * t1032;
	t994 = t966 * t1020 - t947 * t971;
	t881 = qJD(4) * t900 - t1000 * t980 - t954 * t1026 + t929 * t976;
	t993 = qJD(6) * t900 + t899 * t1031 - t881 * t979;
	t883 = qJD(4) * t902 + t956 * t1026 + t932 * t976 - t980 * t999;
	t992 = qJD(6) * t902 + t901 * t1031 - t883 * t979;
	t895 = qJD(4) * t911 + t949 * t976 + (-qJD(2) * t1021 + t1032 * t996) * t980;
	t991 = qJD(6) * t911 + t910 * t1031 - t895 * t979;
	t948 = (-t964 * t1035 - t1037) * t1032;
	t940 = (t1005 * t966 + t1022) * t1032;
	t939 = t971 * t1020 + t947 * t966;
	t934 = t955 * t1046 + t956 * t969;
	t931 = -t953 * t1046 - t954 * t969;
	t917 = -t955 * t1042 - t933 * t966;
	t916 = -t1011 * t966 - t956 * t1042;
	t915 = t953 * t1042 - t930 * t966;
	t914 = -t1012 * t966 + t954 * t1042;
	t904 = t1055 * qJD(4) + t948 * t980 + t994 * t976;
	t903 = qJD(4) * t921 + t948 * t976 - t980 * t994;
	t888 = t1054 * qJD(4) - t1015 * t976 + t934 * t980;
	t887 = qJD(4) * t908 + t1015 * t980 + t934 * t976;
	t886 = t1053 * qJD(4) - t1016 * t976 + t931 * t980;
	t885 = qJD(4) * t906 + t1016 * t980 + t931 * t976;
	t880 = t904 * t979 + t939 * t975 + (-t921 * t975 + t941 * t979) * qJD(5);
	t879 = qJD(5) * t897 + t896 * t979 + t940 * t975;
	t878 = -qJD(5) * t898 - t896 * t975 + t940 * t979;
	t877 = t888 * t979 + t917 * t975 + (-t908 * t975 + t919 * t979) * qJD(5);
	t876 = t886 * t979 + t915 * t975 + (-t906 * t975 + t918 * t979) * qJD(5);
	t875 = t891 * qJD(5) + t884 * t979 + t916 * t975;
	t874 = -t892 * qJD(5) - t884 * t975 + t916 * t979;
	t873 = t889 * qJD(5) + t882 * t979 + t914 * t975;
	t872 = -t890 * qJD(5) - t882 * t975 + t914 * t979;
	t1 = [0, t877 * t978 + t887 * t974 + (-t1054 * t978 - t894 * t974) * qJD(6), 0, t1018 * t974 + t992 * t978, -t891 * t1029 + t874 * t978, -t875 * t974 + t883 * t978 + (-t892 * t978 - t901 * t974) * qJD(6); 0, t876 * t978 + t885 * t974 + (-t1053 * t978 - t893 * t974) * qJD(6), 0, t1019 * t974 + t993 * t978, -t889 * t1029 + t872 * t978, -t873 * t974 + t881 * t978 + (-t890 * t978 - t899 * t974) * qJD(6); 0, t880 * t978 + t903 * t974 + (-t1055 * t978 - t909 * t974) * qJD(6), 0, t1017 * t974 + t991 * t978, -t897 * t1029 + t878 * t978, -t879 * t974 + t895 * t978 + (-t898 * t978 - t910 * t974) * qJD(6); 0, -t877 * t974 + t887 * t978 + (t1054 * t974 - t894 * t978) * qJD(6), 0, t1018 * t978 - t992 * t974, -t891 * t1028 - t874 * t974, -t875 * t978 - t883 * t974 + (t892 * t974 - t901 * t978) * qJD(6); 0, -t876 * t974 + t885 * t978 + (t1053 * t974 - t893 * t978) * qJD(6), 0, t1019 * t978 - t993 * t974, -t889 * t1028 - t872 * t974, -t873 * t978 - t881 * t974 + (t890 * t974 - t899 * t978) * qJD(6); 0, -t880 * t974 + t903 * t978 + (t1055 * t974 - t909 * t978) * qJD(6), 0, t1017 * t978 - t991 * t974, -t897 * t1028 - t878 * t974, -t879 * t978 - t895 * t974 + (t898 * t974 - t910 * t978) * qJD(6); 0, qJD(5) * t894 + t888 * t975 - t917 * t979, 0, -t901 * t1030 - t883 * t975, t875, 0; 0, t893 * qJD(5) + t886 * t975 - t915 * t979, 0, -t899 * t1030 - t881 * t975, t873, 0; 0, qJD(5) * t909 + t904 * t975 - t939 * t979, 0, -t910 * t1030 - t895 * t975, t879, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end