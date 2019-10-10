% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PPPRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiRD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (62->26), mult. (202->64), div. (0->0), fcn. (252->14), ass. (0->33)
	t220 = sin(pkin(12));
	t229 = cos(pkin(6));
	t242 = t220 * t229;
	t222 = sin(pkin(7));
	t223 = sin(pkin(6));
	t241 = t222 * t223;
	t240 = t222 * t229;
	t228 = cos(pkin(7));
	t239 = t223 * t228;
	t225 = cos(pkin(13));
	t238 = t225 * t228;
	t226 = cos(pkin(12));
	t237 = t226 * t229;
	t219 = sin(pkin(13));
	t214 = -t220 * t219 + t225 * t237;
	t215 = t219 * t237 + t220 * t225;
	t218 = sin(pkin(14));
	t221 = sin(pkin(8));
	t224 = cos(pkin(14));
	t227 = cos(pkin(8));
	t233 = t214 * t228 - t226 * t241;
	t236 = -(-t215 * t218 + t233 * t224) * t227 - (-t214 * t222 - t226 * t239) * t221;
	t216 = -t226 * t219 - t225 * t242;
	t217 = -t219 * t242 + t226 * t225;
	t232 = t216 * t228 + t220 * t241;
	t235 = -(-t217 * t218 + t232 * t224) * t227 - (-t216 * t222 + t220 * t239) * t221;
	t234 = -(t224 * t240 + (-t218 * t219 + t224 * t238) * t223) * t227 - (-t225 * t241 + t229 * t228) * t221;
	t231 = cos(qJ(4));
	t230 = sin(qJ(4));
	t210 = t223 * t219 * t224 + (t223 * t238 + t240) * t218;
	t208 = t217 * t224 + t232 * t218;
	t206 = t215 * t224 + t233 * t218;
	t1 = [0, 0, 0, (-t208 * t231 + t235 * t230) * qJD(4), 0, 0; 0, 0, 0, (-t206 * t231 + t236 * t230) * qJD(4), 0, 0; 0, 0, 0, (-t210 * t231 + t234 * t230) * qJD(4), 0, 0; 0, 0, 0, (t208 * t230 + t235 * t231) * qJD(4), 0, 0; 0, 0, 0, (t206 * t230 + t236 * t231) * qJD(4), 0, 0; 0, 0, 0, (t210 * t230 + t234 * t231) * qJD(4), 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:24
	% EndTime: 2019-10-10 08:49:24
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (333->47), mult. (1025->108), div. (0->0), fcn. (1310->16), ass. (0->58)
	t457 = sin(pkin(14));
	t458 = sin(pkin(13));
	t462 = sin(pkin(6));
	t463 = cos(pkin(14));
	t464 = cos(pkin(13));
	t467 = cos(pkin(7));
	t481 = t464 * t467;
	t461 = sin(pkin(7));
	t468 = cos(pkin(6));
	t483 = t461 * t468;
	t449 = t462 * t458 * t463 + (t462 * t481 + t483) * t457;
	t470 = sin(qJ(4));
	t472 = cos(qJ(4));
	t448 = t463 * t483 + (-t457 * t458 + t463 * t481) * t462;
	t484 = t461 * t462;
	t452 = -t464 * t484 + t468 * t467;
	t460 = sin(pkin(8));
	t466 = cos(pkin(8));
	t475 = t448 * t466 + t452 * t460;
	t440 = t449 * t472 + t475 * t470;
	t465 = cos(pkin(12));
	t459 = sin(pkin(12));
	t485 = t459 * t468;
	t456 = -t458 * t485 + t465 * t464;
	t455 = -t465 * t458 - t464 * t485;
	t473 = t455 * t467 + t459 * t484;
	t446 = t456 * t463 + t473 * t457;
	t445 = -t456 * t457 + t473 * t463;
	t482 = t462 * t467;
	t451 = -t455 * t461 + t459 * t482;
	t476 = t445 * t466 + t451 * t460;
	t436 = t446 * t472 + t476 * t470;
	t480 = t465 * t468;
	t454 = t458 * t480 + t459 * t464;
	t453 = -t459 * t458 + t464 * t480;
	t474 = t453 * t467 - t465 * t484;
	t444 = t454 * t463 + t474 * t457;
	t443 = -t454 * t457 + t474 * t463;
	t450 = -t453 * t461 - t465 * t482;
	t477 = t443 * t466 + t450 * t460;
	t434 = t444 * t472 + t477 * t470;
	t469 = sin(qJ(5));
	t479 = qJD(5) * t469;
	t471 = cos(qJ(5));
	t478 = qJD(5) * t471;
	t433 = -t444 * t470 + t477 * t472;
	t435 = -t446 * t470 + t476 * t472;
	t439 = -t449 * t470 + t475 * t472;
	t447 = -t448 * t460 + t452 * t466;
	t442 = -t445 * t460 + t451 * t466;
	t441 = -t443 * t460 + t450 * t466;
	t438 = t440 * qJD(4);
	t437 = t439 * qJD(4);
	t432 = t436 * qJD(4);
	t431 = t435 * qJD(4);
	t430 = t434 * qJD(4);
	t429 = t433 * qJD(4);
	t1 = [0, 0, 0, -t432 * t471 - t435 * t479, -t431 * t469 + (-t436 * t471 - t442 * t469) * qJD(5), 0; 0, 0, 0, -t430 * t471 - t433 * t479, -t429 * t469 + (-t434 * t471 - t441 * t469) * qJD(5), 0; 0, 0, 0, -t438 * t471 - t439 * t479, -t437 * t469 + (-t440 * t471 - t447 * t469) * qJD(5), 0; 0, 0, 0, t432 * t469 - t435 * t478, -t431 * t471 + (t436 * t469 - t442 * t471) * qJD(5), 0; 0, 0, 0, t430 * t469 - t433 * t478, -t429 * t471 + (t434 * t469 - t441 * t471) * qJD(5), 0; 0, 0, 0, t438 * t469 - t439 * t478, -t437 * t471 + (t440 * t469 - t447 * t471) * qJD(5), 0; 0, 0, 0, t431, 0, 0; 0, 0, 0, t429, 0, 0; 0, 0, 0, t437, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:29
	% EndTime: 2019-10-10 08:49:29
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (1185->83), mult. (3564->171), div. (0->0), fcn. (4596->18), ass. (0->81)
	t648 = sin(pkin(8));
	t654 = cos(pkin(8));
	t645 = sin(pkin(14));
	t646 = sin(pkin(13));
	t650 = sin(pkin(6));
	t651 = cos(pkin(14));
	t652 = cos(pkin(13));
	t655 = cos(pkin(7));
	t688 = t652 * t655;
	t649 = sin(pkin(7));
	t656 = cos(pkin(6));
	t690 = t649 * t656;
	t665 = t651 * t690 + (-t645 * t646 + t651 * t688) * t650;
	t691 = t649 * t650;
	t674 = -t652 * t691 + t656 * t655;
	t696 = t674 * t648 + t665 * t654;
	t653 = cos(pkin(12));
	t647 = sin(pkin(12));
	t692 = t647 * t656;
	t644 = -t646 * t692 + t653 * t652;
	t643 = -t653 * t646 - t652 * t692;
	t675 = t643 * t655 + t647 * t691;
	t666 = t644 * t645 - t675 * t651;
	t689 = t650 * t655;
	t676 = -t643 * t649 + t647 * t689;
	t695 = -t676 * t648 + t666 * t654;
	t687 = t653 * t656;
	t642 = t646 * t687 + t647 * t652;
	t641 = -t647 * t646 + t652 * t687;
	t677 = t641 * t655 - t653 * t691;
	t667 = t642 * t645 - t677 * t651;
	t678 = -t641 * t649 - t653 * t689;
	t694 = -t678 * t648 + t667 * t654;
	t693 = cos(qJ(4));
	t658 = sin(qJ(5));
	t686 = qJD(5) * t658;
	t661 = cos(qJ(5));
	t685 = qJD(5) * t661;
	t657 = sin(qJ(6));
	t684 = qJD(6) * t657;
	t660 = cos(qJ(6));
	t683 = qJD(6) * t660;
	t682 = qJD(6) * t661;
	t633 = t642 * t651 + t677 * t645;
	t659 = sin(qJ(4));
	t621 = t633 * t659 + t694 * t693;
	t615 = t621 * qJD(4);
	t681 = t621 * t682 - t615;
	t634 = t644 * t651 + t675 * t645;
	t623 = t634 * t659 + t695 * t693;
	t617 = t623 * qJD(4);
	t680 = t623 * t682 - t617;
	t639 = t650 * t646 * t651 + (t650 * t688 + t690) * t645;
	t627 = t639 * t659 - t696 * t693;
	t625 = t627 * qJD(4);
	t679 = t627 * t682 - t625;
	t622 = t633 * t693 - t694 * t659;
	t629 = t667 * t648 + t678 * t654;
	t612 = t622 * t661 + t629 * t658;
	t611 = -t622 * t658 + t629 * t661;
	t624 = t634 * t693 - t695 * t659;
	t630 = t666 * t648 + t676 * t654;
	t614 = t624 * t661 + t630 * t658;
	t613 = -t624 * t658 + t630 * t661;
	t628 = t639 * t693 + t696 * t659;
	t635 = -t665 * t648 + t674 * t654;
	t620 = t628 * t661 + t635 * t658;
	t619 = -t628 * t658 + t635 * t661;
	t616 = t622 * qJD(4);
	t670 = qJD(6) * t622 - t616 * t661 + t621 * t686;
	t618 = t624 * qJD(4);
	t669 = qJD(6) * t624 - t618 * t661 + t623 * t686;
	t626 = t628 * qJD(4);
	t668 = qJD(6) * t628 - t626 * t661 + t627 * t686;
	t610 = t619 * qJD(5) - t625 * t661;
	t609 = -t620 * qJD(5) + t625 * t658;
	t608 = t613 * qJD(5) - t617 * t661;
	t607 = -t614 * qJD(5) + t617 * t658;
	t606 = t611 * qJD(5) - t615 * t661;
	t605 = -t612 * qJD(5) + t615 * t658;
	t1 = [0, 0, 0, t680 * t657 + t669 * t660, t607 * t660 - t613 * t684, -t608 * t657 + t618 * t660 + (-t614 * t660 - t623 * t657) * qJD(6); 0, 0, 0, t681 * t657 + t670 * t660, t605 * t660 - t611 * t684, -t606 * t657 + t616 * t660 + (-t612 * t660 - t621 * t657) * qJD(6); 0, 0, 0, t679 * t657 + t668 * t660, t609 * t660 - t619 * t684, -t610 * t657 + t626 * t660 + (-t620 * t660 - t627 * t657) * qJD(6); 0, 0, 0, -t669 * t657 + t680 * t660, -t607 * t657 - t613 * t683, -t608 * t660 - t618 * t657 + (t614 * t657 - t623 * t660) * qJD(6); 0, 0, 0, -t670 * t657 + t681 * t660, -t605 * t657 - t611 * t683, -t606 * t660 - t616 * t657 + (t612 * t657 - t621 * t660) * qJD(6); 0, 0, 0, -t668 * t657 + t679 * t660, -t609 * t657 - t619 * t683, -t610 * t660 - t626 * t657 + (t620 * t657 - t627 * t660) * qJD(6); 0, 0, 0, -t618 * t658 - t623 * t685, t608, 0; 0, 0, 0, -t616 * t658 - t621 * t685, t606, 0; 0, 0, 0, -t626 * t658 - t627 * t685, t610, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end