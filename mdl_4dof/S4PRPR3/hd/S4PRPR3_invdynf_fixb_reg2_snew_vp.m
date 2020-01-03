% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:05
% EndTime: 2019-12-31 16:21:06
% DurationCPUTime: 0.91s
% Computational Cost: add. (2073->116), mult. (4769->181), div. (0->0), fcn. (3389->8), ass. (0->97)
t793 = qJD(2) ^ 2;
t784 = sin(pkin(7));
t779 = t784 ^ 2;
t786 = cos(pkin(7));
t780 = t786 ^ 2;
t804 = t779 + t780;
t763 = t804 * t793;
t788 = sin(qJ(4));
t790 = cos(qJ(4));
t796 = t784 * t790 + t786 * t788;
t752 = t796 * qJDD(2);
t753 = (t784 * t788 - t786 * t790) * qJD(2);
t811 = t753 ^ 2;
t755 = t796 * qJD(2);
t810 = t755 ^ 2;
t809 = t755 * t753;
t808 = t779 * t793;
t807 = t780 * t793;
t806 = t786 * t793;
t785 = sin(pkin(6));
t787 = cos(pkin(6));
t766 = -t787 * g(1) - t785 * g(2);
t789 = sin(qJ(2));
t791 = cos(qJ(2));
t798 = t785 * g(1) - t787 * g(2);
t744 = t791 * t766 + t789 * t798;
t782 = -g(3) + qJDD(1);
t803 = qJD(2) * qJD(3);
t805 = t786 * t782 - 0.2e1 * t784 * t803;
t802 = t784 * qJDD(2);
t777 = t786 * qJDD(2);
t801 = t789 * qJDD(2);
t800 = t791 * qJDD(2);
t740 = -t793 * pkin(2) + qJDD(2) * qJ(3) + t744;
t731 = t784 * t782 + (t740 + 0.2e1 * t803) * t786;
t764 = -t789 * t793 + t800;
t765 = -t791 * t793 - t801;
t799 = -t785 * t764 + t787 * t765;
t743 = -t789 * t766 + t791 * t798;
t729 = t790 * t777 - t788 * t802;
t797 = t787 * t764 + t785 * t765;
t737 = -qJDD(2) * pkin(2) - t793 * qJ(3) + qJDD(3) - t743;
t792 = qJD(4) ^ 2;
t767 = t784 * t806;
t762 = t804 * qJDD(2);
t758 = t786 * t763;
t757 = t784 * t763;
t749 = -t792 - t810;
t748 = -t791 * t758 - t786 * t801;
t747 = t791 * t757 + t784 * t801;
t746 = -t789 * t758 + t786 * t800;
t745 = t789 * t757 - t784 * t800;
t742 = t791 * t762 - t789 * t763;
t741 = t789 * t762 + t791 * t763;
t739 = -0.2e1 * t753 * qJD(4) + t752;
t738 = 0.2e1 * t755 * qJD(4) - t729;
t736 = -qJDD(4) - t809;
t735 = qJDD(4) - t809;
t734 = -t792 - t811;
t732 = -t810 - t811;
t730 = -t784 * t740 + t805;
t728 = -pkin(3) * t777 + t737 + (-t807 - t808) * pkin(5);
t727 = -t789 * t743 + t791 * t744;
t726 = t791 * t743 + t789 * t744;
t725 = t790 * t736 - t788 * t749;
t724 = t788 * t736 + t790 * t749;
t723 = -pkin(3) * t807 + pkin(5) * t777 + t731;
t722 = (pkin(3) * t806 - pkin(5) * qJDD(2) - t740) * t784 + t805;
t721 = t790 * t729 + t788 * t752;
t720 = t788 * t729 - t790 * t752;
t719 = t790 * t734 - t788 * t735;
t718 = t788 * t734 + t790 * t735;
t717 = -t784 * t730 + t786 * t731;
t716 = t786 * t730 + t784 * t731;
t715 = -t784 * t724 + t786 * t725;
t714 = t786 * t724 + t784 * t725;
t713 = t788 * t722 + t790 * t723;
t712 = t790 * t722 - t788 * t723;
t711 = t791 * t717 + t789 * t737;
t710 = t789 * t717 - t791 * t737;
t709 = -t784 * t720 + t786 * t721;
t708 = t786 * t720 + t784 * t721;
t707 = -t784 * t718 + t786 * t719;
t706 = t786 * t718 + t784 * t719;
t705 = t791 * t715 + t789 * t739;
t704 = t789 * t715 - t791 * t739;
t703 = t791 * t707 + t789 * t738;
t702 = t789 * t707 - t791 * t738;
t701 = t791 * t709 + t789 * t732;
t700 = t789 * t709 - t791 * t732;
t699 = -t788 * t712 + t790 * t713;
t698 = t790 * t712 + t788 * t713;
t697 = -t784 * t698 + t786 * t699;
t696 = t786 * t698 + t784 * t699;
t695 = t791 * t697 + t789 * t728;
t694 = t789 * t697 - t791 * t728;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t787 * t766 - t785 * t798, 0, 0, 0, 0, 0, 0, t799, -t797, 0, -t785 * t726 + t787 * t727, 0, 0, 0, 0, 0, 0, -t785 * t746 + t787 * t748, -t785 * t745 + t787 * t747, -t785 * t741 + t787 * t742, -t785 * t710 + t787 * t711, 0, 0, 0, 0, 0, 0, -t785 * t702 + t787 * t703, -t785 * t704 + t787 * t705, -t785 * t700 + t787 * t701, -t785 * t694 + t787 * t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t785 * t766 + t787 * t798, 0, 0, 0, 0, 0, 0, t797, t799, 0, t787 * t726 + t785 * t727, 0, 0, 0, 0, 0, 0, t787 * t746 + t785 * t748, t787 * t745 + t785 * t747, t787 * t741 + t785 * t742, t787 * t710 + t785 * t711, 0, 0, 0, 0, 0, 0, t787 * t702 + t785 * t703, t787 * t704 + t785 * t705, t787 * t700 + t785 * t701, t787 * t694 + t785 * t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t782, 0, 0, 0, 0, 0, 0, 0, 0, 0, t782, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, 0, 0, 0, 0, 0, 0, t706, t714, t708, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t766, 0, 0, 0, 0, 0, 0, t765, -t764, 0, t727, 0, 0, 0, 0, 0, 0, t748, t747, t742, t711, 0, 0, 0, 0, 0, 0, t703, t705, t701, t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t798, 0, 0, 0, 0, 0, 0, t764, t765, 0, t726, 0, 0, 0, 0, 0, 0, t746, t745, t741, t710, 0, 0, 0, 0, 0, 0, t702, t704, t700, t694; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t782, 0, 0, 0, 0, 0, 0, 0, 0, 0, t782, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, 0, 0, 0, 0, 0, 0, t706, t714, t708, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t793, -qJDD(2), 0, t744, 0, 0, 0, 0, 0, 0, -t758, t757, t762, t717, 0, 0, 0, 0, 0, 0, t707, t715, t709, t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t793, 0, t743, 0, 0, 0, 0, 0, 0, t777, -t802, t763, -t737, 0, 0, 0, 0, 0, 0, -t738, -t739, -t732, -t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t782, 0, 0, 0, 0, 0, 0, 0, 0, 0, t716, 0, 0, 0, 0, 0, 0, t706, t714, t708, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t807, t767, t777, t731, 0, 0, 0, 0, 0, 0, t719, t725, t721, t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t767, -t808, -t802, t730, 0, 0, 0, 0, 0, 0, t718, t724, t720, t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t777, t802, -t763, t737, 0, 0, 0, 0, 0, 0, t738, t739, t732, t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t734, t736, t729, t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t735, t749, -t752, t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t738, t739, t732, t728;];
f_new_reg = t1;
