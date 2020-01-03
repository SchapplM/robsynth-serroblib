% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:25
% EndTime: 2019-12-31 17:14:26
% DurationCPUTime: 0.92s
% Computational Cost: add. (1589->119), mult. (2329->132), div. (0->0), fcn. (1283->6), ass. (0->75)
t795 = qJD(1) + qJD(2);
t793 = t795 ^ 2;
t800 = sin(qJ(3));
t803 = cos(qJ(3));
t785 = t803 * t793 * t800;
t779 = qJDD(3) - t785;
t797 = t800 ^ 2;
t806 = qJD(3) ^ 2;
t782 = t797 * t793 + t806;
t761 = t803 * t779 - t800 * t782;
t812 = qJD(3) * t795;
t794 = qJDD(1) + qJDD(2);
t817 = t800 * t794;
t766 = 0.2e1 * t803 * t812 + t817;
t801 = sin(qJ(2));
t804 = cos(qJ(2));
t743 = t801 * t761 + t804 * t766;
t746 = t804 * t761 - t801 * t766;
t802 = sin(qJ(1));
t805 = cos(qJ(1));
t825 = t805 * t743 + t802 * t746;
t824 = t802 * t743 - t805 * t746;
t774 = t801 * t793 - t804 * t794;
t809 = -t804 * t793 - t801 * t794;
t823 = t802 * t774 + t805 * t809;
t822 = t805 * t774 - t802 * t809;
t819 = t803 * g(3);
t818 = (-pkin(3) * t803 - qJ(4) * t800) * t793;
t815 = t803 * t794;
t787 = -t805 * g(1) - t802 * g(2);
t807 = qJD(1) ^ 2;
t777 = -t807 * pkin(1) + t787;
t786 = t802 * g(1) - t805 * g(2);
t808 = qJDD(1) * pkin(1) + t786;
t756 = t804 * t777 + t801 * t808;
t798 = t803 ^ 2;
t813 = t797 + t798;
t811 = t800 * t812;
t752 = -t793 * pkin(2) + t794 * pkin(6) + t756;
t749 = -t800 * g(3) + t803 * t752;
t755 = -t801 * t777 + t804 * t808;
t758 = t800 * t779 + t803 * t782;
t751 = -t794 * pkin(2) - t793 * pkin(6) - t755;
t783 = -t798 * t793 - t806;
t781 = -t802 * qJDD(1) - t805 * t807;
t780 = t805 * qJDD(1) - t802 * t807;
t778 = qJDD(3) + t785;
t776 = t813 * t793;
t771 = t813 * t794;
t767 = -0.2e1 * t811 + t815;
t760 = -t800 * t778 + t803 * t783;
t757 = t803 * t778 + t800 * t783;
t754 = t804 * t771 - t801 * t776;
t753 = t801 * t771 + t804 * t776;
t748 = -t800 * t752 - t819;
t745 = t804 * t760 - t801 * t767;
t742 = t801 * t760 + t804 * t767;
t741 = -t801 * t755 + t804 * t756;
t740 = t804 * t755 + t801 * t756;
t739 = -t802 * t753 + t805 * t754;
t738 = t805 * t753 + t802 * t754;
t737 = t819 + qJDD(4) - t806 * qJ(4) - qJDD(3) * pkin(3) + (t752 + t818) * t800;
t736 = -t806 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t803 * t818 + t749;
t735 = -(-t811 + t815) * pkin(3) + (pkin(3) * qJD(3) - (2 * qJD(4))) * t800 * t795 + t751 - t766 * qJ(4);
t734 = -t800 * t748 + t803 * t749;
t733 = t803 * t748 + t800 * t749;
t732 = -t802 * t742 + t805 * t745;
t731 = t805 * t742 + t802 * t745;
t730 = t804 * t734 + t801 * t751;
t729 = t801 * t734 - t804 * t751;
t728 = t803 * t736 + t800 * t737;
t727 = t800 * t736 - t803 * t737;
t726 = t804 * t728 + t801 * t735;
t725 = t801 * t728 - t804 * t735;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t781, -t780, 0, -t802 * t786 + t805 * t787, 0, 0, 0, 0, 0, 0, t823, t822, 0, -t802 * t740 + t805 * t741, 0, 0, 0, 0, 0, 0, t732, t824, t739, -t802 * t729 + t805 * t730, 0, 0, 0, 0, 0, 0, t732, t739, -t824, -t802 * t725 + t805 * t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t780, t781, 0, t805 * t786 + t802 * t787, 0, 0, 0, 0, 0, 0, -t822, t823, 0, t805 * t740 + t802 * t741, 0, 0, 0, 0, 0, 0, t731, -t825, t738, t805 * t729 + t802 * t730, 0, 0, 0, 0, 0, 0, t731, t738, t825, t805 * t725 + t802 * t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t757, -t758, 0, t733, 0, 0, 0, 0, 0, 0, t757, 0, t758, t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t807, -qJDD(1), 0, t787, 0, 0, 0, 0, 0, 0, t809, t774, 0, t741, 0, 0, 0, 0, 0, 0, t745, -t746, t754, t730, 0, 0, 0, 0, 0, 0, t745, t754, t746, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t807, 0, t786, 0, 0, 0, 0, 0, 0, -t774, t809, 0, t740, 0, 0, 0, 0, 0, 0, t742, -t743, t753, t729, 0, 0, 0, 0, 0, 0, t742, t753, t743, t725; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t757, -t758, 0, t733, 0, 0, 0, 0, 0, 0, t757, 0, t758, t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t793, -t794, 0, t756, 0, 0, 0, 0, 0, 0, t760, -t761, t771, t734, 0, 0, 0, 0, 0, 0, t760, t771, t761, t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t794, -t793, 0, t755, 0, 0, 0, 0, 0, 0, t767, -t766, t776, -t751, 0, 0, 0, 0, 0, 0, t767, t776, t766, -t735; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t757, -t758, 0, t733, 0, 0, 0, 0, 0, 0, t757, 0, t758, t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t783, -t779, t815, t749, 0, 0, 0, 0, 0, 0, t783, t815, t779, t736; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t778, -t782, -t817, t748, 0, 0, 0, 0, 0, 0, t778, -t817, t782, -t737; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t767, t766, -t776, t751, 0, 0, 0, 0, 0, 0, -t767, -t776, -t766, t735; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t783, t815, t779, t736; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t767, -t776, -t766, t735; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t778, t817, -t782, t737;];
f_new_reg = t1;
