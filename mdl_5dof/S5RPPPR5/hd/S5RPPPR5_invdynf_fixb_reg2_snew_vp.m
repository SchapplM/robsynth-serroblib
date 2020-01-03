% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPPPR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:46
% EndTime: 2019-12-31 17:46:47
% DurationCPUTime: 1.51s
% Computational Cost: add. (4123->179), mult. (8094->197), div. (0->0), fcn. (4299->8), ass. (0->108)
t864 = qJD(1) ^ 2;
t855 = sin(pkin(8));
t850 = t855 ^ 2;
t857 = cos(pkin(8));
t851 = t857 ^ 2;
t878 = t850 + t851;
t832 = t878 * t864;
t859 = sin(qJ(5));
t861 = cos(qJ(5));
t823 = (t855 * t859 - t857 * t861) * qJD(1);
t887 = t823 ^ 2;
t871 = t855 * t861 + t857 * t859;
t824 = t871 * qJD(1);
t886 = t824 ^ 2;
t885 = -2 * qJD(5);
t884 = -pkin(1) - pkin(2);
t883 = t823 * t824;
t882 = t850 * t864;
t881 = t851 * t864;
t880 = t857 * t864;
t860 = sin(qJ(1));
t862 = cos(qJ(1));
t836 = -t862 * g(1) - t860 * g(2);
t869 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t836;
t816 = t884 * t864 + t869;
t856 = sin(pkin(7));
t858 = cos(pkin(7));
t835 = t860 * g(1) - t862 * g(2);
t868 = -t864 * qJ(2) + qJDD(2) - t835;
t867 = t884 * qJDD(1) + t868;
t800 = t858 * t816 + t856 * t867;
t854 = g(3) + qJDD(3);
t877 = qJD(1) * qJD(4);
t879 = t857 * t854 + 0.2e1 * t855 * t877;
t847 = t855 * qJDD(1);
t876 = t856 * qJDD(1);
t875 = t857 * qJDD(1);
t874 = t858 * qJDD(1);
t797 = -t864 * pkin(3) - qJDD(1) * qJ(4) + t800;
t792 = t855 * t854 + (t797 - 0.2e1 * t877) * t857;
t799 = -t856 * t816 + t858 * t867;
t830 = t858 * t864 - t876;
t831 = t856 * t864 + t874;
t873 = -t860 * t830 + t862 * t831;
t798 = t859 * t847 - t861 * t875;
t872 = t862 * t830 + t860 * t831;
t870 = -t864 * qJ(4) + qJDD(4) - t799;
t822 = t871 * qJDD(1);
t863 = qJD(5) ^ 2;
t838 = t855 * t880;
t834 = t862 * qJDD(1) - t860 * t864;
t833 = t860 * qJDD(1) + t862 * t864;
t829 = t878 * qJDD(1);
t828 = t857 * t832;
t827 = t855 * t832;
t821 = qJDD(1) * pkin(1) - t868;
t818 = -t864 * pkin(1) + t869;
t817 = -t863 - t886;
t812 = -t858 * t828 + t856 * t875;
t811 = t858 * t827 - t855 * t876;
t810 = -t856 * t828 - t857 * t874;
t809 = t856 * t827 + t855 * t874;
t808 = -t858 * t829 - t856 * t832;
t807 = -t856 * t829 + t858 * t832;
t806 = t823 * t885 + t822;
t805 = -t824 * t885 + t798;
t804 = -qJDD(5) - t883;
t803 = qJDD(5) - t883;
t802 = -t863 - t887;
t801 = -t886 - t887;
t796 = qJDD(1) * pkin(3) + t870;
t794 = t861 * t804 - t859 * t817;
t793 = t859 * t804 + t861 * t817;
t791 = -t855 * t797 + t879;
t790 = (pkin(4) * t857 + pkin(3)) * qJDD(1) + t870 + (-t882 - t881) * pkin(6);
t789 = t861 * t798 - t859 * t822;
t788 = t859 * t798 + t861 * t822;
t787 = t861 * t802 - t859 * t803;
t786 = t859 * t802 + t861 * t803;
t785 = -pkin(4) * t881 - pkin(6) * t875 + t792;
t784 = (pkin(4) * t880 + pkin(6) * qJDD(1) - t797) * t855 + t879;
t783 = -t856 * t799 + t858 * t800;
t782 = t858 * t799 + t856 * t800;
t781 = -t855 * t793 + t857 * t794;
t780 = t857 * t793 + t855 * t794;
t779 = -t855 * t791 + t857 * t792;
t778 = t857 * t791 + t855 * t792;
t777 = -t855 * t788 + t857 * t789;
t776 = t857 * t788 + t855 * t789;
t775 = -t855 * t786 + t857 * t787;
t774 = t857 * t786 + t855 * t787;
t773 = t858 * t781 - t856 * t806;
t772 = t856 * t781 + t858 * t806;
t771 = t859 * t784 + t861 * t785;
t770 = t861 * t784 - t859 * t785;
t769 = t858 * t775 - t856 * t805;
t768 = t856 * t775 + t858 * t805;
t767 = t858 * t777 + t856 * t801;
t766 = t856 * t777 - t858 * t801;
t765 = t858 * t779 + t856 * t796;
t764 = t856 * t779 - t858 * t796;
t763 = -t859 * t770 + t861 * t771;
t762 = t861 * t770 + t859 * t771;
t761 = -t855 * t762 + t857 * t763;
t760 = t857 * t762 + t855 * t763;
t759 = t858 * t761 + t856 * t790;
t758 = t856 * t761 - t858 * t790;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t833, -t834, 0, -t860 * t835 + t862 * t836, 0, 0, 0, 0, 0, 0, -t833, 0, t834, t862 * t818 - t860 * t821, 0, 0, 0, 0, 0, 0, -t872, t873, 0, t860 * t782 + t862 * t783, 0, 0, 0, 0, 0, 0, t860 * t810 + t862 * t812, t860 * t809 + t862 * t811, t860 * t807 + t862 * t808, t860 * t764 + t862 * t765, 0, 0, 0, 0, 0, 0, t860 * t768 + t862 * t769, t860 * t772 + t862 * t773, t860 * t766 + t862 * t767, t860 * t758 + t862 * t759; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t834, -t833, 0, t862 * t835 + t860 * t836, 0, 0, 0, 0, 0, 0, t834, 0, t833, t860 * t818 + t862 * t821, 0, 0, 0, 0, 0, 0, t873, t872, 0, -t862 * t782 + t860 * t783, 0, 0, 0, 0, 0, 0, -t862 * t810 + t860 * t812, -t862 * t809 + t860 * t811, -t862 * t807 + t860 * t808, -t862 * t764 + t860 * t765, 0, 0, 0, 0, 0, 0, -t862 * t768 + t860 * t769, -t862 * t772 + t860 * t773, -t862 * t766 + t860 * t767, -t862 * t758 + t860 * t759; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t854, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t778, 0, 0, 0, 0, 0, 0, -t774, -t780, -t776, -t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t864, -qJDD(1), 0, t836, 0, 0, 0, 0, 0, 0, -t864, 0, qJDD(1), t818, 0, 0, 0, 0, 0, 0, -t830, t831, 0, t783, 0, 0, 0, 0, 0, 0, t812, t811, t808, t765, 0, 0, 0, 0, 0, 0, t769, t773, t767, t759; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t864, 0, t835, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t864, t821, 0, 0, 0, 0, 0, 0, t831, t830, 0, -t782, 0, 0, 0, 0, 0, 0, -t810, -t809, -t807, -t764, 0, 0, 0, 0, 0, 0, -t768, -t772, -t766, -t758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t854, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t778, 0, 0, 0, 0, 0, 0, -t774, -t780, -t776, -t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t864, 0, qJDD(1), t818, 0, 0, 0, 0, 0, 0, -t830, t831, 0, t783, 0, 0, 0, 0, 0, 0, t812, t811, t808, t765, 0, 0, 0, 0, 0, 0, t769, t773, t767, t759; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t854, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t778, 0, 0, 0, 0, 0, 0, -t774, -t780, -t776, -t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t864, -t821, 0, 0, 0, 0, 0, 0, -t831, -t830, 0, t782, 0, 0, 0, 0, 0, 0, t810, t809, t807, t764, 0, 0, 0, 0, 0, 0, t768, t772, t766, t758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t864, qJDD(1), 0, t800, 0, 0, 0, 0, 0, 0, -t828, t827, -t829, t779, 0, 0, 0, 0, 0, 0, t775, t781, t777, t761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t864, 0, t799, 0, 0, 0, 0, 0, 0, -t875, t847, t832, -t796, 0, 0, 0, 0, 0, 0, t805, t806, -t801, -t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t854, 0, 0, 0, 0, 0, 0, 0, 0, 0, t778, 0, 0, 0, 0, 0, 0, t774, t780, t776, t760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t881, t838, -t875, t792, 0, 0, 0, 0, 0, 0, t787, t794, t789, t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t838, -t882, t847, t791, 0, 0, 0, 0, 0, 0, t786, t793, t788, t762; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t875, -t847, -t832, t796, 0, 0, 0, 0, 0, 0, -t805, -t806, t801, t790; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t802, t804, t798, t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t803, t817, t822, t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t805, -t806, t801, t790;];
f_new_reg = t1;
