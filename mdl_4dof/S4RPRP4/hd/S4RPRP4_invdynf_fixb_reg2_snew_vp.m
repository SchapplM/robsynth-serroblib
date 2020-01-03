% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP4
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:02
% EndTime: 2019-12-31 16:44:03
% DurationCPUTime: 0.92s
% Computational Cost: add. (1119->115), mult. (2329->131), div. (0->0), fcn. (1283->6), ass. (0->73)
t773 = sin(qJ(3));
t775 = cos(qJ(3));
t778 = qJD(1) ^ 2;
t757 = t775 * t778 * t773;
t751 = qJDD(3) - t757;
t765 = t773 ^ 2;
t777 = qJD(3) ^ 2;
t754 = t765 * t778 + t777;
t732 = t775 * t751 - t773 * t754;
t785 = t773 * qJDD(1);
t786 = qJD(1) * qJD(3);
t742 = 0.2e1 * t775 * t786 + t785;
t770 = sin(pkin(6));
t771 = cos(pkin(6));
t716 = t770 * t732 + t771 * t742;
t719 = t771 * t732 - t770 * t742;
t774 = sin(qJ(1));
t776 = cos(qJ(1));
t792 = t776 * t716 + t774 * t719;
t791 = t774 * t716 - t776 * t719;
t753 = -t776 * g(1) - t774 * g(2);
t741 = -t778 * pkin(1) + t753;
t752 = t774 * g(1) - t776 * g(2);
t779 = qJDD(1) * pkin(1) + t752;
t725 = t771 * t741 + t770 * t779;
t723 = -t778 * pkin(2) + qJDD(1) * pkin(5) + t725;
t767 = -g(3) + qJDD(2);
t714 = t775 * t723 + t773 * t767;
t766 = t775 ^ 2;
t788 = t765 + t766;
t787 = t778 * (-pkin(3) * t775 - qJ(4) * t773);
t784 = t775 * qJDD(1);
t783 = t773 * t786;
t724 = -t770 * t741 + t771 * t779;
t744 = -t770 * qJDD(1) - t771 * t778;
t745 = t771 * qJDD(1) - t770 * t778;
t781 = t776 * t744 - t774 * t745;
t780 = t774 * t744 + t776 * t745;
t729 = t773 * t751 + t775 * t754;
t722 = -qJDD(1) * pkin(2) - t778 * pkin(5) - t724;
t762 = t775 * t767;
t755 = -t766 * t778 - t777;
t750 = qJDD(3) + t757;
t749 = t788 * t778;
t748 = -t774 * qJDD(1) - t776 * t778;
t747 = t776 * qJDD(1) - t774 * t778;
t746 = t788 * qJDD(1);
t743 = -0.2e1 * t783 + t784;
t731 = -t773 * t750 + t775 * t755;
t728 = t775 * t750 + t773 * t755;
t727 = t771 * t746 - t770 * t749;
t726 = t770 * t746 + t771 * t749;
t718 = t771 * t731 - t770 * t743;
t715 = t770 * t731 + t771 * t743;
t713 = -t773 * t723 + t762;
t712 = -t774 * t726 + t776 * t727;
t711 = t776 * t726 + t774 * t727;
t710 = -t770 * t724 + t771 * t725;
t709 = t771 * t724 + t770 * t725;
t708 = qJDD(4) - t762 - t777 * qJ(4) - qJDD(3) * pkin(3) + (t723 + t787) * t773;
t707 = -t777 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t775 * t787 + t714;
t706 = -(-t783 + t784) * pkin(3) + (pkin(3) * qJD(3) - (2 * qJD(4))) * t773 * qJD(1) + t722 - t742 * qJ(4);
t705 = -t774 * t715 + t776 * t718;
t704 = t776 * t715 + t774 * t718;
t703 = -t773 * t713 + t775 * t714;
t702 = t775 * t713 + t773 * t714;
t701 = t771 * t703 + t770 * t722;
t700 = t770 * t703 - t771 * t722;
t699 = t775 * t707 + t773 * t708;
t698 = t773 * t707 - t775 * t708;
t697 = t771 * t699 + t770 * t706;
t696 = t770 * t699 - t771 * t706;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t748, -t747, 0, -t774 * t752 + t776 * t753, 0, 0, 0, 0, 0, 0, t781, -t780, 0, -t774 * t709 + t776 * t710, 0, 0, 0, 0, 0, 0, t705, t791, t712, -t774 * t700 + t776 * t701, 0, 0, 0, 0, 0, 0, t705, t712, -t791, -t774 * t696 + t776 * t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t747, t748, 0, t776 * t752 + t774 * t753, 0, 0, 0, 0, 0, 0, t780, t781, 0, t776 * t709 + t774 * t710, 0, 0, 0, 0, 0, 0, t704, -t792, t711, t776 * t700 + t774 * t701, 0, 0, 0, 0, 0, 0, t704, t711, t792, t776 * t696 + t774 * t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t767, 0, 0, 0, 0, 0, 0, t728, -t729, 0, t702, 0, 0, 0, 0, 0, 0, t728, 0, t729, t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t778, -qJDD(1), 0, t753, 0, 0, 0, 0, 0, 0, t744, -t745, 0, t710, 0, 0, 0, 0, 0, 0, t718, -t719, t727, t701, 0, 0, 0, 0, 0, 0, t718, t727, t719, t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t778, 0, t752, 0, 0, 0, 0, 0, 0, t745, t744, 0, t709, 0, 0, 0, 0, 0, 0, t715, -t716, t726, t700, 0, 0, 0, 0, 0, 0, t715, t726, t716, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t767, 0, 0, 0, 0, 0, 0, t728, -t729, 0, t702, 0, 0, 0, 0, 0, 0, t728, 0, t729, t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t778, -qJDD(1), 0, t725, 0, 0, 0, 0, 0, 0, t731, -t732, t746, t703, 0, 0, 0, 0, 0, 0, t731, t746, t732, t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t778, 0, t724, 0, 0, 0, 0, 0, 0, t743, -t742, t749, -t722, 0, 0, 0, 0, 0, 0, t743, t749, t742, -t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t767, 0, 0, 0, 0, 0, 0, t728, -t729, 0, t702, 0, 0, 0, 0, 0, 0, t728, 0, t729, t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t755, -t751, t784, t714, 0, 0, 0, 0, 0, 0, t755, t784, t751, t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t750, -t754, -t785, t713, 0, 0, 0, 0, 0, 0, t750, -t785, t754, -t708; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t743, t742, -t749, t722, 0, 0, 0, 0, 0, 0, -t743, -t749, -t742, t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t755, t784, t751, t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t743, -t749, -t742, t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t750, t785, -t754, t708;];
f_new_reg = t1;
