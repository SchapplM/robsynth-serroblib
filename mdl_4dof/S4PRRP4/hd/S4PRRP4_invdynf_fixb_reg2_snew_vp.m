% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRP4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:06
% EndTime: 2019-12-31 16:28:07
% DurationCPUTime: 0.85s
% Computational Cost: add. (958->105), mult. (2084->125), div. (0->0), fcn. (1275->6), ass. (0->69)
t766 = sin(qJ(3));
t768 = cos(qJ(3));
t771 = qJD(2) ^ 2;
t750 = t766 * t771 * t768;
t746 = qJDD(3) - t750;
t758 = t766 ^ 2;
t770 = qJD(3) ^ 2;
t747 = t758 * t771 + t770;
t729 = t746 * t768 - t747 * t766;
t778 = t766 * qJDD(2);
t779 = qJD(2) * qJD(3);
t738 = 0.2e1 * t768 * t779 + t778;
t767 = sin(qJ(2));
t769 = cos(qJ(2));
t712 = t729 * t767 + t738 * t769;
t716 = t729 * t769 - t738 * t767;
t763 = sin(pkin(6));
t764 = cos(pkin(6));
t785 = t712 * t764 + t716 * t763;
t784 = t712 * t763 - t716 * t764;
t744 = -g(1) * t764 - g(2) * t763;
t773 = g(1) * t763 - g(2) * t764;
t724 = t744 * t769 + t767 * t773;
t720 = -pkin(2) * t771 + qJDD(2) * pkin(5) + t724;
t760 = -g(3) + qJDD(1);
t714 = t720 * t768 + t760 * t766;
t759 = t768 ^ 2;
t781 = t758 + t759;
t780 = t771 * (-pkin(3) * t768 - qJ(4) * t766);
t777 = t768 * qJDD(2);
t776 = t766 * t779;
t741 = qJDD(2) * t769 - t767 * t771;
t742 = -qJDD(2) * t767 - t769 * t771;
t774 = -t741 * t763 + t742 * t764;
t723 = -t767 * t744 + t769 * t773;
t772 = t741 * t764 + t742 * t763;
t726 = t746 * t766 + t747 * t768;
t719 = -qJDD(2) * pkin(2) - t771 * pkin(5) - t723;
t755 = t768 * t760;
t748 = -t759 * t771 - t770;
t745 = qJDD(3) + t750;
t743 = t781 * t771;
t740 = t781 * qJDD(2);
t739 = -0.2e1 * t776 + t777;
t728 = -t745 * t766 + t748 * t768;
t725 = t745 * t768 + t748 * t766;
t722 = t740 * t769 - t743 * t767;
t721 = t740 * t767 + t743 * t769;
t715 = t728 * t769 - t739 * t767;
t711 = t728 * t767 + t739 * t769;
t710 = -t720 * t766 + t755;
t709 = -t723 * t767 + t724 * t769;
t708 = t723 * t769 + t724 * t767;
t707 = -t721 * t763 + t722 * t764;
t706 = t721 * t764 + t722 * t763;
t705 = qJDD(4) - t755 - t770 * qJ(4) - qJDD(3) * pkin(3) + (t720 + t780) * t766;
t704 = -pkin(3) * t770 + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t768 * t780 + t714;
t703 = -(-t776 + t777) * pkin(3) + (pkin(3) * qJD(3) - (2 * qJD(4))) * t766 * qJD(2) + t719 - t738 * qJ(4);
t702 = -t710 * t766 + t714 * t768;
t701 = t710 * t768 + t714 * t766;
t700 = -t711 * t763 + t715 * t764;
t699 = t711 * t764 + t715 * t763;
t698 = t702 * t769 + t719 * t767;
t697 = t702 * t767 - t719 * t769;
t696 = t704 * t768 + t705 * t766;
t695 = t704 * t766 - t705 * t768;
t694 = t696 * t769 + t703 * t767;
t693 = t696 * t767 - t703 * t769;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t764 * t744 - t763 * t773, 0, 0, 0, 0, 0, 0, t774, -t772, 0, -t708 * t763 + t709 * t764, 0, 0, 0, 0, 0, 0, t700, t784, t707, -t697 * t763 + t698 * t764, 0, 0, 0, 0, 0, 0, t700, t707, -t784, -t693 * t763 + t694 * t764; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t763 * t744 + t764 * t773, 0, 0, 0, 0, 0, 0, t772, t774, 0, t708 * t764 + t709 * t763, 0, 0, 0, 0, 0, 0, t699, -t785, t706, t697 * t764 + t698 * t763, 0, 0, 0, 0, 0, 0, t699, t706, t785, t693 * t764 + t694 * t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t760, 0, 0, 0, 0, 0, 0, 0, 0, 0, t760, 0, 0, 0, 0, 0, 0, t725, -t726, 0, t701, 0, 0, 0, 0, 0, 0, t725, 0, t726, t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t744, 0, 0, 0, 0, 0, 0, t742, -t741, 0, t709, 0, 0, 0, 0, 0, 0, t715, -t716, t722, t698, 0, 0, 0, 0, 0, 0, t715, t722, t716, t694; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, 0, 0, 0, 0, 0, 0, t741, t742, 0, t708, 0, 0, 0, 0, 0, 0, t711, -t712, t721, t697, 0, 0, 0, 0, 0, 0, t711, t721, t712, t693; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t760, 0, 0, 0, 0, 0, 0, 0, 0, 0, t760, 0, 0, 0, 0, 0, 0, t725, -t726, 0, t701, 0, 0, 0, 0, 0, 0, t725, 0, t726, t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t771, -qJDD(2), 0, t724, 0, 0, 0, 0, 0, 0, t728, -t729, t740, t702, 0, 0, 0, 0, 0, 0, t728, t740, t729, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t771, 0, t723, 0, 0, 0, 0, 0, 0, t739, -t738, t743, -t719, 0, 0, 0, 0, 0, 0, t739, t743, t738, -t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t760, 0, 0, 0, 0, 0, 0, t725, -t726, 0, t701, 0, 0, 0, 0, 0, 0, t725, 0, t726, t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t748, -t746, t777, t714, 0, 0, 0, 0, 0, 0, t748, t777, t746, t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t745, -t747, -t778, t710, 0, 0, 0, 0, 0, 0, t745, -t778, t747, -t705; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t739, t738, -t743, t719, 0, 0, 0, 0, 0, 0, -t739, -t743, -t738, t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t748, t777, t746, t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t739, -t743, -t738, t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t745, t778, -t747, t705;];
f_new_reg = t1;
