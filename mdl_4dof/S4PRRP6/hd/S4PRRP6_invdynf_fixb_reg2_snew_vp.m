% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRP6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:51
% EndTime: 2019-12-31 16:30:51
% DurationCPUTime: 0.80s
% Computational Cost: add. (909->107), mult. (1923->121), div. (0->0), fcn. (1152->6), ass. (0->69)
t728 = sin(qJ(3));
t730 = cos(qJ(3));
t733 = qJD(2) ^ 2;
t713 = t728 * t733 * t730;
t709 = qJDD(3) - t713;
t721 = t728 ^ 2;
t732 = qJD(3) ^ 2;
t710 = t721 * t733 + t732;
t689 = t730 * t709 - t728 * t710;
t737 = t728 * qJDD(2);
t738 = qJD(2) * qJD(3);
t700 = 0.2e1 * t730 * t738 + t737;
t729 = sin(qJ(2));
t731 = cos(qJ(2));
t676 = t731 * t689 - t729 * t700;
t686 = t728 * t709 + t730 * t710;
t725 = sin(pkin(6));
t726 = cos(pkin(6));
t748 = t725 * t676 - t726 * t686;
t747 = t726 * t676 + t725 * t686;
t672 = t729 * t689 + t731 * t700;
t706 = t725 * g(1) - t726 * g(2);
t744 = t725 * t706;
t741 = -g(3) + qJDD(1);
t707 = -t726 * g(1) - t725 * g(2);
t692 = t731 * t707 + t729 * t741;
t684 = -t733 * pkin(2) + qJDD(2) * pkin(5) + t692;
t674 = t730 * t684 - t728 * t706;
t722 = t730 ^ 2;
t740 = t721 + t722;
t739 = t733 * (-pkin(3) * t730 - qJ(4) * t728);
t736 = t730 * qJDD(2);
t735 = t728 * t738;
t691 = -t729 * t707 + t731 * t741;
t683 = -qJDD(2) * pkin(2) - t733 * pkin(5) - t691;
t711 = -t722 * t733 - t732;
t708 = qJDD(3) + t713;
t705 = t740 * t733;
t704 = t731 * qJDD(2) - t729 * t733;
t703 = -t729 * qJDD(2) - t731 * t733;
t702 = t740 * qJDD(2);
t701 = -0.2e1 * t735 + t736;
t697 = t730 * t706;
t694 = t726 * t706;
t688 = -t728 * t708 + t730 * t711;
t685 = t730 * t708 + t728 * t711;
t682 = t731 * t702 - t729 * t705;
t681 = t729 * t702 + t731 * t705;
t679 = t726 * t682;
t678 = t725 * t682;
t675 = t731 * t688 - t729 * t701;
t671 = t729 * t688 + t731 * t701;
t670 = -t728 * t684 - t697;
t669 = -t729 * t691 + t731 * t692;
t668 = t731 * t691 + t729 * t692;
t667 = qJDD(4) + t697 - t732 * qJ(4) - qJDD(3) * pkin(3) + (t684 + t739) * t728;
t666 = -t732 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t730 * t739 + t674;
t665 = -(-t735 + t736) * pkin(3) + (pkin(3) * qJD(3) - (2 * qJD(4))) * t728 * qJD(2) + t683 - t700 * qJ(4);
t664 = t726 * t675 + t725 * t685;
t663 = t725 * t675 - t726 * t685;
t662 = -t728 * t670 + t730 * t674;
t661 = t730 * t670 + t728 * t674;
t660 = t731 * t662 + t729 * t683;
t659 = t729 * t662 - t731 * t683;
t658 = t730 * t666 + t728 * t667;
t657 = t728 * t666 - t730 * t667;
t656 = t731 * t658 + t729 * t665;
t655 = t729 * t658 - t731 * t665;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t726 * t707 - t744, 0, 0, 0, 0, 0, 0, t726 * t703, -t726 * t704, 0, t726 * t669 - t744, 0, 0, 0, 0, 0, 0, t664, -t747, t679, t726 * t660 + t725 * t661, 0, 0, 0, 0, 0, 0, t664, t679, t747, t726 * t656 + t725 * t657; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t725 * t707 + t694, 0, 0, 0, 0, 0, 0, t725 * t703, -t725 * t704, 0, t725 * t669 + t694, 0, 0, 0, 0, 0, 0, t663, -t748, t678, t725 * t660 - t726 * t661, 0, 0, 0, 0, 0, 0, t663, t678, t748, t725 * t656 - t726 * t657; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t741, 0, 0, 0, 0, 0, 0, t704, t703, 0, t668, 0, 0, 0, 0, 0, 0, t671, -t672, t681, t659, 0, 0, 0, 0, 0, 0, t671, t681, t672, t655; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, 0, 0, 0, 0, 0, 0, t703, -t704, 0, t669, 0, 0, 0, 0, 0, 0, t675, -t676, t682, t660, 0, 0, 0, 0, 0, 0, t675, t682, t676, t656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t706, 0, 0, 0, 0, 0, 0, 0, 0, 0, t706, 0, 0, 0, 0, 0, 0, -t685, t686, 0, -t661, 0, 0, 0, 0, 0, 0, -t685, 0, -t686, -t657; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t741, 0, 0, 0, 0, 0, 0, t704, t703, 0, t668, 0, 0, 0, 0, 0, 0, t671, -t672, t681, t659, 0, 0, 0, 0, 0, 0, t671, t681, t672, t655; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t733, -qJDD(2), 0, t692, 0, 0, 0, 0, 0, 0, t688, -t689, t702, t662, 0, 0, 0, 0, 0, 0, t688, t702, t689, t658; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t733, 0, t691, 0, 0, 0, 0, 0, 0, t701, -t700, t705, -t683, 0, 0, 0, 0, 0, 0, t701, t705, t700, -t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t706, 0, 0, 0, 0, 0, 0, t685, -t686, 0, t661, 0, 0, 0, 0, 0, 0, t685, 0, t686, t657; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t711, -t709, t736, t674, 0, 0, 0, 0, 0, 0, t711, t736, t709, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t708, -t710, -t737, t670, 0, 0, 0, 0, 0, 0, t708, -t737, t710, -t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t701, t700, -t705, t683, 0, 0, 0, 0, 0, 0, -t701, -t705, -t700, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t711, t736, t709, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t701, -t705, -t700, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t708, t737, -t710, t667;];
f_new_reg = t1;
