% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:21
% EndTime: 2019-12-31 19:36:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (161->37), mult. (412->79), div. (0->0), fcn. (253->6), ass. (0->35)
t152 = pkin(3) + pkin(7);
t141 = qJD(1) ^ 2;
t151 = t141 / 0.2e1;
t150 = pkin(6) + qJ(3);
t140 = cos(qJ(2));
t149 = t140 * t141;
t138 = sin(qJ(2));
t148 = qJD(1) * t138;
t131 = qJD(2) * pkin(2) - t150 * t148;
t147 = qJD(1) * t140;
t132 = t150 * t147;
t135 = sin(pkin(8));
t136 = cos(pkin(8));
t121 = t135 * t131 + t136 * t132;
t146 = qJD(1) * qJD(2);
t145 = t138 * t146;
t144 = t140 * t146;
t120 = t136 * t131 - t135 * t132;
t119 = -qJD(2) * qJ(4) - t121;
t143 = qJD(4) - t120;
t133 = qJD(3) + (-pkin(2) * t140 - pkin(1)) * qJD(1);
t128 = (t135 * t140 + t136 * t138) * qJD(1);
t142 = -t128 * qJ(4) + t133;
t139 = cos(qJ(5));
t137 = sin(qJ(5));
t127 = t135 * t148 - t136 * t147;
t126 = qJD(5) + t128;
t123 = t139 * qJD(2) + t137 * t127;
t122 = t137 * qJD(2) - t139 * t127;
t118 = -qJD(2) * pkin(3) + t143;
t117 = t127 * pkin(3) + t142;
t116 = -t127 * pkin(4) - t119;
t115 = t128 * pkin(4) - t152 * qJD(2) + t143;
t114 = t152 * t127 + t142;
t1 = [t151, 0, 0, t138 ^ 2 * t151, t138 * t149, t145, t144, qJD(2) ^ 2 / 0.2e1, pkin(1) * t149 - pkin(6) * t145, -t141 * pkin(1) * t138 - pkin(6) * t144, -t120 * t128 - t121 * t127, t121 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1, t118 * t128 + t119 * t127, t118 * qJD(2) - t117 * t127, -t119 * qJD(2) - t117 * t128, t117 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1, t123 ^ 2 / 0.2e1, -t123 * t122, t123 * t126, -t122 * t126, t126 ^ 2 / 0.2e1, (-t137 * t114 + t139 * t115) * t126 + t116 * t122, -(t139 * t114 + t137 * t115) * t126 + t116 * t123;];
T_reg = t1;
