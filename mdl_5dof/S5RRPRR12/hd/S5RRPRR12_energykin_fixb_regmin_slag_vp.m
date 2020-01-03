% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:36
% EndTime: 2019-12-31 20:30:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (141->39), mult. (322->86), div. (0->0), fcn. (183->6), ass. (0->33)
t149 = qJD(1) ^ 2;
t159 = t149 / 0.2e1;
t148 = cos(qJ(2));
t158 = t148 * t149;
t145 = sin(qJ(2));
t156 = qJD(1) * t145;
t154 = pkin(6) * t156 + qJD(3);
t124 = -pkin(7) * t156 + (-pkin(2) - pkin(3)) * qJD(2) + t154;
t155 = qJD(1) * t148;
t132 = pkin(6) * t155 + qJD(2) * qJ(3);
t129 = -pkin(7) * t155 + t132;
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t157 = t144 * t124 + t147 * t129;
t153 = qJD(1) * qJD(2);
t130 = -qJD(1) * pkin(1) - pkin(2) * t155 - qJ(3) * t156;
t152 = t145 * t153;
t151 = t148 * t153;
t123 = pkin(3) * t155 - t130;
t150 = t147 * t124 - t144 * t129;
t127 = (t144 * t145 + t147 * t148) * qJD(1);
t146 = cos(qJ(5));
t143 = sin(qJ(5));
t139 = qJD(2) - qJD(4);
t131 = -qJD(2) * pkin(2) + t154;
t128 = (-t144 * t148 + t145 * t147) * qJD(1);
t126 = qJD(5) + t127;
t121 = t146 * t128 - t143 * t139;
t120 = t143 * t128 + t146 * t139;
t119 = -t139 * pkin(8) + t157;
t118 = t139 * pkin(4) - t150;
t117 = t127 * pkin(4) - t128 * pkin(8) + t123;
t1 = [t159, 0, 0, t145 ^ 2 * t159, t145 * t158, t152, t151, qJD(2) ^ 2 / 0.2e1, pkin(1) * t158 - pkin(6) * t152, -t149 * pkin(1) * t145 - pkin(6) * t151, -t131 * qJD(2) - t130 * t155, (t131 * t145 + t132 * t148) * qJD(1), t132 * qJD(2) - t130 * t156, t132 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1, t128 ^ 2 / 0.2e1, -t128 * t127, -t128 * t139, t127 * t139, t139 ^ 2 / 0.2e1, t123 * t127 - t150 * t139, t123 * t128 + t157 * t139, t121 ^ 2 / 0.2e1, -t121 * t120, t121 * t126, -t120 * t126, t126 ^ 2 / 0.2e1, (t146 * t117 - t143 * t119) * t126 + t118 * t120, -(t143 * t117 + t146 * t119) * t126 + t118 * t121;];
T_reg = t1;
