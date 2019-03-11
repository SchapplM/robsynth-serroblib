% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:10
% EndTime: 2019-03-09 01:56:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (210->42), mult. (445->88), div. (0->0), fcn. (259->6), ass. (0->36)
t156 = pkin(4) + pkin(8);
t147 = qJD(1) ^ 2;
t155 = t147 / 0.2e1;
t141 = sin(pkin(9));
t133 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t151 = -pkin(7) * qJD(1) + t133;
t126 = t151 * t141;
t142 = cos(pkin(9));
t127 = t151 * t142;
t144 = sin(qJ(4));
t146 = cos(qJ(4));
t154 = t146 * t126 + t144 * t127;
t153 = qJD(1) * t141;
t152 = qJD(1) * t142;
t136 = qJD(1) * qJ(2) + qJD(3);
t131 = pkin(3) * t153 + t136;
t150 = -t144 * t126 + t146 * t127;
t130 = -t144 * t153 + t146 * t152;
t118 = -qJD(4) * qJ(5) - t154;
t149 = qJD(5) - t150;
t148 = -t130 * qJ(5) + t131;
t145 = cos(qJ(6));
t143 = sin(qJ(6));
t139 = t142 ^ 2;
t138 = t141 ^ 2;
t137 = -qJD(1) * pkin(1) + qJD(2);
t129 = (t141 * t146 + t142 * t144) * qJD(1);
t128 = qJD(6) + t130;
t121 = t145 * qJD(4) + t143 * t129;
t120 = t143 * qJD(4) - t145 * t129;
t119 = t129 * pkin(4) + t148;
t117 = -qJD(4) * pkin(4) + t149;
t116 = t156 * t129 + t148;
t115 = -t129 * pkin(5) - t118;
t114 = t130 * pkin(5) - t156 * qJD(4) + t149;
t1 = [t155, 0, 0, t137 * qJD(1), t147 * qJ(2), qJ(2) ^ 2 * t155 + t137 ^ 2 / 0.2e1, t136 * t153, t136 * t152 (-t138 - t139) * t133 * qJD(1), t136 ^ 2 / 0.2e1 + (t138 / 0.2e1 + t139 / 0.2e1) * t133 ^ 2, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * qJD(4), -t129 * qJD(4), qJD(4) ^ 2 / 0.2e1, t150 * qJD(4) + t131 * t129, -t154 * qJD(4) + t131 * t130, t117 * t130 + t118 * t129, t117 * qJD(4) - t119 * t129, -t118 * qJD(4) - t119 * t130, t119 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1, t121 ^ 2 / 0.2e1, -t121 * t120, t121 * t128, -t120 * t128, t128 ^ 2 / 0.2e1 (t145 * t114 - t143 * t116) * t128 + t115 * t120 -(t143 * t114 + t145 * t116) * t128 + t115 * t121;];
T_reg  = t1;
